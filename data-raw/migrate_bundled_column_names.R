# =============================================================================
# Migrate the bundled datasets to the canonical workbook column names
# =============================================================================
#
# The 4.10.0 column renames (Q_prior -> Catchability_init, Comp_loglike ->
# Comp_distribution, ...) were applied to BS2017SS, BS2017MS and GOA2018SS but
# not to the other bundled datasets. Those still shipped the legacy spellings,
# so simply loading one -- GOApollock, Atka2022, ... -- printed a dozen
# deprecation messages before the model did anything.
#
# The rename goes through the package's own schema-driven upgraders, so the
# mapping cannot drift from `.rce_column_schema()`'s `aliases` field. Only
# NAMES change; every value is carried over untouched, which the verification
# at the bottom asserts.
#
# Re-run with:  Rscript data-raw/migrate_bundled_column_names.R
# =============================================================================

suppressMessages(pkgload::load_all(".", quiet = TRUE))

upgrade_fc   <- Rceattle:::.rce_upgrade_fleet_control_aliases
upgrade_list <- Rceattle:::.rce_upgrade_data_list_aliases

datasets <- sub("[.]rda$", "", list.files("data", pattern = "[.]rda$"))

for (ds in datasets) {
  f <- file.path("data", paste0(ds, ".rda"))
  e <- new.env()
  load(f, envir = e)
  if (!ds %in% ls(e)) next
  before <- get(ds, envir = e)
  if (!is.list(before)) next

  after <- suppressMessages(upgrade_list(before))
  if (!is.null(after$fleet_control)) {
    after$fleet_control <- suppressMessages(upgrade_fc(after$fleet_control))
  }

  # `est_M1` -> `M1_model` is handled by a bespoke helper rather than the schema
  # aliases, so the upgraders above leave it alone. Mirror .alias_est_M1()'s
  # exact semantics: an all-NA est_M1 is a placeholder meaning "not specified"
  # and must be dropped so M1_model still falls through to its default; a real
  # value is promoted only when M1_model is unset.
  if (!is.null(after$est_M1)) {
    if (is.null(after$M1_model) && !all(is.na(after$est_M1))) {
      after$M1_model <- after$est_M1
    }
    after$est_M1 <- NULL
  }

  if (identical(before, after)) next

  # Values must be preserved exactly -- this is a rename, not a data edit.
  # The upgrader also reorders the renamed columns into canonical schema order,
  # so compare BY NAME (old -> canonical via the schema aliases), not by
  # position.
  fb <- before$fleet_control
  fa <- after$fleet_control
  if (!is.null(fb)) {
    stopifnot(nrow(fb) == nrow(fa), ncol(fb) == ncol(fa))
    canon_of <- function(nm) {
      for (row in Rceattle:::.rce_column_schema()) {
        if (identical(row$name, nm)) return(nm)
        if (!is.null(row$aliases) && nm %in% row$aliases) return(row$name)
      }
      nm
    }
    for (nm in colnames(fb)) {
      target <- canon_of(nm)
      stopifnot(target %in% colnames(fa),
                identical(unname(fb[[nm]]), unname(fa[[target]])))
    }
  }

  assign(ds, after, envir = e)
  save(list = ds, file = f, envir = e, compress = "xz")
  cat(sprintf("migrated %-22s (%d fleet_control cols)\n", ds,
              if (is.null(fa)) 0L else ncol(fa)))
}

# ---- Verify no legacy spelling survives -------------------------------------
legacy_fc <- c("Time_varying_sel_sd_prior", "Sel_norm_bin1", "Sel_norm_bin2",
               "Comp_loglike", "CAAL_loglike", "Weight1_Numbers2", "Q_index",
               "Q_prior", "Q_sd_prior", "Time_varying_q_sd_prior",
               "Index_sd_prior", "Catch_sd_prior")
legacy_top <- c("proj_F_prop", "sigma_rec_prior", "Diet_loglike", "est_M1")

left <- 0L
for (ds in datasets) {
  e <- new.env(); load(file.path("data", paste0(ds, ".rda")), envir = e)
  if (!ds %in% ls(e)) next
  x <- get(ds, envir = e)
  if (!is.list(x)) next
  hits <- c(intersect(legacy_fc, colnames(x$fleet_control)),
            intersect(legacy_top, names(x)))
  if (length(hits)) {
    left <- left + length(hits)
    cat(sprintf("  STILL LEGACY %-20s %s\n", ds, paste(hits, collapse = ", ")))
  }
}
cat("remaining legacy names across all bundled datasets:", left, "\n")
