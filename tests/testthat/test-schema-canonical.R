# Canonical workbook-column schema (R/0-column_schema.R): faithful-reproduction
# and integrity guards. These prove the single-source schema reproduces the
# hand-maintained copies BEFORE any consumer is switched over to it.

# Compare descriptions on their meaningful text: collapse CR/LF and runs of
# spaces so a drift guard flags wording changes, not cosmetic whitespace.
norm_desc <- function(x) {
  x <- gsub("[\r\n]+", "\n", x)
  x <- gsub("[ \t]+", " ", x)
  x <- gsub(" ?\n ?", "\n", x)
  trimws(x)
}

test_that("schema is structurally sound", {
  schema <- .rce_column_schema()
  nms <- vapply(schema, function(r) r$name, character(1))

  # Unique canonical names; names() keyed correctly.
  expect_equal(anyDuplicated(nms), 0L)
  expect_equal(names(schema), unname(nms))

  # Valid sheets.
  sheets <- vapply(schema, function(r) r$sheet, character(1))
  expect_true(all(sheets %in% c("control", "fleet_control", "bioenergetics_control")))

  # No alias collides with a canonical name (would make a rename ambiguous).
  aliases <- unlist(lapply(schema, function(r) r$aliases))
  expect_length(intersect(aliases, nms), 0)

  # No alias is shared by two columns (would make an upgrade ambiguous).
  expect_equal(anyDuplicated(aliases), 0L)

  # The orphan set is pinned: these are documented historically but read
  # nowhere and slated for removal. Promoting one to "live" (or adding a new
  # orphan) must be a deliberate, visible change -- pinning catches it, so the
  # file's "explicit, testable record of the drop" claim actually holds.
  orphan <- vapply(Filter(function(r) r$status == "orphan", schema),
                   function(r) r$name, character(1))
  expect_setequal(unname(orphan), c("R_sexr", "Catch_units", "Sex"))

  # Orphans never leak into the live accessors used by the generators.
  expect_length(intersect(.rce_schema_names("fleet_control"), orphan), 0)
  expect_length(intersect(.rce_schema_names("control"), orphan), 0)
})

test_that("committed meta_data_names.xlsx matches the schema generator exactly", {
  # The bundled xlsx is generated from the schema (data-raw/regenerate_meta_xlsx.R);
  # this drift guard fails if the two diverge (i.e. the schema was edited without
  # re-running the generator, or the asset was hand-edited).
  gen <- .rce_build_meta_data_df()
  xlsx <- as.data.frame(readxl::read_xlsx(
    system.file("extdata", "meta_data_names.xlsx", package = "Rceattle")))

  expect_equal(dim(xlsx), dim(gen))
  expect_equal(colnames(xlsx), colnames(gen))

  blank <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
  expect_equal(blank(xlsx[["Sheet name"]]), blank(gen[["Sheet name"]]))
  expect_equal(blank(xlsx[["Column/row name"]]), blank(gen[["Column/row name"]]))
  # Descriptions compared on meaningful text (whitespace-normalized), row-aligned.
  # USE.NAMES = FALSE so the comparison is on the normalized text, not names
  # carrying the raw strings (which would defeat the whitespace tolerance).
  expect_equal(vapply(blank(xlsx$Description), norm_desc, character(1), USE.NAMES = FALSE),
               vapply(blank(gen$Description), norm_desc, character(1), USE.NAMES = FALSE))
})

test_that("schema pins the known legacy -> canonical column aliases", {
  # These are the deprecated names accepted on read and upgraded in place
  # (read_data / switch_check). Pinning them here guards the alias table the
  # generator and the rename cascades rely on.
  expected <- list(
    Estimate_index_sd   = "Estimate_survey_sd",
    Index_sd            = c("Survey_sd_prior", "Index_sd_prior"),
    Catch_sd            = "Catch_sd_prior",
    Catchability_index    = "Q_index",
    Catchability_init     = c("Q_prior", "Q_init"),
    Catchability_prior_sd = "Q_sd_prior",
    Comp_distribution     = "Comp_loglike",
    CAAL_distribution     = "CAAL_loglike",
    Index_distribution    = "Index_loglike",
    N_sel_bins          = "Nselages",
    Time_varying_sel_sd = c("Sel_sd_prior", "Time_varying_sel_sd_prior"),
    Time_varying_q_sd   = "Time_varying_q_sd_prior",
    Catchability        = "Estimate_q",
    Bin_first_selected  = "Age_first_selected",
    Sel_norm_bin1       = "Age_max_selected",
    Sel_norm_bin2       = "Age_max_selected_upper",
    Sel_pen_first_bin   = "Sel_pen_first_age",
    Sel_pen_last_bin    = "Sel_pen_last_age",
    Sel_cap_bin         = "Sel_cap_age"
  )
  schema <- .rce_column_schema()
  for (canon in names(expected))
    expect_setequal(schema[[canon]]$aliases, expected[[canon]])
})
