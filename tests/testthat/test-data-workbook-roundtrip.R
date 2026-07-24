# write_data() sources the control-sheet and bioenergetics-sheet object names
# and order from the canonical schema (R/0-column_schema.R). These guards assert
# the write -> read round-trip preserves every control and bioenergetics value on
# the bundled datasets, so the schema-driven assembly is faithful.

# Control scalars and bioenergetics scalars are single values recycled across
# species (or per-species vectors); read_data returns them at nspp width, so
# compare the original recycled to the read-back length.
expect_roundtrip_field <- function(orig, back, key) {
  a <- orig[[key]]
  if (is.null(a)) return(invisible())
  b <- back[[key]]
  expect_equal(as.numeric(b), rep(as.numeric(a), length.out = length(b)),
               info = key)
}

roundtrip <- function(d) {
  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f))
  suppressWarnings(suppressMessages(write_data(d, f)))
  suppressWarnings(suppressMessages(read_data(f)))
}

CONTROL_KEYS <- c("nspp", "styr", "endyr", "projyr", "nsex", "spawn_month",
                  "nages", "minage", "nlengths", "pop_wt_index", "ssb_wt_index",
                  "alpha_wt_len", "beta_wt_len", "pop_age_transition_index",
                  "sigma_rec", "other_food", "estDynamics")
BIO_KEYS <- c("Ceq", "Cindex", "Pvalue", "fday", "CA", "CB", "Qc", "Tco", "Tcm",
              "Tcl", "CK1", "CK4", "Diet_distribution", "Diet_comp_weights")

# GOA2018MS is not bundled (the multispecies GOA reference reuses GOA2018SS), so
# only the three shipped data lists are exercised.
for (nm in Filter(exists, c("BS2017SS", "BS2017MS", "GOA2018SS"))) {
  test_that(paste("control + bioenergetics round-trip on", nm), {
    d <- get(nm)
    back <- roundtrip(d)
    for (k in c(CONTROL_KEYS, BIO_KEYS)) expect_roundtrip_field(d, back, k)
    # The fleet_control sheet is unchanged by this step but must still survive.
    expect_equal(nrow(back$fleet_control), nrow(d$fleet_control))
    expect_true(all(d$fleet_control$Fleet_code %in% back$fleet_control$Fleet_code))
  })
}

test_that("Diet_distribution / Diet_comp_weights round-trip (bundled data leave them empty)", {
  # Both are NULL in the shipped datasets, so set them explicitly to exercise the
  # two bioenergetics rows added to the schema for the write_data label vector.
  d <- BS2017MS
  d$Diet_distribution <- c(0, 1, 0)
  d$Diet_comp_weights <- c(2, 3, 4)
  back <- roundtrip(d)
  expect_equal(as.numeric(back$Diet_distribution), c(0, 1, 0))
  expect_equal(as.numeric(back$Diet_comp_weights), c(2, 3, 4))
})

test_that("read_data parses the control sheet by name, not position", {
  # The schema (hence write_data) orders alpha/beta before the age-transition
  # index; read_data must recover the same values regardless of row order. Write,
  # reverse the non-dimension control rows, and confirm the read-back is
  # unchanged -- the property that makes the schema-driven row order safe.
  d <- BS2017MS
  f1 <- tempfile(fileext = ".xlsx"); on.exit(unlink(f1), add = TRUE)
  suppressWarnings(suppressMessages(write_data(d, f1)))
  sheets <- readxl::excel_sheets(f1)
  book <- lapply(sheets, function(s) as.data.frame(readxl::read_xlsx(f1, sheet = s)))
  names(book) <- sheets
  ctrl <- book$control
  scal <- ctrl[1:4, , drop = FALSE]                    # nspp/styr/endyr/projyr
  rest <- ctrl[nrow(ctrl):5, , drop = FALSE]           # reverse rows 5+
  book$control <- rbind(scal, rest)
  f2 <- tempfile(fileext = ".xlsx"); on.exit(unlink(f2), add = TRUE)
  writexl::write_xlsx(book, f2)
  d_ord <- suppressWarnings(suppressMessages(read_data(f1)))
  d_rev <- suppressWarnings(suppressMessages(read_data(f2)))
  for (k in c(CONTROL_KEYS, BIO_KEYS))
    if (!is.null(d_ord[[k]]))
      expect_equal(as.numeric(d_rev[[k]]), as.numeric(d_ord[[k]]), info = k)
})
