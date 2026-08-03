# read_data() robustness: required sheets error clearly, optional sheets are
# skipped when absent, and a non-numeric control/bioenergetics cell errors by
# name instead of silently becoming NA.

# Read every sheet of a workbook into a named list (for patch-and-rewrite).
.read_all_sheets <- function(f) {
  sheets <- readxl::excel_sheets(f)
  stats::setNames(lapply(sheets, function(s)
    as.data.frame(readxl::read_xlsx(f, sheet = s))), sheets)
}

testthat::test_that("a missing required sheet errors clearly", {
  testthat::skip_if_not_installed("writexl")
  testthat::skip_if_not_installed("readxl")

  f <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(control = data.frame(
    Object = c("nspp", "styr", "endyr", "projyr"), sp1 = c(1, 1, 2, 3))), f)
  testthat::expect_error(Rceattle::read_data(f),
                         "Required sheet 'fleet_control' not found")
})

testthat::test_that("optional sheets are skipped when absent (minimal workbook reads)", {
  testthat::skip_if_not_installed("writexl")
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not_installed("Rceattle")

  f <- withr::local_tempfile(fileext = ".xlsx")
  suppressMessages(suppressWarnings(Rceattle::write_template(f)))
  wb <- .read_all_sheets(f)

  # Drop the optional predation / environment sheets a single-species workbook
  # does not need.
  for (s in c("bioenergetics_control", "env_data", "diet_data", "ration_data",
              "caal_data"))
    wb[[s]] <- NULL
  f2 <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(wb, f2)

  back <- suppressMessages(suppressWarnings(Rceattle::read_data(f2)))
  # Reads cleanly; the dropped optional elements are simply absent.
  testthat::expect_equal(back$nspp, 1)
  testthat::expect_true(!is.null(back$fleet_control))
  testthat::expect_null(back$env_data)
})

testthat::test_that("a non-numeric control cell errors by name instead of becoming NA", {
  testthat::skip_if_not_installed("writexl")
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not_installed("Rceattle")

  f <- withr::local_tempfile(fileext = ".xlsx")
  suppressMessages(suppressWarnings(Rceattle::write_template(f)))
  wb <- .read_all_sheets(f)

  ctl <- wb[["control"]]
  # Corrupt the sigma_rec value with a typo (a stray non-number).
  row <- which(ctl$Object == "sigma_rec")
  testthat::skip_if(length(row) == 0, "no sigma_rec row to corrupt")
  ctl[row, 2] <- "0.O"          # letter O instead of zero
  wb[["control"]] <- ctl
  f2 <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(wb, f2)

  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::read_data(f2))),
    "Non-numeric value.*sigma_rec")
})

testthat::test_that("a genuinely blank control cell stays NA and does not error", {
  testthat::skip_if_not_installed("writexl")
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not_installed("Rceattle")

  f <- withr::local_tempfile(fileext = ".xlsx")
  suppressMessages(suppressWarnings(Rceattle::write_template(f)))
  wb <- .read_all_sheets(f)

  ctl <- wb[["control"]]
  row <- which(ctl$Object == "other_food")
  testthat::skip_if(length(row) == 0, "no other_food row to blank")
  ctl[row, 2] <- NA          # a genuinely empty per-species cell
  wb[["control"]] <- ctl
  f2 <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(wb, f2)

  # Reads without error; the blanked scalar comes back NA (not a stop()).
  back <- suppressMessages(suppressWarnings(Rceattle::read_data(f2)))
  testthat::expect_true(is.na(back$other_food))
})
