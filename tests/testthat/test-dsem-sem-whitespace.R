# A sem written with its columns aligned must mean what it reads as.
#
# dsem::make_dsem_ram() decides whether a variable already has a variance of its
# own by comparing the path string it would add -- paste(v, "<->", v),
# single-spaced -- against the paths already in the model, with no whitespace
# normalization. So `recdevs1  <->  recdevs1` failed to match, dsem appended a
# second variance row for the same node, and MakeADFun() then died inside the
# dsem DLL: a SEGFAULT on a multi-variable sem, an error on a single-variable
# one. Two of the three live pollock DSEM scripts were written that way and
# could not run at all.
#
# Nothing in this package could catch it: every sem in tests/, examples/ and
# vignettes/dsem.Rmd happens to be single-spaced, so they all take the safe path.
# These fixtures are deliberately aligned.

.ws_data <- function(nspp = 1L) {
  d <- list(styr = 1990L, endyr = 2010L, projyr = 2015L, nspp = nspp,
            spnames = letters[seq_len(nspp)], sigma_rec = 1, random_rec = TRUE,
            proj_mean_rec = TRUE)
  ed <- data.frame(Year = 1990:2015)
  ed$BT <- as.numeric(seq_along(ed$Year))
  ed$CP <- as.numeric(rev(seq_along(ed$Year)))
  d$env_data <- ed
  d
}

testthat::test_that("the arrow normalizer puts one space around every arrow", {
  n <- Rceattle:::.dsem_normalize_arrows
  testthat::expect_equal(n("recdevs1  <->  recdevs1, 0, sigmaR1, 1"),
                         "recdevs1 <-> recdevs1, 0, sigmaR1, 1")
  testthat::expect_equal(n("BT           ->  recdevs1, 1, a, 0"),
                         "BT -> recdevs1, 1, a, 0")
  testthat::expect_equal(n("a<->b, 0, p, 1"), "a <-> b, 0, p, 1")
  testthat::expect_equal(n("a <- b, 0, p, 1"), "a <- b, 0, p, 1")
  # `<->` must be protected before `->` and `<-` are touched, or the three
  # arrows rewrite each other's halves.
  testthat::expect_equal(n("a<->b, 0, p, 1\nc->d, 1, q, 0"),
                         "a <-> b, 0, p, 1\nc -> d, 1, q, 0")
  testthat::expect_null(n(NULL))
})

testthat::test_that("an aligned sem builds the same RAM as a compact one", {
  testthat::skip_if_not_installed("dsem")
  d <- .ws_data()
  compact <- "BT -> BT, 1, AR_BT, 0
BT -> recdevs1, 1, BT_to_R, 0
recdevs1 <-> recdevs1, 0, sigmaR1, 1"
  aligned <- "
  BT         ->   BT,         1,   AR_BT,     0
  BT         ->   recdevs1,   1,   BT_to_R,   0
  recdevs1  <->   recdevs1,   0,   sigmaR1,   1
"
  b1 <- suppressWarnings(Rceattle::build_dsem_objects(
    Rceattle::build_DSEM(sem = compact, family = "fixed"), data_list = d))
  b2 <- suppressWarnings(Rceattle::build_dsem_objects(
    Rceattle::build_DSEM(sem = aligned, family = "fixed"), data_list = d))

  cols <- c("name", "parameter", "first", "second", "direction", "lag")
  testthat::expect_equal(b2$sem_full[, cols], b1$sem_full[, cols])
  # The specific symptom: exactly ONE variance term on recdevs1, not two.
  self <- b2$sem_full[b2$sem_full$first == "recdevs1" &
                        b2$sem_full$second == "recdevs1" &
                        abs(as.numeric(b2$sem_full$direction)) == 2, ]
  testthat::expect_equal(nrow(self), 1L)
  testthat::expect_equal(as.character(self$name), "sigmaR1")
  testthat::expect_equal(length(b2$tmb_inputs$parameters$beta_z),
                         length(b1$tmb_inputs$parameters$beta_z))
})

testthat::test_that("a multi-variable aligned sem builds -- the segfault case", {
  testthat::skip_if_not_installed("dsem")
  # Verbatim shape of Rceattle-models/GOA pollock/2025/dsem.R, which took the
  # whole R session down before the normalization.
  d <- .ws_data()
  sem <- "
  BT  ->  BT,        1,   BT_AR1,    0
  CP  ->  CP,        1,   CP_AR1,    0
  recdevs1  <->  recdevs1,   0,   sigmaR1,   1
"
  b <- suppressWarnings(Rceattle::build_dsem_objects(
    Rceattle::build_DSEM(sem = sem, family = "fixed"), data_list = d))
  self <- b$sem_full[b$sem_full$first == b$sem_full$second &
                       abs(as.numeric(b$sem_full$direction)) == 2 &
                       as.numeric(b$sem_full$lag) == 0, ]
  testthat::expect_false(any(duplicated(as.character(self$first))))
})

testthat::test_that("a genuinely duplicated variance term is named, not fatal", {
  testthat::skip_if_not_installed("dsem")
  d <- .ws_data()
  sem <- "recdevs1 <-> recdevs1, 0, sigmaR1, 1
recdevs1 <-> recdevs1, 0, sigmaR1b, 1"
  testthat::expect_error(
    suppressWarnings(Rceattle::build_dsem_objects(
      Rceattle::build_DSEM(sem = sem, family = "fixed"), data_list = d)),
    "two variance terms")
})

testthat::test_that("check_dsem_spec() refuses a specification instead of passing it", {
  d <- .ws_data()
  spec <- Rceattle::build_DSEM(sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                               family = "fixed")
  # The natural mistake: `dsem` is what fit_mod() takes, so a user passes the
  # specification. Reporting "OK, no checks run" on a spec nothing has looked at
  # is worse than useless.
  testthat::expect_error(Rceattle::check_dsem_spec(d, spec),
                         "build_dsem_objects")
})
