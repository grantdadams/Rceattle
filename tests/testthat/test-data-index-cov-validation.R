# data_check() validation of the covariance (MVN) survey index: Sigma must be
# usable by the likelihood, and a Sigma that no fleet will read must not pass
# silently. Pure validation (no fit), so unguarded/fast.
#
# Sigma already had presence / squareness / dimension / symmetry checks. The two
# added here close the cases that used to get through:
#   * symmetric but not positive definite -- failed later inside the TMB
#     objective, in a message naming neither the fleet nor Sigma;
#   * supplied for a fleet that is not MVN/MVNORM -- ignored by both data_check()
#     and .align_index_cov(), so the fleet quietly fit its default lognormal.

# data_check() runs inside fit_mod() after switches are filled and the HCR
# fields are set from build_hcr(), so reproduce that minimum rather than calling
# it on a raw fixture.
prep <- function(d) {
  d <- suppressMessages(Rceattle::switch_check(d))
  d$HCR <- "NoFishing"
  d
}

# Equicorrelated Sigma. rho < -1/(n-1) makes it symmetric but indefinite, and
# rho = 1 makes it singular -- the two cases that used to reach TMB.
equicor <- function(n, rho, sd = 20) {
  R <- matrix(rho, n, n)
  diag(R) <- 1
  diag(rep(sd, n)) %*% R %*% diag(rep(sd, n))
}

mvn_data <- function(nyrs = 8, rho = 0.3) {
  d <- make_test_data(nyrs = nyrs, nages = 5, seed = 42)
  srv <- d$fleet_control$Fleet_name == "Survey"
  d$fleet_control$Index_distribution[srv] <- "MVN"
  d$index_cov <- list(Survey = equicor(nyrs, rho))
  prep(d)
}

testthat::test_that("a positive definite Sigma passes", {
  testthat::expect_no_error(data_check(mvn_data(rho = 0.3)))
})

testthat::test_that("a symmetric but indefinite Sigma is rejected, naming the fleet", {
  # rho = -0.5 with n = 8 gives eigenvalues 1 + 7*(-0.5) = -2.5 and 1.5:
  # symmetric, correctly sized, and not positive definite.
  d <- mvn_data(rho = -0.5)
  testthat::expect_error(data_check(d), "not positive definite")
  testthat::expect_error(data_check(d), "Survey")
  # The smallest eigenvalue is reported so the user can see how far off it is.
  testthat::expect_error(data_check(d), "smallest eigenvalue")
})

testthat::test_that("a singular Sigma is rejected too", {
  # rho = 1 exactly: rank 1, positive SEMI-definite. Cholesky fails, and the
  # likelihood could not invert it either.
  testthat::expect_error(data_check(mvn_data(rho = 1)), "not positive definite")
})

testthat::test_that("index_cov for a fleet that is not MVN warns instead of being ignored", {
  # Index_distribution left at its default: the covariance would be dropped and
  # the fleet would fit lognormal, with nothing said.
  d <- make_test_data(nyrs = 8, nages = 5, seed = 42)
  d$index_cov <- list(Survey = equicor(8, 0.3))
  d <- prep(d)
  testthat::expect_warning(data_check(d), "not using Index_distribution")
  testthat::expect_warning(data_check(d), "Survey")
})

testthat::test_that("a covariance keyed by Fleet_code does not trigger the warning", {
  # index_cov may be keyed by Fleet_name OR Fleet_code, so the stray check has to
  # accept both or it would fire on a perfectly good model.
  d <- make_test_data(nyrs = 8, nages = 5, seed = 42)
  srv <- d$fleet_control$Fleet_name == "Survey"
  d$fleet_control$Index_distribution[srv] <- "MVN"
  code <- d$fleet_control$Fleet_code[srv]
  d$index_cov <- setNames(list(equicor(8, 0.3)), as.character(code))
  testthat::expect_no_warning(data_check(prep(d)))
})

# Note: the branch guarding a missing Index_distribution column is not reachable
# through the package's own path -- switch_check() always supplies the column
# before data_check() sees it. The initialisation of mvn_flts is kept because
# data_check() is also callable directly on a hand-built list.

testthat::test_that(".is_positive_definite accepts only usable matrices", {
  testthat::expect_true(.is_positive_definite(diag(3)))
  testthat::expect_true(.is_positive_definite(equicor(5, 0.3)))
  testthat::expect_false(.is_positive_definite(equicor(5, -0.5)))  # indefinite
  testthat::expect_false(.is_positive_definite(equicor(5, 1)))     # singular
  testthat::expect_false(.is_positive_definite(matrix(NA_real_, 2, 2)))
  testthat::expect_false(.is_positive_definite(matrix(1:6, 2, 3))) # not square
  testthat::expect_false(.is_positive_definite(matrix(numeric(0), 0, 0)))
})

testthat::test_that("the helper reads only the upper triangle, so symmetry must be checked first", {
  # chol() ignores the lower triangle, so an asymmetric matrix whose upper
  # triangle implies a positive definite one passes the helper. data_check()
  # is safe because it rejects on symmetry in the branch BEFORE the positive
  # definiteness branch -- this pins that dependency, so the ordering is not
  # rearranged and the helper is not "fixed" without moving the symmetry check.
  M <- matrix(c(4, 99, 1, 4), 2, 2)          # upper triangle implies [[4,1],[1,4]]
  testthat::expect_false(isTRUE(all.equal(M, t(M))))
  testthat::expect_true(.is_positive_definite(M))

  # Behaviourally: an asymmetric Sigma must be reported as asymmetric, not
  # waved through by the positive definiteness branch.
  n <- 8
  A <- equicor(n, 0.3)
  A[2, 1] <- A[2, 1] + 5000               # break symmetry, leave upper tri PD
  d <- mvn_data()
  d$index_cov <- list(Survey = A)
  testthat::expect_error(data_check(d), "not symmetric")
})
