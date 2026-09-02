# What the composition draw is scaled and dispersed by: Comp_distribution
# (multinomial against Dirichlet-multinomial) and Comp_weights, which multiplies
# the likelihood and is therefore an effective sample size.
#
# Split out of test-functions-sim-mod-comp.R, which was 27% of the coverage
# suite in one file. Shared fixtures: helpers-comp-sim.R.
testthat::skip_on_cran()

testthat::test_that("the Dirichlet-multinomial composition is over-dispersed", {
  testthat::skip_if_not_installed("TMB")

  # theta is exp(Comp_weights) and Comp_weights is ESTIMATED, so a fitted model
  # drives it large, where the Dirichlet-multinomial is a multinomial to machine
  # precision and its expected over-dispersion is 1.0000 -- a check run there
  # passes whatever the simulator does. Hold every parameter at one fit and vary
  # only theta, exactly as the diet harness has to.
  base <- .comp_fixture(dist = "Multinomial", weight = 1)
  hold <- function(dist, w) {
    d <- base$data_list
    d$fleet_control$Comp_distribution <- dist
    inits <- base$estimated_params
    inits$comp_weights[] <- w
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, inits = inits, file = NULL, estimateMode = 3, msmMode = 0,
      random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE))))
  }
  mn <- hold("Multinomial", 1)            # weight 1: draw at the nominal N
  dm <- hold("DirichletMultinomial", 0)   # theta = exp(0) = 1

  cols <- grep("^Comp_", names(mn$data_list$comp_data))
  hat  <- mn$quantities$comp_hat[, seq_along(cols), drop = FALSE]
  N    <- mn$data_list$comp_data$Sample_size[1]
  big  <- which(hat[1, ] > 0.05)
  testthat::expect_gt(length(big), 1)

  sd_mn <- apply(t(.comp_draws(mn, 600)[1, , ])[, big, drop = FALSE], 2, stats::sd)
  sd_dm <- apply(t(.comp_draws(dm, 600)[1, , ])[, big, drop = FALSE], 2, stats::sd)

  expected <- sqrt((N + N * 1) / (1 + N * 1))   # ~1.41 at N = 200, theta = 1
  testthat::expect_gt(expected, 1.3)            # the check has power
  testthat::expect_equal(mean(sd_dm / sd_mn), expected, tolerance = 0.2)
})


testthat::test_that("Comp_weights scales the multinomial draw", {
  testthat::skip_if_not_installed("TMB")

  # The weight multiplies the multinomial log-likelihood, which makes it an
  # effective sample size. Drawing at the nominal N regardless would hand the
  # estimator data of the wrong precision.
  w <- 4
  f1 <- .comp_fixture(weight = 1)
  fw <- .comp_fixture(weight = w)

  cols <- grep("^Comp_", names(f1$data_list$comp_data))
  N    <- f1$data_list$comp_data$Sample_size[1]

  r1 <- t(.comp_draws(f1, 400)[1, , ])
  rw <- t(.comp_draws(fw, 400)[1, , ])
  testthat::expect_equal(unname(rowSums(r1)[1]), N)
  testthat::expect_equal(unname(rowSums(rw)[1]), N * w)

  hat <- f1$quantities$comp_hat[1, seq_along(cols)]
  big <- which(hat > 0.05)
  # Proportions, so the extra sample size shows as a sqrt(w) tighter spread.
  testthat::expect_equal(
    mean(apply(rw[, big, drop = FALSE] / (N * w), 2, stats::sd) /
         apply(r1[, big, drop = FALSE] / N, 2, stats::sd)),
    1 / sqrt(w), tolerance = 0.15)
})
