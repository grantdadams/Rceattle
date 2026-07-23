# Fleets sharing a Q_index estimate one catchability and one deviate vector, so
# the q prior and the deviate penalties must be accumulated once (flt_q_lead),
# not once per sharing fleet.

testthat::skip_on_cran()

q_prior_dat <- function(share) {
  data("BS2017SS", envir = environment())
  d <- BS2017SS
  srv <- which(d$fleet_control$Fleet_type == 2)[1:2]
  d$fleet_control$Catchability[srv] <- "Estimated-with-prior"
  d$fleet_control$Q_init[srv]      <- 1
  d$fleet_control$Q_sd_prior[srv]   <- 0.2
  if (share) d$fleet_control$Q_index[srv] <- 99
  list(data = d, codes = d$fleet_control$Fleet_code[srv])
}

q_prior_fit <- function(share) {
  s <- q_prior_dat(share)
  fit <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(data_list = s$data, inits = NULL, estimateMode = 3, msmMode = 0,
                      fit_control = fit_control(phase = FALSE, verbose = 0))))
  list(prior = fit$quantities$jnll_comp[7, s$codes],
       qmap  = fit$map$mapList$index_log_q[s$codes],
       lead  = fit$obj$env$data$flt_q_lead[s$codes])
}

testthat::test_that("a shared Q_index counts the catchability prior once", {
  testthat::skip_if_not_installed("TMB")

  sep <- q_prior_fit(FALSE)
  shr <- q_prior_fit(TRUE)

  # Separate: two free q parameters, each carrying its own prior.
  testthat::expect_false(sep$qmap[[1]] == sep$qmap[[2]])
  testthat::expect_equal(unname(sep$lead), c(1L, 1L))
  testthat::expect_equal(sum(sep$prior), 2 * sep$prior[[1]])

  # Shared: ONE free q parameter, so the prior is accumulated on the lead only.
  testthat::expect_equal(shr$qmap[[1]], shr$qmap[[2]])
  testthat::expect_equal(unname(shr$lead), c(1L, 0L))
  testthat::expect_equal(shr$prior[[2]], 0)
  testthat::expect_equal(sum(shr$prior), sep$prior[[1]])
})

testthat::test_that("PowerEquation catchability is rejected as not implemented", {
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS", envir = environment())
  d <- BS2017SS
  d$fleet_control$Catchability[which(d$fleet_control$Fleet_type == 2)[1]] <- "PowerEquation"
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(data_list = d, inits = NULL, estimateMode = 3,
                                       msmMode = 0,
                                       fit_control = fit_control(phase = FALSE, verbose = 0))),
    regexp = "PowerEquation"
  )
})
