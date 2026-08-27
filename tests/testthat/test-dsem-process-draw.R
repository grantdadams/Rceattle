# Redrawing a DSEM's process, and its covariate observations.
#
# PROVENANCE. Three defects sat behind one report -- "self_test(process = TRUE)
# does not simulate process error":
#
#   1. calculate_dsem() carried its own do_simulate branch, `x_tj = draw_tj`,
#      which assigned the WHOLE latent field and so discarded the conditional
#      draw sim_mod() had just substituted. It consulted neither the map nor
#      dsem_cond_k, so under family = "fixed" it also overwrote the covariate
#      columns -- the env_data series -- and each replicate was generated under
#      an environment the refit never saw (measured: a scaled series with sd 1
#      moved by up to 2.96). The draw is now taken only in R.
#
#   2. The template's rec_dev_drawn_sim mask is written inside its own IID draw,
#      which is gated off under a DSEM, so it stayed zero however much was
#      redrawn. .sim_process_truth() reads that mask, so attr(, "process_sim")
#      came back NULL on exactly the models whose process HAD been redrawn --
#      indistinguishable from nothing happening, and enough to stop
#      compare_sim() warning that it was measuring bias against the wrong thing.
#
#   3. init_dev was gated on dsem_on == 0, so the initial age structure kept its
#      fitted values while the hindcast deviations were re-realized.
#
# Plus one gap: a covariate under any family but "fixed" is a latent state
# OBSERVED WITH ERROR, and neither its state nor its observation was drawn.

.dsem_fit <- function(sem = NULL, family = "fixed", env = NULL, mode = 1) {
  d <- Rceattle::BS2017SS
  d$env_data <- if (is.null(env)) {
    data.frame(Year = d$styr:d$endyr, BT = 0)
  } else env(d)
  spec <- if (is.null(sem)) Rceattle::build_DSEM(family = family) else
    Rceattle::build_DSEM(sem = sem, family = family)
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = mode,
    random_rec = TRUE, msmMode = 0, dsem = spec,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
}

testthat::test_that("a DSEM process draw reaches the dynamics and reports itself", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  fit <- .dsem_fit()
  set.seed(3)
  sim <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE, process = TRUE)))
  truth <- attr(sim, "process_sim")

  testthat::expect_false(is.null(truth))
  # The latent field is the process a DSEM self test should be scored against;
  # rec_dev is what every existing recovery statistic already reads. Both.
  testthat::expect_true(all(c("rec_dev", "rec_dev_drawn", "dsem_x_tj",
                              "dsem_x_tj_drawn", "init_dev", "init_dev_drawn")
                            %in% names(truth)))
  testthat::expect_equal(dim(truth$dsem_x_tj),
                         dim(as.matrix(fit$estimated_params$dsem_x_tj)))

  # The field spans the hindcast (estimate_projection = FALSE), so exactly the
  # hindcast cells of rec_dev are marked and the projection cells are not.
  nyrs_hind <- fit$data_list$endyr - fit$data_list$styr + 1L
  testthat::expect_equal(sum(truth$rec_dev_drawn),
                         as.integer(fit$data_list$nspp * nyrs_hind))
  testthat::expect_true(all(truth$rec_dev_drawn[, seq_len(nyrs_hind)]))
  testthat::expect_false(any(truth$rec_dev_drawn[, -seq_len(nyrs_hind)]))

  # Drawn cells MOVED, measured against the deviations the fit actually used.
  # Not against estimated_params$rec_dev: under a DSEM that block is mapped out
  # and sits at zero, so comparing with it would only assert "the draw is not
  # zero" and would pass on a draw that ignored the fit entirely. The fitted
  # deviations are the ones a process = FALSE draw reports, since they are
  # derived from the latent field on every evaluation.
  fitted_dev <- Rceattle:::.sim_draw(
    Rceattle:::.sim_obj(fit), state = Rceattle:::.sim_state_codes(FALSE))$rec_dev_sim
  testthat::expect_gt(
    max(abs(truth$rec_dev[truth$rec_dev_drawn] -
              fitted_dev[truth$rec_dev_drawn])), 1e-6)
  # ...and the cells the draw did NOT touch are untouched.
  testthat::expect_equal(truth$rec_dev[!truth$rec_dev_drawn],
                         fitted_dev[!truth$rec_dev_drawn])

  # init_dev is drawn under a DSEM too. Its density is its own
  # N(-bias*R_sd^2/2, R_sd) whatever dsem_on is, and leaving it alone gave a
  # replicate whose initial age structure and hindcast came from different
  # realizations.
  testthat::expect_gt(sum(truth$init_dev_drawn), 0L)

  # process = FALSE draws nothing and says nothing.
  set.seed(3)
  plain <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE, process = FALSE)))
  testthat::expect_null(attr(plain, "process_sim"))
})

testthat::test_that("the substituted field is what the template uses", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # The regression test for defect 1. A second draw inside calculate_dsem()
  # overwrote the substituted field, so the R-side draw was inert -- adding 100
  # to every latent state changed nothing at all. Comparing two seeds cannot see
  # this, because the observation draw makes any two calls differ.
  fit <- .dsem_fit(mode = 3)
  obj <- Rceattle:::.sim_obj(fit)
  par <- obj$env$last.par.best
  bumped <- par
  bumped[names(par) == "dsem_x_tj"] <- par[names(par) == "dsem_x_tj"] + 100
  st <- Rceattle:::.sim_state_codes(TRUE)

  set.seed(1); a <- Rceattle:::.sim_draw(obj, state = st, par = par)
  set.seed(1); b <- Rceattle:::.sim_draw(obj, state = st, par = bumped)
  testthat::expect_gt(max(abs(a$rec_dev_sim - b$rec_dev_sim)), 50)
  testthat::expect_false(isTRUE(all.equal(a$index_hat, b$index_hat)))
})

testthat::test_that("a DSEM refuses the draw where the standard path does", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # srr_fun = 0 with srr_pred_fun > 0 scores the deviations a second time
  # through the stock-recruit penalty, so there is no single distribution to
  # draw from. That is true whether or not the first density is a GMRF, and the
  # SEM draw has to respect the same gate the template's own draw does.
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    recFun = Rceattle::build_srr(srr_fun = 0, srr_pred_fun = "BevertonHolt",
                                 srr_est_mode = 1),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  testthat::expect_equal(as.numeric(fit$quantities$rec_srr_single_density), 0)

  # "recruitment" rather than TRUE so the only warning is the gate under test;
  # TRUE would also warn that this model gives M and growth no distribution,
  # which is true and unrelated. Both spellings reach the same gate.
  testthat::expect_warning(
    sim <- suppressMessages(
      Rceattle::sim_mod(fit, simulate = TRUE, process = "recruitment")),
    "AMAK/Ianelli")
  testthat::expect_null(attr(sim, "process_sim"))

  # Nothing moved, rather than "nothing was reported".
  obj <- Rceattle:::.sim_obj(fit)
  off <- Rceattle:::.sim_draw(obj, state = Rceattle:::.sim_state_codes(FALSE))
  on  <- Rceattle:::.sim_draw(obj, state = Rceattle:::.sim_state_codes(TRUE))
  testthat::expect_equal(off$rec_dev_sim, on$rec_dev_sim)
  testthat::expect_equal(off$init_dev_sim, on$init_dev_sim)
})

testthat::test_that("a fixed covariate is data and is held", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  sem <- "
    recdevs1 -> recdevs1, 1, rho_R,     0.3
    temp     -> recdevs1, 0, temp_to_R, 0.2
    recdevs1 <-> recdevs1, 0, sigmaR1,  0.6
    recdevs2 <-> recdevs2, 0, sigmaR2,  0.6
    recdevs3 <-> recdevs3, 0, sigmaR3,  0.6
  "
  env <- function(d) {
    set.seed(1)
    yrs <- d$styr:d$endyr
    data.frame(Year = yrs,
               temp = as.numeric(scale(cumsum(stats::rnorm(length(yrs))))))
  }
  fit <- .dsem_fit(sem = sem, family = "fixed", env = env)

  set.seed(5)
  sim <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE, process = TRUE)))

  # family = "fixed" means the covariate column IS env_data: familycode 0, no
  # measurement density, pinned in the map. It must come back untouched --
  # a redrawn covariate beside the ORIGINAL env_data handed to the refit is a
  # model/data mismatch in every replicate.
  testthat::expect_equal(sim$env_data$temp, fit$data_list$env_data$temp)

  drawn <- attr(sim, "process_sim")$dsem_x_tj_drawn
  rc <- fit$dsem$tmb_inputs$data$rec_dev_col + 1L
  cov_col <- setdiff(seq_len(ncol(drawn)), rc)
  testthat::expect_length(cov_col, 1L)
  testthat::expect_false(any(drawn[, cov_col]))
  testthat::expect_true(all(drawn[, rc]))
})

testthat::test_that("a covariate observed with error is drawn, state and observation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # Any family but "fixed" makes the covariate a latent state observed with
  # error. Both halves then belong in a simulation: the STATE is process (drawn
  # with the rest of the field) and the OBSERVATION is data (drawn whenever the
  # observations are, process error or not, and written back into env_data so
  # the refit is fitted to it).
  sem <- "
    recdevs1 -> recdevs1, 1, rho_R,     0.3
    temp     -> temp,     1, rho_T,     0.4
    temp     -> recdevs1, 0, temp_to_R, 0.2
    temp     <-> temp,    0, sigmaT,    0.5
    recdevs1 <-> recdevs1, 0, sigmaR1,  0.6
    recdevs2 <-> recdevs2, 0, sigmaR2,  0.6
    recdevs3 <-> recdevs3, 0, sigmaR3,  0.6
  "
  env <- function(d) {
    set.seed(1)
    yrs <- d$styr:d$endyr
    data.frame(Year = yrs,
               temp = as.numeric(scale(cumsum(stats::rnorm(length(yrs))))))
  }
  fit <- .dsem_fit(sem = sem, family = "normal", env = env, mode = 3)

  set.seed(5)
  sim <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE, process = TRUE)))
  testthat::expect_false(isTRUE(all.equal(sim$env_data$temp,
                                          fit$data_list$env_data$temp)))
  # Written back in place: same years, same rows, no NA introduced.
  testthat::expect_identical(sim$env_data$Year, fit$data_list$env_data$Year)
  testthat::expect_false(anyNA(sim$env_data$temp))

  drawn <- attr(sim, "process_sim")$dsem_x_tj_drawn
  testthat::expect_true(all(drawn))   # every column is latent here

  # The observation is data, so it is redrawn without process error too.
  set.seed(5)
  obs_only <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE, process = FALSE)))
  testthat::expect_false(isTRUE(all.equal(obs_only$env_data$temp,
                                          fit$data_list$env_data$temp)))
  testthat::expect_null(attr(obs_only, "process_sim"))
})

testthat::test_that("the covariate observation draw follows its measurement SD", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # The draw has to be the density that scores it, or a self test measures
  # recovery against a different data-generating process. Checked as a moment:
  # y_sim - mu over many draws is centred on 0 with the fitted measurement SD.
  sem <- "
    recdevs1 -> recdevs1, 1, rho_R,     0.3
    temp     -> temp,     1, rho_T,     0.4
    temp     -> recdevs1, 0, temp_to_R, 0.2
    temp     <-> temp,    0, sigmaT,    0.5
    recdevs1 <-> recdevs1, 0, sigmaR1,  0.6
    recdevs2 <-> recdevs2, 0, sigmaR2,  0.6
    recdevs3 <-> recdevs3, 0, sigmaR3,  0.6
  "
  env <- function(d) {
    set.seed(1)
    yrs <- d$styr:d$endyr
    data.frame(Year = yrs,
               temp = as.numeric(scale(cumsum(stats::rnorm(length(yrs))))))
  }
  fit <- .dsem_fit(sem = sem, family = "normal", env = env, mode = 3)
  sigma <- exp(as.numeric(fit$estimated_params$dsem_lnsigma_z))[1]

  obj <- Rceattle:::.sim_obj(fit)
  none <- Rceattle:::.sim_state_codes(FALSE)
  cov_col <- ncol(as.matrix(fit$estimated_params$dsem_x_tj))
  mu <- as.matrix(Rceattle:::.sim_draw(obj, state = none)$dsem_z_tj)[, cov_col]

  set.seed(11)
  resid <- as.numeric(replicate(200, {
    d <- Rceattle:::.sim_draw(obj, state = none)
    as.matrix(d$dsem_y_tj_sim)[, cov_col] - mu
  }))
  testthat::expect_equal(stats::sd(resid), sigma, tolerance = 0.08)
  testthat::expect_lt(abs(mean(resid)), 4 * sigma / sqrt(length(resid)))
})
