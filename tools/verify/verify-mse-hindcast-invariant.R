# tools/verify/verify-mse-hindcast-invariant.R
# Correctness invariant (independent of the .refit_like refactor): an MSE must
# not perturb the hindcast. With regenerate_past = FALSE, the operating model's
# spawning biomass over styr:endyr is fixed history and must be reproduced by the
# advanced OM_use, regardless of simulate_data / sample_rec (those add data only
# in the projection). For multispecies this only holds if the predation
# suitability window is pinned (the site-9 suit_styr/suit_endyr = om$ pins).
#
# Reports max |OM_use.ssb - OM.ssb| over the hindcast for SS and MS, under both
# (simulate=FALSE, sample=FALSE) and (simulate=TRUE, sample=TRUE). EM_use is
# reported for information: an EM legitimately updates its hindcast estimate as
# the assessment advances, so only the OM figure is the strict invariant.
#
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-mse-hindcast-invariant.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# First hind_nyrs columns of a fit's ssb (column 1 == styr), as a matrix.
hind_ssb <- function(fit, hind_nyrs) {
  s <- fit$quantities$ssb
  if (is.null(dim(s))) s <- matrix(s, nrow = 1)
  hind_nyrs <- min(hind_nyrs, ncol(s))
  s[, seq_len(hind_nyrs), drop = FALSE]
}

check_mse <- function(label, om, em) {
  hind_nyrs <- om$data_list$endyr - om$data_list$styr + 1
  om_ref <- hind_ssb(om, hind_nyrs)
  em_ref <- hind_ssb(em, hind_nyrs)
  for (flags in list(c(FALSE, FALSE), c(TRUE, TRUE))) {
    m <- tryCatch(
      Rceattle::run_mse(om = om, em = em, nsim = 1, assessment_period = 1,
                        sampling_period = 1, simulate_data = flags[1],
                        sample_rec = flags[2], regenerate_past = FALSE, seed = 99),
      error = function(e) e)
    if (inherits(m, "error")) {
      cat(sprintf("  [%s sim=%-5s samp=%-5s] ERROR: %s\n",
                  label, flags[1], flags[2], conditionMessage(m)))
      next
    }
    for (nm in names(m)) {
      s <- m[[nm]]
      d_om <- tryCatch(max(abs(hind_ssb(s$OM, hind_nyrs) - om_ref)), error = function(e) NA_real_)
      d_em <- tryCatch(max(abs(hind_ssb(s$EM, hind_nyrs) - em_ref)), error = function(e) NA_real_)
      cat(sprintf("  [%s sim=%-5s samp=%-5s %s] OM hindcast SSB max|dev| = %.3e   (EM %.3e, informational)\n",
                  label, flags[1], flags[2], nm, d_om, d_em))
    }
  }
}

# --- single species ----------------------------------------------------------
data(BS2017SS)
BS2017SS$projyr <- 2020
BS2017SS$fleet_control$Proj_F_proportion <- rep(1, 7)
ss <- Rceattle::fit_mod(data_list = BS2017SS, inits = NULL, file = NULL,
  estimateMode = 1, random_rec = FALSE, msmMode = 0,
  fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0))
ss_em <- suppressWarnings(Rceattle::fit_mod(
  data_list = BS2017SS, inits = ss$estimated_params, estimateMode = 2,
  HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2, Alpha = 0.05),
  msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0)))
cat("Single species (BS2017SS):\n")
check_mse("SS", ss, ss_em)

# --- multispecies (the case that needs the suitability pin) ------------------
data(BS2017MS)
BS2017MS$projyr <- 2020
BS2017MS$fleet_control$Proj_F_proportion <- rep(1, nrow(BS2017MS$fleet_control))
ms <- suppressWarnings(Rceattle::fit_mod(
  data_list = BS2017MS, inits = ss$estimated_params, file = NULL,
  estimateMode = 1, random_rec = FALSE, msmMode = 1, suitMode = 0, niter = 5,
  fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
ms_em <- suppressWarnings(Rceattle::fit_mod(
  data_list = BS2017MS, inits = ms$estimated_params, estimateMode = 2,
  HCR = build_hcr(HCR = 2, Ftarget = rep(0.15, ms$data_list$nspp)),
  msmMode = 1, suitMode = 0, niter = 5,
  fit_control = fit_control(getsd = FALSE, verbose = 0)))
cat("Multispecies (BS2017MS):\n")
check_mse("MS", ms, ms_em)
cat("\nDone. OM hindcast max|dev| ~ 0 confirms the hindcast is not perturbed.\n")
