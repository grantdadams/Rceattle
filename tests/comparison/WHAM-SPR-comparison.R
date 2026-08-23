# Spawning biomass per recruit: Rceattle vs WHAM ------------------------------
#
# Cross-check src/TMB/spr.hpp against an independent implementation. WHAM's
# `get_SPR()` (R/wham_extra.R, unexported) is the reference: it is the routine
# behind WHAM's own F40% search, written by a different author from a different
# spec, so agreement is evidence about the formula rather than about our
# transcription of it.
#
# The two use the same conventions, which is what makes the comparison exact
# rather than approximate:
#
#   * plus group   WHAM does n[A] <- n[A-1]*exp(-Z[A-1]) then n[A]/(1-exp(-Z[A]));
#                  spr.hpp's per_recruit_survivors() is the same expression.
#   * spawn timing WHAM discounts by exp(-fracyrssb * Z[a]); Rceattle by
#                  exp(-Z[a] * spawn_month/12), so fracyrssb = spawn_month/12.
#   * F at age     WHAM takes a scalar F and multiplies by selectivity.
#                  Rceattle builds Flimit_at_age as Flimit times the sum over
#                  fishery fleets of sel * Proj_F_proportion, so that sum is the
#                  selectivity WHAM's scalar F multiplies -- one fleet or three.
#
# The one difference is sex: WHAM is female-only and folds that into its
# spawning weight, while Rceattle carries an explicit sex_ratio at age. It is
# passed into WHAM's maturity argument here.
#
# Not part of test_check() -- it needs a WHAM install (see the other scripts in
# this directory), and get_SPR is unexported, so calling it from tests/testthat
# would put a ::: on another package into R CMD check.
#
# Result, Rceattle 5.12.0 vs WHAM 1.0.7.9000: 48 comparisons -- two datasets,
# two spawn months, three species, four quantities -- worst relative difference
# 2.22e-16, one ulp. GOA2018SS is the one that earns its place: nages differs
# per species (10/21/12), species 3 has THREE fishery fleets, and species 1
# ships a fractional spawn_month of 2.52.

library(Rceattle)
stopifnot(requireNamespace("wham", quietly = TRUE))
wham_spr <- get("get_SPR", asNamespace("wham"))

run_one <- function(ds, spawn_month) {
  data(list = ds); d <- get(ds)
  # Proportions must sum to 1 across each species' fishery fleets.
  for (sp in seq_len(d$nspp)) {
    f <- which(d$fleet_control$Fleet_type == 1 & d$fleet_control$Species == sp)
    d$fleet_control$Proj_F_proportion[f] <- 1 / length(f)
  }
  if (!is.null(spawn_month)) d$spawn_month <- rep(spawn_month, d$nspp)

  # initMode 3 is required: Finit is 0 in every other mode, and SPRFinit is
  # then identically SPR0.
  hcr <- build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2, Alpha = 0.05)
  base <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, HCR = hcr,
    initMode = "FishedNonEquilibrium",
    fit_control = fit_control(getsd = FALSE, verbose = 0))))
  # Three different rates, so the four quantities are distinct.
  inits <- base$estimated_params
  inits$log_Flimit  <- rep(log(0.25), length(inits$log_Flimit))
  inits$log_Ftarget <- rep(log(0.15), length(inits$log_Ftarget))
  inits$log_Finit   <- rep(log(0.10), length(inits$log_Finit))
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = inits, estimateMode = 3, msmMode = 0, HCR = hcr,
    initMode = "FishedNonEquilibrium",
    fit_control = fit_control(getsd = FALSE, verbose = 0))))

  q <- fit$quantities; dl <- fit$obj$env$data
  nyrs_hind <- dl$endyr - dl$styr + 1L
  Finit <- exp(fit$estimated_params$log_Finit)
  frac  <- dl$spawn_month / 12
  worst <- 0

  for (sp in seq_along(dl$nages)) {
    na <- dl$nages[sp]; ages <- seq_len(na)
    wt_idx <- 2 * (sp - 1) + 2                     # C++ wt_idx_ssb = 2*sp + 1, 0-based
    flts <- which(dl$flt_type == 1 & dl$flt_spp + 1L == sp)
    sel <- rowSums(vapply(flts, function(f)
      q$sel_at_age[f, 1, ages, nyrs_hind] * d$fleet_control$Proj_F_proportion[f],
      numeric(na)))

    M_end <- q$M_at_age[sp, 1, ages, nyrs_hind]
    M_one <- q$M_at_age[sp, 1, ages, 1]
    wt_end <- q$weight_hat[wt_idx, 1, ages, nyrs_hind]
    wt_one <- q$weight_hat[wt_idx, 1, ages, 1]
    mat <- dl$maturity[sp, ages] * dl$sex_ratio[sp, ages]   # WHAM has no sex term

    cmp <- c(
      SPR0      = q$SPR0[sp]      / wham_spr(0,              M_end, sel, mat, wt_end, frac[sp]),
      SPRlimit  = q$SPRlimit[sp]  / wham_spr(q$Flimit[sp],   M_end, sel, mat, wt_end, frac[sp]),
      SPRtarget = q$SPRtarget[sp] / wham_spr(q$Ftarget[sp],  M_end, sel, mat, wt_end, frac[sp]),
      # Finit is flat across ages, not selectivity-scaled, and reads year 1.
      SPRFinit  = q$SPRFinit[sp]  / wham_spr(Finit[sp], M_one, rep(1, na), mat, wt_one, frac[sp]))
    worst <- max(worst, abs(cmp - 1))
  }
  cat(sprintf("%-10s nages=%-10s fishery fleets=%-8s spawn=%-14s worst rel diff %.3e\n",
              ds, paste(dl$nages, collapse = "/"),
              paste(tabulate(dl$flt_spp[dl$flt_type == 1] + 1L, dl$nspp), collapse = "/"),
              paste(signif(dl$spawn_month, 3), collapse = "/"), worst))
  worst
}

worst_all <- 0
for (ds in c("BS2017SS", "GOA2018SS")) {
  for (sm in list(NULL, 5)) worst_all <- max(worst_all, run_one(ds, sm))
}

cat(sprintf("\nworst relative difference across all 48 comparisons: %.3e\n", worst_all))
stopifnot(worst_all < 1e-10)
cat("Rceattle and WHAM agree on all four per-recruit quantities.\n")
