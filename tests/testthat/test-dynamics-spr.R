# Spawning biomass per recruit, section 6.2 of the template.
#
# Four per-recruit quantities are built from the same survival recursion and
# differ only in the mortality applied and the year the rates are read from:
#
#   SPR0       M only                            terminal hindcast year
#   SPRlimit   M + Flimit_at_age  (at age)       terminal hindcast year
#   SPRtarget  M + Ftarget_at_age (at age)       terminal hindcast year
#   SPRFinit   M + Finit          (flat scalar)  FIRST year
#
# SPRFinit differing on both counts is the part that is easy to lose in a
# refactor, so it is asserted explicitly here rather than inferred.
#
# `test-dynamics-brps.R` already checks SPR0 against a closed form, but on a
# fixture with constant M, weight 1, maturity 1 and spawn_month 0 -- so it
# cannot catch a dropped spawn-timing term, a maturity or sex-ratio
# misindexing, or a wrong weight year. SPRlimit, SPRtarget and SPRFinit had no
# direct assertion at all, and those three drive the reference-point penalties
# and R_init. These tests exist so the section can be moved without silently
# changing what it computes.

# Rebuild SPR in R from the reported rates, with every index written out.
.spr_by_hand <- function(fit, dat) {
  q <- fit$quantities
  d <- fit$obj$env$data
  nyrs_hind <- d$endyr - d$styr + 1L
  term <- nyrs_hind                       # 1-based terminal hindcast year
  Finit <- exp(fit$estimated_params$log_Finit)

  # Flimit/Ftarget at age are selectivity-weighted sums over fishery fleets,
  # the way section 5.12 builds them.
  f_at_age <- function(sp, yr, scale) {
    out <- numeric(d$nages[sp])
    for (flt in seq_along(d$flt_type)) {
      if (d$flt_type[flt] != 1) next
      if (d$flt_spp[flt] + 1L != sp) next
      prop <- dat$fleet_control$Proj_F_proportion[flt]
      out <- out + q$sel_at_age[flt, 1, seq_len(d$nages[sp]), yr] * prop * scale[sp]
    }
    out
  }

  out <- list()
  for (sp in seq_along(d$nages)) {
    na <- d$nages[sp]
    wt_idx <- 2 * (sp - 1) + 2          # C++ wt_idx_ssb = 2*sp + 1, 0-based

    M_term <- q$M_at_age[sp, 1, seq_len(na), term]
    M_yr1  <- q$M_at_age[sp, 1, seq_len(na), 1]
    Z <- list(
      SPR0      = M_term,
      SPRlimit  = M_term + f_at_age(sp, term, q$Flimit),
      SPRtarget = M_term + f_at_age(sp, term, q$Ftarget),
      SPRFinit  = M_yr1  + Finit[sp]
    )
    wt <- list(
      SPR0      = q$weight_hat[wt_idx, 1, seq_len(na), term],
      SPRlimit  = q$weight_hat[wt_idx, 1, seq_len(na), term],
      SPRtarget = q$weight_hat[wt_idx, 1, seq_len(na), term],
      SPRFinit  = q$weight_hat[wt_idx, 1, seq_len(na), 1]
    )

    for (nm in names(Z)) {
      z <- Z[[nm]]
      n <- numeric(na)
      n[1] <- 1
      for (a in 2:(na - 1)) n[a] <- n[a - 1] * exp(-z[a - 1])
      # Plus group: survivors into it, then the infinite tail at its own Z.
      n[na] <- n[na - 1] * exp(-z[na - 1]) / (1 - exp(-z[na]))
      out[[nm]] <- c(out[[nm]],
                     sum(n * wt[[nm]] * d$maturity[sp, seq_len(na)] *
                           d$sex_ratio[sp, seq_len(na)] *
                           exp(-z * d$spawn_month[sp] / 12)))
    }
  }
  out
}

testthat::test_that("all four SPRs match an independent recomputation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  d <- BS2017SS
  d$fleet_control$Proj_F_proportion <- rep(1, nrow(d$fleet_control))

  # initMode "FishedNonEquilibrium" (3) is required: section 5.4 sets Finit to 0
  # for every other mode, under which SPRFinit is identically SPR0. That is why
  # a broken SPRFinit could go unnoticed -- no default configuration separates
  # the two.
  hcr  <- Rceattle::build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35,
                              Plimit = 0.2, Alpha = 0.05)
  base <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0, HCR = hcr,
                      initMode = "FishedNonEquilibrium",
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))

  # Distinct rates, so the four quantities cannot pass by coinciding. The
  # shipped starting values are Flimit = Ftarget = 1 and Finit = exp(-10),
  # under which SPRlimit == SPRtarget and SPRFinit is SPR0 to seven decimals.
  inits <- base$estimated_params
  inits$log_Flimit  <- rep(log(0.25), length(inits$log_Flimit))
  inits$log_Ftarget <- rep(log(0.15), length(inits$log_Ftarget))
  inits$log_Finit   <- rep(log(0.10), length(inits$log_Finit))

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, inits = inits, estimateMode = 3, msmMode = 0,
                      HCR = hcr, initMode = "FishedNonEquilibrium",
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))

  want <- .spr_by_hand(fit, d)
  for (nm in c("SPR0", "SPRlimit", "SPRtarget", "SPRFinit")) {
    testthat::expect_equal(as.numeric(fit$quantities[[nm]]), want[[nm]],
                           tolerance = 1e-8, info = nm)
  }

  # The four really are distinct under these rates, so the loop above compared
  # four different numbers rather than one repeated.
  got <- vapply(c("SPR0", "SPRlimit", "SPRtarget", "SPRFinit"),
                function(nm) as.numeric(fit$quantities[[nm]])[1], numeric(1))
  testthat::expect_equal(length(unique(signif(got, 8))), 4L)
  # More fishing removes spawning output: SPR0 is the largest.
  testthat::expect_true(all(got[-1] < got[["SPR0"]]))
})

testthat::test_that("SPR carries spawn timing, maturity and the weight year", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # A closed form on a controlled fixture: constant M, so the recursion is
  # writable by hand, but a non-zero spawn month and a real maturity ogive, so
  # dropping either term fails. `test-dynamics-brps.R` sets both to their
  # neutral values and cannot see them.
  data("GOA2018SS")
  g <- GOA2018SS
  g$M1_base[, 3:32] <- 0.2
  g$weight[, 6:35]  <- 1
  g$sex_ratio[, -1] <- 0.5
  g$spawn_month     <- rep(4, g$nspp)      # spawning a third of the way through
  mat <- c(0, 0.1, 0.5, 0.9, rep(1, 30))   # knife-ish ogive, not all ones
  for (sp in seq_len(g$nspp)) g$maturity[sp, -1] <- mat[seq_len(ncol(g$maturity) - 1)]

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = g, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))

  M <- 0.2
  for (sp in seq_len(g$nspp)) {
    na <- fit$obj$env$data$nages[sp]
    n <- numeric(na)
    n[1] <- 1
    for (a in 2:(na - 1)) n[a] <- n[a - 1] * exp(-M)
    n[na] <- n[na - 1] * exp(-M) / (1 - exp(-M))
    want <- sum(n * 1 * mat[seq_len(na)] * 0.5 * exp(-M * 4 / 12))
    testthat::expect_equal(as.numeric(fit$quantities$SPR0)[sp], want,
                           tolerance = 1e-6, info = paste("species", sp))
  }
})

testthat::test_that("NbyageSPR is zero where it is never filled", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Two holes the recursion never writes: ages past a species' own nages when
  # another species is longer-lived, and every element under multispecies mode,
  # where section 6.2 does not run at all. Both are reported, so both are read.
  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  nages <- ss$obj$env$data$nages
  max_age <- dim(ss$quantities$NbyageSPR)[3]
  testthat::expect_gt(max_age, min(nages))          # the ragged case exists here
  for (sp in seq_along(nages)) {
    if (nages[sp] >= max_age) next
    tail_ages <- (nages[sp] + 1):max_age
    testthat::expect_true(all(ss$quantities$NbyageSPR[, sp, tail_ages] == 0),
                          info = paste("species", sp))
  }

  data("BS2017MS")
  ms <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017MS, estimateMode = 3, msmMode = 1,
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  testthat::expect_true(all(ms$quantities$NbyageSPR == 0))
  testthat::expect_true(all(ms$quantities$SPR0 == 0))
})

testthat::test_that("a stock-recruit curve that needs SPR is refused under predation", {
  testthat::skip_if_not_installed("TMB")

  # Section 6.2 does not compute SPR under predation, on purpose: total
  # mortality carries M2, so per-recruit spawning output is not a property of
  # the prey stock alone. Section 6.3 read it anyway.
  #
  #   srr_fun >= 2                        SPR0 = 0, so R0 = (alpha - 1/0)/beta
  #                                       is -Inf and the objective is NaN.
  #   srr_fun < 2, srr_pred_fun >= 2      the Ianelli configuration, and the
  #                                       dangerous one: the fit RUNS, with
  #                                       steepness 0, R_hat -Inf and a finite
  #                                       objective.
  data("BS2017MS")

  for (cfg in list(list(f = "BevertonHolt", p = "BevertonHolt"),
                   list(f = "Ricker",       p = "Ricker"),
                   list(f = "mean",         p = "BevertonHolt"))) {
    srr <- suppressWarnings(Rceattle::build_srr(
      srr_fun = cfg$f, srr_pred_fun = cfg$p, proj_mean_rec = TRUE,
      srr_est_mode = "Fixed", srr_prior = 0.8))
    testthat::expect_error(
      suppressMessages(suppressWarnings(
        Rceattle::fit_mod(data_list = BS2017MS, estimateMode = 3, msmMode = 1,
                          recFun = srr,
                          fit_control = Rceattle::fit_control(verbose = 0)))),
      "spawning biomass per recruit",
      info = paste(cfg$f, cfg$p))
  }

  # Mean recruitment is unaffected, and so is the same curve single-species.
  data("BS2017SS")
  testthat::expect_no_error(suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017MS, estimateMode = 3, msmMode = 1,
                      fit_control = Rceattle::fit_control(verbose = 0)))))

  # srr_fun can still be the canonical string here. fit_mod() coerces it --
  # with recFun = NULL it rebuilds the default and the string never survives --
  # but data_check() is documented as the way to validate a data_list on its
  # own, so it has to read the string form correctly. This is the case that
  # matters: R compares "mean" >= 2 as characters and calls it TRUE, so a guard
  # reading the raw value refuses a good mean-recruitment multispecies model.
  # switch_check() first, as fit_mod() does -- it normalises the fleet_control
  # switches but leaves srr_fun alone, which is exactly how a string reaches
  # data_check().
  d_str <- suppressMessages(suppressWarnings(switch_check(BS2017MS)))
  d_str$msmMode <- 1
  d_str$srr_fun <- d_str$srr_pred_fun <- "mean"
  testthat::expect_no_error(suppressMessages(suppressWarnings(data_check(d_str))))

  d_bh <- d_str
  d_bh$srr_fun <- d_bh$srr_pred_fun <- "BevertonHolt"
  testthat::expect_error(suppressMessages(suppressWarnings(data_check(d_bh))),
                         "spawning biomass per recruit")
  srr_bh <- suppressWarnings(Rceattle::build_srr(
    srr_fun = "BevertonHolt", srr_pred_fun = "BevertonHolt", proj_mean_rec = TRUE,
    srr_est_mode = "Fixed", srr_prior = 0.8))
  testthat::expect_no_error(suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      recFun = srr_bh,
                      fit_control = Rceattle::fit_control(verbose = 0)))))
})
