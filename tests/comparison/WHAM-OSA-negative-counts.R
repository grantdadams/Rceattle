# Issue #108, point 1 -- negative expected composition counts, reproduced in WHAM.
#
# The report against Rceattle was that `osa_residuals()` sometimes returns a
# negative `predicted`, which cannot happen for an expected bin count. This
# script shows WHAM doing the same thing, on WHAM's own example, unmodified.
#
# The comparison is exact rather than analogous:
#
#   * Both hold composition observations in `obsvec` as COUNTS. WHAM multiplies
#     proportions-at-age by Neff inside prepare_wham_input(); Rceattle stores
#     (proportion + comp_offset) * Sample_size. Support is [0, N] either way.
#   * Both default to TMB::oneStepPredict(method = "oneStepGaussianOffMode").
#   * Rceattle's comp_osa.hpp was ported from WHAM's age_comp_osa.hpp, and the
#     conditional-binomial / beta-binomial decomposition matches (ddirmultinom
#     is byte-identical, and Rceattle's osa_squeeze is TMB's squeeze).
#
# So `OSA.agecomp$mean` in WHAM and `predicted` in Rceattle are the same
# quantity on the same scale, and both go below zero.
#
# Requires wham (tested against 1.0.7.9000).

library(wham)

# ---- 1. Fit WHAM's example 1, untouched -----------------------------------
# NOTE: do.osa = TRUE is silently skipped unless do.sdrep = TRUE -- fit_wham()
# gates on `do.osa & mod$is_sdrep` and only warns. Without it you get a fitted
# model with no $osa at all, which looks like "no residuals" rather than
# "never computed".
dat <- read_asap3_dat(file.path(system.file("extdata", package = "wham"),
                                "ex1_SNEMAYT.dat"))
input <- prepare_wham_input(
  dat, model_name = "ex1_osa",
  selectivity = list(model = rep("logistic", 3),
                     initial_pars = list(c(3, 3), c(3, 3), c(3, 3))),
  NAA_re = list(sigma = "rec", cor = "iid"))

fit <- fit_wham(input, do.osa = TRUE, do.sdrep = TRUE,
                do.retro = FALSE, do.proj = FALSE, MakeADFun.silent = TRUE,
                osa.opts = list(method = "oneStepGaussianOffMode",
                                parallel = FALSE))

# ---- 2. The symptom -------------------------------------------------------
osa <- fit$OSA.agecomp        # observation, nll, grad, mean, residual

sum(osa$mean < 0)             # 17
nrow(osa)                     # 792
min(osa$mean)                 # -3.0457
range(osa$observation)        # counts, not proportions: 0 .. 90.35

# `mean` is an expected count, so every one of those rows is impossible.

# ---- 3. It tracks the bin count, not the sample size -----------------------
# Neff is 100 for the fleet and 50 for the indices throughout, so nothing here
# is a small-sample problem -- the same point the AFSC assessments make, where
# 69 of 353 negative rows occurred at Sample_size > 100.
#
# The negatives concentrate in low-count bins (median observed count 0.95
# against 5.94 over all rows), but a small count is not on its own enough: only
# 3 of the 69 rows with an observed count of exactly 0 go negative, and some
# negative rows carry counts as high as 20.5. That is what a conditional-mean
# artifact looks like -- the mean for one observation depends on the whole
# preceding conditioning sequence, not on that bin's count alone. So read this
# as "sparse bins are where it shows up", not as a threshold rule.
neg <- osa$mean < 0
cat(sprintf("observed count -- negative rows: median %.4f, max %.4f\n",
            median(osa$observation[neg]), max(osa$observation[neg])))
cat(sprintf("observed count -- all rows     : median %.4f, max %.4f\n",
            median(osa$observation), max(osa$observation)))
cat(sprintf("rows with observed count == 0  : %d of %d, of which %d have mean < 0\n",
            sum(osa$observation == 0), nrow(osa),
            sum(osa$observation == 0 & neg)))

# ---- 4. Changing Neff after prepare_wham_input() does nothing --------------
# Worth recording because it looks like it should work. The effective sample
# size is baked into obsvec when the input is built, so a later assignment to
# input$data$catch_Neff leaves the likelihood untouched. To vary it, set the
# ESS in the ASAP object before prepare_wham_input().
input2 <- input
input2$data$catch_Neff[] <- 2
input2$data$index_Neff[] <- 2
stopifnot(identical(input$data$obsvec, input2$data$obsvec))   # unchanged

# ---- What this means for issue #108 ---------------------------------------
# The negative expected counts are a property of the one-step-ahead machinery
# applied to bounded, discrete composition observations, not of Rceattle's
# implementation of it. TMB's Gaussian OSA methods return a conditional mean
# from a continuous approximation, which can fall outside the support when a
# bin holds almost no fish. Switching Rceattle to method = "oneStepGeneric"
# makes it worse, not better, which points the same way.
#
# The residual is affected too, not just the mean.
# The residual standardises the observation against that same mean, so pushing
# the mean below the support biases it positive -- here every one of the 17 rows
# has a positive residual (median 1.33, against 0.00 elsewhere), and Rceattle
# behaves the same way (median 0.66, against -0.04, also 100% positive).
#
# The effect is bounded, though: those rows are 17 of 792, only 1 of the 15
# largest absolute residuals has a negative mean, and dropping them barely moves
# the residual distribution (Shapiro p 1.5e-16 either way -- these residuals are
# strongly non-normal for reasons that have nothing to do with this).
neg <- osa$mean < 0
c(neg = median(osa$residual[neg]), other = median(osa$residual[!neg]))
mean(osa$residual[neg] > 0)     # 1 -- all positive
#
# The Rceattle side is pinned by tests/testthat/test-likelihood-osa-residuals.R,
# "a negative expected composition count warns": `predicted` goes negative at
# N <= 2 and is clean by N = 5, and osa_residuals() warns rather than returning
# the residual silently.
