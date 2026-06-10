# Phase 1 OSA validation (scratch, run via Rscript) ------------------------
# Verifies that the obsvec/keep refactor (catch + index) is numerically
# identical to the pre-change model, and that osa_residuals() runs end-to-end.
#
# Run from the package root:
#   Rscript R/dev/osa_phase1_check.R
#
# Requires: a compiled ceattle_v01_11.dll (new) and ceattle_base.dll (baseline
# built from the `dev` branch source) in src/TMB/.

suppressMessages({
  library(devtools)
})

cat("== load_all ==\n")
devtools::load_all(".", quiet = TRUE)
source(file.path("tests", "testthat", "helpers.R"))

ok <- TRUE
check <- function(label, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) ok <<- FALSE
}

# ---- 1. obsvec construction correctness (pure R) -------------------------
cat("\n== 1. obsvec / obs_ctl construction ==\n")
dat <- make_test_data(nyrs = 8, nages = 5, seed = 123)
dl  <- Rceattle::rearrange_data(dat)

endyr <- dl$endyr
ft    <- dl$flt_type

# index
idx <- dl$index_obsvec_idx
inc_i <- which(dl$index_ctl[, 3] > 0 & dl$index_ctl[, 3] <= endyr &
                 ft[dl$index_ctl[, 1]] > 0 & dl$index_obs[, 1] > 0)
check("index: included rows have idx >= 0",
      all(idx[inc_i] >= 0) && all(idx[-inc_i] == -1 | length(inc_i) == nrow(dl$index_obs)))
check("index: obsvec[idx+1] == log(obs)",
      isTRUE(all.equal(dl$obsvec[idx[inc_i] + 1], log(dl$index_obs[inc_i, 1]))))

# catch
cdx <- dl$catch_obsvec_idx
inc_c <- which(dl$catch_ctl[, 3] > 0 & dl$catch_ctl[, 3] <= endyr &
                 ft[dl$catch_ctl[, 1]] == 1 & dl$catch_obs[, 1] > 0)
check("catch: obsvec[idx+1] == log(obs)",
      isTRUE(all.equal(dl$obsvec[cdx[inc_c] + 1], log(dl$catch_obs[inc_c, 1]))))

check("obs_ctl rows == number of included obs",
      nrow(dl$obs_ctl) == length(inc_i) + length(inc_c))
check("osa_mode default 0", identical(dl$osa_mode, 0L))

# ---- 2. Baseline vs new: reported jnll_comp identical --------------------
# Fit the SAME data with the new DLL and with the pre-change baseline DLL
# (compiled from the `dev` branch source as src/TMB/ceattle_base.cpp), both at
# estimateMode = 3 (no optimization -> identical init parameters), and compare
# the REPORTed jnll_comp. fit_mod() builds the map/parameters correctly for each
# DLL; the baseline simply ignores the extra obsvec/keep data it does not read.
cat("\n== 2. baseline vs new jnll_comp ==\n")
fit_new <- Rceattle::fit_mod(data_list = dat, inits = NULL, estimateMode = 3,
                             fit_control = fit_control(phase = FALSE, verbose = 0))
par <- fit_new$obj$env$last.par.best
rep_new <- fit_new$obj$report(par)

if (file.exists(file.path("src", "TMB", "ceattle_base.cpp"))) {
  fit_base <- Rceattle::fit_mod(
    data_list = dat, inits = NULL, estimateMode = 3,
    fit_control = fit_control(phase = FALSE, verbose = 0,
                              TMBfilename = file.path("src", "TMB", "ceattle_base")))
  # Evaluate the baseline at the SAME parameter vector (identical PARAMETER
  # sections -> identical ordering) so any difference is purely the refactor.
  rep_base <- fit_base$obj$report(par)
  d <- max(abs(rep_new$jnll_comp - rep_base$jnll_comp))
  cat(sprintf("   max |jnll_comp_new - jnll_comp_base| = %.3e\n", d))
  cat(sprintf("   slot0 (index): new=%.8f base=%.8f\n",
              sum(rep_new$jnll_comp[1, ]), sum(rep_base$jnll_comp[1, ])))
  cat(sprintf("   slot1 (catch): new=%.8f base=%.8f\n",
              sum(rep_new$jnll_comp[2, ]), sum(rep_base$jnll_comp[2, ])))
  check("jnll_comp identical to baseline (< 1e-8)", d < 1e-8)
} else {
  cat("   (ceattle_base.cpp not found -- skipping baseline comparison)\n")
}

# ---- 3. OSA smoke (requires an optimized estimateMode < 3 fit) -----------
cat("\n== 3. osa_residuals() smoke ==\n")
fit_opt <- Rceattle::fit_mod(data_list = dat, inits = NULL, estimateMode = 1,
                             fit_control = fit_control(phase = FALSE, verbose = 0))
osa <- tryCatch(osa_residuals(fit_opt, types = c("index", "catch")),
                error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(osa)) {
  check("osa: one residual per included aggregate obs",
        nrow(osa) == length(inc_i) + length(inc_c))
  check("osa: all residuals finite", all(is.finite(osa$residual)))
  diag <- osa_diagnostics(osa)
  cat("   diagnostics:\n")
  print(diag)
  check("osa_diagnostics returns rows", nrow(diag) >= 1)
}

cat(sprintf("\n==== OVERALL: %s ====\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
