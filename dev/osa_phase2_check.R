# Phase 2 OSA validation (composition: comp + caal). Run via:
#   Rscript R/dev/osa_phase2_check.R
# Requires a baseline src/TMB/ceattle_base.cpp (from the dev branch) for the
# regression comparison.

suppressMessages(library(devtools))
devtools::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "helpers.R"))

ok <- TRUE
check <- function(label, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) ok <<- FALSE
}

# ---- 1. Regression: osa_mode = 0 jnll_comp identical to baseline ----------
# make_test_data carries comp data, so this exercises the comp/caal fitting
# path (which is the original code wrapped in `if (osa_mode == 0)`).
cat("== 1. baseline vs new jnll_comp (with composition data) ==\n")
dat <- make_test_data(nyrs = 8, nages = 5, seed = 123)
fit_new <- Rceattle::fit_mod(dat, inits = NULL, estimateMode = 3,
                             fit_control = fit_control(phase = FALSE, verbose = 0))
par <- fit_new$obj$env$last.par.best
rep_new <- fit_new$obj$report(par)

if (file.exists(file.path("src", "TMB", "ceattle_base.cpp"))) {
  fit_base <- Rceattle::fit_mod(
    dat, inits = NULL, estimateMode = 3,
    fit_control = fit_control(phase = FALSE, verbose = 0,
                              TMBfilename = file.path("src", "TMB", "ceattle_base")))
  rep_base <- fit_base$obj$report(par)
  d <- max(abs(rep_new$jnll_comp - rep_base$jnll_comp))
  cat(sprintf("   max |jnll_comp diff| = %.3e\n", d))
  cat(sprintf("   comp slot (2): new=%.8f base=%.8f\n",
              sum(rep_new$jnll_comp[3, ]), sum(rep_base$jnll_comp[3, ])))
  check("jnll_comp identical to baseline (<1e-8)", d < 1e-8)
} else {
  cat("   (no baseline; skipping)\n")
}

# ---- 2. Composition OSA end-to-end on GOApollock --------------------------
cat("\n== 2. comp OSA on GOApollock ==\n")
data("GOApollock")
fit <- Rceattle::fit_mod(GOApollock, inits = NULL, file = NULL, estimateMode = 1,
                         random_rec = FALSE, msmMode = 0,
                         fit_control = fit_control(verbose = 0, phase = TRUE))
cat("   obs_ctl types:\n"); print(table(fit$obs_ctl$type))

osa <- osa_residuals(fit, source = c("index", "catch", "comp"))
cat("   total residuals:", nrow(osa),
    " comp residuals:", sum(osa$type == "comp"), "\n")
check("comp residuals present", sum(osa$type == "comp") > 0)
check("comp residuals finite", all(is.finite(osa$residual[osa$type == "comp"])))
check("aggregate residuals still finite",
      all(is.finite(osa$residual[osa$type %in% c("index", "catch")])))

osa2 <- osa_residuals(fit, source = c("index", "catch", "comp"))
check("deterministic re-run", isTRUE(all.equal(osa$residual, osa2$residual)))

cat("   diagnostics (head):\n")
print(utils::head(osa_diagnostics(osa), 12))

cat(sprintf("\n==== OVERALL: %s ====\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
