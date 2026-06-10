# Phase 4 process-residuals check. Run via: Rscript R/dev/osa_phase4_check.R
suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))

ok <- TRUE
check <- function(label, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) ok <<- FALSE
}

data("GOApollock")
# random_rec = FALSE: rec_dev is a penalized fixed effect (same N(R_sd^2/2,
# R_sd^2) prior) and the model converges cleanly, so process residuals come from
# the fixed-effect covariance. (random_rec = TRUE would use the joint precision.)
fit <- Rceattle::fit_mod(
  data_list = GOApollock, inits = NULL, file = NULL, estimateMode = 1,
  random_rec = FALSE, msmMode = 0,
  fit_control = fit_control(verbose = 0, phase = TRUE, getsd = TRUE))

mg <- tryCatch(max(abs(fit$obj$gr(fit$obj$env$last.par.best))),
               error = function(e) NA)
cat("obj:", fit$opt$objective, " maxgrad:", mg, "\n")

pr <- process_residuals(fit, process = "recruitment")
cat("recruitment process residuals:", nrow(pr), "\n")
check("one residual per hindcast year",
      nrow(pr) == length(GOApollock$styr:GOApollock$endyr))
check("residuals finite", all(is.finite(pr$residual)))
cat(sprintf("mean = %.3f  sd = %.3f  (expect ~0, ~1 for a good process)\n",
            mean(pr$residual), stats::sd(pr$residual)))
cat("year range:", range(pr$year), "\n")
print(utils::head(pr))

pr2 <- process_residuals(fit, process = "recruitment")
check("deterministic (same seed)", isTRUE(all.equal(pr$residual, pr2$residual)))

cat("diagnostics:\n"); print(osa_diagnostics(pr))

r <- residuals(fit, type = "process")
check("residuals(type='process') schema",
      all(c("Source", "Species", "Year", "Residual") %in% names(r)))

cat("\n-- process = 'all' --\n")
pr_all <- process_residuals(fit, process = "all")
cat("processes present:", paste(unique(pr_all$type), collapse = ", "),
    " (n =", nrow(pr_all), ")\n")
check("'all' returns >= recruitment", "recruitment" %in% pr_all$type)
check("'all' residuals finite", all(is.finite(pr_all$residual)))

cat(sprintf("\n==== OVERALL: %s ====\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
