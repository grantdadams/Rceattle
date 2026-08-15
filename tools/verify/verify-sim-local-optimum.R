# tools/verify/verify-sim-local-optimum.R
# When a self_test() recovers badly, is that basin the MLE for the simulated
# data, or did the optimizer miss? Evaluate the operating model's TRUE
# parameters on the SAME simulated data (estimateMode = 3 builds the real
# objective without optimizing) and compare objectives:
#
#   refit obj < truth obj -> the data really do support it (estimation; nothing
#                            to fix in the simulator or the optimizer)
#   refit obj > truth obj -> the optimizer missed the optimum (actionable, e.g.
#                            phasing or starting values)
#
# Diagnostic, not a gate. Reach for it before concluding a simulator change
# moved recovery.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-sim-local-optimum.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, e)

d <- e$make_test_data(nyrs = 20, nages = 5, seed = 123)
d$index_data$Log_sd[d$index_data$Fleet_name == "Survey"] <- 0.2
om <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
  d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
  fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                      phase = FALSE))))
nyrs <- om$data_list$endyr - om$data_list$styr + 1L
truth <- as.numeric(om$quantities$ssb[1, nyrs])
cat("truth terminal SSB:", signif(truth, 6), "\n\n")

cat(sprintf("%4s %10s %12s %12s %10s\n",
            "sim", "SSB %err", "refit obj", "truth obj", "verdict"))
for (k in 1:6) {
  set.seed(123 + k)
  sim <- Rceattle::sim_mod(om, simulate = TRUE)

  ## (a) refit as self_test does
  f <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    sim, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    inits = om$initial_params,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))

  ## (b) the OM's true parameters evaluated on the SAME simulated data
  ##     (estimateMode = 3 builds the objective without optimizing)
  t3 <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    sim, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    inits = om$estimated_params,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  truth_obj <- t3$obj$fn(t3$obj$par)

  got <- as.numeric(f$quantities$ssb[1, nyrs])
  verdict <- if (f$opt$objective < truth_obj - 1e-6) "refit BETTER" else
    if (f$opt$objective > truth_obj + 1e-6) "refit WORSE" else "tie"
  cat(sprintf("%4d %9.1f%% %12.4f %12.4f  %s\n", k,
              100 * (got - truth) / truth, f$opt$objective, truth_obj, verdict))
}

## Does phasing rescue it?
cat("\n--- same sims, phase = TRUE ---\n")
for (k in 1:6) {
  set.seed(123 + k)
  sim <- Rceattle::sim_mod(om, simulate = TRUE)
  f <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    sim, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    inits = om$initial_params,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = TRUE))))
  got <- as.numeric(f$quantities$ssb[1, nyrs])
  cat(sprintf("  sim %d  SSB %+7.1f%%  obj %11.4f\n", k,
              100 * (got - truth) / truth, f$opt$objective))
}
