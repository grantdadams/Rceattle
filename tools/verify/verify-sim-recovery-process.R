# End-to-end self-consistency: draw process error with a known sigma, refit the
# simulated data, and see what comes back. Closes the loop the other harnesses
# leave open -- they check the draw's moments and the density's structure
# separately, which does not show that a refit recovers what was simulated.
#
# GATED on deviation recovery, not on sigma. make_test_data() cannot identify an
# M random-effect variance (sigma goes to ~1e-20 with no simulation involved), so
# a sigma gate would test the fixture rather than the code. The correlation still
# catches the failure this exists for: if the draw and the density disagreed, the
# estimated field would bear no relation to the simulated one.
#
# For sigma on data that do inform M, see verify-sim-recovery-M.R.
#
# Run: export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript tools/verify/verify-sim-recovery-process.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, e)

CTL   <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
NREP  <- 12L
SIGMA <- 0.5
fails <- character(0)

# A longer hindcast than the default fixture: one variance estimated from eight
# deviations is not identified, and a "failure" there would be the fixture.
d <- e$make_test_data(nyrs = 40)

m1_linkage <- function() Rceattle::build_M1(
  M1_model = 1,
  linkages = list(M1 = Rceattle::linkage_spec(~ (1 | Year), by = ~ species)))

cat("Fitting the operating model...\n")
om <- suppressWarnings(Rceattle::fit_mod(
  data_list = d, estimateMode = 1, msmMode = 0, fit_control = CTL,
  M1Fun = m1_linkage()))

# Pin the generating sigma. .sim_obj() hands back this same object, so writing
# last.par.best here is what the draw will use.
om$obj$fn(om$obj$par)
par <- om$obj$env$last.par.best
stopifnot(any(names(par) == "log_sigma_linkage"))
par[names(par) == "log_sigma_linkage"] <- log(SIGMA)
om$obj$env$last.par.best <- par

sig_hat <- rep(NA_real_, NREP)
re_cor  <- rep(NA_real_, NREP)

for (r in seq_len(NREP)) {
  set.seed(1000 + r)
  sim <- suppressWarnings(
    Rceattle::sim_mod(om, simulate = TRUE, process = "M"))
  truth <- attr(sim, "process_sim")$beta_linkage_re
  if (is.null(truth)) stop("no linkage truth returned -- the draw did not happen")

  ref <- try(suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = sim, estimateMode = 1, msmMode = 0, fit_control = CTL,
    M1Fun = m1_linkage()))), silent = TRUE)
  if (inherits(ref, "try-error")) { cat(sprintf("  rep %2d: refit failed\n", r)); next }

  sig_hat[r] <- exp(as.numeric(ref$estimated_params$log_sigma_linkage)[1])
  est <- as.numeric(ref$estimated_params$beta_linkage_re)
  n   <- min(length(est), length(truth))
  re_cor[r] <- suppressWarnings(stats::cor(est[seq_len(n)], as.numeric(truth)[seq_len(n)]))
  cat(sprintf("  rep %2d: sigma_hat %.3f   cor(est, true RE) %+.3f\n",
              r, sig_hat[r], re_cor[r]))
}

ok     <- sum(!is.na(sig_hat))
usable <- sum(!is.na(re_cor))
cat(sprintf("\n%d/%d replicates converged; %d gave a usable correlation\n",
            ok, NREP, usable))
# A replicate whose estimated field came back constant has an undefined
# correlation. Reported rather than quietly dropped by na.rm: if most of them
# went that way the gate below would be averaging over almost nothing.
if (usable < ok) {
  cat(sprintf("  (%d replicate(s) estimated a constant field -- correlation undefined)\n",
              ok - usable))
}
if (usable < 3L) {
  cat("FAIL: too few usable correlations to judge recovery\n"); quit(status = 1)
}
if (ok < 3L) {
  cat("FAIL: too few replicates to judge recovery\n"); quit(status = 1)
}

m_sig <- mean(sig_hat, na.rm = TRUE)
se    <- stats::sd(sig_hat, na.rm = TRUE) / sqrt(ok)
m_cor <- mean(re_cor, na.rm = TRUE)

cat(sprintf("INFORMATIONAL  sigma: true %.3f, mean estimate %.4g (se %.3g)\n",
            SIGMA, m_sig, se))
if (m_sig < 0.1 * SIGMA) {
  cat("               (collapsed toward zero -- this fixture does not identify\n",
      "               an M random-effect variance, with or without simulation)\n", sep = "")
}
cat(sprintf("GATED          mean cor(estimated RE, true RE): %+.3f\n", m_cor))

# The gate. A draw taken from a different distribution than the one being fitted
# would leave the estimated field unrelated to the simulated one, so a clearly
# positive correlation across replicates is the evidence that draw and density
# agree. It is not a precision claim: shrinkage under a collapsing sigma caps
# how high this can go.
if (!is.finite(m_cor) || m_cor < 0.3) {
  fails <- c(fails, sprintf("deviation recovery: mean cor %.3f (expected > 0.3)", m_cor))
}
# Every replicate should point the same way; a sign flip would mean the field
# is being recovered mirrored rather than reproduced.
if (any(re_cor < 0, na.rm = TRUE)) {
  fails <- c(fails, sprintf("%d replicate(s) with negative correlation",
                            sum(re_cor < 0, na.rm = TRUE)))
}

cat("\n")
if (length(fails)) {
  cat("FAIL:\n"); for (f in fails) cat("  - ", f, "\n", sep = ""); quit(status = 1)
}
cat("PASS: a refit recovers the simulated deviation field (draw and density agree)\n")
