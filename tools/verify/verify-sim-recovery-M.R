# Can the M random-effect variance be recovered from simulated process error, on
# data that actually inform M? verify-sim-recovery-process.R asks the weaker
# question (does the deviation FIELD come back) on a fixture that cannot identify
# a variance at all; this uses BS2017SS, a real assessment.
#
# A diagnostic, not a gate: a negative answer describes what these data support,
# not a defect, and should be reported rather than made to pass.
#
# RECORDED RESULT (BS2017SS, M1_re = 2, sigma = 0.3, 8 replicates):
#   deviations  mean cor(estimated, true) = +0.52, positive every replicate
#   sigma       true 0.300, mean estimate 0.089 -- 70% low, 3 of 8 collapsed to 0
#
# The field is recovered; the variance is not. That is estimation, not the draw:
# the draw's marginal SD is exact to four figures in verify-sim-process-error.R
# and the density orientation is pinned against an MVN reference. A variance
# component with little information is biased low and piles up at zero, and a
# year-varying M trades off against selectivity and catchability.
#
# So read self_test(process = "M") as testing deviation recovery, not variance
# recovery -- a recovered sigma is not evidence it works, a collapsed one not
# evidence it is broken.
#
# Run: export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript tools/verify/verify-sim-recovery-M.R [nrep]

suppressMessages(pkgload::load_all(".", quiet = TRUE))

args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args)) as.integer(args[1]) else 8L
SIGMA <- 0.3          # generating SD of the M deviations
CTL   <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)

# M1_re = 2: one deviation per YEAR, shared across ages. The smallest random
# structure that still has a variance to estimate -- an age x year field on
# these data is confounded with selectivity, which would confound the question
# being asked with a separate identifiability problem.
m_fun <- function() Rceattle::build_M1(M1_model = 1, M1_re = 2)

cat("Fitting the operating model on BS2017SS (M1_re = 2)...\n")
om <- suppressWarnings(Rceattle::fit_mod(
  data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
  estimateMode = 1, random_rec = FALSE, msmMode = 0,
  M1Fun = m_fun(), fit_control = CTL))
cat(sprintf("  objective %.4f\n", om$opt$objective))

# Pin the generating sigma. M1_dev_log_sd is the INNOVATION sd; with rho left at
# 0 the marginal equals it, so SIGMA is what the draws should show.
om$obj$fn(om$obj$par)
par <- om$obj$env$last.par.best
stopifnot(any(names(par) == "M1_dev_log_sd"))
par[names(par) == "M1_dev_log_sd"] <- log(SIGMA)
par[names(par) == "M1_rho"] <- 0
om$obj$env$last.par.best <- par

sig_hat <- rep(NA_real_, NREP)
dev_cor <- rep(NA_real_, NREP)

for (r in seq_len(NREP)) {
  set.seed(20000 + r)
  sim <- suppressWarnings(Rceattle::sim_mod(om, simulate = TRUE, process = "M"))
  truth <- attr(sim, "process_sim")$log_M1_dev
  if (is.null(truth)) stop("no M truth returned -- the draw did not happen")

  ref <- try(suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = sim, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0, M1Fun = m_fun(), fit_control = CTL))),
    silent = TRUE)
  if (inherits(ref, "try-error")) { cat(sprintf("  rep %2d: refit failed\n", r)); next }

  sig_hat[r] <- exp(as.numeric(ref$estimated_params$M1_dev_log_sd)[1])
  # Compare the year-varying deviation actually used: species 1, sex 1, age 1,
  # over the hindcast (modes 2/5 hold the deviation constant across ages).
  nyr <- om$data_list$endyr - om$data_list$styr + 1L
  est <- as.numeric(ref$estimated_params$log_M1_dev[1, 1, 1, seq_len(nyr)])
  tru <- as.numeric(truth[1, 1, 1, seq_len(nyr)])
  dev_cor[r] <- suppressWarnings(stats::cor(est, tru))
  cat(sprintf("  rep %2d: sigma_hat %.4f   cor(est, true dev) %+.3f\n",
              r, sig_hat[r], dev_cor[r]))
}

ok     <- sum(!is.na(sig_hat))
usable <- sum(!is.na(dev_cor))
cat(sprintf("\n%d/%d refits converged; %d gave a usable correlation\n",
            ok, NREP, usable))
if (ok < 3L) { cat("INCONCLUSIVE: too few refits converged\n"); quit(status = 0) }

m_sig <- mean(sig_hat, na.rm = TRUE)
se    <- stats::sd(sig_hat, na.rm = TRUE) / sqrt(ok)
m_cor <- mean(dev_cor, na.rm = TRUE)

cat(sprintf("\nsigma:  true %.3f   mean estimate %.4f   se %.4f\n", SIGMA, m_sig, se))
cat(sprintf("        relative bias %+.1f%%\n", 100 * (m_sig - SIGMA) / SIGMA))
cat(sprintf("deviations: mean cor(estimated, true) %+.3f\n", m_cor))

cat("\nREAD THIS AS:\n")
if (m_sig < 0.25 * SIGMA) {
  cat("  sigma COLLAPSES on these data -- the M random-effect variance is not\n")
  cat("  identified even with real index and composition data. A self_test()\n")
  cat("  that redraws M measures deviation recovery only; do not read a\n")
  cat("  recovered sigma from it.\n")
} else if (abs(m_sig - SIGMA) <= max(3 * se, 0.2 * SIGMA)) {
  cat("  sigma is RECOVERED within sampling error. Redrawing M in a self-test\n")
  cat("  is meaningful on data of this quality.\n")
} else {
  cat("  sigma is estimable but BIASED at this sample size. Worth understanding\n")
  cat("  before reading a self_test() that redraws M -- check whether the bias\n")
  cat("  shrinks as the series lengthens before suspecting the simulator.\n")
}
cat(sprintf("  deviation recovery cor = %+.3f (0 would mean the draw and the\n", m_cor))
cat("  density disagree, which the other harnesses would also have caught).\n")
