# dev/verify-mse-om-horizon.R
# Equivalence + timing harness for refitting the MSE operating model on a
# shortened projection horizon (R/10-run_mse.R).
#
# Two scenarios, because they gate different things:
#
#   sim=FALSE  sim_mod() draws no random numbers at all (it copies the expected
#              values), so the shortened-horizon run must be BIT-IDENTICAL to
#              the full-horizon one. This is the real gate on the truncate /
#              restore logic: parameters, data tables and the catch quantities
#              all have to come back intact.
#
#   sim=TRUE   sim_mod() draws per observation row, and the operating model
#              carries fewer rows on a shortened horizon, so the random stream
#              advances differently. Results legitimately differ here; the digest
#              is recorded to show HOW MUCH, not to require equality.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript dev/verify-mse-om-horizon.R dev/mse-horizon-before.rds
#   NOT_CRAN=true Rscript dev/verify-mse-om-horizon.R dev/mse-horizon-after.rds compare dev/mse-horizon-before.rds

# Fourth argument selects the operating model: "ss" (single-species, cheap) or
# "ms" (multispecies, where the tape cost per model year is large enough for the
# shortened horizon to matter).
args         <- commandArgs(trailingOnly = TRUE)
out_path     <- if (length(args) >= 1) args[1] else "dev/mse-horizon-digest.rds"
do_compare   <- length(args) >= 2 && args[2] == "compare"
compare_path <- if (length(args) >= 3) args[3] else NULL
model        <- if (length(args) >= 4) args[4] else "ss"

suppressMessages(pkgload::load_all(".", quiet = TRUE))
data(BS2017SS); data(BS2017MS)

BS2017SS$projyr <- 2040
BS2017SS$fleet_control$Proj_F_proportion <- rep(1, 7)

om_ss <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
  data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
  random_rec = FALSE, msmMode = 0,
  fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))

if (identical(model, "ms")) {
  BS2017MS$projyr <- 2040
  BS2017MS$fleet_control$Proj_F_proportion <- rep(1, 7)
  # Warm-started from the single-species MLEs, as everywhere else in the repo:
  # the predation likelihood is non-convex, so the start point picks the optimum.
  om <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = BS2017MS, inits = om_ss$estimated_params, file = NULL,
    estimateMode = 1, niter = 5, random_rec = FALSE, msmMode = 1, suitMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0))))
} else {
  om <- om_ss
}

em <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
  data_list = BS2017SS, inits = om_ss$estimated_params, estimateMode = 2,
  HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2,
                  Alpha = 0.05),
  msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0))))

# Digest: the operating model trajectory plus the catch advice the estimation
# models produced, which together determine everything an MSE reports.
digest_one <- function(simulate) {
  t0 <- proc.time()[["elapsed"]]
  mse <- suppressWarnings(suppressMessages(Rceattle::run_mse(
    om = om, em = em, nsim = 1, assessment_period = 1, sampling_period = 1,
    simulate_data = simulate, sample_rec = simulate, seed = 42)))
  sec <- proc.time()[["elapsed"]] - t0
  sim <- mse$Sim_1
  list(
    sec       = sec,
    om_ssb    = sim$OM$quantities$ssb,
    om_R      = sim$OM$quantities$R,
    om_catch  = sim$OM$data_list$catch_data$Catch,
    em_ssb    = lapply(sim$EM, function(x) x$quantities$ssb),
    em_ftarget = lapply(sim$EM, function(x) x$quantities$Ftarget),
    use_sim   = sim$use_sim)
}

digest <- list(sim_off = digest_one(FALSE), sim_on = digest_one(TRUE))

saveRDS(digest, out_path)
cat("wrote digest ->", out_path, "\n")
for (nm in names(digest)) {
  cat(sprintf("  %-8s run_mse %6.1fs  use_sim=%s\n",
              nm, digest[[nm]]$sec, digest[[nm]]$use_sim))
}

if (do_compare && !is.null(compare_path) && file.exists(compare_path)) {
  ref <- readRDS(compare_path)
  cat("\n=== COMPARE vs", compare_path, "===\n")

  # Timing is expected to move; drop it before testing equality.
  strip <- function(x) { x$sec <- NULL; x }

  off_ok <- identical(strip(digest$sim_off), strip(ref$sim_off))
  cat(if (off_ok) "PASS" else "FAIL",
      ": sim=FALSE bit-identical (no RNG consumed; must match exactly)\n", sep = "")
  if (!off_ok) {
    for (f in names(strip(digest$sim_off))) {
      if (!identical(digest$sim_off[[f]], ref$sim_off[[f]])) cat("   DIFF:", f, "\n")
    }
  }

  # sim=TRUE: report the magnitude rather than demanding equality.
  d <- suppressWarnings(max(abs(as.numeric(digest$sim_on$om_ssb) -
                                as.numeric(ref$sim_on$om_ssb)), na.rm = TRUE))
  rel <- d / max(abs(as.numeric(ref$sim_on$om_ssb)), na.rm = TRUE)
  cat(sprintf("INFO: sim=TRUE OM SSB max|dev| = %.4g (%.2f%% of peak SSB) -- ",
              d, 100 * rel),
      "expected: the observation-error stream differs.\n", sep = "")

  cat(sprintf("\nSpeed: sim=FALSE %6.1fs -> %6.1fs (%+.1f%%);  sim=TRUE %6.1fs -> %6.1fs (%+.1f%%)\n",
              ref$sim_off$sec, digest$sim_off$sec,
              100 * (digest$sim_off$sec / ref$sim_off$sec - 1),
              ref$sim_on$sec, digest$sim_on$sec,
              100 * (digest$sim_on$sec / ref$sim_on$sec - 1)))
}
