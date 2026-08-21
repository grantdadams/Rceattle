# tools/verify/verify-sim-recovery.R
# Which observation type costs a self_test() its recovery? Deterministic
# recovery (simulate = FALSE) is exact, so any loss is in the random draws.
# Re-simulate ONE type at a time, holding the rest at their expected values, and
# refit -- that localises the loss to a data type instead of leaving one number
# for the whole simulator.
#
# NOT a gate. 8 replicates is nowhere near enough to judge recovery from, and
# make_test_data()'s defaults are not a recovery benchmark in the first place:
# 49 free parameters against 20 index points and 2 comp rows at Sample_size = 1
# is under-determined, and the refits legitimately find a BETTER objective than
# the truth (by the expected npar/2). A -50% terminal SSB there is correct
# estimation of an under-determined model, not a simulator bug. Use
# Sample_size = 200 or the plain CV 0.1 fixture when recovery is the thing being
# judged, and verify-sim-local-optimum.R to tell estimation from optimization.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-sim-recovery.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, e)

dat <- e$make_test_data(nyrs = 20, nages = 5, seed = 123)
om <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
  dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
  fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                      phase = FALSE))))
nyrs <- om$data_list$endyr - om$data_list$styr + 1L
truth <- as.numeric(om$quantities$ssb[1, nyrs])
cat("truth terminal SSB:", signif(truth, 6), "\n")
cat("index Log_sd :", unique(om$data_list$index_data$Log_sd), "\n")
cat("catch Log_sd :", unique(om$data_list$catch_data$Log_sd), "\n")
cat("comp Sample_size:", unique(om$data_list$comp_data$Sample_size), "\n\n")

exp_dat <- Rceattle::sim_mod(om, simulate = FALSE)

one <- function(label, mutate, nrep = 8) {
  got <- numeric(0)
  for (k in seq_len(nrep)) {
    set.seed(100 + k)
    sim <- Rceattle::sim_mod(om, simulate = TRUE)   # everything redrawn
    d <- exp_dat                                    # start from expected
    d <- mutate(d, sim)                             # swap in ONE random piece
    f <- tryCatch(suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
      inits = om$estimated_params,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE)))),
      error = function(err) NULL)
    if (!is.null(f)) got <- c(got, as.numeric(f$quantities$ssb[1, nyrs]))
  }
  cat(sprintf("%-34s n=%d  median %8.3f  (%+7.1f%%)  range [%+.1f%%, %+.1f%%]\n",
              label, length(got), stats::median(got),
              100 * (stats::median(got) - truth) / truth,
              100 * (min(got) - truth) / truth,
              100 * (max(got) - truth) / truth))
}

one("index only", function(d, s) {
  d$index_data$Observation <- s$index_data$Observation; d })
one("catch only", function(d, s) {
  d$catch_data$Catch <- s$catch_data$Catch; d })
one("comp only", function(d, s) {
  d$comp_data <- s$comp_data; d })
one("everything (= self_test)", function(d, s) s)
