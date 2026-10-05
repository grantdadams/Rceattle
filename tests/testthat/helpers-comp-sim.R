# Shared fixtures for the composition-simulation tests
# (test-functions-sim-mod-comp*.R).
#
# These live in a helper rather than in one test file because that file was
# 1906 s = 31.8 min of a 117.5-minute coverage suite -- 27% of it in one file,
# which set the wall-clock floor no matter how many shards it was split across
# (tools/ci/coverage-costs.tsv). The assertions are unchanged; they are spread
# over several files so the shard planner can put them on different runners.

# One fitted single-species model whose composition switches the caller cares
# about are set before the fit. Sample_size is per row, so it is set on the
# whole comp_data table.
.comp_fixture <- function(nages = 6, N = 200, yng = NULL, old = NULL,
                          dist = NULL, weight = NULL, seed = 123) {
  d <- make_test_data(nyrs = 20, nages = nages, seed = seed)
  d$comp_data$Sample_size <- N
  if (!is.null(yng))    d$fleet_control$Comp_accum_young <- yng
  if (!is.null(old))    d$fleet_control$Comp_accum_old   <- old
  if (!is.null(dist))   d$fleet_control$Comp_distribution <- dist
  if (!is.null(weight)) d$fleet_control$Comp_weights <- weight
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
}

# nrep composition draws from one fit, as a bins x replicates array per row.
.comp_draws <- function(fit, nrep, seed = 9) {
  cols <- grep("^Comp_", names(fit$data_list$comp_data))
  set.seed(seed)
  replicate(nrep, suppressWarnings(
    as.matrix(Rceattle::sim_mod(fit, simulate = TRUE)$comp_data[, cols])))
}
