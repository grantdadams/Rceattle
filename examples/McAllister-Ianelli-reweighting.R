
# Composition data weighting by the McAllister-Ianelli method.
#
# The weight multiplies a fleet's input sample size, and the weight the fit
# implies depends on the fit itself, so it has to be tuned iteratively: refit,
# read the implied weights, refit again, until they settle.

# Load
library(Rceattle)

# Run the 2017 single species assessment for the Bering Sea
data(BS2017SS)
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            fit_control = fit_control(
                              phase = TRUE,
                              verbose = 1))


# reweight_comps() runs the tuning loop: it refits the model the way it was
# originally fitted, updating Comp_weights each pass until every fleet's weight
# moves less than `tol` or `n_iter` is reached.
ss_run_reweighted <- Rceattle::reweight_comps(ss_run, n_iter = 10, tol = 0.01)

# One row per fleet per iteration; the last rows are the weights the returned
# model was fitted with.
ss_run_reweighted$reweight$history
ss_run_reweighted$reweight$converged

# Tune only selected fleets. `fleets` takes Fleet_code values, not names
# (unlike linkage_spec(fleet = ), which accepts either):
# Rceattle::reweight_comps(ss_run, fleets = c(1, 2))


# Doing one pass by hand, for reference. Note the weight has to be set on the
# parameter the fit reads, not only the fleet_control column: since 5.3.0 a
# model warm-started from `inits` takes its composition weights from there.
# Building afresh (inits = NULL, below) reads the column.
BS2017SS_weighted <- BS2017SS
BS2017SS_weighted$fleet_control$Comp_weights <-
  ss_run$data_list$fleet_control$Comp_weights_mcallister

ss_run_onepass <- Rceattle::fit_mod(data_list = BS2017SS_weighted,
                            inits = NULL, # Build afresh so the column is read
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            fit_control = fit_control(
                              phase = TRUE,
                              verbose = 1))
