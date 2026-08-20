# tools/

Ad-hoc verification harnesses. These are not part of `test_check()`: they run
real fits, take minutes to hours, and answer questions a unit test cannot. Run
them by hand (`Rscript tools/verify/<script>.R`) after a change to the
simulation, the likelihood, or the MSE loop.

Each script prints a pass/fail digest and leaves its numbers on stdout.

## verify/

**Equivalence — must stay bit-identical**

| Script | Proves |
|---|---|
| `verify-refit-like.R` | The `.refit_like()` refit path is unchanged across all nine entry points. Golden covers none of them. |
| `verify-mse-hindcast-invariant.R` | An MSE does not perturb the hindcast, independent of the refit refactor. |
| `verify-mse-om-horizon.R` | Refitting the operating model on a shortened projection horizon gives the same answer, and how much time it saves. |
| `verify-mse-repro.R` | Full-MSE reproduction digest, single- and multispecies. |

**Simulation — is the draw the distribution the likelihood assumes?**

| Script | Proves |
|---|---|
| `verify-sim-centering.R` | Simulated observations are centred where the likelihood expects. An sd check alone passes a mis-centred draw. |
| `verify-sim-index-families.R` | Each survey is drawn from its own `Index_distribution`, not from lognormal for everything. |
| `verify-sim-diet.R` | Stomach contents are drawn from the distribution the diet likelihood assumes. |
| `verify-sim-linkage-re.R` | Recruitment, M and growth written as a random linkage are drawn from the covariance they were given. |
| `verify-sim-process-error.R` | The process draws, and the two `M1_re` defects fixed alongside them. |

**Recovery — does a refit get back what was simulated?**

| Script | Proves |
|---|---|
| `verify-sim-recovery.R` | Which observation type costs a `self_test()` its recovery, starting from deterministic recovery. |
| `verify-sim-recovery-process.R` | Simulate a process, refit, and see what comes back. |
| `verify-sim-recovery-M.R` | Whether the data actually inform M. Reproduces the documented limit: deviations recover, their SD does not. |
| `verify-sim-local-optimum.R` | When recovery is poor, whether that basin is the MLE for the simulated data or a local optimum. |
