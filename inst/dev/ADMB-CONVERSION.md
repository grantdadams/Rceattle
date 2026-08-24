# Converting an ADMB model to Rceattle

## `dev_vector` sum-to-zero cannot be replicated in TMB

ADMB deviation vectors (`dev_vector`) carry a built-in *sum-to-zero* constraint. TMB has no
direct equivalent, so a straight port leaves the vector unidentified: the deviations and the
parameter they deviate from trade off freely, and the optimizer wanders along that ridge.

**Fix: turn off estimation of the first element** of the dev vector — map it to `NA` in
`R/3-build_map.R`. That pins the level and recovers identifiability.

**Exception:** if the dev vector already carries an additional likelihood penalty (a normal
penalty, or a random-effects density), the penalty pins it and every element can stay
estimated. Turning off the first element *as well* would then over-constrain it.

## Reference conversions in this repo

`src/TMB/selectivity.hpp` carries the AMAK conventions explicitly — e.g. AMAK evaluates the
logistic at `age_vector(j) = j + 0.5` (mid-age), which a port must reproduce rather than
"correct". `NonParametricPM` (`sel_map` 9) and `LogisticPM` (11) exist specifically to match
ADMB AMAK's "pm" parameterizations, penalties included.

The literature citations scattered through `src/TMB/` (AMAK, Ianelli, Punt, Holsman, Francis,
Methot, Kinzey & Punt) are the specification for those blocks, not historical notes. Keep them.
