# Making `NonParametric` + `Time_varying_sel = "IID"` integrable

Status: **proposed, not implemented.** 5.25.0 refuses `random_sel = TRUE` on this
combination. The inline markers in `src/TMB/ceattle.cpp` (`flt_sel_type == 2`
penalty block, search "PROPOSED") say where each piece goes.

## What is wrong now

Under IID the realized curve is `sel_coff + sel_coff_dev(yr)`, mean-centred.
Terms 1, 2 and 4 of the penalty block are charged on that realized curve, every
year, so all three are functions of the random effects. Two separate defects
follow, and 5.25.0's `NEWS` and `.rce_np_unintegrable_fleets()` only name the
first.

**(a) The kink.** Term 1 is `w * max(d, 0)^2`, whose second derivative is a step
at `d = 0`. The Laplace correction is `0.5 * log det H(u)`, so a step in `H`
steps the marginal objective. `Atka2022` halts at a maximum gradient of 6.8 and
reports a deviation sd 27% from what the same model reaches with
`Sel_curve_pen1 = 0` (maximum gradient 4e-4).

**(b) The unnormalized prior.** None of terms 1, 2 or 4 scales with
`sel_dev_sd`. The density actually integrated is

```
N(u; 0, sd) * exp(-w * P(u)) / Z(sd, w)
```

and `log Z(sd, w)` is nowhere in the objective. `Z` shrinks as `sd` grows, so
`sel_dev_sd` absorbs a term belonging to the normalizer. Terms 2 (`d^2`) and 4
(`avg_sel^2`) are already smooth: they contribute no kink and all bias. This is
the effect `vignette("model-options-and-functionality")` measures for term 4.

**Smoothing the hinge is not the fix.** Replacing `CppAD::abs(d)` with
`sqrt(d*d + eps*eps)` gives a C-infinity penalty and cures (a), but leaves (b)
untouched -- it buys a converged fit with a wrong variance, which is worse than
the current refusal because nothing announces it.

## The fix

Under IID the model is a base curve plus independent annual noise. The shape
prior belongs to the shape, and the shape is `sel_coff`, which -- unlike
`RandomWalk` -- `build_map()` leaves estimated. Charge terms 1, 2 and 4 once on
`sel_coff`; leave term 3 per year.

Then `H(u) = sd^-2 * I + d2(data)/du2`: smooth, no step, and the random-effect
prior is exactly `N(0, sd^2)` with its own normalizing constant. Both defects go,
and `Sel_curve_pen1 = 0` stops being a precondition. It also makes type 2
consistent with type 9, whose `avgsel` penalty already reads `sel_coff`.

## Three things to get right

1. **Gate on the mode, not on `random_sel`.** Under `RandomWalk` the base curve
   is mapped off, so there is nothing to hold a shape prior; leave that branch
   exactly as it is. `RandomWalk` stays non-integrable for its own reason -- the
   per-year renormalization leaves the level of each year's coefficients
   improper -- and no penalty change reaches that.
2. **Replicate the tail fill.** Term 1 runs to `nbins - 1`. On the realized
   curve, bins at or past `flt_n_sel_bins` are copies of the last estimated bin,
   so `d = 0` there and they contribute nothing. `sel_coff` holds 0 in those
   bins instead, so reading it unbounded would charge a penalty against a cliff
   the fitted selectivity does not have. Bound the loop at
   `flt_n_sel_bins(flt) - 1`, or fill forward before differencing.
3. **It moves the penalized fit too.** Same code path, so `random_sel = FALSE`
   changes with it, and the AMAK/Atka ADMB formulation charges the shape penalty
   on each annual curve. Needs `/golden-check` plus a before/after on `Atka2022`.
   No bundled dataset combined type 2 with IID before 5.25.0, so golden should be
   untouched -- but confirm the `nyrs_tmp` change does not reach a `RandomWalk`
   fleet such as `GOAcod`.

## The alternative, and why not

Keeping the penalty on the realized curve means putting `log Z(sd, w)` in the
objective. For the curvature term alone that is closed form -- `N(0, sd^2 I)`
times `exp(-w u'Du)` is Gaussian, so `log Z = -0.5 * log det(I + 2 w sd^2 D)` --
but term 1's hinge has no closed-form normalizer, so it would need the smoothing
of (a) *and* a numerical constant. Not worth it for a mode whose whole premise is
that the base curve carries the shape.
