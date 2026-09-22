# TODO: follow-ups found while preparing 5.34.0 (PR #144)

Status: **open except the `zero_N_pen` item, which shipped in 5.35.0.** Re-checked
2026-09-21 against 5.41.0. Pre-existing defects seen during the projection and MSE work.
None was fixed there because each changes an objective or a result on its own.

## PFMC's `Ftarget` is assigned after the hindcast years of the F loop

`ceattle.cpp` section 5.12 sets `Ftarget(sp) = Flimit(sp) + QnormHCR(sp)` inside the
fleet/sex/age/year loop, at the first projection year, so the hindcast years of `Ftarget_at_age`
for the first fishery, first sex and age `minage` carry `exp(log_Ftarget)`, the start value of a
parameter PFMC never estimates. `SPRtarget` reads a hindcast year, and the single-species
reference-point penalty moves `Flimit` with it: 0.4812, 0.5853 and 0.7379 from starts of 0, 3
and -999 on `make_test_data()` (`TRAPS.md`, "Silent-wrong-number traps"). Assign `Ftarget`
from `Flimit` before the loop, then re-run every single-species PFMC fit (the hake EM
`ss_run_DM_hcr_CSL` is one) and record the change. Results change.

## An estimated `log_Ftarget` is a single unbounded `nlminb` start

The single-species projection starts `log_Ftarget` at its `inits` value with no bounds
(`R/4-build_parameter_bounds.R`). 5.34.0 resets a no-fishing start (-999 or non-finite) to 0,
where the gradient is not degenerate, but a large start also sticks: NPFMC on `make_test_data()`
from `log_Ftarget = 3` (F = 20) stays at 20.08 with objective 83886.99 against 83877.42 from
0. A stored `CMSY` or `ConstantFSSB` fit with a large estimated `Ftarget` handed on as `inits`
is the realistic route. Bound the parameter, or start every estimated `log_Ftarget` at 0 as
the multispecies loop does.

## ~~`zero_N_pen` over-counts~~ -- FIXED in 5.35.0 (`c20f2ec7`)

`penalty` was reset per section rather than per cell, so `zero_N_pen(sp) += penalty` added a
running total that also carried earlier species'. Every floor site now declares its own
accumulator (`{ Type pen_N = 0; ... zero_N_pen(sp) += pen_N; }`, nine sites in `ceattle.cpp`)
and `check_convergence()` reports the row. **Do not re-do this**: the prescription here was to
reset the penalty, re-run golden and refit the live assessments, which is a day's work against
a fix that shipped.

## `run_mse(regenerate_past = TRUE)`'s average-F refit never runs

`R/10-run_mse.R:647` tests `em$data_list$HCR == 2`, but the rule reaches it under its name
(`"ConstantF"`) on the string path, so the branch is dead there. `.normalize_hcr()`'s own
comment says either form may be stored depending on the processing path, so it is dead on that
path rather than universally. `.normalize_hcr()` in `R/10-mse_summary.R` is the
existing way to compare either spelling. If the branch is revived, `Ftarget` needs a full
per-species vector: `avg_F$avg_F` covers only species with a fleet in `fleet_control`, and
`extend_length()` stops on any other length.

## `goa_ss` has a second local minimum 52.9 units higher

`TRAPS.md`, "Coverage gaps in the golden check", has the account. The golden recipe needs a
robustness fix so a one-ULP gradient change cannot move it: a warm start from the reference
`par`, or a second start keeping the lower negative log-likelihood. Until then a `goa_ss` delta
of 52.9 with the other three bit-identical is the known second minimum.
