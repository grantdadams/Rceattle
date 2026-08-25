# TODO: get the MSE operating-model horizon saving back

Status: **open**, deliberately deferred. Opened 2026-08-24, alongside the change
that dropped the shortened horizon in 5.18.0.

## What was removed, and why

`run_mse()` used to refit the operating model only as far as the next assessment
year rather than over the whole projection. `projyr` sizes the AD tape at
roughly a fixed cost per model year, so this was worth **9% on a Bering Sea
multispecies MSE projecting to 2040, growing with the length of the projection**
(single-species models saw little change). It was introduced in 5.14.0 and
removed in 5.18.0.

It was removed because the shortened horizon set how far the random stream
advanced. `sim_mod()` draws once per observation row; `clean_data()` filters the
data frames to `styr:projyr`; so a shorter horizon carried fewer rows and
consumed fewer draws. The horizon was `assess_yrs[k + 1]` — **the next
assessment year** — so the number of draws taken at assessment `k` depended on
when the following assessment fell.

Two runs differing only in their assessment schedule therefore drew different
observation error from the first assessment onward, before the schedules had
diverged in any other way. Measured on BS2017SS (endyr 2017, projyr 2020,
`nsim = 1`, `seed = 666`), dropping the 2019 assessment moved 2019 catch by
**2.1%** — in a year whose advice is identical by construction, since both
schedules set it from the same 2018 assessment.

That is precisely the comparison an MSE of two assessment schedules exists to
make, so the saving lost to correctness. Per-year re-seeding (also added in
5.18.0) took the contamination to 0.53%; removing the horizon shortening took it
to 0.003%, which is numerical rather than stream divergence.

## The design that should get it back

The horizon serves two purposes at once, and they need different lengths:

| purpose | horizon needed | depends on the schedule? |
|---|---|---|
| simulate observations | `assess_yrs[k]` — later rows are discarded | **no** |
| `max_catch_hat` for the next iteration's cap | `assess_yrs[k + 1]` | yes |

Only the second needs the look-ahead, and `max_catch_hat` is the **only**
operating-model quantity the assessment loop reads (`R/10-run_mse.R`, the
exploitable-biomass limit in section 1). Everything else it consumes comes from
the estimation model. So:

```
per assessment k:
  refit the OM to assess_yrs[k]                   # SHORTER than 5.14.0-5.17.0
  simulate from it                                # horizon set by k alone
  project (estimateMode = 2) to assess_yrs[k + 1] # only to fill max_catch_hat
```

The simulation horizon then depends only on the current assessment year, which
two schedules share wherever they both assess — so the draws agree wherever the
two runs are genuinely in the same state, which is the most any single random
stream can give. And each object's quantities align with its own data frames, so
nothing has to trim a fit's per-observation quantities to match a trimmed data
frame. That was the hazard that made the obvious version of this dangerous:
`sim_mod()` matches quantities to data frames positionally.

## What is not known

Whether it is actually faster. The extension pass still builds an AD tape at
`assess_yrs[k + 1]`, and tape size is what the 9% was about — so per assessment
it trades a longer optimization for a shorter optimization plus one extra tape
build. That could net out faster, neutral or slower. **Measure before
committing to it.**

## Recovering the deleted code

The 5.14.0 truncate/restore helpers are the starting point for the refit-short
half, and the `estimateMode = 2` extension is new. They were deleted in the
5.18.0 commit on `mse-common-random-numbers`, along with
`tools/verify/verify-mse-om-horizon.R`, which was their equivalence and timing
harness:

```
git show <5.18.0-commit>^:R/10-run_mse.R
git show <5.18.0-commit>^:tools/verify/verify-mse-om-horizon.R
```

`.mse_trim_proj_params()` / `.mse_slice_year_dim()` are directly reusable;
`.mse_restore_om_horizon()` probably is not, since the new design never needs to
put a shortened fit back on the full projection.

## The other half: common random numbers are still not complete

Related, and worth doing in the same pass, because both are about which year a
draw belongs to.

Seeding is keyed to the **assessment** year, but one `sim_mod()` call draws
every year in the assessment interval. So a year sitting inside a longer
interval is drawn under that interval's seed, while the same year in a schedule
that assessed it directly is drawn under its own. On the verify fixture,
schedules `c(2018, 2019, 2020)` and `c(2018, 2020)` reach 2019 in the same state
— 2018's advice sets 2019 in both, and the catch agrees to 3e-5 — yet 2019's
simulated index is drawn under `seed_2019` in the first and `seed_2020` in the
second. That difference then enters the 2020 assessment as noise rather than as
the schedule.

Closing it needs the draw keyed to the year of the **observation**, not the year
of the assessment that collects it — per-observation-year streams inside the
`SIMULATE` blocks, or one `sim_mod()` call per projection year (much more
expensive). Neither is a small change, which is why 5.18.0 documents the gap in
`?run_mse` and `vignette("hcrs-and-mses")` rather than claiming more than it
delivers.

Measured effect of what 5.18.0 *did* fix, on the same fixture: 2.1% → 0.53%
(per-year seeding) → 0.003% (also refitting over the whole projection).

## The gate

`tools/verify/verify-mse-schedule.R` measures the contamination directly: the
`gap_crn_clean` check asserts that a year whose advice is identical under two
schedules comes out identical. Any horizon work has to keep that check passing,
and `verify-mse-repro.R` has to stay bit-identical for a single schedule.
