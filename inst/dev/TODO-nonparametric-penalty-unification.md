# TODO: share one penalty implementation between NonParametric (2) and NonParametricPM (9)

Status: not started. Raised 2026-09-01 while adding `Sel_penalty_form` (5.25.0).

## The observation

`Sel_penalty_form = "AMAK"` makes type 2 score its penalties the way type 9 already
does. Type 9's random walk is already the bare sum of squares, and its level weight
already comes from `Sel_avgsel_pen`, so the switch adds to type 2 exactly what type 9
already had. That is an argument that the penalties, not the selectivity types, are
the wrong unit of organization.

What genuinely separates the two types is how the curve is BUILT, and that has to
stay a `Selectivity` type: type 2's `sel_coff` are log-selectivity directly, while
type 9's deviates are random-walk increments accumulated by cumulative sum
(`selectivity.hpp`, cases 2 and 9). Everything after "here is `log_non_par_sel`" is
penalty algebra that the two blocks compute independently, in
`ceattle.cpp` sections 1a (type 2) and the `flt_sel_type(flt) == 9` block below it.

## How they have drifted

| term | NonParametric (2) | NonParametricPM (9) |
| --- | --- | --- |
| decreasing / shape | whole range, decreasing only | `Sel_pen_first_bin..Sel_pen_last_bin`, directional or smooth |
| curvature | whole range, `JNLL_SEL_NONPARAM` | skips the base year when `Sel_start_year > styr`, `JNLL_SEL_DEV` |
| random walk | `dnorm`, or bare SSQ under `Sel_penalty_form = "AMAK"` | always bare SSQ |
| dev magnitude | absent | `Sel_curve_pen3 * sel_coff_dev^2` |
| level | `avg_sel(yr)^2` summed over years | `log_mean_exp(base coffs)^2`, once per sex |

Two of those are defects rather than choices:

* **The curvature term lands in a different `jnll_comp` row for the two types**, so the
  same penalty is reported under two different names. `.JNLL_ROW_AXIS` and the display
  names in `R/6-rename_output.R` cannot distinguish them.
* **The level penalty is computed on a different quantity** -- a per-year centering
  constant summed over years, against the base curve once.

## Measured against ADMB (2026-09-01)

Bridging `Rceattle-models/BSAI atka mackerel` against `Data/mod23/amak.tpl`, with the
penalty weights taken from `Data/mod23/input.log` and mapped onto the columns Rceattle
reads, the fishery selectivity likelihood comes to **98.5212** against ADMB's
**98.9469**, and the survey to **4.7809** against **4.4366**. Both fleets reconcile
exactly against a hand reconstruction of the package's own algebra, so the residual is
range, not arithmetic:

* the walk runs over `nbins - 1` increments (`ceattle.cpp`, the type-2 branch) where
  ADMB's `norm2` (`amak.tpl:2524`) runs over all `nbins`, double-counting the
  plus-group increment because `log_sel` is flat above `nselages`;
* the decreasing penalty covers every adjacent pair where ADMB starts at
  `seldecage = int(nages/2)` (`amak.tpl:499, 2527`) -- i.e. `Sel_pen_first_bin`. On the
  atka curve this one happens to cost nothing, because the decreasing violations all
  sit above that bin, but it is a real difference on any other curve;
* ADMB charges these only at `yrs_sel_ch` (selectivity change years) rather than every
  year. For this model `n_sel_ch_fsh = 45`, every year, so it does not bite here.

Type 9 already uses the whole-`nbins` walk bound (`ceattle.cpp`, the type-9 branch),
which is a second reason to share one implementation.

## The live trap

Type 2 silently ignores `Sel_pen_first_bin`, `Sel_pen_last_bin`, `Sel_shape_mode`,
`Sel_start_year` and `Sel_curve_pen3`; type 9 honours all five. `data_check()`
(`R/1-data_check.R`, the `Sel_pen_first_bin` / `Sel_pen_last_bin` block) VALIDATES the
first two and will reject an out-of-range value on a `NonParametric` fleet, after which
the model discards the setting. Someone giving a type-2 fleet a penalty range has no
way to learn it did nothing.

## The standard-deviation columns already exist -- finish the job

`Sel_shape_sd`, `Sel_curvature_sd` and `Sel_devmag_sd` (`R/0-column_schema.R`,
applied in `R/0-switches.R`) already express the three penalties as standard
deviations and set the corresponding `Sel_curve_pen1/2/3` to `1/(2*sd^2)`. So the
intuitive form is there; what remains is that it is an alternative rather than the
convention, and that the sd still cannot be estimated.

Worth recording, because it is what makes the ADMB bridge legible: **an AMAK
control file's raw inputs ARE these sds, one to one, with no arithmetic.**
`amak.tpl` reads a curvature input and stores `1/(2*input^2)` (line 948), and
squares the decreasing input into a variance (line 615) that it then divides
`0.5*d^2` by -- both exactly Rceattle's `1/(2*sd^2)`. Verified on
`Data/mod23/input.log`:

| ADMB input | Rceattle column | weight |
| --- | --- | --- |
| `curv_pen` 0.5 / 0.946 | `Sel_curvature_sd` | 2 / 0.5587117 |
| `seldec_pen` 0.16 / 106.09 | `Sel_shape_sd` | 19.53125 / 4.442435e-05 |

What is left:

* **`Sel_curve_pen1/2/3` remain the primary form** -- the schema, the C++
  (`sel_curve_pen`) and every bundled workbook carry weights, and the sd columns
  are converted away in `switch_check()`. A fitted object reports the weight, so
  a user who supplied an sd cannot read one back.
* **The sd is still not estimable**, which is the substantive half. Without a
  `log(sigma)` term a free sigma runs to infinity, so these can only ever be
  fixed constants. Adding the constant changes the objective by a constant and
  breaks every pinned reference number, so it belongs with the unification below
  rather than as a standalone change.
* **`Sel_devmag_sd` sets `Sel_curve_pen3`, which type 2 ignores** -- see the live
  trap above. Supplying it on a `NonParametric` fleet does nothing.

## Sequence

1. **Factor the shared terms into `selectivity.hpp`**, each block calling them with
   exactly the arguments it uses today. No numeric change by construction; the golden
   suite confirms it.
2. **Make type 2 read the five columns it ignores**, defaulting to today's behaviour
   (whole range, directional, `styr`, weight 0). Numeric only for a model that sets
   them, and no shipped model does -- setting them is currently a no-op. Until this
   lands, `data_check()` should warn that they are ignored on type 2.
3. **Unify the curvature row and the level quantity.** This one moves existing type-9
   numbers, so it needs the GOA pollock bridge re-run, not only `/golden-check`.

Steps 1 and 2 are safe. Step 3 needs a real refit to justify.

## Related

* `Sel_penalty_form` (5.25.0) is the switch that exposed this; its `"AMAK"` branch is
  the code steps 1-3 would absorb.
* `inst/dev/TRAPS.md` records that `TMB::compile()` does not track `.hpp` headers, which
  matters for step 1 -- moving code into `selectivity.hpp` will not rebuild on its own.
