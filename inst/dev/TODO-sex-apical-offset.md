# A per-sex apical selectivity offset

Status: **proposed, not implemented.**

## What is missing

Selectivity is the only route to a sex difference in fishing mortality:
`F_flt_age(flt, sex, age, yr) = sel_at_age(flt, sex, age, yr) * exp(log_F(flt, yr))`
(`src/TMB/ceattle.cpp:1270`) has one `log_F` per fleet-year, shared across sexes,
and `index_q` has no sex dimension at all.

But **no selectivity form has an apical height parameter**, so the level a sex
reaches is never something you can set or estimate directly:

| Form | Per-sex level |
|---|---|
| `Logistic`, `DescendingLogistic`, `LogisticPM` | none — the curve asymptotes to 1 for every sex |
| `NonParametric`, `NonParametricPM` | none — each sex is re-centred to mean 1, every year |
| `Hake` | none — normalizes within sex by its own max |
| `DoubleNormal` | tail plateau only (`right_floor` is sex-indexed); both sexes peak at exactly 1 |
| `DoubleLogistic`, `2DAR1`, `3DAR1` | emergent, not a parameter you can offset |

Per-sex offsets on the **shape** already work — `sex` is a valid `by` term, and
`linkage_spec(by = ~ fleet + sex)` reaches `slp_asc`, `slp_desc`, `inf_asc`,
`inf_desc`, `sigma_asc`, `sigma_desc`, `right_floor`. Only the height is absent.

**Naming trap.** Rceattle's `peak` is an *alias of `inf_asc`* — the ascending
inflection **age**, i.e. where the peak sits (`R/0-build_selectivity.R:31-33`).
Stock Synthesis's "peak selectivity" is how **high** it is. Do not describe the
existing `peak` linkage as the thing this TODO adds; a user wiring it up gets a
shifted curve, not a lowered one.

## Precedent

Stock Synthesis, `SS_selex.tpl`:

```c
Apical_Selex = sp(Maleselparm(f) + 4);            // 5th male offset parameter
asc = point1 + (Apical_Selex - point1) * (...);   // scales the male curve height
dsc = Apical_Selex + (point2 - Apical_Selex) * (...);
...
temp1 = max(max(sel_l(y,f,1)), max(sel_l(y,f,2)));
sel_l(y, f) /= temp1;                             // joint rescale over BOTH sexes
```

`seltype(f,3) = 3` makes males offsets from females, `4` the reverse — **the
user picks the reference sex**, which Rick's own description stresses: you need
to know which is the lesser-selected sex. The joint rescale on the last two
lines is already our `Sel_norm_scope = "AcrossSexes"`. We lack only the scalar.

## Design

- New `PARAMETER_ARRAY(log_sel_apical)`, `[n_flt, max_sex]`, init 0. Fully mapped
  off unless a linkage asks for it, so `exp(0) = 1` and a model without one is
  bit-identical.
- **Expose it through the linkage grammar**, as param `apical` with base array
  `log_sel_apical` — exactly parallel to `inf_asc` → `sel_inf`. No new
  `fleet_control` column: the grammar already carries per-sex strata, priors,
  bounds and phases, and gives "any parameter as an offset" for free.
- **Reference sex** is whichever sex the `sex =` filter does not name. That
  reproduces SS3's option 3 / 4 with no new switch. Note the existing usability
  trap: `sex =` is a **silent no-op unless `sex` appears in `by`**
  (`R/0-build_linkage.R:1066-1069`) — worth a warning while here.
- **Only one sex may carry it.** `log_F` is shared, so the common level is
  confounded with F and only the male:female *ratio* is identified. Estimating
  both sexes' apical would add a pure ridge; `data_check()` should refuse it.
- **Where it is applied:** after the form's `switch` in `selectivity.hpp`, before
  `normalize_and_project_selectivity()`. That placement is what makes it work for
  **every** form — it lands after `NonParametric`'s per-sex mean-centring
  (`:434-441`) and after `Hake`'s own internal normalization (`:490-505`), both of
  which would otherwise cancel it.
- **Normalization interaction.** The offset survives under `AcrossSexes`, and
  under `Sel_norm_bin = "Off"` where nothing is divided. Under `"WithinSex"` it is
  divided straight back out, exactly as the Hake column is today —
  `data_check()` must error on `apical` + `"WithinSex"` rather than fit a model
  whose parameter does nothing.
- **Mirrored fleets** index by fleet like every other selectivity block;
  `adjust_map_shared_params()` copies lead → mirror after the linkage masking, so
  a group shares one offset.

## What this subsumes

Applied after the form switch, `apical` gives **`Hake` a per-sex level without
touching its normalization**, so it takes over the level half of
`inst/dev/TODO-hake-sel-norm-scope.md`. That TODO is narrowed to what is left:
`Sel_norm_scope` silently does nothing on a Hake fleet, which needs a
`data_check()` warning whether or not this ships. Do not implement both the
apical offset and a freed first `sel_coff` bin — they are two parameters for one
quantity.

More generally, this is the answer for the whole class: one mechanism gives
every form a per-sex level, instead of fixing each form's normalization quirk
one at a time.

## Identifiability — say this in the docs, not just here

The ratio is informed **only** by joint composition (`comp_data$Sex = 3`), which
stacks both sexes into one vector summing to 1. With `Sex = 1`/`2` each sex is
normalized independently and carries no sex-ratio information, so the offset
would rest entirely on its prior. Stock Synthesis's manual says the same of its
own data types.

Measured on `GOAatf` (fleet 3, joint comps) with `DoubleLogistic` +
`AcrossSexes`, which is today's nearest approximation: female max 0.5080, male
1.0000, ratio 1.9687, PD Hessian, jnll 348.979. But the male descending limb sits
at a boundary — `log_sel_slp` SE 5391.6, `sel_inf` SE 8878.7 — and
`hessian_conditioning` puts **96% of the least-determined direction on `log_F`**.
That is the confounding this parameter has to live inside, so ship it with the
rule: quote the ratio from the sex whose parameters are determined, and treat a
standard error in the thousands as a boundary rather than an estimate.

## Verification

- **Default-off must be bit-identical.** `/golden-check`; no bundled dataset
  carries a linkage table, so the reference fits should not move at all.
- A new parameter block **requires** a `parameter_dictionary()` entry —
  `test-schema-parameter-index.R` asserts every block in `parList()` is in the
  dictionary, so that guard fires automatically if it is forgotten.
- Add the arms from `inst/dev/two-sex-selectivity-example.R` §8 with an
  `apical` offset on a `Logistic` fleet: today that form cannot produce a sex
  level difference at all, so a ratio ≠ 1 there is the proof the parameter works.
- `inits`: no saved `estimated_params` objects are carried across versions here,
  so the missing-block `stop()` at `R/6-fit_mod.R:612` needs no change. Revisit
  if that ever stops being true — adding a parameter block is otherwise a hard
  error for any warm start from an older fit.

Files: `src/TMB/ceattle.cpp`, `src/TMB/selectivity.hpp`, `R/2-build_params.R`,
`R/3-build_map.R`, `R/4-build_bounds.R`, `R/0-build_selectivity.R`,
`R/0-linkage_encode.R` + `src/TMB/linkage.hpp` (in lockstep),
`R/0-parameter_dictionary.R`, `R/1-data_check.R`.
