# N5 — Reformulate `Sel_curve_pen1/2/3` as purpose-named SD columns

## Why

`Sel_curve_pen1/2/3` are the worst-named public surface in Rceattle: three
columns whose meaning **changes with the selectivity form**, whose units are
opaque penalty weights, and one of which encodes a *direction in its sign*. The
CIE review and Grant both flagged them. Grant's ask: "reformulate them in terms
of standard deviations if we rename them (all the penalties are technically
squared/normal, right?)." — yes, they are (proof below), so a penalty weight
maps exactly to an SD.

## What the three columns actually mean (source: `ceattle_v01_11.cpp`)

Every penalty term is `weight · (·)²`, i.e. a Gaussian SSQ. For a weight `w`,
`w·x² = x²/(2σ²)` ⟹ **`σ = 1/√(2w)`**, **`w = 1/(2σ²)`**.

**Non-parametric forms (Selectivity 2 = NonParametric, 9 = NonParametricPM):**
- **col 0** = shape/smoothness penalty weight ([cpp:2949-2957](../src/TMB/ceattle_v01_11.cpp#L2949)).
  A *separate* flag `flt_sel_shape_mode` (0 = directional, 1 = smooth/two-sided)
  already selects the mode; within directional mode the **sign of col 0** picks
  the direction (`≥0` penalize decreasing, `<0` penalize increasing), and `|col0|`
  is the weight. Also used at [cpp:3076](../src/TMB/ceattle_v01_11.cpp#L3076) (type 9).
- **col 1** = curvature (2nd-difference) penalty weight ([cpp:2967](../src/TMB/ceattle_v01_11.cpp#L2967)).
- **col 2** = dev-magnitude penalty weight on the raw coefficient increments
  ([cpp:2984](../src/TMB/ceattle_v01_11.cpp#L2984), [cpp:3080](../src/TMB/ceattle_v01_11.cpp#L3080)).
- (The RW step penalty at [cpp:2974](../src/TMB/ceattle_v01_11.cpp#L2974) already
  uses `sel_dev_sd`, an SD — untouched here.)

**AR1 forms (Selectivity 6 = 2DAR1, 7 = 3DAR1):** the same columns hold **AR1
correlations**, not penalties: col 0 = `rho_age`, col 1 = `rho_year`
([cpp:3111-3112](../src/TMB/ceattle_v01_11.cpp#L3111)), col 2 = third rho / passed
natural-scale to `construct_Q` ([cpp:3143](../src/TMB/ceattle_v01_11.cpp#L3143)).
These are **not** SD-reformulatable.

## Design — R-side only, C++ untouched

**Keep the C++ contract fixed:** the template keeps reading the `sel_curve_pen`
`[n_flt × 3]` matrix exactly as today. All change is in R: `build_params()`
constructs `sel_curve_pen` from new, purpose-named `fleet_control` columns,
dispatching on the selectivity form. This guarantees the fit is **bit-identical**
to a model whose old columns held the equivalent weights, and confines risk to a
single, testable R construction.

**New `fleet_control` columns** (all optional; a fleet uses only those its form
needs):

| New column | Replaces | Meaning | C++ target |
|---|---|---|---|
| `Sel_shape_sd` | `Sel_curve_pen1` (mag) | SD of the shape/smoothness penalty | col0 = `±1/(2·sd²)` |
| `Sel_shape_dir` | `Sel_curve_pen1` (sign) | `"Decreasing"`/`"Increasing"` (directional mode only; ignored in smooth mode) | sign of col0 |
| `Sel_curvature_sd` | `Sel_curve_pen2` | SD of the curvature penalty | col1 = `1/(2·sd²)` |
| `Sel_devmag_sd` | `Sel_curve_pen3` | SD of the dev-magnitude penalty | col2 = `1/(2·sd²)` |
| `Sel_ar1_rho_age` | `Sel_curve_pen1` (AR1 forms) | AR1 age correlation | col0 (rho) |
| `Sel_ar1_rho_year` | `Sel_curve_pen2` (AR1 forms) | AR1 year correlation | col1 (rho) |
| `Sel_ar1_rho_3` | `Sel_curve_pen3` (AR1 forms) | third rho (3DAR1) | col2 (rho) |

**Conversion (`build_params`)**, per fleet, dispatching on `Selectivity`:
- Non-parametric (2, 9): `col0 = dir_sign · 1/(2·Sel_shape_sd²)` where
  `dir_sign = +1` for `"Decreasing"` / `-1` for `"Increasing"` (and irrelevant in
  smooth mode); `col1 = 1/(2·Sel_curvature_sd²)`; `col2 = 1/(2·Sel_devmag_sd²)`.
  A `Sel_*_sd` of `Inf`/`NA` ⇒ weight `0` (no penalty); guard `sd > 0`.
- AR1 (6, 7): `col0/1/2 = Sel_ar1_rho_age/year/3` (pass-through, no transform).
- Other forms: no `sel_curve_pen` (unchanged).

**Back-compat (the alias, released package):** if a fleet supplies the legacy
`Sel_curve_pen1/2/3`, use them **directly** as today (weights / rhos), emit one
deprecation message, and ignore the new columns for that fleet. If it supplies the
new columns, convert. This mirrors the `_sd` / `Q_init` alias contract and keeps
every existing model (and the 27 R / 52 test / 57 sibling references) working. The
bundled `.rda` (do any carry non-trivial `Sel_curve_pen`? — verify) are handled by
the pass-through, and can optionally be regenerated to the new columns in a
follow-up once equivalence is proven.

## Sub-commits (each compiles where relevant, `/golden-check`, suite green)

- **P1 — construction refactor, no interface change.** Move the `sel_curve_pen`
  matrix construction in [2-build_params.R:348](../R/2-build_params.R#L348) into a
  small, documented, form-dispatching helper that *today* still reads
  `Sel_curve_pen1/2/3`. Pure refactor. *Gate:* golden + GOA2018SS bit-identical.
- **P2 — new columns + conversion + alias.** Add the seven columns, the
  SD→weight / rho pass-through conversion, the legacy pass-through + deprecation.
  Defaults so an all-legacy or all-new model both work. *Gate (the real bar):*
  a fixture with active non-parametric penalties fit with the **new SD columns**
  is **bit-identical** to the same fit with the equivalent **legacy weights**
  (`w = 1/(2σ²)`), across a shape (both directions) + curvature + dev-magnitude
  case; an AR1-sel fixture with `Sel_ar1_rho_*` matches the legacy `Sel_curve_pen`
  rho encoding; golden unaffected.
- **P3 — docs + tests + NEWS.** `R/data.R` `@format` (replace the
  `Sel_curve_pen1/2/3` items), the selectivity vignette, committed conversion +
  alias tests, `NEWS.md` deprecation bullet. Optionally regenerate any bundled
  `.rda` that carry these columns.

## Hard gate

The **GOA2018SS** (and/or a constructed non-parametric + AR1 fixture) must fit
**bit-identically** through both the legacy and the reformulated columns — the
`w ↔ σ` conversion is where a subtle error would silently reweight composition
data and move the quota. Adversarial review of the conversion + the form dispatch
+ the direction-sign handling after P2, per the federal-rigor standard.

## Findings from digging into the code (scope-affecting)

- **`Sel_curve_pen1/2/3` is overloaded FOUR ways, not two.** By selectivity form:
  NonParametric (2/9) = penalty *weights* (SD-reformulatable); **LogisticPM (11)** =
  random-walk weights on slope/inflection/age-1 ([data_check.R:485-497](../R/1-data_check.R#L485));
  **2DAR1/3DAR1 (6/7)** = *logit-scale* AR1 correlations ([data_check.R:498-502](../R/1-data_check.R#L498));
  and a separate `Sel_shape_mode` column ("Directional"/"Smooth") already exists
  ([0-switches.R:386](../R/0-switches.R#L386)) with the *sign* of `Sel_curve_pen1`
  choosing decreasing-vs-increasing within directional mode. So an SD reformulation
  is clean ONLY for the NonParametric penalty weights; LogisticPM and AR1 uses must
  stay on the legacy columns (or get their own separate fields).
- **The σ↔w conversion is NOT floating-point exact.** `w = 1/(2σ²)`, `σ = 1/√(2w)`:
  `w=12.5→σ=0.2→12.5` and `w=200→0.05→200` fail `all.equal(tol=0)` (round to ~1e-15
  off). So migrating a model from weights to SD columns is equivalent only to
  ~machine precision — a strict bit-identity gate cannot use the *new-column* path.
  **Mitigation:** leave the legacy `Sel_curve_pen` path untouched (golden = exactly
  bit-identical), add SD columns as an ALTERNATIVE input converted in `switch_check`,
  and test the new-column path against the legacy with tolerance ~1e-10.

**Recommended scope:** NonParametric penalties only (`Sel_shape_sd` + `Sel_shape_dir`,
`Sel_curvature_sd`, `Sel_devmag_sd`), converted to `Sel_curve_pen1/2/3` in
`switch_check`, legacy untouched. LogisticPM / 2D-3D-AR1 keep the legacy columns.

## Honest scope note

This is the highest-effort, highest-risk item in the naming sweep: it changes a
*numeric input interface* (not just a name) for the minority of models using
non-parametric or AR1 selectivity, and the penalty↔SD conversion must be exact.
Because the C++ is untouched and the legacy columns pass through unchanged, the
blast radius is contained and every existing model keeps its exact numerics — but
P2's bit-identity gate is non-negotiable before it lands.
