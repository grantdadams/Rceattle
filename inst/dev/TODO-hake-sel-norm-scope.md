# `Hake` selectivity ignores `Sel_norm_scope`

Status: **proposed, not implemented.** The inline marker in
`src/TMB/selectivity.hpp` (case 5, search "PROPOSED") says where the fix goes.

## What is wrong now

`Hake` (type 5) takes its normalization reference **inside** the sex loop
(`src/TMB/selectivity.hpp:490-505`), so the reference is always that sex's own
maximum. It is also excluded from the shared normalizer —
`normalize_and_project_selectivity()` gates on `sel_type != 5`
(`selectivity.hpp:73`) — so the `Sel_norm_scope` logic never runs for it.

The column is therefore **silently inert** on a Hake fleet. Measured on `GOAatf`
fleet 3 (`GOA_atf_fishery`, `nsex = 2`, `Time_varying_sel = "IID"`,
`Sel_norm_bin = 0`), the two settings are identical to every digit reported:

| `Sel_norm_scope` | female max | male max | ratio | jnll |
|---|---|---|---|---|
| `"AcrossSexes"` (default) | 1.0000 | 1.0000 | 1.0000 | 46.4710 |
| `"WithinSex"` | 1.0000 | 1.0000 | 1.0000 | 46.4710 |

Two consequences:

1. A user who sets `Sel_norm_scope` on a Hake fleet gets no error, no warning
   and no change — the one failure mode this package cannot afford.
2. Hake cannot carry a sex difference in selectivity **level**. Since F is
   shared across sexes (`F_flt_age = sel_at_age * exp(log_F)`,
   `src/TMB/ceattle.cpp:1270`), that means a Hake fleet cannot express
   sex-specific fishing mortality at all, whatever the data say.

A second, smaller gap sits in the same block: Hake honours a single named bin
(`sel_norm_bin1 >= 0 && sel_norm_bin2 < 0`) but falls through to the maximum
when a bin **range** is given, so `Sel_norm_bin_upper` is ignored too. The
shared normalizer averages over the range (`selectivity.hpp:84-88`).

## Does across- or within-sex even make sense here?

Across-sex, but pooling the reference is **not sufficient on its own** — and the
reason is specific to this form.

Hake's curve is a cumulative sum of `sel_coff` on the log scale, and the first
coefficient is mapped off at 0 for **both** sexes — `bins_on <-
(bin_first_selected + 1):N_sel_bins`, "first parameter is not-identifiable and
is not estimated" (`R/3-build_map.R:1100-1107`; `sel_coff` inits to 0,
`R/2-build_params.R:306`). So each sex's raw log-curve starts at exactly 0, and
the normalization decides **where the two sexes are welded together**:

| | reference subtracted | sexes forced equal at | free |
|---|---|---|---|
| `"WithinSex"` (today) | each sex's own max | the **peak** — both are 1 | the ratio at young ages |
| `"AcrossSexes"` | one pooled max `M` | the **first selected bin** — both `exp(-M)` | the ratio at the peak |

Neither frees the male:female level outright. Pooling only relocates the pin.

**Across-sex relocates it to the better place.** For a dimorphic stock the sexes
are similar in size where they first recruit to the gear and diverge as females
grow, so equal selectivity at the first selected bin is defensible; equal
selectivity at the peak asserts they are equally available exactly where they
differ most. Measured on `GOAatf` under today's within-sex rule, age-1
selectivity is 1.83e-4 for females and **2.25e-47** for males — with the peak
welded to 1, a near-step curve is the only way left for the model to say males
are less available overall. Suggestive rather than proof, but it is the shape
the constraint predicts.

**The complete fix frees one more parameter.** Under across-sex pooling, the
first coefficient stops being non-identifiable *for the second sex*: a common
shift of both sexes still cancels against `M`, but a shift of one sex relative
to the other does not. So estimate `sel_coff(flt, 2, bin_first_selected)` while
holding sex 1 at 0 as the reference. That single parameter **is** the
male:female level offset — the direct analogue of Stock Synthesis's
`Apical_Selex` (`SS_selex.tpl`, male parameter 5), and the same quantity
`inst/dev/` records as the deferred per-sex apical scalar. Keep it mapped off
under `"WithinSex"`, where it genuinely is not identified.

Identified only where the composition is joint (`comp_data$Sex = 3`); with
`Sex = 1`/`2` each sex's row is normalized independently and carries no
sex-ratio information, so the offset would rest on its prior.

## What it should do — narrowed

**The level question now belongs to `inst/dev/TODO-sex-apical-offset.md`.** A
per-sex apical offset gives Hake a sex level difference directly and
identifiably, applied after this form's own normalization. Freeing
`sel_coff`'s first bin for sex 2, sketched below, would be a *second* mechanism
for that same quantity — two grammars for one number, which is what we avoid.

So the remaining defect is the narrow one: **a column that silently does
nothing**. A user sets `Sel_norm_scope` on a Hake fleet and gets no error, no
warning and no change. Cheapest honest fix, and the one to do first:

- `data_check()` reports that `Sel_norm_scope` has no effect on a `Hake` fleet,
  the way it already reports other per-fleet settings a form ignores.
- Say the same in `R/0-column_schema.R`'s `Sel_norm_scope` description and in
  the vignette's normalization table.

Implementing the pooling itself is optional after that, and buys only
consistency of the column's meaning across forms — not a capability. The sketch
below is kept for whoever decides it is worth it.

## Fix sketch (only if the pooling itself is wanted)

Hake normalizes on the **log scale** (`exp(val - max_sel)`, `:503-505`) because
its pre-normalization curve is a cumulative sum of `sel_coff`. The shared
normalizer divides on the natural scale. `exp(a - m) == exp(a) / exp(m)`, so the
pooled reference is `max2()` over sexes of the per-sex **log** reference — not a
natural-scale divisor.

So: **do not simply drop `sel_type != 5` from the `:73` gate.** That would run
the natural-scale division on top of Hake's own log-scale normalization and
normalize twice.

Instead, move Hake's normalization out of the `switch` into a block after the
sex loop closes but still inside the year loop (`:596-597`), where both sexes'
un-normalized cumulative curves are already in `sel_at_age` / `sel_at_length`:

1. per sex, compute `max_sel(sex)` exactly as now (named bin, else max over bins);
2. if `sel_norm_scope(flt) == 1`, replace both with `max2(max_sel(0), max_sel(1))`;
3. then apply `exp(val - max_sel(sex))` per sex.

Do **not** also free `sel_coff`'s first bin for sex 2 to carry the level — that
is `TODO-sex-apical-offset.md`'s job, and doing both gives two parameters for
one quantity. Pooling here only relocates where the sexes are welded together;
it is not a substitute for the offset.

Fold rather than branch on the AD type, as the existing code does at `:501`.
While there, decide whether to honour `sel_norm_bin2` as a range mean — on the
log scale that is the mean of the logs, i.e. the geometric mean of the curve,
which is **not** what the shared normalizer computes; state whichever is chosen
in `R/0-column_schema.R`.

## Verification

- No bundled golden model uses Hake — `BS2017SS`/`BS2017MS` are `Selectivity`
  2/2/2/1 and `GOA2018SS` has no 5 — so `/golden-check` should be
  bit-identical. Confirm that rather than assume it.
- The fit **does** move for any two-sex Hake fleet left on the default
  `"AcrossSexes"`, which is a behaviour change: `NEWS.md` + `DESCRIPTION`
  `Version:` + the `Sel_norm_scope` row in
  `vignettes/model-options-and-functionality.Rmd`, which currently states the
  limitation as permanent ("**`Selectivity = "Hake"` always normalizes within
  sex**") and would become wrong.
- Add the two-scope pair above to `tests/testthat/` as a regression: today they
  are identical, and after the fix `"AcrossSexes"` must differ from
  `"WithinSex"` on `GOAatf` fleet 3.

See [[TRAPS.md]] for why `GOAatf` is the fixture for anything sex-structured.
