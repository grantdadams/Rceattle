# Maturity: what is open

Opened 2026-09-29, from the PR #178 review. Maturity gained its first length-based option in
5.46.0 (`L50_mat_len` / `slope_mat_len`), which raised the question in Open 1 and surfaced the
data anomaly in Open 2.

---

## Open 1 — should `maturity` carry a bin column and an age/length flag?

**Status: proposed, not implemented. Its own PR.** Raised by Grant: the composition sheet says
which axis its bins are on with `Age0_Length1`, so should `maturity` do the same instead of
being an age sheet with optional length scalars beside it?

The instinct is right — there is a real inconsistency — but the two are not interchangeable, and
the flag costs more here than it does for compositions.

### What the two sheets look like now

| | `comp_data` | `maturity` |
|---|---|---|
| shape | one row per observation | one row per species (wide) |
| bin columns | generic `Comp_1 .. Comp_n` | literally `Age1 .. Age30` |
| axis named by | `Age0_Length1` | nothing; the column names are the convention |
| read by | `comp_ctl` / `comp_obs` | `select(contains("Age"))`, `R/5-rearrange_data.R:864` |
| fitted? | yes, it has a likelihood | no, it is a fixed input |

So "age" is hard-coded in **two** places on the maturity side: the column names and the selector
that reads them. Adopting the flag means renaming the bin columns to something generic, which the
schema's `aliases` mechanism supports (rule 1, deprecate never delete), and fixing the selector.

### Why the flag is heavier than it looks

`Age0_Length1` on a composition row only picks which dimension of an existing prediction to
compare against. A length-based maturity row is not a reinterpretation:

- it has to be integrated over `P(L|age)` to give the maturity-at-age that SSB is built from,
  so it requires estimated growth (`growth_model > 0`) — the same prerequisite `data_check()`
  already enforces for `L50_mat_len` (`R/1-data_check.R:1924`);
- its bins land on the **population** length grid, not the data length bins. That is the
  data-ordinal-versus-population-ordinal problem 5.46.0 spent its review on, and 5.46.0's answer
  to that whole class was **refuse, don't translate** (`R/1-data_check.R:1973`). Designing a
  length-based maturity row means answering it for real, not refusing it.

### It is a superset, not a substitute

- A flagged sheet buys an **arbitrary empirical ogive at length**. That is a genuine new
  capability and nothing else provides it.
- `L50_mat_len` / `slope_mat_len` is **parametric**, and it is exactly SS3's maturity option 1,
  which is what the AI and GOA cod bridges have to reproduce. Binning a logistic throws away the
  two parameters.
- SS3 carries both for that reason. So this is an addition; the logistic shipped in 5.46.0 stays.

### Recommendation

Do it, in its own PR, and **reuse an existing spelling for the axis**. There are already three
ways to say "age or length" in this package:

- `Age0_Length1` — `comp_data`, `caal_data`
- `Selectivity_dimension` (`"Age"` / `"Length"`) — `fleet_control`
- maturity's implicit convention — an age sheet, plus length scalars on the control sheet

A fourth spelling would be a second grammar for something two existing columns already express.
Pick `Age0_Length1` (it is the one on a data sheet, which is what `maturity` is) and say so in the
schema `doc`.

Owed before it ships: the column-name migration through `aliases`, the `contains("Age")` selector,
a `write_data()`/`read_data()` round-trip test (a `data_list` element without that support
round-trips to nothing), and the population-grid bin decision above.

---

## Open 2 — `GOA2018SS` Cod maturity is 2.0 at ages 1-12, and nothing checks the range

**Status: DIAGNOSED 2026-10-12 — it is the first branch below (the data is wrong), not a
convention. Documented, not fixed; the remedy is Grant's call.** Full evidence, with the measured
two-sex doubling and the golden blindness, is in `TRAPS.md` under "Silent-wrong-number traps".
Pre-existing; noticed while answering Open 1, not introduced by 5.46.0.

Measured on the bundled dataset, 2026-09-29, **over every Age column**:

```
species 1 (Pollock)             min 0        max 0.9931
species 2 (Arrowtooth flounder) min 0        max 1
species 3 (Cod)                 min 0.9983   max 2       <- 2.0 for ages 1-12
```

Two corrections to those figures (2026-10-12). Read over the ages the model actually **uses**
(`nages = c(10, 21, 12)`), Cod is **min 2.0, max 2.0**: the 0.9983 is the Age13-20 padding, which
is beyond Cod's 12 ages and never read — and that padding is an ogive tail on the *undoubled*
scale, which is itself evidence the used cells were overwritten. And when reading the sheet,
drop the `Species` ID column first (`R/5-rearrange_data.R` does, via `select(contains("Age"))`);
including it makes species 2 and 3 report maxima of 2.0 and 3.0 in *every* dataset, which is a
trap of its own.

`maturity` is read straight into the template as a proportion: `mature_females(sp, age) =
maturity(sp, age) * sex_ratio(sp, age)` for a one-sex species (`ceattle.cpp:1015`), and Cod is
`nsex = 1` with `estDynamics = 0`, so the value feeds `spawn_output` and therefore SSB, SB0,
the depletions and the SPR reference points.

`data_check()` validates maturity **age coverage and gaps but not its range** — there is no
`[0, 1]` check anywhere (grepped `R/1-data_check.R`, 2026-09-29; the maturity rules there are
column count at `:451` and the NA/gap rules at `:479`).

Two possibilities were open. **The first is the answer** (2026-10-12):

- ~~the bundled GOA2018SS Cod ogive is wrong~~ — **this one.** 14 one-sex species rows across the
  11 bundled datasets carry maturity <= 1 and exactly 2 carry 2.0, both GOA Pacific cod. The
  control is decisive: BS2017SS/MS Cod is the same species, also `nsex = 1`, also
  `sex_ratio = 0.5`, and carries an ordinary ogive topping at 0.99728.
- ~~there is a convention for `maturity > 1`~~ — no. Greps over `R/`, `src/TMB/`, `vignettes/` and
  `inst/dev/` find no statement of one, and the per-dataset pattern a convention would produce is
  absent.

Three things the earlier version of this item did not have. **Cod's SSB is not in fact inflated**:
`ceattle.cpp:1025` multiplies by `sex_ratio` for a one-sex species, so `2.0 x 0.5 = 1.0`. The
exposure is a **two-sex** species, where that multiply does not apply and a 2.0 doubles spawning
output — measured exactly x2 on `GOA2018SS` Arrowtooth, objective unchanged. **Golden cannot see
it either way**, because `GOA2018SS` sets no `srr_fun`/`srr_pred_fun` and maturity never enters the
objective there (`dObj = 0`, against working positive controls), so correcting the data leaves both
GOA references bit-unchanged and no re-pin is owed. And the anomaly is **not confined to the stale
bundled copy**: 70 of 528 maturity rows across the ecosystem exceed 1, all exactly 2.0, all cod,
including the live GOA 2025 workbook.

**Still do not "fix" the data.** What Cod's maturity-at-age should be has to come from the
length-based ogive the SS3 model uses (`maturity_option: 1`), not from a guess, and the live
assessment workbooks would need the same correction. The range check in `data_check()` is still
owed — the upper half is already written and commented out at `R/1-data_check.R` (`# | mat_vals > 1`,
born commented in `ab9872bf`) — and enabling it as an **error** would refuse two bundled datasets
and 70 live rows, so it wants a warning first or a data fix alongside. Grant's call, 2026-10-12.
