# Next steps from the IPHC halibut trial

Status: **open.** Four items from an external IPHC user building a two-sex,
random-effects assessment entirely in R, 2026-09-04 to 09-10. Two have their
own design notes already; this file is the index and holds the two that do not.

The user's model: two sexes, logistic selectivity with RE deviations on the
survey, non-parametric on the fishery, and **sex-specific composition for only
part of the fishery time series** — which is why they want one offset parameter
rather than two free curves.

## Already settled — do not redo

- **Per-sex selectivity linkage fixed both sexes.** Fixed in 5.28.1
  (`564475b8`); `test-linkage-selectivity-per-sex.R` guards it.
- **`random_sel` confusion.** The user set `random_sel = TRUE` believing linkage
  REs needed it. They do not: `random_sel` gates `log_sel_slp_dev` /
  `sel_inf_dev` / `sel_coff_dev` (`fit_mod()`, `if (random_sel)`), while linkage
  REs enter through `beta_linkage_re`, gated independently on the map. Verified
  by fitting both on `GOAatf`: with `random_sel = FALSE` the survey's
  `rw(1 | Year)` still integrates (57 random effects, all `beta_linkage_re`).
- **Which forms can separate the sexes by level.** Documented in
  `vignette("model-options-and-functionality")` and
  `inst/dev/two-sex-selectivity-example.R` §8.

---

## 1. Sex-specific apex offset — the one the user needs

**See `inst/dev/TODO-sex-apical-offset.md` for the design.** This section records
only what the email thread added.

The request, verbatim: *"If there was an ability to specify which logistic (male
or female) was treated as an offset and then estimate that offset rather than set
the logistic to go to 1.0 that might solve all my problems."*

Three things follow that the design note should be read against:

- **It must work on a plain `Logistic`.** The request is not to change form. The
  double-logistic route works only because that form's peak height is a free
  function of its parameters; a logistic tends to 1 for either sex whatever the
  parameters, so no normalization setting reaches it.
- **A prior on `inf_desc` is not this.** Suggested as a stopgap during the
  thread, but it puts a prior on the male *descending inflection age* of a
  double logistic and lets the peak fall out of the shape indirectly. Measured on
  `GOAatf`: free per-sex parameters give a male:female max ratio of 1.9687, the
  tight `normal(8, 0.1)` prior 1.9294 — the prior barely moved the ratio, it
  constrained the shape. A loose `normal(8, 1)` reached 2.80 with **no
  sdreport**. It is a shape prior standing in for an apex offset and will not
  behave predictably on their data.
- **Partial sex-specific comps make the offset parameterization more valuable,
  not less.** One parameter informed by the sex-structured subset, with the whole
  series informing a shared shape, is the case Stock Synthesis's male-offset
  option exists for. Two free curves each see a fraction of the data — the
  user's stated problem with the two-sex non-parametric.

Identifiability is unchanged: the ratio is informed only by joint
(`comp_data$Sex = 3`) composition, so on their fishery it is informed only over
the years that have it. Say so when it ships.

---

## 2. Integrating non-parametric deviations

**See `inst/dev/TODO-nonparametric-iid-integrable.md` for the design.** One gap
found afterwards that the note does not cover:

**The guard is necessary but not sufficient.** `.rce_np_unintegrable_fleets()`
(`R/3-build_map.R`) tests `Sel_curve_pen1` only. Zeroing it removes the kink
and the fit converges — verified on `GOAatf`, `random_sel = TRUE` with both
shape penalties zeroed fits with 2,261 random effects across `beta_linkage_re`
and `sel_coff_dev`. But the avgsel term is hardcoded:

```cpp
jnll_comp(JNLL_SEL_NONPARAM, flt) += 2.0 * square(avg_sel(flt, sex, yr));  // ceattle.cpp
```

Weight 2.0, no column disables it, charged every year on the realized curve. It
is smooth, so it adds no kink and — in the design note's words — "all bias". The
same holds for `Sel_curve_pen2` if left non-zero. So a fit that passes the guard
still reports a deviation SD that maximizes a marginal missing its
`log Z(sd, w)` term, with nothing announcing it: the failure mode the note calls
"worse than the current refusal".

Either widen the guard to cover `Sel_curve_pen2` and the avgsel term, or fix the
underlying charge. Do not ship the `Sel_curve_pen1 = 0` workaround as advice
without saying the SD is not trustworthy.

Timing: the fix moves the **penalized** fit too, since it is the same code path,
so `random_sel = FALSE` models change with it. It reaches the GOA ATF assessment
this cycle — hold until after, then `/golden-check` plus a before/after on
`Atka2022`.

---

## 3. Too many switches, split across two places

In the user's words: *"it is challenging to navigate setting up switches in both
the fleetcontrol and model run sections… I am pretty much just doing all my work
in R rather than excel, so combining the fleet control and fit_mod sections would
really streamline things."*

Not a lift worth taking whole, but two pieces are cheap and both remove a switch
rather than adding one:

- **`random_sel` is derivable.** It is one global flag whose only job is deciding
  whether the `Time_varying_sel` deviations are integrated — and whether a fleet's
  deviations *can* be integrated is already computed per fleet from
  `fleet_control` alone (`.rce_np_unintegrable_fleets(fleet_control)`,
  `R/3-build_map.R`). Integrating every fleet that can be and penalizing the
  rest would delete the flag and the class of error the user hit. Behaviour
  change: needs a deprecation path, since `random_sel = FALSE` on an integrable
  fleet is a legitimate choice today.
- **Most of the merge already exists; the user has not been shown it.**
  `build_data()` takes `fleet_control` and `model_config` (linkages included) as
  named blocks, so one R object carries every switch, and `run_config()` /
  `save_config()` put the estimation controls in one YAML. Verified on 5.29.0,
  `GOAatf` with a survey `rw(1 | Year)` linkage: the pattern below integrates 57
  linkage REs, and the YAML round-trips the linkage exactly.

  ```r
  d   <- build_data(base = my_data, fleet_control = fc,
                    model_config = model_config(selFun = sel))
  cfg <- run_config(d, estimateMode = "Hindcast", random_sel = FALSE,
                    fit_control = fit_control(phase = TRUE))
  save_config(cfg, "halibut.yaml")
  fit <- fit_mod(d, config = cfg)
  ```

  What is left:
  - **The `config =` trap** (`TRAPS.md`). A config built from `model_config()`
    rather than from `d` silently replaces `d`'s linkages (`R/6-fit_mod.R:318`):
    57 REs become 0, no message. Warn when the two differ, or skip the overwrite
    when the config's `model_config` is the default.
  - **Fix the `random_sel` / `random_q` docs** (`R/0-save_config.R:305-306`).
    "Estimate time-varying selectivity as random effects" is what prompted the
    question whether `random_sel = FALSE` turns off the linkage REs. It does not:
    they are integrated unless the spec sets `integrate = FALSE`. `fit_mod()`'s
    `@param`s (`R/6-fit_mod.R:25-26`) name the `Time_varying_*_sd` columns but
    should say the same.
  - **`save_config()` writes no `fleet_control`,** so the YAML is not the whole
    model; the `fleet_control` still lives in the data object or workbook.
    Serializing it is the remaining step toward one file.
  - **`random_sel` cannot live on the data object** — `fit_mod()` overwrites
    `data_list$random_sel` from its argument (`R/6-fit_mod.R:373`). Moot if the
    first bullet derives it.

---

## 4. OSA outlier symbols scale with sample size

In the user's words: *"the number of extreme values depends on the total number
of residuals."* That is right, and it is a fixed threshold:

```r
osa$shape <- ifelse(abs(osa$residual) > 3, "outlier", "normal")   # R/7-plot_osa.R:354
```

Under the null, the expected count grows linearly with the panel's sample size:

| residuals in panel | expected \|resid\| > 3 |
|---|---|
| 100 | 0.3 |
| 500 | 1.3 |
| 2,000 | 5.4 |
| 5,000 | 13.5 |

So a length-composition panel always looks worse than an age panel on the same
model. Flag on the expected order statistic instead — a Bonferroni threshold
`qnorm(1 - 0.05 / (2n))` expects ~0.05 marked points per panel regardless of `n`
(3.48 at n = 100, 4.42 at n = 5,000). Cheap: one line plus the `@details` note at
`R/7-plot_osa.R:346`, which documents the fixed 3.

The user has raised this twice and calls it trivial; it is also the only one of
the four that is a half-hour job.

---

## Suggested order

1. **(4)** OSA threshold — smallest, and closes a point raised twice.
2. **(1)** apex offset — what the user needs, and the only one that unblocks
   their assessment.
3. **(2)** widen the guard now; the underlying fix waits on the ATF cycle.
4. **(3)** send the user the `build_data()` + `run_config(d)` pattern now, and
   fix the `config =` trap and the `random_sel` wording with it. The
   `random_sel` derivation waits on whether the deprecation path is acceptable.
