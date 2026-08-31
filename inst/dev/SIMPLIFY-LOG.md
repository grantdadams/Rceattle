# Simplification log

Opportunities to clarify or simplify a workflow, an API or a document, found
while working in the code. **Add here; do not act on an entry unasked.** Comment
and roxygen tightening in a file a change already edits goes into that change;
anything that alters an API, a switch, a default or a workflow waits for Grant
to pick it from this list. One line each: where, what, the simpler form, and
whether behaviour changes.

| Where | What | Simpler form | Behaviour change |
|---|---|---|---|
| `fit_mod()` / `run_config()` | `random_sel` and `random_q` are global flags whose only job is a per-fleet decision `fleet_control` already determines (`.rce_np_unintegrable_fleets()`) | derive per fleet; deprecate the flags | yes, needs a deprecation path |
| `save_config()` | writes no `fleet_control`, so the YAML is not the whole model | serialise it | additive |
| `Catchability` / `Time_varying_q` | `Time_varying_q` is a mode, or comma-separated `env_data` column indices under `Catchability = "Environmental"` (`R/3-build_map.R`, `# FIXME: use formula`) | the q linkage, one grammar | yes, needs a shim |
| `fit_mod(initMode =)` | overwrites `data_list$initMode` unconditionally, so a switch stored on the data is never read (`BS2017MS$initMode = 1` never reached the golden fit) | read the data's value when the argument is missing, as `model_config` does | yes; moves `ms` golden |
| `estDynamics` | under `msmMode = 0` codes 1 and 2 estimate nothing, so two codes describe one model | one code in single-species mode, or estimate the scalar there | yes |
| `.map_switch()` vs `.conv()` | two switch normalisers; PR 1 gave `.map_switch()` the same factor and numeric-string handling | one should call the other | internal |
| `switch_check()` `.bad_ed` | the `estDynamics` range check repeats `.map_switch()`'s error, because `.map_switch()` passes numeric codes through for every switch (`msmMode`, `suitMode`, `Catchability`, the sd switches) | validate numerics in `.map_switch()` once each map's legacy numeric codes are known | yes, may reject old numeric codes |
| `build_map_fixed_natage()` | three blocks fix everything a fixed-numbers species would estimate; the linkage rows are refused in `fit_mod()` instead | one header comment saying it is the single place for "what a fixed species estimates" | internal |
| `.osa_bubble_plot()` | one function serves OSA and Pearson residuals, whose null distributions differ | a `type` argument carrying the outlier rule (PR 2) | internal |
| `linkage_spec(species =)` | a spec with no `species =` expands to one row per species, including species whose recruitment is input (`estDynamics > 0`), and is then refused; the estimated species have to be listed by hand | drop the fixed species from the expansion with a message | yes |
| `Sel_norm_scope`, `Sel_norm_bin` | inert unless the fleet is two-sex, normalized at a bin or the max, and not Hake/LogisticPM; three conditions the reader has to assemble | one sentence in the schema (PR 1) and a `data_check()` note (PR 1) | done in PR 1 |
| `R/0-osa_data.R:79-85` | five-line comment plus a FIXME for a one-line fact (`comp_offset` has three fill sites) | one line (PR 2) | doc |
| `CONTRIBUTOR-EXPERIENCE.md` item F | proposes work the code already does (`revert_switches()`) | delete (PR 6) | doc |
| `R/10-run_mse.R:900` | every deviation array but `log_M1_dev` is carried into the operating-model projection | carry the terminal year, as `index_q_dev` is | yes, `M1_re` operating models only |
| `materialize_linkage()` | the filter warnings (5.36.0) fire on every build, so a retrospective, jitter or MSE repeats them once per refit | hoist the check into the `build_*()` validators, or message once per fit | internal |
