# Rceattle configuration switch maps

Internal lookup tables that translate the human-readable switch strings
a user supplies (e.g. `"Logistic"`, `"NPFMC"`) into the integer codes
the TMB template consumes (and back again, via `revert_switches()`).
They are the single source of truth for the set of valid switch values;
`validate_switches()` errors on anything not listed here. See the
lifecycle note at the top of `R/0-switches.R` for how the maps are
applied.

## Details

- `sel_map`:

  Selectivity form (`fleet_control$Selectivity`).

- `tv_sel_map`:

  Time-varying selectivity structure (`fleet_control$Time_varying_sel`).

- `q_map`:

  Catchability form (`fleet_control$Catchability`).

- `tv_q_map`:

  Time-varying catchability structure (`fleet_control$Time_varying_q`).

- `comp_loglike_map`:

  Composition likelihood (`fleet_control$Comp_loglike` and
  `CAAL_loglike`).

- `fleet_map`:

  Fleet type (`fleet_control$Fleet_type`).

- `initMode_map`:

  Initial age-structure mode (`data_list$initMode`; see
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)).

- `suitMode_map`:

  Predator-prey suitability mode, per predator (`data_list$suitMode`).

- `hcr_map`:

  Harvest control rule (see
  [`build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)).
