# Package index

## Data input / output

Read, write, build, and manipulate Rceattle data lists.

- [`build_data()`](https://grantdadams.github.io/Rceattle/reference/build_data.md)
  : Build an Rceattle data list in R
- [`data_requirements()`](https://grantdadams.github.io/Rceattle/reference/data_requirements.md)
  : Report which data inputs a model configuration requires
- [`model_config()`](https://grantdadams.github.io/Rceattle/reference/model_config.md)
  : Build a model-configuration slot for a data list
- [`read_data()`](https://grantdadams.github.io/Rceattle/reference/read_data.md)
  : Read a CEATTLE excel data file
- [`write_data()`](https://grantdadams.github.io/Rceattle/reference/write_data.md)
  : Write a data list to a CEATTLE xlsx workbook
- [`write_template()`](https://grantdadams.github.io/Rceattle/reference/write_template.md)
  : Write a minimal starter CEATTLE data workbook
- [`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md)
  : Default-fill and tidy a data list before fitting
- [`combine_data()`](https://grantdadams.github.io/Rceattle/reference/combine_data.md)
  : Combine two Rceattle data lists
- [`rearrange_data()`](https://grantdadams.github.io/Rceattle/reference/rearrange_data.md)
  [`rearrange_dat()`](https://grantdadams.github.io/Rceattle/reference/rearrange_data.md)
  : Rearrange a data_list for TMB

## Model setup

Build parameter lists, map objects, bounds, and biological inputs.

- [`build_params()`](https://grantdadams.github.io/Rceattle/reference/build_params.md)
  : Build parameter list from cpp file

- [`build_map()`](https://grantdadams.github.io/Rceattle/reference/build_map.md)
  : Main function to construct the TMB map argument for CEATTLE

- [`build_bounds()`](https://grantdadams.github.io/Rceattle/reference/build_bounds.md)
  : Build parameter bounds

- [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)
  : Specify the residual natural mortality (M1) model for Rceattle

- [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
  : Specify the growth model for Rceattle

- [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
  : Specify the stock-recruit relationship (SRR) for Rceattle

- [`build_selectivity()`](https://grantdadams.github.io/Rceattle/reference/build_selectivity.md)
  : Selectivity specification

- [`build_catchability()`](https://grantdadams.github.io/Rceattle/reference/build_catchability.md)
  : Catchability specification

- [`build_composition()`](https://grantdadams.github.io/Rceattle/reference/build_composition.md)
  : Composition-weighting specification

- [`build_hcr()`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
  : Specify the harvest control rule (HCR) used for Rceattle

- [`build_hcr_map()`](https://grantdadams.github.io/Rceattle/reference/build_hcr_map.md)
  : Function to construct the TMB map argument for CEATTLE for
  projecting under alternative harvest control rules

- [`fit_control()`](https://grantdadams.github.io/Rceattle/reference/fit_control.md)
  :

  Bundle the optimizer / sdreport / phasing controls for
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)

- [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  : Capture a linkage specification

## Run configuration

Save, load, and inspect the full model + estimation setup as a YAML run
configuration.

- [`save_config()`](https://grantdadams.github.io/Rceattle/reference/save_config.md)
  : Save a model run configuration to a documented YAML file
- [`load_config()`](https://grantdadams.github.io/Rceattle/reference/load_config.md)
  : Load a model run configuration from a YAML file
- [`run_config()`](https://grantdadams.github.io/Rceattle/reference/run_config.md)
  : Extract the run configuration from a fit, data list, or config
  object

## Map helpers

Low-level helpers called by
[`build_map()`](https://grantdadams.github.io/Rceattle/reference/build_map.md)
to configure parameter estimation.

- [`adjust_map_shared_params()`](https://grantdadams.github.io/Rceattle/reference/adjust_map_shared_params.md)
  : Helper to adjust map for shared catchability/selectivity indices
- [`build_map_catchability()`](https://grantdadams.github.io/Rceattle/reference/build_map_catchability.md)
  : Helper to set map for Catchability parameters
- [`build_map_debug()`](https://grantdadams.github.io/Rceattle/reference/build_map_debug.md)
  : Helper to set map for debug mode
- [`build_map_f_and_data_weights()`](https://grantdadams.github.io/Rceattle/reference/build_map_f_and_data_weights.md)
  : Helper to set map for Fishing Mortality and Data Weights
- [`build_map_fixed_natage()`](https://grantdadams.github.io/Rceattle/reference/build_map_fixed_natage.md)
  : Helper to set map for Fixed N-at-Age models
- [`build_map_growth()`](https://grantdadams.github.io/Rceattle/reference/build_map_growth.md)
  : Helper to set map for growth parameters
- [`build_map_m1()`](https://grantdadams.github.io/Rceattle/reference/build_map_m1.md)
  : Helper to set map for Natural Mortality (M1) parameters
- [`build_map_predation()`](https://grantdadams.github.io/Rceattle/reference/build_map_predation.md)
  : Helper to set map for Predation Mortality (M2) parameters
- [`build_map_recruitment()`](https://grantdadams.github.io/Rceattle/reference/build_map_recruitment.md)
  : Helper to set map for Recruitment parameters
- [`build_map_selectivity()`](https://grantdadams.github.io/Rceattle/reference/build_map_selectivity.md)
  : Helper to set map for Selectivity parameters

## Model fitting

Fit single- or multi-species CEATTLE models and manage output.

- [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  : Fit the CEATTLE assessment model
- [`reweight_comps()`](https://grantdadams.github.io/Rceattle/reference/reweight_comps.md)
  : Iteratively reweight composition data (McAllister-Ianelli)
- [`set_phases()`](https://grantdadams.github.io/Rceattle/reference/set_phases.md)
  : Default phasing order for CEATTLE parameters
- [`TMBphase()`](https://grantdadams.github.io/Rceattle/reference/TMBphase.md)
  : Run TMB with ADMB-style phased estimation
- [`TMBAIC()`](https://grantdadams.github.io/Rceattle/reference/TMBAIC.md)
  : Calculate marginal AIC for a fitted model
- [`model_average()`](https://grantdadams.github.io/Rceattle/reference/model_average.md)
  : Model average of derived quantities
- [`rename_output()`](https://grantdadams.github.io/Rceattle/reference/rename_output.md)
  : Label the derived quantities reported by a CEATTLE fit

## Projections and MSE

Forward projections under alternative harvest strategies and closed-loop
management strategy evaluation.

- [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  : Run a management strategy evaluation
- [`load_mse()`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)
  : Load saved MSE simulation runs
- [`check_mse()`](https://grantdadams.github.io/Rceattle/reference/check_mse.md)
  : Check which saved MSE simulation files can be loaded
- [`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md)
  : Management strategy evaluation performance metric summary
- [`remove_F()`](https://grantdadams.github.io/Rceattle/reference/remove_F.md)
  : Rerun with F = 0.
- [`sample_rec()`](https://grantdadams.github.io/Rceattle/reference/sample_rec.md)
  : Sample historical recruitment deviations into the projection

## Diagnostics and simulation

Retrospective analysis, jitter testing, residual diagnostics, and
simulation.

- [`convergence_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/convergence_diagnostics.md)
  : Convergence diagnostics for a fitted Rceattle model
- [`osa_residuals()`](https://grantdadams.github.io/Rceattle/reference/osa_residuals.md)
  : One-step-ahead (OSA) residuals for an Rceattle model
- [`osa_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/osa_diagnostics.md)
  : Statistical diagnostics for OSA residuals
- [`process_residuals()`](https://grantdadams.github.io/Rceattle/reference/process_residuals.md)
  : Process residuals for an Rceattle model's random-effect processes
- [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md)
  : Retrospective peels
- [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md)
  : Jitter analysis
- [`profile(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/profile.Rceattle.md)
  : Likelihood profile across one or more parameter cells
- [`profile_components()`](https://grantdadams.github.io/Rceattle/reference/profile_components.md)
  : Likelihood components along a profile
- [`self_test()`](https://grantdadams.github.io/Rceattle/reference/self_test.md)
  : Self-test simulation analysis
- [`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md)
  : Simulate Rceattle data
- [`simulate(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/simulate.Rceattle.md)
  : Simulate data from a fitted Rceattle model
- [`compare_sim()`](https://grantdadams.github.io/Rceattle/reference/compare_sim.md)
  : Evaluate simulation performance
- [`print(`*`<Rceattle_retro>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_retro.md)
  : Print method for a retrospective analysis
- [`print(`*`<Rceattle_jitter>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_jitter.md)
  : Print method for a jitter analysis
- [`print(`*`<Rceattle_selftest>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_selftest.md)
  : Print method for a self-test
- [`print(`*`<Rceattle_profile>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_profile.md)
  : Print method for a likelihood profile
- [`print(`*`<Rceattle_mse_summary>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_mse_summary.md)
  : Print method for an MSE summary

## Diagnostics — shared arguments

The common vocabulary retrospective(), jitter() and self_test() share:
object, cores and fit_control, and what fit_control() reaches through a
refit.

- [`rceattle-refit-args`](https://grantdadams.github.io/Rceattle/reference/rceattle-refit-args.md)
  : Shared refit-diagnostic arguments

## Plotting — shared arguments

The common vocabulary the timeseries, predation and selectivity plotters
share: line_col, lwd, lty, alpha, species, spnames, minyr/maxyr,
incl_proj, incl_mean, add_ci and model_names.

- [`rceattle-plot-args`](https://grantdadams.github.io/Rceattle/reference/rceattle-plot-args.md)
  : Shared plotting arguments

## Plotting — population dynamics

Visualize biomass, recruitment, SSB, and depletion trajectories.

- [`plot_timeseries()`](https://grantdadams.github.io/Rceattle/reference/plot_timeseries.md)
  : Plot time-series
- [`plot_biomass()`](https://grantdadams.github.io/Rceattle/reference/plot_biomass.md)
  : Plot biomass
- [`plot_ssb()`](https://grantdadams.github.io/Rceattle/reference/plot_ssb.md)
  : Plot spawning stock biomass (SSB)
- [`plot_recruitment()`](https://grantdadams.github.io/Rceattle/reference/plot_recruitment.md)
  : Plot recruitment
- [`plot_depletion()`](https://grantdadams.github.io/Rceattle/reference/plot_depletion.md)
  : Plot biomass depletion
- [`plot_depletionSSB()`](https://grantdadams.github.io/Rceattle/reference/plot_depletionSSB.md)
  : Plot SSB depletion
- [`plot_exploitable_biomass()`](https://grantdadams.github.io/Rceattle/reference/plot_exploitable_biomass.md)
  : Plot exploitable biomass

## Plotting — fishing and natural mortality

Visualize fishing mortality, catch, and total mortality.

- [`plot_catch()`](https://grantdadams.github.io/Rceattle/reference/plot_catch.md)
  : Fishery catch fits
- [`plot_f()`](https://grantdadams.github.io/Rceattle/reference/plot_f.md)
  : plot F
- [`plot_mortality()`](https://grantdadams.github.io/Rceattle/reference/plot_mortality.md)
  : Plot M1 or M2 at age
- [`plot_m_at_age()`](https://grantdadams.github.io/Rceattle/reference/plot_m_at_age.md)
  : Plot natural mortality by age
- [`plot_m2_at_age_prop()`](https://grantdadams.github.io/Rceattle/reference/plot_m2_at_age_prop.md)
  : Plot predation mortality by age and predator

## Plotting — surveys and fit

Visualize survey indices, fits, and residuals.

- [`plot_index()`](https://grantdadams.github.io/Rceattle/reference/plot_index.md)
  : Survey index fits
- [`plot_indexresidual()`](https://grantdadams.github.io/Rceattle/reference/plot_indexresidual.md)
  : Survey index residuals
- [`plot_comp()`](https://grantdadams.github.io/Rceattle/reference/plot_comp.md)
  : Plot composition fits and residuals
- [`plot_profile()`](https://grantdadams.github.io/Rceattle/reference/plot_profile.md)
  : Plot the likelihood components across a profile
- [`plot(`*`<Rceattle_profile>`*`)`](https://grantdadams.github.io/Rceattle/reference/plot.Rceattle_profile.md)
  : Plot method for a likelihood profile
- [`plot(`*`<rceattle_osa>`*`)`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md)
  : Plot one-step-ahead (OSA) residual diagnostics
- [`plot_data()`](https://grantdadams.github.io/Rceattle/reference/plot_data.md)
  : Timeline of data used in the model likelihoods
- [`plot_form()`](https://grantdadams.github.io/Rceattle/reference/plot_form.md)
  : Plot functional form

## Plotting — predation and diet

Visualize predation mortality, rations, and diet composition.

- [`plot_b_eaten()`](https://grantdadams.github.io/Rceattle/reference/plot_b_eaten.md)
  : Plot biomass eaten by predation
- [`plot_b_eaten_prop()`](https://grantdadams.github.io/Rceattle/reference/plot_b_eaten_prop.md)
  : Plot biomass consumed of each prey species by predator
- [`plot_diet_comp()`](https://grantdadams.github.io/Rceattle/reference/plot_diet_comp.md)
  : Plot diet composition fits
- [`plot_diet_comp1()`](https://grantdadams.github.io/Rceattle/reference/plot_diet_comp1.md)
  : Plot diet composition fits (bubble/grid diagnostic)
- [`plot_diet_comp2()`](https://grantdadams.github.io/Rceattle/reference/plot_diet_comp2.md)
  : Plot diet composition fits (aggregation-aware)
- [`plot_ration()`](https://grantdadams.github.io/Rceattle/reference/plot_ration.md)
  : Plot ration

## Plotting — biology

Visualize selectivity, maturity, and stock-recruitment.

- [`plot_selectivity()`](https://grantdadams.github.io/Rceattle/reference/plot_selectivity.md)
  : Plot selectivity
- [`plot_selectivity_vs_maturity()`](https://grantdadams.github.io/Rceattle/reference/plot_selectivity_vs_maturity.md)
  : Plot fishery selectivity and maturity
- [`plot_maturity()`](https://grantdadams.github.io/Rceattle/reference/plot_maturity.md)
  : Plot maturity
- [`plot_stock_recruit()`](https://grantdadams.github.io/Rceattle/reference/plot_stock_recruit.md)
  : Plot stock recruit function

## Example datasets

Built-in data lists for the Bering Sea (BS), Gulf of Alaska (GOA),
Georges Bank, and other applications.

- [`BS2017SS`](https://grantdadams.github.io/Rceattle/reference/BS2017SS.md)
  : Data inputs for single species CEATTLE of the Bering Sea from 1979
  to 2017
- [`BS2017MS`](https://grantdadams.github.io/Rceattle/reference/BS2017MS.md)
  : Data inputs for multispecies CEATTLE of the Bering Sea from 1979 to
  2017
- [`EBS_ss_run`](https://grantdadams.github.io/Rceattle/reference/EBS_ss_run.md)
  : Fitted single-species CEATTLE model for the Eastern Bering Sea
- [`EBS_ss_M_run`](https://grantdadams.github.io/Rceattle/reference/EBS_ss_M_run.md)
  : Fitted single-species CEATTLE model with estimated M for the Eastern
  Bering Sea
- [`EBS_ms_run`](https://grantdadams.github.io/Rceattle/reference/EBS_ms_run.md)
  : Fitted multispecies CEATTLE model for the Eastern Bering Sea
- [`GOA2018SS`](https://grantdadams.github.io/Rceattle/reference/GOA2018SS.md)
  : Data inputs for a single-species Gulf of Alaska CEATTLE model (2018)
- [`GOApollock`](https://grantdadams.github.io/Rceattle/reference/GOApollock.md)
  : Data inputs for Gulf of Alaska walleye pollock CEATTLE model
- [`GOAsafe2018`](https://grantdadams.github.io/Rceattle/reference/GOAsafe2018.md)
  : Gulf of Alaska 2018 SAFE report reference values
- [`GOAatf`](https://grantdadams.github.io/Rceattle/reference/GOAatf.md)
  : Data inputs for Gulf of Alaska arrowtooth flounder CEATTLE model
- [`GOAatf2023`](https://grantdadams.github.io/Rceattle/reference/GOAatf2023.md)
  : Data inputs for Gulf of Alaska arrowtooth flounder CEATTLE model
  (2023)
- [`GOAcod`](https://grantdadams.github.io/Rceattle/reference/GOAcod.md)
  : Data inputs for Gulf of Alaska Pacific cod CEATTLE model
- [`GeorgesBank3spp`](https://grantdadams.github.io/Rceattle/reference/GeorgesBank3spp.md)
  : Data inputs for a three-species Georges Bank CEATTLE model
- [`NorthernRockfish2022`](https://grantdadams.github.io/Rceattle/reference/NorthernRockfish2022.md)
  : Data inputs for Northern Rockfish CEATTLE model (2022)
- [`Atka2022`](https://grantdadams.github.io/Rceattle/reference/Atka2022.md)
  : Data inputs for Atka mackerel CEATTLE model (2022)
- [`whamGrowthData`](https://grantdadams.github.io/Rceattle/reference/whamGrowthData.md)
  : Data inputs for CEATTLE model with WHAM-estimated growth

## Miscellaneous

S3 methods, utility functions, and internal helpers.

- [`print(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle.md)
  : Print method for fitted Rceattle models
- [`summary(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/summary.Rceattle.md)
  : Compact summary method for Rceattle fits
- [`print(`*`<summary.Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.summary.Rceattle.md)
  : Print method for an Rceattle model summary
- [`print(`*`<Rceattle_data>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_data.md)
  [`summary(`*`<Rceattle_data>`*`)`](https://grantdadams.github.io/Rceattle/reference/print.Rceattle_data.md)
  : Print method for an Rceattle data list
- [`switch_check()`](https://grantdadams.github.io/Rceattle/reference/switch_check.md)
  : Function to check for missing switches for map and parameter
  functions
- [`get_growth_matrix_r()`](https://grantdadams.github.io/Rceattle/reference/get_growth_matrix_r.md)
  : Generate Length-at-Age Transition Matrix
- [`get_weight_at_age_r()`](https://grantdadams.github.io/Rceattle/reference/get_weight_at_age_r.md)
  : Calculate Predicted Weight-at-Age
- [`rich.colors.short()`](https://grantdadams.github.io/Rceattle/reference/rich.colors.short.md)
  : Make a vector of colors.
- [`as.data.frame(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/as.data.frame.Rceattle.md)
  : Tidy long-format derived quantities from an Rceattle fit
- [`coef(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/coef.Rceattle.md)
  : Extract estimated parameters from an Rceattle fit
- [`logLik(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/logLik.Rceattle.md)
  : Log-likelihood of an Rceattle fit
- [`plot(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/plot.Rceattle.md)
  : Plot method for fitted Rceattle models
- [`residuals(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/residuals.Rceattle.md)
  : Residuals from an Rceattle fit
- [`vcov(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/vcov.Rceattle.md)
  : Variance-covariance matrix for an Rceattle fit
- [`prior_beta()`](https://grantdadams.github.io/Rceattle/reference/prior_beta.md)
  : Beta prior on a linkage coefficient
- [`prior_gamma()`](https://grantdadams.github.io/Rceattle/reference/prior_gamma.md)
  : Gamma prior on a linkage coefficient
- [`prior_lognormal()`](https://grantdadams.github.io/Rceattle/reference/prior_lognormal.md)
  : Lognormal prior on a linkage coefficient
- [`prior_normal()`](https://grantdadams.github.io/Rceattle/reference/prior_normal.md)
  : Normal prior on a linkage coefficient

## Deprecated

Superseded functions kept for backward compatibility. They still work
but emit a deprecation message; use the named replacement in new code.

- [`plot_logindex()`](https://grantdadams.github.io/Rceattle/reference/plot_logindex.md)
  : Survey index fits on the log scale (deprecated)
- [`plot_ssb_depletion()`](https://grantdadams.github.io/Rceattle/reference/plot_ssb_depletion.md)
  : Plot SSB depletion (deprecated name)
