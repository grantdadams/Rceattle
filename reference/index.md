# Package index

## Data input / output

Read, write, and manipulate Rceattle data lists.

- [`read_data()`](https://grantdadams.github.io/Rceattle/reference/read_data.md)
  : Read a CEATTLE excel data file
- [`write_data()`](https://grantdadams.github.io/Rceattle/reference/write_data.md)
  : Write data file
- [`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md)
  : Function to clean data prior to Rceattle runs
- [`combine_data()`](https://grantdadams.github.io/Rceattle/reference/combine_data.md)
  : Combine data sets. Will use the env_data data set from data_set1 and
  diet data will have to be updated.
- [`rearrange_dat()`](https://grantdadams.github.io/Rceattle/reference/rearrange_dat.md)
  : Rearrange
- [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
  : Function to check data for errors. Does not update the data set!

## Model setup

Build parameter lists, map objects, bounds, and biological inputs.

- [`build_params()`](https://grantdadams.github.io/Rceattle/reference/build_params.md)
  : Build parameter list from cpp file

- [`build_map()`](https://grantdadams.github.io/Rceattle/reference/build_map.md)
  : Main function to construct the TMB map argument for CEATTLE

- [`build_bounds()`](https://grantdadams.github.io/Rceattle/reference/build_bounds.md)
  : Build parameter bounds

- [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)
  : Define M1 specifications

- [`build_growth()`](https://grantdadams.github.io/Rceattle/reference/build_growth.md)
  : Specify the growth model for Rceattle

- [`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
  : Specify the stock-recruit relationship (SRR) for Rceattle

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
  : This function runs CEATTLE
- [`set_phases()`](https://grantdadams.github.io/Rceattle/reference/set_phases.md)
  : Function to set phasing order
- [`TMBphase()`](https://grantdadams.github.io/Rceattle/reference/TMBphase.md)
  : Run TMB using phases
- [`TMBAIC()`](https://grantdadams.github.io/Rceattle/reference/TMBAIC.md)
  : Calculate marginal AIC for a fitted model
- [`model_average()`](https://grantdadams.github.io/Rceattle/reference/model_average.md)
  : Model average of derived quantities
- [`rename_output()`](https://grantdadams.github.io/Rceattle/reference/rename_output.md)
  : Function to rename derived quantities from Rceattle

## Projections and MSE

Forward projections under alternative harvest strategies and closed-loop
management strategy evaluation.

- [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  : Run a management strategy evaluation
- [`load_mse()`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)
  : Function to load .RDs files from MSE runs
- [`check_mse()`](https://grantdadams.github.io/Rceattle/reference/check_mse.md)
  : Function to load .RDs files from MSE runs
- [`mse_summary()`](https://grantdadams.github.io/Rceattle/reference/mse_summary.md)
  : Management strategy evaluation performance metric summary
- [`remove_F()`](https://grantdadams.github.io/Rceattle/reference/remove_F.md)
  : Rerun with F = 0.
- [`sample_rec()`](https://grantdadams.github.io/Rceattle/reference/sample_rec.md)
  : Sample historical recruitment deviates and place in the projection

## Diagnostics and simulation

Retrospective analysis, jitter testing, and simulation.

- [`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md)
  : Retrospective peels
- [`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md)
  : Jitter analysis
- [`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md)
  : Simulate Rceattle data
- [`compare_sim()`](https://grantdadams.github.io/Rceattle/reference/compare_sim.md)
  : Evaluate simulation performance

## Data weighting

Composition data reweighting methods.

- [`calc_mcall_ianelli()`](https://grantdadams.github.io/Rceattle/reference/calc_mcall_ianelli.md)
  : Function to calculate McAllister-Ianelli weights
- [`calc_mcall_ianelli_diet()`](https://grantdadams.github.io/Rceattle/reference/calc_mcall_ianelli_diet.md)
  : Function to calculate McAllister-Ianelli weights for diet data
- [`check_composition_data()`](https://grantdadams.github.io/Rceattle/reference/check_composition_data.md)
  : Check and Clean Composition Data
- [`check_caal_data()`](https://grantdadams.github.io/Rceattle/reference/check_caal_data.md)
  : Check and Clean CAAL Data

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
- [`plot_ssb_depletion()`](https://grantdadams.github.io/Rceattle/reference/plot_ssb_depletion.md)
  : Plot SSB depletion (deprecated name)

## Plotting — fishing and natural mortality

Visualize fishing mortality, catch, and total mortality.

- [`plot_catch()`](https://grantdadams.github.io/Rceattle/reference/plot_catch.md)
  : Landings fits
- [`plot_f()`](https://grantdadams.github.io/Rceattle/reference/plot_f.md)
  : plot F
- [`plot_mortality()`](https://grantdadams.github.io/Rceattle/reference/plot_mortality.md)
  : Plot M1 + M2
- [`plot_m_at_age()`](https://grantdadams.github.io/Rceattle/reference/plot_m_at_age.md)
  : Plot natural mortality by age
- [`plot_m2_at_age_prop()`](https://grantdadams.github.io/Rceattle/reference/plot_m2_at_age_prop.md)
  : Plot predation mortality by age and predator

## Plotting — surveys and fit

Visualize survey indices, fits, and residuals.

- [`plot_index()`](https://grantdadams.github.io/Rceattle/reference/plot_index.md)
  : CPUE fits
- [`plot_logindex()`](https://grantdadams.github.io/Rceattle/reference/plot_logindex.md)
  : log(CPUE) fits
- [`plot_indexresidual()`](https://grantdadams.github.io/Rceattle/reference/plot_indexresidual.md)
  : CPUE residual
- [`plot_comp()`](https://grantdadams.github.io/Rceattle/reference/plot_comp.md)
  : Plot time series of comp data
- [`plot_data()`](https://grantdadams.github.io/Rceattle/reference/plot_data.md)
  : Timeline of data used in the model likelihoods
- [`plot_form()`](https://grantdadams.github.io/Rceattle/reference/plot_form.md)
  : Plot functional form

## Plotting — predation and diet

Visualize predation mortality, rations, and diet composition.

- [`plot_b_eaten()`](https://grantdadams.github.io/Rceattle/reference/plot_b_eaten.md)
  : Plot biomass eaten
- [`plot_b_eaten_prop()`](https://grantdadams.github.io/Rceattle/reference/plot_b_eaten_prop.md)
  : Plot biomass consumed of each prey species by predator
- [`plot_diet_comp()`](https://grantdadams.github.io/Rceattle/reference/plot_diet_comp.md)
  : Plot diet composition fits
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
- [`switch_check()`](https://grantdadams.github.io/Rceattle/reference/switch_check.md)
  : Function to check for missing switches for map and parameter
  functions
- [`convert_switches()`](https://grantdadams.github.io/Rceattle/reference/convert_switches.md)
  : Convert intuitive text strings to integer switches for TMB
- [`revert_switches()`](https://grantdadams.github.io/Rceattle/reference/revert_switches.md)
  : Convert integer switches to intuitive text strings. Maintains
  backwards compatability.
- [`get_growth_matrix_r()`](https://grantdadams.github.io/Rceattle/reference/get_growth_matrix_r.md)
  : Generate Length-at-Age Transition Matrix
- [`get_weight_at_age_r()`](https://grantdadams.github.io/Rceattle/reference/get_weight_at_age_r.md)
  : Calculate Predicted Weight-at-Age
- [`rich.colors.short()`](https://grantdadams.github.io/Rceattle/reference/rich.colors.short.md)
  : Make a vector of colors.
- [`t_col()`](https://grantdadams.github.io/Rceattle/reference/t_col.md)
  : \#https://www.dataanalytics.org.uk/make-transparent-colors-in-r/
- [`as.data.frame(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/as.data.frame.Rceattle.md)
  : Tidy long-format derived quantities from an Rceattle fit
- [`coef(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/coef.Rceattle.md)
  : Extract estimated parameters from an Rceattle fit
- [`logLik(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/logLik.Rceattle.md)
  : Log-likelihood of an Rceattle fit
- [`plot(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/plot.Rceattle.md)
  : Plot method for fitted Rceattle models
- [`residuals(`*`<Rceattle>`*`)`](https://grantdadams.github.io/Rceattle/reference/residuals.Rceattle.md)
  : Observed-vs-fitted residuals from an Rceattle fit
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
- [`validate_switches()`](https://grantdadams.github.io/Rceattle/reference/validate_switches.md)
  : Validates switches are correct
