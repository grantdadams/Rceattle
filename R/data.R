#' Data inputs for single species CEATTLE of the Bering Sea from 1979 to 2017
#'
#' A dataset containing the inputs used for CEATTLE
#'
#' @format
#' Control
#' \describe{
#' \item{nspp}{Number of species included in CEATTLE}
#' \item{styr}{Start year of the hindcast}
#' \item{endyr}{End year of the hindcast}
#' \item{projyr}{End year of the forecast}
#' \item{nsex}{Number of sexes to model in the population (1 = combined/1sex, 2 = models both female/male)}
#' \item{spawn_month}{Spawning month of the population to adjust the numbers spawning}
#' \item{R_sexr}{Percent of recruitment that is female (ignored if nsex = 1)}
#' \item{nages}{Number of ages of each species included in the hindcast}
#' \item{minage}{Minimum age for each population (i.e.does recruitment correspond to age 0, 1, 2?)}
#' \item{nlengths}{Number of lengths of each species included in the hindcast}
#' \item{pop_wt_index}{Weight-at-age (wt) index to use for calculation of each species population derived quantities (SSB, Consumption/Ration, Suitability, etc)}
#' \item{ssb_wt_index}{Weight-at-age (wt) index to use for calculation of each species spawning biomass}
#' \item{pop_age_transition_index}{Age transition matrix (e.g. growth trajectory) index to use for derived quantities of the population to convert age to length (also used in length-based predation estimation)}
#' \item{sigma_rec_prior}{Standard deviation to use for recruitment}
#' \item{other_food}{Other food in the ecosystem for each species}
#' \item{estDynamics}{Estimate or fix numbers-at-age: 0 = estimate dynamics, 1 = use input numbers-at-age in NbyageFixed, 2 = multiply input numbers-at-age (NbyageFixed) by a single scaling coefficient, 3 = multiply input numbers-at-age (NbyageFixed) by age specific scaling coefficient.}
#' \item{est_sex_ratio}{Is sex ration F/(M+F) to be included in the likelihood (assumed normal); 0 = no, 1 = use annual average across ages (uses 2nd age in propF data), 2 = age, and year specific (TBD)}
#' \item{sex_ratio_sigma}{Initial value or fixed value for sd of normal likelihood for sex ration. Not yet able to estimate.}
#' \item{est_M1}{Estimate residual (multi-species mode) or total natural mortality (single-species mode). 0 = use fixed natural mortality from M1_base, 1 = estimate sex- and age-invariant M1, 2 = sex-specific (two-sex model), age-invariant M1, 3 =   estimate sex- and age-specific M1.}
#' \item{fleet_control}{Survey and fishery data specifications}
#' \item{srv_biom}{Survey index in weight (kg) or numbers data}
#' \item{fsh_biom}{Total catch in weight (kg) or numbers data}
#' \item{comp_data}{Survey/fishery age or length composition data. Note if sex is 3, put female composition data then male composition data (similar to SS).}
#' \item{emp_sel}{Empirical/fixed selectivity for surveys and fisheries (leave empty if not used)}
#' \item{age_trans_matrix}{Age transition matrix (e.g. growth trajectory) used to convert age to length for length comp data. Can have multiple matrices for a species specified by Age_transition_index.}
#' \item{age_error}{Aging error matrices. Can have only one per species.}
#' \item{wt}{Weight-at-age (wt) to use for calculation of derived quantities (SSB, Consumption/Ration, Suitability, Total Catch, Survey Biomass, etc). Can have multiple weight-at-age data-sets for each species, but must be full for all years of the hindcast.}
#' \item{pmature}{Maturity-at-age for each species}
#' \item{propF}{Percent female at age for each species}
#' \item{M1_base}{Residual natural mortality for each species}
#' \item{Mn_LatAge}{Mean length-at-age for each species. Used when estimating time-invariant length-based gamma suitability (suitMode = 1) or time invariant length-based lognormal suitability (suitMode = 4)}
#' \item{aLW}{Parameters for weight-at-length power function for each species. . Used when estimating time-variant length-based gamma suitability (suitMode = 2) or time-variant length-based lognormal suitability (suitMode = 5)}
#' \item{Ceq}{Which bioenergetics equation to use for each species for ft to scale max consumtion: 1 = Exponential (Stewart et al 1983), 2 = Temperature-dependendence for warm-water species (Kitchell et al 1977; sensu Holsman et al 2015), 3 = temperature dependence for cool and cold-water species (Thornton and Lessem 1979)}
#' \item{Cindex}{Which environmental index in env_data to use to drive bioenergetics}
#' \item{Pvalue}{This scales the maximum consumption used for ration for each species; Pvalue is in Cmax*fT*Pvalue*Pyrs}
#' \item{fday}{Number of foraging days per year for each species}
#' \item{CA}{Intercept of allometric mass function for calculating maximum consumption: CA * Weight ^ CB}
#' \item{CB}{Slope of allometric mass function for calculating maximum consumption: CA * Weight ^ CB}
#' \item{Qc}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tco}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tcm}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tcl}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{CK1}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{CK4}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{env_data}{Environmental indices such as bottom temperature data to incorporate into ration equation specificed by Ceq and Cindex. Also used to drive catchability if Estimate_q = 5. Will use the mean for missing years. Temperature should be in celcius.}
#' \item{Pyrs}{Annual relative foraging rate by age. Multiplied by pvalue and fday to scale maximum consumption to the number of days in a year that foraging occurs.}
#' \item{UobsAge}{Stomach proportion by numbers for each predator, prey, predator age, prey age combination. Can also be year specific by including the column, "Year"}
#' \item{UobsWtAge}{Stomach proportion by weight for each predator, prey, predator age, prey age combination. Can also be year specific by including the column, "Year"}
#' }
#'
#' fleet_control: controls for survey and fisheries data
#'\describe{
#' \item{Fleet_name}{Name of survey or fishery}
#' \item{Fleet_code}{Index of survey/fishery ACROSS species}
#' \item{Fleet_type}{0 = Do not estimate; 1 = Fishery; 2 = Survey}
#' \item{Species}{Species number}
#' \item{Selectivity_index}{index to use if selectivities of different surveys are to be the same}
#' \item{Selectivity}{Selectivity to use for the species: 0 = empirical selectivity provided in srv_emp_sel; 1 = logistic selectivity; 2 = non-parametric selecitivty sensu Ianelli et al 2018; 3 = double logistic; 4 = descending logistic}
#' \item{Nselages}{Number of ages to estimate non-parametric selectivity for Selectivity = 2. Not used otherwise}
#' \item{Time_varying_sel}{Wether a time-varying selectivity should be estimated for logistic, double logistic selectivity, or descending logistic. 0 = no, 1 = penalized deviates given sel_sd_prior, 2 = random effect, 3 = time blocks with no penality, 4 = random walk following Dorn, 5 = random walk on ascending portion of double logistic only. If selectivity is set to type = 2 (non-parametric) this value will be the 1st penalty on selectivity.}
#' \item{Sel_sd_prior}{The sd to use for the random walk of time varying selectivity if set to 1. If selectivity is set to type = 2 (non-parametric) this value will be the 2nd penalty on selectivity.}
#' \item{Age_first_selected}{Age at which selectivity is non-zero}
#' \item{Acuumulation_age_lower}{Ages below this will be grouped to this age for composition data. For example, if set to 2, comp data for age 2 will include 1 and 2 year olds.}
#' \item{Acuumulation_age_upper}{Ages above this will be grouped to this age for composition data. For example, if set to 9 for a species with 10 ages, comp data for age 9 will include 9 and 10 year olds.}
#' \item{weight1_Numbers2}{Is the observation in weight (kg) set as 1, if the observation is in numbers caught, set as 2}
#' \item{Weight_index}{Weight-at-age (wt) index to use for calculation of derived quantities}
#' \item{Age_transition_index}{Age transition matrix (e.g. growth trajectory) index to use for derived quantities to convert age to length}
#' \item{Q_index}{index to use if catchability coefficients are to be set the same}
#' \item{Estimate_q}{Estimate catchability? (0 = fixed at prior; - 1 = Estimate single parameter; - 2 = Estimate single parameter with prior; - 3 = Estimate analytical q  from Ludwig and Walters 1994;  - 4 = Estimate power equation; - 5 - Linear equation log(q_y) = q_mu + beta * index_y)}
#' \item{Q_prior}{Starting value or fixed value for catchability}
#' \item{Q_sd_prior}{Variance of q prior: dnorm (log_q, log_q_prior, q_sd_prior)}
#' \item{Time_varying_q}{Wether a time-varying q should be estimated. 0 = no, 1 = penalized deviate, 2 = random effect, 3 = time blocks with no penalty; 4 = random walk from mean following Dorn 2018 (dnorm(q_y - q_y-1, 0, sigma). If Estimate_q = 5, this determines the environmental index to be used in the equation log(q_y) = q_mu + beta * index_y}
#' \item{Time_varying_q_sd_prior}{The sd to use for the random walk of time varying q if set to 1}
#' \item{Estimate_survey_sd}{Estimate survey variance (0 = use CV from srv_biom, 1 = yes, 2 = analytically estimate following (Ludwig and Walters 1994)}
#' \item{Survey_sd_prior}{Starting value to be used if Estimate_sigma_index = 1}
#' \item{Estimate_catch_sd}{Estimate fishery variance (0 = use CV from srv_biom, 1 = yes, 2 = analytically estimate following (Ludwig and Walters 1994)}
#' \item{Catch_sd_prior}{Starting value to be used if Estimate_sigma_catch = 1}
#' \item{Comp_weights}{Composition weights to be used for multinomial likelihood. These are multiplied. After running model, these will update to McAllister & Ianelli 1997 weights using the harmonic mean.}
#' \item{Catch_units}{Units used for survey: 1 = kg; 2 = numbers}
#' \item{proj_F_prop}{The proportion of future fishing mortality assigned to this fleet}
#' \item{Sex}{sex codes: 0=combined; 1=use female only; 2=use male only; 3 = joint female and male}
#'}
"BS2017SS"


#' Data inputs for multi-species CEATTLE of the Bering Sea from 1979 to 2017
#'
#' A dataset containing the inputs used for CEATTLE
#'
#' @format
#' Control
#' \describe{
#' \item{nspp}{Number of species included in CEATTLE}
#' \item{styr}{Start year of the hindcast}
#' \item{endyr}{End year of the hindcast}
#' \item{projyr}{End year of the forecast}
#' \item{nsex}{Number of sexes to model in the population (1 = combined/1sex, 2 = models both female/male)}
#' \item{spawn_month}{Spawning month of the population to adjust the numbers spawning}
#' \item{R_sexr}{Percent of recruitment that is female (ignored if nsex = 1)}
#' \item{nages}{Number of ages of each species included in the hindcast}
#' \item{minage}{Minimum age for each population (i.e.does recruitment correspond to age 0, 1, 2?)}
#' \item{nlengths}{Number of lengths of each species included in the hindcast}
#' \item{pop_wt_index}{Weight-at-age (wt) index to use for calculation of each species population derived quantities (SSB, Consumption/Ration, Suitability, etc)}
#' \item{ssb_wt_index}{Weight-at-age (wt) index to use for calculation of each species spawning biomass}
#' \item{pop_age_transition_index}{Age transition matrix (e.g. growth trajectory) index to use for derived quantities of the population to convert age to length (also used in length-based predation estimation)}
#' \item{sigma_rec_prior}{Standard deviation to use for recruitment}
#' \item{other_food}{Other food in the ecosystem for each species}
#' \item{estDynamics}{Estimate or fix numbers-at-age: 0 = estimate dynamics, 1 = use input numbers-at-age in NbyageFixed, 2 = multiply input numbers-at-age (NbyageFixed) by a single scaling coefficient, 3 = multiply input numbers-at-age (NbyageFixed) by age specific scaling coefficient.}
#' \item{proj_F}{Unused}
#' \item{est_sex_ratio}{Is sex ration F/(M+F) to be included in the likelihood (assumed normal); 0 = no, 1 = use annual average across ages (uses 2nd age in propF data), 2 = age, and year specific (TBD)}
#' \item{sex_ratio_sigma}{Initial value or fixed value for sd of normal likelihood for sex ration. Not yet able to estimate.}
#' \item{est_M1}{Estimate residual (multi-species mode) or total natural mortality (single-species mode). 0 = use fixed natural mortality from M1_base, 1 = estimate sex- and age-invariant M1, 2 = sex-specific (two-sex model), age-invariant M1, 3 =   estimate sex- and age-specific M1.}
#' \item{fleet_control}{Survey and fishery data specifications}
#' \item{srv_biom}{Survey index in weight (kg) or numbers data}
#' \item{fsh_biom}{Total catch in weight (kg) or numbers data}
#' \item{comp_data}{Survey/fishery age or length composition data. Note if sex is 3, put female composition data then male composition data (similar to SS).}
#' \item{emp_sel}{Empirical/fixed selectivity for surveys and fisheries (leave empty if not used)}
#' \item{age_trans_matrix}{Age transition matrix (e.g. growth trajectory) used to convert age to length for length comp data. Can have multiple matrices for a species specified by Age_transition_index.}
#' \item{age_error}{Aging error matrices. Can have only one per species.}
#' \item{wt}{Weight-at-age (wt) to use for calculation of derived quantities (SSB, Consumption/Ration, Suitability, Total Catch, Survey Biomass, etc). Can have multiple weight-at-age data-sets for each species, but must be full for all years of the hindcast.}
#' \item{pmature}{Maturity-at-age for each species}
#' \item{propF}{Percent female at age for each species}
#' \item{M1_base}{Residual natural mortality for each species}
#' \item{Mn_LatAge}{Mean length-at-age for each species. Used when estimating time-invariant length-based gamma suitability (suitMode = 1) or time invariant length-based lognormal suitability (suitMode = 4)}
#' \item{aLW}{Parameters for weight-at-length power function for each species. . Used when estimating time-variant length-based gamma suitability (suitMode = 2) or time-variant length-based lognormal suitability (suitMode = 5)}
#' \item{Ceq}{Which bioenergetics equation to use for each species for ft to scale max consumtion: 1 = Exponential (Stewart et al 1983), 2 = Temperature-dependendence for warm-water species (Kitchell et al 1977; sensu Holsman et al 2015), 3 = temperature dependence for cool and cold-water species (Thornton and Lessem 1979)}
#' \item{Cindex}{Which environmental index in env_data to use to drive bioenergetics}
#' \item{Pvalue}{This scales the maximum consumption used for ration for each species; Pvalue is in Cmax*fT*Pvalue*Pyrs}
#' \item{fday}{Number of foraging days per year for each species}
#' \item{CA}{Intercept of allometric mass function for calculating maximum consumption: CA * Weight ^ CB}
#' \item{CB}{Slope of allometric mass function for calculating maximum consumption: CA * Weight ^ CB}
#' \item{Qc}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tco}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tcm}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tcl}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{CK1}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{CK4}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{env_data}{Environmental indices such as bottom temperature data to incorporate into ration equation specificed by Ceq and Cindex. Also used to drive catchability if Estimate_q = 5. Will use the mean for missing years. Temperature should be in celcius.}
#' \item{Pyrs}{Annual relative foraging rate by age. Multiplied by pvalue and fday to scale maximum consumption to the number of days in a year that foraging occurs.}
#' \item{UobsAge}{Stomach proportion by numbers for each predator, prey, predator age, prey age combination. Can also be year specific by including the column, "Year"}
#' \item{UobsWtAge}{Stomach proportion by weight for each predator, prey, predator age, prey age combination. Can also be year specific by including the column, "Year"}
#' }
#'
#' fleet_control: controls for survey and fisheries data
#'\describe{
#' \item{Fleet_name}{Name of survey or fishery}
#' \item{Fleet_code}{Index of survey/fishery ACROSS species}
#' \item{Fleet_type}{0 = Do not estimate; 1 = Fishery; 2 = Survey}
#' \item{Species}{Species number}
#' \item{Selectivity_index}{index to use if selectivities of different surveys are to be the same}
#' \item{Selectivity}{Selectivity to use for the species: 0 = empirical selectivity provided in srv_emp_sel; 1 = logistic selectivity; 2 = non-parametric selecitivty sensu Ianelli et al 2018; 3 = double logistic; 4 = descending logistic}
#' \item{Nselages}{Number of ages to estimate non-parametric selectivity for Selectivity = 2. Not used otherwise}
#' \item{Time_varying_sel}{Wether a time-varying selectivity should be estimated for logistic, double logistic selectivity, or descending logistic. 0 = no, 1 = penalized deviates given sel_sd_prior, 2 = random effect, 3 = time blocks with no penality, 4 = random walk following Dorn, 5 = random walk on ascending portion of double logistic only. If selectivity is set to type = 2 (non-parametric) this value will be the 1st penalty on selectivity.}
#' \item{Sel_sd_prior}{The sd to use for the random walk of time varying selectivity if set to 1. If selectivity is set to type = 2 (non-parametric) this value will be the 2nd penalty on selectivity.}
#' \item{Age_first_selected}{Age at which selectivity is non-zero}
#' \item{Acuumulation_age_lower}{Ages below this will be grouped to this age for composition data. For example, if set to 2, comp data for age 2 will include 1 and 2 year olds.}
#' \item{Acuumulation_age_upper}{Ages above this will be grouped to this age for composition data. For example, if set to 9 for a species with 10 ages, comp data for age 9 will include 9 and 10 year olds.}
#' \item{weight1_Numbers2}{Is the observation in weight (kg) set as 1, if the observation is in numbers caught, set as 2}
#' \item{Weight_index}{Weight-at-age (wt) index to use for calculation of derived quantities}
#' \item{Age_transition_index}{Age transition matrix (e.g. growth trajectory) index to use for derived quantities to convert age to length}
#' \item{Q_index}{index to use if catchability coefficients are to be set the same}
#' \item{Estimate_q}{Estimate catchability? (0 = fixed at prior; - 1 = Estimate single parameter; - 2 = Estimate single parameter with prior; - 3 = Estimate analytical q  from Ludwig and Walters 1994;  - 4 = Estimate power equation; - 5 - Linear equation log(q_y) = q_mu + beta * index_y)}
#' \item{Q_prior}{Starting value or fixed value for catchability}
#' \item{Q_sd_prior}{Variance of q prior: dnorm (log_q, log_q_prior, q_sd_prior)}
#' \item{Time_varying_q}{Wether a time-varying q should be estimated. 0 = no, 1 = penalized deviate, 2 = random effect, 3 = time blocks with no penalty; 4 = random walk from mean following Dorn 2018 (dnorm(q_y - q_y-1, 0, sigma). If Estimate_q = 5, this determines the environmental index to be used in the equation log(q_y) = q_mu + beta * index_y}
#' \item{Time_varying_q_sd_prior}{The sd to use for the random walk of time varying q if set to 1}
#' \item{Estimate_survey_sd}{Estimate survey variance (0 = use CV from srv_biom, 1 = yes, 2 = analytically estimate following (Ludwig and Walters 1994)}
#' \item{Survey_sd_prior}{Starting value to be used if Estimate_sigma_index = 1}
#' \item{Estimate_catch_sd}{Estimate fishery variance (0 = use CV from srv_biom, 1 = yes, 2 = analytically estimate following (Ludwig and Walters 1994)}
#' \item{Catch_sd_prior}{Starting value to be used if Estimate_sigma_catch = 1}
#' \item{Comp_weights}{Composition weights to be used for multinomial likelihood. These are multiplied. After running model, these will update to McAllister & Ianelli 1997 weights using the harmonic mean.}
#' \item{Catch_units}{Units used for survey: 1 = kg; 2 = numbers}
#' \item{proj_F_prop}{The proportion of future fishing mortality assigned to this fleet}
#' \item{Sex}{sex codes: 0=combined; 1=use female only; 2=use male only; 3 = joint female and male}
#'}
"BS2017MS"


