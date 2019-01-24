#' Data inputs for single species CEATTLE of the Bering Sea from 1979 to 2017
#'
#' A dataset containing the inputs used for CEATTLE
#'
#' @format A list of different integers and vectors, specifying model switches:
#' MODEL CONFIGURATION
#' \describe{
#' \item{debug}{Integer = Logical to debug (1) or not (0)}
#' \item{msmMode}{Integer: Single species (0) or multi-species (1) mode}
#' \item{est_diet}{Integer: include diet data in likelihood. 0 = no, 1 = yes}
#' \item{suitMode}{Integer: Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = gamma selectivity from Kinzey and Punt (2009)}
#' \item{avgnMode}{Integer: Average numbers-at-age to be used in predation function. 0 = average numbers-at-age, 1 = \eqn{= N * e^{-Z/2} }, 2 = \eqn{N}.
#' \item{random_rec}{Integer: Switch to estimate recruitment deviations as random effects. 0 = no, 1 = yes}
#' \item{niter}{Integer: Number of loops around population/predation equations}
#' \item{nspp}{Integer: number of species in model}
#' \item{styr}{Integer: Start year}
#' \item{nyrs}{Integer: Number of estimation years}
#' \item{logist_sel_phase}{Ivector: selectivity type. 0 = logistic, 1 = non-parametric selectivity estimated for each age; length = nspp:}
#' \item{nselages}{Ivector: if non-parametric selectivity, the number of ages to estimate selectivity; length = nspp:}
#' }
#'
#'
#'
#'
#' DATA INPUTS: A list of different integers, vectors, matrices, and arrays. NOTE: Matrices and arrays will be filled with NA to the largest size. For example if nages of species 1 is 12 and nages of species 2 is 24}
#'
#' \describe{
#' \item{nages}{Ivector: Number of species ages; length = nspp}
#' \item{1. FISHERY COMPONENTS}{}
#' \item{nyrs_tc_biom}{Vector: Number of years with total observed catch; length = nspp}
#' \item{yrs_tc_biom}{Matrix: Years with total observed catch; dim = [nspp, nyrs_tc_biom]}
#' \item{tcb_obs}{Matrix: Observed total yield (kg}{ dim = [nspp, nyrs_tc_biom]}
#' \item{ nyrs_fsh_comp }{Ivector: Number of years in the fishery sp_age composition data; dim = [nspp]}
#' \item{ yrs_fsh_comp }{Imatrix: Years for the fishery sp_age composition data; dim = [nspp, nyrs_fsh_comp]}
#' \item{ fsh_age_type }{Ivector: Composition type for fishery catch (1 = age based, 2 = length based}{ dim = [nspp]}
#' \item{ fsh_age_bins }{Ivector: Number of bins for fishery age/length composition data; dim = [nspp]}
#' \item{ obs_catch }{Array: Observed fishery catch-at-age or catch-at-length; dim = [nspp, fsh_age_bins, nyrs_fsh_comp]}
#' \item{ 2. SURVEY COMPONENTS}
#' \item{ nyrs_srv_biom }{Ivector: Number of years of survey biomass data; dim = [nspp]}
#' \item{ yrs_srv_biom }{Imatrix: Years of survey biomass data; dim = [nspp, nyrs_srv_biom]}
#' \item{ srv_biom }{Matrix: Observed BT survey biomass (kg}{ dim = [nspp, nyrs]}
#' \item{ srv_biom_se }{Matrix: Observed annual biomass error (SE}{ dim = [nspp, nyrs_srv_biom]}
#' \item{ nyrs_srv_age }{Ivector: Number of years of survey age/length composition; length = [nspp]}
#' \item{ yrs_srv_age }{Imatrix: Years for the survey age/length composition data; dim = [nspp, nyrs_srv_age]}
#' \item{ srv_age_type }{Ivector: Type of compisition (1 = age; 2 = length}{ length = [nspp]}
#' \item{ srv_age_bins }{Ivector: Number of size binds for the age/length comps; dim = [nspp]}
#' \item{ srv_age_n }{Matrix: Sample size for the multinomial; dim = [nspp, nyrs_srv_age]}
#' \item{ srv_age_obs }{Array: Observed BT age comp; dim = [nspp, nages, nyrs]}
#' \item{ age_trans_matrix}{Array: observed sp_age/size compositions; dim = [nspp, nages, srv_age_bins]}
#' \item{ n_eit}{Integer: Number of years with EIT data; dim = [1]}
#' \item{ yrs_eit}{Ivector: Years for available EIT data; dim = [n_eit]}
#' \item{ eit_age_n}{Vector: Numberof EIT Hauls for multinomial; dim = [yrs_eit]}
#' \item{ obs_eit_age}{Matrix: Observed EIT catch-at-age; dim = [n_eit, nyrs]}
#' \item{ obs_eit}{Vector: Observed EIT survey biomass (kg}{ length = [nyrs]}
#' \item{ eit_sel}{Matrix: Observed EIT survey selectivity; dim = [eit_age, nyrs_eit_sel]}
#' \item{3. RATION AND CONSUMPTION PARAMETERS}{}
#' \item{ fday }{Vector: number of foraging days for each predator; length = [nspp]}
#' \item{ Pyrs }{Array: : dim = [nspp, nyrs+1, nages]}
#' \item{ Uobs }{Array: pred, prey, predL, preyL U matrix (mean number of prey in each pred}{ dim = [nspp, nspp, maxL, maxL]}
#' \item{ UobsWt }{Array: pred, prey, predL, preyL U matrix (mean wt_hat of prey in each pred}{ dim = [nspp, nspp, maxL, maxL] }
#' \item{ UobsAge }{Array: pred, prey, predA, preyA U matrix (mean number of prey in each pred age}{ dim = [nspp, nspp, max_age, max_age]}
#' \item{ Mn_LatAge }{Matrix: Mean length-at-age; dim = [nspp, nages]}
#' \item{ nTyrs }{Integer: Number of temperature years; dim = [1] }
#' \item{ Tyrs }{Ivector: Years of hindcast data; length = [nTyrs]}
#' \item{ BTempC_retro }{Vector: bottom temperature; length = [nTyrs ]}
#' \item{ other_food }{Vector: Biomass of other prey (kg}{ length = [nspp]}
#' \item{ C_model }{Ivector: f == 1, the use Cmax*fT*P; length = [nspp]}
#' \item{ Pvalue }{Vector: This scales the pvalue used if C_model ==1 , proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge; length = [nspp]}
#' \item{ Ceq }{Ivector: Ceq: which Comsumption equation to use; length = [nspp]; Currently all sp = 1}
#' \item{ CA }{Vector: Wt specific intercept of Cmax=CA*W^CB; length = [nspp]}
#' \item{ CB }{Vector: Wt specific slope of Cmax=CA*W^CB; length = [nspp]}
#' \item{ Qc }{Vector: used in fT, QC value; length = [nspp]}
#' \item{ Tco }{Vector: used in fT, thermal optimum; length = [nspp]}
#' \item{ Tcm }{Vector: used in fT, thermal max; length = [nspp]}
#' \item{ Tcl }{Vector: used in fT eq 3, limit; length = [nspp]}
#' \item{ CK1 }{Vector: used in fT eq 3, limit where C is .98 max (ascending}{ length = [nspp]}
#' \item{ CK4 }{Vector: used in fT eq 3, temp where C is .98 max (descending}{ length = [nspp]}
#' \item{ S_a }{Matrix: S_a, S_b, S_b2, S_b3, S_b4, S_b5: a,L,L^2,L^3,L^4,L^5 (rows)coef for mean S=a+b*L+b2*L*L, whith a cap at 80cm for each pred spp(cols}{ dim = [6, nspp]}
#' \item{3. WEIGHT, MATURITY, MORTALITY}{}
#' \item{ wt }{Array: Weight-at-age by year; dim = [nyrs, nages, nspp]}
#' \item{ aLW }{Matrix: LW a&b regression coefs for W=a*L^b; dim = [2, nspp]}
#' \item{ M1_base }{Matrix: Residual natural mortality; dim = [nspp, nages]}
#' \item{ mf_type }{Ivector: Sex specific mort and weight at age? : 1 = same for both, 2 = seperate wt at sp_age for each sex; length = nspp}
#' \item{ propMorF }{Matrix: Proportion-at-age of females of population; dim = [nspp, nages]}
#' \item{ pmature }{Matrix: Proportion of mature females at age; dim = [nspp, nages]}
#' ...
#' }
"BS2017SS"


#' Data inputs for multi-species CEATTLE of the Bering Sea from 1979 to 2017
#'
#' A dataset containing the inputs used for CEATTLE
#'
#' @format A list of different integers and vectors, specifying model switches:
#' MODEL CONFIGURATION
#' \describe{
#' \item{debug}{Integer = Logical to debug (1) or not (0)}
#' \item{msmMode}{Integer: Single species (0) or multi-species (1) mode}
#' \item{est_diet}{Integer: include diet data in likelihood. 0 = no, 1 = yes}
#' \item{suitMode}{Integer: Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = gamma selectivity from Kinzey and Punt (2009)}
#' \item{avgnMode}{Integer: Average numbers-at-age to be used in predation function. 0 = average numbers-at-age, 1 = \eqn{= N * e^{-Z/2} }, 2 = \eqn{N}.
#' \item{random_rec}{Integer: Switch to estimate recruitment deviations as random effects. 0 = no, 1 = yes}
#' \item{niter}{Integer: Number of loops around population/predation equations}
#' \item{nspp}{Integer: number of species in model}
#' \item{styr}{Integer: Start year}
#' \item{nyrs}{Integer: Number of estimation years}
#' \item{logist_sel_phase}{Ivector: selectivity type. 0 = logistic, 1 = non-parametric selectivity estimated for each age; length = nspp:}
#' \item{nselages}{Ivector: if non-parametric selectivity, the number of ages to estimate selectivity; length = nspp:}
#' }
#'
#'
#'
#'
#' DATA INPUTS: A list of different integers, vectors, matrices, and arrays. NOTE: Matrices and arrays will be filled with NA to the largest size. For example if nages of species 1 is 12 and nages of species 2 is 24}
#'
#' \describe{
#' \item{nages}{Ivector: Number of species ages; length = nspp}
#' \item{1. FISHERY COMPONENTS}{}
#' \item{nyrs_tc_biom}{Vector: Number of years with total observed catch; length = nspp}
#' \item{yrs_tc_biom}{Matrix: Years with total observed catch; dim = [nspp, nyrs_tc_biom]}
#' \item{tcb_obs}{Matrix: Observed total yield (kg}{ dim = [nspp, nyrs_tc_biom]}
#' \item{ nyrs_fsh_comp }{Ivector: Number of years in the fishery sp_age composition data; dim = [nspp]}
#' \item{ yrs_fsh_comp }{Imatrix: Years for the fishery sp_age composition data; dim = [nspp, nyrs_fsh_comp]}
#' \item{ fsh_age_type }{Ivector: Composition type for fishery catch (1 = age based, 2 = length based}{ dim = [nspp]}
#' \item{ fsh_age_bins }{Ivector: Number of bins for fishery age/length composition data; dim = [nspp]}
#' \item{ obs_catch }{Array: Observed fishery catch-at-age or catch-at-length; dim = [nspp, fsh_age_bins, nyrs_fsh_comp]}
#' \item{ 2. SURVEY COMPONENTS}
#' \item{ nyrs_srv_biom }{Ivector: Number of years of survey biomass data; dim = [nspp]}
#' \item{ yrs_srv_biom }{Imatrix: Years of survey biomass data; dim = [nspp, nyrs_srv_biom]}
#' \item{ srv_biom }{Matrix: Observed BT survey biomass (kg}{ dim = [nspp, nyrs]}
#' \item{ srv_biom_se }{Matrix: Observed annual biomass error (SE}{ dim = [nspp, nyrs_srv_biom]}
#' \item{ nyrs_srv_age }{Ivector: Number of years of survey age/length composition; length = [nspp]}
#' \item{ yrs_srv_age }{Imatrix: Years for the survey age/length composition data; dim = [nspp, nyrs_srv_age]}
#' \item{ srv_age_type }{Ivector: Type of compisition (1 = age; 2 = length}{ length = [nspp]}
#' \item{ srv_age_bins }{Ivector: Number of size binds for the age/length comps; dim = [nspp]}
#' \item{ srv_age_n }{Matrix: Sample size for the multinomial; dim = [nspp, nyrs_srv_age]}
#' \item{ srv_age_obs }{Array: Observed BT age comp; dim = [nspp, nages, nyrs]}
#' \item{ age_trans_matrix}{Array: observed sp_age/size compositions; dim = [nspp, nages, srv_age_bins]}
#' \item{ n_eit}{Integer: Number of years with EIT data; dim = [1]}
#' \item{ yrs_eit}{Ivector: Years for available EIT data; dim = [n_eit]}
#' \item{ eit_age_n}{Vector: Numberof EIT Hauls for multinomial; dim = [yrs_eit]}
#' \item{ obs_eit_age}{Matrix: Observed EIT catch-at-age; dim = [n_eit, nyrs]}
#' \item{ obs_eit}{Vector: Observed EIT survey biomass (kg}{ length = [nyrs]}
#' \item{ eit_sel}{Matrix: Observed EIT survey selectivity; dim = [eit_age, nyrs_eit_sel]}
#' \item{3. RATION AND CONSUMPTION PARAMETERS}{}
#' \item{ fday }{Vector: number of foraging days for each predator; length = [nspp]}
#' \item{ Pyrs }{Array: : dim = [nspp, nyrs+1, nages]}
#' \item{ Uobs }{Array: pred, prey, predL, preyL U matrix (mean number of prey in each pred}{ dim = [nspp, nspp, maxL, maxL]}
#' \item{ UobsWt }{Array: pred, prey, predL, preyL U matrix (mean wt_hat of prey in each pred}{ dim = [nspp, nspp, maxL, maxL] }
#' \item{ UobsAge }{Array: pred, prey, predA, preyA U matrix (mean number of prey in each pred age}{ dim = [nspp, nspp, max_age, max_age]}
#' \item{ Mn_LatAge }{Matrix: Mean length-at-age; dim = [nspp, nages]}
#' \item{ nTyrs }{Integer: Number of temperature years; dim = [1] }
#' \item{ Tyrs }{Ivector: Years of hindcast data; length = [nTyrs]}
#' \item{ BTempC_retro }{Vector: bottom temperature; length = [nTyrs ]}
#' \item{ other_food }{Vector: Biomass of other prey (kg}{ length = [nspp]}
#' \item{ C_model }{Ivector: f == 1, the use Cmax*fT*P; length = [nspp]}
#' \item{ Pvalue }{Vector: This scales the pvalue used if C_model ==1 , proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge; length = [nspp]}
#' \item{ Ceq }{Ivector: Ceq: which Comsumption equation to use; length = [nspp]; Currently all sp = 1}
#' \item{ CA }{Vector: Wt specific intercept of Cmax=CA*W^CB; length = [nspp]}
#' \item{ CB }{Vector: Wt specific slope of Cmax=CA*W^CB; length = [nspp]}
#' \item{ Qc }{Vector: used in fT, QC value; length = [nspp]}
#' \item{ Tco }{Vector: used in fT, thermal optimum; length = [nspp]}
#' \item{ Tcm }{Vector: used in fT, thermal max; length = [nspp]}
#' \item{ Tcl }{Vector: used in fT eq 3, limit; length = [nspp]}
#' \item{ CK1 }{Vector: used in fT eq 3, limit where C is .98 max (ascending}{ length = [nspp]}
#' \item{ CK4 }{Vector: used in fT eq 3, temp where C is .98 max (descending}{ length = [nspp]}
#' \item{ S_a }{Matrix: S_a, S_b, S_b2, S_b3, S_b4, S_b5: a,L,L^2,L^3,L^4,L^5 (rows)coef for mean S=a+b*L+b2*L*L, whith a cap at 80cm for each pred spp(cols}{ dim = [6, nspp]}
#' \item{3. WEIGHT, MATURITY, MORTALITY}{}
#' \item{ wt }{Array: Weight-at-age by year; dim = [nyrs, nages, nspp]}
#' \item{ aLW }{Matrix: LW a&b regression coefs for W=a*L^b; dim = [2, nspp]}
#' \item{ M1_base }{Matrix: Residual natural mortality; dim = [nspp, nages]}
#' \item{ mf_type }{Ivector: Sex specific mort and weight at age? : 1 = same for both, 2 = seperate wt at sp_age for each sex; length = nspp}
#' \item{ propMorF }{Matrix: Proportion-at-age of females of population; dim = [nspp, nages]}
#' \item{ pmature }{Matrix: Proportion of mature females at age; dim = [nspp, nages]}
#' ...
#' }
"BS2017MS"


