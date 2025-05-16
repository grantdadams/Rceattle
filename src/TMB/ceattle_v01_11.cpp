#include <TMB.hpp>
// ------------------------------------------------------------------------- //
//                 CEATTLE version 3.1.2                                     //
//                  Template Model Builder                                   //
//               Multispecies Statistical Model                              //
//          Bioenergetic-based Assessment for Understanding                  //
//              Biomass Linkages To The Environment                          //
// AUTHORS:   Grant Adams                                                    //
// CITATIONS:                                                                //
// 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
// 2. Adams, G. D., Holsman, K. K., Barbeaux, S_at_age. J., Dorn, M_at_age. W., Ianelli, J. N., Spies, I., ... & Punt, A. E. (2022). An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fisheries Research, 251, 106303.
// ------------------------------------------------------------------------- //
//
// Best attempt to use google C++ style guide https://google.github.io/styleguide/cppguide.html#General_Naming_Rules
//
//  INDEX:
//  0. Load dependencies
//  1. Model configuration
//  2. Model inputs
//  3. Model parameters
//  4. Derived quantities
//  5. Initial calculations
//  6. Population dynamics
//  7. Ration equations
//  8. Predation mortality equations
//  -- 8.0. Suitability equations
//  -- 8.1. Holsman et al Predation mortality
//  -- 8.2. Kinzey and Punt predation mortality
//  9. Survey components
//  10. Fishery components
//  11. Likelihood components
//  12. Report section
//  13. Model return/end
//  14. Change log
//
//
// ------------------------------------------------------------------------- //
// 0. LOAD DEPENDENCIES                                                      //
// ------------------------------------------------------------------------- //

// Function for getting max of an IVECTOR and return an int
template <class Type>
Type imax(const vector<Type> &x)
{
  int res = x[0];
  for(int i=0; i < x.size(); i++){
    int res2 = x[i];
    if(res < res2){
      res = res2;
    }
  }
  return res;
}

// Function to maintain > limit
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);   // tune a new candidate x given eps
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * (x - eps) * (x - eps));  // penalize via x vs eps
  return ans;
}


// Dirichlet multinomial (from WHAM)
template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type phi=sum(alpha);
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + alpha(a)) - lgamma(alpha(a));
  if(do_log == 1) return ll;
  else return exp(ll);
}

// From TMB examples, modified name, constrains between -1 and 1
template <class Type>
Type rho_trans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}


// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


// Function for detecting Inf
template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}


//  Function to getting sqaure
template <class Type> Type square(Type x){return x*x;}


// Function for mean of a vector
template <class Type>
Type mean(const vector<Type> &x)
{
  Type mean = x.sum()/x.size();
  return mean;
}

// Function to calculate the first difference
template <class Type>
vector<Type> first_difference(const vector<Type> &x)
{
  int length = x.size() - 1;

  if(length <= 0){
    std::cerr << "Error -- vector size too small in first_difference(const dvector&)" << std::endl;
    return(0);
  }

  vector<Type> tmp(length);

  for(int i = 0; i < length; i++){
    tmp(i) = x(i+1) - x(i);
  }

  return tmp;
}


// Function to get max of two values (differentiable)
// https://groups.google.com/g/tmb-users/c/QAVmgj66OC0
template <class Type>
Type max2(Type x, Type y){
  Type ans = 0.5*(CppAD::abs(x-y)+x+y);
  return ans;
}


template<class Type>
Type objective_function<Type>::operator() () {
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  // ------------------------------------------------------------------------- //
  // 1. MODEL CONFIGURATION                                                    //
  // ------------------------------------------------------------------------- //
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(estimateMode);             // Logical to debug or not
  DATA_INTEGER(msmMode);
  //    0 = run in single species mode
  //    1 = run in Type II MSVPA based (sensu Holsman et al (2015))
  //    2 = run in Type III MSVPA based
  //    3 = Holling Type I (linear)
  //    4 = Holling Type II
  //    5 = Holling Type III
  //    6 = Predator interference
  //    7 = Predator preemption
  //    8 = Hassell-Varley
  //    9 = Ecosim
  // DATA_INTEGER(est_diet);              // Include diet data in the likelihood
  DATA_IVECTOR(suitMode);                 // Estimate suitability
  // DATA_INTEGER(avgnMode);              // N used for predation function
  //    0 = avgN_at_age
  //    1 = N*exp(-Z / 2))
  //    2 = N
  DATA_INTEGER(initMode);                 // How the age-structure is initialized
  DATA_IVECTOR(forecast);                 // Run the model in forecast mo


  // 1.2. Temporal dimensions
  DATA_INTEGER( styr );                   // Start year
  DATA_INTEGER( endyr );                  // End of estimation years
  DATA_INTEGER( projyr );                 // End year of projection

  DATA_INTEGER( srr_meanyr );             // The last year used to calculate average recruitment. Used for MSE runs.
  DATA_INTEGER( suit_styr );              // The first year used to calculate suitability averages.
  DATA_INTEGER( suit_endyr );             // The last year used to calculate suitability averages.
  DATA_INTEGER( srr_hat_styr );           // The first year used to calculate stock-recuitment penalties or env-rec relationship.
  DATA_INTEGER( srr_hat_endyr );          // The last year used to calculate stock-recuitment penalties or env-rec relationship.

  int nyrs = projyr - styr + 1;
  int nyrs_hind = endyr - styr + 1;

  suit_endyr = suit_endyr - styr;
  suit_styr = suit_styr - styr;
  int nyrs_suit = suit_endyr - suit_styr + 1;
  int nyrs_srrmean = srr_meanyr - styr + 1;

  srr_hat_styr = srr_hat_styr - styr;
  srr_hat_endyr = srr_hat_endyr - styr;
  if(nyrs_srrmean > nyrs_hind){nyrs_srrmean = nyrs_hind;}
  if(srr_hat_styr == 0){srr_hat_styr = 1;} // R_hat starts at year 1.

  // 1.3. Number of species
  DATA_INTEGER( nspp );                   // Number of species (prey)
  DATA_IVECTOR( pop_wt_index );           // Dim 3 of weight to use for population dynamics
  DATA_IVECTOR( ssb_wt_index );           // Dim 3 of weight to use for spawning stock biomass calculation
  DATA_IVECTOR( pop_age_transition_index );// Dim 3 of weight to use for age_transition_matrix
  DATA_IVECTOR( estDynamics );            // Index indicating wether the population parameters are estimated (0), numbers-at-age are provided (1), or an index of numbers-at-age multiplied by an estimated scalar is used (2)
  pop_wt_index -= 1;                      // Indexing starts at 0
  ssb_wt_index -= 1;                      // Indexing starts at 0
  pop_age_transition_index -= 1;          // Indexing starts at 0


  // 1.4. RECRUITMENT SETTINGS
  DATA_INTEGER(srr_fun);                  // Stock recruit relationship for hindcast estimation
  DATA_INTEGER(srr_pred_fun);             // Stock recruit relationship for projection/brps/penalty
  DATA_INTEGER(proj_mean_rec);
  //    0 = project recruitment using ln_R0 and rec devs
  //    1 = project recruitment using mean rec (can also have adjusted rec devs)
  DATA_INTEGER(srr_est_mode);             // Logical of wether to add normal prior to stock recruit-relationship
  DATA_VECTOR(srr_prior);                 // Prior mean for stock recruit relationship parameter
  DATA_VECTOR(srr_prior_sd);              // Prior sd for stock recruit relationship parameter
  DATA_INTEGER( niter );                  // Number of loops for MSM mode
  DATA_VECTOR( Bmsy_lim );                // Upper limit for Bmsy in ricker function. Will add penalty if 1/beta > lim


  // 1.5. HARVEST CONTROL RULE (HCR) SETTINGS
  DATA_INTEGER(HCR);                      // Function to be used for harvest control rule
  DATA_INTEGER(DynamicHCR);               // TRUE/FALSE. Wether to use static or dynamic reference points (default = FALSE)
  DATA_VECTOR(Ptarget);                   // Target spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0
  DATA_VECTOR(Plimit);                    // Limit spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0
  DATA_VECTOR(Ftarget_percent);           // Percentage of spawning-stock biomass per recruit or SB0 at F = 0 used to find the target F
  DATA_VECTOR(Flimit_percent);            // Percentage of spawning-stock biomass per recruit or SB0 at F = 0 used to find the limit F
  DATA_VECTOR(Alpha);                     // Parameter used in NPFMC Tier 3 HCR
  DATA_VECTOR(Fmult);                     // Multiplier for target fishing mortality (Fmult * Ftarget). For example 0.75 * F40%.
  DATA_VECTOR(QnormHCR);                  // Add to Flimit to set Ftarget based on the Pstar approach: Flimit + qnorm(Pstar, 0, Sigma)


  // 1.6. MODEL OBJECTS
  // 1.5.1. LOOPING INDICES -- k = observation, sp = species, sex = sex (0 = combined; 1 = females; 2 = males)
  // age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length
  int sp, sex, age, ln, ksp, k_sex, k_age, k_ln, yr, rsp, r_sex, r_age, r_ln; // k
  int index, flt;                                                         // Survey and fishery indices
  int flt_yr, flt_sex, comp_type;
  int flt_ind, fsh_ind, index_ind, comp_ind, yr_ind;                      // Indices for survey sets
  int sel_type, nselages;
  Type mo = 0;                                                            // Month float
  if(msmMode == 0) { niter = 1; }                                        // Number of iterations for SS mode


  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                           //
  // ------------------------------------------------------------------------- //

  // -- 2.1. Species attributes
  DATA_IVECTOR( nsex );                   // Number of sexes to be modelled; 1 = sexes combined/single sex, 2 = 2 sexes
  int max_sex = imax(nsex);               // Integer of maximum sexes to make the arrays
  DATA_VECTOR( spawn_month );             // Month of spawning to adjust mortality
  DATA_IVECTOR( nages );                  // Number of species (prey) ages
  DATA_IVECTOR( minage );                 // Minimum age of each species
  DATA_IVECTOR( nlengths );               // Number of species (prey) lengths
  DATA_ARRAY( NByageFixed );              // Provided estimates of numbers- or index-at-age to be multiplied (or not) by pop_scalar to get N_at_age
  int max_age = imax(nages);              // Integer of maximum nages to make the arrays
  DATA_VECTOR( MSSB0 );                   // SB0 from projecting the model forward in multi-species mode under no fishing
  DATA_VECTOR( MSB0 );                    // B0 from projecting the model forward in multi-species mode under no fishing

  // -- 2.2. M1_at_age specifications
  DATA_IVECTOR(M1_model);
  DATA_IVECTOR(M1_re);
  DATA_IVECTOR(M1_use_prior);
  DATA_IVECTOR(M2_use_prior);
  DATA_VECTOR(M_prior);
  DATA_VECTOR(M_prior_sd);

  // -- 2.3. Data controls (i.e. how to assign data to objects)
  DATA_IMATRIX( fleet_control );          // Fleet specifications

  // -- 2.3.1 Fishery Components
  DATA_IMATRIX( catch_ctl );              // Info for fishery biomass; columns = Fishery_name, Fishery_code, Species, Year
  DATA_IMATRIX( catch_n );                // Info for fishery biomass; columns = Month
  DATA_MATRIX( catch_obs ); DATA_UPDATE( catch_obs ) // Observed fishery catch biomass (kg) and log_sd; n = [nobs_fish_biom, 2]; columns = Observation, Error

    // -- 2.3.2 Survey components
    DATA_IMATRIX( index_ctl );            // Info for survey biomass; columns = Survey_name, Survey_code, Species, Year
  DATA_MATRIX( index_n );                 // Info for survey biomass; columns = Month
  DATA_MATRIX( index_obs );               // Observed survey biomass (kg) and log_sd; columns = Observation, Error
  DATA_VECTOR( index_ln_q_prior );        // Prior mean for survey catchability

  // -- 2.3.3. Composition data
  DATA_IMATRIX( comp_ctl );               // Info on observed survey age/length comp; columns = Survey_name, Survey_code, Species, Year
  DATA_MATRIX( comp_n );                  // Month and sample size on observed survey age/length comp; columns = Month, Sample size
  DATA_MATRIX( comp_obs );                // Observed survey age/length comp; cols = Comp_1, Comp_2, etc. can be proportion

  // -- 2.3.4 Age and selectivity
  DATA_IMATRIX( emp_sel_ctl );            // Info on empirical fishery selectivity; columns =  Fishery_name, Fishery_code, Species, Year
  DATA_MATRIX( emp_sel_obs );             // Observed emprical fishery selectivity; columns = Compe_1, Comp_2, etc.
  DATA_ARRAY( age_trans_matrix);          // observed sp_age/size compositions; n = [nspp, nages, index_age_bins]
  DATA_ARRAY( age_error );                // Array of aging error matrices for each species; n = [nspp, nages, nages]

  // -- 2.3.5. Growth
  DATA_ARRAY( weight );                   // Weight-at-age by year; n = [nweight, sex, nages, nyrs]
  // DATA_ARRAY( laa );                      // Length-at-age by year; n = [nspp, sex, nages, nyrs]

  // 2.3.6. Diet data
  DATA_VECTOR( fday );                    // number of foraging days for each predator
  DATA_ARRAY( Pyrs );                     // Relative-foraging rate
  DATA_MATRIX( diet_obs );                // pred, prey, predA, preyA U observations (mean wt_hat of prey in each pred age)
  DATA_IMATRIX( diet_ctl );               // Info on pred, prey, predA, preyA U matrix (mean wt_hat of prey in each pred age)

  // 2.3.7. Environmental data
  DATA_MATRIX( env_index );               // Matrix of environmental predictors such as bottom temperature

  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( other_food );              // Biomass of other prey (kg)
  DATA_VECTOR( Pvalue );                  // This scales the pvalue used, proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge
  DATA_IVECTOR( Ceq );                    // Ceq: which Comsumption equation to use; Currently all sp = 1
  DATA_IVECTOR( Cindex );                 // Cindex, which environmental index in env_index should drive bioenergetics.
  DATA_VECTOR( CA );                      // Wt specific intercept of Cmax=CA*W^CB
  DATA_VECTOR( CB );                      // Wt specific slope of Cmax=CA*W^CB
  DATA_VECTOR( Qc );                      // used in fT, QC value
  DATA_VECTOR( Tco );                     // used in fT, thermal optimum
  DATA_VECTOR( Tcm );                     // used in fT, thermal max
  DATA_VECTOR( Tcl );                     // used in fT eq 3, limit
  DATA_VECTOR( CK1 );                     // used in fT eq 3, limit where C is .98 max (ascending)
  DATA_VECTOR( CK4 );                     // used in fT eq 3, temp where C is .98 max (descending)

  // -- 2.4.3. Others
  DATA_MATRIX( sex_ratio );               // Proportion-at-age of females of population
  DATA_MATRIX( maturity );                // Proportion of mature females at age; [nspp, nages]


  // ------------------------------------------------------------------------- //
  // 3. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //

  PARAMETER( dummy );                             // Variable to test derived quantities given input parameters; n = [1]
  PARAMETER_MATRIX( ln_pop_scalar );              // Scalar to multiply supplied numbers at age by

  // -- 3.1. Recruitment parameters
  PARAMETER_MATRIX( rec_pars );                   // Stock-recruit parameters: col1 = mean rec, col2 = SRR alpha, col3 = SRR beta
  PARAMETER_MATRIX( beta_rec_pars );              // Regression parameters for environmental linkage to stock-recruit function row is spp
  PARAMETER_VECTOR( R_ln_sd );                    // Standard deviation of recruitment deviations
  PARAMETER_MATRIX( rec_dev );                    // Annual recruitment deviation; n = [nspp, nyrs]
  PARAMETER_MATRIX( init_dev );                   // Initial abundance-at-age # NOTE: Need to figure out how to best vectorize this

  // -- 3.2. Natural mortality (M1)
  PARAMETER_ARRAY( ln_M1 );                       // Natural mortality (residual if multispecies mode or total if single species mode); n = [nspp, nsex, nages]
  PARAMETER_ARRAY( ln_M1_dev );                   // Natural mortality annual deviate; n = [nspp, nsex, nyrs]
  PARAMETER_ARRAY( M1_beta );                     // Regression coefficients for environmnetally linked M1; n = [nspp, nsex, n indices]
  PARAMETER_ARRAY( M1_rho );                      // Correlation for AR1 random effects on age and year; n = [nspp, nsex, 2]
  PARAMETER_ARRAY( M1_dev_ln_sd );                // Standard deviation of random effects on age and year; n = [nspp, nsex, 2]

  // -- 3.3. Fishing mortality parameters
  PARAMETER_VECTOR( ln_Flimit );                  // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_VECTOR( ln_Ftarget );                 // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_VECTOR( ln_Finit );                   // Fishing mortality for initial population to induce non-equilibrium; n = [nspp]
  PARAMETER_VECTOR( proj_F_prop );                // Proportion of fishing mortality from each fleet for projections; n = [n_fsh]
  PARAMETER_MATRIX( ln_F );                       // Annual fishing mortality; n = [n_fsh, nyrs] # NOTE: The size of this will likely change

  // -- 3.4. Survey catchability parameters
  PARAMETER_VECTOR( index_ln_q );                 // Survey catchability; n = [n_index]
  PARAMETER_VECTOR( index_q_rho );                // Correlation parameter for AR1 on natural scale; n = [n_index]
  PARAMETER_MATRIX( index_q_beta );               // Survey catchability regression coefficient and rho parameters
  // PARAMETER_VECTOR( index_q_pow );             // Survey catchability power coefficient q * B ^ q_pow or beta ln(q_y) = q_mut + beta * index_y; n = [n_index]
  PARAMETER_MATRIX( index_q_dev );                // Annual survey catchability deviates; n = [n_index, nyrs_hind]
  PARAMETER_VECTOR( index_q_ln_sd );              // Log standard deviation of prior on survey catchability; n = [1, n_index]
  PARAMETER_VECTOR( index_q_dev_ln_sd );// Log standard deviation of time varying survey catchability; n = [1, n_index]

  // -- 3.5. Selectivity parameters
  PARAMETER_ARRAY( sel_coff );                    // selectivity parameters for non-parametric; n = [n_selectivities, nsex, nselages]
  PARAMETER_ARRAY( sel_coff_dev );                // Annual deviates for non-parametric selectivity parameters; n = [n_selectivities, nsex, nselages]
  PARAMETER_ARRAY( ln_sel_slp );                  // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( sel_inf );                     // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( ln_sel_slp_dev );              // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( sel_inf_dev );                 // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_VECTOR( sel_dev_ln_sd );               // Log standard deviation of selectivity; n = [1, n_selectivities]
  PARAMETER_MATRIX( sel_curve_pen );              // Selectivity penalty for non-parametric selectivity, 2nd column is for monotonic bit

  // -- 3.6. Data variance
  //FIXME: remove ln_sd terms below
  PARAMETER_VECTOR( index_ln_sd );                // Log standard deviation of survey index time-series; n = [1, n_index]
  PARAMETER_VECTOR( catch_ln_sd );                // Log standard deviation of fishery catch time-series; n = [1, n_fsh]
  PARAMETER_VECTOR( comp_weights );               // Weights for composition data
  vector<Type>  DM_pars = exp(comp_weights);      // Dirichlet-multinomial scalars

  PARAMETER_VECTOR( diet_comp_weights );          // Weights for diet composition data
  // vector<Type>  DM_diet_pars = exp(diet_comp_weights);// Dirichlet-multinomial scalars

  // -- 3.7. Kinzery predation function parameters
  /*
   PARAMETER_MATRIX(logH_1);                       // Predation functional form; n = [nspp, nspp2];
   PARAMETER_VECTOR(logH_1a);                      // Age adjustment to H_1; // FIXME: make matrix
   PARAMETER_VECTOR(logH_1b);                      // Age adjustment to H_1; // FIXME: make matrix
   PARAMETER_MATRIX(logH_2);                       // Predation functional form; n = [nspp, nspp]
   PARAMETER_MATRIX(logH_3);                       // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
   PARAMETER_MATRIX(H_4);                          // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;
   */

  // 3.8. Gamma selectivity parameters
  PARAMETER_VECTOR( log_gam_a );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  PARAMETER_VECTOR( log_gam_b );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10

  // 3.9. Preference
  PARAMETER_MATRIX( log_phi );                    // Species preference coefficient; n = [nspp, nspp]


  // ------------------------------------------------------------------------- //
  // 4. DERIVED QUANTITIES SECTION                                             //
  // ------------------------------------------------------------------------- //

  // 4.1. Derived indices
  // int max_bin = imax( nlengths );                                                   // Integer of maximum number of length/age bins.
  int n_flt = fleet_control.rows();
  vector<int> joint_adjust(comp_obs.rows()); joint_adjust.setZero();
  Type penalty = 0.0;
  Type ricker_intercept = 0.0;

  // -- 4.2. Estimated population quantities
  matrix<Type>  pop_scalar = ln_pop_scalar;  pop_scalar = exp(ln_pop_scalar.array());// Fixed n-at-age scaling coefficient
  vector<Type>  avg_R(nspp); avg_R.setZero();                                       // Mean recruitment of hindcast
  matrix<Type>  R_hat(nspp, nyrs); R_hat.setZero();                                 // Expected recruitment given SR curve
  matrix<Type>  mort_sum(nspp, max_age); mort_sum.setZero();
  vector<Type>  R0(nspp); R0.setZero();                                             // Equilibrium recruitment at F = 0.
  vector<Type>  R_init(nspp); R_init.setZero();                                     // Equilibrium recruitment at F = Finit (non-equilibrium).
  Type srr_alpha = 0.0;
  matrix<Type>  R(nspp, nyrs); R.setZero();                                         // Estimated recruitment (n)
  vector<Type>  steepness(nspp); steepness.setZero();                               // Expected % of R0 at 20% SSB0.
  array<Type>   biomass_at_age(nspp, max_sex, max_age, nyrs); biomass_at_age.setZero();// Estimated biomass-at-age (kg)
  matrix<Type>  biomass(nspp, nyrs); biomass.setZero();                             // Estimated biomass (kg)
  matrix<Type>  exploitable_biomass(nspp, nyrs); exploitable_biomass.setZero();     // Estimated exploitable biomass (kg)
  matrix<Type>  ssb(nspp, nyrs); ssb.setZero();                                     // Estimated spawning stock biomass (kg)
  matrix<Type>  biomass_depletion(nspp, nyrs); biomass_depletion.setZero();         // Estimated biomass biomass_depletion
  matrix<Type>  ssb_depletion(nspp, nyrs); ssb_depletion.setZero();                 // Estimated biomass_depletion of spawning stock biomass
  // array<Type>   ssb_at_age(nspp, max_age, nyrs); ssb_at_age.setZero();           // Spawning biomass at age (kg)
  array<Type>   M_at_age(nspp, max_sex, max_age, nyrs); M_at_age.setZero();               // Total natural mortality at age
  array<Type>   M1_at_age(nspp, max_sex, max_age, nyrs); M1_at_age.setZero();             // Residual or total natural mortality at age
  array<Type>   N_at_age(nspp, max_sex, max_age, nyrs); N_at_age.setZero();               // Numbers at age
  array<Type>   avgN_at_age(nspp, max_sex, max_age, nyrs); avgN_at_age.setZero();         // Average numbers-at-age
  // array<Type>   S_at_age(nspp, max_sex, max_age, nyrs); S_at_age.setZero();            // Survival at age
  array<Type>   Z_at_age(nspp, max_sex, max_age, nyrs); Z_at_age.setZero();               // Total mortality at age
  vector<Type>  R_sd(nspp); R_sd.setZero();                                         // Standard deviation of recruitment variation
  vector<Type>  zero_N_pen(nspp); zero_N_pen.setZero();                             // Additional penalty to add to likelihood if n-at-age goes < 0

  // -- 4.3. Selectivity parameters
  array<Type>   sel(n_flt, max_sex, max_age, nyrs); sel.setZero();                        // Estimated selectivity at age
  matrix<Type>  avg_sel(n_flt, max_sex); avg_sel.setZero();                               // Average selectivity
  array<Type>   non_par_sel(n_flt, max_sex, max_age); non_par_sel.setZero();              // Temporary saved selectivity at age for estimated bits
  vector<Type>  sel_dev_sd(n_flt); sel_dev_sd.setZero();                            // Standard deviation of selectivity deviates

  // -- 4.4. Fishery components
  matrix<Type>  F_spp(nspp, nyrs); F_spp.setZero();                                 // Fully selected fishing mortality by species
  matrix<Type>  F_flt(n_flt, nyrs); F_flt.setZero();                                // Fully selected fishing mortality by fleet
  array<Type>   F_flt_age(n_flt, max_sex, max_age, nyrs); F_flt_age.setZero();            // Estimated fishing mortality-at-age/sex for each fishery
  array<Type>   Flimit_age_spp(nspp, max_sex, max_age, nyrs); Flimit_age_spp.setZero();   // Estimated target fishing mortality-at-age/sex for each species
  array<Type>   Ftarget_age_spp(nspp, max_sex, max_age, nyrs); Ftarget_age_spp.setZero(); // Estimated limit fishing mortality-at-age/sex for each species
  array<Type>   F_spp_at_age(nspp, max_sex, max_age, nyrs); F_spp_at_age.setZero();       // Sum of annual estimated fishing mortalities for each species-at-age
  vector<Type>  catch_hat(catch_obs.rows()); catch_hat.setZero();                   // Estimated fishery yield/numbers (kg)
  vector<Type>  max_catch_hat(catch_obs.rows()); max_catch_hat.setZero();           // Estimated exploitable biomass/numbers by fleet (kg)
  vector<Type>  ln_catch_sd(catch_obs.rows()); ln_catch_sd.setZero();               // Estimated/fixed fishery log_sd (kg)

  // -- 4.5. Biological reference points
  array<Type>   NByage0(nspp, max_sex, max_age, nyrs); NByage0.setZero();                 // Numbers at age at mean recruitment and F = 0
  array<Type>   NByageF(nspp, max_sex, max_age, nyrs); NByageF.setZero();                 // Numbers at age at mean recruitment and F = Flimit
  array<Type>   DynamicNByage0(nspp, max_sex, max_age, nyrs); DynamicNByage0.setZero();   // Numbers at age at F = 0 (accounts for annual recruitment)
  array<Type>   DynamicNByageF(nspp, max_sex, max_age, nyrs); DynamicNByageF.setZero();   // Female numbers at age at F = Ftarget (accounts for annual recruitment)
  matrix<Type>  DynamicSB0(nspp, nyrs); DynamicSB0.setZero();                       // Estimated dynamic spawning biomass at F = 0 (accounts for S_at_age-R curve)
  matrix<Type>  DynamicB0(nspp, nyrs); DynamicB0.setZero();                         // Estimated dynamic  biomass at F = 0 (accounts for S_at_age-R curve)
  matrix<Type>  DynamicSBF(nspp, nyrs); DynamicSBF.setZero();                       // Estimated dynamic spawning biomass at F = Ftarget (accounts for S_at_age-R curve)
  array<Type>   NbyageSPR(4, nspp, max_age);                                        // Estimated numbers at age for spawning biomass per recruit reference points
  vector<Type>  SPRlimit(nspp); SPRlimit.setZero();                                 // Estimated Plimit SPR
  vector<Type>  SPRtarget(nspp); SPRtarget.setZero();                               // Estimated Ptarget SPR
  vector<Type>  SPR0(nspp); SPR0.setZero();                                         // Estimated spawning biomass per recruit at F = 0
  vector<Type>  SPRFinit(nspp); SPRFinit.setZero();                                 // Estimated spawning biomass per recruit at Finit
  matrix<Type>  SB0(nspp, nyrs); SB0.setZero();                                     // Estimated spawning stock biomass at F = 0 (Accounts for S_at_age-R)
  matrix<Type>  SBF(nspp, nyrs); SBF.setZero();                                     // Estimated spawning stock biomass at F = target (Accounts for S_at_age-R)
  matrix<Type>  B0(nspp, nyrs); B0.setZero();                                       // Estimated biomass at F = 0 (Accounts for S_at_age-R)
  vector<Type>  Flimit = exp(ln_Flimit);                                            // Target F parameter on natural scale
  vector<Type>  Ftarget = exp(ln_Ftarget);                                          // Limit F parameter on natural scale
  vector<Type>  Finit = exp(ln_Finit);                                              // Initial F for non-equilibrium age-structure
  vector<Type>  srr_mult(beta_rec_pars.cols()); srr_mult.setZero();                 // Environmental design matrix for rec
  vector<Type>  index_q_mult(index_q_beta.cols()); index_q_mult.setZero();          // Environmental design matrix for q
  vector<Type>  M1_mult(M1_beta.dim(2)); M1_mult.setZero();                         // Environmental design matrix for M1
  vector<Type>  beta_rec_tmp(beta_rec_pars.cols()); beta_rec_tmp.setZero();         // Temporary vector to store beta parameters by species for matrix mult
  vector<Type>  env_rec_tmp(beta_rec_pars.cols()); env_rec_tmp.setZero();           // Temporary vector to store env data by year for matrix mult
  vector<Type>  beta_q_tmp(index_q_beta.cols()); beta_q_tmp.setZero();              // Temporary vector to store Q beta parameters by species for matrix mult
  vector<Type>  env_q_tmp(index_q_beta.cols()); env_q_tmp.setZero();                // Temporary vector to store Q env data by year for matrix mult
  vector<Type>  beta_M1_tmp(M1_beta.dim(2)); beta_M1_tmp.setZero();                 // Temporary vector to store beta parameters by species for matrix mult
  vector<Type>  env_M1_tmp(M1_beta.dim(2)); env_M1_tmp.setZero();                   // Temporary vector to store env data by year for matrix mult
  matrix<Type>  proj_F(nspp, nyrs); proj_F.setZero();                               // Projected F (Fabc/Ftac/etc) using harvest control rule


  // -- 4.6. Survey components
  vector<Type>  index_q_sd(n_flt); index_q_sd.setZero();                            // Vector of standard deviation of survey catchability prior
  vector<Type>  index_q_dev_sd(n_flt); index_q_dev_sd.setZero();                    // Vector of standard deviation of time-varying survey catchability deviation
  Type avgsel_tmp = 0;                                                              // Temporary object for average selectivity across all ages
  vector<Type>  index_hat(index_obs.rows()); index_hat.setZero();                   // Estimated survey biomass (kg)
  vector<Type>  ln_index_sd(index_obs.rows()); ln_index_sd.setZero();               // Estimated/fixed log index sd (kg)
  vector<Type>  ln_index_analytical_sd(n_flt); ln_index_analytical_sd.setZero();    // Temporary vector to save analytical sd follow Ludwig and Walters 1994
  vector<Type>  index_q_analytical(n_flt); index_q_analytical.setZero();            // Temporary vector to save analytical sd follow Ludwig and Walters 1994
  matrix<Type>  index_q(n_flt, nyrs_hind); index_q.setZero();                       // Estimated survey catchability
  vector<Type>  index_n_obs(n_flt); index_n_obs.setZero();                          // Vector to save the number of observations for each survey time series

  // -- 4.7. Composition data - FIXME: will blow up in nlengths is less than nages
  vector<Type>  n_hat(comp_obs.rows()) ; n_hat.setZero() ;                          // Estimated catch (numbers)
  matrix<Type>  age_hat = comp_obs; age_hat.setZero();                              // Estimated catch at true age
  matrix<Type>  true_age_comp_hat = comp_obs; true_age_comp_hat.setZero();          // True estimated age composition
  matrix<Type>  age_obs_hat = comp_obs; age_obs_hat.setZero();                      // Estimated catch at observed age (accounts for ageing error)
  matrix<Type>  comp_hat = comp_obs; comp_hat.setZero();                            // Estimated comp

  // -- 4.8. Ration components
  array<Type>   consumption_at_age( nspp, max_sex, max_age, nyrs ); consumption_at_age.setZero();           // Pre-allocated indiviudal consumption in grams per predator-age
  matrix<Type>  fT( nspp, nyrs ); fT.setZero();                                     // Pre-allocation of temperature function of consumption
  // matrix<Type>  mnWt_obs( nspp, max_age ); mnWt_obs.setZero();                   // Mean observed weight at age (across years)
  array<Type>   ration( nspp, max_sex, max_age, nyrs ); ration.setZero();                 // Annual ration at age (kg/yr)

  // -- 4.9. Diet components
  array<Type>   diet_prop(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); diet_prop.setZero();             // Stomach proportion by weight U
  array<Type>   diet_prop_sum(nspp, max_sex, max_age, nyrs); diet_prop_sum.setZero();                               // Sum of stomach proportion for a predator
  array<Type>   diet_prop_hat(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); diet_prop_hat.setZero();     // Predicted stomach proportion by weight U
  array<Type>   other_food_diet_prop(nspp, max_sex, max_age, nyrs); other_food_diet_prop.setZero();                 // Other food diet proportion by weight
  matrix<Type>  diet_hat = diet_obs; diet_hat.setZero();                                                            // Estimated stomach proportion by weight U (formated following data input)

  // -- 4.10. Suitability components
  array<Type>   avail_food(nspp, max_sex, max_age, nyrs); avail_food.setZero();                               // Available food to predator
  array<Type>   stom_div_bio(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); stom_div_bio.setZero(); // Stomach proportion over biomass; U/ (W * N)
  array<Type>   suitability(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); suitability.setZero();   // Suitability/gamma selectivity of predator age u on prey age a
  array<Type>   suit_other(nspp, max_sex, max_age, nyrs); suit_other.setZero();                               // Suitability not accounted for by the included prey
  array<Type>   suma_suit(nspp, max_sex, max_age, nyrs); suma_suit.setZero();                                 // Sum of suitabilities

  // -- 4.11. Suitability parameters
  vector<Type> gam_a = exp(log_gam_a);                                    // Predator size-selectivity: shape parameter for gamma suitability, mean for normal of logs
  vector<Type> gam_b = exp(log_gam_b);                                    // Predator size-selectivity: scale parameter for gamma suitability, sd for normal of logs
  vector<Type> sum_phi(nspp); sum_phi.setZero();                          // Sum of predator-prey preference coefficients for multinomial transformation
  matrix<Type> vulnerability(nspp, nspp); vulnerability.setZero();        // Predator-prey preference coefficients
  vector<Type> vulnerability_other(nspp); vulnerability_other.setZero();  // Preference for other food


  // -- 4.12. Predation components
  array<Type>   M2_at_age(nspp, max_sex, max_age, nyrs); M2_at_age.setZero();                        // Predation mortality at age
  // array<Type>   M2_prop(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); M2_prop.setZero();     // Relative predation mortality at age from each species at age
  array<Type>   B_eaten(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); B_eaten.setZero();        // Biomass of prey eaten via predation by a predator at age
  array<Type>   B_eaten_as_prey(nspp, max_sex, max_age, nyrs); B_eaten_as_prey.setZero();            // Biomass eaten as prey via predation
  // array<Type>   B_eaten_as_pred(nspp, max_sex, max_age, nyrs); B_eaten_as_pred.setZero();         // Biomass eaten as predator via predation (used for Kinzey and Punt)
  // array<Type>   N_eaten(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); N_eaten.setZero();     // Number of prey of age a eaten by predator age u

  // -- 4.13. Kinzey Functional response parameters
  /*
   matrix<Type> H_1(nspp, nspp + 1); H_1 = exp(logH_1.array());
   vector<Type> H_1a(nspp); H_1a = exp(logH_1a);
   vector<Type> H_1b(nspp); H_1b = exp(logH_1b);
   matrix<Type> H_2(nspp, nspp); H_2 = exp(logH_2.array());
   matrix<Type> H_3(nspp, nspp); H_3 = exp(logH_3.array());

   array<Type>  N_pred_yrs(nspp, max_sex, max_age, nyrs); N_pred_yrs.setZero();                // Effective numbers of predators for each age of prey FIXME: should be avgN_at_age?
   array<Type>  N_prey_yrs(nspp, max_sex, max_age, nyrs); N_prey_yrs.setZero();                // Effective numbers of prey for each age of predator
   array<Type>  N_pred_eq(nspp, max_sex, max_age); N_pred_eq.setZero();                        // Effective numbers of predators for each age of prey (styr_pred)
   array<Type>  N_prey_eq(nspp, max_sex, max_age); N_prey_eq.setZero();                        // Effective numbers of prey for each age of predator

   array<Type>  pred_resp(nspp * max_sex, (nspp * max_sex)+1, max_age, max_age, nyrs); pred_resp.setZero();// Predator functional response +1 for other species
   array<Type>  Pred_r(nspp, max_sex, max_age, nyrs); Pred_r.setZero();                          // save Pred_ratio values
   array<Type>  Prey_r(nspp, max_sex, max_age, nyrs); Prey_r.setZero();                          // save Prey_ratio values

   array<Type> ration_hat(nspp, max_sex, max_age, nyrs); ration_hat.setZero();                   // Annual ration by predator age each year
   array<Type> ration_hat_ave(nspp, max_sex, max_age); ration_hat_ave.setZero();                 // Annual ration by predator age averaged over years
   */

  // ------------------------------------------------------------------------- //
  // 5. INITIAL CALCULATIONS                                                   //
  // ------------------------------------------------------------------------- //


  // 5.1. MATURITY AND SEX RATIO
  matrix<Type> mature_females = maturity;
  for( sp = 0; sp < nspp ; sp++) {

    // Sex ratio for SSB derivation
    for( age = 0 ; age < nages(sp); age++ ) {
      if(nsex(sp) == 1){
        mature_females( sp, age ) = maturity( sp, age ) * sex_ratio(sp, age); // Mulitply sex_ratio and maturity for 1 sex models
      }
    }

    // Sex ratio at recruitment (1-sex model gets all R)
    if(nsex(sp)  == 1){
      sex_ratio(sp, 0) = 1.0;
    }
  }


  // 5.4. DATA VARIANCE TERMS
  R_sd = exp(R_ln_sd); // Convert log sd to natural scale
  sel_dev_sd = exp(sel_dev_ln_sd) ;
  index_q_sd = exp(index_q_ln_sd) ;
  index_q_dev_sd = exp(index_q_dev_ln_sd) ;
  Cindex -=1; // Subtract 1 from Cindex to deal with indexing start at 0


  // 5.5. SURVEY CONTROL SWITCHES
  matrix<Type> flt_q(n_flt, nyrs_hind); flt_q.setZero();                        // Vector to save q on natural scale
  vector<int> flt_sel_ind(n_flt); flt_sel_ind.setZero();                        // Vector to store survey index
  vector<int> flt_type(n_flt); flt_type.setZero();                              // Index wether the data are included in the likelihood or not (0 = no, 1 = yes)
  vector<int> flt_sel_type(n_flt); flt_sel_type.setZero();                      // Vector to save survey selectivity type
  vector<int> flt_nselages(n_flt); flt_nselages.setZero();                      // Vector to save number of ages to estimate non-parametric selectivity (1 = age, 2 = length)
  vector<int> flt_varying_sel(n_flt); flt_varying_sel.setZero();                // Vector storing information on wether time-varying selectivity is estimated (0 = no, 1 = random walk with fixed variance, 2 = random effect)
  vector<int> flt_spp(n_flt); flt_spp.setZero();                                // Vector to save survey species
  vector<int> flt_sel_age(n_flt); flt_sel_age.setZero();                        // Vector to save age first selected (selectivity below this age = 0)
  vector<int> flt_sel_maxage(n_flt); flt_sel_maxage.setZero();                  // Vector to save age of max selectivity for normalization (if NA not used)
  vector<int> flt_sel_maxage_upper(n_flt); flt_sel_maxage_upper.setZero();      // Vector to save upper age of max selectivity for normalization (if NA not used)
  vector<int> comp_ll_type(n_flt); comp_ll_type.setZero();                      // Vector to save composition type
  vector<int> flt_units(n_flt); flt_units.setZero();                            // Vector to save survey units (1 = weight, 2 = numbers)
  vector<int> flt_wt_index(n_flt); flt_wt_index.setZero();                      // Vector to save 1st dim of weight to use for weight-at-age
  vector<int> flt_age_transition_index(n_flt); flt_age_transition_index.setZero(); // Vector to save 3rd dim of age_trans_matrix to use for ALK
  vector<int> flt_q_ind(n_flt); flt_q_ind.setZero();                            // Vector storing index of survey q for mapping
  vector<int> est_index_q(n_flt); est_index_q.setZero();                            // Vector to save wether or not analytical q is used
  vector<int> index_varying_q(n_flt); index_varying_q.setZero();                    // Vector storing information on wether time-varying q is estimated (0 = no, 1 = random walk with fixed variance, 2 = random effect)
  vector<int> est_sigma_index(n_flt); est_sigma_index.setZero();                    // Vector to save wether sigma survey is estimated
  vector<int> est_sigma_fsh(n_flt); est_sigma_fsh.setZero();                    // Vector to save wether sigma fishery is estimated


  for(flt_ind = 0; flt_ind < n_flt; flt_ind++){
    flt = fleet_control(flt_ind, 1) - 1;                     // Temporary survey index
    flt_type(flt) = fleet_control(flt_ind, 2);               // Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
    flt_spp(flt) = fleet_control(flt_ind, 3) - 1;            // Species
    flt_sel_ind(flt) = fleet_control(flt_ind, 4);            // Survey selectivity index
    flt_sel_type(flt) = fleet_control(flt_ind, 5);           // Selectivity type
    flt_nselages(flt) = fleet_control(flt_ind, 6);           // Non-parametric selectivity ages
    flt_varying_sel(flt) = fleet_control(flt_ind, 7);        // Time-varying selectivity type.
    flt_sel_age(flt) = fleet_control(flt_ind, 8) - minage(flt_spp(flt));                 // First age selected
    flt_sel_maxage(flt) = fleet_control(flt_ind, 9) - minage(flt_spp(flt));              // Age of max selectivity (used for normalization). If NA, does not normalize
    flt_sel_maxage_upper(flt) = fleet_control(flt_ind, 10) - minage(flt_spp(flt));       // Upper age of max selectivity (used for normalization). If NA, does not normalize
    comp_ll_type(flt) = fleet_control(flt_ind, 11);          // Index indicating wether to do dirichlet multinomial for a fleet's composition data (0 = multinomial; 1 = dirichlet-multinomial)
    flt_units(flt) = fleet_control(flt_ind, 12);             // Survey units
    flt_wt_index(flt) = fleet_control(flt_ind, 13) - 1;      // Dim1 of weight
    flt_age_transition_index(flt) = fleet_control(flt_ind, 14) - 1;     // Dim3 of age transition matrix
    flt_q_ind(flt) = fleet_control(flt_ind, 15) - 1;         // Index of survey q
    est_index_q(flt) = fleet_control(flt_ind, 16);           // Estimate analytical q?
    index_varying_q(flt) = fleet_control(flt_ind, 17);       // Time varying q type (or if AR1, what environmental index to fit to)
    est_sigma_index(flt) = fleet_control(flt_ind, 18);       // Wether to estimate standard deviation of survey time series
    est_sigma_fsh(flt) = fleet_control(flt_ind, 19);         // Wether to estimate standard deviation of fishery time series
  }


  // 5.6. CATCHABILITY
  for(flt = 0; flt < n_flt; flt++){
    for(yr = 0; yr < nyrs_hind; yr++){
      index_q(flt, yr) = exp(index_ln_q(flt) + index_q_dev(flt, yr));                 // Exponentiate

      // Q as a function of environmental index
      if(est_index_q(flt) == 5){
        beta_q_tmp = index_q_beta.row(flt);
        env_q_tmp = env_index.row(yr) ;
        index_q_mult =  env_q_tmp * beta_q_tmp;
        index_q(flt, yr) = exp(index_ln_q(flt) + (index_q_mult).sum());
      }

      // QAR1 deviates fit to environmental index (sensu Rogers et al 2024; 10.1093/icesjms/fsae005)
      if(est_index_q(flt) == 6){
        index_q(flt, yr) = exp(index_ln_q(flt) + index_q_beta(flt, 0) * index_q_dev(flt, yr));
      }
    }
  }


  matrix<int> r_sexes(diet_obs.rows(), 2); r_sexes.setZero();
  matrix<int> k_sexes(diet_obs.rows(), 2); k_sexes.setZero();


  // 5.7. SELECTIVITY
  // -- 5.7.1 EMPIRICAL SELECTIVITY
  sel.setZero();
  for(int sel_ind = 0; sel_ind < emp_sel_obs.rows(); sel_ind++){

    // Fishery
    flt = emp_sel_ctl(sel_ind, 0) - 1;        // Temporary index
    sp = emp_sel_ctl(sel_ind, 1) - 1;         // Temporary index of species
    sex = emp_sel_ctl(sel_ind, 2);            // Temporary index for sex
    yr = emp_sel_ctl(sel_ind, 3) - styr;      // Temporary index for years of data

    // Switch to make sure it doesnt fill out selectivity
    if(flt_sel_type(flt) == 0){

      // 1 sex model
      vector<int> sexes(1); sexes(0) = 0;

      // 2 sex model and emp_sel is for both sex
      if( (nsex(sp) == 2) & (sex == 0) ){
        vector<int> sexes(2); sexes(0) = 0; sexes(1) = 1;
      }
      // 2 sex model and emp_sel is for 1 sex
      if( (nsex(sp) == 2) & (sex > 0) ){
        sexes(0) = sex - 1;
      }

      if(yr < nyrs_hind){
        for(sex = 0; sex < sexes.size(); sex++){
          for(age = 0; age < nages(sp); age++) {
            if(!isNA(emp_sel_obs(sel_ind, age))){
              sel(flt, sexes(sex), age, yr) = emp_sel_obs(sel_ind, age);
            }
          }
        }
      }
    }
  }


  // -- 5.7.2 ESTIMATED SELECTIVITY
  avg_sel.setZero();
  non_par_sel.setZero();
  Type max_sel = 0;
  for(flt = 0; flt < n_flt; flt++) {

    // Temporay indices
    sp = flt_spp(flt);             // Temporary index of species
    sel_type = flt_sel_type(flt);
    nselages = flt_nselages(flt);


    switch(sel_type){

    case 1:  // 6.1.1. Logisitic selectivity
      for(age = 0; age < nages(sp); age++){
        for(yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex++){
            // Random walk and block
            sel(flt, sex, age, yr) = 1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr)) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr )))));

          }
        }
      }
      break;


    case 2:  // 6.1.2. Non-parametric selectivity fit to age ranges (Ianelli 20??). NOTE: This can likely be improved
      for(sex = 0; sex < nsex(sp); sex++){
        for(age = 0; age < nselages; age++) {
          non_par_sel(flt, sex, age) = sel_coff(flt, sex, age);
          avg_sel(flt, sex) +=  exp(sel_coff(flt, sex, age));
        }
        //  Average selectivity up to nselages
        avg_sel(flt, sex) = log(avg_sel(flt, sex) / nselages);

        // Plus group selectivity
        for(age = nselages; age < nages(sp); age++) {
          non_par_sel(flt, sex, age) = non_par_sel(flt, sex, nselages - 1);
        }

        // Average selectivity across all ages
        avgsel_tmp = 0; // Temporary object for average selectivity across all ages
        for(age = 0; age < nages(sp); age++) {
          avgsel_tmp += exp(non_par_sel(flt, sex, age));
        }
        avgsel_tmp = log(avgsel_tmp / nages(sp));

        // Standardize selectivity
        for(age = 0; age < nages(sp); age++) {
          non_par_sel(flt, sex, age) -= avgsel_tmp;
          non_par_sel(flt, sex, age) = exp(non_par_sel(flt, sex, age));
        }
      }

      // Move to rest of years
      for(age = 0; age < nages(sp); age++){
        for(sex = 0; sex < nsex(sp); sex++){
          for(yr = 0; yr < nyrs_hind; yr++) {
            sel(flt, sex, age, yr) = non_par_sel(flt, sex, age);
          }
        }
      }
      break;


    case 3: // 6.1.3. Double logistic (Dorn and Methot 1990)
      for(age = 0; age < nages(sp); age++){
        for(yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex++){

            // Random walk and block
            sel(flt, sex, age, yr) = (1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr)) * ((age + 1) - ( sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr )))))) * // Upper slope
              (1 - (1 / (1 + exp( -exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr )))))));  // Downward slope;
          }
        }
      }
      break;


    case 4: // 6.1.4. Descending logistic (Dorn and Methot 1990)
      for(age = 0; age < nages(sp); age++){
        for(yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex++){

            // Random walk and block
            sel(flt, sex, age, yr) = (1 - (1 / (1 + exp( - exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr )))))));  // Downward slope;
          }
        }
      }
      break;


    case 5: // 6.1.5. Non-parametric selectivity fit to age ranges. Hake version (Taylor et al 2014)
      // -- For each age, sum coefficients from first age selected to age
      // -- Setzero
      for(yr = 0; yr < nyrs_hind; yr++) {
        for(sex = 0; sex < nsex(sp); sex++){
          for(age = 0; age < nages(sp); age++) {
            sel(flt, sex, age, yr) = 0;
          }
        }
      }

      for(yr = 0; yr < nyrs_hind; yr++) {
        for(sex = 0; sex < nsex(sp); sex++){
          for(age = flt_sel_age(flt); age < nselages; age++) {
            for(int age_tmp = flt_sel_age(flt); age_tmp <= age; age_tmp++) {
              sel(flt, sex, age, yr) += sel_coff(flt, sex, age_tmp) + sel_coff_dev(flt, sex, age_tmp, yr);
            }
          }
        }
      }

      // -- For each sex & year, subtract max sel across ages from sel and take exp
      for(yr = 0; yr < nyrs_hind; yr++) {
        for(sex = 0; sex < nsex(sp); sex++){

          // Single-normalization age
          if(flt_sel_maxage(flt) >= 0){
            max_sel = 0.001;
            max_sel = max2(max_sel, sel(flt, sex, flt_sel_maxage(flt), yr)); // Get max sel by sex/year (split or otherwise, divides by 1 for ages > maxselage)
          }

          // Normalize by max
          //FIXME: non-differentiable
          if(flt_sel_maxage(flt) < 0){
            max_sel = 0;
            for(age = flt_sel_age(flt); age < nselages; age++) {
              if(sel(flt, sex, age, yr) > max_sel){
                max_sel = sel(flt, sex, age, yr);
              }
            }
          }

          // Normalize by age rage between max lower and max upper
          if((flt_sel_maxage(flt) >= 0) & (flt_sel_maxage_upper(flt) >= 0)){
            max_sel = 0;
            for(age = flt_sel_maxage(flt); age < flt_sel_maxage_upper(flt); age++) {
              max_sel += sel(flt, sex, age, yr)/(flt_sel_maxage_upper(flt) - flt_sel_maxage(flt) + 1);
            }
          }


          for(age = flt_sel_age(flt); age < nselages; age++) {
            sel(flt, sex, age, yr) = exp(sel(flt, sex, age, yr) - max_sel);
          }

          // Fill in rest of ages
          for(age = nselages-1; age < nages(sp); age++) {
            sel(flt, sex, age, yr) = sel(flt, sex, nselages-1, yr);
          }
        }
      }
      break;
    } // End selectivity switch


    // 6.1.6. Normalize and account for unselected ages
    if(sel_type > 0) {


      // Ages not selected
      for(yr = 0; yr < nyrs_hind; yr++) {
        for(age = 0; age < nages(sp); age++){
          for(sex = 0; sex < nsex(sp); sex++){
            if(age < flt_sel_age(flt)){
              sel(flt, sex, age, yr) = 0;
            }
          }
        }
      }


      // Normalize by selectivity by specific age
      if((flt_sel_maxage(flt) >= 0) & (sel_type < 5)) {
        for(yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex++){

            // Single-normalization age
            if(flt_sel_maxage(flt) >= 0){
              max_sel = 0.001;
              max_sel = max2(max_sel, sel(flt, sex, flt_sel_maxage(flt), yr)); // Get max sel by sex/year (split or otherwise, divides by 1 for ages > maxselage)
            }

            //Normalize by age rage between max lower and max upper
            if((flt_sel_maxage(flt) >= 0) & (flt_sel_maxage_upper(flt) >= 0)){
              max_sel = 0;
              for(age = flt_sel_maxage(flt); age < flt_sel_maxage_upper(flt); age++) {
                max_sel += sel(flt, sex, age, yr)/(flt_sel_maxage_upper(flt) - flt_sel_maxage(flt) + 1);
              }
            }

            // Normalize by max
            for(age = 0; age < nages(sp); age++){
              sel(flt, sex, age, yr) /= max_sel;
            }
          }
        }
      }

      // Find max for each fishery and year across ages, and sexes
      if((sel_type < 5) & (flt_sel_maxage(flt) < 0)) {
        for(yr = 0; yr < nyrs_hind; yr++) {
          max_sel = 0;
          for(age = 0; age < nages(sp); age++){
            for(sex = 0; sex < nsex(sp); sex++){


              // Normalize by max
              if(sel(flt, sex, age, yr) > max_sel){
                max_sel = sel(flt, sex, age, yr);
              }
            }
          }

          // Normalize selectivity
          for(age = 0; age < nages(sp); age++){
            for(sex = 0; sex < nsex(sp); sex++){
              sel(flt, sex, age, yr) /= max_sel;
            }
          }
        }
      }
    }

    // Project outwards assuming final years selectivity
    for(yr = nyrs_hind; yr < nyrs; yr++) {
      for(age = 0; age < nages(sp); age++){
        for(sex = 0; sex < nsex(sp); sex++){
          sel(flt, sex, age, yr) = sel(flt, sex, age, nyrs_hind - 1);
        }
      }
    }
  }

  // ------------------------------------------------------------------------- //
  // 6. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0
  // Start iterations for multi-species convergence
  for(int iter = 0; iter < niter; iter++) {


    // 6.1. FISHING MORTALITY and FSPRs
    F_spp.setZero();
    F_flt_age.setZero();
    F_spp_at_age.setZero();
    Ftarget_age_spp.setZero();
    Flimit_age_spp.setZero();
    for(flt = 0; flt < n_flt; flt++) {

      sp = flt_spp(flt);  // Temporary index of fishery survey

      if(flt_type(flt) == 1){
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            for(yr = 0; yr < nyrs; yr++) {

              // Hindcast
              if( yr < nyrs_hind){
                F_flt_age(flt, sex, age, yr) = sel(flt, sex, age, yr) * exp(ln_F(flt, yr));
              }


              // Forecast
              if( yr >= nyrs_hind){
                // -- Apply HCRs
                switch(HCR){
                case 0: // No fishing
                  proj_F(sp, yr) = 0;
                  break;

                case 1: // CMSY
                  proj_F(sp, yr) = Ftarget(sp);
                  break;

                case 2: // Constant F
                  proj_F(sp, yr) = Ftarget(sp);
                  break;

                case 3: // Constant F that acheives X% of SSB0
                  proj_F(sp, yr) = Ftarget(sp);
                  break;

                case 4: // Constant Fspr
                  proj_F(sp, yr) = Ftarget(sp) * Fmult(sp);
                  break;

                case 5: // NPFMC Tier 3 HCR
                  proj_F(sp, yr) = Ftarget(sp); // Used Fabc of Ftarget_age%
                  break;

                case 6: // PFMC Category 1 HCR
                  proj_F(sp, yr) = Ftarget(sp) = Flimit(sp) + QnormHCR(sp);
                  break;

                case 7: // SESSF Tier 1 HCR
                  proj_F(sp, yr) = Ftarget(sp); // Used Fabc of Ftarget_age%
                  break;
                }

                // Set F to zero if not running forecast
                if(forecast(sp) == 0){
                  proj_F(sp, yr) = 0;
                }


                F_flt_age(flt, sex, age, yr) = sel(flt, sex, age, yr) * proj_F_prop(flt) * proj_F(sp, yr);
              }

              // -- Sum F across fleets
              F_spp_at_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);

              // -- Calculate F target by age and sex for reference points
              Flimit_age_spp(sp, sex, age, yr) += sel(flt, sex, age, yr) * proj_F_prop(flt) * Flimit(sp); // account for time-varying sel
              Ftarget_age_spp(sp, sex, age, yr) += sel(flt, sex, age, yr) * proj_F_prop(flt) * Ftarget(sp); // account for time-varying sel
            }
          }
        }

        // F across fleets or species
        for(yr = 0; yr < nyrs; yr++) {
          // Hindcast
          if( yr < nyrs_hind){
            F_flt(flt, yr) = exp(ln_F(flt, yr));
            F_spp(sp, yr) += exp(ln_F(flt, yr)); // Fully selected fishing mortality
          }

          // Forecast
          if( yr >= nyrs_hind){
            F_flt(flt, yr) = proj_F_prop(flt) * proj_F(sp, yr);
            F_spp(sp, yr) +=  proj_F_prop(flt) * proj_F(sp, yr);
          }
        }
      }
    }


    // 6.2. TOTAL MORTALITY-AT-AGE
    for(sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for(int i = 0; i < M1_beta.dim(2); i++){
          beta_M1_tmp = M1_beta(sp, sex, i);
        }
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs; yr++) {
            // Matrix multiplication from sliced arrays doesn't work
            env_M1_tmp = env_index.row(yr);
            M1_mult = env_M1_tmp * beta_M1_tmp;

            M1_at_age(sp, sex, age, yr) = exp(ln_M1(sp, sex, age) + ln_M1_dev(sp, sex, age, yr) + M1_mult.sum());
            M_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            Z_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + F_spp_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            // S_at_age(sp, sex, age, yr) = exp(-Z_at_age(sp, sex, age, yr));
          }
        }
      }
    }


    // 6.3. SPR BASED REFERENCE POINTS
    //FIXME - make time-varying?
    SPR0.setZero();
    SPRFinit.setZero();
    SPRlimit.setZero();
    SPRtarget.setZero();
    for(sp = 0; sp < nspp; sp++) {

      if(initMode != 3){
        Finit(sp) = 0; // If population starts out at equilibrium set Finit to 0 (R_init and R0 will be the same)
      }

      //FIXME: set to 1
      NbyageSPR(0, sp, 0) = 1.0; // F = 0
      NbyageSPR(1, sp, 0) = 1.0; // F = Flimit
      NbyageSPR(2, sp, 0) = 1.0; // F = Ftarget
      NbyageSPR(3, sp, 0) = 1.0; // F = Finit

      for(age = 1; age < nages(sp)-1; age++) {
        NbyageSPR(0, sp, age) =  NbyageSPR(0, sp, age-1) * exp(-M_at_age(sp, 0, age-1, nyrs_hind - 1));
        NbyageSPR(1, sp, age) =  NbyageSPR(1, sp, age-1) * exp(-(M_at_age(sp, 0, age-1, nyrs_hind - 1) + Flimit_age_spp(sp, 0, age-1, nyrs_hind - 1))); //FIXME: time-vary sel in the forecast
        NbyageSPR(2, sp, age) =  NbyageSPR(2, sp, age-1) * exp(-(M_at_age(sp, 0, age-1, nyrs_hind - 1) + Ftarget_age_spp(sp, 0, age-1, nyrs_hind - 1)));
        NbyageSPR(3, sp, age) =  NbyageSPR(3, sp, age-1) * exp(-(M1_at_age(sp, 0, age-1, nyrs_hind - 1) + Finit(sp)));

      }

      // Plus group
      NbyageSPR(0, sp, nages(sp) - 1) = NbyageSPR(0, sp, nages(sp) - 2) * exp(-M_at_age(sp, 0, nages(sp) - 2, nyrs_hind - 1)) / (1 - exp(-M_at_age(sp, 0, nages(sp) - 1, nyrs_hind - 1)));
      NbyageSPR(1, sp, nages(sp) - 1) = NbyageSPR(1, sp, nages(sp) - 2) * exp(-(M_at_age(sp, 0,  nages(sp) - 2, nyrs_hind - 1) + Flimit_age_spp(sp, 0,  nages(sp) - 2, nyrs_hind - 1))) / (1 - exp(-(M_at_age(sp, 0,  nages(sp) - 1, nyrs_hind - 1) + Flimit_age_spp(sp, 0,  nages(sp) - 1, nyrs_hind - 1))));
      NbyageSPR(2, sp, nages(sp) - 1) = NbyageSPR(2, sp, nages(sp) - 2) * exp(-(M_at_age(sp, 0,  nages(sp) - 2, nyrs_hind - 1) + Ftarget_age_spp(sp, 0,  nages(sp) - 2, nyrs_hind - 1))) / (1 - exp(-(M_at_age(sp, 0,  nages(sp) - 1, nyrs_hind - 1) + Ftarget_age_spp(sp, 0,  nages(sp) - 1, nyrs_hind - 1))));
      NbyageSPR(3, sp, nages(sp) - 1) = NbyageSPR(3, sp, nages(sp) - 2) * exp(-(M1_at_age(sp, 0,  nages(sp) - 2, nyrs_hind - 1) + Finit(sp))) / (1 - exp(-(M1_at_age(sp, 0,  nages(sp) - 1, nyrs_hind - 1) + Finit(sp))));

      // Calulcate SPR
      for(age = 0; age < nages(sp); age++) {
        SPR0(sp) +=  NbyageSPR(0, sp, age) *  weight( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * mature_females(sp, age) * exp(-M_at_age(sp, 0,  age, nyrs_hind - 1) * spawn_month(sp)/12.0);
        SPRlimit(sp) +=  NbyageSPR(1, sp, age) *  weight( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * mature_females(sp, age) * exp(-(M_at_age(sp, 0,  age, nyrs_hind - 1) + Flimit_age_spp(sp, 0,  age, nyrs_hind - 1)) * spawn_month(sp)/12.0);
        SPRtarget(sp) +=  NbyageSPR(2, sp, age) *  weight( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * mature_females(sp, age) * exp(-(M_at_age(sp, 0,  age, nyrs_hind - 1) + Ftarget_age_spp(sp, 0,  age, nyrs_hind - 1)) * spawn_month(sp)/12.0);
        SPRFinit(sp) +=  NbyageSPR(3, sp, age) *  weight( ssb_wt_index(sp), 0, age, 0) * mature_females(sp, age) * exp(-(M1_at_age(sp, 0,  age, nyrs_hind - 1) + Finit(sp)) * spawn_month(sp)/12.0);
      }
    }



    // 6.4. STOCK-RECRUIT PARAMETERS
    // -- For beverton-holt, steepness and R0 are derived from SPR0
    penalty = 0.0;
    zero_N_pen.setZero();
    for( sp = 0; sp < nspp ; sp++) {
      switch(srr_fun){
      case 0: // Random about mean (e.g. Alaska)
        steepness(sp) = 0.99;
        R_init(sp) = R0(sp) = exp(rec_pars(sp, 0));
        break;

      case 1: // Random about mean with environmental linkage
        steepness(sp) = 0.99;
        beta_rec_tmp = beta_rec_pars.row(sp);
        env_rec_tmp = env_index.row(0);
        srr_mult = env_rec_tmp * beta_rec_tmp;
        R_init(sp) = R0(sp) = exp(rec_pars(sp, 0) + srr_mult.sum());
        break;

      case 2: // Beverton-Holt
        steepness(sp) = exp(rec_pars(sp, 1)) * SPR0(sp)/(4.0 + exp(rec_pars(sp, 1)) * SPR0(sp));
        R0(sp) = (exp(rec_pars(sp, 1))-1.0/SPR0(sp)) / exp(rec_pars(sp, 2)); // (Alpha-1/SPR0)/beta
        R_init(sp) = (exp(rec_pars(sp, 1))-1/SPRFinit(sp)) / exp(rec_pars(sp, 2)); // (Alpha-1/SPR0)/beta
        break;

      case 3: // Beverton-Holt with environmental impacts on alpha
        //FIXME make time-varying
        beta_rec_tmp = beta_rec_pars.row(sp);
        env_rec_tmp = env_index.row(0);
        srr_mult = env_rec_tmp * beta_rec_tmp;
        srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
        steepness(sp) = srr_alpha * SPR0(sp)/(4.0 + srr_alpha * SPR0(sp));
        R0(sp) = (srr_alpha-1.0/SPR0(sp)) / exp(rec_pars(sp, 2)); // (Alpha-1/SPR0)/beta
        R_init(sp) = (srr_alpha-1.0/SPRFinit(sp)) / exp(rec_pars(sp, 2)); // (Alpha-1/SPR0)/beta
        break;

      case 4: // Ricker
        steepness(sp) = 0.2 * exp(0.8*log(exp(rec_pars(sp, 1)) * SPR0(sp))); //

        // - R at F0
        ricker_intercept = exp(rec_pars(sp, 1)) * SPR0(sp) - 1.0;
        ricker_intercept =  posfun(ricker_intercept, Type(0.001), penalty) + 1.0;

        R0(sp) = log(ricker_intercept)/(exp(rec_pars(sp, 2)) * SPR0(sp)/1000000.0); // FIXME - make time-varying

        // R at equilibrium F
        ricker_intercept = exp(rec_pars(sp, 1)) * SPRFinit(sp) - 1.0;
        ricker_intercept =  posfun(ricker_intercept, Type(0.001), penalty) + 1.0;

        R_init(sp) = log(ricker_intercept)/(exp(rec_pars(sp, 2)) * SPRFinit(sp)/1000000.0); // FIXME - make time-varying

        zero_N_pen(sp) += penalty;
        break;

      case 5: // Ricker with environmental impacts on alpha
        beta_rec_tmp = beta_rec_pars.row(sp);
        env_rec_tmp = env_index.row(0);
        srr_mult = env_rec_tmp * beta_rec_tmp;
        srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
        steepness(sp) = 0.2 * exp(0.8*log(srr_alpha * SPR0(sp))); //

        ricker_intercept = srr_alpha * SPR0(sp) - 1.0;
        ricker_intercept =  posfun(ricker_intercept, Type(0.001), penalty) + 1.0;
        R0(sp) = log(ricker_intercept)/(exp(rec_pars(sp, 2)) * SPR0(sp)/1000000.0); // FIXME - make time-varying

        ricker_intercept = srr_alpha * SPRFinit(sp) - 1.0;
        ricker_intercept =  posfun(ricker_intercept, Type(0.001), penalty) + 1.0;
        R_init(sp) = log(ricker_intercept)/(exp(rec_pars(sp, 2)) * SPRFinit(sp)/1000000.0); // FIXME - make time-varying
        zero_N_pen(sp) += penalty;
        break;

      default:
        error("Invalid 'srr_fun'");
      }
    }


    // 6.5. INITIAL ABUNDANCE AT AGE, BIOMASS, AND ssb (YEAR 1)
    biomass.setZero();
    ssb.setZero();
    // ssb_at_age.setZero();
    for(sp = 0; sp < nspp; sp++) {
      for(age = 0; age < nages(sp); age++){
        for(sex = 0; sex < nsex(sp); sex ++){


          switch(estDynamics(sp)){
          case 0: // Estimated

            // - Estimate as free parameters
            if(initMode == 0){
              R(sp, 0) = exp(init_dev(sp, 0));
              N_at_age(sp, 0, age, 0) = exp(init_dev(sp, age)) * sex_ratio(sp, 0);
              N_at_age(sp, 1, age, 0) = exp(init_dev(sp, age)) * (1-sex_ratio(sp, 0));
            }

            // - Equilibrium or non-equilibrium estimated as function of R0, Finit, mortality, and init devs
            // Finit is set to 0 when initMode != 2
            if(initMode > 0){
              // -- 6.5.1. Amin (i.e. recruitment)
              if(age == 0){
                R(sp, 0) = R_init(sp) * exp(rec_dev(sp, 0));
                N_at_age(sp, 0, 0, 0) = R(sp, 0) * sex_ratio(sp, 0);
                N_at_age(sp, 1, 0, 0) = R(sp, 0) * (1-sex_ratio(sp, 0));
              }

              // Sum M1 until age - 1
              mort_sum(sp, age) = 0;
              for(int age_tmp = 0; age_tmp < age; age_tmp++){
                mort_sum(sp, age) += M1_at_age(sp, sex, age_tmp, 0) + Finit(sp);
              }

              // -- 6.5.2. Age Amin+1:Amax-1 (initial abundance)
              if((age > 0) & (age < nages(sp) - 1)) {

                if(sex == 0){
                  N_at_age(sp, 0, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1)) * sex_ratio(sp, 0);
                }
                if(sex == 1){
                  N_at_age(sp, 1, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1)) * (1-sex_ratio(sp, 0));
                }
              }

              // -- 6.5.3. Amax
              if(age == (nages(sp) - 1)) {

                if(sex == 0){// NOTE: This solves for the geometric series
                  N_at_age(sp, 0, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1)) / (1 - exp(-M1_at_age(sp, sex, nages(sp) - 1, 0))) * sex_ratio(sp, 0);
                }

                if(sex == 1){
                  N_at_age(sp, 1, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1)) / (1 - exp(-M1_at_age(sp, sex, nages(sp) - 1, 0))) * (1-sex_ratio(sp, 0));
                }
              }
            }
            break;

          case 1: // Fixed numbers-at-age - fixed scalar
            N_at_age(sp, sex, age, 0) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, 0);
            break;

          case 2: // Fixed numbers-at-age age-independent scalar
            N_at_age(sp, sex, age, 0) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, 0);
            break;

          case 3: // Fixed numbers-at-age age-dependent scalar
            N_at_age(sp, sex, age, 0) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, 0);
            break;

          default:
            error("Invalid 'estDynamics'");
          }

          // -- 6.5.3. Estimate total biomass in year 1
          biomass_at_age(sp, sex, age, 0) = N_at_age(sp, sex, age, 0) * weight( pop_wt_index(sp), sex, age, 0);
          biomass(sp, 0) += N_at_age(sp, sex, age, 0) * weight( pop_wt_index(sp), sex, age, 0);
        }

        // -- 6.5.4. Estimated initial female SSB
        // ssb_at_age(sp, age, 0) = N_at_age(sp, 0, age, 0) * pow(S(sp, 0, age, 0), spawn_month(sp)/12) * weight(ssb_wt_index(sp), 0, age, 0 ) * mature_females(sp, age); // 6.6.
        ssb(sp, 0) += N_at_age(sp, 0, age, 0) * exp(-Z_at_age(sp, 0, age, 0) * spawn_month(sp)/12) * weight(ssb_wt_index(sp), 0, age, 0 ) * mature_females(sp, age); // 6.6. ssb_at_age(sp, age, 0);
      }
    }


    // 6.6. HINDCAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and ssb-AT-AGE (kg)
    penalty = 0.0;
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 1; yr < nyrs_hind; yr++) {

        // Switch for srr (for MSEs if using Ianelli SRR method)
        int srr_switch = srr_fun;
        if(yr >= nyrs_srrmean){
          srr_switch = srr_pred_fun;
        }

        // -- 6.6.1. Recruitment
        switch(srr_switch){
        case 0: // Random about mean (e.g. Alaska)
          R(sp, yr) = R0(sp) * exp(rec_dev(sp, yr));
          break;

        case 1: // Random about mean with environmental effects
          beta_rec_tmp = beta_rec_pars.row(sp);
          env_rec_tmp = env_index.row(yr);
          srr_mult = env_rec_tmp * beta_rec_tmp;
          R(sp, yr) = R0(sp) * exp(rec_dev(sp, yr) + srr_mult.sum());
          break;

        case 2: // Beverton-Holt
          R(sp, yr) = exp(rec_pars(sp, 1)) * ssb(sp, yr-minage(sp)) * exp(rec_dev(sp, yr)) / (1.0+exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp)));
          break;

        case 3: // Beverton-Holt with environmental impacts on alpha
          beta_rec_tmp = beta_rec_pars.row(sp);
          env_rec_tmp = env_index.row(yr);
          srr_mult = env_rec_tmp * beta_rec_tmp;
          srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
          R(sp, yr) = srr_alpha * ssb(sp, yr-minage(sp)) * exp(rec_dev(sp, yr)) / (1.0+exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp)));
          break;

        case 4: // Ricker: a * ssb * exp(-beta * ssb). Beta is divided by 1,000,000 for estimation
          R(sp, yr) = exp(rec_pars(sp, 1)) * ssb(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp))/1000000.0) * exp(rec_dev(sp, yr));
          break;

        case 5: // Ricker with environmental impacts on alpha
          beta_rec_tmp = beta_rec_pars.row(sp);
          env_rec_tmp = env_index.row(yr);
          srr_mult = env_rec_tmp * beta_rec_tmp;
          srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
          R(sp, yr) = srr_alpha * ssb(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp))/1000000.0) * exp(rec_dev(sp, yr));
          break;

        default:
          error("Invalid 'srr_fun'");
        }

        N_at_age(sp, 0, 0, yr) = R(sp, yr) * sex_ratio(sp, 0);
        N_at_age(sp, 1, 0, yr) = R(sp, yr) * (1.0-sex_ratio(sp, 0));


        // -- 6.6.2. Ages beyond recruitment
        for(age = 0; age < nages(sp); age++){
          for(sex = 0; sex < nsex(sp); sex ++){

            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- Where Amin < age < Amax
              if(age < (nages(sp) - 1)) {
                N_at_age(sp, sex, age + 1, yr) = N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));// S_at_age(sp, sex, age, yr - 1);
              }

              // -- Plus group where age = Amax.
              if(age == (nages(sp) - 1)) {
                N_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age - 1, yr - 1) * exp(-Z_at_age(sp, sex, age-1, yr-1)) + N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1)); // S_at_age(sp, sex, age, yr - 1);
              }
              break;

            case 1: // Fixed numbers-at-age - fixed scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;

            case 2: // Fixed numbers-at-age age-independent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;
            case 3: // Fixed numbers-at-age age-dependent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, yr);
              break;

            default:
              error("Invalid 'estDynamics'");
            }

            N_at_age(sp, sex, age, yr) = posfun(N_at_age(sp, sex, age, yr), Type(0.001), penalty);
            zero_N_pen(sp) += penalty;

            // -- 6.6.3. Estimate total biomass

            biomass_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * weight( pop_wt_index(sp), sex, age, yr );
            biomass(sp, yr) += biomass_at_age(sp, sex, age, yr);
          }

          // -- 6.6.4. Estimated female ssb
          /*
           ssb_at_age(sp, age, yr) = N_at_age(sp, 0, age, yr) * pow(S_at_age(sp, 0, age, yr), spawn_month(sp)/12.0) * weight(ssb_wt_index(sp), 0, age, yr ) * mature_females(sp, age); // 6.6.
           ssb(sp, yr) += ssb_at_age(sp, age, yr);
           */
          ssb(sp, yr) += N_at_age(sp, 0, age, yr) * exp(-Z_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0) * weight(ssb_wt_index(sp), 0, age, yr ) * mature_females(sp, age); // 6.6.
        }
      }
    }


    // 6.7. DEPLETION BASED BIOMASS REFERENCE POINTS (i.e. SB0 and dynamic SB0)
    // -- calculate mean recruitment
    avg_R.setZero();
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 0; yr < nyrs_srrmean; yr++) {
        avg_R(sp) += R(sp, yr)/Type(nyrs_srrmean); // Update mean rec
      }
    }

    // -- Calc biomass_depletion based BRPs
    NByage0.setZero();
    NByageF.setZero();
    DynamicNByageF.setZero();
    DynamicNByage0.setZero();
    B0.setZero();
    SB0.setZero();
    SBF.setZero();
    DynamicB0.setZero();
    DynamicSB0.setZero();
    DynamicSBF.setZero();
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 0; yr < nyrs; yr++) {

        // - Year 1
        // -- Initial abundance-at-age is the same as the hindcast
        if(yr == 0){
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 0; age < nages(sp); age++){
              NByage0(sp, sex, age , 0) = NByageF(sp, sex, age , 0) = DynamicNByageF(sp, sex, age, 0) = DynamicNByage0(sp, sex, age, 0) = N_at_age(sp,sex,age,0);
            }
          }
        }

        // - Recruitment Year > 1
        if(yr > 0){
          // - Option 1a: Use mean rec
          if((proj_mean_rec == 1) & (srr_pred_fun < 2)){
            // Equilibrium RPs used mean rec across all years
            NByage0(sp, 0, 0, yr) = avg_R(sp);
            NByageF(sp, 0, 0, yr) = avg_R(sp);

            // Dynamic RPs use mean rec for forecast and observed rec for hindcast
            if(yr < nyrs_hind){
              DynamicNByageF(sp, 0, 0, yr) = DynamicNByage0(sp, 0, 0, yr) = R(sp, yr); // Hindcast use observed R
            }
            if(yr >= nyrs_hind){
              DynamicNByageF(sp, 0, 0, yr) = DynamicNByage0(sp, 0, 0, yr) = exp(log(avg_R(sp)) + rec_dev(sp, yr)); // Projections use mean R given bias in R0
            }
          }

          // - Option 1b: Use mean rec and env
          if((proj_mean_rec == 1) & (srr_pred_fun > 1)){
            beta_rec_tmp = beta_rec_pars.row(sp);
            env_rec_tmp = env_index.row(yr);
            srr_mult = env_rec_tmp * beta_rec_tmp;
            NByage0(sp, 0, 0, yr) = avg_R(sp) * exp(srr_mult.sum());

            // Dynamic RPs use mean rec for forecast and observed rec for hindcast
            if(yr < nyrs_hind){
              DynamicNByageF(sp, 0, 0, yr) = DynamicNByage0(sp, 0, 0, yr) = R(sp, yr); // Hindcast use observed R
            }
            if(yr >= nyrs_hind){
              DynamicNByageF(sp, 0, 0, yr) = DynamicNByage0(sp, 0, 0, yr) = exp(log(avg_R(sp)) + rec_dev(sp, yr) + srr_mult.sum()); // Projections use mean R given bias in R0
            }
          }


          // - Option 2: Use SRR
          if(proj_mean_rec == 0){

            // Recruitment
            switch(srr_pred_fun){
            case 0: // Random about mean (e.g. Alaska)
              NByage0(sp, 0, 0, yr) = R0(sp);
              NByageF(sp, 0, 0, yr) = R0(sp);

              DynamicNByage0(sp, 0, 0, yr) = R0(sp) * exp(rec_dev(sp, yr));
              DynamicNByageF(sp, 0, 0, yr) = R0(sp) * exp(rec_dev(sp, yr));
              break;

            case 1: // Random about mean with environmental covariates (e.g. Alaska)
              NByage0(sp, 0, 0, yr) = R0(sp);
              NByageF(sp, 0, 0, yr) = R0(sp);

              beta_rec_tmp = beta_rec_pars.row(sp);
              env_rec_tmp = env_index.row(yr);
              srr_mult = env_rec_tmp * beta_rec_tmp;

              DynamicNByage0(sp, 0, 0, yr) = R0(sp) * exp(rec_dev(sp, yr) + srr_mult.sum());
              DynamicNByageF(sp, 0, 0, yr) = R0(sp) * exp(rec_dev(sp, yr) + srr_mult.sum());
              break;

            case 2: // Beverton-Holt
              // FIXME: error if minage > 1
              NByage0(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * SB0(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * SB0(sp, yr-minage(sp)));
              NByageF(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * SBF(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * SBF(sp, yr-minage(sp)));

              DynamicNByage0(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * DynamicSB0(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * DynamicSB0(sp, yr-minage(sp))) * exp(rec_dev(sp, yr));
              DynamicNByageF(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * DynamicSBF(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * DynamicSBF(sp, yr-minage(sp))) * exp(rec_dev(sp, yr));
              break;

            case 3: // Beverton-Holt with environmental impacts on alpha
              // FIXME: error if minage > 1
              beta_rec_tmp = beta_rec_pars.row(sp);
              env_rec_tmp = env_index.row(yr);
              srr_mult = env_rec_tmp * beta_rec_tmp;
              srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
              NByage0(sp, 0, 0, yr) = srr_alpha * SB0(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * SB0(sp, yr-minage(sp)));
              NByageF(sp, 0, 0, yr) = srr_alpha * SBF(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * SBF(sp, yr-minage(sp)));

              DynamicNByage0(sp, 0, 0, yr) = srr_alpha * DynamicSB0(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * DynamicSB0(sp, yr-minage(sp))) * exp(rec_dev(sp, yr));
              DynamicNByageF(sp, 0, 0, yr) = srr_alpha * DynamicSBF(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * DynamicSBF(sp, yr-minage(sp))) * exp(rec_dev(sp, yr));
              break;

            case 4: // Ricker
              // FIXME: error if minage > 1
              NByage0(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * SB0(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * SB0(sp, yr-minage(sp))/1000000.0);
              NByageF(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * SBF(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * SBF(sp, yr-minage(sp))/1000000.0);

              DynamicNByage0(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * DynamicSB0(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * DynamicSB0(sp, yr-minage(sp))/1000000.0) * exp(rec_dev(sp, yr));
              DynamicNByageF(sp, 0, 0, yr) = exp(rec_pars(sp, 1)) * DynamicSBF(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * DynamicSBF(sp, yr-minage(sp))/1000000.0) * exp(rec_dev(sp, yr));
              break;

            case 5: // Ricker with environmental impacts on alpha
              // FIXME: error if minage > 1
              beta_rec_tmp = beta_rec_pars.row(sp);
              env_rec_tmp = env_index.row(yr);
              srr_mult = env_rec_tmp * beta_rec_tmp;
              srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
              NByage0(sp, 0, 0, yr) = srr_alpha * SB0(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * SB0(sp, yr-minage(sp))/1000000.0);
              NByageF(sp, 0, 0, yr) = srr_alpha * SBF(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * SBF(sp, yr-minage(sp))/1000000.0);

              DynamicNByage0(sp, 0, 0, yr) = srr_alpha * DynamicSB0(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * DynamicSB0(sp, yr-minage(sp))/1000000.0) * exp(rec_dev(sp, yr));
              DynamicNByageF(sp, 0, 0, yr) = srr_alpha * DynamicSBF(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * DynamicSBF(sp, yr-minage(sp))/1000000.0) * exp(rec_dev(sp, yr));
              break;

            default:
              error("Invalid 'srr_fun'");
            }
          } // End recruitment switch

          // Account for sex ratio
          NByage0(sp, 0, 0, yr) = NByage0(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          NByage0(sp, 1, 0, yr) = NByage0(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0)); // Males (divide by sex_r because we multiplied in the line before)

          NByageF(sp, 0, 0, yr) = NByageF(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          NByageF(sp, 1, 0, yr) = NByageF(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0)); // Males (divide by sex_r because we multiplied in the line before)

          DynamicNByage0(sp, 0, 0, yr) = DynamicNByage0(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          DynamicNByage0(sp, 1, 0, yr) = DynamicNByage0(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0)); // Males (divide by sex_r because we multiplied in the line before)

          DynamicNByageF(sp, 0, 0, yr) = DynamicNByageF(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          DynamicNByageF(sp, 1, 0, yr) = DynamicNByageF(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0)); // Males (divide by sex_r because we multiplied in the line before)


          // N-at-age year > 1
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 1; age < nages(sp)-1; age++) {
              NByage0(sp, sex, age, yr) =  NByage0(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr-1)); // F = 0

              NByageF(sp, sex, age, yr) =  NByageF(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr-1) - Ftarget_age_spp(sp, sex, age-1, yr-1)); // F = target

              DynamicNByage0(sp, sex, age, yr) =  DynamicNByage0(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr - 1)); // F = 0

              DynamicNByageF(sp, sex, age, yr) =  DynamicNByageF(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr - 1) - Ftarget_age_spp(sp, sex, age-1, yr-1)); // F = Ftarget
            }

            // Plus group  -- No fishing
            NByage0(sp, sex, nages(sp)-1, yr) = NByage0(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1))  + NByage0(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1));

            NByageF(sp, sex, nages(sp)-1, yr) = NByageF(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1) - Ftarget_age_spp(sp, sex, nages(sp)-2, yr-1))  + NByageF(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1) - Ftarget_age_spp(sp, sex, nages(sp)-1, yr-1));

            DynamicNByage0(sp, sex, nages(sp)-1, yr) = DynamicNByage0(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1))  + DynamicNByage0(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1));

            DynamicNByageF(sp, sex, nages(sp)-1, yr) = DynamicNByageF(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1) - Ftarget_age_spp(sp, sex, nages(sp)-2, yr-1))  + DynamicNByageF(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1) - Ftarget_age_spp(sp, sex, nages(sp)-1, yr-1));

          }
        }


        // Calculate Dynamic SB0 and SB at F target

        if(yr < nyrs_hind){
          yr_ind = yr;
        }
        if(yr >= nyrs_hind){
          yr_ind = nyrs_hind - 1;
        }

        for(age = 0; age < nages(sp); age++) {
          SB0(sp, yr) +=  NByage0(sp, 0, age, yr) *  weight( ssb_wt_index(sp), 0, age, nyrs_hind - 1 ) * mature_females(sp, age) * exp(-M_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0);
          SBF(sp, yr) +=  NByageF(sp, 0, age, yr) *  weight( ssb_wt_index(sp), 0, age, nyrs_hind - 1 ) * mature_females(sp, age) * exp(-(M_at_age(sp, 0, age, yr) + Ftarget_age_spp(sp, 0, age, yr)) * spawn_month(sp)/12.0);
          DynamicSB0(sp, yr) +=  DynamicNByage0(sp, 0, age, yr) *  weight( ssb_wt_index(sp), 0, age, yr_ind ) * mature_females(sp, age) * exp(-M_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0);
          DynamicSBF(sp, yr) +=  DynamicNByageF(sp, 0, age, yr) *  weight( ssb_wt_index(sp), 0, age, yr_ind ) * mature_females(sp, age) * exp(-(M_at_age(sp, 0, age, yr) + Ftarget_age_spp(sp, 0, age, yr)) * spawn_month(sp)/12.0);

          for(sex = 0; sex < nsex(sp); sex ++){
            B0(sp, yr) +=  NByage0(sp, sex, age, yr) *  weight( ssb_wt_index(sp), sex, age, nyrs_hind - 1 );
            DynamicB0(sp, yr) +=  DynamicNByage0(sp, sex,age, yr) *  weight( ssb_wt_index(sp), sex, age, yr_ind );
          }
        }


        // Input SB0 (if running in multi-species mode)
        if(msmMode > 0){
          SB0(sp, yr) = MSSB0(sp);
          B0(sp, yr) = MSB0(sp);
        }
      }
    }

    // 6.8-6.9. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and ssb-AT-AGE (kg)
    // Includes Harvest Control Rules
    for(sp = 0; sp < nspp; sp++) {
      for(yr = nyrs_hind; yr < nyrs; yr++){

        // 6.8. HARVEST CONTROL RULES FOR PROJECTION (i.e. SB0 and dynamic SB0)
        // -- Equilibrium Harvest Control Rules
        if(DynamicHCR == 0){
          switch(HCR){
          case 0: // No fishing
            proj_F(sp, yr) = 0.0;
            break;

          case 1: // CMSY
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 2: // Constant F
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 3: // Constant F to acheive X% of SB0
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 4: // Constant Fspr with multiplier
            proj_F(sp, yr) = proj_F(sp, yr) * Fmult(sp);
            break;

          case 5: // NPFMC Tier 3 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(ssb(sp, yr-1) < SBF(sp, nyrs-1)){
              proj_F(sp, yr) = Ftarget(sp) * (((ssb(sp, yr-1)/SBF(sp, nyrs-1))-Alpha(sp))/(1-Alpha(sp))); // Used Fabc of FtargetSPR%
            }
            if((ssb(sp, nyrs_hind-1) < SB0(sp, nyrs-1) * Plimit(sp)) | (ssb(sp, yr-1) / SBF(sp, nyrs-1) < Alpha(sp))){ // If overfished
              proj_F(sp, yr) = 0.0;
            }
            break;

          case 6: // PFMC Category 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(ssb(sp, yr-1) < SB0(sp, nyrs-1) * Ptarget(sp)){
              proj_F(sp, yr) = (Flimit(sp) + QnormHCR(sp)) * (SB0(sp, nyrs-1) * Ptarget(sp) * (ssb(sp, yr-1) - SB0(sp, nyrs-1) * Plimit(sp))) / (ssb(sp, yr-1) * (SB0(sp, nyrs-1) * (Ptarget(sp) - Plimit(sp))));
            }
            if(ssb(sp, yr-1) < SB0(sp, nyrs-1) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) = 0.0;
            }
            break;

          case 7: // SESSF Tier 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(ssb(sp, yr-1) < SB0(sp, nyrs-1) * Ptarget(sp)){
              proj_F(sp, yr) = Ftarget(sp) * ((ssb(sp, yr-1)/(SB0(sp, nyrs-1) * Plimit(sp)))-1); // Used Fabc of FtargetSPR%
            }
            if(ssb(sp, yr-1) < SB0(sp, nyrs-1) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) = 0.0;
            }
            break;
          }
        }

        // Dynamic Harvest Control Rules
        if(DynamicHCR == 1){
          switch(HCR){
          case 0: // No fishing
            proj_F(sp, yr) = 0.0;
            break;

          case 1: // CMSY
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 2: // Constant F
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 3: // Constant F to acheive X% of SB0
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 4: // Constant Fspr with multiplier
            proj_F(sp, yr) = proj_F(sp, yr) * Fmult(sp);
            break;

          case 5: // NPFMC Tier 3 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(ssb(sp, yr-1) < DynamicSBF(sp, yr-1)){
              proj_F(sp, yr) = Ftarget(sp) * (((ssb(sp, yr-1)/(DynamicSBF(sp, yr-1)))-Alpha(sp))/(1-Alpha(sp))); // Used Fabc of FtargetSPR%
            }
            if((ssb(sp, yr-1) < DynamicSB0(sp, yr-1) * Plimit(sp)) | (ssb(sp, yr-1) / (DynamicSBF(sp, yr-1)) < Alpha(sp))){ // If overfished
              proj_F(sp, yr) = 0.0;
            }
            break;

          case 6: // PFMC Category 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(ssb(sp, yr-1) < DynamicSB0(sp, yr-1) * Ptarget(sp)){
              proj_F(sp, yr) = (Flimit(sp) + QnormHCR(sp)) * (DynamicSB0(sp, yr-1) * Ptarget(sp) * (ssb(sp, yr-1) - DynamicSB0(sp, yr-1) * Plimit(sp))) / (ssb(sp, yr-1) * (DynamicSB0(sp, yr-1) * (Ptarget(sp) - Plimit(sp))));
            }
            if(ssb(sp, yr-1) < DynamicSB0(sp, yr-1) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) = 0.0;
            }
            break;

          case 7: // SESSF Tier 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(ssb(sp, yr-1) < DynamicSB0(sp, yr-1) * Ptarget(sp)){
              proj_F(sp, yr) = Ftarget(sp) * ((ssb(sp, yr-1)/(DynamicSB0(sp, yr-1) * Plimit(sp)))-1); // Used Fabc of FtargetSPR%
            }
            if(ssb(sp, yr-1) < DynamicSB0(sp, yr-1) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) =  0.0;
            }
            break;
          }
        }

        // Set F to 0 if not forecast
        if(forecast(sp) == 0){
          proj_F(sp, yr) =  0.0;
        }

        // Adjust F*selex
        // -- 6.8.3. Update F for the projection (account for selectivity and fleets)
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            F_spp_at_age(sp, sex, age, yr) = 0.0;
          }
        }

        // -- Multiply F from HCR by selectivity and fleet proportion
        F_spp(sp, yr) = proj_F(sp, yr);
        for(flt = 0; flt < n_flt; flt++) {
          if(sp == flt_spp(flt)){
            F_flt(sp, yr) = proj_F_prop(flt) * proj_F(sp, yr);
            for(age = 0; age < nages(sp); age++) {
              for(sex = 0; sex < nsex(sp); sex ++){
                F_flt_age(flt, sex, age, yr) = sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * proj_F(sp, yr); // FIXME using last year of selectivity
                if(flt_type(flt) == 1){
                  F_spp_at_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);
                }
              }
            }
          }
        }


        // -- 6.8.4. Update mortality for forcast
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            M_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            Z_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + F_spp_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            //S_at_age(sp, sex, age, yr) = exp(-Z_at_age(sp, sex, age, yr));
          }
        }


        // 6.9. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and ssb-AT-AGE (kg)
        // -- 6.9.1. Forecasted recruitment
        // - Option 1: Use mean rec
        if((proj_mean_rec == 1) & (srr_pred_fun > 1)){
          R(sp, yr) = exp(log(avg_R(sp)) + rec_dev(sp, yr)); //  // Projections use mean R given bias in R0
        }

        // - Mean rec and environment
        if((proj_mean_rec == 1) & (srr_pred_fun < 2)){
          beta_rec_tmp = beta_rec_pars.row(sp);
          env_rec_tmp = env_index.row(yr);
          srr_mult = env_rec_tmp * beta_rec_tmp;
          R(sp, yr) = exp(log(avg_R(sp)) + rec_dev(sp, yr)) * exp(srr_mult.sum());
        }

        // - Option 2: Use SRR and rec devs
        if(proj_mean_rec == 0){
          switch(srr_pred_fun){
          case 0: // Random about mean (e.g. Alaska)
            R(sp, yr) = R0(sp) * exp(rec_dev(sp, yr));
            break;

          case 1: // Random about mean with environmental effects
            beta_rec_tmp = beta_rec_pars.row(sp);
            env_rec_tmp = env_index.row(yr);
            srr_mult = env_rec_tmp * beta_rec_tmp;
            R(sp, yr) = R0(sp) * exp(rec_dev(sp, yr) + srr_mult.sum());
            break;

          case 2: // Beverton-Holt
            R(sp, yr) = exp(rec_pars(sp, 1)) * ssb(sp, yr-minage(sp)) * exp(rec_dev(sp, yr)) / (1+exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp)));
            break;

          case 3: // Beverton-Holt with environmental impacts on alpha
            beta_rec_tmp = beta_rec_pars.row(sp);
            env_rec_tmp = env_index.row(yr);
            srr_mult = env_rec_tmp * beta_rec_tmp;
            srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
            R(sp, yr) = srr_alpha * ssb(sp, yr-minage(sp)) * exp(rec_dev(sp, yr)) / (1+exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp)));
            break;

          case 4: // Ricker
            R(sp, yr) = exp(rec_pars(sp, 1)) * ssb(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp))/1000000.0)* exp(rec_dev(sp, yr));
            break;

          case 5: // Ricker with environmental impacts on alpha
            beta_rec_tmp = beta_rec_pars.row(sp);
            env_rec_tmp = env_index.row(yr);
            srr_mult = env_rec_tmp * beta_rec_tmp;
            srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
            R(sp, yr) = srr_alpha * ssb(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp))/1000000.0)* exp(rec_dev(sp, yr));
            break;

          default:
            error("Invalid 'srr_pred_fun'");
          }
        }

        N_at_age(sp, 0, 0 , yr) = R(sp, yr) * sex_ratio(sp, 0);
        N_at_age(sp, 1, 0 , yr) = R(sp, yr) * (1 - sex_ratio(sp, 0));


        // -- 6.9.2. Ages > recruitment
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- 6.9.2.  Where minage + 1 <= age < Ai
              if(age < (nages(sp) - 1)) {
                N_at_age(sp, sex, age + 1, yr) = N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));// S_at_age(sp, sex, age, yr - 1);
              }

              // -- 6.9.3. Plus group where age > Ai.
              if(age == (nages(sp) - 1)) {
                N_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age - 1, yr - 1) * exp(-Z_at_age(sp, sex, age-1, yr-1)) + N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));// S_at_age(sp, sex, age, yr - 1);
              }
              break;

            case 1: // Fixed numbers-at-age - fixed scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;

            case 2: // Fixed numbers-at-age age-independent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;

            case 3: // Fixed numbers-at-age age-dependent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, yr);
              break;

            default: // Wrong estDynamics
              error("Invalid 'estDynamics'");
            }

            //  Constraint to reduce population collapse
            // N_at_age(sp, sex, age, yr) = posfun(N_at_age(sp, sex, age, yr), Type(0.001), penalty);
            // zero_N_pen(sp) += penalty;

            // -- 6.9.4. FORECAST BIOMASS
            biomass_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * weight( pop_wt_index(sp), sex, age, nyrs_hind-1 ); // 6.5.
            biomass(sp, yr) += biomass_at_age(sp, sex, age, yr);
          } // End sex loop

          // -- 6.9.5. FORECAST SSB (SUM ACROSS AGES)
          /*
           ssb_at_age(sp, age, yr) = N_at_age(sp, 0, age, yr) * pow(S_at_age(sp, 0, age, yr), spawn_month(sp)/12.0) * weight(ssb_wt_index(sp), 0, age, nyrs_hind-1 ) * mature_females(sp, age); // 6.6.
           ssb(sp, yr) += ssb_at_age(sp, age, yr);
           */
          ssb(sp, yr) += N_at_age(sp, 0, age, yr) * exp(-Z_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0) * weight(ssb_wt_index(sp), 0, age, nyrs_hind-1 ) * mature_females(sp, age); // 6.6.
        }
      }
    }


    // 6.10. EXPECTED RECRUITMENT
    for(sp = 0; sp < nspp; sp++) {

      // Year 1 (arent fit in likelihood)
      switch(srr_pred_fun){
      case 0: // Random about mean (e.g. Alaska)
        R_hat(sp, 0)  = R0(sp);
        break;

      case 1: // Random about mean with environmental effects
        beta_rec_tmp = beta_rec_pars.row(sp);
        env_rec_tmp = env_index.row(0);
        srr_mult = env_rec_tmp * beta_rec_tmp;
        R_hat(sp, 0) = R0(sp) * exp(srr_mult.sum());
        break;

      case 2: // Beverton-Holt
        R_hat(sp, 0) = (exp(rec_pars(sp, 1))-1/SPRFinit(sp)) / exp(rec_pars(sp, 2)); // (Alpha-1/SPR0)/beta
        break;

      case 3: // Beverton-Holt with environmental impacts on alpha
        beta_rec_tmp = beta_rec_pars.row(sp);
        env_rec_tmp = env_index.row(0);
        srr_mult = env_rec_tmp * beta_rec_tmp;
        srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
        R_hat(sp, 0) = (srr_alpha-1/SPRFinit(sp)) / exp(rec_pars(sp, 2)); // (Alpha-1/SPR0)/beta
        break;

      case 4: // Ricker
        R_hat(sp, 0) = log(exp(rec_pars(sp, 1)) * SPRFinit(sp))/(exp(rec_pars(sp, 2)) * SPRFinit(sp)/1000000.0); // FIXME - make time-varying
        break;

      case 5: // Ricker with environmental impacts on alpha
        beta_rec_tmp = beta_rec_pars.row(sp);
        env_rec_tmp = env_index.row(0);
        srr_mult = env_rec_tmp * beta_rec_tmp;
        srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
        R_hat(sp, 0) = log(srr_alpha * SPRFinit(sp))/(exp(rec_pars(sp, 2)) * SPRFinit(sp)/1000000.0); // FIXME - make time-varying
        break;
      default:
        error("Invalid 'srr_pred_fun'");
      }

      // Year 1+
      for(yr = 1; yr < nyrs; yr++){
        switch(srr_pred_fun){
        case 0: // Random about mean (e.g. Alaska)
          R_hat(sp, yr)  = R0(sp);
          break;

        case 1: // Random about mean with environmental effects
          beta_rec_tmp = beta_rec_pars.row(sp);
          env_rec_tmp = env_index.row(yr);
          srr_mult = env_rec_tmp * beta_rec_tmp;
          R_hat(sp, 0) = R0(sp) * exp(srr_mult.sum());
          break;

        case 2: // Beverton-Holt
          R_hat(sp, yr) = exp(rec_pars(sp, 1)) * ssb(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp)));
          break;

        case 3: // Beverton-Holt with environmental impacts on alpha
          beta_rec_tmp = beta_rec_pars.row(sp);
          env_rec_tmp = env_index.row(yr);
          srr_mult = env_rec_tmp * beta_rec_tmp;
          srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
          R_hat(sp, yr) = srr_alpha * ssb(sp, yr-minage(sp)) / (1+exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp)));
          break;

        case 4: // Ricker
          R_hat(sp, yr) = exp(rec_pars(sp, 1)) * ssb(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp))/1000000.0); // Divide by 1e6 for numerical stability
          break;

        case 5: // Ricker with environmental impacts on alpha
          beta_rec_tmp = beta_rec_pars.row(sp);
          env_rec_tmp = env_index.row(yr);
          srr_mult = env_rec_tmp * beta_rec_tmp;
          srr_alpha = exp(rec_pars(sp, 1) + srr_mult.sum());
          R_hat(sp, yr) = srr_alpha * ssb(sp, yr-minage(sp)) * exp(-exp(rec_pars(sp, 2)) * ssb(sp, yr-minage(sp))/1000000.0);
          break;

        default:
          error("Invalid 'srr_pred_fun'");
        }
      }
    }



    // 6.11. ESTIMATE AVERAGE NUMBERS AT AGE
    for(sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs; yr++) {
            // switch(avgnMode){
            // case 0: // MSVPA approach
            avgN_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) / Z_at_age(sp, sex, age, yr);
            /*
             break;
          case 1: // Kinzey and Punt (2009) approximation
             avgN_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * exp(- Z_at_age(sp, sex, age, yr) / 2);
             break;
          case 2: // Van Kirk et al (2010) approximation
             avgN_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr);
             break;
          default:
             error("Invalid 'avgnMode'");
             }
             */
          }
        }
      }
    }


    // ------------------------------------------------------------------------- //
    // 7. RATION EQUATIONS                                                       //
    // ------------------------------------------------------------------------- //
    // NOTE -- LOOPING INDICES -- sp = species, age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length

    // 7.2. Calculate temperature function of consumption
    Type Yc = 0.0;
    Type Zc = 0.0;
    Type Vc = 0.0;
    Type Xc = 0.0;
    Type G2 = 0.0;
    Type L2 = 0.0;
    Type G1 = 0.0;
    Type L1 = 0.0;
    Type Ka = 0.0;
    Type Kb = 0.0;
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 0; yr < nyrs; yr++) {

        switch(Ceq(sp)){
        case 1:// Exponential function from Stewart et al. 1983
          fT(sp, yr) = exp(Qc(sp) * env_index(yr, Cindex(sp)));
          break;

        case 2:// Temperature dependence for warm-water-species from Kitchell et al. 1977
          Yc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp) + 2.0);
          Zc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp));
          Vc = (Tcm(sp) - env_index(yr, Cindex(sp))) / (Tcm(sp) - Tco(sp));
          Xc = pow(Zc, 2) * pow((1.0 + pow((1.0 + 40.0 / Yc), 0.5)), 2) / 400.0;
          fT(sp, yr) = pow(Vc, Xc) * exp(Xc * (1.0 - Vc));
          break;

        case 3:// Temperature dependence for cool and cold-water species from Thornton and Lessem 1979
          G2 = (1.0 / (Tcl(sp) - Tcm(sp))) * log((0.98 * (1.0 - CK4(sp))) / (CK4(sp) * 0.02));
          L2 = exp(G2 * (Tcl( sp ) -  env_index(yr, Cindex(sp))));
          Kb = (CK4(sp) * L2) / (1.0 + CK4(sp) * (L2 - 1.0));
          G1 = (1.0 / (Tco(sp) - Qc(sp))) * log((0.98 * (1.0 - CK1(sp))) / (CK1(sp) * 0.02));
          L1 = exp(G1 * (env_index(yr, Cindex(sp)) - Qc(sp)));
          Ka = (CK1(sp) * L1) / (1.0 + CK1(sp) * (L1 - 1.0));
          fT(sp, yr) = Ka * Kb;
          break;

        case 4:
          fT(sp, yr) = 1.0;
          break;
        }
      }
    }


    // 7.3. Calculate historic ration
    for(sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex++) {
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs; yr++) {
            // p = proportion of maximum consumption
            // f(T) = temperature dependence function
            // CA = intercept of allometric mass function
            // CB = slope of allometric mass function
            // fday = number of forageing days per year

            // Hindcast
            if(yr < nyrs_hind){
              consumption_at_age(sp, sex, age, yr) = CA(sp) * pow(weight( pop_wt_index(sp), sex, age, yr ) * Type(1000.0), 1 + CB( sp )) //  C_max = CA * W ^ 1+CB; where C_max is grams consumed per grams of predator per day
              * fT(sp, yr) * fday( sp );                           //  C_max * f(T) * weight * fday g/pred.yr
              consumption_at_age(sp, sex, age, yr) = consumption_at_age(sp, sex, age, yr) * Pvalue(sp) * Pyrs(sp, sex, age, yr); //
            }

            // Projection (wt indexing)
            if(yr >= nyrs_hind){
              consumption_at_age(sp, sex, age, yr) = CA(sp) * pow(weight( pop_wt_index(sp), sex, age, (nyrs_hind - 1) ) * Type(1000.0), 1 + CB( sp ))  //  C_max = CA * W ^ 1+CB; where C_max is grams consumed per grams of predator per day
              * fT(sp, yr) * fday( sp );                            //  C_max * f(T) * weight * fday g/pred.yr
              consumption_at_age(sp, sex, age, yr) = consumption_at_age(sp, sex, age, yr) * Pvalue(sp) * Pyrs(sp, sex, age, (nyrs_hind-1)); //
            }

            ration(sp, sex, age, yr) = consumption_at_age(sp, sex, age, yr) / 1000.0;      // Annual ration kg/yr
          }
        }
      }
    }


    // 7.4. Reorganize diet_obs content
    for(int stom_ind = 0; stom_ind < diet_obs.rows(); stom_ind++){
      rsp = diet_ctl(stom_ind, 0) - 1; // Index of pred
      ksp = diet_ctl(stom_ind, 1) - 1; // Index of prey
      r_sex = diet_ctl(stom_ind, 2); // Index of pred sex
      k_sex = diet_ctl(stom_ind, 3); // Index of prey sex
      r_age = diet_ctl(stom_ind, 4) - minage(rsp); // Index of pred age
      k_age = diet_ctl(stom_ind, 5) - minage(ksp); // Index of prey age
      flt_yr = diet_ctl(stom_ind, 6); // Index of year

      // Predator
      // 1 sex model
      r_sexes(stom_ind, 0) = 0; r_sexes(stom_ind, 1) = 0;
      k_sexes(stom_ind, 0) = 0; k_sexes(stom_ind, 1) = 0;

      // 2 sex model
      // This is to account for situations where nsex = 2, but r_sex or k_sex = 0
      if(nsex(rsp) == 2){
        // But k_sex = 0 (indicating diet data is for both sexes)
        r_sexes(stom_ind, 0) = 0; r_sexes(stom_ind, 1) = 1;

        if(r_sex > 0){
          r_sexes(stom_ind, 0) = r_sex - 1;  r_sexes(stom_ind, 1) = r_sex - 1;
        }
      }

      if(nsex(ksp) == 2){
        // But k_sex = 0 (indicating diet data is for both sexes)
        k_sexes(stom_ind, 0) = 0; k_sexes(stom_ind, 1) = 1;

        // Sex-specific diet data
        if(k_sex > 0){
          k_sexes(stom_ind, 0) = k_sex - 1;  k_sexes(stom_ind, 1) = k_sex - 1;
        }
      }

      // Diet proportion of prey-at-at in predator-at-age
      // if K_age < 0, data are diet proportion of prey-spp in predator-at-age (summed across prey ages)
      // and are not used for MSVPA based suitability
      if(k_age >= 0){
        for(int j = 0; j < 2; j ++){
          for(int k = 0; k < 2; k ++){

            // Annual data
            if(flt_yr > 0){
              yr = flt_yr - styr;
              if(yr < nyrs_hind){
                diet_prop(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr) = diet_obs(stom_ind, 1);
              }
            }

            // Average of years
            if(flt_yr == 0){
              for(yr = 0; yr < nyrs; yr++) {
                diet_prop(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr) = diet_obs(stom_ind, 1);
              }
            }
          }
        }
      }
    }


    // 7.5. Calculate other food stomach content
    other_food_diet_prop.setZero();
    for(yr = 0; yr < nyrs; yr++) {                           // Year loop
      for(rsp = 0; rsp < nspp; rsp++) {                      // Predator species loop
        for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
          for(r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
            other_food_diet_prop(rsp, r_sex, r_age, yr) = Type( 1 );             // Initialize other suitability
            for(ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for(k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                  other_food_diet_prop(rsp, r_sex, r_age, yr) -= diet_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr);
                }
              }
            }
            if(other_food(rsp) > 0){
              other_food_diet_prop(rsp, r_sex, r_age, yr) /= other_food(rsp); // Penalize this
            }
            if(other_food(rsp) == 0){
              other_food_diet_prop(rsp, r_sex, r_age, yr) = 0.0;
            }
          }
        }
      }
    }

    // START PREDATION
    if(msmMode > 0) {

      // ------------------------------------------------------------------------- //
      // 8.1. SUITABILITY EQUATIONS                                                //
      // ------------------------------------------------------------------------- //

      suitability.setZero();
      suma_suit.setZero();
      sum_phi.setZero();
      diet_prop_sum.setZero();

      for(rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop

        // 8.1.1. MSVPA based suitability // FIXME - not flexible for interannual variation
        if(suitMode(rsp) == 0) {

          // 8.1.1.1. Calculate stomach proportion over biomass; U/ (W * N)
          for(yr = 0; yr < nyrs; yr++) {                         // Year loop
            for(ksp = 0; ksp < nspp; ksp++) {                    // Prey species loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                  for(r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                    for(k_age = 0; k_age < nages(ksp); k_age++) {  // Prey age loop

                      stom_div_bio(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 0.0;

                      if(yr < nyrs_hind){
                        yr_ind = yr;
                      }
                      if(yr >= nyrs_hind){
                        yr_ind = nyrs_hind - 1;
                      }

                      stom_div_bio(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = diet_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / (avgN_at_age(ksp, k_sex, k_age, yr));

                      // Make into Type 3 MSVPA
                      if(msmMode == 2){
                        stom_div_bio(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) /= avgN_at_age(ksp, k_sex, k_age, yr);
                      }

                      // Make sure it is a real number, if not set to 0
                      stom_div_bio(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 0.0;


                      if(weight(pop_wt_index(ksp), k_sex, k_age, yr_ind ) != 0) {
                        stom_div_bio(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) /= weight( pop_wt_index(ksp), k_sex, k_age, yr_ind );
                        suma_suit(rsp, r_sex, r_age, yr) += stom_div_bio(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Calculate sum of stom_div_bio across prey and  prey age for each predator, predator age, and year
                        diet_prop_sum(rsp, r_sex, r_age, yr) += diet_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Sum diet proportion across predator
                      }
                    }
                  }
                }
              }
            }
          }
        }


        // 8.1.1.2. Calculate suitability
        for(r_sex = 0; r_sex < nsex(rsp); r_sex++){               // Predator sex
          for(r_age = 0; r_age < nages(rsp); r_age++) {          // Predator age loop
            for(yr = 0; yr < nyrs; yr++){
              suit_other(rsp, r_sex, r_age, yr) = 1; // Initialize suitability for other prey at 1
            }
            for(ksp = 0; ksp < nspp; ksp++) {                    // Prey species loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){         // Prey sex loop
                for(k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                  for(yr = suit_styr; yr <= suit_endyr; yr++) {  // Suit year loop (over specific years)

                    // Average suitability across years
                    if(diet_prop_sum(rsp, r_sex, r_age, yr) > 0){ // If the predator has diet data
                      suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0) += stom_div_bio(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / (suma_suit(rsp, r_sex, r_age, yr ) + other_food_diet_prop(rsp, r_sex, r_age, yr));
                    }
                  }       // End year loop

                  // FIXME - Add in interannual variability here?
                  suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0) /= nyrs_suit;

                  // Fill in years
                  for(yr = 1; yr < nyrs; yr++) {                 // Year loop
                    suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0);
                  }

                  // Other suitabilitity
                  for(yr = 0; yr < nyrs; yr++) {
                    suit_other(rsp, r_sex, r_age, yr) -= suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // FIXME - include overlap indices
                  }
                }
              }
            }
          }
        } // End Holsman/MSVPA suitability





        // 8.1.2. Estimate suitability
        if(suitMode(rsp) > 0){

          // -- Transform predator-prey preference parameters
          // Adopted from https://github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp (Trijoulet et al 2020)
          // Suitability for other food = 1-sum(predator-prey preference)
          // Criteria for predator-prey preference:
          // 1. predator-prey preference > 0 (hence logs)
          // 2. sum(predator-prey preference) + vuln_other = 1
          // 3. 0 <= sum(predator-prey preference) <= 1 (hence logit transformation)
          for(ksp = 0; ksp < nspp; ksp++) {                                   // Prey loop
            sum_phi(rsp) += exp(log_phi(rsp, ksp));
          }
          for(ksp = 0; ksp < nspp; ksp++) {                                   // Prey loop
            vulnerability(rsp, ksp) = exp(log_phi(rsp, ksp))/(1+sum_phi(rsp));// multinomial logistic transformation
          }
          vulnerability_other(rsp) = 1 - vulnerability.row(rsp).sum();        // vulnerability-other=1-sum-vulnerability but transform so sum_vuln+vuln_other=1
        }



        // 8.1.3. GAMMA suitability
        if((suitMode(rsp) == 1) | (suitMode(rsp) == 2)){
          Type log_size_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))

          for(r_age = 0; r_age < nages(rsp); r_age++) {             // Pred age
            for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
              for(ksp = 0; ksp < nspp; ksp++) {                     // Prey loop
                for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                  for(k_age = 0; k_age < nages(ksp); k_age++) {     // Prey age
                    for(yr = 0; yr < nyrs; yr++) {                  // Year loop

                      suit_other(rsp, r_sex, r_age, yr) = vulnerability_other(rsp);

                      // Year ind
                      if(yr < nyrs_hind){
                        yr_ind = yr;
                      }
                      if(yr >= nyrs_hind){
                        yr_ind = nyrs_hind - 1;
                      }

                      switch(suitMode(rsp)){
                      case 1: // Length-based GAMMA suitability
                        // log_size_ratio = log(laa( rsp, r_sex, r_age, yr_ind) / laa( ksp, k_sex, k_age, yr_ind) ); // Log ratio of lengths
                        break;
                      case 2: // Weight-based GAMMA suitability
                        log_size_ratio = log(weight(pop_wt_index(rsp), r_sex, r_age, yr_ind) / weight(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // Log ratio of weights
                        break;
                      }
                      if(log_size_ratio > 0){
                        // See line 452 from https://github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp
                        suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = vulnerability(rsp, ksp) * dgamma(log_size_ratio, gam_a( rsp ), gam_b(rsp)) / dgamma((gam_a(rsp)-1) * gam_b(rsp), gam_a(rsp), gam_b(rsp)); // Scale to 0,1 by dividing by max
                      }
                    }
                  }
                }
              }
            }
          }
        } // End gamma suitability


        // 8.1.4. Lognormal suitability
        if((suitMode(rsp) == 3) | (suitMode(rsp) == 4) | (suitMode(rsp) == 5)| (suitMode(rsp) == 6)){
          Type log_size_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
          for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
            for(r_age = 0; r_age < nages(rsp); r_age++) {                  // Pred age
              for(ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                  for(k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                    for(yr = 0; yr < nyrs; yr++){

                      suit_other(rsp, r_sex, r_age, yr) = vulnerability_other(rsp);

                      // Year ind
                      if(yr < nyrs_hind){
                        yr_ind = yr;
                      }
                      if(yr >= nyrs_hind){
                        yr_ind = nyrs_hind - 1;
                      }

                      switch(suitMode(rsp)){
                      case 3: // Length-based lognormal suitability
                        // log_size_ratio = log(laa( rsp, r_sex, r_age, yr_ind) / laa( ksp, k_sex, k_age, yr_ind) ); // Log ratio of lengths
                        break;
                      case 4: // Weight-based lognormal suitability
                        log_size_ratio = log(weight(pop_wt_index(rsp), r_sex, r_age, yr_ind) / weight(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // Log ratio of weights
                        break;

                      case 5: // Length-based normal suitability
                        // log_size_ratio = laa( rsp, r_sex, r_age, yr_ind) / laa( ksp, k_sex, k_age, yr_ind); // Log ratio of lengths
                        break;
                      case 6: // Weight-based normal suitability
                        log_size_ratio = weight(pop_wt_index(rsp), r_sex, r_age, yr_ind) / weight(pop_wt_index(ksp), k_sex, k_age, yr_ind); // Log ratio of weights
                        break;
                      }
                      // if prey are smaller than predator:
                      if(log_size_ratio > 0){
                        suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = vulnerability(rsp, ksp) * dnorm(log_size_ratio, gam_a(rsp), gam_b(rsp)) / dnorm(gam_a(rsp), gam_a(rsp), gam_b(rsp)); // Divide by mode to scale to 1
                      }
                    }
                  }
                }
              }
            }
          }
        }// End lognormal selectivity
      } // End  suitability estimation

      // ------------------------------------------------------------------------- //
      // 8. PREDATION MORTALITY EQUATIONS                                          //
      // ------------------------------------------------------------------------- //
      // -- 8.1. HOLSMAN PREDATION MORTALITY
      if((msmMode == 1) | (msmMode == 2)) {

        // 8.1.3. Calculate available food
        avail_food.setZero();
        for(rsp = 0; rsp < nspp; rsp++) {                        // Predator species loop
          for(r_sex = 0; r_sex < nsex(rsp); r_sex++){             // Predator sex loop
            for(r_age = 0; r_age < nages(rsp); r_age++) {          // Predator age loop
              for(yr = 0; yr < nyrs; yr++) {                       // Year loop

                if(yr < nyrs_hind){
                  yr_ind = yr;
                }
                if(yr >= nyrs_hind){
                  yr_ind = nyrs_hind - 1;
                }

                for(ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                    for(k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                      avail_food(rsp, r_sex, r_age, yr) += suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * pow(avgN_at_age(ksp, k_sex, k_age, yr), msmMode) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) ; // FIXME - include overlap indices: FIXME - mn_wt_stom?
                    }
                  }
                }
                // Other food
                avail_food(rsp, r_sex, r_age, yr) += other_food(rsp) * suit_other(rsp, r_sex, r_age, yr);
              }
            }
          }
        }


        // 8.1.3. Calculate predation mortality, diet proportion, and biomass easten
        M2_at_age.setZero();
        B_eaten_as_prey.setZero();
        for(ksp = 0; ksp < nspp; ksp++) {                        // Prey species loop
          for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
            for(k_age = 0; k_age < nages(ksp); k_age++) {          // Prey age loop
              for(rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop
                for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                  for(r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                    for(yr = 0; yr < nyrs; yr++) {                       // Year loop

                      if(yr < nyrs_hind){
                        yr_ind = yr;
                      }
                      if(yr >= nyrs_hind){
                        yr_ind = nyrs_hind - 1;
                      }

                      if(avail_food(rsp, r_sex, r_age, yr) > 0){

                        //switch(msmMode){
                        //case 1: // Type 2 MSVPA
                        M2_at_age(ksp, k_sex, k_age, yr) += (avgN_at_age(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr); // #FIXME - include indices of overlap
                        //M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = (avgN_at_age(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr);
                        B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = avgN_at_age(ksp, k_sex, k_age, yr) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) * avgN_at_age(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                        B_eaten_as_prey(ksp, k_sex, k_age, yr) += B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr);
                        diet_prop_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = avgN_at_age(ksp, k_sex, k_age, yr) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)/ avail_food(rsp, r_sex, r_age, yr);
                        //break;

                        /*
                      case 2: // Type 3 MSVPA
                         M2_at_age(ksp, k_sex, k_age, yr) += pow(avgN_at_age(ksp, k_sex, k_age, yr) , msmMode) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (avgN_at_age(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr)
                         * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / (avgN_at_age(ksp, k_sex, k_age, yr) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // #FIXME - include indices of overlap
                         M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pow(avgN_at_age(ksp, k_sex, k_age, yr), msmMode) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (avgN_at_age(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / (avgN_at_age(ksp, k_sex, k_age, yr) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind));
                         B_eaten_as_prey(ksp, k_sex, k_age, yr) += pow(avgN_at_age(ksp, k_sex, k_age, yr), msmMode) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) * avgN_at_age(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                         B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pow(avgN_at_age(ksp, k_sex, k_age, yr), msmMode) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) *  avgN_at_age(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                         diet_prop_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pow(avgN_at_age(ksp, k_sex, k_age, yr) , msmMode) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr) / (avgN_at_age(ksp, k_sex, k_age, yr) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind));
                         break;
                         }
                         */
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } // End 8..1. Holsman predation mortality


      // 8.2. KINZEY PREDATION EQUATIONS
      /*
       if(msmMode > 2) {

       // 8.2.3. Initialize counters
       Type Pred_ratio = 0.0;          // Predator ratio
       Type Prey_ratio = 0.0;          // Prey ratio
       Type NS_Z = 0.0;                // N(k,yr,a) * survival/Z = 0.0;
       Type Tmort = 0.0;               // Mortality on other
       Type Q_ksum_l = 0.0;            // Diet sum
       Type Term = 0.0;                // Linear adjustment for predation


       // 8.2.4. Calculate equilibrium N predators and prey in styr_pred for each species X age: FIXME: May want to have this be the final year of a projection!
       N_pred_eq.setZero();
       N_prey_eq.setZero();
       for(rsp = 0; rsp < nspp; rsp++) {
       for(ksp = 0; ksp < nspp; ksp++) {
       for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
       for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
       for(r_age = 0; r_age < nages(rsp); r_age++) {
       for(k_age = 0; k_age < nages(ksp); k_age++) {
       N_pred_eq(rsp, r_sex, r_age) += N_at_age(rsp, r_sex, r_age, 0) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0); // Denominator of Eq. 17 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
       N_prey_eq(ksp, k_sex, k_age) += N_at_age(ksp, k_sex, k_age, 0) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0); // Denominator of Eq. 16 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
       }
       }
       }
       }
       }
       }


       // 8.2.5. Calculate available prey and predator for each year
       N_pred_yrs.setZero();
       N_prey_yrs.setZero();
       for(yr = 0; yr < nyrs; yr++) {
       for(rsp = 0; rsp < nspp; rsp++) {
       for(ksp = 0; ksp < nspp; ksp++) {
       for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
       for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
       for(r_age = 0; r_age < nages(rsp); r_age++) {
       for(k_age = 0; k_age < nages(ksp); k_age++) {
       N_pred_yrs(rsp, r_sex, r_age, yr) += N_at_age(rsp, r_sex, r_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Numerator of Eq. 17 Kinzey and Punt (2009) 1st year // FIXME: Use averageN?
       N_prey_yrs(ksp, k_sex, k_age, yr) += N_at_age(ksp, k_sex, k_age, yr) * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Numerator of Eq. 16 Kinzey and Punt (2009) 1st year
       }
       }
       }
       }
       }
       }
       }


       // 8.2.6. Calculate predator functional response (Table 1 Kinzey and Punt (2009))
       for(rsp = 0; rsp < nspp; rsp++) {                          // Predator loop
       for(ksp = 0; ksp < (nspp + 1); ksp++) {                  // Prey loop
       Term = 1.0e-10 + H_1(rsp, ksp) * (Type(1) + H_1a(rsp) * H_1b(rsp) / (Type(r_age) + H_1b(rsp) + Type(1.0e-10))); // Eq. 15 Kinzey and Punt (2009)
       for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
       for(r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
       for(yr = 0; yr < nyrs; yr++) {                         // Year loop

       // Observed species
       if(ksp < nspp){
       for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
       for(k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop

       // Predator-prey ratios
       Pred_ratio = (N_pred_yrs(rsp, r_sex, r_age, yr) + Type(1.0e-10)) / (N_pred_eq(rsp, r_sex, r_age) + Type(1.0e-10)); // Eq. 17 Kinzey and Punt (2009): Predator biomass relative to equilibrium
       Prey_ratio = (N_prey_yrs(ksp, k_sex, k_age, yr) + Type(1.0e-10)) / (N_prey_eq(ksp, k_sex, k_age) + Type(1.0e-10)); // Eq. 16 Prey Kinzey and Punt (2009): biomass relative to equilibrium
       Pred_r(rsp, r_sex, r_age, yr) = Pred_ratio;
       Prey_r(ksp, k_sex, k_age, yr) = Prey_ratio;

       switch (msmMode) {
    case 3: // Holling Type I (linear)
       pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 1.0e-10 + Term;
       break;
    case 4: // Holling Type II
       pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
       ( 1 + H_2(rsp, ksp) * Prey_ratio + 1.0e-10);
       break;
    case 5: // Holling Type III
       pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 1.0e-10 +
       Term * (1 + H_2(rsp, ksp)) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) /
       (1 + H_2(rsp, ksp) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) + 1.0e-10 );
       break;
    case 6: // predator interference
       pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
       ( 1 + H_2(rsp, ksp) * Prey_ratio + H_3(rsp, ksp) * (Pred_ratio - 1) + 1.0e-10);
       break;
    case 7: // predator preemption
       pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
       ( (1 + H_2(rsp, ksp) * Prey_ratio) * (1 + H_3(rsp, ksp) * (Pred_ratio - 1)) + Type(1.0e-10));
       break;
    case 8: // Hassell-Varley
       pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 1.0e-10 + Term * (2 + H_2(rsp, ksp) + 1.0e-10) /
       (1.0 + H_2(rsp, ksp) * Prey_ratio + pow((Prey_ratio + Type(1.0e-10)), H_4(rsp, ksp)) + 1.0e-10 );
       break;
    case 9: // Ecosim
       pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 1.0e-10 + Term /
       (1 + H_3(rsp, ksp) * (Pred_ratio - 1 + 1.0e-10));
       break;
    default:
       error("Invalid 'msmMode'");
       }
       }
       }
       }
       // "other food" is linear
       else{
       pred_resp(rsp + (nspp * r_sex), (nspp * 2), r_age, 0, yr)  = 1.0e-10 + Term;
       } // end of r_ages, k_ages loop
       }   // =========================
       }
       }
       }
       }


       // 8.2.7  Predation mortality: Eq. 6 & 7 Kinzey and Punt (2009)
       M2_at_age.setZero();
       for(yr = 0; yr < nyrs; yr++) {
       for(rsp = 0; rsp  <  nspp; rsp++) {
       for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
       for(ksp = 0; ksp  <  nspp; ksp++) {
       for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
       for(r_age = 0; r_age < nages(rsp); r_age++) {
       for(k_age = 0; k_age < nages(ksp); k_age++) {
       M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)
       * suitability(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * N_at_age(rsp, r_sex, r_age, yr);
       M2_at_age(ksp, k_sex, k_age, yr) += M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr);
       }
       }
       }
       }
       }
       }
       }

       // 8.2.9. Numbers and mass eaten (of modeled prey species and "other prey"); Equations 8 and 9 from Kinzey and Punt 2009
       for(yr = 0; yr < nyrs; yr++) {
       for(rsp = 0; rsp  <  nspp; rsp++) {
       for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
       for(r_age = 0; r_age  <  nages(rsp); r_age++) {
       for(ksp = 0; ksp  <  nspp+1; ksp++) {
       // Species included
       if(ksp < nspp) {
       for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
       for(k_age = 0; k_age < nages(ksp); k_age++) {
       // Year indices
       if(yr < nyrs_hind){
       yr_ind = yr;
       }
       if(yr >= nyrs_hind){
       yr_ind = nyrs_hind - 1;
       }
       N_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * avgN_at_age(ksp, k_sex, k_age, yr); // Eq. 8 Kinzey and Punt (2009)
       B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * avgN_at_age(ksp, k_sex, k_age, yr) * weight(pop_wt_index(ksp), k_sex, k_age, yr_ind);
       B_eaten_as_pred(rsp, r_sex, r_age, yr) += B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr);
       B_eaten_as_prey(ksp, k_sex, k_age, yr) += B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Results by all species: Eq. 11 Kinzey and Punt (2009)
       }
       }
       // Other food
       }else{
       B_eaten_as_pred(rsp, r_sex, r_age, yr)  += other_food(rsp) * (Type(1) - exp(-pred_resp(rsp + (nspp * r_sex), nspp*2, r_age, 0, yr)  * N_at_age(rsp, r_sex, r_age, yr)));      // Eq. 11 Kinzey and Punt (2009)
       }
       }
       }
       }
       }
       }


       // 8.2.11. Predicted diet proportion as weight of prey-at-age in diet of predator-at-age (Eqn 15)
       // NOTE: Only including prey species in the model. Does not include other food
       for(rsp = 0; rsp  <  nspp; rsp++) {
       for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
       for(r_age = 0; r_age  <  nages(rsp); r_age++) {
       for(ksp = 0; ksp  <  nspp; ksp++) {
       for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
       for(k_age = 0; k_age < nages(ksp); k_age++) {
       for(yr = 0; yr < nyrs_hind; yr++) {
       diet_prop_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)/B_eaten_as_pred(rsp, r_sex, r_age, yr);
       }
       }
       }
       }
       }
       }
       }
       }


       // 8.2.12. Predict ration Eq. 13 Kinzey and Punt (2009)
       Type numer;
       Type denom;

       ration_hat_ave.setZero();
       ration_hat.setZero();
       for(rsp = 0; rsp < nspp; rsp++) {
       for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
       for(r_age = 0; r_age < nages(rsp); r_age++) {
       numer = 0;
 denom = 0;
 // Calculate year-specific values
 for(yr = 0; yr < nyrs; yr++) {
 if(avgN_at_age(rsp, r_sex, r_age, yr) > 0){
 ration_hat(rsp, r_sex, r_age, yr) = B_eaten_as_pred(rsp, r_sex, r_age, yr)/avgN_at_age(rsp, r_sex, r_age, yr); // NOTE: Divide by 365 to make into daily ration
 }
 }

 // Find average over hindcast// FIXME: suit_years?
 for(yr = 0; yr < nyrs_hind; yr++) {
 numer += B_eaten_as_pred(rsp, r_sex, r_age, yr); // NOTE: Divide by 365 to make into daily ration
 denom += avgN_at_age(rsp, r_sex, r_age, yr);
 }
 ration_hat_ave(rsp, r_sex, r_age) = numer / denom;
       }
       }
       }
 // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
 } // End 8.2. Kinzey predation
 // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
 */
    } // End 8. Predation mortality
    // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //



    // ------------------------------------------------------------------------- //
    // 9. INDEX COMPONENTS EQUATIONS                                             //
    // ------------------------------------------------------------------------- //

    // -- 9.1. Index of abundance/biomass
    for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

      index = index_ctl(index_ind, 0) - 1;            // Temporary survey index
      sp = index_ctl(index_ind, 1) - 1;             // Temporary index of species
      flt_yr = index_ctl(index_ind, 2);             // Temporary index for years of data

      mo = index_n(index_ind, 0);                   // Temporary index for month

      index_hat(index_ind) = 0.0;                    // Initialize

      if(flt_yr > 0){
        flt_yr = flt_yr - styr;
      }
      if(flt_yr < 0){
        flt_yr = -flt_yr - styr;
      }

      // Hindcast
      if(flt_yr < nyrs_hind){
        yr_ind = flt_yr;
      }

      // Projection
      if(flt_yr >= nyrs_hind){
        yr_ind = nyrs_hind - 1;
      }

      for(age = 0; age < nages(sp); age++) {
        for(sex = 0; sex < nsex(sp); sex++){
          // Weight
          if(flt_units(index) == 1){
            index_hat(index_ind) += N_at_age(sp, sex, age, flt_yr) * exp( - (mo/12.0) * Z_at_age(sp, sex, age, flt_yr)) * sel(index, sex, age, yr_ind) * weight( flt_wt_index(index), sex, age, yr_ind );
          }
          // Numbers
          if(flt_units(index) == 2){
            index_hat(index_ind) += N_at_age(sp, sex, age, flt_yr) * exp( - (mo/12.0) * Z_at_age(sp, sex, age, flt_yr)) * sel(index, sex, age, yr_ind);
          }
        }
      }
    }


    // -- 9.2. Analytical survey q following Ludwig and Martell 1994
    index_n_obs.setZero();
    index_q_analytical.setZero();
    for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){



      index = index_ctl(index_ind, 0) - 1;            // Temporary survey index
      sp = index_ctl(index_ind, 1) - 1;             // Temporary index of species
      flt_yr = index_ctl(index_ind, 2);             // Temporary index for years of data

      if(flt_yr > 0){
        flt_yr = flt_yr - styr;


        mo = index_n(index_ind, 0);                    // Temporary index for month
        if(flt_yr < nyrs_hind){

          // If etimated standard deviation or analytical sigma (non-time-varying)
          if(est_sigma_index(index) > 0) {
            index_n_obs(index) += 1; // Add one if survey is used
            index_q_analytical(index) += log(index_obs(index_ind, 0) / index_hat(index_ind));
          }

          // If time-varying sigma
          if(est_sigma_index(index) == 0 ) {
            index_n_obs(index) += 1 / square(index_obs(index_ind, 1));
            index_q_analytical(index) += log(index_obs(index_ind, 0) / index_hat(index_ind)) / square(index_obs(index_ind, 1));
          }
        }
      }
    }

    // Take average
    for(index = 0 ; index < n_flt; index ++){
      index_q_analytical(index) = exp(index_q_analytical(index) / index_n_obs(index));

      // Set index_q to analytical if used
      if(est_index_q(index) == 3){
        for(yr = 0; yr < nyrs_hind; yr++){
          index_q(index, yr) = index_q_analytical(index);
        }
      }
    }


    // -- 9.3. Survey Biomass - multiply by q
    for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

      index = index_ctl(index_ind, 0) - 1;            // Temporary survey index
      flt_yr = index_ctl(index_ind, 2);      // Temporary index for years of data

      if(flt_yr > 0){
        flt_yr = flt_yr - styr;
      }
      if(flt_yr < 0){
        flt_yr = -flt_yr - styr;
      }

      // Hindcast
      if(flt_yr < nyrs_hind){
        yr_ind = flt_yr;
      }

      // Projection
      if(flt_yr >= nyrs_hind){
        yr_ind = nyrs_hind - 1;
      }

      index_hat(index_ind) = index_q(index, yr_ind) * index_hat(index_ind); // pow(index_hat(index_ind), (1 + index_q_pow(index)));

    }



    // -- 9.4. Calculate analytical sigma following Ludwig and Walters 1994
    index_n_obs.setZero();
    ln_index_analytical_sd.setZero();
    for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

      index = index_ctl(index_ind, 0) - 1;            // Temporary survey index
      flt_yr = index_ctl(index_ind, 2);      // Temporary index for years of data

      if(flt_yr > 0){
        flt_yr = flt_yr - styr;

        if(flt_yr < nyrs_hind){
          index_n_obs(index) += 1; // Add one if survey is used
          ln_index_analytical_sd(index) += square( log(index_obs(index_ind, 0)) - log(index_hat(index_ind)));
        }
      }
    }

    for(index = 0 ; index < n_flt; index ++){
      ln_index_analytical_sd(index) = sqrt(ln_index_analytical_sd(index) / index_n_obs(index));
    }


    // ------------------------------------------------------------------------- //
    // 10. FISHERY COMPONENTS EQUATIONS                                          //
    // ------------------------------------------------------------------------- //
    // 10.1. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
    for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){

      flt = catch_ctl(fsh_ind, 0) - 1;     // Temporary fishery index
      sp = catch_ctl(fsh_ind, 1) - 1;      // Temporary index of species
      flt_yr = catch_ctl(fsh_ind, 2);      // Temporary index for years of data
      mo = catch_n(fsh_ind, 0);            // Temporary index for month

      catch_hat(fsh_ind) = 0.0;  // Initialize
      max_catch_hat(fsh_ind) = 0.0; // Initialize

      if(flt_yr > 0){
        flt_yr = flt_yr - styr;
      }
      if(flt_yr < 0){
        flt_yr = -flt_yr - styr;
      }

      // Hindcast
      if(flt_yr < nyrs_hind){
        yr_ind = flt_yr;
      }

      // Projection
      if(flt_yr >= nyrs_hind){
        yr_ind = nyrs_hind - 1;
      }

      for(sex = 0; sex < nsex(sp); sex ++){
        for(age = 0; age < nages(sp); age++) {
          // FIXME: add wt_yr_ind for MSE

          // By weight
          if(flt_units(flt) == 1){
            catch_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Z_at_age(sp, sex, age, flt_yr) * (1.0 - exp(-Z_at_age(sp, sex, age, flt_yr))) * N_at_age(sp, sex, age, flt_yr) * weight( flt_wt_index(flt), sex, age, yr_ind ); // 5.5.
            max_catch_hat(fsh_ind) += N_at_age(sp, sex, age, flt_yr) * weight( flt_wt_index(flt), sex, age, yr_ind ) * sel(flt, sex, age, yr_ind) * proj_F_prop(flt); // FIXME using last year of selectivity;
          }

          // By numbers
          if(flt_units(flt) == 2){
            catch_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Z_at_age(sp, sex, age, flt_yr) * (1.0 - exp(-Z_at_age(sp, sex, age, flt_yr))) * N_at_age(sp, sex, age, flt_yr);
            max_catch_hat(fsh_ind) += N_at_age(sp, sex, age, flt_yr) * sel(flt, sex, age, yr_ind) * proj_F_prop(flt); // FIXME using last year of selectivity;
          }
        }
      }
    }

    // 10.2 Exploitable biomass ----
    exploitable_biomass.setZero();

    for(flt = 0; flt < n_flt; flt++) {

      sp = flt_spp(flt);

      if(flt_type(flt) == 1){ //Fishery only
        for(yr = 0; yr < nyrs; yr++) {

          // Hindcast
          if(yr < nyrs_hind){
            yr_ind = yr;
          }

          // Projection
          if(yr >= nyrs_hind){
            yr_ind = nyrs_hind - 1;
          }

          for(age = 0; age < nages(sp); age++) {
            for(sex = 0; sex < nsex(sp); sex ++){
              exploitable_biomass(sp, yr) += N_at_age(sp, sex, age, yr) * weight( flt_wt_index(flt), sex, age, yr_ind ) * sel(flt, sex, age, yr_ind) * proj_F_prop(flt); // FIXME using last year of selectivity;
            }
          }
        }
      }
    }


    // ------------------------------------------------------------------------- //
    // 11. COMPOSITION EQUATIONS                                                  //
    // ------------------------------------------------------------------------- //

    // -- 11.1. Composition
    age_obs_hat.setZero();
    comp_hat.setZero();
    age_hat.setZero();
    for(comp_ind = 0; comp_ind < comp_hat.rows(); comp_ind++){

      flt = comp_ctl(comp_ind, 0) - 1;            // Temporary fishery index
      sp = comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
      flt_sex = comp_ctl(comp_ind, 2);            // Temporary index for comp sex (0 = combined, 1 = female, 2 = male)
      comp_type = comp_ctl(comp_ind, 3);          // Temporary index for comp type (0 = age, 1 = length, 2 = CAAL)
      yr = comp_ctl(comp_ind, 4);                 // Temporary index for years of data
      mo = comp_n(comp_ind, 0);                   // Temporary index for month

      n_hat(comp_ind) = 0.0;                          // Initialize

      // Likelihood (yr > 0) vs prediction (yr < 0)
      if(yr > 0){
        yr = yr - styr;
      }
      if(yr < 0){
        yr = -yr - styr;
      }

      // Hindcast
      if(yr < nyrs_hind){
        yr_ind = yr;
      }

      // Projection
      if(yr >= nyrs_hind){
        yr_ind = nyrs_hind - 1;
      }

      // Get catch-at-age
      for(age = 0; age < nages(sp); age++) {


        switch(flt_type(flt)){
        case 1: // - Fishery
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              // Catch-at-age
              age_hat(comp_ind, age ) += F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
              // Total numbers
              n_hat( comp_ind ) += F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            // Catch-at-age
            age_hat(comp_ind, age ) = F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            // Total numbers
            n_hat( comp_ind ) += F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
              // Total numbers
              n_hat(comp_ind) += F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            }
            break;
          }
          break;


        case 2: // - Survey
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age ) += N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
              // Total numbers
              n_hat(comp_ind) += N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));  // FIXME - put power in for catchability
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            // Survey catch-at-age
            age_hat(comp_ind, age ) = N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            // Total numbers
            n_hat(comp_ind) += N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
              // Total numbers
              n_hat(comp_ind) += N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            }
            break;
          }
        }
      }

      // Adjustment for joint sex composition data
      joint_adjust(comp_ind) = 1;
      if(flt_sex == 3){
        joint_adjust(comp_ind) = 2;
      }

      // True age comp standardize to sum to 1
      for(age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
        true_age_comp_hat(comp_ind, age ) = age_hat(comp_ind, age ) / posfun(n_hat(comp_ind), Type(0.0001), penalty);
      }


      // Adjust for aging error
      for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
        for(int true_age = 0; true_age < nages(sp); true_age++) {
          age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(sp, true_age, obs_age);
        }
      }

      // Adjust for aging error for joint data
      if(flt_sex == 3){
        for(int obs_age = nages(sp); obs_age < nages(sp) * 2; obs_age++) {
          for(int true_age = nages(sp); true_age < nages(sp) * 2; true_age++) {

            // Adjust indexing for joint age/length comp
            int true_age_tmp = true_age - nages(sp);
            int obs_age_tmp = obs_age - nages(sp);

            age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(sp, true_age_tmp, obs_age_tmp);
          }
        }
      }


      //  Observed age comp standardize to sum to 1
      if(comp_type == 0) {
        for(age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
          comp_hat(comp_ind, age) = age_obs_hat(comp_ind, age ) / posfun(n_hat(comp_ind), Type(0.0001), penalty);
        }
      }


      // Convert from catch-at-age to catch-at-length
      if( comp_type == 1) {
        for(ln = 0; ln < nlengths(sp); ln++) {
          for(age = 0; age < nages(sp); age++) {

            sex = 0;

            // Adjust sex for males/females
            if((flt_sex > 0) & (flt_sex < 3)){
              sex = flt_sex - 1;
            }

            comp_hat(comp_ind, ln ) += age_obs_hat(comp_ind, age) * age_trans_matrix(flt_age_transition_index(flt), sex, age, ln );
          }
        }


        // Convert from catch-at-age to catch-at-length for joint comp data
        if( flt_sex == 3) {
          for(ln = nlengths(sp); ln < nlengths(sp) * 2; ln++) {
            for(age = nages(sp); age < nages(sp) * 2; age++) {

              // Adjust indexing for joint age/length comp
              int obs_ln_tmp = ln;
              int obs_age_tmp = age;

              obs_age_tmp = age - nages(sp);
              obs_ln_tmp = ln - nlengths(sp);
              sex = 1;

              comp_hat(comp_ind, ln ) += age_obs_hat(comp_ind, age) * age_trans_matrix(flt_age_transition_index(flt), sex, obs_age_tmp, obs_ln_tmp );
            }
          }
        }


        // Standardize to sum to 1
        for(ln = 0; ln < nlengths(sp) * joint_adjust(comp_ind); ln++) {
          comp_hat(comp_ind, ln ) = comp_hat(comp_ind, ln) / posfun(n_hat(comp_ind), Type(0.0001), penalty);
        }
      }
    }



    // ------------------------------------------------------------------------- //
    // 12. PREDICTED STOMACH CONTENT                                             //
    // ------------------------------------------------------------------------- //

    // Predict stomach content
    // 12. Reorganize diet_hat content
    for(int stom_ind = 0; stom_ind < diet_obs.rows(); stom_ind++){
      rsp = diet_ctl(stom_ind, 0) - 1;             // Index of pred
      ksp = diet_ctl(stom_ind, 1) - 1;             // Index of prey
      r_sex = diet_ctl(stom_ind, 2);               // Index of pred sex
      k_sex = diet_ctl(stom_ind, 3);               // Index of prey sex
      r_age = diet_ctl(stom_ind, 4) - minage(rsp); // Index of pred age
      k_age = diet_ctl(stom_ind, 5) - minage(ksp); // Index of prey age
      flt_yr = diet_ctl(stom_ind, 6);              // Index of year

      // 1 sex model
      r_sexes(stom_ind, 0) = 0; r_sexes(stom_ind, 1) = 0;
      k_sexes(stom_ind, 0) = 0; k_sexes(stom_ind, 1) = 0;

      // 2 sex model
      // This is to account for situations where nsex = 2, but r_sex or k_sex = 0 should be weighted average
      if(nsex(rsp) == 2){
        // But k_sex = 0 (indicating diet data is for both sexes)
        r_sexes(stom_ind, 0) = 0; r_sexes(stom_ind, 1) = 1;

        if(r_sex > 0){
          r_sexes(stom_ind, 0) = r_sex - 1;  r_sexes(stom_ind, 1) = r_sex - 1;
        }
      }

      if(nsex(ksp) == 2){
        // But k_sex = 0 (indicating diet data is for both sexes)
        k_sexes(stom_ind, 0) = 0; k_sexes(stom_ind, 1) = 1;

        // Sex-specific diet data
        if(k_sex > 0){
          k_sexes(stom_ind, 0) = k_sex - 1;  k_sexes(stom_ind, 1) = k_sex - 1;
        }
      }

      // Initialize
      diet_hat(stom_ind, 1) = 0;

      for(int j = 0; j < 2; j ++){
        for(int k = 0; k < 2; k ++){

          if(flt_yr > 0){
            yr = flt_yr - styr;

            // Annual data
            if(yr < nyrs_hind){
              // Diet proportion of prey-at-age in predator-at-age
              if(k_age >= 0){
                diet_hat(stom_ind, 1) += diet_prop_hat(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr)/4; //FIXME: take weighted average for two-sex models?
              }


              // Diet proportion of prey-spp in predator-at-age (sum across prey ages)
              if(k_age < 0){
                for(k_age = 0; k_age < nages(ksp); k_age++){
                  diet_hat(stom_ind, 1) += diet_prop_hat(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr)/4; //FIXME: take weighted average for two-sex models?
                }
              }
            }
          }

          // Average of years
          if(flt_yr == 0){
            for(yr = suit_styr; yr <= suit_endyr; yr++) {  // Suit year loop (over specific years)

              // Diet proportion of prey-at-age in predator-at-age
              if(k_age >= 0){
                diet_hat(stom_ind, 1) += diet_prop_hat(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr)/4/nyrs_suit; //FIXME: take weighted average for two-sex models?
              }


              // Diet proportion of prey-spp in predator-at-age (sum across prey ages)
              if(k_age < 0){
                for(k_age = 0; k_age < nages(ksp); k_age++){
                  diet_hat(stom_ind, 1) += diet_prop_hat(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr)/4/nyrs_suit; //FIXME: take weighted average for two-sex models?
                }
              }
            }
          }
        }
      }
    }
    // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
  } // End population dynamics iterations
  // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //




  // ------------------------------------------------------------------------- //
  // 13. DERIVED QUANTITIES
  // ------------------------------------------------------------------------- //


  // 13.1. Depletion
  biomass_depletion.setZero();
  ssb_depletion.setZero();
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs; yr++){

      if(DynamicHCR == 0){
        biomass_depletion(sp, yr) = biomass(sp, yr)/B0(sp, nyrs-1);
        ssb_depletion(sp, yr) = ssb(sp, yr)/SB0(sp, nyrs-1);
      }

      if(DynamicHCR == 1){
        biomass_depletion(sp, yr) = biomass(sp, yr)/DynamicB0(sp, yr);
        ssb_depletion(sp, yr) = ssb(sp, yr)/DynamicSB0(sp, yr);
      }

      // Multi-species and no HCR (MSSB0 is input otherwises)
      if((HCR == 0) & (msmMode > 0)){
        biomass_depletion(sp, yr) = biomass(sp, yr)/biomass(sp, nyrs-1);
        ssb_depletion(sp, yr) = ssb(sp, yr)/ssb(sp, nyrs-1);
      }
    }
  }




  // ------------------------------------------------------------------------- //
  // 14. LIKELIHOOD EQUATIONS
  // ------------------------------------------------------------------------- //

  // 14.0. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(19, n_flt); jnll_comp.setZero();  // matrix of negative log-likelihood components

  // -- Data likelihood components (Fleet specific)
  // Slot 0 -- Survey biomass
  // Slot 1 -- Total catch (kg)
  // Slot 2 -- Age/length composition
  // Slot 3 -- Non-parametric selectivity
  // Slot 4 -- Selectivity annual deviates
  // Slot 5 -- Survey catchability prior
  // Slot 6 -- Survey catchability annual deviates
  // Slot 7 -- EMPTY
  // -- Species-specific priors/penalties
  // Slot 8 -- Stock recruitment parameter prior
  // Slot 9 -- init_dev -- Initial abundance-at-age
  // Slot 10 -- Annual recruitment deviation
  // Slot 11 -- Stock-recruit penalty
  // Slot 12 -- Reference point penalities
  // Slot 13 -- N-at-age < 0 penalty
  // Slot 14 -- M_at_age prior
  // Slot 15 -- M_at_age random effects
  // -- M2_at_age likelihood components
  // Slot 16 -- Ration likelihood
  // Slot 17 -- Ration penalties
  // Slot 18 -- Diet proportion by weight likelihood


  // 14.1. FIT OBJECTIVE FUNCTION
  // Slot 0 -- Index of abundance biomass/numbers --
  Type index_std_dev = 0;
  for(index_ind = 0; index_ind < index_obs.rows(); index_ind++) {

    index = index_ctl(index_ind, 0) - 1;            // Temporary survey index
    sp = index_ctl(index_ind, 1) - 1;             // Species is the 2nd column
    flt_yr = index_ctl(index_ind, 2);             // Temporary index for years of data

    if(flt_yr > 0){
      // Set up variance
      switch (est_sigma_index(index)) {
      case 0:     // Provided standard deviation
        index_std_dev = index_obs(index_ind, 1);
        break;
      case 1:     // Estimated standard deviation
        index_std_dev = exp(index_ln_sd(index));
        break;
      case 2:     // Analytical
        index_std_dev = ln_index_analytical_sd(index);
        break;
      default:
        error("Invalid 'Estimate_sigma_index'");
      }

      ln_index_sd(index_ind) = index_std_dev;


      // Only include years from hindcast
      if(flt_type(index) > 0){
        if(flt_yr <= endyr){
          if(index_obs(index_ind) > 0){
            jnll_comp(0, index) -= dnorm(log(index_obs(index_ind, 0)), log(index_hat(index_ind)) - square(index_std_dev)/2.0, index_std_dev, true);
          }
        }
      }
    }
  }


  // Slot 1 -- Total catch --
  Type fsh_std_dev = 0;
  for(fsh_ind = 0; fsh_ind < catch_obs.rows(); fsh_ind++) {

    flt = catch_ctl(fsh_ind, 0) - 1;            // Temporary fishery index
    sp = catch_ctl(fsh_ind, 1) - 1;             // Species is the column 3
    flt_yr = catch_ctl(fsh_ind, 2);             // Temporary index for years of data
    yr = flt_yr - styr;                            // Temporary index of years. Start at 0.

    if(flt_yr > 0){

      // Set up variance
      switch (est_sigma_fsh(flt)) {
      case 0: // Provided standard deviation
        fsh_std_dev = catch_obs(fsh_ind, 1);
        break;
      case 1:     // Estimated standard deviation
        fsh_std_dev = exp(catch_ln_sd(flt));
        break;
      default:
        error("Invalid 'Estimate_sigma_catch'");
      }

      ln_catch_sd(fsh_ind) = fsh_std_dev; // Save estimated log_sd

      // Add only years from hindcast
      if(flt_type(flt) == 1){
        if(flt_yr <= endyr){
          if(catch_obs(fsh_ind, 0) > 0){
            jnll_comp(1, flt) -= dnorm(log(catch_obs(fsh_ind, 0)), log(catch_hat(fsh_ind)) - square(fsh_std_dev)/2.0, fsh_std_dev, true) ;

            // Martin's
            // jnll_comp(1, flt)+= 0.5*square((log(catch_obs(fsh_ind, 0))-log(catch_hat(fsh_ind)))/fsh_std_dev);
          }
        }
      }
    }
  }



  // Slot 2 -- Age/length composition
  for(comp_ind = 0; comp_ind < comp_obs.rows(); comp_ind++) {

    flt = comp_ctl(comp_ind, 0) - 1;        // Temporary fleet index
    sp = comp_ctl(comp_ind, 1) - 1;         // Temporary index of species
    flt_sex = comp_ctl(comp_ind, 2);        // Temporary index for comp sex (0 = combined, 1 = female, 2 = male, 3 = joint)
    comp_type = comp_ctl(comp_ind, 3);      // Temporary index for comp type (0 = age, 1 = length)
    yr = comp_ctl(comp_ind, 4);             // Temporary index for years of data

    // Adjustment for joint sex composition data
    joint_adjust(comp_ind) = 1;
    if(flt_sex == 3){
      joint_adjust(comp_ind) = 2;
    }

    // Number of ages
    int n_comp = 0;
    if(comp_type == 0){
      n_comp = nages(sp) * joint_adjust(comp_ind); // For sex = 3 (joint sex comp data)

    }
    if(comp_type == 1){
      n_comp = nlengths(sp) * joint_adjust(comp_ind);
    }

    // Select sections
    vector<Type> comp_obs_tmp = comp_obs.row(comp_ind).segment(0, n_comp); // Observed proportion
    vector<Type> comp_hat_tmp = comp_hat.row(comp_ind).segment(0, n_comp); // Expected proportion

    // Add offset (for some reason can't do above in single line....)
    comp_obs_tmp += 0.00001;
    comp_hat_tmp += 0.00001;

    // Convert observed prop to observed numbers
    comp_obs_tmp *= comp_n(comp_ind, 1);
    vector<Type> alphas = sum(comp_obs_tmp) * comp_hat_tmp * DM_pars(flt); // DM alpha

    // Only use years wanted
    if((yr <= endyr) & (yr > 0) & (flt_type(flt) > 0)){

      switch(comp_ll_type(flt)){

      case -1:
        for(ln = 0; ln < n_comp; ln++) {
          // Martin's
          jnll_comp(2, flt) -= comp_weights(flt) * Type(comp_n(comp_ind, 1)) * (comp_obs(comp_ind, ln) + 0.00001) * log((comp_hat(comp_ind, ln)+0.00001) / (comp_obs(comp_ind, ln) + 0.00001)) ;
        }
        break;

      case 0:  // Full multinomial
        jnll_comp(2, flt) -= comp_weights(flt) * dmultinom(comp_obs_tmp, comp_hat_tmp, true);
        break;

      case 1:  // Dirichlet-multinomial
        jnll_comp(2, flt) -= ddirmultinom(comp_obs_tmp, alphas,  true);
        break;

      default:
        error("Invalid 'comp_ll_type'");
      }
    }
  }


  // Slot 3-4 -- Selectivity
  for(flt = 0; flt < n_flt; flt++){ // Loop around surveys
    jnll_comp(3, flt) = 0;
    jnll_comp(4, flt) = 0;
    sp = flt_spp(flt);

    // If estimating survey or fishery
    if(flt_type(flt) > 0){

      // Ianelli non-parametic selectivity penalties
      // - using non-normalized selectivities following the arrowtooth ADMB model
      // - updated to make differentiable using abs to only penalize when sel_ratio_tmp > 0 (decreasing sel)
      if(flt_sel_type(flt) == 2) {

        for(sex = 0; sex < nsex(sp); sex++){
          for(age = 0; age < (nages(sp) - 1); age++) {
            Type sel_ratio_tmp = log(non_par_sel(flt, sex, age) / non_par_sel(flt, sex, age + 1) );
            jnll_comp(3, flt) += sel_curve_pen(flt, 0) * square( (CppAD::abs(sel_ratio_tmp) + sel_ratio_tmp)/2.0);
          }
        }

        for(sex = 0; sex < nsex(sp); sex++){
          // Extract only the selectivities we want
          vector<Type> sel_tmp(nages(sp)); sel_tmp.setZero();

          for(age = 0; age < nages(sp); age++) {
            sel_tmp(age) = log(non_par_sel(flt, sex, age));
          }

          for(age = 0; age < nages(sp) - 2; age++) {
            sel_tmp(age) = first_difference( first_difference( sel_tmp ) )(age);
            jnll_comp(3, flt) += sel_curve_pen(flt, 1) * pow( sel_tmp(age) , 2);
          }
        }

        // Survey selectivity normalization (non-parametric)
        for(sex = 0; sex < nsex(sp); sex++){
          jnll_comp(3, flt) += 1.0 * square(avg_sel(flt, sex));
        }
      }



      // Penalized/random effect likelihood time-varying logistic/double-logistic selectivity deviates
      if(((flt_varying_sel(flt) == 1)|(flt_varying_sel(flt) == 2)) & (flt_sel_type(flt) != 2) & (flt_sel_type(flt) != 5)){
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 0; yr < nyrs_hind; yr++){

            jnll_comp(4, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr), Type(0.0), sel_dev_sd(flt), true);
            jnll_comp(4, flt) -= dnorm(ln_sel_slp_dev(0, flt, sex, yr), Type(0.0), 4 * sel_dev_sd(flt), true);

            // Double logistic deviates
            if(flt_sel_type(flt) == 3){
              jnll_comp(4, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(4, flt) -= dnorm(ln_sel_slp_dev(1, flt, sex, yr), Type(0.0), 4 * sel_dev_sd(flt), true);
            }
          }
        }
      }

      // Penalized/random effect likelihood time-varying non-parametric (Taylor et al 2014) selectivity deviates
      if(((flt_varying_sel(flt) == 1) | (flt_varying_sel(flt) == 2)) & (flt_sel_type(flt) == 5)){
        for(age = 0; age < flt_nselages(flt); age++){ //NOTE: extends beyond selectivity age range, but should be mapped to 0 in map function
          for(sex = 0; sex < nsex(sp); sex++){
            for(yr = 0; yr < nyrs_hind; yr++){
              jnll_comp(4, flt) -= dnorm(sel_coff_dev(flt, sex, age, yr), Type(0.0), sel_dev_sd(flt), true);
            }
          }
        }
      }


      // Random walk: Type 4 = random walk on ascending and descending for double logistic; Type 5 = ascending only for double logistics
      if(((flt_varying_sel(flt) == 4)|(flt_varying_sel(flt) == 5)) & (flt_sel_type(flt) != 2) & (flt_sel_type(flt) != 5)){
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 1; yr < nyrs_hind; yr++){ // Start at second year

            // Logistic deviates
            jnll_comp(4, flt) -= dnorm(ln_sel_slp_dev(0, flt, sex, yr) - ln_sel_slp_dev(0, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt), true);
            jnll_comp(4, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr) - sel_inf_dev(0, flt, sex, yr-1), Type(0.0), 4 * sel_dev_sd(flt), true);

            // Double logistic deviates
            if((flt_sel_type(flt) == 3) & (flt_varying_sel(flt) == 4)){
              jnll_comp(4, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr) - sel_inf_dev(1, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(4, flt) -= dnorm(ln_sel_slp_dev(1, flt, sex, yr) - ln_sel_slp_dev(1, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt) * 4, true);
            }

          }
        }
      }
    }
  } // End selectivity loop


  // Slot 5-6 -- Survey catchability deviates
  for(flt = 0; flt < n_flt; flt++){

    // Prior on catchability
    if( est_index_q(flt) == 2){
      jnll_comp(5, flt) -= dnorm(index_ln_q(flt), index_ln_q_prior(flt), index_q_sd(flt), true);
    }

    // QAR1 deviates fit to environmental index (sensu Rogers et al 2024; 10.1093/icesjms/fsae005)
    if(est_index_q(flt) == 6){

      // AR1 process error
      Type rho=rho_trans(index_q_rho(flt));
      vector<Type> index_q_dev_tmp = index_q_dev.row(flt);
      jnll_comp(5, flt) = SCALE(AR1(rho), index_q_dev_sd(flt))(index_q_dev_tmp);

      // Observation error
      // - Fit to environmental index
      int q_index = index_varying_q(flt) - 1;
      for(yr = 0; yr < nyrs_hind; yr++){
        jnll_comp(6, flt) -= dnorm(env_index(yr, q_index), index_q_dev(flt, yr), index_q_sd(flt), true); //FIXME: index by env-year
      }
    }

    // Penalized/random deviate likelihood
    // - Estimate_q = 1 or 2
    // - Time_varying_q  = 1 or 2
    if(((index_varying_q(flt) == 1) | (index_varying_q(flt) == 2)) & (flt_type(flt) == 2) & ((est_index_q(flt) == 1) | (est_index_q(flt) == 2))){
      for(yr = 0; yr < nyrs_hind; yr++){
        jnll_comp(6, flt) -= dnorm(index_q_dev(flt, yr), Type(0.0), index_q_dev_sd(flt), true );
      }
    }

    // Random walk
    // - Estimate_q = 1 or 2
    // - Time_varying_q  = 4
    if((index_varying_q(flt) == 4) & (flt_type(flt) == 2) & ((est_index_q(flt) == 1) | (est_index_q(flt) == 2))){
      for(yr = 1; yr < nyrs_hind; yr++){
        jnll_comp(6, flt) -= dnorm(index_q_dev(flt, yr) - index_q_dev(flt, yr-1), Type(0.0), index_q_dev_sd(flt), true );
      }
    }
  } // End q loop


  // Slots 8-11 -- RECRUITMENT PARAMETERS
  for(sp = 0; sp < nspp; sp++) {
    penalty = 0.0;
    // Slot 9 -- stock-recruit prior for Beverton
    // -- Lognormal
    if((srr_est_mode == 2) & ((srr_pred_fun == 2) | (srr_pred_fun == 3))){
      jnll_comp(8, sp) -= dnorm(log(steepness(sp)), log(srr_prior(sp)) + square(srr_prior_sd(sp))/2.0, srr_prior_sd(sp), true);
    }

    // -- Beta
    if((srr_est_mode == 3) & ((srr_pred_fun == 2) | (srr_pred_fun == 3))){
      // Convert mean and SD to beta params
      Type beta_alpha = ((1 - srr_prior(sp))/ square(srr_prior_sd(sp)) - 1/srr_prior(sp)) * square(srr_prior(sp));
      Type beta_beta = beta_alpha * (1/srr_prior(sp) - 1);
      jnll_comp(8, sp) -= dbeta(steepness(sp), beta_alpha, beta_beta, true);
    }

    // Slot 9 -- stock-recruit prior for Ricker
    if((srr_est_mode == 2) & ((srr_pred_fun == 4) | (srr_pred_fun == 5))){
      jnll_comp(8, sp) -= dnorm((rec_pars(sp, 1)), log(srr_prior(sp)), srr_prior_sd(sp), true);
    }

    // Slot 9 -- penalty for Bmsy > Bmsy_lim for Ricker
    if(!isNA(Bmsy_lim(sp)) & ((srr_pred_fun == 4) | (srr_pred_fun == 5))){ // Using pred_fun in case ianelli method is used
      Type bmsy = 1.0/exp(rec_pars(sp, 2));
      bmsy =  posfun(Bmsy_lim(sp)/Type(1000000.0) - bmsy, Type(0.001), penalty);
      jnll_comp(8, sp) += 100 * penalty;
    }


    // Slot 9 -- init_dev -- Initial abundance-at-age
    if(initMode > 0){
      for(age = 1; age < nages(sp); age++) {
        jnll_comp(9, sp) -= dnorm( init_dev(sp, age - 1), square(R_sd(sp))/2.0, R_sd(sp), true);
      }
    }

    // Slot 10 -- Tau -- Annual recruitment deviation
    for(yr = 0; yr < nyrs_hind; yr++) {
      jnll_comp(10, sp) -= dnorm( rec_dev(sp, yr),  square(R_sd(sp))/2.0, R_sd(sp), true);    // Recruitment deviation using random effects.
    }

    // Slot 11 -- Additional penalty for SRR curve (sensu AMAK/Ianelli)
    if((srr_fun == 0) & (srr_pred_fun  > 0)){
      for(yr = srr_hat_styr; yr <= srr_hat_endyr; yr++) {
        jnll_comp(11, sp) -= dnorm( log(R(sp, yr)), log(R_hat(sp,yr)), R_sd(sp), true);
      }
    }
  }


  // Slot 12 -- Reference point penalties
  // -- CMSY
  Type CMSY = 0;

  // --- Sum terminal catch across species
  if(HCR == 1){

    // -- Loop through catch data
    for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){

      flt = catch_ctl(fsh_ind, 0) - 1;
      sp = catch_ctl(fsh_ind, 1) - 1;
      flt_yr = catch_ctl(fsh_ind, 2);

      // Add fishery data from terminal year
      if((flt_type(flt) == 1) & (forecast(sp) == 1) & (estDynamics(sp) == 0)){
        if(flt_yr == projyr){
          CMSY  += catch_hat(fsh_ind);
        }
      }
    }

    jnll_comp(12, 0) = - square(CMSY/1000000.0); // CMSY is ll


    // --- Add biomass_depletion constraint
    for(sp = 0; sp < nspp; sp++) {
      if((forecast(sp) == 1) & (estDynamics(sp) == 0)){
        penalty = 0.0;
        Type nothing_useful =  posfun( (ssb_depletion(sp, nyrs-1) - Plimit(sp)), Type(0.0001), penalty);
        jnll_comp(12, sp) += 500.0 * square(CMSY/1000.0) * penalty; // CMSY
      }
    }
  }


  for(sp = 0; sp < nspp; sp++) {
    // -- Single-species static reference points
    if((DynamicHCR == 0) & (forecast(sp) == 1) & (msmMode == 0) & (estDynamics(sp) == 0)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(12, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
      }

      // F that acheives \code{Ftarget}% of SSB0 in the end of the projection
      if(HCR == 3){
        // Using ssb rather than SBF because of multi-species interactions arent in SBF
        jnll_comp(12, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Ftarget_percent(sp));
      }

      // -- SPR
      if(HCR > 3){
        jnll_comp(12, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
        jnll_comp(12, sp)  += 200*square((SPRtarget(sp)/SPR0(sp))-Ftarget_percent(sp));
      }
    }

    // -- Single-species dynamic reference points
    if((DynamicHCR == 1) & (forecast(sp) == 1) & (msmMode == 0) & (estDynamics(sp) == 0)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(12, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
      }

      for(yr = 1; yr < nyrs; yr++){ // No initial abundance
        // F that acheives Ftarget% of SSB0y
        if(HCR == 3){
          jnll_comp(12, sp)  += 200*square((DynamicSBF(sp, yr)/DynamicSB0(sp, yr))-Ftarget_percent(sp));
        }
      }

      // -- SPR
      if(HCR > 3){
        jnll_comp(12, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
        jnll_comp(12, sp)  += 200*square((SPRtarget(sp)/SPR0(sp))-Ftarget_percent(sp));
      }
    }


    // -- Multi-species static reference points (all biomass_depletion based)
    if((forecast(sp) == 1) & (msmMode > 0) & (estDynamics(sp) == 0)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(12, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Flimit_percent(sp));
      }

      // F that achieves \code{Ftarget}% of SSB0 in the end of the projection
      if(HCR == 3){
        jnll_comp(12, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Ftarget_percent(sp));
      }

      // Tiered HCRs with Limit and Targets
      if(HCR == 6){
        jnll_comp(12, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Flimit_percent(sp));
      }
    }
  }


  // Slot 13 -- N-at-age < 0 penalty. See posfun
  for(sp = 0; sp < nspp; sp++){
    if(estDynamics(sp) == 0){
      jnll_comp(13, sp) += zero_N_pen(sp);
    }
  }


  // Slots 14-15 -- M_at_age prior
  // M1_model 0 = use fixed natural mortality from M1_base, 1 = estimate sex- and age-invariant M1_at_age, 2 = sex-specific (two-sex model), age-invariant M1_at_age, 3 =   estimate sex- and age-specific M1_at_age.
  for(sp = 0; sp < nspp; sp++) {
    // PRIORS
    // Prior on M1_at_age only and using species specific M1_at_age
    if( (M1_model(sp) == 1) & (M1_use_prior(sp) == 1) & (M2_use_prior(sp) == 0)){
      jnll_comp(15, sp) -= dnorm(ln_M1(sp, 0, 0), log(M_prior(sp)) + square(M_prior_sd(sp))/2, M_prior_sd(sp), true);
    }

    // Prior on M1_at_age only and using species and sex specific M1_at_age
    if( (M1_model(sp) == 2) & (M1_use_prior(sp) == 1) & (M2_use_prior(sp) == 0)){
      for(sex = 0; sex < nsex(sp); sex ++){
        jnll_comp(15, sp) -= dnorm(ln_M1(sp, sex, 0), log(M_prior(sp)) + square(M_prior_sd(sp))/2, M_prior_sd(sp), true);
      }
    }

    // Prior on M1_at_age only and using species, sex, and age specific M1_at_age
    if( (M1_model(sp) == 3) & (M1_use_prior(sp) == 1) & (M2_use_prior(sp) == 0)){
      for(sex = 0; sex < nsex(sp); sex ++){
        for(age = 0; age < nages(sp); age++) {
          jnll_comp(15, sp) -= dnorm(ln_M1(sp, sex, age), log(M_prior(sp)) + square(M_prior_sd(sp))/2, M_prior_sd(sp), true);
        }
      }
    }

    // Prior on total M_at_age (M1_at_age and M2_at_age)
    if( M2_use_prior(sp) == 1){
      for(sex = 0; sex < nsex(sp); sex ++){
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs_hind; yr++) {
            jnll_comp(15, sp) -= dnorm(log(M_at_age(sp, sex, age, yr)), log(M_prior(sp)) + square(M_prior_sd(sp))/2, M_prior_sd(sp), true);
          }
        }
      }
    }


    // RANDOM EFFECTS LIKELIHOOD
    // M1 random effects are applied to each species if M1_model = 1 or each species AND sex if M1_model = 2.
    // Variance and correlation coefficients are species-specific, but sex-invariant.
    // - M1_re = 0: No random effects (default).
    // - M1_re = 1: Random effects varies by age, but uncorrelated (IID) and constant over years.
    // - M1_re = 2: Random effects varies by year, but uncorrelated (IID) and constant over ages.
    // - M1_re = 3: Random effects varies by year and age, but uncorrelated (IID).
    // - M1_re = 4: Correlated AR1 random effects varies by age, but constant over years.
    // - M1_re = 5: Correlated AR1 random effects varies by year, but constant over ages.
    // - M1_re = 6: Correlated 2D-AR1 random effects varies by year and age.

    // M1_re = 1/4: Random effects varies by age (IID or AR1) and constant over years.
    if((M1_re(sp) == 1) | (M1_re(sp) == 4)){
      Type sigma_M = exp(M1_dev_ln_sd(sp, 0));
      Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
      Type Sigma_M = pow(pow(sigma_M, 2) / (1 - pow(rho_M_a, 2)), 0.5);


      // likelihood of M deviations
      vector<Type> M_re_age(nages(sp)); M_re_age.setZero();
      for(age = 0; age < nages(sp); age++) {
        M_re_age(age) = ln_M1_dev(sp, 0, age, 0);
      }
      jnll_comp(16, sp) += SCALE(AR1(rho_M_a), Sigma_M)(M_re_age);

      // Add males for 2-sex M1
      if(M1_model(sp) == 2){
        for(age = 0; age < nages(sp); age++) {
          M_re_age(age) = ln_M1_dev(sp, 1, age, 0);
        }
        jnll_comp(16, sp) += SCALE(AR1(rho_M_a), Sigma_M)(M_re_age);
      }
    }

    // M1_re = 2/5: Random effects varies by year (IID or AR1) and constant over ages.
    if((M1_re(sp) == 2) | (M1_re(sp) == 5)){

      Type sigma_M = exp(M1_dev_ln_sd(sp, 0));
      Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
      Type Sigma_M = pow(pow(sigma_M, 2) / (1 - pow(rho_M_y, 2)), 0.5);


      // likelihood of M deviations
      vector<Type> M_re_yr(nyrs_hind); M_re_yr.setZero();
      for(yr = 0; yr < nyrs_hind; yr++) {
        M_re_yr(yr) = ln_M1_dev(sp, 0, 0, yr);
      }
      jnll_comp(16, sp) += SCALE(AR1(rho_M_y), Sigma_M)(M_re_yr);

      // Add males for 2-sex M1
      if(M1_model(sp) == 2){
        for(yr = 0; yr < nyrs_hind; yr++) {
          M_re_yr(yr) = ln_M1_dev(sp, 1, 0, yr);
        }
        jnll_comp(16, sp) += SCALE(AR1(rho_M_y), Sigma_M)(M_re_yr);
      }
    }


    // M1_re = 3/6: Random effects varies by age and year (IID or 2D-AR1)
    if((M1_re(sp) == 3) | (M1_re(sp) == 6)){

      Type sigma_M = exp(M1_dev_ln_sd(sp, 0));
      Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
      Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
      Type Sigma_M = pow(pow(sigma_M,2) / ((1-pow(rho_M_y,2))*(1-pow(rho_M_a,2))),0.5);


      array<Type> M_re_a_yr(nages(sp), nyrs_hind);
      for(age = 0; age < nages(sp); age++) {
        for(yr = 0; yr < nyrs_hind; yr++) {
          M_re_a_yr(age, yr) = ln_M1_dev(sp, 0, age, yr);
        }
      }

      jnll_comp(16, sp) += SCALE(SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)), Sigma_M)(M_re_a_yr); // must be array, not matrix!

      // Add males for 2-sex M1
      if(M1_model(sp) == 2){
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs_hind; yr++) {
            M_re_a_yr(age, yr) = ln_M1_dev(sp, 1, age, yr);
          }
        }

        jnll_comp(16, sp) += SCALE(SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)), Sigma_M)(M_re_a_yr); // must be array, not matrix!
      }
    }
  }


  // 14.3. Diet likelihood components
  // -- MSVPA-based predation with estimated suitability AND Kinzey and Punt
  if((msmMode > 2) | (imax(suitMode) > 0)) {
    // Slot 15 -- Diet weight likelihood
    for(int stom_ind = 0; stom_ind < diet_obs.rows(); stom_ind++){

      rsp = diet_ctl(stom_ind, 0) - 1; // Index of predator

      if((msmMode > 2) | (suitMode(rsp) > 0)) { // If estimating
        jnll_comp(18, rsp) -= diet_comp_weights(rsp) * Type(diet_obs(stom_ind, 0)) * (diet_obs(stom_ind, 1) + 0.00001) * log((diet_hat(stom_ind, 1) + 0.00001)/(diet_obs(stom_ind, 1) + 0.00001));
      }
    }
  } // End diet proportion by weight component


  // 14.4. Diet likelihood components for Kinzey and Punt predation
  /*
   if(msmMode > 2){
   // Slot 15 -- Ration likelihood
   for(yr = 0; yr < nyrs_hind; yr++) {
   for(sp = 0; sp < nspp; sp++) {
   for(sex = 0; sex < nsex(sp); sex ++){
   for(age = 1; age < nages(sp); age++) { // don't include age zero in likelihood
   if(ration(sp, sex, age, yr) > 0){
   if(ration_hat(sp, sex, age, yr) > 0){
   jnll_comp(16, sp) -= dnorm(log(ration(sp, sex, age, yr))  - square(sd_ration) / 2,  log( ration_hat(sp, sex, age, yr)), sd_ration, true);
   }
   }
   }
   }
   }
   }


   // Slot 16 -- Ration penalalties FIXME: not sure if necessary: talk to Andre
   for(sp = 0; sp < nspp; sp++) {
   for(sex = 0; sex < nsex(sp); sex ++){
   for(age = 0; age < nages(sp); age++) {
   for(yr = 0; yr < nyrs_hind; yr++) {
   if(ration_hat(sp, sex, age, yr) > 0){
   //jnll_comp(17, sp) += 20 *  pow(ration_hat(sp, sex, age, yr) - ration_hat_ave(sp, sex, age), 2);
   }
   }
   }
   }
   }
   } // End if statement for Kinzey diet likelihood
   */


  // ------------------------------------------------------------------------- //
  // 12. REPORT SECTION                                                        //
  // ------------------------------------------------------------------------- //

  // 12.0 Report indices
  /*
   REPORT( flt_type );
   REPORT( flt_spp );
   REPORT( flt_sel_ind );
   REPORT( flt_sel_type );
   REPORT( flt_nselages );
   REPORT( flt_varying_sel );
   REPORT( flt_units );
   REPORT( flt_wt_index );
   REPORT( flt_age_transition_index );
   REPORT( flt_q_ind );
   REPORT( est_index_q );
   REPORT( index_varying_q );
   REPORT( est_sigma_index );
   REPORT( est_sigma_fsh );
   REPORT( diet_ctl );
   REPORT( r_sexes );
   REPORT( k_sexes );
   REPORT( joint_adjust );
   */


  // 12.1. Population components
  REPORT( pop_scalar );
  ADREPORT( pop_scalar );
  REPORT( avg_R );
  //REPORT( mature_females );
  REPORT( Z_at_age );
  REPORT( N_at_age );
  REPORT( avgN_at_age );
  //REPORT( sex_ratio_hat );
  //REPORT( avg_sex_ratio_hat );
  //REPORT( S_at_age );
  REPORT( biomass_at_age );
  //REPORT( ssb_at_age );
  REPORT( biomass_depletion );
  REPORT( ssb_depletion );
  REPORT( biomass );
  REPORT( ssb );
  REPORT( exploitable_biomass );
  //REPORT(R_sd);
  REPORT( R0 );
  REPORT( R_init );
  REPORT( R );
  REPORT( R_hat );
  REPORT( M_at_age );

  // ADREPORT( B_eaten_as_prey );
  // ADREPORT( M_at_age );
  // ADREPORT( Z_at_age );
  // ADREPORT( biomass_depletion );
  // ADREPORT( ssb_depletion );
  ADREPORT( biomass );
  ADREPORT( ssb );
  matrix<Type>  log_ssb = ssb;  log_ssb = log(ssb.array());// Fixed n-at-age scaling coefficient
  REPORT( log_ssb );
  ADREPORT( log_ssb );
  // ADREPORT( R_sd );
  ADREPORT( R );


  // -- 12.2. Biological reference points
  REPORT( NByageF );
  REPORT( NByage0);
  //REPORT( DynamicNByage0);
  //REPORT( DynamicNByageF);
  //REPORT( NbyageSPR);
  REPORT( B0 );
  REPORT( SB0 );
  REPORT( SBF );
  REPORT( DynamicB0 );
  REPORT( DynamicSB0 );
  REPORT( DynamicSBF );
  REPORT( SPR0 );
  REPORT( SPRFinit );
  REPORT( SPRlimit );
  REPORT( SPRtarget );
  REPORT( steepness );
  REPORT( proj_F );
  REPORT( Ftarget );
  REPORT( Flimit );
  //REPORT( Flimit_age_spp );
  //REPORT( Ftarget_age_spp );


  // -- 12.3. Selectivity
  REPORT( sel );
  /*
   REPORT( avg_sel );
   REPORT( emp_sel_obs );
   REPORT( sel_tmp );
   REPORT( sel_dev_sd );
   REPORT( sel_curve_pen );
   */


  // -- 12.4. Survey components
  REPORT( index_hat );
  REPORT( ln_index_sd );
  REPORT( index_q );
  vector<Type>  log_index_hat = index_hat;  log_index_hat = log(index_hat.array());// Fixed n-at-age scaling coefficient
  REPORT( log_index_hat );
  ADREPORT( log_index_hat );
  /*
   REPORT( index_q_analytical );
   REPORT( index_q_sd );
   REPORT( index_q_dev_sd );
   REPORT( index_q_dev );
   */


  // -- 12.5. Fishery components
  REPORT( F_spp );
  REPORT( F_flt );
  REPORT( F_flt_age );
  REPORT( F_spp_at_age );
  REPORT( catch_hat );
  REPORT( max_catch_hat );
  REPORT( ln_catch_sd );


  // 12.6. Age/length composition
  REPORT( age_obs_hat );
  REPORT( age_hat );
  REPORT( comp_obs );
  REPORT( comp_hat );
  REPORT( true_age_comp_hat );
  // REPORT( n_hat );
  REPORT( comp_n );


  // -- 12.7. Likelihood components
  REPORT( jnll_comp );


  // -- 12.8. Ration components
  REPORT( consumption_at_age );
  // REPORT( mnWt_obs );
  REPORT( fT );
  REPORT( ration );

  // -- 12.9. Predation components
  /*
   REPORT( suma_suit );
   REPORT( suitability );
   REPORT( suit_other );
   REPORT( stom_div_bio );
   REPORT( diet_prop );
   REPORT( avail_food );
   REPORT( other_food_diet_prop );
   REPORT( M2_prop );
   REPORT( diet_prop_hat );
   REPORT( B_eaten_as_pred );
   REPORT( N_eaten );
   */

  REPORT( suitability );
  REPORT( M1_at_age );
  REPORT( M2_at_age );
  REPORT( B_eaten );
  REPORT( B_eaten_as_prey );
  REPORT(vulnerability);
  REPORT(vulnerability_other);
  REPORT( gam_a );
  REPORT( gam_b );
  REPORT( diet_hat );


  // -- 12.10. Kinzey predation functions
  /*
   REPORT( H_1 );
   REPORT( H_1a );
   REPORT( H_1b );
   REPORT( H_2 );
   REPORT( H_3 );

   REPORT( N_pred_yrs );
   REPORT( N_prey_yrs );
   REPORT( N_pred_eq );
   REPORT( N_prey_eq );

   REPORT( pred_resp );
   REPORT( Pred_r );
   REPORT( Prey_r );

   REPORT( ration_hat );
   REPORT( ration_hat_ave );
   */


  // ------------------------------------------------------------------------- //
  // 13. END MODEL                                                             //
  // ------------------------------------------------------------------------- //

  Type jnll = 0;

  // Estimation mode
  if(estimateMode < 3) {
    jnll = jnll_comp.sum();
  }

  // Debug mode
  if(estimateMode > 2) {
    jnll = dummy * dummy;
  }


  REPORT( jnll );
  return jnll;
}

// ------------------------------------------------------------------------- //
// 14. CHANGE LOG                                                            //
// ------------------------------------------------------------------------- //
// Changes from 2017 assessment
// 1. Normalize age-transition matrix prior to use (PCod rows did not sum to 1)
// 2. Added log-normal recruitment deviation bias correction
// 3. Added month to survey index calculation. EBS specific. BT survey age/length composition estimated for month 6 rather than 0, similar to BT survey biomass
// 4. Fixed mis-specification of multinomial for fishery composition
// 5. Normalized survey and fishery selectivity so that max = 1
// 6. Fixed initialization of population
// 7. Fixed estimation routine of suitability coefficients (N used to caclulate suitability is no-longer different)
// 8. Changed ration calculation over nyrs rather than nTyrs
// 9. Allowed retrospective estimation
// 10. Added month to comp caclulation. EBS specific: Acoustic comp and biomass data are assumed to come from month 6 rather than 0, given they survey in the summer
// 11. Added time varying selectivity and catchability
// 12. Can make survey selectivity and/or catchability equal across surveys
// 13. Added cabability to have a two sex model
// 14. Added in spawning stock biomass weight
// 15. Added in multiple weights for surveys
// 16. Added in spawning month mortality adjustment
// 17. Removed constant 0.0001 added to M1_at_age
// 18. Had the model estimate diet_obs
// 19. Added flexibility to fixed n-at-age to have age, sex specific scaling paramters and estimate selectivity
//  Fixme: denominator is zero somewhere. Log of negative number. Check suitability. Make other prey a very large number.
//  Look at M2_at_age: suitability: and consumption. Make sure positive.
// 20. Added analytical q for time-varying survey sigma inputs
// 21. All random effect selectivity and catchability deviates are commented out
// 22. Estimate M1_at_age for sex/species
// 23. Fixed suitability estimation (use hindcast only)
// 24. Added in non-parametric time varying selectivity similar to Hake
// 25. Added in dynamic reference points
// 26. Added biomass_depletion and F
// 27. Added stock-recruitment relationships
// 28. Added dirichlet multinomial for age/length composition data
// 29. Removed ln_mean_F + F_dev, converted to ln_F
