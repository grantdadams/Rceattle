#include <TMB.hpp>
// ------------------------------------------------------------------------- //
//                 CEATTLE version 3.1.2                                     //
//                  Template Model Builder                                   //
//               Multispecies Statistical Model                              //
//          Bioenergetic-based Assessment for Understanding                  //
//              Biomass Linkages To The Environment                          //
//                  in the Bering Sea                                        //
// AUTHORS:   Kirstin Holsman, Jim Ianelli                                   //
//            Modified by Grant Adams                                        //
// CITATIONS:                                                                //
// 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
// ------------------------------------------------------------------------- //
//
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
// 17. Removed constant 0.0001 added to M1
// 18. Had the model estimate UobsWtAge
// 19. Added flexibility to fixed n-at-age to have age, sex specific scaling paramters and estimate selectivity
//  Fixme: denominator is zero somewhere. Log of negative number. Check suitability. Make other prey a very large number.
//  Look at M2: suitability: and consumption. Make sure positive.
// 20. Added analytical q for time-varying survey sigma inputs
// 21. All random effect selectivity and catchability deviates are commented out
// 22. Estimate M1 for sex/species
// 23. Fixed suitability estimation (use hindcast only)
// 24. Added in non-parametric time varying selectivity similar to Hake
// 25. Added in dynamic reference points
// 26. Added depletion and F
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
//  13. End
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

// Function for taking the maximum of a matrix
template<class Type>
Type max_matrix(matrix<Type> m1)
{
  Type res = m1(0,0);
  for(int i=0; i < m1.rows(); i++){
    for(int j=0; j < m1.cols(); j++){

      Type res2 = m1(i,j);
      if(res < res2){
        res = res2;
      }
    }
  }
  return res;
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


// Positive function
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}


//  Function to getting sqaure
template <class Type> Type square(Type x){return x*x;}

// Function for getting min of an IVECTOR and return an int
template <class Type>
Type imin(const vector<Type> &x)
{
  int res = x[0];
  for(int i=0; i < x.size(); i++){
    int res2 = x[i];
    if(res > res2){
      res = res2;
    }
  }
  return res;
}

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

  if (length <= 0){
    std::cerr << "Error -- vector size too small in first_difference(const dvector&)" << std::endl;
    return(0);
  }

  vector<Type> tmp(length);

  for (int i = 0; i < length; i++){
    tmp(i) = x(i+1) - x(i);
  }

  return tmp;
}

// Function for elementwise division
template <class Type>
matrix<Type> elem_div(matrix<Type> m1, matrix<Type> m2){

  int m1r = m1.rows();
  int m2r = m2.rows();

  int m1c = m1.cols();
  int m2c = m1.cols();

  if (m1r != m2r){
    std::cerr << "Error -- number of rows in matrices does not match" << std::endl;
    return(0);
  }

  if (m1c != m2c){
    std::cerr << "Error -- number of columns in matrices does not match" << std::endl;
    return(0);
  }

  matrix<Type> m3(m1r, m1c);

  // Elementwise division
  for(int r = 0; r < m1r; r++){
    for(int c = 0; c < m1c; c++){
      m3(r, c) = m1(r, c) / m2(r, c);
    }
  }
  return m3;
}

// Function for elementwise matrix exponential functions
template <class Type>
matrix<Type> elem_pow(matrix<Type> m1, Type exponent){

  int nrow = m1.rows();
  int ncol = m1.cols();

  matrix<Type> m2(nrow, ncol);

  // Elementwise division
  for(int r = 0; r < nrow; r++){
    for(int c = 0; c < ncol; c++){
      m2(r, c) = pow( m1(r, c), exponent);
    }
  }
  return m2;
}

// Function for to extract a layer from an array
template<class Type>
matrix<Type> matrix_from_array_d3(array<Type> a1, int sheet){
  vector<int> a1_dim = a1.dim;
  matrix<Type> m1(a1_dim(0), a1_dim(1));

  // If array is array
  if(a1_dim.size() == 3){
    for(int row = 0; row < a1_dim(0); row++){
      for(int col = 0; col < a1_dim(1); col++){
        m1(row, col) = a1(row, col, sheet);
      }
    }
  }

  // If array is matrix
  if(a1_dim.size() == 2){
    for(int row = 0; row < a1_dim(0); row++){
      for(int col = 0; col < a1_dim(1); col++){
        m1(row, col) = a1(row, col);
      }
    }
  }
  return m1;
}

// Function for to extract a layer from an array
template<class Type>
matrix<Type> matrix_from_array_d1(array<Type> a1, int sheet){
  vector<int> a1_dim = a1.dim;
  matrix<Type> m1(a1_dim(1), a1_dim(2));

  // If array is array
  if(a1_dim.size() == 3){
    for(int row = 0; row < a1_dim(1); row++){
      for(int col = 0; col < a1_dim(2); col++){
        m1(row, col) = a1(sheet, row, col);
      }
    }
  }

  // If array is matrix
  if(a1_dim.size() == 2){
    for(int row = 0; row < a1_dim(0); row++){
      for(int col = 0; col < a1_dim(1); col++){
        m1(row, col) = a1(row, col);
      }
    }
  }
  return m1;
}

// Function for the multiplaction of a vector and matrix: in r "v %*% m"
template <class Type>
vector<Type> vec_mat_prod(vector<Type> v1, matrix<Type> m1){

  if (v1.size() != m1.rows()){
    std::cerr << "Error -- length of vector does not equal the number of rows in the matrix" << std::endl;
    return(0);
  }

  vector<Type> v2(m1.rows());
  vector<Type> v3(m1.cols());

  // Vector * matrix
  for(int c = 0; c < m1.cols(); c++){
    v2 = m1.col(c);
    v3(c) = (v1 * v2).sum();
  }
  return v3;
}


// function to get row vector from array
template<class Type>
vector<Type> col_from_3D_array(array<Type> a1, int col, int sheet){
  vector<int> a1_dim = a1.dim;
  vector<Type> v1(a1_dim(0));

  // If array is array
  if(a1_dim.size() == 3){
    for(int row = 0; row < a1_dim(0); row++){
      v1(row) = a1(row, col, sheet);
    }
  }

  // If array is matrix
  if(a1_dim.size() == 2){
    for(int row = 0; row < a1_dim(0); row++){
      v1(row) = a1(row, col);
    }
  }
  return v1;
}


// Function to trim matrix
template<class Type>
matrix<Type> trim_matrix(matrix<Type> m1, int nrow, int ncol){
  // dim is a DATA_IMATRIX with structure nrow, ncol, ID # of object we wish to interogate
  matrix<Type> m2(nrow, ncol);

  // If array is matrix
  for(int row = 0; row < nrow; row++){
    for(int col = 0; col < ncol; col++){
      m2(row, col) = m1(row, col);
    }
  }
  return m2;
}

// Function to trim vector
template<class Type>
matrix<Type> trim_vector(vector<Type> v1, int length){
  // dim is a DATA_IMATRIX with structure nrow, ncol, ID # of object we wish to interogate
  vector<Type> v2(length);

  // If array is matrix
  for(int i = 0; i < length; i++){
    v2(i) = v1(i);
  }
  return v2;
}


template<class Type>
Type objective_function<Type>::operator() () {
  // ------------------------------------------------------------------------- //
  // 1. MODEL CONFIGURATION                                                    //
  // ------------------------------------------------------------------------- //
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(estimateMode);             // Logical to debug or not
  DATA_SCALAR( minNByage );               // Hard constraint that the lowest a numbers at age can be is 1
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
  DATA_INTEGER(suitMode);                 // Estimate suitability
  DATA_INTEGER(avgnMode);                 // N used for predation function
  //    0 = AvgN
  //    1 = N*exp(-Z / 2))
  //    2 = N
  DATA_INTEGER( niter );                  // Number of loops for MSM mode
  DATA_INTEGER(initMode);                 // How the age-structure is initialized


  // Columns =  Species, Survey, Selectivity, Estimate_q, log_q_prior
  // Selectivity 0 = empirical, 1 = logistic, 2 = non-parametric, 3 = double logistic

  // 1.2. Temporal dimensions
  DATA_INTEGER( styr );                   // Start year
  DATA_INTEGER( endyr );                  // End of estimation years
  DATA_INTEGER( meanyr );                 // The last year used to calculate population averages.
  DATA_INTEGER( projyr );                 // End year of projection
  int nyrs = projyr - styr + 1;
  int nyrs_mean = meanyr - styr + 1;
  int nyrs_hind = endyr - styr + 1;

  // 1.3. Number of species
  DATA_INTEGER( nspp );                   // Number of species (prey)
  DATA_IVECTOR( pop_wt_index );           // Dim 3 of wt to use for population dynamics
  DATA_IVECTOR( ssb_wt_index );           // Dim 3 of wt to use for spawning stock biomass calculation
  DATA_IVECTOR( pop_age_transition_index );// Dim 3 of wt to use for age_transition_matrix
  DATA_IVECTOR( estDynamics );            // Index indicating wether the population parameters are estimated (0), numbers-at-age are provided (1), or an index of numbers-at-age multiplied by an estimated scalar is used (2)
  DATA_IVECTOR( est_sex_ratio );          // Is sex ration F/(M+F) to be included in the likelihood; 0 = no, 1 = use annual average across ages (uses 2nd age in sex_ratio data), 2 = age, and year specific (TBD)
  pop_wt_index -= 1;                      // Indexing starts at 0
  ssb_wt_index -= 1;                      // Indexing starts at 0
  pop_age_transition_index -= 1;          // Indexing starts at 0

  // 1.4. RECRUITMENT SETTINGS
  DATA_INTEGER(srr_fun);                  // Stock recruit relationship
  DATA_INTEGER(proj_mean_rec);
  //    0 = project recruitment using ln_R0 and rec devs
  //    1 = project recruitment using mean rec (can also have adjusted rec devs)
  DATA_INTEGER(srr_use_prior);            // Logical of wether to add normal prior to stock recruit-relationship
  DATA_VECTOR(srr_prior_mean);           // Prior mean for stock recruit relationship parameter
  DATA_VECTOR(srr_prior_sd);             // Prior sd for stock recruit relationship parameter

  // 1.5. HCR and projection settings
  DATA_INTEGER(forecast);                 // Switch, wether or not we are estimating the hindcast or forecast
  DATA_INTEGER(HCR);                      // Function to be used for harvest control rule
  DATA_INTEGER(DynamicHCR);               // TRUE/FALSE. Wether to use static or dynamic reference points (default = FALSE)
  DATA_VECTOR(Ptarget);                   // Target spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0
  DATA_VECTOR(Plimit);                    // Limit spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0
  DATA_VECTOR(FsprTarget);                // Percentage of spawning-stock biomass per recruit at F = 0 used to find the target F-spr
  DATA_VECTOR(FsprLimit);                 // Percentage of spawning-stock biomass per recruit at F = 0 used to find the limit F-spr
  DATA_VECTOR(Alpha);                     // Parameter used in NPFMC Tier 3 HCR
  DATA_VECTOR(Fmult);                     // Multiplier for target fishing mortality (Fmult * Ftarget). For example 0.75 * F40%.
  DATA_VECTOR(QnormHCR);                  // Add to Flimit to set Ftarget based on the Pstar approach: Flimit + qnorm(Pstar, 0, Sigma)


  // 1.5. MODEL OBJECTS
  // 1.5.1. LOOPING INDICES -- k = observation, sp = species, sex = sex (0 = combines; 1 = females; 2 = males), age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length
  int sp, sex, age, ln, ksp, k_sex, k_age, k_ln, yr, rsp, r_sex, r_age, r_ln; // k
  int srv, flt;                                                         // Survey and fishery indices
  int flt_yr, flt_sex, comp_type;
  int flt_ind, fsh_ind, srv_ind, comp_ind, yr_ind;                      // Indices for survey sets
  int sel_type, sel_varying, nselages;
  int ncnt;    // Pointers
  Type mo = 0;                                                          // Month float
  if (msmMode == 0) { niter = 1; }                                      // Number of iterations for SS mode
  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                           //
  // ------------------------------------------------------------------------- //

  // 2.1. FIXED VALUES
  Type sd_ration = 0.05;                  // SD of ration likelihood


  // 2.2. DIMENSIONS

  // -- 2.2.2. Species attributes
  DATA_IVECTOR( nsex );                   // Number of sexes to be modelled; n = [1, nspp]; 1 = sexes combined/single sex, 2 = 2 sexes
  DATA_VECTOR( spawn_month );             // Month of spawning to adjust mortality; n = [1, nspp]
  DATA_VECTOR( R_sexr );                  // Sex ratio of recruitment; n = [1, nspp]; R_sexr = R to females/ R to females and males
  DATA_IVECTOR( nages );                  // Number of species (prey) ages; n = [1, nspp]
  DATA_IVECTOR( minage );                 // Minimum age of each species; n = [1, nspp]
  DATA_IVECTOR( nlengths );               // Number of species (prey) lengths; n = [1, nspp]
  DATA_ARRAY( NByageFixed );              // Provided estimates of numbers- or index-at-age to be multiplied (or not) by pop_scalar to get Nbyage
  Type stom_tau = 20;                     // Stomach sample size: FIXME - have as input
  int max_age = imax(nages);              // Integer of maximum nages to make the arrays; n = [1]
  DATA_VECTOR( MSSB0 );                   // SB0 from projecting the model forward in multi-species mode under no fishing

  // 2.3. Data controls (i.e. how to assign data to objects)
  DATA_IMATRIX( fleet_control );          // Fleet specifications

  // -- 2.3.1 Fishery Components
  DATA_IMATRIX( fsh_biom_ctl );           // Info for fishery biomass; columns = Fishery_name, Fishery_code, Species, Year
  DATA_IMATRIX( fsh_biom_n );             // Info for fishery biomass; columns = Month
  DATA_MATRIX( fsh_biom_obs );            // Observed fishery catch biomass (kg) and log_sd; n = [nobs_fish_biom, 2]; columns = Observation, Error

  // -- 2.3.2 Survey components
  DATA_IMATRIX( srv_biom_ctl );           // Info for survey biomass; columns = Survey_name, Survey_code, Species, Year
  DATA_MATRIX( srv_biom_n );              // Info for survey biomass; columns = Month
  DATA_MATRIX( srv_biom_obs );            // Observed survey biomass (kg) and log_sd; n = [nobs_srv_biom, 2]; columns = Observation, Error
  DATA_VECTOR( ln_srv_q_prior );          // Prior mean for survey catchability; n = [nflt]

  // -- 2.3.3. Composition data
  DATA_IMATRIX( comp_ctl );               // Info on observed survey age/length comp; columns = Survey_name, Survey_code, Species, Year
  DATA_MATRIX( comp_n );                  // Month and sample size on observed survey age/length comp; columns = Month, Sample size
  DATA_MATRIX( comp_obs );                // Observed survey age/length comp; cols = Comp_1, Comp_2, etc. can be proportion

  // -- 2.3.4 Age and selectivity
  DATA_IMATRIX( emp_sel_ctl );            // Info on empirical fishery selectivity; columns =  Fishery_name, Fishery_code, Species, Year
  DATA_MATRIX( emp_sel_obs );             // Observed emprical fishery selectivity; columns = Compe_1, Comp_2, etc.
  DATA_ARRAY( age_trans_matrix);          // observed sp_age/size compositions; n = [nspp, nages, srv_age_bins]
  DATA_ARRAY( age_error );                // Array of aging error matrices for each species; n = [nspp, nages, nages]

  // -- 2.3.5. Weight-at-age
  DATA_ARRAY( wt );                       // Weight-at-age by year; n = [nweight, sex, nages, nyrs]: FIXME: Change nyrs to nyrs_wt_at_age if data don't span entire bit

  // 2.3.6. Diet data
  DATA_VECTOR( fday );                    // number of foraging days for each predator; n = [1, nspp]
  DATA_ARRAY( Pyrs );                     // Relative-foraging rate;  n = [nspp, nyrs+1, nages]: #FIXME - Assuming this is the same as Pby_yr?
  DATA_MATRIX( UobsWtAge );               // pred, prey, predA, preyA U observations (mean wt_hat of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_IMATRIX( UobsWtAge_ctl );          // Info on pred, prey, predA, preyA U matrix (mean wt_hat of prey in each pred age); n = [nspp, nspp, max_age, max_age]

  // 2.3.7. Environmental data
  DATA_IVECTOR( env_yrs );                // Years of hindcast data; n = [1, nTyrs] #FIXME - changed the name of this in retro_data2017_asssmnt.dat
  DATA_MATRIX( env_index );               // Matrix o environmental predictors such as bottom temperature; n = [1,  nTyrs ]
  int nTyrs = env_yrs.size();             // Number of temperature years; n = [1] #FIXME - changed the name of this in retro_data2017_asssmnt.dat

  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( other_food );              // Biomass of other prey (kg); n = [1, nspp]
  DATA_VECTOR( Pvalue );                  // This scales the pvalue used, proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge; n = [1, nspp]
  DATA_IVECTOR( Ceq );                    // Ceq: which Comsumption equation to use; n = [1, nspp]; Currently all sp = 1
  DATA_IVECTOR( Cindex );                 // Cindex, which environmental index in env_index should drive bioenergetics.
  DATA_VECTOR( CA );                      // Wt specific intercept of Cmax=CA*W^CB; n = [1, nspp]
  DATA_VECTOR( CB );                      // Wt specific slope of Cmax=CA*W^CB; n = [1, nspp]
  DATA_VECTOR( Qc );                      // used in fT, QC value; n = [1, nspp]
  DATA_VECTOR( Tco );                     // used in fT, thermal optimum; n = [1, nspp]
  DATA_VECTOR( Tcm );                     // used in fT, thermal max; n = [1, nspp]
  DATA_VECTOR( Tcl );                     // used in fT eq 3, limit; n = [1, nspp]
  DATA_VECTOR( CK1 );                     // used in fT eq 3, limit where C is .98 max (ascending); n = [1, nspp]
  DATA_VECTOR( CK4 );                     // used in fT eq 3, temp where C is .98 max (descending); n = [1, nspp]

  // -- 2.4.2. NOT USED von Bertalannfy growth function (VBGF): This is used to calculate future weight-at-age: NOT YET IMPLEMENTED

  // -- 2.4.3. Weight-at-length parameters
  DATA_MATRIX( aLW );                     // LW a&b regression coefs for W=a*L^b; n = [2, nspp]

  // -- 2.4.4. Others
  DATA_MATRIX( sex_ratio );               // Proportion-at-age of females of population; n = [nspp, nages]
  DATA_MATRIX( pmature );                 // Proportion of mature females at age; [nspp, nages]


  // ------------------------------------------------------------------------- //
  // 3. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //

  PARAMETER( dummy );                             // Variable to test derived quantities given input parameters; n = [1]
  PARAMETER_MATRIX( ln_pop_scalar );              // Scalar to multiply supplied numbers at age by; n = [nspp, nages]

  // -- 3.1. Recruitment parameters
  PARAMETER_MATRIX( rec_pars );                   // Recruitment parameters; n = [1, nspp]
  PARAMETER_VECTOR( ln_rec_sigma );               // Standard deviation of recruitment deviations; n = [1, nspp]
  PARAMETER_MATRIX( rec_dev );                    // Annual recruitment deviation; n = [nspp, nyrs]
  PARAMETER_MATRIX( init_dev );                   // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this

  // -- 3.2. Population parameters
  PARAMETER_ARRAY( ln_M1 );                       // Natural mortality (residual if multispecies mode or total if single species mode); n = [nspp, nsex, nages]
  PARAMETER_VECTOR( ln_sex_ratio_sigma );         // Variance for sex ratio to be used; n = [nspp]


  // -- 3.3. fishing mortality parameters
  PARAMETER_VECTOR( ln_mean_F );                  // Log mean fishing mortality; n = [1, n_fsh]
  PARAMETER_MATRIX( ln_Flimit );                  // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_MATRIX( ln_Ftarget );                 // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_VECTOR( proj_F_prop );                // Proportion of fishing mortality from each fleet for projections; n = [n_fsh]
  PARAMETER_MATRIX( F_dev );                      // Annual fishing mortality deviations; n = [n_fsh, nyrs] # NOTE: The size of this will likely change

  // -- 3.4. Survey catchability parameters
  PARAMETER_VECTOR( ln_srv_q );                   // Survey catchability; n = [n_srv]
  // PARAMETER_VECTOR( srv_q_pow );                  // Survey catchability power coefficient q * B ^ q_pow or beta ln(q_y) = q_mut + beta * index_y; n = [n_srv]
  PARAMETER_MATRIX( ln_srv_q_dev );               // Annual survey catchability deviates; n = [n_srv, nyrs_hind]
  PARAMETER_VECTOR( ln_sigma_srv_q );             // Log standard deviation of prior on survey catchability; n = [1, n_srv]
  PARAMETER_VECTOR( ln_sigma_time_varying_srv_q );// Log standard deviation of time varying survey catchability; n = [1, n_srv]

  // -- 3.5. Selectivity parameters
  PARAMETER_ARRAY( sel_coff );                    // selectivity parameters for non-parametric; n = [n_selectivities, nsex, nselages]
  PARAMETER_ARRAY( sel_coff_dev );                // Annual deviates for non-parametric selectivity parameters; n = [n_selectivities, nsex, nselages]
  PARAMETER_ARRAY( ln_sel_slp );                  // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( sel_inf );                     // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( ln_sel_slp_dev );              // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( sel_inf_dev );                 // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_VECTOR( ln_sigma_sel );               // Log standard deviation of selectivity; n = [1, n_selectivities]
  PARAMETER_MATRIX( sel_curve_pen );              // Selectivity penalty for non-parametric selectivity, 2nd column is for monotonic bit

  // -- 3.6. Variance of survey and fishery time series
  PARAMETER_VECTOR( ln_sigma_srv_index );         // Log standard deviation of survey index time-series; n = [1, n_srv]
  PARAMETER_VECTOR( ln_sigma_fsh_catch );         // Log standard deviation of fishery catch time-series; n = [1, n_fsh]
  PARAMETER_VECTOR( comp_weights );               // MacCallister-Ianelli weights for fisheries data

  // -- 3.7. Kinzery predation function parameters
  PARAMETER_MATRIX(logH_1);                       // Predation functional form; n = [nspp, nspp2];
  PARAMETER_VECTOR(logH_1a);                      // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix
  PARAMETER_VECTOR(logH_1b);                      // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix
  PARAMETER_MATRIX(logH_2);                       // Predation functional form; n = [nspp, nspp]
  PARAMETER_MATRIX(logH_3);                       // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
  PARAMETER_MATRIX(H_4);                          // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;

  // 3.8. Gamma selectivity parameters
  PARAMETER_VECTOR( log_gam_a );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  PARAMETER_VECTOR( log_gam_b );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10

  // 3.9. Preference
  PARAMETER_MATRIX( log_phi );                    // Species preference coefficient; n = [nspp, nspp]


  // ------------------------------------------------------------------------- //
  // 4. DERIVED QUANTITIES SECTION                                             //
  // ------------------------------------------------------------------------- //

  // 4.1. Derived indices
  int max_bin = imax( nlengths );                                                   // Integer of maximum number of length/age bins.
  int n_flt = fleet_control.rows();
  vector<int> joint_adjust(comp_obs.rows()); joint_adjust.setZero();

  // -- 4.2. Estimated population quantities
  matrix<Type>  pop_scalar = ln_pop_scalar;  pop_scalar = exp(ln_pop_scalar.array());// Fixed n-at-age scaling coefficient; n = [nspp, nages]
  vector<Type>  mean_rec(nspp); mean_rec.setZero();                                     // Mean recruitment of hindcast; n = [nspp]
  vector<Type>  R0(nspp); R0.setZero();                                             // Expected recruitment at F = 0.
  vector<Type>  Rinit(nspp); Rinit.setZero();                                       // Expected recruitment at F = Finit (non-equilibrium).
  vector<Type>  Steepness(nspp); Steepness.setZero();                               // Expected % of R0 at 20% SSB0.
  array<Type>   biomassByage(nspp, 2, max_age, nyrs); biomassByage.setZero();          // Estimated biomass-at-age (kg)
  matrix<Type>  biomass(nspp, nyrs); biomass.setZero();                             // Estimated biomass (kg)
  matrix<Type>  biomassSSB(nspp, nyrs); biomassSSB.setZero();                       // Estimated spawning stock biomass (kg)
  matrix<Type>  depletion(nspp, nyrs); depletion.setZero();                         // Estimated biomass depletion
  matrix<Type>  depletionSSB(nspp, nyrs); depletionSSB.setZero();                   // Estimated depletion of spawning stock biomass
  array<Type>   biomassSSBByage(nspp, max_age, nyrs); biomassSSBByage.setZero();    // Spawning biomass at age (kg)
  array<Type>   M(nspp, 2, max_age, nyrs); M.setZero();                             // Total natural mortality at age
  array<Type>   M1(ln_M1.dim); M1 = exp(ln_M1);                                     // Residual or total natural mortality at age
  array<Type>   NByage(nspp, 2, max_age, nyrs); NByage.setZero();                   // Numbers at age
  array<Type>   AvgN(nspp, 2, max_age, nyrs); AvgN.setZero();                       // Average numbers-at-age
  array<Type>   sex_ratio_hat(nspp, max_age, nyrs); sex_ratio_hat.setZero();        // Estimated age-specific sex ratin
  matrix<Type>  sex_ratio_mean_hat(nspp, nyrs); sex_ratio_mean_hat.setZero();       // Estimated sex ratio across all ages
  matrix<Type>  R(nspp, nyrs); R.setZero();                                         // Estimated recruitment (n)
  array<Type>   S(nspp, 2, max_age, nyrs); S.setZero();                             // Survival at age
  array<Type>   Zed(nspp, 2, max_age, nyrs); Zed.setZero();                         // Total mortality at age
  vector<Type>  r_sigma(nspp); r_sigma.setZero();                                   // Standard deviation of recruitment variation
  vector<Type>  zero_pop_pen(nspp); zero_pop_pen.setZero();                         // Additional penalty to add to likelihood if n-at-age goes < 0
  vector<Type>  sex_ratio_sigma(nspp); sex_ratio_sigma.setZero();                   // Variance of sex ratio

  // -- 4.3. Selectivity parameters
  matrix<Type>  avgsel(n_flt, 2); avgsel.setZero();                                 // Average selectivity
  array<Type>   sel(n_flt, 2, max_age, nyrs_hind); sel.setZero();                   // Estimated selectivity at age
  array<Type>   sel_tmp(n_flt, 2, max_age); sel_tmp.setZero();                      // Temporary saved selectivity at age for estimated bits
  vector<Type>  sigma_sel(n_flt); sigma_sel.setZero();                              // Standard deviation of selectivity deviates

  // -- 4.4. Fishery components
  vector<Type>  sigma_fsh_catch(n_flt); sigma_fsh_catch.setZero();                  // Standard deviation of fishery time-series
  matrix<Type>  F_spp(nspp, nyrs); F_spp.setZero();                                 // Fully selected fishing mortality by species
  matrix<Type>  F_flt(n_flt, nyrs); F_flt.setZero();                                // Fully selected fishing mortality by fleet
  array<Type>   F_flt_age(n_flt, 2, max_age, nyrs); F_flt_age.setZero();            // Estimated fishing mortality for each fishery
  array<Type>   FlimitSPR(nspp, 2, max_age); FlimitSPR.setZero();                   // Estimated fishing mortality for each species that leads to SPRlimit
  array<Type>   FtargetSPR(nspp, 2, max_age); FtargetSPR.setZero();                 // Estimated fishing mortality for each species that leads to SPRtarget
  array<Type>   F_spp_age(nspp, 2, max_age, nyrs+1); F_spp_age.setZero();           // Sum of annual estimated fishing mortalities for each species-at-age
  vector<Type>  fsh_bio_hat(fsh_biom_obs.rows()); fsh_bio_hat.setZero();            // Estimated fishery yield (kg)
  vector<Type>  fsh_log_sd_hat(fsh_biom_obs.rows()); fsh_log_sd_hat.setZero();      // Estimated/fixed fishery log_sd (kg)

  // -- 4.5. Biological reference points
  array<Type>   NByage0(nspp, 2, max_age, nyrs); NByage0.setZero();                 // Numbers at age at mean recruitment and F = 0
  array<Type>   DynamicNByage0(nspp, 2, max_age, nyrs); DynamicNByage0.setZero();   // Numbers at age at F = 0 (accounts for annual recruitment)
  array<Type>   DynamicNByageF(nyrs, nspp, max_age, nyrs); DynamicNByageF.setZero(); // Female numbers at age at F = Ftarget (accounts for annual recruitment)
  matrix<Type>  DynamicSB0(nspp, nyrs); DynamicSB0.setZero();                       // Estimated dynamic spawning biomass at F = 0 (accounts for S-R curve)
  matrix<Type>  DynamicB0(nspp, nyrs); DynamicB0.setZero();                         // Estimated dynamic  biomass at F = 0 (accounts for S-R curve)
  matrix<Type>  DynamicSBF(nspp, nyrs); DynamicSBF.setZero();                       // Estimated dynamic spawning biomass at F = Ftarget (accounts for S-R curve)
  array<Type>   NbyageSPR(3, nspp, max_age);                                        // Estimated numbers at age for spawning biomass per recruit reference points
  vector<Type>  SPRlimit(nspp); SPRlimit.setZero();                                 // Estimated Plimit SPR
  vector<Type>  SPRtarget(nspp); SPRtarget.setZero();                               // Estimated Ptarget SPR
  vector<Type>  SPR0(nspp); SPR0.setZero();                                         // Estimated spawning biomass per recruit at F = 0
  vector<Type>  SB0(nspp); SB0.setZero();                                           // Estimated spawning stock biomass at F = 0 (Accounts for S-R)
  vector<Type>  B0(nspp); B0.setZero();                                             // Estimated biomass at F = 0 (Accounts for S-R)
  matrix<Type>  Flimit = exp(ln_Flimit.array());                                    // Target F parameter on natural scale
  matrix<Type>  Ftarget = exp(ln_Ftarget.array());                                  // Limit F parameter on natural scale
  matrix<Type>  proj_F(nspp, nyrs); proj_F.setZero();                               // Projected F using harvest control rule


  // -- 4.6. Survey components
  vector<Type>  sigma_srv_index(n_flt); sigma_srv_index.setZero();                  // Vector of standard deviation of survey index
  vector<Type>  sigma_srv_q(n_flt); sigma_srv_q.setZero();                          // Vector of standard deviation of survey catchability prior
  vector<Type>  time_varying_sigma_srv_q(n_flt); time_varying_sigma_srv_q.setZero();// Vector of standard deviation of time-varying survey catchability deviation
  Type avgsel_tmp = 0;                                                              // Temporary object for average selectivity across all ages
  vector<Type>  srv_bio_hat(srv_biom_obs.rows()); srv_bio_hat.setZero();            // Estimated survey biomass (kg)
  vector<Type>  srv_log_sd_hat(srv_biom_obs.rows()); srv_log_sd_hat.setZero();      // Estimated/fixed survey log_sd (kg)
  vector<Type>  sigma_srv_analytical(n_flt); sigma_srv_analytical.setZero();        // Temporary vector to save analytical sigma follow Ludwig and Walters 1994
  vector<Type>  srv_q_analytical(n_flt); srv_q_analytical.setZero();                // Temporary vector to save analytical sigma follow Ludwig and Walters 1994
  matrix<Type>  srv_q(n_flt, nyrs_hind); srv_q.setZero();                           // Estimated survey catchability
  vector<Type>  srv_n_obs(n_flt); srv_n_obs.setZero();                              // Vector to save the number of observations for each survey time series

  // -- 4.7. Composition data - FIXME: will blow up in nlengths is less than nages
  vector<Type>  n_hat(comp_obs.rows()) ; n_hat.setZero() ;                          // Estimated catch (numbers)
  matrix<Type>  age_hat = comp_obs; age_hat.setZero();                              // Estimated catch at true age
  matrix<Type>  true_age_comp_hat = comp_obs; true_age_comp_hat.setZero();          // True estimated age composition
  matrix<Type>  age_obs_hat = comp_obs; age_obs_hat.setZero();                      // Estimated catch at observed age (accounts for ageing error)
  matrix<Type>  comp_hat = comp_obs; comp_hat.setZero();                            // Estimated comp

  // -- 4.8. Ration components
  array<Type>   ConsumAge( nspp, 2, max_age, nyrs ); ConsumAge.setZero();           // Pre-allocated indiviudal consumption in grams per predator-age
  matrix<Type>  fT( nspp, nyrs ); fT.setZero();                                     // Pre-allocation of temperature function of consumption
  array<Type>   LbyAge( nspp, 2, max_age, nyrs ); LbyAge.setZero();                 // Length by age from LW regression
  matrix<Type>  mnWt_obs( nspp, max_age ); mnWt_obs.setZero();                      // Mean observed weight at age (across years)
  array<Type>   ration( nspp, 2, max_age, nyrs ); ration.setZero();                 // Annual ration at age (kg/yr)
  matrix<Type>  env_index_hat( nyrs, env_index.cols() ); env_index_hat.setZero();   // Environmental indices like bottom temperature

  // -- 4.9. Diet components
  array<Type>   diet_prop_weight(nspp * 2, nspp * 2, max_age, max_age, nyrs); diet_prop_weight.setZero();           // Stomach proportion by weight U
  array<Type>   diet_prop_weight_hat(nspp * 2, nspp * 2, max_age, max_age, nyrs); diet_prop_weight_hat.setZero();   // Predicted stomach proportion by weight U
  array<Type>   diet_prop_weight_ave_hat(nspp * 2, nspp * 2, max_age, max_age); diet_prop_weight_ave_hat.setZero();// Average predicted stomach proportion by weight across years U
  array<Type>   other_food_diet_prop_weight(nspp, 2, max_age, nyrs); other_food_diet_prop_weight.setZero();         // Other food diet proportion by weight
  matrix<Type>  UobsWtAge_hat = UobsWtAge; UobsWtAge_hat.setZero();                                                 // Estimated stomach proportion by weight U

  // -- 4.10. Suitability components
  array<Type>   avail_food(nspp, 2, max_age, nyrs); avail_food.setZero();                                           // Available food to predator
  array<Type>   stom_div_bio2(nspp * 2, nspp * 2, max_age, max_age, nyrs); stom_div_bio2.setZero();                 // Stomach proportion over biomass; U/ (W * N)
  array<Type>   suit_tmp(nspp * 2, nspp * 2, max_age, max_age, nyrs); suit_tmp.setZero();                           // Temporary suitability storage U
  array<Type>   suit_main(nspp * 2, nspp * 2, max_age, max_age, nyrs); suit_main.setZero();                         // Suitability/gamma selectivity of predator age u on prey age a
  array<Type>   suit_other(nspp, 2, max_age, nyrs); suit_other.setZero();                                           // Suitability not accounted for by the included prey
  array<Type>   suma_suit(nspp, 2, max_age, nyrs); suma_suit.setZero();                                             // Sum of suitabilities

  // -- 4.11. Suitability parameters
  vector<Type> gam_a = exp(log_gam_a);                                    // Predator size-selectivity: shape parameter for gamma suitability, mean for normal of logs
  vector<Type> gam_b = exp(log_gam_b);                                    // Predator size-selectivity: scale parameter for gamma suitability, sd for normal of logs
  vector<Type> sum_phi(nspp); sum_phi.setZero();                          // Sum of predator-prey preference coefficients for multinomial transformation
  matrix<Type> vulnerability(nspp, nspp); vulnerability.setZero();        // Predator-prey preference coefficients
  vector<Type> vulnerability_other(nspp); vulnerability_other.setZero();  // Preference for other food

  // -- 4.12. Predation components
  array<Type>   M2(nspp, 2, max_age, nyrs); M2.setZero();                                   // Predation mortality at age
  array<Type>   M2_prop(nspp * 2, nspp * 2, max_age, max_age, nyrs); M2_prop.setZero();     // Relative predation mortality at age from each species at age
  array<Type>   B_eaten_as_prey(nspp, 2, max_age, nyrs); B_eaten_as_prey.setZero();         // Biomass eaten as prey via predation
  array<Type>   B_eaten_as_pred(nspp, 2, max_age, nyrs); B_eaten_as_pred.setZero();         // Biomass eaten as predator via predation
  array<Type>   B_eaten(nspp * 2, nspp * 2, max_age, max_age, nyrs); B_eaten.setZero();     // Biomass of prey eaten via predation by a predator at age
  array<Type>   N_eaten(nspp * 2, nspp * 2, max_age, max_age, nyrs); N_eaten.setZero();     // Number of prey of age a eaten by predator age u

  // -- 4.13. Kinzey Functional response parameters
  matrix<Type> H_1(nspp, nspp + 1); H_1 = exp(logH_1.array());
  vector<Type> H_1a(nspp); H_1a = exp(logH_1a);
  vector<Type> H_1b(nspp); H_1b = exp(logH_1b);
  matrix<Type> H_2(nspp, nspp); H_2 = exp(logH_2.array());
  matrix<Type> H_3(nspp, nspp); H_3 = exp(logH_3.array());

  array<Type>  N_pred_yrs(nspp, 2, max_age, nyrs); N_pred_yrs.setZero();                // Effective numbers of predators for each age of prey FIXME: should be AvgN?
  array<Type>  N_prey_yrs(nspp, 2, max_age, nyrs); N_prey_yrs.setZero();                // Effective numbers of prey for each age of predator
  array<Type>  N_pred_eq(nspp, 2, max_age); N_pred_eq.setZero();                        // Effective numbers of predators for each age of prey (styr_pred)
  array<Type>  N_prey_eq(nspp, 2, max_age); N_prey_eq.setZero();                        // Effective numbers of prey for each age of predator

  array<Type>  pred_resp(nspp * 2, (nspp * 2)+1, max_age, max_age, nyrs); pred_resp.setZero();// Predator functional response +1 for other species
  array<Type>  Pred_r(nspp, 2, max_age, nyrs); Pred_r.setZero();                            // save Pred_ratio values
  array<Type>  Prey_r(nspp, 2, max_age, nyrs); Prey_r.setZero();                            // save Prey_ratio values

  array<Type> ration_hat(nspp, 2, max_age, nyrs); ration_hat.setZero();                   // Annual ration by predator age each year
  array<Type> ration_hat_ave(nspp, 2, max_age); ration_hat_ave.setZero();                 // Annual ration by predator age averaged over years

  // ------------------------------------------------------------------------- //
  // 5. INITIAL CALCULATIONS                                                   //
  // ------------------------------------------------------------------------- //

  // 5.5. Maturity and sex ratio
  for ( sp = 0; sp < nspp ; sp++) {
    for ( age = 0 ; age < nages(sp); age++ ) {
      if(nsex(sp) == 1){
        pmature( sp, age ) = pmature( sp, age ) * sex_ratio(sp, age); // Mulitply sex_ratio and pmature for 1 sex models
      }
    }
  }

  // 5.6. Calculate temperature to use
  for(int i = 0; i < env_index.cols(); i++){
    Type env_mu = env_index.col(i).sum() / nTyrs; // Fill with average bottom temperature
    for(yr = 0; yr < nyrs; yr++){
      env_index_hat(yr, i) = env_mu; // Fill in missing years with average
    }
  }

  for(int i = 0; i < env_index.cols(); i++){
    yr_ind = 0;
    for (yr = 0; yr < nTyrs; yr++) {
      yr_ind = env_yrs( yr ) - styr;
      if ((yr_ind >= 0) & (yr_ind < nyrs)) {
        env_index_hat(yr_ind, i) = env_index( yr , i );
      }
    }
  }


  // 5.7. Calculate length-at-age
  for (sp = 0; sp < nspp; sp++) {
    for (age = 0; age < nages(sp); age++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for (yr = 0; yr < nyrs; yr++) {
          // Hindcast
          if(yr < nyrs_hind){
            LbyAge( sp, sex, age, yr) = ( pow( ( 1 / aLW(sp, 0) ), (1 / aLW(sp, 1) ) ) )  * pow( ( wt(pop_wt_index(sp), sex, age, yr) * 1000), (1 / aLW(sp, 1))); // W = a L ^ b is the same as (W/a)^(1/b)
          }
          if(yr >= nyrs_hind){
            LbyAge( sp, sex, age, yr) = ( pow( ( 1 / aLW(sp, 0) ), (1 / aLW(sp, 1) ) ) )  * pow( ( wt( pop_wt_index(sp), sex, age, (nyrs_hind - 1) ) * 1000), (1 / aLW(sp, 1))); // W = a L ^ b is the same as (W/a)^(1/b)
          }
        }
      }
    }
  }


  // 5.8. Parameter transform
  r_sigma = exp(ln_rec_sigma); // Convert log sd to natural scale
  sigma_srv_index = exp(ln_sigma_srv_index ) ;
  sigma_fsh_catch = exp( ln_sigma_fsh_catch) ;
  sigma_sel = exp(ln_sigma_sel) ;
  sigma_srv_q = exp(ln_sigma_srv_q) ;
  time_varying_sigma_srv_q = exp(ln_sigma_time_varying_srv_q) ;
  sex_ratio_sigma = exp(ln_sex_ratio_sigma);
  Cindex -=1; // Subtract 1 from Cindex to deal with indexing start at 0


  // 5.1. Reorganize survey control bits
  matrix<Type> flt_q(n_flt, nyrs_hind); flt_q.setZero();                        // Vector to save q on natural scale
  vector<int> flt_sel_ind(n_flt); flt_sel_ind.setZero();                        // Vector to store survey index
  vector<int> flt_type(n_flt); flt_type.setZero();                              // Index wether the data are included in the likelihood or not (0 = no, 1 = yes)
  vector<int> flt_sel_type(n_flt); flt_sel_type.setZero();                      // Vector to save survey selectivity type
  vector<int> flt_nselages(n_flt); flt_nselages.setZero();                      // Vector to save number of ages to estimate non-parametric selectivity (1 = age, 2 = length)
  vector<int> flt_varying_sel(n_flt); flt_varying_sel.setZero();                // Vector storing information on wether time-varying selectivity is estimated (0 = no, 1 = random walk with fixed variance, 2 = random effect)
  vector<int> flt_spp(n_flt); flt_spp.setZero();                                // Vector to save survey species
  vector<int> flt_sel_age(n_flt); flt_sel_age.setZero();                        // Vector to save age first selected (selectivity below this age = 0)
  vector<int> flt_accum_age_lower(n_flt); flt_accum_age_lower.setZero();        // Vector to save lower accumulation age
  vector<int> flt_accum_age_upper(n_flt); flt_accum_age_upper.setZero();        // Vector to save upper accumulation age
  vector<int> flt_units(n_flt); flt_units.setZero();                            // Vector to save survey units (1 = weight, 2 = numbers)
  vector<int> flt_wt_index(n_flt); flt_wt_index.setZero();                      // Vector to save 1st dim of wt to use for weight-at-age
  vector<int> flt_age_transition_index(n_flt); flt_age_transition_index.setZero(); // Vector to save 3rd dim of age_trans_matrix to use for ALK
  vector<int> flt_q_ind(n_flt); flt_q_ind.setZero();                            // Vector storing index of survey q for mapping
  vector<int> est_srv_q(n_flt); est_srv_q.setZero();                            // Vector to save wether or not analytical q is used
  vector<int> srv_varying_q(n_flt); srv_varying_q.setZero();                    // Vector storing information on wether time-varying q is estimated (0 = no, 1 = random walk with fixed variance, 2 = random effect)
  vector<int> est_sigma_srv(n_flt); est_sigma_srv.setZero();                    // Vector to save wether sigma survey is estimated
  vector<int> est_sigma_fsh(n_flt); est_sigma_fsh.setZero();                    // Vector to save wether sigma fishery is estimated
  vector<int> sel_norm_age(n_flt); sel_norm_age.setZero();                      // Vector to save age to normalize selectivty


  for (flt_ind = 0; flt_ind < n_flt; flt_ind++){
    flt = fleet_control(flt_ind, 1) - 1;                     // Temporary survey index
    flt_type(flt) = fleet_control(flt_ind, 2);               // Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
    flt_spp(flt) = fleet_control(flt_ind, 3) - 1;            // Species
    flt_sel_ind(flt) = fleet_control(flt_ind, 4);            // Survey selectivity index
    flt_sel_type(flt) = fleet_control(flt_ind, 5);           // Selectivity type
    flt_nselages(flt) = fleet_control(flt_ind, 6);           // Non-parametric selectivity ages
    flt_varying_sel(flt) = fleet_control(flt_ind, 7);        // Time-varying selectivity type.
    flt_sel_age(flt) = fleet_control(flt_ind, 8) - minage(flt_spp(flt));              // First age selected
    flt_accum_age_lower(flt) = fleet_control(flt_ind, 9) - minage(flt_spp(flt));       // Dim1 of wt
    flt_accum_age_upper(flt) = fleet_control(flt_ind, 10) - minage(flt_spp(flt));     // Dim3 of age transition matrix
    flt_units(flt) = fleet_control(flt_ind, 11);             // Survey units
    flt_wt_index(flt) = fleet_control(flt_ind, 12) - 1;      // Dim1 of wt
    flt_age_transition_index(flt) = fleet_control(flt_ind, 13) - 1;     // Dim3 of age transition matrix
    flt_q_ind(flt) = fleet_control(flt_ind, 14) - 1;         // Index of survey q
    est_srv_q(flt) = fleet_control(flt_ind, 15);             // Estimate analytical q?
    srv_varying_q(flt) = fleet_control(flt_ind, 16);         // Time varying q type
    est_sigma_srv(flt) = fleet_control(flt_ind, 17);         // Wether to estimate standard deviation of survey time series
    est_sigma_fsh(flt) = fleet_control(flt_ind, 18);         // Wether to estimate standard deviation of fishery time series
    sel_norm_age(flt) =  flt_nselages(flt) - minage(flt_spp(flt));  // Age to normalize logistic selectivities
  }

  // Set up survey q
  for(flt = 0; flt < n_flt; flt++){
    for(yr = 0; yr < nyrs_hind; yr++){
      // Random walk, penalized deviate
      if(srv_varying_q(flt) != 2){
        srv_q(flt, yr) = exp(ln_srv_q(flt) + ln_srv_q_dev(flt, yr));                 // Exponentiate
      }

      // Random effect
      if(srv_varying_q(flt) == 2){
        // srv_q(flt, yr) = exp(ln_srv_q(flt) + ln_srv_q_dev_re(flt, yr));                 // Exponentiate
      }

      // Q as a function of environmental index
      if(est_srv_q(flt) == 5){
        // srv_q(flt, yr) = exp(ln_srv_q(flt) + srv_q_pow(flt) * env_index_hat(yr, srv_varying_q(flt)));
      }
    }
  }


  matrix<int> r_sexes(UobsWtAge.rows(), 2); r_sexes.setZero();
  matrix<int> k_sexes(UobsWtAge.rows(), 2); k_sexes.setZero();
  // Good above here
  // ------------------------------------------------------------------------- //
  // 6. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0
  // Start iterations


  for (int iter = 0; iter < niter; iter++) {
    // 6.0. EMPIRICAL SELECTIVITY
    sel.setZero();
    for (int sel_ind = 0; sel_ind < emp_sel_obs.rows(); sel_ind++){

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
            for (age = 0; age < nages(sp); age++) {
              if(!isNA(emp_sel_obs(sel_ind, age))){
                sel(flt, sexes(sex), age, yr) = emp_sel_obs(sel_ind, age);
              }
            }
          }
        }
      }
    }      // FIXME - set all empirical selectivities after nyrs_hind to the terminal year

    // 6.1. ESTIMATE SELECTIVITY
    avgsel.setZero();
    sel_tmp.setZero();
    for (flt = 0; flt < n_flt; flt++) {

      // Temporay indices
      sp = flt_spp(flt);             // Temporary index of species
      sel_type = flt_sel_type(flt);
      sel_varying = flt_varying_sel(flt);
      nselages = flt_nselages(flt);


      switch(sel_type){

      case 1:  // 6.1.1. Logisitic selectivity
        for (age = 0; age < nages(sp); age++){
          for (yr = 0; yr < nyrs_hind; yr++) {
            for(sex = 0; sex < nsex(sp); sex++){
              // Random walk and block
              if(sel_varying != 2){
                sel(flt, sex, age, yr) = 1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr)) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr )))));
              }

              // Random effect
              if(sel_varying == 2){
                // sel(flt, sex, age, yr) = 1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev_re(0, flt, sex, yr)) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev_re(0, flt, sex, yr )))));
              }
            }
          }
        }
        break;


      case 2:  // 6.1.2. Non-parametric selectivity fit to age ranges (Ianelli 20??). NOTE: This can likely be improved
        for(sex = 0; sex < nsex(sp); sex++){
          for (age = 0; age < nselages; age++) {
            sel_tmp(flt, sex, age) = sel_coff(flt, sex, age);
            avgsel(flt, sex) +=  exp(sel_coff(flt, sex, age));
          }
          //  Average selectivity up to nselages
          avgsel(flt, sex) = log(avgsel(flt, sex) / nselages);

          // Plus group selectivity
          for (age = nselages; age < nages(sp); age++) {
            sel_tmp(flt, sex, age) = sel_tmp(flt, sex, nselages - 1);
          }

          // Average selectivity across all ages
          avgsel_tmp = 0; // Temporary object for average selectivity across all ages
          for (age = 0; age < nages(sp); age++) {
            avgsel_tmp += exp(sel_tmp(flt, sex, age));
          }
          avgsel_tmp = log(avgsel_tmp / nages(sp));

          // Standardize selectivity
          for (age = 0; age < nages(sp); age++) {
            sel_tmp(flt, sex, age) -= avgsel_tmp;
            sel_tmp(flt, sex, age) = exp(sel_tmp(flt, sex, age));
          }
        }

        // Move to rest of years
        for (age = 0; age < nages(sp); age++){
          for(sex = 0; sex < nsex(sp); sex++){
            for (yr = 0; yr < nyrs_hind; yr++) {
              sel(flt, sex, age, yr) = sel_tmp(flt, sex, age);
            }
          }
        }
        break;


      case 3: // 6.1.3. Double logistic (Dorn and Methot 1990)
        for (age = 0; age < nages(sp); age++){
          for (yr = 0; yr < nyrs_hind; yr++) {
            for(sex = 0; sex < nsex(sp); sex++){

              // Random walk and block
              if(sel_varying != 2){
                sel(flt, sex, age, yr) = 1 / (1 + exp( -exp((ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr))) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr ))))) * // Upper slope
                  (1 - (1 / (1 + exp( -(exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr))) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr )))))));  // Downward slope;
              }

              // Random effect
              if(sel_varying == 2){
                // sel(flt, sex, age, yr) = 1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev_re(0, flt, sex, yr)) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev_re(0, flt, sex, yr ))))) * // Upper slope
                //(1 - (1 / (1 + exp( -exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev_re(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev_re(1, flt, sex, yr )))))));  // Downward slope;
              }
            }
          }
        }
        break;


      case 4: // 6.1.4. Descending logistic (Dorn and Methot 1990)
        for (age = 0; age < nages(sp); age++){
          for (yr = 0; yr < nyrs_hind; yr++) {
            for(sex = 0; sex < nsex(sp); sex++){

              // Random walk and block
              if(sel_varying != 2){
                sel(flt, sex, age, yr) = (1 - (1 / (1 + exp( - exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr )))))));  // Downward slope;
              }

              // Random effect
              if(sel_varying == 2){
                // sel(flt, sex, age, yr) = (1 - (1 / (1 + exp( - exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev_re(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev_re(1, flt, sex, yr )))))));  // Downward slope;
              }
            }
          }
        }
        break;


      case 5: // 6.1.5. Non-parametric selectivity fit to age ranges. Hake version (Taylor et al 2014)
        // -- For each age, sum coefficients from first age selected to age
        for (yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex++){
            for (age = 0; age < nages(sp); age++) {
              sel(flt, sex, age, yr) = 0;
            }
          }
        }

        for (yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex++){
            for (age = flt_sel_age(flt); age < nselages; age++) {
              for (int age_tmp = flt_sel_age(flt); age_tmp <= age; age_tmp++) {
                sel(flt, sex, age, yr) += sel_coff(flt, sex, age_tmp) + sel_coff_dev(flt, sex, age_tmp, yr);
              }
            }
          }
        }

        // -- For each sex & year, subtract max sel across ages from sel and take exp
        for (yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex++){
            Type max_sel = 0;
            for (age = flt_sel_age(flt); age < nselages; age++) {
              // Max
              if(sel(flt, sex, age, yr) > max_sel){
                max_sel = sel(flt, sex, age, yr);
              }
            }

            for (age = flt_sel_age(flt); age < nselages; age++) {
              sel(flt, sex, age, yr) = exp(sel(flt, sex, age, yr) - max_sel);
            }

            // Fill in rest of ages
            for (age = nselages-1; age < nages(sp); age++) {
              sel(flt, sex, age, yr) = sel(flt, sex, nselages-1, yr);
            }
          }
        }
        break;
      } // End selectivity switch


      // 6.1.6. Normalize
      if (sel_type > 0) {

        for (yr = 0; yr < nyrs_hind; yr++) {
          for (age = 0; age < nages(sp); age++){
            for(sex = 0; sex < nsex(sp); sex++){
              if(age < flt_sel_age(flt)){
                sel(flt, sex, age, yr) = 0;
              }
            }
          }
        }



        // Find max for each fishery and year across ages, and sexes
        // FIXME: may not be necessary
        if (sel_type < 5) {
          for (yr = 0; yr < nyrs_hind; yr++) {
            Type max_sel = 0;
            for (age = 0; age < nages(sp); age++){
              for(sex = 0; sex < nsex(sp); sex++){

                // Normalize by specific age
                if((nselages > 0) & (sel_type != 2)){
                  max_sel = sel(flt, sex, sel_norm_age(flt), yr);
                }


                // Normalize by max
                if(sel(flt, sex, age, yr) > max_sel){
                  max_sel = sel(flt, sex, age, yr);
                }
              }
            }

            // Normalize selectivity
            for (age = 0; age < nages(sp); age++){
              for(sex = 0; sex < nsex(sp); sex++){
                sel(flt, sex, age, yr) /= max_sel;
              }
            }
          }
        }
      }
    }



    // 6.1. ESTIMATE HINDCAST FISHING MORTALITY and FSPRs
    F_spp.setZero();
    F_flt_age.setZero();
    F_spp_age.setZero();
    FtargetSPR.setZero();
    FlimitSPR.setZero();
    for (flt = 0; flt < n_flt; flt++) {

      sp = flt_spp(flt);  // Temporary index of fishery survey

      if(flt_type(flt) == 1){
        for (age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            for (yr = 0; yr < nyrs; yr++) {

              // Hindcast
              if( yr < nyrs_hind){
                F_flt_age(flt, sex, age, yr) = sel(flt, sex, age, yr) * exp(ln_mean_F(flt) + F_dev(flt, yr));
              }

              // Forecast
              if( yr >= nyrs_hind){
                // -- Apply HCRs
                switch(HCR){
                case 0: // No fishing
                  proj_F(sp, yr) = 0;
                  break;

                case 2: // Constant F
                  proj_F(sp, yr) = Ftarget(sp, yr);
                  break;

                case 3: // Constant F that acheives X% of SSB0
                  proj_F(sp, yr) = Ftarget(sp, yr);
                  break;

                case 4: // Constant Fspr
                  proj_F(sp, yr) = Ftarget(sp, yr) * Fmult(sp);
                  break;

                case 5: // NPFMC Tier 3 HCR
                  proj_F(sp, yr) = Ftarget(sp, yr); // Used Fabc of FtargetSPR%
                  break;

                case 6: // PFMC Category 1 HCR
                  proj_F(sp, yr) = Flimit(sp, yr) + QnormHCR(sp);
                  Ftarget(sp, yr) = Flimit(sp, yr) + QnormHCR(sp);
                  break;

                case 7: // SESSF Tier 1 HCR
                  proj_F(sp, yr) = Ftarget(sp, yr); // Used Fabc of FtargetSPR%
                  break;
                }


                F_flt_age(flt, sex, age, yr) = sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * proj_F(sp, yr); // FIXME using last year of selectivity
              }

              // -- Sum F across fleets
              F_spp_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);
            }

            // -- Calculate static reference points
            FlimitSPR(sp, sex, age) += sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * Flimit(sp, 0); // FlimitSPR%
            FtargetSPR(sp, sex, age) += sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * Ftarget(sp, 0); // FtargetSPR%
          }
        }

        // F across fleets or species
        for (yr = 0; yr < nyrs; yr++) {
          // Hindcast
          if( yr < nyrs_hind){
            F_flt(flt, yr) = exp(ln_mean_F(flt) + F_dev(flt, yr));
            F_spp(sp, yr) += exp(ln_mean_F(flt) + F_dev(flt, yr)); // Fully selected fishing mortality
          }

          // Forecast
          if( yr >= nyrs_hind){
            F_flt(flt, yr) = proj_F_prop(flt) * proj_F(sp, yr);
            F_spp(sp, yr) +=  proj_F_prop(flt) * proj_F(sp, yr);
          }
        }
      }
    }


    // 6.2. Estimate total mortality at age NOTE: May need to go above population dynamics
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            M(sp, sex, age, yr) = M1(sp, sex, age) + M2(sp, sex, age, yr);
            Zed(sp, sex, age, yr) = M1(sp, sex, age) + F_spp_age(sp, sex, age, yr) + M2(sp, sex, age, yr);
            S(sp, sex, age, yr) = exp(-Zed(sp, sex, age, yr));
          }
        }
      }
    }


    // 6.4. STOCK-RECRUIT PARAMETERS
    // -- For beverton-holt, steepness and R0 are derived from SPR0
    for ( sp = 0; sp < nspp ; sp++) {
      switch(srr_fun){
      case 0: // Random about mean (e.g. Alaska)
        Steepness(sp) = 0.99;
        Rinit(sp) = R0(sp) = exp(rec_pars(sp, 0));
        break;

      case 1: // Beverton-Holt
        Steepness(sp) = exp(rec_pars(sp, 0)) * SPR0(sp)/(4.0 + exp(rec_pars(sp, 0)) * SPR0(sp));
        R0(sp) = (exp(rec_pars(sp, 0))-1/SPR0(sp)) / exp(rec_pars(sp, 1)); // (Alpha-1/SPR0)/beta
        // Rinit(sp) = (exp(rec_pars(sp, 0))-1/SPRFinit(sp)) / exp(rec_pars(sp, 1)); // (Alpha-1/SPR0)/beta
        break;

      case 2: // Beverton-Holt with environmental impacts on alpha
        //FIXME make time-varying
        Steepness(sp) = exp(rec_pars(sp, 0)) * SPR0(sp)/(4.0 + exp(rec_pars(sp, 0)) * SPR0(sp));
        R0(sp) = (exp(rec_pars(sp, 0))-1/SPR0(sp)) / exp(rec_pars(sp, 1)); // (Alpha-1/SPR0)/beta
        // Rinit(sp) = (exp(rec_pars(sp, 0))-1/SPRFinit(sp)) / exp(rec_pars(sp, 1)); // (Alpha-1/SPR0)/beta
        break;

      case 3: // Ricker
        Steepness(sp) = 0.2 * exp(0.8*log(exp(rec_pars(sp, 0)) * SPR0(sp))); //
        R0(sp) = log(exp(rec_pars(sp, 0)) * SPR0(sp))/(exp(rec_pars(sp, 1)) * SPR0(sp)); // FIXME - make time-varying
        // Rinit(sp) = log(exp(rec_pars(sp, 0)) * SPRFinit(sp))/(exp(rec_pars(sp, 1)) * SPRFinit(sp)); // FIXME - make time-varying
        break;

      default:
        error("Invalid 'srr_fun'");
      }
    }


    // 6.5. INITIAL ABUNDANCE AT AGE, BIOMASS, AND SSB (YEAR 1)
    biomass.setZero();
    biomassSSB.setZero();
    biomassByage.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++){
        for(sex = 0; sex < nsex(sp); sex ++){

          switch(estDynamics(sp)){
          case 0: // Estimated

            // - Estimate as free parameters
            if(initMode == 0){
              R(sp, 0) = exp(init_dev(sp, 0));
              NByage(sp, 0, age, 0) = exp(init_dev(sp, age)) * R_sexr(sp);
              NByage(sp, 1, age, 0) = exp(init_dev(sp, age)) * (1-R_sexr(sp));
            }

            // - Equilibrium or non-equilibrium estimated as function of R0, Finit, mortality, and init devs
            // Finit is set to 0 when initMode != 2
            if(initMode > 0){
              // -- 6.5.1. Amin (i.e. recruitment)
              if(age == 0){
                R(sp, 0) = Rinit(sp) * exp(rec_dev(sp, 0));
                NByage(sp, 0, 0, 0) = R(sp, 0) * R_sexr(sp);
                NByage(sp, 1, 0, 0) = R(sp, 0) * (1-R_sexr(sp));
              }

              // Sum M1 until age - 1
              Type mort_sum = 0;
              for(int age_tmp = 0; age_tmp < age; age_tmp++){
                mort_sum += M1(sp, sex, age_tmp) ;// + Finit(sp);
              }

              // -- 6.5.2. Age Amin+1:Amax-1 (initial abundance)
              if ((age > 0) & (age < nages(sp) - 1)) {

                if(sex == 0){
                  NByage(sp, 0, age, 0) = R0(sp) * exp( - mort_sum + init_dev(sp, age - 1)) * R_sexr(sp);
                }
                if(sex == 1){
                  NByage(sp, 1, age, 0) = R0(sp) * exp( - mort_sum + init_dev(sp, age - 1)) * (1-R_sexr(sp));
                }
              }

              // -- 6.5.3. Amax
              if (age == (nages(sp) - 1)) {

                if(sex ==0){// NOTE: This solves for the geometric series
                  NByage(sp, 0, age, 0) = R0(sp) * exp( - mort_sum + init_dev(sp, age - 1)) / (1 - exp(-M1(sp, sex, nages(sp) - 1))) * R_sexr(sp);
                }

                if(sex == 1){
                  NByage(sp, 1, age, 0) = R0(sp) * exp( - mort_sum + init_dev(sp, age - 1)) / (1 - exp(-M1(sp, sex, nages(sp) - 1))) * (1-R_sexr(sp));
                }
              }
            }
            break;

          case 1: // Fixed numbers-at-age - fixed scalar
            NByage(sp, sex, age, 0) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, 0);
            break;

          case 2: // Fixed numbers-at-age age-independent scalar
            NByage(sp, sex, age, 0) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, 0);
            break;

          case 3: // Fixed numbers-at-age age-dependent scalar
            NByage(sp, sex, age, 0) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, 0);
            break;

          default:
            error("Invalid 'estDynamics'");
          }

          // -- 6.5.3. Estimate total biomass in year 1
          biomassByage(sp, sex, age, 0) = NByage(sp, sex, age, 0) * wt( pop_wt_index(sp), sex, age, 0);
          biomass(sp, 0) += biomassByage(sp, sex, age, 0);
        }

        // -- 6.5.4. Estimated initial female SSB
        biomassSSBByage(sp, age, 0) = NByage(sp, 0, age, 0) * pow(S(sp, 0, age, 0), spawn_month(sp)/12) * wt(ssb_wt_index(sp), 0, age, 0 ) * pmature(sp, age); // 6.6.
        biomassSSB(sp, 0) += biomassSSBByage(sp, age, 0);
      }
    }





    // 6.4. ESTIMATE HINDCAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
    biomass.setZero();
    biomassSSB.setZero();
    biomassByage.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 1; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- 6.6.1. Amin (i.e. recruitment)
              R(sp, yr) = R0(sp) * exp(rec_dev(sp, yr));
              NByage(sp, 0, 0, yr) = R(sp, yr) * R_sexr(sp);
              NByage(sp, 1, 0, yr) = R(sp, yr) * (1-R_sexr(sp));

              // -- 6.6.2.  Where Amin < age < Amax
              if (age < (nages(sp) - 1)) {
                NByage(sp, sex, age + 1, yr) = NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }

              // -- 6.6.3. Plus group where age > Amax. NOTE: This is not the same as T1.3 because sp used age = A_i rather than age > A_i.
              if (age == (nages(sp) - 1)) {
                NByage(sp, sex, age, yr) = NByage(sp, sex, age - 1, yr - 1) * S(sp, sex, age - 1, yr - 1) + NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }
              break;
            case 1: // Fixed numbers-at-age - fixed scalar
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;
            case 2: // Fixed numbers-at-age age-independent scalar
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;
            case 3: // Fixed numbers-at-age age-dependent scalar
              NByage(sp, sex, age, yr) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, yr);
              break;
            default:
              error("Invalid 'estDynamics'");
            }
          }
        }

        // -- 6.6.4. Estimate Biomass and SSB
        for (yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            biomassByage(sp, sex, age, yr) = NByage(sp, sex, age, yr) * wt( pop_wt_index(sp), sex, age, yr ); // 6.5.
            biomassSSBByage(sp, age, yr) = NByage(sp, 0, age, yr) * pow(S(sp, 0, age, yr), spawn_month(sp)/12) * wt(ssb_wt_index(sp), 0, age, yr ) * pmature(sp, age); // 6.6.
            biomass(sp, yr) += biomassByage(sp, sex, age, yr);
          }

          biomassSSB(sp, yr) += biomassSSBByage(sp, age, yr);
        }
      }
    }


    // 6.5. Static SPR based reference points
    // -- calculate mean recruitment
    mean_rec.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs_mean; yr++) {
        mean_rec(sp) += R(sp, yr)/nyrs_mean; // Update mean rec
      }
    }

    SPR0.setZero();
    SPRlimit.setZero();
    SPRtarget.setZero();
    for (sp = 0; sp < nspp; sp++) {

      if(nsex(sp)  == 1){
        R_sexr(sp) = 1; // We multiply by pmature later on which is sex-ratio * maturity for 1 sex and maturity for 2-sex
      }

      NbyageSPR(0, sp, 0) = mean_rec(sp) * R_sexr(sp); // F = 0
      NbyageSPR(1, sp, 0) = mean_rec(sp) * R_sexr(sp); // F = Flimit
      NbyageSPR(2, sp, 0) = mean_rec(sp) * R_sexr(sp); // F = Ftarget

      for (age = 1; age < nages(sp)-1; age++) {
        NbyageSPR(0, sp, age) =  NbyageSPR(0, sp, age-1) * exp(-M(sp, 0, age-1, nyrs_hind - 1));
        NbyageSPR(1, sp, age) =  NbyageSPR(1, sp, age-1) * exp(-(M(sp, 0, age-1, nyrs_hind - 1) + FlimitSPR(sp, 0, age-1)));
        NbyageSPR(2, sp, age) =  NbyageSPR(2, sp, age-1) * exp(-(M(sp, 0, age-1, nyrs_hind - 1) + FtargetSPR(sp, 0, age-1)));
      }

      // Plus group
      NbyageSPR(0, sp, nages(sp) - 1) = NbyageSPR(0, sp, nages(sp) - 2) * exp(-M(sp, 0, nages(sp) - 2, nyrs_hind - 1)) / (1 - exp(-M(sp, 0, nages(sp) - 1, nyrs_hind - 1)));
      NbyageSPR(1, sp, nages(sp) - 1) = NbyageSPR(1, sp, nages(sp) - 2) * exp(-(M(sp, 0,  nages(sp) - 2, nyrs_hind - 1) + FlimitSPR(sp, 0,  nages(sp) - 2))) / (1 - exp(-(M(sp, 0,  nages(sp) - 1, nyrs_hind - 1) + FlimitSPR(sp, 0,  nages(sp) - 1))));
      NbyageSPR(2, sp, nages(sp) - 1) = NbyageSPR(2, sp, nages(sp) - 2) * exp(-(M(sp, 0,  nages(sp) - 2, nyrs_hind - 1) + FtargetSPR(sp, 0,  nages(sp) - 2))) / (1 - exp(-(M(sp, 0,  nages(sp) - 1, nyrs_hind - 1) + FtargetSPR(sp, 0,  nages(sp) - 1))));

      // Calulcate SPR
      for (age = 0; age < nages(sp); age++) {
        SPR0(sp) +=  NbyageSPR(0, sp, age) *  wt( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * pmature(sp, age) * exp(-M(sp, 0, age, nyrs_hind - 1) * spawn_month(sp)/12);
        SPRlimit(sp) +=  NbyageSPR(1, sp, age) *  wt( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * pmature(sp, age) * exp(-(M(sp, 0,  age, nyrs_hind - 1) + FlimitSPR(sp, 0,  age)) * spawn_month(sp)/12);
        SPRtarget(sp) +=  NbyageSPR(2, sp, age) *  wt( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * pmature(sp, age) * exp(-(M(sp, 0,  age, nyrs_hind - 1) + FtargetSPR(sp, 0,  age)) * spawn_month(sp)/12);
      }
    }


    // 6.6. Dynamic SPR based reference points and SB0
    B0.setZero();
    SB0.setZero();
    DynamicB0.setZero();
    DynamicSB0.setZero();
    DynamicSBF.setZero();
    for (sp = 0; sp < nspp; sp++) {
      // Year 1
      // Initial abundance-at-age is the same as the hindcast
      for(sex = 0; sex < nsex(sp); sex ++){
        for(age = 0; age < nages(sp); age++){
          NByage0(sp, sex, age , 0) = DynamicNByage0(sp, sex, age, 0) = NByage(sp,sex,age,0);
        }
      }

      // -- Nbyage0
      for (yr = 1; yr < nyrs; yr++) {

        // Recruitment - SB0 and Dynamic SB0
        // FIXME account for S-R relationship down the line
        NByage0(sp, 0, 0 , yr) = mean_rec(sp) * R_sexr(sp); // Using mean rec because R0 is biased
        NByage0(sp, 1, 0 , yr) = mean_rec(sp) * (1 - R_sexr(sp)); // Using mean rec because R0 is biased

        // Recruitment - Dynamic SPR
        if(yr < nyrs_hind){ // Hindcast using R0
          DynamicNByage0(sp, 0, 0, yr) = R0(sp) * exp( rec_dev(sp, yr)) * R_sexr(sp); // Using R0 because mse functions adjust rec_devs to account for bias
          DynamicNByage0(sp, 1, 0, yr) = R0(sp) * exp( rec_dev(sp, yr)) * (1-R_sexr(sp)); // Using R0 because mse functions adjust rec_devs to account for bias
        }
        if(yr >= nyrs_hind){ // Forecast
          // - Use mean rec
          if(proj_mean_rec == 1){
            DynamicNByage0(sp, 0, 0, yr) = exp(log(mean_rec(sp)) + rec_dev(sp, yr)) * R_sexr(sp); // Using R0 because mse functions adjust rec_devs to account for bias
            DynamicNByage0(sp, 1, 0, yr) = exp(log(mean_rec(sp)) + rec_dev(sp, yr)) * (1-R_sexr(sp)); // Using R0 because mse functions adjust rec_devs to account for bias
          }
          // -- Use rec dev
          if(proj_mean_rec == 0){
            DynamicNByage0(sp, 0, 0 , yr) = R0(sp) * exp( rec_dev(sp, yr)) * R_sexr(sp);
            DynamicNByage0(sp, 1, 0 , yr) = R0(sp) * exp( rec_dev(sp, yr)) * (1 - R_sexr(sp));
          }
        }

        // N-at-age -- No fishing
        for(sex = 0; sex < nsex(sp); sex ++){
          for (age = 1; age < nages(sp)-1; age++) {
            NByage0(sp, sex, age, yr) =  NByage0(sp, sex, age-1, yr-1) * exp(-M(sp, sex, age-1, yr-1)); // F = 0
            DynamicNByage0(sp, sex, age, yr) =  DynamicNByage0(sp, sex, age-1, yr-1) * exp(-M(sp, sex, age-1, yr - 1)); // F = 0
          }

          // Plus group  -- No fishing
          NByage0(sp, sex, nages(sp)-1, yr) = NByage0(sp, sex, nages(sp)-2, yr - 1) * exp(-M(sp, sex, nages(sp)-2, yr - 1))  + NByage0(sp, sex, nages(sp)-1, yr - 1) * exp(-M(sp, sex, nages(sp)-1, yr - 1));
          DynamicNByage0(sp, sex, nages(sp)-1, yr) = DynamicNByage0(sp, sex, nages(sp)-2, yr - 1) * exp(-M(sp, sex, nages(sp)-2, yr - 1))  + DynamicNByage0(sp, sex, nages(sp)-1, yr - 1) * exp(-M(sp, sex, nages(sp)-1, yr - 1));

        }
      }


      // -- NbyageF
      for (int F_yr = 1; F_yr < nyrs; F_yr++) { // Loop across years for F
        for (yr = 0; yr <= F_yr; yr++) { // Loop across all years prior to dynamic F year

          if(yr < nyrs_hind){
            yr_ind = yr;
          }
          if(yr >= nyrs_hind){
            yr_ind = nyrs_hind - 1;
          }

          if(yr == 0){
            DynamicNByageF(F_yr, sp, age, 0) = NByage(sp,0,age,0);
          }

          if(yr > 0){
            // Recruitment - Dynamic SBF
            // FIXME account for S-R relationship down the line
            if(yr < nyrs_hind){ // Hindcast using R0
              DynamicNByageF(F_yr, sp, 0, yr) = R0(sp) * exp( rec_dev(sp, yr)) * R_sexr(sp); // Using R0 because mse_functions account for bias

            }
            if(yr >= nyrs_hind){ // Forecast using mean R
              if(proj_mean_rec == 1){
                DynamicNByageF(F_yr, sp, 0, yr) = exp(log(mean_rec(sp)) + rec_dev(sp, yr)) * R_sexr(sp); // Using R0 because mse_functions account for bias
              }
              if(proj_mean_rec == 0){
                DynamicNByageF(F_yr, sp, 0, yr) = R0(sp) * exp( rec_dev(sp, yr)) * R_sexr(sp);
              }
            }

            // N-at-age -- With fishing
            for (flt = 0; flt < n_flt; flt++) {
              if(flt_spp(flt) == sp){
                for (age = 1; age < nages(sp)-1; age++) {
                  DynamicNByageF(F_yr, sp, age, yr) =  DynamicNByageF(F_yr, sp, age-1, yr-1) * exp(-M(sp, 0, age-1, yr - 1) - sel(flt, 0, age-1, yr_ind-1) * proj_F_prop(flt) * Ftarget(sp, F_yr)); // F = 0
                }

                // Plus group  -- No fishing
                DynamicNByageF(F_yr, sp, nages(sp)-1, yr) = DynamicNByageF(F_yr, sp, nages(sp)-2, yr - 1) * exp(-M(sp, 0, nages(sp)-2, yr - 1) - sel(flt, 0, nages(sp)-2, yr_ind - 1) * proj_F_prop(flt) * Ftarget(sp, F_yr))
                  + DynamicNByageF(F_yr, sp, nages(sp)-1, yr - 1) * exp(-M(sp, 0, nages(sp)-1, yr - 1) - sel(flt, 0, nages(sp)-1, yr_ind - 1) * proj_F_prop(flt) * Ftarget(sp, F_yr));
              }
            }
          }
        }
      }


      // Calulcate Dynamic SB0 and Dynamic SBF
      for (yr = 1; yr < nyrs; yr++) { // No initial
        if(yr < nyrs_hind){
          yr_ind = yr;
        }
        if(yr >= nyrs_hind){
          yr_ind = nyrs_hind - 1;
        }
        for (age = 0; age < nages(sp); age++) {
          DynamicSB0(sp, yr) +=  DynamicNByage0(sp, 0, age, yr) *  wt( ssb_wt_index(sp), 0, age, yr_ind ) * pmature(sp, age) * exp(-M(sp, 0, age, yr) * spawn_month(sp)/12);

          for(sex = 0; sex < nsex(sp); sex ++){
            DynamicB0(sp, yr) +=  DynamicNByage0(sp, sex,age, yr) *  wt( ssb_wt_index(sp), sex, age, yr_ind );
          }

          // Dynamic SB with F (loop across fleets)
          for (flt = 0; flt < n_flt; flt++) {
            if(flt_spp(flt) == sp){
              DynamicSBF(sp, yr) +=  DynamicNByageF(yr, sp, age, yr) *  wt( ssb_wt_index(sp), 0, age, yr_ind ) * pmature(sp, age) * exp(-(M(sp, 0, age, yr) + sel(flt, 0, age, yr_ind) * proj_F_prop(flt) * Ftarget(sp, yr)) * spawn_month(sp)/12); // Fixme - add F
            }
          }
        }
      }

      // Calculate SB0
      for (age = 0; age < nages(sp); age++) {
        SB0(sp) +=  NByage0(sp, 0, age, nyrs-1) *  wt( ssb_wt_index(sp), 0, age, nyrs_hind - 1 ) * pmature(sp, age) * exp(-M(sp, 0, age, nyrs-1) * spawn_month(sp)/12);
        for(sex = 0; sex < nsex(sp); sex ++){
          B0(sp) +=  NByage0(sp, sex, age, nyrs-1) *  wt( pop_wt_index(sp), sex, age, nyrs_hind - 1 );
        }
      }

      // Input SB0 (if running in multi-species mode)
      if(msmMode > 0){
        SB0(sp) = MSSB0(sp);
      }
    }


    // 6.7. HCR and FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
    // Includes Harvest Control Rules
    for (sp = 0; sp < nspp; sp++) {
      for (yr = nyrs_hind; yr < nyrs; yr++){

        // -- Apply HCRs
        // Static Harvest Control Rules
        if(DynamicHCR == 0){
          switch(HCR){
          case 0: // No fishing
            proj_F(sp, yr) = 0;
            break;

          case 2: // Constant F
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 3: // Constant F to acheive X% of SB0
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 4: // Constant Fspr
            proj_F(sp, yr) = proj_F(sp, yr) * Fmult(sp);
            break;

          case 5: // NPFMC Tier 3 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(biomassSSB(sp, nyrs_hind - 1) < SPRtarget(sp)){
              proj_F(sp, yr) = Ftarget(sp, 0) * (((biomassSSB(sp, nyrs_hind-1)/SPRtarget(sp))-Alpha(sp))/(1-Alpha(sp))); // Used Fabc of FtargetSPR%
            }
            if((biomassSSB(sp, nyrs_hind-1) < SB0(sp) * Plimit(sp)) | (biomassSSB(sp, nyrs_hind-1) / SPRtarget(sp) < Alpha(sp))){ // If overfished
              proj_F(sp, yr) =  0;
            }
            break;

          case 6: // PFMC Category 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(biomassSSB(sp, nyrs_hind-1) < SB0(sp) * Ptarget(sp)){
              proj_F(sp, yr) = (Flimit(sp, 0) + QnormHCR(sp)) * (SB0(sp) * Ptarget(sp) * (biomassSSB(sp, nyrs_hind-1) - SB0(sp) * Plimit(sp))) / (biomassSSB(sp, nyrs_hind-1) * (SB0(sp) * (Ptarget(sp) - Plimit(sp))));
            }
            if(biomassSSB(sp, nyrs_hind-1) < SB0(sp) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) =  0;
            }
            break;

          case 7: // SESSF Tier 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(biomassSSB(sp, nyrs_hind-1) < SB0(sp) * Ptarget(sp)){
              proj_F(sp, yr) = Ftarget(sp, 0) * ((biomassSSB(sp, nyrs_hind-1)/(SB0(sp) * Plimit(sp)))-1); // Used Fabc of FtargetSPR%
            }
            if(biomassSSB(sp, nyrs_hind-1) < SB0(sp) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) =  0;
            }
            break;
          }
        }

        // Dynamic Harvest Control Rules
        if(DynamicHCR == 1){
          switch(HCR){
          case 0: // No fishing
            proj_F(sp, yr) = 0;
            break;

          case 2: // Constant F
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 3: // Constant F to acheive X% of SB0
            proj_F(sp, yr) = proj_F(sp, yr);
            break;

          case 4: // Constant Fspr
            proj_F(sp, yr) = proj_F(sp, yr) * Fmult(sp);
            break;

          case 5: // NPFMC Tier 3 HCR
            //FIXME: Using DynamicSB0 * Ptarget(sp) rather than DynamicSPRtarget because they are the same with no S-R curve
            proj_F(sp, yr) = proj_F(sp, yr);
            if(biomassSSB(sp, nyrs_hind-1) < DynamicSB0(sp, nyrs_hind-1) * FsprTarget(sp)){
              proj_F(sp, yr) = Ftarget(sp, yr) * (((biomassSSB(sp, nyrs_hind-1)/(DynamicSB0(sp, nyrs_hind-1) * Ptarget(sp)))-Alpha(sp))/(1-Alpha(sp))); // Used Fabc of FtargetSPR%
            }
            if((biomassSSB(sp, nyrs_hind-1) < DynamicSB0(sp, nyrs_hind-1) * Plimit(sp)) | (biomassSSB(sp, nyrs_hind-1) / (DynamicSB0(sp, nyrs_hind-1) * Ptarget(sp)) < Alpha(sp))){ // If overfished
              proj_F(sp, yr) =  0;
            }
            break;

          case 6: // PFMC Category 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(biomassSSB(sp, nyrs_hind-1) < DynamicSB0(sp, nyrs_hind-1) * Ptarget(sp)){
              proj_F(sp, yr) = (Flimit(sp, yr) + QnormHCR(sp)) * (DynamicSB0(sp, nyrs_hind-1) * Ptarget(sp) * (biomassSSB(sp, nyrs_hind-1) - DynamicSB0(sp, nyrs_hind-1) * Plimit(sp))) / (biomassSSB(sp, nyrs_hind-1) * (DynamicSB0(sp, nyrs_hind-1) * (Ptarget(sp) - Plimit(sp))));
            }
            if(biomassSSB(sp, nyrs_hind-1) < DynamicSB0(sp, nyrs_hind-1) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) =  0;
            }
            break;

          case 7: // SESSF Tier 1 HCR
            proj_F(sp, yr) = proj_F(sp, yr);
            if(biomassSSB(sp, nyrs_hind-1) < DynamicSB0(sp, nyrs_hind-1) * Ptarget(sp)){
              proj_F(sp, yr) = Ftarget(sp, yr) * ((biomassSSB(sp, nyrs_hind-1)/(DynamicSB0(sp, nyrs_hind-1) * Plimit(sp)))-1); // Used Fabc of FtargetSPR%
            }
            if(biomassSSB(sp, nyrs_hind-1) < DynamicSB0(sp, nyrs_hind-1) * Plimit(sp)){ // If overfished
              proj_F(sp, yr) =  0;
            }
            break;
          }
        }


        // -- Update F for the projection
        for (age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            F_spp_age(sp, sex, age, yr) = 0;
          }
        }

        F_spp(sp, yr) = proj_F(sp, yr);
        for (flt = 0; flt < n_flt; flt++) {
          if(sp == flt_spp(flt)){
            F_flt(sp, yr) = proj_F_prop(flt) * proj_F(sp, yr);
            for (age = 0; age < nages(sp); age++) {
              for(sex = 0; sex < nsex(sp); sex ++){
                F_flt_age(flt, sex, age, yr) = sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * proj_F(sp, yr); // FIXME using last year of selectivity
                if(flt_type(flt) == 1){
                  F_spp_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);
                }
              }
            }
          }
        }

        for (age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            biomassByage(sp, sex, age, yr) = 0;
            M(sp, sex, age, yr) = M1(sp, sex, age) + M2(sp, sex, age, yr);
            Zed(sp, sex, age, yr) = M1(sp, sex, age) + F_spp_age(sp, sex, age, yr) + M2(sp, sex, age, yr);
            S(sp, sex, age, yr) = exp(-Zed(sp, sex, age, yr));
          }
        }


        // -- Estimate NbyAge and B/SSB for the projection given the HCR
        biomass(sp, yr) = 0;
        biomassSSB(sp, yr) = 0;
        for (age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- 6.7.1. Amin (i.e. recruitment)
              // - Use mean rec
              if(proj_mean_rec == 1){
                R(sp, yr) = exp(log(mean_rec(sp)) + rec_dev(sp, yr)); //  exp(+ rec_dev(sp, yr)); // Projections use mean R given bias in R0
                NByage(sp, 0, 0, yr) = R(sp, yr) * R_sexr(sp);
                NByage(sp, 1, 0, yr) = R(sp, yr) * (1-R_sexr(sp));
              }

              // - Use S-R/rec devs
              if(proj_mean_rec == 0){
                R(sp, yr) = R0(sp) * exp(rec_dev(sp, yr));
                NByage(sp, 0, 0 , yr) = R0(sp) * exp(rec_dev(sp, yr)) * R_sexr(sp);
                NByage(sp, 1, 0 , yr) = R0(sp) * exp(rec_dev(sp, yr)) * (1 - R_sexr(sp));
              }

              // -- 6.7.2.  Where Amin < age < Amax
              if (age < (nages(sp) - 1)) {
                NByage(sp, sex, age + 1, yr) = NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }

              // -- 6.7.3. Plus group where age > Ai. NOTE: This is not the same as T1.3 because sp used age = A_i rather than age > A_i.
              if (age == (nages(sp) - 1)) {
                NByage(sp, sex, age, yr) = NByage(sp, sex, age - 1, yr - 1) * S(sp, sex, age - 1, yr - 1) + NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }
              break;
            case 1: // Fixed numbers-at-age - fixed scalar
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;
            case 2:// Fixed numbers-at-age age-independent scalar
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
              break;
            case 3:// Fixed numbers-at-age age-dependent scalar
              NByage(sp, sex, age, yr) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, yr);
              break;
            default:
              error("Invalid 'estDynamics'");
            }

            // Hard constraint to reduce population collapse
            if(NByage(sp, sex, age, yr) < minNByage){
              NByage(sp, sex, age, yr) = minNByage;
            }

            // -- 6.3.3. FORECAST Biomass
            biomassByage(sp, sex, age, yr) = NByage(sp, sex, age, yr) * wt( pop_wt_index(sp), sex, age, nyrs_hind-1 ); // 6.5.
            biomass(sp, yr) += biomassByage(sp, sex, age, yr);
          }

          // -- 6.3.3. FORECAST SSB
          biomassSSBByage(sp, age, yr) = NByage(sp, 0, age, yr) * pow(S(sp, 0, age, yr), spawn_month(sp)/12) * wt(ssb_wt_index(sp), 0, age, nyrs_hind-1 ) * pmature(sp, age); // 6.6.
          biomassSSB(sp, yr) += biomassSSBByage(sp, age, yr);
        }
      }
    }


    // 6.8. ESTIMATE SEX RATIO
    sex_ratio_hat.setZero();
    Type sexr_denom = 0; // Denominator for sex ration across ages
    Type sexr_nom = 0; // Nominator for sex ratio across ages
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs; yr++) {
        // Re-initialize
        sexr_denom = 0;
        sexr_nom = 0;
        for (age = 0; age < nages(sp); age++) {
          sex_ratio_hat(sp, age, yr) = NByage(sp, 0, age, yr) / (NByage(sp, 0, age, yr) + NByage(sp, 1, age, yr));
          sexr_nom +=  NByage(sp, 0, age, yr); // Females
          sexr_denom += (NByage(sp, 0, age, yr) + NByage(sp, 1, age, yr)); // Males
        }
        sex_ratio_mean_hat(sp, yr) = sexr_nom/sexr_denom;
      }
    }



    // 6.9. ESTIMATE AVERAGE NUMBERS AT AGE
    for (sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for (age = 0; age < nages(sp); age++) {
          for (yr = 0; yr < nyrs; yr++) {
            switch(avgnMode){
            case 0: // MSVPA approach
              AvgN(sp, sex, age, yr) = NByage(sp, sex, age, yr) * (1 - S(sp, sex, age, yr)) / Zed(sp, sex, age, yr);
              break;
            case 1: // Kinzey and Punt (2009) approximation
              AvgN(sp, sex, age, yr) = NByage(sp, sex, age, yr) * exp(- Zed(sp, sex, age, yr) / 2);
              break;
            case 2: // Van Kirk et al (2010) approximation
              AvgN(sp, sex, age, yr) = NByage(sp, sex, age, yr);
              break;
            default:
              error("Invalid 'avgnMode'");
            }
          }
        }
      }
    }


    // ------------------------------------------------------------------------- //
    // 7. RATION EQUATIONS                                                       //
    // ------------------------------------------------------------------------- //
    // NOTE -- LOOPING INDICES -- sp = species, age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length

    // 7.2. Calculate temperature function of consumption
    Type Yc = 0;
    Type Zc = 0;
    Type Vc = 0;
    Type Xc = 0;
    Type G2 = 0;
    Type L2 = 0;
    Type G1 = 0;
    Type L1 = 0;
    Type Ka = 0;
    Type Kb = 0;
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs; yr++) {
        switch(Ceq(sp)){
        case 1:// Exponential function from Stewart et al. 1983
          fT(sp, yr) = exp(Qc(sp) * env_index_hat(yr, Cindex(sp)));
          break;

        case 2:// Temperature dependence for warm-water-species from Kitchell et al. 1977
          Yc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp) + 2);
          Zc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp));
          Vc = (Tcm(sp) - env_index_hat(yr, Cindex(sp))) / (Tcm(sp) - Tco(sp));
          Xc = pow(Zc, 2) * pow((1 + pow((1 + 40 / Yc), 0.5)), 2) / 400;
          fT(sp, yr) = pow(Vc, Xc) * exp(Xc * (1 - Vc));
          break;

        case 3:// Temperature dependence for cool and cold-water species from Thornton and Lessem 1979
          G2 = (1 / (Tcl(sp) - Tcm(sp))) * log((0.98 * (1 - CK4(sp))) / (CK4(sp) * 0.02));
          L2 = exp(G2 * (Tcl( sp ) -  env_index_hat(yr, Cindex(sp))));
          Kb = (CK4(sp) * L2) / (1 + CK4(sp) * (L2 - 1));
          G1 = (1 / (Tco(sp) - Qc(sp))) * log((0.98 * (1 - CK1(sp))) / (CK1(sp) * 0.02));
          L1 = exp(G1 * (env_index_hat(yr, Cindex(sp)) - Qc(sp)));
          Ka = (CK1(sp) * L1) / (1 + CK1(sp) * (L1 - 1));
          fT(sp, yr) = Ka * Kb;
          break;

        case 4:
          fT(sp, yr) = 1;
          break;
        }
      }
    }


    // 7.3. Calculate historic ration
    for (sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex++) {
        for (age = 0; age < nages(sp); age++) {
          for (yr = 0; yr < nyrs; yr++) {
            // p = proportion of maximum consumption
            // f(T) = temperature dependence function
            // CA = intercept of allometric mass function
            // CB = slope of allometric mass function
            // fday = number of forageing days per year

            // Hindcast
            if(yr < nyrs_hind){
              ConsumAge(sp, sex, age, yr) = CA(sp) * pow(wt( pop_wt_index(sp), sex, age, yr ) * Type(1000), 1 + CB( sp )) //  C_max = CA * W ^ 1+CB; where C_max is grams consumed per grams of predator per day
              * fT(sp, yr) * fday( sp );                           //  C_max * f(T) * wt * fday g/pred.yr
              ConsumAge(sp, sex, age, yr) = ConsumAge(sp, sex, age, yr) * Pvalue(sp) * Pyrs(sp, sex, age, yr); //
            }

            // Projection
            if(yr >= nyrs_hind){
              ConsumAge(sp, sex, age, yr) = CA(sp) * pow(wt( pop_wt_index(sp), sex, age, (nyrs_hind - 1) ) * Type(1000), 1 + CB( sp ))  //  C_max = CA * W ^ 1+CB; where C_max is grams consumed per grams of predator per day
              * fT(sp, yr) * fday( sp );                            //  C_max * f(T) * wt * fday g/pred.yr
              ConsumAge(sp, sex, age, yr) = ConsumAge(sp, sex, age, yr) * Pvalue(sp) * Pyrs(sp, sex, age, (nyrs_hind-1)); //
            }

            ration(sp, sex, age, yr) = ConsumAge(sp, sex, age, yr) / 1000;      // Annual ration kg/yr //aLW(predd)*pow(lengths(predd,age),bLW(predd));//mnwt_bin(predd,age);
          }
        }
      }
    }


    // 7.4. Reorganize UobsWTAge content
    for(int stom_ind = 0; stom_ind < UobsWtAge.rows(); stom_ind++){
      rsp = UobsWtAge_ctl(stom_ind, 0) - 1; // Index of pred
      ksp = UobsWtAge_ctl(stom_ind, 1) - 1; // Index of prey
      r_sex = UobsWtAge_ctl(stom_ind, 2); // Index of pred sex
      k_sex = UobsWtAge_ctl(stom_ind, 3); // Index of prey sex
      r_age = UobsWtAge_ctl(stom_ind, 4) - minage(rsp); // Index of pred age
      k_age = UobsWtAge_ctl(stom_ind, 5) - minage(ksp); // Index of prey age
      flt_yr = UobsWtAge_ctl(stom_ind, 6); // Index of year

      // Predator
      // 1 sex model
      r_sexes(stom_ind, 0) = 0; r_sexes(stom_ind, 1) = 1;
      k_sexes(stom_ind, 0) = 0; k_sexes(stom_ind, 1) = 1;

      if(r_sex > 0){
        r_sexes(stom_ind, 0) = r_sex - 1;  r_sexes(stom_ind, 1) = r_sex - 1;
      }

      if(k_sex > 0){
        k_sexes(stom_ind, 0) = k_sex - 1;  k_sexes(stom_ind, 1) = k_sex - 1;
      }


      for(int j = 0; j < 2; j ++){
        for(int k = 0; k < 2; k ++){

          if(flt_yr > 0){
            yr = flt_yr - styr;

            if(yr < nyrs_hind){
              diet_prop_weight(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr) = UobsWtAge(stom_ind, 1);
            }
          }

          // Average of years
          if(flt_yr == 0){
            for (yr = 0; yr < nyrs; yr++) {
              diet_prop_weight(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr) = UobsWtAge(stom_ind, 1);
            }
          }
        }
      }
    }


    // 7.5. Calculate other food stomach content
    other_food_diet_prop_weight.setZero();
    for (yr = 0; yr < nyrs; yr++) {                           // Year loop
      for (rsp = 0; rsp < nspp; rsp++) {                      // Predator species loop
        for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
          for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
            other_food_diet_prop_weight(rsp, r_sex, r_age, yr) = Type( 1 );             // Initialize other suitability
            for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                  other_food_diet_prop_weight(rsp, r_sex, r_age, yr) -= diet_prop_weight(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr);
                }
              }
            }
            if(other_food(rsp) > 0){
              other_food_diet_prop_weight(rsp, r_sex, r_age, yr) /= other_food(rsp); // Penalize this
            }
            if(other_food(rsp) == 0){
              other_food_diet_prop_weight(rsp, r_sex, r_age, yr) = 0;
            }
          }
        }
      }
    }

    // START PREDATION
    if (msmMode > 0) {

      // ------------------------------------------------------------------------- //
      // 8.1. SUITABILITY EQUATIONS                                                //
      // ------------------------------------------------------------------------- //

      // 8.1.1. Holsman and MSVPA based suitability // FIXME - not flexible for interannual variation
      if (suitMode == 0) {
        // 8.1.1.1. Calculate stomach proportion over biomass; U/ (W * N)
        suma_suit.setZero();
        for (yr = 0; yr < nyrs; yr++) {                         // Year loop
          for (ksp = 0; ksp < nspp; ksp++) {                    // Prey species loop
            for (rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                  for (r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                    for (k_age = 0; k_age < nages(ksp); k_age++) {  // Prey age loop

                      suit_tmp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 0;

                      if(yr < nyrs_hind){
                        yr_ind = yr;
                      }
                      if(yr >= nyrs_hind){
                        yr_ind = nyrs_hind - 1;
                      }


                      if(AvgN(ksp, k_sex, k_age, yr) > 0){
                        suit_tmp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = diet_prop_weight(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / (AvgN(ksp, k_sex, k_age, yr));

                        // Make into Type 3 MSVPA
                        // U =
                        if(msmMode == 2){
                          suit_tmp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) /= AvgN(ksp, k_sex, k_age, yr);
                        }
                      }


                      // Make sure it is a real number, if not set to 0
                      if(!isFinite(suit_tmp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr))){
                        suit_tmp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = 0;
                      }



                      if (wt(pop_wt_index(ksp), k_sex, k_age, yr_ind ) != 0) {
                        stom_div_bio2(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = suit_tmp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / wt( pop_wt_index(ksp), k_sex, k_age, yr_ind );
                        suma_suit(rsp, r_sex, r_age, yr ) += stom_div_bio2(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Calculate sum of stom_div_bio2 across prey and  prey age for each predator, predator age, and year
                      }
                    }
                  }
                }
              }
            }
          }
        }


        // 8.1.1.2. Calculate suitability
        suit_main.setZero();
        suit_other.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                          // Predator species loop
          for(r_sex = 0; r_sex < nsex(rsp); r_sex++){               // Predator sex
            for (r_age = 0; r_age < nages(rsp); r_age++) {          // Predator age loop
              for (ksp = 0; ksp < nspp; ksp++) {                    // Prey species loop
                for(k_sex = 0; k_sex < nsex(ksp); k_sex++){         // Prey sex loop
                  for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                    for (yr = 0; yr < nyrs_mean; yr++) {            // Suit year loop

                      // Average suitability across years
                      if( suma_suit(rsp, r_sex, r_age, yr ) + other_food_diet_prop_weight(rsp, r_sex, r_age, yr) > 0){
                        suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0) += stom_div_bio2(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / ( suma_suit(rsp, r_sex, r_age, yr ) + other_food_diet_prop_weight(rsp, r_sex, r_age, yr) );
                      }
                    }       // End year loop

                    // FIXME - Add in interannual variability here
                    suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0) /= nyrs_mean;

                    // Remove NAs from crashing populations
                    if(!isFinite(suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0))){
                      suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0) = 0;
                    }

                    // Fill in years
                    for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                      suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0);
                    }

                    // Other suitabilitity
                    for (yr = 0; yr < nyrs; yr++) {
                      suit_other(rsp, r_sex, r_age, yr) += suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // FIXME - include overlap indices
                    }
                  }
                }
              }
            }
          }
        }
      } // End Holsman/MSVPA suitability



      // 8.1.2. Estimate suitability
      if(suitMode > 0){

        // -- Transform predator-prey preference parameters
        // Adopted from https://github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp (Trijoulet et al 2020)
        // Suitability for other food = 1-sum(predator-prey preference)
        // Criteria for predator-prey preference:
        // 1. predator-prey preference > 0 (hence logs)
        // 2. sum(predator-prey preference) + vuln_other = 1
        // 3. 0 <= sum(predator-prey preference) <= 1 (hence logit transformation)
        sum_phi.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          for (ksp = 0; ksp < nspp; ksp++) {              // Prey loop
            sum_phi(rsp) += exp(log_phi(rsp, ksp));
          }
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            vulnerability(rsp,ksp) = exp(log_phi(rsp,ksp))/(1+sum_phi(rsp)); // multinomial logistic transformation
          }
          vulnerability_other(rsp) = vulnerability.row(rsp).sum(); // vulnerability-other=1-sum-vulnerability but transform so sum_vuln+vuln_other=1
        }

        // 8.1.3. GAMMA suitability
        if((suitMode == 1) | (suitMode == 2)){
          Type log_size_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))

          suit_main.setZero();
          for (rsp = 0; rsp < nspp; rsp++) {                           // Pred loop
            for (r_age = 0; r_age < nages(rsp); r_age++) {             // Pred age
              for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                for (ksp = 0; ksp < nspp; ksp++) {                     // Prey loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {     // Prey age
                      for (yr = 0; yr < nyrs; yr++) {                  // Year loop

                        if(yr < nyrs_hind){
                          yr_ind = yr;
                        }
                        if(yr >= nyrs_hind){
                          yr_ind = nyrs_hind - 1;
                        }

                        suit_other(rsp, r_sex, r_age, yr) = vulnerability_other(rsp);

                        switch(suitMode){
                        case 1: // Length-based GAMMA suitability
                          log_size_ratio = log(LbyAge( rsp, r_sex, r_age, yr) / LbyAge( ksp, k_sex, k_age, yr) ); // Log ratio of lengths
                          break;
                        case 2: // Weight-based GAMMA suitability
                          log_size_ratio = log(wt(pop_wt_index(rsp), r_sex, r_age, yr_ind) / wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // Log ratio of weights
                          break;
                        }
                        if(log_size_ratio > 0){
                          // See line 452 from https://github.com/vtrijoulet/Multisp_model_JAE/blob/master/MS_SSM.cpp
                          suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = vulnerability(rsp, ksp) * dgamma(log_size_ratio, gam_a( rsp ), gam_b(rsp)) / dgamma((gam_a(rsp)-1) * gam_b(rsp), gam_a(rsp), gam_b(rsp)); // Scale to 0,1 by dividing by max
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        } // End gamma suitability


        // 8.1.4. Lognormal suitability
        if((suitMode == 3) | (suitMode == 4)){
          Type log_size_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
          suit_main.setZero();
          for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
            for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
              for (r_age = 0; r_age < nages(rsp); r_age++) {                  // Pred age
                for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                      for(yr = 0; yr < nyrs; yr++){

                        if(yr < nyrs_hind){
                          yr_ind = yr;
                        }
                        if(yr >= nyrs_hind){
                          yr_ind = nyrs_hind - 1;
                        }

                        suit_other(rsp, r_sex, r_age, yr) = vulnerability_other(rsp);

                        switch(suitMode){
                        case 3: // Length-based lognormal suitability
                          log_size_ratio = log(LbyAge( rsp, r_sex, r_age, yr) / LbyAge( ksp, k_sex, k_age, yr) ); // Log ratio of lengths
                          break;
                        case 4: // Weight-based lognormal suitability
                          log_size_ratio = log(wt(pop_wt_index(rsp), r_sex, r_age, yr_ind) / wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // Log ratio of weights
                          break;
                        }
                        // if prey are smaller than predator:
                        if(log_size_ratio > 0){
                          suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = vulnerability(rsp, ksp) * dnorm(log_size_ratio, gam_a(rsp), gam_b(rsp)) / dnorm(gam_a(rsp), gam_a(rsp), gam_b(rsp)); // Divide by mode to scale to 1
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        } // End lognormal selectivity
      } // End suitability estimation


      // ------------------------------------------------------------------------- //
      // 8. PREDATION MORTALITY EQUATIONS                                          //
      // ------------------------------------------------------------------------- //
      // -- 8.1. HOLSMAN PREDATION MORTALITY
      if ((msmMode == 1) | (msmMode == 2)) {

        // 8.1.3. Calculate available food
        avail_food.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                        // Predator species loop
          for(r_sex = 0; r_sex < nsex(rsp); r_sex++){             // Predator sex loop
            for (r_age = 0; r_age < nages(rsp); r_age++) {          // Predator age loop
              for (yr = 0; yr < nyrs; yr++) {                       // Year loop

                if(yr < nyrs_hind){
                  yr_ind = yr;
                }
                if(yr >= nyrs_hind){
                  yr_ind = nyrs_hind - 1;
                }

                for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                      avail_food(rsp, r_sex, r_age, yr) += suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) ; // FIXME - include overlap indices: FIXME - mn_wt_stom?
                    }
                  }
                }
                // Other food
                avail_food(rsp, r_sex, r_age, yr) += other_food(rsp) * (Type(1) - suit_other(rsp, r_sex, r_age, yr));
              }
            }
          }
        }


        // 8.1.3. Calculate predation mortality, diet proportion, and biomass easten
        M2.setZero();
        B_eaten_as_prey.setZero();
        diet_prop_weight_ave_hat.setZero();
        for (ksp = 0; ksp < nspp; ksp++) {                        // Prey species loop
          for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
            for (k_age = 0; k_age < nages(ksp); k_age++) {          // Prey age loop
              for (rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop
                for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                  for (r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                    for (yr = 0; yr < nyrs; yr++) {                       // Year loop

                      if(yr < nyrs_hind){
                        yr_ind = yr;
                      }
                      if(yr >= nyrs_hind){
                        yr_ind = nyrs_hind - 1;
                      }

                      if(avail_food(rsp, r_sex, r_age, yr) > 0){

                        switch(msmMode){
                        case 1: // Type 2 MSVPA
                          M2(ksp, k_sex, k_age, yr) += (AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr); // #FIXME - include indices of overlap
                          M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = (AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr);
                          B_eaten_as_prey(ksp, k_sex, k_age, yr) += AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)/ avail_food(rsp, r_sex, r_age, yr);
                          B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                          diet_prop_weight_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)/ avail_food(rsp, r_sex, r_age, yr);
                          break;


                        case 2: // Type 3 MSVPA
                          M2(ksp, k_sex, k_age, yr) += pow(AvgN(ksp, k_sex, k_age, yr) , msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr)
                                                                                                                                                     * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / (AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // #FIXME - include indices of overlap
                          M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / (AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind));
                          B_eaten_as_prey(ksp, k_sex, k_age, yr) += pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                          B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) *  AvgN(rsp, r_sex, r_age, yr) * ration(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                          diet_prop_weight_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pow(AvgN(ksp, k_sex, k_age, yr) , msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr) / (AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind));
                          break;
                        }
                      }
                    }
                    // Average diet proportion
                    for (yr = 0; yr < nyrs_hind; yr++) {                       // Year loop
                      diet_prop_weight_ave_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age) += diet_prop_weight_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)/nyrs_hind;
                    }
                  }
                }
              }
            }
          }
        }
      } // End 8..1. Holsman predation mortality


      // 8.2. KINZEY PREDATION EQUATIONS
      if (msmMode > 2) {

        // 8.2.3. Initialize counters
        Type Pred_ratio = 0;          // Predator ratio
        Type Prey_ratio = 0;          // Prey ratio
        Type NS_Z = 0;                // N(k,yr,a) * survival/Z = 0;
        Type Tmort = 0;               // Mortality on other
        Type Q_ksum_l = 0;            // Diet sum
        Type Term = 0;                // Linear adjustment for predation


        // 8.2.4. Calculate equilibrium N predators and prey in styr_pred for each species X age: FIXME: May want to have this be the final year of a projection!
        N_pred_eq.setZero();
        N_prey_eq.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {
          for (ksp = 0; ksp < nspp; ksp++) {
            for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for (r_age = 0; r_age < nages(rsp); r_age++) {
                  for (k_age = 0; k_age < nages(ksp); k_age++) {
                    N_pred_eq(rsp, r_sex, r_age) += NByage(rsp, r_sex, r_age, 0) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0); // Denominator of Eq. 17 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
                    N_prey_eq(ksp, k_sex, k_age) += NByage(ksp, k_sex, k_age, 0) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, 0); // Denominator of Eq. 16 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
                  }
                }
              }
            }
          }
        }


        // 8.2.5. Calculate available prey and predator for each year
        N_pred_yrs.setZero();
        N_prey_yrs.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp < nspp; rsp++) {
            for (ksp = 0; ksp < nspp; ksp++) {
              for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                  for (r_age = 0; r_age < nages(rsp); r_age++) {
                    for (k_age = 0; k_age < nages(ksp); k_age++) {
                      N_pred_yrs(rsp, r_sex, r_age, yr) += NByage(rsp, r_sex, r_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Numerator of Eq. 17 Kinzey and Punt (2009) 1st year // FIXME: Use averageN?
                      N_prey_yrs(ksp, k_sex, k_age, yr) += NByage(ksp, k_sex, k_age, yr) * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Numerator of Eq. 16 Kinzey and Punt (2009) 1st year
                    }
                  }
                }
              }
            }
          }
        }


        // 8.2.6. Calculate predator functional response (Table 1 Kinzey and Punt (2009))
        for (rsp = 0; rsp < nspp; rsp++) {                          // Predator loop
          for (ksp = 0; ksp < (nspp + 1); ksp++) {                  // Prey loop
            Term = 1.0e-10 + H_1(rsp, ksp) * (Type(1) + H_1a(rsp) * H_1b(rsp) / (Type(r_age) + H_1b(rsp) + Type(1.0e-10))); // Eq. 15 Kinzey and Punt (2009)
            for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
              for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
                for (yr = 0; yr < nyrs; yr++) {                         // Year loop

                  // Observed species
                  if(ksp < nspp){
                    for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                      for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop

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
        M2.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
              for (ksp = 0; ksp  <  nspp; ksp++) {
                for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                  for (r_age = 0; r_age < nages(rsp); r_age++) {
                    for (k_age = 0; k_age < nages(ksp); k_age++) {
                      M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = pred_resp(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)
                      * suit_main(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * NByage(rsp, r_sex, r_age, yr);
                      M2(ksp, k_sex, k_age, yr) += M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr);
                    }
                  }
                }
              }
            }
          }
        }

        // 8.2.9. Numbers and mass eaten (of modeled prey species and "other prey"); Equations 8 and 9 from Kinzey and Punt 2009
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
              for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                for (ksp = 0; ksp  <  nspp+1; ksp++) {
                  // Species included
                  if (ksp < nspp) {
                    for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                      for (k_age = 0; k_age < nages(ksp); k_age++) {
                        // Year indices
                        if(yr < nyrs_hind){
                          yr_ind = yr;
                        }
                        if(yr >= nyrs_hind){
                          yr_ind = nyrs_hind - 1;
                        }
                        N_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * AvgN(ksp, k_sex, k_age, yr); // Eq. 8 Kinzey and Punt (2009)
                        B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = M2_prop(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) * AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind);
                        B_eaten_as_pred(rsp, r_sex, r_age, yr) += B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr);
                        B_eaten_as_prey(ksp, k_sex, k_age, yr) += B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr); // Results by all species: Eq. 11 Kinzey and Punt (2009)
                      }
                    }
                    // Other food
                  }else{
                    B_eaten_as_pred(rsp, r_sex, r_age, yr)  += other_food(rsp) * (Type(1) - exp(-pred_resp(rsp + (nspp * r_sex), nspp*2, r_age, 0, yr)  * NByage(rsp, r_sex, r_age, yr)));      // Eq. 11 Kinzey and Punt (2009)
                  }
                }
              }
            }
          }
        }


        // 8.2.11. Predicted diet proportion as weight of prey-at-age in diet of predator-at-age (Eqn 15)
        // NOTE: Only including prey species in the model. Does not include other food
        diet_prop_weight_ave_hat.setZero();
        for (yr = 0; yr < nyrs_hind; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
              for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                for (ksp = 0; ksp  <  nspp; ksp++) {
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {
                      diet_prop_weight_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr) = B_eaten(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)/B_eaten_as_pred(rsp, r_sex, r_age, yr);
                      diet_prop_weight_ave_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age) += diet_prop_weight_hat(rsp + (nspp * r_sex), ksp + (nspp * k_sex), r_age, k_age, yr)/nyrs_hind; // FIXME: Only use years from the diet data
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
        for (rsp = 0; rsp < nspp; rsp++) {
          for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
            for (r_age = 0; r_age < nages(rsp); r_age++) {
              numer = 0;
              denom = 0;
              // Calculate year-specific values
              for (yr = 0; yr < nyrs; yr++) {
                if(AvgN(rsp, r_sex, r_age, yr) > 0){
                  ration_hat(rsp, r_sex, r_age, yr) = B_eaten_as_pred(rsp, r_sex, r_age, yr)/AvgN(rsp, r_sex, r_age, yr); // NOTE: Divide by 365 to make into daily ration
                }
              }

              // Find average over hindcast// FIXME: suit_years?
              for (yr = 0; yr < nyrs_hind; yr++) {
                numer += B_eaten_as_pred(rsp, r_sex, r_age, yr); // NOTE: Divide by 365 to make into daily ration
                denom += AvgN(rsp, r_sex, r_age, yr);
              }
              ration_hat_ave(rsp, r_sex, r_age) = numer / denom;
            }
          }
        }
        // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
      } // End 8.2. Kinzey predation
      // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
    } // End 8. Predation mortality
    // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //


    // ------------------------------------------------------------------------- //
    // 9. INDEX COMPONENTS EQUATIONS                                             //
    // ------------------------------------------------------------------------- //

    // -- 9.1. Index of abundance/biomass
    for(srv_ind = 0; srv_ind < srv_biom_ctl.rows(); srv_ind++){

      srv = srv_biom_ctl(srv_ind, 0) - 1;            // Temporary survey index
      sp = srv_biom_ctl(srv_ind, 1) - 1;             // Temporary index of species
      flt_yr = srv_biom_ctl(srv_ind, 2);              // Temporary index for years of data

      mo = srv_biom_n(srv_ind, 0);                    // Temporary index for month

      srv_bio_hat(srv_ind) = 0;                      // Initialize

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

      for (age = 0; age < nages(sp); age++) {
        for(sex = 0; sex < nsex(sp); sex++){
          // Weight
          if(flt_units(srv) == 1){
            srv_bio_hat(srv_ind) += NByage(sp, sex, age, flt_yr) * exp( - (mo/12) * Zed(sp, sex, age, flt_yr)) * sel(srv, sex, age, yr_ind) * wt( flt_wt_index(srv), sex, age, yr_ind );
          }
          // Numbers
          if(flt_units(srv) == 2){
            srv_bio_hat(srv_ind) += NByage(sp, sex, age, flt_yr) * exp( - (mo/12) * Zed(sp, sex, age, flt_yr)) * sel(srv, sex, age, yr_ind);
          }
        }
      }
    }


    // -- 9.2. Analytical survey q following Ludwig and Martell 1994
    srv_n_obs.setZero();
    srv_q_analytical.setZero();
    for(srv_ind = 0; srv_ind < srv_biom_ctl.rows(); srv_ind++){



      srv = srv_biom_ctl(srv_ind, 0) - 1;            // Temporary survey index
      sp = srv_biom_ctl(srv_ind, 1) - 1;             // Temporary index of species
      flt_yr = srv_biom_ctl(srv_ind, 2);             // Temporary index for years of data

      if(flt_yr > 0){
        flt_yr = flt_yr - styr;


        mo = srv_biom_n(srv_ind, 0);                    // Temporary index for month
        if(flt_yr < nyrs_hind){

          // If etimated standard deviation or analytical sigma (non-time-varying)
          if (est_sigma_srv(srv) > 0) {
            srv_n_obs(srv) += 1; // Add one if survey is used
            srv_q_analytical(srv) += log(srv_biom_obs(srv_ind, 0) / srv_bio_hat(srv_ind));
          }

          // If time-varying sigma
          if (est_sigma_srv(srv) == 0 ) {
            srv_n_obs(srv) += 1 / square(srv_biom_obs(srv_ind, 1));
            srv_q_analytical(srv) += log(srv_biom_obs(srv_ind, 0) / srv_bio_hat(srv_ind)) / square(srv_biom_obs(srv_ind, 1));
          }
        }
      }
    }

    // Take average
    for(srv = 0 ; srv < n_flt; srv ++){
      srv_q_analytical(srv) = exp(srv_q_analytical(srv) / srv_n_obs(srv));

      // Set srv_q to analytical if used
      if(est_srv_q(srv) == 3){
        for(yr = 0; yr < nyrs_hind; yr++){
          srv_q(srv, yr) = srv_q_analytical(srv);
        }
      }
    }


    // -- 9.3. Survey Biomass - multiply by q
    for(srv_ind = 0; srv_ind < srv_biom_ctl.rows(); srv_ind++){

      srv = srv_biom_ctl(srv_ind, 0) - 1;            // Temporary survey index
      flt_yr = srv_biom_ctl(srv_ind, 2);      // Temporary index for years of data

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

      srv_bio_hat(srv_ind) = srv_q(srv, yr_ind) * srv_bio_hat(srv_ind); // pow(srv_bio_hat(srv_ind), (1 + srv_q_pow(srv)));

    }



    // -- 9.4. Calculate analytical sigma following Ludwig and Walters 1994
    srv_n_obs.setZero();
    sigma_srv_analytical.setZero();
    for(srv_ind = 0; srv_ind < srv_biom_ctl.rows(); srv_ind++){

      srv = srv_biom_ctl(srv_ind, 0) - 1;            // Temporary survey index
      flt_yr = srv_biom_ctl(srv_ind, 2);      // Temporary index for years of data

      if(flt_yr > 0){
        flt_yr = flt_yr - styr;

        if(flt_yr < nyrs_hind){
          srv_n_obs(srv) += 1; // Add one if survey is used
          sigma_srv_analytical(srv) += square( log(srv_biom_obs(srv_ind, 0)) - log(srv_bio_hat(srv_ind)));
        }
      }
    }

    for(srv = 0 ; srv < n_flt; srv ++){
      sigma_srv_analytical(srv) = sqrt(sigma_srv_analytical(srv) / srv_n_obs(srv));
    }




    // ------------------------------------------------------------------------- //
    // 10. FISHERY COMPONENTS EQUATIONS                                          //
    // ------------------------------------------------------------------------- //
    // 10.5. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
    for(fsh_ind = 0; fsh_ind < fsh_biom_ctl.rows(); fsh_ind++){

      flt = fsh_biom_ctl(fsh_ind, 0) - 1;            // Temporary fishery index
      sp = fsh_biom_ctl(fsh_ind, 1) - 1;             // Temporary index of species
      flt_yr = fsh_biom_ctl(fsh_ind, 2);      // Temporary index for years of data
      mo = fsh_biom_n(fsh_ind, 0);                   // Temporary index for month

      fsh_bio_hat(fsh_ind) = 0;                      // Initialize

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
        for (age = 0; age < nages(sp); age++) {
          // FIXME: add wt_yr_ind for MSE

          // By weight
          if(flt_units(flt) == 1){
            fsh_bio_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Zed(sp, sex, age, flt_yr) * (1 - exp(-Zed(sp, sex, age, flt_yr))) * NByage(sp, sex, age, flt_yr) * wt( flt_wt_index(flt), sex, age, yr_ind ); // 5.5.
          }

          // By numbers
          if(flt_units(flt) == 2){
            fsh_bio_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Zed(sp, sex, age, flt_yr) * (1 - exp(-Zed(sp, sex, age, flt_yr))) * NByage(sp, sex, age, flt_yr);
          }
        }
      }
    }


    // ------------------------------------------------------------------------- //
    // 11. COMPOSITION EQUATIONS                                                  //
    // ------------------------------------------------------------------------- //

    // -- 10.3. Composition
    age_obs_hat.setZero();
    comp_hat.setZero();
    age_hat.setZero();
    for(comp_ind = 0; comp_ind < comp_hat.rows(); comp_ind++){

      flt = comp_ctl(comp_ind, 0) - 1;            // Temporary fishery index
      sp = comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
      flt_sex = comp_ctl(comp_ind, 2);            // Temporary index for comp sex (0 = combined, 1 = female, 2 = male)
      comp_type = comp_ctl(comp_ind, 3);          // Temporary index for comp type (0 = age, 1 = length)
      yr = comp_ctl(comp_ind, 4);                 // Temporary index for years of data
      mo = comp_n(comp_ind, 0);                   // Temporary index for month

      n_hat(comp_ind) = 0;                          // Initialize

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
      for (age = 0; age < nages(sp); age++) {


        switch(flt_type(flt)){
        case 1: // - Fishery
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              age_hat(comp_ind, age ) += F_flt_age(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr); // 5.4.
              n_hat( comp_ind ) += F_flt_age(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            age_hat(comp_ind, age ) = F_flt_age(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
            n_hat( comp_ind ) += F_flt_age(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = F_flt_age(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
              // Total numbers
              n_hat(comp_ind) += F_flt_age(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
            }
            break;
          }
          break;


        case 2: // - Survey
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age ) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
              // Total numbers
              n_hat(comp_ind) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));   // Total numbers
              // FIXME - put power in for catchability
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            // Survey catch-at-age
            age_hat(comp_ind, age ) = NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
            // Total numbers
            n_hat(comp_ind) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
              // Total numbers
              n_hat(comp_ind) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
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

      // Get true age comp
      for (age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
        if(n_hat(comp_ind) > 0){ //FIXME - not differentiable
          true_age_comp_hat(comp_ind, age ) = age_hat(comp_ind, age ) / n_hat(comp_ind);
        }
      }


      // Adjust for aging error
      for (int obs_age = 0; obs_age < nages(sp); obs_age++) {
        for (int true_age = 0; true_age < nages(sp); true_age++) {
          age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(sp, true_age, obs_age);
        }
      }

      // Adjust for aging error for joint data
      if(flt_sex == 3){
        for (int obs_age = nages(sp); obs_age < nages(sp) * 2; obs_age++) {
          for (int true_age = nages(sp); true_age < nages(sp) * 2; true_age++) {

            // Adjust indexing for joint age/length comp
            int true_age_tmp = true_age - nages(sp);
            int obs_age_tmp = obs_age - nages(sp);

            age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(sp, true_age_tmp, obs_age_tmp);
          }
        }
      }


      //  Survey catch-at-age - standardize to sum to 1
      if (comp_type == 0) {
        for (age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
          if(n_hat(comp_ind) > 0){ //FIXME - not differentiable
            comp_hat(comp_ind, age) = age_obs_hat(comp_ind, age ) / n_hat(comp_ind);
          }
        }
      }


      // Convert from catch-at-age to catch-at-length
      if ( comp_type == 1) {
        for (ln = 0; ln < nlengths(sp); ln++) {
          for (age = 0; age < nages(sp); age++) {

            sex = 0;

            // Adjust sex for males/females
            if((flt_sex > 0) & (flt_sex < 3)){
              sex = flt_sex - 1;
            }

            comp_hat(comp_ind, ln ) += age_obs_hat(comp_ind, age) * age_trans_matrix(flt_age_transition_index(flt), sex, age, ln );
          }
        }


        // Convert from catch-at-age to catch-at-length for joint comp data
        if ( flt_sex == 3) {
          for (ln = nlengths(sp); ln < nlengths(sp) * 2; ln++) {
            for (age = nages(sp); age < nages(sp) * 2; age++) {

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
        for (ln = 0; ln < nlengths(sp) * joint_adjust(comp_ind); ln++) {
          if(n_hat(comp_ind) > 0){ //FIXME - not differentiable
            comp_hat(comp_ind, ln ) = comp_hat(comp_ind, ln) / n_hat(comp_ind);
          }
        }
      }
    }


    // ------------------------------------------------------------------------- //
    // 12. PREDICTED STOMACH CONTENT                                             //
    // ------------------------------------------------------------------------- //

    // Predict stomach content
    // 7.4. Reorganize UobsWTAge content
    for(int stom_ind = 0; stom_ind < UobsWtAge.rows(); stom_ind++){
      rsp = UobsWtAge_ctl(stom_ind, 0) - 1; // Index of pred
      ksp = UobsWtAge_ctl(stom_ind, 1) - 1; // Index of prey
      r_sex = UobsWtAge_ctl(stom_ind, 2); // Index of pred sex
      k_sex = UobsWtAge_ctl(stom_ind, 3); // Index of prey sex
      r_age = UobsWtAge_ctl(stom_ind, 4) - minage(rsp); // Index of pred age
      k_age = UobsWtAge_ctl(stom_ind, 5) - minage(ksp); // Index of prey age
      flt_yr = UobsWtAge_ctl(stom_ind, 6); // Index of year

      // Predator
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

        if(k_sex > 0){
          k_sexes(stom_ind, 0) = k_sex - 1;  k_sexes(stom_ind, 1) = k_sex - 1;
        }
      }

      // Initialize
      UobsWtAge_hat(stom_ind, 1) = 0;

      for(int j = 0; j < 2; j ++){
        for(int k = 0; k < 2; k ++){

          if(flt_yr > 0){
            yr = flt_yr - styr;

            if(yr < nyrs_hind){
              UobsWtAge_hat(stom_ind, 1) += diet_prop_weight_hat(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age, yr)/4;; //FIXME: take weighted average
            }
          }

          // Average of years
          if(flt_yr == 0){
            for (yr = 0; yr < nyrs_hind; yr++) {
              UobsWtAge_hat(stom_ind, 1) += diet_prop_weight_ave_hat(rsp + (nspp * r_sexes(stom_ind, j)), ksp + (nspp * k_sexes(stom_ind, k)), r_age, k_age)/4; //FIXME: take weighted average
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
  depletion.setZero();
  depletionSSB.setZero();
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs; yr++){
      if(DynamicHCR == 0){
        depletion(sp, yr) = biomass(sp, yr)/B0(sp);
        depletionSSB(sp, yr) = biomassSSB(sp, yr)/SB0(sp);
      }

      if(DynamicHCR == 1){
        depletion(sp, yr) = biomass(sp, yr)/DynamicB0(sp, yr);
        depletionSSB(sp, yr) = biomassSSB(sp, yr)/DynamicSB0(sp, yr);
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 14. LIKELIHOOD EQUATIONS
  // ------------------------------------------------------------------------- //

  // 14.0. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(18, n_flt); jnll_comp.setZero();  // matrix of negative log-likelihood components

  // -- Data likelihood components
  // Slot 0 -- Survey biomass
  // Slot 1 -- Total catch (kg)
  // Slot 2 -- Age/length composition
  // Slot 3 -- Sex ratio likelihood
  // Slot 4 -- Selectivity
  // Slot 5 -- Selectivity annual deviates
  // Slot 6 -- Survey selectivity normalization
  // Slot 7 -- Survey catchability prior
  // Slot 8 -- Survey catchability annual deviates
  // -- Priors/penalties
  // Slot 9 -- Alpha prior for stock recruit relationship
  // Slot 10 -- Tau -- Annual recruitment deviation
  // Slot 11 -- init_dev -- Initial abundance-at-age
  // Slot 12 -- Epsilon -- Annual fishing mortality deviation
  // Slot 13 -- SPR penalities
  // Slot 14 -- N-at-age < 0 penalty
  // -- M2 likelihood components
  // Slot 15 -- Ration likelihood
  // Slot 16 -- Ration penalties
  // Slot 17 -- Diet proportion by weight likelihood


  // 14.1. FIT OBJECTIVE FUNCTION
  // Slot 0 -- Survey biomass -- NFMS annual BT survey and EIT survey
  Type srv_std_dev = 0;
  for (srv_ind = 0; srv_ind < srv_biom_obs.rows(); srv_ind++) {

    srv = srv_biom_ctl(srv_ind, 0) - 1;            // Temporary survey index
    sp = srv_biom_ctl(srv_ind, 1) - 1;             // Species is the 2nd column
    flt_yr = srv_biom_ctl(srv_ind, 2);             // Temporary index for years of data

    if(flt_yr > 0){
      // Set up variance
      switch (est_sigma_srv(srv)) {
      case 0:     // Provided standard deviation
        srv_std_dev = srv_biom_obs(srv_ind, 1);
        break;
      case 1:     // Estimated standard deviation
        srv_std_dev = sigma_srv_index(srv);
        break;
      case 2:     // Analytical
        srv_std_dev = sigma_srv_analytical(srv);
        sigma_srv_index(srv) = sigma_srv_analytical(srv);
        break;
      default:
        error("Invalid 'Estimate_sigma_index'");
      }

      srv_log_sd_hat(srv_ind) = srv_std_dev;


      // Only include years from hindcast
      if(flt_type(srv) > 0){
        if(flt_yr <= endyr){
          if(srv_bio_hat(srv_ind) > 0){
            jnll_comp(0, srv) -= dnorm(log(srv_biom_obs(srv_ind, 0)) - square(srv_std_dev)/2, log(srv_bio_hat(srv_ind)), srv_std_dev, true);  // pow(log(srv_biom_obs(srv_ind, 0)) - log(srv_bio_hat(srv_ind)), 2) / (2 * square( srv_std_dev )); // NOTE: This is not quite the lognormal and biohat will be the median.

            // Martin's
            // jnll_comp(0, srv)+= 0.5*square((log(srv_biom_obs(srv_ind, 0))-log(srv_bio_hat(srv_ind))+square(srv_std_dev)/2.0)/srv_std_dev);


            SIMULATE {
              srv_biom_obs(srv_ind, 0) = rnorm(log(srv_bio_hat(srv_ind)), srv_std_dev);  // Simulate response
              srv_biom_obs(srv_ind, 0) = exp(srv_biom_obs(srv_ind, 0));
            }
          }
        }
      }
    }
  }


  // Slot 1 -- Total catch -- Fishery observer dat
  Type fsh_std_dev = 0;
  for (fsh_ind = 0; fsh_ind < fsh_biom_obs.rows(); fsh_ind++) {

    flt = fsh_biom_ctl(fsh_ind, 0) - 1;            // Temporary fishery index
    sp = fsh_biom_ctl(fsh_ind, 1) - 1;             // Species is the column 3
    flt_yr = fsh_biom_ctl(fsh_ind, 2);             // Temporary index for years of data
    yr = flt_yr - styr;                            // Temporary index of years. Start at 0.

    if(flt_yr > 0){

      // Set up variance
      switch (est_sigma_fsh(flt)) {
      case 0: // Provided standard deviation
        fsh_std_dev = fsh_biom_obs(fsh_ind, 1);
        break;
      case 1:     // Estimated standard deviation
        fsh_std_dev = sigma_fsh_catch(flt);
        break;
      default:
        error("Invalid 'Estimate_sigma_catch'");
      }

      fsh_log_sd_hat(fsh_ind) = fsh_std_dev; // Save estimated log_sd

      // Add only years from hindcast
      if(flt_type(flt) == 1){
        if(flt_yr <= endyr){
          if(fsh_biom_obs(fsh_ind, 0) > 0){
            jnll_comp(1, flt) -= dnorm(log(fsh_biom_obs(fsh_ind, 0)) - square(fsh_std_dev)/2, log(fsh_bio_hat(fsh_ind)), fsh_std_dev, true) ;

            // Martin's
            // jnll_comp(1, flt)+= 0.5*square((log(fsh_biom_obs(fsh_ind, 0))-log(fsh_bio_hat(fsh_ind)))/fsh_std_dev);


            SIMULATE {
              fsh_biom_obs(fsh_ind, 0) = rnorm(log(fsh_bio_hat(fsh_ind)), fsh_std_dev);  // Simulate response
              fsh_biom_obs(fsh_ind, 0) = exp(fsh_biom_obs(fsh_ind, 0));
            }


            // Slot 12 -- Epsilon -- Annual fishing mortality deviation
            jnll_comp(12, flt) += square(F_dev(flt, yr));      // Fishing mortality deviation using penalized likelihood.
          }
        }
      }
    }
  }



  // Slot 2 -- Age composition
  for (comp_ind = 0; comp_ind < comp_obs.rows(); comp_ind++) {

    flt = comp_ctl(comp_ind, 0) - 1;        // Temporary survey index
    sp = comp_ctl(comp_ind, 1) - 1;         // Temporary index of species
    flt_sex = comp_ctl(comp_ind, 2);        // Temporary index for comp sex (0 = combined, 1 = female, 2 = male, 3 = joint)
    comp_type = comp_ctl(comp_ind, 3);      // Temporary index for comp type (0 = age, 1 = length)
    yr = comp_ctl(comp_ind, 4);             // Temporary index for years of data

    if(yr > 0){

      // Adjustment for joint sex composition data
      joint_adjust(comp_ind) = 1;
      if(flt_sex == 3){
        joint_adjust(comp_ind) = 2;
      }

      // Number of ages/lengths
      Type n_comp = 0;
      if(comp_type == 0){
        n_comp = nages(sp) * joint_adjust(comp_ind); // For sex = 3 (joint sex comp data)

        // Accumulation age
        for(age = 0; age < n_comp; age++){
          int age_test = age;
          if(age >= nages(sp)){
            age_test = age - nages(sp);
          }

          // Lower
          if(age_test < flt_accum_age_lower(flt)){
            comp_obs(comp_ind, flt_accum_age_lower(flt) + (nages(sp) * (joint_adjust(comp_ind)-1))) += comp_obs(comp_ind, age);
            comp_obs(comp_ind, age) = 0;

            comp_hat(comp_ind, flt_accum_age_lower(flt) + (nages(sp) * (joint_adjust(comp_ind)-1))) += comp_hat(comp_ind, age);
            comp_hat(comp_ind, age) = 0;
          }
          // Upper
          if(age_test > flt_accum_age_upper(flt)){
            comp_obs(comp_ind, flt_accum_age_upper(flt) + (nages(sp) * (joint_adjust(comp_ind)-1))) += comp_obs(comp_ind, age);
            comp_obs(comp_ind, age) = 0;

            comp_hat(comp_ind, flt_accum_age_upper(flt) + (nages(sp) * (joint_adjust(comp_ind)-1))) += comp_hat(comp_ind, age);
            comp_hat(comp_ind, age) = 0;
          }
        }
      }
      if(comp_type == 1){
        n_comp = nlengths(sp) * joint_adjust(comp_ind);
      }

      // Only use years wanted
      if(flt_type(flt) > 0){
        if(yr <= endyr){
          for (ln = 0; ln < n_comp; ln++) {
            if(!isNA( comp_obs(comp_ind, ln) )){
              if(comp_hat(comp_ind, ln) > 0){
                // jnll_comp(2, flt) -= comp_weights(flt) * Type(comp_n(comp_ind, 1)) * (comp_obs(comp_ind, ln) + 0.00001) * log(comp_hat(comp_ind, ln) + 0.00001) ;

                // Martin's
                jnll_comp(2, flt) -= comp_weights(flt) * Type(comp_n(comp_ind, 1)) * (comp_obs(comp_ind, ln) + 0.00001) * log((comp_hat(comp_ind, ln)+0.00001) / (comp_obs(comp_ind, ln) + 0.00001)) ;

              }
            }
          }
        }
      }
    }
  }


  // Slot 3 -- Observed sex ratio
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs_hind; yr++){
      if((nsex(sp) == 2) & (estDynamics(sp) == 0)){
        if(est_sex_ratio(sp) == 1){
          jnll_comp(3, sp) -= dnorm(sex_ratio_mean_hat(sp, yr), sex_ratio(sp, 1), sex_ratio_sigma(sp)); // Using the 2nd age here, cause recruitment is assumed to be age 1 (c++ 0)
        }

        if(est_sex_ratio(sp) == 2){
          for(age = 1; age < nages(sp); age++){ // Start at age 2 because age 1 is fixed
            jnll_comp(3, sp) -= dnorm(sex_ratio_hat(sp, age, yr), sex_ratio(sp, age), sex_ratio_sigma(sp));
          }
        }
      }
    }
  }



  // Slot 4-6 -- Selectivity
  matrix<Type> sel_slp_dev_ll(n_flt, nyrs_hind);sel_slp_dev_ll.setZero();
  matrix<Type> sel_inf_dev_ll(n_flt, nyrs_hind);sel_inf_dev_ll.setZero();
  for(flt = 0; flt < n_flt; flt++){ // Loop around surveys
    jnll_comp(4, flt) = 0;
    sp = flt_spp(flt);

    if(flt_type(flt) > 0){
      if (flt_sel_type(flt) == 2) {

        for(sex = 0; sex < nsex(sp); sex++){
          for (age = 0; age < (nages(sp) - 1); age++) {
            if ( sel(flt, sex, age, 0) > sel(flt, sex, age + 1, 0)) {
              jnll_comp(4, flt) += sel_curve_pen(flt, 0) * pow( log(sel(flt, sex, age, 0) / sel(flt, sex, age + 1, 0) ), 2);
            }
          }
        }

        for(sex = 0; sex < nsex(sp); sex++){
          // Extract only the selectivities we want
          vector<Type> sel_tmp(nages(sp)); sel_tmp.setZero();

          for (age = 0; age < nages(sp); age++) {
            sel_tmp(age) = log(sel(flt, sex, age, 0));
          }

          for (age = 0; age < nages(sp) - 2; age++) {
            sel_tmp(age) = first_difference( first_difference( sel_tmp ) )(age);
            jnll_comp(4, flt) += sel_curve_pen(flt, 1) * pow( sel_tmp(age) , 2);
          }
        }
      }
    }


    // Penalized/random effect likelihood time-varying logistic/double-logistic selectivity deviates
    if(((flt_varying_sel(flt) == 1)|(flt_varying_sel(flt) == 2)) & (flt_sel_type(flt) != 2) & (flt_sel_type(flt) != 5) & (flt_type(flt) > 0)){
      for(sex = 0; sex < nsex(sp); sex ++){
        for(yr = 0; yr < nyrs_hind; yr++){

          jnll_comp(5, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr), Type(0.0), sigma_sel(flt), true);
          jnll_comp(5, flt) -= dnorm(ln_sel_slp_dev(0, flt, sex, yr), Type(0.0), sigma_sel(flt), true);

          // Double logistic deviates
          if(flt_sel_type(flt) == 3){
            jnll_comp(5, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr), Type(0.0), sigma_sel(flt), true);
            jnll_comp(5, flt) -= dnorm(ln_sel_slp_dev(1, flt, sex, yr), Type(0.0), sigma_sel(flt), true);
          }
        }
      }
    }

    // Penalized/random effect likelihood time-varying non-parametric (Taylor et al 2014) selectivity deviates
    if(((flt_varying_sel(flt) == 1) | (flt_varying_sel(flt) == 2)) & (flt_sel_type(flt) == 5) & (flt_type(flt) > 0)){
      for(age = 0; age < nages(sp); age++){ //NOTE: extends beyond selectivity age range, but should be mapped to 0 in map function
        for(sex = 0; sex < nsex(sp); sex++){
          for(yr = 0; yr < nyrs_hind; yr++){
            jnll_comp(5, flt) -= dnorm(sel_coff_dev(flt, sex, age, yr), Type(0.0), sigma_sel(flt), true);
          }
        }
      }
    }


    // Random walk: Type 4 = random walk on ascending and descending for double logistic; Type 5 = ascending only for double logistics
    if(((flt_varying_sel(flt) == 4)|(flt_varying_sel(flt) == 5)) & (flt_sel_type(flt) != 2) & (flt_sel_type(flt) != 5) & (flt_type(flt) > 0)){
      for(sex = 0; sex < nsex(sp); sex ++){
        for(yr = 1; yr < nyrs_hind; yr++){ // Start at second year

          // Logistic deviates
          sel_inf_dev_ll(flt, yr) = -dnorm(sel_inf_dev(0, flt, sex, yr) - sel_inf_dev(0, flt, sex, yr-1), Type(0.0), sigma_sel(flt), true);
          sel_slp_dev_ll(flt, yr) = -dnorm(ln_sel_slp_dev(0, flt, sex, yr) - ln_sel_slp_dev(0, flt, sex, yr-1), Type(0.0), sigma_sel(flt), true);

          // Logistic deviates
          jnll_comp(5, flt) -= dnorm(ln_sel_slp_dev(0, flt, sex, yr) - ln_sel_slp_dev(0, flt, sex, yr-1), Type(0.0), sigma_sel(flt), true);
          jnll_comp(5, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr) - sel_inf_dev(0, flt, sex, yr-1), Type(0.0), 4 * sigma_sel(flt), true);

          // Double logistic deviates
          if((flt_sel_type(flt) == 3) & (flt_varying_sel(flt) == 4)){
            jnll_comp(5, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr) - sel_inf_dev(1, flt, sex, yr-1), Type(0.0), sigma_sel(flt), true);
            jnll_comp(5, flt) -= dnorm(ln_sel_slp_dev(1, flt, sex, yr) - ln_sel_slp_dev(1, flt, sex, yr-1), Type(0.0), sigma_sel(flt), true);
          }
        }

        // Sum to zero constraint
        Type sel_slp_dev1_sum = 0;
        Type sel_inf_dev1_sum = 0;
        Type sel_slp_dev2_sum = 0;
        Type sel_inf_dev2_sum = 0;

        for(yr = 0; yr < nyrs_hind; yr++){
          sel_slp_dev1_sum += ln_sel_slp_dev(0, flt, sex, yr);
          sel_inf_dev1_sum += sel_inf_dev(0, flt, sex, yr);
          sel_slp_dev2_sum += ln_sel_slp_dev(1, flt, sex, yr);
          sel_inf_dev2_sum += sel_inf_dev(1, flt, sex, yr);
        }


        jnll_comp(6, flt) += square(sel_slp_dev1_sum) * 10000;
        jnll_comp(6, flt) += square(sel_inf_dev1_sum) * 10000;
        jnll_comp(6, flt) += square(sel_slp_dev2_sum) * 10000;
        jnll_comp(6, flt) += square(sel_inf_dev2_sum) * 10000;
      }
    }
  } // End selectivity loop



  // Slot 7 -- Add survey selectivity normalization (non-parametric)
  for(flt = 0; flt < n_flt; flt++){
    sp = flt_spp(flt);
    for(sex = 0; sex < nsex(sp); sex++){
      if(flt_type(flt) > 0){
        jnll_comp(6, flt) += 50 * square(avgsel(flt, sex));
      }
    }
  }


  // Slot 8-9 -- Survey catchability deviates
  for(flt = 0; flt < n_flt; flt++){

    // Prior on catchability
    if( est_srv_q(flt) == 2){
      jnll_comp(7, flt) -= dnorm(ln_srv_q(flt) - square(sigma_srv_q(flt))/2, ln_srv_q_prior(flt), sigma_srv_q(flt), true);
    }

    // Penalized/random deviate likelihood
    if(((srv_varying_q(flt) == 1) | (srv_varying_q(flt) == 2)) & (flt_type(flt) == 2)){
      for(yr = 0; yr < nyrs_hind; yr++){
        jnll_comp(8, flt) -= dnorm(ln_srv_q_dev(flt, yr), Type(0.0), time_varying_sigma_srv_q(flt), true );
      }
    }

    // Random walk
    if((srv_varying_q(flt) == 4) & (flt_type(flt) == 2)){

      // Sum to zero constraint
      Type q_dev_sum = 0;

      for(yr = 0; yr < nyrs_hind; yr++){
        q_dev_sum += ln_srv_q_dev(flt, yr);
      }
      jnll_comp(8, flt) += square(q_dev_sum) * 10000;

      for(yr = 1; yr < nyrs_hind; yr++){
        jnll_comp(8, flt) -= dnorm(ln_srv_q_dev(flt, yr) - ln_srv_q_dev(flt, yr-1), Type(0.0), time_varying_sigma_srv_q(flt), true );
      }
    }
  } // End q loop


  // Slots 10-12 -- RECRUITMENT DEVIATES
  for (sp = 0; sp < nspp; sp++) {
    // Slot 9 -- storck-recruit prior
    if(srr_use_prior == 1){
      jnll_comp(9, sp) -= dnorm(rec_pars(sp, 0), srr_prior_mean(sp), srr_prior_sd(sp), true);
    }

    // Slot 10 -- init_dev -- Initial abundance-at-age
    for (age = 1; age < nages(sp); age++) {
      jnll_comp(11, sp) -= dnorm( init_dev(sp, age - 1), Type(square(r_sigma(sp)) / 2), r_sigma(sp), true);
    }

    for (yr = 0; yr < nyrs_hind; yr++) {
      // Slot 11 -- Tau -- Annual recruitment deviation
      jnll_comp(10, sp) -= dnorm( rec_dev(sp, yr), Type(square(r_sigma(sp)) / 2), r_sigma(sp), true);    // Recruitment deviation using random effects.
    }
  }


  // Slot 13 -- Reference point penalties
  for (sp = 0; sp < nspp; sp++) {
    // -- Static reference points
    if((DynamicHCR == 0) & (forecast == 1)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(13, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-FsprLimit(sp));
      }

      // F that acheives \code{Ftarget}% of SSB0 in the end of the projection
      if(HCR == 3){
        jnll_comp(13, sp)  += 200*square((biomassSSB(sp, nyrs-1)/SB0(sp))-FsprTarget(sp));
      }

      // -- SPR
      if(HCR > 3){
        jnll_comp(13, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-FsprLimit(sp));
        jnll_comp(13, sp)  += 200*square((SPRtarget(sp)/SPR0(sp))-FsprTarget(sp));
      }
    }

    // -- Dynamic referene points
    if((DynamicHCR == 1) & (forecast == 1)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(13, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-FsprLimit(sp));
      }

      for(yr = 1; yr < nyrs; yr++){ // No initial abundance
        // F that acheives Ftarget% of SSB0y
        if(HCR == 3){
          jnll_comp(13, sp)  += 200*square((DynamicSBF(sp, yr)/DynamicSB0(sp, yr))-FsprTarget(sp));
        }
      }

      // -- SPR
      if(HCR > 3){
        jnll_comp(13, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-FsprLimit(sp));
        jnll_comp(13, sp)  += 200*square((SPRtarget(sp)/SPR0(sp))-FsprTarget(sp));
      }
    }
  }


  // Slot 14 -- N-at-age < 0 penalty. See posfun
  for (sp = 0; sp < nspp; sp++){
    if(estDynamics(sp) == 0){
      jnll_comp(14, sp) += zero_pop_pen(sp);
    }
  }


  // 14.3. Diet likelihood components from MSVPA with estimated suitability and Kinzey and Punt
  if ((msmMode > 2) | (suitMode > 0)) {
    // Slot 15 -- Diet weight likelihood
    for(int stom_ind = 0; stom_ind < UobsWtAge.rows(); stom_ind++){

      rsp = UobsWtAge_ctl(stom_ind, 0) - 1; // Index of pred
      if(UobsWtAge_hat(stom_ind, 1) > 0){
        jnll_comp(17, rsp) -= Type(UobsWtAge(stom_ind, 0)) * (UobsWtAge(stom_ind, 1)) * log(UobsWtAge_hat(stom_ind, 1)/UobsWtAge(stom_ind, 1));
      }
    }
  } // End diet proportion by weight component


  // 14.4. Diet likelihood components for Kinzey and Punt predation
  if (msmMode > 2){
    // Slot 15 -- Ration likelihood
    for (yr = 0; yr < nyrs_hind; yr++) {
      for (sp = 0; sp < nspp; sp++) {
        for(sex = 0; sex < nsex(sp); sex ++){
          for (age = 1; age < nages(sp); age++) { // don't include age zero in likelihood
            if(ration(sp, sex, age, yr) > 0){
              if(ration_hat(sp, sex, age, yr) > 0){
                jnll_comp(15, sp) -= dnorm(log(ration(sp, sex, age, yr))  - square(sd_ration) / 2,  log( ration_hat(sp, sex, age, yr)), sd_ration, true);
              }
            }
          }
        }
      }
    }


    // Slot 16 -- Ration penalalties FIXME: not sure if necessary: talk to Andre
    for (sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for (age = 0; age < nages(sp); age++) {
          for (yr = 0; yr < nyrs_hind; yr++) {
            if(ration_hat(sp, sex, age, yr) > 0){
              //jnll_comp(16, sp) += 20 *  pow(ration_hat(sp, sex, age, yr) - ration_hat_ave(sp, sex, age), 2);
            }
          }
        }
      }
    }
  } // End if statement for Kinzey diet likelihood


  // ------------------------------------------------------------------------- //
  // 12. REPORT SECTION                                                        //
  // ------------------------------------------------------------------------- //

  // 12.0 Report indices
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
  REPORT( est_srv_q );
  REPORT( srv_varying_q );
  REPORT( est_sigma_srv );
  REPORT( est_sigma_fsh );
  REPORT( UobsWtAge_ctl );
  REPORT( sel_norm_age );
  REPORT( r_sexes );
  REPORT( k_sexes );
  REPORT( joint_adjust );
  REPORT( sel_slp_dev_ll);
  REPORT( sel_inf_dev_ll);


  // 12.1. Population components
  REPORT( pop_scalar );
  REPORT( mean_rec );
  REPORT( pmature );
  REPORT( Zed );
  REPORT( NByage );
  REPORT( AvgN );
  REPORT( sex_ratio_hat );
  REPORT( sex_ratio_mean_hat );
  REPORT( S );
  REPORT( biomassByage );
  REPORT( biomassSSBByage );
  REPORT( depletion );
  REPORT( depletionSSB );
  REPORT( biomass );
  REPORT( biomassSSB );
  REPORT(r_sigma);
  REPORT( R_sexr );
  REPORT( R );
  REPORT( M );

  ADREPORT( B_eaten_as_prey );
  ADREPORT( M );
  ADREPORT( Zed );
  ADREPORT( depletion );
  ADREPORT( depletionSSB );
  ADREPORT( biomass );
  ADREPORT( biomassSSB );
  ADREPORT( r_sigma );
  ADREPORT( R );
  ADREPORT( pop_scalar );


  // -- 12.2. Biological reference points
  REPORT( NByage0);
  REPORT( DynamicNByage0);
  REPORT( DynamicNByageF);
  REPORT( NbyageSPR);
  REPORT( B0 );
  REPORT( SB0 );
  REPORT( DynamicB0 );
  REPORT( DynamicSB0 );
  REPORT( DynamicSBF );
  REPORT( SPR0 );
  REPORT( SPRlimit );
  REPORT( SPRtarget );
  REPORT( proj_F );
  REPORT( Ftarget );
  REPORT( Flimit );
  REPORT( FlimitSPR );
  REPORT( FtargetSPR );


  // -- 12.3. Selectivity
  REPORT( sel );
  REPORT( avgsel );
  REPORT( emp_sel_obs );
  REPORT( sel_tmp );
  REPORT( sigma_sel );
  REPORT( sel_curve_pen );


  // -- 12.4. Survey components
  REPORT( srv_bio_hat );
  REPORT( srv_log_sd_hat );
  REPORT( sigma_srv_index );
  REPORT( srv_q );
  REPORT( srv_q_analytical );
  REPORT( sigma_srv_q );
  REPORT( time_varying_sigma_srv_q );
  REPORT( ln_srv_q_dev );


  // -- 12.5. Fishery components
  REPORT( F_spp );
  REPORT( F_flt );
  REPORT( F_flt_age );
  REPORT( F_spp_age );
  REPORT( fsh_bio_hat );
  REPORT( fsh_log_sd_hat );


  // 12.6. Age/length composition
  REPORT( age_obs_hat );
  REPORT( age_hat );
  REPORT( comp_obs );
  REPORT( comp_hat );
  REPORT( true_age_comp_hat );
  REPORT( n_hat );
  REPORT( comp_n );


  // -- 12.7. Likelihood components
  REPORT( jnll_comp );


  // -- 12.8. Ration components
  REPORT( ConsumAge );
  REPORT( LbyAge );
  REPORT( mnWt_obs );
  REPORT( fT );
  REPORT( env_index_hat );
  REPORT( ration );


  // 12.9. Suitability components
  REPORT( suma_suit );
  REPORT( suit_main );
  REPORT( suit_other );
  REPORT( stom_div_bio2 );
  REPORT( diet_prop_weight );
  REPORT( avail_food );
  REPORT( other_food_diet_prop_weight );
  REPORT( M1 );
  REPORT( M2 );
  REPORT( M2_prop );
  REPORT( B_eaten_as_prey );
  REPORT( B_eaten_as_pred );
  REPORT( B_eaten );
  REPORT( N_eaten );
  REPORT( UobsWtAge_hat );
  REPORT( flt_accum_age_upper );
  REPORT( flt_accum_age_lower );


  // -- 12.10. Kinzey predation functions
  REPORT( H_1 );
  REPORT( H_1a );
  REPORT( H_1b );
  REPORT( H_2 );
  REPORT( H_3 );
  REPORT( gam_a );
  REPORT( gam_b );

  REPORT( N_pred_yrs );
  REPORT( N_prey_yrs );
  REPORT( N_pred_eq );
  REPORT( N_prey_eq );

  REPORT( pred_resp );
  REPORT( Pred_r );
  REPORT( Prey_r );

  REPORT( ration_hat );
  REPORT( ration_hat_ave );


  // ------------------------------------------------------------------------- //
  // 13. REPORT SIMULATED QUANTITIES SECTION                                   //
  // ------------------------------------------------------------------------- //
  SIMULATE {
    REPORT(srv_biom_obs);          // Report the simulation
    REPORT(fsh_biom_obs);
  }

  // ------------------------------------------------------------------------- //
  // 14. END MODEL                                                             //
  // ------------------------------------------------------------------------- //

  Type jnll = 0;

  // Estimation mode
  if (estimateMode < 3) {
    jnll = jnll_comp.sum();
  }

  // Debug mode
  if (estimateMode > 2) {
    jnll = dummy * dummy;
  }


  REPORT( jnll );
  return jnll;
}
