#include <TMB.hpp>
// ------------------------------------------------------------------------- //
//                 CEATTLE version 3.1.2                                     //
//                  Template Model Builder                                   //
//               Multispecies Statistical Model                              //
//          Bioenergetic-based Assessment for Understanding                  //
//              Biomass Linkages To The Environment                          //
//                  in the Bering Sea                                        //
//                                                                           //
// AUTHORS:   Kirstin Holsman, Jim Ianelli                                   //
//            Modified by Grant Adams                                        //
// CITATIONS:                                                                //
// 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
// ------------------------------------------------------------------------- //
//
// Changes from 2017 assessment
// 1. Normalize age-transition matrix prior to use (PCod rows did not sum to 1)
// 2. Added log-normal recruitment deviation bias correction
// 3. BT survey age/length composition estimated for month 6 rather than 0, similar to BT survey biomass
// 4. Fixed mis-specification of multinomial for fishery composition
// 5. Normalized survey and fishery selectivity so that max = 1
// 6. Fixed initialization of population
// 7. Fixed estimation routine of suitability coefficients (N used to caclulate suitability is no-longer different)
// 8. Changed ration calculation over nyrs rather than nTyrs
// 9. Allowed retrospective estimation
// 10. Acoustic comp and biomass data are assumed to come from month 6 rather than 0, given they survey in the summer
// 11. Added time varying selectivity and catchability
// 12. Can make survey selectivity and/or catchability equal across surveys
// 13. Added cabability to have a two sex model
// 14. Added in spawning stock biomass weight
// 15. Added in multiple weights for surveys
// 16. Added in spawning month mortality adjustment
// 17. Removed constant 0.0001 added to M1
// 18. Had the model estimate Uobs
// 19. Added flexibility to fixed n-at-age to have age, sex specific scaling paramters and estimate selectivity
//  Fixme: denominator is zero somewhere. Log of negative number. Check suitability. Make other prey a very large number.
//  Look at M2: suitability: and consumption. Make sure positive.
// 20. Added analytical q for time-varying survey sigma inputs
// 21. All kinzey predation bits are commented out in v&
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
  DATA_INTEGER(debug);                    // Logical to debug or not
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

  // DATA_INTEGER( random_rec );             // Logical of whether to treate recruitment deviations as random effects
  DATA_INTEGER( niter );                  // Number of loops for MSM mode


  // Columns =  Species, Survey, Selectivity, Estimate_q, log_q_prior
  // Selectivity 0 = empirical, 1 = logistic, 2 = non-parametric, 3 = double logistic

  // 1.2. Temporal dimensions
  DATA_INTEGER( styr );                      // Start year
  DATA_INTEGER( endyr );                     // End of estimation years
  DATA_INTEGER( projyr );                    // End year of projection
  int nyrs = projyr - styr + 1;
  int nyrs_hind = endyr - styr + 1;

  // 1.3. Number of species
  DATA_INTEGER( nspp );                   // Number of species (prey)
  DATA_IVECTOR( pop_wt_index );           // Dim 3 of wt to use for population dynamics
  DATA_IVECTOR( ssb_wt_index );           // Dim 3 of wt to use for spawning stock biomass calculation
  DATA_IVECTOR( pop_age_transition_index );          // Dim 3 of wt to use for age_transition_matrix
  DATA_IVECTOR( estDynamics );            // Index indicating wether the population parameters are estimated (0), numbers-at-age are provided (1), or an index of numbers-at-age multiplied by an estimated scalar is used (2)
  DATA_IVECTOR( est_sex_ratio );              // Is sex ration F/(M+F) to be included in the likelihood; 0 = no, 1 = use annual average across ages (uses 2nd age in sex_ratio data), 2 = age, and year specific (TBD)
  pop_wt_index -= 1;                      // Indexing starts at 0
  ssb_wt_index -= 1;                      // Indexing starts at 0
  pop_age_transition_index -= 1;                     // Indexing starts at 0


  // 1.4. MODEL OBJECTS
  // 1.4.1. LOOPING INDICES -- k = observation, sp = species, sex = sex (0 = combines; 1 = females; 2 = males), age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length
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
  DATA_VECTOR( fday );                    // number of foraging days for each predator; n = [1, nspp] #FIXME - assuming this is the same as fdays
  DATA_ARRAY( Pyrs );                     // n = [nspp, nyrs+1, nages]: #FIXME - Assuming this is the same as Pby_yr?
  DATA_MATRIX( UobsAge );                 // pred, prey, predA, preyA U observations (mean number of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_IMATRIX( UobsAge_ctl );            // Info on pred, prey, predA, preyA U matrix (mean number of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_MATRIX( UobsWtAge );               // pred, prey, predA, preyA U observations (mean wt_hat of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_IMATRIX( UobsWtAge_ctl );          // Info on pred, prey, predA, preyA U matrix (mean wt_hat of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_ARRAY( Mn_LatAge );                // Mean length-at-age; n = [nspp, sex, nages], ALSO: mean_laa in Kinzey

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

  // -- 2.4.2. von Bertalannfy growth function (VBGF): This is used to calculate future weight-at-age: NOT YET IMPLEMENTED

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
  PARAMETER_VECTOR( ln_mean_rec );                // Mean recruitment; n = [1, nspp]
  PARAMETER_VECTOR( ln_rec_sigma );               // Standard deviation of recruitment deviations; n = [1, nspp]
  PARAMETER_MATRIX( rec_dev );                    // Annual recruitment deviation; n = [nspp, nyrs]
  PARAMETER_VECTOR( ln_sex_ratio_sigma );         // Variance for sex ratio to be used; n = [nspp]


  // -- 3.2. Abundance parameters
  PARAMETER_MATRIX( init_dev );                   // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this
  PARAMETER_ARRAY( ln_M1 );                       // Natural mortality (residual if multispecies mode or total if single species mode); n = [nspp, nsex, nages]


  // -- 3.3. fishing mortality parameters
  PARAMETER_VECTOR( ln_mean_F );                  // Log mean fishing mortality; n = [1, n_fsh]
  PARAMETER_MATRIX( ln_FSPR );                    // Fishing mortality for projections; n = [1, nspp]
  PARAMETER_VECTOR( proj_F_prop );                // Proportion of fishing mortality from each fleet for projections; n = [1, n_fsh]
  PARAMETER_MATRIX( F_dev );                      // Annual fishing mortality deviations; n = [n_fsh, nyrs] # NOTE: The size of this will likely change

  // -- 3.4. Survey catchability parameters
  PARAMETER_VECTOR( ln_srv_q );                   // Survey catchability; n = [n_srv]
  PARAMETER_VECTOR( srv_q_pow );                  // Survey catchability power coefficient q * B ^ q_pow or beta ln(q_y) = q_mut + beta * index_y; n = [n_srv]
  PARAMETER_MATRIX( ln_srv_q_dev );               // Annual survey catchability deviates; n = [n_srv, nyrs_hind]
  PARAMETER_MATRIX( ln_srv_q_dev_re );            // Annual survey catchability random effect deviates; n = [n_srv, nyrs_hind]
  PARAMETER_VECTOR( ln_sigma_srv_q );             // Log standard deviation of prior on survey catchability; n = [1, n_srv]
  PARAMETER_VECTOR( ln_sigma_time_varying_srv_q );// Log standard deviation of time varying survey catchability; n = [1, n_srv]


  // -- 3.5. Selectivity parameters
  PARAMETER_ARRAY( sel_coff );                    // selectivity parameters for non-parametric; n = [n_selectivities, nsex, nselages]
  PARAMETER_ARRAY( ln_sel_slp );                  // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( sel_inf );                     // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( ln_sel_slp_dev );              // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( sel_inf_dev );                 // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( ln_sel_slp_dev_re );           // selectivity parameter random effect deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( sel_inf_dev_re );              // selectivity parameter random effect deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_VECTOR( ln_sigma_sel );               // Log standard deviation of selectivity; n = [1, n_selectivities]
  PARAMETER_MATRIX( sel_curve_pen );              // Selectivity penalty for non-parametric selectivity, 2nd column is for monotonic bit


  // -- 3.6. Variance of survey and fishery time series
  PARAMETER_VECTOR( ln_sigma_srv_index );         // Log standard deviation of survey index time-series; n = [1, n_srv]
  PARAMETER_VECTOR( ln_sigma_fsh_catch );         // Log standard deviation of fishery catch time-series; n = [1, n_fsh]
  PARAMETER_VECTOR( comp_weights );               // MacCallister-Ianelli weights for fisheries data

  // -- 3.7. Kinzery predation function parameters
  /* Commenting out kinzey bits for speed
  PARAMETER_MATRIX(logH_1);                       // Predation functional form; n = [nspp, nspp2];
  PARAMETER_VECTOR(logH_1a);                      // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix
  PARAMETER_VECTOR(logH_1b);                      // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix

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
  int max_bin = imax( nlengths );                                                   // Integer of maximum number of length/age bins.
  int n_flt = fleet_control.rows();
  vector<int> joint_adjust(comp_obs.rows()); joint_adjust.setZero();

  // -- 4.2. Estimated population quantities
  matrix<Type>  pop_scalar = ln_pop_scalar;  pop_scalar = exp(ln_pop_scalar.array());// Fixed n-at-age scaling coefficient; n = [nspp, nages]
  vector<Type>  mn_rec = exp(ln_mean_rec);                                          // Mean recruitment; n = [1, nspp]
  array<Type>   biomassByage(nspp, max_age, nyrs); biomassByage.setZero();          // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type>  biomass(nspp, nyrs); biomass.setZero();                             // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type>  biomassSSB(nspp, nyrs); biomassSSB.setZero();                       // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type>   biomassSSBByage(nspp, max_age, nyrs); biomassSSBByage.setZero();    // Spawning biomass at age (kg); n = [nspp, nages, nyrs]
  array<Type>   M(nspp, 2, max_age, nyrs); M.setZero();                             // Total natural mortality at age; n = [nyrs, 2 sexes, nages, nspp]
  array<Type>   M1(ln_M1.dim); M1 = exp(ln_M1);                                     // Residual or total natural mortality at age; n = [nspp, 2 sexes, nages]
  array<Type>   M2(nspp, 2, max_age, nyrs); M2.setZero();                           // Predation mortality at age; n = [nyrs, nages, nspp]
  array<Type>   M2_prop(nspp, nspp, 2, 2, max_age, max_age, nyrs); M2_prop.setZero();     // Relative predation mortality at age from each species at age; n = [nyrs, nages, nspp]
  array<Type>   NByage(nspp, 2, max_age, nyrs); NByage.setZero();                   // Numbers at age; n = [nspp, nages, nyrs]
  array<Type>   AvgN(nspp, 2, max_age, nyrs); AvgN.setZero();                       // Average numbers-at-age; n = [nspp, nages, nyrs]
  array<Type>   sex_ratio_hat(nspp, max_age, nyrs); sex_ratio_hat.setZero();        // Estimated age-specific sex ratin; n = [nspp, nages, nyrs]
  matrix<Type>  sex_ratio_mean_hat(nspp, nyrs); sex_ratio_mean_hat.setZero();       // Estimated sex ratio across all ages; n = [nspp, nyrs]
  matrix<Type>  R(nspp, nyrs); R.setZero();                                         // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type>   S(nspp, 2, max_age, nyrs); S.setZero();                             // Survival at age; n = [nspp, 2 sexes, nages, nyrs]
  array<Type>   Zed(nspp, 2, max_age, nyrs); Zed.setZero();                         // Total mortality at age; n = [nspp, 2 sexes, nages, nyrs]
  vector<Type>  r_sigma(nspp); r_sigma.setZero();                                   // Standard deviation of recruitment variation
  vector<Type>  zero_pop_pen(nspp); zero_pop_pen.setZero();                         // Additional penalty to add to likelihood if n-at-age goes < 0
  vector<Type>  sex_ratio_sigma(nspp); sex_ratio_sigma.setZero();                   // Variance of sex ratio; n = [nspp]

  // -- 4.3. Selectivity parameters
  matrix<Type>  avgsel(n_flt, 2); avgsel.setZero();                                 // Average selectivity; n = [1, nselectivities]
  array<Type>   sel(n_flt, 2, max_age, nyrs_hind); sel.setZero();                   // Estimated selectivity at age; n = [nselectivities, 2 sexes, nages, nyrs]
  array<Type>   sel_tmp(n_flt, 2, max_age); sel_tmp.setZero();                      // Temporary saved selectivity at age for estimated bits; n = [nselectivities, nages, nyrs]
  vector<Type>  sigma_sel(n_flt); sigma_sel.setZero();                              // Standard deviation of selectivity deviates; n = [1, nselectivities]

  // -- 4.4. Fishery components
  vector<Type>  sigma_fsh_catch(n_flt); sigma_fsh_catch.setZero();                  // Standard deviation of fishery time-series; n = [1, n_fsh]
  array<Type>   F(n_flt, 2, max_age, nyrs); F.setZero();                            // Estimated fishing mortality for each fishery; n = [n_flt, 2 sexes, nages, nyrs]
  array<Type>   F35_tot(nspp, 2, max_age); F35_tot.setZero();                       // Estimated fishing mortality for each species that leads to SB35; n = [nspp, 2 sexes, nages]
  array<Type>   F40_tot(nspp, 2, max_age); F40_tot.setZero();                       // Estimated fishing mortality for each species that leads to SB40; n = [nspp, 2 sexes, nages]
  array<Type>   F_tot(nspp, 2, max_age, nyrs+1); F_tot.setZero();                   // Sum of annual estimated fishing mortalities for each species-at-age; n = [nspp, 2 sexes, nages, nyrs]
  vector<Type>  fsh_bio_hat(fsh_biom_obs.rows()); fsh_bio_hat.setZero();            // Estimated fishery yield (kg);
  vector<Type>  fsh_log_sd_hat(fsh_biom_obs.rows()); fsh_log_sd_hat.setZero();      // Estimated/fixed fishery log_sd (kg);
  array<Type>   NbyageSPR(3, nspp, max_age);                                        // Estimated numbers at age for spawning biomass per recruit reference points
  vector<Type>  SB35(nspp); SB35.setZero();                                         // Estimated 35% spawning biomass per recruit
  vector<Type>  SB40(nspp); SB40.setZero();                                         // Estimated 40% spawning biomass per recruit
  vector<Type>  SB0(nspp); SB0.setZero();                                           // Estimated spawning biomass per recruit at F = 0
  matrix<Type>  FSPR(nspp, 2); FSPR = exp(ln_FSPR.array());
  matrix<Type>  proj_FABC(nspp, nyrs); proj_FABC.setZero();                         // Projected FABC using tier 3 harvest control rule
  // matrix<Type>  FSPR(nspp, 2); FSPR = exp(ln_FSPR.array());

  // -- 4.5. Survey components
  vector<Type>  sigma_srv_index(n_flt); sigma_srv_index.setZero();                  // Vector of standard deviation of survey index; n = [1, n_srv]
  vector<Type>  sigma_srv_q(n_flt); sigma_srv_q.setZero();                          // Vector of standard deviation of survey catchability prior; n = [1, n_srv]
  vector<Type>  time_varying_sigma_srv_q(n_flt); time_varying_sigma_srv_q.setZero();// Vector of standard deviation of time-varying survey catchability deviation; n = [1, n_srv]
  Type avgsel_tmp = 0;                                                              // Temporary object for average selectivity across all ages
  vector<Type>  srv_bio_hat(srv_biom_obs.rows()); srv_bio_hat.setZero();            // Estimated survey biomass (kg)
  vector<Type>  srv_log_sd_hat(srv_biom_obs.rows()); srv_log_sd_hat.setZero();      // Estimated/fixed survey log_sd (kg)
  vector<Type>  sigma_srv_analytical(n_flt); sigma_srv_analytical.setZero();        // Temporary vector to save analytical sigma follow Ludwig and Walters 1994
  vector<Type>  srv_q_analytical(n_flt); srv_q_analytical.setZero();                // Temporary vector to save analytical sigma follow Ludwig and Walters 1994
  matrix<Type>  srv_q(n_flt, nyrs_hind); srv_q.setZero();                           // Estimated survey catchability
  vector<Type>  srv_n_obs(n_flt); srv_n_obs.setZero();                              // Vector to save the number of observations for each survey time series

  // -- 4.6. Composition data - FIXME: will blow up in nlengths is less than nages
  vector<Type>  n_hat(comp_obs.rows()) ; n_hat.setZero() ;                          // Estimated catch (n); n = [nspp, nyrs]
  matrix<Type>  age_hat = comp_obs; age_hat.setZero();                              // Estimated catch at true age; n = [nspp, nages, nyrs]
  matrix<Type>  true_age_comp_hat = comp_obs; true_age_comp_hat.setZero();          // True age composition
  matrix<Type>  age_obs_hat = comp_obs; age_obs_hat.setZero();                      // Estimated catch at observed age; n = [nspp, nages, nyrs]
  matrix<Type>  comp_hat = comp_obs; comp_hat.setZero();                            // Estimated comp; n = [nspp, nages, nyrs]

  // -- 4.7. Ration components
  array<Type>   ConsumAge( nspp, 2, max_age, nyrs ); ConsumAge.setZero();           // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
  matrix<Type>  fT( nspp, nyrs ); fT.setZero();                                     // Pre-allocation of temperature function of consumption; n = [nspp, nTyrs]
  array<Type>   LbyAge( nspp, 2, max_age, nyrs ); LbyAge.setZero();                 // Length by age from LW regression
  matrix<Type>  mnWt_obs( nspp, max_age ); mnWt_obs.setZero();                      // Mean observed weight at age (across years); n = [nspp, nages]
  array<Type>   ration2Age( nspp, 2, max_age, nyrs ); ration2Age.setZero();         // Annual ration at age (kg/yr); n = [nyrs, nages, nspp]
  matrix<Type>  env_index_hat( nyrs, env_index.cols() ); env_index_hat.setZero();   // Environmental indices like bottom temperature; n = [1, nTyrs]

  // -- 4.8. Suitability components
  array<Type>   avail_food(nspp, 2, max_age, nyrs); avail_food.setZero();                        // Available food to predator; n = [nspp, 2 sexes, nages, nyrs]
  array<Type>   othersuit(nspp, 2, max_age, nyrs); othersuit.setZero();                          // Sum of available food to predator; n = [nspp, 2 sexes, nages, nyrs]
  array<Type>   B_eaten(nspp, 2, max_age, nyrs); B_eaten.setZero();                              // Biomass of prey eaten via predation; n = [nyrs, nages, nspp]
  array<Type>   B_eaten_prop(nspp, nspp, 2, 2, max_age, max_age, nyrs); B_eaten.setZero();       // Biomass of prey eaten via predation by a predator at age; n = [nyrs, nages, nspp]
  array<Type>   of_stomKir(nspp, 2, max_age, nyrs); of_stomKir.setZero();                        // Other food stomach content; n = [nspp, 2 sexes, nages, nyrs]
  array<Type>   stom_div_bio2(nspp, nspp, 2, 2, max_age, max_age, nyrs); stom_div_bio2.setZero();// Stomach proportion over biomass; U/ (W * N) ; n = [nspp, nspp, 2 sexes, 2 sexes, nages, nages, nyrs]
  array<Type>   stomKir(nspp, nspp, 2, 2, max_age, max_age, nyrs); stomKir.setZero();            // Stomach proportion by numbers U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   suit_tmp(nspp, nspp, 2, 2, max_age, max_age, nyrs); suit_tmp.setZero();       // Temporary suitability storage U; n = [nspp, nspp, 2 sexes, 2 sexes, nages, nages, nyrs]
  array<Type>   stomKirWt(nspp, nspp, 2, 2, max_age, max_age, nyrs); stomKirWt.setZero();     // Stomach proportion by weight U; n = [nspp, nspp, nages, nages, nyrs]
  // array<Type>   diet_w_dat(nspp, nspp, max_bin, max_bin, nyrs); diet_w_dat.setZero();         // Observed stomach contents by weight of prey length j in predator length l
  // array<Type>   diet_w_sum(nspp, nspp, max_bin, nyrs); diet_w_sum.setZero();                  // Observed stomach contentes by weight of prey in predator length j
  array<Type>   suit_main(nspp, nspp,  2, 2, max_age, max_age, nyrs); suit_main.setZero();    // Suitability/gamma selectivity of predator age u on prey age a; n = [nspp, nspp, 2 sexes, 2 sexes, nages, nages]
  array<Type>   suit_other(nspp, 2, max_age); suit_other.setZero();                           // Suitability not accounted for by the included prey; n = [nspp, 2 sexes, nages]
  array<Type>   suma_suit(nspp, 2, max_age, nyrs); suma_suit.setZero();                       // Sum of suitabilities; n = [nyrs, 2 sexes, nages, nspp]
  matrix<Type>  UobsAge_hat = UobsAge; UobsAge.setZero();                                     // Estimated stomach proportion by weight U; n = [nspp, nspp, nages, nages, nyrs]
  matrix<Type>  UobsWtAge_hat = UobsWtAge; UobsWtAge_hat.setZero();                           // Estimated stomach proportion by weight U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   Uobs(nspp, nspp, max_bin, max_bin);Uobs.setZero();                            // pred, prey, predL, preyL U matrix (mean number of prey in each pred); n = [nspp, nspp, maxL, maxL]
  array<Type>   UobsWt(nspp, nspp, max_bin, max_bin);UobsWt.setZero();                         // pred, prey, predL, preyL U matrix (mean wt_hat of prey in each pred); n = [nspp, nspp, maxL, maxL] #FIXME - Changed name in stomach2017.dat

  // -- 4.9. Kinzey selectivity
  vector<Type> gam_a = exp(log_gam_a); // Predator selectivity
  vector<Type> gam_b = exp(log_gam_b); // Predator selectivity

  // -- 4.10. Kinzey Functional response
  /* Commenting out kinzey bits for speed
  matrix<Type> H_1(nspp, nspp + 1); H_1 = exp(logH_1.array());
  vector<Type> H_1a(nspp); H_1a = exp(logH_1a);
  vector<Type> H_1b(nspp); H_1b = exp(logH_1b);
  matrix<Type> H_2(nspp, nspp); H_2 = exp(logH_2.array());
  matrix<Type> H_3(nspp, nspp); H_3 = exp(logH_3.array());

  array<Type>  N_pred_yrs(nspp, 2, max_age, nyrs); N_pred_yrs.setZero();                // Effective numbers of predators for each age of prey
  array<Type>  N_prey_yrs(nspp, 2, max_age, nyrs); N_prey_yrs.setZero();                // Effective numbers of prey for each age of predator
  array<Type>  N_pred_eq(nspp, 2, max_age); N_pred_eq.setZero();                        // Effective numbers of predators for each age of prey (styr_pred)
  array<Type>  N_prey_eq(nspp, 2, max_age); N_prey_eq.setZero();                        // Effective numbers of prey for each age of predator

  array<Type>  pred_resp(nspp, nspp + 1, 2, 2, max_age, max_age, nyrs); pred_resp.setZero();// Predator functional response
  array<Type>  Pred_r(nspp, 2, max_age, nyrs); Pred_r.setZero();                        // save Pred_ratio values
  array<Type>  Prey_r(nspp, 2, max_age, nyrs); Prey_r.setZero();                        // save Prey_ratio values

  array<Type>  Vmort_ua(nspp, nspp, 2, 2, max_age, max_age, nyrs); Vmort_ua.setZero();     // Predation mortality on prey age a by single predator age u
  array<Type>  Pmort_ua(nspp, nspp, 2, 2, max_age, nyrs); Pmort_ua.setZero();              // Predation mortality on prey age a by all predators age u

  // FIXME - add flexible diet components
  array<Type>  eaten_la(nspp, nspp, max_bin, max_age, nyrs); eaten_la.setZero();     // Number of prey of age a eaten by predator length l
  array<Type>  eaten_ua(nspp, nspp, max_age, max_age, nyrs); eaten_ua.setZero();     // Number of prey of age a eaten by predator age u

  array<Type>  Q_mass_l(nspp, nspp + 1, max_bin, nyrs); Q_mass_l.setZero();          // Mass of each prey sp consumed by predator at length // FIXME: make into 4D array
  array<Type>  Q_mass_u(nspp, nspp + 1, max_bin, nyrs); Q_mass_u.setZero();          // Mass of each prey sp consumed by predator at age // FIXME: make into 4D array
  matrix<Type> Q_other_u(nspp, max_age); Q_other_u.setZero();                        // Mass of other prey consumed by predator at age
  array<Type>  Q_hat(nspp, nspp + 1, max_bin, nyrs); Q_hat.setZero();                // Fraction for each prey type of total mass eaten by predator length
  array<Type>  T_hat(nspp, nspp, max_bin, max_bin); T_hat.setZero();                 // Fraction of prey of length m in predator of length l

  array<Type> omega_hat(nspp, max_age, nyrs); omega_hat.setZero();                   // Daily ration by predator age each year
  matrix<Type> omega_hat_ave(nspp, max_age); omega_hat_ave.setZero();                // Daily ration by predator age averaged over years
  */


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
        srv_q(flt, yr) = exp(ln_srv_q(flt) + ln_srv_q_dev_re(flt, yr));                 // Exponentiate
      }

      // Q as a function of environmental index
      if(est_srv_q(flt) == 5){
        srv_q(flt, yr) = exp(ln_srv_q(flt) + srv_q_pow(flt) * env_index_hat(yr, srv_varying_q(flt)));
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



      // 6.1.1. Logisitic selectivity
      if (sel_type == 1) {
        for (age = 0; age < nages(sp); age++){
          for (yr = 0; yr < nyrs_hind; yr++) {
            for(sex = 0; sex < nsex(sp); sex++){
              // Random walk and block
              if(sel_varying != 2){
                sel(flt, sex, age, yr) = 1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr)) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr )))));
              }

              // Random effect
              if(sel_varying == 2){
                sel(flt, sex, age, yr) = 1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev_re(0, flt, sex, yr)) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev_re(0, flt, sex, yr )))));
              }
            }
          }
        }
      }

      // 6.1.2. Non-parametric selectivity fit to age ranges. NOTE: This can likely be improved
      if (sel_type == 2) {
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
              for(sex = 0; sex < nsex(sp); sex++){
                sel(flt, sex, age, yr) = sel_tmp(flt, sex, age);
              }
            }
          }
        }
      }


      // 6.1.3. Double logistic (Dorn and Methot 1990)
      if (sel_type == 3) {
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
                sel(flt, sex, age, yr) = 1 / (1 + exp( -exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev_re(0, flt, sex, yr)) * ((age + 1) -( sel_inf(0, flt, sex) + sel_inf_dev_re(0, flt, sex, yr ))))) * // Upper slope
                  (1 - (1 / (1 + exp( -exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev_re(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev_re(1, flt, sex, yr )))))));  // Downward slope;
              }
            }
          }
        }
      }


      // 6.1.4. Descending logistic (Dorn and Methot 1990)
      if (sel_type == 4) {
        for (age = 0; age < nages(sp); age++){
          for (yr = 0; yr < nyrs_hind; yr++) {
            for(sex = 0; sex < nsex(sp); sex++){

              // Random walk and block
              if(sel_varying != 2){
                sel(flt, sex, age, yr) = (1 - (1 / (1 + exp( - exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr )))))));  // Downward slope;
              }

              // Random effect
              if(sel_varying == 2){
                sel(flt, sex, age, yr) = (1 - (1 / (1 + exp( - exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev_re(1, flt, sex, yr)) * ((age + 1) - (sel_inf(1, flt, sex) + sel_inf_dev_re(1, flt, sex, yr )))))));  // Downward slope;
              }
            }
          }
        }
      }


      // 6.1.4. Normalize
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
        // FIXME: there is probably a more elegant way to do this
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



    // 6.1. ESTIMATE FISHING MORTALITY
    F.setZero();
    F_tot.setZero();
    F40_tot.setZero();
    F35_tot.setZero();
    for (flt = 0; flt < n_flt; flt++) {

      sp = flt_spp(flt);  // Temporary index of fishery survey

if(flt_type(flt) == 1){
      for (age = 0; age < nages(sp); age++) {
        for(sex = 0; sex < nsex(sp); sex ++){
          for (yr = 0; yr < nyrs; yr++) {
            // Hindcast
            if( yr < nyrs_hind){
              F(flt, sex, age, yr) = sel(flt, sex, age, yr) * exp(ln_mean_F(flt) + F_dev(flt, yr));
            }

            // Forecast
            if( yr >= nyrs_hind){
              F(flt, sex, age, yr) = sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * FSPR(sp, 1); // FIXME using last year of selectivity
            }

              F_tot(sp, sex, age, yr) += F(flt, sex, age, yr);
          }

            F35_tot(sp, sex, age) += sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * FSPR(sp, 0); // F35_tot%
            F40_tot(sp, sex, age) += sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * FSPR(sp, 1); // F40_tot%
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
            Zed(sp, sex, age, yr) = M1(sp, sex, age) + F_tot(sp, sex, age, yr) + M2(sp, sex, age, yr);
            S(sp, sex, age, yr) = exp(-Zed(sp, sex, age, yr));
          }
        }
      }
    }


    // 6.3. SPR based reference points
    SB0.setZero();
    SB35.setZero();
    SB40.setZero();
    for (sp = 0; sp < nspp; sp++) {
      NbyageSPR(0, sp, 0) = exp(ln_mean_rec(sp)); // F = 0
      NbyageSPR(1, sp, 0) = exp(ln_mean_rec(sp)); // F35
      NbyageSPR(2, sp, 0) = exp(ln_mean_rec(sp)); // F40

      for (age = 1; age < nages(sp)-1; age++) {
        NbyageSPR(0, sp, age) =  NbyageSPR(0, sp, age-1) * exp(-M(sp, 0, age-1, nyrs_hind - 1));
        NbyageSPR(1, sp, age) =  NbyageSPR(1, sp, age-1) * exp(-(M(sp, 0, age-1, nyrs_hind - 1) + F35_tot(sp, 0, age-1)));
        NbyageSPR(2, sp, age) =  NbyageSPR(2, sp, age-1) * exp(-(M(sp, 0, age-1, nyrs_hind - 1) + F40_tot(sp, 0, age-1)));
      }

      // Plus group
      NbyageSPR(0, sp, nages(sp) - 1) = NbyageSPR(0, sp, nages(sp) - 2) * exp(-M(sp, 0, nages(sp) - 2, nyrs_hind - 1)) / (1 - exp(-M(sp, 0, nages(sp) - 1, nyrs_hind - 1)));
      NbyageSPR(1, sp, nages(sp) - 1) = NbyageSPR(1, sp, nages(sp) - 2) * exp(-(M(sp, 0,  nages(sp) - 2, nyrs_hind - 1) + F35_tot(sp, 0,  nages(sp) - 2))) / (1 - exp(-(M(sp, 0,  nages(sp) - 1, nyrs_hind - 1) + F35_tot(sp, 0,  nages(sp) - 1))));
      NbyageSPR(2, sp, nages(sp) - 1) = NbyageSPR(2, sp, nages(sp) - 2) * exp(-(M(sp, 0,  nages(sp) - 2, nyrs_hind - 1) + F40_tot(sp, 0,  nages(sp) - 2))) / (1 - exp(-(M(sp, 0,  nages(sp) - 1, nyrs_hind - 1) + F40_tot(sp, 0,  nages(sp) - 1))));

      // Calulcate SPR
      for (age = 0; age < nages(sp); age++) {
        SB0(sp) +=  NbyageSPR(0, sp, age) *  wt( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * pmature(sp, age) * exp(-M(sp, 0, nages(sp) - 1, nyrs_hind - 1) * spawn_month(sp)/12);
        SB35(sp) +=  NbyageSPR(1, sp, age) *  wt( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * pmature(sp, age) * exp(-(M(sp, 0,  age, nyrs_hind - 1) + F35_tot(sp, 0,  age)) * spawn_month(sp)/12);
        SB40(sp) +=  NbyageSPR(2, sp, age) *  wt( ssb_wt_index(sp), 0, age, (nyrs_hind - 1) ) * pmature(sp, age) * exp(-(M(sp, 0,  age, nyrs_hind - 1) + F40_tot(sp, 0,  age)) * spawn_month(sp)/12);
      }
    }


    // 6.4. ESTIMATE RECRUITMENT T1.1
    for (sp = 0; sp < nspp; sp++) {

      if(nsex(sp)  == 1){
        R_sexr(sp) = 1;
      }

      for(sex = 0; sex < nsex(sp); sex ++){
        for (yr = 0; yr < nyrs; yr++) {
          R(sp, yr) = exp(ln_mean_rec(sp) + rec_dev(sp, yr));

          // Estimated numbers-at-age
          if(estDynamics(sp) == 0){
            NByage(sp, sex, 0, yr) = R(sp, yr) * R_sexr(sp);
          }
          // Fixed numbers-at-age
          if(estDynamics(sp) > 0){
            NByage(sp, sex, 0, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, 0, yr);
          }
        }
      }
    }


    // 6.5. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
    for (sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for (age = 0; age < nages(sp); age++) {
          // Estimated abundance
          if(estDynamics(sp) == 0){

            if ((age > 0) & (age < nages(sp) - 1)) {

              // Sum M1 until age - 1
              Type mort_sum = 0;
              for(int age_tmp = 0; age_tmp < age; age_tmp++){
                mort_sum += M1(sp, sex, age_tmp);
              }
              NByage(sp, sex, age, 0) = exp(ln_mean_rec(sp) - mort_sum + init_dev(sp, age - 1)) * R_sexr(sp);
            }
            // -- 6.2.2. Where yr = 1 and age > Ai.
            if (age == (nages(sp) - 1)) {
              // Sum M1 until age - 1
              Type mort_sum = 0;
              for(int age_tmp = 0; age_tmp < age; age_tmp++){
                mort_sum += M1(sp, sex, age_tmp);
              }
              NByage(sp, sex, age, 0) = exp(ln_mean_rec(sp) - mort_sum + init_dev(sp, age - 1)) / (1 - exp(-M1(sp, sex, nages(sp) - 1))) * R_sexr(sp); // NOTE: This solves for the geometric series
            }
          }

          // Fixed numbers-at-age - fixed scalar
          if(estDynamics(sp) == 1){
            NByage(sp, sex, age, 0) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, 0);
          }

          // Fixed numbers-at-age age-independent scalar
          if(estDynamics(sp) == 2){
            NByage(sp, sex, age, 0) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, 0);
          }

          // Fixed numbers-at-age age-dependent scalar
          if(estDynamics(sp) == 3){
            NByage(sp, sex, age, 0) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, 0);
          }
        }
      }
    }


    // 6.6. ESTIMATE HINDCAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
    biomass.setZero();
    biomassSSB.setZero();
    biomassByage.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 1; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            // Estimated numbers-at-age
            if(estDynamics(sp) == 0){
              // -- 6.3.1.  Where 1 <= age < Ai
              if (age < (nages(sp) - 1)) {
                NByage(sp, sex, age + 1, yr) = NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }

              // -- 6.3.2. Plus group where age > Ai. NOTE: This is not the same as T1.3 because sp used age = A_i rather than age > A_i.
              if (age == (nages(sp) - 1)) {
                NByage(sp, sex, age, yr) = NByage(sp, sex, age - 1, yr - 1) * S(sp, sex, age - 1, yr - 1) + NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }
            }

            // Fixed numbers-at-age - fixed scalar
            if(estDynamics(sp) == 1){
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
            }

            // Fixed numbers-at-age age-independent scalar
            if(estDynamics(sp) == 2){
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
            }

            // Fixed numbers-at-age age-dependent scalar
            if(estDynamics(sp) == 3){
              NByage(sp, sex, age, yr) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, yr);
            }

            // Constrain to reduce population collapse
            // NByage(sp, sex, age, yr) = posfun(NByage(sp, sex, age, yr), Type(0.001), zero_pop_pen(sp));

          }
        }

        // -- 6.3.3. Estimate Biomass and SSB
        for (yr = 0; yr < nyrs_hind; yr++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            biomassByage(sp, age, yr) += NByage(sp, sex, age, yr) * wt( pop_wt_index(sp), sex, age, yr ); // 6.5.
            biomassSSBByage(sp, age, yr) = NByage(sp, 0, age, yr) * pow(S(sp, 0, age, yr), spawn_month(sp)/12) * wt(ssb_wt_index(sp), 0, age, yr ) * pmature(sp, age); // 6.6.

          }
          biomass(sp, yr) += biomassByage(sp, age, yr);
          biomassSSB(sp, yr) += biomassSSBByage(sp, age, yr);
        }
      }
    }



    // 6.7. ESTIMATE FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
    // Tier 3 Harvest control
    for (sp = 0; sp < nspp; sp++) {
      for (yr = nyrs_hind; yr < nyrs; yr++){
        for (age = 0; age < nages(sp); age++) {

          for(sex = 0; sex < nsex(sp); sex ++){

            // Estimated numbers-at-age
            if(estDynamics(sp) == 0){
              // -- 6.3.1.  Where 1 <= age < Ai
              if (age < (nages(sp) - 1)) {
                NByage(sp, sex, age + 1, yr) = NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }

              // -- 6.3.2. Plus group where age > Ai. NOTE: This is not the same as T1.3 because sp used age = A_i rather than age > A_i.
              if (age == (nages(sp) - 1)) {
                NByage(sp, sex, age, yr) = NByage(sp, sex, age - 1, yr - 1) * S(sp, sex, age - 1, yr - 1) + NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }
            }


            // Fixed numbers-at-age - fixed scalar
            if(estDynamics(sp) == 1){
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
            }

            // Fixed numbers-at-age age-independent scalar
            if(estDynamics(sp) == 2){
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
            }

            // Fixed numbers-at-age age-dependent scalar
            if(estDynamics(sp) == 3){
              NByage(sp, sex, age, yr) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, yr);
            }

            // Hard constraint to reduce population collapse
            if(NByage(sp, sex, age, yr) < minNByage){
              NByage(sp, sex, age, yr) = minNByage;
            }


            // -- 6.3.3. Estimate Biomass and SSB
            biomassByage(sp, age, yr) += NByage(sp, sex, age, yr) * wt( pop_wt_index(sp), sex, age, nyrs_hind-1 ); // 6.5.
            biomassSSBByage(sp, age, yr) = NByage(sp, 0, age, yr) * pow(S(sp, 0, age, yr), spawn_month(sp)/12) * wt(ssb_wt_index(sp), 0, age, nyrs_hind-1 ) * pmature(sp, age); // 6.6.

            biomass(sp, yr) += biomassByage(sp, age, yr);
            biomassSSB(sp, yr) += biomassSSBByage(sp, age, yr);
          }
        }


        // Adjust F based on Tier 3 control rule
        if(biomassSSB(sp, yr) < SB40(sp)){
          proj_FABC(sp, yr) = FSPR(sp, 1) * (((biomassSSB(sp, yr)/SB40(sp))-0.05)/(1-0.05)); // Used Fabc of F40_tot%


          for (age = 0; age < nages(sp); age++) {
            for(sex = 0; sex < nsex(sp); sex ++){
              F_tot(sp, sex, age, yr) = 0;
            }
          }


          for (flt = 0; flt < n_flt; flt++) {
            if(sp == flt_spp(flt)){
              for (age = 0; age < nages(sp); age++) {
                for(sex = 0; sex < nsex(sp); sex ++){

                  F(flt, sex, age, yr) = sel(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * proj_FABC(sp, yr); // FIXME using last year of selectivity

                  if(flt_type(flt) == 1){
                    F_tot(sp, sex, age, yr) += F(flt, sex, age, yr);
                  }
                }
              }
            }
          }
        }

        for (age = 0; age < nages(sp); age++) {
          biomassByage(sp, age, yr) = 0;
          for(sex = 0; sex < nsex(sp); sex ++){
            M(sp, sex, age, yr) = M1(sp, sex, age) + M2(sp, sex, age, yr);
            Zed(sp, sex, age, yr) = M1(sp, sex, age) + F_tot(sp, sex, age, yr) + M2(sp, sex, age, yr);
            S(sp, sex, age, yr) = exp(-Zed(sp, sex, age, yr));
          }
        }

        biomass(sp, yr) = 0;
        biomassSSB(sp, yr) = 0;


        // Rerun projection
        for (age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){

            // Estimated numbers-at-age
            if(estDynamics(sp) == 0){
              // -- 6.3.1.  Where 1 <= age < Ai
              if (age < (nages(sp) - 1)) {
                NByage(sp, sex, age + 1, yr) = NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }

              // -- 6.3.2. Plus group where age > Ai. NOTE: This is not the same as T1.3 because sp used age = A_i rather than age > A_i.
              if (age == (nages(sp) - 1)) {
                NByage(sp, sex, age, yr) = NByage(sp, sex, age - 1, yr - 1) * S(sp, sex, age - 1, yr - 1) + NByage(sp, sex, age, yr - 1) * S(sp, sex, age, yr - 1);
              }
            }


            // Fixed numbers-at-age - fixed scalar
            if(estDynamics(sp) == 1){
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
            }

            // Fixed numbers-at-age age-independent scalar
            if(estDynamics(sp) == 2){
              NByage(sp, sex, age, yr) = pop_scalar(sp, 0) * NByageFixed(sp, sex, age, yr);
            }

            // Fixed numbers-at-age age-dependent scalar
            if(estDynamics(sp) == 3){
              NByage(sp, sex, age, yr) = pop_scalar(sp, age) * NByageFixed(sp, sex, age, yr);
            }

            // Hard constraint to reduce population collapse
            if(NByage(sp, sex, age, yr) < minNByage){
              NByage(sp, sex, age, yr) = minNByage;
            }


            // -- 6.3.3. Estimate Biomass and SSB
            biomassByage(sp, age, yr) += NByage(sp, sex, age, yr) * wt( pop_wt_index(sp), sex, age, nyrs_hind-1 ); // 6.5.
            biomassSSBByage(sp, age, yr) = NByage(sp, 0, age, yr) * pow(S(sp, 0, age, yr), spawn_month(sp)/12) * wt(ssb_wt_index(sp), 0, age, nyrs_hind-1 ) * pmature(sp, age); // 6.6.

            biomass(sp, yr) += biomassByage(sp, age, yr);
            biomassSSB(sp, yr) += biomassSSBByage(sp, age, yr);
          }
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
            if (avgnMode == 0) {
              AvgN(sp, sex, age, yr) = NByage(sp, sex, age, yr) * (1 - S(sp, sex, age, yr)) / Zed(sp, sex, age, yr); // MSVPA approach
            }else if (avgnMode == 1) {
              AvgN(sp, sex, age, yr) = NByage(sp, sex, age, yr) * exp(- Zed(sp, sex, age, yr) / 2); // Kinzey and Punt (2009) approximation
            }else if (avgnMode == 2) {
              AvgN(sp, sex, age, yr) = NByage(sp, sex, age, yr); // Van Kirk et al (2010) approximation
            } else{
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

        // Exponential function from Stewart et al. 1983
        if ( Ceq(sp) == 1) {
          fT(sp, yr) = exp(Qc(sp) * env_index_hat(yr, Cindex(sp)));
        }

        // Temperature dependence for warm-water-species from Kitchell et al. 1977
        if ( Ceq(sp) == 2) {
          Yc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp) + 2);
          Zc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp));
          Vc = (Tcm(sp) - env_index_hat(yr, Cindex(sp))) / (Tcm(sp) - Tco(sp));
          Xc = pow(Zc, 2) * pow((1 + pow((1 + 40 / Yc), 0.5)), 2) / 400;
          fT(sp, yr) = pow(Vc, Xc) * exp(Xc * (1 - Vc));
        }

        // Temperature dependence for cool and cold-water species from Thornton and Lessem 1979
        if (Ceq(sp) == 3) {
          G2 = (1 / (Tcl(sp) - Tcm(sp))) * log((0.98 * (1 - CK4(sp))) / (CK4(sp) * 0.02));
          L2 = exp(G2 * (Tcl( sp ) -  env_index_hat(yr, Cindex(sp))));
          Kb = (CK4(sp) * L2) / (1 + CK4(sp) * (L2 - 1));
          G1 = (1 / (Tco(sp) - Qc(sp))) * log((0.98 * (1 - CK1(sp))) / (CK1(sp) * 0.02));
          L1 = exp(G1 * (env_index_hat(yr, Cindex(sp)) - Qc(sp)));
          Ka = (CK1(sp) * L1) / (1 + CK1(sp) * (L1 - 1));
          fT(sp, yr) = Ka * Kb;
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
              ConsumAge(sp, sex, age, yr) = CA(sp) * pow(wt( pop_wt_index(sp), sex, age, yr ) * Type(1000), CB( sp )) //  C_max = CA * W ^ CB; where C_max is grams consumed per grams of predator per day
              * fT(sp, yr) * fday( sp ) * wt(pop_wt_index(sp), sex, age, yr) * 1000;                           //  C_max * f(T) * wt * fday g/pred.yr
              ConsumAge(sp, sex, age, yr) = ConsumAge(sp, sex, age, yr) * Pvalue(sp) * Pyrs(sp, sex, age, yr); //
            }

            // Projection
            if(yr >= nyrs_hind){
              ConsumAge(sp, sex, age, yr) = CA(sp) * pow(wt( pop_wt_index(sp), sex, age, (nyrs_hind - 1) ) * Type(1000), CB( sp ))  //  C_max = CA * W ^ CB; where C_max is grams consumed per grams of predator per day
              * fT(sp, yr) * fday( sp ) * wt( pop_wt_index(sp), sex, age, (nyrs_hind - 1) ) * 1000;                            //  C_max * f(T) * wt * fday g/pred.yr
              ConsumAge(sp, sex, age, yr) = ConsumAge(sp, sex, age, yr) * Pvalue(sp) * Pyrs(sp, sex, age, (nyrs_hind-1)); //
            }

            ration2Age(sp, sex, age, yr) = ConsumAge(sp, sex, age, yr) / 1000;      // Annual ration kg/yr //aLW(predd)*pow(lengths(predd,age),bLW(predd));//mnwt_bin(predd,age);
          }
        }
      }
    }


    // 7.4. Reorganize UobsAge content
    /*
     for(int stom_ind = 0; stom_ind < UobsAge.rows(); stom_ind++){
     rsp = UobsAge_ctl(stom_ind, 0) - 1; // Index of pred
     ksp = UobsAge_ctl(stom_ind, 1) - 1; // Index of prey
     r_sex = UobsAge_ctl(stom_ind, 2); // Index of pred sex
     k_sex = UobsAge_ctl(stom_ind, 3); // Index of prey sex
     r_age = UobsAge_ctl(stom_ind, 4) - minage(rsp); // Index of pred age
     k_age = UobsAge_ctl(stom_ind, 5) - minage(ksp); // Index of prey age
     yr = UobsAge_ctl(stom_ind, 6) - styr; // Index of year

     // Predator
     // 1 sex model
     vector<int> r_sexes(1); r_sexes(0) = 0;

     // 2 sex model and Uobs is for both sex
     if( (nsex(rsp) == 2) & (r_sex == 0) ){
     vector<int> r_sexes(2); r_sexes(0) = 0; r_sexes(1) = 1;
     }
     // 2 sex model and Uobs is for 1 sex
     if( (nsex(rsp) == 2) & (r_sex > 0) ){
     vector<int> r_sexes(1); r_sexes(0) = r_sex;
     }

     // Prey
     // 1 sex model
     vector<int> k_sexes(1); k_sexes(0) = 0;

     // 2 sex model and Uobs is for both sex
     if( (nsex(ksp) == 2) & (k_sex == 0) ){
     vector<int> k_sexes(2); k_sexes(0) = 0; k_sexes(1) = 1;
     }
     // 2 sex model and Uobs is for 1 sex
     if( (nsex(ksp) == 2) & (k_sex > 0) ){
     vector<int> k_sexes(1); k_sexes(0) = k_sex;
     }

     // Only use hindcast
     if(yr < nyrs_hind){

     for(r_sex = 0; r_sex < r_sexes.size(); r_sex ++){
     for(k_sex = 0; k_sex < k_sexes.size(); k_sex ++){
     // Average of years
     if(yr == -styr){
     for (yr = 0; yr < nyrs; yr++) {
     stomKir(rsp, ksp, r_sexes(r_sex), k_sexes(k_sex), r_age, k_age, yr) = UobsAge(stom_ind, 1);
     }
     } else { // Not average
     stomKir(rsp, ksp, r_sexes(r_sex), k_sexes(k_sex), r_age, k_age, yr) = UobsAge(stom_ind, 1);
     }
     }
     }
     }
     }
     */



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
              stomKirWt(rsp, ksp, r_sexes(stom_ind, j), k_sexes(stom_ind, k), r_age, k_age, yr) = UobsWtAge(stom_ind, 1);
            }
          }

          // Average of years
          if(flt_yr == 0){
            for (yr = 0; yr < nyrs; yr++) {
              stomKirWt(rsp, ksp, r_sexes(stom_ind, j), k_sexes(stom_ind, k), r_age, k_age, yr) = UobsWtAge(stom_ind, 1);
            }
          }
        }
      }
    }



    /* Turning off kinzey bits
    diet_w_sum.setZero();
    for (yr = 0; yr < nyrs; yr++) {                             // Year loop
      for (rsp = 0; rsp < nspp; rsp++) {                        // Predator species loop
        for (ksp = 0; ksp < nspp; ksp++) {                      // Prey species loop
          // By length
          for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) {      // Predator length loop // FIXME: Change to diet length bins
            for (k_ln = 0; k_ln < nlengths(ksp); k_ln++) {    // Prey length loop FIXME: Change to diet length bins
              diet_w_dat(rsp, ksp, r_ln, k_ln, yr) = UobsWt(rsp , ksp , r_ln, k_ln);

              // Use hindcast only data
              // FIXME - this is wrong, need to have annual data
              if(yr < nyrs_hind){
                diet_w_sum(rsp, ksp, r_ln, yr) += UobsWt(rsp , ksp , r_ln, k_ln);
              }
            }
          }
        }
      }
    }
    */ 


    // 7.5. Calculate other food stomach content
    of_stomKir.setZero();
    for (yr = 0; yr < nyrs; yr++) {                           // Year loop
      for (rsp = 0; rsp < nspp; rsp++) {                      // Predator species loop
        for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
          for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
            of_stomKir(rsp, r_sex, r_age, yr) = Type( 1 );             // Initialize other suitability
            for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                  of_stomKir(rsp, r_sex, r_age, yr) -= stomKirWt(rsp, ksp, r_sex, k_sex, r_age, k_age, yr);
                }
              }
            }
            if(other_food(rsp) > 0){
              of_stomKir(rsp, r_sex, r_age, yr) /= other_food(rsp); // Penalize this
            }
            if(other_food(rsp) == 0){
              of_stomKir(rsp, r_sex, r_age, yr) = 0;
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

                      suit_tmp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 0;

                      if(yr < nyrs_hind){
                        yr_ind = yr;
                      }
                      if(yr >= nyrs_hind){
                        yr_ind = nyrs_hind - 1;
                      }


                      if(AvgN(ksp, k_sex, k_age, yr) > 0){
                        suit_tmp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = stomKirWt(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) / (AvgN(ksp, k_sex, k_age, yr));

                        // Make into Type 3 MSVPA
                        // U =
                        if(msmMode == 2){
                          suit_tmp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) /= AvgN(ksp, k_sex, k_age, yr);
                        }
                      }


                      // Make sure it is a real number, if not set to 0
                      if(!isFinite(suit_tmp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr))){
                        suit_tmp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 0;
                      }



                      if (wt(pop_wt_index(ksp), k_sex, k_age, yr_ind ) != 0) {
                        stom_div_bio2(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = suit_tmp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) / wt( pop_wt_index(ksp), k_sex, k_age, yr_ind );
                        suma_suit(rsp, r_sex, r_age, yr ) += stom_div_bio2(rsp, ksp, r_sex, k_sex, r_age, k_age, yr); // Calculate sum of stom_div_bio2 across prey and  prey age for each predator, predator age, and year
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
        for (rsp = 0; rsp < nspp; rsp++) {                      // Predator species loop
          for(r_sex = 0; r_sex < nsex(rsp); r_sex++){           // Predator sex
            for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
              suit_other(rsp, r_sex, r_age) = Type( 1 );                 // Initialize other suitability
              for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
                for(k_sex = 0; k_sex < nsex(ksp); k_sex++){       // Prey sex loop
                  for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                    for (yr = 0; yr < nyrs; yr++) {                 // Year loop

                      // Average suitability across years
                      if( suma_suit(rsp, r_sex, r_age, yr ) + of_stomKir(rsp, r_sex, r_age, yr) > 0){
                        suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0) += stom_div_bio2(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) / ( suma_suit(rsp, r_sex, r_age, yr ) + of_stomKir(rsp, r_sex, r_age, yr) );
                      }
                    }       // End year loop

                    // FIXME - Add in interannual variability here
                    suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0) /= nyrs;
                    suit_other(rsp, r_sex, r_age) -= suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0); // Subtract observed suitability from entire suitability (ksp.e. 1)

                    // Remove NAs from crashing populations
                    if(!isFinite(suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0))){
                      suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0) = 0;
                    }

                    // Fill in years
                    for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                      suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0);
                    }
                  }
                }
              }
            }
          }
        }
      } // End Holsman/MSVPA suitability




      // 8.1.2. Length-based GAMMA suitability
      // -- Turned off if not estimating selectivity
      if(suitMode == 1){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
            for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
                  ncnt = 0;
                  gsum = 1.0e-10;                                           // Initialize
                  for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age

                    // if prey are smaller than predator:
                    if (Mn_LatAge(rsp, r_sex, r_age) > Mn_LatAge(ksp, k_sex, k_age)) {
                      x_l_ratio = log(Mn_LatAge(rsp, r_sex, r_age) / Mn_LatAge(ksp, k_sex, k_age));
                      suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                        (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                      ncnt += 1;
                      gsum += exp( suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0) );
                    }
                    else
                      suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0) = 0;
                  }
                  for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                    // if prey are smaller than predator:
                    if (Mn_LatAge(rsp, r_sex, r_age) > Mn_LatAge(ksp, k_sex, k_age)) {
                      suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                    }

                    // Fill in years
                    for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                      suit_main(rsp, ksp, r_age, k_age, yr) = suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0);
                    }
                  }
                }
              }
            }
          }
        }
      } // End GAMMA selectivity


      // 8.1.2. Time-varying length-based GAMMA suitability
      // -- Turned off if not estimating selectivity
      if(suitMode == 2){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
            for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
                  for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                    ncnt = 0;
                    gsum = 1.0e-10;                                           // Initialize
                    for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                      // if prey are smaller than predator:
                      if ( LbyAge( rsp, r_sex, r_age, yr)  > LbyAge( ksp, k_sex, k_age, yr)) {
                        x_l_ratio = log(LbyAge( rsp, r_sex, r_age, yr) / LbyAge( ksp, k_sex, k_age, yr) ); // Log ratio of lengths
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                          (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                        ncnt += 1;
                        gsum += exp( suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) );
                      }
                      else
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = 0;
                    }
                    for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                      // if prey are smaller than predator:
                      if ( LbyAge( rsp, r_sex, r_age, yr)  > LbyAge( ksp, k_sex, k_age, yr)) {
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } // End GAMMA selectivity


      // 8.1.3. Time-varying weight-based GAMMA suitability
      // -- Turned off if not estimating selectivity
      if(suitMode == 3){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
              for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
                  for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                    ncnt = 0;
                    gsum = 1.0e-10;                                           // Initialize

                    if(yr < nyrs_hind){
                      yr_ind = yr;
                    }
                    if(yr >= nyrs_hind){
                      yr_ind = nyrs_hind - 1;
                    }

                    for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                      // if prey are smaller than predator:
                      if (wt(pop_wt_index(rsp), r_sex, r_age, yr_ind ) > wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)) {
                        x_l_ratio = log(wt(pop_wt_index(rsp), r_sex, r_age, yr_ind) / wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // Log ratio of lengths
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                          (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                        ncnt += 1;
                        gsum += exp( suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) );
                      }
                      else{
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = 0;
                      }
                    }
                    for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                      // if prey are smaller than predator:
                      if (wt(pop_wt_index(rsp), r_sex, r_age, yr_ind) > wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)) {
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } // End GAMMA selectivity


      //-------------------------------------------------------------------------------------------------------------
      // 8.1.4. Length-based lognormal suitability
      // -- Turned off if not estimating selectivity
      if(suitMode == 4){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
          for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
            for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
              gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                  for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                    // if prey are smaller than predator:
                    if (Mn_LatAge(rsp, r_sex, r_age) > Mn_LatAge(ksp, k_sex, k_age)) {
                      x_l_ratio = log(Mn_LatAge(rsp, r_sex, r_age) / Mn_LatAge(ksp, k_sex, k_age)); // Log ratio of lengths
                      suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0) = exp(log_phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                      gsum += suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0);
                    }
                    else{
                      suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0) = 0;
                    }
                  }
                }
              }
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                  for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                    suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0) /= gsum;                // Scale, so it sums to 1.

                    // Fill in years
                    for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                      suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, 0);
                    }
                  }
                }
              }
            }
          }
        }
      } // End lognormal selectivity


      // 8.1.5. Time varying length-based lognormal suitability
      if(suitMode == 5){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
          for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
            for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
              for(yr = 0; yr < nyrs; yr++){
                gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
                for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                      // if prey are smaller than predator:
                      if ( LbyAge( rsp, r_sex, r_age, yr)  > LbyAge( ksp, k_sex, k_age, yr)) {
                        x_l_ratio = log(LbyAge( rsp, r_sex, r_age, yr) / LbyAge( ksp, k_sex, k_age, yr) ); // Log ratio of lengths
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = exp(log_phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                        gsum += suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr);
                      }
                      else{
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = 0;
                      }
                    }
                  }
                }
                for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                      suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) /= gsum;                // Scale, so it sums to 1.
                    }
                  }
                }
              }
            }
          }
        }
      } // End lognormal selectivity


      // 8.1.6. Time varying weight-based lognormal suitability
      if(suitMode == 6){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
          for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
            for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
              for(yr = 0; yr < nyrs; yr++){
                gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1

                if(yr < nyrs_hind){
                  yr_ind = yr;
                }
                if(yr >= nyrs_hind){
                  yr_ind = nyrs_hind - 1;
                }

                for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                      // if prey are smaller than predator:
                      if (wt(pop_wt_index(rsp), r_sex, r_age, yr_ind) > wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)) {
                        x_l_ratio = log(wt(pop_wt_index(rsp), r_sex, r_age, yr_ind) / wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // Log ratio of lengths
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = exp(log_phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                        gsum += suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr);
                      }
                      else{
                        suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) = 0;
                      }
                    }
                  }
                }
                for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                      suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr) /= gsum;                // Scale, so it sums to 1.
                    }
                  }
                }
              }
            }
          }
        }
      } // End lognormal selectivity


      // ------------------------------------------------------------------------- //
      // 8. PREDATION MORTALITY EQUATIONS                                          //
      // ------------------------------------------------------------------------- //
      // -- 8.1. HOLSMAN PREDATION MORTALITY

      if ((msmMode == 1) | (msmMode == 2)) {


        // 8.1.3. Calculate available food
        avail_food.setZero();
        othersuit.setZero();
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
                      avail_food(rsp, r_sex, r_age, yr) += suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) * pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) ; // FIXME - include overlap indices: FIXME - mn_wt_stom?
                      othersuit(rsp, r_sex, r_age, yr) += suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr); // FIXME - include overlap indices
                    }
                  }
                }
                if(suitMode == 0){ // Holsman
                  avail_food(rsp, r_sex, r_age, yr) += other_food(rsp) * (Type(1) - (othersuit(rsp, r_sex, r_age, yr)));
                }
                if(suitMode > 0){ // Log-normal and gamma
                  avail_food(rsp, r_sex, r_age, yr) += other_food(rsp) * Type(1);
                }
              }
            }
          }
        }


        // 8.1.3. Calculate predation mortality
        M2.setZero();
        B_eaten.setZero();
        for (ksp = 0; ksp < nspp; ksp++) {                        // Prey species loop
          for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
            for (k_age = 0; k_age < nages(ksp); k_age++) {          // Prey age loop
              for (yr = 0; yr < nyrs; yr++) {                       // Year loop

                if(yr < nyrs_hind){
                  yr_ind = yr;
                }
                if(yr >= nyrs_hind){
                  yr_ind = nyrs_hind - 1;
                }

                for (rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop
                  for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                    for (r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop


                      if(avail_food(rsp, r_sex, r_age, yr) > 0){

                        /*
                         // Predation mortality
                         if(AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) > 0){
                         M2(ksp, k_sex, k_age, yr) += pow(AvgN(ksp, k_sex, k_age, yr) , msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind); // #FIXME - include indices of overlap
                         M2_prop(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind);
                         }
                         */

                        // Type 2 MSVPA
                        if(msmMode == 1){
                          M2(ksp, k_sex, k_age, yr) += (AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr); // #FIXME - include indices of overlap
                          M2_prop(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = (AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr);
                          B_eaten(ksp, k_sex, k_age, yr) += AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr)/ avail_food(rsp, r_sex, r_age, yr);
                          B_eaten_prop(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                        }


                        // Type 3 MSVPA
                        if(msmMode == 2){
                          M2(ksp, k_sex, k_age, yr) += pow(AvgN(ksp, k_sex, k_age, yr) , msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr)
                                                                                                                             * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / (AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind)); // #FIXME - include indices of overlap
                          M2_prop(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * (AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr)) / avail_food(rsp, r_sex, r_age, yr) / (AvgN(ksp, k_sex, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind));
                          B_eaten(ksp, k_sex, k_age, yr) += pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) * AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                          B_eaten_prop(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = pow(AvgN(ksp, k_sex, k_age, yr), msmMode) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind) *  AvgN(rsp, r_sex, r_age, yr) * ration2Age(rsp, r_sex, r_age, yr) * suit_main(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) / avail_food(rsp, r_sex, r_age, yr);
                        }
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
      /* Removing kinzey bits for speed
      if (msmMode > 2) {

        // 8.2.2. Populate other food
        for (rsp = 0; rsp < nspp; rsp++) {
          for (r_age = 0; r_age < nages(rsp); r_age++) {
            Q_other_u(rsp, r_age) = other_food(rsp); // FIXME: may want to estimate other food exp(Q_other_est(rsp,r_age));
          }
        }


        // 8.2.3. Initialize counters
        Type Pred_ratio = 0;          // Predator ratio
        Type Prey_ratio = 0;          // Prey ratio
        Type pred_effect = 0;         // pred_resp * N(r,stage,yr)
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
                    N_pred_eq(rsp, r_sex, r_age) += NByage(rsp, r_sex, r_age, 0) * suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0); // Denominator of Eq. 17 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
                  }
                }

                for (k_age = 0; k_age < nages(ksp); k_age++) {
                  for (r_age = 0; r_age < nages(rsp); r_age++) {
                    N_prey_eq(ksp, k_sex, k_age) += NByage(ksp, k_sex, k_age, 0) * suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, 0); // Denominator of Eq. 16 Kinzey and Punt (2009) 1st year - FIXME: probably should use 2015
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
                      N_pred_yrs(rsp, r_sex, r_age, yr) += NByage(rsp, r_sex, r_age, yr) * suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr); // Numerator of Eq. 17 Kinzey and Punt (2009) 1st year // FIXME: Use averageN?
                    }
                  }


                  for (k_age = 0; k_age < nages(ksp); k_age++) {
                    for (r_age = 0; r_age < nages(rsp); r_age++) {
                      N_prey_yrs(ksp, k_sex, k_age, yr) += NByage(ksp, k_sex, k_age, yr) * suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr); // Numerator of Eq. 16 Kinzey and Punt (2009) 1st year
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
            for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
              for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                for (yr = 0; yr < nyrs; yr++) {                         // Year loop
                  for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop

                    Term = 1.0e-10 + H_1(rsp, ksp) * (Type(1) + H_1a(rsp) * H_1b(rsp) / (Type(r_age) + H_1b(rsp) + Type(1.0e-10))); // Eq. 15 Kinzey and Punt (2009)

                    // Observed species
                    if(ksp < nspp){
                      for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop

                        // Predator-prey ratios
                        Pred_ratio = (N_pred_yrs(rsp, r_sex, r_age, yr) + Type(1.0e-10)) / (N_pred_eq(rsp, r_sex, r_age) + Type(1.0e-10)); // Eq. 17 Kinzey and Punt (2009): Predator biomass relative to equilibrium
                        Prey_ratio = (N_prey_yrs(ksp, k_sex, k_age, yr) + Type(1.0e-10)) / (N_prey_eq(ksp, k_sex, k_age) + Type(1.0e-10)); // Eq. 16 Prey Kinzey and Punt (2009): biomass relative to equilibrium
                        Pred_r(rsp, r_age, yr) = Pred_ratio;
                        Prey_r(ksp, k_age, yr) = Prey_ratio;

                        switch (msmMode) {
                        case 3: // Holling Type I (linear)
                          pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 1.0e-10 + Term;
                          break;
                        case 4: // Holling Type II
                          pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                            ( 1 + H_2(rsp, ksp) * Prey_ratio + 1.0e-10);
                          break;
                        case 5: // Holling Type III
                          pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 1.0e-10 +
                            Term * (1 + H_2(rsp, ksp)) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) /
                              (1 + H_2(rsp, ksp) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) + 1.0e-10 );
                          break;
                        case 6: // predator interference
                          pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                            ( 1 + H_2(rsp, ksp) * Prey_ratio + H_3(rsp, ksp) * (Pred_ratio - 1) + 1.0e-10);
                          break;
                        case 7: // predator preemption
                          pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                            ( (1 + H_2(rsp, ksp) * Prey_ratio) * (1 + H_3(rsp, ksp) * (Pred_ratio - 1)) + Type(1.0e-10));
                          break;
                        case 8: // Hassell-Varley
                          pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 1.0e-10 + Term * (2 + H_2(rsp, ksp) + 1.0e-10) /
                            (1.0 + H_2(rsp, ksp) * Prey_ratio + pow((Prey_ratio + Type(1.0e-10)), H_4(rsp, ksp)) + 1.0e-10 );
                          break;
                        case 9: // Ecosim
                          pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = 1.0e-10 + Term /
                            (1 + H_3(rsp, ksp) * (Pred_ratio - 1 + 1.0e-10));
                          break;
                        default:
                          error("Invalid 'msmMode'");
                        }
                      }
                    }
                    // "other" is linear
                    else{
                      pred_resp(rsp, ksp, r_sex, 0, r_age, 0, yr) = 1.0e-10 + Term;
                    }
                  } // end of r_ages, k_ages loop
                }   // =========================
              }
            }
          }
        }


        // 8.2.7. Mortality as a function of predator age: Eq. 7 Kinzey and Punt (2009)
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
              for (ksp = 0; ksp  <  nspp; ksp++) {
                for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                  for (r_age = 0; r_age < nages(rsp); r_age++) {
                    for (k_age = 0; k_age < nages(ksp); k_age++) {
                      pred_effect = pred_resp(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) * suit_main(rsp , ksp, r_sex, k_sex, r_age, k_age, yr);
                      Vmort_ua(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) = pred_effect * NByage(rsp, r_sex, r_age, yr);
                    }
                  }
                }
              }
            }
          }
        }


        // 8.2.8 Accumulate predation mortality: Eq. 6 Kinzey and Punt (2009)
        Pmort_ua.setZero();
        M2.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp < nspp; rsp++) {
            for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
              for (ksp = 0; ksp < nspp; ksp++) {
                for(k_sex = 0; k_sex < nsex(ksp); k_sex ++){
                  for (k_age = 0; k_age < nages(ksp); k_age++) {
                    for (r_age = 0; r_age < nages(rsp); r_age++) {
                      Pmort_ua(rsp, ksp, r_sex, k_sex, k_age, yr) += Vmort_ua(rsp, ksp, r_sex, k_sex, r_age, k_age, yr);
                    }
                    M2(ksp, k_sex, k_age, yr) += Pmort_ua(rsp, ksp, r_sex, k_sex, k_age, yr);
                  }
                }
              }
            }
          }
        }


        // 8.2.9. Numbers eaten (of modeled prey species); Equations 8 and 9 from Kinzey and Punt 2009
        eaten_la.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (ksp = 0; ksp  <  nspp; ksp++) {
            for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                // Relative number
                NS_Z = NByage(ksp, k_sex, k_age, yr) * (1 - exp(-Zed(ksp, k_sex, k_age, yr) )) / Zed(ksp, k_sex, k_age, yr);                                  // Eq. 8 Kinzey and Punt (2009) Baranov

                for (rsp = 0; rsp  <  nspp; rsp++) {
                  for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
                    for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                      eaten_ua(rsp, ksp, r_age, k_age, yr) = Vmort_ua(rsp, ksp, r_sex, k_sex, r_age, k_age, yr) * NS_Z;                                  // Eq. 8 Kinzey and Punt (2009)

                      for (r_ln = 0; r_ln  <  nlengths(rsp); r_ln++) { // FIXME: Change to diet length bins
                        eaten_la(rsp, ksp, r_ln, k_age, yr) += eaten_ua(rsp, ksp, r_age, k_age, yr)  * age_trans_matrix(pop_age_transition_index(rsp), r_sex, r_age, r_ln ); // Eq. 9 Kinzey and Punt (2009)
                      }
                    }
                  }
                }
              }
            }
          }
        }


        // 8.2.10. Mass eaten (including "other")
        Q_mass_l.setZero();
        Q_mass_u.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          if(yr < nyrs_hind){
            yr_ind = yr;
          }
          if(yr >= nyrs_hind){
            yr_ind = nyrs_hind - 1;
          }

          for (rsp = 0; rsp  < nspp; rsp++) {
            for (ksp = 0; ksp  <  nspp + 1; ksp++) {


              // Species included
              if (ksp < nspp) {
                // Results by length: Eqn 10 kinda Kinzey and Punt (2009)
                for (r_ln = 0; r_ln  <  nlengths(rsp); r_ln++) { // FIXME: Change to diet bins
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {
                      Q_mass_l(rsp, ksp, r_ln, yr) += eaten_la(rsp, ksp, r_ln, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind);
                    }
                  }
                }
                // Results by k_age: Eq. 11 Kinzey and Punt (2009)
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
                    for (k_age = 0; k_age < nages(ksp); k_age++) {
                      Q_mass_u(rsp, ksp, r_age, yr) += eaten_ua(rsp, ksp, r_age, k_age, yr) * wt(pop_wt_index(ksp), k_sex, k_age, yr_ind);
                    }
                  }
                }
              }

              // Other food
              else {
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  Tmort = 0;
                  for(r_sex = 0; r_sex < nsex(rsp); r_sex++){
                    pred_effect = pred_resp(rsp, ksp, r_sex, 0, r_age, 0, yr); // Other food k_age = 0 because 1 age is used
                    Tmort += pred_effect * NByage(rsp, r_sex, r_age, yr);
                  }
                  // FIXME: add sex specific diet comp if desired. Move Q_mass u to previous sex loop
                  Q_mass_u(rsp, ksp, r_age, yr)  = Q_other_u(rsp, r_age) * (Type(1) - exp(-Tmort));                     // Eq. 11 Kinzey and Punt (2009)
                }
                for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) { // FIXME: Change to diet bins
                  for (r_age = 0; r_age < nages(rsp); r_age++) {
                    Q_mass_l(rsp, ksp, r_ln, yr) += Q_mass_u(rsp, ksp, r_age, yr) * age_trans_matrix(pop_age_transition_index(rsp), r_sex, r_age, r_ln); // Eq. 10 Kinzey and Punt (2009)
                  }
                }
              }
            }
          }
        }


        // 8.2.11. Total up the consumption by each predator and normalize (Eqn 15)
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) { // FIXME: Change to diet length bins
              Q_ksum_l = 0;
              for (ksp = 0; ksp < (nspp + 1); ksp++) {
                Q_ksum_l += Q_mass_l(rsp, ksp, r_ln, yr) + 1.0e-10; // Table 3 - Diet weight proportions Kinzey and Punt (2009)
              }
              for (ksp = 0; ksp < (nspp + 1); ksp++) {
                Q_hat(rsp, ksp, r_ln, yr) = (1.0e-10 + Q_mass_l(rsp, ksp, r_ln, yr) / Q_ksum_l); // changed paranthesis -dhk apr 28 09
              }
            }
          }
        }


        // 8.2.12. Predict ration Eq. 13 Kinzey and Punt (2009)
        Type n_avg, numer, denom;

        omega_hat_ave.setZero();
        omega_hat.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {
          for(r_sex = 0; r_sex < nsex(rsp); r_sex ++){
            for (r_age = 0; r_age < nages(rsp); r_age++) {
              // Calculate year-specific values
              numer = 0;
              denom = 0;
              for (yr = 0; yr < nyrs; yr++) {
                // Average abundance
                n_avg = 1.0e-10 + NByage(rsp, r_sex, r_age, yr) * sqrt(S(rsp, r_sex, r_age, yr));  // FIXME: Add average N added 1.0e-10 -dhk June 24 08. not in NRM tpl -dhk apr 28 09
                // find total consumption by this r_age-class
                for (ksp = 0; ksp < (nspp + 1); ksp++) {
                  omega_hat(rsp, r_age, yr) += Q_mass_u(rsp, ksp, r_age, yr);
                }
                numer += omega_hat(rsp, r_age, yr) / + Type(1.0e-10); // FIXME: Divide by 365 to make into daily ration
                denom += n_avg;
                // normalize
                omega_hat(rsp, r_age, yr) /= (n_avg); // FIXME: Divide by 365 to make into daily ration
              }
              omega_hat_ave(rsp, r_age) = numer / denom;
            }
          }
        }
        // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
      } // End 8.2. Kinzey predation
      */
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

      srv_bio_hat(srv_ind) = srv_q(srv, yr_ind) * pow(srv_bio_hat(srv_ind), (1 + srv_q_pow(srv)));

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
            fsh_bio_hat(fsh_ind) += F(flt, sex, age, flt_yr) / Zed(sp, sex, age, flt_yr) * (1 - exp(-Zed(sp, sex, age, flt_yr))) * NByage(sp, sex, age, flt_yr) * wt( flt_wt_index(flt), sex, age, yr_ind ); // 5.5.
          }

          // By numbers
          if(flt_units(flt) == 2){
            fsh_bio_hat(fsh_ind) += F(flt, sex, age, flt_yr) / Zed(sp, sex, age, flt_yr) * (1 - exp(-Zed(sp, sex, age, flt_yr))) * NByage(sp, sex, age, flt_yr); // 5.5.
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


        // Fishery
        if(flt_type(flt) == 1){
          // Sexes combined or 1 sex assessment
          if((flt_sex == 0) | (nsex(sp) == 1)){
            for(sex = 0; sex < nsex(sp); sex ++){
              age_hat(comp_ind, age ) += F(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr); // 5.4.
              n_hat( comp_ind ) += F(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
            }
          }

          // Sex specific composition data
          if((flt_sex == 1)| (flt_sex == 2)){
            sex = flt_sex - 1;
            age_hat(comp_ind, age ) = F(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
            n_hat( comp_ind ) += F(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
          }

          // Joint composition data
          if(flt_sex == 3){
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = F(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
              // Total numbers
              n_hat(comp_ind) += F(flt, sex, age, yr) / Zed(sp, sex, age, yr) * (1 - exp(-Zed(sp, sex, age, yr))) * NByage(sp, sex, age, yr);
            }
          }
        }


        // Survey
        if(flt_type(flt) == 2){
          // Sexes combined or 1 sex assessment
          if((flt_sex == 0) | (nsex(sp) == 1)){
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age ) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
              // Total numbers
              n_hat(comp_ind) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));   // Total numbers
              // FIXME - put power in for catchability
            }
          }

          // Sex specific composition data
          if((flt_sex == 1) | (flt_sex == 2)){
            sex = flt_sex - 1;
            // Survey catch-at-age
            age_hat(comp_ind, age ) = NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
            // Total numbers
            n_hat(comp_ind) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
          }

          // Joint composition data
          if(flt_sex == 3){
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
              // Total numbers
              n_hat(comp_ind) += NByage(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * srv_q(flt, yr_ind) * exp( - Type(mo/12) * Zed(sp, sex, age, yr));
            }
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
        true_age_comp_hat(comp_ind, age ) = age_hat(comp_ind, age ) / n_hat(comp_ind);
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
          comp_hat(comp_ind, age) = age_obs_hat(comp_ind, age ) / n_hat(comp_ind);
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
          comp_hat(comp_ind, ln ) = comp_hat(comp_ind, ln) / n_hat(comp_ind);
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
      r_sexes(stom_ind, 0) = 0; r_sexes(stom_ind, 1) = 1;
      k_sexes(stom_ind, 0) = 0; k_sexes(stom_ind, 1) = 1;

      if(r_sex > 0){
        r_sexes(stom_ind, 0) = r_sex - 1;  r_sexes(stom_ind, 1) = r_sex - 1;
      }

      if(k_sex > 0){
        k_sexes(stom_ind, 0) = k_sex - 1;  k_sexes(stom_ind, 1) = k_sex - 1;
      }

      // Initialize
      UobsWtAge_hat(stom_ind, 1) = 0;

      for(int j = 0; j < 2; j ++){
        for(int k = 0; k < 2; k ++){

          if(flt_yr > 0){
            yr = flt_yr - styr;

            if(yr < nyrs_hind){
              UobsWtAge_hat(stom_ind, 1) += (AvgN(ksp, k_sexes(stom_ind, k), k_age, yr) * suit_main(rsp, ksp, r_sexes(stom_ind, j), k_sexes(stom_ind, k), r_age, k_age, yr) * wt(pop_wt_index(ksp), k_sexes(stom_ind, k), k_age, yr)) / avail_food(rsp, r_sexes(stom_ind, j), r_age, yr) / 4 ; // NOTE: Divide by species
            }
          }

          // Average of years
          if(flt_yr == 0){
            for (yr = 0; yr < nyrs_hind; yr++) {
              UobsWtAge_hat(stom_ind, 1) += (AvgN(ksp, k_sexes(stom_ind, k), k_age, yr) * suit_main(rsp, ksp, r_sexes(stom_ind, j), k_sexes(stom_ind, k), r_age, k_age, yr) * wt(pop_wt_index(ksp), k_sexes(stom_ind, k), k_age, yr)) / avail_food(rsp, r_sexes(stom_ind, j), r_age, yr) / 4 / nyrs_hind;
            }
          }
        }
      }
    }
    // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
  } // End population dynamics iterations
  // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //


  // ------------------------------------------------------------------------- //
  // 11. LIKELIHOOD EQUATIONS

  // 11.0. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(19, n_flt); jnll_comp.setZero();  // matrix of negative log-likelihood components

  // -- Data likelihood components
  // Slot 0 -- Survey biomass
  // Slot 1 -- Total catch (kg)
  // Slot 2 -- Age/length composition
  // Slot 3 -- Sex ratio likelihood
  // Slot 4 -- Selectivity
  // Slot 5 -- Selectivity random walk deviates
  // Slot 6 -- Selectivity random effect deviates
  // Slot 7 -- Survey selectivity normalization
  // Slot 8 -- Survey catchability random walk deviates
  // Slot 9 -- Survey catchability random effect deviates
  // -- Priors/penalties
  // Slot 10 -- Tau -- Annual recruitment deviation
  // Slot 11 -- init_dev -- Initial abundance-at-age
  // Slot 12 -- Epsilon -- Annual fishing mortality deviation
  // Slot 13 -- SPR penalities
  // Slot 14 -- N-at-age < 0 penalty
  // -- M2 likelihood components
  // Slot 15 -- Ration likelihood
  // Slot 16 -- Ration penalties
  // Slot 17 -- Diet weight likelihood
  // Slot 18 -- Stomach content of prey length ln in predator length a likelihood


  // 11.1. OFFSETS AND PENALTIES
  // 11.1.1 -- Set up offset objects
  vector<Type> offset(n_flt); offset.setZero(); // Offset for multinomial likelihood
  vector<Type> offset_diet_w(nspp); offset_diet_w.setZero(); // Offset for total stomach content weight likelihood
  vector<Type> offset_diet_l(nspp); offset_diet_l.setZero(); // Offset for total stomach content of prey length ln in predator length a
  vector<Type> offset_uobsagewt(nspp); offset_uobsagewt.setZero(); // Offset for stomach proportion by weight likelihood


  // 11.1.1. -- Age/length comp offsets
  for (comp_ind = 0; comp_ind < comp_obs.rows(); comp_ind++) {

    flt = comp_ctl(comp_ind, 0) - 1;            // Temporary fishery index
    sp = comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
    flt_sex = comp_ctl(comp_ind, 2);            // Temporary index for comp sex (0 = combined, 1 = female, 2 = male)
    comp_type = comp_ctl(comp_ind, 3);          // Temporary index for comp type (0 = age, 1 = length)
    yr = comp_ctl(comp_ind, 4) - styr;          // Temporary index for years of data

    if(yr > 0){
      // Adjustment for joint sex composition data
      Type joint_adjust = 1;
      if(flt_sex == 3){
        joint_adjust = 2;
      }

      // Number of ages/lengths
      Type n_comp = 0;
      if(comp_type == 0){
        n_comp = nages(sp) * joint_adjust;
      }
      if(comp_type == 1){
        n_comp = nlengths(sp) * joint_adjust;
      }

      // Add years from hindcast only
      if(yr <= endyr){
        for (ln = 0; ln < n_comp; ln++) {
          if(!isNA( comp_obs(comp_ind, ln ))){
            if(comp_obs(comp_ind, ln) > 0){
              offset(flt) -= Type(comp_n(comp_ind, 1)) * (comp_obs(comp_ind, ln)) * log(comp_obs(comp_ind, ln)) ;
            }
          }
        }
      }
    }
  }

  // 11.1.4. -- Offsets for diet (weights)
  /* Turning off bits for kinzey
  for (rsp = 0; rsp < nspp; rsp++) {
    for (yr = 0; yr < nyrs_hind; yr++) {
      for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) {
        // if (stoms_w_N(rsp,r_ln,yr) > 0){ Sample size
        for (ksp = 0; ksp < nspp; ksp++) { // FIXME: no offset for "other food"
          if ( diet_w_sum(rsp, ksp, r_ln, yr) > 0) {

            offset_diet_w(rsp) -= stom_tau * (diet_w_sum(rsp, ksp, r_ln, yr)) * log(diet_w_sum(rsp, ksp, r_ln, yr)); // FIXME add actual sample size
          }
        }
      }
    }
  }
  */

  // 11.1.5. -- Offsets for diet (lengths) FIXME: add real stomach sample size
  for (rsp = 0; rsp < nspp; rsp++) {
    for (ksp = 0; ksp < nspp; ksp++) {
      for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) {
        for (k_ln = 0; k_ln < nlengths(ksp); k_ln++) {
          if (Uobs(rsp, ksp, r_ln, k_ln) > 0) {
            offset_diet_l(rsp) -= stom_tau * Uobs(rsp, ksp, r_ln, k_ln) * log(Uobs(rsp, ksp, r_ln, k_ln) + 1.0e-10);
          }
        }
      }
    }
  }

  // 11.1.6. -- Offsets for diet proportion by weight
  // FIXME - change years to adjust for missing years of diet data
  // If 5D array
  for(int stom_ind = 0; stom_ind < UobsWtAge.rows(); stom_ind++){
    rsp = UobsWtAge_ctl(stom_ind, 0) - 1; // Index of pred

    if(UobsWtAge(stom_ind, 1) > 0){
      offset_uobsagewt(rsp) -= (UobsWtAge(stom_ind, 0) * (UobsWtAge(stom_ind, 1)) * log(UobsWtAge(stom_ind, 1)));
    }
  }


  // 11.2. FIT OBJECTIVE FUNCTION
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
        n_comp = nages(sp) * joint_adjust(comp_ind);

        // Accumulation age
        for(age = 0; age < n_comp; age++){
          int age_test = age;
          if(age >= nages(sp)){
            age_test = age - nages(sp);
          }

          // Lower
          if(age_test < flt_accum_age_lower(flt)){

            comp_obs(comp_ind, flt_accum_age_lower(flt)) += comp_obs(comp_ind, age);
            comp_obs(comp_ind, age) = 0;

            comp_hat(comp_ind, flt_accum_age_lower(flt)) += comp_hat(comp_ind, age);
            comp_hat(comp_ind, age) = 0;
          }
          // Upper
          if(age_test > flt_accum_age_upper(flt)){
            comp_obs(comp_ind, flt_accum_age_upper(flt)) += comp_obs(comp_ind, age);
            comp_obs(comp_ind, age) = 0;

            comp_hat(comp_ind, flt_accum_age_upper(flt)) += comp_hat(comp_ind, age);
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
                jnll_comp(2, flt) -= comp_weights(flt) * Type(comp_n(comp_ind, 1)) * (comp_obs(comp_ind, ln) + 0.00001) * log(comp_hat(comp_ind, ln) + 0.00001) ;
              }
            }
          }
        }
      }
    }
  }


  // Remove offsets
  for (flt = 0; flt < n_flt; flt++) {
    if(flt_type(flt) > 0){
      jnll_comp(2, flt) -= offset(flt);
    }
  }


  // Slot 3 -- Observed sex ratio
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs_hind; yr++){
      if((nsex(sp) == 2) & (estDynamics(sp) == 0)){
        if(est_sex_ratio(sp) == 1){
          jnll_comp(4, sp) -= dnorm(sex_ratio_mean_hat(sp, yr), sex_ratio(sp, 1), sex_ratio_sigma(sp)); // Using the 2nd age here, cause recruitment is assumed to be age 1 (c++ 0)
        }

        if(est_sex_ratio(sp) == 2){
          for(age = 1; age < nages(sp); age++){ // Start at age 2 because age 1 is fixed
           jnll_comp(4, sp) -= dnorm(sex_ratio_hat(sp, age, yr), sex_ratio(sp, age), sex_ratio_sigma(sp));
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


    // Penalized likelihood time-varying selectivity deviates
    if((flt_varying_sel(flt) == 1) & (flt_sel_type(flt) != 2) & (flt_type(flt) > 0)){
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

    // Penalize time-varying selectivity deviates (Random effects)
    if((flt_varying_sel(flt) == 2) & (flt_sel_type(flt) != 2) & (flt_type(flt) > 0)){
      for(sex = 0; sex < nsex(sp); sex ++){
        for(yr = 0; yr < nyrs_hind; yr++){

          // Logistic deviates
          jnll_comp(6, flt) -= dnorm(sel_inf_dev_re(0, flt, sex, yr), Type(0.0), sigma_sel(flt), true);
          jnll_comp(6, flt) -= dnorm(ln_sel_slp_dev_re(0, flt, sex, yr), Type(0.0), sigma_sel(flt), true);

          // Double logistic deviates
          if(flt_sel_type(flt) == 3){
            jnll_comp(6, flt) -= dnorm(sel_inf_dev_re(1, flt, sex, yr), Type(0.0), sigma_sel(flt), true);
            jnll_comp(6, flt) -= dnorm(ln_sel_slp_dev_re(1, flt, sex, yr), Type(0.0), sigma_sel(flt), true);
          }
        }
      }
    }


    // Random walk: Type 4 = random walk on ascending and descending for double logistic; Type 5 = ascending only for double logistics
    if(((flt_varying_sel(flt) == 4)|(flt_varying_sel(flt) == 5)) & (flt_sel_type(flt) != 2) & (flt_type(flt) > 0)){
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


        jnll_comp(7, flt) += square(sel_slp_dev1_sum) * 10000;
        jnll_comp(7, flt) += square(sel_inf_dev1_sum) * 10000;
        jnll_comp(7, flt) += square(sel_slp_dev2_sum) * 10000;
        jnll_comp(7, flt) += square(sel_inf_dev2_sum) * 10000;
      }
    }
  } // End selectivity loop



  // Slot 7 -- Add survey selectivity normalization
  for(flt = 0; flt < n_flt; flt++){
    sp = flt_spp(flt);
    for(sex = 0; sex < nsex(sp); sex++){
      if(flt_type(flt) > 0){
        jnll_comp(7, flt) += 50 * square(avgsel(flt, sex));
      }
    }
  }


  // Slot 8-9 -- Survey catchability deviates
  for(flt = 0; flt < n_flt; flt++){

    // Prior on catchability
    if( est_srv_q(flt) == 2){
      jnll_comp(8, flt) -= dnorm(ln_srv_q(flt), ln_srv_q_prior(flt), sigma_srv_q(flt), true);
    }

    // Penalized likelihood
    if((srv_varying_q(flt) == 1) & (flt_type(flt) == 2)){
      for(yr = 0; yr < nyrs_hind; yr++){
        jnll_comp(8, flt) -= dnorm(ln_srv_q_dev(flt, yr), Type(0.0), time_varying_sigma_srv_q(flt), true );
      }
    }

    // Random effects
    if((srv_varying_q(flt) == 2) & (flt_type(flt) == 2)){
      for(yr = 0; yr < nyrs_hind; yr++){
        jnll_comp(9, flt) -= dnorm(ln_srv_q_dev_re(flt, yr), Type(0.0), time_varying_sigma_srv_q(flt), true );
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


  // Slots 9-11 -- PRIORS: PUT RANDOM EFFECTS SWITCH HERE
  for (sp = 0; sp < nspp; sp++) {
    // Slot 10 -- init_dev -- Initial abundance-at-age
    for (age = 1; age < nages(sp); age++) {
      jnll_comp(11, sp) -= dnorm( init_dev(sp, age - 1) - square(r_sigma(sp)) / 2, Type(0.0), r_sigma(sp), true);
    }

    for (yr = 0; yr < nyrs_hind; yr++) {
      // Slot 11 -- Tau -- Annual recruitment deviation
      jnll_comp(10, sp) -= dnorm( rec_dev(sp, yr)  - square(r_sigma(sp)) / 2, Type(0.0), r_sigma(sp), true);    // Recruitment deviation using random effects.
    }
  }


  // Slot 13 -- SPR reference point penalties
  for (sp = 0; sp < nspp; sp++) {
    if((msmMode == 0) & (proj_F_prop.sum() > 0)){
      jnll_comp(13, sp)  += 200*square((SB35(sp)/SB0(sp))-0.35);
      jnll_comp(13, sp)  += 200*square((SB40(sp)/SB0(sp))-0.40);
    }
  }


  // Slot 14 -- N-at-age < 0 penalty. See posfun
  for (sp = 0; sp < nspp; sp++){
    if(estDynamics(sp) == 0){
      jnll_comp(14, sp) += zero_pop_pen(sp);
    }
  }


  // 11.3. Diet likelihood components from MSVPA
  if ((msmMode == 1) & (suitMode > 0)) {
    // Slot 15 -- Diet weight likelihood
    for(int stom_ind = 0; stom_ind < UobsWtAge.rows(); stom_ind++){

      rsp = UobsWtAge_ctl(stom_ind, 0) - 1; // Index of pred
      if(UobsWtAge_hat(stom_ind, 1) > 0){
        jnll_comp(17, rsp) -= Type(UobsWtAge(stom_ind, 0)) * (UobsWtAge(stom_ind, 1)) * log(UobsWtAge_hat(stom_ind, 1));
      }
    }

    // Remove offset
    for (rsp = 0; rsp < nspp; rsp++) {
      jnll_comp(17, rsp) -= offset_uobsagewt(rsp);
    }

  } // End diet proportion by weight component

  // 11.4. Diet likelihood components from Kinzey and Punt
/* Turning off Kinzey bits for speed
  if ((msmMode > 2)) {
    // Slot 15 -- Ration likelihood
    for (yr = 0; yr < nyrs_hind; yr++) {
      for (sp = 0; sp < nspp; sp++) {
        for(sex = 0; sex < nsex(sp); sex ++){
          for (age = 1; age < nages(sp); age++) { // don't include age zero in likelihood
            jnll_comp(15, sp) += 0.5 * pow( log( omega_hat(sp, age, yr) + 1.0e-10) -
              log(ration2Age(sp, sex, age, yr)), 2) / (sd_ration * sd_ration); // FIXME: add year indices for ration
          }
        }
      }
    }


    // Slot 16 -- Ration penalalties FIXME: not sure if necessary: talk to Andre
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        Type mean_ohat = 0;
        for (yr = 0; yr < nyrs; yr++) {
          mean_ohat += omega_hat(sp, age, yr) / nyrs;
        }
        for (yr = 0; yr < nyrs; yr++) {
          jnll_comp(16, sp) += 20 *  pow(omega_hat(sp, age, yr) - mean_ohat, 2);
        }
      }
    }


    // Slot 15 -- Diet weight likelihood
    for (rsp = 0; rsp < nspp; rsp++) {
      for (yr = 0; yr < nyrs_hind; yr++) {
        for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) {
          for (ksp = 0; ksp < nspp; ksp++) {                  // FIME: need to add in other food
            if (diet_w_sum(rsp, ksp, r_ln, yr) > 0) { // (rsp, ksp, a, r_ln, yr)
              jnll_comp(17, rsp) -= stom_tau * diet_w_sum(rsp, ksp, r_ln, yr) * log(Q_hat(rsp, ksp, r_ln, yr) + 1.0e-10); // FIXME: Q_hat has some NAs
            }
          }
        }
      }
      jnll_comp(17, rsp) -= offset_diet_w(rsp);
    }


    Type Denom = 0;

    // Calculate the predicted fraction by length-class (Eqn 17)
    T_hat.setZero();
    for (rsp = 0; rsp < nspp; rsp++) {
      for (ksp = 0; ksp < nspp; ksp++) {
        for(k_sex = 0; k_sex < nsex(ksp); k_sex++){
          vector<Type> eaten_lmy(nlengths(ksp)); eaten_lmy.setZero(); // no. of prey@length eaten by a predator length during iyr

          for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) { //
            // This is Equation 17
            for (yr = 0; yr < nyrs_hind; yr++) { // FIXME: loop by stomach year
              // int stom_yr = yrs_stomlns(rsp, ksp, stm_yr); // FIXME: index by stomach year
              for (k_ln = 0; k_ln < nlengths(ksp); k_ln ++) {
                for (k_age = 0; k_age < nages(ksp); k_age++) {
                  eaten_lmy(k_ln) += eaten_la(rsp, ksp, r_ln, k_age, yr) * age_trans_matrix(pop_age_transition_index(ksp), k_sex, k_age, k_ln);
                }
                T_hat(rsp, ksp, r_ln, k_ln) += stom_tau * eaten_lmy(k_ln); // FIXME: add actual stomach content sample size
              }
            }


            Denom = 0;
            for (k_ln = 0; k_ln < nlengths(ksp); k_ln++) {
              Denom += T_hat(rsp, ksp, r_ln, k_ln);
            }

            // Renormalize the eaten vector
            for (k_ln = 0; k_ln < nlengths(ksp); k_ln++) {
              if(Denom > 0){
                T_hat(rsp, ksp, r_ln, k_ln) /= Denom;
              }

              // Likelihood of diet length         / This is equation 16
              if (Uobs(rsp, ksp, r_ln, k_ln) > 0) {
                jnll_comp(18, rsp) -= stom_tau * Uobs(rsp, ksp, r_ln, k_ln) * log(T_hat(rsp, ksp, r_ln, k_ln)  + 1.0e-10);
              }
            }
          }
        }
      }
      jnll_comp(18, rsp) -= offset_diet_l(rsp);
    }
  } // End if statement for Kinzey diet likelihood
  */

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
  REPORT( mn_rec );
  REPORT( pmature );
  REPORT( Zed );
  REPORT( NByage );
  REPORT( AvgN );
  REPORT( sex_ratio_hat );
  REPORT( sex_ratio_mean_hat );
  REPORT( S );
  REPORT( biomassByage );
  REPORT( biomassSSBByage );
  REPORT( biomass );
  REPORT( biomassSSB );
  REPORT(r_sigma);
  REPORT( R_sexr );
  REPORT( R );
  REPORT( M );
  REPORT( NbyageSPR);
  REPORT( SB0 );
  REPORT( SB35 );
  REPORT( SB40 );
  ADREPORT( Zed );
  ADREPORT( biomass );
  ADREPORT( biomassSSB );
  ADREPORT( r_sigma );
  ADREPORT( R );
  ADREPORT( pop_scalar );


  // -- 12.2. Selectivity
  REPORT( sel );
  REPORT( avgsel );
  REPORT( emp_sel_obs );
  REPORT( sel_tmp );
  REPORT( sigma_sel );


  // -- 12.3. Survey components
  REPORT( srv_bio_hat );
  REPORT( srv_log_sd_hat );
  REPORT( sigma_srv_index );
  REPORT( srv_q );
  REPORT( srv_q_analytical );
  REPORT( sigma_srv_q );
  REPORT( time_varying_sigma_srv_q );
  REPORT( ln_srv_q_dev );


  // -- 12.4. Fishery components
  REPORT( F );
  REPORT( F_tot );
  REPORT( fsh_bio_hat );
  REPORT( fsh_log_sd_hat );
  REPORT( proj_FABC );
  REPORT( FSPR );
  REPORT( F35_tot );
  REPORT( F40_tot );


  // 12.5. Age/length composition
  REPORT( age_obs_hat );
  REPORT( age_hat );
  REPORT( comp_obs );
  REPORT( comp_hat );
  REPORT( true_age_comp_hat );
  REPORT( n_hat );
  REPORT( comp_n );


  // -- 12.6. Likelihood components
  REPORT( jnll_comp );
  REPORT( offset );
  REPORT( offset_diet_w );
  REPORT( offset_diet_l );
  REPORT( offset_uobsagewt );


  // -- 12.7. Ration components
  REPORT( ConsumAge );
  REPORT( LbyAge );
  REPORT( mnWt_obs );
  REPORT( fT );
  REPORT( env_index_hat );
  REPORT( ration2Age );


  // 12.8. Suitability components
  REPORT( suma_suit );
  REPORT( suit_main );
  REPORT( suit_other );
  REPORT( stom_div_bio2 );
  REPORT( stomKir );
  REPORT( stomKirWt );
  REPORT( avail_food );
  REPORT( othersuit );
  REPORT( of_stomKir );
  REPORT( M1 );
  REPORT( M2 );
  REPORT( M2_prop );
  REPORT( B_eaten );
  REPORT( B_eaten_prop );
  REPORT( UobsAge_hat );
  REPORT( UobsWtAge_hat );

  // -- 12.9. Kinzey predation functions
  /* Turning off kinzey bits for speed
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
  REPORT( Vmort_ua );
  REPORT( eaten_la );
  REPORT( eaten_ua );

  REPORT( Q_mass_l );
  REPORT( Q_mass_u );
  REPORT( Q_other_u );
  REPORT( Q_hat );
  REPORT( T_hat );

  REPORT( omega_hat );
  REPORT( omega_hat_ave );
  */


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
  if (debug == 0) {
    jnll = jnll_comp.sum();
  }

  // Debug mode
  if (debug > 0) {
    jnll = dummy * dummy;
  }


  REPORT( jnll );
  return jnll;
}
