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
//
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
//  -- 8.1. Holsman et al Predation mortality
//  -- 9.2. Kinzey and Punt predation mortality
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

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
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
matrix<Type> matrix_from_array(array<Type> a1, int sheet){
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
  DATA_INTEGER(msmMode);
  //    0 = run in single species mode
  //    1 = run in MSM mode  Holsman et al (2015) MSVPA based
  // DATA_INTEGER(est_diet);              // Include diet data in the likelihood
  DATA_INTEGER(suitMode);                 // Estimate suitability
  DATA_INTEGER(avgnMode);                 // N used for predation function
  //    0 = AvgN
  //    1 = N*exp(-Z / 2))
  //    2 = N

  DATA_INTEGER( random_rec );             // Logical of whether to treate recruitment deviations as random effects
  DATA_INTEGER( niter );                  // Number of loops for MSM mode


  // Columns =  Species, Survey, Selectivity, Estimate_q, log_q_prior
  // Selectivity 0 = empirical, 1 = logistic, 2 = non-parametric, 3 = double logistic

  // 1.2. Temporal dimensions
  DATA_INTEGER( styr );                      // Start year
  DATA_INTEGER( endyr );                     // End of estimation years
  DATA_INTEGER( projyr );                    // End year of projection
  int nyrs = projyr - styr + 1;
  int nyrs_hind = endyr - styr + 1;
  // DATA_INTEGER( endyr);                   // End year of projection

  // 1.3. Number of species
  DATA_INTEGER( nspp );                   // Number of species (prey)
  DATA_IVECTOR( pop_wt_index );           // Dim 3 of wt to use for population dynamics
  DATA_IVECTOR( pop_alk_index );          // Dim 3 of wt to use for age_transition_matrix
  pop_wt_index = pop_wt_index - 1;
  pop_alk_index -= 1;


  // 1.4. MODEL OBJECTS
  // 1.4.1. LOOPING INDICES -- k = observation, sp = species, age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length
  int sp, age, ln, ksp, k_age, k_ln, yr, rsp, r_age, r_ln; // k
  int srv, fsh, sex;                                                    // Survey and fishery indices
  int fsh_yr, srv_yr, srv_comp_type, fsh_comp_type;                     // Year indices
  int fsh_ind, srv_ind, comp_ind, ctl_ind;                              // Indices for survey data sets
  int ncnt;    // Pointers
  Type mo = 0;                                                          // Month float
  if (msmMode == 0) { niter = 1; }                                      // Number of iterations for SS mode

  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                           //
  // ------------------------------------------------------------------------- //

  // 2.1. FIXED VALUES
  Type MNConst = 0.001;                   // Constant additive for logistic functions
  Type curv_pen_fsh = 12.5;               // Fishery selectivity penalty
  Type sd_ration = 0.05;                  // SD of ration likelihood


  // 2.2. DIMENSIONS
  // int nyrs_proj = proj_yr - styr + 1;        // End year

  // -- 2.2.2. Species attributes
  DATA_IVECTOR( nages );                  // Number of species (prey) ages; n = [1, nspp]
  DATA_IVECTOR( nlengths );               // Number of species (prey) lengths; n = [1, nspp]
  DATA_VECTOR( stom_tau );               // Stomach sample size
  int max_age = imax(nages);              // Integer of maximum nages to make the arrays; n = [1]
  DATA_ARRAY( age_error );                // Array of aging error matrices for each species; n = [nspp, nages, nages]

  // 2.3. DATA INPUTS (i.e. assign data to objects)
  // -- 2.3.1 Fishery Components

  DATA_IMATRIX( fsh_control );            // Survey specifications
  DATA_IMATRIX( fsh_biom_ctl );           // Info for fishery biomass; columns = Fishery_name, Fishery_code, Species, Year
  DATA_IMATRIX( fsh_biom_n );             // Info for fishery biomass; columns = Month
  DATA_MATRIX( fsh_biom_obs );            // Observed fishery catch biomass (kg) and cv; n = [nobs_fish_biom, 2]; columns = Observation, Error
  DATA_IMATRIX( fsh_emp_sel_ctl );        // Info on empirical fishery selectivity; columns =  Fishery_name, Fishery_code, Species, Year

  DATA_MATRIX( fsh_emp_sel_obs );         // Observed emprical fishery selectivity; columns = Compe_1, Comp_2, etc.
  DATA_IMATRIX( fsh_comp_ctl );           // Info on observed fishery age/length comp; columns = Fishery_name, Fishery_code, Species, Year
  DATA_MATRIX( fsh_comp_n );              // Info on month and sample size of observed fishery age/length comp; columns = Month, Sample size
  DATA_MATRIX( fsh_comp_obs );            // Observed fishery age/length comp; cols = Comp_1, Comp_2, etc. can be proportion

  // -- 2.3.2 Survey components
  DATA_IMATRIX( srv_control );            // Survey specifications
  DATA_IMATRIX( srv_biom_ctl );           // Info for survey biomass; columns = Survey_name, Survey_code, Species, Year
  DATA_MATRIX( srv_biom_n );              // Info for survey biomass; columns = Month

  DATA_MATRIX( srv_biom_obs );            // Observed survey biomass (kg) and cv; n = [nobs_srv_biom, 2]; columns = Observation, Error
  DATA_IMATRIX( srv_emp_sel_ctl );        // Info on empirical survey selectivity; columns =  Survey_name, Survey_code, Species, Year
  DATA_MATRIX( srv_emp_sel_obs );         // Observed emprical survey selectivity; columns = Compe_1, Comp_2, etc.
  DATA_IMATRIX( srv_comp_ctl );           // Info on observed survey age/length comp; columns = Survey_name, Survey_code, Species, Year
  DATA_MATRIX( srv_comp_n );             // Month and sample size on observed survey age/length comp; columns = Month, Sample size
  DATA_MATRIX( srv_comp_obs );            // Observed survey age/length comp; cols = Comp_1, Comp_2, etc. can be proportion


  DATA_ARRAY( age_trans_matrix);          // observed sp_age/size compositions; n = [nspp, nages, srv_age_bins]

  // -- 2.3.5. Weight-at-age
  DATA_ARRAY( wt );                       // Weight-at-age by year; n = [nyrs, nages, nspp]: FIXME: Change nyrs to nyrs_wt_at_age if data don't span entire bit

  // 2.3.6. Diet data
  DATA_VECTOR( fday );                    // number of foraging days for each predator; n = [1, nspp] #FIXME - assuming this is the same as fdays
  DATA_ARRAY( Pyrs );                     // n = [nspp, nyrs+1, nages]: #FIXME - Assuming this is the same as Pby_yr?
  DATA_ARRAY( Uobs );                     // pred, prey, predL, preyL U matrix (mean number of prey in each pred); n = [nspp, nspp, maxL, maxL]
  DATA_ARRAY( UobsWt );                   // pred, prey, predL, preyL U matrix (mean wt_hat of prey in each pred); n = [nspp, nspp, maxL, maxL] #FIXME - Changed name in stomach2017.dat
  DATA_ARRAY( UobsAge );                  // pred, prey, predA, preyA U matrix (mean number of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_ARRAY( UobsWtAge );                // pred, prey, predA, preyA U matrix (mean wt_hat of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_MATRIX( Mn_LatAge );               // Mean length-at-age; n = [nspp, nages], ALSO: mean_laa in Kinzey

  // 2.3.8. Environmental data
  DATA_IVECTOR( Tyrs );                   // Years of hindcast data; n = [1, nTyrs] #FIXME - changed the name of this in retro_data2017_asssmnt.dat
  DATA_VECTOR( BTempC );                  // Vector of bottom temperature; n = [1,  nTyrs ]
  int nTyrs = Tyrs.size();                // Number of temperature years; n = [1] #FIXME - changed the name of this in retro_data2017_asssmnt.dat

  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( other_food );              // Biomass of other prey (kg); n = [1, nspp]
  DATA_VECTOR( Pvalue );                  // This scales the pvalue used, proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge; n = [1, nspp]
  DATA_IVECTOR( Ceq );                    // Ceq: which Comsumption equation to use; n = [1, nspp]; Currently all sp = 1
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
  DATA_MATRIX( M1_base );                 // Residual natural mortality; n = [nspp, nages]
  DATA_MATRIX( propF );                   // Proportion-at-age of females of population; n = [nspp, nages]
  DATA_MATRIX( pmature );                 // Proportion of mature females at age; [nspp, nages]

  // -- 2.4.5. F Profile data: NOTUSED


  // ------------------------------------------------------------------------- //
  // 2.7. Debugging with data inputs                                           //
  // ------------------------------------------------------------------------- //
  if (debug == 1) {
    /*
     // -- 2.8.2.1 Check to make sure the first year of survey data are not before start year
     for (sp = 0; sp < nspp; sp++) {
     if (yrs_tc_biom(sp, 0) < styr) {
     std::cerr << "First year of total catch biomass of species " << sp + 1 << " is before specified start year" << std::endl;
     return (0);
     }
     }

     // -- 2.8.2.2 Check to make sure the first year of survey data are not before start year
     for (sp = 0; sp < nspp; sp++) {
     if (yrs_srv_biom(sp, 0) < styr) {
     std::cerr << "First year of survey biomass of species " << sp + 1 << " is before specified start year" << std::endl;
     return (0);
     }
     }

     // -- 2.8.2.3. Check to make sure the years of survey data biomass and age are the same
     for (sp = 0; sp < nspp; sp++) {
     if (nyrs_srv_biom(sp) != nyrs_srv_age(sp)) {
     std::cerr << "Nyrs of survey biomass and age-comp of species " << sp + 1 << " do not match" << std::endl;
     return (0);
     }
     }

     if (yrs_eit(0) < styr) {
     std::cerr << "First year of EIT survey biomass is before specified start year" << std::endl;
     return (0);
     }
     */
  }


  // ------------------------------------------------------------------------- //
  // 3. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //

  PARAMETER(dummy);                               // Variable to test derived quantities given input parameters; n = [1]

  // -- 3.1. Recruitment parameters
  PARAMETER_VECTOR( ln_mn_rec );                  // Mean recruitment; n = [1, nspp]
  PARAMETER_VECTOR( ln_rec_sigma );               // Standard deviation of recruitment deviations; n = [1, nspp]
  PARAMETER_MATRIX( rec_dev );                    // Annual recruitment deviation; n = [nspp, nyrs]


  // -- 3.2. Abundance parameters
  PARAMETER_MATRIX( init_dev );                   // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this

  // -- 3.3. fishing mortality parameters
  PARAMETER_VECTOR( ln_mean_F );                  // Log mean fishing mortality; n = [1, nspp]
  PARAMETER_VECTOR( proj_F );                     // Fishing mortality for projections
  PARAMETER_MATRIX( F_dev );                      // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change

  // -- 3.4. Selectivity parameters
  PARAMETER_MATRIX( srv_sel_coff );               // Survey selectivity parameters; n = [n_srv, nselages]
  PARAMETER_MATRIX( srv_sel_slp );                // Survey selectivity paramaters for logistic; n = [2, n_srv]
  PARAMETER_MATRIX( srv_sel_inf );                // Survey selectivity paramaters for logistic; n = [2, n_srv]
  PARAMETER_VECTOR( log_srv_q );                  // Survey catchability; n = [n_srv]

  PARAMETER_MATRIX( fsh_sel_coff );               // Fishery age selectivity coef; n = [n_srv, nselages]
  PARAMETER_MATRIX( fsh_sel_slp );                // Fishery selectivity paramaters for logistic; n = [2, n_fsh]
  PARAMETER_MATRIX( fsh_sel_inf );                // Fishery selectivity paramaters for logistic; n = [2, n_fsh]

  // FIXME: Create maps
  // -- 3.5. Kinzery predation function parameters
  PARAMETER_MATRIX(logH_1);                       // Predation functional form; n = [nspp, nspp2];
  PARAMETER_VECTOR(logH_1a);                      // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix
  PARAMETER_VECTOR(logH_1b);                      // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix

  PARAMETER_MATRIX(logH_2);                       // Predation functional form; n = [nspp, nspp]
  PARAMETER_MATRIX(logH_3);                       // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
  PARAMETER_MATRIX(H_4);                          // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;

  // 3.6 Gamma selectivity parameters
  PARAMETER_VECTOR( log_gam_a );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  PARAMETER_VECTOR( log_gam_b );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10

  // 3.7. Preference
  PARAMETER_MATRIX( log_phi );                        // Species preference coefficient; n = [nspp, nspp]


  // ------------------------------------------------------------------------- //
  // 4. DERIVED QUANTITIES SECTION                                             //
  // ------------------------------------------------------------------------- //

  // 4.1. Derived indices
  int max_bin = imax( nlengths );                                                   // Integer of maximum number of length/age bins.
  int n_srv = srv_control.rows();
  int n_fsh = fsh_control.rows();

  // -- 4.2. Estimated population parameters
  vector<Type>  mn_rec = exp(ln_mn_rec);                                            // Mean recruitment; n = [1, nspp]
  array<Type>   AvgN(nspp, max_age, nyrs); AvgN.setZero();                          // Average numbers-at-age; n = [nspp, nages, nyrs]
  array<Type>   biomassByage(nspp, max_age, nyrs); biomassByage.setZero();          // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type>  biomass(nspp, nyrs); biomass.setZero();                             // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type>  biomassSSB(nspp, nyrs); biomassSSB.setZero();                       // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type>   biomassSSBByage(nspp, max_age, nyrs); biomassSSBByage.setZero();    // Spawning biomass at age (kg); n = [nspp, nages, nyrs]
  array<Type>   M(nspp, max_age, nyrs); M.setZero();                                // Total natural mortality at age; n = [nyrs, nages, nspp]
  matrix<Type>  M1(nspp, max_age); M1.setZero();                                    // Base natural mortality; n = [nspp, nages]
  array<Type>   M2(nspp, max_age, nyrs); M2.setZero();                              // Predation mortality at age; n = [nyrs, nages, nspp]
  array<Type>   NByage(nspp, max_age, nyrs); NByage.setZero();                      // Numbers at age; n = [nspp, nages, nyrs]
  matrix<Type>  R(nspp, nyrs); R.setZero();                                         // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type>   S(nspp, max_age, nyrs); S.setZero();                                // Survival at age; n = [nspp, nages, nyrs]
  array<Type>   Zed(nspp, max_age, nyrs); Zed.setZero();                            // Total mortality at age; n = [nspp, nages, nyrs]
  vector<Type>  r_sigma(nspp); r_sigma.setZero();                                   // Standard deviation of recruitment variation

  // -- 4.3. Fishery observations
  vector<Type>  avgsel_fsh(n_fsh); avgsel_fsh.setZero();                            // Average fishery selectivity; n = [1, nspp]
  array<Type>   fsh_sel(n_fsh, max_age, nyrs); fsh_sel.setZero();                   // Estimated fishery selectivity at age; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_sel2(n_fsh, max_age); fsh_sel2.setZero();                       // Temporary saved fishery selectivity at age for estimated bits; n = [nspp, nages, nyrs]

  array<Type>   catch_hat(nspp, max_age, nyrs); catch_hat.setZero();                // Estimated catch-at-age (n); n = [nspp, nages, nyrs]
  array<Type>   F(n_fsh, max_age, nyrs); F.setZero();                               // Estimated fishing mortality for each fishery; n = [n_fsh, nages, nyrs]
  array<Type>   F_tot(nspp, max_age, nyrs); F_tot.setZero();                        // Sum of annual estimated fishing mortalities for each species-at-age; n = [nspp, nages, nyrs]
  vector<Type>  fsh_bio_hat(fsh_biom_obs.rows()); fsh_bio_hat.setZero();            // Estimated fishery yield (kg); columns = Species, Index, Year, Month, Estimate
  vector<Type>  fsh_hat(fsh_comp_obs.rows()) ; fsh_hat.setZero() ;                  // Estimated fishery catch (n); n = [nspp, nyrs]
  matrix<Type>  fsh_age_hat = fsh_comp_obs; fsh_age_hat.setZero();                  // Estimated fishery catch at true age; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_age_obs_hat = fsh_comp_obs; fsh_age_obs_hat.setZero();          // Estimated fishery catch at observed age; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_comp_hat = fsh_comp_obs; fsh_comp_hat.setZero();                // Estimated fishery comp; n = [nspp, nages, nyrs]

  // -- 4.4. Survey components
  vector<Type>  avgsel_srv(n_srv); avgsel_srv.setZero();                          // Average survey selectivity; n = [1, nspp]
  array<Type>   srv_sel(n_srv, max_age, nyrs); srv_sel.setZero();                 // Estimated survey selectivity at age; n = [nspp, nages, nyrs]
  matrix<Type>  srv_sel_tmp(n_srv, max_age); srv_sel_tmp.setZero();               // Temporary saved survey selectivity at age for estimated bits; n = [nspp, nages, nyrs]

  Type avgsel_tmp = 0;                                                            // Temporary object for average selectivity across all ages
  vector<Type>  srv_bio_hat(srv_biom_obs.rows()); srv_bio_hat.setZero();          // Estimated survey biomass (kg); columns = Species, Index, Year, Month, Estimate
  vector<Type>  srv_hat( srv_comp_obs.rows() ); srv_hat.setZero();                // Estimated survey total abundance (n); n = [nspp, nyrs]
  matrix<Type>  srv_age_hat = srv_comp_obs; srv_age_hat.setZero();                // Estimated survey abundance at true age; columns = Comp_1, Comp_2, etc.
  matrix<Type>  srv_age_obs_hat = srv_comp_obs; srv_age_obs_hat.setZero();        // Estimated survey abundance at observed age; columns = Comp_1, Comp_2, etc.
  matrix<Type>  srv_comp_hat = srv_comp_obs; srv_comp_hat.setZero();              // Estimated survey comp; columns = Comp_1, Comp_2, etc.

  // -- 4.6. Ration components
  array<Type>   ConsumAge( nspp, max_age, nyrs ); ConsumAge.setZero();                        // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
  matrix<Type>  fT( nspp, nyrs ); fT.setZero();                                              // Pre-allocation of temperature function of consumption; n = [nspp, nTyrs]
  array<Type>   LbyAge( nspp, max_age, nyrs ); LbyAge.setZero();                              // Length by age from LW regression
  matrix<Type>  mnWt_obs( nspp, max_age ); mnWt_obs.setZero();                                // Mean observed weight at age (across years); n = [nspp, nages]
  array<Type>   ration2Age( nspp, max_age, nyrs ); ration2Age.setZero();                      // Annual ration at age (kg/yr); n = [nyrs, nages, nspp]
  vector<Type>  TempC( nyrs ); TempC.setZero();                                              // Bottom temperature; n = [1, nTyrs]

  // -- 4.7. Suitability components
  Type tmp_othersuit = 0  ;
  Type suit_tmp = 0;                                                                          //  Temporary storage variable
  array<Type>   avail_food(nspp, max_age, nyrs); avail_food.setZero();                        // Available food to predator; n = [nyrs, nages, nspp]
  array<Type>   B_eaten(nspp, max_age, nyrs); B_eaten.setZero();                              // Biomass of prey eaten via predation; n = [nyrs, nages, nspp]
  array<Type>   of_stomKir(nspp, max_age, nyrs); of_stomKir.setZero();                        // Other food stomach content; n = [nyrs, nages, nspp]
  array<Type>   stom_div_bio2(nspp, nspp, max_age, max_age, nyrs); stom_div_bio2.setZero();   // Stomach proportion over biomass; U/ (W * N) ; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   stomKir(nspp, nspp, max_age, max_age, nyrs); stomKir.setZero();               // Stomach proportion by numbers U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   stomKirWt(nspp, nspp, max_age, max_age, nyrs); stomKirWt.setZero();           // Stomach proportion by weight U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   diet_w_dat(nspp, nspp, max_bin, max_bin, nyrs); diet_w_dat.setZero();         // Observed stomach contents by weight of prey length j in predator length l
  array<Type>   diet_w_sum(nspp, nspp, max_bin, nyrs); diet_w_sum.setZero();                  // Observed stomach contentes by weight of prey in predator length j
  array<Type>   suit_main(nspp, nspp, max_age, max_age, nyrs); suit_main.setZero();           // Suitability/gamma selectivity of predator age u on prey age a; n = [nspp, nspp, nages, nages]
  matrix<Type>  suit_other(nspp, max_age); suit_other.setZero();                              // Suitability not accounted for by the included prey; n = [nspp, nages]
  array<Type>   suma_suit(nspp, max_age, nyrs); suma_suit.setZero();                          // Sum of suitabilities; n = [nyrs, nages, nspp]
  array<Type>   UobsWtAge_hat(nspp, nspp, max_age, max_age, nyrs); UobsWtAge_hat.setZero();   // Estimated stomach proportion by weight U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   mn_UobsWtAge_hat(nspp, nspp, max_age, max_age); mn_UobsWtAge_hat.setZero();   // Average estimated stomach proportion by weight U; n = [nspp, nspp, nages, nages]

  // -- 4.8. Kinzey selectivity
  vector<Type> gam_a = exp(log_gam_a); // Predator selectivity
  vector<Type> gam_b = exp(log_gam_b); // Predator selectivity

  // -- 4.9. Kinzey Functional response
  matrix<Type> H_1(nspp, nspp + 1); H_1 = exp(logH_1.array());                       // FIXME: make matrix
  vector<Type> H_1a(nspp); H_1a = exp(logH_1a);                                      // FIXME: make matrix
  vector<Type> H_1b(nspp); H_1b = exp(logH_1b);                                      // FIXME: make matrix
  matrix<Type> H_2(nspp, nspp); H_2 = exp(logH_2.array());
  matrix<Type> H_3(nspp, nspp); H_3 = exp(logH_3.array());

  array<Type>  N_pred_yrs(nspp, max_age, nyrs); N_pred_yrs.setZero();                // Effective numbers of predators for each age of prey
  array<Type>  N_prey_yrs(nspp, max_age, nyrs); N_prey_yrs.setZero();                // Effective numbers of prey for each age of predator
  matrix<Type> N_pred_eq(nspp, max_age); N_pred_eq.setZero();                        // Effective numbers of predators for each age of prey (styr_pred)
  matrix<Type> N_prey_eq(nspp, max_age); N_prey_eq.setZero();                        // Effective numbers of prey for each age of predator

  array<Type>  pred_resp(nspp, nspp + 1, max_age, max_age, nyrs); pred_resp.setZero();     // Predator functional response
  array<Type>  Pred_r(nspp, max_age, nyrs); Pred_r.setZero();                        // save Pred_ratio values
  array<Type>  Prey_r(nspp, max_age, nyrs); Prey_r.setZero();                        // save Prey_ratio values

  array<Type>  Vmort_ua(nspp, nspp, max_age, max_age, nyrs); Vmort_ua.setZero();     // Predation mortality on prey age a by single predator age u
  array<Type>  Pmort_ua(nspp, nspp, max_age, nyrs); Pmort_ua.setZero();              // Predation mortality on prey age a by all predators age u

  array<Type>  eaten_la(nspp, nspp, max_bin, max_age, nyrs); eaten_la.setZero();     // Number of prey of age a eaten by predator length l
  array<Type>  eaten_ua(nspp, nspp, max_age, max_age, nyrs); eaten_ua.setZero();     // Number of prey of age a eaten by predator age u

  array<Type>  Q_mass_l(nspp, nspp + 1, max_bin, nyrs); Q_mass_l.setZero();          // Mass of each prey sp consumed by predator at length // FIXME: make into 4D array
  array<Type>  Q_mass_u(nspp, nspp + 1, max_bin, nyrs); Q_mass_u.setZero();          // Mass of each prey sp consumed by predator at age // FIXME: make into 4D array
  matrix<Type> Q_other_u(nspp, max_age); Q_other_u.setZero();                        // Mass of other prey consumed by predator at age
  array<Type>  Q_hat(nspp, nspp + 1, max_bin, nyrs); Q_hat.setZero();                // Fraction for each prey type of total mass eaten by predator length
  array<Type>  T_hat(nspp, nspp, max_bin, max_bin); T_hat.setZero();                 // Fraction of prey of length m in predator of length l

  array<Type> omega_hat(nspp, max_age, nyrs); omega_hat.setZero();                   // Daily ration by predator age each year
  matrix<Type> omega_hat_ave(nspp, max_age); omega_hat_ave.setZero();                // Daily ration by predator age averaged over years


  // ------------------------------------------------------------------------- //
  // 5. INITIAL CALCULATIONS                                                   //
  // ------------------------------------------------------------------------- //

  // 5.5. Maturity and sex ratio
  M1 = M1_base.array() + Type(0.0001);
  for ( sp = 0; sp < nspp ; sp++) {
    for ( age = 0 ; age < nages(sp); age++ ) {
      pmature( sp, age ) = pmature( sp, age ) * propF(sp, age);
    }
  }

  // 5.6. Calculate temperature to use
  TempC = BTempC.sum() / nTyrs; // Fill with average bottom temperature

  int yr_ind = 0;
  for (yr = 0; yr < nTyrs; yr++) {
    yr_ind = Tyrs( yr ) - styr;
    if ((yr_ind >= 0) & (yr_ind < nyrs)) {
      TempC(yr_ind) = BTempC( yr );
    }
  }


  // 5.7. Calculate length-at-age
  for (sp = 0; sp < nspp; sp++) {
    for (age = 0; age < nages(sp); age++) {
      for (yr = 0; yr < nyrs; yr++) {
        // Hindcast
        if(yr < nyrs_hind){
          LbyAge( sp, age, yr) = ( pow( ( 1 / aLW(0, sp) ), (1 / aLW(1, sp) ) ) )  * pow( ( wt(yr, age, pop_wt_index(sp)) * 1000), (1 / aLW(1, sp))); // W = a L ^ b is the same as (W/a)^(1/b)
        }
        if(yr >= nyrs_hind){
          LbyAge( sp, age, yr) = ( pow( ( 1 / aLW(0, sp) ), (1 / aLW(1, sp) ) ) )  * pow( ( wt((nyrs_hind - 1), age, pop_wt_index(sp)) * 1000), (1 / aLW(1, sp))); // W = a L ^ b is the same as (W/a)^(1/b)
        }
        }
    }
  }


  // 5.8. Parameter transform
  r_sigma = exp(ln_rec_sigma); // Convert log sd to natural scale

  // 5.1. Reorganize survey control bits
  vector<Type> srv_q(n_srv); srv_q.setZero();                                   // Vector to save q on natural scale
  vector<int> srv_sel_type(n_srv); srv_sel_type.setZero();                      // Vector to save survey selectivity type
  vector<int> srv_nselages(n_srv); srv_nselages.setZero();                      // Vector to save number of ages to estimate non-parametric selectivity (1 = age, 2 = length)
  vector<int> srv_spp(n_srv); srv_spp.setZero();                                // Vector to save survey species
  vector<int> srv_units(n_srv); srv_units.setZero();                            // Vector to save survey units (1 = weight, 2 = numbers)
  vector<int> srv_wt_index(n_srv); srv_wt_index.setZero();                      // Vector to save 3rd dim of wt to use for weight-at-age
  vector<int> srv_alk_index(n_srv); srv_alk_index.setZero();                    // Vector to save 3rd dim of age_trans_matrix to use for ALK

  for (srv_ind = 0; srv_ind < n_srv; srv_ind++){
    srv = srv_control(srv_ind, 1) - 1;                    // Temporary survey index
    srv_q(srv) = exp(log_srv_q(srv_ind));                 // Exponentiate
    srv_spp(srv) = srv_control(srv_ind, 2) - 1;           // Species
    srv_sel_type(srv) = srv_control(srv_ind, 3);          // Selectivity type
    srv_nselages(srv) = srv_control(srv_ind, 4);          // Non-parametric selectivity ages
    srv_units(srv) = srv_control(srv_ind, 5);             // Survey units
    srv_wt_index(srv) = srv_control(srv_ind, 6) - 1;      // Dim3 of wt
    srv_alk_index(srv) = srv_control(srv_ind, 7) - 1;     // Dim3 of ALK
  }


  // 5.2. Reorganize fishery control bits
  vector<int> fsh_sel_type(n_fsh); fsh_sel_type.setZero();                      // Vector to save fishery selectivity type
  vector<int> fsh_nselages(n_fsh); fsh_nselages.setZero();                      // Vector to save number of ages to estimate non-parametric selectivity (1 = age, 2 = length)
  vector<int> fsh_spp(n_fsh); fsh_spp.setZero();                                // Vector to save specives of fishery
  vector<int> fsh_units(n_fsh); fsh_units.setZero();                            // Vector to save survey units (1 = weight, 2 = numbers)
  vector<int> fsh_wt_index(n_fsh); fsh_wt_index.setZero();                      // Vector to save 3rd dim of wt to use for weight-at-age
  vector<int> fsh_alk_index(n_fsh); fsh_alk_index.setZero();                    // Vector to save 3rd dim of age_trans_matrix to use for ALK

  for (fsh_ind = 0; fsh_ind < n_fsh; fsh_ind++){
    fsh = fsh_control(fsh_ind, 1) - 1;                    // Temporary survey index
    fsh_spp(fsh) = fsh_control(fsh_ind, 2) - 1;           // Species
    fsh_sel_type(fsh) = fsh_control(fsh_ind, 3);          // Selectivity type
    fsh_nselages(fsh) = fsh_control(fsh_ind, 4);          // Non-parametric selectivity ages
    fsh_units(fsh) = fsh_control(fsh_ind, 5);             // Fishery units
    fsh_wt_index(fsh) = fsh_control(fsh_ind, 6) - 1;      // Dim3 of wt
    fsh_alk_index(fsh) =  fsh_control(fsh_ind, 7) - 1;    // Dim3 of ALK
  }


  // 5.3. Normalize fishery composition
  for(comp_ind = 0; comp_ind < fsh_comp_obs.rows(); comp_ind++){
    fsh = fsh_comp_ctl(comp_ind, 0)-1;
    // Get sum of composition
    Type sum_temp = 0;
    for(int col = 0; col < fsh_comp_obs.cols(); col ++){
      if(!isNA(fsh_comp_obs(comp_ind, col))){
        sum_temp += fsh_comp_obs(comp_ind, col);
      }
    }
    // Divide each row by sum
    fsh_comp_obs.row(comp_ind) = fsh_comp_obs.row(comp_ind) / sum_temp;
  }


  // 5.4. Normalize survey composition
  for(comp_ind = 0; comp_ind < srv_comp_obs.rows(); comp_ind++){
    // Get sum of composition
    Type sum_temp = 0;
    for(int col = 0; col < srv_comp_obs.cols(); col ++){
      if(!isNA(srv_comp_obs(comp_ind, col))){
        sum_temp += srv_comp_obs(comp_ind, col);
      }
    }
    // Divide each row by sum
    srv_comp_obs.row(comp_ind) = srv_comp_obs.row(comp_ind) / sum_temp;
  }

  /*
   // 5.5. Normalize age-transition matrix
   for(sp = 0; sp < nspp; sp ++){
   Type sum_temp = 0;
   for (age = 0; age < nages(sp); age++) {
   for (ln = 0; ln < nlengths(sp); ln++) {
   sum_temp +=  age_trans_matrix(age, ln, sp);
   }
   for (ln = 0; ln < nlengths(sp); ln++) {
   age_trans_matrix(age, ln, sp) = age_trans_matrix(age, ln, sp) / sum_temp;
   }
   }
   }
   */

  // Dimensions of diet matrix
  vector<int> a1_dim = UobsAge.dim; // dimension of a1
  vector<int> a2_dim = UobsWtAge.dim; // dimension of a1

  // Good above here
  // ------------------------------------------------------------------------- //
  // 6. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0
  // Start iterations
  for (int iter = 0; iter < niter; iter++) {


    // 6.0. EMPIRICAL FISHERY SELECTIVITY
    avgsel_fsh.setZero();
    fsh_sel2.setZero();
    fsh_sel.setZero();
    for (int sel_ind = 0; sel_ind < fsh_emp_sel_obs.rows(); sel_ind++){

      fsh = fsh_emp_sel_ctl(sel_ind, 0) - 1;            // Temporary survey index
      sp = fsh_emp_sel_ctl(sel_ind, 1) - 1;             // Temporary index of species
      fsh_yr = fsh_emp_sel_ctl(sel_ind, 2) - styr;      // Temporary index for years of data

      if(fsh_yr < nyrs){
        for (age = 0; age < fsh_emp_sel_obs.cols(); age++) {
          if(!isNA(fsh_emp_sel_obs(sel_ind, age))){
            fsh_sel(fsh, age, fsh_yr) = fsh_emp_sel_obs(sel_ind, age);
          }
        }
      }
      // FIXME - set all fishing selectivities after nyrs_hind to the terminal year
    }

    // 6.1. ESTIMATE FISHERY SELECTIVITY
    for (fsh = 0; fsh < n_fsh; fsh++) {

      sp = fsh_spp(fsh);             // Temporary index of species

      // 6.1.1. Logisitic selectivity
      if (fsh_sel_type(fsh) == 1) {
        for (age = 0; age < nages(sp); age++){
          fsh_sel2(fsh, age) = 1 / (1 + exp( -fsh_sel_slp(0, fsh) * ((age + 1) - fsh_sel_inf(0, fsh))));
        }
      }

      // 6.1.2. Non-parametric selectivity fit to age ranges. NOTE: This can likely be improved
      if (fsh_sel_type(fsh) == 2) {
        for (age = 0; age < fsh_nselages(fsh); age++) {
          fsh_sel2(fsh, age) = fsh_sel_coff(fsh, age);
          avgsel_fsh(fsh) +=  exp(fsh_sel_coff(fsh, age));
        }
        //  Average selectivity up to nselages
        avgsel_fsh(fsh) = log(avgsel_fsh(fsh) / fsh_nselages(fsh));

        // Plus group selectivity
        for (age = fsh_nselages(fsh); age < nages(sp); age++) {
          fsh_sel2(fsh, age) = fsh_sel2(fsh, fsh_nselages(fsh) - 1);
        }

        // Average selectivity across all ages
        avgsel_tmp = 0; // Temporary object for average selectivity across all ages
        for (age = 0; age < nages(sp); age++) {
          avgsel_tmp += exp(fsh_sel2(fsh, age));
        }
        avgsel_tmp = log(avgsel_tmp / nages(sp));

        // Standardize selectivity
        for (age = 0; age < nages(sp); age++) {
          fsh_sel2(fsh, age) -= avgsel_tmp;
          fsh_sel2(fsh, age) = exp(fsh_sel2(fsh, age));
        }
      }


      // 6.1.3. Double logistic (Dorn and Methot 1990)
      if (fsh_sel_type(fsh) == 3) {
        for (age = 0; age <  nages(sp); age++){
          fsh_sel2(fsh, age) = (1 / (1 + exp( -fsh_sel_slp(0, fsh) * ((age + 1) - fsh_sel_inf(0, fsh))))) * // Upper slop
            (1 - (1 / (1 + exp( -fsh_sel_slp(1, fsh) * ((age + 1) - fsh_sel_inf(1, fsh))))));  // Downward slope;
        }
      }

      // 6.1.4. Normalize
      if (fsh_sel_type(fsh) > 0) {
        vector<Type> sel_sp = fsh_sel2.row(fsh); // Can't max a matrix....
        fsh_sel2.row(fsh) /= max(sel_sp); // Standardize so max sel = 1 for each species
      }

      for (age = 0; age < nages(sp); age++){
        for (yr = 0; yr < nyrs; yr++) {
          fsh_sel(fsh, age, yr) = fsh_sel2(fsh, age);
        }
      }
    }


    // 6.1. ESTIMATE FISHING MORTALITY
    F.setZero();
    F_tot.setZero();
    for (fsh_ind = 0; fsh_ind < n_fsh; fsh_ind++) {

      sp = fsh_spp(fsh_ind);  // Temporary index of fishery survey

      for (yr = 0; yr < nyrs; yr++) {
        for (age = 0; age < nages(sp); age++) {
          // Hindcast
          if( yr < nyrs_hind){
F(fsh_ind, age, yr) = fsh_sel(fsh_ind, age, yr) * exp(ln_mean_F(fsh_ind) + F_dev(fsh_ind, yr));
          }
if( yr >= nyrs_hind){
F(fsh_ind, age, yr) = fsh_sel(fsh_ind, age, yr) * proj_F(fsh_ind);
          }
          
          // F(fsh_ind, age, yr) = fsh_sel(fsh_ind, age, yr) * exp(ln_mean_F(fsh_ind) + F_dev(fsh_ind, yr));
          F_tot(sp, age, yr) += F(fsh_ind, age, yr);

          // FIXME - make f zero for projections
        }
      }
    }

    // Good above here


    // 6.2. Estimate total mortality at age NOTE: May need to go above population dynamics
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          M(sp, age, yr) = M1(sp, age) + M2(sp, age, yr);
          Zed(sp, age, yr) = M1(sp, age) + F_tot(sp, age, yr) + M2(sp, age, yr);
          S(sp, age, yr) = exp(-Zed(sp, age, yr));
        }
      }
    }


    // 6.3. ESTIMATE RECRUITMENT T1.1
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs; yr++) {

        // Hindcast
        if(yr < nyrs_hind){
                  R(sp, yr) = exp(ln_mn_rec(sp) + rec_dev(sp, yr));
        } 

        // Project under mean recruitment
        if(yr >= nyrs_hind){
          R(sp, yr) = exp(ln_mn_rec(sp));
        }

        NByage(sp, 0, yr) = R(sp, yr);
      }
    }

    // Good above here


    // 6.4. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        if ((age > 0) & (age < nages(sp) - 1)) {
          Type mort_sum = M1.row(sp).segment(0, age).sum(); // Sum M1 until age - 1
          NByage(sp, age, 0) = exp(ln_mn_rec(sp) - mort_sum + init_dev(sp, age - 1));
        }
        // -- 6.2.2. Where yr = 1 and age > Ai.
        if (age == (nages(sp) - 1)) {
          Type mort_sum = M1.row(sp).segment(0, age).sum(); // Sum M1 until age - 1
          NByage(sp, age, 0) = exp(ln_mn_rec(sp) - mort_sum + init_dev(sp, age - 1)) / (1 - exp(-M1(sp, nages(sp) - 1))); // NOTE: This solves for the geometric series
        }
      }
    }


    // 6.5. ESTIMATE NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
    biomass.setZero();
    biomassSSB.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 1; yr < nyrs; yr++) {
          // -- 6.3.1.  Where 1 <= age < Ai
          if (age < (nages(sp) - 1)) {
            NByage(sp, age + 1, yr) = NByage(sp, age, yr - 1) * S(sp, age, yr - 1);
          }

          // -- 6.3.2. Plus group where age > Ai. NOTE: This is not the same as T1.3 because sp used age = A_i rather than age > A_i.
          if (age == (nages(sp) - 1)) {
            NByage(sp, nages(sp) - 1, yr) = NByage(sp, nages(sp) - 2, yr - 1) * S(sp, nages(sp) - 2, yr - 1) + NByage(sp, nages(sp) - 1, yr - 1) * S(sp, nages(sp) - 1, yr - 1);
          }
        }

        // -- 6.3.3. Estimate Biomass and SSB
        for (yr = 0; yr < nyrs; yr++) {
          // Hindcast
          if(yr < nyrs_hind){
                      biomassByage(sp, age, yr) = NByage(sp, age, yr) * wt(yr, age, pop_wt_index(sp)); // 6.5.
          }
          // Projection - Use last year of wt
          if(yr >= nyrs_hind){
                      biomassByage(sp, age, yr) = NByage(sp, age, yr) * wt((nyrs_hind - 1), age, pop_wt_index(sp)); // 6.5.
          }

          biomassSSBByage(sp, age, yr) = biomassByage(sp, age, yr) * pmature(sp, age); // 6.6.

          biomass(sp, yr) += biomassByage(sp, age, yr);
          biomassSSB(sp, yr) += biomassSSBByage(sp, age, yr);
        }
      }
    }


    // 6.6. ESTIMATE AVERAGE NUMBERS AT AGE
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          if (avgnMode == 0) {
            AvgN(sp, age, yr) = NByage(sp, age, yr) * (1 - S(sp, age, yr)) / Zed(sp, age, yr); // MSVPA approach
          }else if (avgnMode == 1) {
            AvgN(sp, age, yr) = NByage(sp, age, yr) * exp(- Zed(sp, age, yr) / 2); // Kinzey and Punt (2009) approximation
          }else if (avgnMode == 2) {
            AvgN(sp, age, yr) = NByage(sp, age, yr); // Van Kirk et al (2010) approximation
          } else{
            error("Invalid 'avgnMode'");
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
          fT(sp, yr) = exp(Qc(sp) * TempC(yr));
        }

        // Temperature dependence for warm-water-species from Kitchell et al. 1977
        if ( Ceq(sp) == 2) {
          Yc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp) + 2);
          Zc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp));
          Vc = (Tcm(sp) - TempC(yr)) / (Tcm(sp) - Tco(sp));
          Xc = pow(Zc, 2) * pow((1 + pow((1 + 40 / Yc), 0.5)), 2) / 400;
          fT(sp, yr) = pow(Vc, Xc) * exp(Xc * (1 - Vc));
        }

        // Temperature dependence for cool and cold-water species from Thornton and Lessem 1979
        if (Ceq(sp) == 3) {
          G2 = (1 / (Tcl(sp) - Tcm(sp))) * log((0.98 * (1 - CK4(sp))) / (CK4(sp) * 0.02));
          L2 = exp(G2 * (Tcl( sp ) - TempC( yr )));
          Kb = (CK4(sp) * L2) / (1 + CK4(sp) * (L2 - 1));
          G1 = (1 / (Tco(sp) - Qc(sp))) * log((0.98 * (1 - CK1(sp))) / (CK1(sp) * 0.02));
          L1 = exp(G1 * (TempC(yr) - Qc(sp)));
          Ka = (CK1(sp) * L1) / (1 + CK1(sp) * (L1 - 1));
          fT(sp, yr) = Ka * Kb;
        }
    }
  }


    // 7.3. Calculate historic ration
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          // p = proportion of maximum consumption
          // f(T) = temperature dependence function
          // CA = intercept of allometric mass function
          // CB = slope of allometric mass function
          // fday = number of forageing days per year

// Hindcast
if(yr < nyrs_hind){
          ConsumAge(sp, age, yr) = CA(sp) * pow(wt(yr, age, pop_wt_index(sp)) * Type(1000), CB( sp )) // C_max = CA * W ^ CB; where C_max is grams consumed per grams of predator per day
          * fT(sp, yr) * fday( sp ) * wt(yr, age, pop_wt_index(sp)) * 1000; //  C_max * f(T) * wt * fday g/pred.yr
          ConsumAge(sp, age, yr) = ConsumAge(sp, age, yr) * Pvalue(sp) * Pyrs(yr, age, sp); //
}
// Projection
if(yr >= nyrs_hind){
          ConsumAge(sp, age, yr) = CA(sp) * pow(wt((nyrs_hind - 1), age, pop_wt_index(sp)) * Type(1000), CB( sp )) // C_max = CA * W ^ CB; where C_max is grams consumed per grams of predator per day
          * fT(sp, yr) * fday( sp ) * wt((nyrs_hind - 1), age, pop_wt_index(sp)) * 1000; //  C_max * f(T) * wt * fday g/pred.yr
          ConsumAge(sp, age, yr) = ConsumAge(sp, age, yr) * Pvalue(sp) * Pyrs((nyrs_hind-1), age, sp); //
}

          ration2Age(sp, age, yr) = ConsumAge(sp, age, yr) / 1000;      // Annual ration kg/yr //aLW(predd)*pow(lengths(predd,age),bLW(predd));//mnwt_bin(predd,age);
        }
      }
    }


    // 7.4. Calculate stomach content
    diet_w_sum.setZero();
    for (yr = 0; yr < nyrs; yr++) {                             // Year loop
      for (rsp = 0; rsp < nspp; rsp++) {                        // Predator species loop
        for (ksp = 0; ksp < nspp; ksp++) {                      // Prey species loop
          for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
            for (k_age = 0; k_age < nages(ksp); k_age++) {      // Prey age loop
              // Change to 5D array if 4D
              if(a1_dim.size() == 4){
                stomKir(rsp, ksp, r_age, k_age, yr) = UobsAge(rsp , ksp , r_age, k_age);
              }
              if(a1_dim.size() == 5){
                // Hindcast
                if(yr < nyrs_hind){
                                 stomKir(rsp, ksp, r_age, k_age, yr) = UobsAge(rsp , ksp , r_age, k_age, yr); 
                }

                // Projection - use last year
                if(yr >= nyrs_hind){
                                 stomKir(rsp, ksp, r_age, k_age, yr) = UobsAge(rsp , ksp , r_age, k_age, nyrs_hind - 1); 
                }

              }

              if(a2_dim.size() == 4){
                stomKirWt(rsp, ksp, r_age, k_age, yr) = UobsWtAge(rsp, ksp, r_age, k_age);
              }
              if(a2_dim.size() == 5){
                // Hindcast
                if(yr < nyrs_hind){
                stomKirWt(rsp, ksp, r_age, k_age, yr) = UobsWtAge(rsp, ksp, r_age, k_age, yr);
              }

// Projection - use last year
                if(yr >= nyrs_hind){
stomKirWt(rsp, ksp, r_age, k_age, yr) = UobsWtAge(rsp, ksp, r_age, k_age, nyrs_hind - 1);
              }
              }
            }
          }
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


    // 7.5. Calculate other food stomach content
    of_stomKir.setZero();
    for (yr = 0; yr < nyrs; yr++) {                           // Year loop
      for (rsp = 0; rsp < nspp; rsp++) {                      // Predator species loop
        for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
          of_stomKir(rsp, r_age, yr) = Type( 1 );             // Initialize other suitability
          for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
            for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
              of_stomKir(rsp, r_age, yr) -= stomKir(rsp, ksp, r_age, k_age, yr);
            }
          }
        }
      }
    }
    of_stomKir /= other_food(0);

    // Good above here

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
              for (r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {  // Prey age loop
                  suit_tmp = stomKir(rsp, ksp, r_age, k_age, yr) / (AvgN(ksp, k_age, yr));

if(yr < nyrs_hind){
  int yr_ind = yr;
}
if(yr >= nyrs_hind){
  int yr_ind = nyrs_hind - 1;
}

                  if (wt(yr_ind, k_age, pop_wt_index(ksp)) != 0) {
                    stom_div_bio2(rsp, ksp, r_age, k_age, yr) = suit_tmp / wt(yr_ind, k_age, pop_wt_index(ksp));
                    suma_suit(rsp, r_age, yr ) += stom_div_bio2(rsp, ksp, r_age, k_age, yr); // Calculate sum of stom_div_bio2 across prey and  prey age for each predator, predator age, and year
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
          for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
            suit_other(rsp, r_age) = Type( 1 );                 // Initialize other suitability
            for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
              for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, 0) += stom_div_bio2(rsp, ksp, r_age, k_age, yr) / ( suma_suit(rsp, r_age, yr ) + of_stomKir(rsp, r_age, yr) );
                }                                 // End year loop
                suit_main(rsp, ksp, r_age, k_age, 0) /= nyrs;
                suit_other(rsp, r_age) -= suit_main(rsp, ksp, r_age, k_age, 0); // Subtract observed suitability from entire suitability (ksp.e. 1)

                // Fill in years
                for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, yr) = suit_main(rsp, ksp, r_age, k_age, 0);
                }
              }
            }
          }
        }
      } // End Holsman/MSVPA suitability




      // 8.1.2. Length-based GAMMA suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 1){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
              ncnt = 0;
              gsum = 1.0e-10;                                           // Initialize
              for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                // if prey are smaller than predator:
                if (Mn_LatAge(rsp, r_age) > Mn_LatAge(ksp, k_age)) {
                  x_l_ratio = log(Mn_LatAge(rsp, r_age) / Mn_LatAge(ksp, k_age));
                  suit_main(rsp , ksp, r_age, k_age, 0) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                    (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                  ncnt += 1;
                  gsum += exp( suit_main(rsp , ksp, r_age, k_age, 0) );
                }
                else
                  suit_main(rsp , ksp, r_age, k_age, 0) = 0;
              }
              for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                // if prey are smaller than predator:
                if (Mn_LatAge(rsp, r_age) > Mn_LatAge(ksp, k_age)) {
                  suit_main(rsp , ksp, r_age, k_age, 0) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_age, k_age, 0) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                }

                // Fill in years
                for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, yr) = suit_main(rsp, ksp, r_age, k_age, 0);
                }
              }
            }
          }
        }
      } // End GAMMA selectivity


      // 8.1.2. Time-varying length-based GAMMA suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 2){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
              for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                ncnt = 0;
                gsum = 1.0e-10;                                           // Initialize
                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if ( LbyAge( rsp, r_age, yr)  > LbyAge( ksp, k_age, yr)) {
                    x_l_ratio = log(LbyAge( rsp, r_age, yr) / LbyAge( ksp, k_age, yr) ); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                      (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                    ncnt += 1;
                    gsum += exp( suit_main(rsp , ksp, r_age, k_age, yr) );
                  }
                  else
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                }
                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if ( LbyAge( rsp, r_age, yr)  > LbyAge( ksp, k_age, yr)) {
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_age, k_age, yr) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                  }
                }
              }
            }
          }
        }
      } // End GAMMA selectivity


      // 8.1.3. Time-varying weight-based GAMMA suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 3){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
              for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                ncnt = 0;
                gsum = 1.0e-10;                                           // Initialize

                if(yr < nyrs_hind){
  int yr_ind = yr;
}
if(yr >= nyrs_hind){
  int yr_ind = nyrs_hind - 1;
}

                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if (wt(yr_ind, r_age, pop_wt_index(rsp)) > wt(yr_ind, k_age, pop_wt_index(ksp))) {
                    x_l_ratio = log(wt(yr_ind, r_age, pop_wt_index(rsp)) / wt(yr_ind, k_age, pop_wt_index(ksp))); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                      (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                    ncnt += 1;
                    gsum += exp( suit_main(rsp , ksp, r_age, k_age, yr) );
                  }
                  else
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                }
                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if (wt(yr_ind, r_age, pop_wt_index(rsp)) > wt(yr_ind, k_age, pop_wt_index(ksp))) {
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_age, k_age, yr) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                  }
                }
              }
            }
          }
        }
      } // End GAMMA selectivity


      //-------------------------------------------------------------------------------------------------------------
      // 8.1.4. Length-based lognormal suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 4){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
          for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
            gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
            for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
              for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                // if prey are smaller than predator:
                if (Mn_LatAge(rsp, r_age) > Mn_LatAge(ksp, k_age)) {
                  x_l_ratio = log(Mn_LatAge(rsp, r_age) / Mn_LatAge(ksp, k_age)); // Log ratio of lengths
                  suit_main(rsp , ksp, r_age, k_age, 0) = exp(log_phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                  gsum += suit_main(rsp , ksp, r_age, k_age, 0);
                }
                else
                  suit_main(rsp , ksp, r_age, k_age, 0) = 0;
              }
            }
            for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
              for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                suit_main(rsp , ksp, r_age, k_age, 0) /= gsum;                // Scale, so it sums to 1.

                // Fill in years
                for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, yr) = suit_main(rsp, ksp, r_age, k_age, 0);
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
          for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
            for(yr = 0; yr < nyrs; yr++){
              gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  // if prey are smaller than predator:
                  if ( LbyAge( rsp, r_age, yr)  > LbyAge( ksp, k_age, yr)) {
                    x_l_ratio = log(LbyAge( rsp, r_age, yr) / LbyAge( ksp, k_age, yr) ); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = exp(log_phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                    gsum += suit_main(rsp , ksp, r_age, k_age, yr);
                  }
                  else{
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                  }
                }
              }
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  suit_main(rsp , ksp, r_age, k_age, yr) /= gsum;                // Scale, so it sums to 1.
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
          for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
            for(yr = 0; yr < nyrs; yr++){
              gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
              
                              if(yr < nyrs_hind){
  int yr_ind = yr;
}
if(yr >= nyrs_hind){
  int yr_ind = nyrs_hind - 1;
}

              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  // if prey are smaller than predator:
                  if (wt(yr_ind, r_age, pop_wt_index(rsp)) > wt(yr_ind, k_age, pop_wt_index(ksp))) {
                    x_l_ratio = log(wt(yr_ind, r_age, pop_wt_index(rsp)) / wt(yr_ind, k_age, pop_wt_index(ksp))); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = exp(log_phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                    gsum += suit_main(rsp , ksp, r_age, k_age, yr);
                  }
                  else{
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                  }
                }
              }
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  suit_main(rsp , ksp, r_age, k_age, yr) /= gsum;                // Scale, so it sums to 1.
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
      if (msmMode == 1) {
        // 8.1.3. Calculate available food
        avail_food.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                        // Predator species loop
          for (r_age = 0; r_age < nages(rsp); r_age++) {          // Predator age loop
            for (yr = 0; yr < nyrs; yr++) {                       // Year loop
              tmp_othersuit = 0.;

                              if(yr < nyrs_hind){
  int yr_ind = yr;
}
if(yr >= nyrs_hind){
  int yr_ind = nyrs_hind - 1;
}

              for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                  avail_food(rsp, r_age, yr) += suit_main(rsp, ksp, r_age, k_age, yr) * AvgN(ksp, k_age, yr) * wt(yr_ind, k_age, pop_wt_index(ksp)); // FIXME - include overlap indices: FIXME - mn_wt_stom?
                  tmp_othersuit += suit_main(rsp, ksp, r_age, k_age, yr); // FIXME - include overlap indices
                }
              }
              if(suitMode == 0){ // Holsman
                avail_food(rsp, r_age, yr) += other_food(rsp) * (Type(1) - (tmp_othersuit));
              }
              if(suitMode > 0){ // Log-normal and gamma
                avail_food(rsp, r_age, yr) += other_food(rsp) * Type(1);
              }
            }
          }
        }


        // 8.1.3. Calculate predation mortality
        M2.setZero();
        B_eaten.setZero();
        mn_UobsWtAge_hat.setZero();
        for (ksp = 0; ksp < nspp; ksp++) {                        // Prey species loop
          for (k_age = 0; k_age < nages(ksp); k_age++) {          // Prey age loop
            for (yr = 0; yr < nyrs; yr++) {                       // Year loop

                              if(yr < nyrs_hind){
  int yr_ind = yr;
}
if(yr >= nyrs_hind){
  int yr_ind = nyrs_hind - 1;
}

              for (rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop
                for (r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                  M2(ksp, k_age, yr) += (AvgN(rsp, r_age, yr) * ration2Age(rsp, r_age, yr) * suit_main(rsp , ksp , r_age, k_age, yr)) / avail_food(rsp, r_age, yr); // #FIXME - include indices of overlap
                  B_eaten(ksp, k_age, yr) += AvgN(rsp, r_age, yr) * ration2Age(rsp, r_age, yr) * suit_main(rsp , ksp , r_age, k_age, yr);

                  // Estimated stomach proportion
                  UobsWtAge_hat(rsp, ksp, r_age, k_age, yr) = (AvgN(ksp, k_age, yr) * suit_main(rsp , ksp , r_age, k_age, yr) * wt(yr_ind, k_age, pop_wt_index(ksp))) / avail_food(rsp, r_age, yr);
                  mn_UobsWtAge_hat(rsp, ksp, r_age, k_age) += UobsWtAge_hat(rsp, ksp, r_age, k_age, yr) / ( endyr - styr + 1) ;
                }
              }
            }
          }
        }
      } // End 8..1. Holsman predation mortality


      // 8.2. KINZEY PREDATION EQUATIONS
      if (msmMode > 1) {
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
            for (r_age = 0; r_age < nages(rsp); r_age++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                N_pred_eq(rsp, r_age) += NByage(rsp, r_age, 0) * suit_main(rsp , ksp, r_age, k_age, yr); // Denominator of Eq. 17 Kinzey and Punt (2009) 1st year
              }
            }
            for (k_age = 0; k_age < nages(ksp); k_age++) {
              for (r_age = 0; r_age < nages(rsp); r_age++) {
                N_prey_eq(ksp, k_age) += NByage(ksp, k_age, 0) * suit_main(rsp , ksp, r_age, k_age, yr); // Denominator of Eq. 16 Kinzey and Punt (2009) 1st year
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
              for (r_age = 0; r_age < nages(rsp); r_age++) {
                for (k_age = 0; k_age < nages(ksp); k_age++) {
                  N_pred_yrs(rsp, r_age, yr) += NByage(rsp, r_age, yr) * suit_main(rsp , ksp, r_age, k_age, yr); // Numerator of Eq. 17 Kinzey and Punt (2009) 1st year // FIXME: Use averageN?
                }
              }
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                for (r_age = 0; r_age < nages(rsp); r_age++) {
                  N_prey_yrs(ksp, k_age, yr) += NByage(ksp, k_age, yr) * suit_main(rsp , ksp, r_age, k_age, yr); // Numerator of Eq. 16 Kinzey and Punt (2009) 1st year
                }
              }
            }
          }
        }


        // 8.2.6. Calculate predator functional response (Table 1 Kinzey and Punt (2009))
        for (rsp = 0; rsp < nspp; rsp++) {                          // Predator loop
          for (ksp = 0; ksp < (nspp + 1); ksp++) {                  // Prey loop
            for (yr = 0; yr < nyrs; yr++) {                         // Year loop
              for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop

                Term = 1.0e-10 + H_1(rsp, ksp) * (Type(1) + H_1a(rsp) * H_1b(rsp) / (Type(r_age) + H_1b(rsp) + Type(1.0e-10))); // Eq. 15 Kinzey and Punt (2009)

                // Observed species
                if(ksp < nspp){
                  for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop

                    // Predator-prey ratios
                    Pred_ratio = (N_pred_yrs(rsp, r_age, yr) + Type(1.0e-10)) / (N_pred_eq(rsp, r_age) + Type(1.0e-10)); // Eq. 17 Kinzey and Punt (2009): Predator biomass relative to equilibrium
                    Prey_ratio = (N_prey_yrs(ksp, k_age, yr) + Type(1.0e-10)) / (N_prey_eq(ksp, k_age) + Type(1.0e-10)); // Eq. 16 Prey Kinzey and Punt (2009): biomass relative to equilibrium
                    Pred_r(rsp, r_age, yr) = Pred_ratio;
                    Prey_r(ksp, k_age, yr) = Prey_ratio;

                    switch (msmMode) {
                    case 2: // Holling Type I (linear)
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term;
                      break;
                    case 3: // Holling Type II
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                        ( 1 + H_2(rsp, ksp) * Prey_ratio + 1.0e-10);
                      break;
                    case 4: // Holling Type III
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 +
                        Term * (1 + H_2(rsp, ksp)) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) /
                          (1 + H_2(rsp, ksp) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) + 1.0e-10 );
                      break;
                    case 5: // predator interference
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                        ( 1 + H_2(rsp, ksp) * Prey_ratio + H_3(rsp, ksp) * (Pred_ratio - 1) + 1.0e-10);
                      break;
                    case 6: // predator preemption
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                        ( (1 + H_2(rsp, ksp) * Prey_ratio) * (1 + H_3(rsp, ksp) * (Pred_ratio - 1)) + Type(1.0e-10));
                      break;
                    case 7: // Hassell-Varley
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (2 + H_2(rsp, ksp) + 1.0e-10) /
                        (1.0 + H_2(rsp, ksp) * Prey_ratio + pow((Prey_ratio + Type(1.0e-10)), H_4(rsp, ksp)) + 1.0e-10 );
                      break;
                    case 8: // Ecosim
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term /
                        (1 + H_3(rsp, ksp) * (Pred_ratio - 1 + 1.0e-10));
                      break;
                    default:
                      error("Invalid 'msmMode'");
                    }
                  }
                }
                // "other" is linear
                else{
                  pred_resp(rsp, ksp, r_age, 0, yr) = 1.0e-10 + Term;
                }
              } // end of r_ages, k_ages loop
            }   // =========================
          }
        }


        // 8.2.7. Mortality as a function of predator age: Eq. 7 Kinzey and Punt (2009)
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for (ksp = 0; ksp  <  nspp; ksp++) {
              for (r_age = 0; r_age < nages(rsp); r_age++) {
                for (k_age = 0; k_age < nages(ksp); k_age++) {
                  pred_effect = pred_resp(rsp, ksp, r_age, k_age, yr) * suit_main(rsp , ksp, r_age, k_age, yr);
                  Vmort_ua(rsp, ksp, r_age, k_age, yr) = pred_effect * NByage(rsp, r_age, yr);
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
            for (ksp = 0; ksp < nspp; ksp++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                for (r_age = 0; r_age < nages(rsp); r_age++) {
                  Pmort_ua(rsp, ksp, k_age, yr) += Vmort_ua(rsp, ksp, r_age, k_age, yr);
                }
                M2(ksp, k_age, yr) += Pmort_ua(rsp, ksp, k_age, yr);
              }
            }
          }
        }


        // 8.2.9. Numbers eaten (of modeled prey species); Equations 8 and 9 from Kinzey and Punt 2009
        eaten_la.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (ksp = 0; ksp  <  nspp; ksp++) {
            for (k_age = 0; k_age < nages(ksp); k_age++) {
              // Relative number
              NS_Z = NByage(ksp, k_age, yr) * (1 - exp(-Zed(ksp, k_age, yr) )) / Zed(ksp, k_age, yr);                                  // Eq. 8 Kinzey and Punt (2009) Baranov

              for (rsp = 0; rsp  <  nspp; rsp++) {
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  eaten_ua(rsp, ksp, r_age, k_age, yr) = Vmort_ua(rsp, ksp, r_age, k_age, yr) * NS_Z;                                  // Eq. 8 Kinzey and Punt (2009)

                  for (r_ln = 0; r_ln  <  nlengths(rsp); r_ln++) { // FIXME: Change to diet length bins
                    eaten_la(rsp, ksp, r_ln, k_age, yr) += eaten_ua(rsp, ksp, r_age, k_age, yr)  * age_trans_matrix(r_age, r_ln, pop_alk_index(rsp)); // Eq. 9 Kinzey and Punt (2009)
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
  int yr_ind = yr;
}
if(yr >= nyrs_hind){
  int yr_ind = nyrs_hind - 1;
}

          for (rsp = 0; rsp  < nspp; rsp++) {
            for (ksp = 0; ksp  <  nspp + 1; ksp++) {

              // Species included
              if (ksp < nspp) {
                // Results by length: Eqn 10 kinda Kinzey and Punt (2009)
                for (r_ln = 0; r_ln  <  nlengths(rsp); r_ln++) { // FIXME: Change to diet bins
                  for (k_age = 0; k_age < nages(ksp); k_age++) {
                    Q_mass_l(rsp, ksp, r_ln, yr) += eaten_la(rsp, ksp, r_ln, k_age, yr) * wt(yr_ind, k_age, pop_wt_index(ksp));
                  }
                }
                // Results by k_age: Eq. 11 Kinzey and Punt (2009)
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  for (k_age = 0; k_age < nages(ksp); k_age++) {
                    Q_mass_u(rsp, ksp, r_age, yr) += eaten_ua(rsp, ksp, r_age, k_age, yr) * wt(yr_ind, k_age, pop_wt_index(ksp));
                  }
                }
              }

              // Other food
              else {
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  pred_effect = pred_resp(rsp, ksp, r_age, 0, yr); // Other food k_age = 0 because 1 age is used
                  Tmort = pred_effect * NByage(rsp, r_age, yr);
                  Q_mass_u(rsp, ksp, r_age, yr)  = Q_other_u(rsp, r_age) * (Type(1) - exp(-Tmort));                     // Eq. 11 Kinzey and Punt (2009)
                }
                for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) { // FIXME: Change to diet bins
                  for (r_age = 0; r_age < nages(rsp); r_age++) {
                    Q_mass_l(rsp, ksp, r_ln, yr) += Q_mass_u(rsp, ksp, r_age, yr) * age_trans_matrix(r_age, r_ln, pop_alk_index(rsp)); // Eq. 10 Kinzey and Punt (2009)
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
          for (r_age = 0; r_age < nages(rsp); r_age++) {
            // Calculate year-specific values
            numer = 0;
            denom = 0;
            for (yr = 0; yr < nyrs; yr++) {
              // Average abundance
              n_avg = 1.0e-10 + NByage(rsp, r_age, yr) * sqrt(S(rsp, r_age, yr));  // FIXME: Add average N added 1.0e-10 -dhk June 24 08. not in NRM tpl -dhk apr 28 09
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

        // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
      } // End 8.2. Kinzey predation
      // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
    } // End 8. Predation mortality
    // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //

    // ------------------------------------------------------------------------- //
    // 9. SURVEY COMPONENTS EQUATIONS                                            //
    // ------------------------------------------------------------------------- //
    // Good above here

    // 9.0. Empirical selectivity
    srv_sel.setZero();
    for (int sel_ind = 0; sel_ind < srv_emp_sel_obs.rows(); sel_ind++){

      srv = srv_emp_sel_ctl(sel_ind, 0) - 1;            // Temporary survey index
      sp = srv_emp_sel_ctl(sel_ind, 1) - 1;             // Temporary index of species
      srv_yr = srv_emp_sel_ctl(sel_ind, 2) - styr;      // Temporary index for years of data

      if(srv_yr < nyrs_hind){

        for (age = 0; age < srv_emp_sel_obs.cols(); age++) {
          srv_sel(srv, age, srv_yr) = srv_emp_sel_obs(sel_ind, age);
        }
      }

      // FIXME - have survey selectivity after nyrs_hind filled with last year
    }
    // Good above here


    // 9.1. Estimated survey selectivity
    avgsel_srv.setZero();
    srv_sel_tmp.setZero();
    for(ctl_ind = 0; ctl_ind < n_srv; ctl_ind++){

      srv = srv_control(ctl_ind, 1) - 1;            // Temporary survey index
      sp = srv_control(ctl_ind, 2) - 1;             // Temporary index of species

      // 9.1.1. Logisitic selectivity
      if (srv_sel_type(srv) == 1) {
        for (age = 0; age < nages(sp); age++){
          srv_sel_tmp(ctl_ind, age) = 1 / (1 + exp( -srv_sel_slp(0, ctl_ind) * ((age + 1) - srv_sel_inf(0, ctl_ind))));
        }
      }

      // 9.1.2. Non-parametric selectivity fit to age ranges. NOTE: This can likely be improved
      if (srv_sel_type(srv) == 2) {
        for (age = 0; age < srv_nselages(srv); age++) {
          srv_sel_tmp(ctl_ind, age) = srv_sel_coff(ctl_ind, age);
          avgsel_srv(ctl_ind) +=  exp(srv_sel_coff(ctl_ind, age));
        }
        // 9.1.3 Average selectivity up to nselages(sp)
        avgsel_srv(ctl_ind) = log(avgsel_srv(ctl_ind) /srv_nselages(srv));

        // 9.1.4. Plus group selectivity
        for (age = srv_nselages(srv); age < nages(sp); age++) {
          srv_sel_tmp(ctl_ind, age) = srv_sel_tmp(ctl_ind, srv_nselages(srv) - 1);
        }

        // 9.1.5. Average selectivity across all ages
        avgsel_tmp = 0; // set to zero
        for (age = 0; age < nages(sp); age++) {
          avgsel_tmp += exp(srv_sel_tmp(ctl_ind, age));
        }
        avgsel_tmp = log(avgsel_tmp / nages(sp));

        // 9.1.6. Standardize selectivity
        for (age = 0; age < nages(sp); age++) {
          srv_sel_tmp(ctl_ind, age) -=  avgsel_tmp;
          srv_sel_tmp(ctl_ind, age) = exp(srv_sel_tmp(ctl_ind, age));
        }
      }


      // 9.1.2. Double logistic (Dorn and Methot 1990)
      if (srv_sel_type(srv) == 3) {
        for (age = 0; age < nages(sp); age++){
          srv_sel_tmp(ctl_ind, age) = (1 / (1 + exp( -srv_sel_slp(0, ctl_ind) * ((age + 1) - srv_sel_inf(0, ctl_ind))))) * // Upper slop
            (1 - (1 / (1 + exp( -srv_sel_slp(1, ctl_ind) * ((age + 1) - srv_sel_inf(1, ctl_ind))))));  // Downward slope;
        }
      }

      // Standardize estimated selectivities and add to vector
      if (srv_sel_type(srv) > 0) {
        vector<Type> sel_sp = srv_sel_tmp.row(ctl_ind); // Can't max a matrix....
        srv_sel_tmp.row(ctl_ind) /= max(sel_sp); // Standardize so max sel = 1 for each species

        for (age = 0; age < nages(sp); age++){
          for (yr = 0; yr < nyrs; yr++) {
            srv_sel(srv, age, yr) = srv_sel_tmp(ctl_ind, age);
          }
        }
      }
    } // End loop
    // Good above here


    // -- 9.2. Survey Biomass
    for(srv_ind = 0; srv_ind < srv_biom_ctl.rows(); srv_ind++){

      srv = srv_biom_ctl(srv_ind, 0) - 1;            // Temporary survey index
      sp = srv_biom_ctl(srv_ind, 1) - 1;             // Temporary index of species
      sex = srv_biom_ctl(srv_ind, 2);                // Temporary index for years of data
      srv_yr = srv_biom_ctl(srv_ind, 3) - styr;      // Temporary index for years of data

      mo = srv_biom_n(srv_ind, 0);                    // Temporary index for month

      srv_bio_hat(srv_ind) = 0;                      // Initialize

      if(srv_yr < nyrs_hind){

        for (age = 0; age < nages(sp); age++) {
          // Weight
          if(srv_units(srv) == 1){
            srv_bio_hat(srv_ind) += NByage(sp, age, srv_yr) * exp( - (mo/12) * Zed(sp, age, srv_yr)) * srv_sel(srv, age, srv_yr) * srv_q(srv) * wt(srv_yr, age, srv_wt_index(srv));
          }
          // Numbers
          if(srv_units(srv) == 2){
            srv_bio_hat(srv_ind) += NByage(sp, age, srv_yr) * exp( - (mo/12) * Zed(sp, age, srv_yr)) * srv_sel(srv, age, srv_yr) * srv_q(srv);
          }
        }
      }
    }


    // -- 9.3. Survey composition
    srv_age_obs_hat.setZero();
    srv_comp_hat.setZero();
    for(comp_ind = 0; comp_ind < srv_comp_hat.rows(); comp_ind++){

      srv = srv_comp_ctl(comp_ind, 0) - 1;            // Temporary survey index
      sp = srv_comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
      sex = srv_comp_ctl(comp_ind, 2);                // Temporary index for comp sex (0 = combined, 1 = female, 2 = male)
      srv_comp_type = srv_comp_ctl(comp_ind, 3);      // Temporary index for comp type (0 = age, 1 = length)
      srv_yr = srv_comp_ctl(comp_ind, 4) - styr;      // Temporary index for years of data
      mo = srv_comp_n(comp_ind, 0);                   // Temporary index for month

      if(srv_yr < nyrs_hind){

        srv_hat(comp_ind) = 0;                       // Initialize

        // Total numbers
        for (age = 0; age < nages(sp); age++) {
          srv_age_hat(comp_ind, age ) = NByage(sp, age, srv_yr)  * srv_sel(srv, age, srv_yr) * srv_q(srv) * exp( - Type(mo/12) * Zed(sp, age, srv_yr));
          srv_hat(comp_ind) += srv_age_hat(comp_ind, age );   // Total numbers
        }

        // Add in aging error
        for (int obs_age = 0; obs_age < nages(sp); obs_age++) {
          for (int true_age = 0; true_age < nages(sp); true_age++) {
            srv_age_obs_hat(comp_ind, obs_age) += srv_age_hat(comp_ind, true_age ) * age_error(sp, true_age, obs_age);
          }
        }


        //  Normalize survey catch-at-age
        if (srv_comp_type == 0) {
          for (age = 0; age < nages(sp); age++) {
            srv_comp_hat(comp_ind, age) = srv_age_obs_hat(comp_ind, age) / srv_hat(comp_ind);
          }
        }

        // Convert from catch-at-age to catch-at-length and normalize
        if ( srv_comp_type == 1) {
          for (ln = 0; ln < nlengths(sp); ln++) {
            for (age = 0; age < nages(sp); age++) {
              srv_comp_hat(comp_ind, ln ) += srv_age_obs_hat(comp_ind, age) * age_trans_matrix(age, ln, srv_alk_index(sp));
            }
          }

          // Normalize
          for (ln = 0; ln < nlengths(sp); ln++) {
            srv_comp_hat(comp_ind, ln) = srv_comp_hat(comp_ind, ln) / srv_hat(comp_ind);
          }
        }
      }
    }

    // Good above here

    // ------------------------------------------------------------------------- //
    // 10. FISHERY COMPONENTS EQUATIONS                                          //
    // ------------------------------------------------------------------------- //
    // 10.5. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
    for(fsh_ind = 0; fsh_ind < fsh_biom_ctl.rows(); fsh_ind++){

      fsh = fsh_biom_ctl(fsh_ind, 0) - 1;            // Temporary fishery index
      sp = fsh_biom_ctl(fsh_ind, 1) - 1;             // Temporary index of species
      sex = fsh_biom_ctl(fsh_ind, 2);                // Temporary index for years of data
      fsh_yr = fsh_biom_ctl(fsh_ind, 3) - styr;      // Temporary index for years of data
      mo = fsh_biom_n(fsh_ind, 0);                   // Temporary index for month

      fsh_bio_hat(fsh_ind) = 0;                  // Initialize

      if(fsh_yr < nyrs_hind){

        for (age = 0; age < nages(sp); age++) {
          // By weight
          if(fsh_units(fsh) == 1){
            fsh_bio_hat(fsh_ind) += F(fsh, age, fsh_yr) / Zed(sp, age, fsh_yr) * (1 - exp(-Zed(sp, age, fsh_yr))) * NByage(sp, age, fsh_yr) * wt(fsh_yr, age, fsh_wt_index(fsh)); // 5.5.
          }

          // By numbers
          if(fsh_units(fsh) == 2){
            fsh_bio_hat(fsh_ind) += F(fsh, age, fsh_yr) / Zed(sp, age, fsh_yr) * (1 - exp(-Zed(sp, age, fsh_yr))) * NByage(sp, age, fsh_yr); // 5.5.
          }
        }
      }
    }
    // Good above here


    // -- 9.3. Fishery composition
    fsh_age_obs_hat.setZero();
    fsh_comp_hat.setZero();
    for(comp_ind = 0; comp_ind < fsh_comp_hat.rows(); comp_ind++){

      fsh = fsh_comp_ctl(comp_ind, 0) - 1;            // Temporary fishery index
      sp = fsh_comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
      sex = fsh_comp_ctl(comp_ind, 2);                // Temporary index for comp sex (0 = combined, 1 = female, 2 = male)
      fsh_comp_type = fsh_comp_ctl(comp_ind, 3);      // Temporary index for comp type (0 = age, 1 = length)
      fsh_yr = fsh_comp_ctl(comp_ind, 4) - styr;      // Temporary index for years of data
      mo = fsh_comp_n(comp_ind, 0);                   // Temporary index for month

      fsh_hat(comp_ind) = 0;                          // Initialize

      if(fsh_yr < nyrs_hind){

        // Total numbers
        for (age = 0; age < nages(sp); age++) {
          fsh_age_hat(comp_ind, age ) = F(fsh, age, fsh_yr) / Zed(sp, age, fsh_yr) * (1 - exp(-Zed(sp, age, fsh_yr))) * NByage(sp, age, fsh_yr); // 5.4.
          fsh_hat( comp_ind ) += fsh_age_hat(comp_ind, age );   // Total numbers
        }


        // Adjust for aging error
        for (int obs_age = 0; obs_age < nages(sp); obs_age++) {
          for (int true_age = 0; true_age < nages(sp); true_age++) {

            fsh_age_obs_hat(comp_ind, obs_age) += fsh_age_hat(comp_ind, true_age ) * age_error(sp, true_age, obs_age);
          }
        }

        //  Survey catch-at-age
        if (fsh_comp_type == 0) {
          for (age = 0; age < nages(sp); age++) {
            fsh_comp_hat(comp_ind, age) = fsh_age_obs_hat(comp_ind, age ) / fsh_hat(comp_ind);
          }
        }

        // Convert from catch-at-age to catch-at-length
        if ( fsh_comp_type == 1) {
          for (ln = 0; ln < nlengths(sp); ln++) {
            for (age = 0; age < nages(sp); age++) {
              fsh_comp_hat(comp_ind, ln ) += fsh_age_obs_hat(comp_ind, age) * age_trans_matrix(age, ln, fsh_alk_index(fsh));
            }
          }

          for (ln = 0; ln < nlengths(sp); ln++) {
            fsh_comp_hat(comp_ind, ln ) = fsh_comp_hat(comp_ind, ln) / fsh_hat(comp_ind);
          }
        }
      }
    }
    // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //
  } // End iterations
  // - END LOOP - END LOOP - END LOOP - END LOOP - END LOOP - //


  // ------------------------------------------------------------------------- //
  // 11. LIKELIHOOD EQUATIONS                                                  //
  // Get max number of fisheries or surveys
  vector<int> n_tot_vec(2);
  n_tot_vec(0) = n_srv; n_tot_vec(1) = n_fsh;
  int n_tot = imax(n_tot_vec);

  // 11.0. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(17, n_tot); jnll_comp.setZero();  // matrix of negative log-likelihood components


  // -- Survey components
  // Slot 0 -- Survey biomass
  // Slot 1 -- Survey age/length composition
  // Slot 2 -- Survey selectivity
  // Slot 3 -- Survey selectivity normalization
  // -- Fishery components
  // Slot 4 -- Total catch (kg)
  // Slot 5 -- Fishery age/length composition
  // Slot 6 -- Fishery selectivity
  // Slot 7 -- Fishery selectivity normalization
  // -- Priors/penalties
  // Slot 8 -- Tau -- Annual recruitment deviation
  // Slot 9 -- init_dev -- Initial abundance-at-age
  // Slot 10 -- Epsilon -- Annual fishing mortality deviation
  // -- M2 likelihood components
  // Slot 13 -- Ration likelihood
  // Slot 14 -- Ration penalties
  // Slot 15 -- Diet weight likelihood
  // Slot 16 -- Stomach content of prey length ln in predator length a likelihood


  // 11.1. OFFSETS AND PENALTIES
  // 11.1.1 -- Set up offset objects
  matrix<Type> offset(n_tot, 2); offset.setZero(); // Offset for multinomial likelihood

  vector<Type> offset_diet_w(nspp); offset_diet_w.setZero(); // Offset for total stomach content weight likelihood
  vector<Type> offset_diet_l(nspp); offset_diet_l.setZero(); // Offset for total stomach content of prey length ln in predator length a
  vector<Type> offset_uobsagewt(nspp); offset_uobsagewt.setZero(); // Offset for stomach proportion by weight likelihood


  // 11.1.1. -- Fishery age comp offsets
  for (comp_ind = 0; comp_ind < fsh_comp_obs.rows(); comp_ind++) {

    fsh = fsh_comp_ctl(comp_ind, 0) - 1;            // Temporary survey index
    sp = fsh_comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
    fsh_yr = fsh_comp_ctl(comp_ind, 4);             // Temporary index for years of data

    // Add years from hindcast only
    if(fsh_yr <= endyr){
      for (ln = 0; ln < fsh_comp_obs.cols(); ln++) {
        if(!isNA( fsh_comp_obs(comp_ind, ln ))){
          offset(fsh, 1) -= Type(fsh_comp_n(comp_ind, 1)) * (fsh_comp_obs(comp_ind, ln) + MNConst) * log(fsh_comp_obs(comp_ind, ln) + MNConst ) ;
        }
      }
    }
  }


  // 11.1.2. -- Survey age comp offsets
  for (comp_ind = 0; comp_ind < srv_comp_obs.rows(); comp_ind++) {

    srv = srv_comp_ctl(comp_ind, 0) - 1;            // Temporary survey index
    sp = srv_comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
    srv_yr = srv_comp_ctl(comp_ind, 4);             // Temporary index for years of data


    // Only use years from hindcast only
    if(srv_yr <= endyr){
      for (ln = 0; ln < srv_comp_obs.cols(); ln++) {
        if(!isNA( srv_comp_obs(comp_ind, ln ))){
          offset(srv, 0) -= Type(srv_comp_n(comp_ind, 1)) * (srv_comp_obs(comp_ind, ln) + MNConst) * log(srv_comp_obs(comp_ind, ln) + MNConst ) ;
        }
      }
    }
  }


  // 11.1.4. -- Offsets for diet (weights)
  for (rsp = 0; rsp < nspp; rsp++) {
    for (yr = 0; yr < nyrs_hind; yr++) {
      for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) {
        // if (stoms_w_N(rsp,r_ln,yr) > 0){ Sample size
        for (ksp = 0; ksp < nspp; ksp++) { // FIXME: no offset for "other food"
          if ( diet_w_sum(rsp, ksp, r_ln, yr) > 0) {
            offset_diet_w(rsp) -= stom_tau(rsp) * (diet_w_sum(rsp, ksp, r_ln, yr) + MNConst) * log(diet_w_sum(rsp, ksp, r_ln, yr) + MNConst); // FIXME add actual sample size
          }
        }
      }
    }
  }

  // 11.1.5. -- Offsets for diet (lengths) FIXME: add real stomach sample size
  for (rsp = 0; rsp < nspp; rsp++) {
    for (ksp = 0; ksp < nspp; ksp++) {
      for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) {
        for (k_ln = 0; k_ln < nlengths(ksp); k_ln++) {
          if (Uobs(rsp, ksp, r_ln, k_ln) > 0) {
            offset_diet_l(rsp) -= stom_tau(rsp) * Uobs(rsp, ksp, r_ln, k_ln) * log(Uobs(rsp, ksp, r_ln, k_ln) + 1.0e-10);
          }
        }
      }
    }
  }

  // 11.1.6. -- Offsets for diet proportion by weight
  // FIXME - change years to adjust for missing years of diet data
  // If 5D array
  if(a2_dim.size() == 5){
    for (rsp = 0; rsp < nspp; rsp++) {
      for (yr = 0; yr < nyrs_hind; yr++) {
        for (r_age = 0; r_age < nages(rsp); r_age++) {
          for (ksp = 0; ksp < nspp; ksp++) {
            for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
              //if (UobsWtAge(rsp, ksp, r_age, k_age, yr) > 0) { // (rsp, ksp, a, r_ln, yr)
              offset_uobsagewt(rsp) -= stom_tau(rsp) * (UobsWtAge(rsp, ksp, r_age, k_age, yr) + MNConst) * log(UobsWtAge(rsp, ksp, r_age, k_age, yr) + MNConst);
              //}
            }
          }
        }
      }
    }
  }

  // If 4D array
  if(a2_dim.size() == 4){
    for (rsp = 0; rsp < nspp; rsp++) {
      for (r_age = 0; r_age < nages(rsp); r_age++) {
        for (ksp = 0; ksp < nspp; ksp++) {
          for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
            //if (UobsWtAge(rsp, ksp, r_age, k_age) > 0) { // (rsp, ksp, a, r_ln, yr)
            offset_uobsagewt(rsp) -= stom_tau(rsp) * (UobsWtAge(rsp, ksp, r_age, k_age) + MNConst) * log(UobsWtAge(rsp, ksp, r_age, k_age) + MNConst);
            //}
          }
        }
      }
    }
  }


  // 11.2. FIT OBJECTIVE FUNCTION
  // Slot 0 -- Survey biomass -- NFMS annual BT survey and EIT survey
  for (srv_ind = 0; srv_ind < srv_biom_obs.rows(); srv_ind++) {

    srv = srv_biom_ctl(srv_ind, 0) - 1;            // Temporary survey index
    sp = srv_biom_ctl(srv_ind, 1) - 1;             // Species is the 2nd column
    srv_yr = srv_biom_ctl(srv_ind, 3);             // Temporary index for years of data

    // Add years from hindcast
    if(srv_yr <= endyr){
      jnll_comp(0, srv) += pow(log(srv_biom_obs(srv_ind, 0)) - log(srv_bio_hat(srv_ind)), 2) / (2 * square(srv_biom_obs(srv_ind, 1))); // NOTE: This is not quite the lognormal and biohat will be the median.
    }
  }



  // Slot 1 -- Survey age composition
  for (comp_ind = 0; comp_ind < srv_comp_obs.rows(); comp_ind++) {

    srv = srv_comp_ctl(comp_ind, 0) - 1;            // Temporary survey index
    sp = srv_comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
    srv_yr = srv_comp_ctl(comp_ind, 4);             // Temporary index for years of data


    // Only use years wanted
    if(srv_yr <= endyr){
      for (ln = 0; ln < srv_comp_obs.cols(); ln++) {
        if(!isNA( srv_comp_obs(comp_ind, ln) )){
          jnll_comp(1, srv) -= Type(srv_comp_n(comp_ind, 1)) * (srv_comp_obs(comp_ind, ln) + MNConst) * log(srv_comp_hat(comp_ind, ln) + MNConst ) ;
        }
      }
    }
  }

  // Remove offsets
  for (srv = 0; srv < n_srv; srv++) {
    // srv_yr = yrs_srv_age(sp, yr) - styr;
    jnll_comp(1, srv) -= offset(srv, 0);
  }


  // Slot 2 -- Survey selectivity
  for(srv = 0; srv < n_srv; srv++){ // Loop around surveys
    jnll_comp(2, srv) = 0; // FIXME: Likeliy redundant
    sp = srv_spp(srv);

    if (srv_sel_type(srv) == 2) {

      for (age = 0; age < (nages(sp) - 1); age++) {
        if ( srv_sel(srv, age, 0) > srv_sel(srv, age + 1, 0)) {
          jnll_comp(2, sp) += 20 * pow( log(srv_sel(srv, age, 0) / srv_sel(srv, age + 1, 0) ), 2);
        }
      }

      // Extract only the selectivities we want
      vector<Type> sel_tmp(nages(sp));

      for (age = 0; age < nages(sp); age++) {
        sel_tmp(age) = log(srv_sel(srv, age, 0));  sel_tmp.setZero();
      }

      for (age = 0; age < nages(sp) - 2; age++) {

        sel_tmp(age) = first_difference( first_difference( sel_tmp ) )(age);
        jnll_comp(2, srv) += curv_pen_fsh * pow( sel_tmp(age) , 2); // FIX

      }
    }
  }


  // Slot 3 -- Add survey selectivity normalization
  for(srv = 0; srv < n_srv; srv++){
    jnll_comp(3, srv) += 50 * square(avgsel_srv(srv));
  }


  // Slot 4 -- Total catch -- Fishery observer dat
  for (fsh_ind = 0; fsh_ind < fsh_biom_obs.rows(); fsh_ind++) {

    fsh = fsh_biom_ctl(fsh_ind, 0) - 1;            // Temporary fishery index
    sp = fsh_biom_ctl(fsh_ind, 1) - 1;             // Species is the column 3
    fsh_yr = fsh_biom_ctl(fsh_ind, 3);             // Temporary index for years of data

    // Ad only years from hindcast
    if(fsh_yr <= endyr){
      jnll_comp(4, fsh) += pow(log(fsh_biom_obs(fsh_ind, 0)) - log(fsh_bio_hat(fsh_ind)), 2) / (2 * square(fsh_biom_obs(fsh_ind, 1))); // NOTE: This is not quite the log  normal and biohat will be the median.
    }
  }



  // Slot 5 -- Fishery age composition -- Fishery observer data
  for (comp_ind = 0; comp_ind < fsh_comp_obs.rows(); comp_ind++) {

    fsh = fsh_comp_ctl(comp_ind, 0) - 1;            // Temporary survey index
    sp = fsh_comp_ctl(comp_ind, 1) - 1;             // Temporary index of species
    fsh_yr = fsh_comp_ctl(comp_ind, 4);             // Temporary index for years of data

    // Add years from hindcast only
    if(fsh_yr <= endyr){
      for (ln = 0; ln < fsh_comp_obs.cols(); ln++) {
        if(!isNA( fsh_comp_obs(comp_ind, ln) )){
          jnll_comp(5, fsh) -= Type(fsh_comp_n(comp_ind, 1)) * (fsh_comp_obs(comp_ind, ln) + MNConst) * log(fsh_comp_hat(comp_ind, ln) + MNConst ) ;
        }
      }
    }
  }

  // Remove offsets
  for (fsh = 0; fsh < n_fsh; fsh++) {
    // srv_yr = yrs_srv_age(sp, yr) - styr;
    jnll_comp(5, fsh) -= offset(fsh, 1);
  }



  // Slot 6 -- Fishery selectivity
  for(fsh = 0; fsh < n_fsh; fsh++){ // Loop around surveys
    jnll_comp(6, fsh) = 0; // FIXME: Likeliy redundant
    sp = fsh_spp(fsh);

    if (fsh_sel_type(fsh) == 2) {

      for (age = 0; age < (nages(sp) - 1); age++) {
        //if ( fsh_sel(fsh, age, 0) > fsh_sel(fsh, age + 1, 0)) {
        if ( fsh_sel(fsh, age, 0) > fsh_sel(fsh, age + 1, 0)) {
          //jnll_comp(6, sp) += 20 * pow( log(fsh_sel(fsh, age, 0) / fsh_sel(fsh, age + 1, 0) ), 2);
          jnll_comp(6, sp) += 20 * pow( log(fsh_sel(fsh, age, 0) / fsh_sel(fsh, age + 1, 0) ), 2);
        }
      }

      // Extract only the selectivities we want
      vector<Type> sel_tmp(nages(sp)); sel_tmp.setZero();

      for (age = 0; age < nages(sp); age++) {
        sel_tmp(age) = log(fsh_sel(fsh, age, 0));
      }

      for (age = 0; age < nages(sp) - 2; age++) {

        sel_tmp(age) = first_difference( first_difference( sel_tmp ) )(age);
        jnll_comp(6, fsh) += curv_pen_fsh * pow( sel_tmp(age) , 2); // FIX

      }
    }
  }


  // Slot 7 -- Add fishery selectivity normalization
  for (fsh = 0; fsh < n_fsh; fsh++) {
    jnll_comp(7, fsh) = 50 * square(avgsel_fsh(fsh));
  }

  // Slots 8-10 -- PRIORS: PUT RANDOM EFFECTS SWITCH HERE
  for (sp = 0; sp < nspp; sp++) {
    // Slot 10 -- init_dev -- Initial abundance-at-age
    for (age = 1; age < nages(sp); age++) {

      if (random_rec == 0) {
        jnll_comp(8, sp) += pow( init_dev(sp, age - 1) - Type(0.25), 2);
      }

      if (random_rec == 1) {
        jnll_comp(8, sp) -= dnorm( init_dev(sp, age - 1) - square(r_sigma(sp)) / 2, Type(0.0), r_sigma(sp), true);
      }

    }

    for (yr = 0; yr < nyrs_hind; yr++) {

      // Slot 11 -- Tau -- Annual recruitment deviation
      if (random_rec == 0) {
        jnll_comp(9, sp) += pow( rec_dev(sp, yr) - Type(0.25), 2);    // Recruitment deviation using penalized likelihood. ADDED lognormal bias correction
      }
      if (random_rec == 1) {
        jnll_comp(9, sp) -= dnorm( rec_dev(sp, yr) - square(r_sigma(sp)) / 2, Type(0.0), r_sigma(sp), true);    // Recruitment deviation using random effects.
      }

      // Slot 12 -- Epsilon -- Annual fishing mortality deviation
      jnll_comp(10, sp) += pow( F_dev(sp, yr), 2);      // Fishing mortality deviation using penalized likelihood.
    }
  }



  // 11.3. Diet likelihood components from MSVPA
  if ((msmMode == 1) & (suitMode > 0)) {
    // Slot 14 -- Diet weight likelihood
    // If 5D array
    if(a2_dim.size() == 5){
      for (rsp = 0; rsp < nspp; rsp++) {
        jnll_comp(15, rsp) = 0;
        // FIXME - change diet years for years with missing data
        for (yr = 0; yr < nyrs_hind; yr++) {
          for (r_age = 0; r_age < nages(rsp); r_age++) {
            for (ksp = 0; ksp < nspp; ksp++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
                //if (UobsWtAge(rsp, ksp, r_age, k_age, yr) > 0) { // (rsp, ksp, a, r_ln, yr)
                jnll_comp(15, rsp) -= stom_tau(rsp) * (UobsWtAge(rsp, ksp, r_age, k_age, yr) + MNConst) * (log(UobsWtAge_hat(rsp, ksp, r_age, k_age, yr) + MNConst));
                // }
              }
            }
          }
        }
        jnll_comp(15, rsp) -= offset_uobsagewt(rsp);
      }
    }

    // If 4D array
    if(a2_dim.size() == 4){
      for (rsp = 0; rsp < nspp; rsp++) {
        jnll_comp(15, rsp) = 0;
        for (r_age = 0; r_age < nages(rsp); r_age++) {
          for (ksp = 0; ksp < nspp; ksp++) {
            for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
              //if (UobsWtAge(rsp, ksp, r_age, k_age) > 0) { // (rsp, ksp, a, r_ln, yr)
              jnll_comp(15, rsp) -= stom_tau(rsp) * (UobsWtAge(rsp, ksp, r_age, k_age) + MNConst) * (log(mn_UobsWtAge_hat(rsp, ksp, r_age, k_age) + MNConst));
              //}
            }
          }
        }
        jnll_comp(15, rsp) -= offset_uobsagewt(rsp);
      }
    }
  } // End diet proportion by weight component

  // 11.4. Diet likelihood components from Kinzey and Punt
  if ((msmMode > 1)) {
    // Slot 13 -- Ration likelihood
    for (yr = 0; yr < nyrs_hind; yr++) {
      for (sp = 0; sp < nspp; sp++) {
        for (age = 1; age < nages(sp); age++) { // don't include age zero in likelihood
          jnll_comp(13, sp) += 0.5 * pow( log( omega_hat(sp, age, yr) + 1.0e-10) -
            log(ration2Age(sp, age, yr)), 2) / (sd_ration * sd_ration); // FIXME: add year indices for ration
        }
      }
    }


    // Slot 14 -- Ration penalalties FIXME: not sure if necessary: talk to Andre
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        Type mean_ohat = 0;
        for (yr = 0; yr < nyrs; yr++) {
          mean_ohat += omega_hat(sp, age, yr) / nyrs;
        }
        for (yr = 0; yr < nyrs; yr++) {
          jnll_comp(14, sp) += 20 *  pow(omega_hat(sp, age, yr) - mean_ohat, 2);
        }
      }
    }


    // Slot 14 -- Diet weight likelihood
    for (rsp = 0; rsp < nspp; rsp++) {
      for (yr = 0; yr < nyrs_hind; yr++) {
        for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) {
          for (ksp = 0; ksp < nspp; ksp++) {                  // FIME: need to add in other food
            if (diet_w_sum(rsp, ksp, r_ln, yr) > 0) { // (rsp, ksp, a, r_ln, yr)
              jnll_comp(15, rsp) -= stom_tau(rsp) * diet_w_sum(rsp, ksp, r_ln, yr) * log(Q_hat(rsp, ksp, r_ln, yr) + 1.0e-10);
            }
          }
        }
      }
      jnll_comp(15, rsp) -= offset_diet_w(rsp);
    }


    Type Denom = 0;

    // Calculate the predicted fraction by length-class (Eqn 17)
    T_hat.setZero();
    for (rsp = 0; rsp < nspp; rsp++) {
      for (ksp = 0; ksp < nspp; ksp++) {
        vector<Type> eaten_lmy(nlengths(ksp)); eaten_lmy.setZero(); // no. of prey@length eaten by a predator length during iyr

        for (r_ln = 0; r_ln < nlengths(rsp); r_ln++) { //
          // This is Equation 17
          for (yr = 0; yr < nyrs_hind; yr++) { // FIXME: loop by stomach year
            // int stom_yr = yrs_stomlns(rsp, ksp, stm_yr); // FIXME: index by stomach year
            for (k_ln = 0; k_ln < nlengths(ksp); k_ln ++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                eaten_lmy(k_ln) += eaten_la(rsp, ksp, r_ln, k_age, yr) * age_trans_matrix(k_age, k_ln, pop_alk_index(ksp));
              }
              T_hat(rsp, ksp, r_ln, k_ln) += stom_tau(rsp) * eaten_lmy(k_ln); // FIXME: add actual stomach content sample size
            }
          }


          Denom = 0;
          for (k_ln = 0; k_ln < nlengths(ksp); k_ln++) {
            Denom += T_hat(rsp, ksp, r_ln, k_ln);
          }

          // Renormalize the eaten vector
          for (k_ln = 0; k_ln < nlengths(ksp); k_ln++) {
            T_hat(rsp, ksp, r_ln, k_ln) /= Denom;

            // Likelihood of diet length         / This is equation 16
            if (Uobs(rsp, ksp, r_ln, k_ln) > 0) {
              jnll_comp(16, rsp) -= stom_tau(rsp) * Uobs(rsp, ksp, r_ln, k_ln) * log(T_hat(rsp, ksp, r_ln, k_ln)  + 1.0e-10);
            }
          }
        }
      }
      jnll_comp(16, rsp) -= offset_diet_l(rsp);
    }
  } // End if statement for diet likelihood


  // ------------------------------------------------------------------------- //
  // 12. REPORT SECTION                                                        //
  // ------------------------------------------------------------------------- //

  // 12.1. Population components
  REPORT( mn_rec );
  REPORT( Zed );
  REPORT( NByage );
  REPORT( AvgN );
  REPORT( S );
  REPORT( biomassByage );
  REPORT( biomassSSBByage );
  REPORT( biomass );
  REPORT( biomassSSB );
  REPORT( pmature );
  REPORT(r_sigma);
  REPORT( R );
  REPORT( M1 );
  REPORT( M );

  ADREPORT( Zed );
  ADREPORT( biomass );
  ADREPORT( biomassSSB );
  ADREPORT( r_sigma );
  ADREPORT( R );
  ADREPORT( M1);

  // -- 12.2. Fishery components
  REPORT( avgsel_fsh );
  REPORT( fsh_sel );
  REPORT( fsh_sel2 );
  REPORT( catch_hat );
  REPORT( F );
  REPORT( F_tot );
  REPORT( fsh_bio_hat );
  REPORT( fsh_hat );
  REPORT( fsh_age_hat );
  REPORT( fsh_comp_hat );
  REPORT( fsh_comp_obs );
  REPORT( fsh_sel_type );

  // -- 12.3. Survey components
  REPORT( avgsel_srv );
  REPORT( srv_sel );
  REPORT( srv_sel_tmp );
  REPORT(  srv_bio_hat );
  REPORT( srv_hat );
  REPORT( srv_age_hat );
  REPORT( srv_comp_hat );


  // -- 12.4. Likelihood components
  REPORT( jnll_comp );
  REPORT( offset );
  REPORT( offset_diet_w );
  REPORT( offset_diet_l );
  REPORT( offset_uobsagewt );


  // -- 12.5. Ration components
  REPORT( ConsumAge );
  REPORT( LbyAge );
  REPORT( mnWt_obs );
  REPORT( fT );
  REPORT( TempC );
  REPORT( ration2Age );


  // 12.6. Suitability components
  REPORT( suma_suit );
  REPORT( suit_main );
  REPORT( suit_other );
  REPORT( stom_div_bio2 );
  REPORT( stomKir );
  REPORT( avail_food );
  REPORT( of_stomKir );
  REPORT( M2 );
  REPORT( B_eaten );
  REPORT( UobsWtAge_hat );
  REPORT( mn_UobsWtAge_hat );

  // -- 12.7. Kinzey predation functions
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

  // ------------------------------------------------------------------------- //
  // 13. END MODEL                                                             //
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
