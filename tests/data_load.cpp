#include <TMB.hpp>
// ------------------------------------------------------------------------- //
//                 CEATTLE version 3.1.2                                     //
//                  Template Model Builder                                   //
//               Multispecies Statistical Model                              //
//          Bioenergetic-based Assessment for Understanding                  //
//              Biomass Linkages To The Environment                          //
//                  in the Bering Sea                                        //
//                   Oct. 2017                                               //
//                                                                           //
// AUTHORS:   Kirstin Holsman, Jim Ianelli                                   //
//            Modified by Grant Adams                                        //
// CITATIONS:                                                                //
// 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
// 0. LOAD DEPENDENCIES                                                      //
// ------------------------------------------------------------------------- //
#include "../src/include/functions.hpp"


template<class Type>
Type objective_function<Type>::operator() (){
  // ------------------------------------------------------------------------- //
  // 1. MODEL CONFIGURATION                                                    //
  // ------------------------------------------------------------------------- //
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(debug);      // Logical vector to debug or not
  DATA_INTEGER(msmMode);
  //    0 = run in single species mode                                         
  //    1 = run in MSM mode 
                                                     
  DATA_IVECTOR(logist_sel_phase); // Selectivity for BT survey 
  //    0 = fit to data
  //    1 = logistic

  // 1.2. Temporal dimensions
  DATA_INTEGER( nyrs );       // Number of estimation years
  DATA_INTEGER( styr );       // Start year

  // 1.3. Number of species
  DATA_INTEGER( nspp );       // Number of species (prey)


  // 1.4. MODEL OBJECTS
  // 1.4.1. LOOPING INDICES -- k = observation, i = species/prey, j = age/prey age (yr), y = year, p = predator, a = predator age (yr)
  int  i, j, y; //, k, p, a;

  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                           //
  // ------------------------------------------------------------------------- //

  // 2.1. FIXED VALUES
  int tau = 200; // Fishery age composition sample size
  Type MNConst = 0.001;      // Constant additive for logistic functions
  int nselages = 8;
  Type curv_pen_fsh = 12.5; // Fishery selectivity penalty
  Type sigma_catch = 0.05;      // SD of catch


  // 2.2. DIMENSIONS OF DATA
  int endyr = nyrs + styr;   // End year

  // -- 2.2.2. Species attributes
  DATA_IVECTOR( nages );      // Number of species (prey) ages
  // DATA_INTEGER( n_pred);  // Number of predator species
  // DATA_IVECTOR( n_age_pred); // Number of predator ages

  // 2.3. DATA INPUTS (i.e. assign data to objects)
  // -- 2.3.1 Fishery Components
  DATA_IVECTOR( nyrs_tc_biom );// Number of years with total observed catch; n = [nspp]
  DATA_IMATRIX( yrs_tc_biom ); // Years with total observed catch; n = [nspp, nyrs_tc_biom]
  DATA_MATRIX( tcb_obs );      // Observed total yield (kg); n = [nspp, nyrs_tc_biom]

  DATA_IVECTOR( nyrs_fsh_comp );  // Number of years in the fishery sp_age composition data; n = [nspp]
  DATA_IMATRIX( yrs_fsh_comp );   // Years for the fishery sp_age composition data; n = [nspp, nyrs_fsh_comp]
  DATA_IVECTOR( fsh_age_type );   // Which method of calculating fishery age hat (2 = ATF); n = [nspp]
  DATA_IVECTOR( fsh_age_bins );   // Bins for fishery age composition data; n = [nspp]
  DATA_ARRAY( obs_catch );        // Observed fishery catch-at-age or catch-at-length; n = [nspp, fsh_age_bins, nyrs_fsh_comp]

  // Calculate in house
  array<Type>  fsh_age_obs(imax(nages), nyrs, nspp);  // Observed fishery age comp; n = [nspp, fsh_age_bins, nyrs_fsh_comp]
  matrix<Type> tc_obs(nspp, nyrs); // Observed total catch (n); n = [nspp, nyrs] NOTE: This may not be necessary if loading data from tmp

  // -- 2.3.2 BT Survey Components
  DATA_IVECTOR( nyrs_srv_biom );  // Number of years of survey biomass data; n = [nspp]
  DATA_IMATRIX( yrs_srv_biom );   // Years of survey biomass data; n = [nspp, nyrs_srv_biom]
  DATA_MATRIX( srv_biom );         // Observed BT survey biomass (kg); n = [nspp, nyrs]
  DATA_MATRIX( srv_biom_se );     // Observed annual biomass error (SE); n = [nspp, nyrs_srv_biom]
  matrix<Type> srv_biom_lse(nspp, imax(nyrs_srv_biom)); // Observed annual biomass CV; n = [nspp, nyrs_srv_biom]
  srv_biom_lse = srv_biom_se.array()/ srv_biom.array();          // CV estimation
  srv_biom_lse = pow( log( ( pow( srv_biom_lse.array(), Type(2) ).array() + 1).array()).array(), Type(0.5));

  DATA_IVECTOR( nyrs_srv_age );   // Number of years of survey age/length composition; n = [nspp]
  DATA_IMATRIX( yrs_srv_age );    // Years for the survey age/length composition data; n = [nspp, nyrs_srv_age]
  DATA_IVECTOR( srv_age_type );  // Type of compisition (1 = age; 2 = length); n = [nspp]
  DATA_IVECTOR( srv_age_bins );   // Number of size binds for the age/length comps; n = [nspp]
  DATA_MATRIX( srv_age_n );       // Sample size for the multinomial; n = [nspp, nyrs_srv_age]
  DATA_ARRAY( srv_age_obs) ;      // Observed BT age comp; n = [nspp, nages, nyrs]
  DATA_MATRIX( srv_age_sizes );  // Observed size composition
  DATA_ARRAY( age_trans_matrix); // observed sp_age/size compositions; n = [nspp, nages, srv_age_bins]

  // -- 2.3.3 EIT Survey components
  DATA_INTEGER( n_eit);          // Number of years with EIT data; n = [1]
  DATA_IVECTOR( yrs_eit);        // Years for available EIT data; n = [n_eit]
  DATA_VECTOR( eit_age_n);       // Number  of  EIT Hauls for multinomial; n = [yrs_eit]
  DATA_MATRIX( obs_eit_age);     // Observed EIT age comp; n = [1, nages, nyrs] NOTE: may need to change this for future
  DATA_VECTOR( obs_eit);         // Observed EIT survey biomass (kg); n = [1, nyrs]
  DATA_MATRIX( eit_sel);         // Observed EIT survey selectivity; n = [eit_age, nyrs_eit_sel]

  // -- 2.3.4 Other
  // DATA_VECTOR( TempC );           // Bottom temperature (degrees C); n = [1, nyrs] # NOTE: Need to figure out how to make it flexible for alternative environmental predictors
  // DATA_ARRAY( Diet_Mat );         // Annual gravimetric proportion of prey in predator stomach; n = [n_pred, n_age_pred, nspp, nages, nyrs]
  // DATA_INTEGER( other_food );     // Biomass of other prey (kg); n = [nyrs, n_pred] # QUESTION: Is this year specific?


  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  // DATA_VECTOR( phi_p_bp);        // Annual relative foraging rate (d yr^-1)
  // DATA_VECTOR( CA );             //  Intercept of the allometric maximum consumption function (g g^-1 yr^-1); n = [1, nspp]
  // DATA_VECTOR( CB );             //  Allometric slope of maximum consumption; n = [1, nspp]
  DATA_VECTOR( Tcm );               //  Consumption maximum physiological temperature (degree C); n = [1, n_pred]
  DATA_VECTOR( Tco );               //  Consumption optimum physiological temperature (degree C); n = [1, n_pred]
  // DATA_VECTOR( Qc );             //  Max consumption parameter; n = [1, n_pred]

  // 2.6. DERIVED QUANTITIES # Calculate these in the model
  // DATA_MATRIX(d);               // VBGF allometric slope of consumption (d); n = [nspp, nyrs]
  // DATA_MATRIX(Winf);            // VBGF max asymptoptic weight; n = [nspp, nyrs]

  // 2.7. FIXED PARAMETERS
  // -- 2.7.1. von Bertalannfy growth function (VBGF)
  // DATA_VECTOR(d0_vbgf);      // VBGF intercept for d parameter; n = [1, nspp]
  // DATA_MATRIX(d_dev_vbgf);   // Annual deviation for VBGF d parameter; n = [nspp, nyrs] # NOTE: Need to figure out how to best vectorize this
  // DATA_VECTOR(Beta_d_vbgf);  // Temperature covariate for VBGF d parameter; n = [1, nspp]
  // DATA_VECTOR( logK );       // VBGF energy loss constant (kg kg^-1 yr^-1); n[1, nspp]
  // DATA_VECTOR( logH );       // VBGF assimilation constant (kg kg^-1 yr^-1); n[1, nspp]
  // DATA_VECTOR( t0 );         // VBGF age at Weight 0 (yr); n[1, nspp]

  // -- 2.7.2. Others
  DATA_MATRIX( M1_base );         // Residual natural mortality; n = [nspp, nages]
  DATA_IVECTOR( mf_type );        // Sex specific mort and weight at age? : 1 = same for both, 2 = seperate wt at sp_age for each sex
  DATA_MATRIX( propMorF );        // Proportion-at-age of females of population; n = [nspp, nages]
  DATA_MATRIX( pmature );         // Proportion of mature females at age; [nspp, nages]


  // ------------------------------------------------------------------------- //
  // 2.8. Debugging with data inputs                                           //
  // ------------------------------------------------------------------------- //
  if(debug == 1){
    // -- 2.8.2.1 Check to make sure the first year of survey data are not before start year
    for(int i = 0; i<nspp; i++){
      if(yrs_tc_biom(i,0) < styr){
        std::cerr<<"First year of total catch biomass of species "<< i + 1 << " is before specified start year"<<std::endl;
        return(0);
      }
    }

    // -- 2.8.2.2 Check to make sure the first year of survey data are not before start year
    for(int i = 0; i<nspp; i++){
      if(yrs_srv_biom(i,0) < styr){
        std::cerr<<"First year of survey biomass of species "<< i + 1 << " is before specified start year"<<std::endl;
        return(0);
      }
    }

    // -- 2.8.2.3. Check to make sure the years of survey data biomass and age are the same
    for(int i = 0; i<nspp; i++){
      if(nyrs_srv_biom(i) != nyrs_srv_age(i)){
        std::cerr<<"Nyrs of survey biomass and age-comp of species "<< i + 1 << " do not match"<<std::endl;
        return(0);
      }
    }

    if(yrs_eit(0) < styr){
      std::cerr<<"First year of EIT survey biomass is before specified start year"<<std::endl;
      return(0);
    }
  }


  // ------------------------------------------------------------------------- //
  // 3. INITIAL CALCULATIONS                                                   //
  // ------------------------------------------------------------------------- //
  for ( i = 1; i < nspp ; i++){
    for( j = 0 ; j < nages(i); j++ ){
    pmature( i, j ) = pmature( i, j ) * propMorF(i + (mf_type(i) - 1), j);
  }
  }


  // ------------------------------------------------------------------------- //
  // 4. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //

  // 4.1. PARAMETERS (assign parameters to objects)
  // -- 4.1.1 Recruitment parameters
  PARAMETER_VECTOR( ln_mn_rec );       // Mean recruitment; n = [1, nspp]
  PARAMETER_MATRIX( rec_dev );         // Annual recruitment deviation; n = [nspp, nyrs]
  // PARAMETER(sigma_rec);             // Standard deviation of recruitment variation # NOTE: Have this estimated if using random effects.
  // -- 4.1.2. Abundance parameters
  PARAMETER_MATRIX( init_dev );             // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this
  // -- 4.1.3. fishing mortality parameters
  PARAMETER_VECTOR( ln_mean_F );      // Log mean fishing mortality; n = [1, nspp]
  PARAMETER_MATRIX( F_dev );          // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change
  // PARAMETER_VECTOR(sigma_F); // SD of fishing mortality deviations; n = [1, nspp] # NOTE: Have this estimated if using random effects.
  // -- 4.1.4. Selectivity parameters
  PARAMETER_MATRIX( srv_sel_coff );    // Survey selectivity parameters; n = [nspp, nselages]
  PARAMETER_MATRIX( fsh_sel_coff );    // Fishery age selectivity coef; n = [nspp, nselages]
  PARAMETER( log_eit_q );              // EIT Catchability; n = [1]
  PARAMETER_VECTOR( log_srv_q );       // BT Survey catchability; n = [1, nspp]


  // 4.2. DERIVED QUANTITIES
  int max_age = imax(nages);    // Integer of maximum nages to make the arrays.
  int max_bin = imax(srv_age_bins); // Integer of maximum number of length/age bins.

  // -- 4.2.1. Fishery observations
  matrix<Type>  tc_biom_hat(nspp, nyrs);              // Estimated total yield (kg); n = [nspp, nyrs]
  array<Type>   catch_hat(nyrs, max_age, nspp);       // Estimated catch-at-age (n); n = [nspp, nages, nyrs]
  array<Type>   fsh_age_hat(nyrs, max_age, nspp);     // Estimated fishery age comp; n = [nspp, nages, nyrs]
  matrix<Type>  tc_hat(nspp, nyrs);                   // Estimated total catch (n); n = [nspp, nyrs]
  array<Type>   F(nyrs, max_age, nspp);               // Estimated fishing mortality; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_sel(nspp, max_age); fsh_sel.setZero(); // Log estimated fishing selectivity
  vector<Type>  avgsel_fsh(nspp); avgsel_fsh.setZero();    // Average fishery selectivity

  // -- 4.2.2. BT Survey components
  matrix<Type>  srv_bio_hat(nspp, imax(nyrs_srv_biom));              // Estimated BT survey biomass (kg); n = [nspp, nyrs]
  array<Type>   srv_age_hat(imax(nyrs_srv_age), max_bin, nspp);     // Estimated BT age comp; n = [nspp, nages, nyrs]
  matrix<Type>  srv_hat(nspp, imax(nyrs_srv_biom));                  // Estimated BT survey total abundance (n); n = [nspp, nyrs]
  matrix<Type>  srv_sel(nspp, max_age); srv_sel.setZero(); // Estimated survey selectivity at age; n = [nspp, nyrs]
  vector<Type>  avgsel_srv(nspp); avgsel_srv.setZero();    // Average survey selectivity

  // -- 4.2.3. EIT Survey Components
  array<Type>   Weight_at_Age(nyrs, max_age, nspp);   // Estimated weight-at-age; n = [nspp, nages, nyrs]
  matrix<Type>  eit_age_hat(12, n_eit);               // Estimated EIT age comp; n = [12 ages, nyrs]
  vector<Type>  eit_hat(n_eit);                       // Estimated EIT survey biomass (kg); n = [nyrs]

  // -- 4.2.4. Estimated population parameters
  matrix<Type>  R(nspp, nyrs);                        // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type>   NByage(nyrs, max_age, nspp);          // Numbers at age; n = [nspp, nages, nyrs]
  array<Type>   biomassByage(nyrs, max_age, nspp);    // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type>  biomass(nspp, nyrs);                  // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type>  biomassSSB(nspp, nyrs);               // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type>   biomassSSBByage(nyrs, max_age, nspp); // Spawning biomass at age (kg); n = [nspp, nages, nyrs]

  // -- 4.3. Parameter transformations
  Type eit_q = exp(log_eit_q);                        // EIT Catchability
  array<Type>   Zed(nyrs, max_age, nspp);             // Total mortality at age; n = [nspp, nages, nyrs]


  // ------------------------------------------------------------------------- //
  // 5. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0

  // 5.1. ESTIMATE RECRUITMENT T1.1
  for(i=0; i< nspp; i++){
    for(y=0; y < nyrs; y++){
      R(i,y) = exp(ln_mn_rec(i) + rec_dev(i,y));
      NByage(y, 0, i) = R(i, y);
    }
  }
  REPORT( R );


  // 5.2. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
  matrix<Type> M1( nspp, max_age); M1.setZero();
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      M1(i,j) = M1_base(i,j) + Type(0.0001);
    
      if((j > 0) & (j < nages(i) - 1)){
        NByage(0, j, i) = exp(ln_mn_rec(i) - (j+1)*M1(i,j) + init_dev(i, j));
      }
      // -- 5.2.2. Where y = 1 and j > Ai.
      if(j == (nages(i)-1)){
        NByage(0, j, i) = exp(ln_mn_rec(i) - (j+1)*M1(i,j) + init_dev(i, j))/ (1-exp(-M1(i,nages(i)))); // NOTE: This solves for the geometric series
      }
    }
  }
  REPORT( NByage );
  REPORT( M1 );

return 0;
}

