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

  // -- 2.3.2 BT Survey Components
  DATA_IVECTOR( nyrs_srv_biom );  // Number of years of survey biomass data; n = [nspp]
  DATA_IMATRIX( yrs_srv_biom );   // Years of survey biomass data; n = [nspp, nyrs_srv_biom]
  DATA_MATRIX( srv_biom );         // Observed BT survey biomass (kg); n = [nspp, nyrs]
  DATA_MATRIX( srv_biom_se );     // Observed annual biomass error (SE); n = [nspp, nyrs_srv_biom]

  // -- 2.3.3. BT Survey age components
  DATA_IVECTOR( nyrs_srv_age );   // Number of years of survey age/length composition; n = [1, nspp]
  DATA_IMATRIX( yrs_srv_age );    // Years for the survey age/length composition data; n = [nspp, nyrs_srv_age]
  DATA_IVECTOR( srv_age_type );  // Type of compisition (1 = age; 2 = length); n = [1, nspp]
  DATA_IVECTOR( srv_age_bins );   // Number of size binds for the age/length comps; n = [nspp]
  DATA_MATRIX( srv_age_n );       // Sample size for the multinomial; n = [nspp, nyrs_srv_age]
  DATA_ARRAY( srv_age_obs) ;      // Observed BT age comp; n = [nspp, nages, nyrs]
  DATA_MATRIX( srv_age_sizes );  // Observed size composition
  DATA_ARRAY( age_trans_matrix); // observed sp_age/size compositions; n = [nspp, nages, srv_age_bins]

  // -- 2.3.4 EIT Survey components
  DATA_INTEGER( n_eit);          // Number of years with EIT data; n = [1]
  DATA_IVECTOR( yrs_eit);        // Years for available EIT data; n = [n_eit]
  DATA_VECTOR( eit_age_n);       // Number  of  EIT Hauls for multinomial; n = [yrs_eit]
  DATA_MATRIX( obs_eit_age);     // Observed EIT age comp; n = [1, nages, nyrs] NOTE: may need to change this for future
  DATA_VECTOR( obs_eit);         // Observed EIT survey biomass (kg); n = [1, nyrs]
  DATA_MATRIX( eit_sel);         // Observed EIT survey selectivity; n = [eit_age, nyrs_eit_sel]

  // -- 2.3.5. Weight-at-age
  DATA_IVECTOR( nyrs_wt_at_age );  // Number of years of weight at age data; n = [nspp]
  DATA_IMATRIX( yrs_wt_at_age );   // Years of weight-at-age; data n = [nspp, nyrs_wt_at_age]
  DATA_ARRAY( wt );                // Weight-at-age by year; n = [nyrs_wt_at_age, nages, nspp]

  // -- 2.3.6 Other
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

  // 2.5. DERIVED QUANTITIES # Calculate these in the model
  // DATA_MATRIX(d);               // VBGF allometric slope of consumption (d); n = [nspp, nyrs]
  // DATA_MATRIX(Winf);            // VBGF max asymptoptic weight; n = [nspp, nyrs]

  // 2.6. FIXED PARAMETERS
  // -- 2.6.1. von Bertalannfy growth function (VBGF)
  // DATA_VECTOR(d0_vbgf);      // VBGF intercept for d parameter; n = [1, nspp]
  // DATA_MATRIX(d_dev_vbgf);   // Annual deviation for VBGF d parameter; n = [nspp, nyrs] # NOTE: Need to figure out how to best vectorize this
  // DATA_VECTOR(Beta_d_vbgf);  // Temperature covariate for VBGF d parameter; n = [1, nspp]
  // DATA_VECTOR( logK );       // VBGF energy loss constant (kg kg^-1 yr^-1); n[1, nspp]
  // DATA_VECTOR( logH );       // VBGF assimilation constant (kg kg^-1 yr^-1); n[1, nspp]
  // DATA_VECTOR( t0 );         // VBGF age at Weight 0 (yr); n[1, nspp]

  // -- 2.6.2. Others
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

  int max_age = imax(nages);    // Integer of maximum nages to make the arrays.
  int max_bin = imax(srv_age_bins); // Integer of maximum number of length/age bins.

  array<Type>  fsh_age_obs(nyrs, max_bin, nspp);  // Observed fishery age comp; n = [nyrs_fsh_comp, fsh_age_bins, nspp]
  matrix<Type> tc_obs(nspp, nyrs); // Observed total catch (n); n = [nspp, nyrs] NOTE: This may not be necessary if loading data from tmp

  matrix<Type> srv_biom_lse(nspp, nyrs); // Observed annual biomass CV; n = [nspp, nyrs_srv_biom]
  srv_biom_lse = srv_biom_se.array()/ srv_biom.array();          // CV estimation
  srv_biom_lse = pow( log( ( pow( srv_biom_lse.array(), Type(2) ).array() + 1).array()).array(), Type(0.5));

  matrix<Type> M1( nspp, max_age); M1.setZero();
  M1 = M1_base.array() + Type(0.0001);

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
  // -- 4.2.1. Fishery observations
  matrix<Type>  tc_biom_hat(nspp, nyrs);              // Estimated total yield (kg); n = [nspp, nyrs]
  array<Type>   catch_hat(nyrs, max_age, nspp);       // Estimated catch-at-age (n); n = [nspp, nages, nyrs]
  array<Type>   fsh_age_hat(nyrs, max_bin, nspp);     // Estimated fishery age comp; n = [nspp, nages, nyrs]
  matrix<Type>  tc_hat(nspp, nyrs);                   // Estimated total catch (n); n = [nspp, nyrs]
  array<Type>   F(nyrs, max_age, nspp);               // Estimated fishing mortality; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_sel(nspp, max_age); fsh_sel.setZero(); // Log estimated fishing selectivity
  vector<Type>  avgsel_fsh(nspp); avgsel_fsh.setZero();    // Average fishery selectivity

  // -- 4.2.2. BT Survey components
  matrix<Type>  srv_bio_hat(nspp, nyrs);              // Estimated BT survey biomass (kg); n = [nspp, nyrs]
  array<Type>   srv_age_hat(nyrs, max_bin, nspp);     // Estimated BT age comp; n = [nspp, nages, nyrs]
  matrix<Type>  srv_hat(nspp, nyrs);                  // Estimated BT survey total abundance (n); n = [nspp, nyrs]
  matrix<Type>  srv_sel(nspp, max_age); srv_sel.setZero(); // Estimated survey selectivity at age; n = [nspp, nyrs]
  vector<Type>  avgsel_srv(nspp); avgsel_srv.setZero();    // Average survey selectivity

  // -- 4.2.3. EIT Survey Components
  matrix<Type>  eit_age_hat(nyrs, srv_age_bins(0));               // Estimated EIT age comp; n = [nyrs, 12 ages]
  vector<Type>  eit_hat(nyrs);                       // Estimated EIT survey biomass (kg); n = [nyrs]

  // -- 4.2.4. Estimated population parameters
  matrix<Type>  R(nspp, nyrs);                        // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type>   NByage(nyrs + 1, max_age, nspp);          // Numbers at age; n = [nspp, nages, nyrs]
  array<Type>   biomassByage(nyrs, max_age, nspp);    // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type>  biomass(nspp, nyrs);                  // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type>  biomassSSB(nspp, nyrs);               // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type>   biomassSSBByage(nyrs, max_age, nspp); // Spawning biomass at age (kg); n = [nspp, nages, nyrs]

  // -- 4.3. Parameter transformations
  Type eit_q = exp(log_eit_q);                        // EIT Catchability
  array<Type>   Zed(nyrs, max_age, nspp);             // Total mortality at age; n = [nspp, nages, nyrs]
  array<Type>   S(nyrs, max_age, nspp);               // Survival at age; n = [nspp, nages, nyrs]


  // ------------------------------------------------------------------------- //
  // 5. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0

    // 7.1. ESTIMATE FISHERY SELECTIVITY
  for (i = 0; i <nspp; i++){
    for(j=0; j < nselages; j++){
      fsh_sel(i, j) = fsh_sel_coff(i, j);
      avgsel_fsh(i) +=  exp(fsh_sel_coff(i, j));
    }
    // 7.1.3 Average selectivity up to nselages
    avgsel_fsh(i) = log(avgsel_fsh(i) / nselages);

    // 7.1.4. Plus group selectivity
    for(j = nselages; j < nages(i); j++){
      fsh_sel(i, j) = fsh_sel(i, nselages - 1);
    }

    // 7.1.5. Average selectivity across all ages
    Type avgsel_tmp = 0; // Temporary object for average selectivity across all ages
    for(j = 0; j < nages(i); j++){
      avgsel_tmp += exp(fsh_sel(i, j));
    }
    avgsel_tmp = log(avgsel_tmp / nages(i));

    // 7.1.6. Standardize selectivity
    for(j = 0; j < nages(i); j++){
      fsh_sel(i, j) -=  avgsel_tmp;
      fsh_sel(i, j) = exp(fsh_sel(i, j));
    }
  }


  // 7.2. ESTIMATE FISHING MORTALITY
  for(i=0; i<nspp; i++){
    for(y=0; y<nyrs; y++){
      for(j=0; j<nages(i); j++){
        F(y, j, i) = fsh_sel(i, j) * exp(ln_mean_F(i) + F_dev(i, y));
      }
    }
  }


  // 7.4. Estimate total mortality at age NOTE: May need to go above population dynamics
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        Zed(y, j, i) = M1(i, j) + F(y, j, i);
        S(y, j, i) = exp(-Zed(y, j, i));
      }
    }
  }

  // 5.1. ESTIMATE RECRUITMENT T1.1
  for(i=0; i< nspp; i++){
    for(y=0; y < nyrs; y++){
      R(i,y) = exp(ln_mn_rec(i) + rec_dev(i,y));
      NByage(y, 0, i) = R(i, y);
    }
    NByage(nyrs, 0, i) = R(i, nyrs-1);
  }


  // 5.2. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
   
      if((j > 0) & (j < nages(i) - 1)){
        NByage(0, j, i) = exp(ln_mn_rec(i) - (j)*M1(i,j) + init_dev(i, j));
      }
      // -- 5.2.2. Where y = 1 and j > Ai.
      if(j == (nages(i)-1)){
        NByage(0, j, i) = exp(ln_mn_rec(i) - (j)*M1(i,j) + init_dev(i, j))/ (1-exp(-M1(i, nages(i)-1))); // NOTE: This solves for the geometric series
      }
    }
  }


  // 5.3. ESTIMATE NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
  biomass.setZero();  // Initialize Biomass
  biomassSSB.setZero(); // Initialize SSB
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=1; y < nyrs+1; y++){
        // -- 5.3.1.  Where 1 <= j < Ai
        if(j < (nages(i)-1)){
          NByage(y, j+1, i) = NByage(y-1, j, i) * S(y-1, j, i);
        }

        // -- 5.3.2. Plus group where j > Ai. NOTE: This is not the same as T1.3 because I used j = A_i rather than j > A_i.
        if(j == (nages(i)-1)){ 
          NByage(y, nages(i)-1, i) = NByage(y-1, nages(i)-2, i) * S(y-1, nages(i)-2, i) + NByage(y-1, nages(i)-1, i) * S(y-1, nages(i)-1, i);
        }
        // NOTE: The "-1" is because the NByage has the initial population numbers
        biomassByage(y-1, j, i) = NByage(y, j, i) * wt(y-1, j, i); // 5.5.
        biomassSSBByage(y-1, j, i) = biomassByage(y-1, j, i) * pmature(i, j); // 5.6.

        // -- 5.3.3. Estimate Biomass and SSB
        biomass(i, y-1) += biomassByage(y-1, j, i);
        biomassSSB(i, y-1) += biomassSSBByage(y-1, j, i);
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 6. SURVEY COMPONENTS EQUATIONS                                            //
  // ------------------------------------------------------------------------- //
  // 6.1. Survey selectivity
  for (i=0; i<nspp; i++){
    // 6.1.1. Logisitic selectivity
    if(logist_sel_phase(i) > 0){
      for (j=0; j < nages(i); j++)
        srv_sel(i, j) = 1/ (1 + exp( -srv_sel_coff(i, 0) * ((j+1) - srv_sel_coff(i, 1))));
    }

    // 6.1.2. Selectivity fit to age ranges. NOTE: This can likely be improved
    if(logist_sel_phase(i) < 0){
      for(j=0; j < nselages; j++){
        srv_sel(i, j) = srv_sel_coff(i, j);
        avgsel_srv(i) +=  exp(srv_sel_coff(i, j));
      }
      // 6.1.3 Average selectivity up to nselages
      avgsel_srv(i) = log(avgsel_srv(i) / nselages);

      // 6.1.4. Plus group selectivity
      for(j = nselages; j < nages(i); j++){
        srv_sel(i, j) = srv_sel(i, nselages - 1);
      }

      // 6.1.5. Average selectivity across all ages
      Type avgsel_tmp = 0; // Temporary object for average selectivity across all ages
      for(j = 0; j < nages(i); j++){
        avgsel_tmp += exp(srv_sel(i, j));
      }
      avgsel_tmp = log(avgsel_tmp / nages(i));

      // 6.1.6. Standardize selectivity
      for(j = 0; j < nages(i); j++){
        srv_sel(i, j) -=  avgsel_tmp;
        srv_sel(i, j) = exp(srv_sel(i, j));
      }
    }
  }

  // 6.2 EIT Components
  // -- 6.2.1 EIT Survey Biomass
  int eit_yr_ind;
  eit_age_hat.setZero();
  eit_hat.setZero();
  for(j=0; j < nages(0); j++){
    for(y=0; y < n_eit; y++){
      eit_yr_ind = yrs_eit(y) - styr;

      eit_age_hat(eit_yr_ind, j) = NByage(eit_yr_ind+1, j, 0) * eit_sel(y, j) * eit_q; // Remove the mid-year trawl?
      eit_hat(eit_yr_ind) += eit_age_hat(eit_yr_ind, j) * wt(eit_yr_ind, j, 0);  //
    }
  }


  // -- 6.2.2 EIT Survey Age Composition
  for(j=0; j < nages(0); j++){
    for (y=0; y < n_eit; y++){
      eit_yr_ind = yrs_eit(y) - styr;
      eit_age_hat(eit_yr_ind, j) = eit_age_hat(eit_yr_ind, j) / eit_age_hat.row(eit_yr_ind).sum(); // Divide numbers at age by total numbers for each year
    }
  }


  // 6.3 BT Components
  // -- 6.3.1 BT Survey Biomass
  int srv_yr_ind;
  srv_age_hat.setZero();
  srv_hat.setZero();
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs_srv_biom(i); y++){

        srv_yr_ind = yrs_srv_biom(i, y) - styr; // Temporary index for years of data

        srv_age_hat(srv_yr_ind, j, i) = NByage(srv_yr_ind+1, j, i) * exp(-0.5 * Zed(srv_yr_ind, j, i)) * srv_sel(i, j) * exp(log_srv_q(i));
        srv_hat(i, srv_yr_ind) += srv_age_hat(srv_yr_ind, j, i);   // Total numbers
        srv_bio_hat(i, srv_yr_ind) += srv_age_hat(srv_yr_ind, j, i) * wt(srv_yr_ind, j, i);  //
      }
    }
  }


  // -- 6.4.2 BT Survey Age Composition: NOTE: will have to alter if age comp data are not the same length as survey biomass data
  vector<Type> srv_age_tmp( imax(nages)); // Temporary vector of survey-catch-at-age for matrix multiplication
  
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_srv_age(i); y++){

      srv_yr_ind = yrs_srv_age(i, y) - styr; // Temporary index for years of data

      // 6.4.2.1 -- BT Survey catch-at-age
      if(srv_age_type(i)==1){
        for(j=0; j < nages(i); j++){
          srv_age_hat(srv_yr_ind, j, i) = srv_age_hat(srv_yr_ind, j, i) / srv_hat(i, srv_yr_ind);
        }
      }
      // 6.4.2.2 -- Convert from catch-at-age to catch-at-length: NOTE: There has got to be a better way
      if(srv_age_type(i)!=1){
        for(j=0; j < nages(i); j++){
          srv_age_tmp(j) = srv_age_hat(srv_yr_ind, j, i);
        }

        matrix<Type> ALK = trim_matrix( matrix_from_array(age_trans_matrix, i), nages(i), srv_age_bins(i) );
        vector<Type> srv_age_tmp_trimmed = trim_vector(srv_age_tmp, nages(i) );
        vector<Type> srv_len_tmp = vec_mat_prod( srv_age_tmp_trimmed , ALK ); // Multiply the ALK for species i against the survey catch-at-age for year y

        for(j=0; j < srv_age_bins(i); j++){
          srv_age_hat(srv_yr_ind, j, i) = srv_len_tmp(j) / srv_hat(i, srv_yr_ind) ; // * age_trans_matrix.col().col(i)) / srv_hat(i, y); // # NOTE: Double check the matrix algebra here
        }
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 7. FISHERY COMPONENTS EQUATIONS                                           //
  // ------------------------------------------------------------------------- //


  // NOTE: The above may need to be before the population dynamics


  // 7.5. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
  int fsh_yr_ind;
  tc_hat.setZero();
  tc_biom_hat.setZero(); // Initialize tc_biom_hat
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs_tc_biom(i); y++){

        fsh_yr_ind = yrs_tc_biom(i, y) - styr; // Temporary index for years of data

        catch_hat(fsh_yr_ind, j, i) = F(fsh_yr_ind, j, i)/Zed(fsh_yr_ind, j, i) * (1 - exp(-Zed(fsh_yr_ind, j, i))) * NByage(fsh_yr_ind + 1, j, i); // 5.4.
        tc_hat(i, fsh_yr_ind) += catch_hat(fsh_yr_ind, j, i); // Estimate catch in numbers
        tc_biom_hat(i, fsh_yr_ind) += catch_hat(fsh_yr_ind, j, i) * wt(fsh_yr_ind, j, i); // 5.5.
      }
    }
  }


  // 7.6. ESTIMATE FISHERY AGE COMPOSITION
  vector<Type> fsh_age_tmp( imax(nages)); // Temporary vector of survey-catch-at-age for matrix multiplication
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_fsh_comp(i); y++){

      fsh_yr_ind = yrs_fsh_comp(i, y) - styr; // Temporary index for years of data

      /// 7.7.2.1 -- Estimate age composition of the fishery
      if(fsh_age_type(i)==1){
        for(j=0; j < nages(i); j++){
          fsh_age_hat(fsh_yr_ind, j, i) = catch_hat(fsh_yr_ind, j, i) / tc_hat(i, fsh_yr_ind);
        }
      }

      // 7.7.2.1 -- Convert from catch-at-age to catch-at-length: NOTE: There has got to be a better way
      if(fsh_age_type(i)!=1){
        for(j=0; j < nages(i); j++){
          fsh_age_tmp(j) = catch_hat(fsh_yr_ind, j, i);
        }

        matrix<Type> ALK = trim_matrix( matrix_from_array(age_trans_matrix, i), nages(i), fsh_age_bins(i) );
        vector<Type> fsh_age_tmp_trimmed = trim_vector(fsh_age_tmp, nages(i) );
        vector<Type> fsh_len_tmp = vec_mat_prod( fsh_age_tmp_trimmed , ALK ); // Multiply the ALK for species i against the survey catch-at-age for year y

        for(j=0; j < fsh_age_bins(i); j++){
          fsh_age_hat(fsh_yr_ind, j, i) = fsh_len_tmp(j) / tc_hat(i, fsh_yr_ind) ; // * age_trans_matrix.col().col(i)) / srv_hat(i, y); // # NOTE: Double check the matrix algebra here
        }
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 8. LIKELIHOOD EQUATIONS                                                   //
  // ------------------------------------------------------------------------- //
  // 8.0. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(13,nspp); // matrix of negative log-likelihood components
  jnll_comp.setZero();
  Type jnll = 0;

  // -- Data components
  // Slot 0 -- BT survey biomass -- NFMS annual BT survey
  // Slot 1 -- BT survey age composition -- NFMS annual BT survey
  // Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
  // Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
  // Slot 4 -- Total catch -- Fishery observer data
  // Slot 5 -- Fishery age composition -- Fishery observer data
  // -- Likelihood penalties
  // Slot 6 -- Fishery selectivity
  // Slot 7 -- Fishery selectivity normalization
  // Slot 8 -- Survey selectivity
  // Slot 9 -- Survey selectivity normalization
  // -- Priors
  // Slot 10 -- Tau -- Annual recruitment deviation
  // Slot 11 -- init_dev -- Initial abundance-at-age
  // Slot 12 -- Epsilon -- Annual fishing mortality deviation


  // 8.1. OFFSETS AND PENALTIES
  // 8.1.1 -- Set up offset objects
  vector<Type> offset_srv(nspp); offset_srv.setZero(); // Offset for multinomial likelihood
  vector<Type> offset_fsh(nspp); offset_srv.setZero(); // Offset for multinomial likelihood
  Type offset_eit = 0; // Offset for multinomial likelihood
  tc_obs.setZero(); // Set total catch to 0 to initialize

  for(i=0; i < nspp; i++){

    // 8.1.1. -- Fishery age comp offsets
    for(y = 0; y < nyrs_fsh_comp(i); y++){
        for(j=0; j < fsh_age_bins(i); j++){
        tc_obs(i, y) += obs_catch(y, j, i);
        }
        for(j=0; j < fsh_age_bins(i); j++){
         fsh_age_obs(y, j, i) = obs_catch(y, j, i)/(tc_obs(i, y) + 0.01); // Calculate age comp
         offset_fsh( i ) -= tau * (fsh_age_obs(y, j, i) + MNConst) * log(fsh_age_obs(y, j, i) + MNConst);
        }
      }


    for (y = 0; y < nyrs_srv_age(i); y++){
    // 8.1.2. -- Survey age comp offsets
        for(j = 0; j < srv_age_bins(i); j++){
          matrix<Type> srv_age_tmp = trim_matrix( matrix_from_array(srv_age_obs, i), nyrs_srv_age(i), srv_age_bins(i) ); // Get srv_age_ob from array for species i
          srv_age_obs(y, j, i) = srv_age_obs(y, j, i) / srv_age_tmp.row(y).sum() ; // Convert numbers-at-age to age comp
      offset_srv(i) -= srv_age_n(i, y)*(srv_age_obs(y, j, i) + MNConst) * log(srv_age_obs(y, j, i) + MNConst ) ;
    }
  }
  }


  // 8.1.3. -- Offsets for acoustic survey
  for(y = 0; y < n_eit; y++){
    for(j = 0; j < srv_age_bins(0); j++){
      obs_eit_age(y, j) = obs_eit_age.row(y).sum(); // Convert from catch-at-age to age comp
      offset_eit -= eit_age_n(y) * (obs_eit_age(y, j) + MNConst) * log(obs_eit_age(y, j) + MNConst );
    }
  }


  // 8.2. FIT OBJECTIVE FUNCTION
  // Slot 0 -- BT survey biomass -- NFMS annual BT survey
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_srv_biom(i); y++){
      srv_yr_ind = yrs_srv_biom(i, y) - styr;
      jnll_comp(0, i) += pow(log(srv_biom(i, y)) - log(srv_bio_hat(i, srv_yr_ind)), 2) / (2 * pow(srv_biom_lse(i, y), 2)); // NOTE: This is not quite the lognormal and biohat will be the median.
    }
  }


  // Slot 1 -- BT survey age composition -- NFMS annual BT survey
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for (y=0; y < nyrs_srv_age(i); y++){
        srv_yr_ind = yrs_srv_age(i, y) - styr;
        jnll_comp(1, i) -= srv_age_n(i, y) * (srv_age_obs(y, j, i) + MNConst) * log(srv_age_hat(srv_yr_ind, j, i) + MNConst); // Should srv_age_obs  be in log space?
      }
    }
    jnll_comp(1, i) -= offset_srv(i);
  }


  // Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
  for (y=0; y < n_eit; y++){
    eit_yr_ind = yrs_eit(y) - styr;
    jnll_comp(2, 0) += 12.5 * pow(log(obs_eit(y)) - log(eit_hat(eit_yr_ind) + 1.e-04), 2);
  }


  // Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
  for(j=0; j < nages(0); j++){
    for (y=0; y < n_eit; y++){
      eit_yr_ind = yrs_eit(y) - styr;
      jnll_comp(3, 0) -= eit_age_n(y) *  (obs_eit_age(y, j) + MNConst) * log(eit_age_hat(eit_yr_ind, j) + MNConst);
    }
  }
  jnll_comp(3, 0) -= offset_eit;


  // Slot 4 -- Total catch -- Fishery observer data
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_tc_biom(i); y++){
      fsh_yr_ind = yrs_tc_biom(i, y) - styr;
      jnll_comp(4,i) += pow((log(tc_biom_hat(i, fsh_yr_ind) + Type(1.e-4)) - log(tcb_obs(i, y) + Type(1.e-4))), 2) / (2 * pow(sigma_catch, 2)); // T.4.5
    }
  }


  // Slot 5 -- Fishery age composition -- Fishery observer data
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for (y=0; y < nyrs_fsh_comp(i); y++){
        fsh_yr_ind = yrs_fsh_comp(i, y) - styr; // Temporary index for years of data
        // int yr_ind  = yrs_fsh_comp(i, y) - styr;
        jnll_comp(5, i) -= tau * (fsh_age_obs(y, j, i) + MNConst) * log(fsh_age_hat(fsh_yr_ind, j, i) + MNConst );
      }
    }
    jnll_comp(5, i) -= offset_fsh(i);
  }


  // Slot 6 -- Fishery selectivity
  for(i=0; i < nspp; i++){
    for(j=0; j < (nages(i) - 1); j++){
      if( fsh_sel(i, j) > fsh_sel( i, j + 1 )){
        jnll_comp(6, i) -= 20 * pow( log( fsh_sel(i, j) / fsh_sel(i, j + 1 ) ), 2);
      }
    }

    // Extract only the selectivities we want
    vector<Type> sel_tmp(nages(i));
    for(j=0; j < nages(i); j++){
      sel_tmp(j) = log(fsh_sel(i, j));
    }

    for(j=0; j < nages(i) - 2; j++){
      sel_tmp(j) = first_difference( first_difference( sel_tmp ) )(j);
      jnll_comp(6, i) += curv_pen_fsh * pow( sel_tmp(j) , 2); // FIX
    }
  }


  // Slot 7 -- Fishery selectivity normalization
  for(i=0; i < nspp; i++){
    jnll_comp(7, i) += 50 * pow(avgsel_fsh(i), 2);
  }


  // Slot 8 -- Survey selectivity
  for(i=0; i < nspp; i++){
    if (logist_sel_phase(i) < 0){
      // Extract only the selectivities we want
      vector<Type> sel_tmp(nages(i));
      for(j=0; j < nages(i); j++){
        sel_tmp(j) = log(srv_sel(i, j));
      }
      for(j=0; j < nages(i) - 2; j++){
        sel_tmp(j) = first_difference( first_difference( sel_tmp ) )(j);
        jnll_comp(8, i) += curv_pen_fsh * pow( sel_tmp(j) , 2); // FIX
      }
    }
  }


  // Slot 9 -- Survey selectivity normalization
  for(i=0; i < nspp; i++){
    jnll_comp(9, i) += 50 * pow(avgsel_srv(i), 2);
  }


  // Slots 10-12 -- PRIORS: PUT RANDOM EFFECTS SWITCH HERE
  for(i=0; i < nspp; i++){
    // Slot 9 -- init_dev -- Initial abundance-at-age
    for(j=0; j < nages(i); j++){
      jnll_comp(11, i) += pow( init_dev(i,j), 2);
    }

    // Slot 8 -- Tau -- Annual recruitment deviation
    // Slot 10 -- Epsilon -- Annual fishing mortality deviation
    for (y=0; y < nyrs; y++){
      jnll_comp(10, i) += pow( rec_dev(i,y), 2);     // Recruitment deviation using penalized likelihood.
      jnll_comp(12, i) += pow( F_dev(i,y), 2);       // Fishing mortality deviation using penalized likelihood.
    }
  }

  // ------------------------------------------------------------------------- //
  // 9. SIMULATION SECTION                                                     //
  // ------------------------------------------------------------------------- //


  // ------------------------------------------------------------------------- //
  // 10. REPORT SECTION                                                        //
  // ------------------------------------------------------------------------- //
  REPORT( pmature );
  REPORT( R );
  REPORT( M1 );
  REPORT( avgsel_srv );
  REPORT( srv_sel );
  REPORT( avgsel_fsh );
  REPORT( fsh_sel );
  REPORT( eit_hat );
  REPORT( eit_age_hat );
  REPORT( NByage );
  REPORT( S );
  REPORT( biomassByage );
  REPORT( biomassSSBByage );
  REPORT( biomass );
  REPORT( biomassSSB );
  REPORT( srv_bio_hat );
  REPORT( srv_hat );
  REPORT( srv_age_hat );
  REPORT( F );
  REPORT( F_dev );
  REPORT( tc_biom_hat );
  REPORT( catch_hat );
  REPORT( tc_hat );
  REPORT( fsh_age_hat );
  REPORT( Zed );
  REPORT( srv_age_obs );

  // ------------------------------------------------------------------------- //
  // END MODEL                                                                 //
  // ------------------------------------------------------------------------- //
  jnll = jnll_comp.sum();
  REPORT( jnll );
  REPORT( jnll_comp );
  return jnll;
}
