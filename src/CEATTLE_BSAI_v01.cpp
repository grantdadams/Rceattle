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
// 0. FUNCTIONS                                                              //
// ------------------------------------------------------------------------- //
// Function for getting max of an IVECTOR
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


template<class Type>
Type objective_function<Type>::operator() (){
  // ------------------------------------------------------------------------- //
  // 1. MODEL CONFIGURATION                                                    //
  // ------------------------------------------------------------------------- //
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(debug); // Logical vector to debug or not


  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                           //
  // ------------------------------------------------------------------------- //
  // 2.1. FIXED VALUES
  int tau = 200; // Fishery age composition sample size
  Type MNConst = 0.001;     // Constant additive for logistic functions
  int nselages = 8;


  // 2.2. DIMENSIONS OF DATA
  // -- 2.2.1. Temporal dimensions
  DATA_INTEGER( nyrs);       // Number of estimation years
  DATA_INTEGER( styr);       // Start year
  int endyr = nyrs + styr;  // End year

  // -- 2.2.2. Number of species
  DATA_INTEGER( nspp);       // Number of species (prey)
  DATA_IVECTOR( nages);      // Number of species (prey) ages
  // DATA_INTEGER( n_pred);  // Number of predator species
  // DATA_IVECTOR( n_age_pred); // Number of predator ages

  // 2.3. DATA INPUTS (i.e. assign data to objects)
  // -- 2.3.1 Fishery Components
  DATA_IVECTOR( nyrs_tc_biom_obs );// Number of years with total observed catch; n = [nspp]
  DATA_IMATRIX( yrs_tc_biom_obs ); // Years with total observed catch; n = [nspp, nyrs_tc_biom_obs]
  DATA_MATRIX( tc_biom_obs );      // Observed total yield (kg); n = [nspp, nyrs_tc_biom_obs]

  DATA_IVECTOR( nyrs_fsh_comp );  // Number of years in the fishery sp_age composition data; n = [nspp]
  DATA_IMATRIX( yrs_fsh_comp );   // Years for the fishery sp_age composition data; n = [nspp, nyrs_fsh_comp]
  DATA_IVECTOR( fsh_age_type );   // Which method of calculating fishery age hat (2 = ATF); n = [nspp]
  DATA_IVECTOR( fsh_age_bins );   // Bins for fishery age composition data; n = [nspp]
  DATA_ARRAY( obs_catch );        // Observed fishery catch-at-age or catch-at-length; n = [nspp, fsh_age_bins, nyrs_fsh_comp]

  // May need to calculate in house
  DATA_ARRAY( fsh_age_obs);      // Observed fishery age comp; n = [nspp, fsh_age_bins, nyrs_fsh_comp]
  // matrix<Type> tc_obs(nspp, nyrs); // Observed total catch (n); n = [nspp, nyrs] NOTE: This may not be necessary if loading data from tmp

  // -- 2.3.2 BT Survey Components
  DATA_IVECTOR( nyrs_srv_biom );  // Number of years of survey biomass data; n = [nspp]
  DATA_IMATRIX( yrs_srv_biom );   // Years of survey biomass data; n = [nspp, nyrs_srv_biom]
  DATA_MATRIX( srv_bio );         // Observed BT survey biomass (kg); n = [nspp, nyrs]
  DATA_MATRIX( srv_biom_se );     // Observed annual biomass error (SE); n = [nspp, nyrs_srv_biom]

  DATA_IVECTOR( nyrs_srv_age);   // Number of years of survey age/length composition; n = [nspp]
  DATA_IMATRIX( yrs_srv_age);    // Years for the survey age/length composition data; n = [nspp, nyrs_srv_age]
  DATA_IVECTOR( srv_age_type );  // Type of compisition (1 = age; 2 = length); n = [nspp]
  DATA_IVECTOR( srv_age_bins);   // Number of size binds for the age/length comps; n = [nspp]
  DATA_MATRIX( srv_age_n);       // Sample size for the multinomial
  DATA_ARRAY( srv_age_obs);      // Observed BT age comp; n = [nspp, nages, nyrs]
  DATA_MATRIX( srv_age_sizes );  // Observed size composition
  DATA_ARRAY( age_trans_matrix); // observed sp_age/size compositions; n = [nspp, nages, srv_age_bins]

  // -- 2.3.3 EIT Survey components
  DATA_INTEGER( n_eit);          // Number of years with EIT data; n = [1]
  DATA_IVECTOR( yrs_eit);        // Years for available EIT data; n = [n_eit]
  DATA_VECTOR( eit_age_n);       // Number  of  EIT Hauls;
  DATA_MATRIX( obs_eit_age);     // Observed EIT age comp; n = [1, nages, nyrs] NOTE: may need to change this for future
  DATA_VECTOR( obs_eit);         // Observed EIT survey biomass (kg); n = [1, nyrs]
  DATA_MATRIX( eit_sel);         // Observed EIT survey selectivity; n = [eit_age, nyrs_eit_sel]
  // -- 2.3.4 Other
  DATA_VECTOR( TempC);           // Bottom temperature (degrees C); n = [1, nyrs] # NOTE: Need to figure out how to make it flexible for alternative environmental predictors
  DATA_ARRAY( Diet_Mat);         // Annual gravimetric proportion of prey in predator stomach; n = [n_pred, n_age_pred, nspp, nages, nyrs]
  DATA_INTEGER( other_food);     // Biomass of other prey (kg); n = [nyrs, n_pred] # QUESTION: Is this year specific?

  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( phi_p_bp);       // Annual relative foraging rate (d yr^-1)
  DATA_VECTOR( aLW);            //  Intercept of the allometric maximum consumption function (g g^-1 yr^-1); n = [1, nspp]
  DATA_VECTOR( bLW);            //  Allometric slope of maximum consumption; n = [1, nspp]
  DATA_VECTOR( Tcm);            //  Consumption maximum physiological temperature (degree C); n = [1, n_pred]
  DATA_VECTOR( Tco);            //  Consumption optimum physiological temperature (degree C); n = [1, n_pred]
  DATA_VECTOR( Qc);             //  Max consumption parameter; n = [1, n_pred]

  // 2.6. DERIVED QUANTITIES # Calculate these in the model
  DATA_MATRIX(d);               // VBGF allometric slope of consumption (d); n = [nspp, nyrs]
  DATA_MATRIX(Winf);            // VBGF max asymptoptic weight; n = [nspp, nyrs]
  DATA_MATRIX(Prop_Mat);        // Proportion of mature females at age; [nspp, nages]

  // 2.7. FIXED PARAMETERS
  // -- 2.7.1. von Bertalannfy growth function (VBGF)
  PARAMETER_VECTOR(d0_vbgf);      // VBGF intercept for d parameter; n = [1, nspp]
  PARAMETER_MATRIX(d_dev_vbgf);   // Annual deviation for VBGF d parameter; n = [nspp, nyrs] # NOTE: Need to figure out how to best vectorize this
  PARAMETER_VECTOR(Beta_d_vbgf);  // Temperature covariate for VBGF d parameter; n = [1, nspp]
  PARAMETER_VECTOR(K_vbgf);       // VBGF energy loss constant (kg kg^-1 yr^-1); n[1, nspp]
  PARAMETER_VECTOR(H_vbgf);       // VBGF assimilation constant (kg kg^-1 yr^-1); n[1, nspp]
  PARAMETER_VECTOR(t0_vbgf);      // VBGF age at Weight 0 (yr); n[1, nspp]

  // -- 2.7.2. Others
  PARAMETER_MATRIX(M1);         // Residual natural mortality; n = [nspp, nages]
  PARAMETER_VECTOR(M_f);        // Female natural mortality; n = [1, nspp]
  PARAMETER_VECTOR(M_m);        // Male natural mortality; n = [1, nspp]
  PARAMETER_MATRIX(Prop_age_f); // Proportion-at-age of females of population; n = [nspp, nages]
  PARAMETER_MATRIX(Prop_mat);   // Age-specific maturity; n = [nspp, nages]
  Type sigma_catch = 0.05;      // SD of catch

  // ------------------------------------------------------------------------- //
  // 3. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //
  // 3.1. PARAMETERS (assign parameters to objects)
  // -- 3.1.1 Recruitment parameters
  PARAMETER_VECTOR(ln_mn_rec);      // Mean recruitment; n = [1, nspp]
  PARAMETER_MATRIX(rec_dev);        // Annual recruitment deviation; n = [nspp, nyrs]
  // PARAMETER(sigma_rec);          // Standard deviation of recruitment variation # NOTE: Have this estimated if using random effects.
  // -- 3.1.2. Abundance parameters
  PARAMETER_ARRAY(init_dev);             // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this
  // -- 3.1.3. fishing mortality parameters
  PARAMETER_VECTOR(ln_mean_F);      // Log mean fishing mortality; n = [1, nspp]
  PARAMETER_MATRIX(F_dev);          // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change
  // PARAMETER_VECTOR(sigma_F); // SD of fishing mortality deviations; n = [1, nspp] # NOTE: Have this estimated if using random effects.
  // -- 3.1.4. Selectivity parameters
  PARAMETER_ARRAY(srv_sel_coff);    // Survey selectivity parameters; n = [nspp,nselages,sel_phase]
  PARAMETER_MATRIX(n_fish_sel);     // Fishery age selectivity coef; n = [nspp, 8]
  PARAMETER_VECTOR(b_surv_sel);     // Survey age selectivity slope ; n = [1, nspp]
  PARAMETER_VECTOR(a_surv_sel);     // Survey age selectivity limit ; n = [1, nspp]
  PARAMETER(log_eit_q);             // EIT Catchability
  PARAMETER_VECTOR(log_srv_q);      // BT Survey catchability; n = [1, nspp]

  // 3.2. DERIVED QUANTITIES
  int max_age = imax(nages);    // Integer of maximum nages to make the arrays.
  // -- 3.2.1. Fishery observations
  matrix<Type>  tc_biom_est(nspp, nyrs);              // Estimated total yield (kg); n = [nspp, nyrs]
  array<Type>   obs_catch_hat(nspp, max_age, nyrs);   // Estimated catch-at-age (n); n = [nspp, nages, nyrs]
  array<Type>   fsh_age_hat(nspp, max_age, nyrs);     // Estimated fishery age comp; n = [nspp, nages, nyrs]
  matrix<Type>  tc_hat(nspp, nyrs);                   // Estimated total catch (n); n = [nspp, nyrs]
  array<Type>   F(nspp, max_age, nyrs);               // Estimated fishing mortality; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_sel(nspp, max_age);               // Log estimated fishing selectivity
  // -- 3.2.2. BT Survey components
  vector<Type>  srv_bio_hat(nspp, nyrs);              // Estimated BT survey biomass (kg); n = [nspp, nyrs]
  array<Type>   srv_age_hat(nspp, max_age, nyrs);     // Estimated BT age comp; n = [nspp, nages, nyrs]
  matrix<Type>  srv_hat(nspp, nyrs);                  // Estimated BT survey total abundance (n); n = [nspp, nyrs]
  matrix<Type>  srv_sel(nspp, max_age);               // Estimated survey selectivity at age; n = [nspp, nyrs]
  // -- 3.2.3. EIT Survey Components
  array<Type>   Weight_at_Age(nspp, max_age, nyrs);   // Estimated weight-at-age; n = [nspp, nages, nyrs]
  matrix<Type>  eit_age_hat(12, n_eit);               // Estimated EIT age comp; n = [12 ages, nyrs]
  vector<Type>  eit_hat(n_eit);                       // Estimated EIT survey biomass (kg); n = [nyrs]

  // -- 3.3. ESTIMATED POPULATION PARAMETERS
  matrix<Type>  R(nspp, nyrs);                        // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type>   NByage(nspp, max_age, nyrs);          // Numbers at age; n = [nspp, nages, nyrs]
  array<Type>   biomassByage(nspp, max_age, nyrs);    // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type>  biomass(nspp, nyrs);                  // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type>  biomassSSB(nspp, nyrs);               // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type>   biomassSSBByage(nspp, max_age, nyrs); // Spawning biomass at age (kg); n = [nspp, nages, nyrs]

  // -- 3.4. Parameter transformations
  Type eit_q = exp(log_eit_q);                        // EIT Catchability

  array<Type>   Zed(nspp, max_age, nyrs);             // Total mortality at age; n = [nspp, nages, nyrs]

  // ------------------------------------------------------------------------- //
  // 4. MODEL OBJECTS                                                          //
  // ------------------------------------------------------------------------- //
  // 4.1. LOOPING INDICES -- k = observation, i = species/prey, j = age/prey age (yr), y = year, p = predator, a = predator age (yr)
  int k, i, j, y, p, a;

  // ------------------------------------------------------------------------- //
  // 5. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0
  //
  // 5.1. ESTIMATE RECRUITMENT T1.1
  for(i=0; i<nspp; y++){
    for(y=0; y<nyrs; y++){
      R(i,y) = exp(ln_mn_rec(i) + rec_dev(i,y));

      for(j=0; j<nages(i); j++){
        F(i,j,y) = exp(fsh_sel(i,j)) * exp(ln_mean_F(i) + F_dev(i,y));
      }
    }
  }

  // Debug it
  if(debug == 1){
    error("Unknown fleet code");
    return 0;
  }

  // 5.2. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      // -- 5.2.1. Plus group where y = 1 and 1 < j <= Ai
      if(j > 0 & j <= nages(i)){
        NByage(i,j,0) = ln_mn_rec(i) * exp(-(j+1)*M1(i,j)) * init_dev(i,j);
      }
      // -- 5.2.2. Where y = 1 and j > Ai.
      if(j>nages(i)){
        NByage(i,j,0) = ln_mn_rec(i) * exp(-(j+1)*M1(i,nages(i))) * init_dev(i,nages(i))/ (1-exp(-(j+1)*M1(i,nages(i)))); // NOTE: This solves for the geometric series
      }
    }
  }
  //
  // 5.3. ESTIMATE NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
  biomass.setZero();  // Initialize Biomass
  biomassSSB.setZero(); // Initialize SSB
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        // -- 5.3.1.  Where 1 <= j < Ai
        if(j >= 0 & j < (nages(i)-1)){
          NByage(i,j+1,y+1) = NByage(i,j,y) * exp(-Zed(i,j,y));
          // -- 5.3.2. Plus group where j > Ai. NOTE: This is not the same as T1.3 because I used j = A_i rather than j > A_i.
        }
        if(j == (nages(i)-1)){ // # NOTE: May need to increase j loop to "nages(i) + 1" .
          NByage(i,nages(i),y+1) = NByage(i,nages(i)-1,y) * exp(-Zed(i,nages(i)-1,y)) + NByage(i,nages(i),y) * exp(-Zed(i,nages(i),y));
        }
        biomassByage(i, j, y) = NByage(i, j, y) * Weight_at_Age(i, j, y); // 5.5.
        biomassSSBByage(i, j, y) = biomassByage(i, j, y) * Prop_Mat(i, j); // 5.6.

        // -- 5.3.3. Estimate Biomass and SSB
        biomass(i, y) += biomassByage(i, j, y);
        biomassSSB(i, y) += biomassSSBByage(i, j, y);
      }
    }
  }

  // ------------------------------------------------------------------------- //
  // 6. FISHERY COMPONENTS EQUATIONS                                           //
  // ------------------------------------------------------------------------- //
  // 6.1. ESTIMATE FISHING MORTALITY
  for(i=0; i<nspp; y++){
    for(y=0; y<nyrs; y++){
      for(j=0; j<nages(i); j++){
        F(i,j,y) = exp(fsh_sel(i,j)) * exp(ln_mean_F(i) + F_dev(i,y));
      }
    }
  }

  // 6.1. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
  tc_hat.setZero();
  tc_biom_est.setZero(); // Initialize tc_biom_est
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        obs_catch_hat(i, j, y) = F(i, j, y)/Zed(i, j, y) * (1 - exp(-Zed(i, j, y))) * NByage(i,j,y); // 5.4.
        tc_hat(i, y) += obs_catch_hat(i, j, y); // Estimate catch in numbers
        tc_biom_est(i, y) += obs_catch_hat(i, j, y) * Weight_at_Age(i, j, y); // 5.5.
      }
    }
  }

  // 5.5. ESTIMATE FISHERY AGE COMPOSITION
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        // -- 5.5.1 Estimate age composition of the fishery
        if(fsh_age_type(sp)==1){
          fsh_age_hat(i, j, y) = obs_catch_hat(i, j, y) / tc_hat(i, y);
        }
        if(fsh_age_type(sp)!=1){
          // ADD AGE TRANSITION MATRIX TO INPUT
          fsh_age_hat(i, j, y)  = obs_catch_hat(i, j, y)*age_trans_matrix.col().col(i) / tc_hat(sp,yr); // NOTE: double check indexing here of the age_trans_matrix
        }
      }
    }
  }

  // 5.8. Estimate total mortality at age
  // 5.9. Estimate fishing mortality at age
  // 5.17. Estimate fishery selectivity
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        Zed(i,j,y) = M1(i,j) + F(i, j, y);
      }
    }
  }

  // ------------------------------------------------------------------------- //
  // 7. SURVEY COMPONENTS EQUATIONS                                            //
  // ------------------------------------------------------------------------- //
  // 576 EIT Components
  // -- 7.6.1 EIT Survey Biomass
  eit_age_hat.setZero();
  eit_hat.setZero();
  for(j=0; j < nages(0); j++){
    for(y=0; y < n_eit; y++){
      eit_age_hat(j,y) = NByage(0,j,y) * exp(-Zed(0,j,y)) * eit_sel(j, y) * eit_q; // Remove the mid-year trawl?
      eit_hat(y) += eit_age_hat(j,y) * Weight_at_Age(0, j, y);  //
    }
  }
  // -- 7.6.1 EIT Survey Age Composition
  for(j=0; j < nages(0); j++){
    for (y=0; y < n_eit; y++){
      eit_age_hat(j, y) = eit_age_hat(j, y) / eit_age_hat.col(y).sum();
    }
  }

  // 7.7 BT Components
  // -- 7.7.1 BT Survey Biomass
  srv_age_hat.setZero();
  srv_hat.setZero();
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < n_eit; y++){
        srv_age_hat(i, j, y) = NByage(i,j,y) * exp(-0.5 * Zed(i,j,y)) * srv_sel(i, j) * srv_q(i);
        srv_hat(i, y) += srv_age_hat(i, j, y)   //
          srv_bio_hat(i, y) += srv_age_hat(i, j, y) * Weight_at_Age(i, j, y);  //
      }
    }
  }

  // -- 7.7.2 BT Survey Age Composition
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for (y=0; y < n_eit; y++){
        if(srv_age_type(i)==1){
          srv_age_hat(i, j, y) = srv_age_hat(i, j, y) / srv_hat(i, y)
        }
        if(srv_age_type(i)!=1){
          srv_age_hat(i, j, y) = (srv_age_hat(i, j, y) * *age_trans_matrix.col().col(i)) / srv_hat(i, y) // # NOTE: Double check the matrix algebra here
        }
      }
    }
  }

  // ------------------------------------------------------------------------- //
  // 8. LIKELIHOOD EQUATIONS                                                   //
  // ------------------------------------------------------------------------- //
  // 8.1. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(10,nspp); // matrix of negative log-likelihood components
  // -- Data components
  // Slot 0 -- BT survey biomass -- NFMS annual BT survey
  // Slot 1 -- BT survey age composition -- NFMS annual BT survey
  // Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
  // Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
  // Slot 4 -- Total catch -- Fishery observer data
  // Slot 5 -- Fishery age composition -- Fishery observer data
  // -- Likelihood penalties
  // Slot 6 -- Fishery selectivity
  // -- Priors
  // Slot 7 -- Tau -- Annual recruitment deviation
  // Slot 8 -- init_dev -- Initial abundance-at-age
  // Slot 9 -- Epsilon -- Annual fishing mortality deviation
  jnll_comp.setZero();
  Type jnll = 0;

  // 8.2. Fit the objective function
  // Slot 0 -- BT survey biomass -- NFMS annual BT survey
  for(i=0; i < nspp; i++){
    for (y=0; y < n_eit; y++){
      jnll_comp(0, i) += 12.5 * pow(log(srv_bio(i, y)) - log(srv_bio_hat(i, y)  + 1.e-04)), 2);
    }
  }

  // Slot 1 -- BT survey age composition -- NFMS annual BT survey
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for (y=0; y < nyrs_srv_biom(i); y++){
        jnll_comp(1, i) += eit_age_n(y) *  (srv_age_obs(i, j, y) + MNConst) * log(srv_age_hat(i, j, y) + MNConst);
      }
    }
  }

  // Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
  for (y=0; y < n_eit; y++){
    jnll_comp(2, 0) += 12.5 * pow(log(obs_eit(y)) - log(eit_hat(y) + 1.e-04)), 2);
  }

  // Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
  for(j=0; j < nages(0); j++){
    for (y=0; y < n_eit; y++){
      jnll_comp(3, 0) += eit_age_n(y) *  (obs_eit_age(j, y) + MNConst) * log(eit_age_hat(j, y) + MNConst);
    }
  }

  // Slot 4 -- Total catch -- Fishery observer data
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs; y++){
      jnll_comp(4,i) = pow((log(tc_biom_est(i,y)) - log(tc_biom_obs(i, y))), 2) / (2 * pow(sigma_catch, 2)); // T.4.5
    }
  }

  // Slot 5 -- Fishery age composition -- Fishery observer data
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for (y=0; y < n_eit; y++){
        jnll_comp(5, i) -= tau * (fsh_age_obs(i, j, y) + MNConst)* log(fsh_age_hat(i, j, y) + MNConst);
      }
    }
  }

  // Slot 6 -- Fishery selectivity

  // PRIORS: PUT RANDOM EFFECTS SWITCH HERE
  for(i=0; i < nspp; i++){
    // Slot 8 -- init_dev -- Initial abundance-at-age
    for(j=0; j < nages(i); j++){
      jnll_comp(8,i) += pow(init_dev(i,j), 2);
    }

    // Slot 7 -- Tau -- Annual recruitment deviation
    // Slot 9 -- Epsilon -- Annual fishing mortality deviation
    for (y=0; y < nyrs; y++){
      jnll_comp(7,i) += pow(rec_dev(i,y), 2);     // Recruitment deviation using penalized likelihood.
      jnll_comp(9,i) += pow(F_dev(i,y), 2);       // Fishing mortality deviation using penalized likelihood.
    }
  }

  // ------------------------------------------------------------------------- //
  // 9. REPORT SECTION                                                         //
  // ------------------------------------------------------------------------- //

  jnll = jnll_comp.sum();
  REPORT(jnll);
  return jnll;

  // ------------------------------------------------------------------------- //
  // END MODEL                                                                 //
  // ------------------------------------------------------------------------- //
}
