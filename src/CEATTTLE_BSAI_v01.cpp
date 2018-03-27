#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
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
  // 1. MODEL CONFIGURATION                                                    //
  // ------------------------------------------------------------------------- //
  // 1.1. CONFIGURE MODEL (this section sets up the switches)

  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                    //
  // ------------------------------------------------------------------------- //
  // 2.1. FIXED VALUES
  int tau = 200; // Fishery age composition sample size
  Type MNConst = 0.001; // Constant additive for logistic functions


  // 2.2. DIMENSIONS OF DATA
  // -- 2.2.1. Temporal dimensions
  DATA_INTEGER(nyrs); //Number of simulation years
  DATA_INTEGER(styr); // Start year
  
  // -- 2.2.2. Number of species
  DATA_INTEGER(nspp); // Number of species (prey)
  DATA_IVECTOR(nages); // Number of species (prey) ages
  DATA_INTEGER(n_pred); // Number of predator species
  DATA_IVECTOR(n_age_pred); // Number of predator ages
  
  // 2.3. DATA INPUTS (i.e. assign data to objects)
  DATA_MATRIX(tc_biom_obs); // Observed total yield (kg); n = [nspp, nyrs]
  DATA_ARRAY(obs_catch); // Observed fishery catch-at-age; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  DATA_ARRAY(fsh_age_obs); // Observed fishery age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  DATA_ARRAY(srv_age_obs); // Observed BT age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  DATA_MATRIX(obs_eit_age); // Observed EIT age comp; n = [1, nages, nyrs] 
  DATA_MATRIX(srv_biom); // Observed BT survey biomass (kg); n = [nspp, nyrs] 
  DATA_VECTOR(obs_eit); // Observed EIT survey biomass (kg); n = [1, nyrs] 
  DATA_VECTOR(TempC); // Bottom temperature (degrees C); n = [1, nyrs] # NOTE: Need to figure out how to make it flexible for alternative environmental predictors
  DATA_ARRAY(Diet_Mat); // Annual gravimetric proportion of prey in predator stomach; n = [n_pred, n_age_pred, nspp, nages, nyrs]
  DATA_INTEGER(other_food); // Biomass of other prey (kg); n = [nyrs, n_pred] # QUESTION: Is this year specific?
  
  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR(phi_p_bp); // Annual relative foraging rate (d yr^-1)
  DATA_VECTOR(aLW); //  Intercept of the allometric maximum consumption function (g g^-1 yr^-1); n = [1, nspp]
  DATA_VECTOR(bLW); //  Allometric slope of maximum consumption; n = [1, nspp]
  DATA_VECTOR(Tcm); //  Consumption maximum physiological temperature (degree C); n = [1, n_pred]
  DATA_VECTOR(Tco); //  Consumption optimum physiological temperature (degree C); n = [1, n_pred]
  DATA_VECTOR(Qc); //  Max consumption parameter; n = [1, n_pred]
  
  // 2.6. DERIVED QUANTITIES # Calculate these in the model
  DATA_MATRIX(d); // VBGF allometric slope of consumption (d); n = [nspp, nyrs]
  DATA_MATRIX(Winf); // VBGF max asymptoptic weight; n = [nspp, nyrs]
  DATA_MATRIX(Prop_Mat); // Proportion of mature females at age; [nspp, nages]
  
  // 2.7. FIXED PARAMETERS
  // -- 2.7.1. von Bertalannfy growth function (VBGF)
  PARAMETER_VECTOR(d0_vbgf); // VBGF intercept for d parameter; n = [1, nspp]
  PARAMETER_MATRIX(d_dev_vbgf); // Annual deviation for VBGF d parameter; n = [nspp, nyrs] # NOTE: Need to figure out how to best vectorize this
  PARAMETER_VECTOR(Beta_d_vbgf); // Temperature covariate for VBGF d parameter; n = [1, nspp]
  PARAMETER_VECTOR(K_vbgf); // VBGF energy loss constant (kg kg^-1 yr^-1); n[1, nspp]
  PARAMETER_VECTOR(H_vbgf); // VBGF assimilation constant (kg kg^-1 yr^-1); n[1, nspp]
  PARAMETER_VECTOR(t0_vbgf); // VBGF age at Weight 0 (yr); n[1, nspp]
  
  // -- 2.7.2. Others
  PARAMETER_MATRIX(M1); // Residual natural mortality; n = [nspp, nages]
  PARAMETER_VECTOR(S_EIT); // EIT survey selectivity # QUESTION: What size is this?
  PARAMETER_VECTOR(M_f); // Female natural mortality; n = [1, nspp]
  PARAMETER_VECTOR(M_m); // Male natural mortality; n = [1, nspp]
  PARAMETER_MATRIX(Prop_age_f); // Proportion-at-age of females of population; n = [nspp, nages]
  PARAMETER_MATRIX(Prop_mat); // Age-specific maturity; n = [nspp, nages]
  int sigma_catch = 0.05; // SD of catch
  
  // ------------------------------------------------------------------------- //
  // 3. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //
  // 3.1. PARAMETERS (assign parameters to objects)
  // -- 3.1.1 Recruitment parameters
  PARAMETER_VECTOR(ln_mn_rec);  // Mean recruitment; n = [1, nspp]
  PARAMETER_MATRIX(rec_dev); // Annual recruitment deviation; n = [nspp, nyrs]
  // PARAMETER(sigma_rec); // Standard deviation of recruitment variation # NOTE: Have this estimated if using random effects.
  // -- 3.1.2. Abundance parameters
  PARAMETER_ARRAY(N_0); // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this
  // -- 3.1.3. fishing mortality parameters
  PARAMETER_VECTOR(F_0); // Mean fishing mortality; n = [1, nspp]
  PARAMETER_ARRAY(F_dev); // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change
  // PARAMETER_VECTOR(sigma_F); // SD of fishing mortality deviations; n = [1, nspp] # NOTE: Have this estimated if using random effects.
  // -- 3.1.4. Selectivity parameters
  PARAMETER_MATRIX(n_fish_sel); // Fishery age selectivity coef; n = [nspp, 8]
  PARAMETER_VECTOR(b_surv_sel); // Survey age selectivity slope ; n = [1, nspp]
  PARAMETER_VECTOR(a_surv_sel); // Survey age selectivity limit ; n = [1, nspp]
  
  // 3.2. DERIVED QUANTITIES
  // -- 3.2.1. Estimated observations
  matrix<Type> tc_biom_est(nspp, nyrs); // Estimated total yield (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  array<Type> obs_catch_hat(nspp, nages, nyrs); // Estimated catch-at-age; n = [nspp, nages, nyrs]
  matrix<Type> tc_est(nspp, nyrs); // Estimated catch (n); n = [nspp, nyrs]
  array<Type> Weight_at_Age(nspp, nages, nyrs); // Estimated weight-at-age; n = [nspp, nages, nyrs]
  array<Type> fsh_age_hat(nspp, nages, nyrs); // Estimated fishery age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  array<Type> Est_BT_Age_Comp(nspp, nages, nyrs); // Estimated BT age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  array<Type> Est_EIT_Age_Comp(nspp, nages, nyrs); // Estimated EIT age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  vector<Type> Est_BT_IOA(nspp, nyrs); // Estimated BT survey biomass (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  vector<Type> Est_EIT_IOA(nspp, nyrs); // Estimated BT survey biomass (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.

  // -- 3.2.2. Estimated Population Parameters 
  matrix<Type> R(nspp, nyrs); // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type> NByage(nspp, nages, nyrs); // Numbers at age; n = [nspp, nages, nyrs]
  array<Type> biomassByage(nspp, nages, nyrs); // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type> biomass(nspp, nyrs); // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type> biomassSSB(nspp, nyrs); // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type> biomassSSBByage(nspp, nages, nyrs); // Spawning biomass at age (kg); n = [nspp, nages, nyrs]

  array<Type> Zed(nspp, nages, nyrs); // Total mortality at age; n = [nspp, nages, nyrs]
  array<Type> F_Tot(nspp, nages, nyrs); // Fishing mortality at age; n = [nspp, nages, nyrs]
  matrix<Type> R(nspp, nyrs); // Recruitment; n = [nspp, nyrs]

  // ------------------------------------------------------------------------- //
  // 4. MODEL OBJECTS AND LIKELIHOOD                                                      //
  // ------------------------------------------------------------------------- //
  // 4.1. LOOPING INDICES -- k = observation, i = species/prey, j = age/prey age (yr), y = year, p = predator, a = predator age (yr)
  int k, i, j, y, p, a;
  
  // 4.2. OBJECTIVE FUNCTION
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
  // Slot 8 -- N_0 -- Initial abundance-at-age
  // Slot 9 -- Epsilon -- Annual fishing mortality deviation
  jnll_comp.setZero();
  Type jnll = 0;
    
  
  // ------------------------------------------------------------------------- //
  // 5. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0
  //
  // 5.1. ESTIMATE RECRUITMENT: T1.1
  for(i=0; i<nspp; y++){
    for(y=0; y<n_yrs; y++){
      R(i,y) = exp(ln_mn_rec(i) + rec_dev(i,y));
      
      // Likelihood 
      jnll_comp(7,i) += rec_dev(i,y); // Recruitment deviation using penalized likelihood. QUESTION: Do we want recruitment deviations for each species to have seperate sigmas?
    }
  }
  
  //
  // 5.2. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      // -- 5.2.1. Plus group where y = 1 and 1 < j <= Ai
      if(j > 0 & j <= nages(i)){
        NByage(i,j,0) = ln_mn_rec(i) * exp(-(j+1)*M1(i,j)) * N_0(i,j);
      }
      // -- 5.2.2. Where y = 1 and j > Ai.
      if(j>nages(i)){
        NByage(i,j,0) = ln_mn_rec(i) * exp(-(j+1)*M1(i,nages(i))) * N_0(i,nages(i))/ (1-exp(-(j+1)*M1(i,nages(i)))); // NOTE: This solves for the geometric series
      }
    }
  }
  //
  // 5.3. ESTIMATE NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
  Biomass.setZero();  // Initialize Biomass 
  BiomassSSB.setZero(); // Initialize SSB 
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
        Biomass(i, y) += biomassByage(i, j, y);
        BiomassSSB(i, y) += BiomassSSB(i, j, y);
      }
    }
  }

  // 5.4. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
  tc_est.setZero();
  tc_biom_est.setZero(); // Initialize tc_biom_est
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        obs_catch_hat(i, j, y) = F_Tot(i, j, y)/Zed(i, j, y) * (1 - exp(-Zed(i, j, y))) * NByage(i,j,y); // 5.4.
        tc_est(i, y) += obs_catch_hat(i, j, y); // Estimate catch in numbers
        tc_biom_est(i, y) += obs_catch_hat(i, j, y) * Weight_at_Age(i, j, y); // 5.5.
        
        // Total Catch Likelihood 
        jnll_comp(4,i) = log(tc_biom_est(i,y) - log(tc_biom_obs(i, y))^2 / (2 * sigma_catch ^ 2); // T.4.5
      }
    }
  }

  // 5.5. ESTIMATE FISHERY AGE COMPOSITION
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
      // -- 5.5.1 Estimate age composition of the fishery
      fsh_age_hat(i, j, y) = obs_catch_hat(i, j, y) / tc_est(i, y);

      jnll_comp(5, i) -= tau * (fsh_age_obs(i, j, y) + MNConst)* log(fsh_age_hat(i, j, y) + MNConst);
    }
  }
  
  // 5.8. Estimate total mortality at age
  // 5.9. Estimate fishing mortality at age
  // 5.17. Estimate fishery selectivity
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      M1(i,j) = 
      
      for(y=0; y < nyrs; y++){
        F_Tot(i, j, y) = F_0(i) * exp(F_dev(i,y));
        Zed(i,j,y) = M1(i,j) + F_Tot(i, j, y);
        
        // Likelihood 
        jnll_comp(9, i) += F_dev(i,y); // Penalized likelihood for Fishing mortality deviation.
      }
    }
  }

  jnll = jnll_comp.sum()
  // ------------------------------------------------------------------------- //
  // END MODEL                                                                 //
  // ------------------------------------------------------------------------- //
}
