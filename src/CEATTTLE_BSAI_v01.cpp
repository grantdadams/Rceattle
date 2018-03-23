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
  
  // 1.2. DIMENSIONS OF DATA
  // -- 1.2.1. Temporal dimensions
  DATA_INTEGER(nyrs); //Number of simulation years
  DATA_INTEGER(styr); // Start year
  
  // -- 1.2.2. Number of species
  DATA_INTEGER(nspp); // Number of species (prey)
  DATA_IVECTOR(nages); // Number of species (prey) ages
  DATA_INTEGER(n_pred); // Number of predator species
  DATA_IVECTOR(n_age_pred); // Number of predator ages
  
  // 1.3. DATA VECTORS (i.e. assign data to objects)
  DATA_MATRIX(Tot_Yield); // Observed total yield (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  DATA_ARRAY(F_Age_Comp); // Observed fishery age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  DATA_ARRAY(BT_Age_Comp); // Observed BT age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  DATA_ARRAY(EIT_Age_Comp); // Observed EIT age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  DATA_VECTOR(BT_IOA); // Observed BT survey biomass (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  DATA_VECTOR(EIT_IOA); // Observed EIT survey biomass (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  DATA_VECTOR(BOT_TEMP); // Bottom temperature (degrees C); n = [1, nyrs] # NOTE: Need to figure out how to make it flexible for alternative environmental predictors
  DATA_ARRAY(Diet_Mat); // Annual gravimetric proportion of prey in predator stomach; n = [n_pred, n_age_pred, nspp, nages, nyrs]
  DATA_MATRIX(B_Prey); // Biomass of other prey (kg); n = [nyrs, n_pred] # QUESTION: Is this year specific?
  
  // 1.4. INPUT PARAMETERS
  // -- 1.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR(phi_p_bp); // Annual relative foraging rate (d yr^-1)
  DATA_VECTOR(aLW); //  Intercept of the allometric maximum consumption function (g g^-1 yr^-1); n = [1, nspp]
  DATA_VECTOR(bLW); //  Allometric slope of maximum consumption; n = [1, nspp]
  DATA_VECTOR(Tcm); //  Consumption maximum physiological temperature (degree C); n = [1, n_pred]
  DATA_VECTOR(Tco); //  Consumption optimum physiological temperature (degree C); n = [1, n_pred]
  DATA_VECTOR(Qc); //  Max consumption parameter; n = [1, n_pred]
  
  // 1.5. ESTIMATED PARAMETERS (assign parameters to objects)
  // -- 1.5.1 Recruitment parameters
  PARAMETER_VECTOR(R_0);  // Mean recruitment; n = [1, nspp]
  PARAMETER_MATRIX(R_dev); // Annual recruitment deviation; n = [nspp, nyrs]
  PARAMETER(sigma_rec); // Standard deviation of recruitment variation
  // -- 1.5.2. Abundance parameters
  PARAMETER_ARRAY(N_0); // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this
  // -- 1.5.3. fishing mortality parameters
  PARAMETER_VECTOR(F_0); // Mean fishing mortality; n = [1, nspp]
  PARAMETER_ARRAY(F_dev); // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change
  PARAMETER_VECTOR(sigma_F); // SD of fishing mortality deviations; n = [1, nspp]
  // -- 1.5.4. Selectivity parameters
  PARAMETER_MATRIX(n_fish_sel); // Fishery age selectivity coef; n = [nspp, 8]
  PARAMETER_VECTOR(b_surv_sel); // Survey age selectivity slope ; n = [1, nspp]
  PARAMETER_VECTOR(a_surv_sel); // Survey age selectivity limit ; n = [1, nspp]
  
  // 1.6. DERIVED QUANTITIES # Calculate these in the model
  DATA_MATRIX(d); // VBGF allometric slope of consumption (d); n = [nspp, nyrs]
  DATA_MATRIX(Winf); // VBGF max asymptoptic weight; n = [nspp, nyrs]
  DATA_MATRIX(Prop_Mat); // Proportion of mature females at age; [nspp, nages]
  
  // 1.7. FIXED PARAMETERS
  // -- 1.7.1. von Bertalannfy growth function (VBGF)
  PARAMETER_VECTOR(d0_vbgf); // VBGF intercept for d parameter; n = [1, nspp]
  PARAMETER_MATRIX(d_dev_vbgf); // Annual deviation for VBGF d parameter; n = [nspp, nyrs] # NOTE: Need to figure out how to best vectorize this
  PARAMETER_VECTOR(Beta_d_vbgf); // Temperature covariate for VBGF d parameter; n = [1, nspp]
  PARAMETER_VECTOR(K_vbgf); // VBGF energy loss constant (kg kg^-1 yr^-1); n[1, nspp]
  PARAMETER_VECTOR(H_vbgf); // VBGF assimilation constant (kg kg^-1 yr^-1); n[1, nspp]
  PARAMETER_VECTOR(t0_vbgf); // VBGF age at Weight 0 (yr); n[1, nspp]
  
  // -- 1.7.2. Others
  PARAMETER_MATRIX(M1); // Residual natural mortality; n = [nspp, nages]
  PARAMETER_VECTOR(S_EIT); // EIT survey selectivity # QUESTION: What size is this?
  PARAMETER_VECTOR(M_f); // Female natural mortality; n = [1, nspp]
  PARAMETER_VECTOR(M_m); // Male natural mortality; n = [1, nspp]
  PARAMETER_MATRIX(Prop_age_f); // Proportion-at-age of females of population; n = [nspp, nages]
  PARAMETER_MATRIX(Prop_mat); // Age-specific maturity; n = [nspp, nages]
  int sigma_catch = 0.05; // SD of catch
  
  // ------------------------------------------------------------------------- //
  // 2. ESTIMATED OBJECT CONFIGURATION                                         //
  // ------------------------------------------------------------------------- //
  // 2.1. LOOPING INDICES -- k = observation, i = species/prey, j = age/prey age (yr), y = year, p = predator, a = predator age (yr)
  int k, i, j, y, p, a;
  
  // 2.2. OBJECTIVE FUNCTION
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
  Type jnll = 0; // Joint nll
  
  // 2.3. ESTIMATED QUANTITIES
  // -- 2.3.2. Estimated observations
  matrix<Type> Est_Tot_Yield(nspp, nyrs); // Estimated total yield (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  array<Type> Est_Catch_at_Age(nspp, nages, nyrs); // Estimated catch-at-age; n = [nspp, nages, nyrs]
  array<Type> Weight_at_Age(nspp, nages, nyrs); // Estimated weight-at-age; n = [nspp, nages, nyrs]
  array<Type> Est_B_at_Age(nspp, nages, nyrs); // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  array<Type> Est_F_Age_Comp(nspp, nages, nyrs); // Estimated fishery age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of fleets. Maybe look at ASAP code.
  array<Type> Est_BT_Age_Comp(nspp, nages, nyrs); // Estimated BT age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  array<Type> Est_EIT_Age_Comp(nspp, nages, nyrs); // Estimated EIT age comp; n = [nspp, nages, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  vector<Type> Est_BT_IOA(nspp, nyrs); // Estimated BT survey biomass (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  vector<Type> Est_EIT_IOA(nspp, nyrs); // Estimated BT survey biomass (kg); n = [nspp, nyrs] # NOTE: Need to figure out how to make the age comps flexible for varying number of surveys.
  // -- 2.3.2. Estimated Population Parameters 
  array<Type> N_at_age(nspp, nages, nyrs); // Numbers at age; n = [nspp, nages, nyrs]
  array<Type> SSB_at_Age(nspp, nages, nyrs); // Spawning biomass at age (kg); n = [nspp, nages, nyrs]
  array<Type> Zed(nspp, nages, nyrs); // Total mortality at age; n = [nspp, nages, nyrs]
  array<Type> F_Tot(nspp, nages, nyrs); // Fishing mortality at age; n = [nspp, nages, nyrs]
  
  
  // ------------------------------------------------------------------------- //
  // 3. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0
  //
  // 3.1. ESTIMATE RECRUITMENT: T1.1
  for(i=0; i<nspp; y++){
    for(y=0; y<n_yrs; y++){
      N_at_age(i,0,y) = R_0(i) * exp(R_dev(i,y));
      
      // Likelihood 
      jnll_comp(7,nspp) -= R_dev(i,y); // Recruitment deviation using penalized likelihood. QUESTION: Do we want recruitment deviations for each species to have seperate sigmas?
    }
  }
  
  //
  // 3.2. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      // -- 3.2.1. Plus group where y = 1 and 1 < j <= Ai
      if(j > 0 & j <= nages(i)){
        N_at_age(i,j,0) = R_0(i) * exp(-(j+1)*M1(i,j)) * N_0(i,j);
      }
      // -- 3.2.2. Where y = 1 and j > Ai.
      if(j>nages(i)){
        N_at_age(i,j,0) = R_0(i) * exp(-(j+1)*M1(i,nages(i))) * N_0(i,nages(i))/ (1-exp(-(j+1)*M1(i,nages(i)))); // NOTE: This solves for the geometric series
      }
    }
  }
  //
  // 3.3. ESTIMATE NUMBERS AT AGE
  // 3.6. ESTIMATE BIOMASS-AT-AGE (kg)
  // 3.7. ESTIMATE SSB-AT-AGE (kg)
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        // -- 3.3.1.  Where 1 <= j < Ai
        if(j >= 0 & j < (nages(i)-1)){
          N_at_age(i,j+1,y+1) = N_at_age(i,j,y) * exp(-Zed(i,j,y))
          // -- 3.3.2. Plus group where j > Ai. NOTE: This is not the same as T1.3 because I used j = A_i rather than j > A_i.
        }
        if(j == (nages(i)-1)){ // # NOTE: May need to increase j loop to "nages(i) + 1" .
          N_at_age(i,nages(i),y+1) = N_at_age(i,nages(i)-1,y) * exp(-Zed(i,nages(i)-1,y)) + N_at_age(i,nages(i),y) * exp(-Zed(i,nages(i),y))
        }
        Est_B_at_Age(i, j, y) = Est_Catch_at_Age(i, j, y) * Weight_at_Age(i, j, y); // 3.5.
        SSB_at_Age(i, j, y) = Est_B_at_Age(i, j, y) * Prop_Mat(i, j); // 3.6.
      }
    }
  }
  //
  // 3.4. ESTIMATE CATCH-AT-AGE
  // 3.5. ESTIMATE TOTAL YIELD (kg)
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        Est_Catch_at_Age(i, j, y) = F_Tot(i, j, y)/Zed(i, j, y) * (1 - exp(-Zed(i, j, y))) * N_at_age(i,j,y); // 3.4.
        Est_Tot_Yield(i, y) += Est_Catch_at_Age(i, j, y) * Weight_at_Age(i, j, y); // 3.5.
        
        // Likelihood 
        jnll_comp(4,nspp) -= log(Est_Tot_Yield(i,y) - log(Tot_Yield(i, y))^2 / (2 * sigma_catch ^ 2); // T.4.5
      }
    }
  }
  //
  // 3.8. Estimate total mortality at age
  // 3.9. Estimate fishing mortality at age
  // 3.17. Estimate fishery selectivity
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      M1(i,j) = 
      
      for(y=0; y < nyrs; y++){
        F_Tot(i, j, y) = F_0(i) * exp(F_dev(i,y));
        Zed(i,j,y) = M1(i,j) + F_Tot(i, j, y);
        
        // Likelihood 
        jnll_comp(9,nspp) -= F_dev(i,y); // Penalized likelihood for Fishing mortality deviation.
      }
    }
  }
  // ------------------------------------------------------------------------- //
  // END MODEL                                                                 //
  // ------------------------------------------------------------------------- //
}
