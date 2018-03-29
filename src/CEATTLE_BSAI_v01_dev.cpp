#include <TMB.hpp>

// Function for getting max of an Ivector
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
  int nselages = 8;


  // 2.2. DIMENSIONS OF DATA
  // -- 2.2.1. Temporal dimensions
  DATA_INTEGER(nyrs); //Number of simulation years
  DATA_INTEGER(styr); // Start year

  // -- 2.2.2. Number of species
  DATA_INTEGER(nspp); // Number of species (prey)
  DATA_IVECTOR(nages); // Number of species (prey) ages
  // DATA_INTEGER(n_pred); // Number of predator species
  // DATA_IVECTOR(n_age_pred); // Number of predator ages

  // 2.3. DATA INPUTS (i.e. assign data to objects)
  // -- 2.3.1 Fishery Components
  DATA_MATRIX(tc_biom_obs); // Observed total yield (kg); n = [nspp, nyrs]
  DATA_ARRAY(obs_catch); // Observed fishery catch-at-age; n = [nspp, nages, nyrs]
  DATA_ARRAY(fsh_age_obs); // Observed fishery age comp; n = [nspp, nages, nyrs]
  // -- 2.3.2 BT Survey Components
  DATA_ARRAY(srv_age_obs); // Observed BT age comp; n = [nspp, nages, nyrs]
  DATA_MATRIX(srv_biom); // Observed BT survey biomass (kg); n = [nspp, nyrs]
  // -- 2.3.3 EIT Survey components
  DATA_INTEGER(n_eit); // Number of years with EIT data
  DATA_IVECTOR(yrs_eit); // Years for available EIT data
  DATA_VECTOR(eit_age_n); // Number  of  EIT Hauls
  DATA_MATRIX(obs_eit_age); // Observed EIT age comp; n = [1, nages, nyrs]
  DATA_VECTOR(obs_eit); // Observed EIT survey biomass (kg); n = [1, nyrs]
  DATA_MATRIX(eit_sel); // Observed EIT survey selectivity; n = [eit_age, nyrs_eit_sel]
  // -- 2.3.4 Other
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
  PARAMETER_VECTOR(M_f); // Female natural mortality; n = [1, nspp]
  PARAMETER_VECTOR(M_m); // Male natural mortality; n = [1, nspp]
  PARAMETER_MATRIX(Prop_age_f); // Proportion-at-age of females of population; n = [nspp, nages]
  PARAMETER_MATRIX(Prop_mat); // Age-specific maturity; n = [nspp, nages]
  Type sigma_catch = 0.05; // SD of catch

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
  PARAMETER_VECTOR(ln_mean_F); // Log mean fishing mortality; n = [1, nspp]
  PARAMETER_MATRIX(F_dev); // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change
  // PARAMETER_VECTOR(sigma_F); // SD of fishing mortality deviations; n = [1, nspp] # NOTE: Have this estimated if using random effects.
  // -- 3.1.4. Selectivity parameters
  PARAMETER_MATRIX(n_fish_sel); // Fishery age selectivity coef; n = [nspp, 8]
  PARAMETER_VECTOR(b_surv_sel); // Survey age selectivity slope ; n = [1, nspp]
  PARAMETER_VECTOR(a_surv_sel); // Survey age selectivity limit ; n = [1, nspp]
  PARAMETER(log_eit_q); // EIT Catchability


  // 3.2. DERIVED QUANTITIES
  int max_age = imax(nages); // Integer of maximum nages to make the arrays.
  // -- 3.2.1. Fishery observations
  matrix<Type> tc_biom_est(nspp, nyrs); // Estimated total yield (kg); n = [nspp, nyrs]
  array<Type> obs_catch_hat(nspp, max_age, nyrs); // Estimated catch-at-age (n); n = [nspp, nages, nyrs]
  array<Type> fsh_age_hat(nspp, max_age, nyrs); // Estimated fishery age comp; n = [nspp, nages, nyrs]
  matrix<Type> tc_est(nspp, nyrs); // Estimated total catch (n); n = [nspp, nyrs]
  array<Type> F(nspp, max_age, nyrs); // Estimated fishing mortality; n = [nspp, nages, nyrs]
  matrix<Type> fsh_sel(nspp, max_age); // Estimated fishing selectivity
  // -- 3.2.2. BT Survey components
  vector<Type> Est_BT_IOA(nspp, nyrs); // Estimated BT survey biomass (kg); n = [nspp, nyrs]
  array<Type> Est_BT_Age_Comp(nspp, max_age, nyrs); // Estimated BT age comp; n = [nspp, nages, nyrs]
  // -- 3.2.3. EIT Survey Components
  array<Type> Weight_at_Age(nspp, max_age, nyrs); // Estimated weight-at-age; n = [nspp, nages, nyrs]
  matrix<Type> eit_age_hat(12, n_eit); // Estimated EIT age comp; n = [12 ages, nyrs]
  vector<Type> eit_hat(n_eit); // Estimated EIT survey biomass (kg); n = [nyrs]

  // -- 3.2.2. Estimated Population Parameters
  matrix<Type> R(nspp, nyrs); // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type> NByage(nspp, max_age, nyrs); // Numbers at age; n = [nspp, nages, nyrs]
  array<Type> biomassByage(nspp, max_age, nyrs); // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type> biomass(nspp, nyrs); // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type> biomassSSB(nspp, nyrs); // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type> biomassSSBByage(nspp, max_age, nyrs); // Spawning biomass at age (kg); n = [nspp, nages, nyrs]

  // -- 3.2.3. Parameter transformations
  Type eit_q = exp(log_eit_q); // EIT Catchability
  array<Type> Zed(nspp, max_age, nyrs); // Total mortality at age; n = [nspp, nages, nyrs]

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






  ///////////////////////
  DATA_VECTOR(Y);
  DATA_ARRAY(x);
  PARAMETER(ab);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));

  int n_data = Y.size();

  vector<Type> mu(n_data);
  vector<Type> x2(n_data);

  if(x2(1) == 1){
    std::cout << "x2 is " << sum(x2) << " and the code is broken \n";
    return EXIT_FAILURE;
  }
  x2.setZero();
  x2 += x.col(0).col(0);
  mu = (x2 * x.col(1).col(1));
  std::cout << "max mu is " << max(x2) << "\n";
  return 1;
  Type nll = -sum(dnorm(Y, mu, exp(logSigma), true));

  std::cout << "The nll is " << nll << "\n";
  ADREPORT(nll);
  return nll;
}
