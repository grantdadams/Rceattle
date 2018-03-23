tmb_ceattle_model  
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  
  // Weight-at-Age Data
  DATA_INTEGER(N_SPP); // Number of species
  DATA_IVECTOR(N_YEARS); // Number of years of available weight-at-age data
  DATA_IVECTOR(YEAR_MIN); // First year of available weight-at-age data
  DATA_IVECTOR(n_weight_at_age); // Number of weight-at-age observations
  DATA_ARRAY(WEIGHT); // Weight in units?????
  DATA_ARRAY(AGE); // Age of 
  DATA_ARRAY(TEMP); // Temperature
  DATA_IVECTOR(YEAR); // Year of weight-at-length observation
  
  // VBGF Weight-at-Age Parameters
  PARAMETER_VECTOR(D_BASE_SCALAR); // Mean consumption intercept
  PARAMETER_VECTOR(D_BETA); // Coefficient of the residual effect of temperature on consumption
  PARAMETER_VECTOR(K_SCALAR); // Energy loss constant
  PARAMETER_VECTOR(H_SCALAR); // Assimilation constant
  PARAMETER_VECTOR(T0_SCALAR); // Age at weight 0 
  PARAMETER_VECTOR(sigma_VBGF); // SD of weight-at-age
  array<Type> D_MAT(N_YEARS, 1, N_SPP)
  
  // Estimate VBGF Params
  for(int i = 0; i < N_SPP; i++){
    for(int y = 0; y < N_YEARS(i); y++){
      D_MAT(y, 1, i) = exp(D_BASE_SCALAR(i) + D_BETA(i) * TEMP(y,1,i)); Estimate the year specific temperate allometric slop of consumption
      W_INF_VEC(y,i) = (H_SCALAR(i)/K_SCALAR(i)) ^ (1./(1. - D_MAT(y, 1, i))); //Derive the year specific asymptoptic weight
    }
  }
  
  // Fit VBGF
  vector<Type> vbgf_nll(N_SPP); // Initialize VBGF negative log-likelihood
  array<Type> pred_log_weight(n_weight_at_age,1,N_SPP); // Predicted weight
  int ind_year; // Indexing variable for year specific VBGF yarameters
  
  for(int i = 0; i < N_SPP; i++){
    for(int ind = 0; ind < n_weight_at_age(i); ind++){ // Loop through weight-at-age observations
      
      ind_year = YEAR(ind, 1, i) - YEAR_MIN(i); // Create indexing variable for VBGF Parameters
      
      pred_log_weight(ind, 1, i) = log(W_INF_VEC(ind_year, i) * (1 - exp(-K_SCALAR(i) * (1 - D_MAT(ind_year, 1, i))  * (AGE(ind) - T0_SCALAR(i)))) ^ (1 / (1 - D_MAT(ind_year, 1, i))));
      
      // Calculate negative log-likelihood
      vbgf_nll -= dnorm(log(WEIGHT(ind,1,i)), pred_log_weight(ind, 1, i), sigma_VBGF, true);
    }
  }
  
  // Report VBGF Parameters
  ADREPORT(K_SCALAR);
  ADREPORT(H_SCALAR);
  ADREPORT(T0_SCALAR);
  ADREPORT(D_BASE_SCALAR);
  ADREPORT(D_BETA);
  ADREPORT(vbgf_nll); //Report the likelihood component of weight-at-age
  
  }