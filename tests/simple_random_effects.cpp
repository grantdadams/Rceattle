#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() (){
  // Set-up
  DATA_VECTOR(x);
  DATA_IVECTOR(group);
  DATA_INTEGER(n_group);

  PARAMETER_VECTOR(re);
  PARAMETER(log_sigma_re);
  PARAMETER(log_sigma);
  PARAMETER(meanR);

  // Model objects
  Type sigma_re = exp(log_sigma_re);
  Type sigma = exp(log_sigma);
  Type jnll = 0;

  // Likelihood
  for(int j = 0; j < n_group; j++){
    jnll -= dnorm(re[j], Type(0), sigma_re, true);
  }

  for(int i = 0; i < x.size(); i++){
    jnll -= dnorm( x[i], meanR + re[group[i]-1], sigma, true );
  }

  vector<Type> rec = meanR + re;

  REPORT(sigma_re);
  REPORT(sigma);
  REPORT(meanR);
  REPORT(rec);

  ADREPORT(rec);
  ADREPORT(sigma);
  ADREPORT(sigma_re);

  return jnll;
}
