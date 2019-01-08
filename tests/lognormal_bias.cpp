#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() (){
  // Set-up
  DATA_VECTOR(x);
  DATA_IVECTOR(group);
  DATA_INTEGER(n_group);
  DATA_INTEGER( mod );

  PARAMETER_VECTOR(log_re);
  PARAMETER(log_sigma_re);
  PARAMETER(log_sigma);
  PARAMETER(log_meanR);

  // Model objects
  Type sigma_re = exp(log_sigma_re);
  Type sigma = exp(log_sigma);
  Type jnll = 0;

  // Random effects likelihood
  for(int j = 0; j < n_group; j++){
    if(mod == 1){
      jnll -= dnorm(log_re[j], Type(0), sigma_re, true);
    }
    if(mod == 2){
      jnll -= dnorm(log_re[j] - (sigma_re * sigma_re)/2, Type(0), sigma_re, true);
    }

  }

  vector<Type> rec(group.size());
  for(int i = 0; i < x.size(); i++){
    if(mod == 1){
      rec[i] = exp(log_meanR + log_re[group[i]] - (sigma_re * sigma_re) / 2);
    }

    if(mod == 2){
      rec[i] = exp(log_meanR + log_re[group[i]]);
    }
  }

  // Likelihood
  for(int i = 0; i < x.size(); i++){
    jnll -= dnorm( x[i], rec[i], sigma, true );
  }

  REPORT(sigma_re);
  REPORT(sigma);
  REPORT(log_meanR);
  REPORT(rec);

  ADREPORT(rec);
  ADREPORT(sigma);
  ADREPORT(sigma_re);

  return jnll;
}
