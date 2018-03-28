#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() (){
  DATA_VECTOR(Y);
  DATA_ARRAY(x);
  PARAMETER(a);
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
  mu = x2 * x.col(1).col(1);
  std::cout << "mu is " << mu(1) << "\n";
  return 1;
  Type nll = -sum(dnorm(Y, mu, exp(logSigma), true));

  std::cout << "The nll is " << nll << "\n";
  ADREPORT(nll);
  return nll;
}
