#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() (){
  DATA_VECTOR(Y);
  DATA_ARRAY(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));

  vector<Type> mu;
  vector<Type> x2;
  x2 = x.col(0).col(0);
  mu = x2 * b + a;

  Type nll = -sum(dnorm(Y, mu, exp(logSigma), true));

  std::cout << "The nll is " << nll << "\n";
  ADREPORT(nll);
  return nll;
}
