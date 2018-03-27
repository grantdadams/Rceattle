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

  x2 = 0;
  std::cout << "x2 is " << x2 << "\n";
  x2 += x.col(0).col(0);
  mu = x2 * b + a;

  Type nll = -sum(dnorm(Y, mu, exp(logSigma), true));

  std::cout << "The nll is " << nll << "\n";
  ADREPORT(nll);
  return nll;
}
