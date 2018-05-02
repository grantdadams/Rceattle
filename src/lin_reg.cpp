#include <TMB.hpp>
// Function for getting max of an Ivector

#include "../inst/include/functions.hpp"

template<class Type>
  Type objective_function<Type>::operator() (){
  DATA_VECTOR(Y);
  DATA_MATRIX(m1);
  DATA_ARRAY(a1);
  PARAMETER_VECTOR(a);
  PARAMETER(logSigma);


  int n_data = Y.size();
  vector<Type> mu(n_data);
  matrix<Type> m2(n_data, m1.cols());
  matrix<Type> m3(n_data, m1.cols());
  vector<Type> x1(n_data);

  m2 = pow(m1.array(), Type(2.0));
  m2 = m2.array() + 1;

  x1 = m1.col(0);

  mu = a(0) + x1 * a(1);

  Type nll = -sum(dnorm(Y, mu, exp(logSigma), true));

 // Get 2nd sheet of array a1
   m3 = array_to_matrix(a1, 1);

   // Check vector * matrix multiplication
   vector<Type> v1(2);
   v1 = Y * m3;

  std::cout << "The a1 ncols is " << a1.dim << "\n";


  ADREPORT(nll);
  ADREPORT(exp(2*logSigma));
  REPORT(m3);
  REPORT(v1);
  return nll;
}
