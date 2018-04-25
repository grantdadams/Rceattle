#include <TMB.hpp>
// Function for getting max of an Ivector
template <class Type>
matrix<Type> pow_mat(matrix<Type> m1, Type exponent){

  int m1r = m1.rows();
  int m1c = m1.cols();

  matrix<Type> m2(m1r, m1c);

  // Elementwise division
  for(int r = 0; r < m1r; r++){
    for(int c = 0; c < m1c; c++){
      m2(r, c) = pow( m1(r, c), exponent);
    }
  }
  return m2;
}

template<class Type>
  Type objective_function<Type>::operator() (){
  DATA_VECTOR(Y);
  DATA_IVECTOR(ints);
  DATA_MATRIX(x);
  PARAMETER_VECTOR(a);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));

  int n_data = Y.size();
  vector<Type> mu(n_data);
  matrix<Type> m2()
  vector<Type> x2(n_data);

  x2 = x.col(0);

  mu = a(0) + x2 * a(1);

  Type nll = -sum(dnorm(Y, mu, exp(logSigma), true));

  std::cout << "The nll is " << nll << "\n";
  ADREPORT(nll);
  return nll;
}
