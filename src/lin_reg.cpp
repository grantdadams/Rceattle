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
  DATA_MATRIX(m1);
  PARAMETER_VECTOR(a);
  PARAMETER(logSigma);


  int n_data = Y.size();
  vector<Type> mu(n_data);
  matrix<Type> m2(n_data, m1.cols());
  vector<Type> x1(n_data);

  m2 = pow(m1.array(), Type(2.0));
  m2 = m2.array() + 1;

  x1 = m1.col(0);

  mu = a(0) + x1 * a(1);

  Type nll = -sum(dnorm(Y, mu, exp(logSigma), true));

  std::cout << "The matrix 1,1 are " << m2(n_data-1,0) << "\n";


  ADREPORT(nll);
  ADREPORT(exp(2*logSigma));
  REPORT(m2);
  return nll;
}
