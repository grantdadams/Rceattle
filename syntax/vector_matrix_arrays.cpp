// Overview of customized vector, matrix and array operations in TMB.

#include <TMB.hpp>
#include "../src/include/functions.hpp"

template<class Type>
  Type objective_function<Type>::operator() (){
  DATA_VECTOR(v1);
  DATA_IVECTOR(iv1);
  DATA_MATRIX(m1);
  DATA_ARRAY(a1);

  // Vectors -----------------------------------


  vector<Type> v2(v1.size()-1); v2 = first_difference(v1); // First difference of vector i.e. c(v1(0) - v1(1), v1(1) - v1(2), ..., v1(n-1) - c1(n))
  REPORT(v2);


  vector<Type> v3(m1.cols()); v3 = vec_mat_prod(v1, m1); // Multiplication of a vector and matrix: "v %*% m" in r
  REPORT(v3);


  // IVectors ----------------------------------


  int i1 = imax(iv1); // Max of an ivector
  REPORT(i1);


  int i2 = imin(iv1); // Min of an ivector
  REPORT(i2);


  // Matrices ----------------------------------


  matrix<Type> m2(m1.rows(), m1.cols()); m2 = elem_pow(m1, Type(2.0)); // Elementwise power of matrix: pow(m1.array(), 2) also works
  REPORT(m2);


  //matrix<Type> m3(m1.rows(), m1.cols()); m3 = elem_div(m1, m1); // Elementwise division of two matrices: m1.array()/m1.array() is equivalent
  //REPORT(m3);


  // Arrays ------------------------------------


  matrix<Type> m4(a1.dim(0), a1.dim(1)); m4 = matrix_from_array(a1, 1); // Select sheet of array
  REPORT(m4);


  return 0;
}
