#include <TMB.hpp>
// Function for getting max of an Ivector
template <class Type>
Type max2(const vector<Type> &x)
{
  int res = x[0];
  for(int i=0; i < x.size(); i++){
    int res2 = x[i];
    if(res < res2){
      res = res2;
      }
  }
  return res;
}

template <class Type>
Type mean_vec(const vector<Type> &x)
{
  Type mean = x.sum()/x.size();
  return mean;
}

template<class Type>
  Type objective_function<Type>::operator() (){
  DATA_VECTOR(Y);
  DATA_IVECTOR(ints);
  DATA_ARRAY(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));

  int n_data = Y.size();


  vector<Type> mu(n_data);
  vector<Type> x2(n_data);

  x2.setZero();
  x2 += x.col(0).col(0);

  mu = a + x2 * b;

  if(x2(1) < 1){

  }

  std::cout << "max x2 is " << mean_vec(ints) << "\n";
  return 1;
  Type nll = -sum(dnorm(Y, log(mu), exp(logSigma), true));

  std::cout << "The nll is " << nll << "\n";
  ADREPORT(nll);
  return nll;
}