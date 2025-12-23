
// Function for getting max of an IVECTOR and return an int
template <class Type>
Type imax(const vector<Type> &x)
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

// Function to maintain > limit
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);   // tune a new candidate x given eps
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * (x - eps) * (x - eps));  // penalize via x vs eps
  return ans;
}


// Dirichlet multinomial (from WHAM)
template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type phi=sum(alpha);
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + alpha(a)) - lgamma(alpha(a));
  if(do_log == 1) return ll;
  else return exp(ll);
}

// From TMB examples, modified name, constrains between -1 and 1
template <class Type>
Type rho_trans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}


// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


// Function for detecting Inf
template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}


//  Function to getting sqaure
template <class Type> Type square(Type x){return x*x;}


// Function for mean of a vector
template <class Type>
Type mean(const vector<Type> &x)
{
  Type mean = x.sum()/x.size();
  return mean;
}

// Function to calculate the first difference
template <class Type>
vector<Type> first_difference(const vector<Type> &x)
{
  int length = x.size() - 1;

  if(length <= 0){
    std::cerr << "Error -- vector size too small in first_difference(const dvector&)" << std::endl;
    return(0);
  }

  vector<Type> tmp(length);

  for(int i = 0; i < length; i++){
    tmp(i) = x(i+1) - x(i);
  }

  return tmp;
}


// Function to get max of two values (differentiable)
// https://groups.google.com/g/tmb-users/c/QAVmgj66OC0
template <class Type>
Type max2(Type x, Type y){
  Type ans = 0.5*(CppAD::abs(x-y)+x+y);
  return ans;
}
