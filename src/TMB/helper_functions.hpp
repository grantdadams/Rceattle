
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


// Function to assemble sparse precision matrix
template<class Type>
// @description: Function that constructs a precision matrix, separable along the
// year, age, and cohort axis. Var_Param allows users to switch between conditional
// variance, and marginal variance.
Eigen::SparseMatrix<Type> construct_Q(int n_years, // Integer of years
                                      int n_ages, // Integer of ages
                                      matrix<Type> ay_Index, // Index matrix to construct
                                      Type rho_y, // Partial correlation by years
                                      Type rho_a, // Partial correlation by ages
                                      Type rho_c, // Partial correlation by cohort
                                      Type log_sigma2, // Variance parameter governing GMRF
                                      int Var_Param // Parameterization of Variance ==0 (Conditional), == 1(Marginal)
) {

  // Dimension to construct matrices
  int total_n = n_years * n_ages;

  // Construct matrices for precision matrix
  Eigen::SparseMatrix<Type> B(total_n,total_n); // B matrix
  Eigen::SparseMatrix<Type> I(total_n,total_n); // Identity matrix
  I.setIdentity(); // Set I to identity matrix
  Eigen::SparseMatrix<Type> Omega(total_n,total_n); // Omega matrix (variances)
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix

  for(int n = 0; n < total_n; n++) {

    // Define year and age objects
    Type age = ay_Index(n,0);
    Type year = ay_Index(n,1);

    // Constructing B matrix to determine where the correlation pars should go
    if(age > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year)
          B.coeffRef(n, n1) = rho_y;
      } // n1 loop

    } // end age > 1

    if(year > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n, n1) = rho_a;
      } // n1 loop

    } // if year > 1

    if(year > 1 && age > 1) {

      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year - 1)
          B.coeffRef(n,n1) = rho_c; // correlation by cohort
      } // n1 loop

    } // if both year and age > 1

  } // end n loop

  // Fill in Omega matrix here (variances)
  if(Var_Param == 0) { // Conditional variance

    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/exp(log_sigma2);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if conditional variance

  if(Var_Param == 1) { // Marginal Variance

    // Construct container objects
    matrix<Type> L(total_n, total_n); // L Matrix
    matrix<Type> tmp_I_B = I-B; // Temporary Matrix to store I-B
    L =  tmp_I_B.inverse(); // Invert to get L
    vector<Type> d(total_n); // Store variance calculations

    for(int n = 0; n < total_n; n++) {
      if(n == 0) {
        d(n) = exp(log_sigma2); // marginal variance parameter
      } else{

        Type cumvar = 0; // Cumulative Variance Container

        for(int n1 = 0; n1 < n; n1++) {
          cumvar += L(n, n1) * d(n1) * L(n, n1);
        } // n1 loop

        // Calculate diagonal values for omega
        d(n) = (exp(log_sigma2) - cumvar) / pow(L(n, n), 2);

      } // else loop
    } // n loop

    // Now fill in our diagonals for Omega
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/d(i);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop

  } // end if marginal variance

  // Now, do calculations to construct (Q = (I - t(B)) %*% Omega %*% (I-B))
  Eigen::SparseMatrix<Type> B_transpose = B.transpose(); // transpose B matrix

  // Calculate Precision Matrix
  Q_sparse = (I - B_transpose) * Omega * (I-B);

  return(Q_sparse);

} // end construct_Q function

