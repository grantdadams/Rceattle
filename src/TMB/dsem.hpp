/** * @brief Dynamic Structural Equation Model (DSEM) Module
 * * This module assembles the precision matrix for a Gaussian Markov Random Field (GMRF)
 * using the Reticular Action Model (RAM) specification. It handles SEM path
 * coefficients, sparsity-optimized inversion, and multiple likelihood families.
 * * @param jnll [ref] Joint negative log-likelihood to be updated.
 * @param RAM Reticular Action Model specification matrix.
 * @param RAMstart Starting values or fixed values for RAM paths.
 * @param familycode_j Likelihood family codes for each variable j.
 * @param y_tj Observed data array (time t, variable j).
 * @param x_tj Latent state array (random effects).
 * @param beta_z Estimated path coefficients.
 * @param lnsigma_j Log-standard deviations for likelihood families.
 * @param mu_j Mean offsets for each variable j.
 * @param delta0_j Initial condition offsets.
 * @param options Vector of configuration flags for rank and variance scaling.
 */
template<class Type>
void calculate_dsem(
    Type &jnll, // Modified by reference
    matrix<int> RAM,
    vector<Type> RAMstart,
    vector<int> familycode_j,
    array<Type> y_tj,
    array<Type> x_tj,
    vector<Type> beta_z,
    vector<Type> lnsigma_j,
    vector<Type> mu_j,
    vector<Type> delta0_j,
    vector<int> options
) {
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  using namespace Eigen;

  // Copied from https://github.com/James-Thorson-NOAA/dsem/blob/main/src/dsem.cpp
  // Thorson, J. T., Andrews, A. G., Essington, T., & Large, S. (2024). Dynamic structural equation models synthesize ecosystem dynamics constrained by ecological mechanisms. Methods in Ecology and Evolution 15(4): 744-755. https://doi.org/10.1111/2041-210X.14289

  // Indices and Dimensions
  int n_t = y_tj.rows();
  int n_j = y_tj.cols();
  int n_k = n_t * n_j;

  Type jnll_gmrf = 0;
  matrix<Type> loglik_tj(n_t, n_j); loglik_tj.setZero();
  vector<Type> sigma_j = exp(lnsigma_j);

  // Assemble Precision Matrices
  Eigen::SparseMatrix<Type> Rho_kk(n_k, n_k); Rho_kk.setZero();
  Eigen::SparseMatrix<Type> Gamma_kk(n_k, n_k); Gamma_kk.setZero();
  Eigen::SparseMatrix<Type> I_kk(n_k, n_k); I_kk.setIdentity();

  for (int r = 0; r < RAM.rows(); r++) {
    // Extract estimated or fixed value
    Type tmp = (RAM(r, 3) >= 1) ? beta_z(RAM(r, 3) - 1) : RAMstart(r);
    if (RAM(r, 0) == 1) Rho_kk.coeffRef(RAM(r, 1) - 1, RAM(r, 2) - 1) = tmp;
    if (RAM(r, 0) == 2) Gamma_kk.coeffRef(RAM(r, 1) - 1, RAM(r, 2) - 1) = tmp; // Cholesky of covariance, so -Inf to Inf;
  }

  // Compute inverse LU-decomposition
  Eigen::SparseMatrix<Type> IminusRho_kk = I_kk - Rho_kk;
  Eigen::SparseLU<Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int>> inverseIminusRho_kk;
  inverseIminusRho_kk.compute(IminusRho_kk);

  // Marginal Variance Rescaling Logic
  if (options(1) == 1 || options(1) == 2) {
    Eigen::SparseMatrix<Type> invIminusRho_kk = inverseIminusRho_kk.solve(I_kk); // WORKS:  Based on: https://github.com/kaskr/adcomp/issues/74

    // Hadamard squared LU-decomposition
    // See: https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
    Eigen::SparseMatrix<Type> squared_invIminusRho_kk = invIminusRho_kk.cwiseProduct(invIminusRho_kk);
    Eigen::SparseLU<Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int>> invsquared_invIminusRho_kk;
    invsquared_invIminusRho_kk.compute(squared_invIminusRho_kk);

    if (options(1) == 1) {
      matrix<Type> ones_k1(n_k, 1); ones_k1.setOnes();

      // Calculate diag( t(Gamma) * Gamma )
      Eigen::SparseMatrix<Type> squared_Gamma_kk = Gamma_kk.cwiseProduct(Gamma_kk);
      matrix<Type> sigma2_k1 = squared_Gamma_kk.transpose() * ones_k1;
      matrix<Type> margvar_k1 = invsquared_invIminusRho_kk.solve(sigma2_k1);

      // Rescale IminusRho_kk and Gamma
      Eigen::SparseMatrix<Type> invmargsd_kk(n_k, n_k);
      Eigen::SparseMatrix<Type> invsigma_kk(n_k, n_k);
      for (int k = 0; k < n_k; k++) {
        invmargsd_kk.coeffRef(k, k) = pow(margvar_k1(k, 0), -0.5);
        invsigma_kk.coeffRef(k, k) = pow(sigma2_k1(k, 0), -0.5);
      }
      IminusRho_kk = invmargsd_kk * IminusRho_kk;
      Gamma_kk = invsigma_kk * Gamma_kk;

      // Recompute inverse LU-decomposition
      inverseIminusRho_kk.compute(IminusRho_kk);
    } else {
      // calculate diag(Gamma)^2
      matrix<Type> targetvar_k1(n_k, 1);
      for (int k = 0; k < n_k; k++) targetvar_k1(k, 0) = pow(Gamma_kk.coeffRef(k, k), 2);

      // Rescale Gamma
      matrix<Type> margvar_k1 = invsquared_invIminusRho_kk.solve(targetvar_k1);
      for (int k = 0; k < n_k; k++) Gamma_kk.coeffRef(k, k) = pow(margvar_k1(k, 0), 0.5);
    }
  }

  // Initial Conditions (Delta)
  vector<Type> delta_k(n_k); delta_k.setZero();
  if (delta0_j.size() > 0) {
    matrix<Type> delta0_k1(n_k, 1); delta0_k1.setZero();
    for (int j = 0; j < n_j; j++) delta0_k1(j * n_t, 0) = delta0_j(j);
    matrix<Type> x = inverseIminusRho_kk.solve(delta0_k1);
    delta_k = x.array();
  }

  // GMRF Density Calculation
  array<Type> xhat_tj(n_t, n_j);
  array<Type> delta_tj(n_t, n_j);
  for (int j = 0; j < n_j; j++) {
    for (int t = 0; t < n_t; t++) {
      xhat_tj(t, j) = mu_j(j);
      delta_tj(t, j) = delta_k(j * n_t + t);
    }
  }

  // Apply GMRF
  array<Type> z_tj(n_t, n_j);
  if (options(0) == 0) { // Full Rank
    Eigen::SparseMatrix<Type> V_kk = Gamma_kk.transpose() * Gamma_kk;
    matrix<Type> Vinv_kk = invertSparseMatrix(V_kk);
    Eigen::SparseMatrix<Type> Q_kk = IminusRho_kk.transpose() * asSparseMatrix(Vinv_kk) * IminusRho_kk;

    // Centered GMRF
    jnll_gmrf = GMRF(Q_kk)(x_tj - xhat_tj - delta_tj);
    z_tj = x_tj;
  } else if (options(0) == 1) { // Rank Deficient
    jnll_gmrf = GMRF(I_kk)(x_tj);

    // Forward-format matrix
    matrix<Type> z_k1(n_k, 1);
    for (int j = 0; j < n_j; j++){
      for (int t = 0; t < n_t; t++){
        z_k1(j * n_t + t, 0) = x_tj(t, j);
      }
    }

    // (I-Rho)^{-1} * Gamma * Epsilon
    matrix<Type> z3_k1 = inverseIminusRho_kk.solve(matrix<Type>(Gamma_kk * z_k1));

    // Back-format vector
    for (int j = 0; j < n_j; j++){
      for (int t = 0; t < n_t; t++){
        z_tj(t, j) = z3_k1(j * n_t + t, 0);
      }
    }

    // Add back mean and deviation
    z_tj += xhat_tj + delta_tj;
  }

  // Likelihoods for family types
  array<Type> mu_tj( n_t, n_j );
  for (int t = 0; t < n_t; t++) {
    for (int j = 0; j < n_j; j++) {
      mu_tj(t, j) = (familycode_j(j) == 2) ? invlogit(z_tj(t, j)) :
      (familycode_j(j) == 3 || familycode_j(j) == 4) ? exp(z_tj(t, j)) : z_tj(t, j);

      if (!R_IsNA(asDouble(y_tj(t, j))) && familycode_j(j) != 0) {
        if (familycode_j(j) == 1) loglik_tj(t, j) = dnorm(y_tj(t, j), mu_tj(t, j), sigma_j(j), true);
        if (familycode_j(j) == 2) loglik_tj(t, j) = dbinom(y_tj(t, j), Type(1.0), mu_tj(t, j), true);
        if (familycode_j(j) == 3) loglik_tj(t, j) = dpois(y_tj(t, j), mu_tj(t, j), true);
        if (familycode_j(j) == 4) loglik_tj(t, j) = dgamma(y_tj(t, j), pow(sigma_j(j), -2), mu_tj(t, j) * pow(sigma_j(j), 2), true);
      }
    }
  }
  jnll -= loglik_tj.sum();
  jnll += jnll_gmrf;

  // Reporting
  // REPORT( xhat_tj ); // needed to simulate new GMRF in R
  // REPORT( delta_k ); // FIXME>  Eliminate in simulate.dsem
  // REPORT( delta_tj ); // needed to simulate new GMRF in R
  // REPORT( Rho_kk );
  // REPORT( Gamma_kk );
  // //REPORT( mu_tj );
  // //REPORT( devresid_tj );
  // REPORT( IminusRho_kk );
  // REPORT( jnll_dsem );
  // REPORT( loglik_tj_dsem );
  // REPORT( jnll_gmrf_dsem );
  // SIMULATE{
  //   REPORT( y_tj );
  // }
  // REPORT( z_tj );
  // ADREPORT( z_tj );
}
