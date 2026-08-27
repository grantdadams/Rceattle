// Standalone TMB model wrapping src/TMB/dsem.hpp's calculate_dsem().
//
// The vendored DSEM code has no numerical verification: every check so far has
// been compilation only. This exposes calculate_dsem() as an objective function
// on its own, so its jnll can be compared against dsem::dsem(run_model = TRUE)
// on identical inputs, and so the input assertions can be shown to fire.
//
// Used by tools/verify/verify-dsem-equivalence.R. Not part of the package build.
#include <TMB.hpp>
#include "dsem.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_IVECTOR( options );
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );
  DATA_IVECTOR( familycode_j );
  DATA_IVECTOR( linkcode_j );
  DATA_IVECTOR( sigmastart_j );
  DATA_ARRAY( eps_tj );
  DATA_ARRAY( y_tj );
  DATA_IVECTOR( obs_idx );
  DATA_IVECTOR( unobs_idx );
  DATA_IVECTOR( cond_k );        // 1 -> this CELL is given, k = j * n_t + t

  PARAMETER_VECTOR( beta_z );
  PARAMETER_VECTOR( lnsigma_z );
  PARAMETER_VECTOR( mu_j );
  PARAMETER_VECTOR( delta0_j );
  PARAMETER_ARRAY( x_tj );

  Type jnll = 0;
  array<Type> z_tj( y_tj.rows(), y_tj.cols() ); z_tj.setZero();
  // Variance of each latent state given the cells cond_k marks known. Reported
  // so it can be checked against a closed form -- for a first-order self-path
  // with nothing known it must approach sigma^2 / (1 - rho^2) away from the
  // series edges, and against diag(solve(Q[unknown, unknown])) when cells are
  // known. The recruitment bias correction in ceattle.cpp is built on this.
  array<Type> margvar_tj( y_tj.rows(), y_tj.cols() ); margvar_tj.setZero();
  matrix<Type> Q(0, 0);
  array<Type> xhat_tj( y_tj.rows(), y_tj.cols() ); xhat_tj.setZero();
  array<Type> delta_tj( y_tj.rows(), y_tj.cols() ); delta_tj.setZero();
  // Mean and SD of the measurement density, handed back so ceattle.cpp can draw
  // a covariate observation from the density that scores it. Unused here beyond
  // keeping the signature honest, but REPORTed so the equivalence harness can
  // check the link is applied where it should be.
  array<Type> mu_tj( y_tj.rows(), y_tj.cols() ); mu_tj.setZero();
  vector<Type> sigma_z( lnsigma_z.size() ); sigma_z.setZero();
  calculate_dsem( jnll, options, RAM, RAMstart, familycode_j, linkcode_j,
                  sigmastart_j, eps_tj, y_tj, obs_idx, unobs_idx,
                  beta_z, lnsigma_z, mu_j, delta0_j, x_tj, z_tj,
                  Q, xhat_tj, delta_tj, 0, 1, cond_k, margvar_tj,
                  mu_tj, sigma_z );
  REPORT( Q );
  REPORT( z_tj );
  REPORT( margvar_tj );
  REPORT( mu_tj );
  REPORT( sigma_z );
  return jnll;
}
