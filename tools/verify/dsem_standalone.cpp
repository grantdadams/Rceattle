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
  DATA_IVECTOR( cond_j );        // 1 -> condition on this column

  PARAMETER_VECTOR( beta_z );
  PARAMETER_VECTOR( lnsigma_z );
  PARAMETER_VECTOR( mu_j );
  PARAMETER_VECTOR( delta0_j );
  PARAMETER_ARRAY( x_tj );

  Type jnll = 0;
  array<Type> z_tj( y_tj.rows(), y_tj.cols() ); z_tj.setZero();
  // Marginal variance of each latent state. Reported so it can be checked
  // against a closed form -- for a first-order self-path it must approach
  // sigma^2 / (1 - rho^2) away from the series edges. The recruitment bias
  // correction in ceattle.cpp is built on this, and nothing else tests it.
  array<Type> margvar_tj( y_tj.rows(), y_tj.cols() ); margvar_tj.setZero();
  calculate_dsem( jnll, options, RAM, RAMstart, familycode_j, linkcode_j,
                  sigmastart_j, eps_tj, y_tj, obs_idx, unobs_idx,
                  beta_z, lnsigma_z, mu_j, delta0_j, x_tj, z_tj,
                  1, cond_j, margvar_tj );
  REPORT( z_tj );
  REPORT( margvar_tj );
  return jnll;
}
