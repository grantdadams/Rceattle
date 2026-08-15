/**
 * @file comp_sim.hpp
 * @brief Random draws for composition-style observations, used by SIMULATE.
 *
 * Every composition likelihood in the model -- age/length comps, conditional
 * age-at-length, and stomach contents -- is either a multinomial or a
 * Dirichlet-multinomial over a vector of bins. Simulation needs the matching
 * random draws, and TMB supplies neither: `distributions_R.hpp` stops at
 * univariate generators (`rbinom`, `rgamma`, `runif`, ...). These are the two
 * multivariate draws built from them, kept beside each other so a likelihood
 * cannot gain a family without its simulator being visibly absent.
 *
 * The draws are exact, not approximations, and both are the standard
 * constructions:
 *
 * - A multinomial is drawn as a sequence of conditional binomials. Having
 *   placed the first `a` bins, the count in bin `a+1` is binomial in the
 *   observations still unplaced, with probability rescaled by the probability
 *   still unspent. This costs one `rbinom` per bin. Drawing one categorical
 *   variate per observation instead (as WHAM's `rmultinom` does,
 *   src/age_comp_sim.hpp) is equally valid but costs one uniform per
 *   observation, which stomach sample sizes make expensive.
 * - A Dirichlet-multinomial is a multinomial whose probabilities are themselves
 *   a Dirichlet draw, and a Dirichlet is a vector of independent gammas divided
 *   by its own sum.
 *
 * Both return COUNTS summing to `N`, matching what the densities in
 * `helper_functions.hpp` and `comp_osa.hpp` read. The caller converts back to
 * whatever scale its data object stores.
 *
 * These run only inside `SIMULATE{}`, where `Type` is always double and there is
 * no AD tape, so `CppAD::Integer()` and the R generators underneath are safe.
 */

#ifndef RCEATTLE_COMP_SIM_HPP
#define RCEATTLE_COMP_SIM_HPP

/**
 * @brief Draw multinomial counts by sequential conditional binomials.
 *
 * @param N Number of observations to place. Rounded to a whole number, since
 *   the binomial size must be one; effective sample sizes reaching here are
 *   whole counts (stomachs sampled, fish aged) or a weighting of them.
 * @param p Bin probabilities. Rescaled internally, so they need not sum to 1,
 *   but must be non-negative.
 * @return Counts, one per bin, summing to the rounded `N`.
 */
template<class Type>
vector<Type> rmultinom_rce(Type N, vector<Type> p)
{
  int dim = p.size();
  vector<Type> x(dim);
  x.setZero();
  if(dim == 0) return x;

  // Round to a whole number for rbinom. Clamp first: CppAD::Integer saturates on
  // some platforms and wraps to INT_MIN on others, which would make n_left
  // negative and silently draw nothing.
  double n_req = asDouble(N) + 0.5;
  if(!(n_req > 0.0)) return x;
  if(n_req > 2147483000.0) n_req = 2147483000.0;
  Type n_left = Type(CppAD::Integer(n_req));
  Type p_left = p.sum();

  for(int a = 0; a < dim - 1; a++){
    if(!(n_left > Type(0))) break;          // everything already placed
    // Probability of landing in bin a given it did not land in bins 0..a-1.
    Type pa = (p_left > Type(0)) ? p(a) / p_left : Type(0);
    if(pa < Type(0)) pa = Type(0);
    if(pa > Type(1)) pa = Type(1);          // guards rounding at the last bin
    x(a) = rbinom(n_left, pa);
    n_left -= x(a);
    p_left -= p(a);
  }
  // The final bin is whatever is left: fixed by the sum-to-N constraint, never
  // drawn, exactly as the density treats it.
  if(n_left > Type(0)) x(dim - 1) = n_left;
  return x;
}


/**
 * @brief Draw Dirichlet proportions as normalized independent gammas.
 *
 * With a very small concentration every gamma can underflow to zero in double
 * precision, leaving nothing to normalize. The limit of a Dirichlet as
 * `alpha -> 0` is degenerate -- all the mass in ONE bin, chosen with probability
 * proportional to `alpha` -- so that is what is returned. Returning a flat
 * vector instead would be the opposite of the near-degenerate composition the
 * parameters describe, and would quietly turn a concentrated diet into a diffuse
 * one. Onset is around `sum(alpha)` below 1e-3.
 *
 * @param alpha Concentration parameters, one per bin, non-negative.
 * @return Proportions summing to 1.
 */
template<class Type>
vector<Type> rdirichlet_rce(vector<Type> alpha)
{
  int dim = alpha.size();
  vector<Type> x(dim);
  x.setZero();
  if(dim == 0) return x;

  for(int a = 0; a < dim; a++){
    x(a) = (alpha(a) > Type(0)) ? rgamma(alpha(a), Type(1)) : Type(0);
  }
  Type tot = x.sum();
  if(!(tot > Type(0))){
    // Everything underflowed: fall back to the degenerate limit, picking the
    // bin in proportion to alpha.
    Type a_tot = alpha.sum();
    x.setZero();
    if(a_tot > Type(0)){
      Type u = runif(Type(0), Type(1)) * a_tot;
      Type acc = 0.0;
      for(int a = 0; a < dim; a++){
        acc += alpha(a);
        if(u <= acc){ x(a) = Type(1); return x; }
      }
    }
    x(dim - 1) = Type(1);
    return x;
  }
  return x / tot;
}


/**
 * @brief Draw Dirichlet-multinomial counts.
 *
 * @param N Number of observations to place.
 * @param alpha Concentration parameters, one per bin. In this model these are
 *   `predicted proportion * N * theta`. The draw reproduces the MEAN the density
 *   predicts; it is not identical to the density in every moment, because the
 *   density is evaluated on offset data while the draw is taken on the raw
 *   scale.
 * @return Counts, one per bin, summing to the rounded `N`.
 */
template<class Type>
vector<Type> rdirmultinom_rce(Type N, vector<Type> alpha)
{
  vector<Type> p = rdirichlet_rce(alpha);
  return rmultinom_rce(N, p);
}

#endif
