/**
 * @file comp_osa.hpp
 * @brief One-step-ahead (OSA) residual helpers for composition data.
 *
 * TMB's oneStepPredict() requires every observation entering the likelihood to
 * live in a flat `obsvec` with a `keep` indicator, and each likelihood term to
 * be multiplied by its keep element. Multinomial and Dirichlet-multinomial
 * composition data are multivariate and discrete, so oneStepPredict cannot use
 * them directly. Following Trijoulet et al. (2023) and WHAM
 * (src/age_comp_osa.hpp), a composition of A bins is decomposed into A-1
 * univariate conditional distributions -- binomials for the multinomial,
 * beta-binomials for the Dirichlet-multinomial -- each gated by its keep
 * element. The final bin is fixed by the sum-to-N constraint and contributes
 * nothing (its residual is not defined and is reported as NA on the R side).
 *
 * Key property: with keep == 1 the sum of the conditional log-densities equals
 * the joint multinomial / Dirichlet-multinomial log-density. Turning OSA on
 * therefore changes how the likelihood is *decomposed*, not its value.
 *
 * These functions are named with an `_osa` suffix and take a data_indicator so
 * they never collide with TMB's built-in dmultinom() or the plain
 * ddirmultinom() in helper_functions.hpp used during ordinary fitting.
 */

// Reorder categories so kept (keep > 0.5) entries come first, preserving the
// original order within the kept and dropped groups.
template<class Type>
vector<int> osa_order(vector<Type> k){
  int n = k.size();
  vector<int> ord(n);
  ord.setZero();
  int at = -1;
  for(int i = 0; i < n; ++i){
    if(k(i) > 0.5){ ord(++at) = i; }
  }
  at = n;
  for(int i = n - 1; i >= 0; --i){
    if(k(i) < 0.5){ ord(--at) = i; }
  }
  return ord;
}

// Map the closed [0, 1] interval onto the open interval for numerical stability.
// Identical to TMB's convenience.hpp squeeze() (eps = machine epsilon), which is
// what WHAM uses in its OSA composition likelihood, so the two implementations
// agree to the last digit.
template<class Type>
Type osa_squeeze(Type u){
  Type eps = std::numeric_limits<double>::epsilon();
  return (Type(1.0) - eps) * (u - Type(0.5)) + Type(0.5);
}

/**
 * @brief Keep-aware multinomial log-density via conditional binomials.
 *
 * With `do_osa`, decomposes the A-bin multinomial into A-1 conditional binomials
 * (each gated by its keep element) so oneStepPredict can consume it; the sum
 * equals the joint multinomial log-density when all keeps are 1.
 *
 * @param x Observed counts (length A; may be non-integer, e.g. Neff * prop).
 * @param p Predicted proportions (length A; normalized internally).
 * @param keep Data indicator for these A bins (from keep.segment(...)).
 * @param do_osa 1 = conditional-binomial decomposition; 0 = plain multinomial.
 * @param give_log 1 = return the log-density.
 */
template <class Type>
Type dmultinom_osa(vector<Type> x, vector<Type> p,
                   data_indicator<vector<Type>, Type> keep,
                   int do_osa, int give_log)
{
  Type logres = 0;
  // Normalize predicted proportions. Guard the denominator so a length/age bin
  // with all-zero predicted mass (sum(p) == 0, e.g. a growth-transition row that
  // is exactly zero in the tails) yields p_x == 0 rather than 0/0 = NaN; the
  // conditional binomials below then read each squeezed-to-eps probability,
  // matching WHAM's behaviour for vanishing predicted proportions.
  vector<Type> p_x = p / (sum(p) + Type(1e-300));

  if(do_osa){
    vector<Type> k = keep;
    vector<int> o = osa_order(k);
    x = x(o); p_x = p_x(o); k = k(o);     // kept bins first
    Type nUnused = asDouble(sum(x));
    Type pUsed = 0;
    for(int i = 0; i < x.size(); ++i){
      if(i != (x.size() - 1)){
        // Bin i versus all remaining bins: a binomial written as a 2-category
        // multinomial. Conditional success probability p_i / (1 - pUsed).
        vector<Type> x2(2), p2(2);
        x2(0) = x(i);
        x2(1) = nUnused - x(i);
        p2(0) = osa_squeeze(p_x(i)) / osa_squeeze(Type(1) - pUsed);
        p2(0) = osa_squeeze(p2(0));
        // When one bin holds ~all remaining mass the separately-squeezed
        // numerator/denominator can push the ratio marginally above 1; clamp so
        // p2(1) = 1 - p2(0) stays non-negative (avoids dmultinom log(negative) =
        // NaN). A no-op for well-behaved compositions.
        Type eps_clamp = std::numeric_limits<double>::epsilon();
        p2(0) = CppAD::CondExpGt(p2(0), Type(1) - eps_clamp, Type(1) - eps_clamp, p2(0));
        p2(1) = Type(1) - p2(0);
        logres += k(i) * dmultinom(x2, p2, true);
        nUnused -= x(i);
        pUsed   += p_x(i);
        // Keep pUsed inside (0, 1) so 1 - pUsed never reaches 0 or goes
        // negative. WHAM does NOT squeeze pUsed itself; it relies on the
        // squeeze of 1 - pUsed above, which is also applied here. This extra
        // squeeze is therefore redundant protection and contracts pUsed toward
        // 0.5 by one machine epsilon per bin, so a composition of A bins sits
        // O(A * eps) from WHAM's value -- ~1e-15 relative, far below the
        // reference tolerances, but not bit-identical. Removing it would match
        // WHAM exactly at the cost of re-pinning the OSA references.
        pUsed    = osa_squeeze(pUsed);
      } else {
        logres += k(i) * Type(0);         // last bin: fixed by sum-to-N
      }
    }
  } else {
    logres = dmultinom(x, p_x, true);
  }
  return give_log ? logres : exp(logres);
}

/**
 * @brief Keep-aware Dirichlet-multinomial log-density via conditional beta-binomials.
 *
 * With `do_osa`, decomposes the A-bin Dirichlet-multinomial into A-1 conditional
 * beta-binomials (each gated by its keep element); the sum equals the joint D-M
 * log-density when all keeps are 1.
 *
 * @param obs Observed counts (length A).
 * @param alpha Dirichlet-multinomial concentration parameters (length A).
 * @param keep Data indicator for these A bins.
 * @param do_osa 1 = conditional beta-binomial decomposition; 0 = plain D-M.
 * @param give_log 1 = return the log-density.
 */
template<class Type>
Type ddirmultinom_osa(vector<Type> obs, vector<Type> alpha,
                      data_indicator<vector<Type>, Type> keep,
                      int do_osa, int give_log)
{
  Type ll = 0.0;
  if(do_osa){
    vector<Type> k = keep;
    vector<int> o = osa_order(k);
    obs = obs(o); alpha = alpha(o); k = k(o);   // kept bins first
    int dim = obs.size();
    Type alp_sum = alpha.sum();
    Type obs_sum = asDouble(obs.sum());
    for(int a = 0; a < dim - 1; ++a){
      // Bin a versus all remaining bins: a beta-binomial written as a
      // 2-category Dirichlet-multinomial (alpha_a vs sum of remaining alphas).
      obs_sum -= obs(a);
      alp_sum -= alpha(a);
      vector<Type> obs_a(2), alphas_a(2);
      obs_a(0)    = obs(a);   obs_a(1)    = obs_sum;
      alphas_a(0) = alpha(a); alphas_a(1) = alp_sum;
      ll += k(a) * ddirmultinom(obs_a, alphas_a, 1);
    }
  } else {
    ll = ddirmultinom(obs, alpha, 1);
  }
  return give_log ? ll : exp(ll);
}
