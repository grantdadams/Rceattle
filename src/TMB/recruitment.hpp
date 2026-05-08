#ifndef RECRUITMENT_HPP
#define RECRUITMENT_HPP

/**
 * @brief Calculates the recruitment for a given species and year.
 *
 * This function supports multiple stock-recruit relationships (SRR), including
 * random about mean, Beverton-Holt, and Ricker, with optional environmental covariates.
 *
 * @param srr_fun      Integer switch for the SRR type:
 * 0 = Random about mean (e.g., Alaska default)
 * 1 = Random about mean with environmental linkage
 * 2 = Beverton-Holt
 * 3 = Beverton-Holt with environmental impacts on alpha
 * 4 = Ricker
 * 5 = Ricker with environmental impacts on alpha
 * @param R0           Equilibrium recruitment (e.g., unfished or initial recruitment).
 * @param ssb          Spawning stock biomass (must be pre-lagged by minage before passing).
 * @param rec_par_1    First SRR parameter (e.g., alpha/productivity parameter).
 * @param rec_par_2    Second SRR parameter (e.g., beta/density-dependence parameter).
 * @param rec_dev      Annual log recruitment deviation.
 * @param srr_env_mult Pre-calculated environmental multiplier (dot product of env_index and betas).
 * Pass 0.0 if not using an environmentally linked SRR.
 * * @return             Calculated expected recruitment (Type).
 */
template <class Type>
Type calculate_recruitment(int srr_fun,
                           Type R0,
                           Type ssb,
                           Type rec_par_1,
                           Type rec_par_2,
                           Type rec_dev,
                           Type srr_env_mult) {
  Type R = 0;
  Type srr_alpha = 0;

  switch(srr_fun) {
  case 0: // Random about mean
    R = R0 * exp(rec_dev);
    break;

  case 1: // Random about mean with environmental effects
    R = R0 * exp(rec_dev + srr_env_mult);
    break;

  case 2: // Beverton-Holt
    R = exp(rec_par_1) * ssb * exp(rec_dev) / (Type(1.0) + exp(rec_par_2) * ssb);
    break;

  case 3: // Beverton-Holt with environmental impacts on alpha
    srr_alpha = exp(rec_par_1 + srr_env_mult);
    R = srr_alpha * ssb * exp(rec_dev) / (Type(1.0) + exp(rec_par_2) * ssb);
    break;

  case 4: // Ricker
    // Beta is divided by 1,000,000 for estimation stability
    R = exp(rec_par_1) * ssb * exp(-exp(rec_par_2) * ssb / Type(1000000.0)) * exp(rec_dev);
    break;

  case 5: // Ricker with environmental impacts on alpha
    srr_alpha = exp(rec_par_1 + srr_env_mult);
    R = srr_alpha * ssb * exp(-exp(rec_par_2) * ssb / Type(1000000.0)) * exp(rec_dev);
    break;

  default:
    error("Invalid 'srr_fun'");
  }

  return R;
}

#endif
