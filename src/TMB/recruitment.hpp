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
 * @param srr_alpha    First SRR parameter (e.g., alpha/productivity parameter).
 * @param srr_beta     Second SRR parameter (e.g., beta/density-dependence parameter).
 * @param rec_dev      Annual log recruitment deviation.
 * @param spr0         Equilibrium SPR at F = 0 (assume M in yr 1)
 * * @return             Calculated expected recruitment (Type).
 */
template <class Type>
Type calculate_recruitment(int srr_fun,
                           Type R0,
                           Type ssb,
                           Type srr_alpha,
                           Type srr_beta,
                           Type rec_dev,
                           Type spr0) {
  Type R = 0;

  switch(srr_fun) {
  case 0: // Random about mean
    R = R0 * exp(rec_dev);
    break;

  case 1: // Random about mean with environmental effects
    R = R0 * exp(rec_dev);
    break;

  case 2: // Beverton-Holt
    R = srr_alpha * ssb * exp(rec_dev) / (Type(1.0) + srr_beta * ssb);
    // Recruitment depends only on alpha, beta and SSB. The implied unfished
    // recruitment, (alpha - 1/spr0)/beta, is derived by the caller.
    break;

  case 3: // Beverton-Holt with environmental impacts on alpha
    R = srr_alpha * ssb * exp(rec_dev) / (Type(1.0) + srr_beta * ssb);
    // Recruitment depends only on alpha, beta and SSB. The implied unfished
    // recruitment, (alpha - 1/spr0)/beta, is derived by the caller.
    break;

  case 4: // Ricker
    // Beta is divided by 1,000,000 for estimation stability
    R = srr_alpha * ssb * exp(-srr_beta * ssb / Type(1000000.0)) * exp(rec_dev);
    // Recruitment depends only on alpha, beta and SSB. The implied unfished
    // recruitment is derived by the caller, which guards the log() with posfun().
    break;

  case 5: // Ricker with environmental impacts on alpha
    R = srr_alpha * ssb * exp(-srr_beta * ssb / Type(1000000.0)) * exp(rec_dev);
    // Recruitment depends only on alpha, beta and SSB. The implied unfished
    // recruitment is derived by the caller, which guards the log() with posfun().
    break;

  default:
    error("Invalid 'srr_fun'");
  }

  return R;
}

#endif
