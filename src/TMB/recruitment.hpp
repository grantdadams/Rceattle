#ifndef RECRUITMENT_HPP
#define RECRUITMENT_HPP

/**
 * @brief Calculates the recruitment for a given species and year.
 *
 * @param srr_fun   Stock-recruit relationship (SRR) type: 0 = random about mean (the
 *                  default), 2 = Beverton-Holt, 4 = Ricker. Codes 1, 3 and 5 are retired
 *                  aliases of 0, 2 and 4 that compute the same curve here; build_srr()
 *                  refuses them, because the environmental effect now comes from a
 *                  recruitment linkage on R0 (code 1) or alpha (3 and 5).
 * @param R0        Equilibrium or mean recruitment level (thousands of fish); read by the
 *                  random-about-mean forms only.
 * @param ssb       Female spawning-stock biomass (mt), lagged by minage by the caller.
 * @param srr_alpha SRR productivity parameter.
 * @param srr_beta  SRR density-dependence parameter.
 * @param rec_dev   Annual recruitment deviation, log scale.
 * @param spr0      Equilibrium spawning biomass per recruit at F = 0 (M of year 1).
 *                  Unused here: both curves are functions of srr_alpha, srr_beta and ssb,
 *                  and the caller derives the unfished recruitment they imply -- (alpha -
 *                  1/spr0)/beta for Beverton-Holt, a posfun()-guarded log for Ricker.
 * @return          Expected recruitment (thousands of fish).
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

  case 1: // Retired code: the same mean form as 0
    R = R0 * exp(rec_dev);
    break;

  case 2: // Beverton-Holt
    R = srr_alpha * ssb * exp(rec_dev) / (Type(1.0) + srr_beta * ssb);
    break;

  case 3: // Retired code: the same Beverton-Holt curve as 2
    R = srr_alpha * ssb * exp(rec_dev) / (Type(1.0) + srr_beta * ssb);
    break;

  case 4: // Ricker
    // Beta is estimated divided by 1,000,000 for estimation stability, so the fitted beta
    // is 1e6 times the curve's coefficient.
    R = srr_alpha * ssb * exp(-srr_beta * ssb / Type(1000000.0)) * exp(rec_dev);
    break;

  case 5: // Retired code: the same Ricker curve as 4
    R = srr_alpha * ssb * exp(-srr_beta * ssb / Type(1000000.0)) * exp(rec_dev);
    break;

  default:
    error("Invalid 'srr_fun'");
  }

  return R;
}

#endif
