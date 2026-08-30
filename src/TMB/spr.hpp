#ifndef SPR_HPP
#define SPR_HPP

/**
 * @file spr.hpp
 * @brief Spawning biomass per recruit, for the Amendment-56 SPR proxies.
 *
 * F40% and F35% are the F that leave a target fraction of unfished spawning
 * output per recruit; the template searches for them by penalising
 * SPR(F)/SPR(0) against that fraction (section 13, JNLL_REFPT_PENALTY).
 *
 * Both functions work on the female schedule (sex index 0), and report spawning
 * output per TOTAL recruit. The caller supplies the mature-female schedule with
 * the female fraction already in it -- the age-varying ratio for a one-sex
 * species (via `mature_females`, section 5.4), the recruitment split alone for a
 * two-sex species, whose schedule is female already.
 *
 * Not defined under predation: total mortality then carries M2, which scales
 * with predator abundance, so spawning output per recruit is not a property of
 * the prey stock alone. data_check() refuses the switch combinations that would
 * read it.
 */


/**
 * @brief Survivors at age from one recruit under a fixed mortality schedule.
 *
 * The oldest bin is a plus group, so it also accumulates the tail of survivors
 * that stay in it: `/ (1 - exp(-Z_plus))`.
 *
 * Two unchecked preconditions, both inherited from the code this replaced and
 * both left out of the AD path deliberately: `Z.size() >= 2`, since the plus
 * group reads the age below it, and a plus-group `Z` above zero, or that tail
 * divides by zero. No age-structured species reaches either.
 *
 * @param Z  Total instantaneous mortality at age (yr^-1), length nages >= 2.
 * @return   Numbers at age per recruit.
 */
template <class Type>
vector<Type> per_recruit_survivors(const vector<Type>& Z)
{
  int nages = Z.size();
  vector<Type> n(nages);
  n.setZero();
  n(0) = 1.0;

  for(int age = 1; age < nages - 1; age++){
    n(age) = n(age - 1) * exp(-Z(age - 1));
  }
  n(nages - 1) = n(nages - 2) * exp(-Z(nages - 2)) / (1.0 - exp(-Z(nages - 1)));

  return n;
}


/**
 * @brief Female spawning biomass per recruit, given a survivor schedule.
 *
 * @param n            Numbers at age per recruit, from per_recruit_survivors().
 * @param Z            The same mortality schedule that produced `n`, reused for
 *                     the mortality served before spawning within the year.
 * @param weight       Spawning weight at age (kg).
 * @param mature_female Proportion mature AND female at age -- `mature_females`
 *                     from section 5.4, which is maturity times the sex ratio
 *                     for a one-sex species and maturity alone for a two-sex
 *                     one, whose sex-0 schedule is already female.
 * @param spawn_month  Month of spawning, 0-12; 0 spawns at the start of the year.
 * @return             Spawning biomass per recruit (kg per recruit).
 */
template <class Type>
Type spawning_biomass_per_recruit(const vector<Type>& n,
                                  const vector<Type>& Z,
                                  const vector<Type>& weight,
                                  const vector<Type>& mature_female,
                                  Type spawn_month)
{
  Type spr = 0.0;
  for(int age = 0; age < n.size(); age++){
    spr += n(age) * weight(age) * mature_female(age) *
      exp(-Z(age) * spawn_month / 12.0);
  }
  return spr;
}

#endif
