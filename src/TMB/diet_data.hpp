#ifndef DIET_CONTENT_HPP
#define DIET_CONTENT_HPP

/**
 * @brief Parses empirical diet observations and maps them into a 5D proportion array.
 *
 * @details Maps the observed proportions in `diet_obs` onto the 5D `diet_prop` array
 * used for MSVPA based suitability, handling 1-sex and 2-sex indexing. A row with
 * `flt_yr == 0` is written to every modeled year. Rows aggregated over predator or
 * prey ages are skipped, since only diet-data type 1 (a prey age in a predator age)
 * can be used for MSVPA suitability.
 *
 * @param nspp Number of modeled species.
 * @param nyrs Total number of years in the model (hindcast + projection).
 * @param nyrs_hind Number of years in the hindcast period.
 * @param styr The starting year of the model.
 * @param minage Vector containing the minimum modeled age for each species.
 * @param nsex Vector containing the number of sexes modeled for each species (1 or 2).
 * @param diet_obs Empirical diet observations; column 0 sample size, column 1 the proportion.
 * @param diet_ctl Control matrix mapping each observation row to predator, prey, sexes, ages, and year.
 * @param diet_prop [Output] 5D array (pred_idx, prey_idx, pred_age, prey_age, yr) of observed diet proportions. Modified by reference.
 */
template <class Type>
void organize_diet_obs(
    const int& nspp,
    const int& nyrs,
    const int& nyrs_hind,
    const int& styr,
    const vector<int>& minage,
    const vector<int>& nsex,
    matrix<Type>& diet_obs,
    matrix<int>& diet_ctl,
    array<Type>& diet_prop // Modified by reference
) {
  for(int stom_ind = 0; stom_ind < diet_obs.rows(); stom_ind++) {
    int rsp = diet_ctl(stom_ind, 0) - 1;             // Index of pred
    int ksp = diet_ctl(stom_ind, 1) - 1;             // Index of prey
    int r_sex = diet_ctl(stom_ind, 2);               // Pred sex: 0 combined, 1 female, 2 male
    int k_sex = diet_ctl(stom_ind, 3);               // Prey sex, same convention
    int r_age = diet_ctl(stom_ind, 4) - minage(rsp); // Index of pred age
    int k_age = diet_ctl(stom_ind, 5) - minage(ksp); // Index of prey age
    int flt_yr = diet_ctl(stom_ind, 6);              // Year of the row, or 0 for all years
    // A negative age marks a row aggregated over prey ages (k_age) or over predator
    // ages (r_age). Those are not diet-data type 1, so they carry no MSVPA
    // suitability; predict_stomach_content() still predicts them for the likelihood.
    if((k_age < 0) || (r_age < 0)) continue;

    // Sex looping bounds for predator; a combined-sex row is written to both sexes.
    int r_start = 0;
    int r_end = 0;
    if(nsex(rsp) == 2) {
      if(r_sex == 0) { r_end = 1; }
      else { r_start = r_sex - 1; r_end = r_sex - 1; }
    }

    // Sex looping bounds for prey
    int k_start = 0;
    int k_end = 0;
    if(nsex(ksp) == 2) {
      if(k_sex == 0) { k_end = 1; }
      else { k_start = k_sex - 1; k_end = k_sex - 1; }
    }

    Type obs_val = diet_obs(stom_ind, 1);

    for(int j = r_start; j <= r_end; j++) {
      int r_idx = rsp + (nspp * j);

      for(int k = k_start; k <= k_end; k++) {
        int k_idx = ksp + (nspp * k);

        // Only flt_yr >= 0 is loaded here: one year, dropped if outside the hindcast,
        // or year 0 for every year. A negative flt_yr loads nothing.
        if(flt_yr > 0) {
          int yr = flt_yr - styr;
          if(yr >= 0 && yr < nyrs_hind) {
            diet_prop(r_idx, k_idx, r_age, k_age, yr) = obs_val;
          }
        }
        else if(flt_yr == 0) {
          for(int yr = 0; yr < nyrs; yr++) {
            diet_prop(r_idx, k_idx, r_age, k_age, yr) = obs_val;
          }
        }
      }
    }
  }
}


/**
 * @brief Calculates the residual stomach content representing "other food" not explicitly modeled.
 *
 * Dividing the unattributed diet share by other_food (mt) puts it on the same per-biomass
 * scale as stom_div_bio, which calculate_msvpa_suitability() sums over prey and adds it to.
 */
template <class Type>
void calculate_other_food_diet_prop(
    const int& nspp,
    const int& nyrs,
    const vector<int>& nsex,
    const vector<int>& nages,
    const vector<Type>& other_food,
    array<Type>& diet_prop,
    array<Type>& other_food_diet_prop // Modified by reference
) {
  other_food_diet_prop.setZero();
  for(int yr = 0; yr < nyrs; yr++) {
    for(int rsp = 0; rsp < nspp; rsp++) {
      for(int r_sex = 0; r_sex < nsex(rsp); r_sex++) {
        int r_idx = rsp + (nspp * r_sex);
        for(int r_age = 0; r_age < nages(rsp); r_age++) {

          other_food_diet_prop(rsp, r_sex, r_age, yr) = Type(1.0);

          for(int ksp = 0; ksp < nspp; ksp++) {
            for(int k_sex = 0; k_sex < nsex(ksp); k_sex++) {
              int k_idx = ksp + (nspp * k_sex);
              for(int k_age = 0; k_age < nages(ksp); k_age++) {
                other_food_diet_prop(rsp, r_sex, r_age, yr) -= diet_prop(r_idx, k_idx, r_age, k_age, yr);
              }
            }
          }

          if(other_food(rsp) > 0) {
            other_food_diet_prop(rsp, r_sex, r_age, yr) /= other_food(rsp);
          } else {
            other_food_diet_prop(rsp, r_sex, r_age, yr) = 0.0;
          }
        }
      }
    }
  }
}


/**
 * @brief Predicts the expected stomach content proportions given the modeled population state.
 *
 * @details Maps `diet_prop_hat` back onto the rows of the empirical `diet_obs`, which are
 * often aggregated over ages or sexes. The four diet-data types, as the ages encode them:
 * - **Type 1:** Specific predator age and specific prey age (weighted across sexes if combined).
 * - **Type 2:** Specific predator age, summed across all prey ages.
 * - **Type 3:** Unweighted mean across predator ages, summed across prey ages (r_age > -500).
 * - **Type 4:** As Type 3, but weighted by predator numbers-at-age (r_age < -500).
 *
 * @param nspp Number of modeled species.
 * @param nyrs_hind Number of years in the hindcast period.
 * @param styr The starting year of the model.
 * @param suit_styr Start year index (per predator species) for the suitability averaging period.
 * @param suit_endyr End year index (per predator species) for the suitability averaging period.
 * @param minage Vector containing the minimum modeled age for each species.
 * @param nsex Vector containing the number of sexes modeled for each species.
 * @param nages Vector containing the number of age bins for each species.
 * @param diet_obs Matrix of empirical diet observations (used here for row-matching).
 * @param diet_ctl Control matrix mapping each observation row to specific aggregation levels.
 * @param avgN_at_age 4D array of average numbers-at-age (thousands of fish), used for weighting.
 * @param diet_prop_hat 5D array of the model's predicted diet proportions for specific ages/sexes.
 * @param diet_hat [Output] Matrix of predicted diet proportions formatted to match `diet_obs`. Modified by reference.
 */
template <class Type>
void predict_stomach_content(
    const int& nspp,
    const int& nyrs_hind,
    const int& styr,
    const vector<int>& suit_styr,
    const vector<int>& suit_endyr,
    const vector<int>& minage,
    const vector<int>& nsex,
    const vector<int>& nages,
    matrix<Type>& diet_obs,
    matrix<int>& diet_ctl,
    array<Type>& avgN_at_age,
    array<Type>& diet_prop_hat,
    matrix<Type>& diet_hat // Modified by reference
) {
  for(int stom_ind = 0; stom_ind < diet_obs.rows(); stom_ind++) {
    int rsp = diet_ctl(stom_ind, 0) - 1;
    int ksp = diet_ctl(stom_ind, 1) - 1;
    int r_sex = diet_ctl(stom_ind, 2);
    int k_sex = diet_ctl(stom_ind, 3);
    int r_age = diet_ctl(stom_ind, 4) - minage(rsp);
    int k_age = diet_ctl(stom_ind, 5) - minage(ksp);
    int flt_yr = diet_ctl(stom_ind, 6);

    // Sex looping bounds (0 = combined, 1 = female, 2 = male), predator then prey.
    int r_start = 0;
    int r_end = 0;
    if(nsex(rsp) == 2) {
      if(r_sex == 0) { r_end = 1; }
      else { r_start = r_sex - 1; r_end = r_sex - 1; }
    }
    int num_r_sexes = r_end - r_start + 1; // Used for averaging in TYPE 3

    int k_start = 0;
    int k_end = 0;
    if(nsex(ksp) == 2) {
      if(k_sex == 0) { k_end = 1; }
      else { k_start = k_sex - 1; k_end = k_sex - 1; }
    }

    // flt_yr > 0 is that year, 0 the suitability window, negative the year |flt_yr|.
    int yr_start = 0, yr_end = 0;
    if(flt_yr > 0) {
      yr_start = flt_yr - styr;
      yr_end = flt_yr - styr;
    } else if (flt_yr == 0) {
      yr_start = suit_styr(rsp);
      yr_end = suit_endyr(rsp);
    } else {
      yr_start = -flt_yr - styr;
      yr_end = -flt_yr - styr;
    }

    // Only evaluate if within hindcast limits
    if(yr_start < nyrs_hind && yr_end < nyrs_hind) {
      Type yr_accum = 0;
      int n_yrs_evaluated = yr_end - yr_start + 1;

      for(int yr = yr_start; yr <= yr_end; yr++) {
        Type type_accum = 0;

        // TYPE 1: Specific pred age, specific prey age
        if((k_age >= 0) & (r_age >= 0)) {
          Type weighted_sum = 0;
          Type total_numbers = 0;
          for(int r_s = r_start; r_s <= r_end; r_s++) {
            Type pred_numbers = avgN_at_age(rsp, r_s, r_age, yr);
            total_numbers += pred_numbers;
            for(int k_s = k_start; k_s <= k_end; k_s++) {
              weighted_sum += diet_prop_hat(rsp + (nspp * r_s), ksp + (nspp * k_s), r_age, k_age, yr) * pred_numbers;
            }
          }
          if(total_numbers > 0) type_accum += weighted_sum / total_numbers;
        }

        // TYPE 2: Specific pred age, sum across prey ages
        else if((k_age < 0) & (r_age >= 0)) {
          Type weighted_sum = 0;
          Type total_numbers = 0;
          for(int r_s = r_start; r_s <= r_end; r_s++) {
            Type pred_numbers = avgN_at_age(rsp, r_s, r_age, yr);
            total_numbers += pred_numbers;
            for(int k_s = k_start; k_s <= k_end; k_s++) {
              for(int k_age_tmp = 0; k_age_tmp < nages(ksp); k_age_tmp++) {
                weighted_sum += diet_prop_hat(rsp + (nspp * r_s), ksp + (nspp * k_s), r_age, k_age_tmp, yr) * pred_numbers;
              }
            }
          }
          if(total_numbers > 0) type_accum += weighted_sum / total_numbers;
        }

        // TYPE 3: Mean across pred ages, sum across prey ages (Unweighted)
        else if((k_age < 0) & (r_age < 0) & (r_age > -500)) {
          Type unweighted_sum = 0;
          for(int r_s = r_start; r_s <= r_end; r_s++) {
            for(int r_age_tmp = 0; r_age_tmp < nages(rsp); r_age_tmp++) {
              for(int k_s = k_start; k_s <= k_end; k_s++) {
                for(int k_age_tmp = 0; k_age_tmp < nages(ksp); k_age_tmp++) {
                  unweighted_sum += diet_prop_hat(rsp + (nspp * r_s), ksp + (nspp * k_s), r_age_tmp, k_age_tmp, yr);
                }
              }
            }
          }
          type_accum += unweighted_sum / (num_r_sexes * nages(rsp));
        }

        // TYPE 4: Weighted mean across pred ages, sum across prey ages
        else if((k_age < 0) & (r_age < -500)) {
          Type weighted_sum = 0;
          Type total_numbers = 0;
          for(int r_s = r_start; r_s <= r_end; r_s++) {
            for(int r_age_tmp = 0; r_age_tmp < nages(rsp); r_age_tmp++) {
              Type pred_numbers = avgN_at_age(rsp, r_s, r_age_tmp, yr);
              total_numbers += pred_numbers;
              for(int k_s = k_start; k_s <= k_end; k_s++) {
                for(int k_age_tmp = 0; k_age_tmp < nages(ksp); k_age_tmp++) {
                  weighted_sum += diet_prop_hat(rsp + (nspp * r_s), ksp + (nspp * k_s), r_age_tmp, k_age_tmp, yr) * pred_numbers;
                }
              }
            }
          }
          if(total_numbers > 0) type_accum += weighted_sum / total_numbers;
        }
        yr_accum += type_accum;
      }
      diet_hat(stom_ind, 1) = yr_accum / n_yrs_evaluated;
    }
  }
}

#endif
