#ifndef DIET_CONTENT_HPP
#define DIET_CONTENT_HPP

/**
 * @brief Parses empirical diet observations and maps them into a 5D proportion array.
 * * @details This function takes the raw diet observation matrix (`diet_obs`) and its
 * corresponding control matrix (`diet_ctl`) and maps the observed proportions into
 * the structured 5D `diet_prop` array for MSVPA based suitability.
 * It dynamically handles indexing for 1-sex
 * versus 2-sex models. If the control matrix specifies a year `flt_yr == 0`, the
 * observed proportion is applied as an average across all modeled years. Data
 * rows missing specific predator or prey ages are skipped during this step.
 * * @param nspp Number of modeled species.
 * @param nyrs Total number of years in the model (hindcast + projection).
 * @param nyrs_hind Number of years in the hindcast period.
 * @param styr The starting year of the model.
 * @param minage Vector containing the minimum modeled age for each species.
 * @param nsex Vector containing the number of sexes modeled for each species (1 or 2).
 * @param diet_obs Matrix of empirical diet observations (column 1 contains the proportions).
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
    int r_sex = diet_ctl(stom_ind, 2);               // Index of pred sex
    int k_sex = diet_ctl(stom_ind, 3);               // Index of prey sex
    int r_age = diet_ctl(stom_ind, 4) - minage(rsp); // Index of pred age
    int k_age = diet_ctl(stom_ind, 5) - minage(ksp); // Index of prey age
    int flt_yr = diet_ctl(stom_ind, 6);              // Index of year
    // Diet proportion of prey-at-at in predator-at-age
    // if K_age < 0, data are diet proportion of prey-spp in predator-at-age (summed across prey ages) and are not used for MSVPA based suitability
    // if k_age < 0 & r_age < 0, data are diet proportion of prey-spp in predator-spp (summed across prey ages and averaged across predator ages)
    // and are not used for MSVPA based suitability

    // Only process if specific ages are provided for both predator and prey
    if((k_age < 0) || (r_age < 0)) continue;

    // Determine sex looping bounds for predator (0 = combined, 1 = female, 2 = male)
    int r_start = 0;
    int r_end = 0;
    if(nsex(rsp) == 2) {
      if(r_sex == 0) { r_end = 1; } // Loop over both sexes
      else { r_start = r_sex - 1; r_end = r_sex - 1; } // Specific sex
    }

    // Determine sex looping bounds for prey
    int k_start = 0;
    int k_end = 0;
    if(nsex(ksp) == 2) {
      if(k_sex == 0) { k_end = 1; } // Loop over both sexes
      else { k_start = k_sex - 1; k_end = k_sex - 1; } // Specific sex
    }

    // Value to assign
    Type obs_val = diet_obs(stom_ind, 1);

    // Assign Diet proportion
    for(int j = r_start; j <= r_end; j++) {
      int r_idx = rsp + (nspp * j);

      for(int k = k_start; k <= k_end; k++) {
        int k_idx = ksp + (nspp * k);

        // Specific year data
        if(flt_yr > 0) {
          int yr = flt_yr - styr;
          if(yr >= 0 && yr < nyrs_hind) {
            diet_prop(r_idx, k_idx, r_age, k_age, yr) = obs_val;
          }
        }
        // Average across all years
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

          other_food_diet_prop(rsp, r_sex, r_age, yr) = Type(1.0); // Initialize

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
 * * @details This function maps the model's internal prediction of diet proportions (`diet_prop_hat`)
 * back to the structure of the original empirical data (`diet_obs`). Because empirical
 * stomach data is often aggregated across ages or sexes, this function dynamically
 * calculates the appropriate expected value using abundance-weighted averages (`avgN_at_age`).
 * * It handles four types of data aggregation:
 * - **Type 1:** Specific predator age and specific prey age (weighted across sexes if combined).
 * - **Type 2:** Specific predator age, summed across all prey ages.
 * - **Type 3:** Unweighted arithmetic mean across predator ages, summed across prey ages.
 * - **Type 4:** Abundance-weighted mean across predator ages, summed across prey ages.
 * * @param nspp Number of modeled species.
 * @param nyrs_hind Number of years in the hindcast period.
 * @param styr The starting year of the model.
 * @param suit_styr The starting year used for suitability averaging.
 * @param suit_endyr The ending year used for suitability averaging.
 * @param minage Vector containing the minimum modeled age for each species.
 * @param nsex Vector containing the number of sexes modeled for each species.
 * @param nages Vector containing the maximum number of age bins for each species.
 * @param diet_obs Matrix of empirical diet observations (used here for row-matching).
 * @param diet_ctl Control matrix mapping each observation row to specific aggregation levels.
 * @param avgN_at_age 4D array of average population numbers-at-age (used for weighting).
 * @param diet_prop_hat 5D array of the model's predicted diet proportions for specific ages/sexes.
 * @param diet_hat [Output] Matrix of predicted diet proportions formatted to match `diet_obs`. Modified by reference.
 */
template <class Type>
void predict_stomach_content(
    const int& nspp,
    const int& nyrs_hind,
    const int& styr,
    const int& suit_styr,
    const int& suit_endyr,
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

    // Determine sex looping bounds for predator
    int r_start = 0;
    int r_end = 0;
    if(nsex(rsp) == 2) {
      if(r_sex == 0) { r_end = 1; } // Loop over both sexes
      else { r_start = r_sex - 1; r_end = r_sex - 1; } // Specific sex
    }
    int num_r_sexes = r_end - r_start + 1; // Used for averaging in TYPE 3

    // Determine sex looping bounds for prey
    int k_start = 0;
    int k_end = 0;
    if(nsex(ksp) == 2) {
      if(k_sex == 0) { k_end = 1; } // Loop over both sexes
      else { k_start = k_sex - 1; k_end = k_sex - 1; } // Specific sex
    }

    // Determine evaluation years
    int yr_start = 0, yr_end = 0;
    if(flt_yr > 0) {
      yr_start = flt_yr - styr;
      yr_end = flt_yr - styr;
    } else if (flt_yr == 0) {
      yr_start = suit_styr;
      yr_end = suit_endyr;
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
