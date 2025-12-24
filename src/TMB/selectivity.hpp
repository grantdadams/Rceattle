/**
 * @brief Calculates fishery and survey selectivity across fleets, sexes, ages, and years.
 * * Supports empirical data lookup and multiple parametric/non-parametric functional forms.
 * * @section functional_forms Selectivity Models:
 * - Case 0: Empirical (Input data)
 * - Case 1: Logistic (2-parameter: slope and inf)
 * - Case 2: Non-parametric (Age-specific coefficients with smoothing)
 * - Case 3: Double Logistic (Dorn and Methot 1990)
 * - Case 4: Descending Logistic
 * - Case 5: Non-parametric (Hake/Taylor style cumulative coefficients)
 * * @param nspp Number of species.
 * @param n_flt Number of fleets (fisheries and surveys).
 * @param nyrs Total years (including hindcast and projection).
 * @param nyrs_hind Number of hindcast years.
 * @param styr Start year (for indexing).
 * @param nsex Vector of number of sexes per species.
 * @param nages Vector of max age per species.
 * @param flt_spp Vector mapping fleet to species index.
 * @param flt_sel_type Functional form index per fleet.
 * @param flt_sel_age Minimum age selected.
 * @param nselages Maximum age with estimated coefficients.
 * @param flt_sel_maxage Age used for normalization.
 * @param flt_sel_maxage_upper Upper age range for mean normalization.
 * @param emp_sel_obs Matrix of empirical selectivity values.
 * @param emp_sel_ctl Control matrix for empirical selectivity indexing.
 * @param ln_sel_slp Matrix of log-scale logistic slope parameters.
 * @param ln_sel_slp_dev 4D array of yearly slope deviations.
 * @param sel_inf Matrix of age-at-inflection parameters.
 * @param sel_inf_dev 4D array of yearly inflection deviations.
 * @param sel_coff 3D array of age-specific coefficients.
 * @param sel_coff_dev 4D array of yearly coefficient deviations.
 * @param non_par_sel [Output] 4D array of unnormalized non-parametric selectivity (fleet, sex, age, year) to be filled.
 * @param sel [Output] 4D array (fleet, sex, age, year) to be filled.
 */
template<class Type>
void calculate_age_selectivity(
    int nspp, int n_flt,
    int nyrs, int nyrs_hind, int styr,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_sel_age,
    const vector<int>&  flt_nselages,
    const vector<int>&  flt_sel_maxage,
    const vector<int>&  flt_sel_maxage_upper,
    matrix<Type> emp_sel_obs,
    const matrix<int>& emp_sel_ctl,
    array<Type> ln_sel_slp,
    array<Type> ln_sel_slp_dev,
    array<Type> sel_inf,
    array<Type> sel_inf_dev,
    array<Type> sel_coff,
    array<Type> sel_coff_dev,
    array<Type> &avg_sel,
    array<Type> &non_par_sel,
    array<Type> &sel
) {
  sel.setZero();
  Type max_sel = 0;
  Type avgsel_tmp = 0; // Temporary object for average selectivity across all ages

  // --- 1. EMPIRICAL SELECTIVITY ---
  for (int i = 0; i < emp_sel_obs.rows(); i++) {
    int flt = emp_sel_ctl(i, 0) - 1;
    int sp  = emp_sel_ctl(i, 1) - 1;
    int s_idx = emp_sel_ctl(i, 2);
    int yr = emp_sel_ctl(i, 3) - styr;

    if (flt_sel_type(flt) == 0 && yr < nyrs_hind) {
      // Map sex index (0=both, 1=female, 2=male)
      int s_start = (s_idx == 0) ? 0 : s_idx - 1;
      int s_end   = (s_idx == 0 && nsex(sp) == 2) ? 1 : s_start;

      for (int sex = s_start; sex <= s_end; sex++) {
        for (int age = 0; age < nages(sp); age++) {
          if (!isNA(emp_sel_obs(i, age))) {
            sel(flt, sex, age, yr) = emp_sel_obs(i, age);
          }
        }
      }
    }
  }

  // --- 2. ESTIMATED SELECTIVITY ---
  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    int sel_type = flt_sel_type(flt);
    int nselages = flt_nselages(flt);
    if (sel_type == 0) continue;

    for (int yr = 0; yr < nyrs_hind; yr++) {
      for (int sex = 0; sex < nsex(sp); sex++) {

        switch (sel_type) {
        case 1: // Logistic
          for (int age = 0; age < nages(sp); age++) {
            Type slope = exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr));
            Type inf   = sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr);
            sel(flt, sex, age, yr) = 1.0 / (1.0 + exp(-slope * (Type(age + 1) - inf)));
          }
          break;

        case 2: // Non-parametric (Ianelli style)
          // {
          //   Type avgsel_tmp = 0;
          //   for (int age = 0; age < nages(sp); age++) {
          //     int a_idx = (age < nselages) ? age : nselages - 1;
          //     Type val = exp(sel_coff(flt, sex, a_idx) + sel_coff_dev(flt, sex, a_idx, yr));
          //     sel(flt, sex, age, yr) = val;
          //     avgsel_tmp += val;
          //   }
          //   avgsel_tmp /= Type(nages(sp));
          //   for (int age = 0; age < nages(sp); age++) sel(flt, sex, age, yr) /= avgsel_tmp;
          // }

          for(int age = 0; age < nselages; age++) {
            non_par_sel(flt, sex, age, yr) = sel_coff(flt, sex, age) + sel_coff_dev(flt, sex, age, yr);
            avg_sel(flt, sex, yr) +=  exp(sel_coff(flt, sex, age) + sel_coff_dev(flt, sex, age, yr));
          }
          //  Average selectivity (up to nselages)
          avg_sel(flt, sex, yr) = log(avg_sel(flt, sex, yr) / nselages);

          // Fill in plus group
          for(int age = nselages; age < nages(sp); age++) {
            non_par_sel(flt, sex, age, yr) = non_par_sel(flt, sex, nselages - 1, yr);
          }

          // Average selectivity (ALL ages)

          avgsel_tmp = 0;
          for(int age = 0; age < nages(sp); age++) {
            avgsel_tmp += exp(non_par_sel(flt, sex, age, yr));
          }
          avgsel_tmp = log(avgsel_tmp / nages(sp));

          // Standardize selectivity
          for(int age = 0; age < nages(sp); age++) {
            non_par_sel(flt, sex, age, yr) -= avgsel_tmp;
            non_par_sel(flt, sex, age, yr) = exp(non_par_sel(flt, sex, age, yr));
            sel(flt, sex, age, yr) = non_par_sel(flt, sex, age, yr);
          }
          break;

        case 3: // Double Logistic
          for (int age = 0; age < nages(sp); age++) {
            Type slp1 = exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr));
            Type inf1 = sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr);
            Type slp2 = exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr));
            Type inf2 = sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr);
            sel(flt, sex, age, yr) = (1.0 / (1.0 + exp(-slp1 * (Type(age + 1) - inf1)))) * (1.0 - (1.0 / (1.0 + exp(-slp2 * (Type(age + 1) - inf2)))));
          }
          break;

        case 4: // Descending Logistic
          for (int age = 0; age < nages(sp); age++) {
            Type slp2 = exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr));
            Type inf2 = sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr);
            sel(flt, sex, age, yr) = (1.0 - (1.0 / (1.0 + exp(-slp2 * (Type(age + 1) - inf2)))));
          }
          break;

        case 5: // Hake/Taylor Non-parametric
          {
            Type cum_sum = 0;
            for (int age = flt_sel_age(flt); age < nages(sp); age++) {
              if (age < nselages) {
                cum_sum += sel_coff(flt, sex, age) + sel_coff_dev(flt, sex, age, yr);
              }
              sel(flt, sex, age, yr) = cum_sum;
            }
            // Normalize inside year/sex block
            max_sel = -1e10; // Use log-space for normalization
            if (flt_sel_maxage(flt) >= 0 && flt_sel_maxage_upper(flt) < 0)
              max_sel = sel(flt, sex, flt_sel_maxage(flt), yr);
            else {
              for(int age=0; age<nages(sp); age++) if(sel(flt,sex,age,yr) > max_sel) max_sel = sel(flt,sex,age,yr);
            }
            for (int age = 0; age < nages(sp); age++) sel(flt, sex, age, yr) = exp(sel(flt, sex, age, yr) - max_sel);
          }
          break;
        }
      }
    }
  }
}

/**
 * @brief Converts length-based selectivity to age-based selectivity and stores it in a 4D array.
 * * @details This function maps selectivity from length bins to age classes by
 * integrating over a growth transition matrix (Prob(Length | Age)). The resulting
 * age-based selectivity is stored directly into the provided 'sel' array for the
 * specified fleet, sex, and year using pass by reference.
 *
 * @tparam Type The numeric type (e.g., double, TMB::adouble).
 * @param flt Index of the fleet/fishery.
 * @param sp Index of the species.
 * @param sex Index of the sex.
 * @param yr Index of the year.
 * @param wtind Index for the weight-at-age/growth transition matrix type.
 * @param nages Reference to vector containing number of ages per species.
 * @param nlengths Reference to vector containing number of length bins per species.
 * @param sel_length Reference to 4D array of selectivity-at-length [flt, sex, ln, yr].
 * @param growth_matrix Reference to 5D transition matrix [wtind, sex, age, ln, yr].
 * @param sel Reference to 4D output array [flt, sex, age, yr] to be updated in-place.
 */
template <class Type>
void convert_length_selectivity(
    int flt,
    int sp,
    int sex,
    int yr,
    int wtind,
    const vector<int>& nages,
    const vector<int>& nlengths,
    array<Type>& sel_length,
    array<Type>& growth_matrix,
    array<Type>& sel
) {
  // Iterate through ages for the specific species
  for(int age = 0; age < nages(sp); age++) {

    // Initialize the specific cell in the 4D array to zero before accumulating
    sel(flt, sex, age, yr) = 0.0;

    for(int ln = 0; ln < nlengths(sp); ln++) {
      // sel_at_age = sum over lengths of (Prob(Length|Age) * sel_at_length)
      sel(flt, sex, age, yr) += growth_matrix(wtind, sex, age, ln, yr) * sel_length(flt, sex, ln, yr);
    }
  }
}

/**
 * @brief Calculates fishery and survey length-based selectivity across fleets, sexes, length, and years.
 * * Supports multiple parametric/non-parametric functional forms.
 * * @section functional_forms Selectivity Models:
 * - Case 1: Logistic (2-parameter: slope and inf)
 * - Case 2: Non-parametric (Age-specific coefficients with smoothing)
 * - Case 3: Double Logistic (Dorn and Methot 1990)
 * - Case 4: Descending Logistic
 * - Case 5: Non-parametric (Hake/Taylor style cumulative coefficients)
 * * @param nspp Number of species.
 * @param n_flt Number of fleets (fisheries and surveys).
 * @param nsex Vector of number of sexes per species.
 * @param nlengths Vector of max length per species.
 * @param nyrs Total years (including hindcast and projection).
 * @param nyrs_hind Number of hindcast years.
 * @param styr Start year (for indexing).
 * @param flt_spp Vector mapping fleet to species index.
 * @param flt_sel_type Functional form index per fleet.
 * @param flt_sel_age Minimum age selected.
 * @param nselages Maximum length with estimated coefficients. //FIXME
 * @param flt_sel_maxage Age used for normalization.
 * @param flt_sel_maxage_upper Upper age range for mean normalization.
 * @param emp_sel_obs Matrix of empirical selectivity values.
 * @param emp_sel_ctl Control matrix for empirical selectivity indexing.
 * @param ln_sel_slp Matrix of log-scale logistic slope parameters.
 * @param ln_sel_slp_dev 4D array of yearly slope deviations.
 * @param sel_inf Matrix of length-at-inflection parameters.
 * @param sel_inf_dev 4D array of yearly inflection deviations.
 * @param sel_coff 3D array of length-specific coefficients.
 * @param sel_coff_dev 4D array of yearly coefficient deviations.
 * @param non_par_sel [Output] 4D array of unnormalized non-parametric selectivity (fleet, sex, age, year) to be filled.
 * @param sel [Output] 4D array (fleet, sex, age, year) to be filled.
 * @param growth_matrix growth matrix defining the p(age|length) calculated in "calculate_growth" in caal.hpp
 *
 */
template<class Type>
void calculate_length_selectivity(
    int& nspp,
    int& n_flt,
    int& nyrs,
    int& nyrs_hind,
    int& styr,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_sel_age,
    const vector<int>&  flt_nselages,
    const vector<int>&  flt_sel_maxage,
    const vector<int>&  flt_sel_maxage_upper,
    array<Type>& ln_sel_slp,
    array<Type>& ln_sel_slp_dev,
    array<Type>& sel_inf,
    array<Type>& sel_inf_dev,
    array<Type>& sel_coff,
    array<Type>& sel_coff_dev,
    array<Type> &avg_sel,
    array<Type> &non_par_sel, // Pass by reference
    array<Type> &sel, // Pass by reference
    array<Type>& growth_matrix
) {

  Type max_sel = 0;
  Type avgsel_tmp = 0; // Temporary object for average selectivity across all ages

  array<Type>   length_sel(n_flt, imax(nsex), imax(nlengths), nyrs);

  length_sel.setZero();


  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    int sel_type = flt_sel_type(flt);
    int nselages = flt_nselages(flt);
    if (sel_type == 0) continue;

    for (int yr = 0; yr < nyrs_hind; yr++) {
      for (int sex = 0; sex < nsex(sp); sex++) {

        switch (sel_type) {
        case 6: // Logistic ----
          for (int ln = 0; ln < nlengths(sp); ln++) {
            Type slope = exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr));
            Type inf   = sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr);
            length_sel(flt, sex, ln, yr) = 1.0 / (1.0 + exp(-slope * (Type(ln + 1) - inf)));
          }
          break;

        case 7: // Non-parametric (Ianelli style) ----
          // {
          //   Type avgsel_tmp = 0;
          //   for (int ln = 0; ln < nlengths(sp); ln++) {
          //     int a_idx = (ln < nselages) ? ln : nselages - 1;
          //     Type val = exp(sel_coff(flt, sex, a_idx) + sel_coff_dev(flt, sex, a_idx, yr));
          //     length_sel(flt, sex, ln, yr) = val;
          //     avgsel_tmp += val;
          //   }
          //   avgsel_tmp /= Type(nlengths(sp));
          //   for (int ln = 0; ln < nlengths(sp); ln++) length_sel(flt, sex, ln, yr) /= avgsel_tmp;
          // }
          for(int ln = 0; ln < nselages; ln++) {
            non_par_sel(flt, sex, ln, yr) = sel_coff(flt, sex, ln) + sel_coff_dev(flt, sex, ln, yr);
            avg_sel(flt, sex, yr) +=  exp(sel_coff(flt, sex, ln) + sel_coff_dev(flt, sex, ln, yr));
          }
          //  Average selectivity (up to nselages)
          avg_sel(flt, sex, yr) = log(avg_sel(flt, sex, yr) / nselages);

          // Fill in plus group
          for(int ln = nselages; ln < nlengths(sp); ln++) {
            non_par_sel(flt, sex, ln, yr) = non_par_sel(flt, sex, nselages - 1, yr);
          }

          // Average selectivity (ALL ages)
          avgsel_tmp = 0; // Temporary object for average selectivity across all ages
          for(int ln = 0; ln < nlengths(sp); ln++) {
            avgsel_tmp += exp(non_par_sel(flt, sex, ln, yr));
          }
          avgsel_tmp = log(avgsel_tmp / nlengths(sp));

          // Standardize selectivity
          for(int ln = 0; ln < nlengths(sp); ln++) {
            non_par_sel(flt, sex, ln, yr) -= avgsel_tmp;
            non_par_sel(flt, sex, ln, yr) = exp(non_par_sel(flt, sex, ln, yr));
            sel(flt, sex, ln, yr) = non_par_sel(flt, sex, ln, yr);
          }
          break;

        case 8: // Double Logistic ----
          for (int ln = 0; ln < nlengths(sp); ln++) {
            Type slp1 = exp(ln_sel_slp(0, flt, sex) + ln_sel_slp_dev(0, flt, sex, yr));
            Type inf1 = sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr);
            Type slp2 = exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr));
            Type inf2 = sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr);
            length_sel(flt, sex, ln, yr) = (1.0 / (1.0 + exp(-slp1 * (Type(ln + 1) - inf1)))) * (1.0 - (1.0 / (1.0 + exp(-slp2 * (Type(ln + 1) - inf2)))));
          }
          break;

        case 9: // Descending Logistic ----
          for (int ln = 0; ln < nlengths(sp); ln++) {
            Type slp2 = exp(ln_sel_slp(1, flt, sex) + ln_sel_slp_dev(1, flt, sex, yr));
            Type inf2 = sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr);
            length_sel(flt, sex, ln, yr) = (1.0 - (1.0 / (1.0 + exp(-slp2 * (Type(ln + 1) - inf2)))));
          }
          break;

        case 10: // Hake/Taylor Non-parametric ----
          {
            Type cum_sum = 0;
            for (int ln = flt_sel_age(flt); ln < nlengths(sp); ln++) {
              if (ln < nselages) {
                cum_sum += sel_coff(flt, sex, ln) + sel_coff_dev(flt, sex, ln, yr);
              }
              length_sel(flt, sex, ln, yr) = cum_sum;
            }
            // Normalize inside year/sex block
            max_sel = -1e10; // Use log-space for normalization
            if (flt_sel_maxage(flt) >= 0 && flt_sel_maxage_upper(flt) < 0)
              max_sel = length_sel(flt, sex, flt_sel_maxage(flt), yr);
            else {
              for(int ln=0; ln<nlengths(sp); ln++) if(length_sel(flt,sex,ln,yr) > max_sel) max_sel = length_sel(flt,sex,ln,yr);
            }
            for (int ln = 0; ln < nlengths(sp); ln++) length_sel(flt, sex, ln, yr) = exp(length_sel(flt, sex, ln, yr) - max_sel);
          }
          break;
        }


        // Convert length-based selectivity to age-based
        int wtind = nspp * 2 + flt;
        convert_length_selectivity(
          flt,
          sp,
          sex,
          yr,
          wtind,
          nages,
          nlengths,
          length_sel,
          growth_matrix,
          sel
        );
      } // End sex
    } // End yr
  } // End flt
}


/**
 * @brief Performs final normalization and projection of fishery selectivity across all fleets.
 * * @details This function iterates through all fleets to process selectivity data.
 * It handles:
 * 1. Zeroing out selectivity for ages below flt_sel_age.
 * 2. Normalizing selectivity values based on a single age (flt_sel_maxage >= 0),
 * an age range average (flt_sel_maxage_upper >= 0), or the maximum value
 * across all ages/sexes (flt_sel_maxage between -1 and -499).
 * 3. Projecting the final year of hindcast selectivity into future projection years.
 *
 * @param n_flt Total number of fleets to iterate through.
 * @param nyrs_hind Total number of years in the hindcast period.
 * @param nyrs Total number of years including projection.
 * @param flt_spp Function/Array mapping fleet index to species index.
 * @param flt_sel_type Function/Array mapping fleet index to selectivity type.
 * @param flt_sel_age Array/function returning the minimum age of selection for a fleet.
 * @param nages Array/function returning the number of ages for a species.
 * @param nsex Array/function returning the number of sexes for a species.
 * @param flt_sel_maxage Array/function returning the normalization age or control flag.
 * @param flt_sel_maxage_upper Array/function returning the upper bound for age-range normalization.
 * @param sel 4D container (fleet, sex, age, year) modified in-place.
 */

template<class Type>
void normalize_and_project_selectivity(
    int n_flt,
    int nyrs_hind,
    int nyrs,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_sel_age,
    const vector<int>&  nages,
    const vector<int>&  nsex,
    const vector<int>&  flt_sel_maxage,
    const vector<int>&  flt_sel_maxage_upper,
    array<Type> &sel
) {
  Type max_sel; // Local declaration for safety

  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    int sel_type = flt_sel_type(flt);
    // int nselages = flt_nselages(flt); // Note: defined but currently unused in logic

    if(sel_type > 0) {
      // Ages not selected
      for(int yr = 0; yr < nyrs_hind; yr++) {
        for(int age = 0; age < nages(sp); age++){
          for(int sex = 0; sex < nsex(sp); sex++){
            if(age < flt_sel_age(flt)){
              sel(flt, sex, age, yr) = 0;
            }
          }
        }
      }

      // Normalize selectivity
      if(flt_sel_maxage(flt) > -500){

        // 1. Normalize by selectivity by specific age or age-range
        if((flt_sel_maxage(flt) >= 0) && (sel_type < 5)) {
          for(int yr = 0; yr < nyrs_hind; yr++) {
            for(int sex = 0; sex < nsex(sp); sex++){

              // Single-normalization age
              if((flt_sel_maxage(flt) >= 0) && (flt_sel_maxage_upper(flt) < 0)){
                max_sel = 0.001;
                max_sel = max2(max_sel, sel(flt, sex, flt_sel_maxage(flt), yr));
              }

              // Normalize by age range between max lower and max upper
              if((flt_sel_maxage(flt) >= 0) && (flt_sel_maxage_upper(flt) >= 0)){
                max_sel = 0;
                for(int age = flt_sel_maxage(flt); age <= flt_sel_maxage_upper(flt); age++) {
                  max_sel += sel(flt, sex, age, yr)/(flt_sel_maxage_upper(flt) - flt_sel_maxage(flt) + 1);
                }
              }

              // Normalize
              for(int age = 0; age < nages(sp); age++){
                sel(flt, sex, age, yr) /= max_sel;
              }
            }
          }
        }

        // 2. Normalize by max for each fishery and year across ages, and sexes
        if((sel_type < 5) && (flt_sel_maxage(flt) < 0) && (flt_sel_maxage(flt) > -500)) {
          for(int yr = 0; yr < nyrs_hind; yr++) {
            max_sel = 0;
            for(int age = 0; age < nages(sp); age++){
              for(int sex = 0; sex < nsex(sp); sex++){
                // Find max
                if(sel(flt, sex, age, yr) > max_sel){
                  max_sel = sel(flt, sex, age, yr);
                }
              }
            }

            // Normalize selectivity
            for(int age = 0; age < nages(sp); age++){
              for(int sex = 0; sex < nsex(sp); sex++){
                sel(flt, sex, age, yr) /= max_sel;
              }
            }
          }
        }
      }
    }

    // Project outwards assuming final years selectivity
    for(int yr = nyrs_hind; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++){
        for(int sex = 0; sex < nsex(sp); sex++){
          sel(flt, sex, age, yr) = sel(flt, sex, age, nyrs_hind - 1);
        }
      }
    }
  }
}
