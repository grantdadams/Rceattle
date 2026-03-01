#ifndef BIOENERGETICS_HPP
#define BIOENERGETICS_HPP

/**
 * @brief Calculates the temperature-dependent function of consumption (fT) for bioenergetics.
 *
 * Updates the fT matrix in place. Supports multiple bioenergetics models:
 * 1 = Exponential function (Stewart et al. 1983)
 * 2 = Warm-water species temperature dependence (Kitchell et al. 1977)
 * 3 = Cool/cold-water species temperature dependence (Thornton and Lessem 1979)
 * 4 = Constant (fT = 1.0)
 *
 * @param fT         Matrix of temperature functions to be updated [nspp, nyrs].
 * @param nspp       Number of species.
 * @param nyrs       Number of years.
 * @param Ceq        Vector specifying which consumption equation to use for each species.
 * @param Cindex     Vector specifying which environmental index drives each species.
 * @param Qc         Bioenergetics parameter (QC value / intercept).
 * @param Tcm        Bioenergetics parameter (Thermal max).
 * @param Tco        Bioenergetics parameter (Thermal optimum).
 * @param Tcl        Bioenergetics parameter (Thermal limit for Thornton/Lessem).
 * @param CK1        Bioenergetics parameter (Ascending limit for Thornton/Lessem).
 * @param CK4        Bioenergetics parameter (Descending limit for Thornton/Lessem).
 * @param env_index  Matrix of environmental predictors (e.g., bottom temperature) [nyrs, n_indices].
 */
template <class Type>
void calculate_temperature_function(matrix<Type> &fT,
                                    int nspp,
                                    int nyrs,
                                    const vector<int> &Ceq,
                                    const vector<int> &Cindex,
                                    const vector<Type> &Qc,
                                    const vector<Type> &Tcm,
                                    const vector<Type> &Tco,
                                    const vector<Type> &Tcl,
                                    const vector<Type> &CK1,
                                    const vector<Type> &CK4,
                                    const matrix<Type> &env_index) {

  // Pre-allocate temporary variables used in Kitchell and Thornton/Lessem equations
  // Hoisted out of the loops to optimize memory allocation on the AD tape
  Type Yc = 0.0;
  Type Zc = 0.0;
  Type Vc = 0.0;
  Type Xc = 0.0;
  Type G2 = 0.0;
  Type L2 = 0.0;
  Type G1 = 0.0;
  Type L1 = 0.0;
  Type Ka = 0.0;
  Type Kb = 0.0;

  for(int sp = 0; sp < nspp; sp++) {
    for(int yr = 0; yr < nyrs; yr++) {

      switch(Ceq(sp)) {
      case 1: // Exponential function from Stewart et al. 1983
        fT(sp, yr) = exp(Qc(sp) * env_index(yr, Cindex(sp)));
        break;

      case 2: // Temperature dependence for warm-water species from Kitchell et al. 1977
        Yc = log(Qc(sp)) * (Tcm(sp) - Tco(sp) + 2.0);
        Zc = log(Qc(sp)) * (Tcm(sp) - Tco(sp));
        Vc = (Tcm(sp) - env_index(yr, Cindex(sp))) / (Tcm(sp) - Tco(sp));
        Xc = pow(Zc, 2) * pow((1.0 + pow((1.0 + 40.0 / Yc), 0.5)), 2) / 400.0;
        fT(sp, yr) = pow(Vc, Xc) * exp(Xc * (1.0 - Vc));
        break;

      case 3: // Temperature dependence for cool and cold-water species from Thornton and Lessem 1979
        G2 = (1.0 / (Tcl(sp) - Tcm(sp))) * log((0.98 * (1.0 - CK4(sp))) / (CK4(sp) * 0.02));
        L2 = exp(G2 * (Tcl(sp) - env_index(yr, Cindex(sp))));
        Kb = (CK4(sp) * L2) / (1.0 + CK4(sp) * (L2 - 1.0));

        G1 = (1.0 / (Tco(sp) - Qc(sp))) * log((0.98 * (1.0 - CK1(sp))) / (CK1(sp) * 0.02));
        L1 = exp(G1 * (env_index(yr, Cindex(sp)) - Qc(sp)));
        Ka = (CK1(sp) * L1) / (1.0 + CK1(sp) * (L1 - 1.0));

        fT(sp, yr) = Ka * Kb;
        break;

      case 4: // Constant
        fT(sp, yr) = 1.0;
        break;
      }
    }
  }
}


/**
 * @brief Calculates the historic daily and annual ration for predators.
 *
 * This function calculates the amount of food consumed by a predator of a specific
 * age and sex in a given year. It updates both the `consumption_at_age` and `ration`
 * arrays in place. It supports either full bioenergetics equations (Ceq < 4) or
 * direct empirical input of annual ration (Ceq == 4).
 *
 * @param consumption_at_age Array updated with individual consumption (kg/yr) [nspp, max_sex, max_age, nyrs].
 * @param nspp               Total number of species in the model.
 * @param nyrs               Total number of years (hindcast + forecast).
 * @param nyrs_hind          Number of years in the hindcast period.
 * @param nsex               Vector containing the number of sexes modeled for each species.
 * @param nages              Vector containing the number of age classes for each species.
 * @param Ceq                Vector specifying the consumption equation type for each species.
 * @param CA                 Vector of intercepts for the allometric mass function (Cmax = CA * W^(1+CB)).
 * @param CB                 Vector of slopes for the allometric mass function.
 * @param fday               Vector of the number of foraging days per year for each species.
 * @param Pvalue             Vector of the proportion of maximum consumption (scales Cmax).
 * @param fT                 Matrix of the temperature-dependence function [nspp, nyrs].
 * @param ration_data        Array of empirical relative-foraging rates or input ration limits.
 * @param weight_hat         Array of estimated weight-at-age [n_wt_indices, max_sex, max_age, nyrs].
 */
template <class Type>
void calculate_ration(array<Type> &consumption_at_age,
                      int nspp,
                      int nyrs,
                      int nyrs_hind,
                      const vector<int> &nsex,
                      const vector<int> &nages,
                      const vector<int> &Ceq,
                      const vector<Type> &CA,
                      const vector<Type> &CB,
                      const vector<Type> &fday,
                      const vector<Type> &Pvalue,
                      matrix<Type> &fT, // "const", but doesnt work in TMB for matrices and arrays
                      array<Type> &ration_data,
                      array<Type> &weight_hat) {

  for(int sp = 0; sp < nspp; sp++) {
    for(int sex = 0; sex < nsex(sp); sex++) {
      for(int age = 0; age < nages(sp); age++) {
        for(int yr = 0; yr < nyrs; yr++) {

          // Use the last hindcast year's ration data for forecast years
          int yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);

          if(Ceq(sp) < 4){
            int wt_idx_pop = 2 * sp;

            // Step 1: Calculate max consumption (C_max * f(T) * weight * fday) in grams/predator/year
            consumption_at_age(sp, sex, age, yr) = CA(sp) * pow(weight_hat(wt_idx_pop, sex, age, yr) * Type(1000.0), 1.0 + CB(sp)) * fT(sp, yr) * fday(sp);

            // Step 2: Scale by Pvalue and ration_data, then convert to kg/yr
            consumption_at_age(sp, sex, age, yr) = consumption_at_age(sp, sex, age, yr) * Pvalue(sp) * ration_data(sp, sex, age, yr_ind) / 1000.0;
          }

          // If Ceq == 4, bypass bioenergetics and use empirical inputs directly
          if(Ceq(sp) == 4){
            consumption_at_age(sp, sex, age, yr) = Pvalue(sp) * ration_data(sp, sex, age, yr_ind);
          }
        }
      }
    }
  }
}


#endif
