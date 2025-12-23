/**
 * @brief Integrated Predator Growth, Size-Transition, and Weight-at-Age Module.
 * * This function handles the core biological calculations for a single species,
 * integrating growth dynamics with size-structured transitions and allometric
 * weight relationships.
 *
 * @section math_models Mathematical Models:
 * 1. Mean Length-at-Age ($L_a$):
 * - Von Bertalanffy (Model 1): $L_a = L_{\infty} (1 - e^{-K(a - a_0)})$
 * - Richards (Model 2): $L_a = L_{\infty} [1 + e^{-K(a - a_0)}]^{-1/m}$
 * - Note: Linear growth is applied for ages < minage to ensure model stability.
 * * 2. Weight-at-Age ($W_a$):
 * Calculated via integration across the length distribution to account for Jensen's Inequality:
 * $W_a = \sum_{ln} P(ln | a) \times (\alpha \cdot L_{mid}^{ \beta })$
 * * @section logic Biological Logic:
 * - **Temporal Resolution**: Incorporates 'fracyr' to allow for within-year (seasonal)
 * growth and differentiability for time-varying parameters.
 * - **Plus-Group (SS Style)**: Oldest age class is adjusted using a weighted mean
 * based on an assumed mortality decay ($\exp(-0.2)$).
 * - **Size Transition**: Converts mean length and SD into a probability density matrix
 * ($P(Length | Age)$) using a cumulative normal distribution (pnorm).
 *
 * @param sp Species index.
 * @param fracyr Fraction of the year (0.0 - 1.0) for within-year growth estimation.
 * @param max_sex/max_age/max_nlengths Global array dimensions for safe memory allocation.
 * @param nyrs Number of years in the simulation/model.
 * @param nsex/nages/nlengths Vectors containing species-specific counts for each dimension.
 * @param lengths Matrix of length bin boundaries (nspp x length_bins).
 * @param growth_parameters 4D array of time-varying growth parameters (K, L1, Linf, m).
 * @param growth_ln_sd 3D array of log-scale standard deviations for length (youngest and oldest ages).
 * @param LW_par 4D array of length-weight relationship parameters ($\alpha$, $\beta$).
 *
 * *@note This version uses pass-by-reference for outputs to optimize memory.
 * @param length_hat [Output] 3D Array to be filled with mean length-at-age.
 * @param growth_matrix [Output] 4D Array to be filled with the growth transition matrix.
 * @param weight_hat [Output] 3D Array to be filled with integrated weight-at-age.
 */
template<class Type>
void calculate_growth(int sp, Type fracyr, int nspp, int nyrs,
                       vector<int> nsex, vector<int> nages, matrix<int> lengths,
                       vector<int> nlengths, vector<int> minage, vector<int> maxage,
                       array<Type> growth_parameters, array<Type> growth_ln_sd,
                       array<Type> weight_length_pars,
                       vector<int> growth_model,
                       array<Type> &length_hat,     // Pass by reference
                       array<Type> &growth_matrix,  // Pass by reference
                       array<Type> &weight_hat      // Pass by reference
) {

  // Initialize output and temporary storage
  array<Type> length_sd(nsex(sp), nages(sp), nyrs); length_sd.setZero();         // SD in length-at-age


  // Calculate mean-length, SD, and growth matrix, for all years:
  // lengths is vector with lengths mm (2, 4, 6, 8, etc)
  Type Fac1, Fac2, Slope, b_len, last_linear, current_age;

  Type Lmin_sp = lengths(sp, 0);
  Type Lmax_sp = lengths(sp, nlengths(sp) - 1);
  Type age_L1 = minage(sp);
  Type age_L1_ceil = maxage(sp);

  for(int sex = 0; sex < nsex(sp); sex++) {
    for(int yr = 0; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++) {

        // Parameters for parametric growth
        Type kappa = growth_parameters(sp, sex, yr, 0);
        Type l1 = growth_parameters(sp, sex, yr, 1);
        Type linf = growth_parameters(sp, sex, yr, 2);
        Type m = growth_parameters(sp, sex, yr, 3);
        current_age = age + 1.0 + fracyr;

        // 1. Calculate Mean Length at Age ---
        switch(growth_model(sp)) {
        case 1: // Von Bertalanffy

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1;

          // Age < minage
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = linf + (l1 - linf) * exp(-kappa * (current_age - age_L1));
            } else { // Yr > 0
              if((age + 1) == age_L1_ceil) {
                last_linear = Lmin_sp + b_len * age_L1;
                length_hat(wtind,  sex, age, yr) = last_linear + (last_linear - linf) * (exp(-kappa * (current_age - age_L1)) - 1.0);
              } else {
                // Lag 1-year parameters
                Type lagkappa = growth_parameters(sp, sex, yr - 1, 0);
                Type laglinf = growth_parameters(sp, sex, yr - 1, 2);
                length_hat(wtind,  sex, age, yr) = length_hat(wtind,  sex, age - 1, yr - 1) + (length_hat(wtind,  sex, age - 1, yr - 1) - laglinf) * (exp(-lagkappa) - 1.0);
              }
            }
          }
          break;

        case 2: // Richards

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1;

          if((current_age) <= age_L1) {
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = pow(pow(l1, m) + (pow(linf, m) - pow(l1, m)) * exp(-kappa * (current_age - age_L1)), 1 / m);
            } else { // Yr > 0
              if((age + 1) == age_L1_ceil) {
                last_linear = Lmin_sp + b_len * age_L1;
                length_hat(wtind,  sex, age, yr) = pow(pow(last_linear, m) + (pow(last_linear, m) - pow(l1, m)) * (exp(-kappa * (current_age - age_L1)) - 1.0), 1 / m);
              } else {
                // Lag 1-year parameters
                Type lagkappa = growth_parameters(sp, sex, yr - 1, 0);
                Type laglinf = growth_parameters(sp, sex, yr - 1, 2);
                Type lagm = growth_parameters(sp, sex, yr - 1, 3);
                length_hat(wtind,  sex, age, yr) = pow(pow(length_hat(wtind,  sex, age - 1, yr - 1), lagm) + (pow(length_hat(wtind,  sex, age - 1, yr - 1), lagm) - pow(laglinf, lagm)) * (exp(-lagkappa) - 1.0), 1 / lagm);
              }
            }
          }
          break;

        case 3: // Non-parametric (Free parameters)
          error("Non-parametric growth not yet implemented");
          // length_hat(wtind,  sex, age, yr) = exp(length_par(sp, sex, age) + length_par_re(sp, sex, age, yr));
          break;

        default:
          error("Invalid 'growth_model");
        } // Growth_model switch

        // 2. Plus-Group Correction (Oldest Age Only) ---
        if(growth_model(sp) < 3 && age == (nages(sp) - 1)) {
          Type temp_n = 0, temp_sum = 0, weight_a = 1.0;
          Type diff = linf - length_hat(wtind,  sex, age, yr);
          for(int a = 0; a < nages(sp); a++) {
            temp_sum += weight_a * (length_hat(wtind,  sex, nages(sp) - 1, yr) + (Type(a) / Type(nages(sp))) * diff);
            temp_n += weight_a;
            weight_a *= exp(-0.2); //FIXME: update mortality?
          }
          length_hat(wtind,  sex, age, yr) = temp_sum / temp_n;
        }

        // 3. Calculate SD (Integrated) ---
        if(growth_model(sp) < 3) {
          if((current_age) < age_L1) {
            length_sd(sex, age, yr) = exp(growth_ln_sd(sp, sex, 0));
          } else if(age == (nages(sp) - 1)) {
            length_sd(sex, age, yr) = exp(growth_ln_sd(sp, sex, 1));
          } else {
            Slope = (exp(growth_ln_sd(sp, sex, 1)) - exp(growth_ln_sd(sp, sex, 0))) / (linf - l1);
            length_sd(sex, age, yr) = exp(growth_ln_sd(sp, sex, 0) + Slope * (length_hat(wtind,  sex, age, yr) - l1));
          }

          // Free parameters
          if(growth_model(sp) == 3) {
            // Slope = (exp(growth_ln_sd(sp, sex, 1)) - exp(growth_ln_sd(sp, sex, 0)))/(length_hat(wtind,  sex, nages(sp)-1, yr) - length_hat(wtind,  sex, 0, yr));
            // length_sd(sex, age, yr) = exp(growth_ln_sd(sp, sex, 0) + Slope * (length_hat(wtind,  sex, age, yr) - length_hat(wtind,  sex, 0, yr));
          }
        }

        // 4. Build Growth Matrix & Weight-at-Age simultaneously ---
        Type expected_weight = 0.0;

        for(int ln = 0; ln < nlengths(sp); ln++) {
          if(ln == 0) {
            Fac1 = (Lmin_sp + lengths(sp, 1) - lengths(sp, 0) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            growth_matrix(wtind,  sex, ln, age, yr) = pnorm(Fac1);
          } else if(ln == (nlengths(sp) - 1)) {
            Fac1 = (Lmax_sp - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            growth_matrix(wtind,  sex, ln, age, yr) = 1.0 - pnorm(Fac1);
          } else {
            Fac1 = (lengths(sp, ln + 1) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            Fac2 = (lengths(sp, ln) - length_hat(wtind,  sex, age, yr)) / length_sd(sex, age, yr);
            growth_matrix(wtind,  sex, ln, age, yr) = pnorm(Fac1) - pnorm(Fac2);
          }

          // Midpoint calculation for Weight-at-Length
          Type lenmid = (lengths(sp, 1) - lengths(sp, 0))/2;

          // Weighted sum for Weight-at-Age
          expected_weight += growth_matrix(wtind, sex, age, ln, yr) * weight_length_pars(sp, 0) * pow(lengths(sp, ln) + lenmid, weight_length_pars(sp, 1));
        }
        weight_hat(wtind, sex, age, yr) = expected_weight;
      } // age
    } // yr
  } // sex
}



// ------------------------------------------------------------------------- //
// TODO                                                           //
// ------------------------------------------------------------------------- //
// ADD MONTH TO FLEET CONTROL
// ADD WEIGHT LENGTH PARAMETERS TO DATA
// add function for weight-at-age from frac-yr
// update wtind across model
// rename "nselages"
// rename "flt_sel_maxage" and "flt_sel_maxage_upper"
// add non_par_sel to selectivity functions



/**
 * @brief Calculates population and fleet-specific weight-at-age.
 * * This function populates the weight_hat array based on either empirical data
 * (growth_model == 0) or estimated growth parameters (growth_model > 0).
 * It handles both hindcast and projection years by carrying over the last
 * hindcast year's empirical data.
 * * @param weight_hat [ref] 4D array to store calculated weights (wt_index, sex, age, year)
 * @param length_hat [ref] 4D array for estimated lengths
 * @param growth_matrix [ref] 5D array for age-length transition matrices
 * @param weight Empirical weight data array
 * @param growth_model Integer vector indicating growth type (0=empirical, >0=estimated)
 * @param nspp Number of species
 * @param nyrs Total number of years (hindcast + projection)
 * @param nyrs_hind Number of hindcast years
 * @param n_flt Number of fleets
 * @param flt_spp Vector mapping fleet to species index
 * @param flt_month Vector mapping fleet to month of operation
 * @param nsex Vector of number of sexes per species
 * @param nages Vector of number of ages per species
 * @param pop_wt_index Index for population biomass weights
 * @param ssb_wt_index Index for spawning stock biomass weights
 * @param flt_wt_index Index for fleet-specific weights
 * @param spawn_month Vector of spawning months per species
 * [Other parameters for calculate_growth: lengths, nlengths, minage, maxage,
 * growth_parameters, growth_ln_sd, weight_length_pars]
 */
template <class Type>
void calculate_weight(
    array<Type> &weight_hat,
    array<Type> &length_hat,
    array<Type> &growth_matrix,
    const array<Type> &weight,
    const ivector &growth_model,
    int nspp,
    int nyrs,
    int nyrs_hind,
    int n_flt,
    const ivector &flt_spp,
    const vector<Type> &flt_month,
    const ivector &nsex,
    const ivector &nages,
    const ivector &pop_wt_index,
    const ivector &ssb_wt_index,
    const ivector &flt_wt_index,
    const vector<Type> &spawn_month,
    const vector<Type> &lengths,
    const ivector &nlengths,
    const ivector &minage,
    const ivector &maxage,
    const array<Type> &growth_parameters,
    const array<Type> &growth_ln_sd,
    const array<Type> &weight_length_pars
) {
  int yr_ind;
  int wt_idx_pop;
  int wt_idx_ssb;
  int wt_idx_flt;

  // 1. POPULATION WEIGHT-AT-AGE
  for (int sp = 0; sp < nspp; sp++) {

    wt_idx_pop = (nspp - 1) * 2 * sp;
    wt_idx_ssb = (nspp - 1) * 2 * sp + 1;

    // -- 1.1. Empirical weight-at-age
    if (growth_model(sp) == 0) {
      for (int sex = 0; sex < nsex(sp); sex++) {
        for (int age = 0; age < nages(sp); age++) {
          for (int yr = 0; yr < nyrs; yr++) {

            // Handle projection logic
            yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);

            // Biomass weight
            weight_hat(wt_idx_pop, sex, age, yr) = weight(pop_wt_index(sp), sex, age, yr_ind);

            // SSB weight
            weight_hat(wt_idx_ssb, sex, age, yr) = weight(ssb_wt_index(sp), sex, age, yr_ind);
          }
        }
      }
    }

    // -- 1.2. Estimated growth
    if (growth_model(sp) > 0) {
      // Biomass weight (beginning of year / month 0)
      calculate_growth(wt_idx_pop, sp, Type(0.0), nspp, nyrs,
                       nsex, nages, lengths, nlengths, minage, maxage,
                       growth_parameters, growth_ln_sd, weight_length_pars,
                       growth_model, length_hat, growth_matrix, weight_hat);

      // SSB weight (at month of spawning)
      calculate_growth(wt_idx_ssb, sp, spawn_month(sp) / Type(12.0), nspp, nyrs,
                       nsex, nages, lengths, nlengths, minage, maxage,
                       growth_parameters, growth_ln_sd, weight_length_pars,
                       growth_model, length_hat, growth_matrix, weight_hat);
    }
  }

  // 2. FLEET WEIGHT-AT-AGE
  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    Type mo = flt_month(flt);
    wt_idx_flt = nspp * 2 + flt;

    // -- 2.1. Empirical weight-at-age
    if (growth_model(sp) == 0) {
      for (int sex = 0; sex < nsex(sp); sex++) {
        for (int age = 0; age < nages(sp); age++) {
          for (int yr = 0; yr < nyrs; yr++) {

            yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);
            weight_hat(wt_idx_flt, sex, age, yr) = weight(flt_wt_index(flt), sex, age, yr_ind);
          }
        }
      }
    }

    // -- 2.2. Estimated growth
    if (growth_model(sp) > 0) {
      calculate_growth(wt_idx_flt, sp, mo / Type(12.0), nspp, nyrs,
                       nsex, nages, lengths, nlengths, minage, maxage,
                       growth_parameters, growth_ln_sd, weight_length_pars,
                       growth_model, length_hat, growth_matrix, weight_hat);
    }
  }
}



