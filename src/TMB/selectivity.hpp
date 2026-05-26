#ifndef SELECTIVITY_HPP
#define SELECTIVITY_HPP

/**
 * @brief Performs final normalization and projection of fishery selectivity (age/length) across all fleets.
 * * @details This function iterates through all fleets to process selectivity data.
 * It handles:
 * 1. Zeroing out selectivity for ages below bin_first_selected.
 * 2. Normalizing selectivity values based on a single age (sel_norm_bin1 >= 0),
 * an age range average (sel_norm_bin2 >= 0), or the maximum value
 * across all ages/sexes (sel_norm_bin1 between -1 and -499).
 * 3. Projecting the final year of hindcast selectivity into future projection years.
 *
 * @param flt fleet to iterate through.
 * @param nyrs_hind Total number of years in the hindcast period.
 * @param nyrs Total number of years including projection.
 * @param flt_spp Function/Array mapping fleet index to species index.
 * @param flt_sel_type Function/Array mapping fleet index to selectivity type.
 * @param flt_sel_dim Age or length based selectivity.
 * @param bin_first_selected Array/function returning the minimum age/length bin of selection for a fleet.
 * @param nages Array/function returning the number of ages for a species.
 * @param nlengths Array/function returning the number of length bins for a species.
 * @param nsex Array/function returning the number of sexes for a species.
 * @param sel_norm_bin1 Array/function returning the normalization age or control flag.
 * @param sel_norm_bin2 Array/function returning the upper bound for age-range normalization.
 * @param sel_at_age 4D container (fleet, sex, age, year) modified in-place.
 */
template<class Type>
void normalize_and_project_selectivity(
    const int& flt,
    const int& nyrs_hind,
    const int& nyrs,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_sel_dim,
    const vector<int>&  bin_first_selected,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  nsex,
    const vector<int>&  sel_norm_bin1,
    const vector<int>&  sel_norm_bin2,
    array<Type> &selectivity // Modified by reference
) {
  Type max_sel; // Local declaration for safety

  int sp = flt_spp(flt);
  int sel_type = flt_sel_type(flt);
  int nbins = (flt_sel_dim(flt) == 0) ? nages(sp) : nlengths(sp);

  // Ages not selected
  for(int yr = 0; yr < nyrs_hind; yr++) {
    for(int bin = 0; bin < nbins; bin++){
      for(int sex = 0; sex < nsex(sp); sex++){
        if(bin < bin_first_selected(flt)){
          selectivity(flt, sex, bin, yr) = 0;
        }
      }
    }
  }

  // Normalize selectivity
  if(sel_norm_bin1(flt) > -500){

    // 1. Normalize by selectivity by specific bin or bin-range
    if((sel_norm_bin1(flt) >= 0) && (sel_type != 5) && (sel_type != 12)) { // Dont normalize hake type selex
      for(int yr = 0; yr < nyrs_hind; yr++) {
        for(int sex = 0; sex < nsex(sp); sex++){

          // Single-normalization bin
          if((sel_norm_bin1(flt) >= 0) && (sel_norm_bin2(flt) < 0)){
            max_sel = 0.001;
            max_sel = max2(max_sel, selectivity(flt, sex, sel_norm_bin1(flt), yr));
          }

          // Normalize by bin range between max lower and max upper
          if((sel_norm_bin1(flt) >= 0) && (sel_norm_bin2(flt) >= 0)){
            max_sel = 0;
            for(int bin = sel_norm_bin1(flt); bin <= sel_norm_bin2(flt); bin++) {
              max_sel += selectivity(flt, sex, bin, yr)/(sel_norm_bin2(flt) - sel_norm_bin1(flt) + 1);
            }
          }

          // Normalize
          for(int bin = 0; bin < nbins; bin++){
            selectivity(flt, sex, bin, yr) /= max_sel;
          }
        }
      }
    }

    // 2. Normalize by max for each fishery and year across bins, and sexes
    // - Don't for hake non-parametric
    if((sel_type != 5) && (sel_type != 12) && (sel_norm_bin1(flt) < 0) && (sel_norm_bin1(flt) > -500)) {
      for(int yr = 0; yr < nyrs_hind; yr++) {
        max_sel = 0;
        for(int bin = 0; bin < nbins; bin++){
          for(int sex = 0; sex < nsex(sp); sex++){
            // Find max
            if(selectivity(flt, sex, bin, yr) > max_sel){
              max_sel = selectivity(flt, sex, bin, yr);
            }
          }
        }

        // Normalize selectivity
        for(int bin = 0; bin < nbins; bin++){
          for(int sex = 0; sex < nsex(sp); sex++){
            selectivity(flt, sex, bin, yr) /= max_sel;
          }
        }
      }
    }
  }

  // Project outwards assuming final years selectivity
  for(int yr = nyrs_hind; yr < nyrs; yr++) {
    for(int bin = 0; bin < nbins; bin++){
      for(int sex = 0; sex < nsex(sp); sex++){
        selectivity(flt, sex, bin, yr) = selectivity(flt, sex, bin, nyrs_hind - 1);
      }
    }
  }
}


/**
 * @brief Converts length-based selectivity to age-based selectivity and stores it in a 4D array.
 * * @details This function maps selectivity from length bins to age classes by
 * integrating over a growth transition matrix (Prob(Length | Age)). The resulting
 * age-based selectivity is stored directly into the provided 'sel_at_age' array for the
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
 * @param sel_at_age Reference to 4D output array [flt, sex, age, yr] to be updated in-place.
 */
template <class Type>
void convert_length_selectivity(
    const int& flt,
    const int& sp,
    const int& sex,
    const int& yr,
    const int& wtind,
    const vector<int>& nages,
    const vector<int>& nlengths_pop,
    array<Type> sel_length_pop,
    array<Type> growth_matrix_pop,
    array<Type>& sel_at_age // Modified by reference
) {
  // SS3-parity convolution: integrate sel(L) over P(L|age) on the FINE
  // pop grid. Evaluating on data bins under-resolves the sel curve where
  // it changes rapidly (peak / asc limb) and biases sel-at-age.
  for(int age = 0; age < nages(sp); age++) {
    sel_at_age(flt, sex, age, yr) = 0.0;
    for(int ln = 0; ln < nlengths_pop(sp); ln++) {
      sel_at_age(flt, sex, age, yr) +=
        growth_matrix_pop(wtind, sex, age, ln, yr) * sel_length_pop(flt, sex, ln, yr);
    }
  }
}


/**
 * @brief Unified function to calculate selectivity across fleets, sexes, bins (age or length), and years.
 * * @details This function calculates selectivity for the hindcast period,
 * normalizes it, projects it into the future, and—if the base calculation was length-based—uses a
 * growth transition matrix to output the equivalent age-based selectivity.
 * * @section functional_forms Selectivity Models:
 * The dimension (age vs length) is determined by `flt_sel_dim` (0 = age, 1 = length) and is shared
 * across all parametric cases; the age/length distinction is handled inside each case, not by a
 * separate case number. After calculation, length-based selectivity is converted to age-based via
 * `convert_length_selectivity()`, and all selectivities are normalized/projected via
 * `normalize_and_project_selectivity()`.
 *
 * - Case 0: Empirical — fixed selectivity read from `emp_sel` input; assumed age-based.
 *
 * - Case 1: Logistic — 2-parameter ascending logistic [Age or Length]
 *     sel(x) = 1 / (1 + exp(-slope * (x - inflection)))
 *     Parameters: log_sel_slp(0) = log(slope); sel_inf(0) = inflection point
 *
 * - Case 2: Non-parametric (Ianelli/AMAK style) [Age or Length]
 *     Bin-specific log-scale coefficients, normalized by mean; penalized for monotonicity and curvature.
 *     Parameters: sel_coff[flt, sex, bin] for bins 0..N_sel_bins-1
 *
 * - Case 3: Double Logistic (Dorn and Methot 1990) [Age or Length]
 *     sel(x) = asc(x) * (1 - desc(x)), where each limb is a 2-parameter logistic.
 *     Parameters: log_sel_slp(0/1) = log(ascending/descending slope); sel_inf(0/1) = ascending/descending inflection
 *
 * - Case 4: Descending Logistic [Age or Length]
 *     sel(x) = 1 - 1/(1 + exp(-slope * (x - inflection)))
 *     Parameters: log_sel_slp(1) = log(slope); sel_inf(1) = inflection point
 *
 * - Case 5: Non-parametric (Hake/Taylor style) [Age or Length]
 *     Cumulative sum of bin coefficients, exponentiated and normalized to max = 1.
 *     Parameters: sel_coff[flt, sex, bin] for bins bin_first_selected..N_sel_bins-1
 *
 * - Case 6: 2D AR1 (age x year) [Age or Length]
 *     Logistic-transformed bin-year deviates with AR1 correlation across age and year.
 *
 * - Case 7: 3D AR1 (age x cohort x year, Cheng et al 2024) [Age or Length]
 *     Conditional-variance parameterisation of a 3D random field across age, cohort and year.
 *
 * - Case 8: Double Normal [Age or Length]
 *     Gaussian ascending limb left of peak, Gaussian descending limb right of peak, blended
 *     smoothly via a steep logistic weight so the function is twice differentiable everywhere.
 *       sel(x) = (1-w(x)) * exp(-0.5*((x-peak)/sigma_asc)^2)
 *              +    w(x)  * exp(-0.5*((x-peak)/sigma_desc)^2)
 *       w(x)  = 1 / (1 + exp(-20*(x - peak)))   [~0 left of peak, ~1 right]
 *     Parameters (reuses existing arrays):
 *       sel_inf(0)     = peak bin (mode); TV deviate: sel_inf_dev(0)
 *       log_sel_slp(0) = log(sigma_ascending);  TV deviate: log_sel_slp_dev(0)
 *       log_sel_slp(1) = log(sigma_descending); TV deviate: log_sel_slp_dev(1)
 *       sel_inf(1) is intentionally unused (fixed at initial value via map).
 * * @param nspp Number of species.
 * @param n_flt Number of fleets (fisheries and surveys).
 * @param nyrs Total years (including hindcast and projection).
 * @param nyrs_hind Number of hindcast years.
 * @param styr Start year (for indexing empirical data).
 * @param nsex Vector of number of sexes per species.
 * @param nages Vector of max age bins per species.
 * @param nlengths Vector of max length bins per species.
 * @param lengths Matrix of length bin boundaries (nspp x max_length_bins).
 * @param flt_spp Vector mapping fleet index to species index.
 * @param flt_sel_type Vector of selectivity form.
 * @param flt_sel_dim Age or length based selectivity
 * @param bin_first_selected Vector of minimum age/length bins selected per fleet.
 * @param flt_n_sel_bins Vector of maximum bins with estimated coefficients per fleet.
 * @param sel_norm_bin1 Vector of bins used for normalization (or control flags) per fleet.
 * @param sel_norm_bin2 Vector of upper bins for mean normalization per fleet.
 * @param emp_sel_obs Matrix of empirical selectivity observation values.
 * @param emp_sel_ctl Control matrix for empirical selectivity indexing.
 * @param log_sel_slp Matrix of log-scale logistic slope parameters.
 * @param log_sel_slp_dev 4D array of yearly slope deviations.
 * @param sel_inf Matrix of bin-at-inflection parameters.
 * @param sel_inf_dev 4D array of yearly inflection deviations.
 * @param sel_coff 3D array of bin-specific coefficients for non-parametric models.
 * @param sel_coff_dev 4D array of yearly coefficient deviations.
 * @param avg_sel [Output] 3D array of average selectivity per fleet/sex/year (modified by reference).
 * @param non_par_sel [Output] 4D array of unnormalized non-parametric selectivity (modified by reference).
 * @param sel_at_length [Output] 4D array of final length-based selectivity (modified by reference).
 * @param sel_at_age [Output] 4D array of final age-based selectivity (modified by reference).
 * @param growth_matrix 5D array defining the length-at-age probability transition matrix.
 */
template<class Type>
void calculate_selectivity(
    const int& nspp,
    const int& n_flt,
    const int& nyrs,
    const int& nyrs_hind,
    const int& styr,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  nlengths_pop,
    matrix<Type> lengths,
    matrix<Type> lengths_pop,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_sel_dim,
    const vector<int>&  bin_first_selected,
    const vector<int>&  flt_n_sel_bins,
    const vector<int>&  sel_norm_bin1,
    const vector<int>&  sel_norm_bin2,
    matrix<Type> emp_sel_obs,
    matrix<int>& emp_sel_ctl,
    array<Type> log_sel_slp,
    array<Type> log_sel_slp_dev,
    array<Type> sel_inf,
    array<Type> sel_inf_dev,
    array<Type> sel_coff,
    array<Type> sel_coff_dev,
    array<Type> &avg_sel,
    array<Type> &non_par_sel,
    array<Type> &sel_at_length,
    array<Type> &sel_at_length_pop,
    array<Type> &sel_at_age,
    array<Type> growth_matrix,
    array<Type> growth_matrix_pop
) {
  sel_at_age.setZero();
  sel_at_length.setZero();
  sel_at_length_pop.setZero();
  Type max_sel = 0;
  Type avgsel_tmp = 0;

  // --- 1. EMPIRICAL SELECTIVITY (Assumed Age-based from input) ---
  for (int i = 0; i < emp_sel_obs.rows(); i++) {
    int flt = emp_sel_ctl(i, 0) - 1;
    int sp  = emp_sel_ctl(i, 1) - 1;
    int s_idx = emp_sel_ctl(i, 2);
    int yr = emp_sel_ctl(i, 3) - styr;

    if (flt_sel_type(flt) == 0 && yr < nyrs_hind) {
      int s_start = (s_idx == 0) ? 0 : s_idx - 1;
      int s_end   = (s_idx == 0 && nsex(sp) == 2) ? 1 : s_start;

      for (int sex = s_start; sex <= s_end; sex++) {
        for (int age = 0; age < nages(sp); age++) {
          if (!isNA(emp_sel_obs(i, age))) {
            sel_at_age(flt, sex, age, yr) = emp_sel_obs(i, age);
          }
        }
      }
    }
  }

  // --- 2. ESTIMATED SELECTIVITY (Unified Age & Length) ---
  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    int sel_type = flt_sel_type(flt);
    if (sel_type == 0) continue;

    bool is_length_based = flt_sel_dim(flt) == 1;
    int nbins =  is_length_based? nlengths(sp) : nages(sp);
    int n_sel_bins = flt_n_sel_bins(flt);
    Type binwidth = is_length_based ? (lengths(sp, 1) - lengths(sp, 0)) : Type(1.0);

    for (int yr = 0; yr < nyrs_hind; yr++) {
      for (int sex = 0; sex < nsex(sp); sex++) {

        switch (sel_type) {
        case 1: // Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1);
            Type slope = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr));
            Type inf   = sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr);
            Type val = 1.0 / (1.0 + exp(-slope * (x_val - inf)));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 2: // Non-parametric (Ianelli style)
          for(int bin = 0; bin < n_sel_bins; bin++) {
            non_par_sel(flt, sex, bin, yr) = sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr);
            avg_sel(flt, sex, yr) += exp(sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr));
          }
          avg_sel(flt, sex, yr) = log(avg_sel(flt, sex, yr) / n_sel_bins);

          for(int bin = n_sel_bins; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, n_sel_bins - 1, yr);
          }

          avgsel_tmp = 0;
          for(int bin = 0; bin < nbins; bin++) {
            avgsel_tmp += exp(non_par_sel(flt, sex, bin, yr));
          }
          avgsel_tmp = log(avgsel_tmp / nbins);

          for(int bin = 0; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) -= avgsel_tmp;
            non_par_sel(flt, sex, bin, yr) = exp(non_par_sel(flt, sex, bin, yr));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            else sel_at_age(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
          }
          break;

        case 3: // Double Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1);
            Type slp1 = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr));
            Type inf1 = sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr);
            Type slp2 = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr));
            Type inf2 = sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr);
            Type val = (1.0 / (1.0 + exp(-slp1 * (x_val - inf1)))) * (1.0 - (1.0 / (1.0 + exp(-slp2 * (x_val - inf2)))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 4: // Descending Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? (lengths(sp, bin) + 0.5 * binwidth) : Type(bin + 1);
            Type slp2 = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr));
            Type inf2 = sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr);
            Type val = (1.0 - (1.0 / (1.0 + exp(-slp2 * (x_val - inf2)))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 5: // Hake/Taylor Non-parametric
          {
            Type cum_sum = 0;
            for (int bin = bin_first_selected(flt); bin < nbins; bin++) {
              if (bin < n_sel_bins) {
                cum_sum += sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr);
              }
              if (is_length_based) sel_at_length(flt, sex, bin, yr) = cum_sum;
              else sel_at_age(flt, sex, bin, yr) = cum_sum;
            }

            // Normalize inside year/sex block
            max_sel = -1e10;
            if (sel_norm_bin1(flt) >= 0 && sel_norm_bin2(flt) < 0) {
              max_sel = is_length_based ? sel_at_length(flt, sex, sel_norm_bin1(flt), yr) : sel_at_age(flt, sex, sel_norm_bin1(flt), yr);
            } else {
              for(int bin=0; bin<nbins; bin++) {
                Type val = is_length_based ? sel_at_length(flt, sex, bin, yr) : sel_at_age(flt, sex, bin, yr);
                if(val > max_sel) max_sel = val;
              }
            }
            for (int bin = 0; bin < nbins; bin++) {
              if (is_length_based) sel_at_length(flt, sex, bin, yr) = exp(sel_at_length(flt, sex, bin, yr) - max_sel);
              else sel_at_age(flt, sex, bin, yr) = exp(sel_at_age(flt, sex, bin, yr) - max_sel);
            }
          }
          break;

        case 8: { // SS3 pattern 24 -- length-based Double Normal with plateau
          // 6 parameters per fleet/sex (SS3 Pn ↔ Rceattle slot):
          //   sel_inf(0)     = P1 peak (mode start of plateau)
          //   sel_inf(1)     = P6 final/right-tail floor (on logit scale)
          //   sel_inf(2)     = P5 init/left-tail floor   (on logit scale;
          //                     -inf -> 0 = full Gaussian ascending limb)
          //   log_sel_slp(0) = P3 log(sigma_asc)
          //   log_sel_slp(1) = P4 log(sigma_desc)
          //   log_sel_slp(2) = P2 top-width logit, controls plateau width:
          //                     peak2 = peak + binwidth + (xmax - peak - binwidth)
          //                             / (1 + exp(-P2))
          // SS3 formula (selectivity.tpl pattern 24):
          //   asc(L) = exp(-(L-peak)^2 / exp(P3))
          //   asc_scaled(L) = init + (1 - init) * (asc(L) - asc(xmin))/(1 - asc(xmin))
          //   dsc(L) = exp(-(L-peak2)^2 / exp(P4))
          //   dsc_scaled(L) = 1 + (final - 1) * (dsc(L) - 1)/(dsc(xmax) - 1)
          //   join1(L) = 1 / (1 + exp(-20/(1+|L-peak|)  * (peak - L)))
          //   join2(L) = 1 / (1 + exp(-20/(1+|L-peak2|) * (peak2 - L)))
          //   sel(L) = asc_scaled*(1-join1) + join1*((1-join2) + dsc_scaled*join2)
          // Note the SS3 Gaussian denom is exp(P3) NOT exp(2*P3): SS3 uses
          // (L-peak)^2 / exp(P3) where exp(P3) plays the role of 2*sigma^2.
          Type peak     = sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr);
          Type final_lt = sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr);
          Type init_lt  = sel_inf(2, flt, sex) + sel_inf_dev(2, flt, sex, yr);
          Type upselex   = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr));
          Type downselex = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr));
          Type top_lt    = log_sel_slp(2, flt, sex) + log_sel_slp_dev(2, flt, sex, yr);
          Type init      = Type(1.0) / (Type(1.0) + exp(-init_lt));
          Type finalv    = Type(1.0) / (Type(1.0) + exp(-final_lt));

          // xmin / xmax = first / last bin midpoint on the integration grid
          // (length midpoints when length-based; age 1..nages when age-based).
          Type xmin = is_length_based ? (lengths(sp, 0)         + Type(0.5) * binwidth) : Type(1);
          Type xmax = is_length_based ? (lengths(sp, nbins - 1) + Type(0.5) * binwidth) : Type(nbins);
          Type peak2 = peak + binwidth + (xmax - peak - binwidth) / (Type(1.0) + exp(-top_lt));

          // Normalization anchors so asc -> init at xmin, dsc -> final at xmax
          Type t1min = exp(-(xmin - peak)  * (xmin - peak)  / upselex);
          Type t2min = exp(-(xmax - peak2) * (xmax - peak2) / downselex);

          // Helper lambda to evaluate the curve at a length L. Pulls in
          // peak, peak2, init, finalv, upselex, downselex, t1min, t2min
          // by capture. Used for both the DATA-bin and POP-bin loops below.
          auto eval_p24 = [&](Type L) {
            Type asc  = exp(-(L - peak)  * (L - peak)  / upselex);
            Type dsc  = exp(-(L - peak2) * (L - peak2) / downselex);
            Type asc_scaled = init + (Type(1.0) - init) * (asc - t1min) / (Type(1.0) - t1min);
            Type dsc_scaled = Type(1.0) + (finalv - Type(1.0)) * (dsc - Type(1.0)) / (t2min - Type(1.0));
            // SS3 location-dependent join sigmoids. Sign convention from
            // SS3 selectivity.tpl: join1=0 LEFT of peak, =1 RIGHT;
            // join2=0 LEFT of peak2, =1 RIGHT.
            Type denom1 = Type(1.0) + ((L - peak)  >= Type(0) ? (L - peak)  : (peak  - L));
            Type denom2 = Type(1.0) + ((L - peak2) >= Type(0) ? (L - peak2) : (peak2 - L));
            Type join1  = Type(1.0) / (Type(1.0) + exp(-(Type(20.0) / denom1) * (L - peak )));
            Type join2  = Type(1.0) / (Type(1.0) + exp(-(Type(20.0) / denom2) * (L - peak2)));
            return asc_scaled * (Type(1.0) - join1)
                   + join1 * ((Type(1.0) - join2) + dsc_scaled * join2);
          };

          for (int bin = 0; bin < nbins; bin++) {
            Type L = is_length_based ? (lengths(sp, bin) + Type(0.5) * binwidth) : Type(bin + 1);
            Type val = eval_p24(L);
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else                 sel_at_age(flt, sex, bin, yr) = val;
          }

          // For length-based pattern 24, also evaluate at the FINE POP-grid
          // midpoints and store in sel_at_length_pop. The sel-to-age
          // convolution then happens on the pop grid (SS3 parity), which is
          // critical when data bins are coarser than pop bins.
          if (is_length_based) {
            Type pop_bw = lengths_pop(sp, 1) - lengths_pop(sp, 0);
            for (int lp = 0; lp < nlengths_pop(sp); lp++) {
              Type Lp = lengths_pop(sp, lp) + Type(0.5) * pop_bw;
              sel_at_length_pop(flt, sex, lp, yr) = eval_p24(Lp);
            }
          }
          break;
        }

        case 6: // 2D AR1-age x year
        case 7: // 3D AR1 conditional variance
          for (int bin = 0; bin < nbins; bin++) {

            // Safely cap the bin index to prevent out-of-bounds access
            int active_bin = (bin < n_sel_bins) ? bin : (n_sel_bins - 1);
            Type val = 1.0 / (1.0 + exp(-(sel_coff(flt, sex, active_bin) + sel_coff_dev(flt, sex, active_bin, yr))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;
        }

      } // End sex
    } // End year

    // --- 3. NORMALIZATION & PROJECTION ---
    normalize_and_project_selectivity(
      flt, nyrs_hind, nyrs, flt_spp, flt_sel_type, flt_sel_dim, bin_first_selected, nages, nlengths, nsex,
      sel_norm_bin1, sel_norm_bin2,
      is_length_based ? sel_at_length : sel_at_age
    );

    // --- 4. CONVERT LENGTH TO AGE (If necessary) ---
    if (is_length_based) {
      for (int yr = 0; yr < nyrs; yr++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          int wtind = nspp * 2 + flt;
          convert_length_selectivity(
            flt, sp, sex, yr, wtind, nages, nlengths,
            sel_at_length, growth_matrix, sel_at_age
          );
        }
      }
    }
  } // End fleet
}

#endif
