#ifndef SELECTIVITY_HPP
#define SELECTIVITY_HPP

/**
 * @brief Zeroes, normalizes and projects ONE fleet's selectivity (age or length based).
 *
 * @details Bins below bin_first_selected are zeroed; the last hindcast year carries
 * across the projection. Normalization runs only when Sel_norm_bin is set (default Off),
 * never for forms 5, 11 and 15: scaling to 1 at sel_norm_bin1, to a mean of 1 over
 * [bin1, bin2], or to the maximum when bin1 is negative. Under Sel_norm_scope =
 * AcrossSexes the less-selected sex deliberately peaks below 1.
 *
 * @param flt Fleet index.
 * @param nyrs_hind Years in the hindcast.
 * @param nyrs Years, hindcast plus projection.
 * @param flt_spp Fleet index to species index.
 * @param flt_sel_type Selectivity form per fleet.
 * @param flt_sel_dim Age (0) or length (1) based selectivity.
 * @param bin_first_selected First selected age/length bin per fleet, 0-based.
 * @param nages Age bins per species.
 * @param nlengths_pop Population length bins per species.
 * @param nsex Sexes modelled per species: 1 = combined, 2 = both.
 * @param sel_norm_bin1 Normalization bin, or a control flag when negative.
 * @param sel_norm_bin2 Upper bin for mean-over-range normalization.
 * @param sel_norm_scope Per fleet: 0 = WithinSex (each sex divided by its own reference,
 *   so both reach 1), 1 = AcrossSexes (one pooled reference, so the less-selected sex
 *   stays below 1). Orthogonal to sel_norm_bin1/2, which say WHERE it is taken.
 * @param selectivity 4D (fleet, sex, age or length bin, year), modified in place.
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
    const vector<int>&  nlengths_pop,
    const vector<int>&  nsex,
    const vector<int>&  sel_norm_bin1,
    const vector<int>&  sel_norm_bin2,
    const vector<int>&  sel_norm_scope,
    array<Type> &selectivity // Modified by reference
) {
  int sp = flt_spp(flt);
  int sel_type = flt_sel_type(flt);
  // A length-based curve lives on the POPULATION bins, so everything below runs on them.
  int nbins = (flt_sel_dim(flt) == 0) ? nages(sp) : nlengths_pop(sp);

  // Bins below the first selected bin
  for(int yr = 0; yr < nyrs_hind; yr++) {
    for(int bin = 0; bin < nbins; bin++){
      for(int sex = 0; sex < nsex(sp); sex++){
        if(bin < bin_first_selected(flt)){
          selectivity(flt, sex, bin, yr) = 0;
        }
      }
    }
  }

  // Excluded forms: Hake (5) normalizes per sex inside the switch, LogisticPM (11) reuses
  // sel_norm_bin1/2 as a penalty bin-range, and DoubleNormalSS3 (15) is SS3's
  // unnormalized pattern 24, whose plateau is 1 by construction. Codes 10, 12
  // and 14 are unused (R/0-switches.R).
  if((sel_norm_bin1(flt) > -500) && (sel_type != 5) && (sel_type != 12) && (sel_type != 11) && (sel_type != 15)) {
    bool at_bin     = (sel_norm_bin1(flt) >= 0);
    bool over_range = at_bin && (sel_norm_bin2(flt) >= 0);
    bool across     = (sel_norm_scope(flt) == 1);

    vector<Type> ref(nsex(sp));

    for(int yr = 0; yr < nyrs_hind; yr++) {

      // Reference value for each sex.
      for(int sex = 0; sex < nsex(sp); sex++){
        if(over_range){
          // Mean selectivity over the bin range
          ref(sex) = 0;
          for(int bin = sel_norm_bin1(flt); bin <= sel_norm_bin2(flt); bin++) {
            ref(sex) += selectivity(flt, sex, bin, yr)/(sel_norm_bin2(flt) - sel_norm_bin1(flt) + 1);
          }
        } else if(at_bin){
          // Single normalization bin, floored so the division stays finite
          ref(sex) = max2(Type(0.001), selectivity(flt, sex, sel_norm_bin1(flt), yr));
        } else {
          // Maximum over bins. Fold rather than branch: TMB tapes once, so an `if`
          // on an AD Type would freeze the argmax bin at its initial position.
          ref(sex) = 0;
          for(int bin = 0; bin < nbins; bin++){
            ref(sex) = max2(ref(sex), selectivity(flt, sex, bin, yr));
          }
        }
      }

      // One divisor for every sex, so their relative levels survive.
      if(across){
        Type pooled = ref(0);
        for(int sex = 1; sex < nsex(sp); sex++) pooled = max2(pooled, ref(sex));
        for(int sex = 0; sex < nsex(sp); sex++) ref(sex) = pooled;
      }

      for(int bin = 0; bin < nbins; bin++){
        for(int sex = 0; sex < nsex(sp); sex++){
          selectivity(flt, sex, bin, yr) /= ref(sex);
        }
      }
    }
  }

  // Projection years repeat the last hindcast year's selectivity
  for(int yr = nyrs_hind; yr < nyrs; yr++) {
    for(int bin = 0; bin < nbins; bin++){
      for(int sex = 0; sex < nsex(sp); sex++){
        selectivity(flt, sex, bin, yr) = selectivity(flt, sex, bin, nyrs_hind - 1);
      }
    }
  }
}


/**
 * @brief Converts one fleet/sex/year of length-based selectivity to age-based.
 *
 * @details Integrates selectivity-at-length over the growth transition matrix
 * Prob(Length | Age) on the POPULATION length grid both arrays live on, so the product is
 * formed at the resolution the curve is defined on, not on a bin average.
 *
 * @tparam Type The numeric type (e.g., double, TMB::adouble).
 * @param flt Fleet index.
 * @param sp Species index.
 * @param sex Sex index.
 * @param yr Year index.
 * @param wtind Index for the weight-at-age/growth transition matrix type.
 * @param nages Age bins per species.
 * @param nlengths_pop Population length bins per species.
 * @param sel_length 4D array of selectivity-at-length [flt, sex, pop bin, yr].
 * @param growth_matrix 5D transition matrix [wtind, sex, age, pop bin, yr].
 * @param sel_at_age 4D output array [flt, sex, age, yr] updated in-place.
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
    array<Type>& sel_length,
    array<Type>& growth_matrix,
    array<Type>& sel_at_age // Modified by reference
) {
  for(int age = 0; age < nages(sp); age++) {

    sel_at_age(flt, sex, age, yr) = 0.0;

    for(int lp = 0; lp < nlengths_pop(sp); lp++) {
      sel_at_age(flt, sex, age, yr) += growth_matrix(wtind, sex, age, lp, yr) * sel_length(flt, sex, lp, yr);
    }
  }
}


/**
 * @brief Calculates selectivity for every fleet, sex, bin (age or length) and year.
 *
 * @details Builds each fleet's hindcast curve from its Selectivity form, applies the
 * per-sex apical height, then normalizes and projects via
 * `normalize_and_project_selectivity()`; `convert_length_selectivity()` carries a
 * length-based curve to ages. `flt_sel_dim` (0 = age, 1 = length) picks the grid per
 * case, so a length-based curve has no case of its own.
 *
 * @section functional_forms Selectivity forms
 * The codes are fleet_control$Selectivity; R/0-switches.R is the authority. The
 * parametric forms (1, 3, 4, 8, 11, 15) take linkage offsets in the slots named at their
 * case; 15 has no deviate array and varies only through linkages. The coefficient forms
 * (2, 5, 6, 7, 9, 13) take annual deviates but no linkage offset -- build_selectivity()
 * refuses one, since mean-centring cancels a per-year offset. x below is the length-bin
 * midpoint (cm) or, on an age-based fleet, the 1-based age-bin ordinal bin + 1
 * (LogisticPM uses bin + 1.5); minage is NOT added, so an inflection or peak is on the
 * bin scale, not the age scale.
 *
 * - 0 Fixed: read from emp_sel_obs; age-based.
 * - 1 Logistic: ascending logistic on log_sel_slp(0) and sel_inf(0).
 * - 2 NonParametric (Ianelli et al. 2018, AMAK): one log-scale coefficient per bin,
 *     sel_coff bins 0..N_sel_bins-1 with the rest plateauing, mean-centred, and penalized
 *     for monotonicity and curvature.
 * - 3 DoubleLogistic (Dorn and Methot 1990): ascending logistic x (1 - descending).
 * - 4 DescendingLogistic: the descending limb alone, on log_sel_slp(1)/sel_inf(1).
 * - 5 Hake (Taylor et al. 2014): cumulative sum of sel_coff from bin_first_selected up,
 *     exponentiated and normalized per sex to 1 at one sel_norm_bin1, else at its
 *     maximum.
 * - 6 2DAR1: inverse-logit of sel_coff + sel_coff_dev, the deviates carrying an AR1
 *     correlation across bin and year.
 * - 7 3DAR1 (Cheng et al 2024): conditional-variance parameterisation of a 3D random
 *     field across age, cohort and year.
 * - 8 DoubleNormal: four-parameter dome, Gaussian limbs blended by a steep logistic so
 *     the curve is twice differentiable, with a logit-scale right-tail floor analogous to
 *     SS3 pattern 24's P6.
 * - 9 NonParametricPM (AMAK "pm"): random walk on the per-year-centred log-selectivity,
 *     held flat above flt_sel_cap_bin; its increments are penalized, not integrated, so
 *     it is not integrable under random_sel.
 * - 11 LogisticPM (ADMB AMAK "pm" bottom-trawl survey): logistic with MULTIPLICATIVE
 *     slope and inflection deviates, plus a free first-bin log-selectivity.
 * - 13 NonParametricIntegrable: the NonParametric curve whose deviates carry a proper
 *     density, so random_sel can integrate them; Time_varying_sel picks Off, IID or
 *     RandomWalk.
 * - 15 DoubleNormalSS3: Stock Synthesis size pattern 24 (age pattern 20) on SS3's own six
 *     parameters -- P1 peak (cm or age), P2 top width (logit), P3/P4 ascending and
 *     descending widths (log), P5/P6 initial and final selectivity (logit) -- held in
 *     their own array sel_dn6. Unnormalized: the curve peaks at 1 by construction.
 *
 * @param nspp Number of species.
 * @param n_flt Number of fleets.
 * @param nyrs Years, hindcast plus projection.
 * @param nyrs_hind Years in the hindcast.
 * @param styr Start year (for indexing empirical data).
 * @param nsex Sexes modelled per species: 1 = combined, 2 = both.
 * @param nages Age bins per species.
 * @param nlengths Observation length bins per species.
 * @param lengths Data length bins per species; only the first is read, to locate SS3's
 *   startbin for DoubleNormalSS3.
 * @param flt_spp,flt_sel_type,flt_sel_dim,bin_first_selected,sel_norm_bin1,sel_norm_bin2
 *   Per fleet, as in normalize_and_project_selectivity() above.
 * @param flt_n_sel_bins Bins with estimated coefficients per fleet (N_sel_bins).
 * @param emp_sel_obs Empirical selectivity-at-age observations, one row per record.
 * @param emp_sel_ctl Fleet, species, sex and year index for each emp_sel_obs row.
 * @param log_sel_slp Log-scale logistic slope, two slots per fleet and sex.
 * @param log_sel_slp_dev Annual slope deviations.
 * @param sel_inf Inflection bin, two slots per fleet and sex.
 * @param sel_inf_dev Annual inflection deviations.
 * @param sel_coff Per-bin coefficients for the non-parametric and AR1 forms.
 * @param sel_coff_dev Annual coefficient deviations.
 * @param avg_sel [Output] Log mean selectivity over the estimated bins, per
 *   fleet/sex/year. Written by forms 2 and 13 only; forms 5 and 9 leave it at 0.
 * @param non_par_sel [Output] Non-parametric selectivity on the natural scale.
 * @param log_non_par_sel [Output] The same curve on the log scale, the scale the
 *   non-parametric selectivity penalties are written on.
 * @param sel_at_length,sel_at_age [Output] Final length- and age-based selectivity.
 * @param growth_matrix Length-at-age probability transition matrix.
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
    matrix<Type>& lengths,
    // A length-based curve is built on the POPULATION bins, as SS3 does; with no
    // population grid these equal nlengths / lengths. The curve stays on that grid,
    // and ceattle.cpp sums each data bin's run of population bins where it needs them.
    const vector<int>&  nlengths_pop,
    matrix<Type>& lengths_pop,
    const vector<int>&  flt_spp,
    const vector<int>&  flt_sel_type,
    const vector<int>&  flt_varying_sel,
    const vector<int>&  flt_sel_dim,
    const vector<int>&  bin_first_selected,
    const vector<int>&  flt_n_sel_bins,
    const vector<int>&  flt_sel_cap_bin,
    const vector<int>&  sel_norm_bin1,
    const vector<int>&  sel_norm_bin2,
    const vector<int>&  sel_norm_scope,
    const vector<int>&  flt_sel_start_yr,
    matrix<Type>& emp_sel_obs,
    matrix<int>& emp_sel_ctl,
    array<Type>& log_sel_slp,
    array<Type>& log_sel_slp_dev,
    array<Type>& sel_inf,
    array<Type>& sel_inf_dev,
    array<Type>& sel_coff,
    array<Type>& sel_coff_dev,
    array<Type> &avg_sel,
    array<Type> &non_par_sel,
    array<Type> &log_non_par_sel,
    array<Type> &sel_at_length,
    array<Type> &sel_at_age,
    array<Type>& growth_matrix,
    // Environmental-linkage offsets, one log-scale and one natural-scale tensor per
    // parameter family; all zero unless a selectivity linkage was supplied. Read, never
    // written, here -- non-const only because TMB's array operator() is not const.
    array<Type>& sel_slp_off,      array<Type>& sel_slp_off_nat,
    array<Type>& sel_inf_off,      array<Type>& sel_inf_off_nat,
    array<Type>& sel_coff_off,     array<Type>& sel_coff_off_nat,
    // Per-sex apical height, applied in section 2b below only when a selectivity linkage
    // names `apical` (sel_apical_on), so every other model's AD tape is unchanged.
    const int& sel_apical_on,
    array<Type>& log_sel_apical,   // [n_flt, max_sex]
    array<Type>& sel_apical_off,   array<Type>& sel_apical_off_nat,  // [n_flt, max_sex, nyrs]
    // DoubleNormalSS3 (case 15): SS3 pattern-24 parameters on SS3's scales, their
    // linkage offsets, and per fleet whether each end's scaling is used
    // (sel_dn6_ends(flt, 0) = initial, (flt, 1) = final; 1 = used, 0 = SS3 -999).
    array<Type>& sel_dn6,          // [6, n_flt, max_sex]
    array<Type>& sel_dn6_off,      array<Type>& sel_dn6_off_nat,     // [6, n_flt, max_sex, nyrs]
    matrix<int>& sel_dn6_ends      // [n_flt, 2]
) {
  sel_at_age.setZero();
  sel_at_length.setZero();
  Type max_sel = 0;
  Type avgsel_tmp = 0;

  // --- 1. FIXED SELECTIVITY (code 0, read from emp_sel_obs; assumed age-based) ---
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

    // NonParametricIntegrable (13) builds its curve two ways -- a random walk on the
    // base coefficients under Time_varying_sel = "RandomWalk" (4), otherwise the shared
    // non-parametric branch -- so the switch dispatches on both columns. The negative
    // case label is internal; it is not a Selectivity value.
    const int SEL_CASE_NP_INTEGRABLE_WALK = -13;
    int sel_case = (sel_type == 13 && flt_varying_sel(flt) == 4)
                     ? SEL_CASE_NP_INTEGRABLE_WALK : sel_type;

    bool is_length_based = flt_sel_dim(flt) == 1;
    // A length-based curve runs over the POPULATION bins, an age-based one over ages.
    int nbins =  is_length_based? nlengths_pop(sp) : nages(sp);
    int n_sel_bins = flt_n_sel_bins(flt);
    // Bin midpoints, SS3's len_bins_m. A population length grid need not be uniform --
    // SS3 models routinely widen the tail bins -- so each bin takes its OWN width,
    // keeping the curve on the same midpoints as the age-length key it multiplies.
    // A single length bin has no neighbour to take a width from, so it keeps its edge.
    vector<Type> xmid(nbins);
    for (int b = 0; b < nbins; b++) {
      if (!is_length_based) {
        xmid(b) = Type(b + 1);
      } else if (nbins < 2) {
        xmid(b) = lengths_pop(sp, b);
      } else {
        Type w = (b + 1 < nbins) ? (lengths_pop(sp, b + 1) - lengths_pop(sp, b))
                                 : (lengths_pop(sp, b) - lengths_pop(sp, b - 1));
        xmid(b) = lengths_pop(sp, b) + 0.5 * w;
      }
    }
    // One scalar width, as SS3's binwidth2: the population width, or the middle bin's
    // where the grid is not uniform (SS_readdata_330.tpl:1644). Read only by the
    // DoubleNormalSS3 peak2, as SS_selex.tpl:153 does.
    Type binwidth = (is_length_based && nbins > 1)
      ? (lengths_pop(sp, nbins / 2) - lengths_pop(sp, nbins / 2 - 1)) : Type(1.0);

    // Running log-scale walk state, carried across years by the NonParametricPM (9)
    // and NonParametricIntegrable (13) random walks.
    array<Type> np_unc(nsex(sp), nbins, nyrs_hind); np_unc.setZero();

    for (int yr = 0; yr < nyrs_hind; yr++) {
      for (int sex = 0; sex < nsex(sp); sex++) {

        switch (sel_case) {
        case 1: // Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? xmid(bin) : Type(bin + 1);
            // Linkage convention for all forms below: a log-scale parameter takes its
            // log-link offset inside the exp, its natural offset after; a natural-scale
            // one takes the identity offset additively, the log one multiplicatively.
            Type slope = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                             + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
            Type inf   = (sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr)
                          + sel_inf_off_nat(0, flt, sex, yr)) * exp(sel_inf_off(0, flt, sex, yr));
            Type val = 1.0 / (1.0 + exp(-slope * (x_val - inf)));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        // NonParametricIntegrable (13) under Time_varying_sel = "RandomWalk": the Ianelli
        // curve on the base coefficients plus the running sum of increments (years after
        // the fleet's start year). Only the realized curve is mean-centred -- the walk
        // carries the raw sum -- so with no increments every year is NonParametric's
        // curve. No cap; the increments are scored in ceattle.cpp.
        case SEL_CASE_NP_INTEGRABLE_WALK: {
          for(int bin = 0; bin < n_sel_bins; bin++){
            Type inc = (yr > flt_sel_start_yr(flt)) ? sel_coff_dev(flt, sex, bin, yr) : Type(0.0);
            np_unc(sex, bin, yr) = (yr > 0 ? np_unc(sex, bin, yr - 1) : Type(0.0)) + inc;
            non_par_sel(flt, sex, bin, yr) = sel_coff(flt, sex, bin) + np_unc(sex, bin, yr);
          }
          { vector<Type> base(n_sel_bins);
            for(int bin = 0; bin < n_sel_bins; bin++) base(bin) = non_par_sel(flt, sex, bin, yr);
            avg_sel(flt, sex, yr) = log_mean_exp(base); }
          for(int bin = n_sel_bins; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, n_sel_bins - 1, yr);
          }
          { vector<Type> ls(nbins);
            for(int bin = 0; bin < nbins; bin++) ls(bin) = non_par_sel(flt, sex, bin, yr);
            avgsel_tmp = log_mean_exp(ls); }
          for(int bin = 0; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) -= avgsel_tmp;
            log_non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            non_par_sel(flt, sex, bin, yr) = exp(non_par_sel(flt, sex, bin, yr));
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            else                 sel_at_age(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
          }
          break;
        }

        case 9: { // NonParametricPM (AMAK "pm"): random walk on the per-year-
                  // renormalized log-selectivity, then a flat bin-cap.
          // sel_coff = base coefficients (ages 0..n_sel_bins-1, the rest plateauing);
          // sel_coff_dev = RAW increments at the year they apply. The UNCAPPED centred
          // curve (np_unc) carries the walk; the realized curve is that held flat above
          // flt_sel_cap_bin and re-centred (mean(exp) = 1).
          // At or before the fleet's selectivity start year the curve is built from the
          // base coefficients alone rather than the running walk -- the AMAK convention
          // of setting a survey's base selectivity once at its start year.
          bool from_base = (yr <= flt_sel_start_yr(flt));
          for(int bin = 0; bin < nbins; bin++){
            Type prev = from_base ? sel_coff(flt, sex, (bin < n_sel_bins ? bin : n_sel_bins - 1))
                                  : np_unc(sex, (bin < n_sel_bins ? bin : n_sel_bins - 1), yr - 1);
            Type inc  = (bin < n_sel_bins) ? sel_coff_dev(flt, sex, bin, yr) : Type(0.0);
            np_unc(sex, bin, yr) = prev + inc;
          }
          for(int bin = n_sel_bins; bin < nbins; bin++) np_unc(sex, bin, yr) = np_unc(sex, n_sel_bins - 1, yr);
          // Bins below bin_first_selected (e.g. the acoustic-survey age-1) are held at 0
          // before each year's mean-centring, the AMAK convention for a fixed first-bin
          // log-selectivity. This keeps the excluded bin out of the shared normalization
          // of the selected ages and out of the deviate/curvature penalties.
          for(int bin = 0; bin < bin_first_selected(flt); bin++) np_unc(sex, bin, yr) = 0.0;
          { vector<Type> u(nbins);
            for(int bin = 0; bin < nbins; bin++) u(bin) = np_unc(sex, bin, yr);
            Type m = log_mean_exp(u);
            for(int bin = 0; bin < nbins; bin++) np_unc(sex, bin, yr) -= m; }
          { vector<Type> cl(nbins);
            for(int bin = 0; bin < nbins; bin++) cl(bin) = np_unc(sex, bin, yr);
            int cap = flt_sel_cap_bin(flt);
            if(cap >= 0) for(int bin = cap + 1; bin < nbins; bin++) cl(bin) = cl(cap);
            Type m2 = log_mean_exp(cl);
            for(int bin = 0; bin < nbins; bin++){
              log_non_par_sel(flt, sex, bin, yr) = cl(bin) - m2;
              non_par_sel(flt, sex, bin, yr) = exp(cl(bin) - m2);
              if (is_length_based) sel_at_length(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
              else                 sel_at_age(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            }
          }
          break;
        }

        case 13: // NonParametricIntegrable under "Off"/"IID"; penalties on the base
        case 2: // Non-parametric (Ianelli style)
          for(int bin = 0; bin < n_sel_bins; bin++) {
            non_par_sel(flt, sex, bin, yr) = sel_coff(flt, sex, bin) + sel_coff_dev(flt, sex, bin, yr);
          }
          // AMAK avgsel: the level of the estimated coefficients, which the mean-centring
          // below removes from the curve. The penalty on it is the only thing holding
          // that level (ceattle.cpp, JNLL_SEL_NONPARAM).
          { vector<Type> base(n_sel_bins);
            for(int bin = 0; bin < n_sel_bins; bin++) base(bin) = non_par_sel(flt, sex, bin, yr);
            avg_sel(flt, sex, yr) = log_mean_exp(base); }

          for(int bin = n_sel_bins; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, n_sel_bins - 1, yr);
          }

          { vector<Type> ls(nbins);
            for(int bin = 0; bin < nbins; bin++) ls(bin) = non_par_sel(flt, sex, bin, yr);
            avgsel_tmp = log_mean_exp(ls); }

          for(int bin = 0; bin < nbins; bin++) {
            non_par_sel(flt, sex, bin, yr) -= avgsel_tmp;
            log_non_par_sel(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            non_par_sel(flt, sex, bin, yr) = exp(non_par_sel(flt, sex, bin, yr));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
            else sel_at_age(flt, sex, bin, yr) = non_par_sel(flt, sex, bin, yr);
          }
          break;

        case 3: // Double Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? xmid(bin) : Type(bin + 1);
            Type slp1 = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                            + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
            Type inf1 = (sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr)
                         + sel_inf_off_nat(0, flt, sex, yr)) * exp(sel_inf_off(0, flt, sex, yr));
            Type slp2 = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr)
                            + sel_slp_off(1, flt, sex, yr)) + sel_slp_off_nat(1, flt, sex, yr);
            Type inf2 = (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr)
                         + sel_inf_off_nat(1, flt, sex, yr)) * exp(sel_inf_off(1, flt, sex, yr));
            Type val = (1.0 / (1.0 + exp(-slp1 * (x_val - inf1)))) * (1.0 - (1.0 / (1.0 + exp(-slp2 * (x_val - inf2)))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 4: // Descending Logistic
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? xmid(bin) : Type(bin + 1);
            Type slp2 = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr)
                            + sel_slp_off(1, flt, sex, yr)) + sel_slp_off_nat(1, flt, sex, yr);
            Type inf2 = (sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr)
                         + sel_inf_off_nat(1, flt, sex, yr)) * exp(sel_inf_off(1, flt, sex, yr));
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

            // Normalized per sex here, not in normalize_and_project_selectivity(), so
            // both sexes reach 1: Sel_norm_scope is inert on Hake, which cannot carry a
            // sex difference in level (issue #236).
            max_sel = -1e10;
            if (sel_norm_bin1(flt) >= 0 && sel_norm_bin2(flt) < 0) {
              max_sel = is_length_based ? sel_at_length(flt, sex, sel_norm_bin1(flt), yr) : sel_at_age(flt, sex, sel_norm_bin1(flt), yr);
            } else {
              for(int bin=0; bin<nbins; bin++) {
                Type val = is_length_based ? sel_at_length(flt, sex, bin, yr) : sel_at_age(flt, sex, bin, yr);
                max_sel = max2(max_sel, val);   // fold, don't branch on an AD Type
              }
            }
            for (int bin = 0; bin < nbins; bin++) {
              if (is_length_based) sel_at_length(flt, sex, bin, yr) = exp(sel_at_length(flt, sex, bin, yr) - max_sel);
              else sel_at_age(flt, sex, bin, yr) = exp(sel_at_age(flt, sex, bin, yr) - max_sel);
            }
          }
          break;

        case 8: { // DoubleNormal (4-param dome with a right-tail floor)
          // Reuses the logistic slots, each with its own annual deviate:
          //   sel_inf(0)     = peak bin (mode of selectivity)
          //   log_sel_slp(0) = log(sigma_ascending), the ascending limb SD
          //   log_sel_slp(1) = log(sigma_descending), the descending limb SD
          //   sel_inf(1)     = logit(right_floor), analogous to SS3 P6 (end_logit); a
          //                    floor of 0 is a full dome, 1 a logistic (ascending only).
          Type peak        = (sel_inf(0, flt, sex) + sel_inf_dev(0, flt, sex, yr)
                              + sel_inf_off_nat(0, flt, sex, yr)) * exp(sel_inf_off(0, flt, sex, yr));
          Type sigma_asc   = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                              + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
          Type sigma_desc  = exp(log_sel_slp(1, flt, sex) + log_sel_slp_dev(1, flt, sex, yr)
                              + sel_slp_off(1, flt, sex, yr)) + sel_slp_off_nat(1, flt, sex, yr);
          // The linkage offsets act on the LOGIT scale, inside the logistic, matching the
          // other inflection slots. That keeps right_floor in [0, 1] for any offset,
          // where scaling the probability itself could push it above 1.
          Type right_floor = 1.0 / (1.0 + exp(-((sel_inf(1, flt, sex) + sel_inf_dev(1, flt, sex, yr)
                              + sel_inf_off_nat(1, flt, sex, yr)) * exp(sel_inf_off(1, flt, sex, yr)))));
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val      = is_length_based ? xmid(bin) : Type(bin + 1);
            // Smooth logistic blend: ~0 left of peak, ~1 right of peak.
            Type w          = 1.0 / (1.0 + exp(-20.0 * (x_val - peak)));
            Type asc_gauss  = exp(-0.5 * pow((x_val - peak) / sigma_asc,  2.0));
            Type desc_gauss = exp(-0.5 * pow((x_val - peak) / sigma_desc, 2.0));
            Type val_desc   = right_floor + (1.0 - right_floor) * desc_gauss;
            Type val        = (1.0 - w) * asc_gauss + w * val_desc;
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else                 sel_at_age(flt, sex, bin, yr) = val;
          }
          break;
        }

        case 15: { // DoubleNormalSS3: SS3 size pattern 24 (see the form list above)
          // Each of the six takes P = (base + identity offset) x exp(log offset), so an
          // SS3 block replacement (identity link on a block factor) followed by
          // multiplicative annual deviates (log link) reproduces SS3's own order.
          vector<Type> P(6);
          for (int k = 0; k < 6; k++) {
            P(k) = (sel_dn6(k, flt, sex) + sel_dn6_off_nat(k, flt, sex, yr)) * exp(sel_dn6_off(k, flt, sex, yr));
          }
          Type peak      = P(0);
          Type upselex   = exp(P(2));
          Type downselex = exp(P(3));
          // SS3 anchors the ASCENDING limb at `startbin`: the first population bin that
          // reaches the first composition DATA bin, not the first population bin
          // (SS_readdata_330.tpl:1642-1646, used at SS_selex.tpl:142). With one grid
          // startbin is 0.
          int startbin = 0;
          if (is_length_based) {
            while (startbin < nbins - 1 && lengths_pop(sp, startbin) < lengths(sp, 0)) startbin++;
          }
          Type x_first   = is_length_based ? xmid(startbin) : Type(1);
          Type x_last    = is_length_based ? xmid(nbins - 1) : Type(nbins);
          Type peak2     = peak + binwidth + (0.99 * x_last - peak - binwidth) / (1.0 + exp(-P(1)));
          Type point1    = 1.0 / (1.0 + exp(-P(4)));
          Type point2    = 1.0 / (1.0 + exp(-P(5)));
          Type t1min     = exp(-square(x_first - peak) / upselex);
          Type t2min     = exp(-square(x_last - peak2) / downselex);
          for (int bin = startbin; bin < nbins; bin++) {
            Type x_val = is_length_based ? xmid(bin) : Type(bin + 1);
            Type t1    = x_val - peak;
            Type t2    = x_val - peak2;
            Type join1 = 1.0 / (1.0 + exp(-(20.0 * t1 / (1.0 + CppAD::abs(t1)))));
            Type join2 = 1.0 / (1.0 + exp(-(20.0 * t2 / (1.0 + CppAD::abs(t2)))));
            Type asc   = (sel_dn6_ends(flt, 0) == 1)
                         ? point1 + (1.0 - point1) * (exp(-square(t1) / upselex) - t1min) / (1.0 - t1min)
                         : exp(-square(t1) / upselex);
            Type dsc   = (sel_dn6_ends(flt, 1) == 1)
                         ? 1.0 + (point2 - 1.0) * (exp(-square(t2) / downselex) - 1.0) / (t2min - 1.0)
                         : exp(-square(t2) / downselex);
            Type val   = asc * (1.0 - join1) + join1 * ((1.0 - join2) + dsc * join2);
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else                 sel_at_age(flt, sex, bin, yr) = val;
          }
          // Below startbin SS3 does not leave the curve at zero: it scales the value at
          // startbin by (length / length_at_startbin)^2, a quadratic ramp
          // (SS_selex.tpl:197-201). Reached only where the data bins start above the
          // population bins.
          if (is_length_based && startbin > 0) {
            Type sel_start = sel_at_length(flt, sex, startbin, yr);
            Type x_start   = xmid(startbin);
            for (int bin = 0; bin < startbin; bin++) {
              Type x_val = xmid(bin);
              sel_at_length(flt, sex, bin, yr) = square(x_val / x_start) * sel_start;
            }
          }
          break;
        }

        case 6: // 2DAR1 (bin x year)
        case 7: // 3DAR1 (bin x year, plus a cohort term; conditional variance)
          for (int bin = 0; bin < nbins; bin++) {

            // Bins past N_sel_bins plateau on the last estimated coefficient
            int active_bin = (bin < n_sel_bins) ? bin : (n_sel_bins - 1);
            Type val = 1.0 / (1.0 + exp(-(sel_coff(flt, sex, active_bin) + sel_coff_dev(flt, sex, active_bin, yr))));

            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else sel_at_age(flt, sex, bin, yr) = val;
          }
          break;

        case 11: { // LogisticPM (ADMB AMAK "pm" bottom-trawl survey form)
          // 2-parameter logistic with MULTIPLICATIVE deviations on slope and inflection,
          // matching AMAK's
          //   sel = 1/(1 + exp(-slp*exp(slp_dev) * (age - a50*exp(a50_dev)))),
          // plus a FREE first-bin (age-1) log-selectivity independent of the logistic
          // (AMAK sel_age_one_bts*exp(sel_age_one_bts_dev)), held in the unused
          // descending-limb slots sel_inf(1)/sel_inf_dev(1). Never normalized: AMAK does
          // not renormalize the BTS curve and age-1 may exceed 1.
          // AMAK evaluates the logistic at age_vector(j) = j + 0.5 (mid-age), so the
          // age-based x is bin + 1.5, NOT bin + 1 as in the standard Logistic (case 1).
          // The 0.5 shift cannot be folded into a50 because the inflection deviate is
          // multiplicative; a log-link offset therefore rides inside that same exp, so
          // log-link still means multiplicative-on-natural as it does elsewhere.
          Type slope = exp(log_sel_slp(0, flt, sex) + log_sel_slp_dev(0, flt, sex, yr)
                           + sel_slp_off(0, flt, sex, yr)) + sel_slp_off_nat(0, flt, sex, yr);
          Type inf   = (sel_inf(0, flt, sex) + sel_inf_off_nat(0, flt, sex, yr))
                         * exp(sel_inf_dev(0, flt, sex, yr) + sel_inf_off(0, flt, sex, yr));
          for (int bin = 0; bin < nbins; bin++) {
            Type x_val = is_length_based ? xmid(bin) : Type(bin + 1.5);
            Type val = 1.0 / (1.0 + exp(-slope * (x_val - inf)));
            if (is_length_based) sel_at_length(flt, sex, bin, yr) = val;
            else                 sel_at_age(flt, sex, bin, yr) = val;
          }
          // Free first-bin (age-1) log-selectivity override
          Type log_s1 = (sel_inf(1, flt, sex) + sel_inf_off_nat(1, flt, sex, yr))
                          * exp(sel_inf_dev(1, flt, sex, yr) + sel_inf_off(1, flt, sex, yr));
          if (is_length_based) sel_at_length(flt, sex, 0, yr) = exp(log_s1);
          else                 sel_at_age(flt, sex, 0, yr) = exp(log_s1);
          break;
        }
        }

      } // End sex
    } // End year

    // --- 2b. PER-SEX APICAL HEIGHT ---
    // One sex's whole curve times exp(log_sel_apical + log-link offset), plus the
    // identity-link offset, the convention the log-scale slopes use. Applied after each
    // form's own centring or normalization and before the shared normalizer, which must
    // pool across sexes (AcrossSexes) or be off for the offset to survive;
    // .check_sel_apical_rows() refuses the rest.
    if (sel_apical_on == 1) {
      for (int yr = 0; yr < nyrs_hind; yr++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          Type scale = exp(log_sel_apical(flt, sex) + sel_apical_off(flt, sex, yr))
                       + sel_apical_off_nat(flt, sex, yr);
          for (int bin = 0; bin < nbins; bin++) {
            if (is_length_based) sel_at_length(flt, sex, bin, yr) *= scale;
            else                 sel_at_age(flt, sex, bin, yr) *= scale;
          }
        }
      }
    }

    // --- 3. NORMALIZATION & PROJECTION ---
    normalize_and_project_selectivity(
      flt, nyrs_hind, nyrs, flt_spp, flt_sel_type, flt_sel_dim, bin_first_selected, nages, nlengths_pop, nsex,
      sel_norm_bin1, sel_norm_bin2, sel_norm_scope,
      is_length_based ? sel_at_length : sel_at_age
    );

    // --- 4. CONVERT LENGTH TO AGE (If necessary) ---
    if (is_length_based) {
      for (int yr = 0; yr < nyrs; yr++) {
        for (int sex = 0; sex < nsex(sp); sex++) {
          int wtind = nspp * 2 + flt;
          convert_length_selectivity(
            flt, sp, sex, yr, wtind, nages, nlengths_pop,
            sel_at_length, growth_matrix, sel_at_age
          );
        }
      }
    }
  } // End fleet
}

#endif
