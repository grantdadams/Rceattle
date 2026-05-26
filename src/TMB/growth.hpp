#ifndef GROWTH_HPP
#define GROWTH_HPP

/**
 * @brief Integrated Growth, Size-Transition, and Weight-at-Age Module for Month = 0.
 *
 * Computes Jan-1 mean length-at-age, the age->length probability matrix, and
 * integrated weight-at-age for a single species.
 *
 * @section math_models Mathematical Models:
 * 1. Mean Length-at-Age ($L_a$) — anchored at $L(a_{L1}) = l_1$:
 *    - Von Bertalanffy (Model 1): $L_a = L_{\infty} + (l_1 - L_{\infty}) \cdot e^{-K(a - a_{L1})}$
 *    - Richards (Model 2): $L_a = [L_{\infty}^m + (l_1^m - L_{\infty}^m) \cdot e^{-K(a - a_{L1})}]^{1/m}$
 *    - For year > 0, a cohort recursion advances $L_{a, y}$ from $L_{a-1, y-1}$
 *      using lag-year parameters. At the cohort boundary (current_age ==
 *      age_L1_ceil) the closed-form anchor at $l_1$ is used; `age_L1_safe`
 *      keeps this anchor at $l_1$ for both minage > 0 and minage = 0.
 *    - Linear ramp from Lmin_sp at age 0 to $l_1$ at age_L1 applies when
 *      current_age <= age_L1 (only reachable for minage > 0).
 *
 * 2. Weight-at-Age ($W_a$):
 *    Integrated across length to account for Jensen's Inequality:
 *    $W_a = \sum_{ln} P(ln | a) \cdot \alpha \cdot L_{mid}^{\beta}$
 *    Bin midpoints are computed per-bin to support non-uniform length bins;
 *    the length plus-group is extended by half the final interior bin width.
 *
 * @section logic Biological Logic:
 * - **Temporal Resolution**: Jan-1 (month = 0).
 * - **Plus-Age Group**: Oldest age class is corrected via a static
 *   $\exp(-0.2 \cdot a)$-weighted mean of $[current\_size, ..., L_{\infty}]$
 *   over $a = 0..nages$. This is the WHAM static analogue of SS3's
 *   N-at-age-weighted recruitment correction at the season transition.
 * - **SD-at-Age**: For current_age <= age_L1, SD = $e^{sd_0}$. Otherwise
 *   linear interpolation in length between SD($l_1$) = $e^{sd_0}$ and
 *   SD($L_{\infty}$) = $e^{sd_1}$. The plus group uses the same
 *   interpolation (no special-case) for continuity.
 * - **Size Transition**: Converts mean length and SD into a probability
 *   matrix $P(\text{Length} | \text{Age})$ via `pnorm`. First length bin is
 *   a minus-group; last length bin is a plus-group.
 *
 * @param wtind Weight index slot to write into.
 * @param sp Species index.
 * @param nyrs Number of years in the simulation/model.
 * @param nsex/nages/nlengths Vectors of species-specific counts per dimension.
 * @param minage Vector of minimum modeled ages (0 supported).
 * @param growth_model Per-species growth model selector (1 = VBGF, 2 = Richards).
 * @param lengths Matrix of length bin boundaries (nspp x length_bins). Bin
 *                widths may be non-uniform.
 * @param growth_parameters 4D array of time-varying growth parameters (K, L1, Linf, m).
 * @param growth_log_sd 3D array of log-scale SDs (index 0 = SD at L1, index 1 = SD at Linf).
 * @param weight_length_pars Matrix of length-weight parameters ($\alpha$, $\beta$).
 *
 * @note Outputs are passed by reference.
 * @param length_hat [Output] 4D array filled with mean length-at-age.
 * @param growth_matrix [Output] 5D array filled with the growth transition matrix.
 * @param weight_hat [Output] 4D array filled with integrated weight-at-age.
 */
template<class Type>
void estimate_growth(
    int wtind,
    int sp,
    int nspp,
    int nyrs,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  nlengths_pop,
    const vector<int>&  minage,
    const vector<Type>& growth_age_L1,
    const vector<int>&  growth_model,
    matrix<Type> lengths,
    matrix<Type> lengths_pop,
    array<Type>& growth_parameters,
    array<Type>& growth_log_sd,
    matrix<Type> weight_length_pars,
    array<Type> &length_hat,         // Modified by reference
    array<Type> &growth_matrix,      // DATA-bin ALK (length-comp / CAAL)
    array<Type> &growth_matrix_pop,  // POP-bin ALK (sel-to-age, SS3 parity)
    array<Type> &weight_hat          // Modified by reference
) {

  // Initialize output and temporary storage
  array<Type> length_sd(nsex(sp), nages(sp), nyrs); length_sd.setZero();         // SD in length-at-age


  // Calculate mean-length, SD, and growth matrix, for all years:
  // lengths_pop is the FINE integration grid (SS3 parity: 1cm pop bins).
  // lengths is the DATA grid (e.g. 5cm) used for length-comp likelihood.
  // ALK + WAA are computed on lengths_pop, then ALK is aggregated to
  // lengths for the growth_matrix output that downstream length-comp /
  // CAAL likelihood paths consume.
  Type Fac1, Fac2, Slope, b_len, last_linear, current_age;

  // Anchor Lmin/Lmax on the POP grid (SS3 parity for the age-0 floor and
  // the length plus-group at the top of the pop bins). When the user has
  // not supplied a pop grid, fit_mod() falls back to lengths_pop = lengths
  // so this preserves prior behavior.
  Type Lmin_sp = lengths_pop(sp, 0);
  Type Lmax_sp = lengths_pop(sp, nlengths_pop(sp) - 1);
  // age_L1 is the VB anchor age (= age at which `l1` is the length). Read
  // from data_list$growth_age_L1[sp] (= SS3's Growth_Age_for_L1 ctl input).
  // R-side fit_mod() resolves the default to max(0.5, minage[sp]) so models
  // with minage = 0 get an SS3-style half-year anchor and minage >= 1 stays
  // backwards-compatible.
  Type age_L1 = growth_age_L1(sp);
  // age_L1_ceil is compared to current_age (the slot age in years) to detect
  // the youngest VB-relevant slot. That slot needs the closed-form boundary
  // formula (not the cohort recursion which would index age-1 = -1 at slot 0
  // when minage = 0, or use a stale value at slot 1 when the anchor is at
  // age 0.5). For minage = 0, the youngest VB slot is C++ age = 1
  // (current_age = 1, anchor at 0.5). For minage >= 1, slot 0 itself
  // (current_age = minage). Either way: age_L1_ceil = max(1, minage).
  Type age_L1_ceil = (Type(minage(sp)) >= Type(1)) ? Type(minage(sp)) : Type(1);

  // When minage = 0, the linear ramp from Lmin_sp at age 0 to l1 at age_L1
  // becomes degenerate (divide by zero in b_len) AND the ramp branch should
  // never execute (current_age >= 1 always > age_L1 = 0). Replace the ramp
  // anchor with a small positive value so b_len is finite; the value is
  // unused because the if(current_age <= age_L1) branch is unreachable.
  Type age_L1_safe = age_L1;
  if (age_L1 <= Type(0)) age_L1_safe = Type(1);  // safe denominator only

  for(int sex = 0; sex < nsex(sp); sex++) {
    for(int yr = 0; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++) {

        // Parameters for parametric growth
        Type kappa = growth_parameters(sp, sex, yr, 0);
        Type l1 = growth_parameters(sp, sex, yr, 1);
        Type linf = growth_parameters(sp, sex, yr, 2);
        Type m = growth_parameters(sp, sex, yr, 3);
        // current_age is the slot's actual age (start of year). Slot index
        // age (0-based) holds the (age + minage) cohort. The historical
        // `age + 1.0` shifted everything by one year at any minage != 1
        // (gave end-of-year length instead of start). Note this is the SLOT
        // age, NOT the VB anchor (which is age_L1, defined above).
        current_age = Type(age) + Type(minage(sp));

        // 1. Calculate Mean Length at Age ---
        switch(growth_model(sp)) {
        case 1: // Von Bertalanffy

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          // Age < minage
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = linf + (l1 - linf) * exp(-kappa * (current_age - age_L1));
            } else { // Yr > 0
              if(current_age == age_L1_ceil) {
                // age_L1_safe ensures last_linear == l1 in both minage > 0
                // (where age_L1_safe == age_L1) and minage = 0 (where
                // age_L1 = 0 would otherwise collapse this to Lmin_sp).
                last_linear = Lmin_sp + b_len * age_L1_safe;
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
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          if((current_age) <= age_L1) {
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else {
            if(yr == 0) {
              length_hat(wtind,  sex, age, yr) = pow(pow(linf, m) + (pow(l1, m) - pow(linf, m)) * exp(-kappa * (current_age - age_L1)), 1 / m);
            } else { // Yr > 0
              if(current_age == age_L1_ceil) {
                // age_L1_safe ensures last_linear == l1 in both minage > 0
                // and minage = 0; see VBGF branch above for full reasoning.
                last_linear = Lmin_sp + b_len * age_L1_safe;
                length_hat(wtind,  sex, age, yr) = pow(pow(last_linear, m) + (pow(last_linear, m) - pow(linf, m)) * (exp(-kappa * (current_age - age_L1)) - 1.0), 1 / m);
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
        // CURRENT (WHAM-style static approximation): assume the plus-group
        // age distribution follows exp(-0.2 * a) and linearly interpolate
        // length from the youngest plus-group cohort (current_size) up to
        // Linf as fish age 0..nages within the plus group. Cheap, no N
        // dependence, decoupled from the N propagation -- but ~0.2-1%
        // off SS3 in tested models because (a) the decay rate is fixed,
        // (b) lengths are interpolated linearly instead of advanced via
        // VB, and (c) the cohort weights don't reflect the actual
        // hindcast N-at-plus-group.
        if(growth_model(sp) < 3 && age == (nages(sp) - 1)) {
          Type current_size = length_hat(wtind,  sex, age, yr);
          Type temp_n = 0, temp_sum = 0, weight_a = 1.0;
          Type diff = linf - current_size;
          for(int a = 0; a <= nages(sp); a++) {
            temp_sum += weight_a * (current_size + (Type(a) / Type(nages(sp))) * diff);
            temp_n += weight_a;
            weight_a *= exp(-0.2); //FIXME: update mortality?
          }
          length_hat(wtind,  sex, age, yr) = temp_sum / temp_n;
        }

        // SS3 reference (dynamic N-weighted plus-group correction) ----
        // SS3's morph-tracking would compute the plus-group length as the
        // N-weighted mean over the cohorts that have been promoted into
        // the plus group. Each year, the N at plus-age picks up the cohort
        // freshly aged in (at the just-derived VB length L_new) plus the
        // surviving previous plus-group (at the advanced previous-year
        // length L_old). Both are weighted by their N and survival to
        // produce the mean Jan-1 length used here.
        //
        // To enable below, calculate_weight() needs N_at_age and S
        // (survival) passed in; and the loop must run with knowledge of
        // the previous year's plus-group length (which we already store
        // in length_hat at (wtind, sex, age, yr - 1)).
        //
        // if(growth_model(sp) < 3 && age == (nages(sp) - 1) && yr > 0) {
        //   // VB length the just-promoted cohort enters the plus group at
        //   Type L_new = linf - (linf - l1) * exp(-kappa * (current_age - age_L1));
        //   // Previous-year plus-group length advanced one year of VB
        //   Type L_old_prev = length_hat(wtind, sex, age, yr - 1);
        //   Type L_old_adv  = linf + (L_old_prev - linf) * exp(-kappa);
        //   // N entering vs surviving (read from N_at_age + survival arrays)
        //   Type N_in  = N_at_age(sp, sex, age - 1, yr - 1) *
        //                exp(-Z_at_age(sp, sex, age - 1, yr - 1));
        //   Type N_old = N_at_age(sp, sex, age,     yr - 1) *
        //                exp(-Z_at_age(sp, sex, age,     yr - 1));
        //   Type N_tot = N_in + N_old;
        //   length_hat(wtind, sex, age, yr) =
        //     (N_in * L_new + N_old * L_old_adv) / N_tot;
        // }
        // For yr == 0 (initial), SS3 iterates the recursion to equilibrium
        // under M = M1; in single-species mode that converges quickly. A
        // static approximation: use the VB closed-form at age nages-1 plus
        // a few extra ages weighted by exp(-M).

        // 3. Calculate SD (Integrated) ---
        // Length-based linear interpolation: SD(l1) = sd0, SD(linf) = sd1.
        // Applied uniformly above age_L1 (including the plus group) to avoid
        // the discontinuity that a special-case plus-group rule would create.
        if(growth_model(sp) < 3) {
          if((current_age) <= age_L1) {
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0));
          } else {
            Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0))) / (linf - l1);
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0)) + Slope * (length_hat(wtind,  sex, age, yr) - l1);
          }

          // Free parameters
          if(growth_model(sp) == 3) {
            // Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0)))/(length_hat(wtind,  sex, nages(sp)-1, yr) - length_hat(wtind,  sex, 0, yr));
            // length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0) + Slope * (length_hat(wtind,  sex, age, yr) - length_hat(wtind,  sex, 0, yr));
          }
        }

        // 4. Build Growth Matrix & Weight-at-Age ---
        // ALK + WAA integrated on the POP grid (SS3 parity); the pop-bin
        // probabilities are aggregated to the DATA grid and stored in
        // growth_matrix (consumed by length-comp / CAAL likelihood code
        // downstream). See estimate_growth_within_yr() for the same
        // pattern at non-zero fracyr.
        Type expected_weight = 0.0;
        for(int ln = 0; ln < nlengths_pop(sp); ln++) {
          Type prob;
          if(ln == 0) {
            Fac1 = (Lmin_sp + lengths_pop(sp, 1) - lengths_pop(sp, 0) - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1);
          } else if(ln == (nlengths_pop(sp) - 1)) {
            Fac1 = (Lmax_sp - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            prob = 1.0 - pnorm(Fac1);
          } else {
            Fac1 = (lengths_pop(sp, ln + 1) - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            Fac2 = (lengths_pop(sp, ln) - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1) - pnorm(Fac2);
          }
          // Store pop-bin ALK (consumed by SS3-parity sel-to-age convolution)
          growth_matrix_pop(wtind, sex, age, ln, yr) = prob;

          // Pop-bin midpoint for WAA Jensen's integration; plus-group
          // extends by half the last pop-bin width (supports non-uniform).
          Type lenmid;
          if(ln < nlengths_pop(sp) - 1) {
            lenmid = (lengths_pop(sp, ln) + lengths_pop(sp, ln + 1)) / Type(2.0);
          } else {
            Type last_width = lengths_pop(sp, nlengths_pop(sp) - 1) - lengths_pop(sp, nlengths_pop(sp) - 2);
            lenmid = lengths_pop(sp, nlengths_pop(sp) - 1) + last_width / Type(2.0);
          }
          expected_weight += prob * weight_length_pars(sp, 0) * pow(lenmid, weight_length_pars(sp, 1));
        }
        weight_hat(wtind, sex, age, yr) = expected_weight;

        // Aggregate pop-bin probs to data-bin growth_matrix (SS3 rule:
        // minus-group below first data edge, plus-group at/above last
        // edge; interior by [edge_k, edge_{k+1}) containing pop_lo).
        for(int ln = 0; ln < nlengths(sp); ln++) growth_matrix(wtind, sex, age, ln, yr) = Type(0.0);
        Type data_lo0 = lengths(sp, 0);
        for(int lp = 0; lp < nlengths_pop(sp); lp++) {
          Type pop_lo = lengths_pop(sp, lp);
          int target;
          if(pop_lo < data_lo0) {
            target = 0;
          } else {
            target = 0;
            for(int k = 0; k < nlengths(sp); k++) {
              if(lengths(sp, k) <= pop_lo) target = k;
              else break;
            }
          }
          growth_matrix(wtind, sex, age, target, yr) += growth_matrix_pop(wtind, sex, age, lp, yr);
        }
      } // age
    } // yr
  } // sex
}



/**
 * @brief Integrated Predator Growth, Size-Transition, and Weight-at-Age Module at Month X.
 *
 * @section math_models Mathematical Models:
 * 1. Mean Length-at-Age ($L_a$):
 *    Advances the Jan-1 length stored at id_pop forward by fracyr of growth.
 *    - Von Bertalanffy (Model 1): $L_a = L_{\infty} + (L_{a, jan1} - L_{\infty}) \cdot e^{-K \cdot fracyr}$
 *    - Richards (Model 2): $L_a = [L_{\infty}^m + (L_{a, jan1}^m - L_{\infty}^m) \cdot e^{-K \cdot fracyr}]^{1/m}$
 *    - Linear growth is applied for ages < minage (current_age <= age_L1).
 *
 * 2. Weight-at-Age ($W_a$):
 *    Calculated via integration across the length distribution to account for Jensen's Inequality:
 *    $W_a = \sum_{ln} P(ln | a) \cdot \alpha \cdot L_{mid}^{\beta}$
 *
 * @section logic Biological Logic:
 * - **Temporal Resolution**: Incorporates `fracyr` to allow within-year (seasonal)
 *   growth and differentiability for time-varying parameters.
 * - **Plus Group**: Advanced by within-year growth identically to other ages
 *   (SS3 convention). The recruitment-weighted-mean correction for fish promoted
 *   into the plus group is applied at every year boundary by `estimate_growth()`
 *   (month 0), so `id_pop` already carries the corrected Jan-1 length.
 * - **SD-at-Age**: Length-based linear interpolation between `SD(l1)` and
 *   `SD(linf)`, applied uniformly above `age_L1` (including the plus group).
 * - **Size Transition**: Converts mean length and SD into a probability matrix
 *   $P(\text{Length} | \text{Age})$ using a cumulative normal distribution (`pnorm`).
 *   The first length bin is a minus-group on length; the last is a plus-group.
 *
 * @param wtind Weight index for population/fleet.
 * @param id_pop Index for population (Jan-1) weight-at-age; the within-year
 *               growth advance is anchored on `length_hat(id_pop, ...)`.
 * @param sp Species index.
 * @param fracyr Fraction of the year (0.0 - 1.0) for within-year growth.
 * @param nyrs Number of years in the simulation/model.
 * @param nsex/nages/nlengths Vectors of species-specific counts per dimension.
 * @param lengths Matrix of length bin boundaries (nspp x length_bins). Bin
 *                widths may be non-uniform.
 * @param growth_parameters 4D array of time-varying growth parameters (K, L1, Linf, m).
 * @param growth_log_sd 3D array of log-scale SDs (index 0 = SD at L1, index 1 = SD at Linf).
 * @param weight_length_pars Matrix of length-weight parameters ($\alpha$, $\beta$).
 *
 * @note Outputs are passed by reference.
 * @param length_hat [Output] 4D array filled with mean length-at-age.
 * @param growth_matrix [Output] 5D array filled with the growth transition matrix.
 * @param weight_hat [Output] 4D array filled with integrated weight-at-age.
 */
template<class Type>
void estimate_growth_within_yr(
    int wtind,
    int id_pop,
    int sp,
    Type fracyr,
    int nspp,
    int nyrs,
    const vector<int>&  nsex,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  nlengths_pop,
    const vector<int>&  minage,
    const vector<Type>& growth_age_L1,
    const vector<int>&  growth_model,
    matrix<Type> lengths,
    matrix<Type> lengths_pop,
    array<Type> growth_parameters,
    array<Type> growth_log_sd,
    matrix<Type> weight_length_pars,
    array<Type> &length_hat,         // Modified by reference
    array<Type> &growth_matrix,      // DATA-bin ALK
    array<Type> &growth_matrix_pop,  // POP-bin ALK (SS3 parity)
    array<Type> &weight_hat          // Modified by reference
) {

  // Initialize output and temporary storage
  array<Type> length_sd(nsex(sp), nages(sp), nyrs); length_sd.setZero();         // SD in length-at-age


  // Calculate mean-length, SD, and growth matrix, for all years:
  // ALK + WAA computed on lengths_pop (fine SS3-style grid); growth_matrix
  // is filled on the data-bin grid by aggregation -- see estimate_growth()
  // for the same pattern at Jan 1.
  Type Fac1, Fac2, Slope, b_len, current_age;

  Type Lmin_sp = lengths_pop(sp, 0);
  Type Lmax_sp = lengths_pop(sp, nlengths_pop(sp) - 1);
  // age_L1 is the VB anchor age (= age at which `l1` is the length). Read
  // from data_list$growth_age_L1[sp] (= SS3's Growth_Age_for_L1 ctl input).
  // R-side fit_mod() resolves the default to max(0.5, minage[sp]) so models
  // with minage = 0 get an SS3-style half-year anchor and minage >= 1 stays
  // backwards-compatible.
  Type age_L1 = growth_age_L1(sp);

  // Safe denominator for the linear ramp slope when minage = 0. See month=0
  // estimate_growth() for full context.
  Type age_L1_safe = age_L1;
  if (age_L1 <= Type(0)) age_L1_safe = Type(1);

  for(int sex = 0; sex < nsex(sp); sex++) {
    for(int yr = 0; yr < nyrs; yr++) {
      for(int age = 0; age < nages(sp); age++) {

        // Parameters for parametric growth
        Type kappa = growth_parameters(sp, sex, yr, 0);
        Type l1 = growth_parameters(sp, sex, yr, 1);
        Type linf = growth_parameters(sp, sex, yr, 2);
        Type m = growth_parameters(sp, sex, yr, 3);
        // See month=0 overload: slot k = age (k - 1 + minage); current_age is
        // slot age (NOT the VB anchor) offset by fracyr. Was
        // `age + 1.0 + fracyr`, which is off by one when minage != 1.
        current_age = Type(age) + Type(minage(sp)) + fracyr;

        // 1. Calculate Mean Length at Age ---
        switch(growth_model(sp)) {
        case 1: // Von Bertalanffy

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          // Age < minage
          // The plus group is advanced through within-year growth like every
          // other age (SS3 convention). The recruitment-weighted-mean
          // correction for fish promoted into the plus group is reapplied at
          // every year boundary by estimate_growth(),
          // so id_pop already carries the corrected Jan-1 length.
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else if((current_age - fracyr) < age_L1) {
            // Jan-1 length for this slot was in the ramp (not yet on the
            // VB curve), so the cohort-style fracyr advance would carry
            // forward a ramp value as if it were on VB and overshoot. Use
            // closed-form VB anchored at (l1, age_L1) instead. Only fires
            // when minage = 0 and the slot's integer age is below age_L1.
            length_hat(wtind, sex, age, yr) = linf + (l1 - linf) * exp(-kappa * (current_age - age_L1));
          } else { // Growth curve (incl. plus group)
            length_hat(wtind,  sex, age, yr) = length_hat(id_pop,  sex, age, yr) + (length_hat(id_pop,  sex, age, yr) - linf) * (exp(-kappa * fracyr) - 1.0); // Add fracyr growth
          }
          break;

        case 2: // Richards

          // Slope from Lmin to L1
          b_len = (l1 - Lmin_sp) / age_L1_safe;  // age_L1_safe == age_L1 when minage > 0; == 1 when minage = 0 (ramp branch unreachable)

          // Age < minage (plus group advances via Richards growth — see VBGF
          // branch above for the year-boundary correction reasoning).
          if((current_age) <= age_L1){
            length_hat(wtind,  sex, age, yr) = Lmin_sp + b_len * (current_age);
          } else if((current_age - fracyr) < age_L1) {
            // See VBGF branch above: ramp Jan-1 cannot be cohort-advanced.
            length_hat(wtind, sex, age, yr) = pow(pow(linf, m) + (pow(l1, m) - pow(linf, m)) * exp(-kappa * (current_age - age_L1)), 1 / m);
          } else {
            length_hat(wtind,  sex, age, yr) = pow(pow(length_hat(id_pop,  sex, age, yr), m) + (pow(length_hat(id_pop,  sex, age, yr), m) - pow(linf, m)) * (exp(-kappa * fracyr) - 1.0), 1 / m); // Add fracyr growth
          }
          break;

        case 3: // Non-parametric (Free parameters)
          error("Non-parametric growth not yet implemented");
          // length_hat(wtind,  sex, age, yr) = exp(length_par(sp, sex, age) + length_par_re(sp, sex, age, yr));
          break;

        default:
          error("Invalid 'growth_model");
        } // Growth_model switch

        // 2. Calculate SD (Integrated) ---
        // Length-based linear interpolation: SD(l1) = sd0, SD(linf) = sd1.
        // Applied uniformly above age_L1 (including the plus group).
        if(growth_model(sp) < 3) {
          if((current_age) <= age_L1) {
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0));
          } else {
            Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0))) / (linf - l1);
            length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0)) + Slope * (length_hat(wtind,  sex, age, yr) - l1);
          }

          // Free parameters
          if(growth_model(sp) == 3) {
            // Slope = (exp(growth_log_sd(sp, sex, 1)) - exp(growth_log_sd(sp, sex, 0)))/(length_hat(wtind,  sex, nages(sp)-1, yr) - length_hat(wtind,  sex, 0, yr));
            // length_sd(sex, age, yr) = exp(growth_log_sd(sp, sex, 0) + Slope * (length_hat(wtind,  sex, age, yr) - length_hat(wtind,  sex, 0, yr));
          }
        }

        // 4. Build Growth Matrix & Weight-at-Age (pop-grid integration,
        //    aggregated to data-bin growth_matrix). See estimate_growth()
        //    for the algorithm.
        Type expected_weight = 0.0;
        for(int ln = 0; ln < nlengths_pop(sp); ln++) {
          Type prob;
          if(ln == 0) {
            Fac1 = (Lmin_sp + lengths_pop(sp, 1) - lengths_pop(sp, 0) - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1);
          } else if(ln == (nlengths_pop(sp) - 1)) {
            Fac1 = (Lmax_sp - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            prob = 1.0 - pnorm(Fac1);
          } else {
            Fac1 = (lengths_pop(sp, ln + 1) - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            Fac2 = (lengths_pop(sp, ln) - length_hat(wtind, sex, age, yr)) / length_sd(sex, age, yr);
            prob = pnorm(Fac1) - pnorm(Fac2);
          }
          growth_matrix_pop(wtind, sex, age, ln, yr) = prob;

          Type lenmid;
          if(ln < nlengths_pop(sp) - 1) {
            lenmid = (lengths_pop(sp, ln) + lengths_pop(sp, ln + 1)) / Type(2.0);
          } else {
            Type last_width = lengths_pop(sp, nlengths_pop(sp) - 1) - lengths_pop(sp, nlengths_pop(sp) - 2);
            lenmid = lengths_pop(sp, nlengths_pop(sp) - 1) + last_width / Type(2.0);
          }
          expected_weight += prob * weight_length_pars(sp, 0) * pow(lenmid, weight_length_pars(sp, 1));
        }
        weight_hat(wtind, sex, age, yr) = expected_weight;

        // Aggregate pop-bin probs to data-bin growth_matrix.
        for(int ln = 0; ln < nlengths(sp); ln++) growth_matrix(wtind, sex, age, ln, yr) = Type(0.0);
        Type data_lo0 = lengths(sp, 0);
        for(int lp = 0; lp < nlengths_pop(sp); lp++) {
          Type pop_lo = lengths_pop(sp, lp);
          int target;
          if(pop_lo < data_lo0) {
            target = 0;
          } else {
            target = 0;
            for(int k = 0; k < nlengths(sp); k++) {
              if(lengths(sp, k) <= pop_lo) target = k;
              else break;
            }
          }
          growth_matrix(wtind, sex, age, target, yr) += growth_matrix_pop(wtind, sex, age, lp, yr);
        }
      } // age
    } // yr
  } // sex
}


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
 * [Other parameters for estimate_growth: lengths, nlengths, minage,
 * growth_parameters, growth_log_sd, weight_length_pars]
 */
template <class Type>
void calculate_weight(
    array<Type> &weight_hat,         // Modified by reference
    array<Type> &length_hat,         // Modified by reference
    array<Type> &growth_matrix,      // DATA-bin ALK
    array<Type> &growth_matrix_pop,  // POP-bin ALK (SS3 parity)
    array<Type> weight_obs,
    const vector<int>&  growth_model,
    const int& nspp,
    const int& nyrs,
    const int& nyrs_hind,
    const int& n_flt,
    const vector<int>&  flt_spp,
    vector<Type> flt_month,
    const vector<int>&  nsex,
    const vector<int>&  minage,
    const vector<Type>& growth_age_L1,
    const vector<int>&  nages,
    const vector<int>&  nlengths,
    const vector<int>&  nlengths_pop,
    const vector<int>&  pop_wt_index,
    const vector<int>&  ssb_wt_index,
    const vector<int>&  flt_wt_index,
    vector<Type> spawn_month,
    matrix<Type> lengths,
    matrix<Type> lengths_pop,
    array<Type> growth_parameters,
    array<Type> growth_log_sd,
    matrix<Type> weight_length_pars
) {
  int yr_ind;
  int wt_idx_pop;
  int wt_idx_ssb;
  int wt_idx_flt;

  // 1. POPULATION WEIGHT-AT-AGE
  for (int sp = 0; sp < nspp; sp++) {

    wt_idx_pop = 2 * sp;
    wt_idx_ssb = 2 * sp + 1;

    // -- 1.1. Empirical weight-at-age
    if (growth_model(sp) == 0) {
      for (int sex = 0; sex < nsex(sp); sex++) {
        for (int age = 0; age < nages(sp); age++) {
          for (int yr = 0; yr < nyrs; yr++) {

            // Handle projection logic
            yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);

            // Biomass weight
            weight_hat(wt_idx_pop, sex, age, yr) = weight_obs(pop_wt_index(sp), sex, age, yr_ind);

            // SSB weight
            weight_hat(wt_idx_ssb, sex, age, yr) = weight_obs(ssb_wt_index(sp), sex, age, yr_ind);
          }
        }
      }
    }

    // -- 1.2. Estimated growth
    if (growth_model(sp) > 0) {
      // Biomass weight (beginning of year / month 0)
      estimate_growth(
        wt_idx_pop,
        sp,
        nspp,
        nyrs,
        nsex,
        nages,
        nlengths,
        nlengths_pop,
        minage,
        growth_age_L1,
        growth_model,
        lengths,
        lengths_pop,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        growth_matrix_pop,  // Pass by reference (POP-bin)
        weight_hat      // Pass by reference
      );

      // SSB weight (at month of spawning)
      estimate_growth_within_yr(
        wt_idx_ssb,
        wt_idx_pop,
        sp,
        spawn_month(sp) / Type(12.0),
        nspp,
        nyrs,
        nsex,
        nages,
        nlengths,
        nlengths_pop,
        minage,
        growth_age_L1,
        growth_model,
        lengths,
        lengths_pop,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        growth_matrix_pop,  // Pass by reference (POP-bin)
        weight_hat      // Pass by reference
      );
    }
  }

  // 2. FLEET WEIGHT-AT-AGE
  for (int flt = 0; flt < n_flt; flt++) {
    int sp = flt_spp(flt);
    Type mo = flt_month(flt);
    wt_idx_pop = 2 * sp;
    wt_idx_flt = nspp * 2 + flt;

    // -- 2.1. Empirical weight-at-age
    if (growth_model(sp) == 0) {
      for (int sex = 0; sex < nsex(sp); sex++) {
        for (int age = 0; age < nages(sp); age++) {
          for (int yr = 0; yr < nyrs; yr++) {

            yr_ind = (yr < nyrs_hind) ? yr : (nyrs_hind - 1);
            weight_hat(wt_idx_flt, sex, age, yr) = weight_obs(flt_wt_index(flt), sex, age, yr_ind);
          }
        }
      }
    }

    // -- 2.2. Estimated growth
    if (growth_model(sp) > 0) {

      estimate_growth_within_yr(
        wt_idx_flt,
        wt_idx_pop,
        sp,
        mo / Type(12.0),
        nspp,
        nyrs,
        nsex,
        nages,
        nlengths,
        nlengths_pop,
        minage,
        growth_age_L1,
        growth_model,
        lengths,
        lengths_pop,
        growth_parameters,
        growth_log_sd,
        weight_length_pars,
        length_hat,     // Pass by reference
        growth_matrix,  // Pass by reference
        growth_matrix_pop,  // Pass by reference (POP-bin)
        weight_hat      // Pass by reference
      );
    }
  }
}



// ------------------------------------------------------------------------- //
// TODO                                                                      //
// ------------------------------------------------------------------------- //
// add growth parameters to build param and build map (AR1 and variance)
// define age_L1 and age_L1 ceiling

#endif

