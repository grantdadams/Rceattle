#include <TMB.hpp>


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
void get_growth_matrix(int sp, Type fracyr, int nspp, int nyrs,
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
          expected_weight += growth_matrix(wtind,  sex, age, ln, yr) * weight_length_pars(sp, 0) * pow(lengths(sp, ln) + lenmid, weight_length_pars(sp, 1));
        }
        weight_hat(wtind, , sex, age, yr) = expected_weight;
      } // age
    } // yr
  } // sex
}



// ------------------------------------------------------------------------- //
// 1. DATA SECTION                                                           //
// ------------------------------------------------------------------------- //
DATA_MATRIX(lengths);       // Length bins for each species [sp, nbins]
DATA_IVECTOR(growth_model); // 0: "input", 1: "vB-classic", 2: "Richards", 3: "nonparametric LAA" [sp]
// RENAME WEIGHT TO WEIGHT_OBS

// ------------------------------------------------------------------------- //
// 2. PARAMETER SECTION                                                      //
// ------------------------------------------------------------------------- //
PARAMETER_ARRAY(mu_growth_pars);  // Mean growth curve parameters [sp, sex, par]
PARAMETER_ARRAY(re_growth_pars);  // Annual random effects for growth curve parameters [sp, sex, year, par]
PARAMETER_ARRAY(growth_ln_sd);    // Log standard deviation of length-at- min and max age [sp, sex, 2]
PARAMETER_MATRIX(weight_length_pars);         // Length-weight parameters [sp, (alpha, beta)]

for(sp = 0; sp < nspp; sp++){
  for(sex = 0; sex < nsex(sp); sex ++){
    for(yr = 0; yr < nyrs; yr++){
      for(int par = 0; par < 3; par++){
        growth_parameters(sp, sex, yr, par) = exp(mu_growth_pars(sp, sex, par) + re_growth_pars(sp, sex, yr, par));
      }
    }
  }
}



// ------------------------------------------------------------------------- //
// 3. DERIVED QUANTITIES                                                     //
// ------------------------------------------------------------------------- //
arry<Type> growth_matrix(nspp * 2 + nflt, max_sex, max_age, max_nlengths, nyrs); growth_matrix.setZero(); // growth transition matrix for each fleet and each species derived quantity (biomass and ssb)
arry<Type> weight_hat(nspp * 2 + nflt, max_sex, max_age, nyrs); weight_hat.setZero(); // Estimated weight-at-age for each fleet and each species derived quantity (biomass and ssb)
arry<Type> length_hat(nspp * 2 + nflt, max_sex, max_age, nyrs); length_hat.setZero(); // Estimated length-at-age for each fleet and each species derived quantity (biomass and ssb)


// 1. POPULATION WEIGHT-AT-AGE
for(sp = 0; sp < nspp; sp++){

  // -- 1.1. Empirical weight-at-age
  if(growth_model(sp) == 0){

    for(int sex = 0; sex < nsex(sp); sex++) {
      for(int age = 0; age < nages(sp); age++) {
        for(int yr = 0; yr < nyrs; yr++) {

          // Hindcast
          if(yr < nyrs_hind){
            yr_ind = yr;
          }

          // Projection
          if(yr >= nyrs_hind){
            yr_ind = nyrs_hind - 1;
          }

          // Biomass weight
          int wtind = (nspp - 1) * 2 * sp;
          weight_hat(wtind, sex, age, yr) = weight( pop_wt_index(sp), sex, age, yr_ind);

          // SSB weight
          int wtind = (nspp - 1) * 2 * sp + 1;
          weight_hat(wtind, sex, age, yr) = weight( ssb_wt_index(sp), sex, age, yr_ind);
        }
      }
    }
  }

  // -- 1.2. Estimated growth
  if(growth_model(sp) > 0){

    // Biomass weight
    int wtind = (nspp - 1) * 2 * sp;
    get_growth_matrix(wtind, sp, 0.0, nspp, nyrs,
                      nsex, nages, lengths,
                      nlengths, minage, maxage,
                      growth_parameters, growth_ln_sd,
                      weight_length_pars,
                      growth_model,
                      length_hat,     // Pass by reference
                      growth_matrix,  // Pass by reference
                      weight_hat      // Pass by reference
    );

    // SSB weight (at mongth of spawning)
    int wtind = (nspp - 1) * 2 * sp + 1;
    get_growth_matrix(wtind, sp, spawn_month(sp)/12.0, nspp, nyrs,
                      nsex, nages, lengths,
                      nlengths, minage, maxage,
                      growth_parameters, growth_ln_sd,
                      weight_length_pars,
                      growth_model,
                      length_hat,     // Pass by reference
                      growth_matrix,  // Pass by reference
                      weight_hat      // Pass by reference
    );
  }
}

// 2. FLEET WEIGHT-AT-AGE
for(flt_ind = 0; flt_ind < n_flt; flt_ind++){

  sp = flt_spp(flt);
  mo = flt_month(flt);

  // -- 1.1. Empirical weight-at-age
  if(growth_model(sp) == 0){

    for(int sex = 0; sex < nsex(sp); sex++) {
      for(int age = 0; age < nages(sp); age++) {
        for(int yr = 0; yr < nyrs; yr++) {

          // Hindcast
          if(yr < nyrs_hind){
            yr_ind = yr;
          }

          // Projection
          if(yr >= nyrs_hind){
            yr_ind = nyrs_hind - 1;
          }

          // Assign weight
          int wtind = nspp * 2 + flt;
          weight_hat(wtind, sex, age, yr) = weight( flt_wt_index(flt), sex, age, yr_ind);
        }
      }
    }
  }

  // -- 1.2. Estimated growth
  if(growth_model(sp) > 0){

    // Estimated weight
    int wtind = nspp * 2 + flt;
    get_growth_matrix(wtind, sp, mo / 12.0, nspp, nyrs,
                      nsex, nages, lengths,
                      nlengths, minage, maxage,
                      growth_parameters, growth_ln_sd,
                      weight_length_pars,
                      growth_model,
                      length_hat,     // Pass by reference
                      growth_matrix,  // Pass by reference
                      weight_hat      // Pass by reference
    );
  }
}





// -- 11.1. CAAL
matrix<Type>  caal_hat = caal_obs; caal_hat.setZero();                            // Estimated CAAL

caal_hat.setZero();
for(int caal_ind = 0; caal_ind < caal_ctl.rows(); caal_ind++){

  flt = caal_ctl(caal_ind, 0) - 1;            // Temporary fishery index
  sp = caal_ctl(caal_ind, 1) - 1;             // Temporary index of species
  sex = caal_ctl(caal_ind, 2);                // Temporary index for caal sex (0 = combined, 1 = female, 2 = male)
  yr = caal_ctl(caal_ind, 4);                 // Temporary index for years of data
  mo = caal_n(caal_ind, 0);                   // Temporary index for month
  int wtind = nspp * 2 + flt;

  // Hindcast
  if(yr < nyrs_hind){
    yr_ind = yr;
  }

  // Projection
  if(yr >= nyrs_hind){
    yr_ind = nyrs_hind - 1;
  }

  for(age = 0; age < nages(sp); age++) {
    for(ln = 0; ln < nlengths(sp); ln++) {

      // - Fishery
      if(flt_type(flt) == 1){
        // FIXME: we can use either age or length selectivity
        pred_CAAL(flt, sex, ln, age, yr) = F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr) * growth_matrix(wtind,  sex, age, ln, yr); // * F_flt_ln(flt, sex, ln, yr)
      }

      // - Survey
      if(flt_type(flt) == 2){
        pred_CAAL(flt, sex, age, ln, yr_ind) = N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr)) * growth_matrix(wtind,  sex, age, ln, yr); //TODO length-based selectivity
        // * sel_ln(flt, sex, ln, yr_ind)
      }

      // // Marginals
      // pred_CAA(flt, sex, age, yr) += pred_CAAL(flt, sex, ln, age, yr); // Predicted catch-at-age (numbers)
      // pred_CAL(flt, sex, ln, yr) += pred_CAAL(flt, sex, ln, age, yr);  // Predicted catch-at-length
    }
  }
}


// Adjustment for joint sex composition data
joint_adjust(caal_ind) = 1;
if(flt_sex == 3){
  joint_adjust(caal_ind) = 2;
}

// Get true age comp
for(age = 0; age < nages(sp) * joint_adjust(caal_ind); age++) {
  if(n_hat(caal_ind) > 0){ //FIXME - not differentiable
    true_age_caal_hat(caal_ind, age ) = age_hat(caal_ind, age ) / n_hat(caal_ind);
  }
}


// Adjust for aging error
for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
  for(int true_age = 0; true_age < nages(sp); true_age++) {
    age_obs_hat(caal_ind, obs_age) += age_hat(caal_ind, true_age ) * age_error(sp, true_age, obs_age);
  }
}

// Adjust for aging error for joint data
if(flt_sex == 3){
  for(int obs_age = nages(sp); obs_age < nages(sp) * 2; obs_age++) {
    for(int true_age = nages(sp); true_age < nages(sp) * 2; true_age++) {

      // Adjust indexing for joint age/length comp
      int true_age_tmp = true_age - nages(sp);
      int obs_age_tmp = obs_age - nages(sp);

      age_obs_hat(caal_ind, obs_age) += age_hat(caal_ind, true_age ) * age_error(sp, true_age_tmp, obs_age_tmp);
    }
  }
}


;
// * sel_l(flt, sex, ln, yr_ind)
for(int f = 0; f < n_fleets; f++)
{
  lsum.setZero();
  asum.setZero();
  pred_catch(y,f) = 0.0;
  for(int a = 0; a < n_ages; a++){
    for(int l = 0; l < n_lengths; l++) {
      // for model years, we can use either age or len selex: F vector goes until n_model_years
      if(y < n_years_model) pred_CAAL(flt, sex, ln, age, yr) = selAA(selblock_pointer_fleets(usey,f)-1)(usey,a)*selLL(selblock_pointer_fleets(usey,f)-1)(usey,l)*catch_phi_mat(y,l,a)*NAA(y,a)*F(y,f)*(1-exp(-ZAA(y,a)))/ZAA(y,a);
      // for projection years, we can only use selectivity at age as calculated above. TODO: implement len-selex for projection years.
      if(y > n_years_model-1) pred_CAAL(flt, sex, ln, age, yr) = catch_phi_mat(y,l,a)*NAA(y,a)*FAA(y,f,a)*(1-exp(-ZAA(y,a)))/ZAA(y,a);
      lsum(l) += pred_CAAL(flt, sex, ln, age, yr);
      asum(a) += pred_CAAL(flt, sex, ln, age, yr);
    }
  }
  for(int l = 0; l < n_lengths; l++) pred_CAL(y,f,l) = lsum(l); // predicted catch-at-length
  for(int a = 0; a < n_ages; a++) pred_CAA(y,f,a) = asum(a); // predicted catch-at-age (numbers)

  // -- CAAL
  vector<Type> paa_obs_y(nages(sp));
  paa_obs_y.setZero();
  if(use_index_aging_error(i) == 1) { // use aging error
    for(int age = 0; age < nages(sp); age++){
      for(int a2 = 0; a2 < nages(sp); a2++) tmp_aging(a2, age) = pred_CAAL(yr, i, ln, age) * age_error(sp, true_age, obs_age);
    }
    tmp_agecomps = tmp_aging.rowwise().sum();
    for(int age = 0; age < nages(sp); age++) {
      pred_index_caal(yr, i, ln, age) = tmp_agecomps(age)/lsum(ln); // this object will contain the paa with aging error
      t_pred_paa(age) = pred_index_caal(yr, i, ln, age);
    }
  } else { // not use aging error
    for(int age = 0; age < nages(sp); age++){
      pred_index_caal(yr, i, ln, age) = pred_CAAL(yr, i, ln, age)/lsum(ln);
      t_pred_paa(age) = pred_index_caal(yr, i, ln, age);
    }
  }
