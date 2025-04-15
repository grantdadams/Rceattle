#include <TMB.hpp>

template <class Type>
array<Type> pred_length(matrix<Type> mLAA_jan1, int n_yrs, int n_years_model, vector<Type> expSD,
                        array<Type> growth_pars, vector<Type> lengths, vector<Type> fracyr_vec, int growth_model, Type age_L1){

  // ------------------------------------------------------------------------- //
  // 1. DATA SECTION                                                           //
  // ------------------------------------------------------------------------- //
  DATA_MATRIX(lengths); // Length bins for each species

  // ------------------------------------------------------------------------- //
  // 2. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //
  PARAMETER_ARRAY(mu_growth_pars);  // Mean growth curve parameters [sp, sex, par]
  PARAMETER_ARRAY(re_growth_pars);  // Annual random effects for growth curve parameters [sp, sex, year, par]
  PARAMETER_ARRAY(length_ln_sd);    // Log standard deviation of length-at- min and max age [sp, sex, 2]


  // ------------------------------------------------------------------------- //
  // 3. CALCULATE MEAN GROWTH AT TIME-STEP 0 (i.e. Jan 1)                      //
  // ------------------------------------------------------------------------- //
  // Calculate mean-length, SD, and growth matrix, for all years:
  // lengths is vector with lengths mm (2, 4, 6, 8, etc)
  Type Fac1 = 0.0;
  Type Fac2 = 0.0;
  Type Slope = 0.0;
  Type b_len = 0.0;
  Type last_linear = 0.0;

  // required for mean length at age plus group:
  Type temp = 0.0;
  Type temp1 = 0.0;
  Type temp3 = 0.0;
  Type temp4 = 0.0;
  Type div_age = 0.0;

  array<Type> length_sd(nspp, max_sex, max_nages, nyrs); length_sd.setZero();     // SD in length-at-age
  array<Type> length(nspp, max_sex, max_nages, nyrs); length.setZero();           // Length-at-age
  array<Type> weight(nspp, max_sex, max_nages, nyrs); weight.setZero();           // Wength-at-age
  array<Type> growth_matrix(nspp, max_sex, max_nlengths, max_nages, nyrs);        // Growth matrix

  // Loop through species, sex, ages, and years
  for(sp = 0; sp < nspp; sp++){

    Type Lmin_sp = lengths(sp, 0);
    Type Lmax_sp = lengths(sp, nlengths(sp)-1);
    Type age_L1 = minage(sp); //FIXME: double check
    Type age_L1_ceil = maxage(sp); //FIXME: double check

    for(sex = 0; sex < nsex(sp); sex ++){
      for(age = 0; age < nages(sp); age++){
        for(yr = 0; yr < nyrs; yr++){

          switch(growth_model){

          case 1:  // Von Bertalanffy Growth
            b_len = (L1(sp, age, yr) - Lmin_sp)/age_L1; // Slope from Lmin_sp to L1 (reference age)

            // Year 1
            if(yr == 0){
              if((age + 1.0) <= age_L1) { // Linear growth
                length(sp, sex, age, 0) = Lmin_sp + b_len*(age+1.0);
              } else { // use growth equation
                length(sp, sex, age, 0) = growth_pars(sp, sex, yr, 1) + (growth_pars(sp, sex, yr, 2) - growth_pars(sp, sex, yr, 1)) * exp(-growth_pars(sp, sex, yr, 0)*(age+1.0-age_L1));
              }
            }

            // All other years
            if(yr > 0){
              // First age
              if((age + 1.0) < age_L1) { // linear growth
                length(sp, sex, age, yr) = Lmin_sp + b_len * (age+1.0);
              } else { // use growth equation
                if((age + 1) == age_L1_ceil) { // linear growth + growth equation
                  last_linear = Lmin_sp + b_len * age_L1; // last age (cont) with linear growth
                  length(sp, sex, age, yr) = last_linear + (last_linear - growth_pars(sp, sex, yr, 1))*(exp(-growth_pars(sp, sex, yr, 0)*(a+1.0-age_L1)) - 1.0); // use growth parameters y
                } else { // only growth curve
                  length(sp, sex, age, yr) = length(sp, sex, age-1, yr-1) + (length(sp, sex, age-1, yr-1) - growth_pars(sp, sex, yr-1, 1))*(exp(-growth_pars(sp, sex, yr-1, 0)) - 1.0);// use growth parameters y-1 and a-1 because it is jan1
                }
              }
            }
            break;


          case 2: // Richard's growth model
            b_len = (L1(sp, age, yr) - Lmin_sp)/age_L1; // Slope from Lmin_sp to L1 (reference age)

            // Year 1
            if(yr == 0) { //yr = 0
              if((age + 1.0) <= age_L1) { // linear growth
                length(sp, sex, age, yr) = Lmin_sp + b_len * (age+1.0);
              } else { // use growth equation
                length(sp, sex, age, yr) = pow(pow(growth_pars(sp, sex, yr, 1), growth_pars(sp, sex, yr, 3)) + (pow(growth_pars(sp, sex, yr, 2), growth_pars(sp, sex, yr, 3)) - pow(growth_pars(sp, sex, yr, 1), growth_pars(sp, sex, yr, 3))) * exp(-growth_pars(sp, sex, yr, 0) * (a+1.0-age_L1)), 1/growth_pars(sp, sex, yr, 3));
              }
            }

            // All other years
            if(yr > 0){
              if((age + 1.0) < age_L1) { // linear growth
                length(sp, sex, age, yr) = Lmin_sp + b_len*(age+1.0);
              } else { // use growth equation
                if((age+1) == age_L1_ceil) { // linear growth + growth equation
                  last_linear = Lmin_sp + b_len * age_L1; // last age (cont) with linear growth
                  length(sp, sex, age, yr) = pow(pow(last_linear, growth_pars(sp, sex, yr, 3)) + (pow(last_linear,growth_pars(sp, sex, yr, 3)) - pow(growth_pars(sp, sex, yr, 1),growth_pars(sp, sex, yr, 3)))*(exp(-growth_pars(sp, sex, yr, 0)*(a+1.0-age_L1)) - 1.0), 1/growth_pars(sp, sex, yr, 3)); // use growth parameters y
                } else { // only growth curve
                  length(sp, sex, age, yr) = pow(pow(length(sp, sex, age-1, yr-1), growth_pars(sp, sex, yr-1,3)) + (pow(length(sp, sex, age-1, yr-1), growth_pars(sp, sex, yr-1, 3)) - pow(growth_pars(sp, sex, yr-1, 1), growth_pars(sp, sex, yr-1, 3)))*(exp(-growth_pars(sp, sex, yr-1,0)) - 1.0), 1/growth_pars(sp, sex, yr-1, 3));
                }
              }
            }
            break;

          case 3: // Non-parametric (Free parameters)
            // length(sp, sex, age, yr) = exp(length_par(sp, sex, age) + length_par_re(sp, sex, age, yr));
            break;

          default:
            error("Invalid 'growth_model");
          } // Growth_model switch

          // Correction for oldest age (as in SS)
          // - parametric growth only
          if(growth_model < 3 & age == (nages(sp)-1)) {
            temp = 0;
            temp1 = 0;
            temp3 = growth_pars(sp, sex, yr, 1) - length(sp, sex, age, yr);
            temp4 = 1;
            for(int age = 0; age < nages(sp); age++) {
              div_age = (age+0.0)/(nages(sp)+0.0);
              temp += temp4 * (length(sp, sex, nages(sp)-1, yr) + div_age * temp3); //  accumulate weight for mean size calculation
              temp1 += temp4; //  accumulate numbers to create denominator for mean size calculation
              temp4 *= exp(-0.2); //  decay numbers at age by exp(-0.xxx)
            }
            length(sp, sex, age, yr) = temp/temp1; // oldest age
          }
        } // Year
      }  // Age
    } // Sex
  } // Species


  // SD calculation:
  for(sex = 0; sex < nsex(sp); sex ++){
    for(age = 0; age < nages(sp); age++){
      for(yr = 0; yr < nyrs; yr++){

        // Parametric models
        if(growth_model < 3) {

          if((age + 1.0) < age_L1) { // same as SD1
            length_sd(sp, sex, age, yr) = exp(length_ln_sd(sp, sex, 0);
          } else {
            if(age == (nages(sp)-1)) { // same as SDA
              length_sd(sp, sex, age, yr) = exp(length_ln_sd(sp, sex, 1);
            } else { // linear interpolation
              Slope = (exp(length_ln_sd(sp, sex, 1) - exp(length_ln_sd(sp, sex, 0))/(growth_pars(sp, sex, yr, 1) - growth_pars(sp, sex, yr, 2));
                             length_sd(sp, sex, age, yr) = exp(length_ln_sd(sp, sex, 0) + Slope * (length(sp, sex, age, yr) - growth_pars(sp, sex, yr, 2));
            }
          }

          // Free parameters
          if(growth_model == 3) {
            // Slope = (exp(length_ln_sd(sp, sex, 1) - exp(length_ln_sd(sp, sex, 0))/(length(sp, sex, nages(sp)-1, yr) - length(sp, sex, 0, yr));
            // length_sd(sp, sex, age, yr) = exp(length_ln_sd(sp, sex, 0) + Slope * (length(sp, sex, age, yr) - length(sp, sex, 0, yr));
          }
        }
      }
    }

    // Derive growth matrix
    for(sex = 0; sex < nsex(sp); sex ++){
      for(age = 0; age < nages(sp); age++){
        for(yr = 0; yr < nyrs; yr++){
          for(ln = 0; ln < nlengths(sp); ln++) {

            if(ln == 0) {
              Fac1 = (Lmin_sp + lengths(sp, 1) - lengths(sp, 0) - length(sp, sex, age, yr)) / length_sd(sp, sex, age, yr); // upper limit smallest len bin, important colsums = 0
              growth_matrix(sp, sex, ln, age, yr) = pnorm(Fac1);
            } else {
              if(ln == (nlengths(sp)-1)) {
                Fac1 = (Lmax_sp - length(sp, sex, age, yr)) / length_sd(sp, sex, age, yr);
                growth_matrix(sp, sex, ln, age, yr) = 1.0 - pnorm(Fac1);
              } else {
                Fac1 = (lengths(sp, ln + 1) - length(sp, sex, age, yr)) / length_sd(sp, sex, age, yr);
                Fac2 = (lengths(sp, ln) - length(sp, sex, age, yr)) / length_sd(sp, sex, age, yr);
                growth_matrix(sp, sex, ln, age, yr) = pnorm(Fac1) - pnorm(Fac2);
              }
            }
          }
        }
      }
    }
  }


  return growth_matrix;
}



pred_CAAL(flt, sex, age, ln, yr_ind) = N_at_age(sp, sex, age, yr)  * sel(flt, sex, age, yr_ind) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr)) * growth_mat(sp, sex, age, ln, yr); //TODO length-based selectivity

;
// * sel_l(flt, sex, ln, yr_ind)


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
