#include <TMB.hpp>

template <class Type>
array<Type> pred_length(matrix<Type> mLAA_jan1, int n_yrs, int n_years_model, vector<Type> expSD,
                        array<Type> growth_pars, vector<Type> lengths, vector<Type> fracyr_vec, int growth_model, Type age_L1){

  // ------------------------------------------------------------------------- //
  // 1. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //

  PARAMETER_MATRIX(ln_mean_Linf);   // Mean asymptotic length [Sp, sex]
  PARAMETER_MATRIX(ln_mean_K);      // Mean growth coefficient [Sp, sex]
  PARAMETER_MATRIX(ln_mean_L1);     // Mean length at reference age-a [Sp, sex]

  PARAMETER_ARRAY(Linf_dev);        // Time-varying deviate of asymptotic length [Sp, sex, year]
  PARAMETER_ARRAY(K_dev);           // Time-varying deviate of growth coefficient [Sp, sex, year]
  PARAMETER_ARRAY(L1_dev);          // Time-varying deviate of mean length at reference age-a [Sp, sex, year]

  PARAMETER_VECTOR(Linf_dev_ln_sd); // Log standard deviation of time varying deviate of asymptotic length [Sp, sex, year]
  PARAMETER_VECTOR(K_dev_ln_sd);    // Log standard deviation of time varying deviate of growth coefficient [Sp, sex, year]
  PARAMETER_VECTOR(L1_dev_ln_sd);   // Log standard deviation of time varying deviate of mean length at reference age-a [Sp, sex, year]


  // ------------------------------------------------------------------------- //
  // 2. CALCULATE MEAN GROWTH AT TIME-STEP 0 (i.e. Jan 1)                      //
  // ------------------------------------------------------------------------- //
  // Calculate mean-LAA, SDAA, and transition matrix, for all years:
  // lengths is vector with lengths mm (2, 4, 6, 8, etc)
  Type len_bin = lengths(1) - lengths(0);
  Type Fac1 = 0.0;
  Type Fac2 = 0.0;
  Type Ll1p = 0.0;
  Type Llp = 0.0;
  Type Slope = 0.0;
  Type b_len = 0.0;
  Type last_linear = 0.0;

  // required for mean length at age plus group:
  Type current_size = 0.0;
  Type temp = 0.0;
  Type temp1 = 0.0;
  Type temp3 = 0.0;
  Type temp4 = 0.0;
  Type div_age = 0.0;

  array<Type>   length_sd(nspp, max_sex, max_nages, nyrs); length_sd.setZero();             // SD in length-at-age
  array<Type>   length(nspp, max_sex, max_nlengths, nyrs); length.setZero();            // Length-at-age
  array<Type>   weight(nspp, max_sex, max_nages, nyrs); weight.setZero();               // Wength-at-age

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
            length(sp, sex, age, yr) = exp(length_par(sp, sex, age) + length_par_re(sp, sex, age, yr));
            break;

          default:
            error("Invalid 'growth_model");
          } // Growth_model switch

          // Correction for oldest age (as in SS)
          // - parametric growth only
          if(growth_model < 3 & age == (nages(sp)-1)) {
              current_size = length(sp, sex, age, yr);
              temp = 0;
              temp1 = 0;
              temp3 = growth_pars(sp, sex, yr, 1) - current_size;
              temp4 = 1;
              for(int age = 0; age < nages(sp); age++) {
                div_age = (age+0.0)/(nages(sp)+0.0);
                temp += temp4 * (current_size + div_age*temp3);
                temp1 += temp4;
                temp4 *= exp(-0.2);
              }
              length(sp, sex, age, yr) = temp/temp1; // oldest age
            }
          }
        } // Year
      } // Age
    } // Sex
  } // Species



      // SD calculation:
      if(growth_model < 3) { // for parametric approach
        for(int age = 0; age < nages(sp); age++) {
          if((age + 1.0) < age_L1) { // same as SD1
            length_sd(sp, sex, age, yr) = exp(length_ln_sd(0, sp));
          }
        } else {

        }
      }


      if(growth_model == 3) { // only works for yr effect

        for(int age = 0; age < nages(sp); age++) {

        }


        Slope = (exp(length_ln_sd(1, sp)) - exp(length_ln_sd(0, sp)))/(length(yr, nages(sp)-1)-length(yr, 0));
        length_sd(sp, sex, age, yr) = exp(length_ln_sd(0, sp)) + Slope*(length(sp, sex, age, yr)-length(yr, 0));
      }

      for(int ln = 0; ln < nlengths(sp); ln++) {

        if(ln == 0) {
          Fac1 = (Lmin_sp + len_bin - length(sp, sex, age, yr))/length_sd(sp, sex, age, yr); // upper limit smallest len bin, important colsums = 0
          growth_matrix(sp, sex, age, ln, yr) = pnorm(Fac1);
        } else {
          if(ln == (nlengths(sp)-1)) {
            Fac1 = (Lmax_sp - length(sp, sex, age, yr))/length_sd(sp, sex, age, yr);
            growth_matrix(sp, sex, age, ln, yr) = 1.0 - pnorm(Fac1);
          } else {
            Ll1p = lengths(ln+1);
            Llp = lengths(ln);
            Fac1 = (Ll1p - length(sp, sex, age, yr))/length_sd(sp, sex, age, yr);
            Fac2 = (Llp - length(sp, sex, age, yr))/length_sd(sp, sex, age, yr);
            growth_matrix(sp, sex, age, ln, yr) = pnorm(Fac1) - pnorm(Fac2);
          }
        }

      }// loop length
    } // End age loop
  } // End yr loop
} // End sex loop
} // End sp loop

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
