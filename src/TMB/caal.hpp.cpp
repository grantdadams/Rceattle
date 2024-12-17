
template <class Type>
array<Type> pred_LAA(matrix<Type> mLAA_jan1, int n_yrs, int n_years_model, vector<Type> expSD,
                     array<Type> GW_par, vector<Type> lengths, vector<Type> fracyr_vec, int growth_model, Type age_L1){

  // --------------------------------------------------------------------------
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
  Type temp2 = exp(-0.2);
  Type temp3 = 0.0;
  Type temp4 = 0.0;
  Type div_age = 0.0;

  array<Type>   SDAA(nspp, 2, max_nages, nyrs); SDAA.setZero();             // SD in length-at-age
  array<Type>   LAA(nspp, 2, max_nlengths, nyrs); LAA.setZero();            // Length-at-age
  array<Type>   WAA(nspp, 2, max_nages, nyrs); WAA.setZero();               // Wength-at-age

  // Loop through species, sex, ages, and years
  for(sp = 0; sp < nspp; sp++){

    Type Lmin_sp = min(lengths(sp));
    Type Lmax_sp = max(lengths);
    Type age_L1 = minage(sp)ln

    for(sex = 0; sex < nsex(sp); sex ++){
      for(yr = 0; yr < nyrs; yr++){

        // Mean length-at-age
        for(age = 0; age < nages(sp); age++){

          //if(growth_model == 1) { // vB classic growth model

          b_len = (L1(sp, age, yr) - Lmin_sp)/age_L1; // Slope from Lmin_sp to L1 (reference age)

          if(yr == 0) { //yr = 0
            if((age + 1.0) <= age_L1) { // linear growth
              LAA(sp, sex, age, yr) = Lmin_sp + b_len*(age+1.0);
            } else { // use growth equation
              LAA(sp, sex, age, yr) = Linf(sp, age, yr) + (L1(sp, age, yr) - Linf(sp, age, yr)) * exp(-K(sp, age, yr)*(age+1.0-age_L1));
            }


            if(yr > 0) { // yr > 0
              if((age + 1.0) < age_L1) { // linear growth
                LAA(sp, sex, age, yr) = Lmin_sp + b_len*(age+1.0);
              } else { // use growth equation

                if((age+1) == age_L1_ceil) { // linear growth + growth equation
                  last_linear = Lmin_sp + b_len * age_L1; // last age (cont) with linear growth
                  LAA(sp, sex, age, yr) = last_linear + (last_linear - Linf(sp, age, yr))*(exp(-K(sp, age, yr)*(age+1.0-age_L1)) - 1.0); // use growth parameters yr

                } else { // only growth curve
                  LAA(sp, sex, age, yr) = LAA(sp, sex, age-1, yr-1) + (LAA(sp, sex, age-1, yr-1) - Linf(sp, yr-1, age-1))*(exp(-K(sp, yr-1, age-1)) - 1.0);// use growth parameters yr-1 and age-1 because it is jan1
                }
              }
            }

            // correction for oldest age (as in SS)
            if(age == (nages(sp)-1)) {
              current_size = LAA(sp, sex, age, yr);
              temp = 0;
              temp1 = 0;
              temp3 = Linf(sp, age, yr) - current_size;
              temp4 = 1;

              for(int age = 0; age <= nages(sp); age++) {
                div_age = (age+0.0)/(nages(sp)+0.0);
                temp += temp4*(current_size + div_age*temp3);
                temp1 += temp4;
                temp4 *= temp2;
              }
              LAA(sp, sex, age, yr) = temp/temp1; // oldest age
            }
          }
          /*
           if(growth_model == 2) { // Richards growth model
           b_len = (L1(sp, age, yr) - Lmin_sp)/age_L1; // Slope from Lmin to L1
           if(yr == 0) { //yr = 0
           if((age + 1.0) <= age_L1) { // linear growth
           LAA(sp, sex, age, yr) = Lmin_sp + b_len*(age+1.0);
           } else { // use growth equation
           LAA(sp, sex, age, yr) = pow(pow(Linf(sp, age, yr),GW_par(yr, age, 3)) + (pow(L1(sp, age, yr),GW_par(yr, age, 3)) - pow(Linf(sp, age, yr),GW_par(yr, age, 3))) * exp(-K(sp, age, yr)*(age+1.0-age_L1)),1/GW_par(yr, age, 3));
           }
           } else { // yr > 0
           if((age + 1.0) < age_L1) { // linear growth
           LAA(sp, sex, age, yr) = Lmin_sp + b_len*(age+1.0);
           } else { // use growth equation
           if((age+1) == age_L1_ceil) { // linear growth + growth equation
           last_linear = Lmin_sp + b_len*age_L1; // last age (cont) with linear growth
           LAA(sp, sex, age, yr) = pow(pow(last_linear, GW_par(yr, age, 3)) + (pow(last_linear, GW_par(yr, age, 3)) - pow(Linf(sp, age, yr),GW_par(yr, age, 3)))*(exp(-K(sp, age, yr)*(age+1.0-age_L1)) - 1.0),1/GW_par(yr, age, 3)); // use growth parameters yr
           } else { // only growth curve
           LAA(sp, sex, age, yr) = pow(pow(LAA(sp, sex, age-1, yr-1),GW_par(yr-1, age-1, 3)) + (pow(LAA(sp, sex, age-1, yr-1),GW_par(yr-1, age-1, 3)) - pow(Linf(sp, yr-1, age-1),GW_par(yr-1, age-1, 3)))*(exp(-K(sp, yr-1, age-1)) - 1.0),1/GW_par(yr-1, age-1, 3));
           }
           }
           }
           // correction for oldest age (as in SS)
           if(age == (nages(sp)-1)) {
           current_size = LAA(sp, sex, age, yr);
           temp = 0;
           temp1 = 0;
           temp3 = Linf(sp, age, yr) - current_size;
           temp4 = 1;
           for(int age = 0; age <= nages(sp); age++) {
           div_age = (age+0.0)/(nages(sp)+0.0);
           temp += temp4*(current_size + div_age*temp3);
           temp1 += temp4;
           temp4 *= temp2;
           }
           LAA(sp, sex, age, yr) = temp/temp1; // oldest age
           }
           }

           if(growth_model == 3) LAA(sp, sex, age, yr) = LAA_par(yr, age); // LAA model
           */
        } // End age loop


        // SD calculation:
        for(int age = 0; age < nages(sp); age++) {
          // if(growth_model < 3) { // for parametric approach
          if((age + 1.0) < age_L1) { // same as SD1
            SDAA(sp, sex, age, yr) = exp(laa_ln_sd(0, sp));
          } else {
            if(age == (nages(sp)-1)) { // same as SDA
              SDAA(sp, sex, age, yr) = exp(laa_ln_sd(1, sp));
            } else { // linear interpolation
              Slope = (exp(laa_ln_sd(1, sp)) - exp(laa_ln_sd(0, sp)))/(Linf(sp, age, yr)-L1(sp, age, yr));
              SDAA(sp, sex, age, yr) = exp(laa_ln_sd(0, sp)) + Slope*(LAA(sp, sex, age, yr)-L1(sp, age, yr));
            }
          }
          //}
          // if(growth_model == 3) { // only works for yr effect
          //   Slope = (exp(laa_ln_sd(1, sp)) - exp(laa_ln_sd(0, sp)))/(LAA(yr, nages(sp)-1)-LAA(yr, 0));
          //   SDAA(sp, sex, age, yr) = exp(laa_ln_sd(0, sp)) + Slope*(LAA(sp, sex, age, yr)-LAA(yr, 0));
          // }

          for(int ln = 0; ln < nlengths(sp); ln++) {

            if(ln == 0) {
              Fac1 = (Lmin_sp + len_bin - LAA(sp, sex, age, yr))/SDAA(sp, sex, age, yr); // upper limit smallest len bin, important colsums = 0
              growth_matrix(sp, sex, age, ln, yr) = pnorm(Fac1);
            } else {
              if(ln == (nlengths(sp)-1)) {
                Fac1 = (Lmax_sp - LAA(sp, sex, age, yr))/SDAA(sp, sex, age, yr);
                growth_matrix(sp, sex, age, ln, yr) = 1.0 - pnorm(Fac1);
              } else {
                Ll1p = lengths(ln+1);
                Llp = lengths(ln);
                Fac1 = (Ll1p - LAA(sp, sex, age, yr))/SDAA(sp, sex, age, yr);
                Fac2 = (Llp - LAA(sp, sex, age, yr))/SDAA(sp, sex, age, yr);
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

comp_hat(comp_ind, ln )
  pred_CAAL(flt, sex, age, ln, yr_ind) = index_q(flt, yr_ind)  * sel(flt, sex, age, yr_ind) * N_at_age(sp, sex, age, yr) * exp(-Type(mo/12.0) * Z_at_age(sp, sex, age, yr)) * growth_mat(sp, sex, age, ln, yr); //TODO length-based selectivity
// * sel_l(flt, sex, ln, yr_ind)


// -- CAAL
vector<Type> paa_obs_y(nages(sp));
paa_obs_y.setZero();
if(use_index_aging_error(i) == 1) { // use aging error
  for(int age = 0; age < nages(sp); age++){
    for(int a2 = 0; a2 < nages(sp); a2++) tmp_aging(a2, age) = pred_CAAL(yr, i, ln, age)*index_aging_error(i, a2, age);
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
if(yr < n_years_model) if(use_index_caal(yr, i, ln) == 1) {
  for(int age = 0; age < nages(sp); age++) paa_obs_y(age) = index_caal(i, yr, ln, age);
  //NB: indexing in obsvec MUST be: keep_Ipaa(i, yr, 0),...,keep_Ipaa(i, yr, 0) + keep_Ipaa(i, yr, 1) - 1
  //keep_Ipaa(i, yr, 0) is first val, keep_Ipaa(i, yr, 1) is the length of the vector
  vector<Type> tf_paa_obs = obsvec.segment(keep_Icaal(i, yr, ln, 0), keep_Icaal(i, yr, ln, 1));
  vector<int> ages_obs_y = agesvec.segment(keep_Icaal(i, yr, ln, 0), keep_Icaal(i, yr, ln, 1));
  nll_index_caal(yr, i, ln) -= get_acomp_ll(tf_paa_obs, t_pred_paa,
