#if defined(__clang__)
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <TMB.hpp>
#if defined(__clang__)
# pragma clang diagnostic pop
#endif

// GCC/Eigen emits false -Warray-bounds positives where its vectorized reductions
// inline TMB's multi-dimensional indexing (e.g. the 5-D suitability() access in
// predation.hpp). A source pragma rather than a -Wno-array-bounds compiler flag,
// so CRAN's compilation-flags check stays clean and other diagnostics still show.
#if defined(__GNUC__)
# pragma GCC diagnostic ignored "-Warray-bounds"
#endif

#include "helper_functions.hpp"
#include "comp_osa.hpp"
#include "comp_sim.hpp"
#include "growth.hpp"
#include "selectivity.hpp"
#include "recruitment.hpp"
#include "spr.hpp"
#include "bioenergetics.hpp"
#include "predation.hpp"
#include "diet_data.hpp"
#include "linkage.hpp"
#include "dsem.hpp"

/**
 * @file ceattle.cpp
 * @brief The CEATTLE objective function: configuration, population dynamics,
 *   predation, likelihood and report. The numbered section index is in the
 *   banner below; the process equations live in the headers this file includes.
 */

// Reads an R list() of matrices into a vector<matrix<Type>>: the per-fleet
// survey-index covariance matrices Sigma supplied when Index_loglike is
// "MVN"/"MVNORM" (the AMAK/ebswp DoCovBTS covariance survey likelihood; sections
// 8.2b and 13.1). Other fleets pass a 1x1 inert dummy, so the list is always
// length n_flt and indexes by fleet code.
template<class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x) {  // x = R list of matrices
    (*this).resize(LENGTH(x));
    for(int i = 0; i < LENGTH(x); i++){
      SEXP m = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(m);
    }
  }
};

/* ------------------------------------------------------------------------- //
 * CEATTLE: Bioenergetic-based Assessment for Understanding Biomass Linkages To
 * The Environment -- a multispecies statistical model in Template Model
 * Builder. Version: see DESCRIPTION.
 *
 * Follows the Google C++ style guide
 * (https://google.github.io/styleguide/cppguide.html#General_Naming_Rules),
 * except that free functions stay snake_case.
 *
 * CITATIONS:
 * 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
 * 2. Adams, G. D., Holsman, K. K., Barbeaux, S. J., Dorn, M. W., Ianelli, J. N., Spies, I., ... & Punt, A. E. (2022). An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fisheries Research, 251, 106303.
 * 3. Wassermann, S. N., Adams, G. D., Haltuch, M. A., Kaplan, I. C., Marshall, K. N., & Punt, A. E. (2025). Even low levels of cannibalism can bias population estimates for Pacific hake. ICES Journal of Marine Science, 82(1), fsae064.
 * ------------------------------------------------------------------------- //
 *  INDEX:
 *  0. Load dependencies
 *  1. Model configuration
 *  2. Model inputs
 *  3. Model parameters
 *  4. Derived quantities
 *  5. Initial calculations
 *  6. Population dynamics
 *  7. Predation mortality equations
 *  8. Index components equations
 *  9. Fishery components equations
 *  10. Composition equations
 *  11. Predicted stomach content
 *  12. Derived quantities
 *  13. Likelihood components
 *  14. Report section
 *  15. Model return/end
 *
 */


template<class Type>
Type objective_function<Type>::operator() () {
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  /** ------------------------------------------------------------------------ //
   * 1. MODEL CONFIGURATION                                                    //
   * ------------------------------------------------------------------------- */
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(estimateMode);             // 0 = hindcast + HCR projection, 1 = hindcast, 2 = projection from inits, 3 = build only, 4 = smoke test
  DATA_INTEGER(msmMode);
  //    0 = single species (no predation mortality)
  //    1 = Type II MSVPA based (sensu Holsman et al 2015)
  //    2 = Type III MSVPA based
  //    3-9 = not yet implemented (Holling forms, predator interference and
  //          preemption, Hassell-Varley, Ecosim)
  DATA_IVECTOR(suitMode);                 // Estimate suitability
  // DATA_INTEGER(est_diet);              // Include diet data in the likelihood
  // DATA_INTEGER(avgnMode);              // Not read; predation uses avgN_at_age
  DATA_INTEGER(initMode);                 // How the age-structure is initialized
  DATA_IVECTOR(forecast);                 // Run the model in forecast mode


  // 1.2. Temporal dimensions
  DATA_INTEGER( styr );                   // Start year
  DATA_INTEGER( endyr );                  // End of estimation years
  DATA_INTEGER( projyr );                 // End year of projection

  DATA_INTEGER( srr_mse_switchyr );             // Last year of the average-recruitment window; for MSE runs.
  DATA_IVECTOR( suit_styr );              // First year of the suitability average, per predator species.
  DATA_IVECTOR( suit_endyr );             // Last year of the suitability average, per predator species.
  DATA_INTEGER( srr_hat_styr );           // First year of the stock-recruitment-penalty / env-rec window.
  DATA_INTEGER( srr_hat_endyr );          // Last year of the stock-recruitment-penalty / env-rec window.

  int nyrs = projyr - styr + 1;
  int nyrs_hind = endyr - styr + 1;

  // suit_styr / suit_endyr are per-predator (offset to model-year indices below,
  // once `nspp` is available).
  int nyrs_srrmean = srr_mse_switchyr - styr + 1;

  srr_hat_styr = srr_hat_styr - styr;
  srr_hat_endyr = srr_hat_endyr - styr;
  if(nyrs_srrmean > nyrs_hind){nyrs_srrmean = nyrs_hind;}
  if(srr_hat_styr == 0){srr_hat_styr = 1;} // R_hat starts at year 1.

  // 1.3. Number of species
  DATA_INTEGER( nspp );                   // Number of species (prey)

  // Offset per-predator suitability reference years to model-year indices and
  // count the averaging window length for each predator species.
  vector<int> nyrs_suit(nspp);
  for(int sp = 0; sp < nspp; sp++){
    suit_endyr(sp) = suit_endyr(sp) - styr;
    suit_styr(sp)  = suit_styr(sp) - styr;
    nyrs_suit(sp)  = suit_endyr(sp) - suit_styr(sp) + 1;
  }

  DATA_IVECTOR( pop_wt_index );           // Dim 1 of weight_obs to use for population dynamics
  DATA_IVECTOR( ssb_wt_index );           // Dim 1 of weight_obs to use for spawning stock biomass
  DATA_IVECTOR( pop_age_transition_index );// Dim 1 of age_trans_matrix to use for population dynamics
  DATA_IVECTOR( estDynamics );            // 0 = population parameters estimated, 1 = numbers-at-age provided, 2 = provided numbers-at-age times an estimated scalar
  pop_wt_index -= 1;                      // Indexing starts at 0
  ssb_wt_index -= 1;                      // Indexing starts at 0
  pop_age_transition_index -= 1;          // Indexing starts at 0


  // 1.4. RECRUITMENT SETTINGS
  DATA_INTEGER(srr_fun);                  // Stock recruit relationship for hindcast estimation
  DATA_INTEGER(srr_pred_fun);             // Stock recruit relationship for projection/brps/penalty
  DATA_INTEGER(proj_mean_rec);
  //    0 = project recruitment using log_R0 and rec devs
  //    1 = project recruitment using mean rec (can also have adjusted rec devs)
  DATA_INTEGER(srr_est_mode);             // Logical of wether to add normal prior to stock recruit-relationship
  DATA_VECTOR(srr_prior);                 // Prior mean for stock recruit relationship parameter
  DATA_VECTOR(srr_prior_sd);              // Prior sd for stock recruit relationship parameter
  DATA_INTEGER( niter );                  // Number of loops for MSM mode
  DATA_VECTOR( Bmsy_lim );                // Upper limit for Bmsy in ricker function. Will add penalty if 1/beta > lim


  // 1.5. HARVEST CONTROL RULE (HCR) SETTINGS
  DATA_INTEGER(HCR);                      // Function to be used for harvest control rule
  DATA_INTEGER(DynamicHCR);               // TRUE/FALSE. Whether to use static or dynamic reference points (default = FALSE)
  DATA_VECTOR(Ptarget);                   // Target SSB as a PROPORTION of static or dynamic SSB at F = 0 (B40% is 0.40)
  DATA_VECTOR(Plimit);                    // Limit SSB as a PROPORTION of static or dynamic SSB at F = 0 (B10% is 0.10)
  DATA_VECTOR(Ftarget_percent);           // Proportion of SSB-per-recruit or SB0 at F = 0 that sets the target F (F40% is 0.40)
  DATA_VECTOR(Flimit_percent);            // Proportion of SSB-per-recruit or SB0 at F = 0 that sets the limit F (F35% is 0.35)
  DATA_VECTOR(Alpha);                     // Parameter used in NPFMC Tier 3 HCR
  DATA_VECTOR(Fmult);                     // Multiplier for target fishing mortality (Fmult * Ftarget). For example 0.75 * F40%.
  DATA_VECTOR(QnormHCR);                  // Add to Flimit to set Ftarget based on the Pstar approach: Flimit + qnorm(Pstar, 0, Sigma)


  // 1.6. MODEL OBJECTS
  // 1.6.1. LOOPING INDICES -- k = observation, sp = species, sex = sex (0 = combined; 1 = females; 2 = males)
  // age = age, ln = length, yr = year
  int sp, sex, age, ln, yr;
  int index, flt;                                                         // Survey and fishery indices
  int flt_yr, flt_sex, comp_type;
  int fsh_ind, index_ind, comp_ind, yr_ind;                       // Indices for survey sets
  int wt_idx_pop, wt_idx_ssb, wt_idx_flt; // Indices for weight indices
  Type mo = 0;                                                           // Month float
  if(msmMode == 0) { niter = 1; }                                        // Number of iterations for SS mode


  /** ------------------------------------------------------------------------ //
   * 2. MODEL INPUTS                                                           //
   * ------------------------------------------------------------------------- */

  // -- 2.1. Species attributes
  DATA_IVECTOR( nsex );                   // Sexes modelled per species: 1 = combined, 2 = both.
  int max_sex = imax(nsex);               // Max sexes, for array dimensions
  DATA_VECTOR( spawn_month );             // Month of spawning to adjust mortality
  DATA_IVECTOR( nages );                  // Age bins per species.
  DATA_IVECTOR( minage );                 // Minimum age of each species
  DATA_IVECTOR( nlengths );               // Observation length bins per species.
  DATA_MATRIX(lengths);                   // Length bins for each species [sp, nlengths]
  DATA_IVECTOR( nlengths_pop );           // Population length bins per species. Equals nlengths when no population grid is given
  DATA_MATRIX( lengths_pop );             // Population length bins, lower edges (cm) [sp, nlengths_pop]: the grid the age-length key, weight and maturity are integrated on
  DATA_IMATRIX( pop_to_data_bin );        // Data length bin (0-based) that each population bin sums into [sp, nlengths_pop]
  DATA_ARRAY( NByageFixed );              // Input numbers- or index-at-age, multiplied (or not) by pop_scalar to get N_at_age
  DATA_VECTOR( MSSB0 );                   // SB0 from projecting forward in multi-species mode under no fishing
  DATA_VECTOR( MSB0 );                    // B0 from projecting forward in multi-species mode under no fishing

  int max_nlengths = imax(nlengths);      // Max nlengths, for array dimensions
  // The age-length key and selectivity-at-length live on the POPULATION length bins,
  // which are the data bins whenever no pop_lengths is supplied.
  int max_nlengths_pop = imax(nlengths_pop);

  // -- 2.3c. First and last population bin in each data bin, per species. SS3 forms
  // sel(L) * P(L|age) at population resolution and bins the result, so pred_CAAL and the
  // length compositions sum each data bin's run rather than applying one bin-average
  // selectivity across it. The map is monotone, so the runs are contiguous; with one grid
  // every run is one bin.
  matrix<int> pop_bin_lo(nspp, max_nlengths); pop_bin_lo.setZero();
  matrix<int> pop_bin_hi(nspp, max_nlengths); pop_bin_hi.setZero();
  for(int sp_i = 0; sp_i < nspp; sp_i++){
    for(int d = 0; d < nlengths(sp_i); d++){ pop_bin_lo(sp_i, d) = -1; pop_bin_hi(sp_i, d) = -2; }
    for(int lp = 0; lp < nlengths_pop(sp_i); lp++){
      int d = pop_to_data_bin(sp_i, lp);
      if(pop_bin_lo(sp_i, d) < 0) pop_bin_lo(sp_i, d) = lp;
      pop_bin_hi(sp_i, d) = lp;
    }
  }
  int max_age = imax(nages);              // Max nages, for array dimensions
  // Sized on the POPULATION bins: the non-parametric curves are built over nlengths_pop,
  // so a finer grid would otherwise write past non_par_sel.
  int max_bin = (max_age > max_nlengths_pop) ? max_age : max_nlengths_pop;

  // -- 2.2. M1_at_age specifications
  DATA_IVECTOR(M1_model);
  DATA_IVECTOR(M1_re);
  DATA_IVECTOR(M1_use_prior);
  DATA_IVECTOR(M2_use_prior);
  DATA_VECTOR(M_prior);
  DATA_VECTOR(M_prior_sd);

  // -- 2.3. Growth model specifications
  DATA_IVECTOR(growth_model); // 0: "input", 1: "vB-classic", 2: "Richards", 3: "nonparametric LAA" [sp]
  DATA_IVECTOR(growth_sd_style); // Plus-group SD-at-age treatment [sp]: 1 = WHAM (pin to exp(sd_Linf)), 2 = interpolate by length
  DATA_IVECTOR(growth_sd_form);  // Growth variability endpoints [sp]: 1 = SDs in cm (SS3 CV_Growth_Pattern 2), 2 = CVs, SD = CV * L (SS3 pattern 0)
  DATA_IVECTOR(growth_plus_length); // Plus-group mean length [sp]: 1 = M1-weighted, 2 = none (SS3 Linf_decay -998), 3 = SS3.24 (-999), 4 = SS3 decay rate
  DATA_VECTOR(plus_group_decay); // Decay rate (per year) for growth_plus_length == 4 [sp]
  DATA_IMATRIX(sel_dn6_ends);    // DoubleNormalSS3, per fleet [flt, 2]: 1 = scale the initial / final end by P5 / P6, 0 = SS3's -999 (unscaled)
  DATA_IVECTOR(mat_len_use);     // 1 = maturity-at-length replaces the age-based maturity in spawning output [sp]
  DATA_MATRIX(mat_len_pars);     // Maturity-at-length [sp, 2]: L50 (cm), logistic slope (per cm)
  DATA_VECTOR(growth_age_L1); // VB anchor age (= SS3 Growth_Age_for_L1) per sp; defaults to max(0.5, minage[sp]) in R-side fit_mod()

  // -- 2.3b. Long-format linkage table (see R/0-linkage_encode.R). Each row is one
  //          estimated coefficient connecting a process parameter (growth/M/recruit/...)
  //          to a column of the shared design matrix `linkage_X`. Empty when no
  //          build_*() supplied a `linkages` argument.
  DATA_IVECTOR(linkage_process);       // process code (recruitment = 0, M = 1, growth = 2, q = 3, sel = 4)
  DATA_IVECTOR(linkage_param);         // per-process parameter code
  DATA_IVECTOR(linkage_species);       // 1-based sp id; 0 = all
  DATA_IVECTOR(linkage_sex);           // 1-based sex id; 0 = all
  DATA_IVECTOR(linkage_age_bin);       // 1-based age id; 0 = all
  DATA_IVECTOR(linkage_fleet);         // 1-based Fleet_code; 0 = all
  DATA_IVECTOR(linkage_X_col);         // 0-based column of linkage_X
  DATA_IVECTOR(linkage_link);          // identity=0, log=1, logit=2
  // Only an identity-link recruitment offset can make the curve non-positive, so
  // the recruitment floors (sections 6.3, 6.5, 6.6 and 6.9) run only then; every
  // other model's AD tape is left exactly as it was.
  int rec_floor_on = 0;
  for (int i = 0; i < linkage_link.size(); ++i) {
    if (linkage_process(i) == RCEATTLE_PROC_RECRUIT && linkage_link(i) == 0) rec_floor_on = 1;
  }
  DATA_IVECTOR(linkage_re_index);      // -1 = fixed row; else 0-based slot in beta_linkage_re
  DATA_IVECTOR(linkage_re_sigma);      // per RE slot: its 0-based log_sigma_linkage group
  DATA_IVECTOR(linkage_re_integrate);  // per RE slot: 1 = Laplace-integrated, 0 = penalized
  DATA_IVECTOR(linkage_re_slot);       // per RE slot: 0-based position within its own parameter vector
  DATA_IVECTOR(linkage_re_struct);     // per RE group: covariance structure (0=us/IID, 1=rw, 2=ar1)
  DATA_IVECTOR(linkage_re_rho);        // per RE group: 0-based slot in trans_rho_linkage (ar1 only; -1 otherwise)
  DATA_IVECTOR(linkage_re_sigma_prior_family); // per RE group: prior on the SD (0=none,1=normal,2=lognormal,3=gamma,4=beta)
  DATA_IVECTOR(linkage_re_rho_prior_family);   // per RE group: prior on rho (natural (-1,1) scale)
  DATA_VECTOR(linkage_re_rho_prior_p1);        // per RE group: rho prior param 1
  DATA_VECTOR(linkage_re_rho_prior_p2);        // per RE group: rho prior param 2
  DATA_IVECTOR(linkage_re_obs);        // per RE group: 0-based slot in beta_linkage_obs (Rogers QAR1 observed ar1; -1 otherwise)
  DATA_VECTOR(linkage_re_obs_sd);      // per RE group: fixed measurement SD (0 if unobserved)
  DATA_VECTOR(linkage_re_obs_value);   // per beta_linkage_re slot: observed covariate value (0 if unobserved)
  DATA_IVECTOR(linkage_re_obs_mask);   // per beta_linkage_re slot: 1 if the covariate is observed that year, 0 otherwise
  DATA_IVECTOR(linkage_re_obsvec_idx); // per beta_linkage_re slot: 0-based obsvec position for OSA (-1 if unobserved / not built)
  DATA_VECTOR(linkage_re_sigma_prior_p1);      // per RE group: prior param 1
  DATA_VECTOR(linkage_re_sigma_prior_p2);      // per RE group: prior param 2
  DATA_IVECTOR(linkage_is_intercept);  // 1 if design_col == "(Intercept)", 0 otherwise
  DATA_IVECTOR(linkage_prior_family);  // none=0, normal=1, lognormal=2, gamma=3, beta=4
  DATA_VECTOR(linkage_prior_p1);       // family-specific prior param 1
  DATA_VECTOR(linkage_prior_p2);       // family-specific prior param 2
  DATA_MATRIX(linkage_X);              // dense design matrix [nyrs, n_design_cols]

  // -- 2.3.8. DSEM on the recruitment deviations
  // A dynamic structural equation model makes the recruitment deviations the
  // latent states of a GMRF, so environmental covariates and lagged or
  // simultaneous paths drive recruitment. Everything here is inert when
  // dsem_on = 0: build_dsem_objects() supplies zero-length stand-ins, and the
  // single block in section 5 that reads them is guarded on the switch, so a
  // model without a DSEM follows exactly the same path it did before.
  // Assembled R-side from dsem::dsem(run_model = FALSE); see R/0-build_DSEM.R.
  DATA_INTEGER(dsem_on);               // 1 if a DSEM is attached at all
  // 1 if that DSEM owns the recruitment deviations. Separate from dsem_on
  // because a SEM can drive one process without driving another: the GMRF
  // density and its conditioning are gated on dsem_on, while whether rec_dev
  // is a derived state or an estimated deviation is gated on this. Supplied as
  // 0 when there is no DSEM, so the recruitment density below is reached on the
  // same branch a non-DSEM model has always taken.
  DATA_INTEGER(dsem_drives_rec);
  DATA_INTEGER(nyrs_dsem);             // DSEM time steps; may be < nyrs (hindcast only)
  DATA_IVECTOR(dsem_options);          // dsem parameterization switches; see dsem.hpp
  DATA_IMATRIX(dsem_RAM);              // path structure [n_path, 6]
  DATA_VECTOR(dsem_RAMstart);          // fixed value per RAM row (used where parameter == 0)
  DATA_IVECTOR(dsem_familycode_j);     // per variable: measurement distribution
  DATA_IVECTOR(dsem_linkcode_j);       // per variable: link function
  DATA_IVECTOR(dsem_sigmastart_j);     // per variable: index into lnsigma_z
  DATA_ARRAY(dsem_eps_tj);             // per variable: known measurement SD where used
  DATA_ARRAY(dsem_y_tj);               // observations; latent columns are NA
  DATA_IVECTOR(dsem_obs_idx);          // stacked-space partition, projecting parameterizations
  DATA_IVECTOR(dsem_unobs_idx);
  DATA_IVECTOR(rec_dev_col);           // 0-based x_tj column carrying each species' recdevs, -1 where the SEM does not drive recruitment
  // Growth latents the SEM drives, one entry per named column: its 0-based x_tj
  // column, the 0-based species and sex it drives, and the 0-based slot in
  // growth_linkage_offset's last dimension (K 0, L1 1, Linf 2, Richards m 3).
  // All four are zero-length unless a sem names a growth latent.
  DATA_IVECTOR(growth_dev_col);
  DATA_IVECTOR(growth_dev_sp);
  DATA_IVECTOR(growth_dev_sex);
  DATA_IVECTOR(growth_dev_par);
  // Natural-mortality latents the SEM drives, one entry per named column: its
  // 0-based x_tj column and the 0-based species and sex it drives. There is no
  // parameter slot -- M_linkage_offset is dimensioned by AGE instead, and one
  // latent applies to every age of that species and sex. Zero-length unless a
  // sem names an M latent.
  DATA_IVECTOR(mdev_col);
  DATA_IVECTOR(mdev_sp);
  DATA_IVECTOR(mdev_sex);
  DATA_IVECTOR(rec_sd_idx);            // 1-based beta_z index of each species' recruitment SD (0 if fixed)
  DATA_VECTOR(rec_sd_fixed);           // fixed recruitment SD, used where rec_sd_idx == 0
  DATA_VECTOR(rec_sd_prior);           // prior centre for the recruitment SD
  DATA_VECTOR(rec_sd_prior_sd);        // prior log-scale SD (0 where no prior)
  DATA_IVECTOR(rec_sd_use_prior);      // 1 if the recruitment-SD prior applies to this species

  // -- 2.4. Fleet controls (i.e. how to assign data to objects)
  DATA_IVECTOR(flt_type);                 // Fleet_type: 0 = Off, 1 = Fishery, 2 = Survey (NOT a boolean)
  DATA_VECTOR(flt_month);
  DATA_IVECTOR(flt_sel_type);             // Fleet selectivity parameterization
  DATA_IVECTOR(flt_sel_dim);              // Fleet selectivity dimension (0 = age, 1 = length)
  DATA_IVECTOR(flt_n_sel_bins);           // Age/length bins for non-parametric selectivity
  DATA_IVECTOR(flt_sel_cap_bin);          // NonParametricPM (type 9): first bin (0-based, age or length per flt_sel_dim) at/after which realized selectivity is held flat (RTMB cap_old_age). < 0 = no cap.
  DATA_IVECTOR(flt_varying_sel);          // Whether time-varying selectivity is estimated
  DATA_IVECTOR(flt_spp);                  // Fleet species
  DATA_IVECTOR(bin_first_selected);       // First selected age/length BIN, 0-based (selectivity below it = 0)
  DATA_IVECTOR(sel_norm_bin1);            // Normalization BIN, 0-based, or a control flag when negative (-999 = do not normalize, -99 = normalize by the maximum); rearrange_data() converts, so NA never arrives. LogisticPM (type 11): lower bin of the selectivity-penalty range.
  DATA_IVECTOR(sel_norm_bin2);            // Upper normalization BIN, same encoding as sel_norm_bin1. LogisticPM (type 11): upper bin of the selectivity-penalty range.
  DATA_IVECTOR(sel_norm_scope);           // Whether normalization pools its reference across sexes: 0 = WithinSex (each sex divided by its own reference; both reach 1), 1 = AcrossSexes (one pooled reference; the lower sex stays below 1). Orthogonal to sel_norm_bin1/2, which say WHERE the reference is taken.
  DATA_IVECTOR(flt_sel_start_yr);         // Per-fleet selectivity start year (0-based from styr); selectivity penalties start the year after this. Default 0 (= styr).
  DATA_IVECTOR(flt_sel_pen_first_bin);    // Per-fleet first bin (0-based, age or length per flt_sel_dim) for the non-parametric shape/monotonicity penalty. < 0 -> defaults to bin_first_selected. Lets the shape constraint span a narrower range than the (possibly non-zero) first selected bin (e.g. ATS mina_ats > first selected age).
  DATA_IVECTOR(flt_sel_lead);             // 1 if this fleet's selectivity penalty should be accumulated; 0 if it mirrors an earlier fleet's selectivity (same Selectivity_index + type) so the shared penalty is counted once.
  DATA_IVECTOR(flt_q_lead);               // As flt_sel_lead, for catchability: 1 if this fleet carries the q prior / deviate penalties, 0 if it shares an earlier fleet's Q_index so the shared block is counted once.
  DATA_IVECTOR(flt_sel_pen_last_bin);     // Per-fleet last bin (0-based, = left bin of the last adjacent pair) for the non-parametric shape penalty. < 0 -> defaults to nbins-2 (whole range).
  DATA_IVECTOR(flt_sel_shape_mode);       // Non-parametric shape-penalty mode: 0 = directional (sign of Sel_curve_pen1 -> penalize decreasing/increasing, one-sided, ADMB/AMAK); 1 = smooth (two-sided d^2 over adjacent ages, RTMB "rpm").
  DATA_VECTOR(flt_sel_avgsel_pen);        // Per-fleet weight on the AMAK "avgsel" base-level penalty (type 9 only); the form is penalty (5) in section 13. 0 = off (default).
  DATA_IVECTOR(comp_ll_type);             // Composition log-likelihood type
  DATA_IVECTOR(comp_accum_young);         // Per-fleet age/length-composition young-tail accumulation bin (1-based ordinal on the comp dimension); bins below it fold into it. 1 (or <1) = no young accumulation (AFSC ac_yng).
  DATA_IVECTOR(comp_accum_old);           // Per-fleet age/length-composition old-tail accumulation bin (1-based); bins above it fold into it. 0/NA or >= nbins = no old accumulation (AFSC ac_old).
  DATA_IVECTOR(caal_ll_type);             // CAAL composition log-likelihood type
  DATA_IVECTOR(flt_units);                // Fleet units (1 = weight, 2 = numbers)
  DATA_IVECTOR(flt_wt_index);             // Vector to save 1st dim of weight to use for weight-at-age
  DATA_IVECTOR(flt_age_transition_index); // Dim 1 of age_trans_matrix each fleet reads for its ALK
  DATA_IVECTOR(flt_ageing_error_index);   // Ageing error matrix each fleet reads (0-based); the fleet's own species unless fleet_control gives one
  DATA_IVECTOR(est_index_q);              // Catchability form per fleet; see q_map and sections 8.2/8.2b
  DATA_IVECTOR(index_varying_q);          // Whether time-varying q is estimated
  DATA_IVECTOR(est_sigma_index);          // Whether survey sigma is estimated
  DATA_IVECTOR(est_sigma_fsh);            // Whether fishery sigma is estimated

  // -- 2.4.1 Fishery Components
  DATA_IMATRIX( catch_ctl );              // Fishery catch metadata; columns = Fleet_code, Species, Year
  DATA_IMATRIX( catch_n );                // Info for fishery biomass; columns = Month
  DATA_MATRIX( catch_obs );               // Observed fishery catch (mt or thousands of fish per the fleet's Observation_units) and log_sd; columns = Observation, Error

  // -- 2.4.2 Survey components
  DATA_IMATRIX( index_ctl );              // Index metadata; columns = Fleet_code, Species, Year
  DATA_MATRIX( index_n );                 // Info for index; columns = Month
  DATA_MATRIX( index_obs );               // Observed index and log_sd; columns = Observation, Error
  DATA_VECTOR( index_log_q_prior );        // Prior mean for catchability
  DATA_IVECTOR( index_ll_type );          // Survey index likelihood family per fleet (0 = lognormal IID, 1 = MVN bare quadratic form, 2 = MVNORM full density, 3 = natural-scale normal with absolute sd)
  DATA_STRUCT( index_cov_mat, LOM_t );    // Per-fleet covariance matrices Sigma for the MVN/MVNORM survey likelihood (1x1 dummy if unused)
  DATA_VECTOR( index_cov_const );         // Per-fleet 0.5*(logdet(Sigma) + n*log(2pi)); subtracted for "MVN" so it reports the bare quadratic form

  // -- 2.4.2b One-step-ahead (OSA) residual support
  // `obsvec` is every observation entering the likelihood, flat: log catch, log index,
  // each comp/caal row's bin counts, each stomach's diet composition. `keep` is the
  // indicator TMB::oneStepPredict() toggles element by element; it defaults to all ones,
  // so ordinary fitting is unchanged. Each `*_obsvec_idx` gives its `*_obs` row's 0-based
  // position in `obsvec` (for a composition, where the row's bins start), or -1 when the
  // row is not fitted. osa_mode 1 switches the comp/caal/diet branches to a proper,
  // unweighted, keep-gated density; mode 2 adds the conditional CDF terms
  // oneStepPredict(method = "cdf") reads -- log F(x) gated by keep.cdf_lower,
  // log(1 - F(x)) by keep.cdf_upper -- computed only in mode 2, so no other path pays
  // for the pbeta/pnorm. The Dirichlet-multinomial families have no CDF term
  // (comp_osa.hpp) and osa_residuals() keeps those fleets off this method.
  DATA_VECTOR( obsvec );                    // Flat observations for OSA residuals
  DATA_VECTOR_INDICATOR( keep, obsvec );    // oneStepPredict keep indicator (defaults to 1 when fitting)
  DATA_IVECTOR( catch_obsvec_idx );         // obsvec position for each catch_obs row (-1 = excluded)
  DATA_IVECTOR( index_obsvec_idx );         // obsvec position for each index_obs row (-1 = excluded)
  DATA_IVECTOR( comp_obsvec_idx );          // obsvec start position for each comp_obs row's bins (-1 = excluded)
  DATA_IVECTOR( caal_obsvec_idx );          // obsvec start position for each caal_obs row's bins (-1 = excluded)
  DATA_IVECTOR( diet_obsvec_idx );          // obsvec start position for each stomach's prey bins (incl. "other prey"); length n_stomach_obs (-1 = excluded)
  DATA_INTEGER( osa_mode );                 // 0 = normal fitting (default); 1 = OSA build (unweighted keep-gated comp/caal/diet densities); 2 = OSA build + conditional CDF terms for method = "cdf"
  DATA_INTEGER( adreport_sel );             // 0 = no selectivity standard errors (default); 1 = ADREPORT log_sel_at_age; set via fit_control(selectivity_se =)
  DATA_SCALAR( comp_offset );               // proportion offset added to comp/caal/diet obs & pred before the multinomial; set via rearrange_data()/fit_control()

  // -- 2.4.2c. Simulation switches (sim_mod(simulate = TRUE))
  // Observations are always drawn; process error is a choice, so it is switched.
  //
  // simulate_state: which process. Slots 0 recruitment, 1 M, 2 growth, 3 catchability,
  // 4 selectivity -- the codes linkage_process uses, so a random linkage gates on
  // simulate_state(linkage_process(i)) untranslated. Default off: redrawing a process
  // changes what self_test() measures, from recovering parameters to recovering a
  // process. Recruitment covers init_dev too, the initial age structure being recruitment
  // pre-styr.
  //
  // simulate_period: slot 0 the fitted window (0 < Year <= endyr), slot 1 outside it.
  // Only slot 0 is read, and only by the PROCESS draws, whose densities cover the
  // hindcast alone. Observations are drawn in every period deliberately, because
  // run_mse() splices the negative-Year rows back in as the next assessment's data. Slot
  // 1 is reserved and inert; sim_mod() exposes neither slot.
  DATA_IVECTOR( simulate_state );
  DATA_IVECTOR( simulate_period );

  // -- 2.4.3. Composition data
  DATA_IMATRIX( comp_ctl );               // Age/length comp metadata; columns = Fleet_code, Species, Sex, Age0_Length1, Year
  DATA_MATRIX( comp_n );                  // Month and sample size on observed age/length comp; columns = Month, Sample size
  DATA_MATRIX( comp_obs );                // Observed age/length comp; cols = Comp_1, Comp_2, etc. can be proportion
  DATA_IMATRIX( caal_ctl );               // CAAL metadata; columns = Fleet_code, Species, Sex, Year, Length_bin
  DATA_MATRIX( caal_n );                  // Sample size on CAAL; ONE column = Sample size. Unlike comp_n there is no month: a CAAL observation is placed at its fleet's Month, and the age-length key is annual.

  // -- 2.4.4. Initial equilibrium catch: the catch yielded under the initial F,
  // in the year before the hindcast. 0 rows unless the data supply one.
  DATA_IMATRIX( equil_catch_ctl );        // columns = Fleet_code, Species
  DATA_MATRIX( equil_catch_obs );         // columns = Catch, Log_sd
  DATA_MATRIX( caal_obs );                // Observed CAAL; cols = Comp_1, Comp_2, etc. can be proportion

  // -- 2.4.5 Age and selectivity
  DATA_IMATRIX( emp_sel_ctl );            // Empirical selectivity metadata; columns = Fleet_code, Species, Sex, Year
  DATA_MATRIX( emp_sel_obs );             // Observed emprical fishery selectivity; columns = Compe_1, Comp_2, etc.
  DATA_ARRAY( age_trans_matrix);          // observed sp_age/size compositions; n = [nspp, nages, index_age_bins]
  DATA_ARRAY( age_error );                // Aging error matrices; n = [n_matrices, nages, nages]. One per species unless fleet_control$Ageing_error_index selects otherwise.

  // -- 2.3.5. Growth
  DATA_ARRAY( weight_obs );               // Weight-at-age by year; n = [nweight, sex, nages, nyrs]

  // 2.3.6. Diet data
  DATA_IVECTOR(diet_ll_type);             // Diet composition log-likelihood type
  DATA_VECTOR( fday );                    // number of foraging days for each predator
  DATA_ARRAY( ration_data );              // Relative-foraging rate
  DATA_MATRIX( diet_obs );                // Pred, prey, pred-age, prey-age for diet matrix (weight of prey in pred stomach)
  DATA_IMATRIX( diet_ctl );               // Per diet_obs row: predator, prey, both sexes, both ages, year
  DATA_INTEGER(n_stomach_obs);            // The total number of unique stomach samples (groups)
  DATA_IVECTOR(stomach_id);               // A vector mapping each diet data row to a stomach ID

  // 2.3.7. Environmental data
  DATA_MATRIX( env_index );               // Matrix of environmental predictors such as bottom temperature

  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( other_food );              // Biomass of other prey (mt)
  DATA_VECTOR( Pvalue );                  // Proportion of Cmax; the P in Cmax*fT*Pvalue*PAge
  DATA_IVECTOR( Ceq );                    // Which consumption equation to use; currently all sp = 1
  DATA_IVECTOR( Cindex );                 // Which env_index column drives bioenergetics.
  DATA_VECTOR( CA );                      // Wt specific intercept of Cmax=CA*W^CB
  DATA_VECTOR( CB );                      // Wt specific slope of Cmax=CA*W^CB
  DATA_VECTOR( Qc );                      // used in fT, QC value
  DATA_VECTOR( Tco );                     // used in fT, thermal optimum
  DATA_VECTOR( Tcm );                     // used in fT, thermal max
  DATA_VECTOR( Tcl );                     // used in fT eq 3, limit
  DATA_VECTOR( CK1 );                     // used in fT eq 3, limit where C is .98 max (ascending)
  DATA_VECTOR( CK4 );                     // used in fT eq 3, temp where C is .98 max (descending)

  // -- 2.4.3. Others
  DATA_MATRIX( sex_ratio );               // Proportion-at-age of females of population
  DATA_MATRIX( maturity );                // Proportion of mature females at age; [nspp, nages]
  DATA_SCALAR(bias_adjust_obs);           // Whether to apply bias adjustment (-sigma^2/2) in index likelihoods
  DATA_SCALAR(bias_adjust_proc);           // Whether to apply bias adjustment (-sigma^2/2) in process likelihoods


  /** ------------------------------------------------------------------------ //
   * 3. PARAMETER SECTION                                                      //
   * ------------------------------------------------------------------------- */

  PARAMETER( dummy );                             // Tests derived quantities given input parameters; n = [1]
  PARAMETER_VECTOR( log_pop_scalar );              // Log multiplier on input numbers-at-age, per species (estDynamics = 2)

  // -- 3.1. Recruitment parameters
  PARAMETER_MATRIX( rec_pars );                   // Stock-recruit parameters: col1 = mean rec, col2 = SRR alpha, col3 = SRR beta
  PARAMETER_VECTOR( R_log_sd );                    // Standard deviation of recruitment deviations
  PARAMETER_MATRIX( rec_dev );                    // Annual recruitment deviation; n = [nspp, nyrs]

  // DSEM parameters. All zero-length and mapped out when dsem_on = 0, so they
  // add nothing to the objective. rec_dev and R_log_sd above stay parameters
  // either way -- under a DSEM the deviations are overwritten from the latent
  // states rather than estimated, which is what keeps the non-DSEM path intact.
  PARAMETER_VECTOR( dsem_beta_z );                 // SEM path coefficients
  PARAMETER_VECTOR( dsem_lnsigma_z );              // log measurement SDs
  PARAMETER_VECTOR( dsem_mu_j );                   // per-variable mean
  PARAMETER_VECTOR( dsem_delta0_j );               // initial-condition offset (may be empty)
  PARAMETER_ARRAY( dsem_x_tj );                    // latent states [nyrs_dsem, n_var]
  PARAMETER_MATRIX( init_dev );                   // Log initial numbers-at-age under initMode 0, log deviations off equilibrium under 2/3/4/6; n = [nspp, max_age]

  // -- 3.2. Natural mortality (M1)
  PARAMETER_ARRAY( log_M1 );                       // Natural mortality: residual under multispecies mode, total under single species; n = [nspp, nsex, nages]
  PARAMETER_ARRAY( log_M1_dev );                   // Natural mortality annual deviate; n = [nspp, nsex, nyrs]
  PARAMETER_ARRAY( M1_beta );                     // Regression coefficients for environmnetally linked M1; n = [nspp, nsex, n indices]
  PARAMETER_ARRAY( M1_rho );                      // Correlation for AR1 random effects on age and year; n = [nspp, nsex, 2]
  PARAMETER_ARRAY( M1_dev_log_sd );                // Innovation sd of the M1 random effects, one per species/sex; n = [nspp, nsex]

  // -- 3.3. Growth
  PARAMETER_ARRAY(log_growth_pars);                // Mean growth curve parameters [sp, sex, par]
  // Time-varying growth comes from growth_linkage_offset below.
  PARAMETER_ARRAY(growth_log_sd);                  // Log standard deviation of length-at- min and max age [sp, sex, 2]
  PARAMETER_MATRIX(weight_length_pars);           // Length-weight parameters [sp, (alpha, beta)]

  // -- 3.3b. Linkage-table coefficients (see linkage.hpp). Aligned row-for-row with the
  //          linkage_* DATA vectors above. May be length 0 (no linkages supplied).
  PARAMETER_VECTOR(beta_linkage);

  // Random-effect linkage coefficients and their hyperparameters. All length 0 unless a
  // `~ (1|group)` / rw() / ar1() linkage was supplied; when present, beta_linkage_re
  // enters the Laplace approximation (added to random=) and its density goes to jnll_comp
  // row 20.
  PARAMETER_VECTOR(beta_linkage_re);     // RE deviation coefficients (Laplace-integrated)
  // Penalized RE deviations: same row-20 density as beta_linkage_re, but NOT in the
  // `random` set, so the density acts as a plain penalty instead of being integrated out.
  // Two vectors rather than one flag because TMB selects random effects by parameter
  // NAME, and one model may mix both. Permitted only with a fixed SD -- estimating
  // deviations and their SD jointly as fixed effects is degenerate (see
  // .re_group_table()).
  PARAMETER_VECTOR(beta_linkage_re_pen);
  PARAMETER_VECTOR(log_sigma_linkage);   // one log-SD per RE group
  PARAMETER_VECTOR(trans_rho_linkage);   // one transformed rho per AR1 group
  PARAMETER_VECTOR(beta_linkage_obs);    // Rogers QAR1 effect size: one per observed ar1 group
  PARAMETER_VECTOR(log_obs_sd_linkage);  // Rogers QAR1 observation log-SD: one per observed ar1 group (estimated)

  // -- 3.4. Fishing mortality parameters
  PARAMETER_VECTOR( log_Flimit );                  // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_VECTOR( log_Ftarget );                 // Target fishing mortality for projections on log scale; n = [nspp, nyrs]
  PARAMETER_VECTOR( log_Finit );                   // Fishing mortality for initial population to induce non-equilibrium; n = [nspp]
  PARAMETER_VECTOR( proj_F_prop );                // Proportion of fishing mortality from each fleet for projections; n = [n_fsh]
  PARAMETER_MATRIX( log_F );                       // Annual fishing mortality; n = [n_fsh, nyrs]

  // -- 3.4. Survey catchability parameters
  PARAMETER_VECTOR( index_log_q );                 // Survey catchability; n = [n_index]
  PARAMETER_MATRIX( index_q_beta );               // Survey catchability regression coefficients on env_data columns
  PARAMETER_MATRIX( index_q_dev );                // Annual survey catchability deviates; n = [n_index, nyrs_hind]
  PARAMETER_VECTOR( index_q_log_sd );              // Log standard deviation of prior on survey catchability; n = [1, n_index]
  PARAMETER_VECTOR( index_q_dev_log_sd );// Log standard deviation of time varying survey catchability; n = [1, n_index]

  // -- 3.5. Selectivity parameters
  PARAMETER_ARRAY( sel_coff );                    // selectivity parameters for non-parametric; n = [n_selectivities, nsex, n_sel_bins]
  PARAMETER_ARRAY( sel_coff_dev );                // Annual deviates for non-parametric selectivity parameters; n = [n_selectivities, nsex, n_sel_bins]
  PARAMETER_ARRAY( log_sel_slp );                  // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( sel_inf );                     // selectivity paramaters for logistic; n = [2, n_selectivities, nsex]
  PARAMETER_ARRAY( log_sel_slp_dev );              // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( sel_inf_dev );                 // selectivity parameter deviate for logistic; n = [2, n_selectivities, nsex, n_sel_blocks]
  PARAMETER_ARRAY( log_sel_apical );              // per-sex log multiplier on the whole curve, after the form and before normalization; n = [n_selectivities, nsex]
  PARAMETER_ARRAY( sel_dn6 );                     // DoubleNormalSS3 (SS3 pattern 24) parameters on SS3's scales: peak, logit top, log asc, log desc, logit init, logit final; n = [6, n_selectivities, nsex]
  PARAMETER_VECTOR( sel_dev_log_sd );              // Log standard deviation of selectivity; n = [1, n_selectivities]
  PARAMETER_MATRIX( sel_curve_pen );              // Non-parametric selectivity penalty columns (Sel_curve_pen1-3), read differently per form -- see the penalty blocks in section 13.1

  // -- 3.6. Data variance
  //FIXME: remove log_sd terms below
  PARAMETER_VECTOR( index_log_sd );                // Log standard deviation of survey index time-series; n = [1, n_index]
  PARAMETER_VECTOR( catch_log_sd );                // Log standard deviation of fishery catch time-series; n = [1, n_fsh]
  PARAMETER_VECTOR( comp_weights );               // Weights for composition data
  PARAMETER_VECTOR( caal_weights );               // Weights for CAAL data
  vector<Type>  DM_pars_comp = exp(comp_weights); // Dirichlet-multinomial scalars
  vector<Type>  DM_pars_caal = exp(caal_weights); // Dirichlet-multinomial scalars

  PARAMETER_VECTOR( diet_comp_weights );          // Weights for diet composition data
  vector<Type>  DM_diet_pars = exp(diet_comp_weights);// Dirichlet-multinomial scalars

  // 3.8. Gamma selectivity parameters
  PARAMETER_VECTOR( log_gam_a );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  PARAMETER_VECTOR( log_gam_b );                  // Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10

  // 3.9. Preference
  PARAMETER_MATRIX( log_phi );                    // Species preference coefficient; n = [nspp, nspp]


  /** ------------------------------------------------------------------------ //
   * 4. DERIVED QUANTITIES SECTION ----                                        //
   * ------------------------------------------------------------------------- */

  // 4.1. Derived indices
  int n_flt = flt_type.size();
  vector<int> joint_adjust(comp_obs.rows()); joint_adjust.setZero();
  Type penalty = 0.0;
  Type ricker_intercept = 0.0;

  // -- 4.2. Growth
  array<Type> growth_matrix(nspp * 2 + n_flt, max_sex, max_age, max_nlengths_pop, nyrs); growth_matrix.setZero(); // Growth transition matrix on the POPULATION length bins, per fleet and per species derived quantity (biomass, ssb)
  array<Type> weight_hat(nspp * 2 + n_flt, max_sex, max_age, nyrs); weight_hat.setZero(); // Estimated weight-at-age, per fleet and per species derived quantity (biomass, ssb)
  array<Type> mat_weight_hat(nspp * 2 + n_flt, max_sex, max_age, nyrs); mat_weight_hat.setZero(); // Mature weight-at-age (kg), maturity-at-length species only
  array<Type> length_hat(nspp * 2 + n_flt, max_sex, max_age, nyrs); length_hat.setZero(); // Estimated length-at-age, per fleet and per species derived quantity (biomass, ssb)

  // -- 4.3. Estimated population quantities
  vector<Type>  pop_scalar = exp(log_pop_scalar);                                   // Multiplier on input numbers-at-age, per species
  vector<Type>  avg_R(nspp); avg_R.setZero();                                       // Mean recruitment of hindcast
  matrix<Type>  R_hat(nspp, nyrs); R_hat.setZero();                                 // Expected recruitment given SR curve
  matrix<Type>  mort_sum(nspp, max_age); mort_sum.setZero();
  // Fishery selectivity at age used by initMode 6 to spread Finit over ages;
  // filled in section 5.12.1, and 1 everywhere for every other mode.
  array<Type>   sel_init(nspp, max_sex, max_age); sel_init.setZero();
  // Initial age structure WITHOUT the initial deviates: the equilibrium they are applied
  // to. Only the equilibrium catch reads it, and that is the quantity it is defined on --
  // SS3 computes the equilibrium catch from the deviation-free equilibrium, not from the
  // year-1 numbers. Filled alongside N_at_age in section 6.5 so the two cannot drift
  // apart.
  array<Type>   N_eq(nspp, max_sex, max_age); N_eq.setZero();
  matrix<Type>  R0(nspp, nyrs); R0.setZero();                                       // Equilibrium recruitment at F = 0.
  matrix<Type>  alpha(nspp, nyrs); alpha.setZero();                                 // Stock recruit alpha
  matrix<Type>  Beta(nspp, nyrs); Beta.setZero();                                   // Stock recruit beta
  matrix<Type>  steepness(nspp, nyrs); steepness.setZero();                         // Expected % of R0 at 20% SSB0.
  vector<Type>  R_init(nspp); R_init.setZero();                                     // Equilibrium recruitment at F = Finit (non-equilibrium).
  matrix<Type>  R(nspp, nyrs); R.setZero();                                         // Estimated recruitment (thousands of fish; n-at-age is in thousands)
  array<Type>   biomass_at_age(nspp, max_sex, max_age, nyrs); biomass_at_age.setZero();// Estimated biomass-at-age (mt; thousands of fish x kg)
  matrix<Type>  biomass(nspp, nyrs); biomass.setZero();                             // Estimated biomass (mt)
  matrix<Type>  exploitable_biomass(nspp, nyrs); exploitable_biomass.setZero();     // Estimated exploitable biomass (mt)
  matrix<Type>  ssb(nspp, nyrs); ssb.setZero();                                     // Estimated spawning stock biomass (mt)
  matrix<Type>  biomass_depletion(nspp, nyrs); biomass_depletion.setZero();         // Estimated biomass biomass_depletion
  matrix<Type>  ssb_depletion(nspp, nyrs); ssb_depletion.setZero();                 // Estimated biomass_depletion of spawning stock biomass
  array<Type>   M_at_age(nspp, max_sex, max_age, nyrs); M_at_age.setZero();         // Total natural mortality at age
  array<Type>   M_at_age_dB0(nspp, max_sex, max_age, nyrs); M_at_age_dB0.setZero(); // Total natural mortality at age (dynamic B0)
  array<Type>   M_at_age_dBF(nspp, max_sex, max_age, nyrs); M_at_age_dBF.setZero(); // Total natural mortality at age (dynamic BF)
  array<Type>   M1_at_age(nspp, max_sex, max_age, nyrs); M1_at_age.setZero();       // Residual or total natural mortality at age
  array<Type>   N_at_age(nspp, max_sex, max_age, nyrs); N_at_age.setZero();         // Numbers at age
  array<Type>   avgN_at_age(nspp, max_sex, max_age, nyrs); avgN_at_age.setZero();   // Average numbers-at-age
  array<Type>   avgN_at_age_dB0(nspp, max_sex, max_age, nyrs); avgN_at_age_dB0.setZero();   // Average numbers-at-age (dynamic B0)
  array<Type>   avgN_at_age_dBF(nspp, max_sex, max_age, nyrs); avgN_at_age_dBF.setZero();   // Average numbers-at-age (dynamic BF)
  array<Type>   Z_at_age(nspp, max_sex, max_age, nyrs); Z_at_age.setZero();         // Total mortality at age
  vector<Type>  R_sd(nspp); R_sd.setZero();                                         // Standard deviation of recruitment variation
  vector<Type>  zero_N_pen(nspp); zero_N_pen.setZero();                             // Additional penalty to add to likelihood if n-at-age goes < 0

  // -- 4.4. Selectivity parameters
  array<Type>   sel_at_age(n_flt, max_sex, max_age, nyrs); sel_at_age.setZero();    // Estimated selectivity at age
  array<Type>   sel_at_length(n_flt, max_sex, max_nlengths_pop, nyrs); sel_at_length.setZero();// Estimated selectivity at length, on the POPULATION length bins
  array<Type>   avg_sel(n_flt, max_sex, nyrs_hind); avg_sel.setZero();              // Average selectivity for non-parametric up to n_sel_bins
  array<Type>   non_par_sel(n_flt, max_sex, max_bin, nyrs); non_par_sel.setZero();  // Estimated selectivity for AMAK non-parametric (pre-normalization)
  array<Type>   log_non_par_sel(n_flt, max_sex, max_bin, nyrs); log_non_par_sel.setZero(); // Same curve on the log scale, where the penalties are written (exp() underflows below about -745, so it is carried, not recovered)
  vector<Type>  sel_dev_sd(n_flt); sel_dev_sd.setZero();                            // Standard deviation of selectivity deviates

  // -- 4.5. Fishery components
  matrix<Type>  F_spp(nspp, nyrs); F_spp.setZero();                                 // Fully selected fishing mortality by species
  matrix<Type>  F_flt(n_flt, nyrs); F_flt.setZero();                                // Fully selected fishing mortality by fleet
  array<Type>   F_flt_age(n_flt, max_sex, max_age, nyrs); F_flt_age.setZero();            // Estimated fishing mortality-at-age/sex for each fishery
  array<Type>   Flimit_at_age(nspp, max_sex, max_age, nyrs); Flimit_at_age.setZero();   // Estimated limit fishing mortality-at-age/sex for each species (FOFL proxy)
  array<Type>   Ftarget_at_age(nspp, max_sex, max_age, nyrs); Ftarget_at_age.setZero(); // Estimated target fishing mortality-at-age/sex for each species (max FABC proxy)
  array<Type>   F_at_age(nspp, max_sex, max_age, nyrs); F_at_age.setZero();       // Sum of annual estimated fishing mortalities for each species-at-age
  vector<Type>  catch_hat(catch_obs.rows()); catch_hat.setZero();                   // Estimated fishery yield/numbers (mt, or thousands of fish)
  vector<Type>  max_catch_hat(catch_obs.rows()); max_catch_hat.setZero();           // Estimated exploitable biomass/numbers by fleet (mt, or thousands of fish)
  // As index_sd below, for the fishery. Catch is lognormal throughout, so this is
  // always a log-scale sd (unitless), not a log of anything. rename_output()
  // also writes it as log_catch_sd, the name fits saved before 5.9.0 carry.
  vector<Type>  catch_sd(catch_obs.rows()); catch_sd.setZero();
  // Concentrated (analytical) catch sd used when Estimate_catch_sd = 2; a
  // log-scale sd (unitless), not a log. Mirrors index_analytical_sd below.
  vector<Type>  catch_analytical_sd(n_flt); catch_analytical_sd.setZero();
  vector<Type>  catch_n_obs(n_flt); catch_n_obs.setZero();  // fitted catch rows per fleet

  // -- 4.6. Biological reference points
  array<Type>   NByage0(nspp, max_sex, max_age, nyrs); NByage0.setZero();                 // Numbers at age at mean recruitment and F = 0
  array<Type>   NByageF(nspp, max_sex, max_age, nyrs); NByageF.setZero();                 // Numbers at age at mean recruitment and F = Ftarget; the age structure behind SBF
  array<Type>   N_at_age_dB0(nspp, max_sex, max_age, nyrs); N_at_age_dB0.setZero();   // Numbers at age at F = 0 (accounts for annual recruitment)
  array<Type>   N_at_age_dBF(nspp, max_sex, max_age, nyrs); N_at_age_dBF.setZero();   // Female numbers at age at F = Ftarget (accounts for annual recruitment)
  matrix<Type>  DynamicSB0(nspp, nyrs); DynamicSB0.setZero();                       // Estimated dynamic spawning biomass at F = 0 (accounts for S-R curve)
  matrix<Type>  DynamicB0(nspp, nyrs); DynamicB0.setZero();                         // Estimated dynamic  biomass at F = 0 (accounts for S-R curve)
  matrix<Type>  DynamicSBF(nspp, nyrs); DynamicSBF.setZero();                       // Estimated dynamic spawning biomass at F = Ftarget (accounts for S-R curve)
  // Zeroed explicitly, like every other array here: the recursion writes only
  // ages below a species' own nages, and nothing at all under predation, so the
  // rest is read from the REPORT as whatever it was initialised to.
  array<Type>   NbyageSPR(4, nspp, max_age); NbyageSPR.setZero();                   // Numbers at age per recruit for the SPR reference points; slots 0-3 are F = 0, Flimit, Ftarget, Finit
  vector<Type>  SPRlimit(nspp); SPRlimit.setZero();                                 // Estimated Plimit SPR
  vector<Type>  SPRtarget(nspp); SPRtarget.setZero();                               // Estimated Ptarget SPR
  vector<Type>  SPR0(nspp); SPR0.setZero();                                         // Estimated spawning biomass per recruit at F = 0
  vector<Type>  SPRFinit(nspp); SPRFinit.setZero();                                 // Estimated spawning biomass per recruit at Finit
  matrix<Type>  SB0(nspp, nyrs); SB0.setZero();                                     // Estimated spawning stock biomass at F = 0 (Accounts for S-R)
  matrix<Type>  SBF(nspp, nyrs); SBF.setZero();                                     // Estimated spawning stock biomass at F = target (Accounts for S-R)
  matrix<Type>  B0(nspp, nyrs); B0.setZero();                                       // Estimated biomass at F = 0 (Accounts for S-R)
  vector<Type>  Flimit = exp(log_Flimit);                                            // Limit F (FOFL proxy, e.g. F35%) on natural scale
  vector<Type>  Ftarget = exp(log_Ftarget);                                          // Target F (max FABC proxy, e.g. F40%) on natural scale
  vector<Type>  Finit = exp(log_Finit);                                              // Initial F for non-equilibrium age-structure
  vector<Type>  index_q_mult(index_q_beta.cols()); index_q_mult.setZero();          // Environmental design matrix for q
  vector<Type>  M1_mult(M1_beta.dim(2)); M1_mult.setZero();                         // Environmental design matrix for M1
  vector<Type>  beta_q_tmp(index_q_beta.cols()); beta_q_tmp.setZero();              // Q beta parameters by species, for the matrix multiply
  vector<Type>  env_q_tmp(index_q_beta.cols()); env_q_tmp.setZero();                // Q env data by year, for the matrix multiply
  vector<Type>  beta_M1_tmp(M1_beta.dim(2)); beta_M1_tmp.setZero();                 // Beta parameters by species, for the matrix multiply
  vector<Type>  env_M1_tmp(M1_beta.dim(2)); env_M1_tmp.setZero();                   // Env data by year, for the matrix multiply
  matrix<Type>  proj_F(nspp, nyrs); proj_F.setZero();                               // Projected F (Fabc/Ftac/etc) using harvest control rule

  // -- 4.7. Survey components
  vector<Type>  index_q_sd(n_flt); index_q_sd.setZero();                            // Vector of standard deviation of survey catchability prior
  vector<Type>  index_q_dev_sd(n_flt); index_q_dev_sd.setZero();                    // Vector of standard deviation of time-varying survey catchability deviation
  vector<Type>  index_hat(index_obs.rows()); index_hat.setZero();                   // Estimated survey index; mt or thousands of fish per the fleet's Observation_units
  // Rejection bookkeeping for the natural-scale survey draws (sim_mod only), per
  // index_obs row: a non-zero tries count marks the rows a rejection-capable branch drew
  // -- the untruncated "Normal" family and MVN/MVNORM, which must redraw a non-positive
  // index. "TruncatedNormal" draws by inverse CDF and never enters here. sim_mod() sizes
  // the draw/density gap from index_hat and index_sd, not from these counts.
  vector<Type>  index_trunc_tries_sim(index_obs.rows());   index_trunc_tries_sim.setZero();
  vector<Type>  index_trunc_rejects_sim(index_obs.rows()); index_trunc_rejects_sim.setZero();
  // The marginal sd the draw itself used, so sim_mod() can size the truncated mass as
  // Phi(-index_hat/sd). Separate from index_sd because a covariance fleet draws from
  // index_cov_mat, leaving its index_sd whatever the (unused) Index_sd column carried.
  vector<Type>  index_trunc_sd_sim(index_obs.rows());      index_trunc_sd_sim.setZero();
  // The redraw budget ACTUALLY applied to each row, which is what sim_mod() quotes: the
  // constant below on the univariate branch, and that scaled by the fleet's row count on
  // the correlated one, which rejects the whole vector at once, so its joint rejection
  // probability climbs with the row count.
  vector<Type>  index_trunc_budget_sim(index_obs.rows()); index_trunc_budget_sim.setZero();
  const int index_trunc_max_tries = 100;
  // The observation sd the index likelihood actually used for each row, whichever of the
  // three est_sigma_index routes supplied it. NOT a log: the sd itself, on the scale that
  // fleet's Index_distribution works on -- a log-scale sd (unitless) for Lognormal, an
  // ABSOLUTE sd in the units of the index for the natural-scale families. rename_output()
  // also writes it as log_index_sd, the name fits saved before 5.9.0 carry.
  vector<Type>  index_sd(index_obs.rows()); index_sd.setZero();
  vector<Type>  index_analytical_sd(n_flt); index_analytical_sd.setZero();          // Concentrated (analytical) index sd used when Estimate_index_sd = 2; a log-scale sd, not a log
  vector<Type>  index_q_analytical(n_flt); index_q_analytical.setZero();            // Analytical catchability (not an sd) used when Catchability = 3; "AnalyticalArith" (7) has its own block at 8.2b
  matrix<Type>  index_q(n_flt, nyrs_hind); index_q.setZero();                       // Estimated survey catchability //FIXME: extend out to full time-series
  vector<Type>  index_n_obs(n_flt); index_n_obs.setZero();                          // Observations per survey time series

  // -- 4.8. Composition data
  // age_hat / age_obs_hat are indexed by AGE whatever dimension the fleet's composition
  // data are on -- a length composition is built at age, then converted through the
  // age-length transition -- so the width is the widest age index any row writes:
  // nages(sp), or nages(sp) * 2 for a joint-sex row (comp_ctl column 2 == 3), never
  // narrower than comp_obs. Sizing from comp_obs alone gives only the workbook's Comp_
  // columns, i.e. nlengths on a length-only model, and Eigen does not bounds-check in a
  // release build, so nlengths < nages would write silently into adjacent memory. Species
  // is range-checked because nothing between rearrange_data() and MakeADFun validates
  // comp_data$Species, and a non-numeric one arrives as NA_integer_.
  int max_age_cols = static_cast<int>(comp_obs.cols());
  for(int comp_row = 0; comp_row < comp_ctl.rows(); comp_row++){
    int comp_row_sp = comp_ctl(comp_row, 1) - 1;
    if((comp_row_sp < 0) | (comp_row_sp >= nspp)){
      error("comp_data Species on row %d is not in 1..nspp", comp_row + 1);
    }
    int comp_row_cols = (comp_ctl(comp_row, 2) == 3) ? nages(comp_row_sp) * 2 : nages(comp_row_sp);
    if(comp_row_cols > max_age_cols){ max_age_cols = comp_row_cols; }
  }
  vector<Type>  n_hat(comp_obs.rows()) ; n_hat.setZero() ;                          // Estimated catch (numbers)
  matrix<Type>  age_hat(comp_obs.rows(), max_age_cols); age_hat.setZero();          // Estimated catch at true age
  matrix<Type>  age_obs_hat(comp_obs.rows(), max_age_cols); age_obs_hat.setZero();  // Estimated catch at observed age (accounts for ageing error)
  matrix<Type>  comp_hat = comp_obs; comp_hat.setZero();                            // Estimated comp
  matrix<Type>  caal_hat = caal_obs; caal_hat.setZero();                            // Estimated CAAL
  array<Type>   pred_CAAL(n_flt, max_sex, max_age, max_nlengths, nyrs); pred_CAAL.setZero(); // Predicted CAAL for each fleet

  // -- 4.9. Ration components
  array<Type>   consumption_at_age( nspp, max_sex, max_age, nyrs ); consumption_at_age.setZero(); // Annual ration at age (kg/yr)
  matrix<Type>  fT( nspp, nyrs ); fT.setZero();                                     // Pre-allocation of temperature function of consumption

  // -- 4.10. Diet components
  array<Type>   diet_prop(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); diet_prop.setZero();             // Stomach proportion by weight U
  array<Type>   diet_prop_hat(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); diet_prop_hat.setZero();     // Predicted stomach proportion by weight U
  array<Type>   other_food_diet_prop(nspp, max_sex, max_age, nyrs); other_food_diet_prop.setZero();                 // Other food diet proportion by weight
  matrix<Type>  diet_hat = diet_obs; diet_hat.setZero();                                                            // Estimated stomach proportion by weight U (formated following data input)

  // -- 4.11. Suitability components
  array<Type>   avail_food(nspp, max_sex, max_age, nyrs); avail_food.setZero();                               // Available food to predator
  array<Type>   avail_food_dB0(nspp, max_sex, max_age, nyrs); avail_food_dB0.setZero();                       // Available food to predator (dynamic B0)
  array<Type>   avail_food_dBF(nspp, max_sex, max_age, nyrs); avail_food_dBF.setZero();                       // Available food to predator (dynamic BF)
  array<Type>   stom_div_bio(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); stom_div_bio.setZero(); // Stomach proportion over biomass; U/ (W * N)
  array<Type>   suitability(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); suitability.setZero();   // Suitability/gamma selectivity of predator age u on prey age a
  array<Type>   suit_other(nspp, max_sex, max_age, nyrs); suit_other.setZero();                               // Suitability not accounted for by the included prey
  array<Type>   other_diet_prop_hat(nspp, max_sex, max_age, nyrs); other_diet_prop_hat.setZero();             // Diet of prey not included in the model

  // -- 4.12. Suitability parameters
  vector<Type> gam_a = exp(log_gam_a);                                    // Predator size-selectivity: shape parameter for gamma suitability, mean for normal of logs
  vector<Type> gam_b = exp(log_gam_b);                                    // Predator size-selectivity: scale parameter for gamma suitability, sd for normal of logs
  matrix<Type> vulnerability(nspp, nspp); vulnerability.setZero();        // Predator-prey preference coefficients
  vector<Type> vulnerability_other(nspp); vulnerability_other.setZero();  // Preference for other food

  // -- 4.13. Predation components
  array<Type>   M2_at_age(nspp, max_sex, max_age, nyrs); M2_at_age.setZero();                        // Predation mortality at age
  array<Type>   M2_at_age_dB0(nspp, max_sex, max_age, nyrs); M2_at_age_dB0.setZero();                // Predation mortality at age (dynamic B0)
  array<Type>   M2_at_age_dBF(nspp, max_sex, max_age, nyrs); M2_at_age_dBF.setZero();                // Predation mortality at age (dynamic BF)
  array<Type>   M2_prop(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); M2_prop.setZero();  // Each predator-at-age's contribution to prey M2 (yr^-1), not a share of it
  array<Type>   B_eaten(nspp * max_sex, nspp * max_sex, max_age, max_age, nyrs); B_eaten.setZero();  // Biomass of prey eaten via predation by a predator at age
  array<Type>   B_eaten_as_prey(nspp, max_sex, max_age, nyrs); B_eaten_as_prey.setZero();            // Biomass eaten as prey via predation

  /** ------------------------------------------------------------------------ //
   * 5. INITIAL CALCULATIONS                                                   //
   * ------------------------------------------------------------------------- */

  // 5.1. DATA VARIANCE TERMS
  R_sd = exp(R_log_sd); // Convert log sd to natural scale
  sel_dev_sd = exp(sel_dev_log_sd) ;
  index_q_sd = exp(index_q_log_sd) ;
  index_q_dev_sd = exp(index_q_dev_log_sd) ;
  Cindex -=1; // Subtract 1 from Cindex to deal with indexing start at 0


  // -- 5.2. SIMULATE LINKAGE RANDOM EFFECTS (sim_mod(simulate = TRUE))
  // Time-varying recruitment, M, growth, q and selectivity are all written as a random
  // linkage, so all five draw here, before the processes read them. Each structure draws
  // from the density that scores it in section 14.6b, slots being in ascending time order
  // so "next slot" is "next year": a rw group leaves slot 0 alone, a level the density
  // never sees and the map pins, and an ar1 draw is standardized then scaled by the
  // MARGINAL SD sigma.
  //
  // beta_linkage_re_drawn_sim marks which slots this call wrote. The reported deviation
  // vector spans EVERY random linkage, drawn or not, so without the mask a caller would
  // score recovery on slots left at their fitted values.
  SIMULATE {
    vector<Type> beta_linkage_re_drawn_sim(linkage_re_sigma.size());
    beta_linkage_re_drawn_sim.setZero();
    if (log_sigma_linkage.size() > 0 && simulate_period(0) == 1) {
      int n_slot = linkage_re_sigma.size();
      int n_grp  = log_sigma_linkage.size();

      // Drawn only if EVERY row using the group names a process the caller asked for.
      // Groups are per (process, param) in practice, so this matters only for a sigma
      // shared across processes, where redrawing for one would quietly move the other.
      vector<int> grp_used(n_grp), grp_draw(n_grp);
      grp_used.setZero();
      grp_draw.fill(1);
      for (int i = 0; i < linkage_re_index.size(); ++i) {
        if (linkage_re_index(i) < 0) continue;
        int grp = linkage_re_sigma(linkage_re_index(i));
        grp_used(grp) = 1;
        // linkage_process also codes composition (5), which simulate_state does not cover
        // because composition linkages carry priors, not random effects. Treat anything
        // outside the five as "not asked for" rather than reading past simulate_state.
        int proc = linkage_process(i);
        if (proc < 0 || proc >= simulate_state.size() || simulate_state(proc) != 1) {
          grp_draw(grp) = 0;
        }
      }

      for (int grp = 0; grp < n_grp; ++grp) {
        if (grp_used(grp) == 0 || grp_draw(grp) == 0) continue;

        // Rogers QAR1: the latent is measured by an observed covariate series, so
        // redrawing it alone would leave the two describing different histories.
        // Left as fitted; sim_mod() warns that it was.
        if (linkage_re_obs(grp) >= 0) continue;

        int len = 0;
        for (int s = 0; s < n_slot; ++s) if (linkage_re_sigma(s) == grp) len++;
        if (len == 0) continue;

        vector<int> slot_of(len);
        int j = 0;
        for (int s = 0; s < n_slot; ++s) if (linkage_re_sigma(s) == grp) slot_of(j++) = s;

        Type sigma = exp(log_sigma_linkage(grp));
        vector<Type> re(len);
        int st = linkage_re_struct(grp);

        if (st == 1) {                                   // rw / RandomWalk
          re(0) = linkage_re_integrate(slot_of(0))
              ? beta_linkage_re(linkage_re_slot(slot_of(0)))
              : beta_linkage_re_pen(linkage_re_slot(slot_of(0)));
          for (int t = 1; t < len; ++t) re(t) = re(t - 1) + rnorm(Type(0), sigma);
        } else if (st == 2) {                            // ar1
          Type rho = rho_trans(trans_rho_linkage(linkage_re_rho(grp)));
          AR1(rho).simulate(re);
          re *= sigma;
        } else {                                         // us / IID
          for (int t = 0; t < len; ++t) re(t) = rnorm(Type(0), sigma);
        }

        for (int t = 0; t < len; ++t) {
          if (linkage_re_integrate(slot_of(t))) {
            beta_linkage_re(linkage_re_slot(slot_of(t))) = re(t);
          } else {
            beta_linkage_re_pen(linkage_re_slot(slot_of(t))) = re(t);
          }
          // A `rw` group keeps slot 0 at its fitted level -- the density never
          // sees it -- so that slot is written but not drawn.
          if (!(st == 1 && t == 0)) beta_linkage_re_drawn_sim(slot_of(t)) = 1;
        }
      }
    }
    REPORT(beta_linkage_re_drawn_sim);
  }

  // Every deviation in one vector, indexed by RE slot: integrated deviations in
  // beta_linkage_re, penalized ones in beta_linkage_re_pen, with linkage_re_slot giving
  // the position within whichever holds it. With no penalized group it copies
  // beta_linkage_re.
  vector<Type> beta_linkage_re_all(linkage_re_sigma.size());
  for (int s = 0; s < beta_linkage_re_all.size(); ++s) {
    beta_linkage_re_all(s) = linkage_re_integrate(s)
        ? beta_linkage_re(linkage_re_slot(s))
        : beta_linkage_re_pen(linkage_re_slot(s));
  }

  // THE _sim NAMING RULE, stated once here and referred to below. TMB never clears the
  // report environment, so anything REPORTed inside a SIMULATE block stays visible in
  // obj$report() for this object's life -- under the observed object's own name, a random
  // replicate would read as the data. Every draw is reported under a name ending _sim,
  // and the composition, CAAL and diet draws go into separate copies (comp_sim, caal_sim,
  // diet_sim) for the same reason.
  SIMULATE {
    vector<Type> beta_linkage_re_sim = beta_linkage_re_all;
    REPORT(beta_linkage_re_sim);
  }

  // Effective linkage coefficient per table row. A fixed row uses its beta_linkage(i); a
  // random-effect row (linkage_re_index >= 0) takes its deviation from the slot-space
  // vector, whose density is jnll_comp row 20. The accumulators below see one beta vector
  // either way; with no RE rows this copies beta_linkage.
  vector<Type> beta_linkage_eff = beta_linkage;
  for (int i = 0; i < beta_linkage_eff.size(); ++i) {
    if (linkage_re_index(i) >= 0) {
      Type z = beta_linkage_re_all(linkage_re_index(i));
      // Rogers QAR1: an observed ar1 latent enters the target scaled by an
      // estimated effect size beta. Unobserved groups (the common case) keep
      // the deviate as-is, so those fits stay bit-identical.
      int grp = linkage_re_sigma(linkage_re_index(i));
      if (linkage_re_obs(grp) >= 0) z *= beta_linkage_obs(linkage_re_obs(grp));
      beta_linkage_eff(i) = z;
    }
  }


  // 5.3. CATCHABILITY
  // Environmental linkage offsets, one pass per link scale. Zero unless a
  // q linkage was supplied, so models without one are unaffected.
  matrix<Type> q_linkage_offset(n_flt, nyrs_hind);     q_linkage_offset.setZero();
  matrix<Type> q_linkage_offset_nat(n_flt, nyrs_hind); q_linkage_offset_nat.setZero();
  matrix<Type> q_linkage_log_mult(n_flt, nyrs_hind);   q_linkage_log_mult.setZero();

  rceattle_apply_q_linkages(
    q_linkage_offset,
    /*link_code=*/ 1,   // log-link rows -> log-scale tensor
    linkage_process, linkage_param, linkage_species, linkage_sex,
    linkage_age_bin, linkage_fleet, linkage_X_col, linkage_link,
    linkage_X, beta_linkage_eff, n_flt, nyrs_hind);

  rceattle_apply_q_linkages(
    q_linkage_offset_nat,
    /*link_code=*/ 0,   // identity-link rows -> natural-scale tensor
    linkage_process, linkage_param, linkage_species, linkage_sex,
    linkage_age_bin, linkage_fleet, linkage_X_col, linkage_link,
    linkage_X, beta_linkage_eff, n_flt, nyrs_hind);

  rceattle_apply_q_linkages(
    q_linkage_log_mult,
    /*link_code=*/ 3,   // exponential rows -> log-scale multiplier tensor
    linkage_process, linkage_param, linkage_species, linkage_sex,
    linkage_age_bin, linkage_fleet, linkage_X_col, linkage_link,
    linkage_X, beta_linkage_eff, n_flt, nyrs_hind);

  REPORT(q_linkage_offset);
  REPORT(q_linkage_offset_nat);
  REPORT(q_linkage_log_mult);

  for(flt = 0; flt < n_flt; flt++){
    for(yr = 0; yr < nyrs_hind; yr++){
      // Order follows SS3: an `exponential` row multiplies log q, and the
      // index_q_dev deviate is added after that product, as SS3 applies its dev
      // blocks after the env link. Inside it the deviate's effective sd would be
      // sd * exp(beta * x), year-varying, against a constant index_q_dev_sd.
      index_q(flt, yr) = exp((index_log_q(flt) + q_linkage_offset(flt, yr))
                               * exp(q_linkage_log_mult(flt, yr))
                             + index_q_dev(flt, yr))
                           + q_linkage_offset_nat(flt, yr);

      // Q as a function of environmental index
      if(est_index_q(flt) == 5){
        beta_q_tmp = index_q_beta.row(flt);
        env_q_tmp = env_index.row(yr) ;
        index_q_mult =  env_q_tmp * beta_q_tmp;
        index_q(flt, yr) = exp(index_log_q(flt) + (index_q_mult).sum());
      }
    }
  }

  matrix<int> r_sexes(diet_obs.rows(), 2); r_sexes.setZero();
  matrix<int> k_sexes(diet_obs.rows(), 2); k_sexes.setZero();


  // 5.4. MATURITY, Finit, AND SEX RATIO
  // -- for SSB derivation, not SPR
  matrix<Type> mature_females = maturity;
  for( sp = 0; sp < nspp ; sp++) {

    if(initMode < 3 || initMode == 5){
      Finit(sp) = 0; // An equilibrium start takes Finit = 0, so R_init and R0 agree. initMode 5 (OffsetEquilibrium) is an F = 0 equilibrium seeded by first-year recruitment.
    }

    // Sex ratio for SSB derivation
    for( age = 0 ; age < nages(sp); age++ ) {
      if(nsex(sp) == 1){
        mature_females( sp, age ) = maturity( sp, age ) * sex_ratio(sp, age); // Multiply sex_ratio and maturity for 1 sex models
      }
    }
  }

  // 5.5. OFFSETS
  // Each process gets two parallel offset tensors, combined at the consume site
  // as param_nat_yr = exp(log_base + log_offset) + nat_offset:
  //   `*_linkage_offset`      log-scale, from log-link rows (linkfn == 1, the
  //                           default); added to log_<param> before exp.
  //   `*_linkage_offset_nat`  natural-scale, from identity-link rows
  //                           (linkfn == 0); added after exp.
  // Catchability has a third, `q_linkage_log_mult`, from exponential-link rows
  // (linkfn == 3), which MULTIPLIES the log rather than shifting it (see 5.3).
  // With no linkages every tensor stays at zero.
  //
  // - RECRUITMENT OFFSETS
  array<Type> recruitment_linkage_offset(nspp,
                                         int(RCEATTLE_N_REC_PARAMS),
                                         nyrs);
  array<Type> recruitment_linkage_offset_nat(nspp,
                                             int(RCEATTLE_N_REC_PARAMS),
                                             nyrs);
  recruitment_linkage_offset.setZero();
  recruitment_linkage_offset_nat.setZero();
  rceattle_apply_recruitment_linkages(
    recruitment_linkage_offset,
    /*link_code=*/ 1,   // log-link rows -> log-scale tensor
    linkage_process,
    linkage_param,
    linkage_species,
    linkage_sex,
    linkage_age_bin,
    linkage_X_col,
    linkage_link,
    linkage_X,
    beta_linkage_eff,
    nspp,
    nyrs
  );
  rceattle_apply_recruitment_linkages(
    recruitment_linkage_offset_nat,
    /*link_code=*/ 0,   // identity-link rows -> natural-scale tensor
    linkage_process,
    linkage_param,
    linkage_species,
    linkage_sex,
    linkage_age_bin,
    linkage_X_col,
    linkage_link,
    linkage_X,
    beta_linkage_eff,
    nspp,
    nyrs
  );
  REPORT(recruitment_linkage_offset);
  REPORT(recruitment_linkage_offset_nat);

  // - M OFFSETS
  array<Type> M_linkage_offset(nspp, max_sex, max_age, nyrs);
  array<Type> M_linkage_offset_nat(nspp, max_sex, max_age, nyrs);
  M_linkage_offset.setZero();
  M_linkage_offset_nat.setZero();
  rceattle_apply_M_linkages(
    M_linkage_offset,
    /*link_code=*/ 1,   // log-link rows -> log-scale tensor
    linkage_process,
    linkage_param,
    linkage_species,
    linkage_sex,
    linkage_age_bin,
    linkage_X_col,
    linkage_link,
    linkage_X,
    beta_linkage_eff,
    nspp,
    nsex,
    nages,
    nyrs
  );
  rceattle_apply_M_linkages(
    M_linkage_offset_nat,
    /*link_code=*/ 0,   // identity-link rows -> natural-scale tensor
    linkage_process,
    linkage_param,
    linkage_species,
    linkage_sex,
    linkage_age_bin,
    linkage_X_col,
    linkage_link,
    linkage_X,
    beta_linkage_eff,
    nspp,
    nsex,
    nages,
    nyrs
  );

  // - GROWTH OFFSETS
  array<Type> growth_linkage_offset(nspp, max_sex, nyrs, RCEATTLE_N_GROWTH_PARAMS);
  array<Type> growth_linkage_offset_nat(nspp, max_sex, nyrs, RCEATTLE_N_GROWTH_PARAMS);
  growth_linkage_offset.setZero();
  growth_linkage_offset_nat.setZero();
  rceattle_apply_growth_linkages(
    growth_linkage_offset,
    /*link_code=*/ 1,   // log-link rows -> log-scale tensor
    linkage_process,
    linkage_param,
    linkage_species,
    linkage_sex,
    linkage_age_bin,
    linkage_X_col,
    linkage_link,
    linkage_X,
    beta_linkage_eff,
    nspp,
    nsex,
    nyrs
  );
  rceattle_apply_growth_linkages(
    growth_linkage_offset_nat,
    /*link_code=*/ 0,   // identity-link rows -> natural-scale tensor
    linkage_process,
    linkage_param,
    linkage_species,
    linkage_sex,
    linkage_age_bin,
    linkage_X_col,
    linkage_link,
    linkage_X,
    beta_linkage_eff,
    nspp,
    nsex,
    nyrs
  );

  // 5.5b. DSEM ON THE RECRUITMENT DEVIATIONS
  // Where a DSEM is supplied, the recruitment deviations are the latent states
  // of a GMRF rather than free parameters: calculate_dsem() supplies the
  // density, and rec_dev is OVERWRITTEN from those states. Everything
  // downstream -- the recruitment build below, initial abundance, projections --
  // reads rec_dev and R_sd exactly as before, so nothing else needs to know.
  //
  // Placed here, before rec_dev and R_sd are first read, and guarded on
  // dsem_on: with no DSEM this block does not execute and the standard
  // parameterization is untouched.
  Type jnll_dsem = 0;
  array<Type> dsem_z_tj(dsem_y_tj.rows(), dsem_y_tj.cols()); dsem_z_tj.setZero();
  // Marginal variance of each latent state, for the lognormal bias correction
  // below. Only computed when the correction is actually on -- it costs an LU
  // solve with n_k right-hand sides per objective evaluation.
  array<Type> dsem_margvar_tj(dsem_y_tj.rows(), dsem_y_tj.cols()); dsem_margvar_tj.setZero();
  int dsem_want_margvar = (bias_adjust_proc > 0) ? 1 : 0;
  // Precision and mean of the latent field, handed back so R can draw from the
  // SAME density that scores it -- the mechanism dsem's own simulate() uses
  // (rmvnorm_prec(xhat + delta, Q)). Without these a caller has to reproduce
  // the recursion, which is a second implementation of the process and free to
  // drift from the first.
  matrix<Type> dsem_Q(0, 0);
  array<Type> dsem_xhat_tj(dsem_y_tj.rows(), dsem_y_tj.cols()); dsem_xhat_tj.setZero();
  array<Type> dsem_delta_tj(dsem_y_tj.rows(), dsem_y_tj.cols()); dsem_delta_tj.setZero();
  // Mean and SD of the covariate MEASUREMENT density, for the observation draw
  // below. A covariate is data only under family = "fixed"; every other family
  // makes it a latent state observed with error, and that observation is drawn
  // like any other.
  array<Type> dsem_mu_tj(dsem_y_tj.rows(), dsem_y_tj.cols()); dsem_mu_tj.setZero();
  vector<Type> dsem_sigma_z(dsem_lnsigma_z.size()); dsem_sigma_z.setZero();
  // The covariate observations as simulated, and a mask of the cells the draw
  // touched. Declared out here so both are REPORTed whether or not this model
  // carries a DSEM; with dsem_on = 0 they are the (empty) inputs, unchanged.
  array<Type> dsem_y_tj_sim = dsem_y_tj;
  array<Type> dsem_y_tj_drawn_sim(dsem_y_tj.rows(), dsem_y_tj.cols());
  dsem_y_tj_drawn_sim.setZero();
  // The DSEM process draw is taken in R (sim_mod()), not here, and this flag
  // stays 0. calculate_dsem()'s own do_simulate branch assigns the WHOLE latent
  // field -- `x_tj = draw_tj` in dsem.hpp -- so it redraws the covariate columns
  // as well as the recruitment ones and consults neither the map nor dsem_cond_k.
  // Under family = "fixed" those covariate columns ARE the env_data series, so
  // it generated each replicate under an environment the refit is never shown
  // (measured on a scaled series with sd 1: the drawn column moved by up to
  // 2.96). Drawing in R instead draws the recruitment columns CONDITIONAL on
  // everything the model was given, from this same reported precision, which is
  // the density that scored the fit -- so there is still exactly one
  // implementation of the process.
  int dsem_do_sim = 0;
  // Which latent cells the model is GIVEN, so the bias correction below uses
  // the variance of recruitment conditional on them. A cell is given when
  // env_data supplied a finite value for that variable in that year.
  //
  // Per CELL: build_dsem_objects() pads env_data to the model's full span, so a
  // covariate that starts late, has a gap, or carries no projection scenario is
  // a free latent state in those years. Treating the whole column as known
  // under-corrects there -- 17.1% high on projected recruitment measured on
  // BS2017SS with a hindcast-only covariate, and 48.8% once the covariate is
  // autoregressive.
  //
  // Two exclusions:
  //
  //   Recruitment-deviation columns, outright. They are latent by construction,
  //   so a stray value in one cannot zero that species' correction.
  //
  //   dsem_unobs_idx, the cells the projecting parameterizations compute rather
  //   than estimate. calculate_dsem() overwrites z_tj at those cells from the
  //   deterministic projection and never reads y there -- scaling such a column
  //   leaves the objective bit-identical -- so they are not values the model is
  //   given. They also carry no exogenous variance, which would make the
  //   conditioning solve singular.
  //
  // A variable with a measurement likelihood (family other than "fixed") is
  // conditioned on its OBSERVATION rather than on a pinned state, which
  // understates the correction by the measurement variance: exact at a small
  // observation SD (0.2% at SD = 0.05) and up to 8% on recruitment at SD = 1.
  // Treating it as unknown instead would be the joint prior variance, which is
  // wrong by more and in the same direction as the defect above.
  // The two switches have to be a coherent pair, and an incoherent one is
  // silent rather than loud: dsem_drives_rec = 1 with no SEM leaves rec_dev
  // free, in `random`, and scored by nothing -- a flat improper prior that
  // still converges -- and reads rec_sd_use_prior, which is zero-length
  // without a SEM, out of bounds. Anything other than 0/1 lands in the same
  // state, because the R side tests == 1 while the template tests == 0. The
  // builders only ever emit (0,0) or (1,1); a caller assembling a data list by
  // hand, which the regression test does deliberately, is what reaches this.
  if(dsem_drives_rec != 0 && dsem_drives_rec != 1)
    Rf_error("dsem_drives_rec is %d; it must be 0 or 1", (int)dsem_drives_rec);
  if(dsem_drives_rec == 1 && dsem_on == 0)
    Rf_error("dsem_drives_rec = 1 with dsem_on = 0: no SEM is attached, so nothing can drive the recruitment deviations");

  vector<int> dsem_cond_k(dsem_y_tj.rows() * dsem_y_tj.cols());
  dsem_cond_k.setZero();
  // dsem_on, NOT dsem_drives_rec: the conditioning count describes the latent
  // field's observed cells, which exist whatever process the SEM drives.
  if(dsem_on == 1){
    int n_t_dsem = dsem_y_tj.rows();
    for(int j2 = 0; j2 < dsem_y_tj.cols(); j2++){
      bool is_rec = false;
      for(int sp2 = 0; sp2 < rec_dev_col.size(); sp2++)
        if(rec_dev_col(sp2) == j2) is_rec = true;
      if(is_rec) continue;
      for(int t2 = 0; t2 < n_t_dsem; t2++)
        if(R_FINITE(asDouble(dsem_y_tj(t2, j2))))
          dsem_cond_k(j2 * n_t_dsem + t2) = 1;
    }
    for(int u = 0; u < dsem_unobs_idx.size(); u++){
      if(dsem_unobs_idx(u) >= 0 && dsem_unobs_idx(u) < dsem_cond_k.size())
        dsem_cond_k(dsem_unobs_idx(u)) = 0;
    }
  }
  // dsem_on, NOT dsem_drives_rec: this builds the latent field and its
  // density, which every SEM needs whatever process it drives.
  if(dsem_on == 1){
    calculate_dsem(jnll_dsem, dsem_options, dsem_RAM, dsem_RAMstart,
                   dsem_familycode_j, dsem_linkcode_j, dsem_sigmastart_j,
                   dsem_eps_tj, dsem_y_tj, dsem_obs_idx, dsem_unobs_idx,
                   dsem_beta_z, dsem_lnsigma_z, dsem_mu_j, dsem_delta0_j,
                   dsem_x_tj, dsem_z_tj, dsem_Q, dsem_xhat_tj, dsem_delta_tj,
                   dsem_do_sim, dsem_want_margvar, dsem_cond_k,
                   dsem_margvar_tj, dsem_mu_tj, dsem_sigma_z);

    // Dimension contract for THIS call site. calculate_dsem() checks its own
    // inputs against each other; nyrs_dsem and rec_dev_col are ours.
    if(nyrs_dsem > dsem_z_tj.rows())
      Rf_error("DSEM spans %d years but the latent states have %d rows",
               (int)nyrs_dsem, (int)dsem_z_tj.rows());

    // Against the MODEL's year span, not just the DSEM's own. Every tensor a
    // latent is written into is dimensioned nyrs, and the rec_dev.cols() check
    // below is inside dsem_drives_rec == 1 -- so a growth-only SEM reached the
    // growth writer with no year bound at all. A dsem built for a longer model
    // and handed to fit_mod() through `dsem = <built objects>` is the reached
    // case.
    if(nyrs_dsem > nyrs)
      Rf_error("DSEM spans %d years but the model spans %d -- the DSEM was built for a different model",
               (int)nyrs_dsem, (int)nyrs);

    // The recruitment takeover: R_sd from the SEM's variance path and rec_dev
    // from the latent states. Gated on dsem_drives_rec, not on dsem_on -- a SEM
    // driving some other process must leave both alone, or rec_dev would be
    // scored by JNLL_REC_DEV and then overwritten from the field (a density
    // evaluated away from its own parameter) while R_log_sd sat estimated with
    // zero gradient, a flat direction and a singular Hessian.
    if(dsem_drives_rec == 1){
    if(nyrs_dsem > rec_dev.cols())
      Rf_error("DSEM spans %d years but rec_dev has %d columns -- the DSEM was built for a different model",
               (int)nyrs_dsem, (int)rec_dev.cols());
    if(rec_dev_col.size() != nspp)
      Rf_error("rec_dev_col has %d entries, expected %d (one per species)",
               (int)rec_dev_col.size(), nspp);
    // The value, not only the count. Where the SEM does not drive recruitment
    // the R side sets -1, a deliberate non-column; reading dsem_z_tj at it would
    // be out of bounds. Reached by hand-setting dsem_drives_rec = 1 on a SEM
    // that does not name a recruitment column, which the test suite does.
    for(sp = 0; sp < nspp; sp++)
      if(rec_dev_col(sp) < 0 || rec_dev_col(sp) >= dsem_z_tj.cols())
        Rf_error("dsem_drives_rec = 1 but rec_dev_col[%d] is %d; this SEM does not name a recruitment latent column",
                 (int)sp, (int)rec_dev_col(sp));

    for(sp = 0; sp < nspp; sp++){
      // R_sd FIRST -- the bias correction below uses it. The recruitment SD is
      // the SEM's two-headed self-loop on this species' recdevs; its sign is not
      // identified (it is a Cholesky factor), hence the absolute value.
      if(rec_sd_idx(sp) >= 1){
        R_sd(sp) = sqrt(square( dsem_beta_z(rec_sd_idx(sp) - 1) ));
      } else {
        R_sd(sp) = rec_sd_fixed(sp);
      }

      // The latent states span nyrs_dsem, which is the hindcast only unless
      // build_DSEM(estimate_projection = TRUE). Copy just that span -- a
      // whole-row assignment would be a dimension error, or worse would
      // misalign if the lengths happened to match.
      //
      // Lognormal bias correction. The GMRF centres its latent states at 0,
      // while the standard density centres the deviations at -sigma^2/2 so that
      // E[R] = R0. Without this a DSEM has a different shrinkage target from the
      // equivalent non-DSEM model: SSB absorbs the offset, but R0 does not --
      // measured 20-51% low on BS2017SS with an IID sem -- and dynamic B0 and
      // the B40% proxy are keyed to R0.
      //
      // The correction is -Var/2 with Var the MARGINAL variance of this
      // species' recruitment-deviation column, which is what makes E[exp(dev)]
      // = 1. It is NOT -R_sd^2/2: R_sd is the GMRF's conditional (innovation)
      // SD, and the two coincide only for a state with no incoming lagged
      // paths. For a first-order self-path they differ by 1/(1-rho^2), so using
      // R_sd would under-correct exactly the lagged sems a DSEM exists to fit.
      // calculate_dsem() returns the marginal variance per state, taken from
      // diag((I-Rho)^-1 Gamma^T Gamma (I-Rho)^-T) after any constant_variance
      // rescaling, so it is correct term by term and year by year rather than
      // relying on a stationary approximation.
      for(yr = 0; yr < nyrs_dsem; yr++){
        rec_dev(sp, yr) = dsem_z_tj(yr, rec_dev_col(sp))
                          - bias_adjust_proc * dsem_margvar_tj(yr, rec_dev_col(sp)) / 2.0;
      }
    }
    } // end dsem_drives_rec == 1

    // 5.5b-growth. GROWTH LATENTS INTO THE GROWTH OFFSET
    // A SEM may drive mean growth per species and sex: L1 (cm at the anchor
    // age), Linf (cm) or the Richards shape m. The latent is a log-scale
    // multiplier on the parameter, which is what growth_linkage_offset already
    // is, so it lands in the same tensor section 5.7 reads at
    // exp(log_growth_pars + growth_linkage_offset).
    //
    // `+=`, not `=`: a fixed environmental linkage may already have written this
    // cell, and the two compose -- the covariate gives the explained part and
    // the latent the residual. Overwriting would silently delete the linkage.
    //
    // Bias-corrected on the same terms as the recruitment deviations above and
    // for the same reason: the offset is exponentiated, so without
    // -margvar/2 the mean parameter drifts up with the latent's variance and
    // mean weight-at-age -- the weight in SSB, B40% and the ABC -- goes with it.
    // The correction uses the GMRF's MARGINAL variance per state, not the
    // innovation SD.
    //
    // nyrs_dsem, not nyrs: with estimate_projection = FALSE -- the DEFAULT --
    // the latent matrix is shorter than the model's year span, and the
    // projection years keep whatever the linkage fill left them (zero unless a
    // covariate linkage wrote there). This applies to section 5.5b-M below too.
    // The four maps are PARALLEL -- entry g of each describes one latent -- but
    // the R side supplies them as four independent fields, so a hand-built or
    // hand-edited $dsem list can carry them at different lengths. Checked
    // before the loop because the value guards inside it cannot see this: a
    // SHORT sp/sex/par vector is read out of bounds (measured: a segfault, or
    // two latents silently landing on one species, halving one species' offset
    // and doubling another's), and a LONG col vector writes phantom entries.
    if(growth_dev_sp.size() != growth_dev_col.size() ||
       growth_dev_sex.size() != growth_dev_col.size() ||
       growth_dev_par.size() != growth_dev_col.size())
      Rf_error("growth_dev_col/sp/sex/par have %d/%d/%d/%d entries; they are parallel and must agree",
               (int)growth_dev_col.size(), (int)growth_dev_sp.size(),
               (int)growth_dev_sex.size(), (int)growth_dev_par.size());

    for(int g = 0; g < growth_dev_col.size(); g++){
      // Dimension contract for this call site, as for rec_dev_col above: the R
      // side builds these from the latent registry, but a hand-built or
      // hand-edited $dsem list reaches here too, and an out-of-range index
      // would write into a neighbouring species' growth.
      if(growth_dev_col(g) < 0 || growth_dev_col(g) >= dsem_z_tj.cols())
        Rf_error("growth_dev_col[%d] is %d; it must index a column of dsem_z_tj (0-%d)",
                 g, (int)growth_dev_col(g), (int)dsem_z_tj.cols() - 1);
      if(growth_dev_sp(g) < 0 || growth_dev_sp(g) >= nspp)
        Rf_error("growth_dev_sp[%d] is %d; it must be a 0-based species (0-%d)",
                 g, (int)growth_dev_sp(g), (int)nspp - 1);
      if(growth_dev_sex(g) < 0 || growth_dev_sex(g) >= nsex(growth_dev_sp(g)))
        Rf_error("growth_dev_sex[%d] is %d, but species %d has %d sex(es)",
                 g, (int)growth_dev_sex(g), (int)growth_dev_sp(g) + 1,
                 (int)nsex(growth_dev_sp(g)));
      if(growth_dev_par(g) < 0 || growth_dev_par(g) >= RCEATTLE_N_GROWTH_PARAMS)
        Rf_error("growth_dev_par[%d] is %d; it must be a mean-growth slot (0-%d)",
                 g, (int)growth_dev_par(g), RCEATTLE_N_GROWTH_PARAMS - 1);

      for(yr = 0; yr < nyrs_dsem; yr++){
        growth_linkage_offset(growth_dev_sp(g), growth_dev_sex(g), yr,
                              growth_dev_par(g)) +=
          dsem_z_tj(yr, growth_dev_col(g))
          - bias_adjust_proc * dsem_margvar_tj(yr, growth_dev_col(g)) / 2.0;
      }
    }

    // 5.5b-M. WRITE THE NATURAL-MORTALITY LATENTS INTO M_linkage_offset
    // The latent is a log-scale multiplier on residual natural mortality M1,
    // which is what M_linkage_offset already is: section 6.1 reads
    // exp(log_M1 + log_M1_dev + M1_mult.sum() + M_linkage_offset). One latent
    // per species and sex, broadcast FLAT over ages -- a year in which
    // mortality was high was high for the whole stock. `M1_re` writes the same
    // exponent through log_M1_dev, and its modes 2/5 have exactly this shape,
    // so a sem naming an M latent for a species whose M1_re is on is refused on
    // the R side.
    //
    // `+=` and the -margvar/2 bias correction for the same reasons as growth
    // above: a fixed environmental linkage may have written this cell and the
    // two compose, and the offset is exponentiated, so without the correction
    // mean M1 drifts up with the latent's variance. The correction leaves mean
    // M1 unbiased, not mean survivorship: exp(-M) is convex, so SB0 and B40%
    // keep a Jensen bias, as mean weight-at-age does under growth.
    // Parallel maps, same contract and same hazard as growth above.
    if(mdev_sp.size() != mdev_col.size() ||
       mdev_sex.size() != mdev_col.size())
      Rf_error("mdev_col/sp/sex have %d/%d/%d entries; they are parallel and must agree",
               (int)mdev_col.size(), (int)mdev_sp.size(), (int)mdev_sex.size());

    for(int m = 0; m < mdev_col.size(); m++){
      if(mdev_col(m) < 0 || mdev_col(m) >= dsem_z_tj.cols())
        Rf_error("mdev_col[%d] is %d; it must index a column of dsem_z_tj (0-%d)",
                 m, (int)mdev_col(m), (int)dsem_z_tj.cols() - 1);
      if(mdev_sp(m) < 0 || mdev_sp(m) >= nspp)
        Rf_error("mdev_sp[%d] is %d; it must be a 0-based species (0-%d)",
                 m, (int)mdev_sp(m), (int)nspp - 1);
      if(mdev_sex(m) < 0 || mdev_sex(m) >= nsex(mdev_sp(m)))
        Rf_error("mdev_sex[%d] is %d, but species %d has %d sex(es)",
                 m, (int)mdev_sex(m), (int)mdev_sp(m) + 1,
                 (int)nsex(mdev_sp(m)));

      for(yr = 0; yr < nyrs_dsem; yr++){
        Type mdev_yr = dsem_z_tj(yr, mdev_col(m))
          - bias_adjust_proc * dsem_margvar_tj(yr, mdev_col(m)) / 2.0;
        // Ages 0 .. nages(sp)-1, the span section 6.1 reads. The tensor is
        // sized to the longest-lived species, so a shorter-lived one's trailing
        // cells are its own padding and stay zero.
        for(age = 0; age < nages(mdev_sp(m)); age++){
          M_linkage_offset(mdev_sp(m), mdev_sex(m), age, yr) += mdev_yr;
        }
      }
    }

    // 5.5c. SIMULATE THE COVARIATE OBSERVATIONS
    // Under family = "fixed" a covariate column IS the environmental series:
    // familycode 0, no measurement density, pinned in the map, nothing to draw.
    // Every other family makes the column a latent state OBSERVED WITH ERROR,
    // and that observation is data like a survey index -- so it is redrawn from
    // the density in dsem.hpp that scores it, using the mean and SD that density
    // used. Without this a replicate carried the ORIGINAL covariate series
    // beside freshly drawn recruitment, so the covariate looked more informative
    // than the model says it is and a self test overstated how well the SEM's
    // paths are recovered.
    //
    // NOT gated on simulate_state: this is an observation, and sim_mod() draws
    // every observation whenever simulate = TRUE. Process error decides whether
    // the LATENT state moves, not whether its observation is redrawn.
    //
    // A cell with no observation stays unobserved -- drawing there would hand
    // the refit data the assessment never had.
    SIMULATE {
      for(int t2 = 0; t2 < dsem_y_tj.rows(); t2++){
        for(int j2 = 0; j2 < dsem_y_tj.cols(); j2++){
          int fam = dsem_familycode_j(j2);
          if(fam == 0) continue;
          if(!R_FINITE(asDouble(dsem_y_tj(t2, j2)))) continue;
          // dsem.hpp fills mu_tj through a four-way link switch with no else,
          // so a link code outside 0-3 leaves the cell UNINITIALIZED. Scoring
          // reads it too, but a draw would turn that into data. Refuse rather
          // than draw from whatever was on the stack.
          if(dsem_linkcode_j(j2) < 0 || dsem_linkcode_j(j2) > 3) continue;
          // Tweedie reads a second sigma slot for its power parameter, so both
          // have to exist before either is read. A family whose SDs are not
          // where it expects them is a build error, not something to draw
          // through: leave the observation alone and let the mask say so.
          int s0 = dsem_sigmastart_j(j2);
          int need = (fam == 7) ? s0 + 1 : s0;
          if(s0 < 0 || need >= dsem_sigma_z.size()) continue;
          Type mu = dsem_mu_tj(t2, j2);
          Type sd = dsem_sigma_z(s0);
          if(fam == 1){        // normal
            dsem_y_tj_sim(t2, j2) = rnorm(mu, sd);
          } else if(fam == 2){ // Bernoulli
            dsem_y_tj_sim(t2, j2) = rbinom(Type(1.0), mu);
          } else if(fam == 3){ // Poisson
            dsem_y_tj_sim(t2, j2) = rpois(mu);
          } else if(fam == 4){ // Gamma; shape = 1/CV^2, scale = mean*CV^2
            dsem_y_tj_sim(t2, j2) = rgamma(pow(sd, Type(-2.0)), mu * pow(sd, Type(2.0)));
          } else if(fam == 5){ // normal with known SD
            dsem_y_tj_sim(t2, j2) = rnorm(mu, dsem_eps_tj(t2, j2));
          } else if(fam == 6){ // lognormal
            dsem_y_tj_sim(t2, j2) = exp(rnorm(log(mu), sd));
          } else if(fam == 7){ // Tweedie; phi and p exactly as the density reads them
            dsem_y_tj_sim(t2, j2) = rtweedie(mu, exp(dsem_sigma_z(s0)),
                                             Type(1.0) + invlogit(dsem_sigma_z(s0 + 1)));
          } else {
            continue;
          }
          dsem_y_tj_drawn_sim(t2, j2) = 1;
        }
      }
    }
  }

  // Position matters: REPORT() copies the tensor at the call site, and sections
  // 5.5b-growth and 5.5b-M add a DSEM's latent states to growth_linkage_offset
  // and M_linkage_offset above. Reported before that, either offset would be a
  // pre-DSEM snapshot contradicting the growth parameters and M-at-age that
  // sections 5.7 and 6.1 actually used.
  REPORT(M_linkage_offset);
  REPORT(M_linkage_offset_nat);
  REPORT(growth_linkage_offset);
  REPORT(growth_linkage_offset_nat);

  // 5.6. RECRUITMENT PARAMETERS
  // Linkage offsets combine log-link (multiplicative) and identity-link
  // (natural-scale additive) contributions:
  //   R0(yr) = exp(rec_pars(sp,0) + log_offset) + nat_offset.
  // TODO: a negative nat_offset can make R0, alpha or Beta non-positive; only hindcast R, R_hat and the penalty curve are floored (issue #234).
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs; yr++){
      R0(sp, yr)    = exp(rec_pars(sp, 0) + recruitment_linkage_offset(sp, RCEATTLE_REC_R0,    yr))
                    + recruitment_linkage_offset_nat(sp, RCEATTLE_REC_R0,    yr);
      alpha(sp, yr) = exp(rec_pars(sp, 1) + recruitment_linkage_offset(sp, RCEATTLE_REC_ALPHA, yr))
                    + recruitment_linkage_offset_nat(sp, RCEATTLE_REC_ALPHA, yr);
      Beta(sp, yr)  = exp(rec_pars(sp, 2) + recruitment_linkage_offset(sp, RCEATTLE_REC_BETA,  yr))
                    + recruitment_linkage_offset_nat(sp, RCEATTLE_REC_BETA,  yr);
    }
  }


  // 5.7. GROWTH
  // -- Rearange growth parameters. Log-link offset enters additively
  //    on the log scale (multiplicative on natural K/L1/Linf/m);
  //    identity-link offset enters additively on the natural scale.
  array<Type> growth_parameters(nspp, max_sex, nyrs, RCEATTLE_N_GROWTH_PARAMS); growth_parameters.setZero(); // K, L1, Linf, m
  for(sp = 0; sp < nspp; sp++){
    for(sex = 0; sex < nsex(sp); sex ++){
      for(yr = 0; yr < nyrs; yr++){
        for(int par = 0; par < 4; par++){
          growth_parameters(sp, sex, yr, par) = exp(
            log_growth_pars(sp, sex, par)
          + growth_linkage_offset(sp, sex, yr, par)
          ) + growth_linkage_offset_nat(sp, sex, yr, par);
        }
      }
    }
  }
  REPORT(growth_parameters);

  // -- Calculate weight
  calculate_weight(
    weight_hat,
    length_hat,
    growth_matrix,
    mat_weight_hat,
    weight_obs,
    growth_model,
    growth_sd_style,
    growth_sd_form,
    growth_plus_length,
    plus_group_decay,
    nspp,
    nyrs,
    nyrs_hind,
    n_flt,
    flt_spp,
    flt_month,
    nsex,
    minage,
    growth_age_L1,
    nages,
    nlengths,
    pop_wt_index,
    ssb_wt_index,
    flt_wt_index,
    spawn_month,
    lengths,
    nlengths_pop,
    lengths_pop,
    growth_parameters,
    growth_log_sd,
    weight_length_pars,
    mat_len_use,
    mat_len_pars,
    log_M1
  );

  // -- Spawning output per fish (kg) at spawning time, female-only for a
  //    two-sex species. Weight x maturity-at-age x sex ratio by default; with
  //    maturity-at-length it is the mature weight integrated over the length
  //    distribution at spawning, so maturity and weight vary together with length.
  array<Type> spawn_output(nspp, max_age, nyrs); spawn_output.setZero();
  for(sp = 0; sp < nspp; sp++){
    wt_idx_ssb = 2 * sp + 1;
    for(age = 0; age < nages(sp); age++){
      for(yr = 0; yr < nyrs; yr++){
        if(mat_len_use(sp) == 1){
          spawn_output(sp, age, yr) = mat_weight_hat(wt_idx_ssb, 0, age, yr) * ((nsex(sp) == 1) ? sex_ratio(sp, age) : Type(1.0));
        } else {
          spawn_output(sp, age, yr) = weight_hat(wt_idx_ssb, 0, age, yr) * mature_females(sp, age);
        }
      }
    }
  }
  REPORT(spawn_output);


  // 5.8. SELECTIVITY
  // Environmental-linkage offsets. Two tensors per parameter family (log and
  // natural scale), all zero unless a selectivity linkage was supplied.
  array<Type> sel_slp_off (2, n_flt, max_sex, nyrs);            sel_slp_off.setZero();
  array<Type> sel_slp_off_nat (2, n_flt, max_sex, nyrs);        sel_slp_off_nat.setZero();
  array<Type> sel_inf_off (2, n_flt, max_sex, nyrs);            sel_inf_off.setZero();
  array<Type> sel_inf_off_nat (2, n_flt, max_sex, nyrs);        sel_inf_off_nat.setZero();
  array<Type> sel_coff_off (n_flt, max_sex, max_bin, nyrs);     sel_coff_off.setZero();
  array<Type> sel_coff_off_nat (n_flt, max_sex, max_bin, nyrs); sel_coff_off_nat.setZero();
  array<Type> sel_apical_off (n_flt, max_sex, nyrs);            sel_apical_off.setZero();
  array<Type> sel_apical_off_nat (n_flt, max_sex, nyrs);        sel_apical_off_nat.setZero();
  array<Type> sel_dn6_off        (6, n_flt, max_sex, nyrs);     sel_dn6_off.setZero();
  array<Type> sel_dn6_off_nat    (6, n_flt, max_sex, nyrs);     sel_dn6_off_nat.setZero();

  // The per-sex apical multiplier is compiled in only when a selectivity
  // linkage names it (param 5); otherwise the curve is untouched and the AD
  // tape of every existing model is unchanged.
  int sel_apical_on = 0;
  for (int i = 0; i < linkage_process.size(); i++) {
    if (linkage_process(i) == RCEATTLE_PROC_SEL && linkage_param(i) == 5) sel_apical_on = 1;
  }

  rceattle_apply_sel_linkages(
    sel_slp_off, sel_inf_off, sel_coff_off, sel_apical_off, sel_dn6_off,
    /*link_code=*/ 1,   // log-link rows -> log-scale tensors
    linkage_process, linkage_param, linkage_species, linkage_sex,
    linkage_age_bin, linkage_fleet, linkage_X_col, linkage_link,
    linkage_X, beta_linkage_eff, n_flt, max_sex, max_bin, nyrs);

  rceattle_apply_sel_linkages(
    sel_slp_off_nat, sel_inf_off_nat, sel_coff_off_nat, sel_apical_off_nat, sel_dn6_off_nat,
    /*link_code=*/ 0,   // identity-link rows -> natural-scale tensors
    linkage_process, linkage_param, linkage_species, linkage_sex,
    linkage_age_bin, linkage_fleet, linkage_X_col, linkage_link,
    linkage_X, beta_linkage_eff, n_flt, max_sex, max_bin, nyrs);

  // "sel_at_age" and "sel_at_length" modified via pass-by-reference
  // when "sel_at_length" is used, it is converted to "sel_at_age" using the growth matrix
  calculate_selectivity(
    nspp,                 // Number of species.
    n_flt,                // Number of fleets.
    nyrs,                 // Years, hindcast plus projection.
    nyrs_hind,            // Years in the hindcast.
    styr,                 // Start year
    nsex,                 // Sexes modelled per species: 1 = combined, 2 = both.
    nages,                // Age bins per species.
    nlengths,             // Observation length bins per species.
    lengths,              // Length bin boundaries matrix
    nlengths_pop,         // Population length bins per species (= nlengths when no pop grid)
    lengths_pop,          // Population length bin lower edges
    flt_spp,              // Fleet to species mapping
    flt_sel_type,         // Selectivity model type per fleet
    flt_varying_sel,      // Time_varying_sel per fleet (picks NonParametricIntegrable's construction)
    flt_sel_dim,          // Age or length based
    bin_first_selected,   // Min bin selected per fleet
    flt_n_sel_bins,       // Max estimated bins per fleet
    flt_sel_cap_bin,      // Bin (0-based) at/after which realized non-par sel is capped flat (NonParametricPM)
    sel_norm_bin1,        // Normalization control/bin 1
    sel_norm_bin2,        // Normalization control/bin 2
    sel_norm_scope,       // Normalization reference pooled across sexes?
    flt_sel_start_yr,     // Per-fleet selectivity start year (0-based)
    emp_sel_obs,          // Empirical observations matrix
    emp_sel_ctl,          // Empirical control matrix
    log_sel_slp,           // Logistic slope parameters
    log_sel_slp_dev,       // Slope deviations
    sel_inf,              // Inflection parameters
    sel_inf_dev,          // Inflection deviations
    sel_coff,             // Non-parametric coefficients
    sel_coff_dev,         // Coefficient deviations
    avg_sel,              // [Modified] Average selectivity
    non_par_sel,          // [Modified] Unnormalized non-parametric selectivity
    log_non_par_sel,      // [Modified] The same curve on the log scale, for the penalties
    sel_at_length,        // [Modified] Final length-based selectivity
    sel_at_age,           // [Modified] Final age-based selectivity
    growth_matrix,        // Length to age transition matrix
    sel_slp_off, sel_slp_off_nat,   // selectivity linkage offsets (log / natural)
    sel_inf_off, sel_inf_off_nat,
    sel_coff_off, sel_coff_off_nat,
    sel_apical_on, log_sel_apical,  // per-sex apical height, only when a linkage names it
    sel_apical_off, sel_apical_off_nat,
    sel_dn6, sel_dn6_off, sel_dn6_off_nat, sel_dn6_ends   // DoubleNormalSS3
  );

  // -- Selected body weight (kg): sum_l P(l|a) s(l) w(l) / sum_l P(l|a) s(l),
  //    SS3's "bodywt". A length-selective fleet takes the larger or smaller fish
  //    of an age class, so its catch and survey biomass weigh that class by the
  //    fish it selects, not by the class mean. Estimated growth only.
  for(flt = 0; flt < n_flt; flt++){
    sp = flt_spp(flt);
    if(growth_model(sp) == 0 || flt_sel_dim(flt) != 1 || flt_sel_type(flt) == 0) continue;
    int wt_idx_sel = nspp * 2 + flt;
    for(sex = 0; sex < nsex(sp); sex++){
      for(age = 0; age < nages(sp); age++){
        for(yr = 0; yr < nyrs; yr++){
          Type num = 0, den = 0;
          // Integrated over POPULATION bins: both the key and the curve live
          // there, so the selectivity weighting is applied at the resolution
          // the curve is defined on rather than to a bin average.
          for(int ln = 0; ln < nlengths_pop(sp); ln++){
            // The last bin takes the previous bin's width; a single-bin grid has
            // no neighbour to take one from, so it keeps the edge itself.
            Type lenmid = (ln < nlengths_pop(sp) - 1)
              ? (lengths_pop(sp, ln) + lengths_pop(sp, ln + 1)) / Type(2.0)
              : (nlengths_pop(sp) < 2
                   ? lengths_pop(sp, ln)
                   : lengths_pop(sp, ln) + (lengths_pop(sp, ln) - lengths_pop(sp, ln - 1)) / Type(2.0));
            Type p_sel = growth_matrix(wt_idx_sel, sex, age, ln, yr) * sel_at_length(flt, sex, ln, yr);
            num += p_sel * weight_length_pars(sp, 0) * pow(lenmid, weight_length_pars(sp, 1));
            den += p_sel;
          }
          // An age class the fleet does not select contributes nothing either way.
          weight_hat(wt_idx_sel, sex, age, yr) = num / (den + Type(1e-30));
        }
      }
    }
  }


  // 5.9. BIOENERGETICS AND CONSUMPTION
  // - Calculate temperature function of consumption
  calculate_temperature_function(fT, nspp, nyrs, Ceq, Cindex, Qc, Tcm, Tco, Tcl, CK1, CK4, env_index);

  // - Calculate historic ration
  calculate_ration(consumption_at_age, nspp, nyrs, nyrs_hind, nsex, nages, Ceq, CA, CB, fday, Pvalue, fT, ration_data, weight_hat);


  // 5.10. PARAMETRIC SUITABILITY
  // - Estimate vulnerability (Parametric Suitability)
  calculate_vulnerability(vulnerability, vulnerability_other, nspp, suitMode, log_phi);

  // - Parametric Suitability (Gamma, Lognormal, Normal)
  calculate_parametric_suitability(suitability, suit_other, nspp, nyrs, nsex, nages, suitMode,
                                   length_hat, weight_hat, vulnerability, vulnerability_other, gam_a, gam_b);


  // 5.11. REORGANIZE MSVPA DIET DATA
  // - Reorganize diet_obs content
  organize_diet_obs(nspp, nyrs, nyrs_hind, styr, minage, nsex,
                    diet_obs, diet_ctl, diet_prop);

  // - Calculate other food stomach content
  calculate_other_food_diet_prop(nspp, nyrs, nsex, nages, other_food,
                                 diet_prop, other_food_diet_prop);

  // 5.12. FISHING MORTALITY and FSPRs
  // Outside the multi-species iteration loop: hindcast F needs no M. The
  // forecast arm below is provisional -- section 6.7 recomputes it once the
  // HCR has its reference points.
  F_spp.setZero();
  F_flt_age.setZero();
  F_at_age.setZero();
  Ftarget_at_age.setZero();
  Flimit_at_age.setZero();
  for(flt = 0; flt < n_flt; flt++) {

    sp = flt_spp(flt);  // Species this fleet fishes or surveys

    if(flt_type(flt) == 1){
      for(age = 0; age < nages(sp); age++) {
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 0; yr < nyrs; yr++) {

            // Hindcast
            if( yr < nyrs_hind){
              F_flt_age(flt, sex, age, yr) = sel_at_age(flt, sex, age, yr) * exp(log_F(flt, yr));
            }


            // Forecast
            if( yr >= nyrs_hind){
              // -- Apply HCRs
              switch(HCR){
              case 0: // No fishing
                proj_F(sp, yr) = 0;
                break;

              case 1: // CMSY
                proj_F(sp, yr) = Ftarget(sp);
                break;

              case 2: // Constant F
                proj_F(sp, yr) = Ftarget(sp);
                break;

              case 3: // Constant F that acheives X% of SSB0
                proj_F(sp, yr) = Ftarget(sp);
                break;

              case 4: // Constant Fspr
                proj_F(sp, yr) = Ftarget(sp) * Fmult(sp);
                break;

              case 5: // NPFMC Tier 3 HCR
                proj_F(sp, yr) = Ftarget(sp); // Used Fabc of Ftarget_age%
                break;

              case 6: // PFMC Category 1 HCR
                proj_F(sp, yr) = Ftarget(sp) = Flimit(sp) + QnormHCR(sp);
                break;

              case 7: // SESSF Tier 1 HCR
                proj_F(sp, yr) = Ftarget(sp); // Used Fabc of Ftarget_age%
                break;
              }

              // No projected F when the forecast is off or numbers-at-age are input (estDynamics > 0).
              if((forecast(sp) == 0) | (estDynamics(sp) > 0)){
                proj_F(sp, yr) = 0;
              }
              F_flt_age(flt, sex, age, yr) = sel_at_age(flt, sex, age, yr) * proj_F_prop(flt) * proj_F(sp, yr);
            }

            // -- Sum F across fleets
            F_at_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);

            // -- Calculate F target by age and sex for reference points
            Flimit_at_age(sp, sex, age, yr) += sel_at_age(flt, sex, age, yr) * proj_F_prop(flt) * Flimit(sp); // account for time-varying sel
            Ftarget_at_age(sp, sex, age, yr) += sel_at_age(flt, sex, age, yr) * proj_F_prop(flt) * Ftarget(sp); // account for time-varying sel
          }
        }
      }

      // F across fleets or species
      for(yr = 0; yr < nyrs; yr++) {
        // Hindcast
        if( yr < nyrs_hind){
          F_flt(flt, yr) = exp(log_F(flt, yr));
          F_spp(sp, yr) += exp(log_F(flt, yr)); // Fully selected fishing mortality
        }

        // Forecast
        if( yr >= nyrs_hind){
          F_flt(flt, yr) = proj_F_prop(flt) * proj_F(sp, yr);
          F_spp(sp, yr) +=  proj_F_prop(flt) * proj_F(sp, yr);
        }
      }
    }
  }


  // 5.12.1. FISHERY SELECTIVITY AT AGE FOR THE INITIAL STATE (initMode 6 only)
  // SS3's InitF convention: the initial equilibrium decays at Finit weighted by fishery
  // selectivity, so a size-selective fishery barely touches the young ages. Mean over the
  // species' fishery fleets in year 1, and 1 everywhere for every other mode. Filled
  // before the SPR reference points build Z_init, and so SPRFinit, from it.
  for(sp = 0; sp < nspp; sp++){
    for(sex = 0; sex < nsex(sp); sex++){
      for(age = 0; age < nages(sp); age++){
        sel_init(sp, sex, age) = 1.0;
      }
    }
  }
  if(initMode == 6){
    for(sp = 0; sp < nspp; sp++){
      int n_fsh = 0;
      for(flt = 0; flt < n_flt; flt++){
        if((flt_type(flt) == 1) && (flt_spp(flt) == sp)) n_fsh++;
      }
      if(n_fsh > 0){
        for(sex = 0; sex < nsex(sp); sex++){
          for(age = 0; age < nages(sp); age++){
            Type acc = 0.0;
            for(flt = 0; flt < n_flt; flt++){
              if((flt_type(flt) == 1) && (flt_spp(flt) == sp)){
                acc += sel_at_age(flt, sex, age, 0);
              }
            }
            sel_init(sp, sex, age) = acc / Type(n_fsh);
          }
        }
      }
    }
  }


  // 5.13. SIMULATE PROCESS ERROR (sim_mod(simulate = TRUE), simulate_state)
  // Drawn before the dynamics consume the deviations -- the initial-numbers block below
  // is the first reader of rec_dev and init_dev -- so the dynamics and every observation
  // draw downstream are consistent with the simulated process. No second copy of the
  // population dynamics is needed (WHAM's sim_pop), these deviations and their SDs all
  // being parameters. Each draw uses the density that scores it in section 13, bias
  // correction included, which is what makes simulated recruitment mean-unbiased.
  //
  // True if exactly one density scores rec_dev. False for the AMAK/Ianelli configuration
  // (srr_fun == 0 with srr_pred_fun > 0), where the stock-recruit curve is fitted as a
  // second penalty on the same deviation. REPORTed so sim_mod() can warn without
  // re-deriving it.
  int rec_srr_single_density = !((srr_fun == 0) && (srr_pred_fun > 0));
  REPORT(rec_srr_single_density);

  SIMULATE {
    // Masks marking which cells this call actually drew. Every deviation array
    // is REPORTed whole -- hindcast and projection alike -- so without these a
    // caller cannot tell a simulated cell from one left at its fitted value, and
    // a recovery statistic taken over the whole array scores both.
    matrix<Type> rec_dev_drawn_sim(rec_dev.rows(), rec_dev.cols());
    matrix<Type> init_dev_drawn_sim(init_dev.rows(), init_dev.cols());
    array<Type>  log_M1_dev_drawn_sim(log_M1_dev.dim);
    rec_dev_drawn_sim.setZero();
    init_dev_drawn_sim.setZero();
    log_M1_dev_drawn_sim.setZero();

    for(sp = 0; sp < nspp; sp++){

      // The density covers the hindcast only -- projection recruitment comes from the
      // harvest control rule or sample_rec() -- so the fitted window is the only period
      // this can honestly redraw. srr_fun == 0 with srr_pred_fun > 0
      // (rec_srr_single_density, above) is excluded: in that AMAK/Ianelli configuration
      // section 13 scores rec_dev through JNLL_REC_DEV and again through
      // JNLL_SRR_PENALTY, so there is no single distribution to draw from. Nothing is
      // drawn, sim_mod() says why, and the full argument is in
      // vignette("model-diagnostics").
      //
      // Excluded under a DSEM (dsem_on == 1). There, rec_dev is not an
      // independent deviation at all: section 5.5b OVERWRITES it from the GMRF
      // latent states, so this draw would run first and be discarded -- except
      // that it runs AFTER 5.5b in the same pass, so it would instead discard
      // the DSEM. Drawing IID here would silently throw away the SEM's lagged
      // paths, covariate effects and cross-species structure, and a self test
      // would report recovery of a DSEM from data generated by an IID process.
      // The scale would be wrong too: under a DSEM R_sd is the GMRF's
      // CONDITIONAL (innovation) SD, so an IID draw using it understates the
      // marginal variance by 1/(1-rho^2). The correct draw is x_tj ~ GMRF(Q)
      // followed by recomputing rec_dev, and that is what happens instead:
      // sim_mod() draws the field from the reported precision and writes it into
      // dsem_x_tj, so section 5.5b derives rec_dev from the drawn states and
      // applies the bias correction once, where it always is.
      // dsem_drives_rec, not dsem_on: where the SEM drives some other process,
      // rec_dev is an ordinary estimated deviation and this is the density that
      // scores it, so it is redrawn exactly as in a model with no SEM.
      if(simulate_state(0) == 1 && simulate_period(0) == 1 && rec_srr_single_density &&
         dsem_drives_rec == 0){
        for(yr = 0; yr < nyrs_hind; yr++){
          rec_dev(sp, yr) = rnorm(-bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp));
          rec_dev_drawn_sim(sp, yr) = 1;
        }
      }

      // Initial abundance-at-age: the same process, for years before styr. The
      // gate is the density's own -- the equilibrium modes and OffsetEquilibrium
      // (5) fix init_dev and carry no penalty, so there is nothing to draw -- and
      // the rec_srr_single_density gate as well, a fresh initial age structure on
      // top of FITTED hindcast deviations being no history the model generated.
      //
      // Drawn under a DSEM too. init_dev is NOT a latent state of the GMRF --
      // the field spans styr..endyr, and section 13 scores init_dev with its own
      // N(-bias*R_sd^2/2, R_sd) whatever dsem_on is, with R_sd taken from the
      // SEM's variance path -- so this is still the density that scores it. The
      // DSEM redraws the hindcast deviations from R, so leaving init_dev alone
      // is what would break "redrawn whole or not at all": the replicate's
      // initial age structure would come from one realization and its hindcast
      // from another.
      if(simulate_state(0) == 1 && simulate_period(0) == 1 && rec_srr_single_density &&
         (initMode > 1) && (initMode != 5)){
        for(age = 1; age < nages(sp); age++){
          init_dev(sp, age - 1) = rnorm(-bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp));
          init_dev_drawn_sim(sp, age - 1) = 1;
        }
      }
    }

    // Natural mortality random effects. The IID modes are the AR1 code at rho = 0, so one
    // construction covers all six: draw a standardized AR1 and scale by the marginal SD
    // afterwards. M1_dev_log_sd is the INNOVATION sd; Sigma_M is the marginal.
    //
    // The draw is BROADCAST along whichever dimension the mode holds constant (modes 1/4
    // over years, modes 2/5 over ages) and across the sexes sharing one parameter --
    // build_map() separates the sexes only at M1_model = 2, so num_re_sexes is the number
    // of INDEPENDENT draws. This block writes the expanded array, so writing only the
    // cell the density reads would leave the rest at fitted values and the deviation
    // would not be constant.
    if(simulate_state(1) == 1 && simulate_period(0) == 1){
      for(sp = 0; sp < nspp; sp++){
        int num_re_sexes = (M1_model(sp) == 2 && nsex(sp) > 1) ? 2 : 1;

        if((M1_re(sp) == 1) || (M1_re(sp) == 4)){        // by age, constant over years
          Type sigma_M = exp(M1_dev_log_sd(sp, 0));
          Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
          Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_a, 2)), 0.5);
          for(int s = 0; s < num_re_sexes; s++){
            vector<Type> M_re_age(nages(sp));
            AR1(rho_M_a).simulate(M_re_age);
            M_re_age *= Sigma_M;
            int sex_end = (num_re_sexes == 1) ? nsex(sp) : s + 1;
            for(int sex = s; sex < sex_end; sex++){
              for(age = 0; age < nages(sp); age++){
                for(yr = 0; yr < nyrs; yr++){
                  log_M1_dev(sp, sex, age, yr) = M_re_age(age);
                  log_M1_dev_drawn_sim(sp, sex, age, yr) = 1;
                }
              }
            }
          }
        }

        if((M1_re(sp) == 2) || (M1_re(sp) == 5)){        // by year, constant over ages
          Type sigma_M = exp(M1_dev_log_sd(sp, 0));
          Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
          Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_y, 2)), 0.5);
          for(int s = 0; s < num_re_sexes; s++){
            vector<Type> M_re_yr(nyrs_hind);
            AR1(rho_M_y).simulate(M_re_yr);
            M_re_yr *= Sigma_M;
            int sex_end = (num_re_sexes == 1) ? nsex(sp) : s + 1;
            for(int sex = s; sex < sex_end; sex++){
              for(yr = 0; yr < nyrs_hind; yr++){
                for(age = 0; age < nages(sp); age++){
                  log_M1_dev(sp, sex, age, yr) = M_re_yr(yr);
                  log_M1_dev_drawn_sim(sp, sex, age, yr) = 1;
                }
              }
            }
          }
        }

        if((M1_re(sp) == 3) || (M1_re(sp) == 6)){        // by age and year
          Type sigma_M = exp(M1_dev_log_sd(sp, 0));
          Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
          Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
          Type Sigma_M = pow(pow(sigma_M, 2) / ((1.0 - pow(rho_M_y, 2)) * (1.0 - pow(rho_M_a, 2))), 0.5);
          for(int s = 0; s < num_re_sexes; s++){
            array<Type> M_re_a_yr(nages(sp), nyrs_hind);
            // Same argument order as the density in section 14 -- year on the
            // outermost dimension, age on the fastest-running one.
            SEPARABLE(AR1(rho_M_y), AR1(rho_M_a)).simulate(M_re_a_yr);
            int sex_end = (num_re_sexes == 1) ? nsex(sp) : s + 1;
            for(int sex = s; sex < sex_end; sex++){
              for(age = 0; age < nages(sp); age++){
                for(yr = 0; yr < nyrs_hind; yr++){
                  log_M1_dev(sp, sex, age, yr) = Sigma_M * M_re_a_yr(age, yr);
                  log_M1_dev_drawn_sim(sp, sex, age, yr) = 1;
                }
              }
            }
          }
        }
      }
    }

    // Reported so the draw can be checked directly. Inferring it from the data
    // it produces cannot tell "drew correctly" from "drew nothing, and the
    // observation error moved". Each goes out with its drawn-cell mask.
    matrix<Type> rec_dev_sim = rec_dev;
    matrix<Type> init_dev_sim = init_dev;
    array<Type>  log_M1_dev_sim = log_M1_dev;
    REPORT(rec_dev_sim);
    REPORT(init_dev_sim);
    REPORT(log_M1_dev_sim);
    REPORT(rec_dev_drawn_sim);
    REPORT(init_dev_drawn_sim);
    REPORT(log_M1_dev_drawn_sim);
  }

  /** ------------------------------------------------------------------------ //
   * 6. POPULATION DYNAMICS EQUATIONS                                          //
   * ------------------------------------------------------------------------- */
  // Start iterations for multi-species convergence (anything dependent on M)
  for(int iter = 0; iter < niter; iter++) {

    // 6.1. TOTAL MORTALITY-AT-AGE
    for(sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        // Element-wise: assigning a scalar to a vector would broadcast.
        for(int i = 0; i < M1_beta.dim(2); i++){
          beta_M1_tmp(i) = M1_beta(sp, sex, i);
        }
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs; yr++) {
            // Matrix multiplication from sliced arrays doesn't work
            env_M1_tmp = env_index.row(yr);
            M1_mult = env_M1_tmp * beta_M1_tmp;

            // M1_mult.sum() is LEGACY (scheduled removal: v4.5.0). It carries
            // the environmental effect of the soft-deprecated M1_indices path,
            // so fits built that way still reproduce; `linkages` is the current
            // route. See "Scheduled removal" in NEWS.md (4.1.0 section).
            M1_at_age(sp, sex, age, yr) = exp(
              log_M1(sp, sex, age)
              + log_M1_dev(sp, sex, age, yr)
              + M1_mult.sum()
              + M_linkage_offset(sp, sex, age, yr)
            ) + M_linkage_offset_nat(sp, sex, age, yr);
            M_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            M_at_age_dB0(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age_dB0(sp, sex, age, yr);
            M_at_age_dBF(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age_dBF(sp, sex, age, yr);
            Z_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + F_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
          }
        }
      }
    }


    // 6.2. SPR BASED REFERENCE POINTS
    //FIXME - make time-varying?
    SPR0.setZero();
    SPRFinit.setZero();
    SPRlimit.setZero();
    SPRtarget.setZero();
    // Not defined under predation: M carries M2, which scales with predator abundance, so
    // spawning output per recruit would not be a property of the prey stock alone.
    // data_check() refuses the switch combinations that read these downstream.
    if(msmMode == 0){
      for(sp = 0; sp < nspp; sp++) {

        int na = nages(sp);
        int term_yr = nyrs_hind - 1;      // rates for a reference point are the terminal hindcast year's
        wt_idx_ssb = 2 * sp + 1;

        // The four schedules. Flimit and Ftarget are selectivity-scaled and read the
        // terminal hindcast year; Finit is the flat F that generated the initial age
        // structure, so it is constant across ages and reads year 0, weight included.
        // FIXME: time-vary sel in the forecast
        vector<Type> Z_unfished(na), Z_limit(na), Z_target(na), Z_init(na);
        vector<Type> wt_term(na), wt_first(na), mature_at_age(na);
        Type female_split = (nsex(sp) == 1) ? Type(1.0) : sex_ratio(sp, 0);

        for(age = 0; age < na; age++){
          Z_unfished(age) = M_at_age(sp, 0, age, term_yr);
          Z_limit(age)    = M_at_age(sp, 0, age, term_yr) + Flimit_at_age(sp, 0, age, term_yr);
          Z_target(age)   = M_at_age(sp, 0, age, term_yr) + Ftarget_at_age(sp, 0, age, term_yr);
          Z_init(age)     = M_at_age(sp, 0, age, 0) + Finit(sp) * sel_init(sp, 0, age);

          // Spawning output per TOTAL recruit, so the female fraction enters once:
          // spawn_output (5.7) carries maturity and the age-varying ratio for a
          // one-sex species, female_split the recruitment split (6.6) for a two-sex one.
          wt_term(age)      = spawn_output(sp, age, term_yr);
          wt_first(age)     = spawn_output(sp, age, 0);
          mature_at_age(age) = female_split;
        }

        vector<Type> n_unfished = per_recruit_survivors(Z_unfished);
        vector<Type> n_limit    = per_recruit_survivors(Z_limit);
        vector<Type> n_target   = per_recruit_survivors(Z_target);
        vector<Type> n_init     = per_recruit_survivors(Z_init);

        SPR0(sp)      = spawning_biomass_per_recruit(n_unfished, Z_unfished, wt_term,  mature_at_age, spawn_month(sp));
        SPRlimit(sp)  = spawning_biomass_per_recruit(n_limit,    Z_limit,    wt_term,  mature_at_age, spawn_month(sp));
        SPRtarget(sp) = spawning_biomass_per_recruit(n_target,   Z_target,   wt_term,  mature_at_age, spawn_month(sp));
        SPRFinit(sp)  = spawning_biomass_per_recruit(n_init,     Z_init,     wt_first, mature_at_age, spawn_month(sp));

        // Reported for inspection; slot order matches the four calls above.
        for(age = 0; age < na; age++){
          NbyageSPR(0, sp, age) = n_unfished(age);
          NbyageSPR(1, sp, age) = n_limit(age);
          NbyageSPR(2, sp, age) = n_target(age);
          NbyageSPR(3, sp, age) = n_init(age);
        }
      }
    }


    // 6.3. INITIAL RECRUITMENT
    // -- For beverton-holt, steepness and R0 are derived from SPR0
    penalty = 0.0;
    zero_N_pen.setZero();
    for( sp = 0; sp < nspp ; sp++) {
      // Under predation SPR is undefined: a hindcast curve takes R_init from the
      // free level R0(sp, 0) and reports steepness 0, skipping the cases below.
      if((msmMode > 0) & (srr_fun > 1)){
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0; }
        R_init(sp) = R0(sp, 0);
      } else switch(srr_fun){
      case 0: // Random about mean (e.g. Alaska)
        // No compensation, so steepness is constant across years.
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0.99; }
        {
          R_init(sp) = R0(sp, 0);
        }
        break;

      case 1: // Random about mean with environmental linkage
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0.99; }
        {
          R_init(sp) = R0(sp, 0);
        }
        break;

      case 2: // Beverton-Holt
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            steepness(sp, yr) = alpha(sp, yr) * SPR0(sp)/(4.0 + alpha(sp, yr) * SPR0(sp));
          }
          // Unfished recruitment implied by the curve, first year only. Filling every
          // year would make R0 an AD function of alpha and Beta across the whole series
          // -- ~6x the fit time for a value that is reported rather than used,
          // Beverton-Holt recruitment being a function of alpha, Beta and SSB, not of R0.
          R0(sp, 0) = (alpha(sp, 0) - 1.0/SPR0(sp)) / Beta(sp, 0); // (Alpha-1/SPR0)/beta
          R_init(sp) = (alpha(sp, 0) - 1.0/SPRFinit(sp)) / Beta(sp, 0); // (Alpha-1/SPRFinit)/beta
          // Same floor as R in 6.6.1: an identity-link alpha offset can make both negative.
          if (rec_floor_on) { Type pen_R = 0; R0(sp, 0) = posfun(R0(sp, 0), Type(1e-3), pen_R); R_init(sp) = posfun(R_init(sp), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
        }
        break;

      case 3: // Retired alias of 2; build_srr() refuses it (alpha moves via a linkage)
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            steepness(sp, yr) = alpha(sp, yr) * SPR0(sp)/(4.0 + alpha(sp, yr) * SPR0(sp));
          }
          // Unfished recruitment implied by the curve, first year only -- see
          // case 2 for why it is not filled for every year.
          R0(sp, 0) = (alpha(sp, 0) - 1.0/SPR0(sp)) / Beta(sp, 0); // (Alpha-1/SPR0)/beta
          R_init(sp) = (alpha(sp, 0) - 1.0/SPRFinit(sp)) / Beta(sp, 0); // (Alpha-1/SPRFinit)/beta
          // Same floor as R in 6.6.1: an identity-link alpha offset can make both negative.
          if (rec_floor_on) { Type pen_R = 0; R0(sp, 0) = posfun(R0(sp, 0), Type(1e-3), pen_R); R_init(sp) = posfun(R_init(sp), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
        }
        break;

      case 4: // Ricker
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            { Type aS = alpha(sp, yr) * SPR0(sp);   // kept positive under an identity-link alpha offset
              if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
              steepness(sp, yr) = 0.2 * exp(0.8*log(aS)); }
          }

          // - R at F0
          Type pen_rk = 0;   // this species' floor penalties alone
          ricker_intercept = alpha(sp, 0) * SPR0(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          // Year 0 only: posfun() keeps the Ricker intercept positive, charging the
          // excursion to the objective. The Ianelli configuration (srr_fun < 2)
          // never reaches this case, so its R0 keeps its estimated value.
          R0(sp, 0) = log(ricker_intercept)/(Beta(sp, 0) * SPR0(sp)/1000000.0);

          // R at equilibrium F
          ricker_intercept = alpha(sp, 0) * SPRFinit(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          R_init(sp) = log(ricker_intercept)/(Beta(sp, 0) * SPRFinit(sp)/1000000.0);
          zero_N_pen(sp) += pen_rk;
        }
        break;

      case 5: // Ricker with environmental impacts on alpha
        {
          // Steepness for every year: alpha may be time-varying through a linkage.
          for(yr = 0; yr < nyrs; yr++){
            { Type aS = alpha(sp, yr) * SPR0(sp);   // kept positive under an identity-link alpha offset
              if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
              steepness(sp, yr) = 0.2 * exp(0.8*log(aS)); }
          }

          // - R at F0
          Type pen_rk = 0;   // this species' floor penalties alone
          ricker_intercept = alpha(sp, 0) * SPR0(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          // Year 0 only, as case 4.
          R0(sp, 0) = log(ricker_intercept)/(Beta(sp, 0) * SPR0(sp)/1000000.0);

          // R at equilibrium F
          ricker_intercept = alpha(sp, 0) * SPRFinit(sp) - 1.0;
          ricker_intercept =  posfun(ricker_intercept, Type(0.001), pen_rk) + 1.0;
          R_init(sp) = log(ricker_intercept)/(Beta(sp, 0) * SPRFinit(sp)/1000000.0);
          zero_N_pen(sp) += pen_rk;
        }
        break;

      default:
        error("Invalid 'srr_fun'");
      }

      // Under the Ianelli configuration (srr_fun 0 or 1, srr_pred_fun > 1) the curve is
      // estimated as a recruitment penalty rather than driving the hindcast, and the
      // switch above -- keyed on srr_fun -- would leave steepness at the mean-recruitment
      // constant. Derive it from the penalty curve's alpha instead, since the
      // stock-recruit prior (JNLL_SRR_PRIOR) is evaluated on it. R0 is deliberately NOT
      // re-derived: recruitment here is R0 * exp(rec_dev), so R0 keeps its estimated
      // value.
      if((srr_fun < 2) & (srr_pred_fun > 1) & (msmMode == 0)){
        for(yr = 0; yr < nyrs; yr++){
          if((srr_pred_fun == 2) | (srr_pred_fun == 3)){
            steepness(sp, yr) = alpha(sp, yr) * SPR0(sp)/(4.0 + alpha(sp, yr) * SPR0(sp));
          }
          if((srr_pred_fun == 4) | (srr_pred_fun == 5)){
            { Type aS = alpha(sp, yr) * SPR0(sp);   // kept positive under an identity-link alpha offset
              if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
              steepness(sp, yr) = 0.2 * exp(0.8*log(aS)); }
          }
        }
      }
      // SPR0 is undefined under predation, so steepness is reported as 0;
      // data_check() refuses a steepness prior there.
      if((srr_fun < 2) & (srr_pred_fun > 1) & (msmMode > 0)){
        for(yr = 0; yr < nyrs; yr++){ steepness(sp, yr) = 0; }
      }
    }


    // 6.4. INITIAL ABUNDANCE AT AGE, BIOMASS, AND SSB (YEAR 1)
    biomass.setZero();
    ssb.setZero();
    for(sp = 0; sp < nspp; sp++) {

      // Sex ratio at recruitment (1-sex model gets all R)
      if(nsex(sp)  == 1){
        sex_ratio(sp, 0) = 1.0;
      }

      // A free initial recruitment level, log scale; 0 without an R_init linkage.
      // Read at year 0 because the initial state is one year, and per species
      // because that is the stratum it is estimated on. Only the estimated branch
      // reads it: input numbers-at-age have no initial state to scale.
      Type init_level = recruitment_linkage_offset(sp, RCEATTLE_REC_R_INIT, 0);

      for(age = 0; age < nages(sp); age++){
        for(sex = 0; sex < nsex(sp); sex ++){

          switch(estDynamics(sp)){
          case 0: // Estimated

            // - Amin (i.e. recruitment).
            if(age == 0){
              // Year-1 recruitment is a recruitment year with its own deviate,
              // not part of the pre-hindcast level, so it stays at R_init.
              R(sp, 0) = R_init(sp) * exp(rec_dev(sp, 0));
              N_at_age(sp, 0, 0, 0) = R(sp, 0) * sex_ratio(sp, 0);
              // Guarded because the male slot exists only where some species is
              // two-sex: the value would be 0 anyway (sex_ratio 1, set above),
              // but the index is out of range when max_sex is 1.
              if(nsex(sp) > 1){
                N_at_age(sp, 1, 0, 0) = R(sp, 0) * (1-sex_ratio(sp, 0));
              }
              // The equilibrium carries no recruitment deviation: it is the mean
              // level the deviates depart from. It DOES carry the initial level,
              // on every age including this one, so the equilibrium that
              // equil_catch_hat integrates sits at one recruitment level.
              N_eq(sp, 0, 0) = R_init(sp) * exp(init_level) * sex_ratio(sp, 0);
              if(nsex(sp) > 1){
                N_eq(sp, 1, 0) = R_init(sp) * exp(init_level) * (1-sex_ratio(sp, 0));
              }
            }

            // - Estimate  as free parameters
            if(initMode == 0){
              if(age > 0){
                N_at_age(sp, 0, age, 0) = exp(init_dev(sp, age-1)) * sex_ratio(sp, 0);
                if(nsex(sp) > 1){
                  N_at_age(sp, 1, age, 0) = exp(init_dev(sp, age-1)) * (1-sex_ratio(sp, 0));
                }
                // No equilibrium underlies free initial numbers; Finit is 0 here
                // and an equilibrium catch is refused, so this is never fitted.
                N_eq(sp, 0, age) = N_at_age(sp, 0, age, 0);
                if(nsex(sp) > 1){
                  N_eq(sp, 1, age) = N_at_age(sp, 1, age, 0);
                }
              }
            }

            // - Equilibrium or non-equilibrium from Rinit, Finit, mortality and init devs
            // Finit is 0 for initMode 0, 1, 2 and 5; it is estimated for the
            // fished non-equilibrium modes 3, 4 and 6 (see section 6.1).
            if(initMode > 0){

              // How far the initial age-structure sits off the mean-recruitment
              // equilibrium R0, on the log scale; 0 under neither source below. Two
              // sources, never both (an R_init linkage is refused under mode 5): the free
              // R_init level -- SS3's SR_regime block on the year before the hindcast,
              // which SS3 shrinks at sigmaR/ave_age and Rceattle leaves free, so put a
              // prior on the coefficient to match SS3 -- or, under OffsetEquilibrium (5),
              // the first-year recruitment deviation, seeding the structure off
              // exp(rec_pars + rec_dev(sp, 0)) with init devs off (Cole Monnahan / AFSC
              // GOA pollock convention).
              Type init_log_scalar = init_level;
              if(initMode == 5){
                init_log_scalar += rec_dev(sp, 0);
              }

              // Sum M1 until age - 1. OffsetEquilibrium (5) uses the same
              // standard departing-age cumulative-M decay as the equilibrium
              // modes (Finit is 0 for mode 5).
              if((initMode == 1) | (initMode == 2) | (initMode == 3) | (initMode == 5)){
                mort_sum(sp, age) = 0;
                for(int age_tmp = 0; age_tmp < age; age_tmp++){
                  mort_sum(sp, age) += M1_at_age(sp, sex, age_tmp, 0) + Finit(sp);
                }
              }

              // SS3's InitF convention: the initial F is spread over ages by
              // the fishery's selectivity before it accumulates, so an
              // unselected young age carries almost none of it.
              if(initMode == 6){
                mort_sum(sp, age) = 0;
                for(int age_tmp = 0; age_tmp < age; age_tmp++){
                  mort_sum(sp, age) += M1_at_age(sp, sex, age_tmp, 0)
                                     + Finit(sp) * sel_init(sp, sex, age_tmp);
                }
              }

              if(initMode == 4){
                mort_sum(sp, age) = 0;
                for(int age_tmp = 0; age_tmp < age; age_tmp++){
                  mort_sum(sp, age) += M1_at_age(sp, sex, age_tmp, 0);
                }
                mort_sum(sp, age) += Finit(sp);
              }

              // -- 6.5.2. Age Amin+1:Amax-1 (initial abundance)
              if((age > 0) & (age < nages(sp) - 1)) {

                if(sex == 0){
                  N_at_age(sp, 0, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) * sex_ratio(sp, 0);
                  N_eq(sp, 0, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) * sex_ratio(sp, 0);
                }
                if(sex == 1){
                  N_at_age(sp, 1, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) * (1-sex_ratio(sp, 0));
                  N_eq(sp, 1, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) * (1-sex_ratio(sp, 0));
                }
              }

              // -- 6.5.3. Amax
              if(age == (nages(sp) - 1)) {

                // The plus group accumulates at its own total mortality. Under initMode
                // 6 the initial F reaches an age through that age's selectivity, so the
                // oldest age decays at M1 + Finit * sel (SS3's equilibrium convention);
                // every other mode charges the full Finit, as in mort_sum above.
                Type Z_plus = M1_at_age(sp, sex, nages(sp) - 1, 0) + Finit(sp) *
                  ((initMode == 6) ? sel_init(sp, sex, nages(sp) - 1) : Type(1.0));
                Type plus_grp = 1 - exp(-Z_plus);

                if(sex == 0){// NOTE: This solves for the geometric series
                  N_at_age(sp, 0, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) / plus_grp * sex_ratio(sp, 0);
                  N_eq(sp, 0, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) / plus_grp * sex_ratio(sp, 0);
                }

                if(sex == 1){
                  N_at_age(sp, 1, age, 0) = R_init(sp) * exp( - mort_sum(sp, age) + init_dev(sp, age - 1) + init_log_scalar) / plus_grp * (1-sex_ratio(sp, 0));
                  N_eq(sp, 1, age)        = R_init(sp) * exp( - mort_sum(sp, age) + init_log_scalar) / plus_grp * (1-sex_ratio(sp, 0));
                }
              }
            }
            break;

          case 1: // Numbers-at-age fixed exactly to NByageFixed (pop_scalar mapped to
            // NA in build_map so log_pop_scalar = 0 -> pop_scalar = 1.0)
            N_at_age(sp, sex, age, 0) = pop_scalar(sp) * NByageFixed(sp, sex, age, 0);
            break;

          case 2: // Numbers-at-age scaled by a single estimated age-independent scalar
            N_at_age(sp, sex, age, 0) = pop_scalar(sp) * NByageFixed(sp, sex, age, 0);
            break;

          default:
            error("Invalid 'estDynamics'");
          }

          // -- 6.5.3. Estimate total biomass in year 1
          wt_idx_pop = 2 * sp;
          biomass_at_age(sp, sex, age, 0) = N_at_age(sp, sex, age, 0) * weight_hat( wt_idx_pop, sex, age, 0);
          biomass(sp, 0) += N_at_age(sp, sex, age, 0) * weight_hat( wt_idx_pop, sex, age, 0);
        }

        // -- 6.5.4. Estimated initial female SSB
        wt_idx_ssb = 2 * sp + 1;
        ssb(sp, 0) += N_at_age(sp, 0, age, 0) * exp(-Z_at_age(sp, 0, age, 0) * spawn_month(sp)/12) * spawn_output(sp, age, 0); // 6.6.
      }
    }


    // 6.5. HINDCAST NUMBERS AT AGE, BIOMASS-AT-AGE (mt), and ssb-AT-AGE (mt)
    penalty = 0.0;
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 1; yr < nyrs_hind; yr++) {

        // Switch for srr (for MSEs if using Ianelli SRR method)
        int srr_switch = srr_fun;
        if(yr >= nyrs_srrmean){
          srr_switch = srr_pred_fun;
        }

        // -- 6.6.1. Recruitment

        // Recruits arriving at age minage in year yr were spawned in yr - minage,
        // before styr while yr < minage. With no modelled SSB those years take
        // R_init * exp(rec_dev) -- equilibrium recruitment at F = Finit, what
        // year 0 gets, so ages seeded from the same equilibrium agree. Stock
        // Synthesis treats the pre-start period the same way; WHAM avoids it by
        // fixing the lag at one year. Cannot fire at minage = 1.
        int spawn_yr = yr - minage(sp);
        int srr_use = (spawn_yr < 0) ? 0 : srr_switch;
        Type ssb_tmp = (spawn_yr < 0) ? Type(0.0) : ssb(sp, spawn_yr);
        // Pre-styr years take R_init, not R0(sp, yr): under a hindcast curve R0 is
        // set only at year 0 (derived, or the free level under predation).
        Type rec_mean = (spawn_yr < 0) ? R_init(sp) : R0(sp, yr);

        R(sp, yr) = calculate_recruitment(srr_use, rec_mean, ssb_tmp, alpha(sp, yr), Beta(sp, yr), rec_dev(sp, yr), SPR0(sp));
        // Under an identity-link recruitment offset the curve can go non-positive: keep R
        // positive (posfun, 0.001 barrier), charging the excursion to the zero-N row.
        if (rec_floor_on) { Type pen_R = 0; R(sp, yr) = posfun(R(sp, yr), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }

        N_at_age(sp, 0, 0, yr) = R(sp, yr) * sex_ratio(sp, 0);
        if(nsex(sp) > 1){
          N_at_age(sp, 1, 0, yr) = R(sp, yr) * (1.0-sex_ratio(sp, 0));
        }


        // -- 6.6.2. Ages beyond recruitment
        for(age = 0; age < nages(sp); age++){
          for(sex = 0; sex < nsex(sp); sex ++){

            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- Where Amin < age < Amax
              if(age < (nages(sp) - 1)) {
                N_at_age(sp, sex, age + 1, yr) = N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));// S_at_age(sp, sex, age, yr - 1);
              }

              // -- Plus group where age = Amax.
              if(age == (nages(sp) - 1)) {
                N_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age - 1, yr - 1) * exp(-Z_at_age(sp, sex, age-1, yr-1)) + N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1)); // S_at_age(sp, sex, age, yr - 1);
              }
              break;

            case 1: // Numbers-at-age fixed exactly to NByageFixed (pop_scalar = 1.0 via map)
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            case 2: // Numbers-at-age scaled by a single estimated age-independent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            default:
              error("Invalid 'estDynamics'");
            }

            // A fresh accumulator per cell: the penalty is this cell's excursion alone.
            { Type pen_N = 0; N_at_age(sp, sex, age, yr) = posfun(N_at_age(sp, sex, age, yr), Type(0.001), pen_N); zero_N_pen(sp) += pen_N; }

            // -- 6.6.3. Estimate total biomass
            wt_idx_pop = 2 * sp ;
            biomass_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * weight_hat( wt_idx_pop, sex, age, yr );
            biomass(sp, yr) += biomass_at_age(sp, sex, age, yr);
          }

          // -- 6.6.4. Estimated female ssb
          wt_idx_ssb = 2 * sp + 1;
          ssb(sp, yr) += N_at_age(sp, 0, age, yr) * exp(-Z_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0) * spawn_output(sp, age, yr); // 6.6.
        }
      }
    }


    // 6.6. DEPLETION BASED BIOMASS REFERENCE POINTS (i.e. SB0 and dynamic SB0)
    // -- calculate mean recruitment
    // FIXME: may want specified year ranges
    avg_R.setZero();
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 0; yr < nyrs_srrmean; yr++) {
        avg_R(sp) += R(sp, yr)/Type(nyrs_srrmean);
      }
    }

    // -- Calc biomass_depletion based BRPs
    NByage0.setZero();
    NByageF.setZero();
    N_at_age_dBF.setZero();
    N_at_age_dB0.setZero();
    B0.setZero();
    SB0.setZero();
    SBF.setZero();
    DynamicB0.setZero();
    DynamicSB0.setZero();
    DynamicSBF.setZero();
    for(sp = 0; sp < nspp; sp++) {
      for(yr = 0; yr < nyrs; yr++) {

        // - Year 1
        // -- Initial abundance-at-age is the same as the hindcast
        if(yr == 0){
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 0; age < nages(sp); age++){
              NByage0(sp, sex, age , 0) = NByageF(sp, sex, age , 0) = N_at_age_dBF(sp, sex, age, 0) = N_at_age_dB0(sp, sex, age, 0) = N_at_age(sp,sex,age,0);
            }
          }
        }

        // - Recruitment Year > 1
        if(yr > 0){
          // - Option 1a: Use mean rec
          if((proj_mean_rec == 1) & (srr_pred_fun < 2)){
            // Equilibrium RPs used mean rec across all years
            NByage0(sp, 0, 0, yr) = avg_R(sp);
            NByageF(sp, 0, 0, yr) = avg_R(sp);

            // Dynamic RPs use mean rec for the forecast and observed rec for the
            // hindcast. LEGACY (scheduled removal: v4.5.0): prefer the
            // stock-recruit arm below. This split feeds SB0(sp, nyrs-1), so it
            // wants a deliberate decision rather than a sweep.
            if(yr < nyrs_hind){
              N_at_age_dBF(sp, 0, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) = R(sp, yr); // Hindcast use observed R
            }
            if(yr >= nyrs_hind){
              N_at_age_dBF(sp, 0, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) = exp(log(avg_R(sp)) + rec_dev(sp, yr)); // Projections use mean R given bias in R0
            }
          }


          // - Option 2: Use SRR
          // Reached whenever a curve exists, not only when the projection runs on it. The
          // breadth is deliberate: Option 1a excludes srr_pred_fun >= 2, so proj_mean_rec
          // = 1 (the build_srr() default) with an estimated curve would match NEITHER
          // arm, leaving reference-point recruitment 0 after the first year. SB0 would
          // decay towards zero, and SB0(sp, nyrs-1) is the depletion reference HCRs 5 and
          // 6 read. Projection recruitment is unaffected: 6.8 reads that switch
          // separately, so mean-recruitment projections stay mean.
          if((proj_mean_rec == 0) | (srr_pred_fun >= 2)){

            // - Equilibrium reference points (No recruitment deviation: pass Type(0.0))
            // Equilibrium spawning biomass exists every year; yr < minage reads the
            // first year's value so the lag stays in bounds.
            int rp_yr = yr - minage(sp);
            if(rp_yr < 0){ rp_yr = 0; }

            // TODO: not floored under an identity-link recruitment offset, so SB0, SBF and B0 can go negative (issue #234).
            NByage0(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), SB0(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));
            NByageF(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), SBF(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));

            // -  Dynamic reference points (Includes annual recruitment deviation: pass rdev)
            // Cohorts spawned before the first year take the hindcast's realized R
            // (R_init with its deviation), so both runs start from the same cohorts.
            Type rdev = rec_dev(sp, yr);
            if(yr < minage(sp)){
              Type R_early = (yr < nyrs_hind) ? R(sp, yr) : R_init(sp) * exp(rdev);
              N_at_age_dB0(sp, 0, 0, yr) = R_early;
              N_at_age_dBF(sp, 0, 0, yr) = R_early;
            } else {
              // Under the Ianelli form the hindcast is on mean recruitment until
              // srr_mse_switchyr, so rec_dev there is measured from R0; the curve takes
              // the realized deviation from it, log R - log R_hat.
              if((srr_fun != srr_pred_fun) & (yr < nyrs_srrmean)){
                Type R_curve = calculate_recruitment(srr_pred_fun, R0(sp, yr), ssb(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));
                // Same floor as R in 6.6.1, for an identity-link offset that makes the curve non-positive.
                if (rec_floor_on) { Type pen_R = 0; R_curve = posfun(R_curve, Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
                rdev = log(R(sp, yr)) - log(R_curve);
              }
              // TODO: dynamic-B0 recruitment is not floored under an identity-link offset (issue #234).
              N_at_age_dB0(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), DynamicSB0(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), rdev, SPR0(sp));
              N_at_age_dBF(sp, 0, 0, yr) = calculate_recruitment(srr_pred_fun, R0(sp, yr), DynamicSBF(sp, rp_yr), alpha(sp, yr), Beta(sp, yr), rdev, SPR0(sp));
            }

          } // End recruitment switch

          // Account for sex ratio. The male block is guarded because that slot
          // exists only where some species is two-sex; the value would be 0
          // anyway (sex_ratio 1, set in 6.4).
          NByage0(sp, 0, 0, yr) = NByage0(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          NByageF(sp, 0, 0, yr) = NByageF(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          N_at_age_dB0(sp, 0, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females
          N_at_age_dBF(sp, 0, 0, yr) = N_at_age_dBF(sp, 0, 0, yr) * sex_ratio(sp, 0); // Females

          if(nsex(sp) > 1){
            // Males (divide by sex_ratio because we multiplied in the lines above)
            NByage0(sp, 1, 0, yr) = NByage0(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
            NByageF(sp, 1, 0, yr) = NByageF(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
            N_at_age_dB0(sp, 1, 0, yr) = N_at_age_dB0(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
            N_at_age_dBF(sp, 1, 0, yr) = N_at_age_dBF(sp, 0, 0, yr) / sex_ratio(sp, 0) * (1.0-sex_ratio(sp, 0));
          }


          // N-at-age year > 1
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 1; age < nages(sp)-1; age++) {
              NByage0(sp, sex, age, yr) =  NByage0(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr-1)); // F = 0

              NByageF(sp, sex, age, yr) =  NByageF(sp, sex, age-1, yr-1) * exp(-M_at_age(sp, sex, age-1, yr-1) - Ftarget_at_age(sp, sex, age-1, yr-1)); // F = target

              // No 0.001 floor here: the hindcast's is a numerical guard, and dynamic B0
              // must not assume a stock that fell to it would not have crashed without
              // fishing.
              N_at_age_dB0(sp, sex, age, yr) =  N_at_age_dB0(sp, sex, age-1, yr-1) * exp(-M_at_age_dB0(sp, sex, age-1, yr - 1)); // F = 0

              N_at_age_dBF(sp, sex, age, yr) =  N_at_age_dBF(sp, sex, age-1, yr-1) * exp(-M_at_age_dBF(sp, sex, age-1, yr - 1) - Ftarget_at_age(sp, sex, age-1, yr-1)); // F = Ftarget
            }

            // Plus group  -- No fishing
            NByage0(sp, sex, nages(sp)-1, yr) = NByage0(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1))  + NByage0(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1));

            NByageF(sp, sex, nages(sp)-1, yr) = NByageF(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-2, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-2, yr-1))  + NByageF(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age(sp, sex, nages(sp)-1, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-1, yr-1));

            N_at_age_dB0(sp, sex, nages(sp)-1, yr) = N_at_age_dB0(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age_dB0(sp, sex, nages(sp)-2, yr - 1))  + N_at_age_dB0(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age_dB0(sp, sex, nages(sp)-1, yr - 1));

            N_at_age_dBF(sp, sex, nages(sp)-1, yr) = N_at_age_dBF(sp, sex, nages(sp)-2, yr - 1) * exp(-M_at_age_dBF(sp, sex, nages(sp)-2, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-2, yr-1))  + N_at_age_dBF(sp, sex, nages(sp)-1, yr - 1) * exp(-M_at_age_dBF(sp, sex, nages(sp)-1, yr - 1) - Ftarget_at_age(sp, sex, nages(sp)-1, yr-1));

          }
        }

        // A species with input numbers-at-age (estDynamics > 0) keeps them in the
        // dynamic runs; projection years are read directly, as 6.8 fills N later.
        if(estDynamics(sp) > 0){
          for(sex = 0; sex < nsex(sp); sex ++){
            for(age = 0; age < nages(sp); age++){
              Type n_fixed = (yr < nyrs_hind) ? N_at_age(sp, sex, age, yr) :
                pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              N_at_age_dB0(sp, sex, age, yr) = n_fixed;
              N_at_age_dBF(sp, sex, age, yr) = n_fixed;
            }
          }
        }


        // Calculate Dynamic SB0 and SB at F target
        for(age = 0; age < nages(sp); age++) {

          wt_idx_ssb = 2 * sp + 1;
          // Multispecies: M_at_age carries the projection's realized M2. The one rule
          // that reads SBF, NPFMC (HCR 5), is refused there; SB0 is replaced by MSSB0
          // below.
          SB0(sp, yr) +=  NByage0(sp, 0, age, yr) *  spawn_output(sp, age, nyrs_hind - 1) * exp(-M_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0);
          SBF(sp, yr) +=  NByageF(sp, 0, age, yr) *  spawn_output(sp, age, nyrs_hind - 1) * exp(-(M_at_age(sp, 0, age, yr) + Ftarget_at_age(sp, 0, age, yr)) * spawn_month(sp)/12.0);
          DynamicSB0(sp, yr) +=  N_at_age_dB0(sp, 0, age, yr) *  spawn_output(sp, age, yr) * exp(-M_at_age_dB0(sp, 0, age, yr) * spawn_month(sp)/12.0);
          DynamicSBF(sp, yr) +=  N_at_age_dBF(sp, 0, age, yr) *  spawn_output(sp, age, yr) * exp(-(M_at_age_dBF(sp, 0, age, yr) + Ftarget_at_age(sp, 0, age, yr)) * spawn_month(sp)/12.0);

          for(sex = 0; sex < nsex(sp); sex ++){

            wt_idx_pop = 2 * sp;
            B0(sp, yr) +=  NByage0(sp, sex, age, yr) *  weight_hat( wt_idx_pop, sex, age, nyrs_hind - 1 );
            DynamicB0(sp, yr) +=  N_at_age_dB0(sp, sex,age, yr) *  weight_hat( wt_idx_pop, sex, age, yr );
          }
        }


        // Input SB0 (if running in multi-species mode)
        if(msmMode > 0){
          SB0(sp, yr) = MSSB0(sp);
          B0(sp, yr) = MSB0(sp);
        }
      }
    }

    // 6.7-6.9. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (mt), and ssb-AT-AGE (mt)
    // Includes Harvest Control Rules
    for(sp = 0; sp < nspp; sp++) {
      for(yr = nyrs_hind; yr < nyrs; yr++){

        // 6.7. HARVEST CONTROL RULES FOR PROJECTION
        // Equilibrium or dynamic reference points (i.e. SB0 or dynamic SB0)
        Type ref_SB0 = (DynamicHCR == 1) ? DynamicSB0(sp, yr-1) : SB0(sp, nyrs-1);
        Type ref_SBF = (DynamicHCR == 1) ? DynamicSBF(sp, yr-1) : SBF(sp, nyrs-1);

        switch(HCR){
        case 0: // No fishing
          proj_F(sp, yr) = 0.0;
          break;

        case 1: // CMSY
          // proj_F(sp, yr) retains its incoming value
          break;

        case 2: // Constant F
          // proj_F(sp, yr) retains its incoming value
          break;

        case 3: // Constant F to acheive X% of SB0
          // proj_F(sp, yr) retains its incoming value
          break;

        case 4: // Constant Fspr with multiplier
          // proj_F(sp, yr) retains its incoming value
          break;

        case 5: // NPFMC Tier 3 HCR
          if(ssb(sp, yr-1) < ref_SBF){
            proj_F(sp, yr) = Ftarget(sp) * (((ssb(sp, yr-1)/ref_SBF)-Alpha(sp))/(1-Alpha(sp))); // Used Fabc of FtargetSPR%
          }
          if((ssb(sp, nyrs_hind-1) < ref_SB0 * Plimit(sp)) || (ssb(sp, yr-1) / ref_SBF < Alpha(sp))){ // If overfished
            proj_F(sp, yr) = 0.0;
          }
          break;

        case 6: // PFMC Category 1 HCR
          if(ssb(sp, yr-1) < ref_SB0 * Ptarget(sp)){
            proj_F(sp, yr) = (Flimit(sp) + QnormHCR(sp)) * (ref_SB0 * Ptarget(sp) * (ssb(sp, yr-1) - ref_SB0 * Plimit(sp))) / (ssb(sp, yr-1) * (ref_SB0 * (Ptarget(sp) - Plimit(sp))));
          }
          if(ssb(sp, yr-1) < ref_SB0 * Plimit(sp)){ // If overfished
            proj_F(sp, yr) = 0.0;
          }
          break;

        case 7: // SESSF Tier 1 HCR
          if(ssb(sp, yr-1) < ref_SB0 * Ptarget(sp)){
            proj_F(sp, yr) = Ftarget(sp) * ((ssb(sp, yr-1)/(ref_SB0 * Plimit(sp)))-1); // Used Fabc of FtargetSPR%
          }
          if(ssb(sp, yr-1) < ref_SB0 * Plimit(sp)){ // If overfished
            proj_F(sp, yr) = 0.0;
          }
          break;
        }


        // No projected F when the forecast is off or numbers-at-age are input (estDynamics > 0).
        if((forecast(sp) == 0) | (estDynamics(sp) > 0)){
          proj_F(sp, yr) =  0.0;
        }

        // -- Update F for the projection (account for selectivity and fleets)
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            F_at_age(sp, sex, age, yr) = 0.0;
          }
        }

        // -- Multiply F from HCR by selectivity and fleet proportion
        F_spp(sp, yr) = proj_F(sp, yr);
        for(flt = 0; flt < n_flt; flt++) {
          if((sp == flt_spp(flt)) & (flt_type(flt) == 1)){ // Fisheries only, as in 5.12
            F_flt(flt, yr) = proj_F_prop(flt) * proj_F(sp, yr);
            for(age = 0; age < nages(sp); age++) {
              for(sex = 0; sex < nsex(sp); sex ++){
                F_flt_age(flt, sex, age, yr) = sel_at_age(flt, sex, age, nyrs_hind - 1) * proj_F_prop(flt) * proj_F(sp, yr); // Projection F uses the terminal hindcast year's selectivity
                if(flt_type(flt) == 1){
                  F_at_age(sp, sex, age, yr) += F_flt_age(flt, sex, age, yr);
                }
              }
            }
          }
        }


        // -- Update mortality for forecast with new F
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            M_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
            Z_at_age(sp, sex, age, yr) = M1_at_age(sp, sex, age, yr) + F_at_age(sp, sex, age, yr) + M2_at_age(sp, sex, age, yr);
          }
        }


        // 6.8. FORECAST NUMBERS AT AGE, BIOMASS-AT-AGE (mt), and ssb-AT-AGE (mt)
        // -- 6.8.1. Forecasted recruitment
        // - Option 1: Use mean rec
        if(proj_mean_rec == 1){
          R(sp, yr) = avg_R(sp); //  // Projections use mean R given bias in R0
        }

        // - Option 2: Use SRR and rec devs
        if(proj_mean_rec == 0){

          // - Calculate recruitment (with linkage offsets pre-added).
          // Pre-styr spawning year takes R_init, as in the hindcast block;
          // reachable only if the hindcast is shorter than minage.
          int proj_spawn_yr = yr - minage(sp);
          int proj_srr_use = (proj_spawn_yr < 0) ? 0 : srr_pred_fun;
          Type ssb_tmp = (proj_spawn_yr < 0) ? Type(0.0) : ssb(sp, proj_spawn_yr);
          Type proj_rec_mean = (proj_spawn_yr < 0) ? R_init(sp) : R0(sp, yr);
          // TODO: projected recruitment is not floored under an identity-link offset (issue #234).
          R(sp, yr) = calculate_recruitment(proj_srr_use, proj_rec_mean, ssb_tmp, alpha(sp, yr), Beta(sp, yr), rec_dev(sp, yr), SPR0(sp));
        }

        N_at_age(sp, 0, 0 , yr) = R(sp, yr) * sex_ratio(sp, 0);
        if(nsex(sp) > 1){
          N_at_age(sp, 1, 0 , yr) = R(sp, yr) * (1 - sex_ratio(sp, 0));
        }


        // -- 6.8.2. Ages > recruitment
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            switch(estDynamics(sp)){
            case 0: // Estimated numbers-at-age

              // -- 6.8.2.  Where minage + 1 <= age < Ai
              if(age < (nages(sp) - 1)) {
                N_at_age(sp, sex, age + 1, yr) = N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));
              }

              // -- 6.8.3. Plus group where age > Ai.
              if(age == (nages(sp) - 1)) {
                N_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age - 1, yr - 1) * exp(-Z_at_age(sp, sex, age-1, yr-1)) + N_at_age(sp, sex, age, yr - 1) * exp(-Z_at_age(sp, sex, age, yr-1));
              }
              break;

            case 1: // Fixed numbers-at-age - fixed scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            case 2: // Fixed numbers-at-age age-independent scalar
              N_at_age(sp, sex, age, yr) = pop_scalar(sp) * NByageFixed(sp, sex, age, yr);
              break;

            default: // Wrong estDynamics
              error("Invalid 'estDynamics'");
            }

            // -- 6.8.4. FORECAST BIOMASS
            wt_idx_pop = 2 * sp;
            biomass_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * weight_hat( wt_idx_pop, sex, age, nyrs_hind-1 ); // 6.5.
            biomass(sp, yr) += biomass_at_age(sp, sex, age, yr);
          } // End sex loop

          // -- 6.8.5. FORECAST SSB (SUM ACROSS AGES)
          wt_idx_ssb = 2 * sp + 1;
          ssb(sp, yr) += N_at_age(sp, 0, age, yr) * exp(-Z_at_age(sp, 0, age, yr) * spawn_month(sp)/12.0) * spawn_output(sp, age, nyrs_hind-1); // 6.6.
        }
      }
    }


    // 6.9. EXPECTED RECRUITMENT
    // The function-scope `yr` is left at `nyrs` by the forecast loop above, so
    // the year-0 block uses the explicit `first_yr` constant: reading
    // R0/alpha/Beta at `yr` would run one column past the end and segfault
    // inside MakeADFun. The assertion pins those matrices to the year horizon.
    const int first_yr = 0;
    if (R0.cols() != nyrs || alpha.cols() != nyrs || Beta.cols() != nyrs) {
      error("Section 6.9: R0/alpha/Beta column count does not equal nyrs; "
            "recruitment-parameter matrices were sized inconsistently.");
    }
    for(sp = 0; sp < nspp; sp++) {

      // Year 1 (arent fit in likelihood)
      // Under predation there is no SPRFinit (6.2) to place the first year on
      // the curve, so it takes R_init, the anchor its realised recruitment uses.
      if((msmMode > 0) & (srr_pred_fun > 1)){
        R_hat(sp, first_yr) = R_init(sp);
      } else {
        switch(srr_pred_fun){
        case 0: // Random about mean (e.g. Alaska)
          R_hat(sp, first_yr) = R0(sp, first_yr);
          break;

        case 1: // Random about mean with environmental effects
          R_hat(sp, first_yr) = R0(sp, first_yr);
          break;

        case 2: // Beverton-Holt
          R_hat(sp, first_yr) = (alpha(sp, first_yr) - 1/SPRFinit(sp)) / Beta(sp, first_yr); // (Alpha-1/SPR0)/beta
          break;

        case 3: // Retired alias of 2; build_srr() refuses it (alpha moves via a linkage)
          R_hat(sp, first_yr) = (alpha(sp, first_yr) - 1/SPRFinit(sp)) / Beta(sp, first_yr); // (Alpha-1/SPR0)/beta
          break;

        case 4: // Ricker
          { Type aS = alpha(sp, first_yr) * SPRFinit(sp);   // kept positive so the log is finite; R_hat itself is floored below
            if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
            R_hat(sp, first_yr) = log(aS) / (Beta(sp, first_yr) * SPRFinit(sp)/1000000.0); }
          break;

        case 5: // Ricker with environmental impacts on alpha
          { Type aS = alpha(sp, first_yr) * SPRFinit(sp);   // kept positive so the log is finite; R_hat itself is floored below
            if (rec_floor_on) { Type pen_d = 0; aS = posfun(aS, Type(1e-3), pen_d); }
            R_hat(sp, first_yr) = log(aS) / (Beta(sp, first_yr) * SPRFinit(sp)/1000000.0); }
          break;
        default:
          error("Invalid 'srr_pred_fun'");
        }
      }
      // Same floor as R in 6.6.1: the first year's R_hat is reported and read by retrospective().
      if (rec_floor_on) { Type pen_R = 0; R_hat(sp, first_yr) = posfun(R_hat(sp, first_yr), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }

      // Year 1+
      for(yr = 1; yr < nyrs; yr++){
        // Pre-styr spawning year takes R_init, the SAME anchor the realised recruitment
        // above uses, so expected and realised differ only by rec_dev and the
        // stock-recruit penalty (13.1) scores the deviation alone -- these years have no
        // spawning biomass to carry information about the curve's level. Anchoring on the
        // curve instead (R_hat's own first year) leaves a mean-versus-curve level gap the
        // penalty reads as signal: under the Ianelli configuration (srr_fun mean,
        // srr_pred_fun BevertonHolt) a fixed 0.83 nats per guarded year, pulling rec_pars
        // against alpha and Beta from years with no data. Identical for
        // BevertonHolt/BevertonHolt and Ricker/Ricker, where R_init is already R_hat's
        // first-year value.
        int hat_spawn_yr = yr - minage(sp);
        int hat_srr_use = (hat_spawn_yr < 0) ? 0 : srr_pred_fun;
        Type ssb_tmp = (hat_spawn_yr < 0) ? Type(0.0) : ssb(sp, hat_spawn_yr);
        Type hat_rec_mean = (hat_spawn_yr < 0) ? R_init(sp) : R0(sp, yr);

        // Expected recruitment carries no deviation, hence Type(0.0).
        R_hat(sp, yr) = calculate_recruitment(hat_srr_use, hat_rec_mean, ssb_tmp, alpha(sp, yr), Beta(sp, yr), Type(0.0), SPR0(sp));
        // Same floor as R in 6.6.1: log(R_hat) enters the stock-recruit penalty.
        if (rec_floor_on) { Type pen_R = 0; R_hat(sp, yr) = posfun(R_hat(sp, yr), Type(1e-3), pen_R); zero_N_pen(sp) += pen_R; }
      }
    }



    // 6.11. ESTIMATE AVERAGE NUMBERS AT AGE
    for(sp = 0; sp < nspp; sp++) {
      for(sex = 0; sex < nsex(sp); sex ++){
        for(age = 0; age < nages(sp); age++) {
          for(yr = 0; yr < nyrs; yr++) {
          // MSVPA average N over the year.
            avgN_at_age(sp, sex, age, yr) = N_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) / Z_at_age(sp, sex, age, yr);
            avgN_at_age_dB0(sp, sex, age, yr) = N_at_age_dB0(sp, sex, age, yr) * (1 - exp(-M_at_age_dB0(sp, sex, age, yr))) / M_at_age_dB0(sp, sex, age, yr);
            avgN_at_age_dBF(sp, sex, age, yr) = N_at_age_dBF(sp, sex, age, yr) * (1 - exp(-M_at_age_dBF(sp, sex, age, yr) - Ftarget_at_age(sp, sex, age, yr))) / (M_at_age_dBF(sp, sex, age, yr) + Ftarget_at_age(sp, sex, age, yr));
          }
        }
      }
    }


    /** ------------------------------------------------------------------------ //
     // 7. MSVPA PREDATION MORTALITY EQUATIONS                                   //
     * ------------------------------------------------------------------------- */

    // START PREDATION
    if(msmMode > 0) {

      // NOTE -- LOOPING INDICES -- sp = species, age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length

      // 7.1. Empirical MSVPA Suitability
      calculate_msvpa_suitability(stom_div_bio, suitability, suit_other,
                                  nspp, nyrs, suit_styr, suit_endyr, nyrs_suit, msmMode,
                                  nsex, nages, suitMode, diet_prop, avgN_at_age, weight_hat, other_food_diet_prop);


      // 7.2. MSVPA PREDATION MORTALITY
      if((msmMode == 1) | (msmMode == 2)) {
        calculate_msvpa_predation(avail_food, avail_food_dB0, avail_food_dBF,
                                  M2_at_age, M2_at_age_dB0, M2_at_age_dBF, M2_prop,
                                  B_eaten, B_eaten_as_prey, diet_prop_hat, other_diet_prop_hat,
                                  nspp, nyrs, nsex, nages, msmMode,
                                  avgN_at_age, avgN_at_age_dB0, avgN_at_age_dBF,
                                  suitability, suit_other, weight_hat, consumption_at_age, other_food);
      }
    } // End 7. Predation mortality
  } // End population dynamics iterations



  /** ------------------------------------------------------------------------ //
   // 8. INDEX COMPONENTS EQUATIONS                                             //
   * ------------------------------------------------------------------------- */

  // -- 8.1. Index of abundance/biomass
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

    index = index_ctl(index_ind, 0) - 1;          // Survey index
    sp = index_ctl(index_ind, 1) - 1;             // Species
    flt_yr = index_ctl(index_ind, 2);             // Data year
    wt_idx_flt = nspp * 2 + index;                // Fleet weight index
    mo = index_n(index_ind, 0);                   // Month

    index_hat(index_ind) = 0.0;

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
    }
    else if(flt_yr < 0){
      flt_yr = -flt_yr - styr;
    }

    for(age = 0; age < nages(sp); age++) {
      for(sex = 0; sex < nsex(sp); sex++){

        // Index by fleet type. Month is not used for a fishery: a year-average
        // has no instant to be taken at. SS3 splits the same way on its
        // per-fleet survey timing (SS_expval.tpl: timing >= 0 -> exp(-Z*timing);
        // timing < 0, required of every fishing fleet -> (1 - exp(-Z))/Z).
        //
        // A seasonal fishery takes the same form. Its exact predictor over [t1, t1+D] is
        //     N_a * exp(-Zpre_a * t1) * (1 - exp(-Zwin_a * D)) / (Zwin_a * D)
        //     Zpre_a = M_a + other fleets' F_a      (this fishery has yet to fish)
        //     Zwin_a = M_a + F_a/D + other fleets' F_a
        // Since Zwin*D = M*D + F, the window factor tracks F almost as the annual
        // one does, and q absorbs the F-free exp(-M*t1) offset. Trend error
        // against the exact form is ~1.5% over F 0.05-0.8, against -29% to +33%
        // for the snapshot. TODO: fit the window form directly if a fleet needs it.
        Type n_index;
        if(flt_type(index) == 1){
          n_index = N_at_age(sp, sex, age, flt_yr) * (1.0 - exp( - Z_at_age(sp, sex, age, flt_yr))) / Z_at_age(sp, sex, age, flt_yr);
        } else {
          n_index = N_at_age(sp, sex, age, flt_yr) * exp( - (mo/12.0) * Z_at_age(sp, sex, age, flt_yr));
        }

        // Weight
        if(flt_units(index) == 1){
          index_hat(index_ind) += n_index * sel_at_age(index, sex, age, flt_yr) * weight_hat(wt_idx_flt, sex, age, flt_yr );
        }
        // Numbers
        if(flt_units(index) == 2){
          index_hat(index_ind) += n_index * sel_at_age(index, sex, age, flt_yr);
        }
      }
    }
  }


  // -- 8.2. Analytical survey q following Walters and Ludwig 1994
  index_n_obs.setZero();
  index_q_analytical.setZero();
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){



    index = index_ctl(index_ind, 0) - 1;            // Survey index
    sp = index_ctl(index_ind, 1) - 1;             // Species
    flt_yr = index_ctl(index_ind, 2);             // Data year

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;


      mo = index_n(index_ind, 0);                    // Month
      if(flt_yr < nyrs_hind){

        // If etimated standard deviation or analytical sigma (non-time-varying)
        if(est_sigma_index(index) > 0) {
          index_n_obs(index) += 1;
          index_q_analytical(index) += log(index_obs(index_ind, 0) / index_hat(index_ind));
        }

        // If time-varying sigma
        if(est_sigma_index(index) == 0 ) {
          index_n_obs(index) += 1 / square(index_obs(index_ind, 1));
          index_q_analytical(index) += log(index_obs(index_ind, 0) / index_hat(index_ind)) / square(index_obs(index_ind, 1));
        }
      }
    }
  }

  // Take average
  for(index = 0 ; index < n_flt; index ++){
    index_q_analytical(index) = exp(index_q_analytical(index) / index_n_obs(index));

    // Set index_q to analytical if used
    if(est_index_q(index) == 3){
      for(yr = 0; yr < nyrs_hind; yr++){
        index_q(index, yr) = index_q_analytical(index);
      }
    }
  }


  // -- 8.2b. Arithmetic-mean analytical q (AMAK/ebswp): q = sum(obs)/sum(pred)
  //    over the fitted survey years (est_index_q == 7, "AnalyticalArith"). This
  //    reproduces the AMAK BTS q_bts = mean(ob_bts)/mean(eb_bts) and is the q
  //    form paired with the MVN covariance survey likelihood (Index_loglike =
  //    "MVN"). Uses index_hat before it is scaled by q just below.
  vector<Type> index_obs_sum(n_flt); index_obs_sum.setZero();
  vector<Type> index_hat_sum(n_flt); index_hat_sum.setZero();
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){
    index = index_ctl(index_ind, 0) - 1;
    flt_yr = index_ctl(index_ind, 2);
    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
      if(flt_yr < nyrs_hind){
        index_obs_sum(index) += index_obs(index_ind, 0);
        index_hat_sum(index) += index_hat(index_ind);
      }
    }
  }
  for(index = 0; index < n_flt; index++){
    if(est_index_q(index) == 7){
      Type q_arith = index_obs_sum(index) / index_hat_sum(index);
      for(yr = 0; yr < nyrs_hind; yr++){
        index_q(index, yr) = q_arith;
      }
    }
  }


  // -- 8.3. Survey Biomass - multiply by q
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

    index = index_ctl(index_ind, 0) - 1;            // Survey index
    flt_yr = index_ctl(index_ind, 2);      // Data year

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
    }
    else if(flt_yr < 0){
      flt_yr = -flt_yr - styr;
    }

    // Hindcast or projection indexing for q
    int yr_ind = (flt_yr < nyrs_hind) ? flt_yr : (nyrs_hind - 1);

    index_hat(index_ind) = index_q(index, yr_ind) * index_hat(index_ind); // pow(index_hat(index_ind), (1 + index_q_pow(index)));
  }



  // --8.4. Calculate analytical sigma following Walters and Ludwig 1994
  // The positive-observation guard matches the catch estimator below; it is
  // belt-and-braces, data_check() having already refused a non-positive index.
  index_n_obs.setZero();
  index_analytical_sd.setZero();
  for(index_ind = 0; index_ind < index_ctl.rows(); index_ind++){

    index = index_ctl(index_ind, 0) - 1;            // Survey index
    flt_yr = index_ctl(index_ind, 2);      // Data year

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;

      if((flt_yr < nyrs_hind) && (index_obs(index_ind, 0) > 0)){
        index_n_obs(index) += 1;
        index_analytical_sd(index) += square( log(index_obs(index_ind, 0)) - log(index_hat(index_ind)));
      }
    }
  }

  for(index = 0 ; index < n_flt; index ++){
    if(index_n_obs(index) > 0){
      index_analytical_sd(index) = concentrated_lognormal_sd(index_analytical_sd(index) / index_n_obs(index), bias_adjust_obs);
    }
  }


  /** ------------------------------------------------------------------------ //
   // 9. FISHERY COMPONENTS EQUATIONS                                          //
   * ------------------------------------------------------------------------- */
  // 9.1. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (mt or thousands of fish)
  for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){

    flt = catch_ctl(fsh_ind, 0) - 1;     // Fleet index (0-based)
    sp = catch_ctl(fsh_ind, 1) - 1;      // Species
    flt_yr = catch_ctl(fsh_ind, 2);      // Data year
    mo = catch_n(fsh_ind, 0);            // Month
    wt_idx_flt = nspp * 2 + flt;         // Weight index for fleet

    catch_hat(fsh_ind) = 0.0;
    max_catch_hat(fsh_ind) = 0.0;

    if(flt_yr > 0){
      flt_yr = flt_yr - styr;
    }
    else if(flt_yr < 0){
      flt_yr = -flt_yr - styr;
    }

    for(sex = 0; sex < nsex(sp); sex ++){
      for(age = 0; age < nages(sp); age++) {

        // By weight
        if(flt_units(flt) == 1){
          catch_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Z_at_age(sp, sex, age, flt_yr) * (1.0 - exp(-Z_at_age(sp, sex, age, flt_yr))) * N_at_age(sp, sex, age, flt_yr) * weight_hat( wt_idx_flt, sex, age, flt_yr ); // 5.5.
          max_catch_hat(fsh_ind) += N_at_age(sp, sex, age, flt_yr) * weight_hat( wt_idx_flt, sex, age, flt_yr ) * sel_at_age(flt, sex, age, flt_yr) * proj_F_prop(flt);
        }

        // By numbers
        if(flt_units(flt) == 2){
          catch_hat(fsh_ind) += F_flt_age(flt, sex, age, flt_yr) / Z_at_age(sp, sex, age, flt_yr) * (1.0 - exp(-Z_at_age(sp, sex, age, flt_yr))) * N_at_age(sp, sex, age, flt_yr);
          max_catch_hat(fsh_ind) += N_at_age(sp, sex, age, flt_yr) * sel_at_age(flt, sex, age, flt_yr) * proj_F_prop(flt);
        }
      }
    }
  }

  // -- 9.1a. Initial equilibrium catch, in weight or numbers per the fleet's
  //    Observation_units, as catch_hat is:
  //   C_eq = sum_a Finit * s_a * w_a * N_eq_a * (1 - exp(-Z_a)) / Z_a
  //   Z_a  = M1_a + Finit * s_a
  // Baranov on the deviation-free equilibrium age structure, at the fleet's own
  // selectivity and weight in the first hindcast year; SS3's Equil_catch
  // (SS_popdyn.tpl, Do_Equil_Calc). R supplies these rows only under initMode 6,
  // the one mode that builds N_eq at this same Finit * selectivity.
  vector<Type> equil_catch_hat(equil_catch_ctl.rows()); equil_catch_hat.setZero();
  for(int eq_ind = 0; eq_ind < equil_catch_ctl.rows(); eq_ind++){
    flt = equil_catch_ctl(eq_ind, 0) - 1;
    sp  = equil_catch_ctl(eq_ind, 1) - 1;
    wt_idx_flt = nspp * 2 + flt;
    for(sex = 0; sex < nsex(sp); sex++){
      for(age = 0; age < nages(sp); age++){
        Type s_a = sel_at_age(flt, sex, age, 0);
        Type Z_a = M1_at_age(sp, sex, age, 0) + Finit(sp) * s_a;
        Type baranov = Finit(sp) * s_a / Z_a * (1.0 - exp(-Z_a)) * N_eq(sp, sex, age);
        if(flt_units(flt) == 1){          // by weight
          equil_catch_hat(eq_ind) += baranov * weight_hat( wt_idx_flt, sex, age, 0 );
        }
        if(flt_units(flt) == 2){          // by numbers
          equil_catch_hat(eq_ind) += baranov;
        }
      }
    }
  }

  // -- 9.1b. Analytical catch sigma, following Walters and Ludwig (1994)
  // The sd that minimises the catch density: squared log residuals accumulated
  // here, concentrated in concentrated_lognormal_sd(). The row predicate matches
  // the likelihood's own -- a fitted year, a fishery, a positive observation --
  // so the sd is estimated over exactly the rows it is applied to. A zero-catch
  // year is legal input and log(0) is -inf, which would poison the fleet's sd,
  // hence the positive-observation guard.
  //
  // It reads catch_obs, not obsvec, so oneStepPredict()'s candidate value never
  // reaches it: the sd enters an OSA sweep as a fixed full-data number, as an
  // estimated sd held at its MLE would.
  catch_n_obs.setZero();
  catch_analytical_sd.setZero();
  for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){

    flt = catch_ctl(fsh_ind, 0) - 1;
    flt_yr = catch_ctl(fsh_ind, 2);

    if((flt_yr > 0) && (flt_yr <= endyr) && (flt_type(flt) == 1) && (catch_obs(fsh_ind, 0) > 0)){
      catch_n_obs(flt) += 1;
      catch_analytical_sd(flt) += square( log(catch_obs(fsh_ind, 0)) - log(catch_hat(fsh_ind)) );
    }
  }

  for(flt = 0; flt < n_flt; flt++){
    if(catch_n_obs(flt) > 0){
      catch_analytical_sd(flt) = concentrated_lognormal_sd(catch_analytical_sd(flt) / catch_n_obs(flt), bias_adjust_obs);
    }
  }


  // 9.2 Exploitable biomass ----
  exploitable_biomass.setZero();

  for(flt = 0; flt < n_flt; flt++) {

    sp = flt_spp(flt);
    wt_idx_flt = nspp * 2 + flt;

    if(flt_type(flt) == 1){ //Fishery only
      for(yr = 0; yr < nyrs; yr++) {
        for(age = 0; age < nages(sp); age++) {
          for(sex = 0; sex < nsex(sp); sex ++){
            exploitable_biomass(sp, yr) += N_at_age(sp, sex, age, yr) * weight_hat( wt_idx_flt, sex, age, yr ) * sel_at_age(flt, sex, age, yr) * proj_F_prop(flt);
          }
        }
      }
    }
  }


  /** ------------------------------------------------------------------------ //
   * 10. COMPOSITION EQUATIONS                                                  //
   * ------------------------------------------------------------------------- */

  // -- 10.1. Marginal age/length composition
  age_hat.setZero();     // Age-composition
  age_obs_hat.setZero(); // Age-composition with ageing error
  comp_hat.setZero();    // Age- or length-comp following data
  for(comp_ind = 0; comp_ind < comp_hat.rows(); comp_ind++){

    flt = comp_ctl(comp_ind, 0) - 1;            // Fleet index (0-based)
    sp = comp_ctl(comp_ind, 1) - 1;             // Species
    flt_sex = comp_ctl(comp_ind, 2);            // Comp sex: 0 = combined, 1 = female, 2 = male
    comp_type = comp_ctl(comp_ind, 3);          // Comp type: 0 = age, 1 = length, 2 = CAAL
    yr = comp_ctl(comp_ind, 4);                 // Data year
    mo = (growth_model(sp) > 0) ? flt_month(flt) : comp_n(comp_ind, 0); // Month (use fleet month if estimating growth, otherwise use from data)
    int wtind = nspp * 2 + flt;
    n_hat(comp_ind) = 0.0;

    // Likelihood (yr > 0) vs prediction (yr < 0)
    if(yr > 0){
      yr = yr - styr;
    }
    else if(yr < 0){
      yr = -yr - styr;
    }

    // Hindcast
    Type Frate = 0.0;
    if(yr < nyrs_hind){
      yr_ind = yr;
      Frate = exp(log_F(flt, yr));
    }

    // Projection
    if(yr >= nyrs_hind){
      yr_ind = nyrs_hind - 1;
      Frate = proj_F_prop(flt) * proj_F(sp, yr);
    }

    // 10.1.1. Estimated growth
    if(growth_model(sp) > 0){
      for(age = 0; age < nages(sp); age++) {
        for(ln = 0; ln < nlengths(sp); ln++) {
          for(sex = 0; sex < nsex(sp); sex ++){

            switch(flt_type(flt)){
            case 1: // - Fishery
              if(flt_sel_dim(flt) == 1){ // Length based
                {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (Frate / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr));
                }
              }
              break;


            case 2: // - Survey
              if(flt_sel_dim(flt) == 1){ // Length based
                {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (N_at_age(sp, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr)));
                }
              }
              break;
            }
          }
        }
      }

      // Assign CAAL to marginal age or length composition
      for(age = 0; age < nages(sp); age++) {
        for(ln = 0; ln < nlengths(sp); ln++) {

          // Age composition
          // - Catch-at-age
          if(comp_type == 0){
            switch(flt_sex){
            case 0: // Sexes combined or 1 sex assessment
              for(sex = 0; sex < nsex(sp); sex ++){
                age_hat(comp_ind, age ) += pred_CAAL(flt, sex, age, ln, yr);
              }
              break;

            case 1: case 2: // Sex-specific composition data
              sex = flt_sex - 1;
              age_hat(comp_ind, age ) += pred_CAAL(flt, sex, age, ln, yr);
              break;

            case 3: // Joint composition data
              for(sex = 0; sex < nsex(sp); sex ++){
                age_hat(comp_ind, age + nages(sp) * sex ) += pred_CAAL(flt, sex, age, ln, yr);
              }
              break;
            }
          }


          // Length composition
          // - Catch-at-length
          if(comp_type == 1){
            switch(flt_sex){
            case 0: // Sexes combined or 1 sex assessment
              for(sex = 0; sex < nsex(sp); sex ++){
                comp_hat(comp_ind, ln ) += pred_CAAL(flt, sex, age, ln, yr);
              }
              break;

            case 1: case 2: // Sex-specific composition data
              sex = flt_sex - 1;
              comp_hat(comp_ind, ln ) += pred_CAAL(flt, sex, age, ln, yr);
              break;

            case 3: // Joint composition data
              for(sex = 0; sex < nsex(sp); sex ++){
                comp_hat(comp_ind, ln + nlengths(sp) * sex ) += pred_CAAL(flt, sex, age, ln, yr) ;
              }
              break;
            }
          }
        }
      }

      // Adjust for aging error
      for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
        for(int true_age = 0; true_age < nages(sp); true_age++) {
          age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age, obs_age);
        }
      }

      // Adjust for aging error for joint data
      if(flt_sex == 3){
        for(int obs_age = nages(sp); obs_age < nages(sp) * 2; obs_age++) {
          for(int true_age = nages(sp); true_age < nages(sp) * 2; true_age++) {

            // Adjust indexing for joint age/length comp
            int true_age_tmp = true_age - nages(sp);
            int obs_age_tmp = obs_age - nages(sp);

            age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age_tmp, obs_age_tmp);
          }
        }
      }

      // Assign age-comp to estimated comp
      if( comp_type == 0) {
        joint_adjust(comp_ind) = 1;
        if(flt_sex == 3){
          joint_adjust(comp_ind) = 2;
        }
        for(age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
          comp_hat(comp_ind, age ) = age_obs_hat(comp_ind, age );
        }
      }
    } // End estimated growth loop


    // 10.1.2. Empirical weight-at-age
    if(growth_model(sp) == 0){
      for(age = 0; age < nages(sp); age++) {

        switch(flt_type(flt)){
        case 1: // - Fishery
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              // Catch-at-age
              age_hat(comp_ind, age ) += F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            // Catch-at-age
            age_hat(comp_ind, age ) = F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = F_flt_age(flt, sex, age, yr) / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr);
            }
            break;
          }
          break;


        case 2: // - Survey
          switch(flt_sex){
          case 0: // Sexes combined or 1 sex assessment
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age ) += N_at_age(sp, sex, age, yr)  * sel_at_age(flt, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            }
            break;

          case 1: case 2: // Sex-specific composition data
            sex = flt_sex - 1;
            // Survey catch-at-age
            age_hat(comp_ind, age ) = N_at_age(sp, sex, age, yr)  * sel_at_age(flt, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            break;

          case 3: // Joint composition data
            for(sex = 0; sex < nsex(sp); sex ++){
              // Survey catch-at-age
              age_hat(comp_ind, age + nages(sp) * sex ) = N_at_age(sp, sex, age, yr)  * sel_at_age(flt, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr));
            }
            break;
          }
        }
      }


      // Adjust for aging error
      for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
        for(int true_age = 0; true_age < nages(sp); true_age++) {
          age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age, obs_age);
        }
      }

      // Adjust for aging error for joint data
      if(flt_sex == 3){
        for(int obs_age = nages(sp); obs_age < nages(sp) * 2; obs_age++) {
          for(int true_age = nages(sp); true_age < nages(sp) * 2; true_age++) {

            // Adjust indexing for joint age/length comp
            int true_age_tmp = true_age - nages(sp);
            int obs_age_tmp = obs_age - nages(sp);

            age_obs_hat(comp_ind, obs_age) += age_hat(comp_ind, true_age ) * age_error(flt_ageing_error_index(flt), true_age_tmp, obs_age_tmp);
          }
        }
      }

      // Assign age-comp to estimated comp
      if( comp_type == 0) {
        joint_adjust(comp_ind) = 1;
        if(flt_sex == 3){
          joint_adjust(comp_ind) = 2;
        }
        for(age = 0; age < nages(sp) * joint_adjust(comp_ind); age++) {
          comp_hat(comp_ind, age ) = age_obs_hat(comp_ind, age );
        }
      }

      // Convert from catch-at-age to catch-at-length, reading age_hat, the TRUE ages, not
      // age_obs_hat: the age-length key is P(length | true age), and a measured length
      // carries no otolith reading, so ageing error has no part in a length composition.
      // SS3 likewise applies its age_age matrix only to age and CAAL data.
      if( comp_type == 1) {
        for(ln = 0; ln < nlengths(sp); ln++) {
          for(age = 0; age < nages(sp); age++) {

            sex = 0;

            // Adjust sex for males/females
            if((flt_sex > 0) & (flt_sex < 3)){
              sex = flt_sex - 1;
            }
            comp_hat(comp_ind, ln ) += age_hat(comp_ind, age) * age_trans_matrix(flt_age_transition_index(flt), sex, age, ln );
          }
        }


        // Convert from catch-at-age to catch-at-length for joint comp data
        if( flt_sex == 3) {
          for(ln = nlengths(sp); ln < nlengths(sp) * 2; ln++) {
            for(age = nages(sp); age < nages(sp) * 2; age++) {

              // Adjust indexing for joint age/length comp
              int obs_log_tmp = ln;
              int obs_age_tmp = age;

              obs_age_tmp = age - nages(sp);
              obs_log_tmp = ln - nlengths(sp);
              sex = 1;

              comp_hat(comp_ind, ln ) += age_hat(comp_ind, age) * age_trans_matrix(flt_age_transition_index(flt), sex, obs_age_tmp, obs_log_tmp );
            }
          }
        }
      }
    } // End empirical weight-at-age loop

    // Standardize to sum to 1
    Type penalty = 0.0;
    comp_hat.row(comp_ind) /= posfun(comp_hat.row(comp_ind).sum(), Type(0.001), penalty);
    (void) penalty;
  }


  // -- 10.2. CAAL
  // -- length-based selectivity was converted above
  caal_hat.setZero();
  for(int caal_ind = 0; caal_ind < caal_hat.rows(); caal_ind++){

    flt = caal_ctl(caal_ind, 0) - 1;            // Fleet index (0-based)
    sp = caal_ctl(caal_ind, 1) - 1;             // Species
    flt_sex = caal_ctl(caal_ind, 2);            // CAAL sex: 0 = combined, 1 = female, 2 = male
    yr = caal_ctl(caal_ind, 3);                 // Data year
    ln = caal_ctl(caal_ind, 4) - 1;             // Length bin
    mo = flt_month(flt);                        // Month (use fleet month b/c estimating growth)

    if(growth_model(sp) == 0) continue; // Skip empirical weight-at-age species

    if(flt_sex > 0){
      sex = flt_sex - 1;
    } else {
      sex = flt_sex;
    }

    int wtind = nspp * 2 + flt;

    // Likelihood (yr > 0) vs prediction (yr < 0)
    if(yr > 0){
      yr = yr - styr;
    }
    else if(yr < 0){
      yr = -yr - styr;
    }
    // Hindcast
    Type Frate = 0.0;
    if(yr < nyrs_hind){
      yr_ind = yr;
      Frate = exp(log_F(flt, yr));
    }

    // Projection
    if(yr >= nyrs_hind){
      yr_ind = nyrs_hind - 1;
      Frate = proj_F_prop(flt) * proj_F(sp, yr);
    }

    // 10.2. Calculate CAAL (predicted conditional age-at-length)
    for(age = 0; age < nages(sp); age++) {

      switch(flt_type(flt)){
      case 1: // - Fishery
        if(flt_sel_dim(flt) == 1){
          {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (Frate / Z_at_age(sp, sex, age, yr) * (1 - exp(-Z_at_age(sp, sex, age, yr))) * N_at_age(sp, sex, age, yr));
                }
        }
        break;


      case 2: // - Survey
        if(flt_sel_dim(flt) == 1){
          {
                  // sel(L) * P(L|age) over this data bin's population bins (2.3c).
                  Type ps_ln = 0.0;
                  for(int lp = pop_bin_lo(sp, ln); lp <= pop_bin_hi(sp, ln); lp++){
                    ps_ln += sel_at_length(flt, sex, lp, yr) * growth_matrix(wtind, sex, age, lp, yr);
                  }
                  pred_CAAL(flt, sex, age, ln, yr) = ps_ln * (N_at_age(sp, sex, age, yr) * index_q(flt, yr_ind) * exp( - Type(mo/12.0) * Z_at_age(sp, sex, age, yr)) );
                }
        }
        break;
      }
    }

    // Adjust for aging error
    for(int obs_age = 0; obs_age < nages(sp); obs_age++) {
      for(int true_age = 0; true_age < nages(sp); true_age++) {
        caal_hat(caal_ind, obs_age) += pred_CAAL(flt, sex, true_age, ln, yr) * age_error(flt_ageing_error_index(flt), true_age, obs_age);
      }
    }

    //  Observed CAAL standardize to sum to 1
    Type penalty = 0.0;
    caal_hat.row(caal_ind) /= posfun(caal_hat.row(caal_ind).sum(), Type(0.0001), penalty);
    (void) penalty;
  }



  /** ------------------------------------------------------------------------ //
   * 11. PREDICTED STOMACH CONTENT                                             //
   * ------------------------------------------------------------------------- */

  // Predict stomach content
  if((msmMode > 2) || (imax(suitMode) > 0)) {
    predict_stomach_content(nspp, nyrs_hind, styr, suit_styr, suit_endyr,
                            minage, nsex, nages, diet_obs, diet_ctl,
                            avgN_at_age, diet_prop_hat, diet_hat);
  }


  /** ------------------------------------------------------------------------ //
   // 12. DERIVED QUANTITIES
   * ------------------------------------------------------------------------- */


  // 12.1. Depletion
  biomass_depletion.setZero();
  ssb_depletion.setZero();
  for(sp = 0; sp < nspp; sp++){
    for(yr = 0; yr < nyrs; yr++){

      if(DynamicHCR == 0){
        biomass_depletion(sp, yr) = biomass(sp, yr)/B0(sp, nyrs-1);
        ssb_depletion(sp, yr) = ssb(sp, yr)/SB0(sp, nyrs-1);
      }

      if(DynamicHCR == 1){
        biomass_depletion(sp, yr) = biomass(sp, yr)/DynamicB0(sp, yr);
        ssb_depletion(sp, yr) = ssb(sp, yr)/DynamicSB0(sp, yr);
      }

      // Multi-species, no HCR: the projection runs unfished, so SSB in the last
      // projection year is multispecies SB0 once `projyr` is far enough out to
      // equilibrate. Excluded under DynamicHCR, where DynamicSB0 is the
      // reference and this would overwrite it.
      if((HCR == 0) & (msmMode > 0) & (DynamicHCR == 0)){
        biomass_depletion(sp, yr) = biomass(sp, yr)/biomass(sp, nyrs-1);
        ssb_depletion(sp, yr) = ssb(sp, yr)/ssb(sp, nyrs-1);
      }
    }
  }




  /** ------------------------------------------------------------------------ //
   // 13. LIKELIHOOD EQUATIONS
   * ------------------------------------------------------------------------- */

  // 13.0. OBJECTIVE FUNCTION
  int n_col = std::max(n_flt, nspp);
  // Named rows of jnll_comp / unweighted_jnll_comp. These codes MUST stay in
  // sync with the row labels in R/6-rename_output.R (rownames of
  // quantities$jnll_comp): add or reorder a row here and there in lockstep.
  enum JnllRow {
    JNLL_INDEX          = 0,   // Index data
    JNLL_CATCH          = 1,   // Catch data
    JNLL_COMP           = 2,   // Composition data
    JNLL_CAAL           = 3,   // CAAL data
    JNLL_SEL_NONPARAM   = 4,   // Non-parametric selectivity
    JNLL_SEL_DEV        = 5,   // Selectivity deviates
    JNLL_Q_PRIOR        = 6,   // Catchability prior
    JNLL_Q_DEV          = 7,   // Catchability deviates
    JNLL_SRR_PRIOR      = 8,   // Stock-recruit prior
    JNLL_INIT_DEV       = 9,   // Initial abundance deviates
    JNLL_REC_DEV        = 10,  // Recruitment deviates
    JNLL_SRR_PENALTY    = 11,  // Stock-recruit penalty
    JNLL_REFPT_PENALTY  = 12,  // Reference point penalties
    JNLL_ZERO_N_PENALTY = 13,  // Zero n-at-age penalty
    JNLL_M_PRIOR        = 14,  // M prior
    JNLL_M_RE           = 15,  // M random effects
    JNLL_RATION         = 16,  // Ration
    JNLL_RATION_PENALTY = 17,  // Ration penalties
    JNLL_STOMACH        = 18,  // Stomach content data
    JNLL_LINKAGE_PRIOR  = 19,  // Linkage-table priors (per-row)
    JNLL_LINKAGE_RE     = 20,  // Linkage random effects
    JNLL_EQUIL_CATCH    = 21,  // Initial equilibrium catch
    JNLL_DSEM           = 22,  // DSEM latent-state GMRF (recruitment)
    JNLL_N_ROWS         = 23   // total row count (for dimensioning)
  };
  matrix<Type> jnll_comp(JNLL_N_ROWS, n_col); jnll_comp.setZero();  // negative log-likelihood components
  matrix<Type> unweighted_jnll_comp(JNLL_N_ROWS, n_col); unweighted_jnll_comp.setZero();  // same, without likelihood weights

  // What a COLUMN means changes from row to row: fleets on rows 0-7 (data,
  // selectivity, catchability), species on rows 8-19, neither on row 20 (one
  // model-wide linkage density), fleets again on row 21, and row 22 is mixed --
  // the recruitment-SD prior writes to each species' column, the joint GMRF
  // density to column 0. Column j is therefore a fleet on one row and a species
  // on the next, so summing down a column mixes the two. .JNLL_ROW_AXIS in
  // R/9-profile.R is the registry of that axis.


  // 13.1. FIT OBJECTIVE FUNCTION
  // Slot 0 -- Index of abundance biomass/numbers --
  Type index_std_dev = 0;
  for(index_ind = 0; index_ind < index_obs.rows(); index_ind++) {

    index = index_ctl(index_ind, 0) - 1;            // Survey index
    sp = index_ctl(index_ind, 1) - 1;             // Species is the 2nd column
    flt_yr = index_ctl(index_ind, 2);             // Data year

    // Set up variance
    switch (est_sigma_index(index)) {
    case 0:     // Provided standard deviation
      index_std_dev = index_obs(index_ind, 1);
      break;
    case 1:     // Estimated standard deviation
      index_std_dev = exp(index_log_sd(index));
      break;
    case 2:     // Analytical (Walters and Ludwig 1994); see section 8.4
      index_std_dev = index_analytical_sd(index);
      break;
    default:
      error("Invalid 'Estimate_sigma_index'");
    }

    index_sd(index_ind) = index_std_dev;

    // Only include years from hindcast. Lognormal IID branch (index_ll_type == 0);
    // MVN covariance fleets (index_ll_type == 1) are handled in the per-fleet loop
    // below and contribute nothing here.
    if((flt_yr > 0) && (flt_yr <= endyr) && (flt_type(index) > 0) && (index_ll_type(index) == 0)){
      if(index_obs(index_ind) > 0){
        // Read the (log) observation from obsvec and gate it with keep so
        // oneStepPredict() can compute OSA residuals. With keep == 1 (ordinary
        // fitting) and obsvec(pos) == log(index_obs) this is the lognormal
        // likelihood exactly. The `pos >= 0` guard matches the comp/caal/diet
        // and MVN/Normal OSA branches; build_osa_data()'s inclusion set is this
        // same predicate, so it guards only against a future divergence.
        int pos = index_obsvec_idx(index_ind);
        if(pos >= 0){
          Type index_mu = log(index_hat(index_ind)) - bias_adjust_obs*square(index_std_dev)/2.0;
          jnll_comp(JNLL_INDEX, index) -= keep(pos) * dnorm(obsvec(pos), index_mu, index_std_dev, true);
          // method = "cdf": the same lognormal read as a CDF on the log scale.
          if(osa_mode == 2){
            jnll_comp(JNLL_INDEX, index) -= osa_norm_cdf_terms(
              obsvec(pos), index_mu, index_std_dev, keep.cdf_lower(pos), keep.cdf_upper(pos));
          }
        }
      }
    }

    // Natural-scale normal: the residual (obs - q*pred) is normal with an ABSOLUTE sd in
    // the units of the index, not a log-scale CV, so there is no lognormal bias term. Two
    // families share this branch, differing only in whether the density is renormalized
    // over the positive half line:
    //
    //   3  "Normal"           plain normal on (-inf, inf). Matches the AMAK
    //                         avo_like/cpue_like = 0.5*(obs - q*pred)^2 / sd^2 term for
    //                         term, which is what the ADMB bridges compare against.
    //   4  "TruncatedNormal"  normal left-truncated at zero. An index cannot be negative
    //                         and data_check() will not accept one, so the support is
    //                         (0, inf) and the density is renormalized:
    //                         log f = log phi(x; mu, sd) - log(Phi(mu/sd)).
    //
    // A modelling choice, not a correction: only family 4 draws from the same
    // distribution the likelihood scores (see the draw below). Family 3's draw is exact
    // for the untruncated normal it is fitted under but must reject the non-positive
    // draws data_check() would refuse, so draw and density differ where the absolute sd
    // is large relative to the index.
    //
    // mu = q * (selected biomass) >= 0 and sd > 0, so mu/sd >= 0 and the truncation
    // constant lies in [log(0.5), 0]: bounded, no underflow.
    if((flt_yr > 0) && (flt_yr <= endyr) && (flt_type(index) > 0) &&
       ((index_ll_type(index) == 3) || (index_ll_type(index) == 4))){
      if(index_obs(index_ind, 0) > 0){
        Type log_Z = (index_ll_type(index) == 4) ? log(pnorm(index_hat(index_ind) / index_std_dev)) : Type(0.0);
        if(osa_mode == 0){
          jnll_comp(JNLL_INDEX, index) -= dnorm(index_obs(index_ind, 0), index_hat(index_ind), index_std_dev, true);
          jnll_comp(JNLL_INDEX, index) += log_Z;
        } else {
          // OSA: the natural-scale observation from obsvec (build_osa_data()
          // stores it untransformed for both families), keep-gated, so
          // oneStepPredict() residualizes against the density that was fit.
          int pos = index_obsvec_idx(index_ind);
          if(pos >= 0){
            jnll_comp(JNLL_INDEX, index) -= keep(pos) * dnorm(obsvec(pos), index_hat(index_ind), index_std_dev, true);
            jnll_comp(JNLL_INDEX, index) += keep(pos) * log_Z;
            // method = "cdf". "TruncatedNormal" gets its own CDF, renormalized over (0,
            // inf) -- the correction the Gaussian methods cannot see -- so it needs no
            // oneStepGeneric call and none of its integration failures.
            if(osa_mode == 2){
              if(index_ll_type(index) == 4){
                jnll_comp(JNLL_INDEX, index) -= osa_trunc_norm_cdf_terms(
                  obsvec(pos), index_hat(index_ind), index_std_dev,
                  keep.cdf_lower(pos), keep.cdf_upper(pos));
              } else {
                jnll_comp(JNLL_INDEX, index) -= osa_norm_cdf_terms(
                  obsvec(pos), index_hat(index_ind), index_std_dev,
                  keep.cdf_lower(pos), keep.cdf_upper(pos));
              }
            }
          }
        }
      }
    }

    // -- Simulate the survey observation, for sim_mod(simulate = TRUE) --
    // Families 0 (lognormal), 3 (normal) and 4 (truncated normal), each drawing
    // what its own density assumes; the correlated families are drawn per fleet
    // below, where their covariance lives.
    //
    // Outside the hindcast gate, deliberately: run_mse() splices the
    // negative-Year rows in sign-flipped as the next assessment's data, so a
    // gated draw would leave the projection carrying unsimulated observations.
    SIMULATE {
      if((flt_type(index) > 0) && (index_ll_type(index) == 0)){
        index_obs(index_ind, 0) = exp(rnorm(log(index_hat(index_ind)) - bias_adjust_obs*square(index_std_dev)/2.0, index_std_dev));
        // The density reads obsvec for this family; keep the two in step.
        int sim_pos = index_obsvec_idx(index_ind);
        if(sim_pos >= 0){
          obsvec(sim_pos) = log(index_obs(index_ind, 0));
        }
      }
      if((flt_type(index) > 0) && (index_ll_type(index) == 3)){
        // "Normal": the density is untruncated, so the exact draw is a plain
        // normal, which can come back non-positive for data_check() to reject.
        // Redraw, and count the rejections: a redrawn row follows the normal
        // truncated at zero, not the untruncated one the likelihood scores, and
        // sim_mod() warns when that is doing enough of the work to matter. Use
        // "TruncatedNormal" to remove the gap rather than measure it.
        Type sim_draw = rnorm(index_hat(index_ind), index_std_dev);
        int  tn_tries = 0;
        while((sim_draw <= Type(0)) && (tn_tries < index_trunc_max_tries)){
          index_trunc_rejects_sim(index_ind) += 1;
          sim_draw = rnorm(index_hat(index_ind), index_std_dev);
          tn_tries++;
        }
        index_trunc_tries_sim(index_ind) += tn_tries + 1;
        index_trunc_sd_sim(index_ind) = index_std_dev;
        index_trunc_budget_sim(index_ind) = index_trunc_max_tries;
        index_obs(index_ind, 0) = sim_draw;
      }
      if((flt_type(index) > 0) && (index_ll_type(index) == 4)){
        // "TruncatedNormal": left-truncated at zero by inverse CDF, matching the
        // density above. With a = Phi(-mu/sd), drawing u uniformly on (a, 1) and
        // returning mu + sd*Phi^-1(u) is exact for the normal restricted to
        // (0, inf): no rejection loop, no retry budget, no non-positive draw.
        Type a = pnorm(-index_hat(index_ind) / index_std_dev);
        Type u = a + (Type(1.0) - a) * runif(Type(0.0), Type(1.0));
        index_obs(index_ind, 0) = index_hat(index_ind) + index_std_dev * qnorm(u);
      }
      // Keep obsvec in step for both natural-scale families, in whichever layout
      // it is carrying: build_osa_data() stores the UNTRANSFORMED observation
      // for these families, which the OSA branch above reads, while the ordinary
      // fitting path keeps the log(obs) layout used for every fleet. Writing the
      // wrong one would leave obsvec_sim describing different data from
      // index_obs_sim. Same reasoning as the MVN block below.
      if((flt_type(index) > 0) && ((index_ll_type(index) == 3) || (index_ll_type(index) == 4))){
        int sim_pos = index_obsvec_idx(index_ind);
        if(sim_pos >= 0){
          if(osa_mode == 0){
            if(index_obs(index_ind, 0) > Type(0)) obsvec(sim_pos) = log(index_obs(index_ind, 0));
          } else {
            obsvec(sim_pos) = index_obs(index_ind, 0);
          }
        }
      }
    }
  }

  // MVN survey biomass likelihood (Index_loglike == "MVN"): the AMAK/ebswp DoCovBTS
  // covariance survey likelihood 0.5 * r' Sigma^-1 r on each MVN fleet's fitted residual
  // vector r = obs - q*pred (arithmetic, natural scale; pair with est_index_q == 7 for
  // the AMAK arithmetic-mean q). Sigma is the per-fleet index_cov_mat(index), and r is
  // assembled in index_obs row order to match its rows and columns. Reuses jnll_comp
  // row 0.
  //
  // OSA (osa_mode == 1): the correlated block is whitened into independent standard
  // normals so oneStepPredict() can residualize each observation. With the lower Cholesky
  // Sigma = L L', build_osa_data() stores z = L^-1 obs in obsvec and the whitened
  // mean m = L^-1 index_hat is formed here with the SAME L; z_k ~ N(m_k, 1) are
  // independent and z_k - m_k is the multivariate-Gaussian one-step-ahead residual
  // (Thygesen et al. 2017, the SAM/TMB construction). Sigma is fixed data, so L is
  // constant and the whitened density differs from the fitting density only by
  // 0.5*logdet(Sigma) (family 1 also drops n/2*log(2*pi)), not moving the residuals.
  for(index = 0; index < n_flt; index++){
    if((flt_type(index) > 0) && (index_ll_type(index) == 1 || index_ll_type(index) == 2)){   // 1 = MVN (bare quadratic form), 2 = MVNORM (full density) ONLY -- other families (0 lognormal, 3 normal) carry a 1x1 dummy Sigma and must not enter here
      int n_mvn = 0;
      for(index_ind = 0; index_ind < index_obs.rows(); index_ind++){
        if((index_ctl(index_ind, 0) - 1) == index){
          flt_yr = index_ctl(index_ind, 2);
          if((flt_yr > 0) && (flt_yr <= endyr) && (index_obs(index_ind, 0) > 0)) n_mvn++;
        }
      }
      if(n_mvn > 0){
        vector<Type> resid(n_mvn);
        vector<Type> mu(n_mvn);
        vector<int>  posv(n_mvn);
        vector<int>  rowv(n_mvn);   // index_obs row per Sigma position, for SIMULATE
        int k = 0;
        for(index_ind = 0; index_ind < index_obs.rows(); index_ind++){
          if((index_ctl(index_ind, 0) - 1) == index){
            flt_yr = index_ctl(index_ind, 2);
            if((flt_yr > 0) && (flt_yr <= endyr) && (index_obs(index_ind, 0) > 0)){
              resid(k) = index_obs(index_ind, 0) - index_hat(index_ind);  // arithmetic residual (obs - q*pred)
              mu(k)    = index_hat(index_ind);
              posv(k)  = index_obsvec_idx(index_ind);
              rowv(k)  = index_ind;
              k++;
            }
          }
        }
        if(osa_mode == 0){
          // TMB-native multivariate normal: density::MVNORM(Sigma) factorizes Sigma
          // internally (robust; no explicit inverse). MVNORM(Sigma)(r) returns the full
          // negative log-density 0.5*(r' Sigma^-1 r + logdet(Sigma) + n*log(2*pi)).
          Type dens = MVNORM(index_cov_mat(index))(resid);
          // "MVN" (1) reports the bare quadratic form 0.5 r' Sigma^-1 r (the AMAK/ebswp
          // value) by removing the fixed normalizing constant; "MVNORM" (2) keeps the
          // full density. Both give an identical fit (the constant has zero gradient).
          jnll_comp(JNLL_INDEX, index) += (index_ll_type(index) == 2) ? dens : (dens - index_cov_const(index));
        } else {
          // OSA: whiten with the lower Cholesky of Sigma (same unique factor
          // build_osa_data() uses to whiten the observations), then score each
          // whitened observation as an independent standard normal.
          Eigen::LLT<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > llt(index_cov_mat(index));
          matrix<Type> Lmat = llt.matrixL();
          matrix<Type> muc(n_mvn, 1);
          for(k = 0; k < n_mvn; k++) muc(k, 0) = mu(k);
          matrix<Type> mwhite = Lmat.template triangularView<Eigen::Lower>().solve(muc);  // m = L^-1 mu
          for(k = 0; k < n_mvn; k++){
            int pos = posv(k);
            if(pos >= 0){
              jnll_comp(JNLL_INDEX, index) -= keep(pos) * dnorm(obsvec(pos), mwhite(k, 0), Type(1.0), true);
              // method = "cdf": the whitened block is independent standard
              // normals, so each row's CDF is pnorm of its own innovation.
              if(osa_mode == 2){
                jnll_comp(JNLL_INDEX, index) -= osa_norm_cdf_terms(
                  obsvec(pos), mwhite(k, 0), Type(1.0),
                  keep.cdf_lower(pos), keep.cdf_upper(pos));
              }
            }
          }
        }

        // -- Simulate the correlated survey block --
        // MVN and MVNORM differ only by a constant in the density, so they
        // simulate identically; MVNORM_t draws mean-zero, hence adding mu back.
        // Only the FITTED rows are drawn -- Sigma is dimensioned to exactly
        // those, in the assembly order above -- so unlike the independent
        // families, rows outside the fitted set keep the values they came with.
        SIMULATE {
          // Redraw while any row is non-positive. The target is the joint normal
          // truncated to the positive orthant, which has no closed-form inverse CDF, so
          // rejection is the exact sampler here, not a fallback; truncating each margin
          // separately would break the correlation this likelihood models. There is no
          // truncated family to switch to, so the draw/density gap stays. The budget
          // scales with the row count (index_trunc_max_tries * n_mvn) because a vector
          // is rejected if ANY row is non-positive, and the acceptance test sits at the
          // TOP of the loop, so an exhausted budget writes a non-positive index that has
          // been examined and sim_mod() reports it as any other.
          vector<Type> sim_dev = MVNORM(index_cov_mat(index)).simulate();
          int mvn_tries = 0;
          while(true){
            bool any_nonpos = false;
            for(k = 0; k < n_mvn; k++){
              if(mu(k) + sim_dev(k) <= Type(0)) any_nonpos = true;
            }
            if(!any_nonpos) break;
            if(mvn_tries >= index_trunc_max_tries * n_mvn) break;   // budget spent; the non-positive draw stands and is reported
            for(k = 0; k < n_mvn; k++){
              if(mu(k) + sim_dev(k) <= Type(0)) index_trunc_rejects_sim(rowv(k)) += 1;
            }
            sim_dev = MVNORM(index_cov_mat(index)).simulate();
            mvn_tries++;
          }
          for(k = 0; k < n_mvn; k++){
            index_trunc_tries_sim(rowv(k)) += mvn_tries + 1;
            index_trunc_sd_sim(rowv(k)) = sqrt(index_cov_mat(index)(k, k));
            index_trunc_budget_sim(rowv(k)) = index_trunc_max_tries * n_mvn;
            index_obs(rowv(k), 0) = mu(k) + sim_dev(k);
            // obsvec holds log(obs) for this family on the ordinary fitting path
            // (build_osa_data() only whitens it when the OSA data are built), so
            // keep it in step on that path and leave the whitened layout alone.
            int sim_pos = posv(k);
            if((sim_pos >= 0) && (osa_mode == 0) && (index_obs(rowv(k), 0) > Type(0))){
              obsvec(sim_pos) = log(index_obs(rowv(k), 0));
            }
          }
        }
      }
    }
  }

  // Reported under a _sim name -- see the naming rule in section 5.2.
  SIMULATE {
    matrix<Type> index_obs_sim = index_obs;
    REPORT(index_obs_sim);
    REPORT(index_trunc_tries_sim);
    REPORT(index_trunc_rejects_sim);
    REPORT(index_trunc_sd_sim);
    REPORT(index_trunc_budget_sim);
  }


  // Slot 1 -- Total catch --
  Type fsh_std_dev = 0;
  for(fsh_ind = 0; fsh_ind < catch_obs.rows(); fsh_ind++) {

    flt = catch_ctl(fsh_ind, 0) - 1;            // Fleet index (0-based)
    sp = catch_ctl(fsh_ind, 1) - 1;             // Species is the column 3
    flt_yr = catch_ctl(fsh_ind, 2);             // Data year
    yr = flt_yr - styr;                         // Year index, 0-based

    // Set up variance
    switch (est_sigma_fsh(flt)) {
    case 0: // Provided standard deviation
      fsh_std_dev = catch_obs(fsh_ind, 1);
      break;
    case 1:     // Estimated standard deviation
      fsh_std_dev = exp(catch_log_sd(flt));
      break;
    case 2:     // Analytical (Walters and Ludwig 1994); see section 9.1b
      fsh_std_dev = catch_analytical_sd(flt);
      break;
    default:
      error("Invalid 'Estimate_sigma_catch'");
    }

    catch_sd(fsh_ind) = fsh_std_dev; // the sd this row was actually fitted with

    // Add only years from hindcast
    if((flt_yr > 0) && (flt_yr <= endyr) && (flt_type(flt) == 1)){
      if(catch_obs(fsh_ind, 0) > 0){
        // The (log) observation from obsvec, keep-gated for OSA residuals (see
        // the index slot above). With keep == 1 and obsvec(pos) == log(catch_obs)
        // this is the lognormal likelihood exactly; the `pos >= 0` guard is the
        // same defensive parity as the index/comp/caal/diet reads.
        int pos = catch_obsvec_idx(fsh_ind);
        if(pos >= 0){
          Type catch_mu = log(catch_hat(fsh_ind)) - bias_adjust_obs*square(fsh_std_dev)/2.0;
          jnll_comp(JNLL_CATCH, flt) -= keep(pos) * dnorm(obsvec(pos), catch_mu, fsh_std_dev, true) ;
          // method = "cdf": the same lognormal read as a CDF on the log scale.
          if(osa_mode == 2){
            jnll_comp(JNLL_CATCH, flt) -= osa_norm_cdf_terms(
              obsvec(pos), catch_mu, fsh_std_dev, keep.cdf_lower(pos), keep.cdf_upper(pos));
          }
        }
      }
    }

    // -- Simulate the catch observation, for sim_mod(simulate = TRUE) --
    // Same mean as the dnorm above, bias term included: the estimator fits to that mean,
    // so a draw centred elsewhere biases every self-test.
    //
    // Outside the hindcast gate, deliberately. clean_data() appends a catch row per
    // fishery per projection year with Catch = NA (99 of BS2017SS's 216) that sim_mod()
    // fills; a gated draw leaves them NA and data_check() rejects the refit. catch_hat is
    // 0 there, so the mean is log(0) = -inf and rnorm() returns it undrawn, giving
    // exp(-inf) = 0 -- safe only because SIMULATE runs on doubles. Do not share this mean
    // with the dnorm above, where log(0) puts a NaN adjoint on the tape, and do not add a
    // catch_hat > 0 guard, which would change how many random numbers the draw consumes.
    SIMULATE {
      catch_obs(fsh_ind, 0) = exp(rnorm(log(catch_hat(fsh_ind)) - bias_adjust_obs*square(fsh_std_dev)/2.0, fsh_std_dev));
      // The density reads obsvec (log scale), the write-back catch_obs (natural
      // scale); keep both in step, the index being -1 outside the fitted window.
      // Both writes last only for this evaluation -- DATA_ objects are re-read
      // each call -- so a simulate call still returns the original data's jnll.
      // Fit simulated data by building a new model on them.
      int sim_pos = catch_obsvec_idx(fsh_ind);
      if(sim_pos >= 0){
        obsvec(sim_pos) = log(catch_obs(fsh_ind, 0));
      }
    }
  }

  // -- Initial equilibrium catch likelihood
  // The lognormal the hindcast catch uses, on the row's own Log_sd, booked to
  // its own jnll row as SS3 does (equ_catch_like). Its own loop, so it still
  // runs for a fleet with no hindcast catch rows. No OSA slot: one observation
  // has no residual sequence to condition on.
  for(int eq_ind = 0; eq_ind < equil_catch_ctl.rows(); eq_ind++){
    int eq_flt = equil_catch_ctl(eq_ind, 0) - 1;
    if((equil_catch_obs(eq_ind, 0) > 0) && (flt_type(eq_flt) == 1)){
      Type eq_sd = equil_catch_obs(eq_ind, 1);
      Type eq_mu = log(equil_catch_hat(eq_ind)) - bias_adjust_obs*square(eq_sd)/2.0;
      jnll_comp(JNLL_EQUIL_CATCH, eq_flt) -=
        dnorm(log(equil_catch_obs(eq_ind, 0)), eq_mu, eq_sd, true);

      // -- Simulate the equilibrium catch, for sim_mod(simulate = TRUE) --
      // Same mean as the dnorm above, bias term included. This observation is the ONLY
      // one informing Finit, and therefore the initial age-structure and the SSB scale
      // under initMode 6, so leaving it at its real value would condition every
      // self_test() replicate on data its own operating model had not generated.
      //
      // Inside the gate, unlike the hindcast catch: every equilibrium row is a real
      // observation data_check() requires to be positive, so there is no projection row
      // to fill.
      SIMULATE {
        equil_catch_obs(eq_ind, 0) = exp(rnorm(eq_mu, eq_sd));
      }
    }
  }

  // Reported under _sim names -- see the naming rule in section 5.2. obsvec_sim holds the
  // catch and survey-index entries redrawn and every other entry as observed, and no
  // more: the composition, CAAL and diet draws below write their own matrices and never
  // touch obsvec, their densities reading it only in OSA mode. It is for inspection, not
  // a simulated data set -- sim_mod() assembles the data_list from the *_obs_sim
  // matrices.
  SIMULATE {
    matrix<Type> catch_obs_sim = catch_obs;
    matrix<Type> equil_catch_obs_sim = equil_catch_obs;
    vector<Type> obsvec_sim = obsvec;
    REPORT(catch_obs_sim);
    REPORT(equil_catch_obs_sim);
    REPORT(obsvec_sim);
  }


  // Slot 2 -- Age/length composition
  // Small proportion offset added to obs and predicted comps before the
  // multinomial to keep log(0) out. Set from R through
  // rearrange_data()/fit_control(), default 1e-5, and the OSA obsvec is built
  // with the same offset so fitting and OSA residuals stay consistent.
  Type comp_prop_offset = comp_offset;
  // FIXME: case 0 fitting routes through dmultinom_osa(), which renormalizes p,
  // so the *reported* case-0 multinomial NLL carries a per-row additive constant
  // Neff*log(1 + n_comp*offset) that an un-renormalized multinomial does not.
  // Gradient and MLE are unchanged. Subtract the constant if the reported NLL
  // has to be comparable with a version fitted without it.
  //
  // Simulated compositions go to this copy, not into the REPORTed comp_obs --
  // see the _sim naming rule in section 5.2. Rows the draw does not touch keep
  // the observed values they were copied with.
  matrix<Type> comp_sim = comp_obs;

  for(comp_ind = 0; comp_ind < comp_obs.rows(); comp_ind++) {

    flt = comp_ctl(comp_ind, 0) - 1;        // Fleet index
    sp = comp_ctl(comp_ind, 1) - 1;         // Species
    flt_sex = comp_ctl(comp_ind, 2);        // Comp sex: 0 = combined, 1 = female, 2 = male, 3 = joint
    comp_type = comp_ctl(comp_ind, 3);      // Comp type: 0 = age, 1 = length
    yr = comp_ctl(comp_ind, 4);             // Data year

    // Adjustment for joint sex composition data
    joint_adjust(comp_ind) = 1;
    if(flt_sex == 3){
      joint_adjust(comp_ind) = 2;
    }

    // Number of ages
    int n_comp = 0;
    if(comp_type == 0){
      n_comp = nages(sp) * joint_adjust(comp_ind); // For sex = 3 (joint sex comp data)

    }
    if(comp_type == 1){
      n_comp = nlengths(sp) * joint_adjust(comp_ind);
    }

    // Select sections
    vector<Type> comp_obs_tmp = comp_obs.row(comp_ind).segment(0, n_comp); // Observed proportion
    vector<Type> comp_hat_tmp = comp_hat.row(comp_ind).segment(0, n_comp); // Expected proportion

    // Composition young/old accumulation (AFSC ac_yng / ac_old): fold the bins below
    // `yng` into the `yng` bin and above `old` into the `old` bin, then restrict the
    // likelihood to [yng, old]. Done PER SEX BLOCK for joint-sex comps
    // (joint_adjust == 2) so a fold never crosses the sex boundary. yng and old are
    // 1-based bin ordinals on the fleet's comp dimension (age or length), and the default
    // (yng = 1, old >= nbins) is a no-op. Applied on the OSA path too, since OSA
    // residuals must correspond to the folded model that was fit: the branch below reads
    // a folded obsvec, built identically in build_osa_data(), against this folded
    // comp_hat_tmp / n_comp.
    {
      int nbins_blk = (comp_type == 0) ? nages(sp) : nlengths(sp);
      int yng = comp_accum_young(flt);
      int old = comp_accum_old(flt);
      if (yng < 1) yng = 1;
      if (old < 1 || old > nbins_blk) old = nbins_blk;
      if (yng > old) yng = old;   // defensive: keep 1 <= yng <= old (data_check() validates)
      if (yng > 1 || old < nbins_blk) {
        int nkeep = old - yng + 1;
        int nblk  = joint_adjust(comp_ind);
        vector<Type> obs_acc(nblk * nkeep); obs_acc.setZero();
        vector<Type> hat_acc(nblk * nkeep); hat_acc.setZero();
        for (int b = 0; b < nblk; b++) {
          for (int j = 0; j < nbins_blk; j++) {
            int tgt = j;
            if (j < yng - 1) tgt = yng - 1;   // fold young tail into the yng bin
            if (j > old - 1) tgt = old - 1;   // fold old tail into the old bin
            obs_acc(b * nkeep + tgt - (yng - 1)) += comp_obs_tmp(b * nbins_blk + j);
            hat_acc(b * nkeep + tgt - (yng - 1)) += comp_hat_tmp(b * nbins_blk + j);
          }
        }
        comp_obs_tmp = obs_acc;
        comp_hat_tmp = hat_acc;
        n_comp = nblk * nkeep;
      }
    }

    // Folded proportions (pre-offset, pre-numbers) for the AFSC pseudo-likelihood
    // (comp_ll_type == -1), which reads proportions directly rather than the
    // numbers vector built below. Copied after the fold, so case -1 accumulates
    // the tails as cases 0 and 1 do.
    vector<Type> comp_obs_prop = comp_obs_tmp;
    vector<Type> comp_hat_prop = comp_hat_tmp;

    // Add offset (for some reason can't do above in single line....)
    comp_obs_tmp += comp_prop_offset;
    comp_hat_tmp += comp_prop_offset;

    // Convert observed prop to observed numbers
    comp_obs_tmp *= comp_n(comp_ind, 1);
    // Dirichlet-multinomial concentration. The observed-count total comp_obs_tmp.sum() =
    // N*(1 + offset*nbins) is the effective sample size, matching the N ddirmultinom()
    // reads from obs.sum(), the CAAL DM, and the AFSC linear-DM parameterization.
    vector<Type> alphas = comp_obs_tmp.sum() * comp_hat_tmp * DM_pars_comp(flt); // DM alpha
    vector<Type> unweighted_alphas = comp_obs_tmp.sum() * comp_hat_tmp;          // DM alpha

    // Only use years wanted
    if((yr <= endyr) && (yr > 0) && (flt_type(flt) > 0) && (comp_n(comp_ind, 1) > 0)){

      if(osa_mode == 0){
        // ---- Normal fitting: weighted density read from comp_obs (unchanged) ----
        switch(comp_ll_type(flt)){

        case -1:
          for(ln = 0; ln < n_comp; ln++) {
            // Martin's -- read the folded proportions (comp_*_prop), so an active
            // Comp_accum_young/old actually accumulates the tails here too.
            jnll_comp(JNLL_COMP, flt) -= comp_weights(flt) * Type(comp_n(comp_ind, 1)) * (comp_obs_prop(ln) + comp_prop_offset) * log((comp_hat_prop(ln)+comp_prop_offset) / (comp_obs_prop(ln) + comp_prop_offset)) ;
            unweighted_jnll_comp(JNLL_COMP, flt) -= Type(comp_n(comp_ind, 1)) * (comp_obs_prop(ln) + comp_prop_offset) * log((comp_hat_prop(ln)+comp_prop_offset) / (comp_obs_prop(ln) + comp_prop_offset));
          }
          break;

        case 0: {  // Full multinomial -- via the OSA conditional-binomial decomposition (keep == 1)
          data_indicator<vector<Type>, Type> keep_ones(comp_obs_tmp, true);
          jnll_comp(JNLL_COMP, flt) -= comp_weights(flt) * dmultinom_osa(comp_obs_tmp, comp_hat_tmp, keep_ones, 1, 1);
          unweighted_jnll_comp(JNLL_COMP, flt) -= dmultinom_osa(comp_obs_tmp, comp_hat_tmp, keep_ones, 1, 1);
          break;
        }

        case 1:  // Dirichlet-multinomial
          jnll_comp(JNLL_COMP, flt) -= ddirmultinom(comp_obs_tmp, alphas,  true);
          unweighted_jnll_comp(JNLL_COMP, flt) -= ddirmultinom(comp_obs_tmp, unweighted_alphas,  true);
          break;
        default:
          error("Invalid 'comp_ll_type'");
        }
      } else {
        // ---- OSA build: unweighted, keep-gated conditional density read from
        // obsvec so oneStepPredict() can residualize each bin. The AFSC
        // pseudo-likelihood (-1) has no proper density, so it is residualized
        // under the full multinomial. ----
        int start = comp_obsvec_idx(comp_ind);
        if(start >= 0){
          vector<Type> osa_x = obsvec.segment(start, n_comp);
          if(comp_ll_type(flt) == 1){     // Dirichlet-multinomial (uses fitted DM par)
            // No CDF term: the beta-binomial has no elementary one (comp_osa.hpp).
            jnll_comp(JNLL_COMP, flt) -= ddirmultinom_osa(osa_x, alphas, keep.segment(start, n_comp), 1, 1);
          } else {                        // multinomial (cases 0 and -1)
            jnll_comp(JNLL_COMP, flt) -= dmultinom_osa(osa_x, comp_hat_tmp, keep.segment(start, n_comp), 1, 1, osa_mode == 2);
          }
        }
      }
    }

    // -- Simulate the composition, for sim_mod(simulate = TRUE) --
    // Drawn in RAW bin space, from comp_hat before the transforms above. Tail
    // accumulation folds bins many-to-one before the density, so a draw taken there has
    // no way back; drawing raw and letting the refit fold is exact, both families being
    // closed under merging categories. comp_prop_offset is affine on obs and prediction
    // alike, so drawing from the un-offset prediction keeps the round trip unbiased.
    //
    // Outside the year gate, deliberately: run_mse() reveals the negative-Year rows to
    // the estimation model as new composition data.
    SIMULATE {
      int nbin_raw = ((comp_type == 0) ? nages(sp) : nlengths(sp)) * joint_adjust(comp_ind);
      vector<Type> sim_p = comp_hat.row(comp_ind).segment(0, nbin_raw);
      Type sim_tot = sim_p.sum();
      Type n_nom = comp_n(comp_ind, 1);
      // A row the draw cannot be taken for comes back EMPTY, not carried over and not
      // filled in from the prediction. Two ways to get there:
      //   sim_tot == 0  nothing is predicted -- a fleet with no catch that year, or one
      //                 switched off.
      //   n_nom  == 0   nothing was sampled -- no sample size to draw at.
      // comp_hat rows are normalized to sum to one, so falling back to the prediction
      // would hand back noise-free proportions summing to one: a row that was never
      // sampled, indistinguishable from a real full-weight composition. run_mse() decides
      // "no catch, so no composition" by whether the row sums above zero
      // (R/10-run_mse.R), zeroes Sample_size for empty rows and feeds them back in, so
      // the second case is reachable on the next assessment; sim_mod() warns on an empty
      // row the same way.
      vector<Type> sim_out(nbin_raw); sim_out.setZero();

      if((sim_tot > Type(0)) && (n_nom > Type(0))){
        vector<Type> p = sim_p / sim_tot;
        if(comp_ll_type(flt) == 1){
          // Dirichlet-multinomial: the concentration already carries the weight
          // (DM_pars_comp = exp(comp_weights)), so the draw is at the nominal N.
          vector<Type> sim_alphas = p * n_nom * DM_pars_comp(flt);
          sim_out = rdirmultinom_rce(n_nom, sim_alphas);
        } else {
          // Multinomial (and the AFSC pseudo-likelihood, which has no proper
          // density): comp_weights multiplies the log-likelihood, so it is an
          // effective sample size and the draw must match it, or the data are
          // more informative than the estimator treats them as being.
          sim_out = rmultinom_rce(n_nom * comp_weights(flt), p);
        }
      }
      for(int b = 0; b < nbin_raw; b++){
        comp_sim(comp_ind, b) = sim_out(b);
      }
      // Clear the ragged tail. comp_obs is sized to the widest row, so a row with fewer
      // bins has padding columns belonging to no bin; rearrange_data() normalizes across
      // ALL columns on the next fit, so leaving the padding would dilute every real bin.
      for(int b = nbin_raw; b < comp_sim.cols(); b++){
        comp_sim(comp_ind, b) = Type(0);
      }
    }
  }

  SIMULATE {
    matrix<Type> comp_obs_sim = comp_sim;
    REPORT(comp_obs_sim);
  }


  // Slot 3 -- Conditional age at length (CAAL)
  // Simulated CAAL goes to this copy, not into caal_obs, for the reason comp_sim
  // and diet_sim exist: caal_obs is REPORTed further down.
  matrix<Type> caal_sim = caal_obs;

  for(int caal_ind = 0; caal_ind < caal_hat.rows(); caal_ind++){

    flt = caal_ctl(caal_ind, 0) - 1;            // Fleet index (0-based)
    sp = caal_ctl(caal_ind, 1) - 1;             // Species
    flt_sex = caal_ctl(caal_ind, 2);            // CAAL sex: 0 = combined, 1 = female, 2 = male, 3 = joint
    yr = caal_ctl(caal_ind, 3);                 // Data year
    ln = caal_ctl(caal_ind, 4) - 1;             // Length bin

    // Number of ages
    int n_caal = nages(sp); // * joint_adjust(caal_ind); // For sex = 3 (joint sex caal data)

    // Select sections
    vector<Type> caal_obs_tmp = caal_obs.row(caal_ind).segment(0, n_caal); // Observed proportion
    vector<Type> caal_hat_tmp = caal_hat.row(caal_ind).segment(0, n_caal); // Expected proportion

    // Proportions kept pre-offset and pre-numbers for the AFSC pseudo-likelihood
    // (caal_ll_type == -1), which reads proportions directly. Mirrors the
    // comp_obs_prop / comp_hat_prop copies in the age-composition slot above.
    vector<Type> caal_obs_prop = caal_obs_tmp;
    vector<Type> caal_hat_prop = caal_hat_tmp;

    // Add offset (for some reason can't do above in single line....)
    caal_obs_tmp += comp_prop_offset;
    caal_hat_tmp += comp_prop_offset;

    // Convert observed prop to observed numbers
    caal_obs_tmp *= caal_n(caal_ind, 0);
    vector<Type> alphas = sum(caal_obs_tmp) * caal_hat_tmp * DM_pars_caal(flt); // DM alpha
    vector<Type> unweighted_alphas = sum(caal_obs_tmp) * caal_hat_tmp;          // DM alpha

    // Only use years wanted
    if((yr <= endyr) && (yr > 0) && (flt_type(flt) > 0) && (caal_n(caal_ind, 0) > 0)){

      if(osa_mode == 0){
        // ---- Normal fitting ----
        // Two parallel accumulators per CAAL fleet: jnll_comp, weighted by
        // caal_weights(flt) and entering the objective, and unweighted_jnll_comp,
        // reported so Francis and McAllister-Ianelli reweighting can read the raw
        // likelihood. Both address JNLL_CAAL, not JNLL_COMP: age composition and
        // conditional age-at-length are separate likelihood rows (see the JnllRow enum).
        switch(caal_ll_type(flt)){

        case -1:  // AFSC/AMAK multinomial pseudo-likelihood, as for age comps
          // Same form as comp_ll_type == -1: sum_a N * (obs_a + off) *
          // log((hat_a + off) / (obs_a + off)), read off the proportions rather
          // than the numbers vector. It drops the multinomial normalizing
          // constant, so it is a pseudo-likelihood and not a density -- hence the
          // OSA branch residualizes it under the full multinomial, as comps do.
          for(int age_ind = 0; age_ind < n_caal; age_ind++){
            jnll_comp(JNLL_CAAL, flt) -= caal_weights(flt) * Type(caal_n(caal_ind, 0)) * (caal_obs_prop(age_ind) + comp_prop_offset) * log((caal_hat_prop(age_ind) + comp_prop_offset) / (caal_obs_prop(age_ind) + comp_prop_offset));
            unweighted_jnll_comp(JNLL_CAAL, flt) -= Type(caal_n(caal_ind, 0)) * (caal_obs_prop(age_ind) + comp_prop_offset) * log((caal_hat_prop(age_ind) + comp_prop_offset) / (caal_obs_prop(age_ind) + comp_prop_offset));
          }
          break;

        case 0: {  // Full multinomial -- via the OSA conditional-binomial decomposition (keep == 1)
          data_indicator<vector<Type>, Type> keep_ones(caal_obs_tmp, true);
          jnll_comp(JNLL_CAAL, flt) -= caal_weights(flt) * dmultinom_osa(caal_obs_tmp, caal_hat_tmp, keep_ones, 1, 1);
          unweighted_jnll_comp(JNLL_CAAL, flt) -= dmultinom_osa(caal_obs_tmp, caal_hat_tmp, keep_ones, 1, 1);
          break;
        }

        case 1:  // Dirichlet-multinomial
          jnll_comp(JNLL_CAAL, flt) -= ddirmultinom(caal_obs_tmp, alphas,  true);
          unweighted_jnll_comp(JNLL_CAAL, flt) -= ddirmultinom(caal_obs_tmp, unweighted_alphas,  true);
          break;
        default:
          error("Invalid 'caal_ll_type'");
        }
      } else {
        // ---- OSA build: unweighted, keep-gated conditional density read from
        // obsvec (alpha total uses the fixed caal_obs counts, not obsvec). ----
        int start = caal_obsvec_idx(caal_ind);
        if(start >= 0){
          vector<Type> osa_x = obsvec.segment(start, n_caal);
          if(caal_ll_type(flt) == 1){     // Dirichlet-multinomial
            // No CDF term: the beta-binomial has no elementary one (comp_osa.hpp).
            jnll_comp(JNLL_CAAL, flt) -= ddirmultinom_osa(osa_x, alphas, keep.segment(start, n_caal), 1, 1);
          } else {                        // multinomial
            jnll_comp(JNLL_CAAL, flt) -= dmultinom_osa(osa_x, caal_hat_tmp, keep.segment(start, n_caal), 1, 1, osa_mode == 2);
          }
        }
      }
    }

    // -- Simulate the conditional age-at-length, for sim_mod(simulate = TRUE) --
    // Same construction as the marginal composition above and for the same
    // reasons: drawn from caal_hat before comp_prop_offset is added, and outside
    // the year gate. CAAL has no tail accumulation, so no fold to undo.
    SIMULATE {
      vector<Type> sim_p = caal_hat.row(caal_ind).segment(0, n_caal);
      Type sim_tot = sim_p.sum();
      Type n_nom = caal_n(caal_ind, 0);
      vector<Type> sim_out(n_caal); sim_out.setZero();

      if((sim_tot > Type(0)) && (n_nom > Type(0))){
        vector<Type> p = sim_p / sim_tot;
        if(caal_ll_type(flt) == 1){
          vector<Type> sim_alphas = p * n_nom * DM_pars_caal(flt);
          sim_out = rdirmultinom_rce(n_nom, sim_alphas);
        } else {
          sim_out = rmultinom_rce(n_nom * caal_weights(flt), p);
        }
      }
      // As for the marginal composition: a row with nothing predicted or nothing sampled
      // comes back empty rather than filled in from the prediction, and the ragged tail
      // is cleared so a later normalization cannot dilute the real bins.
      for(int b = 0; b < n_caal; b++){
        caal_sim(caal_ind, b) = sim_out(b);
      }
      for(int b = n_caal; b < caal_sim.cols(); b++){
        caal_sim(caal_ind, b) = Type(0);
      }
    }
  }

  SIMULATE {
    matrix<Type> caal_obs_sim = caal_sim;
    REPORT(caal_obs_sim);
  }


  // Slot 4-5 -- Selectivity
  // FIXME: penalize every selectivity deviation rather than a sub-range. The
  // range controls (Sel_pen_first_bin, Sel_pen_last_bin, Sel_cap_bin,
  // Sel_start_year and the matching build_map masking) exist only because these
  // penalties cover a sub-range; penalizing all deviations would pin the
  // unidentified directions and remove the need to index by bin and year.
  for(flt = 0; flt < n_flt; flt++){ // Every fleet; the gate below picks the lead ones
    jnll_comp(JNLL_SEL_NONPARAM, flt) = 0;
    jnll_comp(JNLL_SEL_DEV, flt) = 0;
    sp = flt_spp(flt);

    // Non-parametric penalties act over the fleet's selectivity dimension:
    // nbins = nages for age-based, nlengths for length-based selectivity.
    bool sel_is_length = (flt_sel_dim(flt) == 1);
    // sel_at_length is indexed on the POPULATION bins, so the penalty's bin
    // range is counted there. Identical to the data bins with one grid, and
    // data_check() refuses a time-varying length-based selectivity on a coarser
    // data grid, where Sel_pen_first_bin / Sel_pen_last_bin would be ambiguous.
    int  nbins = sel_is_length ? nlengths_pop(sp) : nages(sp);

    // If estimating survey or fishery (and not a selectivity mirror of an earlier
    // fleet - the shared penalty is accumulated once, on the lead fleet).
    if(flt_type(flt) > 0 && flt_sel_lead(flt) == 1){

      // 1a) Ianelli/AMAK non-parametic selectivity penalties
      // - using non-normalized selectivities following the arrowtooth ADMB model
      // - differentiable via abs: only penalize sel_ratio_tmp > 0 (decreasing sel_at_age)
      // - Added time-varying component following Atka mackerel
      if(flt_sel_type(flt) == 2) {

        // Time-invariant selectivity has one curve, so the shape penalties are charged
        // once; both time-varying modes give each year its own curve, charged each year.
        //
        // PROPOSED, to make IID integrable under random_sel = TRUE: charge terms 1, 2 and
        // 4 once on the base curve sel_coff instead of per year on the realized curve,
        // term 1 over bins below flt_n_sel_bins - 1 (past which sel_coff holds 0 rather
        // than the last estimated bin the realized curve repeats) and term 4 over
        // [bin_first_selected, flt_n_sel_bins). Term 3 would then be the only term on the
        // deviates and WOULD NEED its own for(y < nyrs_hind) loop here, since nyrs_tmp
        // becomes 1 for IID. The reason to move them: on the deviates, term 1 steps the
        // Laplace log-determinant (max(d, 0)^2 has a step second derivative), and terms 2
        // and 4 do not scale with sel_dev_sd, so they integrate an unnormalized prior and
        // bias it. Type 13 below is the same idea as its own form, but extends the base
        // by repeating the last estimated bin and takes its level over [0, n_sel_bins).
        int nyrs_tmp = 1;
        if((flt_varying_sel(flt) == 1) || (flt_varying_sel(flt) == 4)){
          nyrs_tmp = nyrs_hind;
        }

        for(yr = 0; yr < nyrs_tmp; yr++){

          // 1. Decreasing selectivity penalty (over the fleet's own bins)
          // FIXME: AMAK starts at nbins/2
          for(sex = 0; sex < nsex(sp); sex++){
            for(age = 0; age < (nbins - 1); age++) {
              Type sel_ratio_tmp = log_non_par_sel(flt, sex, age, yr) - log_non_par_sel(flt, sex, age + 1, yr); // Positive if decreasing
              jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0) * square( (CppAD::abs(sel_ratio_tmp) + sel_ratio_tmp)/2.0);
            }
          }

          // 2. Curvature penalty
          for(sex = 0; sex < nsex(sp); sex++){
            // Extract only the selectivities we want
            vector<Type> sel_tmp(nbins); sel_tmp.setZero();

            for(age = 0; age < nbins; age++) {
              sel_tmp(age) = log_non_par_sel(flt, sex, age, yr);
            }

            // Second difference computed once (matches the type-9 branch).
            vector<Type> sel_d2 = first_difference( first_difference( sel_tmp ) );
            for(int a2 = 0; a2 < sel_d2.size(); a2++) {
              jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 1) * sel_d2(a2) * sel_d2(a2);
            }
          }

          // 3. Time-varying penalty
          //
          // RandomWalk (4) scores the year-to-year change in the realized
          // log-selectivity (Atka mackerel ADMB model), which is renormalized to mean 1
          // each year; the level is held by the avgsel penalty at (4). IID (1) scores
          // each deviation itself, over exactly the coefficient bins build_map()
          // estimates. Bins outside that range are fixed at 0 and left unscored: each
          // would add log(sel_dev_sd) + 0.5*log(2*pi), constant in the deviates but
          // rising with the sd, so scoring them would drag the estimate toward zero.
          if(flt_varying_sel(flt) == 4){
            if(yr > 0){
              for(sex = 0; sex < nsex(sp); sex++){
                for(age = 0; age < (nbins - 1); age++) {
                  jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_non_par_sel(flt, sex, age, yr), log_non_par_sel(flt, sex, age, yr - 1), sel_dev_sd(flt), true);
                }
              }
            }
          } else if(flt_varying_sel(flt) == 1){
            for(sex = 0; sex < nsex(sp); sex++){
              for(int bin = bin_first_selected(flt); bin < flt_n_sel_bins(flt); bin++) {
                jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_coff_dev(flt, sex, bin, yr), Type(0.0), sel_dev_sd(flt), true);
              }
            }
          }

          // 4. Survey selectivity normalization (non-parametric). Reads the
          //    per-year avg_sel, which includes the deviates. Type 13 takes the
          //    level from sel_coff with log_mean_exp(); type 9 does too, but only
          //    when Sel_avgsel_pen is set (0 = off by default).
          for(sex = 0; sex < nsex(sp); sex++){
            jnll_comp(JNLL_SEL_NONPARAM, flt) += 2.0 * square(avg_sel(flt, sex, yr));
          }
        }
      }

      // 1b) Non-parametric with the ADMB ("pm"/AMAK) selectivity penalty
      //     (NonParametricPM, type 9). Construction is identical to type 2; only
      //     the penalty differs, to match ADMB's Selectivity_Likelihood:
      //       jnll_comp(4) = ADMB sel_like(1): decreasing-only penalty, weight
      //         sel_curve_pen(flt,0) (= ctrl_flag(13)).
      //       jnll_comp(5) = ADMB sel_like_dev(1): curvature, weight
      //         sel_curve_pen(flt,1) (= ctrl_flag(11)/nch) on the 2nd difference
      //         of log-selectivity every year; a random walk, bare Gaussian SSQ
      //         of the year-to-year change in log-selectivity / (2*sd^2) over all
      //         ages, with NO dnorm normalizing constant; and dev magnitude,
      //         weight sel_curve_pen(flt,2) (= ctrl_flag(10)/group_num) on norm2
      //         of the year-to-year coefficient increments.
      if(flt_sel_type(flt) == 9) {

        int nyrs_tmp = 1;
        if(flt_varying_sel(flt) == 4){
          nyrs_tmp = nyrs_hind;
        }

        // Per-fleet selectivity start year (0-based). Penalties over realized
        // selectivity skip pre-survey years, and the survey base-year curvature
        // with them, matching ADMB's styr-anchored base term, which is 0 for a
        // survey starting after styr. A fleet starting at styr has start_yr = 0.
        int start_yr = flt_sel_start_yr(flt);
        int shape_a0 = flt_sel_pen_first_bin(flt);              // first (left) bin of the shape-penalty pairs
        if(shape_a0 < 0) shape_a0 = bin_first_selected(flt);    // default: first selected bin
        int shape_a1 = flt_sel_pen_last_bin(flt);               // last (left) bin of the shape-penalty pairs
        if(shape_a1 < 0) shape_a1 = nbins - 2;                  // default: whole range (pairs up to (nbins-2, nbins-1))
        // The pair is (age, age + 1), so the left bin stops one short of the
        // last. data_check() bounds Sel_pen_last_bin by the number of bins, so
        // without this a final-bin setting reads bin `nbins` -- allocated to the
        // widest species, never written, scoring as a phantom selected bin.
        if(shape_a1 > nbins - 2) shape_a1 = nbins - 2;
        int shape_mode = flt_sel_shape_mode(flt);               // 0 = directional (sign of pen), 1 = smooth (two-sided d^2)
        for(sex = 0; sex < nsex(sp); sex++){
          for(yr = start_yr; yr < nyrs_tmp; yr++){

            // (1) Shape penalty over adjacent ages [shape_a0 .. shape_a1], all active years.
            //     mode 0 (directional, ADMB/AMAK): SIGN of sel_curve_pen(flt,0) sets
            //       direction (>=0 penalize DECREASING, <0 penalize INCREASING), one-sided,
            //       differentiable via (|d|+d)/2 = max(d,0) and (|d|-d)/2 = max(-d,0).
            //     mode 1 (smooth, RTMB "rpm"): two-sided weight * d^2 smoothness.
            for(age = shape_a0; age <= shape_a1; age++) {
              Type d = log_non_par_sel(flt, sex, age, yr) - log_non_par_sel(flt, sex, age + 1, yr); // > 0 if decreasing
              if(shape_mode == 1)
                jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0) * d * d;                                 // two-sided smoothness
              else if(sel_curve_pen(flt, 0) >= 0)
                jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0)  * square( (CppAD::abs(d) + d)/2.0 );    // penalize decreasing
              else
                jnll_comp(JNLL_SEL_NONPARAM, flt) += -sel_curve_pen(flt, 0) * square( (CppAD::abs(d) - d)/2.0 );    // penalize increasing
            }

            // (2) Curvature (2nd-difference) penalty  [ADMB term 2/3]. Includes the
            //     base year only when the fleet starts at styr (start_yr == 0);
            //     otherwise the survey base-year curvature is excluded (ADMB anchors
            //     the base curvature term at styr, which is 0 pre-survey).
            if((yr > start_yr) || (start_yr == 0)){
              vector<Type> ls(nbins); ls.setZero();
              for(age = 0; age < nbins; age++) ls(age) = log_non_par_sel(flt, sex, age, yr);
              vector<Type> d2 = first_difference( first_difference( ls ) );
              for(int a2 = 0; a2 < d2.size(); a2++) jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 1) * d2(a2) * d2(a2);
            }

            // (3) Random-walk penalty (bare SSQ, no normalizing constant)  [ADMB term 4]
            if(yr > start_yr){
              for(age = 0; age < nbins; age++) {
                Type dd = log_non_par_sel(flt, sex, age, yr) - log_non_par_sel(flt, sex, age, yr - 1);
                jnll_comp(JNLL_SEL_DEV, flt) += dd * dd / (2.0 * sel_dev_sd(flt) * sel_dev_sd(flt));
              }
            }
          }

          // (4) Dev-magnitude penalty: norm2 of the RAW per-year increments
          //     (sel_coff_dev IS the random-walk increment for NonParametricPM;
          //     = RTMB norm2(sel_devs)). Increments are 0 at non-change years.
          for(int bin = 0; bin < flt_n_sel_bins(flt); bin++){
            for(yr = start_yr; yr < nyrs_tmp; yr++){
              jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 2) * sel_coff_dev(flt, sex, bin, yr) * sel_coff_dev(flt, sex, bin, yr);
            }
          }

          // (5) AMAK "avgsel" base-level penalty: weight * (log(mean(exp(base
          //     coffs))))^2 over the estimated coefficient bins
          //     [bin_first_selected, n_sel_bins-1]. A mild regulariser on the overall
          //     level of the base coefficients, which the per-year mean-centering leaves
          //     unconstrained; equivalent to AMAK's 10*square(avgsel_*).
          //     flt_sel_avgsel_pen defaults to 0 (off); opt in via Sel_avgsel_pen (e.g.
          //     10 to match AMAK). The range holds at least one bin, data_check()
          //     requiring Bin_first_selected <= N_sel_bins.
          if(flt_sel_avgsel_pen(flt) > 0){
            int n_base = flt_n_sel_bins(flt) - bin_first_selected(flt);
            vector<Type> base(n_base);
            for(int bin = 0; bin < n_base; bin++){
              base(bin) = sel_coff(flt, sex, bin + bin_first_selected(flt));
            }
            Type avgsel = log_mean_exp(base);
            jnll_comp(JNLL_SEL_NONPARAM, flt) += flt_sel_avgsel_pen(flt) * avgsel * avgsel;
          }
        }
      }


      // 1c) NonParametricIntegrable (13): the Ianelli shape priors (decreasing,
      //     curvature, average selectivity) are charged ONCE on the base curve
      //     sel_coff, and the deviates carry a proper Gaussian density with
      //     sel_dev_sd, so under random_sel = TRUE the Laplace approximation
      //     integrates a density whose normalizing constant is complete. With no
      //     deviates the objective equals NonParametric's.
      if(flt_sel_type(flt) == 13) {
        int n_sel_bins = flt_n_sel_bins(flt);
        int start_yr   = flt_sel_start_yr(flt);
        for(sex = 0; sex < nsex(sp); sex++){

          // The base curve as selectivity.hpp builds a year without deviates:
          // the estimated bins, the tail held at the last estimated bin, the
          // whole vector centred to mean exp() of 1. Bins below the first
          // selected bin read the 0 that build_map() leaves in sel_coff there,
          // exactly as the realized curve does before the normalizer zeroes them.
          vector<Type> base(nbins);
          for(int bin = 0; bin < nbins; bin++){
            base(bin) = sel_coff(flt, sex, (bin < n_sel_bins ? bin : n_sel_bins - 1));
          }
          vector<Type> est(n_sel_bins);
          for(int bin = 0; bin < n_sel_bins; bin++) est(bin) = base(bin);
          Type avg_base = log_mean_exp(est);
          Type ctr = log_mean_exp(base);
          for(int bin = 0; bin < nbins; bin++) base(bin) -= ctr;

          // 1. Decreasing-selectivity penalty on the base.
          for(int bin = 0; bin < (nbins - 1); bin++) {
            Type d = base(bin) - base(bin + 1);
            jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 0) * square( (CppAD::abs(d) + d)/2.0 );
          }
          // 2. Curvature penalty on the base.
          vector<Type> d2 = first_difference( first_difference( base ) );
          for(int a2 = 0; a2 < d2.size(); a2++) {
            jnll_comp(JNLL_SEL_NONPARAM, flt) += sel_curve_pen(flt, 1) * d2(a2) * d2(a2);
          }
          // 4. Average-selectivity level of the base coefficients.
          jnll_comp(JNLL_SEL_NONPARAM, flt) += 2.0 * square(avg_base);

          // 3. The deviates: iid about the base over every hindcast year, or random-walk
          //    increments from the year after the fleet's start year (the start-year
          //    increment is fixed at 0). Only estimated coefficient bins are scored; a
          //    bin held at 0 would add a constant rising with the sd and pull it toward
          //    zero. With Time_varying_sel = "Off" there are no deviates and no density.
          bool scored = (flt_varying_sel(flt) == 1) || (flt_varying_sel(flt) == 4);
          int yr_lo = (flt_varying_sel(flt) == 4) ? start_yr + 1 : 0;
          if(scored){
            for(yr = yr_lo; yr < nyrs_hind; yr++){
              for(int bin = bin_first_selected(flt); bin < n_sel_bins; bin++) {
                jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_coff_dev(flt, sex, bin, yr), Type(0.0), sel_dev_sd(flt), true);
              }
            }
          }
        }
      }


      // 2) Logistic selectivity penalties
      // Penalized/random-effect density: time-varying logistic/double-logistic deviates
      if(((flt_varying_sel(flt) == 1)||(flt_varying_sel(flt) == 2)) && (flt_sel_type(flt) != 2) && (flt_sel_type(flt) != 5) && (flt_sel_type(flt) != 11) && (flt_sel_type(flt) != 13)){
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 0; yr < nyrs_hind; yr++){

            // Logistic / ascending-limb deviates (types 1, 3, 8)
            // For DoubleNormal (8): sel_inf_dev(0) = peak deviate, log_sel_slp_dev(0) = ascending-SD deviate
            if((flt_sel_type(flt) == 1) || (flt_sel_type(flt) == 3) || (flt_sel_type(flt) == 8)){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(0, flt, sex, yr), Type(0.0), 4 * sel_dev_sd(flt), true);
            }

            // Double logistic / descending-limb deviates (types 3, 4, 8)
            // For DoubleNormal (8): sel_inf_dev(1) = right-floor logit deviate; log_sel_slp_dev(1) = descending-SD deviate
            if((flt_sel_type(flt) == 3) || (flt_sel_type(flt) == 4) || (flt_sel_type(flt) == 8)){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(1, flt, sex, yr), Type(0.0), 4 * sel_dev_sd(flt), true);
            }
          }
        }
      }

      // Random walk:
      // - Type 4 = random walk on ascending and descending for double logistic
      // - Type 5 = ascending only for double logistics
      if(((flt_varying_sel(flt) == 4)||(flt_varying_sel(flt) == 5)) && (flt_sel_type(flt) != 2) && (flt_sel_type(flt) != 5) && (flt_sel_type(flt) != 11) && (flt_sel_type(flt) != 13)){
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = 1; yr < nyrs_hind; yr++){ // Start at second year

            // Logistic / ascending-limb random walk (types 1, 3, 8)
            if((flt_sel_type(flt) == 1) || (flt_sel_type(flt) == 3) || (flt_sel_type(flt) == 8)){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(0, flt, sex, yr) - log_sel_slp_dev(0, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(0, flt, sex, yr) - sel_inf_dev(0, flt, sex, yr-1), Type(0.0), 4 * sel_dev_sd(flt), true);
            }

            // Double logistic / descending-limb random walk (types 3, 4, 8).
            // Gated on mode 4 as well as the type: mode 5 (RandomWalkAscending)
            // varies the ascending limb ONLY, and build_map() estimates no
            // descending deviate under it for any type -- type 3 restricts to
            // j = 1, types 4 and 8 exclude mode 5 outright -- so those deviates
            // sit at their init of 0. Under random_sel = FALSE the penalty would
            // then be a pure constant, shifting the reported objective and no
            // gradient; under random_sel = TRUE sel_dev_log_sd is estimated, so
            // the term is 2n log(sigma) + const and biases that SD downward.
            if((flt_varying_sel(flt) == 4) &&
               ((flt_sel_type(flt) == 3) || (flt_sel_type(flt) == 4) || (flt_sel_type(flt) == 8))){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_inf_dev(1, flt, sex, yr) - sel_inf_dev(1, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt), true);
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(log_sel_slp_dev(1, flt, sex, yr) - log_sel_slp_dev(1, flt, sex, yr-1), Type(0.0), sel_dev_sd(flt) * 4, true);
            }
          }
        }
      }

      // 2b) LogisticPM (type 11): ADMB AMAK ("pm") bottom-trawl-survey penalties,
      //     ctrl_flag(19) > 0 branch. Two random-walk (norm2 of first-difference)
      //     terms, BOTH starting the year AFTER flt_sel_start_yr so pre-survey
      //     years and the start-year boundary jump are excluded (ADMB's
      //     dev_vector / flat pre-survey selectivity contribute ~0 there):
      //       (1) RW on the REALIZED log-selectivity-at-age over the penalty
      //           age-range [sel_norm_bin1, sel_norm_bin2] (= ADMB ctrl_flag(26) *
      //           sum_{q_amin..q_amax-1} norm2(first_difference(log_sel_bts))),
      //           weight sel_curve_pen(flt,0);
      //       (2) RW on the free age-1 parameter deviates (= ADMB
      //           8 * norm2(first_difference(sel_age_one_bts_dev))),
      //           weight sel_curve_pen(flt,2).
      //     sel_curve_pen(flt,1) is unused in this branch.
      if(flt_sel_type(flt) == 11){
        int start_yr = flt_sel_start_yr(flt);                 // 0-based; first first-difference is start_yr+1
        int alo = sel_norm_bin1(flt); int ahi = sel_norm_bin2(flt);
        if(alo < 0) alo = bin_first_selected(flt);            // default: whole selected range
        if(ahi < 0) ahi = nbins - 1;
        for(sex = 0; sex < nsex(sp); sex ++){
          for(yr = start_yr + 1; yr < nyrs_hind; yr++){
            // (1) realized log-selectivity random walk over the penalty bin-range
            for(age = alo; age <= ahi; age++){
              Type s_now  = sel_is_length ? sel_at_length(flt, sex, age, yr)     : sel_at_age(flt, sex, age, yr);
              Type s_prev = sel_is_length ? sel_at_length(flt, sex, age, yr - 1) : sel_at_age(flt, sex, age, yr - 1);
              Type d = log(s_now) - log(s_prev);
              jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 0) * d * d;
            }
            // (2) free age-1 parameter-deviate random walk
            Type da1 = sel_inf_dev(1, flt, sex, yr) - sel_inf_dev(1, flt, sex, yr - 1);
            jnll_comp(JNLL_SEL_DEV, flt) += sel_curve_pen(flt, 2) * da1 * da1;
          }
        }
      }

      // 3) Hake style non-parametric
      // Penalized/random-effect time-varying sel deviates (Taylor et al 2014)
      if(((flt_varying_sel(flt) == 1) || (flt_varying_sel(flt) == 2)) && (flt_sel_type(flt) == 5)){
        for(int bin = 0; bin < flt_n_sel_bins(flt); bin++){ //NOTE: extends beyond selectivity age range, but should be mapped to 0 in map function
          for(sex = 0; sex < nsex(sp); sex++){
            for(yr = 0; yr < nyrs_hind; yr++){
              jnll_comp(JNLL_SEL_DEV, flt) -= dnorm(sel_coff_dev(flt, sex, bin, yr), Type(0.0), sel_dev_sd(flt), true);
            }
          }
        }
      }


      // 4) 2D AR1
      if(flt_sel_type(flt) == 6){
        int n_sel_bins = flt_n_sel_bins(flt);
        array<Type> tmp_AR2(n_sel_bins, nyrs_hind);


        for(sex = 0; sex < nsex(sp); sex ++){
          for (int bin = 0; bin < n_sel_bins; bin++) {
            for (int yr = 0; yr < nyrs_hind; yr++) {
              tmp_AR2(bin, yr) =  sel_coff_dev(flt, sex, bin, yr);
            }
          }

          Type rho_a = rho_trans(sel_curve_pen(flt, 0)); // Scale from -1 to 1
          Type rho_y = rho_trans(sel_curve_pen(flt, 1));

          Type Sigma_sig_sel = pow(pow(sel_dev_sd(flt),2) / ((1-pow(rho_y,2))*(1-pow(rho_a,2))),0.5);
          // As in the M 2D-AR1 above: SEPARABLE(f, g) puts f on the outermost
          // dimension. tmp_AR2 is (bin, year), so year leads and bin follows.
          jnll_comp(JNLL_SEL_DEV, flt) += SCALE(SEPARABLE(AR1(rho_y),AR1(rho_a)), Sigma_sig_sel)(tmp_AR2);
        } // end sex loop
      }

      // 5) 3D AR1
      if(flt_sel_type(flt) == 7){
        // Build bin-year index
        int n_sel_bins = flt_n_sel_bins(flt);
        int num_rows = n_sel_bins * nyrs_hind;
        matrix<Type> ay_index(num_rows, 2);

        for(sex = 0; sex < nsex(sp); sex ++){

          int row = 0;
          for (int year = 1; year <= nyrs_hind; ++year) {
            for (int bin = 1; bin <= n_sel_bins; ++bin) {
              ay_index(row, 0) = bin;
              ay_index(row, 1) = year;
              row++;
            }
          }

          // -- Define precision matrix for GMRF
          Eigen::SparseMatrix<Type> Q_sparse(num_rows, num_rows); // Precision matrix

          // -- Construct precision matrix
          // Argument order is bin, year, cohort. construct_Q()'s own parameter
          // names read the other way round, but its `age > 1` branch (previous
          // bin, same year) multiplies the FIRST argument, and this ay_index is
          // filled (bin, year) -- the mirror of WHAM's -- so the two inversions
          // cancel. Do not "correct" this call to match those names.
          Q_sparse = construct_Q(nyrs_hind, n_sel_bins, ay_index,
                                 sel_curve_pen(flt, 0), sel_curve_pen(flt, 1), sel_curve_pen(flt, 2), // natural scale
                                 log(square(sel_dev_sd(flt))), 0); // Conditional variance


          // --- Derive likelihood
          array<Type> tmp_AR2(n_sel_bins, nyrs_hind);
          for (int bin = 0; bin < n_sel_bins; bin++) {
            for (int yr = 0; yr < nyrs_hind; yr++) {
              tmp_AR2(bin, yr) =  sel_coff_dev(flt, sex, bin, yr);
            }
          }
          jnll_comp(JNLL_SEL_DEV, flt) += GMRF(Q_sparse)(tmp_AR2);
        }
      }
    }
  } // End selectivity loop


  // Slot 6-7 -- Catchability
  // Fleets sharing a Q_index estimate one catchability and one deviate vector,
  // so only the lead fleet accumulates the prior and deviate penalties.
  for(flt = 0; flt < n_flt; flt++){
   if(flt_q_lead(flt) == 1){

    // Lognormal prior on catchability: Catchability_init is the mean of q when
    // bias_adjust_proc = 1, the median when 0.
    if( est_index_q(flt) == 2){
      jnll_comp(JNLL_Q_PRIOR, flt) -= dnorm(index_log_q(flt), index_log_q_prior(flt) - bias_adjust_proc*square(index_q_sd(flt))/2.0, index_q_sd(flt), true);
    }

    // Penalized/random deviate likelihood
    if(((index_varying_q(flt) == 1) || (index_varying_q(flt) == 2))  // - Time_varying_q = 1 (IID) or 2 (AR1, retired -- scored as IID)
         && (flt_type(flt) > 0) &&                                    // - If survey or fishery CPUE
           ((est_index_q(flt) == 1) || (est_index_q(flt) == 2))){        // - Catchability = 1 (Estimated) or 2 (Estimated-with-prior)
      for(yr = 0; yr < nyrs_hind; yr++){
        jnll_comp(JNLL_Q_DEV, flt) -= dnorm(index_q_dev(flt, yr), Type(0.0), index_q_dev_sd(flt), true );
      }
    }

    // Random walk
    if((index_varying_q(flt) == 4) &&                          // - Time_varying_q = 4 (RandomWalk)
       (flt_type(flt) > 0) &&                                  // - If survey or fishery CPUE
       ((est_index_q(flt) == 1) || (est_index_q(flt) == 2)))   // - Catchability = 1 (Estimated) or 2 (Estimated-with-prior)
    {
      for(yr = 1; yr < nyrs_hind; yr++){
        jnll_comp(JNLL_Q_DEV, flt) -= dnorm(index_q_dev(flt, yr) - index_q_dev(flt, yr-1), Type(0.0), index_q_dev_sd(flt), true );
      }
    }
   } // End q lead gate
  } // End q loop


  // Slots 8-11 -- Recruitment
  for(sp = 0; sp < nspp; sp++) {
    penalty = 0.0;
    // Input numbers-at-age (estDynamics > 0) leave R a placeholder, so no SRR
    // prior, Bmsy penalty or curve penalty.
    int srr_terms_on = (estDynamics(sp) == 0);

    // Slot 8 -- stock-recruit prior for Beverton
    // -- Lognormal on steepness. srr_prior is its mean when bias_adjust_proc = 1
    //    (centred at -sigma^2/2, as rec_dev), its median when 0.
    if(srr_terms_on && (srr_est_mode == 2) & ((srr_pred_fun == 2) | (srr_pred_fun == 3))){
      jnll_comp(JNLL_SRR_PRIOR, sp) -= dnorm(log(steepness(sp, 0)), log(srr_prior(sp)) - bias_adjust_proc*square(srr_prior_sd(sp))/2.0, srr_prior_sd(sp), true);
    }

    // -- Beta
    if(srr_terms_on && (srr_est_mode == 3) & ((srr_pred_fun == 2) | (srr_pred_fun == 3))){
      // Convert mean and SD to beta params
      Type beta_alpha = ((1 - srr_prior(sp))/ square(srr_prior_sd(sp)) - 1/srr_prior(sp)) * square(srr_prior(sp));
      Type beta_beta = beta_alpha * (1/srr_prior(sp) - 1);
      jnll_comp(JNLL_SRR_PRIOR, sp) -= dbeta(steepness(sp, 0), beta_alpha, beta_beta, true);
    }

    // Slot 8 -- stock-recruit prior for Ricker, lognormal on alpha. srr_prior is alpha's
    //    mean when bias_adjust_proc = 1, its median when 0.
    if(srr_terms_on && (srr_est_mode == 2) & ((srr_pred_fun == 4) | (srr_pred_fun == 5))){
      jnll_comp(JNLL_SRR_PRIOR, sp) -= dnorm(rec_pars(sp, 1), log(srr_prior(sp)) - bias_adjust_proc*square(srr_prior_sd(sp))/2.0, srr_prior_sd(sp), true);
    }

    // Slot 8 -- penalty for Bmsy > Bmsy_lim for Ricker
    if(srr_terms_on && (!isNA(Bmsy_lim(sp))) && ((srr_pred_fun == 4) || (srr_pred_fun == 5))){ // Using pred_fun in case ianelli method is used
      Type bmsy = 1.0/exp(rec_pars(sp, 2));
      bmsy =  posfun(Bmsy_lim(sp)/Type(1000000.0) - bmsy, Type(0.001), penalty);
      jnll_comp(JNLL_SRR_PRIOR, sp) += 100 * penalty;
    }


    // Slot 9 -- init_dev -- Initial abundance-at-age
    // Lognormal bias correction: dev ~ N(-sigma^2/2, sigma) so E[N_init] = deterministic
    // equilibrium. initMode 5 (OffsetEquilibrium) fixes init_dev at 0 (off), like the
    // equilibrium modes, so it carries no init_dev penalty.
    if(initMode > 1 && initMode != 5){
      for(age = 1; age < nages(sp); age++) {
        jnll_comp(JNLL_INIT_DEV, sp) -= dnorm( init_dev(sp, age - 1), -bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp), true);
      }
    }

    // Slot 10 -- Tau -- Annual recruitment deviation
    // Lognormal bias correction: dev ~ N(-sigma^2/2, sigma) so E[R] = R0 (mean-unbiased).
    // Skipped where the SEM drives recruitment: the GMRF in section 5.5b is
    // already the density for these deviations, and applying both would count
    // them twice. Gated on dsem_drives_rec rather than dsem_on, so a SEM that
    // drives only some other process leaves this density in place -- without
    // that, attaching such a SEM would delete the recruitment process while
    // the fit still converged and still reported an ABC.
    if(dsem_drives_rec == 0){
      for(yr = 0; yr < nyrs_hind; yr++) {
        jnll_comp(JNLL_REC_DEV, sp) -= dnorm( rec_dev(sp, yr),  -bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp), true);    // Recruitment deviation using random effects.
      }
    } else if(rec_sd_use_prior(sp) == 1){
      // Optional lognormal prior on the estimated recruitment SD, centred on the
      // assessment value. Regularizes R_sd away from the 1/sigma^2 collapse that
      // occurs when the covariates over-explain the deviations.
      jnll_comp(JNLL_DSEM, sp) -= dnorm( log(R_sd(sp)), log(rec_sd_prior(sp)),
                                         rec_sd_prior_sd(sp), true );
    }

    // Slot 11 -- Additional penalty for SRR curve (sensu AMAK/Ianelli). R_hat is the mean of R
    //    when bias_adjust_proc = 1 (Ianelli's pm.tpl, Dorn 2002), the median (AMAK) when 0.
    if(srr_terms_on && (srr_fun == 0) & (srr_pred_fun  > 0)){
      for(yr = srr_hat_styr; yr <= srr_hat_endyr; yr++) {
        jnll_comp(JNLL_SRR_PENALTY, sp) -= dnorm( log(R(sp, yr)), log(R_hat(sp,yr)) - bias_adjust_proc*square(R_sd(sp))/2.0, R_sd(sp), true);
      }
    }
  }


  // Slot 12 -- Reference point penalties
  // -- CMSY
  Type CMSY = 0;

  // --- Sum terminal catch across species
  if(HCR == 1){

    // -- Loop through catch data
    for(fsh_ind = 0; fsh_ind < catch_ctl.rows(); fsh_ind++){

      flt = catch_ctl(fsh_ind, 0) - 1;
      sp = catch_ctl(fsh_ind, 1) - 1;
      flt_yr = catch_ctl(fsh_ind, 2);

      // Add fishery data from terminal year
      if((flt_type(flt) == 1) && (forecast(sp) == 1) && (estDynamics(sp) == 0)){
        if(flt_yr == projyr){
          CMSY  += catch_hat(fsh_ind);
        }
      }
    }

    jnll_comp(JNLL_REFPT_PENALTY, 0) = - square(CMSY/1000000.0); // CMSY is ll


    // --- Add biomass_depletion constraint
    for(sp = 0; sp < nspp; sp++) {
      if((forecast(sp) == 1) && (estDynamics(sp) == 0)){
        penalty = 0.0;
        Type nothing_useful =  posfun( (ssb_depletion(sp, nyrs-1) - Plimit(sp)), Type(0.0001), penalty); (void) nothing_useful;
        jnll_comp(JNLL_REFPT_PENALTY, sp) += 500.0 * square(CMSY/1000.0) * penalty; // CMSY
      }
    }
  }


  for(sp = 0; sp < nspp; sp++) {
    // -- Single-species static reference points
    if((DynamicHCR == 0) && (forecast(sp) == 1) && (msmMode == 0) && (estDynamics(sp) == 0)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
      }

      // F that acheives \code{Ftarget}% of SSB0 in the end of the projection
      if(HCR == 3){
        // Depletion target: realized SSB in the final projection year against SB0.
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Ftarget_percent(sp));
      }

      // -- SPR
      if(HCR > 3){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRtarget(sp)/SPR0(sp))-Ftarget_percent(sp));
      }
    }

    // -- Single-species dynamic reference points
    if((DynamicHCR == 1) && (forecast(sp) == 1) && (msmMode == 0) && (estDynamics(sp) == 0)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
      }

      for(yr = 1; yr < nyrs; yr++){ // No initial abundance
        // F that acheives Ftarget% of SSB0y
        if(HCR == 3){
          jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((DynamicSBF(sp, yr)/DynamicSB0(sp, yr))-Ftarget_percent(sp));
        }
      }

      // -- SPR
      if(HCR > 3){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRlimit(sp)/SPR0(sp))-Flimit_percent(sp));
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((SPRtarget(sp)/SPR0(sp))-Ftarget_percent(sp));
      }
    }


    // -- Multi-species static reference points (all biomass_depletion based)
    if((forecast(sp) == 1) && (msmMode > 0) && (estDynamics(sp) == 0)){

      // -- Avg F (have F limit)
      if(HCR == 2){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Flimit_percent(sp));
      }

      // F that achieves \code{Ftarget}% of SSB0 in the end of the projection
      if(HCR == 3){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Ftarget_percent(sp));
      }

      // Tiered HCRs with Limit and Targets
      if(HCR == 6){
        jnll_comp(JNLL_REFPT_PENALTY, sp)  += 200*square((ssb(sp, nyrs-1)/SB0(sp, nyrs-1))-Flimit_percent(sp));
      }
    }
  }


  // Slot 13 -- N-at-age < 0 penalty. See "posfun"
  for(sp = 0; sp < nspp; sp++){
    if(estDynamics(sp) == 0){
      jnll_comp(JNLL_ZERO_N_PENALTY, sp) += zero_N_pen(sp);
    }
  }


  // Slots 14-15 -- M PRIORS AND RANDOM EFFECTS
  for(sp = 0; sp < nspp; sp++) {

    // -- Slot 14 -- M priors
    // M1_model 0 = fixed natural mortality from M1_base; 1 = estimate sex- and
    // age-invariant M1_at_age; 2 = sex-specific (two-sex model), age-invariant;
    // 3 = sex- and age-specific.

    // Lognormal prior on M: M_prior is its mean when bias_adjust_proc = 1, its median
    // when 0. The log-scale centre is computed once per species, outside the age loop.
    Type M_prior_mean = log(M_prior(sp)) - bias_adjust_proc*square(M_prior_sd(sp)) / 2.0;

    // Prior on M1_at_age only
    if( (M1_use_prior(sp) == 1) && (M2_use_prior(sp) == 0) ) {
      // Determine loop bounds dynamically based on M1_model
      int nsex_tmp = (M1_model(sp) > 1) ? nsex(sp) : 1;
      int nage_tmp = (M1_model(sp) == 3) ? nages(sp) : 1;

      for(int sex = 0; sex < nsex_tmp; sex++) {
        for(int age = 0; age < nage_tmp; age++) {
          jnll_comp(JNLL_M_PRIOR, sp) -= dnorm(log_M1(sp, sex, age), M_prior_mean, M_prior_sd(sp), true);
        }
      }
    }

    // Prior on total M_at_age (M1_at_age and M2_at_age)
    // FIXME: sex-specific prior?
    if( M2_use_prior(sp) == 1 ) {
      for(int sex = 0; sex < nsex(sp); sex++) {
        for(int age = 0; age < nages(sp); age++) {
          for(int yr = 0; yr < nyrs_hind; yr++) {
            jnll_comp(JNLL_M_PRIOR, sp) -= dnorm(log(M_at_age(sp, sex, age, yr)), M_prior_mean, M_prior_sd(sp), true);
          }
        }
      }
    }


    // -- Slot 15 -- M random effects
    // Applied per species at M1_model = 1, per species AND sex at M1_model = 2.
    // Variance and correlation coefficients are species-specific, sex-invariant.
    // - M1_re = 0: none (default).
    // - M1_re = 1: IID by age, constant over years.
    // - M1_re = 2: IID by year, constant over ages.
    // - M1_re = 3: IID by year and age.
    // - M1_re = 4: AR1 by age, constant over years.
    // - M1_re = 5: AR1 by year, constant over ages.
    // - M1_re = 6: 2D-AR1 by year and age.

    // Determine how many sexes to run REs for (only add 2nd sex if M1_model == 2)
    int num_re_sexes = (M1_model(sp) == 2 && nsex(sp) > 1) ? 2 : 1;

    // M1_re = 1/4: Random effects varies by age (IID or AR1) and constant over years.
    if( (M1_re(sp) == 1) || (M1_re(sp) == 4) ) {
      Type sigma_M = exp(M1_dev_log_sd(sp, 0));
      Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
      Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_a, 2)), 0.5);
      vector<Type> M_re_age(nages(sp));

      for(int sex = 0; sex < num_re_sexes; sex++) {
        for(int age = 0; age < nages(sp); age++) {
          M_re_age(age) = log_M1_dev(sp, sex, age, 0);
        }
        jnll_comp(JNLL_M_RE, sp) += SCALE(AR1(rho_M_a), Sigma_M)(M_re_age);
      }
    }

    // M1_re = 2/5: Random effects varies by year (IID or AR1) and constant over ages.
    if( (M1_re(sp) == 2) || (M1_re(sp) == 5) ) {
      Type sigma_M = exp(M1_dev_log_sd(sp, 0));
      Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
      Type Sigma_M = pow(pow(sigma_M, 2) / (1.0 - pow(rho_M_y, 2)), 0.5);
      vector<Type> M_re_yr(nyrs_hind);

      for(int sex = 0; sex < num_re_sexes; sex++) {
        for(int yr = 0; yr < nyrs_hind; yr++) {
          M_re_yr(yr) = log_M1_dev(sp, sex, 0, yr);
        }
        jnll_comp(JNLL_M_RE, sp) += SCALE(AR1(rho_M_y), Sigma_M)(M_re_yr);
      }
    }

    // M1_re = 3/6: Random effects varies by age and year (IID or 2D-AR1)
    if( (M1_re(sp) == 3) || (M1_re(sp) == 6) ) {
      Type sigma_M = exp(M1_dev_log_sd(sp, 0));
      Type rho_M_a = rho_trans(M1_rho(sp, 0, 0));
      Type rho_M_y = rho_trans(M1_rho(sp, 0, 1));
      Type Sigma_M = pow(pow(sigma_M, 2) / ((1.0 - pow(rho_M_y, 2)) * (1.0 - pow(rho_M_a, 2))), 0.5);
      array<Type> M_re_a_yr(nages(sp), nyrs_hind);

      for(int sex = 0; sex < num_re_sexes; sex++) {
        for(int age = 0; age < nages(sp); age++) {
          for(int yr = 0; yr < nyrs_hind; yr++) {
            M_re_a_yr(age, yr) = log_M1_dev(sp, sex, age, yr);
          }
        }
        // SEPARABLE(f, g) applies f to the OUTERMOST array dimension and g to
        // the fastest-running one (TMB density.hpp). The array is
        // (age, year), so the year correlation goes first and the age
        // correlation second -- written the other way round, rho_M_a would
        // correlate over years and rho_M_y over ages.
        jnll_comp(JNLL_M_RE, sp) += SCALE(SEPARABLE(AR1(rho_M_y), AR1(rho_M_a)), Sigma_M)(M_re_a_yr);
      }
    }
  }


  // Slot 19 -- Linkage-table priors on the natural scale.
  // A prior is specified on the natural-scale parameter it names (R0, alpha, M1,
  // K, ...), so an (Intercept) row standing in for a log-scale base parameter is
  // transformed first. On a slope row b_nat == b, the coefficient already being
  // on the natural / linear scale. Families:
  //   0 = none
  //   1 = normal    -- dnorm(b_nat, p1, p2)
  //   2 = lognormal -- dnorm(log(b_nat), p1 - bias_adjust_proc*p2^2/2, p2);
  //                    exp(p1) = mean if bias_adjust_proc = 1, median if 0
  //   3 = gamma     -- dgamma(b_nat, p1, 1/p2)
  //   4 = beta      -- dbeta(b_nat, p1, p2)
  //
  // Most (Intercept) rows are mapped out (beta_linkage(i) stays at 0), so their
  // prior is evaluated against the *base parameter* (`rec_pars`, `log_M1`,
  // `log_growth_pars`) instead. The recruitment `R_init` intercept has no base
  // parameter and stays estimable, so its prior is read on beta_linkage(i).
  for (int i = 0; i < beta_linkage.size(); ++i) {
    int fam = linkage_prior_family(i);
    if (fam == 0) continue;
    Type p1  = linkage_prior_p1(i);
    Type p2  = linkage_prior_p2(i);
    int slot_col = linkage_species(i) - 1;
    if (slot_col < 0) slot_col = 0;        // shared/all => column 0
    if (slot_col >= n_col) slot_col = 0;

    Type b = beta_linkage(i);
    // Whether the intercept's base parameter lives on the log scale (so b_nat = exp(b)).
    // True for every log-scale base (rec_pars, log_M1, log_growth_pars, index_log_q,
    // comp_weights); set false for a natural-scale base (the selectivity inflection
    // sel_inf), whose prior is read on the natural scale directly.
    bool base_is_log = true;
    if (linkage_is_intercept(i) == 1) {
      // Re-target the prior to the base parameter that this intercept row stands in for.
      // Sentinel 0 in species/sex/age_bin means "all levels"; pick the first concrete
      // cell as the prior target (in practice intercepts are stratified per
      // species/sex/age and there is a single base-parameter scalar per row).
      int sp_in = linkage_species(i);
      int sx_in = linkage_sex(i);
      int ab_in = linkage_age_bin(i);
      int sp_idx = (sp_in == 0) ? 0 : (sp_in - 1);
      int sx_idx = (sx_in == 0) ? 0 : (sx_in - 1);
      int ab_idx = (ab_in == 0) ? 0 : (ab_in - 1);

      int proc  = linkage_process(i);
      int param = linkage_param(i);
      if (proc == RCEATTLE_PROC_RECRUIT) {
        // recruitment params: R0=0, alpha=1, beta=2 -> rec_pars cols 0..2 (stored as
        // log_R0, log_alpha, log_beta on the log scale). R_init (3) has NO rec_pars
        // column -- its intercept IS the level, so the prior stays on beta_linkage(i)
        // and reading rec_pars here would run off the end of the matrix.
        if (param != RCEATTLE_REC_R_INIT) {
          b = rec_pars(sp_idx, param);
        }
      } else if (proc == RCEATTLE_PROC_M) {
        b = log_M1(sp_idx, sx_idx, ab_idx);
      } else if (proc == RCEATTLE_PROC_GROWTH) {
        // growth params: K=0, L1=1, Linf=2, m=3 -> log_growth_pars(sp, sx, param)
        // SD endpoints: sd_L1=4, sd_Linf=5 -> growth_log_sd(sp, sx, param - N_GROWTH)
        if (param < RCEATTLE_N_GROWTH_PARAMS) {
          b = log_growth_pars(sp_idx, sx_idx, param);
        } else {
          b = growth_log_sd(sp_idx, sx_idx, param - RCEATTLE_N_GROWTH_PARAMS);
        }
      } else if (proc == RCEATTLE_PROC_Q) {
        // Catchability is fleet- not species-indexed. The intercept
        // stands in for the base log-catchability index_log_q(fleet).
        int fl_in  = linkage_fleet(i);
        int fl_idx = (fl_in == 0) ? 0 : (fl_in - 1);
        b = index_log_q(fl_idx);
      } else if (proc == RCEATTLE_PROC_COMP) {
        // Dirichlet-multinomial overdispersion. comp/caal weights are fleet-indexed,
        // diet weights predator(species)-indexed; all stored on the log scale
        // (theta = exp(weight)), so b_nat = exp(b) is the natural DM scalar the prior
        // targets. param: theta_comp=0, theta_caal=1, theta_diet=2.
        int fl_in  = linkage_fleet(i);
        int fl_idx = (fl_in == 0) ? 0 : (fl_in - 1);
        if (param == 0)      b = comp_weights(fl_idx);
        else if (param == 1) b = caal_weights(fl_idx);
        else if (param == 2) b = diet_comp_weights(sp_idx);
      } else if (proc == RCEATTLE_PROC_SEL) {
        // Parametric selectivity base parameters are (slot, fleet, sex), with asc -> slot
        // 0 and desc -> slot 1. param: slp_asc=0 / slp_desc=1 -> log_sel_slp (log scale,
        // exp -> slope); inf_asc=2 / inf_desc=3 -> sel_inf (natural scale, the inflection
        // age / peak). The (Intercept) beta_linkage row is pinned at 0, so it adds no
        // offset to selectivity: the base parameter carries the mean and this prior
        // regularizes it (Cole Monnahan's GOA pollock priors -- lognormal on the
        // log-slope, normal on the natural inflection). coff (param 4) has no scalar base
        // parameter and cannot carry a linkage (.check_sel_linkage_support).
        int fl_in  = linkage_fleet(i);
        int fl_idx = (fl_in == 0) ? 0 : (fl_in - 1);
        int slot   = param % 2;   // asc -> 0, desc -> 1
        if (param == 0 || param == 1) {
          b = log_sel_slp(slot, fl_idx, sx_idx);
        } else if (param == 2 || param == 3) {
          b = sel_inf(slot, fl_idx, sx_idx);
          base_is_log = false;    // inflection is natural-scale, not log
        } else if (param == 5) {
          // apical: the per-sex multiplier, stored logged, so the prior reads
          // on the multiplier itself (lognormal centred on 1 = no offset).
          b = log_sel_apical(fl_idx, sx_idx);
        } else if (param >= 6 && param <= 11) {
          // DoubleNormalSS3: the prior reads the parameter as stored, on SS3's
          // own scale (peak in cm, P3/P4 log widths, P2/P5/P6 logits), which
          // is where an SS3 control file puts its priors.
          b = sel_dn6(param - 6, fl_idx, sx_idx);
          base_is_log = false;
        }
      }
    }

    // Back-transform to natural scale when this is an intercept row on a log-scale base
    // parameter (b holds the log-scale parameter). For all slope rows, and for a
    // natural-scale intercept base (sel_inf), b_nat == b already.
    Type b_nat = (linkage_is_intercept(i) == 1 && base_is_log) ? exp(b) : b;

    if (fam == 1) {                         // normal(p1, p2) on natural scale
      jnll_comp(JNLL_LINKAGE_PRIOR, slot_col)            -= dnorm(b_nat, p1, p2, true);
      unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, slot_col) -= dnorm(b_nat, p1, p2, true);
    } else if (fam == 2) {                  // lognormal on log(b_nat); exp(p1) = mean when bias_adjust_proc = 1, median when 0
      // For a log-scale intercept base: log(b_nat) = b (efficient form avoids
      // log(exp(b))). A natural-scale intercept base (sel_inf) takes log(b_nat).
      Type log_b_nat = (linkage_is_intercept(i) == 1 && base_is_log) ? b : log(b_nat);
      Type mu_log    = p1 - bias_adjust_proc * square(p2) / 2.0;
      jnll_comp(JNLL_LINKAGE_PRIOR, slot_col)            -= dnorm(log_b_nat, mu_log, p2, true);
      unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, slot_col) -= dnorm(log_b_nat, mu_log, p2, true);
    } else if (fam == 3) {                  // gamma(p1=shape, p2=rate) on natural scale
      jnll_comp(JNLL_LINKAGE_PRIOR, slot_col)            -= dgamma(b_nat, p1, Type(1.0)/p2, true);
      unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, slot_col) -= dgamma(b_nat, p1, Type(1.0)/p2, true);
    } else if (fam == 4) {                  // beta(p1=shape1, p2=shape2) on natural scale
      jnll_comp(JNLL_LINKAGE_PRIOR, slot_col)            -= dbeta(b_nat, p1, p2, true);
      unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, slot_col) -= dbeta(b_nat, p1, p2, true);
    }
  }


  // 13.2. Diet likelihood components
  // Simulated stomach contents go to this copy, not into the REPORTed diet_obs
  // -- see the _sim naming rule in section 5.2. Stomachs the model does not fit
  // keep the observed values they were copied with.
  matrix<Type> diet_sim = diet_obs;

  if((msmMode > 2) || (imax(suitMode) > 0)) {

    int current_j = 0; // Track position in diet_ctl

    for (int i = 0; i < n_stomach_obs; ++i) {

      // Take this stomach's prey as the run of diet_ctl rows with stomach_id == i. The
      // cursor only moves forward, so the ids must be sorted 0, 1, 2, ... with no gaps:
      // clean_data() rebuilds stomach_id from the pred/sex/age/year stratum and sorts by
      // it on every pass, and data_check() rejects a table that breaks it. Out of that
      // order the cursor runs past whole stomachs, which then drop out of the likelihood
      // silently, and build_osa_data(), matching on which(stomach_id == i), would
      // disagree with this scan.
      int start_j = current_j;
      while((current_j < diet_ctl.rows()) && (stomach_id(current_j) == i)) {
        current_j++;
      }

      // --- Process the predator if suitability is estimated and data are available
      // Test the empty run first: start_j can sit one past the last row, where
      // reading diet_ctl / diet_obs would be out of bounds.
      int n_prey = current_j - start_j;
      if (n_prey == 0) continue;

      int rsp = diet_ctl(start_j, 0) - 1;
      Type N_s = diet_obs(start_j, 0);

      if (suitMode(rsp) <= 0) continue;

      // -- Pre-allocate TMB vectors with space for "Other prey" (+1)
      vector<Type> obs_diet_prop(n_prey + 1); obs_diet_prop.setZero();
      vector<Type> pred_diet_prop(n_prey + 1); pred_diet_prop.setZero();

      for (int k = 0; k < n_prey; ++k) {
        obs_diet_prop(k) = diet_obs(start_j + k, 1);
        pred_diet_prop(k) = diet_hat(start_j + k, 1);
      }

      // --- Add in "Other prey" ---
      Type sum_obs_p = obs_diet_prop.head(n_prey).sum();
      if (sum_obs_p > 1.0) sum_obs_p = 1.0;
      obs_diet_prop(n_prey) = 1.0 - sum_obs_p;

      Type sum_est_p = pred_diet_prop.head(n_prey).sum();
      pred_diet_prop(n_prey) = posfun(1.0 - sum_est_p, Type(0.00001), penalty); // Making it differentiable (cant do if statement)

      // Vectorized Offset & Normalization. Uses the same configurable proportion
      // offset (DATA_SCALAR comp_offset, default 1e-5) as the age/length comps,
      // so fitting and the diet OSA obsvec stay consistent.
      obs_diet_prop += comp_offset;
      pred_diet_prop += comp_offset;

      obs_diet_prop /= obs_diet_prop.sum();
      pred_diet_prop /= pred_diet_prop.sum();

      // Convert observed prop to observed numbers
      vector<Type> obs_diet_content = obs_diet_prop * N_s;

      Type stomach_log_likelihood = 0;
      Type unweighted_stomach_log_likelihood = 0;
      Type DM_diet_par = DM_diet_pars(rsp);

      // Calculate alpha parameters.
      vector<Type> diet_alphas = pred_diet_prop * N_s * DM_diet_par;
      vector<Type> unweighted_diet_alphas = pred_diet_prop * N_s; // For "unweighted" likelihood (DM_diet_par = 1)

      // Likelihood
      if(osa_mode == 0){
        // ---- Normal fitting: weighted density read from diet_obs (unchanged) ----
        switch(diet_ll_type(rsp)){

        case 0:  // Full multinomial
          stomach_log_likelihood = dmultinom(obs_diet_content, pred_diet_prop, true);

          unweighted_jnll_comp(JNLL_STOMACH, rsp) -= stomach_log_likelihood;
          jnll_comp(JNLL_STOMACH, rsp) -= diet_comp_weights(rsp) * stomach_log_likelihood;
          break;
        case 1:  // Dirichlet-multinomial
          // Calculate the log-likelihood
          stomach_log_likelihood = ddirmultinom(obs_diet_content, diet_alphas, true);
          unweighted_stomach_log_likelihood = ddirmultinom(obs_diet_content, unweighted_diet_alphas, true);

          unweighted_jnll_comp(JNLL_STOMACH, rsp) -= unweighted_stomach_log_likelihood;
          jnll_comp(JNLL_STOMACH, rsp) -= stomach_log_likelihood;
          break;

        default:
          error("Invalid 'diet_ll_type'");
        }
      } else {
        // ---- OSA build: unweighted, keep-gated conditional density read from
        // obsvec. The diet composition for a stomach is its prey items plus an
        // "other prey" category (the last bin, fixed by sum-to-1 and dropped).
        // diet_obsvec_idx gives the obsvec start for this stomach. ----
        int start = diet_obsvec_idx(i);
        if(start >= 0){
          vector<Type> osa_x = obsvec.segment(start, n_prey + 1);
          if(diet_ll_type(rsp) == 1){   // Dirichlet-multinomial (fitted DM par)
            // No CDF term: the beta-binomial has no elementary one (comp_osa.hpp).
            jnll_comp(JNLL_STOMACH, rsp) -= ddirmultinom_osa(osa_x, diet_alphas, keep.segment(start, n_prey + 1), 1, 1);
          } else {                      // multinomial
            jnll_comp(JNLL_STOMACH, rsp) -= dmultinom_osa(osa_x, pred_diet_prop, keep.segment(start, n_prey + 1), 1, 1, osa_mode == 2);
          }
        }
      }

      // -- Simulate the stomach contents, for sim_mod(simulate = TRUE) --
      // Drawn under the same family and from the same predicted proportions and
      // concentration the density above uses. The composition is this stomach's
      // prey plus the "other prey" balance, so the draw covers n_prey + 1 bins;
      // only the prey bins are written back, "other prey" being whatever is left
      // over and recomputed on the next fit, so storing it would double-count.
      SIMULATE {
        // Draw from the composition BEFORE the offset, rebuilt here from diet_hat so
        // ordinary fitting pays nothing for it. The offset keeps log(0) out of the
        // density rather than describing the stomachs, and the transform is affine, so
        // drawing from the un-offset prediction keeps a simulate-then-refit round trip
        // unbiased.
        vector<Type> sim_p(n_prey + 1);
        Type prey_sum = 0.0;
        for(int k = 0; k < n_prey; k++){
          sim_p(k) = diet_hat(start_j + k, 1);
          prey_sum += sim_p(k);
        }
        sim_p(n_prey) = (prey_sum < Type(1.0)) ? Type(1.0) - prey_sum : Type(0.0);
        Type sim_tot = sim_p.sum();
        if(sim_tot > Type(0)) sim_p /= sim_tot;

        if((N_s > Type(0)) && (sim_tot > Type(0))){
          vector<Type> sim_counts(n_prey + 1);
          switch(diet_ll_type(rsp)){
          case 0: {
            // diet_comp_weights multiplies the multinomial density, making it an
            // effective sample size: weighting a multinomial log-likelihood by w is, up
            // to a constant, a multinomial at w*N. Draw at that size, or a down-weighted
            // fleet gets stomachs the likelihood treats as w times less informative than
            // they are. For the Dirichlet-multinomial the same parameter is log theta and
            // already enters through the concentration below.
            Type n_eff = N_s * diet_comp_weights(rsp);
            sim_counts = rmultinom_rce(n_eff, sim_p);
            break;
          }
          case 1: {
            vector<Type> sim_alphas = sim_p * N_s * DM_diet_par;
            sim_counts = rdirmultinom_rce(N_s, sim_alphas);
            break;
          }
          default:
            error("Invalid 'diet_ll_type'");
          }
          // Divide by the number of observations actually PLACED, not by N_s:
          // rmultinom_rce rounds N_s to a whole number, so a fractional sample
          // size would leave the prey proportions summing to round(N_s)/N_s > 1,
          // which data_check() rejects and self_test() then reports as a failed
          // replicate rather than as bad data.
          //
          // A stomach whose weighted sample size rounds to zero places nothing
          // and comes back EMPTY rather than carrying its observed proportions
          // over, the same rule the slot 2 composition draw follows: keeping
          // them would hand the refit the real diet under a weight saying the
          // stomach is worth less than one observation.
          Type n_drawn = sim_counts.sum();
          for(int k = 0; k < n_prey; k++){
            diet_sim(start_j + k, 1) = (n_drawn > Type(0))
              ? sim_counts(k) / n_drawn
              : Type(0);
          }
        }
      }
    }
  }

  // Reported under a _sim name -- see the naming rule in section 5.2. Equal to
  // diet_obs for a model with no diet likelihood, where nothing above is drawn.
  SIMULATE {
    matrix<Type> diet_obs_sim = diet_sim;
    REPORT(diet_obs_sim);
  }
  // Paste unweighted likelihood parts over
  // - Only comp and diet comp use data weights!
  for(flt = 0; flt < n_col; flt++){
    for(int ind = 0; ind < 19; ind++){
      if((ind!= 2) & (ind != 18)){
        unweighted_jnll_comp(ind, flt) = jnll_comp(ind, flt);
      }
    }
  }


  /** ------------------------------------------------------------------------ //
   // 14. REPORT SECTION                                                        //
   * ------------------------------------------------------------------------- */

  // 14.1. Population components
  REPORT( pop_scalar );
  ADREPORT( pop_scalar );
  REPORT( avg_R );
  REPORT( Z_at_age );
  REPORT( N_at_age );
  REPORT( avgN_at_age );
  REPORT( biomass_at_age );
  REPORT( biomass_depletion );
  REPORT( ssb_depletion );
  REPORT( biomass );
  REPORT( ssb );
  REPORT( exploitable_biomass );
  REPORT(R_sd);
  REPORT( R0 );
  REPORT( R_init );
  REPORT( R );
  REPORT( R_hat );
  REPORT( M_at_age );
  //REPORT( mature_females );
  //REPORT( sex_ratio_hat );
  //REPORT( N_at_age_dB0);
  //REPORT( N_at_age_dBF);
  // ADREPORT( M_at_age );            // no SEs: see ?plot_mortality
  REPORT( weight_hat );
  REPORT( growth_matrix );
  REPORT( length_hat );

  ADREPORT( biomass );
  ADREPORT( ssb );

  // -- Secondary display series, so add_ci = TRUE can draw an interval at all.
  // Natural scale ONLY: exploitable_biomass sums over fisheries alone, so it is
  // exactly 0 for a survey-only species or proj_F_prop = 0, and the depletions
  // divide by B0 / SB0. log() of either puts -Inf on the tape and NaNs the whole
  // sdreport; the plotters recover sd(log x) = sd(x)/x instead.
  ADREPORT( exploitable_biomass );
  ADREPORT( biomass_depletion );
  ADREPORT( ssb_depletion );

  // -- Log-scale derived series.
  // sdreport linearizes once about the MLE, so its SD is exact only for a linear function
  // of the parameters. These are built multiplicatively (R = R0 * exp(rec_dev); n-at-age
  // is a product of survivals), so log(x) is near-linear and exp(log(x) +/- z * sd) is
  // right-skewed and cannot cross zero. The log SD is also the CV, the form the ABC / OFL
  // buffer wants. The natural-scale ADREPORTs stay so existing sdrep$value callers work.
  matrix<Type>  log_biomass = biomass;  log_biomass = log(biomass.array());
  matrix<Type>  log_ssb     = ssb;      log_ssb     = log(ssb.array());
  matrix<Type>  log_R       = R;        log_R       = log(R.array());
  REPORT( log_biomass );
  REPORT( log_ssb );
  REPORT( log_R );
  ADREPORT( log_biomass );
  ADREPORT( log_ssb );
  ADREPORT( log_R );
  ADREPORT( R_sd );
  ADREPORT( R );

  // -- Log selectivity-at-age, for a delta-method interval on the curve. Log and not
  // logit, because the non-parametric forms normalize to a MEAN of one, so sel_at_age
  // routinely exceeds one. Gated because the delta method costs every reported value
  // against every parameter; see ?fit_control.
  //
  // Only estimated, age-based lead fleets, from their first selected bin: those cells
  // cannot be a structural zero, and one log(0) = -Inf on the tape turns EVERY sdreport
  // quantity to NaN, biomass and SSB included. The four routes to a structural zero are
  // enumerated in ?fit_control and all read from the data, so the reported SET never
  // depends on a parameter value. A reported VALUE can still underflow -- a narrow dome
  // puts exp(-x^2/2) below 1e-308 many bins from its peak -- so the log is floored with
  // max2(), exact above the floor. log_sel_at_age_index carries each row's fleet, sex,
  // absolute age and calendar year, since which cells are reported depends on the model.
  if(adreport_sel == 1){
    // First bin and bin count per fleet. data_check() bounds Bin_first_selected
    // by nages, but a 0 in the workbook still arrives here as -1, and the model
    // builds safebounds = FALSE, so both ends are clamped before indexing.
    vector<int> sel_first_bin(n_flt); sel_first_bin.setZero();
    vector<int> sel_rep_bins(n_flt);  sel_rep_bins.setZero();
    for(int flt = 0; flt < n_flt; flt++){
      if((flt_type(flt) <= 0) || (flt_sel_lead(flt) != 1) ||
         (flt_sel_type(flt) == 0) || (flt_sel_dim(flt) == 1)) continue;
      int n_age = nages(flt_spp(flt));
      int first = bin_first_selected(flt);
      if(first < 0) first = 0;
      if(first > n_age) first = n_age;
      sel_first_bin(flt) = first;
      sel_rep_bins(flt)  = n_age - first;
    }

    int n_sel = 0;
    for(int flt = 0; flt < n_flt; flt++){
      n_sel += nsex(flt_spp(flt)) * sel_rep_bins(flt) * nyrs_hind;
    }

    vector<Type> log_sel_at_age(n_sel); log_sel_at_age.setZero();
    matrix<Type> log_sel_at_age_index(n_sel, 4); log_sel_at_age_index.setZero();

    int i_sel = 0;
    for(int flt = 0; flt < n_flt; flt++){
      if(sel_rep_bins(flt) == 0) continue;
      int sel_sp = flt_spp(flt);
      for(int sel_sex = 0; sel_sex < nsex(sel_sp); sel_sex++){
        for(int sel_age = sel_first_bin(flt); sel_age < nages(sel_sp); sel_age++){
          for(int sel_yr = 0; sel_yr < nyrs_hind; sel_yr++){
            log_sel_at_age(i_sel) = log(max2(sel_at_age(flt, sel_sex, sel_age, sel_yr),
                                             Type(1e-300)));
            log_sel_at_age_index(i_sel, 0) = Type(flt + 1);
            log_sel_at_age_index(i_sel, 1) = Type(sel_sex + 1);
            log_sel_at_age_index(i_sel, 2) = Type(minage(sel_sp) + sel_age);
            log_sel_at_age_index(i_sel, 3) = Type(styr + sel_yr);
            i_sel++;
          }
        }
      }
    }

    REPORT( log_sel_at_age );
    REPORT( log_sel_at_age_index );
    ADREPORT( log_sel_at_age );
  }


  // -- 14.2. Biological reference points
  REPORT( NByageF );
  REPORT( NByage0);
  REPORT( NbyageSPR);
  REPORT( B0 );
  REPORT( SB0 );
  REPORT( SBF );
  REPORT( DynamicB0 );
  REPORT( DynamicSB0 );
  REPORT( DynamicSBF );
  REPORT( SPR0 );
  REPORT( SPRFinit );
  REPORT( SPRlimit );
  REPORT( SPRtarget );
  REPORT( steepness );
  REPORT( proj_F );
  REPORT( Ftarget );
  REPORT( Flimit );

  // -- 14.3. Selectivity
  REPORT( sel_at_age );
  REPORT( sel_at_length );


  // -- 14.4. Survey components
  REPORT( index_hat );
  REPORT( index_sd );
  REPORT( index_q );
  vector<Type>  log_index_hat = index_hat;  log_index_hat = log(index_hat.array());// Log predicted index
  REPORT( log_index_hat );
  ADREPORT( log_index_hat );

  // -- 14.4b. Formula-linkage effect sizes: the estimated coefficients, the random-effect
  // deviations and the Rogers (2024) QAR1 effect size, so a user can read the *effect*
  // and, through ADREPORT, its uncertainty -- not just the realized parameter (index_q
  // above). All may be length 0 (no linkages), in which case they report nothing.
  REPORT( beta_linkage );
  ADREPORT( beta_linkage );
  REPORT( beta_linkage_re );
  REPORT( beta_linkage_re_pen );
  // Slot-ordered deviations across both vectors, so a caller can read a group's
  // walk without knowing how it is stored.
  REPORT( beta_linkage_re_all );
  REPORT( beta_linkage_obs );
  ADREPORT( beta_linkage_obs );


  // -- 14.5. Fishery components
  REPORT( F_spp );
  REPORT( F_flt );
  REPORT( F_flt_age );
  REPORT( F_at_age );
  REPORT( catch_hat );
  REPORT( max_catch_hat );
  REPORT( catch_sd );
  REPORT( catch_analytical_sd );


  // 14.6. Age/length composition
  REPORT( age_obs_hat );
  REPORT( age_hat );
  REPORT( comp_obs );
  REPORT( comp_hat );
  REPORT( caal_hat );
  REPORT( equil_catch_hat );
  REPORT( caal_obs );
  REPORT( pred_CAAL );


  // -- 14.6a. Random-effect SD priors (jnll_comp row 19). A prior on a group's deviation
  // SD, routed from linkage_spec(priors = list(sigma = ...)), once per group because the
  // loop is over groups. Shares the linkage-prior row and the same families as the
  // fixed-beta prior loop above, on the natural-scale SD exp(log_sigma).
  // FIXME(jacobian): the normal/gamma/beta families place the prior on the natural-scale
  // SD while the free parameter is log_sigma, WITHOUT the change-of-variables Jacobian
  // d(sd)/d(log_sigma) = sd -- a penalized-likelihood prior, not a proper Bayesian one,
  // matching the fixed-beta loop's convention. `lognormal` is exempt, dnorm(log(sd), ..)
  // having a constant Jacobian in log_sigma. Add `+ log_sigma` for a proper density.
  if (log_sigma_linkage.size() > 0) {
    for (int g = 0; g < log_sigma_linkage.size(); ++g) {
      int fam = linkage_re_sigma_prior_family(g);
      if (fam == 0) continue;
      Type sd = exp(log_sigma_linkage(g));
      Type p1 = linkage_re_sigma_prior_p1(g);
      Type p2 = linkage_re_sigma_prior_p2(g);
      if (fam == 1) {                         // normal(p1, p2) on the SD
        jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= dnorm(sd, p1, p2, true);
        unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= dnorm(sd, p1, p2, true);
      } else if (fam == 2) {                  // lognormal on the SD; exp(p1) = mean when bias_adjust_proc = 1, median when 0
        Type mu_log = p1 - bias_adjust_proc * square(p2) / 2.0;
        jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= dnorm(log(sd), mu_log, p2, true);
        unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= dnorm(log(sd), mu_log, p2, true);
      } else if (fam == 3) {                  // gamma(shape, rate)
        jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= dgamma(sd, p1, Type(1.0)/p2, true);
        unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= dgamma(sd, p1, Type(1.0)/p2, true);
      } else if (fam == 4) {                  // beta(shape1, shape2)
        jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= dbeta(sd, p1, p2, true);
        unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= dbeta(sd, p1, p2, true);
      }
    }
  }

  // -- 14.6a'. Random-effect AR1 correlation priors (jnll_comp row 19), routed from
  // linkage_spec(priors = list(rho = ...)) on an ar1 group. On the natural (-1, 1)
  // correlation rho = rho_trans(trans_rho_linkage): normal(p1, p2) directly, or
  // beta(p1, p2) on (rho + 1) / 2. Once per group.
  // FIXME(jacobian): as with the sigma priors above, the prior is on the natural rho
  // while the free parameter is trans_rho_linkage, without the rho_trans Jacobian -- a
  // penalized-likelihood prior, not a proper density on rho.
  if (trans_rho_linkage.size() > 0) {
    for (int g = 0; g < linkage_re_rho.size(); ++g) {
      int fam = linkage_re_rho_prior_family(g);
      if (fam == 0) continue;
      Type rho = rho_trans(trans_rho_linkage(linkage_re_rho(g)));
      Type p1 = linkage_re_rho_prior_p1(g);
      Type p2 = linkage_re_rho_prior_p2(g);
      if (fam == 1) {                         // normal(p1, p2) on rho
        jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= dnorm(rho, p1, p2, true);
        unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= dnorm(rho, p1, p2, true);
      } else if (fam == 4) {                  // beta(p1, p2) on (rho + 1) / 2
        jnll_comp(JNLL_LINKAGE_PRIOR, 0)            -= dbeta((rho + Type(1)) / Type(2), p1, p2, true);
        unweighted_jnll_comp(JNLL_LINKAGE_PRIOR, 0) -= dbeta((rho + Type(1)) / Type(2), p1, p2, true);
      }
    }
  }

  // -- 14.6b. Random-effect linkage density (jnll_comp row 20)
  // Group-oriented so each covariance structure couples its deviations correctly. Per RE
  // group, gather its slot-space deviations in ascending slot order -- which the registry
  // assigns in elapsed-time order -- then dispatch on the group's structure:
  //   us  (IID)        : sum of N(0, sigma) densities (the index_q_dev idiom);
  //   rw  (RandomWalk) : N(0, sigma) on successive first differences, the first deviate
  //                      pinned at 0 by the map (as index_q_dev);
  //   ar1 (AR1)        : stationary AR1 with MARGINAL SD sigma and correlation
  //                      rho = rho_trans(trans_rho_linkage), the glmmTMB / Rogers-QAR1
  //                      convention SCALE(AR1(rho), sigma).
  // Guarded on size so row 20 stays exactly 0 without a random linkage, and placed before
  // REPORT(jnll_comp) so the reported matrix carries it.

  // The DSEM GMRF density, computed back in section 5.5b. Routed through
  // jnll_comp rather than added straight to jnll so it appears in the reported
  // breakdown like every other component; JNLL_DSEM stays exactly 0 without one.
  // Column 0 because the density is joint over all species, not separable.
  //
  // dsem_on, NOT dsem_drives_rec: jnll_dsem is the SEM's whole negative log
  // likelihood -- the GMRF over the latent field AND the measurement density on
  // any observed covariate cell -- which a SEM has whatever process it drives.
  // Gating it on dsem_drives_rec would drop both from the objective for a SEM
  // that drives something else, leaving the field declared random with nothing
  // scoring it.
  if (dsem_on == 1) {
    jnll_comp(JNLL_DSEM, 0) += jnll_dsem;
    // The unweighted copy is built from rows < 19 only, so mirror this one too:
    // a diagnostic reading unweighted_jnll_comp as "the objective without data
    // weights" would otherwise be short by the whole DSEM density.
    unweighted_jnll_comp(JNLL_DSEM, 0) += jnll_dsem;
  }
  REPORT(jnll_dsem);
  // Reported so the bias correction can be read, inverted and tested from R.
  // retrospective() needs it to write a target recruitment deviation into
  // dsem_x_tj: the template subtracts margvar/2 when it derives rec_dev, so a
  // caller writing the deviation directly lands that much low. A quota-relevant
  // correction that cannot be inspected from R cannot be checked.
  REPORT(dsem_margvar_tj);
  REPORT(dsem_Q);
  REPORT(dsem_xhat_tj);
  REPORT(dsem_delta_tj);
  // The latent field the model actually used, after any rank-reduced or
  // conditional-kriging option has mapped the parameter block onto it. rec_dev
  // is derived from this, so it is what says whether a substituted draw reached
  // the dynamics; without it a caller can only infer that from rec_dev.
  REPORT(dsem_z_tj);
  // The covariate observations as simulated, and the cells the draw touched.
  // Empty for a model with no DSEM, and all-zero mask for one whose covariates
  // are family = "fixed" -- those carry no measurement density to draw from.
  REPORT(dsem_y_tj_sim);
  REPORT(dsem_y_tj_drawn_sim);

  if (log_sigma_linkage.size() > 0) {
    int n_re = beta_linkage_re_all.size();
    for (int grp = 0; grp < log_sigma_linkage.size(); ++grp) {
      Type sigma = exp(log_sigma_linkage(grp));
      // Collect this group's slots in slot (= time) order.
      int len = 0;
      for (int g = 0; g < n_re; ++g) if (linkage_re_sigma(g) == grp) len++;
      if (len == 0) continue;
      vector<Type> re(len), obs(len);
      vector<int> obs_mask(len), obs_slot(len);
      int j = 0;
      for (int g = 0; g < n_re; ++g) if (linkage_re_sigma(g) == grp) {
        obs(j) = linkage_re_obs_value(g); re(j) = beta_linkage_re_all(g);
        obs_mask(j) = linkage_re_obs_mask(g); obs_slot(j) = g; j++;
      }

      int st = linkage_re_struct(grp);
      if (st == 1) {                                   // rw / RandomWalk
        for (int t = 1; t < len; ++t) {
          jnll_comp(JNLL_LINKAGE_RE, 0)            -= dnorm(re(t) - re(t - 1), Type(0), sigma, true);
          unweighted_jnll_comp(JNLL_LINKAGE_RE, 0) -= dnorm(re(t) - re(t - 1), Type(0), sigma, true);
        }
      } else if (st == 2) {                            // ar1 / first-order AR
        // Stationary AR1 with marginal SD sigma and correlation rho. SCALE(...)
        // returns the negative log density, so it is ADDED (unlike the -= dnorm
        // forms). Reduces to the IID sum at rho = 0.
        Type rho = rho_trans(trans_rho_linkage(linkage_re_rho(grp)));
        jnll_comp(JNLL_LINKAGE_RE, 0)            += SCALE(AR1(rho), sigma)(re);
        unweighted_jnll_comp(JNLL_LINKAGE_RE, 0) += SCALE(AR1(rho), sigma)(re);
      } else {                                         // us / IID (default)
        jnll_comp(JNLL_LINKAGE_RE, 0)            -= sum(dnorm(re, Type(0), sigma, true));
        unweighted_jnll_comp(JNLL_LINKAGE_RE, 0) -= sum(dnorm(re, Type(0), sigma, true));
      }

      // Rogers QAR1 observation: the ar1 latent re(t) is measured as
      // linkage_re_obs_value(t) with SD exp(log_obs_sd_linkage(obs_slot)), one ESTIMATED
      // SD per observed group (matching the reference Estimate_q = 6 / GOApollock),
      // started from the spec's obs_sd. This term pins the latent to the observed env
      // series and, with the beta scaling in beta_linkage_eff, identifies the effect
      // size.
      if (linkage_re_obs(grp) >= 0) {
        Type osd = exp(log_obs_sd_linkage(linkage_re_obs(grp)));
        for (int t = 0; t < len; ++t) {
          if (obs_mask(t) == 0) continue;   // year absent from env_data: latent exists, but no observation
          if (osa_mode == 0) {
            jnll_comp(JNLL_LINKAGE_RE, 0)            -= dnorm(obs(t), re(t), osd, true);
            unweighted_jnll_comp(JNLL_LINKAGE_RE, 0) -= dnorm(obs(t), re(t), osd, true);
          } else {
            // OSA: read the covariate observation from obsvec, keep-gated, so
            // oneStepPredict() residualizes it against the latent AR1 state re(t)
            // (WHAM's Ecov OSA). build_osa_data() lays each observed slot in.
            int pos = linkage_re_obsvec_idx(obs_slot(t));
            if (pos >= 0) {
              jnll_comp(JNLL_LINKAGE_RE, 0)            -= keep(pos) * dnorm(obsvec(pos), re(t), osd, true);
              unweighted_jnll_comp(JNLL_LINKAGE_RE, 0) -= keep(pos) * dnorm(obsvec(pos), re(t), osd, true);
              // method = "cdf": the covariate measurement is normal about the latent
              // state, so its CDF is pnorm of the standardized residual. Only jnll_comp,
              // as at every other CDF site: unweighted_jnll_comp is reported as an
              // unweighted LIKELIHOOD (Francis and McAllister-Ianelli read it), and a CDF
              // term is not one. Rows 19 and 20 keep their own values, the wholesale copy
              // in section 13.2 running `ind < 19` and executing before this block.
              if (osa_mode == 2) {
                jnll_comp(JNLL_LINKAGE_RE, 0) -= osa_norm_cdf_terms(
                  obsvec(pos), re(t), osd, keep.cdf_lower(pos), keep.cdf_upper(pos));
              }
            }
          }
        }
      }
    }
  }

  // -- 14.7. Likelihood components
  REPORT( jnll_comp );
  REPORT( unweighted_jnll_comp );


  // -- 14.8. Ration components
  REPORT( consumption_at_age );
  REPORT( fT );

  // -- 14.9. Predation components
  REPORT( M2_prop );
  REPORT( diet_obs );
  REPORT( suitability );
  REPORT( M1_at_age );
  REPORT( M2_at_age );
  REPORT( B_eaten );
  REPORT( B_eaten_as_prey );
  // ADREPORT( B_eaten_as_prey );    // no SEs: see ?plot_b_eaten_prop
  REPORT( vulnerability );
  REPORT( vulnerability_other );
  REPORT( gam_a );
  REPORT( gam_b );
  REPORT( diet_hat );
  REPORT( suit_other );
  REPORT( stom_div_bio );
  REPORT( avail_food );


  // -- 14.10. Initial-state cumulative mortality (feeds the initial age structure)
  REPORT(mort_sum);

  /** ------------------------------------------------------------------------ //
   // 15. END MODEL                                                             //
   * ------------------------------------------------------------------------- */

  Type jnll = 0;

  // Modes 0-3 return the real objective. Mode 3 builds the TMB object without
  // running the optimizer, so obj$fn()/obj$gr() are usable for diagnostics.
  if(estimateMode < 4) {
    jnll = jnll_comp.sum();
  }

  // Mode 4: build_map() maps out every hindcast parameter, leaving `dummy` as
  // the only free parameter. Placeholder objective so nlminb has something
  // valid to minimize -- a plumbing smoke test, not a likelihood.
  if(estimateMode > 3) {
    jnll = dummy * dummy;
  }


  REPORT( jnll );
  return jnll;
}

