// =====================================================================
// linkage.hpp
//
// Turns each row of the Rceattle linkage table -- encoded by R into parallel
// IVECTORs/VECTORs and the design matrix `linkage_X` -- into a tensor of
// year-varying offsets that the recruitment, mortality, growth, catchability
// and selectivity modules add to their parameters.
//
// Process and parameter codes are in lockstep with R/0-linkage_encode.R
// (LINKAGE_PROCESS_CODES / LINKAGE_PARAM_CODES):
//   PROC: recruitment=0, M=1, growth=2, q=3, sel=4, comp=5
//   PARAM is a per-process namespace: growth K=0, L1=1, Linf=2, m=3
//     (natural-scale names) plus the SD endpoints 4/5; M1=0 for M.
//
// A stratum id of 0 in `linkage_species`, `linkage_sex`, `linkage_age_bin` or
// `linkage_fleet` applies the row to every level of that stratum; otherwise it
// is a 1-based index into that stratum's levels.
// =====================================================================
#ifndef RCEATTLE_LINKAGE_HPP
#define RCEATTLE_LINKAGE_HPP

// Process codes (must match R/0-linkage_encode.R::LINKAGE_PROCESS_CODES)
#define RCEATTLE_PROC_RECRUIT 0
#define RCEATTLE_PROC_M       1
#define RCEATTLE_PROC_GROWTH  2
#define RCEATTLE_PROC_Q       3
#define RCEATTLE_PROC_SEL     4
#define RCEATTLE_PROC_COMP    5   // prior-only: DM composition-weighting overdispersion

// Growth parameters in the offset tensor: log_K, log_L1, log_Linf, log_m.
#define RCEATTLE_N_GROWTH_PARAMS 4

// Recruitment parameters in the offset tensor: R0, alpha and beta, held on the log scale
// in `rec_pars`, plus R_init, the initial recruitment level, read at year 0 only.
#define RCEATTLE_N_REC_PARAMS 4
#define RCEATTLE_REC_R0    0
#define RCEATTLE_REC_ALPHA 1
#define RCEATTLE_REC_BETA  2
#define RCEATTLE_REC_R_INIT 3

/**
 * @brief Expand a stratum sentinel id into the half-open iteration range [lo, hi).
 *
 * Clamped to the levels this species actually has, so a sex-2 row on a single-sex stock
 * resolves to an empty range.
 *
 * @param id 1-based stratum id, or 0 for "all levels".
 * @param n_levels Number of levels this stratum has for the current species.
 * @param lo,hi [out] The range bounds: `lo` inclusive, `hi` exclusive.
 */
inline void rceattle_stratum_range(int id, int n_levels, int& lo, int& hi) {
  lo = (id == 0) ? 0 : (id - 1);
  hi = (id == 0) ? n_levels : id;
  if (hi > n_levels) hi = n_levels;
}


// ---------------------------------------------------------------------
// Every accumulator below has one shape: for each linkage row of this process on this
// link scale, expand the stratum sentinels and add beta * X(yr, col) into the offset
// tensor for every (stratum, year) the row applies to. Three arguments are shared and
// documented only here: `linkage_X`, the environmental covariate matrix whose rows are
// years; `beta`, the per-row effect sizes (0-length = no-op); and `link_code`, which
// picks the rows a call consumes -- 1 = log (added inside the exp), 0 = identity (added
// to the natural-scale value after it), 3 = exponential (MULTIPLIES the log, raising the
// parameter to exp(beta * x); catchability only, where index_q = exp((index_log_q +
// log_offset) * exp(log_mult) + dev) + nat_offset). Years past linkage_X.rows() keep a
// zero offset, so env_data need not span the projection. See
// vignette("environmental-linkages-and-priors").
// ---------------------------------------------------------------------


/**
 * @brief Accumulate growth offsets into `growth_offset[sp, sex, yr, param]`.
 *
 * Not age-stratified, so the age_bin sentinel is ignored. The growth SD codes (sd_L1 = 4,
 * sd_Linf = 5) are skipped by the param-range guard below: they live on `growth_log_sd`
 * and consume intercept rows only.
 *
 * @param growth_offset [in,out] Offset tensor [nspp, max_sex, nyrs, n_growth_params].
 * @param link_code,linkage_X,beta Shared arguments; see the note above this group.
 */
template<class Type>
void rceattle_apply_growth_linkages(
    array<Type>&          growth_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   nspp,
    const vector<int>&    nsex,
    int                   nyrs)
{
  (void)linkage_age_bin;               // growth is not age-stratified
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_GROWTH) continue;
    if (linkage_link(i)    != link_code)            continue;

    int param = linkage_param(i);
    if (param < 0 || param >= RCEATTLE_N_GROWTH_PARAMS) continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int sp_lo, sp_hi;
    rceattle_stratum_range(linkage_species(i), nspp, sp_lo, sp_hi);

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      int sx_lo, sx_hi;
      rceattle_stratum_range(linkage_sex(i), nsex(sp), sx_lo, sx_hi);

      for (int sx = sx_lo; sx < sx_hi; ++sx) {
        for (int yr = 0; yr < yr_hi; ++yr) {
          growth_offset(sp, sx, yr, param) += b * linkage_X(yr, xc);
        }
      }
    }
  }
}


/**
 * @brief Accumulate natural-mortality offsets into `M_offset[sp, sex, age, yr]`.
 *
 * Age-stratified, unlike growth; `log_M1` (M1 = 0) is the only param code.
 *
 * @param M_offset [in,out] Offset tensor [nspp, max_sex, max_age, nyrs].
 * @param link_code,linkage_X,beta Shared arguments; see the note above this group.
 */
template<class Type>
void rceattle_apply_M_linkages(
    array<Type>&          M_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   nspp,
    const vector<int>&    nsex,
    const vector<int>&    nages,
    int                   nyrs)
{
  (void)linkage_param;                 // only log_M1 is exposed so far
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_M) continue;
    if (linkage_link(i)    != link_code)       continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int sp_lo, sp_hi;
    rceattle_stratum_range(linkage_species(i), nspp, sp_lo, sp_hi);

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      int sx_lo, sx_hi, ab_lo, ab_hi;
      rceattle_stratum_range(linkage_sex(i),     nsex(sp),  sx_lo, sx_hi);
      rceattle_stratum_range(linkage_age_bin(i), nages(sp), ab_lo, ab_hi);

      for (int sx = sx_lo; sx < sx_hi; ++sx) {
        for (int ab = ab_lo; ab < ab_hi; ++ab) {
          for (int yr = 0; yr < yr_hi; ++yr) {
            M_offset(sp, sx, ab, yr) += b * linkage_X(yr, xc);
          }
        }
      }
    }
  }
}


/**
 * @brief Accumulate recruitment offsets into `rec_offset[sp, param, yr]`.
 *
 * Recruitment happens at age 0 with no sex stratification on the parameter itself
 * (`sex_ratio` splits recruits downstream), so the sex and age_bin sentinels are ignored.
 * Param codes R0 = 0, alpha = 1, beta = 2, R_init = 3.
 *
 * @param rec_offset [in,out] Offset tensor [nspp, n_rec_params, nyrs].
 * @param link_code,linkage_X,beta Shared arguments; see the note above this group.
 */
template<class Type>
void rceattle_apply_recruitment_linkages(
    array<Type>&          rec_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   nspp,
    int                   nyrs)
{
  (void)linkage_sex;
  (void)linkage_age_bin;
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_RECRUIT) continue;
    if (linkage_link(i)    != link_code)             continue;

    int param = linkage_param(i);
    if (param < 0 || param >= RCEATTLE_N_REC_PARAMS) continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int sp_lo, sp_hi;
    rceattle_stratum_range(linkage_species(i), nspp, sp_lo, sp_hi);

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      for (int yr = 0; yr < yr_hi; ++yr) {
        rec_offset(sp, param, yr) += b * linkage_X(yr, xc);
      }
    }
  }
}


/**
 * @brief Accumulate catchability offsets into `q_offset[fleet, yr]`.
 *
 * Indexed by fleet rather than species/sex/age -- a fleet already implies its species
 * through `fleet_control` -- so those sentinels are ignored, and `q` (q = 0) is the only
 * param code.
 *
 * @param q_offset [in,out] Offset matrix [n_flt, nyrs].
 * @param linkage_fleet Fleet stratum id per row (0 = all fleets).
 * @param link_code,linkage_X,beta Shared arguments; see the note above this group.
 */
template<class Type>
void rceattle_apply_q_linkages(
    matrix<Type>&         q_offset,
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_fleet,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   n_flt,
    int                   nyrs)
{
  (void)linkage_param;                 // only q itself is exposed
  (void)linkage_species;
  (void)linkage_sex;
  (void)linkage_age_bin;
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_Q) continue;
    if (linkage_link(i)    != link_code)       continue;

    int xc = linkage_X_col(i);
    Type b = beta(i);

    int fl_lo, fl_hi;
    rceattle_stratum_range(linkage_fleet(i), n_flt, fl_lo, fl_hi);

    for (int flt = fl_lo; flt < fl_hi; ++flt) {
      for (int yr = 0; yr < yr_hi; ++yr) {
        q_offset(flt, yr) += b * linkage_X(yr, xc);
      }
    }
  }
}


/**
 * @brief Accumulate selectivity offsets into the five parameter-family tensors.
 *
 * Routed by the row's param code, for every (fleet, sex, year) the row applies to. The
 * tensors are dimensioned per fleet and sex; fleets sharing a `Selectivity_index` are
 * collapsed onto one parameter block by the TMB map.
 *
 * Param codes: `0`/`1` -> log_sel_slp[asc/desc] (`slp_offset`);
 *              `2`/`3` -> sel_inf[asc/desc]      (`inf_offset`);
 *              `4`     -> sel_coff, all bins     (`coff_offset`);
 *              `5`     -> log_sel_apical         (`apical_offset`);
 *              `6`-`11` -> sel_dn6[0..5], the six DoubleNormalSS3 (SS3
 *                         pattern 24) parameters (`dn6_offset`).
 *
 * @param slp_offset [in,out] Slope offsets [2, n_flt, max_sex, nyrs].
 * @param inf_offset [in,out] Inflection offsets [2, n_flt, max_sex, nyrs].
 * @param coff_offset [in,out] Nonparametric offsets [n_flt, max_sex, n_sel_bins, nyrs].
 * @param apical_offset [in,out] Per-sex apical-height offsets [n_flt, max_sex, nyrs].
 * @param dn6_offset [in,out] DoubleNormalSS3 parameter offsets [6, n_flt, max_sex, nyrs].
 * @param link_code,linkage_X,beta Shared arguments; see the note above this group.
 */
template<class Type>
void rceattle_apply_sel_linkages(
    array<Type>&          slp_offset,   // [2, n_flt, max_sex, nyrs]
    array<Type>&          inf_offset,   // [2, n_flt, max_sex, nyrs]
    array<Type>&          coff_offset,  // [n_flt, max_sex, n_sel_bins, nyrs]
    array<Type>&          apical_offset,// [n_flt, max_sex, nyrs]
    array<Type>&          dn6_offset,   // [6, n_flt, max_sex, nyrs]
    int                   link_code,
    const vector<int>&    linkage_process,
    const vector<int>&    linkage_param,
    const vector<int>&    linkage_species,
    const vector<int>&    linkage_sex,
    const vector<int>&    linkage_age_bin,
    const vector<int>&    linkage_fleet,
    const vector<int>&    linkage_X_col,
    const vector<int>&    linkage_link,
    const matrix<Type>&   linkage_X,
    const vector<Type>&   beta,
    int                   n_flt,
    int                   max_sex,
    int                   n_sel_bins,
    int                   nyrs)
{
  (void)linkage_species;
  (void)linkage_age_bin;
  int n = beta.size();
  if (n == 0) return;

  int yr_hi = std::min(nyrs, (int)linkage_X.rows());

  for (int i = 0; i < n; ++i) {
    if (linkage_process(i) != RCEATTLE_PROC_SEL) continue;
    if (linkage_link(i)    != link_code)         continue;

    int param = linkage_param(i);
    int xc    = linkage_X_col(i);
    Type b    = beta(i);

    int fl_lo, fl_hi, sx_lo, sx_hi;
    rceattle_stratum_range(linkage_fleet(i), n_flt,   fl_lo, fl_hi);
    rceattle_stratum_range(linkage_sex(i),   max_sex, sx_lo, sx_hi);

    for (int flt = fl_lo; flt < fl_hi; ++flt) {
      for (int sx = sx_lo; sx < sx_hi; ++sx) {
        for (int yr = 0; yr < yr_hi; ++yr) {
          Type v = b * linkage_X(yr, xc);
          if (param == 0 || param == 1) {
            slp_offset(param, flt, sx, yr) += v;
          } else if (param == 2 || param == 3) {
            inf_offset(param - 2, flt, sx, yr) += v;
          } else if (param == 4) {
            for (int bin = 0; bin < n_sel_bins; ++bin) {
              coff_offset(flt, sx, bin, yr) += v;
            }
          } else if (param == 5) {
            apical_offset(flt, sx, yr) += v;
          } else if (param == 6) {
            dn6_offset(0, flt, sx, yr) += v;   // P1 peak
          } else if (param == 7) {
            dn6_offset(1, flt, sx, yr) += v;   // P2 top width (logit)
          } else if (param == 8) {
            dn6_offset(2, flt, sx, yr) += v;   // P3 ascending width (log)
          } else if (param == 9) {
            dn6_offset(3, flt, sx, yr) += v;   // P4 descending width (log)
          } else if (param == 10) {
            dn6_offset(4, flt, sx, yr) += v;   // P5 initial selectivity (logit)
          } else if (param == 11) {
            dn6_offset(5, flt, sx, yr) += v;   // P6 final selectivity (logit)
          }
        }
      }
    }
  }
}

#endif
