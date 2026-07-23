// =====================================================================
// linkage.hpp
//
// Generic accumulator for the Rceattle long-format linkage table. The
// R side encodes the pooled `Rceattle_linkage_table` into parallel
// IVECTORs/VECTORs and a dense design matrix `linkage_X`. This header
// reads those, evaluates one (yearly) coefficient contribution per
// linkage row, and accumulates the result into a per-process tensor of
// year-varying offsets that the growth / mortality / recruitment
// modules can add to their linear-predictor parameters.
//
// Process and parameter codes are kept in lockstep with R/0-linkage_encode.R:
//   PROC: recruitment=0, M=1, growth=2, q=3, sel=4
//   PARAM (growth): K=0, L1=1, Linf=2, m=3 (natural-scale names)
//   PARAM (M): M1=0
//   ...
//
// Stratum sentinels: a 0 in `linkage_species`, `linkage_sex`, or
// `linkage_age_bin` means "applies to every level of that stratum".
// Otherwise the value is a 1-based index into 1..nspp / 1..nsex(sp) /
// 1..nages(sp).
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

// Number of growth parameters tracked in the offset tensor
// (log_K, log_L1, log_Linf, log_m).
#define RCEATTLE_N_GROWTH_PARAMS 4

// Number of recruitment parameters tracked in the offset tensor
// (R0, alpha, beta -- natural-scale names; stored as log_R0 etc. in TMB).
#define RCEATTLE_N_REC_PARAMS 3
#define RCEATTLE_REC_R0    0
#define RCEATTLE_REC_ALPHA 1
#define RCEATTLE_REC_BETA  2

// ---------------------------------------------------------------------
// Sentinel expansion.
//
// A stratum id of 0 means "applies to every level of this stratum";
// otherwise it is a 1-based index. Resolve one id into the half-open
// range [lo, hi) to iterate, clamped to the levels this species
// actually has (e.g. a sex-2 row applied to a single-sex stock).
// ---------------------------------------------------------------------
inline void rceattle_stratum_range(int id, int n_levels, int& lo, int& hi) {
  lo = (id == 0) ? 0 : (id - 1);
  hi = (id == 0) ? n_levels : id;
  if (hi > n_levels) hi = n_levels;
}


// ---------------------------------------------------------------------
// The three accumulators below share one shape:
//
//   for each linkage row of this process, on this link scale:
//     expand the stratum sentinels, then add beta * X(yr, col) into the
//     offset tensor for every (stratum, year) the row applies to.
//
// `link_code` selects which rows a call consumes and therefore which
// tensor it fills: 1 = log (added inside the exp at the consume site),
// 0 = identity (added to the natural-scale value afterwards). The
// consumer combines them as
//   value_yr = exp(log_base + log_offset(yr)) + nat_offset(yr)
// so each process is called twice, once per scale, with the matching
// tensor. Years beyond `linkage_X.rows()` keep a zero offset: env_data
// need not span the projection horizon.
// ---------------------------------------------------------------------


// Growth: offset tensor is [nspp, max_nsex, nyrs, n_growth_params].
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


// Natural mortality: offset tensor is [nspp, max_nsex, max_age, nyrs].
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


// Recruitment: offset tensor is [nspp, n_rec_params, nyrs]. Recruitment
// happens at age 0 with no sex stratification on the parameter itself
// (sex_ratio splits the recruits downstream), so those sentinels are
// accepted and ignored.
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


// Catchability: offset tensor is [n_flt, nyrs]. Catchability is indexed by
// fleet rather than by species/sex/age, so those sentinels are accepted and
// ignored -- a fleet already implies its species through fleet_control.
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


// Selectivity: fills three offset tensors, one per parameter family, routed by
// the row's param code. Selectivity is indexed by fleet and sex (mirroring is
// collapsed by the TMB map, so the parameter arrays are effectively per-fleet).
// Param codes: 0/1 -> log_sel_slp[asc/desc] (slp_offset),
//              2/3 -> sel_inf[asc/desc]      (inf_offset),
//              4   -> sel_coff (all bins)    (coff_offset).
// As with the other processes, called once per link scale; the caller passes
// the log-scale or the natural-scale trio of tensors.
template<class Type>
void rceattle_apply_sel_linkages(
    array<Type>&          slp_offset,   // [2, n_flt, max_sex, nyrs]
    array<Type>&          inf_offset,   // [2, n_flt, max_sex, nyrs]
    array<Type>&          coff_offset,  // [n_flt, max_sex, n_sel_bins, nyrs]
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
          }
        }
      }
    }
  }
}

#endif
