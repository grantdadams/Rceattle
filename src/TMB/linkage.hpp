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
//   PARAM (growth): log_K=0, log_L1=1, log_Linf=2, log_m=3
//   PARAM (M): log_M=0
//   ...
//
// Stratum sentinels: a 0 in `linkage_species`, `linkage_sex`, or
// `linkage_age_bin` means "applies to every level of that stratum".
// Otherwise the value is a 1-based index into 1..nspp / 1..nsex(sp) /
// 1..nages(sp).
//
// Step 3B note: this header currently only exposes the accumulator
// against the *growth* parameters and routes its output into a
// 4D array `growth_linkage_offset(sp, sex, growth_param, year)`.
// Subsequent steps will widen the API to recruitment / M / q.
// =====================================================================
#ifndef RCEATTLE_LINKAGE_HPP
#define RCEATTLE_LINKAGE_HPP

// Process codes (must match R/0-linkage_encode.R::LINKAGE_PROCESS_CODES)
#define RCEATTLE_PROC_RECRUIT 0
#define RCEATTLE_PROC_M       1
#define RCEATTLE_PROC_GROWTH  2
#define RCEATTLE_PROC_Q       3
#define RCEATTLE_PROC_SEL     4

// Number of growth parameters tracked in the offset tensor
// (log_K, log_L1, log_Linf, log_m).
#define RCEATTLE_N_GROWTH_PARAMS 4


// Build the per-(species, sex, growth_param, year) offset tensor.
//
// @tparam Type TMB scalar type.
// @param[out] growth_offset 4D array shaped [nspp, max_nsex, n_params, nyrs].
//   Caller is expected to allocate and zero-initialize it.
// @param[in] linkage_process,linkage_param,linkage_species,linkage_sex,
//   linkage_age_bin,linkage_X_col,linkage_link the encoded table columns.
// @param[in] linkage_X dense design matrix, [nyrs x n_design_cols].
// @param[in] beta the parameter vector aligned with the linkage rows.
// @param[in] nspp,max_nsex,nyrs dimensions of the offset tensor.
// @param[in] nsex per-species sex count for sentinel expansion.
template<class Type>
void rceattle_apply_linkages(
    array<Type>& growth_offset,
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
  // age_bin is reserved for parameters that vary by age (growth dispatch
  // does not yet read it; mortality will). Reference once to silence
  // -Wunused-parameter while we wait for that wiring.
  (void)linkage_age_bin;

  int n = beta.size();
  if (n == 0) return;

  for (int i = 0; i < n; ++i) {
    int proc = linkage_process(i);
    if (proc != RCEATTLE_PROC_GROWTH) continue;

    int param = linkage_param(i);
    if (param < 0 || param >= RCEATTLE_N_GROWTH_PARAMS) continue;

    int sp_in  = linkage_species(i);
    int sx_in  = linkage_sex(i);
    int xc     = linkage_X_col(i);
    int linkfn = linkage_link(i);
    Type b     = beta(i);

    int sp_lo = (sp_in == 0) ? 0 : (sp_in - 1);
    int sp_hi = (sp_in == 0) ? nspp : sp_in;

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      int max_sx_for_sp = nsex(sp);
      int sx_lo = (sx_in == 0) ? 0 : (sx_in - 1);
      int sx_hi = (sx_in == 0) ? max_sx_for_sp : sx_in;
      // Skip rows whose requested sex stratum is not represented for
      // this species (e.g. a sex-2 row applied to a single-sex stock).
      if (sx_hi > max_sx_for_sp) sx_hi = max_sx_for_sp;

      for (int sx = sx_lo; sx < sx_hi; ++sx) {
        // The growth params we accumulate for (log_K, log_L1, log_Linf,
        // log_m) all live on a log scale already, so an "identity" link
        // and a "log" link both contribute additively to the linear
        // predictor. The branch is a placeholder for future processes
        // (e.g. logit for steepness on a 0-1 scale).
        (void)linkfn;
        // Clamp to the rows actually present in the design matrix:
        // env_data may not span the full projection horizon. Years
        // beyond linkage_X.rows() retain a zero offset, i.e. the
        // unperturbed value of the underlying growth parameter.
        int yr_hi = std::min(nyrs, (int)linkage_X.rows());
        for (int yr = 0; yr < yr_hi; ++yr) {
          growth_offset(sp, sx, param, yr) += b * linkage_X(yr, xc);
        }
      }
    }
  }
}


// =====================================================================
// M (natural mortality) accumulator.
//
// Builds the per-(species, sex, age, year) offset tensor that is
// added (additively, on the log scale) to ln_M1 inside the M1_at_age
// compute in ceattle_v01_11.cpp. Iterates the same encoded linkage
// table but only consumes rows whose process == RCEATTLE_PROC_M.
//
// Stratum sentinels: a 0 in `linkage_species`, `linkage_sex`, or
// `linkage_age_bin` expands to "every level of that stratum". A
// linkage row with all three at 0 thus broadcasts a single beta
// across the entire (species, sex, age) cube.
// =====================================================================
template<class Type>
void rceattle_apply_M_linkages(
    array<Type>&          M_offset,        // [nspp, max_sex, max_age, nyrs]
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
  int n = beta.size();
  if (n == 0) return;

  for (int i = 0; i < n; ++i) {
    int proc = linkage_process(i);
    if (proc != RCEATTLE_PROC_M) continue;

    // Currently only one M parameter is exposed (log_M).
    // linkage_param is reserved for future M parameters but is
    // accepted as-is here.
    (void)linkage_param;

    int sp_in  = linkage_species(i);
    int sx_in  = linkage_sex(i);
    int ab_in  = linkage_age_bin(i);
    int xc     = linkage_X_col(i);
    int linkfn = linkage_link(i);
    Type b     = beta(i);
    (void)linkfn;

    int sp_lo = (sp_in == 0) ? 0 : (sp_in - 1);
    int sp_hi = (sp_in == 0) ? nspp : sp_in;

    for (int sp = sp_lo; sp < sp_hi; ++sp) {
      int max_sx_for_sp = nsex(sp);
      int max_age_for_sp = nages(sp);
      int sx_lo = (sx_in == 0) ? 0 : (sx_in - 1);
      int sx_hi = (sx_in == 0) ? max_sx_for_sp : sx_in;
      if (sx_hi > max_sx_for_sp) sx_hi = max_sx_for_sp;
      int ab_lo = (ab_in == 0) ? 0 : (ab_in - 1);
      int ab_hi = (ab_in == 0) ? max_age_for_sp : ab_in;
      if (ab_hi > max_age_for_sp) ab_hi = max_age_for_sp;

      // Clamp years to the design matrix length: env_data may not
      // span the full projection horizon. Years beyond linkage_X.rows()
      // retain a zero offset.
      int yr_hi = std::min(nyrs, (int)linkage_X.rows());

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

#endif  // RCEATTLE_LINKAGE_HPP
