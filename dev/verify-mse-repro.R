# dev/verify-mse-repro.R
# Full-MSE reproduction digest. verify-refit-like.R exercises a single tiny
# seeded refit per entry function; it does NOT capture a multi-year MSE
# TRAJECTORY (OM advanced + EM re-assessed over several assessment years, with
# data simulation and recruitment resampling active). This script does: a
# fixed-seed short MSE on BS2017SS (single-species), BS2017MS (multispecies,
# exercises the pinned predation-suitability window) and GOA2018SS (larger fleet
# count), digested to a numeric fingerprint of every fit node in the returned
# trajectory. It is the bit-identical gate for any future MSE optimization
# (AD-object reuse, invariant-work hoisting, C++ micro-ops) -- capture on the
# clean baseline, require identical afterwards.
#
# Determinism: cores = 1 (no scheduler reorder) and each sim self-seeds from
# seed + sim, so simulate_data / sample_rec stay reproducible.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript dev/verify-mse-repro.R dev/mse-repro-before.rds
#   NOT_CRAN=true Rscript dev/verify-mse-repro.R dev/mse-repro-after.rds compare dev/mse-repro-before.rds

args         <- commandArgs(trailingOnly = TRUE)
out_path     <- if (length(args) >= 1) args[1] else "dev/mse-repro-digest.rds"
do_compare   <- length(args) >= 2 && args[2] == "compare"
compare_path <- if (length(args) >= 3) args[3] else NULL

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# Same walker as verify-refit-like.R: recurse any structure, fingerprint every
# node that has $quantities (both optimized fits and estimateMode = 3 build-only
# objects), stop descending there. Deterministic (list order preserved). Also
# capture the MSE keep-quantities so the digest reflects everything a summary
# would consume, not just ssb/R.
KEEP <- Rceattle:::.mse_keep_quantities
collect <- function(x, acc = list()) {
  if (is.list(x) && !is.null(x[["quantities"]])) {
    q <- x$quantities
    acc[[length(acc) + 1L]] <- list(
      obj  = tryCatch(as.numeric(x$opt$objective), error = function(e) NA_real_),
      ssb  = tryCatch(as.numeric(q$ssb), error = function(e) NA_real_),
      R    = tryCatch(as.numeric(q$R),   error = function(e) NA_real_),
      keep = lapply(KEEP, function(nm) tryCatch(as.numeric(q[[nm]]), error = function(e) NA_real_)))
    return(acc)
  }
  if (is.list(x)) for (el in x) acc <- collect(el, acc)
  acc
}
section <- function(expr) {
  r <- tryCatch(expr, error = function(e) structure(list(.err = conditionMessage(e)),
                                                    class = "harness_err"))
  if (inherits(r, "harness_err")) return(r$.err)
  collect(r)
}

# Activate the projection: Proj_F_proportion must sum to 1 across each species'
# FISHERY fleets (Fleet_type == 1); surveys/other get 0. Equal split is fine for
# a repro fixture. Single-fishery species (all of BS) get 1 -- numerically the
# same as the rep(1, nfleets) form used in verify-refit-like.R, since projection
# F only reads fishery proportions.
activate_proj <- function(dl) {
  fc <- dl$fleet_control
  p <- rep(0, nrow(fc))
  is_fish <- fc$Fleet_type == 1
  for (s in unique(fc$Species[is_fish])) {
    idx <- which(is_fish & fc$Species == s)
    p[idx] <- 1 / length(idx)
  }
  dl$fleet_control$Proj_F_proportion <- p
  dl
}

digest <- list()

# --- BS2017SS single-species: Tier-3 EM (F40/F35) --------------------------
digest$bs_ss <- section({
  data(BS2017SS)
  BS2017SS$projyr <- 2019
  BS2017SS <- activate_proj(BS2017SS)
  om <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))
  em <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = om$estimated_params, estimateMode = 2,
    HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2, Alpha = 0.05),
    msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0)))
  Rceattle::run_mse(
    om = om, em = em, nsim = 2, assessment_period = 1, sampling_period = 1,
    simulate_data = TRUE, sample_rec = TRUE, regenerate_past = FALSE,
    seed = 666, cores = 1)
})

# --- BS2017MS multispecies: ConstantF EM, warm-started (pins suitability) ---
digest$bs_ms <- section({
  data(BS2017SS); data(BS2017MS)
  ss <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))
  BS2017MS$projyr <- 2019
  BS2017MS <- activate_proj(BS2017MS)
  om <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017MS, inits = ss$estimated_params, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 1, suitMode = 0, niter = 5,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
  em <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017MS, inits = om$estimated_params, estimateMode = 2,
    HCR = build_hcr(HCR = 2, Ftarget = rep(0.15, om$data_list$nspp)),
    msmMode = 1, suitMode = 0, niter = 5,
    fit_control = fit_control(getsd = FALSE, verbose = 0)))
  Rceattle::run_mse(
    om = om, em = em, nsim = 2, assessment_period = 1, sampling_period = 1,
    simulate_data = TRUE, sample_rec = TRUE, regenerate_past = FALSE,
    seed = 666, cores = 1)
})

# --- GOA2018SS single-species: larger fleet count (16 fleets) ---------------
digest$goa_ss <- section({
  data(GOA2018SS)
  GOA2018SS$projyr <- 2020
  GOA2018SS <- activate_proj(GOA2018SS)
  om <- suppressWarnings(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))
  em <- suppressWarnings(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = om$estimated_params, estimateMode = 2,
    HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2, Alpha = 0.05),
    msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0)))
  Rceattle::run_mse(
    om = om, em = em, nsim = 2, assessment_period = 1, sampling_period = 1,
    simulate_data = TRUE, sample_rec = TRUE, regenerate_past = FALSE,
    seed = 666, cores = 1)
})

saveRDS(digest, out_path)
cat("wrote digest ->", out_path, "\n")
for (nm in names(digest)) {
  v <- digest[[nm]]
  n <- if (is.character(v)) paste0("ERR(", substr(v, 1, 50), ")") else paste0(length(v), " fit node(s)")
  cat(sprintf("  %-8s %s\n", nm, n))
}

if (do_compare && !is.null(compare_path)) {
  before <- readRDS(compare_path)
  cat("\n=== COMPARE vs", compare_path, "===\n")
  if (identical(before, digest)) {
    cat("PASS: bit-identical across all sections.\n")
  } else {
    cat("FAIL: sections differ:\n")
    for (nm in union(names(before), names(digest))) {
      if (!identical(before[[nm]], digest[[nm]])) cat("  DIFF:", nm, "\n")
    }
  }
}
