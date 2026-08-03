# dev/bench-perf.R
# Read-only runtime + scaling harness for Rceattle. Measures where wall-clock
# goes across three cost centers (single fit, MSE iteration, OSA residuals) at
# several model sizes, and fits scaling exponents so super-linear phases are
# visible. Touches NO package code -- it only wraps public/internal entry points
# with proc.time()/Rprof() and, for the objective, microbenchmarks obj$fn/obj$gr
# directly (estimateMode = 3 keeps the AD object alive and returns the real
# objective). Emphasis is the general-fit / C++ per-objective-eval axis.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript dev/bench-perf.R                       # -> dev/bench-perf-baseline.rds
#   NOT_CRAN=true Rscript dev/bench-perf.R dev/bench-perf-x.rds  # custom out path
#   NOT_CRAN=true Rscript dev/bench-perf.R dev/bench-perf-x.rds compare dev/bench-perf-baseline.rds
#
# Model recipes mirror .claude/commands/golden-check.md EXACTLY (multispecies
# warm-started from the single-species MLEs -- a fresh-init MS fit lands ~37
# units off the pinned optimum and is not comparable).

args         <- commandArgs(trailingOnly = TRUE)
out_path     <- if (length(args) >= 1) args[1] else "dev/bench-perf-baseline.rds"
do_compare   <- length(args) >= 2 && args[2] == "compare"
compare_path <- if (length(args) >= 3) args[3] else NULL

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Keep the harness single-threaded so timings are attributable and reproducible.
Sys.setenv(OMP_NUM_THREADS = "1")
# R-only harness; the TMB DLL is loaded as-is (compile it once beforehand with a
# plain load_all if the .cpp changed). compile = FALSE avoids a recompile race.
suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = FALSE))
set.seed(666)

# ---------------------------------------------------------------------------
# Small timing helpers
# ---------------------------------------------------------------------------

# Elapsed wall-clock (seconds) of an expression, plus its value.
timed <- function(expr) {
  t <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(sec = proc.time()[["elapsed"]] - t, val = val)
}

# Run expr under Rprof and return summed SELF seconds for each named bucket of
# function-name regexes, plus total profiled self time. Buckets are matched
# against summaryRprof()$by.self rownames; a frame is counted in the FIRST
# bucket it matches (order matters), else "other".
prof_buckets <- function(expr, buckets, interval = 0.01) {
  pf <- tempfile(fileext = ".Rprof")
  Rprof(pf, interval = interval, line.profiling = FALSE)
  val <- tryCatch(force(expr), finally = Rprof(NULL))
  s <- tryCatch(summaryRprof(pf), error = function(e) NULL)
  unlink(pf)
  if (is.null(s) || nrow(s$by.self) == 0)
    return(list(val = val, buckets = setNames(rep(NA_real_, length(buckets)), names(buckets)),
                total_self = NA_real_))
  self <- s$by.self$self.time
  names(self) <- rownames(s$by.self)
  assigned <- rep("other", length(self))
  for (bn in names(buckets)) {
    pat <- paste(buckets[[bn]], collapse = "|")
    hit <- assigned == "other" & grepl(pat, names(self))
    assigned[hit] <- bn
  }
  out <- sapply(c(names(buckets), "other"),
                function(bn) sum(self[assigned == bn]))
  list(val = val, buckets = out, total_self = sum(self))
}

# Median seconds per call of f() over `reps` (after `warm` warmups).
microtime <- function(f, reps = 25L, warm = 3L) {
  for (i in seq_len(warm)) f()
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t <- proc.time()[["elapsed"]]; f(); ts[i] <- proc.time()[["elapsed"]] - t
  }
  c(median = median(ts), mean = mean(ts), min = min(ts), n = reps)
}

# ---------------------------------------------------------------------------
# Size metrics for a fitted model
# ---------------------------------------------------------------------------
size_metrics <- function(fit, niter) {
  dl <- fit$data_list
  n_par <- tryCatch(length(unlist(fit$estimated_params)), error = function(e) NA_integer_)
  n_free <- tryCatch(length(fit$obj$par), error = function(e) NA_integer_)
  n_flt <- tryCatch(nrow(dl$fleet_control), error = function(e) NA_integer_)
  n_idx <- tryCatch(nrow(dl$index_data), error = function(e) NA_integer_)
  n_cat <- tryCatch(nrow(dl$catch_data), error = function(e) NA_integer_)
  n_comp <- tryCatch(nrow(dl$comp_data), error = function(e) NA_integer_)
  nyrs <- tryCatch(dl$endyr - dl$styr + 1, error = function(e) NA_integer_)
  c(n_params = n_par, n_free = n_free, n_fleets = n_flt, n_index_obs = n_idx,
    n_catch_obs = n_cat, n_comp_obs = n_comp, n_yrs = nyrs, nspp = dl$nspp,
    niter = niter, total_obs = sum(c(n_idx, n_cat, n_comp), na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# Model ladder (recipes = golden-check)
# ---------------------------------------------------------------------------
data(BS2017SS); data(BS2017MS); data(GOA2018SS)

results <- list(single = list(), objeval = list(), mse = list(), osa = list(), sizes = list())

# ===========================================================================
# INSTRUMENT 1  --  single fit: wall-clock + Rprof buckets (getsd = TRUE)
# ===========================================================================
fit_buckets <- list(
  data_prep = c("rearrange_data", "clean_data", "data_check", "switch_check",
                "build_params", "build_map", "build_bounds", "pool_linkages",
                "_upgrade_", "column_schema"),
  build     = c("MakeADFun"),
  objective = c("^fn$", "^gr$", "^he$", "MakeADFun.*env", "\\.Call", "ADREPORT",
                "EvalADFunObject", "getValues", "TransformADFunObject"),
  optimize  = c("nlminb", "fit_tmb", "newton", "optimHess"),
  sdreport  = c("sdreport", "getJointPrecision", "solve", "Cholesky"),
  post      = c("rename_output", "calc_mcall", "report", "convergence"))

# Fit each ladder point at estimateMode = 0 (hindcast + projection, getsd = TRUE)
# under Rprof, and record wall-clock, bucketed self-time, and optimizer counts.
run_single <- function(label, fit_fun, niter) {
  cat(sprintf("[single] %s ...\n", label))
  wall <- proc.time()[["elapsed"]]
  pb <- prof_buckets(fit_fun(), fit_buckets)
  wall <- proc.time()[["elapsed"]] - wall
  fit <- pb$val
  opt <- fit$opt
  list(
    label = label,
    wall_sec = wall,
    buckets = pb$buckets,
    total_self = pb$total_self,
    nlminb_iter = tryCatch(opt$iterations, error = function(e) NA_integer_),
    nlminb_eval = tryCatch(paste(opt$evaluations, collapse = "/"), error = function(e) NA),
    objective = tryCatch(as.numeric(opt$objective), error = function(e) NA_real_),
    max_grad = tryCatch(max(abs(fit$quantities$mgc), na.rm = TRUE), error = function(e) NA_real_),
    size = size_metrics(fit, niter))
}

single_specs <- list(
  BS_SS = list(niter = 1, fun = function() Rceattle::fit_mod(
    data_list = BS2017SS, file = NULL, inits = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = TRUE, verbose = 0))),
  GOA_SS = list(niter = 1, fun = function() Rceattle::fit_mod(
    data_list = GOA2018SS, file = NULL, inits = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = TRUE, verbose = 0)))
)

# MS points need SS MLEs as warm-start inits -> fit SS first (getsd = FALSE, fast)
cat("[single] warm-start SS fits for MS inits ...\n")
ss_bs  <- suppressMessages(Rceattle::fit_mod(
  data_list = BS2017SS, file = NULL, inits = NULL, estimateMode = 0,
  random_rec = FALSE, msmMode = 0, fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))
ss_goa <- suppressMessages(Rceattle::fit_mod(
  data_list = GOA2018SS, file = NULL, inits = NULL, estimateMode = 0,
  random_rec = FALSE, msmMode = 0, fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))

single_specs$BS_MS <- list(niter = 5, fun = function() suppressMessages(Rceattle::fit_mod(
  data_list = BS2017MS, inits = ss_bs$estimated_params, file = NULL, estimateMode = 0,
  niter = 5, random_rec = FALSE, msmMode = 1, suitMode = 0,
  fit_control = fit_control(verbose = 0))))
single_specs$GOA_MS <- list(niter = 3, fun = function() suppressMessages(Rceattle::fit_mod(
  data_list = GOA2018SS, inits = ss_goa$estimated_params, file = NULL, estimateMode = 0,
  niter = 3, random_rec = FALSE, msmMode = 1, suitMode = 0,
  fit_control = fit_control(phase = TRUE, verbose = 0))))

for (lab in names(single_specs)) {
  sp <- single_specs[[lab]]
  results$single[[lab]] <- tryCatch(run_single(lab, sp$fun, sp$niter),
    error = function(e) list(label = lab, error = conditionMessage(e)))
}

# ===========================================================================
# INSTRUMENT 2 (PRIORITY)  --  objective-eval microbench via estimateMode = 3
# Isolates the C++ per-evaluation cost (obj$fn / obj$gr) from the optimizer, and
# its scaling with niter (Section-6 population dynamics runs niter x per eval)
# and with fleets/obs (Section-13 likelihood). Build-only obj is not freed.
# ===========================================================================
objeval_specs <- list(
  BS_SS  = list(niter = 1, data = quote(BS2017SS),  inits = NULL,                    msm = 0),
  BS_MS  = list(niter = 5, data = quote(BS2017MS),  inits = quote(ss_bs$estimated_params),  msm = 1),
  GOA_SS = list(niter = 1, data = quote(GOA2018SS), inits = NULL,                    msm = 0),
  GOA_MS = list(niter = 3, data = quote(GOA2018SS), inits = quote(ss_goa$estimated_params), msm = 1)
)
# Extra niter sweep on BS-MS (same data) to isolate the Section-6 exponent in niter.
niter_sweep <- c(1, 2, 3, 5, 8)

build_mode3 <- function(dl, inits, msm, niter) {
  suppressMessages(Rceattle::fit_mod(
    data_list = dl, inits = inits, file = NULL, estimateMode = 3,
    niter = niter, random_rec = FALSE, msmMode = msm, suitMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0)))
}
objeval_one <- function(fit) {
  obj <- fit$obj
  p <- obj$par
  # TMB caches fn/gr on the last par vector; calling with an identical par
  # returns the cached value (~0 ms) and hides the true eval cost. Alternate
  # between two nearby par points so every call forces a fresh sweep.
  eps <- 1e-6 * (abs(p) + 1)
  i <- 0L
  perturbed <- function() { i <<- i + 1L; p + (as.numeric(i %% 2L) - 0.5) * eps }
  fn <- microtime(function() obj$fn(perturbed()))
  gr <- microtime(function() obj$gr(perturbed()))
  list(fn = fn, gr = gr, n_free = length(p))
}
for (lab in names(objeval_specs)) {
  sp <- objeval_specs[[lab]]
  cat(sprintf("[objeval] %s ...\n", lab))
  results$objeval[[lab]] <- tryCatch({
    dl <- eval(sp$data); inits <- if (is.null(sp$inits)) NULL else eval(sp$inits)
    fit <- build_mode3(dl, inits, sp$msm, sp$niter)
    c(objeval_one(fit), list(size = size_metrics(fit, sp$niter)))
  }, error = function(e) list(label = lab, error = conditionMessage(e)))
}
# niter sweep (BS multispecies, shared data)
cat("[objeval] niter sweep on BS-MS ...\n")
results$objeval$niter_sweep <- tryCatch({
  lapply(niter_sweep, function(ni) {
    fit <- build_mode3(BS2017MS, ss_bs$estimated_params, 1, ni)
    ev <- objeval_one(fit)
    list(niter = ni, fn_med = ev$fn[["median"]], gr_med = ev$gr[["median"]], n_free = ev$n_free)
  })
}, error = function(e) list(error = conditionMessage(e)))

# ===========================================================================
# INSTRUMENT 3  --  MSE iteration: run_mse serial, Rprof buckets, per-iteration
# ===========================================================================
mse_buckets <- list(
  refit     = c("refit_like", "fit_mod", "MakeADFun", "\\.fit_tmb", "nlminb", "newton"),
  data_prep = c("rearrange_data", "clean_data", "data_check", "switch_check",
                "build_params", "build_map", "_upgrade_"),
  simulate  = c("sim_mod", "sample_rec", "simulate"),
  reshape   = c("arrange", "bind_rows", "rbind", "filter", "slice", "mutate", "\\[\\.data"),
  objective = c("^fn$", "^gr$", "\\.Call", "EvalADFunObject"))

run_mse_bench <- function(label, om, em, n_assess = 3) {
  cat(sprintf("[mse] %s ...\n", label))
  wall <- proc.time()[["elapsed"]]
  pb <- prof_buckets(suppressMessages(Rceattle::run_mse(
    om = om, em = em, nsim = 1, assessment_period = 1, sampling_period = 1,
    simulate_data = TRUE, sample_rec = TRUE, regenerate_past = FALSE,
    seed = 42, cores = 1)), mse_buckets)
  wall <- proc.time()[["elapsed"]] - wall
  list(label = label, wall_sec = wall, n_assess = n_assess,
       per_assess_sec = wall / n_assess, buckets = pb$buckets, total_self = pb$total_self)
}

# SS MSE: Tier-3 EM (F40/F35). OM/EM fit with Proj_F_proportion activated and
# projyr set so the projection has years to advance through (mirrors
# verify-refit-like.R). n_assess = projyr - endyr assessment years.
results$mse$BS_SS <- tryCatch({
  dl <- BS2017SS
  dl$projyr <- dl$endyr + 3
  dl$fleet_control$Proj_F_proportion <- rep(1, nrow(dl$fleet_control))
  om <- suppressMessages(Rceattle::fit_mod(
    data_list = dl, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))
  em <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dl, inits = om$estimated_params, estimateMode = 2,
    HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2, Alpha = 0.05),
    msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0))))
  run_mse_bench("BS_SS", om, em, 3)
}, error = function(e) list(label = "BS_SS", error = conditionMessage(e)))

# MS MSE: ConstantF EM (multispecies supports HCR 0/1/2/3/6), warm-started.
results$mse$BS_MS <- tryCatch({
  dl <- BS2017MS
  dl$projyr <- dl$endyr + 3
  dl$fleet_control$Proj_F_proportion <- rep(1, nrow(dl$fleet_control))
  om <- suppressMessages(Rceattle::fit_mod(
    data_list = dl, inits = ss_bs$estimated_params, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 1, suitMode = 0, niter = 5,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
  em <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dl, inits = om$estimated_params, estimateMode = 2,
    HCR = build_hcr(HCR = 2, Ftarget = rep(0.15, om$data_list$nspp)),
    msmMode = 1, suitMode = 0, niter = 5,
    fit_control = fit_control(getsd = FALSE, verbose = 0))))
  run_mse_bench("BS_MS", om, em, 3)
}, error = function(e) list(label = "BS_MS", error = conditionMessage(e)))

# ===========================================================================
# INSTRUMENT 4  --  OSA residuals: overall + per source, on a converged fit
# ===========================================================================
osa_bench <- function(label, fit) {
  cat(sprintf("[osa] %s ...\n", label))
  out <- list(label = label)
  # component timings via internal entry points (single retape reused)
  bd <- tryCatch(timed(Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)),
                 error = function(e) list(sec = NA, val = NULL, err = conditionMessage(e)))
  out$build_osa_data_sec <- bd$sec
  if (!is.null(bd$val)) {
    ob <- tryCatch(timed(Rceattle:::.osa_build_obj(fit, bd$val)),
                   error = function(e) list(sec = NA, err = conditionMessage(e)))
    out$osa_build_obj_sec <- ob$sec
  }
  # per-source end-to-end (each call is a full osa_residuals for that source)
  srcs <- c("index", "catch", "comp", "caal")
  per <- list()
  for (s in srcs) {
    r <- tryCatch(timed(Rceattle::osa_residuals(fit, source = s, seed = 123)),
                  error = function(e) list(sec = NA, val = NULL, err = conditionMessage(e)))
    per[[s]] <- list(sec = r$sec, n_resid = tryCatch(nrow(r$val), error = function(e) NA_integer_),
                     err = if (is.null(r$val)) r$err else NULL)
  }
  out$per_source <- per
  # all-at-once
  ra <- tryCatch(timed(Rceattle::osa_residuals(fit, seed = 123)),
                 error = function(e) list(sec = NA, val = NULL, err = conditionMessage(e)))
  out$all_sec <- ra$sec
  out$all_n_resid <- tryCatch(nrow(ra$val), error = function(e) NA_integer_)
  out
}

# Converged fits for OSA (hindcast only -> obj available at MLE). getsd not needed.
results$osa$BS_SS <- tryCatch({
  f <- suppressMessages(Rceattle::fit_mod(
    data_list = BS2017SS, file = NULL, inits = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0, fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))
  osa_bench("BS_SS", f)
}, error = function(e) list(label = "BS_SS", error = conditionMessage(e)))

results$osa$BS_MS <- tryCatch({
  f <- suppressMessages(Rceattle::fit_mod(
    data_list = BS2017MS, inits = ss_bs$estimated_params, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 1, suitMode = 0, niter = 5,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
  osa_bench("BS_MS", f)
}, error = function(e) list(label = "BS_MS", error = conditionMessage(e)))

# ===========================================================================
# SCALING TABLE  --  fit log(time) ~ log(size) per phase where >= 3 usable points
# ===========================================================================
scaling_exponent <- function(size, time) {
  ok <- is.finite(size) & is.finite(time) & size > 0 & time > 0
  if (sum(ok) < 3) return(c(k = NA_real_, r2 = NA_real_, n = sum(ok)))
  fit <- lm(log(time[ok]) ~ log(size[ok]))
  c(k = unname(coef(fit)[2]), r2 = summary(fit)$r.squared, n = sum(ok))
}

# Objective-eval scaling in niter (Section-6): from the BS-MS niter sweep.
ns <- results$objeval$niter_sweep
if (is.list(ns) && is.null(ns$error)) {
  niter_v <- sapply(ns, `[[`, "niter")
  fn_v    <- sapply(ns, `[[`, "fn_med")
  gr_v    <- sapply(ns, `[[`, "gr_med")
  results$scaling$objeval_fn_vs_niter <- scaling_exponent(niter_v, fn_v)
  results$scaling$objeval_gr_vs_niter <- scaling_exponent(niter_v, gr_v)
}
# Single-fit wall-clock scaling vs n_params across the ladder.
sg <- results$single
lab_ok <- names(sg)[sapply(sg, function(x) is.null(x$error))]
if (length(lab_ok) >= 3) {
  npar <- sapply(lab_ok, function(l) sg[[l]]$size[["n_free"]])
  wall <- sapply(lab_ok, function(l) sg[[l]]$wall_sec)
  results$scaling$single_wall_vs_nparams <- scaling_exponent(npar, wall)
  for (bk in c("data_prep","build","objective","optimize","sdreport")) {
    bt <- sapply(lab_ok, function(l) sg[[l]]$buckets[[bk]])
    results$scaling[[paste0("single_", bk, "_vs_nparams")]] <- scaling_exponent(npar, bt)
  }
}

# ---------------------------------------------------------------------------
# Persist + print summary
# ---------------------------------------------------------------------------
results$meta <- list(when = as.character(Sys.time()),
                     branch = tryCatch(system("git rev-parse --abbrev-ref HEAD", intern = TRUE),
                                       error = function(e) NA),
                     sha = tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                                    error = function(e) NA))
saveRDS(results, out_path)
cat("\nwrote ->", out_path, "\n")

cat("\n=== SINGLE FIT (wall sec | data_prep/build/objective/optimize/sdreport) ===\n")
for (l in names(results$single)) {
  x <- results$single[[l]]
  if (!is.null(x$error)) { cat(sprintf("  %-7s ERROR %s\n", l, x$error)); next }
  b <- x$buckets
  cat(sprintf("  %-7s wall=%6.1fs  dp=%4.1f bld=%4.1f obj=%5.1f opt=%5.1f sd=%5.1f  npar=%d niter=%g iters=%s\n",
              l, x$wall_sec, b[["data_prep"]], b[["build"]], b[["objective"]],
              b[["optimize"]], b[["sdreport"]], x$size[["n_free"]], x$size[["niter"]], x$nlminb_iter))
}
cat("\n=== OBJECTIVE EVAL (median ms per fn / gr) ===\n")
for (l in setdiff(names(results$objeval), "niter_sweep")) {
  x <- results$objeval[[l]]
  if (!is.null(x$error)) { cat(sprintf("  %-7s ERROR %s\n", l, x$error)); next }
  cat(sprintf("  %-7s fn=%7.2fms gr=%7.2fms  nfree=%d niter=%g\n",
              l, 1000*x$fn[["median"]], 1000*x$gr[["median"]], x$n_free, x$size[["niter"]]))
}
if (is.list(ns) && is.null(ns$error)) {
  cat("  niter sweep (BS-MS): ")
  cat(paste(sapply(ns, function(z) sprintf("n%d fn=%.1fms", z$niter, 1000*z$fn_med)), collapse = "  "), "\n")
}
cat("\n=== MSE (wall | per-assessment-year | refit/prep/sim/reshape) ===\n")
for (l in names(results$mse)) {
  x <- results$mse[[l]]
  if (!is.null(x$error)) { cat(sprintf("  %-7s ERROR %s\n", l, x$error)); next }
  b <- x$buckets
  cat(sprintf("  %-7s wall=%6.1fs  /yr=%5.1fs  refit=%5.1f prep=%4.1f sim=%4.1f reshape=%4.1f\n",
              l, x$wall_sec, x$per_assess_sec, b[["refit"]], b[["data_prep"]], b[["simulate"]], b[["reshape"]]))
}
cat("\n=== OSA (sec: build_data/build_obj | per-source | all) ===\n")
for (l in names(results$osa)) {
  x <- results$osa[[l]]
  if (!is.null(x$error)) { cat(sprintf("  %-7s ERROR %s\n", l, x$error)); next }
  ps <- x$per_source
  cat(sprintf("  %-7s bd=%4.1f bo=%4.1f | idx=%4.1f cat=%4.1f comp=%5.1f caal=%4.1f | all=%5.1f (n=%s)\n",
              l, x$build_osa_data_sec, x$osa_build_obj_sec %||% NA,
              ps$index$sec, ps$catch$sec, ps$comp$sec, ps$caal$sec, x$all_sec, x$all_n_resid))
}
cat("\n=== SCALING EXPONENTS (k = d log time / d log size) ===\n")
for (nm in names(results$scaling)) {
  s <- results$scaling[[nm]]
  cat(sprintf("  %-34s k=%6.2f  r2=%4.2f  n=%d\n", nm, s[["k"]], s[["r2"]], s[["n"]]))
}

if (do_compare && !is.null(compare_path)) {
  before <- readRDS(compare_path)
  cat("\n=== COMPARE vs", compare_path, "(wall-clock deltas) ===\n")
  cmp <- function(sec, a, b) {
    for (l in intersect(names(a), names(b))) {
      av <- a[[l]][[sec]]; bv <- b[[l]][[sec]]
      if (is.numeric(av) && is.numeric(bv))
        cat(sprintf("  %-8s %-8s %6.1fs -> %6.1fs  (%+5.1f%%)\n",
                    sec, l, bv, av, 100*(av-bv)/bv))
    }
  }
  cmp("wall_sec", results$single, before$single)
  cmp("wall_sec", results$mse, before$mse)
  cmp("all_sec", results$osa, before$osa)
}
