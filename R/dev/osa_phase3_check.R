# Phase 3 OSA validation (diet). Run via: Rscript R/dev/osa_phase3_check.R
suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))
source(file.path("tests", "testthat", "helpers-make-msm-data.R"))
source(file.path("tests", "testthat", "helpers.R"))

ok <- TRUE
check <- function(label, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) ok <<- FALSE
}

# ---- 1. build_osa_data diet block (direct unit test) ----------------------
# A minimal data_list with two stomachs of predator species 1 (suitMode > 0).
empty_i <- function(nc) matrix(integer(0), 0, nc)
empty_n <- function(nc) matrix(numeric(0), 0, nc)
dl <- list(
  endyr = 5, flt_type = c(1L, 2L), nages = c(3L, 3L), nlengths = c(3L, 3L),
  index_ctl = empty_i(3), index_obs = empty_n(2),
  catch_ctl = empty_i(3), catch_obs = empty_n(2),
  comp_ctl = empty_i(5), comp_obs = empty_n(3), comp_n = empty_n(2),
  caal_ctl = empty_i(5), caal_obs = empty_n(3), caal_n = empty_n(1),
  n_stomach_obs = 2L,
  stomach_id = c(0L, 0L, 1L),                       # stomach 0: 2 prey, stomach 1: 1 prey
  diet_ctl = matrix(c(1, 1, 0, 0, 1, 1, 3,
                      1, 2, 0, 0, 1, 1, 3,
                      1, 1, 0, 0, 2, 1, 4),
                    nrow = 3, byrow = TRUE),
  diet_obs = matrix(c(50, 0.6,
                      50, 0.3,
                      40, 0.7),
                    nrow = 3, byrow = TRUE),
  suitMode = c(2L, 0L))

dl <- build_osa_data(dl)

check("diet_obsvec_idx length == n_stomach", length(dl$diet_obsvec_idx) == 2)
check("both stomachs included", all(dl$diet_obsvec_idx >= 0))

diet_rows <- dl$obs_ctl[dl$obs_ctl$type == "diet", ]
# stomach 0: 2 prey + other = 3 bins; stomach 1: 1 prey + other = 2 bins
check("diet obs_ctl rows == 5", nrow(diet_rows) == 5)
check("one last-bin (other prey) per stomach",
      sum(diet_rows$is_last_bin) == 2)
check("stomach_id recorded", setequal(unique(diet_rows$stomach_id), c(0, 1)))

# Recompute stomach 0 counts: props (0.6, 0.3) + other (0.1), offset, normalize, * N_s
v <- c(0.6, 0.3, 1 - 0.9) + 0.00001
expected0 <- v / sum(v) * 50
start0 <- dl$diet_obsvec_idx[1]
check("stomach 0 counts match", isTRUE(all.equal(dl$obsvec[start0 + 1:3], expected0)))

# ---- 2. Regression: jnll_comp identical to baseline (diet ACTIVE) ----------
# msmMode = 1 + suitMode = 4 estimates suitability so the diet likelihood
# (slot 18 == jnll_comp[19, ]) is active, exercising the diet osa_mode == 0 path.
cat("\n== 2. baseline vs new jnll_comp (diet active) ==\n")
set.seed(123)
sim <- make_msm_test_data(years = 1:10)
simData <- sim$data_list
fit_new <- Rceattle::fit_mod(simData, inits = NULL, estimateMode = 3,
                             msmMode = 1, suitMode = 4,
                             fit_control = fit_control(phase = FALSE, verbose = 0))
par <- fit_new$obj$env$last.par.best
rep_new <- fit_new$obj$report(par)
cat("   diet slot (18) nll:", sum(rep_new$jnll_comp[19, ]), "\n")
if (file.exists(file.path("src", "TMB", "ceattle_base.cpp"))) {
  fit_base <- Rceattle::fit_mod(
    simData, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 4,
    fit_control = fit_control(phase = FALSE, verbose = 0,
                              TMBfilename = file.path("src", "TMB", "ceattle_base")))
  rep_base <- fit_base$obj$report(par)
  d <- max(abs(rep_new$jnll_comp - rep_base$jnll_comp))
  cat(sprintf("   max |jnll_comp diff| = %.3e\n", d))
  check("jnll_comp identical to baseline (<1e-8)", d < 1e-8)
}

# ---- 3. Diet OSA runtime (structural) -------------------------------------
cat("\n== 3. diet OSA runtime ==\n")
fit <- Rceattle::fit_mod(simData, inits = NULL, estimateMode = 1,
                         msmMode = 1, suitMode = 4, random_rec = FALSE,
                         fit_control = fit_control(phase = TRUE, verbose = 0))
mg <- tryCatch(max(abs(fit$obj$gr(fit$obj$env$last.par.best))),
               error = function(e) NA)
cat("   obj:", fit$opt$objective, " maxgrad:", mg, "\n")
cat("   obs_ctl types:\n"); print(table(fit$obs_ctl$type))
check("diet positions exist", sum(fit$obs_ctl$type == "diet") > 0)

osa <- suppressWarnings(osa_residuals(fit, types = "diet"))
cat("   diet residuals:", nrow(osa),
    " non-finite:", sum(!is.finite(osa$residual)), "\n")
check("diet residuals present", nrow(osa) > 0)
osa2 <- suppressWarnings(osa_residuals(fit, types = "diet"))
check("deterministic", isTRUE(all.equal(osa$residual, osa2$residual)))
if (isTRUE(mg < 0.1)) check("diet residuals finite (converged)",
                            all(is.finite(osa$residual)))

cat(sprintf("\n==== OVERALL: %s ====\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
