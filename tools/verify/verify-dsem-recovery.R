# tools/verify/verify-dsem-recovery.R
#
# Does the DSEM's process density RECOVER a known process, and is the marginal
# variance it reports correct?
#
# Everything else on this branch checks EQUIVALENCE -- that a naive (IID) DSEM
# reproduces the non-DSEM model. Those checks are blind to a lagged sem by
# construction: a lag is exactly where the two models stop being the same. So
# nothing has ever verified that the GMRF estimates rho, a covariate effect, or
# sigma correctly. `sim_mod(process = TRUE)` does redraw a DSEM's field, but it
# draws from the fit's OWN reported precision, so a recovery harness built on it
# could only ever confirm the model agrees with itself. This drives the density
# from a process specified independently, which is what makes it a check.
#
# This drives calculate_dsem() directly through tools/verify/dsem_standalone.cpp:
#
#   1. RECOVERY. Simulate latent states from a known AR1-with-covariate process,
#      hold them fixed, and estimate the path coefficients. rho, the covariate
#      effect and sigma must come back.
#   2. MARGINAL VARIANCE. The value the recruitment bias correction is built on.
#      For a first-order self-path it must approach sigma^2 / (1 - rho^2) away
#      from the series edges, and equal sigma^2 exactly when rho = 0.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-dsem-recovery.R

stopifnot(requireNamespace("dsem", quietly = TRUE))
suppressMessages(pkgload::load_all(".", quiet = TRUE))

ok <- TRUE
report <- function(label, pass, detail = "") {
  cat(sprintf("  %-44s %s %s\n", label,
              if (is.na(pass)) "SKIP" else if (pass) "PASS" else "FAIL", detail))
  if (identical(pass, FALSE)) ok <<- FALSE
}

dsem_data <- function(d) {
  list(options      = as.integer(d$options),
       RAM          = matrix(as.integer(as.matrix(d$RAM)), ncol = 6),
       RAMstart     = as.numeric(d$RAMstart),
       familycode_j = as.integer(d$familycode_j),
       linkcode_j   = as.integer(d$linkcode_j),
       sigmastart_j = as.integer(d$sigmastart_j),
       eps_tj       = array(as.numeric(d$eps_tj), dim = dim(d$eps_tj)),
       y_tj         = array(as.numeric(d$y_tj),   dim = dim(d$y_tj)),
       obs_idx      = if (is.null(d$obs_idx))   integer(0) else as.integer(d$obs_idx),
       unobs_idx    = if (is.null(d$unobs_idx)) integer(0) else as.integer(d$unobs_idx),
       # The cells the model is given, matching what ceattle.cpp builds: a
       # finite env_data value, per cell, so a covariate with a gap is known
       # where it was observed and random where it was not. The recruitment
       # columns and the deterministic unobs_idx cells are excluded, as there.
       cond_k       = { y  <- as.matrix(d$y_tj)
                        ck <- as.integer(is.finite(as.vector(y)))
                        ck[as.vector(col(y)) %in%
                           (as.integer(d$rec_dev_col) + 1L)] <- 0L
                        if (length(d$unobs_idx))
                          ck[as.integer(d$unobs_idx) + 1L] <- 0L
                        ck })
}

# ---- 1. recovery of a lagged, covariate-driven process ----------------------
# Inputs come from Rceattle's own build_dsem_objects(), not from a bare
# dsem::dsem() call: it is the path a real fit takes, and calling dsem::dsem()
# directly on a long simulated series segfaults inside its MakeADFun in this
# environment, which is a dsem-package problem and not what this is testing.
cat("\n=== 1. parameter recovery, AR1 + covariate ===\n")
set.seed(42)
n_t   <- 300                       # long series: this is about bias, not variance
RHO   <- 0.60
BETA  <- 0.45
SIGMA <- 0.55

styr <- 1L
temp <- as.numeric(scale(stats::rnorm(n_t)))
x <- numeric(n_t)
x[1] <- stats::rnorm(1, 0, SIGMA / sqrt(1 - RHO^2))     # stationary start
for (t in 2:n_t) x[t] <- RHO * x[t - 1] + BETA * temp[t] + stats::rnorm(1, 0, SIGMA)

# Minimal single-species data_list: only the DSEM builder is exercised.
dl <- list(nspp = 1L, styr = styr, endyr = styr + n_t - 1L,
           projyr = styr + n_t - 1L, spnames = "sp1",
           # sigma_rec is only a starting value here -- the sem names sigma
           # explicitly, and the point of the test is that it is ESTIMATED back
           # to SIGMA from a deliberately wrong start.
           sigma_rec = 1, random_rec = TRUE,
           env_data = data.frame(Year = styr:(styr + n_t - 1L), temp = temp))
sem <- "
  recdevs1 -> recdevs1, 1, rho,   0.1
  temp     -> recdevs1, 0, beta,  0.1
  recdevs1 <-> recdevs1, 0, sigma, 1
"
# Named arguments: the signature is (dsem_settings, debug, data_list), so a
# positional call passes the data list as the settings and fails deep inside.
built <- Rceattle:::build_dsem_objects(
  dsem_settings = Rceattle::build_DSEM(sem = sem, family = "fixed"),
  data_list = dl)
dat  <- dsem_data(built$tmb_inputs$data)
pars <- built$tmb_inputs$parameters
sem_full <- built$sem_full

dir.create("dev", showWarnings = FALSE)
file.copy("tools/verify/dsem_standalone.cpp", "dev/dsem_standalone.cpp", overwrite = TRUE)
file.copy("src/TMB/dsem.hpp", "dev/dsem.hpp", overwrite = TRUE)
cat("compiling the standalone DSEM model...\n")
TMB::compile("dev/dsem_standalone.cpp", framework = "TMBad", flags = "-O1")
dyn.load(TMB::dynlib("dev/dsem_standalone"))

# Hold the latent states at the simulated truth and estimate only the paths.
# That isolates the density: any failure is the GMRF, not the assessment around
# it. Column order is the builder's: recdevs first, then the covariates.
rec_col <- built$tmb_inputs$data$rec_dev_col[1] + 1L
xt <- pars$x_tj
xt[, rec_col] <- x
for (j in seq_len(ncol(xt))) if (j != rec_col) xt[, j] <- temp
pars$x_tj <- xt
map <- list(x_tj = factor(rep(NA, length(pars$x_tj))))
for (nm in c("mu_j", "delta0_j", "lnsigma_z")) {
  if (length(pars[[nm]])) map[[nm]] <- factor(rep(NA, length(pars[[nm]])))
}

obj <- TMB::MakeADFun(data = dat, parameters = pars, map = map,
                      DLL = "dsem_standalone", silent = TRUE)
opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                     control = list(iter.max = 1000, eval.max = 1000))

est <- setNames(rep(NA_real_, 3), c("rho", "beta", "sigma"))
for (p in names(est)) {
  k <- sem_full$parameter[sem_full$name == p]
  if (length(k) == 1L && k >= 1L) est[[p]] <- opt$par[k]
}
truth <- c(rho = RHO, beta = BETA, sigma = SIGMA)
cat("  estimate:", paste(sprintf("%s=%.4f", names(est), est), collapse = "  "), "\n")
cat("  truth   :", paste(sprintf("%s=%.4f", names(truth), truth), collapse = "  "), "\n")
for (p in names(truth)) {
  e <- abs(est[[p]]) ; t0 <- truth[[p]]        # sigma's sign is not identified
  report(paste0("recovers ", p), !is.na(e) && abs(e - t0) / t0 < 0.10,
         sprintf("%.4f vs %.4f (%+.1f%%)", e, t0, 100 * (e / t0 - 1)))
}
report("optimizer converged", opt$convergence == 0)

# ---- 2. the marginal variance the bias correction uses ----------------------
cat("\n=== 2. marginal variance vs closed form ===\n")
set_beta <- function(vals) {
  p2 <- pars
  for (nm in names(vals)) {
    k <- sem_full$parameter[sem_full$name == nm]
    if (length(k) == 1L && k >= 1L) p2$beta_z[k] <- vals[[nm]]
  }
  o <- TMB::MakeADFun(data = dat, parameters = p2, map = map,
                      DLL = "dsem_standalone", silent = TRUE)
  o$fn(o$par); o$report()
}
mid <- round(n_t / 2)                                   # away from the edges
r1 <- set_beta(list(rho = RHO, beta = BETA, sigma = SIGMA))
analytic <- SIGMA^2 / (1 - RHO^2)
cat(sprintf("  mid-series margvar = %.6f   sigma^2/(1-rho^2) = %.6f\n",
            r1$margvar_tj[mid, rec_col], analytic))
report("AR1 marginal variance matches closed form",
       abs(r1$margvar_tj[mid, rec_col] - analytic) / analytic < 0.01,
       sprintf("%+.2f%%", 100 * (r1$margvar_tj[mid, rec_col] / analytic - 1)))

# rho = 0 must collapse to sigma^2 exactly: that is the case where the naive
# -sigma^2/2 correction is right, and where a DSEM must equal the standard path.
r0 <- set_beta(list(rho = 0, beta = 0, sigma = SIGMA))
report("rho = 0 gives margvar = sigma^2 exactly",
       abs(r0$margvar_tj[mid, rec_col] - SIGMA^2) < 1e-8,
       sprintf("%.8f vs %.8f", r0$margvar_tj[mid, rec_col], SIGMA^2))

cat("\n", if (ok) "ALL CHECKS PASSED" else "SOME CHECKS FAILED", "\n", sep = "")
quit(status = if (ok) 0L else 1L)
