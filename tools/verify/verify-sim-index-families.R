# tools/verify/verify-sim-index-families.R
# Does sim_mod() draw each survey from its OWN index likelihood? Every fleet was
# once drawn as an independent lognormal whatever Index_distribution said, so an
# MVN fleet was simulated on the wrong scale, with the wrong spread, and with no
# correlation -- and nothing errored; self_test() simply reported recovery
# against a data-generating process the likelihood never assumed. This is the
# harness that caught it, kept as the gate for the index stage of the migration
# of simulation into the template's SIMULATE blocks. TMB's MVNORM_t simulates
# natively but WHAM has no covariance survey likelihood, so there is no worked
# example to check against -- these moments are the check.
#
# Correlation is measured ACROSS replicates, between years. Measuring it within
# one replicate subtracts the sample mean, which under compound symmetry absorbs
# the common factor and returns about -1/(n-1) whatever rho really is; that
# produced a convincing -0.06 against a true 0.6 and nearly sent a correct fix
# back for rework.
#
# Expected: Lognormal log-sd 0.1 / corr 0; Normal natural sd 20 / corr 0;
#           MVN and MVNORM natural sd 20 / corr 0.6.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-sim-index-families.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, e)

nyrs <- 20
rho <- 0.6
sdc <- 20

build <- function(dist, absolute_sd = FALSE) {
  dat <- e$make_test_data(nyrs = nyrs, nages = 5, seed = 123)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- dist
  dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
  if (dist %in% c("MVN", "MVNORM")) {
    Rho <- matrix(rho, nyrs, nyrs)
    diag(Rho) <- 1
    dat$index_cov <- list(
      Survey = diag(rep(sdc, nyrs)) %*% Rho %*% diag(rep(sdc, nyrs)))
  }
  if (absolute_sd) {
    dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sdc
  }
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
}

for (dist in c("Lognormal", "MVN", "MVNORM", "Normal")) {
  abs_sd <- dist == "Normal"
  fit <- build(dist, absolute_sd = abs_sd)
  srv <- fit$data_list$index_data$Fleet_name == "Survey"
  hat <- fit$quantities$index_hat[srv]

  set.seed(42)
  reps <- replicate(600, {
    s <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))
    s$index_data$Observation[srv]
  })

  ## Correlation must be measured ACROSS replicates, between years. Measuring it
  ## within one replicate subtracts the sample mean, which for compound symmetry
  ## absorbs the common factor and returns about -1/(n-1) whatever rho is.
  nat <- dist %in% c("MVN", "MVNORM", "Normal")
  dev <- if (nat) reps - hat else log(reps) - log(hat)
  C <- stats::cor(t(dev))                      # year x year, over replicates
  off <- mean(C[upper.tri(C)])

  cat("\n=========== Index_distribution =", dist, "===========\n")
  if (nat) {
    cat("  expected  : natural-scale sd", sdc, " pairwise correlation",
        if (dist == "Normal") 0 else rho, "\n")
  } else {
    cat("  expected  : log-scale sd",
        unique(signif(fit$quantities$log_index_sd[srv], 4))[1],
        " pairwise correlation 0\n")
  }
  cat("  simulated : sd", signif(mean(apply(dev, 1, sd)), 4),
      " pairwise correlation", signif(off, 3), "\n")
}
