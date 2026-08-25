# Verify the linkage random-effect simulator (ceattle.cpp section 5.12b).
#
# Time-varying recruitment, M and growth are all written as a random linkage on
# the process, so one block draws all of them. The draw switches on the group's
# covariance structure and the gate switches on its process, and those are
# independent -- so the structures are checked where they build cheaply (M) and
# the processes are checked separately.
#
# Two things are deliberate:
#  * sigma and rho are PINNED in the parameter vector. These fixtures hold
#    little information about a variance, so a fitted sigma collapses toward
#    zero and a moment check on it would pass whatever the simulator did.
#  * Models are built (estimateMode = 3), not optimized. The simulator draws
#    from the density at the parameters it is given, so an optimum buys nothing
#    -- and a year-varying R0 is confounded with rec_dev, so that fit grinds.
#
# Run: export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript tools/verify/verify-sim-linkage-re.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, e)

CTL   <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
NDRAW <- 400L
fails <- character(0)

report <- function(label, got, want, tol) {
  ok <- is.finite(got) && abs(got - want) <= tol
  cat(sprintf("  %-40s %9.4f  (target %8.4f, tol %6.4f)  %s\n",
              label, got, want, tol, if (ok) "ok" else "FAIL"))
  if (!ok) fails <<- c(fails, label)
  invisible(ok)
}

# Correlation between adjacent slots taken ACROSS draws, not along the series.
# The within-series estimate is biased by -1/(n-1), which is -0.14 at these
# fixture lengths -- big enough to swamp what is being measured.
adjacent_cor <- function(m) {
  mean(vapply(2:ncol(m), function(t) cor(m[, t], m[, t - 1]), numeric(1)))
}

# The deviations a draw must leave alone, read from the parameter vector the
# simulate pass starts from. NOT from $estimated_params: obj$fn() runs the inner
# Laplace step, so the random effects sit at their conditional modes, which is
# what an undrawn group comes back as. The `process = FALSE` check below is what
# confirms this is the right reference.
undrawn_re <- function(fit, sigma) {
  obj <- Rceattle:::.sim_obj(fit)
  obj$fn(obj$par)
  par <- obj$env$last.par.best
  par[names(par) == "log_sigma_linkage"] <- log(sigma)
  as.numeric(par[names(par) == "beta_linkage_re"])
}

# Draw the linkage REs NDRAW times with sigma/rho pinned. Returns a
# draw x slot matrix of deviations in slot (= time) order.
draw_re <- function(fit, state, sigma = NULL, rho = NULL, period = c(1, 0),
                    n = NDRAW, seed = 1L) {
  obj <- Rceattle:::.sim_obj(fit)
  obj$fn(obj$par)                        # populate last.par.best (never optimized)
  par <- obj$env$last.par.best
  if (!is.null(sigma)) par[names(par) == "log_sigma_linkage"] <- log(sigma)
  if (!is.null(rho))   par[names(par) == "trans_rho_linkage"] <- atanh(rho)
  obj$env$last.par.best <- par
  set.seed(seed)
  out <- NULL
  for (i in seq_len(n)) {
    v <- as.numeric(Rceattle:::.sim_draw(obj, state = state, period = period)$beta_linkage_re_sim)
    if (is.null(out)) out <- matrix(NA_real_, n, length(v))
    out[i, ] <- v
  }
  out
}

state_of <- function(...) Rceattle:::.sim_state_codes(c(...))

m_fit <- function(form) {
  Rceattle::fit_mod(
    data_list = e$make_test_data(), estimateMode = 3, msmMode = 0, fit_control = CTL,
    M1Fun = Rceattle::build_M1(
      M1_model = 1,
      linkages = list(M1 = Rceattle::linkage_spec(form, by = ~ species))))
}

# ---- structures ------------------------------------------------------------

# IID: sum of dnorm(re, 0, sigma) -- every slot independent.
cat("\n[1] IID `(1 | Year)`  -- independent N(0, sigma) per year\n")
SIG <- 0.4
d1 <- draw_re(suppressWarnings(m_fit(~ (1 | Year))), state_of("M"), sigma = SIG)
report("realised sd",               sd(as.numeric(d1)),   SIG, 0.03)
report("realised mean",             mean(as.numeric(d1)), 0,   0.03)
report("adjacent-slot correlation", adjacent_cor(d1),     0,   0.09)

# rw: dnorm on FIRST DIFFERENCES. Slot 0 is a level the density never sees, so
# it is not identified, the map pins it, and the draw must leave it alone.
cat("\n[2] rw `rw(1 | Year)`  -- N(0, sigma) on first differences\n")
SIG_RW <- 0.3
fit_rw <- suppressWarnings(m_fit(~ rw(1 | Year)))
d2 <- draw_re(fit_rw, state_of("M"), sigma = SIG_RW)
dd <- t(apply(d2, 1, diff))
report("sd of first differences",   sd(as.numeric(dd)),   SIG_RW, 0.03)
report("mean of first differences", mean(as.numeric(dd)), 0,      0.03)
report("slot 0 spread across draws", sd(d2[, 1]),         0,      1e-12)
# A random walk's variance grows linearly in the step index; an IID draw's does
# not. This is what catches drawing the wrong structure.
report("var(last - first) / (n-1)",
       var(d2[, ncol(d2)] - d2[, 1]) / (ncol(d2) - 1), SIG_RW^2, 0.02)

# ar1: SCALE(AR1(rho), sigma), so sigma is the MARGINAL sd. Reading it as the
# innovation sd instead inflates the draw by 1/sqrt(1-rho^2).
cat("\n[3] ar1 `ar1(1 | Year)` -- marginal sd sigma, correlation rho\n")
SIG_A <- 0.5; RHO_A <- 0.6
d3 <- draw_re(suppressWarnings(m_fit(~ ar1(1 | Year))), state_of("M"),
              sigma = SIG_A, rho = RHO_A)
report("realised MARGINAL sd",      sd(as.numeric(d3)),   SIG_A, 0.04)
report("realised mean",             mean(as.numeric(d3)), 0,     0.05)
report("adjacent-slot correlation", adjacent_cor(d3),     RHO_A, 0.09)
cat(sprintf("  (reading sigma as the innovation sd would give sd %.4f)\n",
            SIG_A / sqrt(1 - RHO_A^2)))

# ---- processes -------------------------------------------------------------
# Growth is the process with no working random effect before this: the legacy
# log_growth_par_devs array carries no density, so a random linkage is the only
# way to make growth time-varying.
cat("\n[4] Each population-dynamics process draws\n")
SIG_P <- 0.35

fit_g <- Rceattle::fit_mod(
  data_list = e$make_test_data(growth = "vonBertalanffy"),
  estimateMode = 3, msmMode = 0, fit_control = CTL,
  growthFun = Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(K = Rceattle::linkage_spec(~ (1 | Year), by = ~ species))))

fit_r <- Rceattle::fit_mod(
  data_list = e$make_test_data(), estimateMode = 3, msmMode = 0, fit_control = CTL,
  recFun = Rceattle::build_srr(
    linkages = list(R0 = Rceattle::linkage_spec(~ (1 | Year), by = ~ species))))

for (p in list(list("growth", fit_g), list("recruitment", fit_r),
               list("M", suppressWarnings(m_fit(~ (1 | Year)))))) {
  dp <- draw_re(p[[2]], state_of(p[[1]]), sigma = SIG_P)
  report(paste0(p[[1]], ": realised sd"),   sd(as.numeric(dp)),   SIG_P, 0.04)
  report(paste0(p[[1]], ": realised mean"), mean(as.numeric(dp)), 0,     0.04)
}

# ---- gating ----------------------------------------------------------------
# simulate_state slots ARE the linkage process codes. Slots 3/4 (q, selectivity)
# were crossed against them in an earlier revision; this is the check for that.
cat("\n[5] A draw is confined to the process asked for (growth fixture)\n")
ref <- undrawn_re(fit_g, SIG_P)

check_gate <- function(state, label, should_move, period = c(1, 0)) {
  mx <- max(abs(sweep(draw_re(fit_g, state, sigma = SIG_P, period = period, n = 20L),
                      2, ref)))
  ok <- if (should_move) mx > 0 else mx == 0
  cat(sprintf("  %-40s max|draw - fitted| = %.3e  %s\n", label, mx,
              if (ok) "ok" else "FAIL"))
  if (!ok) fails <<- c(fails, label)
}

check_gate(state_of(),               "process = FALSE",         FALSE)
check_gate(state_of("M"),            "process = 'M'",           FALSE)
check_gate(state_of("recruitment"),  "process = 'recruitment'", FALSE)
check_gate(state_of("catchability"), "process = 'catchability'",FALSE)
check_gate(state_of("selectivity"),  "process = 'selectivity'", FALSE)
check_gate(state_of("growth"),       "process = 'growth' DOES draw", TRUE)
# The density covers the fitted window, so the projection period draws nothing.
check_gate(state_of("growth"), "simulate_period(0) = 0", FALSE, period = c(0, 0))

cat("\n")
if (length(fails)) {
  cat("FAIL (", length(fails), "):\n", sep = "")
  for (f in unique(fails)) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("PASS: linkage random effects simulate correctly for recruitment, M and growth\n")
