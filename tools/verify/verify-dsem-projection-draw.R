# Does sample_rec() draw a DSEM's projection from the process that was fitted?
#
# tests/testthat/test-dsem-projection-draw.R pins the conditional-draw
# mathematics against precisions built by hand. This script is the other half:
# it checks that the MODEL'S OWN reported precision and parameter map reach that
# code correctly, on a real fit, in the configuration the feature exists for --
# a climate covariate whose values are KNOWN over the projection.
#
# That configuration is the one that fails silently. Under family = "fixed" a
# covariate column of x_tj is the environmental data pinned by the map; if the
# draw treats those known future values as unknown and draws them, projected
# recruitment stops responding to the environment and quietly returns the
# no-covariate answer. Nothing errors, and the numbers look plausible.
#
# Run: export PATH=/usr/bin:$PATH; Rscript tools/verify/verify-dsem-projection-draw.R

suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = FALSE))

ok <- function(label, pass, detail = "") {
  cat(sprintf("%-4s %s%s\n", if (isTRUE(pass)) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  [", detail, "]") else ""))
  invisible(isTRUE(pass))
}
results <- c()

# ---- a fit whose covariate is known through the projection --------------------
d <- Rceattle::BS2017SS
d$projyr <- d$endyr + 8
all_yrs <- d$styr:d$projyr
set.seed(1)
# A covariate with a clear projection-period swing, so a projection that
# responds to it is visibly different from one that does not.
env <- as.numeric(scale(sin(seq_along(all_yrs) / 3)))
d$env_data <- data.frame(Year = all_yrs, BT = env)

sem <- "recdevs1 <-> recdevs1, 0, sigmaR1, 0.6
recdevs2 <-> recdevs2, 0, sigmaR2, 0.6
recdevs3 <-> recdevs3, 0, sigmaR3, 0.6
recdevs1 -> recdevs1, 1, rho1, 0.4
BT <-> BT, 0, sdBT, 1
BT -> recdevs1, 0, bBT, 0.2"
fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
fit_one <- function(dat, proj_mean_rec = FALSE) {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    dat, inits = NULL, file = NULL, estimateMode = 0, random_rec = TRUE,
    msmMode = 0, recFun = Rceattle::build_srr(proj_mean_rec = proj_mean_rec),
    dsem = Rceattle::build_DSEM(sem = sem, family = "fixed",
                                estimate_projection = TRUE),
    fit_control = fc)))
}
fit <- fit_one(d)

n_t <- nrow(fit$estimated_params$dsem_x_tj)
n_j <- ncol(fit$estimated_params$dsem_x_tj)
n_h <- fit$data_list$endyr - fit$data_list$styr + 1L
pj  <- (n_h + 1L):n_t
rcol <- fit$dsem$tmb_inputs$data$rec_dev_col[1] + 1L
mp <- matrix(is.na(as.numeric(fit$map$mapList$dsem_x_tj)), n_t, n_j)
ecol <- setdiff(which(colSums(mp) > 0), rcol)

results <- c(results, ok(
  "the covariate column is pinned over the projection",
  length(ecol) == 1L && all(mp[pj, ecol]) && !any(mp[pj, rcol]),
  sprintf("env col %s fixed in %d/%d proj yrs; recdev col %d free",
          paste(ecol, collapse = ","), sum(mp[pj, ecol]), length(pj), rcol)))

# ---- 1. the known future environment is preserved, never redrawn -------------
X0 <- as.matrix(fit$estimated_params$dsem_x_tj)
set.seed(11)
draws <- lapply(1:25, function(i)
  as.matrix(suppressWarnings(Rceattle::sample_rec(
    fit, update_model = FALSE))$estimated_params$dsem_x_tj))
env_moved <- max(vapply(draws, function(m) max(abs(m[, ecol] - X0[, ecol])), 0))
rec_moved <- stats::sd(vapply(draws, function(m) m[n_t, rcol], 0))
results <- c(results, ok(
  "draws leave the known covariate alone but move recruitment",
  env_moved == 0 && rec_moved > 0.1,
  sprintf("env max|d| = %.2e, terminal recdev state sd = %.3f",
          env_moved, rec_moved)))

# ---- 2. projected recruitment RESPONDS to the future covariate ---------------
# The decisive check, and it holds the FIT fixed: same Q, same paths, same
# hindcast -- only the stored future covariate changes. A projection that
# conditions on those values must move with them; one that draws them away
# cannot see the change at all, so the buggy version returns an identical
# conditional mean. Refitting instead would confound this with the parameters
# moving.
mean_state <- function(f) {
  m <- .dsem_draw_projection(f, sample = FALSE)
  as.numeric(m[pj, rcol])
}
s_lo <- mean_state(fit)
fit_hi <- fit
fit_hi$estimated_params$dsem_x_tj[pj, ecol] <-
  fit$estimated_params$dsem_x_tj[pj, ecol] + 3
s_hi <- mean_state(fit_hi)
b <- fit$estimated_params$dsem_beta_z[
  fit$dsem$sem_full$parameter[fit$dsem$sem_full$name == "bBT"]]
rho <- fit$estimated_params$dsem_beta_z[
  fit$dsem$sem_full$parameter[fit$dsem$sem_full$name == "rho1"]]
# The lag-0 path gives an immediate response of b * 3 in the first projection
# year; the lag-1 self-path then accumulates it toward b * 3 / (1 - rho).
results <- c(results, ok(
  "a shifted future covariate shifts the projected states",
  abs(s_hi[1] - s_lo[1] - b * 3) < 1e-6 &&
    abs(s_hi[length(s_hi)] - s_lo[length(s_lo)]) > abs(b * 3),
  sprintf("bBT = %.3f, rho = %.3f; year 1 shift = %.4f (expected %.4f), ",
          b, rho, s_hi[1] - s_lo[1], b * 3)))

# ---- 3. sample_rec = FALSE is MEAN recruitment, not the median ---------------
# The standard path projects at the arithmetic mean of hindcast recruitment, so
# the DSEM path has to mean the same thing: R = R0 * exp(conditional mean). If
# the bias correction is not added back, R lands at the lognormal MEDIAN --
# roughly exp(margvar/2) low, ~26% on a typical fit, in a deterministic
# projection that sets the ABC.
det <- suppressWarnings(Rceattle::sample_rec(fit, sample_rec = FALSE,
                                             update_model = FALSE))
mv <- fit$quantities$dsem_margvar_tj
bias <- as.numeric(fit$data_list$bias_adjust_proc)
got <- as.numeric(det$estimated_params$rec_dev[1, pj])
want <- s_lo
results <- c(results, ok(
  "sample_rec = FALSE returns the conditional mean deviation",
  max(abs(got - want)) < 1e-8,
  sprintf("max|d| = %.2e (median convention would be off by %.3f)",
          max(abs(got - want)), bias * mean(mv[pj, rcol]) / 2)))

# ---- 4. rec_dev matches what the template will derive ------------------------
drawn <- suppressWarnings(Rceattle::sample_rec(fit, update_model = FALSE))
expect <- as.numeric(drawn$estimated_params$dsem_x_tj[pj, rcol]) -
  bias * as.numeric(mv[pj, rcol]) / 2
results <- c(results, ok(
  "rec_dev mirrors the template's x - bias * margvar / 2",
  max(abs(as.numeric(drawn$estimated_params$rec_dev[1, pj]) - expect)) < 1e-12))

# ---- 5. the draw survives the rebuild ---------------------------------------
set.seed(3); a <- suppressWarnings(Rceattle::sample_rec(fit, update_model = FALSE))
set.seed(3); b2 <- suppressWarnings(Rceattle::sample_rec(fit, update_model = TRUE))
results <- c(results, ok(
  "update_model = TRUE keeps the drawn states",
  max(abs(as.matrix(b2$estimated_params$dsem_x_tj)[pj, ] -
          as.matrix(a$estimated_params$dsem_x_tj)[pj, ])) < 1e-12 &&
    !is.null(b2$quantities$R)))

# ---- 6. the two contradictory-by-default settings are refused ----------------
fit_mean <- fit_one(d, proj_mean_rec = TRUE)
e1 <- tryCatch({Rceattle::sample_rec(fit_mean, update_model = FALSE); ""},
               error = conditionMessage)
results <- c(results, ok(
  "proj_mean_rec = TRUE is refused, not silently deterministic",
  grepl("proj_mean_rec = FALSE", e1)))

d_noproj <- d
fit_np <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
  d_noproj, inits = NULL, file = NULL, estimateMode = 0, random_rec = TRUE,
  msmMode = 0, recFun = Rceattle::build_srr(proj_mean_rec = FALSE),
  dsem = Rceattle::build_DSEM(sem = sem, family = "fixed",
                              estimate_projection = FALSE),
  fit_control = fc)))
e2 <- tryCatch({Rceattle::sample_rec(fit_np, update_model = FALSE); ""},
               error = conditionMessage)
results <- c(results, ok(
  "estimate_projection = FALSE is refused, not extrapolated",
  grepl("estimate_projection = TRUE", e2)))

cat(sprintf("\n%d/%d checks passed\n", sum(results), length(results)))
if (!all(results)) quit(status = 1)
