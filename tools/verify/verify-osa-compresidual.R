# tools/verify/verify-osa-compresidual.R
# Do Rceattle's composition OSA residuals agree with an independent
# implementation? Cross-checks against compResidual (Trijoulet et al. 2023,
# Fish. Res. 257:106487), the randomized-quantile implementation by the
# authors of the method, on BS2017SS.
#
# THE RULE THIS HARNESS EXISTS TO ENFORCE: report the RATIO to the ceiling, not
# the bare correlation. compResidual's residuals carry a random draw, so the
# yardstick is their own seed-to-seed correlation -- the best score any method
# could achieve against them. Reading the bare number instead produced a false
# "fleet 2 diverges" conclusion on the first run.
#
# Measured 2026-10-12 (estimateMode = 1, random_rec = FALSE, phase = FALSE):
#
#   fleet  K   ceiling        ours-vs-theirs  ratio  sparse frac
#   1      12  0.973 (0.005)  0.971           1.00   0.10
#   2      25  0.661 (0.029)  0.648           0.98   0.41
#   3      25  NaN            NaN             --     0.19
#   4      12  0.974 (0.005)  0.972           1.00   0.06
#
# So fleet 2's 0.648 is 0.98 OF ITS CEILING, not a divergence -- the ceiling is
# low there because 41% of its cells hold less than one fish, and a sparse
# multinomial is where the random draw dominates. Fleet 3 returns NaN from
# resMulti(), an instability on their side, not a disagreement.
#
# Reports, per fleet: the seed-to-seed ceiling with its sd, our correlation
# against theirs, the ratio of the two, and the fraction of cells below 1.
#
# Caveats. compResidual is external, so it sees one fleet's composition matrix
# in isolation while our OSA conditions on the whole obsvec; and counts must be
# rounded to integers for resMulti(), which is the fractional-count problem in
# miniature. Needs compResidual installed; it is not a package dependency.

suppressMessages(pkgload::load_all(".", quiet = TRUE)); suppressMessages(library(compResidual))
fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
  Rceattle::BS2017SS, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
  fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
osa <- suppressWarnings(Rceattle::osa_residuals(fit, source = "comp")); pear <- attr(osa, "pearson")
mats <- function(flt) {
  po <- pear[pear$fleet == flt, ]
  yrs <- sort(unique(po$year)); bins <- sort(unique(po$age_length_bin))
  o <- p <- matrix(NA_real_, length(bins), length(yrs))
  for (j in seq_along(yrs)) { s <- po[po$year == yrs[j], ]; s <- s[match(bins, s$age_length_bin), ]
    o[, j] <- s$observed * s$sample_size; p[, j] <- s$predicted }
  ok <- colSums(is.na(o)) == 0 & colSums(is.na(p)) == 0
  list(o = round(o[, ok, drop = FALSE]), p = p[, ok, drop = FALSE], yrs = yrs[ok], bins = bins)
}
ours_vec <- function(flt, m) {
  v <- numeric(0)
  for (y in m$yrs) { s <- osa[osa$fleet == flt & osa$year == y, ]; v <- c(v, s$residual[order(s$age_length_bin)]) }
  v
}
for (flt in c(1, 2, 3, 4)) {
  m <- mats(flt)
  ceil <- sapply(1:5, function(s) {
    set.seed(s);       x <- as.numeric(suppressWarnings(compResidual::resMulti(m$o, m$p)))
    set.seed(s + 100); y <- as.numeric(suppressWarnings(compResidual::resMulti(m$o, m$p)))
    suppressWarnings(stats::cor(x, y, use = "complete.obs")) })
  set.seed(7); theirs <- as.numeric(suppressWarnings(compResidual::resMulti(m$o, m$p)))
  ours <- ours_vec(flt, m)
  n <- min(length(theirs), length(ours))
  got <- suppressWarnings(stats::cor(theirs[1:n], ours[1:n], use = "complete.obs"))
  cat(sprintf("C fleet %d (K=%d): ceiling %.3f (sd %.3f) | ours-vs-theirs %.3f | ratio %.2f | sparse frac %.2f\n",
      flt, length(m$bins), mean(ceil, na.rm=TRUE), stats::sd(ceil, na.rm=TRUE), got,
      got / mean(ceil, na.rm=TRUE), mean(m$o < 1)))
}
