# Verify the population-dynamics process-error draws (ceattle.cpp section 5.13)
# and the two M random-effect defects fixed alongside them.
#
# None of this is covered by /golden-check: no golden model uses M1_re, so those
# fits are bit-identical either way. These are the regression net.
#
# Run: export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript tools/verify/verify-sim-process-error.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, e)

CTL   <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
fails <- character(0)

report <- function(label, got, want, tol) {
  ok <- is.finite(got) && abs(got - want) <= tol
  cat(sprintf("  %-44s %9.4f  (target %8.4f, tol %6.4f)  %s\n",
              label, got, want, tol, if (ok) "ok" else "FAIL"))
  if (!ok) fails <<- c(fails, label)
}
ok_flag <- function(label, ok, detail = "") {
  cat(sprintf("  %-44s %-22s %s\n", label, detail, if (ok) "ok" else "FAIL"))
  if (!ok) fails <<- c(fails, label)
}

# Draw n times with named parameters pinned; returns a list of reported arrays.
proc_draws <- function(fit, state, what, n, pin = list(), period = c(1, 0), seed = 1L) {
  obj <- Rceattle:::.sim_obj(fit)
  obj$fn(obj$par)
  par <- obj$env$last.par.best
  for (nm in names(pin)) {
    # last.par.best carries only the FREE parameters, so a mapped-off name
    # matches nothing and the assignment is a silent no-op -- which reads as a
    # passing moment check against whatever the model already had.
    if (!any(names(par) == nm))
      stop("cannot pin '", nm, "': not a free parameter of this fixture", call. = FALSE)
    par[names(par) == nm] <- pin[[nm]]
  }
  obj$env$last.par.best <- par
  set.seed(seed)
  lapply(seq_len(n), function(i)
    Rceattle:::.sim_draw(obj, state = state, period = period)[[what]])
}

goa_m_fit <- function(m1_re) {
  # Species 2 (arrowtooth) is the two-sex stock; M1_model = 1 makes M
  # sex-invariant, so build_map() maps the male deviations onto the female ones.
  suppressWarnings(Rceattle::fit_mod(
    data_list = Rceattle::GOA2018SS, estimateMode = 3, msmMode = 0, fit_control = CTL,
    M1Fun = Rceattle::build_M1(M1_model = 1, M1_re = m1_re)))
}

ST_M <- Rceattle:::.sim_state_codes("M")

# ---- 1. sexes sharing one M deviation get the same draw ---------------------
cat("\n[1] Sex-invariant M (M1_model = 1): both sexes must get the same draw\n")
for (mode in 1:6) {
  fit <- goa_m_fit(c(0, mode, 0))
  d <- proc_draws(fit, ST_M, "log_M1_dev_sim", n = 3L,
                  pin = list(M1_dev_log_sd = log(0.4)))[[1]]
  sp <- 2                                    # arrowtooth
  gap <- max(abs(d[sp, 1, , ] - d[sp, 2, , ]))
  moved <- max(abs(d[sp, 1, , ]))
  ok_flag(sprintf("M1_re = %d: male == female", mode),
          gap == 0 && moved > 0,
          sprintf("gap %.2e, drew %.3f", gap, moved))
}

# ---- 2. M1_re = 3/6 has one free deviation per age-year cell ----------------
cat("\n[2] 2D M random effects: one free parameter per age x year cell\n")
for (mode in c(3L, 6L)) {
  fit <- goa_m_fit(c(0, mode, 0))
  # $map is list(mapFactor, mapList) -- fit$map$log_M1_dev is NULL, and reading
  # it that way makes any is.na() check vacuously true.
  mp <- fit$map$mapList$log_M1_dev
  stopifnot(!is.null(mp))
  # Species 2 is the one carrying M1_re in this configuration.
  nages <- Rceattle::GOA2018SS$nages[2]
  nyrs  <- Rceattle::GOA2018SS$endyr - Rceattle::GOA2018SS$styr + 1L
  slice <- mp[2, 1, 1:nages, 1:nyrs]
  n_free <- length(unique(stats::na.omit(as.integer(slice))))
  ok_flag(sprintf("M1_re = %d: distinct free deviations", mode),
          n_free == nages * nyrs,
          sprintf("%d of %d cells", n_free, nages * nyrs))
}

# ---- 3. the 2D-AR1 correlations act on the dimensions they name -------------
cat("\n[3] 2D-AR1 M: rho_age correlates ages, rho_year correlates years\n")
# SEPARABLE(f, g) puts f on the OUTERMOST array dimension and g on the fastest.
# The array is (age, year), so year must be passed first. Written the other way
# round the two correlations are silently exchanged.
RHO_A <- 0.8; RHO_Y <- 0.2; SIG <- 0.3
fit6 <- goa_m_fit(c(0, 6L, 0))
dd <- proc_draws(fit6, ST_M, "log_M1_dev_sim", n = 200L,
                 pin = list(M1_dev_log_sd = log(SIG),
                            M1_rho = c(atanh(RHO_A), atanh(RHO_Y))))
nages <- Rceattle::GOA2018SS$nages[2]
nyrs  <- Rceattle::GOA2018SS$endyr - Rceattle::GOA2018SS$styr + 1L
# Correlate adjacent cells across independent draws (unbiased at this n).
cor_along <- function(along) {
  v <- vapply(dd, function(x) as.numeric(x[2, 1, 1:nages, 1:nyrs]), numeric(nages * nyrs))
  m <- array(v, c(nages, nyrs, length(dd)))
  if (along == "age") {
    mean(vapply(2:nages, function(a) mean(vapply(1:nyrs, function(y)
      cor(m[a, y, ], m[a - 1, y, ]), numeric(1))), numeric(1)))
  } else {
    mean(vapply(2:nyrs, function(y) mean(vapply(1:nages, function(a)
      cor(m[a, y, ], m[a, y - 1, ]), numeric(1))), numeric(1)))
  }
}
report("adjacent-AGE correlation",  cor_along("age"),  RHO_A, 0.08)
report("adjacent-YEAR correlation", cor_along("year"), RHO_Y, 0.08)

# ---- 4. rec_dev / init_dev moments -----------------------------------------
cat("\n[4] Recruitment and initial-age deviations\n")
d <- e$make_test_data(); d$initMode <- 2L      # NonEquilibrium: init_dev penalized
fit_r <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0, fit_control = CTL)
# R_log_sd is fixed in this fixture, so take the target from the model rather
# than pinning it -- the density uses whatever R_sd the model holds.
R_SD <- as.numeric(fit_r$quantities$R_sd)[1]
cat(sprintf("  (model R_sd = %.4f)\n", R_SD))
ST_R <- Rceattle:::.sim_state_codes("recruitment")
# Hindcast columns only. rec_dev spans styr:projyr, but the density -- and so
# the draw -- covers the fitted window; projection recruitment comes from the
# harvest control rule or sample_rec(). Averaging the undrawn projection years
# in pulls the sd toward zero and looks like a biased draw.
n_hind <- d$endyr - d$styr + 1L
rd <- proc_draws(fit_r, ST_R, "rec_dev_sim", n = 400L)
rv <- unlist(lapply(rd, function(x) x[1, seq_len(n_hind)]))
report("rec_dev sd",   stats::sd(rv), R_SD,        0.05)
report("rec_dev mean", mean(rv),      -R_SD^2 / 2, 0.05)   # lognormal bias correction
# The projection years must come back untouched.
proj_moved <- max(vapply(rd, function(x)
  max(abs(x[1, (n_hind + 1):ncol(x)] - rd[[1]][1, (n_hind + 1):ncol(x)])), numeric(1)))
ok_flag("projection years not drawn", proj_moved == 0,
        sprintf("spread %.2e", proj_moved))
# init_dev is dimensioned nages but only ages 1..nages-1 are drawn; the spare
# slot holds a -999 sentinel that would swamp any moment taken over the array.
n_init <- d$nages[1] - 1L
iv <- unlist(lapply(proc_draws(fit_r, ST_R, "init_dev_sim", n = 400L),
                    function(x) x[1, seq_len(n_init)]))
report("init_dev sd",   stats::sd(iv), R_SD,        0.08)
report("init_dev mean", mean(iv),      -R_SD^2 / 2, 0.08)

# ---- 5. sim_mod() hands back the truth it drew ------------------------------
cat("\n[5] sim_mod() returns the deviations that generated the data\n")
fit_s <- Rceattle::fit_mod(data_list = e$make_test_data(), estimateMode = 1,
                           msmMode = 0, fit_control = CTL)
set.seed(42)
plain <- Rceattle::sim_mod(fit_s, simulate = TRUE)
ok_flag("process = FALSE attaches nothing", is.null(attr(plain, "process_sim")))

set.seed(42)
withp <- Rceattle::sim_mod(fit_s, simulate = TRUE, process = "recruitment")
tr <- attr(withp, "process_sim")
ok_flag("process = 'recruitment' attaches truth", !is.null(tr),
        paste(names(tr), collapse = ", "))
ok_flag("truth holds rec_dev and init_dev",
        all(c("rec_dev", "init_dev") %in% names(tr)))
ok_flag("truth is not the fitted deviation",
        !is.null(tr) && max(abs(tr$rec_dev - fit_s$estimated_params$rec_dev)) > 0)
ok_flag("M not drawn, so no log_M1_dev in truth",
        is.null(tr$log_M1_dev))
ok_flag("returned object is still a plain data_list",
        identical(class(plain), class(fit_s$data_list)))

cat("\n")
if (length(fails)) {
  cat("FAIL (", length(fails), "):\n", sep = "")
  for (f in unique(fails)) cat("  - ", f, "\n", sep = "")
  quit(status = 1)
}
cat("PASS: process-error draws, sex sharing, 2D map and 2D-AR1 orientation\n")
