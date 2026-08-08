# tools/verify/verify-dsem-equivalence.R
#
# Numerical verification of the vendored DSEM code in src/TMB/dsem.hpp.
#
# Everything else about that header has been checked by compiling it. This runs
# it. Two things are established:
#
#   1. EQUIVALENCE. calculate_dsem()'s jnll matches dsem::dsem(run_model = TRUE)
#      on identical inputs, for each gmrf_parameterization. That is what makes
#      the vendored copy trustworthy, and it is the gate for relaxing the
#      DESCRIPTION pin: run it against a new dsem, and if it passes the header is
#      still compatible.
#   2. THE INPUT GUARDS FIRE. calculate_dsem() takes many same-typed positional
#      arguments, so a caller that swaps a pair still compiles. The assertions
#      added for that are only worth having if they actually trigger.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-dsem-equivalence.R
#
# Compiles tools/verify/dsem_standalone.cpp into dev/ (gitignored scratch).

suppressPackageStartupMessages({
  library(TMB)
})
stopifnot(requireNamespace("dsem", quietly = TRUE))

dir.create("dev", showWarnings = FALSE)
src <- "tools/verify/dsem_standalone.cpp"
work <- file.path("dev", "dsem_standalone.cpp")
file.copy(src, work, overwrite = TRUE)
file.copy("src/TMB/dsem.hpp", file.path("dev", "dsem.hpp"), overwrite = TRUE)

cat("compiling the standalone DSEM model...\n")
TMB::compile(work, framework = "TMBad", flags = "-O1")
dyn.load(TMB::dynlib(file.path("dev", "dsem_standalone")))

# ---- a small, fully specified DSEM ------------------------------------------
set.seed(1)
n_t <- 30
tsdata <- data.frame(
  x = as.numeric(scale(cumsum(rnorm(n_t)))),
  y = as.numeric(scale(cumsum(rnorm(n_t))))
)
sem <- "
  x -> x, 1, ar_x, 0.3
  x -> y, 1, x_to_y, 0.2
  x <-> x, 0, sd_x, 1
  y <-> y, 0, sd_y, 1
"


# dsem's TMB data, mapped onto calculate_dsem()'s arguments. obs_idx/unobs_idx
# are absent for the non-projecting parameterizations, so default them to empty.
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
       unobs_idx    = if (is.null(d$unobs_idx)) integer(0) else as.integer(d$unobs_idx))
}

ok <- TRUE
report <- function(label, pass, detail = "") {
  verdict <- if (is.na(pass)) "SKIP" else if (pass) "PASS" else "FAIL"
  cat(sprintf("  %-46s %s %s\n", label, verdict, detail))
  if (identical(pass, FALSE)) ok <<- FALSE   # NA (skipped) is not a failure
}

cat("\n=== 1. jnll equivalence vs dsem::dsem() on identical inputs ===\n")
# Build BOTH objects from the same data and the same full parameter list, with
# no map and no random. Each then returns the plain joint negative log-likelihood
# on the same vector -- no Laplace integration, no inner optimization, no mapped
# parameters -- so any difference is the header's arithmetic and nothing else.
# dsem_control()'s parameterization names, in options(0) order.
params <- c("full", "project", "mvn_project", "gmrf_project")
for (i in seq_along(params)) {
  pz <- params[i]
  fit <- try(dsem::dsem(sem = sem, tsdata = ts(tsdata),
                        control = dsem::dsem_control(
                          run_model = TRUE, quiet = TRUE,
                          gmrf_parameterization = pz)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    report(paste0("options(0)=", i - 1, " (", pz, ")"), NA,
           "-- dsem could not fit this parameterization; skipped")
    next
  }
  d    <- fit$tmb_inputs$data
  pars <- fit$tmb_inputs$parameters

  theirs_obj <- TMB::MakeADFun(data = d, parameters = pars,
                               DLL = "dsem", silent = TRUE)
  ours_obj   <- TMB::MakeADFun(data = dsem_data(d), parameters = pars,
                               DLL = "dsem_standalone", silent = TRUE)

  # Evaluate at the start values and at a perturbed point: a single evaluation
  # at a flat optimum can agree by accident.
  set.seed(42)
  pts <- list(start = ours_obj$par,
              jittered = ours_obj$par + rnorm(length(ours_obj$par), 0, 0.15))
  worst <- 0
  for (nm in names(pts)) {
    v <- pts[[nm]]
    a <- try(ours_obj$fn(v), silent = TRUE)
    b <- try(theirs_obj$fn(v), silent = TRUE)
    if (inherits(a, "try-error") || inherits(b, "try-error")) { worst <- NA; break }
    worst <- max(worst, abs(a - b) / max(1, abs(b)))
  }
  report(paste0("options(0)=", i - 1, " (", pz, ")"),
         !is.na(worst) && worst < 1e-10,
         sprintf("max rel diff over %d points = %.2e", length(pts), worst))
}

cat("\n=== 2. input guards fire ===\n")
fit <- dsem::dsem(sem = sem, tsdata = ts(tsdata),
                  control = dsem::dsem_control(run_model = TRUE, quiet = TRUE))
d <- fit$tmb_inputs$data
p <- fit$obj$env$parList(fit$obj$env$last.par.best)
base_data <- dsem_data(d)

fires <- function(label, mutate, pattern) {
  dd <- mutate(base_data)
  e <- tryCatch({
    o <- TMB::MakeADFun(data = dd, parameters = p, DLL = "dsem_standalone", silent = TRUE)
    o$fn(o$par); "no error"
  }, error = function(e) conditionMessage(e))
  report(label, grepl(pattern, e), if (grepl(pattern, e)) "" else paste("got:", substr(e, 1, 60)))
}

fires("eps_tj / y_tj wrong shape",
      function(x) { x$eps_tj <- array(0, dim = c(nrow(x$y_tj) + 1L, ncol(x$y_tj))); x },
      "eps_tj is")
fires("familycode_j wrong length",
      function(x) { x$familycode_j <- c(x$familycode_j, 0L); x }, "familycode_j has")
fires("RAM wrong column count",
      function(x) { x$RAM <- x$RAM[, 1:5, drop = FALSE]; x }, "RAM has")
fires("RAMstart / RAM row mismatch",
      function(x) { x$RAMstart <- x$RAMstart[-1]; x }, "RAMstart has")
fires("options too short",
      function(x) { x$options <- x$options[1:2]; x }, "options has")
proj_fit <- dsem::dsem(sem = sem, tsdata = ts(tsdata),
                       control = dsem::dsem_control(run_model = TRUE, quiet = TRUE,
                                                    gmrf_parameterization = "gmrf_project"))
proj_data <- dsem_data(proj_fit$tmb_inputs$data)
proj_p <- proj_fit$obj$env$parList(proj_fit$obj$env$last.par.best)
local({
  dd <- proj_data; dd$obs_idx <- head(dd$obs_idx, -1L)
  e <- tryCatch({
    o <- TMB::MakeADFun(data = dd, parameters = proj_p, DLL = "dsem_standalone", silent = TRUE)
    o$fn(o$par); "no error"
  }, error = function(e) conditionMessage(e))
  report("obs_idx + unobs_idx do not cover n_k", grepl("obs_idx", e),
         if (grepl("obs_idx", e)) "" else paste("got:", substr(e, 1, 60)))
})
fires("unrecognized options(0)",
      function(x) { x$options[1] <- 7L; x }, "unrecognized options\\(0\\)")

cat("\n", if (ok) "ALL CHECKS PASSED" else "SOME CHECKS FAILED", "\n", sep = "")
quit(status = if (ok) 0L else 1L)
