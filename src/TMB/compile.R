tmb_name <- "ceattle"
tmb_flags <- commandArgs(trailingOnly = TRUE)
if(length(tmb_flags) == 0) tmb_flags <- ""

# NOTE: The spurious GCC/Eigen -Warray-bounds false positives (from Eigen's
# vectorized reductions inlining TMB's multi-dimensional array indexing, e.g.
# the 5-D suitability(...) access in predation.hpp) are suppressed at the SOURCE
# level via a `#pragma GCC diagnostic ignored "-Warray-bounds"` near the top of
# ceattle.cpp. We deliberately do NOT add a -Wno-array-bounds compiler
# flag here: CRAN's "checking compilation flags used" rejects flags that
# suppress warnings, whereas a targeted source pragma is permitted.
tmb_flags <- paste(tmb_flags, collapse = " ")

# Under covr (covr::package_coverage()), the model would otherwise be built with
# covr's coverage flags, which force -O0 -- making every fit_mod() optimization ~10x
# slower than the -O2 production build (see dev/PERF-findings.md). covr injects those
# flags through a temporary `R_MAKEVARS_USER` (`CXXFLAGS += -O0 --coverage`) whose
# -O0 lands in CXXFLAGS, after any PKG_CXXFLAGS, so it cannot be overridden on the
# compile line. covr reports R-line coverage, NOT C++ coverage, so instrumenting this
# TMB unit buys nothing -- so drop covr's Makevars for THIS compile, letting the
# model build at the normal -O2. The coverage run's fits then run at production speed
# and R-line coverage is unaffected. This file runs in a dedicated Rscript subprocess
# (see ../Makevars), so unsetting the variable here scopes to this compile and its
# child R CMD SHLIB only -- the parent covr process is untouched, no restore needed.
.mkv <- Sys.getenv("R_MAKEVARS_USER")
if (nzchar(.mkv) && file.exists(.mkv) &&
    any(grepl("--coverage|-fprofile", readLines(.mkv, warn = FALSE)))) {
  Sys.unsetenv("R_MAKEVARS_USER")
}

if(file.exists(paste0(tmb_name, ".cpp"))) {
  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               framework = "TMBad",
               safebounds = FALSE, safeunload = FALSE)
  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}

# cleanup done in ../Makevars[.win]
