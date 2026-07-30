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

if(file.exists(paste0(tmb_name, ".cpp"))) {
  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               framework = "TMBad",
               safebounds = FALSE, safeunload = FALSE)
  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}

# cleanup done in ../Makevars[.win]
