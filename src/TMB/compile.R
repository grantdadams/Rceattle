tmb_name <- "ceattle_v01_11"
tmb_flags <- commandArgs(trailingOnly = TRUE)
if(length(tmb_flags) == 0) tmb_flags <- ""

# Silence the spurious -Warray-bounds diagnostics ("array subscript ... partly
# outside array bounds", "return *__P", "... allocated by malloc") that GCC
# 12-14 emit from Eigen's vectorized reductions when they inline TMB's multi-
# dimensional array indexing (e.g. the 5-D suitability(...) access in
# predation.hpp). The index math is correct -- this is a well-known GCC/Eigen
# false positive, not a bug in this package. We disable ONLY -Warray-bounds, so
# real diagnostics (e.g. -Wmaybe-uninitialized) still surface. Harmless on
# compilers that do not emit it (Clang accepts and ignores -Wno-array-bounds).
tmb_flags <- paste("-Wno-array-bounds", paste(tmb_flags, collapse = " "))

if(file.exists(paste0(tmb_name, ".cpp"))) {
  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               framework = "TMBad",
               safebounds = FALSE, safeunload = FALSE)
  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}

# cleanup done in ../Makevars[.win]
