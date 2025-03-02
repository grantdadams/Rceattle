tmb_name <- "ceattle_v01_10"
tmb_flags <- commandArgs(trailingOnly = TRUE)

if(file.exists(paste0(tmb_name, ".cpp"))) {
  if(length(tmb_flags) == 0) tmb_flags <- ""
  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               framework = "TMBad",
               safebounds = FALSE, safeunload = FALSE)
  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}

# cleanup done in ../Makevars[.win]
