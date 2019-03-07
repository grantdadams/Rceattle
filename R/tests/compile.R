library(TMB)
TMBfilename = "ceattle_v01_04"

version <- TMBfilename
cpp_directory <- "inst/executables"
cpp_file <- paste0(cpp_directory, "/", version)

TMB::compile(paste0(cpp_file, ".cpp"))
dyn.load(TMB::dynlib(paste0(cpp_file)))
