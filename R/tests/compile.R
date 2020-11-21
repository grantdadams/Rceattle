library(TMB)
TMBfilename = "ceattle_v01_06"

version <- TMBfilename
cpp_directory <- "inst/executables"
cpp_file <- paste0(cpp_directory, "/", version)
dyn.unload(TMB::dynlib(paste0(cpp_file)))
TMB::compile(paste0(cpp_file, ".cpp"), CPPFLAGS="-Wno-ignored-attributes")
dyn.load(TMB::dynlib(paste0(cpp_file)))


