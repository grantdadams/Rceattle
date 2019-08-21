# Steps to install Rceattle from the source file
# Grant D. Adams adamsgd@uw.edu

# Step 1 - install compilers
# Windows install rtools: https://cran.r-project.org/bin/windows/Rtools/
# Mac install clang and gfortran https://cran.r-project.org/bin/macosx/tools/
# Macs will also need to update makevars file

# Step 2 - Install TMB
# Instructions can be found here for non-pc: https://github.com/kaskr/adcomp/wiki/Download
install.packages('TMB', type = 'source')
# Try "TMB::runExample(all = TRUE)" to see if TMB works

# Step 3 - Install dependencies
install.packages(c("TMB", "testthat", "reshape2", "oce", "plyr", "readxl", "writexl", "tidyr"))

# Step 4 - Download Rceattle_0.0.0.9000.tar.gz file from:
# https://drive.google.com/drive/folders/1sHe_KxvKZi7UyjnWB4Oz9HmKY6CWGWoM

# Step 5 - set "file_directory" to the directory where the download is
file_directory <- "your_download_directory"
path_to_file <- paste(file_directory, "/", "Rceattle_0.0.0.9000.tar.gz", sep = "")

# Step 6 - Install Rceattle from source
install.packages(path_to_file, repos = NULL, type="source", dependencies = TRUE)
