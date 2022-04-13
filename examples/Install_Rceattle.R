# Steps to install Rceattle
# Grant D. Adams adamsgd@uw.edu

################################################
# Install from GitHub - Recommended
################################################
# Install devtools if you don't already have it
install.packages("devtools")
# Install TMB and rtools https://cran.r-project.org/bin/windows/Rtools/
# Instructions can be found here for non-pc: https://github.com/kaskr/adcomp/wiki/Download

install.packages('TMB', type = 'source')
# Try "TMB::runExample(all = TRUE)" to see if TMB works

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper") # install tmb helper

# Install Rceattle
devtools::install_github("grantdadams/Rceattle")


################################################
# Install from source
################################################
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
path_to_file <- paste(file_directory, "/", "Rceattle_1.0.0.0000.tar.gz", sep = "")

# Step 6 - Install Rceattle from source
install.packages(path_to_file, repos = NULL, type="source", dependencies = TRUE)


