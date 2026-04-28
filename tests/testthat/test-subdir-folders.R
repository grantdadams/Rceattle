
testthat::skip_on_cran()
testthat::skip_if(nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                  "Subdirectory test_dir() runs require manual sourcing of helpers; run interactively.")
source(testthat::test_path("helpers-make-msm-data.R")) # for cran checks
source(testthat::test_path("helpers.R"))
testthat::test_dir("tests-Dynamics/")
testthat::test_dir("tests-Mortality/")
testthat::test_dir("tests-Likelihoods/")
testthat::test_dir("tests-Data-processing/")
testthat::test_dir("tests-Selectivity/")
testthat::test_dir("tests-Functions/")
