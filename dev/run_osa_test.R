# Run the Phase 1 OSA testthat file directly (it lives in a subdir that
# devtools::test(filter=) does not scan). Usage: Rscript R/dev/run_osa_test.R
suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))
library(testthat)
e <- new.env(parent = environment())
source("tests/testthat/helpers.R", local = e)
source("tests/testthat/helpers-make-msm-data.R", local = e)
testthat::test_file("tests/testthat/tests-Likelihoods/test-osa-residuals.R", env = e)
