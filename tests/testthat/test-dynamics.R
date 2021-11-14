test_that("single-species population dynamics is correct", {
  data(BS2017SS) # ?BS2017SS for more information on the data
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              debug = 1, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              verbose = 0)

  expect_equal(round(ss_run$quantities$biomassSSB[1,10],4), 164.4717)
  expect_equal(round(ss_run$quantities$biomassSSB[2,10],4), 553.5195)
  expect_equal(round(ss_run$quantities$biomassSSB[3,10],4), 578.1277)
})
#> Test passed 😀


test_that("multi-species population dynamics is correct", {
  data(BS2017SS) # ?BS2017SS for more information on the data
  ms_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              debug = 1, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 1, # Single species mode
                              niter = 3,
                              verbose = 0)

  expect_equal(round(ms_run$quantities$biomassSSB[1,2],4), 84.7607)
  expect_equal(round(ms_run$quantities$biomassSSB[2,2],4), 6434.2652)
  expect_equal(round(ms_run$quantities$biomassSSB[3,2],4), 11153.5178)
})
#> Test passed 😀
