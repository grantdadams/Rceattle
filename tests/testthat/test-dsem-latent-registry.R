# The DSEM latent-process registry is the single place that decides what the
# SEM's latent state variables are called. Six sites depend on it: the default
# sem, the tsdata column injection, the `fixed`-family regex, the recruitment-SD
# self-loop lookup, the rec_dev_col map handed to the C++, and check_dsem_spec().
# These tests pin the contract so that adding growth/M stays a one-line edit and
# so a future change cannot silently desynchronize those six.

testthat::test_that("the registry names one latent series per species, in model order", {
  d <- list(nspp = 3L)

  reg <- .dsem_latent_registry(d)
  testthat::expect_named(reg, "recruitment")
  testthat::expect_equal(reg$recruitment$prefix, "recdevs")
  testthat::expect_equal(reg$recruitment$n, 3L)
  testthat::expect_equal(reg$recruitment$columns,
                         c("recdevs1", "recdevs2", "recdevs3"))
  testthat::expect_equal(.dsem_latent_columns(d), reg$recruitment$columns)

  # Order is positional -- rec_dev_col indexes x_tj by it -- so a double-digit
  # species count must not sort lexically.
  testthat::expect_equal(.dsem_latent_columns(list(nspp = 11L))[10:11],
                         c("recdevs10", "recdevs11"))
})

testthat::test_that("the latent pattern matches exactly the latent column names", {
  p <- .dsem_latent_pattern()

  testthat::expect_true(all(grepl(p, c("recdevs1", "recdevs10", "recdevs001"))))
  # Near-misses must NOT match: these are ordinary env-data variables.
  testthat::expect_false(any(grepl(p, c("recdevs", "recdevsA", "recdevs_anom",
                                        "recdevs1_lag", "xrecdevs1", "BT"))))
})

testthat::test_that("the pattern needs no data_list, so check_dsem_spec keeps its contract", {
  # Regression: folding the name stems into the per-model registry made
  # .dsem_latent_pattern() require data_list$nspp, and check_dsem_spec() -- whose
  # documented contract asks only for env_data/styr/endyr -- began erroring on
  # inputs it used to accept.
  testthat::expect_silent(.dsem_latent_pattern())
  testthat::expect_equal(length(formals(.dsem_latent_pattern)), 0L)

  sf <- data.frame(path = "BT -> recdevs1", lag = 1, name = "b", start = 0,
                   parameter = 1L, first = "BT", second = "recdevs1",
                   direction = 1, stringsAsFactors = FALSE)
  d_no_nspp <- list(env_data = data.frame(Year = 1990:2000, BT = 0),
                    styr = 1990, endyr = 2000)
  out <- check_dsem_spec(d_no_nspp, list(sem_full = sf))
  testthat::expect_s3_class(out, "Rceattle_convergence")
})

testthat::test_that("latent columns are injected without clobbering an env_data column", {
  testthat::skip_if_not_installed("dsem")

  # Regression: the injection used a fixed temporary column name, so an env_data
  # column of that name was silently overwritten.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L, spnames = "a",
            sigma_rec_prior = 1, random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, .dsem_latent = 1, BT = 0,
                                  check.names = FALSE))

  obj <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = "BT -> recdevs1, 1, b, 0
                                      recdevs1 <-> recdevs1, 0, sigmaR1, 1"),
    data_list = d)

  # The latent column leads, the referenced env variable survives, and the
  # unreferenced .dsem_latent column is dropped as usual (not silently reused).
  testthat::expect_equal(colnames(obj$tmb_inputs$data$y_tj)[1], "recdevs1")
  testthat::expect_true("BT" %in% colnames(obj$tmb_inputs$data$y_tj))
  testthat::expect_equal(obj$tmb_inputs$data$rec_dev_col, 0L)
})
