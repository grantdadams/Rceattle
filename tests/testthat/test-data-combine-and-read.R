# Regression tests for defects in the data-assembly helpers.
#
# All three were silent: they produced a wrong object rather than an error,
# so nothing downstream complained until much later (or never).

testthat::test_that("combine_data() keeps merged env_data and adds no junk elements", {
  d1 <- make_test_data(seed = 1)
  d2 <- make_test_data(seed = 2)

  # Give the two data sets distinct, partially overlapping env series so a
  # successful merge is observable.
  d1$env_data <- data.frame(Year = 1:8,  temp = seq(0, 1, length.out = 8))
  d2$env_data <- data.frame(Year = 5:12, pdo  = seq(-1, 1, length.out = 8))

  out <- Rceattle::combine_data(d1, d2)

  # The merge used to be assigned to the *input* list and then discarded.
  testthat::expect_false(is.null(out$env_data))
  testthat::expect_true(all(c("Year", "temp", "pdo") %in% names(out$env_data)))
  testthat::expect_setequal(out$env_data$Year, 1:12)

  # The matrix loop used to run over `mat_names[1:17]` for a 15-element
  # vector; `l[[NA_character_]] <- v` succeeds silently, so two NA-named
  # elements were appended to the returned list.
  testthat::expect_false(anyNA(names(out)))
  testthat::expect_true(all(nzchar(names(out))))
})


testthat::test_that("read_data() defaults a missing Month but preserves a supplied one", {
  testthat::skip_if_not_installed("writexl")
  testthat::skip_if_not_installed("readxl")

  dat <- make_test_data(seed = 1)
  dat$fleet_control$Month <- 6  # user-supplied survey month

  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f), add = TRUE)
  Rceattle::write_data(dat, file = f)
  back <- suppressMessages(Rceattle::read_data(file = f))

  # The guard was inverted relative to its own message, so a supplied Month
  # column was overwritten with 0.
  testthat::expect_true(all(back$fleet_control$Month == 6))
})


testthat::test_that("set_phases() has no duplicated parameter names", {
  # `log_M1 = 4` was declared twice. TMBphase() pairs phases to parameters
  # positionally across two differently-ordered lists, so a duplicate
  # silently misaligns every parameter after it.
  ph <- Rceattle:::set_phases()
  testthat::expect_equal(anyDuplicated(names(ph)), 0L)
})


testthat::test_that("build_hcr_map() disables F only when ALL fisheries are inactive", {
  # Build the map once from a valid data set, then vary only the
  # Proj_F_proportion column so the any()-vs-all() branch is isolated.
  dat <- make_test_data(seed = 1)
  dat <- utils::modifyList(dat, Rceattle::build_growth())  # growth switches
  dat <- utils::modifyList(dat, Rceattle::build_hcr(        # HCR switches
    HCR = "NPFMC", DynamicHCR = FALSE,
    Ptarget = 0.4, Plimit = 0.2, Alpha = 0.05))
  dat <- Rceattle::switch_check(dat)
  pars <- suppressWarnings(Rceattle::build_params(data_list = dat))
  map  <- suppressMessages(
    Rceattle::build_map(data_list = dat, params = pars))

  fsh <- which(dat$fleet_control$Fleet_type == "Fishery")
  testthat::skip_if(length(fsh) < 1)

  # Add a second fishery for the same species so "one of several inactive"
  # is representable. build_hcr_map() reads only Species, Fleet_type and
  # Proj_F_proportion, so duplicating the row is sufficient here.
  dat_two <- dat
  dat_two$fleet_control <- rbind(dat_two$fleet_control,
                                 dat_two$fleet_control[fsh[1], ])
  fsh2 <- which(dat_two$fleet_control$Fleet_type == "Fishery")

  # One active, one inactive -> F must stay estimable.
  dat_two$fleet_control$Proj_F_proportion[fsh2] <- c(1, 0)
  mixed <- suppressMessages(Rceattle::build_hcr_map(dat_two, map))
  testthat::expect_false(all(is.na(mixed$mapList$log_Ftarget)))
  testthat::expect_false(all(is.na(mixed$mapList$log_Flimit)))

  # All inactive -> F is correctly turned off (the intended behaviour).
  dat_two$fleet_control$Proj_F_proportion[fsh2] <- 0
  none <- suppressMessages(Rceattle::build_hcr_map(dat_two, map))
  testthat::expect_true(all(is.na(none$mapList$log_Ftarget)))
  testthat::expect_true(all(is.na(none$mapList$log_Flimit)))
})
