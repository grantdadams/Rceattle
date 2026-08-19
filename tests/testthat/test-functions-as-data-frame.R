# Build a minimal stand-in `Rceattle` object so we can test the long-format
# extraction without paying the cost of a full TMB fit. The contract under
# test is purely a reshape from data_list/quantities/sdrep -> long form.
make_fake_rceattle <- function(nspp = 2, nyrs = 4, max_age = 3, max_sex = 2,
                               with_sdrep = TRUE) {
  styr <- 2001L
  projyr <- styr + nyrs - 1L

  nages <- c(3, 2)[seq_len(nspp)]
  nsex  <- c(2, 1)[seq_len(nspp)]

  data_list <- list(
    nspp    = nspp,
    styr    = styr,
    endyr   = projyr,
    projyr  = projyr,
    spnames = paste0("sp", seq_len(nspp)),
    nages   = nages,
    nsex    = nsex,
    minage  = c(1L, 2L)[seq_len(nspp)]
  )

  # Helper: deterministic biomass-like matrix
  fill_sy <- function(scale) {
    out <- matrix(NA_real_, nspp, nyrs)
    for (sp in seq_len(nspp)) {
      out[sp, ] <- scale * sp + seq_len(nyrs)
    }
    out
  }

  # Helper: deterministic age/sex/year array (NA for padded cells so we can
  # confirm the function actually drops them rather than keeping NAs)
  fill_ssay <- function() {
    out <- array(NA_real_, dim = c(nspp, max_sex, max_age, nyrs))
    for (sp in seq_len(nspp))
      for (s in seq_len(nsex[sp]))
        for (a in seq_len(nages[sp]))
          for (y in seq_len(nyrs))
            out[sp, s, a, y] <- 1000 * sp + 100 * s + 10 * a + y
    out
  }

  biomass        <- fill_sy(100)
  ssb            <- fill_sy(50)
  R              <- fill_sy(20)
  N_at_age       <- fill_ssay()

  quantities <- list(
    biomass             = biomass,
    ssb                 = ssb,
    R                   = R,
    biomass_depletion   = biomass / max(biomass),
    ssb_depletion       = ssb / max(ssb),
    F_spp               = fill_sy(0.05),
    N_at_age            = N_at_age,
    M_at_age            = N_at_age * 0  # all-zero, valid array
  )

  sdrep <- NULL
  if (with_sdrep) {
    # column-major flattening of the 2D ADREPORT'd quantities, with constant
    # SD per quantity so we can assert the band is symmetric around `value`.
    val <- c(as.numeric(biomass), as.numeric(ssb), as.numeric(R))
    sd  <- c(rep(1.0, length(biomass)),
             rep(0.5, length(ssb)),
             rep(2.0, length(R)))
    nm  <- c(rep("biomass", length(biomass)),
             rep("ssb",     length(ssb)),
             rep("R",       length(R)))
    names(val) <- nm
    sdrep <- list(value = val, sd = sd)
  }

  fit <- list(
    data_list  = data_list,
    quantities = quantities,
    sdrep      = sdrep
  )
  class(fit) <- "Rceattle"
  fit
}


testthat::test_that("as.data.frame.Rceattle returns expected columns and shape", {
  fit <- make_fake_rceattle()
  df  <- as.data.frame(fit)

  testthat::expect_named(
    df,
    c("year", "species", "sex", "age", "quantity", "value", "se", "lwr", "upr")
  )

  # default `which` = 6 species-by-year quantities, two species, four years
  testthat::expect_equal(nrow(df), 6 * 2 * 4)

  # species/year shapes: sex and age are NA
  testthat::expect_true(all(is.na(df$sex)))
  testthat::expect_true(all(is.na(df$age)))

  # year column matches styr:projyr
  testthat::expect_setequal(df$year, 2001:2004)

  # quantity ordering follows the requested `which`
  testthat::expect_equal(unique(df$quantity),
                         c("biomass", "ssb", "R",
                           "biomass_depletion", "ssb_depletion", "F_spp"))
})


testthat::test_that("CIs come from sdrep for ADREPORT'd quantities", {
  fit <- make_fake_rceattle()
  df  <- as.data.frame(fit, which = c("biomass", "ssb", "biomass_depletion"))

  # A strictly positive series takes its interval on the log scale --
  # exp(log(x) +/- z * sd(x)/x) -- so it is right-skewed rather than symmetric.
  # This is the same construction plot_timeseries() draws; the two agreeing is
  # the point (see .rce_ci_bounds()).
  bio <- df[df$quantity == "biomass", ]
  z95 <- stats::qnorm(0.975)
  testthat::expect_equal(bio$upr, exp(log(bio$value) + z95 * 1.0 / bio$value),
                         tolerance = 1e-10)
  testthat::expect_equal(bio$lwr, exp(log(bio$value) - z95 * 1.0 / bio$value),
                         tolerance = 1e-10)
  # Skewed away from the estimate on the high side, and never at zero.
  testthat::expect_true(all((bio$upr - bio$value) > (bio$value - bio$lwr)))
  testthat::expect_true(all(bio$lwr > 0))

  ssb <- df[df$quantity == "ssb", ]
  testthat::expect_equal(ssb$upr, exp(log(ssb$value) + z95 * 0.5 / ssb$value),
                         tolerance = 1e-10)

  # Quantities without an sdreport entry get NA bands.
  bd <- df[df$quantity == "biomass_depletion", ]
  testthat::expect_true(all(is.na(bd$lwr)))
  testthat::expect_true(all(is.na(bd$upr)))
})


testthat::test_that("ci_level controls the band width", {
  fit <- make_fake_rceattle()
  df80 <- as.data.frame(fit, which = "biomass", ci_level = 0.80)
  df95 <- as.data.frame(fit, which = "biomass", ci_level = 0.95)

  # Bounds are on the log scale, so ci_level enters through the quantile rather
  # than as a fixed offset from the estimate.
  testthat::expect_equal(df80$upr,
                         exp(log(df80$value) + stats::qnorm(0.90) / df80$value),
                         tolerance = 1e-10)
  testthat::expect_equal(df95$upr,
                         exp(log(df95$value) + stats::qnorm(0.975) / df95$value),
                         tolerance = 1e-10)
  # ...and a wider level is still a wider band.
  testthat::expect_true(all(df95$upr > df80$upr))
  testthat::expect_true(all(df95$lwr < df80$lwr))
})


testthat::test_that("4D quantities drop sex/age padding for species with fewer", {
  # sp1 has nsex=2, nages=3; sp2 has nsex=1, nages=2.
  # Per year we expect: sp1 = 2*3 = 6 cells, sp2 = 1*2 = 2 cells, total 8/yr.
  fit <- make_fake_rceattle()
  df  <- as.data.frame(fit, which = "N_at_age")

  testthat::expect_equal(nrow(df), 8 * 4)
  testthat::expect_true(all(!is.na(df$value)))

  sp1 <- df[df$species == "sp1", ]
  sp2 <- df[df$species == "sp2", ]
  testthat::expect_setequal(unique(sp1$sex), c(1L, 2L))
  testthat::expect_setequal(unique(sp2$sex), 1L)

  # `age` column is biological age (array index + minage - 1)
  testthat::expect_setequal(unique(sp1$age), c(1L, 2L, 3L))   # minage=1
  testthat::expect_setequal(unique(sp2$age), c(2L, 3L))       # minage=2

  # Spot-check the value at sp1, sex2, age3 (array idx 3), year idx 4.
  cell <- sp1[sp1$sex == 2 & sp1$age == 3 & sp1$year == 2004, ]
  testthat::expect_equal(cell$value, 1000 * 1 + 100 * 2 + 10 * 3 + 4)
})


testthat::test_that("sdrep absent -> all NA bands, no error", {
  fit <- make_fake_rceattle(with_sdrep = FALSE)
  df  <- as.data.frame(fit, which = c("biomass", "ssb"))
  testthat::expect_true(all(is.na(df$lwr)))
  testthat::expect_true(all(is.na(df$upr)))
})


testthat::test_that("which='all' returns every known quantity present", {
  fit <- make_fake_rceattle()
  df  <- as.data.frame(fit, which = "all")
  # Only the quantities we put on the fake fit should appear.
  testthat::expect_setequal(
    unique(df$quantity),
    c("biomass", "ssb", "R", "biomass_depletion", "ssb_depletion",
      "F_spp", "N_at_age", "M_at_age")
  )
})


testthat::test_that("unknown quantity name errors loudly", {
  fit <- make_fake_rceattle()
  testthat::expect_error(as.data.frame(fit, which = "not_a_thing"),
                         "Unknown quantity")
})


testthat::test_that("invalid ci_level errors", {
  fit <- make_fake_rceattle()
  testthat::expect_error(as.data.frame(fit, ci_level = 1.5))
  testthat::expect_error(as.data.frame(fit, ci_level = c(0.9, 0.95)))
})
