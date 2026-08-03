# The diet Dirichlet-multinomial weight (diet_comp_weights / theta_diet) is only
# informed where the diet composition likelihood is fit: the template enters that
# block for `msmMode > 2 || max(suitMode) > 0` and skips any predator with
# `suitMode <= 0`. These tests pin (a) that build_map() frees the weight on
# exactly that set, and (b) that a theta_diet prior aimed at a fixed weight is
# ignored rather than adding a constant to the objective.

msm_dm_data <- function(nyrs = 30) {
  sim <- make_msm_test_data(years = seq_len(nyrs))
  d <- sim$data_list
  d$Diet_distribution <- rep(1L, d$nspp)   # DirichletMultinomial diet comps
  d
}

# make_msm_test_data() is a parametric-growth fixture (CAAL populated), so the
# fit has to be told to use von Bertalanffy growth.
build_only <- function(data_list, msmMode, suitMode, ...) {
  suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = data_list,
                      estimateMode = 3,          # build, do not optimize
                      msmMode = msmMode,
                      suitMode = suitMode,
                      niter = 3,
                      growthFun = build_growth(fun = "vonBertalanffy"),
                      fit_control = fit_control(getsd = FALSE, verbose = 0),
                      ...)))
}

diet_map <- function(data_list, msmMode, suitMode) {
  build_only(data_list, msmMode, suitMode)$map$mapList$diet_comp_weights
}


test_that("the diet DM weight is estimated only where the diet likelihood is fit", {
  d <- msm_dm_data()

  # Single species: no predation, no diet likelihood.
  expect_true(all(is.na(diet_map(d, msmMode = 0, suitMode = 0))))

  # Multispecies with empirical suitability: diet proportions are taken as
  # given, so there is still nothing to estimate the overdispersion from.
  expect_true(all(is.na(diet_map(d, msmMode = 1, suitMode = 0))))

  # Multispecies with a parametric (weight-based lognormal) suitability: the
  # diet composition IS fit, so the weight is free.
  expect_true(all(!is.na(diet_map(d, msmMode = 1, suitMode = 4))))
})


test_that("the diet DM weight follows suitMode per predator", {
  d <- msm_dm_data()
  m <- diet_map(d, msmMode = 1, suitMode = c(0, 4))
  expect_true(is.na(m[1]))
  expect_false(is.na(m[2]))
})


test_that("a multinomial diet keeps its weight fixed even under parametric suitability", {
  d <- msm_dm_data()
  d$Diet_distribution <- rep(0L, d$nspp)
  expect_true(all(is.na(diet_map(d, msmMode = 1, suitMode = 4))))
})


test_that("a theta_diet prior on a fixed weight is ignored, not added as a constant", {
  d <- msm_dm_data()

  compFun <- build_composition(linkages = list(
    theta_diet = linkage_spec(~ 1, by = ~ species,
                              species = seq_len(d$nspp),
                              priors = list(`(Intercept)` = prior_lognormal(0, 2)))))

  # msmMode = 0: no predation, so the diet weight is fixed.
  build <- function(cf) build_only(d, msmMode = 0, suitMode = 0, compFun = cf)

  with_prior <- build(compFun)
  no_prior   <- build(build_composition())

  # The prior is neutralized, so it contributes nothing anywhere.
  expect_equal(sum(with_prior$quantities$jnll_comp["Linkage-table priors", ]), 0)
  expect_equal(with_prior$quantities$jnll, no_prior$quantities$jnll)

  # The rows survive (beta_linkage keeps its length, so `inits` stay portable
  # between a single-species and a multispecies fit sharing one compFun) but
  # carry no prior family.
  tab <- with_prior$data_list$linkage_table
  expect_equal(nrow(tab), d$nspp)
  expect_true(all(tab$prior_family == "none"))
})


test_that("a theta_diet prior is kept where the weight is estimated", {
  d <- msm_dm_data()
  d$suitMode <- rep(4, d$nspp)

  compFun <- build_composition(linkages = list(
    theta_diet = linkage_spec(~ 1, by = ~ species,
                              species = seq_len(d$nspp),
                              priors = list(`(Intercept)` = prior_lognormal(0, 2)))))

  fit <- build_only(d, msmMode = 1, suitMode = 4, compFun = compFun)

  expect_true(all(fit$data_list$linkage_table$prior_family == "lognormal"))
  expect_gt(sum(fit$quantities$jnll_comp["Linkage-table priors", ]), 0)
})
