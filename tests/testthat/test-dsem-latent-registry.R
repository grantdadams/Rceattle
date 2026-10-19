# The DSEM latent-process registry is the single place that decides what the
# SEM's latent state variables are called. Six sites depend on it: the default
# sem, the tsdata column injection, the `fixed`-family regex, the recruitment-SD
# self-loop lookup, the rec_dev_col / growth_dev_* maps handed to the C++, and
# check_dsem_spec(). These tests pin the contract so a future change cannot
# silently desynchronize those six.
#
# One entry per drivable TARGET, not per process: growth contributes three
# (`L1`, `Linf`, Richards `m`). Growth `K` is deliberately NOT drivable: the
# Jan-1 recursion reads the PREVIOUS year's K, so the terminal K latent reaches
# a hindcast observation only through the within-year advance, whose factor
# `exp(-kappa * fracyr) - 1` is exactly zero at a Jan-1 fleet/spawning month.
# Measured there: bit-identical objective under perturbation. Not structural --
# at a mid-year month it is informed, though weakly.
# Natural mortality keeps `mdevs` for when it lands, which is why Richards `m`
# is `richdevs`; and a trailing `m` has to stay readable as the sex token.

testthat::test_that("the registry names one latent series per target, in model order", {
  d <- list(nspp = 3L)

  reg <- .dsem_latent_registry(d)
  testthat::expect_named(reg, c("recruitment", "growth_L1", "growth_Linf",
                                "growth_m"))
  testthat::expect_equal(reg$recruitment$prefix, "recdevs")
  testthat::expect_equal(reg$recruitment$n, 3L)
  testthat::expect_equal(reg$recruitment$columns,
                         c("recdevs1", "recdevs2", "recdevs3"))
  testthat::expect_equal(reg$growth_L1$columns,
                         c("l1devs1", "l1devs2", "l1devs3"))
  testthat::expect_equal(reg$growth_Linf$columns,
                         c("linfdevs1", "linfdevs2", "linfdevs3"))
  testthat::expect_equal(reg$growth_m$columns,
                         c("richdevs1", "richdevs2", "richdevs3"))
  # Every declared column, recruitment first then growth, in registry order.
  testthat::expect_equal(.dsem_latent_columns(d),
                         c(reg$recruitment$columns, reg$growth_L1$columns,
                           reg$growth_Linf$columns, reg$growth_m$columns))

  # Order is positional -- rec_dev_col indexes x_tj by it -- so a double-digit
  # species count must not sort lexically.
  testthat::expect_equal(.dsem_latent_columns(list(nspp = 11L))[10:11],
                         c("recdevs10", "recdevs11"))
})

testthat::test_that("growth targets carry the growth_linkage_offset slot, from one source", {
  # `param` indexes growth_linkage_offset's last dimension. It is read from
  # LINKAGE_PARAM_CODES rather than written out again, so the sem names and the
  # linkage grammar cannot disagree about which slot `Linf` is.
  reg <- .dsem_latent_registry(list(nspp = 1L))
  testthat::expect_equal(reg$growth_L1$param,   LINKAGE_PARAM_CODES$growth[["L1"]])
  testthat::expect_equal(reg$growth_Linf$param, LINKAGE_PARAM_CODES$growth[["Linf"]])
  testthat::expect_equal(reg$growth_m$param,    LINKAGE_PARAM_CODES$growth[["m"]])
  # Recruitment drives no growth slot.
  testthat::expect_true(all(is.na(reg$recruitment$param)))
  # K is not drivable, so no entry may claim its slot.
  testthat::expect_false(LINKAGE_PARAM_CODES$growth[["K"]] %in%
                           unlist(lapply(reg, `[[`, "param")))
})

testthat::test_that("growth latents are per sex, and a one-sex species carries no token", {
  # The package convention everywhere else: a one-sex species has no sex token
  # (parameter_index() reports NA, data_check() uses ""). Growth parameters are
  # sex-specific, so a two-sex species needs one latent per sex. The f/m token
  # is the spelling linkage_spec() already accepts on input; the package's own
  # OUTPUT labels are female/male, so this is an input alias, not a rename.
  reg <- .dsem_latent_registry(list(nspp = 2L, nsex = c(1L, 2L)))
  testthat::expect_equal(reg$growth_L1$columns,
                         c("l1devs1", "l1devs2f", "l1devs2m"))
  testthat::expect_equal(reg$growth_L1$species, c(1L, 2L, 2L))
  testthat::expect_equal(reg$growth_L1$sex,     c(1L, 1L, 2L))
  # Recruitment is not sex-specific whatever nsex says.
  testthat::expect_equal(reg$recruitment$columns, c("recdevs1", "recdevs2"))

  # nsex absent on a partially built data_list reads as one sex; data_check()
  # is what enforces its presence.
  testthat::expect_equal(.dsem_latent_registry(list(nspp = 1L))$growth_L1$columns,
                         "l1devs1")
})

testthat::test_that("the latent pattern matches exactly the latent column names", {
  p <- .dsem_latent_pattern()

  testthat::expect_true(all(grepl(p, c("recdevs1", "recdevs10", "recdevs001",
                                       "l1devs1", "linfdevs2", "richdevs1",
                                       "l1devs2f", "l1devs2m", "linfdevs10m"))))
  # Near-misses must NOT match: these are ordinary env-data variables.
  testthat::expect_false(any(grepl(p, c("recdevs", "recdevsA", "recdevs_anom",
                                        "recdevs1_lag", "xrecdevs1", "BT",
                                        "l1devs1x", "l1devs1fm", "l1devsf"))))
  # `mdevs` belongs to natural mortality and `kdevs` to a parameter that is not
  # drivable; neither may be recognised as a latent today, or an env column by
  # that name would be silently shadowed.
  testthat::expect_false(any(grepl(p, c("mdevs1", "mdevs2f", "kdevs1"))))
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

testthat::test_that("a recruitment-only sem is untouched by the growth targets", {
  testthat::skip_if_not_installed("dsem")

  # The registry DECLARES growth latents for every model, so the thing that
  # keeps existing fits on the existing code path is the column drop: a latent
  # column the sem does not name never reaches the tsdata. This is the FRESH
  # build path; the saved-object refit path is .dsem_tmb_inputs(), pinned
  # separately below.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L, spnames = "a",
            sigma_rec_prior = 1, random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, BT = 0))

  obj <- build_dsem_objects(
    dsem_settings = build_DSEM(),          # default: IID recruitment only
    data_list = d)

  cols <- colnames(obj$tmb_inputs$data$y_tj)
  testthat::expect_equal(cols[grepl(.dsem_latent_pattern(), cols)], "recdevs1")
  # The growth index map is empty, which is the branch the C++ reads as "no
  # growth latent".
  testthat::expect_length(obj$tmb_inputs$data$growth_dev_col, 0L)
  testthat::expect_equal(obj$tmb_inputs$data$drives_rec, 1L)
})

testthat::test_that("a sem may drive growth without driving recruitment", {
  testthat::skip_if_not_installed("dsem")

  # Naming NO recruitment column is a growth-only SEM, which must build with
  # drives_rec = 0 so the template keeps scoring the recruitment deviations on
  # their own density.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L, spnames = "a",
            sigma_rec_prior = 1, random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, BT = 0))

  obj <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = "l1devs1 <-> l1devs1, 0, sigL1, 0.3"),
    data_list = d)

  testthat::expect_equal(obj$tmb_inputs$data$drives_rec, 0L)
  # The sentinel is -1, NOT 0. Several R consumers index x_tj as
  # rec_dev_col[sp] + 1L, and 0 would silently alias the first latent column --
  # a growth field -- where -1 gives index 0 and errors. NA is not an option
  # either: it reaches a DATA_IVECTOR as INT_MIN.
  testthat::expect_equal(obj$tmb_inputs$data$rec_dev_col, -1L)
  testthat::expect_false(anyNA(obj$tmb_inputs$data$rec_dev_col))
  testthat::expect_equal(obj$tmb_inputs$data$growth_dev_par,
                         LINKAGE_PARAM_CODES$growth[["L1"]])
  # The VALUE, not just the length: review removed the `- 1L` that makes this
  # 0-based and every other assertion in this file still passed. l1devs1 is the
  # only latent column in this sem, so it is x_tj column 0.
  testthat::expect_equal(obj$tmb_inputs$data$growth_dev_col, 0L)
  testthat::expect_equal(obj$tmb_inputs$data$growth_dev_sp,  0L)
  testthat::expect_equal(obj$tmb_inputs$data$growth_dev_sex, 0L)
})

testthat::test_that("the growth index map decodes species and sex, 0-based", {
  testthat::skip_if_not_installed("dsem")

  # Two species, the second two-sex, with L1 driven for the MALE of species 2
  # only. The decode is the one thing standing between a sem variable name and
  # which (species, sex, growth parameter) slot the template will write, so the
  # 0-based conversion is pinned on a case where every index is non-zero.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 2L,
            spnames = c("a", "b"), nsex = c(1L, 2L),
            sigma_rec_prior = c(1, 1), random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, BT = 0))

  obj <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = paste0(
      "recdevs1 <-> recdevs1, 0, sigmaR1, 1\n",
      "recdevs2 <-> recdevs2, 0, sigmaR2, 1\n",
      "l1devs2m <-> l1devs2m, 0, sigL1, 0.3")),
    data_list = d)

  dt <- obj$tmb_inputs$data
  testthat::expect_length(dt$growth_dev_col, 1L)
  testthat::expect_equal(dt$growth_dev_sp,  1L)   # species 2, 0-based
  testthat::expect_equal(dt$growth_dev_sex, 1L)   # male, 0-based
  testthat::expect_equal(dt$growth_dev_par, LINKAGE_PARAM_CODES$growth[["L1"]])
  # The column index must point at l1devs2m's own column of the latent matrix.
  testthat::expect_equal(dt$growth_dev_col,
                         which(colnames(dt$y_tj) == "l1devs2m") - 1L)
})

testthat::test_that("a saved pre-growth dsem object refits with no growth latents", {
  # Around 30 saved $dsem objects in the assessment repos refit through
  # .dsem_tmb_inputs(), NOT through build_dsem_objects(), and they carry no
  # growth_dev_* fields at all. The %||% defaults are what keeps them loading;
  # without them the refit errors on a missing field.
  saved <- list(tmb_inputs = list(
    data = list(drives_rec = 1L, nyrs_dsem = 11L, options = integer(4),
                RAM = matrix(0L, 0, 6), RAMstart = numeric(0),
                familycode_j = 0L, linkcode_j = 0L, sigmastart_j = 0L,
                eps_tj = array(0, c(11, 1)), y_tj = array(0, c(11, 1)),
                obs_idx = integer(0), unobs_idx = integer(0),
                rec_dev_col = 0L, rec_sd_idx = 0L, rec_sd_fixed = 0,
                rec_sd_prior = 0, rec_sd_prior_sd = 0, rec_sd_use_prior = 0L),
    parameters = list(beta_z = numeric(0), lnsigma_z = numeric(0),
                      mu_j = numeric(0), delta0_j = numeric(0),
                      x_tj = array(0, c(11, 1)))))

  got <- .dsem_tmb_inputs(saved)$data
  testthat::expect_length(got$growth_dev_col, 0L)
  testthat::expect_length(got$growth_dev_sp,  0L)
  testthat::expect_length(got$growth_dev_sex, 0L)
  testthat::expect_length(got$growth_dev_par, 0L)
  # A SEM built before the flag existed could only have driven recruitment.
  testthat::expect_equal(got$dsem_drives_rec, 1L)
  # The null and real field sets must stay identical, or TMB will be handed a
  # different data list depending on which path built it.
  testthat::expect_setequal(names(got), names(.dsem_null_inputs()$data))
})

testthat::test_that("driving recruitment for only SOME species is refused", {
  testthat::skip_if_not_installed("dsem")

  # A subset is a typo in a species index, not a design: it would delete the
  # recruitment process for the unnamed species while keeping it for the rest.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 2L,
            spnames = c("a", "b"), sigma_rec_prior = c(1, 1), random_rec = TRUE,
            proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, BT = 0))

  testthat::expect_error(
    build_dsem_objects(
      dsem_settings = build_DSEM(sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1"),
      data_list = d),
    "1 of 2 species")
})

testthat::test_that("a sem that drives no process at all is refused", {
  testthat::skip_if_not_installed("dsem")

  # Before growth became drivable, the all-species recruitment refusal caught
  # this. Such a model still fits a GMRF to its covariates and still adds
  # jnll_dsem to the objective while nothing in the population model reads it,
  # so it converges and reports an ABC whose objective is not comparable to the
  # same model without a SEM.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L, spnames = "a",
            sigma_rec_prior = 1, random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, BT = 0))

  testthat::expect_error(
    build_dsem_objects(dsem_settings = build_DSEM(sem = "BT <-> BT, 0, sdBT, 1"),
                       data_list = d),
    "drives nothing")
})

testthat::test_that(".dsem_drives() answers per growth parameter, not just recruitment", {
  testthat::skip_if_not_installed("dsem")

  # Regression: widening .DSEM_LATENT_PREFIXES widened the set of processes
  # .dsem_drives() ACCEPTS, while its switch() still had only a recruitment
  # branch -- so a growth question got a silent FALSE where it used to error.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L, spnames = "a",
            sigma_rec_prior = 1, random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, BT = 0))

  g <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = "l1devs1 <-> l1devs1, 0, sigL1, 0.3"),
    data_list = d)
  testthat::expect_false(.dsem_drives(g, "recruitment"))
  testthat::expect_true(.dsem_drives(g, "growth_L1"))
  testthat::expect_false(.dsem_drives(g, "growth_Linf"))
  testthat::expect_false(.dsem_drives(g, "growth_m"))

  r <- build_dsem_objects(dsem_settings = build_DSEM(), data_list = d)
  testthat::expect_true(.dsem_drives(r, "recruitment"))
  testthat::expect_false(.dsem_drives(r, "growth_L1"))
})

testthat::test_that("a growth latent the species' growth model never reads is refused", {
  testthat::skip_if_not_installed("dsem")

  # The silent-wrong-model case this guards. Empirical weight-at-age estimates
  # no growth parameters, and von Bertalanffy has no Richards `m` (log_m is
  # mapped out at R/3-build_map.R, npar = 3). A latent on either is a free field
  # carrying its own GMRF density and driving nothing, so the model converges
  # and reports an ABC for a process the user believes is time-varying.
  # Measured before the guard: `richdevs1` on a vB species built with
  # growth_dev_par = 3.
  mk <- function(gm, nspp = 1L) {
    list(styr = 1990, endyr = 2000, projyr = 2000, nspp = nspp,
         spnames = letters[seq_len(nspp)], nsex = rep(1L, nspp),
         growth_model = gm, sigma_rec_prior = rep(1, nspp),
         random_rec = TRUE, proj_mean_rec = TRUE,
         env_data = data.frame(Year = 1990:2000, BT = 0))
  }
  R1 <- "recdevs1 <-> recdevs1, 0, sR, 1\n"

  testthat::expect_error(
    build_dsem_objects(dsem_settings = build_DSEM(
      sem = paste0(R1, "l1devs1 <-> l1devs1, 0, sL, .3\n")), data_list = mk(0L)),
    "empirical weight-at-age")
  testthat::expect_error(
    build_dsem_objects(dsem_settings = build_DSEM(
      sem = paste0(R1, "richdevs1 <-> richdevs1, 0, sm, .3\n")),
      data_list = mk(1L)),
    "Richards shape parameter")

  # The legal combinations still build.
  for (s in c("l1devs1 <-> l1devs1, 0, sL, .3\n",
              "linfdevs1 <-> linfdevs1, 0, sI, .3\n")) {
    testthat::expect_s3_class(
      build_dsem_objects(dsem_settings = build_DSEM(sem = paste0(R1, s)),
                         data_list = mk(1L)), "dsem")
  }
  testthat::expect_equal(
    build_dsem_objects(dsem_settings = build_DSEM(
      sem = paste0(R1, "richdevs1 <-> richdevs1, 0, sm, .3\n")),
      data_list = mk(2L))$tmb_inputs$data$growth_dev_par,
    LINKAGE_PARAM_CODES$growth[["m"]])

  # Per species, not model-wide: sp1 vB, sp2 Richards.
  R2 <- paste0(R1, "recdevs2 <-> recdevs2, 0, sR2, 1\n")
  d2 <- mk(c(1L, 2L), nspp = 2L)
  testthat::expect_equal(
    build_dsem_objects(dsem_settings = build_DSEM(
      sem = paste0(R2, "richdevs2 <-> richdevs2, 0, sm, .3\n")),
      data_list = d2)$tmb_inputs$data$growth_dev_sp, 1L)
  testthat::expect_error(
    build_dsem_objects(dsem_settings = build_DSEM(
      sem = paste0(R2, "richdevs1 <-> richdevs1, 0, sm, .3\n")),
      data_list = d2),
    "Richards shape parameter")

  # A direct caller may pass no growth_model; there is nothing to check then.
  d3 <- mk(1L); d3$growth_model <- NULL
  testthat::expect_s3_class(
    build_dsem_objects(dsem_settings = build_DSEM(
      sem = paste0(R1, "l1devs1 <-> l1devs1, 0, sL, .3\n")), data_list = d3),
    "dsem")
})

testthat::test_that("an unrecognised latent-looking sem variable already fails loudly", {
  testthat::skip_if_not_installed("dsem")

  # Pins that NO Rceattle-side refusal is needed for these: dsem itself rejects
  # a sem variable absent from the tsdata. Recorded because one review called
  # them "silently intersect()ed away", which measurement contradicted -- and
  # because if dsem ever stopped erroring, `kdevs1` would become a free latent.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L, spnames = "a",
            nsex = 1L, growth_model = 1L, sigma_rec_prior = 1,
            random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, BT = 0))
  R1 <- "recdevs1 <-> recdevs1, 0, sR, 1\n"
  for (nm in c("kdevs1", "L1devs1", "l1dev1")) {
    testthat::expect_error(
      build_dsem_objects(dsem_settings = build_DSEM(
        sem = paste0(R1, nm, " <-> ", nm, ", 0, s, .3\n")), data_list = d),
      "missing from tsdata")
  }
})

testthat::test_that("an env_data column named like a latent warns before being replaced", {
  testthat::skip_if_not_installed("dsem")

  # The injection is unconditional and precedes the sem-variable select, so a
  # user's series under a latent's name is OVERWRITTEN with NA -- data deletion,
  # not just shadowing in the family regex. Measured before the warning: a
  # supplied 10..20 series came back all-NA and the column became a free latent.
  d <- list(styr = 1990, endyr = 2000, projyr = 2000, nspp = 1L, spnames = "a",
            nsex = 1L, growth_model = 1L, sigma_rec_prior = 1,
            random_rec = TRUE, proj_mean_rec = TRUE,
            env_data = data.frame(Year = 1990:2000, l1devs1 = 10:20))

  testthat::expect_warning(
    build_dsem_objects(
      dsem_settings = build_DSEM(sem = paste0(
        "recdevs1 <-> recdevs1, 0, sigmaR1, 1\n",
        "l1devs1 <-> l1devs1, 0, sigL1, 0.3\n")),
      data_list = d),
    "share a name with a DSEM latent")

  # A non-colliding env column is silent on this count.
  d2 <- d; d2$env_data <- data.frame(Year = 1990:2000, BT = 0)
  testthat::expect_no_warning(
    build_dsem_objects(
      dsem_settings = build_DSEM(sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1\n"),
      data_list = d2))
})
