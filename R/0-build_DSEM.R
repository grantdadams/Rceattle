
#' Function to fit a dynamic structural equation model related to recruitment
#'
#' @param sem Specification for the time-series structural equation model,
#'   including lagged and simultaneous effects; see \code{dsem::make_dsem_ram}.
#'   Every variable is either a column of \code{env_data} or one of the LATENT
#'   process columns below, which the SEM supplies and the population model then
#'   reads. The default is IID recruitment deviations.
#' @details
#' Latent process columns, one per species \code{<sp>} in \code{1:nspp}, with
#' \code{f}/\code{m} appended for a two-sex species and no token for a one-sex
#' species:
#'
#' \describe{
#'   \item{\code{recdevs<sp>}}{recruitment deviations}
#'   \item{\code{l1devs<sp>}}{\code{L1}, mean length at the growth anchor age}
#'   \item{\code{linfdevs<sp>}}{\code{Linf}, asymptotic mean length}
#'   \item{\code{richdevs<sp>}}{\code{m}, the Richards shape parameter}
#'   \item{\code{mdevs<sp>}}{\code{M1}, residual natural mortality}
#' }
#'
#' Every latent is on the LOG scale -- a multiplier on its target, not an
#' additive offset in the target's units. A growth latent of 0.2 multiplies the
#' parameter by \code{exp(0.2)} = 1.22; it does not add 0.2 cm.
#'
#' The growth latents are log-scale multipliers on the parameter, applied through
#' \code{growth_linkage_offset}, so they compose with a fixed environmental
#' growth linkage rather than replacing it. They are bias-corrected, which makes
#' \code{E[exp(latent)] = 1} and so leaves the PARAMETER's mean unbiased rather
#' than its median. Note what that does and does not buy: weight-at-age is a
#' convex power of length (\code{beta} near 3), so an unbiased mean parameter
#' still leaves mean weight-at-age -- the weight in SSB and B40\% -- biased
#' upward by Jensen's inequality. The correction removes the drift in the
#' parameter, not in the weight.
#'
#' The natural-mortality latent is the same kind of log-scale multiplier,
#' applied through \code{M_linkage_offset} and bias-corrected the same way, so
#' mean \code{M1} is unbiased rather than its median. The same Jensen caveat as
#' growth applies: survivorship is \code{exp(-M)}, convex in M, so an unbiased
#' mean M still leaves mean survivorship -- and so SB0 and B40\% -- biased. It
#' is FLAT over ages: one value per species, sex and year. With an
#' identity-link M linkage the latent multiplies the exponentiated part only,
#' since that linkage adds its offset outside the \code{exp()}. No code writes it into
#' predation mortality \code{M2}, but under \code{msmMode > 0} M2 moves with it
#' anyway, because M1 sets survivorship, survivorship sets the average numbers
#' predators eat from, and M2 is predator demand over available prey.
#'
#' Name a recruitment column for every species to have the SEM drive
#' recruitment, or none of them to leave recruitment on its own density; naming
#' only some is an error. A sem that names no latent column at all is refused,
#' since it would contribute a density to the objective while driving nothing.
#'
#' Naming one SEX of a two-sex species' growth and not the other is allowed, and
#' the asymmetry with the species rule above is deliberate. A species with no
#' recruitment column has no recruitment deviations at all -- the process is
#' gone -- whereas a sex with no growth column simply has time-invariant growth,
#' which is what every model without a growth latent has. The sexes' mean-growth
#' parameters are separate blocks (only the length-at-age SD is shared), so
#' driving the better-sampled sex alone is a coherent model rather than a
#' half-specified one.
#'
#' Growth \code{K} is NOT drivable. The Jan-1 recursion reads the previous
#' year's K, so a terminal-year K latent reaches a hindcast observation only
#' through the within-year advance, whose factor is zero at a Jan-1 fleet and
#' spawning month -- measured there, perturbing it leaves the objective
#' bit-identical. A growth latent is also refused where the species' growth model
#' never reads that parameter: any of them under empirical weight-at-age, and
#' \code{richdevs} under von Bertalanffy. An \code{mdevs} column is refused for
#' a species whose \code{M1_re} is on: both are deviations on the same
#' \code{M1}, so a flat annual multiplier beside an \code{M1_re} field asks the
#' model to split one quantity two ways.
#'
#' One consequence of the latent naming: a name matching the latent pattern is
#' treated as a latent even where the process never generates one, so an
#' \code{env_data} column called e.g. \code{recdevs1f} would be shadowed.
#' @param family Distribution for the measurement error of each \code{env_data}
#'   variable used in the \code{sem}. The latent process columns are never
#'   observed and are always treated as \code{fixed} (measured exactly), so
#'   \code{family} describes only the env-data variables retained in the
#'   \code{sem}. Supply a single string (recycled across all env variables) or a
#'   character vector with one element per env variable, each one of: \code{fixed}
#'   (default; no measurement error / identity link), \code{normal} (Gaussian,
#'   identity link), \code{gamma} (gamma, fixed CV and log link), \code{bernoulli}
#'   (logit link), \code{poisson} (log link), \code{lognormal} (log link), or
#'   \code{tweedie}. Translated internally to the named list of \code{family}
#'   objects that \code{dsem::dsem} (>= 3.0.0) expects. You may also pass
#'   \code{dsem} \code{family} objects directly (e.g.
#'   \code{dsem::gaussian_fixed_sd("identity", 0.1)} to fix an observation-error
#'   SD): a single object is recycled across env variables, and a vector or list
#'   mixing strings and objects is matched to env variables by name when named,
#'   otherwise by position.
#' @param sigmaR_prior_sd Optional log-scale SD of a lognormal prior placed on
#'   each species' estimated recruitment SD, centered at the assessment value
#'   \code{data_list$sigma_rec}. Default \code{NA} = no prior (recruitment
#'   SD estimated freely). Supplying a finite value (e.g. 0.5--1) regularizes the
#'   recruitment SD away from the \eqn{1/\sigma_R^2} collapse that occurs when
#'   environmental covariates over-explain the recruitment deviations, while
#'   still estimating it. The prior applies only where the SD is estimated (it is
#'   ignored for species whose SD is fixed in the \code{sem}).
#' @param estimate_projection Whether the SEM covers the projection period.
#'   Default \code{FALSE} builds the latent states over the hindcast
#'   (\code{styr:endyr}) only, leaving the projection years out of the DSEM
#'   likelihood; projection recruitment then comes from the usual method (mean
#'   recruitment or the stock-recruit relationship). \code{TRUE} extends them to
#'   \code{projyr} and projects recruitment through the SEM. This overrides
#'   \code{build_srr(proj_mean_rec = )}: under mean-recruitment projection the
#'   projected recruitment is \code{avg_R} and the SEM's projection states
#'   would reach only the dynamic B0/BF reference series, so \code{fit_mod()}
#'   turns mean-recruitment projection off and says so.
#'   \code{FALSE} leaves the projection years out, rather than holding them
#'   fixed: fixed states still sit in the GMRF, where lagged paths tie them to
#'   the last hindcast years and pull on the terminal recruitment deviations.
#' @param constant_variance What a \code{<->} term means when the sem carries
#'   lagged paths. Passed to \code{dsem::dsem_control}; the three settings are
#'   identical for a sem with no lags.
#'   \describe{
#'     \item{\code{"conditional"}}{the default, and \code{dsem}'s. A
#'       \code{<->} term is the DISTURBANCE (innovation) variance, in the
#'       structural-equation sense the RAM notation comes from. The marginal
#'       variance then grows along the causal graph -- for a first-order
#'       self-path it is \eqn{\sigma^2 / (1 - \rho^2)}, so it diverges as
#'       \eqn{\rho \to 1}.}
#'     \item{\code{"marginal"}}{rescales \code{Gamma} and \code{(I - Rho)} so the
#'       marginal variance is constant. This is the time-series convention used
#'       by \code{ar1(1 | Year)} elsewhere in Rceattle, and by WHAM's
#'       \code{Ecov} AR1 and \code{glmmTMB}'s \code{ar1()}, where the estimated
#'       SD is the marginal one.}
#'     \item{\code{"diagonal"}}{also holds the marginal variance constant, but
#'       preserves the conditional correlation rather than the conditional
#'       variance.}
#'   }
#'   The choice matters for the recruitment SD reported by \code{summary()} and
#'   for what \code{sigmaR_prior_sd} centres on: under \code{"conditional"} the
#'   prior lands on the innovation SD, while \code{data_list$sigma_rec} is an
#'   assessment's marginal recruitment SD. Those coincide for an IID sem, and
#'   differ by \eqn{1/(1 - \rho^2)} once a self-path is added.
#' @param bound_sd Whether to give the sem's standard-deviation paths a lower
#'   bound of 0. Default \code{TRUE}. A lag-0 two-headed self path
#'   (\code{x <-> x}) is a diagonal of the Cholesky factor of the exogenous
#'   covariance; the likelihood sees only \eqn{\Gamma'\Gamma}, so its sign is
#'   not identified and the surface is exactly symmetric about 0. That is a
#'   harmless pair of mirrored optima for an MLE and fatal for MCMC, where the
#'   posterior is bimodal by construction. Bounding the diagonal below at 0 is
#'   the standard identifying restriction for a Cholesky factor and rules out
#'   nothing the likelihood could distinguish, so leave it on unless you
#'   specifically want the unbounded parameterization -- reproducing a fit made
#'   before the bound existed, parameter for parameter, or demonstrating the
#'   symmetry. Only paths that qualify are bounded; \code{\link{build_bounds}}
#'   says which and why, and \code{fit_mod()} names any variable it had to
#'   skip. With \code{FALSE}, a starting value at the negative root is also left
#'   where it is rather than flipped.
#'
#' @description
#' The code links dynamic structural equation models to recruitment within Rceattle. The internals of \code{dsem} were copy and pasted into Rceattle. See \code{??dsem} for more description.
#'
#' @export
#'
build_DSEM <- function(sem = NULL,
                       family = "fixed",
                       sigmaR_prior_sd = NA,
                       estimate_projection = FALSE,
                       constant_variance = c("conditional", "marginal", "diagonal"),
                       bound_sd = TRUE
){
  constant_variance <- match.arg(constant_variance)
  if (!is.logical(bound_sd) || length(bound_sd) != 1L || is.na(bound_sd)) {
    stop("`bound_sd` must be TRUE or FALSE.", call. = FALSE)
  }
  return(list(sem = sem, family = family, sigmaR_prior_sd = sigmaR_prior_sd,
              estimate_projection = estimate_projection,
              constant_variance = constant_variance,
              bound_sd = bound_sd))
}


#' Parse a single SEM path specification
#'
#' Ported from \code{dsem::parse_path} (originally \code{sem::parse.path}). Used
#' by \code{\link{build_dsem_objects}} to collect the variables referenced in a
#' \code{sem}. dsem is distributed under GPL-3; \code{parse_path} originates from
#' package \code{sem} (GPL >= 2, used with permission from John Fox).
#'
#' @param path text to parse
#' @return tagged list with \code{first}, \code{second}, \code{direction}
#' @keywords internal
#' @noRd
parse_path <- function( path ){
  path.1 <- gsub("-", "", gsub(" ", "", path))
  direction <- if(regexpr("<>", path.1) > 0){
    2
  }else if(regexpr("<", path.1) > 0){
    -1
  }else if(regexpr(">", path.1) > 0){
    1
  }else{
    stop(paste("ill-formed path:", path))
  }
  path.1 <- strsplit(path.1, "[<>]")[[1]]
  list(first = path.1[1], second = path.1[length(path.1)], direction = direction)
}


#' Put one space either side of every arrow in a sem
#'
#' @description
#' `dsem::make_dsem_ram()` decides whether a variable already has a variance of
#' its own by comparing the path string it would add -- `paste(v, "<->", v)`,
#' single-spaced -- against the paths already in the model, with no whitespace
#' normalization. A sem whose columns are aligned (`recdevs1  <->  recdevs1`)
#' therefore fails to match, dsem appends a second variance row for the same
#' node, and `MakeADFun()` then dies inside the dsem DLL -- a segfault on a
#' multi-variable sem, an error on a single-variable one, either way with a
#' traceback that says nothing about the spacing. Normalizing on the way in
#' makes an aligned sem mean exactly what it reads as.
#'
#' `<->` is protected before `->` and `<-` are touched, so the three arrows
#' cannot rewrite each other's halves.
#'
#' @param sem the sem text, or `NULL`.
#' @return `sem` with each arrow surrounded by exactly one space.
#' @noRd
.dsem_normalize_arrows <- function(sem) {
  if (is.null(sem) || !length(sem)) return(sem)
  x <- gsub("[ \t]*<->[ \t]*", "\001", sem)
  x <- gsub("[ \t]*->[ \t]*", " -> ", x)
  x <- gsub("[ \t]*<-[ \t]*", " <- ", x)
  gsub("\001", " <-> ", x, fixed = TRUE)
}


#' Map an Rceattle DSEM family string to a dsem 3.0.0 family object
#'
#' dsem >= 3.0.0 expects glm-style \code{family} objects rather than the
#' character strings Rceattle's API uses. Returns the input unchanged if it is
#' already a \code{family} object.
#'
#' @param fam a single family string (or a \code{family} object)
#' @return a \code{family} object understood by \code{dsem::dsem}
#' @keywords internal
#' @noRd
dsem_family_object <- function( fam ){
  if(inherits(fam, "family")) return(fam)
  switch(as.character(fam),
         "fixed"     = dsem::fixed(),
         "normal"    = stats::gaussian(),
         "gaussian"  = stats::gaussian(),
         "bernoulli" = stats::binomial(),
         "binomial"  = stats::binomial(),
         "poisson"   = stats::poisson(),
         "gamma"     = stats::Gamma(link = "log"),
         "Gamma"     = stats::Gamma(link = "log"),
         "lognormal" = dsem::lognormal(),
         "tweedie"   = dsem::tweedie(),
         stop("Unsupported DSEM 'family': '", fam, "'. Use one of ",
              "fixed, normal, gamma, bernoulli, poisson, lognormal, tweedie.",
              call. = FALSE))
}


#' Registry of the CEATTLE process deviations DSEM can drive
#'
#' DSEM works by making a process's deviations the latent states of the SEM, so
#' several places need to agree on what those states are *called*: the default
#' sem, the columns injected into the tsdata, the regex that gives them the
#' `fixed` family, the SD self-loop lookup, the column-index map handed to the
#' C++, and the pre-fit spec check. This function is the single place that
#' decides, so adding a process is one edit here rather than a hunt through six
#' `paste0()`s and regexes.
#'
#' Each entry carries the sem variable-name `prefix` users write (e.g.
#' `recdevs` in `recdevs1 <-> recdevs1, 0, sigmaR1, 1`), the number of latent
#' series `n` that process contributes, and the resulting `columns` in model
#' order, plus the `species`, `sex` and growth-parameter `param` each column
#' drives. The template reads the growth maps through
#' `DATA_IVECTOR(growth_dev_*)` into `growth_linkage_offset` in section
#' 5.5b-growth, and the M map through `DATA_IVECTOR(mdev_*)` into
#' `M_linkage_offset` in section 5.5b-M.
#'
#' One entry per drivable TARGET, not per process: a sem names one latent column
#' per target, and growth has three (`L1`, `Linf`, Richards `m`). Growth
#' parameters are sex-specific, so a two-sex species carries one column per sex
#' with the `f`/`m` token [linkage_spec()] already accepts on input, while a
#' one-sex species carries no token -- the convention the package follows
#' everywhere (`parameter_index()` reports NA, `data_check()` uses ""). `param` is the 0-based slot in `growth_linkage_offset`'s
#' last dimension, taken from [LINKAGE_PARAM_CODES] so the two cannot drift.
#'
#' Growth `K` is deliberately absent, on a measured identification argument
#' rather than a structural one. In the Jan-1 cohort recursion year *t*'s length
#' reads the PREVIOUS year's K (`growth.hpp` `lagkappa`), so the terminal year's
#' K reaches a hindcast observation only through the within-year advance, whose
#' factor is `exp(-kappa * fracyr) - 1` -- exactly zero at `fracyr = 0`.
#' Measured there, perturbing the terminal K latent leaves the objective
#' *bit-identical*, against a 71-nat response one year earlier.
#'
#' That measurement holds under three conditions, and not otherwise: a Jan-1
#' fleet and spawning month; an INTEGER `growth_age_L1 - minage` (a fractional
#' anchor puts one age slot on the `last_linear` branch, which reads the current
#' year's kappa every year); and `plus_group_length` other than `"decay"`, whose
#' plus-group loop also reads the current year's kappa. At a mid-year month, a
#' fractional anchor, or `"decay"`, the terminal K is informed -- weakly.
#'
#' `L1` and `Linf` are drivable because each has a current-year channel that K
#' lacks, and they are different channels: `L1` enters the linear ramp, `Linf`
#' the plus group and the length-at-age SD denominator. Neither is universal
#' either -- `Linf`'s terminal information vanishes under
#' `plus_group_length = "none"` with equal SD endpoints.
#' A `kdevs<sp>` column is therefore not a latent name; it is currently dropped
#' as an unreferenced variable rather than refused, which is owed.
#' Natural mortality is `mdevs`; Richards `m` is `richdevs` so the two cannot
#' collide, and so that a trailing `m` is unambiguously the sex token. M's
#' `param` slot is `NA`: `M_linkage_offset` is dimensioned by AGE rather than by
#' growth parameter, so there is no slot to claim, and the latent is broadcast
#' flat over that age dimension.
#'
#' Declaring a target costs nothing: a latent column the sem does not name is
#' dropped before the model is built, so a recruitment-only sem is unaffected.
#'
#' This registry settles NAMING, not the per-process statistics, and the
#' recruitment-SD machinery below is still recruitment-specific: a growth
#' latent's `<->` variance takes its start value from the sem alone, gets no
#' `sigmaR_prior_sd` analogue, and has no `random_rec` map handling. Generalising
#' that is owed.
#'
#' @param data_list a data_list; only `nspp` and `nsex` are used. Declaration is
#'   unconditional: a `richdevs` column is declared even where the species has no
#'   Richards `m` to drive, and is dropped unless the sem names it.
#' @noRd
.dsem_latent_registry <- function(data_list) {
  nspp <- data_list[["nspp"]]
  spp  <- seq_len(nspp)
  reg <- list(
    recruitment = list(
      prefix  = .DSEM_LATENT_PREFIXES[["recruitment"]],
      n       = nspp,
      columns = paste0(.DSEM_LATENT_PREFIXES[["recruitment"]], spp),
      species = spp,
      sex     = rep(0L, nspp),          # recruitment is not sex-specific
      param   = rep(NA_integer_, nspp)
    )
  )

  # nsex may be absent on a partially built data_list; one sex is the safe
  # reading, and data_check() is what enforces its presence.
  nsex <- data_list[["nsex"]]
  if (is.null(nsex)) nsex <- rep(1L, nspp)
  nsex <- rep_len(as.integer(nsex), nspp)

  # One latent per (species, sex), with the one-sex rule: no sex token where the
  # species has one sex. Growth and natural mortality share this shape, which is
  # what earns it a helper rather than a second copy.
  by_sex <- function(prefix, param) {
    cols <- character(0); sp_i <- integer(0); sx_i <- integer(0)
    for (sp in spp) {
      toks <- if (nsex[sp] > 1L) c("f", "m") else ""
      cols <- c(cols, paste0(prefix, sp, toks))
      sp_i <- c(sp_i, rep(sp, length(toks)))
      sx_i <- c(sx_i, seq_along(toks))
    }
    list(prefix = prefix, n = length(cols), columns = cols,
         species = sp_i, sex = sx_i, param = rep(param, length(cols)))
  }

  for (nm in names(.DSEM_GROWTH_TARGETS)) {
    reg[[nm]] <- by_sex(
      .DSEM_GROWTH_TARGETS[[nm]],
      unname(LINKAGE_PARAM_CODES$growth[[sub("^growth_", "", nm)]]))
  }

  # Natural mortality. `param` is NA, not LINKAGE_PARAM_CODES$M's 0: M's offset
  # tensor is `M_linkage_offset(nspp, max_sex, max_age, nyrs)` -- there is no
  # parameter dimension to index, because M1 is the only drivable M parameter and
  # the AGE dimension takes that slot's place. An M latent is one series per
  # (species, sex) BROADCAST over ages, so the C++ needs species and sex and
  # nothing else.
  #
  # The sex rule follows growth's rather than `M1_re`'s: M random effects run per
  # sex only at `M1_model == 2` (ceattle.cpp 13.1, jnll slot 15), but that switch governs
  # whether the BASELINE M1 is estimated per sex, not whether an offset can be.
  # `log_M1` and `M_linkage_offset` are both sex-dimensioned regardless, so a
  # shared baseline with sex-specific time variation is a coherent model, and one
  # naming rule across processes is what a user actually needs.
  reg[["natural_mortality"]] <-
    by_sex(.DSEM_LATENT_PREFIXES[["natural_mortality"]], NA_integer_)
  reg
}

# TMB requires every DATA_ and PARAMETER_ macro in the template to be supplied,
# whether or not the model uses it. These are the stand-ins for a fit with no
# DSEM: zero-length everywhere, dsem_on = 0. That keeps fit_mod() on one code
# path -- the switch is the only branch -- and keeps the non-DSEM objective
# untouched, since nothing zero-length can contribute to it.
#' @noRd
.dsem_null_inputs <- function() {
  list(
    data = list(
      dsem_on           = 0L,
      # 0, so the recruitment-deviation density in section 13 is reached on the
      # branch a model with no DSEM has always taken.
      dsem_drives_rec   = 0L,
      nyrs_dsem         = 0L,
      dsem_options      = integer(4),
      dsem_RAM          = matrix(0L, 0, 6),
      dsem_RAMstart     = numeric(0),
      dsem_familycode_j = integer(0),
      dsem_linkcode_j   = integer(0),
      dsem_sigmastart_j = integer(0),
      dsem_eps_tj       = array(0, dim = c(0, 0)),
      dsem_y_tj         = array(0, dim = c(0, 0)),
      dsem_obs_idx      = integer(0),
      dsem_unobs_idx    = integer(0),
      rec_dev_col       = integer(0),
      # Growth latent index map, empty for a model with no growth latent.
      # Supplied here because TMB requires every DATA_ macro in the template,
      # and section 5.5b-growth declares these; the null and real input field
      # sets have to stay identical or the two paths hand TMB different lists.
      growth_dev_col    = integer(0),
      growth_dev_sp     = integer(0),
      growth_dev_sex    = integer(0),
      growth_dev_par    = integer(0),
      mdev_col          = integer(0),
      mdev_sp           = integer(0),
      mdev_sex          = integer(0),
      rec_sd_idx        = integer(0),
      rec_sd_fixed      = numeric(0),
      rec_sd_prior      = numeric(0),
      rec_sd_prior_sd   = numeric(0),
      rec_sd_use_prior  = integer(0)
    ),
    parameters = list(
      dsem_beta_z    = numeric(0),
      dsem_lnsigma_z = numeric(0),
      dsem_mu_j      = numeric(0),
      dsem_delta0_j  = numeric(0),
      dsem_x_tj      = array(0, dim = c(0, 0))
    )
  )
}

# The DSEM parameter blocks in a shape merge_dsem_params() can consume: named
# for the template, with the built starting values.
#' @noRd
.dsem_par_shape <- function(dsem) {
  list(tmb_inputs = list(parameters = .dsem_par_template(dsem)))
}

# A warm start must match the DSEM this fit is building. inits from another fit
# -- a different year span, a different sem -- would otherwise be reinterpreted
# silently against the new dimensions.
#' @noRd
.dsem_check_shape <- function(param_list, dsem) {
  if (is.null(dsem)) return(param_list)
  want <- .dsem_par_template(dsem)
  for (nm in names(want)) {
    got <- param_list[[nm]]
    # Absent, or present-but-empty (a warm start from a non-DSEM fit carries the
    # zero-length blocks): nothing to conflict with, take the template's.
    if (is.null(got) || length(got) == 0L) {
      param_list[[nm]] <- want[[nm]]
      next
    }
    if (!identical(dim(got), dim(want[[nm]])) || length(got) != length(want[[nm]])) {
      stop("`inits` carries a `", nm, "` of ",
           paste(dim(got) %||% length(got), collapse = "x"),
           " but this DSEM needs ",
           paste(dim(want[[nm]]) %||% length(want[[nm]]), collapse = "x"),
           ". Rebuild the starting values for this model (inits = NULL).",
           call. = FALSE)
    }
  }
  param_list
}

# DSEM map entries, renamed onto the template's dsem_-prefixed parameter names.
# build_dsem_objects() names them as dsem::dsem() does (beta_z, x_tj, ...); the
# CEATTLE template prefixes them so they cannot collide with anything already in
# the model. The rename has to happen in one place or the map silently fails to
# match the parameters.
#' @noRd
.dsem_map_entries <- function(dsem) {
  if (is.null(dsem)) return(list(mapList = list(), mapFactor = list()))
  ren <- function(x) { names(x) <- paste0("dsem_", names(x)); x }
  list(mapList   = ren(dsem$mapList),
       mapFactor = ren(dsem$tmb_inputs$map))
}

# The DSEM parameter blocks for a fit, in the template's naming. Zero-length
# when there is no DSEM.
#' @noRd
.dsem_par_template <- function(dsem) {
  if (is.null(dsem)) .dsem_null_inputs()$parameters else .dsem_tmb_inputs(dsem)$parameters
}

# The same objects for a fit that HAS a DSEM, mapped from build_dsem_objects()
# onto the template's dsem_-prefixed names.
#' @noRd
.dsem_tmb_inputs <- function(dsem) {
  d <- dsem$tmb_inputs$data
  p <- dsem$tmb_inputs$parameters
  list(
    data = list(
      dsem_on           = 1L,
      # Defaults to 1 for a `dsem` object built before the flag existed: such a
      # SEM could only have driven recruitment, and every saved $dsem in the
      # assessment repos refits through here. .dsem_drives() takes the same
      # default, or a saved fit would refit as a different model than it reports.
      dsem_drives_rec   = as.integer(d$drives_rec %||% 1L),
      nyrs_dsem         = as.integer(d$nyrs_dsem),
      dsem_options      = as.integer(d$options),
      dsem_RAM          = matrix(as.integer(as.matrix(d$RAM)), ncol = 6),
      dsem_RAMstart     = as.numeric(d$RAMstart),
      dsem_familycode_j = as.integer(d$familycode_j),
      dsem_linkcode_j   = as.integer(d$linkcode_j),
      dsem_sigmastart_j = as.integer(d$sigmastart_j),
      dsem_eps_tj       = array(as.numeric(d$eps_tj), dim = dim(d$eps_tj)),
      dsem_y_tj         = array(as.numeric(d$y_tj),   dim = dim(d$y_tj)),
      dsem_obs_idx      = if (is.null(d$obs_idx))   integer(0) else as.integer(d$obs_idx),
      dsem_unobs_idx    = if (is.null(d$unobs_idx)) integer(0) else as.integer(d$unobs_idx),
      rec_dev_col       = as.integer(d$rec_dev_col),
      # Empty for a `dsem` object built before growth was drivable: such a SEM
      # drove recruitment only, so it has no growth latents. Every saved $dsem
      # in the assessment repos refits through here, so the default has to be
      # the no-growth case rather than an error.
      growth_dev_col    = as.integer(d$growth_dev_col %||% integer(0)),
      growth_dev_sp     = as.integer(d$growth_dev_sp  %||% integer(0)),
      growth_dev_sex    = as.integer(d$growth_dev_sex %||% integer(0)),
      growth_dev_par    = as.integer(d$growth_dev_par %||% integer(0)),
      # Empty for any `dsem` saved before M was drivable.
      mdev_col          = as.integer(d[["mdev_col"]] %||% integer(0)),
      mdev_sp           = as.integer(d[["mdev_sp"]] %||% integer(0)),
      mdev_sex          = as.integer(d[["mdev_sex"]] %||% integer(0)),
      rec_sd_idx        = as.integer(d$rec_sd_idx),
      rec_sd_fixed      = as.numeric(d$rec_sd_fixed),
      rec_sd_prior      = as.numeric(d$rec_sd_prior),
      rec_sd_prior_sd   = as.numeric(d$rec_sd_prior_sd),
      rec_sd_use_prior  = as.integer(d$rec_sd_use_prior)
    ),
    parameters = list(
      dsem_beta_z    = as.numeric(p$beta_z),
      dsem_lnsigma_z = as.numeric(p$lnsigma_z),
      dsem_mu_j      = as.numeric(p$mu_j),
      dsem_delta0_j  = as.numeric(p$delta0_j),
      dsem_x_tj      = array(as.numeric(p$x_tj), dim = dim(p$x_tj))
    )
  )
}

# The DSEM parameter block names, as dsem::dsem() emits them. Used by fit_mod()
# to strip a DSEM warm start when the current fit has no DSEM: these are absent
# from build_params(), so the inits shape guard passes them through as unflagged
# extras and MakeADFun then rejects parameters the template never declares.
# Kept as a constant rather than derived from a built object, because the scrub
# has to happen precisely when no DSEM object exists.
.DSEM_PARAM_NAMES <- paste0("dsem_", c("beta_z", "lnsigma_z", "mu_j",
                                       "delta0_j", "x_tj"))

# The sem variable-name stems, one per drivable process. Deliberately separate
# from the registry above: the stems are a property of the grammar, not of any
# particular model, so the pattern below can be built without a data_list.
# Folding them together made check_dsem_spec() -- which needs only the stems --
# start requiring data_list$nspp, and error on inputs it used to accept.
.DSEM_LATENT_PREFIXES <- c(recruitment       = "recdevs",
                           growth_L1         = "l1devs",
                           growth_Linf       = "linfdevs",
                           growth_m          = "richdevs",
                           natural_mortality = "mdevs")

# The growth subset, keyed by registry entry name. Richards `m` takes the stem
# `richdevs`, not `mdevs`: `mdevs` belongs to natural mortality, and a trailing
# `m` has to stay readable as the sex token.
.DSEM_GROWTH_TARGETS <- .DSEM_LATENT_PREFIXES[c("growth_L1", "growth_Linf",
                                                "growth_m")]


# The first RANDOM linkage spec under `keys` (every key when NULL) that covers
# any latent at 1-based species `sp` and sex `sex`, with `hit` over `sp`; NULL
# if there is none. A spec's `species` holds ids or names, matched to `spnames`.
#' @noRd
.random_linkage_cover <- function(linkages, sp, spnames, keys = NULL, sex = NULL) {
  groups <- if (is.null(keys)) linkages else
    linkages[intersect(keys, names(linkages))]
  for (specs in groups) {
    # A key holds one spec, or a list of them (the per-species-formula form).
    if (inherits(specs, "Rceattle_linkage_spec")) specs <- list(specs)
    for (spec in specs) {
      if (!inherits(spec, "Rceattle_linkage_spec")) {
        stop("A linkage in `data_list` is not a linkage_spec() object. Build ",
             "linkages with linkage_spec() inside build_M1() or build_growth().",
             call. = FALSE)
      }
      if (!length(.parse_linkage_formula(spec[["formula"]])$re_terms)) next
      # A spec with no `species` applies to every species.
      cover <- spec[["species"]]
      if (is.character(cover)) cover <- match(trimws(cover), trimws(spnames))
      hit <- if (is.null(cover)) rep(TRUE, length(sp)) else sp %in% cover
      # `sex` restricts a spec only where it is stratified by sex (`by`).
      if (!is.null(sex) && !is.null(spec[["sex"]]) &&
          "sex" %in% all.vars(spec[["by"]])) {
        hit <- hit & sex %in% spec[["sex"]]
      }
      if (any(hit)) return(list(spec = spec, hit = hit))
    }
  }
  NULL
}


# Run the per-process refusals on a built DSEM's own index maps, so they cover
# objects from build_dsem_objects() handed to fit_mod() as well as a spec.
#' @param data_list the data_list the model is fitted with.
#' @param obj built DSEM objects carrying `$tmb_inputs`.
#' @noRd
.check_dsem_support <- function(data_list, obj) {
  dt  <- obj[["tmb_inputs"]][["data"]]
  nm  <- colnames(dt[["y_tj"]])
  # Maps of disagreeing lengths are the template's to refuse, with its own
  # message (ceattle.cpp sections 5.5b-growth and 5.5b-M), so a process whose
  # maps do not line up is left to it.
  aligned <- function(f) length(unique(lengths(lapply(f, function(x) dt[[x]])))) == 1L
  g <- dt[["growth_dev_col"]] %||% integer(0)
  if (aligned(c("growth_dev_col", "growth_dev_sp", "growth_dev_sex", "growth_dev_par"))) {
    .check_dsem_growth_model(data_list, nm[g + 1L], dt[["growth_dev_sp"]] + 1L,
                             dt[["growth_dev_par"]], dt[["growth_dev_sex"]] + 1L)
  }
  m <- dt[["mdev_col"]] %||% integer(0)
  if (aligned(c("mdev_col", "mdev_sp", "mdev_sex"))) {
    .check_dsem_m_support(data_list, nm[m + 1L], dt[["mdev_sp"]] + 1L,
                          dt[["mdev_sex"]] + 1L)
  }
}


# Refuse a growth latent the species' growth model never reads. Empirical
# weight-at-age (`growth_model` 0) has no growth parameters at all, and von
# Bertalanffy (1) has no Richards `m` -- `R/3-build_map.R` maps `log_m` out at
# npar = 3. Such a latent is a free field carrying its own GMRF density and
# driving nothing: the model converges and reports an ABC for a process the user
# believes is time-varying and is not. Refused rather than warned because the
# user named the column explicitly.
#
# Also refuse one beside a RANDOM growth linkage on the same parameter: two
# random fields on it, each scored by its own density. A FIXED linkage composes.
#
# `growth_model` is set by fit_mod() from `growthFun`; a direct caller of
# build_dsem_objects() may not have it, and there is nothing to check against then.
#' @param data_list carries `growth_model`, `growth_linkages` and `spnames`.
#' @param cols,sp,par,sex The named growth latents, their 1-based species, their
#'   `LINKAGE_PARAM_CODES$growth` slot, and their 1-based sex (NULL for any).
#' @noRd
.check_dsem_growth_model <- function(data_list, cols, sp, par, sex = NULL) {
  if (!length(cols)) return(invisible(NULL))
  codes <- LINKAGE_PARAM_CODES$growth
  for (p in unique(par)) {
    on_p <- par == p
    key  <- names(codes)[codes == p]
    rnd  <- .random_linkage_cover(data_list[["growth_linkages"]], sp[on_p],
                                  data_list[["spnames"]], keys = key,
                                  sex = sex[on_p])
    if (is.null(rnd)) next
    stop("The sem drives growth (", paste(cols[on_p][rnd$hit], collapse = ", "),
         ") for species that also carry a RANDOM ", key, " growth linkage (",
         deparse(rnd$spec[["formula"]]), "). Each is a random field on ", key,
         ", scored by its own density, so ", key, " would carry two with ",
         "nothing to tell them apart. Use the sem or the linkage, not both. A ",
         "FIXED growth linkage does compose with the latent and is not refused.",
         call. = FALSE)
  }

  growth_model <- data_list[["growth_model"]]
  if (is.null(growth_model)) return(invisible(NULL))
  gm <- rep_len(as.integer(growth_model), max(sp))
  m_slot <- unname(LINKAGE_PARAM_CODES$growth[["m"]])
  bad_empirical <- gm[sp] == 0L
  bad_richards  <- gm[sp] == 1L & par == m_slot
  if (any(bad_empirical)) {
    stop("The sem drives growth (", paste(cols[bad_empirical], collapse = ", "),
         ") for species with empirical weight-at-age (growth_model = 0), which ",
         "estimates no growth parameters. Those latents would be fitted and ",
         "drive nothing. Set a parametric growth model for those species ",
         "(build_growth(fun = \"vonBertalanffy\")), or drop the columns from ",
         "the sem.", call. = FALSE)
  }
  if (any(bad_richards)) {
    stop("The sem drives the Richards shape parameter (",
         paste(cols[bad_richards], collapse = ", "),
         ") for species whose growth model is von Bertalanffy, which has no ",
         "`m` -- it is mapped out. That latent would be fitted and drive ",
         "nothing. Use build_growth(fun = \"Richards\") for those species, or ",
         "drop the columns from the sem.", call. = FALSE)
  }
  invisible(NULL)
}

# Refuse a SEM M latent beside another time-varying M process for the same
# species: M1_re, which fills `log_M1_dev`, or a RANDOM M1 linkage. Each enters
# M1_at_age (ceattle.cpp section 6.1) beside the latent and is scored by its own
# density, so two fields compete for the same signal.
#
# Modes 2 and 5 (by year, constant over ages) are the DSEM latent's shape
# exactly; modes 3 and 6 vary by age AND year, so they contain that shape. Both
# must stay refused. Modes 1 and 4 (by age, constant over years) are orthogonal
# in time -- `M1(a) x M_year(t)` is separable and arguably legitimate -- so they
# are refused for caution rather than for a conflict, and the message says which
# is which so relaxing it later is a documented step.
# There is deliberately NO estDynamics refusal here, where
# .check_dsem_growth_model() has one for growth. The two cases look alike and are
# not: a growth latent under growth_model = 0 drives nothing because no growth
# parameter exists, whereas a fixed-numbers species still has mortality. Its
# numbers-at-age are input, but the Baranov catch equation reads Z, so the catch
# likelihood still informs the latent -- weaker than for an estimated-N species,
# not inert. Measured figures in inst/dev/TRAPS.md.
#' @param data_list carries `M1_re`, `M1_linkages` and `spnames`; absent means
#'   no M random effects.
#' @param cols,sp,sex the named M latents, their 1-based species, and their
#'   1-based sex (NULL for any).
#' @noRd
.check_dsem_m_support <- function(data_list, cols, sp, sex = NULL) {
  if (!length(cols)) return(invisible(NULL))

  # A FIXED M1 linkage is the intended pairing -- the covariate carries the
  # explained part, the latent the residual -- and composes via `+=`.
  rnd <- .random_linkage_cover(data_list[["M1_linkages"]], sp,
                               data_list[["spnames"]], sex = sex)
  if (!is.null(rnd)) {
    stop("The sem drives natural mortality (", paste(cols[rnd$hit], collapse = ", "),
         ") for species that also carry a RANDOM M1 linkage (",
         deparse(rnd$spec[["formula"]]), "). Each is a random field on M1, ",
         "scored by its own density, so M would carry two with nothing to tell ",
         "them apart. Use the sem or the linkage, not both. A FIXED M1 linkage ",
         "does compose with the latent and is not refused.", call. = FALSE)
  }

  m1_re <- data_list[["M1_re"]]
  if (is.null(m1_re)) return(invisible(NULL))
  # build_M1() accepts the string form ("iid_year"), and data_check() tests for
  # it, so a data_list in the wild carries it. as.integer() would give NA and the
  # `any(bad)` below would then error on a missing value with no useful message.
  m1_re <- rep_len(.coerce_M1_arg(m1_re, .M1_RES, "M1_re"), max(sp))
  bad <- m1_re[sp] > 0L
  if (!any(bad)) return(invisible(NULL))
  yr_modes <- m1_re[sp] %in% unname(.M1_RES[grepl("year", names(.M1_RES))])
  stop("The sem drives natural mortality (", paste(cols[bad], collapse = ", "),
       ") for species whose `M1_re` is already non-zero (",
       paste(m1_re[sp][bad], collapse = ", "), "). Both write the same ",
       "exponent -- M1_re through log_M1_dev and the SEM through ",
       "M_linkage_offset -- so M would be driven twice, each part scored by its ",
       "own density. ",
       if (any(yr_modes & bad))
         paste0("M1_re 2/5 vary M by year alone, which is this latent's shape ",
                "exactly, and 3/6 vary M by age and year, which contains it. ") else
         paste0("M1_re 1/4 vary M by age only, so they are separable from a ",
                "year latent in principle; they are refused here for caution, ",
                "not for a measured conflict. "),
       "Set M1_re = 0 to drive M through the SEM.", call. = FALSE)
}



# Latent column names across the whole registry, in model order.
#' @noRd
.dsem_latent_columns <- function(data_list) {
  unlist(lapply(.dsem_latent_registry(data_list), `[[`, "columns"),
         use.names = FALSE)
}

# Regex matching any latent column name. These are states, never observations,
# so they take dsem's `fixed` family regardless of what the user asked for.
# Takes no arguments by design -- see the note on .DSEM_LATENT_PREFIXES.
#' @noRd
.dsem_latent_pattern <- function() {
  # Optional trailing f/m is the sex token a two-sex species carries. No prefix
  # is a prefix of another, so the alternation is unambiguous; the suffix class
  # stays narrow enough to keep rejecting near-misses like `recdevs_anom`,
  # `recdevs1_lag` or `xrecdevs1`. One accepted consequence: `recdevs1f` is now
  # RECOGNISED as a latent name though recruitment never generates one, so an
  # env column by that name would be shadowed.
  paste0("^(", paste(unique(.DSEM_LATENT_PREFIXES), collapse = "|"),
         ")[0-9]+[fm]?$")
}


#' Function to build the map and parameter objects for DSEM recruitment linkages
#'
#' @param dsem_settings dsem specifications from \code{\link{build_DSEM}}.
#' @param debug Runs the model without estimating parameters to get derived quantities given initial parameter values. If TRUE, sets all map values to NA
#' @param data_list a data_list read in via \code{\link{read_data}} or built directly in R.
#'
#' @export
build_dsem_objects <- function(dsem_settings = NULL, debug = FALSE, data_list = NULL){

  latent_cols <- .dsem_latent_columns(data_list)
  rec_cols    <- .dsem_latent_registry(data_list)$recruitment$columns

  # Build IID sem if NULL
  if(is.null(dsem_settings$sem)){
    sem = c()
    for(sp in seq_len(data_list$nspp)){
      nm <- rec_cols[sp]
      sem <- c(sem, paste0(nm, " <-> ", nm, ", 0, sigmaR", sp,", 1\n")) # No space after
    }
    sem <- paste0(sem, collapse = " ")
    dsem_settings$sem <- sem
  }

  # Put one space either side of every arrow, so a sem written with its columns
  # aligned means the same thing as one written compactly.
  dsem_settings$sem <- .dsem_normalize_arrows(dsem_settings$sem)

  # DSEM time span. With estimate_projection = FALSE the state space is built
  # over the HINDCAST ONLY, so projection years are absent from the GMRF
  # entirely rather than present-but-pinned.
  #
  # This matters: the density is a quadratic form over the stacked [t, j] state
  # space, and lagged RAM paths couple each row to its predecessors. Merely
  # mapping the projection rows off (the previous behavior) fixes those states
  # at their initial value but leaves them in the quadratic form, so their
  # cross-terms with the terminal hindcast years are still not constant in the
  # estimated states -- the terminal recruitment deviations pick up an extra
  # shrinkage pull toward whatever makes the pinned projection states likely.
  # Terminal recdevs feed terminal SSB and hence the ABC, so this is not a
  # cosmetic difference. Truncating the span removes the coupling at the source.
  dsem_endyr <- if (isTRUE(dsem_settings$estimate_projection)) {
    data_list$projyr
  } else {
    data_list$endyr
  }

  # DSEM data
  dsem_data <- data_list$env_data %>%
    # Adding NA in missing years (match assessment begining)
    dplyr::full_join(data.frame(Year=c(data_list$styr:dsem_endyr)), by = dplyr::join_by(Year)) %>%
    dplyr::filter(Year >= data_list$styr & Year <= dsem_endyr) %>% # FIXME: if including init devs, adjust
    dplyr::arrange(Year)

  # - Prepend one all-NA latent column per registry entry, in model order. All-NA
  #   is what marks them as states rather than observations downstream.
  # Base R rather than mutate(!!nm := ): the tidy-eval `:=` is not an object this
  # package imports, so R CMD check reports it as an undefined global. Column
  # ORDER is load-bearing -- rec_dev_col indexes x_tj by position -- so each
  # latent column goes to the front, and iterating in reverse leaves them in
  # model order.
  # An env_data column sharing a latent's name is OVERWRITTEN with NA here, not
  # merely shadowed in the family regex -- the user's series is deleted and the
  # column becomes a free latent state. Warn before doing it: growth added three
  # more name families (l1devs/linfdevs/richdevs) plus the f/m sex forms, so the
  # chance of a collision with a real covariate went up. No sibling script is
  # affected today (swept), which is why this warns rather than refuses.
  .clobbered <- intersect(latent_cols, names(dsem_data))
  if (length(.clobbered)) {
    warning("env_data column(s) ", paste(.clobbered, collapse = ", "),
            " share a name with a DSEM latent process column, so the supplied ",
            "series are replaced by latent states and the data in them is not ",
            "used. Rename those env_data columns.", call. = FALSE)
  }
  for(nm in rev(latent_cols)){
    dsem_data[[nm]] <- NA_real_
    dsem_data <- dsem_data[, c(nm, setdiff(names(dsem_data), nm)), drop = FALSE]
  }

  # - Drop the Year column once all recdev columns are added
  dsem_data <- dsem_data %>%
    dplyr::select(-Year)

  # - Keep only variables referenced in the sem.
  # Mirror make_dsem_ram(): scan the sem to a model table, then parse_path each
  # path to collect the variables on both heads of every arrow.
  sem_model <- scan(
    text         = dsem_settings$sem,
    what         = list(path = "", lag = 1, par = "", start = 1, dump = ""),
    sep          = ",",
    strip.white  = TRUE,
    comment.char = "#",
    fill         = TRUE,
    quiet        = TRUE
  )
  sem_vars <- unique(unlist(lapply(sem_model$path, function(p){
    pp <- parse_path(p)
    c(pp$first, pp$second)
  })))

  # An ESTIMATED lag-0 two-headed self path starting at exactly 0 puts the AD
  # tape's start on a singular exogenous covariance, and MakeADFun() inside
  # dsem::dsem() below does not error -- it SEGFAULTS, and R aborts with no
  # condition for anything to catch. check_dsem_spec() cannot help either: it
  # screens this function's OUTPUT, so the crash happens before it is reachable.
  # Refuse here, while there is still a session to report into.
  #
  # The parameter NAME is what separates this from the legitimate case, not the
  # start value. `x <-> x, 0, NA, 0` pins the variance at zero deliberately: the
  # variable becomes deterministic, dsem handles it, and check_dsem_spec()'s
  # covariate_variance record is what reports it (see test-dsem-spec-variance.R).
  # `x <-> x, 0, sigma, 0` instead asks to ESTIMATE a variance from a start of
  # zero, which is the shape that crashes. Measured against dsem 3.0.0: named +
  # 0 aborts R; NA + 0 returns a usable model on both a covariate and a recdev
  # column. A missing start is fine (dsem supplies its own), a zero CROSS
  # covariance is a legitimate "these two do not covary", and a lag-1 two-headed
  # path dsem already refuses with a clear message.
  .fixed_par <- function(p) {
    is.na(p) || !nzchar(trimws(p)) || toupper(trimws(p)) %in% c("NA", "NAN")
  }
  .zero_var <- vapply(seq_along(sem_model$path), function(i) {
    pp <- parse_path(sem_model$path[i])
    isTRUE(pp$direction == 2) && isTRUE(pp$first == pp$second) &&
      isTRUE(sem_model$lag[i] == 0) &&
      !.fixed_par(sem_model$par[i]) &&
      isTRUE(is.finite(sem_model$start[i])) && isTRUE(sem_model$start[i] == 0)
  }, logical(1))
  if (any(.zero_var)) {
    bad <- unique(vapply(which(.zero_var),
                         function(i) parse_path(sem_model$path[i])$first,
                         character(1)))
    stop("sem asks to estimate a two-headed self path (`<->`) for ",
         paste(sprintf("'%s'", bad), collapse = ", "),
         " from a start value of 0. A variance cannot be estimated from zero: ",
         "the covariance is singular there, and dsem crashes rather than ",
         "returning an error. Give it a positive start value. To pin the ",
         "variance at zero deliberately -- making the variable deterministic -- ",
         "write the parameter name as NA instead, which check_dsem_spec() then ",
         "reports.", call. = FALSE)
  }

  # Preserve dsem_data's column order (recdevs<sp> first, then env vars)
  dsem_data <- dsem_data %>%
    dplyr::select(dplyr::any_of(intersect(colnames(dsem_data), sem_vars)))


  # DSEM family. dsem >= 3.0.0 wants a named list of `family` objects keyed by
  # tsdata column. The latent recdevs<sp> columns are never observed, so they
  # are always `fixed`; the user `family` applies only to the env-data variables
  # retained in the sem. `family` may be:
  #   * a single string or `family` object  -> recycled across all env variables
  #   * a vector/list of strings and/or `family` objects, one per env variable,
  #     matched to env columns by name when named, otherwise by position.
  var_names <- colnames(dsem_data)
  is_recdev <- grepl(.dsem_latent_pattern(), var_names)
  env_names <- var_names[!is_recdev]

  fam_in <- dsem_settings$family
  # A bare `family` object is itself a list; wrap it so it counts as one entry.
  if(inherits(fam_in, "family")) fam_in <- list(fam_in)
  # Normalize so character vectors and lists share one code path.
  fam_in  <- as.list(fam_in)
  fam_nms <- names(fam_in)

  if(length(fam_in) == 1L && is.null(fam_nms)){
    # Single unnamed family -> recycle across every env variable.
    env_family <- stats::setNames(rep(fam_in, length(env_names)), env_names)
  } else if(!is.null(fam_nms) && all(nzchar(fam_nms))){
    # Fully named -> match by env-column name (order-independent).
    unknown <- setdiff(fam_nms, env_names)
    if(length(unknown) > 0){
      stop("'family' names are not env-data variables in the sem (",
           paste(env_names, collapse = ", "), "): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    missing_fam <- setdiff(env_names, fam_nms)
    if(length(missing_fam) > 0){
      stop("'family' is missing an entry for env-data variable(s): ",
           paste(missing_fam, collapse = ", "), ".", call. = FALSE)
    }
    env_family <- fam_in[env_names]
  } else {
    # Unnamed vector/list -> matched to env columns by position.
    if(length(fam_in) != length(env_names)){
      stop("Length of 'family' must be 1 or the number of env-data variables in ",
           "the sem (", length(env_names), "), but got ", length(fam_in), ".",
           call. = FALSE)
    }
    env_family <- stats::setNames(fam_in, env_names)
  }

  dsem_family <- stats::setNames(vector("list", length(var_names)), var_names)
  dsem_family[is_recdev] <- list(dsem::fixed())
  dsem_family[env_names] <- lapply(env_family[env_names], dsem_family_object)

  # The compiled TMB model (src/TMB/dsem.hpp) are
  # matched to dsem 3.0.0. Other dsem versions can emit a different
  # `tmb_inputs$data` layout (e.g. missing `obs_idx`/`unobs_idx` under the
  # gmrf_project parameterization), which surfaces downstream as the cryptic TMB
  # error "Error when reading the variable: 'obs_idx'". Fail fast instead.
  if (!requireNamespace("dsem", quietly = TRUE)) {
    stop("The 'dsem' package (version 3.0.0) is required to build DSEM inputs. ",
         "Install the matching version with:\n",
         "  remotes::install_version(\"dsem\", version = \"3.0.0\")",
         call. = FALSE)
  }
  if (utils::packageVersion("dsem") != "3.0.0") {
    stop("Rceattle's DSEM module requires dsem 3.0.0 (the version src/TMB/dsem.hpp ",
         "is matched to), but dsem ", utils::packageVersion("dsem"),
         " is installed.\n",
         "Other versions can produce an incompatible data layout and a downstream ",
         "TMB error such as \"Error when reading the variable: 'obs_idx'\".\n",
         "Install the matching version with:\n",
         "  remotes::install_version(\"dsem\", version = \"3.0.0\")",
         call. = FALSE)
  }

  # Two variance terms on one node is the failure the arrow normalization above
  # exists to prevent, and it is fatal inside the dsem DLL rather than reported.
  # Catch it here, on the model table dsem itself builds, so any OTHER route to
  # it -- a sem that genuinely names the same self-loop twice, or a future dsem
  # that adds a default row for a different reason -- is named rather than
  # crashing the session.
  .ram <- dsem::make_dsem_ram(dsem_settings$sem, times = seq_len(nrow(dsem_data)),
                              variables = colnames(dsem_data), quiet = TRUE)
  .self <- .ram$model[abs(as.numeric(.ram$model$direction)) == 2 &
                        .ram$model$first == .ram$model$second &
                        as.numeric(.ram$model$lag) == 0, , drop = FALSE]
  .dup <- unique(.self$first[duplicated(.self$first)])
  if (length(.dup) > 0) {
    stop("Variable(s) ", paste(.dup, collapse = ", "), " carry two variance ",
         "terms (`<->` self-loops at lag 0) in this sem: ",
         paste(sprintf("'%s'", .self$name[.self$first %in% .dup]), collapse = ", "),
         ". One variable takes one variance. Remove the duplicate row, or, if ",
         "you did not write one, leave the variance out entirely and dsem will ",
         "estimate it.", call. = FALSE)
  }

  # constant_variance decides what a `<->` term means once the sem carries lagged
  # paths -- the disturbance variance (dsem's default, and the SEM reading) or
  # the marginal one. Defaulted here as well as in build_DSEM(), so a spec
  # assembled by hand or restored from a config written before this argument
  # existed keeps dsem's default rather than becoming NULL.
  .cv <- dsem_settings$constant_variance
  if (is.null(.cv) || !length(.cv)) .cv <- "conditional"
  .cv <- match.arg(as.character(.cv)[1],
                   c("conditional", "marginal", "diagonal"))

  fit_dsem <- dsem::dsem(sem = dsem_settings$sem,
                         tsdata = stats::ts(dsem_data),
                         family = dsem_family,
                         control = dsem::dsem_control(use_REML = FALSE,
                                                 quiet = TRUE,
                                                 run_model = FALSE,
                                                 constant_variance = .cv))

  # Extract dsem map and parameter objects
  # - Create mapList object
  mapList <- sapply(fit_dsem$tmb_inputs$parameters, function(x) replace(x, values = c(1:length(x))))

  # - Copy dsem map-factor to map-list
  for(i in 1:length(fit_dsem$tmb_inputs$map)){
    parname <- names(fit_dsem$tmb_inputs$map)[i]
    mapList[[parname]] <- replace(fit_dsem$tmb_inputs$parameters[[parname]], values = as.numeric(fit_dsem$tmb_inputs$map[[i]]))
  }

  # Debug mode
  if(debug){
    mapList <- sapply(mapList, function(x) replace(x, values = rep(NA, length=length(x))))
  }

  # Recruitment variance.
  # Locate the beta_z index of each species' recruitment SD: the two-headed
  # self-loop on recdevs[sp] (recdevs[sp] <-> recdevs[sp]).
  sf       <- fit_dsem$sem_full
  sf_start <- suppressWarnings(as.numeric(sf$start))
  sf_par   <- suppressWarnings(as.numeric(sf$parameter))
  sf_dir   <- abs(suppressWarnings(as.numeric(sf$direction)))
  # v5.0 renamed this control element sigma_rec_prior -> sigma_rec (the schema
  # keeps the old spelling as an alias and upgrades it on entry to the pipeline),
  # so reading the old name yields NULL and rep_len() then fails. Accept either,
  # preferring the canonical one.
  sigma_rec_in <- data_list$sigma_rec
  if (is.null(sigma_rec_in)) sigma_rec_in <- data_list$sigma_rec_prior
  if (is.null(sigma_rec_in)) {
    stop("`data_list` carries neither `sigma_rec` nor `sigma_rec_prior`; the ",
         "DSEM needs a recruitment-SD starting value.", call. = FALSE)
  }
  sigma_rec_prior <- rep_len(sigma_rec_in, data_list$nspp)

  rec_sd_idx   <- integer(data_list$nspp)   # 1-based beta_z index; 0 if SD is fixed in the sem
  rec_sd_fixed <- numeric(data_list$nspp)   # fixed SD value used when rec_sd_idx == 0
  for(sp in seq_len(data_list$nspp)){
    nm   <- rec_cols[sp]
    rows <- which(sf$first == nm & sf$second == nm & sf_dir == 2)
    pn   <- unique(sf_par[rows]); pn <- pn[!is.na(pn) & pn > 0]
    if(length(pn) > 0){
      rec_sd_idx[sp] <- pn[1]
      # Initialize the recruitment SD at sigma_rec_prior; fix it if !random_rec
      fit_dsem$tmb_inputs$parameters$beta_z[rec_sd_idx[sp]] <- sigma_rec_prior[sp]
      if(!data_list$random_rec){
        mapList$beta_z[rec_sd_idx[sp]] <- NA
      }
    } else {
      # Recruitment SD fixed in the sem (NA parameter name): use its start value
      val <- sf_start[rows]; val <- val[!is.na(val)]
      rec_sd_fixed[sp] <- if(length(val) > 0) val[1] else sigma_rec_prior[sp]
    }
  }
  fit_dsem$tmb_inputs$data$rec_sd_idx   <- as.integer(rec_sd_idx)
  fit_dsem$tmb_inputs$data$rec_sd_fixed <- as.numeric(rec_sd_fixed)

  # Optional lognormal prior on the estimated recruitment SD, centered at the
  # assessment value (sigma_rec_prior). Regularizes R_sd away from the
  # 1/R_sd^2 collapse when env covariates over-explain recruitment. Enabled per
  # species only when the SD is estimated (rec_sd_idx >= 1) and a finite
  # sigmaR_prior_sd is supplied via build_DSEM(); off (NA) by default.
  prior_sd_in <- dsem_settings$sigmaR_prior_sd
  if(is.null(prior_sd_in)) prior_sd_in <- NA
  prior_sd <- rep_len(prior_sd_in, data_list$nspp)
  rec_sd_use_prior <- as.integer((rec_sd_idx >= 1) & is.finite(prior_sd))
  fit_dsem$tmb_inputs$data$rec_sd_prior     <- as.numeric(sigma_rec_prior)
  fit_dsem$tmb_inputs$data$rec_sd_prior_sd  <- as.numeric(ifelse(is.finite(prior_sd), prior_sd, 0))
  fit_dsem$tmb_inputs$data$rec_sd_use_prior <- as.integer(rec_sd_use_prior)
  # The recruitment-SD prior is applied per species, which is what a sem with a
  # sigmaR<sp> per species gives. A sem that SHARES one SD parameter across
  # species is different: the prior then lands on that one parameter once per
  # species (tightening it by sqrt(n)), and the per-species prior centres and
  # initial values collide. Warn rather than let either happen silently.
  .shared <- rec_sd_idx[rec_sd_idx >= 1]
  if (anyDuplicated(.shared)) {
    warning("Species share a recruitment-SD parameter in the sem. The ",
            "recruitment-SD prior is applied per species, so it is applied ",
            "once per species to that one parameter; and only one species' ",
            "`sigma_rec` sets its starting value. Give each species its own ",
            "sigmaR term if that is not what you intend.", call. = FALSE)
  }

  # x_tj column of each species' recdevs (0-based for the cpp).
  # Growth latents the sem actually names, decoded against the registry so the
  # names and the index map cannot drift. One row per named column: its x_tj
  # column (0-based), the species and sex it drives, and the 0-based slot in
  # growth_linkage_offset's last dimension. Empty for a recruitment-only sem,
  # which is what keeps those fits on the existing code path.
  g_reg <- .dsem_latent_registry(data_list)[names(.DSEM_GROWTH_TARGETS)]
  g_col <- unlist(lapply(g_reg, `[[`, "columns"), use.names = FALSE)
  g_hit <- match(g_col, colnames(dsem_data))
  keep  <- !is.na(g_hit)
  fit_dsem$tmb_inputs$data$growth_dev_col <- as.integer(g_hit[keep] - 1L)
  fit_dsem$tmb_inputs$data$growth_dev_sp  <-
    as.integer(unlist(lapply(g_reg, `[[`, "species"), use.names = FALSE)[keep] - 1L)
  fit_dsem$tmb_inputs$data$growth_dev_sex <-
    as.integer(unlist(lapply(g_reg, `[[`, "sex"), use.names = FALSE)[keep] - 1L)
  fit_dsem$tmb_inputs$data$growth_dev_par <-
    as.integer(unlist(lapply(g_reg, `[[`, "param"), use.names = FALSE)[keep])

  # Natural mortality latents the sem names. No param slot: M_linkage_offset is
  # (nspp, max_sex, max_age, nyrs), so the template needs species and sex and
  # broadcasts over ages.
  m_reg <- .dsem_latent_registry(data_list)[["natural_mortality"]]
  m_hit <- match(m_reg$columns, colnames(dsem_data))
  m_keep <- !is.na(m_hit)
  fit_dsem$tmb_inputs$data$mdev_col <- as.integer(m_hit[m_keep] - 1L)
  fit_dsem$tmb_inputs$data$mdev_sp  <- as.integer(m_reg$species[m_keep] - 1L)
  fit_dsem$tmb_inputs$data$mdev_sex <- as.integer(m_reg$sex[m_keep] - 1L)
  .check_dsem_support(data_list, fit_dsem)

  rec_dev_col <- match(rec_cols, colnames(dsem_data)) - 1L
  n_rec_named <- sum(!is.na(rec_dev_col))
  # A sem may drive growth and leave recruitment alone, so naming NO recruitment
  # column is allowed. Naming SOME but not all stays an error: that is a typo in
  # a species index, not a design, and it would silently delete the recruitment
  # process for the species whose column is missing while keeping it for the
  # others.
  if (n_rec_named > 0L && n_rec_named < length(rec_cols)) {
    stop("The sem names a recruitment latent column for ", n_rec_named, " of ",
         length(rec_cols), " species. Name one for every species (",
         paste(rec_cols, collapse = ", "), ") so the SEM drives recruitment, ",
         "or none of them so it does not. Driving a subset would leave the ",
         "remaining species with no recruitment deviations at all.",
         call. = FALSE)
  }
  # ... but a sem that drives NOTHING is a different thing, and before growth
  # became drivable the refusal above caught it. Such a model still fits a GMRF
  # to its covariates and still adds jnll_dsem to the objective, while no part
  # of the population model reads any of it -- so it converges, reports an ABC,
  # and its objective is not comparable to the same model without the SEM.
  # Every drivable process has to be counted here.
  if (n_rec_named == 0L && !any(keep) && !any(m_keep)) {
    stop("The sem names no latent process column, so it drives nothing in the ",
         "model. Name a recruitment column for every species (",
         paste(rec_cols, collapse = ", "), ") and/or one of this model's growth ",
         "columns (", paste(g_col, collapse = ", "),
         ") and/or a natural-mortality column (",
         paste(m_reg$columns, collapse = ", "),
         "). A sem over covariates alone still contributes its density to the ",
         "objective without driving any process.", call. = FALSE)
  }
  # -1, not 0: unread when drives_rec = 0, but several R consumers index x_tj as
  # rec_dev_col[sp] + 1L, and 0 would alias the FIRST latent column -- a growth
  # latent. -1 does NOT make such a read error: measured, index 0 silently
  # no-ops an assignment, returns a 0-column matrix, and inverts a mask built by
  # negative subscripting. What makes it loud is the `.col >= 1L` assertion at
  # each consumer; the sentinel's job is only to be a value no column can have,
  # so a missed assertion cannot quietly hit a real column. NA is not an option:
  # it reaches a DATA_IVECTOR as INT_MIN.
  if (n_rec_named == 0L) rec_dev_col <- rep(-1L, length(rec_cols))
  fit_dsem$tmb_inputs$data$rec_dev_col <- as.integer(rec_dev_col)

  # Whether this SEM owns the recruitment deviations, which the template reads
  # as dsem_drives_rec. Read off the columns the sem actually names, not off
  # rec_dev_col: match() returns NA rather than a shorter vector, so
  # length(rec_dev_col) is nspp whatever the sem says and would always test
  # TRUE.
  fit_dsem$tmb_inputs$data$drives_rec <-
    as.integer(any(rec_cols %in% colnames(dsem_data)))

  # estimate_projection (DSEM latent projection states) overlaps with
  # proj_mean_rec (SRR projection method): under mean-recruitment projection
  # the template sets R = avg_R past endyr, so the SEM's projection states
  # reach only the dynamic B0/BF reference series (which does set the ABC under
  # DynamicHCR = 1) and the states are otherwise fitted and then ignored.
  # fit_mod() resolves this by letting estimate_projection win, so this only
  # fires for a direct caller that built the objects itself.
  if( isTRUE(dsem_settings$estimate_projection) &&
      isTRUE(as.logical(data_list$proj_mean_rec)) ){
    warning("`estimate_projection = TRUE` with `proj_mean_rec = TRUE`: ",
            "projected recruitment is set to avg_R, so the SEM's projection ",
            "states do not drive it. Set `proj_mean_rec = FALSE` to project ",
            "recruitment through the SEM, or `estimate_projection = FALSE` to ",
            "leave the projection years out of the DSEM. fit_mod() applies the ",
            "first of those for you.", call. = FALSE)
  }

  # Number of DSEM time steps, and whether they cover the projection. The C++
  # call site needs this: with estimate_projection = FALSE the latent state
  # matrix is SHORTER than nyrs, so the recruitment-deviation copy must be
  # rec_dev.row(sp).head(nyrs_dsem), not a whole-row assignment.
  # Take the span from the BUILT data, not the year range: a duplicated year in
  # env_data makes them differ, and the C++ then copies latent states into the
  # wrong assessment years with no error.
  if (nrow(dsem_data) != length(data_list$styr:dsem_endyr)) {
    stop("env_data resolves to ", nrow(dsem_data), " rows over ",
         data_list$styr, ":", dsem_endyr, " (", length(data_list$styr:dsem_endyr),
         " years) -- check for duplicated or missing years.", call. = FALSE)
  }
  fit_dsem$tmb_inputs$data$nyrs_dsem <- as.integer(nrow(dsem_data))
  fit_dsem$covers_projection <- isTRUE(dsem_settings$estimate_projection)

  # Return
  fit_dsem$tmb_inputs$map <- sapply(mapList, function(x) factor(x))
  fit_dsem$mapList <- mapList
  fit_dsem$sem <- dsem_settings$sem
  # Carry the specification the objects were built from. fit_mod() records it on
  # data_list$dsem_settings, which is what .refit_like() forwards and what the
  # diagnostic guards look for -- so a fit handed pre-built objects would
  # otherwise look DSEM-free to both.
  fit_dsem$dsem_settings <- dsem_settings
  return(fit_dsem = .as_rceattle_dsem(fit_dsem))
}


#' Subclass a dsem object so Rceattle's methods take precedence
#'
#' Prepending keeps `inherits(x, "dsem")` TRUE for anything that tests it.
#' `oldClass()`, not `class()`: an unclassed list would otherwise gain a literal
#' `"list"` class, which `NextMethod()` then falls through to.
#'
#' @param x An object from [build_dsem_objects()].
#' @return `x` with `"Rceattle_dsem"` prepended to its class.
#' @noRd
.as_rceattle_dsem <- function(x) {
  if (inherits(x, "Rceattle_dsem")) return(x)
  class(x) <- c("Rceattle_dsem", oldClass(x))
  x
}


# -----------------------------------------------------------------------------
# DSEM <-> CEATTLE integration helpers
# -----------------------------------------------------------------------------
# build_dsem_objects() produces TMB inputs that fit_mod() merges into the
# CEATTLE parameter list, map, random-effects vector, and data list at four
# points in its pipeline. The contract of what DSEM contributes lives here so
# those merges stay consistent:
#   * parameters: beta_z, lnsigma_z, mu_j, delta0_j, x_tj
#   * map:        mapList + mapFactor entries for the above
#   * random:     x_tj (when random_rec = TRUE)
#   * data:       options, RAM, RAMstart, familycode_j, linkcode_j, sigmastart_j,
#                 eps_tj, y_tj, obs_idx, unobs_idx

#' Names of the DSEM parameters contributed to the CEATTLE parameter list
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
dsem_param_names <- function(dsem) {
  names(dsem$tmb_inputs$parameters)
}

#' Merge DSEM parameters into a CEATTLE parameter list
#'
#' Adds any DSEM parameters absent from \code{param_list}, pulling defaults from
#' the built DSEM objects. A no-op for parameters already present, so it is safe
#' for both a fresh \code{build_params()} list (none present -> all added) and a
#' prior fit's \code{estimated_params} (all present -> unchanged). Nothing is
#' overwritten or duplicated.
#'
#' @param param_list CEATTLE parameter list (from build_params or inits)
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
merge_dsem_params <- function(param_list, dsem) {
  dsem_par <- dsem$tmb_inputs$parameters
  missing_par <- setdiff(names(dsem_par), names(param_list))
  if (length(missing_par) > 0) {
    param_list <- c(param_list, dsem_par[missing_par])
  }
  param_list
}

#' Merge DSEM map entries into a CEATTLE map object
#'
#' DSEM contributes both \code{mapList} (named numeric/NA vectors) and
#' \code{mapFactor} (the TMB factor map). Merge them together so the two halves
#' never drift apart.
#'
#' @param map CEATTLE map object (from \code{build_map}) with \code{mapList} and
#'   \code{mapFactor}
#' @param dsem object returned by \code{\link{build_dsem_objects}}
#' @keywords internal
#' @noRd
merge_dsem_map <- function(map, dsem) {
  map$mapList   <- c(map$mapList, dsem$mapList)
  map$mapFactor <- c(map$mapFactor, dsem$tmb_inputs$map)
  map
}


# ---- Pre-fit DSEM spec screen ----
# Moved from R/0-convergence.R during the dev-data-workflow (v5.0) merge so it
# lives with the other DSEM builders. Depends on the .conv_record/.conv_overall
# convergence internals -- reconcile with v5.0 R/0-convergence.R in Tier 3.
#' Pre-fit screen of a DSEM recruitment specification
#'
#' Runs the Tier-1 (\code{"spec"}) convergence checks on a built DSEM and the
#' model's \code{env_data}, returning the same \code{"Rceattle_convergence"}
#' object as [convergence_diagnostics()]. Call it directly to screen a
#' specification before fitting; \code{fit_mod()} does **not** run it for you.
#'
#' Checks (tier \code{"spec"}):
#' \itemize{
#'   \item \code{rec_predictor_observability} -- recruitment predictors observed
#'     in a low fraction of hindcast years (free latents may absorb recruitment
#'     variance -> sigmaR collapse).
#'   \item \code{rec_design_conditioning} -- condition number / max pairwise
#'     correlation of the lag-aligned recruitment design matrix (collinear
#'     predictors).
#'   \item \code{growth_predictor_observability},
#'     \code{growth_design_conditioning} -- the same two checks for paths into a
#'     growth latent (\code{l1devs}, \code{linfdevs}, \code{richdevs}), reported
#'     separately so a growth predictor's coverage or collinearity is not
#'     attributed to recruitment.
#'   \item \code{mdev_predictor_observability},
#'     \code{mdev_design_conditioning} -- the same two for paths into a natural
#'     mortality latent (\code{mdevs}).
#'   \item \code{covariate_scale} -- covariates of ANY driven process spanning
#'     orders of magnitude in SD (ill-conditioning; consider standardizing).
#'     Pooled across processes deliberately: they share one latent field and one
#'     precision, so scale heterogeneity is a property of the joint field.
#'   \item \code{covariate_variance} -- an observed recruitment predictor with
#'     no exogenous variance of its own (no two-headed path, or one fixed at
#'     zero). Such a variable is deterministic: the model computes it rather
#'     than reading \code{env_data}, so the series supplied has no effect.
#' }
#'
#' @param data_list A cleaned data list (must carry \code{env_data},
#'   \code{styr}, \code{endyr}).
#' @param dsem A built DSEM object (from [build_dsem_objects()]) carrying
#'   \code{sem_full}.
#'
#' @return An object of class \code{"Rceattle_convergence"}.
#' @export
check_dsem_spec <- function(data_list, dsem) {
  empty <- structure(list(status = "OK", checks = list()),
                     class = "Rceattle_convergence")
  sf <- tryCatch(dsem$sem_full, error = function(e) NULL)
  ed <- data_list$env_data
  # A build_DSEM() SPECIFICATION carries `sem` but no `sem_full`, and passing one
  # here is the natural mistake -- the argument is named `dsem` and that is what
  # fit_mod(dsem = ) takes. Returning an empty pass would report OK on a
  # specification nothing has looked at. Say what to call instead.
  if (is.null(sf) && is.list(dsem) && "sem" %in% names(dsem)) {
    stop("`dsem` is a build_DSEM() specification, which carries no `sem_full`, ",
         "so there is nothing to check. Build it first:\n",
         "  check_dsem_spec(data_list, build_dsem_objects(dsem, data_list = data_list))",
         call. = FALSE)
  }
  if (is.null(sf)) return(empty)
  out <- list()
  finish <- function(out) structure(list(status = .conv_overall(out),
                                         checks = out),
                                    class = "Rceattle_convergence")

  # T1.4 -- a variable with no exogenous variance
  # dsem's projecting parameterizations compute such a variable from its
  # incoming paths instead of estimating it, so the env_data supplied for it is
  # never read (scaling that column leaves the objective bit-identical), it
  # drops out of the reported precision that sample_rec() draws from, and the
  # recruitment bias correction cannot condition on it.
  #
  # The rule is dsem's own: a lag-0 two-headed self-path whose parameter AND
  # start are both zero. Not "has no two-headed path" -- make_dsem_ram() adds an
  # estimated V[x] row for any variable that lacks one, so a missing line is the
  # ordinary way to ask for a freely estimated variance.
  #
  # Over every variable in the sem rather than the recruitment predictors alone,
  # because a variable two steps upstream of recruitment does the same thing.
  # Runs before the env_data checks below, which it does not need.
  sf_dir <- abs(suppressWarnings(as.numeric(sf$direction)))
  sf_par <- suppressWarnings(as.numeric(sf$parameter))
  sf_st  <- suppressWarnings(as.numeric(sf$start))
  sf_lag <- suppressWarnings(as.numeric(sf$lag))
  all_vars <- unique(c(as.character(sf$first), as.character(sf$second)))
  novar <- all_vars[vapply(all_vars, function(p) {
    rows <- which(sf$first == p & sf$second == p & sf_dir == 2 &
                  (is.na(sf_lag) | sf_lag == 0))
    if (length(rows) == 0) return(FALSE)
    all((is.na(sf_par[rows]) | sf_par[rows] == 0) &
        (is.na(sf_st[rows])  | sf_st[rows]  == 0))
  }, logical(1))]
  if (length(novar) > 0) {
    out$covariate_variance <- .conv_record(
      "covariate_variance", "spec", "WARN",
      sprintf(paste0(
        "Variable(s) with no exogenous variance: %s. Each is deterministic, so ",
        "the env_data supplied for it is not read and it drops out of the ",
        "precision sample_rec() draws from. Give it a non-zero variance (e.g. ",
        "`%s <-> %s, 0, sd%s, 1`)."),
        paste(novar, collapse = ", "), novar[1], novar[1], novar[1]),
      data.frame(variable = novar, row.names = NULL))
  }

  # Everything below reads env_data.
  if (is.null(ed) || is.null(ed$Year)) return(finish(out))

  # Paths into a latent node, one-headed (direction 1), SPLIT BY THE PROCESS the
  # target drives, so each process is reported under its own records and never
  # shares a design matrix with another. The prefix is enough to classify -- no
  # data_list needed, which keeps this function's contract.
  tgt <- as.character(sf$second)
  # NA-safe: `direction` is user-supplied on a hand-built sem_full, and
  # `any(NA)` would make the guard below throw.
  # Every neighbouring read of this column is equally defensive.
  dir1 <- suppressWarnings(as.numeric(sf$direction))
  dir1 <- !is.na(dir1) & dir1 == 1
  into_latent <- grepl(.dsem_latent_pattern(), tgt) & dir1
  # The suffix is already validated by the pattern above, so the stem alone
  # classifies. startsWith() rather than a fourth copy of `[0-9]+[fm]?$`.
  is_rec <- into_latent & startsWith(tgt, .DSEM_LATENT_PREFIXES[["recruitment"]])
  is_gro <- into_latent &
    Reduce(`|`, lapply(.DSEM_GROWTH_TARGETS, function(st) startsWith(tgt, st)))
  is_m <- into_latent &
    startsWith(tgt, .DSEM_LATENT_PREFIXES[["natural_mortality"]])
  # `into_latent` is the UNION of every registered prefix while the three groups
  # above enumerate them one at a time, so a process added to the registry would
  # otherwise be checked by none of them and say nothing. Loud instead.
  if (any(into_latent & !(is_rec | is_gro | is_m))) {
    stop("No specification check covers paths into ",
         paste(unique(tgt[into_latent & !(is_rec | is_gro | is_m)]), collapse = ", "),
         ". Add the process to check_dsem_spec()'s classification.", call. = FALSE)
  }

  if (!any(into_latent)) return(finish(out))

  styr <- data_list$styr; endyr <- data_list$endyr
  yrs  <- styr:endyr
  nyr  <- length(yrs)
  edh  <- ed[ed$Year %in% yrs, , drop = FALSE]

  # Predictors of ANY driven process, for the scale check below, which is about
  # the covariates rather than about one process.
  preds <- intersect(unique(as.character(sf$first[into_latent])), colnames(ed))
  if (length(preds) == 0) return(finish(out))

  # T1.1 / T1.2 per process. Two callers, so this earns its place as a helper;
  # the recruitment wording and record names are unchanged from before the
  # split, because a script may read them.
  design_checks <- function(paths, noun, sd_name, obs_name, cond_name) {
    res <- list()
    pr <- intersect(unique(as.character(paths$first)), colnames(ed))
    if (!nrow(paths) || !length(pr)) return(res)

    # T1.1 -- latent observability of this process's predictors
    cov_frac <- vapply(pr, function(p) sum(!is.na(edh[[p]])) / nyr, numeric(1))
    low <- pr[cov_frac < 0.5]
    if (length(low) > 0) {
      sev <- if (any(cov_frac[low] < 0.25)) "WARN" else "NOTE"
      res[[obs_name]] <- .conv_record(
        obs_name, "spec", sev,
        sprintf(paste0(
          "%s predictor(s) observed in <50%% of hindcast years: %s. ",
          "Sparse latents can absorb %s variance and drive %s to 0."),
          noun, paste(sprintf("%s (%.0f%%)", low, 100 * cov_frac[low]),
                      collapse = ", "), tolower(noun), sd_name),
        data.frame(predictor = pr, coverage = round(cov_frac, 3),
                   row.names = NULL))
    }

    # T1.2 -- lag-aligned design conditioning, PER LATENT EQUATION.
    # Collinearity is a property of one target's design matrix: each path
    # coefficient is identified by its own target latent's residual, so
    # correlation ACROSS equations says nothing about either. Pooling targets
    # makes this fire by construction -- one predictor driving two latents
    # appears as two identical columns, which make.unique() renames and cor()
    # scores at 1.00, reporting a condition number of ~6e16 for a design with no
    # collinearity at all. That is the canonical growth sem (one covariate, L1
    # and Linf) and the canonical multispecies recruitment sem (one covariate,
    # recdevs1 and recdevs2), so the pooled form was wrong far more often than
    # it was right. Worst severity across equations wins, and the record names
    # the equation.
    worst <- NULL
    for (tg in unique(as.character(paths$second))) {
      pe <- paths[as.character(paths$second) == tg, , drop = FALSE]
      cols <- lapply(seq_len(nrow(pe)), function(i) {
        p <- as.character(pe$first[i]); lag <- suppressWarnings(as.numeric(pe$lag[i]))
        if (!p %in% colnames(ed) || is.na(lag)) return(NULL)
        ed[[p]][match(yrs - lag, ed$Year)]   # predictor at the lag the SEM uses
      })
      names(cols) <- make.unique(paste0(pe$first, "_lag", pe$lag))
      cols <- cols[!vapply(cols, is.null, logical(1))]
      if (length(cols) < 2) next
      X  <- do.call(cbind, cols)
      Xc <- X[stats::complete.cases(X), , drop = FALSE]
      if (nrow(Xc) <= ncol(Xc)) next
      keep <- apply(Xc, 2, function(z) stats::sd(z) > 0)
      Xc   <- Xc[, keep, drop = FALSE]
      if (ncol(Xc) < 2) next
      R    <- stats::cor(Xc)
      kap  <- tryCatch(kappa(R, exact = TRUE), error = function(e) NA_real_)
      maxr <- max(abs(R[upper.tri(R)]))
      sev  <- if ((is.finite(kap) && kap > 100) || maxr > 0.9) "WARN"
              else if (maxr > 0.8) "NOTE" else "OK"
      if (sev == "OK") next
      rank <- c(NOTE = 1L, WARN = 2L)
      if (is.null(worst) || rank[[sev]] > rank[[worst$sev]])
        worst <- list(sev = sev, tg = tg, kap = kap, maxr = maxr, n = nrow(Xc),
                      R = R)
    }
    if (!is.null(worst)) {
      res[[cond_name]] <- .conv_record(
        cond_name, "spec", worst$sev,
        sprintf(paste0("Collinear lag-aligned %s design for %s ",
          "(condition number = %.3g, max |r| = %.2f over %d years)."),
          tolower(noun), worst$tg, worst$kap, worst$maxr, worst$n),
        list(target = worst$tg, condition_number = worst$kap,
             max_abs_cor = worst$maxr, R = worst$R))
    }
    res
  }

  # "sigmaR" is the recruitment wording from before the split, kept byte-for-byte
  # so a script matching on the message still matches.
  out <- c(out, design_checks(sf[is_rec, , drop = FALSE], "Recruitment",
                              "sigmaR",
                              "rec_predictor_observability",
                              "rec_design_conditioning"))
  out <- c(out, design_checks(sf[is_gro, , drop = FALSE], "Growth",
                              "the growth latent's SD",
                              "growth_predictor_observability",
                              "growth_design_conditioning"))
  out <- c(out, design_checks(sf[is_m, , drop = FALSE], "Natural mortality",
                              "the M latent's SD",
                              "mdev_predictor_observability",
                              "mdev_design_conditioning"))

  # T1.5 -- covariate scale heterogeneity
  sds <- vapply(preds, function(p) stats::sd(ed[[p]], na.rm = TRUE), numeric(1))
  sds <- sds[is.finite(sds) & sds > 0]
  if (length(sds) >= 2) {
    ratio <- max(sds) / min(sds)
    if (ratio > 20) {
      out$covariate_scale <- .conv_record(
        "covariate_scale", "spec", "NOTE",
        sprintf(paste0("DSEM covariates span %.0fx in SD (%s: %.3g; ",
          "%s: %.3g); consider standardizing env_data."),
          ratio, names(which.max(sds)), max(sds),
          names(which.min(sds)), min(sds)),
        data.frame(predictor = names(sds), sd = signif(sds, 3),
                   row.names = NULL))
    }
  }

  finish(out)
}

# Draw the DSEM's UNKNOWN projection states from the precision that scores them.
#
# This is the conditional draw of the unknown projection nodes given everything
# the fit already knows, taken from the GMRF the model was fitted with -- the
# mechanism dsem's own simulate() uses, rmvnorm_prec(xhat + delta, Q).
# Partitioning the field into drawn (P) and known (H):
#
#   x_P | x_H  ~  N( mu_P - Q_PP^-1 Q_PH (x_H - mu_H),  Q_PP^-1 )
#
# so a projected recruitment deviation inherits whatever the sem says: a
# self-path makes it persist, covariate paths make it respond to the
# environment, and a cross-species path couples the stocks. Drawing iid instead
# would project a different process from the one estimated.
#
# WHAT GOES IN P IS THE WHOLE POINT. Not every projection node is unknown. Under
# family = "fixed" a covariate column of x_tj IS the environmental data, pinned
# by the map; where env_data extends past endyr those projection nodes are KNOWN
# FUTURE VALUES -- a climate scenario the projection is meant to respond to.
# Drawing them instead of conditioning on them replaces the scenario with noise
# and returns the no-covariate answer, which is the opposite of what a
# climate-linked assessment is asking for. It also inflates the draw variance,
# against a bias correction the template computes CONDITIONAL on the
# environment (that is what dsem_cond_k is for in ceattle.cpp section 5.5b, and
# it is built from the same per-cell rule). So P is read off the MAP: draw what
# is estimated, condition on what is fixed, wherever it sits.
#
# Requires the latent field to SPAN the projection -- build_DSEM(
# estimate_projection = TRUE). With estimate_projection = FALSE the GMRF stops
# at endyr, the likelihood says nothing about projection years, and any draw
# there is an extrapolation rather than a draw from the model.
#' @noRd
.dsem_draw_projection <- function(fit, sample = TRUE) {
  q <- fit$quantities
  Q <- q$dsem_Q
  if (is.null(Q) || !length(Q)) {
    stop("This model reports no DSEM precision, so its projection states ",
         "cannot be drawn from the fitted process. Either it carries no DSEM, ",
         "or it predates the reporting of `dsem_Q` and must be refitted.",
         call. = FALSE)
  }
  X   <- as.matrix(fit$estimated_params$dsem_x_tj)
  mu  <- as.matrix(q$dsem_xhat_tj) + as.matrix(q$dsem_delta_tj)
  n_t <- nrow(X); n_j <- ncol(X)
  if (nrow(Q) != n_t * n_j) {
    # Under the projecting parameterizations the template reports the precision
    # of the OBSERVED block only, indexed by position within obs_idx rather than
    # by the stacked cell k, so it is short whenever the sem carries a
    # deterministic variable. Name that as the usual cause without claiming it
    # is the only one -- a mis-wired precision reaches here too.
    stop("The reported DSEM precision is ", nrow(Q), "x", ncol(Q),
         " but the latent field is ", n_t, "x", n_j,
         "; they cannot be aligned. Usually this means a variable in the sem ",
         "has no exogenous variance -- its two-headed path is pinned at zero -- ",
         "so the model computes that variable instead of estimating it. Give it ",
         "a non-zero variance (e.g. `BT <-> BT, 0, sdBT, 1`); ",
         "check_dsem_spec() names the variable.", call. = FALSE)
  }

  n_hind <- fit$data_list$endyr - fit$data_list$styr + 1L
  if (n_t <= n_hind) {
    stop("The DSEM's latent states stop at endyr, so the model says nothing ",
         "about the projection years and they cannot be drawn from it. ",
         "Rebuild with build_DSEM(estimate_projection = TRUE).", call. = FALSE)
  }

  # k = j * n_t + t in the template's 0-based indexing, which is column-major
  # over x_tj -- so as.numeric(X)[k] is the node the GMRF scores at k.
  k_of <- function(t, j) (j - 1L) * n_t + t
  proj <- as.integer(unlist(lapply(seq_len(n_j), function(j)
    k_of((n_hind + 1L):n_t, j))))

  # NA in the map means the node is FIXED -- a known value to condition on, not
  # something to draw. An absent map means nothing is fixed.
  mp <- fit$map$mapList$dsem_x_tj
  estimated <- if (is.null(mp)) rep(TRUE, n_t * n_j) else !is.na(as.numeric(mp))
  if (length(estimated) != n_t * n_j) {
    stop("The DSEM's parameter map covers ", length(estimated), " states but ",
         "the latent field has ", n_t * n_j, "; they cannot be aligned.",
         call. = FALSE)
  }

  P <- proj[estimated[proj]]
  if (!length(P)) {
    stop("Every projection state of this DSEM is fixed, so there is nothing ",
         "to draw. The recruitment deviations must be estimated over the ",
         "projection for the process to supply them.", call. = FALSE)
  }
  H <- setdiff(seq_len(n_t * n_j), P)

  Qpp <- Q[P, P, drop = FALSE]
  Qph <- Q[P, H, drop = FALSE]
  dev_h <- as.numeric(X)[H] - as.numeric(mu)[H]
  cond_mu <- as.numeric(mu)[P] - solve(Qpp, Qph %*% dev_h)

  out <- X
  if (isTRUE(sample)) {
    # Draw N(0, Qpp^-1) from the precision, as dsem::simulate() does.
    L <- chol(Qpp)              # Qpp = L'L, so backsolve(L, z) ~ N(0, Qpp^-1)
    out[P] <- as.numeric(cond_mu) + backsolve(L, stats::rnorm(length(P)))
  } else {
    # The conditional MEAN -- the deterministic counterpart of the draw, and
    # what the process expects given everything already known.
    out[P] <- as.numeric(cond_mu)
  }
  out
}


#' Does this DSEM want its standard-deviation paths bounded at 0?
#'
#' Reads `build_DSEM(bound_sd = )` off whichever object the caller has -- a
#' built DSEM carries the settings it was built from in `$dsem_settings`, and a
#' bare specification is the settings list itself.
#'
#' Defaults to TRUE when the field is absent, which is what a model fitted
#' before the argument existed carries. That matches the argument's own default,
#' so an old fit and a new one are bounded the same way.
#'
#' @param dsem A built DSEM or a `build_DSEM()` specification.
#' @return `TRUE` or `FALSE`.
#' @noRd
.dsem_bound_sd <- function(dsem) {
  v <- dsem$dsem_settings$bound_sd
  if (is.null(v)) v <- dsem$bound_sd
  if (is.null(v)) return(TRUE)
  isTRUE(v)
}


#' Which `dsem_beta_z` entries are standard deviations, not path coefficients
#'
#' `dsem_beta_z` holds every free SEM path in one vector: one-headed paths
#' (`->`) are regression coefficients whose sign is meaningful, two-headed paths
#' (`<->`) are entries of the Cholesky factor of the exogenous covariance.
#'
#' A lag-0 two-headed SELF path is a diagonal of that factor -- a standard
#' deviation. Its sign is not identified: the likelihood sees only
#' \eqn{\Gamma'\Gamma}, and the template reads it as
#' `R_sd = sqrt(square(beta_z))`. So the surface is exactly symmetric about 0.
#' For an MLE that is a harmless pair of mirrored optima; for MCMC it is fatal --
#' the posterior is bimodal by construction, chains sit in whichever mode they
#' started in or wander between them, and R-hat and any interval on sigma are
#' meaningless. Constraining the diagonal to be non-negative is the standard
#' identifying restriction for a Cholesky factor and discards nothing.
#'
#' Only diagonals that stand ALONE in their row qualify. If the variable also
#' carries a cross-covariance (`A <-> B`) or a lagged two-headed path, its row
#' of \eqn{\Gamma} has off-diagonal entries, the sign flip is no longer
#' invariant term by term, and identifying it means flipping the whole row --
#' which this does not attempt. Those entries keep their unbounded support, and
#' `build_bounds()` says so.
#'
#' An index is returned only when EVERY sem row sharing it qualifies: two paths
#' constrained to one parameter by name share a `beta_z` entry.
#'
#' @param dsem A built DSEM (`mod_objects$dsem`), carrying `sem_full`.
#' @return A list of `sd` (1-based indices into `dsem_beta_z` that are
#'   standalone SDs) and `entangled` (names of variables whose SD could not be
#'   constrained), or NULL if the fit carries no usable `sem_full`.
#' @noRd
.dsem_sd_indices <- function(dsem) {
  sf <- dsem$sem_full
  if (is.null(sf) || !nrow(sf)) return(NULL)
  need <- c("direction", "first", "second", "lag", "parameter")
  if (!all(need %in% names(sf))) return(NULL)

  two <- sf$direction == 2
  # A row of Gamma is clean only if this variable carries no other two-headed
  # path -- no cross-covariance, no lagged variance.
  messy <- two & (sf$first != sf$second | sf$lag != 0)
  entangled <- unique(c(sf$first[messy], sf$second[messy]))

  qualifies <- two & sf$first == sf$second & sf$lag == 0 &
    !(sf$first %in% entangled)
  idx <- unique(sf$parameter[qualifies & sf$parameter > 0])
  # Shared by a path that does not qualify -> not an SD on its own terms.
  shared <- unique(sf$parameter[!qualifies & sf$parameter > 0])
  idx <- setdiff(idx, shared)

  blocked <- unique(sf$first[two & sf$first == sf$second & sf$lag == 0 &
                               sf$first %in% entangled])
  list(sd = sort(as.integer(idx)), entangled = blocked)
}


#' Draw a DSEM's whole latent field from the density that scored it
#'
#' @description
#' The process-error draw for a DSEM. Where the standard recruitment process is
#' redrawn one year at a time from `N(-bias*sigma^2/2, sigma)`, a DSEM's
#' deviations are the latent states of a GMRF -- a self-path makes them
#' autocorrelated, covariate paths make them respond to the environment, and a
#' cross-species path couples the stocks -- so the whole field is drawn at once
#' from the same precision the density uses, exactly as `dsem::simulate()` does.
#'
#' Every estimated cell OF A REQUESTED PROCESS is redrawn, conditional both on
#' the cells the map pins and on the estimated cells of processes not requested
#' -- those are held at their fitted values, which is what lets a caller redraw
#' growth without redrawing recruitment. See `state`. The
#' map is where the model says which is which. Under `family = "fixed"` a
#' covariate column *is* the `env_data` series and is pinned, so it is held --
#' redrawing it would simulate a different environment rather than a different
#' recruitment process. Under any other family the covariate is a latent state
#' observed with error: the state belongs to the process and is drawn here,
#' while its OBSERVATION is redrawn by the template and written back into
#' `env_data` by `sim_mod()`.
#'
#' The caller writes the result back into `dsem_x_tj`; the template then derives
#' `rec_dev` from it, so the lognormal bias correction is applied once, where it
#' always is, rather than reproduced here.
#'
#' Under `build_DSEM(estimate_projection = TRUE)` the field spans the projection
#' and those years are drawn too. That is wider than the standard path, which
#' stops at `endyr` -- but it is the same rule, draw what the density scores,
#' because with a projected field the GMRF *is* the density over those years.
#' With the default `estimate_projection = FALSE` the field stops at `endyr` and
#' the two agree.
#'
#' @param fit a fitted `Rceattle` model carrying a DSEM.
#' @param state `.sim_state_codes()`'s five-slot request vector. Slot 1 selects
#'   the recruitment latent columns, slot 2 the natural-mortality ones and slot 3
#'   the growth ones; a column left out is held at its fitted value and
#'   CONDITIONED on, not ignored, so asking for one process does not redraw
#'   another.
#' @return A matrix the shape of `dsem_x_tj`, carrying a logical matrix of the
#'   same shape in `attr(, "drawn")` marking the cells this call redrew. The
#'   mask is written by the draw itself so a caller cannot re-derive it wrongly;
#'   `sim_mod()` maps it onto `rec_dev` for `attr(x, "process_sim")`.
#' @noRd
.dsem_draw_process <- function(fit, state = rep(1L, 5L)) {
  q  <- fit$quantities
  Q  <- q$dsem_Q
  if (is.null(Q) || !length(Q)) {
    stop("This model reports no DSEM precision, so its latent process ",
         "cannot be redrawn from the density that scored it. Refit with the ",
         "current version.", call. = FALSE)
  }
  Q  <- as.matrix(Q)
  X  <- as.matrix(fit$estimated_params$dsem_x_tj)
  mu <- as.matrix(q$dsem_xhat_tj) + as.matrix(q$dsem_delta_tj)
  n_k <- nrow(X) * ncol(X)
  if (nrow(Q) != n_k) {
    stop("The reported DSEM precision is ", nrow(Q), "x", ncol(Q),
         " but the latent field is ", nrow(X), "x", ncol(X),
         "; they cannot be aligned. Usually this means a variable in the sem ",
         "has no exogenous variance -- its two-headed path is pinned at zero -- ",
         "so the model computes that variable instead of estimating it. ",
         "check_dsem_spec() names the variable.", call. = FALSE)
  }

  dsd <- fit$dsem$tmb_inputs$data
  # Columns by role, 1-based. A SEM that does not drive recruitment carries the
  # -1 sentinel, which must map to NO column: `(0 - 1) * nrow + seq_len(nrow)`
  # is a NEGATIVE subscript, and assigning TRUE through it inverts the mask
  # rather than erroring.
  rc <- dsd$rec_dev_col + 1L                       # 0-based in the C++
  rc <- rc[!is.na(rc) & rc >= 1L]
  gc_ <- unique(dsd$growth_dev_col + 1L)
  gc_ <- gc_[!is.na(gc_) & gc_ >= 1L]
  mc <- unique(dsd[["mdev_col"]] + 1L)
  mc <- mc[!is.na(mc) & mc >= 1L]
  # Every drivable process has to be counted here: a process missing from this
  # list falls through to the latent-covariate branch below, which draws
  # unconditionally.
  if (!length(rc) && !length(gc_) && !length(mc)) {
    stop("This fit records no latent process columns, so there is nothing to ",
         "redraw. Refit before simulating.", call. = FALSE)
  }

  cells <- function(j) if (!length(j)) integer(0) else
    unlist(lapply(j, function(k) (k - 1L) * nrow(X) + seq_len(nrow(X))))

  # Which cells this call MAY draw. Drawing is joint over the GMRF, so a cell
  # left out is not ignored -- it is CONDITIONED on, held at its fitted value,
  # which is what lets a caller ask for one process without redrawing another.
  # Latent covariate columns belong to no process slot and are always drawn:
  # under family = "fixed" a covariate IS the env_data series and the map pins
  # it (so `est` excludes it), while under any other family it is a latent state
  # observed with error, so holding it at a FITTED value would be neither data
  # nor a draw.
  want <- rep(FALSE, n_k)
  if (isTRUE(state[1] == 1L)) want[cells(rc)]  <- TRUE
  if (isTRUE(state[2] == 1L)) want[cells(mc)]  <- TRUE
  if (isTRUE(state[3] == 1L)) want[cells(gc_)] <- TRUE
  want[setdiff(seq_len(n_k), cells(c(rc, gc_, mc)))] <- TRUE

  mp  <- fit$map$mapList$dsem_x_tj
  est <- if (is.null(mp)) rep(TRUE, n_k) else !is.na(as.numeric(mp))
  E <- which(est & want); F_ <- setdiff(seq_len(n_k), E)
  # The guard has to test the REQUESTED PROCESS's cells, not E. E always
  # includes the latent covariate columns, so `!length(E)` cannot fire while any
  # covariate is free -- the draw would then "succeed" having redrawn no process
  # state, hand back a non-empty drawn mask, and sim_mod() would attach a
  # process_sim truth for a process that was never redrawn. compare_sim() then
  # refuses to compare against the operating model, pointing at a truth the draw
  # did not produce.
  proc_want <- rep(FALSE, n_k)
  if (isTRUE(state[1] == 1L)) proc_want[cells(rc)]  <- TRUE
  if (isTRUE(state[2] == 1L)) proc_want[cells(mc)]  <- TRUE
  if (isTRUE(state[3] == 1L)) proc_want[cells(gc_)] <- TRUE
  if (any(proc_want) && !any(est & proc_want)) {
    stop("Every latent state of the process(es) asked for is fixed in this ",
         "model, so there is nothing to redraw.", call. = FALSE)
  }
  if (!length(E)) {
    stop("Every latent state this draw covers is fixed in this model, so there ",
         "is no process to redraw.", call. = FALSE)
  }

  # A latent node with no exogenous variance is deterministic given the others,
  # so its row of the precision is degenerate and both the conditional mean and
  # the Cholesky fail on it. Say which node rather than pass LAPACK's "leading
  # minor is not positive definite" up to an assessment author.
  Qee <- Q[E, E, drop = FALSE]
  degenerate <- function(e) {
    stop("The DSEM precision is not positive definite over the states this ",
         "draw covers, so the process cannot be redrawn (", conditionMessage(e),
         "). This happens when a latent variable carries no exogenous variance ",
         "-- its two-headed path is fixed at zero -- so the model computes it ",
         "rather than estimating it. check_dsem_spec() names the variable.",
         call. = FALSE)
  }

  m <- as.numeric(mu)[E]
  if (length(F_)) {
    m <- m - tryCatch(solve(Qee, Q[E, F_, drop = FALSE] %*%
                              (as.numeric(X)[F_] - as.numeric(mu)[F_])),
                      error = degenerate)
  }
  # Qee = U'U from chol(), so backsolve(U, z) ~ N(0, Qee^-1).
  U <- tryCatch(chol(Qee), error = degenerate)
  out <- X
  out[E] <- as.numeric(m) + backsolve(U, stats::rnorm(length(E)))
  drawn <- array(FALSE, dim = dim(X))
  drawn[E] <- TRUE
  attr(out, "drawn") <- drawn
  out
}


#' Carry a DSEM process draw across to the shapes `sim_mod()` reports
#'
#' The template REPORTs `rec_dev_drawn_sim` from inside its own IID draw, which
#' is gated off under a DSEM, so it stays zero however much of the field was
#' redrawn. This maps what `.dsem_draw_process()` actually drew onto the two
#' shapes a caller needs: the latent field itself, which is the process a DSEM
#' self test should be scored against, and a `rec_dev`-shaped mask, since
#' `rec_dev` is what every recovery statistic already reads.
#'
#' Row `t` of the latent field is model year `t` counting from `styr`, so it is
#' column `t` of `rec_dev`; the species' column comes from the `rec_dev_col`
#' registry the fit carries. `rec_dev` runs through the projection while the
#' field may stop at `endyr`, so the mask is padded with `FALSE`.
#'
#' @param fit a fitted `Rceattle` model carrying a DSEM.
#' @param drawn the return of `.dsem_draw_process()`, carrying its `"drawn"`
#'   attribute.
#' @param rec_dim `dim()` of the template's `rec_dev`, to pad the mask to.
#' @return A list of `x` (the drawn field), `x_drawn` (its mask) and `rec` (the
#'   `rec_dev`-shaped mask), or NULL if the draw carried no mask.
#' @noRd
.dsem_rec_truth <- function(fit, drawn, rec_dim = NULL) {
  mask <- attr(drawn, "drawn")
  if (is.null(mask)) return(NULL)
  # Drop the -1 sentinel a SEM that does not drive recruitment carries: rc = 0
  # would index zero columns and yield an empty `rec` mask rather than saying so.
  rc <- fit$dsem$tmb_inputs$data$rec_dev_col + 1L      # 0-based in the C++
  rc <- rc[!is.na(rc) & rc >= 1L]
  rec <- if (!length(rc)) NULL else t(mask[, rc, drop = FALSE])  # species x yr
  if (!is.null(rec) && !is.null(rec_dim)) {
    full <- array(FALSE, dim = rec_dim)
    nyr <- min(ncol(rec), rec_dim[2])
    full[seq_len(nrow(rec)), seq_len(nyr)] <- rec[, seq_len(nyr), drop = FALSE]
    rec <- full
  }
  list(x = matrix(as.numeric(drawn), nrow(mask), ncol(mask)),
       x_drawn = mask, rec = rec)
}
