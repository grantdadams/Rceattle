# Rejection budget for a non-positive natural-scale index draw, as a multiple of
# the fleet's row count. Correlated fleets reject the whole vector whenever any
# single row is non-positive, so the joint rejection probability climbs with the
# number of rows; scaling the budget the same way on both branches stops a wide
# fleet from exhausting it while every row on its own is fine.
.SIM_INDEX_MAX_TRIES <- 100L

# Per-row rejection rate above which truncation is doing enough of the work that
# the simulated data no longer follow the likelihood's own normal. Compared
# against the WORST row, not the fleet mean -- truncation bites one marginal row
# at a time, and a fleet average hides it. Set well below the rate that would
# matter: rejecting even a twentieth of a row's draws already shifts that row's
# mean by several percent.
.SIM_INDEX_WARN_FRAC <- 0.02


#' Resolve `Index_distribution` to its integer family code
#'
#' Accepts either spelling a `fleet_control` column can hold -- the name
#' (`"MVN"`) or the code (`1`) -- because `index_ll_type` is built inside
#' `rearrange_data()` and is not carried on a fitted model's `data_list`.
#' Anything unrecognized falls back to lognormal, matching the column default.
#'
#' @param x The `Index_distribution` column.
#' @return Integer vector of codes, one per fleet. See `index_distribution_map`.
#' @noRd
.index_family_codes <- function(x) {
  if (is.null(x)) return(integer(0))
  chr <- trimws(as.character(x))
  num <- suppressWarnings(as.integer(chr))
  nmd <- as.integer(index_distribution_map[chr])
  out <- ifelse(!is.na(num), num, nmd)
  out[is.na(out)] <- 0L
  as.integer(out)
}


#' Warn where a random linkage the caller asked for was left at its fitted values
#'
#' Observed AR1 (Rogers QAR1) is the one random linkage the template leaves
#' alone: a new latent state paired with the original covariate observations
#' would describe two different histories. Say so, rather than let a process the
#' caller asked for come back unchanged without explanation.
#'
#' @param obj The simulated TMB object.
#' @param state Integer switch vector from `.sim_state_codes()`.
#' @noRd
.sim_warn_linkage_qar1 <- function(obj, state) {
  d <- obj$env$data
  obs <- d$linkage_re_obs
  if (is.null(obs) || !length(obs) || all(obs < 0)) return(invisible(FALSE))

  procs <- c("recruitment", "M", "growth", "catchability", "selectivity")
  hit <- character(0)
  for (i in seq_along(d$linkage_re_index)) {
    slot <- d$linkage_re_index[i]
    if (slot < 0) next
    grp <- d$linkage_re_sigma[slot + 1L]
    if (obs[grp + 1L] < 0) next
    # Composition is a linkage process (code 5) with no simulate_state slot, so
    # guard the index the way the template does rather than compare against NA.
    p <- d$linkage_process[i]
    if (is.na(p) || p < 0 || p >= length(state)) next
    if (state[p + 1L] == 1L) hit <- c(hit, procs[p + 1L])
  }
  if (length(hit)) {
    warning("Random linkage(s) on ", paste(unique(hit), collapse = ", "),
            " are observed AR1 (Rogers QAR1), whose latent state is measured by ",
            "an observed covariate series. Those were NOT redrawn and keep their ",
            "fitted values, because a new latent state paired with the original ",
            "covariate observations would describe two different histories. ",
            "Every other random linkage on those processes was drawn normally.",
            call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}


# The legacy time-varying catchability and selectivity deviations
# (index_q_dev, log_sel_slp_dev, sel_inf_dev, sel_coff_dev) carry densities in
# the template but have no SIMULATE block yet -- only the linkage-grammar
# versions are drawn. Handing back the fitted values without a word is how a
# switch that does nothing goes unnoticed, which is exactly what growth_re did.
.sim_warn_process_unsimulated <- function(data_list, state) {
  fc <- data_list$fleet_control
  if (is.null(fc)) return(invisible(FALSE))
  code <- function(x, map) {
    if (is.null(x)) return(integer(0))
    if (is.character(x)) unname(map[x]) else suppressWarnings(as.integer(x))
  }
  hit <- character(0)

  # Block (3) is a time block with no penalty -- fixed structure, not process
  # error -- so it is not something a draw would ever redraw.
  if (state[4] == 1L) {
    tvq <- code(fc$Time_varying_q, tv_q_map)
    qq  <- code(fc$Catchability, q_map)
    # Catchability 5/6 read Time_varying_q as an environmental index instead;
    # 6 (Rogers AR1) is itself a random process on the deviations.
    stochastic <- (tvq %in% c(1L, 2L, 4L) & !(qq %in% c(5L, 6L))) | qq %in% 6L
    if (any(stochastic, na.rm = TRUE)) hit <- c(hit, "catchability")
  }
  if (state[5] == 1L) {
    tvs <- code(fc$Time_varying_sel, tv_sel_map)
    # The AR1 selectivity FORMS (2DAR1 = 6, 3DAR1 = 7) carry sel_coff_dev random
    # effects of their own and ignore Time_varying_sel entirely, so keying only
    # off that column would leave them frozen without a word -- the loudest case,
    # since an AR1 selectivity field is exactly what someone redrawing
    # selectivity means to test.
    sel <- code(fc$Selectivity, sel_map)
    if (any(tvs %in% c(1L, 2L, 4L, 5L), na.rm = TRUE) ||
        any(sel %in% c(6L, 7L), na.rm = TRUE)) hit <- c(hit, "selectivity")
  }

  if (length(hit)) {
    warning("Process error was requested for ", paste(hit, collapse = " and "),
            ", but this model expresses it with the switch-based deviations ",
            "(Time_varying_q / Time_varying_sel, or an AR1 selectivity form), ",
            "which are not simulated yet -- only the equivalent random linkages ",
            "are. Those deviations keep their fitted values, so a self-test on ",
            "this model measures parameter recovery, not recovery of the ",
            paste(hit, collapse = " or "), " process. Express it as a random ",
            "linkage (see vignette('environmental-linkages-and-priors')) to ",
            "have it drawn.", call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}


# A process the model gives no distribution to cannot be redrawn: M with
# M1_re = 0, or growth with no random linkage on a growth parameter. Nothing is
# drawn and nothing is wrong, but silence there is indistinguishable from a draw
# that worked -- which is how a switch that does nothing survives.
.sim_warn_process_absent <- function(obj, state) {
  procs <- c("recruitment", "M", "growth", "catchability", "selectivity")
  hit <- character(0)

  m1_re <- suppressWarnings(as.integer(obj$env$data$M1_re))
  if (state[2] == 1L && !any(m1_re > 0, na.rm = TRUE) &&
      !.sim_linkage_drawn(obj, c(0L, 1L, 0L, 0L, 0L))) {
    hit <- c(hit, "M (M1_re = 0 and no random linkage on M)")
  }
  if (state[3] == 1L && !.sim_linkage_drawn(obj, c(0L, 0L, 1L, 0L, 0L))) {
    hit <- c(hit, "growth (no random linkage on a growth parameter)")
  }

  if (length(hit)) {
    warning("Process error was requested for ", paste(hit, collapse = "; "),
            ", but this model gives ", if (length(hit) > 1) "those" else "that",
            " no distribution to draw from, so nothing was redrawn and the ",
            "fitted values stand. A self-test on this model measures parameter ",
            "recovery, not recovery of the process.", call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}


# Recruitment scored by two densities has no single distribution to draw from.
# Under srr_fun = 0 with srr_pred_fun > 0 -- the AMAK/Ianelli configuration --
# recruitment is R = exp(R0 + rec_dev) and the stock-recruit curve is fitted as a
# PENALTY on the same deviation, so rec_dev is scored by JNLL_REC_DEV and again
# by JNLL_SRR_PENALTY. Drawing at R_sd from the first alone ignores the second,
# and the second is not a density on rec_dev alone either -- it couples the
# deviations to the stock-recruit parameters, which are estimated too -- so there
# is no closed-form conditional to fall back on. The template draws nothing
# (ceattle.cpp section 5.13); this says so.
.sim_warn_rec_srr_penalty <- function(obj, state) {
  if (state[1] != 1L) return(invisible(FALSE))
  d <- obj$env$data
  if (!isTRUE(d$srr_fun == 0 && d$srr_pred_fun > 0)) return(invisible(FALSE))

  warning("Process error was requested for recruitment, but this model fits the ",
          "stock-recruit curve as a penalty on the recruitment deviations ",
          "(srr_fun = 0 with srr_pred_fun > 0, the AMAK/Ianelli configuration). ",
          "The deviations are then scored by two terms at once, which do not ",
          "compose into a single distribution to draw from, so the recruitment ",
          "deviations were not redrawn and their fitted values stand. Any random ",
          "linkage on a recruitment parameter is a separate latent and was drawn ",
          "normally. Set srr_fun to the stock-recruit function itself to have the ",
          "deviations simulated too.",
          call. = FALSE)
  invisible(TRUE)
}


# Truncation diagnostics for the natural-scale survey draws. An index cannot be
# negative and data_check() rejects one, so a draw from an untruncated normal has
# to be redrawn when it comes back non-positive; the template counts the attempts
# and rejections per row (see ceattle.cpp) and this reads them.
#
# Two families land here, for the same reason: `Normal` (drawn row by row) and
# `MVN`/`MVNORM` (the whole vector redrawn, since truncating each margin
# separately would break the correlation the likelihood exists to model). In both
# the redrawn data follow the normal truncated at zero while the likelihood
# scores the untruncated one, so a self-test built on them measures a different
# data-generating process. That gap is what these warnings are about.
#
# `TruncatedNormal` never enters here: it is fitted AND drawn as a normal
# left-truncated at zero, by inverse CDF, so its draw already follows its own
# likelihood. It is the fix for a `Normal` fleet that warns; there is no
# equivalent for MVN, which has no closed-form truncated sampler.
.sim_warn_index_truncated <- function(sim_rep, data_list) {
  tries   <- as.numeric(sim_rep$index_trunc_tries_sim)
  rejects <- as.numeric(sim_rep$index_trunc_rejects_sim)
  if (!length(tries) || !any(tries > 0)) return(invisible(FALSE))

  idx <- data_list$index_data
  obs <- idx$Observation
  drawn <- tries > 0
  # Per-ROW rate, worst row rather than the fleet mean: truncation bites one
  # marginal row at a time and an average over a wide fleet hides exactly the row
  # being resampled.
  frac <- ifelse(drawn, rejects / pmax(tries, 1), 0)

  fleet_of <- function(sel) unique(as.character(idx$Fleet_name[sel]))
  nonpos <- fleet_of(drawn & is.finite(obs) & obs <= 0)
  heavy  <- setdiff(fleet_of(drawn & frac > .SIM_INDEX_WARN_FRAC), nonpos)

  # The remedy differs by family, so name it: a `Normal` fleet has a truncated
  # counterpart to switch to, a covariance fleet does not.
  fam <- .index_family_codes(data_list$fleet_control$Index_distribution)
  names(fam) <- as.character(data_list$fleet_control$Fleet_name)
  remedy <- function(flts) {
    codes <- fam[flts]
    if (any(codes %in% 3L, na.rm = TRUE)) {
      paste0("An absolute sd that large relative to the index means the ",
             "untruncated normal is a poor observation model here. Use ",
             "Index_distribution = \"TruncatedNormal\" to fit and simulate the ",
             "same truncated density, or reduce the observation sd.")
    } else {
      paste0("Check the covariance against the scale of the index. MVN has no ",
             "closed-form truncated sampler, so there is no family to switch to.")
    }
  }

  if (length(nonpos)) {
    warning("Simulated a non-positive survey index for fleet(s) ",
            paste(nonpos, collapse = ", "), ", after up to ", .SIM_INDEX_MAX_TRIES,
            " redraws per row. data_check() requires Observation > 0, so ",
            "refitting this data set fails and self_test() counts the run as ",
            "not converged. ", remedy(nonpos), call. = FALSE)
  } else if (length(heavy)) {
    warning("Simulating the survey index for fleet(s) ",
            paste(heavy, collapse = ", "), " needed a non-positive draw ",
            "redrawn more than ", round(100 * .SIM_INDEX_WARN_FRAC),
            "% of the time. The draws are positive, but they come from the ",
            "normal truncated at zero rather than the untruncated normal the ",
            "likelihood scores, so a self_test() built on them tests a ",
            "different data-generating process. ", remedy(heavy), call. = FALSE)
  }
  invisible(length(nonpos) > 0 || length(heavy) > 0)
}


#' Warn where a covariance survey fleet carries rows it cannot simulate
#'
#' The `MVN`/`MVNORM` draw is a correlated block, and the covariance is supplied
#' for the fleet's FITTED observations only -- there is no covariance entry for a
#' year the model does not fit, so those rows keep the values they came in with.
#' The independent families have no such limit and redraw every row.
#'
#' That matters because `run_mse()` reveals the operating model's negative-`Year`
#' rows to the estimation model as each assessment's new data. For a covariance
#' fleet those rows were never drawn, so the estimation model would be handed
#' survey observations carried over unchanged, and a management-strategy
#' evaluation would report performance against a survey that never varied.
#'
#' @param data_list Data list being simulated into.
#' @noRd
.sim_warn_index_unsimulated <- function(data_list) {
  fc  <- data_list$fleet_control
  idx <- data_list$index_data
  if (is.null(fc) || is.null(idx) || !nrow(idx)) return(invisible(FALSE))
  fam <- .index_family_codes(fc$Index_distribution)
  if (!length(fam)) return(invisible(FALSE))

  hit <- character(0)
  for (i in seq_len(nrow(fc))) {
    if (!(fam[i] %in% c(1L, 2L))) next
    rows <- which(idx$Fleet_code == fc$Fleet_code[i])
    if (!length(rows)) next
    outside <- idx$Year[rows] <= 0 | idx$Year[rows] > data_list$endyr
    if (any(outside)) hit <- c(hit, as.character(fc$Fleet_name[i]))
  }
  if (length(hit)) {
    warning("Fleet(s) ", paste(unique(hit), collapse = ", "),
            " use a covariance (MVN/MVNORM) survey likelihood, whose covariance ",
            "covers only the years the model fits, so their observations ",
            "outside that window were NOT simulated and carry their original ",
            "values. run_mse() reveals those rows to the estimation model as new ",
            "survey data, so an MSE on this model evaluates against a survey ",
            "that does not vary.", call. = FALSE)
  }
  invisible(length(hit) > 0)
}


#' Check a simulated matrix has the columns the data frame expects
#'
#' The write-back takes the first `n` columns positionally. `comp_obs` is built
#' with `dplyr::select(contains("Comp_"))` while [.composition_cols()] matches
#' `^Comp_`, so a column named like `Total_Comp_n` would put the two out of step
#' and every bin would land one place over -- silently, since the copy would
#' still be the right shape. The per-row assignment this replaced would have
#' errored on the length mismatch; assert it instead.
#'
#' @param n_sim Columns the template returned.
#' @param n_dat Bin columns in the data frame.
#' @param what Data type name.
#' @noRd
.sim_check_cols <- function(n_sim, n_dat, what) {
  if (n_sim < n_dat) {
    stop("sim_mod(): the fitted model simulated ", n_sim, " ", what,
         " bin columns but its data_list holds ", n_dat,
         ". The two must describe the same bins.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Warn where a composition row came back with nothing in it
#'
#' A row sums to zero when the model predicts no composition for it -- a fleet
#' with no catch that year, or one switched off -- and also when the effective
#' sample size rounds away, which `Sample_size * weight` can do for a heavily
#' down-weighted fleet. The first is meaningful and expected; the second is a
#' fleet silently dropping out of the refit, so name both rather than let a
#' `self_test()` score a replicate whose data quietly went missing.
#'
#' @param x Simulated bin matrix.
#' @param fleet Fleet code per row.
#' @param what Data type name.
#' @noRd
.sim_warn_empty_comp <- function(x, fleet, what) {
  empty <- rowSums(x, na.rm = TRUE) <= 0
  if (any(empty)) {
    warning(sum(empty), " ", what, " row(s) came back empty for fleet(s) ",
            paste(unique(fleet[empty]), collapse = ", "),
            ". The model predicts no composition for them, or their sample size ",
            "times the composition weight rounded to zero. Those rows carry no ",
            "information into a refit.", call. = FALSE)
  }
  invisible(any(empty))
}


#' Resolve `process` to the template's `simulate_state` vector
#'
#' Process error is off unless asked for, because redrawing a process changes
#' what a self-test measures -- from "can the estimator recover these
#' parameters" to "can it recover this process". The slots follow
#' `ceattle.cpp`: recruitment, M and growth are the population dynamics;
#' selectivity and catchability are the observation process.
#'
#' `"recruitment"` covers the initial age structure as well as the annual
#' deviations. The initial numbers-at-age are recruitment from the years before
#' `styr` -- `init_dev` and `rec_dev` share `R_sd` and the same bias correction,
#' and `rec_dev` at year 0 feeds the initial scalar -- so they are one process.
#'
#' @param process `FALSE` / `"none"` for none, `TRUE` / `"all"` for every
#'   process, `"dynamics"` for the population dynamics, `"observation"` for the
#'   observation process, or any subset of `"recruitment"`, `"M"`, `"growth"`,
#'   `"selectivity"`, `"catchability"`.
#' @return Integer vector of length 5.
#' @noRd
.sim_state_codes <- function(process) {
  # Slot order is the linkage process code (see .RCE_LINKAGE_PROCESS in
  # 0-linkage_encode.R), so the template can index simulate_state by a
  # linkage row's process directly. Keep the two in step.
  slots <- c("recruitment", "M", "growth", "catchability", "selectivity")
  out <- rep(0L, length(slots))
  if (is.null(process) || isFALSE(process)) return(out)
  if (isTRUE(process)) return(rep(1L, length(slots)))

  chr <- as.character(process)
  if (length(chr) == 1L && chr %in% c("none", "all", "dynamics", "observation")) {
    if (chr == "none") return(out)
    if (chr == "all")  return(rep(1L, length(slots)))
    chr <- if (chr == "dynamics") slots[1:3] else slots[4:5]
  }
  bad <- setdiff(chr, slots)
  if (length(bad)) {
    stop("sim_mod(): unknown process(es) ", paste(bad, collapse = ", "),
         ". Use FALSE, TRUE, \"dynamics\", \"observation\", or any of ",
         paste(slots, collapse = ", "), ".", call. = FALSE)
  }
  out[match(chr, slots)] <- 1L
  out
}


#' Get a TMB object that can run the template's `SIMULATE` blocks
#'
#' The draws live beside their densities in `ceattle.cpp`, so simulating means
#' evaluating the compiled model. Usually the fit's own `TMB::MakeADFun` object
#' does it, including one saved and reloaded (TMB re-tapes from the stored
#' `obj$env$data`).
#'
#' A fit whose `$obj` was dropped to save space still carries its `data_list` and
#' its estimates, which is all `.refit_like()` needs to rebuild the model without
#' re-estimating it. The rebuild is the same model rather than an approximation,
#' but that is checked rather than assumed -- see [.sim_check_rebuild()].
#'
#' Rebuild at the mode the source was fitted under. `estimateMode = 3` builds the
#' objective and stops, which leaves `catch_hat` at 0 on every projection row, so
#' a fit that ran a projection needs one too -- mode 2 runs it from the supplied
#' estimates without re-optimizing the hindcast.
#'
#' `model_average()` output cannot be rebuilt and is not meant to be. It drops
#' the estimates along with `$obj` and replaces `quantities` with an average over
#' models, so no parameter vector produced those expected values and there is
#' nothing for a likelihood to draw around.
#'
#' @param Rceattle A fitted `Rceattle` model.
#' @return A `TMB::MakeADFun` object at the model's estimates.
#' @noRd
.sim_obj <- function(Rceattle) {
  if (!is.null(Rceattle$obj)) return(Rceattle$obj)

  if (is.null(Rceattle$estimated_params) || is.null(Rceattle$data_list)) {
    stop("sim_mod(simulate = TRUE) needs a fitted model, and this one carries ",
         "neither a TMB object nor the estimates to rebuild one. Averaged ",
         "models (model_average()) are the usual case: their quantities are an ",
         "average over models rather than any one model's fit, so there are no ",
         "parameters to simulate around. Simulate from one of the fits instead.",
         call. = FALSE)
  }

  src_mode <- Rceattle$data_list$estimateMode
  mode <- if (!is.null(src_mode) && src_mode %in% c(0, 2)) 2 else 3

  rebuilt <- suppressMessages(.refit_like(
    data_list    = Rceattle$data_list,
    inits        = Rceattle$estimated_params,
    estimateMode = mode))

  .sim_check_rebuild(Rceattle, rebuilt)
  rebuilt$obj
}


#' Confirm a rebuilt model is the one whose estimates it was given
#'
#' `.refit_like()` reconstructs the model specification from the `data_list`, so
#' a `data_list` edited since the fit -- or one that no longer records everything
#' the fit used -- rebuilds something else, which would then simulate quite
#' happily around the wrong expected values.
#'
#' Check what the draws are made of: the predicted catch, the observation
#' standard deviation and the bias-adjustment multiplier (`ceattle.cpp`, catch
#' slot), plus the predicted diet composition. Extend the list as later stages
#' add their draws.
#'
#' `catch_hat` and `log_catch_sd` are checked against the fit's own stored
#' values, so an edited `data_list` shows up as a mismatch. `bias_adjust_obs`
#' cannot be checked that way: it enters neither of them and nothing outside the
#' `data_list` records it, so the rebuild simply honours whatever the
#' `data_list` now says. What is worth catching is its ABSENCE --
#' `.refit_bias_adjust()` falls back to the `fit_control()` default when the
#' field is missing, so a fit made with the correction off would rebuild with it
#' on and every simulated catch would be mis-centred by `exp(sigma^2 / 2)`, with
#' nothing to reveal it.
#'
#' A quantity that is absent is not a pass: a rebuild that cannot be checked
#' cannot be trusted.
#'
#' @param Rceattle The fit being simulated from.
#' @param rebuilt The model returned by `.refit_like()`.
#' @noRd
.sim_check_rebuild <- function(Rceattle, rebuilt) {
  for (q in c("catch_hat", "log_catch_sd", "diet_hat")) {
    want <- Rceattle$quantities[[q]]
    got  <- rebuilt$quantities[[q]]
    if (is.null(want) || is.null(got)) {
      stop("sim_mod(): this model has no TMB object, and the rebuild cannot be ",
           "checked because '", q, "' is missing, so there is no way to tell ",
           "whether it is the same model. Refit with fit_mod() and simulate ",
           "from that fit.", call. = FALSE)
    }
    if (length(want) != length(got) ||
        !isTRUE(all.equal(as.numeric(want), as.numeric(got), tolerance = 1e-8))) {
      stop("sim_mod(): this model has no TMB object, and rebuilding one from ",
           "its data_list gives a different '", q, "', so the two no longer ",
           "describe the same model. This happens when the data_list has been ",
           "edited since the fit, or when the fit used a setting the data_list ",
           "does not record. Refit with fit_mod() and simulate from that fit.",
           call. = FALSE)
    }
  }

  ba <- Rceattle$data_list$bias_adjust_obs
  if (is.null(ba) || !length(ba) || anyNA(ba)) {
    stop("sim_mod(): this model has no TMB object, and its data_list does not ",
         "record 'bias_adjust_obs', so a rebuild cannot reproduce how the ",
         "original draw was centred. Refit with fit_mod() and simulate from ",
         "that fit.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Run the template's `SIMULATE` blocks once
#'
#' Called once per `sim_mod()` call. Each `obj$simulate()` re-runs the whole
#' model and draws every simulated quantity, so calling it per data type would
#' give each type a draw from a different replicate.
#'
#' `par` is passed explicitly because TMB defaults to `obj$env$last.par`, the
#' last parameter vector evaluated, which after an optimization can be a step the
#' optimizer rejected. `last.par.best` is the MLE, and is also where
#' `fit$quantities` was reported from, so the draws and the expected values agree.
#'
#' The draw lasts only for this call: TMB re-reads `DATA_` objects each
#' evaluation, so the object's data and objective are unchanged, and the `jnll`
#' returned is still the original data's. Fit simulated data by building a new
#' model on them, as `self_test()` and `run_mse()` do.
#'
#' @param obj A `TMB::MakeADFun` object from [.sim_obj()].
#' @return The report list, including the simulated `*_obs_sim` matrices.
#' @noRd
.sim_draw <- function(obj, state = NULL, period = NULL) {
  if (!is.null(state) || !is.null(period)) {
    # obj$env is by reference, so an override has to be undone or it would
    # follow the caller's fitted object around for the rest of the session.
    #
    # Written as DOUBLE, not integer. fit_mod() sanitizes every DATA_ element to
    # double before MakeADFun, and TMB re-reads the stored list on each
    # evaluation, so handing back an integer where it stored a double fails with
    # "Error when reading the variable" -- a DATA_IVECTOR reads a double vector
    # perfectly well, it is the change of storage type that breaks it.
    old <- obj$env$data[c("simulate_state", "simulate_period")]
    on.exit(obj$env$data[names(old)] <- old, add = TRUE)
    if (!is.null(state))  obj$env$data$simulate_state  <- as.double(state)
    if (!is.null(period)) obj$env$data$simulate_period <- as.double(period)
  }
  tryCatch(
    obj$simulate(par = obj$env$last.par.best),
    error = function(e) {
      # TMB re-tapes a stored model against the template loaded now. If that
      # template wants an input the stored data list lacks, the failure is a bare
      # "Error when reading the variable" with nothing to act on. Only name that
      # cause when it is the one that happened -- any other error is passed
      # through as itself rather than given a diagnosis it may not have.
      msg <- conditionMessage(e)
      if (grepl("reading the variable", msg, fixed = TRUE)) {
        stop("sim_mod(): the model's stored data does not match the TMB ",
             "template now loaded (", msg, "). This happens when the model was ",
             "fitted with a different version of Rceattle. Refit it with this ",
             "version and simulate from that fit.", call. = FALSE)
      }
      stop(e)
    })
}


#' Check a simulated matrix lines up with the data frame it is written into
#'
#' The write-back copies by row position: `rearrange_data()` builds each `*_obs`
#' matrix straight from its data frame, so row `i` of one is row `i` of the
#' other. Editing `data_list` after the fit breaks that. The old R draw handed
#' the mismatched vectors to `rnorm()`, which recycles the shorter one silently
#' and returns a full-length wrong answer; fail instead.
#'
#' @param n_sim Rows the template returned.
#' @param n_dat Rows of the data frame being written.
#' @param what Name of the data type, for the message.
#' @noRd
.sim_check_rows <- function(n_sim, n_dat, what) {
  if (n_sim != n_dat) {
    stop("sim_mod(): the fitted model simulated ", n_sim, " ", what,
         " rows but its data_list holds ", n_dat,
         ". The two must line up row for row. This happens when data_list is ",
         "edited after the fit -- simulate from the model as it was fitted.",
         call. = FALSE)
  }
  invisible(TRUE)
}


#' Pull a simulated observation matrix out of the report list
#'
#' A template with no `SIMULATE` block for this data type returns a report
#' lacking the entry, and the write-back would put `NULL` into the data frame.
#' Since TMB always uses the template currently loaded, this really only fires on
#' a stale shared object in a development session.
#'
#' @param rep Report list from [.sim_draw()].
#' @param name Name of the simulated object, e.g. `"catch_obs_sim"`.
#' @return The matrix.
#' @noRd
.sim_report_obs <- function(rep, name) {
  out <- rep[[name]]
  if (is.null(out)) {
    stop("sim_mod(): the loaded TMB template returned no simulated '", name,
         "', so it has no simulator for this data type. Recompile the package ",
         "(pkgload::load_all() or R CMD INSTALL) so the loaded model matches ",
         "the installed source.", call. = FALSE)
  }
  out
}


#' Warn when a simulated observation is one the model cannot be refit on
#'
#' A non-finite or negative draw is not quietly dropped: `data_check()` rejects
#' it, the refit errors, and `self_test()` counts the replicate as not converged
#' -- which reads as a convergence problem rather than a data one. The usual
#' cause is an observation standard deviation that never got a value
#' (`Estimate_catch_sd = 1` with `Catch_sd` blank gives `exp(log(NA))`), or, for a
#' natural-scale survey draw, observation error large relative to the index.
#'
#' @param x Simulated observations.
#' @param fleet Fleet code per element, for the message.
#' @param what Data type name.
#' @param strictly_positive `TRUE` where zero is invalid too. `data_check()`
#'   requires a survey index above zero, but accepts a catch of exactly zero --
#'   and a fishery closed in a projection year legitimately draws one.
#' @noRd
.sim_warn_unusable <- function(x, fleet, what, strictly_positive = FALSE) {
  bad <- if (strictly_positive) !is.finite(x) | x <= 0 else !is.finite(x) | x < 0
  if (any(bad)) {
    warning("Simulated an unusable ", what, " for fleet(s) ",
            paste(unique(fleet[bad]), collapse = ", "),
            ". data_check() rejects those observations, so refitting this data ",
            "set fails and self_test() counts the run as not converged. Check ",
            "the fleet's observation standard deviation against the scale of ",
            "the data.", call. = FALSE)
  }
  invisible(x)
}


#' Simulate Rceattle data
#'
#' @description Simulates the data an Rceattle model would have produced, either
#' as expected values or as a random draw. Every observation type is covered:
#' survey biomass (under the fleet's own \code{Index_distribution} -- lognormal,
#' natural-scale normal, or the correlated MVN/MVNORM draw from its covariance),
#' total catch (lognormal), age/length composition and conditional
#' age-at-length (multinomial or Dirichlet-multinomial), and stomach contents
#' (multinomial or Dirichlet-multinomial).
#'
#' @details
#' Every draw is taken by the TMB model itself, in a \code{SIMULATE} block beside
#' the likelihood that defines it, so the simulated data and the density that
#' will be fitted to them cannot drift apart. Two copies of one observation model
#' drift silently, and a \code{\link{self_test}} then reports recovery against a
#' process the likelihood never assumed.
#'
#' Consequently \code{simulate = TRUE} needs a model to evaluate. A model loaded
#' from an \code{.Rdata}/\code{.rds} file has one, and a fit whose \code{$obj}
#' was dropped to save space is rebuilt from its \code{data_list} and estimates,
#' provided the rebuild reproduces the fit's own expected values.
#' \code{\link{model_average}} output cannot be simulated from at all: its
#' quantities are an average over models rather than any one model's fit, so no
#' parameters produced them. \code{simulate = FALSE} reads only
#' \code{$quantities} and draws no random numbers, so it works on any model.
#'
#' Rows the model predicts nothing for are not invented. A composition for a
#' fleet with no catch that year comes back empty, a stomach whose predator has
#' an empirical suitability keeps its observed proportions, and a covariance
#' survey fleet's observations outside its fitted window keep theirs; each of
#' those is reported by a warning rather than passed over in silence.
#'
#' Simulating leaves the drawn values in the object's report environment, under
#' names ending \code{_sim}. The estimates, the data and the objective function
#' are untouched, so a later \code{osa_residuals()} or \code{vcov()} on the same
#' model is unaffected.
#'
#' @param Rceattle A CEATTLE model object exported from \code{Rceattle}.
#' @param simulate Logical. If \code{TRUE}, simulates data from distributions.
#'   If \code{FALSE}, returns the expected values (hats).
#' @param process Which process error to redraw alongside the observations.
#'   \code{FALSE} (default) or \code{"none"} keeps the fitted deviations;
#'   \code{TRUE} or \code{"all"} redraws every process; \code{"dynamics"} covers
#'   recruitment, natural mortality and growth, \code{"observation"} covers
#'   catchability and selectivity. Any subset of \code{"recruitment"},
#'   \code{"M"}, \code{"growth"}, \code{"catchability"} and
#'   \code{"selectivity"} may be given instead. Ignored when
#'   \code{simulate = FALSE}.
#'
#' @return A \code{data_list} object containing the simulated or expected data
#'   values, formatted for use in \code{Rceattle}. When \code{process} redrew
#'   something, the deviations that generated the data are attached as
#'   \code{attr(x, "process_sim")} -- a named list holding whichever of
#'   \code{rec_dev}, \code{init_dev}, \code{log_M1_dev} and
#'   \code{beta_linkage_re} were drawn. Those are the truth a refit has to
#'   recover; without them the only comparison available is against the original
#'   fitted deviations, which are no longer the values that generated the data.
#' @export
#'
sim_mod <- function(Rceattle, simulate = FALSE, process = FALSE) {
  dat_sim <- Rceattle$data_list
  quantities <- Rceattle$quantities

  # Indices of abundance/biomass ----
  index_hat <- quantities$index_hat

  if (simulate) {
    # Every simulated observation in one call. obj$simulate() re-runs the whole
    # model and draws every type the template covers, so this is the only place
    # it is called; the catch and diet blocks below read the same report.
    #
    # The survey draw follows each fleet's own Index_distribution (ceattle.cpp,
    # slot 0): lognormal, natural-scale normal, or a correlated draw from the
    # fleet's covariance. Drawing every fleet as lognormal, as this once did,
    # simulates from an observation model the likelihood does not assume -- and a
    # self-test then reports recovery as though it had.
    sim_obj <- .sim_obj(Rceattle)
    sim_state <- .sim_state_codes(process)
    sim_rep <- .sim_draw(sim_obj, state = sim_state)
    .sim_warn_linkage_qar1(sim_obj, sim_state)
    .sim_warn_process_unsimulated(dat_sim, sim_state)
    .sim_warn_process_absent(sim_obj, sim_state)
    .sim_warn_rec_srr_penalty(sim_obj, sim_state)
    index_sim <- .sim_report_obs(sim_rep, "index_obs_sim")
    .sim_check_rows(nrow(index_sim), nrow(dat_sim$index_data), "index")
    dat_sim$index_data$Observation <-
      .sim_warn_unusable(as.numeric(index_sim[, 1]),
                         dat_sim$index_data$Fleet_code, "survey index",
                         strictly_positive = TRUE)
    .sim_warn_index_unsimulated(dat_sim)
    .sim_warn_index_truncated(sim_rep, dat_sim)
  } else {
    # Expected value
    dat_sim$index_data$Observation <- index_hat
  }


  # Age/Length composition ----
  # Drawn by the template, in RAW bin space (ceattle.cpp, slot 2). Tail
  # accumulation folds bins before the density and the fold has no inverse, so a
  # draw taken there could not be written back; drawing raw and letting the refit
  # fold again is exact, because both families are closed under merging
  # categories.
  #
  # Counts are stored, as the R draw stored them, and rearrange_data() normalizes
  # each row on the next fit. A row the model predicts nothing for -- a fleet
  # with no catch that year, or one switched off -- comes back as the prediction,
  # which is zero, NOT as the values it went in with. run_mse() relies on that:
  # it reads a zero row as "not sampled" and drops the sample size with it.
  if (nrow(dat_sim$comp_data) > 0) {
    comp_cols <- .composition_cols(dat_sim$comp_data, "Comp_")
    if (simulate) {
      comp_sim <- .sim_report_obs(sim_rep, "comp_obs_sim")
      .sim_check_rows(nrow(comp_sim), nrow(dat_sim$comp_data), "composition")
      .sim_check_cols(ncol(comp_sim), length(comp_cols), "composition")
      dat_sim$comp_data[, comp_cols] <- comp_sim[, seq_along(comp_cols), drop = FALSE]
      .sim_warn_empty_comp(comp_sim, dat_sim$comp_data$Fleet_code, "composition")
    } else {
      dat_sim$comp_data[, comp_cols] <-
        quantities$comp_hat[, seq_along(comp_cols), drop = FALSE]
    }
  }


  # CAAL ----
  # As for the marginal composition above; CAAL has no tail accumulation.
  if (nrow(dat_sim$caal_data) > 0) {
    caal_cols <- .composition_cols(dat_sim$caal_data, "CAAL_")
    if (simulate) {
      caal_sim <- .sim_report_obs(sim_rep, "caal_obs_sim")
      .sim_check_rows(nrow(caal_sim), nrow(dat_sim$caal_data), "CAAL")
      .sim_check_cols(ncol(caal_sim), length(caal_cols), "CAAL")
      dat_sim$caal_data[, caal_cols] <- caal_sim[, seq_along(caal_cols), drop = FALSE]
      .sim_warn_empty_comp(caal_sim, dat_sim$caal_data$Fleet_code, "CAAL")
    } else {
      dat_sim$caal_data[, caal_cols] <-
        quantities$caal_hat[, seq_along(caal_cols), drop = FALSE]
    }
  }


  # Catch ----
  catch_hat <- quantities$catch_hat

  if (simulate) {
    # Drawn by the template's SIMULATE block, beside the catch density that
    # defines it (ceattle.cpp, slot 1), rather than re-derived here. Two copies
    # of one observation model drift apart silently, and a self-test then reports
    # recovery against a process the likelihood never assumed.
    #
    # Read from the draw taken in the index block above -- one obj$simulate() per
    # sim_mod() call. Calling it again here would give catch a draw from a
    # different replicate than the index, and consume twice the random numbers.
    catch_sim <- .sim_report_obs(sim_rep, "catch_obs_sim")
    .sim_check_rows(nrow(catch_sim), nrow(dat_sim$catch_data), "catch")
    # Column 1 is the observation; the template writes the natural scale there
    # (obsvec holds its log). Column 2 is the supplied sd, untouched.
    dat_sim$catch_data$Catch <-
      .sim_warn_unusable(as.numeric(catch_sim[, 1]),
                         dat_sim$catch_data$Fleet_code, "catch")
  } else {
    # Expected values
    dat_sim$catch_data$Catch <- catch_hat
  }


  # Diet (stomach content) ----
  # Drawn by the template alongside the catch, under each predator's own
  # Diet_distribution (ceattle.cpp, section 13.2). Until this was added, a
  # multispecies self_test() resampled every other data type and refit against
  # the same stomachs every replicate, so suitability was recovered from data
  # that never varied and the test read better than it was.
  #
  # Only rows belonging to a stomach the model actually fits are redrawn: the
  # template skips a predator whose suitability is not estimated, and those rows
  # come back carrying the values they went in with. The "other prey" balance is
  # not stored -- it is recomputed from the prey proportions on the next fit.
  if (simulate && !is.null(dat_sim$diet_data) && nrow(dat_sim$diet_data) > 0) {
    diet_sim <- .sim_report_obs(sim_rep, "diet_obs_sim")
    .sim_check_rows(nrow(diet_sim), nrow(dat_sim$diet_data), "diet")
    # Column 2 is Stomach_proportion_by_weight; column 1 is Sample_size.
    kept <- as.numeric(diet_sim[, 2])
    was  <- dat_sim$diet_data$Stomach_proportion_by_weight
    dat_sim$diet_data$Stomach_proportion_by_weight <- kept

    # Warn per predator, not once for the whole table. A model that estimates
    # suitability for some predators and not others redraws only the former, and
    # an aggregate test cannot see that: on BS2017MS with suitMode = c(4, 0, 4)
    # the middle predator's stomachs are frozen while the table as a whole
    # changes. Rows can also come back untouched because the stomach's sample
    # size rounds to zero.
    #
    # This only matters where predation is modelled. Under empirical suitability
    # the stomach proportions set suitability directly (predation.hpp,
    # calculate_msvpa_suitability) and hence predation mortality, so a
    # self_test() that holds them fixed makes recovery of predation look better
    # than it is. In a single-species model the diet rows are inert and there is
    # nothing to say.
    if (!is.null(dat_sim$msmMode) && any(dat_sim$msmMode > 0)) {
      pred_sp <- dat_sim$diet_data$Pred
      frozen <- vapply(split(seq_along(kept), pred_sp),
                       function(ix) all(abs(kept[ix] - was[ix]) < 1e-14),
                       logical(1))
      if (any(frozen)) {
        sp <- names(frozen)[frozen]
        nm <- if (!is.null(dat_sim$spnames)) {
          paste(dat_sim$spnames[as.integer(sp)], collapse = ", ")
        } else {
          paste(sp, collapse = ", ")
        }
        warning("sim_mod() did not simulate the diet data for predator(s) ", nm,
                ": the model predicts no stomach composition for them, usually ",
                "because their suitability is empirical (suitMode = 0) rather ",
                "than estimated. Those stomach proportions still set suitability ",
                "and predation mortality, so a self_test() holds them fixed and ",
                "recovery of predation is optimistic.", call. = FALSE)
      }
    }
  }

  # When process error was redrawn, the deviations that generated these data are
  # the truth a self-test has to recover. They are attached rather than added as
  # a list element so the return value is still a plain data_list -- every
  # existing caller keeps working, and fit_mod() ignores the attribute.
  if (simulate && any(sim_state == 1L)) {
    attr(dat_sim, "process_sim") <- .sim_process_truth(sim_rep, sim_state, sim_obj)
  }
  return(dat_sim)
}


#' The process deviations a simulated data set was generated from
#'
#' Only the processes that were actually redrawn are returned: handing back a
#' deviation that was NOT drawn would look like a truth to recover, when it is
#' just the fitted value the refit is starting from.
#'
#' A process is only listed if the model gives it a distribution to draw from.
#' Asking for M on a model with `M1_re = 0` draws nothing, so returning
#' `log_M1_dev` there would hand back the fitted values dressed as a truth, and a
#' self-test would report perfect recovery of a process it never simulated.
#'
#' @param sim_rep The report from one `obj$simulate()` call.
#' @param state Integer switch vector from `.sim_state_codes()`.
#' @param obj The simulated object, for the model's own switches.
#' @return Named list of the drawn deviations, or NULL if none were drawn.
#' @noRd
.sim_process_truth <- function(sim_rep, state, obj) {
  out <- list()
  # rec_srr_single_density is the template's own gate on the recruitment draw
  # (ceattle.cpp section 5.13): FALSE under the AMAK/Ianelli configuration, where
  # the stock-recruit penalty scores rec_dev a second time and there is no single
  # distribution to draw from. Nothing was drawn there, so nothing is returned.
  if (state[1] == 1L && isTRUE(as.logical(sim_rep$rec_srr_single_density))) {
    out$rec_dev  <- sim_rep$rec_dev_sim
    out$init_dev <- sim_rep$init_dev_sim
  }
  m1_re <- suppressWarnings(as.integer(obj$env$data$M1_re))
  if (state[2] == 1L && any(m1_re > 0, na.rm = TRUE)) {
    out$log_M1_dev <- sim_rep$log_M1_dev_sim
  }
  # One vector covering every random linkage, in the registry's slot order; the
  # linkage table says which process and parameter each slot belongs to. Only
  # attached when a group belonging to a requested process was actually drawn --
  # otherwise `!is.null(attr(x, "process_sim"))` would report process error on a
  # model that has none.
  re <- sim_rep$beta_linkage_re_sim
  if (length(re) && .sim_linkage_drawn(obj, state)) out$beta_linkage_re <- re
  if (!length(out)) NULL else out
}


#' Does this model hold a random linkage on a process the caller asked for?
#'
#' Mirrors the gate in ceattle.cpp section 5.12b, including its skip of observed
#' AR1 (QAR1) groups, so R and the template agree on what was drawn.
#'
#' @param obj The simulated object.
#' @param state Integer switch vector from `.sim_state_codes()`.
#' @noRd
.sim_linkage_drawn <- function(obj, state) {
  d <- obj$env$data
  idx <- d$linkage_re_index
  if (is.null(idx) || !length(idx)) return(FALSE)
  for (i in seq_along(idx)) {
    slot <- idx[i]
    if (is.na(slot) || slot < 0) next
    p <- d$linkage_process[i]
    if (is.na(p) || p < 0 || p >= length(state) || state[p + 1L] != 1L) next
    grp <- d$linkage_re_sigma[slot + 1L]
    if (!is.null(d$linkage_re_obs) && d$linkage_re_obs[grp + 1L] >= 0) next
    return(TRUE)
  }
  FALSE
}

#' Sample historical recruitment deviations into the projection
#'
#' @param Rceattle CEATTLE model object exported from \code{Rceattle}
#' @param sample_rec Include resampled recruitment deviations from the hindcast in the OM projection. Resampled deviations are used rather than drawing from N(0, sigmaR) because the initial deviations bias R0 low. If FALSE, uses the mean recruitment deviation.
#' @param update_model Update model dynamics. Default = TRUE
#' @param rec_trend Linear increase or decrease in mean recruitment from \code{endyr} to \code{projyr}. This is the terminal multiplier \code{mean rec * (1 + (rec_trend/projection years) * 1:projection years)}. Can be of length 1 or of length nspp. If length 1, all species get the same trend.
#'
#' @returns Rceattle model
#' @export
#'
sample_rec <- function(Rceattle, sample_rec = TRUE, update_model = TRUE, rec_trend = 0){

  # Years for simulations
  hind_yrs <- (Rceattle$data_list$styr) : Rceattle$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  proj_yrs <- (Rceattle$data_list$endyr + 1) : Rceattle$data_list$projyr
  proj_nyrs <- length(proj_yrs)

  # - Adjust rec trend
  if(length(rec_trend)==1){
    rec_trend = rep(rec_trend, Rceattle$data_list$nspp)
  }

  # Replace future rec devs ----
  #FIXME - update non-sample rec for stock recruit relationship
  for(sp in 1:Rceattle$data_list$nspp){

    # -- where SR curve is estimated directly
    if(Rceattle$data_list$srr_fun == Rceattle$data_list$srr_pred_fun){
      if(sample_rec){ # Sample devs from hindcast
        rec_dev <- sample(x = Rceattle$estimated_params$rec_dev[sp, 1:hind_nyrs], size = proj_nyrs, replace = TRUE) + log((1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
      } else{ # Set to mean rec otherwise
        rec_dev <- log(mean(Rceattle$quantities$R[sp,1:hind_nyrs]) * (1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs))  - log(Rceattle$quantities$R0[sp]) # - Scale mean rec for rec trend
      }
    }

    # -- OMs where SR curve is estimated as penalty (sensu Ianelli)
    if(Rceattle$data_list$srr_fun != Rceattle$data_list$srr_pred_fun){
      if(sample_rec){ # Sample devs from hindcast
        rec_dev <- sample(x = (log(Rceattle$quantities$R) - log(Rceattle$quantities$R_hat))[sp, 1:hind_nyrs],
                          size = proj_nyrs, replace = TRUE) + log((1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
      } else{ # Set to mean rec otherwise
        # `log(R) - log(R_hat)` is already a log-scale deviation centred near
        # zero, so its mean is routinely negative and the outer log() this
        # used to take returned NaN. Take the mean deviation directly and add
        # the log trend, mirroring the sampling branch above.
        rec_dev <- mean((log(Rceattle$quantities$R) - log(Rceattle$quantities$R_hat))[sp, 1:hind_nyrs]) +
          log((1+(rec_trend[sp]/proj_nyrs) * 1:proj_nyrs)) # - Scale mean rec for rec trend
      }
    }

    # - Update OM with devs
    Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1] <- replace(
      Rceattle$estimated_params$rec_dev[sp,proj_yrs - Rceattle$data_list$styr + 1],
      values =  rec_dev)
  }

  if(update_model){
    # * Update fit ----
    estMode <- Rceattle$data_list$estimateMode
    Rceattle <-
      suppressWarnings(
        suppressMessages(
          # Rebuild so the resampled projection recruitment propagates into the
          # reported dynamics. .refit_like() carries the source model's whole
          # specification across the refit -- HCR, stock-recruit, M, growth, and
          # the catchability / selectivity / composition linkages.
          .refit_like(
            data_list    = Rceattle$data_list,
            inits        = Rceattle$estimated_params,
            estimateMode = 3,
            getsd        = TRUE)
        )
      )
    Rceattle$data_list$estimateMode <- estMode
  }

  return(Rceattle)
}

#' Evaluate simulation performance
#'
#' @description Function to evaluate the simulation performance with regard to bias using the median relative error (MRE) and precision using the coefficient of variation.
#'
#' @param operating_mod CEATTLE model object exported from \code{Rceattle} to be used as the operating model
#' @param simulation_mods List of CEATTLE model objects exported from \code{Rceattle} fit to simulated data
#' @param object character string specifying which part of the model to compare (default = "quantities")
#' @return A data frame summarizing simulation performance metrics
#' @note Compares against \code{operating_mod} only, so it is valid when the
#'   replicates redrew the observations alone. If they were produced with
#'   \code{sim_mod(process = )} / \code{self_test(process = )}, the operating
#'   model's deviations are no longer what generated the data and the bias this
#'   reports is an artefact. Compare against
#'   \code{attr(sims, "process_sim")} in that case.
#' @export
compare_sim <- function(operating_mod, simulation_mods, object = "quantities") {
  # TODO update

  # Get differences
  sim_mre <- list()
  sim_mse <- list()
  sim_mean <- list()
  sim_median <- list()
  sim_sd <- list()
  sim_cv <- list()
  sim_params <- list()

  for (j in 1:length(names(operating_mod[[object]]))) {
    param <- names(operating_mod[[object]])[j]

    sim_mre[[param]] <- list()
    sim_mse[[param]] <- list()
    sim_mean[[param]] <- list()
    sim_sd[[param]] <- list()
    sim_cv[[param]] <- list()
    sim_params[[param]] <- list()

    om_params <- operating_mod[[object]][[param]]

    for (i in 1:length(simulation_mods)) {

      sm_params <- simulation_mods[[i]][[object]][[param]]

      sim_params[[param]][[i]] <- sm_params
      sim_mre[[param]][[i]] <- (sm_params - om_params)/om_params
      sim_mse[[param]][[i]] <- (sm_params - om_params)^2
    }

    param_dim <- length(dim(om_params))

    # If 1 value
    if (param_dim == 0) {
      sim_mean[[param]] <- mean(unlist(sim_params[[param]]))
      sim_sd[[param]] <- sd(unlist(sim_params[[param]]))
      sim_cv[[param]] <- sim_sd[[param]]/sim_mean[[param]]

      sim_mre[[param]] <- median(unlist(sim_mre[[param]]))
      sim_mse[[param]] <- mean(unlist(sim_mse[[param]]))
    }

    # If multiple values
    if (param_dim > 0) {
      # Get mean, sd, and CV
      sim_mean[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, mean)
      sim_median[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, median)
      sim_sd[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, sd)
      sim_cv[[param]] <- sim_sd[[param]]/sim_mean[[param]]

      sim_mre[[param]] <- apply(simplify2array(sim_mre[[param]]), 1:param_dim, median)
      sim_mse[[param]] <- apply(simplify2array(sim_mse[[param]]), 1:param_dim, mean)
    }
  }


  result_list <- list(Mean = sim_mean, Median = sim_median, SD = sim_sd, CV = sim_cv, MRE = sim_mre, MSE = sim_mse)
  return(result_list)
}


#' Generate Length-at-Age Transition Matrix
#'
#' This function calculates a probability transition matrix that defines the
#' probability of a fish of a given age belonging to specific length bins.
#' It supports Von Bertalanffy and Richards growth models and includes
#' a Stock Synthesis (SS) style plus-group correction.
#'
#' @param fracyr Numeric. Fraction of the year (0 = Jan 1st).
#' @param nsex_sp Integer. Number of sexes for the species.
#' @param nages_sp Integer. Number of age classes.
#' @param nlengths_sp Integer. Number of length bins.
#' @param nyrs Integer. Number of years in the simulation.
#' @param lengths_sp Vector. Boundaries of the length bins.
#' @param minage_sp Numeric. The reference age (L1) for growth estimation.
#' @param maxage_sp Numeric. The age at which growth enters the asymptotic phase.
#' @param growth_params_sp Array. Dimensions (sex, yr, 4).
#'   Params: K, L1, Linf, Richards m.
#' @param growth_log_sd_sp Array. Dimensions (sex, 2).
#'   Log-SD of length: 1st param is SD at minage, 2nd param is SD at maxage.
#' @param growth_model_sp Integer. 1 = Von Bertalanffy, 2 = Richards.
#'
#' @return A 4D array of probabilities with dimensions (sex, age, length, year).
get_growth_matrix_r <- function(fracyr, nsex_sp, nages_sp, nlengths_sp, nyrs,
                                lengths_sp, minage_sp, maxage_sp,
                                growth_params_sp, growth_log_sd_sp, growth_model_sp) {

  # Define names for the dimensions
  dim_names <- list(
    sex    = paste0("Sex_", 1:nsex_sp),
    age    = paste0("Age_", 1:nages_sp),
    length = paste0("Len_", lengths_sp),
    year   = paste0("Year_", 1:nyrs)
  )

  # Initialize Output: (sex, age, ln, yr)
  growth_matrix <- array(0, dim = c(nsex_sp, nages_sp, nlengths_sp, nyrs),
                         dimnames = dim_names)
  length_at_age <- array(0, dim = c(nsex_sp, nages_sp, nyrs),
                         dimnames =   list(
                           sex    = paste0("Sex_", 1:nsex_sp),
                           age    = paste0("Age_", 0:(nages_sp - 1)),
                           year   = paste0("Year_", 1:nyrs)
                         ))
  length_sd     <- array(0, dim = c(nsex_sp, nages_sp, nyrs))

  l_min <- lengths_sp[1]
  l_max <- lengths_sp[nlengths_sp]

  for(s in 1:nsex_sp) {
    for(y in 1:nyrs) {
      # --- 1. Calculate Mean Length at Age ---
      # Params: 1:K, 2:L1, 3:Linf, 4:Richards_m
      k    <- growth_params_sp[s, y, 1]
      l1   <- growth_params_sp[s, y, 2]
      linf <- growth_params_sp[s, y, 3]

      b_len <- (l1 - l_min) / minage_sp

      for(a in 1:nages_sp) {
        current_age <- a + fracyr

        if (growth_model_sp == 1) { # VB
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            if(y == 1){
              length_at_age[s, a, y] <- linf + (l1 - linf) * (exp(-k * (current_age - minage_sp)))
            } else{
              if(a == nages_sp){ # linear growth + growth equation
                last_linear = l_min + b_len * minage_sp # last age (cont) with linear growth

                length_at_age[s, a, y] = last_linear + (last_linear - linf) * (exp(-k * (current_age - minage_sp)) - 1.0)
              }else{
                length_at_age[s, a, y] <- length_at_age[s, a-1, y-1] + (length_at_age[s, a-1, y-1] - growth_params_sp[s, y-1, 3]) * (exp(-growth_params_sp[s, y-1, 1]) - 1)
              }
            }
          }
        } else if (growth_model_sp == 2) { # Richards
          m <- growth_params_sp[s, y, 4]
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            if(y == 1){
              length_at_age[s, a, y] <- (linf^m + (l1^m - linf^m) * (exp(-k * (current_age - minage_sp))))^(1/m)
            } else{
              if(a == nages_sp){ # linear growth + growth equation
                last_linear = l_min + b_len * minage_sp # last age (cont) with linear growth

                length_at_age[s, a, y] = (last_linear^m + (last_linear^m - linf^m) * (exp(-k * (current_age - minage_sp)) - 1.0))^(1/m)
              }else{
                lagk <- growth_params_sp[s, y-1, 2]
                lagm <- growth_params_sp[s, y-1, 4]
                laglinf <-  growth_params_sp[s, y-1, 3]
                length_at_age[s, a, y] <- (length_at_age[s, a-1, y-1]^lagm + (length_at_age[s, a-1, y-1]^lagm - laglinf^lagm) * (exp(-lagk) - 1))^1/lagm
              }
            }
          }
        }

        # --- 2. Plus Group Correction (SS Style) ---
        if(a == nages_sp) {
          diff <- growth_params_sp[s, y, 3] - length_at_age[s, a, y] # Linf - current size
          ages <- 0:(nages_sp)
          weight_a <- exp(-0.2 * ages)
          vals <- length_at_age[s, a, y] + (ages / nages_sp) * diff
          length_at_age[s, a, y] <- sum(vals * weight_a) / sum(weight_a)
        }

        # --- 3. SD Calculation ---
        sd1 <- exp(growth_log_sd_sp[s, 1])
        sda <- exp(growth_log_sd_sp[s, 2])

        if(current_age <= minage_sp) {
          length_sd[s, a, y] <- sd1
        } else if(a == nages_sp) {
          length_sd[s, a, y] <- sda
        } else {
          slope <- (sda - sd1) / (linf - l1) # Match C++ interpolation
          length_sd[s, a, y] <- sd1 + slope * (length_at_age[s, a, y] - l1)
        }

        # --- 4. Matrix Distribution ---
        for(l in 1:nlengths_sp) {
          if(l == 1) {
            fac1 <- (lengths_sp[2] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- stats::pnorm(fac1)
          } else if(l == nlengths_sp) {
            fac1 <- (lengths_sp[nlengths_sp] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- 1 - stats::pnorm(fac1)
          } else {
            fac1 <- (lengths_sp[l+1] - length_at_age[s, a, y]) / length_sd[s, a, y]
            fac2 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- stats::pnorm(fac1) - stats::pnorm(fac2)
          }
        }
      }
    }
  }

  return(list(length_at_age = length_at_age, growth_matrix = growth_matrix))
}


#' Calculate Predicted Weight-at-Age
#'
#' Converts a growth matrix (length-at-age probabilities) into mean weight-at-age
#' using a length-weight relationship (W = a * L^b).
#'
#' @param nsex_sp Integer. Number of sexes.
#' @param nages_sp Integer. Number of age classes.
#' @param nlengths_sp Integer. Number of length bins.
#' @param nyrs Integer. Number of years.
#' @param lengths_sp Vector. Boundaries of the length bins.
#' @param length_at_age Array. Mean length at age from get_growth_matrix_r.
#' @param growth_matrix Array. 4D array (sex, age, length, year) from get_growth_matrix_r.
#' @param lw_params Array. Dimensions (sex, yr, 2).
#'   Params: 1st is alpha (a), 2nd is beta (b).
#'
#' @details The function calculates midpoints for length bins to avoid bias.
#' For the first bin, it assumes the width is equal to the second bin's width.
#' The final weight-at-age is the expected value across all length bins for that age.
#'
#' @return A 3D array of mean weights with dimensions (sex, age, year).
get_weight_at_age_r <- function(nsex_sp, nages_sp, nlengths_sp, nyrs,
                                lengths_sp, length_at_age, growth_matrix, lw_params) {
  # Define names for the dimensions
  dim_names <- list(
    sex  = paste0("Sex_", 1:nsex_sp),
    age  = paste0("Age_", 1:nages_sp),
    year = paste0("Year_", 1:nyrs)
  )

  # Output: (sex, age, yr)
  waa <- array(0, dim = c(nsex_sp, nages_sp, nyrs),
               dimnames = dim_names)


  for(s in 1:nsex_sp) {
    for(y in 1:nyrs) {
      alpha <- lw_params[s, y, 1]
      beta  <- lw_params[s, y, 2]

      # Weight at length for all bins
      wal <- alpha * (lengths_sp + (lengths_sp[2] - lengths_sp[1])/2) ^ beta

      for(a in 1:nages_sp) {
        # Matrix multiply: Prob(length | age) * Weight(length)
        waa[s, a, y] <- sum(growth_matrix[s, a, , y] * wal)
      }
    }
  }
  return(waa)
}


