#' Forecast skill across retrospective peels
#'
#' Scores how well each retrospective peel's projection recovers the full
#' time-series model's estimate over the years the peel did not see, given the
#' known catch. The full model is taken as the best available estimate of what
#' happened, so this measures forecast skill rather than the estimation
#' consistency [retrospective()] reports through Mohn's rho.
#'
#' @details
#' Each peel drops the last `i` years of data, refits, and then solves fishing
#' mortality against the catch that was actually taken in those years -- so the
#' projection differs from the full model only in what it could not know:
#' the recruitment, and anything downstream of it. That makes this the natural
#' comparison for recruitment-projection assumptions, e.g.
#' `proj_mean_rec = TRUE` against `FALSE` against a DSEM.
#'
#' Skill is reported as the mean absolute scaled error (MASE), following the
#' hindcast cross-validation of Kell et al. (2016) but scored against the full
#' model's estimate rather than against held-out observations:
#'
#' \deqn{MASE = \frac{\overline{|forecast - reference|}}{\overline{|naive - reference|}}}
#'
#' The naive forecast is the peel's own terminal estimate carried forward,
#' i.e. "assume nothing changes". So:
#' \itemize{
#'   \item `MASE < 1` -- the projection beats assuming no change.
#'   \item `MASE > 1` -- it is worse than assuming no change, which usually means
#'     the recruitment assumption is pulling the projection the wrong way.
#'   \item `MASE = 1` -- no better than persistence.
#' }
#'
#' A MASE is undefined when the naive forecast happens to be exactly right
#' (a zero denominator); those rows come back `NA` rather than `Inf`.
#'
#' Read short horizons with care. At `years_ahead = 1` the denominator is a
#' single \eqn{|naive - reference|}, so a peel whose persistence forecast lands
#' close to the reference by chance produces a very large MASE that says more
#' about the denominator than about the projection -- values in the tens appear
#' routinely at one-year horizons and settle by three. Compare the deeper peels,
#' or the median across peels, rather than any single number. This is a property
#' of the scaling, not of the model, and it is shared with the `ss3diags`
#' implementation the definition matches.
#'
#' @param Rceattle A fitted Rceattle model (the full time series).
#' @param peels Number of peels. Passed to [retrospective()].
#' @param quantity Quantities to score against the full model. Any of `"ssb"`,
#'   `"biomass"`, `"R"`. Ignored when `reference = "observed"`.
#' @param reference What to score the projection against.
#'   `"model"` (default) compares each peel's projection to the full time-series
#'   model's estimate, taking that as the best available account of what
#'   happened. `"observed"` is classic hindcast cross-validation: it compares the
#'   peel's PREDICTED SURVEY INDEX to the index values actually observed in the
#'   held-out years. Both can be asked for. `"model"` scores the quantity a
#'   projection is used for (SSB, recruitment); `"observed"` scores the only
#'   thing that was really measured, and is the version comparable to MASE
#'   values published for other assessment platforms.
#' @param forecast_rec How the peeled years get their recruitment.
#'   `"model"` (default here) uses the model's own projection rule, in
#'   precedence order: `proj_mean_rec = TRUE` projects at mean recruitment,
#'   whatever process the model carries, because that is what the setting means;
#'   otherwise the LATENT STATES supply it where the deviations are random
#'   effects (`random_rec = TRUE`, or a DSEM), so an AR1's autocorrelation or a
#'   DSEM's lagged and covariate paths propagate into the forecast; otherwise
#'   recruitment comes off the stock-recruit curve. `"mean"` forces the
#'   historical mean for all of them, which is [retrospective()]'s default and
#'   the convention Mohn's rho is computed under.
#'
#'   The default differs from [retrospective()]'s on purpose. A peel's forecast
#'   years are hindcast years, so `proj_mean_rec` cannot reach them by itself,
#'   and with `"mean"` every model forecasts identically -- which makes this
#'   function unable to tell projection methods apart, the one thing it exists
#'   to do. Measured on the GOA arrowtooth model, `proj_mean_rec = TRUE` and
#'   `FALSE` returned byte-identical MASE under `"mean"`.
#' @param retro Optionally an already-computed [retrospective()] result, to
#'   avoid refitting when both are wanted. Note a `retro` computed with
#'   [retrospective()]'s defaults carries `forecast_rec = "mean"`.
#' @param ... Passed to [retrospective()] (`cores`, `getsd`, `rescale`).
#'
#' @return A list with
#'   \describe{
#'     \item{`mase`}{one row per peel x species x quantity, with `mase`,
#'       `mae_forecast`, `mae_naive` and the number of forecast years scored.}
#'     \item{`by_year`}{the underlying series: peel, species, quantity, year,
#'       years-ahead, `forecast`, `reference`, `naive`.}
#'   }
#'
#' @references
#' Kell, L.T., Kimoto, A., Kitakado, T. (2016) Evaluation of the prediction
#' skill of stock assessment using hindcasting. \emph{Fisheries Research} 183,
#' 119-127.
#'
#' @export
hindcast_skill <- function(Rceattle = NULL, peels = 5,
                           quantity = c("ssb", "biomass", "R"),
                           reference = c("model", "observed"),
                           forecast_rec = c("model", "mean"),
                           retro = NULL, ...) {
  if (!inherits(Rceattle, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }
  quantity     <- match.arg(quantity, several.ok = TRUE)
  reference    <- match.arg(reference, several.ok = TRUE)
  forecast_rec <- match.arg(forecast_rec)

  for (q in quantity) {
    if (is.null(Rceattle$quantities[[q]])) {
      stop("`quantity` '", q, "' is not reported by this fit.", call. = FALSE)
    }
  }

  if (is.null(retro)) {
    retro <- retrospective(Rceattle, peels = peels,
                           forecast_rec = forecast_rec, ...)
  }
  # retrospective() returns rev(c(list(Rceattle), peels)), so the LAST element
  # is the input model and the peels run deepest-first. Score only the peels;
  # comparing the parent to itself would give MASE 0 and mean nothing.
  ml <- retro$Rceattle_list
  if (length(ml) < 2L) {
    stop("No peel survived, so there is nothing to score. Inspect ",
         "retrospective(..., peels = 1) and read its $convergence.",
         call. = FALSE)
  }
  peel_models <- ml[-length(ml)]

  styr   <- Rceattle$data_list$styr
  endyr  <- Rceattle$data_list$endyr
  nspp   <- Rceattle$data_list$nspp
  spnames <- if (!is.null(Rceattle$data_list$spnames)) {
    Rceattle$data_list$spnames
  } else as.character(seq_len(nspp))

  by_year <- list()
  if ("model" %in% reference) for (m in peel_models) {
    ep <- m$data_list$endyr_peel
    # A peel with no forecast years (endyr_peel == endyr) scores nothing.
    if (is.null(ep) || ep >= endyr) next
    fyrs <- (ep + 1):endyr
    cols <- fyrs - styr + 1L
    last <- ep - styr + 1L

    for (q in quantity) {
      ref <- Rceattle$quantities[[q]]
      fc  <- m$quantities[[q]]
      if (is.null(fc) || ncol(fc) < max(cols)) next
      for (sp in seq_len(nspp)) {
        by_year[[length(by_year) + 1L]] <- data.frame(
          peel        = endyr - ep,
          species     = spnames[sp],
          quantity    = q,
          year        = fyrs,
          years_ahead = seq_along(fyrs),
          forecast    = as.numeric(fc[sp, cols]),
          reference   = as.numeric(ref[sp, cols]),
          # "Assume nothing changes": the peel's own terminal estimate held flat.
          naive       = as.numeric(fc[sp, last]),
          stringsAsFactors = FALSE)
      }
    }
  }
  # --- reference = "observed": classic hindcast cross-validation ----------
  # A peel's index_data is filtered to endyr_peel, so its index_hat has no rows
  # for the held-out years. Rebuild each peel at estimateMode = 3 -- build, do
  # not optimize -- against the FULL index_data, with the peel's own parameters
  # and map. That yields the peel's PREDICTION at observations it never fitted,
  # which is the whole point, without letting it re-estimate anything.
  if ("observed" %in% reference) {
    idx_full <- Rceattle$data_list$index_data
    for (m in peel_models) {
      ep <- m$data_list$endyr_peel
      if (is.null(ep) || ep >= endyr) next
      dl <- Rceattle$data_list
      dl$endyr_peel <- ep
      pred <- try(suppressWarnings(suppressMessages(.refit_like(
        data_list = dl, inits = m$estimated_params, map = m$map,
        estimateMode = 3, phase = FALSE, getsd = FALSE))), silent = TRUE)
      if (inherits(pred, "try-error") || is.null(pred$quantities$index_hat)) next
      ih <- as.numeric(pred$quantities$index_hat)
      if (length(ih) != nrow(idx_full)) next

      held <- which(idx_full$Year > ep & idx_full$Year <= endyr)
      for (flt in unique(idx_full$Fleet_code[held])) {
        rows <- held[idx_full$Fleet_code[held] == flt]
        if (!length(rows)) next
        # Naive = the last index value this fleet was observed at before the
        # peel, carried forward. Undefined if the fleet has no pre-peel
        # observation, in which case there is no persistence baseline.
        prior <- which(idx_full$Fleet_code == flt & idx_full$Year <= ep)
        if (!length(prior)) next
        naive <- idx_full$Observation[prior[which.max(idx_full$Year[prior])]]
        sp <- idx_full$Species[rows][1]
        by_year[[length(by_year) + 1L]] <- data.frame(
          peel        = endyr - ep,
          species     = if (is.numeric(sp) && sp >= 1 && sp <= nspp) spnames[sp] else as.character(sp),
          quantity    = paste0("index_fleet_", flt),
          year        = idx_full$Year[rows],
          years_ahead = idx_full$Year[rows] - ep,
          forecast    = ih[rows],
          reference   = idx_full$Observation[rows],
          naive       = as.numeric(naive),
          stringsAsFactors = FALSE)
      }
    }
  }

  if (!length(by_year)) {
    stop("No peel had forecast years to score.", call. = FALSE)
  }
  by_year <- do.call(rbind, by_year)

  key <- interaction(by_year$peel, by_year$species, by_year$quantity,
                     drop = TRUE, lex.order = TRUE)
  mase <- do.call(rbind, lapply(split(by_year, key), function(z) {
    mae_f <- mean(abs(z$forecast - z$reference))
    mae_n <- mean(abs(z$naive    - z$reference))
    data.frame(
      peel       = z$peel[1], species = z$species[1], quantity = z$quantity[1],
      n_years    = nrow(z),
      mae_forecast = mae_f,
      mae_naive    = mae_n,
      # NA rather than Inf when persistence was exactly right: an undefined
      # ratio is not infinitely bad skill, and Inf would poison any mean.
      mase       = if (mae_n > 0) mae_f / mae_n else NA_real_,
      stringsAsFactors = FALSE)
  }))
  rownames(mase) <- NULL
  mase <- mase[order(mase$quantity, mase$species, mase$peel), ]

  list(mase = mase, by_year = by_year)
}
