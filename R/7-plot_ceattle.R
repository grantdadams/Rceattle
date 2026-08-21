#' Make a vector of colors.
#'
#' A subset of rich.colors by Arni Magnusson from the gplots package, with the
#' addition of alpha transparency (which is now available in the gplots version
#' as well)
#'
#'
#' @param n Number of colors to generate.
#' @param alpha Alpha transparency value for all colors in vector. Value is
#' passed to rgb function.
#' @author Arni Magnusson, Ian Taylor
rich.colors.short <- function(n,alpha=1){
  x <- seq(0, 1, length = n)
  r <- 1/(1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
  rgb.m <- matrix(c(r, g, b), ncol = 3)
  apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha = alpha))
}

#' Plot time-series
#'
#' @description Function that plots the time-series (SSB/B/R/Depletion) 95% CI trends as estimated from Rceattle
#'
#' @details
#' # Units
#'
#' The model carries numbers-at-age in **thousands** and weight-at-age in
#' **kg**, so every biomass series (`biomass`, `ssb`, `exploitable_biomass`)
#' comes out of the model in **mt** and recruitment comes out in **thousands of
#' fish**. For display these are divided by 1e6 (million mt) and 1e3 (millions
#' of recruits) respectively; depletion is a ratio and is not rescaled. Supply
#' the model's inputs on that convention -- catch and index in mt, weight-at-age
#' in kg -- or the axis labels will not describe what is plotted.
#'
#' # Confidence intervals
#'
#' `add_ci = TRUE` draws a 95% interval. Every strictly positive series takes it
#' on the log scale as `exp(log(x) +/- qnorm(0.975) * sd_log)`. These quantities are
#' built multiplicatively, so `log(x)` is close to linear in the estimated
#' parameters where `x` is exponential in them: the interval is both a better
#' linearization and right-skewed, and it cannot cross zero the way a symmetric
#' natural-scale interval does for weak year classes and depleted stocks.
#'
#' `sd_log` comes from the model's own `log_biomass` / `log_ssb` / `log_R` where
#' those are reported, and is otherwise recovered as `sd(x) / x` -- the delta
#' method's own identity, which matches the reported values to machine
#' precision. That covers `exploitable_biomass` and the two depletions, which
#' cannot be reported on the log scale (`exploitable_biomass` is identically 0
#' without projection F), and models fit before the `log_*` series existed. A
#' non-positive or unreported value keeps the symmetric interval.
#'
#' @inheritParams rceattle-plot-args
#' @param output derived quantity of interest: recruitment, biomass, ssb, depletion, or ssb_depletion. Uses same name as ".cpp" file.
#' @param ylab Y-axis label. `NULL` (default) derives one from `output` and the
#'   model's `minage`.
#' @param save Write the plotted series to CSV alongside the figure.
#' @param mse Is an MSE object from \code{\link{load_mse}} or \code{\link{run_mse}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#' @param mod_avg TRUE/FALSE
#' @param zero_y Anchor the y-axis at zero. TRUE for the absolute series
#'   (biomass, SSB, recruitment), where a truncated axis exaggerates change;
#'   FALSE for the depletions, which are already on a relative scale.
#' @param reference Reference model, drawn in black at 1.5x `lwd`.
#'
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_timeseries <- function(Rceattle,
                            output = "biomass",
                            ylab = NULL,
                            file = NULL,
                            model_names = NULL,
                            line_col = NULL,
                            species = NULL,
                            spnames = NULL,
                            add_ci = FALSE,
                            lwd = 3,
                            save = FALSE,
                            legend.pos = "topright",
                            right_adj = 0,
                            width = 7,
                            height = 6.5,
                            minyr = NULL,
                            maxyr = NULL,
                            incl_proj = FALSE,
                            mod_cex = 1,
                            lty = rep(1, length(Rceattle)),
                            alpha = 0.4,
                            mse = FALSE,
                            OM = TRUE,
                            reference = NULL,
                            zero_y = FALSE,
                            mod_avg = rep(FALSE, length(Rceattle))) {

  ## Model object manipulation ----
  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    } else {
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    nmse = length(Rceattle)
    add_ci = TRUE
    incl_proj = TRUE
  }


  # Convert single one into a list
  if(inherits(Rceattle, "Rceattle")){
    Rceattle <- list(Rceattle)
  }


  ## Line width / colour ----
  lwd <- rep(lwd, length(Rceattle))
  # line_col stays NULL unless the user supplies colours, in which case it falls
  # through to the default colourblind-safe scale (.rceattle_scale) at render
  # time. When supplied, recycle to one colour per model.
  if (!is.null(line_col)) {
    line_col <- rep(line_col, length.out = length(Rceattle))
  }

  ## Add reference model
  if(!is.null(reference)){
    Rceattle <- c(Rceattle, list(reference))
    mod_avg = c(mod_avg, FALSE)
    lty = c(lty, 1)
    if (!is.null(line_col)) line_col <- c(line_col, "black")
    lwd <- c(lwd, lwd[1] * 1.5)
    # Give the reference its own label so .model_labels() does not recycle the
    # user's model_names (which would collapse the reference into the first
    # model's factor level and drop its distinct colour / width / type).
    if (!is.null(model_names)) model_names <- c(model_names, "Reference")
  }

  # Extract years
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }
  if(!is.null(maxyr)){
    years <- lapply(Rceattle, function(x) x$data_list$styr:min(c(maxyr, x$data_list$projyr)))
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(years, length)
  nyrs <- max(nyrs_vec)
  if(is.null(minyr)){minyr <- min((sapply(years, min)))}
  if(is.null(maxyr)){maxyr <- max((sapply(years, max)))}

  # Extract species information
  nspp <- Rceattle[[1]]$data_list$nspp
  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics

  # `species` selects (indices, names, a logical mask, or "all"); `spnames`
  # labels. Shared with every other plotter so the pair means one thing.
  sp_sel  <- .resolve_species(Rceattle, species, spnames)
  spp     <- sp_sel$index
  spnames <- sp_sel$spnames

  # Built here rather than in .ts_wrapper() because the age it names is
  # per-species data. An explicit `ylab` always wins.
  if (is.null(ylab)) ylab <- .rce_ts_ylab(output, minage)


  ## Get quantity of interest ----
  # Save objects
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  # Delta-method SD of log(quantity), when the model ADREPORTed a matching
  # `log_<output>` series. Distinct from log_quantity_sd below, which is the
  # across-sample SD used by the model-averaging path.
  quantity_log_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_mu <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  ptarget = matrix(NA, nrow = length(Rceattle), ncol = nspp)
  plimit = matrix(NA, nrow = length(Rceattle), ncol = nspp)


  for (i in 1:length(Rceattle)) {
    # - Get quantities
    ptarget[i,] <- Rceattle[[i]]$data_list$Ptarget
    plimit[i,] <- Rceattle[[i]]$data_list$Plimit
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities[[output]][,1:nyrs_vec[i]]

    # Get SD of quantity
    if(add_ci & !mse){
      # The series is ADREPORTed over the whole modelled period, flattened
      # column-major with species varying fastest, so the first
      # nspp * nyrs_hindcast values are exactly the years being plotted.
      # n_total states that full shape, so the slice is checked rather than
      # assumed -- an unexpected block length draws no interval instead of one
      # belonging to different cells.
      n_need  <- nyrs_vec[i] * nspp
      n_total <- ncol(Rceattle[[i]]$quantities[[output]]) * nspp
      sd_vals <- .rce_series_sd(Rceattle[[i]], output, n_need, n_total)
      if (is.null(sd_vals)) {
        # Leave the SD array as NA and draw no ribbon. Only warn for a series
        # the model could have reported: F_spp and the other REPORT()-only
        # quantities never carry standard errors, so warning on them would fire
        # on every fit.
        .rce_no_ci(TRUE, paste0("`", output, "` in model ", i),
                   paste("this fit carries no standard errors for it",
                         "(not reported by the model, or getsd = FALSE)"),
                   warn = .rce_quantity_adreport(output))
      } else {
        quantity_sd[, 1:nyrs_vec[i], i] <-
          replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_vals)

        # Log-scale SD where the model reported one; the rest are recovered
        # from sd/x below.
        log_vals <- .rce_series_sd(Rceattle[[i]], paste0("log_", output),
                                   n_need, n_total)
        if (!is.null(log_vals)) {
          quantity_log_sd[, 1:nyrs_vec[i], i] <-
            replace(quantity_log_sd[, 1:nyrs_vec[i], i], values = log_vals)
        }
      }
    }

    # - Model average
    if(mod_avg[i]){
      log_quantity_sd[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples[quantity][,1:nyrs_vec[i],], c(1,2), function(x) sd(as.vector(log(x))))
      log_quantity_mu[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples[quantity][,1:nyrs_vec[i],], c(1,2), function(x) mean(as.vector(log(x))))
    }
  }

  ## Get confidence intervals ----
  if(!mse){
    # - Single model. .rce_ci_bounds() is shared with as.data.frame.Rceattle(),
    #   so the figure and the table cannot disagree: a strictly positive series
    #   takes exp(log(x) +/- z * sd_log), recovering sd(log x) = sd(x)/x where
    #   the model reported no log_* row, and everything else stays symmetric.
    ci95 <- .rce_ci_bounds(quantity, quantity_sd, quantity_log_sd,
                           z = stats::qnorm(0.975))
    ci50 <- .rce_ci_bounds(quantity, quantity_sd, quantity_log_sd,
                           z = stats::qnorm(0.75))

    quantity_upper95 <- array(ci95$upr, dim = dim(quantity))
    quantity_lower95 <- array(ci95$lwr, dim = dim(quantity))
    quantity_upper50 <- array(ci50$upr, dim = dim(quantity))
    quantity_lower50 <- array(ci50$lwr, dim = dim(quantity))
  } else {
    # - MSE objects
    ptarget <- ptarget[1,]
    plimit <- plimit[1,]

    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.975) )
    quantity_lower95 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.025) )
    quantity_upper50 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.75) )
    quantity_lower50 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.25) )

    # -- Put back in array for indexing below
    if(is.null(reference)){
      quantity <- array(apply( quantity[,,1:nmse], c(1,2), mean ), dim = c(nspp, nyrs,  1))
      quantity_upper95 <- array(quantity_upper95, dim = c(nspp, nyrs,  1))
      quantity_lower95 <- array(quantity_lower95, dim = c(nspp, nyrs,  1))
      quantity_upper50 <- array(quantity_upper50, dim = c(nspp, nyrs,  1))
      quantity_lower50<- array(quantity_lower50, dim = c(nspp, nyrs,  1))
    } else {
      quantity_upper95 <- array(c(quantity_upper95, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity_lower95 <- array(c(quantity_lower95, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity_upper50 <- array(c(quantity_upper50, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity_lower50<- array(c(quantity_lower50, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity <- array(c(apply( quantity[,,1:nmse], c(1,2), mean ), quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
    }
  }

  # - Model Average
  for (i in 1:length(Rceattle)) {
    if(mod_avg[i]){
      quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
    }
  }

  ## Rescale ----
  # Numbers-at-age are in thousands and weight-at-age in kg, so biomass is in mt
  # and recruitment in thousands of fish. Depletion is a ratio and is not
  # rescaled.
  rescale <- .RCE_TS_RESCALE[output]
  if (!is.na(rescale)) {
    quantity <- quantity / rescale
    quantity_upper95 <- quantity_upper95 / rescale
    quantity_lower95 <- quantity_lower95 / rescale
    quantity_upper50 <- quantity_upper50 / rescale
    quantity_lower50 <- quantity_lower50 / rescale
  }


  ## Save ----
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      if (ncol(dat) > 1) {
        for (j in 2:ncol(dat)) {
          dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
          colnames(dat_new2) <- rep(model_names[j], 3)
          dat_new <- cbind(dat_new, dat_new2)
        }
      }

      write.csv(dat_new, file = paste0(file, "_", output,"_species_", i, ".csv"))
    }
  }


  ## Assemble tidy data ----
  # Build one row per (species, year, model) holding exactly the plotted
  # quantities (the same arrays, reshaped) so the returned ggplot's `$data` can
  # be tested against the model's source quantities.
  model_names_use <- .model_labels(Rceattle, model_names)
  df_list <- list()
  for (k in seq_len(dim(quantity)[3])) {
    yk  <- years[[k]]
    nyk <- length(yk)
    for (sp in spp) {
      df_list[[length(df_list) + 1L]] <- data.frame(
        Species = spnames[sp],
        Year    = yk,
        Model   = model_names_use[k],
        value   = quantity[sp, seq_len(nyk), k],
        lower95 = quantity_lower95[sp, seq_len(nyk), k],
        upper95 = quantity_upper95[sp, seq_len(nyk), k],
        lower50 = quantity_lower50[sp, seq_len(nyk), k],
        upper50 = quantity_upper50[sp, seq_len(nyk), k],
        show_ci = (estDynamics[sp] == 0),
        stringsAsFactors = FALSE
      )
    }
  }
  plot_df <- do.call(rbind, df_list)
  plot_df$Species <- factor(plot_df$Species, levels = spnames[spp])
  plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))

  ## Render ----
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$Year, y = .data$value,
                 colour = .data$Model, fill = .data$Model))

  # Credible-interval ribbons (only for estimated species)
  if (add_ci) {
    ci_df <- plot_df[plot_df$show_ci, , drop = FALSE]
    if (nrow(ci_df) > 0) {
      p <- p + ggplot2::geom_ribbon(
        data = ci_df,
        ggplot2::aes(ymin = .data$lower95, ymax = .data$upper95),
        alpha = alpha / 2, colour = NA)
      if (mse) {
        p <- p + ggplot2::geom_ribbon(
          data = ci_df,
          ggplot2::aes(ymin = .data$lower50, ymax = .data$upper50),
          alpha = alpha, colour = NA)
      }
    }
  }

  # Line width / type: honour per-model lwd / lty, keyed to Model when they
  # vary so the uniform default adds no extra legend. The lwd / 3 bridge lives
  # in the helper, which every other plotter now shares.
  line_p <- .rce_line_params(lwd = lwd, lty = lty,
                             lwd_by = "Model", lty_by = "Model")
  p <- .rce_add_line(p, line_p)

  # Projection divider at the terminal hindcast year
  p <- p + .rce_proj_divider(Rceattle, incl_proj)

  # Depletion reference points (per species)
  if (output %in% "ssb_depletion") {
    ref_df <- data.frame(
      Species = factor(spnames[spp], levels = spnames[spp]),
      target  = ptarget[spp],
      limit   = plimit[spp])
    p <- p +
      ggplot2::geom_hline(data = ref_df, inherit.aes = FALSE,
                          ggplot2::aes(yintercept = .data$target),
                          colour = "blue", linetype = 2) +
      ggplot2::geom_hline(data = ref_df, inherit.aes = FALSE,
                          ggplot2::aes(yintercept = .data$limit),
                          colour = "red", linetype = 2)
  }

  # A single-species model needs no species strip: the facet label would name
  # the only thing on the plot. Faceting is what carries that label, so drop it
  # rather than blanking the strip and leaving its gap.
  if (length(spp) > 1L) {
    p <- p + ggplot2::facet_wrap(~ Species, scales = "free_y", ncol = 1,
                                 strip.position = "top")
  }
  p <- p +
    .rce_year_limits(minyr, maxyr) +
    ggplot2::labs(x = "Year", y = ylab) +
    .rceattle_theme()
  # A quantity on an absolute scale is read against zero, so an axis starting at
  # the series minimum overstates how much it moves. Applied per panel, so
  # `scales = "free_y"` still fits each species.
  if (isTRUE(zero_y)) p <- p + ggplot2::expand_limits(y = 0)
  # `line_col = NULL` keeps the default Okabe-Ito palette; supplied colours are
  # applied to the model levels in plotting order.
  p <- .rceattle_scale(p, line_col = line_col)

  # A single model needs no colour legend
  if (nlevels(plot_df$Model) < 2L) {
    p <- p + ggplot2::guides(colour = "none", fill = "none",
                             linewidth = "none", linetype = "none")
  }

  .save_ggplot(p, file = file, suffix = paste0(output, "_trajectory"),
               width = width, height = height)
}


# Display divisor per derived series; a series absent here is plotted unscaled.
# Any entry added needs its matching unit in .rce_ts_ylab(): the two drifting
# apart is how recruitment came to be plotted in thousands under a "million"
# label.
.RCE_TS_RESCALE <- c(
  biomass             = 1e6,   # mt              -> million mt
  ssb                 = 1e6,   # mt              -> million mt
  exploitable_biomass = 1e6,   # mt              -> million mt
  R                   = 1e3    # thousands of fish -> millions of recruits
)


# Default y-axis label for a derived series. The minimum age is per-species
# data, so it is interpolated at plot time; where species disagree on minage the
# age prefix is dropped, since one shared y axis cannot name two ages.
.rce_ts_ylab <- function(output, minage) {
  # Only recruitment names an age -- it is an age class. SSB is mature females,
  # not the minage+ stock, so its old "Age-1+" prefix was wrong as well as noisy.
  ages <- unique(as.integer(minage[!is.na(minage)]))
  age  <- if (length(ages) == 1L) paste0("Age-", ages, " ") else ""
  switch(output,
         biomass             = "Biomass (million mt)",
         ssb                 = "SSB (million mt)",
         exploitable_biomass = "Max exploitable biomass (million mt)",
         R                   = paste0(age, "recruits (millions)"),
         ssb_depletion       = "SSB depletion",
         biomass_depletion   = "Biomass depletion",
         output)
}


# Build a plot_timeseries() wrapper that pins the derived series (`output`) for
# one quantity while exposing plot_timeseries()'s full argument list unchanged.
# The six timeseries plotters below differ only in that string, so they share one
# body through this factory. The y-axis label is resolved inside
# plot_timeseries() from .rce_ts_ylab(), which needs the model's minage.
.ts_wrapper <- function(output, zero_y = FALSE) {
  force(output); force(zero_y)
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           spnames = NULL,
           add_ci = FALSE,
           lwd = 3,
           save = FALSE,
           right_adj = 0,
           legend.pos = "topright",
           width = 7,
           height = 6.5,
           minyr = NULL,
           maxyr = NULL,
           incl_proj = FALSE,
           mod_cex = 1,
           lty = rep(1, length(Rceattle)),
           alpha = 0.4,
           mod_avg = rep(FALSE, length(Rceattle)),
           mse = FALSE,
           OM = TRUE,
           reference = NULL,
           ylab = NULL) {

    plot_timeseries(Rceattle,
                    output = output,
                    ylab = ylab,
                    zero_y = zero_y,
                    file = file,
                    model_names = model_names,
                    line_col = line_col,
                    species = species,
                    spnames = spnames,
                    add_ci = add_ci,
                    lwd = lwd,
                    save = save,
                    right_adj = right_adj,
                    legend.pos = legend.pos,
                    width = width,
                    height = height,
                    minyr = minyr,
                    maxyr = maxyr,
                    incl_proj = incl_proj,
                    mod_cex = mod_cex,
                    lty = lty,
                    alpha = alpha,
                    mod_avg = mod_avg,
                    mse = mse,
                    OM = OM,
                    reference = reference)
  }
}

#' Plot biomass
#'
#' @description Plots the mean minage+ biomass (million mt) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the biomass trajectory.
plot_biomass <- .ts_wrapper("biomass", zero_y = TRUE)


#' Plot recruitment
#'
#' @description Plots the mean minage recruitment (millions) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the recruitment trajectory.
plot_recruitment <- .ts_wrapper("R", zero_y = TRUE)


#' Plot spawning stock biomass (SSB)
#'
#' @description Plots the mean SSB (million mt) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the SSB trajectory.
plot_ssb <- .ts_wrapper("ssb", zero_y = TRUE)


#' Plot exploitable biomass
#'
#' @description Plots the mean exploitable biomass (million mt) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the exploitable biomass trajectory.
plot_exploitable_biomass <- .ts_wrapper("exploitable_biomass", zero_y = TRUE)

#' Plot SSB depletion
#'
#' @description Plots the mean SSB depletion and 95\% CI trends as estimated from Rceattle.
#'   Depletion reference lines for Ptarget and Plimit are drawn in blue and red respectively.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the SSB depletion trajectory.
plot_depletionSSB <- .ts_wrapper("ssb_depletion", zero_y = TRUE)

#' Plot SSB depletion (deprecated name)
#'
#' @description Deprecated alias for \code{\link{plot_depletionSSB}}. Please use
#'   \code{plot_depletionSSB()} instead.
#'
#' @inheritParams plot_timeseries
#' @export
plot_ssb_depletion <- plot_depletionSSB

#' Plot biomass depletion
#'
#' @description Plots the mean biomass depletion and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the biomass depletion trajectory.
plot_depletion <- .ts_wrapper("biomass_depletion", zero_y = TRUE)


#' Plot selectivity
#'
#' @description Function that plots the fishery and survey selectivity as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param width Figure width in inches
#' @param height Figure height in inches
#'
#' @export
plot_selectivity <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           width = 7,
           height = 6.5,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3) {

    Rceattle <- .as_model_list(Rceattle)
    fc     <- Rceattle[[1]]$data_list$fleet_control
    nflt   <- nrow(fc)
    nages  <- Rceattle[[1]]$data_list$nages
    minage <- Rceattle[[1]]$data_list$minage
    nsex   <- Rceattle[[1]]$data_list$nsex
    styr   <- Rceattle[[1]]$data_list$styr
    endyr  <- Rceattle[[1]]$data_list$endyr
    hindyears <- styr:endyr

    # Selectivity-at-age over hindcast years for the first model:
    # sel_at_age[fleet, sex, age, year] (year index 1 == styr).
    sel <- Rceattle[[1]]$quantities$sel_at_age

    df_list <- list()
    for (flt in seq_len(nflt)) {
      sp <- fc$Species[fc$Fleet_code == flt]
      if (length(sp) == 0) next
      for (sex in seq_len(nsex[sp])) {
        sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
        for (yr in seq_along(hindyears)) {
          df_list[[length(df_list) + 1L]] <- data.frame(
            Fleet = as.character(fc$Fleet_name[fc$Fleet_code == flt]),
            Sex   = sex_lab,
            Age   = seq_len(nages[sp]) - 1 + minage[sp],
            Year  = hindyears[yr],
            Selectivity = as.numeric(sel[flt, sex, seq_len(nages[sp]), yr]),
            stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Fleet <- factor(plot_df$Fleet, levels = unique(fc$Fleet_name))

    # Year-coloured selectivity-at-age curves, faceted by fleet. Time-varying
    # selectivity shows as a colour fan; time-invariant collapses to one line.
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Age, y = .data$Selectivity,
                   group = interaction(.data$Year, .data$Sex),
                   colour = .data$Year, linetype = .data$Sex)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ Fleet) +
      ggplot2::scale_colour_viridis_c("Year") +
      ggplot2::labs(x = "Age", y = "Selectivity") +
      .rceattle_theme()
    if (length(unique(plot_df$Sex)) < 2L) {
      p <- p + ggplot2::guides(linetype = "none")
    }

    .save_ggplot(p, file = file, suffix = "selectivity",
                 width = width, height = height)
  }




#' Plot functional form
#'
#' @description Function to plot the functional form estimated or specified by \code{Rceattle}
#'
#' @param params Parameter list object from \code{\link{build_params}} or \code{Rceattle}
#' @param pred Predator index
#' @param pred_age Predator age
#' @param prey Prey index
#' @param msmMode Multispecies mode integer specifying functional form
#' @export
plot_form <- function( params = NULL, pred = 1, pred_age = 1, prey = 1, msmMode = 3){
  # The Kinzey & Punt (2009) functional-response parameters are commented out
  # in the TMB template, build_params(), build_map() and build_bounds(), and
  # data_check() blocks msmMode 3:9. Body retained below so this can be revived
  # alongside the C++.
  stop(paste0(
    "plot_form() plots the Kinzey & Punt (2009) functional responses ",
    "(msmMode 3-9), which are not available in this version of Rceattle: ",
    "the functional-response parameters are not estimated, and data_check() ",
    "blocks these modes. Use msmMode = 1 (Holsman et al. 2015 MSVPA) or ",
    "2 (Holling Type III MSVPA)."), call. = FALSE)

  # ---- retained for future development ----
  #
  #   # Get indices
  #   rsp = pred
  #   ksp = prey
  #
  #   # Get parameter values
  #   H_1 <- exp(params$logH_1)
  #   H_1a <- exp(params$logH_1a)
  #   H_1b <- exp(params$logH_1b)
  #
  #   H_2 <- exp(params$logH_2)
  #   H_3 <- exp(params$logH_3)
  #
  #   H_4 <- params$H_4
  #
  #   # Set up ratios
  #   Pred_r <- seq(from = 0.001, to = 5, length.out = 100) # Pred biomass relative to equilibrium
  #   Prey_r <- seq(from = 0.001, to = 5, length.out = 100) # Prey biomass relative to equilibrium
  #
  #   # Calculate functional form
  #   Term = H_1[rsp, ksp] * (1 + H_1a[rsp] * H_1b[rsp] / (pred_age + H_1b[rsp]))
  #   response <- matrix(NA, ncol = length(Prey_r), nrow = length(Pred_r))
  #   rownames(response) <- Pred_r
  #   colnames(response) <- Prey_r
  #
  #   for(i in 1:length(Pred_r)){
  #     for(j in 1:length(Prey_r)){
  #       response[i, j] <- Prey_r[j] * switch (
  #         as.character(msmMode),
  #         "2" = { # Holling Type I (linear)
  #           Term},
  #         "3" = { #  Holling Type II
  #           Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] )
  #         },
  #         "4" = { #  Holling Type III
  #           Term * (1 + H_2[rsp, ksp]) * ((Prey_r[j] ) ^ H_4[rsp, ksp]) / (1 + H_2[rsp, ksp] * ((Prey_r[j] ) ^ H_4[rsp, ksp])  )},
  #         "5" = { #  predator interference
  #           Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] + H_3[rsp, ksp] * (Pred_r[i] - 1) )
  #         },
  #         "6" = { # predator preemption
  #           Term * (1 + H_2[rsp, ksp] ) / ( (1 + H_2[rsp, ksp] * Prey_r[j]) * (1 + H_3[rsp, ksp] * (Pred_r[i] - 1)) )
  #         },
  #         "7" = { # Hassell-Varley
  #           Term * (2 + H_2[rsp, ksp] ) / (1 + H_2[rsp, ksp] * Prey_r[j] + ((Prey_r[j] ) ^ H_4[rsp, ksp]))
  #         },
  #         "8" = { #  Ecosim
  #           Term / (1 + H_3[rsp, ksp] * (Pred_r[i] - 1 ))},
  #         {
  #           stop("msmMode not implemented: ", msmMode)
  #         }
  #       )
  #     }
  #   }
  #
  #   df <- reshape2::melt(response)
  #   colnames(df) <- c("Pred_r", "Prey_r", "Response")
  #   df$Pred_r <- as.numeric(as.character(df$Pred_r))
  #   df$Prey_r <- as.numeric(as.character(df$Prey_r))
  #
  #   if (msmMode %in% c(2, 3, 4, 7)) {
  #     # functional response depends on the prey ratio only
  #     d <- stats::aggregate(Response ~ Prey_r, data = df, FUN = mean)
  #     p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Prey_r, y = .data$Response)) +
  #       ggplot2::geom_line(linewidth = 1) +
  #       ggplot2::labs(x = "Prey ratio", y = "Functional response")
  #   } else if (msmMode == 8) {
  #     d <- stats::aggregate(Response ~ Pred_r, data = df, FUN = mean)
  #     p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Pred_r, y = .data$Response)) +
  #       ggplot2::geom_line(linewidth = 1) +
  #       ggplot2::labs(x = "Predator ratio", y = "Functional response")
  #   } else if (msmMode %in% c(5, 6)) {
  #     p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Prey_r, y = .data$Pred_r,
  #                                           fill = .data$Response)) +
  #       ggplot2::geom_raster() +
  #       ggplot2::scale_fill_viridis_c("Functional\nresponse") +
  #       ggplot2::labs(x = "Prey ratio", y = "Predator ratio")
  #   } else {
  #     stop("msmMode not implemented: ", msmMode)
  #   }
  #   p + .rceattle_theme()
}





#' Plot M1 + M2
#'
#' @description Function that plots the M1 and M2 as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param incl_proj Include the projection years (TRUE/FALSE)
#' @param zlim zlim for M1 + M2 plots. Character - use max range across species in model. NULL - use species specific ranges. Vector of two.
#' @param type 0 = Tiles, 1 = contour, 2 = facet lines, 3 = persp
#' @param width Plot width when saved "inches"
#' @param height Plot height when saved "inches"
#' @param title Additional title to add. Will also add species names if not NULL
#' @param title_cex Font size for title
#' @param species Species to plot. Plots all if null.
#' @param log TRUE/FALSE use log M1 + M2
#' @param minyr First year to plot
#' @param theta theta for persp plot
#' @param maxage Plot up to this age. Plots all ages if NULL
#' @param M2 TRUE/FALSE Use M2 only (True) or total M (False)
#'
#' @export
plot_mortality <-
  function(Rceattle,
           file = NULL,
           incl_proj = FALSE,
           zlim = NULL,
           type = "lines",
           width = 8,
           height = 5.5,
           title = NULL,
           log = FALSE,
           minyr = NULL,
           theta = 155,
           species = NULL,
           maxage = NULL,
           title_cex = 10,
           M2 = TRUE) {

    Rceattle <- .as_model_list(Rceattle)
    if (length(Rceattle) > 1) stop("Can only plot one model")
    dl <- Rceattle[[1]]$data_list

    if (is.null(minyr)) minyr <- dl$styr
    endyr <- if (incl_proj) dl$projyr else dl$endyr
    years <- minyr:endyr
    nyrs  <- length(years)

    nspp        <- dl$nspp
    spnames     <- dl$spnames
    estdynamics <- dl$estDynamics
    nages       <- dl$nages
    if (!is.null(maxage)) nages <- pmin(nages, maxage)
    minage <- dl$minage
    nsex   <- dl$nsex
    if (is.null(species)) species <- seq_len(nspp)
    spp <- intersect(species, seq_len(nspp))

    # M-at-age over the selected years: M1 only (M2 = FALSE) or M1 + M2.
    yr_idx <- seq_len(nyrs) + (minyr - dl$styr)
    Msrc <- if (M2) Rceattle[[1]]$quantities$M2_at_age
            else    Rceattle[[1]]$quantities$M1_at_age

    df_list <- list()
    for (sp in spp) {
      if (estdynamics[sp] != 0) next
      ages <- seq_len(nages[sp]) - 1 + minage[sp]
      for (sex in seq_len(nsex[sp])) {
        sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
        for (yi in seq_along(years)) {
          val <- as.numeric(Msrc[sp, sex, seq_len(nages[sp]), yr_idx[yi]])
          df_list[[length(df_list) + 1L]] <- data.frame(
            Species = spnames[sp], Sex = sex_lab,
            Age = ages, Year = years[yi],
            M = if (log) log(val) else val,
            stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Species <- factor(plot_df$Species, levels = spnames[spp])
    # Msrc is M2_at_age when M2 = TRUE, else M1_at_age; label matches the source.
    ylab <- if (M2) "M2" else "M1"
    if (log) ylab <- paste0("log(", ylab, ")")

    if (identical(type, "heatmap") || identical(type, 0)) {
      # Age x year heatmap, faceted by species
      p <- ggplot2::ggplot(plot_df,
                           ggplot2::aes(x = .data$Year, y = .data$Age,
                                        fill = .data$M)) +
        ggplot2::geom_tile() +
        .facet_species(plot_df) +
        ggplot2::scale_fill_viridis_c(ylab) +
        ggplot2::labs(x = "Year", y = "Age") +
        .rceattle_theme()
    } else {
      # Year-coloured M-at-age curves (standard), faceted by species
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$Age, y = .data$M,
                     group = interaction(.data$Year, .data$Sex),
                     colour = .data$Year, linetype = .data$Sex)) +
        ggplot2::geom_line() +
        .facet_species(plot_df, scales = "free_y") +
        ggplot2::scale_colour_viridis_c("Year") +
        ggplot2::labs(x = "Age", y = ylab) +
        .rceattle_theme()
      if (length(unique(plot_df$Sex)) < 2L) {
        p <- p + ggplot2::guides(linetype = "none")
      }
    }
    if (!is.null(title)) p <- p + ggplot2::ggtitle(title)

    .save_ggplot(p, file = file, suffix = "mortality_at_age",
                 width = width, height = height)
  }







#' Plot maturity
#'
#' @description Function that plots the maturity of each species
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param width Figure width in inches
#' @param height Figure height in inches
#'
#' @export
plot_maturity <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           width = 4,
           height = 5.5,
           lwd = 3) {

    Rceattle <- .as_model_list(Rceattle)
    model_names_use <- .model_labels(Rceattle, model_names)
    nspp  <- Rceattle[[1]]$data_list$nspp
    nages <- Rceattle[[1]]$data_list$nages
    if (is.null(species)) species <- Rceattle[[1]]$data_list$spnames

    # Tidy maturity-at-age (data_list$maturity, species column dropped)
    df_list <- list()
    for (k in seq_along(Rceattle)) {
      mat <- Rceattle[[k]]$data_list$maturity[, -1, drop = FALSE]
      for (sp in seq_len(nspp)) {
        df_list[[length(df_list) + 1L]] <- data.frame(
          Model    = model_names_use[k],
          Species  = species[sp],
          Age      = seq_len(nages[sp]),
          Maturity = as.numeric(mat[sp, seq_len(nages[sp])]),
          stringsAsFactors = FALSE)
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
    plot_df$Species <- factor(plot_df$Species, levels = species)

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$Age, y = .data$Maturity,
                                      colour = .data$Model)) +
      ggplot2::geom_line(linewidth = 1) +
      .facet_species(plot_df, ncol = 1) +
      ggplot2::coord_cartesian(ylim = c(0, 1.1)) +
      ggplot2::labs(x = "Age", y = "Maturity")
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
    if (nlevels(plot_df$Model) < 2L) p <- p + ggplot2::guides(colour = "none")

    .save_ggplot(p, file = file, suffix = "maturity",
                 width = width, height = height)
  }


#' Plot biomass eaten
#'
#' @description Function that plots the biomass consumed trends as estimated from Rceattle. Returns and saves a figure with the biomass eaten trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param save Save figure to file?
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param minyr first year to plot
#' @param incl_proj TRUE/FALSE include projections years
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#' @param mod_cex Cex of text for model name legend
#' @param alpha Shading for confidence intervals
#' @param mod_avg TRUE/FALSE
#' @param mse Is an MSE object from \code{\link{load_mse}} or \code{\link{run_mse}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#'
#' @export
#'
plot_b_eaten <- function(Rceattle,
                         file = NULL,
                         model_names = NULL,
                         line_col = NULL,
                         species = NULL,
                         spnames = NULL,
                         add_ci = FALSE,
                         lwd = 3,
                         save = FALSE,
                         right_adj = 0,
                         width = 7,
                         height = 6.5,
                         minyr = NULL,
                         incl_proj = FALSE,
                         mod_cex = 1,
                         alpha = 0.4,
                         mod_avg = rep(FALSE, length(Rceattle)),
                         mse = FALSE,
                         OM = TRUE) {

  models <- .as_model_list(Rceattle, mse = mse, OM = OM)
  if (mse) incl_proj <- TRUE
  model_names_use <- .model_labels(models, model_names)
  if (is.null(spnames)) spnames <- models[[1]]$data_list$spnames
  nspp <- models[[1]]$data_list$nspp
  if (is.null(species)) species <- seq_len(nspp)

  # Total biomass eaten as prey per species/year: sum of B_eaten_as_prey
  # over sex and age.
  df_list <- list()
  for (k in seq_along(models)) {
    dl  <- models[[k]]$data_list
    yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
    be  <- models[[k]]$quantities$B_eaten_as_prey
    tot <- apply(be[, , , seq_along(yrs), drop = FALSE], c(1, 4), sum)
    for (sp in species) {
      df_list[[length(df_list) + 1L]] <- data.frame(
        Model = model_names_use[k], Species = spnames[sp],
        Year = yrs, value = as.numeric(tot[sp, ]),
        stringsAsFactors = FALSE)
    }
  }
  plot_df <- do.call(rbind, df_list)
  plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
  plot_df$Species <- factor(plot_df$Species, levels = spnames[species])

  if (mse) {
    agg <- stats::aggregate(value ~ Species + Year, plot_df,
      FUN = function(v) c(m = mean(v),
                          l95 = stats::quantile(v, 0.025, names = FALSE),
                          u95 = stats::quantile(v, 0.975, names = FALSE)))
    mdf <- data.frame(Species = agg$Species, Year = agg$Year, agg$value)
    p <- ggplot2::ggplot(mdf, ggplot2::aes(x = .data$Year)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$l95, ymax = .data$u95),
                           alpha = 0.3, fill = "grey40") +
      ggplot2::geom_line(ggplot2::aes(y = .data$m), linewidth = 1) +
      .facet_species(mdf, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "Biomass eaten as prey") +
      .rceattle_theme()
    return(.save_ggplot(p, file = file, suffix = "biomass_eaten",
                        width = width, height = height))
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$Year, y = .data$value, colour = .data$Model)) +
    ggplot2::geom_line(linewidth = 1) +
    .facet_species(plot_df, scales = "free_y") +
    ggplot2::labs(x = "Year", y = "Biomass eaten as prey")
  if (incl_proj) {
    p <- p + ggplot2::geom_vline(
      xintercept = models[[length(models)]]$data_list$endyr,
      linetype = 2, colour = "grey50")
  }
  p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
  if (nlevels(plot_df$Model) < 2L) p <- p + ggplot2::guides(colour = "none")
  .save_ggplot(p, file = file, suffix = "biomass_eaten",
               width = width, height = height)
}



#' Plot biomass consumed of each prey species by predator
#'
#' @description Function that plots the biomass consumed trends as estimated from Rceattle. Returns and saves a figure with the biomass eaten trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param top_adj Adjustment for top margin
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param minyr first year to plot
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#' @param incl_mean TRUE/FALSE include horizontal long term mean
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#' @param mod_cex Cex of text for model name legend
#'
#' @export
#'
plot_b_eaten_prop <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           right_adj = 0,
           top_adj = 0.15,
           minyr = NULL,
           mohns = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE,
           mod_cex = 1) {

    models <- .as_model_list(Rceattle)
    if (is.null(spnames)) spnames <- models[[1]]$data_list$spnames
    model_names_use <- .model_labels(models, model_names)
    nspp    <- models[[1]]$data_list$nspp
    max_sex <- max(models[[1]]$data_list$nsex)
    if (is.null(species)) species <- 1:nspp

    # Biomass of prey ksp eaten by predator rsp (summed over sex/age), in
    # million mt. B_eaten indexing kept verbatim from the original.
    df_list <- list()
    for (k in seq_along(models)) {
      dl  <- models[[k]]$data_list
      yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
      Be  <- models[[k]]$quantities$B_eaten
      for (rsp in 1:nspp) {
        for (ksp in species) {
          val <- vapply(seq_along(yrs), function(yr)
            sum(Be[c(rsp, (rsp + nspp) * (max_sex - 1)),
                   c(ksp, (ksp + nspp) * (max_sex - 1)), , , yr, drop = FALSE]),
            numeric(1))
          df_list[[length(df_list) + 1L]] <- data.frame(
            Model = model_names_use[k], Prey = spnames[ksp],
            Predator = spnames[rsp], Year = yrs, value = val / 1e6,
            stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Predator <- factor(plot_df$Predator, levels = spnames)
    plot_df$Prey     <- factor(plot_df$Prey, levels = spnames[species])

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$value,
                   colour = .data$Predator, linetype = .data$Model)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::facet_wrap(~ Prey, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "Biomass eaten (million mt)")
    if (incl_proj) {
      p <- p + ggplot2::geom_vline(
        xintercept = models[[length(models)]]$data_list$endyr,
        linetype = 2, colour = "grey50")
    }
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
    if (length(unique(plot_df$Model)) < 2L) p <- p + ggplot2::guides(linetype = "none")
    .save_ggplot(p, file = file, suffix = "biomass_eaten_by_predator",
                 width = width, height = height)
  }




#' Plot natural mortality by age
#'
#' @description Function that plots the natural mortality at age (M1 + M2) as estimated from Rceattle. Returns and saves a figure with the M-at-age trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param age Age to plot M at age
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param lty Line type
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param minyr first year to plot
#' @param incl_proj TRUE/FALSE include projections years
#' @param incl_mean TRUE/FALSE include time series mean as horizontal line
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @export
#'
plot_m_at_age <-
  function(Rceattle,
           file = NULL,
           age = 1,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           lty = 1,
           right_adj = 0,
           minyr = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE) {

    Rceattle <- .as_model_list(Rceattle)
    model_names_use <- .model_labels(Rceattle, model_names)
    if (is.null(spnames)) spnames <- Rceattle[[1]]$data_list$spnames
    nspp <- Rceattle[[1]]$data_list$nspp
    nsex <- Rceattle[[1]]$data_list$nsex
    if (is.null(species)) species <- seq_len(nspp)

    # Total natural mortality (M1 + M2) at a single age, as a time series.
    df_list <- list()
    for (k in seq_along(Rceattle)) {
      dl  <- Rceattle[[k]]$data_list
      yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
      Ma  <- Rceattle[[k]]$quantities$M_at_age
      for (sp in species) {
        for (sex in seq_len(nsex[sp])) {
          sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
          df_list[[length(df_list) + 1L]] <- data.frame(
            Model   = model_names_use[k],
            Species = spnames[sp],
            Sex     = sex_lab,
            Year    = yrs,
            M       = as.numeric(Ma[sp, sex, age, seq_along(yrs)]),
            stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
    plot_df$Species <- factor(plot_df$Species, levels = spnames[species])

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$M,
                   colour = .data$Model, linetype = .data$Sex)) +
      ggplot2::geom_line(linewidth = 1) +
      .facet_species(plot_df, scales = "free_y") +
      ggplot2::labs(x = "Year", y = paste0("M at age ", age))
    if (incl_proj) {
      p <- p + ggplot2::geom_vline(
        xintercept = Rceattle[[length(Rceattle)]]$data_list$endyr,
        linetype = 2, colour = "grey50")
    }
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
    if (length(unique(plot_df$Sex)) < 2L) p <- p + ggplot2::guides(linetype = "none")
    if (nlevels(plot_df$Model) < 2L)      p <- p + ggplot2::guides(colour = "none")

    .save_ggplot(p, file = file, suffix = paste0("m_at_age", age),
                 width = width, height = height)
  }


#' Plot predation mortality by age and predator
#'
#' @description Function that plots the predation mortality at age (M2) by predator as estimated from Rceattle. Returns and saves a figure with the M-at-age trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param age Age to plot M at age
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param top_adj Adjustment for top margin
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param minyr first year to plot
#' @param incl_proj TRUE/FALSE include projections years
#' @param incl_mean TRUE/FALSE include time series mean as horizontal line
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @export
#'
plot_m2_at_age_prop <-
  function(Rceattle,
           file = NULL,
           age = 1,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           right_adj = 0,
           top_adj = 0.15,
           minyr = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE) {

    Rceattle <- .as_model_list(Rceattle)
    if (is.null(spnames)) spnames <- Rceattle[[1]]$data_list$spnames
    model_names_use <- .model_labels(Rceattle, model_names)
    years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    if (incl_proj) years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
    nyrs_vec <- sapply(years, length)
    nyrs <- max(nyrs_vec)
    nsex <- Rceattle[[1]]$data_list$nsex
    nspp <- Rceattle[[1]]$data_list$nspp
    if (is.null(species)) species <- 1:nspp

    # Proportion of M2 (predation mortality) on each prey age attributable to
    # each predator species (kept verbatim from the original extraction).
    m2_at_age_prop <- array(NA, dim = c(nspp, nspp, 2, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      for (ksp in 1:nspp) {
        for (k_sex in 1:nsex[ksp]) {
          for (rsp in 1:nspp) {
            for (yr in 1:nyrs_vec[i]) {
              m2_at_age_prop[rsp, ksp, k_sex, yr, i] <- sum(Rceattle[[i]]$quantities$M2_prop[c(rsp, (rsp + nspp)*(max(nsex)-1)), ksp + (nspp*(k_sex-1)),,age,yr])
            }
          }
        }
      }
    }

    # Tidy: one row per (model, prey species/sex panel, predator, year)
    df_list <- list()
    for (k in seq_along(Rceattle)) {
      yrs <- years[[k]]
      for (ksp in species) {            # prey
        for (sex in seq_len(nsex[ksp])) {
          panel <- if (nsex[ksp] == 1) spnames[ksp]
                   else paste(spnames[ksp], c("Female", "Male")[sex])
          for (rsp in seq_len(nspp)) {  # predator
            df_list[[length(df_list) + 1L]] <- data.frame(
              Model     = model_names_use[k],
              Prey      = panel,
              Predator  = spnames[rsp],
              Year      = yrs,
              Proportion = as.numeric(m2_at_age_prop[rsp, ksp, sex, seq_along(yrs), k]),
              stringsAsFactors = FALSE)
          }
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Predator <- factor(plot_df$Predator, levels = spnames)

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$Proportion,
                   colour = .data$Predator, linetype = .data$Model)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::facet_wrap(~ Prey, scales = "free_y") +
      ggplot2::labs(x = "Year", y = paste0("M2 proportion at age ", age))
    if (incl_proj) {
      p <- p + ggplot2::geom_vline(
        xintercept = Rceattle[[length(Rceattle)]]$data_list$endyr,
        linetype = 2, colour = "grey50")
    }
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
    if (length(unique(plot_df$Model)) < 2L) p <- p + ggplot2::guides(linetype = "none")

    .save_ggplot(p, file = file, suffix = paste0("m2_at_age_prop", age),
                 width = width, height = height)
  }


#' plot F
#'
#' @description Function that plots the F time series per species from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param spnames Species names for legend
#' @param add_ci NOT WORKING If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param minyr First year to plot
#' @param height plot height
#' @param width plot width
#' @param incl_proj TRUE/FALSE, include projection years
#' @param mod_cex Cex of text for model name legend
#' @param mse Is an MSE object from \code{\link{load_mse}} or \code{\link{run_mse}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#' @param maxyr max year to plot
#' @param alpha shading for confidence intervals
#' @param mod_avg is the list a model average? (DEPRECATED)
#' @importFrom stats quantile
#' @importFrom grDevices png dev.off adjustcolor
#' @importFrom graphics layout par plot.new abline legend polygon lines plot
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_f <- function(Rceattle,
                   file = NULL,
                   model_names = NULL,
                   line_col = NULL,
                   species = NULL,
                   spnames = NULL,
                   add_ci = FALSE,
                   lwd = 3,
                   right_adj = 0,
                   width = 7,
                   height = 6.5,
                   minyr = NULL,
                   maxyr = NULL,
                   incl_proj = FALSE,
                   mod_cex = 1,
                   alpha = 0.4,
                   mod_avg = rep(FALSE, length(Rceattle)),
                   mse = FALSE,
                   OM = TRUE) {

  # Fishing mortality shares the time-series machinery (output = "F_spp");
  # only the Ftarget / Flimit reference lines are F-specific. Build the plot
  # without saving, add the reference lines, then save so the F lines are
  # included in the figure file.
  p <- plot_timeseries(Rceattle, output = "F_spp",
                       ylab = "Fishing mortality (F)",
                       file = NULL, model_names = model_names,
                       line_col = line_col, species = species,
                       spnames = spnames, add_ci = add_ci, lwd = lwd,
                       right_adj = right_adj, width = width, height = height,
                       minyr = minyr, maxyr = maxyr, incl_proj = incl_proj,
                       mod_cex = mod_cex, alpha = alpha, mod_avg = mod_avg,
                       mse = mse, OM = OM)

  # F reference points: target (blue) and limit (red), per species facet.
  models <- .as_model_list(Rceattle, mse = mse, OM = OM)
  sp_all  <- if (is.null(spnames)) models[[1]]$data_list$spnames else spnames
  if (is.null(species)) species <- seq_along(sp_all)
  ftarget <- models[[1]]$quantities$Ftarget
  flimit  <- models[[1]]$quantities$Flimit
  if (!is.null(ftarget) && !is.null(flimit)) {
    ref_df <- data.frame(
      Species = factor(sp_all[species], levels = levels(p$data$Species)),
      Ftarget = ftarget[species],
      Flimit  = flimit[species])
    p <- p +
      ggplot2::geom_hline(data = ref_df, inherit.aes = FALSE,
                          ggplot2::aes(yintercept = .data$Ftarget),
                          colour = "blue", linetype = 2) +
      ggplot2::geom_hline(data = ref_df, inherit.aes = FALSE,
                          ggplot2::aes(yintercept = .data$Flimit),
                          colour = "red", linetype = 2)
  }

  .save_ggplot(p, file = file, suffix = "f_trajectory",
               width = width, height = height)
}




#' Plot ration
#'
#' @description Function that plots the ration across ages (minage:nages) as estimated from Rceattle. Returns and saves a figure with the ration trajectory. Ration is multiplied by biomass-at-age/sex to get population level estimates
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param minage minage to plot ration (i.e. age "minage"+)
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param lty Line type
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param minyr first year to plot
#' @param incl_proj TRUE/FALSE include projections years
#' @param incl_mean TRUE/FALSE include time series mean as horizontal line
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @export
#'
plot_ration <-
  function(Rceattle,
           file = NULL,
           minage = 1,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           lty = 1,
           right_adj = 0,
           minyr = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE) {

    models <- .as_model_list(Rceattle)
    if (is.null(spnames)) spnames <- models[[1]]$data_list$spnames
    model_names_use <- .model_labels(models, model_names)
    nspp  <- models[[1]]$data_list$nspp
    nsex  <- models[[1]]$data_list$nsex
    nages <- models[[1]]$data_list$nages
    if (is.null(species)) species <- 1:nspp

    # Population-level consumption: ration (consumption_at_age) x
    # biomass-at-age, summed over age minage+, in million mt.
    df_list <- list()
    for (k in seq_along(models)) {
      dl  <- models[[k]]$data_list
      yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
      Caa <- models[[k]]$quantities$consumption_at_age
      Baa <- models[[k]]$quantities$biomass_at_age
      for (sp in species) {
        ages <- minage:nages[sp]
        for (sex in seq_len(nsex[sp])) {
          sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
          val <- vapply(seq_along(yrs), function(yr)
            sum(Caa[sp, sex, ages, yr] * Baa[sp, sex, ages, yr]),
            numeric(1)) / 1e6
          df_list[[length(df_list) + 1L]] <- data.frame(
            Model = model_names_use[k], Species = spnames[sp], Sex = sex_lab,
            Year = yrs, value = val, stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
    plot_df$Species <- factor(plot_df$Species, levels = spnames[species])

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$value,
                   colour = .data$Model, linetype = .data$Sex)) +
      ggplot2::geom_line(linewidth = 1) +
      .facet_species(plot_df, scales = "free_y") +
      ggplot2::labs(x = "Year",
                    y = paste0("Consumption (million mt), age ", minage, "+"))
    if (incl_proj) {
      p <- p + ggplot2::geom_vline(
        xintercept = models[[length(models)]]$data_list$endyr,
        linetype = 2, colour = "grey50")
    }
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
    if (length(unique(plot_df$Sex)) < 2L) p <- p + ggplot2::guides(linetype = "none")
    if (nlevels(plot_df$Model) < 2L)      p <- p + ggplot2::guides(colour = "none")
    .save_ggplot(p, file = file, suffix = paste0("ration", minage, "plus"),
                 width = width, height = height)
  }

