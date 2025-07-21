#' Plot stock recruit function
#'
#' @description Function the plots the stock recruit function as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param spnames Species names for legend
#' @param incl_proj TRUE/FALSE, include projection years for environmental relationship
#'
#' @param lwd Line width as specified by user
#'
#' @export
plot_stock_recruit <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           width = 7,
           height = 6.5,
           species = NULL,
           spnames = NULL,
           lwd = 3,
           lty = 1,
           incl_proj = FALSE,
           plot_env = FALSE,
           mod_cex = 1) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Setup ----
    # * Extract dimensions ----
    years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
    nyrs <- max(sapply(years, length))

    hindyears <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    nyrshind <- sapply(hindyears, length)

    srr_fun <- sapply(Rceattle, function(x) x$data_list$srr_fun)
    srr_pred_fun <- sapply(Rceattle, function(x) x$data_list$srr_pred_fun)
    nspp <- Rceattle[[1]]$data_list$nspp
    minage <- Rceattle[[1]]$data_list$minage
    est_dynamics <- Rceattle[[1]]$data_list$estDynamics

    # Species names
    if(is.null(spnames)){
      spnames =  Rceattle[[1]]$data_list$spnames
    }

    # Species included
    if(is.null(species)){
      species <- 1:nspp
    }

    # * Extract quantities ----
    ssb_array <-
      array(NA, dim = c(nspp, max(nyrshind), length(Rceattle)))
    rec_array <-
      array(NA, dim = c(nspp, max(nyrshind), length(Rceattle)))
    rec_pars <-
      array(NA, dim = c(nspp, 3, length(Rceattle)))
    beta_rec_pars <- list()
    env_data <- list()


    for (mod in 1:length(Rceattle)) {
      ssb_array[,1:nyrshind[mod],mod] <- Rceattle[[mod]]$quantities$ssb[,1:nyrshind[mod]]/1000000
      rec_array[,1:nyrshind[mod],mod] <- Rceattle[[mod]]$quantities$R[,1:nyrshind[mod]]/1000000
      rec_pars[,,mod] <- Rceattle[[mod]]$estimated_params$rec_pars
      beta_rec_pars[[mod]] <- Rceattle[[mod]]$estimated_params$beta_rec_pars
      env_data[[mod]] <-  Rceattle::rearrange_dat(Rceattle[[mod]]$data_list)$env_index
      if(!incl_proj){
        env_data[[mod]] <- env_data[[mod]][1:nyrshind[mod],]
      }
    }

    # * Plot limits ----
    ymax <- c()
    ymin <- c()

    xmax <- c()
    xmin <- c()
    for (sp in 1:nspp) {
      xmax[sp] <- max(c(ssb_array[sp,,], 0), na.rm = T)
      xmin[sp] <- min(c(ssb_array[sp,,], 0), na.rm = T)


      ymax[sp] <- max(c(rec_array[sp,,], 0), na.rm = T)
      ymin[sp] <- min(c(rec_array[sp,,], 0), na.rm = T)
    }

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }

    if(length(lty) == 1){
      lty <- rep(lty, length(Rceattle))
    }

    # Plot stock recruit ----
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "stock_recruit_function", ".png")
        png(
          file = filename ,
          width = width,# 169 / 25.4,
          height = height,# 150 / 25.4,
          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(length(species) + 2), nrow = (length(species) + 2)), heights = c(0.2, rep(1, length(species)), 0.2))
      par(
        mar = c(1, 3 , 0.5 , 1) ,
        oma = c(0.5, 0 , 0.5 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      # 1) Species loop ----
      for (sp in 1:length(species)) {

        # * Base plot ----
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[species[sp]], ymax[species[sp]]),
          xlim = c(xmin[species[sp]], xmax[species[sp]]),
          xlab = NA,
          ylab = NA
        )

        # 2) Model loop ----
        for (mod in 1:length(Rceattle)) {

          # * No SRR ----
          if(srr_pred_fun[mod] == 0){
            abline(h = exp(rec_pars[species[sp], 1, mod]),
                   lty = 1, lwd = lwd, col = line_col[mod])
          }

          # * Beverton and holt ----
          if(srr_pred_fun[mod] %in% c(2,3)){
            curve(exp(rec_pars[species[sp], 2, mod]) * x / (1+exp(rec_pars[species[sp], 3, mod]) * x * 1000000),
                  from = 0, to = xmax[species[sp]],
                  lty = 1, lwd = lwd, col = line_col[mod],
                  add = TRUE)
          }

          # * Ricker ----
          if(srr_pred_fun[mod] %in% c(4,5)){
            curve(exp(rec_pars[species[sp], 2, mod]) * x * exp(-exp(rec_pars[species[sp], 3, mod]) * x),
                  from = 0, to = xmax[species[sp]],
                  lty = lty[mod], lwd = lwd, col = line_col[mod],
                  add = TRUE)
          }

          # * Ricker w/ environment ----
          if(srr_pred_fun[mod] == 5 & plot_env){

            # - Env curves
            for(env in 1:nrow(env_data[[mod]])){
              curve(exp(rec_pars[species[sp], 2, mod]) * exp(env_data[[mod]][env,] %*% beta_rec_pars[[mod]][species[sp],]) * x * exp(-exp(rec_pars[species[sp], 3, mod]) * x),
                    from = 0, to = xmax[species[sp]],
                    lty = lty[mod], lwd = lwd/2, col = t_col(line_col[mod], percent = 85),
                    add = TRUE)
            }

            # - Base curve
            curve(exp(rec_pars[species[sp], 2, mod]) * x * exp(-exp(rec_pars[species[sp], 3, mod]) * x),
                  from = 0, to = xmax[species[sp]],
                  lty = lty[mod], lwd = lwd, col = line_col[mod],
                  add = TRUE)

            # Updated curves
          }

          # * Recruitment points ----
          points(
            x = ssb_array[species[sp], 1:nyrshind[mod], mod],
            y = rec_array[species[sp], 1:nyrshind[mod], mod],
            pch = 21,
            lwd = 0.5,
            bg = t_col(line_col[mod], percent = 50)
          ) # Median
        }

        # Labels and legends ----
        # - Axis labels
        if(sp == 2){
          mtext("Recruitment (millions)", side = 2, line = 1.7, cex = 0.9)
        }
        if(sp == length(species)){
          mtext("Spawning stock biomass (million mt)", side = 1, line = 2, cex = 0.9)
        }

        # - Legends
        legend("topleft",
               legend = spnames[species[sp]],
               bty = "n",
               cex = 1)

        if (species[sp] == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              lty = lty,
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = mod_cex
            )
          }
        }
      }

      if (i == 2) {
        dev.off()
      }
    }
  }



#' #https://www.dataanalytics.org.uk/make-transparent-colors-in-r/
#'
#' @param color color name
#' @param percent % transparency
#' @param name an optional name for the color
#'
#'
t_col <- function(color, percent = 50, name = NULL) {

  ## Get RGB values for named color
  rgb.val <- col2rgb(color)

  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)

  ## Save the color
  invisible(t.col)
}
