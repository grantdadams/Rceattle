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
           mod_cex = 1) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects ----
    years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
    hindyears <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    nyrs <- max(sapply(years, length))
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
    spp <- species

    # Extract derived quantities ----
    ssb_array <-
      array(NA, dim = c(nspp, max(nyrshind), length(Rceattle)))
    rec_array <-
      array(NA, dim = c(nspp, max(nyrshind), length(Rceattle)))
    rec_pars <-
      array(NA, dim = c(nspp, 3, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      ssb_array[,1:nyrshind[i],i] <- Rceattle[[i]]$quantities$ssb[,1:nyrshind[i]]/1000000
      rec_array[,1:nyrshind[i],i] <- Rceattle[[i]]$quantities$R[,1:nyrshind[i]]/1000000
      rec_pars[,,i] <- Rceattle[[i]]$estimated_params$rec_pars
    }

    # Plot limits
    ymax <- c()
    ymin <- c()

    xmax <- c()
    xmin <- c()
    for (i in 1:nspp) {
      xmax[i] <- max(c(ssb_array[i,,], 0), na.rm = T)
      xmin[i] <- min(c(ssb_array[i,,], 0), na.rm = T)


      ymax[i] <- max(c(rec_array[i,,], 0), na.rm = T)
      ymin[i] <- min(c(rec_array[i,,], 0), na.rm = T)
    }

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot stock recruit function ----
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
      layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.2, rep(1, length(spp)), 0.2))
      par(
        mar = c(1, 3 , 0.5 , 1) ,
        oma = c(0.5, 0 , 0.5 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      for (j in 1:length(spp)) {
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[spp[j]], ymax[spp[j]]),
          xlim = c(xmin[spp[j]], xmax[spp[j]]),
          xlab = NA,
          ylab = NA
        )

        # Curves
        for (k in 1:length(Rceattle)) {

          # Raw points
          points(
            x = ssb_array[spp[j], 1:nyrshind[k], k],
            y = rec_array[spp[j], 1:nyrshind[k], k],
            pch = 16,
            col = t_col(line_col[k], percent = 50)
          ) # Median

          # - No SRR
          if(srr_pred_fun[k] == 0){
            abline(h = exp(rec_pars[spp[j], 1, k]),
                  lty = 1, lwd = lwd, col = line_col[k])
          }

          # - Beverton and holt
          if(srr_pred_fun[k] %in% c(2,3)){
            curve(exp(rec_pars[spp[j], 2, k]) * x / (1+exp(rec_pars[spp[j], 3, k]) * x * 1000000),
                  from = 0, to = xmax[spp[j]],
                  lty = 1, lwd = lwd, col = line_col[k],
                  add = TRUE)
          }

          # - Ricker
          if(srr_pred_fun[k] %in% c(4,5)){
            curve(exp(rec_pars[spp[j], 2, k]) * x * exp(-exp(rec_pars[spp[j], 3, k]) * x),
                  from = 0, to = xmax[spp[j]],
            lty = 1, lwd = lwd, col = line_col[k],
            add = TRUE)
          }
        }

        # Axis labels
        if(j == 2){
          mtext("Recruitment (millions)", side = 2, line = 1.7, cex = 0.9)
        }
        if(j == length(spp)){
          mtext("Spawning stock biomass (million mt)", side = 1, line = 2, cex = 0.9)
        }

        # Legends
        legend("topleft",
               legend = spnames[spp[j]],
               bty = "n",
               cex = 1)

        if (spp[j] == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              lty = rep(1, length(line_col)),
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
