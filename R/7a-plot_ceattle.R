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
#' @export
rich.colors.short <- function(n,alpha=1){
  x <- seq(0, 1, length = n)
  r <- 1/(1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
  rgb.m <- matrix(c(r, g, b), ncol = 3)
  rich.vector <- apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha=alpha))
}

#' plot_biomass
#'
#' @description Function the plots the mean biomass and 95% CI trends as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param add_ci If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param minyr First year to plot
#' @param height
#' @param width
#' @param save Save biomass?
#' @param incl_proj TRUE/FALSE, include projection years
#' @param mod_cex Cex of text for model name legend
#' @param mse Is if an MSE object from \code{\link{load_mse}} or \code{\link{mse_run}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#'
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_biomass <- function(Rceattle,
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
                         OM = TRUE,
                         reference = NULL) {

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    nmse = length(Rceattle)
    add_ci = TRUE
    incl_proj = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }

  # Add reference model
  if(!is.null(reference)){
    Rceattle <- c(Rceattle, list(reference))
  }

  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  if(is.null(minyr)){minyr <- min((sapply(Years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp

  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
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
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$biomass[,1:nyrs_vec[i]]

    # Get SD of quantity
    if(add_ci & !mse){
      sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "biomass")
      sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
      quantity_sd[,  1:nyrs_vec[i], i] <-
        replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_temp[1:(nyrs_vec[i] * nspp)])
    }

    # - Model average
    if(mod_avg[i]){
      log_quantity_sd[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples$biomass[,1:nyrs_vec[i],], c(1,2), function(x) sd(as.vector(log(x))))
      log_quantity_mu[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples$biomass[,1:nyrs_vec[i],], c(1,2), function(x) mean(as.vector(log(x))))
    }
  }

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){
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

  # - Rescale
  quantity <- quantity / 1000000
  quantity_upper95 <- quantity_upper95 / 1000000
  quantity_lower95 <- quantity_lower95 / 1000000
  quantity_upper50 <- quantity_upper50 / 1000000
  quantity_lower50 <- quantity_lower50 / 1000000


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_biomass_species_", i, ".csv"))
    }
  }


  ## Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * 1.2

  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }
  if(!is.null(reference)){
    line_col <- c(line_col, 1)
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_biomass_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = NA,
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      if(j == 2){
        mtext("Age-1+ biomass (million mt)", side = 2, line = 1.7, cex = 0.9)
      }

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
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


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(Years[[k]], rev(Years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(Years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = Years[[k]],
          y = quantity[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd,
          col = line_col[k]
        ) # Median
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}




#' plot_recruitment
#'
#' @description Function the plots the mean recruitment and 95% CI trends as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param spnames Species names for legend
#' @param add_ci If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param minyr First year to plot
#' @param height
#' @param width
#' @param save Save biomass?
#' @param incl_proj TRUE/FALSE, include projection years
#' @param mod_cex Cex of text for model name legend
#' @param mse Is if an MSE object from \code{\link{load_mse}} or \code{\link{mse_run}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#'
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_recruitment <- function(Rceattle,
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

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    add_ci = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }


  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  if(is.null(minyr)){minyr <- min((sapply(Years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp

  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
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
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$R[,1:nyrs_vec[i]]

    # Get SD of quantity
    if(add_ci & !mse){
      sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "R")
      sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
      quantity_sd[,  1:nyrs_vec[i], i] <-
        replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_temp[1:(nyrs_vec[i] * nspp)])
    }

    # - Model average
    if(mod_avg[i]){
      log_quantity_sd[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples$R[,1:nyrs_vec[i],], c(1,2), function(x) sd(as.vector(log(x))))
      log_quantity_mu[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples$R[,1:nyrs_vec[i],], c(1,2), function(x) mean(as.vector(log(x))))
    }
  }

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){
    ptarget <- ptarget[1,]
    plimit <- plimit[1,]

    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.975) )
    quantity_lower95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.025) )
    quantity_upper50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.75) )
    quantity_lower50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.25) )
    quantity <- apply( quantity, c(1,2), mean ) # Get mean quantity

    # -- Put back in array for indexing below
    quantity <- array(quantity, dim = c(nspp, nyrs,  1))
    quantity_upper95 <- array(quantity_upper95, dim = c(nspp, nyrs,  1))
    quantity_lower95 <- array(quantity_lower95, dim = c(nspp, nyrs,  1))
    quantity_upper50 <- array(quantity_upper50, dim = c(nspp, nyrs,  1))
    quantity_lower50<- array(quantity_lower50, dim = c(nspp, nyrs,  1))
  }

  # - Model Average
  for (i in 1:length(Rceattle)) {
    if(mod_avg[i]){
      quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
    }
  }

  # - Rescale
  quantity <- quantity / 1000000
  quantity_upper95 <- quantity_upper95 / 1000000
  quantity_lower95 <- quantity_lower95 / 1000000
  quantity_upper50 <- quantity_upper50 / 1000000
  quantity_lower50 <- quantity_lower50 / 1000000


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_recruitment_species_", i, ".csv"))
    }
  }


  ## Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * 1.2

  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_recruitment_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = NA,
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      if(j == 2){
        mtext("Age-1 recruits (million)", side = 2, line = 1.7, cex = 0.9)
      }

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
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


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(Years[[k]], rev(Years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(Years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = Years[[k]],
          y = quantity[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd,
          col = line_col[k]
        ) # Median
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}




#' Plot selectivity
#'
#' @description Function the plots the fishery and survey selectivity as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
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

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp
    fleet_control <- Rceattle[[1]]$data_list$fleet_control
    nflt <- nrow(Rceattle[[1]]$data_list$fleet_control)
    nages <- Rceattle[[1]]$data_list$nages
    minage <- Rceattle[[1]]$data_list$minage
    nsex <- Rceattle[[1]]$data_list$nsex
    est_dynamics <- Rceattle[[1]]$data_list$estDynamics

    # Get biomass
    selectivity_array <-
      array(NA, dim = c(nflt, 2, max(nages), nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      selectivity_array[, , , ,i] <- Rceattle[[i]]$quantities$sel[,,,]
    }

    # Plot limits
    ymax_sel <- c()
    ymin_sel <- c()
    for (i in 1:dim(selectivity_array)[1]) {
      ymax_sel[i] <- max(c(selectivity_array[i,,,,], 0), na.rm = T)
      ymin_sel[i] <- min(c(selectivity_array[i,,,,], 0), na.rm = T)
    }

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }

    max_obj <- nflt


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)

    #################################
    # Selectivity time series
    #################################
    for(j in 1:nflt){
      for (i in 1:loops) {
        # Species
        sp <- fleet_control$Species[which(fleet_control$Fleet_code == j)]

        for(sex in 1:nsex[sp]){

          # Get sex for legend
          legend_sex = sex
          legend_sex2 = ifelse(sex == 1, "Female", "Male")
          if(nsex[sp] == 1){
            legend_sex <- 0
            legend_sex2 = "Combined"
          }

          if (i == 2) {
            filename <- paste0(file, "time-varying_selectivity_fleet",j,"_sex",legend_sex, ".png")
            png(
              file = filename ,
              width = width,
              height = height,
              units = "in",
              res = 300
            )
          }
          sel_subset <- (selectivity_array[j, sex, 1:nages[sp], 1:nyrs, 1])

          par(
            mar = c(3.2, 3.2 , 1 , 0.5) ,
            oma = c(0 , 0 , 0 , 0),
            tcl = -0.8,
            mgp = c(10, 0.6, 0)
          )
          persp(y = Years[[1]], x =  (1:nages[sp]) - 1 + minage[sp], z = sel_subset, col="white",xlab = "Age",ylab= "\n\nYear", zlab= "\n\nSelectivity",expand=0.5,box=TRUE,ticktype="detailed",phi=35,theta=-19, main = NA, zlim = c(ymin_sel[i], ymax_sel[j]))
          mtext(text = paste(legend_sex2, as.character(fleet_control$Fleet_name[j])), side = 3, cex = 0.8, line = -.5)

          if (i == 2) {
            dev.off()
          }
        }
      }
    }

    #################################
    # Terminal selectivity
    #################################
    for(sp in 1:nspp){
      if(est_dynamics[sp] == 0){
        for(sex in 1:nsex[sp]){

          # Get sex for legend
          legend_sex = sex
          legend_sex2 = ifelse(sex == 1, "Female", "Male")
          if(nsex[sp] == 1){
            legend_sex <- 0
            legend_sex2 = "Combined"
          }

          for (i in 1:loops) {
            if (i == 2) {

              filename <- paste0(file, "_terminal_selectivity_species",sp,"_sex",legend_sex, ".png")
              png(
                file = filename ,
                width = width,
                height = height,
                units = "in",
                res = 300
              )
            }

            fleets <- fleet_control$Fleet_code[which(fleet_control$Species == sp)]
            flt_colors <- rich.colors.short(length(fleets))

            # Plot configuration
            par(
              mar = c(3.2, 3.2 , 0.5 , 0.5) ,
              oma = c(0 , 0 , 0 , 0),
              tcl = -0.35,
              mgp = c(1.75, 0.5, 0)
            )

            plot(
              y = NA,
              x = NA,
              ylim = c(min(ymin_sel[fleets]), max(ymax_sel[fleets])),
              xlim = c(min(0),  max(nages[sp], na.rm = TRUE)),
              xlab = "Age",
              ylab = "Terminal selectivity"
            )

            # Mean selectivity
            for (flt in 1:length(fleets)) {
              lines(
                x = 1:nages[sp],
                y = selectivity_array[fleets[flt], sex, 1:nages[sp], nyrs, 1],
                lty = 1,
                lwd = lwd,
                col = flt_colors[flt]
              )
            }

            # Legends
            legend("bottomright", paste(legend_sex2, as.character(fleet_control$Fleet_name[fleets])), col = flt_colors, bty = "n", lty = rep(1, length(fleets)), lwd = rep(2, length(fleets)), cex = 0.8)

            # Save plot
            if (i == 2) {
              dev.off()
            }
          }
        }
      }
    }
  }




#' Plot functional form
#'
#' @description Function to plot the functional form estimated or specified by \code{\link{Rceattle}}
#'
#' @param params Parameter list object from \code{\link{build_params}} or \code{\link{Rceattle}}
#' @param pred Predator index
#' @param pred_age Predator age
#' @param prey Prey index
#' @export
plot_form <- function( params = NULL, pred = 1, pred_age = 1, prey = 1, msmMode = 3){

  # Get indices
  rsp = pred
  ksp = prey

  # Get parameter values
  H_1 <- exp(params$logH_1)
  H_1a <- exp(params$logH_1a)
  H_1b <- exp(params$logH_1b)

  H_2 <- exp(params$logH_2)
  H_3 <- exp(params$logH_3)

  H_4 <- params$H_4

  # Set up ratios
  Pred_r <- seq(from = 0.001, to = 5, length.out = 100) # Pred biomass relative to equilibrium
  Prey_r <- seq(from = 0.001, to = 5, length.out = 100) # Prey biomass relative to equilibrium

  # Calculate functional form
  Term = H_1[rsp, ksp] * (1 + H_1a[rsp] * H_1b[rsp] / (pred_age + H_1b[rsp]))
  response <- matrix(NA, ncol = length(Prey_r), nrow = length(Pred_r))
  rownames(response) <- Pred_r
  colnames(response) <- Prey_r

  for(i in 1:length(Pred_r)){
    for(j in 1:length(Prey_r)){
      response[i, j] <- Prey_r[j] * switch (
        as.character(msmMode),
        "2" = { # Holling Type I (linear)
          Term},
        "3" = { #  Holling Type II
          Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] )
        },
        "4" = { #  Holling Type III
          Term * (1 + H_2[rsp, ksp]) * ((Prey_r[j] ) ^ H_4[rsp, ksp]) / (1 + H_2[rsp, ksp] * ((Prey_r[j] ) ^ H_4[rsp, ksp])  )},
        "5" = { #  predator interference
          Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] + H_3[rsp, ksp] * (Pred_r[i] - 1) )
        },
        "6" = { # predator preemption
          Term * (1 + H_2[rsp, ksp] ) / ( (1 + H_2[rsp, ksp] * Prey_r[j]) * (1 + H_3[rsp, ksp] * (Pred_r[i] - 1)) )
        },
        "7" = { # Hassell-Varley
          Term * (2 + H_2[rsp, ksp] ) / (1 + H_2[rsp, ksp] * Prey_r[j] + ((Prey_r[j] ) ^ H_4[rsp, ksp]))
        },
        "8" = { #  Ecosim
          Term / (1 + H_3[rsp, ksp] * (Pred_r[i] - 1 ))},
        {
          print("msmMode not implemented")
        }
      )
    }
  }

  response_reshape <- reshape2::melt(response, id.vars = c("Pred_r", "Prey_r"), measure.vars = "Response")
  colnames(response_reshape) <- c("Pred_r", "Prey_r", "Response")

  if(msmMode %in% c(2, 3, 4, 7)){
    plot(x = response_reshape$Prey_r, y = response_reshape$Response, xlab = "Prey ratio", ylab = "Functional response", type = "l")
  }

  if(msmMode %in% c(8)){
    plot(x = response_reshape$Pred_r, y = response_reshape$Response, xlab = "Pred ratio", ylab = "Functional response", type = "l")
  }

  if(msmMode %in% c(5, 6)){
    filled.contour(response, xlab = "Prey ratio", ylab = "Predator ratio", col = oce.colorsDensity(30))
  }
}





#' Plot M1 + M2
#'
#' @description Function the plots the M1 and M2 as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param incl_proj Include the projection years (TRUE/FALSE)
#' @param zlim zlim for M1 + M2 plots. Character - use max range across species in model. NULL - use species specific ranges. Vector of two.
#' @param type 0 = Tiles, 1 = contour, 2 = facet lines, 3 = persp
#' @param width Plot width when saved "inches"
#' @param height Plot height when saved "inches"
#' @param title Additional title to add. Will also add species names if not NULL
#' @param title_cex Font size for title
#' @param spp Species to plot. Plots all if null.
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
           type = 0,
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

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    if(length(Rceattle) > 1){
      stop("Can only plot one model")
    }

    # Extract data objects
    if(is.null(minyr)){ minyr <- Rceattle[[1]]$data_list$styr}

    Years <- minyr:Rceattle[[1]]$data_list$endyr
    if(incl_proj){
      Years <- minyr:Rceattle[[1]]$data_list$projyr
    }
    nyrs_vec <- length(Years)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))

    nspp <- Rceattle[[1]]$data_list$nspp
    spnames <- Rceattle[[1]]$data_list$spnames
    estdynamics <- Rceattle[[1]]$data_list$estDynamics
    nages <- Rceattle[[1]]$data_list$nages

    if(!is.null(maxage)){
      nages <- sapply(nages, function(x) ifelse(x > maxage, maxage, x))
    }


    minage <- Rceattle[[1]]$data_list$minage
    nsex <- Rceattle[[1]]$data_list$nsex

    # Get M
    M_array <-
      array(NA, dim = c(nspp, 2, max(nages), nyrs, length(Rceattle)))
    M1_array <-
      array(NA, dim = c(nspp, 2, max(nages), length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      M1_array[, , ,i] <- Rceattle[[i]]$quantities$M1[,,1:max(nages)]
      if(!M2){
        M_array[, , , ,i] <- Rceattle[[i]]$quantities$M[,,1:max(nages),(1:nyrs)+(minyr - Rceattle[[1]]$data_list$styr)]
      }
      if(M2){
        M_array[, , , ,i] <- Rceattle[[i]]$quantities$M2[,,1:max(nages),(1:nyrs)+(minyr - Rceattle[[1]]$data_list$styr)]
      }
    }

    if(log){
      M1_array = log(M1_array)
      M_array = log(M_array)
    }

    # Plot limits
    zmax <- c()
    zmin <- c()
    for (i in 1:dim(M_array)[1]) {
      zmax[i] <- max(c(M_array[i,,,,], 0), na.rm = T)
      zmin[i] <- min(c(M_array[i,,,,], 0), na.rm = T)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)

    #################################
    # Mortality time series
    #################################
    if(is.null(species)){
      species <- 1:nspp
    }
    spp <- species


    # Species
    for(j in 1:nspp){
      sp <- j

      if(estdynamics[j] == 0 & sp %in% spp){

        # Sexes
        for(sex in 1:nsex[sp]){

          # Get sex for legend
          legend_sex = sex
          legend_sex2 = ifelse(sex == 1, "Female", "Male")
          if(nsex[sp] == 1){
            legend_sex <- 0
            legend_sex2 = "Combined"
          }

          # Save
          for (i in 1:loops) {
            if (i == 2) {
              filename <- paste0(file, "predation_and_residual_mortality_spp_",sp,"_sex_",legend_sex2,".png")
              png(
                file = filename ,
                width = width,
                height = height,
                units = "in",
                res = 300
              )
            }

            # Subset mortality data
            m_subset <- (M_array[j, sex, (1:nages[sp]), 1:nyrs, 1])

            # Get ages
            ages <- (1:(nages[sp])) - 1 + minage[sp]

            # Rearrange data
            data <- data.frame(Year = rep(Years, each = length(ages)), Age = rep(ages, length(Years)), M = c(m_subset))

            # Plot limits
            if(is.null(zlim)){
              zlim <- c(zmin[sp], zmax[sp])
            }

            if(is.character(zlim)){
              zlim <- c(min(zmin), max(zmax))
            }

            # Plot as tiles
            if(type == 0){
              p = ggplot2::ggplot(data, aes(y = Age, x = Year, zmin = zlim[1], zmax = zlim[2])) + geom_tile(aes(fill = M))  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) + coord_equal() +  scale_x_continuous(expand = c(0, 0))+ theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
              if(!is.null(title)){
                p = p + ggtitle(paste0(title,": ",spnames[j] )) + theme(plot.title = element_text(size = title_cex))
              }
              if(log){
                p = p + scale_fill_viridis_c("log(M1 + M2)", limits = c(zlim[1], zlim[2]))
              } else {
                p = p + scale_fill_viridis_c("M1 + M2", limits = c(zlim[1], zlim[2]))
              }
              print(p)
            }


            # Plot as contours
            if(type == 1){
              print(ggplot2::ggplot(data, aes(y = Age, x = Year, z = M, zmin = zlim[1], zmax = zlim[2])) + geom_contour(colour = 1, size = 0.5) + geom_contour_filled()  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) +  scale_x_continuous(expand = c(0, 0)) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_fill_viridis_d("M1 + M2"))
            }

            # Plot as facets
            if(type == 2){
              p = ggplot(data=data, aes(x=Year, y = M, colour = Age, group = Age)) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + geom_line(size = 2) + scale_color_viridis_c("Age")
              print(p)
            }

            # Plot as persp
            if(type == 3){
              par( mar=c(1 , 2 , 1 , 1) , tcl=-.25 , mgp=c(2 ,  1 ,  0) ,  oma=c(0 , 2 , 0 , 0))
              pmat = persp(y = Years, x = ages, z = m_subset, zlab = NA, zlim = zlim, xlab = "Age", ylab = "Year", theta = theta, ticktype = "detailed")
              mtext(ifelse(M2, "M2", "M"), side = 2, line = 0.5, at = 0)
              if(M2){
                text(-0.25,.15, labels = paste0("M1 = ",round((M1_array[j, sex, 1, 1]), 3)))
              }


              if(nsex[sp] == 1){
                mtext(paste0(title,": ",spnames[j]), side = 3, line = -2, at = 0)
              }
              if(nsex[sp] == 2){
                mtext(paste0(title,": ",spnames[j], " ",legend_sex2), side = 3, line = -2, at = 0)
              }
            }

            if (i == 2) {
              dev.off()
            }
          }
        }
      }
    }
  }







#' Plot maturity
#'
#' @description Function the plots the maturity of each species
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
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

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    maturity <- list()
    for(i in 1:length(Rceattle)){
      maturity[[i]] <- Rceattle[[i]]$data_list$pmature[,-1] # Remove species column
    }

    nspp <- Rceattle[[1]]$data_list$nspp
    nages <- Rceattle[[1]]$data_list$nages

    # Line colors
    if (is.null(line_col)) {
      line_col <- oce::oce.colorsViridis(length(Rceattle))
    }


    # Species names
    if(is.null(species)){
      species =  Rceattle[[1]]$data_list$spnames
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_maturity", ".png")
        png(
          file = filename ,
          width = width,
          height = height,
          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:((nspp + 2)), nrow = (nspp + 2), ncol = 1, byrow = FALSE), heights = c(0.2, rep(1, nspp), 0.3))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )

      plot.new()

      # Survey selectivity
      for (j in 1:nspp) {
        plot(
          y = NA,
          x = NA,
          ylim = c(0, 1.1),
          xlim = c(min(0), max(nages, na.rm = TRUE)),
          xlab = "Age",
          ylab = "Maturity",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        if(j == nspp){
          mtext(side = 1, "Age", cex  = 0.75, line = 2)
        }

        # Mean maturity
        for (k in 1:length(maturity)) {
          lines(
            x = 1:nages[j],
            y = maturity[[k]][j, 1:nages[j]],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          )
        }

        # Species legends
        legend("topleft", species[j], bty = "n", cex = 1)

        # Model name legends
        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "bottomright",
              legend = model_names,
              lty = rep(1, length(line_col)),
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = 0.72
            )
          }
        }
      }

      if (i == 2) {
        dev.off()
      }
    }
  }


#' plot_ssb
#'
#' @description Function the plots the mean ssb and 95% CI trends as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param spnames Species names for legend
#' @param add_ci If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param minyr First year to plot
#' @param height
#' @param width
#' @param save Save ssb?
#' @param incl_proj TRUE/FALSE, include projection years
#' @param mod_cex Cex of text for model name legend
#' @param mod_avg Vector of length Rceattle denoting if it is a model average object
#' @param mse Is if an MSE object from \code{\link{load_mse}} or \code{\link{mse_run}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#'
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_ssb <- function(Rceattle,
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
                     OM = TRUE,
                     reference = NULL) {

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    nmse = length(Rceattle)
    add_ci = TRUE
    incl_proj = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }

  # Add reference model
  if(!is.null(reference)){
    Rceattle <- c(Rceattle, list(reference))
  }

  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  if(is.null(minyr)){minyr <- min((sapply(Years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp

  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
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
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$biomassSSB[,1:nyrs_vec[i]]

    # Get SD of quantity
    if(add_ci & !mse){
      sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "biomassSSB")
      sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
      quantity_sd[,  1:nyrs_vec[i], i] <-
        replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_temp[1:(nyrs_vec[i] * nspp)])
    }

    # - Model average
    if(mod_avg[i]){
      log_quantity_sd[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples$biomassSSB[,1:nyrs_vec[i],], c(1,2), function(x) sd(as.vector(log(x))))
      log_quantity_mu[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples$biomassSSB[,1:nyrs_vec[i],], c(1,2), function(x) mean(as.vector(log(x))))
    }
  }

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){
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

  # - Rescale
  quantity <- quantity / 1000000
  quantity_upper95 <- quantity_upper95 / 1000000
  quantity_lower95 <- quantity_lower95 / 1000000
  quantity_upper50 <- quantity_upper50 / 1000000
  quantity_lower50 <- quantity_lower50 / 1000000


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_ssb_species_", i, ".csv"))
    }
  }


  # - Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * 1.2


  # - Line colors
  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }
  if(!is.null(reference)){
    line_col <- c(line_col, 1)
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_ssb_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = NA,
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      if(j == 2){
        mtext("Spawning stock biomass (million mt)", side = 2, line = 1.7, cex = 0.9)
      }

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
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


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(Years[[k]], rev(Years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(Years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = Years[[k]],
          y = quantity[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd,
          col = line_col[k]
        ) # Median
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}


#' Plot biomass eaten
#'
#' @description Function the plots the biomass consumed trends as estimated from Rceattle. Returns and saves a figure with the biomass eaten trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param minyr first year to plot
#' @param incl_proj TRUE/FALSE include projections years
#' @param incl_mean TRUE/FALSE include time series mean as horizontal line
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @export
#'
plot_b_eaten <-  function(Rceattle,
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

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    add_ci = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }


  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  if(is.null(minyr)){minyr <- min((sapply(Years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp

  minage <- Rceattle[[1]]$data_list$minage
  maxage <- max(Rceattle[[1]]$data_list$nages)
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
    array(0, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_mu <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))

  for (i in 1:length(Rceattle)) {

    # - Get quantities
    quantity[,1:nyrs_vec[i] , i] <- apply(Rceattle[[i]]$quantities$B_eaten_as_prey[,,,1:nyrs_vec[i]], c(1,4), sum)

    # # Get SD of quantity
    # # NOTE: No uncertainty estimates currently
    # if(add_ci & !mse){
    #   sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "B_eaten_as_prey")
    #   sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
    #   quantity_sd[,  1:nyrs_vec[i], i] <-
    #     replace(quantity_sd[,,,, i], values = sd_temp)
    # }

    # - Model average
    if(mod_avg[i]){
      log_quantity_sd[,1:nyrs_vec[i], i] <- apply(
        apply(Rceattle[[i]]$asymptotic_samples$B_eaten_as_prey[,,,1:nyrs_vec[i],], c(1,4), function(x) sum), # Sum across age-sex
        c(1,2), sd(as.vector(log(x)))) # SD across samples
      log_quantity_mu[,1:nyrs_vec[i], i] <- apply(
        apply(Rceattle[[i]]$asymptotic_samples$B_eaten_as_prey[,,,1:nyrs_vec[i],], c(1,4), function(x) sum), # Sum across age-sex
        c(1,2), mean(as.vector(log(x)))) # Mean across samples
    }
  }

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){

    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.975) )
    quantity_lower95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.025) )
    quantity_upper50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.75) )
    quantity_lower50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.25) )
    quantity <- apply( quantity, c(1,2), mean ) # Get mean quantity

    # -- Put back in array for indexing below
    quantity <- array(quantity, dim = c(nspp, nyrs,  1))
    quantity_upper95 <- array(quantity_upper95, dim = c(nspp, nyrs,  1))
    quantity_lower95 <- array(quantity_lower95, dim = c(nspp, nyrs,  1))
    quantity_upper50 <- array(quantity_upper50, dim = c(nspp, nyrs,  1))
    quantity_lower50<- array(quantity_lower50, dim = c(nspp, nyrs,  1))
  }

  # - Model Average
  for (i in 1:length(Rceattle)) {
    if(mod_avg[i]){
      quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
    }
  }

  # - Rescale
  quantity <- quantity / 1000000
  quantity_upper95 <- quantity_upper95 / 1000000
  quantity_lower95 <- quantity_lower95 / 1000000
  quantity_upper50 <- quantity_upper50 / 1000000
  quantity_lower50 <- quantity_lower50 / 1000000


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_b_eaten_as_prey_trajectory", i, ".csv"))
    }
  }


  ## Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * 1.2

  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_b_eaten_as_prey_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = "Biomass consumed (million t)",
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
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


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(Years[[k]], rev(Years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(Years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = Years[[k]],
          y = quantity[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd,
          col = line_col[k]
        ) # Median
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}



#' Plot biomass consumbed of each prey species by predator
#'
#' @description Function the plots the biomass consumed trends as estimated from Rceattle. Returns and saves a figure with the biomass eaten trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param minyr first year to plot
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#' @param incl_mean TRUE/FALSE include horizontal long term mean
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
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

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Species names
    if(is.null(spnames)){
      spnames =  Rceattle[[1]]$data_list$spnames
    }


    # Extract data objects
    Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    if(incl_proj){
      Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
    }

    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    if(is.null(minyr)){minyr <- min((sapply(Years, min)))}
    max_age <- max(Rceattle[[1]]$data_list$nages)

    nspp <- Rceattle[[1]]$data_list$nspp

    if(is.null(species)){
      species <- 1:nspp
    }


    # Get B_eaten
    B_eaten <-
      array(NA, dim = c(nspp * 2, nspp * 2, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      for(rsp in 1:(nspp)){
        for(ksp in 1:(nspp)){
          for(yr in 1:nyrs_vec[i]){
            B_eaten[rsp, ksp,yr,i] <- sum(Rceattle[[i]]$quantities$B_eaten[c(rsp, rsp + nspp),c(ksp, ksp + nspp),,,yr])
          }
        }
      }
    }

    ind = 1
    B_eaten <- B_eaten/1000000

    # Plot limits
    ymax <- c()
    ymin <- c()
    for (ksp in 1:dim(B_eaten)[2]) { # Loop through prey
      ymax[ksp] <- max(c(B_eaten[,ksp, , ], 0), na.rm = T)
      ymin[ksp] <- min(c(B_eaten[,ksp, , ], 0), na.rm = T)
    }
    ymax <- ymax + top_adj * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_B_eaten_trajectory", ".png")
        png(
          file = filename ,
          width = width,
          height = height,
          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(length(species) + 2), nrow = (length(species) + 2)), heights = c(0.1, rep(1, length(species)), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      for (j in 1:length(species)) {
        spp <- species[j]
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[spp], ymax[spp]),
          xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
          xlab = "Year",
          ylab = NA,
          xaxt = c(rep("n", length(species) - 1), "s")[j]
        )

        if(j == 2){
          mtext("Biomass consumed (million t)", side = 2, line = 1.6, cex = 0.8)
        }


        # Horizontal line
        if(incl_proj){
          abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", spnames[spp], bty = "n", cex = 1)

        if(!is.null(mohns)){
          legend("top", paste0("B_eaten Rho = ",  round(mohns[2,spp+1], 2) ), bty = "n", cex = 0.8) # B_eaten rho
        }

        if (j == 1) {
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

        if (j == 2 | length(species) == 1) {
          legend(
            "topright",
            legend = c("Predator:", spnames[species]),
            lty = c(NA, 1:nspp),
            lwd = lwd,
            col = c(0, rep(1, nspp)),
            bty = "n",
            cex = 1
          )
        }



        # Mean B_eaten
        for (k in 1:length(species)) {
          pred <- species[k]
          for (mod in 1:dim(B_eaten)[4]) {
            lines(
              x = Years[[mod]],
              y = B_eaten[pred, spp, 1:length(Years[[mod]]), mod],
              lwd = lwd,
              lty = k,
              col = line_col[mod]) # Median
          }

          # Average across time
          if(incl_mean){
            abline(h = mean(B_eaten[pred, spp, 1:length(Years[[1]]), ]), lwd  = lwd, col = "grey", lty = pred)
          }
        }
      }


      if (i == 2) {
        dev.off()
      }
    }
  }



#' Plot natural mortality by age
#'
#' @description Function the plots the natural mortality at age (M1 + M2) as estimated from Rceattle. Returns and saves a figure with the M-at-age trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param age Age to plot M at age
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param minyr first year to plot
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
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
           right_adj = 0,
           minyr = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Species names
    if(is.null(spnames)){
      spnames =  Rceattle[[1]]$data_list$spnames
    }


    # Extract data objects
    Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    if(incl_proj){
      Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
    }

    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    if(is.null(minyr)){minyr <- min((sapply(Years, min)))}
    nsex <- Rceattle[[1]]$data_list$nsex
    nspp <- Rceattle[[1]]$data_list$nspp

    if(is.null(species)){
      species <- 1:nspp
    }

    # Get B_eaten
    m_at_age <-
      array(NA, dim = c(nspp, 2, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      for(sp in 1:nspp){
        for(yr in 1:nyrs_vec[i]){
          m_at_age[sp, , yr, i] <- Rceattle[[i]]$quantities$M[sp,,age,yr]
        }
      }
    }

    # Plot limits
    ymax <- matrix(0, nrow = nspp, ncol = 2)
    ymin <- matrix(0, nrow = nspp, ncol = 2)
    for (i in 1:dim(m_at_age)[1]) {
      for(sex in 1:nsex[sp]){
        ymax[i,sex] <- max(c(m_at_age[i,sex,, ],0), na.rm = T)
        ymin[i,sex] <- min(c(m_at_age[i,sex,, ]), na.rm = T)
      }
    }
    ymax <- ymax + 0.15 * ymax
    ymin <- ymin - 0.15 * ymin

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_m_at_age",age,"_trajectory", ".png")
        png(
          file = filename ,
          width = width,
          height = height,
          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(sum(nsex[species]) + 2), nrow = (sum(nsex[species]) + 2)), heights = c(0.1, rep(1, sum(nsex[species])), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()
      ind = 0

      for (j in 1:length(species)) {

        spp = species[j]

        for(sex in 1:nsex[spp]){
          ind = ind+1

          # Get sex for legend
          legend_sex = sex
          legend_sex2 = ifelse(sex == 1, "female", "male")
          if(nsex[spp] == 1){
            legend_sex <- 0
            legend_sex2 = "combined"
          }

          plot(
            y = NA,
            x = NA,
            ylim = c(ymin[spp,sex], ymax[spp,sex]),
            xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
            xlab = "Year",
            ylab = paste0("M-at-age-",age),
            xaxt = c(rep("n", sum(nsex[species]) - 1), "s")[ind]
          )


          # Horizontal line
          if(incl_proj){
            abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
          }

          # Species legends
          legend("topleft", paste0(spnames[spp], " ", legend_sex2), bty = "n", cex = 1)

          # Model names legends
          if (ind == 1) {
            if(!is.null(model_names)){
              legend(
                "topright",
                legend = model_names,
                lty = rep(1, length(line_col)),
                lwd = lwd,
                col = line_col,
                bty = "n",
                cex = 0.72
              )
            }
          }

          # M-at-age
          for (k in 1:dim(m_at_age)[4]) {
            lines(
              x = Years[[k]],
              y = m_at_age[spp, sex, 1:length(Years[[k]]), k],
              lty = 1,
              lwd = lwd,
              col = line_col[k]
            ) # Median
          }


          # Average across time
          if(incl_mean){
            abline(h = mean(m_at_age[spp, sex, 1:length(Years[[1]]), ]), lwd  = lwd, col = "grey", lty = 1)
          }
        }
      }


      if (i == 2) {
        dev.off()
      }
    }
  }


#' Plot predation mortality by age and predator
#'
#' @description Function the plots the predation mortality at age (M2) by predator as estimated from Rceattle. Returns and saves a figure with the M-at-age trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param age Age to plot M at age
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param minyr first year to plot
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
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

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Species names
    if(is.null(spnames)){
      spnames =  Rceattle[[1]]$data_list$spnames
    }


    # Extract data objects
    Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    if(incl_proj){
      Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
    }

    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    if(is.null(minyr)){minyr <- min((sapply(Years, min)))}
    nsex <- Rceattle[[1]]$data_list$nsex
    nspp <- Rceattle[[1]]$data_list$nspp

    if(is.null(species)){
      species <- 1:nspp
    }

    # Get B_eaten
    m2_at_age_prop <-
      array(NA, dim = c(nspp, nspp, 2, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      for(ksp in 1:nspp){
        for(k_sex in 1:nsex[ksp]){
          for(rsp in 1:nspp){
            for(yr in 1:nyrs_vec[i]){
              m2_at_age_prop[rsp, ksp, k_sex, yr, i] <- sum(Rceattle[[i]]$quantities$M2_prop[c(rsp, rsp + nspp),ksp + (nspp*(k_sex-1)),,age,yr])
            }
          }
        }
      }
    }

    # Plot limits
    ymax <- matrix(0, nrow = nspp, ncol = 2)
    ymin <- matrix(0, nrow = nspp, ncol = 2)
    for (i in 1:nspp) {
      for(sex in 1:nsex[i]){
        ymax[i,sex] <- max(c(m2_at_age_prop[,i,sex,, ], 0), na.rm = T)
        ymin[i,sex] <- min(c(m2_at_age_prop[,i,sex,, ], 0), na.rm = T)
      }
    }
    ymax <- ymax + top_adj * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_m2_at_age_prop",age,"_trajectory", ".png")
        png(
          file = filename ,
          width = width,
          height = height,
          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(sum(nsex[species]) + 2), nrow = (sum(nsex[species]) + 2)), heights = c(0.1, rep(1, sum(nsex[species])), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()
      ind = 0

      for (j in 1:length(species)) {

        spp = species[j]

        for(sex in 1:nsex[spp]){

          ind = ind+1

          # Get sex for legend
          legend_sex = sex
          legend_sex2 = ifelse(sex == 1, "Female", "Male")
          if(nsex[spp] == 1){
            legend_sex <- 0
            legend_sex2 = "Combined"
          }

          plot(
            y = NA,
            x = NA,
            ylim = c(ymin[spp,sex], ymax[spp,sex]),
            xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
            xlab = "Year",
            ylab = paste0("M-at-age-",age),
            xaxt = c(rep("n", sum(nsex[species]) - 1), "s")[ind]
          )

          # Horizontal line
          if(incl_proj){
            abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
          }

          # Species legends
          legend("topleft", paste0(spnames[spp], " ", legend_sex2), bty = "n", cex = 1)

          # Model names legends
          if (ind == 1) {
            if(!is.null(model_names)){
              legend(
                "topright",
                legend = model_names,
                lty = rep(1, length(line_col)),
                lwd = lwd,
                col = line_col,
                bty = "n",
                cex = 0.72
              )
            }
          }

          # Predator legend
          if (ind == 2 | length(species) == 1) {
            legend(
              "topright",
              legend = c("Predator:", spnames),
              lty = c(NA, 1:nspp),
              lwd = lwd,
              col = c(0, rep(1, nspp)),
              bty = "n",
              cex = 1
            )
          }

          # M-at-age
          for (k in 1:dim(m2_at_age_prop)[5]) {
            for(rsp in 1:nspp){
              lines(
                x = Years[[k]],
                y = m2_at_age_prop[rsp, spp, sex, 1:length(Years[[k]]), k],
                lty = rsp,
                lwd = lwd,
                col = line_col[k]
              ) # Median

              # Average across time
              if(incl_mean){
                abline(h = mean(m2_at_age_prop[rsp, spp, sex, 1:length(Years[[1]]), ]), lwd  = lwd, col = "grey", lty = rsp)
              }
            }
          }
        }
      }


      if (i == 2) {
        dev.off()
      }
    }
  }

#' plot_depletionSSB
#'
#' @description Function the plots the mean depletion of spawning stock biomass (SSB) as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param spnames Species names for legend
#' @param add_ci NOT WORKING If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param top_adj Multiplier for y-axis to add to the top side of the figure for fitting the legend.
#' @param minyr First year to plot
#' @param height
#' @param width
#' @param save Save output?
#' @param incl_proj TRUE/FALSE, include projection years
#' @param mod_cex Cex of text for model name legend
#' @param mse Is if an MSE object from \code{\link{load_mse}} or \code{\link{mse_run}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#'
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_depletionSSB <- function(Rceattle,
                              file = NULL,
                              model_names = NULL,
                              line_col = NULL,
                              species = NULL,
                              spnames = NULL,
                              add_ci = FALSE,
                              lwd = 3,
                              save = FALSE,
                              right_adj = 0,
                              top_adj = 1.2,
                              width = 7,
                              height = 6.5,
                              minyr = NULL,
                              incl_proj = FALSE,
                              mod_cex = 1,
                              alpha = 0.4,
                              mod_avg = rep(FALSE, length(Rceattle)),
                              mse = FALSE,
                              OM = TRUE,
                              reference = NULL) {

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
      incl_proj = TRUE
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
      incl_proj = FALSE
    }
    nmse = length(Rceattle)
    add_ci = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }

  # Add reference model
  if(!is.null(reference)){
    Rceattle <- c(Rceattle, list(reference))
  }


  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  if(is.null(minyr)){minyr <- min((sapply(Years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp

  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_mu <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  ptarget = matrix(NA, nrow = length(Rceattle), ncol = nspp)
  plimit = matrix(NA, nrow = length(Rceattle), ncol = nspp)

  for (i in 1:length(Rceattle)) {
    ptarget[i,] <- Rceattle[[i]]$data_list$Ptarget
    plimit[i,] <- Rceattle[[i]]$data_list$Plimit
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$depletionSSB[,1:nyrs_vec[i]]


    # Get SD of quantity
    if(add_ci & !mse){
      sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "depletionSSB")
      sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
      quantity_sd[,  1:nyrs_vec[i], i] <-
        replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_temp[1:(nyrs_vec[i] * nspp)])
    }
  }
  quantity[is.infinite(quantity)] <- NA # for dynamicB0

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){
    ptarget <- ptarget[1,]
    plimit <- plimit[1,]

    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.975, na.rm = TRUE) )
    quantity_lower95 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.025, na.rm = TRUE) )
    quantity_upper50 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.75, na.rm = TRUE) )
    quantity_lower50 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.25, na.rm = TRUE) )

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

  # # - Model Average
  # for (i in 1:length(Rceattle)) {
  #   if(mod_avg[i]){
  #     quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #     quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #     quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #   }
  # }


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_depletionssb_species_", i, ".csv"))
    }
  }


  ## Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * top_adj

  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }
  if(!is.null(reference)){
    line_col <- c(line_col, 1)
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_depletionssb_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = "SSB depletion",
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
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


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(Years[[k]], rev(Years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(Years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = Years[[k]],
          y = quantity[spp[j], 1:length(Years[[k]]), k],
          lty = ifelse(is.null(reference) == FALSE & mse == TRUE, c(1,2)[k], 1),
          lwd = lwd,
          col = line_col[k]
        ) # Median

        # Ptarget and plimit
        abline(h= ptarget[spp[j]], lwd = lwd/2, col = "blue")
        abline(h= plimit[spp[j]], lwd = lwd/2, col = "red")
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}


#' plot_depletion
#'
#' @description Function the plots the mean depletion of total biomass as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param spnames Species names for legend
#' @param add_ci NOT WORKING If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param minyr First year to plot
#' @param height
#' @param width
#' @param save Save depletion?
#' @param incl_proj TRUE/FALSE, include projection years
#' @param mod_cex Cex of text for model name legend
#'
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_depletion <- function(Rceattle,
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

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    add_ci = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }


  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  if(is.null(minyr)){minyr <- min((sapply(Years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp

  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_mu <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  ptarget = matrix(NA, nrow = length(Rceattle), ncol = nspp)
  plimit = matrix(NA, nrow = length(Rceattle), ncol = nspp)

  for (i in 1:length(Rceattle)) {
    ptarget[i,] <- Rceattle[[i]]$data_list$Ptarget
    plimit[i,] <- Rceattle[[i]]$data_list$Plimit
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$depletion[,1:nyrs_vec[i]]


    # Get SD of quantity
    if(add_ci & !mse){
      sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "depletion")
      sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
      quantity_sd[,  1:nyrs_vec[i], i] <-
        replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_temp[1:(nyrs_vec[i] * nspp)])
    }
  }

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){
    ptarget <- ptarget[1,]
    plimit <- plimit[1,]

    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.975) )
    quantity_lower95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.025) )
    quantity_upper50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.75) )
    quantity_lower50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.25) )
    quantity <- apply( quantity, c(1,2), mean ) # Get mean quantity

    # -- Put back in array for indexing below
    quantity <- array(quantity, dim = c(nspp, nyrs,  1))
    quantity_upper95 <- array(quantity_upper95, dim = c(nspp, nyrs,  1))
    quantity_lower95 <- array(quantity_lower95, dim = c(nspp, nyrs,  1))
    quantity_upper50 <- array(quantity_upper50, dim = c(nspp, nyrs,  1))
    quantity_lower50<- array(quantity_lower50, dim = c(nspp, nyrs,  1))
  }

  # # - Model Average
  # for (i in 1:length(Rceattle)) {
  #   if(mod_avg[i]){
  #     quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #     quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #     quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #   }
  # }


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_depletion_species_", i, ".csv"))
    }
  }


  ## Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * 1.2

  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_depletion_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = "Biomass depletion",
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
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


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(Years[[k]], rev(Years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(Years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = Years[[k]],
          y = quantity[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd,
          col = line_col[k]
        ) # Median

        # Ptarget and plimit
        abline(h= ptarget[spp[j]], lwd = lwd/2, col = "blue")
        abline(h= plimit[spp[j]], lwd = lwd/2, col = "red")
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}


#' plot F
#'
#' @description Function the plots the F time series per species from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param spnames Species names for legend
#' @param add_ci NOT WORKING If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj Multiplier for to add to the right side of the figure for fitting the legend.
#' @param minyr First year to plot
#' @param height
#' @param width
#' @param save Save output?
#' @param incl_proj TRUE/FALSE, include projection years
#' @param mod_cex Cex of text for model name legend
#' @param mse Is if an MSE object from \code{\link{load_mse}} or \code{\link{mse_run}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#'
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

  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    }
    if(!OM){
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    add_ci = TRUE
  }


  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }


  # Species names
  if(is.null(spnames)){
    spnames =  Rceattle[[1]]$data_list$spnames
  }

  # Extract data objects
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  if(is.null(minyr)){minyr <- min((sapply(Years, min)))}

  nspp <- Rceattle[[1]]$data_list$nspp
  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics


  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species


  # Get depletion
  quantity <- quantity_sd <- log_quantity_sd <- log_quantity_mu <-  ftarget <- flimit <- array(NA, dim = c(nspp, nyrs,  length(Rceattle)))

  for (i in 1:length(Rceattle)) {
    ftarget[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$Ftarget[,1:nyrs_vec[i]]
    flimit[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$Flimit[,1:nyrs_vec[i]]
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities$F_spp[,1:nyrs_vec[i]]


    # Get SD of quantity
    if(add_ci & !mse){
      sd_temp <- which(names(Rceattle[[i]]$sdrep$value) == "F_spp")
      sd_temp <- Rceattle[[i]]$sdrep$sd[sd_temp]
      quantity_sd[,  1:nyrs_vec[i], i] <-
        replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_temp[1:(nyrs_vec[i] * nspp)])
    }
  }

  ## Get confidence intervals
  # - Single model
  if(!mse){
    quantity_upper95 <- quantity + quantity_sd * 1.92
    quantity_lower95 <- quantity - quantity_sd * 1.92

    quantity_upper50 <- quantity + quantity_sd * 0.674
    quantity_lower50 <- quantity - quantity_sd * 0.674
  }

  # - MSE objects
  if(mse){

    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.975) )
    quantity_lower95 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.025) )
    quantity_upper50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.75) )
    quantity_lower50 <- apply( quantity, c(1,2), function(x) quantile(x, probs = 0.25) )
    quantity <- apply( quantity, c(1,2), mean ) # Get mean quantity

    # -- Put back in array for indexing below
    quantity <- array(quantity, dim = c(nspp, nyrs,  1))
    quantity_upper95 <- array(quantity_upper95, dim = c(nspp, nyrs,  1))
    quantity_lower95 <- array(quantity_lower95, dim = c(nspp, nyrs,  1))
    quantity_upper50 <- array(quantity_upper50, dim = c(nspp, nyrs,  1))
    quantity_lower50<- array(quantity_lower50, dim = c(nspp, nyrs,  1))
  }

  # # - Model Average
  # for (i in 1:length(Rceattle)) {
  #   if(mod_avg[i]){
  #     quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #     quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #     quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i]) / 1000000
  #   }
  # }


  ## Save
  if (save) {
    for (i in 1:nspp) {
      dat <- data.frame(quantity[i, , ])
      datup <- data.frame(quantity_upper95[i, , ])
      datlow <- data.frame(quantity_lower95[i, , ])

      dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
      colnames(dat_new) <- rep(model_names[1], 3)

      for (j in 2:ncol(dat)) {
        dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
        colnames(dat_new2) <- rep(model_names[j], 3)
        dat_new <- cbind(dat_new, dat_new2)

      }

      write.csv(dat_new, file = paste0(file, "_f_species_", i, ".csv"))
    }
  }


  ## Plot limits
  ymax <- c()
  ymin <- c()
  for (sp in 1:nspp) {
    if (add_ci & (estDynamics[sp] == 0)) {
      ymax[sp] <- max(c(quantity_upper95[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity_upper95[sp, , ], 0), na.rm = T)
    } else{
      ymax[sp] <- max(c(quantity[sp, , ], 0), na.rm = T)
      ymin[sp] <- min(c(quantity[sp, , ], 0), na.rm = T)
    }
  }
  ymax <- ymax * 1.2

  if (is.null(line_col)) {
    if(!mse){
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }
    if(mse){
      line_col <- 1
    }
  }


  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (i in 1:loops) {
    if (i == 2) {
      filename <- paste0(file, "_f_trajectory", ".png")
      png(
        file = filename ,
        width = width,# 169 / 25.4,
        height = height,# 150 / 25.4,
        units = "in",
        res = 300
      )
    }

    # Plot configuration
    layout(matrix(1:(length(spp) + 2), nrow = (length(spp) + 2)), heights = c(0.1, rep(1, length(spp)), 0.2))
    par(
      mar = c(0, 3 , 0 , 1) ,
      oma = c(0 , 0 , 0 , 0),
      tcl = -0.35,
      mgp = c(1.75, 0.5, 0)
    )
    plot.new()

    for (j in 1:length(spp)) {
      plot(
        y = NA,
        x = NA,
        ylim = c(ymin[spp[j]], ymax[spp[j]]),
        xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
        xlab = "Year",
        ylab = "Fishing mortality",
        xaxt = c(rep("n", length(spp) - 1), "s")[j]
      )

      # Horizontal line at end yr
      if(incl_proj){
        abline(v = Rceattle[[length(Rceattle)]]$data_list$meanyr, lwd  = lwd, col = "grey", lty = 2)
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


      # Credible interval
      if(estDynamics[spp[j]] == 0){
        if (add_ci) {
          for (k in 1:dim(quantity)[3]) {
            # - 95% CI
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(quantity_upper95[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower95[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = alpha/2),
              border = NA
            )

            # - 50% CI
            if(mse){
              polygon(
                x = c(Years[[k]], rev(Years[[k]])),
                y = c(quantity_upper50[spp[j], 1:length(Years[[k]]), k], rev(quantity_lower50[spp[j], 1:length(Years[[k]]), k])),
                col = adjustcolor( line_col[k], alpha.f = alpha),
                border = NA
              )
            }
          }
        }
      }

      # Mean quantity
      for (k in 1:dim(quantity)[3]) {
        lines(
          x = Years[[k]],
          y = quantity[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd,
          col = line_col[k]
        ) # Median

        # - Ftarget
        lines(
          x = Years[[k]],
          y = ftarget[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd/2,
          col = "blue"
        )

        # - Flimit
        lines(
          x = Years[[k]],
          y = flimit[spp[j], 1:length(Years[[k]]), k],
          lty = 1,
          lwd = lwd/2,
          col = "red"
        )
      }
    }


    if (i == 2) {
      dev.off()
    }
  }
}
