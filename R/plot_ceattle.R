
#' plot_biomass
#'
#' Function the plots the biomass and spawning stock biomass trends as estimated from Rceattle
#'
#' @param file_name name of a file to identified the files exported by the
#'   function.
#' @param ceattle_list List of CEATTLE model objects
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param include_srv Boolian of whether to include survey biomass estimates and 95% CI
#'
#' @return Returns and saves a figure with the population trajectory.
plot_biomass <- function(ceattle_list, file_name = "NULL", model_names = NULL, line_col = NULL, species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"), lwd = 3, include_srv = FALSE) {
  library(extrafont)

  # Extract data objects
  nyrs <- ceattle_list[[1]]$data_list$nyrs
  Years <- ceattle_list[[1]]$data_list$styr : (ceattle_list[[1]]$data_list$styr + nyrs - 1)
  nspp <- ceattle_list[[1]]$data_list$nspp

  # Get biomass
  Biomass <- array(NA, dim = c(nspp, nyrs, length(ceattle_list) ))
  for(i in 1:length(ceattle_list)){
    Biomass[,,i] <- ceattle_list[[i]]$quantities$biomass
  }
  Biomass <- Biomass / 1000000

  # Get SSB
  SSB <- array(NA, dim = c(nspp, nyrs, length(ceattle_list) ))
  for(i in 1:length(ceattle_list)){
    SSB[,,i] <- ceattle_list[[i]]$quantities$biomassSSB
  }
  SSB <- SSB / 1000000

  # Plot limits
  ymax <- c()
  for(i in 1:dim(Biomass)[1]){
    ymax[i] <- max(Biomass[i,,], na.rm = T)
  }
  ymax <- ymax + 0.15 * ymax
  ymin <- 0

  if(is.null(line_col)){
    line_col <- 1:length(ceattle_list)
  }


  # Plot trajectory
  for(i in 1:2){
    if(i == 1){
      filename <- paste0(file_name, "_biomass_trajectory", ".tiff")
      tiff( file = filename , width=169 / 25.4, height = 150 / 25.4, family = "Helvetica" , units = "in", res = 300)
    }

    # Plot configuration
    layout(matrix(1:(nspp+2), nrow = (nspp+2)), heights = c(0.1, rep(1, nspp), 0.2))
    par(mar=c(0, 3 , 0 , 1) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
    plot.new()

    for(j in 1:nspp){
      plot(y = NA, x = NA,
           ylim = c(ymin, ymax[j]),
           xlim = c(min(Years), max(Years)),
           xlab = "Year", ylab = "Biomass (million t)",
           xaxt = c(rep("n", nspp - 1), "s")[j])

      # Legends
      legend("topleft", species[j], bty = "n", cex = 1.4)

      if(j == 1){
        legend("topright", model_names, lty = rep(1, length(line_col)), lwd = lwd, col = line_col, bty = "n", cex = 1.175)
      }

      if(j == 2){
        legend("topright", c("Biomass", "SSB"), lty = c(1, 2), lwd = lwd, col = c(1, 1), bty = "n", cex = 1.175)
      }

      # Survey data
      if(include_srv){
        srv_yrs <- ceattle_list[[1]]$data_list$yrs_srv_biom[j,]
        srv_biom <- ceattle_list[[1]]$data_list$srv_biom[j,]
        points( x = srv_yrs,
                y =  srv_biom / 1000000,
                pch = 16, cex = 1, col = "#696773")


        # Get 95% CI and range
        Upper95 <- (ceattle_list[[1]]$data_list$srv_biom[j,] + ceattle_list[[1]]$data_list$srv_biom_se[j,] * 1.92)/1000000
        Lower95 <- (ceattle_list[[1]]$data_list$srv_biom[j,] - ceattle_list[[1]]$data_list$srv_biom_se[j,] * 1.92)/1000000


        arrows( x0 = srv_yrs,
                y0 = Upper95,
                x1 = srv_yrs,
                y1 = Lower95,
                length=0.05, angle=90, code=3, lwd = 2, col = "#696773")


        # Mean biomass
        for(k in 1:dim(Biomass)[3]){
          lines( x = Years, y = Biomass[j,,k], lty = 1, lwd = lwd, col = line_col[k]) # Median
          lines( x = Years, y = SSB[j,,k], lty = 2, lwd = lwd, col = line_col[k]) # Median
        }

        # # Credible interval
        # polygon(
        #   x = c(Years,rev(Years)),
        #   y = c(output_summary[3, ],rev(output_summary[4, ])),
        #   col = "Grey80", border = NA) # 95% CI
        # polygon( x = c(Years,rev(Years)),
        #          y = c(output_summary[5, ], rev(output_summary[6, ])),
        #          col = "Grey60", border = NA) # 90% CI
      }
    }


    if(i == 1){ dev.off()}
  }
}


#' plot_recruitment
#'
#' Function the plots the mean recruitment and 95% CI trends as estimated from Rceattle
#'
#' @param file_name name of a file to identified the files exported by the
#'   function.
#' @param ceattle_list List of CEATTLE model objects
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param ci_col Colors to be used for CI color
#' @param lwd Line width as specified by user
#'
#' @return Returns and saves a figure with the population trajectory.
plot_recruitment <- function(ceattle_list, file_name = "NULL", model_names = NULL, line_col = NULL, species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"), ci_col = NULL, lwd = 3) {


  library(extrafont)

  # Extract data objects
  nyrs <- ceattle_list[[1]]$data_list$nyrs
  Years <- ceattle_list[[1]]$data_list$styr : (ceattle_list[[1]]$data_list$styr + nyrs - 1)
  nspp <- ceattle_list[[1]]$data_list$nspp


  # Get biomass
  recruitment <- array(NA, dim = c(nspp, nyrs, length(ceattle_list) ))
  recruitment_sd <- array(NA, dim = c(nspp, nyrs, length(ceattle_list) ))
  for(i in 1:length(ceattle_list)){
    recruitment[,,i] <- t(ceattle_list[[i]]$quantities$NByage[,1,])

    # Get SD of rec
    sd_rec <- which(names(ceattle_list[[i]]$sdrep$value) == "R")
    sd_rec <- ceattle_list[[i]]$sdrep$sd[sd_rec]
    recruitment_sd[,,i] <- replace(recruitment_sd[,,i], values = sd_rec)
  }


  recruitment <- recruitment / 1000000
  recruitment_sd <- recruitment_sd / 1000000
  recruitment_upper <- recruitment + recruitment_sd * 1.92
  recruitment_lower <- recruitment - recruitment_sd * 1.92

  # Plot limits
  ymax <- c()
  for(i in 1:dim(recruitment)[1]){
    ymax[i] <- max(recruitment_upper[i,,], na.rm = T)
  }
  ymax <- ymax + 0.2 * ymax
  ymin <- 0

  if(is.null(line_col)){
    line_col <- 1:length(ceattle_list)
  }

  if(is.null(ci_col)){
    ci_col <- rep(NA, length(ceattle_list))
  }


  # Plot trajectory
  for(i in 1:2){
    if(i == 1){
      filename <- paste0(file_name, "_recruitment_trajectory", ".tiff")
      tiff( file = filename , width=169 / 25.4, height = 150 / 25.4, family = "Helvetica", units = "in", res = 300)
    }

    # Plot configuration
    layout(matrix(1:(nspp+2), nrow = (nspp+2)), heights = c(0.1, rep(1, nspp), 0.2))
    par(mar=c(0, 3 , 0 , 1) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
    plot.new()

    for(j in 1:nspp){
      plot(y = NA, x = NA,
           ylim = c(ymin, ymax[j]),
           xlim = c(min(Years), max(Years)),
           xlab = "Year", ylab = "Recruitment (millions)", xaxt = c(rep("n", nspp - 1), "s")[j])

      # Legends
      legend("topleft", species[j], bty = "n", cex = 1.4)
      if(j == 1){
        legend("topright", model_names, lty = rep(1, length(line_col)), lwd = lwd, col = line_col, bty = "n", cex = 1.175)
      }


      # Credible interval
      for(k in 1:dim(recruitment)[3]){
        polygon(
          x = c(Years,rev(Years)),
          y = c(recruitment_upper[j,,k], rev(recruitment_lower[j,,k])),
          col = ci_col[k], border = NA, density = 15 * k, angle = 22 * k) # 95% CI
      }

      # Mean biomass
      for(k in 1:dim(recruitment)[3]){
        lines( x = Years, y = recruitment[j,,k], lty = 1, lwd = lwd, col = line_col[k]) # Median
      }
    }


    if(i == 1){ dev.off()}
  }
}



#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated indices of abundance a SIR  model including: median, 95%
#' credible interval, 90% credible interval, catch, and observed indices of abundance abundance.
#'
#' @param SIR A fit SIR model
#' @param file_name name of a file to identified the files exported by the
#'   function.
#'
#' @return Returns and saves a figure with the IOA trajectories.
plot_ioa <- function(SIR, file_name = "NULL"){

  rel.abundance <- SIR$inputs$rel.abundance
  row_names <- c("mean", "median",
                 "2.5%PI", "97.5%PI",
                 "5%PI", "95%PI",
                 "min", "max", "n")

  if(is.null(rel.abundance)){
    stop("SIR model did not include an IOA")
  }

  rel.abundance$Upper95 <- qlnorm(0.975, log(rel.abundance$IA.obs), rel.abundance$Sigma)
  rel.abundance$Lower95 <- qlnorm(0.025, log(rel.abundance$IA.obs), rel.abundance$Sigma)

  # Plot IOA

  IA.yrs <- rel.abundance$Year
  IA.yr.range <- c((min(IA.yrs) - 1):(max(IA.yrs) + 1)) # Range +- 1 of IOA years

  # Predict IOA
  N_hat <- SIR$resamples_trajectories[, paste0("N_", IA.yr.range)] # Estimates of N within IOA years
  q_cols <- grep("q_IA", colnames(SIR$resamples_output)) # Columns of resample Q estimates
  q_est <- SIR$resamples_output[, q_cols]


  IA_pread <- list()
  IA_summary <- list()
  for(i in 1:length(q_cols)){
    # Predict
    IA_pread[[i]] <- N_hat * q_est

    # Summarize
    IA_summary[[i]] <-  matrix(nrow = length(row_names), ncol = dim(IA_pread[[i]])[2])
    IA_summary[[i]][1, ] <- sapply(IA_pread[[i]], mean)
    IA_summary[[i]][2:6, ] <- sapply(IA_pread[[i]], quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
    IA_summary[[i]][7, ] <- sapply(IA_pread[[i]], min)
    IA_summary[[i]][8, ] <- sapply(IA_pread[[i]], max)
    IA_summary[[i]][9, ] <- sapply(IA_pread[[i]], length)



    IA_summary[[i]] <- data.frame(IA_summary[[i]])
    names(IA_summary[[i]]) <- paste0("IA", i, "_", IA.yr.range)
    row.names(IA_summary[[i]]) <- row_names
  }



  ymax <- max(c(sapply(IA_summary, function(x) max(x[2:6,])), rel.abundance$Lower95, rel.abundance$Upper95))
  ymin <- 0


  for(j in 1:2){
    if(j == 1){
      filename <- paste0(file_name, "_IOA_fit", ".tiff")
      tiff( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
    }

    par(mfrow = c(1,length(IA_summary)))

    # Loop through inices
    for(i in 1:length(IA_summary)){
      rel.abundance.sub <- rel.abundance[which(rel.abundance$Index == i),]

      # Plot configuration
      par( mar=c(3, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
      plot(y = NA, x = NA,
           ylim = c(ymin, ymax),
           xlim = c(min(IA.yr.range), max(IA.yr.range)),
           xlab = "Year", ylab = "Relative abundance")


      # Credible interval
      polygon(
        x = c(IA.yr.range, rev(IA.yr.range)),
        y = c(IA_summary[[i]][3, ],rev(IA_summary[[i]][4, ])),
        col = "Grey80", border = NA) # 95% CI

      polygon( x = c(IA.yr.range, rev(IA.yr.range)),
               y = c(IA_summary[[i]][5, ], rev(IA_summary[[i]][6, ])),
               col = "Grey60", border = NA) # 90% CI


      # Median and catch series
      lines( x = IA.yr.range, y = IA_summary[[i]][2, ], lwd = 3) # Median


      # Relative abundance
      points( x = rel.abundance.sub$Year,
              y = rel.abundance.sub$IA.obs,
              col = 1, pch = 16, cex = 2)
      arrows( x0 = rel.abundance.sub$Year,
              y0 = rel.abundance.sub$Lower95,
              x1 = rel.abundance.sub$Year,
              y1 = rel.abundance.sub$Upper95,
              length=0.05, angle=90, code=3, lwd = 3, col = 1)

    }
    if(j == 1){ dev.off()}
  }
}


#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated posterior densities of parameters from  SIR  model.
#'
#' @param SIR A fit SIR model or list of SIR models. Plots in the order provided.
#' @param file_name name of a file to identified the files exported by the
#'   function.
#' @param multiple_sirs Logical whether or not multiple SIRS are provided as a list.
#'
#' @return Returns and saves a figure with the posterior densities of parameters.
plot_density <- function(SIR, file_name = "NULL", multiple_sirs = FALSE){

  if(multiple_sirs == FALSE){
    sir_list <- list(SIR)
  }
  if(multiple_sirs == TRUE){
    sir_list <- SIR
  }

  # Vars of interest
  years <- c( sir_list[[1]]$inputs$target.Yr, sir_list[[1]]$inputs$output.Years)
  vars <- c("r_max", "K", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))
  vars_latex <- c("$r_{max}$", "$K$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years))

  # Plot
  for(j in 1:2){
    if(j == 1){
      filename <- paste0(file_name, "_posterior_density", ".tiff")
      tiff( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
    }

    par(mfrow = c(2,length(vars)/2))
    par( mar=c(3, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

    # Loop through vars
    for(i in 1:length(vars)){

      # Extract posterio densities
      posterior_dens <- list()
      for(k in 1:length(sir_list)){
        posterior_dens[[k]] <- density(sir_list[[k]]$resamples_output[,vars[i]])
      }
      posteriors_lwd <- rep(3, length(posterior_dens))
      posteriors_lty <- c(1:length(posterior_dens))

      # Extract prior densities
      # prior_dens <- list()
      # for(k in 1:length(sir_list)){
      #     prior_dens[[k]] <- density(sir_list[[k]]$resamples_output[,vars[i]])
      # }
      # priors_lwd <- rep(3, length(prior_dens))
      # priors_lty <- c(1:length(prior_dens))

      # Plot them
      plot(NA,
           xlim = quantile(sapply(posterior_dens, "[", "x")$x, probs= c(0.02, 0.95)),
           ylim = range(sapply(posterior_dens, "[", "y")),
           ylab = "Density", xlab = latex2exp::TeX(vars_latex[i]))
      mapply(lines, posterior_dens, lwd = posteriors_lwd, lty = posteriors_lty)
    }

    if(j == 1){ dev.off()}
  }
}
