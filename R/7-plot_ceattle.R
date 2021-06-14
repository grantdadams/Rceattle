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
#' @description Function the plots the biomass and spawning stock biomass trends as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param species Which species to plot e.g. c(1,4). Default = NULL plots them all
#' @param lwd Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#'
#' @return Returns and saves a figure with the population trajectory.
#' @export
plot_biomass <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           spnames = NULL,
           lwd = 3,
           right_adj = 0,
           mohns = NULL,
           incl_proj = FALSE) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Species names
    if(is.null(spnames)){
      spnames =  Rceattle[[1]]$data_list$spnames
    }


    # Extract data objects
    Years <- list()
    Endyrs <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    if(is.null(species)){
      nspp <- Rceattle[[1]]$data_list$nspp
      species <- 1:nspp
    }

    nspp <- length(species)

    # Get biomass
    Biomass <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      Biomass[species - min(species) + 1, 1:length(Years[[i]]), i] <- Rceattle[[i]]$quantities$biomass[species,1:nyrs_vec[i]]
    }
    Biomass <- Biomass / 1000000

    # Get SSB
    SSB <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      SSB[species - min(species) + 1, 1:length(Years[[i]]), i] <- Rceattle[[i]]$quantities$biomassSSB[species,1:nyrs_vec[i]]
    }

    SSB <- SSB / 1000000

    # Plot limits
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(Biomass)[1]) {
      ymax[i] <- max(c(Biomass[i, , ], SSB[i, , ], 0), na.rm = T)
      ymin[i] <- min(c(Biomass[i, , ], SSB[i, , ], 0), na.rm = T)
    }
    ymax <- ymax + 0.15 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_biomass_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(nspp + 2), nrow = (nspp + 2)), heights = c(0.1, rep(1, nspp), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      for (j in 1:nspp) {
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[j], ymax[j]),
          xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
          xlab = "Year",
          ylab = "Biomass (million t)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", spnames[species[j]], bty = "n", cex = 0.8)

        if(!is.null(mohns)){
          legend("top", paste0("B Rho = ", round(mohns[1,j+1], 2), "; SSB Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 0.8) # Biomass rho
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
              cex = 0.8
            )
          }
        }

        if (j == 2) {
          legend(
            "topright",
            legend = c("Biomass", "SSB"),
            lty = c(1, 2),
            lwd = lwd,
            col = c(1, 1),
            bty = "n",
            cex = 0.8
          )
        }


        # Mean biomass
        for (k in 1:dim(Biomass)[3]) {
          lines(
            x = Years[[k]],
            y = Biomass[j, 1:length(Years[[k]]), k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
          lines(
            x = Years[[k]],
            y = SSB[j, 1:length(Years[[k]]), k],
            lty = 2,
            lwd = lwd,
            col = line_col[k]
          ) # Median
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
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE, include projection years
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_recruitment <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           spnames = NULL,
           add_ci = FALSE,
           lwd = 3,
           save_rec = FALSE,
           right_adj = 0,
           mohns = NULL,
           incl_proj = FALSE) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Species names
    if(is.null(spnames)){
      spnames =  Rceattle[[1]]$data_list$spnames
    }

    # Extract data objects
    Years <- list()
    Endyrs <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    spp <- which(Rceattle[[1]]$data_list$estDynamics == 0)
    nspp <- Rceattle[[1]]$data_list$nspp
    minage <- Rceattle[[1]]$data_list$minage

    if(is.null(species)){

      species <- 1:nspp
    }

    spp <- spp[which(spp %in% species)]


    # Get biomass
    recruitment <-
      array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
    recruitment_sd <-
      array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      recruitment[, 1:length(Years[[i]]) , i] <- Rceattle[[i]]$quantities$R[,1:nyrs_vec[i]]

      # Get SD of rec
      if (add_ci) {
        sd_rec <- which(names(Rceattle[[i]]$sdrep$value) == "R")
        sd_rec <- Rceattle[[i]]$sdrep$sd[sd_rec]
        recruitment_sd[, , i] <-
          replace(recruitment_sd[, , i], values = sd_rec[1:(nyrs_vec[i] * nspp)])
      }
    }


    recruitment <- recruitment / 1000000
    recruitment_sd <- recruitment_sd / 1000000
    recruitment_upper <- recruitment + recruitment_sd * 1.92
    recruitment_lower <- recruitment - recruitment_sd * 1.92

    if (save_rec) {
      for (i in 1:nspp) {
        dat <- data.frame(recruitment[i, , ])
        datup <- data.frame(recruitment_upper[i, , ])
        datlow <- data.frame(recruitment_lower[i, , ])

        dat_new <- cbind(dat[, 1], datlow[, 1], datup[, 1])
        colnames(dat_new) <- rep(model_names[1], 3)

        for (j in 2:ncol(dat)) {
          dat_new2 <- cbind(dat[, j], datlow[, j], datup[, j])
          colnames(dat_new2) <- rep(model_names[j], 3)
          dat_new <- cbind(dat_new, dat_new2)

        }


        filename <-
          paste0(file, "_recruitment_species_", i, ".csv")
        write.csv(dat_new, file = filename)
      }
    }


    # Plot limits
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(recruitment)[1]) {
      if (add_ci) {
        ymax[i] <- max(c(recruitment_upper[i, , ], 0), na.rm = T)
        ymin[i] <- min(c(recruitment_upper[i, , ], 0), na.rm = T)
      } else{
        ymax[i] <- max(c(recruitment[i, , ], 0), na.rm = T)
        ymin[i] <- min(c(recruitment[i, , ], 0), na.rm = T)
      }
    }
    ymax <- ymax + 0.2 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_recruitment_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,
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
          ylab = "Recruitment (millions)",
          xaxt = c(rep("n", length(spp) - 1), "s")[j]
        )

        # Horizontal line at end yr
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft",
               legend = spnames[spp[j]],
               bty = "n",
               cex = 0.8)

        if(!is.null(mohns)){
          legend("top", paste0("Rho = ", round(mohns[3,spp[j]+1], 2) ), bty = "n", cex = 0.8) # Biomass rho
        }

        if (spp[j] == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              lty = rep(1, length(line_col)),
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = 0.8
            )
          }

        }


        # Credible interval
        if (add_ci) {
          for (k in 1:dim(recruitment)[3]) {
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(recruitment_upper[spp[j], 1:length(Years[[k]]), k], rev(recruitment_lower[spp[j], 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = 0.4),
              border = NA
            ) # 95% CI
          }
        }

        # Mean recruitment
        for (k in 1:dim(recruitment)[3]) {
          lines(
            x = Years[[k]],
            y = recruitment[spp[j], 1:length(Years[[k]]), k],
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
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    Years <- list()
    for(i in 1:length(Rceattle)){
      Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
    }
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
              width = 7,# 169 / 25.4,
              height = 6.5, # 150 / 25.4,

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
                width = 7,# 169 / 25.4,
                height = 6.5, # 150 / 25.4,

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
#' @param contour If plot it to be done as contours rather than tiles.
#' @param species Species names for legend
#' @param width Plot width when saved "inches"
#' @param height Plot height when saved "inches"
#' @param title Additional title to add. Will also add species names if not NULL
#' @param title_cex Font size for title
#' @param spp Species to plot. Plots all if null.
#' @param log TRUE/FALSE use log M1 + M2
#' @param maxage Plot up to this age. Plots all ages if NULL
#'
#' @export
plot_mortality <-
  function(Rceattle,
           file = NULL,
           incl_proj = FALSE,
           zlim = NULL,
           contour = FALSE,
           width = NULL,
           height = NULL,
           title = NULL,
           log = FALSE,
           spp = NULL,
           maxage = NULL,
           title_cex = 10) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    if(length(Rceattle) > 1){
      stop("Can only plot one model")
    }

    # Extract data objects
    Years <- list()
    for(i in 1:length(Rceattle)){
      Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp
    spnames <- Rceattle[[1]]$data_list$spnames
    estdynamics <- Rceattle[[1]]$data_list$estDynamics
    nages <- Rceattle[[1]]$data_list$nages
    for(i in 1:length(nages)){
      nages[i] <- ifelse(nages[i] > maxage, maxage, nages[i] )
    }


    minage <- Rceattle[[1]]$data_list$minage
    nsex <- Rceattle[[1]]$data_list$nsex

    # Get biomass
    m_array <-
      array(NA, dim = c(nspp, 2, max(nages), nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      m_array[, , , ,i] <- Rceattle[[i]]$quantities$M[,,1:maxage,1:nyrs]
    }

    if(log){
      m_array = log(m_array)
    }

    # Plot limits
    zmax <- c()
    zmin <- c()
    for (i in 1:dim(m_array)[1]) {
      zmax[i] <- max(c(m_array[i,,,,], 0), na.rm = T)
      zmin[i] <- min(c(m_array[i,,,,], 0), na.rm = T)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)

    #################################
    # Mortality time series
    #################################
    if(is.null(spp)){
      spp <- 1:nspp
    }


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
                width = ifelse(is.null(width), 8, width),
                height = ifelse(is.null(height), 5.5, height),

                units = "in",
                res = 300
              )
            }

            # Subset mortality data
            m_subset <- (m_array[j, sex, (1:nages[sp]), 1:nyrs, 1])

            # Get ages
            ages <- (1:(nages[sp])) - 1 + minage[sp]

            # Rearrange data
            data <- data.frame(Year = rep(Years[[1]], each = length(ages)), Age = rep(ages, length(Years[[1]])), M = c(m_subset))

            # Plot limits
            if(is.null(zlim)){
              zlim <- c(zmin[sp], zmax[sp])
            }

            if(is.character(zlim)){
              zlim <- c(min(zmin), max(zmax))
            }

            # Plot as contours
            if(contour){
              print(ggplot2::ggplot(data, aes(y = Age, x = Year, z = M, zmin = zlim[1], zmax = zlim[2])) + geom_contour(colour = 1, size = 0.5) + geom_contour_filled()  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) +  scale_x_continuous(expand = c(0, 0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size = 6)) + scale_fill_viridis_d("M1 + M2"))
            }

            # Plot as tiles
            if(contour == FALSE){
              p = ggplot2::ggplot(data, aes(y = Age, x = Year, zmin = zlim[1], zmax = zlim[2])) + geom_tile(aes(fill = M))  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) + coord_equal() +  scale_x_continuous(expand = c(0, 0))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size = 6))
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
           lwd = 3) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    maturity <- list()
    for(i in 1:length(Rceattle)){
      maturity[[i]] <- Rceattle[[i]]$data_list$pmature
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
          width = 4,# 169 / 25.4,
          height = 6.5, # 150 / 25.4,
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
        legend("topleft", species[j], bty = "n", cex = 0.8)

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
              cex = 0.8
            )
          }
        }
      }

      if (i == 2) {
        dev.off()
      }
    }
  }


#' Plot SSB
#'
#' @description Function the plots the spawning stock biomass trends as estimated from Rceattle. Returns and saves a figure with the ssb trajectory.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param spnames Species names for legend
#' @param lwd Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @export
#'
plot_ssb <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           spnames = NULL,
           lwd = 3,
           right_adj = 0,
           mohns = NULL,
           incl_proj = FALSE,
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
    Years <- list()
    Endyrs <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))
    nspp <- Rceattle[[1]]$data_list$nspp

    if(is.null(species)){
      species <- 1:nspp
    }


    # Get SSB
    SSB <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle)))
    ssb_sd <- array(NA, dim = c(nspp, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      SSB[, 1:length(Years[[i]]), i] <- Rceattle[[i]]$quantities$biomassSSB[,1:nyrs_vec[i]]

      # Get SD of rec
      if (add_ci) {
        ssb_sd_sub <- which(names(Rceattle[[i]]$sdrep$value) == "biomassSSB")
        ssb_sd_sub <- Rceattle[[i]]$sdrep$sd[ssb_sd_sub]
        ssb_sd[, , i] <-
          replace(ssb_sd[, , i], values = ssb_sd_sub[1:(nyrs_vec[i] * nspp)])
      }
    }

    SSB <- SSB / 1000000
    ssb_sd <- ssb_sd / 1000000
    SSB_upper <- SSB + ssb_sd * 1.92
    SSB_lower <- SSB - ssb_sd * 1.92

    # Plot limits
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(SSB)[1]) {
      ymax[i] <- max(c(SSB[i, , ], SSB_lower[i, , ], SSB_upper[i, , ], 0), na.rm = T)
      ymin[i] <- min(c(SSB[i, , ], SSB_lower[i, , ], SSB_upper[i, , ], 0), na.rm = T)
    }
    ymax <- ymax + 0.15 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_ssb_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

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

      for (j in species) {
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[j], ymax[j]),
          xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
          xlab = "Year",
          ylab = "SSB (million t)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", spnames[j], bty = "n", cex = 0.8)

        if(!is.null(mohns)){
          legend("top", paste0("SSB Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 0.8) # SSB rho
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
              cex = 0.8
            )
          }
        }



        # Mean SSB
        for (k in 1:dim(SSB)[3]) {
          lines(
            x = Years[[k]],
            y = SSB[j, 1:length(Years[[k]]), k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
        }

        # Credible interval
        if (add_ci) {
          for (k in 1:dim(SSB)[3]) {
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(SSB_upper[j, 1:length(Years[[k]]), k], rev(SSB_lower[j, 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = 0.4),
              border = NA
            )
          }
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
#' @param lwd Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @export
#'
plot_b_eaten <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           lwd = 3,
           right_adj = 0,
           mohns = NULL,
           incl_proj = FALSE,
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
    Years <- list()
    Endyrs <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }

    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp


    # Get B_eaten
    B_eaten <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle)))
    ssb_sd <- array(NA, dim = c(nspp, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      for(sp in 1:nspp){
        for(yr in 1:nyrs_vec[i]){
          B_eaten[sp, yr, i] <- sum(Rceattle[[i]]$quantities$B_eaten[sp,,,yr])


          # # Get SD of ssb
          # if (add_ci) {
          #   ssb_sd_sub <- which(names(Rceattle[[i]]$sdrep$value) == "biomassSSB")
          #   ssb_sd_sub <- Rceattle[[i]]$sdrep$sd[ssb_sd_sub]
          #   ssb_sd[, , i] <-
          #     replace(ssb_sd[, , i], values = ssb_sd_sub[1:(nyrs_vec[i] * nspp)])
          # }
        }
      }
    }




    ind = 1

    B_eaten <- B_eaten / 1000000
    ssb_sd <- ssb_sd / 1000000
    SSB_upper <- B_eaten + ssb_sd * 1.92
    SSB_lower <- B_eaten - ssb_sd * 1.92

    # Plot limits
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(B_eaten)[1]) {
      ymax[i] <- max(c(B_eaten[i, , ], SSB_lower[i, , ], SSB_upper[i, , ], 0), na.rm = T)
      ymin[i] <- min(c(B_eaten[i, , ], SSB_lower[i, , ], SSB_upper[i, , ], 0), na.rm = T)
    }
    ymax <- ymax + 0.15 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_b_eaten_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(nspp + 2), nrow = (nspp + 2)), heights = c(0.1, rep(1, nspp), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      for (j in 1:nspp) {
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[j], ymax[j]),
          xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
          xlab = "Year",
          ylab = "Biomass consumed (million t)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", spnames[j], bty = "n", cex = 0.8)

        if(!is.null(mohns)){
          legend("top", paste0("B_eaten Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 0.8) # B_eaten rho
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
              cex = 0.8
            )
          }
        }



        # Mean B_eaten
        for (k in 1:dim(B_eaten)[3]) {
          lines(
            x = Years[[k]],
            y = B_eaten[j, 1:length(Years[[k]]), k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
        }

        # Credible interval
        if (add_ci) {
          for (k in 1:dim(B_eaten)[3]) {
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(SSB_upper[j, 1:length(Years[[k]]), k], rev(SSB_lower[j, 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = 0.4),
              border = NA
            )
          }
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
#' @param lwd Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
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
           lwd = 3,
           right_adj = 0,
           mohns = NULL,
           incl_proj = FALSE,
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
    Years <- list()
    Endyrs <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }

    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))
    max_age <- max(Rceattle[[1]]$data_list$nages)

    nspp <- Rceattle[[1]]$data_list$nspp


    # Get B_eaten
    B_eaten_prop <-
      array(NA, dim = c(nspp * 2, nspp * 2, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      for(rsp in 1:(nspp)){
        for(ksp in 1:(nspp)){
          for(yr in 1:nyrs_vec[i]){
            B_eaten_prop[rsp, ksp,yr,i] <- sum(Rceattle[[i]]$quantities$B_eaten_prop[c(rsp, rsp + nspp),c(ksp, ksp + nspp),,,yr])
          }
        }
      }
    }

    ind = 1
    B_eaten_prop <- B_eaten_prop/1000000

    # Plot limits
    ymax <- c()
    ymin <- c()
    for (ksp in 1:dim(B_eaten_prop)[2]) { # Loop through prey
      ymax[ksp] <- max(c(B_eaten_prop[,ksp, , ], 0), na.rm = T)
      ymin[ksp] <- min(c(B_eaten_prop[,ksp, , ], 0), na.rm = T)
    }
    ymax <- ymax + 0.15 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_b_eaten_prop_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(nspp + 2), nrow = (nspp + 2)), heights = c(0.1, rep(1, nspp), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      for (j in 1:nspp) {
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[j], ymax[j]),
          xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj),
          xlab = "Year",
          ylab = "Biomass consumed (million t)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", spnames[j], bty = "n", cex = 0.8)

        if(!is.null(mohns)){
          legend("top", paste0("B_eaten_prop Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 0.8) # B_eaten_prop rho
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
              cex = 0.8
            )
          }
        }

        if (j == 2) {
          legend(
            "topright",
            legend = c("Predator:", spnames),
            lty = c(NA, 1:nspp),
            lwd = lwd,
            col = c(0, rep(1, nspp)),
            bty = "n",
            cex = 0.8
          )
        }



        # Mean B_eaten_prop
        for (mod in 1:dim(B_eaten_prop)[4]) {
          for(pred in 1:nspp){
            lines(
              x = Years[[k]],
              y = B_eaten_prop[pred, j, 1:length(Years[[k]]), mod],
              lwd = lwd,
              lty = pred,
              col = line_col[mod]) # Median
          }
        }
      }

      # # Credible interval
      # if (add_ci) {
      #   for (k in 1:dim(B_eaten)[3]) {
      #     polygon(
      #       x = c(Years[[k]], rev(Years[[k]])),
      #       y = c(SSB_upper[j, 1:length(Years[[k]]), k], rev(SSB_lower[j, 1:length(Years[[k]]), k])),
      #       col = adjustcolor( line_col[k], alpha.f = 0.4),
      #       border = NA
      #     )
      #   }
      # }



      if (i == 2) {
        dev.off()
      }
    }
  }
