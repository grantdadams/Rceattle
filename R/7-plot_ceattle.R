#' plot_biomass
#'
#' @description Function the plots the biomass and spawning stock biomass trends as estimated from Rceattle
#'
#' @param file_name name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param include_srv Boolian of whether to include survey biomass estimates and 95 CI
#'
#' @return Returns and saves a figure with the population trajectory.
#' @export
plot_biomass <-
  function(Rceattle,
           tmp_list = NULL,
           file_name = NULL,
           model_names = NULL,
           line_col = NULL,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3,
           include_srv = FALSE) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    nyrs <- Rceattle[[1]]$data_list$nyrs
    Years <-
      Rceattle[[1]]$data_list$styr:(Rceattle[[1]]$data_list$styr + nyrs - 1)
    nspp <- Rceattle[[1]]$data_list$nspp

    # Get biomass
    Biomass <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      Biomass[, , i] <- Rceattle[[i]]$quantities$biomass
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(Rceattle) + 1):(length(Rceattle) + length(tmp_list))) {
        for (k in 1:nspp) {
          Biomass[k, , i] <- tmp_list[[ind]][[paste0("Biomass_", k)]]
        }
        ind = ind + 1
      }
    }
    Biomass <- Biomass / 1000000

    # Get SSB
    SSB <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      SSB[, , i] <- Rceattle[[i]]$quantities$biomassSSB
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(Rceattle) + 1):(length(Rceattle) + length(tmp_list))) {
        for (k in 1:nspp) {
          SSB[k, , i] <- tmp_list[[ind]][[paste0("BiomassSSB_", k)]]
        }
        ind = ind + 1
      }
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
      line_col <- 1:length(Rceattle)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file_name), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file_name, "_biomass_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,
          family = "Helvetica" ,
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
          xlim = c(min(Years), max(Years)),
          xlab = "Year",
          ylab = "Biomass (million t)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Legends
        legend("topleft", species[j], bty = "n", cex = 1.4)

        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              lty = rep(1, length(line_col)),
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = 1.175
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
            cex = 1.175
          )
        }

        # Survey data
        if (include_srv) {
          srv_yrs <- Rceattle[[1]]$data_list$yrs_srv_biom[j, ]
          srv_biom <- Rceattle[[1]]$data_list$srv_biom[j, ]
          points(
            x = srv_yrs,
            y =  srv_biom / 1000000,
            pch = 16,
            cex = 1,
            col = "#696773"
          )


          # Get 95% CI and range
          Upper95 <-
            (
              Rceattle[[1]]$data_list$srv_biom[j, ] + Rceattle[[1]]$data_list$srv_biom_se[j, ] * 1.92
            ) / 1000000
          Lower95 <-
            (
              Rceattle[[1]]$data_list$srv_biom[j, ] - Rceattle[[1]]$data_list$srv_biom_se[j, ] * 1.92
            ) / 1000000


          arrows(
            x0 = srv_yrs,
            y0 = Upper95,
            x1 = srv_yrs,
            y1 = Lower95,
            length = 0.05,
            angle = 90,
            code = 3,
            lwd = 2,
            col = "#696773"
          )
        }


        # Mean biomass
        for (k in 1:dim(Biomass)[3]) {
          lines(
            x = Years,
            y = Biomass[j, , k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
          lines(
            x = Years,
            y = SSB[j, , k],
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
#' @param file_name name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param ci_col Colors to be used for CI color
#' @param lwd Line width as specified by user
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_recruitment <-
  function(Rceattle,
           tmp_list = NULL,
           file_name = NULL,
           model_names = NULL,
           line_col = NULL,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           ci_col = NULL,
           lwd = 3,
           save_rec = FALSE) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    nyrs <- Rceattle[[1]]$data_list$nyrs
    Years <-
      Rceattle[[1]]$data_list$styr:(Rceattle[[1]]$data_list$styr + nyrs - 1)
    nspp <- Rceattle[[1]]$data_list$nspp


    # Get biomass
    recruitment <-
      array(NA, dim = c(nspp, nyrs,  length(Rceattle) + length(tmp_list)))
    recruitment_sd <-
      array(NA, dim = c(nspp, nyrs,  length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      recruitment[, , i] <- Rceattle[[i]]$quantities$R[, ]

      # Get SD of rec
      if (!is.null(ci_col)) {
        sd_rec <- which(names(Rceattle[[i]]$sdrep$value) == "R")
        sd_rec <- Rceattle[[i]]$sdrep$sd[sd_rec]
        recruitment_sd[, , i] <-
          replace(recruitment_sd[, , i], values = sd_rec)
      }
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(Rceattle) + 1):(length(Rceattle) + length(tmp_list))) {
        for (k in 1:nspp) {
          recruitment[k, , i] <- tmp_list[[ind]][[paste0("R_", k)]]
          recruitment_sd[k, , i] <-
            replace(recruitment_sd[k, , i], values = rep(NA, length(recruitment_sd[k, , i])))
        }
        ind = ind + 1
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
          paste0(file_name, "_recruitment_species_", i, ".csv")
        write.csv(dat_new, file = filename)
      }
    }


    # Plot limits
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(recruitment)[1]) {
      if (!is.null(ci_col)) {
        ymax[i] <- max(c(recruitment_upper[i, , ], 0), na.rm = T)
        ymin[i] <- min(c(recruitment_upper[i, , ], 0), na.rm = T)
      } else{
        ymax[i] <- max(c(recruitment[i, , ], 0), na.rm = T)
        ymin[i] <- min(c(recruitment[i, , ], 0), na.rm = T)
      }
    }
    ymax <- ymax + 0.2 * ymax

    if (is.null(line_col)) {
      line_col <- 1:length(Rceattle)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file_name), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file_name, "_recruitment_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,
          family = "Helvetica",
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
          xlim = c(min(Years), max(Years)),
          xlab = "Year",
          ylab = "Recruitment (millions)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Legends
        legend("topleft",
               legend = species[j],
               bty = "n",
               cex = 1.4)
        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              lty = rep(1, length(line_col)),
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = 1.175
            )
          }

        }


        # Credible interval
        if (!is.null(ci_col)) {
          for (k in 1:dim(recruitment)[3]) {
            polygon(
              x = c(Years, rev(Years)),
              y = c(recruitment_upper[j, , k], rev(recruitment_lower[j, , k])),
              col = adjustcolor( ci_col[k], alpha.f = 0.2),
              border = NA
            ) # 95% CI
          }
        }

        # Mean recruitment
        for (k in 1:dim(recruitment)[3]) {
          lines(
            x = Years,
            y = recruitment[j, , k],
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
#' @param file_name name of a file to identified the files exported by the
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
           tmp_list = NULL,
           file_name = NULL,
           model_names = NULL,
           line_col = NULL,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    nyrs <- Rceattle[[1]]$data_list$nyrs
    Years <-
      Rceattle[[1]]$data_list$styr:(Rceattle[[1]]$data_list$styr + nyrs - 1)
    nspp <- Rceattle[[1]]$data_list$nspp
    nages <- Rceattle[[1]]$data_list$nages

    # Get biomass
    srv_selectivity <-
      array(NA, dim = c(nspp, max(nages), length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      srv_selectivity[, , i] <- Rceattle[[i]]$quantities$srv_sel
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(Rceattle) + 1):(length(Rceattle) + length(tmp_list))) {
        for (k in 1:nspp) {
          srv_selectivity[k,1:nages[k], i] <- tmp_list[[ind]][[paste0("srv_sel_", k)]]
        }
        ind = ind + 1
      }
    }

    # Get SSB
    fsh_selectivity <-
      array(NA, dim = c(nspp, max(nages), length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      fsh_selectivity[, , i] <- Rceattle[[i]]$quantities$fsh_sel
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(Rceattle) + 1):(length(Rceattle) + length(tmp_list))) {
        for (k in 1:nspp) {
          fsh_selectivity[k, , i] <- tmp_list[[ind]][[paste0("fsh_sel_", k)]]
        }
        ind = ind + 1
      }
    }

    # Plot limits
    ymax_srv <- c()
    ymin_srv <- c()

    ymax_fsh <- c()
    ymin_fsh <- c()
    for (i in 1:dim(srv_selectivity)[1]) {
      ymax_srv[i] <- max(c(srv_selectivity[i,1:nages[i], ], 0), na.rm = T)
      ymin_srv[i] <- min(c(srv_selectivity[i,1:nages[i], ], 0), na.rm = T)

      ymax_fsh[i] <- max(c(fsh_selectivity[i,1:nages[i], ], 0), na.rm = T)
      ymin_fsh[i] <- min(c(fsh_selectivity[i,1:nages[i], ], 0), na.rm = T)
    }
    ymax_srv <- ymax_srv + 0.15 * ymax_srv
    ymax_fsh <- ymax_fsh + 0.15 * ymax_fsh

    if (is.null(line_col)) {
      line_col <- 1:length(Rceattle)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file_name), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file_name, "_selectivity", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5, # 150 / 25.4,
          family = "Helvetica" ,
          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:((nspp + 2)*2), nrow = (nspp + 2), ncol = 2, byrow = TRUE), heights = c(0.1, rep(1, nspp), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()
      plot.new()

      for (j in 1:nspp) {

        # Survey selectivity
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin_srv[j], ymax_srv[j]),
          xlim = c(min(0), nages[j]),
          xlab = "Age",
          ylab = "Survey selectivity",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Legends
        legend("topleft", species[j], bty = "n", cex = 1.4)

        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "bottomright",
              legend = model_names,
              lty = rep(1, length(line_col)),
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = 1.175
            )
          }
        }


        # Mean selectivity
        for (k in 1:dim(srv_selectivity)[3]) {
          lines(
            x = 1:nages[j],
            y = srv_selectivity[j, 1:nages[j], k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          )
        }


        # Fishery selectivity
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin_fsh[j], ymax_fsh[j]),
          xlim = c(min(0), nages[j]),
          xlab = "Age",
          ylab = "Fishery selectivity",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "bottomright",
              legend = model_names,
              lty = rep(1, length(line_col)),
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = 1.175
            )
          }
        }


        # Mean selectivity
        for (k in 1:dim(fsh_selectivity)[3]) {
          lines(
            x = 1:nages[j],
            y = fsh_selectivity[j, 1:nages[j], k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          )
        }
      }


      if (i == 2) {
        dev.off()
      }
    }
  }


#' Plot functional form
#'
#' @describtion Function to plot the functional form estimated or specified by \code{\link{Rceattle}}
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


#' plot_M2
#'
#' @description Function the plots the predation mortality trends as estimated from Rceattle
#'
#' @param file_name name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param age Age specified
#'
#' @return Returns and saves a figure with the population trajectory.
#' @export
plot_mort <-
  function(Rceattle,
           tmp_list = NULL,
           file_name = NULL,
           model_names = NULL,
           line_col = NULL,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3,
           age = 3) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    nyrs <- Rceattle[[1]]$data_list$nyrs
    Years <-
      Rceattle[[1]]$data_list$styr:(Rceattle[[1]]$data_list$styr + nyrs - 1)
    nspp <- Rceattle[[1]]$data_list$nspp
    max_age <- max(Rceattle[[1]]$data_list$nages)

    # Get M2
    M2 <-
      array(NA, dim = c(nspp, max_age, nyrs, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      M2[, , ,i] <- Rceattle[[i]]$quantities$M2
    }

    # ind = 1
    # if (!is.null(tmp_list)) {
    #   for (i in (length(Rceattle) + 1):(length(Rceattle) + length(tmp_list))) {
    #     for (k in 1:nspp) {
    #       M2[k, , i] <- tmp_list[[ind]][[paste0("M2_", k)]]
    #     }
    #     ind = ind + 1
    #   }
    # }

    # Get f_mat
    f_mat <-
      array(NA, dim = c(nspp, max_age, nyrs, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      f_mat[, , ,i] <- Rceattle[[i]]$quantities$F
    }

    m_mat <-
      array(NA, dim = c(nspp, max_age, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      m_mat[, ,i] <- Rceattle[[i]]$quantities$M1
    }

    # ind = 1
    # if (!is.null(tmp_list)) {
    #   for (i in (length(Rceattle) + 1):(length(Rceattle) + length(tmp_list))) {
    #     for (k in 1:nspp) {
    #       f_mat[k, , i] <- tmp_list[[ind]][[paste0("F_", k)]]
    #     }
    #     ind = ind + 1
    #   }
    # }

    # Plot limits
    #if(fishing){
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(M2)[1]) {
      ymax[i] <- max(c(M2[i, age, ,], f_mat[i, age, ,], m_mat[i, age,], 0), na.rm = T)
      ymin[i] <- min(c(M2[i, age, ,], f_mat[i, age, ,], m_mat[i, age,], 0), na.rm = T)
    }
    ymax <- ymax + 0.15 * ymax
    # }
    # if(fishing == FALSE){
    #   ymax <- c()
    #   ymin <- c()
    #   for (i in 1:dim(M2)[1]) {
    #     ymax[i] <- max(c(M2[i, , ,], 0), na.rm = T)
    #     ymin[i] <- min(c(M2[i, , ,], 0), na.rm = T)
    #   }
    #   ymax <- ymax + 0.15 * ymax
    # }

    if (is.null(line_col)) {
      line_col <- 1:length(Rceattle)
    }

    alphas <- seq(1, 0.1, length.out = nages)

    # Plot trajectory
    loops <- ifelse(is.null(file_name), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file_name,"_age",age, "_mortality_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,
          family = "Helvetica" ,
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
          xlim = c(min(Years), max(Years)),
          xlab = "Year",
          ylab = "Mortality",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Legends
        legend("topleft", species[j], bty = "n", cex = 1.4)

        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              lty = rep(1, length(line_col)),
              lwd = lwd,
              col = line_col,
              bty = "n",
              cex = 1.175
            )
          }
        }

        if (j == 2) {
          legend(
            "topright",
            legend = c("M1", "M2", "F"),
            lty = c(3, 1, 2),
            lwd = lwd,
            col = c(1, 1, 1),
            bty = "n",
            cex = 1.175
          )
        }



        # Mean M2
        for (k in 1:dim(M2)[4]) {
          #if(fishing == FALSE){
          lines(
            x = Years,
            y = M2[j, age, , k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
          #}
          #if(fishing){
          lines(
            x = Years,
            y = f_mat[j, age, , k],
            lty = 2,
            lwd = lwd,
            col = line_col[k]
          ) # Median

          # M
          abline(
            h = m_mat[j, age, k],
            lty = 3,
            lwd = lwd,
            col = line_col[k]
          ) # Median
          #}
        }
      }


      if (i == 2) {
        dev.off()
      }
    }
  }


