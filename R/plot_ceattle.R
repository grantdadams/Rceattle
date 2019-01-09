#' plot_biomass
#'
#' @description Function the plots the biomass and spawning stock biomass trends as estimated from Rceattle
#'
#' @param file_name name of a file to identified the files exported by the
#'   function.
#' @param ceattle_list List of CEATTLE model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param include_srv Boolian of whether to include survey biomass estimates and 95% CI
#'
#' @return Returns and saves a figure with the population trajectory.
plot_biomass <-
  function(ceattle_list,
           tmp_list = NULL,
           file_name = NULL,
           model_names = NULL,
           line_col = NULL,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3,
           include_srv = FALSE) {
    library(extrafont)

    # Extract data objects
    nyrs <- ceattle_list[[1]]$data_list$nyrs
    Years <-
      ceattle_list[[1]]$data_list$styr:(ceattle_list[[1]]$data_list$styr + nyrs - 1)
    nspp <- ceattle_list[[1]]$data_list$nspp

    # Get biomass
    Biomass <-
      array(NA, dim = c(nspp, nyrs, length(ceattle_list) + length(tmp_list)))
    for (i in 1:length(ceattle_list)) {
      Biomass[, , i] <- ceattle_list[[i]]$quantities$biomass
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(ceattle_list) + 1):(length(ceattle_list) + length(tmp_list))) {
        for (k in 1:nspp) {
          Biomass[k, , i] <- tmp_list[[ind]][[paste0("Biomass_", k)]]
        }
        ind = ind + 1
      }
    }
    Biomass <- Biomass / 1000000

    # Get SSB
    SSB <-
      array(NA, dim = c(nspp, nyrs, length(ceattle_list) + length(tmp_list)))
    for (i in 1:length(ceattle_list)) {
      SSB[, , i] <- ceattle_list[[i]]$quantities$biomassSSB
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(ceattle_list) + 1):(length(ceattle_list) + length(tmp_list))) {
        for (k in 1:nspp) {
          SSB[k, , i] <- tmp_list[[ind]][[paste0("BiomassSSB_", k)]]
        }
        ind = ind + 1
      }
    }

    SSB <- SSB / 1000000

    # Plot limits
    ymax <- c()
    for (i in 1:dim(Biomass)[1]) {
      ymax[i] <- max(Biomass[i, , ], na.rm = T)
    }
    ymax <- ymax + 0.15 * ymax
    ymin <- 0

    if (is.null(line_col)) {
      line_col <- 1:length(ceattle_list)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file_name), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file_name, "_biomass_trajectory", ".tiff")
        tiff(
          file = filename ,
          width = 169 / 25.4,
          height = 150 / 25.4,
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
          ylim = c(ymin, ymax[j]),
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
          srv_yrs <- ceattle_list[[1]]$data_list$yrs_srv_biom[j, ]
          srv_biom <- ceattle_list[[1]]$data_list$srv_biom[j, ]
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
              ceattle_list[[1]]$data_list$srv_biom[j, ] + ceattle_list[[1]]$data_list$srv_biom_se[j, ] * 1.92
            ) / 1000000
          Lower95 <-
            (
              ceattle_list[[1]]$data_list$srv_biom[j, ] - ceattle_list[[1]]$data_list$srv_biom_se[j, ] * 1.92
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
#' @param ceattle_list List of CEATTLE model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param ci_col Colors to be used for CI color
#' @param lwd Line width as specified by user
#'
#' @return Returns and saves a figure with the population trajectory.
plot_recruitment <-
  function(ceattle_list,
           tmp_list = NULL,
           file_name = NULL,
           model_names = NULL,
           line_col = NULL,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           ci_col = NULL,
           lwd = 3,
           save_rec = FALSE) {
    library(extrafont)

    # Extract data objects
    nyrs <- ceattle_list[[1]]$data_list$nyrs
    Years <-
      ceattle_list[[1]]$data_list$styr:(ceattle_list[[1]]$data_list$styr + nyrs - 1)
    nspp <- ceattle_list[[1]]$data_list$nspp


    # Get biomass
    recruitment <-
      array(NA, dim = c(nspp, nyrs,  length(ceattle_list) + length(tmp_list)))
    recruitment_sd <-
      array(NA, dim = c(nspp, nyrs,  length(ceattle_list) + length(tmp_list)))
    for (i in 1:length(ceattle_list)) {
      recruitment[, , i] <- ceattle_list[[i]]$quantities$R[, ]

      # Get SD of rec
      sd_rec <- which(names(ceattle_list[[i]]$sdrep$value) == "R")
      sd_rec <- ceattle_list[[i]]$sdrep$sd[sd_rec]
      recruitment_sd[, , i] <-
        replace(recruitment_sd[, , i], values = sd_rec)
    }

    ind = 1
    if (!is.null(tmp_list)) {
      for (i in (length(ceattle_list) + 1):(length(ceattle_list) + length(tmp_list))) {
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
    for (i in 1:dim(recruitment)[1]) {
      ymax[i] <- max(recruitment_upper[i, , ], na.rm = T)
    }
    ymax <- ymax + 0.2 * ymax
    ymin <- 0

    if (is.null(line_col)) {
      line_col <- 1:length(ceattle_list)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file_name), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file_name, "_recruitment_trajectory", ".tiff")
        tiff(
          file = filename ,
          width = 169 / 25.4,
          height = 150 / 25.4,
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
          ylim = c(ymin, ymax[j]),
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
