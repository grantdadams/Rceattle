#' Plot time series of comp data
#'
#' @description Function the plots the comp data as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param species Species names for legend
#' @param cex Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#'
#' @return Returns and saves a figure
#' @export
plot_comp <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           species = NULL,
           cex = 3,
           lwd = 3,
           right_adj = 0) {

    # Make sure we are using only one model
    if(class(Rceattle) != "Rceattle"){
      stop("Please only use one Rceattle model")
    }

    # Species names
    if(is.null(species)){
      species =  Rceattle$data_list$spnames
    }


    # Get colors
    colvec=c("red", "blue", "black")

    # Extract data objects
    # Get observed
    comp_data <- Rceattle$data_list$comp_data
    # Normalize
    comp_data[,grep("Comp_", colnames(comp_data))] = comp_data[,grep("Comp_", colnames(comp_data))]/rowSums(comp_data[,grep("Comp_", colnames(comp_data))], na.rm = TRUE)

    fleet_control <- Rceattle$data_list$fleet_control


    # Get estimated
    comp_hat <- Rceattle$data_list$comp_data
    comp_hat[,grep("Comp_", colnames(comp_data))] = Rceattle$quantities$comp_hat


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # # Plot comps Type 1 ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # for(comp_type in c(0,1)){
    #
    # TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
    ## devtools::install_github("fishfollower/compResidual/compResidual")
    # o <- round(Neff*obs/rowSums(obs),0); p=exp/rowSums(exp)
    # ## default output
    # res <-  compResidual::resMulti(t(o), t(p))

    #
    #   srv <- unique(comp_data$Fleet_code[which(comp_data$Age0_Length1 == comp_type)])
    #   nsrv <- length(srv)
    #
    #   loops <- ifelse(is.null(file), 1, 2)
    #   for (i in 1:loops) {
    #     if (i == 2) {
    #       filename <- paste0(file, c("_survey_age_comps_1", "_survey_length_comps_1")[comp_type + 1], ".png")
    #       png(
    #         file = filename ,
    #         width = 7,# 169 / 25.4,
    #         height = 6.5,# 150 / 25.4,
    #
    #         units = "in",
    #         res = 300
    #       )
    #     }
    #
    #     # Plot configuration
    #     layout(matrix(1:(nsrv + 2), nrow = (nsrv + 2)), heights = c(0.1, rep(1, nsrv), 0.2))
    #     par(
    #       mar = c(0, 3 , 0 , 1) ,
    #       oma = c(0 , 0 , 0 , 0),
    #       tcl = -0.35,
    #       mgp = c(1.75, 0.5, 0)
    #     )
    #     plot.new()
    #
    #     for (j in 1:nsrv) {
    #
    #
    #       # Extract comps
    #       srv_ind <- which(comp_data$Fleet_code == srv[j] & comp_data$Age0_Length1 == comp_type)
    #       comp_tmp <- comp_data[srv_ind,]
    #       comp_hat_tmp <- comp_hat[srv_ind, ]
    #
    #       # Reorganize and clean
    #       comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
    #       comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))
    #       comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]
    #       comp_tmp <- comp_tmp[which(comp_tmp$comp > 0),]
    #
    #       comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
    #       comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))
    #       comp_hat_tmp <- comp_hat_tmp[which(!is.na(comp_hat_tmp$comp)),]
    #       comp_hat_tmp <- comp_hat_tmp[which(comp_hat_tmp$comp > 0),]
    #
    #       sp <- fleet_control$Species[which(fleet_control$Fleet_code == srv[j])]
    #       nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]
    #
    #       plot(
    #         y = NA,
    #         x = NA,
    #         ylim = c(0, nages[sp] * 1.20),
    #         xlim = c(min(comp_data$Year, na.rm = TRUE), max(comp_data$Year, na.rm = TRUE) + right_adj),
    #         xlab = "Year",
    #         ylab = c("Survey age comp", "Survey length comp")[comp_type + 1],
    #         xaxt = c(rep("n", nsrv - 1), "s")[j]
    #       )
    #
    #
    #       # Horizontal line
    #       if(incl_proj){
    #         abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
    #       }
    #
    #       # Legends
    #       legend("topleft", as.character(fleet_control$Fleet_name[srv[j]]), bty = "n", cex = 1.4)
    #
    #
    #       # Type
    #       if (j == 1) {
    #         legend(
    #           "topright",
    #           legend = c("Observed", "Estimated"),
    #           pch = c(16, 16),
    #           cex = 1.125,
    #           col = line_col,
    #           bty = "n"
    #         )
    #
    #         x_loc <- c(mean(comp_data$Year, na.rm = T) - 2, mean(comp_data$Year, na.rm = T), mean(comp_data$Year, na.rm = T) +2)
    #         symbols( x = x_loc , y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.07, 3) , circle = c(0.1, 0.25, 0.5), inches=0.10,add=T, fg = line_col[2])
    #         text(x = x_loc, y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.16, 3), labels = c(0.1, 0.25, 0.5))
    #       }
    #
    #
    #       if(nrow(comp_tmp) > 0){
    #         # Observed
    #         symbols( x = comp_tmp$Year , y = comp_tmp$age , circle = comp_tmp$comp, inches=0.10,add=T, fg = line_col[1])
    #
    #         # Estimated
    #         symbols( x = comp_hat_tmp$Year , y = comp_hat_tmp$age , circle = comp_hat_tmp$comp, inches=0.10,add=T, fg = line_col[2])
    #       }
    #     }
    #
    #
    #     if (i == 2) {
    #       dev.off()
    #     }
    #   }
    # }
    #
    #

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Plot comps Type 2 - Pearson residual ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    for(comp_type in c(0,1)){ # Age0, Length 1

      srv <- unique(comp_data$Fleet_code[which(comp_data$Age0_Length1 == comp_type)])
      nsrv <- length(srv)


      if(nsrv > 0){
        loops <- ifelse(is.null(file), 1, 2)
        for (i in 1:loops) {

          # Plot configuration
          par(
            mar = c(3, 3 , 1 , 1) ,
            oma = c(0 , 0 , 0 , 0),
            tcl = -0.35,
            mgp = c(1.75, 0.5, 0)
          )

          # Loop around fleets with comp data
          for (j in 1:nsrv) {

            #  Plot name
            if (i == 2) {
              filename <- paste0(file, paste0(c("_comps_pearson_residual_", "_comps_pearson_residual_"), "fleet_code_",srv[j] )[comp_type + 1], ".png")
              png(
                file = filename ,
                width = 7,# 169 / 25.4,
                height = 6.5,# 150 / 25.4,

                units = "in",
                res = 300
              )
            }


            # Extract comps
            srv_ind <- which(comp_data$Fleet_code == srv[j] & comp_data$Age0_Length1 == comp_type)
            comp_tmp <- comp_data[srv_ind,]
            comp_hat_tmp <- comp_hat[srv_ind, ]

            # Reorganize and clean
            comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
            comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))

            comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
            comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))

            # Calculate pearson residual
            comp_tmp$comp_hat <- comp_hat_tmp$comp
            comp_tmp <- comp_tmp[which(comp_tmp$comp > 0),]
            comp_tmp <- comp_tmp[which(comp_tmp$comp_hat > 0),]

            # Calculate pearson
            comp_tmp$pearson <- (comp_tmp$comp - comp_tmp$comp_hat) / sqrt( ( comp_tmp$comp_hat * (1 - comp_tmp$comp_hat)) / comp_tmp$Sample_size)
            max_pearson <- max(abs(comp_tmp$pearson), na.rm = TRUE)


            sp <- fleet_control$Species[which(fleet_control$Fleet_code == srv[j])]
            nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

            plot(
              y = NA,
              x = NA,
              ylim = c(0, nages[sp] * 1.25),
              xlim = c(min(abs(comp_tmp$Year), na.rm = TRUE), max(abs(comp_tmp$Year), na.rm = TRUE) + right_adj),
              xlab = "Year",
              ylab = c("Pearson residual: Age comp", "Pearson residual: Length comp")[comp_type + 1],
              xaxt = "s"
            )


            # Legends
            legend("topleft", as.character(fleet_control$Fleet_name[srv[j]]), bty = "n", cex = 1)



            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            # Legend ----
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            round = 0
            if(sum(seq(from = 1, to = max_pearson, length.out = 3) < 3) == 3){
              round = 1
            }

            # Positive
            x_loc <- c(max(abs(comp_tmp$Year), na.rm = T) - 6, max(abs(comp_tmp$Year), na.rm = T) - 3.5, max(abs(comp_tmp$Year), na.rm = T) - 1)
            symbols( x = x_loc , y = rep(nages[sp] * 1.1, 3) , circle = round(seq(from = 1, to = max_pearson, length.out = 3) , round), inches=0.20,add=T, bg = colvec[3])
            text(x = x_loc, y = rep(nages[sp] * 1.23, 3), labels = round(seq(from = 1, to = max_pearson, length.out = 3) , round))

            # Negative
            x_loc <- c(max(abs(comp_tmp$Year), na.rm = T) - 8.5, max(abs(comp_tmp$Year), na.rm = T) - 11, max(abs(comp_tmp$Year), na.rm = T) - 13.5)
            symbols( x = x_loc , y = rep(nages[sp] * 1.1, 3) , circle = -round(seq(from = -1, to = -max_pearson, length.out = 3) , round), inches=0.20,add=T, bg = NA)
            text(x = x_loc, y = rep(nages[sp] * 1.23, 3), labels = round(seq(from = -1, to = -max_pearson, length.out = 3) , round) )


            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            # Plot comp pearson residuals
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            if(nrow(comp_tmp) > 0){

              # Get colors
              comp_tmp$colors <- NA
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 0, colvec[3], comp_tmp$colors) # Combined sex / 1 sex model
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 1, colvec[1], comp_tmp$colors) # Females
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 2, colvec[2], comp_tmp$colors) # Males

              # Joint sex
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age <= nages[sp], colvec[1], comp_tmp$colors) # Females
              comp_tmp$colors <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age > nages[sp], colvec[2], comp_tmp$colors) # Males

              # Adjust years for joint
              comp_tmp$Year <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age > nages[sp], comp_tmp$Year + 0.2, comp_tmp$Year) # Males
              comp_tmp$Year <- ifelse(comp_tmp$Sex == 3 & comp_tmp$age <= nages[sp], comp_tmp$Year - 0.2, comp_tmp$Year) # Females

              # Adjust age for joint data
              comp_tmp$age <- ifelse(comp_tmp$age > nages[sp], comp_tmp$age - nages[sp], comp_tmp$age)

              # Background colors
              comp_tmp$bg_colors <- ifelse(comp_tmp$pearson > 0, comp_tmp$colors, NA)

              # Plot
              symbols( x = comp_tmp$Year , y = comp_tmp$age , circle = abs(comp_tmp$pearson), inches=0.2,add=T, bg = comp_tmp$bg_colors, fg = comp_tmp$colors)
            }



            if (i == 2) {
              dev.off()
            }
          }
        }
      }
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Plot comps Type 3 - Histograms ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    for(comp_type in c(0,1)){

      srv <- unique(comp_data$Fleet_code[which(comp_data$Age0_Length1 == comp_type)])
      nsrv <- length(srv)

      if(nsrv > 0){

        loops <- ifelse(is.null(file), 1, 2)
        for (i in 1:loops) {
          # Loop around fleets with comp data
          for (j in 1:nsrv) {

            #  Plot name
            if (i == 2) {
              filename <- paste0(file, paste0(c("_age_comps_histograms_", "_length_comps_histograms_"), "fleet_code_",srv[j] )[comp_type + 1], ".png")
              png(
                file = filename ,
                width = 7,# 169 / 25.4,
                height = 6.5,# 150 / 25.4,

                units = "in",
                res = 300
              )
            }

            # Extract comps
            srv_ind <- which(comp_data$Fleet_code == srv[j] & comp_data$Age0_Length1 == comp_type)
            comp_tmp <- comp_data[srv_ind,]
            comp_hat_tmp <- comp_hat[srv_ind, ]

            # Reorganize and clean
            comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
            comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))

            comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
            comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))

            # Combine
            comp_tmp$comp_hat <- comp_hat_tmp$comp
            comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]

            # Get comp dims
            sp <- fleet_control$Species[which(fleet_control$Fleet_code == srv[j])]
            nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

            yrs <- sort(unique(c(abs(comp_hat_tmp$Year), abs(comp_tmp$Year))))
            nyrs <- length(yrs)

            # Min and max
            max_comp <- max(c(comp_tmp$comp, comp_tmp$comp_hat))
            min_comp <- min(c(comp_tmp$comp, comp_tmp$comp_hat))

            # Plot configuration
            plot_rows <- ceiling(nyrs/4) + 1
            layout(matrix(1:(plot_rows * 5), ncol = 5, byrow = FALSE), widths = c(0.25,1,1,1,1))
            par(
              mar = c(0, 0 , 0 , 0) ,
              oma = c(0 , 0 , 0 , 0),
              tcl = -0.35,
              mgp = c(1.75, 0.5, 0)
            )

            # Plot black row
            for(k in 1:plot_rows){
              plot.new()
            }


            row_cnt = 0
            # Plot each comp for each year
            for(yr in 1:nyrs){
              row_cnt = row_cnt + 1

              # Subset year for observed and predicted comp
              comp_tmp_yr <- comp_tmp[which(abs(comp_tmp$Year) == yrs[yr]),]

              # Get lower bound for joint sex
              lowercomp <- 0
              sex <- unique(comp_tmp_yr$Sex)

              if(sex == 3){
                lowercomp <- - max_comp * 1.10
              }

              # Set plot up
              plot(
                y = NA,
                x = NA,
                ylim = c(lowercomp, max_comp * 1.10),
                xlim = c(0, nages[sp]),
                xlab = NA,
                ylab = NA,
                xaxt = "n",
                yaxt = "n"
              )

              # x-axis
              if(row_cnt == (plot_rows - 1) | yr == nyrs){
                axis(side = 1)
                mtext(text =  paste(as.character(fleet_control$Fleet_name[srv[j]]), c("age", "length"))[comp_type + 1],
                      side = 1, line = 2, cex = 0.7)
              }


              # y-axis
              if(yr < plot_rows){
                # Single sex or combined
                if(sex != 3){
                  y_axis <- round(seq(0, max_comp * 1.15, length.out = 4)[1:3],1)
                  axis(side = 2, at = y_axis, labels = c(0,y_axis[2:3]))
                }

                # Joint sex comp
                if(sex == 3){
                  y_axis <- round(seq(-max_comp * 1.15, max_comp * 1.15, length.out = 5)[1:5],1)
                  axis(side = 2, at = y_axis, labels = y_axis)
                }

                mtext(text =  paste0(c("Age comp", "Length comp"))[comp_type + 1],
                      side = 2, line = 2, cex = 0.7)
              }

              # Legends
              legend("topleft", legend = yrs[yr], bty = "n", cex = 1.2)


              # Plot observed and predicted comp
              # Single sex or combined
              if(sex < 3){
                if(sex == 1){
                  line_col <- colvec[1]
                }
                if(sex == 2){
                  line_col <- colvec[2]
                }
                if(sex == 0){
                  line_col <- colvec[3]
                }
                polygon(c(0,comp_tmp_yr$age, max(comp_tmp_yr$age) + 1), c(0, comp_tmp_yr$comp, 0),col='grey80',border=NA)
                lines(c(0,comp_tmp_yr$age, max(comp_tmp_yr$age) + 1), c(0, comp_tmp_yr$comp_hat, 0),col = line_col, lwd = lwd)
              }

              # Joint sex
              if(sex == 3){
                # Subset males and females and adjust ages
                comp_tmp_yr_males <- comp_tmp_yr[which(comp_tmp_yr$age > nages[sp]),]
                comp_tmp_yr_males$age <- comp_tmp_yr_males$age - nages[sp]
                comp_tmp_yr_females <- comp_tmp_yr[which(comp_tmp_yr$age <= nages[sp]),]

                # Plot females
                polygon(x = c(0,comp_tmp_yr_females$age, nages[sp] + 1), y = c(0, comp_tmp_yr_females$comp, 0),col='grey80',border=NA)
                lines(x = c(0,comp_tmp_yr_females$age, nages[sp] + 1), y =  c(0, comp_tmp_yr_females$comp_hat, 0),col = colvec[1], lwd = lwd, lty = ifelse(comp_tmp_yr$Year > 0, 1, 2))

                # Plot males
                polygon(x = c(0,comp_tmp_yr_males$age, nages[sp] + 1), y = c(0, -comp_tmp_yr_males$comp, 0),col='grey80',border=NA)
                lines(x = c(0,comp_tmp_yr_males$age, nages[sp] + 1), y = c(0, -comp_tmp_yr_males$comp_hat, 0),col = colvec[2], lwd = lwd, lty = ifelse(comp_tmp_yr$Year > 0, 1, 2))

                # Middle line
                abline( h = 0, col = 1)
              }


              # Make bottom row of empty plots
              if(row_cnt == (plot_rows - 1)){
                plot.new()
                row_cnt = 0
              }
            }

            if (i == 2) {
              dev.off()
            }
          }
        }
      }
    }



    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Plot comps Type 4 - Histograms of aggregated comps ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    for(comp_type in c(0,1)){

      srv <- unique(comp_data$Fleet_code[which(comp_data$Age0_Length1 == comp_type)])
      nsrv <- length(srv)

      if(nsrv > 0){

        loops <- ifelse(is.null(file), 1, 2)
        for (i in 1:loops) {
          # Loop around fleets with comp data
          for (j in 1:nsrv) {

            #  Plot name
            if (i == 2) {
              filename <- paste0(file, paste0(c("_aggregated_age_comps_histograms_", "_aggregated_length_comps_histograms_"), "fleet_code_",srv[j] )[comp_type + 1], ".png")
              png(
                file = filename ,
                width = 7,# 169 / 25.4,
                height = 6.5,# 150 / 25.4,

                units = "in",
                res = 300
              )
            }

            par(
              mfrow = c(1,1),
              mar = c(3.5, 3.5 , 0.5 , 0.5) ,
              oma = c(0 , 0 , 0 , 0),
              tcl = -0.35,
              mgp = c(1.75, 0.5, 0)
            )

            # Extract comps
            srv_ind <- which(comp_data$Fleet_code == srv[j] & comp_data$Age0_Length1 == comp_type)
            comp_tmp <- comp_data[srv_ind,]
            comp_hat_tmp <- comp_hat[srv_ind, ]

            # Reorganize and clean
            comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
            comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))
            comp_tmp <- comp_tmp[which(comp_tmp$Year > 0),]

            comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
            comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))
            comp_hat_tmp <- comp_hat_tmp[which(comp_hat_tmp$Year > 0),]

            # Combine
            comp_tmp$comp_hat <- comp_hat_tmp$comp
            comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]

            # Combine across year for observed and predicted comp
            comp_tmp_sum <- aggregate(comp_tmp$comp, by=list(Category=comp_tmp$age), FUN=sum)
            comp_tmp_sum$x <- comp_tmp_sum$x / sum(comp_tmp_sum$x)

            comp_hat_tmp_sum <- aggregate(comp_tmp$comp_hat, by=list(Category=comp_tmp$age), FUN=sum)
            comp_hat_tmp_sum$x <- comp_hat_tmp_sum$x / sum(comp_hat_tmp_sum$x)


            # Get comp dims
            sp <- fleet_control$Species[which(fleet_control$Fleet_code == srv[j])]
            nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

            yrs <- sort(unique(c(comp_hat_tmp$Year, comp_tmp$Year)))
            nyrs <- length(yrs)


            # Min and max
            max_comp <- max(c(comp_hat_tmp_sum$x, comp_tmp_sum$x))
            min_comp <- min(c(comp_hat_tmp_sum$x, comp_tmp_sum$x))

            sex <- unique(comp_tmp$Sex)
            lowercomp <- 0
            if(sex == 3){
              lowercomp <- - max_comp * 1.10
            }

            # Set plot up
            plot(
              y = NA,
              x = NA,
              ylim = c(lowercomp, max_comp * 1.10),
              xlim = c(0, nages[sp]),
              xlab = NA,
              ylab = NA,
              xaxt = "n",
              yaxt = "n"
            )

            # x-axisyrs){
            axis(side = 1)
            mtext(text =  paste(as.character(fleet_control$Fleet_name[srv[j]]), c("age", "length"))[comp_type + 1],
                  side = 1, line = 2, cex = 0.7)



            # y-axis
            # Single sex or combined
            if(sex != 3){
              y_axis <- round(seq(0, max_comp * 1.15, length.out = 4)[1:3],2)
              axis(side = 2, at = y_axis, labels = c(0,y_axis[2:3]))
            }

            # Joint sex comp
            if(sex == 3){
              y_axis <- round(seq(-max_comp * 1.15, max_comp * 1.15, length.out = 5)[1:5],2)
              axis(side = 2, at = y_axis, labels = y_axis)
            }

            mtext(text =  paste0(c("Aggregated age comp", "Aggregated length comp"))[comp_type + 1],
                  side = 2, line = 2, cex = 0.7)



            # Plot observed and predicted comp
            # Single sex or combined
            if(sex < 3){
              if(sex == 1){
                line_col <- colvec[1]
              }
              if(sex == 2){
                line_col <- colvec[2]
              }
              if(sex == 0){
                line_col <- colvec[3]
              }
              polygon(c(0,comp_tmp_sum$Category, nages[sp] + 1), c(0, comp_tmp_sum$x, 0),col='grey80',border=NA)
              lines(c(0,comp_hat_tmp_sum$Category, nages[sp] + 1), c(0, comp_hat_tmp_sum$x, 0),col = line_col, lwd = lwd)
            }

            # Joint sex
            if(sex == 3){
              # Subset males and females and adjust ages
              comp_tmp_sum_males <- comp_tmp_sum[which(comp_tmp_sum$Category > nages[sp]),]
              comp_tmp_sum_males$Category <- comp_tmp_sum_males$Category - nages[sp]
              comp_tmp_sum_females <- comp_tmp_sum[which(comp_tmp_sum$Category <= nages[sp]),]

              comp_hat_tmp_sum_males <- comp_hat_tmp_sum[which(comp_hat_tmp_sum$Category > nages[sp]),]
              comp_hat_tmp_sum_males$Category <- comp_hat_tmp_sum_males$Category - nages[sp]
              comp_hat_tmp_sum_females <- comp_hat_tmp_sum[which(comp_hat_tmp_sum$Category <= nages[sp]),]

              # Plot females
              polygon(x = c(0,comp_tmp_sum_females$Category, nages[sp] + 1), y = c(0, comp_tmp_sum_females$x, 0),col='grey80',border=NA)
              lines(x = c(0,comp_hat_tmp_sum_females$Category, nages[sp] + 1), y =  c(0, comp_hat_tmp_sum_females$x, 0),col = colvec[1], lwd = lwd)

              # Plot males
              polygon(x = c(0,comp_tmp_sum_males$Category, nages[sp] + 1), y = c(0, -comp_tmp_sum_males$x, 0),col='grey80',border=NA)
              lines(x = c(0,comp_hat_tmp_sum_males$Category, nages[sp] + 1), y = c(0, -comp_hat_tmp_sum_males$x, 0),col = colvec[2], lwd = lwd)

              # Middle line
              abline( h = 0, col = 1)
            }



            if (i == 2) {
              dev.off()
            }
          }
        }
      }
    }
  }







#' Plot diet composition fits
#'
#' @description
#' If year == 0, diet data are averaged from suit_styr to suit_endyr
#' If prey_age >= 0 diet data are diet proportion of prey-at-age in predator-at-age
#' If prey_age < 0 diet data are diet proportion of prey-spp in predator-at-age (sum across prey ages)
#' If prey_age < 0 and pred_age < 0, diet data are mean diet proportion of prey-spp in predator-spp (sum across prey ages and take mean across predator ages)
#' If prey_age < 0 and pred_age < -500, diet data are weighted mean diet proportion of prey-spp in predator-spp (sum across prey ages and take weighted mean across predator ages)
#'
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param species Species names for legend
#'
#' @return Returns and saves a figure
#' @export
plot_diet_comp <-
  function(Rceattle,
           file = NULL,
           species = NULL) {

    # Make sure we are using only one model
    if(class(Rceattle) != "Rceattle"){
      stop("Please only use one Rceattle model")
    }
    data_list <- Rceattle$data_list

    # Species names
    if(is.null(species)){
      species =  Rceattle$data_list$spnames
    }

    # Get colors
    colvec=c("red", "blue", "black")

    # * Extract data objects ----
    # - Get observed
    comp_data <- Rceattle$data_list$diet_data
    # - Get estimated
    comp_data$Est = Rceattle$quantities$diet_hat[,2]

    comp_data <- comp_data %>%
      dplyr::mutate(pearson = (Stomach_proportion_by_weight - Est)/ sqrt( ( Est * (1 - Est)) / Sample_size))

    # If year == 0, diet data are averaged from suit_styr to suit_endyr
    # If prey_age >= 0 diet data are diet proportion of prey-at-age in predator-at-age
    # If prey_age < 0 diet data are diet proportion of prey-spp in predator-at-age (sum across prey ages)
    # If prey_age < 0 and pred_age < 0, diet data are mean diet proportion of prey-spp in predator-spp (sum across prey ages and take mean across predator ages)
    # If prey_age < 0 and pred_age < -500, diet data are weighted mean diet proportion of prey-spp in predator-spp (sum across prey ages and take weighted mean across predator ages)


    # Loop around predators ----
    for(pred in 1:data_list$nspp) {
      for(pred_sex in 1:data_list$nsex[pred]){
        for(prey in 1:data_list$nspp) {

          # * Get sex for legend ----
          if(data_list$nsex[pred] > 1){
            pred_legend <- paste("Pred-", species[pred], ifelse(sex == 1, "female", "male"))
            pred_sex = pred_sex - 1
          } else{
            pred_legend <- paste("Pred-", species[pred])
            pred_sex = 0
          }

          # * Extract comps ----
          comp_tmp <- comp_data %>%
            dplyr::filter(Pred == pred & Pred_sex == pred_sex & Prey == prey)

          # - Years
          yrs <- sort(unique(comp_tmp$Year))
          nyrs <- length(yrs)

          # - Min and max
          range_comp <- range(c(comp_tmp$Stomach_proportion_by_weight, comp_tmp$Est))
          range_pearson <- range(comp_tmp$pearson)

          # * Plot annual comps ----
          for(yr in 1:nyrs){

            # Subset year for observed and predicted comp
            comp_tmp_yr <- comp_tmp %>%
              dplyr::filter(Year == yrs[yr] ) %>%
              dplyr::mutate(Prey_age = ifelse(Prey_sex == 2, -Prey_age, Prey_age))


            plot_obs <- comp_tmp_yr %>%
              dplyr::filter(Stomach_proportion_by_weight > 0) %>%
              ggplot2::ggplot(ggplot2::aes(x = Pred_age, y = Prey_age, size = Stomach_proportion_by_weight)) +
              ggplot2::geom_point(alpha = 1) +
              # scale_size(range = c(range_comp[1], range_comp[2]), name="Population (M)") +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Observed diet: year", yrs[yr])) +
              theme(legend.position = c(0.25, 0.7),
                    legend.title = element_blank())

            plot_est <- comp_tmp_yr %>%
              dplyr::filter(Stomach_proportion_by_weight > 0) %>%
              ggplot2::ggplot(ggplot2::aes(x = Pred_age, y = Prey_age, size = Est)) +
              ggplot2::geom_point(alpha = 1) +
              # scale_size(range = c(range_comp[1], range_comp[2]), name="Population (M)") +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Estimated diet: year", yrs[yr])) +
              theme(legend.position = "none")

            plot_pear <- comp_tmp_yr %>%
              dplyr::filter(Stomach_proportion_by_weight > 0) %>%
              ggplot2::ggplot(ggplot2::aes(x = Pred_age, y = Prey_age, size = abs(pearson), color = pearson < 0)) +
              ggplot2::geom_point(alpha = 1) +
              # scale_size(range = c(range_comp[1], range_comp[2]), name="Population (M)") +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Pearson residual: year", yrs[yr])) +
              theme(legend.position = c(0.25, 0.8))

            p1 <- cowplot::plot_grid(plot_obs, plot_est, plot_pear, nrow = 1)
            print(p1)

            # Save ----
            if (!is.null(file)) {
              filename <- paste0(file, "_aggregated_diet_comps_histograms_year", yr,"_", pred_legend, "prey", species[prey],".png")
              ggplot2::ggsave(filename = filename,
                              plot = p1,
                              width = 10,
                              height = 6.5,
                              units = "in",
                              dpi = 300
              )
            }
          }
        }
      }
    }

    # End
  }

