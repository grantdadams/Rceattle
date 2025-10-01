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



#' Plot diet composition fits (Final Version)
#'
#' @description Creates diagnostic plots for diet composition fits. It automatically
#' detects the data aggregation level for each predator-prey interaction and
#' generates the most appropriate plot type.
#'
#' @param Rceattle A single Rceattle model object.
#' @param file Optional file path prefix for saving plots.
#' @param species Optional character vector of species names.
#'
#' @return Prints plots and invisibly returns a list of plot objects.
#' @export
plot_diet_comp <- function(Rceattle, file = NULL, species = NULL) {

  # 1. SETUP & DATA PREPARATION ----
  if(!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for this function.")
  if(!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required for this function.")
  if(!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required for this function.")
  if(!requireNamespace("cowplot", quietly = TRUE)) stop("cowplot is required for this function.")
  if (!inherits(Rceattle, "Rceattle")) stop("Input 'Rceattle' must be a single object of class Rceattle.")
  if (is.null(species)) species <- Rceattle$data_list$spnames

  plot_data <- match_diet_preds(
    data_list = Rceattle$data_list,
    quantities = Rceattle$quantities
  )
  if (is.null(plot_data) || nrow(plot_data) == 0) { message("No diet data to plot."); return(invisible(NULL)) }

  plot_data <- plot_data %>%
    dplyr::rename(Observed = Stomach_proportion_by_weight, Est = Est_diet) %>%
    dplyr::mutate(
      Est_clipped = pmin(0.9999, pmax(0.0001, Est)),
      Pearson = (Observed - Est) / sqrt((Est_clipped * (1 - Est_clipped)) / Sample_size),
      AbsPearson = abs(Pearson)
    )

  # 2. PLOTTING LOGIC ----
  plot_list <- list()

  # Loop through each predator-prey interaction
  for (pred_ind in 1:Rceattle$data_list$nspp) {
    for (prey_ind in 1:Rceattle$data_list$nspp) {

      subset_data <- plot_data %>% dplyr::filter(Pred == pred_ind, Prey == prey_ind)
      if(nrow(subset_data) == 0) next

      is_prey_age_agg <- any(subset_data$Prey_age < 0)
      is_pred_age_agg <- any(subset_data$Pred_age < 0)

      pred_legend <- paste("Predator:", species[pred_ind])

      # --- PATHWAY FOR PREY-AGE AGGREGATED (Line plot) ---
      if(is_prey_age_agg && !is_pred_age_agg) {

        message(paste("Generating line plot for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))
        long_data <- subset_data %>% tidyr::pivot_longer(cols = c(Observed, Est), names_to = "Source", values_to = "Proportion")

        p <- ggplot2::ggplot(long_data, ggplot2::aes(x = Pred_age, y = Proportion, color = Source, linetype = Source, alpha = Source)) +
          ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 2.5) +
          ggplot2::facet_wrap(~ Year, scales = "free_y", labeller = ggplot2::labeller(Year = ~paste("Year:", .))) +
          ggplot2::scale_color_manual(name = "Source", values = c("Observed" = "black", "Est" = "darkred")) +
          ggplot2::scale_linetype_manual(name = "Source", values = c("Observed" = "dashed", "Est" = "solid")) +
          ggplot2::scale_alpha_manual(name = "Source", values = c("Observed" = 1.0, "Est" = 0.7), guide = "none") +
          ggplot2::labs(x = "Predator Age", y = "Diet Proportion", title = paste("Diet of", species[pred_ind], "on", species[prey_ind])) +
          ggplot2::theme_bw()

        print(p)
        plot_list[[length(plot_list) + 1]] <- p

        # --- PATHWAY FOR PRED & PREY-AGE AGGREGATED (Bar plot) ---
      } else if (is_prey_age_agg && is_pred_age_agg) {

        message(paste("Generating bar plot for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))

        plot_data_long <- subset_data %>%
          tidyr::pivot_longer(cols = c(Observed, Est), names_to = "Source", values_to = "Proportion") %>%
          mutate(Prey_name = species[Prey])

        p_fit <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = factor(Year), y = Proportion, fill = Source)) +
          ggplot2::geom_bar(stat = "identity", position = "dodge") +
          ggplot2::scale_fill_manual(name = "Source", values = c("Observed" = "grey50", "Est" = "red")) +
          ggplot2::labs(x = "Year", y = "Diet Proportion", title = paste("Fit to Aggregated Diet:", species[pred_ind], "on", species[prey_ind])) +
          ggplot2::theme_bw()

        print(p_fit)
        plot_list[[length(plot_list) + 1]] <- p_fit

        # --- PATHWAY FOR AGE-DISAGGREGATED (Bubble plot) ---
      } else {

        message(paste("Generating bubble plots for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))
        yrs <- sort(unique(subset_data$Year))

        for(i in 1:length(yrs)) {
          current_yr <- yrs[i]
          comp_tmp_yr <- subset_data %>% dplyr::filter(Year == current_yr)
          if(sum(comp_tmp_yr$Observed, na.rm = TRUE) == 0) next

          # FIX: Use the correct loop variable `prey_ind` instead of `prey`
          title <- paste(pred_legend, "- Prey:", species[prey_ind], "- Year:", current_yr)
          if(current_yr == 0) title <- paste(pred_legend, "- Prey:", species[prey_ind], "(Avg over Years)")

          # Bubble plot logic
          p_obs <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = Pred_age, y = Prey_age, size = Observed)) +
            ggplot2::geom_point(alpha=0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Observed", size = "Prop.")

          p_est <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = Pred_age, y = Prey_age, size = Est)) +
            ggplot2::geom_point(alpha=0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Estimated", size = "Prop.")

          p_pear <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = Pred_age, y = Prey_age, size = AbsPearson, color = Pearson < 0)) +
            ggplot2::geom_point(alpha=0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Pearson Residuals", size = "Abs(Resid)")

          p1 <- cowplot::plot_grid(p_obs, p_est, p_pear, nrow = 1)
          # FIX: Use the correct loop variable `prey_ind`
          p1 <- cowplot::ggdraw(p1) + cowplot::draw_label(title, x = 0.5, y = 0.98)

          print(p1)
          plot_list[[length(plot_list) + 1]] <- p1

          if (!is.null(file)) {
            # FIX: Use the correct loop variables `pred_ind` and `prey_ind`
            # Also removed pred_sex_ind as it is not defined in this loop
            ggplot2::ggsave(paste0(file, "_diet_bubble_Pred", pred_ind, "_Prey", prey_ind, "_Yr", current_yr, ".png"), p1, width = 12, height = 4)
          }
        }
      }
    }
  }
  return(invisible(plot_list))
}
