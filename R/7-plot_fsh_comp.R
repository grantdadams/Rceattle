#' Plot time series of fishery comp data
#'
#' @description Function the plots the fishery comp data as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param cex Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#'
#' @return Returns and saves a figure with the catch trajectory.
#' @export
plot_fsh_comp <-
  function(Rceattle,
           tmp_list = NULL,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           cex = 3,
           lwd = 3,
           right_adj = 0,
           mohns = NULL,
           incl_proj = FALSE) {

    # Make sure we are using only one model
    if(class(Rceattle) != "Rceattle"){
      stop("Please only use one Rceattle model")
    }

    # Species names
    if(is.null(species)){
      species =  Rceattle$data_list$spnames
    }


    # Extract data objects
    # Get observed
    fsh_list <- Rceattle$data_list$fsh_comp
    fsh_list[,grep("Comp_", colnames(fsh_list))] = fsh_list[,grep("Comp_", colnames(fsh_list))]/rowSums(fsh_list[,grep("Comp_", colnames(fsh_list))], na.rm = TRUE)
    fsh_control <- Rceattle$data_list$fsh_control


    # Get estimated
    fsh_hat_list <- Rceattle$data_list$fsh_comp
    fsh_hat_list[,grep("Comp_", colnames(fsh_list))] = Rceattle$quantities$fsh_comp_hat


    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(3))[2:3]
    }

    # #############################################
    # # Plot comps Type 1
    # #############################################
    # for(comp_type in c(0,1)){
    #
    #   fsh <- unique(fsh_list$Fishery_code[which(fsh_list$Age0_Length1 == comp_type)])
    #   nfsh <- length(fsh)
    #
    #   loops <- ifelse(is.null(file), 1, 2)
    #   for (i in 1:loops) {
    #     if (i == 2) {
    #       filename <- paste0(file, c("_fishery_age_comps_1", "_fishery_length_comps_1")[comp_type + 1], ".png")
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
    #     layout(matrix(1:(nfsh + 2), nrow = (nfsh + 2)), heights = c(0.1, rep(1, nfsh), 0.2))
    #     par(
    #       mar = c(0, 3 , 0 , 1) ,
    #       oma = c(0 , 0 , 0 , 0),
    #       tcl = -0.35,
    #       mgp = c(1.75, 0.5, 0)
    #     )
    #     plot.new()
    #
    #     for (j in 1:nfsh) {
    #
    #
    #       # Extract comps
    #       fsh_ind <- which(fsh_list$Fishery_code == fsh[j] & fsh_list$Age0_Length1 == comp_type)
    #       comp_tmp <- fsh_list[fsh_ind,]
    #       comp_hat_tmp <- fsh_hat_list[fsh_ind, ]
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
    #       sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
    #       nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]
    #
    #       plot(
    #         y = NA,
    #         x = NA,
    #         ylim = c(0, nages[sp] * 1.20),
    #         xlim = c(min(fsh_list$Year, na.rm = TRUE), max(fsh_list$Year, na.rm = TRUE) + right_adj),
    #         xlab = "Year",
    #         ylab = c("Fishery age comp", "Fishery length comp")[comp_type + 1],
    #         xaxt = c(rep("n", nfsh - 1), "s")[j]
    #       )
    #
    #
    #       # Horizontal line
    #       if(incl_proj){
    #         abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
    #       }
    #
    #       # Legends
    #       legend("topleft", as.character(fsh_control$Fishery_name[fsh[j]]), bty = "n", cex = 1.4)
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
    #         x_loc <- c(mean(fsh_list$Year, na.rm = T) - 2, mean(fsh_list$Year, na.rm = T), mean(fsh_list$Year, na.rm = T) +2)
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



    #############################################
    # Plot comps Type 2 - Pearson residual
    #############################################
    for(comp_type in c(0,1)){

      fsh <- unique(fsh_list$Fishery_code[which(fsh_list$Age0_Length1 == comp_type)])
      nfsh <- length(fsh)

      loops <- ifelse(is.null(file), 1, 2)
      for (i in 1:loops) {
        if (i == 2) {
          filename <- paste0(file, c("_fishery_age_comps_pearson_residual_1", "_fishery_length_comps_pearson_residual_1")[comp_type + 1], ".png")
          png(
            file = filename ,
            width = 7,# 169 / 25.4,
            height = 6.5,# 150 / 25.4,

            units = "in",
            res = 300
          )
        }

        # Plot configuration
        layout(matrix(1:(nfsh + 2), nrow = (nfsh + 2)), heights = c(0.1, rep(1, nfsh), 0.2))
        par(
          mar = c(0, 3 , 0 , 1) ,
          oma = c(0 , 0 , 0 , 0),
          tcl = -0.35,
          mgp = c(1.75, 0.5, 0)
        )
        plot.new()

        for (j in 1:nfsh) {


          # Extract comps
          fsh_ind <- which(fsh_list$Fishery_code == fsh[j] & fsh_list$Age0_Length1 == comp_type)
          comp_tmp <- fsh_list[fsh_ind,]
          comp_hat_tmp <- fsh_hat_list[fsh_ind, ]

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


          sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
          nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

          plot(
            y = NA,
            x = NA,
            ylim = c(0, nages[sp] * 1.25),
            xlim = c(min(fsh_list$Year, na.rm = TRUE), max(fsh_list$Year, na.rm = TRUE) + right_adj),
            xlab = "Year",
            ylab = c("Fishery age comp", "Fishery length comp")[comp_type + 1],
            xaxt = c(rep("n", nfsh - 1), "s")[j]
          )


          # Horizontal line
          if(incl_proj){
            abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
          }

          # Legends
          legend("topleft", as.character(fsh_control$Fishery_name[fsh[j]]), bty = "n", cex = 1.4)


          if(j == 1){
            legend("topright", "Pearson residual", bty = "n", cex = 1)
          }


          # Positive
          x_loc <- c(mean(fsh_list$Year, na.rm = T) + 1, mean(fsh_list$Year, na.rm = T) + 3.5, mean(fsh_list$Year, na.rm = T) + 6)
          symbols( x = x_loc , y = rep(max(comp_tmp$age, na.rm = TRUE) * 1.1, 3) , circle = round(seq(from = 0.5, to = max_pearson, length.out = 3) , 1), inches=0.10,add=T, bg = line_col[2])
          text(x = x_loc, y = rep(max(comp_tmp$age, na.rm = TRUE) * 1.23, 3), labels = round(seq(from = 0.5, to = max_pearson, length.out = 3) , 1))

          # Negative
          x_loc <- c(mean(fsh_list$Year, na.rm = T) - 1, mean(fsh_list$Year, na.rm = T) - 3.5, mean(fsh_list$Year, na.rm = T) - 6)
          symbols( x = x_loc , y = rep(max(comp_tmp$age, na.rm = TRUE) * 1.1, 3) , circle = -round(seq(from = -0.5, to = -max_pearson, length.out = 3) , 1), inches=0.10,add=T, bg = line_col[1])
          text(x = x_loc, y = rep(max(comp_tmp$age, na.rm = TRUE) * 1.23, 3), labels = round(seq(from = -0.5, to = -max_pearson, length.out = 3) , 1) )



          if(nrow(comp_tmp) > 0){
            # Positive
            symbols( x = comp_tmp$Year , y = comp_tmp$age , circle = abs(comp_tmp$pearson), inches=0.1,add=T, bg = ifelse(comp_tmp$pearson > 0, line_col[2], line_col[1]))
          }
        }


        if (i == 2) {
          dev.off()
        }
      }
    }


    #############################################
    # Plot comps Type 3 - Histograms
    #############################################
    for(comp_type in c(0,1)){

      fsh <- unique(fsh_list$Fishery_code[which(fsh_list$Age0_Length1 == comp_type)])
      nfsh <- length(fsh)

      loops <- ifelse(is.null(file), 1, 2)

      for (j in 1:nfsh) {
        for (i in 1:loops) {
          if (i == 2) {
            filename <- paste0(file,"_",as.character(fsh_control$Fishery_name[fsh[j]]), "_", c("fsh_age_comps_histograms", "fsh_length_comps_histograms")[comp_type + 1], ".png")
            png(
              file = filename ,
              width = 7.5,# 169 / 25.4,
              height = 10,# 150 / 25.4,

              units = "in",
              res = 300
            )
          }


          # Extract comps
          fsh_ind <- which(fsh_list$Fishery_code == fsh[j] & fsh_list$Age0_Length1 == comp_type)
          comp_tmp <- fsh_list[fsh_ind,]
          comp_hat_tmp <- fsh_hat_list[fsh_ind, ]

          # Reorganize and clean
          comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
          comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))

          comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
          comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))

          # Combine
          comp_tmp$comp_hat <- comp_hat_tmp$comp
          comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]
          # comp_tmp <- comp_tmp[which(comp_tmp$comp_hat > 0),]

          # Get comp dims
          sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
          nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

          yrs <- sort(unique(c(comp_hat_tmp$Year, comp_tmp$Year)))
          nyrs <- length(yrs)

          # Min and max
          max_comp <- max(c(comp_tmp$comp, comp_tmp$comp_hat))
          min_comp <- min(c(comp_tmp$comp, comp_tmp$comp_hat))

          # Plot configuration
          plot_rows <- ceiling(nyrs/4) + 1
          layout(matrix(1:(plot_rows * 5), ncol = 5, byrow = FALSE), widths = c(0.2,1,1,1,1))
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

            # Set plot up
            plot(
              y = NA,
              x = NA,
              ylim = c(0, max_comp * 1.10),
              xlim = c(0, nages[sp]),
              xlab = NA,
              ylab = NA,
              xaxt = "n",
              yaxt = "n"
            )

            # x-axis
            if(row_cnt == (plot_rows - 1) | yr == nyrs){
              axis(side = 1)
              mtext(text =  paste(as.character(fsh_control$Fishery_name[fsh[j]]), c("age", "length"))[comp_type + 1],
                    side = 1, line = 2, cex = 0.7)
            }


            # y-axis
            if(yr < plot_rows){
              y_axis <- round(seq(0, max_comp * 1.15, length.out = 4)[1:3],1)
              axis(side = 2, at = y_axis, labels = c(0,y_axis[2:3]))
              mtext(text =  "fsh comp",
                    side = 2, line = 2, cex = 0.7)
            }

            # Legends
            legend("topleft", legend = yrs[yr], bty = "n", cex = 1.2)

            if(yr == 1){
              legend("topright", legend = c("Observed", "Estimated"), lty = 1, lwd = 2, col = c("grey80", 1), bty = "n")
            }

            # Subset year for observed and predicted comp
            comp_tmp_yr <- comp_tmp[which(comp_tmp$Year == yrs[yr]),]

            # Plot observed and predicted comp
            polygon(c(0,comp_tmp_yr$age, max(comp_tmp_yr$age) + 1), c(0, comp_tmp_yr$comp, 0),col='grey80',border=NA)
            lines(c(0,comp_tmp_yr$age, max(comp_tmp_yr$age) + 1), c(0, comp_tmp_yr$comp_hat, 0),col=1, lwd = lwd)


            # Make bottom row of empty plots
            if(row_cnt == (plot_rows - 1)){
              plot.new()
              row_cnt = 0
            }
          }
        }


        if (i == 2) {
          dev.off()
        }
      }
    }



    #############################################
    # Plot comps Type 4 - Histograms of aggregated comps
    #############################################
    for(comp_type in c(0,1)){

      fsh <- unique(fsh_list$Fishery_code[which(fsh_list$Age0_Length1 == comp_type)])
      nfsh <- length(fsh)

      loops <- ifelse(is.null(file), 1, 2)


      for (i in 1:loops) {
        if (i == 2) {
          filename <- paste0(file,"_", c("fsh_age_agg_comps_histograms", "fsh_length_agg_comps_histograms")[comp_type + 1], ".png")
          png(
            file = filename ,
            width = 7.5,# 169 / 25.4,
            height = 10,# 150 / 25.4,

            units = "in",
            res = 300
          )
        }



        # Plot configuration
        if(nfsh < 4){
          layout(matrix(1:(nfsh + 2), nrow = (nfsh + 2), byrow = TRUE), heights = c(0.2, rep(1, nfsh), 0.2))
          par(
            mar = c(2, 3 , 0 , 1) ,
            oma = c(0 , 0 , 0 , 0),
            tcl = -0.35,
            mgp = c(1.75, 0.5, 0)
          )
          plot.new()
          nrows <- nfsh
        }

        if(nfsh >= 4){
          nrows <- ceiling(nfsh/2)
          layout(matrix(1:(((nrows+2) *2)), nrow = (nrows + 2), byrow = TRUE), heights = c(0.2, rep(1, nrows), 0.2))
          par(
            mar = c(2, 3 , 0 , 1) ,
            oma = c(0 , 0 , 0 , 0),
            tcl = -0.35,
            mgp = c(1.75, 0.5, 0)
          )
          plot.new()
          plot.new()
        }

        for (j in 1:nfsh) {

          # Extract comps
          fsh_ind <- which(fsh_list$Fishery_code == fsh[j] & fsh_list$Age0_Length1 == comp_type)
          comp_tmp <- fsh_list[fsh_ind,]
          comp_hat_tmp <- fsh_hat_list[fsh_ind, ]

          # Reorganize and clean
          comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
          comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))

          comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
          comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))

          # Combine
          comp_tmp$comp_hat <- comp_hat_tmp$comp
          comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]


          # Combine across year for observed and predicted comp
          comp_tmp_sum <- aggregate(comp_tmp$comp, by=list(Category=comp_tmp$age), FUN=sum)
          comp_tmp_sum$x <- comp_tmp_sum$x / sum(comp_tmp_sum$x)

          comp_hat_tmp_sum <- aggregate(comp_tmp$comp_hat, by=list(Category=comp_tmp$age), FUN=sum)
          comp_hat_tmp_sum$x <- comp_hat_tmp_sum$x / sum(comp_hat_tmp_sum$x)


          # Get comp dims
          sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
          nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

          yrs <- sort(unique(c(comp_hat_tmp$Year, comp_tmp$Year)))
          nyrs <- length(yrs)


          # Min and max
          max_comp <- max(c(comp_hat_tmp_sum$x, comp_tmp_sum$x))
          min_comp <- min(c(comp_hat_tmp_sum$x, comp_tmp_sum$x))

          # Set plot up
          plot(
            y = NA,
            x = NA,
            ylim = c(0, max_comp * 1.10),
            xlim = c(0, nages[sp]),
            ylab = "Comp",
            xlab = NA
          )


          # Legend
          legend("topleft", legend = as.character(fsh_control$Fishery_name[fsh[j]]), bty = "n")

          if(j == 1){
            legend("topright", legend = c("Observed", "Estimated"), lty = 1, lwd = 2, col = c("grey80", 1), bty = "n")
          }

          # x-axis
          if(j == nfsh){
            axis(side = 1)
            mtext(text =  paste(c("Age", "Length"))[comp_type + 1],
                  side = 1, line = 2, cex = 0.7)
          }


          # Plot observed and predicted comp
          polygon(c(0, comp_tmp_sum$Category, max(comp_tmp_sum$Category) + 1), c(0, comp_tmp_sum$x, 0),col='grey80',border=NA)
          lines(c(0, comp_hat_tmp_sum$Category, max(comp_tmp_sum$Category) + 1), c(0, comp_hat_tmp_sum$x, 0),col=1, lwd = lwd)
        }



        if (i == 2) {
          dev.off()
        }
      }
    }
  }
