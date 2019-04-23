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

    #############################################
    # Plot comps Type 1
    #############################################
    for(comp_type in c(0,1)){

      fsh <- unique(fsh_list$Fishery_code[which(fsh_list$Age0_Length1 == comp_type)])
      nfsh <- length(fsh)

      loops <- ifelse(is.null(file), 1, 2)
      for (i in 1:loops) {
        if (i == 2) {
          filename <- paste0(file, c("_fishery_age_comps_1", "_fishery_length_comps_1")[comp_type + 1], ".png")
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
          comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]
          comp_tmp <- comp_tmp[which(comp_tmp$comp > 0),]

          comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
          comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))
          comp_hat_tmp <- comp_hat_tmp[which(!is.na(comp_tmp$comp)),]
          comp_hat_tmp <- comp_hat_tmp[which(comp_hat_tmp$comp > 0),]

          sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
          nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

          plot(
            y = NA,
            x = NA,
            ylim = c(0, nages[sp] * 1.20),
            xlim = c(min(comp_tmp$Year, na.rm = TRUE), max(comp_tmp$Year, na.rm = TRUE) + right_adj),
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


          # Type
          if (j == 1) {
            legend(
              "topright",
              legend = c("Observed", "Estimated"),
              pch = c(16, 16),
              cex = 1.125,
              col = line_col,
              bty = "n"
            )

            x_loc <- c(mean(comp_tmp$Year, na.rm = T) - 2, mean(comp_tmp$Year, na.rm = T), mean(comp_tmp$Year, na.rm = T) +2)
            symbols( x = x_loc , y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.07, 3) , circle = c(0.1, 0.25, 0.5), inches=0.10,add=T, fg = line_col[2])
            text(x = x_loc, y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.16, 3), labels = c(0.1, 0.25, 0.5))
          }


          if(nrow(comp_tmp) > 0){
            # Observed
            symbols( x = comp_tmp$Year , y = comp_tmp$age , circle = comp_tmp$comp, inches=0.10,add=T, fg = line_col[1])

            # Estimated
            symbols( x = comp_hat_tmp$Year , y = comp_hat_tmp$age , circle = comp_hat_tmp$comp, inches=0.10,add=T, fg = line_col[2])
          }
        }


        if (i == 2) {
          dev.off()
        }
      }
    }


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
          comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]
          # comp_tmp <- comp_tmp[which(comp_tmp$comp > 0),]

          comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
          comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))
          comp_hat_tmp <- comp_hat_tmp[which(!is.na(comp_tmp$comp)),]
          # comp_hat_tmp <- comp_hat_tmp[which(comp_hat_tmp$comp > 0),]

          # Calculate pearson residual
          comp_tmp$comp_hat <- comp_hat_tmp$comp

          comp_tmp$pearson <- (comp_tmp$comp - comp_tmp$comp_hat) / sqrt( ( comp_tmp$comp_hat * (1 - comp_tmp$comp_hat)) / comp_tmp$Sample_size)

          max_pearson <- max(abs(comp_tmp$pearson), na.rm = TRUE)


          sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
          nages <- list(Rceattle$data_list$nages, Rceattle$data_list$nlengths)[[comp_type + 1]]

          plot(
            y = NA,
            x = NA,
            ylim = c(0, nages[sp] * 1.20),
            xlim = c(min(comp_tmp$Year, na.rm = TRUE), max(comp_tmp$Year, na.rm = TRUE) + right_adj),
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
          x_loc <- c(mean(comp_tmp$Year, na.rm = T) + 1, mean(comp_tmp$Year, na.rm = T) + 3, mean(comp_tmp$Year, na.rm = T) + 5)
          symbols( x = x_loc , y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.07, 3) , circle = c(0.1, 0.25, 0.5), inches=0.10,add=T, bg = line_col[2])
          text(x = x_loc, y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.16, 3), labels = round(seq(from = 0.5, to = max_pearson, length.out = 3) , 2))

          # Negative
          x_loc <- c(mean(comp_tmp$Year, na.rm = T) - 1, mean(comp_tmp$Year, na.rm = T) - 3, mean(comp_tmp$Year, na.rm = T) - 5)
          symbols( x = x_loc , y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.07, 3) , circle = c(0.1, 0.25, 0.5), inches=0.10,add=T, bg = line_col[1])
          text(x = x_loc, y = rep(max(comp_hat_tmp$age, na.rm = TRUE) * 1.16, 3), labels = round(seq(from = -0.5, to = -max_pearson, length.out = 3) , 2) )



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
  }
