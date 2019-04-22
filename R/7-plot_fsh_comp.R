#' Plot time series of fishery comp
#'
#' @description Function the plots the fishery catch as estimated from Rceattle
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
    Endyrs <- Rceattle$data_list$endyr
    if(incl_proj == FALSE){
      Years <- Rceattle$data_list$styr:Rceattle$data_list$endyr
    }
    if(incl_proj){
      Years <- Rceattle$data_list$styr:Rceattle$data_list$projyr
    }

    # Get observed
    fsh_list <- Rceattle$data_list$fsh_comp
    fsh_list[,grep("Comp_", colnames(fsh_list))] = fsh_list[,grep("Comp_", colnames(fsh_list))]/rowSums(fsh_list[,grep("Comp_", colnames(fsh_list))], na.rm = TRUE)


    # Get estimated
    fsh_hat_list <- Rceattle$data_list$fsh_comp
    fsh_hat_list[,grep("Comp_", colnames(fsh_list))] = Rceattle$quantities$fsh_comp_hat

    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle$data_list$nspp



    # Plot z-limits
    nfsh <- nrow(Rceattle$data_list$fsh_control)
    fsh_control <- (Rceattle$data_list$fsh_control)
    ymax <- rep(0, nfsh)
    ymin <- rep(0, nfsh)
    for(fsh in 1:nrow(Rceattle$data_list$fsh_control)){
      fsh_ind <- which(fsh_list$Fishery_code == fsh & fsh_list$Age0_Length1 == 0)
      comp_tmp <- as.matrix(fsh_list[fsh_ind, grep("Comp_", colnames(fsh_list))])
      comp_hat_tmp <- as.matrix(fsh_hat_list[fsh_ind, grep("Comp_", colnames(fsh_hat_list))])
      ymax[fsh] <- max(c(comp_tmp, comp_hat_tmp, ymax[fsh]), na.rm = T)
      ymin[fsh] <- min(c(comp_tmp, comp_hat_tmp, ymin[fsh]), na.rm = T)
    }

    ymax <- ymax + 0.1 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(3))[2:3]
    }

    #############################################
    # Plot fishery age comps Type 1
    #############################################
    fsh <- unique(fsh_list$Fishery_code[which(fsh_list$Age0_Length1 == 0)])
    nfsh <- length(fsh)

    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_fishery_age_comps_1", ".png")
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

        sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
        nages <- Rceattle$data_list$nages

        plot(
          y = NA,
          x = NA,
          ylim = c(0, nages[sp] * 1.20),
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Fishery age comp",
          xaxt = c(rep("n", nfsh - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", fsh_control$Fishery_name[fsh[j]], bty = "n", cex = 1.4)


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
        }

        # Extract comps
        fsh_ind <- which(fsh_list$Fishery_code == fsh[j] & fsh_list$Age0_Length1 == 0)
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


    #############################################
    # Plot fishery length comps Type 1
    #############################################
    fsh <- unique(fsh_list$Fishery_code[which(fsh_list$Age0_Length1 == 1)])
    nfsh <- length(fsh)

    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_fishery_length_comps_1", ".png")
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

        sp <- fsh_control$Species[which(fsh_control$Fishery_code == fsh[j])]
        nages <- Rceattle$data_list$nlengths

        plot(
          y = NA,
          x = NA,
          ylim = c(0, nages[sp] * 1.20),
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Fishery length comp",
          xaxt = c(rep("n", nfsh - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", fsh_control$Fishery_name[fsh[j]], bty = "n", cex = 1.4)


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
        }

        # Extract comps
        fsh_ind <- which(fsh_list$Fishery_code == fsh[j] & fsh_list$Age0_Length1 == 1)
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
