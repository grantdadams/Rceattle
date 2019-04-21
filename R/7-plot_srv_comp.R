#' Plot time series of survey age/length composition
#'
#' @description Function the plots the survey composition as estimated from Rceattle
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
plot_srv_comp <-
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
    srv_list <- Rceattle$data_list$srv_comp
    srv_list[,grep("Comp_", colnames(srv_list))] = srv_list[,grep("Comp_", colnames(srv_list))]/rowSums(srv_list[,grep("Comp_", colnames(srv_list))], na.rm = TRUE)


    # Get estimated
    srv_hat_list <- Rceattle$data_list$srv_comp
    srv_hat_list[,grep("Comp_", colnames(srv_list))] = Rceattle$quantities$srv_comp_hat

    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle$data_list$nspp



    # Plot z-limits
    nsrv <- nrow(Rceattle$data_list$srv_control)
    srv_control <- (Rceattle$data_list$srv_control)
    ymax <- rep(0, nsrv)
    ymin <- rep(0, nsrv)
    for(srv in 1:nrow(Rceattle$data_list$srv_control)){
      srv_ind <- which(srv_list$Survey_code == srv & srv_list$Age0_Length1 == 0)
      comp_tmp <- as.matrix(srv_list[srv_ind, grep("Comp_", colnames(srv_list))])
      comp_hat_tmp <- as.matrix(srv_hat_list[srv_ind, grep("Comp_", colnames(srv_hat_list))])
      ymax[srv] <- max(c(comp_tmp, comp_hat_tmp, ymax[srv]), na.rm = T)
      ymin[srv] <- min(c(comp_tmp, comp_hat_tmp, ymin[srv]), na.rm = T)
    }

    ymax <- ymax + 0.1 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(3))[2:3]
    }

    #############################################
    # Plot Survey age comps Type 1
    #############################################
    srv <- unique(srv_list$Survey_code[which(srv_list$Age0_Length1 == 0)])
    nsrv <- length(srv)

    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_survey_age_comps_1", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(nsrv + 2), nrow = (nsrv + 2)), heights = c(0.1, rep(1, nsrv), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      for (j in 1:nsrv) {

        sp <- srv_control$Species[which(srv_control$Survey_code == srv[j])]
        nages <- Rceattle$data_list$nages

        plot(
          y = NA,
          x = NA,
          ylim = c(0, nages[sp] * 1.10),
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Survey age comp",
          xaxt = c(rep("n", nsrv - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", srv_control$Survey_name[sp], bty = "n", cex = 1.4)


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
        srv_ind <- which(srv_list$Survey_code == srv[j] & srv_list$Age0_Length1 == 0)
        comp_tmp <- srv_list[srv_ind,]
        comp_hat_tmp <- srv_hat_list[srv_ind, ]

        # Reorganize and clean
        comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
        comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))
        comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]

        comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
        comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))
        comp_hat_tmp <- comp_hat_tmp[which(!is.na(comp_tmp$comp)),]


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
    # Plot Survey length comps Type 1
    #############################################
    srv <- unique(srv_list$Survey_code[which(srv_list$Age0_Length1 == 1)])
    nsrv <- length(srv)

    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_survey_length_comps_1", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:(nsrv + 2), nrow = (nsrv + 2)), heights = c(0.1, rep(1, nsrv), 0.2))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )
      plot.new()

      for (j in 1:nsrv) {

        sp <- srv_control$Species[which(srv_control$Survey_code == srv[j])]
        nages <- Rceattle$data_list$nlengths

        plot(
          y = NA,
          x = NA,
          ylim = c(0, nages[sp] * 1.10),
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Survey length comp",
          xaxt = c(rep("n", nsrv - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", srv_control$Survey_name[sp], bty = "n", cex = 1.4)


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
        srv_ind <- which(srv_list$Survey_code == srv[j] & srv_list$Age0_Length1 == 1)
        comp_tmp <- srv_list[srv_ind,]
        comp_hat_tmp <- srv_hat_list[srv_ind, ]

        # Reorganize and clean
        comp_tmp <- tidyr::gather(comp_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_tmp)))
        comp_tmp$age <- as.numeric(gsub("Comp_", "", comp_tmp$age))
        comp_tmp <- comp_tmp[which(!is.na(comp_tmp$comp)),]

        comp_hat_tmp <- tidyr::gather(comp_hat_tmp, key = "age", value = "comp", grep("Comp_", colnames(comp_hat_tmp)))
        comp_hat_tmp$age <- as.numeric(gsub("Comp_", "", comp_hat_tmp$age))
        comp_hat_tmp <- comp_hat_tmp[which(!is.na(comp_tmp$comp)),]


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
