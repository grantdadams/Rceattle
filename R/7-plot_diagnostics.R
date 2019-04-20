
#' Plot time series of IOA
#'
#' @description Function the plots the survey indices as estimated from Rceattle
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
#' @return Returns and saves a figure with the index trajectory.
#' @export
plot_index <-
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

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Species names
    if(is.null(species)){
      species =  Rceattle[[1]]$data_list$spnames
    }


    # Extract data objects
    Years <- list()
    Endyrs <- list()
    Srv_list <- list()
    Srv_hat_list <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }

      # Get observed
      Srv_list[[i]] <- Rceattle[[i]]$data_list$srv_biom
      Srv_list[[i]]$CV <- Rceattle[[i]]$quantities$srv_cv_hat
      Srv_list[[i]]$Upper95 <- qlnorm(0.975, meanlog = log(Srv_list[[i]]$Observation), sdlog = Srv_list[[i]]$CV)
      Srv_list[[i]]$Lower95 <- qlnorm(0.025, meanlog = log(Srv_list[[i]]$Observation), sdlog = Srv_list[[i]]$CV)


      # Get estimated
      Srv_hat_list[[i]] <- Rceattle[[i]]$data_list$srv_biom
      Srv_hat_list[[i]]$Observation <- Rceattle[[i]]$quantities$srv_bio_hat
      Srv_hat_list[[i]]$CV <- Rceattle[[i]]$quantities$srv_cv_hat
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp



    # Plot limits
    nsrv <- nrow(Rceattle[[1]]$data_list$srv_control)
    srv_control <- (Rceattle[[1]]$data_list$srv_control)
    ymax <- rep(0, nsrv)
    ymin <- rep(0, nsrv)
    for(i in 1:length(Rceattle)){
      for(srv in 1:nrow(Rceattle[[i]]$data_list$srv_control)){
        srv_ind <- which(Srv_list[[i]]$Survey_code == srv)
        ymax[srv] <- max(c(Srv_list[[i]]$Upper95[srv_ind], Srv_hat_list[[i]]$Observation[srv_ind], ymax[srv]), na.rm = T)
        ymin[srv] <- min(c(Srv_list[[i]]$Lower95[srv_ind], Srv_hat_list[[i]]$Observation[srv_ind], ymin[srv]), na.rm = T)
      }
    }
    ymax <- ymax + 0.1 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_survey_index", ".png")
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
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[j], ymax[j]),
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Relative abundance",
          xaxt = c(rep("n", nsrv - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", srv_control$Survey_name[j], bty = "n", cex = 1.4)

        # if(!is.null(mohns)){
        #   legend("top", paste0("B Rho = ", round(mohns[1,j+1], 2), "; SSB Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 1) # Biomass rho
        # }

        # Model names
        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              pch = rep(16, length(line_col)),
              cex = 1.125,
              col = line_col,
              bty = "n"
            )
          }
        }

        # Type
        if (j == 2) {
            legend(
              "topright",
              legend = c("Observed", "Estimated"),
              pch = c(16, NA),
              lty = c(NA, 1),
              cex = 1.125,
              lwd = lwd,
              col = 1,
              bty = "n"
            )
        }

        # Survey data


        # Mean biomass
        for (k in 1:length(Rceattle)) {
          srv_tmp <- Srv_list[[k]][which(Srv_list[[k]]$Survey_code == j),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_list[[k]]$Survey_code == j),]

          # Observed
          points(
            x = srv_tmp$Year,
            y = srv_tmp$Observation,
            pch = 16,
            cex = cex,
            col = line_col[k]
          ) # Median


          # 95% CI
          arrows(srv_tmp$Year, srv_tmp$Lower95, srv_tmp$Year, srv_tmp$Upper95, length=0.05, angle=90, code=3, col = line_col[k])

          # Estimated
          lines(
            x = srv_hat_tmp$Year,
            y = srv_hat_tmp$Observation,
            pch = 18,
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





#' Plot time series of fishery catch
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
plot_catch <-
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

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Species names
    if(is.null(species)){
      species =  Rceattle[[1]]$data_list$spnames
    }


    # Extract data objects
    Years <- list()
    Endyrs <- list()
    Srv_list <- list()
    Srv_hat_list <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }

      # Get observed
      Srv_list[[i]] <- Rceattle[[i]]$data_list$fsh_biom
      Srv_list[[i]]$CV <- Rceattle[[i]]$quantities$fsh_cv_hat

      no_zero <- which(Srv_list[[i]]$Catch_kg > 0)
      Srv_list[[i]]$Lower95 <- NA
      Srv_list[[i]]$Upper95 <- NA

      Srv_list[[i]]$Upper95[no_zero] <- qlnorm(0.975, meanlog = log(Srv_list[[i]]$Catch_kg[no_zero]), sdlog = Srv_list[[i]]$CV[no_zero])
      Srv_list[[i]]$Lower95[no_zero] <- qlnorm(0.025, meanlog = log(Srv_list[[i]]$Catch_kg[no_zero]), sdlog = Srv_list[[i]]$CV[no_zero])


      # Get estimated
      Srv_hat_list[[i]] <- Rceattle[[i]]$data_list$fsh_biom
      Srv_hat_list[[i]]$Catch_kg <- Rceattle[[i]]$quantities$fsh_bio_hat
      Srv_hat_list[[i]]$CV <- Rceattle[[i]]$quantities$fsh_cv_hat
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp



    # Plot limits
    nsrv <- nrow(Rceattle[[1]]$data_list$fsh_control)
    fsh_control <- (Rceattle[[1]]$data_list$fsh_control)
    ymax <- rep(0, nsrv)
    ymin <- rep(0, nsrv)
    for(i in 1:length(Rceattle)){
      for(srv in 1:nrow(Rceattle[[i]]$data_list$fsh_control)){
        srv_ind <- which(Srv_list[[i]]$Fishery_code == srv)
        ymax[srv] <- max(c(Srv_list[[i]]$Upper95[srv_ind], Srv_hat_list[[i]]$Catch_kg[srv_ind], ymax[srv]), na.rm = T)
        ymin[srv] <- min(c(Srv_list[[i]]$Lower95[srv_ind], Srv_hat_list[[i]]$Catch_kg[srv_ind], ymin[srv]), na.rm = T)
      }
    }
    ymax <- ymax + 0.1 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_total_catch", ".png")
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
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin[j], ymax[j]),
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Total catch",
          xaxt = c(rep("n", nsrv - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = 3, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", fsh_control$Fishery_name[j], bty = "n", cex = 1.4)

        # if(!is.null(mohns)){
        #   legend("top", paste0("B Rho = ", round(mohns[1,j+1], 2), "; SSB Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 1) # Biomass rho
        # }

        # Model names
        if (j == 1) {
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              pch = rep(16, length(line_col)),
              cex = 1.125,
              col = line_col,
              bty = "n"
            )
          }
        }

        # Type
        if (j == 2) {
          legend(
            "topright",
            legend = c("Observed", "Estimated"),
            pch = c(16, NA),
            lty = c(NA, 1),
            cex = 1.125,
            lwd = lwd,
            col = 1,
            bty = "n"
          )
        }

        # Survey data


        # Mean biomass
        for (k in 1:length(Rceattle)) {
          srv_tmp <- Srv_list[[k]][which(Srv_list[[k]]$Fishery_code == j),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_list[[k]]$Fishery_code == j),]

          # Observed
          points(
            x = srv_tmp$Year,
            y = srv_tmp$Catch_kg,
            pch = 16,
            cex = cex,
            col = line_col[k]
          ) # Median


          # 95% CI
          arrows(srv_tmp$Year, srv_tmp$Lower95, srv_tmp$Year, srv_tmp$Upper95, length=0.05, angle=90, code=3, col = line_col[k])

          # Estimated
          lines(
            x = srv_hat_tmp$Year,
            y = srv_hat_tmp$Catch_kg,
            pch = 18,
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
