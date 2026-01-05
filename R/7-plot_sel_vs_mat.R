#' Plot fishery selectivity and maturity
#'
#' @description Function the plots the fishery selectivity and input maturity. Useful for debugging SPR based reference points.
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
plot_selectivity_vs_maturity <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           width = 7,
           height = 6.5,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
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
      selectivity_array[, , , ,i] <- Rceattle[[i]]$quantities$sel_at_age[,,,]
    }

    # Extract data objects
    maturity <- list()
    for(i in 1:length(Rceattle)){
      maturity[[i]] <- Rceattle[[i]]$data_list$pmature[,-1] # Remove species column
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

    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)

    #################################
    # Selectivity time series
    #################################
    for(j in 1:nflt){
      if(fleet_control$Fleet_type[j] == 1){
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
              filename <- paste0(file, "time-varying_selectivity_fleet",j,"_sex",legend_sex, "_w_maturity.png")
              png(
                file = filename ,
                width = width,
                height = height,
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
            res1 <- persp(y = Years[[1]], x =  (1:nages[sp]) - 1 + minage[sp], z = sel_subset, col="white",xlab = "Age",ylab= "\n\nYear", zlab= "\n\nSelectivity",expand=0.5,box=TRUE,ticktype="detailed",phi=35,theta=-19, main = NA, zlim = c(ymin_sel[i], ymax_sel[j]))
            mtext(text = paste(legend_sex2, as.character(fleet_control$Fleet_name[j])), side = 3, cex = 0.8, line = -.5)

            check <- matrix(as.numeric(maturity[[1]][sp, 1:nages[sp]]), ncol = length(Years[[1]]), nrow = nages[sp])
            par(new = TRUE)
            persp(y = Years[[1]], x =  (1:nages[sp]) - 1 + minage[sp], z = check, col=2,xlab = "Age",ylab= "\n\nYear", zlab= "\n\nSelectivity",expand=0.5,box=TRUE,ticktype="detailed",phi=35,theta=-19, main = NA, zlim = c(ymin_sel[i], ymax_sel[j]))
            par(new = TRUE)
            persp(y = Years[[1]], x =  (1:nages[sp]) - 1 + minage[sp], z = sel_subset, col="white",xlab = "Age",ylab= "\n\nYear", zlab= "\n\nSelectivity",expand=0.5,box=TRUE,ticktype="detailed",phi=35,theta=-19, main = NA, zlim = c(ymin_sel[i], ymax_sel[j]))



            lines(
              x = 1:nages[sp],
              y = maturity[[1]][sp, 1:nages[sp]],
              lty = 1,
              lwd = 2,
              col = 1
            )

            if (i == 2) {
              dev.off()
            }
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

              filename <- paste0(file, "_terminal_selectivity_w_maturity_species",sp,"_sex",legend_sex, ".png")
              png(
                file = filename ,
                width = width,
                height = height,
                units = "in",
                res = 300
              )
            }

            fleets <- fleet_control$Fleet_code[which(fleet_control$Species == sp & fleet_control$Fleet_type == 1)]
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


            # Mean maturity
              lines(
                x = 1:nages[sp],
                y = maturity[[1]][sp, 1:nages[sp]],
                lty = 1,
                lwd = lwd,
                col = 1
              )


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
