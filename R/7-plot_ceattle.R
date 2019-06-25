#' plot_biomass
#'
#' @description Function the plots the biomass and spawning stock biomass trends as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#'
#' @return Returns and saves a figure with the population trajectory.
#' @export
plot_biomass <-
  function(Rceattle,
           tmp_list = NULL,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
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
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp

    # Get biomass
    Biomass <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      Biomass[, 1:length(Years[[i]]), i] <- Rceattle[[i]]$quantities$biomass[,1:nyrs_vec[i]]
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
      SSB[, 1:length(Years[[i]]), i] <- Rceattle[[i]]$quantities$biomassSSB[,1:nyrs_vec[i]]
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
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_biomass_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

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
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Biomass (million t)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", species[j], bty = "n", cex = 1.4)

        if(!is.null(mohns)){
          legend("top", paste0("B Rho = ", round(mohns[1,j+1], 2), "; SSB Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 1) # Biomass rho
        }

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


        # Mean biomass
        for (k in 1:dim(Biomass)[3]) {
          lines(
            x = Years[[k]],
            y = Biomass[j, 1:length(Years[[k]]), k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
          lines(
            x = Years[[k]],
            y = SSB[j, 1:length(Years[[k]]), k],
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
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param add_ci If the confidence interval is to be added
#' @param lwd Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE, include projection years
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_recruitment <-
  function(Rceattle,
           tmp_list = NULL,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           add_ci = FALSE,
           lwd = 3,
           save_rec = FALSE,
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
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp
    minage <- Rceattle[[1]]$data_list$minage


    # Get biomass
    recruitment <-
      array(NA, dim = c(nspp, nyrs,  length(Rceattle) + length(tmp_list)))
    recruitment_sd <-
      array(NA, dim = c(nspp, nyrs,  length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      recruitment[, 1:length(Years[[i]]) , i] <- Rceattle[[i]]$quantities$R[,1:nyrs_vec[i]]

      # Get SD of rec
      if (add_ci) {
        sd_rec <- which(names(Rceattle[[i]]$sdrep$value) == "R")
        sd_rec <- Rceattle[[i]]$sdrep$sd[sd_rec]
        recruitment_sd[, , i] <-
          replace(recruitment_sd[, , i], values = sd_rec[1:(nyrs_vec[i] * nspp)])
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
          paste0(file, "_recruitment_species_", i, ".csv")
        write.csv(dat_new, file = filename)
      }
    }


    # Plot limits
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(recruitment)[1]) {
      if (add_ci) {
        ymax[i] <- max(c(recruitment_upper[i, , ], 0), na.rm = T)
        ymin[i] <- min(c(recruitment_upper[i, , ], 0), na.rm = T)
      } else{
        ymax[i] <- max(c(recruitment[i, , ], 0), na.rm = T)
        ymin[i] <- min(c(recruitment[i, , ], 0), na.rm = T)
      }
    }
    ymax <- ymax + 0.2 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_recruitment_trajectory", ".png")
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
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Recruitment (millions)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line at end yr
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft",
               legend = species[j],
               bty = "n",
               cex = 1.4)

        if(!is.null(mohns)){
          legend("top", paste0("Rho = ", round(mohns[3,j+1], 2) ), bty = "n", cex = 1.4) # Biomass rho
        }

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
        if (add_ci) {
          for (k in 1:dim(recruitment)[3]) {
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(recruitment_upper[j, 1:length(Years[[k]]), k], rev(recruitment_lower[j, 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = 0.4),
              border = NA
            ) # 95% CI
          }
        }

        # Mean recruitment
        for (k in 1:dim(recruitment)[3]) {
          lines(
            x = Years[[k]],
            y = recruitment[j, 1:length(Years[[k]]), k],
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
#' @param file name of a file to identified the files exported by the
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
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = c("Walleye pollock", "Pacific cod", "Arrowtooth flounder"),
           lwd = 3) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    Years <- list()
    for(i in 1:length(Rceattle)){
      Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
    }
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp
    srv_control <- Rceattle[[1]]$data_list$srv_control
    fsh_control <- Rceattle[[1]]$data_list$fsh_control
    nsrv <- nrow(Rceattle[[1]]$data_list$srv_control)
    nfsh <- nrow(Rceattle[[1]]$data_list$fsh_control)
    nages <- Rceattle[[1]]$data_list$nages

    # Get biomass
    srv_selectivity <-
      array(NA, dim = c(nsrv, max(nages), length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      srv_selectivity[, , i] <- Rceattle[[i]]$quantities$srv_sel[,,1]
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
      array(NA, dim = c(nfsh, max(nages), length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      fsh_selectivity[, , i] <- Rceattle[[i]]$quantities$fsh_sel[,,1]
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
      ymax_srv[i] <- max(c(srv_selectivity[i,,], 0), na.rm = T)
      ymin_srv[i] <- min(c(srv_selectivity[i,,], 0), na.rm = T)
    }

    for (i in 1:dim(fsh_selectivity)[1]) {
      ymax_fsh[i] <- max(c(fsh_selectivity[i,, ], 0), na.rm = T)
      ymin_fsh[i] <- min(c(fsh_selectivity[i,, ], 0), na.rm = T)
    }
    ymax_srv <- ymax_srv + 0.15 * ymax_srv
    ymax_fsh <- ymax_fsh + 0.15 * ymax_fsh

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    max_obj <- max(c(nsrv, nfsh))

    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_selectivity", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5, # 150 / 25.4,

          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:((max_obj + 2)*2), nrow = (max_obj + 2), ncol = 2, byrow = FALSE), heights = c(0.1, rep(1, max_obj), 0.3))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )

      plot.new()

      # Survey selectivity
      for(j in 1:nsrv){
        sp <- srv_control$Species[which(srv_control$Survey_code == j)]
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin_srv[j], ymax_srv[j]),
          xlim = c(min(0),  max(nages, na.rm = TRUE)),
          xlab = "Age",
          ylab = "Survey selectivity",
          xaxt = c(rep("n", nsrv - 1), "s")[j]
        )

        if(j == nsrv){
          mtext(side = 1, "Age", cex  = 0.75, line = 2)
        }

        # Mean selectivity
        for (k in 1:dim(srv_selectivity)[3]) {
          lines(
            x = 1:nages[sp],
            y = srv_selectivity[j, 1:nages[sp], k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          )
        }

        # Legends
        legend("bottomright", as.character(srv_control$Survey_name[j]), bty = "n", cex = 1.4)
      }

      # Add empty plots
      if(nsrv >= nfsh){
        plot.new()
      }

      if(nsrv < nfsh){
        for(i in 1:(nfsh - nsrv)){
          plot.new()
        }
      }


      # Fishery selectivity
      plot.new()
      for (j in 1:nfsh) {
        sp <- fsh_control$Species[which(fsh_control$Fishery_code == j)]
        plot(
          y = NA,
          x = NA,
          ylim = c(ymin_fsh[j], ymax_fsh[j]),
          xlim = c(min(0), max(nages, na.rm = TRUE)),
          xlab = "Age",
          ylab = "Fishery selectivity",
          xaxt = c(rep("n", nfsh - 1), "s")[j]
        )

        if(j == nfsh){
          mtext(side = 1, "Age", cex  = 0.75, line = 2)
        }

        # Mean selectivity
        for (k in 1:dim(fsh_selectivity)[3]) {
          lines(
            x = 1:nages[sp],
            y = fsh_selectivity[j, 1:nages[sp], k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          )
        }

        # Legends
        legend("bottomright", as.character(fsh_control$Fishery_name[j]), bty = "n", cex = 1.4)

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
      }



      if (i == 2) {
        dev.off()
      }
    }
  }


#' Plot functional form
#'
#' @description Function to plot the functional form estimated or specified by \code{\link{Rceattle}}
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
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param age Age specified
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param M2_only TRUE/FALSE plot only M1 and M2
#' @param incl_proj TRUE/FALSE include projection years?
#'
#'
#' @return Returns and saves a figure with the population trajectory.
#' @export
plot_mort <-
  function(Rceattle,
           tmp_list = NULL,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           lwd = 3,
           age = 3,
           right_adj = 0,
           M2_only = FALSE,
           incl_proj = FALSE){

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    Years <- list()
    Endyrs <- list()
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp
    max_age <- max(Rceattle[[1]]$data_list$nages)

    # Species names
    if(is.null(species)){
      species =  Rceattle[[1]]$data_list$spnames
    }

    # Get M2
    M2 <-
      array(NA, dim = c(nspp, max_age, nyrs, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      M2[, , 1:length(Years[[i]]),i] <- Rceattle[[i]]$quantities$M2[,,1:nyrs_vec[i]]
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
      f_mat[, , 1:length(Years[[i]]) ,i] <- Rceattle[[i]]$quantities$F[,,1:nyrs_vec[i]]
    }

    m_mat <-
      array(NA, dim = c(nspp, max_age, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      m_mat[, ,i] <- Rceattle[[i]]$quantities$M1[,1:max_age]
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
      if(M2_only){
        ymax[i] <- max(c(M2[i, age, ,], m_mat[i, age,], 0), na.rm = T)
        ymin[i] <- min(c(M2[i, age, ,], m_mat[i, age,], 0), na.rm = T)
      } else{
        ymax[i] <- max(c(M2[i, age, ,], f_mat[i, age, ,], m_mat[i, age,], 0), na.rm = T)
        ymin[i] <- min(c(M2[i, age, ,], f_mat[i, age, ,], m_mat[i, age,], 0), na.rm = T)
      }
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
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file,"_age",age, "_mortality_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

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
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "Mortality",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line at end yr
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

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
          if(M2_only){
            legend(
              "topright",
              legend = c("M1", "M2"),
              lty = c(3, 1),
              lwd = lwd,
              col = c(1, 1),
              bty = "n",
              cex = 1.175
            )
          } else{
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
        }



        # Mean M2
        for (k in 1:dim(M2)[4]) {
          #if(fishing == FALSE){
          lines(
            x = Years[[k]],
            y = M2[j, age, 1:length(Years[[k]]) , k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
          #}
          #if(fishing){
          if(M2_only == FALSE){
            lines(
              x = Years[[k]],
              y = f_mat[j, age, 1:length(Years[[k]]), k],
              lty = 2,
              lwd = lwd,
              col = line_col[k]
            ) # Median
          }
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


#' Plot maturity
#'
#' @description Function the plots the maturity of each species
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
plot_maturity <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           lwd = 3) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    # Extract data objects
    maturity <- list()
    for(i in 1:length(Rceattle)){
      maturity[[i]] <- Rceattle[[i]]$data_list$pmature
    }

    nspp <- Rceattle[[1]]$data_list$nspp
    nages <- Rceattle[[1]]$data_list$nages

    # Line colors
    if (is.null(line_col)) {
      line_col <- oce::oce.colorsViridis(length(Rceattle))
    }


    # Species names
    if(is.null(species)){
      species =  Rceattle[[1]]$data_list$spnames
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_maturity", ".png")
        png(
          file = filename ,
          width = 4,# 169 / 25.4,
          height = 6.5, # 150 / 25.4,
          family = "Helvetica" ,
          units = "in",
          res = 300
        )
      }

      # Plot configuration
      layout(matrix(1:((nspp + 2)), nrow = (nspp + 2), ncol = 1, byrow = FALSE), heights = c(0.2, rep(1, nspp), 0.3))
      par(
        mar = c(0, 3 , 0 , 1) ,
        oma = c(0 , 0 , 0 , 0),
        tcl = -0.35,
        mgp = c(1.75, 0.5, 0)
      )

      plot.new()

      # Survey selectivity
      for (j in 1:nspp) {
        plot(
          y = NA,
          x = NA,
          ylim = c(0, 1.1),
          xlim = c(min(0), max(nages, na.rm = TRUE)),
          xlab = "Age",
          ylab = "Maturity",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        if(j == nspp){
          mtext(side = 1, "Age", cex  = 0.75, line = 2)
        }

        # Mean maturity
        for (k in 1:length(maturity)) {
          lines(
            x = 1:nages[j],
            y = maturity[[k]][j, 1:nages[j]],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          )
        }

        # Species legends
        legend("topleft", species[j], bty = "n", cex = 1.4)

        # Model name legends
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
      }

      if (i == 2) {
        dev.off()
      }
    }
  }


#' plot_ssb
#'
#' @description Function the plots the spawning stock biomass trends as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param right_adj How many units of the x-axis to add to the right side of the figure for fitting the legend.
#' @param mohns data.frame of mohn's rows extracted from \code{\link{retrospective}}
#' @param incl_proj TRUE/FALSE include projections years
#' @param add_ci TRUE/FALSE, includes 95 percent confidence interval
#'
#' @return Returns and saves a figure with the ssb trajectory.
#' @export
#'
plot_ssb <-
  function(Rceattle,
           tmp_list = NULL,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           lwd = 3,
           right_adj = 0,
           mohns = NULL,
           incl_proj = FALSE,
           add_ci = FALSE) {

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
    for(i in 1:length(Rceattle)){
      Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
      if(incl_proj == FALSE){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      }
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp


    # Get SSB
    SSB <-
      array(NA, dim = c(nspp, nyrs, length(Rceattle) + length(tmp_list)))
    ssb_sd <- array(NA, dim = c(nspp, nyrs, length(Rceattle) + length(tmp_list)))
    for (i in 1:length(Rceattle)) {
      SSB[, 1:length(Years[[i]]), i] <- Rceattle[[i]]$quantities$biomassSSB[,1:nyrs_vec[i]]

      # Get SD of rec
      if (add_ci) {
        ssb_sd_sub <- which(names(Rceattle[[i]]$sdrep$value) == "biomassSSB")
        ssb_sd_sub <- Rceattle[[i]]$sdrep$sd[ssb_sd_sub]
        ssb_sd[, , i] <-
          replace(ssb_sd[, , i], values = ssb_sd_sub[1:(nyrs_vec[i] * nspp)])
      }
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
    ssb_sd <- ssb_sd / 1000000
    SSB_upper <- SSB + ssb_sd * 1.92
    SSB_lower <- SSB - ssb_sd * 1.92

    # Plot limits
    ymax <- c()
    ymin <- c()
    for (i in 1:dim(SSB)[1]) {
      ymax[i] <- max(c(SSB[i, , ], SSB_lower[i, , ], SSB_upper[i, , ], 0), na.rm = T)
      ymin[i] <- min(c(SSB[i, , ], SSB_lower[i, , ], SSB_upper[i, , ], 0), na.rm = T)
    }
    ymax <- ymax + 0.15 * ymax

    if (is.null(line_col)) {
      line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)
    for (i in 1:loops) {
      if (i == 2) {
        filename <- paste0(file, "_ssb_trajectory", ".png")
        png(
          file = filename ,
          width = 7,# 169 / 25.4,
          height = 6.5,# 150 / 25.4,

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
          xlim = c(minyr, maxyr + right_adj),
          xlab = "Year",
          ylab = "SSB (million t)",
          xaxt = c(rep("n", nspp - 1), "s")[j]
        )

        # Horizontal line
        if(incl_proj){
          abline(v = max_endyr, lwd  = lwd, col = "grey", lty = 2)
        }

        # Legends
        legend("topleft", species[j], bty = "n", cex = 1.4)

        if(!is.null(mohns)){
          legend("top", paste0("SSB Rho = ",  round(mohns[2,j+1], 2) ), bty = "n", cex = 1) # SSB rho
        }

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



        # Mean SSB
        for (k in 1:dim(SSB)[3]) {
          lines(
            x = Years[[k]],
            y = SSB[j, 1:length(Years[[k]]), k],
            lty = 1,
            lwd = lwd,
            col = line_col[k]
          ) # Median
        }

        # Credible interval
        if (add_ci) {
          for (k in 1:dim(SSB)[3]) {
            polygon(
              x = c(Years[[k]], rev(Years[[k]])),
              y = c(SSB_upper[j, 1:length(Years[[k]]), k], rev(SSB_lower[j, 1:length(Years[[k]]), k])),
              col = adjustcolor( line_col[k], alpha.f = 0.4),
              border = NA
            )
          }
        }

      }


      if (i == 2) {
        dev.off()
      }
    }
  }

