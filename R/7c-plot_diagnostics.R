#' CPUE fits
#'
#' Plot of fitted CPUE indices on natural-scale (r4ss-style)
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param right_adj How much right side of the x-axis for fitting the legend. As percentage.
#' @param top_adj How much top side of the y-axis for fitting the legend. As percentage.
#' @param incl_proj TRUE/FALSE include projections years
#' @param width plot width
#' @param height plot hight
#' @export

plot_index <- function(Rceattle,
                       file = NULL,
                       model_names = NULL,
                       line_col = NULL,
                       species = NULL,
                       right_adj = 0,
                       top_adj = 0.05,
                       incl_proj = FALSE,
                       single.plots=FALSE,
                       width=NULL,
                       height=NULL){


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
    Srv_list[[i]]$Log_sd <- Rceattle[[i]]$quantities$srv_log_sd_hat
    Srv_list[[i]]$Upper95 <- qlnorm(0.975, meanlog = log(Srv_list[[i]]$Observation), sdlog = Srv_list[[i]]$Log_sd)
    Srv_list[[i]]$Lower95 <- qlnorm(0.025, meanlog = log(Srv_list[[i]]$Observation), sdlog = Srv_list[[i]]$Log_sd)


    # Get estimated
    Srv_hat_list[[i]] <- Rceattle[[i]]$data_list$srv_biom
    Srv_hat_list[[i]]$Observation <- Rceattle[[i]]$quantities$srv_bio_hat
    Srv_hat_list[[i]]$Log_sd <- Rceattle[[i]]$quantities$srv_log_sd_hat
  }
  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)

  nspp <- Rceattle[[1]]$data_list$nspp


  # Plot limits
  fleet_control <- (Rceattle[[1]]$data_list$fleet_control)
  srv_biom <- (Rceattle[[1]]$data_list$srv_biom)
  srvs <- sort(unique(srv_biom$Fleet_code))
  srvs <- srvs[which(srvs %in% fleet_control$Fleet_code[which(fleet_control$Fleet_type == 2)])] # Only use surveys that are estimates
  #FIXME assumes all surveys are the same across models
  nsrv <- length(srvs)

  ymax <- c()
  ymin <- c()

  for(srv in 1:nsrv){
    for(i in 1:length(Rceattle)){
      srv_ind <- which(Srv_list[[i]]$Fleet_code == srvs[srv])
      ymax[srv] <- max(c(Srv_list[[i]]$Upper95[srv_ind], Srv_hat_list[[i]]$Observation[srv_ind], ymax[srv]), na.rm = T)
      ymin[srv] <- min(c(Srv_list[[i]]$Lower95[srv_ind], Srv_hat_list[[i]]$Observation[srv_ind], ymin[srv]), na.rm = T)
    }
  }
  ymax <- ymax + top_adj * (ymax-ymin)

  # Assume colors if not provided
  if (is.null(line_col)) {
    line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
  }

  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (j in 1:loops) {

    # Plot/save each survey individually
    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(srv in 1:nsrv){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

        # Save
        if(j == 2){
          filename <- paste0(file, "fleet",srvs[j]," ",as.character(fleet_control$Fleet_name[srvs[srv]]), "_survey_index", ".png")
          png(file = filename, width = width, height = height, res = 200, units = "in")
        }


        minyr <- min(sapply(Srv_list, function(x) min(x[which(x$Fleet_code == srvs[srv] & x$Year > 0),]$Year)))
        maxyr <- max(sapply(Srv_list, function(x) max(x[which(x$Fleet_code == srvs[srv]),]$Year)))

        par(Par)
        plot(NA, NA, ylab="Index", xlab="Year", ylim = c((ymin[srv]), (ymax[srv])), xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj), type='n', xaxt="n", yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Loop through models
        for (k in 1:length(Rceattle)) {

          # Subset data by fleet and model
          srv_tmp <- Srv_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]

          # Plot predicted CPUE
          lines(abs(srv_hat_tmp$Year), (srv_hat_tmp$Observation),lwd=2,col=line_col[k])

          # Plot observed CPUE
          gplots::plotCI(srv_tmp$Year, (srv_tmp$Observation), ui=(srv_tmp$Upper95), li=(srv_tmp$Lower95),add=T,gap=0,pch=21,xaxt="n",yaxt="n",pt.bg = "white")
        }

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[srvs[srv]]),bty="n",y.intersp = -0.2,cex=0.8)

        # Model names
        if(!is.null(model_names)){
          legend(
            "topright",
            legend = model_names,
            pch = rep(16, length(line_col)), cex=0.8,
            col = line_col,
            bty = "n"
          )

        }
        # Save plot
        if(j == 2){dev.off()}
      }
    }

    # Plot/save each survey individually
    if(single.plots==FALSE){

      # Set heights of plot
      if(is.null(width)) width = 7
      if(is.null(height)) height = ifelse(nsrv==1,5,ifelse(nsrv==2,3.,2.5))*round(nsrv/2+0.01,0)


      Par = list(mfrow=c(round(nsrv/2+0.01,0),ifelse(nsrv==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)

      # Save
      if(j == 2){
        filename <- paste0(file,"_survey_indices", ".png")
        png(file = filename, width = width, height = height, res = 200, units = "in")
      }
      par(Par)


      for(srv in 1:nsrv){
        minyr <- min(sapply(Srv_list, function(x) min(x[which(x$Fleet_code == srvs[srv] & x$Year > 0),]$Year)))
        maxyr <- max(sapply(Srv_list, function(x) max(x[which(x$Fleet_code == srvs[srv]),]$Year)))

        xlim <- c(minyr, maxyr)
        if(srv == 1){
          xlim <- c(minyr, maxyr + (maxyr - minyr) * right_adj)
        }

        plot(NA, NA, ylab="", xlab="", ylim = c((ymin[srv]), (ymax[srv])), xlim = xlim, type='n', xaxt="n", yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[srvs[srv]]),bty="n",y.intersp = -0.2,cex=0.8)

        # Model names
        if(srv == 1){
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              pch = rep(16, length(line_col)), cex=0.8,
              col = line_col,
              bty = "n"
            )
          }
        }

        # Loop through models
        for (k in 1:length(Rceattle)) {

          # Subset data by fleet and model
          srv_tmp <- Srv_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]

          # Plot predicted CPUE
          lines(abs(srv_hat_tmp$Year), (srv_hat_tmp$Observation),lwd=2,col=line_col[k])

          # Plot observed CPUE
          gplots::plotCI(srv_tmp$Year, (srv_tmp$Observation), ui=(srv_tmp$Upper95), li=(srv_tmp$Lower95),add=T,gap=0,pch=21,xaxt="n",yaxt="n",pt.bg = "white")
        }
      }
      mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      mtext(paste("Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(j == 2){dev.off()}
    }
  }
} # End of fit


#' Landings fits
#'
#' Plot of fitted landings data on natural-scale (r4ss-style)
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param right_adj How much right side of the x-axis for fitting the legend. As percentage.
#' @param top_adj How much top side of the y-axis for fitting the legend. As percentage.
#' @param incl_proj TRUE/FALSE include projections years
#' @param width plot width
#' @param height plot hight
#' @param mse Is if an MSE object from \code{\link{load_mse}} or \code{\link{mse_run}}. Will plot data from OMs.
#' @export

plot_catch <- function(Rceattle,
                       file = NULL,
                       model_names = NULL,
                       line_col = NULL,
                       species = NULL,
                       right_adj = 0,
                       top_adj = 0.05,
                       incl_proj = FALSE,
                       single.plots=FALSE,
                       width=NULL,
                       height=NULL,
                       alpha = 0.4,
                       mse = FALSE){

  # Convert mse object to Rceattle list
  if(mse){
    Rceattle <- lapply(Rceattle, function(x) x$OM)
    incl_proj = TRUE
  }

  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }

  # Species names
  if(is.null(species)){
    species =  Rceattle[[1]]$data_list$spnames
  }


  # Extract data objects
  if(incl_proj == FALSE){
    Years <- lapply(Rceattle, function(x) x$data_list$styr: x$data_list$endyr)
  }
  if(incl_proj){
    Years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }
  Endyrs <- lapply(Rceattle, function(x) x$data_list$endyr)
  meanyrs <- lapply(Rceattle, function(x) x$data_list$meanyr)
  fsh_list <- list()
  fsh_hat_list <- list()
  nmods = length(Rceattle)

  for(i in 1:length(Rceattle)){
    # Get observed
    fsh_list[[i]] <- Rceattle[[i]]$data_list$fsh_biom[which(Rceattle[[i]]$data_list$fsh_biom$Year %in% Years[[i]] ),]
    fsh_list[[i]]$Log_sd <- Rceattle[[i]]$quantities$fsh_log_sd_hat[which(Rceattle[[i]]$data_list$fsh_biom$Year %in% Years[[i]] )]

    no_zero <- which(fsh_list[[i]]$Catch > 0)
    fsh_list[[i]]$Lower95 <- 0
    fsh_list[[i]]$Upper95 <- 0
    fsh_list[[i]]$Upper95[no_zero]  <- qlnorm(0.975, meanlog = log(fsh_list[[i]]$Catch[no_zero]), sdlog = fsh_list[[i]]$Log_sd[no_zero])
    fsh_list[[i]]$Lower95[no_zero]  <- qlnorm(0.025, meanlog = log(fsh_list[[i]]$Catch[no_zero]), sdlog = fsh_list[[i]]$Log_sd[no_zero])


    # Get estimated
    fsh_hat_list[[i]] <- Rceattle[[i]]$data_list$fsh_biom[which(Rceattle[[i]]$data_list$fsh_biom$Year %in% Years[[i]] ),]
    fsh_hat_list[[i]]$Catch <- Rceattle[[i]]$quantities$fsh_bio_hat[which(Rceattle[[i]]$data_list$fsh_biom$Year %in% Years[[i]] )]
    fsh_hat_list[[i]]$Log_sd <- Rceattle[[i]]$quantities$fsh_log_sd_hat[which(Rceattle[[i]]$data_list$fsh_biom$Year %in% Years[[i]] )]
  }


  # - MSE objects
  if(mse){
    # -- Get quantiles and mean across simulations
    catch_list_tmp <- simplify2array(lapply(fsh_hat_list, function(x) x$Catch))
    fsh_hat_list[[1]]$Upper95 <- apply( catch_list_tmp, 1, function(x) quantile(x, probs = 0.975) )
    fsh_hat_list[[1]]$Lower95 <- apply( catch_list_tmp, 1, function(x) quantile(x, probs = 0.025) )
    fsh_hat_list[[1]]$Upper50 <- apply( catch_list_tmp, 1, function(x) quantile(x, probs = 0.75) )
    fsh_hat_list[[1]]$Lower50 <- apply( catch_list_tmp, 1, function(x) quantile(x, probs = 0.25) )
    fsh_hat_list[[1]]$Catch <- apply( catch_list_tmp, 1, function(x) mean(x) ) # Get mean quantity
    nmods = 1
  }

  # Plot
  minyr <- min(unlist(Years), na.rm = TRUE)
  maxyr <- max(unlist(Years), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)

  nspp <- Rceattle[[1]]$data_list$nspp


  # Plot limits
  fleet_control <- (Rceattle[[1]]$data_list$fleet_control)
  fsh_biom <- (Rceattle[[1]]$data_list$fsh_biom)
  fshs <- sort(unique(fsh_biom$Fleet_code))
  nfsh <- length(fshs)

  ymax <- c()
  ymin <- c()

  for(fsh in 1:nfsh){
    for(i in 1:nmods){
      fsh_ind <- which(fsh_list[[i]]$Fleet_code == fshs[fsh])
      ymax[fsh] <- max(c(fsh_list[[i]]$Upper95[fsh_ind], fsh_hat_list[[i]]$Catch[fsh_ind], ymax[fsh]), na.rm = T)
      ymin[fsh] <- min(c(fsh_list[[i]]$Lower95[fsh_ind], fsh_hat_list[[i]]$Catch[fsh_ind], ymin[fsh]), na.rm = T)

      if(mse){
        ymax[fsh] <- max(c(fsh_hat_list[[i]]$Upper95[fsh_ind], ymax[fsh]), na.rm = T)
        ymin[fsh] <- min(c(fsh_hat_list[[i]]$Lower95[fsh_ind], ymin[fsh]), na.rm = T)

      }
    }
  }
  ymax <- ymax + top_adj * (ymax-ymin)

  # Assume colors if not provided
  if (is.null(line_col)) {
    line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
  }

  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (j in 1:loops) {

    # Plot/save each survey individually
    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(fsh in 1:nfsh){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

        # Save
        if(j == 2){
          filename <- paste0(file, "fleet",fshs[j]," ",as.character(fleet_control$Fleet_name[fshs[fsh]]), "_fishery_catch", ".png")
          png(file = filename, width = width, height = height, res = 200, units = "in")
        }

        par(Par)
        plot(NA, NA, ylab="Catch", xlab="Year", ylim = c((ymin[fsh]), (ymax[fsh])), xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj), type='n', xaxt="n", yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Loop through models
        for (k in 1:nmods) {

          # Subset data by fleet and model
          fsh_tmp <- fsh_list[[k]] %>%
            filter(Fleet_code == fshs[fsh])

          if(mse){
            fsh_tmp <- fsh_tmp %>% filter(Year <= meanyrs[k]) # Only show historical catch if MSE models
          }

          fsh_hat_tmp <- fsh_hat_list[[k]] %>%
            filter(Fleet_code == fshs[fsh])


          # - Plot predicted CPUE
          lines(fsh_hat_tmp$Year, (fsh_hat_tmp$Catch),lwd=2,col=line_col[k])

          # - Plot observed CPUE
          gplots::plotCI(fsh_tmp$Year, (fsh_tmp$Catch), ui=(fsh_tmp$Upper95), li=(fsh_tmp$Lower95),add=T,gap=0,pch=21,xaxt="n",yaxt="n",pt.bg = "white")

          # - Plot MSE shading
          if(mse){
            fsh_hat_tmp <- fsh_hat_list[[k]] %>%
              filter(Year > meanyrs[k] & Fleet_code == fshs[fsh])

            # 95% CI
            polygon(
              x = c(fsh_hat_tmp$Year, rev(fsh_hat_tmp$Year)),
              y = c(fsh_hat_tmp$Upper95, rev(fsh_hat_tmp$Lower95)),
              col = adjustcolor( line_col[k], alpha.f = alpha),
              border = NA
            )

            # 50% CI
            polygon(
              x = c(fsh_hat_tmp$Year, rev(fsh_hat_tmp$Year)),
              y = c(fsh_hat_tmp$Upper50, rev(fsh_hat_tmp$Lower50)),
              col = adjustcolor( line_col[k], alpha.f = alpha),
              border = NA
            )

          }
        }

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[fshs[fsh]]),bty="n",y.intersp = -0.2,cex=0.8)

        # Model names
        if(!is.null(model_names)){
          legend(
            "topright",
            legend = model_names,
            pch = rep(16, length(line_col)), cex=0.8,
            col = line_col,
            bty = "n"
          )

        }
        # Save plot
        if(j == 2){dev.off()}
      }
    }

    # Plot/save each survey together
    if(single.plots==FALSE){

      # Set heights of plot
      if(is.null(width)) width = 7
      if(is.null(height)) height = ifelse(nfsh==1,5,ifelse(nfsh==2,3.,2.5))*round(nfsh/2+0.01,0)


      Par = list(mfrow=c(round(nfsh/3+0.01,0),ifelse(nfsh==1,1,3)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)

      # Save
      if(j == 2){
        filename <- paste0(file,"_fishery_catch", ".png")
        png(file = filename, width = width, height = height, res = 200, units = "in")
      }
      par(Par)


      for(fsh in 1:nfsh){

        xlim <- c(minyr, maxyr)
        if(fsh == 1){
          xlim <- c(minyr, maxyr + (maxyr - minyr) * right_adj)
        }

        plot(NA, NA, ylab="", xlab="", ylim = c((ymin[fsh]), (ymax[fsh])), xlim = xlim, type='n', xaxt="n", yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[fshs[fsh]]),bty="n",y.intersp = -0.2,cex=0.8)

        # Model names
        if(fsh == 1){
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              pch = rep(16, length(line_col)), cex=0.8,
              col = line_col,
              bty = "n"
            )
          }
        }

        # Loop through models
        for (k in 1:nmods) {

          # Subset data by fleet and model
          fsh_tmp <- fsh_list[[k]] %>%
            filter(Fleet_code == fshs[fsh])

          if(mse){
            fsh_tmp <- fsh_tmp %>% filter(Year <= meanyrs[k]) # Only show historical catch if MSE models
          }

          fsh_hat_tmp <- fsh_hat_list[[k]] %>%
            filter(Fleet_code == fshs[fsh])


          # - Plot predicted CPUE
          lines(fsh_hat_tmp$Year, (fsh_hat_tmp$Catch),lwd=2,col=line_col[k])

          # - Plot observed CPUE
          gplots::plotCI(fsh_tmp$Year, (fsh_tmp$Catch), ui=(fsh_tmp$Upper95), li=(fsh_tmp$Lower95),add=T,gap=0,pch=21,xaxt="n",yaxt="n",pt.bg = "white")

          # - Plot MSE shading
          if(mse){
            fsh_hat_tmp <- fsh_hat_list[[k]] %>%
              filter(Year > meanyrs[k] & Fleet_code == fshs[fsh])

            # 95% CI
            polygon(
              x = c(fsh_hat_tmp$Year, rev(fsh_hat_tmp$Year)),
              y = c(fsh_hat_tmp$Upper95, rev(fsh_hat_tmp$Lower95)),
              col = adjustcolor( line_col[k], alpha.f = alpha),
              border = NA
            )

            # 50% CI
            polygon(
              x = c(fsh_hat_tmp$Year, rev(fsh_hat_tmp$Year)),
              y = c(fsh_hat_tmp$Upper50, rev(fsh_hat_tmp$Lower50)),
              col = adjustcolor( line_col[k], alpha.f = alpha),
              border = NA
            )

          }
        }
      }
      mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      mtext(paste("Catch"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(j == 2){dev.off()}
    }
  }
} # End of fit






#' log(CPUE) fits
#'
#' Plot of fitted CPUE indices on log-scale (r4ss-style)
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param right_adj How much right side of the x-axis for fitting the legend. As percentage.
#' @param top_adj How much top side of the y-axis for fitting the legend. As percentage.
#' @param incl_proj TRUE/FALSE include projections years
#' @param width plot width
#' @param height plot hight
#' @export

plot_logindex <- function(Rceattle,
                          file = NULL,
                          model_names = NULL,
                          line_col = NULL,
                          species = NULL,
                          right_adj = 0,
                          top_adj = 0.05,
                          incl_proj = FALSE,
                          single.plots=FALSE,
                          width=NULL,
                          height=NULL){


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

    Srv_list[[i]] <- Rceattle[[i]]$data_list$srv_biom[which(Rceattle[[i]]$data_list$srv_biom$Year > 0),]
    Srv_list[[i]]$Log_sd <- Rceattle[[i]]$quantities$srv_log_sd_hat[which(Rceattle[[i]]$data_list$srv_biom$Year > 0)]
    Srv_list[[i]]$Upper95 <- qlnorm(0.975, meanlog = log(Srv_list[[i]]$Observation), sdlog = Srv_list[[i]]$Log_sd)
    Srv_list[[i]]$Lower95 <- qlnorm(0.025, meanlog = log(Srv_list[[i]]$Observation), sdlog = Srv_list[[i]]$Log_sd)


    # Get estimated
    Srv_hat_list[[i]] <- Rceattle[[i]]$data_list$srv_biom[which(Rceattle[[i]]$data_list$srv_biom$Year > 0),]
    Srv_hat_list[[i]]$Observation <- Rceattle[[i]]$quantities$srv_bio_hat[which(Rceattle[[i]]$data_list$srv_biom$Year > 0)]
    Srv_hat_list[[i]]$Log_sd <- Rceattle[[i]]$quantities$srv_log_sd_hat[which(Rceattle[[i]]$data_list$srv_biom$Year > 0)]
  }
  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)

  nspp <- Rceattle[[1]]$data_list$nspp


  # Plot limits
  fleet_control <- (Rceattle[[1]]$data_list$fleet_control)
  srv_biom <- (Rceattle[[1]]$data_list$srv_biom)
  srvs <- sort(unique(srv_biom$Fleet_code))
  srvs <- srvs[which(srvs %in% fleet_control$Fleet_code[which(fleet_control$Fleet_type == 2)])] # Only use surveys that are estimates
  nsrv <- length(srvs)

  ymax <- c()
  ymin <- c()

  for(srv in 1:nsrv){
    for(i in 1:length(Rceattle)){
      srv_ind <- which(Srv_list[[i]]$Fleet_code == srvs[srv])
      ymax[srv] <- max(c(log(c(Srv_list[[i]]$Upper95[srv_ind], Srv_hat_list[[i]]$Observation[srv_ind])), ymax[srv]), na.rm = T)
      ymin[srv] <- min(c(log(c(Srv_list[[i]]$Lower95[srv_ind], Srv_hat_list[[i]]$Observation[srv_ind])), ymin[srv]), na.rm = T)
    }
  }
  ymax <- ymax + top_adj * (ymax-ymin)

  # Assume colors if not provided
  if (is.null(line_col)) {
    line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
  }

  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (j in 1:loops) {

    # Plot/save each survey individually
    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(srv in 1:nsrv){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

        # Save
        if(j == 2){
          filename <- paste0(file, "fleet",srvs[j]," ",as.character(fleet_control$Fleet_name[srvs[srv]]), "_log_survey_index", ".png")
          png(file = filename, width = width, height = height, res = 200, units = "in")
        }


        minyr <- min(sapply(Srv_list, function(x) min(x[which(x$Fleet_code == srvs[srv] & x$Year > 0),]$Year)))
        maxyr <- max(sapply(Srv_list, function(x) max(x[which(x$Fleet_code == srvs[srv]),]$Year)))

        par(Par)
        plot(NA, NA, ylab="log Index", xlab="Year", ylim = c((ymin[srv]), (ymax[srv])), xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj), type='n', xaxt="n", yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Loop through models
        for (k in 1:length(Rceattle)) {

          # Subset data by fleet and model
          srv_tmp <- Srv_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]

          # Plot predicted CPUE
          lines(srv_hat_tmp$Year, log(srv_hat_tmp$Observation),lwd=2,col=line_col[k])

          # Plot observed CPUE
          gplots::plotCI(srv_tmp$Year, log(srv_tmp$Observation), ui=log(srv_tmp$Upper95), li=log(srv_tmp$Lower95),add=T,gap=0,pch=21,xaxt="n",yaxt="n",pt.bg = "white")
        }

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[srvs[srv]]),bty="n",y.intersp = -0.2,cex=0.8)

        # Model names
        if(!is.null(model_names)){
          legend(
            "topright",
            legend = model_names,
            pch = rep(16, length(line_col)), cex=0.8,
            col = line_col, ncol = 2,
            bty = "n"
          )

        }
        # Save plot
        if(j == 2){dev.off()}
      }
    }

    # Plot/save each survey individually
    if(single.plots==FALSE){

      # Set heights of plot
      if(is.null(width)) width = 7
      if(is.null(height)) height = ifelse(nsrv==1,5,ifelse(nsrv==2,3.,2.5))*round(nsrv/2+0.01,0)


      Par = list(mfrow=c(round(nsrv/2+0.01,0),ifelse(nsrv==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)

      # Save
      if(j == 2){
        filename <- paste0(file,"_log_survey_indices", ".png")
        png(file = filename, width = width, height = height, res = 200, units = "in")
      }
      par(Par)


      for(srv in 1:nsrv){
        minyr <- min(sapply(Srv_list, function(x) min(x[which(x$Fleet_code == srvs[srv] & x$Year > 0),]$Year)))
        maxyr <- max(sapply(Srv_list, function(x) max(x[which(x$Fleet_code == srvs[srv]),]$Year)))

        xlim <- c(minyr, maxyr)
        if(srv == 1){
          xlim <- c(minyr, maxyr + (maxyr - minyr) * right_adj)
        }

        plot(NA, NA, ylab="", xlab="", ylim = c((ymin[srv]), (ymax[srv])), xlim = xlim, type='n', xaxt="n", yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[srvs[srv]]),bty="n",y.intersp = -0.2,cex=0.8)


        # Loop through models
        for (k in 1:length(Rceattle)) {

          # Subset data by fleet and model
          srv_tmp <- Srv_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_list[[k]]$Fleet_code == srvs[srv]),]

          # Plot predicted CPUE
          lines(srv_hat_tmp$Year, log(srv_hat_tmp$Observation),lwd=2,col=line_col[k])

          # Plot observed CPUE
          gplots::plotCI(srv_tmp$Year, log(srv_tmp$Observation), ui=log(srv_tmp$Upper95), li=log(srv_tmp$Lower95),add=T,gap=0,pch=21,xaxt="n",yaxt="n",pt.bg = "white")
        }

        # Model names
        if(srv == 1 & (nsrv %% 2) == 0){
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              pch = rep(16, length(line_col)), cex=0.8,
              col = line_col,
              bty = "n"
            )
          }
        }

        if(srv == nsrv & (nsrv %% 2) != 0 & nsrv > 1){
          plot(NA, NA, ylab="", xlab="", ylim = c((ymin[srv]), (ymax[srv])), xlim = xlim, type='n', xaxt="n", yaxt="n", bty = "n")
          if(!is.null(model_names)){
            legend(
              "top",
              legend = model_names,
              pch = rep(16, length(line_col)), cex=0.75,
              col = line_col,
              bty = "n"
            )
          }
        }
      }
      mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(j == 2){dev.off()}
    }
  }
} # End of logfit



#' CPUE residual
#'
#' Plot of residuals CPUE indices on log-scale (r4ss-style)
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param right_adj How much right side of the x-axis for fitting the legend. As percentage.
#' @param top_adj How much top side of the y-axis for fitting the legend. As percentage.
#' @param incl_proj TRUE/FALSE include projections years
#' @param width plot width
#' @param height plot hight
#' @export

plot_indexresidual <- function(Rceattle,
                               file = NULL,
                               model_names = NULL,
                               line_col = NULL,
                               species = NULL,
                               right_adj = 0,
                               top_adj = 0.05,
                               incl_proj = FALSE,
                               single.plots=FALSE,
                               width=NULL,
                               height=NULL){


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
  Srv_hat_list <- list()
  for(i in 1:length(Rceattle)){
    Endyrs[[i]] <- Rceattle[[i]]$data_list$endyr
    if(incl_proj == FALSE){
      Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
    }
    if(incl_proj){
      Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
    }

    # Get estimate and residual
    Srv_hat_list[[i]] <- Rceattle[[i]]$data_list$srv_biom
    Srv_hat_list[[i]]$Estimate <- Rceattle[[i]]$quantities$srv_bio_hat
    Srv_hat_list[[i]]$Residual <- log(Rceattle[[i]]$quantities$srv_bio_hat) - log(Rceattle[[i]]$data_list$srv_biom$Observation)
  }
  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(Years, length)
  nyrs <- max(nyrs_vec)

  nspp <- Rceattle[[1]]$data_list$nspp


  # Plot limits
  fleet_control <- (Rceattle[[1]]$data_list$fleet_control)
  srv_biom <- (Rceattle[[1]]$data_list$srv_biom)
  srvs <- sort(unique(srv_biom$Fleet_code))
  nsrv <- length(srvs)

  ymax <- c()
  ymin <- c()

  for(srv in 1:nsrv){
    for(i in 1:length(Rceattle)){
      srv_ind <- which(Srv_hat_list[[i]]$Fleet_code == srvs[srv])
      ymax[srv] <- max((c(ymax[srv], Srv_hat_list[[i]]$Residual[srv_ind])), na.rm = T)
      ymin[srv] <- min((c(ymin[srv], Srv_hat_list[[i]]$Residual[srv_ind])), na.rm = T)
    }
  }
  ymax <- ymax + top_adj * (ymax-ymin)

  positions=seq(-0.2, 0.2, length.out = length(Rceattle))

  # Assume colors if not provided
  if (is.null(line_col)) {
    line_col <- rev(oce::oce.colorsViridis(length(Rceattle)))
  }

  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  for (j in 1:loops) {

    # Plot/save each survey individually
    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(srv in 1:nsrv){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

        # Save
        if(j == 2){
          filename <- paste0(file, "fleet",srvs[j]," ",as.character(fleet_control$Fleet_name[srvs[srv]]), "_survey_index", ".png")
          png(file = filename, width = width, height = height, res = 200, units = "in")
        }

        minyr <- min(sapply(Srv_hat_list, function(x) min(x[which(x$Fleet_code == srvs[srv] & x$Year > 0),]$Year)))
        maxyr <- max(sapply(Srv_hat_list, function(x) max(x[which(x$Fleet_code == srvs[srv]),]$Year)))

        par(Par)
        plot(NA, NA, ylab="Index residual", xlab="Year", ylim = c((ymin[srv]), (ymax[srv])), xlim = c(minyr, maxyr + (maxyr - minyr) * right_adj), type='n', xaxt="n", yaxt="n")
        abline(h = 0, lty = 2)

        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Loop through models
        for (k in 1:length(Rceattle)) {

          # Subset data by fleet and model
          srv_tmp <- Srv_hat_list[[k]][which(Srv_hat_list[[k]]$Fleet_code == srvs[srv]),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_hat_list[[k]]$Fleet_code == srvs[srv]),]

          # Plot residual
          for(yr in 1:length(srv_hat_tmp$Year)){
            lines(rep((srv_hat_tmp$Year[yr] + positions[k]),2),c(0, srv_hat_tmp$Residual[yr]), col = line_col[k])
          }
          points(srv_hat_tmp$Year + positions[k], srv_hat_tmp$Residual, col=1, pch=21, bg = line_col[k])
        }

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[srvs[srv]]),bty="n",y.intersp = -0.2,cex=0.8)

        # Model names
        if(!is.null(model_names)){
          legend(
            "topright",
            legend = model_names,
            pch = rep(16, length(line_col)), cex=0.8,
            col = line_col,
            bty = "n"
          )

        }
        # Save plot
        if(j == 2){dev.off()}
      }
    }

    # Plot/save each survey individually
    if(single.plots==FALSE){

      # Set heights of plot
      if(is.null(width)) width = 7
      if(is.null(height)) height = ifelse(nsrv==1,5,ifelse(nsrv==2,3.,2.5))*round(nsrv/2+0.01,0)


      Par = list(mfrow=c(round(nsrv/2+0.01,0),ifelse(nsrv==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)

      # Save
      if(j == 2){
        filename <- paste0(file,"_survey_indices", ".png")
        png(file = filename, width = width, height = height, res = 200, units = "in")
      }
      par(Par)


      for(srv in 1:nsrv){
        minyr <- min(sapply(Srv_hat_list, function(x) min(x[which(x$Fleet_code == srvs[srv] & x$Year > 0),]$Year)))
        maxyr <- max(sapply(Srv_hat_list, function(x) max(x[which(x$Fleet_code == srvs[srv]),]$Year)))

        xlim <- c(minyr, maxyr)
        if(srv == 1){
          xlim <- c(minyr, maxyr + (maxyr - minyr) * right_adj)
        }

        plot(NA, NA, ylab="", xlab="", ylim = c((ymin[srv]), (ymax[srv])), xlim = xlim, type='n', xaxt="n", yaxt="n")
        abline(h = 0, lty = 2)
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)

        # Index name
        legend('topleft',as.character(fleet_control$Fleet_name[srvs[srv]]),bty="n",y.intersp = -0.2,cex=0.8)

        # Model names
        if(srv == 1){
          if(!is.null(model_names)){
            legend(
              "topright",
              legend = model_names,
              pch = rep(16, length(line_col)), cex=0.8,
              col = line_col,
              bty = "n"
            )
          }
        }

        # Loop through models
        for (k in 1:length(Rceattle)) {

          # Subset data by fleet and model
          srv_tmp <- Srv_hat_list[[k]][which(Srv_hat_list[[k]]$Fleet_code == srvs[srv]),]
          srv_hat_tmp <- Srv_hat_list[[k]][which(Srv_hat_list[[k]]$Fleet_code == srvs[srv]),]

          # Plot residual
          for(yr in 1:length(srv_hat_tmp$Year)){
            lines(rep((srv_hat_tmp$Year[yr] + positions[k]),2),c(0, srv_hat_tmp$Residual[yr]), col = line_col[k])
          }
          points(srv_hat_tmp$Year + positions[k], srv_hat_tmp$Residual, col=1, pch=21, bg = line_col[k])
        }
      }
    }
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext(paste("Index residual"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
    if(j == 2){dev.off()}
  }
} # End of plot

