
#' Plot M1 + M2
#'
#' @description Function the plots the M1 and M2 as estimated from Rceattle
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param incl_proj Include the projection years (TRUE/FALSE)
#' @param zlim zlim for M1 + M2 plots. Character - use max range across species in model. NULL - use species specific ranges. Vector of two.
#' @param contour If plot it to be done as contours rather than tiles.
#' @param species Species names for legend
#' @param width Plot width when saved "inches"
#' @param height Plot height when saved "inches"
#'
#' @export
plot_mortality <-
  function(Rceattle,
           file = NULL,
           incl_proj = FALSE,
           zlim = NULL,
           contour = FALSE,
           width = NULL,
           height = NULL) {

    # Convert single one into a list
    if(class(Rceattle) == "Rceattle"){
      Rceattle <- list(Rceattle)
    }

    if(length(Rceattle) > 1){
      stop("Can only plot one model")
    }

    # Extract data objects
    Years <- list()
    for(i in 1:length(Rceattle)){
      Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$endyr
      if(incl_proj){
        Years[[i]] <- Rceattle[[i]]$data_list$styr:Rceattle[[i]]$data_list$projyr
      }
    }
    nyrs_vec <- sapply(Years, length)
    nyrs <- max(nyrs_vec)
    maxyr <- max((sapply(Years, max)))
    minyr <- min((sapply(Years, min)))

    nspp <- Rceattle[[1]]$data_list$nspp
    estdynamics <- Rceattle[[1]]$data_list$estDynamics
    nages <- Rceattle[[1]]$data_list$nages
    minage <- Rceattle[[1]]$data_list$minage
    nsex <- Rceattle[[1]]$data_list$nsex

    # Get biomass
    m_array <-
      array(NA, dim = c(nspp, 2, max(nages), nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      m_array[, , , ,i] <- Rceattle[[i]]$quantities$M[,,,1:nyrs]
    }

    # Plot limits
    zmax <- c()
    zmin <- c()
    for (i in 1:dim(m_array)[1]) {
      zmax[i] <- max(c(m_array[i,,,,], 0), na.rm = T)
      zmin[i] <- min(c(m_array[i,,,,], 0), na.rm = T)
    }


    # Plot trajectory
    loops <- ifelse(is.null(file), 1, 2)

    #################################
    # Mortality time series
    #################################


    # Species
    for(j in 1:nspp){
      sp <- j

      if(estdynamics[j] == 0){

        # Sexes
        for(sex in 1:nsex[sp]){

          # Get sex for legend
          legend_sex = sex
          legend_sex2 = ifelse(sex == 1, "Female", "Male")
          if(nsex[sp] == 1){
            legend_sex <- 0
            legend_sex2 = "Combined"
          }

          # Save
          for (i in 1:loops) {
            if (i == 2) {
              filename <- paste0(file, "predation_and_residual_mortality_spp_",sp,"_sex_",legend_sex2,".png")
              png(
                file = filename ,
                width = ifelse(is.null(width), 8, width),
                height = ifelse(is.null(height), 5.5, height),

                units = "in",
                res = 300
              )
            }

            # Subset mortality data
            m_subset <- (m_array[j, sex, (1:nages[sp]), 1:nyrs, 1])

            # Get ages
            ages <- (1:(nages[sp])) - 1 + minage[sp]

            # Rearrange data
            data <- data.frame(Year = rep(Years[[1]], each = length(ages)), Age = rep(ages, length(Years[[1]])), M = c(m_subset))

            # Plot limits
            if(is.null(zlim)){
              zlim <- c(zmin[sp], zmax[sp])
            }

            if(is.character(zlim)){
              zlim <- c(min(zmin), max(zmax))
            }

            # Plot as contours
            if(contour){
              print(ggplot(data, aes(y = Age, x = Year, z = M, zmin = zlim[1], zmax = zlim[2])) + geom_contour(colour = 1, size = 0.5) + geom_contour_filled()  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) +  scale_x_continuous(expand = c(0, 0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size = 6)) + scale_fill_viridis_d("M1 + M2"))
            }

            # Plot as tiles
            if(contour == FALSE){
              print(ggplot(data, aes(y = Age, x = Year, zmin = zlim[1], zmax = zlim[2])) + geom_tile(aes(fill = M))  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) + coord_equal() +  scale_x_continuous(expand = c(0, 0))+ scale_fill_viridis_c("M1 + M2") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size = 6)))
            }
            if (i == 2) {
              dev.off()
            }
          }
        }
      }
    }
  }
