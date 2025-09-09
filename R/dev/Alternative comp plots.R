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
    fleet_control <- Rceattle$data_list$fleet_control

    # Get observed
    comp_data <- Rceattle$data_list$comp_data

    # Normalize
    comp_data[,grep("Comp_", colnames(comp_data))] = comp_data[,grep("Comp_", colnames(comp_data))]/rowSums(comp_data[,grep("Comp_", colnames(comp_data))], na.rm = TRUE)

    comp_data <- comp_data %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("Comp_"),
                          names_to = "Comp",
                          names_prefix = "Comp_",
                          values_to = "Obs")

    # Get estimated
    comp_hat <- Rceattle$data_list$comp_data
    comp_hat[,grep("Comp_", colnames(comp_hat))] = Rceattle$quantities$comp_hat
    comp_hat <- comp_hat %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("Comp_"),
                          names_to = "Comp",
                          names_prefix = "Comp_",
                          values_to = "Est")


    # * Calculate pearson ----
    comp_data <- comp_data %>%
      dplyr::full_join(comp_hat) %>%
      dplyr::mutate(Pearson = (Obs - Est)/ sqrt( ( Est * (1 - Est)) / Sample_size))


    # * Distinct data sets ----
    # - age/length and fleets
    data_sets <- comp_data %>%
      distinct(Fleet_name, Fleet_code, Species, Sex, Age0_Length1)


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Loop around data-sets ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    for (j in 1:nrow(data_sets)) {


      # Species and bins
      sp <- data_sets$Species[j]
      nsex <- Rceattle$data_list$nsex[sp]
      comp_type = data_sets$Age0_Length1[j] + 1 # Age or length
      nbin <- c(Rceattle$data_list$nages[sp], Rceattle$data_list$nlengths[sp])[comp_type]

      # Extract comps
      comp_tmp <- comp_data %>%
        dplyr::inner_join(data_sets[j,]) %>%
        dplyr::filter(!is.na(Obs)) %>%
        dplyr::mutate(Comp = as.numeric(Comp)) %>%
        dplyr::mutate(
          Sex = dplyr::case_when(
            Sex == 0 ~ "Sex combined",
            Sex == 1 ~ "Females",
            Sex == 2 ~ "Males",
            Sex == 3 & Comp <= nbin ~ "Females",
            Sex == 3 & Comp > nbin ~ "Males",
            .default = NA
          ),
          CompNeg = ifelse(Sex == "Males" & Comp > nbin, -1 * (Comp - nbin), Comp)
          Comp = ifelse(Sex == "Males" & Comp > nbin, (Comp - nbin), Comp)) # Adjust males if joint sex


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # Plot type 1 - Pearson residual ----
      plot_pear <- comp_tmp %>%
        dplyr::filter(Obs > 0) %>%
        ggplot2::ggplot(ggplot2::aes(x = Year, y = Comp, size = abs(pearson), color = pearson < 0, group = Sex)) +
        ggplot2::facet_grid(~Sex) +
        ggplot2::geom_point(alpha = 0.8) +
        # scale_size(range = c(range_comp[1], range_comp[2]), name="Population (M)") +
        ggplot2::theme_classic() +
        ggplot2::ylim(range(comp_tmp$Comp)) +
        ggplot2::ylab(c("Age", "Length")[comp_type]) +
        ggplot2::xlim(range(comp_tmp$Year)) +
        ggplot2::xlab("Year") +
        ggplot2::ggtitle(c("Pearson residual: Age comp", "Pearson residual: Length comp")[comp_type])

      print(plot_pear)

      #  - Save
      if (!is.null(file)) {
        ggplot2:ggsave(
          filename = paste0(file, data_sets$Fleet_name[j],c("_age-comps", "_length-comps")[comp_type], "_pearson_residual.png"),
          plot_pear,
          width = 10,
          height = 7,
          units = "in",
          res = 300
        )
      }


      # Plot type 2 - Annual histograms ----
      annual_comp <- comp_tmp %>%
        ggplot(aes(x = Comp, y = Obs, fill = Sex)) +
        facet_wrap(~Year+Sex, labeller =
                     labeller(
                       Sex = label_value,
                       Year = label_value,
                       .multi_line = FALSE
                     )) +
        geom_col() +
        ggplot2::xlab(c("Age", "Length")[comp_type]) +
        theme_classic() +
        theme(legend.position="none") +
        theme(axis.title.y=element_blank()) +
        geom_line(aes(
          x = Comp,
          y = Est
        ))

      print(annual_comp)

      #  - Save
      if (!is.null(file)) {
        filename <- paste0(file, paste0(c("_age_comps_histograms_", "_length_comps_histograms_"), "fleet_code_",srv[j] )[comp_type + 1], ".png")
        ggplot2:ggsave(
          filename = filename,
          annual_comp
          width = 10,
          height = 7,
          units = "in",
          res = 300
        )
      }

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # Plot type 3 - Histograms of aggregated comps ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#



      filename <- paste0(file, paste0(c("_aggregated_age_comps_histograms_", "_aggregated_length_comps_histograms_"), "fleet_code_",srv[j] )[comp_type + 1], ".png")
    }





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







#' Plot diet composition fits
#'
#' @description
#' If year == 0, diet data are averaged from suit_styr to suit_endyr
#' If prey_age >= 0 diet data are diet proportion of prey-at-age in predator-at-age
#' If prey_age < 0 diet data are diet proportion of prey-spp in predator-at-age (sum across prey ages)
#' If prey_age < 0 and pred_age < 0, diet data are mean diet proportion of prey-spp in predator-spp (sum across prey ages and take mean across predator ages)
#' If prey_age < 0 and pred_age < -500, diet data are weighted mean diet proportion of prey-spp in predator-spp (sum across prey ages and take weighted mean across predator ages)
#'
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param species Species names for legend
#'
#' @return Returns and saves a figure
#' @export
plot_diet_comp <-
  function(Rceattle,
           file = NULL,
           species = NULL) {

    # Make sure we are using only one model
    if(class(Rceattle) != "Rceattle"){
      stop("Please only use one Rceattle model")
    }
    data_list <- Rceattle$data_list

    # Species names
    if(is.null(species)){
      species =  Rceattle$data_list$spnames
    }

    # Get colors
    colvec=c("red", "blue", "black")

    # * Extract data objects ----
    # - Get observed
    comp_data <- Rceattle$data_list$diet_data
    # - Get estimated
    comp_data$Est = Rceattle$quantities$diet_hat[,2]

    comp_data <- comp_data %>%
      dplyr::mutate(pearson = (Stomach_proportion_by_weight - Est)/ sqrt( ( Est * (1 - Est)) / Sample_size))

    # If year == 0, diet data are averaged from suit_styr to suit_endyr
    # If prey_age >= 0 diet data are diet proportion of prey-at-age in predator-at-age
    # If prey_age < 0 diet data are diet proportion of prey-spp in predator-at-age (sum across prey ages)
    # If prey_age < 0 and pred_age < 0, diet data are mean diet proportion of prey-spp in predator-spp (sum across prey ages and take mean across predator ages)
    # If prey_age < 0 and pred_age < -500, diet data are weighted mean diet proportion of prey-spp in predator-spp (sum across prey ages and take weighted mean across predator ages)


    # Loop around predators ----
    for(pred in 1:data_list$nspp) {
      for(pred_sex in 1:data_list$nsex[pred]){
        for(prey in 1:data_list$nspp) {

          # * Get sex for legend ----
          if(data_list$nsex[pred] > 1){
            pred_legend <- paste("Pred-", species[pred], ifelse(sex == 1, "female", "male"))
            pred_sex = pred_sex - 1
          } else{
            pred_legend <- paste("Pred-", species[pred])
            pred_sex = 0
          }

          # * Extract comps ----
          comp_tmp <- comp_data %>%
            dplyr::filter(Pred == pred & Pred_sex == pred_sex & Prey == prey)

          # - Years
          yrs <- sort(unique(comp_tmp$Year))
          nyrs <- length(yrs)

          # - Min and max
          range_comp <- range(c(comp_tmp$Stomach_proportion_by_weight, comp_tmp$Est))
          range_pearson <- range(comp_tmp$pearson)

          # * Plot annual comps ----
          for(yr in 1:nyrs){

            # Subset year for observed and predicted comp
            comp_tmp_yr <- comp_tmp %>%
              dplyr::filter(Year == yrs[yr] ) %>%
              dplyr::mutate(Prey_age = ifelse(Prey_sex == 2, -Prey_age, Prey_age))


            plot_obs <- comp_tmp_yr %>%
              dplyr::filter(Stomach_proportion_by_weight > 0) %>%
              ggplot2::ggplot(ggplot2::aes(x = Pred_age, y = Prey_age, size = Stomach_proportion_by_weight)) +
              ggplot2::geom_point(alpha = 1) +
              # scale_size(range = c(range_comp[1], range_comp[2]), name="Population (M)") +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Observed diet: year", yrs[yr])) +
              theme(legend.position = c(0.25, 0.7),
                    legend.title = element_blank())

            plot_est <- comp_tmp_yr %>%
              dplyr::filter(Stomach_proportion_by_weight > 0) %>%
              ggplot2::ggplot(ggplot2::aes(x = Pred_age, y = Prey_age, size = Est)) +
              ggplot2::geom_point(alpha = 1) +
              # scale_size(range = c(range_comp[1], range_comp[2]), name="Population (M)") +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Estimated diet: year", yrs[yr])) +
              theme(legend.position = "none")

            plot_pear <- comp_tmp_yr %>%
              dplyr::filter(Stomach_proportion_by_weight > 0) %>%
              ggplot2::ggplot(ggplot2::aes(x = Pred_age, y = Prey_age, size = abs(pearson), color = pearson < 0)) +
              ggplot2::geom_point(alpha = 1) +
              # scale_size(range = c(range_comp[1], range_comp[2]), name="Population (M)") +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Pearson residual: year", yrs[yr])) +
              theme(legend.position = c(0.25, 0.8))

            p1 <- cowplot::plot_grid(plot_obs, plot_est, plot_pear, nrow = 1)
            print(p1)

            # Save ----
            if (!is.null(file)) {
              filename <- paste0(file, "_aggregated_diet_comps_histograms_year", yr,"_", pred_legend, "prey", species[prey],".png")
              ggplot2::ggsave(filename = filename,
                              plot = p1,
                              width = 10,
                              height = 6.5,
                              units = "in",
                              dpi = 300
              )
            }
          }
        }
      }
    }

    # End
  }

