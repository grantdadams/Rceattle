#' Write excel file of Rceattle results
#'
#' @param Rceattle Single or list of Rceattle model objects exported from \code{\link{Rceattle}}
#' @param file name of a file to identified the files exported by the
#'   function.
#'
#' @export
#'
#' @examples
#'
#'# Load package and data
#'library(Rceattle)
#'data(BS2017SS) # ?BS2017SS for more information on the data
#'
#'# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
#'ss_run <- fit_mod(data_list = BS2017SS,
#'    inits = NULL, # Initial parameters = 0
#'    file = NULL, # Don't save
#'    debug = 0, # Estimate
#'    random_rec = FALSE, # No random recruitment
#'    msmMode = 0, # Single species mode
#'    avgnMode = 0,
#'    silent = TRUE)
#'
#' write_results(Rceattle = BS2017SS, file = 'Rceattle_results.xlsx')
write_results <- function(Rceattle, file = "Rceattle_results.xlsx") {
    "%!in%" <- function(x, y) !(x %in% y)


    # Make sure we are using only one model
    if (class(Rceattle) != "Rceattle") {
        stop("Please only use one Rceattle model")
    }

    ## setup a workbook
    data_names <- names(Rceattle$quantities)
    names_used <- c()

    xcel_list <- list()


    # Save control specifications
    data_list <- Rceattle$data_list
    quantities <- Rceattle$quantities


    # control
    control <- matrix(NA, ncol = data_list$nspp, nrow = 15)
    control[1, 1] <- data_list$nspp
    control[2, 1] <- data_list$styr
    control[3, 1] <- data_list$endyr
    control[4, 1] <- data_list$projyr
    control[5, ] <- data_list$nsex
    control[6, ] <- data_list$spawn_month
    control[7, ] <- data_list$R_sexr
    control[8, ] <- data_list$nages
    control[9, ] <- data_list$minage
    control[10, ] <- data_list$nlengths
    control[11, ] <- data_list$pop_wt_index
    control[12, ] <- data_list$ssb_wt_index
    control[13, ] <- data_list$pop_alk_index
    control[14, ] <- data_list$sigma_rec_prior
    control[15, ] <- data_list$other_food
    control <- as.data.frame(control)
    control <- cbind(c("nspp", "styr", "endyr", "projyr", "nsex", "spawn_month", "R_sexr", "nages", "minage", "nlengths", "pop_wt_index", "ssb_wt_index","pop_alk_index", "sigma_rec_prior",
                       "other_food"), control)
    colnames(control) <- c("Object", data_list$spnames)
    names_used <- c(names_used, as.character(control$Object))

    xcel_list$control <- control


    # Jnll
    xcel_list$jnll_comp <- as.data.frame(quantities$jnll_comp)
    xcel_list$jnll_comp <- cbind(data.frame(JNLL_component = rownames(quantities$jnll_comp)), xcel_list$jnll_comp)

    ############################################################ Survey data srv control
    srv_bits <- c("fleet_control")
    for (i in 1:length(srv_bits)) {
        xcel_list[[srv_bits[i]]] <- data_list[[srv_bits[i]]]
    }
    names_used <- c(names_used, srv_bits)


    # srv_biom_hat
    xcel_list$srv_biom_hat <- data_list$srv_biom
    xcel_list$srv_biom_hat$CV <- NULL
    xcel_list$srv_biom_hat$Estimated_index <- quantities$srv_bio_hat
    xcel_list$srv_biom_hat$CV <- quantities$srv_cv_hat


    # fsh_biom_hat
    xcel_list$fsh_biom_hat <- data_list$fsh_biom
    xcel_list$fsh_biom_hat$CV <- NULL
    xcel_list$fsh_biom_hat$Estimated_catch <- quantities$fsh_bio_hat
    xcel_list$fsh_biom_hat$CV <- quantities$fsh_cv_hat


    # srv_comp_hat
    xcel_list$comp_data_hat <- data_list$comp_data
    xcel_list$comp_data_hat[, grep("Comp_", colnames(data_list$comp_data))] <- quantities$comp_hat
    colnames(xcel_list$comp_data_hat)[grep("Comp_", colnames(xcel_list$comp_data_hat))] <- paste0("Comp_hat_", 1:length(grep("Comp_",
        colnames(xcel_list$comp_data_hat))))


    ############################################################ Population data
    yrs <- data_list$styr:data_list$endyr
    nyrs <- length(yrs)

    # Nbyage
    xcel_list$NByage <- as.data.frame(t(quantities$NByage[1, 1, , 1:nyrs]))
    xcel_list$NByage$Species_name <- data_list$spnames[1]
    xcel_list$NByage$Species = 1
    xcel_list$NByage$Sex = ifelse(data_list$nsex[1] == 1, 0, 1)
    xcel_list$NByage$Year = yrs

    if(data_list$nsex[1] == 2){
      nbyage_tmp <-  as.data.frame(t(quantities$NByage[1, 2, , 1:nyrs] ))
      nbyage_tmp$Species_name <- data_list$spnames[1]
      nbyage_tmp$Species = 1
      nbyage_tmp$Sex = 2
      nbyage_tmp$Year = yrs


      xcel_list$NByage <- rbind(xcel_list$NByage, nbyage_tmp)
    }

    for (i in 2:data_list$nspp) {
      for(sex in 1:data_list$nsex[i]){

        nbyage_tmp <-  as.data.frame(t(quantities$NByage[i, sex, , 1:nyrs] ))
        nbyage_tmp$Species_name <- data_list$spnames[i]
        nbyage_tmp$Species = i
        nbyage_tmp$Sex = ifelse(data_list$nsex[i] == 1, 0, sex)
        nbyage_tmp$Year = yrs


        xcel_list$NByage <- rbind(xcel_list$NByage, nbyage_tmp)
      }
    }
    info_cols <- which(colnames(xcel_list$NByage) %in% c("Year", "Species", "Sex", "Species_name"))
    nbyage_info <- xcel_list$NByage[,info_cols]
    nbyage <- xcel_list$NByage[,-info_cols]
    colnames(nbyage) <- paste0("Age_", 1:ncol(nbyage))
    xcel_list$NByage <- cbind(nbyage_info, nbyage)

    # Recruitment
    xcel_list$Recruitment <- as.data.frame(t(quantities$R[, 1:nyrs]))
    colnames(xcel_list$Recruitment) <- data_list$spnames
    xcel_list$Recruitment <- cbind(data.frame(Year = yrs), xcel_list$Recruitment)

    # Biomass
    xcel_list$Biomass <- as.data.frame(t(quantities$biomass[, 1:nyrs]))
    colnames(xcel_list$Biomass) <- data_list$spnames
    xcel_list$Biomass <- cbind(data.frame(Year = yrs), xcel_list$Biomass)

    # SSB
    xcel_list$SSB <- as.data.frame(t(quantities$biomassSSB[, 1:nyrs]))
    colnames(xcel_list$SSB) <- data_list$spnames
    xcel_list$SSB <- cbind(data.frame(Year = yrs), xcel_list$SSB)

    # Age-1 Mortality
    mortality <- matrix(NA, nrow = nyrs, ncol = sum(data_list$nsex) * 3 + 1)
    mortality[,1] <- yrs
    ind <- 2
    colnames_sex <- c()
    colnames_sp <- c()
    for(sp in 1:data_list$nspp){
      for(sex in 1:data_list$nsex[sp]){
      mortality[,ind] <- quantities$M2[sp, sex ,1,1:nyrs]; ind = ind + 1
      mortality[,ind] <- quantities$F_tot[sp, sex,1,1:nyrs]; ind = ind + 1
      mortality[,ind] <- quantities$Zed[sp, sex,1,1:nyrs]; ind = ind + 1
      colnames_sex <- c(colnames_sex, ifelse(data_list$nsex[sp] == 0, 0, sex))
      colnames_sp <- c(colnames_sp, sp)
      }
    }

    colnames(mortality) <- c("Year", paste0(rep(c("M2", "Ftot", "Z"), length(colnames_sex)), paste0("Sp", colnames_sp), paste0("_Sex",colnames_sex)))
    xcel_list$Age1Mortality <- as.data.frame(mortality)

    ############################################################ write the data
    writexl::write_xlsx(xcel_list, file)
}
