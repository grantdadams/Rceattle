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
    control <- matrix(NA, ncol = data_list$nspp, nrow = 12)
    control[1, 1] <- data_list$nspp
    control[2, 1] <- data_list$styr
    control[3, 1] <- data_list$endyr
    control[4, 1] <- data_list$projyr
    control[5, ] <- data_list$nages
    control[6, ] <- data_list$minage
    control[7, ] <- data_list$nlengths
    control[8, ] <- data_list$pop_wt_index
    control[9, ] <- data_list$pop_alk_index
    control[10, ] <- data_list$sigma_rec_prior
    control[11, ] <- data_list$other_food
    control[12, ] <- data_list$stom_tau
    control <- as.data.frame(control)
    control <- cbind(c("nspp", "styr", "endyr", "projyr", "nages", "minage", "nlengths", "pop_wt_index", "pop_alk_index", "sigma_rec_prior",
        "other_food", "stom_sample_size"), control)
    colnames(control) <- c("Object", data_list$spnames)
    names_used <- c(names_used, as.character(control$Object))

    xcel_list$control <- control


    # Jnll
    xcel_list$jnll_comp <- as.data.frame(quantities$jnll_comp)
    xcel_list$jnll_comp <- cbind(data.frame(JNLL_component = rownames(quantities$jnll_comp)), xcel_list$jnll_comp)

    ############################################################ Survey data srv control
    srv_bits <- c("srv_control")
    for (i in 1:length(srv_bits)) {
        xcel_list[[srv_bits[i]]] <- data_list[[srv_bits[i]]]
    }
    names_used <- c(names_used, srv_bits)


    # srv_biom_hat
    xcel_list$srv_biom_hat <- data_list$srv_biom
    xcel_list$srv_biom_hat$CV <- NULL
    xcel_list$srv_biom_hat$Estimated_index <- quantities$srv_bio_hat
    xcel_list$srv_biom_hat$CV <- quantities$srv_cv_hat


    # srv_comp_hat
    xcel_list$srv_comp_hat <- data_list$srv_comp
    xcel_list$srv_comp_hat[, grep("Comp_", colnames(data_list$srv_comp))] <- quantities$srv_comp_hat
    colnames(xcel_list$srv_comp_hat)[grep("Comp_", colnames(xcel_list$srv_comp_hat))] <- paste0("Comp_hat_", 1:length(grep("Comp_",
        colnames(xcel_list$srv_comp_hat))))

    ############################################################ Fisheries data fsh control
    fsh_bits <- c("fsh_control")
    for (i in 1:length(fsh_bits)) {
        xcel_list[[fsh_bits[i]]] <- data_list[[fsh_bits[i]]]
    }
    names_used <- c(names_used, fsh_bits)

    # fsh_biom_hat
    xcel_list$fsh_biom_hat <- data_list$fsh_biom
    xcel_list$fsh_biom_hat$CV <- NULL
    xcel_list$fsh_biom_hat$Estimated_catch <- quantities$fsh_bio_hat
    xcel_list$fsh_biom_hat$CV <- quantities$fsh_cv_hat

    # fsh_comp_hat
    xcel_list$fsh_comp_hat <- data_list$fsh_comp
    xcel_list$fsh_comp_hat[, grep("Comp_", colnames(data_list$fsh_comp))] <- quantities$fsh_comp_hat
    colnames(xcel_list$fsh_comp_hat)[grep("Comp_", colnames(xcel_list$fsh_comp_hat))] <- paste0("Comp_hat_", 1:length(grep("Comp_",
        colnames(xcel_list$fsh_comp_hat))))

    ############################################################ Population data
    yrs <- data_list$styr:data_list$endyr
    nyrs <- length(yrs)

    # Nbyage
    xcel_list$NByage <- t(quantities$NByage[1, , 1:nyrs])
    for (i in 2:data_list$nspp) {
        xcel_list$NByage <- rbind(xcel_list$NByage, t(quantities$NByage[i, , 1:nyrs]))
    }
    xcel_list$NByage <- as.data.frame(xcel_list$NByage)
    colnames(xcel_list$NByage) <- paste0("Age_", 1:ncol(xcel_list$NByage))
    xcel_list$NByage <- cbind(data.frame(Year = rep(yrs, 3)), xcel_list$NByage)
    xcel_list$NByage <- cbind(data.frame(Species = rep(1:data_list$nspp, each = nyrs)), xcel_list$NByage)
    xcel_list$NByage <- cbind(data.frame(Species_name = rep(data_list$spnames, each = nyrs)), xcel_list$NByage)

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
    mortality <- matrix(NA, nrow = nyrs, ncol = data_list$nspp * 3 + 1)
    mortality[,1] <- yrs
    ind <- 2
    for(sp in 1:data_list$nspp){
      mortality[,ind] <- quantities$M2[sp,1,1:nyrs]; ind = ind + 1
      mortality[,ind] <- quantities$F_tot[sp,1,1:nyrs]; ind = ind + 1
      mortality[,ind] <- quantities$Zed[sp,1,1:nyrs]; ind = ind + 1
    }

    colnames(mortality) <- c("Year", paste0(rep(c("M2", "Ftot", "Z"), 3), rep(1:data_list$nspp, each = 3)))
    xcel_list$Age1Mortality <- as.data.frame(mortality)

    ############################################################ write the data
    writexl::write_xlsx(xcel_list, file)
}
