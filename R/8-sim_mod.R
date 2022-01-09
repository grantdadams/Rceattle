#' Simulate Rceattle data
#'
#' @description  Simulates data used in Rceattle from the expected values etimated from Rceattle. The variances and uncertainty are the same as used in the operating model. The function currently simulates (assumed distribution) the following: survey biomass (log-normal), survey catch-at-length/age (multinomial), EIT biomass (log-normal), EIT catch-at-length/age (multinomial), total catch (kg) (log-normal), and catch-at-length/age.
#'
#' @param Rceattle CEATTLE model object exported from \code{\link{Rceattle}}
#' @param simulate TRUE/FALSE, whether to simulate the data or export the expected value
#'
#' @export
sim_mod <- function(Rceattle, simulate = FALSE) {
    # TODO Options for simulation diet data: multinomial, sqrt-normal, dirichlet, multinomial
    dat_sim <- Rceattle$data_list


    # Slot 0 -- BT survey biomass -- NFMS annual BT survey
    srv_biom_lse = Rceattle$quantities$srv_log_sd_hat

    if (simulate) {
        # Simulate
        values <- exp(rnorm(length(dat_sim$srv_biom$Observation), mean = log(Rceattle$quantities$srv_bio_hat), sd = srv_biom_lse))
    } else {
        # Estimated value
        values <- Rceattle$quantities$srv_bio_hat
    }
    dat_sim$srv_biom$Observation = values


    # Slot 1 -- Age composition
    for (obs in 1:nrow(dat_sim$comp_data)) {

        sp <- dat_sim$comp_data$Species[obs]

        # Age or length?
        if(dat_sim$comp_data$Age0_Length1[obs] == 0){
            age_length = Rceattle$data_list$nages[sp]
        }
        if(dat_sim$comp_data$Age0_Length1[obs] == 1){
            age_length = Rceattle$data_list$nlengths[sp]
        }

        if (simulate & (sum(Rceattle$quantities$comp_hat[obs,], na.rm = TRUE) > 0)) {
            # Simulate
            # Sex (0 = combined, 1 = female, 2 = male, 3 = joint)
            #FIXME add dirichlet multinomial

            # Combined, females, or males
            if(dat_sim$comp_data$Sex[obs] %in% c(0,1,2)){

                # If the model is all combined sex, this will throw an error
                if(Rceattle$data_list$nsex[sp] != 1){
                    Rceattle$quantities$comp_hat[obs, (age_length+1):(2*age_length)] <- 0 # Set males to 0
                }
                values <- rmultinom(n = 1, size = dat_sim$comp_data$Sample_size[obs], prob = Rceattle$quantities$comp_hat[obs,])
            }
            # # Males
            # if(dat_sim$comp_data$Sex[obs] %in% c(2)){
            #     Rceattle$quantities$comp_hat[obs, 1:age_length] <- 0 # Set females to 0
            #     values <- rmultinom(n = 1, size = dat_sim$comp_data$Sample_size[obs], prob = Rceattle$quantities$comp_hat[obs,])
            # }
            # Joint
            if(dat_sim$comp_data$Sex[obs] %in% c(3)){
                values <- rmultinom(n = 1, size = dat_sim$comp_data$Sample_size[obs], prob = Rceattle$quantities$comp_hat[obs,])
            }

        } else {
            # Expected value
            values <- Rceattle$quantities$comp_hat[obs, ]

        }

        dat_sim$comp_data[obs, 9:ncol(dat_sim$comp_data)] = values
    }




    # Slot 2 -- Total catch -- Fishery observer data
    fsh_biom_lse = Rceattle$quantities$fsh_log_sd_hat

    if (simulate) {
        # Simulate
        values <- exp(rnorm(length(dat_sim$fsh_biom$Catch), mean = log(Rceattle$quantities$fsh_bio_hat), sd = fsh_biom_lse))
    } else {
        # simulate value
        values <- Rceattle$quantities$fsh_bio_hat
    }

    dat_sim$fsh_biom$Catch = values


    #TODO
    # # Slot 5 -- Diet composition from lognormal suitability 4D
    # if (length(dim(dat_sim$UobsWtAge)) == 4) {
    #     for (sp in 1:dat_sim$nspp) {
    #         for (r_age in 1:dat_sim$nages[sp]) {
    #             if (Rceattle$data_list$suitMode > 0 & simulate & sum(Rceattle$quantities$mn_UobsWtAge_hat[sp, , r_age,
    #                                                                                                       ] > 0) > 0) {
    #                 values <- rmultinom(n = 1, size = dat_sim$stom_tau[sp], prob = Rceattle$quantities$mn_UobsWtAge_hat[sp,
    #                                                                                                                     , r_age, ])  #FIXME change sample size
    #             } else {
    #                 values <- Rceattle$quantities$mn_UobsWtAge_hat[sp, , r_age, ]
    #             }
    #             dat_sim$UobsWtAge[sp, , r_age, ] <- replace(dat_sim$UobsWtAge[sp, , r_age, ], values = values)
    #         }
    #     }
    # }
    #
    # # 5D
    # if (length(dim(dat_sim$UobsWtAge)) == 5) {
    #     for (sp in 1:dat_sim$nspp) {
    #         for (yr in 1:dat_sim$nyrs_fsh_comp[sp]) {
    #             for (r_age in 1:dat_sim$nages[sp]) {
    #                 if (Rceattle$data_list$suitMode > 0 & simulate & sum(Rceattle$quantities$UobsWtAge_hat[sp, , r_age,
    #                                                                                                        , yr] > 0) > 0) {
    #                     values <- rmultinom(n = 1, size = dat_sim$stom_tau[sp], prob = Rceattle$quantities$UobsWtAge_hat[sp,
    #                                                                                                                      , r_age, , yr])  #FIXME change sample size
    #                 } else {
    #                     values <- Rceattle$quantities$UobsWtAge_hat[sp, , r_age, , yr]
    #                 }
    #                 dat_sim$UobsWtAge[sp, , r_age, , yr] <- replace(dat_sim$UobsWtAge[sp, , r_age, , yr], values = values)
    #             }
    #         }
    #     }
    # }

    return(dat_sim)
}



#' Evaluate simulation performance
#'
#' @description Function to evaluate the simulation performance with regard to bias using the median relative error (MRE) and precision using the coefficient of variation.
#'
#' @param operating_mod CEATTLE model object exported from \code{\link{Rceattle}} to be used as the operating model
#' @param simulation_mods List of CEATTLE model objects exported from \code{\link{Rceattle}} fit to simulated data
compare_sim <- function(operating_mod, simulation_mods, object = "quantities") {
    # TODO update

    # Get differences
    sim_mre <- list()
    sim_mse <- list()
    sim_mean <- list()
    sim_median <- list()
    sim_sd <- list()
    sim_cv <- list()
    sim_params <- list()

    for (j in 1:length(names(operating_mod[[object]]))) {
        param <- names(operating_mod[[object]])[j]

        sim_mre[[param]] <- list()
        sim_mse[[param]] <- list()
        sim_mean[[param]] <- list()
        sim_sd[[param]] <- list()
        sim_cv[[param]] <- list()
        sim_params[[param]] <- list()

        om_params <- operating_mod[[object]][[param]]

        for (i in 1:length(simulation_mods)) {

            sm_params <- simulation_mods[[i]][[object]][[param]]

            sim_params[[param]][[i]] <- sm_params
            sim_mre[[param]][[i]] <- (sm_params - om_params)/om_params
            sim_mse[[param]][[i]] <- (sm_params - om_params)^2
        }

        param_dim <- length(dim(om_params))

        # If 1 value
        if (param_dim == 0) {
            sim_mean[[param]] <- mean(unlist(sim_params[[param]]))
            sim_sd[[param]] <- sd(unlist(sim_params[[param]]))
            sim_cv[[param]] <- sim_sd[[param]]/sim_mean[[param]]

            sim_mre[[param]] <- median(unlist(sim_mre[[param]]))
            sim_mse[[param]] <- mean(unlist(sim_mse[[param]]))
        }

        # If multiple values
        if (param_dim > 0) {
            # Get mean, sd, and CV
            sim_mean[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, mean)
            sim_median[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, median)
            sim_sd[[param]] <- apply(simplify2array(sim_params[[param]]), 1:param_dim, sd)
            sim_cv[[param]] <- sim_sd[[param]]/sim_mean[[param]]

            sim_mre[[param]] <- apply(simplify2array(sim_mre[[param]]), 1:param_dim, median)
            sim_mse[[param]] <- apply(simplify2array(sim_mse[[param]]), 1:param_dim, mean)
        }
    }


    result_list <- list(Mean = sim_mean, Median = sim_median, SD = sim_sd, CV = sim_cv, MRE = sim_mre, MSE = sim_mse)
    return(result_list)
}

