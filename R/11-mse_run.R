


mse_run <- function(operating_model, estimation_model, nsim = 10, assessment_period = 1){

  Rceattle_list <- list()

  for(sim in 1:nsim){

    # Set OM rec-devs

    # Years for simulations
    proj_yrs <- (estimation_model$data_list$endyr + 1) : estimation_model$data_list$projyr
    proj_nyrs <- length(proj_yrs)
    assess_yrs <- seq(from = estimation_model$data_list$endyr + 1, to = estimation_model$data_list$projyr,  by = assessment_period)

    # Run through model
    for(yr in 1:assess_yrs){

      # Get first year of catch form EM
      new_data <- estimation_model$data_list
    }



  }


}
