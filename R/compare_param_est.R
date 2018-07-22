

compare_parameter_estimation <- function(opt, params){

  (opt$par - unlist(params))/unlist(params)

}
