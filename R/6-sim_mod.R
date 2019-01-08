# Options for simulation diet data: multinomial, sqrt-normal, dirichlet, multinomial
sim_mod <- function( ceattle_obj ){

  dat_sim <- ceattle_obj$data_list

  for(sp in 1:ceattle_obj$data_list$nspp){ #FIXME: may need to adjust for multiple surveys

    # Survey biomass
    srv_biom_lse = dat_sim$srv_biom_se / dat_sim$srv_biom
    srv_biom_lse = sqrt( log( (  srv_biom_lse^2 ) + 1));

    dat_sim$srv_biom = replace(dat_sim$srv_biom, values = exp(rnorm(length(dat_sim$srv_biom), log(ceattle_obj$quantities$srv_bio_hat[,1:ncol(dat_sim$srv_biom)]), srv_biom_lse)))


    #

  }

}
