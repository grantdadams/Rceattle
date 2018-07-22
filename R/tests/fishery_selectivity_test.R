

# for (int sp_age=1;sp_age<nages(sp);sp_age++)
#   if (fsh_sel(sp,sp_age)>fsh_sel(sp,sp_age+1))
#     fsh_sel_like(sp) += 20.*square(log(fsh_sel(sp,sp_age)/fsh_sel(sp,sp_age+1)));
#     fsh_sel_like(sp) += curv_pen_fsh(sp)*norm2(first_difference( first_difference(log(fsh_sel(sp)))));
#
species = 2

norm2<-function(v){
  nn<-length(v)
  nu<-rep(0,nn)
  for(i in 1:nn){
    nu[i]<-(v[i])^2
  }
  return(sum(nu))
}

first_difference<-function(v){
  nn<-length(v)
  return( v[2:(nn)]-v[1:(nn-1)])
}

curv_pen_fsh <- 12.5

fsh_sel <- rep$fsh_sel[species,]
#fsh_sel <- tmp[[paste0("fsh_sel_", species)]]

fsh_sel_like <- 0



for(i in 1:(data_list$nages[species]-1)){
  if(fsh_sel[i] > fsh_sel[i+1]){
    fsh_sel_like <- fsh_sel_like + 20 * (log(fsh_sel[i]/fsh_sel[i+1]) ^ 2)
  }
}
  # fsh_sel_like <- fsh_sel_like + curv_pen_fsh * norm2(first_difference( first_difference(log(fsh_sel[1:data_list$nages[species]]))))

  fsh_sel_like
  rep$jnll_comp[7,]

  tmp[[paste0("fsh_sel_like_" , species)]]
