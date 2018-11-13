# # FUNCTION CALC_SELECTIVITY
# #   for (int sp=1;sp<=nspp;sp++)
# #   {
# #     // Fishery...
# #     fsh_sel(sp)(1,nselages)          = fsh_sel_coff(sp);                          // fsh_sel_coff(sp) parameters
# #     fsh_sel(sp)(nselages+1,nages(sp)) = fsh_sel(sp,nselages);
# #     avgsel_fsh(sp)                   = log(mean(mfexp(fsh_sel_coff(sp))));
# #     fsh_sel(sp)                     -= log(mean(mfexp(fsh_sel(sp))));
# #     fsh_sel(sp)                      = mfexp(fsh_sel(sp));
# #
# #     // Surveys...
# #     // Do this only if logistic is selected for this species...
# #     if (logist_sel_phase(sp) >0){// logist_sel_phase(sp)= -2 so this is not done
# #       for (int sp_age=1;sp_age<=nages(sp);sp_age++)
# #         srv_sel(sp,sp_age) = 1/(1+mfexp(-srv_sel_slp(sp)*(double(sp_age)-srv_sel_inf(sp))));
# #     }
# #     else // use coefficient form... these are the actual calculations
# #     {
# #       srv_sel(sp)(1,nselages)          = srv_sel_coff(sp);
# #       srv_sel(sp)(nselages+1,nages(sp)) = srv_sel(sp,nselages);
# #       avgsel_srv(sp)                   = log(mean(mfexp(srv_sel_coff(sp))));
# #       srv_sel(sp)                     -= log(mean(mfexp(srv_sel(sp))));
# #       srv_sel(sp)                      = mfexp(srv_sel(sp));
# #     }
# rm(list=ls())
# load(path.expand("data/BS_SS_Files/CEATTLE_results.Rdata"))
# mod<-tmp
# tmp<-read.csv(path.expand("data/BS_SS_Files/ceattle_est.std"),header=T)
# tt<-strsplit(as.character(tmp[1,]),split=" ")[[1]]
# tt<-tt[tt!=""]
# parms<-data.frame(parN=rep(0,dim(tmp)[1]),name="t",est=0,sd=0)
# parms$name<-as.character(parms$name)
# i<-1
# parms[i,1]<-as.numeric(tt[1])
# parms[i,2]<-as.character(tt[2])
# parms[i,3]<-as.numeric(tt[3])
# parms[i,4]<-as.numeric(tt[4])
#
#
# for(i in 2:dim(parms[1])){
#   tt<-strsplit(as.character(tmp[i,]),split=" ")[[1]]
#   tt<-tt[tt!=""]
#   parms[i,1]<-as.numeric(tt[1])
#   parms[i,2]<-as.character(tt[2])
#   parms[i,3]<-as.numeric(tt[3])
#   parms[i,4]<-as.numeric(tt[4])
#
# }
#
# fsh_sel_coff<-parms[grep("fsh_sel_coff",parms[,2]),3][1:8]
# nages<-12
# nselages<-8 #(set on line 436)
#
#   mod$fsh_sel_1
#   mod$srv_sel_1
#   par(mfrow=c(1,2))
#   norm_fsh_like<-50*(1)^2+  50*(log(5))^2+  50*(1)^2   # try playing around with this by changing 5 to 1, or 0.5
#
#   norm_fsh_like<-0
#   Fsel<-function(fsh_sel_coff,nselages,nages){
#     fsh_sel<-rep(0,nages)
#     fsh_sel[1:nselages]<-fsh_sel_coff                #// fsh_sel_coff(sp) parameters
#     fsh_sel[nselages+1:nages]<-fsh_sel[nselages]    #// set to last sel for those ages with less info
#     plot(exp(fsh_sel),type="b",pch=16)
#     avgsel_fsh<- log(mean(exp(fsh_sel_coff)))  # // average fisheries selectivity - used to normalize selectivity
#     norm_fsh_like = 50 * (avgsel_fsh)^2  #// To normalize selectivities across spp
#
#     fsh_sel<-fsh_sel-log(mean(exp(fsh_sel)))  # this is done to help recenter shaper around mean fsh sel
#     fsh_sel <- exp(fsh_sel)
#     plot((fsh_sel),type="b",pch=16)
#   }
#   norm2<-function(v){
#     nn<-length(v)
#     nu<-rep(0,nn)
#     for(i in 1:nn){
#       nu[i]<-(v[i]*i)^2
#     }
#     return(sum(nu))
#   }
#
#   first_difference<-function(v){
#     nn<-length(v)
#     return( v[2:(nn)]-v[1:(nn-1)])
#   }
#   curv_pen_fsh<- 12.5
#
#   # norm_srv_like(sp) += 50. * square(avgsel_srv(sp));                               // To normalize selectivities
#   # norm_fsh_like(sp) += 50. * square(avgsel_fsh(sp));                               // To normalize selectivities
#   #
#   # fsh_cat_like(sp)  += 200.*norm2(log(tc_biom_obs(sp)+1.e-4)-log(tc_biom_hat(sp)+1.e-4));                     // Errors in total catch estimation
#
#   #// ----Selectivity penalties...for smooth second diff etc
#   #//           Invoke a penalty when the partial F's go down with sp_age
#   #for (int sp_age=1;sp_age<nages(sp);sp_age++)
#   #if (fsh_sel(sp,sp_age)>fsh_sel(sp,sp_age+1))
#   #fsh_sel_like(sp) += 20.*square(log(fsh_sel(sp,sp_age)/fsh_sel(sp,sp_age+1)));
#   #fsh_sel_like(sp) += curv_pen_fsh(sp)*norm2(first_difference( first_difference(log(fsh_sel(sp)))));
#
#   sel_pen<-function(nages,fsh_sel,smoothpen=TRUE){
#     fsh_sel_like<-0
#     for (sp_age in 1:(nages-1)){
#       if(fsh_sel[sp_age]>fsh_sel[sp_age+1]){
#         # if the fsh_sel value changes (up or down) penalize (this helps avoid large jumps)
#         fsh_sel_like<-fsh_sel_like+ 20*(log(fsh_sel[sp_age]/fsh_sel[sp_age+1]))^2
#       }
#     }
#     if(smoothpen)
#       fsh_sel_like<-fsh_sel_like+curv_pen_fsh*norm2(first_difference( first_difference(log(fsh_sel))))
#     return(fsh_sel_like)
#
#   }
#   sel_pen(nages,fsh_sel,smoothpen=F)
#   sel_pen(nages,fsh_sel,smoothpen=T)
#
#
#
#
#
#
