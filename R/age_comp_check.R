# Function to return relative error of estimated survey age comps relative to that estimated in CEATTLE ADMBusing outputs from TMB or ADMB based CEATTLE. Trying to match how srv_age_hat in ADMB is calculated.

# species == Whichs species of interest: 1 = pollock, 2 = cod, 3 = ATF.
# ADMB_TMB == Which model output: 1 = Grant's TMB, 2 = Kirstin's ADMB
# data_list == CEATTLE inputs modified for TMB
# rep == TMB CEATTLE output
# tmp == ADMB CEATTLE output

# Trying to mimic the following lines (lines 2223-2235 in source) of ADMB CEATTLE:
# for (int yr=1;yr<=nyrs_srv_age(sp);yr++){
#   yr_ind  = yrs_srv_age(sp,yr) - styr + 1;// convert years into indices for 1-=nyrs counting purposes
#   tmp_age = N(sp,yr_ind);// need to test for type of sp_age data here (if length bins are different
#   tmp_age = elem_prod(srv_sel(sp) , tmp_age);
#   if ( srv_age_type(sp)==1)
#   srv_age_hat(sp,yr) = tmp_age / sum(tmp_age);
#   else{
#   srv_age_hat(sp,yr)  = tmp_age * age_trans_matrix(sp);
#   srv_age_hat(sp,yr) /= sum(srv_age_hat(sp,yr));
#   }
# }

age_comp_comparison <- function(data_list, species, rep, tmp, ADMB_TMB){

  # TMB DATA
  if(ADMB_TMB == 1){
    NByage <- rep$NByage[,,species]
    F = rep$F[,,species]
    M1 = rep$M1[species,]
    srv_sel = rep$srv_sel[species,]
    ALK <- data_list$age_trans_matrix[,,species]
  }

  # ADMB DATA
  if(ADMB_TMB == 2){
    NByage <- tmp[[paste("NByage", species, sep = "_")]]
    F = tmp[[paste("F", species, sep = "_")]]
    M1 = tmp[[paste("M1", species, sep = "_")]]
    srv_sel = tmp[[paste("srv_sel", species, sep = "_")]]
    ALK <- data_list$age_trans_matrix[,,species]
  }

  # Calculate age comp if pollock
  srv_age_hat <- sweep(NByage, 2, srv_sel, "*") # sweep(NByage * exp(-0.5 * sweep(F, 2, M1, "+")), 2, srv_sel, "*") Half year is not included in ceattle
  if(species == 1){
    srv_age_com <- (srv_age_hat/ rowSums(srv_age_hat, na.rm = T))
    srv_age_com <- srv_age_com[1:36, 1:12]
  }
  # Calculate length comp if cod and ATF
  if(species > 1){
    srv_len_hat <- matrix(NA, nrow = 36, ncol = 25)
    for(i in 1:36){
      srv_len_hat[i,] <- srv_age_hat[i,1:c(12,21)[species-1]] %*% ALK[1:c(12,21)[species-1], 1:25]
    }
    srv_age_com <- (srv_len_hat/ rowSums(srv_len_hat, na.rm = T))
  }

  # Calculate relative error and return
  rel_error <- (srv_age_com - tmp[[paste("srv_age_hat", species, sep = "_")]] )/ tmp[[paste("srv_age_hat", species, sep = "_")]]
  return(rel_error)
}

# Evaluate it!
load("ADMB_TMB_CEATTLE_inputs_and_outputs.RData")
results_tmb <- age_comp_comparison(data_list, species = 1, rep, tmp, ADMB_TMB = 1)
results_admb <- age_comp_comparison(data_list, species = 1, rep, tmp, ADMB_TMB = 2)
