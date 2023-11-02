nages = 10
comp_obs = 1:10
flt_accum_age_lower = 1
flt_accum_age_upper = 10

#Accumulation age
for(age in 1:nages){
   age_test = age
  if(age > nages){
    age_test = age - nages
  }


  # Lower
  if(age_test < flt_accum_age_lower){
    comp_obs[flt_accum_age_lower + nages] =  comp_obs[flt_accum_age_lower + nages] + comp_obs[age]
    comp_obs[age] = 0

  }
  # Upper
  if(age_test > flt_accum_age_upper){
    comp_obs[flt_accum_age_upper + nages] =  comp_obs[flt_accum_age_upper + nages] + comp_obs[age]
    comp_obs[age] = 0
  }
}
