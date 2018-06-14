#START filenames (don't alter this line, for this section no spaces after #)
#File
ceattle
#datafile_name
dat_input_files/01_assesment_2017/test_2Jim.dat
#dietfile_name
#dat_input_files/diet.dat
#diet2file_name 
dat_input_files/01_assesment_2017/diet2.dat
#stomfile_name 
dat_input_files/01_assesment_2017/stomach2017.dat
#Fprofile_datafile
dat_input_files/01_assesment_2017/F_profile.dat
#ATFfile_name 
dat_input_files/01_assesment_2017/ATF_Misc.dat
#Recfile_name    #fits_4_CEATTLE/rs_data4CEATTLE_fallZtot.dat
fits_4_CEATTLE/rs_data4CEATTLE_full_avg.dat
#retrofile_name
dat_input_files/01_assesment_2017/retro_data2017_asssmnt.dat
#futfile_name
dat_input_files/01_assesment_2017/projection_data2017_asssmnt.dat
#catch_in_name
dat_input_files/01_assesment_2017/catch_in.dat
#setC_name 
dat_input_files/01_assesment_2017/set_catch.dat
#setF_name 
dat_input_files/01_assesment_2017/set_FabcFofl_0.dat
#END filenames (don't alter this line, needed for FIT_recruitment.R)
#================================================================================ 
#================================================================================ 
# Model setup for MSM                                         
#================================================================================ 
#================================================================================
#  
#debugg  : debugg the model (0 =NO, 1= dump reps)                                          
0
#msmMode  : model mode                                            
0
# 0 = run in  single  species mode                                    
# 2 = run in  MSM mode  using Juarado #NAME?  approach  but new diet  and stomach data                           
#nspp : number  of  species                                         
3 
#recMode  : recruitment mode                                            
0
# rand_rec : 0 project without error around recruitment, 1 project with randomly drawn errors
0
#harvestMode  : used  for the calcmort  function  in  the projection  / simulation  models  
2
#1= project under no  fishing,  
#2= project under mean  fishing rate  from  hindcast, 
#3= project to  F rate  that  matches Btarget,  
#4= fit to NewCatch.dat
#5= project under radnom seeded fishing rate  from  hindcast
#ntemp_scen  : number  of  future  simulations
1
#simset : set of  future  iterations  to  run                                   
1 # 2 3 4 5 6 7 8 9 10 11 12
#alpha_ABC: ADDED F40 alpha - min biomass alpha
0.05
#alpha_OFL: ADDED F35 alpha - min biomass alpha
0.00
#catch_scen: ADDED scenario for Amanda's catch function
1
#crnum  : number  of  control rules to  run                                   
1
#ncr_scen : number  of  control rule  Scenarios                                     
1
#c_mult : control rule  scenarios    
1 1 #Individual stocks  fished  to  B40% where B0 is set to the B0 from no climate (using B0_set in .ctl file)                                    
#1 1 #Individual stocks  fished  to  B40%;
#1 2 #Individual stocks  fished  to  B40%; constrained such  that  no  B0(spp) <.35  B0(spp)
#1 22 #Individual stocks  fished  to  B40%; constrained such  that  no  B0(spp) <.35  B0(spp); given 2 mT cap & ATF ~ mean F (i.e., plk C<=1.35 mT, and pcod=0.95*ABC) 
#1 3 #Individual stocks  fished  to  B40%; pcod  and atf,  then  plk
#1 4 #Individual stocks  fished  to  B40%; pcod, then  atf,  then  plk
#1 5 #Individual stocks  fished  to  B40%; pcod  and plk,  then  atf
#1 6 # To ADD: Individual stocks  fished  to  B40%; given 2 mT cap & ATF ~ mean F (i.e., plk C<=1.35 mT, and pcod=0.95*ABC) 
#1 7 # To ADD: Individual stocks  fished  to  B40%; given 2 mT cap (i.e., plk C<=1.35 mT, and pcod=0.95*ABC) 
#1 8 # To ADD: Individual stocks  fished  to  B40%; ATF ~ mean F
#4 3  # MSY with min B35%
#5 3  # MEY with min B35%
#2 1 #calculated project using YPR method
#3 1 #sum  F Biomass = 0.35% sum (B0); scaled  given average M's
#3 2 #sum  F Biomass = 0.35% sum (B0); scaled  given average M's,  constrained such  that  no  B0(spp) <.35  B0(spp)
#3 3 #sum  F Biomass = 0.35% sum (B0); scaled  given B0
#3 4 #sum  F Biomass = 0.35% sum (B0); scaled  given B0, constrained such  that  no  B0(spp) <.35  B0(spp)
#4 1 # solve for individual msy  
#4 2 #solve  for system  wide  msy  (key  spp); unconstrained  
#4 3 #solve  for system  wide  msy  (key  spp); constrained such  that  no  B0(spp) <.35  B0(spp)
#4 4 #solve  for system  wide  Bmsy  (key  spp); Bmsy  based on  realistic effort  / value of  stocks
#2  1 #calculated F 35% for M at  B0  for each  stock using SR covariates from MSM_SR model fits to data
#rationMode : ration  mode                                            
2
# 0 = use static  rations (read in  from  data)                           
# 1 = use Elliott & presson temp  rations (temp and weigt specific)
# 2 = use bioenergetic  rations (temp and weight) & p-value
#M2mode :0= straight  ratio,  1=  use ADMB  liklihood to  find  M2, 3=  iterate for M2  , 4 uses  jesus's approach  of  suitability       
0 
#Umode  : U calculation mode; 0 = use Uobs  read  in  from  data, 1 = Calculate Unew  (only under msmMode==3)
0 
#maxphase : maximum number  of  phases  to  run                                   
5 
#RationByAge  : 1 = calculate ration  based mean  length  at  age,  0 = calculate ration  and foraging  based on  lengths,  convert to  ages      
1
#Btarget  : target  spawning  biomass for control rules (ie 0.35  such  that  B/B0  =.35)
0.4 
#B0_set: ADDED target  spawning  biomass for control rules (ie 0.35  such  that  B/B0  =.35)
# single spp B0 based on CR 1.1 and no climate 
5386284.5  447268.5  493190.2
# multispp B0 based on CR 1.1 
#3812316.4  398083.9  460652.7
#repn2: number  of  reps  to  find  M2                                      
30  
#rep_in:  number  of  reps  for starting  up  from  a #NAME?  command line (should be 5 or more in MSM mode)
30
#styr :First  year                                              
1979   
#nyrs :number of  years                                           
39
#nyrs_est :number of years for estimation (will be same as nyrs unless fitting to a subset)                                       
39
#logist_sel_phase : run logistic  selectivity for the survey  -2  is  off,  2 is  on  - opposite  is  the phase sel_phase is  calculated  in  2   
3 3 3      
#FP_in : Input F rates for each species for harvest mode 
.2 .2 .2                                          
#nyrs_fut2  :not used - use nyrs_fut in the futureT.dat file of future temperatures // number of  projection  years
#50
#clt_test_num
12345


