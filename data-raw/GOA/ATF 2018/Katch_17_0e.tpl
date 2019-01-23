//                      atfbsai2.tpl
//           Uses the shelf survey, slope survey and the Aleutian Islands (fits biomass and length comps)
//           modified version of the GOA model with domed shaped selectivity for the shelf survey males,
//           this run tries to estimate the temperature
//           effect on shelf survey catchability  
//			1 is female, 2 is male  in model (not in survey data)   
			//this is how to change the name of the data file
// when you run it add 			./atfbsai_is2014_6 -ind atfbsai_is2014_4.dat
DATA_SECTION
!!CLASS ofstream evalout("atf_gen2.mcmc.out");       //wasgen
  init_int styr         //(1) start year of model
  init_int endyr        //(2) end year
  init_int styr_fut     //(3) start year of projections (endyr+1) 
  init_int endyr_fut    //(4) end year of projections 
  init_int assess;  //(4.5) 1 for BSAI 2 for GOA 3 for Kamchatka   
  init_int nsurv   //(5)     
  !!cout<<"nsurv"<<nsurv<<std::endl;
  init_number median_rec  //(6) median recruit value to use for the last 3 years
 
  init_int first_age //(7) first age to use 
  init_int last_age //(8) 
  init_int first_length//8.1 first length to use
  init_int last_length//8.2 last length to use (for the GOA there are 21 lengths, for BSAI there are 25) for 3-15+ only first was not used previously
  init_int nages_read//(8.3) number of ages read    
  !!cout<<"nages_read"<<nages_read<<std::endl;
  int nages; //(8.4)
  !! nages=last_age-first_age+1;          // # of ages in the model  
  !!cout<<"nages"<<nages<<std::endl;
  init_int nsurv_aged //(9) 
  !!cout<<"nsurv_aged"<<nsurv_aged<<std::endl;    
 //phases
  init_int phase_F40      //(10) phase F40 is estimated 
  init_ivector phase_logistic_sel(1,2*nsurv) //(11)  
  init_ivector phase_fishery_sel(1,2) //(12) phase to do selectivities for fishery  
  init_int phase_alphabeta //(13) phase to estimate alpha and beta
  init_ivector q_Phase(1,nsurv); //(14)  
  init_int phase_selcoffs      //(15) generally set to phase 4 phase for smooth selectivity curve fishery
  !!cout<<"phase_selcoffs"<<phase_selcoffs<<std::endl;         
 //selectivity parameters 
  init_int nselages       //(16) fishery (for asymptotic selectivity) set to 19 selectivity is set to the selectivity at nselages-1 after age nselages
  init_ivector nselages_srv(1,nsurv) //(17)
  init_vector fishsel_LB_f(1,2);//(18)fishery selectivity lower and upper bounds
  init_vector fishsel_LB_m(1,2);//(19)
  init_vector fishsel_UB_f(1,2);//(20)
  init_vector fishsel_UB_m(1,2);//(21) 
  init_vector fishsel_prior_f(1,2);//(22) fishsel_prior_f 
  init_vector fishsel_prior_m(1,2);//(23) fishsel_prior_m 
  init_ivector nsel_params(1,nsurv);//(24) number of selectivity parameters for each survey (either 2 or 4)
  !!cout<<"nsel_params"<<nsel_params<<std::endl;
  init_vector sel_prior_f(1,2*nsurv);//(25) sel_prior_f
  init_vector sel_prior_m(1,2*nsurv);//(26) sel_prior_m
  !!cout<<"sel_prior_f"<<sel_prior_f<<std::endl;  
  init_vector sel_LB_f(1,2*nsurv);//(27)   
  init_vector sel_LB_m(1,2*nsurv);//(28)   
  init_vector sel_UB_f(1,2*nsurv);//(29)   
  init_vector sel_UB_m(1,2*nsurv);//(30) 
 //parameters below are for the descending logistic of the shelf survey (srv1) for the BSAI survey. I just leave them in for the GOA survey but they are not used.
  init_vector sel1_desc_prior_f(1,2);//(31)
  init_vector sel1_desc_prior_m(1,2);//(32)
  init_vector sel1_desc_LB_f(1,2);//(33)
  init_vector sel1_desc_LB_m(1,2);//(34)
  init_vector sel1_desc_UB_f(1,2);//(35)
  init_vector sel1_desc_UB_m(1,2);//(36) 

 //sample size for length comps for weighting likelihoods  
  init_int nlen             //(37) # of length bins
  init_int nobs_fish          //(38) # of years of fishery data
  init_ivector yrs_fish(1,nobs_fish)   //(39) years with fishery data 
  !!cout<<"yrs_fish"<<yrs_fish<<std::endl; 
  init_matrix nsamples_fish(1,2,1,nobs_fish)  //(40) sample size (weights) for each sex and yr of fishery data
  init_ivector nobs_srv(1,nsurv) //(41) # yrs of shelf, slope, AI data
  init_imatrix yrs_srv(1,nsurv,1,nobs_srv) //(42) yrs with shelf, slope, AI survey data
  init_ivector nobs_srv_length(1,nsurv) //(43)  
  !!cout<<"nobs_srv_length"<<nobs_srv_length<<std::endl; 
  init_imatrix yrs_srv_length(1,nsurv,1,nobs_srv_length) //(44) yrs with shelf, slope, AI length data
  !!cout<<"yrs_srv_length"<<yrs_srv_length<<std::endl;
  init_imatrix nsamples_srv_length_fem(1,nsurv,1,nobs_srv_length)   //(45) sample sizes for each length comp by sex and year for shelf, slope, AI survey  
  init_imatrix nsamples_srv_length_mal(1,nsurv,1,nobs_srv_length)   //(46)     
  init_3darray obs_p_fish(1,2,1,nobs_fish,1,nlen)  //(48) fishery length comps 
  !!cout<<"obs_p_fish"<<obs_p_fish(1)<<std::endl;  
  init_3darray obs_p_srv_length_fem(1,nsurv,1,nobs_srv_length,1,nlen)  //(49) survey length comps by survey (shelf, slope, AI, bin, sex and yr)
  init_3darray obs_p_srv_length_mal(1,nsurv,1,nobs_srv_length,1,nlen)  //(50)
  init_vector catch_bio(styr,endyr)    //(51) catch by year
  !!cout<<"catch_bio"<<catch_bio<<std::endl;
  init_imatrix obs_srv(1,nsurv,1,nobs_srv) //(52) survey biomass by year (shelf, slope, AI)
  !!cout<<"obs_srv"<<obs_srv<<std::endl;
  init_imatrix obs_srv_sd(1,nsurv,1,nobs_srv) //(53) survey SE by year    
  init_matrix wt(1,2,1,nages)          //(54) weight-at-age by sex   
  !!cout<<"wt"<<wt<<std::endl;  
  init_vector maturity(1,nages)        //(55) female prop. mature-at-age
  !!cout<<"maturity"<<maturity<<std::endl;
  //length age transition matrix 
  //goa: need length age transition matrix with 21 ages and 21 length bins
  init_3darray lenage(1,2,1,nages,1,nlen)  //(56) length-age transition matrix
  !!cout<<"lenage"<<lenage<<std::endl; 
  int nyrs_temps; 
  !!if(assess==1){nyrs_temps = nobs_srv(1);}  
  !!if(assess==2) {nyrs_temps = 33;}//nobs_srv(1); change back for BSAI assessment
  !!if(assess==3){nyrs_temps = nobs_srv(1);}
  !!cout<<"nyrs_temps"<<nyrs_temps<<std::endl; 
  init_vector bottom_temps(1,nyrs_temps); //nobs_srv(1))    //(57) shelf survey bottom temperatures
  !!cout<<"nyrs_temps"<<nyrs_temps<<std::endl;
  init_int monot_sel     //(58) selectivity smoothing function for fishery 
  !!cout<<"monot_sel"<<monot_sel<<std::endl;
  init_vector wt_like(1,8)    //(59) 
  !!cout<<"wt_like"<<wt_like<<std::endl;               
  init_ivector nobs_srv_age(1,nsurv_aged) //(60) # yrs with survey ages 
  init_imatrix yrs_srv_age(1,nsurv_aged,1,nobs_srv_age) //(61) yrs of shelf, ai survey ages
  init_3darray nsamples_srv_age(1,nsurv_aged,1,2,1,nobs_srv_age) //(62) sample sizes of ages read each year by sex and survey
  init_3darray obs_p_srv_age_fem(1,nsurv_aged,1,nobs_srv_age,1,nages) //(63) survey age comps by sex and year females  
  init_3darray obs_p_srv_age_mal(1,nsurv_aged,1,nobs_srv_age,1,nages) //(64) survey age comps by sex and year males   
  init_vector M(1,2) //(65) female then male natural mortality
  !!cout<<"M"<<M<<std::endl;            
  init_number offset_const //(66) a constant to offset zero values
  init_vector q_Lower_bound(1,nsurv); //(67)
  init_vector q_Upper_bound(1,nsurv); //(68)  

  init_vector q_surv_prior_mean(1,nsurv);  //(69) 
  init_ivector nparams_srv(1,nsurv); //(70) this tells you whether you have 2 or 4 parameters for each survey
  !!cout<<"q_surv_prior_mean"<<q_surv_prior_mean<<std::endl;  
  init_int mean_log_rec_prior; //(72)  
  init_int log_avg_fmort_prior; //(73)   
  init_vector like_wght(1,7)    //(74) 
  init_vector fpen_mult(1,2) //(75) 
  init_number catch_err;//(76) 
  init_number fmort_boundL //(77) 
  init_number fmort_boundH //(78) 
  init_matrix agerr_matrix(1,nages_read,1,nages_read) //agerr_matrix 
  int styr_rec;  
  number catch_err_like; 
  matrix cv_srv(1,nsurv,1,nobs_srv);  //matrix to hold CVs for surveys
//not all consistent throughout - could work on.
//year
  int i
//age
  int j
//sex
  int k
//
  int ii
  int m  


 LOCAL_CALCS
   styr_rec=styr-nages+1;
   if(nselages>nages) 
   {nselages=nages;  
   cout<<"Warning selectivity: is set to be estimated on more ages than are in the model."<<std::endl;  }
   for (i=1; i<= nsurv; i++){
   if(nselages_srv(i)>nages) nselages_srv(i)=nages;
   }
   //calculate cv for surveys
   for (int j=1;j<=nsurv;j++){
   for (i=1;i<=nobs_srv(j);i++){ 
   cv_srv(j,i)=obs_srv_sd(j,i)/(double)obs_srv(j,i); }} 
   //change weights to tons
   wt=wt*.001;
  catch_err_like=.5/(catch_err*catch_err);

 END_CALCS
   
  vector obs_sexr(1,nobs_fish)  // prop. males in fishery length data
  matrix obs_sexr_srv_2(1,nsurv,1,nobs_srv_length)  //proportion males in survey data 
  number obs_mean_sexr    //average proportion of males in shelf survey population estimates
  number obs_SD_sexr      //standard deviation from male prop. in shelf survey population estimates
  vector pred_sexr(styr,endyr)   //proportion of males in num at age matrix to be calculated
//  vector fishsel_params(1,2)     //so far always 4 params for fishery data
    
INITIALIZATION_SECTION
  //can have different mortality for males and females  
  F40 .20
  F35 .21
  F30 .23
  mean_log_rec mean_log_rec_prior
  log_avg_fmort log_avg_fmort_prior  
  //proportion in each region constrained with catchability so it does not add to 1. Expect to add to less than 1?
  q_surv q_surv_prior_mean  
//  fishsel_params_f fishsel_params_prior 
  fmort_dev 0.00001   
//note: you can initialize things you do not use. 
//selectivity parameter vectors initialize
  fishsel_params_f fishsel_prior_f
  fishsel_params_m fishsel_prior_m
  srv_params_f sel_prior_f 
  srv_params_m sel_prior_m
  srv1desc_params_f sel1_desc_prior_f 
  srv1desc_params_m sel1_desc_prior_m  
  alpha 1.
  beta 0.  


PARAMETER_SECTION
 //parameters to be estimated are all ones that begin with init_ and have a positive
 //phase, negative phase means are fixed.
 //phase of 8 is greater than last phase so does q1 in last phase  
  // init_bounded_number q1(.5,2,8)
 //fix q1 to be 1 otherwise it went to lower bound of .5

  init_bounded_number_vector q_surv(1,nsurv,q_Lower_bound,q_Upper_bound,q_Phase)
  init_bounded_number_vector fishsel_params_f(1,2,fishsel_LB_f,fishsel_UB_f,phase_fishery_sel)
  init_bounded_number_vector fishsel_params_m(1,2,fishsel_LB_m,fishsel_UB_m,phase_fishery_sel) 
  init_bounded_number_vector srv_params_f(1,2*nsurv,sel_LB_f,sel_UB_f,phase_logistic_sel);
  init_bounded_number_vector srv_params_m(1,2*nsurv,sel_LB_m,sel_UB_m,phase_logistic_sel);    
  init_bounded_number_vector srv1desc_params_f(1,2,sel1_desc_LB_f,sel1_desc_UB_f,phase_fishery_sel); //phase is just 2
  init_bounded_number_vector srv1desc_params_m(1,2,sel1_desc_LB_m,sel1_desc_UB_m,phase_fishery_sel);
//  init_bounded_number_matrix srv_params_m(1,nsurv,1,nsel_params,sel_LB_m,sel_UB_m)   
  init_number alpha(phase_alphabeta)       // used to estimate temperature effect on shelf survey catchability
  init_number beta(phase_alphabeta)  // used to estimate temperature effect on shelf survey catchability
 //phase of -1 means M is fixed   
  init_number mean_log_rec(1)
  init_bounded_dev_vector rec_dev(styr_rec,endyr,-15,15,2) //JNI    
  init_number log_avg_fmort(1)
  init_bounded_dev_vector fmort_dev(styr,endyr,fmort_boundL,fmort_boundH,1)  //was -3,3,1 change IS 2015 
//  Selectivity parameters from the GOA version of the model

  init_matrix log_selcoffs_fish(1,2,1,nselages,phase_selcoffs)     
  init_bounded_number sexr_param_fish(1.0,1.0,-5)  //this was hitting bound of 1.0 so fixed it - should free up to check

// Parameters for computing SPR rates 
  init_bounded_number F40(0.01,1.,phase_F40)
  init_bounded_number F35(0.01,1.,phase_F40)
  init_bounded_number F30(0.01,1.,phase_F40)
  matrix log_sel_fish(1,2,1,nages)
  matrix sel(1,2,1,nages) //fishery selectivity 
  3darray sel_srv(1,2,1,nsurv,1,nages) //try 3d array here
  vector avgsel_fish(1,2)
  matrix popn(1,2,styr,endyr)  
  3darray totn_srv(1,nsurv,1,2,styr,endyr)  //new matrix combine total numbers over 3 surveys 
  vector temp1(1,nages)
  vector temp2(1,nages)
  vector explbiom(styr,endyr)
  vector pred_bio(styr,endyr)
  sdreport_vector fspbio(styr,endyr) 
  matrix pred_srv(1,nsurv,styr,endyr) //matrix combine pred_srv into 3 surveys
  3darray pred_p_fish(1,2,styr,endyr,1,nlen)
  3darray pred_p_srv_age_fem(1,nsurv_aged,1,nobs_srv_age,1,nages)//pred_p_srv_age for males and females for each survey
  3darray pred_p_srv_age_mal(1,nsurv_aged,1,nobs_srv_age,1,nages)//same but males
  3darray pred_p_srv_len_fem(1,nsurv,1,nobs_srv_length,1,nlen)//pred_p_srv_length for males and females for each survey 
  3darray pred_p_srv_len_mal(1,nsurv,1,nobs_srv_length,1,nlen)  //same but males
  vector pred_catch(styr,endyr)
  3darray natage(1,2,styr,endyr,1,nages) 
  sdreport_vector totalbiomass(styr,endyr)
  3darray catage(1,2,styr,endyr,1,nages)
  3darray Z(1,2,styr,endyr,1,nages)
  3darray F(1,2,styr,endyr,1,nages)
  3darray S(1,2,styr,endyr,1,nages)
  vector fmort(styr,endyr)
  number rbar
  vector surv(1,2)  //survival for each sex
  vector offset(1,7) 
  number rec_like
  number catch_like
  number sexr_like   
  vector age_like(1,nsurv_aged) //really only need shelf and AI but may need other elements later
  vector length_like(1,nsurv+1)
  vector sel_like(1,4)
  number fpen 
  vector surv_like(1,nsurv) //survey likelihood for each survey   
  sdreport_number endbiom
  sdreport_number depletion
  objective_function_value obj_fun
  number tmp
  vector pred_sexr(styr,endyr)
 // Stuff for SPR and yield projections
  number sigmar
  number ftmp
  number SB0
  number SBF40
  number SBF35
  number SBF30
  number sprpen
  matrix Nspr(1,4,1,nages)
  3darray nage_future(1,2,styr_fut,endyr_fut,1,nages)
  sdreport_matrix fspbiom_fut(1,4,styr_fut,endyr_fut)
  3darray F_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray Z_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray S_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray catage_future(1,2,styr_fut,endyr_fut,1,nages)
  number avg_rec_dev_future
  vector avg_F_future(1,4)
  sdreport_matrix catch_future(1,3,styr_fut,endyr_fut) // Note, don't project for F=0 (it 
  sdreport_matrix future_biomass(1,4,styr_fut,endyr_fut)
  vector explbiom_fut(styr_fut,endyr_fut)
  number maxsel_fish
  vector maxsel_srv(1,nsurv)
  number mlike
  number qlike
  number flike
  vector qtime(styr,endyr) 
 
//  vector obs_sexr(1,nobs_fish);  

PRELIMINARY_CALCS_SECTION  

  obs_mean_sexr=0.34;  //initial value for avg proportion of male population estimated from shelf surveys; calculated below
  obs_SD_sexr=0.0485;  //initial value for standard deviation of mean male population proportion: calculated below
//sex ratio in the fishery 
  
  for(i=1; i<=nobs_fish;i++)
  {
    obs_sexr(i) = sum(obs_p_fish(1,i))/sum(obs_p_fish(1,i) + obs_p_fish(2,i)); 
  }  

//length obs sex ratio in surveys (all combined); proportion of males     
  for(i=1;i<=nsurv;i++)
  {
	for (j=1;j<=nobs_srv_length(i);j++)
	{    
    	obs_sexr_srv_2(i,j)=sum(obs_p_srv_length_mal(i,j)/
    	      (sum(obs_p_srv_length_mal(i,j))+sum(obs_p_srv_length_fem(i,j))));
	}
  }
  obs_mean_sexr=mean(obs_sexr_srv_2); //previously was just estimated from shelf survey data so kept that here.
  obs_SD_sexr=std_dev(obs_sexr_srv_2(1)); 

 //Compute offset for multinomial and length bin proportions
 // offset is a constant nplog(p) is added to the likelihood     
 // magnitude depends on nsamples(sample size) and p's_

  offset.initialize();

//fishery offset
  if(assess<3){
  for (i=1; i <= nobs_fish; i++)
  {
  double sumtot ;  
  sumtot = sum(obs_p_fish(1,i)+obs_p_fish(2,i)); 
  obs_p_fish(1,i) = obs_p_fish(1,i) / sumtot; 
  obs_p_fish(2,i) = obs_p_fish(2,i) / sumtot; 
  for(k=1; k<=2;k++) 
  {
    offset(1) -= nsamples_fish(k,i)*obs_p_fish(k,i) * log(obs_p_fish(k,i)+offset_const); //this multiplies elements together then sums them up.
  } 
  } 
  }

  if(assess==3)
  {
	for (i=1; i <= nobs_fish; i++)
  {
  double sumtot ;
  sumtot = sum(obs_p_fish(1,i)+obs_p_fish(2,i));
  obs_p_fish(1,i) /= sum(obs_p_fish(1,i));
  obs_p_fish(2,i) /= sum(obs_p_fish(2,i));
  for(k=1; k<=2;k++)
    offset(1) -= nsamples_fish(k,i)*obs_p_fish(k,i) * log(obs_p_fish(k,i)+.0001);
  }
  }   
 
//survey length offset and bin proportions 
  //this loops over all surveys and makes sure all proportions sum to 1.
  for(i=1;i<=nsurv;i++){
	for(j=1;j<=nobs_srv_length(i);j++){    
		double sumtot;
		sumtot=sum(obs_p_srv_length_fem(i,j)+obs_p_srv_length_mal(i,j));
        obs_p_srv_length_mal(i,j)=obs_p_srv_length_mal(i,j)/sumtot;  //changing these to proportions rather than numbers
        obs_p_srv_length_fem(i,j)=obs_p_srv_length_fem(i,j)/sumtot;
        offset(i+1)-= nsamples_srv_length_fem(i,j)*obs_p_srv_length_fem(i,j)*log(obs_p_srv_length_fem(i,j)+offset_const)
                   +nsamples_srv_length_mal(i,j)*obs_p_srv_length_mal(i,j)*log(obs_p_srv_length_mal(i,j)+offset_const); 
	}
  } 

  //survey age offsets
  for (i=1;i<=nsurv_aged;i++)
  {
    for (j=1;j<=nobs_srv_age(i);j++)
    {
	double sumtot;
	sumtot=sum(obs_p_srv_age_fem(i,j)+obs_p_srv_age_mal(i,j));
	obs_p_srv_age_fem(i,j)=(obs_p_srv_age_fem(i,j)/sumtot)*agerr_matrix;
	obs_p_srv_age_mal(i,j)=(obs_p_srv_age_mal(i,j)/sumtot)*agerr_matrix;
	offset(i+nsurv+1)-=nsamples_srv_age(i,1,j)*obs_p_srv_age_fem(i,j)*log(obs_p_srv_age_fem(i,j)+offset_const)+
	             nsamples_srv_age(i,2,j)*obs_p_srv_age_mal(i,j)*log(obs_p_srv_age_mal(i,j)+offset_const);
    }  
  } 
                                                 
PROCEDURE_SECTION
//this is for bootstraping where qrun is a vector of q's from bootstrap irun is the 
//run number.  sets the q (=q1) for the run.
//   q1=qrun(irun);
   get_selectivity();     
   get_mortality();  
    surv(1)=mfexp(-1.0* M(1));
    surv(2)=mfexp(-1.0* M(2));   
   get_numbers_at_age();
   get_catch_at_age();

  if (active(F40))
    compute_spr_rates();
  if (last_phase())
  {
    Future_projections();      
  }
  if (sd_phase() || mceval_phase()) 
   { 
	Do_depend();
    if (mceval_phase())
    {
      evalout << obj_fun << " " ;
      for (i=styr;i<=endyr;i++)
      evalout<<  fspbio(i) << " ";
      for (i=styr;i<=endyr;i++)    
      evalout<< natage(1,i)*wt(1) + natage(2,i)*wt(2) <<" ";
      for (i=styr;i<=endyr;i++)
      evalout << 2*natage(1,i,1) <<" ";
      evalout <<  endl;
    }
  }

    evaluate_the_objective_function();  


FUNCTION get_selectivity

  if(active(log_selcoffs_fish))// 
 {           
    for(k=1;k<=2;k++)
    {
     for (j=1;j<=nselages;j++)
      {
        log_sel_fish(k,j)=log_selcoffs_fish(k,j);
      }
    }
   //sets selectivity of ages older than nselages to selectivity at nselages  
   for(k=1;k<=2;k++)
   {
     for (j=nselages+1;j<=nages;j++)
     {
       log_sel_fish(k,j)=log_sel_fish(k,j-1);
     }
   }

 for(k=1;k<=2;k++)
   {
     avgsel_fish(k)=log(mean(mfexp(log_selcoffs_fish(k))));
   }
 //vector=vector-scalar same as  vector-=scalar  
 //scaling selectivities by subracting the mean so exp(mean(s))=1.   
 //selectivities can be greater than 1 but mean is 1. 
 
  for(k=1;k<=2;k++)
    {
      log_sel_fish(k)-=log(mean(mfexp(log_sel_fish(k))));
      sel(k)=mfexp(log_sel_fish(k));
    }     
 }//  end if(active(log_selcoffs_fish))
  else
    { 
	   sel(1)=get_sel(fishsel_params_f(1),fishsel_params_f(2));  
       sel(2)=get_sel(fishsel_params_m(1),fishsel_params_m(2));  
    
     //logistic selectivity curve
          for (j=1;j<=nages;j++)
          { 
            if(j<=nselages)
             {
               sel(1,j)=1./(1.+mfexp(-1.*fishsel_params_f(1)*(double(j)-fishsel_params_f(2))));
               sel(2,j)=1./(1.+mfexp(-1.*fishsel_params_m(1)*(double(j)-fishsel_params_m(2))));
             }
            else
            {
             sel(1,j)=sel(1,j-1);
             sel(2,j)=sel(2,j-1);
            }    				
          } 
     }

  for(i=1;i<=nsurv;i++)   
  {
     if(nsel_params(i)==4)
     {              
     sel_srv(1,i) = get_sel(srv_params_f(i),srv_params_f(i+1),srv1desc_params_f(1),srv1desc_params_f(2));
     sel_srv(2,i) = get_sel(srv_params_m(i),srv_params_m(i+1),srv1desc_params_m(1),srv1desc_params_m(2));         
     }
     if(nsel_params(i)==2)
       {
       sel_srv(1,i) = get_sel(srv_params_f((2*i)-1),srv_params_f(2*i));
       sel_srv(2,i) = get_sel(srv_params_m((2*i)-1),srv_params_m(2*i)); 
       }        
  }  
  for (j=1;j<=nages;j++)   
  {           
	for (i=1;i<=nsurv;i++) 
	{
    if (j>nselages_srv(i))
    {
	sel_srv(1,i,j)=sel_srv(1,i,j-1); 
	sel_srv(2,i,j)=sel_srv(2,i,j-1); 
    } 
    } 
  }

//     logistic selectivity curves, asymptotic for fishery, slope survey and the Aleutian Islands but domed shape for shelf survey    
  
FUNCTION dvar_vector get_sel(const dvariable& slp, const dvariable& a50)
   {
	dvar_vector sel_tmp(1,nages);
    for (j=1;j<=nages;j++)  //this is selectivity for the surveys
    {  
    sel_tmp(j)=1./(1.+mfexp(-1.*slp*(double(j)-a50))); 
    }          
    return(sel_tmp);
   }          
         
FUNCTION dvar_vector get_sel(const dvariable& slp, const dvariable& a50, const dvariable& dslp, const dvariable& d50)
   {
	dvar_vector sel_tmp(1,nages);
   for (j=1;j<=nages;j++)  //this is selectivity for the surveys         
   {
	  sel_tmp(j) = 1./(1.+mfexp(-1.*slp*(double(j)-a50)));           
      sel_tmp(j) *= 1./(1.+mfexp(dslp*(double(j)-d50)));
   }
 	return(sel_tmp);
  }          

FUNCTION get_mortality 
  maxsel_fish=max(sel(1));     //1 is females
  if(maxsel_fish<max(sel(2)))  //if highest female selectivity is > male selectivity, make maxsel_fish=male high selectivity
      maxsel_fish=max(sel(2));

  	  if (assess==3)   //added to match Kamchatka assessment - consider revising this
	  {
	  maxsel_fish=1.0;
      }

  fmort = mfexp(log_avg_fmort+fmort_dev); 
  for(k=1;k<=2;k++)
  {
    for (i=styr;i<=endyr;i++)
    {
      F(k,i)=(sel(k)/maxsel_fish)*fmort(i);
      Z(k,i)=F(k,i) + M(k); 
    }
  } 
  S = mfexp(-1.0*Z); 

FUNCTION get_numbers_at_age

  for(i=1;i<=nsurv;i++)
  {
   maxsel_srv(i)=max(sel_srv(1,i));
   if(maxsel_srv(i)<max(sel_srv(2,i)))
   maxsel_srv(i)=max(sel_srv(2,i));
  }     
    if (assess==3)   //added to match Kamchatka assessment - consider revising this
  {    
  for(i=2;i<=nsurv;i++)  
  maxsel_srv(i)=1.0;  
  }   
  if (assess==3)   //added to match Kamchatka assessment - consider removing this
  {
  maxsel_srv(1)=max(sel_srv(1,1));
  sel_srv(1,1) /= maxsel_srv(1); // shelf survey normalized by 4th age class 
  maxsel_srv(1)=max(sel_srv(2,1));
  sel_srv(2,1) /= maxsel_srv(1); // shelf survey normalized by 4th age class 
  }
  int itmp;

 //calc initial population  
  for (j=1;j<nages;j++)
    {
      itmp=styr+1-j;
      natage(1,styr,j)=mfexp(mean_log_rec-(M(1)*double(j-1))+rec_dev(itmp));
      natage(2,styr,j)=mfexp(mean_log_rec-(M(2)*double(j-1))+rec_dev(itmp));
    }
    itmp=styr+1-nages;
  //last age    
    natage(1,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(1)*(nages-1)))/(1.- surv(1));
    natage(2,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(2)*(nages-1)))/(1.- surv(2));

 // Now do for next several years----------------------------------
  if(assess<3)
  {
  for (i=styr+1;i<=endyr;i++)
  {
    //for age 1 recruits in the last year use value read in from data file
    if(i<=(endyr-1))
    {
      natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
      natage(2,i,1)=natage(1,i,1);
    }
    else
    {
      natage(1,i,1)=median_rec;
      natage(2,i,1)=natage(1,i,1);
    }
  }
  }

  if(assess==3)
  {
   for (i=styr+1;i<=endyr;i++)
   {
  //for age 1 recruits in the last year use value read in from data file
   natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
   natage(2,i,1)=natage(1,i,1);
   }  
  }
 //numbers at age 

  for(k=1;k<=2;k++)
  {
    for (i=styr;i< endyr;i++)
    {
      for(j=1;j<nages;j++)
      {
        natage(k,i+1,j+1)=natage(k,i,j)*S(k,i,j); 
      }
      natage(k,i+1,nages)+=natage(k,i,nages)*S(k,i,nages);
      popn(k,i)= natage(k,i)*sel(k);
    }
    popn(k,endyr)=natage(k,endyr)*sel(k);
  }

  for (i=styr;i<=endyr;i++)
  {
      pred_sexr(i)=sum(natage(2,i))/(sum((natage(1,i)+natage(2,i))));  //calculation of prop. of males in pred. population 
  }

  //predicted survey values
  fspbio.initialize(); 
  qtime=q_surv(1); 
  
  for (j=1;j<=nsurv;j++)
 {
    for (i=styr;i<=endyr;i++)
   {
  fspbio(i) = natage(1,i)*elem_prod(wt(1),maturity);
  explbiom(i)=0.;
  pred_bio(i)=0.; 
  pred_srv(j,i)=0.;
  //catchability calculation for survey years
  if (assess==1 && (i>=1982) && (i-1981 <= nobs_srv(1)))      //JNI catchability calculation for survey years    
   {   
   qtime(i)=q_surv(1)*mfexp(-alpha+beta*bottom_temps(i-1981));
   }
  for(k=1;k<=2;k++)
    {
    if (j==1 && assess==1)
      {             
    pred_srv(j,i) += qtime(i)*(natage(k,i)*elem_prod(sel_srv(k,j)/maxsel_srv(j),wt(k)));maxsel_srv(j);   //shelf survey, dividing by the maxsel constrains female selectivity to be 1.0
      } 
    else 
      {
    pred_srv(j,i) += q_surv(j)*(natage(k,i)*elem_prod(sel_srv(k,j),wt(k)));///maxsel_srv(j);         //slope survey JNI  do not need to divide by maxsel_srv if it is logistic but does not hurt
      } 
     
       //Aleutian Islands survey JNI

    //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
    // are set different in the tpl file the program will take to value from the bin file and use that 
    explbiom(i)+=natage(k,i)*elem_prod(sel(k),wt(k))/maxsel_fish;
    pred_bio(i)+=natage(k,i)*wt(k);
      }
   }  
 }     

  //Fit survey length compositions
  for (i=1;i<=nsurv;i++)
  {
	for (j=1;j<=nobs_srv_length(i);j++)
	{
			ii=yrs_srv_length(i,j);
			pred_p_srv_len_fem(i,j)=q_surv(i)*elem_prod(sel_srv(1,i),natage(1,ii))*lenage(1);
			pred_p_srv_len_mal(i,j)=q_surv(i)*elem_prod(sel_srv(2,i),natage(2,ii))*lenage(2);
			dvariable sum_tot=sum(pred_p_srv_len_fem(i,j)+pred_p_srv_len_mal(i,j));
			pred_p_srv_len_fem(i,j)/=sum_tot;
			pred_p_srv_len_mal(i,j)/=sum_tot;
	 }
   }
     
  //Fit survey age composition
   for (i=1;i<=nsurv_aged;i++)
  {
	for (j=1;j<=nobs_srv_age(i);j++)
	{
			ii=yrs_srv_age(i,j);
			pred_p_srv_age_fem(i,j)=q_surv(i)*elem_prod(sel_srv(1,i),natage(1,ii));
			pred_p_srv_age_mal(i,j)=q_surv(i)*elem_prod(sel_srv(2,i),natage(2,ii));
			dvariable sum_tot=sum(pred_p_srv_age_fem(i,j)+pred_p_srv_age_mal(i,j));
			pred_p_srv_age_fem(i,j)/=sum_tot;
			pred_p_srv_age_mal(i,j)/=sum_tot;
	 }
   }

  depletion=pred_bio(endyr)/pred_bio(styr);
  endbiom=pred_bio(endyr);
 
FUNCTION get_catch_at_age
  for (i=styr; i<=endyr; i++)
  {
    pred_catch(i)=0.;
    for(k=1;k<=2;k++)
    {      
      //--Baranov's equation here-----------------------------------
      for (j = 1 ; j<= nages; j++)
      {
        catage(k,i,j) = natage(k,i,j)*F(k,i,j)*(1.-S(k,i,j))/Z(k,i,j);
        pred_catch(i) += catage(k,i,j)*wt(k,j);
      }
      pred_p_fish(k,i)=elem_prod(sel(k),natage(k,i))*lenage(k)/(popn(1,i)+popn(2,i));
    }
  }   
    
FUNCTION Future_projections
  for(k=1;k<=2;k++)
  {
    nage_future(k,styr_fut)(2,nages)=++elem_prod(natage(k,endyr)(1,nages-1),S(k,endyr)(1,nages-1));
    nage_future(k,styr_fut,nages)+=natage(k,endyr,nages)*S(k,endyr,nages);
   }
    future_biomass.initialize();
    catch_future.initialize();
    for (int l=1;l<=4;l++)
    {
      switch (l)
      {
        case 1:
          ftmp=F40;
          break;
        case 2:
          ftmp=F35;
          break;
        case 3:
          ftmp=F30;
          break;
        case 4:
          ftmp.initialize();
          break;
      }

// Get future F's
     for(k=1;k<=2;k++)
     {
      for (i=endyr+1;i<=endyr_fut;i++)
      {
        for (j=1;j<=nages;j++)
        {
          F_future(k,i,j) = (sel(k,j)/maxsel_fish)*ftmp;
          Z_future(k,i,j) = F_future(k,i,j)+M(k);
          S_future(k,i,j) = exp(-1.*Z_future(k,i,j));
        }
      }
// Future Recruitment (and spawners)
      for (i=styr_fut;i<endyr_fut;i++)
      {
        nage_future(k,i,1)  = median_rec;
       // Now graduate for the next year....
        nage_future(k,i+1)(2,nages) = ++elem_prod(nage_future(k,i)(1,nages-1),S_future(k,i)(1,nages-1));
        nage_future(k,i+1,nages)   += nage_future(k,i,nages)*S_future(k,i,nages);
      }
      nage_future(k,endyr_fut,1)  = median_rec;
      // Now get catch at future ages
      for (i=styr_fut; i<=endyr_fut; i++)
      {
        for (j = 1 ; j<= nages; j++)
        {
          catage_future(k,i,j) = nage_future(k,i,j) * F_future(k,i,j) * ( 1.- S_future(k,i,j) ) / Z_future(k,i,j);
         if(k==1)
          {
          fspbiom_fut(l,i) += nage_future(1,i,j)*wt(1,j)*maturity(j);
          }
        }
        if (l!=4) catch_future(l,i)   += catage_future(k,i)*wt(k);
        future_biomass(l,i) += nage_future(k,i)*wt(k);
 
      }   //end loop over future years
     }   //end loop over sex
     fspbiom_fut(l)=0.;
     for(i=styr_fut;i<=endyr_fut;i++)
       fspbiom_fut(l,i) = elem_prod(nage_future(1,i),wt(1)) * maturity;
    }   //End of loop over F's

FUNCTION compute_spr_rates
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0.initialize();
  SBF40.initialize();
  SBF35.initialize();
  SBF30.initialize();

  // Initialize the recruit (1) for each F  (F40 etc)
  for (i=1;i<=3;i++)
    Nspr(i,1)=1.;

  for (j=2;j<nages;j++)
  {
    Nspr(1,j)=Nspr(1,j-1)*exp(-1.*M(1));
    Nspr(2,j)=Nspr(2,j-1)*exp(-1.*(M(1)+F40*sel(1,j-1)/maxsel_fish));
    Nspr(3,j)=Nspr(3,j-1)*exp(-1.*(M(1)+F35*sel(1,j-1)/maxsel_fish));
    Nspr(4,j)=Nspr(4,j-1)*exp(-1.*(M(1)+F30*sel(1,j-1)/maxsel_fish));
  }

 // Now do plus group
  Nspr(1,nages)=Nspr(1,nages-1)*exp(-1.*M(1))/(1.-exp(-1.*M(1)));
  Nspr(2,nages)=Nspr(2,nages-1)*exp(-1.* (M(1)+F40*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F40*sel(1,nages)/maxsel_fish)));
  Nspr(3,nages)=Nspr(3,nages-1)*exp(-1.* (M(1)+F35*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F35*sel(1,nages)/maxsel_fish)));
  Nspr(4,nages)=Nspr(4,nages-1)*exp(-1.* (M(1)+F30*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F30*sel(1,nages)/maxsel_fish)));

  for (j=1;j<=nages;j++)
  {
   // Kill them off till april (0.25) atf spawn in winter so put in 0.0
   //         Number   ProportMat  Wt    Amount die off prior to spawning (within that year)
    SB0    += Nspr(1,j)*maturity(j)*wt(1,j)*exp(-0.0*M(1));
    SBF40  += Nspr(2,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F40*sel(1,j)/maxsel_fish));
    SBF35  += Nspr(3,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F35*sel(1,j)/maxsel_fish));
    SBF30  += Nspr(4,j)*maturity(j)*wt(1,j)*exp(-0.0*(M(1)+F30*sel(1,j)/maxsel_fish));
  }
  sprpen    = 200.*square((SBF40/SB0)-0.4);
  sprpen   += 200.*square((SBF35/SB0)-0.35);
  sprpen   += 200.*square((SBF30/SB0)-0.30);

FUNCTION Do_depend
  for (i=styr;  i<=endyr;  i++) 
  totalbiomass(i)=natage(1,i)*wt(1) + natage(2,i)*wt(2);
  obj_fun += 1.*sexr_like;             // male proportion prior, emphasis factor = 1

FUNCTION evaluate_the_objective_function
  length_like.initialize();
  age_like.initialize();   
  sel_like.initialize(); 
  fpen.initialize();
  rec_like.initialize();
  surv_like.initialize();
  catch_like.initialize();
  sexr_like.initialize();
  obj_fun.initialize(); 

  if (active(rec_dev))
  {
  length_like.initialize();
  int ii;
    //recruitment likelihood - norm2 is sum of square values   
    rec_like = norm2(rec_dev);
   
    //length likelihood
    for(k=1;k<=2;k++)
    {
      for (i=1; i <= nobs_fish; i++)
      {
        ii=yrs_fish(i);
        //fishery length likelihood fitting
          length_like(1) -= nsamples_fish(k,i)*(offset_const+obs_p_fish(k,i))*log(pred_p_fish(k,ii)+offset_const);
      }
    }
    //add the offset to the likelihood   
    length_like(1)-=offset(1);

  //survey length composition fitting

   for (i=1;i<=nsurv;i++)
   {
     for (j=1;j<=nobs_srv_length(i);j++) 
     {    
	   length_like(i+1)-=((nsamples_srv_length_fem(i,j)*(offset_const+obs_p_srv_length_fem(i,j))*log(pred_p_srv_len_fem(i,j)+offset_const))
	                      +(nsamples_srv_length_mal(i,j)*(offset_const+obs_p_srv_length_mal(i,j))*log(pred_p_srv_len_mal(i,j)+offset_const)));
	 } 
	  length_like(i+1)-=offset(i+1); 
    } 
 
//survey age composition fitting 
  for (i=1;i<=nsurv_aged;i++)
  {
	for (j=1;j<=nobs_srv_age(i);j++)
	{
   age_like(i)-=nsamples_srv_age(i,1,j)*(offset_const+obs_p_srv_age_fem(i,j))*log(pred_p_srv_age_fem(i,j)+offset_const)+
                  nsamples_srv_age(i,2,j)*(offset_const+obs_p_srv_age_mal(i,j))*log(pred_p_srv_age_mal(i,j)+offset_const); 	
	}	
	age_like(i)-=offset(i+nsurv+1);
  }
  //end of if(active (rec_dev))
  }
  // Fit to indices (lognormal)  
  //weight each years estimate by 1/(2*variance) - use cv as an approx to s.d. of log(biomass) 

  for (i=1;i<=nsurv;i++)
  {   
  surv_like(i) = norm2(elem_div(log(obs_srv(i))-log(pred_srv(i)(yrs_srv(i))),sqrt(2)*cv_srv(i)));   //yrs_srv(i) is an index here.
  }  
   
   catch_like=norm2(log(catch_bio+offset_const)-log(pred_catch+offset_const));

   // sex ratio likelihood
   sexr_like=0.5*norm2((obs_mean_sexr-pred_sexr)/obs_SD_sexr); 
 //selectivity likelihood is penalty on how smooth selectivities are   
 //here are taking the sum of squares of the second differences  

  if(active(log_selcoffs_fish))
  {  
    sel_like(1)=wt_like(1)*norm2(first_difference(first_difference(log_sel_fish(1)))); //fishery females
    sel_like(3)=wt_like(3)*norm2(first_difference(first_difference(log_sel_fish(2)))); //fishery males 

   for (j=1;j<nages;j++)
   {
    if(monot_sel==1)
    { 
        if (log_sel_fish(1,j)>log_sel_fish(1,j+1))
        sel_like(1)+=wt_like(5)*square(log_sel_fish(1,j)-log_sel_fish(1,j+1));   //monotonicity constraint fishery females
        if (log_sel_fish(2,j)>log_sel_fish(2,j+1))
        sel_like(3)+=wt_like(6)*square(log_sel_fish(2,j)-log_sel_fish(2,j+1));  //monotonicity constraing fishery males
    }
   }  

    obj_fun+=1.*sum(sel_like);    
    obj_fun+=1.*square(avgsel_fish(1));
    obj_fun+=1.*square(avgsel_fish(2));    
//obj_fun += sexr_like; note: the sex ratio like was never in the likelihood in GOA or BSAI model
  } //end if active(log_selcoffs_fish)

  // Phases less than 3, penalize High F's
    if (current_phase()<2)
    {
       //F's are low for arrowtooth changed the value to compare from .2 to .001
       //don't know if makes any difference since the penalty is reduced at the end
       fpen=10.*norm2(mfexp(fmort_dev+log_avg_fmort)-offset_const);
    }
    else
    {
      fpen=fpen_mult(1)*norm2(mfexp(fmort_dev+log_avg_fmort)-offset_const);      //0.001 for BSAI and 0.01 for GOA see fpen_mult values  
    }

    if (active(fmort_dev))
    {
      fpen+=fpen_mult(2)*norm2(fmort_dev);    //0.01 for BSAI and 1 for GOA
    }
 
  obj_fun += rec_like;
  obj_fun += like_wght(1)*length_like(1);    //this is fishery length likelihood
  for(i=2;i<=nsurv+1;i++)
   {
   obj_fun+=like_wght(6)*length_like(i);     //survey likelihood
   } 
  obj_fun+= like_wght(7)*sum(age_like); 
  for(i=1;i<=nsurv;i++)
  {
  obj_fun += like_wght(i+1)*surv_like(i);  
  }     
  obj_fun += catch_err_like*like_wght(5)*catch_like;        // large emphasis to fit observed catch 
  //obj_fun += fpen;   
  obj_fun += sprpen;     

REPORT_SECTION     
  report << "Styr" << endl<<styr<<endl;  
  report << "Endyr" << endl<<endyr<<endl; 
  report << "nages_read" <<endl<<nages_read<<endl;  
  report << "nobs_srv" <<endl<<nobs_srv<<endl; 
  report << "nobs_fish" <<endl<<nobs_fish<<endl; 
  report << "yrs_fish" <<endl<<yrs_fish<<endl; 
  report << "yrs_srv" <<endl<<yrs_srv<<endl; 
  report << "nobs_srv_length" <<endl<<nobs_srv_length<<endl; 
  report << "nobs_srv_age"<<endl<<nobs_srv_age<<endl;
  report << "assessment 1 is BSAI ATF, 2 is GOA ATF, 3 is Kamchatka"<<endl<<assess<<endl; 
  report << "Numbers_fem" << endl<<natage(1)<<endl;
  report << "Numbers_mal" << endl<<natage(2)<<endl;
  report << "Catch_est_fem"<<endl<<catage(1)<<endl;   
  report << "Catch_est_mal"<<endl<<catage(2)<<endl; 
  report << "Mort_est_fem"<<endl<<F(1)<<endl; 
  report << "Mort_est_mal"<<endl<<F(2)<<endl; 
  report << "Fishsel_fem"<<endl<<sel(1)/maxsel_fish<<endl; 
  report << "Fishsel_mal"<<endl<<sel(2)/maxsel_fish<<endl; 
  report << "Survsel_fem"<<endl<<sel_srv(1)/maxsel_srv(1)<<endl; 
  report << "Survsel_mal"<<endl<<sel_srv(2)/maxsel_srv(1)<<endl; 
  report << "Obs_srv_biomass"<<endl<<obs_srv<<endl; 
  report << "Pred_srv_biomass"<<endl<<pred_srv<<endl;
  report << "obs_srv_sd" <<endl <<obs_srv_sd<<endl;  
  report << "Obs_srv_lengthcomp_fem"<<endl<<obs_p_srv_length_fem<<endl; 
  report << "Pred_srv_lengthcomp_fem"<<endl<<pred_p_srv_len_fem<<endl; 
  report << "Obs_srv_lengthcomp_mal"<<endl<<obs_p_srv_length_mal<<endl; 
  report << "Pred_srv_lengthcomp_mal"<<endl<<pred_p_srv_len_mal<<endl;
  report << "Yrs_srv_lengthcomp"<<endl<<yrs_srv_length<<endl; 
  report << "Yrs_srv_age"<<endl<<yrs_srv_age<<endl;
  report << "Obs_srv_agecomp_fem"<<endl<<obs_p_srv_age_fem<<endl; 
  report << "Pred_srv_agecomp_fem"<<endl<<pred_p_srv_age_fem<<endl; 
  report << "Obs_srv_agecomp_mal"<<endl<<obs_p_srv_age_mal<<endl;
  report << "Pred_srv_agecomp_mal"<<endl<<pred_p_srv_age_mal<<endl;
  report << "obs_p_fish" << endl<<obs_p_fish << endl;
  report << "pred_p_fish" << endl<<pred_p_fish << endl;
  report << "Obs_catch"<<endl<<catch_bio<<endl; 
  report << "Pred_catch"<<endl<<pred_catch<<endl; 
  report << "FSB"<<endl<<fspbio<<endl; 
  report << "fmort"<<endl<<fmort<<endl;
  report << "F40"<<endl<<F40<<endl;
  report << "F35"<<endl<<F35<<endl;
  report << "F30"<<endl<<F30<<endl;
  report << "SBF40"<<endl<<SBF40<<endl;
  report << "SBF35"<<endl<<SBF35<<endl;
  report << "SBF30"<<endl<<SBF30<<endl;
  report << "SB0"<<endl<<SB0<<endl; 
  report << "surv_like"<<endl<<surv_like<<endl;
  report << "length_like_fishery"<<endl<<length_like(1)<<endl; 
  report << "length_like_survey"<<endl<<length_like(2)<<endl; 
  report << "age_like_survey1"<<endl<<age_like<<endl;  
  report << "fpen" <<endl<<fpen<<endl;   
  report << "catch_like"<<endl<<catch_like<<endl;   
  report << "rec_like"<<endl<<rec_like<<endl; 
  report << "total_likelihood"<<endl<<obj_fun<<endl; 
  report << "sexr_like"<<endl<<sexr_like<<endl;  
  report << "sel_like" <<endl<<sel_like<<endl;
  report <<"offset"<<endl<<offset<<endl; 
  report <<"explbiom"<<endl<<explbiom<<endl; 
  report <<"fspbio"<<endl<<fspbio<<endl; 
  report <<"pred_bio"<<endl<<pred_bio<<endl; 
  report <<"lenage"<<endl<<lenage<<endl;  
  report <<"fishmort"<<endl<<mfexp(log_avg_fmort+fmort_dev)<<endl;
  report <<"wt"<<endl<<wt<<endl; 
  report << "sprpen"<<endl<<sprpen<<endl;
  report << "nsamples_srv_length_fem"<<endl<<nsamples_srv_length_fem<<endl; 
  report << "nsamples_srv_length_mal"<<endl<<nsamples_srv_length_mal<<endl; 
  report << "fspbiom_fut"<<endl<<fspbiom_fut<<endl;
  report << "q_surv"<<endl<<q_surv<<endl;
  report << "M"<<endl<<M<<endl;
  report << "q_time_catchability"<<endl<<qtime<<endl;
  report << "pred_sexr"<<endl<<pred_sexr<<endl; 
  report << "obs_mean_sexr"<<endl<<obs_mean_sexr<<endl;
  report << "obs_SD_sexr"<<endl<<obs_SD_sexr<<endl; 
  report << "alpha"<<endl<<alpha<<endl; 
  report << "beta"<<endl<<beta<<endl; 
  report << "obs_srv_sd"<<endl<<obs_srv_sd<<endl; 
  report << "mean_log_rec"<<endl<<mean_log_rec<<endl; 
  report << "rec_dev"<<endl<<rec_dev<<endl; 
  report<<"Num_parameters_Estimated "<<endl<<initial_params::nvarcalc()<<endl;
    report << "obj_fun"<<endl<<obj_fun<<endl;
  report << "like_wght"<<endl<<like_wght<<endl;
  report << "nsamples_fish"<<endl<<nsamples_fish<<endl;
  report<<"nsamples_srv_age"<<endl<<nsamples_srv_age<<endl; 

  report << "SARA form for Angie Grieg" << endl;

  report << "ATF        # stock  " << endl;
  report << "BSAI       # region     (AI AK BOG BSAI EBS GOA SEO WCWYK)" << endl;
  report << "2013       # ASSESS_YEAR - year assessment is presented to the SSC" << endl;
  report << "3a         # TIER  (1a 1b 2a 2b 3a 3b 4 5 6) " << endl;
  report << "none       # TIER2  if mixed (none 1a 1b 2a 2b 3a 3b 4 5 6)" << endl;
  report << "partial    # UPDATE (new benchmark full partial)" << endl;
  report << "2          # LIFE_HIST - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "2          # ASSES_FREQ - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "5          # ASSES_LEV - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "5          # CATCH_DAT - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "3          # ABUND_DAT - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "567000     # Minimum B  Lower 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "665000     # Maximum B  Upper 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "202138     # BMSY  is equilibrium spawning biomass at MSY (Tiers 1-2) or 7/8 x B40% (Tier 3)" << endl;
  report << "ADMB       # MODEL - Required only if NMFS toolbox software used; optional otherwise " << endl;
  report << "           # VERSION - Required only if NMFS toolbox software used; optional otherwise" << endl;
  report << "2          # number of sexes  if 1 sex=ALL elseif 2 sex=(FEMALE, MALE) " << endl;
  report << "1          # number of fisheries" << endl;
  report << "1          # multiplier for recruitment, N at age, and survey number (1,1000,1000000)" << endl;
  report << "1          # recruitment age used by model or size" << endl;
  report << "1          # age+ or mmCW+ used for biomass estimate" << endl;
  report << "Single age        # Fishing mortality type such as \"Single age\" or \"exploitation rate\"" << endl;
  report << "Age model         # Fishing mortality source such as \"Model\" or \"(total catch (t))/(survey biomass (t))\"" << endl;
  report << "Age of maximum F  # Fishing mortality range such as \"Age of maximum F\"" << endl; 
  report << "#FISHERYDESC -list of fisheries (ALL TWL LGL POT FIX FOR DOM TWLJAN LGLMAY POTAUG ...)" << endl; 
  report << "ALL" << endl; 

  report <<"#FISHERYYEAR - list years used in the model " << endl;
    for (i=styr;  i<=endyr; i++)
       report << i << "	";
       report<<endl;  

  report<<"#AGE - list of ages used in the model"<<endl;
    for (i=1; i<=21;i++)
       report << i << "	";
       report<<endl;    

  report <<"#RECRUITMENT - Number of recruits by year " << endl;
    for (i=styr;  i<=endyr;  i++)
	   report  << 2*natage(1,i,1) << "	";
	   report<<endl;     

  report <<"#SPAWNBIOMASS - Spawning biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << fspbio(i) << "	";
       report<<endl;  

  report <<"#TOTALBIOMASS - Total biomass by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << natage(1,i)*wt(1) + natage(2,i)*wt(2) << "	";
       report<<endl;

  report <<"#TOTFSHRYMORT - Fishing mortality rate by year " << endl;
	for (i=styr;  i<=endyr;  i++)
	   report  << (F(1,i,13)+ F(2,i,13))/2<< "	";
	   report<<endl;

  report <<"#TOTALCATCH - Total catch by year in metric tons " << endl;
    for (i=styr;  i<=endyr;  i++)
       report  << catch_bio(i) << "	";
       report<<endl;

 report <<"#MATURITY - Maturity ratio by age (females only)" << endl;  
       for (i=1;  i<=13;  i++) 
       report  << maturity(i) <<"	";
       report<< endl; 

 report <<"#SPAWNWT - Average spawning weight (in kg) by age"<< endl; 
       report <<"0.019	0.041	0.113	0.224	0.376	0.566	0.784	1.028	1.292	1.569	1.855	2.142	2.417	2.667	2.881	3.057	3.198	3.308	3.393"<<endl;                              
       report<<endl;

 report <<"#NATMORT - Natural mortality rate for females then males"<< endl; 
 for (i=1;  i<=13;  i++) 
 report  << 0.2 <<"	";
 report<< endl;   
 for (i=1;  i<=13;  i++) 
 report  << 0.35 <<"	";
 report<< endl;

 report << "#N_AT_AGE - Estimated numbers of female (first) then male (second) fish at age " << endl;
   for (i=styr; i<=endyr;i++)
     report <<natage(1,i)<< "	";
     report<<endl;

   for (i=styr; i<=endyr;i++)
     report <<natage(2,i)<< "	";
     report<<endl;

 report <<"#FSHRY_WT_KG - Fishery weight at age (in kg) females (first) males (second), only one fishery"<< endl;   
    report <<wt(1)*1000  << "	";
    report<<endl; //1 is females        

    report <<wt(2)*1000  << "	";
    report<<endl; //2 is males


  report << "#SELECTIVITY - Estimated fishery selectivity for females (first) males (second) at age " << endl;
    for (j=1; j<=nages;j++)
      report <<" " <<sel(1,j)<< "	";
      report<<endl;

    for (j=1; j<=nages;j++)
      report <<" "  <<sel(2,j)<< "	";
      report<<endl;
 report << "#SURVEYDESC"<<endl;
 report<<"EBS_trawl_survey BS_slope_trawl_survey AI_trawl_survey"<<endl;

 report<<"SURVEYMULT"<<endl;
 report<<"1 1 1"<<endl;
 report << "#GOA_trawl_survey - Gulf of Alaska survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
    for (i=1; i<=nsurv;i++)
      report << yrs_srv(i) << "	";
      report<<endl;
    for (i=1; i<=nsurv;i++) 
      report<< obs_srv(i)<< "	";
      report<< endl;
    for (i=1; i<=nsurv;i++){     
    for (j=1; j<=nobs_srv(i);j++){
      report<< pred_srv(i,yrs_srv(i,j))<< "	"; }}
      report<< endl;

 report<<"#STOCKNOTES"<<endl;
 report<<"SAFE report indicates that this stock was not subjected to overfishing in 2012 and is neither overfished nor approaching a condition of being overfished in 2013."<<endl;
 
RUNTIME_SECTION
  maximum_function_evaluations 4000
  convergence_criteria 1e-3 1e-4 1e-7

TOP_OF_MAIN_SECTION
  arrmblsize = 20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(300);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);

