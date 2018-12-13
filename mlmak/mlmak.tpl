//////////////////////////////////////////////////////////////////
// Multispecies Age-length Amak
// Sep 15 from Sep 2 2009 
// from /sep14/H2_H3
//
// Extending part of generic model 
// (Alaska Fisheries Toolbox Assessment Model
// Version 0.1, May 2003) to handle
// multiple species, diet data, length compositions and 
// predation mortality
//
// H2 is unbounded for types I-VI
// Type VII has bounded H_3 instead of unbounded H_2
// init_vector logH_2(1,nspp_sq,PhasePredH2); // H2 is unbounded 
// (same as Sep 6 from Sep 2 2009 except ration_like)
//  !! LowerBoundH3 = -30.0;
//  !! UpperBoundH3 =  -0.000001;
//  !! LowerBoundH4 =  -0.1;
//  !! UpperBoundH4 =   20.0;
//  !! if (resp_type == 7) UpperBoundH2 = -0.00001;
//  !! if (resp_type == 6) UpperBoundH2 = 0.99999;

  // ============================
DATA_SECTION
  // ============================
  !! cout << "BEGIN DATA_SECTION" << endl << endl;
  !!CLASS ofstream mceval("mceval.dat")
  !!CLASS ofstream McFile1("SSBMCMC.Out");
  !!CLASS ofstream McFile2("RECMCMC.Out");
  !!CLASS ofstream TempFile("Temp.Out");
  int mcmcmode
  int mcflag
  int count_iters
  !! mcmcmode = 0;
  !! mcflag   = 1;
  !! count_iters = 0;

  // Declare constants 
  // ============================
  int nfsh;                                      // Total number of fleets
  int nsrv;                                      // Total number of surveys
  int tot_yr_fsh_comp                            // Total number of years of fishery age-data
  int tot_yr_srv                                 // Total number of years of survey index data
  int tot_yr_srv_age                             // Total number of years of survey age-data 
  int tot_fsh_yr                                 // Total number of catch data
  int tot_srv_yr                                 // Total number of survey data
  int isp;                                       // Counter
  int iyr;                                       // Counter
  int ifsh;                                      // Counter
  int isrv;                                      // Counter
  int icmp;                                      // Counter for report
  int i_lage;                                    // Counter for report
  int iage;                                      // Counter
  int ipnt;                                      // Pointer
  int istart;                                    // Placeholder
  int IyrPred;                                   // Year for Compute_Predation
  int rk_sp;                                     // Predation sp X sp
  int r_age;                                     // Predator age
  int rsp;                                       // Predator species
  int rksp;                                      // PredPrey index
  int ksp;                                       // Prey species
  int rln;                                       // Pred length
  int kln;                                       // Prey length
  int ru;                                        // Pred age
  int stm_yr;                                    // Year with diet length data
  int k_age;                                     // Prey age
  
  // Debug facility
  init_int DebugOut
  init_int Set_from_pin_file                      // Default is zero
  !! if (Set_from_pin_file != 0)
  !!  cout << "Set from PIN File" << endl;
  !! else
  !!  cout << "Set from DAT File" << endl;
  init_int ResetPhasesToZero                      // Set to 1 to override phases
  init_int Disc_first_phase;                      // Make First phase easier
  init_int Disc_any_phases;                       // Set this 0 to operate in last_phase
  init_int Initial_phase;                         // Initial phase (usually 1)
  init_int Terminal_phase;                        // Last phase for the model

  // predation on:off
  init_int With_Pred                              // 0:No predation; otherwise yes

  // type of predator functional response         // 1:linear; 2:Holling II; 3:Holling III
  init_int resp_type

  // Begin reading fishery data
  // =========================
  init_int      styr;                             // First model year
  init_int      endyr;                            // Last model year
  int FinalYr
  !! FinalYr = endyr;
  init_int      nspp;                             // Number of species
  init_vector   comp_type(1,nspp);                // 1=age, 2= length data for compositions
  init_ivector  oldest_age(1,nspp);               // Last age-class
  init_ivector  l_bins(1,nspp);                   // Number of length bins
  ivector       nages(1,nspp);                    // Age-range
  !!            nages = oldest_age + 1;           // Number of ages (by species)
  init_ivector  nfsh_spp(1,nspp)                  // Number of fishing fleets per species
  !! nfsh = sum(nfsh_spp);                        // Total number of fleets, all species combinations
  init_ivector  spp_fsh(1,nfsh);                  // Link between fleet and species
  ivector       ncomps_fsh(1,nfsh)                // Number of age or length compositions per fishery
  !! for(ifsh=1;ifsh<=nfsh;ifsh++)
  !!  {
  !!   isp = spp_fsh(ifsh);
  !!   if (comp_type(isp) == 1) ncomps_fsh(ifsh)=nages(isp);
  !!   else ncomps_fsh(ifsh) = l_bins(isp);
  !!  }
  !! if (DebugOut >1 ) cout << "ncomps_fsh = " << ncomps_fsh << endl;
  ivector       nages_fsh(1,nfsh);                // Number of ages per fishery
  !! for(ifsh=1;ifsh<=nfsh;ifsh++)
  !!  {
  !!   isp = spp_fsh(ifsh);
  !!   nages_fsh(ifsh) = nages(isp);
  !!  }
  !! if (DebugOut > 1) cout << "nages_fsh = " << nages_fsh << endl;
  ivector       styr_rec(1,nspp);                 // Start recruitment for fished cohorts
  ivector       endyr_sp(1,nspp);                 // Last year
  ivector       styr_sp(1,nspp);                  // First year of spawning biomass (start catches)
  int first_rec_est;
  !! first_rec_est = endyr;
  !! for(isp=1;isp<=nspp;isp++)
  !!  {
  !!   styr_rec(isp) = (styr - nages(isp)) + 1;
  !!   if (styr_rec(isp) < first_rec_est) first_rec_est = styr_rec(isp);
  !!   endyr_sp(isp) = endyr - 1;
  !!  }
  !! styr_sp  = styr_rec ;
  !! if (DebugOut > 1) cout << "styr_rec = " << styr_rec << endl;
  !! if (DebugOut > 1) cout << "styr_sp = " << styr_sp << endl;
  !! if (DebugOut > 1) cout << "endyr_sp = " << endyr_sp << endl;
  !! tot_fsh_yr = (endyr-styr+1)*nfsh;           // Total number of year*fleet*species combinations
  int nyrs;
  !! nyrs = endyr-styr+1;  

  // Read catch biomasses
  // ============================
  init_matrix catch_bio(1,nfsh,styr,endyr)            // Catch biomass
  !! if (DebugOut > 2) cout << "Catch_bio = " << endl << catch_bio << endl;
  init_3darray wt_fsh(1,nfsh,styr,endyr,1,nages_fsh)  // Weight-at-age in the catch
  !! if (DebugOut > 2) cout << "wt_fsh = " << endl << wt_fsh << endl;

  // Read catch age/length-compositions
  // ====================================
  init_ivector nyrs_fsh_comp(1,nfsh)                  // Number of years of age data (by fleet)
  !! tot_yr_fsh_comp = sum(nyrs_fsh_comp);            // Number of years of age data (over fleets)
  init_imatrix yrs_fsh_comp(1,nfsh,1,nyrs_fsh_comp)   // Years with age data
  !! if (DebugOut > 2) cout << "yrs_fsh_comp = " << endl << yrs_fsh_comp << endl;
  init_matrix nsmpl_fsh(1,nfsh,1,nyrs_fsh_comp)       // Effective sample sizes
  !! if (DebugOut > 2) cout << "nsmpl_fsh = " << endl << nsmpl_fsh << endl;
  init_3darray oc_fsh(1,nfsh,1,nyrs_fsh_comp,1,ncomps_fsh); // Fishery age or length composition data
  !! if (DebugOut > 2) cout << "oc_fsh = " << endl << oc_fsh << endl;
  vector offset_fsh(1,nfsh)                           // Fleet-specific offsets

  // Read survey-index data
  // =========================                  
  init_ivector nsrv_spp(1,nspp)                       // Number of surveys per species
  !! nsrv = sum(nsrv_spp);                            // Total number of survey-species combinations
  init_ivector spp_srv(1,nsrv);                       // Link between survey and species
  !! if (DebugOut > 2) cout << "spp_srv = " << spp_srv << endl;
  !! tot_srv_yr = (endyr-styr+1)*nsrv;                // Total number of year*survey*species combinations
  ivector       ncomps_srv(1,nsrv)                    // Number of age or length compositions per survey
  !! for(isrv=1;isrv<=nsrv;isrv++)
  !!  {
  !!   isp = spp_srv(isrv);
  !!   if (comp_type(isp) == 1) ncomps_srv(isrv)=nages(isp);
  !!   else ncomps_srv(isrv) = l_bins(isp);
  !!  }
  ivector       nages_srv(1,nsrv);                    // Number of ages per survey
  !! for(isrv=1;isrv<=nsrv;isrv++)
  !!  {
  !!   isp = spp_srv(isrv);
  !!   nages_srv(isrv) = nages(isp);
  !!  }
  init_ivector nyrs_srv(1,nsrv)                       // Number of years of index data for surveys
  !! tot_yr_srv = sum(nyrs_srv);                      // Number of years of idnex data for surveys (all species)
  init_imatrix yrs_srv(1,nsrv,1,nyrs_srv)             // Years with survey data
  init_vector mo_srv(1,nsrv)                          // Timing of the survey
  init_matrix obs_srv(1,nsrv,1,nyrs_srv)              // The survey indices
  !! if (DebugOut > 2) cout << "obs_srv = " << endl << obs_srv << endl;
  init_matrix obs_se_srv(1,nsrv,1,nyrs_srv)           // The survey standard errors

  //Read survey age data
  // ============================
  init_ivector nyrs_srv_comp(1,nsrv)                   // Number of years of age data for surveys
  !! tot_yr_srv_age = sum(nyrs_srv_comp);              // Number of years of age data for surveys (all species)
  init_imatrix yrs_srv_comp(1,nsrv,1,nyrs_srv_comp)    // Years with survey data
  init_matrix nsmpl_srv(1,nsrv,1,nyrs_srv_comp)        // Effective sample sizes
  matrix age_matrix(1,nspp,1,nages);                   // Ages for use in the model
  !! for (isp=1; isp <= nspp; isp++)
  !!  for (int a=1; a <= nages(isp); a++)
  !!   age_matrix(isp,a) = double(a);

  // Declare the vector of maximum age for each survey proportion-at-age vector
  // ============================
  ivector max_srv_age(1,tot_yr_srv_age);                  
  !! ipnt = 0;
  !! for (isrv = 1; isrv <= nsrv; isrv++)
  !!  {
  !!   for (iyr = 1; iyr <= nyrs_srv_comp(isrv); iyr++) 
  !!    { ipnt +=1; max_srv_age(ipnt) = nages(isrv); }
  !!  }
  init_3darray oc_srv(1,nsrv,1,nyrs_srv_comp,1,ncomps_srv);// Survey age or length composition data
  init_3darray  wt_srv(1,nsrv,styr,endyr,1,nages)       // Weight-at-age (survey)
  vector offset_srv(1,nsrv)                             // Offsets for the likelihood

  // Read other biological parameters
  // ============================
  init_matrix wt_pop(1,nspp,1,nages);                   // Population weight-at-age
  init_matrix maturity(1,nspp,1,nages);                 // Population maturity-at-age
  init_vector spawnmo(1,nspp);                          // Spawning month
  vector spmo_frac(1,nspp);
  !! for (isp = 1; isp<=nspp; isp++)
  !!  {
  !!   if (max(maturity(isp)) >.9 ) maturity(isp)/=2.0; // Adjust maturity
  !!   spmo_frac(isp) = (spawnmo(isp)-1)/12.;           // Adjustment for mortality before spawning
  !!  }

  // Read age-length keys
  // ============================
  init_3darray al_key(1,nspp,1,nages,1,l_bins);         // Age-length transition matrix
  matrix mean_laa(1,nspp,1,nages);                      // Mean length-at-age for predator selectivity
  !! if (DebugOut > 1) cout << "al_key = " << endl << al_key << endl;

  // end of non-predation Malmak inputs
  //================================= 
  // new predator-prey inputs:

  int nspp_sq;
  !! nspp_sq = nspp*nspp; // number pred X prey
  int nspp_sq2;
  !! nspp_sq2 = nspp*(nspp+1); // number pred X (prey + "other")
  // indices using l_bins or nages for each pred-prey combination:
  ivector  r_lens(1,nspp_sq);
  ivector  k_lens(1,nspp_sq);
  ivector  r_ages(1,nspp_sq);
  ivector  k_ages(1,nspp_sq);
  ivector  rr_lens(1,nspp_sq2);
  ivector  rr_ages(1,nspp_sq2);
  ivector  kk_ages(1,nspp_sq2);
  !!rk_sp = 0;
  !!for (rsp=1; rsp<=nspp; rsp++)
  !! for (ksp=1; ksp<=nspp; ksp++)
  !!  {
  !!   rk_sp = rk_sp+1;
  !!   r_lens(rk_sp) = l_bins(rsp);
  !!   k_lens(rk_sp) = l_bins(ksp);
  !!   r_ages(rk_sp) = nages(rsp);
  !!   k_ages(rk_sp) = nages(ksp);
  !!  } 
  !!rk_sp = 0;
  !!for (rsp=1; rsp<=nspp; rsp++)
  !! for (ksp=1; ksp<=nspp+1; ksp++)
  !!  {
  !!   rk_sp = rk_sp+1;
  !!   rr_lens(rk_sp) = l_bins(rsp);
  !!   rr_ages(rk_sp) = nages(rsp);
  !!   if (ksp <=nspp) kk_ages(rk_sp) = nages(ksp);
  !!   else kk_ages(rk_sp) = 1;
  !!  } 
  init_matrix pred_l_bin(1,nspp,1,l_bins);           // mid-points of 5cm length bins
  init_matrix omega_vB(1,nspp,1,nages);              // von_Bertalanffy ration
  init_vector omega_sigma(1,nspp);                   // von-Bertalanffy sigma (data weight)
  init_ivector nyrs_stomwts(1,nspp);                 // Number of years with predator stomach weight samples
  init_ivector nyrs_stomlns(1,nspp_sq);              // Number of years with predator stomach length samples
  init_imatrix yrs_stomwts(1,nspp,1,nyrs_stomwts);   // Years with predator stomach weights
  init_imatrix yrs_stomlns(1,nspp_sq,1,nyrs_stomlns); // Years with predator stomach lengths
  init_3darray stoms_w_N(1,nspp,1,l_bins,1,nyrs_stomwts); // Number of stomachsXpred_lenXyear with prey weights
  init_3darray  stoms_l_N(1,nspp_sq,1,r_lens,1,nyrs_stomlns);  // Number of stomachsXpred_lenXyear with prey lengths
  init_vector min_SS_w(1,nspp);                      // minimum sample size stomach weights
  init_vector max_SS_w(1,nspp);                      // maximum sample size stomach weights
  init_vector min_SS_l(1,nspp_sq);                   // minimum sample size stomach lengths
  init_vector max_SS_l(1,nspp_sq);                   // maximum sample size stomach lengths
  init_ivector i_wt_yrs_all(1,nspp_sq2);             // nyrs_stomwts for each predator species 
  init_3darray diet_w_dat(1,nspp_sq2,1,rr_lens,1,i_wt_yrs_all); // Prey weight fraction data
  init_3darray diet_l_dat(1,nspp_sq,1,r_lens,1,k_lens); // Prey length fraction data
  number offset_diet_w;                              // Scales multinomial likelihood to 1
  number offset_diet_l;                              // Scales multinomial likelihood to 1
  !! if (DebugOut > 1) cout << "nspp_sq = " << nspp_sq << endl;
  !! if (DebugOut > 1) cout << "pred_l_bin = " << endl << pred_l_bin << endl;
  !! if (DebugOut > 1) cout << "omega_vB" << endl << omega_vB << endl;
  !! if (DebugOut > 1) cout << "omega_sigma" << endl << omega_sigma << endl;
  !! if (DebugOut > 1) cout << "nyrs_stomwts" << endl << nyrs_stomwts << endl;
  !! if (DebugOut > 1) cout << "yrs_stomwts" << endl << yrs_stomwts << endl;
  !! if (DebugOut > 1) cout << "stoms_w_N" << endl << stoms_w_N << endl;
  !! if (DebugOut > 1) cout << "nyrs_stomlns" << endl << nyrs_stomlns << endl;
  !! if (DebugOut > 1) cout << "yrs_stomlns" << endl << yrs_stomlns << endl;
  !! if (DebugOut > 1) cout << "yrs_stomwts = " << endl << yrs_stomwts << endl;
  !! if (DebugOut > 1) cout << "diet_l_dat" << endl << diet_l_dat << endl;

  // Check that .DAT is finished
  // ============================
  init_int check;
  !! cout << ".DAT Check = " << check << endl << endl;

  // Now read modeling options (.CTR values)
  // ============================
  init_int SrType                                    // Stock-Recruit type: 2 Bholt, 1 Ricker
  init_vector steepnessprior(1,nspp)                 // Prior for steepness
  init_vector cvsteepnessprior(1,nspp)
  init_int phase_srec                                // Phase for steepness
  init_vector sigmarprior(1,nspp)                    // Prior for sigma-R
  vector log_sigmarprior(1,nspp)
  !! log_sigmarprior = log(sigmarprior);
  init_vector cvsigmarprior(1,nspp)
  init_int phase_sigmar                              // Phase for sigma-R
  init_ivector styr_rec_est(1,nspp)                  // Start year for recruitment estimate
  init_ivector endyr_rec_est(1,nspp)                 // End year for recruitment estimate
  ivector nrecs_est(1,nspp);                         // Number of estimated recruitments
  int n_est_recs;
  !! n_est_recs = 0;
  !! for(isp=1; isp<=nspp; isp++)
  !!  {
  !!   nrecs_est(isp) = endyr_rec_est(isp)-styr_rec_est(isp)+1;
  !!   n_est_recs += endyr_rec_est(isp)- styr_rec(isp)+1;
  !!  cout << "styr_rec_est(" << isp <<") = " <<styr_rec_est(isp) << endl;
  !!  }
  !! cout << "n_est_recs = " << n_est_recs << endl;
  init_vector natmortprior(1,nspp)                   // Prior for M
  init_vector cvnatmortprior(1,nspp)
  init_ivector natmortphase2(1,nspp)                 // Allows some Ms to be estimated and others not
  int NEstNat                                        // Number of species for which M is estimated
  !! NEstNat = 0;
  !! for (isp=1;isp<=nspp;isp++) if (natmortphase2(isp) > 0) NEstNat += 1;
  init_vector qprior(1,nsrv)                         // Prior for survey-q
  vector log_qprior(1,nsrv)
  !! log_qprior = log(qprior);
  init_matrix cvqprior(1,nspp,1,nsrv_spp)
  init_ivector phase_q(1,nsrv)                       // Phase for survey-q
  init_ivector q_age_min(1,nsrv)                     // Minimum age for variable survey selectivity
  init_ivector q_age_max(1,nsrv)                     // Maximum age for variable survey selectivity
  init_vector cv_catchbiomass(1,nspp)                // Penalty on catch biomass
  init_vector sd_ration(1,nspp)                      // Denominator in ration_like
  init_int phase_M                                   // Phase for MEst
  init_int phase_Rzero;                              // Phase for Rzero
  init_int phase_fmort;                              // Phase for fmortality
  init_int phase_fmort1;                             // Phase for fmortality (average)
  init_int phase_LogRec;                             // Phase for mean_log_rec
  init_int phase_RecDev;                             // Phase for rec_dev
  init_int phase_SelFshCoff;                         // Phase for log_selcoffs_fsh
  init_int phase_SelSrvCoff;                         // Phase for log_selcoffs_srv
  int phase_SelSrvCoff2;
  init_int PhasePred1;                               // Phase for predator selectivity log_gam_a, log_gam_b
  init_int PhasePred2;                               // Phase for predator selectivity H_1, H_2, H_3
  init_int PhasePred3;                               // Phase for predator selectivity Q_other
  int PhasePredH1a;
  int PhasePredH2;
  int PhasePredH3;
  int PhasePredH4;
  !! { PhasePredH1a= -1; PhasePredH2 = -1; PhasePredH3 = -1; PhasePredH4 = -1; }
  !! if (phase_M > 0)
  !!  { phase_M = phase_M - Initial_phase+1; if (phase_M < 1) phase_M = 1; }
  !! if (phase_Rzero > 0) 
  !!  { phase_Rzero = phase_Rzero - Initial_phase+1; if (phase_Rzero < 1) phase_Rzero = 1; }
  !! if (phase_fmort > 0)
  !!  { phase_fmort = phase_fmort - Initial_phase+1; if (phase_fmort < 1) phase_fmort = 1; }
  !! if (phase_fmort1 > 0)
  !!  { phase_fmort1 = phase_fmort1 - Initial_phase+1; if (phase_fmort1 < 1) phase_fmort1 = 1; }
  !! if (phase_LogRec > 0)
  !!  { phase_LogRec = phase_LogRec - Initial_phase+1; if (phase_LogRec < 1) phase_LogRec = 1; }
  !! if (phase_RecDev > 0)
  !!  { phase_RecDev = phase_RecDev - Initial_phase+1; if (phase_RecDev < 1) phase_RecDev = 1; }
  !! if (phase_SelFshCoff > 0)
  !!  { phase_SelFshCoff = phase_SelFshCoff - Initial_phase+1; if (phase_SelFshCoff < 1) phase_SelFshCoff = 1; }
  !! if (phase_SelSrvCoff > 0)
  !!  { phase_SelSrvCoff = phase_SelSrvCoff - Initial_phase+1; if (phase_SelSrvCoff < 1) phase_SelSrvCoff = 1; }
  !! if (PhasePred1 > 0)
  !!  { PhasePred1 = PhasePred1 - Initial_phase+1; if (PhasePred1 < 1) PhasePred1 = 1; }
  !! if (PhasePred2 > 0)
  !!  { PhasePred2 = PhasePred2 - Initial_phase+1; if (PhasePred2 < 1) PhasePred2 = 1; }
  !! if (PhasePred3 > 0)
  !!  { PhasePred2 = PhasePred3 - Initial_phase+1; if (PhasePred3 < 1) PhasePred3 = 1; }
  
  // Deal with phasePreds
  !! if (resp_type > 0) PhasePredH1a = PhasePred2;
  !! if (1 < resp_type < 7) PhasePredH2 = PhasePred2 + 1;
  !! if (resp_type > 5)  PhasePredH3 = PhasePred2 + 1;
  !! if (resp_type == 3 || resp_type == 6) PhasePredH4 = PhasePred2 + 2;
  !! if (resp_type == 4 || resp_type == 5) PhasePredH3 = PhasePred2 + 2;
  !! phase_SelSrvCoff2 = phase_SelSrvCoff;
  !! cout << PhasePred1 << " " << PhasePredH1a << " " << PhasePredH2 << " " << PhasePredH3 << " " << PhasePredH4 << endl;


  vector catchbiomass_pen(1,nspp)                    // Convert cv_catchbiomass to penalty
  !! for(isp=1; isp<=nspp; isp++)
  !!  catchbiomass_pen(isp)= 1.0/(2*cv_catchbiomass(isp)*cv_catchbiomass(isp));

  // Declare Fishery Selectivies
  // ============================
  init_ivector fsh_sel_opt(1,nfsh)                   // Options for fishery selectivity
  init_vector nselages_in_fsh(1,nfsh)                // Number of age classes with selectivities by fishery
  init_ivector phase_sel_fsh(1,nfsh)                 // Phase for selectivity parameters
  init_vector curv_pen_fsh(1,nfsh)                   // Penalty on the curvature of selectivity
  init_vector seldec_pen_fsh(1,nfsh)                 // Penalty on declining selectivity with age
  ivector seldecage(1,nfsh);                         // midpoint of selected ages
  !!for(ifsh=1; ifsh<=nfsh; ifsh++)  
  !! { 
  !!  isp =  spp_fsh(ifsh);
  !!  seldecage(ifsh) = int(nages(isp)/2);
  !! }  

  init_matrix sel_change_in_fsh(1,nfsh,styr,endyr)   // Changes in fishery selectivity 
  ivector n_sel_ch_fsh(1,nfsh);                      // Number of years of fishery selectivity changes
  imatrix yrs_sel_ch_tmp(1,nfsh,1,nyrs);             // The years of fishery selectivity changes
  ivector phase_selcoff_fsh(1,nfsh);                 // Phase for fishery selectivity
  !! phase_selcoff_fsh = phase_sel_fsh;
  int fsh_sel_opt_tmp;                               // Switch for type of fishery selectivity

  // Declare Survey Selectivies 
  // ==========================
  init_ivector srv_sel_opt(1,nsrv)                   // Options for survey selectivity
  init_matrix sel_change_in_srv(1,nsrv,styr,endyr)   // Changes survey selectivity
  init_ivector phase_sel_srv(1,nsrv)                 // Phase for survey selectivity
  vector sel_slp_in_srv(1,nsrv)                      // Survey selectivity slope
  vector sel_inf_in_srv(1,nsrv)                      // Survey selectivity inflection
  ivector nselages_in_srv(1,nsrv)                    //  Number of age classes with selectivities by survey
  vector curv_pen_srv(1,nsrv)                        // Penalty of curvature of survey selectivity
  vector seldec_pen_srv(1,nsrv)                      // Penalty on declining selectivity with age
  vector logsel_slp_in_srv(1,nsrv);                  // Slope of logistic selectivity for survey
  vector sel_dinf_in_srv(1,nsrv);                    // maximum logistic selectivity for survey
  ivector n_sel_ch_srv(1,nsrv);                      // Number of years of survey selectivity changes
  imatrix yrs_sel_ch_tsrv(1,nsrv,1,nyrs);            // The years of survey selectivity changes
  ivector phase_selcoff_srv(1,nsrv);                 // Phase for survey selectivity
  !! phase_selcoff_srv = phase_sel_srv;               
  int srv_sel_opt_tmp;                               // Switch for form of survey selectivity

  // Initialize nselages for fisheries and surveys
  // ============================
  !! nselages_in_srv = nages-1;

  // Read survey selectivity-option dependent variables
  // ============================
 LOCAL_CALCS 
  for(isrv = 1; isrv <= nsrv; isrv++)
   {
    curv_pen_srv(isrv) = 0;    // initialize with zeroes for couts
    seldec_pen_srv(isrv) = 0;
    if(srv_sel_opt(isrv) == 1) // mackerel
     {
      *(ad_comm::global_datafile) >> nselages_in_srv(isrv);
      *(ad_comm::global_datafile) >> curv_pen_srv(isrv);
      *(ad_comm::global_datafile) >> seldec_pen_srv(isrv);
      phase_selcoff_srv(isrv) = phase_sel_srv(isrv);
      logsel_slp_in_srv(isrv) = 0.0;
      sel_inf_in_srv(isrv)    = 0.0;
      sel_dinf_in_srv(isrv)    = 0.0;
     }
    if (srv_sel_opt(isrv)==2) // pollock
     {
      *(ad_comm::global_datafile) >> sel_slp_in_srv(isrv);
      *(ad_comm::global_datafile) >> sel_inf_in_srv(isrv);
      phase_selcoff_srv(isrv) = -1;
      logsel_slp_in_srv(isrv) = log(sel_slp_in_srv(isrv)) ;
      sel_dinf_in_srv(isrv)    = 0.0;
     }
    if (phase_selcoff_srv(isrv)>0) curv_pen_srv(isrv) = 1./ (square(curv_pen_srv(isrv))*2);
   }
 END_CALCS  //  ===================================== 

  !! if (DebugOut > 1) cout << "srv_sel_opt = " << srv_sel_opt << endl;
  !! if (DebugOut > 1) cout << "phase_sel_srv = " << phase_sel_srv << endl;
  !! if (DebugOut > 1) cout << "curv_pen_srv = " << curv_pen_srv << endl;
  !! if (DebugOut > 1) cout << "seldec_pen_srv = "  << seldec_pen_srv << endl;
  !! if (DebugOut > 1) cout << "nselages_in_srv = " << nselages_in_srv << endl;

  // .DAT and .CTR values have now been read in.
  // Next declare variables whose values will be
  // assigned in the PRELIMINARY_CALCS_SECTION.
  // ============================
  number Steepness_UB;                           // Upper bound for Steepness
  !! Steepness_UB =   .999;                               
  vector R_guess(1,nspp);                        // Initial values for R0

  // ============================
  // Specify how many changing selectivities there are
  // Used to calculate number of selected age-lengths (nselages)
 LOCAL_CALCS

  // fishery selectivity
  // ============================
  int j;
  for (ifsh=1; ifsh<=nfsh; ifsh++)
   {
    sel_change_in_fsh(ifsh,styr) = 1.0;
    n_sel_ch_fsh(ifsh) = 0;
    j = 1;
    yrs_sel_ch_tmp(ifsh,j) = styr;
    for (iyr=styr+1; iyr<=endyr; iyr++)           // Count the number of changes in selectivity
     {
      if (sel_change_in_fsh(ifsh,iyr) > 0)
       { j++; yrs_sel_ch_tmp(ifsh,j) = iyr; }
     }
    n_sel_ch_fsh(ifsh) = j;
   }
  if (DebugOut > 1) cout << "sel_change_in_fsh = " << endl << sel_change_in_fsh << endl;
  if (DebugOut > 1) cout << "yrs_sel_ch_tmp = " << endl << yrs_sel_ch_tmp << endl;
  if (DebugOut > 1) cout << "n_sel_ch_fsh = " << endl << n_sel_ch_fsh << endl;

  // survey selectivity
  // ============================
  for (isrv=1; isrv<=nspp; isrv++)
   {
    sel_change_in_srv(isrv,styr) = 1.0;
    n_sel_ch_srv(isrv) = 0;
    j = 1;
    yrs_sel_ch_tsrv(isrv,j) = styr;
    for (iyr=styr+1; iyr<=endyr; iyr++)           // Count the number of changes in selectivity
     {
      if (sel_change_in_srv(isrv,iyr) > 0)
       { j++; yrs_sel_ch_tsrv(isrv,j) = iyr; }
     }
    n_sel_ch_srv(isrv) = j;
   }
  if (DebugOut > 1) cout << "sel_change_in_srv = " << endl << sel_change_in_srv << endl;
  if (DebugOut > 1) cout << "yrs_sel_ch_tsrv = " << endl << yrs_sel_ch_tsrv << endl;
  if (DebugOut > 1) cout << "n_sel_ch_srv = " << endl << n_sel_ch_srv << endl;
 END_CALCS  //=====================================

  // Declare years and ages of changing selectivities (Fishery and survey)
  // ==============================================
  imatrix  yrs_sel_ch_fsh(1,nfsh,1,n_sel_ch_fsh);
  imatrix  nselages_fsh(1,nfsh,1,n_sel_ch_fsh);
  imatrix  yrs_sel_ch_srv(1,nsrv,1,n_sel_ch_srv);
  imatrix  nselages_srv(1,nsrv,1,n_sel_ch_srv);  
  // Calculate the number of ages for each fleet / survey
  // ==============================================
 LOCAL_CALCS 
  for (ifsh=1;ifsh<=nfsh;ifsh++) nselages_fsh(ifsh) = nselages_in_fsh(ifsh);
  for (isrv=1;isrv<=nsrv;isrv++) nselages_srv(isrv) = nselages_in_srv(isrv);
 END_CALCS  //=====================================

  ivector endyr_all(1,nspp);
  !! endyr_all = endyr_sp + 1; 
  !! cout << endl << "END DATA_SECTION" << endl << endl;
  
  // Check that .DAT is finished
  // ============================
  init_int check2;
  !! cout << ".CTR Check = " << check2 << endl << endl;
  
  int styr_pred
  !! styr_pred = 1960; // 1945
  
  int nyrs_pred;
  !! nyrs_pred = endyr-styr_pred+1;  

  int PhaseDummy;
  !! PhaseDummy = -1;
  !! if (With_Pred == 0)
  !! {PhasePred1 = -3; PhasePred2 = -1; PhasePred3 = -1; }
  !! if (ResetPhasesToZero==1) 
  !!  {
  !!   cout << "Resetting all phases" << endl;
  !!   PhaseDummy = 1;
  !!   PhasePred1 = -99;
  !!   PhasePred2 = -99;
  !!   PhasePred3 = -99;
  !!   PhasePredH1a= -99; // not in NRM tpl -dhk april 28 09
  !!   PhasePredH2 = -99;
  !!   PhasePredH3 = -99;
  !!   PhasePredH4 = -99;
  !!   phase_M = -99;
  !!   phase_srec = -99;
  !!   phase_Rzero = -99;
  !!   phase_fmort = -99;
  !!   phase_fmort1 = -99;
  !!   for (ifsh=1;ifsh<=nfsh;ifsh++) phase_selcoff_fsh(ifsh) = -99;
  !!   for (isp=1;isp<=nspp;isp++) phase_selcoff_srv(isp) = -99;
  !!   phase_LogRec = -99;
  !!   phase_RecDev = -99;
  !!   phase_SelFshCoff = -99;
  !!   phase_SelSrvCoff = -99;
  !!   phase_SelSrvCoff2 = -99; // not in NRM tpl -dhk april 28 09
  !!  }
  
  // Phase tidy up
  !! if (phase_SelFshCoff == -99)
  !!  for (ifsh=1;ifsh<=nfsh;ifsh++) phase_selcoff_fsh(ifsh) = -99;

  //number LowerBoundH2;
  //number UpperBoundH2;
  number LowerBoundH3;
  number UpperBoundH3;
  number LowerBoundH4;
  number UpperBoundH4;
  //!! LowerBoundH2 = -30.0;
  //!! UpperBoundH2 =   20.0;
  !! LowerBoundH3 = -30.0;
  !! UpperBoundH3 =  -0.000001;
  !! LowerBoundH4 =  -0.1;
  !! UpperBoundH4 =   20.0;
  !! if (resp_type == 7) UpperBoundH3 = -0.00001;
  //!! if (resp_type == 6) UpperBoundH3 = 0.99999;

   
// =========================================================================
// =========================================================================

PARAMETER_SECTION

  !! cout << "BEGIN PARAMETER_SECTION" << endl << endl;

  // Fishery selectivity parameters combined over species
  // ============================
  !! int Nselfshpars = 0;
  !! for (ifsh=1;ifsh<=nfsh;ifsh++)
  !!  for (iyr=1;iyr<=n_sel_ch_fsh(ifsh);iyr++)
  !!   for (iage=1;iage<=nselages_fsh(ifsh,iyr);iage++)
  !!    Nselfshpars += 1;
  !! cout << "Nselfshpars = " << Nselfshpars << endl;

  // Survey selectivity parameters combined over species 
  // =============================
  !! int Nselsrvpars = 0;
  !! int Nselsrvlogs = 0;
  !! for (isrv=1;isrv<=nsrv;isrv++)
  !!  {
  !!   if (srv_sel_opt(isrv) == 1)
  !!    for (iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
  !!     for (iage=1;iage<=nselages_srv(isrv,iyr);iage++)
  !!      Nselsrvpars += 1;
  !!   if (srv_sel_opt(isrv) == 2)
  !!      Nselsrvlogs += 1;
  !!  }
  !! cout << "Nselsrvpars = " << Nselsrvpars << endl;
  !! cout << "Nselsrvlogs = " << Nselsrvlogs << endl;
  !! if (Nselsrvlogs == 0)
  !!  { Nselsrvlogs = -1; phase_SelSrvCoff2 = -1; }
  
  !! cout << "PHASE VALUES: " << endl;
  !! cout << "phase_M = " << phase_M << endl;
  !! cout << "phasePred1 = " << PhasePred1 << endl;
  !! cout << "phasePred2 = " << PhasePred2 << endl;
  !! cout << "phasePred3 = " << PhasePred3 << endl;
  !! cout << "phasePredH1a = " << PhasePredH1a << endl;
  !! cout << "phasePredH2 = " << PhasePredH2 << endl;
  !! cout << "phasePredH3 = " << PhasePredH3 << endl;
  !! cout << "phasePredH4 = " << PhasePredH4 << endl;
  !! cout << "phase_LogRec = " << phase_LogRec << endl;
  !! cout << "phase_srec = " << phase_srec << endl;
  !! cout << "phase_Rzero = " << phase_Rzero << endl;
  !! cout << "phase_RecDev = " << phase_RecDev << endl;
  !! cout << "phase_sigmar = " << phase_sigmar << endl;
  !! cout << "phase_SelFshCoff = " << phase_SelFshCoff << endl;
  !! cout << "phase_q = " << phase_q << endl;
  !! cout << "phase_SelSrvCoff = " << phase_SelSrvCoff << endl;
  !! cout << "phase_SelSrvCoff2 = " << phase_SelSrvCoff2 << endl;
  !! cout << "phase_fmort = " << phase_fmort << endl << endl;
  !! cout << "phase_fmort1 = " << phase_fmort1 << endl << endl;
  !! cout << "phase_selcoff_fsh = " << phase_selcoff_fsh << endl;
  !! cout << "phase_selcoff_srv = " << phase_selcoff_srv << endl;
  
  // Count how many F_devs are needed
  !! int NFdevs = 0;
  !! for (ifsh=1;ifsh<=nfsh;ifsh++)
  !!  for (iyr=styr;iyr<=endyr;iyr++)
  !!   if (catch_bio(ifsh,iyr) > 10e-24) NFdevs ++; // 0 in NRM tpl -dhk apr 28 09
  !! cout << "NFdevs = " << NFdevs << endl;

  // Estimated parameters by species, fishery, or survey
  // =============================
  init_bounded_vector MEst(1,NEstNat,0.02,0.8,phase_M)             // Natural mortality (phase_M: 5, -4)
  init_bounded_vector log_gam_a(1,nspp,1.0e-10,19.9,PhasePred1);   // Predator selectivity
  init_bounded_vector log_gam_b(1,nspp,-5.2,10,PhasePred1);        // Predator selectivity
  init_matrix Q_other_est(1,nspp,1,nages,PhasePred3);              // Other prey biomass in predator diets
  
  init_vector logH_1(1,nspp_sq2,PhasePred2);                       // Predation functional form
  init_vector logH_1a(1,nspp,PhasePredH1a);                        // Age adjustment to H_1
  init_vector logH_1b(1,nspp,PhasePredH1a);                        // Age adjustment to H_1
  !! 
  init_vector logH_2(1,nspp_sq,PhasePredH2);                                    // Predation functional form
  init_bounded_vector logH_3(1,nspp_sq,LowerBoundH3,UpperBoundH3,PhasePredH3);  // Predation functional form
  init_bounded_vector    H_4(1,nspp_sq,LowerBoundH4,UpperBoundH4,PhasePredH4);  // Predation functional form

  vector H_1(1,nspp_sq2);
  vector H_1a(1,nspp);
  vector H_1b(1,nspp);
  vector H_2(1,nspp_sq);
  vector H_3(1,nspp_sq);
  vector gam_a(1,nspp);
  vector gam_b(1,nspp);
  vector M(1,nspp);
  init_vector mean_log_rec(1,nspp,phase_LogRec);                     // Log(mean recruitment) (phase = 1)
  init_bounded_vector steepness(1,nspp,0.21,Steepness_UB,phase_srec) // Steepness (phase_srec: -6, -6)
  init_bounded_vector log_Rzero(1,nspp,-100,100,phase_Rzero)         // Log(R0) (phase_Rzero = 2) : init_vector in NRM -dhk apr 28 09
  init_bounded_vector rec_dev(1,n_est_recs,-15,2,phase_RecDev);      // All species combined (rec devs, phase = 2)
  init_vector log_sigmar(1,nspp,phase_sigmar);                       // Sigma(R) (phase_sigmar = -5, -4)
  init_vector log_selcoffs_fsh(1,Nselfshpars,phase_SelFshCoff)       // Log (selectivity coefficients) (phase_selcoff_fsh: 3, 4)
  init_number_vector log_q_srv(1,nsrv,phase_q)                       // Survey-q (phase_q: -6, -4)
  init_vector log_selcoffs_srv(1,Nselsrvpars,phase_SelSrvCoff);      // Survey selectivity coefficients (phase_selcoff_srv: -1, 5)
  init_vector logsel_slope_srv_par(1,Nselsrvlogs,phase_SelSrvCoff2)  // Selectivity slope ( 3, -1)
  init_vector sel50_srv_par(1,Nselsrvlogs,phase_SelSrvCoff2)         // Length-at-50%-selectivity ( 3, -1)
  matrix logsel_slope_srv(1,nsrv,1,n_sel_ch_srv)                     // Selectivity slope ( 3, -1)
  matrix sel50_srv(1,nsrv,1,n_sel_ch_srv)                            // Length-at-50%-selectivity ( 3, -1)
  init_bounded_vector fmort_dev_est(1,NFdevs,-12,8,phase_fmort)      // Fishing mortality deviations
  init_bounded_vector log_avg_fmort(1,nfsh,-10,2,phase_fmort1)       // Average F (phase = 1)
  init_number dummy(PhaseDummy);
  init_number dummy2(Terminal_phase);

  // Derived parameters
  // =============================
  3darray natage(1,nspp,styr_pred,endyr,1,nages)           // Numbers-at-age
  matrix Sp_Biom(1,nspp,styr_sp,endyr_all)                 // Spawning biomass
  matrix pred_rec(1,nspp,styr_rec,endyr_all)               // Recruitment from s-r relationship
  matrix mod_rec(1,nspp,styr_rec,endyr_all)                // Recruitment as estimated by model
  3darray Z(1,nspp,styr_pred,endyr,1,nages)                // Total mortality
  3darray F(1,nfsh,styr,endyr,1,nages_fsh);                // Fishing mortality
  3darray S(1,nspp,styr_pred,endyr,1,nages)                // Survival from total mortality
  3darray catage(1,nfsh,styr,endyr,1,nages_fsh)            // Catch-at-age
  vector surv(1,nspp)                                      // Survival from natural mortality
  vector natmort(1,nspp)                                   // Natural mortality
  matrix rec_dev_spp(1,nspp,styr_rec,endyr_all);           // Recruitment devs unpacked by spp
  matrix fmort_dev(1,nfsh,styr,endyr)                      // Fishing mortality deviations
  matrix Fmort(1,nspp,styr,endyr);                         // Annual total Fmort
  vector m_sigmarsq(1,nspp)                                // Sigma(R)**2
  vector m_sigmar(1,nspp)                                  // Signa(R)
  vector sigmarsq(1,nspp)                                  // Sigma(R)**2
  vector sigmar(1,nspp)                                    // Sigma(R)
  vector alpha(1,nspp)                                     // S-R alpha
  vector beta(1,nspp)                                      // S-R beta
  vector Bzero(1,nspp)                                     // B0
  vector Rzero(1,nspp)                                     // R0
  vector phizero(1,nspp)                                   // Phi(0)
  vector avg_rec_dev(1,nspp)                               // Average rec dev
  matrix avgsel_fsh(1,nfsh,1,n_sel_ch_fsh);                // Average selectivity
  matrix sel_slope_fsh(1,nfsh,1,n_sel_ch_fsh)              // Selectivity 
  3darray log_sel_fsh(1,nfsh,styr,endyr,1,nages_fsh)       // Log(selectivity)
  3darray sel_fsh(1,nfsh,styr,endyr,1,nages_fsh)           // Selectivity

  // Fishery composition and catch parameters
  // ============================
  3darray eac_fsh(1,nfsh,1,nyrs_fsh_comp,1,nages_fsh);     // Predicted catch-at-age
  3darray ec_fsh(1,nfsh,1,nyrs_fsh_comp,1,ncomps_fsh)      // Predicted catch-at-length or age (depends on composition data)
  matrix pred_catch(1,nfsh,styr,endyr)                     // Total catch

  // Survey composition and selectivity parameters
  // ============================
  matrix sel_slope_srv(1,nsrv,1,n_sel_ch_srv)              // Selectivity slope
  3darray log_sel_srv(1,nsrv,styr,endyr,1,nages)           // Log-survey selectivity
  3darray sel_srv(1,nsrv,styr,endyr,1,nages)               // Survey selectivity
  matrix avgsel_srv(1,nsrv,1,n_sel_ch_srv);                // Average selectivity
  matrix pred_srv(1,nsrv,styr,endyr)                       // Predicted index
  3darray eac_srv(1,nsrv,1,nyrs_srv_comp,1,nages_srv)      // Predicted survey numbers at age
  3darray ec_srv(1,nsrv,1,nyrs_srv_comp,1,ncomps_srv)      // Predicted catch-at-length or age (depends on composition data)
  vector q_srv(1,nsrv)                                     // Survey q

  // Predation variables
  // ===========================
  3darray  gam_ua(1,nspp_sq,1,r_ages,1,k_ages)             // gamma selectivity of predator age u on prey age a
  matrix   N_pred_eq(1,nspp,1,nages)                       // Effective numbers of predators for each age of prey (styr_pred)
  matrix   N_prey_eq(1,nspp,1,nages)                       // Effective numbers of prey for each age of predator
  matrix   N_pred_yr(1,nspp,1,nages)                       // Effective numbers of predators for each age of prey (IyrPred)
  matrix   N_prey_yr(1,nspp,1,nages)                       // Effective numbers of prey for each age of predator
  3darray  N_pred_eqs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_pred_eq for all yrs
  3darray  N_prey_eqs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_prey_eq for all yrs
  3darray  N_pred_yrs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_pred_yr for all yrs
  3darray  N_prey_yrs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_prey_yr for all yrs
  3darray  Pred_r(1,nspp,styr_pred,FinalYr,1,nages)        // save Pred_ratio values
  3darray  Prey_r(1,nspp,styr_pred,FinalYr,1,nages)        // save Prey_ratio values
  4darray  pred_resp(1,nspp_sq2,styr_pred,FinalYr,1,rr_ages,1,kk_ages) // Predator functional response
  3darray  Pmort_ua(1,nspp_sq,styr_pred,FinalYr,1,k_ages)  // Predation mortality on prey age a by all predators age u
  4darray  Vmort_ua(1,nspp_sq,1,r_ages,styr_pred,FinalYr,1,k_ages) // Predation mortality on prey age a by single predator age u
  4darray  eaten_la(1,nspp_sq,1,r_lens,styr_pred,FinalYr,1,k_ages) // Number of prey of age a eaten by predator length l
  4darray  eaten_ua(1,nspp_sq,1,r_ages,styr_pred,FinalYr,1,k_ages) // Number of prey of age a eaten by predator age u
  3darray  Q_mass_l(1,nspp_sq2,styr_pred,endyr,1,rr_lens)  // Mass of prey consumed by length l of predator
  3darray  Q_mass_u(1,nspp_sq2,styr_pred,endyr,1,rr_ages)  // Mass of prey consumed by age u of predator
  3darray  omega_hat(1,nspp,styr_pred,endyr,1,nages)       // Daily ration by predator age each year
  matrix   omega_hat_ave(1,nspp,1,nages)                   // Daily ration by predator age averaged over years
  3darray  Q_hat(1,nspp_sq2,styr_pred,endyr,1,rr_lens)     // Fraction for each prey type of total mass eaten by predator length 
  matrix   Q_other_u(1,nspp,1,nages)                       // Mass of other consumed by age u of predator
  3darray  T_hat(1,nspp_sq,1,r_lens,1,k_lens)              // Fraction of prey of length m in predator of length l
  matrix   Zcurr(1,nspp,1,nages);
  matrix   Zlast(1,nspp,1,nages);

  // Store the number of selected ages
  // ============================
  !! for (ifsh=1;ifsh<=nfsh;ifsh++) nselages_fsh(ifsh)=nselages_in_fsh(ifsh); // Sets all elements of a vector to one scalar value...
  !! for (isrv=1;isrv<=nsrv;isrv++) nselages_srv(isrv)=nselages_in_srv(isrv); // Sets all elements of a vector to one scalar value...
  
  // Likelihood value names         
  // ============================
  vector sigma(1,nspp)                                     // Sigma
  matrix rec_like(1,nspp,1,4)                              // Recruitment 
  vector catch_like(1,nfsh)                                // Catch 
  vector age_like_fsh(1,nfsh)                              // Fishery age-composition
  vector age_like_srv(1,nsrv)                              // Survey age-compostion
  matrix sel_like_fsh(1,nfsh,1,4)                          // Fishery-selectivity penalty
  matrix sel_like_srv(1,nsrv,1,4)                          // Survey-selectivity penalty
  vector surv_like(1,nsrv)                                 // Survey index component
  matrix fpen(1,nspp,1,6)                                  // F-penalty
  vector post_priors(1,4)                                  // Compute_priors penalty
  vector post_priors_srvq(1,nsrv)                          // Compute_priors penalty
  vector ration_like(1,nspp);                              // Ration likelihoods
  number diet_like1;                                       // Prey weights 
  number diet_like2;                                       // Prey lengths 
  vector Zlast_pen(1,nspp);                                // Penalty on stability
  vector obj_comps(1,16)                                   // Values of likelihood components
  number penal_rec_dev;
  //number Misc_penal;
  vector ration_pen(1,nspp);
  number mean_ohat;

  objective_function_value obj_fun;                        // Objective function
  number ObjTemp;                                          // Objective function value from previous iteration
  number SSLow;                                            // Maximum value for ObjTemp

  !! cout << "END OF PARAMETER_SECTION" << endl << endl;

 sdreport_matrix SSBOut(1,nspp,first_rec_est,endyr);

// ========================================================================
PRELIMINARY_CALCS_SECTION
// ========================================================================
  double btmp,ctmp,TotN,LogH1Low;
  int nagestmp, iyrs,Itot,II,iage;
  
  cout << "BEGIN PRELIMINARY_CALCS_SECTION" << endl << endl;

  // Penalty on the curvature of fishery selectivity (only used if opt_fsh_sel=1)
  curv_pen_fsh = 1./ (square(curv_pen_fsh)*2);         

  // Compute an initial guess for Rzero based on exploitation 
  // ====================================================
   for (isp=1; isp<=nspp; isp++)
    {
     nagestmp = nages(isp);
     btmp = 0.0; ctmp=   0.0;
     dvector ntmp(1,nagestmp);
     ntmp(1) = 1.0/exp(-natmortprior(isp)-.05);
     for (int a = 2; a <= nagestmp; a++)
      ntmp(a) = ntmp(a-1)*exp(-natmortprior(isp)-.05);
     btmp = wt_pop(isp) * ntmp;
     ctmp = mean(catch_bio(isp));
     R_guess(isp) = log((ctmp/.05 )/btmp/exp(-natmortprior(isp)) ) ;
    }
  
  cout << "R_guess = " << R_guess << endl;

  // Compute fishery offsets to be used in FUNCTION Age_Like (PRELIMINARY_CALCS_SECTION)
  // ===================================================================================
  offset_fsh.initialize();
  for (ifsh = 1; ifsh <= nfsh; ifsh++)
   for (iyr = 1; iyr <= nyrs_fsh_comp(ifsh); iyr++)
    {
     oc_fsh(ifsh,iyr) /= sum(oc_fsh(ifsh,iyr));
     offset_fsh(ifsh) -= nsmpl_fsh(ifsh,iyr)*(oc_fsh(ifsh,iyr) + 0.001) * log(oc_fsh(ifsh,iyr) + 0.001 ) ;   
    }

  // Compute survey offsets to be used in FUNCTION Age_Like (PRELIMINARY_CALCS_SECTION)
  // ==================================================================================
  offset_srv.initialize();
  for (isrv = 1; isrv <= nsrv; isrv++)
   for (iyr = 1; iyr <= nyrs_srv_comp(isrv); iyr++)
    {
     oc_srv(isrv,iyr) /= sum(oc_srv(isrv,iyr));
     offset_srv(isrv) -= nsmpl_srv(isrv,iyr)*(oc_srv(isrv,iyr) + 0.001) * log(oc_srv(isrv,iyr) + 0.001 ) ;   
    }

  // Find mean length-at-age for gamma selectivity
  for (rsp=1;rsp<=nspp;rsp++)
   for (iage=1;iage<=nages(rsp);iage++)
    {
     mean_laa(rsp,iage) = 0;
     for (rln=1;rln<=l_bins(rsp);rln++)
      mean_laa(rsp,iage) += al_key(rsp,iage,rln)*pred_l_bin(rsp,rln);
    } 

  // Compute years having time-varying selectivities (PRELIMINARY_CALCS_SECTION)
  // ===========================================================================
  for (ifsh = 1; ifsh <= nfsh; ifsh++)
   for (iyr=1;iyr<=n_sel_ch_fsh(ifsh);iyr++)
    yrs_sel_ch_fsh(ifsh,iyr) = yrs_sel_ch_tmp(ifsh,iyr);
  for (isrv = 1; isrv <= nsrv; isrv++)
   for (iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
    yrs_sel_ch_srv(isrv,iyr) = yrs_sel_ch_tsrv(isrv,iyr);
  if (DebugOut > 1) cout << "yrs_sel_ch_fsh: " << yrs_sel_ch_fsh << endl;
  if (DebugOut > 1) cout << "yrs_sel_ch_srv: " << yrs_sel_ch_srv << endl;

  // set min & max sample size for stomach prey wts, lns
  for (rsp=1;rsp <= nspp; rsp++)
    for (rln=1;rln<=l_bins(rsp);rln++)
      for (iyrs=1;iyrs<=nyrs_stomwts(rsp);iyrs++)
       {
        if (stoms_w_N(rsp,rln,iyrs) <= min_SS_w(rsp))
          stoms_w_N(rsp,rln,iyrs) = 0;
        if (stoms_w_N(rsp,rln,iyrs) > max_SS_w(rsp))
          stoms_w_N(rsp,rln,iyrs) = max_SS_w(rsp);          
       }          
  for (rk_sp=1;rk_sp <= nspp_sq; rk_sp++)
    for (rln=1;rln<=r_lens(rk_sp);rln++)
      for (iyrs=1;iyrs<=nyrs_stomlns(rk_sp);iyrs++)
       {
        if (stoms_l_N(rk_sp,rln,iyrs) <= min_SS_l(rk_sp))
          stoms_l_N(rk_sp,rln,iyrs) = 0;
        if (stoms_l_N(rk_sp,rln,iyrs) > max_SS_l(rk_sp))
          stoms_l_N(rk_sp,rln,iyrs) = max_SS_l(rk_sp);
       }

  // Offset for diet (weights)
  offset_diet_w = 0;          
  for (rsp=1;rsp <= nspp; rsp++)
   for (iyrs=1; iyrs<= nyrs_stomwts(rsp); iyrs++)
    for (rln=1;rln<=l_bins(rsp);rln++)
     if (stoms_w_N(rsp,rln,iyrs) > 0)
      {
       for (ksp=1;ksp <=(nspp+1);ksp++)
        {
         rk_sp = (rsp-1)*(nspp+1)+ksp;
         if (diet_w_dat(rk_sp,rln,iyrs) > 0)
          offset_diet_w +=
               -1*stoms_w_N(rsp,rln,iyrs)*diet_w_dat(rk_sp,rln,iyrs) *
                                  log(diet_w_dat(rk_sp,rln,iyrs)+ 1.0e-10);
        }
      }
  
  // Offset for diet (lengths)
  offset_diet_l = 0;
  rk_sp = 0;
  for (rsp = 1; rsp <= nspp; rsp++)
   for (ksp = 1; ksp <= nspp; ksp++)
    {
     rk_sp = rk_sp + 1;
     for (rln = 1; rln <= l_bins(rsp); rln++)
      {
       Itot = int(sum(stoms_l_N(rk_sp,rln))); 
       if (Itot > 0)
        for (kln=1;kln<=l_bins(ksp);kln++)
         if (diet_l_dat(rk_sp,rln,kln) > 0)
          offset_diet_l += -1*Itot*diet_l_dat(rk_sp,rln,kln)*log(diet_l_dat(rk_sp,rln,kln)+ 1.0e-10);
      }
    }

  if (phase_SelFshCoff == -99) 
   for (ifsh=1;ifsh<=nfsh;ifsh++) phase_selcoff_fsh(ifsh) = -99;
  if (phase_SelSrvCoff == -99)
   for (isp=1;isp<=nspp;isp++) phase_selcoff_srv(isp) = -99;
 
  // Initial values for M, steepness, sigmar, R0, etc, (PRELIMINARY_CALCS_SECTION)
  // ==============================================================================
  if (Set_from_pin_file == 0) 
   {
    cout << "Crap" << endl;
    ipnt = 0;
    if (phase_M != -99) 
     for (isp=1;isp<=nspp;isp++)
      if (natmortphase2(isp) > 0) 
       {
        ipnt++;
        MEst(ipnt) = natmortprior(isp); 
       } 
    steepness = steepnessprior;
    log_sigmar = log_sigmarprior;
    if (phase_Rzero != -99) log_Rzero = R_guess;
    if (phase_LogRec != -99) mean_log_rec = 0;
    if (phase_fmort1 != -99) log_avg_fmort = -6.0;
    if (phase_RecDev != -99) rec_dev = 0;
    if (phase_SelFshCoff != -99) log_selcoffs_fsh = 0;
    if (phase_SelSrvCoff != -99) log_selcoffs_srv = 0;
    if (phase_fmort != -99) fmort_dev_est = 0;
    
    for (isrv=1;isrv<=nsrv;isrv++) log_q_srv(isrv) = log_qprior(isrv);
    if ( phase_SelSrvCoff != -99)
     {
      ipnt = 0;
      for (isrv=1;isrv<=nsrv;isrv++) 
       if (srv_sel_opt(isrv) == 2)
        {
         ipnt += 1;
         logsel_slope_srv_par(ipnt) = logsel_slp_in_srv(isrv);
         sel50_srv_par(ipnt) = sel_inf_in_srv(isrv);
        } 
     }   
  
    if (PhasePred3 != -99)
     for (rsp=1; rsp<=nspp; rsp++) Q_other_est = log(10000);
    if (PhasePred2 != -99)
     { 
      if (With_Pred == 0)
       logH_1 = -100;
      else
       logH_1 = -8.5;
      logH_1a = -3; 
      logH_1b = 0; 
      logH_2 = -9;
      logH_3 = -9;
      H_4 = 1;
     } 
    if (PhasePred1 != -99)
     {  
      rk_sp = 0;
      for (rsp = 1; rsp<=nspp; rsp++)
       for (ksp = 1; ksp <= nspp; ksp++)
        {
         rk_sp = rk_sp+1;
         log_gam_a(rsp) = 0.5;
         log_gam_b(rsp) = -0.5;
        }      
     }   
    
    if (phase_SelSrvCoff != -99)
     {
      ipnt = 0;
      for (isrv = 1; isrv <= nsrv; isrv++)
       if (srv_sel_opt(isrv) == 1)
       {
        isp =spp_srv(isrv);
        for (iage=1;iage<=nselages_srv(isrv,1);iage++)
         {
          ipnt += 1;
          log_selcoffs_srv(ipnt) = -log(1.0+mfexp(-log(19)*((double(iage)-8.0)/5.0) ));
         }
       }
     }    
   }
  cout << "END PRELIMINARY_CALCS_SECTION" << endl;
  if (DebugOut > 1) cout << "logH_1 = " << logH_1 << endl;

// ==========================================================================
// ==========================================================================

PROCEDURE_SECTION
 obj_fun.initialize();
 DoAll();


FUNCTION DoAll

  obj_fun.initialize();
 
  // Extract predation parameters from log space
  if (With_Pred == 0) 
   H_1 = 0;
  else 
   H_1 = mfexp(logH_1);
  H_1a = mfexp(logH_1a); 
  H_1b = mfexp(logH_1b); 
  H_2 = mfexp(logH_2);
  H_3 = mfexp(logH_3);
  gam_a = mfexp(log_gam_a);
  gam_b = mfexp(log_gam_b); 

  // Assign natural mortality (= prior if not estimated)
  ipnt = 0;
  for (isp=1;isp<=nspp;isp++)
   if (natmortphase2(isp) <= 0) 
    M(isp) = natmortprior(isp); 
   else
    {
     ipnt +=1;
     M(isp) = MEst(ipnt); 
    } 
  
  if (DebugOut > 1) cout << "logH_1 = " << logH_1 << endl;
  if (DebugOut > 1) cout << "H_1 = " << H_1 << endl;

  Get_Selectivity();
  if (With_Pred != 0) gamma_selectivity();
  Get_Mortality();
  
  Pmort_ua.initialize();
  eaten_ua.initialize();
  eaten_la.initialize();
  Get_Bzero();
  Get_Numbers_at_Age();
  Get_Survey_Predictions();
  Catch_at_Age();
  evaluate_the_objective_function();

  for (isp=1;isp<=nspp;isp++)
   for (iyr=styr_sp(isp);iyr<=endyr_all(isp);iyr++)
    SSBOut(isp,iyr) = Sp_Biom(isp,iyr);
    
  if (mceval_phase())
   {
    for (isp=1;isp<=nspp;isp++)
     for (iyr=styr_sp(isp);iyr<=endyr_all(isp);iyr++)
       McFile1 << Sp_Biom(isp,iyr) << " ";
    McFile1 << endl;  
    for (isp=1;isp<=nspp;isp++)
     for (iyr=styr_pred;iyr<=endyr;iyr++)
       McFile2 << natage(isp,iyr,1) << " ";
    McFile2 << endl;  
   }

  if (DebugOut > 4) for (isp=1; isp <= nspp; isp++)
   {
    cout << "natage(" << isp << ") = " << endl << natage(isp) << endl;
   }
 
 ObjTemp = obj_fun;

  // ============================
FUNCTION gamma_selectivity
  // ============================
  dvariable x_l_ratio;       // Log(mean(predLen@age)/mean(preyLen@age))
  dvariable LenOpt;          // Value of x_l_ratio where selectivity = 1
  dvariable gsum;
  int ncnt;
  int r_age,k_age;
  // -dhk June 26 2009

  if (DebugOut == 1) cout << "begin Gamma_selectivity" << endl;
  gam_ua.initialize();
  rk_sp = 0;
  for (rsp = 1; rsp <= nspp; rsp++)
   {
    LenOpt = 1.0e-10 + (gam_a(rsp)-1)*gam_b(rsp);
    for (ksp = 1; ksp <= nspp; ksp++)
     {
      rk_sp = rk_sp +1;
      for (r_age = 2; r_age <= nages(rsp); r_age++)
       {
        ncnt = 0; gsum = 1.0e-10;
        for (k_age = 1; k_age <= nages(ksp); k_age++)
         {
          // if prey are smaller than predator:
          if(mean_laa(rsp,r_age) > mean_laa(ksp,k_age))
           {
            x_l_ratio = log(mean_laa(rsp,r_age)/mean_laa(ksp,k_age));
            gam_ua(rk_sp,r_age,k_age) = 1.0e-10 +  (1.0e-10 +gam_a(rsp)-1) * log(x_l_ratio/LenOpt+1.0e-10) - 
                                    (1.0e-10 + x_l_ratio-LenOpt)/gam_b(rsp); // -dhk June 26 2009 
            ncnt += 1;
            gsum += mfexp(gam_ua(rk_sp,r_age,k_age));
           }  
          else
           gam_ua(rk_sp,r_age,k_age) = 0; 
         }
        for (k_age = 1; k_age <= nages(ksp); k_age++)
         if(mean_laa(rsp,r_age) > mean_laa(ksp,k_age))
          gam_ua(rk_sp,r_age,k_age) = 1.0e-10 + mfexp(gam_ua(rk_sp,r_age,k_age) - log(1.0e-10 + gsum/double(ncnt)));
                                       // -dhk June 26 2009
       } 
     }
   }  
  if (DebugOut > 5) cout << "gam_ua = " << endl << gam_ua << endl;
  if (DebugOut == 1) cout << "end gamma_selectivity" << endl;

  // ============================
FUNCTION Get_Selectivity
  // ============================
  int iage,max_sel_age;
  
  if (DebugOut == 1) cout << "begin Get_selectivity" << endl;
  
  ipnt = 0;
  for (ifsh = 1; ifsh <= nfsh; ifsh++)
   {
    isp = spp_fsh(ifsh);
    int fsh_sel_opt_tmp = fsh_sel_opt(ifsh);
    switch (fsh_sel_opt_tmp)
     {
      case 1 : // Fishery selectivity coefficients for POLLOCK, MACKEREL, COD 
       {
        if (phase_SelFshCoff > 0 || phase_selcoff_fsh(ifsh) == -99)
         {
          int isel_ch_tmp = 1 ;
          dvar_vector sel_coffs_tmp(1,nselages_fsh(ifsh,isel_ch_tmp));
          for (iyr=styr;iyr<=endyr;iyr++)
           {
            if (iyr==yrs_sel_ch_fsh(ifsh,isel_ch_tmp)) 
             {
              sel_coffs_tmp.initialize();
              for (iage=1;iage<=nselages_fsh(ifsh,isel_ch_tmp);iage++)
                { 
                  ipnt += 1; 
                  sel_coffs_tmp(iage) = log_selcoffs_fsh(ipnt); 
                } 
              avgsel_fsh(ifsh,isel_ch_tmp)  = log(mean(mfexp(sel_coffs_tmp)));
              if (isel_ch_tmp < n_sel_ch_fsh(ifsh)) isel_ch_tmp++;
             }
            max_sel_age = nselages_fsh(ifsh,isel_ch_tmp);
            log_sel_fsh(ifsh,iyr)(1,max_sel_age) = sel_coffs_tmp;
            log_sel_fsh(ifsh,iyr)(max_sel_age,nages(isp)) = log_sel_fsh(ifsh,iyr,max_sel_age);
            log_sel_fsh(ifsh,iyr) -= log(mean(mfexp(log_sel_fsh(ifsh,iyr) )));
           }
         }
       }
      break;
      case 2 : // Fishery asymptotic logistic NOT USED FOR POLLOCK, MACKEREL, COD
               // ===========================
        {
          cout << "case 2 Fishery asymptotic logistic not coded" << endl;
        }
      break;
     }
   }
 
  // Extract the selectivity parameters when there is logistic selectivity
  ipnt = 0;
  for (isrv = 1; isrv <= nsrv; isrv++)
   {
    if (srv_sel_opt(isrv) == 2)
     {
      ipnt += 1;
      logsel_slope_srv(isrv) = logsel_slope_srv_par(ipnt);
      sel50_srv(isrv) = sel50_srv_par(ipnt);
     }
    else
     {
      logsel_slope_srv(isrv) = 0;
      sel50_srv(isrv) = 0;
     } 
   }
  
  ipnt = 0;
  for (isrv = 1; isrv <= nsrv; isrv++)
   {
     isp =spp_srv(isrv);

     switch (srv_sel_opt(isrv))
      {
       case 1 : // Survey selectivity coefficients (mackerel)
                // ===============================
       if (phase_selcoff_srv(isrv) > 0 || phase_selcoff_srv(isrv) == -99)
        {
         int isel_ch_tmp = 1 ;
         dvar_vector sel_coffs_tmp(1,nselages_srv(isrv,isel_ch_tmp));
         for (iyr=styr;iyr<=endyr;iyr++)
           {
            if (iyr==yrs_sel_ch_srv(isrv,isel_ch_tmp)) 
             {
               sel_coffs_tmp.initialize();
               for (iage=1;iage<=nselages_srv(isrv,isel_ch_tmp);iage++)
                { 
                 ipnt += 1;
                 sel_coffs_tmp(iage) = log_selcoffs_srv(ipnt);
                } 
               avgsel_srv(isrv,isel_ch_tmp) = log(mean(mfexp(sel_coffs_tmp)));
               if (isel_ch_tmp < n_sel_ch_srv(isrv)) isel_ch_tmp++;
             }
           max_sel_age = nselages_srv(isrv,isel_ch_tmp);
           log_sel_srv(isrv,iyr)(1,max_sel_age) = sel_coffs_tmp;
           log_sel_srv(isrv,iyr)(max_sel_age,nages(isp)) = log_sel_srv(isrv,iyr,max_sel_age);
           log_sel_srv(isrv,iyr) -= log(mean(mfexp(log_sel_srv(isrv,iyr)(q_age_min(isrv),q_age_max(isrv))))); 
           log_sel_srv(isrv,iyr) -= log(mean(mfexp(log_sel_srv(isrv,iyr))));
          }
        }
       break;
       case 2 : // Survey asymptotic logistic (pollock, cod)
                // ==========================
        {
         //cout << "srv_sel_opt case 2" << endl;
         int isel_ch_tmp = 1 ; // selectivity change pointer can be incremented with n_sel_ch_srv and with srv
                              // in for loop to increment ipnt,isel_ch_tmp for multiple species
         for (iyr=styr;iyr<=endyr;iyr++) // this for loop not used, see comment below
           {
            if (iyr==yrs_sel_ch_srv(isrv,isel_ch_tmp)) 
             if (isel_ch_tmp<n_sel_ch_srv(isrv)) isel_ch_tmp++;
           } // option of incrementing isel_ch_tmp when different logistic survey selectivities 
             // occur at breaks in time (yrs_sel_ch_srv) was not used for cod, pollock, or mackerel

         sel_slope_srv(isrv) = mfexp(logsel_slope_srv(isrv));
         dvariable sel_slope_tmp = sel_slope_srv(isrv,isel_ch_tmp);
         dvariable sel50_tmp     = sel50_srv(isrv,isel_ch_tmp);
         for (iyr=styr;iyr<=endyr;iyr++)
           {
             if (iyr==yrs_sel_ch_srv(isrv,isel_ch_tmp)) // first year of survey only for cod, pollock and mackerel 
              {                                          // so isel_ch_tmp always = 1
               sel_slope_tmp = sel_slope_srv(isrv,isel_ch_tmp);
               sel50_tmp     =     sel50_srv(isrv,isel_ch_tmp);
               if (isel_ch_tmp<n_sel_ch_srv(isrv))      // n_sel_ch_srv always = 1
                 isel_ch_tmp++;                         // so never incremented
              }
             // fill in log_sel_srv values for all the selected age groups
             log_sel_srv(isrv,iyr)(1,nselages_srv(isrv,isel_ch_tmp)) = -1.*log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                                                  ( age_matrix(isp)(1,nselages_srv(isrv,isel_ch_tmp)) - sel50_tmp) ));
             // copy last selected age log_sel_srv value to the remaining older age groups
             log_sel_srv(isrv,iyr)(nselages_srv(isrv,isel_ch_tmp),nages(isp)) = log_sel_srv(isrv,iyr,nselages_srv(isrv,isel_ch_tmp));
           }
         }
        break;
    }        // end of srv_sel_opt switch
   }         // end of isrv loop
  sel_fsh = mfexp(log_sel_fsh);
  sel_srv = mfexp(log_sel_srv);
  if (DebugOut == 1) cout << "end Get_Selectivity" << endl;

  //=============================
FUNCTION Get_Mortality
  //=============================
  int age,ii;
  dvariable Temp,Temp1;

  if (DebugOut == 1) cout << "begin Get_Mortality" << endl;
  
  // Extract the rec_devs
  //int rec_adv.initialize();
  penal_rec_dev.initialize();
  int rec_adv = 0;
  for (isp = 1; isp <= nspp; isp++)
   for (iyr = styr_rec(isp); iyr <= endyr; iyr++)
    {
     rec_adv += 1;
     rec_dev_spp(isp,iyr) = rec_dev(rec_adv);
     Temp = (rec_dev(rec_adv)+6.5)/8.5;
     Temp1 = 10;
     for (ii=1;ii<=10;ii++)
      Temp1 *= Temp;
     penal_rec_dev += Temp1;
    }
  
  // Extract the qs
  for (isrv=1;isrv<=nsrv;isrv++) q_srv(isrv) = mfexp(log_q_srv(isrv) );
  
  // Extract the Fs (only for years WITH data)
  int F_adv = 0; //DHK initialize?
  for (ifsh=1;ifsh<=nfsh;ifsh++)
   for (iyr=styr;iyr<=endyr;iyr++)
    if (catch_bio(ifsh,iyr) > 10e-24) // 0 in NRM tpl -dhk apr 28 09
     {
      F_adv += 1;
      fmort_dev(ifsh,iyr) = fmort_dev_est(F_adv);
     }
    else
     fmort_dev(ifsh,iyr) = -100; // changed from -1000 -dhk Sep 1 09

  // Extract the other prey biomass
  for (rsp = 1; rsp<=nspp; rsp++)
   for (age=1;age<=nages(rsp);age++)
    Q_other_u(rsp,age) = mfexp(Q_other_est(rsp,age));
  
  surv    = mfexp(-1.0 * M);
  natmort = M;
  for (isp = 1; isp <= nspp; isp++) Z(isp) = M(isp);
  Fmort.initialize();
  ipnt = 0;
  for (ifsh = 1; ifsh <= nfsh; ifsh++)
   {
    isp = spp_fsh(ifsh);
    ipnt += 1;
    Fmort(isp) += mfexp(log_avg_fmort(ifsh) + fmort_dev(ifsh))+1.0e-10;
    for (iyr = styr; iyr <= endyr; iyr++)
     {
      F(ipnt,iyr) = mfexp(log_avg_fmort(ifsh) + fmort_dev(ifsh,iyr)) * sel_fsh(ifsh,iyr) + 1.0e-12;
      Z(isp,iyr) += F(ifsh,iyr);
     } 
   }
  if (DebugOut == 1) cout << "end get_mortality" << endl;   
  //=============================
FUNCTION Get_Bzero
  //=============================
  int ages;

  if (DebugOut == 1) cout << "Begin Get_Bzero" << endl;
  
  Rzero =  mfexp(log_Rzero);
  Sp_Biom.initialize();
  Bzero.initialize();
  AltStart();    // Iteratively calculate natage in styr_pred
  if (DebugOut == 1) cout << " Get_Bzero past AltStart" << endl;   

  // Extract all the recruitments
  for (isp=1;isp<=nspp;isp++)
   for (iyr=styr_pred;iyr<=endyr;iyr++)
    if (iyr < styr_rec(isp))
     natage(isp,iyr,1) = Rzero(isp);
    else 
     {
      natage(isp,iyr,1)  = Rzero(isp)*mfexp(rec_dev_spp(isp,iyr) + mean_log_rec(isp))+1.0e-10;//-dhk July 12 09
      mod_rec(isp,iyr) = natage(isp,iyr,1);     
     } 

  // Project forward to "correct" the age-structure
  for (iyr=styr_pred; iyr<styr; iyr++)
   {
    IyrPred = iyr;
    Compute_Predation();
    for (isp=1; isp<=nspp; isp++)
     {
      for (ages=2;ages<=nages(isp);ages++)
        natage(isp,iyr+1,ages) = natage(isp,iyr,ages-1)*S(isp,iyr,ages-1);
      natage(isp,iyr+1,nages(isp))+=natage(isp,iyr,nages(isp))*S(isp,iyr,nages(isp));
     }
   }

  if (DebugOut == 1) cout << " Get_Bzero past Compute_Predation" << endl;   

  // Equilibrium is located - find the parameters of the s-r relationship
  for (isp=1; isp<=nspp; isp++)
   {
    if (DebugOut > 1) cout << "inside Get_Bzero" << endl;
    if (DebugOut > 1) cout << "styr_rec = " << styr_rec << endl;
    if (DebugOut > 1) cout << "elem_prod(wt_pop(isp) , maturity(isp))" << endl;
    if (DebugOut > 1) cout << elem_prod(wt_pop(isp) , maturity(isp)) << endl;
    if (DebugOut > 1) cout << "natage(isp,styr_rec(isp)-1)" << endl;
    if (DebugOut > 1) cout << natage(isp,styr_rec(isp)-1) << endl;
    if (DebugOut > 1) cout << "pow(surv(isp),spmo_frac(isp))*natage(isp,styr_rec(isp)-1)" << endl;
    if (DebugOut > 1) cout << pow(surv(isp),spmo_frac(isp))*natage(isp,styr_rec(isp)-1) << endl;
    if (DebugOut > 1) cout << "pow(surv(isp),spmo_frac(isp))" << endl;
    if (DebugOut > 1) cout << pow(surv(isp),spmo_frac(isp)) << endl;
    
    iyr = styr_rec(isp)-1;
    Bzero(isp) = elem_prod(natage(isp,iyr),pow(S(isp,iyr),spmo_frac(isp))) * 
                 elem_prod(wt_pop(isp),maturity(isp));
    phizero(isp) = Bzero(isp)/Rzero(isp);
    switch (SrType)
     {
      case 1:
        alpha(isp) = log(-4.*steepness(isp)/(steepness(isp)-1.));
        break;
      case 2:
        alpha(isp)  =  Bzero(isp) * (1. - (steepness(isp) - 0.2) / (0.8*steepness(isp)) ) / Rzero(isp);
        beta(isp)   = (5. * steepness(isp) - 1.) / (4. * steepness(isp) * Rzero(isp));
        break;
      case 4:
        beta(isp)  = log(5.*steepness(isp))/(0.8*Bzero(isp)) ;
        alpha(isp) = log(Rzero(isp)/Bzero(isp))+beta(isp)*Bzero(isp);
        break;
     }
    Sp_Biom(isp)(styr_sp(isp),styr_rec(isp)) = Bzero(isp);
   }
  if (DebugOut == 1) cout << "end Get_Bzero" << endl; 

  //=============================
FUNCTION AltStart
  //=============================
  int isp,rsp,ksp,itno,age,ru,rln,ksp_type,rk_sp;
  dvar_matrix NN(1,nspp,1,nages);
  dvariable pred_effect,Term;
  dvariable ParA, ParB, ParC;
   if (DebugOut == 1) cout << " Begin AltStart" << endl;

  //Term.initialize(); // not initialized in NRM tpl -dhk apr 28 09
  
  // Initialize Z to M
  for (isp=1;isp<=nspp;isp++)
   Zlast(isp) = M(isp);
  
  // Calculate Z averaged over 5 iterations
   for (itno = 1; itno <= 25; itno++)
    {  
     //Set up virgin age-structure
    for (isp = 1; isp <= nspp; isp++)
     {
      NN(isp,1) = Rzero(isp);
      for (age=2; age<=nages(isp); age++)
        NN(isp,age) = NN(isp,age-1)*mfexp(-Zlast(isp,age-1));
      NN(isp,nages(isp)) /= (1.-mfexp(-Zlast(isp,nages(isp)))+1.0e-10); // - dhk 28 June 09
     } 
    
    // Mortality as a function of predator AGE (Eqn 3b)
    for (ksp=1;ksp<=nspp;ksp++)
     for (age=1;age<=nages(ksp);age++)
      {
       Zcurr(ksp,age) = M(ksp);
       if (With_Pred > 0)
        for (rsp=1;rsp<=nspp;rsp++)
         {
          ksp_type = (rsp-1)*(nspp+1)+ksp;
          rk_sp = (rsp-1)*nspp+ksp;
          for (ru=1;ru<=nages(rsp);ru++)
           {
            Term = H_1(ksp_type)*(1 + H_1a(rsp)*H_1b(rsp)/(double(ru)+H_1b(rsp)+1.0e-10));
            Zcurr(ksp,age) += Term*NN(rsp,ru)*gam_ua(rk_sp,ru,age) + 1.0e-10;   // - dhk 30 June 09
                                                                                // removed fabs() 1 July 09
           }   
         }
      }
      // added fabs() to calculation of Zcurr -dhk June 24 08

    // returned Zlast_pen back from here - dhk July 1 09

    // Average the Za
    for (isp=1;isp<=nspp;isp++)
     for(age=1;age<=nages(isp);age++)
       Zlast(isp,age) = sqrt(sqrt(Zcurr(isp,age)*Zlast(isp,age)))*sqrt(Zlast(isp,age)); // -dhk 8 July 09
      //Zlast(isp,age) = sqrt(Zcurr(isp,age)*Zlast(isp,age) + 1.0e-10); // -dhk 1 July 09
    }      
   
  // moved Zlast_pen from here -dhk 27 June 2009
  // Zlast penalty (moved back here 1 July 2009 -dhk)
    for (isp=1;isp<=nspp;isp++)
     {
      Zlast_pen(isp) = 0;
      for (age=1;age<=nages(isp);age++)
       Zlast_pen(isp) += 100*square(Zcurr(isp,age)-Zlast(isp,age) + 1.0e-10); // -dhk removed fabs() July 11 08
                                                                              // -dhk added 1.0e-10 July 9 09
     } // Zlast_pen location shouldn't matter now with itno removed -dhk July 2 09
 

  // Copy final NN into natage(,styr_pred,)
  for (isp=1;isp<=nspp;isp++)
   for (age=1;age<=nages(isp);age++)
     natage(isp,styr_pred,age) = NN(isp,age);

  // Calculate equilibrium N predators and prey in styr_pred for each species X age
 N_pred_eq.initialize();
 N_prey_eq.initialize();
 for (rsp=1;rsp<=nspp;rsp++)
  for (ru=1;ru<=nages(rsp);ru++)
    N_pred_yr(rsp,ru) = 1.0e-10; // -dhk 12 June 2009
 rk_sp=0;
  for (rsp=1;rsp<=nspp;rsp++)
   for (ksp=1; ksp<=nspp;ksp++)
    {
     rk_sp = rk_sp+1;
     for (ru=1;ru<=nages(rsp);ru++)
       for (age=1;age<=nages(ksp);age++)
         N_pred_eq(rsp,ru) += natage(rsp,styr_pred,ru) * gam_ua(rk_sp,ru,age);
     for (age=1;age<=nages(ksp);age++)
       for (ru=1;ru<=nages(rsp);ru++)
         N_prey_eq(ksp,age) += natage(ksp,styr_pred,age) * gam_ua(rk_sp,ru,age);
    }

  //=============================
FUNCTION Get_Numbers_at_Age
  //=============================
  int ages; 

  if (DebugOut == 1) cout << "Begin Get_Numbers_at_age" << endl;
 
  // Project ahead
  for (iyr=styr;iyr<endyr;iyr++)
   {
     // Compute the predation (update Z)
     IyrPred = iyr;
     Compute_Predation();
     for (isp=1; isp<=nspp; isp++)
      {
       for (ages=2;ages<=nages(isp);ages++)
        natage(isp,iyr+1,ages) = natage(isp,iyr,ages-1)*S(isp,iyr,ages-1);
       natage(isp,iyr+1,nages(isp))+=natage(isp,iyr,nages(isp))*S(isp,iyr,nages(isp));
      }
   } 
   
  // SSB
  IyrPred = endyr;
  Compute_Predation();
  for (isp=1; isp<=nspp; isp++)
   for (iyr=styr_rec(isp);iyr<=endyr;iyr++)
    Sp_Biom(isp,iyr)  = elem_prod(natage(isp,iyr),pow(S(isp,iyr),spmo_frac(isp))) * 
                        elem_prod(wt_pop(isp),maturity(isp));

  //=============================
FUNCTION Compute_Predation
  //=============================
 dvariable Pred_ratio;          // Predator ratio
 dvariable Prey_ratio;          // Prey ratio
 dvariable pred_effect;         // pred_resp * N(r,stage,y)
 dvariable NS_Z;                // N(k,y,a) * survival/Z;
 dvariable Tmort;               // Mortality on other
 dvariable Q_ksum_l;            // Diet sum
 dvariable Term;                // Linear adjustment for predation
 dvariable ParA, ParB, ParC;    // Parameters of H model
 int age,ksp_type, kall_type;   // Pointer
 
 // Only continue if predation is on; Pmort is zero otherwise
 if (With_Pred !=0) 
  { 
  
  //Term.initialize(); // not initialized in NRM tpl -dhk apr 28 09
  // Calculate N predators and prey in IyrPred for each species X age
 N_pred_yr.initialize();
 N_prey_yr.initialize();
 for (rsp=1;rsp<=nspp;rsp++)
  for (ru=1;ru<=nages(rsp);ru++)
    N_pred_yr(rsp,ru) = 1.0e-10; // -dhk 12 June 2009

 rk_sp=0;
  for (rsp=1;rsp<=nspp;rsp++)
   for (ksp=1; ksp<=nspp;ksp++)
    {
     rk_sp = rk_sp+1;
     for (ru=1;ru<=nages(rsp);ru++){
       for (age=1;age<=nages(ksp);age++)
         N_pred_yr(rsp,ru) += natage(rsp,IyrPred,ru) * gam_ua(rk_sp,ru,age);
       N_pred_yrs(rsp,IyrPred,ru) = N_pred_yr(rsp,ru);
       }
         
     for (age=1;age<=nages(ksp);age++){
       for (ru=1;ru<=nages(rsp);ru++)
         N_prey_yr(ksp,age) += natage(ksp,IyrPred,age) * gam_ua(rk_sp,ru,age);
       N_prey_yrs(ksp,IyrPred,age) = N_prey_yr(ksp,age);
      }         
    }

  // Calculate predator functional response
  rk_sp = 0;
  rksp = 0;

  for (rsp=1;rsp<=nspp;rsp++) 
   for (ksp=1;ksp<=(nspp+1);ksp++)
    {
     rk_sp += 1;
     if (ksp <= nspp)
       rksp += 1;
     for (r_age=1;r_age<=nages(rsp);r_age++)
      for (k_age=1; k_age<=kk_ages(ksp); k_age++)
      {
       Term = 1.0e-10 + H_1(rk_sp)*(1 + H_1a(rsp)*H_1b(rsp)/(double(r_age)+H_1b(rsp)+1.0e-10));

     N_pred_eqs(rsp,IyrPred,r_age) = N_pred_eq(rsp,r_age);

       if (ksp <= nspp)
        {
     N_prey_eqs(ksp,IyrPred,k_age) = N_prey_eq(ksp,k_age);

         // Predator-prey ratios
         Pred_ratio = (N_pred_yr(rsp,r_age)+1.0e-10)/(N_pred_eq(rsp,r_age)+1.0e-10);
         Prey_ratio = (N_prey_yr(ksp,k_age)+1.0e-10)/(N_prey_eq(ksp,k_age)+1.0e-10);
         Pred_r(rsp,IyrPred,r_age) = Pred_ratio;
         Prey_r(ksp,IyrPred,k_age) = Prey_ratio;

         

         if (resp_type == 1 || current_phase() < PhasePred2-1)      // Holling Type I (linear)
          pred_resp(rk_sp,IyrPred,r_age,k_age) = 1.0e-10 + Term;
         else 
          if (resp_type == 2) // Holling Type II
           {
            pred_resp(rk_sp,IyrPred,r_age,k_age) = 1.0e-10 + Term*(1+H_2(rksp)+1.0e-10) /
              ( 1 + H_2(rksp) * Prey_ratio + 1.0e-10);
           }
         else 
          if (resp_type == 3) // Holling Type III
           {
            pred_resp(rk_sp,IyrPred,r_age,k_age) = 1.0e-10 +
              Term*(1+H_2(rksp))*pow((Prey_ratio + 1.0e-10),H_4(rksp)) /
              (1 + H_2(rksp) * pow((Prey_ratio + 1.0e-10),H_4(rksp)) + 1.0e-10 );
           }
         else 
          if (resp_type == 4) // predator interference
           {
            pred_resp(rk_sp,IyrPred,r_age,k_age) = 1.0e-10 + Term*(1+H_2(rksp)+1.0e-10) /
              ( 1 + H_2(rksp)*Prey_ratio + H_3(rksp)*(Pred_ratio-1) + 1.0e-10);  
           }
         else 
          if (resp_type == 5) // predator preemption
           {
            pred_resp(rk_sp,IyrPred,r_age,k_age) = 1.0e-10 + Term*(1+H_2(rksp)+1.0e-10) /
              ( (1+H_2(rksp)*Prey_ratio) * (1+H_3(rksp)*(Pred_ratio-1))+1.0e-10);  
           }
         else 
          if (resp_type == 6) // Hassell-Varley
           {
            pred_resp(rk_sp,IyrPred,r_age,k_age) = 1.0e-10 + Term*(2+H_2(rksp)+ 1.0e-10) /
              (1.0+ H_2(rksp)*Prey_ratio + pow((Prey_ratio+1.0e-10),H_4(rksp)) + 1.0e-10 );  
           }
         else 
          if (resp_type == 7) // Ecosim
           {
            pred_resp(rk_sp,IyrPred,r_age,k_age) = 1.0e-10 + Term /
              (1 + H_3(rksp)*(Pred_ratio - 1 + 1.0e-10));  
           }
        }
       else  // "other" is linear
        pred_resp(rk_sp,IyrPred,r_age,1) = 1.0e-10 + Term; 
      }            // end of r_ages, k_ages loop
                   // =========================
   }


   // Mortality as a function of predator AGE (Eqn 3b)
   for (rsp = 1;rsp <= nspp; rsp++)
    for (ksp = 1; ksp <= nspp; ksp++)
     {
      ksp_type = (rsp-1)*(nspp+1)+ksp;
      rk_sp = (rsp-1)*nspp+ksp;
      for (ru=1;ru<=nages(rsp);ru++)
       for (age=1;age<=nages(ksp);age++)
        {
         pred_effect = pred_resp(ksp_type,IyrPred,ru,age)*gam_ua(rk_sp,ru,age);
         Vmort_ua(rk_sp,ru,IyrPred,age) = pred_effect*natage(rsp,IyrPred,ru);
        }
     }

   // Accumulate the total mortality (Equations 2a / 2b)
   rk_sp = 0; 
   for (rsp=1;rsp<=nspp;rsp++)
    for (ksp=1;ksp<=nspp;ksp++)
     {
      rk_sp += 1;
      for (ru=1;ru<=nages(rsp);ru++)
       Pmort_ua(rk_sp,IyrPred) += Vmort_ua(rk_sp,ru,IyrPred);
      Z(ksp,IyrPred) += Pmort_ua(rk_sp,IyrPred);
     }
  }    
 
  
  // Convert from total mortality to survival
  for (isp=1; isp<=nspp; isp++)
   S(isp,IyrPred) = 1.0e-10+mfexp(-1*Z(isp,IyrPred));//1.0e-30 in NRM tpl -dhk apr 28 09

 // Only continue if predation is on; Pmort is zero otherwise
 if (With_Pred !=0) 
 if (IyrPred >= styr)
  {

    // Numbers eaten (of modeled prey species); Equations 7 and 8
   for (ksp = 1; ksp <= nspp; ksp++)
    for (age=1;age<=nages(ksp);age++)
    {
     // Relative number 
     NS_Z = natage(ksp,IyrPred,age)*(1-mfexp(-Z(ksp,IyrPred,age)))/Z(ksp,IyrPred,age);      
     for (rsp = 1; rsp <= nspp; rsp++)
      {
       rk_sp = (rsp-1)*nspp+ksp;
     
       // Numbers eaten by predator age and length (Eqn 7a & 7b)
       for (ru = 1;ru <= nages(rsp); ru++)
        {
         eaten_ua(rk_sp, ru,IyrPred,age) = Vmort_ua(rk_sp,ru,IyrPred,age)*NS_Z;
         for (rln = 1; rln <= l_bins(rsp); rln++)
          eaten_la(rk_sp,rln,IyrPred,age) += eaten_ua(rk_sp,ru,IyrPred,age)*al_key(rsp,ru,rln);
        } 
      }
    } 
     
    // Mass eaten (including "other")
    for (rsp = 1; rsp <= nspp; rsp++)
     for (ksp = 1; ksp <= nspp+1; ksp++)
      {
       // Pointers to locations of data  
       ksp_type = (rsp-1)*(nspp+1)+ksp;
       kall_type = (rsp-1)*nspp+ksp;
       if (ksp <= nspp)
        {
         // Results by length (Eqn 8a)
         for (rln = 1; rln <= l_bins(rsp); rln++)
          {
           Q_mass_l(ksp_type,IyrPred,rln) = 0;
           for (age=1;age<=nages(ksp);age++)
            Q_mass_l(ksp_type,IyrPred,rln) += eaten_la(kall_type,rln,IyrPred,age)*wt_pop(ksp,age);
          }  
         // Results by age (Eqn 8b)
         for (ru = 1;ru <= nages(rsp); ru++)
          {
           Q_mass_u(ksp_type,IyrPred, ru) = 0;
           for (age=1;age<=nages(ksp);age++)
            Q_mass_u(ksp_type,IyrPred, ru) += eaten_ua(kall_type, ru,IyrPred,age)*wt_pop(ksp,age);
          } 
        }
       else
        {
         for (ru = 1; ru <= nages(rsp); ru++)
          {
           pred_effect = pred_resp(ksp_type,IyrPred,ru,1);
           Tmort = pred_effect * natage(rsp,IyrPred,ru); // Eq.3b ======
           Q_mass_u(ksp_type,IyrPred,ru)  = Q_other_u(rsp,ru)*(1.0-mfexp(-Tmort));
          } 
         for (rln=1;rln<=l_bins(rsp);rln++)
          {
           Q_mass_l(ksp_type,IyrPred,rln) = 0;
           for (ru=1;ru<=nages(rsp);ru++)
             Q_mass_l(ksp_type,IyrPred,rln) += Q_mass_u(ksp_type,IyrPred,ru)*al_key(rsp,ru,rln);
          }  
        } 
      }
     
    // Total up the consumption by each predator and normalize (Eqn 15)
    for (rsp = 1; rsp <= nspp; rsp++)
     for (rln=1;rln<=l_bins(rsp);rln++)
      {
       rk_sp = (rsp-1)*(nspp+1);
       Q_ksum_l.initialize();
       for (ksp = 1; ksp <= (nspp+1); ksp++)
        Q_ksum_l += Q_mass_l(rk_sp+ksp,IyrPred,rln) + 1.0e-10; //1.e-20 in NRM tpl -dhk apr 28 09
       for (ksp = 1; ksp <= (nspp+1); ksp++)
       {
        Q_hat(rk_sp+ksp,IyrPred,rln) = (1.0e-10+Q_mass_l(rk_sp+ksp,IyrPred,rln)/Q_ksum_l); // changed paranthesis -dhk apr 28 09
       }
      } 
  }     
  
  //=============================
FUNCTION Get_Survey_Predictions
  //=============================

  if (DebugOut == 1) cout << "Begin Get_Survey_Predictions" << endl;
  
  dvariable sum_tmp;
  sum_tmp.initialize();
  int yy;
  
  for (isrv=1; isrv<=nsrv; isrv++)
   {
    isp = spp_srv(isrv);
    dvar_matrix natage_spp = natage(isp);
    for (iyr=styr; iyr<=endyr; iyr++)
    { 
     pred_srv(isrv,iyr) = 1.0e-10 + q_srv(isrv) * natage(isp,iyr) * elem_prod(sel_srv(isrv,iyr) , wt_srv(isrv,iyr));//-dhk apr 28 09
    }
    for (iyr=1; iyr<=nyrs_srv_comp(isrv); iyr++)
     {
      yy = yrs_srv_comp(isrv,iyr); 
      dvar_vector tmp_n =elem_prod(sel_srv(isrv,yy),natage(isp,yy));  
      sum_tmp = sum(tmp_n);
      eac_srv(isrv,iyr) = tmp_n/sum_tmp;
     }
   } 

  //=============================
FUNCTION Catch_at_Age
  //=============================

  if (DebugOut == 1) cout << "Begin Catch_at_Age" << endl;
  
  for (ifsh=1; ifsh<=nfsh; ifsh++)
   {
    isp = spp_fsh(ifsh);
    for (iyr=styr;iyr<=endyr;iyr++)
     catage(ifsh,iyr) = elem_prod(elem_div(F(ifsh,iyr),Z(isp,iyr)),elem_prod(1.-S(isp,iyr),natage(isp,iyr)));
    dvar_matrix Ctmp = catage(ifsh); // Copy 3darray to matrix for efficiency...
    for (iyr=styr; iyr<=endyr; iyr++)
     pred_catch(ifsh,iyr) = Ctmp(iyr)*wt_fsh(ifsh,iyr);
    for (iyr=1; iyr<=nyrs_fsh_comp(ifsh); iyr++)
     eac_fsh(ifsh,iyr)=Ctmp(yrs_fsh_comp(ifsh,iyr))/(1.0e-10+sum(Ctmp(yrs_fsh_comp(ifsh,iyr))));
   }

  // ============================
FUNCTION evaluate_the_objective_function
  dvariable Temp_obj;
  
  // ============================
  count_iters = count_iters + 1;
  
  Cat_Like();
  Rec_Like();
  Age_Like();
  Srv_Like();
  Sel_Like();
  Fmort_Pen();
  Compute_priors();
  diet_like1.initialize();
  diet_like2.initialize();
  ration_like.initialize();
  ration_pen.initialize();
  if (With_Pred > 0)
   {
    ration();
    ration_Like();
    diet_wt_Like();
    diet_len_Like();
   }
  if (active(log_Rzero))
   for(isp=1; isp<=nspp; isp++)
    obj_fun += .5 * square(mean_log_rec(isp)); // penalty to constrain Rzero when active
  Temp_obj = obj_fun;  

  obj_comps.initialize();
  obj_comps(1) = sum(catch_like);
  if (Disc_any_phases != 0 & current_phase() < 2-Initial_phase+1 & Disc_first_phase > 0) obj_comps(1) *= 0.1;
  obj_comps(2) = sum(age_like_fsh);
  obj_comps(3) = sum(sel_like_fsh);
  if (Disc_any_phases != 0 & current_phase() == 1-Initial_phase+1) obj_comps(3) = 0;
  obj_comps(4) = sum(surv_like);
  if (Disc_any_phases != 0 & current_phase() < 2-Initial_phase+1 & Disc_first_phase > 0) obj_comps(4) *= 0.1;
  obj_comps(5) = sum(age_like_srv);
  obj_comps(6) = sum(sel_like_srv);
  if (Disc_any_phases != 0 & current_phase() == 1-Initial_phase+1) obj_comps(6) = 0;
  obj_comps(7) = sum(rec_like);
  obj_comps(8) = sum(fpen);
  obj_comps(9) = sum(post_priors_srvq);
  obj_comps(10)= sum(post_priors);
  obj_comps(11) = penal_rec_dev;
  obj_comps(12) = sum(Zlast_pen);
  obj_comps(13) = sum(ration_like);
  if (Disc_any_phases != 0 & current_phase() < 2-Initial_phase+1 & Disc_first_phase > 0) obj_comps(13) *= 0.1;
  obj_comps(14) = diet_like1;
  obj_comps(15) = diet_like2;
  obj_comps(16) = sum(ration_pen);
  obj_fun     += sum(obj_comps);
  obj_fun     += dummy2*dummy2;
  if (current_phase() == Terminal_phase)
    cout << "ration_pen = " << ration_pen << endl;
  
  if (DebugOut > 2)
   {
    cout << "obj_comps(1) = catch_like = " << catch_like << endl;
    cout << "obj_comps(2) = age_like_fsh = " << age_like_fsh << endl;
    cout << "obj_comps(3) = sel_like_fsh = " << endl << sel_like_fsh << endl;
    cout << "obj_comps(4) = surv_like = " << surv_like << endl;
    cout << "obj_comps(5) = age_like_srv = " << age_like_srv << endl;
    cout << "obj_comps(6) = sel_like_srv = " << endl << sel_like_srv << endl;
    cout << "obj_comps(7) = rec_like = " << endl << rec_like << endl;
    cout << "obj_comps(8) = fpen = " << endl << fpen << endl;
    cout << "obj_comps(9) = post_priors_srvq = " << post_priors_srvq << endl;
    cout << "obj_comps(10) = post_priors = " << post_priors << endl;
    cout << "obj_comps(12) = Zlast_pen = " << Zlast_pen << endl;
    cout << "obj_comps(13) = ration_like = " << ration_like << endl;
    cout << "obj_comps(14) = diet_like1 = " << diet_like1 << endl;
    cout << "obj_comps(15) = diet_like2 = " << diet_like2 << endl;
    cout << "obj_comps(16) = ration_pen = " << ration_pen << endl;
    cout << endl << "end of evaluate_the_objective_function, iteration " << count_iters << endl;
    cout << "  ======================================================================" << endl;
    cout << "  ======================================================================" << endl << endl;
   } 
  cout << "obj_fun: " << obj_fun << " Iteration: " << count_iters << " Phase: " << current_phase() << endl;
  cout << obj_comps << " " << Temp_obj  << endl;
  
    

  // ============================
FUNCTION Cat_Like
  // ============================
  if (DebugOut == 1) cout << "begin Cat_Like" << endl;
  
  catch_like.initialize();
  for (ifsh=1; ifsh<=nfsh; ifsh++)
   {
    isp = spp_fsh(ifsh);
    catch_like(ifsh) = catchbiomass_pen(isp) * norm2(log(catch_bio(ifsh) +.000001)-log(pred_catch(ifsh) +.000001));
   }

  // ============================
FUNCTION Rec_Like
  // ============================
  if (DebugOut == 1) cout << "Begin Rec_like" << endl;
  
  rec_like.initialize();
  if (active(rec_dev) || phase_RecDev == -99)
  for (isp=1; isp<=nspp; isp++)
   {
    sigmar(isp)     =  mfexp(log_sigmar(isp));
    sigmarsq(isp)   =  square(sigmar(isp));
    dvariable SSQRec;
    SSQRec.initialize();
    if (Disc_any_phases == 0 || current_phase() > 2-Initial_phase+1)
      {
        pred_rec(isp) = SRecruit(Sp_Biom(isp)(styr_rec(isp),endyr).shift(styr_rec(isp))(styr_rec(isp),endyr));
        dvar_vector chi = log(elem_div(mod_rec(isp)(styr_rec_est(isp),endyr_rec_est(isp)) ,
                                  pred_rec(isp)(styr_rec_est(isp),endyr_rec_est(isp))));
        SSQRec   = norm2( chi ) ;
        m_sigmar(isp)   =  sqrt( SSQRec  / nrecs_est(isp));
        m_sigmarsq(isp) =  m_sigmar(isp) * m_sigmar(isp)   ;
        rec_like(isp,1) += norm2(chi + sigmarsq(isp)/2.)/(2*sigmarsq(isp)) + nrecs_est(isp)*log_sigmar(isp);
        // rec_like(isp,1) above changed to match Dorn(2002) -dhk Jul 3 2008. old (wrong) form used in NRM tpl -dhk apr 28 09
      }
     if (Disc_any_phases == 0 || current_phase() >= phase_RecDev)
      {
        // Variance term for the parts not estimated by sr curve
        rec_like(isp,4) += .5*norm2( rec_dev_spp(isp)(styr_rec(isp),styr_rec_est(isp)) )/sigmarsq(isp) + (styr_rec_est(isp)-styr_rec(isp))*log(sigmar(isp)) ; 
        if (endyr>endyr_rec_est(isp))
          rec_like(isp,4) += .5*norm2( rec_dev_spp(isp)(endyr_rec_est(isp),endyr  ) )/sigmarsq(isp) + (endyr-endyr_rec_est(isp))*log(sigmar(isp)) ; 
      }
     else
      {
        rec_like(isp,2) += norm2( rec_dev_spp(isp)( styr_rec_est(isp),endyr_rec_est(isp)) ) ;
      }
     rec_like(isp,2) += norm2( rec_dev_spp(isp)( styr_rec_est(isp),endyr) ) ;
   }

  // ============================
FUNCTION Age_Like
  // ============================
  if (DebugOut == 1) cout << "Begin Age_Like " << endl;
  
  age_like_fsh.initialize();
  for (ifsh=1;ifsh<=nfsh;ifsh++)
   {
    isp = spp_fsh(ifsh);
    if (comp_type(isp) == 1)
     ec_fsh(ifsh) = eac_fsh(ifsh);
    else if (comp_type(isp)==2)
     ec_fsh(ifsh) = eac_fsh(ifsh)*al_key(isp);
    for ( iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
     {
      ec_fsh(ifsh,iyr) /= sum(ec_fsh(ifsh,iyr));
      age_like_fsh(ifsh) -= nsmpl_fsh(ifsh,iyr)*(oc_fsh(ifsh,iyr) + 0.001) * log(ec_fsh(ifsh,iyr) + 0.001 ) ;
     }
   }
  age_like_fsh -= offset_fsh;
  
  age_like_srv.initialize();
  for (isrv=1;isrv<=nsrv;isrv++)
   {
    isp = spp_srv(isrv);
    if (comp_type(isp) == 1)
     ec_srv(isrv) = eac_srv(isrv);
    else if (comp_type(isp) == 2)
     ec_srv(isrv) = eac_srv(isrv)*al_key(isp);
    for (iyr=1;iyr<=nyrs_srv_comp(isrv);iyr++)
     {
      ec_srv(isrv,iyr) /= sum(ec_srv(isrv,iyr));
      age_like_srv(isrv) -= nsmpl_srv(isrv,iyr)*(oc_srv(isrv,iyr) + 0.001) * log(ec_srv(isrv,iyr) + 0.001 ) ;
     } 
   }
  age_like_srv-=offset_srv;

  // ============================
FUNCTION dvar_vector SRecruit(const dvar_vector& Stmp)
  // ============================
  RETURN_ARRAYS_INCREMENT();
  
  int i_sp = isp;
  dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
  switch (SrType)
    {
    case 1:
      RecTmp = elem_prod((Stmp / phizero(i_sp)) , mfexp( alpha(i_sp) * ( 1. - Stmp / Bzero(i_sp) ))) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = elem_prod(Stmp , 1. / ( alpha(i_sp) + beta(i_sp) * Stmp));        //Beverton-Holt form
      break;
    case 3:
      RecTmp = Rzero(isp)*mfexp(mean_log_rec(i_sp));                    //Avg recruitment
      break;
    case 4:
      RecTmp = elem_prod(Stmp , mfexp( alpha(i_sp)  - Stmp * beta(i_sp))) ; //Old Ricker form
      break;
    }
  
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

  // ============================
FUNCTION Srv_Like
  // ============================
  // Fit to indices (Normal) 
  dvariable qest,nest;
  
  surv_like.initialize();
  for (isrv=1;isrv<=nsrv;isrv++)
   {
    if (Disc_any_phases != 0 & current_phase() < 2-Initial_phase+1)
     {
      qest.initialize(); nest.initialize();
      for (iyr=1;iyr<=nyrs_srv(isrv);iyr++)
       { nest += 1; qest += obs_srv(isrv,iyr)/pred_srv(isrv,yrs_srv(isrv,iyr)); }
      qest = qest/nest;
      surv_like(isrv) += square(qest-1)*100000;
     }  
    else
     qest = 1;
   
    for (iyr=1;iyr<=nyrs_srv(isrv);iyr++)
     surv_like(isrv) += square(obs_srv(isrv,iyr) - qest*pred_srv(isrv,yrs_srv(isrv,iyr)) ) / 
                                   (2.*obs_se_srv(isrv,iyr)*obs_se_srv(isrv,iyr));
   } 

  // ============================
FUNCTION Sel_Like 
  // ============================
  if (DebugOut == 1) cout << "Begin Sel_Like" << endl;
  
  sel_like_fsh.initialize();
  sel_like_srv.initialize();
  for (ifsh=1;ifsh<=nfsh;ifsh++)  //FISHERIES
   {                              //=========
     isp = spp_fsh(ifsh);
     if (active(log_selcoffs_fsh) || phase_SelFshCoff == -99)
      {
       for (iyr=1;iyr<=n_sel_ch_fsh(ifsh);iyr++)
        {
         int i_iyr = yrs_sel_ch_fsh(ifsh,iyr) ;
         sel_like_fsh(ifsh,1) += curv_pen_fsh(ifsh)*norm2(first_difference(
                                 first_difference(log_sel_fsh(ifsh,i_iyr ))));
         // This part is the penalty on the change itself--------------
         if (iyr>1)
           {
            dvariable var_tmp = square(sel_change_in_fsh(ifsh,i_iyr ));
            sel_like_fsh(ifsh,2) += .5*norm2( log_sel_fsh(ifsh,i_iyr-1) - log_sel_fsh(ifsh,i_iyr) ) / var_tmp ;
           }
         int nagestmp = nselages_fsh(ifsh,1);
         for (int j=seldecage(isp);j<=nagestmp;j++)
           {
            dvariable difftmp = log_sel_fsh(ifsh,i_iyr,j-1)-log_sel_fsh(ifsh,i_iyr,j) ;
            if (difftmp > 0.)
              sel_like_fsh(ifsh,3)    += .5*square( difftmp ) / seldec_pen_fsh(ifsh);
           }        
         obj_fun += 20 * square(avgsel_fsh(ifsh,iyr)); // To normalize selectivities
        }
      }
    }
  for (isrv=1;isrv<=nsrv;isrv++)  //SURVEYS
   {                              //=======
    isp = spp_fsh(isrv);
    if (phase_selcoff_srv(isrv) > 0 || phase_SelSrvCoff == -99)
     {
       for (iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
        {
          int i_iyr = yrs_sel_ch_srv(isrv,iyr) ;
          sel_like_srv(isrv,1) += curv_pen_srv(isrv)*norm2(first_difference(
                                                 first_difference(log_sel_srv(isrv,i_iyr))));
          // This part is the penalty on the change itself--------------
          if (iyr>1)
            {
             dvariable var_tmp = square(sel_change_in_srv(isrv,i_iyr ));
             sel_like_srv(isrv,2)    += .5*norm2( log_sel_srv(isrv,i_iyr-1) - log_sel_srv(isrv,i_iyr) ) 
                                   / var_tmp ;
            }
          int nagestmp = nselages_srv(isrv,1);
          for (int j=seldecage(isp);j<=nagestmp;j++)
            {
             dvariable difftmp = log_sel_srv(isrv,i_iyr,j-1)-log_sel_srv(isrv,i_iyr,j) ;
             if (difftmp > 0.)
               sel_like_srv(isrv,3)    += .5*square( difftmp ) / seldec_pen_srv(isrv);
            }
          obj_fun += 20. * square(avgsel_srv(isrv,iyr));  // To normalize selectivities
        }
     }
   }

  // ============================
FUNCTION ration
  // ============================
  int jyr,age;
  dvariable n_avg,numer,denom;

  if (DebugOut == 1) cout << "Begin Ration" << endl;

  // Equations 9 and 10
  omega_hat_ave.initialize();
  omega_hat.initialize();
  for (rsp=1;rsp<=nspp;rsp++)
   for (age=1;age<=nages(rsp);age++)
    {
     // Calculate year-specific values
     numer.initialize();
     denom.initialize();
     //n_avg.initialize(); //-dhk June 30 08
     for (iyr=styr;iyr<=endyr;iyr++)
      {
       // Average abundance
       //n_avg =1.0e-10 + natage(rsp,iyr,age) * mfexp(-0.5*Z(rsp,iyr,age));
       n_avg = 1.0e-10 + natage(rsp,iyr,age) * sqrt(S(isp,iyr,age));  // added 1.0e-10 -dhk June 24 08. not in NRM tpl -dhk apr 28 09
       // find total consumption by this age-class
       rk_sp = (rsp-1)*(nspp+1);
       for (ksp=1;ksp<=(nspp+1);ksp++)
        omega_hat(rsp,iyr,age) += Q_mass_u(rk_sp+ksp,iyr,age);
     
       numer += omega_hat(rsp,iyr,age)/365+1.0e-10;
       denom += n_avg;
       
       // normalize
       omega_hat(rsp,iyr,age) /= (n_avg*365);
       
      }
     omega_hat_ave(rsp,age) = numer/denom;
    }
     
  // ============================
FUNCTION ration_Like
  // ============================
  int ages;
  //ration_like.initialize(); // not initialized in NRM tpl -dhk apr 28 09

  if (DebugOut == 1) cout << "Begin Ration_Like" << endl;
  // Likelihood (Eqn 11)
  for (rsp=1;rsp<=nspp; rsp++)
   {
    ration_like(rsp) = 0; // NRM tpl form -dhk apr 28 09
    for (ages=2;ages<=nages(rsp);ages++) // don't include age zero in likelihood
     ration_like(rsp) += 0.5 * square(log(omega_hat_ave(rsp,ages)+1.0e-10) - 
                          log(omega_vB(rsp,ages)))/square(sd_ration(rsp));
   }  
   for(isp=1; isp<=nspp; isp++)
     //ration_pen(isp) = 0;
     {
      for (iage=1;iage<=nages(isp);iage++)
      {
       mean_ohat = 0;
       for(iyr=styr;iyr<=endyr;iyr++)  mean_ohat += omega_hat(isp,iyr,iage);
       mean_ohat /= float(endyr-styr+1);
       for(iyr=styr;iyr<=endyr;iyr++)
         ration_pen(isp) += 20 * square(omega_hat(isp,iyr,iage)-mean_ohat);
      }
     }
                     
  // ============================
FUNCTION diet_wt_Like
  // ============================
  int iyrs;

  //loop_count = 0;
  // Likelihood (Eqn 14)
  diet_like1.initialize();
  for (rsp=1;rsp <= nspp; rsp++)
   for (iyrs=1; iyrs<= nyrs_stomwts(rsp); iyrs++)
    {
     iyr = yrs_stomwts(rsp,iyrs);
     for (rln=1; rln<=l_bins(rsp); rln++)
      if (stoms_w_N(rsp,rln,iyrs) > 0)
       for (ksp=1; ksp <=(nspp+1); ksp++)
        {
         rk_sp = (rsp-1)*(nspp+1)+ksp;
         if (diet_w_dat(rk_sp,rln,iyrs) > 0)
          {
          diet_like1 +=
               -1*stoms_w_N(rsp,rln,iyrs)*diet_w_dat(rk_sp,rln,iyrs) *
                                  log(Q_hat(rk_sp,iyr,rln)+ 1.0e-10);
         }    
        }
     }

  diet_like1 -= offset_diet_w;

  // ============================
FUNCTION diet_len_Like
  // ============================
  dvariable Denom,TotN;
  Denom.initialize(); // -dhk june 30 08. not initialized in NRM tpl -dhk apr 28 09
  TotN.initialize();  // -dhk june 30 08. not initialized in NRM tpl -dhk apr 28 09
  int Itot;
  int loop_count;
  loop_count = 0;

  if (DebugOut == 1) cout << "Begin Diet_len_Like" << endl;

  // Calculate the predicted fraction by length-class (Eqn 17)
  rk_sp=0;
  T_hat.initialize();
  diet_like2.initialize();
  for (rsp = 1; rsp <= nspp; rsp++)
   for (ksp = 1; ksp <= nspp; ksp++)
    {
     dvar_vector eaten_lmy(1,l_bins(ksp)); // no. of prey@length eaten by a predator length during iyr
     rk_sp = rk_sp + 1;
     for (rln = 1; rln <= l_bins(rsp); rln++)
      {
       TotN = sum(stoms_l_N(rk_sp,rln));
       Itot = int(sum(stoms_l_N(rk_sp,rln))); 
       if (Itot > 0)
        {
         // This is Equation 17
         for (stm_yr = 1; stm_yr <= nyrs_stomlns(rk_sp); stm_yr++)
          if (stoms_l_N(rk_sp,rln,stm_yr) > 0)
           {
            iyr = yrs_stomlns(rk_sp,stm_yr);
            eaten_lmy = eaten_la(rk_sp,rln,iyr) * al_key(ksp);
            T_hat(rk_sp,rln) += stoms_l_N(rk_sp,rln,stm_yr)* eaten_lmy;
           }  
         // Renormalize the eaten vector 
         Denom = sum(T_hat(rk_sp,rln))+1.0e-10;
         T_hat(rk_sp,rln) /= Denom;
         // This is equation 16
         for (kln=1;kln<=l_bins(ksp);kln++)
          if (diet_l_dat(rk_sp,rln,kln) > 0)
          {
           //if (T_hat(rk_sp,rln,kln) < 1.0e-10) T_hat(rk_sp,rln,kln) = 1.0e-10; // -dhk Jul 3 08. not in NRM tpl -dhk apr 28 09
           diet_like2 += -TotN*diet_l_dat(rk_sp,rln,kln)*log(T_hat(rk_sp,rln,kln)+1.0e-10);
          }
        }
      }
    }
  diet_like2 -= offset_diet_l;

  // ============================
FUNCTION Fmort_Pen
  // ============================
  dvariable totalN, totalF, TotalG, Temp;
  
  if (DebugOut == 1) cout << "Begin Fmort_Pen" << endl;
  
  fpen.initialize();
  for(isp=1; isp<=nspp; isp++)
   {      
    if (Disc_any_phases != 0 & current_phase()<3-Initial_phase+1) // penalize High F's for beginning phases
     fpen(isp,1) += 10.* norm2(Fmort(isp) - .2);
    else 
     fpen(isp,1) +=.001*norm2(Fmort(isp) - .2);
   }
  for (ifsh = 1; ifsh <= nfsh; ifsh++)
   {
    isp = spp_fsh(ifsh);
    totalF.initialize(); totalN.initialize(); TotalG.initialize();
    for (iyr=styr;iyr<=endyr;iyr++)
     if (catch_bio(ifsh,iyr) >1.0e-24) // changed from 10e-24 to 1.0e-24 where NRM tpl has 0 -dhk apr 28 09
      { totalF += fmort_dev(ifsh,iyr); totalN += 1; 
        Temp = mfexp(log_avg_fmort(ifsh) + fmort_dev(ifsh,iyr));
        TotalG += 100/(1+mfexp(-log(19)*(Temp-2)/0.25)); 
      }
    fpen(isp,2) += 20*TotalG;    
    fpen(isp,2) += 20.*square(totalF/totalN);
   } 

  // ============================
FUNCTION Compute_priors
  // ============================
  post_priors.initialize();
  post_priors_srvq.initialize();
  for (isrv=1;isrv<=nsrv;isrv++)
   if (active(log_q_srv(isrv)))
    post_priors_srvq(isrv) += square(q_srv(isrv)-qprior(isrv))/(2*cvqprior(isrv)*cvqprior(isrv)); 
  for(isp=1; isp<=1; isp++)
   {
    if (active(MEst) || phase_M == -99)
     post_priors(1) += square(M(isp)-natmortprior(isp))/(2*cvnatmortprior(isp)*cvnatmortprior(isp)); 

    if (active(steepness))
     post_priors(2) += square(steepness(isp)-steepnessprior(isp))/(2*cvsteepnessprior(isp)*cvsteepnessprior(isp)); 

    if (active(log_sigmar))
     post_priors(3) += square(sigmar(isp)-sigmarprior(isp))/(2*cvsigmarprior(isp)*cvsigmarprior(isp)); 
   } 


REPORT_SECTION
  // ============================
  // this report file outputs source code for R
  
 // Andre_crap();

  report << " #### INDEX VALUES ##### " << endl;
  report << "nspp <- " << nspp << endl;
  report << "nfsh <- " << nfsh << endl;
  report << "nsrv <- " << nsrv << endl;
  report << "styr <- " << styr << endl;
  report << "endyr <- " << endyr << endl;
  report << "nyrs <- " << nyrs << endl;
  report << "styr_pred <- " << styr_pred << endl;
  report << "nyrs_pred <- " << nyrs_pred << endl;
  report << "styr_rec <- c(";
   for (isp=1;isp<=nspp;isp++)
    {
     report << styr_rec(isp);
     if(isp<nspp)
      {
       report << ", ";
      }
     else
      {
       report << ")" << endl;
      }
    }
  report << "oldest_age <- c(";
   for (isp=1;isp<=nspp;isp++)
    {
     report << oldest_age(isp);
     if(isp<nspp)
      {
       report << ", ";
      }
     else
      {
       report << ")" << endl;
      }
    }
  report << "l_bins <- c("; 
   for (isp=1;isp<=nspp;isp++)
    {
     report << l_bins(isp);
     if(isp<nspp)
      {
       report << ", ";
      }
     else
      {
       report << ")" << endl;
      }
    }
  report << "nages <- c("; 
   for (isp=1;isp<=nspp;isp++)
    {
     report<< nages(isp);
     if(isp<nspp)
      {
       report << ", ";
      }
     else
      {
       report << ")" << endl;
      }
    }
  report << "nages_fsh <- c("; 
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
     report<< nages_fsh(ifsh);
     if(ifsh < nfsh)
      {
       report << ", ";
      }
     else
      {
       report << ")" << endl;
      }
    }
  report << "nfsh_spp <- c("; 
   for (isp=1;isp<=nspp;isp++)
    {
     report << nfsh_spp(isp);
     if(isp<nspp)
      {
      report << ", ";
      }
     else
      {
       report  << ")" << endl;
      }
    }
  report << "comp_type <- c("; 
   for (isp=1;isp<=nspp;isp++)
    {
     report << comp_type(isp);
     if(isp<nspp)
      {
       report << ", ";
      }
     else if(isp==nspp)
      {
       report << ")" << endl;
      }
    }
  report << "nyrs_fsh_comp <- c(";
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
     report << nyrs_fsh_comp(ifsh);
     if(ifsh<nfsh)
      {
        report << ", ";
      }
     else
      {
       report  << ")" << endl;
      }
    }
  report << "spp_fsh <- c(";
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
     report << spp_fsh(ifsh);
     if(ifsh<nfsh)
      {
        report << ", ";
      }
     else
      {
       report  << ")" << endl;
      }
    }
  report << "spp_srv <- c(";
   for (isrv=1;isrv<=nsrv;isrv++)
    {
     report << spp_srv(isrv);
     if(isrv<nsrv)
      {
        report << ", ";
      }
     else
      {
       report  << ")" << endl;
      }
    }

  report << "ncomps_fsh <- c(";
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
     report << ncomps_fsh(ifsh);
     if(ifsh<nfsh)
      {
        report << ", ";
      }
     else
      {
       report  << ")" << endl;
      }
    }

 // yrs_srv(1,nsrv,1,nyrs_srv)
  report << "yrs_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "yrs_srv[[" << isrv << "]] <- c(" << yrs_srv(isrv,1);
      for(iyr=2;iyr<=nyrs_srv(isrv);iyr++)
       {
        report << ", " << yrs_srv(isrv,iyr);
       }
          report   << ")" << endl;
    }
  report << endl;
 //obs_srv
  report << "obs_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "obs_srv[[" << isrv << "]] <- c(" << obs_srv(isrv,1);
      for(iyr=2;iyr<=nyrs_srv(isrv);iyr++)
       {
        report << ", " << obs_srv(isrv,iyr);
       }
          report   << ")" << endl;
    }
  report << endl;
 //obs_se_srv
  report << "obs_se_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "obs_se_srv[[" << isrv << "]] <- c(" << obs_se_srv(isrv,1);
      for(iyr=2;iyr<=nyrs_srv(isrv);iyr++)
       {
        report << ", " << obs_se_srv(isrv,iyr);
       }
          report  << ")" << endl;
    }
  report << endl;


  report << "ncomps_srv <- c(";
   for (isrv=1;isrv<=nsrv;isrv++)
    {
     report << ncomps_srv(isrv);
     if(isrv<nsrv)
      {
        report << ", ";
      }
     else
      {
       report  << ")" << endl;
      }
    }
  report << "nyrs_srv_comp <- c(";
   for (isrv=1;isrv<=nsrv;isrv++)
    {
     report << nyrs_srv_comp(isrv);
     if(isrv<nsrv)
      {
        report << ", ";
      }
     else
      {
       report  << ")" << endl;
      }
    }

  report << endl << "###### DATA ##########" << endl;

//############ al_key ################
//##
  report << endl << "al_key<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "al_key[[" << isp << "]] <- c(" << endl << al_key(isp,1,1);
      for (iage=1;iage<=nages(isp);iage++)  // ages per species
        {
          if (iage == 1)
            istart = 2;
          else
            istart = 1;

          for (icmp=istart;icmp<=l_bins(isp);icmp++)  // lengths per species
            {
              report << ", " << al_key(isp,iage,icmp);        
            }
        }
          report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      isp = spp_fsh(isp);
      report << "dim(al_key[[" << isp << "]])<-c(" << l_bins(isp);
      report << ", " << nages(isp) << ")" << endl;
      report << "al_key[[" << isp << "]]<-t(al_key[[" << isp << "]])" << endl;
    }
//#### oc_fsh matrices ################
//##
  report << endl << "oc_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << "oc_fsh[[" << ifsh << "]] <- c(" << endl << oc_fsh(ifsh,1,1);
      isp = spp_fsh(ifsh);
      for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (icmp=istart;icmp<=ncomps_fsh(isp);icmp++)
            {
              report << ", " << oc_fsh(ifsh,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
    }
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      isp = spp_fsh(ifsh);
      report << "dim(oc_fsh[[" << ifsh << "]])<-c(" << ncomps_fsh(isp);
      report << ", " << nyrs_fsh_comp(ifsh) << ")" << endl;
      report << "oc_fsh[[" << ifsh << "]]<-t(oc_fsh[[" << ifsh << "]])" << endl;
    }
  report << endl;

//#### wt_pop matrices ################
//##
  report << endl << "wt_pop<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "wt_pop[[" << isp << "]] <- c(" << endl <<wt_pop(isp,1);
      for (iage=2;iage<=nages(isp);iage++)
        {
              report << ", " << wt_pop(isp,iage);        
        }
          report  << endl << ")" << endl;
    }
  report << endl;

  report << " #### ESTIMATED BIOLOGICAL PARAMETERS #### " << endl;
  report << "M<-c(";
   for (isp=1;isp<=nspp;isp++)
    {
      report << M(isp);
      if (isp<nspp)   report << ", ";
      else if (isp == nspp) report << ")" << endl;
    }
  report << "mean_log_rec<- c(";
   for (isp=1;isp<=nspp;isp++)
    {
      report << mean_log_rec(isp);
      if (isp<nspp)   report << ", ";
      else if (isp == nspp) report << ")" << endl;
    }  
  report << "steepness<- c(";
   for (isp=1;isp<=nspp;isp++)
    {
      report << steepness(isp);
      if (isp<nspp)   report << ", ";
      else if (isp == nspp) report << ")" << endl;
    }  
  report << "log_Rzero<- c(";
   for (isp=1;isp<=nspp;isp++)
    {
      report << log_Rzero(isp);
      if (isp<nspp)   report << ", ";
      else if (isp == nspp) report << ")" << endl;
    }  
  report << "Bzero<- c(";
   for (isp=1;isp<=nspp;isp++)
    {
      report << Bzero(isp);
      if (isp<nspp)   report << ", ";
      else if (isp == nspp) report << ")" << endl;
    }  
  report << "rec_dev_spp <- list()" << endl;
   for (isp=1;isp<=nspp;isp++)
    {
      report << "rec_dev_spp[[" << isp << "]] <- c(";
      for (iyr=styr_rec(isp);iyr<=endyr;iyr++)
      {
       report << rec_dev_spp(isp,iyr);
       if (iyr < endyr)  report << ", ";
       else if(iyr==endyr) report << ")" << endl; 
      } 
    }
  report << "log_sigmar<- c(";
   for (isp=1;isp<=nspp;isp++)
    {
      report << log_sigmar(isp);
      if (isp<nspp)   report << ", ";
      else if (isp == nspp) report << ")" << endl;
    }  
  report << endl;
  report << " #### ESTIMATED FISHERIES PARAMETERS #### " << endl;
  report << "log_selcoffs_fsh <- list()" << endl;
   ipnt = 0;
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << "log_selcoffs_fsh[[" << ifsh << "]] <- c(";
     for (int iyr=1;iyr<=n_sel_ch_fsh(ifsh);iyr++)
      for (int iage=1;iage<=nselages_fsh(ifsh,iyr);iage++)
      {
        ipnt += 1;
        i_lage =  nselages_fsh(ifsh,iyr);
        report << log_selcoffs_fsh(ipnt);
        if (iage < i_lage)
          report << ", ";
        else if ((iage == i_lage) && (iyr < n_sel_ch_fsh(ifsh)))
          report << ", ";
      }
       report << ")" << endl;
    }
  report << "sel_slope_fsh <- list()" << endl;
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << "sel_slope_fsh[[" << ifsh << "]] <- c(";
      for (int iyr=1;iyr<=n_sel_ch_fsh(ifsh);iyr++)
       {
        report << sel_slope_fsh(ifsh,iyr);
        if (iyr < n_sel_ch_fsh(ifsh))
          report << ", ";
        else if (iyr == n_sel_ch_fsh(ifsh)) report << ")" << endl;
       }
    } 
  report << "log_avg_fmort <- c(";  
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << log_avg_fmort(ifsh);
      if (ifsh<nfsh) report << ", ";
      else if (ifsh == nfsh) report << ")" << endl;
    }
  report << "fmort_dev <- list()" << endl;
   for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << "fmort_dev[[" << ifsh << "]]<-c(";
      for (iyr=styr;iyr<=endyr;iyr++)
        {
          report << fmort_dev(ifsh,iyr);
          if (iyr < endyr) report << ", ";
          else if (iyr == endyr) report << ")" << endl; 
        }
    }    
  report << endl;

  report << " #### ESTIMATED SURVEY PARAMETERS #### " << endl;
  report << "log_q_srv<- c(";
   for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << log_q_srv(isrv);
      if (isrv<nsrv)   report << ", ";
      else if (isrv == nsrv) report << ")" << endl;
    }   
  report << "log_selcoffs_srv <- list()" << endl;
   ipnt = 0;
   for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "log_selcoffs_srv[[" << isrv << "]] <- c(";
     for (int iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
      if (srv_sel_opt(isrv) == 1)
       {
        for (int iage=1;iage<=nselages_srv(isrv,iyr);iage++)
        {
         ipnt += 1;
         i_lage =  nselages_srv(isrv,iyr);
         report << log_selcoffs_srv(ipnt);
         if (iage < i_lage)  report << ", ";
        }
        report << ")" << endl;
       }
      else
       {
        for (int iage=1;iage<=nselages_srv(isrv,iyr);iage++)
        {
         i_lage =  nselages_srv(isrv,iyr);
         report << 0;
         if (iage < i_lage)  report << ", ";
        }
        report << ")" << endl;
       }
        
    }
  report << "logsel_slope_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
    {
     report << "logsel_slope_srv[[" << isrv << "]] <- c(";
     for (int iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
      {
       report << logsel_slope_srv(isrv,iyr);
       if (iyr < n_sel_ch_srv(isrv))
         report << ", ";
       else if (iyr == n_sel_ch_srv(isrv)) report << ")" << endl;
      }
    }   
  report << "sel50_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "sel50_srv[[" << isrv << "]] <- c(";
      for (int iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
       {
        report << sel50_srv(isrv,iyr);
        if (iyr < n_sel_ch_srv(isrv))
          report << ", ";
        else if (iyr == n_sel_ch_srv(isrv)) report << ")" << endl;
       }
    } 
  report << "sel_slope_srv <- list()" << endl;
   for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "sel_slope_srv[[" << isrv << "]] <- c(";
      for (int iyr=1;iyr<=n_sel_ch_srv(isrv);iyr++)
       {
        report << sel_slope_srv(isrv,iyr);
        if (iyr < n_sel_ch_srv(isrv))
          report << ", ";
        else if (iyr == n_sel_ch_srv(isrv)) report << ")" << endl;
       }
    }

  report << endl;

  report << " #### DERIVED PARAMETERS AND DATA MATRICES #### " << endl;
//############ eac_fsh matrices ################
//##
  report << "eac_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << "eac_fsh[[" << ifsh << "]] <- c(" << endl << eac_fsh(ifsh,1,1);
      isp = spp_fsh(ifsh);
      for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (icmp=istart;icmp<=nages(isp);icmp++)
            {
              report << ", " << eac_fsh(ifsh,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
    }
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      isp = spp_fsh(ifsh);
      report << "dim(eac_fsh[[" << ifsh << "]])<-c(" << nages(isp);
      report << ", " << nyrs_fsh_comp(ifsh) << ")" << endl;
      report << "eac_fsh[[" << ifsh << "]]<-t(eac_fsh[[" << ifsh << "]])" << endl;
    }

//############ ec_fsh matrices ################
//##
  report << endl << "ec_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      report << "ec_fsh[[" << ifsh << "]] <- c(" << endl << ec_fsh(ifsh,1,1);
      isp = spp_fsh(ifsh);
      for (iyr=1;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (icmp=istart;icmp<=ncomps_fsh(isp);icmp++)
            {
              report << ", " << ec_fsh(ifsh,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
    }
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      isp = spp_fsh(ifsh);
      report << "dim(ec_fsh[[" << ifsh << "]])<-c(" << ncomps_fsh(isp);
      report << ", " << nyrs_fsh_comp(ifsh) << ")" << endl;
      report << "ec_fsh[[" << ifsh << "]]<-t(ec_fsh[[" << ifsh << "]])" << endl;
    }

//############## log_sel_fsh matrices #####
//##
  report << endl << "log_sel_fsh<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
         {
           report << "log_sel_fsh[[" << ifsh << "]]<- c(" << endl << log_sel_fsh(ifsh,styr,1);
           isp = spp_fsh(ifsh);
           for (iyr=styr;iyr<=endyr;iyr++)
             {
               if (iyr == styr)
                 istart = 2;
               else
                 istart = 1;

               for (icmp=istart;icmp<=nages(isp);icmp++)
                {
                  report << ", " << log_sel_fsh(ifsh,iyr,icmp);        
                }
              }
          report  << endl << ")" << endl;
         }

  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      isp = spp_fsh(ifsh);
      report << "dim(log_sel_fsh[[" << ifsh << "]])<-c(" << nages(isp);
      report << ", " << nyrs << ")" << endl << endl;
      report << "log_sel_fsh[[" << ifsh << "]]<-t(log_sel_fsh[[" << ifsh << "]])" << endl;
    }

//############## log_sel_srv matrices #####
//##
  report << endl << "log_sel_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
         { 
           report << "log_sel_srv[[" << isrv << "]]<- c(" << endl << log_sel_srv(isrv,styr,1);
           isp = spp_srv(isrv);
           for (iyr=styr;iyr<=endyr;iyr++)
             {
               if (iyr == styr)
                 istart = 2;
               else
                 istart = 1;

               for (icmp=istart;icmp<=nages(isp);icmp++)
                 {
                   report << ", " << log_sel_srv(isrv,iyr,icmp);        
                 }
              }
          report  << endl << ")" << endl;
         }
  for (isrv=1;isrv<=nsrv;isrv++)
    {
      isp = spp_srv(isrv);
      report << "dim(log_sel_srv[[" << isrv << "]])<-c(" << nages(isp);
      report << ", " << nyrs << ")" << endl;
      report << "log_sel_srv[[" << isrv << "]]<-t(log_sel_srv[[" << isrv << "]])" << endl;
    }

//############## pred_srv #####
//## matrix pred_srv(1,nsrv,styr,endyr)

  report << endl << "pred_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
         { 
           report << "pred_srv[[" << isrv << "]]<- c(" << endl << pred_srv(isrv,styr);
           isp = spp_srv(isrv);
           for (iyr=styr+1;iyr<=endyr;iyr++)
             {
               report << ", " << pred_srv(isrv,iyr);
             }
           report << ")" << endl;
          }        

//############## Spawning biomass
//   matrix Sp_Biom(1,nspp,styr_sp,endyr_all)  
  report << endl << "Sp_Biom<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
     report << "Sp_Biom[[" << isp << "]]<-c(" << endl << Sp_Biom(isp,styr);
     for (iyr=styr_sp(isp)+1; iyr<=endyr_all(isp); iyr++)
      {
        report << ", " << Sp_Biom(isp,iyr);
      }
      report << ")" << endl;
    }

//############## predicted recruits
// matrix pred_rec(1,nspp,styr_rec,endyr_all)
  report << endl << "pred_rec<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    iyr = styr_rec(isp);
    report << "pred_rec[[" << isp << "]]<-c(" << endl << pred_rec(isp,iyr);
    for (iyr=styr_rec(isp)+1; iyr<=endyr_all(isp); iyr++)
     {
       report << ", " << pred_rec(isp,iyr);
     }
     report << ")" << endl;
    }

//############## numbers at age #####
//##
  report << "natage <- list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "natage[[" << isp << "]] <- c(" << endl << natage(isp,styr,1);
      for (iyr=styr;iyr<=endyr;iyr++)
        {
          if (iyr == styr)
            istart = 2;
          else
            istart = 1;

          for (iage=istart;iage<=nages(isp);iage++)
            {
              report <<  ", " << natage(isp,iyr,iage);
            }
         }
      report  << endl << ")" << endl;      
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "dim(natage[[" << isp << "]])<-c(" << nages(isp) << ", " << nyrs << ")" << endl;
      report << "natage[[" << isp << "]]<-t(natage[[" << isp << "]])" << endl;
    }
//############## total mortality Z #####
//##
  report << "Z <- list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "Z[[" << isp << "]] <- c(" << endl << Z(isp,styr,1);
      for (iyr=styr;iyr<=endyr;iyr++)
        {
          if (iyr == styr)
            istart = 2;
          else
            istart = 1;

          for (iage=istart;iage<=nages(isp);iage++)
            {
              report <<  ", " << Z(isp,iyr,iage);
            }
         }
      report  << endl << ")" << endl;      
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "dim(Z[[" << isp << "]])<-c(" << nages(isp) << ", " << nyrs << ")" << endl;
      report << "Z[[" << isp << "]]<-t(Z[[" << isp << "]])" << endl;
    }

//############## fishing mortality F #####
//##
  report << "F <- list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      if (ifsh < 3)
        isp = ifsh;
      else
        isp = 3;
      report << "F[[" << ifsh << "]] <- c(" << endl << F(ifsh,styr,1);
      for (iyr=styr;iyr<=endyr;iyr++)
        {
          if (iyr == styr)
            istart = 2;
          else
            istart = 1;

          for (iage=istart;iage<=nages(isp);iage++)
            {
              report <<  ", " << F(ifsh,iyr,iage);
            }
         }
      report  << endl << ")" << endl;      
    }
  for (ifsh=1;ifsh<=nfsh;ifsh++)
    {
      if (ifsh < 3)
        isp = ifsh;
      else
        isp = 3;
      report << "dim(F[[" << ifsh << "]])<-c(" << nages(isp) << ", " << nyrs << ")" << endl;
      report << "F[[" << ifsh << "]]<-t(F[[" << ifsh << "]])" << endl;
    }

//############## mortality rate by predator AGE Pmort_ua #####
//## 3darray  Pmort_ua(1,nspp_sq,1,nyrs_pred,1,k_ages)

  report << "Pmort_ua <- list()" << endl;
  for (isp=1; isp<=nspp_sq; isp++)
   {
    report << "Pmort_ua[[" << isp << "]] <- c(" << endl << Pmort_ua(isp,styr_pred,1);
    for (iyr=styr_pred; iyr<=endyr; iyr++)
     {
          if (iyr == styr_pred)
            istart = 2;
          else
            istart = 1;
      for (iage=istart; iage<=k_ages(isp); iage++)
       {
        report <<  ", " << Pmort_ua(isp,iyr,iage);
       }
     }
      report  << endl << ")" << endl;      
   }
  for (isp=1;isp<=nspp_sq;isp++)
    {
      report << "dim(Pmort_ua[[" << isp << "]])<-c(" << k_ages(isp) << ", " << nyrs_pred << ")" << endl;
      report << "Pmort_ua[[" << isp << "]]<-t(Pmort_ua[[" << isp << "]])" << endl;
    }

//############## pred_catch #####
//##
  report << endl << "pred_catch<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
         { 
           report << "pred_catch[[" << ifsh << "]]<- c(" << endl << pred_catch(ifsh,styr);
           for (iyr=styr+1;iyr<=endyr;iyr++)
             {
               report << ", " << pred_catch(ifsh,iyr);
             }
           report << ")" << endl;
          }  
//############## catch_bio #####
//##
  report << endl << "catch_bio<-list()" << endl;
  for (ifsh=1;ifsh<=nfsh;ifsh++)
         { 
           report << "catch_bio[[" << ifsh << "]]<- c(" << endl << catch_bio(ifsh,styr);
           for (iyr=styr+1;iyr<=endyr;iyr++)
             {
               report << ", " << catch_bio(ifsh,iyr);
             }
           report << ")" << endl;
          }        


  //3darray  omega_hat(1,nspp,1,nyrs_pred,1,nages)

  //matrix  omega_hat_ave(1,nspp,1,nages)

//############ omega_hat ################
//##   3darray  omega_hat(1,nspp,1,nyrs_pred,1,nages)
  report << endl << "omega_hat<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    report << "omega_hat[[" << isp << "]] <- c(" << endl << omega_hat(isp,styr_pred,1);
     for (iyr=styr_pred; iyr <= endyr; iyr++)
     {
       for (iage=1;iage<=nages(isp);iage++)  // ages per species
         {
               report << ", " << omega_hat(isp,iyr,iage); 
         }
     }
           report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "omega_hat[[" << isp << "]] <- omega_hat[[" << isp << "]]";
      report << "[2:length(omega_hat[[ " << isp << "]])]" << endl;
      report << "dim(omega_hat[[" << isp << "]])<-c(" << nages(isp);
      report << ", " << nyrs_pred << ")" << endl;
      report << "omega_hat[[" << isp << "]]<-t(omega_hat[[" << isp << "]])" << endl;
    }

//############ omega_vB ################
//##   init_matrix omega_vB(1,nspp,1,nages);
  report << endl << "omega_vB<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "omega_vB[[" << isp << "]] <- c(" << endl << omega_vB(isp,1);
      for (iage=2;iage <= nages(isp);iage++)
        {
              report << ", " << omega_vB(isp,iage);        
        }
          report  << endl << ")" << endl;
    }
  report << endl;

//############ nsmpl_fsh ###################
//   init_matrix nsmpl_fsh(1,nfsh,1,nyrs_fsh_comp)       // Effective sample sizes
  report << endl << "nsmpl_fsh <- list()" << endl;
  for (ifsh =1; ifsh <= nfsh; ifsh++)
    {
      report << "nsmpl_fsh[[" << ifsh << "]] <- c(" << endl <<nsmpl_fsh(ifsh,1);
      for (iyr=2;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
              report << ", " << nsmpl_fsh(ifsh,iyr);        
        }
          report  << endl << ")" << endl;
    }

//############ yrs_fsh_comp ################
//   init_imatrix yrs_fsh_comp(1,nfsh,1,nyrs_fsh_comp)     // Years with age data
  report << endl << "yrs_fsh_comp <- list()" << endl;
  for (ifsh =1; ifsh <= nfsh; ifsh++)
    {
      report << "yrs_fsh_comp[[" << ifsh << "]] <- c(" << endl <<yrs_fsh_comp(ifsh,1);
      for (iyr=2;iyr<=nyrs_fsh_comp(ifsh);iyr++)
        {
              report << ", " << yrs_fsh_comp(ifsh,iyr);        
        }
          report  << endl << ")" << endl;
    }
  report << endl;

//############ eac_srv matrices ################
//##
  report << endl << "eac_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "eac_srv[[" << isrv << "]] <- c(" << endl << eac_srv(isrv,1,1);
      isp = spp_srv(isrv);
      for (iyr=1;iyr<=nyrs_srv_comp(isrv);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (icmp=istart;icmp<=nages(isp);icmp++)
            {
              report << ", " << eac_srv(isrv,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
    }
  for (isrv=1;isrv<=nsrv;isrv++)
    {
      isp = spp_srv(isrv);
      report << "dim(eac_srv[[" << isrv << "]])<-c(" << nages(isp);
      report << ", " << nyrs_srv_comp(isrv) << ")" << endl;
      report << "eac_srv[[" << isrv << "]]<-t(eac_srv[[" << isrv << "]])" << endl;
    }

//############ ec_srv matrices ################
//##
  report << endl << "ec_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "ec_srv[[" << isrv << "]] <- c(" << endl <<ec_srv(isrv,1,1);
      isp = spp_srv(isrv);
      for (iyr=1;iyr<=nyrs_srv_comp(isrv);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (icmp=istart;icmp<=ncomps_srv (isp);icmp++)
            {
              report << ", " << ec_srv(isrv,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
    }
  for (isrv=1;isrv<=nsrv;isrv++)
    {
      isp = spp_srv(isrv);
      report << "dim(ec_srv[[" << isrv << "]])<-c(" <<   ncomps_srv(isp);
      report << ", " << nyrs_srv_comp(isrv) << ")" << endl;
      report << "ec_srv[[" << isrv << "]]<-t(ec_srv[[" << isrv << "]])" << endl;
    }

//############ nsmpl_srv ###################
//   init_matrix nsmpl_srv(1,nsrv,1,nyrs_srv_comp)       // Effective sample sizes
  report << endl << "nsmpl_srv <- list()" << endl;
  for (isrv =1; isrv <= nsrv; isrv++)
    {
      report << "nsmpl_srv[[" << isrv << "]] <- c(" << endl <<nsmpl_srv(isrv,1);
      for (iyr=2;iyr<=nyrs_srv_comp(isrv);iyr++)
        {
              report << ", " << nsmpl_srv(isrv,iyr);        
        }
          report  << endl << ")" << endl;
    }

//############ yrs_srv_comp ################
//   init_imatrix yrs_srv_comp(1,nsrv,1,nyrs_srv_comp)     // Years with age data
  report << endl << "yrs_srv_comp <- list()" << endl;
  for (isrv =1; isrv <= nsrv; isrv++)
    {
      report << "yrs_srv_comp[[" << isrv << "]] <- c(" << endl <<yrs_srv_comp(isrv,1);
      for (iyr=2;iyr<=nyrs_srv_comp(isrv);iyr++)
        {
              report << ", " << yrs_srv_comp(isrv,iyr);        
        }
          report  << endl << ")" << endl;
    }
  report << endl;

//############ oc_srv matrices ################
//##
  report << endl << "oc_srv<-list()" << endl;
  for (isrv=1;isrv<=nsrv;isrv++)
    {
      report << "oc_srv[[" << isrv << "]] <- c(" << endl <<oc_srv(isrv,1,1);
      isp = spp_srv(isrv);
      for (iyr=1;iyr<=nyrs_srv_comp(isrv);iyr++)
        {
          if (iyr == 1)
            istart = 2;
          else
            istart = 1;

          for (icmp=istart;icmp<=ncomps_srv (isp);icmp++)
            {
              report << ", " << oc_srv(isrv,iyr,icmp);        
            }
        }
          report  << endl << ")" << endl;
    }
  for (isrv=1;isrv<=nsrv;isrv++)
    {
      isp = spp_srv(isrv);
      report << "dim(oc_srv[[" << isrv << "]])<-c(" <<   ncomps_srv(isp);
      report << ", " << nyrs_srv_comp(isrv) << ")" << endl;
      report << "oc_srv[[" << isrv << "]]<-t(oc_srv[[" << isrv << "]])" << endl;
    }
  report << endl;

//#### gamma selectivity of predator age u on prey age a
  //  3darray  gam_ua(1,nspp_sq,1,r_ages,1,k_ages)             // gamma selectivity of predator age u on prey age a
  report << "gam_ua<-list()" << endl;
  for (isp=1; isp<=nspp_sq; isp++)
   {
    report << "gam_ua[[" << isp << "]]<-c(";
    for (ru=1; ru<=r_ages(isp); ru++)
     for (k_age=1; k_age<=k_ages(isp); k_age++)
      {
        report << gam_ua(isp,ru,k_age);
        if (ru < r_ages(isp)) report << ", ";
        else if (k_age < k_ages(isp)) report << ", ";
        else report << ")" << endl;
      }
   }
  for (isp =1; isp <= nspp_sq; isp++)
   {
    report << "dim(gam_ua[[" << isp << "]]) <- c(";
    report << k_ages(isp) << ", " << r_ages(isp) << ")" << endl;
    report << "gam_ua[[" << isp << "]] <- t(gam_ua[[" << isp << "]])" << endl;
   }


//#### mass consumed by given predator age u (Q_mass_u) ################
//##
  // 3darray  Q_mass_u(1,nspp_sq2,1,nyrs_pred,1,rr_ages)
  // 3darray  Q_mass_u(1,nspp_sq2,1,nyrs_pred,1,rr_ages)
  report << "Q_mass_u <- list()" << endl;
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << "Q_mass_u[[" << isp << "]] <- c(";
    for (iyr = styr_pred; iyr <= endyr; iyr++)
     {
      for (iage = 1; iage <= rr_ages(isp); iage++)
       {
        report << Q_mass_u(isp,iyr,iage);
        if (iyr < endyr) report << ", ";
        else if (iage < rr_ages(isp)) report << ", ";
        else report << ")" << endl; 
       }
     }
   }
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << "dim(Q_mass_u[[" << isp << "]]) <- c(";
    report << rr_ages(isp) << ", " << nyrs_pred << ")" << endl;
    report << "Q_mass_u[[" << isp << "]] <- t(Q_mass_u[[" << isp << "]])" << endl;
   }
//#### mass consumed by given predator length l (Q_mass_l) ################
//##
  // 3darray  Q_mass_l(1,nspp_sq2,1,nyrs_pred,1,rr_lens)
  report << "Q_mass_l <- list()" << endl;
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << "Q_mass_l[[" << isp << "]] <- c(";
    for (iyr = styr_pred; iyr <= endyr; iyr++)
     {
      for (iage = 1; iage <= rr_lens(isp); iage++)
       {
        report << Q_mass_l(isp,iyr,iage);
        if (iyr < endyr) report << ", ";
        else if (iage < rr_lens(isp)) report << ", ";
        else report << ")" << endl; 
       }
     }
   }
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << "dim(Q_mass_l[[" << isp << "]]) <- c(";
    report << rr_lens(isp) << ", " << nyrs_pred << ")" << endl;
    report << "Q_mass_l[[" << isp << "]] <- t(Q_mass_l[[" << isp << "]])" << endl;
   }
  report << "rr_ages <- c(";
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << rr_ages(isp);
    if (isp < nspp_sq2) report << ", ";
    else report << ")" << endl;
   }
  report << "rr_lens <- c(";
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << rr_lens(isp);
    if (isp < nspp_sq2) report << ", ";
    else report << ")" << endl;
   }

  //3darray  Q_hat(1,nspp_sq2,1,nyrs_pred,1,rr_lens)
  report << "Q_hat<-list()" << endl;
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << "Q_hat[[" << isp << "]]<-c(";
    for (iyr=styr_pred; iyr <= endyr; iyr++)
     {
      for (icmp=1; icmp<= rr_lens(isp); icmp++)
       {
        report << Q_hat(isp,iyr,icmp);
        if (iyr < endyr) report << ",";
        else if (icmp < rr_lens(isp)) report << ", ";
        else report << ")" << endl;
       }
     }
   }
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << "dim(Q_hat[[" << isp << "]]) <- c(";
    report << rr_lens(isp) << ", " << nyrs_pred << ")" << endl;
    //report << "Q_hat[[" << isp << "]] <- t(Q_hat[[" << isp << "]])" << endl;
   }
  //3darray  T_hat(1,nspp_sq,1,r_lens,1,k_lens)
  report << "T_hat<-list()" << endl;
  for (isp =1; isp <= nspp_sq; isp++)
   {
    report << "T_hat[[" << isp << "]]<-c(";
    for (icmp=1; icmp <= r_lens(isp); icmp++)
     {
      for (r_age=1; r_age<= k_lens(isp); r_age++)
       {
        report << T_hat(isp,icmp,r_age);
        if (icmp < r_lens(isp)) report << ", ";
        else if (r_age < k_lens(isp)) report << ", ";
        else report << ")" << endl;
       }
     }
   }
  for (isp =1; isp <= nspp_sq; isp++)
   {
    report << "dim(T_hat[[" << isp << "]]) <- c(";
    report << k_lens(isp) << ", " << r_lens(isp) << ")" << endl;
    report << "T_hat[[" << isp << "]] <- t(T_hat[[" << isp << "]])" << endl;
   }
  //   init_3darray stoms_w_N(1,nspp,1,l_bins,1,nyrs_stomwts);
  report << "stoms_w_N<-list()" << endl;
  for (isp =1; isp <= nspp; isp++)
   {
    report << "stoms_w_N[[" << isp << "]]<-c(";
    for (icmp=1; icmp<= l_bins(isp); icmp++)
     {
      for (iyr=1; iyr <= nyrs_stomwts(isp); iyr++)
       {
        report << stoms_w_N(isp,icmp,iyr);
        if (icmp < l_bins(isp)) report << ", ";
        else if (iyr < nyrs_stomwts(isp)) report << ", ";
        else report << ")" << endl;
       }
     }
   }

 // init_3darray  stoms_l_N(1,nspp_sq,1,r_lens,1,nyrs_stomlns);
  //rk_sp = 0;
  report << "stoms_l_N<-list()" << endl;
  for (rk_sp =1; rk_sp <= nspp_sq; rk_sp++)
   //for (ksp =1; ksp <= nspp; ksp++)
   {
    //rk_sp = rk_sp + 1;
    report << "stoms_l_N[[" << rk_sp << "]]<-c(";
    for (icmp=1; icmp<= r_lens(rk_sp); icmp++)
     {
      for (iyr=1; iyr <= nyrs_stomlns(rk_sp); iyr++)
       {
        report << stoms_l_N(rk_sp,icmp,iyr);
        if (icmp < r_lens(rk_sp)) report << ", ";
        else if (iyr < nyrs_stomlns(rk_sp)) report << ", ";
        else report << ")" << endl;
       }
     }    
   }
  for (isp =1; isp <= nspp; isp++)
   {
    report << "dim(stoms_w_N[[" << isp << "]]) <- c(";
    report << nyrs_stomwts(isp) << ", " << l_bins(isp) << ")" << endl;
    report << "stoms_w_N[[" << isp << "]] <- t(stoms_w_N[[" << isp << "]])" << endl;
   }
  for (rk_sp =1; rk_sp <= nspp_sq; rk_sp++)
   {
    report << "dim(stoms_l_N[[" << rk_sp << "]]) <- c(";
    report << nyrs_stomlns(rk_sp) << ", " << r_lens(rk_sp) << ")" << endl;
    report << "stoms_l_N[[" << rk_sp << "]] <- t(stoms_l_N[[" << rk_sp << "]])" << endl;
   }

 
  report << "logH_1 <- c(";
  for (isp=1; isp<=nspp_sq; isp++)
   {
    report << logH_1(isp);
    if (isp < nspp_sq) report << ", ";
    else report << ")" << endl;
   }

  //init_3darray diet_w_dat(1,nspp_sq2,1,rr_lens,1,i_wt_yrs_all);
  report << "diet_w_dat<-list()" << endl;
  for (isp=1; isp<=nspp_sq2; isp++)
   {
    report << "diet_w_dat[[" << isp << "]]<-c(";
    for (icmp=1; icmp<=rr_lens(isp); icmp++)
     for (r_age=1; r_age<=i_wt_yrs_all(isp); r_age++)
      {
        report << diet_w_dat(isp,icmp,r_age);
        if (icmp < rr_lens(isp)) report << ", ";
        else if (r_age < i_wt_yrs_all(isp)) report << ", ";
        else report << ")" << endl;
      }
   }
  for (isp =1; isp <= nspp_sq2; isp++)
   {
    report << "dim(diet_w_dat[[" << isp << "]]) <- c(";
    report << i_wt_yrs_all(isp) << ", " << rr_lens(isp) << ")" << endl;
    report << "diet_w_dat[[" << isp << "]] <- t(diet_w_dat[[" << isp << "]])" << endl;
   }

  //init_3darray diet_l_dat(1,nspp_sq,1,r_lens,1,k_lens);
  report << "diet_l_dat<-list()" << endl;
  for (isp=1; isp<=nspp_sq; isp++)
   {
    report << "diet_l_dat[[" << isp << "]]<-c(";
    for (icmp=1; icmp<=r_lens(isp); icmp++)
     for (r_age=1; r_age<=k_lens(isp); r_age++)
      {
        report << diet_l_dat(isp,icmp,r_age);
        if (icmp < r_lens(isp)) report << ", ";
        else if (r_age < k_lens(isp)) report << ", ";
        else report << ")" << endl;
      }
   }
  for (isp =1; isp <= nspp_sq; isp++)
   {
    report << "dim(diet_l_dat[[" << isp << "]]) <- c(";
    report << k_lens(isp) << ", " << r_lens(isp) << ")" << endl;
    report << "diet_l_dat[[" << isp << "]] <- t(diet_l_dat[[" << isp << "]])" << endl;
   }

  //   init_imatrix yrs_stomwts(1,nspp,1,nyrs_stomwts);
  report << "yrs_stomwts<-list()" << endl;
  for (isp=1; isp<=nspp; isp++)
   {
    report << "yrs_stomwts[[" << isp << "]]<-c(";
    for (icmp=1; icmp<=nyrs_stomwts(isp); icmp++)
     {
        report << yrs_stomwts(isp,icmp);
        if (icmp < nyrs_stomwts(isp)) report << ", ";
        else report << ")" << endl;
     }
   }

  //   matrix  omega_hat_ave(1,nspp,1,nages)
  report << endl << "omega_hat_ave<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "omega_hat_ave[[" << isp << "]] <- c(" << endl <<omega_hat_ave(isp,1);
      for (iage=2;iage<=nages(isp);iage++)
        {
              report << ", " << omega_hat_ave(isp,iage);        
        }
          report  << endl << ")" << endl;
    }
  report << endl;

  //   matrix  Q_other_u(1,nspp,1,nages)
  report << endl << "Q_other_u<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "Q_other_u[[" << isp << "]] <- c(" << endl << Q_other_u(isp,1);
      for (iage=2;iage<=nages(isp);iage++)
        {
              report << ", " << Q_other_u(isp,iage);        
        }
          report  << endl << ")" << endl;
    }

 // N_pred_eq(1,nspp,1,nages)
  report << "N_pred_eq <-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_pred_eq[[" << isp << "]] <- c(" << endl << N_pred_eq(isp,1);
      for (iage=2;iage<=nages(isp);iage++)
        {
              report << ", " << N_pred_eq(isp,iage);        
        }
          report  << endl << ")" << endl;
    }  
 // N_prey_eq(1,nspp,1,nages)
  report << "N_prey_eq <-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_prey_eq[[" << isp << "]] <- c(" << endl << N_prey_eq(isp,1);
      for (iage=2;iage<=nages(isp);iage++)
        {
              report << ", " << N_prey_eq(isp,iage);        
        }
          report  << endl << ")" << endl;
    }  
// N_pred_yr(1,nspp,1,nages)
  report << "N_pred_yr <-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_pred_yr[[" << isp << "]] <- c(" << endl << N_pred_yr(isp,1);
      for (iage=2;iage<=nages(isp);iage++)
        {
              report << ", " << N_pred_yr(isp,iage);        
        }
          report  << endl << ")" << endl;
    }  
 // N_prey_yr(1,nspp,1,nages)
  report << "N_prey_yr <-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_prey_yr[[" << isp << "]] <- c(" << endl << N_prey_yr(isp,1);
      for (iage=2;iage<=nages(isp);iage++)
        {
              report << ", " << N_prey_yr(isp,iage);        
        }
          report  << endl << ")" << endl;
    } 
 
//############ N_pred_eqs ################
//  3darray  N_pred_eqs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_pred_eq for all yrs
  report << endl << "N_pred_eqs<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    report << "N_pred_eqs[[" << isp << "]] <- c(" << endl << N_pred_eqs(isp,styr_pred,1);
     for (iyr=styr_pred; iyr <= FinalYr; iyr++)
     {
       for (iage=1;iage<=nages(isp);iage++)  // ages per species
         {
               report << ", " << N_pred_eqs(isp,iyr,iage); 
         }
     }
           report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_pred_eqs[[" << isp << "]] <- N_pred_eqs[[" << isp << "]]";
      report << "[2:length(N_pred_eqs[[ " << isp << "]])]" << endl;
      report << "dim(N_pred_eqs[[" << isp << "]])<-c(" << nages(isp);
      report << ", " << nyrs_pred << ")" << endl;
      report << "N_pred_eqs[[" << isp << "]]<-t(N_pred_eqs[[" << isp << "]])" << endl;
    }
//############ N_prey_eqs ################
//  3darray  N_prey_eqs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_pred_eq for all yrs
  report << endl << "N_prey_eqs<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    report << "N_prey_eqs[[" << isp << "]] <- c(" << endl << N_prey_eqs(isp,styr_pred,1);
     for (iyr=styr_pred; iyr <= FinalYr; iyr++)
     {
       for (iage=1;iage<=nages(isp);iage++)  // ages per species
         {
               report << ", " << N_prey_eqs(isp,iyr,iage); 
         }
     }
           report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_prey_eqs[[" << isp << "]] <- N_prey_eqs[[" << isp << "]]";
      report << "[2:length(N_prey_eqs[[ " << isp << "]])]" << endl;
      report << "dim(N_prey_eqs[[" << isp << "]])<-c(" << nages(isp);
      report << ", " << nyrs_pred << ")" << endl;
      report << "N_prey_eqs[[" << isp << "]]<-t(N_prey_eqs[[" << isp << "]])" << endl;
    }
//############ N_pred_yrs ################
//  3darray  N_pred_yrs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_pred_eq for all yrs
  report << endl << "N_pred_yrs<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    report << "N_pred_yrs[[" << isp << "]] <- c(" << endl << N_pred_yrs(isp,styr_pred,1);
     for (iyr=styr_pred; iyr <= FinalYr; iyr++)
     {
       for (iage=1;iage<=nages(isp);iage++)  // ages per species
         {
               report << ", " << N_pred_yrs(isp,iyr,iage); 
         }
     }
           report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_pred_yrs[[" << isp << "]] <- N_pred_yrs[[" << isp << "]]";
      report << "[2:length(N_pred_yrs[[ " << isp << "]])]" << endl;
      report << "dim(N_pred_yrs[[" << isp << "]])<-c(" << nages(isp);
      report << ", " << nyrs_pred << ")" << endl;
      report << "N_pred_yrs[[" << isp << "]]<-t(N_pred_yrs[[" << isp << "]])" << endl;
    }
//############ N_prey_yrs ################
//  3darray  N_prey_yrs(1,nspp,styr_pred,FinalYr,1,nages)    // save N_pred_eq for all yrs
  report << endl << "N_prey_yrs<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    report << "N_prey_yrs[[" << isp << "]] <- c(" << endl << N_prey_yrs(isp,styr_pred,1);
     for (iyr=styr_pred; iyr <= FinalYr; iyr++)
     {
       for (iage=1;iage<=nages(isp);iage++)  // ages per species
         {
               report << ", " << N_prey_yrs(isp,iyr,iage); 
         }
     }
           report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "N_prey_yrs[[" << isp << "]] <- N_prey_yrs[[" << isp << "]]";
      report << "[2:length(N_prey_yrs[[ " << isp << "]])]" << endl;
      report << "dim(N_prey_yrs[[" << isp << "]])<-c(" << nages(isp);
      report << ", " << nyrs_pred << ")" << endl;
      report << "N_prey_yrs[[" << isp << "]]<-t(N_prey_yrs[[" << isp << "]])" << endl;
    }

//############ Pred_r ################
//##   3darray  Pred_r(1,nspp,1,nyrs_pred,1,nages)
  report << endl << "Pred_r<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    report << "Pred_r[[" << isp << "]] <- c(" << endl << Pred_r(isp,styr_pred,1);
     for (iyr=styr_pred; iyr <= endyr; iyr++)
     {
       for (iage=1;iage<=nages(isp);iage++)  // ages per species
         {
               report << ", " << Pred_r(isp,iyr,iage); 
         }
     }
           report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "Pred_r[[" << isp << "]] <- Pred_r[[" << isp << "]]";
      report << "[2:length(Pred_r[[ " << isp << "]])]" << endl;
      report << "dim(Pred_r[[" << isp << "]])<-c(" << nages(isp);
      report << ", " << nyrs_pred << ")" << endl;
      report << "Pred_r[[" << isp << "]]<-t(Pred_r[[" << isp << "]])" << endl;
    }
//############ Prey_r ################
//##   3darray  Prey_r(1,nspp,1,nyrs_pred,1,nages)
  report << endl << "Prey_r<-list()" << endl;
  for (isp=1;isp<=nspp;isp++)
    {
    report << "Prey_r[[" << isp << "]] <- c(" << endl << Prey_r(isp,styr_pred,1);
     for (iyr=styr_pred; iyr <= endyr; iyr++)
     {
       for (iage=1;iage<=nages(isp);iage++)  // ages per species
         {
               report << ", " << Prey_r(isp,iyr,iage); 
         }
     }
           report  << endl << ")" << endl;
    }
  for (isp=1;isp<=nspp;isp++)
    {
      report << "Prey_r[[" << isp << "]] <- Prey_r[[" << isp << "]]";
      report << "[2:length(Prey_r[[ " << isp << "]])]" << endl;
      report << "dim(Prey_r[[" << isp << "]])<-c(" << nages(isp);
      report << ", " << nyrs_pred << ")" << endl;
      report << "Prey_r[[" << isp << "]]<-t(Prey_r[[" << isp << "]])" << endl;
    }

  //obj_comps
  report << endl << "obj_comps<-vector()" << endl;
  report << "obj_comps <- c(" << obj_comps(1);
  for (icmp=2;icmp<=16;icmp++)
    {
     report << ", " << obj_comps(icmp);
    }
   report << ")" << endl;

  report << endl;



TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(6263);  // replaced 1000 with 6263
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(20000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(15000000);
  arrmblsize=500000000;



