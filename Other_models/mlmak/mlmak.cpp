#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <mlmak.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
 cout << "BEGIN DATA_SECTION" << endl << endl;
  pad_mceval = new ofstream("mceval.dat");
  pad_McFile1 = new ofstream("SSBMCMC.Out");;
  pad_McFile2 = new ofstream("RECMCMC.Out");;
  pad_TempFile = new ofstream("Temp.Out");;
 mcmcmode = 0;
 mcflag   = 1;
 count_iters = 0;
  DebugOut.allocate("DebugOut");
  Set_from_pin_file.allocate("Set_from_pin_file");
 if (Set_from_pin_file != 0)
  cout << "Set from PIN File" << endl;
 else
  cout << "Set from DAT File" << endl;
  ResetPhasesToZero.allocate("ResetPhasesToZero");
  Disc_first_phase.allocate("Disc_first_phase");
  Disc_any_phases.allocate("Disc_any_phases");
  Initial_phase.allocate("Initial_phase");
  Terminal_phase.allocate("Terminal_phase");
  With_Pred.allocate("With_Pred");
  resp_type.allocate("resp_type");
  styr.allocate("styr");
  endyr.allocate("endyr");
 FinalYr = endyr;
  nspp.allocate("nspp");
  comp_type.allocate(1,nspp,"comp_type");
  oldest_age.allocate(1,nspp,"oldest_age");
  l_bins.allocate(1,nspp,"l_bins");
  nages.allocate(1,nspp);
            nages = oldest_age + 1;           // Number of ages (by species)
  nfsh_spp.allocate(1,nspp,"nfsh_spp");
 nfsh = sum(nfsh_spp);                        // Total number of fleets, all species combinations
  spp_fsh.allocate(1,nfsh,"spp_fsh");
  ncomps_fsh.allocate(1,nfsh);
 for(ifsh=1;ifsh<=nfsh;ifsh++)
  {
   isp = spp_fsh(ifsh);
   if (comp_type(isp) == 1) ncomps_fsh(ifsh)=nages(isp);
   else ncomps_fsh(ifsh) = l_bins(isp);
  }
 if (DebugOut >1 ) cout << "ncomps_fsh = " << ncomps_fsh << endl;
  nages_fsh.allocate(1,nfsh);
 for(ifsh=1;ifsh<=nfsh;ifsh++)
  {
   isp = spp_fsh(ifsh);
   nages_fsh(ifsh) = nages(isp);
  }
 if (DebugOut > 1) cout << "nages_fsh = " << nages_fsh << endl;
  styr_rec.allocate(1,nspp);
  endyr_sp.allocate(1,nspp);
  styr_sp.allocate(1,nspp);
 first_rec_est = endyr;
 for(isp=1;isp<=nspp;isp++)
  {
   styr_rec(isp) = (styr - nages(isp)) + 1;
   if (styr_rec(isp) < first_rec_est) first_rec_est = styr_rec(isp);
   endyr_sp(isp) = endyr - 1;
  }
 styr_sp  = styr_rec ;
 if (DebugOut > 1) cout << "styr_rec = " << styr_rec << endl;
 if (DebugOut > 1) cout << "styr_sp = " << styr_sp << endl;
 if (DebugOut > 1) cout << "endyr_sp = " << endyr_sp << endl;
 tot_fsh_yr = (endyr-styr+1)*nfsh;           // Total number of year*fleet*species combinations
 nyrs = endyr-styr+1;  
  catch_bio.allocate(1,nfsh,styr,endyr,"catch_bio");
 if (DebugOut > 2) cout << "Catch_bio = " << endl << catch_bio << endl;
  wt_fsh.allocate(1,nfsh,styr,endyr,1,nages_fsh,"wt_fsh");
 if (DebugOut > 2) cout << "wt_fsh = " << endl << wt_fsh << endl;
  nyrs_fsh_comp.allocate(1,nfsh,"nyrs_fsh_comp");
 tot_yr_fsh_comp = sum(nyrs_fsh_comp);            // Number of years of age data (over fleets)
  yrs_fsh_comp.allocate(1,nfsh,1,nyrs_fsh_comp,"yrs_fsh_comp");
 if (DebugOut > 2) cout << "yrs_fsh_comp = " << endl << yrs_fsh_comp << endl;
  nsmpl_fsh.allocate(1,nfsh,1,nyrs_fsh_comp,"nsmpl_fsh");
 if (DebugOut > 2) cout << "nsmpl_fsh = " << endl << nsmpl_fsh << endl;
  oc_fsh.allocate(1,nfsh,1,nyrs_fsh_comp,1,ncomps_fsh,"oc_fsh");
 if (DebugOut > 2) cout << "oc_fsh = " << endl << oc_fsh << endl;
  offset_fsh.allocate(1,nfsh);
  nsrv_spp.allocate(1,nspp,"nsrv_spp");
 nsrv = sum(nsrv_spp);                            // Total number of survey-species combinations
  spp_srv.allocate(1,nsrv,"spp_srv");
 if (DebugOut > 2) cout << "spp_srv = " << spp_srv << endl;
 tot_srv_yr = (endyr-styr+1)*nsrv;                // Total number of year*survey*species combinations
  ncomps_srv.allocate(1,nsrv);
 for(isrv=1;isrv<=nsrv;isrv++)
  {
   isp = spp_srv(isrv);
   if (comp_type(isp) == 1) ncomps_srv(isrv)=nages(isp);
   else ncomps_srv(isrv) = l_bins(isp);
  }
  nages_srv.allocate(1,nsrv);
 for(isrv=1;isrv<=nsrv;isrv++)
  {
   isp = spp_srv(isrv);
   nages_srv(isrv) = nages(isp);
  }
  nyrs_srv.allocate(1,nsrv,"nyrs_srv");
 tot_yr_srv = sum(nyrs_srv);                      // Number of years of idnex data for surveys (all species)
  yrs_srv.allocate(1,nsrv,1,nyrs_srv,"yrs_srv");
  mo_srv.allocate(1,nsrv,"mo_srv");
  obs_srv.allocate(1,nsrv,1,nyrs_srv,"obs_srv");
 if (DebugOut > 2) cout << "obs_srv = " << endl << obs_srv << endl;
  obs_se_srv.allocate(1,nsrv,1,nyrs_srv,"obs_se_srv");
  nyrs_srv_comp.allocate(1,nsrv,"nyrs_srv_comp");
 tot_yr_srv_age = sum(nyrs_srv_comp);              // Number of years of age data for surveys (all species)
  yrs_srv_comp.allocate(1,nsrv,1,nyrs_srv_comp,"yrs_srv_comp");
  nsmpl_srv.allocate(1,nsrv,1,nyrs_srv_comp,"nsmpl_srv");
  age_matrix.allocate(1,nspp,1,nages);
 for (isp=1; isp <= nspp; isp++)
  for (int a=1; a <= nages(isp); a++)
   age_matrix(isp,a) = double(a);
  max_srv_age.allocate(1,tot_yr_srv_age);
 ipnt = 0;
 for (isrv = 1; isrv <= nsrv; isrv++)
  {
   for (iyr = 1; iyr <= nyrs_srv_comp(isrv); iyr++) 
    { ipnt +=1; max_srv_age(ipnt) = nages(isrv); }
  }
  oc_srv.allocate(1,nsrv,1,nyrs_srv_comp,1,ncomps_srv,"oc_srv");
  wt_srv.allocate(1,nsrv,styr,endyr,1,nages,"wt_srv");
  offset_srv.allocate(1,nsrv);
  wt_pop.allocate(1,nspp,1,nages,"wt_pop");
  maturity.allocate(1,nspp,1,nages,"maturity");
  spawnmo.allocate(1,nspp,"spawnmo");
  spmo_frac.allocate(1,nspp);
 for (isp = 1; isp<=nspp; isp++)
  {
   if (max(maturity(isp)) >.9 ) maturity(isp)/=2.0; // Adjust maturity
   spmo_frac(isp) = (spawnmo(isp)-1)/12.;           // Adjustment for mortality before spawning
  }
  al_key.allocate(1,nspp,1,nages,1,l_bins,"al_key");
  mean_laa.allocate(1,nspp,1,nages);
 if (DebugOut > 1) cout << "al_key = " << endl << al_key << endl;
 nspp_sq = nspp*nspp; // number pred X prey
 nspp_sq2 = nspp*(nspp+1); // number pred X (prey + "other")
  r_lens.allocate(1,nspp_sq);
  k_lens.allocate(1,nspp_sq);
  r_ages.allocate(1,nspp_sq);
  k_ages.allocate(1,nspp_sq);
  rr_lens.allocate(1,nspp_sq2);
  rr_ages.allocate(1,nspp_sq2);
  kk_ages.allocate(1,nspp_sq2);
rk_sp = 0;
for (rsp=1; rsp<=nspp; rsp++)
 for (ksp=1; ksp<=nspp; ksp++)
  {
   rk_sp = rk_sp+1;
   r_lens(rk_sp) = l_bins(rsp);
   k_lens(rk_sp) = l_bins(ksp);
   r_ages(rk_sp) = nages(rsp);
   k_ages(rk_sp) = nages(ksp);
  } 
rk_sp = 0;
for (rsp=1; rsp<=nspp; rsp++)
 for (ksp=1; ksp<=nspp+1; ksp++)
  {
   rk_sp = rk_sp+1;
   rr_lens(rk_sp) = l_bins(rsp);
   rr_ages(rk_sp) = nages(rsp);
   if (ksp <=nspp) kk_ages(rk_sp) = nages(ksp);
   else kk_ages(rk_sp) = 1;
  } 
  pred_l_bin.allocate(1,nspp,1,l_bins,"pred_l_bin");
  omega_vB.allocate(1,nspp,1,nages,"omega_vB");
  omega_sigma.allocate(1,nspp,"omega_sigma");
  nyrs_stomwts.allocate(1,nspp,"nyrs_stomwts");
  nyrs_stomlns.allocate(1,nspp_sq,"nyrs_stomlns");
  yrs_stomwts.allocate(1,nspp,1,nyrs_stomwts,"yrs_stomwts");
  yrs_stomlns.allocate(1,nspp_sq,1,nyrs_stomlns,"yrs_stomlns");
  stoms_w_N.allocate(1,nspp,1,l_bins,1,nyrs_stomwts,"stoms_w_N");
  stoms_l_N.allocate(1,nspp_sq,1,r_lens,1,nyrs_stomlns,"stoms_l_N");
  min_SS_w.allocate(1,nspp,"min_SS_w");
  max_SS_w.allocate(1,nspp,"max_SS_w");
  min_SS_l.allocate(1,nspp_sq,"min_SS_l");
  max_SS_l.allocate(1,nspp_sq,"max_SS_l");
  i_wt_yrs_all.allocate(1,nspp_sq2,"i_wt_yrs_all");
  diet_w_dat.allocate(1,nspp_sq2,1,rr_lens,1,i_wt_yrs_all,"diet_w_dat");
  diet_l_dat.allocate(1,nspp_sq,1,r_lens,1,k_lens,"diet_l_dat");
 if (DebugOut > 1) cout << "nspp_sq = " << nspp_sq << endl;
 if (DebugOut > 1) cout << "pred_l_bin = " << endl << pred_l_bin << endl;
 if (DebugOut > 1) cout << "omega_vB" << endl << omega_vB << endl;
 if (DebugOut > 1) cout << "omega_sigma" << endl << omega_sigma << endl;
 if (DebugOut > 1) cout << "nyrs_stomwts" << endl << nyrs_stomwts << endl;
 if (DebugOut > 1) cout << "yrs_stomwts" << endl << yrs_stomwts << endl;
 if (DebugOut > 1) cout << "stoms_w_N" << endl << stoms_w_N << endl;
 if (DebugOut > 1) cout << "nyrs_stomlns" << endl << nyrs_stomlns << endl;
 if (DebugOut > 1) cout << "yrs_stomlns" << endl << yrs_stomlns << endl;
 if (DebugOut > 1) cout << "yrs_stomwts = " << endl << yrs_stomwts << endl;
 if (DebugOut > 1) cout << "diet_l_dat" << endl << diet_l_dat << endl;
  check.allocate("check");
 cout << ".DAT Check = " << check << endl << endl;
  SrType.allocate("SrType");
  steepnessprior.allocate(1,nspp,"steepnessprior");
  cvsteepnessprior.allocate(1,nspp,"cvsteepnessprior");
  phase_srec.allocate("phase_srec");
  sigmarprior.allocate(1,nspp,"sigmarprior");
  log_sigmarprior.allocate(1,nspp);
 log_sigmarprior = log(sigmarprior);
  cvsigmarprior.allocate(1,nspp,"cvsigmarprior");
  phase_sigmar.allocate("phase_sigmar");
  styr_rec_est.allocate(1,nspp,"styr_rec_est");
  endyr_rec_est.allocate(1,nspp,"endyr_rec_est");
  nrecs_est.allocate(1,nspp);
 n_est_recs = 0;
 for(isp=1; isp<=nspp; isp++)
  {
   nrecs_est(isp) = endyr_rec_est(isp)-styr_rec_est(isp)+1;
   n_est_recs += endyr_rec_est(isp)- styr_rec(isp)+1;
  cout << "styr_rec_est(" << isp <<") = " <<styr_rec_est(isp) << endl;
  }
 cout << "n_est_recs = " << n_est_recs << endl;
  natmortprior.allocate(1,nspp,"natmortprior");
  cvnatmortprior.allocate(1,nspp,"cvnatmortprior");
  natmortphase2.allocate(1,nspp,"natmortphase2");
 NEstNat = 0;
 for (isp=1;isp<=nspp;isp++) if (natmortphase2(isp) > 0) NEstNat += 1;
  qprior.allocate(1,nsrv,"qprior");
  log_qprior.allocate(1,nsrv);
 log_qprior = log(qprior);
  cvqprior.allocate(1,nspp,1,nsrv_spp,"cvqprior");
  phase_q.allocate(1,nsrv,"phase_q");
  q_age_min.allocate(1,nsrv,"q_age_min");
  q_age_max.allocate(1,nsrv,"q_age_max");
  cv_catchbiomass.allocate(1,nspp,"cv_catchbiomass");
  sd_ration.allocate(1,nspp,"sd_ration");
  phase_M.allocate("phase_M");
  phase_Rzero.allocate("phase_Rzero");
  phase_fmort.allocate("phase_fmort");
  phase_fmort1.allocate("phase_fmort1");
  phase_LogRec.allocate("phase_LogRec");
  phase_RecDev.allocate("phase_RecDev");
  phase_SelFshCoff.allocate("phase_SelFshCoff");
  phase_SelSrvCoff.allocate("phase_SelSrvCoff");
  PhasePred1.allocate("PhasePred1");
  PhasePred2.allocate("PhasePred2");
  PhasePred3.allocate("PhasePred3");
 { PhasePredH1a= -1; PhasePredH2 = -1; PhasePredH3 = -1; PhasePredH4 = -1; }
 if (phase_M > 0)
  { phase_M = phase_M - Initial_phase+1; if (phase_M < 1) phase_M = 1; }
 if (phase_Rzero > 0) 
  { phase_Rzero = phase_Rzero - Initial_phase+1; if (phase_Rzero < 1) phase_Rzero = 1; }
 if (phase_fmort > 0)
  { phase_fmort = phase_fmort - Initial_phase+1; if (phase_fmort < 1) phase_fmort = 1; }
 if (phase_fmort1 > 0)
  { phase_fmort1 = phase_fmort1 - Initial_phase+1; if (phase_fmort1 < 1) phase_fmort1 = 1; }
 if (phase_LogRec > 0)
  { phase_LogRec = phase_LogRec - Initial_phase+1; if (phase_LogRec < 1) phase_LogRec = 1; }
 if (phase_RecDev > 0)
  { phase_RecDev = phase_RecDev - Initial_phase+1; if (phase_RecDev < 1) phase_RecDev = 1; }
 if (phase_SelFshCoff > 0)
  { phase_SelFshCoff = phase_SelFshCoff - Initial_phase+1; if (phase_SelFshCoff < 1) phase_SelFshCoff = 1; }
 if (phase_SelSrvCoff > 0)
  { phase_SelSrvCoff = phase_SelSrvCoff - Initial_phase+1; if (phase_SelSrvCoff < 1) phase_SelSrvCoff = 1; }
 if (PhasePred1 > 0)
  { PhasePred1 = PhasePred1 - Initial_phase+1; if (PhasePred1 < 1) PhasePred1 = 1; }
 if (PhasePred2 > 0)
  { PhasePred2 = PhasePred2 - Initial_phase+1; if (PhasePred2 < 1) PhasePred2 = 1; }
 if (PhasePred3 > 0)
  { PhasePred2 = PhasePred3 - Initial_phase+1; if (PhasePred3 < 1) PhasePred3 = 1; }
 if (resp_type > 0) PhasePredH1a = PhasePred2;
 if (1 < resp_type < 7) PhasePredH2 = PhasePred2 + 1;
 if (resp_type > 5)  PhasePredH3 = PhasePred2 + 1;
 if (resp_type == 3 || resp_type == 6) PhasePredH4 = PhasePred2 + 2;
 if (resp_type == 4 || resp_type == 5) PhasePredH3 = PhasePred2 + 2;
 phase_SelSrvCoff2 = phase_SelSrvCoff;
 cout << PhasePred1 << " " << PhasePredH1a << " " << PhasePredH2 << " " << PhasePredH3 << " " << PhasePredH4 << endl;
  catchbiomass_pen.allocate(1,nspp);
 for(isp=1; isp<=nspp; isp++)
  catchbiomass_pen(isp)= 1.0/(2*cv_catchbiomass(isp)*cv_catchbiomass(isp));
  fsh_sel_opt.allocate(1,nfsh,"fsh_sel_opt");
  nselages_in_fsh.allocate(1,nfsh,"nselages_in_fsh");
  phase_sel_fsh.allocate(1,nfsh,"phase_sel_fsh");
  curv_pen_fsh.allocate(1,nfsh,"curv_pen_fsh");
  seldec_pen_fsh.allocate(1,nfsh,"seldec_pen_fsh");
  seldecage.allocate(1,nfsh);
for(ifsh=1; ifsh<=nfsh; ifsh++)  
 { 
  isp =  spp_fsh(ifsh);
  seldecage(ifsh) = int(nages(isp)/2);
 }  
  sel_change_in_fsh.allocate(1,nfsh,styr,endyr,"sel_change_in_fsh");
  n_sel_ch_fsh.allocate(1,nfsh);
  yrs_sel_ch_tmp.allocate(1,nfsh,1,nyrs);
  phase_selcoff_fsh.allocate(1,nfsh);
 phase_selcoff_fsh = phase_sel_fsh;
  srv_sel_opt.allocate(1,nsrv,"srv_sel_opt");
  sel_change_in_srv.allocate(1,nsrv,styr,endyr,"sel_change_in_srv");
  phase_sel_srv.allocate(1,nsrv,"phase_sel_srv");
  sel_slp_in_srv.allocate(1,nsrv);
  sel_inf_in_srv.allocate(1,nsrv);
  nselages_in_srv.allocate(1,nsrv);
  curv_pen_srv.allocate(1,nsrv);
  seldec_pen_srv.allocate(1,nsrv);
  logsel_slp_in_srv.allocate(1,nsrv);
  sel_dinf_in_srv.allocate(1,nsrv);
  n_sel_ch_srv.allocate(1,nsrv);
  yrs_sel_ch_tsrv.allocate(1,nsrv,1,nyrs);
  phase_selcoff_srv.allocate(1,nsrv);
 phase_selcoff_srv = phase_sel_srv;               
 nselages_in_srv = nages-1;
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
 if (DebugOut > 1) cout << "srv_sel_opt = " << srv_sel_opt << endl;
 if (DebugOut > 1) cout << "phase_sel_srv = " << phase_sel_srv << endl;
 if (DebugOut > 1) cout << "curv_pen_srv = " << curv_pen_srv << endl;
 if (DebugOut > 1) cout << "seldec_pen_srv = "  << seldec_pen_srv << endl;
 if (DebugOut > 1) cout << "nselages_in_srv = " << nselages_in_srv << endl;
 Steepness_UB =   .999;                               
  R_guess.allocate(1,nspp);
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
  yrs_sel_ch_fsh.allocate(1,nfsh,1,n_sel_ch_fsh);
  nselages_fsh.allocate(1,nfsh,1,n_sel_ch_fsh);
  yrs_sel_ch_srv.allocate(1,nsrv,1,n_sel_ch_srv);
  nselages_srv.allocate(1,nsrv,1,n_sel_ch_srv);
  for (ifsh=1;ifsh<=nfsh;ifsh++) nselages_fsh(ifsh) = nselages_in_fsh(ifsh);
  for (isrv=1;isrv<=nsrv;isrv++) nselages_srv(isrv) = nselages_in_srv(isrv);
  endyr_all.allocate(1,nspp);
 endyr_all = endyr_sp + 1; 
 cout << endl << "END DATA_SECTION" << endl << endl;
  check2.allocate("check2");
 cout << ".CTR Check = " << check2 << endl << endl;
 styr_pred = 1960; // 1945
 nyrs_pred = endyr-styr_pred+1;  
 PhaseDummy = -1;
 if (With_Pred == 0)
 {PhasePred1 = -3; PhasePred2 = -1; PhasePred3 = -1; }
 if (ResetPhasesToZero==1) 
  {
   cout << "Resetting all phases" << endl;
   PhaseDummy = 1;
   PhasePred1 = -99;
   PhasePred2 = -99;
   PhasePred3 = -99;
   PhasePredH1a= -99; // not in NRM tpl -dhk april 28 09
   PhasePredH2 = -99;
   PhasePredH3 = -99;
   PhasePredH4 = -99;
   phase_M = -99;
   phase_srec = -99;
   phase_Rzero = -99;
   phase_fmort = -99;
   phase_fmort1 = -99;
   for (ifsh=1;ifsh<=nfsh;ifsh++) phase_selcoff_fsh(ifsh) = -99;
   for (isp=1;isp<=nspp;isp++) phase_selcoff_srv(isp) = -99;
   phase_LogRec = -99;
   phase_RecDev = -99;
   phase_SelFshCoff = -99;
   phase_SelSrvCoff = -99;
   phase_SelSrvCoff2 = -99; // not in NRM tpl -dhk april 28 09
  }
 if (phase_SelFshCoff == -99)
  for (ifsh=1;ifsh<=nfsh;ifsh++) phase_selcoff_fsh(ifsh) = -99;
 LowerBoundH3 = -30.0;
 UpperBoundH3 =  -0.000001;
 LowerBoundH4 =  -0.1;
 UpperBoundH4 =   20.0;
 if (resp_type == 7) UpperBoundH3 = -0.00001;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
 cout << "BEGIN PARAMETER_SECTION" << endl << endl;
 int Nselfshpars = 0;
 for (ifsh=1;ifsh<=nfsh;ifsh++)
  for (iyr=1;iyr<=n_sel_ch_fsh(ifsh);iyr++)
   for (iage=1;iage<=nselages_fsh(ifsh,iyr);iage++)
    Nselfshpars += 1;
 cout << "Nselfshpars = " << Nselfshpars << endl;
 int Nselsrvpars = 0;
 int Nselsrvlogs = 0;
 for (isrv=1;isrv<=nsrv;isrv++)
  {
   if (srv_sel_opt(isrv) == 1)
    for (iyUpperBoundH4r=1;iyr<=n_sel_ch_srv(isrv);iyr++)
     for (iage=1;iage<=nselages_srv(isrv,iyr);iage++)
      Nselsrvpars += 1;
   if (srv_sel_opt(isrv) == 2)
      Nselsrvlogs += 1;
  }
 cout << "Nselsrvpars = " << Nselsrvpars << endl;
 cout << "Nselsrvlogs = " << Nselsrvlogs << endl;
 if (Nselsrvlogs == 0)
  { Nselsrvlogs = -1; phase_SelSrvCoff2 = -1; }
 cout << "PHASE VALUES: " << endl;
 cout << "phase_M = " << phase_M << endl;
 cout << "phasePred1 = " << PhasePred1 << endl;
 cout << "phasePred2 = " << PhasePred2 << endl;
 cout << "phasePred3 = " << PhasePred3 << endl;
 cout << "phasePredH1a = " << PhasePredH1a << endl;
 cout << "phasePredH2 = " << PhasePredH2 << endl;
 cout << "phasePredH3 = " << PhasePredH3 << endl;
 cout << "phasePredH4 = " << PhasePredH4 << endl;
 cout << "phase_LogRec = " << phase_LogRec << endl;
 cout << "phase_srec = " << phase_srec << endl;
 cout << "phase_Rzero = " << phase_Rzero << endl;
 cout << "phase_RecDev = " << phase_RecDev << endl;
 cout << "phase_sigmar = " << phase_sigmar << endl;
 cout << "phase_SelFshCoff = " << phase_SelFshCoff << endl;
 cout << "phase_q = " << phase_q << endl;
 cout << "phase_SelSrvCoff = " << phase_SelSrvCoff << endl;
 cout << "phase_SelSrvCoff2 = " << phase_SelSrvCoff2 << endl;
 cout << "phase_fmort = " << phase_fmort << endl << endl;
 cout << "phase_fmort1 = " << phase_fmort1 << endl << endl;
 cout << "phase_selcoff_fsh = " << phase_selcoff_fsh << endl;
 cout << "phase_selcoff_srv = " << phase_selcoff_srv << endl;
 int NFdevs = 0;
 for (ifsh=1;ifsh<=nfsh;ifsh++)
  for (iyr=styr;iyr<=endyr;iyr++)
   if (catch_bio(ifsh,iyr) > 10e-24) NFdevs ++; // 0 in NRM tpl -dhk apr 28 09
 cout << "NFdevs = " << NFdevs << endl;
  MEst.allocate(1,NEstNat,0.02,0.8,phase_M,"MEst");
  log_gam_a.allocate(1,nspp,1.0e-10,19.9,PhasePred1,"log_gam_a");
  log_gam_b.allocate(1,nspp,-5.2,10,PhasePred1,"log_gam_b");
  Q_other_est.allocate(1,nspp,1,nages,PhasePred3,"Q_other_est");
  logH_1.allocate(1,nspp_sq2,PhasePred2,"logH_1");
  logH_1a.allocate(1,nspp,PhasePredH1a,"logH_1a");
  logH_1b.allocate(1,nspp,PhasePredH1a,"logH_1b");
 
  logH_2.allocate(1,nspp_sq,PhasePredH2,"logH_2");
  logH_3.allocate(1,nspp_sq,LowerBoundH3,UpperBoundH3,PhasePredH3,"logH_3");
  H_4.allocate(1,nspp_sq,LowerBoundH4,UpperBoundH4,PhasePredH4,"H_4");
  H_1.allocate(1,nspp_sq2,"H_1");
  #ifndef NO_AD_INITIALIZE
    H_1.initialize();
  #endif
  H_1a.allocate(1,nspp,"H_1a");
  #ifndef NO_AD_INITIALIZE
    H_1a.initialize();
  #endif
  H_1b.allocate(1,nspp,"H_1b");
  #ifndef NO_AD_INITIALIZE
    H_1b.initialize();
  #endif
  H_2.allocate(1,nspp_sq,"H_2");
  #ifndef NO_AD_INITIALIZE
    H_2.initialize();
  #endif
  H_3.allocate(1,nspp_sq,"H_3");
  #ifndef NO_AD_INITIALIZE
    H_3.initialize();
  #endif
  gam_a.allocate(1,nspp,"gam_a");
  #ifndef NO_AD_INITIALIZE
    gam_a.initialize();
  #endif
  gam_b.allocate(1,nspp,"gam_b");
  #ifndef NO_AD_INITIALIZE
    gam_b.initialize();
  #endif
  M.allocate(1,nspp,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  mean_log_rec.allocate(1,nspp,phase_LogRec,"mean_log_rec");
  steepness.allocate(1,nspp,0.21,Steepness_UB,phase_srec,"steepness");
  log_Rzero.allocate(1,nspp,-100,100,phase_Rzero,"log_Rzero");
  rec_dev.allocate(1,n_est_recs,-15,2,phase_RecDev,"rec_dev");
  log_sigmar.allocate(1,nspp,phase_sigmar,"log_sigmar");
  log_selcoffs_fsh.allocate(1,Nselfshpars,phase_SelFshCoff,"log_selcoffs_fsh");
  log_q_srv.allocate(1,nsrv,phase_q,"log_q_srv");
  log_selcoffs_srv.allocate(1,Nselsrvpars,phase_SelSrvCoff,"log_selcoffs_srv");
  logsel_slope_srv_par.allocate(1,Nselsrvlogs,phase_SelSrvCoff2,"logsel_slope_srv_par");
  sel50_srv_par.allocate(1,Nselsrvlogs,phase_SelSrvCoff2,"sel50_srv_par");
  logsel_slope_srv.allocate(1,nsrv,1,n_sel_ch_srv,"logsel_slope_srv");
  #ifndef NO_AD_INITIALIZE
    logsel_slope_srv.initialize();
  #endif
  sel50_srv.allocate(1,nsrv,1,n_sel_ch_srv,"sel50_srv");
  #ifndef NO_AD_INITIALIZE
    sel50_srv.initialize();
  #endif
  fmort_dev_est.allocate(1,NFdevs,-12,8,phase_fmort,"fmort_dev_est");
  log_avg_fmort.allocate(1,nfsh,-10,2,phase_fmort1,"log_avg_fmort");
  dummy.allocate(PhaseDummy,"dummy");
  dummy2.allocate(Terminal_phase,"dummy2");
  natage.allocate(1,nspp,styr_pred,endyr,1,nages,"natage");
  #ifndef NO_AD_INITIALIZE
    natage.initialize();
  #endif
  Sp_Biom.allocate(1,nspp,styr_sp,endyr_all,"Sp_Biom");
  #ifndef NO_AD_INITIALIZE
    Sp_Biom.initialize();
  #endif
  pred_rec.allocate(1,nspp,styr_rec,endyr_all,"pred_rec");
  #ifndef NO_AD_INITIALIZE
    pred_rec.initialize();
  #endif
  mod_rec.allocate(1,nspp,styr_rec,endyr_all,"mod_rec");
  #ifndef NO_AD_INITIALIZE
    mod_rec.initialize();
  #endif
  Z.allocate(1,nspp,styr_pred,endyr,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(1,nfsh,styr,endyr,1,nages_fsh,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(1,nspp,styr_pred,endyr,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  catage.allocate(1,nfsh,styr,endyr,1,nages_fsh,"catage");
  #ifndef NO_AD_INITIALIZE
    catage.initialize();
  #endif
  surv.allocate(1,nspp,"surv");
  #ifndef NO_AD_INITIALIZE
    surv.initialize();
  #endif
  natmort.allocate(1,nspp,"natmort");
  #ifndef NO_AD_INITIALIZE
    natmort.initialize();
  #endif
  rec_dev_spp.allocate(1,nspp,styr_rec,endyr_all,"rec_dev_spp");
  #ifndef NO_AD_INITIALIZE
    rec_dev_spp.initialize();
  #endif
  fmort_dev.allocate(1,nfsh,styr,endyr,"fmort_dev");
  #ifndef NO_AD_INITIALIZE
    fmort_dev.initialize();
  #endif
  Fmort.allocate(1,nspp,styr,endyr,"Fmort");
  #ifndef NO_AD_INITIALIZE
    Fmort.initialize();
  #endif
  m_sigmarsq.allocate(1,nspp,"m_sigmarsq");
  #ifndef NO_AD_INITIALIZE
    m_sigmarsq.initialize();
  #endif
  m_sigmar.allocate(1,nspp,"m_sigmar");
  #ifndef NO_AD_INITIALIZE
    m_sigmar.initialize();
  #endif
  sigmarsq.allocate(1,nspp,"sigmarsq");
  #ifndef NO_AD_INITIALIZE
    sigmarsq.initialize();
  #endif
  sigmar.allocate(1,nspp,"sigmar");
  #ifndef NO_AD_INITIALIZE
    sigmar.initialize();
  #endif
  alpha.allocate(1,nspp,"alpha");
  #ifndef NO_AD_INITIALIZE
    alpha.initialize();
  #endif
  beta.allocate(1,nspp,"beta");
  #ifndef NO_AD_INITIALIZE
    beta.initialize();
  #endif
  Bzero.allocate(1,nspp,"Bzero");
  #ifndef NO_AD_INITIALIZE
    Bzero.initialize();
  #endif
  Rzero.allocate(1,nspp,"Rzero");
  #ifndef NO_AD_INITIALIZE
    Rzero.initialize();
  #endif
  phizero.allocate(1,nspp,"phizero");
  #ifndef NO_AD_INITIALIZE
    phizero.initialize();
  #endif
  avg_rec_dev.allocate(1,nspp,"avg_rec_dev");
  #ifndef NO_AD_INITIALIZE
    avg_rec_dev.initialize();
  #endif
  avgsel_fsh.allocate(1,nfsh,1,n_sel_ch_fsh,"avgsel_fsh");
  #ifndef NO_AD_INITIALIZE
    avgsel_fsh.initialize();
  #endif
  sel_slope_fsh.allocate(1,nfsh,1,n_sel_ch_fsh,"sel_slope_fsh");
  #ifndef NO_AD_INITIALIZE
    sel_slope_fsh.initialize();
  #endif
  log_sel_fsh.allocate(1,nfsh,styr,endyr,1,nages_fsh,"log_sel_fsh");
  #ifndef NO_AD_INITIALIZE
    log_sel_fsh.initialize();
  #endif
  sel_fsh.allocate(1,nfsh,styr,endyr,1,nages_fsh,"sel_fsh");
  #ifndef NO_AD_INITIALIZE
    sel_fsh.initialize();
  #endif
  eac_fsh.allocate(1,nfsh,1,nyrs_fsh_comp,1,nages_fsh,"eac_fsh");
  #ifndef NO_AD_INITIALIZE
    eac_fsh.initialize();
  #endif
  ec_fsh.allocate(1,nfsh,1,nyrs_fsh_comp,1,ncomps_fsh,"ec_fsh");
  #ifndef NO_AD_INITIALIZE
    ec_fsh.initialize();
  #endif
  pred_catch.allocate(1,nfsh,styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  sel_slope_srv.allocate(1,nsrv,1,n_sel_ch_srv,"sel_slope_srv");
  #ifndef NO_AD_INITIALIZE
    sel_slope_srv.initialize();
  #endif
  log_sel_srv.allocate(1,nsrv,styr,endyr,1,nages,"log_sel_srv");
  #ifndef NO_AD_INITIALIZE
    log_sel_srv.initialize();
  #endif
  sel_srv.allocate(1,nsrv,styr,endyr,1,nages,"sel_srv");
  #ifndef NO_AD_INITIALIZE
    sel_srv.initialize();
  #endif
  avgsel_srv.allocate(1,nsrv,1,n_sel_ch_srv,"avgsel_srv");
  #ifndef NO_AD_INITIALIZE
    avgsel_srv.initialize();
  #endif
  pred_srv.allocate(1,nsrv,styr,endyr,"pred_srv");
  #ifndef NO_AD_INITIALIZE
    pred_srv.initialize();
  #endif
  eac_srv.allocate(1,nsrv,1,nyrs_srv_comp,1,nages_srv,"eac_srv");
  #ifndef NO_AD_INITIALIZE
    eac_srv.initialize();
  #endif
  ec_srv.allocate(1,nsrv,1,nyrs_srv_comp,1,ncomps_srv,"ec_srv");
  #ifndef NO_AD_INITIALIZE
    ec_srv.initialize();
  #endif
  q_srv.allocate(1,nsrv,"q_srv");
  #ifndef NO_AD_INITIALIZE
    q_srv.initialize();
  #endif
  gam_ua.allocate(1,nspp_sq,1,r_ages,1,k_ages,"gam_ua");
  #ifndef NO_AD_INITIALIZE
    gam_ua.initialize();
  #endif
  N_pred_eq.allocate(1,nspp,1,nages,"N_pred_eq");
  #ifndef NO_AD_INITIALIZE
    N_pred_eq.initialize();
  #endif
  N_prey_eq.allocate(1,nspp,1,nages,"N_prey_eq");
  #ifndef NO_AD_INITIALIZE
    N_prey_eq.initialize();
  #endif
  N_pred_yr.allocate(1,nspp,1,nages,"N_pred_yr");
  #ifndef NO_AD_INITIALIZE
    N_pred_yr.initialize();
  #endif
  N_prey_yr.allocate(1,nspp,1,nages,"N_prey_yr");
  #ifndef NO_AD_INITIALIZE
    N_prey_yr.initialize();
  #endif
  N_pred_eqs.allocate(1,nspp,styr_pred,FinalYr,1,nages,"N_pred_eqs");
  #ifndef NO_AD_INITIALIZE
    N_pred_eqs.initialize();
  #endif
  N_prey_eqs.allocate(1,nspp,styr_pred,FinalYr,1,nages,"N_prey_eqs");
  #ifndef NO_AD_INITIALIZE
    N_prey_eqs.initialize();
  #endif
  N_pred_yrs.allocate(1,nspp,styr_pred,FinalYr,1,nages,"N_pred_yrs");
  #ifndef NO_AD_INITIALIZE
    N_pred_yrs.initialize();
  #endif
  N_prey_yrs.allocate(1,nspp,styr_pred,FinalYr,1,nages,"N_prey_yrs");
  #ifndef NO_AD_INITIALIZE
    N_prey_yrs.initialize();
  #endif
  Pred_r.allocate(1,nspp,styr_pred,FinalYr,1,nages,"Pred_r");
  #ifndef NO_AD_INITIALIZE
    Pred_r.initialize();
  #endif
  Prey_r.allocate(1,nspp,styr_pred,FinalYr,1,nages,"Prey_r");
  #ifndef NO_AD_INITIALIZE
    Prey_r.initialize();
  #endif
  pred_resp.allocate(1,nspp_sq2,styr_pred,FinalYr,1,rr_ages,1,kk_ages,"pred_resp");
  #ifndef NO_AD_INITIALIZE
    pred_resp.initialize();
  #endif
  Pmort_ua.allocate(1,nspp_sq,styr_pred,FinalYr,1,k_ages,"Pmort_ua");
  #ifndef NO_AD_INITIALIZE
    Pmort_ua.initialize();
  #endif
  Vmort_ua.allocate(1,nspp_sq,1,r_ages,styr_pred,FinalYr,1,k_ages,"Vmort_ua");
  #ifndef NO_AD_INITIALIZE
    Vmort_ua.initialize();
  #endif
  eaten_la.allocate(1,nspp_sq,1,r_lens,styr_pred,FinalYr,1,k_ages,"eaten_la");
  #ifndef NO_AD_INITIALIZE
    eaten_la.initialize();
  #endif
  eaten_ua.allocate(1,nspp_sq,1,r_ages,styr_pred,FinalYr,1,k_ages,"eaten_ua");
  #ifndef NO_AD_INITIALIZE
    eaten_ua.initialize();
  #endif
  Q_mass_l.allocate(1,nspp_sq2,styr_pred,endyr,1,rr_lens,"Q_mass_l");
  #ifndef NO_AD_INITIALIZE
    Q_mass_l.initialize();
  #endif
  Q_mass_u.allocate(1,nspp_sq2,styr_pred,endyr,1,rr_ages,"Q_mass_u");
  #ifndef NO_AD_INITIALIZE
    Q_mass_u.initialize();
  #endif
  omega_hat.allocate(1,nspp,styr_pred,endyr,1,nages,"omega_hat");
  #ifndef NO_AD_INITIALIZE
    omega_hat.initialize();
  #endif
  omega_hat_ave.allocate(1,nspp,1,nages,"omega_hat_ave");
  #ifndef NO_AD_INITIALIZE
    omega_hat_ave.initialize();
  #endif
  Q_hat.allocate(1,nspp_sq2,styr_pred,endyr,1,rr_lens,"Q_hat");
  #ifndef NO_AD_INITIALIZE
    Q_hat.initialize();
  #endif
  Q_other_u.allocate(1,nspp,1,nages,"Q_other_u");
  #ifndef NO_AD_INITIALIZE
    Q_other_u.initialize();
  #endif
  T_hat.allocate(1,nspp_sq,1,r_lens,1,k_lens,"T_hat");
  #ifndef NO_AD_INITIALIZE
    T_hat.initialize();
  #endif
  Zcurr.allocate(1,nspp,1,nages,"Zcurr");
  #ifndef NO_AD_INITIALIZE
    Zcurr.initialize();
  #endif
  Zlast.allocate(1,nspp,1,nages,"Zlast");
  #ifndef NO_AD_INITIALIZE
    Zlast.initialize();
  #endif
 for (ifsh=1;ifsh<=nfsh;ifsh++) nselages_fsh(ifsh)=nselages_in_fsh(ifsh); // Sets all elements of a vector to one scalar value...
 for (isrv=1;isrv<=nsrv;isrv++) nselages_srv(isrv)=nselages_in_srv(isrv); // Sets all elements of a vector to one scalar value...
  sigma.allocate(1,nspp,"sigma");
  #ifndef NO_AD_INITIALIZE
    sigma.initialize();
  #endif
  rec_like.allocate(1,nspp,1,4,"rec_like");
  #ifndef NO_AD_INITIALIZE
    rec_like.initialize();
  #endif
  catch_like.allocate(1,nfsh,"catch_like");
  #ifndef NO_AD_INITIALIZE
    catch_like.initialize();
  #endif
  age_like_fsh.allocate(1,nfsh,"age_like_fsh");
  #ifndef NO_AD_INITIALIZE
    age_like_fsh.initialize();
  #endif
  age_like_srv.allocate(1,nsrv,"age_like_srv");
  #ifndef NO_AD_INITIALIZE
    age_like_srv.initialize();
  #endif
  sel_like_fsh.allocate(1,nfsh,1,4,"sel_like_fsh");
  #ifndef NO_AD_INITIALIZE
    sel_like_fsh.initialize();
  #endif
  sel_like_srv.allocate(1,nsrv,1,4,"sel_like_srv");
  #ifndef NO_AD_INITIALIZE
    sel_like_srv.initialize();
  #endif
  surv_like.allocate(1,nsrv,"surv_like");
  #ifndef NO_AD_INITIALIZE
    surv_like.initialize();
  #endif
  fpen.allocate(1,nspp,1,6,"fpen");
  #ifndef NO_AD_INITIALIZE
    fpen.initialize();
  #endif
  post_priors.allocate(1,4,"post_priors");
  #ifndef NO_AD_INITIALIZE
    post_priors.initialize();
  #endif
  post_priors_srvq.allocate(1,nsrv,"post_priors_srvq");
  #ifndef NO_AD_INITIALIZE
    post_priors_srvq.initialize();
  #endif
  ration_like.allocate(1,nspp,"ration_like");
  #ifndef NO_AD_INITIALIZE
    ration_like.initialize();
  #endif
  diet_like1.allocate("diet_like1");
  #ifndef NO_AD_INITIALIZE
  diet_like1.initialize();
  #endif
  diet_like2.allocate("diet_like2");
  #ifndef NO_AD_INITIALIZE
  diet_like2.initialize();
  #endif
  Zlast_pen.allocate(1,nspp,"Zlast_pen");
  #ifndef NO_AD_INITIALIZE
    Zlast_pen.initialize();
  #endif
  obj_comps.allocate(1,16,"obj_comps");
  #ifndef NO_AD_INITIALIZE
    obj_comps.initialize();
  #endif
  penal_rec_dev.allocate("penal_rec_dev");
  #ifndef NO_AD_INITIALIZE
  penal_rec_dev.initialize();
  #endif
  ration_pen.allocate(1,nspp,"ration_pen");
  #ifndef NO_AD_INITIALIZE
    ration_pen.initialize();
  #endif
  mean_ohat.allocate("mean_ohat");
  #ifndef NO_AD_INITIALIZE
  mean_ohat.initialize();
  #endif
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  ObjTemp.allocate("ObjTemp");
  #ifndef NO_AD_INITIALIZE
  ObjTemp.initialize();
  #endif
  SSLow.allocate("SSLow");
  #ifndef NO_AD_INITIALIZE
  SSLow.initialize();
  #endif
 cout << "END OF PARAMETER_SECTION" << endl << endl;
  SSBOut.allocate(1,nspp,first_rec_est,endyr,"SSBOut");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
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
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
 obj_fun.initialize();
 DoAll();
}

void model_parameters::DoAll(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::gamma_selectivity(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Get_Selectivity(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Get_Mortality(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Get_Bzero(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::AltStart(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Get_Numbers_at_Age(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Compute_Predation(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Get_Survey_Predictions(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Catch_at_Age(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::evaluate_the_objective_function(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Cat_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
  // ============================
  if (DebugOut == 1) cout << "begin Cat_Like" << endl;
  catch_like.initialize();
  for (ifsh=1; ifsh<=nfsh; ifsh++)
   {
    isp = spp_fsh(ifsh);
    catch_like(ifsh) = catchbiomass_pen(isp) * norm2(log(catch_bio(ifsh) +.000001)-log(pred_catch(ifsh) +.000001));
   }
  // ============================
}

void model_parameters::Rec_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Age_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

dvar_vector model_parameters::SRecruit(const dvar_vector& Stmp)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Srv_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Sel_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::ration(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::ration_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::diet_wt_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::diet_len_Like(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Fmort_Pen(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::Compute_priors(void)
{
  ofstream& mceval= *pad_mceval;
  ofstream& McFile1= *pad_McFile1;
  ofstream& McFile2= *pad_McFile2;
  ofstream& TempFile= *pad_TempFile;
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_mceval;
  pad_mceval = NULL;
  delete pad_McFile1;
  pad_McFile1 = NULL;
  delete pad_McFile2;
  pad_McFile2 = NULL;
  delete pad_TempFile;
  pad_TempFile = NULL;
}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(6263);  // replaced 1000 with 6263
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(20000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(15000000);
  arrmblsize=500000000;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
