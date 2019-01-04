// ------------------------------------------------------------------------- //
//                 CEATTLE version 3.1                                       //
//               Multispecies Statistical Model                              //
//          Bioenergetic-based Assessment for Understanding                  //
//              Biomass Linkages To The Environment                          //
//                   Feb. 2015                                               //

// cd /Users/kkari/Dropbox/CEATTLE-master/ceattle-master/main/ceattle.tpl    //
// ------------------------------------------------------------------------- //
//                                                                           //
// AUTHORS:   Kirstin Holsman, Jim Ianelli, Kerim Aydin                      //
//                                                                           //
//                                                                           //
//simMode : simulation mode (can use data file or -opr in command line       
//  
//      1= run in simulation/ operational mode                               //
//      0= run in hindcast/estimation mode                                   //
//msmMode : model mode * indicates working....                               //
//                                                                           //
//      0 = run in single species mode                                       //  
//      1 = run in MSM mode using Jesus's approach (defunct)                 //
//      2 = run in MSM mode using Jesus's approach but new diet and stomach data //
//rationMode : ration mode                                                   //
//                                                                           //
//      0 = use static rations (read in from data)                           //
//      1 = use Elliott & Presson temp rations                               //  
//      2 = (default) use bioenergetic rations (temp and weight) & p-value   //
//                                                                           //
//ntemp_scen : the number of future itrs to run, must >0 to run any -fut     //
//                                                                           //
//recMode : recruitment mode                                                 //
//      0 = project under mean  rec   (no RS)
//      1 = mean RS function (mean RS no covariates, ricker; rs_data4CEATTLE_0_0_0)
//      2 = RS function plus SST (rs_data4CEATTLE_SST
//      3 = RS function based on top AIC selected environm. parms for each spp (rs_data4CEATTLE_TOP)
//      4 = RS function based on model with top R2 value (rs_data4CEATTLE_TopR2)
//      5 = RS function based on Recfile_name (above)   
//      6 = RS function plus bottom temp 
//      7 = RS function plus Cold Pool
//      8 = RS function plus full model plus EAT (tot)
// #harvestMode  : used  for the calcmort  function  in  the projection  / simulation  models  
//      0= hindcast
//      1= project under no fishing,
//      2= project under mean F rate from estimation,
//      3= project to PF (F rate) that matches Btarget,
//      4= fit to NewCatch.dat
//      5= project under F+Fdev radnom seeded fishing rate from hindcast
//      6= project using max F from hindcast
//      7= project using mean catch from hindcast
//      8= project using max catch from hindcast
//      9= project using f profiles
//      10=project under set catch
//      11= project for pollock assessment (project to get target harvest for each sp (or 1) in a given year)
//      12= project using maanda's function an a sloping harvest control rule
//      13= project using amanda's function and a sloping harvest control rule (how do these differ?)
// #ntemp_scen  : number  of  future  iterations
//rand_rec :                                                                 //
//      # 0 = project without error around recruitment                       //
//      # 1 project with randomly drawn errors                               //
//                                                                           //
//             k->species                                                    //
//             i->year                                                       //  
//             l->rep of future iterations                                   //
//             j->misc                                                       //
//             m-> misc                                                      //
//             itemp->reps of future iterations                              //
//                                                                           //
// ------------------------------------------------------------------------  //
// ------------------------------------------------------------------------  //
// need to checK :  Stmp              = mfexp(-Z(sp,yr_ind)*.5);                // note: assumes survey occurs in mid-year       
//        srv_biom_hat(sp,yr) = elem_prod(srv_q(sp) * srv_sel(sp) , elem_prod( Stmp ,N(sp,yr_ind)) ) * wt(sp,yr_ind)(1,nages(sp));

DATA_SECTION

//==============================================================================
// Set up switches
//==============================================================================

//==============================================================================
// Define defaults
//==============================================================================
  int            do_fut;       // run future projections
    !!           do_fut=0; log_input(do_fut)
  int            iseed;
    !!           iseed  = 1234; log_input(iseed)
    !!CLASS random_number_generator r_rep(iseed);
  int            simMode                                         
    !!           simMode=0; log_input(simMode)
  ivector        seed(1,6);
    !!           seed=0   ; log_input(seed)
  number         niter                // number of iterations to run for each loop 
    !!           niter=20  ; log_input(niter)
  number         MNConst              // constant additive for logistic functions
    !!           MNConst = 0.001  ; log_input(MNConst)
  number         cutoff
   !!            cutoff= 0.01;
 

//==============================================================================
// Global parms
//==============================================================================
  int            crn
  int            k
  int            i
  int            itemp

//==============================================================================
// READ IN DATA FILES
//==============================================================================

// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in .ctl file data
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 

 LOCAL_CALCS
    if (ad_comm::argc > 1)
    {
      int on=0;
      if ((on=option_match(ad_comm::argc,ad_comm::argv,"-dump_rep"))>-1)
        dump_rep = atoi(ad_comm::argv[on+1]); 
    }  

    if(dump_rep) cout<<"ad_comm::argc   "<<ad_comm::argc<<endl;

    foldername          = adstring("results/");
    flname              = ad_comm::adprogram_name;
    datafile_name       = adstring("dat_input_files/Ceattle2015.dat");
    diet2file_name      = adstring("dat_input_files/diet2.dat");
    stomfile_name       = adstring("dat_input_files/stomach2015.dat");
    ATFfile_name        = adstring("dat_input_files/ATF_Misc.dat");
    Fprofile_datafile   = adstring("dat_input_files/F_profile.dat");
    Recfile_name        = adstring("fits_4_CEATTLE/rs_data4CEATTLE_0_0_0.dat");
    retrofile_name      = adstring("dat_input_files/retro_data2015.dat.dat");
    futfile_name        = adstring("dat_input_files/projection_data2015.dat");
    catch_in_name       = adstring("dat_input_files/catch_in.dat");
    setC_name           = adstring("dat_input_files/set_catch.dat");
    setF_name           = adstring("dat_input_files/setFabcFofl.dat");

    if (ad_comm::argc > 1)
    {
      int on=0;
      if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ind"))>-1)
      {
        *(ad_comm::global_datafile)>>flname; // where global_datafile is a .ctl file with the names of the .dat files  
        simname=foldername+flname+adstring("_R.rep");
        *(ad_comm::global_datafile)>>datafile_name; 
        *(ad_comm::global_datafile)>>diet2file_name;
        *(ad_comm::global_datafile)>>stomfile_name;
        *(ad_comm::global_datafile)>>Fprofile_datafile;
        *(ad_comm::global_datafile)>>ATFfile_name;
        *(ad_comm::global_datafile)>>Recfile_name;
        *(ad_comm::global_datafile)>>retrofile_name;
        *(ad_comm::global_datafile)>>futfile_name;
        *(ad_comm::global_datafile)>>catch_in_name;
        *(ad_comm::global_datafile)>>setC_name;
        *(ad_comm::global_datafile)>>setF_name;
  
      }
      if ((on=option_match(ad_comm::argc,ad_comm::argv,"-fut"))>-1)
      {
        do_fut=1;
        simname=foldername+flname+adstring("_fut_R.rep");
      } 
    }  
      //new uostream((char*)simname);
       if(dump_rep){
          cout<<"flname     =  "<< flname     <<endl;
          cout<<"foldername =  "<< foldername <<endl;
          cout<<"simname    =  "<< simname    <<endl;
          cout<<"datafile_name   =   "<< datafile_name  <<endl;
          cout<<"stomfile_name   =   "<< stomfile_name  <<endl;
          cout<<"Fprofile_datafile   =   "<< Fprofile_datafile  <<endl;
          cout<<"ATFfile_name    =   "<< ATFfile_name   <<endl;
          cout<<"Recfile_name    =   "<< Recfile_name   <<endl;
          cout<<"retrofile_name    =   "<< retrofile_name   <<endl;
          cout<<"futfile_name    =   "<< futfile_name   <<endl;
          cout<<"catch_in_name    =   "<< catch_in_name   <<endl;
          cout<<"setC_name    =   "<< setC_name   <<endl;
          cout<<"setF_name    =   "<< setF_name   <<endl;
        }
    log_input(flname)
    log_input(foldername)
    log_input(simname)
    log_input(datafile_name)
    log_input(stomfile_name)
    log_input(Fprofile_datafile)
    log_input(ATFfile_name)
    log_input(Recfile_name)
    log_input(retrofile_name)    
    log_input(futfile_name)    
    log_input(catch_in_name) 
    log_input(setC_name) 
    log_input(setF_name) 
 END_CALCS   

// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in .ctl file data
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 
 LOCAL_CALCS
  // if(dump_rep){
      cout<<"============= READ IN CONTROL FILE ============="<<endl;
      cout<<"Control file read in begins:   "<<flname<<endl;
  // }  
 END_CALCS 
  init_int            debugg;// if set to 1 will return debug lines
  !!                   if(dump_rep) cout<<"debugg "<<debugg<<endl;  
  int                 dump_rep;    // print text at end of run
  // !!               dump_rep=0;
  !! dump_rep=debugg;

  init_int            msmMode;                    // 0 is single species, 2 is multispecies; overwritten by command line in ceattke_run.sh
  init_int            nspp                        // number of species CHANGE TO 3
  init_int            recMode                     // recruitment mode for projections; see header for def; overwritten by ceattle_run_futNew.sh
  init_int            rand_rec                    // random recruitment switch for MCMC :  0 project without error around recruitment, 1 project with randomly drawn errors   
  init_int            harvestMode                 // used for the calcmort function in the projection / simulation models (1= project under no fishing, 2= project under mean fishing rate from hindcast, 3= project to F rate that matches Btarget
  init_int            ntemp_scen                  // number  of  future  simulations
   !!                   if(dump_rep) cout<<"ntemp_scen "<<ntemp_scen<<endl;
  init_ivector        simset(1,ntemp_scen);       // set of  future  iterations  to  run e.g., 1 2 3 4 5 6 7 8 9 10                
  init_number         alpha_ABC                   // ADDED min biomass alpha for tier 3 HCR ABC
  init_number         alpha_OFL                   // ADDED min biomass alpha for tier 3 HCR OFL
  init_int            catch_scen                  // ADDED scenario number for amanda's catch function
  init_int            crnum                       // number  of  control rules to  run   
  init_int            ncr_scen                    // number of control rule scenarios
  init_imatrix        c_mult(1,ncr_scen,1,2)      // Control rule scenarios          
  init_int            rationMode;                 // # 0 = use static  rations (read in  from  data); 1 = use Elliott & presson temp  rations (temp and weigt specific); 2 = use bioenergetic  rations (temp and weight) & p-value
   !!                   if(dump_rep) cout<<"rationMode "<<rationMode<<endl;
  init_int            M2mode                      // 0= straight  ratio,  1=  use ADMB  liklihood to  find  M2, 3=  iterate for M2  , 4 uses  jesus's approach  of  suitability       
  init_int            Umode                       // presently defunct :U calculation mode; 0 = use Uobs  read  in  from  data, 1 = Calculate Unew  (only under msmMode==3)
  init_int            maxphase                    // maximum number  of  phases  to  run 
  init_int            RationByAge                 // run bionenergetics on mean length at sp_age (1) or based on lengths converted to ages(0)
  init_number         Btarget                     // target biomass ratio (B/B0) for control rule - ie .40
  init_ivector        B0_set(1,nspp);             // target unfished biomass for control rules 
  init_int            repn2                       // number of reps for finding M2
  init_int            rep_in                      // number  of  reps  for starting  up  from  a #NAME?  command line (should be 5 or more in MSM mode)
  init_int            styr                        // starting year
  init_int            nyrs                        // number of years
  init_int            nyrs_est                    // number of years for estimation (=nyrs unless fitting to a subset)
   !!                   if(dump_rep) cout<<"nyrs_est "<<nyrs_est<<endl;
  init_ivector        logist_sel_phase(1,nspp)    // run logistic selectivity for the survey -2 is off, 2 is on - opposite is the phase sel_phase is calculated in 2
  init_vector         FP_in(1,nspp);              // Input Frates for each species for harvest mode 
  init_int            ctl_test_num                // test read number
  !!                  if (ctl_test_num != 12345) {cout<<"Read .ctl file error"<<endl<<ctl_test_num<<endl;exit(1);}
  int                 endyr                       // end year
  !!                   if(dump_rep) cout<<"Control file read complete "<<endl;
  !!                   if(dump_rep) cout<<"------------------------------"<<endl;

 LOCAL_CALCS
  endyr=styr+nyrs;             // ending year of the model

  // set up random seeds
  if (ad_comm::argc > 1)
  {
    int on=0;
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-opr"))>-1)
    {
      simMode   = atoi(ad_comm::argv[on+1]); 
      do_fut    = 1;
      simname   = foldername+flname+adstring("_fut_R.rep");
  //ntemp_scen = atoi(ad_comm::argv[on+1]);
    }
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-seed"))>-1)
    {
      seed(1) = atoi(ad_comm::argv[on+1]);
      seed(2) = atoi(ad_comm::argv[on+2]);
      seed(3) = atoi(ad_comm::argv[on+3]);
      seed(4) = atoi(ad_comm::argv[on+4]);
      seed(5) = atoi(ad_comm::argv[on+5]);
      seed(6) = atoi(ad_comm::argv[on+6]);
    } 
  }  

  if(dump_rep) cout<<"simMode = "<<simMode<<endl;
  // set estimation or projection modes

  if (ad_comm::argc > 1)
  {
    int on=0;
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-msmMode"))>-1)
      msmMode = atoi(ad_comm::argv[on+1]);
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-recMode"))>-1)
      recMode = atoi(ad_comm::argv[on+1]);
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-rand_rec"))>-1)
      rand_rec = atoi(ad_comm::argv[on+1]);
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-harvestMode"))>-1)
      harvestMode = atoi(ad_comm::argv[on+1]);
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-Fprofiles"))>-1)
    {

      FP_in(1) = atof(ad_comm::argv[on+1]);
      FP_in(2) = atof(ad_comm::argv[on+2]);
      FP_in(3) = atof(ad_comm::argv[on+3]);

    }
    // need to make this flexible for varying numbers of species - OKO
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-crnum"))>-1)
      crnum = atoi(ad_comm::argv[on+1]);
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-Btarget"))>-1)
      Btarget = (atof(ad_comm::argv[on+1]))/10; 

    // Btarget = double(atoi(ad_comm::argv[on+1]))/10; 

  }     
 END_CALCS  

// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in fishery and survey data file e.g., dat_input_files/Ceattle2015.dat
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 

  !!                   ad_comm::change_datafile_name(datafile_name);
 LOCAL_CALCS
  if(dump_rep){
      cout<<"simMode = "<<simMode<<endl;
      cout<<"============= READ IN DATA ============="<<endl;
      cout<<"Main data file read in begins:   "<<datafile_name<<endl;
  }  
 END_CALCS 

  !!                   if(dump_rep) cout<<"nspp = "<<nspp<<endl;
  !!                   if(dump_rep) cout<<"nyrs = "<<nyrs<<endl;
  vector               srv_Mean_CV(1,nspp)                              // annual survey biomass error (SE) Jim
  init_ivector         nages(1,nspp)                                    // sp_age classes vector
  !!                   log_input(nages)
  !!                   if(dump_rep) cout<<"nages = "<<nages<<endl;
  
  // Total catch biomass
  init_ivector         nyrs_tc_biom_obs(1,nspp)                            // number of years with total observed catch   
  !!                   log_input(nyrs_tc_biom_obs)
  init_imatrix         yrs_tc_biom_obs(1,nspp,1,nyrs_tc_biom_obs)               // years with total observed catch   
  !!                   log_input(yrs_tc_biom_obs)
  init_matrix          tc_biom_obs(1,nspp,1,nyrs_tc_biom_obs)               // total observed catch  
  !!                   log_input(tc_biom_obs) 
  
  // Age composition of the catch
  init_ivector         nyrs_fsh_comp(1,nspp);                           // number of years in the fishery sp_age composition data
  !!                   log_input(nyrs_fsh_comp)
  init_imatrix         yrs_fsh_comp(1,nspp,1,nyrs_fsh_comp);            // years for the fishery sp_age composition data
  !!                   log_input(yrs_fsh_comp)
  //!!cout<<yrs_fsh_comp<<endl;
  init_ivector         fsh_age_type(1,nspp)                             // which method of calculating fishery age hat (2 = ATF)  1= catch_hat(sp,yr) / tc_hat(sp,yr); !=1 : fsh_age_hat(sp,yr)  = catch_hat(sp,yr)*age_trans_matrix(sp) / tc_hat(sp,yr);
  !!                   log_input(fsh_age_type)
  init_ivector         fsh_age_bins(1,nspp)                             // bins for fishery age composition data
  !!                   if(dump_rep) cout<<"fsh_age_bins "<<fsh_age_bins<<endl;
  !!                   log_input(fsh_age_bins)
  init_3darray         obs_catch(1,nspp,1,nyrs_fsh_comp,1,fsh_age_bins) // sp_age specific observed catch for plk & pcod, length specific catch for ATF
  !!                   log_input(obs_catch)
    !!                   if(dump_rep) cout<<    "obs_catch (sp=nspp,yr=end) = "<<obs_catch(nspp,nyrs_fsh_comp(nspp))<<endl;
  // Weight at age
  init_ivector         nyrs_wt_at_age(1,nspp);                           // number of years in the weight at age data
  !!                   log_input(nyrs_wt_at_age)
  !!                   if(dump_rep) cout<<"nyrs_wt_at_age "<<nyrs_wt_at_age<<endl;
  init_imatrix         yrs_wt_at_age(1,nspp,1,nyrs_wt_at_age);           // years for the weight at age data
  !!                   log_input(yrs_wt_at_age)
  !!                   if(dump_rep) cout<<"yrs_wt_at_age "<<yrs_wt_at_age<<endl;
  3darray              wt(1,nspp,1,nyrs_wt_at_age,1,nages)               // weight at sp_age in each year (ragged form)
 LOCAL_CALCS  
   for (int sp=1;sp<=nspp;sp++) 
     for(int yr=1;yr<=nyrs_wt_at_age(sp);yr++)
       for (int sp_age=1;sp_age<=nages(sp);sp_age++)
         *(ad_comm::global_datafile) >> wt(sp,yr,sp_age);               // Mean weight at sp_age
 END_CALCS 
 !!                   log_input(wt) 
  // proportion mature & residual mort
  init_matrix          pmature(1,nspp,1,nages)                          // proportion mature by sp_age
  !!                   log_input(pmature)
  init_matrix          M1_base(1,nspp,1,nages)                          // residual mortality
  !!                   log_input(M1_base)
  // survey biomass
  init_ivector         nyrs_srv_biom(1,nspp)                            // Number  of  years of  survey  biomass data  (by number  of  species...)  
  !!                   log_input(nyrs_srv_biom)
  !!                   if(dump_rep) cout<<"nyrs_srv_biom "<<nyrs_srv_biom<<endl;
  init_imatrix         yrs_srv_biom(1,nspp,1,nyrs_srv_biom)             // Years of survey biomass 
  !!                   log_input(yrs_srv_biom)
  !!                   if(dump_rep) cout<<"yrs_srv_biom "<<yrs_srv_biom<<endl;
  init_matrix          srv_biom(1,nspp,1,nyrs_srv_biom)                 // annual survey biomass data
  !!                   log_input(srv_biom)
  init_matrix          srv_biom_se(1,nspp,1,nyrs_srv_biom)              // annual survey biomass error (SE)  
  !!                   log_input(srv_biom_se)
  // survey age composition
  init_ivector         nyrs_srv_age(1,nspp)                             // number of years of survey data sp_age or length composition
  !!                   log_input(nyrs_srv_age)
  !!                   if(dump_rep) cout<<"nyrs_srv_age "<<nyrs_srv_age<<endl;
  init_ivector         srv_age_type(1,nspp)                             // type of composition (sp_age or length) 
  !!                   log_input(srv_age_type)
  !!                   if(dump_rep) cout<<"srv_age_type "<<srv_age_type<<endl;
  init_imatrix         yrs_srv_age(1,nspp,1,nyrs_srv_age)               // years for the survey sp_age comp data
  !!                   log_input(yrs_srv_age)
  init_ivector         srv_age_bins(1,nspp)                             // number of size bins for the survey sp_age comp
  !!                   log_input(srv_age_bins)
  init_matrix          srv_age_n(1,nspp,1,nyrs_srv_age)                 // Sample size for multinomial
  !!                   log_input(srv_age_n)
  init_3darray         srv_age_obs(1,nspp,1,nyrs_srv_age,1,srv_age_bins)// observed sp_age/size compositions
  !!                   log_input(srv_age_obs)
  !!                   if(dump_rep) cout<<"srv_age_obs "<<srv_age_obs<<endl;
  init_matrix          srv_age_sizes(1,nspp,1,srv_age_bins)             // observed size compositions    
  !!                   log_input(srv_age_sizes)
  !!                   if(dump_rep) cout<<"srv_age_sizes "<<srv_age_sizes<<endl;
  init_3darray         age_trans_matrix(1,nspp,1,nages,1,srv_age_bins)       // observed sp_age/size compositions
  !!                   log_input(age_trans_matrix)
  // !!                   if(dump_rep) cout<<"age_trans_matrix "<<age_trans_matrix<<endl;
  
  // Acoustic trawl survey
  init_int             n_eit                                            // number of years with eit data
  !!                   log_input(n_eit)
  !!                   if(dump_rep) cout<<    "n_eit = "<<n_eit<<endl;
  init_ivector          yrs_eit(1,n_eit)                                 // years for available eit data
  !!                   log_input(yrs_eit)
  init_vector          obs_eit(1,n_eit)                                 // observed eit data
  !!                   log_input(obs_eit)
  init_vector          eit_age_n(1,n_eit)                               // number  of  EIT Hauls 2 
  !!                   log_input(eit_age_n)
  !! int eit_age=nages(1);
  init_matrix          obs_eit_age(1,n_eit,1,eit_age)                        // sp_age composition data from eit,
  !!                   log_input(obs_eit_age)
      // !!                   if(dump_rep) cout<<    "obs_eit_age = "<<obs_eit_age<<endl;
  init_int             nyrs_eit_sel
  !!                   log_input(nyrs_eit_sel)
  init_ivector         yrs_eit_sel(1,nyrs_eit_sel)
  !!                   log_input(yrs_eit_sel)
  init_matrix          eit_sel(1,nyrs_eit_sel,1,eit_age)                             // eit selectivity
  !!                   log_input(eit_sel)
      // !!                   if(dump_rep) cout<<    "eit_sel = "<<eit_sel<<endl;
    !!                   if(dump_rep) cout<<    "nyrs_eit_sel = "<<nyrs_eit_sel<<endl;
    !!                   if(dump_rep) cout<<    "yrs_eit_sel = "<<yrs_eit_sel<<endl;
    
  // different mort by sex
  init_ivector         mf_type(1,nspp);                                 // sex specific mort and weight at age? : 1 = same for both, 2 = seperate wt at sp_age for each sex
  !!                   log_input(mf_type)
  !!                   if(dump_rep) cout<<    "mf_type = "<<mf_type<<endl;
  init_3darray         propMorF(1,nspp,1,mf_type,1,nages);              // propMorF
  !!                   log_input(propMorF)
  
  // Von Bert. parms for temp specific weight at age
  init_matrix          t0(1,nspp,1,mf_type);                            // t0 parameter of the temp specific VonB for wt
  !!                   log_input(t0)
  init_matrix          log_mean_d(1,nspp,1,mf_type);                    // log mean d parameter of the temp specific von B for wt
  !!                   log_input(log_mean_d)
  init_matrix          logK(1,nspp,1,mf_type);                          // log k parameter of the temp specific VonB for wt
  !!                   log_input(logK)
  init_matrix          logH(1,nspp,1,mf_type);                          // log H parameter of the temp specific VonB for wt
  !!                   log_input(logH)
  init_matrix          Tcoef(1,nspp,1,mf_type);                         // T coefficent of the linear d equations of the temp specific VonB for wt
  !!                   log_input(Tcoef)
  init_matrix          Pcoef(1,nspp,1,mf_type);                         // P-value coefficent of the linear d equations of the temp specific VonB for wt
  !!                   log_input(Pcoef)
  init_number          test_read;                                       // succesful reading number
  !!                   if (test_read != 12345) {cout<<"Read dat_input_files error"<<endl<<test_read<<endl;exit(1);}
  !!                   if(dump_rep) cout<<    "test_read = "<<test_read<<endl;
  // end of data file
  int                  maxA;
  !!                   maxA=max(nages);
// Set up offsets and penalties
// ------------------------------------------------------------------------- 

LOCAL_CALCS
                      for (int ii =1;ii<=nspp;ii++) 
                        srv_Mean_CV(ii) = mean(elem_div(srv_biom_se(ii),srv_biom(ii)));

                      for (int sp=1;sp<=nspp;sp++) 
                        pmature(sp)=elem_prod(pmature(sp),propMorF(sp,mf_type(sp)));
                      //cout<<"pmature"<<pmature<<endl;
                      for (int sp=1;sp<=nspp;sp++) 
                        for (int yr=1;yr<=nyrs_srv_age(sp);yr++) 
                          srv_age_obs(sp,yr) /= sum(srv_age_obs(sp,yr));
 END_CALCS  


  int                  nselages                                                // number of ages for selectivity (8)
    !!                 nselages=8;  
  ivector              sel_phase(1,nspp)                                       // selectivity phase
  !!                   log_input(sel_phase)
  3darray              fsh_age_obs(1,nspp,1,nyrs_fsh_comp,1,fsh_age_bins)      // obs fishery sp_age composition
  !!                   log_input(fsh_age_obs)
  matrix               tc_obs(1,nspp,1,nyrs_fsh_comp)                          // obs total catch
  !!                   log_input(tc_obs)
  matrix               srv_biom_lse(1,nspp,1,nyrs_srv_biom)                    // survey biomass standard deviations 
    !!                 srv_biom_lse = elem_div(srv_biom_se,srv_biom);          // CV estimation
    !!                 srv_biom_lse = sqrt(log(square(srv_biom_lse) + 1.));   
  vector               offset_fsh(1,nspp)                                      // for multinomial likelihood
  vector               offset_srv(1,nspp)                                      // for multinomial likelihood
  number               offset_eit                                              // for multinomial likelihood
  vector               curv_pen_fsh(1,nspp)                                    // Fishery selectivity penalty
  vector               curv_pen_srv(1,nspp)                                    // Survey selectivity penalty
  int                  tau                                                     // Fishery sample size for multinomial
    !!                 curv_pen_fsh = 12.5;
    !!                 curv_pen_srv = 12.5;                                    // srvy selectivity penalty
    // !!             curv_pen_srv(2)=50;                                      // different penalty for pcod? OKO
    !!                 tau = 200;                                              // Fishery sample size
  !!                   log_input(curv_pen_fsh)
  !!                   log_input(curv_pen_srv)
  !!                   log_input(tau)
 LOCAL_CALCS
                      
                      if(simMode>0||do_fut)
                        niter=20; 
                        // niter=10; 
                      for (int sp=1;sp<=nspp;sp++)                              
                        sel_phase(sp) = -1*logist_sel_phase(sp); //-3

                      offset_srv.initialize();
                      offset_fsh.initialize();
                      offset_eit = 0.;

                      for (int sp=1;sp<=nspp;sp++)
                      {
                        for (int yr=1;yr<=nyrs_fsh_comp(sp);yr++)
                        {
                         
                          tc_obs(sp,yr)= sum(obs_catch(sp,yr));
                          fsh_age_obs(sp,yr) = obs_catch(sp,yr)/(tc_obs(sp,yr)+.01);
                          offset_fsh(sp) -= tau*(fsh_age_obs(sp,yr) + MNConst) * log(fsh_age_obs(sp,yr) + MNConst)  ;
                        }
                        for (int yr=1;yr<=nyrs_srv_age(sp);yr++)
                          offset_srv(sp) -= srv_age_n(sp,yr)*(srv_age_obs(sp,yr) + MNConst) * log(srv_age_obs(sp,yr) + MNConst ) ;
                      }
                      for (int yr=1;yr<=n_eit;yr++)
                      {
                        obs_eit_age(yr) /= sum(obs_eit_age(yr));
                        offset_eit      -= eit_age_n(yr)*(obs_eit_age(yr) + MNConst) * log(obs_eit_age(yr) + MNConst ) ;
                      }
                      
                      log_input(offset_fsh)
                      log_input(srv_biom_lse)
                      log_input(offset_srv)
                      log_input(offset_eit)
                      

 END_CALCS 

  // Unused parms
  // _____________________________________________

     // !!             crulename = "crule_"+ adstring(itoa(Btarget*100,10)) + ".rep";  // not used 
 
  // int            repn;                // number of optimization reps defunct
  // !!             repn=30;                              // misc number of reps for M2 optimization (not used yet)
  // int            iter_nn              // erase? this is re-declared below....
  // vector         trend(1,nyrs)            // trend?
  // !!             trend.fill_seqadd(-1,2/(nyrs-1));         // is this being used?
  //  number         x_f  

// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in Q diet data file e.g., dat_input_files/diet2.dat
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 
  !!                 ad_comm::change_datafile_name(diet2file_name);
  !!                 if(dump_rep) cout<<"Diet2 data read in begins:   "<<diet2file_name<<endl;

  // pre-allocate arrays
  5darray            stomKir(1,nyrs,1,nspp,1,maxA,1,nspp,1,maxA)      // this statement helps to read the stomach data
  4darray            mn_wt_stom(1,nspp,1,maxA,1,nspp,1,maxA)          // this statement helps to read the prey weight in 
  3darray            ration(1,nspp,1,nyrs,1,nages)
  
  // read in data
  init_int           npred3                                           // <unused> number of predator spp for the updated diet2file_name dat file
  init_ivector       nprey3(1,nspp)                                   // <unused> number of prey species for each predator
  init_ivector       n_stomyrs(1,nspp)                                // number years of stomach data 
  init_imatrix       stomyrs(1,nspp,1,n_stomyrs)                      // years of stomach data (1,nspp,1,n_stomyrs)
  4darray            mnQ(1,npred3,1,nyrs+1,1,maxA,1,nspp+1)           // meanwt of each prey spp in the stomach of each predator of sp_age a
  4darray            Qse(1,npred3,1,nyrs+1,1,maxA,1,nspp+1)           // SE wt of each prey spp in the stomach of each predator of sp_age a
  imatrix            stomyr_indx(1,nspp,1,n_stomyrs)
 LOCAL_CALCS
                    for (int sp=1;sp<=nspp;sp++)
                      for (int yr=1;yr<=n_stomyrs(sp);yr++)  
                        stomyr_indx(sp,yr) = stomyrs(sp,yr) - styr+1; 
                    mnQ.initialize();
                      for(int pred=1;pred<=npred3;pred++)  
                        for(int yr=1;yr<=n_stomyrs(pred);yr++)
                          for(int pred_age=1;pred_age<=nages(pred);pred_age++)
                            for(int pp=1;pp<=(nspp+1);pp++) 
                              *(ad_comm::global_datafile) >>  mnQ(pred,stomyr_indx(pred,yr),pred_age,pp); 
                    Qse.initialize();
                      for(int pred=1;pred<=npred3;pred++)  
                        for(int yr=1;yr<=n_stomyrs(pred);yr++)
                          for(int pred_age=1;pred_age<=nages(pred);pred_age++)
                            for(int pp=1;pp<=(nspp+1);pp++) 
                              *(ad_comm::global_datafile) >>  Qse(pred,stomyr_indx(pred,yr),pred_age,pp); 
 END_CALCS  
  init_int           diet2_stomach_test
  !!                 if (diet2_stomach_test != 12345) {cout<<"Read file error "<<diet2file_name<<endl<<diet2_stomach_test<<endl;exit(1);}
    
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in stomach data file e.g., dat_input_files/stomach2015.dat
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 
  !!                 if(dump_rep) cout << "Kir's stomach/diet data read in begins:   "<<stomfile_name<< endl;
  !!                 ad_comm::change_datafile_name(stomfile_name);
  init_vector       other_food(1,nspp)       // other food
  init_ivector       useWt(1,nspp)            // assign relative proportion of prey in the diet according to relative biomass in the system.,otherwise the model with use relative proportion by number
  init_ivector       C_model(1,nspp)          // if ==1, the use Cmax*fT*P
  init_vector        Pvalue(1,nspp)           // This scales the pvalue used if C_model ==1 , proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge - 
  init_vector        Ceq(1,nspp)              // Ceq: which Comsumption equation to use
  init_vector        CA(1,nspp)               // Wt specific intercept of Cmax=CA*W^CB
  init_vector        CB(1,nspp)               // Wt specific slope of Cmax=CA*W^CB
  init_vector        Qc(1,nspp)               // used in fT, QC value
  init_vector        Tco(1,nspp)              // used in fT, thermal optimum
  init_vector        Tcm(1,nspp)              // used in fT, thermal max
  init_vector        Tcl(1,nspp)              // used in fT eq 3, limit
  init_vector        CK1(1,nspp)              // used in fT eq 3, limit where C is .98 max (ascending)
  init_vector        CK4(1,nspp)              // used in fT eq 3, temp where C is .98 max (descending)
  init_int           nTyrs_obs                // <not used>number of temperature years
  init_vector        Tyrs_obs(1,nTyrs_obs)    // <not used> actual years
  init_vector        TempC_obs(1,nTyrs_obs)   // <not used> actual observed bottom temperature from the surveys
   !!                 if(dump_rep) cout << "nTyrs_obs"<<nTyrs_obs<< endl;
  init_vector        S_a(1,nspp)              // #S_a : a,L,L^2,L^3,L^4,L^5 (rows)coef for mean S=a+b*L+b2*L*L, whith a cap at 80cm for each pred spp(cols)"- used for Calc U function <defunct>
  init_vector        S_b(1,nspp)              
  init_vector        S_b2(1,nspp)
  init_vector        S_b3(1,nspp)
  init_vector        S_b4(1,nspp)
  init_vector        S_b5(1,nspp)
  init_vector        aLim(1,nspp)             // #aLim and bLim : upper limit of prey size - used for Calc U function <defunct>
  init_vector        bLim(1,nspp)
  init_vector        aLW(1,nspp)              // aLW : LW a&b regression coef for W=a*L^b  
  init_vector        bLW(1,nspp)              // bLW : LW a&b regression coef for W=a*L^b  
  init_matrix        maxK(1,nspp,1,nspp)      // maxK : nsppxnspp matrix of maximum diet proportions for each predd,preyy combo (used for broken stick expon)
  init_vector        klim(1,nspp)
  init_matrix        prefa(1,nspp,1,nspp)         // prefa of the functional response prefa(pred,prey)
  init_matrix        prefb(1,nspp,1,nspp)         // prefb of the functional response prefb(pred,prey)
  init_matrix        l1(1,nspp,1,nspp)            // l1 of the switch function l1(pred,prey)
  init_matrix        l2(1,nspp,1,nspp)            // l2 of the switch function l2(pred,prey)
   !!                 if(dump_rep) cout << "klim"<<klim<< endl;
  // init_matrix     LA_age(1,nspp,1,nages)       // <unused>length to sp_age coversion matrix
  init_matrix        L2A_matrix(1,nspp,1,nages)      // length to sp_age coversion matrix (rounded to nearest integer)
  init_vector        fdays(1,nspp)                   // number of foraging days for each predator
   !!                 if(dump_rep) cout << "fdays"<<fdays<< endl;
  init_ivector       nlengths(1,nspp)                          // number of Lengths for the matrix to convert ages to lenghths
  !!                 maxL=max(nlengths);
  init_matrix        lengths(1,nspp,1,nlengths)                // Lengths for the matrix to convert ages to lenghths
  init_3darray       A2L_matrix(1,nspp,1,nages,1,nlengths)     // A2L_matrix : Matrix to convert ages to lenghths    
  init_3darray       K(1,nspp,1,nspp,1,nlengths)               // K(1) is npreyXmaxL matrix of stomach proportions of predator 1 
  init_3darray       KAge(1,nspp,1,nspp,1,nages)               // K(1) is npreyXmaxL matrix of stomach proportions of predator 1 
  init_matrix        PAge(1,nspp,1,nages)                      // K(1) is npreyXmaxL matrix of stomach proportions of predator 1 
  init_int           nPyrs   // NOT USED
  !!                 if(dump_rep) cout << "nPyrs"<<nPyrs<< endl;
  !!                 if(dump_rep) cout << "nyrs"<<nyrs<< endl;
  init_3darray       Pby_yr(1,nspp,1,nyrs+1,1,nages)
     !!                 if(dump_rep) cout << "Pby_yr"<<Pby_yr<< endl;
  4darray            Uobs(1,nspp,1,nspp,1,maxL,1,maxL)                      // pred, prey, predL, preyL U matrix (mean number of prey in each pred)
  4darray            UobsWt(1,nspp,1,nspp,1,maxL,1,maxL)                    // pred, prey, predL, preyL U matrix (mean wt of prey in each pred)               
  4darray            UobsAge(1,nspp,1,nspp,1,maxA,1,maxA)                   // pred, prey, predA, preyA U matrix (mean number of prey in each pred age)
  4darray            UobsWtAge(1,nspp,1,nspp,1,maxA,1,maxA)                 // pred, prey, predA, preyA U matrix (mean wt of prey in each pred age)
 
 LOCAL_CALCS
                     Uobs.initialize();
                     for(int predd=1;predd<=nspp;predd++) 
                       for(int preyy=1;preyy<=nspp;preyy++) 
                         for(int predL=1;predL<=nlengths(predd);predL++) 
                           for(int preyL=1;preyL<=nlengths(preyy);preyL++) 
                             *(ad_comm::global_datafile) >>  Uobs(predd,preyy,predL,preyL);
                     UobsWt.initialize();
                     for(int predd=1;predd<=nspp;predd++) 
                       for(int preyy=1;preyy<=nspp;preyy++) 
                         for(int predL=1;predL<=nlengths(predd);predL++) 
                           for(int preyL=1;preyL<=nlengths(preyy);preyL++) 
                             *(ad_comm::global_datafile) >>  UobsWt(predd,preyy,predL,preyL);
 END_CALCS 
  init_matrix        Mn_LatAge(1,nspp,1,nages)                               // <unused>  mean length at age
 LOCAL_CALCS
                     UobsAge.initialize();
                     for(int predd=1;predd<=nspp;predd++) 
                        for(int preyy=1;preyy<=nspp;preyy++) 
                          for(int pred_age=1;pred_age<=nages(predd);pred_age++) 
                            for(int prey_age=1;prey_age<=nages(preyy);prey_age++) 
                               *(ad_comm::global_datafile) >>  UobsAge(predd,preyy,pred_age,prey_age);    
                     UobsWtAge.initialize();
                     for(int predd=1;predd<=nspp;predd++) 
                        for(int preyy=1;preyy<=nspp;preyy++) 
                          for(int predL=1;predL<=nages(predd);predL++) 
                            for(int preyL=1;preyL<=nages(preyy);preyL++) 
                              *(ad_comm::global_datafile) >>  UobsWtAge(predd,preyy,predL,preyL); 
 END_CALCS 
  init_int           test_stomach
  !!                 if (test_stomach != 12345) {cout<<"Read file error test_stomach"<<endl<<test_stomach<<endl;exit(1);}
  ivector            prey_nlengths(1,nspp)
  ivector            prey_nages(1,nspp)
  int                maxL
  3darray            other_food_byAge(1,nspp,1,nspp,1,nages)         // <not used?> K(1) is npreyXmaxL matrix of stomach proportions of predator 1 
 LOCAL_CALCS
                     for(int preyy=1;preyy<=nspp;preyy++)
                     {
                      prey_nlengths(preyy)=nlengths(preyy);
                      prey_nages(preyy)=nages(preyy);
                     }
 END_CALCS 

    //init_matrix   Kglm_a(1,nspp,1,nspp)
    //init_matrix   Kglm_b(1,nspp,1,nspp)
    //init_3darray  K_Coef(1,nspp,1,nspp,1,6)         // spp pref curves fit using a gam in R K_Coef(predd,preyy,6) 


// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in Fprofile data file e.g., dat_input_files/F_profile.dat
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------  

  !!                 if(dump_rep) cout << "F_profiledata read in begins:   "<<Fprofile_datafile<< endl;  
  !!                 ad_comm::change_datafile_name(Fprofile_datafile);
  init_int           n_f;
  init_vector        Frates(1,n_f);
  init_int           np;
  init_matrix        Fprofiles(1,np,1,nspp);
  init_int           testnumFp;  
  !!                 if (testnumFp != 12345) {cout<<"Read file error Fprofile test number"<<endl<<testnumFp<<endl;exit(1);}


// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in retrospective hindcast data file e.g., dat_input_files/retro_data2016_scaled_extnd.dat
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------   

  !!                 if(dump_rep) cout<<"read in "<<retrofile_name<<endl;
  !!                 ad_comm::change_datafile_name(retrofile_name);  
  init_int           nspp2
  !! if(nspp2!=nspp){cout<<"number of species does not match"<<endl;exit(1);} 
  init_int           nTyrs                      // number years for the hindcast data 
    !!                 if(dump_rep) cout << "nTyrs= "<<nTyrs<<endl;
    !!                 if(dump_rep) cout << "SSTempC_retro= "<<SSTempC_retro<<endl;
  
  int                nyrs_use                   
  !!                 nyrs_use=min(nTyrs,nyrs);  // minimum number of years between the two - if hindcast is longer use nyrs, and vise versa
  init_vector        Tyrs(1,nTyrs)               // years for the hindcast data , doesn't have to match model years
  !!                 if(dump_rep) cout << "Tyrs= "<<Tyrs<<endl;
  init_int           ncov                       // number of covariates
  init_int           nEatcovs                   // number of covariates for the food limitation configuration of RS (demand-zoop*Eatcov)
  init_ivector       Eat_covs(1,nEatcovs)       // index of covariates that should be used for the function
  matrix             rs_cov(1,ncov,1,nTyrs);    // covariate values (may be scaled)
 LOCAL_CALCS
                     for (int c=1;c<=ncov;c++) 
                      for(int yr=1;yr<=nTyrs;yr++)
                        *(ad_comm::global_datafile) >> rs_cov(c,yr);  // covariate values
 END_CALCS    
 //!! cout<< rs_cov(ncov) <<endl;
  init_vector        BTempC_retro(1,nTyrs)                              // Bottom Temperature (actual values)
  init_vector        SSTempC_retro(1,nTyrs)                             // Sea surface Temperatuer (actual values)

  !!                 if(dump_rep) cout << "SSTempC_retro= "<<SSTempC_retro<<endl;
  3darray            overlap_dat(1,nspp,1,nspp,1,nTyrs)               // pred prey overlap (scenario,pred,prey, yr); pred overlap with prey species (2-6/0-6 area)  
 LOCAL_CALCS
                     for (int pred=1;pred<=nspp;pred++) 
                       for(int prey=1;prey<=nspp;prey++)
                         for(int yr=1;yr<=nTyrs;yr++)
                           *(ad_comm::global_datafile) >> overlap_dat(pred,prey,yr);          // overlap values
 END_CALCS 
  init_int          test_numZ
  !!                if(test_numZ!=12345){cout<<"ERROR RETROSPECTIVE DATA ("<<retrofile_name <<") WAS NOT LOADED CORRECTLY "<<test_numZ<<endl;exit(1);}



 // read in F40 and FABC for each sp, scenario, yr F40(itemp,k,i)
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in projection data file e.g., dat_input_files/projection_data2016_scaled_extnd.dat
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------   
  !!                if(dump_rep) cout << "Projection data read in begins:   "<<futfile_name<< endl;
  !!                ad_comm::change_datafile_name(futfile_name);  
  init_int          nEatcovs_fut                                      // number of covariates for the food limitation configuration of RS (demand-zoop*Eatcov)
  init_ivector      Eat_covs_fut(1,nEatcovs_fut)                      // index of covariates that should be used for the function
  init_int          num_runs                                          // number of future temperature simulations
  init_int          nyrs_fut                                          // (1, num_runs) number years for the projection data
  init_vector       fut_yrs(1,nyrs_fut)                               // fut_years : matrix(num_runs,nyrs_fut) years for the projection data
  init_int          ncov_fut                                          // number of RS covariates
  3darray           rs_cov_fut(1,num_runs,1,ncov_fut,1,nyrs_fut)      // covariate values (may be scaled)
 LOCAL_CALCS
                    if(ncov_fut!=ncov){cout<<"ERROR ncov and ncov_fut do not match!!!"<<endl;exit(1);}
                    for (int c=1;c<=nEatcovs_fut;c++)  
                    {
                       if(Eat_covs_fut(c)!=Eat_covs(c))
                       {
                           cout<<"ERROR Eat_covs("<<c<<") and Eat_covs_fut do not match!!!"<<endl;
                           exit(1);
                       }
                     }
                    for (int c=1;c<=ncov_fut;c++) 
                      for (int cm=1;cm<=num_runs;cm++)
                        for(int yr=1;yr<=nyrs_fut;yr++)
                          *(ad_comm::global_datafile) >> rs_cov_fut(cm,c,yr);       
 END_CALCS   
  init_matrix       BTempC_fut_all(1,num_runs,1,nyrs_fut)              // Bottom temperature for each scenario (actual Temp for bioenergetics)
  init_matrix       SSTempC_fut_all(1,num_runs,1,nyrs_fut)             // Sea surface temperature (actual temperature)
  4darray           overlap_fut(1,num_runs,1,nspp,1,nspp,1,nyrs_fut);  // pred prey overlap (scenario,pred,prey, yr); pred overlap with prey species (2-6/0-6 area)  
 LOCAL_CALCS
                    for (int pred=1;pred<=nspp;pred++) 
                      for (int prey=1;prey<=nspp;prey++) 
                       for (int cm=1;cm<=num_runs;cm++)
                         for(int yr=1;yr<=nyrs_fut;yr++)
                           *(ad_comm::global_datafile) >> overlap_fut(cm,pred,prey,yr);       
 END_CALCS 
  init_int         test_num_projection
  !!               if (test_num_projection != 12345) {cout<<"Read file error projection data"<<endl<<test_num_projection<<endl;exit(1);}

// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 
// If harvest mode =10 (and the catch is set to read in as a constant value)
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
  // seems duplicative with F =11, maybe merge the two?
  !!                ad_comm::change_datafile_name(catch_in_name);     //where catch_in_name is the datafile with the new catch for each spp
  init_matrix             tmpCatch_in(1,2,1,nspp);
  init_int          testreadcatch;
  !!                if (testreadcatch!=12345) {cout<<"Error reading in "<<catch_in_name<<" test read = "<<testreadcatch<<endl;exit(1);}


      //tmpC=(sqrt( norm2((tc_hat(k)-mean(tc_biom_hat(k))))/nyrs )*rmultF(k,i)+mean(tc_hat(k)));
 


// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// If harvestMode =11 Read itarget catches or Frates
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------   

  // !!                 ad_comm::change_datafile_name("dat_input_files/set_catch.dat");
  !!                 ad_comm::change_datafile_name(setC_name);
  !!                cout<<"reading in "<<setC_name<<endl;
      vector set_catch(1,nspp); 
      vector set_val(1,nspp);   
      !!set_val=-99; 
      int test_numsc;
      // ifstream tt_tmp("../../dat_input_files/set_catch.dat");
      // tt_tmp(1,nspp)>>set_catch;
      // tt_tmp(nspp+1,nspp*2)>>set_val;
      // tt_tmp.close();  

 LOCAL_CALCS
   if(harvestMode==11){
      for (int sp=1;sp<=nspp;sp++) 
          *(ad_comm::global_datafile) >> set_catch(sp); 
      for (int sp=1;sp<=nspp;sp++) 
          *(ad_comm::global_datafile) >> set_val(sp); 
      *(ad_comm::global_datafile) >> test_numsc; 
      cout<<"set_val"<<set_val<<endl;
      if(test_numsc!=12345){cout<<"ERROR "<<setC_name<<" did not load correctly "<<test_numsc<<endl;exit(1);}

    }
 END_CALCS  

// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// If harvestMode =12 Read F40 and F35 - read in set F40 value 
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------   
  3darray        F40(1,ntemp_scen,1,nspp,1,nyrs_fut);                          // Future FAbC proxy / target rate
  3darray        F35(1,ntemp_scen,1,nspp,1,nyrs_fut);                          // Future Fofl proxy/ min target
  int  Fnyrs;
  int test_numsf;

     // !!            F40(1)=0.5589386;
   // !!            F40(2)=0.3456143;
   // !!            F40(3)=0.08429076;
 
 LOCAL_CALCS

  // read in target F rate
  if(harvestMode==12|harvestMode==13|harvestMode==15){
    // set the data file to setF_name ; eg setFabcFofl.dat
    ad_comm::change_datafile_name(setF_name);  
     *(ad_comm::global_datafile) >> Fnyrs;

    if(Fnyrs!=nyrs_fut) {cout<<"ERROR, Fset n future years ("<<Fnyrs<<") is not equal to nyrs_fut ("<<nyrs_fut<<") "<<endl;exit(1);}
    for(itemp=1;itemp<=ntemp_scen;itemp++)
      for (int sp=1;sp<=nspp;sp++)
        for (int yr=1;yr<=nyrs_fut;yr++)
            *(ad_comm::global_datafile) >>F40(itemp,sp,yr);
    for(itemp=1;itemp<=ntemp_scen;itemp++)
      for (int sp=1;sp<=nspp;sp++)
        for (int yr=1;yr<=nyrs_fut;yr++)
            *(ad_comm::global_datafile) >>F35(itemp,sp,yr);
      *(ad_comm::global_datafile) >>   test_numsf;
      if(test_numsf!=12345){cout<<"ERROR "<<setF_name<<" did not load correctly "<<test_numsf<<endl;exit(1);}

  }
 END_CALCS  
 
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Read in RS parameter (recfilename) files from the recruitment fitting function (for projections)  e.g.,fits_4_CEATTLE/rs_data4CEATTLE_full_avg.dat
// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------    
  // pre allocate values
  !!                if(dump_rep) cout << "Read in RS parameters"<< endl;
  ivector            ncovs(1,nspp);
  !!                 for (int pred=1;pred<=nspp;pred++) ncovs(pred)=ncov;
  vector             BevH(1,nspp);
  vector             BLM(1,nspp);
  vector             LM(1,nspp);
  vector             logsigma_rs(1,nspp);
  vector             logsigma_rs_std(1,nspp);
  vector             mnRat(1,nspp);
  vector             sdRat(1,nspp);
  vector             log_aa_c(1,nspp);
  vector             log_aa_c_std(1,nspp);
  vector             noRation(1,nspp);     // ? what is this?
  int                rs_test_num;
  vector             log_bb_c(1,nspp);
  vector             log_bb_c_std(1,nspp);
  vector             rs_f(1,nspp);     // ? what is this?
  number             obj_fun_rs;
 LOCAL_CALCS
   if(simMode>0||do_fut)
   {
    if(dump_rep) cout<<"msmMode = "<<msmMode<<endl;
    if(dump_rep) cout<<"rand_rec = "<<rand_rec<<endl;
    if(dump_rep) cout<<"recMode = "<<recMode<<endl;
    if(dump_rep) cout<<"harvestMode = "<<harvestMode<<endl;
    if(dump_rep) cout<<"ntemp_scen = "<<ntemp_scen<<endl;
    switch (recMode)
    {    
         default:
            ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_0_0_0.dat"); 
            cout<<"           Rec_file = fits_4_CEATTLE/rs_data4CEATTLE_0_0_0.dat"<<endl;
          break;
          case 1:
            ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_0_0_0.dat"); 
            cout<<"           Rec_file = fits_4_CEATTLE/rs_data4CEATTLE_0_0_0.dat"<<endl;
            // Ricker RS no covariates
          break;
          case 2:
            ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_SST.dat"); 
            cout<<"           Rec_file = fits_4_CEATTLE/rs_data4CEATTLE_SST.dat"<<endl;
            // SST only
          break;
          case 3:
            ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_TOP.dat"); 
            cout<<"           Rec_file = fits_4_CEATTLE/rs_data4MSM_TOP.dat"<<endl;
            // top selected model based on AIC
          break;
          case 4:
               ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_TopR2.dat"); 
             cout<<"           Rec_file = fits_4_CEATTLE/rs_data4CEATTLE_TopR2.dat"<<endl;
              // top selected model based on AIC
          break;
          case 5:
            ad_comm::change_datafile_name(Recfile_name); 
            cout<<"           Rec_file = "<<Recfile_name<<endl;
          break;
          case 6:
            ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_BT.dat"); 
            cout<<"           Rec_file = fits_4_CEATTLE/rs_data4CEATTLE_BT.dat"<<endl;
            // Ricker RS no covariates
          break;
          case 7:
            ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_ColdPool.dat"); 
            cout<<"           Rec_file = fits_4_CEATTLE/rs_data4CEATTLE_ColdPool.dat"<<endl;
            // Ricker RS no covariates
          break;
          case 8:
            ad_comm::change_datafile_name("fits_4_CEATTLE/rs_data4CEATTLE_fullEAT_tot.dat"); 
            cout<<"           Rec_file = fits_4_CEATTLE/rs_data4CEATTLE_fullEAT_tot.dat"<<endl;
            // Ricker RS no covariates
          break;


    }
          cout<<"########################################################"<<endl;
          cout<<""<<endl;
          cout<<"Setting up the model, please wait...";

     *(ad_comm::global_datafile) >>  ncovs;
     *(ad_comm::global_datafile) >>  BevH;
     *(ad_comm::global_datafile) >>  BLM;
     *(ad_comm::global_datafile) >>  LM;
     *(ad_comm::global_datafile) >>  logsigma_rs;
     *(ad_comm::global_datafile) >>  logsigma_rs_std;
     *(ad_comm::global_datafile) >>  mnRat;
     *(ad_comm::global_datafile) >>  sdRat;
     *(ad_comm::global_datafile) >>  log_aa_c;
     *(ad_comm::global_datafile) >>  log_aa_c_std;
     *(ad_comm::global_datafile) >>  log_bb_c;
     *(ad_comm::global_datafile) >>  log_bb_c_std;
   }  
   for (int pred=1;pred<=nspp;pred++) 
   {   if (ncovs(pred)!=ncov){
        cout<<"ERROR! ncov from Rec file does not match ncov from hindcast "<<endl;
      }
    }

   if (debugg) cout<< "ncov " <<ncov<<endl;
 END_CALCS
  matrix             cov_phase2(1,nspp,1,ncovs);
  matrix             cov_type2(1,nspp,1,ncovs);
  matrix             rs_parm(1,nspp,1,ncovs);
  matrix             rs_parm_std(1,nspp,1,ncovs);
 LOCAL_CALCS
   if(simMode>0||do_fut){
      for (int sp=1;sp<=nspp;sp++)
         for (int c=1;c<=ncovs(sp);c++)
            *(ad_comm::global_datafile) >>  rs_parm(sp,c);  
      for (int sp=1;sp<=nspp;sp++)
         for (int c=1;c<=ncovs(sp);c++)
            *(ad_comm::global_datafile) >>  rs_parm_std(sp,c);  
     for (int sp=1;sp<=nspp;sp++)
         for (int c=1;c<=ncovs(sp);c++)
            *(ad_comm::global_datafile) >>  cov_type2(sp,c);  
      for (int sp=1;sp<=nspp;sp++)
         for (int c=1;c<=ncovs(sp);c++)
            *(ad_comm::global_datafile) >>  cov_phase2(sp,c);      
      *(ad_comm::global_datafile) >>  rs_test_num;
      if(rs_test_num!=12345)cout<<" ERROR DATA WAS NOT LOADED CORRECTLY "<<rs_test_num<<endl;
      if(rs_test_num!=12345)exit(1);
      
      if(dump_rep) cout<<"logsigma_rs "<<logsigma_rs<<endl;
      if(dump_rep) cout<<"noRation "<<noRation<<endl;
      if(dump_rep) cout<<"log_aa_c "<<log_aa_c<<endl;
      if(dump_rep) cout<<"log_bb_c "<<log_bb_c<<endl;
      if(dump_rep) cout<<"ncov "<<ncov<<endl;
      if(dump_rep) cout<<"cov_phase2 "<<cov_phase2<<endl;
      if(dump_rep) cout<<"rs_parm "<<rs_parm<<endl;
   }
 END_CALCS


// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 
// Extra crap for ATF Jim and also Pcod...
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 
  !!                ad_comm::change_datafile_name(ATFfile_name);  
  int               itmp ;
  int               ktmp ;
  !!                itmp = nages(3);
  !!                ktmp = srv_age_bins(3);
  init_vector       ATF_SexRatio(1,itmp);
  init_matrix       ATF_agesize_F(1,itmp,1,ktmp);
  init_matrix       ATF_agesize_M(1,itmp,1,ktmp);
  !!                itmp = nages(2);
  !!                ktmp = srv_age_bins(2);
  init_matrix       PCOD_agesize(1,itmp,1,ktmp)

!!if(dump_rep) cout<<"============= END: READ IN DATA ============="<<endl;  


//==============================================================================
// Now reformat and allocate some data and variables
//==============================================================================

  // pre allocate  Bioenergetics parameters
    matrix            TempC_futUSE(1,ntemp_scen,1,nyrs_fut)         // rows with future temperature scenarios for bottom temp - extrapolated linearly from SST projection
    matrix            fT(1,nspp,1,nTyrs)                            // pre-allocate temperature function of consumption - only used if C_model==1
    3darray           fT_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)        // pre-allocate temperature function of consumption - only used if C_model==1
    vector            TempC(1,nTyrs);
    vector            TempC_fut(1,nyrs_fut);
  // pre allocate U matrix stuff  
    matrix            S2(1,nspp,1,nlengths)                         // pre-allocate mean stomach weight as a function of lengths
    3darray           S2Age(1,nspp,1,nyrs,1,nages)                  // pre-allocate mean stomach weight as a function of sp_age
    matrix            Limit(1,nspp,1,nlengths)                      // pre-allocate
    3darray           LimitAge(1,nspp,1,nyrs,1,nages)               // pre-allocate
    4darray           LimitAge2(1,nspp,1,nspp,1,maxA,1,maxA)        // pre-allocate
    number            totalpreyL
    number            totalpredL
    number            totalpreyA
  
    matrix            WtbyAge(1,nspp,1,srv_age_bins)                 // currently not used
    matrix            WtbyL(1,nspp,1,nlengths)  
    5darray           swtch(1,nyrs,1,nspp,1,nspp,1,maxA,1,maxA)      // swtch(pred,prey,pred_age,prey_age)
    
    vector            years(1,nyrs);
    3darray           LbyAge(1,nspp,1,nyrs,1,nages)                  // Length by age from LW regression
    matrix            LbyAge2(1,nspp,1,nages)                        // Length by age using A2L matrix
    !!                if (debugg) cout<< "maxA "<<maxA <<endl;
 LOCAL_CALCS
    for(int spp=1;spp<=nspp;spp++)
      for(int sp_age=1;sp_age<=nages(spp);sp_age++)
        LbyAge2(spp,sp_age)=sum(elem_prod(lengths(spp),A2L_matrix(spp,sp_age)));
    for(int spp=1;spp<=nspp;spp++){
      for(int yr=1;yr<=nyrs;yr++)
        LbyAge(spp,yr)=(pow((1/aLW(spp)),(1/bLW(spp))))*pow((wt(spp,yr)*1000),(1/bLW(spp)));
        double AvgWt_tmp; 
        AvgWt_tmp=0.;      
        for(int LL=1;LL<=nlengths(spp);LL++)
          WtbyL(spp,LL)=aLW(spp)*(pow(lengths(spp,LL),bLW(spp))); //in g ?      
        for(int LL=1;LL<=nages(spp);LL++)
         WtbyAge(spp,LL)=aLW(spp)*(pow(L2A_matrix(spp,LL),bLW(spp))); //in g divide by 1000 to match wt // not used delete?
    }
    for (int predd=1;predd<=nspp;predd++)  
      totalpredL+=nlengths(predd);
    for (int preyy=1;preyy<=nspp;preyy++)
    {
      totalpreyL+=nlengths(preyy);
      totalpreyA+=nages(preyy);
    }

    for (int predd=1;predd<=nspp;predd++)
    {
      for (int sp_length=1;sp_length<=nlengths(predd);sp_length++)
      {
        S2(predd,sp_length)=S_a(predd)+(S_b(predd)*lengths(predd,sp_length))+(S_b2(predd)*pow(lengths(predd,sp_length),2))+(S_b3(predd)*pow(lengths(predd,sp_length),3))+(S_b4(predd)*pow(lengths(predd,sp_length),4))+(S_b5(predd)*pow(lengths(predd,sp_length),5));  
        if(lengths(predd,sp_length)>80)
         S2(predd,sp_length)=S_a(predd)+(S_b(predd)*80)+(S_b2(predd)*pow(80,2))+(S_b3(predd)*pow(80,3))+(S_b4(predd)*pow(80,4))+(S_b5(predd)*pow(80,5)); // set everything above 80 to 80
      }
      for(int yr=1;yr<=nyrs;yr++)
      {
        for (int sp_age=1;sp_age<=nages(predd);sp_age++)
        {
          S2Age(predd,yr,sp_age)=S_a(predd)+(S_b(predd)*LbyAge(predd,yr,sp_age))+(S_b2(predd)*pow(LbyAge(predd,yr,sp_age),2))+(S_b3(predd)*pow(LbyAge(predd,yr,sp_age),3))+(S_b4(predd)*pow(LbyAge(predd,yr,sp_age),4))+(S_b5(predd)*pow(LbyAge(predd,yr,sp_age),5));  
          if(LbyAge(predd,yr,sp_age)>80)
            S2Age(predd,yr,sp_age)=S_a(predd)+(S_b(predd)*80)+(S_b2(predd)*pow(80,2))+(S_b3(predd)*pow(80,3))+(S_b4(predd)*pow(80,4))+(S_b5(predd)*pow(80,5)); // set everything above 80 to 80
        }
      }
    } 
    for(int predd=1;predd<=nspp;predd++)
      Limit(predd)=aLim(predd)+bLim(predd)*lengths(predd);    // upper limit of prey size that each pred of each size can consume
    for(int predd=1;predd<=nspp;predd++)
      for(int yr=1;yr<=nyrs;yr++)
        LimitAge(predd,yr)=aLim(predd)+bLim(predd)*LbyAge(predd,yr);    // upper limit of prey size that each pred of each size can consume  
    for(int pred=1;pred<=nspp;pred++)
    {
      for(int prey=1;prey<=nspp;prey++)
      {
        for(int pred_age=1;pred_age<=nages(pred);pred_age++)
        {
          dvector tmp(1,nlengths(prey));
          for(int preyL=1;preyL<=nlengths(prey);preyL++)
          {
            tmp(preyL)=aLim(pred)+bLim(pred)*sum(elem_prod(lengths(pred),A2L_matrix(pred,pred_age)))-lengths(prey,preyL);
            if(tmp(preyL)<0)
              tmp(preyL)=0.;
            if(tmp(preyL)>1)
              tmp(preyL)=1.;
          }// end preyL
          //tmp(preyL)=min(1,max(0,aLim(pred)+bLim(pred)*sum(elem_prod(lengths(pred),A2L_matrix(pred,pred_age)))-lengths(prey,preyL) ));
          for(int prey_age=1;prey_age<=nages(prey);prey_age++)
            LimitAge2(pred,prey,pred_age,prey_age)=sum(elem_prod(tmp,A2L_matrix(prey,prey_age)));
        }// end pred_age
        for(int yr=1;yr<=nyrs;yr++)
        {
          for(int pred_age=1;pred_age<=nages(pred);pred_age++)
          {
            dvector tmp(1,nlengths(prey));
            for(int preyL=1;preyL<=nlengths(prey);preyL++)
            {
              //          tmp(preyL)=min(1,max(0,LimitAge(pred,yr,pred_age)-lengths(prey,preyL) ));
              tmp(preyL)=LimitAge(pred,yr,pred_age)-lengths(prey,preyL);
              if(tmp(preyL)<0)
                tmp(preyL)=0.;
              if(tmp(preyL)>1)
                tmp(preyL)=1.;
            }//end preyL

          for(int prey_age=1;prey_age<=nages(prey);prey_age++)
            swtch(yr,pred,prey,pred_age,prey_age)=sum(elem_prod(tmp,A2L_matrix(prey,prey_age)));
          }// end pred_age
        }// end yr
      } // end prey
    }// end pred  
  // index Tyrs against years

  TempC=sum(BTempC_retro)/nTyrs; 
  int yr_ind;
  for(int yr=1;yr<=nyrs;yr++)
    years(yr)=styr+yr-1; 
  for(int yr=1;yr<=nTyrs;yr++)
  {
    yr_ind = Tyrs(yr)-styr+1;
    if(yr_ind>0)
      if(yr_ind<=nyrs)
        TempC(yr_ind)=BTempC_retro(yr);
  }  
  fT.initialize();
  for(int yr=1;yr<=nTyrs;yr++)
  { 
    for(int predd=1;predd<=nspp;predd++)
    {// for each pred predd
      if (Ceq(predd)==1)
        fT(predd,yr)=mfexp(Qc(predd)*TempC(yr));
      if (Ceq(predd)==2)
      {
        double Yc;
        double Zc;  
        double Vc;
        double Xc;
        Yc=log(Qc(predd))*(Tcm(predd)-Tco(predd)+2);
        Zc=log(Qc(predd))*(Tcm(predd)-Tco(predd));
        Vc=(Tcm(predd)-TempC(yr))/(Tcm(predd)-Tco(predd));
        Xc=(pow(Zc,2))*(pow((1+(pow((1+40/Yc),0.5))),2))/400; 
        fT(predd,yr)=(pow(Vc,Xc))*mfexp(Xc*(1-Vc));    
      }
      if (Ceq(predd)==3)
      {
        double G2;
        double L2;  
        double G1;
        double L1;
        double Ka;
        double Kb;        
        G2=(1/(Tcl(predd)-Tcm(predd)))*log((0.98*(1-CK4(predd)))/(CK4(predd)*0.02));
        L2=mfexp(G2*(Tcl(predd)-TempC(yr)));
        Kb=(CK4(predd)*L2)/(1+CK4(predd)*(L2-1));
        G1=(1/(Tco(predd)-Qc(predd)))*log((0.98*(1-CK1(predd)))/(CK1(predd)*0.02));
        L1=mfexp(G1*(TempC(yr)-Qc(predd)));
        Ka=(CK1(predd)*L1)/(1+CK1(predd)*(L1-1));
        fT(predd,yr)=Ka*Kb;
      }
    }
  }
  fT_fut.initialize();
  for(itemp=1;itemp<=ntemp_scen;itemp++)
  { 
    TempC_fut.initialize();
    int simul;
    simul=simset(itemp);
    TempC_fut=BTempC_fut_all(simul);
    TempC_futUSE(itemp)=BTempC_fut_all(simul);
    for(int yr=1;yr<=nyrs_fut;yr++)
    { 
      for(int predd=1;predd<=nspp;predd++)
      {// for each pred predd
        if (Ceq(predd)==1)
          fT_fut(itemp,predd,yr)=mfexp(Qc(predd)*TempC_fut(yr));
        if (Ceq(predd)==2)
        {
          double Yc;
          double Zc;  
          double Vc;
          double Xc;
          Yc=log(Qc(predd))*(Tcm(predd)-Tco(predd)+2);
          Zc=log(Qc(predd))*(Tcm(predd)-Tco(predd));
          Vc=(Tcm(predd)-TempC_fut(yr))/(Tcm(predd)-Tco(predd));
          Xc=(pow(Zc,2))*(pow((1+(pow((1+40/Yc),0.5))),2))/400; 
          fT_fut(itemp,predd,yr)=(pow(Vc,Xc))*mfexp(Xc*(1-Vc));    
        }
        if (Ceq(predd)==3)
        {
          double G2;
          double L2;  
          double G1;
          double L1;
          double Ka;
          double Kb;        
          G2 =(1/(Tcl(predd)-Tcm(predd)))*log((0.98*(1-CK4(predd)))/(CK4(predd)*0.02));
          L2 =mfexp(G2*(Tcl(predd)-TempC_fut(yr)));
          Kb =(CK4(predd)*L2)/(1+CK4(predd)*(L2-1));
          G1 =(1/(Tco(predd)-Qc(predd)))*log((0.98*(1-CK1(predd)))/(CK1(predd)*0.02));
          L1 =mfexp(G1*(TempC_fut(yr)-Qc(predd)));
          Ka =(CK1(predd)*L1)/(1+CK1(predd)*(L1-1));
          fT_fut(itemp,predd,yr)=Ka*Kb;
        }
      }
    }
  }
 END_CALCS  
   //matrix  rmult(1,ntemp_scen,1,nyrs_fut);                    // random multiplier for random recruitment function
   matrix  rmult(1,nspp,1,nyrs_fut);                    // random multiplier for random recruitment function
   matrix  rmultF(1,nspp,1,nyrs_fut);                    // random multiplier for random recruitment function
   vector  NewCatch(1,nspp);                    //
   int    dietphase;
 LOCAL_CALCS
   for(int yr=1;yr<=nyrs_fut;yr++){
    for(int sp=1;sp<=nspp;sp++){
       rmult(sp,yr)=randn(yr+seed(sp));  //random multiplier
       rmultF(sp,yr)=randn(yr+3+seed(sp+3));  // random F multiplier
    }
   }  
   if(dump_rep) cout<<"_________________________ "<<endl;
   if(dump_rep) cout<<"_________________________ "<<endl;
   if(dump_rep) cout<<"_________________________ "<<endl;
   if(dump_rep) cout<<"seed "<<seed<<endl;
   //cout<<"rmult "<<rmult<<endl;
   if(dump_rep) cout<<"_________________________ "<<endl;
   if(dump_rep) cout<<"_________________________ "<<endl;
   if(dump_rep) cout<<"_________________________ "<<endl;
 END_CALCS 
  // !!                  cout<<"rmult"<<rmult<<endl;

  int    rec_phaseBH;
  int    m2phasenum;
  int    eitphase;
  // int    futphase;
  // int    fNphase1;
  // int    fNphase1_2;
  // int    fNphase2;
  // int    fNphase2_2;
  // int    fNphase3;
  // int    iscenstart;
  // int    agespec_m;
  // int    agespecT;
  
 LOCAL_CALCS
  // iscenstart =   1;  // projection scenarios begin with 1 then switch to 2:2 once scen 1 is complete
  // agespec_m  =  -1;
  // agespecT   =  -1;
  m2phasenum =  -7;
  // fNphase1   =  -7;
  // fNphase1_2 =  -7;
  // fNphase2   =  -7;
  // fNphase2_2 =  -7;
  // fNphase3   =  -7;
  if(recMode==2)
    rec_phaseBH=-7;
  else
    rec_phaseBH=-7; 
    
  eitphase=maxphase+max(max(logist_sel_phase),max(sel_phase))+1;

  // if (do_fut)
  //   futphase=eitphase+2;
 END_CALCS  


// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Set up global Params
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 


// ------------------------------------------------------------------------- 
// -------------------------------------------------------------------------
// Do some initial calcs
// ------------------------------------------------------------------------- 
// ------------------------------------------------------------------------- 

//########################################################################################################################
//########################################################################################################################

INITIALIZATION_SECTION
  srv_sel_inf 4.5
  log_eit_q -6.7025
  srv_sel_slp .9
  ln_mn_rec 9.
  log_srv_q 0.;
  ln_mean_F -.8
  ln_mean_M2 -.8
  h_rec .7


//==============================================================================
// PARAMETER SECTION
//==============================================================================
PARAMETER_SECTION
  !!  if(dump_rep) cout<<"============= PARAMETER SECTION ============="<<endl;  
  // Estimated parameters
  init_vector               ln_mn_rec(1,nspp,1)                             // mean recruitment parameter
  init_vector               ln_mean_F(1,nspp,1)                             // log of mean F    parameter (2 here means estimate in phase 2)  
  init_matrix               ln_mean_M2(1,nspp,1,nages,m2phasenum)           // log of mean F    parameter (2 here means estimate in phase 2)  
  init_bounded_matrix       rec_dev(1,nspp,1,nyrs,-10,10,2)                 // recruitment deviations rec=ln_mn_rec + rec_dev  
  init_bounded_matrix       init_dev(1,nspp,2,nages,-10,10,2)               // initial sp_age structure deviations Na,1979 = ln_mn_rec + init_dev
  init_bounded_matrix       F_dev(1,nspp,1,nyrs,-10,10,2)                   // F_dev   Fa,y = mfexp(ln_mean_F)+F_dev  
  init_bounded_matrix       M1_dev(1,nspp,1,nages,-5,5,-2)                  // different M1 rates for each sp_age and spp. - other mort rates..., if set to-7 don't estimate, just use base M1 rat
  init_bounded_matrix       M2_dev(1,nspp,1,nyrs,-10,10,m2phasenum)
  init_matrix               fsh_sel_coff(1,nspp,1,nselages,2)                        // selectivity parameters
  //  init_bounded_matrix     fsh_sel_coff(1,nspp,1,nselages,0.0001,4,2)    // selectivity parameters
  init_vector_vector        srv_sel_coff(1,nspp,1,nselages,sel_phase)       // survey selectivity parameters
  init_number_vector        srv_sel_inf(1,nspp,logist_sel_phase);           //
  init_number_vector        srv_sel_slp(1,nspp,logist_sel_phase);           //
  init_number               log_eit_q(eitphase)                             // q for eit, phase 4
  init_vector               log_srv_q(1,nspp,-4);                           // q for survey  
  init_matrix               ln_mean_M2_fut(1,nspp,1,nages,-7)               // log of mean F    parameter (2 here means estimate in phase 2)  
  init_vector               ln_sigma_rec(1,nspp,-7) ;                      // mean recruitment parameter 
  init_vector               h_rec(1,nspp,-7) ;

  // init_matrix        ln_fN_mn_alpha(1,nspp,1,nages,fNphase3+1)                // mean attack rate for functional response   of prey consumed
 //  init_vector        ln_fN_mn_alpha(1,nspp,eitphase-2)                // mean attack rate for functional response   of prey consumed
 //  init_matrix        ln_fN_T(1,nspp,1,nspp,eitphase-2)                  // mean handling time for functional response of prey consumed
  // init_matrix        ln_fN_mn_T(1,nspp,1,nages,fNphase2+1)                  // mean handling time for functional response of prey consumed
  //init_bounded_vector    ln_fN_pref(1,nspp,-10,10,eitphase-1)            // devation in attack rate for prey
 //  init_matrix        ln_fN_prefb(1,nspp,1,nages,eitphase-1)            // devation in attack rate for prey (or pref by pred)
 //  init_vector        ln_fN_prefa(1,nspp,eitphase-1)            // devation in attack rate for prey (or pref by pred)
  // init_vector        ln_fN_all_T(1,nspp,fNphase2_2+1)
  // init_matrix        ln_fN_dev_T(1,nspp,1,nages,fNphase2+1)              // mean handling time for functional response of prey consumed
  // init_matrix        ln_fN_dev_alpha(1,nspp,1,nages,fNphase3+1)              // mean handling time for functional response of prey consumed
  // init_matrix        ln_fN_alphaT_other(1,nspp,1,nages,fNphase3+1)          // mean attack rate for functional response   of other prey consumed
  // init_bounded_matrix    ln_fN_m(1,nspp,1,nspp,0,.7,fNphase1_2+1)              // mean handling time for functional response of prey consumed
  // init_bounded_matrix    ln_fN_mn_m(1,nspp,1,nages,0,.3,fNphase1+1)                  // mean handling time for functional response of prey consumed
  // init_bounded_matrix    ln_fN_dev_m(1,nspp,1,nages,0,.4,fNphase1+1)                  // mean handling time for functional response of prey consumed
  // init_vector        oprVal(1,nspp,-1+2*do_fut)                  // parameter used when running the model in operational mode
 
 // Derived parameters
 //==============================================================================
  vector    rec_dev_prior(1,nspp)
  vector    init_dev_prior(1,nspp)
  vector    F_dev_prior(1,nspp)
  vector    PvalueAdj(1,nspp)
  !! PvalueAdj =1;
  matrix    R(1,nspp,1,nyrs)
  3darray   N(1,nspp,1,nyrs+1,1,nages)
  3darray   avail_food(1,nspp,1,nyrs,1,nages)                         // available food 
  3darray   AvgN(1,nspp,1,nyrs,1,nages)
  3darray   F(1,nspp,1,nyrs,1,nages)
  3darray   Z(1,nspp,1,nyrs,1,nages)
  3darray   S(1,nspp,1,nyrs,1,nages) 
  matrix    M1(1,nspp,1,nages)
  3darray   M2(1,nspp,1,nyrs,1,nages)
  3darray   B_eaten(1,nspp,1,nyrs,1,nages)
   // 3darray   of_stom(1,2,1,nspp,1,21)                              // year, pred, pred sp_age? other food stomach content, change nages to 21
  3darray   of_stomKir(1,nspp,1,maxA,1,nyrs)                       // year, pred, pred sp_age? other food stomach content, change nages to 21
  // 3darray   of_stom2(1,nyrs,1,nspp,1,21)                          // year, pred, pred sp_age? other food stomach content, change nages to 21
  matrix    stomtmp(1,nspp,1,nages)                                // year, pred, pred sp_age? other food stomach content, change nages to 21
  3darray   overlap(1,nyrs,1,nspp,1,nspp)
  // 5darray   stom_div_bio(1,nstom,1,nspp,1,21,1,nspp,1,21)
  5darray   stom_div_bio2(1,nyrs,1,nspp,1,21,1,nspp,1,21)
  // 5darray   suit_std(1,nstom,1,nspp,1,21,1,nspp,1,21)
  4darray   suit_main(1,nspp,1,21,1,nspp,1,21)
  4darray   suit_mainKir(1,nspp,1,maxA,1,nspp,1,maxA)
  matrix    suit_other(1,nspp,1,nages)
  matrix    suit_other2(1,nspp,1,nages)
  vector    avgsel_srv(1,nspp)
  vector    avgsel_fsh(1,nspp)
  matrix    tc_biom_hat(1,nspp,1,nyrs)
  vector    mean_F(1,nspp)
  matrix    fsh_sel(1,nspp,1,nages)                                     // selectivity
  3darray   fsh_age_hat(1,nspp,1,nyrs,1,fsh_age_bins)
  matrix    tc_hat(1,nspp,1,nyrs)                  // Total (all ages) catch
  3darray   catch_hat(1,nspp,1,nyrs,1,nages)             // estimated catch
  number    eit_q;
  vector    srv_q(1,nspp);
  matrix    srv_biom_hat(1,nspp,1,nyrs_srv_biom)
  matrix    srv_sel(1,nspp,1,nages)
  3darray   srv_age_hat(1,nspp,1,nyrs_srv_age,1,srv_age_bins)         // estimated sp_age/size compositions
  vector    eit_hat(1,n_eit)
  matrix    eit_age_hat(1,n_eit,1,12)                                 // changed from nages to 12
  
  //Other Parms
  //==============================================================================
  vector    depletion(1,nspp)
  vector    M1_std(1,12)                                       // changed from nages to 12, could be increased as needed
  matrix    prey_consumed(1,nspp,1,nyrs)
  
  //Kir's modification parms
  //==============================================================================  
 
 //By Length
  3darray    L2A_convert(1,nspp,1,nages,1,nlengths)      // conversion matrix for going from lengths to ages
  3darray    AvgN2(1,nspp,1,nyrs,1,nlengths)             // number of obs in each  length matrix for each spp length classes vector
  3darray    N2(1,nspp,1,nyrs+1,1,nlengths)              // number of obs in each  length matrix for each spp length classes vector
  3darray    Consum_living(1,nspp,1,nyrs,1,nlengths)   // pre-allocate, indiviudal consumption in grams per predator
  3darray    ration2(1,nspp,1,nyrs,1,nlengths)
  3darray    pred_E(1,nspp,1,nyrs,1,nspp)           // Total prey of eaten by each pred in each year - grams
  3darray    pred_B(1,nspp,1,nyrs,1,nspp)           // Total prey  eaten by each pred in each year - grams
  3darray    Consum(1,nspp,1,nyrs,1,nlengths)         // pre-allocate, indiviudal consumption in grams per predator
  3darray    E(1,nspp,1,nyrs,1,prey_nlengths)         // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    B4E(1,nspp,1,nyrs,1,prey_nlengths)       // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    availB(1,nspp,1,nyrs,1,prey_nlengths)    // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    Ucheck(1,nspp,1,nyrs,1,prey_nlengths)    // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1 
  3darray    U(1,nyrs,1,totalpreyL,1,totalpredL)        // U is a mega matrix with c(1:nlengths(1),1:nlegnths(2)...1:nlengths(nprey)) as rows, and c(1:nlengths(1),1:nlegnths(2)...1:nlengths(npred)) for column
  4darray    Unew(1,nspp,1,nspp,1,maxL,1,maxL)
  4darray    UnewAVG(1,nspp,1,nspp,1,maxL,1,maxL)
  4darray    Snew(1,nspp,1,nspp,1,maxL,1,maxL)
 
 //By sp_age
  3darray    Consum_livingAge(1,nspp,1,nyrs,1,nages)    // pre-allocate, indiviudal consumption in grams per predator
  3darray    ration2Age(1,nspp,1,nyrs,1,nages)
  3darray    pred_EbyAge(1,nspp,1,nyrs,1,nspp)      // Total prey of eaten by each pred in each year - grams
  3darray    pred_BbyAge(1,nspp,1,nyrs,1,nspp)      // Total prey of eaten by each pred in each year - grams
  3darray    ConsumAge(1,nspp,1,nyrs,1,nages)        // pre-allocate, indiviudal consumption in grams per predator
  3darray    EbyAge(1,nspp,1,nyrs,1,nages)         // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    B4EbyAge(1,nspp,1,nyrs,1,nages)         // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    availBbyAge(1,nspp,1,nyrs,1,nages)       // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    UcheckbyAge(1,nspp,1,nyrs,1,nages)       // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    Eage(1,nspp,1,nyrs,1,nages)          // pre-allocate the biomass of each prey sp_age eaten by all predators
  3darray    Eage_hat(1,nspp,1,nyrs,1,nages)        // pre-allocate the biomass of each prey sp_age eaten by all predators
  3darray    Bage(1,nspp,1,nyrs,1,nages)
  matrix     Bage_other(1,nspp,1,nages)
  3darray    EageN(1,nspp,1,nyrs,1,nages)          // pre-allocate the biomass of each prey sp_age eaten by all predators
  3darray    BageN(1,nspp,1,nyrs,1,nages)
  3darray    Nage(1,nspp,1,nyrs,1,nages) 
  4darray    UnewAge(1,nspp,1,nspp,1,maxA,1,maxA)
  4darray    UnewAVGAge(1,nspp,1,nspp,1,maxA,1,maxA)
  4darray    SnewAge(1,nspp,1,nspp,1,maxA,1,maxA)
  
  matrix     biomassSSB(1,nspp,1,nyrs); // spawning biomass
  3darray    biomassByage(1,nspp,1,nyrs,1,nages);
  matrix     biomass(1,nspp,1,nyrs);
 // 3darray    NByage(1,nspp,1,nyrs,1,nages);  // replace this throughout with N

  // stomach weight parms
  //==============================================================================  

  // 3darray    ration_hat(1,nspp,1,nyrs,1,nages)              // estimated ration
  // 3darray    NEaten(1,nspp,1,nyrs,1,nages)                // number eaten of each prey, prey sp_age across all predators
  // 3darray    NEaten_hat(1,nspp,1,nyrs,1,nages)                // number eaten of each prey, prey sp_age across all predators

  // 4darray    Qhat(1,npred3,1,nyrs+1,1,maxA,1,nspp)            // meanwt of each prey spp in the stomach of each predator of sp_age a
  // 4darray    mn_Uhat(1,nspp,1,nspp,1,maxA,1,maxA)            // V(yr,pred,prey,pred_age,prey_age): meanwt of each prey spp in the stomach of each predator of sp_age a
  // 5darray    Uhat(1,nyrs,1,nspp,1,nspp,1,maxA,1,maxA)          // V(yr,pred,prey,pred_age,prey_age): meanwt of each prey spp in the stomach of each predator of sp_age a
  
  // 5darray    V(1,nyrs,1,nspp,1,nspp,1,maxA,1,maxA)            // V(yr,pred,prey,pred_age,prey_age): meanwt of each prey spp in the stomach of each predator of sp_age a
  //  matrix    fN_mn_alpha(1,nspp,1,nages)                // mean attack rate for functional response   of prey consumed (pred,pred_age)
  //vector    fN_mn_alpha(1,nspp)                      // mean attack rate for functional response   of prey consumed (pred,pred_age)
  //  matrix    fN_mn_T(1,nspp,1,nages)                  // mean handling time for functional response of prey consumed
  //  3darray    fN_pref(1,nspp,1,nspp,1,nages)              // devation in attack rate for prey (prey,prey_age)
  //  vector    fN_pref(1,nspp)                      // devation in attack rate for prey (prey,prey_age)
  //  matrix    fN_dev_T(1,nspp,1,nages)                // mean handling time for functional response of prey consumed
  // matrix    fN_alphaT_other(1,nspp,1,nages)  
  // matrix    qtmp(1,nspp,1,nages)                    // temp file used to calculate q (pred,pred_age)
  // 3darray    other_prey_eaten(1,nspp,1,nyrs,1,nages)            // other prey eaten used in the denominator of the fun response
  
  // 4darray    fN_alpha(1,nspp,1,nspp,1,maxA,1,maxA)            // attack rate for functional response   of prey consumed
  //4darray    fN_T(1,nspp,1,nspp,1,maxA,1,maxA)              // handling time for functional response of prey consumed
  // 4darray    fN_T(1,nspp,1,nspp,1,maxA,1,maxA)              // handling time for functional response of prey consumed
  // 4darray    fN_m(1,nspp,1,nspp,1,maxA,1,maxA)              // handling time for functional response of prey consumed
  // 4darray    KAge_hat(1,nyrs,1,nspp,1,nspp,1,maxA)             // K(1) is npreyXmaxL matrix of stomach proportions of predator 1 
  // 3darray    mnKAge_hat(1,nspp,1,nspp,1,nages)             // K(1) is npreyXmaxL matrix of stomach proportions of predator 1 

  // liklihood parms
  //==============================================================================  

  objective_function_value     obj_fun;
  number    eit_srv_like;  
  vector    srv_sel_like(1,nspp);
  vector    fsh_sel_like(1,nspp);
  vector    norm_srv_like(1,nspp);
  vector    norm_fsh_like(1,nspp);
  vector    srv_bio_like(1,nspp);
  vector    fsh_cat_like(1,nspp);
  vector    srv_age_like(1,nspp);
  vector    fsh_age_like(1,nspp);
  vector    M2_pen(1,nspp);
  number    eit_age_like
  vector    M1_like(1,nspp); 
  vector    M2_like(1,nspp);
  // vector    Q_like(1,nspp);
  // vector    U_like(1,nspp);
  // vector    K_like(1,nspp);
  // vector    ration_like(1,nspp);
  
 //Temp parms  
 //==============================================================================  
  //matrix        Stmp(1,nspp,1,nages);
  //matrix        Ftmp(1,nspp,1,nages);
 // matrix        Ztmp(1,nspp,1,nages);
  //vector        Ctmp(1,nspp);
  //vector        Btmp(1,nspp);  
  //vector        BtmpSSB(1,nspp);  // spawning biomass
  //==============================================================================  
  //==============================================================================  
  // Global parameters
  //==============================================================================  
  //==============================================================================  
  //Future parameters
  //==============================================================================  
  vector     BtargetRep(1,nspp);
  matrix     F_dev_fut(1,nspp,1,nyrs_fut)                 // F_dev   Fa,y = mfexp(ln_mean_F)+F_dev  
  3darray    wt_fut(1,nspp,1,nyrs_fut,1,nages)
  4darray    wt_fut_all(1,ntemp_scen,1,nspp,1,nyrs_fut,1,maxA)
  3darray    R_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)
  3darray    RationScaledKeep_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)
  3darray    RationKeep_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)
  4darray    ration2Age_futKEEP(1,ntemp_scen,1,nspp,1,nyrs_fut,1,maxA)
  3darray    AvgN_futKEEP(1,ntemp_scen,1,nspp,1,nyrs_fut)

  3darray    SSB_RSkeep_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)
  3darray    N_fut(1,nspp,1,nyrs_fut+1,1,nages)
  3darray    av_food_fut(1,nspp,1,nyrs_fut,1,nages)                         // available food 
  3darray    AvgN_fut(1,nspp,1,nyrs_fut,1,nages)
  3darray    F_fut(1,nspp,1,nyrs_fut,1,nages)
  3darray    Z_fut(1,nspp,1,nyrs_fut,1,nages)
  3darray    S_fut(1,nspp,1,nyrs_fut,1,nages) 
  3darray    M2_fut(1,nspp,1,nyrs_fut,1,nages)  
  3darray    B_eaten_fut(1,nspp,1,nyrs_fut,1,nages)  
  4darray    M2_fut_all(1,ntemp_scen,1,nspp,1,nyrs_fut,1,maxA)
  4darray    F_fut_all(1,ntemp_scen,1,nspp,1,nyrs_fut,1,maxA)
  3darray    fsh_age_hat_fut(1,nspp,1,nyrs_fut,1,nages)

  3darray    tc_biom_hat_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)
  3darray    tc_hat_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)                  // Total (all ages) catch
  4darray    catch_hat_fut(1,ntemp_scen,1,nspp,1,nyrs_fut,1,maxA)             // estimated catch
  3darray    Frate_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)

  3darray    ABC_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)
  3darray    ABC_N_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)                  // Total (all ages) catch
  4darray    ABC_N_age_fut(1,ntemp_scen,1,nspp,1,nyrs_fut,1,maxA)             // estimated catch


  // 3darray    overlap_fut(1,nyrs_fut,1,nspp,1,nspp)
  // matrix    srv_biom_hat(1,nspp,1,nyrs_srv_biom)
  3darray    biomassSSB_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut); // spawning biomass
  3darray    biomassSSB0_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut); // spawning biomass
  3darray    R0_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut); // spawning biomass
  3darray    biomassSSB0_hind2fut_CR(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut); // spawning biomass for control rules
  4darray    biomassByage_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut,1,maxA);
  4darray    biomassByage0_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut,1,maxA);
  3darray    biomass_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut);
  3darray    biomass0_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut);
  4darray    NByage_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut,1,maxA);
  4darray    NByage0_hind2fut(1,ntemp_scen,1,nspp,1,nyrs+nyrs_fut,1,maxA);  // consollodate with N_fut_out

  matrix     depletion_fut(1,ntemp_scen,1,nspp)
  matrix     depletion0_fut(1,ntemp_scen,1,nspp)
  3darray    prey_consumed_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)
  3darray    prey_consumed0_fut(1,ntemp_scen,1,nspp,1,nyrs_fut)

  // FUTURE parms
  vector      rec_sigma(1,nspp);                    // random multiplier for random recruitment functionrec_
  matrix      C_fut(1,nspp,1,nyrs_fut);
  //matrix      tmpZ(1,nspp,1,nages);
  //number    tmpF;
  // number     Chat
  3darray    PF(1,ntemp_scen,1,nspp,1,nyrs_fut);                           // Future F rate that is used to calculate catch or is based on catch
  3darray    Fabc(1,ntemp_scen,1,nspp,1,nyrs_fut)                          // Fabc from control rule
  3darray    Fofl(1,ntemp_scen,1,nspp,1,nyrs_fut)                          // Fofl from control rule
  
  4darray    N_fut_out(1,ntemp_scen,1,nspp,1,nyrs_fut+1,1,maxA)                  //

 //By Length
  3darray    AvgN2_fut(1,nspp,1,nyrs_fut,1,nlengths)              // number of obs in each  length matrix for each spp length classes vector
  3darray    N2_fut(1,nspp,1,nyrs_fut+1,1,nlengths)              // number of obs in each  length matrix for each spp length classes vector
  3darray    Consum_living_fut(1,nspp,1,nyrs_fut,1,nlengths)    // pre-allocate, indiviudal consumption in grams per predator
  3darray    ration2_fut(1,nspp,1,nyrs_fut,1,nlengths)
  3darray    pred_E_fut(1,nspp,1,nyrs_fut,1,nspp)        // Total prey of eaten by each pred in each year - grams
  3darray    pred_B_fut(1,nspp,1,nyrs_fut,1,nspp)        // Total prey of eaten by each pred in each year - grams
  3darray    Consum_fut(1,nspp,1,nyrs_fut,1,nlengths)        // pre-allocate, indiviudal consumption in grams per predator
  3darray    E_fut(1,nspp,1,nyrs,1,prey_nlengths)         // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    B4E_fut(1,nspp,1,nyrs,1,prey_nlengths)       // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    availB_fut(1,nspp,1,nyrs_fut,1,prey_nlengths)     // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    Ucheck_fut(1,nspp,1,nyrs_fut,1,prey_nlengths)     // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1 
  //4darray    UnewAVG_fut(1,nspp,1,nspp,1,maxL,1,maxL)
 
 //By sp_age
  3darray    Consum_livingAge_fut(1,nspp,1,nyrs_fut,1,nages)    // pre-allocate, indiviudal consumption in grams per predator
  3darray    ration2Age_fut(1,nspp,1,nyrs_fut,1,nages)
  3darray    pred_EbyAge_fut(1,nspp,1,nyrs_fut,1,nspp)      // Total prey of eaten by each pred in each year - grams
  3darray    pred_BbyAge_fut(1,nspp,1,nyrs_fut,1,nspp)      // Total prey of eaten by each pred in each year - grams
  3darray    ConsumAge_fut(1,nspp,1,nyrs_fut,1,nages)        // pre-allocate, indiviudal consumption in grams per predator
  3darray    EbyAge_fut(1,nspp,1,nyrs_fut,1,nages)         // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    B4EbyAge_fut(1,nspp,1,nyrs_fut,1,nages)         // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    availBbyAge_fut(1,nspp,1,nyrs_fut,1,nages)       // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    UcheckbyAge_fut(1,nspp,1,nyrs_fut,1,nages)       // pre-allocate total grams of a prey eaten by all predators in a given year y, E(1) is a nrys x maxL matrix of grams eaten of prey 1
  3darray    Eage_fut(1,nspp,1,nyrs_fut,1,nages)          // pre-allocate the biomass of each prey sp_age eaten by all predators
  3darray    Eage_hat_fut(1,nspp,1,nyrs_fut,1,nages)        // pre-allocate the biomass of each prey sp_age eaten by all predators
  3darray    Bage_fut(1,nspp,1,nyrs_fut,1,nages)
  3darray    EageN_fut(1,nspp,1,nyrs_fut,1,nages)          // pre-allocate the biomass of each prey sp_age eaten by all predators
  3darray    BageN_fut(1,nspp,1,nyrs_fut,1,nages)
  3darray    Nage_fut(1,nspp,1,nyrs_fut,1,nages) 

  //==============================================================================
// recruitment parms

  matrix      logR_hat(1,nspp,1,nyrs);
  //vector      sigma(1,nspp);
  vector      f_rs(1,nspp);  
  // matrix      logR_obs(1,nspp,1,nyrs);
  vector      nobs(1,nspp);
  
  !!  if(dump_rep) cout<<"============= END: PARAMETER SECTION ============="<<endl;  

   //==============================================================================
   //TRASH Bin
   //==============================================================================
   //4darray    UnewAVGAge_fut(1,nspp,1,nspp,1,maxA,1,maxA)
//==============================================================================
// END OF PARAMETER SECTION
//==============================================================================  
PRELIMINARY_CALCS_SECTION
  if(dump_rep) cout<<"============= PRELIMINARY CALCS SECTION ============="<<endl;  
  cout<<"...";
  CALC_WT_AT_AGE_FUT();
    dietphase=50;  // do not change this or the hindcast iterations will not work
  if(do_fut||simMode){
    dietphase=6;  // do not change this or the hindcast iterations will not work
    nobs=nyrs;
    // for (int sp=1;sp<=nspp;sp++)
      // logR_obs(sp) = (ln_mn_rec_rs(sp) + rec_dev_rs(sp)(1,nyrs) );
  }
  // for (int yr=1;yr<=2;yr++) // year
  //  for (int pred=1;pred<=nspp;pred++) // pred sp
  //   for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age sp
  //    of_stom(yr,pred,pred_age) =(1-sum(stom(yr,pred,pred_age)))/other_food(1); 
  for (int yr=1;yr<=nyrs;yr++) // year
   for (int pred=1;pred<=nspp;pred++) // pred sp
    for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age sp
      for (int prey=1;prey<=nspp;prey++) // pred sp_age sp
        stomKir(yr,pred,pred_age,prey)=UobsAge(pred,prey,pred_age); 
 
 for (int yr=1;yr<=nyrs;yr++) // year
   for (int pred=1;pred<=nspp;pred++) // pred sp
    for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age sp
     of_stomKir(pred,pred_age,yr) =(1-sum(stomKir(yr,pred,pred_age)))/other_food(1); 
  for (int pred=1;pred<=nspp;pred++)
    for (int pred_age=1;pred_age<=nages(pred);pred_age++)
       for (int prey=1;prey<=nspp;prey++) // pred sp_age sp
        for (int prey_age=1;prey_age<=nages(prey);prey_age++)
         stomtmp(pred,pred_age)+=UobsAge(pred,prey,pred_age,prey_age);
  suit_other2=1-stomtmp;
  for (int predd=1;predd<=nspp;predd++) //predator species
    for (int preyy=1;preyy<=nspp;preyy++) //predator species
      for (int pred_age=1;pred_age<=nages(predd);pred_age++) // predator sp 
        other_food_byAge(predd,preyy,pred_age)=(1-KAge(predd,preyy,pred_age))*other_food(predd);
  // to add down the line


  for (int yr=1;yr<=nyrs;yr++) // year
   for (int pred=1;pred<=nspp;pred++) // pred sp
    for (int prey=1;prey<=nspp;prey++) // pred sp_age sp
      overlap(yr,pred,prey)=1.; 
  int tmp_y;
  for (int pred=1;pred<=nspp;pred++){
    for (int prey=1;prey<=nspp;prey++){
      for(int yr=1;yr<=nTyrs;yr++){
        tmp_y = Tyrs(yr)-styr+1;
        if(tmp_y>0)
          if(tmp_y<=nyrs)
          overlap(tmp_y,pred,prey)=overlap_dat(pred,prey,yr);
      }  
    }
  }

  log_input( of_stomKir );
  // cout<<"overlap "<<overlap<<endl;exit(1);
      // overlap(yr,pred,prey)=1.; 


  // for (int yr=1;yr<=nyrs;yr++) // year
  //   overlap(yr,3,1)=PLKoverlap(1,yr);  
  // for (int yr=1;yr<=nyrs_fut;yr++) // year
  //  for (int pred=1;pred<=nspp;pred++) // pred sp
  //   for (int prey=1;prey<=nspp;prey++) // pred sp_age sp
  //     overlap_fut(yr,pred,prey)=sum(PLKoverlap(1))/nyrs;  

  CALC_RATION(0); 
  if (option_match(ad_comm::argc,ad_comm::argv,"-binp")>-1)
  {
      int tmp_dofut=do_fut;
      int tmp_simMode=simMode;
      do_fut=0;
      simMode=0;
      cout<<"...";
      UPDATE_BETWEEN();  
      cout<<"..";
      RUN_ESTIMATION();
      dietphase=0;   // do not change this or the hindcast iterations will not work

      for(int rep=1;rep<=rep_in;rep++)
      {
        cout<<".";
        UPDATE_BETWEEN();  
        RUN_ESTIMATION();
      }
      cout<<".Done"<<endl;
      cout<<""<<endl;
      CALC_PREDICTED_VALS(0); 
      do_fut=tmp_dofut;
      simMode=tmp_simMode;
        if(dump_rep) cout<<"do_fut ="<<do_fut<<"| simMode= "<<simMode<<endl;
  }
  
  if(dump_rep) cout<<"SSB "<< biomassSSB<<endl;
  /*if(do_fut||simMode)
    for (itemp=1;itemp<=ntemp_scen;itemp++)
      for (int i=1;i<=nyrs_fut;i++)
        CALC_REC_FUT(0);  // mean recruitment*/
  dietphase=2;    // do not change this or the hindcast iterations will not work

 if(dump_rep) cout<<"============= END: PRELIMINARY CALCS SECTION ============="<<endl;  
  // if(do_fut)
  //   {
  //     CALC_SUIT(msmMode);
  //     dietphase=0;  // do not change this or the hindcast iterations will not work
  //     if(dump_rep) cout << "In PROCEDURE_SECTION"<<endl;

  //     switch (proj_mode)
  //     {
  //       case 3:
  //         RUN_CONTROL_RULES();
  //       break;
  //       case 9:
  //         RUN_F_PROFILE();
  //       break;
  //       default:
  //         RUN_PROJECTIONS(); 
  //       break;
  //     }
  //   }
PROCEDURE_SECTION
  if(do_fut)
  {
    CALC_SUIT(msmMode);
    dietphase=0;  // do not change this or the hindcast iterations will not work
    if(dump_rep) cout << "IN PROCEDURE_SECTION"<<endl;
    /*if(harvestMode==3)
      RUN_CONTROL_RULES();
    else
      RUN_PROJECTIONS(); */
    switch (harvestMode)
    {
      case 3:
        if(dump_rep) cout << "running control rules; harvestMode = "<<harvestMode<<endl;
        RUN_CONTROL_RULES();
      break;
      case 9:
        if(dump_rep) cout << "running f-profiles; harvestMode = "<<harvestMode<<endl;
        RUN_F_PROFILE();
        exit(1);
      break;
      default:
         if(dump_rep) cout << "running projections under harvestMode = "<<harvestMode<<endl;
        RUN_PROJECTIONS(); 
      break;

    }
  }
  else
  {
    // Note that running hindcast when do_fut=1 results in some errors because suitability is calculated in future mode. Thus it needs to stay inside this if() statement
    if(dump_rep) cout<<"RUN_ESTIMATION begins"<<endl;
    RUN_ESTIMATION();
    if(dump_rep) cout<<"RUN_ESTIMATION ends"<<endl;
    CALC_PREDICTED_VALS(0);    
    if(dump_rep) cout<<"CALC_PREDICTED_VALS ends"<<endl;
  }
  if(dump_rep) cout<<"CALC_OBJ_FUN starts"<<endl;
  CALC_OBJ_FUN();
  if(dump_rep) cout<<"CALC_OBJ_FUN ends"<<endl;

//==============================================================================
// FUNCTIONS
//==============================================================================
FUNCTION UPDATE_BETWEEN
  if(do_fut==0)
  {
    if(msmMode==44)
    {
      for(int predd=1;predd<=nspp;predd++)
        PvalueAdj(predd)=current_phase()/(maxphase-10);
      if((maxphase-10)<=0)
      {
        cout<<"ERROR!!!!! BETWEEN_PHASES_SECTION : maxphase must be larger than 10 - please update .dat file"<<endl;
        exit(1);
      }
      else
      {
        PvalueAdj=1;
      }
    }
  } 
  if(current_phase()>dietphase)
  {
    CALC_SUIT(msmMode);
  }else{  
    if(msmMode>1)
    {
      CALC_UandC(0);
      CALC_SUIT(msmMode);
    }else{
      CALC_SUIT(msmMode);
    }  
  }
  if(dump_rep) cout<<"M2(prey=1,yr=10) "<<M2(1,10)<<endl;
  if(dump_rep) cout<<"M2(prey=2,yr=10) "<<M2(2,10)<<endl;
  if(dump_rep) cout<<"M2(prey=3,yr=10) "<<M2(3,10)<<endl;

FUNCTION CALC_OBJ_FUN
  M1_like.initialize();
  M2_like.initialize();
  srv_sel_like.initialize();
  fsh_sel_like.initialize();
  norm_srv_like.initialize();
  norm_fsh_like.initialize();
  srv_bio_like.initialize();
  fsh_cat_like.initialize();
  srv_age_like.initialize();
  fsh_age_like.initialize();
  eit_age_like.initialize();
  int yr_ind;
 
  for (int sp=1;sp<=nspp;sp++)
  {
    if(debugg) cout<<"Start fsh_age_like"<<endl;
    if(debugg) cout<<"Start nyrs_fsh_comp "<<nyrs_fsh_comp<<endl;
    for (int yr=1;yr<=nyrs_fsh_comp(sp);yr++)
    {
      // need to index years here!!!
      yr_ind  = yrs_fsh_comp(sp,yr) - styr + 1;// convert years into indices for 1-=nyrs counting purposes    
      fsh_age_like(sp)  -= offset_fsh(sp) + tau*sum(elem_prod(fsh_age_obs(sp,yr) + MNConst,log(fsh_age_hat(sp,yr_ind) + MNConst )));    // Multinomial
    }
    if(debugg) cout<<"End fsh_age_like"<<endl;
    norm_srv_like(sp) += 50. * square(avgsel_srv(sp));                               // To normalize selectivity parameters
    norm_fsh_like(sp) += 50. * square(avgsel_fsh(sp));                               // To normalize selectivity parameters
    
    fsh_cat_like(sp)  += 200.*norm2(log(tc_biom_obs(sp)+1.e-4)-log(tc_biom_hat(sp)+1.e-4));                     // Errors in total catch estimation
    // too hevy a weight? should we try 12.5 ? (ie 20% CV)?
    // 5% CV
    
    // M2_like(sp)+=norm2(log(prey_consumed(sp)+1.e-4)-log(Eage(sp)+1.e-4));
    // if (!last_phase()) obj_fun+= 10.*square(log(mean(F(sp))/.15));
   
     // ----Selectivity penalties...for smooth second diff etc
     //           Invoke a penalty when the partial F's change abruptly (unidirectional 1 to .5 is the same as .5 to 1, for example)
    for (int sp_age=1;sp_age<nages(sp);sp_age++)
      if (fsh_sel(sp,sp_age)>fsh_sel(sp,sp_age+1))
        fsh_sel_like(sp) += 20.*square(log(fsh_sel(sp,sp_age)/fsh_sel(sp,sp_age+1)));
    fsh_sel_like(sp) += curv_pen_fsh(sp)*norm2(first_difference( first_difference(log(fsh_sel(sp)))));
     
    if (logist_sel_phase(sp) < 0)
      srv_sel_like(sp) += curv_pen_srv(sp)*norm2(first_difference( first_difference(log(srv_sel(sp)))));
    int tmp;
    //if(current_phase()<last_phase())
    //  tmp=20;
    //else

    tmp=1;
    for (int yr=1;yr<=nyrs_srv_biom(sp);yr++)
      srv_bio_like(sp) += tmp*square(log(srv_biom(sp,yr)) - log(srv_biom_hat(sp,yr)) ) /(2.*srv_biom_lse(sp,yr)*srv_biom_lse(sp,yr));
     
    // Multinomial for bottom trawl survey sp_age-compositions
    // need to index these years! yr != yr here
    for (int yr=1;yr<=nyrs_srv_age(sp);yr++)
      srv_age_like(sp) -= srv_age_n(sp,yr)*sum(elem_prod(srv_age_obs(sp,yr)+MNConst,log(srv_age_hat(sp,yr)+MNConst)));
    srv_age_like(sp) -= offset_srv(sp);
    
    // Priors...
    obj_fun +=  1.* norm2(rec_dev(sp)  );
    obj_fun +=  1.* norm2(init_dev(sp) );
    obj_fun +=  1.* norm2(F_dev(sp)    );

    rec_dev_prior(sp)=1.* norm2( rec_dev(sp)  );
    init_dev_prior(sp)=1.* norm2( init_dev(sp)  );
    F_dev_prior(sp)=1.* norm2( F_dev(sp)  );

    if(msmMode==4)
    {
      if(M2mode==1)
      {
        obj_fun +=  1.* norm2(M2_dev(sp));
        for (int yr=1;yr<=nyrs_srv_biom(sp);yr++)
          M2_like(sp)+= norm2(log(Eage(sp,yr)+1.e-06) - log(Eage_hat(sp,yr)+1.e-06))/ (2*(.2*.2));
      }    
    }

    //      M2_like(sp) += square(log(Eage(sp,i)) - log(Eage_hat(sp,i)) ) / (2*(.05*.05));
    // obj_fun += 100.*(square(mean(F_dev(sp)))) + square(mean(rec_dev(sp))) + square(mean(init_dev(sp)));
    if (last_phase())
      M1_like(sp) += norm2(log(M1) - log(M1_base+1.e-06));
      // for (i=1;i<=n_eit;i++) cout<< eit_age_n(i)*sum(elem_prod(obs_eit_age(i),log(eit_age_hat(i)+1.e-4)))<<endl;
  }//for each spp sp
  
  //obj_fun+=sum(M2_pen);
   // assume 20% var (approx) for eit data(1/(2*.2*.2) = 12.5)
   eit_srv_like=12.5 * norm2(log(obs_eit)-log(eit_hat+1.e-04));
  // obj_fun +=eit_srv_like;

  // shouldn't offset_eit be included here?
  for (int yr=1;yr<=n_eit;yr++)
    eit_age_like -= eit_age_n(yr) * sum( elem_prod(obs_eit_age(yr)+MNConst, log(eit_age_hat(yr)+MNConst)));  
  eit_age_like -= offset_eit;
  
  
  obj_fun += 
       sum(srv_age_like)  + 
       sum(fsh_age_like)  + 
       sum(fsh_cat_like)  + 
       sum(srv_bio_like)  + 
       sum(norm_srv_like) + 
       sum(norm_fsh_like) + 
       sum(srv_sel_like)  +
       sum(fsh_sel_like)  +
       eit_srv_like    + 
       eit_age_like       ;  
  if(active(M1_dev)) // if M1_dev is not estimated in the model
    obj_fun +=sum(M1_like)       ;//adding to CALC_OBJ_FUN
  if(msmMode==4)
    if(M2mode==1)
      obj_fun+=sum(M2_like);  

FUNCTION RUN_ESTIMATION
  for (int sp=1;sp<=nspp;sp++){  
    /*if(active(M1_dev(k)))
        M1(k) = elem_prod((M1_base(k)+1.e-04), mfexp(M1_dev(k)));    
      else*/
    M1(sp) = M1_base(sp)+1.e-04;
  }
  mean_F.initialize();
  M2.initialize();
  eit_q.initialize();
  srv_q.initialize();
  eit_q  = mfexp(log_eit_q);
  srv_q  = mfexp(log_srv_q);
  mean_F = mfexp(ln_mean_F);

  if(msmMode==4)
    if(M2mode==1)
      for(int yr=1;yr<=nyrs;yr++)
        for(int preyy=1;preyy<=nspp;preyy++)
          for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
            M2(preyy,yr,sp_age) = mfexp(ln_mean_M2(preyy,sp_age)+M2_dev(preyy,yr));
  CALC_SELECTIVITY();
  CALC_N_AT_AGE(1);        // Calculate AvgN and N based on pass = 1
 // if (msmMode>0)
  for (int iter=1;iter<=niter;iter++)  
    CALC_N_AT_AGE(2);    // loop three times to get M2 close to answer

FUNCTION RUN_PROJECTIONS

  if(dump_rep) cout<<"============= RUN FUTURE PROJECTIONS | simMode = "<<simMode<<"  ============="<<endl;  
  int iter_nn         = 0;
  prey_consumed_fut.initialize();
  // for projmode=4 only
  if (harvestMode==4){
    for (int lsim=1;lsim<=ntemp_scen;lsim++)
      for (int sp=1;sp<=nspp;sp++)
        PF(lsim,sp) = 0.25;
    }
  R_fut.initialize();
  RationScaledKeep_fut.initialize();
  RationKeep_fut.initialize();
  SSB_RSkeep_fut.initialize();
  ration2Age_futKEEP.initialize();
  AvgN_futKEEP.initialize();
  tc_biom_hat_fut.initialize();
  tc_hat_fut.initialize();
  catch_hat_fut.initialize();
  Frate_fut.initialize();


  for (itemp=1;itemp<=ntemp_scen;itemp++)
  {
    if(itemp==1){
      if(harvestMode==11) WRITE_FUT_PREDATION_REPORT(0);
    }
    cout<<"--- temperature scenario = "<<itemp<<" ---"<<endl;
    iter_nn=iter_nn+1;
    M2_fut.initialize();
    wt_fut.initialize();
    av_food_fut.initialize();
    Z_fut.initialize();
    S_fut.initialize();
    F_fut.initialize();
    N_fut.initialize();
    AvgN_fut.initialize();
    fsh_age_hat_fut.initialize();
    //CAL_RS_PARMS();
    if(dump_rep) cout<<"nyrs_fut "<<nyrs_fut<<endl;
    for (int sp=1;sp<=nspp;sp++)  
      for (int yr=1;yr<=nyrs_fut;yr++)
        wt_fut(sp,yr)=wt_fut_all(itemp,sp,yr)(1,nages(sp));
    CALC_RATION(1);   

    // then run projections under harvestMode=1-3 or 4 if simMode;
    int harvestMode_input = harvestMode;
    iter_nn= 0;
    if(dump_rep) cout<<"harvestMode_input "<<harvestMode_input<<endl;
   //if(harvestMode!=9){
    for (i=1;i<=nyrs_fut;i++)
    {

      // __________________________________________________________________________
      // first run unfishbiomass estimates
      // __________________________________________________________________________

      if(debugg) cout<<"year i =  "<<i<<endl;
      harvestMode=1;                                       
      CALC_N_AT_AGE_FUT(1);                  // Calculate AvgN and N based on pass = 1  
          for (int iter=1;iter<=niter;iter++)
            CALC_N_AT_AGE_FUT(2);
      CALC_PREDICTED_VALS(1);
      if(harvestMode==11)        WRITE_FUT_PREDATION_REPORT(1);
    }//end for each year

    M2_fut.initialize();
    av_food_fut.initialize();
    Z_fut.initialize();
    S_fut.initialize();
    F_fut.initialize();
    N_fut.initialize();
    AvgN_fut.initialize();
    fsh_age_hat_fut.initialize();
    
     
    if(harvestMode==4){
      // Reset all three SS datasets from one script;  This runs all the assessments once with no simulated data added
      if(debugg) cout<<"Now run SS models if needed "<<endl;
      #if defined _WIN32 || defined _WIN64
        system(adstring("reset_all_ss.bat ") + adstring(itoa(endyr-1,10))+"_"+adstring(itoa(itemp,10)));     //  system(adstring("reset_all_ss.bat") );          
      #else 
        system(adstring("sh reset_all_ss.sh ") + adstring(itoa(endyr-1,10))+"_"+adstring(itoa(itemp,10)));     // system("sh reset_all_ss.sh");        
      #endif  
      Get_SS_Catch();       // Get next year's "catch" to initialize the loop
      log_input(NewCatch);
    }

    for (i=1;i<=nyrs_fut;i++)
    {
      // __________________________________________________________________________
      // then run projections under harvestMode=1-3 or 4 if simMode;
      // __________________________________________________________________________

      if(simMode>0)
      {
        // if in simulation mode, execute a script to run all three assessments to get NewCatch.dat
        harvestMode=4;
        cout <<"Catch from assessments: "<<NewCatch<<endl;             
        iter_nn=0;
        CALC_N_AT_AGE_FUT(1);                     // Calculate AvgN and N based on pass = 1  
          for (int iter=1;iter<=niter;iter++)
            CALC_N_AT_AGE_FUT(2);
        CALC_PREDICTED_VALS(1);
          // __________________________________________________________________________
            if(simMode==1)
            {
              // Write this years new data (based on previous year's projected catch)
                      
              ATF_SS();
              POLL_SS();
              PCOD_SS();
              cout<<"Running assmnts, temperature scenario, year "<< itemp<<" "<< endyr+i-1<<" "<<i <<endl;// exit(1);
              #if defined _WIN32 || defined _WIN64
                system(adstring("../Scripts/Shell_scripts/runall.bat ") + adstring(itoa(endyr+i-1,10))+"_"+adstring(itoa(itemp,10)));        
              #else 
                system("ls ../Scripts/Shell_scripts/; exit");
                exit(1);
                system(adstring("sh ../Scripts/Shell_scripts/runall.sh ") + adstring(itoa(endyr+i-1,10))+"_"+adstring(itoa(itemp,10)));        
              #endif  
              Get_SS_Catch();         // Get next year's "catch"  from assessments
              log_input(NewCatch);
            }
          // __________________________________________________________________________

        if(simMode==2){
          // run MSM models
          // defunct
          // create NewCatch.dat file
        }
      } 
      else
      {
        // run the projections according to the projection mode
        harvestMode=harvestMode_input; // reset the harvestMode to the input value
        iter_nn=0;
        if(harvestMode==4)
        {
          Get_SS_Catch();
          cout <<"Catch from assessments2: "<<NewCatch<<endl;
        }
        CALC_N_AT_AGE_FUT(1);                     // Calculate AvgN and N based on pass = 1  
        for (int iter=1;iter<=niter;iter++)
            CALC_N_AT_AGE_FUT(2);
        CALC_PREDICTED_VALS(1);
        if(harvestMode<12|harvestMode>13){
          for (k=1;k<=nspp;k++)
          {
            ABC_fut(itemp,k,i)=tc_biom_hat_fut(itemp,k,i);
            ABC_N_fut(itemp,k,i)=tc_hat_fut(itemp,k,i);
            ABC_N_age_fut(itemp,k,i)=catch_hat_fut(itemp,k,i);
          }
        }
      }// end if simMode pass >0

      if(debugg) cout<<"year (second) i =   "<<i<<endl;
     if (harvestMode==11) WRITE_FUT_PREDATION_REPORT(1);
      
    }//end for each year
  //   if(harvestMode>11){
  //     for (i=1;i<=nyrs_fut;i++)
  //     {
  //       for (k=1;k<=nspp;k++)
  //       {
  //         ABC_fut(itemp,k,i)=tc_biom_hat_fut(itemp,k,i);
  //         ABC_N_fut(itemp,k,i)=tc_hat_fut(itemp,k,i);
  //         ABC_N_age_fut(itemp,k,i)=catch_hat_fut(itemp,k,i);
  //       }
  //     }
  //   }
  // }  // end for ntemp_scen
  // if(harvestMode==3)
  // {
  //   //do not exit
  if(debugg) cout<<"Temp scenario  "<<ntemp_scen<<" Complete"<<endl;
  }  // end for ntemp_scen

  if(harvestMode==3)
  {
    //do not exit

  }else{
    if(dump_rep) cout<<"============= END: RUN FUTURE PROJECTIONS  ============="<<endl; 
    if(dump_rep) cout<<"finish up"<<endl;
    //cout<<"R_fut(1) "<< R_fut(1)<<endl;
    //cout<<"SSB "<< biomassSSB(1)<<endl;
    //cout<<"SSB_fut "<< biomassSSB_hind2fut(1,1)<<endl;
    //cout<<"SSB0_fut "<< biomassSSB_hind2fut(1,1)<<endl;
    if(harvestMode!=9) {
      // write_R();
      if(dump_rep) cout<<"============= END: WRITE R    ========================="<<endl;  
      WRITE_FUT_REPORT(0);
      WRITE_FUT_REPORT(1);
      // oko was WRITE_PROJECTION_REPORT();
      WRITE_PROJECTION_REPORT(0);
      WRITE_PROJECTION_REPORT(1);
            write_R_new();
      exit(1);
    }
  }  

FUNCTION Get_SS_Catch
  // Read next year's catches into 
  // ifstream pollcatch("../poll/poll_catch.dat");
  ifstream pollcatch("../../dat_input_files/SingleSppSub/poll/poll_catch.dat");
  
  pollcatch >> NewCatch(1);
  // NewCatch(1) *= 1000.; // units for pollock change
  pollcatch.close();  
  ifstream pcodcatch("../../dat_input_files/SingleSppSub/pcod_catch.dat");
  pcodcatch >> NewCatch(2);
  pcodcatch.close();  
  ifstream atfcatch("../../dat_input_files/SingleSppSub/atf_catch.dat");
  atfcatch  >> NewCatch(3);
  atfcatch.close();  
  cout <<"Send data to single species models "<<endl;

FUNCTION CAL_RS_PARMS
  //logR_hat.initialize();
  //sigma=mfexp(logsigma);
  //aa_c=mfexp(logaa_c);
  // dvariable Eating;
  // dvariable covars;
  // for (int sp=1;sp<=nspp;sp++)
  // {
  //   logR_hat(sp,1)=logR_obs(sp,1);
  //   for(int yr=2;yr<=nyrs;yr++)
  //   {
  //     //Eating=0.;
  //     //Eating = Total_ration(sp,yr)/fall_zoop(yr);
  //     Eating=0.;
  //     covars=0.;
  //     Eating = Total_ration(sp,yr-1)/rs_cov(fallZ_num,yr-1);
  //     if(noRation(sp))
  //         Eating = rs_cov(fallZ_num,yr-1);
  //       //logR_hat(k,yr)= (log(aa_c(k)*SSB_rs(k,yr-1)) -bb_c(k)*SSB_rs(k,yr-1)+spr_Z(k)*spr_zoop(yr) - fall_Z(k)*Eating+TempCoef(k)*WTC(yr));
  //     for (int c=3;c<=ncov;c++)
  //         covars+=rs_parm(c,sp)*rs_cov(c,yr-1);
  //       logR_hat(sp,yr)= (log(aa_c(sp)*SSB_rs(sp,yr-1)) -bb_c(sp)*SSB_rs(sp,yr-1)+rs_parm(1,sp)*rs_cov(1,yr-1) - rs_parm(fallZ_num,sp)*Eating+covars);
  //     //logR_hat(sp,yr)= (log(aa_c(sp)*SSB_rs(sp,yr-1)) -bb_c(sp)*SSB_rs(sp,yr-1)+spr_Z(sp)*spr_zoop(yr) - fall_Z(sp)*Eating);
  //   } 
  // } 
  // for (int sp=1;sp<=nspp;sp++)
  // {
  //    f_rs(sp) = norm2(logR_hat(sp)-logR_obs(sp));
  //    f_rs(sp)= nobs(sp)*log(sigma(sp)) + f_rs(sp)/(2.0*sigma(sp)*sigma(sp));
  // }
  // obj_fun=sum(f_rs);
  //    //f = regression(size,pred_size);
  //    if(sd_phase())
  //       write_output();
FUNCTION void CALC_PREDICTED_VALS_ol(int fut_pass_number)
  if (fut_pass_number==0)
  {
    prey_consumed.initialize();
    if(dump_rep) cout<<" sp yr calculations begin" <<endl;
    for (int sp=1;sp<=nspp;sp++)
    {
      catch_hat(sp) = elem_prod(elem_div(F(sp),Z(sp)) ,elem_prod(1.-mfexp(-Z(sp)) , N(sp)));
      for (int yr=1;yr<=nyrs;yr++)
      {
        tc_hat(sp,yr)       = sum(catch_hat(sp,yr));
        
        if(fsh_age_type(sp)==1) 
          fsh_age_hat(sp,yr)  = catch_hat(sp,yr) / tc_hat(sp,yr);
        else
          fsh_age_hat(sp,yr)  = catch_hat(sp,yr)*age_trans_matrix(sp) / tc_hat(sp,yr);
        
        tc_biom_hat(sp,yr)  = catch_hat(sp,yr) * wt(sp,yr)(1,nages(sp));// matrix multiplication A%*%B
        biomass(sp,yr)      = N(sp,yr) * wt(sp,yr)(1,nages(sp));  // matrix multiplication A%*%B
        biomassByage(sp,yr) = elem_prod(N(sp,yr), wt(sp,yr)(1,nages(sp)));
       // NByage(sp,yr)       = N(sp,yr);
        biomassSSB(sp,yr)   = elem_prod(N(sp,yr),pmature(sp))*wt(sp,yr)(1,nages(sp));
      }
      if(debugg) cout<<" here 1" <<endl;
      dvar_vector Stmp(1,nages(sp));
      dvar_vector tmp_age(1,nages(sp)); 
            if(debugg) cout<<" here 1.1" <<endl;
      //double yr_ind;
      int yr_ind;
      Stmp.initialize();
      tmp_age.initialize();
      if(debugg) cout<<"nyrs_srv_biom("<<sp<<") "<<nyrs_srv_biom(sp)<<endl;
      for (int yr=1;yr<=nyrs_srv_biom(sp);yr++)
      {      
        yr_ind            = yrs_srv_biom(sp,yr) - styr + 1; // convert years into indices for 1-=nyrs counting purposes 
        //Stmp              = (S(sp,yr_ind)+1.)/2.;// note: assumes survey occurs in mid-year (pow(,.5))      
        Stmp              = mfexp(-Z(sp,yr_ind)*.5);                // note: assumes survey occurs in mid-year       
        srv_biom_hat(sp,yr) = elem_prod(srv_q(sp) * srv_sel(sp) , elem_prod( Stmp ,N(sp,yr_ind)) ) * wt(sp,yr_ind)(1,nages(sp));
      }  
            if(debugg) cout<<" here 2" <<endl;
      for (int yr=1;yr<=nyrs_srv_age(sp);yr++)
      {
        yr_ind  = yrs_srv_age(sp,yr) - styr + 1;// convert years into indices for 1-=nyrs counting purposes    
        tmp_age = N(sp,yr_ind);// need to test for type of sp_age data here (if length bins are different
        tmp_age = elem_prod(srv_sel(sp) , tmp_age);
        if ( srv_age_type(sp)==1)
          srv_age_hat(sp,yr) = tmp_age / sum(tmp_age);
        else
        {
          srv_age_hat(sp,yr)  = tmp_age * age_trans_matrix(sp);
          srv_age_hat(sp,yr) /= sum(srv_age_hat(sp,yr));
        }
      }
            if(debugg) cout<<" here 3" <<endl;
      // Compute time series of biomass consumed by different species (vector x vector = scalar)
      if (sd_phase())
      {   
        depletion(sp) = biomass(sp,nyrs)/(1.e-20 + biomass(sp,1)); 
        if (msmMode>0)
          if (sp<=nspp)
            for (int yr=1;yr<=nyrs;yr++)
              prey_consumed(sp,yr) =sum(elem_prod(wt(sp,yr), elem_prod(elem_div(M2(sp,yr),Z(sp,yr)) ,elem_prod(1.-mfexp(-Z(sp,yr)) , N(sp,yr)))));
              // This makes sense since M2 is already summed over sources of predation...
      }//end if sd_phase()
            if(debugg) cout<<" here 4" <<endl;
    }
    // Hard wired for acoustic survey of pollock only...
          if(debugg) cout<<" here 5" <<endl;
    for (int yr=1;yr<=n_eit;yr++)
    {
      //convert years into indices of counting purposes
      double eit_yrs_ind = yrs_eit(yr) -styr +1;
      eit_age_hat(yr)     = elem_prod((N(1,eit_yrs_ind)),eit_sel(eit_yrs_ind)*eit_q); // No adjustment to mid-year
      eit_hat(yr)         = eit_age_hat(yr) *wt(1,eit_yrs_ind)(1,12);
      eit_age_hat(yr)     /= sum(eit_age_hat(yr));
    }
          if(debugg) cout<<" here 6" <<endl;
  }
  if(fut_pass_number==1)
  {
    for (int sp=1;sp<=nspp;sp++)
    {
      if(i==1)
      {
        for (int hyr=1;hyr<=nyrs;hyr++)
        {
          biomass0_hind2fut(itemp,sp,hyr)                  = biomass(sp,hyr);
          biomass_hind2fut(itemp,sp,hyr)                   = biomass(sp,hyr);
          biomassSSB_hind2fut(itemp,sp,hyr)                = biomassSSB(sp,hyr);
          biomassSSB0_hind2fut(itemp,sp,hyr)               = biomassSSB(sp,hyr);
          R0_hind2fut(itemp,sp,hyr)                        = R(sp,hyr);
          biomassByage_hind2fut(itemp,sp,hyr)(1,nages(sp))  = biomassByage(sp,hyr);
          biomassByage0_hind2fut(itemp,sp,hyr)(1,nages(sp)) = biomassByage(sp,hyr);
          NByage_hind2fut(itemp,sp,hyr)(1,nages(sp))        = N(sp,hyr);
          NByage0_hind2fut(itemp,sp,hyr)(1,nages(sp))       = N(sp,hyr);

        }   
      }
      if(harvestMode==1)
      {
        //unfished biomass
        //        catch_hat_fut(itemp,sp,i)(1,nages(sp))                       = elem_prod(elem_div(F_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)));
        //        tc_hat_fut(itemp,sp,i)                                               = sum(catch_hat_fut(itemp,sp,i)(1,nages(sp)));
        //        fsh_age_hat_fut(sp,i)                                        = catch_hat_fut(itemp,sp,i)(1,nages(sp)) / tc_hat_fut(itemp,sp,i);
        //        tc_biom_hat_fut(itemp,sp,i+nyrs)                                    = catch_hat_fut(itemp,sp,i)(1,nages(sp)) * wt(sp,nyrs)(1,nages(sp));// matrix multiplication A%*%B
        biomass0_hind2fut(itemp,sp,i+nyrs)                  = N_fut(sp,i) * wt_fut(sp,i)(1,nages(sp));  // matrix multiplication A%*%B
        biomassSSB0_hind2fut(itemp,sp,i+nyrs)               = elem_prod(N_fut(sp,i),pmature(sp))*wt_fut(sp,i)(1,nages(sp));   // matrix multiplication A%*%B
        R0_hind2fut(itemp,sp,i+nyrs)                        = R_fut(itemp,sp,i);   // matrix multiplication A%*%B
        
        biomassByage0_hind2fut(itemp,sp,i+nyrs)(1,nages(sp)) = elem_prod(N_fut(sp,i), wt_fut(sp,i)(1,nages(sp)));  // matrix multiplication A%*%B
        NByage0_hind2fut(itemp,sp,i+nyrs)(1,nages(sp))       = N_fut(sp,i)(1,nages(sp));   
        if (i==nyrs_fut)
        {   
           depletion0_fut(itemp,sp) = biomass_hind2fut(itemp,sp,nyrs_fut)/(1.e-20 + biomass_hind2fut(itemp,sp,1)); 
           if (msmMode>0)
             prey_consumed0_fut(itemp,sp,i) =sum(elem_prod(wt_fut(sp,i), elem_prod(elem_div(M2_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)))));// This makes sense since M2 is already summed over sources of predation...
        }//end if 
      }
      else
      {
        //fished biomass
        catch_hat_fut(itemp,sp,i)(1,nages(sp))              = elem_prod(elem_div(F_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)));
        Frate_fut(itemp,sp,i)                               = mean(elem_div(F_fut(sp,i),fsh_sel(sp)));
        tc_hat_fut(itemp,sp,i)                              = sum(catch_hat_fut(itemp,sp,i)(1,nages(sp)));
        fsh_age_hat_fut(sp,i)                               = catch_hat_fut(itemp,sp,i)(1,nages(sp)) / tc_hat_fut(itemp,sp,i);
        tc_biom_hat_fut(itemp,sp,i)                         = catch_hat_fut(itemp,sp,i)(1,nages(sp)) * wt_fut(sp,i)(1,nages(sp));// matrix multiplication A%*%B
        biomass_hind2fut(itemp,sp,i+nyrs)                   = N_fut(sp,i) * wt_fut(sp,i)(1,nages(sp));  // matrix multiplication A%*%B
        biomassSSB_hind2fut(itemp,sp,i+nyrs)                = elem_prod(N_fut(sp,i),pmature(sp))*wt_fut(sp,i)(1,nages(sp));
        biomassByage_hind2fut(itemp,sp,i+nyrs)(1,nages(sp)) = elem_prod(N_fut(sp,i), wt_fut(sp,i)(1,nages(sp)));
        NByage_hind2fut(itemp,sp,i+nyrs)(1,nages(sp))       = N_fut(sp,i)(1,nages(sp));
        if (i==nyrs_fut)
        {   
          depletion_fut(itemp,sp) = biomass_hind2fut(itemp,sp,nyrs_fut)/(1.e-20 + biomass_hind2fut(itemp,sp,1)); 
          if (msmMode>0)
            prey_consumed_fut(itemp,sp,i) =sum(elem_prod(wt_fut(sp,i), elem_prod(elem_div(M2_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)))));// This makes sense since M2 is already summed over sources of predation...
        }//end if 
      } // Compute time series of biomass consumed by different species (vector x vector = scalar)
    }//end sp
  }  // end fut_pass_number



FUNCTION void CALC_N_AT_AGE(int pass_number)
  for (int sp=1;sp<=nspp;sp++)
  { 
    R(sp) = mfexp(ln_mn_rec(sp) + rec_dev(sp)(1,nyrs) );
    for (int yr=1;yr<=nyrs;yr++)
      N(sp,yr,1) = R(sp,yr);       // Recruitment  
    N(sp,nyrs+1,1) = R(sp,nyrs); // Recruitment  in nyrs+1  
    for (int sp_age=2;sp_age<=nages(sp);sp_age++)
      N(sp,1,sp_age) = mfexp(ln_mn_rec(sp) - M1(sp,sp_age)*double(sp_age-1)  + init_dev(sp,sp_age) );  // Initial sp_age composition (first year) OjO this may have a weird effect 
      // Plus group
    N(sp,1,nages(sp)) /= (1.-mfexp(-M1(sp,nages(sp))));   
    switch (pass_number)
    {
      case 1:
        AvgN(sp,1) = N(sp,1);
        break;
      case 2:
        AvgN(sp,1) = elem_div( elem_prod( N(sp,1),(1.-S(sp,1))) ,Z(sp,1));
      default:
        AvgN(sp,1) = elem_div( elem_prod( N(sp,1),(1.-S(sp,1))) ,Z(sp,1));
      break;
    }// end switch pass number
  } // End loop over species
  M2.initialize();
  avail_food.initialize();
  // Qhat.initialize();
  // NEaten.initialize();

  // NEaten_hat.initialize();
  // ration_hat.initialize();
  // other_prey_eaten.initialize();
  // mn_Uhat.initialize();
  // Uhat.initialize();
  // KAge_hat.initialize();
  // mnKAge_hat.initialize();
  
  // for each year calculate the predation mortality
  for (i=1;i<=nyrs;i++)
  {
    if (msmMode>0 && current_phase()>dietphase)
    {
      CALC_AVAIL_FOOD(msmMode);
      // removed A 
      if (pass_number>1)
      {  
        switch (msmMode)
        {
          // CALC_M2(fut_pass_number): value of 0 = estimation mode, value of 1 = future projection mode
          case 2:
            CALC_M2(0);  
          break;
          default:
            CALC_M2(0);
          break;
        }
      }  
    }

    for (k=1;k<=nspp;k++)
    {
      CALC_MORT(0);
      N(k,i+1)(2,nages(k)) = ++elem_prod(N(k,i)(1,nages(k)-1),S(k,i)(1,nages(k)-1)) ;
      N(k,i+1,nages(k))   += S(k,i,nages(k))*N(k,i,nages(k));
      if(i<nyrs){
        switch (pass_number)
        {
          case 1:
            AvgN(k,i+1) = N(k,i+1);
          break;
          case 2:
            AvgN(k,i+1) = elem_div( elem_prod( N(k,i+1),(1.-S(k,i+1))) ,Z(k,i+1));
          default:
            AvgN(k,i+1) = elem_div( elem_prod( N(k,i+1),(1.-S(k,i+1))) ,Z(k,i+1));
          break;
        } 
      } 
    } 
  } 
  //  if(pass_number>1)
  //    CALC_SUIT(msmMode);

FUNCTION void CALC_N_AT_AGE_FUT(int pass_number) 
  //ll (simulation number) is passed to this function
  //i (future year) is passed to this function   
  
  dvariable B2B40;

  // get it started using ending values from estimation for ages 2-nages
  if(i==1&pass_number==1){
    for (int sp=1;sp<=nspp;sp++)
      N_fut(sp,1)= N(sp,nyrs);  
      //N_fut(sp,1)= N(sp,nyrs+1);  // get it started using ending values from estimation for ages 2-nages
  }  
      // now calculate recruitment of N_fut(sp,1,1), need AvgN for this so first pass calc_rec as mean Rec then iterate to find true rec

  switch (pass_number)
  {
    case 1:
     //R_fut(itemp)=CALC_REC_FUT(recMode);
      CALC_REC_FUT(0); // mean recruitment
    break;

    case 2:
     //R_fut(itemp)=CALC_REC_FUT(recMode);
      CALC_REC_FUT(recMode);
    break;

    default:
      //R_fut(itemp)=CALC_REC_FUT(recMode);
      CALC_REC_FUT(recMode);
    break;
  }// end switch pass number

  for (int sp=1;sp<=nspp;sp++)
  {  
    // FUTURE YEAR = 1 :Set first year of N_fut equal to N_nyrs+1 from the estimation mode 
    N_fut(sp,1)= N(sp,nyrs);  // repeated from above?

    // FUTURE YEAR > 1 if not first year of projection...
    if(i>1) N_fut(sp,i,1) = R_fut(itemp,sp,i);                    // Age 1 = Recruitment in >nyrs+1  
    N_fut(sp,nyrs_fut+1,1) =R_fut(itemp,sp,nyrs_fut);             // Age 1 = Recruitment in nyrs_fut+1  

    
    switch (pass_number)
    {
      case 1:
        // if first time through calcs use last est year's S as a starting place
        AvgN_fut(sp,i) = elem_div( elem_prod( N_fut(sp,i),(1.-S(sp,nyrs))) ,Z(sp,nyrs));
      break;
      case 2:
        // if >1 time through calcs use updated S
        AvgN_fut(sp,i) = elem_div( elem_prod( N_fut(sp,i),(1.-S_fut(sp,i))) ,Z_fut(sp,i));
      break;

      default:
        AvgN_fut(sp,i) = elem_div( elem_prod( N_fut(sp,i),(1.-S_fut(sp,i))) ,Z_fut(sp,i));
      break;
    }
  } 


  // If converting ABC to catch with function:
  if(harvestMode==12|harvestMode==13){
    dvariable sumabc=0;
    // For both 12 & 13 calculate the FABC and ABC_fut from the Tier 3 harvest control rule
    for (k=1;k<=nspp;k++) CALC_MORT(harvestMode);
      // assign ABC (will be different than catch)
      for (k=1;k<=nspp;k++){
        double Cbb;
        Cbb=.2;
        if(k==3) Cbb=0.;
        B2B40=0.;
        if(i==1) B2B40=biomassSSB(k,nyrs)/(Btarget*B0_set(k)); 
        if(i>1)  B2B40=biomassSSB_hind2fut(itemp,k,i-1)/(Btarget*B0_set(k)); 

        Fabc(itemp,k,i)=TIER3_CR(F40(itemp,k,i),B2B40, alpha_ABC,Cbb);
        Fofl(itemp,k,i)=TIER3_CR(F35(itemp,k,i),B2B40, alpha_OFL,Cbb);
        PF(itemp,k,i) = Fabc(itemp,k,i);  // if set to 12, PF=ABC, if set to 13, PF will be the 2MT catch below
       
        CALC_MORT(harvestMode);  // will return tc_biom_hat_fut
        ABC_fut(itemp,k,i)=tc_biom_hat_fut(itemp,k,i);
        ABC_N_fut(itemp,k,i)=tc_hat_fut(itemp,k,i);
        ABC_N_age_fut(itemp,k,i)=catch_hat_fut(itemp,k,i);
      } 

      // write ABC to dat file for catch function
          ofstream abc_out("abc_out.dat");
            abc_out << "# catch_scen "<<endl;
            abc_out << catch_scen <<endl;
            abc_out << "#yr"<<endl;
            abc_out << i <<endl;
            abc_out << "#mode"<<endl;
            abc_out << msmMode <<endl;
            abc_out << "#climScn"<<endl;
            abc_out <<itemp <<endl;
            abc_out << "#ABC for each spp"<<endl;
            for (k=1;k<=nspp;k++){
              sumabc += ABC_fut(itemp,k,i);
              abc_out << ABC_fut(itemp,k,i) <<endl;
            }
          abc_out.close(); 
        
        dvariable  killrun;
        ifstream tmpkillrun("killrun.dat");
          tmpkillrun >> killrun;
        tmpkillrun.close(); 

        if(killrun>0.){
          cout<<"Simulations stopped using killrun.dat (value >0 )"<<endl;
          exit(1);
        } 
      // Run the ABC --> Catch function in R
      if(harvestMode==13){
        
        if(sumabc>0){ 
          cout<<" running ABC to CATCH function, mode= "<<msmMode<<", hcr= "<<harvestMode<<", rec ="<<recMode<<", sim="<<itemp<<", catch scenario = "<<catch_scen<<" ,year = "<<i<<endl;
            // run amanda function for this year
          #if defined _WIN32 || defined _WIN64
                system(adstring("Rscript ../Scripts/R_code/abc2c.R") ); //system(adstring("./../Scripts/Shell_scripts/abc2c.bat ") );        
          #else 
                system(adstring("Rscript ../Scripts/R_code/abc2c.R") );
                // vanilla R
                // cpp R
                // littler   TRY this! and install
                //system(adstring("R --no-save --quiet --slave < ../Scripts/R_code/abc2c.R") );   
          #endif  

        }else{
          cout<<" skipped catch function; sum of ABCs = 0, sim="<<itemp<<" ,year = "<<i<<endl;
        }
    // read in the err code
        dvariable  tmpdat;
        ifstream rdtmpdat("tmp.dat");
          rdtmpdat >> tmpdat;
        rdtmpdat.close(); 
        if(tmpdat>10.){
          cout<<"ERROR with Catch Function (did not run/ timed out), rdtmpdat = "<<tmpdat<<endl;
          exit(1);
        } 
    
    // read in the catch
        dvar_vector  tmpCC(1,nspp);
        if(sumabc>0){ 
          //read in catch data from file
          ifstream catchIn("catch_out.dat");
          catchIn >> tmpCC;
          catchIn.close(); 
        }else{
          for (k=1;k<=nspp;k++) tmpCC(k)=0;
        }
    
      // find the corresponding F rate that yeilds the catch above
        for (k=1;k<=nspp;k++)
        {
          if(tmpCC(k)>0){
            PF(itemp,k,i) = SolveF2_kir(N_fut(k,i),(tmpCC(k)),k); // assign effective F rate 
          }else{
            PF(itemp,k,i) =0;
          }
        }

    }
  } 

  // now recalculate mort based on catch from function
  for (k=1;k<=nspp;k++) CALC_MORT(harvestMode);

  if (msmMode>0 && current_phase()>dietphase)
  {
    if (pass_number>1)
    {
      CALC_AVAIL_FOOD_FUT(msmMode);
       // cout<<"av_food_fut (1, "<<i<<") "<<av_food_fut(1,i)<<endl; exit(1);
     // removed A2
      switch (msmMode)
      {
        case 2:
          CALC_M2(1);
        break;
        default:
          CALC_M2(1);
        break;
      }
    }
  }

  // now calculate numbers at age for the nex year
  for (k=1;k<=nspp;k++)
  {
    CALC_MORT(harvestMode);  // will return catch, numbers at age
    N_fut(k,i+1)(2,nages(k))= ++elem_prod(N_fut(k,i)(1,nages(k)-1),S_fut(k,i)(1,nages(k)-1)) ;
    N_fut(k,i+1,nages(k)) += S_fut(k,i,nages(k))*N_fut(k,i,nages(k)); // plus groupd
    if(i<nyrs_fut)
    {
      switch (pass_number)
      {
        case 1:
          AvgN_fut(k,i+1) = elem_div( elem_prod( N_fut(k,i+1),(1.-S(k,nyrs))) ,Z(k,nyrs));//N_fut(k,i+1);
          break;
        case 2:
          //AvgN_fut(k,i+1) = elem_div( elem_prod( N_fut(k,i+1),(1.-S_fut(k,i+1))) ,Z_fut(k,i+1));
          for (int sp_age=1;sp_age<=nages(k);sp_age++)
          {
            if(Z_fut(k,i+1,sp_age)>0)
            {
              AvgN_fut(k,i+1,sp_age) =((N_fut(k,i+1,sp_age))*(1.-S_fut(k,i+1,sp_age)))/Z_fut(k,i+1,sp_age);
            }
            else
            {
              AvgN_fut(k,i+1,sp_age) = N_fut(k,i+1,sp_age);
            }
          }
          break;
        }
      }
    }

FUNCTION void CALC_PREDICTED_VALS(int fut_pass_number)
  if (fut_pass_number==0)
  {
    prey_consumed.initialize();
    if(dump_rep) cout<<" sp yr calculations begin" <<endl;
    for (int sp=1;sp<=nspp;sp++)
    {
      catch_hat(sp) = elem_prod(elem_div(F(sp),Z(sp)) ,elem_prod(1.-mfexp(-Z(sp)) , N(sp)));
      for (int yr=1;yr<=nyrs;yr++)
      {
        tc_hat(sp,yr)       = sum(catch_hat(sp,yr));
        if(fsh_age_type(sp)==1) 
          fsh_age_hat(sp,yr)  = catch_hat(sp,yr) / tc_hat(sp,yr);
        else
          fsh_age_hat(sp,yr)  = catch_hat(sp,yr)*age_trans_matrix(sp) / tc_hat(sp,yr);
        tc_biom_hat(sp,yr)  = catch_hat(sp,yr) * wt(sp,yr)(1,nages(sp));// matrix multiplication A%*%B
        biomass(sp,yr)      = N(sp,yr) * wt(sp,yr)(1,nages(sp));  // matrix multiplication A%*%B
        biomassByage(sp,yr) = elem_prod(N(sp,yr), wt(sp,yr)(1,nages(sp)));
       // NByage(sp,yr)       = N(sp,yr);
        biomassSSB(sp,yr)   = elem_prod(N(sp,yr),pmature(sp))*wt(sp,yr)(1,nages(sp));
      }
      dvar_vector Stmp(1,nages(sp));
      dvar_vector tmp_age(1,nages(sp)); 
      

      int yr_ind;
      Stmp.initialize();
      tmp_age.initialize();
      if(debugg) cout<<"nyrs_srv_biom("<<sp<<") "<<nyrs_srv_biom(sp)<<endl;
      for (int yr=1;yr<=nyrs_srv_biom(sp);yr++)
      {      
        yr_ind            = yrs_srv_biom(sp,yr) - styr + 1; // convert years into indices for 1-=nyrs counting purposes 
        //Stmp              = (S(sp,yr_ind)+1.)/2.;// note: assumes survey occurs in mid-year (pow(,.5))      
        Stmp              = mfexp(-Z(sp,yr_ind)*.5);// note: assumes survey occurs in mid-year      
        srv_biom_hat(sp,yr) = elem_prod(srv_q(sp) * srv_sel(sp) , elem_prod( Stmp ,N(sp,yr_ind)) ) * wt(sp,yr_ind)(1,nages(sp));
      }  
            if(debugg) cout<<" here 2" <<endl;
      for (int yr=1;yr<=nyrs_srv_age(sp);yr++)
      {
        yr_ind  = yrs_srv_age(sp,yr) - styr + 1;// convert years into indices for 1-=nyrs counting purposes    
        tmp_age = N(sp,yr_ind);// need to test for type of sp_age data here (if length bins are different
        tmp_age = elem_prod(srv_sel(sp) , tmp_age);
        if ( srv_age_type(sp)==1)
          srv_age_hat(sp,yr) = tmp_age / sum(tmp_age);
        else
        {
          srv_age_hat(sp,yr)  = tmp_age * age_trans_matrix(sp);
          srv_age_hat(sp,yr) /= sum(srv_age_hat(sp,yr));
        }
      }
            if(debugg) cout<<" here 3" <<endl;
      // Compute time series of biomass consumed by different species (vector x vector = scalar)
      if (sd_phase())
      {   
        depletion(sp) = biomass(sp,nyrs)/(1.e-20 + biomass(sp,1)); 
        if (msmMode>0)
          if (sp<=nspp)
            for (int yr=1;yr<=nyrs;yr++)
              prey_consumed(sp,yr) =sum(elem_prod(wt(sp,yr), elem_prod(elem_div(M2(sp,yr),Z(sp,yr)) ,elem_prod(1.-mfexp(-Z(sp,yr)) , N(sp,yr)))));
              // This makes sense since M2 is already summed over sources of predation...
      }//end if sd_phase()
            if(debugg) cout<<" here 4" <<endl;
    }
    // Hard wired for acoustic survey of pollock only...
          if(debugg) cout<<" here 5" <<endl;
    for (int yr=1;yr<=n_eit;yr++)
    {
      //convert years into indices of counting purposes
      double eit_yrs_ind = yrs_eit(yr) -styr +1;
      eit_age_hat(yr)     = elem_prod((N(1,eit_yrs_ind)),eit_sel(eit_yrs_ind)*eit_q); // No adjustment to mid-year
      eit_hat(yr)         = eit_age_hat(yr) *wt(1,eit_yrs_ind)(1,12);
      eit_age_hat(yr)     /= sum(eit_age_hat(yr));
    }
          if(debugg) cout<<" here 6" <<endl;
  }
  if(fut_pass_number==1)
  {
    for (int sp=1;sp<=nspp;sp++)
    {
      if(i==1)
      {
        // initalize the historical timeseries with est values
        for (int hyr=1;hyr<=nyrs;hyr++)
        {
          biomass0_hind2fut(itemp,sp,hyr)                  = biomass(sp,hyr);
          biomass_hind2fut(itemp,sp,hyr)                   = biomass(sp,hyr);
          biomassSSB_hind2fut(itemp,sp,hyr)                = biomassSSB(sp,hyr);
          biomassSSB0_hind2fut(itemp,sp,hyr)               = biomassSSB(sp,hyr);
          R0_hind2fut(itemp,sp,hyr)                        = R(sp,hyr);
          biomassByage_hind2fut(itemp,sp,hyr)(1,nages(sp))  = biomassByage(sp,hyr);
          biomassByage0_hind2fut(itemp,sp,hyr)(1,nages(sp)) = biomassByage(sp,hyr);
          NByage_hind2fut(itemp,sp,hyr)(1,nages(sp))        = N(sp,hyr);
          NByage0_hind2fut(itemp,sp,hyr)(1,nages(sp))       = N(sp,hyr);
        }

      }
      if(harvestMode==1)
      {
        // No harvest

        //unfished biomass
        //        catch_hat_fut(itemp,sp,i)(1,nages(sp))                       = elem_prod(elem_div(F_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)));
        //        tc_hat_fut(itemp,sp,i)                                               = sum(catch_hat_fut(itemp,sp,i)(1,nages(sp)));
        //        fsh_age_hat_fut(sp,i)                                        = catch_hat_fut(itemp,sp,i)(1,nages(sp)) / tc_hat_fut(itemp,sp,i);
        //        tc_biom_hat_fut(itemp,sp,i+nyrs)                                    = catch_hat_fut(itemp,sp,i)(1,nages(sp)) * wt(sp,nyrs)(1,nages(sp));// matrix multiplication A%*%B
        biomass0_hind2fut(itemp,sp,i+nyrs)                  = N_fut(sp,i) * wt_fut(sp,i)(1,nages(sp));  // matrix multiplication A%*%B
        biomassSSB0_hind2fut(itemp,sp,i+nyrs)               = elem_prod(N_fut(sp,i),pmature(sp))*wt_fut(sp,i)(1,nages(sp));   // matrix multiplication A%*%B
        R0_hind2fut(itemp,sp,i+nyrs)                        = R_fut(itemp,sp,i);   // matrix multiplication A%*%B 
        biomassByage0_hind2fut(itemp,sp,i+nyrs)(1,nages(sp)) = elem_prod(N_fut(sp,i), wt_fut(sp,i)(1,nages(sp)));  // matrix multiplication A%*%B
        NByage0_hind2fut(itemp,sp,i+nyrs)(1,nages(sp))       = N_fut(sp,i)(1,nages(sp));   
        if (i==nyrs_fut)
        {   
           depletion0_fut(itemp,sp) = biomass_hind2fut(itemp,sp,nyrs_fut)/(1.e-20 + biomass_hind2fut(itemp,sp,1)); 
           if (msmMode>0)
             prey_consumed0_fut(itemp,sp,i) =sum(elem_prod(wt_fut(sp,i), elem_prod(elem_div(M2_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)))));// This makes sense since M2 is already summed over sources of predation...
        }//end if 
      }
      else
      {
        //fished biomass
        catch_hat_fut(itemp,sp,i)(1,nages(sp))        = elem_prod(elem_div(F_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)));
        Frate_fut(itemp,sp,i)                         = mean(elem_div(F_fut(sp,i),fsh_sel(sp)));
        tc_hat_fut(itemp,sp,i)                        = sum(catch_hat_fut(itemp,sp,i)(1,nages(sp)));
        fsh_age_hat_fut(sp,i)                         = catch_hat_fut(itemp,sp,i)(1,nages(sp)) / tc_hat_fut(itemp,sp,i);
        tc_biom_hat_fut(itemp,sp,i)                   = catch_hat_fut(itemp,sp,i)(1,nages(sp)) * wt_fut(sp,i)(1,nages(sp));// matrix multiplication A%*%B
        biomass_hind2fut(itemp,sp,i+nyrs)                  = N_fut(sp,i) * wt_fut(sp,i)(1,nages(sp));  // matrix multiplication A%*%B
        biomassSSB_hind2fut(itemp,sp,i+nyrs)               = elem_prod(N_fut(sp,i),pmature(sp))*wt_fut(sp,i)(1,nages(sp));
        biomassByage_hind2fut(itemp,sp,i+nyrs)(1,nages(sp)) = elem_prod(N_fut(sp,i), wt_fut(sp,i)(1,nages(sp)));
        NByage_hind2fut(itemp,sp,i+nyrs)(1,nages(sp))       = N_fut(sp,i)(1,nages(sp));
        if (i==nyrs_fut)
        {   
          depletion_fut(itemp,sp) = biomass_hind2fut(itemp,sp,nyrs_fut)/(1.e-20 + biomass_hind2fut(itemp,sp,1)); 
          if (msmMode>0)
            prey_consumed_fut(itemp,sp,i) =sum(elem_prod(wt_fut(sp,i), elem_prod(elem_div(M2_fut(sp,i),Z_fut(sp,i)) ,elem_prod(1.-mfexp(-Z_fut(sp,i)) , N_fut(sp,i)))));// This makes sense since M2 is already summed over sources of predation...
        }//end if 
      } // Compute time series of biomass consumed by different species (vector x vector = scalar)
    }//end sp
  }  // end fut_pass_number


FUNCTION void CALC_U_AGE(int pass_number)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////based on sp_age Calculate (U) effective prey size freq (Npa/sum(N))*switch -- rows=prey, cols=preds, can be in calcs section if Npa/Np doesn't change each year
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  dvar_matrix preypropN(1,nspp,1,nages);
  dvar_matrix preypropW(1,nspp,1,nages);
  dvar_matrix preypropUSE(1,nspp,1,nages);
  preypropUSE.initialize();
  preypropN.initialize();
  preypropW.initialize();
  UnewAge.initialize();
  for(int predd=1;predd<=nspp;predd++)
  {
    for(int preyy=1;preyy<=nspp;preyy++)
    {
      if(Umode==0){
        UnewAge(predd,preyy)=UobsAge(predd,preyy);
      }
      else
      {
        dvar_vector preyprop(1,nages(preyy));
        dvar_vector pref(1,nages(predd));
        preyprop.initialize();
        pref.initialize();
        pref=1;  
        cout<<"CHANGE THIS TO / mean() from / max()!!!"<<endl;
        if(pass_number==0)
        {
          if (sum(N(preyy,i))>0)
            preyprop=N(preyy,i)/max(N(preyy,i));  
          else
            preyprop(preyy)=0.;
          if(useWt(predd)==1)
          {// if the model is set to forage as a function of biomass....
            preyprop=elem_prod(N(preyy,i),(wt(preyy,i)));
            if (sum(preyprop)>0)
              preyprop=preyprop/max(preyprop);// proportion by weight <- change this to propby weight by year !!!!  ? need to pre-allocate?        
            else
              preyprop(preyy)=0.;      
          }   
        }
        if(pass_number==1){
          if (sum(N_fut(preyy,i))>0)
            preyprop=N_fut(preyy,i)/max(N_fut(preyy,i));  
          else
            preyprop(preyy)=0.;
          if(useWt(predd)==1)
          {// if the model is set to forage as a function of biomass....
            preyprop=elem_prod(N_fut(preyy,i),(wt_fut(preyy,i)));
            if (sum(preyprop)>0)
              preyprop=preyprop/max(preyprop);// proportion by weight <- change this to propby weight by year !!!!  ? need to pre-allocate?        
            else
              preyprop(preyy)=0.;      
          }   
        }
        preypropUSE(preyy)=preyprop;      
        dvar_matrix Utmp(1,nages(predd),1,nages(preyy)); //U(#yrs,pred#,Prey#,preysize,predsize)
        dvar_matrix switch1(1,nages(predd),1,nages(preyy)); //U(#yrs,pred#,Prey#,preysize,predsize)
        dvar_matrix avail(1,nages(predd),1,nages(preyy));
        dvar_matrix part1(1,nages(predd),1,nages(preyy));
        dvar_matrix ttmp(1,nages(predd),1,nages(preyy));
        Utmp.initialize();
        switch1.initialize();
        avail.initialize();
        part1.initialize();
        ttmp.initialize();
        switch1=1; //pre-allocate a vector of 1s 
        part1=1;  
        for(int predL=1;predL<=nages(predd);predL++) // for each predL size of predd
        {  
          ttmp(predL)=(LbyAge(preyy,i)-l2(predd,preyy)*LimitAge(predd,i,predL))/l1(predd,preyy);
          for(int preyL=1;preyL<=nages(preyy);preyL++)  
            if(ttmp(predL,preyL)<0)
              ttmp(predL,preyL)=0.;        
          switch1(predL)(1,nages(preyy))=mfexp(-ttmp(predL));
        }
        for(int predL=1;predL<=nages(predd);predL++) // for each predL size of predd
        {
          avail(predL)=elem_prod(switch1(predL),preypropUSE(preyy));
          if(sum(avail(predL))>0)
            avail(predL)=avail(predL)/sum(avail(predL));
          else
            avail(predL)=0.;
          part1(predL)=0.;
          for(int preyL=1;preyL<=nages(preyy);preyL++) // for each j size of preyy
            if(avail(predL,preyL)>0)
              part1(predL,preyL)=( prefa(predd,preyy)*(pow(avail(predL,preyL),prefb(predd,preyy))) )/(1+(prefa(predd,preyy)*(pow(avail(predL,preyL),prefb(predd,preyy)))) ) ;
          if(sum(part1(predL))>0)
            part1(predL)=part1(predL)/(sum(part1(predL)));// standardize across all prey sizes that the pred predL will eat.        
          else
            part1(predL)=0.;
          Utmp(predL)=KAge(predd,preyy,predL)*part1(predL);
          UnewAge(predd,preyy,predL)(1,nages(preyy))=KAge(predd,preyy,predL)*part1(predL);
          //Utmp is a vector of 1xlengths (predL) of pred lengths, switch is a matrix of prey X pred L, and K is a 1xk vector of diet proportions  
        }//end pred lengths k
      }// end if Umode==0
      for(int predL=1;predL<=nages(predd);predL++)
        if(sum(UnewAge(predd,preyy,predL))>0)
          UnewAge(predd,preyy,predL)=(UnewAge(predd,preyy,predL)/sum(UnewAge(predd,preyy,predL)))*KAge(predd,preyy,predL);
        else
          UnewAge(predd,preyy,predL)=0.;
      UnewAVGAge(predd,preyy)+=UnewAge(predd,preyy);
      if(i==nyrs)
        UnewAVGAge(predd,preyy)=UnewAVGAge(predd,preyy)/nyrs;
    }//end prey  
  }// end pred   
FUNCTION void CALC_CONSUM_AGE(int C_number)
  // Not sure this is used anymore
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// Calculate annual ration2 C=24*0.0134*mfexp(0.0115*TempC)*fdays(predd)*S  units are g of food/g of predator
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  for(i=1;i<=nyrs;i++)
  //  {
  // this appears to be working correctly - Kir has checked the U and E matrices, as well as ration 
  //cout<<"CALC_CONSUM_AGE"<<endl;
 /* if(C_number==0)
  {
      for(int predd=1;predd<=nspp;predd++)
      {// for each pred predd  
      
        ConsumAge(predd,i)=24.*0.0134*mfexp(0.0115*TempC(i))*91.25*elem_prod(S2Age(predd,i),wt(predd,i));//kg/pred.yr
        Consum_livingAge(predd,i)=ConsumAge(predd,i); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        if(C_model(predd)==1){
          ConsumAge(predd,i)=elem_prod(((CA(predd))*pow(wt(predd,i)*1000.,CB(predd))*fT(predd,i)*fdays(predd)),wt(predd,i)*1000.);//g/pred.yr      
          //ConsumAge(predd,i)=elem_prod(ConsumAge(predd,i),Pvalue(predd)*PvalueAdj(predd)*PAge(predd));
          ConsumAge(predd,i)=elem_prod(ConsumAge(predd,i),Pvalue(predd)*PvalueAdj(predd)*Pby_yr(predd,i));
        }
        ration2Age(predd,i)=ConsumAge(predd,i)/1000.;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
      }//end predd  
  }
  if(C_number==1){
    for(int predd=1;predd<=nspp;predd++)
    {// for each pred predd  
      ConsumAge_fut(predd,i)=24.*0.0134*mfexp(0.0115*TempC_futUSE(ll,i))*91.25*elem_prod(S2Age(predd,nyrs),wt(predd,nyrs));//kg/pred.yr
      Consum_livingAge_fut(predd,i)=ConsumAge_fut(predd,i); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
      if(C_model(predd)==1){
        ConsumAge_fut(predd,i)=elem_prod(((CA(predd))*pow(wt(predd,nyrs)*1000.,CB(predd))*fT_fut(ll,predd,i)*fdays(predd)),wt(predd,nyrs)*1000.);//g/pred.yr      
        ConsumAge_fut(predd,i)=elem_prod(ConsumAge_fut(predd,i),Pvalue(predd)*PvalueAdj(predd)*PAge(predd));       
      }
      ration2Age_fut(predd,i)=ConsumAge_fut(predd,i)/1000.;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
    }//end predd
  }*/
  
 CALC_RATION(C_number);
  if(C_number==1){
    for (int preyy=1;preyy<=nspp;preyy++)
      availBbyAge_fut(preyy,i)=elem_prod(AvgN_fut(preyy,i),wt_fut(preyy,i));  //kgB/1000; //kg  
  }else{
    for (int preyy=1;preyy<=nspp;preyy++)
      availBbyAge(preyy,i)=elem_prod(AvgN(preyy,i),wt(preyy,i));  //kgB/1000; //kg  
  }
  dvar_matrix U4bio(1,maxA,1,maxA); // demand
  U4bio.initialize();  
  for (int predd=1;predd<=nspp;predd++)
    for (int preyy=1;preyy<=nspp;preyy++)
      U4bio+=UnewAge(predd,preyy);
  dvar_vector ttmp(1,maxA);
  ttmp.initialize();  
  ttmp=colsum(U4bio);
  for (int preyy=1;preyy<=nspp;preyy++)
  {
    EbyAge(preyy,i).initialize();
    B4EbyAge(preyy,i).initialize();
    pred_EbyAge(preyy,i).initialize();
    //availB(preyy,i).initialize();
    Ucheck(preyy,i).initialize();
    Bage_other.initialize();
    for (int predd=1;predd<=nspp;predd++)
    {    
      dvar_matrix UcolSums(1,nages(predd),1,nages(preyy)); // demand
      dvar_matrix UcolSums2(1,nages(predd),1,nages(preyy)); // demand
      dvar_matrix Etmp3(1,nages(predd),1,nages(preyy)); // demand
      dvar_matrix Btmp3(1,nages(predd),1,nages(preyy)); // available B
      dvar_matrix Utmp3(1,nages(predd),1,nages(preyy)); // available B
      UcolSums.initialize();  
      UcolSums2.initialize();  
      Etmp3.initialize();  
      Btmp3.initialize();  
      Utmp3.initialize();  
      Bage_other(predd)=(1.-colsum(KAge(predd)))*other_food(predd);
      for(int pred_age=1;pred_age<=nages(predd);pred_age++)
        UcolSums2(pred_age)=UnewAge(predd,preyy,pred_age)(1,nages(preyy));
      for(int pred_age=1;pred_age<=nages(predd);pred_age++)
        for(int prey_age=1;prey_age<=nages(preyy);prey_age++)
          if(ttmp(prey_age)>0)
            UcolSums(pred_age,prey_age)=(UcolSums2(pred_age,prey_age)/ttmp(prey_age));
      for(int pred_age=1;pred_age<=nages(predd);pred_age++){
        //Etmp3(pred_age)=UnewAge(predd,preyy,pred_age)(1,nages(preyy))*ration2Age(predd,i,pred_age)*AvgN(predd,i,pred_age);

        switch(C_number)
        {
          default:
            cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;
            cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;    
            cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;
          break;
          case 0:
            Etmp3(pred_age)=UnewAge(predd,preyy,pred_age)(1,nages(preyy))*ration2Age(predd,i,pred_age)*AvgN(predd,i,pred_age);
            Btmp3(pred_age)=elem_prod(UcolSums(pred_age),elem_prod(wt(preyy,i),AvgN(preyy,i)));
            Utmp3(pred_age)=UnewAge(predd,preyy,pred_age)(1,nages(preyy));
          break;
          case 1:
            Etmp3(pred_age)=UnewAge(predd,preyy,pred_age)(1,nages(preyy))*ration2Age_fut(predd,i,pred_age)*AvgN_fut(predd,i,pred_age);
            Btmp3(pred_age)=elem_prod(UcolSums(pred_age),elem_prod(wt_fut(preyy,i),AvgN_fut(preyy,i)));
            Utmp3(pred_age)=UnewAge(predd,preyy,pred_age)(1,nages(preyy));
          break;
        }
      }  
      //EbyAge(preyy,i)+=colsum(Etmp3); //E(preyy,i,prey_age),Unew(i,predd,preyy,pred_age,prey_age) // should be adding across predators?
      //EbyAge(preyy,i)+=elem_prod(colsum(Etmp3),elem_div(availBbyAge(preyy,i),availBbyAge(preyy,i)+funRespA(predd)*other_food(predd))); //E(preyy,i,prey_age),Unew(i,predd,preyy,pred_age,prey_age) // should be adding across predators?
      switch(C_number)
        {
          default:
            cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;
            cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;    
            cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;
          break;
          case 0:
            EbyAge(preyy,i)+=colsum(Etmp3);
            B4EbyAge(preyy,i)+=colsum(Btmp3);
            UcheckbyAge(preyy,i)+=colsum(Utmp3); //E(preyy,i,prey_age),Unew(i,predd,preyy,pred_age,prey_age) // should be adding across predators?
            //availB(preyy,i)+=colsum(Btmp3); //E(preyy,i,prey_age),Unew(i,predd,preyy,pred_age,prey_age) // should be adding across predators?
            pred_EbyAge(preyy,i,predd)=sum(colsum(Etmp3));  
            pred_BbyAge(preyy,i,predd)=sum(colsum(Btmp3));  
          break;
          case 1:
            EbyAge_fut(preyy,i)+=colsum(Etmp3);
            B4EbyAge_fut(preyy,i)+=colsum(Btmp3);
            UcheckbyAge_fut(preyy,i)+=colsum(Utmp3); //E(preyy,i,prey_age),Unew(i,predd,preyy,pred_age,prey_age) // should be adding across predators?
            //availB(preyy,i)+=colsum(Btmp3); //E(preyy,i,prey_age),Unew(i,predd,preyy,pred_age,prey_age) // should be adding across predators?
            pred_EbyAge_fut(preyy,i,predd)=sum(colsum(Etmp3));  
            pred_BbyAge_fut(preyy,i,predd)=sum(colsum(Btmp3));  
          break;
        }
    }//end for pred predd
    
  }// end for preyy
  switch(C_number)
  {
    default:
      cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;
      cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;    
      cout<<"C number must be 0 or 1 !!!!!!!!!+==========================="<<endl;
    break;
    case 0:
      for(int preyy=1;preyy<=nspp;preyy++)
      {
        dvar_vector Eagetmp(1,nages(preyy));
        Eagetmp.initialize();
        //Nage(preyy,i)=AvgN(preyy,i);
        Nage(preyy,i)=N(preyy,i);
        Eage(preyy,i)=EbyAge(preyy,i);
        //Bage(preyy,i)=B4EbyAge(preyy,i);
        Bage(preyy,i)=elem_prod(N(preyy,i),wt(preyy,i));
        EageN(preyy,i)=elem_div(Eage(preyy,i),wt(preyy,i));  // N consumed
        BageN(preyy,i)=elem_div(Bage(preyy,i),wt(preyy,i));  // N available
        Eagetmp=Eage(preyy,i); // temporary Eaten (aka demand) by sp_age matrix for the preyy spp and the year
      }  
    break;
    case 1:
      for(int preyy=1;preyy<=nspp;preyy++)
      {
        dvar_vector Eagetmp(1,nages(preyy));
        Eagetmp.initialize();
        //Nage(preyy,i)=AvgN(preyy,i);
        Nage_fut(preyy,i)=N_fut(preyy,i);
        Eage_fut(preyy,i)=EbyAge_fut(preyy,i);
        //Bage(preyy,i)=B4EbyAge(preyy,i);
        Bage_fut(preyy,i)=elem_prod(N_fut(preyy,i),wt_fut(preyy,i));
        EageN_fut(preyy,i)=elem_div(Eage_fut(preyy,i),wt_fut(preyy,i));  // N consumed
        BageN_fut(preyy,i)=elem_div(Bage_fut(preyy,i),wt_fut(preyy,i));  // N available
        Eagetmp=Eage_fut(preyy,i); // temporary Eaten (aka demand) by sp_age matrix for the preyy spp and the year
      }  
    break;
  }

FUNCTION void CALC_UandC(int pass_number)
  int nyr=nyrs;
  if(pass_number)
    nyr=nyrs_fut;
  for (i=1;i<=nyr;i++){ // year
    if(RationByAge==1){
      CALC_U_AGE(pass_number);        // calculate the proportion of each ration that goes to each prey of each size (pred and pred size specific)
    CALC_CONSUM_AGE(pass_number);      // calculate consumption (ration)
    }else{
      CALC_U(pass_number);        // calculate the proportion of each ration that goes to each prey of each size (pred and pred size specific)
      CALC_CONSUM(pass_number);      // calculate consumption (ration)
    }        
  }
FUNCTION void CALC_SUIT(int msmMode_pass_number)
  if(msmMode_pass_number>0){
   // cout<<" *******FIRST PART  =U/(W*N)*******"<<endl;
    dvariable suit_tmp ;
    stom_div_bio2.initialize();
    if(msmMode_pass_number==1)
    {
      // for (int year=1;year<=2;year++) // year
      // {
      //   for (int pred=1;pred<=nspp;pred++) //predator species
      //   { 
      //     for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age 
      //     {
      //       for (int prey=1;prey<=nspp;prey++) // prey species 
      //       {
      //         for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp_age 
      //         {  
      //          suit_tmp = stom(year,pred,pred_age,prey,prey_age)/(AvgN(prey,stomyrs(year),prey_age)); 
      //          if (mn_wt_stom(pred,pred_age,prey,prey_age)!=0.)
      //            stom_div_bio(year,pred,pred_age,prey,prey_age) = suit_tmp/mn_wt_stom(pred,pred_age,prey,prey_age);    
      //         }  //end prey sp_age 
      //       }  // end  prey sp
      //     }  // end pred sp_age    
      //   } // end  pred species 
      // } // end  year loop  
    }else{
      for (int year=1;year<=nyrs;year++) // year
      {
      for (int pred=1;pred<=nspp;pred++) //predator species
      { 
        for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age 
        {
        for (int prey=1;prey<=nspp;prey++) // prey species 
        {
          for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp_age 
          {  
           suit_tmp = stomKir(year,pred,pred_age,prey,prey_age)/(AvgN(prey,year,prey_age)); 
           if (wt(prey,year,prey_age)!=0.)
             stom_div_bio2(year,pred,pred_age,prey,prey_age) = suit_tmp/wt(prey,year,prey_age);    
  
  //         if(AvgN(prey,year,prey_age)!=0)
  //           stom_div_bio2(year,pred,pred_age,prey,prey_age) = (UobsAge(pred,prey,pred_age,prey_age)*ration2Age(pred,year,pred_age))/(AvgN(prey,year,prey_age)*wt(prey,year,prey_age));    
          }  //end prey sp_age 
        }  // end  prey sp
        }  // end pred sp_age    
      } // end  pred species 
      } // end  year loop  
   }
    //end if kirs_pass_num==0
    // cout<<"************suitabilities calculation*******"<<endl;
    suit_main.initialize(); // sets to zero
    dvariable suma_suit;
    for (int pred=1;pred<=nspp;pred++) //predator species
    {
    for (int pred_age=1;pred_age<=nages(pred);pred_age++) // predator sp 
    {
      for (int prey=1;prey<=nspp;prey++) // prey specie
      {
      for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp
      { 
        switch (msmMode_pass_number)
        {
          case 0:
          break;
          case 1:
            // for (int year=1;year<=nstom;year++) // year
            // {
            //   suma_suit=sum(stom_div_bio(year,pred,pred_age)); // sum of aitemp prey and prey ages in the stom of the pred sp_age j
            //   suit_std(year,pred,pred_age,prey,prey_age) = stom_div_bio(year,pred,pred_age,prey,prey_age)/(suma_suit+of_stom(year,pred,pred_age));  //(sum(stom_div_bio(year,pred,pred_age)+of_stom(year,pred,pred_age)));            
            //   suit_main(pred,pred_age,prey,prey_age)    += suit_std(year,pred,pred_age,prey,prey_age);                                  
            // }// end year
            // suit_main(pred,pred_age,prey,prey_age) /= nstom; 
            cout<<"!!!! THIS MODE NO LONGER USED (MSM MODE = 1) "<<endl;
          break;
          case 2:
            suit_main(pred,pred_age,prey,prey_age)=0;
            for (int year=1;year<=nyrs;year++) // year
            {
              suma_suit=sum(stom_div_bio2(year,pred,pred_age)); // sum of all prey and prey ages in the stom of the pred sp_age j
              suit_main(pred,pred_age,prey,prey_age)    += stom_div_bio2(year,pred,pred_age,prey,prey_age)/(suma_suit+of_stomKir(pred,pred_age,year));                               
            }// end year
            suit_main(pred,pred_age,prey,prey_age) /= nyrs; 
            
          break;
         // REMOVED A3
        }  //end switch              
      }// prey sp_age 
      } // prey
      suit_other(pred,pred_age) = 1. - sum(suit_main(pred,pred_age)); //estimate the other food suitability        
    }// pred sp_age
    }  // pred
  }//if msmMode >0
FUNCTION CALC_SELECTIVITY
  for (int sp=1;sp<=nspp;sp++)
  {
    // Fishery...    
    fsh_sel(sp)(1,nselages)          = fsh_sel_coff(sp);                          // fsh_sel_coff(sp) parameters
    fsh_sel(sp)(nselages+1,nages(sp)) = fsh_sel(sp,nselages);
    avgsel_fsh(sp)                   = log(mean(mfexp(fsh_sel_coff(sp))));
    fsh_sel(sp)                     -= log(mean(mfexp(fsh_sel(sp))));
    fsh_sel(sp)                      = mfexp(fsh_sel(sp));

    // Surveys... 
    // Do this only if logistic is selected for this species...
    if (logist_sel_phase(sp) >0){// logist_sel_phase(sp)= -2 so this is not done
      for (int sp_age=1;sp_age<=nages(sp);sp_age++)
        srv_sel(sp,sp_age) = 1/(1+mfexp(-srv_sel_slp(sp)*(double(sp_age)-srv_sel_inf(sp))));  
    }
    else // use coefficient form... these are the actual calculations
    {
      srv_sel(sp)(1,nselages)          = srv_sel_coff(sp);
      srv_sel(sp)(nselages+1,nages(sp)) = srv_sel(sp,nselages);
      avgsel_srv(sp)                   = log(mean(mfexp(srv_sel_coff(sp))));
      srv_sel(sp)                     -= log(mean(mfexp(srv_sel(sp))));
      srv_sel(sp)                      = mfexp(srv_sel(sp));
    }
  }
FUNCTION void CALC_AVAIL_FOOD(int msmMode_pass_number)
  // cout<<"CALC AVAILABLE FOOD ******************************************************************************"<<endl;
  for (int pred=1; pred<= nspp; pred++)  // predator loop
  {
    for (int pred_age=1;pred_age<=nages(pred);pred_age++)  // predator sp_age loop
    {
      dvariable tmpsuit;
      dvariable tmp_othersuit;
      tmp_othersuit=0.;
      tmpsuit=0.;
      tmpsuit=1-suit_other(pred,pred_age);
      avail_food(pred,i,pred_age)=0.;    
      if(msmMode_pass_number==1){
        for (int prey=1;prey<=nspp;prey++)
        {
          for (int prey_age =1;prey_age<=nages(prey);prey_age++)
            avail_food(pred,i,pred_age) += suit_main(pred,pred_age,prey,prey_age)*overlap(i,pred,prey) *AvgN(prey,i,prey_age) * mn_wt_stom(pred,pred_age,prey,prey_age);   
          tmp_othersuit+=sum(suit_main(pred,pred_age,prey))*overlap(i,pred,prey);
        }
        avail_food(pred,i,pred_age) += other_food(pred)*(1.-(tmp_othersuit));
          //avail_food(pred,i,pred_age) += other_food(pred)* suit_other(pred,pred_age); 
        }else{
          for (int prey=1;prey<=nspp;prey++)
          {
            for (int prey_age =1;prey_age<=nages(prey);prey_age++)
              avail_food(pred,i,pred_age) += suit_main(pred,pred_age,prey,prey_age)*overlap(i,pred,prey) *AvgN(prey,i,prey_age) *wt(prey,i,prey_age); 
            tmp_othersuit+=sum(suit_main(pred,pred_age,prey))*overlap(i,pred,prey);
          }
          avail_food(pred,i,pred_age) += other_food(pred)*(1.-(tmp_othersuit));
          
          //avail_food(pred,i,pred_age) += other_food(pred)* (suit_other(pred,pred_age));     
        }  
      }              // end pred sp_age loop
  }                // end pred loop
FUNCTION void CALC_AVAIL_FOOD_FUT(int msmMode_pass_number)
  // cout<<"CALC AVAILABLE FOOD ******************************************************************************"<<endl;
  for (int pred=1; pred<= nspp; pred++)  // predator loop
  {
    for (int pred_age=1;pred_age<=nages(pred);pred_age++)  // predator sp_age loop
    {
      dvariable tmpsuit;
      dvariable tmp_othersuit;
      tmpsuit=1-suit_other(pred,pred_age);
      av_food_fut(pred,i,pred_age)=0;
        if(msmMode_pass_number==1)
        {
          for (int prey=1;prey<=nspp;prey++)
          {
            for (int prey_age =1;prey_age<=nages(prey);prey_age++)
              av_food_fut(pred,i,pred_age) += suit_main(pred,pred_age,prey,prey_age)*overlap_fut(itemp,pred,prey,i) *AvgN_fut(prey,i,prey_age) * mn_wt_stom(pred,pred_age,prey,prey_age);   
            tmp_othersuit+=sum(suit_main(pred,pred_age,prey)*overlap_fut(itemp,pred,prey,i));
          }
        av_food_fut(pred,i,pred_age) += other_food(pred)*(1.-(tmp_othersuit)); 
        //av_food_fut(pred,i,pred_age) += other_food(pred)* suit_other(pred,pred_age); 
      }else{
        for (int prey=1;prey<=nspp;prey++)
        {
          for (int prey_age =1;prey_age<=nages(prey);prey_age++)
          {
            av_food_fut(pred,i,pred_age) += suit_main(pred,pred_age,prey,prey_age)*overlap_fut(itemp,pred,prey,i)*AvgN_fut(prey,i,prey_age) *wt_fut(prey,i,prey_age); 
            // cout <<pred<<" "<<tmp_othersuit<<" "<<overlap_fut(itemp,pred,prey,i)<<" "<<av_food_fut(pred,i,pred_age)<<" "<<wt_fut(prey,i,prey_age)<<endl;
          }
          tmp_othersuit+=sum(suit_main(pred,pred_age,prey)*overlap_fut(itemp,pred,prey,i));
        }
        av_food_fut(pred,i,pred_age) += other_food(pred)*(1.-(tmp_othersuit));
        //av_food_fut(pred,i,pred_age) += other_food(pred)* suit_other(pred,pred_age);     
      }  
    }// end pred sp_age loop
  }// end pred loop
  
FUNCTION void CALC_RATION(int fut_pass_number)
  int nyrsUSE;
  if(fut_pass_number==0)
    nyrsUSE=nyrs;
  if(fut_pass_number==1)  
    nyrsUSE=nyrs_fut;
  for(int yrr=1;yrr<=nyrsUSE;yrr++)  
  {
    for(int predd=1;predd<=nspp;predd++)
    {// for each pred predd  
      if(fut_pass_number==0){
        ConsumAge(predd,yrr)=24*0.0134*mfexp(0.0115*TempC(yrr))*91.25*elem_prod(S2Age(predd,yrr),wt(predd,yrr));//kg/pred.yr
        Consum_livingAge(predd,yrr)=ConsumAge(predd,yrr); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        if(C_model(predd)==1){
          ConsumAge(predd,yrr)=elem_prod(((CA(predd))*pow(wt(predd,yrr)*1000,CB(predd))*fT(predd,yrr)*fdays(predd)),wt(predd,yrr)*1000);//g/pred.yr      
          //ConsumAge(predd,yrr)=elem_prod(ConsumAge(predd,yrr),Pvalue(predd)*PAge(predd));
          ConsumAge(predd,yrr)=elem_prod(ConsumAge(predd,yrr),Pvalue(predd)*Pby_yr(predd,yrr));
        }
        ration2Age(predd,yrr)=ConsumAge(predd,yrr)/1000;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
      }
      if(fut_pass_number==1){  
        
        ConsumAge_fut(predd,yrr)=24*0.0134*mfexp(0.0115*TempC_futUSE(itemp,yrr))*91.25*elem_prod(S2Age(predd,nyrs),wt_fut(predd,yrr));//kg/pred.yr
        Consum_livingAge_fut(predd,yrr)=ConsumAge_fut(predd,yrr); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        if(C_model(predd)==1){
          ConsumAge_fut(predd,yrr)=elem_prod(((CA(predd))*pow(wt_fut(predd,yrr)*1000,CB(predd))*fT_fut(itemp,predd,yrr)*fdays(predd)),wt_fut(predd,yrr)*1000);//g/pred.yr      
          //ConsumAge_fut(predd,yrr)=elem_prod(ConsumAge_fut(predd,yrr),Pvalue(predd)*PAge(predd));
          ConsumAge_fut(predd,yrr)=elem_prod(ConsumAge_fut(predd,yrr),Pvalue(predd)*Pby_yr(predd,nyrs));
        }
        ration2Age_fut(predd,yrr)=ConsumAge_fut(predd,yrr)/1000;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
        
      }  
    }//end predd
  }
  
FUNCTION void CALC_M2(int fut_pass_number)
  CALC_RATION(fut_pass_number);
  for (int prey=1;prey<=nspp;prey++)  // prey spp loop
  {
    for (int prey_age =1;prey_age<=nages(prey);prey_age++)  // prey sp_age loop
    {
      dvariable Mtmp=0.; 
      dvariable Eatentmp=0.;  
      for (int pred=1;pred<=nspp;pred++)   // pred species loop
      {  
        for (int pred_age=1;pred_age<=nages(pred);pred_age++)  // Pred sp_age loop
        {
          if(fut_pass_number==0)
          {
            if(msmMode==1)
            {
              //Mtmp += (AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/avail_food(pred,i,pred_age);
            }else
            {
              Mtmp += (AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/avail_food(pred,i,pred_age);
              Eatentmp += (AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age));
            }
          }
          if(fut_pass_number==1)
          {
            if(msmMode==1)
            {
              //Mtmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration(pred,nyrs,pred_age)*suit_main(pred,pred_age,prey,prey_age))/av_food_fut(pred,i,pred_age);
            }else{
              Mtmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/avail_food(pred,nyrs,pred_age);
              Eatentmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age));
              // Mtmp += AvgN_fut(pred,i,pred_age)*overlap_fut(i,pred,prey) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/(av_food_fut(pred,i,pred_age));
            }
          }
          /*switch (Kirs)
          {
            case 0:
              if(fut_pass_number==0)
                Mtmp += AvgN(pred,i,pred_age) * ration(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/avail_food(pred,i,pred_age);
              if(fut_pass_number==1)
                Mtmp += AvgN_fut(pred,i,pred_age) * ration(pred,nyrs,pred_age)*suit_main(pred,pred_age,prey,prey_age)/av_food_fut(pred,i,pred_age);  
            break;
            case 1:
            
              if(fut_pass_number==0)
                Mtmp += AvgN(pred,i,pred_age) * ration2Age(pred,i,pred_age)*UobsAge(pred,prey,pred_age,prey_age)/avail_food(pred,i,pred_age);
              if(fut_pass_number==1)
                Mtmp += AvgN_fut(pred,i,pred_age) * ration2Age_fut(pred,i,pred_age)*UobsAge(pred,prey,pred_age,prey_age)/av_food_fut(pred,i,pred_age);
          
            break;
            default:
              if(fut_pass_number==0)
                Mtmp += AvgN(pred,i,pred_age) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/avail_food(pred,i,pred_age);
              if(fut_pass_number==1)
                Mtmp += AvgN_fut(pred,i,pred_age) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/av_food_fut(pred,i,pred_age);
            break;
          }   */     
        }  // end pred sp_age loop
      }  // end pred spp loop
      if(fut_pass_number==0){
        M2(prey,i,prey_age) = Mtmp; //M2_plk(prey,i,prey_age);
        B_eaten(prey,i,prey_age)=Eatentmp;
      } 
      if(fut_pass_number==1){
        M2_fut(prey,i,prey_age) = Mtmp; //M2_plk(prey,i,prey_age);  
        B_eaten_fut(prey,i,prey_age)=Eatentmp;
      }
    }   // end prey sp_age loop  
  }   // end prey spp loop
  
FUNCTION void CALC_MORT(int mort_pass_number)
  dvariable tmpC=0.;
  dvariable chattmp=0.;
  dvar_vector  tmpCC(1,3);
  dvar_vector  tmpF2(1,3);
  dvariable tmpF=0.;
  dvariable B2B40;
  tmpCC=0.;
  switch(mort_pass_number)
  {
    default:
      //hindcast of the model
      F(k,i) = fsh_sel(k) * mean_F(k) * mfexp(F_dev(k,i));
      Z(k,i) = F(k,i) + M1(k)(1,nages(k)) + M2(k,i);
      S(k,i) = mfexp(-Z(k,i));
    break;
      //----------------------------------------------------------------------------------------------
    case 0:
      //hindcast of the model
      F(k,i) = fsh_sel(k) * mean_F(k) * mfexp(F_dev(k,i));
      Z(k,i) = F(k,i) + M1(k)(1,nages(k)) + M2(k,i);
      S(k,i) = mfexp(-Z(k,i));
    break;
      //----------------------------------------------------------------------------------------------
    case 1:
      // forecast without fishing
      F_fut(k,i)=fsh_sel(k) *0;
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
      //----------------------------------------------------------------------------------------------
    case 2:
      // forecast under mean F
      F_fut(k,i) = fsh_sel(k) * mean_F(k);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
      //----------------------------------------------------------------------------------------------
    case 3:
      // forecast under control rules
      F_fut(k,i) = fsh_sel(k) * PF(itemp,k,i);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
      //----------------------------------------------------------------------------------------------
    case 4:
      // forecast under NewCatch.dat
      tmpF = PF(itemp,k,i);
      tmpF = SolveF2(N_fut(k,i), NewCatch(k),k);
      PF(itemp,k,i) = tmpF;
      F_fut(k,i) = fsh_sel(k) *tmpF;
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
      //----------------------------------------------------------------------------------------------
    case 5:
      // forecast under mean fishing with random seed
      tmpF = mfexp(sqrt( norm2(F_dev(k)(1,nyrs))/nyrs )*rmultF(k,i)+ln_mean_F(k));
      F_fut(k,i) = fsh_sel(k) *tmpF;
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
      //----------------------------------------------------------------------------------------------
    case 6:
    // forecast under max fishing from past
      tmpF = mfexp(max(F_dev(k)(1,nyrs))+ln_mean_F(k));
      F_fut(k,i) = fsh_sel(k) *tmpF;
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
    //----------------------------------------------------------------------------------------------    
    case 7:
      // forecast under mean catch
      //tmpC=(sqrt( norm2((tc_hat(k)-mean(tc_biom_hat(k))))/nyrs )*rmultF(k,i)+mean(tc_hat(k)));
      tmpC=mean(tc_biom_obs(k));
      tmpCC(1)=(1226280);
      tmpCC(2)=(191938);
      tmpCC(3)=(13458);
      tmpF = PF(itemp,k,i);
      tmpF = SolveF2_kir(N_fut(k,i),(tmpCC(k)),k);
      PF(itemp,k,i) = tmpF;
      F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
    //----------------------------------------------------------------------------------------------
    case 8:
      // forecast under max catch
      tmpC=(max(tc_biom_obs(k)));
      tmpCC(1)=(1490900);
      tmpCC(2)=(220134);
      tmpCC(3)=(17737);
      tmpF = PF(itemp,k,i);
      tmpF = SolveF2_kir(N_fut(k,i),(tmpCC(k)),k);
      PF(itemp,k,i) = tmpF;
      F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
    //----------------------------------------------------------------------------------------------
    case 9:
      // forecast under f profiles
      PF(itemp,k,i) = FP_in(k);
      F_fut(k,i) = fsh_sel(k) * FP_in(k);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
    //----------------------------------------------------------------------------------------------
    case 10:
      // forecast under set catch
      tmpC=mean(tc_biom_obs(k));
      if (msmMode==0)
        tmpCC=tmpCatch_in(1);
      if (msmMode==2)
        tmpCC=tmpCatch_in(2);
      tmpF = PF(itemp,k,i);
      tmpF = SolveF2_kir(N_fut(k,i),(tmpCC(k)),k);
      PF(itemp,k,i) = tmpF;
      F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
    //----------------------------------------------------------------------------------------------
    case 11:
      // forecast under a set catch (if set_catch(sp)==1); else forecast under a set F
      if(set_catch(k)==1.0){
        tmpF = PF(itemp,k,i);
        tmpF = SolveF2(N_fut(k,i), set_val(k),k);  // solve for the Frate that yields the set catch
        PF(itemp,k,i) = tmpF;
        F_fut(k,i) = fsh_sel(k) *tmpF;
        Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
        S_fut(k,i) = mfexp(-Z_fut(k,i));
        M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
        F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
      }else{
         // project under set fishing rate (PF)
        PF(itemp,k,i)=set_val(k);
        if(PF(itemp,k,i)>100){
          cout<<"ERROR - F rate appears to be Catch; check set_catch.dat values for species "<<k<<endl;
          exit(1);
        } 
        F_fut(k,i) = fsh_sel(k) * PF(itemp,k,i);
        Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
        S_fut(k,i) = mfexp(-Z_fut(k,i));
        M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
        F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
      }
    break;
    //----------------------------------------------------------------------------------------------
    // case 12:
    //     // catulcate Fabc from F40 value from set F values in dat file e.g., dat_input_files/setFabcFofl.dat
    //     double Cbb;
    //     Cbb=.2;
    //     if(k==3) Cbb=0.;
    //     B2B40=0.;
    //     //B2B40=biomassSSB_hind2fut(itemp,k,i)/(Btarget*B0_set(k)); 
    //     if(i==1) B2B40=biomassSSB(k,nyrs)/(Btarget*B0_set(k)); 
    //     if(i>1)  B2B40=biomassSSB_hind2fut(itemp,k,i-1)/(Btarget*B0_set(k)); 

    //     Fabc(itemp,k,i)=TIER3_CR(F40(itemp,k,i),B2B40, alpha_ABC,Cbb);
    //     Fofl(itemp,k,i)=TIER3_CR(F35(itemp,k,i),B2B40, alpha_OFL,Cbb);
    //     if(dump_rep) cout<<"F40(itemp,k= "<<k<<",i= "<<i<<") = "<<F40(itemp,k,i)<<" Fabc ="<<Fabc(itemp,k,i)<<endl;
    //     PF(itemp,k,i) = = Fabc(itemp,k,i);
    //     //PF(itemp,k,i) = tmpF;
    //     //F_fut(k,i) = fsh_sel(k) *tmpF;
    //     F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i) ;
    //     Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
    //     S_fut(k,i) = mfexp(-Z_fut(k,i));
    //     M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
    //     F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);     
    //     catch_hat_fut(itemp,k,i)(1,nages(k))        = elem_prod(elem_div(F_fut(k,i),Z_fut(k,i)) ,elem_prod(1.-mfexp(-Z_fut(k,i)) , N_fut(k,i)));
    //     Frate_fut(itemp,k,i)                         = mean(elem_div(F_fut(k,i),fsh_sel(k)));
    //     tc_hat_fut(itemp,k,i)                        = sum(catch_hat_fut(itemp,k,i)(1,nages(k)));
    //     fsh_age_hat_fut(k,i)                         = catch_hat_fut(itemp,k,i)(1,nages(k)) / tc_hat_fut(itemp,k,i);
    //     tc_biom_hat_fut(itemp,k,i)                   = catch_hat_fut(itemp,k,i)(1,nages(k)) * wt_fut(k,i)(1,nages(k));// matrix multiplication A%*%B
        
    //     ABC_fut(itemp,k,i)=tc_biom_hat_fut(itemp,k,i);
    //     ABC_N_fut(itemp,k,i)=tc_hat_fut(itemp,k,i);
    //     ABC_N_age_fut(itemp,k,i)=catch_hat_fut(itemp,k,i);
    //   } 

    // break;
    //----------------------------------------------------------------------------------------------
    case 12:
      // Don't use amanda's function but otherwise the same as 13 to get catch in each year based on sloping HCR ABC
        // F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i);
        // Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
        // S_fut(k,i) = mfexp(-Z_fut(k,i));
        // M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
        // F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
      // TIER 3 happens before this CALC_mort is called (for case 13)

        F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i);
        Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
        S_fut(k,i) = mfexp(-Z_fut(k,i));
        M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
        F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);     
        catch_hat_fut(itemp,k,i)(1,nages(k))        = elem_prod(elem_div(F_fut(k,i),Z_fut(k,i)) ,elem_prod(1.-mfexp(-Z_fut(k,i)) , N_fut(k,i)));
        Frate_fut(itemp,k,i)                         = mean(elem_div(F_fut(k,i),fsh_sel(k)));
        tc_hat_fut(itemp,k,i)                        = sum(catch_hat_fut(itemp,k,i)(1,nages(k)));
        fsh_age_hat_fut(k,i)                         = catch_hat_fut(itemp,k,i)(1,nages(k)) / tc_hat_fut(itemp,k,i);
        tc_biom_hat_fut(itemp,k,i)                   = catch_hat_fut(itemp,k,i)(1,nages(k)) * wt_fut(k,i)(1,nages(k));// matrix multiplication A%*%B
        

    break;  //----------------------------------------------------------------------------------------------
    case 13:
      // use amanda's function to get catch in each year based on sloping HCR ABC
        // F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i);
        // Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
        // S_fut(k,i) = mfexp(-Z_fut(k,i));
        // M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
        // F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
      // TIER 3 happens before this CALC_mort is called (for case 13)

        F_fut(k,i) = fsh_sel(k) *PF(itemp,k,i);
        Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
        S_fut(k,i) = mfexp(-Z_fut(k,i));
        M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
        F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);     
        catch_hat_fut(itemp,k,i)(1,nages(k))        = elem_prod(elem_div(F_fut(k,i),Z_fut(k,i)) ,elem_prod(1.-mfexp(-Z_fut(k,i)) , N_fut(k,i)));
        Frate_fut(itemp,k,i)                         = mean(elem_div(F_fut(k,i),fsh_sel(k)));
        tc_hat_fut(itemp,k,i)                        = sum(catch_hat_fut(itemp,k,i)(1,nages(k)));
        fsh_age_hat_fut(k,i)                         = catch_hat_fut(itemp,k,i)(1,nages(k)) / tc_hat_fut(itemp,k,i);
        tc_biom_hat_fut(itemp,k,i)                   = catch_hat_fut(itemp,k,i)(1,nages(k)) * wt_fut(k,i)(1,nages(k));// matrix multiplication A%*%B
        

    break;
    //----------------------------------------------------------------------------------------------
    case 14:
      // forecast under set fishing rate (PF)
      F_fut(k,i) = fsh_sel(k) * FP_in(k);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
    //----------------------------------------------------------------------------------------------
    case 15:
      // forecast under set F40
      F_fut(k,i) = fsh_sel(k) * F40(itemp,k,i);
      Z_fut(k,i) = F_fut(k,i) + M1(k)(1,nages(k)) + M2_fut(k,i);
      S_fut(k,i) = mfexp(-Z_fut(k,i));
      M2_fut_all(itemp,k,i)(1,nages(k))=M2_fut(k,i);
      F_fut_all(itemp,k,i)(1,nages(k))=F_fut(k,i);
    break;
    //----------------------------------------------------------------------------------------------
    Frate_fut(itemp,k,i)                         = mean(elem_div(F_fut(k,i),fsh_sel(k)));
  }


FUNCTION dvariable SolveF2_kir(const dvar_vector& N_tmp,dvariable TACin, int k)
  // cout<<"In solver for F "<<endl;
  RETURN_ARRAYS_INCREMENT();
  dvariable dd = 10.;
  dvariable cc = TACin;
  dvariable btmp =  sum(elem_prod(elem_prod(fsh_sel(k),N_tmp), wt_fut(k,i))); // temp fishable biomass
  dvariable ftmp;
  ftmp=0;
  if(btmp>0.1) ftmp = 1.2*TACin/btmp;
  dvar_vector Fatmp = ftmp * fsh_sel(k);
  dvar_vector Z_tmp = Fatmp + M1(k)(1,nages(k)) + M2_fut(k,i);
  dvar_vector S_tmp = exp(-Z_tmp);
  // dvariable flim;
  // flim=10*sum(fsh_sel(k)*M1(k)(1,nages(k)))/sum(fsh_sel(k));

  for (int icount=1;icount<=10;icount++)
  {
    // if(btmp>0.1)    ftmp += (TACin-cc) / btmp;
    if(btmp>0.1)     ftmp = ftmp*mfexp((TACin-cc) / btmp);
    if(btmp<=0.1)    ftmp =0;
     if(ftmp>2.)      ftmp=2.;  
    Fatmp = ftmp * fsh_sel(k);
    Z_tmp = Fatmp + M1(k)(1,nages(k)) + M2_fut(k,i);
    S_tmp = mfexp( -Z_tmp );
    cc = sum(elem_prod(wt_fut(k,i), elem_prod(elem_div(Fatmp,  Z_tmp),elem_prod(1.-S_tmp,N_tmp)))); // Catch equation (vectors)
    dd = cc / TACin - 1.;
  }
  RETURN_ARRAYS_DECREMENT();
  return(ftmp); 

 //----------------------------------------------------------------------------------------------

FUNCTION dvariable SolveF2(const dvar_vector& N_tmp, double  TACin, int k)
  // cout<<"In solver for F "<<endl;
  RETURN_ARRAYS_INCREMENT();
  dvariable dd = 10.;
  dvariable cc = TACin;
  dvariable btmp =  sum(elem_prod(elem_prod(fsh_sel(k),N_tmp), wt_fut(k,i))); // temp fishable biomass
  dvariable ftmp;
  ftmp=0;
  if(btmp>0.1) ftmp = 1.2*TACin/btmp;
  dvar_vector Fatmp = ftmp * fsh_sel(k);
  dvar_vector Z_tmp = Fatmp + M1(k)(1,nages(k)) + M2_fut(k,i);
  dvar_vector S_tmp = exp(-Z_tmp);
  // dvariable flim;
  // flim=10*sum(fsh_sel(k)*M1(k)(1,nages(k)))/sum(fsh_sel(k));


  for (int icount=1;icount<=10;icount++)
  {
    // if(btmp>0.1)    ftmp += (TACin-cc) / btmp;
    if(btmp>0.1)    ftmp = ftmp*mfexp((TACin-cc) / btmp);

    if(btmp<=0.1)    ftmp =0;
    if(ftmp>2.)  ftmp=2.;  // set the upper limit for F mortality
    
    Fatmp = ftmp * fsh_sel(k);
    Z_tmp = Fatmp + M1(k)(1,nages(k)) + M2_fut(k,i);
    S_tmp = mfexp( -Z_tmp );
    cc = sum(elem_prod(wt_fut(k,i), elem_prod(elem_div(Fatmp,  Z_tmp),elem_prod(1.-S_tmp,N_tmp)))); // Catch equation (vectors)
    dd = cc / TACin - 1.;
  }
  RETURN_ARRAYS_DECREMENT();
  return(ftmp); 
//----------------------------------------------------------------------------------------------

FUNCTION void CALC_REC_FUT(int recMode_pass)
 //FUNCTION dvar_matrix CALC_REC_FUT(int recMode_pass)
  //RETURN_ARRAYS_INCREMENT();
  //dvar_matrix    R_fut_tmp(1,nspp,1,nyrs_fut);
  dvariable SSB_rs_fut;
  dvariable covars;
  dvar_vector   rmult_a(1,nyrs_fut);                    // random multiplier for random recruitment function
  dvar_vector   rmult_b(1,nyrs_fut);                    // random multiplier for random recruitment function
  dvar_matrix    rmult_covs(1,nyrs_fut,1,800);                    // random multiplier for random recruitment function
  rmult_a=0.;
  rmult_b=0.;
  rmult_covs=0.;
  dvar_vector Ration_scaled(1,nspp);
  dvar_vector Ration_all(1,nspp);
  Ration_scaled=0.;
  int tmpny;
  tmpny=nyrs-1;
  dvar_vector aa_c(1,nspp);
  dvar_vector bb_c(1,nspp);
  aa_c=0.;
  bb_c=0.;

  for (int sp=1;sp<=nspp;sp++)
  { 

    SSB_rs_fut=0.;
    Ration_scaled=0.;
    Ration_all=0.;
      ration2Age_futKEEP(itemp,sp,i)(1,nages(sp))=0.;
    AvgN_futKEEP(itemp,sp,i)=0.;
    RationScaledKeep_fut(itemp,sp,i)=0.;
    SSB_RSkeep_fut(itemp,sp,i)=0.;
    if(i==1){
      SSB_rs_fut = elem_prod(N(sp,nyrs),pmature(sp))*wt(sp,nyrs)(1,nages(sp));  
      Ration_scaled(sp)=((AvgN(sp,nyrs,1)*ration2Age(sp,nyrs,1))-mnRat(sp))/sdRat(sp);
      Ration_all(sp)=AvgN(sp,nyrs,1)*ration2Age(sp,nyrs,1);
      ration2Age_futKEEP(itemp,sp,i)(1,nages(sp))=ration2Age(sp,nyrs);  //ration2Age_fut(pred,i,pred_age)
      AvgN_futKEEP(itemp,sp,i)=AvgN(sp,nyrs,1);
    }else{
      SSB_rs_fut = elem_prod(N_fut(sp,i-1),pmature(sp))*wt_fut(sp,i-1)(1,nages(sp));
      Ration_scaled(sp)=((AvgN_fut(sp,i-1,1)*ration2Age_fut(sp,i-1,1))-mnRat(sp))/sdRat(sp);
      Ration_all(sp)=AvgN_fut(sp,i,1)*ration2Age_fut(sp,i,1);
      ration2Age_futKEEP(itemp,sp,i)(1,nages(sp))=ration2Age_fut(sp,i); // all ages g/ pred
      AvgN_futKEEP(itemp,sp,i)=AvgN_fut(sp,i,1);
    }
    RationScaledKeep_fut(itemp,sp,i)=Ration_scaled(sp);
    RationKeep_fut(itemp,sp,i)=Ration_all(sp);
    SSB_RSkeep_fut(itemp,sp,i)=SSB_rs_fut;
    if(SSB_rs_fut<=1.0){
      //R_fut_tmp(sp,i)=0.;
      R_fut(itemp,sp,i)=0.;
    }else{
    switch (recMode_pass)
    {
      case 0:  //          0 = project under mean rec
       if(rand_rec){  
         //          1 = project under random rec (based on estimates of sd rec)
          rec_sigma(sp)=sum((rec_dev(sp)(1,nyrs) ))/nyrs;

          //R_fut(itemp,sp,i) = mfexp(sqrt(norm2(rec_dev(sp)(1,nyrs))/nyrs)*rmult(sp,i)+ln_mn_rec(sp));
          R_fut(itemp,sp,i) = mfexp(ln_mn_rec(sp)+rec_sigma(sp)*rmult(sp,i));
          //R_fut_tmp(sp,i) = mfexp(ln_mn_rec(sp)+rec_sigma(sp)*rmult(sp,i));
        }else{
          R_fut(itemp,sp,i) = sum(mfexp(ln_mn_rec(sp) + rec_dev(sp)(1,nyrs) ))/nyrs;
          //R_fut_tmp(sp,i) = sum(mfexp(ln_mn_rec(sp) + rec_dev(sp)(1,nyrs) ))/nyrs;

        }
      break; 
      default:
        covars=0.;
        if(BLM(sp)>0||LM(sp)>0)
        {
          if(BLM(sp))
          {
            covars=0.;
            if(rand_rec){
              rmult_a(i)=randn(rmult(sp,i));
              rmult_b(i)=randn(10+rmult(sp,i));
              aa_c(sp)=log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp);
              bb_c(sp)=log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp);
              for (int c=1;c<=ncovs(sp);c++)
              {
                rmult_covs(i,c)=randn(c+110+rmult(sp,i));
                if(cov_type2(sp,c)){
                  // covars+=rs_parm(c)*(rs_cov(c,i)-Ration_scaled(i));

                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                }else{
                  //covars+=rs_parm(c)*rs_cov(c,i);
                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*rs_cov_fut(itemp,c,i);

                }
              }
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+bb_c(sp)*log(SSB_rs_fut)+covars);
              //R_fut_tmp(sp,i)=mfexp(aa_c(sp)+bb_c(sp)*log(SSB_rs_fut)+covars);
            }else{
              // if not random recruitment
              aa_c(sp)=log_aa_c(sp);
              bb_c(sp)=log_bb_c(sp);
              for (int c=1;c<=ncovs(sp);c++)
              {
                if(cov_type2(sp,c)){
                  covars+=rs_parm(sp,c)*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                  // covars+=rs_parm(c)*(rs_cov(c,i)-Ration_scaled(i));
                }else{
                    covars+=rs_parm(sp,c)*rs_cov_fut(itemp,c,i);
                }
              }
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+bb_c(sp)*log(SSB_rs_fut)+covars);
              //R_fut_tmp(sp,i)=mfexp(aa_c(sp)+bb_c(sp)*log(SSB_rs_fut)+covars);
                   //  cout<<"Ration_scaled(sp) "<<Ration_scaled(sp)<<endl;
                   //  cout<<"Ration(sp) "<<(AvgN(sp,nyrs,1)*ration2Age(sp,nyrs,1))<<endl;
                   //  cout<<"mnRat(sp) "<<mnRat(sp)<<endl;
                   //  cout<<"sdRat(sp) "<<sdRat(sp)<<endl;
                   //  cout<<"cov_type2(sp) "<<cov_type2(sp)<<endl;
                   //  cout<<"rs_parm(sp) "<<rs_parm(sp)<<endl;
                   //  cout<<"aa_c(sp) "<<aa_c(sp)<<endl;
                   //  cout<<"bb_c(sp) "<<bb_c(sp)<<endl;
                   //  cout<<"covars "<<covars<<endl;
                   //  cout<<"rs_cov "<<rs_cov<<endl;
                   //  cout<<"SSB_rs_fut "<<SSB_rs_fut<<endl;
                   //  cout<<"logR "<<aa_c(sp)+bb_c(sp)*SSB_rs_fut+covars<<endl;
                    //cout<<"R_fut(itemp,sp,i) "<<R_fut(itemp,sp,i)<<endl;//exit(1);
            }
          }// if BLM
          if(LM(sp))
          {

            covars=0.;
            if(rand_rec){
              rmult_a(i)=randn(rmult(sp,i));
              rmult_b(i)=randn(10+rmult(sp,i));
              aa_c(sp)=log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp);
              bb_c(sp)=log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp);

              for (int c=1;c<=ncovs(sp);c++)
              {
                rmult_covs(i,c)=randn(c+110+rmult(sp,i));
                if(cov_type2(sp,c)){
                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                }else{
                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*rs_cov_fut(itemp,c,i);
                }
              }
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+covars);
              //R_fut_tmp(sp,i)=mfexp(aa_c(sp)+covars);
            }else{
              aa_c(sp)=log_aa_c(sp);
              bb_c(sp)=log_bb_c(sp);
              for (int c=1;c<=ncovs(sp);c++)
              {
                if(cov_type2(sp,c)){
                  covars+=rs_parm(sp,c)*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                }else{
                  covars+=rs_parm(sp,c)*rs_cov_fut(itemp,c,i);
                }
              }
              
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+covars);
              //R_fut_tmp(sp,i)=mfexp(aa_c(sp)+covars);
            }
          }
        }else{
            //cout<<"RUNNING SR MODEL"<<endl;
            if(BevH(sp)){
              if(rand_rec){
                rmult_a(i)=randn(rmult(sp,i));
                rmult_b(i)=randn(10+rmult(sp,i));
                aa_c(sp)=mfexp(log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp));
                bb_c(sp)=mfexp(log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp));
                for (int c=1;c<=ncovs(sp);c++)
                {
                  rmult_covs(i,c)=randn(c+110+rmult(sp,i));
                  if(cov_type2(sp,c)){
                    covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                  }else{
                    covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*rs_cov_fut(itemp,c,i);
                  }
                }
                R_fut(itemp,sp,i)=(SSB_rs_fut/(aa_c(sp)+(bb_c(sp)*SSB_rs_fut)))+covars;
                //R_hat(i)=(SSB(i)/(aa_c+(bb_c*SSB(i))))+covars;
                //R_fut_tmp(sp,i)=SSB_rs_fut/(aa_c(sp)+(bb_c(sp)*SSB_rs_fut));
              }else{      
                 aa_c(sp)=mfexp(log_aa_c(sp));
                 bb_c(sp)=mfexp(log_bb_c(sp));
                 for (int c=1;c<=ncovs(sp);c++)
                  {
                    if(cov_type2(sp,c)){
                      covars+=rs_parm(sp,c)*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                    }else{
                      covars+=rs_parm(sp,c)*rs_cov_fut(itemp,c,i);
                    }
                  }
                 R_fut(itemp,sp,i)=(SSB_rs_fut/(aa_c(sp)+(bb_c(sp)*SSB_rs_fut)))+covars;
                 //R_fut_tmp(sp,i)=SSB_rs_fut/(aa_c(sp)+(bb_c(sp)*SSB_rs_fut));
              }

            }else{  
              covars=0.;
              if(rand_rec){
                rmult_a(i)=randn(rmult(sp,i));
                rmult_b(i)=randn(10+rmult(sp,i));
                aa_c(sp)=mfexp(log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp));
                bb_c(sp)=mfexp(log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp));
               
                for (int c=1;c<=ncovs(sp);c++)
                {
                  rmult_covs(i,c)=randn(c+110+rmult(sp,i));
                  if(cov_type2(sp,c)){
                    covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                  }else{
                    covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*rs_cov_fut(itemp,c,i);
                  }
                }     
                // R_fut(itemp,sp,i)= mfexp(log(aa_c(sp)*SSB_rs_fut )-bb_c(sp)*SSB_rs_fut+covars);
                 // logR_hat(i)= (log(aa_c*SSB(i)) -bb_c*SSB(i)+covars); // may 2017
               R_fut(itemp,sp,i)=  mfexp(aa_c(sp)-bb_c(sp)*SSB_rs_fut+covars+log(SSB_rs_fut));  //a+b*SSB+sum(covs*betas)+log(SSB)  # mueter 2011 formulation
              

                //R_fut_tmp(sp,i)= mfexp(log(aa_c(sp)*SSB_rs_fut )-bb_c(sp)*SSB_rs_fut+covars);
              }else{
                aa_c(sp)=mfexp(log_aa_c(sp));
                bb_c(sp)=mfexp(log_bb_c(sp));
             
                for (int c=1;c<=ncovs(sp);c++)
                {
                  // cout<<"c "<<c<<endl;
                  // cout<<"cov_type2(sp,c) "<<cov_type2(sp,c)<<endl;
                  // cout<<"Ration_scaled(sp) "<<Ration_scaled(sp)<<endl;
                  // cout<<"rs_cov_fut(itemp,c,i) "<<rs_cov_fut(itemp,c,i)<<endl;
                  // cout<<"covars"<<covars<<endl;
                  // cout<<"rs_parm(sp,c)*(Ration_scaled(sp)/rs_cov_fut(itemp,c,i)) "<<endl;
                  // cout<<rs_parm(sp,c)*(Ration_scaled(sp)/rs_cov_fut(itemp,c,i))<<endl;
                    if(cov_type2(sp,c)){
                      covars+=rs_parm(sp,c)*(rs_cov_fut(itemp,c,i)-Ration_scaled(sp));
                    }else{
                      covars+=rs_parm(sp,c)*rs_cov_fut(itemp,c,i);
                    }
                }
                // R_fut(itemp,sp,i)= mfexp(log(aa_c(sp)*SSB_rs_fut )-bb_c(sp)*SSB_rs_fut+covars);
                 R_fut(itemp,sp,i)=  mfexp(aa_c(sp)-bb_c(sp)*SSB_rs_fut+covars+log(SSB_rs_fut));  //a+b*SSB+sum(covs*betas)+log(SSB)  # mueter 2011 formulation
    
              }
            }  
          } 
        break;
      }
    }
  }
  //RETURN_ARRAYS_DECREMENT();
  //return(R_fut_tmp);  
  //----------------------------------------------------------------------------------------------
FUNCTION void CALC_REC_FUTTRASH(int recMode_pass)

  dvariable SSB_rs_fut;
  dvariable covars;
  dvar_vector   rmult_a(1,nyrs_fut);                    // random multiplier for random recruitment function
  dvar_vector   rmult_b(1,nyrs_fut);                    // random multiplier for random recruitment function
  dvar_matrix    rmult_covs(1,nyrs_fut,1,800);                    // random multiplier for random recruitment function
  rmult_a=0.;
  rmult_b=0.;
  rmult_covs=0.;
  dvar_vector Ration_scaled(1,nspp);
  dvar_vector Ration_all(1,nspp);
  Ration_scaled=0.;
  int tmpny;
  tmpny=nyrs-1;
  dvar_vector aa_c(1,nspp);
  dvar_vector bb_c(1,nspp);
  aa_c=0.;
  bb_c=0.;

  for (int sp=1;sp<=nspp;sp++)
  { 

    SSB_rs_fut=0.;
    Ration_scaled=0.;
    Ration_all=0.;
      ration2Age_futKEEP(itemp,sp,i)(1,nages(sp))=0.;
    AvgN_futKEEP(itemp,sp,i)=0.;
    RationScaledKeep_fut(itemp,sp,i)=0.;
    SSB_RSkeep_fut(itemp,sp,i)=0.;
    if(i==1){
      SSB_rs_fut = elem_prod(N(sp,nyrs),pmature(sp))*wt(sp,nyrs)(1,nages(sp));  
      Ration_scaled(sp)=((AvgN(sp,nyrs,1)*ration2Age(sp,nyrs,1))-mnRat(sp))/sdRat(sp);
      Ration_all(sp)=AvgN(sp,nyrs,1)*ration2Age(sp,nyrs,1);
      ration2Age_futKEEP(itemp,sp,i)(1,nages(sp))=ration2Age(sp,nyrs);
      AvgN_futKEEP(itemp,sp,i)=AvgN(sp,nyrs,1);
    }else{
      SSB_rs_fut = elem_prod(N_fut(sp,i-1),pmature(sp))*wt_fut(sp,i-1)(1,nages(sp));
      Ration_scaled(sp)=((AvgN_fut(sp,i-1,1)*ration2Age_fut(sp,i-1,1))-mnRat(sp))/sdRat(sp);
      Ration_all(sp)=AvgN_fut(sp,i,1)*ration2Age_fut(sp,i,1);
      ration2Age_futKEEP(itemp,sp,i)(1,nages(sp))=ration2Age_fut(sp,i);
      AvgN_futKEEP(itemp,sp,i)=AvgN_fut(sp,i,1);


    }
    RationScaledKeep_fut(itemp,sp,i)=Ration_scaled(sp);
    RationKeep_fut(itemp,sp,i)=Ration_all(sp);
    SSB_RSkeep_fut(itemp,sp,i)=SSB_rs_fut;
    if(SSB_rs_fut<=1.0){
      R_fut(itemp,sp,i)=0.;
    }else{
    switch (recMode_pass)
    {
      case 0:  //          0 = project under mean rec
       if(rand_rec){  
         //          1 = project under random rec (based on estimates of sd rec)
          rec_sigma(sp)=sum((rec_dev(sp)(1,nyrs) ))/nyrs;

          //R_fut(itemp,sp,i) = mfexp(sqrt(norm2(rec_dev(sp)(1,nyrs))/nyrs)*rmult(sp,i)+ln_mn_rec(sp));
          R_fut(itemp,sp,i) = mfexp(ln_mn_rec(sp)+rec_sigma(sp)*rmult(sp,i));
        }else{
          R_fut(itemp,sp,i) = sum(mfexp(ln_mn_rec(sp) + rec_dev(sp)(1,nyrs) ))/nyrs;

        }
      break; 
      default:
        covars=0.;
        if(BLM(sp)>0||LM(sp)>0)
        {
          if(BLM(sp))
          {
            covars=0.;
            if(rand_rec){
              rmult_a(i)=randn(rmult(sp,i));
              rmult_b(i)=randn(10+rmult(sp,i));
              aa_c(sp)=log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp);
              bb_c(sp)=log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp);
              for (int c=1;c<=ncovs(sp);c++)
              {
                rmult_covs(i,c)=randn(c+110+rmult(sp,i));
                if(cov_type2(sp,c)){
                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*(Ration_scaled(sp)/rs_cov(c,i));
                }else{
                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*rs_cov(c,i);
                }
              }
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+bb_c(sp)*log(SSB_rs_fut)+covars);
            }else{
              aa_c(sp)=log_aa_c(sp);
              bb_c(sp)=log_bb_c(sp);
              for (int c=1;c<=ncovs(sp);c++)
              {
                if(cov_type2(sp,c)){
                  covars+=rs_parm(sp,c)*(Ration_scaled(sp)/rs_cov(c,i));
                }else{
                    covars+=rs_parm(sp,c)*rs_cov(c,i);
                }
              }
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+bb_c(sp)*log(SSB_rs_fut)+covars);
                   //  cout<<"Ration_scaled(sp) "<<Ration_scaled(sp)<<endl;
                   //  cout<<"Ration(sp) "<<(AvgN(sp,nyrs,1)*ration2Age(sp,nyrs,1))<<endl;
                   //  cout<<"mnRat(sp) "<<mnRat(sp)<<endl;
                   //  cout<<"sdRat(sp) "<<sdRat(sp)<<endl;
                   //  cout<<"cov_type2(sp) "<<cov_type2(sp)<<endl;
                   //  cout<<"rs_parm(sp) "<<rs_parm(sp)<<endl;
                   //  cout<<"aa_c(sp) "<<aa_c(sp)<<endl;
                   //  cout<<"bb_c(sp) "<<bb_c(sp)<<endl;
                   //  cout<<"covars "<<covars<<endl;
                   //  cout<<"rs_cov "<<rs_cov<<endl;
                   //  cout<<"SSB_rs_fut "<<SSB_rs_fut<<endl;
                   //  cout<<"logR "<<aa_c(sp)+bb_c(sp)*SSB_rs_fut+covars<<endl;
                    //cout<<"R_fut(itemp,sp,i) "<<R_fut(itemp,sp,i)<<endl;//exit(1);
            }
          }
          if(LM(sp))
          {

            covars=0.;
            if(rand_rec){
              rmult_a(i)=randn(rmult(sp,i));
              rmult_b(i)=randn(10+rmult(sp,i));
              aa_c(sp)=log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp);
              bb_c(sp)=log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp);

              for (int c=1;c<=ncovs(sp);c++)
              {
                rmult_covs(i,c)=randn(c+110+rmult(sp,i));
                if(cov_type2(sp,c)){
                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*(Ration_scaled(sp)/rs_cov(c,i));
                }else{
                  covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*rs_cov(c,i);
                }
              }
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+covars);
            }else{
              aa_c(sp)=log_aa_c(sp);
              bb_c(sp)=log_bb_c(sp);
              for (int c=1;c<=ncovs(sp);c++)
              {
                if(cov_type2(sp,c)){
                  covars+=rs_parm(sp,c)*(Ration_scaled(sp)/rs_cov(c,i));
                }else{
                  covars+=rs_parm(sp,c)*rs_cov(c,i);
                }
              }
              
              R_fut(itemp,sp,i)=mfexp(aa_c(sp)+covars);
            }
          }
        }else{
            //cout<<"RUNNING SR MODEL"<<endl;
            if(BevH(sp)){
              if(rand_rec){
                rmult_a(i)=randn(rmult(sp,i));
                rmult_b(i)=randn(10+rmult(sp,i));
                aa_c(sp)=mfexp(log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp));
                bb_c(sp)=mfexp(log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp));
                R_fut(itemp,sp,i)=SSB_rs_fut/(aa_c(sp)+(bb_c(sp)*SSB_rs_fut));
              }else{      
                 aa_c(sp)=mfexp(log_aa_c(sp));
                 bb_c(sp)=mfexp(log_bb_c(sp));
                 R_fut(itemp,sp,i)=SSB_rs_fut/(aa_c(sp)+(bb_c(sp)*SSB_rs_fut));
              }
            }else{  
              covars=0.;
              if(rand_rec){
                rmult_a(i)=randn(rmult(sp,i));
                rmult_b(i)=randn(10+rmult(sp,i));
                aa_c(sp)=mfexp(log_aa_c_std(sp)*rmult_a(i)+log_aa_c(sp));
                bb_c(sp)=mfexp(log_bb_c_std(sp)*rmult_b(i)+log_bb_c(sp));
                for (int c=1;c<=ncovs(sp);c++)
                {
                  rmult_covs(i,c)=randn(c+110+rmult(sp,i));
                  if(cov_type2(sp,c)){
                    covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*(Ration_scaled(sp)/rs_cov(c,i));
                  }else{
                    covars+=(rs_parm_std(sp,c)*rmult_covs(i,c)+rs_parm(sp,c) )*rs_cov(c,i);
                  }
                }         
                R_fut(itemp,sp,i)= mfexp(log(aa_c(sp)*SSB_rs_fut )-bb_c(sp)*SSB_rs_fut+covars);
              }else{
                aa_c(sp)=mfexp(log_aa_c(sp));
                bb_c(sp)=mfexp(log_bb_c(sp));
                for (int c=1;c<=ncovs(sp);c++)
                {
                  // cout<<"c "<<c<<endl;
                  // cout<<"cov_type2(sp,c) "<<cov_type2(sp,c)<<endl;
                  // cout<<"Ration_scaled(sp) "<<Ration_scaled(sp)<<endl;
                  // cout<<"rs_cov(c,i) "<<rs_cov(c,i)<<endl;
                  // cout<<"covars"<<covars<<endl;
                  // cout<<"rs_parm(sp,c)*(Ration_scaled(sp)/rs_cov(c,i)) "<<endl;
                  // cout<<rs_parm(sp,c)*(Ration_scaled(sp)/rs_cov(c,i))<<endl;
                  if(cov_type2(sp,c)){
                    covars+=rs_parm(sp,c)*(Ration_scaled(sp)/rs_cov(c,i));
                  }else{
                    covars+=rs_parm(sp,c)*rs_cov(c,i);
                  }
                }
                R_fut(itemp,sp,i)= mfexp(log(aa_c(sp)*SSB_rs_fut )-bb_c(sp)*SSB_rs_fut+covars);
              }
            }  
          } 
        break;
      }
    }
           
  }
//----------------------------------------------------------------------------------------------
  
FUNCTION  CALC_WT_AT_AGE_FUT
  // ofstream Wage_fut_report("results/Wage_fut.dat");
  int Wage_mode=1;
  dvar_vector tmpage(1,maxA);
  for(int tmp_age=1;tmp_age<=maxA;tmp_age++)
      tmpage(tmp_age)=double(tmp_age);
  Wage_fut_report                    <<
    "Wage_Mode"                      <<" "<<
    "Climate_scenario"                             <<" "<<
    "spp"                            <<" "<<
    "yr"                             <<" "<<
    "wt_fut_all "   <<" "<<
    tmpage                           <<endl;
  for (int itemptmp=1;itemptmp<=ntemp_scen;itemptmp++)
  { 
     if (Wage_mode==0)
     {
        for(int spp=1;spp<=nspp;spp++)
        {
          for (int yr=1;yr<=nyrs_fut;yr++)
         {
            wt_fut(spp,yr)=wt(spp,nyrs);
            wt_fut_all(itemptmp,spp,yr)=wt_fut(spp,yr);
         }
       }
     }
    if(Wage_mode==1)
    {
     dvar_matrix sigma(1,nspp,1, mf_type);
     dvar_matrix tmpK(1,nspp,1, mf_type);
     dvar_matrix tmpH(1,nspp,1, mf_type);
     tmpK  =mfexp(logK);
     tmpH  =mfexp(logH);
      for(int spp=1;spp<=nspp;spp++)
      {
        dvar_matrix d(1,mf_type(spp),1,nyrs_fut);
        dvar_matrix Winf(1,mf_type(spp),1,nyrs_fut);
        
       for (int yr=1;yr<=nyrs_fut;yr++)
       {
        dvar_matrix logwt_fut(1,mf_type(spp),1,nages(spp));
        for(int sp_age=1;sp_age<=nages(spp);sp_age++){
            // same values for both sexes
          for(int mf=1;mf<=mf_type(spp);mf++){
            d(mf,yr)             = mfexp(log_mean_d(spp,mf)+Tcoef(spp,mf)*TempC_futUSE(itemptmp,yr));
            Winf(mf,yr)          = pow((tmpH(spp,mf)/tmpK(spp,mf)),1./(1. - d(mf,yr)) );
            logwt_fut(mf,sp_age) =(  log(Winf(mf,yr)) + (1./(1. - d(mf,yr)))*log(1. - mfexp(-tmpK(spp,mf) * (1. - d(mf,yr)) * (double(sp_age) - t0(spp,mf)))) );
          }  // end mf_type
         }// end sp_age   
         if(mf_type(spp)==1)
            wt_fut(spp,yr)=mfexp(logwt_fut(mf_type(spp)))/1000;
         //Take mean of males and females with their own resid. mort rates 
         if(mf_type(spp)==2)
          wt_fut(spp,yr)=elem_div((elem_prod(mfexp(logwt_fut(1)),propMorF(spp,1))+elem_prod(mfexp(logwt_fut(2)),propMorF(spp,2))),(propMorF(spp,1)+propMorF(spp,2)))/1000;
         wt_fut_all(itemptmp,spp,yr)(1,nages(spp))=wt_fut(spp,yr);
        Wage_fut_report           <<
          Wage_mode               <<" "<<
          itemptmp                   <<" "<<
          spp                     <<" "<<
          yr                      <<" "<<
          wt_fut_all(itemptmp,spp,yr)<<" "<<endl;
       }//end for each fut yr
      }//end for each spp
    }// end Wage_mode
    if(Wage_mode==2){
      // calculate weight at age from bioenergetics function --> need mean Pval, and mean 
      // mean ration from data
      // fit pvalue based on weight at age.
      // 
     dvar_matrix sigma(1,nspp,1, mf_type);
     dvar_matrix tmpK(1,nspp,1, mf_type);
     dvar_matrix tmpH(1,nspp,1, mf_type);
     tmpK  =mfexp(logK);
     tmpH  =mfexp(logH);
      for(int spp=1;spp<=nspp;spp++)
      {
        dvar_matrix d(1,mf_type(spp),1,nyrs_fut);
        dvar_matrix Winf(1,mf_type(spp),1,nyrs_fut);
        
       for (int yr=1;yr<=nyrs_fut;yr++)
       {
        dvar_matrix logwt_fut(1,mf_type(spp),1,nages(spp));
        for(int sp_age=1;sp_age<=nages(spp);sp_age++){
            // same values for both sexes
          for(int mf=1;mf<=mf_type(spp);mf++){
            d(mf,yr)             = mfexp(log_mean_d(spp,mf)+Tcoef(spp,mf)*TempC_futUSE(itemptmp,yr));
            Winf(mf,yr)          = pow((tmpH(spp,mf)/tmpK(spp,mf)),1./(1. - d(mf,yr)) );
            logwt_fut(mf,sp_age) =(  log(Winf(mf,yr)) + (1./(1. - d(mf,yr)))*log(1. - mfexp(-tmpK(spp,mf) * (1. - d(mf,yr)) * (double(sp_age) - t0(spp,mf)))) );
          }  // end mf_type
         }// end sp_age   
         if(mf_type(spp)==1)
            wt_fut(spp,yr)=mfexp(logwt_fut(mf_type(spp)))/1000;
         //Take mean of males and females with their own resid. mort rates 
         if(mf_type(spp)==2)
          wt_fut(spp,yr)=elem_div((elem_prod(mfexp(logwt_fut(1)),propMorF(spp,1))+elem_prod(mfexp(logwt_fut(2)),propMorF(spp,2))),(propMorF(spp,1)+propMorF(spp,2)))/1000;
         wt_fut_all(itemptmp,spp,yr)(1,nages(spp))=wt_fut(spp,yr);
        Wage_fut_report           <<
          Wage_mode               <<" "<<
          itemptmp                   <<" "<<
          spp                     <<" "<<
          yr                      <<" "<<
          wt_fut_all(itemptmp,spp,yr)<<" "<<endl;
       }//end for each fut yr
      }//end for each spp
    }// end Wage_mode
  } // end for each simulation itemptmp
  // Wage_fut_report.close();


FUNCTION dvariable TIER3_CR(const dvariable& Flim, dvariable Bratio, double alpha, double Cbeta)

  RETURN_ARRAYS_INCREMENT();

  dvariable maxFabc;
  if(Bratio>1.){
    maxFabc=Flim;
  }else{
    if(alpha<Bratio){
      maxFabc=Flim*((Bratio-alpha)/(1.-alpha));
    }else{
      maxFabc=0.;
    }
  }
  if(Bratio<=Cbeta) maxFabc=0.; 

  RETURN_ARRAYS_DECREMENT();
  return(maxFabc); 

//----------------------------------------------------------------------------------------------
  
FUNCTION RUN_CONTROL_RULES
  // Profile over Control rule multipliers...c_mult
  ofstream future_report("../results/Future_report.rep");
  WRITE_FUT_REPORT(0);
  WRITE_CR_PROJECTION_REPORT(0);
  WRITE_PROJECTION_REPORT(0);

  for (crn=1;crn<=crnum;crn++)
  {
    int cr=c_mult(crn,1);//ie 3 1 - crn=3, sub = 1
    int sub=c_mult(crn,2);

    cout<<"running control rules: "<<cr<<"."<<sub<<endl;
    for(int spp=1;spp<=nspp;spp++) BtargetRep(spp)=Btarget;
    switch (cr)
    {
      case 1:
        CONTROL_RULE_1(sub);
      break;
      case 2:
        CONTROL_RULE_2(sub);
      break;
      case 3:
        CONTROL_RULE_3(sub);
      break;
      case 4:
        CONTROL_RULE_4(sub);
      break;
      case 5:
        CONTROL_RULE_5(sub);
      break;
      default:
        cout<<"ERROR: Control rule ? sub ? : please check ctl file"<<endl;exit(1);
      break;
    }

    dvariable tmpnum;
    tmpnum=c_mult(crn,1)*10+c_mult(crn,2);
    WRITE_CR_PROJECTION_REPORT(tmpnum);
    WRITE_PROJECTION_REPORT(tmpnum);

    ofstream future_report("../results/Future_report.rep");
    WRITE_FUT_REPORT(crn);
    ofstream Fout("../results/Fproxy.dat");
              
          Fout << "# Control rule F proxy "<<endl;
              Fout<<"#cr "<<"sub "<<"msmMode "<<"harvestMode "<<"recMode "<<"scen "<<"species "<<"Btarget "<<"targetF "<<endl;
          for (itemp=1;itemp<=ntemp_scen;itemp++)
              for (k=1;k<=nspp;k++)
                 Fout << cr<<" "<<sub<<" "<<msmMode<<" "<<harvestMode<<" "<<recMode<<" "<<itemp<<" "<<k<<" "<<Btarget<<" "<<PF(itemp,k,1)<<endl;
          Fout.close();  

     if(dump_rep) cout<<"_______________________________________________"<<endl;
     if(dump_rep) cout<<"write FUTURE REPORT for control rule number = "<<cr<<"."<<sub<<endl;
     if(dump_rep) cout<<"_______________________________________________"<<endl;
  }//end crn
  exit(1);

FUNCTION void CONTROL_RULE_1(int sub)
  int CRrepn;
  CRrepn=10;
  int nyrs_avg;
  int skip = 0;
  nyrs_avg=5;
  dvar_vector est(1,nspp);
  dvar_matrix CRtmpF(1,ntemp_scen,1,nspp);
  dvariable adjn;
  int ntemp_scen_tmp;
  int nn=nyrs+nyrs_fut;
  PF.initialize();
  dvar_matrix run_est(1,ntemp_scen,1,nspp);
  dvariable tmpEst; 
  cout<<"NOW running control rule 1: sub= "<<sub<<endl;
  switch (sub)
  {
    default:
    cout<<"ERROR: Control rule 1 sub ? : please check ctl file"<<endl;exit(1);
    break;
    //-------------------------------------------------------------------------------------------------------
    case 1:
      //  1       1    Individual stocks fished to B35%; min (B/B0 - 0.35) - use mean
      CRtmpF=.5;
      adjn=1e-8;
      skip = 0;
      est.initialize();
      biomassSSB0_hind2fut_CR.initialize();
      CRtmpF=.5;
      adjn=1e-8;
      est.initialize();
      PF.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=2;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      biomassSSB0_hind2fut_CR=biomassSSB0_hind2fut;  // set the SSB0_CR to the SSB0 of all three fished
      
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
        // if(skip==0){
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=nspp;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
            
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
          {
            for (int sp=1;sp<=nspp;sp++)
            {
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              dvariable diff_val;
              // diff_val=est(sp)-Btarget;
              diff_val=est(sp)-Btarget;
              if (fabs(diff_val)<cutoff){
                // cout<<"difference ("<<fabs(diff_val)<<") is less than "<<cutoff<<endl;
                // skip=1;
              }
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn;  
            }// for each spp
          }// end for each simulation
        // }else{
        //   //skip and do nothing
        // }
      }// for each Rep
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      //biomassSSB0_hind2fut=biomassSSB0_hind2fut_CR; // set unfished biomass for reporting
    break;
    //-------------------------------------------------------------------------------------------------------
    case 2:
      //  1 2 Individual stocks fished to Btarget; constrained such that no B0(spp) <.35 B0(spp) 
      CRtmpF=.5;
      adjn=1e-8;
      est.initialize();
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=nspp;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
          {
            for (int sp=1;sp<=nspp;sp++)
            {
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              dvariable diff_val;
              dvariable minval;
              minval=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_fut,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_fut,(nyrs_fut+nyrs))));
              if(minval<0.35)
              {
                diff_val=minval-0.35;
                // If the  min is less than 0.35
                if(CRtmpF(lsim,sp)>0)
                  CRtmpF(lsim,sp)=mfexp(diff_val/0.35)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/0.35)*CRtmpF(lsim,sp);
                else
                  CRtmpF(lsim,sp)=adjn;   
              }else{
                diff_val=est(sp)-Btarget;
                if(CRtmpF(lsim,sp)>0)
                  CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
                else
                  CRtmpF(lsim,sp)=adjn;  
              }
            }// for each spp
          }// end for each simulation
      }// for each Rep
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      //biomassSSB0_hind2fut=biomassSSB0_hind2fut_CR; // set unfished biomass for reporting
    break;
    //-------------------------------------------------------------------------------------------------------
    case 3:
      //  1 3 Individual stocks fished to Btarget%;pcod and atf to Btarget first, then pollock
      biomassSSB0_hind2fut_CR.initialize();
      CRtmpF=.5;
      adjn=1e-8;
      est.initialize();
      PF.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=2;sp<=nspp;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      biomassSSB0_hind2fut_CR=biomassSSB0_hind2fut;  // set the SSB0_CR to the SSB0 of all three fished
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=2;sp<=nspp;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=2;sp<=nspp;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              //est=.5-.4
              dvariable diff_val;
              diff_val=est(sp)-Btarget;
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn;  
          }// for each spp
        }// end for each simulation
      }// for each Rep
      cout<<"Now run CR for plk "<<endl;   
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        biomassSSB0_hind2fut_CR(lsim,1)=biomassSSB_hind2fut(lsim,1); // set unfished biomass for plk to include pred release from pcod & atf   
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=1;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=1;sp<=1;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              dvariable diff_val;
              diff_val=est(sp)-Btarget;
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn;  
          }// for each spp
        }// end for each simulation
      }// for each Rep
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      //biomassSSB0_hind2fut=biomassSSB0_hind2fut_CR; // set unfished biomass for reporting
    break;
    //-------------------------------------------------------------------------------------------------------
    case 4:
      // 1 4 Individual stocks fished to Btarget;pcod to Btarget first, then atf, then pollock
     //pcod
      biomassSSB0_hind2fut_CR.initialize();
      CRtmpF=.5;
      adjn=1e-8;
      est.initialize();
      PF.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=2;sp<=2;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      biomassSSB0_hind2fut_CR=biomassSSB0_hind2fut;  // set the SSB0_CR to the SSB0 of all three fished
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=2;sp<=2;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=2;sp<=2;sp++)
            {
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              //est=.5-.4
              dvariable diff_val;
              diff_val=est(sp)-Btarget;
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn; 
          }// for each spp
        }// end for each simulation
      }// for each Rep
      cout<<"Now run CR for atf "<<endl;
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
      {
        biomassSSB0_hind2fut_CR(lsim,3)=biomassSSB_hind2fut(lsim,3); // set unfished biomass for atf to include release from pcod
      }  
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=3;sp<=3;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=3;sp<=3;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              //est=.5-.4
              dvariable diff_val;
              diff_val=est(sp)-Btarget;
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn; 
          }// for each spp
        }// end for each simulation
      }// for each Rep
      cout<<"Now run CR for plk "<<endl;
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        biomassSSB0_hind2fut_CR(lsim,1)=biomassSSB_hind2fut(lsim,1); // set unfished biomass for plk to include release from pcod & atf
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=1;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=1;sp<=1;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              //est=.5-.4
              dvariable diff_val;
              diff_val=est(sp)-Btarget;
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn; 
          }// for each spp
        }// end for each simulation
      }// for each Rep
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      //biomassSSB0_hind2fut=biomassSSB0_hind2fut_CR; // set unfished biomass for reporting
    break; 
    //-------------------------------------------------------------------------------------------------------
    case 5:
      //  1 3 Individual stocks fished to Btarget%;pcod and atf to Btarget first, then pollock
      biomassSSB0_hind2fut_CR.initialize();
      CRtmpF=.5;
      adjn=1e-8;
      est.initialize();
      PF.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=2;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      biomassSSB0_hind2fut_CR=biomassSSB0_hind2fut;  // set the SSB0_CR to the SSB0 of all three fished
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=2;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=1;sp<=2;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              //est=.5-.4
              dvariable diff_val;
              diff_val=est(sp)-Btarget;
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn;  
          }// for each spp
        }// end for each simulation
      }// for each Rep
      cout<<"Now run CR for atf "<<endl;   
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        biomassSSB0_hind2fut_CR(lsim,3)=biomassSSB_hind2fut(lsim,3); // set unfished biomass for atf to include release from pcod & plk   
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=3;sp<=3;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=3;sp<=3;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              dvariable diff_val;
              diff_val=est(sp)-Btarget;
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn;  
          }// for each spp
        }// end for each simulation
      }// for each Rep
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      //biomassSSB0_hind2fut=biomassSSB0_hind2fut_CR; // set unfished biomass for reporting
    break;
    //-------------------------------------------------------------------------------------------------------
    case 6:
      //   //  1 4 Individual stocks fished to Btarget;pcod to Btarget first, then pollock, then atf
      //such that no B B35%
      biomassSSB0_hind2fut_CR.initialize();
      CRtmpF=0.5;
      adjn=1e-8;
      est.initialize();
      PF.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      biomassSSB0_hind2fut_CR=biomassSSB0_hind2fut;  // set the SSB0_CR to the SSB0 of all three fished
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        PF(lsim,3)=0.5;
      cout<<"Now run CR for atf "<<endl;
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
        harvestMode=3;
        cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=3;sp<=3;sp++)
            PF(lsim,sp)=CRtmpF(lsim,sp);
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++){
          for (int sp=3;sp<=3;sp++){
            est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
            dvariable diff_val;
            diff_val=est(sp)-Btarget;
            if(CRtmpF(lsim,sp)>0)
              CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);
            else
              CRtmpF(lsim,sp)=adjn;  
          }// for each spp
        }// end for each simulation
      }// for each Rep
      
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=0.5;
      cout<<"Now run CR for atf "<<endl;
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
      {
        biomassSSB0_hind2fut_CR(lsim,1)=biomassSSB0_hind2fut(lsim,1); // set unfished biomass for plk to include release from atf
        biomassSSB0_hind2fut_CR(lsim,2)=biomassSSB0_hind2fut(lsim,2); // set unfished biomass for cpdo to include release from atf
      } 
      cout<<"Now for all three with B0 conditioned on ATF at 0"<<endl;
      // set unfished biomass for plk to include release from pcod & atf
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
        harvestMode=3;
        cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp)=CRtmpF(lsim,sp);
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++){
          for (int sp=1;sp<=nspp;sp++){
            est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
            dvariable diff_val;
            diff_val=est(sp)-Btarget;
            if(CRtmpF(lsim,sp)>0)
              CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
            else
              CRtmpF(lsim,sp)=adjn;  
            
          }// for each spp
        }// end for each simulation
      }// for each Rep
      
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      //biomassSSB0_hind2fut=biomassSSB0_hind2fut_CR; // set unfished biomass for reporting
    break; 
    //-------------------------------------------------------------------------------------------------------
    case 7:
      //  1 4 Individual stocks fished to Btarget;atf first (to get B0 for pcod pollock), such that no B B35%
      biomassSSB0_hind2fut_CR.initialize();
      CRtmpF=0.5;
      adjn=1e-8;
      est.initialize();
      PF.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=nspp;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      biomassSSB0_hind2fut_CR=biomassSSB0_hind2fut;  // set the SSB0_CR to the SSB0 of all three fished
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        PF(lsim,3)=0.5;
      cout<<"Now run CR for atf "<<endl;
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=3;sp<=3;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=3;sp<=3;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              //est=.5-.4
              dvariable diff_val;
              dvariable minval;
              minval=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_fut,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-nyrs_fut,(nyrs_fut+nyrs))));
              if(minval<0.35)
              {
                diff_val=minval-0.35;
                // If the  min is less than 0.35
                if(CRtmpF(lsim,sp)>0)
                  CRtmpF(lsim,sp)=mfexp(diff_val/0.35)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/0.35)*CRtmpF(lsim,sp);
                else
                  CRtmpF(lsim,sp)=adjn;   
              }else{
                diff_val=est(sp)-Btarget;
                if(CRtmpF(lsim,sp)>0)
                  CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
                else
                  CRtmpF(lsim,sp)=adjn;  
              }
          }// for each spp
        }// end for each simulation
      }// for each Rep

      for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=nspp;sp++)
              PF(lsim,sp)=0.5;
       cout<<"Now run CR for atf "<<endl;
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
      {
        biomassSSB0_hind2fut_CR(lsim,1)=biomassSSB0_hind2fut(lsim,1); // set unfished biomass for plk to include release from atf
        biomassSSB0_hind2fut_CR(lsim,2)=biomassSSB0_hind2fut(lsim,2); // set unfished biomass for cpdo to include release from atf
      } 
      cout<<"Now for all three with B0 conditioned on ATF at 0"<<endl;
      // for (int lsim=1;lsim<=ntemp_scen;lsim++)
        // biomassSSB0_hind2fut_CR(lsim,1)=biomassSSB_hind2fut(lsim,1); // set unfished biomass for plk to include release from pcod & atf
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=nspp;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++){
            for (int sp=1;sp<=nspp;sp++){
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              //est=.5-.4
              dvariable diff_val;
              dvariable minval;
              minval=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_fut,(nyrs_fut+nyrs)),biomassSSB0_hind2fut_CR(lsim,sp)((nyrs_fut+nyrs)-nyrs_fut,(nyrs_fut+nyrs))));
              
              if(minval<0.35)
              {
                diff_val=minval-0.35;
                // If the  min is less than 0.35
                if(CRtmpF(lsim,sp)>0)
                  CRtmpF(lsim,sp)=mfexp(diff_val/0.35)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/0.35)*CRtmpF(lsim,sp);
                else
                  CRtmpF(lsim,sp)=adjn;   
              }else{
                diff_val=est(sp)-Btarget;
                if(CRtmpF(lsim,sp)>0)
                  CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
                else
                  CRtmpF(lsim,sp)=adjn;  
              }
          }// for each spp
        }// end for each simulation
      }// for each Rep

      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      //biomassSSB0_hind2fut=biomassSSB0_hind2fut_CR; // set unfished biomass for reporting
    break; 
    //-------------------------------------------------------------------------------------------------------
    case 8:
      // use set B0 for estimating B proxy but hold ATF at mean historical F
      // Individual stocks fished to B35%; min (B/B0 - 0.35) - use mean
      CRtmpF=.5;
      adjn=1e-8;
      skip = 0;
      est.initialize();
      for (int CRrep=1;CRrep<=(CRrepn);CRrep++)
      {         
        // if(skip==0){
          harvestMode=3;
          cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            for (int sp=1;sp<=2;sp++)
              PF(lsim,sp)=CRtmpF(lsim,sp);
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
            PF(lsim,3)=mean_F(3);  // SET ATF to mean historical
          
          RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
          for (int lsim=1;lsim<=ntemp_scen;lsim++)
          {
            for (int sp=1;sp<=2;sp++)
            {
              //est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
              est(sp)=mean(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)))/(B0_set(sp));

              dvariable diff_val;
              // diff_val=est(sp)-Btarget;
              diff_val=est(sp)-Btarget;
              if (fabs(diff_val)<cutoff){
                // cout<<"difference ("<<fabs(diff_val)<<") is less than "<<cutoff<<endl;
                // skip=1;
              }
              if(CRtmpF(lsim,sp)>0)
                CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
              else
                CRtmpF(lsim,sp)=adjn;  
            }// for each spp
          }// end for each simulation
      }// for each Rep
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=2;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        PF(lsim,3)=mean_F(3);  // SET ATF to mean historical F
      //Now get F40 for ATF:
      for (int CRrep=1;CRrep<=(CRrepn);CRrep++)
      {         
        // if(skip==0){
        cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=3;sp<=3;sp++)
            PF(lsim,sp)=CRtmpF(lsim,sp);
       
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=3;sp<=3;sp++)
          {
            //est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
            est(sp)=mean(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)))/(B0_set(sp));
            
            dvariable diff_val;
            // diff_val=est(sp)-Btarget;
            diff_val=est(sp)-Btarget;
            if (fabs(diff_val)<cutoff){
              // cout<<"difference ("<<fabs(diff_val)<<") is less than "<<cutoff<<endl;
              // skip=1;
            }
            if(CRtmpF(lsim,sp)>0)
              CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
            else
              CRtmpF(lsim,sp)=adjn;  
          }// for each spp
        }// end for each simulation
      }// for each Rep
     
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()

        cout<<"mean_F= "<<mean_F<<"mfexp(ln_mean_F)= "<<mfexp(ln_mean_F)<<endl;
        cout<<"mfexp(mean(F_dev(1)+ln_mean_F(1)))= "<<exp(mean(F_dev(1)+ln_mean_F(1)))<<endl;
        cout<<""<<endl;
        cout<<"***** Target F for scenario 1  = "<<PF(1,1,1)<<" "<<PF(1,2,1)<<" "<<PF(1,3,1)<<" *****"<<endl;
        cout<<""<<endl;
    break; 
    //-------------------------------------------------------------------------------------------------------
    case 9:
      // use set B0 for estimating B proxy
      CRtmpF=.5;
      adjn=1e-8;
      skip = 0;
      est.initialize();
      for (int CRrep=1;CRrep<=(CRrepn);CRrep++)
      {         
        // if(skip==0){
        harvestMode=3;
        cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp)=CRtmpF(lsim,sp);
        
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            //est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs))));
            est(sp)=mean(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-5,(nyrs_fut+nyrs)))/B0_set(sp);
            dvariable diff_val;
            diff_val=est(sp)-Btarget;
            if (fabs(diff_val)<cutoff){
              // cout<<"difference ("<<fabs(diff_val)<<") is less than "<<cutoff<<endl;
              // skip=1;
            }
            if(CRtmpF(lsim,sp)>0)
              CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
            else
              CRtmpF(lsim,sp)=adjn;  
          }// for each spp
        }// end for each simulation
      }// for each Rep
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      
      cout<<""<<endl;
      cout<<"***** Target F for scenario 1  = "<<PF(1,1,1)<<" "<<PF(1,2,1)<<" "<<PF(1,3,1)<<" *****"<<endl;
      cout<<""<<endl;
    break;
    //-------------------------------------------------------------------------------------------------------
  }//end switch
  

FUNCTION void CONTROL_RULE_2(int sub)
  int CRrepn;
  CRrepn=15;
  int nyrs_avg;
  nyrs_avg=5;
  dvar_vector est(1,nspp);
  dvar_matrix CRtmpF(1,ntemp_scen,1,nspp);
  dvar_matrix SPR_0(1,nspp,1,nages);
  dvar_matrix SPR_F(1,nspp,1,nages);
  dvar_matrix mnM2(1,nspp,1,nages);
  dvar_matrix mnWt(1,nspp,1,nages);
  dvar_matrix NPR(1,nspp,1,nages);
  dvariable adjn;
  int ntemp_scen_tmp;
  int nn=nyrs+nyrs_fut;
  PF.initialize();
  dvar_matrix run_est(1,ntemp_scen,1,nspp);
  dvariable tmpEst; 
  cout<<"NOW running control rule 1: sub= "<<sub<<endl;
  switch (sub)
  {
    default:
    cout<<"ERROR: Control rule 2 sub ? : please check ctl file"<<endl;exit(1);
    break;
    case 1:
      //  2       1    use YPR Individual stocks fished to B35%; min (B/B0 - 0.35) - use mean
      CRtmpF=.5;
      adjn=1e-8;
       harvestMode=3;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp)=0;
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          est.initialize();
          NPR.initialize();
          mnM2.initialize();
          mnWt.initialize();
          SPR_0.initialize();
          SPR_F.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
               for (int sp_age=1;sp_age<=nages(sp);sp_age++)
               {
                  for (int yruse=(nyrs_fut-nyrs_avg);yruse<=nyrs_fut;yruse++)
                  {
                    mnM2(sp,sp_age)+=M2_fut_all(lsim,sp,yruse,sp_age);
                    mnWt(sp,sp_age)+=wt_fut_all(lsim,sp,yruse,sp_age);
                  }  
                  mnM2(sp,sp_age)=mnM2(sp,sp_age)/nyrs_avg;
                  mnWt(sp,sp_age)=mnWt(sp,sp_age)/nyrs_avg;
               }
          
            // now calculate SSB_perR_0 as first year projection:
              
              for (int sp_age=1;sp_age<=nages(sp);sp_age++)
              {
                if(sp_age==1)
                  NPR(sp,sp_age)=1;
                else
                  NPR(sp,sp_age)=NPR(sp,sp_age-1)*mfexp(-0*fsh_sel(sp,sp_age)-M1(sp,sp_age)-mnM2(sp,sp_age));
                SPR_0(sp,sp_age)=pmature(sp,sp_age)*mnWt(sp,sp_age)* NPR(sp,sp_age);
              }      
              CRtmpF(lsim,sp)=.5;            
              for (int CRrep=1;CRrep<=CRrepn;CRrep++)
              {
                cout<<"CR iteration = "<<CRrep<<"  est =" <<est<<endl;
                for (int sp_age=1;sp_age<=nages(sp);sp_age++)
                {
                  if(sp_age==1)
                    NPR(sp,sp_age)=1;
                  else
                    NPR(sp,sp_age)=NPR(sp,sp_age-1)*mfexp(-CRtmpF(lsim,sp)*fsh_sel(sp,sp_age)-M1(sp,sp_age)-mnM2(sp,sp_age));
                  SPR_F(sp,sp_age)=pmature(sp,sp_age)*mnWt(sp,sp_age)* NPR(sp,sp_age);
                }  
                est(sp)=sum(SPR_F(sp))/sum(SPR_0(sp));
                dvariable diff_val;
                diff_val=est(sp)-Btarget;
                if(CRtmpF(lsim,sp)>0)
                  CRtmpF(lsim,sp)=mfexp(diff_val/Btarget)*CRtmpF(lsim,sp);//CRtmpF(lsim,sp)+=(diff_val/Btarget)*CRtmpF(lsim,sp);
                else
                  CRtmpF(lsim,sp)=adjn;  
              } // end for reps
              PF(lsim,sp)=CRtmpF(lsim,sp);
        }// end for each species
    }// end for temp simulations
    RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
    break;   
  }//end switch
FUNCTION void CONTROL_RULE_3(int sub)
  int CRrepn;
  CRrepn=10;
  int nyrs_avg;
  nyrs_avg=5;
  dvar_vector est(1,nspp);
  dvariable totB;
  dvariable totB0;
  dvar_vector minB(1,nspp);
  dvar_matrix CRtmpF(1,ntemp_scen,1,nspp);
  dvar_vector CRtmpF_all(1,ntemp_scen);
  dvariable adjn;
  int ntemp_scen_tmp;
  int nn=nyrs+nyrs_fut;
  PF.initialize();
  dvar_matrix run_est(1,ntemp_scen,1,nspp);
  dvariable tmpEst; 
  dvar_vector propPF(1,nspp);
    
  cout<<"NOW running control rule 3: sub= "<<sub<<endl;

  switch (sub)
  {
    default:
    cout<<"ERROR: Control rule 3 sub ? : please check ctl file"<<endl;exit(1);
    break;
    case 1:
   // 3 1 #sum  F Biomass = 0.35% sum (B0); scaled  given average M's
  //  1 2 Individual stocks fished to Btarget; constrained such that no B0(spp) <.35 B0(spp) 
      CRtmpF=0.5;
      CRtmpF_all=0.5;
      adjn=1e-8;
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      propPF.initialize();
      for (int sp=1;sp<=nspp;sp++)
      {
        for (int yr=1;yr<=nyrs;yr++)
          propPF(sp)+=mean(M1(sp)(1,nages(sp)) + M2(sp,yr));
        propPF(sp)=propPF(sp)/nyrs;
      }  
      propPF=propPF/sum(propPF);
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
        harvestMode=3;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp)=CRtmpF_all(lsim)*propPF(sp);
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          est.initialize();
          totB.initialize();
          totB0.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs))));
              totB+=mean(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
              totB0+=mean(biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
          }  // for each spp
          dvariable diff_val;
          diff_val=(totB/totB0)-Btarget;
          if(CRtmpF_all(lsim)>0)
            CRtmpF_all(lsim)=mfexp(diff_val/Btarget)*CRtmpF_all(lsim);//CRtmpF_all(lsim)+=(diff_val/Btarget)*CRtmpF_all(lsim);
          else
            CRtmpF_all(lsim)=adjn; 
        }// end for each simulation
        cout<<"CR iteration = "<<CRrep<<"  Ball =" <<(totB/totB0)<<endl;
      }// for each Rep
     break; 
    case 2:
      //3 2 #sum  F Biomass = 0.35% sum (B0); scaled  given average M's,  constrained such  that  no  B0(spp) <.35  B0(spp)
      CRtmpF=0.5;
      CRtmpF_all=0.5;
      adjn=1e-8;
      
      
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      propPF.initialize();
      for (int sp=1;sp<=nspp;sp++){
        for (int yr=1;yr<=nyrs;yr++)
          propPF(sp)+=mean(M1(sp)(1,nages(sp)) + M2(sp,yr));
        propPF(sp)=propPF(sp)/nyrs;
      }  
      propPF=propPF/sum(propPF);
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
        harvestMode=3;
        cout<<"CR iteration = "<<CRrep<<"  Ball =" <<(totB/totB0)<<endl;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp)=CRtmpF_all(lsim)*propPF(sp);
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          est.initialize();
          totB.initialize();
          totB0.initialize();
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs))));
              totB+=mean(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
              totB0+=mean(biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
              minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          }  // for each spp
          dvariable diff_val;
          if(min(minB)<0.35)
          {
            diff_val=min(minB)-0.35;
            if(CRtmpF_all(lsim)>0)
              CRtmpF_all(lsim)=mfexp(diff_val/Btarget)*CRtmpF_all(lsim);//CRtmpF_all(lsim)+=(diff_val/Btarget)*CRtmpF_all(lsim);
            else
              CRtmpF_all(lsim)=adjn; 
          }else{ 
            diff_val=(totB/totB0)-Btarget;
            if(CRtmpF_all(lsim)>0)
              CRtmpF_all(lsim)=mfexp(diff_val/Btarget)*CRtmpF_all(lsim);//CRtmpF_all(lsim)+=(diff_val/Btarget)*CRtmpF_all(lsim);
            else
              CRtmpF_all(lsim)=adjn; 
          }
        }// end for each simulation
      }// for each Rep
    break;
    case 3:
      //3 3 #sum  F Biomass = 0.35% sum (B0); scaled  given B0s
      CRtmpF=0.5;
      CRtmpF_all=0.5;
      adjn=1e-8;
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
        harvestMode=3;
        for (int lsim=1;lsim<=ntemp_scen;lsim++){
          for (int sp=1;sp<=nspp;sp++)
            propPF(sp)=mean(biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
          propPF=propPF/sum(propPF);
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp)=CRtmpF_all(lsim)*propPF(sp);
        }  
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          est.initialize();
          totB.initialize();
          totB0.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs))));
              totB+=mean(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
              totB0+=mean(biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
          }  // for each spp
          dvariable diff_val;
          diff_val=(totB/totB0)-Btarget;
          if(CRtmpF_all(lsim)>0)
            CRtmpF_all(lsim)=mfexp(diff_val/Btarget)*CRtmpF_all(lsim);//CRtmpF_all(lsim)+=(diff_val/Btarget)*CRtmpF_all(lsim);
          else
            CRtmpF_all(lsim)=adjn; 
        }// end for each simulation
         cout<<"CR iteration = "<<CRrep<<"  Ball =" <<(totB/totB0)<<endl;
      }// for each Rep
    break;
    case 4:
       //3 4 #sum  F Biomass = 0.35% sum (B0); scaled  given B0, constrained such  that  no  B0(spp) <.35  B0(spp)
    
      CRtmpF=0.5;
      CRtmpF_all=0.5;
      adjn=1e-8;
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
        for (int sp=1;sp<=nspp;sp++)
          PF(lsim,sp)=CRtmpF(lsim,sp);
      RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
      
      for (int CRrep=1;CRrep<=CRrepn;CRrep++)
      {         
        harvestMode=3;
        for (int lsim=1;lsim<=ntemp_scen;lsim++){
          for (int sp=1;sp<=nspp;sp++)
            propPF(sp)=mean(biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
          propPF=propPF/sum(propPF);
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp)=CRtmpF_all(lsim)*propPF(sp);
        }  
        RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          est.initialize();
          totB.initialize();
          totB0.initialize();
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
              est(sp)=mean(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs))));
              totB+=mean(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
              totB0+=mean(biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs_avg,(nyrs_fut+nyrs)));
              minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          }  // for each spp
          dvariable diff_val;
          if(min(minB)<0.35)
          {
            diff_val=min(minB)-0.35;
            if(CRtmpF_all(lsim)>0)
              CRtmpF_all(lsim)=mfexp(diff_val/Btarget)*CRtmpF_all(lsim);//CRtmpF_all(lsim)+=(diff_val/Btarget)*CRtmpF_all(lsim);
            else
              CRtmpF_all(lsim)=adjn; 
          }else{ 
            diff_val=(totB/totB0)-Btarget;
            if(CRtmpF_all(lsim)>0)
              CRtmpF_all(lsim)=mfexp(diff_val/Btarget)*CRtmpF_all(lsim);//CRtmpF_all(lsim)+=(diff_val/Btarget)*CRtmpF_all(lsim);
            else
              CRtmpF_all(lsim)=adjn; 
          }
        }// end for each simulation
      }// for each Rep
   break;

  }//end switch

FUNCTION void CONTROL_RULE_4(int sub)
  // max sustainable yield
  int CRrepn;

  CRrepn=5;
    cout<<"CRrepn = "<<CRrepn<<endl;
  int nyrs_avg;
  nyrs_avg=5;
  dvector est(1,nspp);
  dvar_vector minB(1,nspp);
  dvector maxY(1,nspp);
  // dvar_matrix CRtmpF(1,ntemp_scen,1,nspp);
  dvariable adjn;
  int ntemp_scen_tmp;
  int nn=nyrs+nyrs_fut;
  PF.initialize();
  dmatrix run_est(1,ntemp_scen,1,nspp);
  dvariable tmpEst; 
  dvariable tmpmin; 

  cout<<"NOW running control rule 4: sub= "<<sub<<endl;
  
  dvar_matrix F1(1,ntemp_scen,1,nspp);
  dvar_matrix F2(1,ntemp_scen,1,nspp);
  dvar_matrix F3(1,ntemp_scen,1,nspp);
  dvar_matrix yld1(1,ntemp_scen,1,nspp);
  dvar_matrix yld2(1,ntemp_scen,1,nspp);
  dvar_matrix yld3(1,ntemp_scen,1,nspp);
  dvar_vector yld1_all(1,ntemp_scen);
  dvar_vector yld2_all(1,ntemp_scen);
  dvar_vector yld3_all(1,ntemp_scen);
  dvar_matrix dyldp(1,ntemp_scen,1,nspp);
  dvar_matrix dyld(1,ntemp_scen,1,nspp);
  dvar_vector dyldp_all(1,ntemp_scen);
  dvar_vector dyld_all(1,ntemp_scen);
  imatrix breakout(1,ntemp_scen,1,nspp);
  double df=1.e-05;

  switch (sub)
  {
    default:
    cout<<"ERROR: Control rule 4 sub ? : please check ctl file"<<endl;exit(1);
    break;
    case 1:
      //  Find max yeild of each individual stock 
      F1.initialize();
      breakout.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            F1(lsim,sp) = (0.8*(mean(M1(sp)(1,nages(sp)) + M2(sp,nyrs))));
      breakout=0;
      // Newton Raphson stuff to go here
      for (int ii=1;ii<=8;ii++)
      {
        F2.initialize();
        F3.initialize();
        yld1.initialize();
        yld2.initialize();
        yld3.initialize();
        dyld.initialize();
        dyldp.initialize();
        //if (mceval_phase()&&(F1>5||F1<0.01)) 
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if ((F1(lsim,sp)>5||F1(lsim,sp)<0.01)) 
            {
              ii=8;
              if (F1(lsim,sp)>5) 
                F1(lsim,sp)=5.0; 
              else      
                F1(lsim,sp)=0.001; 
              breakout(lsim,sp)= 1;
            }
          }
        }
        F2     = F1 + df*.5;
        F3     = F2 - df;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld1(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        // run projections with F2
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F2(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld2(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        // run projections with F3
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F3(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld3(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        //yld3   = yield(Fratio,F3);
        dyld   = (yld2 - yld3)/df;                          // First derivative (to find the root of this)
        dyldp  = (yld2 + yld3 - 2.*yld1)/(.25*df*df);       // Second derivative (for Newton Raphson)
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if (breakout(lsim,sp)==0)
            {
              F1(lsim,sp)    -= dyld(lsim,sp)/dyldp(lsim,sp);
            }
            else
            {
              if (F1(lsim,sp)>5) 
                cout<<"Fmsy v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
              else      
                cout<<"Fmsy v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
            }
          }
        }
      }// end for each iteration in 1:8
       for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
    break;
    case 2:
      //  Find max yeild unconditional
      F1.initialize();
      breakout.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            F1(lsim,sp) = (0.8*(mean(M1(sp)(1,nages(sp)) + M2(sp,nyrs))));
      breakout=0;
      // Newton Raphson stuff to go here
      for (int ii=1;ii<=8;ii++)
      {
        F2.initialize();
        F3.initialize();
        yld1.initialize();
        yld2.initialize();
        yld3.initialize();
        dyld_all.initialize();
        dyldp_all.initialize();
        //if (mceval_phase()&&(F1>5||F1<0.01)) 
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if ((F1(lsim,sp)>5||F1(lsim,sp)<0.01)) 
            {
              ii=8;
              if (F1(lsim,sp)>5) 
                F1(lsim,sp)=5.0; 
              else      
                F1(lsim,sp)=0.001; 
              breakout(lsim,sp)= 1;
            }
          }
        }
        F2     = F1 + df*.5;
        F3     = F2 - df;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld1(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        for (int lsim=1;lsim<=ntemp_scen;lsim++)  yld1_all(lsim)=sum(yld1(lsim));
        // run projections with F2
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F2(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld2(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        for (int lsim=1;lsim<=ntemp_scen;lsim++)  yld2_all(lsim)=sum(yld2(lsim));
        // run projections with F3
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F3(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld3(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        for (int lsim=1;lsim<=ntemp_scen;lsim++)  yld3_all(lsim)=sum(yld3(lsim));
        //yld3   = yield(Fratio,F3);
        dyld_all   = (yld2_all - yld3_all)/df;                          // First derivative (to find the root of this)
        dyldp_all  = (yld2_all + yld3_all - 2.*yld1_all)/(.25*df*df);       // Second derivative (for Newton Raphson)
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if (breakout(lsim,sp)==0)
            {
              F1(lsim,sp)    -= dyld_all(lsim)/dyldp_all(lsim);
            }
            else
            {
              if (F1(lsim,sp)>5) 
                cout<<"Fmsy v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
              else      
                cout<<"Fmsy v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
            }
          }
        }
        cout<<"yld1_all(lsim=1) "<<yld1_all(1)<<endl;
      }// end for each iteration in 1:8
       for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
    break;

    case 3:
      //  Find max yeild unconditional, constrained such that none fall below 0.35
      F1.initialize();
      breakout.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            F1(lsim,sp) = (0.8*(mean(M1(sp)(1,nages(sp)) + M2(sp,nyrs))));
      breakout=0;
      // Newton Raphson stuff to go here
      for (int ii=1;ii<=CRrepn;ii++)
      {
        F2.initialize();
        F3.initialize();
        yld1.initialize();
        yld2.initialize();
        yld3.initialize();
        dyld_all.initialize();
        dyldp_all.initialize();
        //if (mceval_phase()&&(F1>5||F1<0.01)) 
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if ((F1(lsim,sp)>5||F1(lsim,sp)<0.01)) 
            {
              ii=8;
              if (F1(lsim,sp)>5) 
                F1(lsim,sp)=5.0; 
              else      
                F1(lsim,sp)=0.001; 
              breakout(lsim,sp)= 1;
            }
          }
        }
        F2     = F1 + df*.5;
        F3     = F2 - df;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
            yld1(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);       
            minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          // here is where the fix needs to occur maybe?

          }
          yld1_all(lsim)=sum(yld1(lsim));
          dvariable diff_val;
          diff_val=min(minB)-0.35;
          if(min(minB)<0.35)
            yld1_all(lsim)=mfexp(diff_val/0.35)*yld1_all(lsim);
        }

        // run projections with F2
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F2(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
            yld2(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);       
            minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          }
          yld2_all(lsim)=sum(yld2(lsim));
          dvariable diff_val;
          diff_val=min(minB)-0.35;
          if(min(minB)<0.35)
            yld2_all(lsim)=mfexp(diff_val/0.35)*yld2_all(lsim);
        }
        // run projections with F3
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F3(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
            yld3(lsim,sp)   = mean(tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);       
            minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          }
          yld3_all(lsim)=sum(yld3(lsim));
          dvariable diff_val;
          diff_val=min(minB)-0.35;
          if(min(minB)<0.35)
            yld3_all(lsim)=mfexp(diff_val/0.35)*yld3_all(lsim);
        }

        dyld_all   = (yld2_all - yld3_all)/df;                          // First derivative (to find the root of this)
        dyldp_all  = (yld2_all + yld3_all - 2.*yld1_all)/(.25*df*df);       // Second derivative (for Newton Raphson)
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if (breakout(lsim,sp)==0)
            {
              F1(lsim,sp)    -= dyld_all(lsim)/dyldp_all(lsim);
            }
            else
            {
              if (F1(lsim,sp)>5) 
                cout<<"Fmsy v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
              else      
                cout<<"Fmsy v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
            }
          }
        }
        cout<<"yld1_all(lsim=1) "<<yld1_all(1)<<endl;
      }// end for each iteration in 1:8
       for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
    break;
  }//end switch

FUNCTION void CONTROL_RULE_5(int sub)
  // max economic yield

  int CRrepn;
  CRrepn=10;
  int nyrs_avg;
  nyrs_avg=5;
  // just to get things going OKO replace with aka values
  dvar_vector price(1,nspp);
  price(1)=1;
  price(2)=1.3;
  price(3)=.2;

  dvar_vector est(1,nspp);
  dvar_vector minB(1,nspp);
  dvar_vector maxY(1,nspp);
  dvar_matrix CRtmpF(1,ntemp_scen,1,nspp);
  dvariable adjn;
  int ntemp_scen_tmp;
  int nn=nyrs+nyrs_fut;
  PF.initialize();
  dvar_matrix run_est(1,ntemp_scen,1,nspp);
  dvariable tmpEst; 
  dvariable tmpmin; 

  cout<<"NOW running control rule 5: sub= "<<sub<<endl;
  
  dvar_matrix F1(1,ntemp_scen,1,nspp);
  dvar_matrix F2(1,ntemp_scen,1,nspp);
  dvar_matrix F3(1,ntemp_scen,1,nspp);
  dvar_matrix yld1(1,ntemp_scen,1,nspp);
  dvar_matrix yld2(1,ntemp_scen,1,nspp);
  dvar_matrix yld3(1,ntemp_scen,1,nspp);
  dvar_vector yld1_all(1,ntemp_scen);
  dvar_vector yld2_all(1,ntemp_scen);
  dvar_vector yld3_all(1,ntemp_scen);
  dvar_matrix dyldp(1,ntemp_scen,1,nspp);
  dvar_matrix dyld(1,ntemp_scen,1,nspp);
  dvar_vector dyldp_all(1,ntemp_scen);
  dvar_vector dyld_all(1,ntemp_scen);
  imatrix breakout(1,ntemp_scen,1,nspp);
  double df=1.e-05;

  switch (sub)
  {
    default:
    cout<<"ERROR: Control rule 5 sub ? : please check ctl file"<<endl;exit(1);
    break;
    case 1:
      //  Find max yeild of each individual stock 
      F1.initialize();
      breakout.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            F1(lsim,sp) = (0.8*(mean(M1(sp)(1,nages(sp)) + M2(sp,nyrs))));
      breakout=0;
      // Newton Raphson stuff to go here
      for (int ii=1;ii<=8;ii++)
      {
        F2.initialize();
        F3.initialize();
        yld1.initialize();
        yld2.initialize();
        yld3.initialize();
        dyld.initialize();
        dyldp.initialize();
        //if (mceval_phase()&&(F1>5||F1<0.01)) 
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if ((F1(lsim,sp)>5||F1(lsim,sp)<0.01)) 
            {
              ii=8;
              if (F1(lsim,sp)>5) 
                F1(lsim,sp)=5.0; 
              else      
                F1(lsim,sp)=0.001; 
              breakout(lsim,sp)= 1;
            }
          }
        }
        F2     = F1 + df*.5;
        F3     = F2 - df;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld1(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        // run projections with F2
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F2(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld2(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        // run projections with F3
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F3(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld3(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        //yld3   = yield(Fratio,F3);
        dyld   = (yld2 - yld3)/df;                          // First derivative (to find the root of this)
        dyldp  = (yld2 + yld3 - 2.*yld1)/(.25*df*df);       // Second derivative (for Newton Raphson)
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if (breakout(lsim,sp)==0)
            {
              F1(lsim,sp)    -= dyld(lsim,sp)/dyldp(lsim,sp);
            }
            else
            {
              if (F1(lsim,sp)>5) 
                cout<<"Fmsy v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
              else      
                cout<<"Fmsy v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
            }
          }
        }
      }// end for each iteration in 1:8
       for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
    break;
    case 2:
      //  Find max yeild unconditional
      F1.initialize();
      breakout.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            F1(lsim,sp) = (0.8*(mean(M1(sp)(1,nages(sp)) + M2(sp,nyrs))));
      breakout=0;
      // Newton Raphson stuff to go here
      for (int ii=1;ii<=8;ii++)
      {
        F2.initialize();
        F3.initialize();
        yld1.initialize();
        yld2.initialize();
        yld3.initialize();
        dyld_all.initialize();
        dyldp_all.initialize();
        //if (mceval_phase()&&(F1>5||F1<0.01)) 
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if ((F1(lsim,sp)>5||F1(lsim,sp)<0.01)) 
            {
              ii=8;
              if (F1(lsim,sp)>5) 
                F1(lsim,sp)=5.0; 
              else      
                F1(lsim,sp)=0.001; 
              breakout(lsim,sp)= 1;
            }
          }
        }
        F2     = F1 + df*.5;
        F3     = F2 - df;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld1(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        for (int lsim=1;lsim<=ntemp_scen;lsim++)  yld1_all(lsim)=sum(yld1(lsim));
        // run projections with F2
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F2(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld2(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        for (int lsim=1;lsim<=ntemp_scen;lsim++)  yld2_all(lsim)=sum(yld2(lsim));
        // run projections with F3
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F3(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            yld3(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);
        for (int lsim=1;lsim<=ntemp_scen;lsim++)  yld3_all(lsim)=sum(yld3(lsim));
        //yld3   = yield(Fratio,F3);
        dyld_all   = (yld2_all - yld3_all)/df;                          // First derivative (to find the root of this)
        dyldp_all  = (yld2_all + yld3_all - 2.*yld1_all)/(.25*df*df);       // Second derivative (for Newton Raphson)
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if (breakout(lsim,sp)==0)
            {
              F1(lsim,sp)    -= dyld_all(lsim)/dyldp_all(lsim);
            }
            else
            {
              if (F1(lsim,sp)>5) 
                cout<<"Fmsy v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
              else      
                cout<<"Fmsy v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
            }
          }
        }
        cout<<"yld1_all(lsim=1) "<<yld1_all(1)<<endl;
      }// end for each iteration in 1:8
       for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
    break;

    case 3:
      //  Find max yeild unconditional, constrained such that none fall below 0.35
      F1.initialize();
      breakout.initialize();
      for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            F1(lsim,sp) = (0.8*(mean(M1(sp)(1,nages(sp)) + M2(sp,nyrs))));
      breakout=0;
      // Newton Raphson stuff to go here
      for (int ii=1;ii<=8;ii++)
      {
        F2.initialize();
        F3.initialize();
        yld1.initialize();
        yld2.initialize();
        yld3.initialize();
        dyld_all.initialize();
        dyldp_all.initialize();
        //if (mceval_phase()&&(F1>5||F1<0.01)) 
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if ((F1(lsim,sp)>5||F1(lsim,sp)<0.01)) 
            {
              ii=8;
              if (F1(lsim,sp)>5) 
                F1(lsim,sp)=5.0; 
              else      
                F1(lsim,sp)=0.001; 
              breakout(lsim,sp)= 1;
            }
          }
        }
        F2     = F1 + df*.5;
        F3     = F2 - df;
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
            yld1(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);       
            minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          }
          yld1_all(lsim)=sum(yld1(lsim));
          dvariable diff_val;
          diff_val=min(minB)-0.35;
          if(min(minB)<0.35)
            yld1_all(lsim)=mfexp(diff_val/0.35)*yld1_all(lsim);
        }

        // run projections with F2
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F2(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
            yld2(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);       
            minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          }
          yld2_all(lsim)=sum(yld2(lsim));
          dvariable diff_val;
          diff_val=min(minB)-0.35;
          if(min(minB)<0.35)
            yld2_all(lsim)=mfexp(diff_val/0.35)*yld2_all(lsim);
        }
        // run projections with F3
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F3(lsim,sp);
        RUN_PROJECTIONS();
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          minB.initialize();
          for (int sp=1;sp<=nspp;sp++)
          {
            yld3(lsim,sp)   = mean(price(sp)*tc_biom_hat_fut(lsim,sp)(nyrs_fut-5,nyrs_fut));//yield(Fratio,F1);       
            minB(sp)=min(elem_div(biomassSSB_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs)),biomassSSB0_hind2fut(lsim,sp)((nyrs_fut+nyrs)-nyrs,(nyrs_fut+nyrs))));
          }
          yld3_all(lsim)=sum(yld3(lsim));
          dvariable diff_val;
          diff_val=min(minB)-0.35;
          if(min(minB)<0.35)
            yld3_all(lsim)=mfexp(diff_val/0.35)*yld3_all(lsim);
        }

        dyld_all   = (yld2_all - yld3_all)/df;                          // First derivative (to find the root of this)
        dyldp_all  = (yld2_all + yld3_all - 2.*yld1_all)/(.25*df*df);       // Second derivative (for Newton Raphson)
        for (int lsim=1;lsim<=ntemp_scen;lsim++)
        {
          for (int sp=1;sp<=nspp;sp++)
          {
            if (breakout(lsim,sp)==0)
            {
              F1(lsim,sp)    -= dyld_all(lsim)/dyldp_all(lsim);
            }
            else
            {
              if (F1(lsim,sp)>5) 
                cout<<"Fmey v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
              else      
                cout<<"Fmey v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
            }
          }
        }
        cout<<"yld1_all(lsim=1) "<<yld1_all(1)<<endl;
      }// end for each iteration in 1:8
       for (int lsim=1;lsim<=ntemp_scen;lsim++)
          for (int sp=1;sp<=nspp;sp++)
            PF(lsim,sp) = F1(lsim,sp);
    break;
  }//end switch

FUNCTION RUN_F_PROFILE
  
  // Profile over Control rule multipliers...c_mult
  //ofstream Fprofile_report("../results/Future_Fprofile_report.rep");
  //ofstream future_report("../results/Future_report.rep");
  // ofstream future_Fprofile_report("results/Future_Fprofile_report.rep");
  WRITE_F_PROFILE_REPORT(0);
  // cout<<"number of fprofiles: np= "<<np<<endl;
  for(int sp=1;sp<=nspp;sp++) BtargetRep(sp)=Btarget;
  //for (crn=1;crn<=crnum;crn++)
  //for (int ff=1;ff<=np;ff++)
  for (int ff=1;ff<=1;ff++)
  {

    //cout<<"running F rate profile for ff = "<<ff<<" of "<<np<<" profiles "<<Fprofiles(ff)<<endl;
    for (int lsim=1;lsim<=ntemp_scen;lsim++)
      for (int sp=1;sp<=nspp;sp++)
        PF(lsim,sp)=FP_in(sp);
    RUN_PROJECTIONS(); // calculate unfished biomass and fished biomass at PF()
    WRITE_F_PROFILE_REPORT(ff); 

     if(dump_rep) cout<<"_______________________________________________"<<endl;
     if(dump_rep) cout<<"write Fprofile_report for ff = "<<ff<<endl;
     if(dump_rep) cout<<"_______________________________________________"<<endl;
  }//end crn
      write_R();
      if(dump_rep) cout<<"============= END: WRITE R    ========================="<<endl;  
      WRITE_FUT_REPORT(0);
      WRITE_FUT_REPORT(1);
      WRITE_PROJECTION_REPORT(0);
      WRITE_PROJECTION_REPORT(1);
      exit(1);

FUNCTION void WRITE_F_PROFILE_REPORT(int fprof_pass)
  
  if(fprof_pass==0)
  {
    cout<<"writing f profile report header"<<endl;
    future_Fprofile_report<<
    "profile_num"             <<" "<<
    "simMode"                 <<" "<<
    "MSMmode"                 <<" "<<
    "recMode"                 <<" "<<
    "fut_simulation"          <<" "<<
    "harvestMode"             <<" "<<
    "species"                 <<" "<<
    "Profile_frate"           <<" "<<
    "B_target"                <<" "<<
    FP_in                     <<" "<<
    "Frate_row_ff"            <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<" "<<
    "objective_fun"           <<" "<<
    "2012_harvest"            <<" "<<
    "2050_harvest"            <<" "<<
    "SSB/SSB0"                <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<" "<<
    "catch(biomass)"          <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<" "<<
    "recruits(numbers)"        <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<endl;
  }
  if(fprof_pass!=0)
  {
    cout<<"writing f profile report"<<endl;
    for(int itmp=1;itmp<=ntemp_scen;itmp++)
    {
      for (int spp1=1;spp1<=nspp;spp1++)
      {
        future_Fprofile_report   <<
        fprof_pass        <<" "<<
        simMode           <<" "<<
        msmMode           <<" "<<
        recMode           <<" "<<
        simset(itmp)      <<" "<<
        harvestMode       <<" "<<
        spp1              <<" "<<
        FP_in(spp1)       <<" "<<
        BtargetRep(spp1)  <<" "<<
        FP_in             <<" "<<
        "F_rate"          <<" "<<
        mean_F(spp1) * mfexp(F_dev(spp1))<<
        Frate_fut(itmp,spp1) <<" "<<
        obj_fun           <<" "<<
        tc_biom_hat_fut(itmp,spp1,1)          <<" "<<
        tc_biom_hat_fut(itmp,spp1,nyrs_fut)   <<" "<<
        "SSB"             <<" "<<
        biomassSSB_hind2fut(itmp,spp1)              <<" "<<
        "catch"           <<" "<<
        tc_biom_hat(spp1)<<tc_biom_hat_fut(itmp,spp1) <<" "<<
        "recruits"           <<" "<<
        R(spp1)<<R_fut(itmp,spp1)<<endl;
      }
    }
  }

FUNCTION write_R_new
  write_Recruitment();
  WRITE_REC_PARS();
  WRITE_REC_COVARS();


  if(dump_rep) cout<<"write_R_new"<<endl;
  R_output << "#nspp"       <<endl;   
  R_output << nspp           <<endl;
  R_output << "#nyrs"       <<endl;   
  R_output << nyrs           <<endl;
  R_output << "#styr"       <<endl;   
  R_output << styr           <<endl;
  R_output << "#fsh_age_type"       <<endl;   
  R_output << fsh_age_type           <<endl;
  R_output << "#srv_age_type"       <<endl;   
  R_output << srv_age_type           <<endl;
  if(do_fut)
  {
    R_output << "#nyrs_fut"       <<endl;   
    R_output << nyrs_fut           <<endl;
  }
  R_output << "#nlengths"       <<endl;   
  R_output << nlengths           <<endl;
  R_output << "#lengths"       <<endl;   
  R_output << lengths           <<endl;
  R_output << "#Tyrs"       <<endl;   
  R_output << Tyrs           <<endl;
  R_output << "#TempC"       <<endl;   
  R_output << BTempC_retro           <<endl;

  R_output << "#fT"       <<endl;   for (int sp=1;sp<=nspp;sp++){R_output << fT;}
  R_output << "#nages"    <<endl;   for (int sp=1;sp<=nspp;sp++){R_output << nages;}
  R_output << "#L2A_convert"    <<endl;   for (int sp=1;sp<=nspp;sp++){R_output << L2A_convert;}
  if (RationByAge==1)
  { 
    repR_4spp(ration2Age);
    repR_4spp(Pby_yr);
    repR_4spp(ConsumAge);
    repR_4spp(Consum_living);
    repR_4spp(Consum_livingAge);
    repR_4spp(pred_EbyAge);
  }
    /*
    else
    {
      R_output << "#ration2_"    <<endl;   for (int sp=1;sp<=nspp;sp++){R_output << ration2);
      R_output << "#Consum_"    <<endl;   for (int sp=1;sp<=nspp;sp++){R_output << Consum);
      R_output << "#Consum_living_"    <<endl;   for (int sp=1;sp<=nspp;sp++){R_output << Consum_living);

         R_output << "Consum_"       <<sp<<endl;     //Consum(1,nspp,1,nyrs,1,nlengths)  
      R_output << Consum(sp)           <<endl;
      R_output << "Consum_living_"       <<sp<<endl;     //Consum(1,nspp,1,nyrs,1,nlengths)  
      R_output << Consum_living(sp)           <<endl;
      R_output << "pred_E_"       <<sp<<endl;
      R_output << pred_E(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths)
    }
    R_output << "Eage_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_output << Eage(sp)           <<endl; 
    R_output << "Eage_hat_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_output << Eage_hat(sp)           <<endl; 
    R_output << "Nage_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_output << Nage(sp)           <<endl; 
    R_output << "Bage_"       <<sp<<endl;    //Bage(1,nspp,1,nyrs,1,prey_nages)
    R_output << Bage(sp)           <<endl; 
    R_output << "EageN_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_output << EageN(sp)           <<endl; 
    R_output << "BageN_"       <<sp<<endl;    //Bage(1,nspp,1,nyrs,1,prey_nages)
    R_output << BageN(sp)           <<endl; 
    R_output << "E_"       <<sp<<endl;
    R_output << E(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 

    R_output << "availB_"       <<sp<<endl;
    R_output << availB(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
    R_output << "Ucheck_"       <<sp<<endl;
    R_output << Ucheck(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
    if(do_fut)
    {
      R_output << "E_fut_"       <<sp<<endl;
      R_output << E_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
      R_output << "availB_fut_"       <<sp<<endl;
      R_output << availB_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
      R_output << "Ucheck_fut_"       <<sp<<endl;
      R_output << Ucheck_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
      R_output << "pred_E_fut_"       <<sp<<endl;
      R_output << pred_E_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
    }
    //cout<<"write_R Kirs end"<<endl;
    R_output << "wt_"    <<sp<<endl;
    R_output << wt(sp)          <<endl;
    R_output << "WtbyAge_"    <<sp<<endl;
    R_output << WtbyAge(sp)          <<endl;
    R_output << "WtbyL_"    <<sp<<endl;
    R_output << WtbyL(sp)          <<endl;
    R_output << "obs_catch_"    <<sp<<endl;
    R_output << tc_biom_obs(sp)          <<endl;
    R_output << "obs_catch_hat_"<<sp<<endl;
    R_output << tc_biom_hat(sp)          <<endl;
    R_output << "wt_in_stom_"   <<sp<<endl;
    R_output << mn_wt_stom(sp)       <<endl;
    R_output << "ration_"       <<sp<<endl;
    R_output << ration(sp)           <<endl;
    R_output << "wt_"           <<sp<<endl;
    R_output << wt(sp)               <<endl;
    R_output << "M1_"           <<sp<<endl;
    R_output << M1(sp)               <<endl;
    //R_output << "alpha_rec_"           <<sp<<endl;
    //R_output << alpha_rec(sp)               <<endl;
    //R_output << "beta_rec_"           <<sp<<endl;
    //R_output << beta_rec(sp)               <<endl;
    R_output << "M2_"           <<sp<<endl;
    R_output <<  M2(sp)              <<endl;
    R_output << "AvgN_"            <<sp<<endl;
    R_output << AvgN(sp)                <<endl;
    R_output << "AvgN2_"            <<sp<<endl;
    R_output << AvgN2(sp)                <<endl;
    R_output << "N_"            <<sp<<endl;
    R_output << N(sp)                <<endl;
    if(do_fut){
      R_output << "M2_fut_"           <<sp<<endl;
      R_output <<  M2_fut(sp)              <<endl;
      R_output << "R_fut_"            <<sp<<endl;
      R_output << R_fut(sp)                <<endl;
      R_output << "AvgN_fut_"            <<sp<<endl;
      R_output << AvgN_fut(sp)                <<endl;
      R_output << "N_fut_"            <<sp<<endl;
      R_output << N_fut(sp)                <<endl;
      R_output << "F_fut_"            <<sp<<endl;
      R_output << F_fut(sp)                <<endl;
    }
    //R_output << "Unfished_N_fut_"            <<sp<<endl;
    //R_output << N_futyr_0(sp)                <<endl;
    R_output << "F_"            <<sp<<endl;
    R_output << F(sp)                <<endl;

    R_output << "R_"            <<sp<<endl;
    R_output << R(sp)                <<endl;
    //R_output << "R_hat_"            <<sp<<endl;
    //R_output << R_hat(sp)                <<endl;
    //R_output << "Rbt_"            <<sp<<endl;
    //R_output << Rbt(sp)                <<endl;
    
    if(do_fut){
      R_output <<"ntemp_scen"<<" "<< "PF_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_output<<fitr<<" " << PF(fitr,sp)<<endl; 
    }       
    R_output << "S_"            <<sp<<endl;
    R_output << S(sp)                <<endl;
    R_output << "Biomass_"      <<sp<<endl;
    R_output << biomass(sp)          <<endl;
    if(do_fut){
      R_output <<"ntemp_scen"<<" "<< "Biomass_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_output<<fitr<<" " << biomass_hind2fut(fitr,sp)<<endl;     
   }
    R_output << "biomassByage_"      <<sp<<endl;
    R_output << biomassByage(sp)          <<endl;
    if(do_fut){
      R_output <<"ntemp_scen"<<" "<< "biomassByage_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_output<<fitr<<" " << biomassByage_hind2fut(fitr,sp)<<endl;     
    }

    R_output << "NByage_"      <<sp<<endl;
    R_output << NByage(sp)          <<endl;
    //R_output << "NByage_0_"      <<sp<<endl;
    //R_output << NByage_0(sp)          <<endl;
    
    R_output << "BiomassSSB_"      <<sp<<endl;
    R_output << biomassSSB(sp)          <<endl;
    if(do_fut){
      R_output <<"ntemp_scen"<<" "<<"BiomassSSB_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
         R_output<<fitr<<" " << biomassSSB_hind2fut(fitr,sp)<<endl;     
      R_output <<"ntemp_scen"<<" "<<"BiomassSSB0_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
        R_output<<fitr<<" " << biomassSSB0_hind2fut(fitr,sp)<<endl;     
    }
   //    R_output << "V_spawn_Biomass_"      <<sp<<endl;
   //    R_output << Vbiomass(sp)          <<endl;
   //    R_output << "spawn_Biomass_fished_"      <<sp<<endl;
   //    R_output << Vbiomass_hat(sp)          <<endl;
    R_output << "tc_hat_"       <<sp<<endl;
    R_output <<  tc_hat(sp)          <<endl;
    R_output << "tc_biom_obs_"       <<sp<<endl;
    R_output <<  tc_biom_obs(sp)          <<endl;
    R_output << "tc_obs_"<<sp      <<endl;
    R_output << tc_obs(sp)          <<endl;
    R_output << "tc_hat_"<<sp      <<endl;
    R_output << tc_hat(sp)          <<endl;
    R_output << "suit_other_"   <<sp<<endl;
    R_output << suit_other(sp)       <<endl;
    R_output << "yrs_fsh_comp_"  <<sp<<endl;
    R_output <<  yrs_fsh_comp(sp)    <<endl;
    R_output << "fsh_age_obs_"  <<sp<<endl;
    R_output <<  fsh_age_obs(sp)     <<endl;
    R_output << "fsh_age_hat_"  <<sp<<endl;
    R_output <<  fsh_age_hat(sp)     <<endl;
  }  
  R_output << "other_food"    <<endl;
  R_output << other_food      <<endl;

  for (int sp=1;sp<=nspp;sp++)
  {
    
    R_output << "avail_food_"    <<sp<<endl;
    R_output << avail_food(sp)      <<endl;
    R_output << "overlap_"      <<sp<<endl;
    dvar_matrix report_tmp(1,nyrs,1,nspp);
    for(int yr=1;yr<=nyrs;yr++)
      report_tmp(yr)=overlap(yr,sp);
    R_output << report_tmp <<endl;
    for (int pred_age=1;pred_age<=nages(sp);pred_age++){
      R_output << "suit_main_"    <<sp<<"_"<<pred_age<<endl;
      R_output << suit_main(sp,pred_age)    <<endl; 
    }  
    //Mtmp += AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/(avail_food(pred,i,pred_age));
    R_output << "srv_sel_"      <<sp<<endl;
    R_output << srv_sel(sp)          <<endl;
    R_output << "fsh_sel_"      <<sp<<endl;
    R_output << fsh_sel(sp)          <<endl;
  }

  
  for (int sp=1;sp<=nspp;sp++)
  {
    R_output << "yrs_srv_biom_" <<sp<<endl;
    R_output << yrs_srv_biom(sp)     <<endl;
    R_output << "srv_bio_"<<sp      <<endl;
    R_output << srv_biom(sp)         <<endl;
    R_output << "srv_bio_se_"<<sp      <<endl;
    R_output << srv_biom_se(sp)         <<endl;
    
    R_output << "srv_bio_hat_"<<sp  <<endl;
    R_output << srv_biom_hat(sp)     <<endl;

    R_output << "srv_age_yrs_"<<sp <<endl;
    R_output << yrs_srv_age(sp)     <<endl;
    R_output << "srv_age_hat_"<<sp <<endl;
    R_output << srv_age_hat(sp)     <<endl;
    R_output << "srv_age_obs_"<<sp <<endl;
    R_output << srv_age_obs(sp)     <<endl;    
    R_output << "M1_like_"<<sp     <<endl;
    R_output << M1_like(sp)         <<endl;
    R_output << "M2_like_"<<sp<<endl;
    R_output << M2_like(sp)    <<endl;
    R_output << "srv_sel_like_"<<sp<<endl;
    R_output << srv_sel_like(sp)    <<endl;
    R_output << "fsh_sel_like_"<<sp<<endl;
    R_output << fsh_sel_like(sp)    <<endl;    
    R_output << "srv_bio_like_"<<sp<<endl;
    R_output << srv_bio_like(sp)    <<endl;
    R_output << "fsh_cat_like_"<<sp<<endl;
    R_output << fsh_cat_like(sp)    <<endl;    
    R_output << "srv_age_like_"<<sp<<endl;
    R_output << srv_age_like(sp)    <<endl;
    R_output << "M2_pen_"<<sp<<endl;
    R_output << M2_pen(sp)    <<endl;
    R_output << "fsh_age_like_"<<sp<<endl;
    R_output << fsh_age_like(sp)    <<endl;
    R_output << "prey_consumed_"<<sp<<endl;
    R_output << prey_consumed(sp)    <<endl;

    if(do_fut){
      R_output <<"ntemp_scen"<<" "<< "prey_consumed_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_output<<fitr<<" " << prey_consumed_fut(fitr,sp)<<endl;    
    }    
  }
  R_output << "obj_fun"      <<endl;
  R_output << obj_fun        <<endl;
  R_output << "PLK_consumed"    <<endl;
  R_output << prey_consumed(1)    <<endl;
  R_output << "PCOD_consumed"     <<endl;
  R_output << prey_consumed(2)    <<endl;
  R_output << "ATF_consumed"       <<endl;
  R_output << prey_consumed(3)    <<endl;
  R_output << "eit_age_like "     <<endl;
  R_output << eit_age_like        <<endl;
  R_output << "eit_srv_like"      <<endl;
  R_output << eit_srv_like        <<endl;
  //R_output << U         <<endl;
  if(dump_rep) cout<<"write_R end"<<endl;
  */

FUNCTION write_RS_files
 // ofstream rs_out("results/ceattle_rs.dat");
  rs_out << nspp <<endl;
  rs_out << nyrs <<endl;
  rs_out << ln_mn_rec <<endl;
  for (int sp=1;sp<=nspp;sp++)
    rs_out << rec_dev(sp) <<endl;
  for (int sp=1;sp<=nspp;sp++)
    rs_out << biomassSSB(sp) <<endl;
  dvar_vector tmpration(1,nyrs);
  for (int sp=1;sp<=nspp;sp++){
    tmpration=0.;
    for (int yr=1;yr<=nyrs;yr++)
      tmpration(yr)=AvgN(sp,yr,1)*ration2Age(sp,yr,1);
    rs_out <<tmpration<<endl;
  } 
  for (int sp=1;sp<=nspp;sp++){
    tmpration=0.;
    for (int yr=1;yr<=nyrs;yr++)
      tmpration(yr)=AvgN(sp,yr)(2,nages(sp))*ration2Age(sp,yr)(2,nages(sp));
    rs_out <<tmpration<<endl;
  }
  rs_out<<12345<<endl;
  // rs_out.close();
FUNCTION write_fut_files
  // ofstream rand_report("results/randn.dat");
  if(dump_rep) cout<<"write future .dat files"<<endl;
  rand_report<<"#rmult"<<endl;
  for(int spp=1;spp<=nspp;spp++)
    rand_report<<rmult(spp)<<endl;
  rand_report<<"#rmultF"<<endl;
  for(int spp=1;spp<=nspp;spp++)
    rand_report<<rmultF(spp)<<endl;
  // rand_report.close();
FUNCTION write_output
  // cout<<"write_output"<<endl;
   // ofstream output_report("results/ceattle_output.dat");
    output_report << "#nselages"       <<endl;   
    output_report << nselages           <<endl;
    output_report << "#Recruitment(1,nspp,1,nyrs)"       <<endl;   
    output_report << R           <<endl;
    output_report << "#avail_food(1,nspp,1,nyrs,1,nages)"       <<endl;   
    output_report << avail_food           <<endl;  
    output_report << "#AvgN(1,nspp,1,nyrs,1,nages)"       <<endl;   
    output_report << AvgN         <<endl;
    output_report << "#N(1,nspp,1,nyrs+1,1,nages)"       <<endl;   
    output_report << N           <<endl;
    output_report << "#S(1,nspp,1,nyrs,1,nages)"       <<endl;   
    output_report << S           <<endl;
    output_report << "#fsh_sel(1,nspp,1,nages)"       <<endl;   
    output_report << fsh_sel           <<endl;
    output_report << "#biomass(1,nspp,1,nyrs)"       <<endl;   
    output_report << biomass           <<endl;
    output_report << "#biomassSSB(1,nspp,1,nyrs)"       <<endl;   
    output_report << biomassSSB           <<endl;
    // output_report.close();

  //output_report << "#Bage1(1,nspp,1,nyrs)"       <<endl;   
  //output_report << Bage1         <<endl; 
  //output_report << "#Eage1(1,nspp,1,nyrs)"       <<endl;   
  //output_report << Eage1         <<endl; 
  //cout<<"write_output end"<<endl;
// FUNCTION write_U
  /*
  cout<<"write_U start"<<endl;
  U_report << "U"       <<endl;   
  U_report << U       <<endl; 
  */
FUNCTION write_Recruitment  
 // ofstream Rec_report("results/ceattle_recruitment.rep");
  if(dump_rep) cout<<"write_Recruitment"<<endl;
  // Rec_report<<"spp"<<" "<<"year"<<" "<<"TempC"<<" "<<"Rec"<<" "<<"SSB"<<" "<<"M2of1YrOlds"<<endl;
  for(int spp=1;spp<=nspp;spp++)
    for(int yr=1;yr<=nyrs;yr++)
      Rec_report<<spp<<" "<<styr+yr-1<<" "<<TempC(yr)<<" "<<R(spp,yr)<<" "<<biomassSSB(spp,yr)<<" "<<M2(spp,yr,1)<<endl;
 Rec_report.close();
  if(do_fut)
  {
    ofstream Rec_fut_report("results/ceattle_fut_recruitment.rep");
    Rec_fut_report<<"spp"<<" "<<"future_itr"<<" "<<"year"<<" "<<"TempC"<<" "<<"Rec"<<" "<<"SSB"<<" "<<"M2of1YrOlds"<<endl;
    for(int spp=1;spp<=nspp;spp++){
      for(int fitr=1;fitr<=ntemp_scen;fitr++){
        for(int yr=1;yr<=nyrs_fut;yr++){
            Rec_fut_report<<spp<<" "<<fitr<<" "<<styr+nyrs+yr-1<<" "<<TempC_futUSE(fitr,yr)<<" "<<R_fut(fitr,spp,yr)<<" "<<biomassSSB_hind2fut(fitr,spp,nyrs+yr)<<" "<<M2_fut_all(fitr,spp,yr,1)<<endl;
        
        }
      }
    }
    // Rec_fut_report.close();
  }



FUNCTION WRITE_REC_COVARS  
  // ofstream covar_report("results/ceattle_recCovars.rep");
  if(dump_rep) cout<<"Write Recruitment Covars"<<endl;
   covar_report<<"Retrospective climate covariates"<<endl;
        covar_report<<"#**************************************************** "<<endl;
        covar_report<<"#**************************************************** "<<endl;
   for(int sp=1;sp<=nspp;sp++){  
        covar_report<<"species "<< sp<<endl;
        covar_report<<"#************************** "<<endl;
        for (int c=1;c<=ncovs(sp);c++)
          covar_report<< rs_cov(c)<<endl; 
    }
    for (itemp=1;itemp<=ntemp_scen;itemp++){
        covar_report<<"Future climate covariates, future scenario "<< itemp<<endl;
        covar_report<<"#**************************************************** "<<endl;
        covar_report<<"#**************************************************** "<<endl;
            covar_report<<"rs_cov"<<endl;
    for(int sp=1;sp<=nspp;sp++){  
        covar_report<<"species "<< sp<<endl;
        covar_report<<"#************************** "<<endl;
        for (int c=1;c<=ncovs(sp);c++)
          covar_report<< rs_cov_fut(itemp,c)<<endl; 
    }
        covar_report<<"rs_parm(sp,c)*rs_cov"<<endl;
    for(int sp=1;sp<=nspp;sp++){  
        covar_report<<"species "<< sp<<endl;
        covar_report<<"#************************** "<<endl;
        for (int c=1;c<=ncovs(sp);c++)
          covar_report<< rs_parm(sp,c)*rs_cov_fut(itemp,c)<<endl;
    }
  }
    // covar_report.close();

FUNCTION WRITE_REC_PARS  
  if(dump_rep) cout<<"Write Recruitment Pars"<<endl;
   // ofstream RecPar_report("results/ceattle_recPar.rep");
    RecPar_report<<"ncov"<<endl;
    RecPar_report<< ncov<<endl;
        RecPar_report<<"BevH"<<endl;
    RecPar_report<<  BevH<<endl;
        RecPar_report<<"BLM"<<endl;
    RecPar_report<< BLM<<endl;
        RecPar_report<<"LM"<<endl;
    RecPar_report<< LM<<endl;
        RecPar_report<<"logsigma_rs"<<endl;
    RecPar_report<< logsigma_rs<<endl;
        RecPar_report<<"logsigma_rs_std"<<endl;
    RecPar_report<< logsigma_rs_std<<endl;
        RecPar_report<<"mnRat"<<endl;
    RecPar_report<< mnRat<<endl;
        RecPar_report<<"sdRat"<<endl;
    RecPar_report<< sdRat<<endl;
        RecPar_report<<"log_aa_c"<<endl;
    RecPar_report<<  log_aa_c<<endl;
        RecPar_report<<"log_aa_c_std"<<endl;
    RecPar_report<< log_aa_c_std<<endl;
        RecPar_report<<"log_bb_c"<<endl;
    RecPar_report<< log_bb_c<<endl;
        RecPar_report<<"log_bb_c_std"<<endl;
    RecPar_report<< log_bb_c_std<<endl;

        RecPar_report<<"rs_parm"<<endl;
      for (int sp=1;sp<=nspp;sp++)
             RecPar_report<< rs_parm(sp)<<endl;  
                   RecPar_report<<"rs_parm_std"<<endl;
      for (int sp=1;sp<=nspp;sp++)
             RecPar_report<< rs_parm_std(sp)<<endl;
                   RecPar_report<<"cov_phase2"<<endl;
      for (int sp=1;sp<=nspp;sp++)
             RecPar_report<< cov_phase2(sp)<<endl;
                    RecPar_report<<"cov_type2"<<endl;
      for (int sp=1;sp<=nspp;sp++)
             RecPar_report<< cov_type2(sp)<<endl;
    // RecPar_report.close();


FUNCTION write_R  
  write_Recruitment();
  WRITE_REC_PARS();
  WRITE_REC_COVARS();
  // ofstream R_report("results/ceattle_R2.rep");

  if(dump_rep) cout<<"write_R"<<endl;
  R_report << "nspp"       <<endl;   
  R_report << nspp           <<endl;
  R_report << "nyrs"       <<endl;   
  R_report << nyrs           <<endl;
  R_report << "styr"       <<endl;   
  R_report << styr           <<endl;
  R_report << "fsh_age_type"       <<endl;   
  R_report << fsh_age_type           <<endl;
  R_report << "srv_age_type"       <<endl;   
  R_report << srv_age_type           <<endl;
  if(do_fut)
  {
    R_report << "nyrs_fut"       <<endl;   
    R_report << nyrs_fut           <<endl;
  }
  R_report << "nlengths"       <<endl;   
  R_report << nlengths           <<endl;
  R_report << "lengths"       <<endl;   
  R_report << lengths           <<endl;
  R_report << "Tyrs"       <<endl;   
  R_report << Tyrs           <<endl;
  R_report << "TempC"       <<endl;   
  R_report << BTempC_retro           <<endl;
 
  for (int sp=1;sp<=nspp;sp++)
  {
    //Kir's
   
    R_report << "pmature_"       <<sp<<endl;   // ration2(1,nspp,1,nyrs,1,nlengths)
    R_report << pmature(sp)           <<endl;
    R_report << "fT_"       <<sp<<endl;   // ration2(1,nspp,1,nyrs,1,nlengths)
    R_report << fT(sp)           <<endl;
    R_report << "L2A_convert_"       <<sp<<endl;     //Consum(1,nspp,1,nyrs,1,nlengths)  
    R_report << L2A_convert(sp)           <<endl;
    R_report << "of_stomKir_"       << sp <<endl;  
    R_report << of_stomKir(sp)           <<endl;
    R_report << "rec_dev_prior_"       <<sp<<endl;   
    R_report << rec_dev_prior(sp)           <<endl;
    R_report << "init_dev_prior_"       <<sp<<endl;  
    R_report << init_dev_prior(sp)           <<endl;
    R_report << "F_dev_prior_"       <<sp<<endl; 
    R_report << F_dev_prior(sp)           <<endl;

      

    if (RationByAge==1)
    {
      
      R_report << "ration2_"       <<sp<<endl;   // ration2(1,nspp,1,nyrs,1,nlengths)
      R_report << ration2Age(sp)           <<endl;
      R_report << "Pvalue_"       <<sp<<endl;   // ration2(1,nspp,1,nyrs,1,nlengths)
      R_report << Pby_yr(sp)           <<endl;

      R_report << "Consum_"       <<sp<<endl;     //Consum(1,nspp,1,nyrs,1,nlengths)  
      R_report << ConsumAge(sp)           <<endl;
      R_report << "Consum_living_"       <<sp<<endl;     //Consum(1,nspp,1,nyrs,1,nlengths)  
      R_report << Consum_livingAge(sp)           <<endl;
      R_report << "pred_E_"       <<sp<<endl;
      R_report << pred_EbyAge(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths)
        
    }
    else
    {
        
      R_report << "ration2_"       <<sp<<endl;   // ration2(1,nspp,1,nyrs,1,nlengths)
      R_report << ration2(sp)           <<endl;
      R_report << "Consum_"       <<sp<<endl;     //Consum(1,nspp,1,nyrs,1,nlengths)  
      R_report << Consum(sp)           <<endl;
      R_report << "Consum_living_"       <<sp<<endl;     //Consum(1,nspp,1,nyrs,1,nlengths)  
      R_report << Consum_living(sp)           <<endl;
      R_report << "pred_E_"       <<sp<<endl;
      R_report << pred_E(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths)
     
    }

    R_report << "Eage_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_report << Eage(sp)           <<endl; 
    R_report << "Eage_hat_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_report << Eage_hat(sp)           <<endl; 
    R_report << "Nage_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_report << Nage(sp)           <<endl; 
    R_report << "Bage_"       <<sp<<endl;    //Bage(1,nspp,1,nyrs,1,prey_nages)
    R_report << Bage(sp)           <<endl; 
    R_report << "EageN_"       <<sp<<endl;    //Eage(1,nspp,1,nyrs,1,prey_nages)  
    R_report << EageN(sp)           <<endl; 
    R_report << "BageN_"       <<sp<<endl;    //Bage(1,nspp,1,nyrs,1,prey_nages)
    R_report << BageN(sp)           <<endl; 
    R_report << "E_"       <<sp<<endl;
    R_report << E(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
  
    R_report << "availB_"       <<sp<<endl;
    R_report << availB(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
    R_report << "Ucheck_"       <<sp<<endl;
    R_report << Ucheck(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
    if(do_fut)
    {
      R_report << "E_fut_"       <<sp<<endl;
      R_report << E_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
      R_report << "availB_fut_"       <<sp<<endl;
      R_report << availB_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
      R_report << "Ucheck_fut_"       <<sp<<endl;
      R_report << Ucheck_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
      R_report << "pred_E_fut_"       <<sp<<endl;
      R_report << pred_E_fut(sp)           <<endl;     //E(1,nspp,1,nyrs,1,prey_nlengths) 
    }
     
    //cout<<"write_R Kirs end"<<endl;
    R_report << "wt_"    <<sp<<endl;
    R_report << wt(sp)          <<endl;
    R_report << "WtbyAge_"    <<sp<<endl;
    R_report << WtbyAge(sp)          <<endl;
    R_report << "WtbyL_"    <<sp<<endl;
    R_report << WtbyL(sp)          <<endl;
    R_report << "obs_catch_"    <<sp<<endl;
    R_report << obs_catch(sp)    <<sp<<endl;
    R_report << "obs_catch_hat_"<<sp<<endl;
    R_report << tc_biom_hat(sp)          <<endl;
    R_report << "wt_in_stom_"   <<sp<<endl;//
    R_report << mn_wt_stom(sp)       <<endl;//
    R_report << "ration_"       <<sp<<endl;//
    R_report << ration(sp)           <<endl;//
    R_report << "wt_"           <<sp<<endl;
    R_report << wt(sp)               <<endl;
    R_report << "M1_"           <<sp<<endl;
    R_report << M1(sp)               <<endl;
        

    //R_report << "alpha_rec_"           <<sp<<endl;
    //R_report << alpha_rec(sp)               <<endl;
    //R_report << "beta_rec_"           <<sp<<endl;
    //R_report << beta_rec(sp)               <<endl;
    R_report << "M2_"           <<sp<<endl;
    R_report <<  M2(sp)              <<endl;
    R_report << "AvgN_"            <<sp<<endl;
    R_report << AvgN(sp)                <<endl;
    R_report << "AvgN2_"            <<sp<<endl;
    R_report << AvgN2(sp)                <<endl;
    R_report << "N_"            <<sp<<endl;
    R_report << N(sp)                <<endl;
          
    if(do_fut){
      R_report << "M2_fut_"           <<sp<<endl;
      R_report <<  M2_fut(sp)              <<endl;
      // R_report << "R_fut_"            <<sp<<endl;
      // R_report << R_fut(itemp,sp)                <<endl;
      R_report << "AvgN_fut_"            <<sp<<endl;
      R_report << AvgN_fut(sp)                <<endl;
      R_report << "N_fut_"            <<sp<<endl;
      R_report << N_fut(sp)                <<endl;
      R_report << "F_fut_"            <<sp<<endl;
      R_report << F_fut(sp)                <<endl;
    }
    //R_report << "Unfished_N_fut_"            <<sp<<endl;
    //R_report << N_futyr_0(sp)                <<endl;
    R_report << "F_"            <<sp<<endl;
    R_report << F(sp)                <<endl;
         


    R_report << "R_"            <<sp<<endl;
    R_report << R(sp)                <<endl;
    //R_report << "R_hat_"            <<sp<<endl;
    //R_report << R_hat(sp)                <<endl;
    //R_report << "Rbt_"            <<sp<<endl;
    //R_report << Rbt(sp)                <<endl;
    
    if(do_fut){
      R_report <<"ntemp_scen"<<" "<< "PF_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_report<<fitr<<" " << PF(fitr,sp)<<endl; 
    }       
         

    R_report << "S_"            <<sp<<endl;
    R_report << S(sp)                <<endl;
    R_report << "Biomass_"      <<sp<<endl;
    R_report << biomass(sp)          <<endl;
    if(do_fut){
      R_report <<"ntemp_scen"<<" "<< "Biomass_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_report<<fitr<<" " << biomass_hind2fut(fitr,sp)<<endl;     
   }
    R_report << "biomassByage_"      <<sp<<endl;
    R_report << biomassByage(sp)          <<endl;
    if(do_fut){
      R_report <<"ntemp_scen"<<" "<< "biomassByage_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_report<<fitr<<" " << biomassByage_hind2fut(fitr,sp)<<endl;     
    }

    R_report << "NByage_"      <<sp<<endl;
    R_report << N(sp)          <<endl;
    //R_report << "NByage_0_"      <<sp<<endl;
    //R_report << NByage_0(sp)          <<endl;
    
    R_report << "BiomassSSB_"      <<sp<<endl;
    R_report << biomassSSB(sp)          <<endl;
    if(do_fut){
      R_report <<"ntemp_scen"<<" "<<"BiomassSSB_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
         R_report<<fitr<<" " << biomassSSB_hind2fut(fitr,sp)<<endl;     
      R_report <<"ntemp_scen"<<" "<<"BiomassSSB0_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
        R_report<<fitr<<" " << biomassSSB0_hind2fut(fitr,sp)<<endl;     
    }
   //    R_report << "V_spawn_Biomass_"      <<sp<<endl;
   //    R_report << Vbiomass(sp)          <<endl;
   //    R_report << "spawn_Biomass_fished_"      <<sp<<endl;
   //    R_report << Vbiomass_hat(sp)          <<endl;
    R_report << "tc_hat_"       <<sp<<endl;
    R_report <<  tc_hat(sp)          <<endl;
    R_report << "tc_biom_obs_"       <<sp<<endl;
    R_report <<  tc_biom_obs(sp)          <<endl;
    R_report << "tc_obs_"<<sp      <<endl;
    R_report << tc_obs(sp)          <<endl;
    R_report << "tc_hat_"<<sp      <<endl;
    R_report << tc_hat(sp)          <<endl;
    R_report << "suit_other_"   <<sp<<endl;
    R_report << suit_other(sp)       <<endl;
    R_report << "yrs_fsh_comp_"  <<sp<<endl;
    R_report <<  yrs_fsh_comp(sp)    <<endl;
    R_report << "fsh_age_obs_"  <<sp<<endl;
    R_report <<  fsh_age_obs(sp)     <<endl;
    R_report << "fsh_age_hat_"  <<sp<<endl;
    R_report <<  fsh_age_hat(sp)     <<endl;
  }  
  R_report << "other_food"    <<endl;
  R_report << other_food      <<endl;
  
  for (int sp=1;sp<=nspp;sp++)
  {
    
    R_report << "avail_food_"    <<sp<<endl;
    R_report << avail_food(sp)      <<endl;
    R_report << "overlap_"      <<sp<<endl;
    dvar_matrix report_tmp(1,nyrs,1,nspp);
    for(int yr=1;yr<=nyrs;yr++)
      report_tmp(yr)=overlap(yr,sp);
    R_report << report_tmp <<endl;
    for (int pred_age=1;pred_age<=nages(sp);pred_age++){
      R_report << "suit_main_"    <<sp<<"_"<<pred_age<<endl;
      R_report << suit_main(sp,pred_age)    <<endl; 

      for(int year = 1; year <= nyrs; year++){
            R_report << "stomKir_"    <<year<<"_"<<sp<<"_"<<pred_age<<endl;
            R_report << stomKir(year, sp, pred_age)    <<endl; 

            R_report << "stom_div_bio2_"    <<year<<"_"<<sp<<"_"<<pred_age<<endl;
            R_report << stom_div_bio2(year, sp, pred_age)    <<endl; 

          }
    }  
    //Mtmp += AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/(avail_food(pred,i,pred_age));
    R_report << "srv_sel_"      <<sp<<endl;
    R_report << srv_sel(sp)          <<endl;
    R_report << "fsh_sel_"      <<sp<<endl;
    R_report << fsh_sel(sp)          <<endl;
  }

   
  for (int sp=1;sp<=nspp;sp++)
  {
    R_report << "yrs_srv_biom_" <<sp<<endl;
    R_report << yrs_srv_biom(sp)     <<endl;
    R_report << "srv_bio_"<<sp      <<endl;
    R_report << srv_biom(sp)         <<endl;
    R_report << "srv_bio_se_"<<sp      <<endl;
    R_report << srv_biom_se(sp)         <<endl;
    
    R_report << "srv_bio_hat_"<<sp  <<endl;
    R_report << srv_biom_hat(sp)     <<endl;

    R_report << "srv_age_yrs_"<<sp <<endl;
    R_report << yrs_srv_age(sp)     <<endl;
    R_report << "srv_age_hat_"<<sp <<endl;
    R_report << srv_age_hat(sp)     <<endl;
    R_report << "srv_age_obs_"<<sp <<endl;
    R_report << srv_age_obs(sp)     <<endl;    
    R_report << "M1_like_"<<sp     <<endl;
    R_report << M1_like(sp)         <<endl;
    R_report << "M2_like_"<<sp<<endl;
    R_report << M2_like(sp)    <<endl;
    R_report << "srv_sel_like_"<<sp<<endl;
    R_report << srv_sel_like(sp)    <<endl;
    R_report << "fsh_sel_like_"<<sp<<endl;
    R_report << fsh_sel_like(sp)    <<endl;    
    R_report << "srv_bio_like_"<<sp<<endl;
    R_report << srv_bio_like(sp)    <<endl;
    R_report << "fsh_cat_like_"<<sp<<endl;
    R_report << fsh_cat_like(sp)    <<endl;    
    R_report << "srv_age_like_"<<sp<<endl;
    R_report << srv_age_like(sp)    <<endl;
    R_report << "M2_pen_"<<sp<<endl;
    R_report << M2_pen(sp)    <<endl;
    R_report << "fsh_age_like_"<<sp<<endl;
    R_report << fsh_age_like(sp)    <<endl;
    R_report << "prey_consumed_"<<sp<<endl;
    R_report << prey_consumed(sp)    <<endl;

    if(do_fut){
      R_report <<"ntemp_scen"<<" "<< "prey_consumed_fut_"      <<sp<<endl;
      for(int fitr=1;fitr<=ntemp_scen;fitr++)
       R_report<<fitr<<" " << prey_consumed_fut(fitr,sp)<<endl;    
    }    
  }
   
  R_report << "obj_fun"      <<endl;
  R_report << obj_fun        <<endl;
  R_report << "PLK_consumed"    <<endl;
  R_report << prey_consumed(1)    <<endl;
  R_report << "PCOD_consumed"     <<endl;
  R_report << prey_consumed(2)    <<endl;
  R_report << "ATF_consumed"       <<endl;
  R_report << prey_consumed(3)    <<endl;
  R_report << "eit_age_like "     <<endl;
  R_report << eit_age_like        <<endl;
  R_report << "eit_srv_like"      <<endl;
  R_report << eit_srv_like        <<endl;
 
  R_report << "obs_eit_age"      <<endl;
  R_report << obs_eit_age        <<endl;
  
  R_report << "obs_eit"      <<endl;
  R_report << obs_eit        <<endl;
  
  R_report << "eit_hat"      <<endl;
  R_report << eit_hat        <<endl;

  // R_report.close();
  //R_report << U         <<endl;
  if(dump_rep) cout<<"write_R end"<<endl;
FUNCTION void WRITE_FUT_REPORT(int crn_pass)
  int cr1=77;
  int sub1=77;
  if(crn_pass==0)
  {
    future_report<<
    "simMode"                 <<" "<<
    "MSMmode"                 <<" "<<
    "recMode"                 <<" "<<
    "fut_simulation"          <<" "<<
    "harvestMode"               <<" "<<
    "species"                 <<" "<<
    "control_rule"            <<" "<<
    "B_target"                <<" "<<
    "deltaB"                  <<" "<<
    "F_rate"                  <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<" "<<
    "objective_fun"           <<" "<<
    "2012_harvest"            <<" "<<
    "2050_harvest"            <<" "<<
    "SSB/SSB0"                <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<" "<<
    "catch(biomass)"          <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<" "<<
    "recruits(numbers)"        <<" "<<
    Tyrs                      <<fut_yrs(1,nyrs_fut) <<endl;
  }else{
    if(harvestMode==3){
      cr1=c_mult(crn_pass,1);//ie 3 1 - crn=3, sub = 1
      sub1=c_mult(crn_pass,2);
    }else{
      cr1=77;
      sub1=77;
    }
    dvar_vector deltaB(1,nspp);
    deltaB.initialize();
    for(int itmp=1;itmp<=ntemp_scen;itmp++)
    {
      for (int spp1=1;spp1<=nspp;spp1++)
      {
        future_report      <<
        simMode            <<" "<<
        msmMode            <<" "<<
        recMode            <<" "<<
        simset(itmp)         <<" "<<
        harvestMode          <<" "<<
        spp1               <<" "<<
        cr1 <<"."<<sub1      <<" "<<
        0                  <<" "<<
        0                  <<" "<<
        "F_rate"           <<" "<<
        mean_F(spp1) * 
          mfexp(F_dev(spp1))<<
        PF(itmp,spp1)*0     <<" "<<
        obj_fun             <<" "<<
        0                   <<" "<<
        0                   <<" "<<
        "SSB0"              <<" "<<
        biomassSSB0_hind2fut(itmp,spp1)    <<" "<<
        "catch"              <<" "<<
        tc_biom_hat(spp1)<<tc_biom_hat_fut(itmp,spp1)*0      <<" "<<
        "recruits"           <<" "<<
        R0_hind2fut(itmp,spp1)<<endl;
      }
      for (int spp1=1;spp1<=nspp;spp1++)
      {
        future_report     <<
        simMode           <<" "<<
        msmMode           <<" "<<
        recMode           <<" "<<
        simset(itmp)      <<" "<<
        harvestMode         <<" "<<
        spp1              <<" "<<
        cr1 <<"."<<sub1   <<" "<<
        BtargetRep(spp1)  <<" "<<
        deltaB(spp1)      <<" "<<
        "F_rate"          <<" "<<
        mean_F(spp1) * mfexp(F_dev(spp1))<<
        Frate_fut(itmp,spp1) <<" "<<
        obj_fun           <<" "<<
        tc_biom_hat_fut(itmp,spp1,1)          <<" "<<
        tc_biom_hat_fut(itmp,spp1,nyrs_fut)   <<" "<<
        "SSB"             <<" "<<
        biomassSSB_hind2fut(itmp,spp1)              <<" "<<
        "catch"           <<" "<<
        tc_biom_hat(spp1)<<tc_biom_hat_fut(itmp,spp1) <<" "<<
        "recruits"           <<" "<<
        R(spp1)<<R_fut(itmp,spp1)<<endl;
      }
    }
  }
//PF(itmp,spp1)     <<" "<<
FUNCTION void WRITE_FUT_PREDATION_REPORT(int crn_pass)
  int cr1=77;
  int sub1=77;
  dvariable Eattmp_age;
  dvariable Eattmp_age_N;
  if(crn_pass==0)
  {

    future_predation_report<<
    "msmMode"                 <<" "<<
    "simMode"                 <<" "<<
    "harvestMode"             <<" "<<
    "recmode"                 <<" "<<
    "year"                    <<" "<<
    "fut_sim"                 <<" "<<
    "control_rule"            <<" "<<
    "Btarget"                 <<" "<<
    "set_val_pred"            <<" "<<
    "pred"                    <<" "<<
    "pred_age"                <<" "<<
    "pred_ration_byage"       <<" "<<
    "suit_pred_prey"          <<" "<<
    "overlap_pred_prey"       <<" "<<
    "N_pred"                  <<" "<<
    "Biomass_pred"            <<" "<<
    "prey"                    <<" "<<
    "prey_age"                <<" "<<
    "N_prey"                  <<" "<<
    "Biomass_prey"            <<" "<<
    "N_eaten"                 <<" "<<
    "B_eaten"                 <<" "<<
    "pred_ration"             <<" "<<             
    "prey_consumed_fut0"      <<" "<<
    "prey_consumed_fut"       <<endl;
  }else{
    if(harvestMode==3){
      cr1=c_mult(crn_pass,1);//ie 3 1 - crn=3, sub = 1
      sub1=c_mult(crn_pass,2);
    }else{
      cr1=77;
      sub1=77;
    }
    
    for (int prey=1;prey<=nspp;prey++)
    {
      for (int prey_age=1;prey_age<=nages(prey);prey_age++)
        {
           for (int pred=1;pred<=nspp;pred++)
           {
             for (int pred_age=1;pred_age<=nages(pred);pred_age++)
              {
                Eattmp_age=0.;
                Eattmp_age_N=0.;
                Eattmp_age=(AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age));
                Eattmp_age_N=(AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/wt_fut(prey,i,prey_age);
                
                future_predation_report<<
                        msmMode                   <<" "<<
                        simMode                   <<" "<<
                        harvestMode               <<" "<<
                        recMode                   <<" "<<
                        i+styr+nyrs-1             <<" "<<
                        simset(itemp)             <<" "<<
                        cr1 <<"."<<sub1           <<" "<<
                        Btarget                   <<" "<<
                        set_val(pred)                   <<" "<<
                        pred                     <<" "<<
                        pred_age                  <<" "<<
                        ration2Age_fut(pred,i,pred_age)       <<" "<<
                        suit_main(pred,pred_age,prey,prey_age)          <<" "<<
                        overlap_fut(itemp,pred,prey,i)       <<" "<<
                        AvgN_fut(pred,i,pred_age) <<" "<<
                        AvgN_fut(pred,i,pred_age) * wt_fut(pred,i,pred_age) <<" "<<
                        prey                      <<" "<<
                        prey_age                  <<" "<<
                        AvgN_fut(prey,i,prey_age) <<" "<<
                        AvgN_fut(prey,i,prey_age) * wt_fut(prey,i,prey_age) <<" "<<
                        Eattmp_age_N                 <<" "<<
                        Eattmp_age                 <<" "<<
                        ( AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age) )            <<" "<<
                        prey_consumed0_fut(itemp,pred,i) <<" "<<
                        prey_consumed_fut(itemp,pred,i) <<endl;

              }
            }
        }
    }

    
  }
//PF(itmp,spp1)     <<" "<<
  
FUNCTION POLL_SS
  // Write to pollmsm.dat
  //ofstream poll_out("../poll/poll_msm.dat");
   ofstream poll_out("../../ceattle_SingleSppSub/poll/poll_msm.dat");
  // Fishery LFRq next
  // store catch at sp_age
  dvector C_tmp(1,15);
  double plusgrp;
  C_tmp(1,nages(1))    = value(catch_hat_fut(itemp,1,i)(1,nages(1)));
  // need to spread the 12 yr plust group in msm to the full 15 ages, apportioned by relative survival Ns (assuming M=.3)
  plusgrp = C_tmp(12);
  C_tmp(12) = plusgrp*.26;
  C_tmp(13) = plusgrp*.19;
  C_tmp(14) = plusgrp*.14;
  C_tmp(15) = plusgrp*.41;
  // poll_out << "# Year: "<<i+endyr<<endl; // Catch by MSM (in kt)
  // poll_out << NewCatch(1)/1000 <<endl; // Catch by MSM (in kt)
  poll_out << NewCatch(1) <<endl; // Catch by MSM (in kt)

  dvector ptmp(1,15);
  ptmp     = C_tmp/sum(C_tmp);
  ptmp     = sim_multin(ptmp, tau);
  poll_out << sum(C_tmp)*ptmp <<endl; // Catch by MSM (in N by sp_age)
  cout << sum(C_tmp)*ptmp <<endl; // Catch by MSM (in N by sp_age)

  // Biomass in bottom trawl survey NOT USED part of assessment
  double biomtmp ;
  double sigmatmp ;
  biomtmp  = value(srv_q(1)*srv_sel(1)*elem_prod(N_fut(1,i),wt_fut(1,i)))/1000;
  sigmatmp = sqrt(log(1+ srv_Mean_CV(1)*srv_Mean_CV(1)));
  biomtmp  *= exp(sigmatmp*randn(r_rep) - sigmatmp*sigmatmp/2.);
  // cout <<  biomtmp <<endl;
  // cout <<  biomtmp * srv_Mean_CV(1) <<endl;
  poll_out <<  biomtmp <<endl;
  poll_out <<  biomtmp * srv_Mean_CV(1) <<endl;
  // Now for numbers (population) of sp_age 2+?
  double   Nsrvtmp;
  Nsrvtmp  = value(srv_q(1)*srv_sel(1)*N_fut(1,i))/1000;
  // population std error
  poll_out <<  Nsrvtmp * srv_Mean_CV(1) <<endl;
  // cout <<  Nsrvtmp * srv_Mean_CV(1) <<endl;

  // BTS sp_age composition
  ptmp(1,12) = value(elem_prod(N_fut(1,i),srv_sel(1)));
  plusgrp = ptmp(12);
  // note that this is the weights for plus group up-scaling assuming (roughly) Z=.35
  ptmp(12) = plusgrp*.26;
  ptmp(13) = plusgrp*.19;
  ptmp(14) = plusgrp*.14;
  ptmp(15) = plusgrp*.41;
  ptmp/= sum(ptmp);
  poll_out << Nsrvtmp*sim_multin(ptmp, srv_age_n(1,nyrs_srv_age(1)))<<endl;
  
  // Produce estimate every year even if it's not used (generated by pm_rw.tpl) if (endyr==(yrs_eit_data(n_eit)+1))
  // obs_eit_data); if (endyr==(yrs_eit_data(n_eit)+1))
  // std_obs_eit_data); if (endyr==(yrs_eit_data(n_eit)+1))
  ptmp(1,12) = value(elem_prod(N_fut(1,i),eit_sel(n_eit)*eit_q));
  // eit_age_hat(i)     = elem_prod((N(1,eit_yrs_ind)),eit_sel(eit_yrs_ind)*eit_q); // No adjustment to mid-year
  plusgrp = ptmp(12);
  ptmp(12) = plusgrp*.26;
  ptmp(13) = plusgrp*.19;
  ptmp(14) = plusgrp*.14;
  ptmp(15) = plusgrp*.41;
  double tottmp;
  tottmp = sum(ptmp);
  ptmp /= tottmp;
  poll_out << tottmp * sim_multin(ptmp, eit_age_n(n_eit))<<endl;
  poll_out.close();

FUNCTION PCOD_SS
  // Write to pcod_msm.dat
  ofstream pcod_out("../../ceattle_SingleSppSub/pcod/pcod_msm.dat");
  // Fishery LFRq next
  // store catch at sp_age
  dvector C_tmp(1,12);
  cout << "pcod catch 1"<<endl;
  C_tmp(1,nages(2))    = value(catch_hat_fut(itemp,2,i)(1,nages(2)));
  // pcod_out << "# Catch "<<endyr+i-1 <<endl; // Catch by MSM (in t)
  pcod_out << NewCatch(2) <<endl; // Catch by MSM (in t)

  dvector ptmp(1,12);
  ptmp     = C_tmp/sum(C_tmp);
  ptmp     = sim_multin(ptmp, tau);
  // pcod_out << "# Catch at sp_age "<<endyr+i-1 <<endl; // Catch by MSM (in t)
  pcod_out << sum(C_tmp)*ptmp <<endl; // Catch by MSM (in N by sp_age)
  cout << "pcod catch "<<endl<< sum(C_tmp)*ptmp <<endl; // Catch by MSM (in N by sp_age)

  // Biomass in bottom trawl survey 
  double biomtmp ;
  double sigmatmp ;
  biomtmp  = value(srv_q(2)*srv_sel(2)*elem_prod(N_fut(2,i),wt_fut(2,i)));
  sigmatmp = sqrt(log(1+ srv_Mean_CV(2)*srv_Mean_CV(2)));
  biomtmp  *= exp(sigmatmp*randn(r_rep) - sigmatmp*sigmatmp/2.);
  cout <<  "Pcod biom "<< biomtmp <<endl;
  cout <<  "Pcod biom cv  "<< biomtmp * srv_Mean_CV(2) <<endl;
  // pcod_out << "# survey biomass "<<endyr+i-1 <<endl; // Catch by MSM (in t)
  pcod_out<<  biomtmp <<endl;
  // pcod_out << "# survey biomass std err "<<endyr+i-1 <<endl; // Catch by MSM (in t)
  pcod_out<<  biomtmp * srv_Mean_CV(2) <<endl;
  biomtmp  = value(srv_q(2)*srv_sel(2)*N_fut(2,i));

  // BTS length composition
  dvector pltmp(1,25) ;
  pltmp = value(elem_prod(N_fut(2,i),srv_sel(2))) * PCOD_agesize;
  pltmp /= sum(pltmp);
  // note that this is the weights for plus group up-scaling assuming (roughly) Z=.35
  pcod_out<< sim_multin(pltmp, srv_age_n(2,nyrs_srv_age(2)))<<endl;
  pcod_out.close();
  

FUNCTION ATF_SS
  // Write to atfmsm.dat
  ofstream atf_out("../../ceattle_SingleSppSub/atf/atf_msm.dat");
  // = elem_prod(elem_div(F_fut(k,i),Z_fut(k,i)) ,elem_prod(1.-mfexp(-Z_fut(k,i)) , N_fut(k,i)));
  // survey length compositions first (Females, then males)
  dvector N_tmp(1,nages(3));
  // set vector to correct begin-year N (ixxxx)
  N_tmp = value(N_fut(3,i));
  // get selected ages
  N_tmp = value(elem_prod(srv_sel(3) , N_tmp));
  // matrix multiplication of N females at sp_age by sp_age-size conversion matrix for females
  dvector ptmp(1,srv_age_bins(3)) ;
  ptmp = elem_prod(N_tmp,  ATF_SexRatio) * ATF_agesize_F;
  ptmp /= sum(ptmp);
  atf_out << sim_multin(ptmp, srv_age_n(3,nyrs_srv_age(3)))<<endl;
  // matrix multiplication of N males at sp_age by sp_age-size conversion matrix for males
  ptmp = elem_prod(N_tmp,  1.-ATF_SexRatio) * ATF_agesize_M;
  ptmp /= sum(ptmp);
  atf_out << sim_multin(ptmp, srv_age_n(3,nyrs_srv_age(3)))<<endl;

  // Fishery LFRq next
  // store catch at sp_age
  dvector C_tmp(1,nages(3));
  C_tmp = value(catch_hat_fut(itemp,3,i)(1,nages(3)));
  // matrix multiplication of C females at sp_age by sp_age-size conversion matrix for females
  ptmp = elem_prod(C_tmp,  ATF_SexRatio) * ATF_agesize_F;
  ptmp /= sum(ptmp);
  atf_out << sim_multin(ptmp, tau) <<endl;

  // matrix multiplication of C males at sp_age by sp_age-size conversion matrix for males
  ptmp = elem_prod(C_tmp,1.-ATF_SexRatio) * ATF_agesize_M;
  ptmp /= sum(ptmp);
  atf_out << sim_multin(ptmp, tau) <<endl;
  atf_out << NewCatch(3) <<endl; // Catch by MSM (in t)
  // atf_out << C_tmp * wt_fut(3,i) <<endl;
  double biomtmp ;
  double sigmatmp ;
  biomtmp = value(srv_q(3)*srv_sel(3)*elem_prod(N_fut(3,i),wt_fut(3,i)));
  sigmatmp = sqrt(log(1+srv_Mean_CV(3)*srv_Mean_CV(3)));
  biomtmp *= exp(sigmatmp*randn(r_rep) - sigmatmp*sigmatmp/2.);
  atf_out <<  biomtmp <<endl;
  atf_out <<  biomtmp * srv_Mean_CV(3) <<endl;
  atf_out.close();

FUNCTION double Get_TAC(const int& nyrstmp, const int& ispp, const int& iyr, const double& tac_lastyr)
  double tac_nextyr;
  dvector xtmp(1,nyrstmp);
  xtmp.fill_seqadd(1,1);
  xtmp -= mean(xtmp);
  dvector ytmp = log(srv_biom(ispp)(iyr-nyrstmp,iyr)) - mean(log(srv_biom(ispp)(iyr-nyrstmp,iyr)));
  double lam = sum(elem_div(elem_prod(xtmp,ytmp),elem_prod(xtmp,xtmp)));
  if (lam<0)
    tac_nextyr = (1.- 1.5*std::abs(lam) )* tac_lastyr;
  else
    tac_nextyr = (1.+ 3.*lam) * tac_lastyr;
  return tac_nextyr;

FUNCTION dvector sim_multin(const dvector& phat, const int& Nsam)
  // Function to do what's expected for multinomial, returns proportions given sample size and phat
  ivector xbin(1,Nsam);
  dvector ff(1,phat.indexmax());
  ff.initialize();
  xbin.fill_multinomial(r_rep,phat);
  for (int ii=1;ii<=Nsam;ii++)
    ff(xbin(ii))++;
  return ff /= sum(ff);
//==============================================================================
//// THESE FUNCTIONS ARE NOT PRESENTLY USED BUT SHOULD NOT BE DELETED
//==============================================================================

FUNCTION void CALC_CONSUM(int C2_number)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// Calculate annual ration2 C=24*0.0134*mfexp(0.0115*TempC)*fdays(predd)*S  units are g of food/g of predator
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  for(i=1;i<=nyrs;i++)
  //  {
  // this appears to be working correctly - Kir has checked the U and E matrices, as well as ration 
 cout<<"PvalueAdj "<<PvalueAdj<<endl;
  if(C2_number==0){
      for(int predd=1;predd<=nspp;predd++)
      {// for each pred predd  
        Consum(predd,i)=24*0.0134*mfexp(0.0115*TempC(i))*91.25*elem_prod(S2(predd),WtbyL(predd));//kg/pred.yr
        Consum_living(predd,i)=Consum(predd,i); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        if(C_model(predd)==1)
          Consum(predd,i)=elem_prod(((CA(predd))*pow(WtbyL(predd),CB(predd))*fT(predd,i)*Pvalue(predd)*PvalueAdj(predd)*fdays(predd)),WtbyL(predd));//g/pred.yr      
        //ration2(predd,i)=elem_prod(Consum(predd,i),AvgN2(predd,i));//annual ration g/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
        ration2(predd,i)=Consum(predd,i)/1000.;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
      }//end predd  
  }
  if(C2_number==1){
      for(int predd=1;predd<=nspp;predd++)
      {// for each pred predd  
        Consum_fut(predd,i)=24.*0.0134*mfexp(0.0115*TempC_futUSE(itemp,i))*91.25*elem_prod(S2(predd),WtbyL(predd));//kg/pred.yr
        Consum_living_fut(predd,i)=Consum_fut(predd,i); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        if(C_model(predd)==1)
          Consum_fut(predd,i)=elem_prod(((CA(predd))*pow(WtbyL(predd),CB(predd))*fT_fut(itemp,predd,i)*Pvalue(predd)*PvalueAdj(predd)*fdays(predd)),WtbyL(predd));//g/pred.yr      
        //ration2(predd,i)=elem_prod(Consum(predd,i),AvgN2(predd,i));//annual ration g/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
        ration2_fut(predd,i)=Consum_fut(predd,i)/1000.;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
      }//end predd

  }  
  //for (int preyy=1;preyy<=nspp;preyy++)
  //  availB(preyy,i)=elem_prod(AvgN2(preyy,i),WtbyL(preyy))/1000;  //kgB/1000; //kg
  for (int preyy=1;preyy<=nspp;preyy++)
    availB(preyy,i)=elem_prod(AvgN2(preyy,i),WtbyL(preyy))/1000.;  //kgB/1000; //kg  
  dvar_matrix U4bio(1,maxL,1,maxL); // demand
  U4bio.initialize();  
  for (int predd=1;predd<=nspp;predd++)
    for (int preyy=1;preyy<=nspp;preyy++)
      U4bio+=Unew(predd,preyy);

  dvar_vector ttmp(1,maxL);
  ttmp.initialize();  
  ttmp=colsum(U4bio);  
  
  for (int preyy=1;preyy<=nspp;preyy++)
  {
    E(preyy,i).initialize();
    B4E(preyy,i).initialize();
    pred_E(preyy,i).initialize();
    //availB(preyy,i).initialize();
    Ucheck(preyy,i).initialize();
    for (int predd=1;predd<=nspp;predd++)
    {    
      dvar_matrix UcolSums(1,nlengths(predd),1,nlengths(preyy)); // demand
      dvar_matrix UcolSums2(1,nlengths(predd),1,nlengths(preyy)); // demand
      dvar_matrix Etmp3(1,nlengths(predd),1,nlengths(preyy)); // demand
      dvar_matrix Btmp3(1,nlengths(predd),1,nlengths(preyy)); // available B
      dvar_matrix Utmp3(1,nlengths(predd),1,nlengths(preyy)); // available B
      UcolSums.initialize();  
      UcolSums2.initialize();  
      Etmp3.initialize();  
      Btmp3.initialize();  
      Utmp3.initialize();  
      for(int predL=1;predL<=nlengths(predd);predL++)
        UcolSums2(predL)=Unew(predd,preyy,predL)(1,nlengths(preyy));
      //cout<<"UcolSums"<<UcolSums<<endl;  
      for(int predL=1;predL<=nlengths(predd);predL++)
        for(int preyL=1;preyL<=nlengths(preyy);preyL++)
          if(ttmp(preyL)>0)
            UcolSums(predL,preyL)=(UcolSums2(predL,preyL)/ttmp(preyL));
      for(int predL=1;predL<=nlengths(predd);predL++){
        if(C2_number==0)
        {
            
            //Etmp3(predL) =Unew(predd,preyy,predL)(1,nlengths(preyy))*ration2(predd,i,predL)*AvgN2(predd,i,predL);
            Etmp3(predL)   =Unew(predd,preyy,predL)(1,nlengths(preyy))*ration2(predd,i,predL)*AvgN2(predd,i,predL);
            //Btmp3(predL) =elem_prod(Unew(predd,preyy,predL)(1,nlengths(preyy)),elem_prod(WtbyL(preyy),AvgN2(preyy,i)));
            //Btmp3(predL) =elem_prod(UcolSums(predL),elem_prod(WtbyL(preyy),AvgN2(preyy,i)));
            Btmp3(predL)   =elem_prod(UcolSums(predL),elem_prod(WtbyL(preyy),AvgN2(preyy,i)));
            //Btmp3(predL) =elem_prod(U4bio(predL)(1,nlengths(preyy)),elem_prod(WtbyL(preyy),AvgN2(preyy,i)));
            Utmp3(predL)   =Unew(predd,preyy,predL)(1,nlengths(preyy));
        }
        if(C2_number==1)
        {
            Etmp3(predL) =Unew(predd,preyy,predL)(1,nlengths(preyy))*ration2_fut(predd,i,predL)*AvgN2_fut(predd,i,predL);
            Btmp3(predL) =elem_prod(UcolSums(predL),elem_prod(WtbyL(preyy),AvgN2_fut(preyy,i)));
            Utmp3(predL) =Unew(predd,preyy,predL)(1,nlengths(preyy));
            
        }  
      }  
      E(preyy,i)            +=colsum(Etmp3); //E(preyy,i,preyL),Unew(i,predd,preyy,predL,preyL) // should be adding across predators?
      B4E(preyy,i)          +=colsum(Btmp3);
      Ucheck(preyy,i)       +=colsum(Utmp3); //E(preyy,i,preyL),Unew(i,predd,preyy,predL,preyL) // should be adding across predators?
      //availB(preyy,i)     +=colsum(Btmp3); //E(preyy,i,preyL),Unew(i,predd,preyy,predL,preyL) // should be adding across predators?
      pred_E(preyy,i,predd) =sum(colsum(Etmp3));  
      pred_B(preyy,i,predd) =sum(colsum(Btmp3));  
    }//end for pred predd
  }
  for(int preyy=1;preyy<=nspp;preyy++)
  {
    dvar_matrix tmpN2(1,nages(preyy),1,nlengths(preyy));
    dvar_matrix tmpB(1,nages(preyy),1,nlengths(preyy));
    dvar_matrix tmpE(1,nages(preyy),1,nlengths(preyy));
    dvar_vector availB2(1,nlengths(preyy));
    dvar_vector Eagetmp(1,nages(preyy));
    Eagetmp.initialize();
    tmpN2.initialize();
    tmpB.initialize();
    tmpE.initialize();
    availB2.initialize();
    if(C2_number==0)
    {
        //availB2=elem_prod(WtbyL(preyy),AvgN2(preyy,i));
        availB2=elem_prod(WtbyL(preyy),N2(preyy,i));
        for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
        {
          //tmpN2(sp_age) =elem_prod(AvgN2(preyy,i),L2A_convert(preyy,sp_age));
          tmpN2(sp_age)   =elem_prod(AvgN2(preyy,i),L2A_convert(preyy,sp_age));
          //tmpB(sp_age)  =elem_prod((availB2),L2A_convert(preyy,sp_age));
          //tmpB(sp_age)  =elem_prod((availB(preyy,i)),L2A_convert(preyy,sp_age));
          tmpB(sp_age)    =elem_prod((B4E(preyy,i)),L2A_convert(preyy,sp_age));
          tmpE(sp_age)    =elem_prod((E(preyy,i)),L2A_convert(preyy,sp_age));
        }  
        Nage(preyy,i)   =rowsum(tmpN2);
        Eage(preyy,i)   =rowsum(tmpE);
        //Bage(preyy,i) =rowsum(tmpB);
        Bage(preyy,i)   =elem_prod(N(preyy,i),wt(preyy,i));
        EageN(preyy,i)  =elem_div(Eage(preyy,i),wt(preyy,i));  // N consumed
        BageN(preyy,i)  =elem_div(Bage(preyy,i),wt(preyy,i));  // N available
        Eagetmp         =Eage(preyy,i); // temporary Eaten (aka demand) by sp_age matrix for the preyy spp and the year
    }
    if(C2_number==1){
        availB2=elem_prod(WtbyL(preyy),N2_fut(preyy,i));
        for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
        {
          tmpN2(sp_age) =elem_prod(AvgN2_fut(preyy,i),L2A_convert(preyy,sp_age));
          tmpB(sp_age)  =elem_prod((B4E(preyy,i)),L2A_convert(preyy,sp_age));
          tmpE(sp_age)  =elem_prod((E(preyy,i)),L2A_convert(preyy,sp_age));
        }  
        Nage_fut(preyy,i)  =rowsum(tmpN2);
        Eage_fut(preyy,i)  =rowsum(tmpE);
        //Bage(preyy,i)    =rowsum(tmpB);
        Bage_fut(preyy,i)  =elem_prod(N_fut(preyy,i),wt_fut(preyy,i));
        EageN_fut(preyy,i) =elem_div(Eage_fut(preyy,i),wt_fut(preyy,i));  // N consumed
        BageN_fut(preyy,i) =elem_div(Bage_fut(preyy,i),wt_fut(preyy,i));  // N available
        Eagetmp            =Eage_fut(preyy,i); // temporary Eaten (aka demand) by sp_age matrix for the preyy spp and the year
    }
  }
  //cout<<"sum(B4E(preyy=1,i))"<<sum(B4E(1,i))<<endl;
  //cout<<"sum(elem_prod(WtbyL(preyy=1),AvgN2(preyy=1,i)))"<<sum(elem_prod(WtbyL(1),AvgN2(1,i)))<<endl;exit(1);
FUNCTION void CALC_U(int pass_number)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////based on sp_age Calculate (U) effective prey size freq (Npa/sum(N))*switch -- rows=prey, cols=preds, can be in calcs section if Npa/Np doesn't change each year
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout<<"CHANGE THIS TO / mean() from / max()!!!"<<endl;
  dvar_matrix preypropN(1,nspp,1,nlengths);
  dvar_matrix preypropW(1,nspp,1,nlengths);
  dvar_matrix preypropUSE(1,nspp,1,nlengths);
  preypropUSE.initialize();
  preypropN.initialize();
  preypropW.initialize();
  Unew.initialize();
  //first calculate AvgN2
  for(int spp=1;spp<=nspp;spp++)
  {
    dvar_matrix tmpN(1,nages(spp),1,nlengths(spp));
    dvar_matrix tmpN2(1,nages(spp),1,nlengths(spp));
    tmpN.initialize();
    tmpN2.initialize();
    N2(spp,i)=0;
    for(int sp_age=1;sp_age<=nages(spp);sp_age++)
      tmpN(sp_age)=N(spp,i,sp_age)*A2L_matrix(spp,sp_age);  
    N2(spp,i)=colsum(tmpN);
    for(int sp_age=1;sp_age<=nages(spp);sp_age++)
    {  
      dvar_vector sumN(1,nlengths(spp));
      sumN.initialize();
      sumN=colsum(tmpN);
      L2A_convert(spp,sp_age)=0;
      for(int sp_L=1;sp_L<=nlengths(spp);sp_L++)
        if(sumN(sp_L)>0)
          L2A_convert(spp,sp_age,sp_L)=(tmpN(sp_age,sp_L)/sumN(sp_L));
    }  
    for(int sp_age=1;sp_age<=nages(spp);sp_age++)
      tmpN2(sp_age)=elem_prod(N2(spp,i),L2A_convert(spp,sp_age));  
    Nage(spp,i)=rowsum(tmpN2);
  }
  for(int predd=1;predd<=nspp;predd++)
  {
    for(int preyy=1;preyy<=nspp;preyy++)
    {
      if(Umode==0){
        Unew(predd,preyy)=Uobs(predd,preyy);
      }else{
        // THIS IS FROM THE OLDER METHOD - remove?
        // dvar_vector preyprop(1,nlengths(preyy));
        // dvar_vector pref(1,nlengths(predd));
        // preyprop.initialize();
        // pref.initialize();
        // pref=1;  
        // if(pass_number==0){
        //   if (sum(N2(preyy,i))>0)
        //     preyprop=N2(preyy,i)/max(N2(preyy,i));  
        //   else
        //     preyprop(preyy)=0;
        //   if(useWt(predd)==1)
        //   {// if the model is set to forage as a function of biomass....
        //     preyprop=elem_prod(N2(preyy,i),(WtbyL(preyy)));
        //     cout<<"CHANGE THIS TO / mean() from / max()!!!"<<endl;
        //     if (sum(preyprop)>0)
        //       preyprop=preyprop/max(preyprop);// proportion by weight <- change this to propby weight by year !!!!  ? need to pre-allocate?        
        //     else
        //       preyprop(preyy)=0;      
        //   } 
        // }
        // if(pass_number==1){
        //   if (sum(N2_fut(preyy,i))>0)
        //     preyprop=N2_fut(preyy,i)/max(N2_fut(preyy,i));  
        //   else
        //     preyprop(preyy)=0;
        //   if(useWt(predd)==1)
        //   {// if the model is set to forage as a function of biomass....
        //     preyprop=elem_prod(N2_fut(preyy,i),(WtbyL(preyy)));
        //     if (sum(preyprop)>0)
        //       preyprop=preyprop/max(preyprop);// proportion by weight <- change this to propby weight by year !!!!  ? need to pre-allocate?        
        //     else
        //       preyprop(preyy)=0;      
        //   }         
        // }
        // preypropUSE(preyy)=preyprop;      
        // dvar_matrix Utmp(1,nlengths(predd),1,nlengths(preyy)); //U(#yrs,pred#,Prey#,preysize,predsize)
        // dvar_matrix switch1(1,nlengths(predd),1,nlengths(preyy)); //U(#yrs,pred#,Prey#,preysize,predsize)
        // dvar_matrix avail(1,nlengths(predd),1,nlengths(preyy));
        // dvar_matrix part1(1,nlengths(predd),1,nlengths(preyy));
        // dvar_matrix ttmp(1,nlengths(predd),1,nlengths(preyy));
        // Utmp.initialize();
        // switch1.initialize();
        // avail.initialize();
        // part1.initialize();
        // ttmp.initialize();
        // switch1=1; //pre-allocate a vector of 1s 
        // part1=1;  
        // for(int predL=1;predL<=nlengths(predd);predL++) // for each predL size of predd
        // {  
        //   ttmp(predL)=(lengths(preyy)-l2(predd,preyy)*Limit(predd,predL))/l1(predd,preyy);
        //   for(int preyL=1;preyL<=nlengths(preyy);preyL++)  
        //     if(ttmp(predL,preyL)<0)
        //       ttmp(predL,preyL)=0;        
        //   switch1(predL)(1,nlengths(preyy))=mfexp(-ttmp(predL));
        // }
        // for(int predL=1;predL<=nlengths(predd);predL++) // for each predL size of predd
        // {
        //   avail(predL)=elem_prod(switch1(predL),preypropUSE(preyy));
        //   if(sum(avail(predL))>0)
        //     avail(predL)=avail(predL)/sum(avail(predL));
        //   else
        //     avail(predL)=0;
        //   part1(predL)=0;
        //   for(int preyL=1;preyL<=nlengths(preyy);preyL++) // for each j size of preyy
        //     if(avail(predL,preyL)>0)
        //       part1(predL,preyL)=( prefa(predd,preyy)*(pow(avail(predL,preyL),prefb(predd,preyy))) )/(1+(prefa(predd,preyy)*(pow(avail(predL,preyL),prefb(predd,preyy)))) ) ;
        //   if(sum(part1(predL))>0)
        //     part1(predL)=part1(predL)/(sum(part1(predL)));// standardize across all prey sizes that the pred predL will eat.        
        //   else
        //     part1(predL)=0;
        //   Utmp(predL)=K(predd,preyy,predL)*part1(predL);
        //   Unew(predd,preyy,predL)(1,nlengths(preyy))=K(predd,preyy,predL)*part1(predL);
        //   //Utmp is a vector of 1xlengths (predL) of pred lengths, switch is a matrix of prey X pred L, and K is a 1xk vector of diet proportions  
        // }//end pred lengths k
      }// end if Umode==0
      // Is this correct below?
      cout<<"Calc_U is here"<<endl;
      for(int predL=1;predL<=nlengths(predd);predL++)
          if(sum(Unew(predd,preyy,predL))>0)
            Unew(predd,preyy,predL)=(Unew(predd,preyy,predL)/sum(Unew(predd,preyy,predL)))*K(predd,preyy,predL);
          else
            Unew(predd,preyy,predL)=0;
      UnewAVG(predd,preyy)+=Unew(predd,preyy);
      if(i==nyrs)
        UnewAVG(predd,preyy)=UnewAVG(predd,preyy)/nyrs;
    }//end prey  
  }// end pred 
FUNCTION void CALC_SUITOLD(int msmMode_pass_number)
  if(msmMode_pass_number>0){
    CALC_RATION(0);
   // cout<<" *******FIRST PART  =U/(W*N)*******"<<endl;
    dvariable suit_tmp ;
    // stom_div_bio.initialize();
    stom_div_bio2.initialize();
    if(msmMode_pass_number==1)
    {
      // for (int year=1;year<=2;year++) // year
      // {
      //   for (int pred=1;pred<=nspp;pred++) //predator species
      //   { 
      //     for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age 
      //     {
      //     for (int prey=1;prey<=nspp;prey++) // prey species 
      //     {
      //       for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp_age 
      //       {  
      //        suit_tmp = stom(year,pred,pred_age,prey,prey_age)/(AvgN(prey,stomyrs(year),prey_age)); 
      //        if (mn_wt_stom(pred,pred_age,prey,prey_age)!=0.)
      //          stom_div_bio(year,pred,pred_age,prey,prey_age) = suit_tmp/mn_wt_stom(pred,pred_age,prey,prey_age);    
      //         }  //end prey sp_age 
      //       }  // end  prey sp
      //     }  // end pred sp_age    
      //   } // end  pred species 
      // } // end  year loop  
    }else{
      for (int year=1;year<=nyrs;year++) // year
      {
      for (int pred=1;pred<=nspp;pred++) //predator species
      { 
        for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age 
        {
        for (int prey=1;prey<=nspp;prey++) // prey species 
        {
          for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp_age 
          {  
           suit_tmp = stomKir(year,pred,pred_age,prey,prey_age)/(AvgN(prey,year,prey_age)); 
           if (wt(prey,year,prey_age)!=0.)
             stom_div_bio2(year,pred,pred_age,prey,prey_age) = suit_tmp/wt(prey,year,prey_age);    
  
  //         if(AvgN(prey,year,prey_age)!=0)
  //           stom_div_bio2(year,pred,pred_age,prey,prey_age) = (UobsAge(pred,prey,pred_age,prey_age)*ration2Age(pred,year,pred_age))/(AvgN(prey,year,prey_age)*wt(prey,year,prey_age));    
          }  //end prey sp_age 
        }  // end  prey sp
        }  // end pred sp_age    
      } // end  pred species 
      } // end  year loop  
   }
    //end if kirs_pass_num==0
    // cout<<"************suitabilities calculation*******"<<endl;
    suit_main.initialize(); // sets to zero
    dvariable suma_suit;
    for (int pred=1;pred<=nspp;pred++) //predator species
    {
    for (int pred_age=1;pred_age<=nages(pred);pred_age++) // predator sp 
    {
      for (int prey=1;prey<=nspp;prey++) // prey specie
      {
      for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp
      { 
        switch (msmMode_pass_number)
        {
          case 0:
          break;
          case 1:
            // for (int year=1;year<=nstom;year++) // year
            // {
            //   suma_suit=sum(stom_div_bio(year,pred,pred_age)); // sum of all prey and prey ages in the stom of the pred sp_age j
            //   suit_std(year,pred,pred_age,prey,prey_age) = stom_div_bio(year,pred,pred_age,prey,prey_age)/(suma_suit+of_stom(year,pred,pred_age));  //(sum(stom_div_bio(year,pred,pred_age)+of_stom(year,pred,pred_age)));            
            //   suit_main(pred,pred_age,prey,prey_age)    += suit_std(year,pred,pred_age,prey,prey_age);                                  
            // }// end year
            // suit_main(pred,pred_age,prey,prey_age) /= nstom; 
            cout<<"!!!! THIS MODE NO LONGER USED (MSM MODE = 1) "<<endl;
          break;
          case 2:
            // add suit_std
            suit_main(pred,pred_age,prey,prey_age)=0;
            for (int year=1;year<=nyrs;year++) // year
            {
              suma_suit=sum(stom_div_bio2(year,pred,pred_age)); // sum of all prey and prey ages in the stom of the pred sp_age j
              suit_main(pred,pred_age,prey,prey_age)    += stom_div_bio2(year,pred,pred_age,prey,prey_age)/(suma_suit+of_stomKir(pred,pred_age,year));                               
            }// end year
            suit_main(pred,pred_age,prey,prey_age) /= nyrs; 
            
          break;
          //removed A4
          default:
            suit_main(pred,pred_age,prey,prey_age)=0;
            for (int year=1;year<=nyrs;year++) // year
            {
              suma_suit=sum(stom_div_bio2(year,pred,pred_age)); // sum of all prey and prey ages in the stom of the pred sp_age j
              suit_main(pred,pred_age,prey,prey_age)    += stom_div_bio2(year,pred,pred_age,prey,prey_age)/(suma_suit+of_stomKir(pred,pred_age,year));                               
            }// end year
            suit_main(pred,pred_age,prey,prey_age) /= nyrs; 
          break;
        }  //end switch              
      }// prey sp_age 
      } // prey
      suit_other(pred,pred_age) = 1. - sum(suit_main(pred,pred_age)); //estimate the other food suitability        
    }// pred sp_age
    }  // pred
  }//if msmMode >0
FUNCTION void CALC_SUITAGE_KIR(int pass_number)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// calculate suitability using jesus's approach and Kir's Umatrix
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  dvar_matrix Spaij(1,nspp,1,maxA); // suitability
  dvar_matrix Upaij(1,nspp,1,maxA); // suitability
  dvariable Spaij_tot;
  for(int predd=1;predd<=nspp;predd++){
    for(int pred_age=1;pred_age<=nages(predd);pred_age++){
      Spaij_tot.initialize();      
      Spaij.initialize();
      Upaij.initialize();
      for(int preyy=1;preyy<=nspp;preyy++){    
        dvar_vector Btmp(1,nages(preyy)); // Prey biomass
        Btmp.initialize();  
        if(pass_number==0)
          Btmp=elem_prod(AvgN(preyy,i),wt(preyy,i));
        if(pass_number==1)
          Btmp=elem_prod(AvgN_fut(preyy,i),wt_fut(preyy,i));
        //cout<<"Btmp"<<Btmp<<endl;
        for(int prey_age=1;prey_age<=nages(preyy);prey_age++){
          if(Btmp(prey_age)>0.){
            Spaij(preyy,prey_age)=UnewAge(predd,preyy,pred_age,prey_age)/Btmp(prey_age);
          }else{
            Spaij(preyy,prey_age)=0.;
          }  
          Spaij_tot+=Spaij(preyy,prey_age);
        }//end prey_age
      }  // end preyy
      //cout<<"Spaij"<<Spaij<<endl;
      if(Spaij_tot>0.){
        for(int preyy=1;preyy<=nspp;preyy++)
          for(int prey_age=1;prey_age<=nages(preyy);prey_age++)
            SnewAge(predd,preyy,pred_age,prey_age)=Spaij(preyy,prey_age)/Spaij_tot;
      }else{
        for(int preyy=1;preyy<=nspp;preyy++)
          for(int prey_age=1;prey_age<=nages(preyy);prey_age++)
            SnewAge(predd,preyy,pred_age,prey_age)=0.;
      }
    }// end
  }  

FUNCTION void CALC_SUIT_KIR(int pass_number)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// calculate suitability using jesus's approach and Kir's Umatrix
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  dvar_matrix Spaij(1,nspp,1,maxL); // suitability
  dvar_matrix Upaij(1,nspp,1,maxL); // suitability
  dvariable Spaij_tot;
  for(int predd=1;predd<=nspp;predd++){
    for(int predL=1;predL<=nlengths(predd);predL++){
      Spaij_tot.initialize();
      for(int preyy=1;preyy<=nspp;preyy++){        
        Spaij.initialize();
        Upaij.initialize();
        for(int preyL=1;preyL<=nlengths(preyy);preyL++){
          if(pass_number==0)
            Spaij(preyy,preyL)=Unew(predd,preyy,predL,preyL)/(AvgN2(preyy,i,preyL)*WtbyL(preyy,preyL));
          if(pass_number==1)
            Spaij(preyy,preyL)=Unew(predd,preyy,predL,preyL)/(AvgN2_fut(preyy,i,preyL)*WtbyL(preyy,preyL));
          Spaij_tot+=Spaij(preyy,preyL);
        }//end prey_age
      }  // end preyy
      for(int preyy=1;preyy<=nspp;preyy++)
        for(int preyL=1;preyL<=nlengths(preyy);preyL++)
          Snew(predd,preyy,predL,preyL)=Spaij(preyy,preyL)/Spaij_tot;
    }// end
  }  

FUNCTION void CALC_M2NEW(int M2new_number)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// LASTLY!! Calculate M2 as E/N for each sp_age
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  M2_pen.initialize();
  for(int preyy=1;preyy<=nspp;preyy++)
  {
    dvar_vector M2tmp(1,nages(preyy));
    dvar_vector availB2(1,nlengths(preyy));
    dvar_matrix num_eaten(1,nages(preyy),1,nlengths(preyy));
    dvar_vector Ftmp2(1,nages(preyy));
    dvar_vector Ztmp2(1,nages(preyy));
    dvar_vector include(1,nages(preyy));
    
    include.initialize();
    Ftmp2.initialize();
    M2tmp.initialize();
    availB2.initialize();
    num_eaten.initialize();
    if(M2new_number==0)
    {  
      Ftmp2=fsh_sel(preyy) * mean_F(preyy) * mfexp(F_dev(preyy,i));
      for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
        if(N(preyy,i,sp_age)>1)
          include(sp_age)=1;
    }
    if(M2new_number==1)
    {  
      Ftmp2=fsh_sel(preyy) * mean_F(preyy) * mfexp(F_dev(preyy,nyrs));
      for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
        if(N_fut(preyy,i,sp_age)>1)
          include(sp_age)=1;
    }
    
    
    //int repn2=20;
    switch(M2mode)
    {
      case 1:
      //M2tmp=M2(preyy,i);
      if(M2new_number==0)
      {
          M2tmp=elem_prod(M2(preyy,i),include);
          Ztmp2.initialize();
          Eage_hat(preyy,i).initialize();
          Ztmp2= M1(preyy) + M2tmp + Ftmp2(preyy) ;
          Eage_hat(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage(preyy,i)));
          Eage_hat(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage(preyy,i))); 
      }
      if(M2new_number ==1)
      {
          M2tmp=elem_prod(M2_fut(preyy,i),include);
          Ztmp2.initialize();
          Eage_hat(preyy,i).initialize();
          Ztmp2= M1(preyy) + M2tmp + Ftmp2(preyy) ;
          Eage_hat_fut(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage_fut(preyy,i)));
          Eage_hat_fut(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage_fut(preyy,i))); 
      }  
    break;
    case 3:
      if(M2new_number==0)
      {
        for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
          M2(preyy,i,sp_age)=pow(double(sp_age)*0.8,-2);
        M2tmp=M2(preyy,i);  
        Ztmp2= M1(preyy) + M2tmp + Ftmp2(preyy) ;
        
        Eage_hat(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage(preyy,i)));
        //if(current_phase()>3){
        for (int rep=1;rep<=repn2;rep++)
        {
          dvar_vector Mhat(1,nages(preyy));
          dvar_vector Mobs(1,nages(preyy));
          Ztmp2.initialize();
          Mhat.initialize();
          Mobs.initialize();
          Eage_hat(preyy,i).initialize();
          Ztmp2= M1(preyy) + M2tmp + Ftmp2(preyy) ;
          Eage_hat(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage(preyy,i)));
          Mhat=elem_div(elem_prod(Eage_hat(preyy,i),include),Bage(preyy,i));
          Mobs=elem_div(Eage(preyy,i),Bage(preyy,i));
          //M2tmp-=elem_prod(M2tmp,elem_div(Mhat-Mobs,Mhat+Mobs));
          for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
          {
            //if(Mobs(sp_age)>2)
            //  Mobs(sp_age)=1.2;
            if((Mhat(sp_age)+Mobs(sp_age))>0.)
              M2tmp(sp_age)+=M2tmp(sp_age)*((Mobs(sp_age)-Mhat(sp_age))/(Mobs(sp_age)+Mhat(sp_age)));   //M2tmp(sp_age)-=M2tmp(sp_age)*(Mhat(sp_age)-Mobs(sp_age));
            else
              M2tmp(sp_age)=0.;    
            if(Bage(preyy,i,sp_age)<=0.)
              M2tmp(sp_age)=0.;  
          }
          M2(preyy,i)=M2tmp;
        }    
      }
      if(M2new_number==1)
      {
        for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
          M2_fut(preyy,i,sp_age)=pow(double(sp_age)*0.8,-2);
        M2tmp=M2(preyy,i);  
        Ztmp2= M1(preyy) + M2tmp + Ftmp2(preyy) ;
        
        Eage_hat_fut(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage_fut(preyy,i)));
        for (int rep=1;rep<=repn2;rep++)
        {
          dvar_vector Mhat(1,nages(preyy));
          dvar_vector Mobs(1,nages(preyy));
          Ztmp2.initialize();
          Mhat.initialize();
          Mobs.initialize();
          Eage_hat(preyy,i).initialize();
          Ztmp2= M1(preyy) + M2tmp + Ftmp2(preyy) ;
          Eage_hat_fut(preyy,i) = elem_prod(elem_div(M2tmp,Ztmp2) ,elem_prod(1.-mfexp(-Ztmp2) , Bage_fut(preyy,i)));
          Mhat=elem_div(elem_prod(Eage_hat_fut(preyy,i),include),Bage_fut(preyy,i));
          Mobs=elem_div(Eage_fut(preyy,i),Bage_fut(preyy,i));
          for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
          {
            if((Mhat(sp_age)+Mobs(sp_age))>0.)
              M2tmp(sp_age)+=M2tmp(sp_age)*((Mobs(sp_age)-Mhat(sp_age))/(Mobs(sp_age)+Mhat(sp_age)));   //M2tmp(sp_age)-=M2tmp(sp_age)*(Mhat(sp_age)-Mobs(sp_age));
            else
              M2tmp(sp_age)=0.;    
            if(Bage(preyy,i,sp_age)<=0)
              M2tmp(sp_age)=0.;  
          }
          M2_fut(preyy,i)=M2tmp;
        }        
      }  //end M2new_number
    break;
    case 0:
      if(M2new_number==0){
        for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
          M2(preyy,i,sp_age)=pow(double(sp_age)*0.8,-2);
        for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
        {
          if(Bage(preyy,i,sp_age)>0.){
            M2(preyy,i,sp_age)=Eage(preyy,i,sp_age)/Bage(preyy,i,sp_age);
          }else{
            M2(preyy,i,sp_age)=0.;
          }
        }
        Eage_hat(preyy,i)=elem_prod(Bage(preyy,i),M2(preyy,i));
      }
      if(M2new_number==1){
          for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
            M2_fut(preyy,i,sp_age)=pow(double(sp_age)*0.8,-2);
          for(int sp_age=1;sp_age<=nages(preyy);sp_age++)
          {
            if(Bage_fut(preyy,i,sp_age)>0.){
              M2_fut(preyy,i,sp_age)=Eage_fut(preyy,i,sp_age)/Bage_fut(preyy,i,sp_age);
            }else{
              M2_fut(preyy,i,sp_age)=0.;
            }
          }
          Eage_hat_fut(preyy,i)=elem_prod(Bage_fut(preyy,i),M2_fut(preyy,i));
      }  
    break;
    default:
    break;
    case 4:
    break;  
    }// end switch(M2mode)
  }//end for each preyy
  
  if(M2mode==4)
  {
    if(RationByAge==0){
      cout<<"ERROR - M2mode=4 can only be used if RationByAge=0!!!! Please adjust .dat file"<<endl;
      exit(1);
    }  
    dvar_matrix sumofSnew (1,nspp,2,nages);
    dvar_matrix suit_otherKir (1,nspp,2,nages);
    suit_otherKir.initialize();
    suit_mainKir.initialize();
    if(M2new_number==0)
    {
      CALC_SUITAGE_KIR(0);
      for (int predd=1;predd<=nspp;predd++) //predator species
        for (int pred_age=1;pred_age<=nages(predd);pred_age++) // predator sp 
          for (int preyy=1;preyy<=nspp;preyy++) // prey specie
            for (int prey_age=1;prey_age<=nages(preyy);prey_age++) // prey sp
              suit_mainKir(predd,pred_age,preyy,prey_age)=SnewAge(predd,preyy,pred_age,prey_age);
      for (int pred=1;pred<=nspp;pred++) //predator species
        for (int pred_age=1;pred_age<=nages(pred);pred_age++) // predator sp 
           suit_otherKir(pred,pred_age) = 1. - sum(suit_mainKir(pred,pred_age)); //estimate the other food suitability        
      
      dvariable M2top;
      dvariable M2bottom;
      for(int preyy=1;preyy<=nspp;preyy++)
        for(int prey_age=1;prey_age<=nages(preyy);prey_age++)
          for(int predd=1;predd<=nspp;predd++)
            for(int pred_age=1;pred_age<=nages(predd);pred_age++){
              M2top.initialize();
              M2bottom.initialize();
              M2top=(AvgN(predd,i,pred_age)*ration2Age(predd,i,pred_age)*SnewAge(predd,preyy,pred_age,prey_age));
              M2bottom=(AvgN(preyy,i,prey_age)*wt(preyy,i,prey_age)*SnewAge(predd,preyy,pred_age,prey_age)+(other_food(predd)*suit_otherKir(predd,pred_age)));
              if(M2bottom>0.)
                M2(preyy,i,prey_age)+=M2top/M2bottom;
            }  
      for(int preyy=1;preyy<=nspp;preyy++){
        dvar_vector Ftmp2(1,nages(preyy));
        dvar_vector Ztmp2(1,nages(preyy));
        dvar_vector M2tmp2(1,nages(preyy));
        dvar_vector Btmp2(1,nages(preyy));
        Ftmp2.initialize();
        Ztmp2.initialize();
        M2tmp2.initialize();
        Btmp2.initialize();
        M2tmp2= M2(preyy,i);
        Btmp2=elem_prod(N(preyy,i),wt(preyy,i));
        Ztmp2= M1(preyy) + M2tmp2 + Ftmp2(preyy) ;
        Eage_hat(preyy,i) = elem_prod(elem_div(M2tmp2,Ztmp2),elem_prod(1.-mfexp(-Ztmp2) ,Btmp2));
      }  
    }
    if(M2new_number==1)
    {
      CALC_SUITAGE_KIR(1);
      
      for (int predd=1;predd<=nspp;predd++) //predator species
        for (int pred_age=1;pred_age<=nages(predd);pred_age++) // predator sp 
          for (int preyy=1;preyy<=nspp;preyy++) // prey specie
            for (int prey_age=1;prey_age<=nages(preyy);prey_age++) // prey sp
              suit_mainKir(predd,pred_age,preyy,prey_age)=SnewAge(predd,preyy,pred_age,prey_age);
      for (int pred=1;pred<=nspp;pred++) //predator species
        for (int pred_age=1;pred_age<=nages(pred);pred_age++) // predator sp 
           suit_otherKir(pred,pred_age) = 1. - sum(suit_mainKir(pred,pred_age)); //estimate the other food suitability        
      
      dvariable M2top;
      dvariable M2bottom;
      for(int preyy=1;preyy<=nspp;preyy++)
        for(int prey_age=1;prey_age<=nages(preyy);prey_age++)
          for(int predd=1;predd<=nspp;predd++)
            for(int pred_age=1;pred_age<=nages(predd);pred_age++){
              M2top.initialize();
              M2bottom.initialize();
              M2top=(AvgN_fut(predd,i,pred_age)*ration2Age_fut(predd,i,pred_age)*SnewAge(predd,preyy,pred_age,prey_age));
              M2bottom=(AvgN_fut(preyy,i,prey_age)*wt_fut(preyy,i,prey_age)*SnewAge(predd,preyy,pred_age,prey_age)+(other_food(predd)*suit_otherKir(predd,pred_age)));
              if(M2bottom>0.)
                M2_fut(preyy,i,prey_age)+=M2top/M2bottom;
            }  
      for(int preyy=1;preyy<=nspp;preyy++){
        dvar_vector Ftmp2(1,nages(preyy));
        dvar_vector Ztmp2(1,nages(preyy));
        dvar_vector M2tmp2(1,nages(preyy));
        dvar_vector Btmp2(1,nages(preyy));
        Ftmp2.initialize();
        Ztmp2.initialize();
        M2tmp2.initialize();
        Btmp2.initialize();
        M2tmp2= M2_fut(preyy,i);
        Btmp2=elem_prod(N_fut(preyy,i),wt_fut(preyy,i));
        Ztmp2= M1(preyy) + M2tmp2 + Ftmp2(preyy) ;
        Eage_hat_fut(preyy,i) = elem_prod(elem_div(M2tmp2,Ztmp2),elem_prod(1.-mfexp(-Ztmp2) ,Btmp2));
      }  
      
    }  // end if(M2new_number==0)
  }// end if M2mode==4
FUNCTION void calc_B0(int output)  
FUNCTION calc_B_fut  
FUNCTION void run_control_rules(int output2)  
FUNCTION void control_rule_1(int sub)
FUNCTION void control_rule_2(int sub)
FUNCTION void control_rule_3(int sub)    
FUNCTION void run_control_rules_old(int output2)  
FUNCTION void Do_Projections(int output)


FUNCTION void WRITE_CR_PROJECTION_REPORT(dvariable cmult_pass1)
    if(cmult_pass1==0){ 
      CR_projection_report << "Scenario "<<"Control_rule "<<"species "<<"future_year "<<"TempC_futUSE(itemp1,i) "<<"AvgN_futKEEP "
            <<"mean_ration2Age_futKEEP "<<"RationKeep_fut "<<"RationScaledKeep_fut "<< "SSB_RSkeep_fut "<<"R_fut"<<endl;
     }else{ 
      for(int itemp1=1;itemp1<=ntemp_scen;itemp1++)
        for(int sp=1;sp<=nspp;sp++)
          for(int yr=1;yr<=nyrs_fut;yr++)
             CR_projection_report<<" "<<itemp1<<" "<<cmult_pass1/10<<" "<<sp<<" "<<yr<<" "<<TempC_futUSE(itemp1,yr)
           <<" "<<AvgN_futKEEP(itemp1,sp,yr)<<" "<<mean(ration2Age_futKEEP(itemp1,sp,yr))<<" "<<RationKeep_fut(itemp1,sp,yr)
           <<" "<<RationScaledKeep_fut(itemp1,sp,yr)<<" "<<SSB_RSkeep_fut(itemp1,sp,yr)<<" "<<R_fut(itemp1,sp,yr)<<endl;
    }
 
FUNCTION void WRITE_PROJECTION_REPORT(dvariable cmult_pass1)
  // write report for age specific Biomass, Numbers, Catch, M2, M1, and catch
  cout<<"write projection report"<<endl;
  if(cmult_pass1==0){ 
     R_projection_report 
     <<"harvestMode "
     << "Scenario "
     <<"Control_rule "
     <<"species "
     <<"age "
     <<"future_year "
     <<"bottomT_C "
     <<"wt_at_age "
     <<"ABC_total_biom "
     <<"Catch_total_biom "
     <<"Catch_Num_byAge "
     <<"pmature "
     <<"B "
     <<"B0 "
     <<"SSB "
     <<"SSB0 "
     <<"SSB_total_biom "
     <<"SSB0_total_biom "
     <<"N "
     <<"N0 "
     <<"M1 "
     <<"M2 "
     <<"F40 "
     <<"F35 "
     <<"Fabc "
     <<"F "
     <<"F_age "
     <<"Ration_g_fish"<<endl;
  }else{ 
    if(harvestMode==3)
    {
      for(int itemp1=1;itemp1<=ntemp_scen;itemp1++)
      {
        for(int sp=1;sp<=nspp;sp++)
        {
          // for estimation years
          for(int yr=1;yr<=nyrs;yr++)
              for(int sp_age=1;sp_age<=nages(sp);sp_age++)
                R_projection_report <<
                  harvestMode<<" "
                  <<itemp1<<" "
                  <<77.77<<" "
                  <<sp<<" "
                  <<sp_age<<" "
                  <<yr<<" "
                  <<TempC(yr)<<" "
                  <<wt(sp,yr,sp_age)<<" "
                  <<0<<" "
                  <<tc_biom_hat(sp,yr)<<" "
                  <<catch_hat(sp,yr,sp_age)<<" "
                  <<pmature(sp,sp_age)<<" "
                  <<biomassByage_hind2fut(itemp1,sp,yr,sp_age)<<" "
                  <<biomassByage0_hind2fut(itemp1,sp,yr,sp_age)<<" "
                  <<biomassByage_hind2fut(itemp1,sp,yr,sp_age)*pmature(sp,sp_age)<<" "
                  <<biomassByage0_hind2fut(itemp1,sp,yr,sp_age)*pmature(sp,sp_age)<<" "
                  <<biomass_hind2fut(itemp1,sp,yr)<<" "
                  <<biomass0_hind2fut(itemp1,sp,yr)<<" "
                  <<NByage_hind2fut(itemp1,sp,yr,sp_age)<<" "
                  <<NByage0_hind2fut(itemp1,sp,yr,sp_age)<<" "
                  <<M1(sp,sp_age)<< " "
                  <<M2(sp,yr,sp_age)<< " "
                  <<0<<" "
                  <<0<<" "
                  <<0<<" "
                  <<mean_F(sp) * mfexp(F_dev(sp,yr))<<" "
                  <<F(sp,yr,sp_age)<<" "
                  <<ration2Age(sp,yr,sp_age)<<endl;

          //now for future years
          for(int yr=1;yr<=nyrs_fut;yr++)
              for(int sp_age=1;sp_age<=nages(sp);sp_age++)
               R_projection_report <<
                  harvestMode<<" "
                  << itemp1<<" "
                  <<(cmult_pass1/10)<<" "
                  <<sp<<" "
                  <<sp_age<<" "
                  <<yr+nyrs<<" "
                  <<TempC_futUSE(itemp1,yr)<<" "
                  <<wt_fut_all(itemp1,sp,yr,sp_age)<<" "
                  <<ABC_fut(itemp1,sp,yr)<<" "
                  <<tc_biom_hat_fut(itemp1,sp,yr)<<" "
                  <<catch_hat_fut(itemp1,sp,yr,sp_age)<<" "
                  <<pmature(sp,sp_age)<<" "
                  <<biomassByage_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                  <<biomassByage0_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                  <<biomassByage_hind2fut(itemp1,sp,yr+nyrs,sp_age)*pmature(sp,sp_age)<<" "
                  <<biomassByage0_hind2fut(itemp1,sp,yr+nyrs,sp_age)*pmature(sp,sp_age)<<" "
                  <<biomass_hind2fut(itemp1,sp,yr+nyrs)<<" "
                  <<biomass0_hind2fut(itemp1,sp,yr+nyrs)<<" "
                  <<NByage_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                  <<NByage0_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                  <<M1(sp,sp_age)<< " "
                  <<M2_fut_all(itemp1,sp,yr,sp_age)<< " "
                  <<F40(itemp1,sp,yr)<<" "
                  <<F35(itemp1,sp,yr)<<" "
                  <<Fabc(itemp1,sp,yr)<<" "
                  <<PF(itemp1,sp,yr)<<" "
                  <<F_fut_all(itemp1,sp,yr,sp_age)<<" "
                  <<ration2Age_futKEEP(itemp1,sp,yr,sp_age)<<endl;
          }
        }
      }
      if(harvestMode!=3)
      {
        for(int itemp1=1;itemp1<=ntemp_scen;itemp1++)
        {
          for(int sp=1;sp<=nspp;sp++)
          {
            // for estimation years
            for(int yr=1;yr<=nyrs;yr++)
                for(int sp_age=1;sp_age<=nages(sp);sp_age++)
                  R_projection_report <<
                    harvestMode<<" "
                    <<itemp1<<" "
                    <<77.77<<" "
                    <<sp<<" "
                    <<sp_age<<" "
                    <<yr<<" "
                    <<TempC(yr)<<" "
                    <<wt(sp,yr,sp_age)<<" "
                    <<0<<" "
                    <<tc_biom_hat(sp,yr)<<" "
                    <<catch_hat(sp,yr,sp_age)<<" "
                    <<pmature(sp,sp_age)<<" "
                    <<biomassByage_hind2fut(itemp1,sp,yr,sp_age)<<" "
                    <<biomassByage0_hind2fut(itemp1,sp,yr,sp_age)<<" "
                    <<biomassByage_hind2fut(itemp1,sp,yr,sp_age)*pmature(sp,sp_age)<<" "
                    <<biomassByage0_hind2fut(itemp1,sp,yr,sp_age)*pmature(sp,sp_age)<<" "
                    <<biomass_hind2fut(itemp1,sp,yr)<<" "
                    <<biomass0_hind2fut(itemp1,sp,yr)<<" "
                    <<NByage_hind2fut(itemp1,sp,yr,sp_age)<<" "
                    <<NByage0_hind2fut(itemp1,sp,yr,sp_age)<<" "
                    <<M1(sp,sp_age)<< " "
                    <<M2(sp,yr,sp_age)<< " "
                    <<0<<" "
                    <<0<<" "
                    <<0<<" "
                    <<mean_F(sp) * mfexp(F_dev(sp,yr))<<" "
                    <<F(sp,yr,sp_age)<<" "
                    <<ration2Age(sp,yr,sp_age)<<endl;

            // for projection years
            for(int yr=1;yr<=nyrs_fut;yr++)
                for(int sp_age=1;sp_age<=nages(sp);sp_age++)
                  R_projection_report <<
                    harvestMode<<" "
                    <<itemp1<<" "
                    <<77.77<<" "
                    <<sp<<" "
                    <<sp_age<<" "
                    <<yr+nyrs<<" "
                    <<TempC_futUSE(itemp1,yr)<<" "
                    <<wt_fut_all(itemp1,sp,yr,sp_age)<<" "
                    <<ABC_fut(itemp1,sp,yr)<<" "
                    <<tc_biom_hat_fut(itemp1,sp,yr)<<" "
                    <<catch_hat_fut(itemp1,sp,yr,sp_age)<<" "
                    <<pmature(sp,sp_age)<<" "
                    <<biomassByage_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                    <<biomassByage0_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                    <<biomassByage_hind2fut(itemp1,sp,yr+nyrs,sp_age)*pmature(sp,sp_age)<<" "
                    <<biomassByage0_hind2fut(itemp1,sp,yr+nyrs,sp_age)*pmature(sp,sp_age)<<" "
                    <<biomass_hind2fut(itemp1,sp,yr+nyrs)<<" "
                    <<biomass0_hind2fut(itemp1,sp,yr+nyrs)<<" "
                    <<NByage_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                    <<NByage0_hind2fut(itemp1,sp,yr+nyrs,sp_age)<<" "
                    <<M1(sp,sp_age)<< " "
                    <<M2_fut_all(itemp1,sp,yr,sp_age)<< " "
                    <<F40(itemp1,sp,yr)<<" "
                    <<F35(itemp1,sp,yr)<<" "
                    <<Fabc(itemp1,sp,yr)<<" "
                    <<PF(itemp1,sp,yr)<<" "
                    <<F_fut_all(itemp1,sp,yr,sp_age)<<" "
                    <<ration2Age_futKEEP(itemp1,sp,yr,sp_age)<<endl;
          }
        }
      }
    }

    // R_projection_report<<" "<<itemp1<<" "<<77.77<<" "<<sp<<" "<<yr<<" "<<TempC_futUSE(itemp1,yr)<<" "<<AvgN_futKEEP(itemp1,sp,yr)<<" "<<ration2Age_futKEEP(itemp1,sp,yr)<<" "<<RationKeep_fut(itemp1,sp,yr)<<" "<<RationScaledKeep_fut(itemp1,sp,yr)<<" "<<SSB_RSkeep_fut(itemp1,sp,yr)<<" "<<R_fut(itemp1,sp,yr)<<endl;
                
    // R_projection_report.close();
        //  age = "all" 
        // tc_biom_hat_fut(itemp,sp,yr)
        // biomassSSB_hind2fut(itemp,sp,yr)
        // biomass_hind2fut(itemp,sp,yr)
        // biomass0_hind2fut(itemp,sp,yr)



// FUNCTION dvar_vector SolveF_Kir(const dvar_vector& Ctarget, int spn)
  /*
    // gradient decent method
    RETURN_ARRAYS_INCREMENT();
    dvar_matrix   tmpZ(1,nspp,1,nages);
    dvariable     adjn;
    dvar_vector   Chatt(1,nspp);
    dvar_vector   tmpFF(1,nspp);
    tmpFF=1.;
    if(i==1){
      tmpFF(spn)=mfexp(F_dev(spn,nyrs)+ln_mean_F(spn));
    }else{
      tmpFF(spn)=PF(itemp,spn,i-1);
    }
    Chatt.initialize();
    adjn=1e-8;
    for(int rep=1;rep<=10;rep++)
    {

        dvar_vector Fatmp=fsh_sel(spn)*tmpFF(spn);
        tmpZ(spn) =Fatmp+M1(spn)(1,nages(spn)) + M2_fut(spn,i);
        Chatt(spn)  = (elem_prod( elem_div( Fatmp,tmpZ(spn)) ,elem_prod( 1.-mfexp(-tmpZ(spn)) , N_fut(spn,i))) )* wt_fut(spn,i)(1,nages(spn));
        if((Chatt(spn)+Ctarget(spn))>0)
         tmpFF(spn)+=tmpFF(spn)*((Ctarget(spn)-Chatt(spn))/(Ctarget(spn)+Chatt(spn))); 
         else
          tmpFF(spn)    =adjn;   

    }
    if(tmpFF(spn)<0)
      tmpFF(spn)=adjn; 
    RETURN_ARRAYS_DECREMENT();
    return(tmpFF); 
  */

//==============================================================================
//// TRASH THESE
//==============================================================================

// FUNCTION void CALC_REC_FUT_OLD(int recMode_pass)

  //   dvariable Eating;
  //   //dvariable noRation;
  //   //noRation=1;
  //   dvariable covars;
  //   dvariable SSB_rs_fut;

  //   for (int sp=1;sp<=nspp;sp++)
  //   { 
  //     switch (recMode_pass)
  //     {
  //       case 0:  //          0 = project under mean rec
  //        if(rand_rec){  
  //          //          1 = project under random rec (based on estimates of sd rec)
  //           rec_sigma(sp)=sum((rec_dev(sp)(1,nyrs) ))/nyrs;
  //           R_fut(itemp,sp,i) = mfexp(sqrt(norm2(rec_dev(sp)(1,nyrs))/nyrs)*rmult(sp,i)+ln_mn_rec(sp));
  //         }else{
  //           R_fut(itemp,sp,i) = sum(mfexp(ln_mn_rec(sp) + rec_dev(sp)(1,nyrs) ))/nyrs;
  //         }
  //       break; 
  //       default:
  //         Eating=0.;
  //         covars=0.;
  //         SSB_rs_fut=0.;
  //         if(i==1){
  //            SSB_rs_fut = elem_prod(N(sp,nyrs),pmature(sp))*wt(sp,nyrs)(1,nages(sp));  
  //            Eating = ((AvgN(sp,nyrs,1)*ration2Age(sp,nyrs,1))/rs_cov(fallZ_num,nyrs)); 
  //           if(noRation(sp))
  //               Eating = rs_cov(fallZ_num,nyrs); 
  //           if(rand_rec){
  //             for (int cc=3;cc<=ncov;cc++)
  //                 covars+=(rs_parm_std(cc,sp)*rmult(sp,i)+rs_parm(cc,sp) )*rs_cov(cc,nyrs); //rs_parm(cc,sp)*rs_cov(cc,nyrs);
  //             R_fut(itemp,sp,i) = mfexp(log(aa_c(sp)*SSB_rs_fut) -bb_c(sp)*SSB_rs_fut+((rs_parm_std(1,sp)*rmult(sp,i)+rs_parm(1,sp) )*rs_cov(1,nyrs))-  (( rs_parm_std(fallZ_num,sp)*rmult(sp,i)+rs_parm(fallZ_num,sp) )*Eating)+covars);
  //           }else{  
  //             for (int cc=3;cc<=ncov;cc++)
  //               covars+=rs_parm(cc,sp)*rs_cov(cc,nyrs);
  //             R_fut(itemp,sp,i) = mfexp(log(aa_c(sp)*SSB_rs_fut) -bb_c(sp)*SSB_rs_fut+rs_parm(1,sp)*rs_cov(1,nyrs) - rs_parm(fallZ_num,sp)*Eating+covars);
  //           }
  //         }else{
  //           SSB_rs_fut = elem_prod(N_fut(sp,i-1),pmature(sp))*wt_fut(sp,i-1)(1,nages(sp));
  //           Eating = ((AvgN_fut(sp,i-1,1)*ration2Age_fut(sp,i-1,1))/rs_cov_fut(itemp,fallZ_num,i-1));
  //           if(noRation(sp))
  //               Eating = rs_cov_fut(itemp,2,i-1); 
  //           if(rand_rec){
  //             for (int cc=3;cc<=ncov;cc++)
  //               covars+=(rs_parm_std(cc,sp)*rmult(sp,i)+rs_parm(cc,sp) )*rs_cov_fut(itemp,cc,i-1); //rs_parm(cc,sp)*rs_cov(cc,nyrs);
  //             R_fut(itemp,sp,i) = mfexp(log(aa_c(sp)*SSB_rs_fut) -bb_c(sp)*SSB_rs_fut+ (( rs_parm_std(1,sp)*rmult(sp,i)+rs_parm(1,sp) )*rs_cov_fut(itemp,1,i-1))- ( (rs_parm_std(fallZ_num,sp)*rmult(sp,i)+rs_parm(fallZ_num,sp)) *Eating)+covars);
  //           }else{  
  //             for (int cc=3;cc<=ncov;cc++)
  //                 covars+=rs_parm(cc,sp)*rs_cov_fut(itemp,cc,i-1);
  //               R_fut(itemp,sp,i) = mfexp(log(aa_c(sp)*SSB_rs_fut) -bb_c(sp)*SSB_rs_fut+rs_parm(1,sp)*rs_cov_fut(itemp,1,i-1) - rs_parm(fallZ_num,sp)*Eating+covars);
  //           }
  //         }
  //       break; 
  //     }
  //   }  


//==============================================================================
// REPORT SECTION
//==============================================================================
REPORT_SECTION
  cout <<endl<<"========End of phase: "<<current_phase()<<" ============"<<endl<<endl;
  //cout<<"in Report section now" <<endl  ;      
  report << "INPUTS"<<endl;
  report << "Catch at sp_age (obs_catch)"<<endl;
  report << obs_catch <<endl;
  report << "weight (wt)"<<endl;
  report << wt <<endl; 
  report << "ration"<<endl;//
  report << ration <<endl;//
  report << " avg weight"<<endl;
  report << wt <<endl;  
  report << "other food"<<endl;
  report << other_food <<endl;
  report << "suit plk as predator, plk as prey, sp_age of prey 1?"<<endl;
  report << suit_main(1) <<endl;
  report << "suit cod as predator, plk as prey, sp_age of prey 1?" << endl;
  report << suit_main(2) <<endl;
  report << "suitabilities other food"<<endl;
  report << suit_other <<endl;
  report << "residual mortality"<<endl;
  report << M1_base <<endl;  

  report << "OUTPUTS"<<endl;
  report << "Numbers at sp_age" <<endl;
  report << N<<endl<<endl;

  report << "Fishing mortality at sp_age" <<endl;
  report << F<<endl<<endl;

  report << "Survival at sp_age" <<endl;
  report << S<<endl<<endl;

  report << "Biomass"<<endl<<biomass<<endl;
  report << "Predicted total catch"<<endl;
  report <<  tc_hat <<endl;

  report << "observed total catch"<<endl;
  report <<  tc_obs <<endl;

  report << "Predicted catch at sp_age"<<endl;
  report <<  fsh_age_hat <<endl;

  report << "observed catch at sp_age"<<endl;
  report <<  fsh_age_obs <<endl;

  report << "Values for M2        "<<endl;
  report << "plk          "<<endl;
  report <<  M2(1) <<endl<<endl;
  report << "pcod          "<<endl;
  report <<  M2(2) <<endl<<endl;


  report << "available food plk"<<endl;
  report << avail_food(1)<< endl;
  report << "available food cod"<<endl;
  report << avail_food(2)<< endl;

  report << "average N plk"<<endl;
  report << AvgN(1)<< endl;
  report << "average N cod"<<endl;
  report << AvgN(2)<< endl;

  report << "======Survey_selectivity-at-sp_age==========="<<endl;
  for (int sp=1;sp<=nspp;sp++)
  {
    report << "==Species: "<<sp<<"----"<<endl;
    report << srv_sel(sp)<<" "<<endl;
  }

  report << "======Survey_Biomass_Fit=================="<<endl;
  for (int sp=1;sp<=nspp;sp++)
  {
    report << "==Species: "<<sp<<"----"<<endl;
    for (int yr=1;yr<=nyrs_srv_biom(sp);yr++)
      report << yrs_srv_biom(sp,yr)<<" "<<srv_biom(sp,yr)<< " "<<srv_biom_hat(sp,yr)<<endl;
  }
  report << "======Survey_age_composition_fits=================="<<endl;
  for (int sp=1;sp<=nspp;sp++)
  {
    report << "==Species: "<<sp<<"----"<<endl;
    for (int yr=1;yr<=nyrs_srv_age(sp);yr++)
      report << yrs_srv_age(sp,yr)<<" "<<srv_age_obs(sp,yr)<< " | | "<<srv_age_hat(sp,yr)<<endl;
  }
  report << "======Fishery_age_composition_fits=================="<<endl;
  for (int sp=1;sp<=nspp;sp++)
  {
    report << "==Species: "<<sp<<"----"<<endl;
    for (int yr=1;yr<=nyrs_fsh_comp(sp);yr++)
      report << styr+yr-1 <<" "<<fsh_age_obs(sp,yr)<< " | | "<<fsh_age_hat(sp,yr)<<endl;
  }
  report << "======Total_catch_fits=================="<<endl;
  for (int sp=1;sp<=nspp;sp++)
  {
    report << "==Species: "<<sp<<"----"<<endl;
    report << tc_obs(sp)<< endl<<tc_hat(sp)<<endl;
  }
  report << "eit_sel"<< eit_sel<<endl;
  report << "eit_sel(1,1)"<< eit_sel(1,1)<<endl;
  report << "eit_sel(1)"<< eit_sel(1)<<endl;
  report << "eit_hat"<< endl;
  report << eit_hat << endl;
    report << "==EIT_Survey: pollock:---"<<endl;
    for (int yr=1;yr<=n_eit;yr++)
      report << yrs_eit(yr)<<" "<<obs_eit_age(yr)<< " | | "<<eit_age_hat(yr)<<endl;
  report<<M1<<endl;
  report<<"predation mortality produce by plk"<<endl;
  report<<M2(1)<<endl;

  report << "obj_fun      "   <<endl;
  report << obj_fun            <<endl;
  report << "M2_like "   <<endl;
  report << M2_like       <<endl;
  report << "srv_sel_like "   <<endl;
  report << srv_sel_like       <<endl;
  report << "fsh_sel_like "   <<endl;
  report << fsh_sel_like       <<endl;

  report << "srv_bio_like "   <<endl;
  report << srv_bio_like       <<endl;
  report << "fsh_cat_like "   <<endl;
  report << fsh_cat_like       <<endl;
  report << "norm_srv_like"   <<endl;
  report << norm_srv_like      <<endl;
  report << "norm_fsh_like"   <<endl;
  report << norm_fsh_like      <<endl;

  report << "srv_age_like "   <<endl;
  report << srv_age_like       <<endl;
  report << "fsh_age_like "   <<endl;
  report << fsh_age_like       <<endl;
  report << "eit_age_like "   <<endl;
  report << eit_age_like       <<endl;


 /*
 FUNCTION  CONTROL_RULE_5(int sub)
// not currently used - was old SR approach for previous CR  dvariable& Btmp)
  Btmp.initialize();
  Rzero    =  mfexp(log_Rzero); 
  dvar_vector Ntmp(1,nages);
  dvar_vector survtmp = mfexp(-natmort);
  Ntmp.initialize();

  Ntmp(1) = Rzero;
  for ( j=1 ; j < nages; j++ )
  {
    Ntmp(j+1)  =   Ntmp(j) * survtmp(j); // Begin numbers in the next year/sp_age class
  }
  Ntmp(nages)  /= (1.- survtmp(nages)); 
  for ( j=1 ; j <= nages; j++ )
    Ntmp(j) *= mfexp(yrfrac*log(survtmp(j)));
    // Ntmp(j) *= pow(survtmp(j),yrfrac);

  Bzero = elem_prod(wt_ssb(endyr_r) , p_mature) * Ntmp ; // p_mature is of Females (half of adults)
  
  switch (SrType)
  {
    case 1:
      alpha = log(-4.*steepness/(steepness-1.)); // Eq. 13
      break;
    case 2:
      alpha  =  Bzero * (1. - (steepness - 0.2) / (0.8*steepness) ) / Rzero;
      beta   = (5. * steepness - 1.) / (4. * steepness * Rzero);
      break;
    case 4:
    //R = S * EXP(alpha - beta * S))
      beta  = log(5.*steepness)/(0.8*Bzero) ;
      alpha = log(Rzero/Bzero)+beta*Bzero;
      break;
  }
  */

BETWEEN_PHASES_SECTION

    if(dump_rep) cout<<" in Between... "<<current_phase()<<endl;
   UPDATE_BETWEEN();   
    if(dump_rep)cout<<" Moving on to phase... "<<current_phase()+1<<endl;
FINAL_SECTION
  cout<<"finish up"<<endl;
  if(dump_rep) {
    cout<<"SSB "<< biomassSSB<<endl;
    cout<<"____________"<<endl;
    cout<<"srv_sel "<<srv_sel<<endl;
    cout<<"____________"<<endl;
    cout<<"fsh_sel "<<fsh_sel<<endl;
    cout<<"____________"<<endl;
    cout<<"N(1) "<< N(1)<<endl;
  }
  if(do_fut || simMode>0){
    WRITE_FUT_REPORT(0);
    WRITE_FUT_REPORT(1);
    //WRITE_PROJECTION_REPORT;
  }  
  write_fut_files();
  write_RS_files();
  write_R();
  write_output();
  //  write_U();
  cout<<"*******************************************"<<endl;
  cout<<"*******************************************"<<endl;
  cout<<"         Model "<<flname<<" complete "<<endl;
  cout<<"*******************************************"<<endl;
  cout<<"*******************************************"<<endl;
RUNTIME_SECTION
  convergence_criteria .1, 0.001,0.00000001
  maximum_function_evaluations   50,  100,  1000,  10000
  //convergence_criteria       .1,  0.01,  0.01,  0.01,  0.01,   0.001,   0.001,  0.00000001
  //maximum_function_evaluations   50,  100,  100,  100,  100,  1000,  1000,  10000
TOP_OF_MAIN_SECTION
  arrmblsize = 100000000; 
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000); 
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);
GLOBALS_SECTION
 #include <math.h>
 #include <cmath>
 #include <admodel.h>
 #include <iostream>
 // adstring crulename;
 adstring flname;
 adstring foldername;
 adstring simname;
 adstring datafile_name;
 // adstring dietfile_name;
 adstring diet2file_name;
 //adstring lengthfile_name;
 adstring stomfile_name; 

 adstring Uobsfile_name;
 //adstring WtAgefile_name;
 adstring ATFfile_name;
 adstring Recfile_name;
 adstring retrofile_name;
 adstring futfile_name;
 adstring Fprofile_datafile;
 adstring catch_in_name;
 adstring setC_name;
 adstring setF_name;

 // currently functional
 ofstream write_input_log("dat_input.log");   
 ofstream R_report("results/ceattle_R2.rep");
 ofstream Wage_fut_report("results/Wage_fut.dat");
 ofstream future_report("results/Future_report.rep");
 ofstream covar_report("results/ceattle_recCovars.rep");
 ofstream Rec_report("results/ceattle_recruitment.rep");
 ofstream RecPar_report("results/ceattle_recPar.rep");
 ofstream CR_projection_report("results/ceattle_ContRule_projection.rep");
 ofstream R_projection_report("results/ceattle_R_projection.rep");
 
 // defunct?
 ofstream rs_out("results/ceattle_rs.dat");
 ofstream rand_report("results/randn.dat");
 ofstream future_predation_report("results/Future_Predation_report.rep");
 ofstream future_Fprofile_report("results/Future_Fprofile_report.rep");

 ofstream output_report("results/ceattle_output.dat");

 // ofstream U_report("results/ceattle_U.rep");
 // ofstream mcout("results/ceattle_mcmc.rep");
 // ofstream cross_val("results/ceattle_cross_val.rep");
   // outfile("", ios::append); // will append the datafile with the most recent run - be sure to add comment at start of file to seperate text string

 
  ofstream R_output("results/R_output.txt");   
  #undef log_input
  #define log_input(object) write_input_log << "# " #object "\n" << object << endl;
  #undef repR_4spp
  #define repR_4spp(object) R_output << "# " #object<<endl; for (int sp=1;sp<=nspp;sp++){          R_output << ration2Age(sp)<<endl;}



