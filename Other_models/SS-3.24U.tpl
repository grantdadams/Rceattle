// Stock Synthesis
// Developed by Richard Methot, NOAA Fisheries

//  ProgLabel_1.0_#  DATA_SECTION
DATA_SECTION

//  ProgLabel_1.1.1  Create string with version info
  !!version_info+="SS-V3.24U-safe;_08/29/2014;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_11.1x64";
  !!version_info_short+="#V3.24U";

//*********COUNTERS*************************
  int z // counters for size (length)
  int z1  // min for z counter
  int z2  // max for z counter
  int  L1  //  for selecting sex specific length data
  int  L2  //  used for l+nlength to get length bin for males
  int  A2  //  used for a+nages+1 to get true age bin for males
  int a1  // use to track a subset of ages
  int f // counter for fleets and surveys
  int fs  //  counter for fleets when looping across size and ageselex; so = f-Ntypes
  int g // counter for gmorph
  int gg  // counter for gender
  int a // counter for ages
  int b // counter for age bins
  int p // counter for area
  int p1
  int p2 // counter for destination area in migration
  int i // counter for observations
  int y // counter for year
  int yz // year, but not allowed to extend past endyr
  int s // counter for seasons
  int s2  // destination season
  int smid  // = s+nseas
  int t // counter for time, combining year and season
  int mo  //  month (1-12), not (0-11)
  int j
  int j1
  int j2
  int k
  int s_off  // offset for male section of vectors
  int Fishon  // whether or not to do fishery catch in equil_calc
  int NP  // number of parameters
  int Ip  // parameter counter
  int firstseas   // used to start season loops at the birthseason
  int t_base;    //
  int niter  // iteration count
  int loop
  int TG_t;  // time counter (in seasons) for tag groups
  int Fcast_catch_start
  int ParCount;
  int N_SC;  // counter for starter comments
  int N_DC;
  int N_CC;
  int N_FC;
  int depletion_fleet;  //  stores fleet(survey) number for the fleet that is defined as "depletion"

  int icycle
  int Ncycle
  int No_Report  //  flag to skip output reports after MCMC and MCeval
  number mcmcFlag
  number temp;
  number temp1;
  
  int Nparm_on_bound;
 !! No_Report=0;
 !! Ncycle=3;

//  ProgLabel_1.1.2   create elements of parameter labels
 LOCAL_CALCS
  adstring_array NumLbl;
  adstring_array GenderLbl;   // gender label
  adstring_array CRLF;   // blank to terminate lines

  CRLF+="";
  GenderLbl+="Fem";
  GenderLbl+="Mal";
  onenum="    ";
  for (i=1;i<=199;i++)
  {
  sprintf(onenum, "%d", i);
  NumLbl+=onenum+CRLF(1);
  }

  adstring sw;
  mcmcFlag = 0;
  for(i=0;i<argc;i++)
  {
    sw = argv[i];
    j=strcmp(sw,"-mcmc");
    if(j==0) {mcmcFlag = 1;}
    j=strcmp(sw,"-mceval");
    if(j==0) {mcmcFlag = 1;}
  }

//  ProgLabel_1.2    read the STARTER.SS file
  ad_comm::change_datafile_name("starter.ss");       //  get filenames
  cout<<" reading from STARTER.SS"<<endl;
  adstring checkchar;
  line_adstring readline;
  checkchar="";
  ifstream Starter_Stream("starter.ss");
   //  this opens a different logical file with a separate pointer from the pointer that ADMB uses when reading using init command to read from global_datafile
  k=0;
  N_SC=0;
  while(k==0)
  {
    Starter_Stream >>  readline;          // reads a single line from input stream
    if(length(readline)>2)
    {
      checkchar=readline(1);
      k=strcmp(checkchar,"#");
      checkchar=readline(1,2);
      j=strcmp(checkchar,"#C");
      if(j==0) {N_SC++; Starter_Comments+=readline;}
    }
  }
 END_CALCS
  !! echoinput<<version_info<<endl;
  !! echoinput<<ctime(&start)<<endl;
  !! warning<<version_info<<endl;
  !! warning<<ctime(&start)<<endl;
  init_adstring datfilename
  !!echoinput<<datfilename<<"  datfilename"<<endl;
  init_adstring ctlfilename
  !!echoinput<<ctlfilename<<"  ctlfilename"<<endl;
  init_int readparfile
  !!echoinput<<readparfile<<"  readparfile"<<endl;
  init_int rundetail
  !!echoinput<<rundetail<<"  rundetail"<<endl;
  init_int reportdetail
 LOCAL_CALCS
  echoinput<<reportdetail<<"  reportdetail"<<endl;
 END_CALCS

  init_int docheckup           // flag for ending dump to "checkup.SS"
   !!echoinput<<docheckup<<"  docheckup"<<endl;
  init_int Do_ParmTrace
   !!echoinput<<Do_ParmTrace<<"  Do_ParmTrace"<<endl;
  init_int Do_CumReport
   !!echoinput<<Do_CumReport<<"  Do_CumReport"<<endl;
  init_int Do_all_priors
   !!echoinput<<Do_all_priors<<"  Do_all_priors"<<endl;
  init_int SoftBound
   !!echoinput<<SoftBound<<"  SoftBound"<<endl;
  init_int N_nudata
   !!echoinput<<N_nudata<<"  N_nudata"<<endl;
  init_int Turn_off_phase
   !!echoinput<<Turn_off_phase<<"  Turn_off_phase"<<endl;

// read in burn and thinning intervals
  init_int burn_intvl
   !!echoinput<<burn_intvl<<"  MCeval burn_intvl"<<endl;
  init_int thin_intvl
   !!echoinput<<thin_intvl<<"  MCeval thin_intvl"<<endl;

  init_number jitter
   !!echoinput<<jitter<<"  jitter fraction for initial parm values"<<endl;

  init_int STD_Yr_min  // min yr for sdreport
   !!echoinput<<STD_Yr_min<<"  STD_Yr_min"<<endl;
  init_int STD_Yr_max  // max yr for sdreport
   !!echoinput<<STD_Yr_max<<"  STD_Yr_max"<<endl;
  init_int N_STD_Yr_RD  // N extra years to read
   !!echoinput<<N_STD_Yr_RD<<"  N extra STD years to read"<<endl;
  int N_STD_Yr
  init_ivector STD_Yr_RD(1,N_STD_Yr_RD)
   !!if(N_STD_Yr_RD>0) echoinput<<STD_Yr_RD<<"  vector of extra STD years"<<endl;
  // wait to process the above until after styr, endyr, N-forecast_yrs are read in data and forecast sections below

  int save_for_report;
  int save_gparm;
  int save_gparm_print;
  int N_warn;
  !! save_for_report=0;
  !! save_gparm=0;
  !! N_warn=0;

// set up the mcmc chain counter
  int mceval_counter;
  int mceval_header;
  !! mceval_counter = 0;
  !! mceval_header = 0;
  int mcmc_counter;
  !! mcmc_counter = 0;
  int done_run;
  !! done_run=0;

// set up the convergence criteria
  vector func_eval(1,50);
  vector func_conv(1,50);
//  number final_conv;
  init_number final_conv;
   !!echoinput<<final_conv<<"  final_conv"<<endl;

  !! func_eval.fill_seqadd(100,0);
  !! func_conv.fill_seqadd(1,0);
  !! func_conv(1)=10.;
  !! func_conv(2)=10.;

  init_int retro_yr;             //  introduce year for retrospective analysis
   !!echoinput<<retro_yr<<"  retro_yr"<<endl;
  int fishery_on_off;
  !! fishery_on_off=1;

  init_int Smry_Age
   !!echoinput<<Smry_Age<<"  Smry_Age"<<endl;

  init_int depletion_basis   // 0=skip; 1=fraction of B0; 2=fraction of Bmsy where fraction is depletion_level 3=rel to styr
   !!echoinput<<depletion_basis<<"  depletion_basis"<<endl;
  init_number depletion_level
   !!echoinput<<depletion_level<<"  depletion_level"<<endl;
  init_int SPR_reporting  // 0=skip; 1=SPR; 2=SPR_MSY; 3=SPR_Btarget; 4=(1-SPR)
   !!echoinput<<SPR_reporting<<"  SPR_reporting"<<endl;
  init_int F_reporting  // 0=skip; 1=exploit(Bio); 2=exploit(Num); 3=sum(frates); 4=true F
 LOCAL_CALCS
  echoinput<<F_reporting<<"  F_reporting"<<endl;
  if(F_reporting==4) {k=2;} else {k=0;}
 END_CALCS
  init_ivector F_reporting_ages_R(1,k);
  //  convert to F_reporting_ages later after nages is read.

 LOCAL_CALCS  
  if(k>0) 
    {
      echoinput<<F_reporting_ages<<"  F_reporting_ages_R"<<endl;
      echoinput<<"Will be checked against maxage later "<<endl;
    }
 END_CALCS

  init_int F_report_basis // 0=raw; 1=rel Fspr; 2=rel Fmsy ; 3=rel Fbtgt; 4=annual F for range of years
  !!echoinput<<F_report_basis<<"  F_report_basis"<<endl;
  !!echoinput<<"For Kobe plot, set depletion_basis=2; depletion_level=1.0; F_reporting=your choose; F_report_basis=2"<<endl;
  init_int finish_starter

 LOCAL_CALCS
   if(finish_starter!=999) {cout<<"CRITICAL error reading starter file "<<finish_starter<<endl; exit(1);}
   echoinput<<"  finish reading starter.ss"<<endl<<endl;
 END_CALCS

  //  end reading  from Starter file

  number pi
  !!  pi=3.14159265358979;

  number neglog19
  !!  neglog19 = -log(19.);

  number NilNumbers           //  used as the minimum for posfun and similar checks
//  !! NilNumbers = 0.0000001;
  !! NilNumbers = 0.000;

//  ProgLabel_1.2.1  Set up a dummy datum for use when max phase = 0
  number dummy_datum
  int dummy_phase
  !! dummy_datum=1.;
  !! if(Turn_off_phase<=0) {dummy_phase=0;} else {dummy_phase=-6;}

  int runnumber
  int N_prof_var;
  int prof_var_cnt
  int prof_junk

 LOCAL_CALCS
//  ProgLabel_1.3   Read runnumber.ss
   ifstream fin1("runnumber.ss", ios::in);
    if (fin1)
    {
      fin1>>runnumber;
      runnumber++;
      fin1.close();
    }
    else
    {
      runnumber=1;
    }
//  ProgLabel_1.3.1   increment runnumber and write to file
    ofstream fin2("runnumber.ss", ios::out);
    fin2 << runnumber;
    fin2.close();

//  ProgLabel_1.4   Read Profilevalues.ss file
   N_prof_var=998;
   ifstream fin3("profilevalues.ss", ios::in);
   fin3>>N_prof_var;   //  if file is null this will not return anything
   if (N_prof_var==998) {N_prof_var=0; prof_junk=0;} else {prof_junk=1;}
   fin3.close();
   if(N_prof_var>0)
   {ad_comm::change_datafile_name("profilevalues.ss");}
   else
   {ad_comm::change_datafile_name("runnumber.ss");}  // just to have something in scope
   prof_var_cnt=(runnumber-1)*N_prof_var+2;
 END_CALCS
   init_vector prof_var(1,prof_junk+runnumber*N_prof_var);


//****************************************************************
//  ProgLabel_2.0  READ DATA FILE

//  ProgLabel_2.1  Read comments and dimension info
//  ProgLabel_2.1.1  Read and save #C comments at top of data file
 LOCAL_CALCS
  ad_comm::change_datafile_name(datfilename);
  cout<<" reading from data file"<<endl;
  ifstream Data_Stream(datfilename);   //  even if the global_datafile name is used, there still is a different logical device created
  k=0;
  N_DC=0;
  while(k==0)
  {
    Data_Stream >>  readline;          // reads the line from input stream
    if(length(readline)>2)
    {
      checkchar=readline(1);
      k=strcmp(checkchar,"#");
      checkchar=readline(1,2);
      j=strcmp(checkchar,"#C");
      if(j==0) {N_DC++; Data_Comments+=readline;}
    }
  }
 END_CALCS

//  ProgLabel_2.1.2  read model time dimensions
  init_int styr  //start year of the model
 !!echoinput<<styr<<" start year "<<endl;

  init_int endyr // end year of the model
 !!echoinput<<endyr<<" end year "<<endl;

  init_int nseas  //  number of seasons
 !!echoinput<<nseas<<" N seasons "<<endl;

  init_vector seasdur(1,nseas) // season duration; enter in units of months, fractions OK; will be rescaled to sum to 1.0 if total is greater than 11.9
 !!echoinput<<seasdur<<" months/seas (fractions OK) "<<endl;

  int TimeMax
  int TimeMax_Fcast_std
  int YrMax;

  int eq_yr;
  int bio_yr;
  number sumseas;

//  ProgLabel_2.1.3  Set up seasons
  vector seasdur_half(1,nseas)   // half a season
  vector azero_seas(1,nseas);   // cumulative time, within year, up until begin of this season
//  ivector Mo_seas(1,12) //  season that a month is within

 LOCAL_CALCS
  sumseas=sum(seasdur);
  if(sumseas>=11.9)
    {seasdur /=sumseas;}
  else
    {seasdur /=12.;}
  seasdur_half = seasdur*0.5;   // half a season
  TimeMax = styr+(endyr-styr)*nseas+nseas-1;
  retro_yr=endyr+retro_yr;

  azero_seas(1)=0.;
  if(nseas>1)
    for(s=2;s<=nseas;s++)
    {azero_seas(s)=sum(seasdur(1,s-1));}

//  calculate the season that the start of a month occurs within
//  s=1;
//  Mo_seas=1;  //  default is all months are in seas 1
//  if(nseas>1)
//  {
//    for (mo=2;mo<=12;mo++)
//    {
//      if( float(mo-1)/12.>=(azero_seas(s)-0.000001) ) {s++;}
//      Mo_seas(mo)=s;
//    }
//  }
//  echoinput<<"Mo_seas: "<<Mo_seas<<" season that a month occurs in "<<endl;

 END_CALCS
 !!echoinput<<seasdur<<" processed season duration "<<endl;

   init_int spawn_seas    //  spawning occurs at beginning of this season
   number spawn_time_ann  //  later this will be entered as a fraction of the year, 
   number spawn_time_seas  //  spawn_seas will be calculated from spawn_time_ann, then spawn_time_seas will be timing within spawn_seas
 LOCAL_CALCS
   spawn_time_ann=0.0;
   spawn_time_seas=0.0;
   echoinput<<spawn_seas<<" spawn_seas "<<endl;
   if(spawn_seas>nseas)
   {
     N_warn++; warning<<" spawn_seas index must be <= nseas "<<endl;
   }
 END_CALCS
//  ProgLabel_2.1.4  Define fleets, surveys and areas
  init_int Nfleet // number of fleets (each one with a different selectivity function)
 !!echoinput<<Nfleet<<" Nfleet "<<endl;

  init_int Nsurvey
  !!echoinput<<Nsurvey<<" Nsurvey "<<endl;

 int      Ntypes
  !! Ntypes=Nfleet+Nsurvey;
  init_int pop   // number of areas (populations)
  !!echoinput<<pop<<" N_areas "<<endl;

  imatrix pfleetname(1,Ntypes,1,2)
  init_adstring fleetnameread;

 LOCAL_CALCS
  if(pop>1 && F_reporting==3)
  {N_warn++; warning<<" F-reporting=3 (sum of full Fs) not advised in multiple area models "<<endl;}
  for(f=1;f<=Ntypes;f++) {pfleetname(f,1)=1; pfleetname(f,2)=1;}    // set whole array to equal 1 in case not enough names are read
  f=1;
  for(i=1;i<=strlen(fleetnameread);i++)
  if(adstring(fleetnameread(i))==adstring("%"))
   {pfleetname(f,2)=i-1; f+=1;  pfleetname(f,1)=i+1;}
  pfleetname(Ntypes,2)=strlen(fleetnameread);
  for(f=1;f<=Ntypes;f++)
  {
    fleetname+=fleetnameread(pfleetname(f,1),pfleetname(f,2))+CRLF(1);
  }
  echoinput<<fleetname<<endl;
 END_CALCS

  init_vector surveytime(1,Ntypes)   // fraction of season (not year) in which survey occurs
  !!echoinput<<surveytime<<" surveytime(fishery+surveys) "<<endl;

  init_ivector fleet_area(1,Ntypes)    // areas in which each fleet/survey operates
  !!echoinput<<fleet_area<<" fleet_area(fishery+surveys) "<<endl;

  init_vector catchunits(1,Nfleet)  // 1=biomass; 2=numbers
  !!echoinput<<catchunits<<" catchunits "<<endl;

  init_vector catch_se_rd(1,Nfleet)  // units are se of log(catch); use -1 to ignore input catch values for discard only fleets
  !!echoinput<<catch_se_rd<<" catch_se "<<endl;

  matrix catch_se(styr-1,TimeMax,1,Nfleet);

//  ProgLabel_2.1.5  define genders and max age
  init_int    gender  //  number of sexes in the model
  !!echoinput<<gender<<" N genders "<<endl;

  init_int    nages // Number of ages
  !!echoinput<<nages<<" nages is maxage "<<endl;

  ivector     age_vector(0,nages)
  vector      r_ages(0,nages)
  vector frac_ages(0,nages)
  ivector     years(styr,endyr) // vector of the years of the model
  vector    r_years(styr,endyr);
  ivector F_reporting_ages(1,2);

 LOCAL_CALCS
  for (a=0;a<=nages;a++) age_vector(a) = a; //age vector
  for (a=0;a<=nages;a++) r_ages(a) = double(a); //age vector
  frac_ages=r_ages/r_ages(nages);
  for (y=styr;y<=endyr;y++) {years(y)=y; r_years(y)=y;}    //year vector
  if (F_reporting==4)
  {
    F_reporting_ages=F_reporting_ages_R;
    if(F_reporting_ages(1)>(nages-2) || F_reporting_ages(1)<0)
    {N_warn++; warning<<" reset lower end of F_reporting_ages to be nages-2  "<<endl; F_reporting_ages(1)=nages-2;}
    if(F_reporting_ages(2)>(nages-2) || F_reporting_ages(2)<0)
    {N_warn++; warning<<" reset upper end of F_reporting_ages to be nages-2  "<<endl; F_reporting_ages(2)=nages-2;}
  }
  else
  {
    F_reporting_ages(1)=nages/2;
    F_reporting_ages(2)=F_reporting_ages(1)+1;
  }
 END_CALCS

//  ProgLabel_2.1.6  create arrays to hold pointer for time/fleet data occurrence
  matrix yr_cr2(1,Ntypes,styr,TimeMax)    // reverse pointers to data occurrence
  matrix yr_disc2(1,Ntypes,styr,TimeMax)
  3darray yr_mnwt2(1,Ntypes,styr,TimeMax,0,2)  // last dimension here is for total, discard, retain
                                               //  store negative value if the observation is mean length, rather than mean weight
  matrix have_data(1,Ntypes,styr,TimeMax)
  !! have_data.initialize();

//  ProgLabel_2.2  Read CATCH amount by fleet
  init_vector obs_equ_catch(1,Nfleet)    //  initial, equilibrium catch.  Annual.
  !!echoinput<<obs_equ_catch<<" obs_equ_catch "<<endl;

  init_int N_ReadCatch;
  !!echoinput<<N_ReadCatch<<" N_ReadCatch "<<endl;

  init_matrix catch_bioT(1,N_ReadCatch,1,Nfleet+2)
  !!echoinput<<" catch as read "<<endl<<catch_bioT<<endl;

  matrix catch_ret_obs(1,Nfleet,styr-nseas,TimeMax+nseas)
  !! k=Nfleet*(TimeMax-styr+1);
  ivector have_catch(1,k);
  3darray catch_seas_area(styr,TimeMax,1,pop,0,Nfleet)
  matrix totcatch_byarea(styr,TimeMax,1,pop)
  vector totcat(styr-1,endyr)
  int first_catch_yr

 LOCAL_CALCS
    catch_ret_obs.initialize();
    t=styr-nseas;
    for(f=1;f<=Nfleet;f++)
      {catch_ret_obs(f,t) = obs_equ_catch(f);}    //  put all init equil catch in first period
  for(k=1;k<=N_ReadCatch;k++)
  {
    y=catch_bioT(k,Nfleet+1); s=catch_bioT(k,Nfleet+2);
    if(y>=styr && y<=endyr)
    {
      if(s>nseas) s=nseas;   // allows for collapsing multiple season catch data down into fewer seasons
                             //  typically to collapse to annual because accumulation will all be in the index "nseas"
      if(s>0)
      {
        t=styr+(y-styr)*nseas+s-1;
        for(f=1;f<=Nfleet;f++) catch_ret_obs(f,t) += catch_bioT(k,f);
      }
      else  // distribute catch equally across seasons
      {
        for(s=1;s<=nseas;s++)
        {
          t=styr+(y-styr)*nseas+s-1;
          for(f=1;f<=Nfleet;f++) catch_ret_obs(f,t) += catch_bioT(k,f)/nseas;
        }
      }
    }
  }
  echoinput<<" processed catch "<<endl<<trans(catch_ret_obs)<<endl;

  for(t=styr-1;t<=TimeMax;t++)
  for(f=1;f<=Nfleet;f++)
   {
    if(catch_se_rd(f)>0.)
      {catch_se(t,f)=catch_se_rd(f);}
      else
      {catch_se(t,f)=0.1;}  //  for bycatch fleets
    }
        

//  calc total catch by year so can calculate the first year with catch and to omit zero catch years from sdreport
  totcat.initialize();
  catch_seas_area.initialize();
  totcatch_byarea.initialize();
  totcat(styr-1)=sum(obs_equ_catch);
  first_catch_yr=0;
  if(totcat(styr-1)>0.0) first_catch_yr=styr-1;
  for(y=styr; y<=endyr; y++)
  {
    for(s=1;s<=nseas;s++)
    {
      t=styr+(y-styr)*nseas+s-1;
      for(p=1;p<=pop;p++)
      for(f=1;f<=Nfleet;f++)
      if(fleet_area(f)==p && catch_ret_obs(f,t) > 0.0)
      {
        totcat(y) += catch_ret_obs(f,t);
        catch_seas_area(t,p,f)=1;
        catch_seas_area(t,p,0)=1;
        totcatch_byarea(t,p)+=catch_ret_obs(f,t);
      }
    }
    if(totcat(y)>0.0 && first_catch_yr==0) first_catch_yr=y;
    if(y==endyr && totcat(y)==0.0)
    {
      N_warn++; warning<<" catch is 0.0 in endyr; this will cause failure in the benchmark and forecast calculations"<<endl;
    }
  }
 END_CALCS


//  ProgLabel_2.3 Read fishery CPUE, effort, and Survey index or abundance
  init_int nobs_cr_rd
  int nobs_cr
 LOCAL_CALCS
  echoinput<<nobs_cr_rd<<" nobs_survey "<<endl;
  if(nobs_cr_rd>0) {k=3; j=Ntypes;} else {k=0; j=Ntypes;}
 END_CALCS
  init_imatrix cr_units_rd(1,Ntypes,1,k)
  ivector cr_units(1,j)   //0=num/1=bio/2=F
  ivector cr_errtype(1,j)  // -1=normal / 0=lognormal / >0=T

 LOCAL_CALCS
  if(k>0)
  {
    echoinput<<"Units:  0=numbers; 1=biomass; 2=F"<<endl;
    echoinput<<"Errtype:  -1=normal; 0=lognormal; >0=T"<<endl;
    echoinput<<"Fleet Units Err_Type"<<endl;
    echoinput<<cr_units_rd<<endl;
    cr_units=column(cr_units_rd,2);
    cr_errtype=column(cr_units_rd,3);
  }
  else
    {
      cr_units=0;
      cr_errtype=0;
    }
 END_CALCS

   init_matrix indexdata(1,nobs_cr_rd,1,5)
  !!if(nobs_cr_rd>0) echoinput<<" indexdata "<<endl<<indexdata<<endl;
  ivector nyr_cr(1,Ntypes)
  int in_superperiod
  ivector N_suprper_cr(1,Ntypes)      // N super_yrs per obs

 LOCAL_CALCS
  //  count the number of observations, exclude those outside the specified year range, count the number of superperiods
  nobs_cr=0;
  nyr_cr=0;
  N_suprper_cr=0;
  if(nobs_cr_rd>0)
  {
    for(i=1;i<=nobs_cr_rd;i++)
    {
      y=indexdata(i,1);
      s=abs(indexdata(i,2));  // because neg season indicates super period
      if(s>nseas)
      {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for survey obs "<<i<<" is > nseas"<<endl; exit(1);}
      if(y>=styr && y<=retro_yr)
      {
        f=abs(indexdata(i,3));  //  negative f turns off observation
        nyr_cr(f)++;
        if(indexdata(i,5)<0) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"cannot use negative se to indicate superperiods in survey data"<<endl; exit(1);}
        if(indexdata(i,2)<0) N_suprper_cr(f)++;  // count the super-periods if seas<0
      }
    }
    nobs_cr=sum(nyr_cr);
    for(f=1;f<=Ntypes;f++)
    if(N_suprper_cr(f)>0)
    {
      j=N_suprper_cr(f)/2;  // because we counted the begin and end
      if(2*j!=N_suprper_cr(f))
      {
        N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"unequal number of starts and ends of survey superperiods "<<endl; exit(1);
      }
      else
      {
        N_suprper_cr(f)=j;
      }
    }
  }
 END_CALCS

  imatrix yr_cr(1,Ntypes,1,nyr_cr)
  imatrix yr_cr_y(1,Ntypes,1,nyr_cr)
  imatrix yr_cr_s(1,Ntypes,1,nyr_cr)
  imatrix yr_cr_use(1,Ntypes,1,nyr_cr)
  matrix  obs_cr(1,Ntypes,1,nyr_cr)
  matrix  Ln_obs_cr(1,Ntypes,1,nyr_cr)
  matrix  se_cr_obs(1,Ntypes,1,nyr_cr)
  matrix vul_bio(1,Ntypes,1,nyr_cr);        // Vulnerable biomass
// arrays for Super-years
  imatrix suprper_cr1(1,Ntypes,1,N_suprper_cr)  //  where N_suprper_cr is a vector
  imatrix suprper_cr2(1,Ntypes,1,N_suprper_cr)
  matrix suprper_cr_sampwt(1,Ntypes,1,nyr_cr)
  
 LOCAL_CALCS
//  ProgLabel_2.3.1  process survey observations, move info into working arrays,create super-periods as needed
    N_suprper_cr=0;
    nyr_cr=0;
    in_superperiod=0;
    if(nobs_cr>0)  // proceed if any data in yr range
    {
      for(i=1;i<=nobs_cr_rd;i++)  // loop all, including those out of yr range
      {
        y=indexdata(i,1);
        if(y>=styr && y<=retro_yr)
        {
          s=abs(indexdata(i,2));  // because neg season indicates super period
          t=styr+(y-styr)*nseas+s-1;
          f=abs(indexdata(i,3));
          nyr_cr(f)++;
          j=nyr_cr(f);
          yr_cr(f,j)=t;
          yr_cr_y(f,j)=y;
          se_cr_obs(f,j)=indexdata(i,5);   // later adjust with varadjust, copy to se_cr_use, then adjust with extra se parameter
          yr_cr_s(f,j)=indexdata(i,2);  // so keeps sign of this value
          if(indexdata(i,3)<0) {yr_cr_use(f,j)=-1;} else {yr_cr_use(f,j)=1;}

          obs_cr(f,j)=indexdata(i,4);
          //  create super_year indexes
          if(indexdata(i,2)<0) // start or stop a super-period;  ALL observations must be continguous in the file
          if(in_superperiod==0) // start superperiod
           {N_suprper_cr(f)++; suprper_cr1(f,N_suprper_cr(f))=j; in_superperiod=1;}
          else
          {
            if(in_superperiod==1)  // end superperiod
            {
              suprper_cr2(f,N_suprper_cr(f))=j;
              in_superperiod=0;
            }
            else
            {
            }
          }
        }
      }
      echoinput<<" processed survey data "<<endl;
      for(f=1;f<=Ntypes;f++)
      {echoinput<<f<<" "<<fleetname(f)<<" "<<obs_cr(f)<<endl;}
    }
    echoinput<<"Successful read of index data; N= "<<nobs_cr<< endl;
    echoinput<<"Number of survey superperiods by fleet: "<<N_suprper_cr<<endl;
 END_CALCS

//  ProgLabel_2.4  read Discard data
   init_int Ndisc_fleets
   int nobs_disc

 LOCAL_CALCS
  echoinput<<" note order of discard read is now: N fleets with disc, then if Ndisc_fleets>0 read:  fleet, disc_units, disc_error(for 1,Ndisc_fleets), then read obs "<<endl;
  echoinput<<Ndisc_fleets<<" N fleets with discard "<<endl;

  if(Ndisc_fleets>0) {j=Ntypes;} else {j=0;}
 END_CALCS
  init_imatrix disc_units_rd(1,Ndisc_fleets,1,3)
  ivector disc_units(1,j)  //  formerly scalar disc_type
  ivector disc_errtype(1,j)  // formerly scalar DF_disc
  vector disc_errtype_r(1,j)  // real version for T-dist
 LOCAL_CALCS
  disc_units.initialize();
  disc_errtype.initialize();
  if(Ndisc_fleets>0)
  {
    echoinput<<"#_discard_units (1=same_as_catchunits(bio/num);2=fraction; 3=numbers)"<<endl;
    echoinput<<"#_discard_error:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal"<<endl;
    echoinput<<"#Fleet Units Err_Type"<<endl;
    echoinput<<disc_units_rd<<endl;
    for(j=1;j<=Ndisc_fleets;j++)
    {
      f=disc_units_rd(j,1);
      disc_units(f)=disc_units_rd(j,2);
      disc_errtype(f)=disc_units_rd(j,3);
      disc_errtype_r(f)=float(disc_errtype(f));
    }
  }
 END_CALCS
   init_int nobs_disc_rd
  !! echoinput<<nobs_disc_rd<<" N discard obs "<<endl;

  init_matrix discdata(1,nobs_disc_rd,1,5)
  !! if(nobs_disc_rd>0) echoinput<<" discarddata "<<endl<<discdata<<endl;
   ivector nyr_disc(1,Ntypes)
   ivector N_suprper_disc(1,Ntypes)      // N super_yrs per obs

 LOCAL_CALCS
  nobs_disc=0;
  nyr_disc=0;
  N_suprper_disc=0;
  if(nobs_disc_rd>0)
  {
    for(i=1;i<=nobs_disc_rd;i++)  // get count of observations in date range
    {
      y=discdata(i,1);
      s=abs(discdata(i,2));
      if(s>nseas)
        {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for discard obs "<<i<<" is > nseas"<<endl; exit(1);}
      if(y>=styr && y<=retro_yr)
      {
        f=abs(discdata(i,3));
        nyr_disc(f)++;
        if(discdata(i,5)<0) {N_warn++; cout<<"EXIT - see warning"<<endl; warning<<"Cannot use negative se as indicator of superperiod in discard data"<<endl;}
        if(discdata(i,2)<0) N_suprper_disc(f)++;  // count the super-periods if seas<0 or se<0
      }
    }
    nobs_disc=sum(nyr_disc);  // sum of obs in the date range
    for(f=1;f<=Ntypes;f++)
    if(N_suprper_disc(f)>0)
    {
      j=N_suprper_disc(f)/2;  // because we counted the begin and end
      if(2*j!=N_suprper_disc(f))
      {
        N_warn++; cout<<"EXIT - see warning"<<endl; warning<<"unequal number of starts and ends of discard superperiods "<<endl; exit(1);
      }
      else
      {
        N_suprper_disc(f)=j;
      }
    }
  }
 END_CALCS

  imatrix yr_disc(1,Ntypes,1,nyr_disc)
  imatrix yr_disc_y(1,Ntypes,1,nyr_disc)
  imatrix yr_disc_s(1,Ntypes,1,nyr_disc)
  imatrix yr_disc_use(1,Ntypes,1,nyr_disc)
  matrix  obs_disc(1,Ntypes,1,nyr_disc)
  matrix  cv_disc(1,Ntypes,1,nyr_disc)
  matrix  sd_disc(1,Ntypes,1,nyr_disc)
// arrays for Super-years
  imatrix suprper_disc1(1,Ntypes,1,N_suprper_disc)
  imatrix suprper_disc2(1,Ntypes,1,N_suprper_disc)
  matrix suprper_disc_sampwt(1,Ntypes,1,nyr_disc)
 LOCAL_CALCS
//  ProgLabel_2.4.1  process discard data and create super periods as needed
    nyr_disc=0;                        // redo the counter to provide pointer for below
    N_suprper_disc=0;
    in_superperiod=0;
    if(nobs_disc>0)
    {
      for(i=1;i<=nobs_disc_rd;i++)
      {
        y=discdata(i,1);
        if(y>=styr && y<=retro_yr)
        {
         s=abs(discdata(i,2));
         f=abs(discdata(i,3));
         t=styr+(y-styr)*nseas+s-1;
         nyr_disc(f)++;
         j=nyr_disc(f);
         yr_disc(f,j)=t;
         yr_disc_y(f,j)=y;
         yr_disc_s(f,j)=discdata(i,2);
         cv_disc(f,j)=discdata(i,5);
         obs_disc(f,j)=discdata(i,4);
         if(discdata(i,4)<0.0) discdata(i,3)=-abs(discdata(i,3));  //  convert to new format using negative fleet
         if(discdata(i,3)<0) {yr_disc_use(f,j)=-1;} else {yr_disc_use(f,j)=1;}
         if(catch_ret_obs(f,t)<=0.0)
         {
           N_warn++; warning<<" discard observation but no corresponding catch, yr, seas, fleet "<<y<<" "<<s<<" "<<f<<endl;
         }

  //  create super_year indexes
         if(discdata(i,2)<0)  // start/stop a super-year  ALL observations must be continguous in the file
         {
         if(in_superperiod==0)  // start a super-year
           {N_suprper_disc(f)++; suprper_disc1(f,N_suprper_disc(f))=j; in_superperiod=1;}
         else if(in_superperiod==1)  // end a super-year
           {suprper_disc2(f,N_suprper_disc(f))=j; in_superperiod=0;}
         }
        }
      }
    }
    echoinput<<"Successful read of discard data, N by fleet = "<< nobs_disc << endl;
    echoinput<<"N superperiods by fleet = "<< N_suprper_disc << endl;
 END_CALCS


//  ProgLabel_2.5  read Mean Body Weight data
   init_int nobs_mnwt_rd
   int nobs_mnwt
  !!echoinput<<nobs_mnwt_rd<<" nobs_mnwt "<<endl;
  init_number DF_bodywt  // DF For meanbodywt T-distribution
  !!echoinput<<DF_bodywt<<" degrees of freedom for bodywt T-distribution "<<endl;

   init_matrix mnwtdata1(1,nobs_mnwt_rd,1,6)
 LOCAL_CALCS
  nobs_mnwt=0;
  if(nobs_mnwt_rd>0) echoinput<<" meanbodywt_data "<<endl<<mnwtdata1<<endl;
  for(i=1;i<=nobs_mnwt_rd;i++)
  {
    y=mnwtdata1(i,1);
    if(y>=styr && y<=retro_yr) nobs_mnwt++;
  }
 END_CALCS
  matrix mnwtdata(1,9,1,nobs_mnwt)
//    yr, seas, type, mkt, obs, CV

 LOCAL_CALCS
  yr_mnwt2=0;
  if(nobs_mnwt>0)
  j=0;
  for(i=1;i<=nobs_mnwt_rd;i++)
  {
    y=mnwtdata1(i,1);
    if(y>=styr && y<=retro_yr)
    {
      j++;
      s=mnwtdata1(i,2);
      if(s<1) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for meanwt obs "<<i<<" is <0; superper is not implemented for meanwt"<<endl; exit(1);}
      if(s>nseas) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for meanwt obs "<<i<<" is > nseas"<<endl; exit(1);}
      t=styr+(y-styr)*nseas+s-1;
      f=abs(mnwtdata1(i,3));
      z=mnwtdata1(i,4);  // z is partition (0, 1, 2);  use (10, 11, 12) as input to indicate use of mean length,, rather than mean weight
      if(z<3)
      {yr_mnwt2(f,t,z)=j;}  //  datum is mean weight
      else
      {yr_mnwt2(f,t,z-10)=-j;}  //  datum is mean length
      have_data(f,t)=1;
      for(k=1;k<=6;k++) mnwtdata(k,j)=mnwtdata1(i,k);
    }
  }
  echoinput<<"Successful read of mean-body-wt data, N= "<< nobs_mnwt <<endl;
 END_CALCS


// **************************************
//  ProgLabel_2.6  setup population Length bins
  number binwidth2  //  width of length bins in population
  number minLread  // input minimum size in population; this is used as the mean size at age 0.00
  number maxLread  //  input maximum size to be considered; should be divisible by binwidth2
  int nlen_bin2  //number of length bins in length comp data doubled for males
  int nlen_binP   //number of length bins in length comp data +1 as needed
  number minL               // minL and maxL store ends of the sizevector and are used as bounds later
  number minL_m  // mean size in first pop bin
  number maxL  // set to the midsize of last population bin for selex calc
  int nlength  // N pop lenbins
  int nlength1  //  +1 as needed
  int nlength2  // doubled for males
  number startbin  // population length bin that matches first data length bin

  init_int LenBin_option  // 1=set to data bins; 2 generate uniform; 3 = read custom
  !!echoinput<<LenBin_option<<" LenBin_option "<<endl;
 LOCAL_CALCS
   if(LenBin_option==1)
   {k=0;}
   else if(LenBin_option==2)
   {k=3;}
   else if(LenBin_option==3)
   {k=1;}
   else
   {N_warn++; warning<<" LenBin_option must be 1, 2 or 3"<<LenBin_option<<endl;}
 END_CALCS

  init_vector PopBin_Read(1,k);
  !!if( k>0) echoinput<<PopBin_Read<<" input for setup of pop length bins "<<endl;
 LOCAL_CALCS
   nlength=0;  // later will be read or calculated
   if(LenBin_option==2)
   {binwidth2=PopBin_Read(1); minLread=PopBin_Read(2); maxLread=PopBin_Read(3);}
   else if(LenBin_option==3)
   {nlength=PopBin_Read(1);}  // number of bins to read
 END_CALCS
  init_vector len_bins_rd(1,nlength)
  !!if(nlength>0) echoinput<<len_bins_rd<<" pop length bins as read "<<endl;

//  ProgLabel_2.7  start length data section
//  ProgLabel_2.7.1  read and process data length bins
  init_number min_tail  //min_proportion_for_compressing_tails_of_observed_composition
  !!echoinput<<min_tail<<" min tail for comps "<<endl;
  init_number min_comp  //  small value added to each composition bins
  !!echoinput<<min_comp<<" value added to comps "<<endl;
  init_int CombGender_l  //  combine genders through this length bin
  !!echoinput<<CombGender_l<<" CombGender_lengths "<<endl;
  init_int nlen_bin //number of length bins in length comp data
 LOCAL_CALCS
  echoinput<<nlen_bin<<" nlen_bin_for_data "<<endl;
  if(CombGender_l>nlen_bin)
  {
    N_warn++; warning<<"Combgender_l cannot be greater than nlen_bin; resetting"<<endl;  CombGender_l=nlen_bin;
  }
  nlen_binP=nlen_bin+1;
  nlen_bin2=gender*nlen_bin;
 END_CALCS

  init_vector len_bins_dat(1,nlen_bin) // length bin lower boundaries
  !!echoinput<<" len_bins_dat "<<endl<<len_bins_dat<<endl;
  vector len_bins_dat2(1,nlen_bin2)  //; doubled for males; for output only
  vector len_bins_dat_m(1,nlen_bin)  //; midbin; for output only

//  ProgLabel_2.7.2  process length bins, create mean length per bin, etc.
 LOCAL_CALCS
   if(LenBin_option==1)
   {
     nlength = nlen_bin;  // set N pop bins same as data bins
   }
   else if(LenBin_option==2)
   {
     nlength = (maxLread-minLread)/binwidth2+1;   // number of population length bins
   }
   else if(LenBin_option==3)
   {
    // nlength was read
   }
   nlength1 = nlength+1;          //  +1 when needed
   nlength2 = gender*nlength;    // doubled for males
 END_CALCS

  vector len_bins(1,nlength)  //vector with lower edge of population length bins
  vector len_bins2(1,nlength2)  //vector with lower edge of population length bins
  vector binwidth(1,nlength2)  //ve
  vector len_bins_m(1,nlength)  //vector with mean size in bin
  vector len_bins_m2(1,nlength2)  //vector with all length bins; doubled for males
  vector len_bins_sq(1,nlength2)  //vector with all length bins; doubled for males
  vector male_offset(1,nlength2)  // used to calculate retained@length as population quantity

 LOCAL_CALCS
  male_offset.initialize();  //  initialize
   if(LenBin_option==1)
   {len_bins=len_bins_dat;}
   else if(LenBin_option==2)
   {
     len_bins(1)=minLread;
     for (z=2;z<=nlength;z++)  {len_bins(z)=len_bins(z-1)+binwidth2;}
   }
   else
   {len_bins=len_bins_rd;}

  if(len_bins(1)==0.0) len_bins(1)=0.001;
  for (z=1;z<=nlength;z++)
  {
    len_bins2(z)=len_bins(z);
    if(z<nlength)
    {
      len_bins_m(z) = (len_bins(z+1)+len_bins(z))/2.;
      binwidth(z)=len_bins(z+1)-len_bins(z);
    }
    else
    {
      len_bins_m(z) = len_bins(z)+binwidth(z-1)/2.;
      binwidth(z)=binwidth(z-1);
    }

    len_bins_m2(z) = len_bins_m(z);     // for use in calc mean size at binned age
    len_bins_sq(z) = len_bins_m2(z)*len_bins_m2(z);        //  for use in calc std dev of size at binned age
    if(gender==2)
    {
      len_bins2(z+nlength)=len_bins(z);
      male_offset(z+nlength)=1.;
      binwidth(z+nlength)=binwidth(z);
      len_bins_m2(z+nlength)=len_bins_m2(z);
      len_bins_sq(z+nlength)=len_bins_sq(z);
    }
  }
  echoinput<<endl<<"Processed Population length bin info "<<endl;

  maxL=len_bins_m(nlength);
  minL=len_bins(1);
  minL_m=len_bins_m(1);
  if(LenBin_option!=2) binwidth2=binwidth(nlength/2);  // set a reasonable value in case LenBin_option !=2

  startbin=1;
  while(len_bins(startbin)<len_bins_dat(1))
  {startbin++;}

  for (z=1;z<=nlen_bin;z++)
  {
    len_bins_dat2(z) = len_bins_dat(z);
    if(gender==2) len_bins_dat2(z+nlen_bin)=len_bins_dat(z);
    if(z<nlen_bin)
    {
      len_bins_dat_m(z)=0.5*(len_bins_dat(z)+len_bins_dat(z+1));  //  this is not gender specific
    }
    else
    {
      len_bins_dat_m(z)=len_bins_dat_m(z-1)+ (len_bins_dat(z)-len_bins_dat(z-1));
    }
  }
  if(len_bins_dat(nlen_bin)>len_bins(nlength))
  {
    N_warn++; cout<<"Critical error, see warning.sso"<<endl; warning<<" Data length bins extend beyond pop len bins "<<len_bins_dat(nlen_bin)<<" "<<len_bins(nlength)<<endl; exit(1);
  }
  echoinput<<endl<<"Processed Data length bin info "<<endl;
 END_CALCS

//  ProgLabel_2.7.3  create conversion of pop length bins to data length bins
  matrix make_len_bin(1,nlen_bin2,1,nlength2);

  int ibin;
  int ibinsave;
  int fini;
  number topbin;
  number botbin;

 LOCAL_CALCS
  make_len_bin.initialize();
  ibin=0;
  topbin=0.;
  botbin=0.;
  for(z=1;z<=nlength;z++)
  {
    if(ibin==nlen_bin)
    {  //echoinput<<" got to last ibin, so put rest of popbins here"<<endl;
      make_len_bin(ibin,z)=1.;
    }
    else
    {
      if(len_bins(z)>=topbin) {ibin++; }  //echoinput<<" incr ibin ";

      if(ibin>1)  {botbin=len_bins_dat(ibin);}
      if(ibin<nlen_bin) {topbin=len_bins_dat(ibin+1);} else {topbin=99999.;}

      if(ibin==nlen_bin)  // echoinput<<" got to last ibin, so put rest of popbins here"<<endl;
      {make_len_bin(ibin,z)=1.;}
      else if(len_bins(z)>=botbin && len_bins(z+1)<=topbin )  //echoinput<<" pop inside dat, put here"<<endl;
      {make_len_bin(ibin,z)=1.;}
      else
      {
      make_len_bin(ibin+1,z)=(len_bins(z+1)-topbin)/(len_bins(z+1)-len_bins(z));
      if(ibin!=1) make_len_bin(ibin,z)=1.-make_len_bin(ibin+1,z);
      }
    }
  }
  if(gender==2)
  {
    for (i=1;i<=nlen_bin;i++)
    for(j=1;j<=nlength;j++)
    make_len_bin(i+nlen_bin,j+nlength)=make_len_bin(i,j);
  }
   if(docheckup==1)
   {
    echoinput<<"make_len_bin "<<len_bins<<endl;
    for(ibin=1;ibin<=nlen_bin;ibin++) echoinput<<len_bins_dat(ibin)<<make_len_bin(ibin)(1,nlength)<<endl;
   }
  echoinput<<endl<<"Processed Population to Data length bin conversion "<<endl;
 END_CALCS


//  ProgLabel_2.7.4 read Length composition data
   init_int nobsl_rd
   int Nobs_l_tot
   !!echoinput<<nobsl_rd<<" N length comp obs "<<endl;
   init_matrix lendata(1,nobsl_rd,1,6+nlen_bin2)
   !!if(nobsl_rd>0) echoinput<<" first lencomp obs "<<endl<<lendata(1)<<endl<<" last obs"<<endl<<lendata(nobsl_rd)<<endl;;

  ivector Nobs_l(1,Ntypes)
  ivector N_suprper_l(1,Ntypes)      // N super_yrs per obs

 LOCAL_CALCS
    Nobs_l=0;                       //  number of observations from each fleet/survey
    N_suprper_l=0;
    if(nobsl_rd>0)
    for(i=1;i<=nobsl_rd;i++)
    {
      y=lendata(i,1);
      if(y>=styr && y<=retro_yr)
      {
      f=abs(lendata(i,3));
      s=abs(lendata(i,2));
      if(lendata(i,6)<0) {N_warn++; cout<<"error in length data "<<endl; warning<<"Error: negative sample size no longer valid as indicator of skip data or superperiods "<<endl; exit(1);}
      if(lendata(i,2)<0) N_suprper_l(f)++;     // count the number of starts and ends of super-periods if seas<0
      Nobs_l(f)++;
      }
    }
  for(f=1;f<=Ntypes;f++) 
  {
    s=N_suprper_l(f)/2.;
    if(s*2==N_suprper_l(f))
    {N_suprper_l(f)=s;}// to get the number of superperiods
    else
    {N_warn++; cout<<"error in length data "<<endl; warning<<"Error: unequal number of length superperiod starts and stops "<<endl; exit(1);}
  }

    echoinput<<"Lendata Nobs by fleet "<<Nobs_l<<endl;
    echoinput<<"Lendata superperiods by fleet "<<N_suprper_l<<endl;
    Nobs_l_tot=sum(Nobs_l);
 END_CALCS

  imatrix yr_l_t(1,Ntypes,1,Nobs_l)
  imatrix yr_l_y(1,Ntypes,1,Nobs_l)
  3darray obs_l(1,Ntypes,1,Nobs_l,1,nlen_bin2)
  3darray obs_l_all(1,4,1,Ntypes,1,nlen_bin)  //  for the sum of all length comp data
  matrix  nsamp_l(1,Ntypes,1,Nobs_l)
  matrix  nsamp_l_read(1,Ntypes,1,Nobs_l)
  imatrix  gen_l(1,Ntypes,1,Nobs_l)
  imatrix  mkt_l(1,Ntypes,1,Nobs_l)
  3darray header_l(1,Ntypes,1,Nobs_l,0,3)
  3darray tails_l(1,Ntypes,1,Nobs_l,1,4)   // min-max bin for females; min-max bin for males
  ivector tails_w(1,4)

// arrays for Super-years
  imatrix suprper_l1(1,Ntypes,1,N_suprper_l)
  imatrix suprper_l2(1,Ntypes,1,N_suprper_l)
  matrix  suprper_l_sampwt(1,Ntypes,1,Nobs_l)  //  will contain calculated weights for obs within super periods
  int floop
  int tloop

//  ProgLabel_2.7.5  process length comps, compress tails, add constant, scale to 1.0
 LOCAL_CALCS
  N_suprper_l=0;
  Nobs_l=0;
  in_superperiod=0;
  suprper_l1.initialize();
  suprper_l2.initialize();
  obs_l_all.initialize();

  if(Nobs_l_tot>0)
  {
    for(floop=1;floop<=Ntypes;floop++)    //  loop fleets
    for(i=1;i<=nobsl_rd;i++)   //  loop all observations to find those for this fleet/time
    {
      y=lendata(i,1);
      if(y>=styr && y<=retro_yr)
      {
        s=abs(lendata(i,2));
        if(s>nseas)
         {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for length obs "<<i<<" is > nseas"<<endl; exit(1);}
        t=styr+(y-styr)*nseas+s-1;
        f=abs(lendata(i,3));
        if(f==floop)
        {
          Nobs_l(f)++;
          j=Nobs_l(f);
//          if(lendata(i,6)==-9999 || lendata(i,6)==-9998) lendata(i,2)=-abs(lendata(i,2));  // convert to new super-period syntax
          if(lendata(i,6)<0.0)
            {N_warn++; warning<<"negative values not allowed for lengthcomp sample size, use -fleet to omit from -logL"<<endl;}
          header_l(f,j,1) = y;
          header_l(f,j,2) = lendata(i,2);  // season
          header_l(f,j,3) = lendata(i,3);   // fleet
          header_l(f,j,0) = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
          yr_l_y(f,j)=y;                   // year
          yr_l_t(f,j)=t;                     // sequential time = year+season
          gen_l(f,j)=lendata(i,4);         // gender 0=combined, 1=female, 2=male, 3=both
          mkt_l(f,j)=lendata(i,5);         // partition: 0=all, 1=discard, 2=retained
          nsamp_l_read(f,j)=lendata(i,6);  // assigned sample size for observation
          nsamp_l(f,j)=nsamp_l_read(f,j);

//  ProgLabel_2.7.6  create super-periods for length compositions
          if(lendata(i,2)<0)  // start/stop a super-period  ALL observations must be continguous in the file
          {
            if(in_superperiod==0)  // start a super-period  ALL observations must be continguous in the file
            {N_suprper_l(f)++; suprper_l1(f,N_suprper_l(f))=j; in_superperiod=1;}
            else if(in_superperiod==1)  // end a super-year
            {suprper_l2(f,N_suprper_l(f))=j; in_superperiod=0;}
          }
          
          for(z=1;z<=nlen_bin2;z++)   // get the composition vector
           {obs_l(f,j,z)=lendata(i,6+z);}

          if(sum(obs_l(f,j))<=0.0) {N_warn++; warning<<" zero fish in size comp (fleet, year) "<<f<<" "<<y<<endl; cout<<" EXIT - see warning "<<endl; exit(1);}
          tails_l(f,j,1)=1;
          tails_l(f,j,2)=nlen_bin;
          tails_l(f,j,3)=nlen_binP;
          tails_l(f,j,4)=nlen_bin2;
          if(gen_l(f,j)==3 && gender==2 && CombGender_l>0)
          {
            for(z=1;z<=CombGender_l;z++)
            {
              obs_l(f,j,z)+=obs_l(f,j,z+nlen_bin);  obs_l(f,j,z+nlen_bin)=0.0;
            }
            tails_l(f,j,3)=nlen_binP+CombGender_l;
          }
          if(gen_l(f,j)==2) obs_l(f,j)(1,nlen_bin) = 0.;   //  zero out females for male-only obs
          if(gen_l(f,j)<=1 && gender==2) obs_l(f,j)(nlen_binP,nlen_bin2) = 0.;   //  zero out males for female-only or combined gender obs
          obs_l(f,j) /= sum(obs_l(f,j));                  // make sum to 1.00

          if(gen_l(f,j)!=2)                      // do females, unless Male-only observation
          {
            k=0;
            temp=sum(obs_l(f,j)(1,nlen_bin));
            for(z=1;z<=nlen_bin;z++)
            if(obs_l(f,j,z)>0.)   //  find Nbins with data
            {k++;}
            if(temp>0.0 && k>1)     // only compress tail if obs exist for this gender and there is more than 1 bin with data
            {
              k=0;
              for(z=1;z<=nlen_bin-1;z++)             // compress Female lower tail until exceeds min_tail
              {
                if(obs_l(f,j,z)<=min_tail && k==0)
                {
                  obs_l(f,j,z+1)+=obs_l(f,j,z);
                  obs_l(f,j,z)=0.00;
                  tails_l(f,j,1)=z+1;
                }
                else
                {k=1;}
              }

              k=0;
              for(z=nlen_bin;z>=tails_l(f,j,1);z--)  // compress Female upper tail until exceeds min_tail
              {
                if(obs_l(f,j,z)<=min_tail && k==0)
                {
                  obs_l(f,j,z-1)+=obs_l(f,j,z);
                  obs_l(f,j,z)=0.00;
                  tails_l(f,j,2)=z-1;
                }
                else
                {k=1;}
              }
            }
            obs_l(f,j)(tails_l(f,j,1),tails_l(f,j,2)) += min_comp;    // add min_comp to bins in range
          }

          if(gen_l(f,j)>=2 && gender==2) // process males
          {
            k=0;
            temp=sum(obs_l(f,j)(nlen_binP,nlen_bin2));
            for(z=nlen_binP;z<=nlen_bin2;z++)
            if(obs_l(f,j,z)>0.)
            {k++;}
            if(temp>0.0 && k>1)     // only compress tail if obs exist for this gender and there is more than 1 bin with data
            {
              k=0;
              k1=tails_l(f,j,3);
              for(z=k1;z<=nlen_bin2-1;z++)
              {
                if(obs_l(f,j,z)<=min_tail && k==0)
                {
                  obs_l(f,j,z+1)+=obs_l(f,j,z);
                  obs_l(f,j,z)=0.00;
                  tails_l(f,j,3)=z+1;
                }
                else
                {k=1;}
              }

              k=0;
              for(z=nlen_bin2;z>=tails_l(f,j,3);z--)  // compress Male upper tail until exceeds min_tail
              {
                if(obs_l(f,j,z)<=min_tail && k==0)
                {
                  obs_l(f,j,z-1)+=obs_l(f,j,z);
                  obs_l(f,j,z)=0.00;
                  tails_l(f,j,4)=z-1;
                }
                else
                {k=1;}
              }
            }
            obs_l(f,j)(tails_l(f,j,3),tails_l(f,j,4)) += min_comp;  // add min_comp to bins in range
          }   // end doing males
          obs_l(f,j) /= sum(obs_l(f,j));                  // make sum to 1.00 again after adding min_comp
          if(gender==1 || gen_l(f,j)!=2) {obs_l_all(1,f)(1,nlen_bin)+=obs_l(f,j)(1,nlen_bin);}  //  females or combined
          if(gender==2)
          {
            if(gen_l(f,j)==1 || gen_l(f,j)==3)  // put females into female only
            {
              obs_l_all(3,f)(1,nlen_bin)+=obs_l(f,j)(1,nlen_bin);
            }
            if(gen_l(f,j)>=2)  //  put males into combined and into male only
            {
              for(z=1;z<=nlen_bin;z++)
              {
                obs_l_all(1,f,z)+=obs_l(f,j,nlen_bin+z);
                obs_l_all(4,f,z)+=obs_l(f,j,nlen_bin+z);
              }
            }
          }
        }
      }
    }
  }

     echoinput<<"Overall_Compositions"<<endl<<"Fleet len_bins "<<len_bins_dat<<endl;
     for(f=1;f<=Ntypes;f++)
     {
       for(j=1;j<=4;j++)
       {
          if(j!=2)
          {
            temp=sum(obs_l_all(j,f));
            if(temp>0.0)
            {
               obs_l_all(j,f)/=temp;
            }
            else
            {
              obs_l_all(j,f)=float(1./nlen_bin);
            }
          }
       }
       
       obs_l_all(2,f,1)=obs_l_all(1,f,1); // first bin
       for(z=2;z<=nlen_bin;z++)
       {
         obs_l_all(2,f,z)=obs_l_all(2,f,z-1)+obs_l_all(1,f,z);
       }
       if(Nobs_l(f)>0)
       {
         echoinput<<f<<" freq"<<obs_l_all(1,f)<<endl;
         echoinput<<f<<" cuml"<<obs_l_all(2,f)<<endl;
         echoinput<<f<<" female"<<obs_l_all(3,f)<<endl;
         echoinput<<f<<" male"<<obs_l_all(4,f)<<endl;
       }
     }

  echoinput<<"Successful processing of length data"<<endl<<endl;
 END_CALCS


//  ProgLabel_2.8  start age composition data section
//  ProgLabel_2.8.1  Read Age bin and ageing error vectors
  init_int n_abins // age classes for data
   !!echoinput<<n_abins<<" N age bins "<<endl;
  int n_abins1;
  int n_abins2;
  int Use_AgeKeyZero;  //  set to ageerr_type for the age data that use parameter approach
  int AgeKeyParm;  //  holds starting parm number for age error parameters
  int AgeKey_StartAge;
  int AgeKey_Linear1;
  int AgeKey_Linear2;

 LOCAL_CALCS
  n_abins1=n_abins+1;
  n_abins2=gender*n_abins;
 END_CALCS

   init_vector age_bins1(1,n_abins) // age classes for data
  !!echoinput << age_bins1<< " agebins "  << endl;
  init_int N_ageerr   // number of ageing error matrices to be calculated
   !!echoinput<<N_ageerr<<" N age error defs "<<endl;
  init_3darray age_err_rd(1,N_ageerr,1,2,0,nages) // ageing imprecision as stddev for each age
 LOCAL_CALCS
  Use_AgeKeyZero=0;
  if(N_ageerr>0)
  {
    for(i=1;i<=N_ageerr;i++)
    {
      if(age_err_rd(i,2,0)<0.) Use_AgeKeyZero=i;  //  set flag for setup of age error parameters
    }
    echoinput<<"set agekey: "<<Use_AgeKeyZero<<" to create key from parameters"<<endl;
  }
 END_CALCS
  vector age_bins(1,n_abins2) // age classes for data
  vector age_bins_mean(1,n_abins2)  //  holds mean age for each data age bin
//  number age;

//  ProgLabel_2.8.2  Read Age data
  init_int nobsa_rd
  int Nobs_a_tot
  !!echoinput << nobsa_rd<< " N ageobs"  << endl;

  init_int Lbin_method  //#_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths
  !!echoinput << Lbin_method<< " Lbin method for defined size ranges "  << endl;

  init_int CombGender_a  //  combine genders through this age bin

 LOCAL_CALCS
  echoinput << CombGender_a<< " merge males&females through this bin number "  << endl;
  if(nobsa_rd>0 && N_ageerr==0)
  {
    N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" must define ageerror vectors because age data exist"<<endl; exit(1);
  }
  if(CombGender_a>n_abins2)
  {
    N_warn++; warning<<"Combgender_a cannot be greater than n_abins; resetting"<<endl;  CombGender_a=n_abins2;
  }
  nlen_binP=nlen_bin+1;
  nlen_bin2=gender*nlen_bin;
 END_CALCS
  init_matrix agedata(1,nobsa_rd,1,9+n_abins2)
   !!if(nobsa_rd>0) echoinput<<" first agecomp obs "<<endl<<agedata(1)<<endl<<" last obs"<<endl<<agedata(nobsa_rd)<<endl;;
  ivector Nobs_a(1,Ntypes)
  ivector N_suprper_a(1,Ntypes)      // N super_yrs per obs

 LOCAL_CALCS
  for (b=1;b<=n_abins;b++)
  {
   age_bins(b) = age_bins1(b);
   if(b<n_abins)
    {
      age_bins_mean(b) =(age_bins1(b)+age_bins1(b+1))*0.5;
    }
    else
    {
      age_bins_mean(b) =age_bins1(b)+0.5*(age_bins1(b)-age_bins1(b-1));
    }
    if(gender==2)
      {
        age_bins(b+n_abins)=age_bins1(b);
        age_bins_mean(b+n_abins)=age_bins_mean(b);
      }
  }
  Nobs_a=0;
  N_suprper_a=0;

  for(i=1;i<=nobsa_rd;i++)
  {
    y=agedata(i,1);
    if(y>=styr && y<=retro_yr)
    {
     f=abs(agedata(i,3));
     if(agedata(i,9)<0) {N_warn++; cout<<"error in age data "<<endl; warning<<"Error: negative sample size no longer valid as indicator of skip data or superperiods "<<endl; exit(1);}
     if(agedata(i,6)==0 || agedata(i,6)>N_ageerr) {N_warn++; cout<<"error in age data "<<endl; warning<<"Error: undefined age_error type: "<<agedata(i,6)<<"  in obs: "<<i<<endl; exit(1);}
     if(agedata(i,2)<0) N_suprper_a(f)++;     // count the number of starts and ends of super-periods if seas<0 or sampsize<0
    
     Nobs_a(f)++;
    }
  }
  for(f=1;f<=Ntypes;f++) 
  {
    s=N_suprper_a(f)/2.;
    if(s*2==N_suprper_a(f))
    {N_suprper_a(f)/=2;}
    else
    {N_warn++; cout<<"error in age data "<<endl; warning<<"Error: unequal number of age superperiod starts and stops "<<endl; exit(1);}
  }
  echoinput<<endl<<"Agedata Nobs by fleet "<<Nobs_a<<endl;
  echoinput<<"Agedata superperiods by fleet "<<N_suprper_a<<endl;
  Nobs_a_tot=sum(Nobs_a);

  for(f=1;f<=Ntypes;f++)
      {if(Nobs_a(f)==0) Nobs_a(f)=1;}  //  why is this needed?
 END_CALCS

  imatrix yr_a_t(1,Ntypes,1,Nobs_a)
  imatrix yr_a_y(1,Ntypes,1,Nobs_a)
  3darray obs_a(1,Ntypes,1,Nobs_a,1,gender*n_abins)
  3darray obs_a_all(1,2,1,Ntypes,1,n_abins)  //  for the sum of all age comp data
  matrix  nsamp_a(1,Ntypes,1,Nobs_a)
  matrix  nsamp_a_read(1,Ntypes,1,Nobs_a)
  imatrix  ageerr_type_a(1,Ntypes,1,Nobs_a)
  imatrix  gen_a(1,Ntypes,1,Nobs_a)
  imatrix  mkt_a(1,Ntypes,1,Nobs_a)
  3darray  Lbin_filter(1,Ntypes,1,Nobs_a,1,nlength2)
  imatrix  use_Lbin_filter(1,Ntypes,1,Nobs_a)
  imatrix  Lbin_lo(1,Ntypes,1,Nobs_a)
  imatrix  Lbin_hi(1,Ntypes,1,Nobs_a)
  3darray tails_a(1,Ntypes,1,Nobs_a,1,4)   // min-max bin for females; min-max bin for males
  3darray header_a(1,Ntypes,1,Nobs_a,0,9)

// arrays for Super-years
  matrix  suprper_a_sampwt(1,Ntypes,1,Nobs_a)  //  will contain calculated weights for obs within super periods
  imatrix suprper_a1(1,Ntypes,1,N_suprper_a)
  imatrix suprper_a2(1,Ntypes,1,N_suprper_a)

//  ProgLabel_2.8.3  pre-process age comps, compress tails, define length bin filters
 LOCAL_CALCS
  Lbin_filter=1.;
  use_Lbin_filter.initialize();     // have to use initialize; imatrix cannot be set to a constant
  suprper_a1.initialize();
  suprper_a2.initialize();
  obs_a_all.initialize();
  N_suprper_a=0;
  Nobs_a=0;
  if(Nobs_a_tot>0)
   {
     echoinput<<"process age comps "<<endl;
     in_superperiod=0;
     for(floop=1;floop<=Ntypes;floop++)
     for(i=1;i<=nobsa_rd;i++)
     {
       y=agedata(i,1);
       if(y>=styr && y<=retro_yr)
       {
         s=abs(agedata(i,2));
         if(s>nseas)
          {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for agecomp obs "<<i<<" is > nseas"<<endl; exit(1);}
         t=styr+(y-styr)*nseas+s-1;
         f=abs(agedata(i,3));
         if(f==floop)
         {
           Nobs_a(f)++;  //  redoing this pointer just to create index j used below
           j=Nobs_a(f);
           header_a(f,j)(1,9)=agedata(i)(1,9);
           yr_a_t(f,j)=t;                       // sequential time = year+season
           yr_a_y(f,j)=y;
           s=abs(agedata(i,2));
//           header_a(f,j,0) = float(y)+float(s-1.)/nseas;
           header_a(f,j,0) = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));

           gen_a(f,j)=agedata(i,4);         // gender 0=combined, 1=female, 2=male, 3=both
           mkt_a(f,j)=agedata(i,5);         // partition: 0=all, 1=discard, 2=retained
           if(agedata(i,6)>N_ageerr)
           {
              N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" ageerror type must be <= "<<N_ageerr<<endl; exit(1);
           }
           ageerr_type_a(f,j)=agedata(i,6);
           nsamp_a_read(f,j)=agedata(i,9);
           nsamp_a(f,j)=nsamp_a_read(f,j);
//  ProgLabel_2.8.4  create super-periods for age compositions
           if(in_superperiod==0 && agedata(i,2)<0)  // start a super-year  ALL observations must be continguous in the file
           {N_suprper_a(f)++; suprper_a1(f,N_suprper_a(f))=j; in_superperiod=1;}
           else if(in_superperiod==1 && agedata(i,2)<0)  // end a super-year
           {suprper_a2(f,N_suprper_a(f))=j; in_superperiod=0;}

           for(b=1;b<=gender*n_abins;b++)   // get the composition vector
           {obs_a(f,j,b)=agedata(i,9+b);}
           if(sum(obs_a(f,j))<=0.0)
           {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" zero fish in age comp "<<header_a(f,j)<<endl;  exit(1);}

           Lbin_lo(f,j)=agedata(i,7);
           Lbin_hi(f,j)=agedata(i,8);
           switch (Lbin_method)   //  here all 3 methods are converted to poplenbins for use internally
           {
             case 1:  // values are population length bin numbers
             {
               if(Lbin_lo(f,j)<=0) Lbin_lo(f,j)=1;
               if(Lbin_hi(f,j)<=0 || Lbin_hi(f,j)>nlength) Lbin_hi(f,j)=nlength;
               break;
             }
             case 2:  // values are data length bin numbers
             {
               if(Lbin_lo(f,j)<=0) Lbin_lo(f,j)=1;
               if(Lbin_hi(f,j)<=0 || Lbin_hi(f,j)>nlen_bin) Lbin_hi(f,j)=nlen_bin;
               s=0;
               for(k=1;k<=nlength;k++)
               {
                 if( len_bins(k)==len_bins_dat(Lbin_lo(f,j)) ) s=k;  //  find poplen bin that matches data len bin
               }
               if(s==0) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" L_bin_lo no match to poplenbins in age comp "<<header_a(f,j)<<endl;  exit(1);}
               Lbin_lo(f,j)=s;

               s=0;
               for(k=1;k<=nlength;k++)
               {
                 if( len_bins(k)==len_bins_dat(Lbin_hi(f,j)) ) s=k;   //  find poplen bin that matches data len bin
               }
               if(s==0) {N_warn++; cout<<" exit - see warning "<<endl; warning<<" L_bin_hi no match to poplenbins in age comp "<<header_a(f,j)<<endl;  exit(1);}
               Lbin_hi(f,j)=s;
               break;
             }
             case 3:   // values are lengths
             {
               if(Lbin_lo(f,j)<=0) Lbin_lo(f,j)=len_bins(1);
               if(Lbin_hi(f,j)<=0 || Lbin_hi(f,j)>len_bins(nlength)) Lbin_hi(f,j)=len_bins(nlength);
               s=0;
               for(k=1;k<=nlength;k++)
               {
                 if( len_bins(k)==Lbin_lo(f,j) ) s=k;  //  find poplen bin that matches input length for lbin_lo
               }
               if(s==0) {N_warn++; cout<<" exit - see warning "<<endl; warning<<"L_bin_lo no match to poplenbins in age comp "<<header_a(f,j)<<endl;  exit(1);}
               Lbin_lo(f,j)=s;

               s=0;
               for(k=1;k<=nlength;k++)
               {
                 if( len_bins(k)==Lbin_hi(f,j) ) s=k;
               }
               if(s==0) {N_warn++; cout<<" exit - see warning "<<endl; warning<<"L_bin_hi no match to poplenbins in age comp "<<header_a(f,j)<<endl;  exit(1);}
               Lbin_hi(f,j)=s;
               break;
             }
           }

           //  lbin_lo and lbin_hi are now in terms of poplenbins; their original values are retained in header_a
           if(Lbin_lo(f,j)>nlength || Lbin_lo(f,j)>Lbin_hi(f,j))
           {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" L_bin_lo is too high in age comp.  Are you using lengths or bin numbers? "<<header_a(f,j)<<endl;  exit(1);}
           if(Lbin_lo(f,j)==1 && Lbin_hi(f,j)==nlength) {use_Lbin_filter(f,j)=0;} else {use_Lbin_filter(f,j)=1;}

           if(use_Lbin_filter(f,j)==1)
           {                                                   // use Lbin_filter for this obs
             Lbin_filter(f,j)=0.;
             Lbin_filter(f,j)(Lbin_lo(f,j),Lbin_hi(f,j)) = 1;
             if(gender==2)
             {
              k=int(Lbin_lo(f,j))+nlength; z=int(Lbin_hi(f,j))+nlength;
              Lbin_filter(f,j)(k,z) = 1.;
             }
           }

           if(gen_a(f,j)==2) obs_a(f,j)(1,n_abins) = 0.;   //  zero out females for male-only obs
           if(gen_a(f,j)<=1&&gender==2) obs_a(f,j)(n_abins1,gender*n_abins) = 0.;   //  zero out males for female-only or combined gender obs
           tails_a(f,j,1)=1;
           tails_a(f,j,2)=n_abins;
           tails_a(f,j,3)=1+(gender-1)*n_abins;
           tails_a(f,j,4)=gender*n_abins;
           if(gender==2)
           {
             if(gen_a(f,j)==3 && CombGender_a>0)
             {
              for(z=1;z<=CombGender_a;z++)
              {
                obs_a(f,j,z)+=obs_a(f,j,z+n_abins);  obs_a(f,j,z+n_abins)=0.0;
              }
              tails_a(f,j,3)=n_abins+1+CombGender_a;
             }
           }

           obs_a(f,j) /= sum(obs_a(f,j));

           if(gen_a(f,j)!=2)                     // do females, unless Male-only observation
           {
             k=0;
             temp=sum(obs_a(f,j)(1,n_abins));
             for(z=1;z<=n_abins;z++)
             if(obs_a(f,j,z)>0.)
             {k++;}
             if(temp>0.0 && k>1)     // only compress tail if obs exist for this gender and there is more than 1 bin with data
             {
               k=0;
               for(z=1;z<=n_abins-1;z++)             // compress Female lower tail until exceeds min_tail
               {
                 if(obs_a(f,j,z)<=min_tail && k==0)
                 {
                   obs_a(f,j,z+1)+=obs_a(f,j,z);
                   obs_a(f,j,z)=0.00;
                   tails_a(f,j,1)=z+1;
                 }
                 else
                 {k=1;}
               }

               k=0;
               for(z=n_abins;z>=tails_a(f,j,1);z--)  // compress Female upper tail until exceeds min_tail
               {
                 if(obs_a(f,j,z)<=min_tail && k==0)
                 {
                   obs_a(f,j,z-1)+=obs_a(f,j,z);
                   obs_a(f,j,z)=0.00;
                   tails_a(f,j,2)=z-1;
                 }
                 else
                 {k=1;}
               }
             }
             obs_a(f,j)(tails_a(f,j,1),tails_a(f,j,2)) += min_comp;    // add min_comp to bins in range
           }                            // done with females

           if(gen_a(f,j)>=2 && gender==2) // compress Male tails until exceeds min_tail
           {
             k=0;
             temp=sum(obs_a(f,j)(n_abins1,n_abins2));
             for(z=n_abins1;z<=n_abins2;z++)
             if(obs_a(f,j,z)>0.)
             {k++;}
             if(temp>0.0 && k>1)     // only compress tail if obs exist for this gender and there is more than 1 bin with data
             {
               k=0;

              for(z=n_abins1;z<=n_abins2-1;z++)
              {
                if(obs_a(f,j,z)<=min_tail && k==0)
                {
                  obs_a(f,j,z+1)+=obs_a(f,j,z);
                  obs_a(f,j,z)=0.00;
                  tails_a(f,j,3)=z+1;
                }
                else
                {k=1;}
              }

              k=0;
              for(z=n_abins2;z>=tails_a(f,j,3);z--)  // compress Male upper tail until exceeds min_tail
              {
                if(obs_a(f,j,z)<=min_tail && k==0)
                {
                  obs_a(f,j,z-1)+=obs_a(f,j,z);
                  obs_a(f,j,z)=0.00;
                  tails_a(f,j,4)=z-1;
                }
                else
                {k=1;}
              }
            }
            obs_a(f,j)(tails_a(f,j,3),tails_a(f,j,4)) += min_comp;  // add min_comp to bins in range
           }
           if(sum(obs_a(f,j))>0.) obs_a(f,j) /= sum(obs_a(f,j));                  // make sum to 1.00 again after adding min_comp
           if(gender==1 || gen_a(f,j)!=2) obs_a_all(1,f)(1,n_abins)+=obs_a(f,j)(1,n_abins);  //  females or combined
           if(gender==2 && gen_a(f,j)>=2)
           {
             for(a=1;a<=n_abins;a++)
             {
               obs_a_all(1,f,a)+=obs_a(f,j,n_abins+a);  //  males
             }
           }
         }
       }
     }

     echoinput<<"Fleet age_bins "<<age_bins<<endl;
     for(f=1;f<=Ntypes;f++)
     {
       if(Nobs_a(f)>0)
       {
         obs_a_all(1,f)/=sum(obs_a_all(1,f));
       }
       else
       {
         obs_a_all(1,f)=float(1./n_abins);
       }
       obs_a_all(2,f,1)=obs_a_all(1,f,1); // first bin
       for(a=2;a<=n_abins;a++)
       {
         obs_a_all(2,f,a)=obs_a_all(2,f,a-1)+obs_a_all(1,f,a);
       }
       echoinput<<f<<" freq "<<obs_a_all(1,f)<<endl;
       echoinput<<f<<" cuml "<<obs_a_all(2,f)<<endl;
     }
   }
   echoinput<<endl<<"Successful processing of age data "<<endl;
 END_CALCS

// ***********************************************************************
//  ProgLabel_2.9  read mean Size_at_Age data

  init_int nobs_ms_rd
  int nobs_ms_tot
  !!echoinput<<nobs_ms_rd<<" N mean size-at-age obs "<<endl;
  init_matrix sizeagedata(1,nobs_ms_rd,1,7+2*n_abins2)
   !!if(nobs_ms_rd>0) echoinput<<" first size-at-age obs "<<endl<<sizeagedata(1)<<endl<<" last obs"<<endl<<sizeagedata(nobs_ms_rd)<<endl;;
  ivector Nobs_ms(1,Ntypes)
  ivector N_suprper_ms(1,Ntypes)      // N super_yrs per obs

 LOCAL_CALCS
  Nobs_ms=0;
  N_suprper_ms=0;
  if(nobs_ms_rd>0)
  for(i=1;i<=nobs_ms_rd;i++)
  {
    y=sizeagedata(i,1);
    if(y>=styr && y<=retro_yr)
    {
      f=abs(sizeagedata(i,3));
      if(sizeagedata(i,7)<0) {N_warn++; cout<<"error in meansize"<<endl; warning<<"error.  cannot use negative sampsize for meansize data ";exit(1);;}
      if(sizeagedata(i,2)<0) N_suprper_ms(f)++;     // count the number of starts and ends of super-periods if seas<0 or sampsize<0
      Nobs_ms(f)++;
    }
  }
  for(f=1;f<=Ntypes;f++) 
  {
    s=N_suprper_ms(f)/2.;
    if(s*2==N_suprper_ms(f))
    {N_suprper_ms(f)/=2;}
    else
    {N_warn++; cout<<"error in meansize data "<<endl; warning<<"Error: unequal number of meansize superperiod starts and stops "<<endl; exit(1);}
  }
  echoinput<<endl<<"meansize data Nobs by fleet "<<Nobs_ms<<endl;
  echoinput<<"meansize superperiods by fleet "<<N_suprper_ms<<endl;

  nobs_ms_tot=sum(Nobs_ms);
 END_CALCS

  imatrix yr_ms_t(1,Ntypes,1,Nobs_ms)
  imatrix yr_ms_y(1,Ntypes,1,Nobs_ms)
  3darray obs_ms(1,Ntypes,1,Nobs_ms,1,n_abins2)
  3darray obs_ms_n(1,Ntypes,1,Nobs_ms,1,n_abins2)
  3darray obs_ms_n_read(1,Ntypes,1,Nobs_ms,1,n_abins2)
  imatrix  ageerr_type_ms(1,Ntypes,1,Nobs_ms)
  imatrix  gen_ms(1,Ntypes,1,Nobs_ms)
  imatrix  mkt_ms(1,Ntypes,1,Nobs_ms)
  imatrix  use_ms(1,Ntypes,1,Nobs_ms)
  3darray header_ms(1,Ntypes,1,Nobs_ms,0,7)
  matrix suprper_ms_sampwt(1,Ntypes,1,Nobs_ms)
  imatrix suprper_ms1(1,Ntypes,1,N_suprper_ms)
  imatrix suprper_ms2(1,Ntypes,1,N_suprper_ms)

 LOCAL_CALCS
   Nobs_ms=0;
   suprper_ms1.initialize();
   suprper_ms2.initialize();
   N_suprper_ms.initialize();
   if(nobs_ms_tot>0)
   {
     in_superperiod=0;
     for(floop=1;floop<=Ntypes;floop++)
     for(i=1;i<=nobs_ms_rd;i++)
     {
       y=sizeagedata(i,1);
       if(y>=styr && y<=retro_yr)
       {
         s=abs(sizeagedata(i,2));
         if(s>nseas)
         {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for mean size-at-age obs "<<i<<" is > nseas"<<endl; exit(1);}
         t=styr+(y-styr)*nseas+s-1;
         f=abs(sizeagedata(i,3));
         if(f==floop)
         {
           Nobs_ms(f)++;
           j=Nobs_ms(f);
           header_ms(f,j)(1,7)=sizeagedata(i)(1,7);
           yr_ms_t(f,j)=t;
           yr_ms_y(f,j)=y;
           s=abs(sizeagedata(i,2));
//           header_ms(f,j,0) = float(y)+float(s-1.)/nseas;
           header_ms(f,j,0) = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
           gen_ms(f,j)=sizeagedata(i,4);
           mkt_ms(f,j)=sizeagedata(i,5);
           if(sizeagedata(i,6)>N_ageerr)
           {
              N_warn++;cout<<" EXIT - see warning "<<endl;
              warning<<" in meansize, ageerror type must be <= "<<N_ageerr<<endl; exit(1);
           }
           ageerr_type_ms(f,j)=sizeagedata(i,6);
           if(sizeagedata(i,3)>0) {use_ms(f,j)=1;} else {use_ms(f,j)=-1;}

//  ProgLabel_2.9.1  create super-periods for meansize data
           if(sizeagedata(i,2)<0)  // start/stop a super-period  ALL observations must be continguous in the file
           {
             if(in_superperiod==0)  //  start superperiod
             {N_suprper_ms(f)++; suprper_ms1(f,N_suprper_ms(f))=j; in_superperiod=1;}
             else
             if(in_superperiod==1)  // end a super-period
             {suprper_ms2(f,N_suprper_ms(f))=j; in_superperiod=0;}
           }

           for(b=1;b<=n_abins2;b++)
           {obs_ms(f,j,b)=sizeagedata(i,7+b);}
           for(b=1;b<=n_abins2;b++)
           {
             obs_ms_n(f,j,b)=sizeagedata(i,7+b+n_abins2);
             obs_ms_n_read(f,j,b)=sizeagedata(i,7+b+n_abins2);
           }
         }
       }
     }
   }

   echoinput<<"Successful read of size-at-age data; N = "<<Nobs_ms<<endl;
 END_CALCS


//  ProgLabel_2.10  read environmental data that will be used to modify processes and expected values
    init_int N_envvar
  !!echoinput<<N_envvar<<" N_envvar "<<endl;
    init_int N_envdata
  !!echoinput<<N_envdata<<" N_envdata "<<endl;

    matrix env_data_RD(styr-1,endyr+100,1,N_envvar)  //  leave enough room for N_Fcast_Yrs which is not yet known
    init_matrix env_temp(1,N_envdata,1,3)
 LOCAL_CALCS
     if(N_envdata>0)
     {
      if(N_envdata>0)echoinput<<" env data "<<endl<<env_temp<<endl;
      env_data_RD=0.;
      for (i=1;i<=N_envdata;i++)
        if(env_temp(i,1)>=(styr-1))
        {env_data_RD(env_temp(i,1), env_temp(i,2) ) = env_temp(i,3);}
     }
 END_CALCS

//  ProgLabel_2.11  start generalized size composition section
//  ProgLabel_2.11.1  read generalized size frequency data (aka wt frequency)
  int SzFreqMethod;
  int iobs;
  init_int SzFreq_Nmeth;                                   // number of sizefreq methods to be read
  !!echoinput<<SzFreq_Nmeth<<" N sizefreq methods to read "<<endl;
  4darray SzFreq_HaveObs(1,Ntypes,1,SzFreq_Nmeth,styr,TimeMax,0,2)
  imatrix SzFreq_HaveObs2(1,SzFreq_Nmeth,styr,TimeMax)
  init_ivector SzFreq_Nbins(1,SzFreq_Nmeth);               //  number of bins for each method
  !!if(SzFreq_Nmeth>0) echoinput<<SzFreq_Nbins<<" Sizefreq N bins per method"<<endl;
  init_ivector SzFreq_units(1,SzFreq_Nmeth);               //  units for proportions (1 = biomass; 2=numbers ) for each method
  !!if(SzFreq_Nmeth>0) echoinput<<SzFreq_units<<" Sizetfreq units(bio/num) per method"<<endl;
  init_ivector SzFreq_scale(1,SzFreq_Nmeth);               //  bin scale (1=kg; 2=lbs; 3=cm; 4=in) for each method
  !!if(SzFreq_Nmeth>0) echoinput<<SzFreq_scale<<" Sizefreq scale(kg/lbs/cm/inches) per method"<<endl;
  init_vector SzFreq_mincomp(1,SzFreq_Nmeth);               //  mincomp to add for each method
  !!if(SzFreq_Nmeth>0) echoinput<<SzFreq_mincomp<<" Sizefreq mincomp per method "<<endl;
  init_ivector SzFreq_nobs(1,SzFreq_Nmeth);
  !!if(SzFreq_Nmeth>0) echoinput<<SzFreq_nobs<<" Sizefreq N obs per method"<<endl;

  ivector SzFreq_Nbins_seas_g(1,SzFreq_Nmeth*nseas);   //  array dimensioner used only for the SzFreqTrans array
  ivector SzFreq_Nbins3(1,SzFreq_Nmeth)        // doubles the Nbins if gender==2
  int SzFreqMethod_seas;

 LOCAL_CALCS
  g=0;
  if(SzFreq_Nmeth>0)
  {
    SzFreq_HaveObs.initialize();
    SzFreq_HaveObs2.initialize();
    for(k=1;k<=SzFreq_Nmeth;k++)
    {
      SzFreq_Nbins3(k)=gender*SzFreq_Nbins(k);
    for(s=1;s<=nseas;s++)
    {
      g++;
      SzFreq_Nbins_seas_g(g)=SzFreq_Nbins(k)*gender;
    }
  }
  }
 END_CALCS

  init_matrix SzFreq_bins1(1,SzFreq_Nmeth,1,SzFreq_Nbins);    // lower edge of wt bins
  !!if(SzFreq_Nmeth>0) echoinput<<" SizeFreq bins-raw "<<endl<<SzFreq_bins1<<endl;
  matrix SzFreq_bins(1,SzFreq_Nmeth,1,SzFreq_Nbins3);     //  szfreq bins as processed and doubled for the males if necessary
  matrix SzFreq_bins2(1,SzFreq_Nmeth,0,SzFreq_Nbins3+1);   //  as above, but one more bin to aid in the search for bin boundaries
  ivector SzFreq_Omit_Small(1,SzFreq_Nmeth);
  int SzFreq_totobs
  int SzFreq_N_Like

//  ProgLabel_2.11.1  size comp bins according to scaling method
 LOCAL_CALCS
  if(SzFreq_Nmeth>0)
  {
  for(k=1;k<=SzFreq_Nmeth;k++)
  {
// set flag for accumulating, or not, fish from small pop len bins up into first SzFreq data bin
    SzFreq_Omit_Small(k)=1;
    if(SzFreq_bins1(k,1)<0)
    {
      SzFreq_Omit_Small(k)=-1;
      SzFreq_bins1(k,1)*=-1;  // make this positive for use in model, then write out as negative in data.ss_new
    }

    SzFreq_bins(k)(1,SzFreq_Nbins(k))=SzFreq_bins1(k)(1,SzFreq_Nbins(k));
    if(gender==2)
    {
      for(j=1;j<=SzFreq_Nbins(k);j++)
      {
        SzFreq_bins(k,j+SzFreq_Nbins(k))=SzFreq_bins1(k,j);
      }
    }
    if(SzFreq_scale(k)==2)  // convert from lbs to kg
    {
      SzFreq_bins(k)*=0.4536;
    }
    else if (SzFreq_scale(k)==4)  // convert from inches to cm
    {
      SzFreq_bins(k)*=2.54;
    }
    SzFreq_bins2(k,0)=0.;
    SzFreq_bins2(k)(1,SzFreq_Nbins(k))=SzFreq_bins(k)(1,SzFreq_Nbins(k));
    if(gender==2)
    {
      SzFreq_bins2(k,SzFreq_Nbins(k)+1)=0.;
      for(j=1;j<=SzFreq_Nbins(k);j++)
      {SzFreq_bins2(k,j+SzFreq_Nbins(k)+1)=SzFreq_bins2(k,j);}
    }
    echoinput<<"Processed_SizeFreqMethod_bins"<<k<<endl<<SzFreq_bins(k)<<endl;;
  }
  SzFreq_totobs=sum(SzFreq_nobs);
  }
 END_CALCS

//  NOTE:  for the szfreq data, which are stored in one list and not by fleet, it is not possible to exclude from the working array on basis of before styr or after retroyr
  ivector SzFreq_Setup(1,SzFreq_totobs);  //  stores the number of bins plus header info to read into ragged array
  ivector SzFreq_Setup2(1,SzFreq_totobs);   //  stores the number of bins for each obs to create the ragged array

 LOCAL_CALCS
  if(SzFreq_Nmeth>0)
  {
  g=0;
  for(k=1;k<=SzFreq_Nmeth;k++)
  for(j=1;j<=SzFreq_nobs(k);j++)
  {g++; SzFreq_Setup(g)=7+gender*SzFreq_Nbins(k); SzFreq_Setup2(g)=gender*SzFreq_Nbins(k);}
  }
 END_CALCS

//  ProgLabel_2.11.2 read size comp observations into a ragged array
// , with the number of elements for each obs stored in sizefreq_setup
//   unlike the size and agecomp, obs from all fleets are in one dimension, rather than having a dimension for fleet
//  to do super-period, obs must be sorted by fleet and time within each method
  init_matrix SzFreq_obs1(1,SzFreq_totobs,1,SzFreq_Setup);
   !!if(SzFreq_totobs>0) echoinput<<" first sizefreq obs "<<endl<<SzFreq_obs1(1)<<endl<<" last obs"<<endl<<SzFreq_obs1(SzFreq_totobs)<<endl;;
  imatrix SzFreq_obs_hdr(1,SzFreq_totobs,1,8);
  // SzFreq_obs1:     Method, Year, season, Fleet, Gender, Partition, SampleSize, <data>
  // SzFreq_obs_hdr:           1=y; 2=s;    3=f;   4=gender; 5=partition; 6=method&skip flag; 7=first bin to use; 8=last bin(e.g. to include males or not)
  vector SzFreq_sampleN(1,SzFreq_totobs);
  vector SzFreq_effN(1,SzFreq_totobs);
  vector SzFreq_eachlike(1,SzFreq_totobs);
  matrix SzFreq_obs(1,SzFreq_totobs,1,SzFreq_Setup2);
  imatrix SzFreq_LikeComponent(1,Ntypes,1,SzFreq_Nmeth)
  number N_suprper_SzFreq   //  no real need to keep track of these by method, so just use a number
 LOCAL_CALCS
  SzFreq_N_Like=0;
  if(SzFreq_Nmeth>0)
  {
    SzFreq_LikeComponent.initialize();
    SzFreq_obs.initialize();
    N_suprper_SzFreq=0;
    iobs=0;
    for(k=1;k<=SzFreq_Nmeth;k++)
    {
      for(j=1;j<=SzFreq_nobs(k);j++)
      {
//       if(y>=styr && y<=retro_yr)  // not used because all obs in one list
        iobs++;
        for(z=1;z<=5;z++)
        {SzFreq_obs_hdr(iobs,z) = int(SzFreq_obs1(iobs,z+1));}
        SzFreq_sampleN(iobs) = SzFreq_obs1(iobs,7);
        if(SzFreq_obs_hdr(iobs,2)<0) N_suprper_SzFreq++;  //  count the number of superperiod start/stops
        if(SzFreq_obs_hdr(iobs,4)==3)  // both genders
        {
          for(z=1;z<=SzFreq_Setup2(iobs);z++) {SzFreq_obs(iobs,z)=SzFreq_obs1(iobs,7+z);}
        }
        else if(SzFreq_obs_hdr(iobs,4)<=1)  // combined gender or female only
        {
          for(z=1;z<=SzFreq_Nbins(k);z++) {SzFreq_obs(iobs,z)=SzFreq_obs1(iobs,7+z);}
        }
        else  // male only
        {
          for(z=SzFreq_Nbins(k)+1;z<=SzFreq_Setup2(iobs);z++) {SzFreq_obs(iobs,z)=SzFreq_obs1(iobs,7+z);}
        }
        if(gender==1) SzFreq_obs_hdr(iobs,4)=1;  // just in case

        SzFreq_obs(iobs)/=sum(SzFreq_obs(iobs));
        SzFreq_obs(iobs)+=SzFreq_mincomp(k);
        SzFreq_obs(iobs)/=sum(SzFreq_obs(iobs));
        y=SzFreq_obs_hdr(iobs,1);
        s=abs(SzFreq_obs_hdr(iobs,2));
          if(s>nseas)
          {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, season for general sizecomp  method, obs "<<k<<" "<<j<<" is > nseas"<<endl; exit(1);}
        t=styr+(y-styr)*nseas+s-1;
        f=abs(SzFreq_obs_hdr(iobs,3));
        if(gender==1) {SzFreq_obs_hdr(iobs,4)=0;}
        z=SzFreq_obs_hdr(iobs,4);  // gender
// get min and max index according to use of 0, 1, 2, 3 gender index
        if(z!=2) {SzFreq_obs_hdr(iobs,7)=1;} else {SzFreq_obs_hdr(iobs,7)=SzFreq_Nbins(k)+1;}
        if(z<=1) {SzFreq_obs_hdr(iobs,8)=SzFreq_Nbins(k);} else {SzFreq_obs_hdr(iobs,8)=2*SzFreq_Nbins(k);}
  //      SzFreq_obs_hdr(iobs,5);  // partition
        SzFreq_obs_hdr(iobs,6)=k;
        if(k!=SzFreq_obs1(iobs,1)) {N_warn++; warning<<" sizefreq ID # doesn't match for obs "<<iobs<<"; Expected: "<<k<<" Read: "<<SzFreq_obs1(iobs,1)<<endl; }

        if(y>=styr && y<=retro_yr)
        {
          SzFreq_LikeComponent(f,k)=1;    // indicates that this combination is being used
          if(SzFreq_HaveObs(f,k,t,1)==0) SzFreq_HaveObs(f,k,t,1)=iobs;  // save first counter in time x fleet locations with data
          SzFreq_HaveObs(f,k,t,2)=iobs;  // saves last pointer to this source of data
          if(SzFreq_HaveObs2(k,t)==0 || f<=SzFreq_HaveObs2(k,t)) SzFreq_HaveObs2(k,t)=f;  // find the smallest numbered f index that uses this method
          have_data(f,t)=1;
          if(SzFreq_obs_hdr(iobs,7)<0) SzFreq_obs_hdr(iobs,3)=-abs(SzFreq_obs_hdr(iobs,3));  //  old method for excluding from logL
        }
        else
        {
          SzFreq_obs_hdr(iobs,3)=-abs(SzFreq_obs_hdr(iobs,3));  //  flag for skipping this obs
        }
      }
    }
    SzFreq_N_Like=sum(SzFreq_LikeComponent);
    if(N_suprper_SzFreq>0)
    {
      j=N_suprper_SzFreq/2;  // because we counted the begin and end
      if(2*j!=N_suprper_SzFreq)
      {
        N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"unequal number of starts and ends of sizefreq superperiods "<<endl; exit(1);
      }
      else
      {
        N_suprper_SzFreq=j;
      }
      echoinput<<"N superperiods for sizecomp "<<N_suprper_SzFreq<<endl;
    }
  }
 END_CALCS

//  ProgLabel_2.11.3  calc logL for a perfect fit to the sizefreq data as an offset
    vector SzFreq_like_base(1,SzFreq_N_Like)  // this is a constant offset, so can be declared in data section
    ivector suprper_SzFreq_start(1,N_suprper_SzFreq)
    ivector suprper_SzFreq_end(1,N_suprper_SzFreq)
    vector suprper_SzFreq_sampwt(1,SzFreq_totobs)  //  will contain calculated weights for obs within super periods

 LOCAL_CALCS
  if(SzFreq_Nmeth>0)
  {
    SzFreq_like_base.initialize();
    suprper_SzFreq_start.initialize();
    suprper_SzFreq_end.initialize();
    suprper_SzFreq_sampwt.initialize();
    N_suprper_SzFreq=0;  // redo this counter so can use the counter
//  count the number of type x methods being used to create vector length for the likelihoods
    g=0;
    for(f=1;f<=Ntypes;f++)
    for(k=1;k<=SzFreq_Nmeth;k++)
    {
      if(SzFreq_LikeComponent(f,k)>0) {g++; SzFreq_LikeComponent(f,k)=g;}  //  so stored value g gives index in list of logL elements
    }
    in_superperiod=0;
    for(iobs=1;iobs<=SzFreq_totobs;iobs++)
    {
      k=SzFreq_obs_hdr(iobs,6);  //  get the method
      f=abs(SzFreq_obs_hdr(iobs,3));
      s=SzFreq_obs_hdr(iobs,2);  // sign used to indicate start/stop of super period
      if(SzFreq_obs_hdr(iobs,3)>0)  // negative for out of range or skip
      {
        z1=SzFreq_obs_hdr(iobs,7);
        z2=SzFreq_obs_hdr(iobs,8);
        g=SzFreq_LikeComponent(f,k);
        SzFreq_like_base(g)-=SzFreq_sampleN(iobs)*SzFreq_obs(iobs)(z1,z2)*log(SzFreq_obs(iobs)(z1,z2));
      }

// identify super-period starts and stops
      if(s<0) // start/stop a super-period  ALL observations must be continguous in the file
      {
        if(in_superperiod==0) 
        {
          N_suprper_SzFreq++;
          suprper_SzFreq_start(N_suprper_SzFreq)=iobs;
          in_superperiod=1;
        }
        else if(in_superperiod==1)  // end a super-period
        {
          suprper_SzFreq_end(N_suprper_SzFreq)=iobs;
          in_superperiod=0;
        }
      }
    }
  }
  echoinput<<" finished processing sizefreq data "<<endl;
  if(N_suprper_SzFreq>0) echoinput<<"sizefreq superperiod start obs: "<<suprper_SzFreq_start<<endl<<"sizefreq superperiod end obs:   "<<suprper_SzFreq_end<<endl;
 END_CALCS

//  ProgLabel_2.12  read tag release and recapture data
  init_int Do_TG
  !!echoinput<<Do_TG<<" Do_TagData(0/1) "<<endl;

  init_vector TG_temp(1,4*Do_TG)
  int TG;
  int N_TG   // N tag groups
  int N_TG2;
  int TG_timestart;
  int N_TG_recap;   //  N recapture events
  int TG_mixperiod; //  First period (seasons) to start comparing obs to expected recoveries; period=0 is the release period
  int TG_maxperiods; //  max number of periods (seasons) to track recoveries; period=0 is the release period
 LOCAL_CALCS
  if(Do_TG>0)
  {
    Do_TG=1;
    N_TG=TG_temp(1);
    N_TG_recap=TG_temp(2);
    TG_mixperiod=TG_temp(3);
    TG_maxperiods=TG_temp(4);
    N_TG2=N_TG;
    TG_timestart=9999;
    echoinput<<N_TG<<" N tag groups "<<endl
    <<N_TG_recap<<" N recapture events"<<endl
    <<TG_mixperiod<<"  Latency period for mixing"<<endl
    <<TG_maxperiods<<" N periods to track recoveries"<<endl;
  }
  else
  {
    N_TG=0;
    N_TG_recap=0;
    TG_mixperiod=0;
    TG_maxperiods=0;
    N_TG2=1;
    TG_timestart=1;
  }
 END_CALCS

  ivector TG_endtime(1,N_TG2)
  init_matrix TG_release(1,N_TG,1,8)
  // TG area  year season tindex gender age N_released
 LOCAL_CALCS
   echoinput<<" Tag Releases "<<endl<<"TG area year seas tindex gender age N_released "<<endl<<TG_release<<endl;
   TG_endtime(1)=0;
   if(N_TG>0)
   {
   for(TG=1;TG<=N_TG;TG++)
   {
     t=styr+int((TG_release(TG,3)-styr)*nseas+TG_release(TG,4)-1);
     TG_release(TG,5)=t;
     if(t<TG_timestart) TG_timestart=t;
     k=TG_maxperiods;
     if((t+TG_maxperiods)>TimeMax) k-=(t+TG_maxperiods-TimeMax);
      TG_endtime(TG)=k;
   }
  }
 END_CALCS

//  ProgLabel_2.12.1  store recapture info by TG group and time to follow it as a cohort
  init_matrix TG_recap_data(1,N_TG_recap,1,5)
  //  TG, year, season, fleet, gender, Number
  3darray TG_recap_obs(1,N_TG2,0,TG_endtime,0,Nfleet);   //  no area index because each fleet is in just one area
 LOCAL_CALCS
   if(N_TG>0)
   {
   echoinput<<"First row of tag-recapture data "<<TG_recap_data(1)<<endl;
   echoinput<<"Last  row of tag-recapture data "<<TG_recap_data(N_TG_recap)<<endl;
   TG_recap_obs.initialize();
   for(j=1;j<=N_TG_recap;j++)
   {
     TG=TG_recap_data(j,1);  // TD the tag group
     t=styr+int((TG_recap_data(j,2)-styr)*nseas+TG_recap_data(j,3)-1) - TG_release(TG,5); // find elapsed time in terms of number of seasons
     if(t>TG_maxperiods) t=TG_maxperiods;
     if(t<0)
     {
       N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" recapture is before tag release for recap: "<<j<<endl;  exit(1);
     }
     TG_recap_obs(TG,t,TG_recap_data(j,4))+=TG_recap_data(j,5);  //   save N recaptures by TG, fleet of recapture, elapsed time
   }
   for(TG=1;TG<=N_TG;TG++)
   {
     for(TG_t=0;TG_t<=TG_endtime(TG);TG_t++)
     {
       TG_recap_obs(TG,TG_t,0) = sum(TG_recap_obs(TG,TG_t)(1,Nfleet));
       if(TG_recap_obs(TG,TG_t,0)>0.) TG_recap_obs(TG,TG_t)(1,Nfleet)/=TG_recap_obs(TG,TG_t,0);
     }
   }
  }
 END_CALCS

//  ProgLabel_2.13  Morphcomposition data
   init_int Do_Morphcomp
  !!echoinput<<Do_Morphcomp<<" Do_Morphcomp(0/1) "<<endl;
   init_vector mc_temp(1,3*Do_Morphcomp);
   int Morphcomp_nobs
   int Morphcomp_nmorph
   number Morphcomp_mincomp
 LOCAL_CALCS
  if(Do_Morphcomp>0)
  {
    Morphcomp_nobs=mc_temp(1);
    Morphcomp_nmorph=mc_temp(2);   // later compare this value to the n morphs in the control file and exit if different
    Morphcomp_mincomp=mc_temp(3);
  echoinput<<Morphcomp_nobs<<" Morphcomp_nobs "<<endl;
  echoinput<<Morphcomp_nmorph<<" Morphcomp_nmorph "<<endl;
  echoinput<<Morphcomp_mincomp<<" Morphcomp_mincomp "<<endl;
  }
  else
  {
    Morphcomp_nobs=0;
  }
 END_CALCS
 init_matrix Morphcomp_obs(1,Morphcomp_nobs,1,5+Morphcomp_nmorph)
//    yr, seas, type, partition, Nsamp, datavector

  3darray Morphcomp_havedata(1,Ntypes*Do_Morphcomp,styr,TimeMax,0,0)    // last dimension is reserved for future use of Partition
 LOCAL_CALCS
  if(Do_Morphcomp>0)
  {
  echoinput<<" morph composition data"<<endl<<"yr seas type partition Nsamp datavector"<<endl<< Morphcomp_obs<<endl;
  Morphcomp_havedata=0;
  for(i=1;i<=Morphcomp_nobs;i++)
  {
    y=Morphcomp_obs(i,1); s=Morphcomp_obs(i,2); t=styr+(y-styr)*nseas+s-1;
    f=Morphcomp_obs(i,3); z=Morphcomp_obs(i,4);   // z not used, partition must be 0 (e.g. combined discard and retained)
    Morphcomp_havedata(f,t,0)=i;
    have_data(f,t)=1;
    if(y>retro_yr) Morphcomp_obs(i,5) = -fabs(Morphcomp_obs(i,5));
    Morphcomp_obs(i)(6,5+Morphcomp_nmorph) /= sum(Morphcomp_obs(i)(6,5+Morphcomp_nmorph));
    Morphcomp_obs(i)(6,5+Morphcomp_nmorph) += Morphcomp_mincomp;
    Morphcomp_obs(i)(6,5+Morphcomp_nmorph) /= sum(Morphcomp_obs(i)(6,5+Morphcomp_nmorph));
  }
  }
 END_CALCS

//  ProgLabel_2.14  End of datafile indicator
  init_int fid
  !! if(fid!=999) {cout<<" final data value in error "<<fid<<endl; exit(1);}
  !! cout<<"Data read sucessful "<<fid<<endl<<endl;
  !!echoinput<<" data read successful"<<endl<<endl;

//**********************************************************************************************************************
//  ProgLabel_3.0   read forecast.ss
 LOCAL_CALCS
  ad_comm::change_datafile_name("forecast.ss");
  cout<<" reading forecast file "<<endl;
  ifstream Forecast_Stream("forecast.ss");   //  even if the global_datafile name is used, there still is a different logical device created
  k=0;
  N_FC=0;
  while(k==0)
  {
    Forecast_Stream >>  readline;          // reads the line from input stream
    if(length(readline)>2)
    {
      checkchar=readline(1);
      k=strcmp(checkchar,"#");
      checkchar=readline(1,2);
      j=strcmp(checkchar,"#C");
      if(j==0) {N_FC++; Forecast_Comments+=readline;}
    }
  }
 END_CALCS

  init_int Do_Benchmark  // 0=skip; do Fspr, Fbtgt, Fmsy
  !!echoinput<<Do_Benchmark<<" Do_Benchmark "<<endl;
  init_int Do_MSY   //  1= set to F(SPR); 2=calc F(MSY); 3=set to F(Btgt); 4=set to F(endyr)
  !!echoinput<<Do_MSY<<" Do_MSY "<<endl;

  int did_MSY;
  int show_MSY;
  init_number SPR_target
  !!echoinput<<SPR_target<<" SPR_target "<<endl;
  init_number BTGT_target
  !!echoinput<<BTGT_target<<" BTGT_target "<<endl;

  ivector Bmark_Yr(1,6)
  ivector Bmark_t(1,2)  //  for range of time values for averaging body size
  init_ivector Bmark_Yr_rd(1,6)
  init_int Bmark_RelF_Basis
 LOCAL_CALCS
  echoinput<<Bmark_Yr_rd<<" Benchmark years as read:  beg-end bio; beg-end selex; beg-end relF"<<endl;
  for(i=1;i<=6;i++)  //  beg-end bio; beg-end selex; beg-end relF
  {
    if(Bmark_Yr_rd(i)==-999)
    {Bmark_Yr(i)=styr;}
    else if(Bmark_Yr_rd(i)>=styr)
    {Bmark_Yr(i)=Bmark_Yr_rd(i);}
    else if(Bmark_Yr_rd(i)<=0)
    {Bmark_Yr(i)=endyr+Bmark_Yr_rd(i);}
    else
    {
      N_warn++;Bmark_Yr(i)=styr;warning<<"benchmark year less than styr; reset to equal styr"<<endl;
    }
  }
  Bmark_t(1)=styr+(Bmark_Yr(1)-styr)*nseas;
  Bmark_t(2)=styr+(Bmark_Yr(2)-styr)*nseas;

  echoinput<<Bmark_Yr<<" Benchmark years as processed:  beg-end bio; beg-end selex; beg-end relF"<<endl;
  echoinput<<Bmark_RelF_Basis<<"  1=use range of years for relF; 2 = set same as forecast relF below"<<endl;
 END_CALCS

  init_int Do_Forecast   //  0=none; 1=F(SPR); 2=F(MSY) 3=F(Btgt); 4=Ave F (enter yrs); 5=read Fmult
  !!echoinput<<Do_Forecast<<" Do_Forecast "<<endl;

  vector Fcast_Input(1,22);

  int N_Fcast_Yrs
  ivector Fcast_yr(1,4)  // yr range for selex, then yr range foreither allocation or for average F
  int Fcast_Sel_yr1
  int Fcast_Sel_yr2
  int Fcast_RelF_yr1
  int Fcast_RelF_yr2
  int Fcast_RelF_Basis  // 1=use year range; 2=read below
  number Fcast_Flevel
  int Do_Rebuilder
  int Rebuild_Ydecl
  int Rebuild_Yinit
  int HarvestPolicy  // 1=west coast adjust catch; 2=AK to adjust F
  number H4010_top
  number H4010_bot
  number H4010_scale
  int Do_Impl_Error
  number Impl_Error_Std
  ivector Fcast_Loop_Control(1,5)
  int N_Fcast_Input_Catches
  int Fcast_InputCatch_Basis  //  2=dead catch; 3=retained catch;  99=F; -1=read fleet/time specific  (biomass vs numbers will match catchunits(fleet)
  int Fcast_Catch_Basis  //  2=dead catch bio, 3=retained catch bio, 5= dead catch numbers 6=retained catch numbers;   Same for all fleets

  int Fcast_Catch_Allocation_Groups;
  int Fcast_Do_Fleet_Cap;
  int Fcast_Do_Area_Cap;
  int Fcast_Cap_FirstYear;

//  matrix Fcast_RelF(1,nseas,1,Nfleet)

 LOCAL_CALCS
  if(Do_Forecast==0)
  {
    k=0;
    echoinput<<"No forecast selected, so rest of forecast file will not be read and can be omitted"<<endl;
    if(Bmark_RelF_Basis==2) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Fatal stop:  no forecast, but bmark set to use fcast"<<endl;  exit(1);}
  }
  else
  {
    k=22;
    echoinput<<"Forecast selected; next 20 input values will be read as a block then parsed and procesed "<<endl;
  }
 END_CALCS
  init_vector Fcast_Input_rd(1,k)

 LOCAL_CALCS
  if(Do_Forecast>0)
  {
    Fcast_Input(1,k)=Fcast_Input_rd(1,k);
  k=0;
  k++; N_Fcast_Yrs=int(Fcast_Input(k));
  echoinput<<N_Fcast_Yrs<<" N_Fcast_Yrs "<<endl;
  if(N_Fcast_Yrs<=0) {N_warn++; cout<<"Critical error in forecast input, see warning"<<endl; warning<<"ERROR: cannot do a forecast of zero years"<<endl; exit(1);}
  YrMax=endyr+N_Fcast_Yrs;
  k++; Fcast_Flevel=Fcast_Input(k);
  echoinput<<Fcast_Flevel<<" Fmult value used only if Do_Forecast==5"<<endl;
  k++; Fcast_yr(1)=int(Fcast_Input(k));
  k++; Fcast_yr(2)=int(Fcast_Input(k));
  k++; Fcast_yr(3)=int(Fcast_Input(k));
  k++; Fcast_yr(4)=int(Fcast_Input(k));
  echoinput<<Fcast_yr<<" Begin-end yrs for average selex; begin-end yrs for allocation"<<endl;
  for(i=1;i<=4;i++)
  {
    if(Fcast_yr(i)==-999)
    {Fcast_yr(i)=styr;}
    else if(Fcast_yr(i)>0)
    {Fcast_yr(i)=Fcast_yr(i);}
    else
    {Fcast_yr(i)=endyr+Fcast_yr(i);}
  }
  Fcast_Sel_yr1=Fcast_yr(1);
  Fcast_Sel_yr2=Fcast_yr(2);
  Fcast_RelF_yr1=Fcast_yr(3);
  Fcast_RelF_yr2=Fcast_yr(4);
  echoinput<<Fcast_yr<<"  After Transformation:  begin-end yrs for average selex; begin-end yrs for rel F"<<endl;

  k++; HarvestPolicy=int(Fcast_Input(k));
  k++; H4010_top=Fcast_Input(k);
  k++; H4010_bot=Fcast_Input(k);
  k++; H4010_scale=Fcast_Input(k);
  echoinput<<HarvestPolicy<<"  HarvestPolicy "<<endl<<
  H4010_top<<"  H4010_top "<<endl<<
  H4010_bot<<"  H4010_bot "<<endl<<
  H4010_scale<<"  H4010_scale "<<endl;
  if(H4010_top<=H4010_bot) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" control rule top: "<<H4010_top<<" must be > control rule bottom "<<H4010_bot<<endl; exit(1);}
  if(H4010_scale>1.0) {N_warn++; warning<<" Sure you want harvest policy scalar > 1.0? "<<H4010_scale<<endl;}

  k++; Fcast_Loop_Control(1)=Fcast_Input(k);
  k++; Fcast_Loop_Control(2)=Fcast_Input(k);
  k++; Fcast_Loop_Control(3)=Fcast_Input(k);
  k++; Fcast_Loop_Control(4)=Fcast_Input(k);
  k++; Fcast_Loop_Control(5)=Fcast_Input(k);
  k++; Fcast_Cap_FirstYear=Fcast_Input(k);
  echoinput<<Fcast_Loop_Control(1)<<" N forecast loops (1-3) (recommend 3)"<<endl;
  echoinput<<Fcast_Loop_Control(2)<<" First forecast loop with stochastic recruitment (recommend 3)"<<endl;
  echoinput<<Fcast_Loop_Control(3)<<" Forecast loop control #3 (reserved for future use) "<<endl;
  echoinput<<Fcast_Loop_Control(4)<<" Forecast loop control #4 (reserved for future use) "<<endl;
  echoinput<<Fcast_Loop_Control(5)<<" Forecast loop control #5 (reserved for future use) "<<endl;
  echoinput<<Fcast_Cap_FirstYear<<"  Fcast_Cap_FirstYear for caps and allocations (should be after any fixed inputs) "<<endl;

  k++; Impl_Error_Std=Fcast_Input(k);
  echoinput<<Impl_Error_Std<<"  Impl_Error_Std "<<endl;
  if(Impl_Error_Std>0.0) Do_Impl_Error=1;

  k++; Do_Rebuilder=int(Fcast_Input(k));
  k++; Rebuild_Ydecl=int(Fcast_Input(k));
  k++; Rebuild_Yinit=int(Fcast_Input(k));
  echoinput<<
  Do_Rebuilder<<"  Do_Rebuilder "<<endl<<
  Rebuild_Ydecl<<"  Rebuild_Ydecl "<<endl<<
  Rebuild_Yinit<<"  Rebuild_Yinit "<<endl;

  k++; Fcast_RelF_Basis=int(Fcast_Input(k));
  echoinput<<Fcast_RelF_Basis<<" fleet relative F (1=use_year_range/2=read_seas_x_fleet) "<<endl;

  k++; Fcast_Catch_Basis=Fcast_Input(k);
    echoinput<<Fcast_Catch_Basis<<"  Fcast_Catch_Basis for caps and allocations "<<endl;
    echoinput<<"2=dead catch bio, 3=retained catch bio, 5= dead catch numbers 6=retained catch numbers;   Same for all fleets"<<endl;
  }
  else
  {
    //  placeholders here until bio_ranges get coded
  Fcast_RelF_Basis=1;
  Fcast_Sel_yr1=endyr;
  Fcast_Sel_yr2=endyr;
  Fcast_RelF_yr1=endyr;
  Fcast_RelF_yr2=endyr;
  N_Fcast_Yrs=0;
  YrMax=endyr;
  Do_Rebuilder=0;
  }
  if(Fcast_RelF_Basis==2)
  {
    if(Do_Forecast==4)
    {
      N_warn++; warning<<"Cannot specify forecast fleet relative F because Do_Forecast==4 specifies relative F directly as F;"<<endl;
      Fcast_RelF_Basis=1;
    }
    z=nseas;
    echoinput<<"Fcast_RelF_Basis==2, so now read seas(row) x fleet (column) array of relative F; will be re-scaled to sum to 1.0"<<endl;
  }
  else
  {z=0;}
 END_CALCS

  init_matrix Fcast_RelF_Input(1,z,1,Nfleet)

 LOCAL_CALCS
  if(Fcast_RelF_Basis==2)
  {
    echoinput<<" fleet relative F by season and fleet as read"<<endl<<Fcast_RelF_Input<<endl;
    // later set Fcast_RelF_Use=Fcast_RelF_Input;
  }
  else
  {
    echoinput<<" do not include relative F by season and fleet, uses endyr instead "<<endl;
  }
  if(Do_Forecast>0) {k=1;} else {k=0;}
  if(Do_Forecast>0)
  {
    echoinput<<"Now read cap for each fleet, then cap for each area (even if only 1 area), then allocation assignment for each fleet"<<endl;
  }
 END_CALCS

  init_vector Fcast_MaxFleetCatch(1,k*Nfleet)
  init_vector Max_Fcast_Catch(1,k*pop)
  init_ivector Allocation_Fleet_Assignments(1,k*Nfleet)

 LOCAL_CALCS
  Fcast_Do_Fleet_Cap=0;
  Fcast_Do_Area_Cap=0;
  Fcast_Catch_Allocation_Groups=0;
  if(k>0)
  {
    for(f=1;f<=Nfleet;f++)
    {if(Fcast_MaxFleetCatch(f)>0.0) Fcast_Do_Fleet_Cap=1;}
    for(p=1;p<=pop;p++)
    {if(Max_Fcast_Catch(p)>0.0) Fcast_Do_Area_Cap=1;}
    Fcast_Catch_Allocation_Groups=max(Allocation_Fleet_Assignments);
  }
 END_CALCS

  init_vector Fcast_Catch_Allocation(1,Fcast_Catch_Allocation_Groups);

 LOCAL_CALCS
  if(Do_Forecast>0)
  {
    k=2;
    echoinput<<" Max totalcatch by fleet "<<endl<<Fcast_MaxFleetCatch<<endl;
    echoinput<<" Max totalcatch by area "<<endl<<Max_Fcast_Catch<<endl;
    echoinput<<" Assign fleets to allocation groups (0 means not in a group) "<<endl<<Allocation_Fleet_Assignments<<endl;
    echoinput<<" calculated number of allocation groups "<<Fcast_Catch_Allocation_Groups<<endl;
    echoinput<<" Allocation among groups (N entries must match number of allocation groups created) "<<Fcast_Catch_Allocation<<endl;
  }
  else
  {k=0;}
 END_CALCS

  init_ivector more_Fcast_input(1,k);
 LOCAL_CALCS
  if(k>0)
  {
    N_Fcast_Input_Catches=more_Fcast_input(1);
    Fcast_InputCatch_Basis=more_Fcast_input(2);
    echoinput<<N_Fcast_Input_Catches<<" N_Fcast_input_catches "<<endl;
    echoinput<<Fcast_InputCatch_Basis<<" # basis for input Fcast catch:  2=dead catch; 3=retained catch; 99=input Hrate(F); -1=read fleet/time specific (bio/num units are from fleetunits; note new codes in SSV3.20)"<<endl;
    k1 = styr+(endyr-styr)*nseas-1 + nseas + 1;
    y=k1+(N_Fcast_Yrs)*nseas-1;
    if(N_Fcast_Input_Catches>0) echoinput<<"Now read "<<N_Fcast_Input_Catches<<" of fixed forecast catches (yr, seas, fleet, (basis,) catch_or_F) "<<endl;
  }
  else
  {
    k1=1;
    y=0;
  }
  if(Fcast_InputCatch_Basis==-1)
    {j=5;}
    else
    {j=4;}
 END_CALCS

  3darray Fcast_InputCatch(k1,y,1,Nfleet,1,2)  //  values and basis to be used
  init_matrix Fcast_InputCatch_rd(1,N_Fcast_Input_Catches,1,j)  //  values to be read:  yr, seas, fleet, value, (basis) 

 LOCAL_CALCS
  if(N_Fcast_Input_Catches>0) echoinput<<" Fcast_catches_input "<<endl<<Fcast_InputCatch_rd<<endl;
  if(Do_Forecast>0)
  {
    for(t=k1;t<=y;t++)
    for(f=1;f<=Nfleet;f++)
    {Fcast_InputCatch(t,f,1)=-1;}
    if(N_Fcast_Input_Catches>0)
    {
      for(i=1;i<=N_Fcast_Input_Catches;i++)
      {
        y=Fcast_InputCatch_rd(i,1); s=Fcast_InputCatch_rd(i,2); f=Fcast_InputCatch_rd(i,3);
        if(y>endyr && y<=YrMax && f<=Nfleet)
        {
          t=styr+(y-styr)*nseas +s-1;
          Fcast_InputCatch(t,f,1)=Fcast_InputCatch_rd(i,4);
          if(y>=Fcast_Cap_FirstYear) {N_warn++;warning<<"Input catches in "<<y<<" can be overridden by caps or allocations"<<endl;}
          if(Fcast_InputCatch_Basis==-1)
          {
            Fcast_InputCatch(t,f,2)=Fcast_InputCatch_rd(i,5);  //  new method
          }
          else
          {
            Fcast_InputCatch(t,f,2)=Fcast_InputCatch_Basis;  //  method before 3.24P
          }
      	}
      }
    }
    echoinput<<"Processed "<<endl<<Fcast_InputCatch<<endl;
  }

  if(Do_Rebuilder>0 && Do_Forecast<=0) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" error: Rebuilder output selected without requesting forecast"<<endl; exit(1);}
  if(Do_Benchmark==0)
  {
    if(Do_Forecast>=1 && Do_Forecast<=3) {Do_Benchmark=1; N_warn++; warning<<" Turn Benchmark on because Forecast needs it"<<endl;}
    if(Do_Forecast==0 && F_report_basis>0) {F_report_basis=0; N_warn++; warning<<" Set F_report_basis=0 because no benchmark or forecast"<<endl;}
    if(depletion_basis==2) {depletion_basis=1; N_warn++; warning<<" Change depletion basis to 1 because benchmarks are off"<<endl;}
    if(SPR_reporting>=1 && SPR_reporting<=3) {SPR_reporting=4; N_warn++; warning<<" Change SPR_reporting to 4 because benchmarks are off"<<endl;}
  }
  else
  {
     if(Do_MSY==0)  {Do_MSY=1; N_warn++; warning<<" Setting Do_MSY=1 because benchmarks are on"<<endl;}
  }
  if(Do_Forecast==2 && Do_MSY!=2) {Do_MSY=2; N_warn++; warning<<" Set MSY option =2 because Forecast option =2"<<endl;}
  if(depletion_basis==2 && Do_MSY!=2) {Do_MSY=2; N_warn++; warning<<" Set MSY option =2 because depletion basis is B_MSY"<<endl;}
  if(SPR_reporting==2 && Do_MSY!=2) {Do_MSY=2; N_warn++; warning<<" Set MSY option =2 because SPR basis is SPR_MSY"<<endl;}
  if(Fcast_Sel_yr1>Fcast_Sel_yr2) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Error, Fcast_Sel_Yr1 must be at or before Fcast_Sel_Yr1"<<endl;  exit(1);}
  if(Fcast_Sel_yr1>endyr || Fcast_Sel_yr1<styr) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Error, Fcast_Sel_Yr1 must be between styr and endyr"<<endl;  exit(1);}
  if(Fcast_Sel_yr2>endyr || Fcast_Sel_yr2<styr) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Error, Fcast_Sel_Yr2 must be between styr and endyr"<<endl;  exit(1);}

  did_MSY=0;
  TimeMax_Fcast_std = styr+(YrMax-styr)*nseas+nseas-1;

 END_CALCS

  init_number fif

 LOCAL_CALCS
  if(Do_Forecast>0 && fif!=999) {cout<<" EXIT, must have 999 to verify end of forecast inputs "<<fif<<endl; exit(1);}
  echoinput<<" done reading forecast "<<endl<<endl;
 END_CALCS

//  ProgLabel_3.2  create complete list of years for STD reporting
  ivector STD_Yr_Reverse(styr-2,YrMax)
  ivector STD_Yr_Reverse_Dep(styr-2,YrMax)
  ivector STD_Yr_Reverse_Ofish(styr-2,YrMax)
  ivector STD_Yr_Reverse_F(styr-2,YrMax)
  int N_STD_Yr_Dep;
  int N_STD_Yr_Ofish;
  int N_STD_Yr_F;
  int N_STD_Mgmt_Quant;

 LOCAL_CALCS
  if(STD_Yr_min<0 || STD_Yr_min<(styr-2) ) STD_Yr_min=styr-2;
  if(STD_Yr_max==-1) STD_Yr_max=endyr;
  if(STD_Yr_max==-2) STD_Yr_max=YrMax;
  if(STD_Yr_max>YrMax) STD_Yr_max=YrMax;
   STD_Yr_Reverse.initialize();
   for(y=STD_Yr_min;y<=STD_Yr_max;y++) {STD_Yr_Reverse(y)=1;}
   for(i=1;i<=N_STD_Yr_RD;i++)
   {if(STD_Yr_RD(i)>=styr && STD_Yr_RD(i)<YrMax) {STD_Yr_Reverse(STD_Yr_RD(i))=1;}}

   N_STD_Yr=sum(STD_Yr_Reverse);
 END_CALCS

 LOCAL_CALCS
  STD_Yr_Reverse_Dep.initialize();
  STD_Yr_Reverse_Ofish.initialize();
  STD_Yr_Reverse_F.initialize();
  j=0;
  N_STD_Yr_Dep=0;
  N_STD_Yr_Ofish=0;
  N_STD_Yr_F=0;

  for(y=styr-2;y<=YrMax;y++)
  {
    if(STD_Yr_Reverse(y)>0)
    {
      j++;
      STD_Yr_Reverse(y)=j;  // use for SPB and recruitment
      if(y>=styr)
      {
      // depletion must start in year AFTER first catch.  It could vary earlier if recdevs happened enough earlier to change spbio, but this is not included
      if((depletion_basis>0 && y>first_catch_yr) || y==endyr) {N_STD_Yr_Dep++; STD_Yr_Reverse_Dep(y) = N_STD_Yr_Dep; }
      if(y<=endyr)
      {
        if((SPR_reporting>0 && totcat(y)>0.0) || y==endyr) {N_STD_Yr_Ofish++; STD_Yr_Reverse_Ofish(y) = N_STD_Yr_Ofish; }
        if((F_reporting>0 && totcat(y)>0.0) || y==endyr) {N_STD_Yr_F++; STD_Yr_Reverse_F(y) = N_STD_Yr_F; }
      }
      else
      {
        if(SPR_reporting>0) {N_STD_Yr_Ofish++; STD_Yr_Reverse_Ofish(y) = N_STD_Yr_Ofish; }
        if(F_reporting>0) {N_STD_Yr_F++; STD_Yr_Reverse_F(y) = N_STD_Yr_F; }
      }

      }
    }
  }
 END_CALCS

//********CONTROLS********************************
//  ProgLabel_4.0   Begin Reading from Control File

 LOCAL_CALCS
  ad_comm::change_datafile_name(ctlfilename);
  echoinput<<endl<<" Begin reading control file "<<endl;
  cout<<" reading from control  file"<<endl;
  ifstream Control_Stream(ctlfilename);   //  even if the global_datafile name is used, there still is a different logical device created

//  ProgLabel_4.1  Read and store comments at top of control file
  k=0;
  N_CC=0;
  while(k==0)
  {
    Control_Stream >>  readline;          // reads the line from input stream
    if(length(readline)>2)
    {
      checkchar=readline(1);
      k=strcmp(checkchar,"#");
      checkchar=readline(1,2);
      j=strcmp(checkchar,"#C");
      if(j==0) {N_CC++; Control_Comments+=readline;}
    }
  }
 END_CALCS

//  ProgLabel_4.2  Read info for growth patterns, gender, birthseasons, sub-morphs
  init_int N_GP  // number of growth patterns
  !!echoinput<<N_GP<<" N growth patterns "<<endl;
  init_int N_platoon  //  number of platoons  1, 3, 5 are best values to use
  !!echoinput<<N_platoon<<"  N platoons (1, 3 or 5)"<<endl;
  number sd_ratio;  // ratio of stddev within platoon to between morphs
  number sd_within_platoon
  number sd_between_platoon
  ivector ishadow(1,N_platoon)
  vector shadow(1,N_platoon)
  vector submorphdist(1,N_platoon);
 LOCAL_CALCS
  if(N_platoon>1)
  {
    *(ad_comm::global_datafile) >> sd_ratio;
    *(ad_comm::global_datafile) >> submorphdist;
  echoinput<<sd_ratio<<"  sd_ratio"<<endl;
  echoinput<<submorphdist<<"  submorphdist"<<endl;
  }
  else
  {
    sd_ratio=1.;
    submorphdist(1)=1.;
    echoinput<<"  do not read sd_ratio or submorphdist"<<endl;
  }
//  ProgLabel_4.2.1  assign distribution among growth platoons if needed
  if(submorphdist(1)<0.)
  {
    if(N_platoon==1)
      {submorphdist(1)=1.;}
    else if (N_platoon==3)
      {submorphdist.fill("{0.15,0.70,0.15}");}
    else if (N_platoon==5)
      {submorphdist.fill("{0.031, 0.237, 0.464, 0.237, 0.031}");}
  }
  submorphdist/=sum(submorphdist);

  if(N_platoon>1)
  {
    sd_within_platoon = sd_ratio * sqrt(1. / (1. + sd_ratio*sd_ratio));
    sd_between_platoon = sqrt(1. / (1. + sd_ratio*sd_ratio));
  }
  else
  {sd_within_platoon=1; sd_between_platoon=0.000001;}

   if(N_platoon==1)
     {ishadow(1)=0; shadow(1)=0.;}
   else if (N_platoon==3)
     {ishadow.fill_seqadd(-1,1); shadow.fill_seqadd(-1.,1.);}
   else if (N_platoon==5)
     {ishadow.fill_seqadd(-2,1); shadow.fill_seqadd(-2.,1.);}
   else
     {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" illegal N submorphs, must be 1, 3 or 5 "<<N_platoon<<endl; cout<<" exit - see warning "<<endl; exit(1);}

 END_CALCS

//  ProgLabel_4.2.1.1  define indexing vectors to keep track of characteristics of each morph
  int gmorph // number of gmorphs
  int gp
  int gp2
  int birthseas;
  !! gmorph = N_GP*N_platoon*gender*nseas;
  ivector sx(1,gmorph) //  define sex for each growth morph
  ivector GP(1,gmorph)    //  index for main gender*gpattern
  ivector GP3(1,gmorph)   // index for main gender*gpattern*birthseason
  ivector GP4(1,gmorph)   // index for main gpattern
  ivector Bseas(1,gmorph)
  ivector GP2(1,gmorph)  // reverse pointer for submorphs
//  following two containers are used to track which morphs are being used
  ivector use_morph(1,gmorph)
  imatrix TG_use_morph(1,N_TG2,1,gmorph)
  imatrix recr_dist_pattern_2(1,N_GP,1,nseas)   // used in get_growth to flag morphs needing to have growth calculated

  vector azero_G(1,gmorph);  //  time since Jan 1 at beginning of birthseason in which "g" was born
  3darray curr_age_beg(1,gmorph,1,nseas,0,nages);  // constants, so can be dvector
  int first_grow_age;
  3darray curr_age_mid(1,gmorph,1,nseas,0,nages);

  3darray lin_grow_beg(1,gmorph,1,nseas,0,nages)
  3darray lin_grow_mid(1,gmorph,1,nseas,0,nages)
  imatrix first_grow_age2(1,gmorph,1,nseas);

//  ProgLabel_4.2.2  define distribution of recruitment among growth patterns, seasons, areas
  3darray recr_dist_pattern(1,N_GP,1,nseas,1,pop);
  int recr_dist_inx
  int recr_dist_read
 LOCAL_CALCS
  recr_dist_pattern.initialize();
  if(N_GP*nseas*pop==1) {j=0;} else {j=2;}
 END_CALCS
  init_ivector recr_dist_input(1,j)

 LOCAL_CALCS
  if(j==0)
    {recr_dist_inx=0;  recr_dist_read=0; echoinput<<" do not need to read any recr distribution info"<<endl;}
    else
    {recr_dist_read=recr_dist_input(1);  recr_dist_inx=recr_dist_input(2);
      echoinput<<recr_dist_read<<" recr_dist_read"<<endl;
      echoinput<<recr_dist_inx<<" recr_dist_interaction"<<endl;
      }
 END_CALCS

  init_imatrix recr_dist_pattern_1(1,recr_dist_read,1,3)
  imatrix GP_finder(1,N_GP,1,gender)
 LOCAL_CALCS
   if(recr_dist_read>0)
    {
      echoinput<<" recr distribution as read "<<endl<<recr_dist_pattern_1<<endl;
      for(i=1;i<=recr_dist_read;i++)
      {
        gp=recr_dist_pattern_1(i,1); birthseas=recr_dist_pattern_1(i,2); p=recr_dist_pattern_1(i,3);
        recr_dist_pattern(gp,birthseas,p)=1;
      }
    }
    else
      {
        recr_dist_pattern=1;
      }
      echoinput<<" recr distribution as processed "<<endl<<recr_dist_pattern<<endl;

//  ProgLabel_4.2.3  set-up arrays and indexing for growth patterns, gender, birthseasons, platoons
   g=0;
   int g3i;
   g3i=0;
   recr_dist_pattern_2.initialize();
   use_morph.initialize();
   TG_use_morph.initialize();
   for(k=1;k<=gender;k++)
   for(gp=1;gp<=N_GP;gp++)
   for(birthseas=1;birthseas<=nseas;birthseas++)
   {
     g3i++;
     for(gp2=1;gp2<=N_platoon;gp2++)
     {
       g++;
       GP3(g)=g3i;  // track counter for main morphs (pattern x gender x bseas)
       Bseas(g)=birthseas;
       sx(g)=k;
       GP(g)=gp+(k-1)*N_GP;   // counter for pattern x gender
       GP2(g)=gp2; //   reverse pointer to platoon counter
       GP4(g)=gp;  //  counter for growth pattern
       azero_G(g)=azero_seas(birthseas);
       for(p=1;p<=pop;p++)
       {
         if(recr_dist_pattern(gp,birthseas,p)>0.)
         {
           recr_dist_pattern_2(gp,birthseas)=1;
           use_morph(g)=1;
         }
       }
       if(use_morph(g)==1)
       {
         if( (N_platoon==1) || (N_platoon==3 && gp2==2) || (N_platoon==5 && gp2==3) ) GP_finder(gp,k)=g;  // finds g for a given GP and gender and last birstseason
       }
     }
   }

   for(g=1;g<=gmorph;g++)
   for(s=1;s<=nseas;s++)
   {
     curr_age_beg(g,s)=r_ages+azero_seas(s)-azero_G(g);
     curr_age_mid(g,s)=curr_age_beg(g,s)+seasdur_half(s);
//     curr_age3(g,s)=curr_age_beg(g,s)+seasdur(s);
     if(azero_G(g)>azero_seas(s))
     {
       a=0;
       while(curr_age_mid(g,s,a)<0.0)
       {curr_age_mid(g,s,a)=0.0; a++;}
       a=0;
       while(curr_age_beg(g,s,a)<0.0)
       {curr_age_beg(g,s,a)=0.0; a++;}
     }
     a=0;
   }

    if(N_TG>0)
    {
      for(TG=1;TG<=N_TG;TG++)
      {
        for(g=1;g<=gmorph;g++)
        {
          if(TG_release(TG,6)>2) {N_warn++; warning<<" gender for tag groups must be 0, 1 or 2 "<<endl;}
          if(use_morph(g)>0 && (TG_release(TG,6)==0 || TG_release(TG,6)==sx(g))) TG_use_morph(TG,g)=1;
        }
      }
    }
 END_CALCS

//  ProgLabel_4.3  define movement between areas
   int do_migration  //  number of explicit movements to define
   number migr_firstage
   matrix migr_start(1,nseas,1,N_GP)
 LOCAL_CALCS
   migr_firstage=0.0;
   do_migration=0;
   if (pop>1)
   {
      *(ad_comm::global_datafile) >> do_migration;
      echoinput<<do_migration<<" N_migration definitions to read"<<endl;
      if(do_migration>0)
      {
        *(ad_comm::global_datafile) >> migr_firstage;
        echoinput<<migr_firstage<<" migr_firstage"<<endl;
      }
    }
    else
    {
      echoinput<<" only 1 area, so no read of do_migration or migr_firstage "<<endl;
    }
 END_CALCS
   init_matrix move_def(1,do_migration,1,6)   // seas morph source dest minage maxge
   4darray move_pattern(1,nseas,1,N_GP,1,pop,1,pop)
   int do_migr2
   ivector firstBseas(1,N_GP)

 LOCAL_CALCS
    move_pattern.initialize();
    do_migr2=0;
    if(do_migration>0)
    {
      echoinput<<" migration setup "<<endl<<move_def<<endl;
      for(k=1;k<=do_migration;k++)
      {
        s=move_def(k,1); gp=move_def(k,2); p=move_def(k,3); p2=move_def(k,4);
        move_pattern(s,gp,p,p2)=k;   // save index for definition of this pattern to find the right parameters
      }

      k=do_migration;
      for(s=1;s<=nseas;s++)
      for(gp=1;gp<=N_GP;gp++)
      for(p=1;p<=pop;p++)
      {
        if(move_pattern(s,gp,p,p)==0) {k++; move_pattern(s,gp,p,p)=k;} //  no explicit migration for staying in this area, so create implicit
      }
      do_migr2=k;
      migr_start.initialize();
      // need to modify so it only does the calc for the first birthseason used for each GP
      for(gp=1;gp<=N_GP;gp++)
      {
        //  use firstBseas so that the start age of migration is calculated only for the first birthseason used for each GP
        firstBseas(gp)=0;
        for(g=1;g<=gmorph;g++)
        if(use_morph(g)>0)
        {
          if(GP4(g)==gp && firstBseas(gp)==0) firstBseas(gp)=Bseas(g);
        }
      }
      for(g=1;g<=gmorph;g++)
      if(use_morph(g)>0 && firstBseas(GP4(g))==Bseas(g))
      {
        for(s=1;s<=nseas;s++)
        {
        a=0;
        while(curr_age_beg(g,s,a)<migr_firstage) {a++;}
        migr_start(s,GP4(g))=a;
        }
      }
    }
 END_CALCS
   matrix move_def2(1,do_migr2,1,6)    //  movement definitions.  First Do_Migration of these are explicit; rest are implicit

 LOCAL_CALCS
    if(do_migration>0)
    {
      for(k=1;k<=do_migration;k++) {move_def2(k)=move_def(k);}
      k=do_migration;
      for(s=1;s<=nseas;s++)
      for(gp=1;gp<=N_GP;gp++)
      for(p=1;p<=pop;p++)
      {
        if(move_pattern(s,gp,p,p)>do_migration)
        {
          k++;
          move_def2(k,1)=s; move_def2(k,2)=gp; move_def2(k,3)=p; move_def2(k,4)=p; move_def2(k,5)=0; move_def2(k,6)=nages;
        }
      }
    }
 END_CALCS


//  ProgLabel_4.4  define the time blocks for time-varying parameters
  int k1
  int k2
  int k3
  init_int N_Block_Designs                      // read N block designs
  !!echoinput<<N_Block_Designs<<" N_Block_Designs"<<endl;
  init_ivector Nblk(1,N_Block_Designs)    // N blocks in each design
 LOCAL_CALCS
  if(N_Block_Designs>0) echoinput<<Nblk<<" N_Blocks_per design"<<endl;
  k1=N_Block_Designs;
  if(k1==0) k1=1;
 END_CALCS

  ivector Nblk2(1,k1)   //  vector to create ragged array of dimensions for block matrix
 LOCAL_CALCS
  Nblk2=2;
  if(N_Block_Designs>0) Nblk2=Nblk + Nblk;
 END_CALCS
  init_imatrix Block_Design(1,N_Block_Designs,1,Nblk2)  // read the ending year for each block

 LOCAL_CALCS
  if(N_Block_Designs>0)
  {
    echoinput<<" read block info "<<endl<<Block_Design<<endl;
    for(j=1;j<=N_Block_Designs;j++)
    {
      a=-1;
      for(k=1;k<=Nblk2(j)/2;k++)
      {
        a+=2;
        if(Block_Design(j,a+1)<Block_Design(j,a)) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Block:"<<j<<" "<<k<<" ends before it starts; fatal error"<<endl; exit(1);}
        if(Block_Design(j,a)<styr) {N_warn++; warning<<"Block:"<<j<<" "<<k<<" starts before styr; resetting"<<endl; Block_Design(j,a)=styr;}
        if(Block_Design(j,a+1)<styr) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Block:"<<j<<" "<<k<<" ends before styr; fatal error"<<endl; exit(1);}
        if(Block_Design(j,a)>retro_yr+1) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Block:"<<j<<" "<<k<<" starts after retroyr+1; fatal error"<<endl; exit(1);}
        if(Block_Design(j,a+1)>retro_yr+1) {N_warn++; warning<<"Block:"<<j<<" "<<k<<" ends after retroyr+1; resetting"<<endl; Block_Design(j,a+1)=retro_yr+1;}

      }
    }
  }
 END_CALCS
//  ProgLabel_4.5  read setup and parameters for natmort, growth, biology, recruitment distribution, and migration
// read setup for natmort parameters:  LO, HI, INIT, PRIOR, PR_type, CV, PHASE, use_env, use_dev, dev_minyr, dev_maxyr, dev_stddev, Block, Block_type
  int N_MGparm
  int N_natMparms
  int N_growparms
  int MGparm_per_def
  int recr_dist_parms
  number natM_amin;
  number natM_amax;
  init_number fracfemale;
  !!echoinput<<fracfemale<<" fracfemale"<<endl;
  !!if(fracfemale>=1.0) fracfemale=0.999999;
  !!if(fracfemale<=0.0) fracfemale=0.000001;

// read natmort setup
  init_int natM_type;  //  0=1Parm; 1=segmented; 2=Lorenzen; 3=agespecific; 4=agespec with seas interpolate
  !!echoinput<<natM_type<<" natM_type"<<endl;
  !! if(natM_type==1 || natM_type==2) {k=1;} else {k=0;}
  init_vector tempvec4(1,k)
 LOCAL_CALCS
  k=0; k1=0;
  if(natM_type==0)
  {N_natMparms=1;}
  else if(natM_type==1)
  {
    N_natMparms=tempvec4(1);  k=N_natMparms;
    echoinput<<N_natMparms<<" N_natMparms for segmented approach"<<endl;
  }
  else if(natM_type==2)
  {
    natM_amin=tempvec4(1);  N_natMparms=1;
    echoinput<<natM_amin<<" natM_A for Lorenzen"<<endl;
  }
  else
  {
    N_natMparms=0;
    if(natM_type>=3) {k1=N_GP*gender;}  // for reading age_natmort
  }
 END_CALCS

  init_vector NatM_break(1,k);  // these breakpoints only get read for natM_type=1
  !!if(k>0) echoinput<<NatM_break<<" NatM_breakages "<<endl;
  init_matrix Age_NatMort(1,k1,0,nages)
  !!if(k1>0) echoinput<<" Age_NatMort "<<Age_NatMort<<endl;

// read growth setup
  init_int Grow_type  // 1=vonbert; 2=Richards; 3=age-specific K;  4=read vector(not implemented)
  !!echoinput<<Grow_type<<" growth model "<<endl;
//  ProgLabel_4.5.1  create time constants for growth
  number AFIX;
  number AFIX2;
  number AFIX2_forCV;
  number AFIX_delta;
  number AFIX_plus;
  !! k=0;
  !! if(Grow_type<=2) {k=2;}  //  AFIX and AFIX2
  !! if (Grow_type==3) {k=3;}  //  add number of age-specific K changes
  init_vector tempvec5(1,k)
  int Age_K_count;
 LOCAL_CALCS
  Age_K_count=0;
  if(k>0) echoinput<<tempvec5<<" growth specifications"<<endl;
  k1=0;
  AFIX=0.;
  AFIX2=999.;  // this value invokes setting Linf equal to the L2 parameter
  if(Grow_type==1)
  {
    N_growparms=5;
    AFIX=tempvec5(1);
    AFIX2=tempvec5(2);
  }
  else if(Grow_type==2)
  {
    N_growparms=6;
    AFIX=tempvec5(1);
    AFIX2=tempvec5(2);
  }
  else if(Grow_type==3)
  {
    AFIX=tempvec5(1);
    AFIX2=tempvec5(2);
    Age_K_count=tempvec5(3);
    N_growparms=5+Age_K_count;;
  }
  else if(Grow_type==4)
  {
    N_growparms=2;  // for the two CV parameters
    k1=N_GP*gender;  // for reading age_natmort
  }
  
  AFIX2_forCV=AFIX2;
  if(AFIX2_forCV>nages) AFIX2_forCV=nages;

  AFIX_delta=AFIX2-AFIX;
  if(AFIX!=0.0)
  {AFIX_plus=AFIX;}
   else
   {AFIX_plus=1.0e-06;}
  MGparm_per_def=N_natMparms+N_growparms;
  lin_grow_beg.initialize();
  lin_grow_mid.initialize();

  for(g=1;g<=gmorph;g++)
  if(use_morph(g)>0)
  {
    first_grow_age=0;
    for(a=0;a<=nages;a++)
    {
      for(s=1;s<=nseas;s++)
      {
        if(a==0 && s<Bseas(g))
          {lin_grow_beg(g,s,a)=0.0;}
        else if(curr_age_beg(g,s,a)<AFIX)
          {lin_grow_beg(g,s,a)=curr_age_beg(g,s,a)/AFIX_plus;}
        else if(curr_age_beg(g,s,a)==AFIX)
          {lin_grow_beg(g,s,a)=1.0;}
        else if (first_grow_age==0)
          {lin_grow_beg(g,s,a)=-1.0; first_grow_age=1;}
        else
          {lin_grow_beg(g,s,a)=-2.0;}

        if(a==0 && s<Bseas(g))
          {lin_grow_mid(g,s,a)=0.0;}
        else if(curr_age_mid(g,s,a)<AFIX)
          {lin_grow_mid(g,s,a)=curr_age_mid(g,s,a)/AFIX_plus;}
        else if(curr_age_mid(g,s,a)==AFIX)
          {lin_grow_mid(g,s,a)=1.0;}
        else if (first_grow_age==0 || first_grow_age==1)
          {lin_grow_mid(g,s,a)=-1.0; first_grow_age=2;}
        else
          {lin_grow_mid(g,s,a)=-2.0;}
      }
    }
  }

 END_CALCS
  init_ivector Age_K_points(1,Age_K_count);  //  points at which age-specific multipliers to K will be applied
  !!if(Age_K_count>0) echoinput<<"Age-specific_K_points"<<Age_K_points<<endl;

  init_matrix Len_At_Age_rd(1,k1,0,nages)
  !!if(k1>0) echoinput<<"  Len_At_Age_rd"<<Len_At_Age_rd<<endl;

  init_number SD_add_to_LAA   // constant added to SD length-at-age (set to 0.1 for compatibility with SS2 V1.x
  !!echoinput<<SD_add_to_LAA<<"  SD_add_to_LAA"<<endl;
  init_int CV_depvar     //  select CV_growth pattern; 0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=f(A)   SS2 V1.x ony had CV=F(LAA)
  !!echoinput<<CV_depvar<<"  CV_depvar"<<endl;
  int CV_depvar_a;
  int CV_depvar_b;
  int Grow_logN
 LOCAL_CALCS

   Grow_logN=0;
   switch (CV_depvar)
   {
     case 0:
     {
       CV_depvar_a=0;
       CV_depvar_b=0;
       break;
     }
     case 1:
     {
       CV_depvar_a=1;
       CV_depvar_b=0;
       break;
     }
     case 2:
     {
       CV_depvar_a=0;
       CV_depvar_b=1;
       break;
     }
     case 3:
     {
       CV_depvar_a=1;
       CV_depvar_b=1;
       break;
     }
     case 4:
     {
       CV_depvar_a=1;
       CV_depvar_b=1;
       Grow_logN=1;
       break;
     }
   }
 END_CALCS

//  ProgLabel_4.5.2  process biology
   init_int Maturity_Option       // 1=length logistic; 2=age logistic; 3=read age-maturity
                                  // 4= read age-fecundity by growth_pattern 5=read all from separate wtatage.ss file
                                  //  6=read length-maturity
 int WTage_rd

 LOCAL_CALCS
  WTage_rd=0;
  echoinput<<Maturity_Option<<"  Maturity_Option"<<endl;
  if(Maturity_Option==3 || Maturity_Option==4)
    {k1=N_GP;}
  else
    {k1=0;}
  if(Maturity_Option==6)
    {k2=N_GP;}
  else
    {k2=0;}
    
  if(Maturity_Option==5)
  {
    echoinput<<" fecundity and weight at age to be read from file:  wtatage.ss"<<endl;
    WTage_rd=1;
  }
 END_CALCS
  init_matrix Age_Maturity(1,k1,0,nages)
  init_matrix Length_Maturity(1,k2,1,nlength)
  !!if(k1>0) echoinput<<"  read Age_Maturity for each GP"<<Age_Maturity<<endl;
  !!if(k2>0) echoinput<<"  read Length_Maturity for each GP"<<Length_Maturity<<endl;

  init_int First_Mature_Age     // first age with non-zero maturity
  !! echoinput<<First_Mature_Age<<"  First_Mature_Age"<<endl;

  init_int Fecund_Option
//   Value=1 means interpret the 2 egg parameters as linear eggs/kg on body weight (current SS default),
//   so eggs = wt * (a+b*wt), so value of a=1, b=0 causes eggs to be equiv to spawning biomass
//   Value=2 sets eggs=a*L^b   so cannot make equal to biomass
//   Value=3 sets eggs=a*W^b, so values of a=1, b=1 causes eggs to be equiv to spawning biomass
//   Value=4 sets eggs=a+b*L
//   Value=5 sets eggs=a+b*W
  !! echoinput<<Fecund_Option<<"  Fecundity option"<<endl;
  !! if(Fecund_Option>5) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Illegal fecundity option:  "<<Fecund_Option<<endl; exit(1);}
  init_int Hermaphro_Option
  int MGparm_Hermaphro  // pointer to start of these parameters
  !! echoinput<<Hermaphro_Option<<"  Hermaphro_Option "<<endl;
  !! MGparm_Hermaphro=0;
  !! k=0;
  !! if(Hermaphro_Option>0) k=2;
  init_ivector Hermaphro_more(1,k);
  int Hermaphro_seas;
  int Hermaphro_maleSPB;
 LOCAL_CALCS
  Hermaphro_seas=0;
  Hermaphro_maleSPB=0;
  if (k>0)
  {
    Hermaphro_seas=Hermaphro_more(1);
    Hermaphro_maleSPB=Hermaphro_more(2);
    echoinput<<Hermaphro_seas<<"  Hermaphro_season "<<endl;
    echoinput<<Hermaphro_maleSPB<<"  Hermaphro_maleSPB "<<endl;
  }
 END_CALCS
// if Hermaphro_Option=1, then read 3 parameters for switch from female to male by age
// FUTURE if Hermaphro_Option=2, then read 3 parameters for switch from female to male by age for each GrowPattern
// FUTURE if Hermaphro_Option=3, then read 3 parameters for switch from female to male by length
// FUTURE if Hermaphro_Option=4, then read 3 parameters for switch from female to male by length for each GrowPattern

   init_int MGparm_def       //  offset approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
   !! echoinput<<MGparm_def<<"  MGparm_def"<<endl;
   init_int MG_adjust_method   //  1=do V1.xx approach to adjustment by env, block or dev; 2=use new logistic approach
   !! echoinput<<MG_adjust_method<<"  MG_adjust_method"<<endl;

  imatrix time_vary_MG(styr-3,YrMax+1,0,6)  // goes to yrmax+1 to allow referencing in forecast, but only endyr+1 is checked
                                            // stores years to calc non-constant MG parms (1=natmort; 2=growth; 3=wtlen & fec; 4=recr_dist; 5=movement; 6=ageerrorkey)
  ivector MG_active(0,6)
  int do_once;
  int doit;
  vector femfrac(1,N_GP*gender)

  int MGP_CGD
  int CGD;  //  switch for cohort growth dev

 LOCAL_CALCS
  femfrac(1,N_GP)=fracfemale;
  if(gender==2) femfrac(N_GP+1,N_GP+N_GP)=1.-fracfemale;

  ParCount=0;

//  ProgLabel_4.5.3  set up indexing and parameter names for MG parameters
  N_MGparm=MGparm_per_def*N_GP*gender;     // dimension to hold the natmort, growth
  for(gg=1;gg<=gender;gg++)
  for(gp=1;gp<=N_GP;gp++)
  {
    for(k=1;k<=N_natMparms;k++)
    {
      ParCount++;
      onenum="    ";
//      _itoa( k, onenum, 10 );
      sprintf(onenum, "%d", k);
      ParmLabel+="NatM_p_"+onenum+"_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
    }
    switch (Grow_type)
    {
      case 1:
      {
        ParmLabel+="L_at_Amin_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParmLabel+="L_at_Amax_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParmLabel+="VonBert_K_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParCount+=3;
        break;
      }
      case 2:
      {
        ParmLabel+="L_at_Amin_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParmLabel+="L_at_Amax_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParmLabel+="VonBert_K_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParmLabel+="Richards_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParCount+=4;
        break;
      }
      case 3:
      {
        ParmLabel+="L_at_Amin_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParmLabel+="L_at_Amax_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParmLabel+="VonBert_K_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
        ParCount+=3;
        for(a=1;a<=Age_K_count;a++)
        {
          ParmLabel+="Age_K_"+GenderLbl(gg)+"_GP_"+NumLbl(gp)+"_a_"+NumLbl(Age_K_points(a));
          ParCount++;
        }
        break;
      }
    }
    ParmLabel+="CV_young_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
    ParmLabel+="CV_old_"+GenderLbl(gg)+"_GP_"+NumLbl(gp);
    ParCount+=2;
  }
  N_MGparm+= 2*gender+2+2; // add for wt-len(by gender), mat-len parms; eggs
  ParmLabel+="Wtlen_1_"+GenderLbl(1);
  ParmLabel+="Wtlen_2_"+GenderLbl(1);
  ParmLabel+="Mat50%_"+GenderLbl(1);
  ParmLabel+="Mat_slope_"+GenderLbl(1);
  if(Fecund_Option==1)
  {
    ParmLabel+="Eggs/kg_inter_"+GenderLbl(1);
    ParmLabel+="Eggs/kg_slope_wt_"+GenderLbl(1);
  }
  else if(Fecund_Option==2)
  {
     ParmLabel+="Eggs_scalar_"+GenderLbl(1);
     ParmLabel+="Eggs_exp_len_"+GenderLbl(1);
  }
  else if(Fecund_Option==3)
  {
     ParmLabel+="Eggs_scalar_"+GenderLbl(1);
     ParmLabel+="Eggs_exp_wt_"+GenderLbl(1);
  }
  else if(Fecund_Option==4)
  {
     ParmLabel+="Eggs_intercept_"+GenderLbl(1);
     ParmLabel+="Eggs_slope_len_"+GenderLbl(1);
  }
  else if(Fecund_Option==5)
  {
     ParmLabel+="Eggs_intercept_"+GenderLbl(1);
     ParmLabel+="Eggs_slope_Wt_"+GenderLbl(1);
  }
  ParCount+=6;

  if(gender==2)
  {
    ParmLabel+="Wtlen_1_"+GenderLbl(2);
    ParmLabel+="Wtlen_2_"+GenderLbl(2);
    ParCount+=2;
  }

  if(Hermaphro_Option==1)
  {
     MGparm_Hermaphro=N_MGparm+1;  // pointer to first parameter
     ParmLabel+="Herm_Infl_age";
     ParmLabel+="Herm_stdev";
     ParmLabel+="Herm_asymptote";
     ParCount+=3;
     N_MGparm+=3;
  }

  recr_dist_parms = N_MGparm;
  N_MGparm+=N_GP+pop+nseas;         // add for the assignment to areas
  for(k=1;k<=N_GP;k++) {ParCount++; ParmLabel+="RecrDist_GP_"+NumLbl(k);}
  for(k=1;k<=pop;k++)  {ParCount++; ParmLabel+="RecrDist_Area_"+NumLbl(k);}
  for(k=1;k<=nseas;k++){ParCount++; ParmLabel+="RecrDist_Seas_"+NumLbl(k);}

  if(recr_dist_inx==1)
    {
      N_MGparm+=N_GP*pop*nseas;  // add for the morph assignments within each area
      for(gp=1;gp<=N_GP;gp++)
      for(s=1;s<=nseas;s++)
      for(p=1;p<=pop;p++)
      {ParCount++; ParmLabel+="RecrDist_interaction_GP_"+NumLbl(gp)+"_seas_"+NumLbl(s)+"_area_"+NumLbl(p);}
    }
  N_MGparm++;
  MGP_CGD=N_MGparm;             // add 1 parameter for cohort-specific growth parameter
  ParCount++;
  ParmLabel+="CohortGrowDev";

  if(do_migration>0)
  {
    N_MGparm+=2*do_migration;
   for(k=1;k<=do_migration;k++)
     {
     s=move_def(k,1); gp=move_def(k,2); p=move_def(k,3); p2=move_def(k,4);
     ParCount++; ParmLabel+="MoveParm_A_seas_"+NumLbl(s)+"_GP_"+NumLbl(gp)+"from_"+NumLbl(p)+"to_"+NumLbl(p2);
     ParCount++; ParmLabel+="MoveParm_B_seas_"+NumLbl(s)+"_GP_"+NumLbl(gp)+"from_"+NumLbl(p)+"to_"+NumLbl(p2);
    }
  }

  if(Use_AgeKeyZero>0)
  {
    AgeKeyParm=N_MGparm+1;
    for(k=1;k<=7;k++)
    {
       N_MGparm++; ParCount++; ParmLabel+="AgeKeyParm"+NumLbl(k);
    }
  }
 END_CALCS

  init_matrix MGparm_1(1,N_MGparm,1,14)   // matrix with natmort and growth parms controls
  ivector MGparm_offset(1,N_MGparm)

 LOCAL_CALCS
  echoinput<<" Biology parameter setup"<<endl;
  for(i=1;i<=N_MGparm;i++)
  echoinput<<i<<" "<<MGparm_1(i)<<" "<<ParmLabel(ParCount-N_MGparm+i)<<endl;

//  find MGparms for which the male parameter value is set equal to the female value
//  only applies for MGparm_def==1 which is direct estimation (no offsets)
//  only for the natmort and growth parameters (not wtlen, fecundity, movement, recr distribution)
  MGparm_offset.initialize();
  if(MGparm_def==1)
  {
    Ip=0;
    for (gp=1;gp<=N_GP*gender;gp++)
    {
      if(gp>1)
      {
        for(j=1;j<=MGparm_per_def;j++)
        {
          if(MGparm_1(Ip+j,3)==0.0 && MGparm_1(Ip+j,7)<0) MGparm_offset(Ip+j)=Ip+j-N_GP*MGparm_per_def;  // value is zero and not estimated
        }
      }
      Ip+=MGparm_per_def;
    }
  }
  echoinput<<"Now read blocks and other adjustments to MGparms "<<endl;
 END_CALCS

//  ProgLabel_4.5.4  set up environmental linkage for MG parms
  int N_MGparm_env                            //  number of MGparms that use env linkage
  int customMGenvsetup  //  0=read one setup (if necessary) and apply to all; 1=read each
  ivector MGparm_env(1,N_MGparm)   // contains the parameter number of the envlink for a
  ivector MGparm_envuse(1,N_MGparm)   // contains the environment data number
  ivector MGparm_envtype(1,N_MGparm)  // 1=multiplicative; 2= additive
  ivector mgp_type(1,N_MGparm)  //  contains category to parameter (1=natmort; 2=growth; 3=wtlen & fec; 4=recr_dist; 5=movement)

 LOCAL_CALCS
   time_vary_MG.initialize();    // stores years to calc non-constant MG parms (1=natmort; 2=growth; 3=wtlen & fec; 4=recr_dist; 5=movement)
   CGD=0;
   k=1;
   for(f=1;f<=N_GP*gender;f++)
   {
     mgp_type(k,k+N_natMparms-1)=1; // natmort parms
     k+=N_natMparms;
//     mgp_type(k,k+4)=2;  // growth parms
//     k=k+5;
     mgp_type(k,k+N_growparms)=2;  // growth parms
     k=k+N_growparms;
   }
   mgp_type(k,k+5)=3;  k=k+6;   // female wtlen and fecundity
   if(gender==2) {mgp_type(k,k+1)=3; k=k+2;}  // male wt-len
   if(Hermaphro_Option>0) {mgp_type(k,k+2)=3; k=k+3;}  //   herma parameters done with wtlen and fecundity
   mgp_type(k,MGP_CGD-1)=4;   // recruit apportionments
   mgp_type(MGP_CGD)=2;   // cohort growth dev
   if(do_migration>0)  mgp_type(MGP_CGD+1,N_MGparm)=5;
   if(Use_AgeKeyZero>0) mgp_type(AgeKeyParm,N_MGparm)=6;

   MGparm_env.initialize();   //  will store the index of environ fxns here
   MGparm_envtype.initialize();
   N_MGparm_env=0;
   for (f=1;f<=N_MGparm;f++)
   {
    if(MGparm_1(f,8)!=0)
    {
     N_MGparm_env ++;  MGparm_env(f)=N_MGparm+N_MGparm_env;
     if(MGparm_1(f,8)>0)
     {
       ParCount++; ParmLabel+=ParmLabel(f)+"_ENV_mult"; MGparm_envtype(f)=1; MGparm_envuse(f)=MGparm_1(f,8);
       if(MG_adjust_method==2) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"multiplicative env effect on MGparm: "<<f
        <<" not allowed because MG_adjust_method==2; STOP"<<endl; exit(1);}
     }
     else if(MGparm_1(f,8)==-999)
     {ParCount++; ParmLabel+=ParmLabel(f)+"_ENV_densdep"; MGparm_envtype(f)=3;  MGparm_envuse(f)=-1;}
     else
     {ParCount++; ParmLabel+=ParmLabel(f)+"_ENV_add"; MGparm_envtype(f)=2; MGparm_envuse(f)=-MGparm_1(f,8);}

     if(f==MGP_CGD) CGD=1;    // cohort growth dev is a fxn of environ, so turn on CGD calculation
     for(y=styr;y<=endyr;y++)
     {
      if(env_data_RD(y,MGparm_envuse(f))!=0.0 || MGparm_envtype(f)==3) {time_vary_MG(y,mgp_type(f))=1; time_vary_MG(y+1,mgp_type(f))=1; }
      //       non-zero data were read    or fxn uses biomass or recruitment
     }
    }
   }

  if(N_MGparm_env>0)
  {
    *(ad_comm::global_datafile) >> customMGenvsetup;
    if(customMGenvsetup==0) {k1=1;} else {k1=N_MGparm_env;}
   echoinput<<customMGenvsetup<<" customMGenvsetup"<<endl;
  }
  else
  {customMGenvsetup=0; k1=0;
   echoinput<<" no mgparm env links, so don't read customMGenvsetup"<<endl;
    }
 END_CALCS
  init_matrix MGparm_env_1(1,k1,1,7)
  !!if(N_MGparm_env>0) echoinput<<" MGparm-env setup "<<endl<<MGparm_env_1<<endl;


//  ProgLabel_4.5.5  set up block for MG parms
  int N_MGparm_blk                            // number of MGparms that use blocks
  imatrix Block_Defs_MG(1,N_MGparm,styr,endyr+1)
  int N_MGparm_trend     //   number of MG parameters using trend or cycle
  int N_MGparm_trend2     //   number of parameters needed to define trends and cycles
  ivector MGparm_trend_point(1,N_MGparm)   //  index of trend parameters associated with each MG parm


 LOCAL_CALCS
  echoinput<<"Process and create labels for the MGparm adjustments"<<endl;
  Block_Defs_MG.initialize();
  N_MGparm_blk=0;  // counter for assigned parms
  for(j=1;j<=N_MGparm;j++)
  {
   z=MGparm_1(j,13);    // specified block definition
   if(z>N_Block_Designs) {N_warn++; warning<<" ERROR, Block > N Blocks "<<z<<" "<<N_Block_Designs<<endl;}
   if(z>0)
   {
     g=1;
     for(a=1;a<=Nblk(z);a++)
     {
      N_MGparm_blk++;
      y=Block_Design(z,g);
      time_vary_MG(y,mgp_type(j))=1;
      sprintf(onenum, "%d", y);
      ParCount++;
      k=int(MGparm_1(j,14));
      switch(k)
      {
        case 0:
        {ParmLabel+=ParmLabel(j)+"_BLK"+NumLbl(z)+"mult_"+onenum+CRLF(1);  break;}
        case 1:
        {ParmLabel+=ParmLabel(j)+"_BLK"+NumLbl(z)+"add_"+onenum+CRLF(1);  break;}
        case 2:
        {ParmLabel+=ParmLabel(j)+"_BLK"+NumLbl(z)+"repl_"+onenum+CRLF(1);  break;}
        case 3:
        {ParmLabel+=ParmLabel(j)+"_BLK"+NumLbl(z)+"delta_"+onenum+CRLF(1);  break;}
      }

      y=Block_Design(z,g+1)+1;  // first year after block
      if(y>endyr+1) y=endyr+1;
      time_vary_MG(y,mgp_type(j))=1;
     for(y=Block_Design(z,g);y<=Block_Design(z,g+1);y++)  // loop years for this block, including yrs past endyr+1
     {
      Block_Defs_MG(j,y)=N_MGparm+N_MGparm_env+N_MGparm_blk;
     }
     g+=2;
    }
    echoinput<<"Block definitions for MGparms"<<endl<<Block_Defs_MG(j)<<endl;
    if(j==MGP_CGD) CGD=1;
   }
   else if(z<0)  // will do parameter trend approach
   {
      for(y=styr;y<=YrMax+1;y++)
      {
        time_vary_MG(y,mgp_type(j))=1;
      }
   }
  }
 END_CALCS

  int customblocksetup_MG  //  0=read one setup and apply to all; 1=read each
 LOCAL_CALCS
  if(N_MGparm_blk>0)
  {
    *(ad_comm::global_datafile) >> customblocksetup_MG;
    if(customblocksetup_MG==0)
    {k1=1;}
    else
    {k1=N_MGparm_blk;}
    echoinput<<customblocksetup_MG<<" customblocksetup_MG"<<endl;
  }
  else
  {
    customblocksetup_MG=0;
    k1=0;
    echoinput<<" no mgparm blocks, so don't read customblocksetup_MG"<<endl;
  }
 END_CALCS
  init_matrix MGparm_blk_1(1,k1,1,7)  // read matrix that defines the block parms
  !!if(N_MGparm_blk>0) echoinput<<" MGparm-blk setup "<<endl<<MGparm_blk_1<<endl;

//  ProgLabel_4.5.6  Setup trends and cycles as alternative to blocks
 LOCAL_CALCS
  if(N_MGparm_blk>0) echoinput<<" MGparm-blk setup "<<endl<<MGparm_blk_1<<endl;
// use negative block as indicator to use time trend
  N_MGparm_trend=0;
  N_MGparm_trend2=0;
  for(j=1;j<=N_MGparm;j++)
  {
    if(MGparm_1(j,13)<0)  //  create timetrend parameter
    {
      N_MGparm_trend++;
      MGparm_trend_point(j)=N_MGparm_trend;
      if(MGparm_1(j,13)==-1)
      {
        ParCount++; ParmLabel+=ParmLabel(j)+"_TrendFinal_Offset"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j)+"_TrendInfl_"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j)+"_TrendWidth_"+CRLF(1);
        N_MGparm_trend2+=3;
      }
      else if(MGparm_1(j,13)==-2)
      {
        ParCount++; ParmLabel+=ParmLabel(j)+"_TrendFinal_"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j)+"_TrendInfl_"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j)+"_TrendWidth_"+CRLF(1);
        N_MGparm_trend2+=3;
      }
      else
      {
        for(icycle=1;icycle<=Ncycle;icycle++)
        {
          ParCount++; ParmLabel+=ParmLabel(j)+"_Cycle_"+NumLbl(icycle)+CRLF(1);
          N_MGparm_trend2++;
        }
      }
    }
  }
 END_CALCS

  init_matrix MGparm_trend_1(1,N_MGparm_trend2,1,7)  // read matrix that defines the parms and trend parms
  !!if(N_MGparm_trend2>0) echoinput<<"MG trend and cycle parameters "<<endl<<MGparm_trend_1<<endl;

  ivector MGparm_trend_rev(1,N_MGparm_trend)
  ivector MGparm_trend_rev_1(1,N_MGparm_trend)
 LOCAL_CALCS
  if(N_MGparm_trend>0)
  {
    k1=0;
    k2=N_MGparm+N_MGparm_env+N_MGparm_blk;
    for(j=1;j<=N_MGparm;j++)
    {
      if(MGparm_1(j,13)<0)  //  timetrend exists
      {
        k1++;
        MGparm_trend_rev(k1)=j;  // reverse pointer from trend to affected parameter
        MGparm_trend_rev_1(k1)=k2;  // pointer to base in list of MGparms (so k2+1 is first parameter used)
        if(MGparm_1(j,13)>=-2)  //  timetrend
        {k2+=3;}
        else
        {k2+=Ncycle;}
      }
    }
  }
 END_CALCS

//  ProgLabel_4.5.7  set up seasonal effects for MG parms
  init_ivector MGparm_seas_effects(1,10)  // femwtlen1, femwtlen2, mat1, mat2, fec1 fec2 Malewtlen1, malewtlen2 L1 K
  int MGparm_doseas
  int N_MGparm_seas                            // number of MGparms that use seasonal effects
 LOCAL_CALCS
   echoinput<<MGparm_seas_effects<<" MGparm_seas_effects"<<endl;
  adstring_array MGseasLbl;
  MGseasLbl+="F-WL1"+CRLF(1);
  MGseasLbl+="F-WL2"+CRLF(1);
  MGseasLbl+="F-Mat1"+CRLF(1);
  MGseasLbl+="F-Mat1"+CRLF(1);
  MGseasLbl+="F-Fec1"+CRLF(1);
  MGseasLbl+="F-Fec1"+CRLF(1);
  MGseasLbl+="M-WL1"+CRLF(1);
  MGseasLbl+="M-WL2"+CRLF(1);
  MGseasLbl+="L1"+CRLF(1);
  MGseasLbl+="VBK"+CRLF(1);
  MGparm_doseas=sum(MGparm_seas_effects);
  N_MGparm_seas=0;  // counter for assigned parms
  if(MGparm_doseas>0)
  {
    for(j=1;j<=10;j++)
    {
      if(MGparm_seas_effects(j)>0)
      {
        MGparm_seas_effects(j)=N_MGparm+N_MGparm_env+N_MGparm_blk+N_MGparm_seas;  // store base parameter count
        for(s=1;s<=nseas;s++)
        {
          N_MGparm_seas++; ParCount++; ParmLabel+=MGseasLbl(j)+"_seas_"+NumLbl(s);
        }
      }
    }
  }
 END_CALCS
  init_matrix MGparm_seas_1(1,N_MGparm_seas,1,7)  // read matrix that defines the seasonal parms
  !!if(N_MGparm_seas>0) echoinput<<" MGparm_seas"<<endl<<MGparm_seas_1<<endl;

//  ProgLabel_4.5.7  Create vectors (e.g. MGparm_PH) to be used to define the actual estimated parameter array
  int N_MGparm2
  !!N_MGparm2=N_MGparm+N_MGparm_env+N_MGparm_blk+N_MGparm_trend2+N_MGparm_seas;
  vector MGparm_LO(1,N_MGparm2)
  vector MGparm_HI(1,N_MGparm2)
  vector MGparm_RD(1,N_MGparm2)
  vector MGparm_PR(1,N_MGparm2)
  ivector MGparm_PRtype(1,N_MGparm2)
  vector MGparm_CV(1,N_MGparm2)
  ivector MGparm_PH(1,N_MGparm2)

 LOCAL_CALCS
   MG_active=0;   // initializes
   for (f=1;f<=N_MGparm;f++)
   {
    MGparm_LO(f)=MGparm_1(f,1);
    MGparm_HI(f)=MGparm_1(f,2);
    MGparm_RD(f)=MGparm_1(f,3);
    MGparm_PR(f)=MGparm_1(f,4);
    MGparm_PRtype(f)=MGparm_1(f,5);
    MGparm_CV(f)=MGparm_1(f,6);
    MGparm_PH(f)=MGparm_1(f,7);
    if(MGparm_PH(f)>0)
    {MG_active(mgp_type(f))=1;}
   }
   if(natM_type==2 && MG_active(2)>0) MG_active(1)=1;  // lorenzen M depends on growth

   j=N_MGparm;
   if(N_MGparm_env>0)
   {
    for (f=1;f<=N_MGparm_env;f++)
    {
     j++;
     if(customMGenvsetup==0) {k=1;}
     else {k=f;}

    MGparm_LO(j)=MGparm_env_1(k,1);
     MGparm_HI(j)=MGparm_env_1(k,2);
     MGparm_RD(j)=MGparm_env_1(k,3);
     MGparm_PR(j)=MGparm_env_1(k,4);
     MGparm_PRtype(j)=MGparm_env_1(k,5);
     MGparm_CV(j)=MGparm_env_1(k,6);
     MGparm_PH(j)=MGparm_env_1(k,7);
    }
   }

   if(N_MGparm_blk>0)
   for (f=1;f<=N_MGparm_blk;f++)
   {
    j++;
    if(customblocksetup_MG==0) k=1;
    else k=f;
    MGparm_LO(j)=MGparm_blk_1(k,1);
    MGparm_HI(j)=MGparm_blk_1(k,2);
    MGparm_RD(j)=MGparm_blk_1(k,3);
    MGparm_PR(j)=MGparm_blk_1(k,4);
    MGparm_PRtype(j)=MGparm_blk_1(k,5);
    MGparm_CV(j)=MGparm_blk_1(k,6);
    MGparm_PH(j)=MGparm_blk_1(k,7);
   }
   if(N_MGparm_trend>0)
   for (f=1;f<=N_MGparm_trend2;f++)
   {
    j++;
    MGparm_LO(j)=MGparm_trend_1(f,1);
    MGparm_HI(j)=MGparm_trend_1(f,2);
    MGparm_RD(j)=MGparm_trend_1(f,3);
    MGparm_PR(j)=MGparm_trend_1(f,4);
    MGparm_PRtype(j)=MGparm_trend_1(f,5);
    MGparm_CV(j)=MGparm_trend_1(f,6);
    MGparm_PH(j)=MGparm_trend_1(f,7);
   }

   if(N_MGparm_seas>0)
   for (f=1;f<=N_MGparm_seas;f++)
   {
    j++;
    MGparm_LO(j)=MGparm_seas_1(f,1);
    MGparm_HI(j)=MGparm_seas_1(f,2);
    MGparm_RD(j)=MGparm_seas_1(f,3);
    MGparm_PR(j)=MGparm_seas_1(f,4);
    MGparm_PRtype(j)=MGparm_seas_1(f,5);
    MGparm_CV(j)=MGparm_seas_1(f,6);
    MGparm_PH(j)=MGparm_seas_1(f,7);
   }
 END_CALCS

//  ProgLabel_4.5.8  set up random deviations for MG parms
  int N_MGparm_dev                            //  number of MGparms that use annual deviations
 LOCAL_CALCS
    N_MGparm_dev=0;
    for(f=1;f<=N_MGparm;f++)
    {
    if(MGparm_1(f,9)>=1) N_MGparm_dev ++;
    }
 END_CALCS

  ivector MGparm_dev_minyr(1,N_MGparm_dev)
  ivector MGparm_dev_maxyr(1,N_MGparm_dev)
  vector  MGparm_dev_stddev(1,N_MGparm_dev)
  ivector MGparm_dev_type(1,N_MGparm_dev)  // contains type of dev:  1 for multiplicative, 2 for additive, 3 for additive randwalk
  ivector MGparm_dev_select(1,N_MGparm)
  int MGparm_dev_PH
 LOCAL_CALCS
   MGparm_dev_select.initialize();
   int N_MGparm_dev_tot=0;
   if(N_MGparm_dev>0)
     {
       j=0;
       for (f=1;f<=N_MGparm;f++)
       {
       if(MGparm_1(f,9)>=1)
       {
           j++;
           if(MG_adjust_method==2 && MGparm_1(f,9)==1)
           {N_warn++; warning<<" cannot use MG_adjust_method==2 and multiplicative devs for parameter "<<f<<endl;}
           MGparm_dev_type(j)=MGparm_1(f,9);  //  1 for multiplicative, 2 for additive, 3 for additive randwalk
           MGparm_dev_select(f)=j;
           MGparm_dev_stddev(j)=MGparm_1(f,12);
           y=MGparm_1(f,10);
           if(y<styr)
           {
             N_warn++; warning<<" reset MGparm_dev start year to styr for MGparm: "<<f<<" "<<y<<endl;
             y=styr;
           }
           MGparm_dev_minyr(j)=y;

           y=MGparm_1(f,11);
           if(y>endyr)
           {
             N_warn++; warning<<" reset MGparm_dev end year to endyr for MGparm: "<<f<<" "<<y<<endl;
             y=endyr;
           }
           MGparm_dev_maxyr(j)=y;
           for(y=MGparm_dev_minyr(j);y<=MGparm_dev_maxyr(j);y++)
           {
             MG_active(mgp_type(f))=1;
             time_vary_MG(y,mgp_type(f))=1;
             if(y<=endyr) time_vary_MG(y+1,mgp_type(f))=1;   // so will recalculate to null value, even for endyr+1
             sprintf(onenum, "%d", y);
             N_MGparm_dev_tot++;
             ParCount++;
             if(MGparm_dev_type(j)==1)
             {ParmLabel+=ParmLabel(f)+"_DEVmult_"+onenum+CRLF(1);}
             else if(MGparm_dev_type(j)==2)
             {ParmLabel+=ParmLabel(f)+"_DEVadd_"+onenum+CRLF(1);}
             else if(MGparm_dev_type(j)==3)
             {ParmLabel+=ParmLabel(f)+"_DEVrwalk_"+onenum+CRLF(1);}
           else
           {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" illegal MGparmdevtype for parm "<<f<<endl; exit(1);}
         }
         if(f==MGP_CGD) CGD=1;
         }
       }
       *(ad_comm::global_datafile) >> MGparm_dev_PH;
       echoinput<<MGparm_dev_PH<<" MGparm_dev_PH"<<endl;
     }
     else
     {
      MGparm_dev_PH=-6;
      echoinput<<" don't read MGparm_dev_PH"<<endl;
    }

//  ProgLabel_4.5.9  populate time/bio_category array defining when biology changes
     k=YrMax+1;
    for(y=styr+1;y<=YrMax;y++)
    {
      if(time_vary_MG(y,2)>0 && y<k)  k=y;
    }
    if(k<YrMax+1)
    {
      for(y=k;y<=YrMax+1;y++)
      {
        time_vary_MG(y,2)=1;
      }
    }
    for(y=styr;y<=YrMax;y++)
    {
      for(f=1;f<=6;f++)
      {
        if(time_vary_MG(y,f)>0)
        {
          MG_active(f)=1;
          time_vary_MG(y,0)=1;  // tracks active status for all MG types
        }
      }
    }
    MG_active(0)=sum(MG_active(1,6));
 END_CALCS

//**************************************
//  ProgLabel_4.6  read setup for Spawner-Recruitment parameters
// read setup for SR parameters:  LO, HI, INIT, PRIOR, PRtype, CV, PHASE
  init_int SR_fxn
  ivector N_SRparm(1,10)
  !!N_SRparm.fill("{0,2,2,2,3,2,3,3,0,0}");
  int N_SRparm2
  !!echoinput<<SR_fxn<<" #_SR_function: 1=null; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=Survival_3Parm; 8=Shepard "<<endl;
  !!N_SRparm2=N_SRparm(SR_fxn)+4;
  init_matrix SR_parm_1(1,N_SRparm2,1,7)
  !!echoinput<<" SR parms "<<endl<<SR_parm_1<<endl;
  init_int SR_env_link
  !!echoinput<<SR_env_link<<" SR_env_link "<<endl;
  init_int SR_env_target_RD   // 0=none; 1=devs; 2=R0; 3=steepness
  !!echoinput<<SR_env_target_RD<<" SR_env_target_RD "<<endl;
  int SR_env_target
  int SR_autocorr;  // will be calculated later

  vector SRvec_LO(1,N_SRparm2)
  vector SRvec_HI(1,N_SRparm2)
  ivector SRvec_PH(1,N_SRparm2)

//  ProgLabel_4.6.1  create S-R parameter labels
 LOCAL_CALCS
   SRvec_LO=column(SR_parm_1,1);
   SRvec_HI=column(SR_parm_1,2);
   SRvec_PH=ivector(column(SR_parm_1,7));
   if(SR_env_link>N_envvar)
   {
     N_warn++;
     warning<<" ERROR:  SR_env_link ( "<<SR_env_link<<" ) was set greater than the highest numbered environmental index ( "<<N_envvar<<" )"<<endl;
     cout<<" EXIT - see warning "<<endl; exit(1);
   }
   SR_env_target=SR_env_target_RD;
   if(SR_env_link==0) SR_env_target=0;
   if(SR_env_link==0 && SR_env_target_RD>0)
   {N_warn++; warning<<" WARNING:  SR_env_target was set, but no SR_env_link selected, SR_env_target set to 0"<<endl;}
//#_SR_function: 1=null; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=Survival_3Parm "<<endl;
  ParmLabel+="SR_LN(R0)";
  switch(SR_fxn)
  {
    case 1: // previous placement for B-H constrained
    {
      N_warn++; cout<<"Critical error:  see warning"<<endl; warning<<"B-H constrained curve is now Spawn-Recr option #6"<<endl; exit(1);
      break;
    }
    case 2:  // Ricker
    {
      ParmLabel+="SR_Ricker";
      break;
    }
    case 3:  // Bev-Holt
    {
      ParmLabel+="SR_BH_steep";
      break;
    }
    case 4:  // SCAA
    {
      ParmLabel+="SR_SCAA_null";
      break;
    }
    case 5:  // Hockey
    {
      ParmLabel+="SR_hockey_infl";
      ParmLabel+="SR_hockey_min_R";
      break;
    }
    case 6:  // Bev-Holt flattop
    {
      ParmLabel+="SR_BH_flat_steep";
      break;
    }
    case 7:  // survival
    {
      ParmLabel+="SR_surv_Sfrac";
      ParmLabel+="SR_surv_Beta";
      break;
    }
    case 8:  // shepard
    {
      ParmLabel+="SR_steepness";
      ParmLabel+="SR_Shepard_c";
      break;
    }
   }
  ParmLabel+="SR_sigmaR";
  ParmLabel+="SR_envlink";
  ParmLabel+="SR_R1_offset";
  ParmLabel+="SR_autocorr";
  ParCount+=N_SRparm2;
 END_CALCS

  init_int do_recdev  //  0=none; 1=devvector; 2=simple deviations
  !!echoinput<<do_recdev<<" do_recdev"<<endl;
  init_int recdev_start;
  !!echoinput<<recdev_start<<" recdev_start"<<endl;
  init_int recdev_end;
  !!echoinput<<recdev_end<<" recdev_end"<<endl;
  init_int recdev_PH_rd;
  !!echoinput<<recdev_PH_rd<<" recdev_PH"<<endl;
  int recdev_PH;
  !! recdev_PH=recdev_PH_rd;
  init_int recdev_adv
  !!echoinput<<recdev_adv<<" recdev_adv"<<endl;

  init_vector recdev_options_rd(1,13*recdev_adv)
  vector recdev_options(1,13)
  int recdev_early_start_rd
  int recdev_early_start
  int recdev_early_end
  int recdev_first
  int recdev_early_PH
  int Fcast_recr_PH
  int Fcast_recr_PH2
  number Fcast_recr_lambda
  vector recdev_adj(1,5)
  int recdev_cycle
  int recdev_do_early
  int recdev_read
  number recdev_LO;
  number recdev_HI;
  ivector recdev_doit(styr-nages,endyr+1)

//  ProgLabel_4.6.2  setup advanced recruitment options
 LOCAL_CALCS
  recdev_doit=0;
  if(recdev_adv>0)
  {
    recdev_options(1,13)=recdev_options_rd(1,13);
    recdev_early_start_rd=recdev_options(1);
    recdev_early_PH=recdev_options(2);
    Fcast_recr_PH=recdev_options(3);
    Fcast_recr_lambda=recdev_options(4);
    recdev_adj(1)=recdev_options(5);
    recdev_adj(2)=recdev_options(6);
    recdev_adj(3)=recdev_options(7);
    recdev_adj(4)=recdev_options(8);
    recdev_adj(5)=recdev_options(9);  // maxbias adj

    recdev_cycle=recdev_options(10);
    recdev_LO=recdev_options(11);
    recdev_HI=recdev_options(12);
    recdev_read=recdev_options(13);
  }
  else
  {
    recdev_early_start_rd=0;   // 0 means no early
    recdev_early_end=-1;
    recdev_early_PH=-4;
    recdev_options(2)=recdev_early_PH;
    Fcast_recr_PH=0;  // so will be reset to maxphase+1
    recdev_options(3)=Fcast_recr_PH;
    Fcast_recr_lambda=1.;
    recdev_adj(1)=double(styr)-1000.;
    recdev_adj(2)=styr-nages;
    recdev_adj(3)=recdev_end;
    recdev_adj(4)=double(endyr)+1.;
    recdev_adj(5)=1.0;
    recdev_cycle=0;
    recdev_LO=-5;
    recdev_HI=5;
    recdev_read=0;
  }

  recdev_early_start=recdev_early_start_rd;
  if(recdev_adv>0)
  {echoinput<<"#_start of advanced SR options"<<endl;}
  else
  {echoinput<<"# advanced options not read;  defaults displayed below"<<endl;}

    echoinput<<recdev_early_start_rd<<" #_recdev_early_start (0=none; neg value makes relative to recdev_start)"<<endl;
    echoinput<<recdev_early_PH<<" #_recdev_early_phase"<<endl;
    echoinput<<Fcast_recr_PH<<" #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)"<<endl;
    echoinput<<Fcast_recr_lambda<<" #_lambda for Fcast_recr_like occurring before endyr+1"<<endl;
    echoinput<<recdev_adj(1)<<" #_last_early_yr_nobias_adj_in_MPD"<<endl;
    echoinput<<recdev_adj(2)<<" #_first_yr_fullbias_adj_in_MPD"<<endl;
    echoinput<<recdev_adj(3)<<" #_last_yr_fullbias_adj_in_MPD"<<endl;
    echoinput<<recdev_adj(4)<<" #_first_recent_yr_nobias_adj_in_MPD"<<endl;
    echoinput<<recdev_adj(5)<<" #_max_bias_adj_in_MPD"<<endl;
    echoinput<<recdev_cycle<<" # period of cycle in recruitment "<<endl;
    echoinput<<recdev_LO<<" #min rec_dev"<<endl;
    echoinput<<recdev_HI<<" #max rec_dev"<<endl;
    echoinput<<recdev_read<<" #_read_recdevs"<<endl;
    echoinput<<"#_end of advanced SR options"<<endl;

 END_CALCS

 LOCAL_CALCS
//  ProgLabel_4.6.3  create parm labels for recruitment cycle parameters
  if(recdev_cycle>0)
  {
    for(y=1;y<=recdev_cycle;y++)
    {
      ParCount++;
      sprintf(onenum, "%d", y);
      ParmLabel+="RecrDev_Cycle_"+onenum+CRLF(1);
    }
  }

//  ProgLabel_4.6.4  setup recruitment deviations and create parm labels for each year
  if(recdev_end>retro_yr) recdev_end=retro_yr;
  if(recdev_start<(styr-nages)) {recdev_start=styr-nages; N_warn++; warning<<" adjusting recdev_start to: "<<recdev_start<<endl;}
  recdev_first=recdev_start;   // stores first recdev, whether from the early period or the standard dev period

  if(recdev_early_start>=recdev_start)
  {
    N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" error, cannot set recdev_early_start after main recdev start"<<endl;
    exit(1);
  }
  else if(recdev_early_start==0)  // do not do early rec devs
  {
    recdev_do_early=0;
    recdev_early_end=-1;
    if(recdev_early_PH>0) recdev_early_PH=-recdev_early_PH;
  }
  else
  {
    if(recdev_early_start<0) recdev_early_start+=recdev_start;  // do relative to start of recdevs
    recdev_do_early=1;
    if(recdev_early_start<(styr-nages))
      {recdev_early_start=styr-nages; N_warn++; warning<<" adjusting recdev_early to: "<<recdev_early_start<<endl;}
    if(recdev_start-recdev_early_start<6)
    {N_warn++; warning<<" Are you sure you want so few early recrdevs? "<<recdev_start-recdev_early_start<<endl;}

    recdev_first=recdev_early_start;  // because this is before recdev_start
    recdev_early_end=recdev_start-1;
    for(y=recdev_early_start;y<=recdev_early_end;y++)
    {
      ParCount++;
      recdev_doit(y)=1;
      if(y>=styr)
      {
        sprintf(onenum, "%d", y);
        ParmLabel+="Early_RecrDev_"+onenum+CRLF(1);
      }
      else
      {
        onenum="    ";
        sprintf(onenum, "%d", styr-y);
        ParmLabel+="Early_InitAge_"+onenum+CRLF(1);
      }
    }
  }

  if(do_recdev>0)
  {
    for(y=recdev_start;y<=recdev_end;y++)
    {
      ParCount++;
      recdev_doit(y)=1;

        if(y>=styr)
        {
        sprintf(onenum, "%d", y);
        ParmLabel+="Main_RecrDev_"+onenum+CRLF(1);
      }
      else
        {
          onenum="    ";
        sprintf(onenum, "%d", styr-y);
        ParmLabel+="Main_InitAge_"+onenum+CRLF(1);
      }
    }
  }

  if(Do_Forecast>0)
  {
    for(y=recdev_end+1;y<=YrMax;y++)
    {
      sprintf(onenum, "%d", y);
      ParCount++;
      if(y>endyr)
      {ParmLabel+="ForeRecr_"+onenum+CRLF(1);}
      else
      {ParmLabel+="Late_RecrDev_"+onenum+CRLF(1);}
    }

    for(y=endyr+1;y<=YrMax;y++)
    {
      sprintf(onenum, "%d", y);
      ParCount++;
      ParmLabel+="Impl_err_"+onenum+CRLF(1);
    }
  }
 END_CALCS

//  ProgLabel_4.6.5  read recdev_cycle parameters and input recruitment deviations if needed
  init_matrix recdev_cycle_parm_RD(1,recdev_cycle,1,14);
  !!k=1;
  !!if(recdev_cycle>0) k=recdev_cycle;
  vector recdev_cycle_LO(1,k);
  vector recdev_cycle_HI(1,k);
  ivector recdev_cycle_PH(1,k);
  !!if(recdev_cycle>0) echoinput<<"recruitment cycle input "<<endl<<recdev_cycle_parm_RD<<endl;

  init_matrix recdev_input(1,recdev_read,1,2);
  !!if(recdev_read>0) echoinput<<"recruitment deviation input "<<endl<<recdev_input<<endl;

//  ProgLabel_4.7  Input F_method setup
  init_number F_ballpark
  !! echoinput<<F_ballpark<<" F ballpark is annual F (=Z-M) for specified year"<<endl;
  init_int F_ballpark_yr
  !! echoinput<<F_ballpark_yr<<" F_ballpark_yr (<0 to ignore)  "<<endl;
  init_int F_Method;           // 1=Pope's; 2=continuouos F; 3=hybrid
  int F_Method_use
  !! echoinput<<F_Method<<" F_Method "<<endl;
  init_number max_harvest_rate
  number Equ_F_joiner

 LOCAL_CALCS
    echoinput<<max_harvest_rate<<" max_harvest_rate "<<endl;
  if(F_Method<1 || F_Method>3)
    {
      N_warn++;
    warning<<" ERROR:  F_Method must be 1 or 2 or 3, value is: "<<F_Method<<endl;
    cout<<" EXIT - see warning "<<endl;
    exit(1);
    }
   if(F_Method==1)
   {
     k=-1;
     j=-1;
     Equ_F_joiner=(log(1./max_harvest_rate -1.))/(max_harvest_rate-0.2);  //  used to spline the harvest rate
     if(max_harvest_rate>0.999)
     {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" max harvest rate must  be <1.0 for F_method 1 "<<max_harvest_rate<<endl; exit(1);}
     if(max_harvest_rate<=0.30)
     {N_warn++; warning<<" unexpectedly small value for max harvest rate for F_method 1:  "<<max_harvest_rate<<endl;}
   }
   else
   {
     if(max_harvest_rate<1.0)
     {N_warn++; warning<<" max harvest rate should be >1.0 for F_method 2 or 3 "<<max_harvest_rate<<endl;}
     if(F_Method==2)
     {
       k=3;
       j=Nfleet*(TimeMax-styr+1);
     }
     else
     {
       k=1;
       j=-1;
     }
   }
 END_CALCS

  init_vector F_setup(1,k)
  vector F_rate_max(1,j)
// setup for F_rate with F_Method=2
// F_setup(1) = overall initial value
// F_setup(2) = overall phase
// F_setup(3) = number of specific initial values and phases to read
  int F_detail
  int F_Tune
 LOCAL_CALCS
  F_detail=-1;
  if(F_Method>1)
  {
    if(F_Method==2)
    {
      echoinput<<F_setup<<" initial F value, F phase, N_detailed Fsetups to read "<<endl;
      F_detail=F_setup(3);
      F_Tune=4;
      F_rate_max=max_harvest_rate;  // used to set upper bound on F_rate parameter
    }
    else if(F_Method==3)
    {
      F_Tune=F_setup(1);
      echoinput<<F_Tune<<" N iterations for tuning hybrid F "<<endl;
    }
  }
 END_CALCS

  init_matrix F_setup2(1,F_detail,1,6)  // fleet, yr, seas, Fvalue, se, phase
  !!echoinput<<" detailed F_setups "<<endl<<F_setup2<<endl;

//  ProgLabel_4.7.1  read setup for init_F parameters and create init_F parameter labels
  init_matrix init_F_parm_1(1,Nfleet,1,7)
  !! echoinput<<" initial equil F parameter setup"<<endl<<init_F_parm_1<<endl;
  vector init_F_LO(1,Nfleet)
  vector init_F_HI(1,Nfleet)
  vector init_F_RD(1,Nfleet)
  vector init_F_PR(1,Nfleet)
  vector init_F_PRtype(1,Nfleet)
  vector init_F_CV(1,Nfleet)
  ivector init_F_PH(1,Nfleet)

 LOCAL_CALCS
   init_F_LO=column(init_F_parm_1,1);
   init_F_HI=column(init_F_parm_1,2);
   init_F_RD=column(init_F_parm_1,3);
   init_F_PR=column(init_F_parm_1,4);
   init_F_PRtype=column(init_F_parm_1,5);
   init_F_CV=column(init_F_parm_1,6);
   init_F_PH=ivector(column(init_F_parm_1,7));
   for (f=1;f<=Nfleet;f++)
     {
      ParCount++; ParmLabel+="InitF_"+NumLbl(f)+fleetname(f);
      if(obs_equ_catch(f)<=0.0)
      {
        if(init_F_RD(f)>0.0)
        {
          N_warn++;
          warning<<f<<" catch: "<<obs_equ_catch(f)<<" initF: "<<init_F_RD(f)<<" initF is reset to be 0.0"<<endl;
        }
        init_F_RD(f)=0.0; init_F_PH(f)=-1;
        }
      if(obs_equ_catch(f)>0.0 && init_F_RD(f)<=0.0)
       {N_warn++; cout<<" EXIT - see warning "<<endl;
        warning<<f<<" catch: "<<obs_equ_catch(f)<<" initF: "<<init_F_RD(f)<<" initF must be >0"<<endl; exit(1);}
     }

//  ProgLabel_4.7.2  create parameter labels for F parameters if F_method==2
  if(F_Method==2)
  {
    for(f=1;f<=Nfleet;f++)
    for(y=styr;y<=endyr;y++)
    for(s=1;s<=nseas;s++)
    {
//      _itoa(y,onenum,10);
      sprintf(onenum, "%d", y);
      ParCount++;
      ParmLabel+="F_fleet_"+NumLbl(f)+"_YR_"+onenum+"_s_"+NumLbl(s)+CRLF(1);
    }
      g=0;
      for(f=1;f<=Nfleet;f++)
      for(y=styr;y<=endyr;y++)
      for(s=1;s<=nseas;s++)
      {
        t=styr+(y-styr)*nseas+s-1;
        g++;
        if(catch_ret_obs(f,t)>0.)
        {
          have_catch(g)=F_setup(2);
        }
        else
        {
          have_catch(g)=-1;
        }
      }
      if(F_detail>0)
      {
        for(k=1;k<=F_detail;k++)
        {
          f=F_setup2(k,1); y=F_setup2(k,2); s=F_setup2(k,3);
          g=(f-1)*(TimeMax-styr+1)+(y-styr)*nseas+s;
          t=styr+(y-styr)*nseas+s-1;
          if(F_setup2(k,6)!=-999) have_catch(g)=F_setup2(k,6);    //   used to setup the phase for F_rate
          if(F_setup2(k,5)!=-999) catch_se(t,f)=F_setup2(k,5);    //    reset the se for this observation
          //  setup of F_rate values occurs later in the parameter section
        }
      }
  }
 END_CALCS


//***********************************
//  ProgLabel_4.8  read catchability (Q) setup
  init_matrix Q_setup(1,Ntypes,1,4)  // do power, env-var,  extra sd, devtype(<0=mirror, 0=float_nobiasadj 1=float_biasadj, 2=parm_nobiasadj, 3=rand, 4=randwalk); num/bio/F, err_type(0=lognormal, >=1 is T-dist-lognormal)
                                        // change to matrix because devstd has real, not integer, values
  int Q_Npar2
  int Q_Npar
  int ask_detail

  imatrix Q_setup_parms(1,Ntypes,1,4)
 LOCAL_CALCS
  echoinput<<" Q setup "<<endl<<Q_setup<<endl;
  echoinput<<"Note that Q parameters are in units of ln(q)"<<endl;
  Q_Npar=0;
  ask_detail=0;
//  ProgLabel_4.8.1  create index to the catchability parameters and create parameter names
  for (f=1;f<=Ntypes;f++)
  {
    Q_setup_parms(f,1)=0;
   if(Q_setup(f,1)>0)
    {
      Q_Npar++; Q_setup_parms(f,1)=Q_Npar;
      ParCount++; ParmLabel+="Q_power_"+NumLbl(f)+"_"+fleetname(f);
      if(Q_setup(f,4)<2) {N_warn++; warning<<" must create base Q parm to use Q_power for fleet: "<<f<<endl;}
    }
  }
  for (f=1;f<=Ntypes;f++)
  {
    Q_setup_parms(f,2)=0;
    if(Q_setup(f,2)!=0)
      {
        Q_Npar++; Q_setup_parms(f,2)=Q_Npar;
        ParCount++; ParmLabel+="Q_envlink_"+NumLbl(f)+"_"+fleetname(f);
        if(Q_setup(f,4)<2) {N_warn++; warning<<" must create base Q parm to use Q_envlink for fleet: "<<f<<endl;}
       }
  }
  for (f=1;f<=Ntypes;f++)
  {
    Q_setup_parms(f,3)=0;
    if(Q_setup(f,3)>0)
    {
      Q_Npar++; Q_setup_parms(f,3)=Q_Npar;
      ParCount++; ParmLabel+="Q_extraSD_"+NumLbl(f)+"_"+fleetname(f);
    }
  }
//  ProgLabel_4.8.2  create time-varying catchability as needed
  Q_Npar2=Q_Npar;
  for (f=1;f<=Ntypes;f++)
  {
    Q_setup_parms(f,4)=0;
    if(Q_setup(f,4)>=2)
    {
      Q_Npar++; Q_Npar2++; Q_setup_parms(f,4)=Q_Npar;
      ParCount++; ParmLabel+="LnQ_base_"+NumLbl(f)+"_"+fleetname(f);
      if(Q_setup(f,4)==3)
      {
        ask_detail=1;
        Q_Npar2++;
        Q_Npar+=nyr_cr(f);
        for(j=1;j<=nyr_cr(f);j++)
        {
          y=yr_cr_y(f,j);
          ParCount++;
          sprintf(onenum, "%d", y);
          onenum+=CRLF(1);
          if(nseas==1)
          {
            ParmLabel+="Q_dev_"+NumLbl(f)+"y_"+onenum;
          }
          else
          {
            ParmLabel+="Q_dev_"+NumLbl(f)+"y_"+onenum+"_s_"+NumLbl(abs(yr_cr_s(f,j)));
          }

        }
      }
      if(Q_setup(f,4)==4)
      {
        ask_detail=1;
        Q_Npar2++;
        Q_Npar+=nyr_cr(f)-1;
        for(j=2;j<=nyr_cr(f);j++)
        {
          y=yr_cr_y(f,j);
          ParCount++;
//          _itoa(y,onenum,10);
          sprintf(onenum, "%d", y);
          onenum+=CRLF(1);
          if(nseas==1)
          {
            ParmLabel+="Q_walk_"+NumLbl(f)+"y_"+onenum;
          }
          else
          {
            ParmLabel+="Q_walk_"+NumLbl(f)+"y_"+onenum+"_s_"+NumLbl(abs(yr_cr_s(f,j)));
          }
        }
      }
    }
    else if(cr_errtype(f)==-1)
    {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Error, cannot use scaling approach to Q if error type is normal "<<endl; exit(1);}
  }

  for (f=1;f<=Ntypes;f++)
  {
    if(cr_units(f)==2)  // effort deviations
    {
      if(cr_errtype(f)>=0)  //  lognormal
      {
        N_warn++;
        warning<<" Lognormal error selected for effort deviations for fleet "<<f<<"; normal error recommended"<<endl;
      }
      if(Q_setup(f,1)>0)  //  density-dependence
      {
        N_warn++;
        warning<<" Do not use Density-dependence for effort deviations (fleet "<<f<<"); "<<endl;
      }
    }
  }

  if(Q_Npar>0)
    {k=Q_Npar;}
  else
    {k=1;}
 END_CALCS

  vector Q_parm_LO(1,k)
  vector Q_parm_HI(1,k)
  ivector Q_parm_PH(1,k)

  int Q_parm_detail
 LOCAL_CALCS
  if(ask_detail>0)
  {
    *(ad_comm::global_datafile) >> Q_parm_detail;
    echoinput<<Q_parm_detail<<" Q_parm detail for time-varying parameters "<<endl;
  }
  else
  {
    Q_parm_detail=0;
    echoinput<<" # No time-varying Q parms, so no q_parm_detail input needed "<<endl;
  }
  if(Q_parm_detail==1) {j=Q_Npar;} else {j=Q_Npar2;}
 END_CALCS

 matrix Q_parm_1(1,Q_Npar,1,7)

//  ProgLabel_4.8.3  read catchability parameters as necessary
  init_matrix Q_parm_2(1,j,1,7)
 LOCAL_CALCS
  if(j>0) echoinput<<" Catchability parameters"<<endl<<Q_parm_2<<endl;
  if(Q_parm_detail==0)
  {
    Q_Npar=0;  Q_Npar2=0;
    for (f=1;f<=Ntypes;f++)
    {
     if(Q_setup(f,1)>0)
      {
        Q_Npar++;
        Q_parm_1(Q_Npar)=Q_parm_2(Q_Npar);
      }
    }
    for (f=1;f<=Ntypes;f++)
    {
      if(Q_setup(f,2)!=0)
      {
        Q_Npar++;
        Q_parm_1(Q_Npar)=Q_parm_2(Q_Npar);
      }
    }
    for (f=1;f<=Ntypes;f++)
    {
     if(Q_setup(f,3)>0)
      {
        Q_Npar++;
        Q_parm_1(Q_Npar)=Q_parm_2(Q_Npar);
      }
    }
    Q_Npar2=Q_Npar;
    for (f=1;f<=Ntypes;f++)
    {
      if(Q_setup(f,4)>=2)
      {
        Q_Npar++; Q_Npar2++;
        Q_parm_1(Q_Npar)=Q_parm_2(Q_Npar2);
        if(Q_setup(f,4)==3)
        {
          Q_Npar2++;
          for(j=1;j<=nyr_cr(f);j++)
          {
            Q_Npar++;
            Q_parm_1(Q_Npar)=Q_parm_2(Q_Npar2);
          }
        }
        if(Q_setup(f,4)==4)
        {
          Q_Npar2++;
          for(j=2;j<=nyr_cr(f);j++)
          {
            Q_Npar++;
            Q_parm_1(Q_Npar)=Q_parm_2(Q_Npar2);
          }
        }
      }
    }
  }
  else
  {
    Q_parm_1=Q_parm_2;
  }
 END_CALCS

  !! if(Q_Npar>0 ) echoinput<<" processed Q parms "<<endl<<Q_parm_1<<endl;

 LOCAL_CALCS
   if(Q_Npar>0)
     {
     for(f=1;f<=Q_Npar;f++)
       {
       Q_parm_LO(f)=Q_parm_1(f,1);
       Q_parm_HI(f)=Q_parm_1(f,2);
       Q_parm_PH(f)=Q_parm_1(f,7);
       }
     }
   else
     {Q_parm_LO=-1.; Q_parm_HI=1.; Q_parm_PH=-4;}
 END_CALCS


//******************************************
//  ProgLabel_4.9  Define Selectivity patterns and N parameters needed per pattern
  ivector seltype_Nparam(0,34)
 LOCAL_CALCS
   seltype_Nparam(0)=0;   // selex=1.0 for all sizes
   seltype_Nparam(1)=2;   // logistic; with 95% width specification
   seltype_Nparam(2)=8;   // double logistic, with defined peak
   seltype_Nparam(3)=6;   // flat middle, power up, power down
   seltype_Nparam(4)=0;   // set size selex=female maturity
   seltype_Nparam(5)=2;   // mirror another selex; PARMS pick the min-max bin to mirror
   seltype_Nparam(6)=2;   // non-parm len selex, additional parm count is in seltype(f,4)
   seltype_Nparam(7)=8;   // New doublelogistic with smooth transitions and constant above Linf option
   seltype_Nparam(8)=8;   // New doublelogistic with smooth transitions and constant above Linf option
   seltype_Nparam(9)=6;   // simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset

   seltype_Nparam(10)=0;   //  First age-selex  selex=1.0 for all ages
   seltype_Nparam(11)=2;   //  pick min-max age
   seltype_Nparam(12)=2;   //   logistic
   seltype_Nparam(13)=8;   //   double logistic
   seltype_Nparam(14)=nages+1;   //   empirical
   seltype_Nparam(15)=0;   //   mirror another selex
   seltype_Nparam(16)=2;   //   Coleraine - Gaussian
   seltype_Nparam(17)=nages+1;   //   empirical as random walk  N parameters to read can be overridden by setting special to non-zero
   seltype_Nparam(18)=8;   //   double logistic - smooth transition
   seltype_Nparam(19)=6;   //   simple 4-parm double logistic with starting age
   seltype_Nparam(20)=6;   //   double_normal,using joiners

   seltype_Nparam(21)=2;   // non-parm len selex, additional parm count is in seltype(f,4), read as pairs of size, then selex
   seltype_Nparam(22)=4;   //   double_normal as in CASAL
   seltype_Nparam(23)=6;   //   double_normal where final value is directly equal to sp(6) so can be >1.0
   seltype_Nparam(24)=6;   //   double_normal with sel(minL) and sel(maxL), using joiners
   seltype_Nparam(25)=3;   //   exponential-logistic in size
   seltype_Nparam(26)=3;   //   exponential-logistic in age
   seltype_Nparam(27)=3;   // cubic spline for selex at length, additional parm count is in seltype(f,4)
//   seltype_Nparam(28)=3;   // cubic spline for selex at age, additional parm count is in seltype(f,4)
   seltype_Nparam(29)=0;   //   undefined
   seltype_Nparam(30)=0;   //   spawning biomass
   seltype_Nparam(31)=0;   //   recruitment dev
   seltype_Nparam(32)=0;   //   pre-recruitment (spawnbio * recrdev)
   seltype_Nparam(33)=0;   //   recruitment
   seltype_Nparam(34)=0;   //   spawning biomass depletion

 END_CALCS

//  ProgLabel_4.9.1  read selectivity definitions
//  do 2*Ntypes to create options for size-selex (first), then age-selex
  init_imatrix seltype(1,2*Ntypes,1,4)    // read selex type for each fleet/survey, Do_retention, Do_male
  !! echoinput<<" selex types "<<endl<<seltype<<endl;
  int N_selparm   // figure out the Total number of selex parameters
  int N_selparm2                 // N selparms plus env links and blocks
  ivector N_selparmvec(1,2*Ntypes)  //  N selparms by type, including extra parms for male selex, retention, etc.
  ivector Maleselparm(1,2*Ntypes)
  ivector RetainParm(1,Ntypes)
  ivector dolen(1,Ntypes)
  int blkparm
  int firstselparm
  int Do_Retain

 LOCAL_CALCS
//  ProgLabel_4.9.2  process selectivity parameter count and create parameter labels
  depletion_fleet=0;
   firstselparm=ParCount;
   N_selparm=0;
   Do_Retain=0;
   for (f=1;f<=Ntypes;f++)
   {
     if(WTage_rd>0 && seltype(f,1)>0)
     {
      N_warn++; warning<<" Use of size selectivity not advised when reading empirical wt-at-age "<<endl;
     }
     N_selparmvec(f)=seltype_Nparam(seltype(f,1));   // N Length selex parms
     if(seltype(f,1)==6) N_selparmvec(f) +=seltype(f,4);  // special setup of N parms
     if(seltype(f,1)==21) N_selparmvec(f) +=2*(seltype(f,4)-1);  // special setup of N parms
     if(seltype(f,1)==27) N_selparmvec(f) +=2*seltype(f,4);  // special setup of N parms for cubic spline
     if(seltype(f,1)>0 && seltype(f,1)<30) {dolen(f)=1;} else {dolen(f)=0;}

     if(seltype(f,1)==27)
     {
         ParCount++; ParmLabel+="SizeSpline_Code_"+fleetname(f)+"_"+NumLbl(f);
         ParCount++; ParmLabel+="SizeSpline_GradLo_"+fleetname(f)+"_"+NumLbl(f);
         ParCount++; ParmLabel+="SizeSpline_GradHi_"+fleetname(f)+"_"+NumLbl(f);
         for(s=1;s<=seltype(f,4);s++)
         {
           ParCount++; ParmLabel+="SizeSpline_Knot_"+NumLbl(s)+"_"+fleetname(f)+"_"+NumLbl(f);
         }
         for(s=1;s<=seltype(f,4);s++)
         {
           ParCount++; ParmLabel+="SizeSpline_Val_"+NumLbl(s)+"_"+fleetname(f)+"_"+NumLbl(f);
         }
     }
     else
     {
       for(j=1;j<=N_selparmvec(f);j++)
       {
         ParCount++; ParmLabel+="SizeSel_"+NumLbl(f)+"P_"+NumLbl(j)+"_"+fleetname(f);
       }
     }

     if(seltype(f,1)==34)  //  special code for depletion, so adjust phases and lambdas
      {
        depletion_fleet=f;
      }
      
     if(seltype(f,2)>=1)
     {
       if(WTage_rd>0)
       {
        N_warn++; warning<<" BEWARE: Retention functions not implemented fully when reading empirical wt-at-age "<<endl;
       }
       Do_Retain=1;
       if(seltype(f,2)==3)
       {RetainParm(f)=0;}  //  no parameters needed
       else
       {
       RetainParm(f)=N_selparmvec(f)+1;
       N_selparmvec(f) +=4*seltype(f,2);          // N retention parms first 4 for retention; next 4 for mortality
       for(j=1;j<=4;j++)
       {
         ParCount++; ParmLabel+="Retain_"+NumLbl(f)+"P_"+NumLbl(j)+"_"+fleetname(f);
       }
       if(seltype(f,2)==2)
       {
         for(j=1;j<=4;j++)
         {
           ParCount++; ParmLabel+="DiscMort_"+NumLbl(f)+"P_"+NumLbl(j)+"_"+fleetname(f);
         }
       }
      }
     }
     if(seltype(f,3)>=1)
      {
        if(gender==1) {N_warn++; cout<<"Critical error"<<endl; warning<<" Male selex cannot be used in one sex model; fleet: "<<f<<endl; exit(1);}
        Maleselparm(f)=N_selparmvec(f)+1;
        if(seltype(f,3)==1 || seltype(f,3)==2)
        {
          N_selparmvec(f)+=4;  // add male parms
          ParCount+=4;
          ParmLabel+="SzSel_"+NumLbl(f)+"MaleDogleg_"+fleetname(f);
          ParmLabel+="SzSel_"+NumLbl(f)+"MaleatZero_"+fleetname(f);
          ParmLabel+="SzSel_"+NumLbl(f)+"MaleatDogleg_"+fleetname(f);
          ParmLabel+="SzSel_"+NumLbl(f)+"MaleatMaxage_"+fleetname(f);
        }
        else if(seltype(f,3)>=3)
        {
          if(seltype(f,3)==3) {anystring="Male_";} else {anystring="Fem_";}
          if(seltype(f,1)==1)
          {
            N_selparmvec(f)++; ParCount++; ParmLabel+="SzSel_"+NumLbl(f)+anystring+"Infl_"+fleetname(f);
            N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f)+anystring+"Slope_"+fleetname(f);
            N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f)+anystring+"Scale_"+fleetname(f);
          }
          else if(seltype(f,1)==24)
          {
            N_selparmvec(f)++; ParCount++; ParmLabel+="SzSel_"+NumLbl(f)+anystring+"Peak_"+fleetname(f);
            N_selparmvec(f)++; ParCount++; ParmLabel+="SzSel_"+NumLbl(f)+anystring+"Ascend_"+fleetname(f);
            N_selparmvec(f)++; ParCount++; ParmLabel+="SzSel_"+NumLbl(f)+anystring+"Descend_"+fleetname(f);
            N_selparmvec(f)++; ParCount++; ParmLabel+="SzSel_"+NumLbl(f)+anystring+"Final_"+fleetname(f);
            N_selparmvec(f)++; ParCount++; ParmLabel+="SzSel_"+NumLbl(f)+anystring+"Scale_"+fleetname(f);
          }
          else
          {
            N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Illegal male selex option selected for fleet "<<f<<endl;  exit(1);
          }
        }
      }

     if(seltype(f,1)==7) {N_warn++; warning<<"ERROR:  selectivity pattern #7 is no longer supported "<<endl;}
     if(seltype(f,1)==23 && F_Method==1) {N_warn++; warning<<"Do not use F_Method = Pope's with selex pattern #23 "<<endl;}
     N_selparm += N_selparmvec(f);
   }
   for (f=Ntypes+1;f<=2*Ntypes;f++)
   {
     if(seltype(f,1)==15) // mirror
     {
       if(seltype(f,4)==0 || seltype(f,4)>=f-Ntypes)
       {
         N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" illegal mirror for age selex fleet "<<f-Ntypes<<endl; exit(1);
       }
       N_selparmvec(f)=0;   // Nunber of Age selex parms
     }
     else if(seltype(f,1)!=17)
     {
       N_selparmvec(f)=seltype_Nparam(seltype(f,1));   // Nunber of Age selex parms
     }
     else if(seltype(f,4)==0)
     {
       N_selparmvec(f)=seltype_Nparam(seltype(f,1));   // this is nages+1
     }
     else
     {
       N_selparmvec(f)=abs(seltype(f,4))+1;   // so reads value for age 0 through this age
     }

     if(seltype(f,1)==27)
     {
       N_selparmvec(f) +=2*seltype(f,4);  // special setup of N parms for cubic spline
       ParCount++; ParmLabel+="AgeSpline_Code_"+fleetname(f-Ntypes)+"_"+NumLbl(f-Ntypes);
       ParCount++; ParmLabel+="AgeSpline_GradLo_"+fleetname(f-Ntypes)+"_"+NumLbl(f-Ntypes);
       ParCount++; ParmLabel+="AgeSpline_GradHi_"+fleetname(f-Ntypes)+"_"+NumLbl(f-Ntypes);
       for(s=1;s<=seltype(f,4);s++)
       {
         ParCount++; ParmLabel+="AgeSpline_Knot_"+NumLbl(s)+"_"+fleetname(f-Ntypes)+"_"+NumLbl(f-Ntypes);
       }
       for(s=1;s<=seltype(f,4);s++)
       {
         ParCount++; ParmLabel+="AgeSpline_Val_"+NumLbl(s)+"_"+fleetname(f-Ntypes)+"_"+NumLbl(f-Ntypes);
       }
     }
     else
     {
       for(j=1;j<=N_selparmvec(f);j++)
       {
         ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+"P_"+NumLbl(j)+"_"+fleetname(f-Ntypes);
       }
     }
     if(seltype(f,3)>=1)
      {
        if(gender==1) {N_warn++; cout<<"Critical error"<<endl; warning<<" Male selex cannot be used in one sex model; fleet: "<<f<<endl; exit(1);}
        Maleselparm(f)=N_selparmvec(f)+1;
        if(seltype(f,3)==1 || seltype(f,3)==2)
        {
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+"MaleDogleg_"+fleetname(f-Ntypes);
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+"MaleatZero_"+fleetname(f-Ntypes);
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+"MaleatDogleg_"+fleetname(f-Ntypes);
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+"MaleatMaxage_"+fleetname(f-Ntypes);
        }
        else if(seltype(f,3)>=3 && seltype(f,1)==20)
        {
          if(seltype(f,3)==3) {anystring="Male_";} else {anystring="Fem_";}
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+anystring+"Peak_"+fleetname(f-Ntypes);
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+anystring+"Ascend_"+fleetname(f-Ntypes);
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+anystring+"Descend_"+fleetname(f-Ntypes);
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+anystring+"Final_"+fleetname(f-Ntypes);
          N_selparmvec(f)++; ParCount++; ParmLabel+="AgeSel_"+NumLbl(f-Ntypes)+anystring+"Scale_"+fleetname(f-Ntypes);
        }
        else
        {
          N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Illegal male selex option selected for fleet "<<f<<endl;  exit(1);
        }
      }
     N_selparm += N_selparmvec(f);
   }
   for(f=1;f<=Nfleet;f++)
     {
     if(nyr_disc(f)>0 && seltype(f,2)==0)
       {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" ERROR:  discard data exist for fleet "<<f<<"  but retention parms not setup "<<endl; exit(1);}
     else if (nyr_disc(f)==0 && seltype(f,2)!=0)
       {N_warn++; warning<<" WARNING:  no discard amount data for fleet "<<f<<"  but retention parms have been defined "<<endl;}
     }
 END_CALCS

//  ProgLabel_4.9.3  read selex parameters
  init_matrix selparm_1(1,N_selparm,1,14)
 LOCAL_CALCS
  echoinput<<" selex base parameters "<<endl;
  for(g=1;g<=N_selparm;g++)
  {
    echoinput<<g<<" ## "<<selparm_1(g)<<" ## "<<ParmLabel(ParCount-N_selparm+g)<<endl;
  }
 END_CALCS

  imatrix time_vary_sel(styr-3,endyr+1,1,2*Ntypes)
  imatrix time_vary_makefishsel(styr-3,endyr+1,1,Nfleet)
  int makefishsel_yr
  int ALK_yr
//  ProgLabel_4.9.4  create and label environmental linkages for selectivity parameters
  int N_selparm_env                            // number of selparms that use env linkage
  int customenvsetup  //  0=read one setup and apply to all; 1=read each
  ivector selparm_env(1,N_selparm)             //  pointer to parameter with env link for each selparm
  ivector selparm_envuse(1,N_selparm)   // contains the environment data number
  ivector selparm_envtype(1,N_selparm)  // 1=multiplicative; 2= additive

 LOCAL_CALCS
  N_selparm_env=0;
  selparm_env.initialize();
  selparm_envtype.initialize();
  selparm_envuse.initialize();

  for(j=1;j<=N_selparm;j++)
  {
    if(selparm_1(j,8)!=0)
    {
      N_selparm_env++; selparm_env(j)=N_selparm+N_selparm_env;
      if(selparm_1(j,8)>0)
      {
        ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_ENV_mult"; selparm_envtype(j)=1; selparm_envuse(j)=selparm_1(j,8);
       }
       else if(selparm_1(j,8)==-999)
       {ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_ENV_densdep"; selparm_envtype(j)=3;  MGparm_envuse(j)=-1;}
       else
       {ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_ENV_add"; selparm_envtype(j)=2; selparm_envuse(j)=-selparm_1(j,8);}
    }
  }

  if(N_selparm_env>0)
  {
    *(ad_comm::global_datafile) >> customenvsetup;
    if(customenvsetup==0) {k1=1;} else {k1=N_selparm_env;}
    echoinput<<customenvsetup<<" customenvsetup"<<endl;
  }
  else
  {customenvsetup=0; k1=0;
    echoinput<<" no envlinks; so don't read customenvsetup"<<endl;
    }
 END_CALCS

  init_matrix selparm_env_1(1,k1,1,7)  // read matrix that sets up the env linkage parms
 LOCAL_CALCS
  if(k1>0)
    {echoinput<<" selex-env parameters "<<endl;
      for(g=1;g<=k1;g++)
      {
        echoinput<<g<<" ## "<<selparm_env_1(g)<<" ## "<<ParmLabel(ParCount-k1+g)<<endl;
      }
    }
 END_CALCS

//  ProgLabel_4.9.5  creae and label block patterns for selectivity parameters
  int N_selparm_blk                            // number of selparms that use blocks
  imatrix Block_Defs_Sel(1,N_selparm,styr,endyr)
  int customblocksetup  //  0=read one setup and apply to all; 1=read each

  int N_selparm_trend     //   number of selex parameters using trend
  int N_selparm_trend2     //   number of parameters needed to define trends and cycles
  ivector selparm_trend_point(1,N_selparm)   //  index of trend parameters associated with each selex parm

 LOCAL_CALCS
  Block_Defs_Sel.initialize();
  N_selparm_blk=0;  // counter for assigned parms
  for(j=1;j<=N_selparm;j++)
  {
    z=selparm_1(j,13);    // specified block definition
    if(z>N_Block_Designs) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" ERROR, Block > N Blocks "<<z<<" "<<N_Block_Designs<<endl; exit(1);}
    if(z>0)
    {
      g=1;
      for(a=1;a<=Nblk(z);a++)
      {
        N_selparm_blk++;
        y=Block_Design(z,g);
        sprintf(onenum, "%d", y);
        ParCount++;
        k=selparm_1(j,14);
        switch(k)
        {
          case 0:
          {ParmLabel+=ParmLabel(j+firstselparm)+"_BLK"+NumLbl(z)+"mult_"+onenum+CRLF(1);  break;}
          case 1:
          {ParmLabel+=ParmLabel(j+firstselparm)+"_BLK"+NumLbl(z)+"add_"+onenum+CRLF(1);  break;}
          case 2:
          {ParmLabel+=ParmLabel(j+firstselparm)+"_BLK"+NumLbl(z)+"repl_"+onenum+CRLF(1);  break;}
          case 3:
          {ParmLabel+=ParmLabel(j+firstselparm)+"_BLK"+NumLbl(z)+"delta_"+onenum+CRLF(1);  break;}
        }
        for(y=Block_Design(z,g);y<=Block_Design(z,g+1);y++)  // loop years for this block
        {
         if(y<=endyr) Block_Defs_Sel(j,y)=N_selparm+N_selparm_env+N_selparm_blk;
        }
        g+=2;
      }
    }
  }

  if(N_selparm_blk>0)
  {
    *(ad_comm::global_datafile) >> customblocksetup;
    if(customblocksetup==0) {k1=1;} else {k1=N_selparm_blk;}
    echoinput<<customblocksetup<<" customblocksetup"<<endl;
  }
  else
  {customblocksetup=0; k1=0;
    echoinput<<" no blocks; so don't read customblocksetup"<<endl;
  }
 END_CALCS

  init_matrix selparm_blk_1(1,k1,1,7);  // read matrix that defines the block parms and trend parms

 LOCAL_CALCS
  if(k1>0)
    {
      echoinput<<" selex-block parameters "<<endl;
      for(g=1;g<=k1;g++)
      {
        echoinput<<g<<" ## "<<selparm_blk_1(g)<<" ## "<<ParmLabel(ParCount-k1+g)<<endl;
      }
    }

// use negative block as indicator to use time trend
//  ProgLabel_4.9.6  create and label parameter trends or cycles as alternative to blocks
  N_selparm_trend=0;
  N_selparm_trend2=0;
  for(j=1;j<=N_selparm;j++)
  {
    if(selparm_1(j,13)<0)  //  create timetrend parameter
    {
      N_selparm_trend++;
      selparm_trend_point(j)=N_selparm_trend;
      if(selparm_1(j,13)==-1)
      {
        ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_TrendFinal_Offset"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_TrendInfl_"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_TrendWidth_"+CRLF(1);
        N_selparm_trend2+=3;
      }
      else if(selparm_1(j,13)==-2)
      {
        ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_TrendFinal_"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_TrendInfl_"+CRLF(1);
        ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_TrendWidth_"+CRLF(1);
        N_selparm_trend2+=3;
      }
      else
      {
        for(icycle=1;icycle<=Ncycle;icycle++)
        {
          ParCount++; ParmLabel+=ParmLabel(j+firstselparm)+"_Cycle_"+NumLbl(icycle)+CRLF(1);
          N_selparm_trend2++;
        }
      }
    }
  }
  if(N_selparm_trend2>0) echoinput<<" Create N selparm_trend "<<N_selparm_trend<<endl;
 END_CALCS

  init_matrix selparm_trend_1(1,N_selparm_trend2,1,7)  // read matrix that defines the parms and trend parms

 LOCAL_CALCS
  if(N_selparm_trend2>0)
    {
      echoinput<<"Selex trend and cycle parameters "<<endl;
      for(g=1;g<=N_selparm_trend2;g++)
      {
        echoinput<<g<<" ## "<<selparm_trend_1(g)<<" ## "<<ParmLabel(ParCount-N_selparm_trend2+g)<<endl;
      }
    }
 END_CALCS

  ivector selparm_trend_rev(1,N_selparm_trend)
  ivector selparm_trend_rev_1(1,N_selparm_trend)

 LOCAL_CALCS
  if(N_selparm_trend>0)
  {
    k1=0;
    k2=N_selparm+N_selparm_env+N_selparm_blk;
    for(j=1;j<=N_selparm;j++)
    {
      if(selparm_1(j,13)<0)  //  create timetrend parameter
      {
        k1++;
        selparm_trend_rev(k1)=j;
        selparm_trend_rev_1(k1)=k2;  // pointer to base in list of MGparms (so k2+1 is first parameter used)
        if(selparm_1(j,13)>=-2)  //  timetrend
        {k2+=3;}
        else
        {k2+=Ncycle;}
      }
    }
  }
 END_CALCS

//  ProgLabel_4.9.7  create and label selectivity parameter annual devs
  int N_selparm_dev   // number of selparms that use random deviations
  int N_selparm_dev_tot   // number of selparms that use random deviations
 LOCAL_CALCS
  N_selparm_dev=0;
  N_selparm_dev_tot=0;
  for(j=1;j<=N_selparm;j++)
  {
    if(selparm_1(j,9)!=0)
    {
      N_selparm_dev++;
      for(y=selparm_1(j,10);y<=selparm_1(j,11);y++)
      {
        N_selparm_dev_tot++;
        sprintf(onenum, "%d", y);
        ParCount++;
        if(selparm_1(j,9)==1)
        {ParmLabel+=ParmLabel(j+firstselparm)+"_DEVmult_"+onenum+CRLF(1);}
        else if(selparm_1(j,9)==2)
        {ParmLabel+=ParmLabel(j+firstselparm)+"_DEVadd_"+onenum+CRLF(1);}
        else if(selparm_1(j,9)==3)
        {ParmLabel+=ParmLabel(j+firstselparm)+"_DEVrwalk_"+onenum+CRLF(1);}
        else
        {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" illegal selparmdevtype for parm "<<j<<endl; exit(1);}
      }
    }
  }
  if(N_selparm_dev>0) j=N_selparm_dev; else j=1;    // for defining non-zero array size
 END_CALCS

  ivector selparm_dev_minyr(1,j)
  ivector selparm_dev_maxyr(1,j)
  vector  selparm_dev_stddev(1,j)
  ivector  selparm_dev_type(1,j)
  ivector  selparm_dev_select(1,N_selparm)
  number selparm_dev_PH
  int selparm_adjust_method   //  1=do V1.xx approach to adjustment by env, block or dev; 2=use new logistic approach; 3=no check

 LOCAL_CALCS
  selparm_dev_select.initialize();
  if(N_selparm_dev==0)
  {
    selparm_dev_PH=-6;
    echoinput<<" No selparm devs selected, so don't read selparm_dev_PH"<<endl;
  }
  else
  {
    *(ad_comm::global_datafile) >> selparm_dev_PH;
    echoinput<<selparm_dev_PH<<" selparm_dev_PH"<<endl;
  }

  if(N_selparm_env+N_selparm_blk+N_selparm_dev > 0)
  {
    *(ad_comm::global_datafile) >> selparm_adjust_method;
    echoinput<<selparm_adjust_method<<" selparm_adjust_method"<<endl;
    if(selparm_adjust_method<1 || selparm_adjust_method>3)
    {
      N_warn++; cout<<" EXIT - see warning "<<endl;
      warning<<" illegal selparm_adjust_method; must be 1 or 2 or 3 "<<endl;  exit(1);
    }
  }
  else
  {
    selparm_adjust_method=0;
    echoinput<<" No selparm adjustments, so don't read selparm_adjust_method"<<endl;
  }

//  ProgLabel_4.9.8  create some indexes need to track selex parms and do some error checking
    j=0;
    for (f=1;f<=N_selparm;f++)
    {
      if(selparm_1(f,9)>0)
      {
        j++;
        selparm_dev_type(j)=selparm_1(f,9);  // 1 for mult; 2 for additive; 3 for additive randwalk
        selparm_dev_select(f)=j;  // pointer to dev vector used by this parameter
        selparm_dev_minyr(j)=selparm_1(f,10);
        selparm_dev_maxyr(j)=selparm_1(f,11);
        selparm_dev_stddev(j)=selparm_1(f,12);
        if(selparm_dev_type(j)==2 && selparm_dev_stddev(j)<0.10*selparm_1(f,3))
        {N_warn++; warning<<" selparm_dev_stddev is small (<10% of parm value) for selparm: "<<f<<endl;}
        if(selparm_adjust_method==2 && selparm_1(f,9)==1)
        {N_warn++; warning<<" cannot use selparm_adjust_method==2 and multiplicative devs for selex parameter; STOP "<<f<<endl; exit(1);}
        if(selparm_dev_minyr(j)<styr) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" reset selparm_dev minyear to styr for selparm: "<<f<<endl;}
        if(selparm_dev_maxyr(j)>endyr) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" reset selparm_dev maxyear to endyr for selparm: "<<f<<endl;}
      }

      if(selparm_adjust_method==2 && selparm_envtype(f)==1)  //  now check on env links for illegal conditions
      {
        {N_warn++; warning<<" cannot use selparm_adjust_method==2 and multiplicative envlink for selex parameter; STOP "<<f<<endl; exit(1);}
      }
    }
 END_CALCS

//  ProgLabel_4.9.9  create arrays for the total set of selex parameters
  !!N_selparm2=N_selparm+N_selparm_env+N_selparm_blk+N_selparm_trend2;
  vector selparm_LO(1,N_selparm2)
  vector selparm_HI(1,N_selparm2)
  vector selparm_RD(1,N_selparm2)
  vector selparm_PR(1,N_selparm2)
  vector selparm_PRtype(1,N_selparm2)
  vector selparm_CV(1,N_selparm2)
  ivector selparm_PH(1,N_selparm2)

 LOCAL_CALCS
//  ProgLabel_4.9.10  special bound checking for size selex parameters
    z=0;  // parameter counter within this section
    for (f=1;f<=Ntypes;f++)
    {
      if(seltype(f,1)==8 || seltype(f,1)==22 || seltype(f,1)==23 || seltype(f,1)==24)
      {
        if(selparm_1(z+1,1)<len_bins_m(2))
        {N_warn++;
          warning<<"Fleet:_"<<f<<" min bound on parameter for size at peak is "<<selparm_1(z+1,1)<<"; should be >= midsize bin 2 ("<<len_bins_m(2)<<")"<<endl;}
        if(selparm_1(z+1,1)<len_bins_dat(1) && seltype(f,1)==24)
        {N_warn++;
          warning<<"Fleet:_"<<f<<" min bound on parameter for size at peak is "<<selparm_1(z+1,1)<<"; which is < min databin ("<<len_bins_dat(1)<<"), so illogical."<<endl;}
        if(selparm_1(z+1,2)>len_bins_m(nlength-1))
        {N_warn++;
          warning<<"Fleet:_"<<f<<" max bound on parameter for size at peak is "<<selparm_1(z+1,2)<<"; should be <= midsize bin N-1 ("<<len_bins_m(nlength-1)<<")"<<endl;}
      }
      z+=N_selparmvec(f);
    }
// end special bound checking

//  ProgLabel_4.9.11  create time/fleet array indicating when changes in selex occcur
  time_vary_sel.initialize();
  time_vary_makefishsel.initialize();
  time_vary_sel(styr)=1;
//  if(Do_Forecast>0) time_vary_sel(endyr+1)=1;
  time_vary_sel(endyr+1)=1;
  time_vary_makefishsel(styr)=1;
  time_vary_makefishsel(styr-3)=1;
//  if(Do_Forecast>0) time_vary_makefishsel(endyr+1)=1;
  time_vary_makefishsel(endyr+1)=1;
  for(y=styr+1;y<=endyr;y++)
  {
    z=0;  // parameter counter within this section
    for (f=1;f<=2*Ntypes;f++)
    {
      if(seltype(f,1)==5 || seltype(f,1)==15)   // mirror
      {
        if(f<=Ntypes) {time_vary_sel(y,f)=time_vary_sel(y,seltype(f,4));} else {time_vary_sel(y,f)=time_vary_sel(y,seltype(f,4)+Ntypes);}
        z+=seltype_Nparam(seltype(f,1));
      }
      else
      {
        if(seltype_Nparam(seltype(f,1))>0 || (seltype(f,2)==1) || (seltype(f,2)==2))      // type has parms, so look for adjustments
        {
          for(j=1;j<=N_selparmvec(f);j++)
          {
            z++;
            if(selparm_envuse(z)!=0)          // env linkage
            {
             if((env_data_RD(y,selparm_envuse(z))!=env_data_RD(y-1,selparm_envuse(z)) || selparm_envtype(z)==3 )) time_vary_sel(y,f)=1;
            }
            if(selparm_1(z,9)>=1)  // dev vector
            {
              s=selparm_1(z,11)+1;
              if(s>endyr) s=endyr;
              if(y>=selparm_1(z,10) && y<=s) time_vary_sel(y,f)=1;
            }

            if(selparm_1(z,13)>0) //   blocks
            {
              if(Block_Defs_Sel(z,y)!=Block_Defs_Sel(z,y-1) ) time_vary_sel(y,f)=1;
            }

            if(selparm_1(z,13)<0) //   trend
            {
              time_vary_sel(y,f)=1;
            }
          }
        }
      }
      if(f<=Ntypes && seltype(f,2)<0)  //  retention is being mirrored
      {
        k=-seltype(f,2);
        if(time_vary_sel(y,k)>0) time_vary_sel(y,f)=1;
      }
    }  // end type

//    time_vary_makefishsel(y)(1,Nfleet)=time_vary_sel(y)(1,Nfleet);  //  error, this will only do size selex
    for(f=1;f<=Nfleet;f++)
    {
      if(time_vary_sel(y,f)>0 || time_vary_sel(y,f+Ntypes)>0) time_vary_makefishsel(y,f)=1;
    }

    if(time_vary_MG(y,2)>0 || time_vary_MG(y,3)>0 || WTage_rd>0)
    {
      time_vary_makefishsel(y)=1;
    }
  } // end years

//  ProgLabel_4.9.12  create vectors, e.g. selparm_PH(), that will be used to create actual array of estimted parameters
   for (f=1;f<=N_selparm;f++)
   {
    selparm_LO(f)=selparm_1(f,1);
    selparm_HI(f)=selparm_1(f,2);
    selparm_RD(f)=selparm_1(f,3);
    selparm_PR(f)=selparm_1(f,4);
    selparm_PRtype(f)=selparm_1(f,5);
    selparm_CV(f)=selparm_1(f,6);
    selparm_PH(f)=selparm_1(f,7);
   }
   j=N_selparm;
   if(N_selparm_env>0)
   for (f=1;f<=N_selparm_env;f++)
   {
    j++;
    if(customenvsetup==0) k=1;
    else k=f;
    selparm_LO(j)=selparm_env_1(k,1);
    selparm_HI(j)=selparm_env_1(k,2);
    selparm_RD(j)=selparm_env_1(k,3);
    selparm_PR(j)=selparm_env_1(k,4);
    selparm_PRtype(j)=selparm_env_1(k,5);
    selparm_CV(j)=selparm_env_1(k,6);
    selparm_PH(j)=selparm_env_1(k,7);
   }

   if(N_selparm_blk>0)
   for (f=1;f<=N_selparm_blk;f++)
   {
    j++;
    if(customblocksetup==0) k=1;
    else k=f;
    selparm_LO(j)=selparm_blk_1(k,1);
    selparm_HI(j)=selparm_blk_1(k,2);
    selparm_RD(j)=selparm_blk_1(k,3);
    selparm_PR(j)=selparm_blk_1(k,4);
    selparm_PRtype(j)=selparm_blk_1(k,5);
    selparm_CV(j)=selparm_blk_1(k,6);
    selparm_PH(j)=selparm_blk_1(k,7);
   }

   if(N_selparm_trend>0)
   for (f=1;f<=N_selparm_trend2;f++)
   {
    j++;
    selparm_LO(j)=selparm_trend_1(f,1);
    selparm_HI(j)=selparm_trend_1(f,2);
    selparm_RD(j)=selparm_trend_1(f,3);
    selparm_PR(j)=selparm_trend_1(f,4);
    selparm_PRtype(j)=selparm_trend_1(f,5);
    selparm_CV(j)=selparm_trend_1(f,6);
    selparm_PH(j)=selparm_trend_1(f,7);
   }

 END_CALCS


//  ProgLabel_4.10  Read tag recapture parameter setup
// if Tags are used, the read parameters for initial tag loss, chronic tag loss, andd
// fleet-specific tag reporting.  Of these, only reporting rate will be allowed to be time-varying
  init_int TG_custom;  // 1=read; 0=create default parameters
  !! echoinput<<TG_custom<<" TG_custom (need to read even if no tag data )"<<endl;
  !! k=TG_custom*Do_TG*(3*N_TG+2*Nfleet);
  init_matrix TG_parm1(1,k,1,14);  // read initial values
  !! if(k>0) echoinput<<" Tag parameters as read "<<endl<<TG_parm1<<endl;
  !! k=Do_TG*(3*N_TG+2*Nfleet);
  matrix TG_parm2(1,k,1,14);
  !!if(Do_TG>0) {k1=k;} else {k1=1;}
  vector TG_parm_LO(1,k1);
  vector TG_parm_HI(1,k1);
  ivector TG_parm_PH(1,k1);
 LOCAL_CALCS
  if(Do_TG>0)
  {
    if(TG_custom==1)
    {
      TG_parm2=TG_parm1;  // assign to the read values
    }
    else
    {
      TG_parm2.initialize();
      onenum="    ";
      for(j=1;j<=N_TG;j++)
      {
        TG_parm2(j,1)=-10;  // min
        TG_parm2(j,2)=10;   // max
        TG_parm2(j,3)=-9.;   // init
        TG_parm2(j,4)=-9.;   // prior
        TG_parm2(j,5)=1.;   // default prior type is symmetric beta
        TG_parm2(j,6)=0.001;  //  prior is quite diffuse
        TG_parm2(j,7)=-4;  // phase
      }
      for(j=1;j<=N_TG;j++)
      {
        TG_parm2(j+N_TG)=TG_parm2(1);  // set chronic tag retention equal to initial tag_retention
      }
      for(j=1;j<=N_TG;j++)  // set overdispersion
      {
        TG_parm2(j+2*N_TG,1)=1;  // min
        TG_parm2(j+2*N_TG,2)=10;   // max
        TG_parm2(j+2*N_TG,3)=2.;   // init
        TG_parm2(j+2*N_TG,4)=2.;   // prior
        TG_parm2(j+2*N_TG,5)=1.;   // default prior type is symmetric beta
        TG_parm2(j+2*N_TG,6)=0.001;  //  prior is quite diffuse
        TG_parm2(j+2*N_TG,7)=-4;  // phase
      }
      for(j=1;j<=Nfleet;j++)
      {
        TG_parm2(j+3*N_TG)=TG_parm2(1);  // set tag reporting equal to near 1.0, as is the tag retention parameters
      }
      // set tag reporting decay to nil decay rate
      for(j=1;j<=Nfleet;j++)
      {
        k=j+3*N_TG+Nfleet;
        TG_parm2(k,1)=-4.;
        TG_parm2(k,2)=0.;
        TG_parm2(k,3)=0.;
        TG_parm2(k,4)=0.;    // prior of zero
        TG_parm2(k,5)=0.;  // default prior is squared dev
        TG_parm2(k,6)=2.;  // sd dev of prior
        TG_parm2(k,7)=-4.;
      }
    }

//  ProgLabel_4.10.1  create parameter count and parameter names for tag parameters
       onenum="    ";
       for(j=1;j<=N_TG;j++)
       {
       sprintf(onenum, "%d", j);
       ParCount++; ParmLabel+="TG_loss_init_"+onenum+CRLF(1);
      }
       for(j=1;j<=N_TG;j++)
      {
       sprintf(onenum, "%d", j);
       ParCount++; ParmLabel+="TG_loss_chronic_"+onenum+CRLF(1);
      }
       for(j=1;j<=N_TG;j++)
      {
       sprintf(onenum, "%d", j);
       ParCount++; ParmLabel+="TG_overdispersion_"+onenum+CRLF(1);
      }
       for(j=1;j<=Nfleet;j++)
      {
       sprintf(onenum, "%d", j);
       ParCount++; ParmLabel+="TG_report_fleet:_"+onenum+CRLF(1);
      }
       for(j=1;j<=Nfleet;j++)
      {
       sprintf(onenum, "%d", j);
       ParCount++; ParmLabel+="TG_rpt_decay_fleet:_"+onenum+CRLF(1);
      }

    TG_parm_LO=column(TG_parm2,1);
    TG_parm_HI=column(TG_parm2,2);
    k=3*N_TG+2*Nfleet;
    for(j=1;j<=k;j++) TG_parm_PH(j)=TG_parm2(j,7);  // write it out due to no typecast available
    echoinput<<" Processed/generated Tag parameters "<<endl<<TG_parm2<<endl;

  }
  else
  {
    TG_parm_LO.initialize();
    TG_parm_HI.initialize();
    TG_parm_PH.initialize();
  }
 END_CALCS


//  ProgLabel_4.11  read variance adjustment and various variance related inputs
  init_int Do_Var_adjust
  init_matrix var_adjust1(1,6*Do_Var_adjust,1,Ntypes)
  matrix var_adjust(1,6,1,Ntypes)
 LOCAL_CALCS
  echoinput<<Do_Var_adjust<<" Do_Var_adjust "<<endl;
  if(Do_Var_adjust>0)
  {
    var_adjust=var_adjust1;
   echoinput<<" Varadjustments as read "<<endl<<var_adjust1<<endl;
  }
  else
  {
    var_adjust(1)=0.;
    var_adjust(2)=0.;
    var_adjust(3)=0.;
    var_adjust(4)=1.;
    var_adjust(5)=1.;
    var_adjust(6)=1.;
  }
 END_CALCS

  init_number max_lambda_phase
  init_number sd_offset

 LOCAL_CALCS
  echoinput<<max_lambda_phase<<" max_lambda_phase "<<endl;
  echoinput<<sd_offset<<" sd_offset (adds log(s)) "<<endl;
  if(sd_offset==0)
  {
    N_warn++; warning<<" With sd_offset set to 0, be sure you are not estimating any variance parameters "<<endl;
  }
  if(depletion_fleet>0 && max_lambda_phase<2)
    {
      max_lambda_phase=2;
      N_warn++; warning<<"Increase max_lambda_phase to 2 because depletion fleet is being used"<<endl;
    }
 END_CALCS

//  ProgLabel_4.11.1  define type/phase arrays for lambdas
  matrix surv_lambda(1,Ntypes,1,max_lambda_phase)
  matrix disc_lambda(1,Ntypes,1,max_lambda_phase)
  matrix mnwt_lambda(1,Ntypes,1,max_lambda_phase)
  matrix length_lambda(1,Ntypes,1,max_lambda_phase)
  matrix age_lambda(1,Ntypes,1,max_lambda_phase)
  matrix sizeage_lambda(1,Ntypes,1,max_lambda_phase)
  vector init_equ_lambda(1,max_lambda_phase)
  matrix catch_lambda(1,Ntypes,1,max_lambda_phase)
  vector recrdev_lambda(1,max_lambda_phase)
  vector parm_prior_lambda(1,max_lambda_phase)
  vector parm_dev_lambda(1,max_lambda_phase)
  vector CrashPen_lambda(1,max_lambda_phase)
  vector Morphcomp_lambda(1,max_lambda_phase)
  matrix SzFreq_lambda(1,SzFreq_N_Like,1,max_lambda_phase)
  matrix TG_lambda1(1,N_TG2,1,max_lambda_phase)
  matrix TG_lambda2(1,N_TG2,1,max_lambda_phase)
  vector F_ballpark_lambda(1,max_lambda_phase)
  
//  ProgLabel_4.11.2  read and process any lambda adjustments
  init_int N_lambda_changes
  init_matrix Lambda_changes(1,N_lambda_changes,1,5)
  int N_changed_lambdas
 LOCAL_CALCS
   echoinput<<N_lambda_changes<<" N lambda changes "<<endl;
   if(N_lambda_changes>0) echoinput<<" lambda changes "<<endl<<Lambda_changes<<endl;
   surv_lambda=1.;  // 1
   disc_lambda=1.;  // 2
   mnwt_lambda=1.;  // 3
   length_lambda=1.; // 4
   age_lambda=1.;  // 5
   SzFreq_lambda=1.;  // 6
   sizeage_lambda=1.; // 7
   catch_lambda=1.; // 8
   init_equ_lambda=1.; // 9
   recrdev_lambda=1.; // 10
   parm_prior_lambda=1.; // 11
   parm_dev_lambda=1.; // 12
   CrashPen_lambda=1.; // 13
   Morphcomp_lambda=1.; // 14
   TG_lambda1=1.; // 15
   TG_lambda2=1.;  //16
   F_ballpark_lambda=1.;  // 17

    if(depletion_fleet>0)
    {
      for(f=1;f<=Ntypes;f++)
      {
        surv_lambda(f,1)=0.0;
        disc_lambda(f,1)=0.0;
        mnwt_lambda(f,1)=0.0;
        length_lambda(f,1)=0.0;
        age_lambda(f,1)=0.0;
        sizeage_lambda(f,1)=0.0;
//        catch_lambda(f,1)=0.0;  //  keep this positive to prevent crashes from bad fit to catch
      }
      if(SzFreq_Nmeth>0)
      {
        for(z=1;z<=SzFreq_N_Like;z++)
        {SzFreq_lambda(z,1)=0.0;}
      }
      if(N_TG2>0)
      {
        for(z=1;z<=N_TG2;z++)
        {
          TG_lambda1(z,1)=0.0;
          TG_lambda2(z,1)=0.0;
        }
      }
      init_equ_lambda(1)=0.0;
      recrdev_lambda(1)=0.0;
      Morphcomp_lambda(1)=0.0;
      F_ballpark_lambda(1)=0.0;

      surv_lambda(depletion_fleet,1)=1.0;
    }

    N_changed_lambdas=0;
    for(j=1;j<=N_lambda_changes;j++)
    {
      k=Lambda_changes(j,1);  // like component
      f=Lambda_changes(j,2);  // fleet
      s=Lambda_changes(j,3);  // phase
      if(k<=14)
      {
        if(f>Ntypes)
        {
          k=0;
          N_warn++;
          warning<<" illegal fleet/survey for lambda change at row: "<<j<<" fleet: "<<f<<" > Ntypes"<<endl;
        }
      }
      else if(k<=16)  // tag data
      {
        if(f>N_TG2)
        {
          k=0;
          N_warn++;
          warning<<" illegal tag group for lambda change at row: "<<j<<" Tag: "<<f<<" > N_taggroups"<<endl;
        }
      }
      else if(k>17)
      {
        k=0;
        N_warn++;
        warning<<" illegal lambda_type for lambda change at row: "<<j<<" Method: "<<k<<" > 17"<<endl;
      }
      if(s>max_lambda_phase)
      {k=0; N_warn++;  warning<<" illegal request for lambda change at row: "<<j<<" phase: "<<s<<" > max_lam_phase: "<<max_lambda_phase<<endl;}
//      if(s>Turn_off_phase) s=max(1,Turn_off_phase);
      temp=Lambda_changes(j,4);  // value
      if(temp!=0.0 && temp!=1.0) N_changed_lambdas++;
      z=Lambda_changes(j,5);   // special for sizefreq
      switch(k)
      {
        case 0:  // do nothing
        {break;}
        case 1:  // survey
          {surv_lambda(f)(s,max_lambda_phase)=temp;  break;}
        case 2:  // discard
          {disc_lambda(f)(s,max_lambda_phase)=temp;  break;}
        case 3:  // meanbodywt
          {mnwt_lambda(f)(s,max_lambda_phase)=temp; break;}
        case 4:  // lengthcomp
          {length_lambda(f)(s,max_lambda_phase)=temp; break;}
        case 5:  // agecomp
        {age_lambda(f)(s,max_lambda_phase)=temp; break;}
        case 6:  // sizefreq comp
        {
          z=Lambda_changes(j,5);  //  sizefreq method
          if(z>SzFreq_Nmeth) {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" reading sizefreq lambda change for method > Nmeth "<<Lambda_changes(j,5)<<endl; exit(1);}
          SzFreq_lambda(SzFreq_LikeComponent(f,z))(s,max_lambda_phase) = temp;
          break;
        }
        case 7:  // size-at-age
          {sizeage_lambda(f)(s,max_lambda_phase)=temp; break;}
        case 8:  // catch
          {catch_lambda(f)(s,max_lambda_phase)=temp; break;}
        case 9:  // init_equ_catch
          {init_equ_lambda(s,max_lambda_phase)=temp; break;}
        case 10:  // recr_dev
          {recrdev_lambda(s,max_lambda_phase)=temp; break;}
        case 11:  // parm_prior
          {parm_prior_lambda(s,max_lambda_phase)=temp; break;}
        case 12:  // parm_dev
          {parm_dev_lambda(s,max_lambda_phase)=temp; break;}
        case 13:  // crash_penalty
          {CrashPen_lambda(s,max_lambda_phase)=temp; break;}
        case 14:  // morphcomp
          {Morphcomp_lambda(s,max_lambda_phase)=temp; break;}
        case 15:  // Tag - multinomial by fleet  where f is now tag group
          {TG_lambda1(f)(s,max_lambda_phase)=temp; break;}
        case 16:  // Tag - total by time where f is now tag group
          {TG_lambda2(f)(s,max_lambda_phase)=temp; break;}
        case 17:  // F ballpark
          {F_ballpark_lambda(s,max_lambda_phase)=temp; break;}
      }
    }
    for(f=1;f<=Ntypes;f++)
    {
      if(nyr_cr(f)==0) surv_lambda(f)=0.;
      if(nyr_disc(f)==0) disc_lambda(f)=0.;
      if(Nobs_l(f)==0) length_lambda(f)=0.;
      if(Nobs_a(f)==0) age_lambda(f)=0.;
      if(Nobs_ms(f)==0) sizeage_lambda(f)=0.;
    }
    if(nobs_mnwt==0) mnwt_lambda=0.;  //  more complicated to turn off for each fleet
 END_CALCS

//  ProgLabel_4.12  read setup for more derived quantities to include in the STD report
  init_int Do_More_Std
  init_ivector More_Std_Input(1,Do_More_Std*9)
 LOCAL_CALCS
  echoinput<<Do_More_Std<<" # read specs for more stddev reporting "<<endl;
  if(Do_More_Std>0)
  {echoinput<<More_Std_Input<<" # vector with selex type, len/age, year, N selex bins, Growth pattern, N growth ages, N_at_age_Area, NatAge_yr, Natage_ages"<<endl;}
  else
  {echoinput<<" # placeholder vector with selex type, len/age, year, N selex bins, Growth pattern, N growth ages"<<endl;}
 END_CALCS

  int Do_Selex_Std;
  int Selex_Std_AL;
  int Selex_Std_Year;
  int Selex_Std_Cnt;
  int Do_Growth_Std;
  int Growth_Std_Cnt;
  int Do_NatAge_Std;
  int NatAge_Std_Year;
  int NatAge_Std_Cnt;
  int Extra_Std_N;   //  dimension for the sdreport vector Selex_Std which also contains the Growth_Std

 LOCAL_CALCS
   if(Do_More_Std==1)
   {
     Do_Selex_Std=More_Std_Input(1);
     Selex_Std_AL=More_Std_Input(2);
     Selex_Std_Year=More_Std_Input(3);
     if(Selex_Std_Year<0) Selex_Std_Year=endyr;
     Selex_Std_Cnt=More_Std_Input(4);
     Do_Growth_Std=More_Std_Input(5);
     if(MG_active(2)==0) Do_Growth_Std=0;
     Growth_Std_Cnt=More_Std_Input(6);
     Do_NatAge_Std=More_Std_Input(7);
     NatAge_Std_Year=More_Std_Input(8);
     if(NatAge_Std_Year<0) NatAge_Std_Year=endyr+1;
     NatAge_Std_Cnt=More_Std_Input(9);
   }
   else
   {
     Do_Selex_Std=0;
     Selex_Std_AL=1;
     Selex_Std_Year=endyr;
     Selex_Std_Cnt=0;
     Do_Growth_Std=0;
     Growth_Std_Cnt=0;
     Do_NatAge_Std=0;
     NatAge_Std_Cnt=0;
     NatAge_Std_Year=endyr;
   }
 END_CALCS

  init_ivector Selex_Std_Pick(1,Selex_Std_Cnt);
  init_ivector Growth_Std_Pick(1,Growth_Std_Cnt);
  init_ivector NatAge_Std_Pick(1,NatAge_Std_Cnt);

 LOCAL_CALCS
  if(Selex_Std_Cnt>0) echoinput<<Selex_Std_Pick<<" # vector with selex std bin picks (-1 in first bin to self-generate)"<<endl;
  if(Growth_Std_Cnt>0) echoinput<<Growth_Std_Pick<<" # vector with growth std bin picks (-1 in first bin to self-generate)"<<endl;
  if(NatAge_Std_Cnt>0) echoinput<<NatAge_Std_Pick<<" # vector with NatAge std bin picks (-1 in first bin to self-generate)"<<endl;

// reset the counter here after using it to dimension the input statement above
  if(Do_Selex_Std<=0) Selex_Std_Cnt=0;
  if(Do_Growth_Std<=0) Growth_Std_Cnt=0;
  if(Do_NatAge_Std==0) NatAge_Std_Cnt=0;

  Extra_Std_N=0;
  if(Do_Selex_Std>0)
  {
    if(Selex_Std_Pick(1)<=0)  //  then self-generate even bin selection
    {
      if(Selex_Std_AL==1)
      {
        j=nlength/(Selex_Std_Cnt-1);
        Selex_Std_Pick(1)=j/2;
        for(i=2;i<=Selex_Std_Cnt-1;i++) Selex_Std_Pick(i)=Selex_Std_Pick(i-1)+j;
        Selex_Std_Pick(Selex_Std_Cnt)=nlength;
      }
      else
      {
        j=nages/(Selex_Std_Cnt-1);
        Selex_Std_Pick(1)=j/2;
        for(i=2;i<=Selex_Std_Cnt-1;i++) Selex_Std_Pick(i)=Selex_Std_Pick(i-1)+j;
        Selex_Std_Pick(Selex_Std_Cnt)=nages;
      }
    }
    Extra_Std_N=gender*Selex_Std_Cnt;
  }

  if(Do_Growth_Std>0)
  {
    if(Growth_Std_Pick(1)<=0)
    {
      Growth_Std_Pick(1)=AFIX;
      Growth_Std_Pick(Growth_Std_Cnt)=nages;
      if(Growth_Std_Cnt>2)
      {
        k=Growth_Std_Cnt/2;
        for(i=2;i<=k;i++) Growth_Std_Pick(i)=Growth_Std_Pick(i-1)+1;
        j=(nages-Growth_Std_Pick(k))/(Growth_Std_Cnt-k);
        for(i=k+1;i<=Growth_Std_Cnt-1;i++) Growth_Std_Pick(i)=Growth_Std_Pick(i-1)+j;
      }
    }
  }
  Extra_Std_N+=gender*Growth_Std_Cnt;

  if(Do_NatAge_Std!=0)
  {
    if(NatAge_Std_Pick(1)<=0)
    {
      NatAge_Std_Pick(1)=1;
      NatAge_Std_Pick(NatAge_Std_Cnt)=nages;
      if(NatAge_Std_Cnt>2)
      {
        k=NatAge_Std_Cnt/2;
        for(i=2;i<=k;i++) NatAge_Std_Pick(i)=NatAge_Std_Pick(i-1)+1;
        j=(nages-NatAge_Std_Pick(k))/(NatAge_Std_Cnt-k);
        for(i=k+1;i<=NatAge_Std_Cnt-1;i++) NatAge_Std_Pick(i)=NatAge_Std_Pick(i-1)+j;
      }
    }
  }
  Extra_Std_N+=gender*NatAge_Std_Cnt;

  if(Extra_Std_N==0) Extra_Std_N=1;   //  assign a minimum length to dimension the sdreport vector Selex_Std
  echoinput<<"After processing"<<endl;
  if(Selex_Std_Cnt>0) echoinput<<Selex_Std_Pick<<" # vector with selex std bin picks (-1 in first bin to self-generate)"<<endl;
  if(Growth_Std_Cnt>0) echoinput<<Growth_Std_Pick<<" # vector with growth std bin picks (-1 in first bin to self-generate)"<<endl;
  if(NatAge_Std_Cnt>0) echoinput<<NatAge_Std_Pick<<" # vector with NatAge std bin picks (-1 in first bin to self-generate)"<<endl;

 END_CALCS

//  ProgLabel_4.13  end of reading from control file
  init_int fim // end of file indicator

 LOCAL_CALCS
  cout<<"If you see 999, we got to the end of the control file successfully! "<<fim<<endl;
  echoinput<<fim<<"  If you see 999, we got to the end of the control file successfully! "<<endl<<endl<<"***********************************"<<endl;
  if(fim!=999) abort();
 END_CALCS

//  ProgLabel_4.14  create count of active parameters and derived quantities
  int CoVar_Count;
  int active_count;    // count the active parameters
  int active_parms;    // count the active parameters
 LOCAL_CALCS
  if(Do_Benchmark>0)
  {
    N_STD_Mgmt_Quant=16;
  }
  else
  {N_STD_Mgmt_Quant=1;}
  Fcast_catch_start=N_STD_Mgmt_Quant;
  if(Do_Forecast>0) {N_STD_Mgmt_Quant+=N_Fcast_Yrs*(1+Do_Retain)+N_Fcast_Yrs;}
  k=ParCount+2*N_STD_Yr+N_STD_Yr_Dep+N_STD_Yr_Ofish+N_STD_Yr_F+N_STD_Mgmt_Quant+gender*Selex_Std_Cnt+gender*Growth_Std_Cnt;
  echoinput<<"N parameters: "<<ParCount<<endl<<"Parameters plus derived quant: "<<k<<endl;
 END_CALCS
  ivector active_parm(1,k)  //  pointer from active list to the element of the full parameter list to get label later


//***********************************************
//  ProgLabel_4.14.1  adjust the phases to negative if beyond turn_off_phase and find resultant max_phase
  int max_phase;

 LOCAL_CALCS
  echoinput<<"Adjust the phases "<<endl;
  max_phase=1;
  active_count=0;
  active_parm(1,ParCount)=0;
  ParCount=0;

  j=MGparm_PH.indexmax();

  for(k=1;k<=j;k++)
  {
    ParCount++;
    if(MGparm_PH(k)==-9999) {MGparm_RD(k)=prof_var(prof_var_cnt); prof_var_cnt+=1;}
    if(depletion_fleet>0 && MGparm_PH(k)>0) MGparm_PH(k)++;  //  add 1 to phase if using depletion fleet
    if(MGparm_PH(k) > Turn_off_phase) MGparm_PH(k) =-1;
    if(MGparm_PH(k) > max_phase) max_phase=MGparm_PH(k);
    if(MGparm_PH(k)>=0)
    {
      active_count++; active_parm(active_count)=ParCount;
    }
  }

  if(depletion_fleet>0 && MGparm_dev_PH>0) MGparm_dev_PH++;  //  add 1 to phase if using depletion fleet
  if(MGparm_dev_PH>Turn_off_phase) MGparm_dev_PH =-1;
  if(MGparm_dev_PH>max_phase) max_phase=MGparm_dev_PH;
  for(k=1;k<=N_MGparm_dev_tot;k++)
  {
    ParCount++;
    if(MGparm_dev_PH>=0)
    {
    active_count++; active_parm(active_count)=ParCount;
    }
  }

  for(j=1;j<=SRvec_PH.indexmax();j++)
  {
    ParCount++;
    if(SRvec_PH(j)==-9999) {SR_parm_1(j,3)=prof_var(prof_var_cnt); prof_var_cnt+=1;}
    if(depletion_fleet>0 && SRvec_PH(j)>0) SRvec_PH(j)++;  //  add 1 to phase if using depletion fleet
    if(depletion_fleet>0 && j==1) SRvec_PH(1)=1;  //
    if(SRvec_PH(j) > Turn_off_phase) SRvec_PH(j) =-1;
    if(SRvec_PH(j) > max_phase) max_phase=SRvec_PH(j);
    if(SRvec_PH(j)>=0)
    {
      active_count++; active_parm(active_count)=ParCount;
    }
  }

  if(recdev_cycle>0)
  {
    for(y=1;y<=recdev_cycle;y++)
    {
      ParCount++;
      recdev_cycle_LO(y)=recdev_cycle_parm_RD(y,1);
      recdev_cycle_HI(y)=recdev_cycle_parm_RD(y,2);
      recdev_cycle_PH(y)=recdev_cycle_parm_RD(y,7);
      if(depletion_fleet>0 && recdev_cycle_PH(y)>0) recdev_cycle_PH(y)++;  //  add 1 to phase if using depletion fleet
      if(recdev_cycle_PH(y) > Turn_off_phase) recdev_cycle_PH(y) =-1;
      if(recdev_cycle_PH(y) > max_phase) max_phase=recdev_cycle_PH(y);
      if(recdev_cycle_PH(y)>=0) {active_count++; active_parm(active_count)=ParCount;}
    }
  }

  if(depletion_fleet>0 && recdev_early_PH>0) recdev_early_PH++;  //  add 1 to phase if using depletion fleet
  if(recdev_early_PH > Turn_off_phase) recdev_early_PH =-1;
  if(recdev_early_PH > max_phase) max_phase=recdev_early_PH;

  if(recdev_do_early>0)
  {
  for(y=recdev_early_start;y<=recdev_early_end;y++)
  {
    ParCount++;
    if(recdev_early_PH>=0) {active_count++; active_parm(active_count)=ParCount;}
  }
  }

  if(depletion_fleet>0 && recdev_PH>0) recdev_PH++;  //  add 1 to phase if using depletion fleet
  if(recdev_PH > Turn_off_phase) recdev_PH =-1;
  if(recdev_PH > max_phase) max_phase=recdev_PH;
  if(do_recdev>0)
  {
  for(y=recdev_start;y<=recdev_end;y++)
  {
    ParCount++;
    if(recdev_PH>=0) {active_count++; active_parm(active_count)=ParCount;}
  }
  }

  Fcast_recr_PH2=max_phase+1;
  if(Do_Forecast>0)
  {
    if(Turn_off_phase>0)
    {
      if(Fcast_recr_PH!=0)  // read value for forecast_PH
      {
        Fcast_recr_PH2=Fcast_recr_PH;
        if(depletion_fleet>0 && Fcast_recr_PH2>0) Fcast_recr_PH2++;
        if(Fcast_recr_PH2 > Turn_off_phase) Fcast_recr_PH2 =-1;
        if(Fcast_recr_PH2 > max_phase) max_phase=Fcast_recr_PH2;
      }
      for(y=recdev_end+1;y<=YrMax;y++)
      {
        ParCount++;
        if(Fcast_recr_PH2>-1) {active_count++; active_parm(active_count)=ParCount;}
      }
    }
    else
      {
        Fcast_recr_PH2=-1;
      }

    for(y=endyr+1;y<=YrMax;y++)
    {
      ParCount++;
      if(Do_Impl_Error>0 && Fcast_recr_PH2>-1)
      {active_count++; active_parm(active_count)=ParCount;}
    }
  }
  else
  {Fcast_recr_PH2=-1;}

  for(f=1;f<=Nfleet;f++)
  {
    ParCount++;
    if(init_F_PH(f)==-9999) {init_F_parm_1(f,3)=prof_var(prof_var_cnt); init_F_RD(f)=init_F_parm_1(f,3);  prof_var_cnt++;}
    if(depletion_fleet>0 && init_F_PH(f)>0) init_F_PH(f)++;
    if(init_F_PH(f) > Turn_off_phase) init_F_PH(f) =-1;
    if(init_F_PH(f) > max_phase) max_phase=init_F_PH(f);
    if(init_F_PH(f)>=0)
    {
      active_count++; active_parm(active_count)=ParCount;
    }
  }

  if(F_Method==2)
  {
    for(g=1;g<=Nfleet*(TimeMax-styr+1);g++)
    {
      ParCount++;
      if(depletion_fleet>0 && have_catch(g)>0) have_catch(g)++;
      if(have_catch(g) > Turn_off_phase) have_catch(g) =-1;
      if(have_catch(g) > max_phase) max_phase=have_catch(g);
      if(have_catch(g)>=0)
      {
        active_count++; active_parm(active_count)=ParCount;
      }
    }
  }

  for(f=1;f<=Q_Npar;f++)
  {
    ParCount++;
    if(Q_parm_PH(f)==-9999) {Q_parm_1(f,3)=prof_var(prof_var_cnt); prof_var_cnt++;}
    if(depletion_fleet>0 && Q_parm_PH(f)>0) Q_parm_PH(f)++;
    if(Q_parm_PH(f) > Turn_off_phase) Q_parm_PH(f) =-1;
    if(Q_parm_PH(f) > max_phase) max_phase=Q_parm_PH(f);
    if(Q_parm_PH(f)>=0)
    {
      active_count++; active_parm(active_count)=ParCount;
    }
  }

//  ProgLabel_4.14.1.1  Auto-generate cubic spline setup while inside this parameter counting loop
  Ip=0;
  int N_knots;
  for (f=1;f<=2*Ntypes;f++)   //  check for cubic spline setup
  {
    if(f<=Ntypes)
    {fs=f;}
    else
    {fs=f-Ntypes;}
    if(seltype(f,1)==27)  //  reset the cubic spline knots for size or age comp
    {
      k=int(selparm_RD(Ip+1));  // setup method
      N_knots=seltype(f,4);  //  number of knots

      if(k==0)
      {}  //  do nothing
      else if(k==1 || k==2)  //  get new knots according to cumulative distribution of data
      {
        echoinput<<"Adjust the cubic_spline setup for fleet: "<<f<<endl;
        s=4;  // counter for which knot is being set
        z=1;  //  counter for  bins in cumulative distribution
        if(N_knots>=3)
        {
          temp=0.025;
          temp1=0.950/float(N_knots-1);  //  increment
        }
        else
        {
          N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"must have at least 3 knots in spline "<<endl;  exit(1);
        }
        if(f<=Ntypes)  // doing size Selex
        {
          while(temp<=0.975001)
          {
            while(obs_l_all(2,f,z)<temp)
            {
              z++;
            }
            //  intermediate knots are calculated from data_length_bins
            if(z>1)
            {selparm_RD(Ip+s)=len_bins_dat(z-1)+(temp-obs_l_all(2,f,z-1))/(obs_l_all(2,f,z)-obs_l_all(2,f,z-1))*(len_bins_dat(z)-len_bins_dat(z-1));}
            else
            {selparm_RD(Ip+s)=len_bins_dat(z);}
            s++;
            temp+=temp1;
          }
          echoinput<<"len_bins_dat: "<<len_bins_dat<<endl<<"Cum_comp: "<<obs_l_all(2,fs)(1,nlen_bin)<<endl<<"Knots: "<<selparm_RD(Ip+3+1,Ip+3+N_knots)<<endl;
        }
        else  //  age selex
        {
          while(temp<=0.975001)
          {
            while(obs_a_all(2,fs,z)<temp)
            {
              z++;
            }
            //  intermediate knots are calculated from data_length_bins
            if(z>1)
            {selparm_RD(Ip+s)=age_bins(z-1)+(temp-obs_a_all(2,fs,z-1))/(obs_a_all(2,fs,z)-obs_a_all(2,fs,z-1))*(age_bins(z)-age_bins(z-1));}
            else
            {selparm_RD(Ip+s)=age_bins(z);}
            s++;
            temp+=temp1;
          }
          echoinput<<"age_bins: "<<age_bins<<endl<<"Cum_comp: "<<obs_a_all(2,fs)(1,n_abins)<<endl<<"Knots: "<<selparm_RD(Ip+3+1,Ip+3+N_knots)<<endl;
        }
        if(k==2)  //  create default bounds, priors, etc.
        {
        echoinput<<"Do complete setup of lo, hi, prior, etc."<<endl;
          for(z=Ip+4;z<=Ip+3+N_knots;z++)
          {
            selparm_LO(z)=age_bins(1);
            selparm_HI(z)=age_bins(n_abins);
            selparm_PR(z)=0.;
            selparm_PRtype(z)=-1;
            selparm_CV(z)=0.;
            selparm_PH(z)=-99;
          }

          if(N_knots==3)
          {p=8;}
          else if (N_knots==4)
          {p=10;}
          else
          {p=3+N_knots+1+0.5*N_knots;}

          echoinput<<" p "<<p<<endl;
          for(z=N_knots+1+3;z<=3+2*N_knots;z++)
          {
            a=Ip+z;
            if(z<=p)
            {selparm_RD(a)=-5. + float(z-(N_knots+4))/float(p-(N_knots+4))*4.;}
            else
            {selparm_RD(a)=0.0;}
            selparm_LO(a)=-9.;
            selparm_HI(a)=7.;
            selparm_PR(a)=0.;
            selparm_PRtype(a)=1;
            selparm_CV(a)=0.001;
            selparm_PH(a)=2;
          }
          selparm_PH(Ip+p)=-99;
          selparm_PRtype(Ip+p)=-1;
          selparm_CV(Ip+p)=0.;

          p=Ip+1;
          selparm_LO(p)=0.;
          selparm_HI(p)=2.;
          selparm_PR(p)=0.;
          selparm_PRtype(p)=-1;
          selparm_CV(p)=0.;
          selparm_PH(p)=-99;
          p++;
          selparm_LO(p)=-0.001;
          selparm_HI(p)=1.;
          selparm_RD(p)=0.1;  // moderate positive gradient at bottom
          selparm_PR(p)=0.;
          selparm_PRtype(p)=1;  // SYMMETRIC BETA
          selparm_CV(p)=0.001;
          selparm_PH(p)=3;
          p++;
          selparm_LO(p)=-1.;
          selparm_HI(p)=0.001;
          if(N_knots>=3)
          {
          selparm_RD(p)=-0.001;  // small negative gradient at top
          selparm_PR(p)=0.;
          selparm_PRtype(p)=1;
          selparm_CV(p)=0.001;
          selparm_PH(p)=3;
          }
          else
          {
          selparm_RD(p)=0.00;
          selparm_PR(p)=0.;
          selparm_PRtype(p)=-1;
          selparm_CV(p)=0.;
          selparm_PH(p)=-99;
          }

          for(z=Ip+1;z<=Ip+3+2*N_knots;z++)
          {
            selparm_1(z,1)=selparm_LO(z);
            selparm_1(z,2)=selparm_HI(z);
            selparm_1(z,3)=selparm_RD(z);
            selparm_1(z,4)=selparm_PR(z);
            selparm_1(z,5)=selparm_PRtype(z);
            selparm_1(z,6)=selparm_CV(z);
            selparm_1(z,7)=selparm_PH(z);
          }
        }
      }
    }
    Ip+=N_selparmvec(f);
  }

   for(k=1;k<=selparm_PH.indexmax();k++)
   {
     ParCount++;
     if(selparm_PH(k)==-9999) {selparm_RD(k)=prof_var(prof_var_cnt); prof_var_cnt++;}
     if(depletion_fleet>0 && selparm_PH(k)>0) selparm_PH(k)++;
     if(selparm_PH(k) > Turn_off_phase) selparm_PH(k) =-1;
     if(selparm_PH(k) > max_phase) max_phase=selparm_PH(k);
     if(selparm_PH(k)>=0)
    {
      active_count++; active_parm(active_count)=ParCount;
    }
   }

   if(depletion_fleet>0 && selparm_dev_PH>0) selparm_dev_PH++;
   if(selparm_dev_PH > Turn_off_phase) selparm_dev_PH =-1;
   if(selparm_dev_PH > max_phase) max_phase=selparm_dev_PH;
  for(k=1;k<=N_selparm_dev_tot;k++)
  {
    ParCount++;
    if(selparm_dev_PH>=0)
    {
    active_count++; active_parm(active_count)=ParCount;
    }
  }

  if(Do_TG>0)
  {
    for(k=1;k<=3*N_TG+2*Nfleet;k++)
    {
      ParCount++;
      if(depletion_fleet>0 && TG_parm_PH(k)>0) TG_parm_PH(k)++;
      if(TG_parm_PH(k) > Turn_off_phase) TG_parm_PH(k) =-1;
      if(TG_parm_PH(k) > max_phase) max_phase=TG_parm_PH(k);
      if(TG_parm_PH(k)>=0)
      {
      active_count++; active_parm(active_count)=ParCount;
      }
    }
  }
  
  if(Do_Forecast>0 && Turn_off_phase>0)
  {
    if(Fcast_recr_PH==0)  // read value for forecast_PH.  This code is repeats earlier code in case other parameters have changed maxphase
    {
      Fcast_recr_PH2=max_phase+1;
    }
  }

  echoinput<<"Active parameters: "<<active_count<<endl<<"Turn_off_phase "<<Turn_off_phase<<endl<<" max_phase "<<max_phase<<endl;
  if(Turn_off_phase<=0)
  {func_eval(1)=1;}
  else
  {
     func_conv(max_phase)=final_conv;  func_eval(max_phase)=10000;
     func_conv(max_phase+1)=final_conv;  func_eval(max_phase+1)=10000;
  }

//  ProgLabel_4.14.2  add count of derived quantities and create labels for these quantities
    j=ParCount;
    active_parms=active_count;
    CoVar_Count=active_count;
  onenum="    ";
  for(y=styr-2;y<=YrMax;y++)
  {
    if(STD_Yr_Reverse(y)>0)
    {
    CoVar_Count++; j++; active_parm(CoVar_Count)=j;
    if(y==styr-2)
    {ParmLabel+="SPB_Virgin";}
    else if(y==styr-1)
    {ParmLabel+="SPB_Initial";}
    else
    {
//      _itoa(y,onenum,10);
      sprintf(onenum, "%d", y);
      ParmLabel+="SPB_"+onenum+CRLF(1);
    }
    }
  }

  for(y=styr-2;y<=YrMax;y++)
  {
    if(STD_Yr_Reverse(y)>0)
    {
    CoVar_Count++; j++; active_parm(CoVar_Count)=j;
    if(y==styr-2)
    {ParmLabel+="Recr_Virgin";
      }
    else if(y==styr-1)
    {ParmLabel+="Recr_Initial";
      }
    else
    {
//      _itoa(y,onenum,10);
     sprintf(onenum, "%d", y);
      ParmLabel+="Recr_"+onenum+CRLF(1);
    }
  }
  }

  for(y=styr;y<=YrMax;y++)
  {
    if(STD_Yr_Reverse_Ofish(y)>0)
    {
      CoVar_Count++; j++; active_parm(CoVar_Count)=j;
//      _itoa(y,onenum,10);
      sprintf(onenum, "%d", y);
      ParmLabel+="SPRratio_"+onenum+CRLF(1);
    }
  }

  //F_std
  for(y=styr;y<=YrMax;y++)
  {
    if(STD_Yr_Reverse_F(y)>0)
    {
      CoVar_Count++; j++; active_parm(CoVar_Count)=j;
//      _itoa(y,onenum,10);
      sprintf(onenum, "%d", y);
      ParmLabel+="F_"+onenum+CRLF(1);
    }
  }

  for(y=styr;y<=YrMax;y++)
  {
    if(STD_Yr_Reverse_Dep(y)>0)
    {
      CoVar_Count++; j++; active_parm(CoVar_Count)=j;
//      _itoa(y,onenum,10);
    sprintf(onenum, "%d", y);
    ParmLabel+="Bratio_"+onenum+CRLF(1);
    }
  }

//  create labels for Mgmt_Quant
  if(Do_Benchmark>0)
    {
      ParmLabel+="SSB_Unfished"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="TotBio_Unfished"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="SmryBio_Unfished"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="Recr_Unfished"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="SSB_Btgt"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="SPR_Btgt"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="Fstd_Btgt"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="TotYield_Btgt"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="SSB_SPRtgt"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="Fstd_SPRtgt"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="TotYield_SPRtgt"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="SSB_MSY"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="SPR_MSY"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="Fstd_MSY"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="TotYield_MSY"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="RetYield_MSY"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
    }
    else
    {
      ParmLabel+="Bzero_again"+CRLF(1); CoVar_Count++; j++; active_parm(CoVar_Count)=j;
    }

    if(Do_Forecast>0)
    {
      for(y=endyr+1;y<=YrMax;y++)
      {
        CoVar_Count++; j++; active_parm(CoVar_Count)=j;
        sprintf(onenum, "%d", y);
        ParmLabel+="ForeCatch_"+onenum+CRLF(1);
      }
      for(y=endyr+1;y<=YrMax;y++)
      {
        CoVar_Count++; j++; active_parm(CoVar_Count)=j;
        sprintf(onenum, "%d", y);
        ParmLabel+="OFLCatch_"+onenum+CRLF(1);
      }
      if(Do_Retain==1)
      {
        for(y=endyr+1;y<=YrMax;y++)
        {
          CoVar_Count++; j++; active_parm(CoVar_Count)=j;
          sprintf(onenum, "%d", y);
          ParmLabel+="ForeCatchret_"+onenum+CRLF(1);
        }
      }
    }

// do labels for Selex_Std
    if(Do_Selex_Std>0)
    {
      for(g=1;g<=gender;g++)
      for(i=1;i<=Selex_Std_Cnt;i++)
      {
        CoVar_Count++; j++; active_parm(CoVar_Count)=j;
        if(Selex_Std_AL==1)
        {
          if(Selex_Std_Pick(i)>nlength)
          {
            N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" cannot select stdev for length bin greater than nlength "<<Selex_Std_Pick(i)<<" > "<<nlength<<endl; exit(1);
          }
          ParmLabel+="Selex_std_"+NumLbl(Do_Selex_Std)+"_"+GenderLbl(g)+"_L_"+NumLbl(len_bins(Selex_Std_Pick(i)))+CRLF(1);
        }
        else
        {
          if(Selex_Std_Pick(i)>nages)
          {
            N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" cannot select stdev for age bin greater than maxage "<<Selex_Std_Pick(i)<<" > "<<nages<<endl; exit(1);
          }
          ParmLabel+="Selex_std_"+NumLbl(Do_Selex_Std)+"_"+GenderLbl(g)+"_A_"+NumLbl(age_vector(Selex_Std_Pick(i)))+CRLF(1);
        }
      }
    }
    if(Do_Growth_Std>0)
    {
      for(g=1;g<=gender;g++)
      for(i=1;i<=Growth_Std_Cnt;i++)
      {
        CoVar_Count++; j++; active_parm(CoVar_Count)=j;
        ParmLabel+="Grow_std_"+NumLbl(Do_Growth_Std)+"_"+GenderLbl(g)+"_A_"+NumLbl(age_vector(Growth_Std_Pick(i)))+CRLF(1);
      }
    }
    if(Do_NatAge_Std!=0)
    {
      for(g=1;g<=gender;g++)
      for(i=1;i<=NatAge_Std_Cnt;i++)
      {
        CoVar_Count++; j++; active_parm(CoVar_Count)=j;
        if(Do_NatAge_Std>0)
        {ParmLabel+="NatAge_std_"+NumLbl(Do_NatAge_Std)+"_"+GenderLbl(g)+"_A_"+NumLbl(age_vector(NatAge_Std_Pick(i)))+CRLF(1);}
        else
        {ParmLabel+="NatAge_std_All_"+GenderLbl(g)+"_A_"+NumLbl(age_vector(NatAge_Std_Pick(i)))+CRLF(1);}
      }
    }
    if(Do_Selex_Std==0 && Do_Growth_Std==0 && Do_NatAge_Std==0)
    {
      CoVar_Count++; j++; active_parm(CoVar_Count)=j;
      ParmLabel+="Bzero_again"+CRLF(1);
    }

   sprintf(onenum, "%d", int(100*depletion_level));
   switch(depletion_basis)
    {
      case 0:
      {
        depletion_basis_label+="no_depletion_basis";
        break;
      }
      case 1:
      {
        depletion_basis_label+=" "+onenum+"%*Virgin_Biomass";
        break;
      }
      case 2:
      {
        depletion_basis_label+=" "+onenum+"%*B_MSY";
        break;
      }
      case 3:
      {
        depletion_basis_label+=" "+onenum+"%*StartYr_Biomass";
        break;
      }
    }

   switch (SPR_reporting)
  {
    case 0:      // keep as raw value
    {
      SPR_report_label+=" raw_SPR";
      break;
    }
    case 1:  // compare to SPR
    {
      sprintf(onenum, "%d", int(100.*SPR_target));
      SPR_report_label+=" (1-SPR)/(1-SPR_"+onenum+"%)";
      break;
    }
    case 2:  // compare to SPR_MSY
    {
      SPR_report_label+=" (1-SPR)/(1-SPR_MSY)";
      break;
    }
    case 3:  // compare to SPR_Btarget
    {
      sprintf(onenum, "%d", int(100.*BTGT_target));
      SPR_report_label+=" (1-SPR)/(1-SPR_at_B"+onenum+"%)";
      break;
    }
    case 4:
    {
      SPR_report_label+=" 1-SPR";
      break;
    }
  }

  switch (F_report_basis)
  {
    case 0:  // raw
    {
      F_report_label="_abs_F";
      break;
    }
    case 1:
    {
      sprintf(onenum, "%d", int(100.*SPR_target));
      F_report_label="(F)/(F"+onenum+"%SPR)";
      break;
    }
    case 2:
    {
      F_report_label="(F)/(Fmsy)";
      break;
    }
    case 3:
    {
      sprintf(onenum, "%d", int(100.*BTGT_target));
      F_report_label="(F)/(F_at_B"+onenum+"%)";
      break;
    }
  }

   switch (F_reporting)
  {
    case 0:      // keep as raw value
    {
      F_report_label+=";_no_F_report";
      break;
    }
    case 1:      // exploitation rate in biomass
    {
      F_report_label+=";_with_F=Exploit(bio)";
      break;
    }
    case 2:      // exploitation rate in numbers
    {
      F_report_label+=";_with_F=Exploit(num)";
      break;
    }
    case 3:      // sum of F mults
    {
      F_report_label+=";_with_F=sum(full_Fs)";
      break;
    }
    case 4:      // F=Z-M for specified ages
    {
      F_report_label+=";_with_F=Z-M;_for_ages_";
      sprintf(onenum, "%d", int(F_reporting_ages(1)));
      F_report_label+=onenum;
      sprintf(onenum, "%d", int(F_reporting_ages(2)));
      F_report_label+="_"+onenum;
      break;
    }
  }
  echoinput<<"Active parameters plus derived quantities:  "<<CoVar_Count<<endl;
 END_CALCS

  !!k=gmorph*(YrMax-styr+1);
//  ProgLabel_4.14.3  create matrix CoVar and set it to receive the covariance output
  matrix save_G_parm(1,k,1,22);
  matrix save_seas_parm(1,nseas,1,10);
  matrix CoVar(1,CoVar_Count,1,CoVar_Count+1);
  !!save_G_parm.initialize();
  !!CoVar.initialize();
  !!set_covariance_matrix(CoVar);

//  ***********************************************
//  ProgLabel_4.15  Read empirical wt-at-age
  int N_WTage_rd
  int N_WTage_maxage
  int y2
 LOCAL_CALCS
   if(WTage_rd>0)
   {
     ad_comm::change_datafile_name("wtatage.ss");
     k1=2;
   }
   else
   {
    k1=0;
    N_WTage_rd=0;
    N_WTage_maxage=0;
   }
 END_CALCS

  init_vector junkvec(1,k1)
 LOCAL_CALCS
  if(k1>0)
  {
    echoinput<<"WT-at-age input"<<junkvec<<endl;
    N_WTage_rd=junkvec(1);
    N_WTage_maxage=junkvec(2);
    k2=TimeMax_Fcast_std+1;
  }
  else
  {
    k2=styr;
    N_WTage_maxage=nages;
  }

 END_CALCS
  init_matrix WTage_in(1,N_WTage_rd,1,7+N_WTage_maxage)
  vector junkvec2(0,nages)
  4darray WTage_emp(styr-3*nseas,k2,1,gender*N_GP*nseas,-2,Ntypes,0,nages)  //  set to begin period for pop (type=0), or mid period for fleet/survey
// read:  yr, seas, gender, growpattern, birthseas, fleet, <age vec> where first value is for age 0!
// if yr=-yr, then fill remaining years for that seas, growpattern, gender, fleet
// fleet 0 contains begin season pop WT
// fleet -1 contains mid season pop WT
// fleet -2 contains maturity*fecundity

 LOCAL_CALCS
  if(k1>0)
  {
    echoinput<<"Wt_age input"<<endl<<WTage_in<<endl<<"end"<<endl;
    WTage_emp.initialize();
    if(N_WTage_maxage>nages) N_WTage_maxage=nages;  //  so extra ages being read will be ignored
    for(i=1;i<=N_WTage_rd;i++)
    {
      y=abs(WTage_in(i,1));
      if(y<styr) y=styr;
      if(WTage_in(i,1)<0) {y2=YrMax;} else {y2=y;}
      s=WTage_in(i,2);
      gg=WTage_in(i,3);
      gp=WTage_in(i,4);
      birthseas=WTage_in(i,5);
      g=(gg-1)*N_GP*nseas + (gp-1)*nseas + birthseas;
      f=WTage_in(i,6);
      if(s<=nseas && gg<=gender && gp<=N_GP && birthseas<=nseas && f<=Ntypes)
      {
        for(j=y;j<=y2;j++)  // loop years
        {
          t=styr+(j-styr)*nseas+s-1;
          for(a=0;a<=N_WTage_maxage;a++) WTage_emp(t,g,f,a)=WTage_in(i,7+a);
          for(a=N_WTage_maxage;a<=nages;a++) WTage_emp(t,g,f,a)=WTage_emp(t,g,f,N_WTage_maxage);  //  fills out remaining ages, if any
          if(j==y) echoinput<<y<<" s "<<s<<" sex "<<gg<<" gp "<<gp<<" bs "<<birthseas<<" morph "<<g<<" pop/fleet "<<f<<" "<<WTage_emp(t,g,f)(0,min(6,nages))<<endl;
        }
      }
      temp=float(Bmark_Yr(2)-Bmark_Yr(1)+1.);  //  get denominator
      for(f=-2;f<=Ntypes;f++)
      for(g=1;g<=gmorph;g++)
      if(use_morph(g)>0)
      {
        for(s=0;s<=nseas-1;s++)
        {
          junkvec2.initialize();
          for(t=Bmark_t(1);t<=Bmark_t(2);t+=nseas) {junkvec2+=WTage_emp(t,GP3(g),f);}
          WTage_emp(styr-3*nseas+s,GP3(g),f)=junkvec2/temp;
        }
      }
    }
  }
 END_CALCS

//****************************************


//  INITIALIZATION_SECTION (not used in SS)
INITIALIZATION_SECTION

//  ProgLabel_5.0   PARAMETER_SECTION
PARAMETER_SECTION

//  ProgLabel_5.0.1  setup convergence critera and max func evaluations
 LOCAL_CALCS
    if(readparfile>=1)
    {cout<<" read parm file"<<endl;
    ad_comm::change_pinfile_name("ss3.par");}
    maximum_function_evaluations.allocate(func_eval.indexmin(),func_eval.indexmax());
    maximum_function_evaluations=func_eval;
    convergence_criteria.allocate(func_conv.indexmin(),func_conv.indexmax());
    convergence_criteria=func_conv;
 END_CALCS

//  ProgLabel_5.0.2  create dummy_parm that will be estimated even if turn_off_phase is set to 0
  init_bounded_number dummy_parm(0,2,dummy_phase)  //  estimate in phase 0

//  ProgLabel_5.1.1  create MGparm vector and associated arrays
  // natural mortality and growth
  init_bounded_number_vector MGparm(1,N_MGparm2,MGparm_LO,MGparm_HI,MGparm_PH)
  matrix MGparm_trend(1,N_MGparm_trend,styr,YrMax+1);
  matrix MGparm_block_val(1,N_MGparm,styr,YrMax+1);
  init_bounded_matrix MGparm_dev(1,N_MGparm_dev,MGparm_dev_minyr,MGparm_dev_maxyr,-10,10,MGparm_dev_PH)
  matrix MGparm_dev_rwalk(1,N_MGparm_dev,MGparm_dev_minyr,MGparm_dev_maxyr);
  vector L_inf(1,N_GP*gender);
  vector Lmax_temp(1,N_GP*gender);
  matrix VBK(1,N_GP*gender,0,nages);
  vector Richards(1,N_GP*gender);

  vector Lmin(1,N_GP*gender);
  vector Lmin_last(1,N_GP*gender);
//  vector natM1(1,N_GP*gender)
//  vector natM2(1,N_GP*gender)
  matrix natMparms(1,N_natMparms,1,N_GP*gender)
  3darray natM(1,nseas,1,N_GP*gender*nseas,0,nages)   //  need nseas to capture differences due to birthseason
  3darray surv1(1,nseas,1,N_GP*gender*nseas,0,nages)
  3darray surv2(1,nseas,1,N_GP*gender*nseas,0,nages)
  vector CVLmin(1,N_GP*gender)
  vector CVLmax(1,N_GP*gender)
  vector CV_const(1,N_GP*gender)
  matrix mgp_save(styr,endyr+1,1,N_MGparm2);
  vector mgp_adj(1,N_MGparm2);
  matrix Cohort_Growth(styr,YrMax,0,nages)
  3darray Cohort_Lmin(1,N_GP*gender,styr,YrMax,0,nages)
  vector VBK_seas(0,nseas);
  matrix wtlen_seas(0,nseas,1,8);
  vector wtlen_p(1,8);

  3darray wt_len(1,nseas,1,gender,1,nlength)  //  stores wt at mid-bin
  matrix wt_len2(1,nseas,1,nlength2)    //  stores wt at midbin; stacked genders
  matrix wt_len2_sq(1,nseas,1,nlength2)    //  stores wt at midbin^2; stacked genders
  matrix wt_len_low(1,nseas,1,nlength2)  //  wt at lower edge of size bin
  matrix wt_len_fd(1,nseas,1,nlength2-1)  //  first diff of wt_len_low
  matrix mat_len(1,N_GP,1,nlength)
  matrix fec_len(1,N_GP,1,nlength)   // fecundity at length
  matrix mat_fec_len(1,N_GP,1,nlength)
  matrix mat_age(1,N_GP,0,nages)
  matrix Hermaphro_val(1,N_GP,0,nages)

 LOCAL_CALCS
   mat_len=1.0;
   mat_age=1.0;
   mat_fec_len=1.0;
   fec_len=1.0;
 END_CALCS
  3darray age_age(0,N_ageerr,1,n_abins2,0,gender*nages+gender-1)
  3darray age_err(1,N_ageerr,1,2,0,nages) // ageing imprecision as stddev for each age

// Age-length keys for each gmorph
  4darray ALK(1,2*nseas,1,gmorph,0,nages,1,nlength)
  matrix exp_AL(0,gender*nages+gender-1,1,nlength2);
  3darray Sd_Size_within(1,2*nseas,1,gmorph,0,nages)  //  2*nseas stacks begin of seas and end of seas
  3darray Sd_Size_between(1,2*nseas,1,gmorph,0,nages)  //  2*nseas stacks begin of seas and end of seas
  4darray Ave_Size(styr-3*nseas,TimeMax_Fcast_std+nseas,0,1,1,gmorph,0,nages)
  3darray CV_G(1,N_GP*gender,1,2*nseas,0,nages);   //  temporary storage of CV enroute to sd of len-at-age
  3darray Save_Wt_Age(styr-3*nseas,TimeMax_Fcast_std+1,1,gmorph,0,nages)
  3darray Wt_Age_beg(1,nseas,1,gmorph,0,nages)
  3darray Wt_Age_mid(1,nseas,1,gmorph,0,nages)

  3darray migrrate(styr-3,endyr+1,1,do_migr2,0,nages)
  3darray recr_dist_calc(1,N_GP*gender,1,nseas,1,pop);
  3darray recr_dist_benchmark(1,N_GP*gender,1,nseas,1,pop);
  3darray recr_dist_forecast(1,N_GP*gender,1,nseas,1,pop);

//  ProgLabel_5.1.2  create SP_parm vector, recruitment vectors
  init_bounded_number_vector SR_parm(1,N_SRparm2,SRvec_LO,SRvec_HI,SRvec_PH)
  number two_sigmaRsq;
  number half_sigmaRsq;
  number sigmaR
  number rho;
  !!echoinput<<"here2"<<endl;
 LOCAL_CALCS
  if(SR_parm(N_SRparm2)!=0.0 || SRvec_PH(N_SRparm2)>0) {SR_autocorr=1;} else {SR_autocorr=0;}  // flag for recruitment autocorrelation
  if(do_recdev==1)
  {k=recdev_start; j=recdev_end; s=1; p=-1;}
  else if(do_recdev==2)
  {s=recdev_start; p=recdev_end; k=1; j=-1;}
  else
  {s=1; p=-1; k=1; j=-1;}

 END_CALCS

  vector biasadj(styr-nages,YrMax)  // biasadj as used; depends on whether a recdev is estimated or not
  vector biasadj_full(styr-nages,YrMax)  //  full time series of biasadj values, only used in defined conditions
  number sd_offset_rec

  init_bounded_number_vector recdev_cycle_parm(1,recdev_cycle,recdev_cycle_LO,recdev_cycle_HI,recdev_cycle_PH)

//  init_bounded_dev_vector recdev_early(recdev_early_start,recdev_early_end,recdev_LO,recdev_HI,recdev_early_PH)
  init_bounded_vector recdev_early(recdev_early_start,recdev_early_end,recdev_LO,recdev_HI,recdev_early_PH)
  init_bounded_dev_vector recdev1(k,j,recdev_LO,recdev_HI,recdev_PH)
  init_bounded_vector recdev2(s,p,recdev_LO,recdev_HI,recdev_PH)
  init_bounded_vector Fcast_recruitments(recdev_end+1,YrMax,recdev_LO,recdev_HI,Fcast_recr_PH2)
  vector recdev(recdev_first,YrMax);

 LOCAL_CALCS
  if(Do_Impl_Error>0)
//  {k=max_phase+1;}
  {k=Fcast_recr_PH2;}
  else
  {k=-1;}
 END_CALCS
  init_bounded_vector Fcast_impl_error(endyr+1,YrMax,-1,1,k)

  number SPB_current;                            // Spawning biomass
  number SPB_vir_LH
  number Recr_virgin
  number SPB_virgin
  number SPR_unf
  number SPR_trial
//  vector S1(0,1);
  3darray SPB_pop_gp(styr-3,YrMax,1,pop,1,N_GP)         //Spawni Output
  vector SPB_yr(styr-3,YrMax)
  vector SPB_B_yr(styr-3,YrMax)  //  mature biomass (no fecundity)
  vector SPB_N_yr(styr-3,YrMax)   //  mature numbers
  number equ_mat_bio
  number equ_mat_num
  !!k=0;
  !!if(Hermaphro_Option>0) k=1;

  3darray MaleSPB(styr-3,YrMax*k,1,pop,1,N_GP)         //Male Spawning biomass

  matrix SPB_equil_pop_gp(1,pop,1,N_GP);
  matrix MaleSPB_equil_pop_gp(1,pop,1,N_GP);
  number SPB_equil;
  number SPR_temp;  //  used to pass quantity into Equil_SpawnRecr
  number Recruits;                            // Age0 Recruits
  matrix Recr(1,pop,styr-3,YrMax)         //Recruitment
  matrix exp_rec(styr-3,YrMax,1,4) //expected value for recruitment: 1=spawner-recr only; 2=with environ and cycle; 3=with bias_adj; 4=with dev
  matrix Nmid(1,gmorph,0,nages);
  matrix Nsurv(1,gmorph,0,nages);
  3darray natage_temp(1,pop,1,gmorph,0,nages)
  vector agetemp(0,gender*nages+gender-1)
  4darray Save_PopLen(styr-3*nseas,TimeMax_Fcast_std+1,1,2*pop,1,gmorph,1,nlength)
  4darray Save_PopWt(styr-3*nseas,TimeMax_Fcast_std+1,1,2*pop,1,gmorph,1,nlength)
  4darray Save_PopAge(styr-3*nseas,TimeMax_Fcast_std+1,1,2*pop,1,gmorph,0,nages)

  number ave_age    //  average age of fish in unfished population; used to weight R1

//  ProgLabel_5.1.3  create F parameters and associated arrays and constants
  init_bounded_number_vector init_F(1,Nfleet,init_F_LO,init_F_HI,init_F_PH)
  vector est_equ_catch(1,Nfleet)

  !!if(Do_Forecast>0) {k=TimeMax_Fcast_std+1;} else {k=TimeMax+nseas;}
  4darray natage(styr-3*nseas,k,1,pop,1,gmorph,0,nages)  //  add +1 year
  4darray catage(styr-nseas,TimeMax,1,Nfleet,1,gmorph,0,nages)
  4darray equ_catage(1,nseas,1,Nfleet,1,gmorph,0,nages)
  4darray equ_numbers(1,nseas,1,pop,1,gmorph,0,3*nages)
  4darray equ_Z(1,nseas,1,pop,1,gmorph,0,nages)
  matrix catage_tot(1,gmorph,0,nages)//sum the catches for all fleets, reuse matrix each year
  matrix Hrate(1,Nfleet,styr-3*nseas,k) //Harvest Rate for each fleet
  3darray catch_fleet(styr-3*nseas,TimeMax_Fcast_std,1,Nfleet,1,6)  //  1=sel_bio, 2=kill_bio; 3=ret_bio; 4=sel_num; 5=kill_num; 6=ret_num
  matrix annual_catch(styr,YrMax,1,6)  //  same six as above
  matrix annual_F(styr,YrMax,1,2)  //  1=sum of hrate (if Pope fmethod) or sum hrate*seasdur if F; 2=Z-M for selected ages
  matrix equ_catch_fleet(1,Nfleet,1,6)

  matrix fec(1,gmorph,0,nages)            //relative fecundity at age, is the maturity times the weight-at-age times eggs/kg for females
  matrix make_mature_bio(1,gmorph,0,nages)  //  mature female weight at age
  matrix make_mature_numbers(1,gmorph,0,nages)  //  mature females at age
  matrix virg_fec(1,gmorph,0,nages)
  number fish_bio;
  number fish_bio_r;
  number fish_bio_e;
  number fish_num_e;
  number fish_num;
  number fish_num_r;
  number vbio;
  number totbio;
  number smrybio;
  number smrynum;
  number smryage;  // mean age of the summary numbers (not accounting for birthseason)
  number catch_mnage;  //  mean age of the catch (not accounting for birthseason or season of the catch)
  number catch_mnage_d;  // total catch numbers for calc of mean age
  number harvest_rate;                        // Harvest rate
  number maxpossF;
  vector Get_EquilCalc(1,2);

  4darray Z_rate(styr-3*nseas,k,1,pop,1,gmorph,0,nages)
  3darray Zrate2(1,pop,1,gmorph,0,nages)

 LOCAL_CALCS
  if(F_Method==2)    // continuous F
    {k=Nfleet*(TimeMax-styr+1);}
  else
    {k=-1;}
 END_CALCS
 init_bounded_number_vector F_rate(1,k,0.,F_rate_max,have_catch)

  vector Nmigr(1,pop);
  number Nsurvive;
  number YPR_tgt_enc;
  number YPR_tgt_dead;
  number YPR_tgt_N_dead;
  number YPR_tgt_ret;
  number YPR_spr; number Vbio_spr; number Vbio1_spr; number SPR_actual;

  number YPR_Btgt_enc;
  number YPR_Btgt_dead;
  number YPR_Btgt_N_dead;
  number YPR_Btgt_ret;
  number YPR_Btgt; number Vbio_Btgt; number Vbio1_Btgt;
  number Btgt; number Btgttgt; number SPR_Btgt; number Btgt_Rec;
  number Bspr; number Bspr_rec;

  number YPR    // variable still used in SPR series
  number MSY
  number Bmsy
  number Recr_msy
  number YPR_msy_enc;
  number YPR_msy_dead;
  number YPR_msy_N_dead;
  number YPR_msy_ret;

  number YPR_enc;
  number YPR_dead;
  number YPR_N_dead;
  number YPR_ret;
  number MSY_Fmult;
  number SPR_Fmult;
  number Btgt_Fmult;

  number caa;
   number Fmult;
   number Fcast_Fmult;
   number Fcurr_Fmult;
   number Fchange;
   number last_calc;
   matrix Fcast_RelF_Use(1,nseas,1,Nfleet);
   matrix Bmark_RelF_Use(1,nseas,1,Nfleet);
   number alpha;
   number beta;
   number MSY_SPR;
   number GenTime;
   vector cumF(1,gmorph);
   vector maxF(1,gmorph);
   number Yield;
   number Adj4010;

//  !!k1 = styr+(endyr-styr)*nseas-1 + nseas + 1;
//  !!y=k1+N_Fcast_Yrs*nseas-1;

//  ProgLabel_5.1.4  create Q_parm and associated arrays
  init_bounded_number_vector Q_parm(1,Q_Npar,Q_parm_LO,Q_parm_HI,Q_parm_PH)

  matrix log_q_cr(1,Ntypes,1,nyr_cr);
  matrix q_cr(1,Ntypes,1,nyr_cr);
  matrix se_cr_use(1,Ntypes,1,nyr_cr)
  matrix exp_cr(1,Ntypes,1,nyr_cr)    //  will store expected survey
  vector surv_like(1,Ntypes) // likelihood of the indices
  matrix Q_dev_like(1,Ntypes,1,2) // likelihood of the Q deviations

  vector disc_like(1,Ntypes) // likelihood of the discard biomass
  vector mnwt_like(1,Ntypes) // likelihood of the mean body wt

  matrix exp_disc(1,Ntypes,1,nyr_disc)
  3darray retain(styr-3,endyr+1,1,Ntypes,1,nlength2)
  vector retain_M(1,nlength)
  3darray discmort(styr-3,endyr+1,1,Ntypes,1,nlength2)
  vector discmort_M(1,nlength)
  vector exp_mnwt(1,nobs_mnwt)

  matrix Morphcomp_exp(1,Morphcomp_nobs,6,5+Morphcomp_nmorph)   // expected value for catch by growthpattern

  3darray SzFreqTrans(1,SzFreq_Nmeth*nseas,1,nlength2,1,SzFreq_Nbins_seas_g);

//  ProgLabel_5.1.5  Selectivity-related parameters
  init_bounded_number_vector selparm(1,N_selparm2,selparm_LO,selparm_HI,selparm_PH)
  matrix selparm_trend(1,N_selparm_trend,styr,endyr);
  matrix selparm_block_val(1,N_selparm,styr,endyr);

  init_bounded_matrix selparm_dev(1,N_selparm_dev,selparm_dev_minyr,selparm_dev_maxyr,-10,10,selparm_dev_PH)
  matrix selparm_dev_rwalk(1,N_selparm_dev,selparm_dev_minyr,selparm_dev_maxyr)
  4darray sel_l(styr-3,endyr+1,1,Ntypes,1,gender,1,nlength)
  4darray sel_l_r(styr-3,endyr+1,1,Ntypes,1,gender,1,nlength)   //  selex x retained
  4darray discmort2(styr-3,endyr+1,1,Ntypes,1,gender,1,nlength)
  4darray sel_a(styr-3,endyr+1,1,Ntypes,1,gender,0,nages)
  vector sel(1,nlength)  //  used to multiply by ALK

//  ProgLabel_5.1.6  create TG_parm and associated arrays
  matrix TG_alive(1,pop,1,gmorph)
  matrix TG_alive_temp(1,pop,1,gmorph)
  3darray TG_recap_exp(1,N_TG2,0,TG_endtime,0,Nfleet)   //  do not need to store POP index because each fleet is in just one area
  vector TG_like1(1,N_TG2)
  vector TG_like2(1,N_TG2)
  number overdisp     // overdispersion

 LOCAL_CALCS
  k=Do_TG*(3*N_TG+2*Nfleet);
 END_CALCS

  init_bounded_number_vector TG_parm(1,k,TG_parm_LO,TG_parm_HI,TG_parm_PH);

 LOCAL_CALCS
  if(Do_Forecast>0)
  {k=TimeMax_Fcast_std+1;}
  else
  {k=TimeMax+nseas;}
 END_CALCS

//  ProgLabel_5.2  Create arrays for storing derived selectivity quantities for use in mortality calculations
  4darray fish_body_wt(styr-3*nseas,k,1,gmorph,1,Nfleet,0,nages);  // wt (adjusted for size selex)
  4darray sel_al_1(1,nseas,1,gmorph,1,Nfleet,0,nages);  // selected * wt
  4darray sel_al_2(1,nseas,1,gmorph,1,Nfleet,0,nages);  // selected * retained * wt
  4darray sel_al_3(1,nseas,1,gmorph,1,Nfleet,0,nages);  // selected numbers
  4darray sel_al_4(1,nseas,1,gmorph,1,Nfleet,0,nages);  // selected * retained numbers
  4darray deadfish(1,nseas,1,gmorph,1,Nfleet,0,nages);  // sel * (retain + (1-retain)*discmort)
  4darray deadfish_B(1,nseas,1,gmorph,1,Nfleet,0,nages);  // sel * (retain + (1-retain)*discmort) * wt

  4darray save_sel_fec(styr-3*nseas,TimeMax_Fcast_std+nseas,1,gmorph,0,Nfleet,0,nages)  //  save sel_al_3 (Asel_2) and save fecundity for output;  +nseas covers no forecast setups

  4darray Sel_for_tag(TG_timestart*Do_TG,TimeMax*Do_TG,1,gmorph*Do_TG,1,Nfleet,0,nages)
  vector TG_report(1,Nfleet*Do_TG);
  vector TG_rep_decay(1,Nfleet*Do_TG);

  3darray save_sp_len(styr,endyr+1,1,2*Ntypes,1,90);     // use to output selex parm values after adjustment

  3darray exp_l(1,Ntypes,1,Nobs_l,1,nlen_bin2)
  matrix neff_l(1,Ntypes,1,Nobs_l)
  vector tempvec_l(1,nlength);
  vector exp_l_temp(1,nlength2);
  vector exp_l_temp_ret(1,nlength2);     // retained lengthcomp
  vector exp_l_temp_dat(1,nlen_bin2);
  vector offset_l(1,Ntypes) // Compute OFFSET for multinomial (i.e, value for the multinonial function
  vector length_like(1,Ntypes)  // likelihood of the length-frequency data

  matrix SzFreq_exp(1,SzFreq_totobs,1,SzFreq_Setup2);
  vector SzFreq_like(1,SzFreq_N_Like)
  3darray exp_a(1,Ntypes,1,Nobs_a,1,n_abins2)
  3darray  exp_meanage(1,Ntypes,1,Nobs_a,1,3)  //  will hold mean age' and 95% range for the range of sizes identified for this age comp observation
  vector exp_a_temp(1,n_abins2)
  vector tempvec_a(0,nages)
  matrix neff_a(1,Ntypes,1,Nobs_a)
  vector offset_a(1,Ntypes) // Compute OFFSET for multinomial (i.e, value for the multinonial function
  vector age_like(1,Ntypes)  // likelihood of the age-frequency data
  vector sizeage_like(1,Ntypes)  // likelihood of the age-frequency data
  3darray exp_ms(1,Ntypes,1,Nobs_ms,1,n_abins2)
  3darray exp_ms_sq(1,Ntypes,1,Nobs_ms,1,n_abins2)

  number Morphcomp_like
  number equ_catch_like
  vector catch_like(1,Ntypes)
  number recr_like
  number Fcast_recr_like
  number parm_like
  number parm_dev_like
  number CrashPen
  number SoftBoundPen
  number Equ_penalty
  number F_ballpark_like

  number R1
  number R1_exp
  number t1
  number t2
  number temp
  number temp1
  number temp2
  number temp3
  number temp4
  number join1
  number join2
  number join3
  number upselex
  number downselex
  number peak
  number peak2
  number point1
  number point2
  number point3
  number point4
  number timing
  number equ_Recr
  number equ_F_std

//  ProgLabel_5.1.7  create matrix called smry to store derived quantities of interest
  matrix Smry_Table(styr-3,YrMax,1,20+2*gmorph);
  // 1=totbio, 2=smrybio, 3=smrynum, 4=enc_catch, 5=dead_catch, 6=ret_catch, 7=spbio, 8=recruit,
  // 9=equ_totbio, 10=equ_smrybio, 11=equ_SPB_virgin, 12=equ_S1, 13=Gentime, 14=YPR, 15=meanage_spawners, 16=meanage_smrynums, 17=meanage_catch
  // 18, 19, 20  not used
  // 21+cumF-bymorph, maxF-by morph

  matrix env_data(styr-1,YrMax,-2,N_envvar)
  matrix TG_save(1,N_TG,1,3+TG_endtime)

//  ProgLabel_5.1.8  create sdreport vectors
  sdreport_vector SPB_std(1,N_STD_Yr);
  sdreport_vector recr_std(1,N_STD_Yr);
  sdreport_vector SPR_std(1,N_STD_Yr_Ofish);
  sdreport_vector F_std(1,N_STD_Yr_F);
  sdreport_vector depletion(1,N_STD_Yr_Dep);
  sdreport_vector Mgmt_quant(1,N_STD_Mgmt_Quant)
  sdreport_vector Extra_Std(1,Extra_Std_N)

//  ProgLabel_5.1.9  create logLikelihood vectors
  vector MGparm_Like(1,N_MGparm2)
  vector init_F_Like(1,Nfleet)
  vector Q_parm_Like(1,Q_Npar)
  vector selparm_Like(1,N_selparm2)
  vector SR_parm_Like(1,N_SRparm2)
  vector recdev_cycle_Like(1,recdev_cycle)
  !! k=Do_TG*(3*N_TG+2*Nfleet);
  vector TG_parm_Like(1,k);
  objective_function_value obj_fun
  number last_objfun
  vector phase_output(1,max_phase+1)
  !!echoinput<<"end of parameter section"<<endl;

//******************************************************************************************
//  ProgLabel_6.0  PRELIMINARY_CALCS_SECTION
PRELIMINARY_CALCS_SECTION
//  ProgLabel_6.1  Some initial housekeeping
//  ProgLabel_6.1.1  create and initialize random number generator
  random_number_generator radm(long(time(&start)));
  if(F_ballpark_yr>retro_yr) F_ballpark_yr=endyr;
  if(F_ballpark_yr<styr) {F_ballpark_lambda=0.;}
  sel_l.initialize();
  sel_a.initialize();
  offset_l.initialize();
  offset_a.initialize();
  save_sp_len.initialize();
  save_sel_fec.initialize();

//  ProgLabel_6.1.2  initialize the dummy parameter as needed
  if(Turn_off_phase<=0) {dummy_parm=0.5;} else {dummy_parm=1.0;}

  Cohort_Growth=1.0;    // adjustment for cohort growth deviations

//  ProgLabel_6.2  Apply input variance adjustments to each data type
//  ProgLabel_6.2.1  do variance adjustment for surveys

  echoinput<<" do variance adjustment for surveys "<<endl;
  for (f=1; f<=Ntypes; f++)
  if(nyr_cr(f)>0)
  {
    for (i=1; i<=nyr_cr(f); i++)
    {
      if(yr_cr_use(f,i)>0)
      {
        Ln_obs_cr(f,i)=obs_cr(f,i);
        if(cr_errtype(f)>=0)  // lognormal or lognormal T_dist
        {
          if(obs_cr(f,i)<=0.0)
          {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Survey obs must be positive for lognormal error"<<endl; exit(1);}
          Ln_obs_cr(f,i)=log(obs_cr(f,i));
          se_cr_obs(f,i)+=var_adjust(1,f);
          if(se_cr_obs(f,i)<=0.0) se_cr_obs(f,i)=0.001;
        }
        else  // normal distribution
        {
          se_cr_obs(f,i)+=var_adjust(1,f);
          if(se_cr_obs(f,i)<=0.0) se_cr_obs(f,i)=0.001;
        }
      }
    }
  }
  echoinput<<" survey stderr has been set-up "<<endl;

//  ProgLabel_6.2.2  set up variance for discard observations
  for (f=1; f<=Ntypes; f++)
  if(nyr_disc(f)>0)
  {
    for (i=1; i<=nyr_disc(f); i++)
    {
      if(yr_disc_use(f,i)>=0.)
      {
        if(cv_disc(f,i)<=0.0) cv_disc(f,i)=0.001;
        if(disc_errtype(f)>=0)
        {sd_disc(f,i)=cv_disc(f,i)*obs_disc(f,i);}
        else
        {sd_disc(f,i)=cv_disc(f,i);}

        sd_disc(f,i)+=var_adjust(2,f);   // note that adjustment is to the sd, not the CV
        if(sd_disc(f,i)<0.001) sd_disc(f,i)=0.001;
      }
    }
  }
  echoinput<<" discard stderr has been set-up "<<endl;

//  ProgLabel_6.2.3  set up variance for mean body wt data, note different reference to array that was read
  for(i=1;i<=nobs_mnwt;i++)
  {
    if(mnwtdata(5,i)>0. && mnwtdata(3,i)>0.)  //  used observation
    {
      mnwtdata(6,i)+=var_adjust(3,mnwtdata(3,i));
      if(mnwtdata(6,i)<=0.0) mnwtdata(6,i)=0.001;
      mnwtdata(7,i)=mnwtdata(5,i)*mnwtdata(6,i);
      mnwtdata(8,i)=DF_bodywt*square(mnwtdata(7,i));
      mnwtdata(9,i)=sd_offset*log(mnwtdata(7,i));
    }
  }
  echoinput<<" mean bodywt stderr has been set-up "<<endl;

//  ProgLabel_6.2.4  do variance adjustment and compute OFFSET for length comp
  if(Nobs_l_tot>0)
  for (f=1; f <= Ntypes; f++)
  for (i=1; i <= Nobs_l(f); i++)
  if(header_l(f,i,3)>0)
  {
    nsamp_l(f,i)*=var_adjust(4,f);  {if(nsamp_l(f,i)<=1.0) nsamp_l(f,i)=1.;}                              //  adjust sample size
    if(gen_l(f,i) !=2) offset_l(f) -= nsamp_l(f,i) *
    obs_l(f,i)(tails_l(f,i,1),tails_l(f,i,2))*log(obs_l(f,i)(tails_l(f,i,1),tails_l(f,i,2)));
    if(gen_l(f,i) >=2 && gender==2) offset_l(f) -= nsamp_l(f,i) *
    obs_l(f,i)(tails_l(f,i,3),tails_l(f,i,4))*log(obs_l(f,i)(tails_l(f,i,3),tails_l(f,i,4)));
  }
  echoinput<<" length_comp offset: "<<offset_l<<endl;

//  ProgLabel_6.2.4.1  //  get sample weights for the super-period components in length comp
//  the combined obs will have a logL sample size equal to the sample size input for the accumulator observation
//  the accumulator observation is assigned a weight of 1.0 (because there is no place to read this from)
//  the obs to be combined with the accumulator get a weight equal to value input in the nsamp_l element
//  so, nsamp_l can no longer have negative observations
  for (f=1; f<=Ntypes; f++)
  {
    if(N_suprper_cr(f)>0)
    {
      echoinput<<"Create superperiod sample weights for survey obs"<<endl<<"Fleet Super Obs SE_input samp_wt"<<endl;
      for(j=1;j<=N_suprper_cr(f);j++)                  // do each super period
      {
        temp=1.0;  //  relative sample weight for time period the accumulator observation
        k=0;  // count of samples with real information
        for(i=suprper_cr1(f,j);i<=suprper_cr2(f,j);i++)  //  loop obs of this super period
        {
          if(yr_cr_use(f,i)<0)  //  so one of the obs to be combined
          {temp+=se_cr_obs(f,i);}  //  add in its weight relative to 1.0 for the observation with real info
          else
          {k++;}
        }
        if(k!=1) {N_warn++; cout<<"EXIT - see warning.sso"; warning<<" must have only 1 sample with real info in survey superperiod "<<j<<endl; exit(1);}
        for(i=suprper_cr1(f,j);i<=suprper_cr2(f,j);i++)
        {
          if(yr_cr_use(f,i)<0)  //  so one of the obs to be combined
          {suprper_cr_sampwt(f,i)=se_cr_obs(f,i)/value(temp);}
          else
          {suprper_cr_sampwt(f,i)=1.0/value(temp);}
          echoinput<<f<<" "<<j<<" "<<i<<" "<<yr_cr_use(f,i)<<" "<<se_cr_obs(f,i)<<" "<<suprper_cr_sampwt(f,i)<<endl;
        }
      }
    }

    if(N_suprper_disc(f)>0)
    {
      echoinput<<"Create superperiod sample weights for discard obs"<<endl<<"fleet Super Obs SE_input samp_wt"<<endl;
      for(j=1;j<=N_suprper_disc(f);j++)                  // do each super period
      {
        temp=1.0;  //  relative sample weight for time period the accumulator observation
        k=0;  // count of samples with real information
        for(i=suprper_disc1(f,j);i<=suprper_disc2(f,j);i++)  //  loop obs of this super period
        {
          if(yr_disc_use(f,i)<0)  //  so one of the obs to be combined
          {temp+=cv_disc(f,i);}  //  add in its weight relative to 1.0 for the observation with real info
          else
          {k++;}
        }
        if(k!=1) {N_warn++; cout<<"EXIT - see warning.sso"; warning<<" must have only 1 sample with real info in survey superperiod "<<j<<endl; exit(1);}
        for(i=suprper_disc1(f,j);i<=suprper_disc2(f,j);i++)
        {
          if(yr_disc_use(f,i)<0)  //  so one of the obs to be combined
          {suprper_disc_sampwt(f,i)=cv_disc(f,i)/value(temp);}
          else
          {suprper_disc_sampwt(f,i)=1.0/value(temp);}
          echoinput<<f<<" "<<j<<" "<<i<<" "<<yr_disc_use(f,i)<<" "<<cv_disc(f,i)<<" "<<suprper_disc_sampwt(f,i)<<endl;
        }
      }
    }

    if(N_suprper_l(f)>0)
    {
      echoinput<<"Create superperiod sample weights for length obs"<<endl<<"Super Obs effN_input samp_wt"<<endl;
      for(j=1;j<=N_suprper_l(f);j++)                  // do each super period
      {
        temp=1.0;  //  relative sample weight for time period the accumulator observation
        k=0;  // count of samples with real information
        for(i=suprper_l1(f,j);i<=suprper_l2(f,j);i++)  //  loop obs of this super period
        {
          if(header_l(f,i,3)<0)  //  so one of the obs to be combined
          {temp+=nsamp_l(f,i);}
          else
          {k++;}
        }
        if(k!=1) {N_warn++; cout<<"error in length data"; warning<<" must have only 1 sample with real info in length superperiod "<<j<<endl; exit(1);}
        for(i=suprper_l1(f,j);i<=suprper_l2(f,j);i++)
        {
          if(header_l(f,i,3)<0)  //  so one of the obs to be combined
          {suprper_l_sampwt(f,i)=nsamp_l(f,i)/value(temp);}
          else
          {suprper_l_sampwt(f,i)=1.0/value(temp);}
          echoinput<<f<<" "<<j<<" "<<i<<" "<<header_l(f,i,3)<<" "<<nsamp_l(f,i)<<" "<<suprper_l_sampwt(f,i)<<endl;
        }
      }
    }

    if(N_suprper_a(f)>0)
    {
      echoinput<<"Create superperiod sample weights for age obs"<<endl<<"Super Obs effN_input samp_wt"<<endl;
      for(j=1;j<=N_suprper_a(f);j++)                  // do each super period
      {
        temp=1.0;  //  relative sample weight for time period the accumulator observation
        k=0;  // count of samples with real information
        for(i=suprper_a1(f,j);i<=suprper_a2(f,j);i++)  //  loop obs of this super period
        {
          if(header_a(f,i,3)<0)  //  so one of the obs to be combined
          {temp+=nsamp_a(f,i);}
          else
          {k++;}
        }
        if(k!=1) {N_warn++; cout<<"error in age data"; warning<<" must have only 1 sample with real info in age superperiod "<<j<<endl; exit(1);}
        for(i=suprper_a1(f,j);i<=suprper_a2(f,j);i++)
        {
          if(header_a(f,i,3)<0)  //  so one of the obs to be combined
          {suprper_a_sampwt(f,i)=nsamp_a(f,i)/value(temp);}
          else
          {suprper_a_sampwt(f,i)=1.0/value(temp);}  //  for the element holding the combined observation
          echoinput<<f<<" "<<j<<" "<<i<<" "<<header_a(f,i,3)<<" "<<nsamp_a(f,i)<<" "<<suprper_a_sampwt(f,i)<<endl;
        }
      }
    }
    if(N_suprper_ms(f)>0)
    {
      echoinput<<"Create superperiod sample weights for meansize obs"<<endl<<"Super Obs effN_input samp_wt"<<endl;
      for(j=1;j<=N_suprper_ms(f);j++)                  // do each super period
      {
        temp=1.0;  //  relative sample weight for time period the accumulator observation
        k=0;  // count of samples with real information
        for(i=suprper_ms1(f,j);i<=suprper_ms2(f,j);i++)  //  loop obs of this super period
        {
          if(header_ms(f,i,3)<0)  //  so one of the obs to be combined
          {temp+=header_ms(f,i,7);}
          else
          {k++;}
        }
        if(k!=1) {N_warn++; cout<<"error in meansize data"; warning<<" must have only 1 sample with real info in meansize superperiod "<<j<<endl; exit(1);}
        for(i=suprper_ms1(f,j);i<=suprper_ms2(f,j);i++)
        {
          if(header_ms(f,i,3)<0)  //  so one of the obs to be combined
          {suprper_ms_sampwt(f,i)=header_ms(f,i,7)/value(temp);}
          else
          {suprper_ms_sampwt(f,i)=1.0/value(temp);}  //  for the element holding the combined observation
          echoinput<<f<<" "<<j<<" "<<i<<" "<<header_ms(f,i,3)<<" "<<header_ms(f,i,7)<<" "<<suprper_ms_sampwt(f,i)<<endl;
        }
      }
    }
  }
  
  if(N_suprper_SzFreq>0)
  {
    echoinput<<"Create superperiod sample weights for sizecomp obs "<<endl<<"Fleet Super OBS Super fleet Sample_N_read samp_wt"<<endl;
    for(j=1;j<=N_suprper_SzFreq;j++)                  // do each super period
    {
      temp=1.0;  //  relative sample weight for time period the accumulator observation
      k=0;  // count of samples with real information
      for(iobs=suprper_SzFreq_start(j);iobs<=suprper_SzFreq_end(j);iobs++)  //  loop obs of this super period
      {
        if(SzFreq_obs_hdr(iobs,3)<0)  //  so one of the obs to be combined
        {temp+=SzFreq_sampleN(iobs);}
        else
        {k++;}  //  so counts the obs that are not just placeholders
      }
      if(k!=1) {N_warn++; cout<<"error in sizecomp data"; warning<<" must have only 1 sample with real info in sizecomp superperiod "<<j<<endl; exit(1);}
      for(iobs=suprper_SzFreq_start(j);iobs<=suprper_SzFreq_end(j);iobs++)
      {
        if(SzFreq_obs_hdr(iobs,3)<0)  //  so one of the obs to be combined
        {suprper_SzFreq_sampwt(iobs)=SzFreq_sampleN(iobs)/value(temp);}
        else
        {suprper_SzFreq_sampwt(iobs)=1.0/value(temp);}  //  for the element holding the combined observation
        echoinput<<SzFreq_obs_hdr(iobs,3)<<" "<<j<<" "<<iobs<<" "<<SzFreq_sampleN(iobs)<<" "<<suprper_SzFreq_sampwt(iobs)<<endl;
      }
    }
  }

//  ProgLabel_6.2.5  do variance adjustment and compute OFFSET for age comp
  if(Nobs_a_tot>0)
  for (f=1; f <= Ntypes; f++)
  for (i=1; i <= Nobs_a(f); i++)
  if(nsamp_a(f,i)>0 && header_a(f,i,3)>0)
  {
    nsamp_a(f,i)*=var_adjust(5,f);
    {if(nsamp_a(f,i)<=1.0) nsamp_a(f,i)=1.;}                                //  adjust sample size
  if(gen_a(f,i) !=2) offset_a(f) -= nsamp_a(f,i) *
    obs_a(f,i)(tails_a(f,i,1),tails_a(f,i,2))*log(obs_a(f,i)(tails_a(f,i,1),tails_a(f,i,2)));
  if(gen_a(f,i) >=2 && gender==2) offset_a(f) -= nsamp_a(f,i) *
    obs_a(f,i)(tails_a(f,i,3),tails_a(f,i,4))*log(obs_a(f,i)(tails_a(f,i,3),tails_a(f,i,4)));
  }
   echoinput<<" agecomp offset "<<offset_a<<endl;

//  ProgLabel_6.2.6  do variance adjustment for mean size-at-age data
  if(nobs_ms_tot>0)
  {for (f=1; f <= Ntypes; f++)
  for (i=1; i <= Nobs_ms(f); i++)
  for(b=1;b<=n_abins2;b++)
  {
    if(obs_ms_n(f,i,b)>0)
    {
      obs_ms_n(f,i,b)=sqrt(var_adjust(6,f)*obs_ms_n(f,i,b));
      if(obs_ms_n(f,i,b)<=1.0) obs_ms_n(f,i,b)=1.;
    }
  }
  }
  echoinput<<" setup stderr for mean size-at-age: "<<endl;

//  ProgLabel_6.2.7  Input variance adjustment not implemented for generalized size comp

//  ProgLabel_6.3  compute the reverse pointers:  array(fleet,time) contains reference to i'th element in data list
  yr_cr2.initialize();
  yr_disc2.initialize();

  for (f=1; f<=Ntypes; f++)
  {
    if(nyr_cr(f)>0)
    for (i=1; i<=nyr_cr(f); i++)
     {
      if(yr_cr_y(f,i)>=styr && yr_cr_y(f,i)<=retro_yr)
       {
         if(yr_cr2(f,yr_cr(f,i))>0)
         {
          cout<<" EXIT - see warning "<<endl;
          warning<<" replicate survey observations are not allowed (fleet, ith obs, year) "<<f<<" "<<i<<" "<<yr_cr(f,i)<<endl; exit(1);
        }
        else
        {
         yr_cr2(f,yr_cr(f,i))=i;have_data(f,yr_cr(f,i))=1;
       }     //  yr_cr(f,i) stores value of t so has period within year info
      }
     }
  echoinput<<"Create reverse pointers to survey observations"<<endl;

   if(nyr_disc(f)>0)
   for (i=1; i<=nyr_disc(f); i++)
       {yr_disc2(f,yr_disc(f,i))=i;
        have_data(f,yr_disc(f,i))=1;
     }
  echoinput<<"Create reverse pointers to discard observations"<<endl;

   if(Nobs_l(f)>0)  //  no reverse pointer here because cannot use it for both superperiods and for replicate obs
   for (i=1; i<=Nobs_l(f); i++)
   {
     have_data(f,yr_l_t(f,i))=1;
   }

   if(Nobs_a(f)>0)
   for (i=1; i <= Nobs_a(f); i++)
   {
     have_data(f,yr_a_t(f,i))=1;
   }

   if(Nobs_ms(f)>0)
   for (i=1; i <= Nobs_ms(f); i++)
   {
    have_data(f,yr_ms_t(f,i))=1;
   }
  }

       for (f=1;f<=Nfleet;f++)
       if(obs_equ_catch(f)<=0.001)
       {init_F(f)=0.0;}

//  ProgLabel_6.4  conditionally copy the initial parameter values read from the "CTL" file into the parameter arrays
//   skip this assignment if the parameters are being read from a "SS2.PAR" file

  if(readparfile==0)
  {
    echoinput<< " set parms to init values in CTL file "<<endl;
    for(i=1;i<=N_MGparm2;i++)
    {MGparm(i) = MGparm_RD(i);}  //  set vector of initial natmort and growth parms
    echoinput<< " MGparms OK "<<endl;

    for(i=1;i<=N_SRparm2;i++)
    {SR_parm(i)=SR_parm_1(i,3);}
    echoinput<< " SRR_parms OK "<<endl;

    if(recdev_cycle>0)
    {
      for(y=1;y<=recdev_cycle;y++)
      {
        recdev_cycle_parm(y)=recdev_cycle_parm_RD(y,3);
      }
    }

    if(recdev_do_early>0) recdev_early.initialize();
    if(Do_Forecast>0) Fcast_recruitments.initialize();
    if(Do_Forecast>0) Fcast_impl_error.initialize();

    if(do_recdev==1)
    {recdev1.initialize();}                // set devs to zero
    else if(do_recdev==2)
    {recdev2.initialize();}                // set devs to zero

    if(recdev_read>0)
    {
      for(j=1;j<=recdev_read;j++)
      {
        y=recdev_input(j,1);
        if(y>=recdev_first && y<=YrMax)
        {
          if(y<recdev_start)
          {
            recdev_early(y)=recdev_input(j,2);
          }
          else if (y<=recdev_end)
          {
            if(do_recdev==1)
            {recdev1(y)=recdev_input(j,2);}
            else if(do_recdev==2)
            {recdev2(y)=recdev_input(j,2);}
          }
          else
          {
            Fcast_recruitments(y)=recdev_input(j,2);
          }
        }
        else
        {
          N_warn++; warning<<" Trying to specify a recdev out of allowable range of years "<<y<<endl;
        }
      }
    }
    echoinput<< " rec_devs OK "<<endl;

// **************************************************
    for(i=1;i<=Q_Npar;i++)
    {Q_parm(i) = Q_parm_1(i,3);}    //  set vector of initial index Q parms
    if(Q_Npar>0) echoinput<< " Q_parms OK "<<endl;

    for(i=1;i<=Nfleet;i++)
    init_F(i) = init_F_RD(i);    //  set vector of initial parms
    echoinput<< " initF_parms OK "<<endl;

    if(F_Method==2)
    {
      if(readparfile==0)
      {
      g=0;
      for(f=1;f<=Nfleet;f++)
      for(y=styr;y<=endyr;y++)
      for(s=1;s<=nseas;s++)
      {
        t=styr+(y-styr)*nseas+s-1;
        g++;
        if(catch_ret_obs(f,t)>0.)  //have_catch
        {
          F_rate(g)=F_setup(1); Hrate(f,t)=F_setup(1);
        }
        else
        {F_rate(g)=0.; Hrate(f,t)=0.0;}
      }
      if(F_detail>0)
      {
        for(k=1;k<=F_detail;k++)
        {
          f=F_setup2(k,1); y=F_setup2(k,2); s=F_setup2(k,3);
          t=styr+(y-styr)*nseas+s-1;
          g=(f-1)*(TimeMax-styr+1)+(y-styr)*nseas+s;
          if(F_setup2(k,4)!=-999) 
            {F_rate(g)=F_setup2(k,4); Hrate(f,t)=F_setup2(k,4);}
        }
      }
       echoinput<< " Fmort_parms have been reset "<<endl;
      }
      else
      {
        echoinput<< " Fmort_parms obtained from ss3.par "<<endl;
//        g=0;
//        for(f=1;f<=Nfleet;f++)
//        for(y=styr;y<=endyr;y++)
//        for(s=1;s<=nseas;s++)
//        {
//          t=styr+(y-styr)*nseas+s-1;
//          g++;
//            Hrate(f,t)=F_rate(g);
//        }
      }
    }

    for(i=1;i<=N_selparm2;i++)
    selparm(i) = selparm_RD(i);    //  set vector of initial selex parms
    echoinput<< " selex_parms OK "<<endl;

    if(Do_TG>0)
    {
      k=Do_TG*(3*N_TG+2*Nfleet);
      for(i=1;i<=k;i++)
      {
        TG_parm(i)=TG_parm2(i,3);
      }
          echoinput<< " Tag_parms OK "<<endl;
    }
  }


//  ProgLabel_6.5  check parameter bounds and do jitter
    echoinput<< " now check bounds and do jitter if requested "<<endl;
    for(i=1;i<=N_MGparm2;i++)
    if(MGparm_PH(i)>0)
    {MGparm(i)=Check_Parm(MGparm_LO(i),MGparm_HI(i), jitter, MGparm(i));}
    echoinput<< " MG_parms OK "<<endl;

    for(i=1;i<=N_SRparm2;i++)
    if(SR_parm_1(i,7)>0)
    {SR_parm(i) = Check_Parm(SR_parm_1(i,1),SR_parm_1(i,2), jitter, SR_parm(i));}
    echoinput<< " SRR_parms OK "<<endl;

    if(recdev_do_early>0 && recdev_early_PH>0)
    {
    for(y=recdev_early_start;y<=recdev_early_end;y++)
      {recdev_early(y) = Check_Parm(recdev_LO, recdev_HI, jitter, recdev_early(y));}
//      recdev_early -=sum(recdev_early)/(recdev_early_end-recdev_early_start+1);
    }

    if(recdev_PH>0 && do_recdev>0)
    {
      if(do_recdev==1)
      {
        for(i=recdev_start;i<=recdev_end;i++)
        {recdev1(i) = Check_Parm(recdev_LO, recdev_HI, jitter, recdev1(i));}
        recdev1 -=sum(recdev1)/(recdev_end-recdev_start+1);
      }
      else
      {
        for(i=recdev_start;i<=recdev_end;i++)
        {recdev2(i) = Check_Parm(recdev_LO, recdev_HI, jitter, recdev2(i));}
//        recdev2 -=sum(recdev2)/(recdev_end-recdev_start+1);
      }
    }
    echoinput<< " rec_devs OK "<<endl;

    if(Q_Npar>0)
    {
      for(i=1;i<=Q_Npar;i++)
      if(Q_parm_1(i,7)>0)
      {Q_parm(i) = Check_Parm(Q_parm_1(i,1),Q_parm_1(i,2), jitter, Q_parm(i));}
      echoinput<< " Q_parms OK "<<endl;
    }

    for(i=1;i<=Nfleet;i++)
      {
      if(init_F_PH(i)>0)
        {init_F(i) = Check_Parm(init_F_LO(i),init_F_HI(i), jitter, init_F(i));}
        echoinput<< " initF_parms OK "<<endl;
      }

    for(i=1;i<=N_selparm2;i++)
    if(selparm_PH(i)>0)
    {selparm(i)=Check_Parm(selparm_LO(i),selparm_HI(i), jitter, selparm(i));}
    echoinput<< " selex_parms OK "<<endl;
    if(Do_TG>0)
    {
      k=Do_TG*(3*N_TG+2*Nfleet);
      for(i=1;i<=k;i++)
      {
      if(TG_parm_PH(i)>0)
        {TG_parm(i)=Check_Parm(TG_parm_LO(i),TG_parm_HI(i), jitter, TG_parm(i));}
      }
      echoinput<< " Tag_parms OK "<<endl;
    }
//  end bound check and jitter

//  ProgLabel_6.6  copy the environmental data as read into the dmatrix environmental data arrays
//  this will allow dynamic derived quantities like biomass and recruitment to be mapped into this same dmatrix
    env_data.initialize();
    if(N_envvar>=1)
      {
      for (y=styr-1;y<=(YrMax);y++)
        for(j=1;j<=N_envvar;j++)
          {env_data(y,j)=env_data_RD(y,j);}
      }

//  ProgLabel_6.7  initialize several items
    if(Rebuild_Ydecl==-1) Rebuild_Ydecl=1999;
    if(Rebuild_Yinit==-1) Rebuild_Yinit=endyr+1;

    if(Rebuild_Ydecl>YrMax) Rebuild_Ydecl=YrMax;
    if(Rebuild_Yinit>YrMax) Rebuild_Yinit=YrMax;

    migrrate.initialize();
    depletion.initialize();
    sel_l.initialize(); sel_a.initialize(); retain.initialize();  discmort.initialize(); discmort2.initialize();

    for(f=1;f<=Ntypes;f++)
    for(y=styr;y<=endyr+1;y++)
    for(gg=1;gg<=gender;gg++)
    {
      discmort2(y,f,gg)=1.0;
      if(y<=endyr+1)
      {
        discmort(y,f)=1.0;
        retain(y,f)=1.0;
      }
    }
    Richards=1.0;

//  ProgLabel_6.8  go thru biological calculations once, with do_once flag=1
    echoinput<<endl<<"**************************************"<<endl<<" ready to do biological and selectivity calculations once "<<endl;
    do_once=1;
    niter=0;
    y=styr;
    yz=styr;
    t_base=styr+(y-styr)*nseas-1;

//  ProgLabel_6.8.1  call fxn get_MGsetup()
    get_MGsetup();
    echoinput<<" initial bio parameters: "<<mgp_adj<<endl;

//  ProgLabel_6.8.2  call fxn get_growth1()
    get_growth1();
    VBK_seas=value(VBK_seas);
    wtlen_seas=value(wtlen_seas);
    CVLmin=value(CVLmin);
    CVLmax=value(CVLmax);
    echoinput<<"  growth factors that cannot be time-varying: "<<endl;
    echoinput<<"seasonal K: "<<VBK_seas<<endl;
    echoinput<<"seasonal wtlen adjuster "<<endl<<wtlen_seas<<endl;
    echoinput<<"CV_at_age_min(for each morph and gender) "<<CVLmin<<endl;
    echoinput<<"CV_at_age_max(for each morph and gender) "<<CVLmax<<endl;

//  ProgLabel_6.8.3  call fxn get_growth2()
    get_growth2();
    if(minL>10.0) {N_warn++; warning<<" Minimum size bin is:_"<<minL<<"; which is >10cm, which is large for use as size-at-age 0.0 recruitment"<<endl;}
    temp=Ave_Size(styr,0,1,nages);
    if(temp>0.95*len_bins(nlength)) {N_warn++; warning<<" Maximum size at age: "<<temp
    <<"; is within 5% of the largest size bin: "<<len_bins(nlength)<<"; Add more bins"<<endl;}

//  ProgLabel_6.8.4  call fxn get_natmort()
    get_natmort();
    natM = value(natM);
    surv1 = value(surv1);
    surv2 = value(surv2);
    echoinput<<" did natmort "<<endl<<natM<<endl;

//  ProgLabel_6.8.5  call fxn get_wtlen()
    get_wtlen();
    wt_len=value(wt_len);
    wt_len2=value(wt_len2);
    wt_len_fd=value(wt_len_fd);
    mat_len=value(mat_len);
    mat_fec_len=value(mat_fec_len);
    mat_age=value(mat_age);
    echoinput<<" did wtlen and maturity"<<endl;

    for(s=1;s<=nseas;s++)
    {
      t = styr+s-1;
      Make_AgeLength_Key();
      ALK(s) = value(ALK(s));
      ALK(smid) = value(ALK(smid));
      Save_Wt_Age(t)=Wt_Age_beg(s);
    }
    Wt_Age_beg=value(Wt_Age_beg);
    fec = value(fec);

//  ProgLabel_6.8.6  call fxn get_recr_distribution()
    get_recr_distribution();
    recr_dist_calc = value(recr_dist_calc);    //  so the just calculated constant values will be used unless its parms are active

//  ProgLabel_6.8.7  call fxn get_migration()
    if(do_migration>0)   // set up migration rates
    {
      get_migration();
      migrrate=value(migrrate);
    }

//  ProgLabel_6.8.8  call fxn get_age_age()
    if(N_ageerr>0)
    {
      AgeKey_StartAge=0;
      AgeKey_Linear1=1;
      AgeKey_Linear2=1;
      for(j=1;j<=N_ageerr;j++)
      {
        if(j!=Use_AgeKeyZero)
        {
          age_err(j)=age_err_rd(j);  //  this is an age err definition that has been read
        }
        else
        {
          AgeKey_StartAge=int(value(mgp_adj(AgeKeyParm)));
          if(mgp_adj(AgeKeyParm+3)==0.0000) {AgeKey_Linear1=1;} else {AgeKey_Linear1=0;}
          if(mgp_adj(AgeKeyParm+6)==0.0000) {AgeKey_Linear2=1;} else {AgeKey_Linear2=0;}
        }
        get_age_age(j,AgeKey_StartAge,AgeKey_Linear1,AgeKey_Linear2);  //  call function to get the age_age key
      }
      age_age=value(age_age);  //   because these are not based on parameters
    }
    echoinput<<" made the age_age' key "<<endl;
//  ProgLabel_6.8.9  calculated values have been set equal to value() to remove derivative info and save space if their parameters are held constant

//  ProgLabel_6.9  set up headers for ParmTrace
    if(Do_ParmTrace>0) ParmTrace<<"Phase Iter ObjFun Change SPB_start SPB_end BiasAdj_st BiasAdj_max BiasAdj_end ";
    if(Do_ParmTrace==1 || Do_ParmTrace==4)
    {
      for(i=1;i<=active_count;i++) {ParmTrace<<" "<<ParmLabel(active_parm(i));}
    }
    else if(Do_ParmTrace>=2)
    {
      for(i=1;i<=ParCount;i++) {ParmTrace<<" "<<ParmLabel(i);}
    }
    ParmTrace<<endl;

//  ProgLabel_6.10  Ready for estimation
    if(Turn_off_phase<0) {cout<<" Requested exit after read when turn_off_phase < 0 "<<endl; exit(1);}
    cout<<endl<<endl<<"Estimating...please wait..."<<endl;
    last_objfun=1.0e30;

// ****************************************************************************************************************
//  ProgLabel_7.0  PROCEDURE_SECTION
PROCEDURE_SECTION
  niter++;
  Mgmt_quant.initialize();
  Extra_Std.initialize();
  CrashPen.initialize();
  if(SzFreq_Nmeth>0) SzFreq_exp.initialize();

  if(mceval_phase() ) mceval_counter ++;   // increment the counter

  if(mcmcFlag==1)  //  so will do mcmc this run or is in mceval
  {
    if(Do_ParmTrace==1) Do_ParmTrace=4;  // to get all iterations
    if(Do_ParmTrace==2) Do_ParmTrace=3;  // to get all iterations
    if(mcmc_counter>10 || mceval_counter>10) Do_ParmTrace=0;
  }

//  ProgLabel_7.1  set up bias_adj vector
  sigmaR=SR_parm(N_SRparm2-3);
  two_sigmaRsq=2.0*sigmaR*sigmaR;
  half_sigmaRsq=0.5*sigmaR*sigmaR;
  if(initial_params::mc_phase==1) mcmc_counter++;

  biasadj.initialize();
  if(mcmcFlag==1)  //  so will do mcmc this run or is in mceval
  {
    biasadj_full=1.0;
  }
  else if(recdev_adj(5)<0.0)
  {
    biasadj_full=1.0;
  }
  else
  {
    for (y=styr-nages; y<=YrMax; y++)
    {
      if(y<recdev_first)  // before start of recrdevs
        {biasadj_full(y)=0.;}
      else if(y<=recdev_adj(1))
        {biasadj_full(y)=0.;}
      else if (y<=recdev_adj(2))
        {biasadj_full(y)=(y-recdev_adj(1)) / (recdev_adj(2)-recdev_adj(1))*recdev_adj(5);}
      else if (y<=recdev_adj(3))
        {biasadj_full(y)=recdev_adj(5);}   // max bias adjustment
      else if (y<=recdev_adj(4))
        {biasadj_full(y)=recdev_adj(5)-(y-recdev_adj(3)) / (recdev_adj(4)-recdev_adj(3))*recdev_adj(5);}
      else
        {biasadj_full(y)=0.;}
    }
  }

  if(SR_fxn==4 || do_recdev==0)
  {
    // keep all at 0.0 if not using SR fxn
  }
  else
  {
    if(recdev_do_early>0 && recdev_options(2)>=0 )    //  do logic on basis of recdev_options(2), which is read, not recdev_PH which can be reset to a neg. value
    {
      for(i=recdev_early_start;i<=recdev_early_end;i++)
      {biasadj(i)=biasadj_full(i);}
    }
    if(do_recdev>0 && recdev_PH_rd>=0 )
    {
      for(i=recdev_start;i<=recdev_end;i++)
      {biasadj(i)=biasadj_full(i);}
    }
    if(Do_Forecast>0 && recdev_options(3)>=0 )
    {
      for(i=recdev_end+1;i<=YrMax;i++)
      {biasadj(i)=biasadj_full(i);}
    }
    if(recdev_read>0)
    {
      for(j=1;j<=recdev_read;j++)
      {
        y=recdev_input(j,1);
        if(y>=recdev_first && y<=YrMax) biasadj(y)=biasadj_full(y);
      }
    }
  }
  if(do_once==1) echoinput<<"bias_adjustment vector: "<<biasadj<<endl;
 
  sd_offset_rec=sum(biasadj)*sd_offset;

//  ProgLabel_7.2  copy recdev parm vectors into full time series vector

  if(recdev_do_early>0) {recdev(recdev_early_start,recdev_early_end)=recdev_early(recdev_early_start,recdev_early_end);}
  if(do_recdev==1)
    {recdev(recdev_start,recdev_end)=recdev1(recdev_start,recdev_end);}
  else if(do_recdev==2)
    {recdev(recdev_start,recdev_end)=recdev2(recdev_start,recdev_end);}
  if(Do_Forecast>0) recdev(recdev_end+1,YrMax)=Fcast_recruitments(recdev_end+1,YrMax);  // only needed here for reporting

//  ProgLabel_7.3  reset Fmethod 2 to Fmethod 3 according to the phase
  if(F_Method==2)
  {
    F_Method_use=2;
    if(current_phase() < F_setup(2)) F_Method_use=3;  // use hybrid
  }
  else
  {
    F_Method_use=F_Method;
  }

//  ProgLabel_7.4  Do the calculations
  if(mceval_counter==0 || (mceval_counter>burn_intvl &&  ((double(mceval_counter)/double(thin_intvl)) - double((mceval_counter/thin_intvl))==0)  )) // check to see if burn in period is over
  {

    y=styr;
//  ProgLabel_7.4.1  call fxn get_initial_conditions()
    get_initial_conditions();
      if(do_once==1) cout<<" OK with initial conditions "<<endl;
//  ProgLabel_7.4.2  call fxn get_time_series()
    get_time_series();
      if(do_once==1) cout<<" OK with time series "<<endl;

//  ProgLabel_7.4.3  call fxn evaluate_the_objective_function()
    evaluate_the_objective_function();
    if(do_once==1) 
    {
      cout<<" OK with obj_func "<<obj_fun<<endl;
      do_once=0;
    }

//  ProgLabel_7.5  get quantities needed for forecasts
    if(Do_Forecast>0)
    {
//      if(save_for_report>0 || last_phase() || current_phase()==max_phase || ((sd_phase() || mceval_phase()) && (initial_params::mc_phase==0)))
      {
//  ProgLabel_7.5.1  calc average selectivity to use in forecast; store in endyr+1
        temp=float(Fcast_Sel_yr2-Fcast_Sel_yr1+1.);
        for(gg=1;gg<=gender;gg++)
        for(f=1;f<=Nfleet;f++)
        {
          tempvec_l.initialize();
          for(y=Fcast_Sel_yr1;y<=Fcast_Sel_yr2;y++) {tempvec_l+=sel_l(y,f,gg);}
          sel_l(endyr+1,f,gg)=tempvec_l/temp;

          tempvec_l.initialize();
          for(y=Fcast_Sel_yr1;y<=Fcast_Sel_yr2;y++) {tempvec_l+=sel_l_r(y,f,gg);}
          sel_l_r(endyr+1,f,gg)=tempvec_l/temp;

          tempvec_l.initialize();
          for(y=Fcast_Sel_yr1;y<=Fcast_Sel_yr2;y++) {tempvec_l+=discmort2(y,f,gg);}
          discmort2(endyr+1,f,gg)=tempvec_l/temp;

          tempvec_a.initialize();
          for(y=Fcast_Sel_yr1;y<=Fcast_Sel_yr2;y++) {tempvec_a+=sel_a(y,f,gg);}
          sel_a(endyr+1,f,gg)=tempvec_a/temp;
        }

//  ProgLabel_7.5.2  set-up relative F among fleets and seasons for forecast
        if(Fcast_RelF_Basis==1)  // set allocation according to range of years
        {
          temp=0.0;
          Fcast_RelF_Use.initialize();
          for(y=Fcast_RelF_yr1;y<=Fcast_RelF_yr2;y++)
          for(f=1;f<=Nfleet;f++)
          for(s=1;s<=nseas;s++)
          {
            t=styr+(y-styr)*nseas+s-1;
            Fcast_RelF_Use(s,f)+=Hrate(f,t);
          }
          temp=sum(Fcast_RelF_Use);
          if(temp==0.0)
          {
            Fcast_RelF_Use(1,1)=1.0;
            Fcurr_Fmult=0.0;
          }
          else
          {
            Fcast_RelF_Use/=temp;
            Fcurr_Fmult=temp/float(Fcast_RelF_yr2-Fcast_RelF_yr1+1);
          }
        }
        else  // Fcast_RelF_Basis==2 so set to values that were read
        {
          temp=0.0;
          for(f=1;f<=Nfleet;f++)
          for(s=1;s<=nseas;s++)
          {
            temp+=Fcast_RelF_Input(s,f);
          }
          Fcast_RelF_Use=Fcast_RelF_Input/temp;
          Fcurr_Fmult=temp;
        }
      }  //  end being in a phase for these calcs
    }  //  end getting quantities for forecasts

//  ProgLabel_7.5.3  calc average selectivity to use in benchmarks; store in styr-3
//  Bmark_Yr(1,6)<<" Benchmark years:  beg-end bio; beg-end selex; beg-end alloc"<<endl;

    if(Do_Benchmark>0)
    {
//      if(save_for_report>0 || last_phase() || current_phase()==max_phase || ((sd_phase() || mceval_phase()) && (initial_params::mc_phase==0)))
      {
    //  calc average body size to use in equil; store in styr-3
        temp=float(Bmark_Yr(2)-Bmark_Yr(1)+1.);  //  get denominator
        for(g=1;g<=gmorph;g++)
        if(use_morph(g)>0)
        {
          for(s=0;s<=nseas-1;s++)
          {
            tempvec_a.initialize();
            for(t=Bmark_t(1);t<=Bmark_t(2);t+=nseas) {tempvec_a+=Ave_Size(t+s,0,g);}
            Ave_Size(styr-3*nseas+s,0,g)=tempvec_a/temp;
            tempvec_a.initialize();
            for(t=Bmark_t(1);t<=Bmark_t(2);t+=nseas) {tempvec_a+=Ave_Size(t+s,1,g);}
            Ave_Size(styr-3*nseas+s,1,g)=tempvec_a/temp;
            for(f=0;f<=Nfleet;f++)
            {
              tempvec_a.initialize();
              for(t=Bmark_t(1);t<=Bmark_t(2);t+=nseas) {tempvec_a+=save_sel_fec(t+s,g,f);}
              save_sel_fec(styr-3*nseas+s,g,f)=tempvec_a/temp;
            }
          }
        }

        if(do_migration>0)
        {
          for(j=1;j<=do_migr2;j++)
          {
            tempvec_a.initialize();
            for(y=Bmark_Yr(1);y<=Bmark_Yr(2);y++) {tempvec_a+=migrrate(y,j);}
            migrrate(styr-3,j)=tempvec_a;
          }
        }

    //  calc average selectivity to use in equil; store in styr-1
        temp=float(Bmark_Yr(4)-Bmark_Yr(3)+1.);  //  get denominator
        for(gg=1;gg<=gender;gg++)
        for(f=1;f<=Nfleet;f++)
        {
          tempvec_l.initialize();
          for(y=Bmark_Yr(3);y<=Bmark_Yr(4);y++) {tempvec_l+=sel_l(y,f,gg);}
          sel_l(styr-3,f,gg)=tempvec_l/temp;

          tempvec_l.initialize();
          for(y=Bmark_Yr(3);y<=Bmark_Yr(4);y++) {tempvec_l+=sel_l_r(y,f,gg);}
          sel_l_r(styr-3,f,gg)=tempvec_l/temp;

          tempvec_l.initialize();
          for(y=Bmark_Yr(3);y<=Bmark_Yr(4);y++) {tempvec_l+=discmort2(y,f,gg);}
          discmort2(styr-3,f,gg)=tempvec_l/temp;

          tempvec_a.initialize();
          for(y=Bmark_Yr(3);y<=Bmark_Yr(4);y++) {tempvec_a+=sel_a(y,f,gg);}
          sel_a(styr-3,f,gg)=tempvec_a/temp;
        }

    //  set-up relative F among fleets and seasons
        if(Bmark_RelF_Basis==1)  // set allocation according to range of years
        {
          temp=0.0;
          Bmark_RelF_Use.initialize();
          for(y=Bmark_Yr(5);y<=Bmark_Yr(6);y++)
          for(f=1;f<=Nfleet;f++)
          for(s=1;s<=nseas;s++)
          {
            t=styr+(y-styr)*nseas+s-1;
            Bmark_RelF_Use(s,f)+=Hrate(f,t);
          }
          temp=sum(Bmark_RelF_Use);
          if(temp==0.0)
          {
            Bmark_RelF_Use(1,1)=1.0;
          }
          else
          {
          Bmark_RelF_Use/=temp;
          }
        }
        else  // Bmark_RelF_Basis==2 so set same as forecast
        {
          Bmark_RelF_Use=Fcast_RelF_Use;
        }
      }  //  end being in a phase for these calcs
    }  //  end getting quantities for benchmarks


//  ProgLabel_7.6  if sdphase or mcevalphase, do benchmarks and forecast and derived quantities
    if( (sd_phase() || mceval_phase()) && (initial_params::mc_phase==0))
    {

//  ProgLabel_7.6.1  call fxn Get_Benchmarks()
      if(Do_Benchmark>0)
      {
        Get_Benchmarks();
        did_MSY=1;
      }
      else
      {Mgmt_quant(1)=SPB_virgin;}

      if(mceval_phase()==0) {show_MSY=1;}

//  ProgLabel_7.6.2  call fxn Get_Forecast()
      if(Do_Forecast>0)
      {
        report5<<"THIS FORECAST FOR PURPOSES OF STD REPORTING"<<endl;
        Get_Forecast();
        did_MSY=1;
      }

//  ProgLabel_7.7  Call fxn Process_STDquant()
      Process_STDquant();
    }  // end of things to do in std_phase

//  ProgLabel_7.9  do screen output of procedure results
    if(current_phase() <= max_phase+1) phase_output(current_phase())=value(obj_fun);
    if(rundetail>1)
      {
       if(nobs_cr>0) cout<<" CPUE " <<surv_like<<endl;
       if(nobs_disc>0) cout<<" Disc " <<disc_like<<endl;
       if(nobs_mnwt>0) cout<<" MnWt " <<mnwt_like<<endl;
       if(Nobs_l_tot>0) cout<<" LEN  " <<length_like<<endl;
       if(Nobs_a_tot>0) cout<<" AGE  " <<age_like<<endl;
       if(nobs_ms_tot>0) cout<<" L-at-A  " <<sizeage_like<<endl;
       if(SzFreq_Nmeth>0) cout<<" sizefreq "<<SzFreq_like<<endl;
       if(Do_TG>0) cout<<" TG-fleetcomp "<<TG_like1<<endl<<" TG-negbin "<<TG_like2<<endl;
       cout<<" Recr " <<recr_like<<endl;
       cout<<" Parm_Priors " <<parm_like<<endl;
       cout<<" Parm_devs " <<parm_dev_like<<endl;
       cout<<" SoftBound "<<SoftBoundPen<<endl;
       cout<<" F_ballpark " <<F_ballpark_like<<endl;
       if(F_Method>1) {cout<<"Catch "<<catch_like;} else {cout<<"  crash "<<CrashPen;}
       cout<<" EQUL_catch " <<equ_catch_like<<endl;
      }
     if(rundetail>0)
     {
       temp=norm2(recdev(recdev_start,recdev_end));
       temp=sqrt((temp+0.0000001)/(double(recdev_end-recdev_start+1)));
     if(mcmc_counter==0 && mceval_counter==0)
     {cout<<current_phase()<<" "<<niter<<" -log(L): "<<obj_fun<<"  Spbio: "<<value(SPB_yr(styr))<<" "<<value(SPB_yr(endyr));}
     else if (mcmc_counter>0)
     {cout<<" MCMC: "<<mcmc_counter<<" -log(L): "<<obj_fun<<"  Spbio: "<<value(SPB_yr(styr))<<" "<<value(SPB_yr(endyr));}
     else if (mceval_counter>0)
     {cout<<" MCeval: "<<mceval_counter<<" -log(L): "<<obj_fun<<"  Spbio: "<<value(SPB_yr(styr))<<" "<<value(SPB_yr(endyr));}
       if(F_Method>1 && sum(catch_like)>0.01) {cout<<" cat "<<sum(catch_like);}
       else if (CrashPen>0.01) {cout<<"  crash "<<CrashPen;}
       cout<<endl;
     }
//  ProgLabel_7.10  write to ParmTrace
      if((Do_ParmTrace==1 && obj_fun<=last_objfun) || Do_ParmTrace==4)
      {
        ParmTrace<<current_phase()<<" "<<niter<<" "<<obj_fun<<" "<<obj_fun-last_objfun
        <<" "<<value(SPB_yr(styr))<<" "<<value(SPB_yr(endyr))<<" "<<biasadj(styr)<<" "<<max(biasadj)<<" "<<biasadj(endyr);
        for(j=1;j<=MGparm_PH.indexmax();j++)
        {
          if(MGparm_PH(j)>=0) {ParmTrace<<" "<<MGparm(j);}
        }
        if(MGparm_dev_PH>0 && N_MGparm_dev>0)
        {
          for(j=1;j<=N_MGparm_dev;j++)
          {ParmTrace<<MGparm_dev(j)<<" ";}
        }
        for(j=1;j<=SRvec_PH.indexmax();j++)
        {
          if(SRvec_PH(j)>=0) {ParmTrace<<" "<<SR_parm(j);}
        }
        if(recdev_cycle>0)
        {
          for(j=1;j<=recdev_cycle;j++)
          {
            if(recdev_cycle_PH(j)>=0) {ParmTrace<<" "<<recdev_cycle_parm(j);}
          }
        }
        if(recdev_early_PH>0) {ParmTrace<<" "<<recdev_early;}
        if(recdev_PH>0)
        {
          if(do_recdev==1) {ParmTrace<<" "<<recdev1;}
          if(do_recdev==2) {ParmTrace<<" "<<recdev2;}
        }
        if(Do_Forecast>0) ParmTrace<<Fcast_recruitments<<" ";
        if(Do_Forecast>0 && Do_Impl_Error>0) ParmTrace<<Fcast_impl_error<<" ";
        for(f=1;f<=Nfleet;f++)
        {
          if(init_F_PH(f)>0) {ParmTrace<<" "<<init_F(f);}
        }
        if(F_Method==2)    // continuous F
        {
          for(k=1;k<=F_rate.indexmax();k++)
          {
            if(have_catch(k)>0) {ParmTrace<<" "<<F_rate(k);}
          }
        }
        for(f=1;f<=Q_Npar;f++)
        {
          if(Q_parm_PH(f)>0) {ParmTrace<<" "<<Q_parm(f);}
        }
        for(k=1;k<=selparm_PH.indexmax();k++)
        {
          if(selparm_PH(k)>0) {ParmTrace<<" "<<selparm(k);}
        }
        if(selparm_dev_PH>0 && N_selparm_dev>0)
        {
          for(j=1;j<=N_selparm_dev;j++)
          {ParmTrace<<selparm_dev(j)<<" ";}
        }
        for(k=1;k<=TG_parm_PH.indexmax();k++)
        {
          if(TG_parm_PH(k)>0) {ParmTrace<<" "<<TG_parm(k);}
        }
        ParmTrace<<endl;
      }
      else if((Do_ParmTrace==2 && obj_fun<=last_objfun) || Do_ParmTrace==3)
      {
        ParmTrace<<current_phase()<<" "<<niter<<" "<<obj_fun<<" "<<obj_fun-last_objfun
        <<" "<<value(SPB_yr(styr))<<" "<<value(SPB_yr(endyr))<<" "<<biasadj(styr)<<" "<<max(biasadj)<<" "<<biasadj(endyr);
        ParmTrace<<" "<<MGparm<<" ";
        if(N_MGparm_dev>0) ParmTrace<<MGparm_dev<<" ";
        if(N_MGparm_dev>0)
        {
          for(j=1;j<=N_MGparm_dev;j++)
          {ParmTrace<<MGparm_dev(j)<<" ";}
        }
        ParmTrace<<SR_parm<<" ";
        if(recdev_cycle>0) ParmTrace<<recdev_cycle_parm;
        if(recdev_do_early>0) ParmTrace<<recdev_early<<" ";
        if(do_recdev==1) {ParmTrace<<recdev1<<" ";}
        if(do_recdev==2) {ParmTrace<<recdev2<<" ";}
        if(Do_Forecast>0) ParmTrace<<Fcast_recruitments<<" "<<Fcast_impl_error<<" ";
        ParmTrace<<init_F<<" ";
        if(F_Method==2) ParmTrace<<" "<<F_rate;
        if(Q_Npar>0) ParmTrace<<Q_parm<<" ";
        ParmTrace<<selparm<<" ";
        if(N_selparm_dev>0)
        {
          for(j=1;j<=N_selparm_dev;j++)
          {ParmTrace<<selparm_dev(j)<<" ";}
        }
        if(Do_TG>0) ParmTrace<<TG_parm<<" ";
        ParmTrace<<endl;
      }
      if(obj_fun<=last_objfun) last_objfun=obj_fun;
     docheckup=0;  // turn off reporting to echoinput.sso

//  ProgLabel_7.11  Call fxn get_posteriors if in mceval_phase
     if(mceval_phase()) get_posteriors();
    }
  if(mceval_phase() || initial_params::mc_phase==1)
  {
    No_Report=1;  //  flag to skip output reports after MCMC and McEVAL
  }
//  ProgLabel_7.12  end of PROCEDURE_SECTION

//********************************************************************
//  ProgLabel_8  RUNTIME_SECTION (not used in SS)
RUNTIME_SECTION
//  maximum_function_evaluations 200, 200, 200, 200, 200, 200, 200, 2000;
//  convergence_criteria 100, 10, 1, 0.1, 1e-4, 1e-4, 1e-4, 1e-4;

//  ProgLabel_9  TOP_OF_MAIN_SECTION
TOP_OF_MAIN_SECTION

//  ProgLabel_9.1  set array and gradient structure space
  arrmblsize = 200000000; // 2e8 = about 0.2 GB.
//  gradient_structure::set_GRADSTACK_BUFFER_SIZE(672647168);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(20000000); // 2e7 gets multiplied by the gradstack size (usually 48), to reach about 0.9 GB.
  //  gradstack is not allocated unless needed by the model, but arrmblsize and cmpdif are allocated immediately
//  gradient_structure::set_CMPDIF_BUFFER_SIZE(500000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000000); // 2e8 = about 0.2GB
  gradient_structure::set_MAX_NVAR_OFFSET(5000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);

//  ProgLabel_9.2  set clock start time
  time(&start); //this is to see how long it takes to run

//  ProgLabel_10.   GLOBALS_SECTION
GLOBALS_SECTION
  #include <admodel.h>
  #include <time.h>
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

//  ProgLabel_10.1  open output files using ofstream
  ofstream warning("warning.sso");
  ofstream echoinput("echoinput.sso");
  ofstream ParmTrace("ParmTrace.sso");
  ofstream report5("Forecast-report.sso");
  ofstream report2("CumReport.sso",ios::app);
  ofstream bodywtout("wtatage.ss_new");
  ofstream SS2out;   // this is just a create

//  ProgLabel_10.2  define some adstring variables
  adstring_array ParmLabel;  // extendable array to hold the parameter labels
  adstring_array fleetname;
  adstring depletion_basis_label;
  adstring F_report_label;
  adstring SPR_report_label;
  adstring onenum(4);
  adstring anystring;
  adstring anystring2;
  adstring version_info;
  adstring version_info_short;
  adstring_array Starter_Comments;
  adstring_array Data_Comments;
  adstring_array Control_Comments;
  adstring_array Forecast_Comments;
  random_number_generator radm(long(time(&start)));

//  ProgLabel_11.  BETWEEN_PHASES_SECTION
BETWEEN_PHASES_SECTION
  int j_phase=current_phase();  // this is the phase to come

//  ProgLabel_11.1 
  if(j_phase>1)
  {
    last_objfun=obj_fun;
  }

//  ProgLabel_11.2 for Fmethod=2, set parameter values (F_rate) equal to Hrate array calculated using hybrid method in previous phase
    if(F_Method==2)
    {
      if(F_setup(2)>1 && j_phase==F_setup(2) && readparfile==0)  //  so now start doing F as parameters
      {
        g=0;
        for(f=1;f<=Nfleet;f++)
        for(t=styr;t<=TimeMax;t++)
        {
          g++;
          F_rate(g)=Hrate(f,t);
        }
      }
    }
    
//  ProgLabel_12.  FINAL_SECTION
FINAL_SECTION
  int jj;
//  ProgLabel_12.1  get run ending time
  time(&finish);
  elapsed_time = difftime(finish,start);
  hour = long(elapsed_time)/3600;
  minute = long(elapsed_time)%3600/60;
  second = (long(elapsed_time)%3600)%60;
  cout<<endl<<"Finish time: "<<ctime(&finish);
  cout<<"Elapsed time: ";
  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds."<<endl;

  if(No_Report==1)
  {
    cout<<"MCMC finished; note: .sso and .ss_new files not produced after MCMC "<<endl;
  }

  else
  {
    cout<<"Final gradient: "<<objective_function_value::pobjfun->gmax << endl<<endl;
    if(objective_function_value::pobjfun->gmax >final_conv)
    {N_warn++; warning<<"Final gradient: "<<objective_function_value::pobjfun->gmax <<" is larger than final_conv: "<<final_conv<<endl;}

//  ProgLabel_12.2  Output the covariance matrix to covar.sso
    ofstream covarout("covar.sso");
    covarout<<version_info_short<<endl;
    covarout<<version_info<<endl;
    covarout<<"start_time: "<<ctime(&start)<<endl;
    covarout<<active_parms<<" "<<CoVar_Count<<endl;
    covarout<<"active-i active-j all-i all-j Par?-i Par?-j label-i label-j corr"<<endl;
    if(CoVar(1,1)==0.00 && CoVar(2,2)==0.0)
    {covarout<<"Variances are 0.0 for first two elements, so do not write "<<endl;}
    else
    {
      for(i=1;i<=CoVar_Count;i++)
      {
        covarout<<i<<" "<<0<<" "<<active_parm(i)<<" "<<active_parm(i);
        if(i<=active_parms) {covarout<<" Par ";} else {covarout<<" Der ";}
        covarout<<" Std "<<ParmLabel(active_parm(i))<<"   _   "<<CoVar(i,1)<<endl;
        for(j=2;j<=i;j++)
        {
          covarout<<i<<" "<<j-1<<" "<<active_parm(i)<<" "<<active_parm(j-1);
          if(i<=active_parms) {covarout<<" Par ";} else {covarout<<" Der ";}
          if((j-1)<=active_parms) {covarout<<" Par ";} else {covarout<<" Der ";}
          covarout<<ParmLabel(active_parm(i))<<" "<<ParmLabel(active_parm(j-1))<<" "<<CoVar(i,j)<<endl;
        }
      }
    }
    cout<<"Finished writing COVAR.SSO"<<endl;

    get_posteriors();
    
//  ProgLabel_12.3  go thru calculations again to get extra output quantities

//  ProgLabel_12.3.1  write out body weights to wtatage.ss_new.  Occurs while doing procedure with save_for_report=2 
    save_for_report=2;
    bodywtout<<1<<"  #_user_must_replace_this_value_with_number_of_lines_with_wtatage_below"<<endl;
    bodywtout<<N_WTage_maxage<<" # maxage"<<endl;
    bodywtout<<"# if yr=-yr, then fill remaining years for that seas, growpattern, gender, fleet"<<endl;
    bodywtout<<"# fleet 0 contains begin season pop WT"<<endl;
    bodywtout<<"# fleet -1 contains mid season pop WT"<<endl;
    bodywtout<<"# fleet -2 contains maturity*fecundity"<<endl;
    bodywtout<<"#yr seas gender growpattern birthseas fleet "<<age_vector<<endl;
    save_gparm=0;
    y=styr;
    CrashPen.initialize();
    get_initial_conditions();
    get_time_series();
    bodywtout.close();
//    evaluate_the_objective_function();

//  ProgLabel_12.3.2  set save_for_report=1 then call initial_conditions and time_series to get other output quantities
    save_for_report=1;
    CrashPen.initialize();
    save_gparm=0;
    y=styr;
    get_initial_conditions();
    get_time_series();
    evaluate_the_objective_function();

//  ProgLabel_12.3.3  do benchmarks and forecast and stdquantities with save_for_report=1
    if(mceval_phase()==0) {show_MSY=1;}
    if(Do_Benchmark>0)
    {
      if(did_MSY==0) Get_Benchmarks();
    }
    else
    {Mgmt_quant(1)=SPB_virgin;}
    if(Do_Forecast>0)
    {
      report5<<"THIS FORECAST FOR PURPOSES OF GETTING DISPLAY QUANTITIES"<<endl;
      Get_Forecast();
    }
    Process_STDquant();

//  ProgLabel_12.4  call fxn write_bigoutput()
    write_bigoutput();
    cout<<" finished big report "<<endl;

//  ProgLabel_12.5  call fxn write_summaryoutput()
    if(Do_CumReport>0) write_summaryoutput();
    cout<<" finished summary report "<<endl;

//  ProgLabel_12.6  call fxn write_rebuilder_output
    if(Do_Rebuilder>0 && mceval_counter<=1) write_rebuilder_output();
    cout<<" finished rebuilder report "<<endl;
//  ProgLabel_12.7  call fxn write_nudata()
    write_nudata();
    cout<<" finished nudata report "<<endl;
//  ProgLabel_12.8  call fxn write_nucontrol()
    write_nucontrol();
    cout<<" finished nucontrol report "<<endl;
//  ProgLabel_12.9  call fxn write_Bzero_output()
    write_Bzero_output();
    warning<<" N warnings: "<<N_warn<<endl;
    if(MG_adjust_method==3) warning<<"time-vary MGparms not bound checked"<<endl;
    if(selparm_adjust_method==3) warning<<"time-vary selparms not bound checked"<<endl;
    cout<<endl<<"!!  Run has completed  !!            ";
    
//  Proglabel_12.10  finish up with final writes to warning.sso
    if(N_changed_lambdas>0)
      {N_warn++; warning<<"Reminder: Number of lamdas !=0.0 and !=1.0:  "<<N_changed_lambdas<<endl; }
    warning<<"Number_of_active_parameters_on_or_near_bounds: "<<Nparm_on_bound<<endl;
    if(Nparm_on_bound>0) N_warn++;
    if(N_warn>0)
    {cout<<"See warning.sso for N warnings: "<<N_warn<<endl;}
    else
    {cout<<"No warnings"<<endl;}
  }

//  ProgLabel_13.  REPORT_SECTION
REPORT_SECTION
//  ProgLabel_13.1  write limited output to SS3.rep
  if(nobs_cr>0) report<<" CPUE " <<surv_like<<endl;
  if(nobs_disc>0) report<<" Disc " <<disc_like<<endl;
  if(nobs_mnwt>0) report<<" MnWt " <<mnwt_like<<endl;
  if(Nobs_l_tot>0) report<<" LEN  " <<length_like<<endl;
  if(Nobs_a_tot>0) report<<" AGE  " <<age_like<<endl;
  if(nobs_ms_tot>0) report<<" L-at-A  " <<sizeage_like<<endl;
  report<<" EQUL " <<equ_catch_like<<endl;
  report<<" Recr " <<recr_like<<endl;
  report<<" Parm " <<parm_like<<endl;
  report<<" F_ballpark " <<F_ballpark_like<<endl;
  if(F_Method>1) {report<<"Catch "<<catch_like<<endl;} else {report<<"  crash "<<CrashPen<<endl;}
  if(SzFreq_Nmeth>0) report<<" sizefreq "<<SzFreq_like<<endl;
  if(Do_TG>0) report<<" TG-fleetcomp "<<TG_like1<<endl<<" TG-negbin "<<TG_like2<<endl;
  report<<" -log(L): "<<obj_fun<<"  Spbio: "<<value(SPB_yr(styr))<<
  " "<<value(SPB_yr(endyr))<<endl;

  report<<endl<<"Year Spbio Recruitment"<<endl;
  report<<"Virg "<<SPB_yr(styr-2)<<" "<<exp_rec(styr-2,4)<<endl;
  report<<"Init "<<SPB_yr(styr-1)<<" "<<exp_rec(styr-1,4)<<endl;
  for(y=styr;y<=endyr;y++) report<<y<<" "<<SPB_yr(y)<<" "<<exp_rec(y,4)<<endl;

  report<<endl<<"EXPLOITATION F_Method: ";
  if(F_Method==1) {report<<" Pope's_approx ";} else {report<<" instantaneous_annual_F ";}
  report<<endl<<"X Catch_Units ";
  for(f=1;f<=Nfleet;f++) if(catchunits(f)==1) {report<<" Bio ";} else {report<<" Num ";}
  report<<endl<<"yr seas"; for(f=1;f<=Nfleet;f++) report<<" "<<f;
  report<<endl<<"init_yr 1 "<<init_F<<endl;
  for(y=styr;y<=endyr;y++)
  for(s=1;s<=nseas;s++)
  {
    t=styr+(y-styr)*nseas+s-1;
    report<<y<<" "<<s<<" "<<column(Hrate,t)<<endl;
  }

  report<<endl<< "LEN_SELEX" << endl;
  report<<"Fleet gender "<<len_bins_m<<endl;
  for (f=1;f<=Ntypes;f++)
  {
    if(seltype(f,1)>0)
    {
      for (gg=1;gg<=gender;gg++) report<<f<<"-"<<fleetname(f)<<gg<<" "<<sel_l(endyr,f,gg)<<endl;
    }
  }

  report<<endl<< "AGE_SELEX" << endl;
  report<<"Fleet gender "<<age_vector<<endl;
  for (f=1;f<=Ntypes;f++)
  {
    if(seltype(f+Ntypes,1)>10)
    {
      for (gg=1;gg<=gender;gg++) report<<f<<"-"<<fleetname(f)<<" "<<gg<<" "<<sel_a(endyr,f,gg)<<endl;
    }
  }


//  ProgLabel_13.2  call fxn write_bigoutput() as last_phase finishes and before doing Hessian
    if(last_phase()>0)
    {
    save_for_report=1;
    CrashPen.initialize();
    save_gparm=0;
    y=styr;
    get_initial_conditions();

    get_time_series();
    evaluate_the_objective_function();
    write_bigoutput();
    save_for_report=0;
    SS2out.close();
    }

//*********************************************************************
 /*  ProgLabel_FUNCTION 14 get_MGsetup  :  apply time-varying factors this year to the MG parameters to create mgp_adj vector */
FUNCTION void get_MGsetup()
  {
  mgp_adj=MGparm;
  int y1;
//  ProgLabel_14.1  calculate any trends that will be needed for any of the MG parameters
  if(N_MGparm_trend>0)
  {
    for(f=1;f<=N_MGparm_trend;f++)
    {
      j=MGparm_trend_rev(f);  //  parameter affected
      k=MGparm_trend_rev_1(f);  // base index for trend parameters
      if(y==styr)
      {
        //  calc endyr value, but use logistic transform to keep with bounds of the base parameter
        if(MGparm_1(j,13)==-1)
        {
          temp=log((MGparm_1(j,2)-MGparm_1(j,1)+0.0000002)/(MGparm(j)-MGparm_1(j,1)+0.0000001)-1.)/(-2.);   // transform the base parameter
          temp+=MGparm(k+1);     //  add the offset  Note that offset value is in the transform space
          temp1=MGparm_1(j,1)+(MGparm_1(j,2)-MGparm_1(j,1))/(1.+mfexp(-2.*temp));   // backtransform
        }
        else if(MGparm_1(j,13)==-2)
        {
          temp1=MGparm(k+1);  // set ending value directly
        }

        if(MGparm_HI(k+2)<=1.1)  // use max bound as switch
        {temp3=r_years(styr)+MGparm(k+2)*(r_years(endyr)-r_years(styr));}  // infl year
        else
        {temp3=MGparm(k+2);}

        temp2=cumd_norm((r_years(styr)-temp3)/MGparm(k+3));     //  cum_norm at styr
        temp=(temp1-MGparm(j)) / (cumd_norm((r_years(endyr)-temp3)/MGparm(k+3))-temp2);   //  delta in cum_norm between styr and endyr
        for(int y1=styr;y1<=YrMax;y1++)
        {
          if(y1<=endyr)
          {MGparm_trend(f,y1)=MGparm(j) + temp * (cumd_norm((r_years(y1)-temp3)/MGparm(k+3) )-temp2);}
          else
          {MGparm_trend(f,y1)=MGparm_trend(f,y1-1);}
        }
      }
      mgp_adj(j)=MGparm_trend(MGparm_trend_point(j),y);
    }
  }

//  ProgLabel_14.2  else create MGparm block values
  else if (N_MGparm_blk>0)
  {
    for(j=1;j<=N_MGparm;j++)
    {
      z=MGparm_1(j,13);    // specified block pattern
      if(z>0)  // uses blocks
      {
        if(y==styr)  // set up the block values time series
        {
          g=1;
          if(MGparm_1(j,14)<3)
          {
            for(a=1;a<=Nblk(z);a++)
            {
              for(int y1=Block_Design(z,g);y1<=Block_Design(z,g+1);y1++)  // loop years for this block
              {
                k=Block_Defs_MG(j,y1);  // identifies parameter that holds the block effect
                MGparm_block_val(j,y1)=MGparm(k);
              }
              g+=2;
            }
          }
          else
          {
            temp=0.0;
            for(a=1;a<=Nblk(z);a++)
            {
              y1=Block_Design(z,g);   // first year of block
              k=Block_Defs_MG(j,y1);  // identifies parameter that holds the block effect
              temp+=MGparm(k);  // increment by the block delta
              for(int y1=Block_Design(z,g);y1<=Block_Design(z,g+1);y1++)  // loop years for this block
              {
                MGparm_block_val(j,y1)=temp;
              }
              g+=2;
            }
          }
        }  // end block setup
      }  // end uses blocks
    }  // end parameter loop
  }  // end block section

//  ProgLabel_14.3  create MGparm dev randwalks if needed
  if(N_MGparm_dev>0 && y==styr)
  {
    for(k=1;k<=N_MGparm_dev;k++)
    {
      if(MGparm_dev_type(k)==3)
      {
        MGparm_dev_rwalk(k,MGparm_dev_minyr(k))=MGparm_dev(k,MGparm_dev_minyr(k));
        j=MGparm_dev_minyr(k);
        for (j=MGparm_dev_minyr(k)+1;j<=MGparm_dev_maxyr(k);j++)
        {
          MGparm_dev_rwalk(k,j)=MGparm_dev_rwalk(k,j-1)+MGparm_dev(k,j);
        }
      }
    }
  }

//  ProgLabel_14.4  switch(MG_adjust_method)
  switch(MG_adjust_method)
  {
    case 3:
    {
      //  no break statement, so will execute case 1 code
    }
//  ProgLabel_14.4.1   standard MG_adjust_method (1 or 3), loop MGparms
    case 1:
    {
      for(f=1;f<=N_MGparm;f++)
      {
//  ProgLabel_14.4.1.1  adjust for blocks
        if(MGparm_1(f,13)>0)   // blocks
        {
          if(Block_Defs_MG(f,yz)>0)
          {
            if(MGparm_1(f,14)==0)
              {mgp_adj(f) *= mfexp(MGparm_block_val(f,yz));}
            else if(MGparm_1(f,14)==1)
              {mgp_adj(f) += MGparm_block_val(f,yz);}
            else if(MGparm_1(f,14)==2)
              {mgp_adj(f) = MGparm_block_val(f,yz);}
            else if(MGparm_1(f,14)==3)  // additive based on delta approach
              {mgp_adj(f) += MGparm_block_val(f,yz);}
          }
        }

//  ProgLabel_14.4.1.2  adjust for env linkage
        if(MGparm_env(f)>0)
        {
          if(MGparm_envtype(f)==1)
          {mgp_adj(f)*=mfexp(MGparm(MGparm_env(f))*env_data(yz,MGparm_envuse(f)));}
          else
          {mgp_adj(f)+=MGparm(MGparm_env(f))*env_data(yz,MGparm_envuse(f));}
        }

//  ProgLabel_14.4.1.3  adjust for Annual deviations
        k=MGparm_dev_select(f);
        if(k>0)
        {
          if(yz>=MGparm_dev_minyr(k) && yz<=MGparm_dev_maxyr(k))
          {
            if(MGparm_dev_type(k)==1)  // multiplicative
            {mgp_adj(f) *= mfexp(MGparm_dev(k,yz));}
            else if(MGparm_dev_type(k)==2)  // additive
            {mgp_adj(f) += MGparm_dev(k,yz);}
            else if(MGparm_dev_type(k)==3)  // additive rwalk
            {mgp_adj(f) += MGparm_dev_rwalk(k,yz);}
          }
        }

//  ProgLabel_14.4.1.4  do bound check if MG_adjust_method=1
        if(MG_adjust_method==1 && (save_for_report==1 || do_once==1))  // so does not check bounds if MG_adjust_method==3
        {
          if(mgp_adj(f)<MGparm_1(f,1) || mgp_adj(f)>MGparm_1(f,2))
          {
            N_warn++;
            warning<<" adjusted MGparm out of bounds (parm#, yr, min, max, base, adj_value) "<<f<<" "<<yz<<" "<<
            MGparm_1(f,1)<<" "<<MGparm_1(f,2)<<" "<<MGparm(f)<<" "<<mgp_adj(f)<<" "<<ParmLabel(f)<<endl;
          }
        }
      }  // end parameter loop (f)
      break;
    }

//  ProgLabel_14.4.2   constrained MG_adjust_method (2), loop MGparms
    case 2:
    {
      for(f=1;f<=N_MGparm;f++)
      {
        j=0;
        temp=log((MGparm_HI(f)-MGparm_LO(f)+0.0000002)/(mgp_adj(f)-MGparm_LO(f)+0.0000001)-1.)/(-2.);   // transform the parameter

//  ProgLabel_14.4.2.1  adjust for blocks
        if(MGparm_1(f,13)>0)   // blocks
        {
          if(Block_Defs_MG(f,yz)>0)
          {
            j=1;  //  change is being made
            if(MGparm_1(f,14)==1)
              {temp+=MGparm_block_val(f,yz);}
            else if(MGparm_1(f,14)==2)  // block as replacement
              {temp=log((MGparm_HI(f)-MGparm_LO(f)+0.0000002)/(MGparm_block_val(f,yz)-MGparm_LO(f)+0.0000001)-1.)/(-2.);}
            else if(MGparm_1(f,14)==3)  // additive based on delta approach
              {temp += MGparm_block_val(f,yz);}
          }
        }

//  ProgLabel_14.4.2.2  adjust for env linkage
        if(MGparm_env(f)>0)  //  do environmental effect;  only additive allowed for adjustment method=2
        {j=1; temp+=MGparm(MGparm_env(f))* env_data(yz,MGparm_envuse(f));}

//  ProgLabel_14.4.2.3  adjust for annual deviations
        k=MGparm_dev_select(f);
        if(k>0)
        {
          if(yz>=MGparm_dev_minyr(k) && yz<=MGparm_dev_maxyr(k))
            {
              j=1;
              if(MGparm_dev_type(k)==2)
              {temp += MGparm_dev(k,yz);}
              else if(MGparm_dev_type(k)==3)
              {temp += MGparm_dev_rwalk(k,yz);}  // note that only additive effect is allowed
            }
        }
        if(j==1) mgp_adj(f)=MGparm_LO(f)+(MGparm_HI(f)-MGparm_LO(f))/(1.+mfexp(-2.*temp));   // backtransform
      }  // end parameter loop (f)
      break;
    }  // end case 2
  }   // end switch method

// if direct estimation (no offsets) is used; then look for M and growth parameters
// for which male value is set equal to female value

//  ProgLabel_14.5  if MGparm method =1 (no offsets), then do direct assignment if parm value is 0.0.
  if(MGparm_def==1)
  {
    for(j=1;j<=N_GP*gender*MGparm_per_def;j++)
    {
      if(MGparm_offset(j)>0) mgp_adj(j) = mgp_adj(MGparm_offset(j));
    }
  }
  if(save_for_report>0) mgp_save(yz)=value(mgp_adj);
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 15 get_growth1 */
FUNCTION void get_growth1()
  {
//  ProgLabel_15.1  create seasonal effects for growth K, and for wt_len parameters
    if(MGparm_doseas>0)
    {
      if(MGparm_seas_effects(10)>0)  // for seasonal K
      {
        VBK_seas(0)=0.0;
        for(s=1;s<=nseas;s++)
        {
          VBK_seas(s)=mfexp(MGparm(MGparm_seas_effects(10)+s));
          VBK_seas(0)+=VBK_seas(s)*seasdur(s);
        }
      }
      else
      {
        VBK_seas=sum(seasdur);  // set vector to null effect
      }

      for(j=1;j<=8;j++)
      {
        if(MGparm_seas_effects(j)>0)
        {
          wtlen_seas(0,j)=0.0;
          for(s=1;s<=nseas;s++)
          {
            wtlen_seas(s,j)=mfexp(MGparm(MGparm_seas_effects(j)+s));
            wtlen_seas(0,j)+=wtlen_seas(s,j)*seasdur(s);
          }
        }
        else
        {
          for(s=0;s<=nseas;s++) {wtlen_seas(s,j)=1.0;}
        }
      }
    }
    else
    {
      VBK_seas=sum(seasdur);  // set vector to null effect
      wtlen_seas=1.0;  // set vector to null effect
    }

//  ProgLabel_15.2  create variability of size-at-age factors using direct assignment or offset approaches
      Ip=1-MGparm_per_def;
      for (gp=1;gp<=N_GP*gender;gp++)
        {
        Ip+=MGparm_per_def;
        j=Ip+MGparm_per_def-2;  // index for CVmin
        k=j+1;  // index for CVmax
        switch(MGparm_def)    // for CV of size-at-age
          {
            case 1:  // direct
            {
            if(MGparm(j)>0)
              {CVLmin(gp)=MGparm(j);} else {CVLmin(gp)=MGparm(MGparm_per_def-1);}
            if(MGparm(k)>0)
              {CVLmax(gp)=MGparm(k);} else {CVLmax(gp)=MGparm(MGparm_per_def);}
            break;
            }
            case 2:  // offset
            {
            if(gp==1)
              {CVLmin(gp)=MGparm(j); CVLmax(gp)=MGparm(k);}
            else
              {CVLmin(gp)=CVLmin(1)*mfexp(MGparm(j)); CVLmax(gp)=CVLmax(1)*mfexp(MGparm(k));}
            break;
            }
            case 3:  // offset like SS2 V1.23
            {
            if(gp==1)
              {CVLmin(gp)=MGparm(j); CVLmax(gp)=CVLmin(1)*mfexp(MGparm(k));}
            else
              {CVLmin(gp)=CVLmin(1)*mfexp(MGparm(j)); CVLmax(gp)=CVLmin(gp)*mfexp(MGparm(k));}
            break;
            }
          }  // end switch
          if((CVLmin(gp)!=CVLmax(gp)) || active(MGparm(MGparm_per_def)) || active(MGparm(k)))
          {CV_const(gp)=1;} else {CV_const(gp)=0;}
        }
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 16 get_growth2 */
FUNCTION void get_growth2()
  {
//  progress mean growth through time series, accounting for seasonality and possible change in parameters
//   get mean size at the beginning and mid-point of the season
    dvariable grow;
    dvariable growhalf;
    dvar_vector CV_delta(1,N_GP*gender);
    dvariable LminR;
    dvariable LmaxR;
    dvariable LinfR;
    dvariable inv_Richards;

    int k2;
//  ProgLabel_16.1  create Cohort_Growth offset for the cohort borne (age 0) this year
      if(CGD>0)   //  cohort specific growth multiplier
      {
        temp=mgp_adj(MGP_CGD);
        k=min(nages,(YrMax-y));
        for(a=0;a<=k;a++) {Cohort_Growth(y+a,a)=temp;}
      }

      Ip=N_natMparms+1;  // first growth parameter
      for (gp=1;gp<=N_GP*gender;gp++)
      {
//  ProgLabel_16.2  loop growth patterns (N_GP * gender)
//  ProgLabel_16.2.1  set Lmin, Lmax, VBK, Richards to this year's values for mgp_adj
        if(MGparm_def>1 && gp>1)   // switch for growth parms
        {
          Lmin(gp)=Lmin(1)*mfexp(mgp_adj(Ip));
          Lmax_temp(gp)=Lmax_temp(1)*mfexp(mgp_adj(Ip+1));
          VBK(gp)=VBK(1)*mfexp(mgp_adj(Ip+2));  //  assigns to all ages
          if(Grow_type==2) Richards(gp)=Richards(1)*mfexp(mgp_adj(Ip+3));
        }
        else
        {
          Lmin(gp)=mgp_adj(Ip);
          Lmax_temp(gp)=mgp_adj(Ip+1);
          VBK(gp)=-mgp_adj(Ip+2);  // because always used as negative; assigns to all ages
          if(Grow_type==2) Richards(gp)=mgp_adj(Ip+3);
        }
        
        if(Grow_type==3)  //  age specific k
        {
          j=1;
          for(a=1;a<=nages;a++)
          {
            if(a==Age_K_points(j))
            {
              VBK(gp,a)=VBK(gp,a-1)*mgp_adj(Ip+2+j);
              if(j<Age_K_count) j++;
            }
            else
            {
              VBK(gp,a)=VBK(gp,a-1);
            }
          }
        }

        k=min(nages,(YrMax-y));
        if(y==styr)
        {
          Cohort_Lmin(gp)=Lmin(gp);
        }
        else
        {
          for(a=0;a<=k;a++) {Cohort_Lmin(gp,y+a,a)=Lmin(gp);}
        }

        if(Grow_type!=2)  //  not Richards
        {
          if(AFIX2==999)
          {L_inf(gp)=Lmax_temp(gp);}
          else
          {
            L_inf(gp)=Lmin(gp)+(Lmax_temp(gp)-Lmin(gp))/(1.-mfexp(VBK(gp,nages)*VBK_seas(0)*(AFIX_delta)));
          }
        }
        else   //  Richards
        {
          inv_Richards=1.0/Richards(gp);
          LminR=pow(Lmin(gp),Richards(gp));
          if(AFIX2==999)
          {
            L_inf(gp)=Lmax_temp(gp);
          }
          else
          {
            LmaxR=pow(Lmax_temp(gp), Richards(gp));
            temp=LminR+(LmaxR-LminR)/(1.-mfexp(VBK(gp,nages)*VBK_seas(0)*(AFIX_delta)));
            L_inf(gp)=pow(temp,inv_Richards);
          }
          LinfR=pow(L_inf(gp),Richards(gp));
        }

//  ProgLabel_16.2.2  set up the delta in growth variability across ages if needed
        if(CV_const(gp)>0)
        {
          if(CV_depvar_a==0)
            {CV_delta(gp)=(CVLmax(gp)-CVLmin(gp))/(Lmax_temp(gp)-Lmin(gp));}
          else
            {CV_delta(gp)=(CVLmax(gp)-CVLmin(gp))/(AFIX2_forCV-AFIX);}
        }
        else
        {
          CV_delta(gp)=0.0;
          CV_G(gp)=CVLmin(gp);  // sets all seasons and whole age range
        }

        g=(gp-1)*nseas*N_platoon+int(N_platoon/2)+1-N_platoon;  // find the mid-morph being processed
//   note that an index may have been created in data section for the above statement
        for(birthseas=1;birthseas<=nseas;birthseas++)
        {
//  ProgLabel_16.2.3  loop birth seasons (really looping more morphs)
          g+=N_platoon;
          if(use_morph(g)>0)
          {

//  ProgLabel_16.2.3.1  if y=styr, get size-at-age in first season of this first year
            if(y==styr)
            {
              if(docheckup==1) echoinput<<" growth in start year for platoon: "<<g<<" "<<gp<<" Lmin: "<<Lmin(gp)<<" Linf: "<<L_inf(gp)<<" VBK: "<<VBK(gp)<<endl<<endl;
              switch(Grow_type)
              {
                default:
                {
//                  Ave_Size(styr,0,g) = L_inf(gp) + (Lmin(gp)-L_inf(gp))*mfexp(VBK(gp,nages)*VBK_seas(0)*(curr_age_beg(g,1)(0,nages)-AFIX));
                  Ave_Size(styr,0,g,0)=L_inf(gp) + (Lmin(gp)-L_inf(gp))*mfexp(VBK(gp,0)*VBK_seas(0)*(curr_age_beg(g,1,0)-AFIX));
                  for(a=1;a<=nages;a++)
                  {
                  //  error below prior to 3.24k.  must be >=-1.0 for the young ages
                    if(lin_grow_beg(g,1,a)>=-1.0)  // younger or equal to AFIX, so use youngest age groups VBK
                    {
                      Ave_Size(styr,0,g,a) = Lmin(gp) + (Lmin(gp)-L_inf(gp))* (mfexp(VBK(gp,0)*VBK_seas(0)*(curr_age_beg(g,1,a)-AFIX))-1.0);
                    }
                    else
                    {
                      Ave_Size(styr,0,g,a) = Ave_Size(styr,0,g,a-1) + (Ave_Size(styr,0,g,a-1)-L_inf(gp))* (mfexp(VBK(gp,a-1)*VBK_seas(0))-1.0);
                    }
                  }  // done ageloop
                  if(birthseas>1) Ave_Size(styr,0,g,0)=len_bins(1);
                  if(docheckup==1) echoinput<<" size_at_age w/o linear ramp overlay: "<<endl<<Ave_Size(styr,0,g)<<endl;
                  break;
                }
                case 2:  //  Richards; age-specific VBK not allowed here
                {
                  for(a=0;a<=nages;a++)
                  {
                    if(lin_grow_beg(g,1,a)<0.0)  // so in growth curve
                    {
                      temp=LinfR + (LminR-LinfR)*mfexp(VBK(gp,nages)*VBK_seas(0)*(curr_age_beg(g,1,a)-AFIX));
                      Ave_Size(styr,0,g,a) = pow(temp,inv_Richards);
                    }
                    else
                    {Ave_Size(styr,0,g,a) = len_bins(1);}  //  will be replaced later when doing the linear portion
                  }
                  break;
                }
              }

//  ProgLabel_16.2.3.1.1  calc approximation to mean size at maxage to account for growth after reaching the maxage (accumulator age)
              temp=0.0;
              temp1=0.0;
              temp2=mfexp(-0.2);  //  cannot use natM or Z because growth is calculated first
              temp3=L_inf(gp)-Ave_Size(styr,0,g,nages);  // delta between linf and the size at nages
              //  frac_ages = age/nages, so is fraction of a lifetime
              temp4=1.0;
              for(a=0;a<=nages;a++)
              {
                temp+=temp4*(Ave_Size(styr,0,g,nages)+frac_ages(a)*temp3);  // so grows linearly from size at nages to size at nages+nages
                temp1+=temp4;   //  accumulate numbers to create denominator for mean size calculation
                temp4*=temp2;  //  decay numbers at age by exp(-0.2)
              }
              Ave_Size(styr,0,g,nages)=temp/temp1;  //  this is weighted mean size at nages
              if(do_once==1) echoinput<<"Adjusted size at maxage: "<<Ave_Size(styr,0,g,nages)<<endl;
            }  //  end initial year calcs

            int add_age;
            for (s=1;s<=nseas;s++)
            {
//  ProgLabel_16.2.3.2  loop seasons
              smid=s+nseas;
              t=t_base+s;
              if(s==nseas) add_age=1; else add_age=0;   //      advance age or not

// growth to mid-period and next period
              switch(Grow_type)
              {
                default:
                {
//  ProgLabel_16.2.3.2.1  standard von Bert growth, loop ages
                  for(a=0;a<=nages;a++)
                  {
//  ProgLabel_16.2.3.2.1.1  calc midseas size-at-age from begin season size-at-age, accounting for transition from linear to von Bert as necessary
                    if(lin_grow_beg(g,s,a)>=0.0)  // in linear phase for begin seas
                    {
                      Ave_Size(t,0,g,a) = len_bins(1)+lin_grow_beg(g,s,a)*(Cohort_Lmin(gp,y,a)-len_bins(1));  // begin of season
                    }
                    if(lin_grow_mid(g,s,a)>=0.0)  // in linear phase for midseas
                    {
                      Ave_Size(t,1,g,a) = len_bins(1)+lin_grow_mid(g,s,a)*(Cohort_Lmin(gp,y,a)-len_bins(1));  // mid season
                    }

                    growhalf = (mfexp(VBK(gp,a)*seasdur_half(s)*VBK_seas(s))-1.0);
// NOTE:  there is no seasonal interpolation, or real age adjustment for age-specific K.  Maybe someday....
                    if(lin_grow_beg(g,s,a)==-1.0)  // first time point beyond AFIX
                    {
                      Ave_Size(t,0,g,a) = Cohort_Lmin(gp,y,a) + (Cohort_Lmin(gp,y,a)-L_inf(gp))* (mfexp(VBK(gp,a)*(curr_age_beg(g,s,a)-AFIX)*VBK_seas(s))-1.0)*Cohort_Growth(y,a);
                    }

                    if(lin_grow_mid(g,s,a)<0.0)  //  so doing growth curve by mid of this seas, so project with growth to begin of next
                    {
                      if(lin_grow_mid(g,s,a)==-1.0)  // first time point beyond AFIX
                      {
                        Ave_Size(t,1,g,a) = Cohort_Lmin(gp,y,a) + (Cohort_Lmin(gp,y,a)-L_inf(gp))* (mfexp(VBK(gp,a)*(curr_age_mid(g,s,a)-AFIX)*VBK_seas(s))-1.0)*Cohort_Growth(y,a);
                      }
                      else  //  calc midseas size starting from size at beginning of this season
                      {
                        t2=Ave_Size(t,0,g,a)-L_inf(gp);  //  remaining growth potential
                        if(time_vary_MG(y,2)>0 && t2>-1.)
                        {
                          join1=1.0/(1.0+mfexp(-(50.*t2/(1.0+fabs(t2)))));  //  note the logit transform is not perfect, so growth near Linf will not be exactly same as with native growth function
                          t2*=(1.-join1);  // trap to prevent decrease in size-at-age
                        }
                        Ave_Size(t,1,g,a) = Ave_Size(t,0,g,a)+growhalf * t2*Cohort_Growth(y,a);
                      }

//  ProgLabel_16.2.3.2.1.2  calc size at end of the season, which will be size at begin of next season using current seasons growth parms
                      if(a<nages) {k2=a+add_age;} else {k2=a;}  // where add_age =1 if s=nseas, else 0
                      t2=Ave_Size(t,1,g,a)-L_inf(gp);  //  remaining growth potential
                      if(time_vary_MG(y,2)>0 && t2>-1.)
                      {
                        join1=1.0/(1.0+mfexp(-(50.*t2/(1.0+fabs(t2)))));  //  note the logit transform is not perfect, so growth near Linf will not be exactly same as with native growth function
                        t2*=(1.-join1);  // trap to prevent decrease in size-at-age
                      }
                      temp1 = Ave_Size(t,1,g,a)+growhalf * t2*Cohort_Growth(y,a);
                      if((a<nages || s<nseas)) {Ave_Size(t+1,0,g,k2) = temp1;}
                    }
                  }  // done ageloop
//  ProgLabel_16.2.3.2.1.3  after age loop, if(s=nseas) get weighted average for size_at_maxage from carryover fish and fish newly moving into this age
                  if(s==nseas)
                  {
                    temp=( (natage(t,1,g,nages-1)+0.01)*Ave_Size(t+1,0,g,nages) + (natage(t,1,g,nages)+0.01)*temp1) / (natage(t,1,g,nages-1)+natage(t,1,g,nages)+0.02);
                    Ave_Size(t+1,0,g,nages)=temp;
                  }
                  break;
                }

                case 2:   //  Richards
                {
                  growhalf = (mfexp(VBK(gp,nages)*seasdur_half(s)*VBK_seas(s))-1.0);
                  for(a=0;a<=nages;a++)
                  {
//  ProgLabel_16.2.3.2.2  von Bert growth with Richards coefficient, loop ages
//  ProgLabel_16.2.3.2.2.1  calc midseas size-at-age from begin season size-at-age, accounting for transition from linear to von Bert as necessary
                    if(lin_grow_beg(g,s,a)>=0.0)  // in linear phase for begin seas
                    {
                      Ave_Size(t,0,g,a) = len_bins(1)+lin_grow_beg(g,s,a)*(Cohort_Lmin(gp,y,a)-len_bins(1));  // begin of season
                    }
                    if(lin_grow_mid(g,s,a)>=0.0)  // in linear phase for midseas
                    {
                      Ave_Size(t,1,g,a) = len_bins(1)+lin_grow_mid(g,s,a)*(Cohort_Lmin(gp,y,a)-len_bins(1));  // mid season
                    }

                    if(lin_grow_beg(g,s,a)==-1.0)  // first time point beyond AFIX
                    {
                      temp=pow(Cohort_Lmin(gp,y,a),Richards(gp));
                      temp1=temp+(temp-LinfR)*(mfexp(VBK(gp,nages)*(curr_age_beg(g,s,a)-AFIX)*VBK_seas(s))-1.0)*Cohort_Growth(y,a);
                      Ave_Size(t,0,g,a) = pow(temp1,inv_Richards);
                    }

                    if(lin_grow_mid(g,s,a)<0.0)  //  so doing growth curve by mid of this seas, so project with growth to begin of next
                    {
                      if(lin_grow_mid(g,s,a)==-1.0)  // first time point beyond AFIX
                      {
                        temp=pow(Cohort_Lmin(gp,y,a),Richards(gp));
                        temp1=temp+(temp-LinfR)*(mfexp(VBK(gp,nages)*(curr_age_mid(g,s,a)-AFIX)*VBK_seas(s))-1.0)*Cohort_Growth(y,a);
                        Ave_Size(t,1,g,a) = pow(temp1,inv_Richards);
                      }
                      else  //  calc midseas size starting from size at beginning of this season
                      {
                        temp=pow(Ave_Size(t,0,g,a),Richards(gp));
                        t2=temp-LinfR;  //  remaining growth potential
                        if(time_vary_MG(y,2)>0 && t2>-1.)
                        {
                          join1=1.0/(1.0+mfexp(-(50.*t2/(1.0+fabs(t2)))));  //  note the logit transform is not perfect, so growth near Linf will not be exactly same as with native growth function
                          t2*=(1.-join1);  // trap to prevent decrease in size-at-age
                        }
                        Ave_Size(t,1,g,a) = pow((temp+growhalf*(t2)*Cohort_Growth(y,a)),inv_Richards);
                      }

//  ProgLabel_16.2.3.2.2.2  calc size at end of the season, which will be size at begin of next season using current seasons growth parms
                    if(a<nages) {k2=a+add_age;} else {k2=a;}  // where add_age =1 if s=nseas, else 0
                      temp=pow(Ave_Size(t,1,g,a),Richards(gp));
                      t2=temp-LinfR;  //  remaining growth potential
                      if(time_vary_MG(y,2)>0 && t2>-1.)
                      {
                        join1=1.0/(1.0+mfexp(-(50.*t2/(1.0+fabs(t2)))));  //  note the logit transform is not perfect, so growth near Linf will not be exactly same as with native growth function
                        t2*=(1.-join1);  // trap to prevent decrease in size-at-age
                      }
                      temp1 = pow((temp+growhalf*t2*Cohort_Growth(y,a)),inv_Richards);
                      if((a<nages || s<nseas)) {Ave_Size(t+1,0,g,k2) = temp1;}
                    }
                  }  // done age loop
//  ProgLabel_16.2.3.2.2.3  // after age loop, if(s=nseas) get weighted average for size_at_maxage from carryover fish and fish newly moving into this age
                  if(s==nseas)
                  {
                    temp=( (natage(t,1,g,nages-1)+0.01)*Ave_Size(t+1,0,g,nages) + (natage(t,1,g,nages)+0.01)*temp1) / (natage(t,1,g,nages-1)+natage(t,1,g,nages)+0.02);
                    Ave_Size(t+1,0,g,nages)=temp;
                  }
                  break;
                }
              }
                  if(docheckup==1) echoinput<<"Size_at_age for season: "<<s<<endl
                  <<" size-st  "<<Ave_Size(t,0,g)(0,min(6,nages))<<" "<<Ave_Size(t,0,g,nages)<<endl
                  <<" size-mid "<<Ave_Size(t,1,g)(0,min(6,nages))<<" "<<Ave_Size(t,1,g,nages)<<endl
                  <<" size+1   "<<Ave_Size(t+1,0,g)(0,min(6,nages))<<" "<<Ave_Size(t+1,0,g,nages)<<endl;

//  ProgLabel_16.2.3.2.3  do calculations related to std.dev. of size-at-age
//  ProgLabel_16.2.3.2.3.1 if (y=styr), calc CV_G(gp,s,a) by interpolation on age or LAA
//  doing this just at y=styr prevents the CV from changing as time-vary growth updates over time
              if(CV_const(gp)>0 && y==styr)
              {
                for (a=0;a<=nages;a++)
                {
                  if(curr_age_beg(g,s,a)<AFIX)
                  {CV_G(gp,s,a)=CVLmin(gp);}
                  else if(curr_age_beg(g,s,a)>=AFIX2_forCV)
                  {CV_G(gp,s,a)=CVLmax(gp);}
                  else if(CV_depvar_a==0)
                  {CV_G(gp,s,a)=CVLmin(gp) + (Ave_Size(t,0,g,a)-Lmin(gp))*CV_delta(gp);}
                  else
                  {CV_G(gp,s,a)=CVLmin(gp) + (curr_age_beg(g,s,a)-AFIX)*CV_delta(gp);}

                  if(curr_age_mid(g,s,a)<AFIX)
                  {CV_G(gp,smid,a)=CVLmin(gp);}
                  else if(curr_age_mid(g,s,a)>=AFIX2_forCV)
                  {CV_G(gp,smid,a)=CVLmax(gp);}
                  else if(CV_depvar_a==0)
                  {CV_G(gp,smid,a)=CVLmin(gp) + (Ave_Size(t,1,g,a)-Lmin(gp))*CV_delta(gp);}
                  else
                  {CV_G(gp,smid,a)=CVLmin(gp) + (curr_age_mid(g,s,a)-AFIX)*CV_delta(gp);}
                }   // end age loop
              }
              else
              {
                //  already set constant to CVLmi
              }

//  ProgLabel_16.2.3.2.3.2  calc stddev of size-at-age from CV_G(gp,s,a) and Ave_Size(t,g,a)
              if(CV_depvar_b==0)
              {
                Sd_Size_within(s,g)=SD_add_to_LAA+elem_prod(CV_G(gp,s),Ave_Size(t,0,g));
                Sd_Size_within(smid,g)=SD_add_to_LAA+elem_prod(CV_G(gp,smid),Ave_Size(t,1,g));
              }
              else
              {
                Sd_Size_within(s,g)=SD_add_to_LAA+CV_G(gp,s);
                Sd_Size_within(smid,g)=SD_add_to_LAA+CV_G(gp,smid);
              }
//  ProgLabel_16.2.3.2.3.3  if platoons being used, calc the stddev between platoons
              if(N_platoon>1)
              {
                Sd_Size_between(s,g)=Sd_Size_within(s,g)*sd_between_platoon;
                Sd_Size_between(smid,g)=Sd_Size_within(smid,g)*sd_between_platoon;
                Sd_Size_within(s,g)*=sd_within_platoon;
                Sd_Size_within(smid,g)*=sd_within_platoon;
              }
              if(docheckup==1)
              {
              echoinput<<y<<" "<<s<<" "<<g<<" "<<gp<<" CV_g_beg "<<CV_G(gp,s)<<endl;
              echoinput<<y<<" "<<s<<" "<<g<<" "<<gp<<" sd_beg "<<Sd_Size_within(s,g)<<endl;
              echoinput<<y<<" "<<s<<" "<<g<<" "<<gp<<" CV_g_mid "<<CV_G(gp,smid)<<endl;
              echoinput<<y<<" "<<s<<" "<<g<<" "<<gp<<" sd_mid "<<Sd_Size_within(smid,g)<<endl;
              }
            }

//  ProgLabel_16.2.3.3  propagate Ave_Size from early years forward until first year that has time-vary growth
            k=y+1;
            j=yz+1;
            while(time_vary_MG(j,2)==0 && k<=YrMax)
            {
              for (s=1;s<=nseas;s++)
              {
                smid=s+nseas;
                t=styr+(k-styr)*nseas+s-1;
                Ave_Size(t,0,g)=Ave_Size(t-nseas,0,g);
                if(s==1 && k<YrMax)
                {
                  Ave_Size(t+nseas,0,g)=Ave_Size(t,0,g);  // prep for time-vary next yr
                }
                Ave_Size(t,1,g)=Ave_Size(t-nseas,1,g);
              }  // end season loop
              k++;
              if(j<endyr+1) j++;

            }
          }  // end need to consider this GP x birthseason combo (usemorph>0)
        }  // end loop of birthseasons
        Ip+=MGparm_per_def;
      }    // end loop of growth patterns, gp
 // end do growth
  }

//*********************************************************************
//  ProgLabel_FUNCTION 17 get_natmort
FUNCTION void get_natmort()
  {
  dvariable Loren_M1;
  dvariable Loren_temp;
  dvariable Loren_temp2;
  dvariable t_age;
  int Do_AveAge;
  Do_AveAge=0;
  t_base=styr+(yz-styr)*nseas-1;
  Ip=-MGparm_per_def;   // start counter for MGparms
//  ProgLabel_17.1  loop growth patterns N_GP*gender
  for (gp=1;gp<=N_GP*gender;gp++)
  {
  	if(N_natMparms>0)
  	{
//  ProgLabel_17.1.1  copy parameter values from mgp_adj to natMparms(gp), doing direct or offset for gp>1
    Ip+=MGparm_per_def;
    for(j=1;j<=N_natMparms;j++) natMparms(j,gp)=mgp_adj(Ip+j);
    switch(MGparm_def)   //  switch for natmort parms
    {
      case 1:  // direct
      {
      	for(j=1;j<=N_natMparms;j++)
      	{
      		if(natMparms(j,gp)<0) natMparms(j,gp)=natMparms(j,1);
      	}
        break;
      }
      case 2:  // offset
      {
        if(gp>1)
        {
          for(j=1;j<=N_natMparms;j++)
          {
            natMparms(j,gp)=natMparms(j,1)*mfexp(natMparms(j,gp));
          }
        }
        break;
      }
      case 3:  // offset like SS2 V1.23
      {
          if(gp>1) natMparms(1,gp)=natMparms(1,1)*mfexp(natMparms(1,gp));
          if(N_natMparms>1)
          {
          for(j=2;j<=N_natMparms;j++)
          {
            natMparms(j,gp)=natMparms(j-1,gp)*mfexp(natMparms(j,gp));
          }
        }
        break;
      }
    }  // end switch
    }  // end have natmort parms

    g=(gp-1)*nseas*N_platoon+int(N_platoon/2)+1-N_platoon;
    for(birthseas=1;birthseas<=nseas;birthseas++)  // loop birth seasons (really looping more morphs)
    {
//  ProgLabel_17.1.2  loop birthseasons
      g+=N_platoon;
      gg=GP3(g);   // gpattern*gender*birthseason
      if(use_morph(g)>0)
      {
        switch(natM_type)
        {
//  ProgLabel_17.1.2.0  case 0:  constant M
          case 0:  // constant M
          {
            for(s=1;s<=nseas;s++)
            {
              natM(s,gg)=natMparms(1,gp);
              surv1(s,gg)=mfexp(-natMparms(1,gp)*seasdur_half(s));   // refers directly to the constant value
              surv2(s,gg)=square(surv1(s,gg));
            }
            break;
          }

//  ProgLabel_17.1.2.1  case 1:  N breakpoints
          case 1:  // breakpoints
          {
            dvariable natM1;
            dvariable natM2;
            for(s=1;s<=nseas;s++)
            {
              if(s>=Bseas(g))
              {a=0; t_age=azero_seas(s)-azero_G(g);}
              else
              {a=1; t_age=1.0+azero_seas(s)-azero_G(g);}
              natM_amax=NatM_break(1);
              natM2=natMparms(1,gp);
              k=a;

              for(loop=1;loop<=N_natMparms+1;loop++)
              {
                natM_amin=natM_amax;
                natM1=natM2;
                if(loop<=N_natMparms)
                {
                  natM_amax=NatM_break(loop);
                  natM2=natMparms(loop,gp);
                }
                else
                {
                  natM_amax=r_ages(nages)+1.;
                }
                if(natM_amax>natM_amin)
                {temp=(natM2-natM1)/(natM_amax-natM_amin);}  //  calc the slope
                else
                {temp=0.0;}
                while(t_age<natM_amax && a<=nages)
                {
                  natM(s,gg,a)=natM1+(t_age-natM_amin)*temp;
                  t_age+=1.0; a++;
                }
              }
              if(k==1) natM(s,gg,0)=natM(s,gg,1);
              surv1(s,gg)=mfexp(-natM(s,gg)*seasdur_half(s));
              surv2(s,gg)=square(surv1(s,gg));
            } // end season
            break;
          }  // end natM_type==1

//  ProgLabel_17.1.2.2  case 2:  lorenzen M
          case 2:  //  Lorenzen M
          {
            Loren_temp2=L_inf(gp)*(mfexp(-VBK(gp,nages)*VBK_seas(0))-1.);   // need to verify use of VBK_seas here
            t=styr+(yz-styr)*nseas+Bseas(g)-1;
            Loren_temp=Ave_Size(styr,1,g,int(natM_amin));  // uses mean size in middle of season 1 for the reference age
            Loren_M1=natMparms(1,gp)/log(Loren_temp/(Loren_temp+Loren_temp2));
            for(s=nseas;s>=1;s--)
            {
              smid=s+nseas;    //  mid-season size-at-age is stored in season s+nseas
              for (a=nages; a>=0;a--)
              {
                if(a==0 && s<Bseas(g))
                {natM(s,gg,a)=natM(s+1,gg,a);}
                else
                {natM(s,gg,a)=log(Ave_Size(t,1,g,a)/(Ave_Size(t,1,g,a)+Loren_temp2))*Loren_M1;}
                surv1(s,gg,a)=mfexp(-natM(s,gg,a)*seasdur_half(s));
                surv2(s,gg,a)=square(surv1(s,gg,a));
              }   // end age loop
            }
            break;
          }
//  ProgLabel_17.1.2.3  case 3:  set to empirical M as read from file, no seasonal interpolation
          case(3):   // read age_natmort as constant
          {
            for(s=1;s<=nseas;s++)
            {
              natM(s,gg)=Age_NatMort(gp);
              surv1(s,gg)=value(mfexp(-natM(s,gg)*seasdur_half(s)));
              surv2(s,gg)=value(square(surv1(s,gg)));
            }
            break;
          }

//  ProgLabel_17.1.2.4  case 4:  read age_natmort as constant and interpolate to seasonal real age
          case(4):
          {
            for(s=1;s<=nseas;s++)
            {
              if(s>=Bseas(g))
              {
                k=0; t_age=azero_seas(s)-azero_G(g);
                for (a=k;a<=nages-1;a++)
                {
                  natM(s,gg,a) = Age_NatMort(gp,a)+t_age*(Age_NatMort(gp,a+1)-Age_NatMort(gp,a));
                } // end age
              }
              else
              {
                k=1; t_age=azero_seas(s)+(1.-azero_G(g));
                for (a=k;a<=nages-1;a++)
                {
                  natM(s,gg,a) = Age_NatMort(gp,a)+t_age*(Age_NatMort(gp,a+1)-Age_NatMort(gp,a));
                } // end age
                natM(s,gg,0)=natM(s,gg,1);
              }
              natM(s,gg,nages)=Age_NatMort(gp,nages);
              surv1(s,gg)=mfexp(-natM(s,gg)*seasdur_half(s));
              surv2(s,gg)=square(surv1(s,gg));
            } // end season
            break;
          }
        }  // end natM_type switch

//  ProgLabel_17.2  calc an ave_age for the first gp as a scaling factor in logL for initial recruitment (R1) deviation
        if(Do_AveAge==0)
        {
          Do_AveAge=1;
          ave_age = 1.0/natM(1,gg,nages/2)-0.5;
        }
      } //  end use of this morph
    } // end birthseason
  }   // end growth pattern x gender loop
// end nat mort
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 18 get_recr_distribution */
FUNCTION void get_recr_distribution()
  {
  dvar_vector recr_dist_parm(1,MGP_CGD-recr_dist_parms-1);
  recr_dist_calc.initialize();
//  ProgLabel_18.1  set rec_dist_parms = exp(mgp_adj) for this year
  Ip=recr_dist_parms;
  for(f=1;f<=MGP_CGD-recr_dist_parms-1;f++)
  {
    recr_dist_parm(f)=mfexp(mgp_adj(Ip+f));
  }
//  ProgLabel_18.2  loop gp * birthseas * area and multiply together the recr_dist_parm values
  for(gp=1;gp<=N_GP;gp++)
  for(s=1;s<=nseas;s++)         //  really looping birthseas
  for(p=1;p<=pop;p++)
  if(recr_dist_pattern(gp,s,p)>0)
  {
    recr_dist_calc(gp,s,p)=femfrac(gp)*recr_dist_parm(gp)*recr_dist_parm(N_GP+p)*recr_dist_parm(N_GP+pop+s);
    if(gender==2) recr_dist_calc(gp+N_GP,s,p)=femfrac(gp+N_GP)*recr_dist_parm(gp)*recr_dist_parm(N_GP+p)*recr_dist_parm(N_GP+pop+s);  //males
  }

//  ProgLabel_18.3  if recr_dist_interaction is chosen, then multiply these in also
  if(recr_dist_inx==1)
  {
    f=N_GP+nseas+pop;
    for(gp=1;gp<=N_GP;gp++)
    for(s=1;s<=nseas;s++)
    for(p=1;p<=pop;p++)
    {
      f++;
      if(recr_dist_pattern(gp,s,p)>0)
      {
        recr_dist_calc(gp,s,p)*=recr_dist_parm(f);
        if(gender==2) recr_dist_calc(gp+N_GP,s,p)*=recr_dist_parm(f);
      }
    }
  }
//  ProgLabel_18.4  scale the recr_dist matrix to sum to 1.0
  recr_dist_calc/=sum(recr_dist_calc);
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 19 get_wtlenet_wtlen, maturity, fecundity, hermaphroditism */
FUNCTION void get_wtlen()
  {

//  ProgLabel_19.1  set wtlen and maturity/fecundity factors equal to annual values from mgp_adj
  Ip=MGparm_per_def*N_GP*gender+1;
  for(f=0;f<=3+2*gender;f++) {wtlen_p(f+1)=mgp_adj(Ip+f);}

  for(s=1;s<=nseas;s++)
  {
//  ProgLabel_19.2   loop seasons
    t=styr+(y-styr)*nseas+s-1;
//  ProgLabel_19.2.1  calc wt_at_length for each season to include seasonal effects on wtlen
    if(MGparm_seas_effects(1)>0 || MGparm_seas_effects(2)>0 )        //  get seasonal effect on fem wtlen parameters
    {
      wt_len(s,1)=(wtlen_p(1)*wtlen_seas(s,1))*pow(len_bins_m(1,nlength),(wtlen_p(2)*wtlen_seas(s,2)));
      wt_len_low(s)(1,nlength)=(wtlen_p(1)*wtlen_seas(s,1))*pow(len_bins2(1,nlength),(wtlen_p(2)*wtlen_seas(s,2)));
    }
    else
    {
      wt_len(s,1) = wtlen_p(1)*pow(len_bins_m(1,nlength),wtlen_p(2));
      wt_len_low(s)(1,nlength) = wtlen_p(1)*pow(len_bins2(1,nlength),wtlen_p(2));
    }

    wt_len2(s)(1,nlength)=wt_len(s,1)(1,nlength);

//  ProgLabel_19.2.2  calculate maturity and fecundity if seas = spawn_seas
//  these calculations are done in spawn_seas, but are not affected by spawn_time within that season
//  so age-specific inputs will assume to be at correct timing already; size-specific will later be adjusted to use size-at-age at the exact correct spawn_time_seas

    if(s==spawn_seas)  // get biology of maturity and fecundity
    {
      for(gp=1;gp<=N_GP;gp++)
      {

        switch(Maturity_Option)
        {
          case 1:  //  Maturity_Option=1  length logistic
          {
            mat_len(gp) = 1./(1. + mfexp(wtlen_p(4)*(len_bins_m(1,nlength)-wtlen_p(3))));
            break;
          }
          case 2:  //  Maturity_Option=2  age logistic
          {
            mat_age(gp) = 1./(1. + mfexp(wtlen_p(4)*(r_ages-wtlen_p(3))));
            break;
          }
          case 3:  //  Maturity_Option=3  read age-maturity
          {
            mat_age(gp)=Age_Maturity(gp);
            break;
          }
          case 4:  //  Maturity_Option=4   read age-fecundity, so no age-maturity
          {
            break;
          }
          case 5:  //  Maturity_Option=5   read age-fecundity from wtatage.ss
          {
            break;
          }
          case 6:  //  Maturity_Option=6   read length-maturity
          {
            mat_len(gp)=Length_Maturity(gp);
            break;
          }
        }

        if(First_Mature_Age>0)
        {mat_age(gp)(0,First_Mature_Age-1)=0.;}
          
        switch (Fecund_Option)
        {
          case 1:    // as eggs/kg (SS original configuration)
          {
            fec_len(gp) = wtlen_p(5)+wtlen_p(6)*wt_len(s,1);
            fec_len(gp) = elem_prod(wt_len(s,1),fec_len(gp));
            break;
          }
          case 2:
          {       // as eggs = f(length)
            fec_len(gp) = wtlen_p(5)*pow(len_bins_m,wtlen_p(6));
            break;
          }
          case 3:
          {       // as eggs = f(body weight)
            fec_len(gp) = wtlen_p(5)*pow(wt_len(s,1),wtlen_p(6));
            break;
          }
          case 4:
          {       // as eggs = a + b*Len
            fec_len(gp) = wtlen_p(5) + wtlen_p(6)*len_bins_m;
            if(wtlen_p(5)<0.0)
            {
              z=1;
              while(fec_len(gp,z)<0.0)
              {
                fec_len(gp,z)=0.0;
                z++;
              }
            }
            break;
          }
          case 5:
          {       // as eggs = a + b*Wt
            fec_len(gp) = wtlen_p(5) + wtlen_p(6)*wt_len(s,1);
            if(wtlen_p(5)<0.0)
            {
              z=1;
              while(fec_len(gp,z)<0.0)
              {
                fec_len(gp,z)=0.0;
                z++;
              }
            }
            break;
          }
        }

//  combine length maturity and fecundity; but will be ignored if reading empirical age-fecundity
        mat_fec_len(gp) = elem_prod(mat_len(gp),fec_len(gp));

      }
    }

//  ProgLabel_19.2.3  calculate male weight_at_length
    if(gender==2)
    {
      if(MGparm_seas_effects(7)>0 || MGparm_seas_effects(8)>0 )        //  get seasonal effect on male wt-len parameters
      {
        wt_len(s,2) = (wtlen_p(7)*wtlen_seas(s,7))*pow(len_bins_m(1,nlength),(wtlen_p(8)*wtlen_seas(s,8)));
        wt_len_low(s)(nlength1,nlength2) = (wtlen_p(7)*wtlen_seas(s,7))*pow(len_bins2(nlength1,nlength2),(wtlen_p(8)*wtlen_seas(s,8)));
      }
      else
      {
        wt_len(s,2) = wtlen_p(7)*pow(len_bins_m(1,nlength),wtlen_p(8));
        wt_len_low(s)(nlength1,nlength2) = wtlen_p(7)*pow(len_bins2(nlength1,nlength2),wtlen_p(8));
      }
      wt_len2(s)(nlength1,nlength2)=wt_len(s,2).shift(nlength1);
      wt_len(s,2).shift(1);
    }

//  ProgLabel_19.2.4  calculate first diff of wt_len for use in generalized sizp comp bin calculations
    wt_len_fd(s)=first_difference(wt_len_low(s));
    wt_len2_sq(s)=elem_prod(wt_len2(s),wt_len2(s));
    if(gender==2) wt_len_fd(s,nlength)=wt_len_fd(s,nlength-1);

//  ProgLabel_19.2.5  Do Hermaphroditism (no seasonality and no gp differences)
//  should build seasonally component here
    if(Hermaphro_Option>0)
    {
        dvariable infl;  // inflection
        dvariable stdev;  // standard deviation
        dvariable maxval;  // max value

      for(gp=1;gp<=N_GP; gp++)
      {
         infl=mgp_adj(MGparm_Hermaphro);  // inflection
         stdev=mgp_adj(MGparm_Hermaphro+1);  // standard deviation
         maxval=mgp_adj(MGparm_Hermaphro+2);  // max value
//      minval is 0.0;
        temp2=cumd_norm((0.0-infl)/stdev);     //  cum_norm at age 0
        temp=maxval / (cumd_norm((r_ages(nages)-infl)/stdev)-temp2);   //  delta in cum_norm between styr and endyr
        for(a=0; a<=nages; a++)
        {
          Hermaphro_val(gp,a)=0.0 + temp * (cumd_norm((r_ages(a)-infl)/stdev)-temp2);
        }
      }
    }
  }
//  end wt-len and fecundity
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 20 get_migration */
FUNCTION void get_migration()
  {
  Ip=MGP_CGD;   // start counter for  movement parms
//  ProgLabel_20.1  loop the needed movement rates
  for(k=1;k<=do_migr2;k++)   //  loop all movement rates for this year (includes seas, morphs)
  {
    t=styr+(yz-styr)*nseas+move_def2(k,1)-1;
    if(k<=do_migration) 
    {
//  ProgLabel_20.1.1  age-specific movement strength based on parameters for selected area pairs
      temp=1./(move_def2(k,6)-move_def2(k,5));
      temp1=temp*(mgp_adj(Ip+2)-mgp_adj(Ip+1));
      for (a=0;a<=nages;a++)
      {
        if(a<=move_def2(k,5)) {migrrate(yz,k,a) = mgp_adj(Ip+1);}
        else if(a>=move_def2(k,6)) {migrrate(yz,k,a) = mgp_adj(Ip+2);}
        else {migrrate(yz,k,a) = mgp_adj(Ip+1) + (r_ages(a)-move_def2(k,5))*temp1;}
      }   // end age loop
      migrrate(yz,k)=mfexp(migrrate(yz,k));
      Ip+=2;
    }
    else
//  ProgLabel_20.1.2  default movement strength =1.0 for other area pairs
    {
      migrrate(yz,k)=1.;
    }
  }

//  ProgLabel_20.2  loop seasons, GP, source areas
  for (s=1;s<=nseas;s++)
  {
    t=styr+(yz-styr)*nseas+s-1;
    for(gp=1;gp<=N_GP;gp++)
    {
      for (p=1;p<=pop;p++)
      {
        tempvec_a.initialize();   // zero out the summation vector
        for (p2=1;p2<=pop;p2++)
        {
//  ProgLabel_20.2.1  for each destination area, adjust movement rate by season duration and sum across all destination areas
          k=move_pattern(s,gp,p,p2);
          if(k>0)
          {
            if(p2!=p && nseas>1) migrrate(yz,k)*=seasdur(move_def2(k,1));  // fraction leaving an area is reduced if the season is short
            tempvec_a+=migrrate(yz,k);          //  sum of all movement weights for the p2 fish
          }
        }   //end destination area
//  ProgLabel_20.2.2  now normalize for all movement from source area p
        for (p2=1;p2<=pop;p2++)
        {
          k=move_pattern(s,gp,p,p2);
          if(k>0)
          {
            migrrate(yz,k)=elem_div(migrrate(yz,k),tempvec_a);
//  ProgLabel_20.2.3  now set rate to 0.0 (or 1.0 for stay rates) below the start age for migration
            if(migr_start(s,gp)>0)
            {
              if(p!=p2)
              {
                migrrate(yz,k)(0,migr_start(s,gp)-1)=0.0;
              }
              else
              {
                migrrate(yz,k)(0,migr_start(s,gp)-1)=1.0;
              }
            }
          }
        }
      }    //  end source areas loop
    }  // end growth pattern
  }  // end season

//  ProgLabel_20.2.4  copy annual migration rates forward until first year with time-varying migration rates
  if(yz<endyr)
  {
    k=yz+1;
    while(time_vary_MG(k,5)==0 && k<=endyr)
    {
      migrrate(k)=migrrate(k-1);  k++;
    }
  }
//  end migration
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 21 get_saveGparm */
FUNCTION void get_saveGparm()
  {
      for (gp=1;gp<=N_GP*gender;gp++)
      {
        g=(gp-1)*nseas*N_platoon+int(N_platoon/2)+1-N_platoon;
        for(birthseas=1;birthseas<=nseas;birthseas++)  // loop birth seasons (really looping more morphs)
        {
          g+=N_platoon;
          save_gparm++;
          save_G_parm(save_gparm,1)=save_gparm;
          save_G_parm(save_gparm,2)=y;
          save_G_parm(save_gparm,3)=g;
          save_G_parm(save_gparm,4)=AFIX;
          save_G_parm(save_gparm,5)=AFIX2;
          save_G_parm(save_gparm,6)=value(Lmin(gp));
          save_G_parm(save_gparm,7)=value(Lmax_temp(gp));
          save_G_parm(save_gparm,8)=value(-VBK(gp,nages)*VBK_seas(0));
          save_G_parm(save_gparm,9)=value( -log(L_inf(gp)/(L_inf(gp)-Lmin(gp))) / (-VBK(gp,nages)*VBK_seas(0)) +AFIX+azero_G(g) );
          save_G_parm(save_gparm,10)=value(L_inf(gp));
          save_G_parm(save_gparm,11)=value(CVLmin(gp));
          save_G_parm(save_gparm,12)=value(CVLmax(gp));
          save_G_parm(save_gparm,13)=natM_amin;
          save_G_parm(save_gparm,14)=natM_amax;
          save_G_parm(save_gparm,15)=value(natM(1,GP3(g),0));
          save_G_parm(save_gparm,16)=value(natM(1,GP3(g),nages));
          if(sx(g)==1)
          {
          for(k=1;k<=6;k++) save_G_parm(save_gparm,16+k)=value(wtlen_p(k));
          }
          else
          {
          for(k=1;k<=2;k++) save_G_parm(save_gparm,16+k)=value(wtlen_p(k+6));
          }
          save_gparm_print=save_gparm;
        }
      }
    if(MGparm_doseas>0)
      {
        for(s=1;s<=nseas;s++)
        {
          for(k=1;k<=8;k++)
          {
          save_seas_parm(s,k)=value(wtlen_p(k)*wtlen_seas(s,k));
          }
          save_seas_parm(s,9)=value(Lmin(1));
          save_seas_parm(s,10)=value(VBK(1,nages)*VBK_seas(s));
        }
      }
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 22 get_selectivity */
FUNCTION void get_selectivity()
  {
  int Ip_env;
  int y1;
  int fs;
  dvariable t1;
  dvariable t2;
  dvariable t3;
  dvariable t4;
  dvariable Apical_Selex;
  dvariable t1min; dvariable t1max; dvariable t1power;
  dvariable t2min; dvariable t2max; dvariable t2power; dvariable final; dvariable sel_maxL;
  dvariable lastsel; dvariable lastSelPoint; dvariable SelPoint; dvariable finalSelPoint;
  dvariable asc;
  dvariable dsc;

  dvar_vector sp(1,199);                 // temporary vector for selex parms

  // begin new section for testing cubic spline
  // define vectors which form the basis for cubic spline selectivity
  // IMPORTANT: these vectors might need to be expanded to fit values for multiple fleets
  dvector splineX(1,200);
  dvar_vector splineY(1,200);
  splineX.initialize();
  splineY.initialize();
  // end section for testing cubic spline

  Ip=0;
//  dvariable temp3;

//  calculate any trends that will be needed for any of the selex parameters
  if(y==styr && N_selparm_trend>0)
  {
    for(f=1;f<=N_selparm_trend;f++)
    {
      j=selparm_trend_rev(f);  //  parameter affected; would include size selex, discard, discmort, ageselex
      k=selparm_trend_rev_1(f);  // base index for trend parameters
      //  calc endyr value, but use logistic transform to keep with bounds of the base parameter
      if(selparm_1(j,13)==-1)
      {
        temp=log((selparm_1(j,2)-selparm_1(j,1)+0.0000002)/(selparm(j)-selparm_1(j,1)+0.0000001)-1.)/(-2.);   // transform the base parameter
        temp+=selparm(k+1);     //  add the offset  Note that offset value is in the transform space
        temp1=selparm_1(j,1)+(selparm_1(j,2)-selparm_1(j,1))/(1.+mfexp(-2.*temp));   // backtransform
      }
      else if(selparm_1(j,13)==-2)
      {
        temp1=selparm(k+1);  // set ending value directly
      }

      if(selparm_HI(k+2)<=1.1)  // use max bound as switch
      {temp3=r_years(styr)+selparm(k+2)*(r_years(endyr)-r_years(styr));}  // infl year
      else
      {temp3=selparm(k+2);}

      temp2=cumd_norm((r_years(styr)-temp3)/selparm(k+3));     //  cum_norm at styr
      temp=(temp1-selparm(j)) / (cumd_norm((r_years(endyr)-temp3)/selparm(k+3))-temp2);   //  delta in cum_norm between styr and endyr
      for(int y1=styr;y1<=endyr;y1++)
      {
        selparm_trend(f,y1)=selparm(j) + temp * (cumd_norm((r_years(y1)-temp3)/selparm(k+3) )-temp2);
      }
    }
  }

  if (N_selparm_blk>0 && y==styr)  // set up the block values time series
  {
    for(j=1;j<=N_selparm;j++)
    {
      z=selparm_1(j,13);    // specified block pattern
      if(z>0)  // uses blocks
      {
//        if(y==styr)
        {
          g=1;
          if(selparm_1(j,14)<3)  //  not as offset from previous block
          {
            for(a=1;a<=Nblk(z);a++)
            {
              for(int y1=Block_Design(z,g);y1<=Block_Design(z,g+1);y1++)  // loop years for this block
              {
                if(y1<=endyr)
                {
                  k=Block_Defs_Sel(j,y1);  // identifies parameter that holds the block effect
                  selparm_block_val(j,y1)=selparm(k);
                }
              }
              g+=2;
            }
          }
          else   // as additive offset to previous block
          {
            temp=0.0;
            for(a=1;a<=Nblk(z);a++)
            {
              y1=Block_Design(z,g);   // first year of block
              k=Block_Defs_Sel(j,y1);  // identifies parameter that holds the block effect
              temp+=selparm(k);  // increment by the block delta
              for(int y1=Block_Design(z,g);y1<=Block_Design(z,g+1);y1++)  // loop years for this block
              {
                if(y1<=endyr) selparm_block_val(j,y1)=temp;
              }
              g+=2;
            }
          }
        }  // end block setup
      }  // end uses blocks
    }  // end parameter loop
  }  // end block section

  if(N_selparm_dev>0 && y==styr)
  {
    for(k=1;k<=N_selparm_dev;k++)
    {
      if(selparm_dev_type(k)==3)
      {
        selparm_dev_rwalk(k,selparm_dev_minyr(k))=selparm_dev(k,selparm_dev_minyr(k));
        j=selparm_dev_minyr(k);
        for (j=selparm_dev_minyr(k)+1;j<=selparm_dev_maxyr(k);j++)
        {
          selparm_dev_rwalk(k,j)=selparm_dev_rwalk(k,j-1)+selparm_dev(k,j);
        }
      }
    }
  }

// **   Do Selectivity
  for (f=1;f<=2*Ntypes;f++)  // loops twice, first for size then for age
  {

    fs=f-Ntypes;  //index for saving age info in the fleet arrays

    if(time_vary_sel(y,f)==1 || save_for_report>0)
    {    // recalculate the selex in this year x type

      if(N_selparmvec(f)>0)      // type has parms, so look for adjustments
      {
        for(j=1;j<=N_selparmvec(f);j++)
        {
          if(selparm_1(Ip+j,13)<0)
          {sp(j)=selparm_trend(selparm_trend_point(Ip+j),y);}
          else
          {sp(j)=selparm(Ip+j);}
        }

        switch(selparm_adjust_method)
        {
          default:
          {
            break;
          }
          case 0:
          {
            break;
          }
          case 3:
          {
            // no break, so will do the case 1 code
          }
          case(1):
          {
            for(j=1;j<=N_selparmvec(f);j++)
            {
              if(selparm_1(Ip+j,13)>0)     //   uses blocks
              {
                blkparm=Block_Defs_Sel(Ip+j,y);
                if(blkparm>0)
                {
                  k=selparm_1(Ip+j,14);
                  if(k==0)  // multiplicative
                  {sp(j) *= mfexp(selparm_block_val(Ip+j,y));}
                  else if(k==1)  // additive
                  {sp(j) += selparm_block_val(Ip+j,y);}
                  else if(k==2)  // replacement
                  {sp(j) = selparm_block_val(Ip+j,y);}
                  else if(k==3)  // additive, but additions are cumulative
                  {sp(j) += selparm_block_val(Ip+j,y);}
                }
              }

              if(selparm_env(Ip+j)>0)       // if env then modify sp
              {
                if(selparm_envtype(Ip+j)==1)
                  {sp(j) *= mfexp(selparm(selparm_env(Ip+j))* env_data(y,selparm_envuse(Ip+j)));}
                else
                  {sp(j) += selparm(selparm_env(Ip+j))* env_data(y,selparm_envuse(Ip+j));}
              }

              k=selparm_dev_select(Ip+j);     // if dev then modify sp
              if(k>0)
              {
                if(y>=selparm_dev_minyr(k) && y<=selparm_dev_maxyr(k))
                {
                  if(selparm_dev_type(k)==1)
                  {sp(j) *= mfexp(selparm_dev(k,y));}
                  else if(selparm_dev_type(k)==2)
                  {sp(j) += selparm_dev(k,y);}
                  else if(selparm_dev_type(k)==3)
                  {sp(j)+=selparm_dev_rwalk(k,y);}
                }
              }
              if(selparm_adjust_method==1 && (save_for_report>0 || do_once==1))  // so does not check bounds if adjust_method==3
              {
                if(sp(j)<selparm_1(Ip+j,1) || sp(j)>selparm_1(Ip+j,2))
                {
                  N_warn++;
                  warning<<" adjusted selparm out of bounds (Parm#, yr, min, max, base, value) "<<
                  Ip+j<<" "<<y<<" "<<selparm_1(Ip+j,1)<<" "<<selparm_1(Ip+j,2)<<" "<<selparm(Ip+j)<<" "<<sp(j)<<endl;
                }
              }
            }  // end j parm loop
            break;
          }
          case(2):
          {
            for(j=1;j<=N_selparmvec(f);j++)
              {

              temp=log((selparm_1(Ip+j,2)-selparm_1(Ip+j,1)+0.0000002)/(sp(j)-selparm_1(Ip+j,1)+0.0000001)-1.)/(-2.);   // transform the parameter
              doit=0;
              if(selparm_1(Ip+j,13)>0)   // blocks
              {
                blkparm=Block_Defs_Sel(Ip+j,y);  // identifies parameter that holds the block effect
                if(blkparm>0)
                {
                  k=selparm_1(Ip+j,14);
                  doit=1;
                  if(k==1)
                    {temp += selparm_block_val(Ip+j,y);}
                  else if (k==2)  //  replacement
                    {temp=log((selparm_1(Ip+j,2)-selparm_1(Ip+j,1)+0.0000002)/(selparm_block_val(Ip+j,y)-selparm_1(Ip+j,1)+0.0000001)-1.)/(-2.);}  // block as replacement
                  else if(k==3)  // additive, but based on cumulative blocks
                    {temp += selparm_block_val(Ip+j,y);}
                  else
                    {N_warn++; warning<<" disabled multiplicative block effect with logistic approach"<<endl;}
                }
              }

              if(selparm_env(Ip+j)>0)  //  do environmental effect  only additive allowed
                {doit=1;temp += selparm(selparm_env(Ip+j))* env_data(y,selparm_envuse(Ip+j));}

              k=selparm_dev_select(Ip+j); //  Annual deviations;  use kth dev series
              if(k>0)
                {
                if(y>=selparm_dev_minyr(k) && y<=selparm_dev_maxyr(k))
                  {
                    doit=1;
                    if(selparm_dev_type(k)==2)
                    {temp += selparm_dev(k,y);}
                    else if(selparm_dev_type(k)==3)
                    {temp += selparm_dev_rwalk(k,y);}
                  }
                }
              if(doit==1) sp(j)=selparm_1(Ip+j,1)+(selparm_1(Ip+j,2)-selparm_1(Ip+j,1))/(1+mfexp(-2.*temp));   // backtransform
              }  // end parameter loop j
            break;
          }
        }
        if(save_for_report>0 || do_once==1)
        {for(j=1;j<=N_selparmvec(f);j++) save_sp_len(y,f,j)=sp(j);}
      }  // end adjustment of parms

      if(f<=Ntypes)  // do size selectivity, retention, discard mort
      {
      if(docheckup==1 && N_selparmvec(f)>0) echoinput<<"length selex parms: "<<sp(1,N_selparmvec(f))<<endl;
      for(gg=1;gg<=gender;gg++)
      {
        if(gg==1 || (gg==2 && seltype(f,3)>=3))
        {
          switch(seltype(f,1))  // select the selectivity pattern
          {
            case 0:   // ***********   constant
             {sel = 1.;break;}

            case 1:                      // ***********   logistic    ProgLabel_055
              {
                if(seltype(f,3)<3 || (gg==1 && seltype(f,3)==3) || (gg==2 && seltype(f,3)==4))  //  do the primary gender
                {sel = 1./(1.+mfexp(neglog19*(len_bins_m-sp(1))/sp(2)));}
                else  //  do the offset gender
                {
                  temp=sp(1)+sp(Maleselparm(f));
                  temp1=sp(2)+sp(Maleselparm(f)+1);
                  sel = sp(Maleselparm(f)+2)/(1.+mfexp(neglog19*(len_bins_m-temp)/temp1));
                }
                break;
              }

            case 2:                  // *******double logistic     ProgLabel_056
              {
                                     // 1=peak, 2=init,  3=infl,  4=slope, 5=final, 6=infl2, 7=slope2
            N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" selex pattern 2 discontinued; use pattern 8 for double logistic "<<endl; exit(1);
           break;
          }    // end double logistic

          case 3:                  // ********   ProgLabel_057
          {
            N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" selex pattern 3 discontinued "<<endl; exit(1);
           break;
          }  // end seltype=3

    //*********
          case 4:               //  use maturity curve as selex  ProgLabel_058
            {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" selex pattern 4 discontinued; use pattern 30 to get spawning biomass "<<endl; exit(1); break;}                   // do this as a numbers survey because wt is included here

    //**********
          case 5:      //  borrow another selex  ProgLabel_059
                                            //  use only the specified bin range
                                           // must refer to a lower numbered type (f)
          {
           i=int(value(sp(1)));  if(i<=0) i=1;
           j=int(value(sp(2)));  if(j<=0) j=nlength;
           if(j>nlength)
           {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" Critical error, size selex mirror length is too large for fleet: "<<f<<endl; exit(1);}
           sel.initialize();
           sel(i,j)=sel_l(y,seltype(f,4),1)(i,j);
           break;
          }

    //**********
          case 15:      //  borrow another selex  ProgLabel_059
          {
           sel.initialize();
           sel=sel_l(y,seltype(f,4),1);
           break;
          }
    //**********
          case 6:       //  non-parametric selex pattern
          {
          lastsel=-10.0;  // log(selex) for first bin;
          lastSelPoint=len_bins_m(1);    //  first size
          finalSelPoint=value(sp(2));  // size beyond which selex is constant
          SelPoint=value(sp(1));   //  first size that will get a parameter.  Value will get incremented by step interval (temp1)
          z=3;  // parameter counter
          temp1 = (finalSelPoint-SelPoint)/(seltype(f,4)-1.0);  // step interval
          for(j=1;j<=nlength;j++)
          {
            if(len_bins_m(j)<SelPoint)
            {
              tempvec_l(j)=lastsel + (len_bins_m(j)-lastSelPoint)/(SelPoint-lastSelPoint) * (sp(z)-lastsel);
            }
            else if(len_bins_m(j)==SelPoint)
            {
              tempvec_l(j)=sp(z);
              lastsel=sp(z);
              lastSelPoint=SelPoint;
              SelPoint+=temp1;
              if(SelPoint<=finalSelPoint)
                {z++;}
              else
                {SelPoint=finalSelPoint;}
            }
            else if(len_bins_m(j)<=finalSelPoint)
            {
              lastsel=sp(z);
              lastSelPoint=SelPoint;
              SelPoint+=temp1;
              if(SelPoint<=finalSelPoint)
                {z++;}
              else
                {SelPoint=finalSelPoint;}
              tempvec_l(j)=lastsel + (len_bins_m(j)-lastSelPoint)/(SelPoint-lastSelPoint) * (sp(z)-lastsel);
            }
            else
            {tempvec_l(j)=sp(z);}

          }
          temp=max(tempvec_l);
          sel = mfexp(tempvec_l-temp);
          break;
          }

          case 7:                  // *******New double logistic
    // 1=peak, 2=init,  3=infl,  4=slope, 5=final, 6=infl2, 7=slope2 8=binwidth;    Mirror=1===const_above_Linf
          {
            N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" selex pattern 7 discontinued; use pattern 8 for double logistic "<<endl; exit(1);
 /*
           t1=minL+(1./(1.+mfexp(-sp(3))))*(sp(1)-minL);    // INFL
           t1min=1./(1.+mfexp(-sp(4)*(minL-t1)))*0.9999;  // asc value at minsize
           t1max=1./(1.+mfexp(-sp(4)*(sp(1)-t1)))*1.00001;  // asc value at peak
           t1power=log(0.5)/log((0.5-t1min)/(t1max-t1min));  // so the parameter will actual correspond to 50% point

           if(seltype(f,4)==0) {sel_maxL=maxL;} else {sel_maxL=Ave_Size(styr,1,1,nages);}
           t2=(sp(1)+sp(8))+(1./(1.+mfexp(-sp(6))))*(sel_maxL-(sp(1)+sp(8)));    // INFL
           t2min=1./(1.+mfexp(-sp(7)*(sp(1)+sp(8)-t2)))*0.9999;  // asc value at peak+
           t2max=1./(1.+mfexp(-sp(7)*(sel_maxL-t2)))*1.00001;  // asc value at maxL
           t2power=log(0.5)/log((0.5-t2min)/(t2max-t2min));
           final=1./(1.+mfexp(-sp(5)));

           for (j=1; j<=nlength; j++)  //calculate the value over length bins
           {sel(j) =
             (
             (
             (sp(2) + (1. - sp(2)) * pow((( 1./(1.+mfexp(-sp(4)*(len_bins_m(j)-t1))) -t1min ) / (t1max-t1min) ),t1power))
              /(1.+mfexp(10.*(len_bins_m(j)-sp(1))))   // scale ascending side
              +
              1./(1.+mfexp(-10.*(len_bins_m(j)-sp(1))))   // flattop, with scaling
              )
              /(1.+mfexp( 10.*(len_bins_m(j)-(sp(1)+sp(8)))))    // scale combo of ascending and flattop
              +
              (1. + (final - 1.) * pow(sqrt(square((( 1./(1.+mfexp(-sp(7)*(len_bins_m(j)-t2))) -t2min ) / (t2max-t2min) ))),t2power))
              /(1.+mfexp( -10.*(len_bins_m(j)-(sp(1)+sp(8)))))    // scale descending
              ) / (1.+mfexp(10.*(len_bins_m(j)-sel_maxL)));       // scale combo of ascend, flattop, descending
             sel(j)+=final/(1.+mfexp(-10.*(len_bins_m(j)-sel_maxL)));  // add scaled portion above Linf
           }   // end size bin loop
  */
           break;
          }    // end New double logistic

    //**********
          case 8:                  // *******New double logistic in simpler code
    // 1=peak, 2=init,  3=infl,  4=slope, 5=final, 6=infl2, 7=slope2 8=binwidth;    Mirror=1===const_above_Linf
          {
           t1=minL+(1./(1.+mfexp(-sp(3))))*(sp(1)-minL);    // INFL
           t1min=1./(1.+mfexp(-mfexp(sp(4))*(minL-t1)))*0.9999;  // asc value at minsize
           t1max=1./(1.+mfexp(-mfexp(sp(4))*(sp(1)-t1)))*1.0001;  // asc value at peak
           t1power=log(0.5)/log((0.5-t1min)/(t1max-t1min));  // so the parameter will actual correspond to 50% point

           if(seltype(f,4)==0) {sel_maxL=maxL;} else {sel_maxL=Ave_Size(styr,1,1,nages);}
           t2=(sp(1)+sp(8))+(1./(1.+mfexp(-sp(6))))*(sel_maxL-(sp(1)+sp(8)));    // INFL
           t2min=1./(1.+mfexp(-mfexp(sp(7))*(sp(1)+sp(8)-t2)))*0.9999;  // asc value at peak+
           t2max=1./(1.+mfexp(-mfexp(sp(7))*(sel_maxL-t2)))*1.0001;  // asc value at maxL
           t2power=log(0.5)/log((0.5-t2min)/(t2max-t2min));
           final=1./(1.+mfexp(-sp(5)));
           for (j=1; j<=nlength; j++)  //calculate the value over length bins
           {join1=1./(1.+mfexp(10.*(len_bins_m(j)-sp(1))));
            join2=1./(1.+mfexp(10.*(len_bins_m(j)-(sp(1)+sp(8)))));
            join3=1./(1.+mfexp(10.*(len_bins_m(j)-sel_maxL)));
            upselex=sp(2) + (1. - sp(2)) * pow((( 1./(1.+mfexp(-mfexp(sp(4))*(len_bins_m(j)-t1)))-t1min ) / (t1max-t1min)),t1power);
            downselex=(1. + (final - 1.) * pow(fabs(((( 1./(1.+mfexp(-mfexp(sp(7))*(len_bins_m(j)-t2))) -t2min ) / (t2max-t2min) ))),t2power));
            sel(j) = ((((upselex*join1)+1.0*(1.0-join1))*join2) + downselex*(1-join2))*join3 + final*(1-join3);
           }   // end size bin loop
           break;
          }    // end New double logistic

    //***********************
          case 9:                      // ***********   old double logistic
          {k1=int(value(sp(5)));
           if(k1>1) sel(1,k1-1) = 0.0;
           sel(k1,nlength) =   elem_prod(  (1/(1+mfexp(-sp(2)*(len_bins_m(k1,nlength)-sp(1)) ))),
                                                (1-1/(1+mfexp(-sp(4)*(len_bins_m(k1,nlength)-(sp(1)*sp(6)+sp(3))) ))) );
           sel += 1.0e-6;
           sel /= max(sel);
           break;
            }

    //***********************************
          case 21:                 // *******New non-parametric
          {
            j=1;
            z=1;
            k=seltype(f,4);  // n points
            lastsel=0.0;
            lastSelPoint=0.0;

            if(do_once==1)
            {
              if(sp(k)>len_bins(nlength))
              {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<"Selex21: cannot have max selpoint > max_pop_lenbin"<<endl;  exit(1);}
              if(sp(k-1)>len_bins(nlength-1))
              {N_warn++; warning<<"Selex21: should not have selpoint(n-1) > pop_lenbin(nlength-1)"<<endl;}
            }

            while(j<=nlength)
            {
              if(len_bins(j)<=sp(z))
              {
                sel(j) = lastsel + (len_bins(j)-lastSelPoint)/(sp(z)-lastSelPoint) * (sp(z+k)-lastsel);
                j++;
              }
              else if(z<=k)
              {
                lastSelPoint=sp(z);
                lastsel=sp(z+k);
                z++;
              }
              else  //  for sizes beyond last point
              {
                sel(j)=sp(k+k);
                j++;
              }
            }
           break;
          }

    // ***************************************
          case 22:                 // *******double_normal_plateau (similar to CASAL)
            {peak2=sp(1)+ (0.99*len_bins(nlength)-sp(1))/(1.+mfexp(-sp(2))); upselex=mfexp(sp(3)); downselex=mfexp(sp(4));
             for(j=1;j<=nlength;j++)
             {
              t1=len_bins_m(j);
              if(t1<sp(1))
                {sel(j)=mfexp(-square(t1-sp(1))/upselex);}
              else if (t1<=peak2)
                {sel(j)=1.0;}
              else
                {sel(j)=mfexp(-square(t1-peak2)/downselex);}
             }
             break;
            }

// ***************************************
//          case 23:                 // *******double_normal_plateau (similar to CASAL)
//          {
//            break;
//          }

// ***************************************
          case 23:                 // *******double_normal_plateau where final value can be greater than 1.0
                                   //  cannot be used with Pope's because can cause selex to be >1.0
          {
//            if(gg==1 || seltype(f,3)<3)
          if(seltype(f,3)<3 || (gg==1 && seltype(f,3)==3) || (gg==2 && seltype(f,3)==4))
            {peak=sp(1); upselex=mfexp(sp(3)); downselex=mfexp(sp(4)); final=sp(6);}
            else
            {   // offset male parameters if seltype(f,3)==3
              peak=sp(1)+sp(Maleselparm(f));
              upselex=mfexp(sp(3)+sp(Maleselparm(f)+1));
              downselex=mfexp(sp(4)+sp(Maleselparm(f)+2));
              if(sp(6)>-999.) final=sp(6)+sp(Maleselparm(f)+3);
            }

            if(sp(5)<-1000.)
            {
              j1=-1001-int(value(sp(5)));      // selex is nil thru bin j1, so set sp(5) equal to first bin with selex (e.g. -1002 to start selex at bin 2)
              sel(1,j1)=1.0e-06;
            }
            else
            {
              j1=startbin-1;                // start selex at bin equal to min sizecomp databin  (=j1+1)
              if(sp(5)>-999)
              {
                point1=1.0/(1.0+mfexp(-sp(5)));
              t1min=mfexp(-(square(len_bins_m(startbin)-peak)/upselex));  // fxn at first bin
              }
            }
            if(sp(6)<-1000.)
            {
              j2=-1000-int(value(sp(6))); // selex is constant beyond this sizebin, so set sp(6) equal to last bin with estimated selex
            }
            else
            {j2=nlength;}
            peak2=peak+binwidth2+ (0.99*len_bins_m(j2)-peak-binwidth2)/(1.+mfexp(-sp(2)));
            if(sp(6)>-999)
            {
              point2=final;
              t2min=mfexp(-(square(len_bins_m(j2)-peak2)/downselex));  // fxn at last bin
            }
            for(j=j1+1;j<=j2;j++)
            {
              t1=len_bins_m(j)-peak;  t2=len_bins_m(j)-peak2;
              join1=1.0/(1.0+mfexp(-(20.*t1/(1.0+fabs(t1)))));  //  note the logit transform on t1 causes range of mfexp to be over -20 to 20
              join2=1.0/(1.0+mfexp(-(20.*t2/(1.0+fabs(t2)))));
              if(sp(5)>-999)
                {asc=point1+(1.0-point1)*(mfexp(-square(t1)/upselex)-t1min)/(1.0-t1min);}
              else
                {asc=mfexp(-square(t1)/upselex);}
              if(sp(6)>-999)
                {dsc=1.0+(point2-1.0)*(mfexp(-square(t2)/downselex)-1.0    )/(t2min-1.0);}
              else
                {dsc=mfexp(-square(t2)/downselex);}
              sel(j)=asc*(1.0-join1)+join1*(1.0-join2+dsc*join2);
            }
            if(startbin>1 && sp(5)>=-1000.)
            {
              for(j=1;j<=startbin-1;j++)
              {
                sel(j)=square(len_bins_m(j)/len_bins_m(startbin))*sel(startbin);
              }
            }

            if(j2<nlength) {sel(j2+1,nlength)=sel(j2);}
            break;
          }


// ***************************************
          case 24:                 // *******double_normal_plateau
          {
//            if(gg==1 || seltype(f,3)<3)
          if(seltype(f,3)<3 || (gg==1 && seltype(f,3)==3) || (gg==2 && seltype(f,3)==4))
            {peak=sp(1); upselex=mfexp(sp(3)); downselex=mfexp(sp(4)); final=sp(6); Apical_Selex=1.;}
            else
            {   // offset male parameters if seltype(f,3)==3, female parameters if seltype(f,3)==4
              peak=sp(1)+sp(Maleselparm(f));
              upselex=mfexp(sp(3)+sp(Maleselparm(f)+1));
              downselex=mfexp(sp(4)+sp(Maleselparm(f)+2));
              if(sp(6)>-999.) final=sp(6)+sp(Maleselparm(f)+3);
              Apical_Selex=sp(Maleselparm(f)+4);
            }

            if(sp(5)<-1000.)
            {
              j1=-1001-int(value(sp(5)));      // selex is nil thru bin j1, so set sp(5) equal to first bin with selex (e.g. -1002 to start selex at bin 2)
              sel(1,j1)=1.0e-06;
            }
            else
            {
              j1=startbin-1;                // start selex at bin equal to min sizecomp databin  (=j1+1)
              if(sp(5)>-999)
              {
                point1=1.0/(1.0+mfexp(-sp(5)));
              t1min=mfexp(-(square(len_bins_m(startbin)-peak)/upselex));  // fxn at first bin
              }
            }
            if(sp(6)<-1000.)
            {
              j2=-1000-int(value(sp(6))); // selex is constant beyond this sizebin, so set sp(6) equal to last bin with estimated selex
            }
            else
            {j2=nlength;}
            peak2=peak+binwidth2+ (0.99*len_bins_m(j2)-peak-binwidth2)/(1.+mfexp(-sp(2)));
            if(sp(6)>-999)
            {
              point2=1.0/(1.0+mfexp(-final));
              t2min=mfexp(-(square(len_bins_m(j2)-peak2)/downselex));  // fxn at last bin
            }
            for(j=j1+1;j<=j2;j++)
            {
              t1=len_bins_m(j)-peak;  t2=len_bins_m(j)-peak2;
              join1=1.0/(1.0+mfexp(-(20.*t1/(1.0+fabs(t1)))));  //  note the logit transform on t1 causes range of mfexp to be over -20 to 20
              join2=1.0/(1.0+mfexp(-(20.*t2/(1.0+fabs(t2)))));
              if(sp(5)>-999)
                {asc=point1+(Apical_Selex-point1)*(mfexp(-square(t1)/upselex)-t1min)/(1.0-t1min);}
              else
                {asc=Apical_Selex*mfexp(-square(t1)/upselex);}
              if(sp(6)>-999)
                {dsc=Apical_Selex+(point2-Apical_Selex)*(mfexp(-square(t2)/downselex)-1.0    )/(t2min-1.0);}
              else
                {dsc=Apical_Selex*mfexp(-square(t2)/downselex);}
              sel(j)=asc*(1.0-join1)+join1*(Apical_Selex*(1.0-join2)+dsc*join2);
            }
            if(startbin>1 && sp(5)>=-1000.)
            {
              for(j=1;j<=startbin-1;j++)
              {
                sel(j)=square(len_bins_m(j)/len_bins_m(startbin))*sel(startbin);
              }
            }

            if(j2<nlength) {sel(j2+1,nlength)=sel(j2);}
            break;
          }

//*******************************************
          case 25:                 // *******exponential-logistic
          {
            peak = len_bins_m(1) + sp(2)*(len_bins_m(nlength)-len_bins_m(1));
            for(j=1;j<=nlength;j++)
              {sel(j) = mfexp(sp(3)*sp(1)*(peak-len_bins_m(j)))/(1.0-sp(3)*(1.0-mfexp(sp(1)*(peak-len_bins_m(j)))));}
            break;
          }

//*******************************************
          case 27: // cubic spline for selex at length
          {
            int j1;
            int j2;

            j=1;
            k=seltype(f,4);  // n points to include in cubic spline
            for(i=1;i<=k;i++)
            {
              splineX(i)=value(sp(i+3)); // "value" required to avoid error, but values should be always fixed anyway
              splineY(i)=sp(i+3+k);
            }
            z=nlength;
            while(len_bins_m(z)>splineX(k)) {z--;}
            j2=z+1;  //  first size bin beyond last node
            vcubic_spline_function splinefn=vcubic_spline_function(splineX(1,k),splineY(1,k),sp(2),sp(3));
            tempvec_l = splinefn(len_bins_m);  // interpolate selectivity at the mid-point of each population size bin
            temp=max(tempvec_l(1,j2));
            tempvec_l-=temp;  // rescale to get max of 0.0
            tempvec_l(j2+1,nlength) = tempvec_l(j2);  //  set constant above last node
            sel = mfexp(tempvec_l);
            break;
          }

          case 30:
          {
            sel=1.0;
            break;
          }
          case 31:
          {
            sel=1.0;
            break;
          }
          case 32:
          {
            sel=1.0;
            break;
          }
          case 33:
          {
            sel=1.0;
            break;
          }
          case 34:
          {
            sel=1.0;
            break;
          }

          default:   //  seltype not found.  But really need this check earlier when the N selex parameters are being processed.
          {
            N_warn++; cout<<"Critical error, see warning"<<endl; warning<<"Size_selex option not valid "<<seltype(f,1)<<endl; exit(1);
            break;
          }

          }
          sel_l(y,f,gg)=sel;    // Store size-selex in year*type array
        }  // end direct calc of selex from parameters

// Do male relative to female selex
        if(gg==2)         // males exist and am now in the male loop
        {
         if(seltype(f,1)==4)
           {sel_l(y,f,2)=0.;}  // set males to zero for spawning biomass
         else if(seltype(f,1)==5)    // set males equal to mirrored males
         {
          i=int(value(sp(1)));  if(i<=0) i=1;
          j=int(value(sp(2)));  if(j<=0) j=nlength;
          sel_l(y,f,2)(i,j)=sel_l(y,seltype(f,4),2)(i,j);
         }
         else if(seltype(f,1)==15)    // set males equal to mirrored males
         {
          sel_l(y,f,2)=sel_l(y,seltype(f,4),2);
         }
         else if(seltype(f,3)==1 || seltype(f,3)==2)   // do gender selex as offset
         {
//           k=seltype_Nparam(seltype(f,1)) + 1;
//           if(seltype(f,2)>0) k+=seltype(f,2)*4;   // first gender offset parm (skip over the retention parameters
//           if(seltype(f,1)==6) k += seltype(f,4);    // for non-parametric in which N parm is stored in special column
           k=Maleselparm(f);
           temp=sp(k);
           temp1=1.;
           switch(seltype(f,3))
           {
             case 1:
             {                            // do males relative to females
               for(j=1;j<=nlength;j++)
               {
                 if(len_bins_m(j)<=temp)
                   {sel(j)*=mfexp(sp(k+1)+(len_bins_m(j)-minL_m)/(temp-minL_m) * (sp(k+2)-sp(k+1)) );}
                 else
                   {sel(j)*=mfexp(sp(k+2)+(len_bins_m(j)-temp) /(maxL-temp)  * (sp(k+3)-sp(k+2)) );}
//                 if(sel(j)>temp1) temp1=sel(j);
               }
               sel_l(y,f,2)=sel;
               tempvec_a(1)=max(sel_l(y,f,1));
               tempvec_a(2)=max(sel_l(y,f,2));
               temp1=max(tempvec_a(1,2));
               sel_l(y,f) /=temp1;
               break;
             }
             case 2:
             {                   //  do females relative to males
               sel_l(y,f,2)=sel;
               for(j=1;j<=nlength;j++)
               {
                 if(len_bins_m(j)<=temp)
                   {sel(j)*=mfexp(sp(k+1)+(len_bins_m(j)-minL_m)/(temp-minL_m) * (sp(k+2)-sp(k+1)) );}
                 else
                   {sel(j)*=mfexp(sp(k+2)+(len_bins_m(j)-temp) /(maxL-temp)  * (sp(k+3)-sp(k+2)) );}
//                 if(sel(j)>temp1) temp1=sel(j);
               }
               sel_l(y,f,1)=sel;
               tempvec_a(1)=max(sel_l(y,f,1));
               tempvec_a(2)=max(sel_l(y,f,2));
               temp1=max(tempvec_a(1,2));
               sel_l(y,f)/=temp1;
               break;
             }
           }  // end switch
         }  // end do gender selex as offset from other gender
         else if(seltype(f,3)!=3 && seltype(f,3)!=4)
         {
           sel_l(y,f,2)=sel;   // set males = females, but also will set male selex to calculated selex if seltype(f,3)==3
         }
        }  // end doing males
      }  // end loop of genders
      if(docheckup==1) echoinput<<"length selex: "<<endl<<sel_l(y,f)<<endl;

  //  SS_Label_Info_22.5.1 #Calculate discmort
  // discmort is the size-specific fraction of discarded fish that die
  //  discmort2 is size-specific fraction that die from being retained or are dead discard
  //   = elem_prod(sel,(retain + (1-retain)*discmort)) */

      if(seltype(f,2)==0)  //  no discard, all retained
      {
        retain(y,f)=1.0;
        sel_l_r(y,f)=sel_l(y,f);
        discmort(y,f)=1.0;
        discmort2(y,f)=sel_l(y,f);
        if(gender==2)
        {
          discmort_M=1.0;
          retain_M=1.0;
        }
      }
      else if(seltype(f,2)==3)  // none retained; all dead
      {
        retain(y,f)=0.0;
        discmort(y,f)=1.0;
        sel_l_r(y,f)=0.0;
        discmort2(y,f)=sel_l(y,f);
        if(gender==2)
        {
          discmort_M=1.0;
          retain_M=0.0;
        }
      }
      else
      {
        if(seltype(f,2)<0)  // mirror
        {
          k=-seltype(f,2);
          retain(y,f)=retain(y,k);
          discmort(y,f)=discmort(y,k);
          if(seltype(k,2)==1)
          {
            discmort2(y,f)=sel_l(y,f);  //  all selected fish are dead;  this statement does both genders implicitly
          }
          else
          {
            discmort2(y,f,1)=elem_prod(sel_l(y,f,1), retain(y,f)(1,nlength) + elem_prod((1.-retain(y,f)(1,nlength)),discmort(y,f)(1,nlength)) );
          }
        }
        else  // use parameters
        {
          k=RetainParm(f);
          temp=1.-sp(k+2);
          temp1=1.-posfun(temp,0.0,CrashPen);
          retain(y,f)=temp1/(1.+mfexp(-(len_bins_m2-(sp(k)+male_offset*sp(k+3)))/sp(k+1)));  // males are at end of vector, so automatically get done
          if(docheckup==1) echoinput<<" retention parms "<<sp(k)<<" "<<sp(k+1)<<" "<<sp(k+3)<<" "<<temp1<<endl<<"maleoff "<<male_offset<<endl;
          if(docheckup==1) echoinput<<"lenbins "<<len_bins_m2<<endl;
          if(docheckup==1) echoinput<<"retention "<<retain(y,f)<<endl;
  
          if(seltype(f,2)==1)  // all discards are dead
          {
            discmort(y,f)=1.0;
            discmort2(y,f)=sel_l(y,f);  //  all selected fish are dead;  this statement does both genders implicitly
          }
          else
          {
            k+=4;  // first discard mortality parm
            temp=1.-sp(k+2);
            temp1=posfun(temp,0.0,CrashPen);
            discmort(y,f)=(1.-temp1/(1+mfexp(-(len_bins_m2-(sp(k)+male_offset*sp(k+3)))/sp(k+1))));  // males are at end of vector, so automatically get done
            if(docheckup==1) echoinput<<"discmort "<<discmort(y,f)<<endl;
            discmort2(y,f,1)=elem_prod(sel_l(y,f,1), retain(y,f)(1,nlength) + elem_prod((1.-retain(y,f)(1,nlength)),discmort(y,f)(1,nlength)) );
          }
        }
  
        sel_l_r(y,f,1)=elem_prod(sel_l(y,f,1),retain(y,f)(1,nlength));
        if(gender==2)
        {
          discmort_M.shift(nlength1)=discmort(y,f)(nlength1,nlength2);
          retain_M.shift(nlength1)=retain(y,f)(nlength1,nlength2);
          sel_l_r(y,f,2)=elem_prod(sel_l(y,f,2),retain_M.shift(1));
          discmort2(y,f,2)=elem_prod(sel_l(y,f,2), retain_M.shift(1) + elem_prod((1.-retain_M.shift(1)),discmort_M.shift(1)) );  // V3.21f
        }
      }
      if(docheckup==1) echoinput<<"sel-len-retained "<<sel_l_r(y,f)<<endl<<" dead "<<discmort2(y,f)<<endl;

      }  //  end loop of fleets for size selex and retention and discard mortality

// ----------------------------------------------------------------------------------
      else   //    do age-selex
      {
        if(docheckup==1 && N_selparmvec(f)>0) echoinput<<" age selex parms: "<<sp(1,N_selparmvec(f))<<endl;
        for(gg=1;gg<=gender;gg++)
        {
          if(gg==1 || (gg==2 && seltype(f,3)>=3))  //  in age selex
          {
            switch(seltype(f,1))
            {

              case 0:   // constant age-specific selex
              {sel_a(y,fs,1)(0,nages)=1.00; break;}

              case 10:   // constant age-specific selex
              {sel_a(y,fs,1)(1,nages)=1.00; break;}                 //  change with version 1.21 to omit age 0

              case 11:   // selex=1.0 within a range of ages
              {
                a=int(value(sp(2)));
                if(a>nages) {a=nages;}
                sel_a(y,fs,1)(int(value(sp(1))),a)=1.;
                break;
              }

              case 12:                      // ***********   logistic    ProgLabel_055
              { sel_a(y,fs,1) = 1/(1+mfexp(neglog19*(r_ages-sp(1))/sp(2))); break;}

              case 13:                 // *******double logistic     ProgLabel_056
                                       // 1=peak, 2=init,  3=infl,  4=slope, 5=final, 6=infl2, 7=slope2
              {
                t1=0.+(1./(1.+mfexp(-sp(3))))*(sp(1)-0.);    // INFL
                t1min=1./(1.+mfexp(-sp(4)*(0.-t1)))*0.9999999;  // asc value at minage
                t1max=1./(1.+mfexp(-sp(4)*(sp(1)-t1)))*1.00001;  // asc value at peak
                t1power=log(0.5)/log((0.5-t1min)/(t1max-t1min));

                t2=(sp(1)+sp(8))+(1./(1.+mfexp(-sp(6))))*(r_ages(nages)-(sp(1)+sp(8)));    // INFL
                t2min=1./(1.+mfexp(-sp(7)*(sp(1)+sp(8)-t2)))*0.9999;  // asc value at peak+
                t2max=1./(1.+mfexp(-sp(7)*(r_ages(nages)-t2)))*1.00001;  // asc value at maxage
                t2power=log(0.5)/log((0.5-t2min)/(t2max-t2min));
                final=1./(1.+mfexp(-sp(5)));
                k1=int(value(sp(1))); k2=int(value(sp(1)+sp(8)));

                for (a=0; a<=nages; a++)  //calculate the value over ages
                {
                  if (a < k1) // ascending limb
                  {
                    sel_a(y,fs,1,a) = sp(2) + (1. - sp(2)) *
                    pow((( 1./(1.+mfexp(-sp(4)*(r_ages(a)-t1))) -t1min ) / (t1max-t1min) ),t1power);
                  }
                  else if (a > k2) // descending limb
                  {
                    sel_a(y,fs,1,a) = 1. + (final - 1.) *
                    pow((( 1./(1.+mfexp(-sp(7)*(r_ages(a)-t2))) -t2min ) / (t2max-t2min) ),t2power);
                  }
                  else // at the peak
                  { sel_a(y,fs,1,a) = 1.0;}
                }   // end age loop
                break;
              }    // end double logistic

              case 14:                  // separate parm for each age
            {
             temp=9.-max(sp(1,nages+1));  //  this forces at least one age to have selex weight equal to 9
             for (a=0;a<=nages;a++)
             {
              if(sp(a+1)>-999)
              {sel_a(y,fs,1,a) = 1./(1.+mfexp(-(sp(a+1)+temp)));}
              else
              {sel_a(y,fs,1,a) = sel_a(y,fs,1,a-1);}
              }
              break;
            }

              case 15:                //  mirror another selex
                                      // must refer to a lower numbered type (f)
            {
              sel_a(y,fs)=sel_a(y,seltype(f,4));
              break;
            }

              case 16:                  // Coleraine - Gaussian
            {
             t1 = 1/(1+mfexp(-sp(1)))*nages;
             for (a=0;a<=nages;a++)
             {
              if(a<t1)
              {sel_a(y,fs,1,a) = mfexp(-square(r_ages(a)-t1)/mfexp(sp(2)));}
              else
              {sel_a(y,fs,1,a)=1.0;}
             }
             break;
            }

              case 17:                  // separate parm for each age as random walk
            {
              lastsel=0.0;  //  value is the change in log(selex);  this is the reference value for age 0
              tempvec_a=-999.;
              tempvec_a(0)=0.0;   //  so do not try to estimate the first value
              int lastage;
              if(seltype(f,4)==0)
              {lastage=nages;}
              else
              {lastage=abs(seltype(f,4));}

              for(a=1;a<=lastage;a++)
              {
                if(sp(a+1)>-999.) {lastsel=sp(a+1);}  //  so, lastsel stays constant until changed, so could create a linear change in ln(selex)
                                                      // use of (a+1) is because the first element, sp(1), is for age zero
                tempvec_a(a)=tempvec_a(a-1)+lastsel;   // cumulative log(selex)
              }
              temp=max(tempvec_a);   //  find max so at least one age will have selex=1.
              sel_a(y,fs,1)=mfexp(tempvec_a-temp);
              a=0;
              while(sp(a+1)==-1000)  //  reset range of young ages to selex=0.0
              {
                sel_a(y,fs,1,a)=0.0;
                a++;
              }
              if(lastage<nages)
              {
                for(a=lastage+1;a<=nages;a++)
                {
                  if(seltype(f,4)>0)
                  {sel_a(y,fs,1,a)=sel_a(y,fs,1,a-1);}
                  else
                  {sel_a(y,fs,1,a)=0.0;}
                }
              }
              break;
            }

              case 18:                 // *******double logistic with smooth transition
                                       // 1=peak, 2=init,  3=infl,  4=slope, 5=final, 6=infl2, 7=slope2
            {
             t1=0.+(1./(1.+mfexp(-sp(3))))*(sp(1)-0.);    // INFL
             t1min=1./(1.+mfexp(-sp(4)*(0.-t1)))*0.9999;  // asc value at minsize
             t1max=1./(1.+mfexp(-sp(4)*(sp(1)-t1)))*1.00001;  // asc value at peak
             t1power=log(0.5)/log((0.5-t1min)/(t1max-t1min));

             t2=(sp(1)+sp(8))+(1./(1.+mfexp(-sp(6))))*(r_ages(nages)-(sp(1)+sp(8)));    // INFL
             t2min=1./(1.+mfexp(-sp(7)*(sp(1)+sp(8)-t2)))*0.9999;  // asc value at peak+
             t2max=1./(1.+mfexp(-sp(7)*(r_ages(nages)-t2)))*1.00001;  // asc value at maxage
             t2power=log(0.5)/log((0.5-t2min)/(t2max-t2min));
             final=1./(1.+mfexp(-sp(5)));
             for (a=0; a<=nages; a++)  //calculate the value over ages
             {
              sel_a(y,fs,1,a) =
                (
                (
                (sp(2) + (1.-sp(2)) *
                 pow((( 1./(1.+mfexp(-sp(4)*(r_ages(a)-t1)))-t1min)/ (t1max-t1min)),t1power))
                /(1.0+mfexp(30.*(r_ages(a)-sp(1))))  // scale ascending side
                +
                1./(1.+mfexp(-30.*(r_ages(a)-sp(1))))   // flattop, with scaling
                )
                /(1.+mfexp( 30.*(r_ages(a)-(sp(1)+sp(8)))))    // scale combo of ascending and flattop
                +
                (1. + (final - 1.) *
                 pow(fabs((( 1./(1.+mfexp(-sp(7)*(r_ages(a)-t2))) -t2min ) / (t2max-t2min) )),t2power))
                /(1.+mfexp( -30.*(r_ages(a)-(sp(1)+sp(8)))))    // scale descending
                );
             }   // end age loop
             break;
            }    // end double logistic with smooth transition

            case 19:                     // ***********   old double logistic   TAG_069
            {
              k1=int(value(sp(5)));
              sel_a(y,fs,1)(k1,nages) =   elem_prod((1./(1.+mfexp(-sp(2)*(r_ages(k1,nages)-sp(1)) ))),
                                                   (1.-1./(1.+mfexp(-sp(4)*(r_ages(k1,nages)-(sp(1)*sp(6)+sp(3))) ))) );
              sel_a(y,fs,1)(k1,nages) /= max(sel_a(y,fs,1)(k1,nages));
              if(k1>0) sel_a(y,fs,1)(0,k1-1)=1.0e-6;
              break;
            }

      // ***************************************
            case 20:                 // *******double_normal_plateau
            {
              if(seltype(f,3)<3 || (gg==1 && seltype(f,3)==3) || (gg==2 && seltype(f,3)==4))
              {peak=sp(1); upselex=mfexp(sp(3)); downselex=mfexp(sp(4)); final=sp(6); Apical_Selex=1.0;}
              else
              {   // offset male parameters if seltype(f,3)==3
                peak=sp(1)+sp(Maleselparm(f));
                upselex=mfexp(sp(3)+sp(Maleselparm(f)+1));
                downselex=mfexp(sp(4)+sp(Maleselparm(f)+2));
                if(sp(6)>-999.) final=sp(6)+sp(Maleselparm(f)+3);
                Apical_Selex=sp(Maleselparm(f)+4);
              }
              if(sp(5)<-1000.)
              {
                j=-1001-int(value(sp(5)));      // selex is nil thru age j, so set sp(5) equal to first age with selex (e.g. -1002 to start selex at age 2)
                sel_a(y,fs,gg)(0,j)=1.0e-06;
              }
              else
              {
                j=-1;                // start selex at age 0
                if(sp(5)>-999)
                {
                  point1=1./(1.+mfexp(-sp(5)));
                  t1min=mfexp(-(square(0.-peak)/upselex));  // fxn at first bin
                }
              }
              if(sp(6)<-1000.)
              {
                j2=-1000-int(value(sp(6))); // selex is constant beyond this age, so set sp(6) equal to last age with estimated selex
                                              //  (e.g. -1008 to be constant beyond age 8)
              }
              else
              {j2=nages;}

              peak2=peak+1.+(0.99*r_ages(j2)-peak-1.)/(1.+mfexp(-sp(2)));        // note, this uses age=j2 as constraint on range of "peak2"
//              peak2=peak+.1+(0.99*r_ages(j2)-peak-.1)/(1.+mfexp(-sp(2)));        // note, this uses age=j2 as constraint on range of "peak2"
              if(sp(6)>-999)
              {
                point2=1./(1.+mfexp(-final));
                t2min=mfexp(-(square(r_ages(nages)-peak2)/downselex));  // fxn at last bin
              }

              for(a=j+1;a<=j2;a++)
              {
                t1=r_ages(a)-peak;  t2=r_ages(a)-peak2;
                join1=1./(1.+mfexp(-(20./(1.+fabs(t1)))*t1));
                join2=1./(1.+mfexp(-(20./(1.+fabs(t2)))*t2));
                if(sp(5)>-999)
                  {asc=point1+(Apical_Selex-point1)*(mfexp(-square(t1)/upselex  )-t1min)/(1.-t1min);}
                else
                  {asc=Apical_Selex*mfexp(-square(t1)/upselex);}
                if(sp(6)>-999)
                  {dsc=Apical_Selex+(point2-Apical_Selex)*(mfexp(-square(t2)/downselex)-1.    )/(t2min-1.);}
                else
                  {dsc=Apical_Selex*mfexp(-square(t2)/downselex);}
                sel_a(y,fs,gg,a)=asc*(1.-join1)+join1*(Apical_Selex*(1.-join2)+dsc*join2);
              }
              if(j2<nages) {sel_a(y,fs,gg)(j2+1,nages)=sel_a(y,fs,gg,j2);}
              break;
            }

      // ***************************************
            case 26:                 // *******exponential-logistic in age
            {
              peak = r_ages(0) + sp(2)*(r_ages(nages)-r_ages(0));
              for(a=0;a<=nages;a++)
                {sel_a(y,fs,1,a) = mfexp(sp(3)*sp(1)*(peak-r_ages(a)))/(1.0-sp(3)*(1.0-mfexp(sp(1)*(peak-r_ages(a)))));}
              break;
            }

//*******************************************
          case 27: // cubic spline for selex at age
          {
            k=seltype(f,4);  // n points to include in cubic spline
            for(i=1;i<=k;i++)
            {
              splineX(i)=value(sp(i+3)); // "value" required to avoid error, but values should be always fixed anyway
              splineY(i)=sp(i+3+k);
            }
            z=nages;
            while(r_ages(z)>splineX(k)) {z--;}
            j2=z+1;  //  first age beyond last node
            vcubic_spline_function splinefn=vcubic_spline_function(splineX(1,k),splineY(1,k),sp(2),sp(3));
            tempvec_a= splinefn(r_ages);  // interpolate selectivity at each age
            temp=max(tempvec_a(0,j2));
            tempvec_a-=temp;  // rescale to get max of 0.0
            tempvec_a(j2+1,nages) = tempvec_a(j2);  //  set constant above last node
            sel_a(y,fs,1)=mfexp(tempvec_a);
            break;
          }

          default:   //  seltype not found.  But really need this check earlier when the N selex parameters are being processed.
          {
            N_warn++; cout<<"Critical error, see warning"<<endl; warning<<"Age_selex option not valid "<<seltype(f,1)<<endl; exit(1);
            break;
          }

            }  // end last age selex pattern
         }  // end direct calc of selex from parameters
  // ***************************************
  // Do male relative to female selex
          if(gg==2)         // males exist
          {
            if(seltype(f,3)==1 || seltype(f,3)==2)   // do gender selex as offset
            {
              k=Maleselparm(f);   // first male parm
              temp=sp(k)-0.00001;
              temp1=1.;
              switch(seltype(f,3))
              {
                case 1:
                {                       // do males relative to females
                  for(a=0;a<=nages;a++)   //
                  {
                    if(r_ages(a)<=temp)
                    {sel_a(y,fs,2,a)=sel_a(y,fs,1,a)*mfexp(sp(k+1)+(r_ages(a)-0.)   /(temp-0.)   * (sp(k+2)-sp(k+1)) );}
                    else
                    {sel_a(y,fs,2,a)=sel_a(y,fs,1,a)*mfexp(sp(k+2)+(r_ages(a)-temp) /(double(nages)-temp) * (sp(k+3)-sp(k+2)) );}
      //              if(sel_a(y,fs,2,a)>temp1) temp1=sel_a(y,fs,2,a);
                  }
                  tempvec_a(1)=max(sel_a(y,fs,1));
                  tempvec_a(2)=max(sel_a(y,fs,2));
                  temp1=max(tempvec_a(1,2));
                  sel_a(y,fs)/=temp1;
                  break;
                }
                case 2:
                {                   //  do females relative to males
                  sel_a(y,fs,2)=sel_a(y,fs,1);
                  for(a=0;a<=nages;a++)   //
                  {
                    if(r_ages(a)<=temp)
                      {sel_a(y,fs,1,a)=sel_a(y,fs,2,a)*mfexp(sp(k+1)+(r_ages(a)-0.)   /(temp-0.)   * (sp(k+2)-sp(k+1)) );}
                    else
                      {sel_a(y,fs,1,a)=sel_a(y,fs,2,a)*mfexp(sp(k+2)+(r_ages(a)-temp) /(double(nages)-temp) * (sp(k+3)-sp(k+2)) );}
      //              if(sel_a(y,fs,1,a)>temp1) temp1=sel_a(y,fs,1,a);
                  }
      //            sel_a(y,fs)/=temp1;
                  tempvec_a(1)=max(sel_a(y,fs,1));
                  tempvec_a(2)=max(sel_a(y,fs,2));
                  temp1=max(tempvec_a(1,2));
                  sel_a(y,fs)/=temp1;

                  break;
                }
              }
            }
            else if(seltype(f,3)!=3 && seltype(f,3)!=4 && seltype(f,1)!=15)
            {
              sel_a(y,fs,2)=sel_a(y,fs,1);
            }   // set males = females
          }
        }  //  end gender
        if(docheckup==1) echoinput<<" age selex: "<<endl<<sel_a(y,fs)<<endl;
      }  // end recalc of the selex
    }

    else
    { // carryover from previous year
      if(f<=Ntypes)
      {
        sel_l(y,f)=sel_l(y-1,f);   // this does both genders
        sel_l_r(y,f)=sel_l_r(y-1,f);
        retain(y,f)=retain(y-1,f);
        discmort(y,f)=discmort(y-1,f);
        discmort2(y,f)=discmort2(y-1,f);
      }
      else  // age
      {
        sel_a(y,fs)=sel_a(y-1,fs);  // does both genders
      }
    }

    if(f<=Ntypes)
    {
      Ip+=N_selparmvec(f);
    }
    else
    {
      Ip+=N_selparmvec(f);
    }
  }  //  end fleet loop for selectivity
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 23 get_initial_conditions */
FUNCTION void get_initial_conditions()
  {
  if(do_once==1) echoinput<<"start initial conditions "<<endl;
  catch_fleet.initialize();
  annual_catch.initialize();
  annual_F.initialize();
    if(F_Method==2)
    {
      if(current_phase()>=F_setup(2) || (readparfile==1 && current_phase()<=1)) //  set Hrate = Frate parameters on first call if readparfile=1, or for advanced phases
      {
        g=0;
        for(f=1;f<=Nfleet;f++)
        for(t=styr;t<=TimeMax;t++)
        {
          g++;
          Hrate(f,t)=F_rate(g);
        }
      }
    }

// Initial size-at-age distribution             ProgLabel_091
  y=styr;
  yz=styr;
//  t_base=styr+(y-styr)*nseas-1;
  t_base=styr-1;
  if(MG_active(0)>0 || save_for_report>0) get_MGsetup();
  if(do_once==1) cout<<" MG setup OK "<<endl;
  if(MG_active(2)>0) get_growth1();   // seasonal effects and CV
  if(do_once==1) cout<<" growth OK"<<endl;
  if(MG_active(2)>0 || save_for_report>0)
  {
    get_growth2();   // do growth before natmort in case lorensen M is used
  }
  if(do_once==1) cout<<" growth-2 OK"<<endl;
  if(MG_active(1)>0) get_natmort();
  if(do_once==1) cout<<" natmort OK"<<endl;
  if(MG_active(3)>0) get_wtlen();
  if(MG_active(4)>0) get_recr_distribution();
  if(MG_active(5)>0) get_migration();
  if(do_once==1) cout<<" migr OK"<<endl;
  if(Use_AgeKeyZero>0)
  {
    if(MG_active(6)>0) get_age_age(Use_AgeKeyZero,AgeKey_StartAge,AgeKey_Linear1,AgeKey_Linear2); //  call function to get the age_age key
    if(do_once==1) cout<<" ageerr_key OK"<<endl;
  }

  if(save_for_report>0) get_saveGparm();
  if(do_once==1) cout<<" growth OK, ready to call selex "<<endl;
  get_selectivity();
  if(do_once==1) cout<<" selex OK, ready to call ALK and fishselex "<<endl;
   for(s=1;s<=nseas;s++)
     {
      smid=s+nseas;
     t = styr+s-1;
     if(MG_active(2)>0 || MG_active(3)>0 || save_for_report>0 || WTage_rd>0)  Make_AgeLength_Key();
     for(g=1;g<=gmorph;g++)
       if(use_morph(g)>0)
       {
       Make_FishSelex();
       }
     }
  virg_fec = fec;

  if(do_once==1) cout<<" ready for virgin age struc "<<endl;
//Virgin age structure                    ProgLabel_092
  eq_yr=styr-2;
  Recr_virgin=mfexp(SR_parm(1));
  exp_rec(eq_yr,1)=Recr_virgin;  //  expected Recr from s-r parms
  exp_rec(eq_yr,2)=Recr_virgin; 
  exp_rec(eq_yr,3)=Recr_virgin; 
  exp_rec(eq_yr,4)=Recr_virgin;
  
   bio_yr=styr;
   Fishon=0;
   equ_Recr=Recr_virgin;

   Do_Equil_Calc();                      //  call function to do equilibrium calculation
   SPB_virgin=SPB_equil;
      Mgmt_quant(1)=SPB_equil;
      if(Do_Benchmark>0)
      {
        Mgmt_quant(2)=totbio;
        Mgmt_quant(3)=smrybio;
        Mgmt_quant(4)=Recr_virgin;
      }

   SPB_pop_gp(eq_yr)=SPB_equil_pop_gp;   // dimensions of pop x N_GP
   if(Hermaphro_Option>0) MaleSPB(eq_yr)=MaleSPB_equil_pop_gp;
   SPB_yr(eq_yr)=SPB_equil;
   SPB_B_yr(eq_yr)=equ_mat_bio;  
   SPB_N_yr(eq_yr)=equ_mat_num;  
   t=styr-2*nseas-1;
   for(p=1;p<=pop;p++)
   for(g=1;g<=gmorph;g++)
   for(s=1;s<=nseas;s++)
     {natage(t+s,p,g)(0,nages)=value(equ_numbers(s,p,g)(0,nages));}


   //****************************************************************************
   //Initial Age Structure; in equilibrium with init_equ_catch    ProgLabel_093
  if(do_once==1) cout<<" ready for initial age struc "<<endl;
  eq_yr=styr-1;
  bio_yr=styr;
  t=styr-nseas;
  if(fishery_on_off==1) {Fishon=1;} else {Fishon=0;}
  for(f=1;f<=Nfleet;f++) {Hrate(f,t) = init_F(f);}

  if(SR_parm_1(N_SRparm2-1,5)>-999)  //  using the PR_type as a flag
  {
//  do initial equilibrium while maintain R1 based on offset from Recr_virgin
   R1_exp=Recr_virgin;
   exp_rec(eq_yr,1)=R1_exp;
   if(SR_env_target==2) {R1_exp*=mfexp(SR_parm(N_SRparm2-2)* env_data(eq_yr,SR_env_link));}
   exp_rec(eq_yr,2)=R1_exp;
   exp_rec(eq_yr,3)=R1_exp;
   R1 = R1_exp*mfexp(SR_parm(N_SRparm2-1));
   exp_rec(eq_yr,4)=R1;

   equ_Recr=R1;   //  value to use in Do_Equil_Calc
   Do_Equil_Calc();
   CrashPen += Equ_penalty;
  }
  else
  {
//  do initial equilibrium with R1 based on offset from spawner-recruitment curve, using same approach as the benchmark calculations
//  first get SPR for this init_F
  equ_Recr=Recr_virgin;
  Do_Equil_Calc();
  CrashPen += Equ_penalty;

  SPR_temp=SPB_equil/equ_Recr;  //  spawners per recruit at initial F
//  next the rquilibrium SSB and recruitment from SPR_temp
  Get_EquilCalc = Equil_Spawn_Recr_Fxn();  //  returns 2 element vector containing equilibrium biomass and recruitment at this SPR
  R1_exp=Get_EquilCalc(2);     //  set the expected recruitment equal to this equilibrium
  exp_rec(eq_yr,1)=R1_exp;
  if(SR_env_target==2) {R1_exp*=mfexp(SR_parm(N_SRparm2-2)* env_data(eq_yr,SR_env_link));}  //  adjust for environment
  exp_rec(eq_yr,2)=R1_exp;
  exp_rec(eq_yr,3)=R1_exp;

  equ_Recr = R1_exp*mfexp(SR_parm(N_SRparm2-1));  // apply R1 offset
  R1=equ_Recr;
  exp_rec(eq_yr,4)=equ_Recr;
  Do_Equil_Calc();  // calculated SPB_equil
  CrashPen += Equ_penalty;
  }

   SPB_pop_gp(eq_yr)=SPB_equil_pop_gp;   // dimensions of pop x N_GP
   if(Hermaphro_Option>0) MaleSPB(eq_yr)=MaleSPB_equil_pop_gp;
   SPB_yr(eq_yr)=SPB_equil;
   SPB_yr(styr)=SPB_equil;
   SPB_B_yr(eq_yr)=equ_mat_bio;  
   SPB_N_yr(eq_yr)=equ_mat_num;  
   SPB_B_yr(styr)=equ_mat_bio;  
   SPB_N_yr(styr)=equ_mat_num;  

   for(f=1;f<=Nfleet;f++)
   {
     if(catchunits(f)==1)
     {
      est_equ_catch(f)=equ_catch_fleet(f,2);
    }
    else
     {
      est_equ_catch(f)=equ_catch_fleet(f,5);
     }
    }

   t=styr-nseas-1;
   if(save_for_report>0)
     {
     for(f=1;f<=Nfleet;f++)
       {
        catch_fleet(t+1,f)=equ_catch_fleet(f);  // gets all 6 elements
       for(g=1;g<=gmorph;g++)
       for(s=1;s<=nseas;s++)
         {catage(t+s,f,g)=value(equ_catage(s,f,g)); }
       }
     }

   for(p=1;p<=pop;p++)
   for(g=1;g<=gmorph;g++)
   for(s=1;s<=nseas;s++)
     {natage(t+s,p,g)(0,nages)=equ_numbers(s,p,g)(0,nages);}

   natage(styr) = natage(t+1);
   if(docheckup==1) echoinput<<" init age comp"<<endl<<natage(styr)<<endl<<endl;


   // if recrdevs start before styr, then use them to adjust the initial agecomp
   //  apply a fraction of the bias adjustment, so bias adjustment gets less linearly as proceed back in time
   if(recdev_first<styr)
   {
     for (p=1;p<=pop;p++)
     for (g=1;g<=gmorph;g++)
     for (a=styr-recdev_first; a>=1; a--)
     {
       j=styr-a;
       natage(styr,p,g,a) *=mfexp(recdev(j)-biasadj(j)*half_sigmaRsq);
     }
   }

   //  note:  the above keeps SPB_pop_gp(styr) = SPB_equil.  It does not adjust for initial agecomp, but probably should
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 24 get_time_series */
FUNCTION void get_time_series()
  {
  dvar_vector pre_AL(1,nlength);
  dvariable crashtemp; dvariable crashtemp1;
  dvariable Z_adjuster2; dvariable Z_adjuster;
  if(Do_Morphcomp) Morphcomp_exp.initialize();

//loop over time
  SPB_current = SPB_yr(styr);  //  need these initial assignments in case recruitment distribution occurs before spawnbio&recruits
  if(recdev_doit(styr-1)>0)
  { Recruits = R1 * mfexp(recdev(styr-1)-biasadj(styr-1)*half_sigmaRsq); }
  else
  { Recruits = R1;}

  for (y=styr;y<=endyr;y++)
  {
    yz=y;
    if(STD_Yr_Reverse_F(y)>0) F_std(STD_Yr_Reverse_F(y))=0.0;
    t_base=styr+(y-styr)*nseas-1;
    if(y>styr)
    {
      if(time_vary_MG(y,0)>0 || save_for_report>0) get_MGsetup();
//      if(time_vary_MG(y,2)>0 || save_for_report>0) get_growth2();  // do growth before natmort in case lorensen M is used
      if(time_vary_MG(y,2)>0)
      {
        get_growth2();  // do growth before natmort in case lorensen M is used
      }
      if(time_vary_MG(y,1)>0) get_natmort();
      if(time_vary_MG(y,3)>0) get_wtlen();
      if(time_vary_MG(y,4)>0) get_recr_distribution();
      if(time_vary_MG(y,5)>0) get_migration();
      if(save_for_report>0)
      {
        if(time_vary_MG(y,1)>0 || time_vary_MG(y,2)>0 || time_vary_MG(y,3)>0)
        {
          get_saveGparm();
        }
      }
      get_selectivity();
    }

//  ProgLabel_  calculate average recruitment distributions to use in benchmark and forecast (V3.24P; Apr 24, 2013)
  if(y>=Bmark_Yr(1) && y<=Bmark_Yr(2))
    {
      if(y==Bmark_Yr(1)) 
      {recr_dist_benchmark.initialize();}
      for(gp=1;gp<=gender*N_GP;gp++)
      for(s=1;s<=nseas;s++)
      for(p=1;p<=pop;p++)
      {
        recr_dist_benchmark(gp,s,p)+=recr_dist_calc(gp,s,p);
      }
      if(y==Bmark_Yr(2))
      {recr_dist_benchmark/=float(Bmark_Yr(2)-Bmark_Yr(1)+1.0);}
    }
  if(y>=Fcast_yr(1) && y<=Fcast_yr(2))
    {
      if(y==Fcast_yr(1)) 
      {recr_dist_forecast.initialize();}
      for(gp=1;gp<=gender*N_GP;gp++)
      for(s=1;s<=nseas;s++)
      for(p=1;p<=pop;p++)
      {
        recr_dist_forecast(gp,s,p)+=recr_dist_calc(gp,s,p);
      }
      if(y==Fcast_yr(2))
      {recr_dist_forecast/=float(Fcast_yr(2)-Fcast_yr(1)+1.0);}
    }

    for (s=1;s<=nseas;s++)
    {
      t = t_base+s;
      smid=s+nseas;
//      if (do_once==1) cout<<" year, seas "<<y<<" "<<s<<" "<<t<<endl;
      if (docheckup==1) echoinput<<endl<<"************************************"<<endl<<" year, seas "<<y<<" "<<s<<endl;
      if(time_vary_MG(y,2)>0 || time_vary_MG(y,3)>0 || save_for_report==1 || WTage_rd>0)         //  ProgLabel_102
      {
        Make_AgeLength_Key();        // indexed by season and gmorph
      }
      if(y>styr)    // because styr is done as part of initial conditions
      {
        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0)
        {Make_FishSelex();}
        if(do_once==1 || MG_active(2)>0 || MG_active(3)>0)
        {
          Save_Wt_Age(t)=Wt_Age_beg(s);
        }
      }


      if(s==1)  //  calc some Smry_Table quantities that could be needed for exploitation rate calculations, but recalc these in the time_series report section
      {
        Smry_Table(y,2)=0.0;
        Smry_Table(y,3)=0.0;
        for(g=1;g<=gmorph;g++)
        if(use_morph(g)>0)
        {
        for(p=1;p<=pop;p++)
        {
          Smry_Table(y,2)+=natage(t,p,g)(Smry_Age,nages)*Save_Wt_Age(t,g)(Smry_Age,nages);
          Smry_Table(y,3)+=sum(natage(t,p,g)(Smry_Age,nages));   //sums to accumulate across submorphs and birthseasons
        }
        }
      }

// Compute spawning biomass
      if(s==spawn_seas && spawn_time_seas<0.0001)    //  compute spawning biomass if spawning at beginning of season so recruits could occur later this season
      {
        SPB_pop_gp(y).initialize();
        SPB_B_yr(y).initialize();
        SPB_N_yr(y).initialize();
        for (p=1;p<=pop;p++)
        {
          for (g=1;g<=gmorph;g++)
          if(sx(g)==1 && use_morph(g)>0)     //  female
          {
            SPB_pop_gp(y,p,GP4(g)) += fec(g)*natage(t,p,g);   // accumulates SSB by area and by growthpattern
            SPB_B_yr(y) += make_mature_bio(GP4(g))*natage(t,p,g);
            SPB_N_yr(y) += make_mature_numbers(GP4(g))*natage(t,p,g);
          }
        }
        SPB_current=sum(SPB_pop_gp(y));
        SPB_yr(y)=SPB_current;

        if(Hermaphro_Option>0)  // get male biomass
        {
          temp=0.0;
          MaleSPB(y).initialize();
          for (p=1;p<=pop;p++)
          {
            for (g=1;g<=gmorph;g++)
            if(sx(g)==2 && use_morph(g)>0)     //  male; all assumed to be mature
            {
              MaleSPB(y,p,GP4(g)) += Save_Wt_Age(t,g)*natage(t,p,g);   // accumulates SSB by area and by growthpattern
              temp += sum(natage(t,p,g));
            }
          }
          if(Hermaphro_maleSPB==1) // add MaleSPB to female SSB
          {
            SPB_current+=sum(MaleSPB(y));
            SPB_yr(y)=SPB_current;
            SPB_B_yr(y) += sum(MaleSPB(y));
            SPB_N_yr(y) += temp;
          }
        }

        Recruits=Spawn_Recr(SPB_current);  // calls to function Spawn_Recr
      }

      else
      {
        //  spawning biomass and total recruits will be calculated later so they can use Z
      }

// distribute Recruitment of age 0 fish among the pops and gmorphs
//  uses the most recent calculation of total recruits
      for (g=1;g<=gmorph;g++)
      if(use_morph(g)>0 && s==Bseas(g))
      {
        for (p=1;p<=pop;p++)
        { natage(t,p,g,0) = Recruits*recr_dist_calc(GP(g),s,p)*submorphdist(GP2(g));}
      }
      for (p=1;p<=pop;p++)     //loop over populations      ProgLabel_104
      {

        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0)
        {
          Nmid(g) = elem_prod(natage(t,p,g),surv1(s,GP3(g)));      //  get numbers-at-age(g,a) surviving to middle of time period
          if(docheckup==1) echoinput<<p<<" "<<g<<" "<<GP3(g)<<" area & morph "<<endl<<"N-at-age "<<natage(t,p,g)(0,min(6,nages))<<endl
           <<"survival "<<surv1(s,GP3(g))(0,min(6,nages))<<endl;
          if(save_for_report==1)
          {
            Save_PopLen(t,p,g)=0.0;
            Save_PopLen(t,p+pop,g)=0.0;  // later put midseason here
            Save_PopWt(t,p,g)=0.0;
            Save_PopWt(t,p+pop,g)=0.0;  // later put midseason here
            Save_PopAge(t,p,g)=0.0;
            Save_PopAge(t,p+pop,g)=0.0;  // later put midseason here
            for (a=0;a<=nages;a++)
            {
              Save_PopLen(t,p,g)+=value(natage(t,p,g,a))*value(ALK(s,g,a));
              Save_PopWt(t,p,g)+=value(natage(t,p,g,a))*value(elem_prod(ALK(s,g,a),wt_len(s,sx(g))));
              Save_PopAge(t,p,g,a)=value(natage(t,p,g,a));
            } // close age loop
          }
        }

        catage_tot.initialize();

        if(catch_seas_area(t,p,0)==1 && fishery_on_off==1)
        {
          switch (F_Method_use)
          {
            case 1:          // F_Method is Pope's approximation
            {
              for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       ProgLabel_105
              if (catch_seas_area(t,p,f)==1)
              {
                dvar_matrix catage_w=catage(t,f);      // do shallow copy

                vbio.initialize();
                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
                  gg=sx(g);              //  ProgLabel_107
                  // use sel_l to get total catch and use sel_l_r to get retained vbio
                  // note that vbio in numbers can be used for both survey abund and fishery available "biomass"
                  // vbio is for retained catch only;  harvest rate = retainedcatch/vbio;
                  // then harvestrate*catage_w = total kill by this fishery for this morph

                  catage_w(g)=elem_prod(Nmid(g),deadfish(s,g,f));     // total kill numbers at age
                  if(docheckup==1) echoinput<<"killrate "<<deadfish(s,g,f)(0,min(6,nages))<<endl;
                  if(catchunits(f)==1)
                  { vbio+=Nmid(g)*sel_al_2(s,g,f);}    // retained catch bio
                  else
                  { vbio+=Nmid(g)*sel_al_4(s,g,f);}  // retained catch numbers

                }  //close gmorph loop
                if(docheckup==1) echoinput<<"fleet vbio obscatch "<<f<<" "<<vbio<<" "<<catch_ret_obs(f,t)<<endl;
                crashtemp1=0.;
                crashtemp=max_harvest_rate-catch_ret_obs(f,t)/(vbio+NilNumbers);
                crashtemp1=posfun(crashtemp,0.000001,CrashPen);
                harvest_rate=max_harvest_rate-crashtemp1;
                if(crashtemp<0.&&rundetail>=2) {cout<<y<<" "<<f<<" crash vbio "<<catch_ret_obs(f,t)/vbio<<" "<<crashtemp<<
                 " "<<crashtemp1<<" "<<CrashPen<<" "<<harvest_rate<<endl;}
                Hrate(f,t) = harvest_rate;

                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
                  catage_w(g) *= harvest_rate;                       // Total catch at age
                  catage_tot(g) += catage_w(g);                      //catch at age for all fleets
                  gg=sx(g);
                  catch_fleet(t,f,2)+=Hrate(f,t)*Nmid(g)*deadfish_B(s,g,f);      // total fishery kill in biomass
                  catch_fleet(t,f,5)+=Hrate(f,t)*Nmid(g)*deadfish(s,g,f);     // total fishery kill in numbers
                  catch_fleet(t,f,1)+=Hrate(f,t)*Nmid(g)*sel_al_1(s,g,f);      //  total fishery encounter in biomass
                  catch_fleet(t,f,3)+=Hrate(f,t)*Nmid(g)*sel_al_2(s,g,f);      // retained fishery kill in biomass
                  catch_fleet(t,f,4)+=Hrate(f,t)*Nmid(g)*sel_al_3(s,g,f);      // encountered numbers
                  catch_fleet(t,f,6)+=Hrate(f,t)*Nmid(g)*sel_al_4(s,g,f);      // retained fishery kill in numbers
                }  // end g loop
              }  // close fishery

              for (g=1;g<=gmorph;g++)
              if(use_morph(g)>0)
              {
                for (a=0;a<=nages;a++)    //  check for negative abundance, starting at age 1
                {
                  crashtemp=max_harvest_rate-catage_tot(g,a)/(Nmid(g,a)+0.0000001);
                  crashtemp1=posfun(crashtemp,0.000001,CrashPen);
                  if(crashtemp<0.&&rundetail>=2) {cout<<" crash age "<<catage_tot(g,a)/(Nmid(g,a)+0.0000001)<<" "<<crashtemp<<
                    " "<<crashtemp1<<" "<<CrashPen<<" "<<(max_harvest_rate-crashtemp1)*Nmid(g,a)<<endl; }
                  if(crashtemp<0.&&docheckup==1) {echoinput<<" crash age "<<catage_tot(g,a)/(Nmid(g,a)+0.0000001)<<" "<<crashtemp<<
                    " "<<crashtemp1<<" "<<CrashPen<<" "<<(max_harvest_rate-crashtemp1)*Nmid(g,a)<<endl; }
                  catage_tot(g,a)=(max_harvest_rate-crashtemp1)*Nmid(g,a);
                }
                if(docheckup==1) echoinput<<"total catch-at-age for morph "<<g<<" "<<catage_tot(g)(0,min(6,nages))<<endl;
              }
              break;
            }   //  end Pope's approx

            case 2:          // continuous F_method
            {
              for (g=1;g<=gmorph;g++)
              if(use_morph(g)>0)
              {
                Z_rate(t,p,g)=natM(s,GP3(g));
                for (f=1;f<=Nfleet;f++)       //loop over fishing fleets to get Z       ProgLabel_105
                if (catch_seas_area(t,p,f)==1)
                {
                  Z_rate(t,p,g)+=deadfish(s,g,f)*Hrate(f,t);
                }
                Zrate2(p,g)=elem_div( (1.-mfexp(-seasdur(s)*Z_rate(t,p,g))), Z_rate(t,p,g));
              }
              for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       ProgLabel_105
              if (catch_seas_area(t,p,f)==1)
              {
                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
                  catch_fleet(t,f,1)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_1(s,g,f))*Zrate2(p,g);
                  catch_fleet(t,f,2)+=Hrate(f,t)*elem_prod(natage(t,p,g),deadfish_B(s,g,f))*Zrate2(p,g);
                  catch_fleet(t,f,3)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_2(s,g,f))*Zrate2(p,g); // retained bio
                  catch_fleet(t,f,4)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_3(s,g,f))*Zrate2(p,g);
                  catch_fleet(t,f,5)+=Hrate(f,t)*elem_prod(natage(t,p,g),deadfish(s,g,f))*Zrate2(p,g);
                  catch_fleet(t,f,6)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_4(s,g,f))*Zrate2(p,g);  // retained numbers
                  catage(t,f,g)=Hrate(f,t)*elem_prod(elem_prod(natage(t,p,g),deadfish(s,g,f)),Zrate2(p,g));
                }  //close gmorph loop
              }  // close fishery
              break;
            }   //  end continuous F method

            case 3:          // hybrid F_method
            {
              // start by doing a Pope's approximation
              for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       TAG_105
              if (catch_seas_area(t,p,f)==1 && catch_se_rd(f)>0.0)  // do exact catch for this fleet
              {
                vbio.initialize();
                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
                  if(catchunits(f)==1)
                    {vbio+=Nmid(g)*sel_al_2(s,g,f);}    // retained catch bio
                  else
                    {vbio+=Nmid(g)*sel_al_4(s,g,f);}  // retained catch numbers
                }  //close gmorph loop
                temp = catch_ret_obs(f,t)/(vbio+0.1*catch_ret_obs(f,t));  //  Pope's rate  robust
                join1=1./(1.+mfexp(30.*(temp-0.95)));  // steep logistic joiner at harvest rate of 0.95
                temp1=join1*temp + (1.-join1)*0.95;
                Hrate(f,t)=-log(1.-temp1)/seasdur(s);  // initial estimate of F (even though labelled as Hrate)
                //  done with starting values from Pope's approximation
              }
              else
              {
                // Hrate(f,t) previously set to zero or set to a parameter value
              }
            	for(int tune_F=1;tune_F<=F_Tune;tune_F++)
              {
                // add F+M to get Z using both exact and parameter based F
                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
                  Z_rate(t,p,g)=natM(s,GP3(g));
                  for (f=1;f<=Nfleet;f++)       //loop over fishing fleets to get Z
                  if (catch_seas_area(t,p,f)!=0)
                  {
                    Z_rate(t,p,g)+=deadfish(s,g,f)*Hrate(f,t);
                  }
                  Zrate2(p,g)=elem_div( (1.-mfexp(-seasdur(s)*Z_rate(t,p,g))), Z_rate(t,p,g));
                }
                //  now calc adjustment to Z based on changes to be made to Hrate
                if(tune_F<F_Tune)
                {
                  Z_adjuster2=0.0;   // this is the expected total catch that would occur with the current Hrates and Z
                  for (f=1;f<=Nfleet;f++)
                  if (catch_seas_area(t,p,f)==1 && catch_se_rd(f)>0.0)
                  {
                    for (g=1;g<=gmorph;g++)
                    if(use_morph(g)>0)
                    {
                      if(catchunits(f)==1)
                      {
                        Z_adjuster2+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_2(s,g,f))*Zrate2(p,g);  // biomass basis
                      }
                      else
                      {
                        Z_adjuster2+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_4(s,g,f))*Zrate2(p,g);  //  numbers basis
                      }
                    }  //close gmorph loop
                  }  // close fishery
                  Z_adjuster = totcatch_byarea(t,p)/(Z_adjuster2+0.0001);   // but this totcatch_by_area needs to exclude fisheries with F from param
                  for (g=1;g<=gmorph;g++)
                  if(use_morph(g)>0)
                  {
                    Z_rate(t,p,g)=natM(s,GP3(g)) + Z_adjuster*(Z_rate(t,p,g)-natM(s,GP3(g)));  // need to modify to only do the exact catches
                    Zrate2(p,g)=elem_div( (1.-mfexp(-seasdur(s)*Z_rate(t,p,g))), Z_rate(t,p,g));
                  }
                  for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       TAG_105
                  if(catch_seas_area(t,p,f)==1 && catch_se_rd(f)>0.0)
                  {
                    Z_adjuster2=0.;  // now use this to calc the available biomass to each fishery with the adjusted Zrate2
                    for (g=1;g<=gmorph;g++)
                    if(use_morph(g)>0)
                    {
                      if(catchunits(f)==1)
                      {
                        Z_adjuster2+=elem_prod(natage(t,p,g),sel_al_2(s,g,f))*Zrate2(p,g);
                      }
                      else
                      {
                        Z_adjuster2+=elem_prod(natage(t,p,g),sel_al_4(s,g,f))*Zrate2(p,g);
                      }
                    }  //close gmorph loop
                    temp=catch_ret_obs(f,t)/(Z_adjuster2+0.0001);
                    join1=1./(1.+mfexp(30.*(temp-0.95*max_harvest_rate)));
                    Hrate(f,t)=join1*temp + (1.-join1)*max_harvest_rate;
                  }  // close fishery
                }
                else
                {
                for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       ProgLabel_105
                if (catch_seas_area(t,p,f)==1)
                {
                  for (g=1;g<=gmorph;g++)
                  if(use_morph(g)>0)
                  {
                    catch_fleet(t,f,1)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_1(s,g,f))*Zrate2(p,g);
                    catch_fleet(t,f,2)+=Hrate(f,t)*elem_prod(natage(t,p,g),deadfish_B(s,g,f))*Zrate2(p,g);
                    catch_fleet(t,f,3)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_2(s,g,f))*Zrate2(p,g);
                    catch_fleet(t,f,4)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_3(s,g,f))*Zrate2(p,g);
                    catch_fleet(t,f,5)+=Hrate(f,t)*elem_prod(natage(t,p,g),deadfish(s,g,f))*Zrate2(p,g);
                    catch_fleet(t,f,6)+=Hrate(f,t)*elem_prod(natage(t,p,g),sel_al_4(s,g,f))*Zrate2(p,g);
                    catage(t,f,g)=Hrate(f,t)*elem_prod(elem_prod(natage(t,p,g),deadfish(s,g,f)),Zrate2(p,g));
                  }  //close gmorph loop
                }  // close fishery
              }
              }
              break;
            }   //  end hybrid F_Method
          }  // end F_Method switch
        }  //  end have some catch in this seas x area
        else
        {
          for (g=1;g<=gmorph;g++)
          if(use_morph(g)>0)
          {Z_rate(t,p,g)=natM(s,GP3(g));}
        }

// Compute spawning biomass if occurs after start of current season
      if(s==spawn_seas && spawn_time_seas>=0.0001)    //  compute spawning biomass
      {
        SPB_pop_gp(y).initialize();
        SPB_B_yr(y).initialize();
        SPB_N_yr(y).initialize();
        for (p=1;p<=pop;p++)
        {
          for (g=1;g<=gmorph;g++)
          if(sx(g)==1 && use_morph(g)>0)     //  female
          {
            SPB_pop_gp(y,p,GP4(g)) += fec(g)*elem_prod(natage(t,p,g),mfexp(-Z_rate(t,p,g)*spawn_time_seas));   // accumulates SSB by area and by growthpattern
            SPB_B_yr(y) += make_mature_bio(GP4(g))*natage(t,p,g);
            SPB_N_yr(y) += make_mature_numbers(GP4(g))*natage(t,p,g);
          }
        }
        SPB_current=sum(SPB_pop_gp(y));
        SPB_yr(y)=SPB_current;

        if(Hermaphro_Option>0)  // get male biomass
        {
          MaleSPB(y).initialize();
          temp=0.0;
          for (p=1;p<=pop;p++)
          {
            for (g=1;g<=gmorph;g++)
            if(sx(g)==2 && use_morph(g)>0)     //  male; all assumed to be mature
            {
              MaleSPB(y,p,GP4(g)) += Save_Wt_Age(t,g)*elem_prod(natage(t,p,g),mfexp(-Z_rate(t,p,g)*spawn_time_seas));   // accumulates SSB by area and by growthpattern
              temp += sum(natage(t,p,g));
            }
          }
          if(Hermaphro_maleSPB==1) // add MaleSPB to female SSB
          {
            SPB_current+=sum(MaleSPB(y));
            SPB_yr(y)=SPB_current;
            SPB_B_yr(y)+=sum(MaleSPB(y));
            SPB_N_yr(y)+=temp;
          }
        }
        Recruits=Spawn_Recr(SPB_current);  // calls to function Spawn_Recr
        //  these will get distributed among morphs at beginning of next season
      }

//  check for data
        for (f=1;f<=Ntypes;f++)       //loop over fleets and surveys with data       ProgLabel_105
        if (fleet_area(f)==p)
        {
          if(have_data(f,t)>0)      // ProgLabel_106  make the age-length available pop summed over all growth morphs
          {
            if(docheckup==1) echoinput<<" have data for fleet "<<f<<" in area: "<<p<<endl;
            exp_AL.initialize();

            timing = seasdur(s)*surveytime(f);
            for (g=1;g<=gmorph;g++)
            if(use_morph(g)>0)
            {
              gg=sx(g);
              if(gg==2)
              { L1=nlength1; L2= nlength2; A2=nages+1;}    //  move over on length dimension to store males
              else
              { L1=1; L2=nlength; A2=0;}
              if(F_Method==1) //  Pope's approximation
              {
                if(timing<0.0) // fishery catch
                {
                  for (a=0;a<=nages;a++)
                  {
                    temp=Nmid(g,a)*sel_a(y,f,gg,a);      //  get numbers-at-age(g,a) surviving to middle of time period
                    if(dolen(f)==1)
                    {pre_AL.shift(1)=temp*elem_prod(ALK(smid,g,a),sel_l(y,f,gg));}
                    else
                    {pre_AL.shift(1)=temp*ALK(smid,g,a);}
                    exp_AL(a+A2)(L1,L2) += pre_AL.shift(L1);  // shifted to store males in right place and accumulated across morphs
                  } // close age loop
                }
                else  // use explicit timing
                {
                  for (a=0;a<=nages;a++)
                  {
                    temp=(natage(t,p,g,a)*mfexp(-natM(s,GP3(g),a)*timing)-surveytime(f)*catage_tot(g,a))*sel_a(y,f,gg,a);
                    if(dolen(f)==1)
                    {pre_AL.shift(1)=temp*elem_prod(ALK(smid,g,a),sel_l(y,f,gg));}
                    else
                    {pre_AL.shift(1)=temp*ALK(smid,g,a);}
                    exp_AL(a+A2)(L1,L2) += pre_AL.shift(L1);  // shifted to store males in right place and accumulated across morphs
                  } // close age loop
                }
              }  // end Pope's approx

              else   //  continuous F
              {
                if(timing<0.0) // mimic fishery catch, but without Hrate so get available numbers
                {
                  for (a=0;a<=nages;a++)
                  {
                    temp=natage(t,p,g,a)*Zrate2(p,g,a)*sel_a(y,f,gg,a);
                    if(dolen(f)==1)
                    {pre_AL.shift(1)=temp*elem_prod(ALK(smid,g,a),sel_l(y,f,gg));}
                    else
                    {pre_AL.shift(1)=temp*ALK(smid,g,a);}
                    exp_AL(a+A2)(L1,L2) += pre_AL.shift(L1);  // shifted to store males in right place and accumulated across morphs
                  } // close age loop
                }
                else  // use explicit timing
                {
                  for (a=0;a<=nages;a++)
                  {
                    temp=natage(t,p,g,a)*mfexp(-Z_rate(t,p,g,a)*timing) *sel_a(y,f,gg,a);
                    if(dolen(f)==1)
                    {pre_AL.shift(1)=temp*elem_prod(ALK(smid,g,a),sel_l(y,f,gg));}
                    else
                    {pre_AL.shift(1)=temp*ALK(smid,g,a);}
                    exp_AL(a+A2)(L1,L2) += pre_AL.shift(L1);  // shifted to store males in right place and accumulated across morphs
                  } // close age loop
                }
              }  // end continuous F

              if(Do_Morphcomp)
              {
                if(Morphcomp_havedata(f,t,0)>0)
                {
                  Morphcomp_exp(Morphcomp_havedata(f,t,0),5+GP4(g))+=sum(exp_AL);     // total catch of this gpattern in this season x area
                }
              }
            } //close gmorph loop

        exp_l_temp=colsum(exp_AL);
        if(seltype(f,2)!=0)
          {exp_l_temp_ret=elem_prod(exp_l_temp,retain(y,f));}
        else
          {exp_l_temp_ret=exp_l_temp;}

//          end creation of selected A-L
// expected discard amount ************************************* ProgLabel_109
        if(yr_disc2(f,t)>0)
          {
            j=yr_disc2(f,t);
            if(catch_ret_obs(f,t)>0.0)
            {
              if(disc_units(f)==3)  // numbers regardless of catchunits for retained catch
              {
                exp_disc(f,j)=catch_fleet(t,f,4)-catch_fleet(t,f,6);
              }
              else if(catchunits(f)==1)  // biomass units for retained and discarded catch
              {
                exp_disc(f,j)=catch_fleet(t,f,1)-catch_fleet(t,f,3);  // discard in biomass
                if(disc_units(f)==2) exp_disc(f,j) /= (catch_fleet(t,f,1) + 0.0000001);
              }
              else   // numbers for retained and discarded catch
              {
                exp_disc(f,j)=catch_fleet(t,f,4)-catch_fleet(t,f,6);   // discard in numbers
                if(disc_units(f)==2) exp_disc(f,j) /= (catch_fleet(t,f,4) + 0.0000001);
              }
            }
            else
            {
              exp_disc(f,j)=-1.;
            }
          }

// expected mean weight *************************************** ProgLabel_109a
        if(nobs_mnwt>0)
        {
          j=yr_mnwt2(f,t,0);  //   sample from total catch
          if(j!=0)
          {
            if(j>0) 
            {exp_mnwt(j) = (exp_l_temp*wt_len2(s)) / sum(exp_l_temp);}
            else
            {exp_mnwt(-j) = (exp_l_temp*len_bins_m2) / sum(exp_l_temp);}
          }
          j=yr_mnwt2(f,t,1);  // sample from discard
          if(j!=0)
          {
            if(j>0) 
            {exp_mnwt(j) = (exp_l_temp-exp_l_temp_ret)*wt_len2(s) / (sum(exp_l_temp)-sum(exp_l_temp_ret));}
            else
            {exp_mnwt(-j) = (exp_l_temp-exp_l_temp_ret)*len_bins_m2 / (sum(exp_l_temp)-sum(exp_l_temp_ret));}
          }
          j=yr_mnwt2(f,t,2);  // sample from retained catch
          if(j!=0)
          {
            if(j>0) 
            {exp_mnwt(j) = (exp_l_temp_ret*wt_len2(s)) / sum(exp_l_temp_ret);}
            else
            {exp_mnwt(-j) = (exp_l_temp_ret*len_bins_m2) / sum(exp_l_temp_ret);}
          }
        }

// expected weight frequency ****************************************ProgLabel_041
    if(SzFreq_Nmeth>0)       //  have some sizefreq data
    {
//     create the transition matrices to convert population length bins to weight freq
      for(SzFreqMethod=1;SzFreqMethod<=SzFreq_Nmeth;SzFreqMethod++)
      {
        SzFreqMethod_seas=nseas*(SzFreqMethod-1)+s;     // index that combines sizefreqmethod and season and used in SzFreqTrans
        if(SzFreq_HaveObs2(SzFreqMethod,t)==f)  // first occurrence of this method at this time is with fleet = f
        {
          if(do_once==1 || (MG_active(3)>0 && (time_vary_MG(y,3)>0 )))  // calc the matrix because it may have changed
          {
            for(gg=1;gg<=gender;gg++)
            {
              if(gg==1)
              {z1=1;z2=nlength;ibin=0; ibinsave=0;}  // female
              else
              {z1=nlength1; z2=nlength2; ibin=0; ibinsave=SzFreq_Nbins(SzFreqMethod);}   // male
              topbin=0.;
              botbin=0.;

              switch(SzFreq_units(SzFreqMethod))    // biomass vs. numbers
              {
                case(1):  // units are biomass, so accumulate body weight into the bins;  Assume that bin demarcations are also in biomass
                {
                  if(SzFreq_Omit_Small(SzFreqMethod)==1)
                  {while(wt_len_low(s,z1+1)<SzFreq_bins(SzFreqMethod,1)) {z1++;}}      // ignore tiny fish

                  if( wt_len_low(s,nlength2) < SzFreq_bins(SzFreqMethod,SzFreq_Nbins(SzFreqMethod)))
                  {
                    N_warn++; cout<<" EXIT - see warning "<<endl;
                    warning<<" error:  max population size "<<wt_len_low(s,nlength2)<<" is less than max data bin "<<
                    SzFreq_bins(SzFreqMethod,SzFreq_Nbins(SzFreqMethod))<<
                    " for SzFreqMethod "<<SzFreqMethod<<endl;
                    exit(1);
                  }

                  for(z=z1;z<=z2;z++)
                  {
                    if(ibin==SzFreq_Nbins(SzFreqMethod))
                    {
                      SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=wt_len2(s,z);
                    }
                    else
                    {
                      if(wt_len_low(s,z)>=topbin)
                      {
                        ibin++; ibinsave++;
                      }
                      if(ibin>1)  {botbin=SzFreq_bins2(SzFreqMethod,ibin);}
                      if(ibin==SzFreq_Nbins(SzFreqMethod))
                      {
                        SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=wt_len2(s,z);
                        topbin=99999.;
                      }
                      else
                      {
                        topbin=SzFreq_bins2(SzFreqMethod,ibin+1);
                        if(wt_len_low(s,z)>=botbin && wt_len_low(s,z+1)<=topbin )
                        {
                          SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=wt_len2(s,z);
                        }
                        else
                        {
                          temp=(wt_len_low(s,z+1)-topbin)/wt_len_fd(s,z);  // frac in pop bin above (data bin +1)
                          temp1=wt_len_low(s,z)+(1.-temp*0.5)*wt_len_fd(s,z);  // approx body wt for these fish
                          temp2=wt_len_low(s,z)+(1.-temp)*0.5*wt_len_fd(s,z);  // approx body wt for  fish below
                          SzFreqTrans(SzFreqMethod_seas,z,ibinsave+1)=temp*temp1;
                          SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=(1.-temp)*temp2;
                        }
                      }
                    }
                  }
                  if(SzFreq_scale(SzFreqMethod)==2 && gg==gender)  // convert to pounds
                  {
                    SzFreqTrans(SzFreqMethod_seas)/=0.4536;
                  }
                  break;
                }  //  end of units in biomass
                // NOTE: even though  the transition matrix is currently in units of biomass distribution, there is no need to
                // normalize to sum to 1.0 here because the normalization will occur after it gets used to create SzFreq_exp

                case(2):   // units are numbers
                {
                  if(SzFreq_scale(SzFreqMethod)<=2)   //  bin demarcations are in weight units (1=kg, 2=lbs), so uses wt_len to compare to bins
                  {
                    if(SzFreq_Omit_Small(SzFreqMethod)==1)
                    {while(wt_len_low(s,z1+1)<SzFreq_bins(SzFreqMethod,1)) {z1++;}}      // ignore tiny fish

                    if( wt_len_low(s,nlength2) < SzFreq_bins(SzFreqMethod,SzFreq_Nbins(SzFreqMethod)))
                    {
                      N_warn++; cout<<" EXIT - see warning "<<endl;
                      warning<<" error:  max population size "<<wt_len_low(s,nlength2)<<" is less than max data bin "<<
                      SzFreq_bins(SzFreqMethod,SzFreq_Nbins(SzFreqMethod))<<
                      " for SzFreqMethod "<<SzFreqMethod<<endl;
                      exit(1);
                    }

                    for(z=z1;z<=z2;z++)
                    {
                      if(ibin==SzFreq_Nbins(SzFreqMethod))
                      {SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.;}  //echoinput<<" got to last ibin, so put rest of popbins here"<<endl;
                      else
                      {
                        if(wt_len_low(s,z)>=topbin) {ibin++; ibinsave++;}  //echoinput<<" incr ibin "<<z<<" "<<ibin<<" "<<len_bins(z)<<" "<<len_bins_dat(ibin);
                        if(ibin>1)  {botbin=SzFreq_bins2(SzFreqMethod,ibin);}
                        if(ibin==SzFreq_Nbins(SzFreqMethod))  // echoinput<<" got to last ibin, so put rest of popbins here"<<endl;
                        {
                          SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.;
                          topbin=99999.;
                        }
                        else
                        {
                          topbin=SzFreq_bins2(SzFreqMethod,ibin+1);
                          if(wt_len_low(s,z)>=botbin && wt_len_low(s,z+1)<=topbin )  //echoinput<<" pop inside dat, put here"<<endl;
                          {SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.;}
                          else       // echoinput<<" overlap"<<endl;
                          {
                            SzFreqTrans(SzFreqMethod_seas,z,ibinsave+1)=(wt_len_low(s,z+1)-topbin)/wt_len_fd(s,z);
                            SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.-SzFreqTrans(SzFreqMethod_seas,z,ibinsave+1);
                          }
                        }
                      }
                    }
                  }

                  else       //  bin demarcations are in length unit (3=cm, 4=inch) so uses population len_bins to compare to data bins
                  {
                    if(SzFreq_Omit_Small(SzFreqMethod)==1)
                    {while(len_bins2(z1+1)<SzFreq_bins(SzFreqMethod,1)) {z1++;}}      // ignore tiny fish
                    for(z=z1;z<=z2;z++)
                    {
                      if(ibin==SzFreq_Nbins(SzFreqMethod))
                      {SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.;} //echoinput<<" got to last ibin, so put rest of popbins here"<<endl;
                      else
                      {
                        if(len_bins2(z)>=topbin) {ibin++; ibinsave++;}  //echoinput<<" incr ibin "<<z<<" "<<ibin<<" "<<len_bins(z)<<" "<<len_bins_dat(ibin);
                        if(ibin>1)  {botbin=SzFreq_bins2(SzFreqMethod,ibin);}
                        if(ibin==SzFreq_Nbins(SzFreqMethod))  // echoinput<<" got to last ibin, so put rest of popbins here"<<endl;
                        {
                          SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.;
                          topbin=99999.;
                        }
                        else
                        {
                          topbin=SzFreq_bins2(SzFreqMethod,ibin+1);
                          if(len_bins2(z)>=botbin && len_bins2(z+1)<=topbin )  //echoinput<<" pop inside dat, put here"<<endl;
                          {SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.;}
                          else       // echoinput<<" overlap"<<endl;
                          {
                            SzFreqTrans(SzFreqMethod_seas,z,ibinsave+1)=(len_bins2(z+1)-topbin)/(len_bins2(z+1)-len_bins2(z));
                            SzFreqTrans(SzFreqMethod_seas,z,ibinsave)=1.-SzFreqTrans(SzFreqMethod_seas,z,ibinsave+1);
                          }
                        }
                      }
                    }
                  }
                  break;
                }  //  end of units in numbers
              }
              if(docheckup==1 && gg==gender) echoinput<<" sizefreq trans_matrix: method/season "<<SzFreqMethod<<" / "<<s<<endl
              <<trans(SzFreqTrans(SzFreqMethod_seas))<<endl;
            }  // end gender loop
          }  //  end needing to calc the matrix because it may have changed
        }  // done calculating the SzFreqTransition matrix for this method

        if(SzFreq_HaveObs(f,SzFreqMethod,t,1)>0)
        {
          for(iobs=SzFreq_HaveObs(f,SzFreqMethod,t,1);iobs<=SzFreq_HaveObs(f,SzFreqMethod,t,2);iobs++)
          {
            switch(SzFreq_obs_hdr(iobs,5))   // discard/retained partition
            {
              case(0):
              {
                SzFreq_exp(iobs)=trans(SzFreqTrans(SzFreqMethod_seas))*exp_l_temp;
                break;
              }
              case(1):
              {
                SzFreq_exp(iobs)=trans(SzFreqTrans(SzFreqMethod_seas))*(exp_l_temp-exp_l_temp_ret);
                break;
              }
              case(2):
              {
                SzFreq_exp(iobs)=trans(SzFreqTrans(SzFreqMethod_seas))*exp_l_temp_ret;
                break;
              }
            }
            if(gender==2)
            {
              k=SzFreq_obs_hdr(iobs,8);  // max bins for this method
              switch(SzFreq_obs_hdr(iobs,4))   //  combine, select or each gender
              {
                case(0):                    // combine genders
                {
                  for(ibin=1;ibin<=k;ibin++) SzFreq_exp(iobs,ibin)+=SzFreq_exp(iobs,k+ibin);
                  SzFreq_exp(iobs)(k+1,2*k)=0.0;
                  SzFreq_exp(iobs)(1,k)/=sum(SzFreq_exp(iobs)(1,k));
                  if(SzFreq_mincomp(SzFreqMethod)>0.0)
                  {
                    SzFreq_exp(iobs)(1,k)+=SzFreq_mincomp(SzFreqMethod);
                    SzFreq_exp(iobs)(1,k)/=sum(SzFreq_exp(iobs)(1,k));
                  }
                  break;
                }
                case(1):     // female only
                {
                  SzFreq_exp(iobs)(k+1,2*k)=0.0;  //  zero out the males so will not interfere with data generation
                  SzFreq_exp(iobs)(1,k)/=sum(SzFreq_exp(iobs)(1,k));
                  if(SzFreq_mincomp(SzFreqMethod)>0.0)
                  {
                    SzFreq_exp(iobs)(1,k)+=SzFreq_mincomp(SzFreqMethod);
                    SzFreq_exp(iobs)(1,k)/=sum(SzFreq_exp(iobs)(1,k));
                  }
                  break;
                }
                case(2):            //   male only
                {
                  ibin=SzFreq_obs_hdr(iobs,7);
                  SzFreq_exp(iobs)(1,ibin-1)=0.0;  //  zero out the females so will not interfere with data generation
                  SzFreq_exp(iobs)(ibin,k)/=sum(SzFreq_exp(iobs)(ibin,k));
                  if(SzFreq_mincomp(SzFreqMethod)>0.0)
                  {
                    SzFreq_exp(iobs)(ibin,k)+=SzFreq_mincomp(SzFreqMethod);
                    SzFreq_exp(iobs)(ibin,k)/=sum(SzFreq_exp(iobs)(ibin,k));
                  }
                  break;
                }
                case(3):           //  each gender
                {
                  SzFreq_exp(iobs)/=sum(SzFreq_exp(iobs));
                  if(SzFreq_mincomp(SzFreqMethod)>0.0)
                  {
                    SzFreq_exp(iobs)+=SzFreq_mincomp(SzFreqMethod);
                    SzFreq_exp(iobs)/=sum(SzFreq_exp(iobs));
                  }
                  break;
                }
              }  //  end gender switch
            }  // end have 2 genders
            else
            {
              k=SzFreq_obs_hdr(iobs,8);  // max bins for this method
              SzFreq_exp(iobs)(1,k)/=sum(SzFreq_exp(iobs)(1,k));
              if(SzFreq_mincomp(SzFreqMethod)>0.0)
              {
                SzFreq_exp(iobs)(1,k)+=SzFreq_mincomp(SzFreqMethod);
                SzFreq_exp(iobs)(1,k)/=sum(SzFreq_exp(iobs)(1,k));
              }
            }
          }  // end loop of obs for fleet = f
        }   //  end having some obs for this method in this fleet
      }  // end loop of sizefreqmethods
    }    //  end use of wt freq data


        // expected biomass index ************************************* ProgLabel_110
// NOTE that the Q scaler is factored in later on
       if(yr_cr2(f,t)>0)
       {
         j=yr_cr2(f,t);
         if (seltype(f,1)>=30)
         {
           if(docheckup==1) echoinput<<" special selex used "<<seltype(f,1)<<endl;
           switch(seltype(f,1))
           {
             case 30:  // spawning biomass  #30
             {
               if(pop==1 || fleet_area(f)==0)
               {
                 vbio=SPB_current;
                }
                else
                {
                  vbio=sum(SPB_pop_gp(y,fleet_area(f)));
                }
               break;
             }
             case 31:  // recruitment deviation  #31
             {
              if(y>=recdev_start && y<=recdev_end) 
              {vbio=mfexp(recdev(y));}
              else
              {vbio=1.0;}
              break;}
             case 32:  // recruitment without density-dependence (for pre-recruit survey) #32
             {
              if(y>=recdev_start && y<=recdev_end) 
              {vbio=SPB_current*mfexp(recdev(y));}
              else
              {vbio=SPB_current;}
              break;}
             case 33:  // recruitment  #33
             {vbio=Recruits; break;}

             case 34:  // spawning biomass depletion
             {
               if(pop==1 || fleet_area(f)==0)
               {
                 vbio=(SPB_current+1.0e-06)/(SPB_virgin+1.0e-06);
                }
                else
                {
                  vbio=(sum(SPB_pop_gp(y,fleet_area(f)))+1.0e-06)/(SPB_virgin+1.0e-06);
                }
               break;
             }
           }
         }
         else
         {
           if(cr_units(f)==1)
           {
             if(WTage_rd==1)  //  using empirical wt-at-age;  note that this cannot use GP specific bodyweights
             {
               if(seltype(f,2)>=1)
               {
                 agetemp = exp_AL * retain(y,f);    // retained only
               }
               else
               {
                 agetemp=rowsum(exp_AL);
               }
               vbio=0.0;
               for(a=0;a<=nages;a++) vbio+=WTage_emp(y,1,f,a)*agetemp(a);
               if(gender==2)
               {
                 for(a=0;a<=nages;a++) vbio+=WTage_emp(y,2,f,a)*agetemp(a+nages+1);
               }
             }
             else
             {vbio=exp_l_temp_ret*wt_len2(s);}   // biomass
           }
           else if(cr_units(f)==0)
             {vbio=sum(exp_l_temp_ret);}              // numbers
           else if(cr_units(f)==2)
           {vbio=Hrate(f,t);}   //  F rate
         }
         vul_bio(f,j)=value(vbio);
         if(Q_setup(f,1)>0) vbio=pow(vbio,1.0+Q_parm(Q_setup_parms(f,1)));
         if(cr_errtype(f)>=0)  //  lognormal
         {exp_cr(f,j)=log(vbio+0.000001);}
         else
         {exp_cr(f,j)=vbio;}
       }

//  Calc expected size comp ************************************* ProgLabel_120

        {
         for(i=1;i<=Nobs_l(f);i++)                          // loop all obs of this type
         if(yr_l_t(f,i)==t)
         {
          if(LenBin_option>1)
          {
          if(mkt_l(f,i)==0) {exp_l(f,i) = make_len_bin*exp_l_temp;}           // expected size comp  MAtrix * vector = vector
          else if(mkt_l(f,i)==1) {exp_l(f,i) = make_len_bin*(exp_l_temp-exp_l_temp_ret);}  // discard sample
          else {exp_l(f,i) = make_len_bin*exp_l_temp_ret;}    // retained only
          }
          else  //  using data_bins same as pop_bins
          {
          if(mkt_l(f,i)==0) {exp_l(f,i) = exp_l_temp;}           // expected size comp  MAtrix * vector = vector
          else if(mkt_l(f,i)==1) {exp_l(f,i) = (exp_l_temp-exp_l_temp_ret);}  // discard sample
          else {exp_l(f,i) = exp_l_temp_ret;}    // retained only
          }
          if(docheckup==1) echoinput<<" len obs "<<mkt_l(f,i)<<" "<<tails_l(f,i)<<endl<<obs_l(f,i)<<endl<<exp_l(f,i)<<endl;
         //  code for tail compression, etc in the likelihood section to allow for superyear combinations                                                                      // mkt=0 Do nothing
         }  // end lengthcomp loop
        }  // end have lengthcomps

//  ProgLabel_130   size@age in the bins and with ageing error
       if(nobs_ms_tot>0)
       {
         for(i=1;i<=Nobs_ms(f);i++)         // loop obs in this type/time  must be contiguous in list
         if(yr_ms_t(f,i)==t)
         {
           k=abs(ageerr_type_ms(f,i));                           //  age-err type  where the sign selects length vs. weight
           if(ageerr_type_ms(f,i)>0)  // values are length at age
           {
             if(mkt_ms(f,i)==0)
             {
               exp_a_temp = age_age(k) * rowsum(exp_AL);             //  numbers at binned age = age_age(bins,age) * sum(age)
               exp_ms(f,i) = age_age(k) * (exp_AL * len_bins_m2);  // numbers * length
               exp_ms_sq(f,i) = age_age(k) * (exp_AL * len_bins_sq);  // numbers * length^2
             }
             if(mkt_ms(f,i)==1)
             {
               exp_a_temp = age_age(k) * (exp_AL * (1-retain(y,f)));             //  numbers at binned age = age_age(bins,age) * sum(age)
               exp_ms(f,i) = age_age(k) * (exp_AL * elem_prod((1-retain(y,f)),len_bins_m2));  // numbers * length
               exp_ms_sq(f,i) = age_age(k) * (exp_AL * elem_prod((1-retain(y,f)),len_bins_sq));  // numbers * length^2
             }
             if(mkt_ms(f,i)==2)
             {
               exp_a_temp = age_age(k) * (exp_AL * retain(y,f) );             //  numbers at binned age = age_age(bins,age) * sum(age)
               exp_ms(f,i) = age_age(k) * (exp_AL * elem_prod((retain(y,f)),len_bins_m2));  // numbers * length
               exp_ms_sq(f,i) = age_age(k) * (exp_AL * elem_prod((retain(y,f)),len_bins_sq));  // numbers * length^2
             }
           }
           else  // values are weight at age
           {
             if(mkt_ms(f,i)==0)
             {
               exp_a_temp = age_age(k) * rowsum(exp_AL);             //  numbers at binned age = age_age(bins,age) * sum(age)
               exp_ms(f,i) = age_age(k) * (exp_AL * wt_len2(s));  // numbers * bodywt
               exp_ms_sq(f,i) = age_age(k) * (exp_AL * wt_len2_sq(s));  // numbers * bodywt^2
             }
             if(mkt_ms(f,i)==1)
             {
               exp_a_temp = age_age(k) * (exp_AL * (1-retain(y,f)));             //  numbers at binned age = age_age(bins,age) * sum(age)
               exp_ms(f,i) = age_age(k) * (exp_AL * elem_prod((1-retain(y,f)),wt_len2(s)));  // numbers * bodywt
               exp_ms_sq(f,i) = age_age(k) * (exp_AL * elem_prod((1-retain(y,f)),wt_len2_sq(s)));  // numbers * bodywt^2
             }
             if(mkt_ms(f,i)==2)
             {
               exp_a_temp = age_age(k) * (exp_AL * retain(y,f) );             //  numbers at binned age = age_age(bins,age) * sum(age)
               exp_ms(f,i) = age_age(k) * (exp_AL * elem_prod((retain(y,f)),wt_len2(s)));  // numbers * bodywt
               exp_ms_sq(f,i) = age_age(k) * (exp_AL * elem_prod((retain(y,f)),wt_len2_sq(s)));  // numbers * bodywt^2
             }
           }
           exp_ms(f,i)+=1.0e-6;
           exp_a_temp+=1.0e-6;
           exp_ms_sq(f,i)+=1.0e-6;
           exp_ms_sq(f,i) = sqrt(
                                 elem_div(
                                          (exp_ms_sq(f,i) - elem_div(elem_prod(exp_ms(f,i),exp_ms(f,i)), exp_a_temp)),
                                          exp_a_temp
                                         )
                                )
                                + 0.000001;    //std.err. of size at binned age = sqrt( (P2-P1*P1/P0) / P0 )
           exp_ms(f,i) = elem_div(exp_ms(f,i), exp_a_temp);   //  mean size at binned age
         }
       }   // end size-at-age

//                              if (do_once==1) cout<<" got LAA  "<<endl;

//  get expected age comp in the bins and with ageing error
        {
         for(i=1;i<=Nobs_a(f);i++)                          // loop all obs of  this type
         if(yr_a_t(f,i)==t)
         {
          k=ageerr_type_a(f,i);                           //  age-err type

          if(use_Lbin_filter(f,i)==0)
           {                                              // sum across all length bins
          if(mkt_a(f,i)==0) agetemp = rowsum(exp_AL);             //  numbers at binned age = age_age(bins,age) * sum(age)
          if(mkt_a(f,i)==1) agetemp = exp_AL * (1.-retain(y,f));  // discard sample
          if(mkt_a(f,i)==2) agetemp = exp_AL * retain(y,f);    // retained only
           }
          else
           {            // only use ages from specified range of size bins
                        // Lbin_filter is a vector with 0 for unselected size bins and 1 for selected bins
          if(mkt_a(f,i)==0) agetemp = exp_AL * Lbin_filter(f,i);             //  numbers at binned age = age_age(bins,age) * sum(age)
          if(mkt_a(f,i)==1) agetemp = exp_AL * elem_prod(Lbin_filter(f,i),(1.-retain(y,f)));  // discard sample
          if(mkt_a(f,i)==2) agetemp = exp_AL * elem_prod(Lbin_filter(f,i),retain(y,f));    // retained only
           }
          exp_a(f,i) = age_age(k) * agetemp;

          if(docheckup==1) echoinput<<" age obs2 "<<mkt_a(f,i)<<endl<<obs_a(f,i)<<endl<<exp_a(f,i)<<endl;

         }  // end agecomp loop within fleet/time
        }  // end have agecomps in this fleet/time
          }  // end have data
        }  //  close fleet loop

          if(docheckup==1) echoinput<<" now do survival "<<endl;

//  survival to next year
        if(s==nseas) {k=1;} else {k=0;}   //      advance age or not
        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0)
        {
          if(F_Method==1)  // pope's
          {
              if(s<nseas) natage(t+1,p,g,0) = natage(t,p,g,0)*surv2(s,GP3(g),0) -catage_tot(g,0)*surv1(s,GP3(g),0);  // advance age zero within year

              for (a=1;a<nages;a++)
              {natage(t+1,p,g,a) = natage(t,p,g,a-k)*surv2(s,GP3(g),a-k) -catage_tot(g,a-k)*surv1(s,GP3(g),a-k);}
              natage(t+1,p,g,nages) = natage(t,p,g,nages)*surv2(s,GP3(g),nages) - catage_tot(g,nages)*surv1(s,GP3(g),nages);   // plus group
              if(s==nseas) natage(t+1,p,g,nages) += natage(t,p,g,nages-1)*surv2(s,GP3(g),nages-1) - catage_tot(g,nages-1)*surv1(s,GP3(g),nages-1);
              if(save_for_report==1)
              {
                j=p+pop;
                for (a=0;a<=nages;a++)
                {
                  Save_PopLen(t,j,g)+=value(Nmid(g,a)-0.5*catage_tot(g,a))*value(ALK(smid,g,a));
                  Save_PopWt(t,j,g)+= value(Nmid(g,a)-0.5*catage_tot(g,a))*value(elem_prod(ALK(smid,g,a),wt_len(s,sx(g))));
                  Save_PopAge(t,j,g,a)=value(Nmid(g,a)-0.5*catage_tot(g,a));
                } // close age loop
              }
          }
          else   // continuous F
          {
            if(s<nseas) natage(t+1,p,g,0) = natage(t,p,g,0)*mfexp(-Z_rate(t,p,g,0)*seasdur(s));  // advance age zero within year
            for (a=1;a<nages;a++) {natage(t+1,p,g,a) = natage(t,p,g,a-k)*mfexp(-Z_rate(t,p,g,a-k)*seasdur(s));}
            natage(t+1,p,g,nages) = natage(t,p,g,nages)*mfexp(-Z_rate(t,p,g,nages)*seasdur(s));   // plus group
            if(s==nseas) natage(t+1,p,g,nages) += natage(t,p,g,nages-1)*mfexp(-Z_rate(t,p,g,nages-1)*seasdur(s));
              if(save_for_report==1)
              {
                j=p+pop;
                for (a=0;a<=nages;a++)
                {
                  Save_PopLen(t,j,g)+=value(natage(t,p,g,a)*mfexp(-Z_rate(t,p,g,a)*0.5*seasdur(s)))*value(ALK(smid,g,a));
                  Save_PopWt(t,j,g)+= value(natage(t,p,g,a)*mfexp(-Z_rate(t,p,g,a)*0.5*seasdur(s)))*value(elem_prod(ALK(smid,g,a),wt_len(s,sx(g))));
                  Save_PopAge(t,j,g,a)=value(natage(t,p,g,a)*mfexp(-Z_rate(t,p,g,a)*0.5*seasdur(s)));
                } // close age loop
              }
          }
          if(docheckup==1)
          {
            echoinput<<"Morph "<<g<<endl<<" natM: "<<natM(s,GP3(g))(0,min(6,nages))<<endl;
            if(F_Method>1) echoinput<<"Zed: "<<Z_rate(t,p,g)(0,min(6,nages))<<endl;
            echoinput<<"Nsurv: "<<natage(t+1,p,g)(0,min(6,nages))<<endl;
          }

        } // close gmorph loop
      } //close area loop

          if(docheckup==1) echoinput<<" finished survival, do hermophroditism and migration "<<endl;

//hermaphroditism
      if(Hermaphro_Option>0)
      {
        if(Hermaphro_seas==-1 || Hermaphro_seas==s)
        {
          k=gmorph/2;
          for (p=1;p<=pop;p++)  //   area
          for (g=1;g<=k;g++)  //  loop females
          if(use_morph(g)>0)
          {
            for (a=1;a<nages;a++)
            {
              natage(t+1,p,g+k,a) += natage(t+1,p,g,a)*Hermaphro_val(GP4(g),a-1); // increment males
              natage(t+1,p,g,a) *= (1.-Hermaphro_val(GP4(g),a-1)); // decrement females
            }
          }
        }
      }

//do migration between populations, for each gmorph and age
      if(do_migration>0)  // movement between areas in time series
      {
        natage_temp=natage(t+1);
        natage(t+1)=0.0;
        for (p=1;p<=pop;p++)  //   source population
        for (p2=1;p2<=pop;p2++)  //  destination population
        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0)
        {
          k=move_pattern(s,GP4(g),p,p2);
          if(k>0) natage(t+1,p2,g) += elem_prod(natage_temp(p,g),migrrate(y,k));}
      }  //  end migration

//  save selectivity*Hrate for tag-recapture
      if(Do_TG>0 && t>=TG_timestart)
      {
        for(g=1;g<=gmorph;g++)
        for(f=1;f<=Nfleet;f++)
        {
          Sel_for_tag(t,g,f) = sel_al_4(s,g,f)*Hrate(f,t);
        }
      }

      if(docheckup==1) echoinput<<" save some annual catch and F quantities "<<endl;

      if( (save_for_report>0) || ((sd_phase() || mceval_phase()) && (initial_params::mc_phase==0)|| (F_ballpark_yr>=styr)) )
      {
        for (f=1;f<=Nfleet;f++)
        {
          for (k=1;k<=6;k++)
          {
            annual_catch(y,k)+=catch_fleet(t,f,k);
          }
          if(F_Method==1)
          {
            annual_F(y,1)+=Hrate(f,t);
          }
          else
          {
            annual_F(y,1)+=Hrate(f,t)*seasdur(s);
          }
        }
        
          if(docheckup==1) echoinput<<" calculate the F reporting statistic "<<endl;

        if(s==nseas)
        {
  //  sum across p and g the number of survivors to end of the year
  //  also project from the initial numbers and M, the number of survivors without F
  //  then F = ln(n+1/n)(M+F) - ln(n+1/n)(M only), but ln(n) cancels out, so only need the ln of the ratio of the two ending quantities
          temp=0.0;
          temp1=0.0;
          temp2=0.0;
          for(g=1;g<=gmorph;g++)
          if(use_morph(g)>0)
          {
            for(p=1;p<=pop;p++)
            {
              for(a=F_reporting_ages(1);a<=F_reporting_ages(2);a++)   //  should not let a go higher than nages-2 because of accumulator
              {
                temp1+=natage(t+1,p,g,a+1);
                if(nseas==1)
                {
                  temp2+=natage(t,p,g,a)*mfexp(-seasdur(s)*natM(s,GP3(g),a));
                }
                else
                {
                  temp3=natage(t-nseas+1,p,g,a);  //  numbers at begin of year
                  for(j=1;j<=nseas;j++) {temp3*=mfexp(-seasdur(j)*natM(j,GP3(g),a));}
                  temp2+=temp3;
                }
              }
            }
          }
          annual_F(y,2) = log(temp2)-log(temp1);
        
          if(docheckup==1) echoinput<<" annual F "<<annual_F(y,2)<<endl;
        
          if(docheckup==1) echoinput<<" do special STD quantities "<<endl;

          if(STD_Yr_Reverse_F(y)>0)  //  save selected std quantity
          {
            if(F_reporting<=1)
            {
              F_std(STD_Yr_Reverse_F(y))=annual_catch(y,2)/Smry_Table(y,2);  // dead catch biomass/summary biomass
                                                                             //  does not exactly correspond to F, which is for total catch
            }
            else if(F_reporting==2)
            {
              F_std(STD_Yr_Reverse_F(y))=annual_catch(y,5)/Smry_Table(y,3);  // dead catch numbers/summary numbers
            }
            else if(F_reporting==3)
            {
              F_std(STD_Yr_Reverse_F(y))=annual_F(y,1);
            }
            else if(F_reporting==4)
            {
              F_std(STD_Yr_Reverse_F(y))=annual_F(y,2);
            }
          }
        }  //  end s==nseas
      }

    } //close season loop


//  calculate the SPR value for
    if( (save_for_report>0) || ((sd_phase() || mceval_phase()) && (initial_params::mc_phase==0)) )
    {
      eq_yr=y; equ_Recr=Recr_virgin; bio_yr=y;
      Fishon=0;
      dvariable SPR_unf;
      Do_Equil_Calc();                      //  call function to do equilibrium calculation with current year's biology
      SPR_unf=SPB_equil;
      Smry_Table(y,11)=SPR_unf;
      Smry_Table(y,13)=GenTime;
      Fishon=1;
      Do_Equil_Calc();                      //  call function to do equilibrium calculation with current year's biology and F
      if(STD_Yr_Reverse_Ofish(y)>0)
      {
        SPR_std(STD_Yr_Reverse_Ofish(y))=SPB_equil/SPR_unf;
      }
      Smry_Table(y,9)=(totbio);
      Smry_Table(y,10)=(smrybio);
      Smry_Table(y,12)=(SPB_equil);
      Smry_Table(y,14)=(YPR_dead);
      for(g=1;g<=gmorph;g++)
      {
        Smry_Table(y,20+g)=(cumF(g));
        Smry_Table(y,20+gmorph+g)=(maxF(g));
      }
    }
          if(docheckup==1) echoinput<<" finished year "<<y<<endl;

  } //close year loop

// ******************************************************************
//  do tag mortality, movement and recapture
  if(Do_TG>0)
  {
    dvariable TG_init_loss;
    dvariable TG_chron_loss;
    TG_recap_exp.initialize();
    for(TG=1;TG<=N_TG;TG++)
    {
      firstseas=int(TG_release(TG,4));  // release season
      t=int(TG_release(TG,5));  // release t index calculated in data section from year and season of release
      p=int(TG_release(TG,2));  // release area
      gg=int(TG_release(TG,6));  // gender (1=fem; 2=male; 0=both
      a1=int(TG_release(TG,7));  // age at release

      TG_alive.initialize();
      if(gg==0)
      {
        for (g=1;g<=gmorph;g++) {TG_alive(p,g) = natage(t,p,g,a1);}   //  gets both genders
      }
      else
      {
        for (g=1;g<=gmorph;g++)
        {
          if(sx(g)==gg) {TG_alive(p,g) = natage(t,p,g,a1);}  //  only does the selected gender
        }
      }
      TG_init_loss = mfexp(TG_parm(TG))/(1.+mfexp(TG_parm(TG)));
      TG_chron_loss = mfexp(TG_parm(TG+N_TG))/(1.+mfexp(TG_parm(TG+N_TG)));
      for(f=1;f<=Nfleet;f++)
      {
        k=3*N_TG+f;
        TG_report(f) = mfexp(TG_parm(k))/(1.0+mfexp(TG_parm(k)));
        TG_rep_decay(f) = TG_parm(k+Nfleet);
      }
      TG_alive /= sum(TG_alive);     // proportions across morphs at age a1 in release area p at time of release t
      TG_alive *= TG_release(TG,8);   //  number released as distributed across morphs
      TG_alive *= (1.-TG_init_loss);  // initial mortality
      if(save_for_report>0)
      {
        TG_save(TG,1)=value(TG_init_loss);
        TG_save(TG,2)=value(TG_chron_loss);
      }

      TG_t=0;
      for(y=TG_release(TG,3);y<=endyr;y++)
      {
        for(s=firstseas;s<=nseas;s++)
        {
          if(save_for_report>0 && TG_t<=TG_endtime(TG))
          {TG_save(TG,3+TG_t)=value(sum(TG_alive));} //  OK to do simple sum because only selected morphs are populated

          for(p=1;p<=pop;p++)
          {
            for (g=1;g<=gmorph;g++)
            if(TG_use_morph(TG,g)>0)
            {
              for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       ProgLabel_105
              if (fleet_area(f)==p)
              {
// calculate recaptures by fleet
// NOTE:  Sel_for_tag(t,g,f,a1) = sel_al_4(s,g,f,a1)*Hrate(f,t)
                if(F_Method==1)
                {
                  TG_recap_exp(TG,TG_t,f)+=TG_alive(p,g)  // tags recaptured
                  *mfexp(-(natM(s,GP3(g),a1)+TG_chron_loss)*seasdur_half(s))
                  *Sel_for_tag(t,g,f,a1)
                  *TG_report(f)
                  *mfexp(TG_t*TG_rep_decay(f));
                }
                else   // use for method 2 and 3
                {
                  TG_recap_exp(TG,TG_t,f)+=TG_alive(p,g)
                  *Sel_for_tag(t,g,f,a1)/(Z_rate(t,p,g,a1)+TG_chron_loss)
                  *(1.-mfexp(-seasdur(s)*(Z_rate(t,p,g,a1)+TG_chron_loss)))
                  *TG_report(f)
                  *mfexp(TG_t*TG_rep_decay(f));
                }
  if(docheckup==1) echoinput<<" TG_"<<TG<<" y_"<<y<<" s_"<<s<<" area_"<<p<<" g_"<<g<<" GP3_"<<GP3(g)<<" f_"<<f<<" a1_"<<a1<<" Sel_"<<Sel_for_tag(t,g,f,a1)<<" TG_alive_"<<TG_alive(p,g)<<" TG_obs_"<<TG_recap_obs(TG,TG_t,f)<<" TG_exp_"<<TG_recap_exp(TG,TG_t,f)<<endl;
              }  // end fleet loop for recaptures

// calculate survival
              if(F_Method==1)
              {
                if(s==nseas) {k=1;} else {k=0;}   //      advance age or not
                if((a1+k)<=(nages-1)) {f=a1-1;} else {f=nages-2;}
                TG_alive(p,g)*=(natage(t+1,p,g,f+k)+1.0e-10)/(natage(t,p,g,f)+1.0e-10)*(1.-TG_chron_loss*seasdur(s));
              }
              else   //  use for method 2 and 3
              {
                TG_alive(p,g)*=mfexp(-seasdur(s)*(Z_rate(t,p,g,a1)+TG_chron_loss));
              }
            }  // end morph loop
          }  // end area loop

          if(Hermaphro_Option>0)
          {
            if(Hermaphro_seas==-1 || Hermaphro_seas==s)
            {
              k=gmorph/2;
              for (p=1;p<=pop;p++)  //   area
              for (g=1;g<=k;g++)  //  loop females
              if(use_morph(g)>0)
              {
                TG_alive(p,g+k) += TG_alive(p,g)*Hermaphro_val(GP4(g),a1); // increment males
                TG_alive(p,g) *= (1.-Hermaphro_val(GP4(g),a1)); // decrement females
              }
            }
          }

          if(do_migration>0)  //  movement between areas of tags
          {
            TG_alive_temp=TG_alive;
            TG_alive=0.0;
            for (g=1;g<=gmorph;g++)
            if(use_morph(g)>0)
            {
              for (p=1;p<=pop;p++)  //   source population
              for (p2=1;p2<=pop;p2++)  //  destination population
              {
                k=move_pattern(s,GP4(g),p,p2);
                if(k>0) TG_alive(p2,g) += TG_alive_temp(p,g)*migrrate(y,k,a1);
              }
            }
            if(docheckup==1) echoinput<<" Tag_alive after survival and movement "<<endl<<TG_alive<<endl;
          }
          t++;         //  increment seasonal time counter
          if(TG_t<TG_endtime(TG)) TG_t++;
          if(s==nseas && a1<nages) a1++;
        }  // end seasons
        firstseas=1;  // so start with season 1 in year following the tag release
      }  // end years
    }  //  end loop of tag groups
  }  // end having tag groups
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 25 evaluate_the_objective_function */
FUNCTION void evaluate_the_objective_function()
  {
  surv_like.initialize();   Q_dev_like.initialize(); disc_like.initialize();   length_like.initialize(); age_like.initialize();
  sizeage_like.initialize(); parm_like.initialize(); parm_dev_like.initialize(); log_q_cr.initialize();
  mnwt_like.initialize(); equ_catch_like.initialize(); recr_like.initialize(); Fcast_recr_like.initialize();
  catch_like.initialize(); Morphcomp_like.initialize(); TG_like1.initialize(); TG_like2.initialize();
  obj_fun=0.0;

    int k_phase=current_phase();
    if(k_phase>max_lambda_phase) k_phase=max_lambda_phase;


// Fit to indices                       ProgLabel_210
    if(nobs_cr>0)
    {
      for(f=1;f<=Ntypes;f++)
      if(surv_lambda(f,k_phase)>0.0 || save_for_report>0)      // skip if zero emphasis
      {
        if(nyr_cr(f)>0)
        {
          se_cr_use(f) = se_cr_obs(f);
          if(Q_setup(f,3)>0) se_cr_use(f)+=Q_parm(Q_setup_parms(f,3));  // add extra stderr

          for(j=1;j<=N_suprper_cr(f);j++)                  // do super years
          {
            temp=0.0;
            for(i=suprper_cr1(f,j);i<=suprper_cr2(f,j);i++) {temp+=exp_cr(f,i)*suprper_cr_sampwt(f,i);} // combine across range of observations
//  sampwt sums to 1.0, so temp contains the weighted average
            for(i=suprper_cr1(f,j);i<=suprper_cr2(f,j);i++) {exp_cr(f,i)=temp;}   // assign average to each obs
          }
          if(Q_setup(f,4)==0 || Q_setup(f,4)==1  || Q_setup(f,4)==5 )   // do float option=automatic scaling      // ProgLabel_211
          {                                       //  NOTE:  cannot use float option if error type is normal
            temp=0.; temp1=0.; temp2=0.;
            for (i=1;i<=nyr_cr(f);i++)
            {
              if(yr_cr_use(f,i) > 0)
              {
                temp2 += (Ln_obs_cr(f,i)-exp_cr(f,i))/square(se_cr_use(f,i));
                temp += 1.0/square(se_cr_use(f,i));
                temp1 += 1.;
              }
            }
            if(Q_setup(f,4)==0)                               // mean q, with nobiasadjustment
            {log_q_cr(f) = temp2/temp;}
            else                  // for value = 1 or 5       // mean q with variance bias adjustment
            {log_q_cr(f) = (temp2 + temp1*0.5)/temp;}
            if(Q_setup(f,4)==5) Q_parm(Q_setup_parms(f,4))=log_q_cr(f,1);   // base Q
          }
          else if(Q_setup(f,4)<=-1)        // mirror Q from lower numbered survey
                                           // because Q is a vector for each observation, the mirror is to the first observation's Q
                                           // so time-varying property cannot be mirrored
          {log_q_cr(f) = log_q_cr(-Q_setup(f,4),1);}

          else                                               //  Q from parameter
          {
            log_q_cr(f) = Q_parm(Q_setup_parms(f,4));   // base Q
                                                        // note that Q is in log space even if normal error structure is used

//  trend or block effect on Q

//  seasonal or cyclic effect on Q

// environmental effect on Q
            if(Q_setup(f,2)>0)    // environ effect on log(q)  multiplicative
            {
              for(i=1;i<=nyr_cr(f);i++)
              {log_q_cr(f,i) += Q_parm(Q_setup_parms(f,2)) * env_data(yr_cr_y(f,i),Q_setup(f,2));}  // note that this environ effect is after the dev effect!
            }
            else if(Q_setup(f,2)<0)    // environ effect on log(q)  additive
            {
              for(i=1;i<=nyr_cr(f);i++)
              {log_q_cr(f,i) += Q_parm(Q_setup_parms(f,2)) + env_data(yr_cr_y(f,i),-Q_setup(f,2));}
            }

// random deviations or random walk
            if(Q_setup(f,4)==3 || Q_setup(f,4)==4 )
            {
              temp=0.0; temp2=0.0; temp1=0.;
              if(Q_setup(f,4)==3)  // random devs
              {
                for(i=1;i<=nyr_cr(f);i++)
                {
                  j=Q_setup_parms(f,4)+i;
                  log_q_cr(f,i)+=Q_parm(j);
                  temp+=Q_parm(j); temp2+=square(Q_parm(j)); temp1+=1.;
                }
              }
              else if(Q_setup(f,4)==4)   // random walk
              {
                for(i=2;i<=nyr_cr(f);i++)
                {
                  j=Q_setup_parms(f,4)+i-1;
                  log_q_cr(f,i)=log_q_cr(f,i-1)+Q_parm(j);
                  temp+=Q_parm(j); temp2+=square(Q_parm(j)); temp1+=1.;
                }
              }
              Q_dev_like(f,1)=square(1.+square(temp))-1.;  // not used for randwalk
              if(temp1>0.0) Q_dev_like(f,2)=sqrt((temp2+0.0000001)/temp1);  // this is calculated but not used because redundant with the prior penalty
            }
          }

          q_cr(f) = mfexp(log_q_cr(f));        // get q in arithmetic space

//  calculate the survey logL
          if(cr_errtype(f)==0)  // lognormal
          {
            for(i=1;i<=nyr_cr(f);i++)
            if(yr_cr_use(f,i)>0)
            {
              surv_like(f) +=0.5*square( ( Ln_obs_cr(f,i)-exp_cr(f,i)-log_q_cr(f,i) ) / se_cr_use(f,i)) + sd_offset*log(se_cr_use(f,i));
//            should add a term for 0.5*s^2 for bias adjustment so that parameter approach will be same as the  biasadjusted scaling approach
            }
          }
          else if (cr_errtype(f)>0)        // t-distribution
          {
            dvariable df = cr_errtype(f);
            for(i=1;i<=nyr_cr(f);i++)
            if(yr_cr_use(f,i)>0)
            {
              surv_like(f) +=((df+1.)/2.)*log((1.+square((Ln_obs_cr(f,i)-exp_cr(f,i)-log_q_cr(f,i) ))/(df*square(se_cr_use(f,i))) )) + sd_offset*log(se_cr_use(f,i));
            }
          }
          else if(cr_errtype(f)==-1)  // normal
          {
            for(i=1;i<=nyr_cr(f);i++)
            {
              if(yr_cr_use(f,i)>0)
              {
                surv_like(f) +=0.5*square( ( Ln_obs_cr(f,i)-exp_cr(f,i)*q_cr(f,i) ) / se_cr_use(f,i)) + sd_offset*log(se_cr_use(f,i));
              }
            }
          }

        }    // end having obs for this survey
      }
       if(do_once==1) cout<<" did survey obj_fun "<<surv_like<<endl;
    }

// Fit to discard biomass (normal)               ProgLabel_215
  if(nobs_disc>0)
  {
    for(f=1;f<=Nfleet;f++)
    if(disc_lambda(f,k_phase)>0.0 || save_for_report>0)
    {
      if(nyr_disc(f)>0)
      {
        for(j=1;j<=N_suprper_disc(f);j++)                  // do super years
        {
          temp=0.0;
          for(i=suprper_disc1(f,j);i<=suprper_disc2(f,j);i++) {temp+=exp_disc(f,i)*suprper_disc_sampwt(f,i);} // combine across range of observations
          for(i=suprper_disc1(f,j);i<=suprper_disc2(f,j);i++) {exp_disc(f,i)=temp;}   // assign back to each obs
        }

        if(disc_errtype(f)>=1)  // T -distribution
        {
          for (i=1;i<=nyr_disc(f);i++)
          if(yr_disc_use(f,i)>=0.0)
          {
            disc_like(f) +=0.5*(disc_errtype(f)+1.)*log((1.+square(obs_disc(f,i)-exp_disc(f,i))/(disc_errtype(f)*square(sd_disc(f,i))) )) + sd_offset*log(sd_disc(f,i));
          }
        }
        else if (disc_errtype(f)==0)  // normal error, with input CV
        {
          for (i=1;i<=nyr_disc(f);i++)
          if(yr_disc_use(f,i)>=0.0)
          {
            disc_like(f) +=0.5*square( (obs_disc(f,i)-exp_disc(f,i)) / sd_disc(f,i)) + sd_offset*log(sd_disc(f,i));
          }
        }
        else if (disc_errtype(f)==-1)  // normal error with input se
        {
          for (i=1;i<=nyr_disc(f);i++)
          if(yr_disc_use(f,i)>=0.0)
          {
            disc_like(f) +=0.5*square( (obs_disc(f,i)-exp_disc(f,i)) / sd_disc(f,i)) + sd_offset*log(sd_disc(f,i));
          }
        }
        else  // lognormal  where input cv_disc must contain se in log space
        {
          for (i=1;i<=nyr_disc(f);i++)
          if(yr_disc_use(f,i)>=0.0)
          {
            disc_like(f) +=0.5*square( log(obs_disc(f,i)/exp_disc(f,i)) / sd_disc(f,i)) + sd_offset*log(sd_disc(f,i));
          }
        }
      }
    }
    if(do_once==1) cout<<" did discard obj_fun "<<disc_like<<endl;
  }

// Fit to mean body wt (normal with student's T)
   if(nobs_mnwt>0)
   {
     for(i=1;i<=nobs_mnwt;i++)
     if(mnwtdata(5,i)>0. && mnwtdata(3,i)>0.)
     {
       mnwt_like(mnwtdata(3,i)) +=0.5*(DF_bodywt+1.)*log(1.+square(mnwtdata(5,i)-exp_mnwt(i))/mnwtdata(8,i))+ mnwtdata(9,i);
     }
      // mnwtdata(6,i)+=var_adjust(3,mnwtdata(3,i));  adjusted input error as a CV
      // if(mnwtdata(6,i)<=0.0) mnwtdata(6,i)=0.001;
      // mnwtdata(7,i)=mnwtdata(5,i)*mnwtdata(6,i);  se = obs*CV
      // mnwtdata(8,i)=DF_bodywt*square(mnwtdata(7,i));   error as T-dist
      // mnwtdata(9,i)=sd_offset*log(mnwtdata(7,i));  for the -log(s) component

     if(do_once==1) cout<<" did meanwt obj_fun "<<mnwt_like<<endl;
   }

//Length comp likelihood (multinomial)        // ProgLabel_220
   if(Nobs_l_tot>0)
   {
   for(f=1;f<=Ntypes;f++)
   if(length_lambda(f,k_phase)>0.0 || save_for_report>0)
    {
    if(Nobs_l(f)>=1)
    {

     length_like(f) = -offset_l(f);
     for(j=1;j<=N_suprper_l(f);j++)                  // do each super period
     {
       exp_l_temp_dat.initialize();
       for(i=suprper_l1(f,j);i<=suprper_l2(f,j);i++) 
       {
         exp_l_temp_dat+=exp_l(f,i)*suprper_l_sampwt(f,i);  // combine across range of observations
       }
//       exp_l_temp_dat/=sum(exp_l_temp_dat);   // normalize not needed because converted to proportions later
       for(i=suprper_l1(f,j);i<=suprper_l2(f,j);i++)
       {
         exp_l(f,i)=exp_l_temp_dat;   // assign back to all obs
       }
     }

     for(i=1;i<=Nobs_l(f);i++)
     {
     if(gender==2)
     {
       if(gen_l(f,i)==0)                        // combined sex observation   ProgLabel_124
       {
         for(z=1;z<=nlen_bin;z++) {exp_l(f,i,z)+=exp_l(f,i,z+nlen_bin);}
         exp_l(f,i)(nlen_binP,nlen_bin2)=0.00;
       }
       else if(gen_l(f,i)==1)                   // female only
       {exp_l(f,i)(nlen_binP,nlen_bin2)=0.00;}
       else if(gen_l(f,i)==2)                   // male only
       {exp_l(f,i)(1,nlen_bin)=0.00;}
       else if(gen_l(f,i)==3 && CombGender_l>0)
       {
         for(z=1;z<=CombGender_l;z++)
         {exp_l(f,i,z)+=exp_l(f,i,z+nlen_bin);  exp_l(f,i,z+nlen_bin)=0.00;}
       }
     }
     exp_l(f,i) /= sum(exp_l(f,i));
     tails_w=ivector(tails_l(f,i));


        if(gen_l(f,i)!=2)                          // ProgLabel_126
        {
         if(tails_w(1)>1)
           {exp_l(f,i,tails_w(1))=sum(exp_l(f,i)(1,tails_w(1)));  exp_l(f,i)(1,tails_w(1)-1)=0.;}
         if(tails_w(2)<nlen_bin)
           {exp_l(f,i,tails_w(2))=sum(exp_l(f,i)(tails_w(2),nlen_bin));  exp_l(f,i)(tails_w(2)+1,nlen_bin)=0.;}
         exp_l(f,i)(tails_w(1),tails_w(2))+= min_comp;
        }

        if(gender==2 && gen_l(f,i)>=2)
        {
         if(tails_w(3)>nlen_binP)
           {exp_l(f,i,tails_w(3))=sum(exp_l(f,i)(nlen_binP,tails_w(3)));  exp_l(f,i)(nlen_binP,tails_w(3)-1)=0.;}
         if(tails_w(4)<nlen_bin2)
           {exp_l(f,i,tails_w(4))=sum(exp_l(f,i)(tails_w(4),nlen_bin2));  exp_l(f,i)(tails_w(4)+1,nlen_bin2)=0.;}
         exp_l(f,i)(tails_w(3),tails_w(4))+= min_comp;
        }
        exp_l(f,i) /= sum(exp_l(f,i));

     if(nsamp_l(f,i)>0 && header_l(f,i,3)>0)
     {
       if(gen_l(f,i) !=2) length_like(f) -= nsamp_l(f,i) *
       obs_l(f,i)(tails_w(1),tails_w(2)) * log(exp_l(f,i)(tails_w(1),tails_w(2)));
       if(gen_l(f,i) >=2 && gender==2) length_like(f) -= nsamp_l(f,i) *
       obs_l(f,i)(tails_w(3),tails_w(4)) * log(exp_l(f,i)(tails_w(3),tails_w(4)));
     }
     }
    }
   }
  if(do_once==1) cout<<" did lencomp obj_fun "<<length_like<<endl;
   }

// age comp likelihood (multinomial)        // ProgLabel_230
  if(Nobs_a_tot>0)
  {
    for(f=1;f<=Ntypes;f++)
    if(age_lambda(f,k_phase)>0.0 || save_for_report>0)
    {
      if(Nobs_a(f)>=1)
      {
        age_like(f) = -offset_a(f);

        for(j=1;j<=N_suprper_a(f);j++)                  // do super years  Max of 20 allowed per type(f)
        {
          exp_a_temp.initialize();
          for(i=suprper_a1(f,j);i<=suprper_a2(f,j);i++) 
          {
            exp_a_temp+=exp_a(f,i)*suprper_a_sampwt(f,i);  // combine across range of observations
          }
//          exp_a_temp/=(1.0e-15+sum(exp_a_temp));                                        // normalize
          for(i=suprper_a1(f,j);i<=suprper_a2(f,j);i++) exp_a(f,i)=exp_a_temp;   // assign back to each original obs
        }

        for(i=1;i<=Nobs_a(f);i++)
        {
          if(gender==2)
          {
            if(gen_a(f,i)==0)                         // combined sex observation
            {
              for(z=1;z<=n_abins;z++) {exp_a(f,i,z)+=exp_a(f,i,z+n_abins);}
              exp_a(f,i)(n_abins1,n_abins2)=0.00;
            }
            else if(gen_a(f,i)==1)                   // female only
            {exp_a(f,i)(n_abins1,n_abins2)=0.00;}
            else if(gen_a(f,i)==2)                   // male only
            {exp_a(f,i)(1,n_abins)=0.00;}
            else if(gen_a(f,i)==3 && CombGender_a>0)
            {
              for(z=1;z<=CombGender_a;z++)
              {exp_a(f,i,z)+=exp_a(f,i,z+n_abins);  exp_a(f,i,z+n_abins)=0.00;}
            }
          }
          exp_a(f,i) /= (1.0e-15+sum(exp_a(f,i)));                      // proportion at binned age

          if(save_for_report==1)
          {
            exp_a_temp=exp_a(f,i);
            exp_meanage(f,i,2)=sum(elem_prod(exp_a_temp,age_bins_mean));  //  get mean age across both sexes before adding min_comp and compressing
            temp=exp_a_temp(1);
            temp1=0.0;
            if(gender==2) temp+=exp_a_temp(1+n_abins);
            z=1;
            while(temp<=0.05)
            {
              temp1=temp;
              z++;
              temp+=exp_a_temp(z);
              if(gender==2) temp+=exp_a_temp(z+n_abins);
            }
            if(z<n_abins) 
              {exp_meanage(f,i,1)=age_bins(z) + (0.05-temp1)/(temp-temp1);}
            else
              {exp_meanage(f,i,1)=age_bins_mean(z);}

            while(temp<=0.95)
            {
              temp1=temp;
              z++;
              temp+=exp_a_temp(z);
              if(gender==2) temp+=exp_a_temp(z+n_abins);
            }
            if(z<n_abins) 
              {exp_meanage(f,i,3)=age_bins(z) + (0.95-temp1)/(temp-temp1);}
            else
              {exp_meanage(f,i,3)=age_bins_mean(z);}
          }

          tails_w=ivector(tails_a(f,i));
          if(gen_a(f,i)!=2)
          {
            if(tails_w(1)>1)
            {exp_a(f,i,tails_w(1))=sum(exp_a(f,i)(1,tails_w(1)));  exp_a(f,i)(1,tails_w(1)-1)=0.;}
            if(tails_w(2)<n_abins)
            {exp_a(f,i,tails_w(2))=sum(exp_a(f,i)(tails_w(2),n_abins));  exp_a(f,i)(tails_w(2)+1,n_abins)=0.;}
            exp_a(f,i)(tails_w(1),tails_w(2))+= min_comp;
          }

          if(gender==2 && gen_a(f,i)>=2)
          {
            if(tails_w(3)>n_abins1)
            {exp_a(f,i,tails_w(3))=sum(exp_a(f,i)(n_abins1,tails_w(3)));  exp_a(f,i)(n_abins1,tails_w(3)-1)=0.;}
            if(tails_w(4)<n_abins2)
            {exp_a(f,i,tails_w(4))=sum(exp_a(f,i)(tails_w(4),n_abins2));  exp_a(f,i)(tails_w(4)+1,n_abins2)=0.;}
            exp_a(f,i)(tails_w(3),tails_w(4))+= min_comp;
          }

          exp_a(f,i) /= (1.0e-15+sum(exp_a(f,i)));

          if(nsamp_a(f,i)>0 && header_a(f,i,3)>0)
          {
            if(gen_a(f,i) !=2) age_like(f) -= nsamp_a(f,i) *
            obs_a(f,i)(tails_w(1),tails_w(2)) * log(exp_a(f,i)(tails_w(1),tails_w(2)));
            if(gen_a(f,i) >=2 && gender==2) age_like(f) -= nsamp_a(f,i) *
            obs_a(f,i)(tails_w(3),tails_w(4)) * log(exp_a(f,i)(tails_w(3),tails_w(4)));
          }
        }
      }
    }
    if(do_once==1) cout<<" did agecomp obj_fun "<<age_like<<endl;
  }

// size@age likelihood  (normal)      // ProgLabel_240
    if(nobs_ms_tot>0)
    {
      for(f=1;f<=Ntypes;f++)
      if(Nobs_ms(f)>0 && sizeage_lambda(f,k_phase)>0.0)
      {
         for(j=1;j<=N_suprper_ms(f);j++)
         {
           exp_a_temp.initialize();
           for(i=suprper_ms1(f,j);i<=suprper_ms2(f,j);i++) {exp_a_temp+=exp_ms(f,i)*suprper_ms_sampwt(f,i);} // combine across range of observations
           for(i=suprper_ms1(f,j);i<=suprper_ms2(f,j);i++) exp_ms(f,i)=exp_a_temp;   // assign back to each original obs
         }

         for(i=1;i<=Nobs_ms(f);i++)
         if(use_ms(f,i)>0)
         {
           for(b=1;b<=n_abins2;b++)
           {
             if(obs_ms_n(f,i,b)>0 && obs_ms(f,i,b)>0)
             {
               sizeage_like(f) += 0.5*square((obs_ms(f,i,b) -exp_ms(f,i,b)) / (exp_ms_sq(f,i,b)/obs_ms_n(f,i,b)))
               + sd_offset*log(exp_ms_sq(f,i,b)/obs_ms_n(f,i,b));
               //  where:        obs_ms_n(f,i,b)=sqrt(var_adjust(6,f)*obs_ms_n(f,i,b));
             }
           }
         }
       }
        if(do_once==1) cout<<" did meanlength obj_fun "<<sizeage_like<<endl;
     }


//  Size Frequency likelihood
    if(SzFreq_Nmeth>0)       //  have some sizefreq data
    {
      // create super-period expected values
      for(j=1;j<=N_suprper_SzFreq;j++)
      {
        a=suprper_SzFreq_start(j);  // get observation index
        SzFreq_exp(a)*=suprper_SzFreq_sampwt(a);  //  start creating weighted average
        for(iobs=a+1;iobs<=suprper_SzFreq_end(j);iobs++) {SzFreq_exp(a)+=SzFreq_exp(iobs)*suprper_SzFreq_sampwt(iobs);}  //  combine into the first obs of this superperiod
        for(iobs=a+1;iobs<=suprper_SzFreq_end(j);iobs++) {SzFreq_exp(iobs)=SzFreq_exp(a);}  //  assign back to all obs
      }

      SzFreq_like=-SzFreq_like_base;  // initializes
      for(iobs=1;iobs<=SzFreq_totobs;iobs++)
      {
        if(SzFreq_obs_hdr(iobs,3)>0)
        {
          k=SzFreq_obs_hdr(iobs,6);
          f=abs(SzFreq_obs_hdr(iobs,3));
          z1=SzFreq_obs_hdr(iobs,7);
          z2=SzFreq_obs_hdr(iobs,8);
          SzFreq_like(SzFreq_LikeComponent(f,k))-=SzFreq_sampleN(iobs)*SzFreq_obs(iobs)(z1,z2)*log(SzFreq_exp(iobs)(z1,z2));
        }
      }

      if(do_once==1) cout<<" did sizefreq obj_fun: "<<SzFreq_like<<endl;
    }

    if(Do_Morphcomp>0)
    {
      for(iobs=1;iobs<=Morphcomp_nobs;iobs++)
      {
        k=5+Morphcomp_nmorph;
        Morphcomp_exp(iobs)(6,k) /= sum(Morphcomp_exp(iobs)(6,k));
        Morphcomp_exp(iobs)(6,k) += Morphcomp_mincomp;
        Morphcomp_exp(iobs)(6,k) /= 1.+Morphcomp_mincomp*Morphcomp_nmorph;
        if(Morphcomp_obs(iobs,5)>0.) Morphcomp_like -= Morphcomp_obs(iobs,5)*Morphcomp_obs(iobs)(6,k) * log(elem_div(Morphcomp_exp(iobs)(6,k),Morphcomp_obs(iobs)(6,k)));
      }
    if(do_once==1) cout<<" did morphcomp obj_fun "<<Morphcomp_like<<endl;
    }

    if(Do_TG>0)
    {
      for(TG=1;TG<=N_TG;TG++)
      {
        overdisp=TG_parm(2*N_TG+TG);
        for(TG_t=TG_mixperiod;TG_t<=TG_endtime(TG);TG_t++)
        {
          TG_recap_exp(TG,TG_t)(1,Nfleet)+=1.0e-6;  // add a tiny amount
          TG_recap_exp(TG,TG_t,0) = sum(TG_recap_exp(TG,TG_t)(1,Nfleet));
          TG_recap_exp(TG,TG_t)(1,Nfleet)/=TG_recap_exp(TG,TG_t,0);
          if(Nfleet>1) TG_like1(TG)-=TG_recap_obs(TG,TG_t,0)* (TG_recap_obs(TG,TG_t)(1,Nfleet) * log(TG_recap_exp(TG,TG_t)(1,Nfleet)));
          TG_like2(TG)-=log_negbinomial_density(TG_recap_obs(TG,TG_t,0),TG_recap_exp(TG,TG_t,0),overdisp);
        }
      }
    if(do_once==1) cout<<" did tag obj_fun "<<TG_like1<<endl<<TG_like2<<endl;
    }

//  initial equilibrium catch  (lognormal)      // ProgLabel_250
    if(init_equ_lambda(k_phase)>0.0)
    {
      for(f=1;f<=Nfleet;f++)
      {
        if(catch_se_rd(f)>0.0 &&  obs_equ_catch(f)>0.0)
//          {equ_catch_like += 0.5*square( (log(obs_equ_catch(f)) -log(est_equ_catch(f)+0.000001)) / catch_se(styr-1,f));}
          {equ_catch_like += 0.5*square( (log(1.1*obs_equ_catch(f)) -log(est_equ_catch(f)+0.1*obs_equ_catch(f))) / catch_se(styr-1,f));}
      }
      if(do_once==1) cout<<" initequ_catch -log(L) "<<equ_catch_like<<endl;
    }

    if(F_Method>1)
    {
      for(f=1;f<=Nfleet;f++)
      {
        if(catchunits(f)==1)
        {i=3;}  //  biomass
        else
        {i=6;}  //  numbers
      
        for(t=styr;t<=TimeMax;t++)
        {
          if(catch_se_rd(f)>0.0 && catch_ret_obs(f,t)>0.0)
          {
//          catch_like(f) += 0.5*square( (log(catch_ret_obs(f,t)) -log(catch_fleet(t,f,i)+0.000001)) / catch_se(t,f));
            catch_like(f) += 0.5*square( (log(1.1*catch_ret_obs(f,t)) -log(catch_fleet(t,f,i)+0.1*catch_ret_obs(f,t))) / catch_se(t,f));
          }
        }
      }
    }

//The recruitment prior is assumed to be a lognormal pdf with expected
// value equal to the deterministic stock-recruitment curve          // ProgLabel_260
//  R1 deviation is weighted by ave_age because R1 represents a time series of recruitments
//  SR_parm(N_SRparm+1) is sigmaR
//  SR_parm(N_SRparm+4) is rho, the autocorrelation coefficient
//  POP code from Ianelli
//  if (rho>0)
//    for (i=styr_rec+1;i<=endyr;i++)
//      rec_like(1) += square((chi(i)- rho*chi(i-1)) /(sqrt(1.-rho*rho))) / (2.*sigmaRsq) + log(sigr);
//  else
//    rec_like(1)    = (norm2( chi +  sigmaRsq/2. ) / (2*sigmaRsq)) / (2.*sigmaRsq) + size_count(chi)*log(sigr);

    if(recrdev_lambda(k_phase)>0.0 && (do_recdev>0 || recdev_do_early>0) )
    {
      recr_like = sd_offset_rec*log(sigmaR);
      // where sd_offset_rec takes account for the number of recruitment years fully estimated
      // this is calculated as the sum of the biasadj vector
      if(SR_autocorr==0)
      {
      recr_like += norm2(recdev(recdev_first,recdev_end))/two_sigmaRsq;
      }
      else
      {
        rho=SR_parm(N_SRparm2);
        recr_like += square(recdev(recdev_first))/two_sigmaRsq;
        for(y=recdev_first+1;y<=recdev_end;y++)
        {
          recr_like += square(recdev(y)-rho*recdev(y-1)) / ((1.0-rho*rho)*two_sigmaRsq);
        }
      }
      recr_like += 0.5 * square( log(R1/R1_exp) / (sigmaR/ave_age) );
      if(do_once==1) cout<<" did recruitdev obj_fun "<<recr_like<<endl;
    }
    if(Do_Forecast>0)
    {
      if(recdev_end<endyr)
      {
        Fcast_recr_like = Fcast_recr_lambda*(norm2(Fcast_recruitments(recdev_end+1,endyr)))/two_sigmaRsq;
//        Fcast_recr_like += sd_offset_fore*log(sigmaR);  this is now part of the recr_liker logL calculated above
      }
      else
      {Fcast_recr_like=0.0;}
      Fcast_recr_like += (norm2(Fcast_recruitments(endyr+1,YrMax)))/two_sigmaRsq;  // ss3
      if(Do_Impl_Error>0) Fcast_recr_like+=(norm2(Fcast_impl_error(endyr+1,YrMax)))/(2.0*Impl_Error_Std*Impl_Error_Std);  // implementation error
    }

//  ____________________________________________________
//  likelihood for the parameter priors  // ProgLabel_270
    dvariable mu; dvariable tau; dvariable Aprior; dvariable Bprior;
    int Ptype;
    dvariable Pconst;
    Pconst=0.0001;

    if(parm_prior_lambda(k_phase)>0.0 || Do_all_priors>0)
    {
      for(i=1;i<=N_MGparm2;i++)
      if(MGparm_PRtype(i)>-1 && (active(MGparm(i))|| Do_all_priors>0))
        {
        MGparm_Like(i)=Get_Prior(MGparm_PRtype(i), MGparm_LO(i), MGparm_HI(i), MGparm_PR(i), MGparm_CV(i), MGparm(i));
        parm_like+=MGparm_Like(i);
        }

      for(i=1;i<=Nfleet;i++)
      if(init_F_PRtype(i)>-1 && (active(init_F(i))|| Do_all_priors>0))
        {
        init_F_Like(i)=Get_Prior(init_F_PRtype(i), init_F_LO(i), init_F_HI(i), init_F_PR(i), init_F_CV(i), init_F(i));
        parm_like+=init_F_Like(i);
        }

      for(i=1;i<=Q_Npar;i++)
      if(Q_parm_1(i,5)>-1 && (active(Q_parm(i))|| Do_all_priors>0))
        {
        Q_parm_Like(i)=Get_Prior(Q_parm_1(i,5), Q_parm_1(i,1), Q_parm_1(i,2), Q_parm_1(i,4), Q_parm_1(i,6), Q_parm(i));
        parm_like+=Q_parm_Like(i);
        }

      for(i=1;i<=N_selparm2;i++)
      if(selparm_PRtype(i)>-1 && (active(selparm(i))|| Do_all_priors>0))
        {
        selparm_Like(i)=Get_Prior(selparm_PRtype(i), selparm_LO(i), selparm_HI(i), selparm_PR(i), selparm_CV(i), selparm(i));
        parm_like+=selparm_Like(i);
        }

    if(Do_TG>0)
    {
      k=3*N_TG+2*Nfleet;
      for(i=1;i<=k;i++)
      if(TG_parm2(i,5)>-1 && (active(TG_parm(i))|| Do_all_priors>0))
      {
        TG_parm_Like(i)=Get_Prior(TG_parm2(i,5), TG_parm_LO(i), TG_parm_HI(i), TG_parm2(i,4), TG_parm2(i,6), TG_parm(i));
        parm_like+=TG_parm_Like(i);
      }
    }

    for(i=1;i<=N_SRparm2;i++)
      if(SR_parm_1(i,5)>-1 && (active(SR_parm(i))|| Do_all_priors>0))
        {
        SR_parm_Like(i)=Get_Prior(SR_parm_1(i,5), SR_parm_1(i,1), SR_parm_1(i,2), SR_parm_1(i,4), SR_parm_1(i,6), SR_parm(i));
        parm_like+=SR_parm_Like(i);
        }
    }

    if(recdev_cycle>0)
    {
      temp=0.0; temp1=0.0;
      for(i=1;i<=recdev_cycle;i++)
      {
        if(recdev_cycle_parm_RD(i,5)>-1 && (active(recdev_cycle_parm(i))|| Do_all_priors>0))
        {
          recdev_cycle_Like(i)=Get_Prior(recdev_cycle_parm_RD(i,5), recdev_cycle_parm_RD(i,1), recdev_cycle_parm_RD(i,2), recdev_cycle_parm_RD(i,4), recdev_cycle_parm_RD(i,6), recdev_cycle_parm(i));
          parm_like+=recdev_cycle_Like(i);
          temp+=mfexp(recdev_cycle_parm(i));
          temp1+=1.0;
        }
      }
      temp-=temp1;
      parm_like+=10000.*temp*temp;  //  similar to ADMB's approach to getting zero-centered dev_vectors
    }
//   likelihood for the random parameter time series  ProgLabel_280
    if(MGparm_dev_PH>0 && parm_dev_lambda(k_phase)>0.0 )
    {
     for(i=1;i<=N_MGparm_dev;i++)
     for(j=MGparm_dev_minyr(i);j<=MGparm_dev_maxyr(i);j++)
     {parm_dev_like += 0.5*square( MGparm_dev(i,j) / MGparm_dev_stddev(i) );}
    }

    for(f=1;f<=Ntypes;f++)
      if(Q_setup(f,4)==3)
      {
      parm_dev_like += Q_dev_like(f,1); // mean component for dev approach (var component is already in the parm priors)
                                        //  do not include for randwalk (Qsetup==4)
      }

    if(selparm_dev_PH>0 && parm_dev_lambda(k_phase)>0.0 )
    {
     for(i=1;i<=N_selparm_dev;i++)
     for(j=selparm_dev_minyr(i);j<=selparm_dev_maxyr(i);j++)
     {parm_dev_like += 0.5*square( selparm_dev(i,j) / selparm_dev_stddev(i) );}
    }

    if(F_ballpark_yr>=styr)
      {
        if(F_Method==1)
        {temp=annual_F(F_ballpark_yr,1);}  //  sum of H rates
        else
        {temp=annual_F(F_ballpark_yr,2);}  //  actual total F for specified age range
//  in future, could allow specification of a range of years for averaging the F statistic        
        F_ballpark_like = 0.5*square( log((F_ballpark+1.0e-6)/(temp+1.0e-6)) / 1.0) + sd_offset*log(1.0);
      }
    else
      {F_ballpark_like=0.0;}

//  *******************************************
//  Penalty for soft boundaries, uses the symmetric beta prior code
  if(SoftBound>0)
  {
    SoftBoundPen=0.0;

      for(i=1;i<=N_selparm2;i++)
      if(active(selparm(i)))
        {SoftBoundPen+=Get_Prior(1, selparm_LO(i), selparm_HI(i), 1., 0.001, selparm(i));}
  }


//   CrashPen = square(1.0+CrashPen)-1.0;   this was used until V3.00L  7/10/2008
     CrashPen = square(1.0+ (1000.*CrashPen/(1000.+CrashPen)))-1.0;
//  sum the likelihood components weighted by lambda factors // ProgLabel_290
//   cout<<" obj_fun start "<<obj_fun<<endl;
   obj_fun = column(surv_lambda,k_phase)*surv_like;
//   cout<<" obj_fun surv "<<obj_fun<<surv_like<<endl;
   obj_fun += column(disc_lambda,k_phase)*disc_like;
//   cout<<" obj_fun disc "<<obj_fun<<endl;
   obj_fun += column(mnwt_lambda,k_phase)*mnwt_like;
//   cout<<" obj_fun mnwt "<<obj_fun<<endl;
   obj_fun += column(length_lambda,k_phase)*length_like;
//   cout<<" obj_fun len "<<obj_fun<<endl;
   obj_fun += column(age_lambda,k_phase)*age_like;
//   cout<<" obj_fun age "<<obj_fun<<endl;
   obj_fun += column(sizeage_lambda,k_phase)*sizeage_like;
//   cout<<" obj_fun ms "<<obj_fun<<endl;

   obj_fun += equ_catch_like*init_equ_lambda(k_phase);
//   cout<<" obj_fun equ_cat "<<obj_fun<<endl;
   obj_fun += column(catch_lambda,k_phase)*catch_like;
//   cout<<" obj_fun catch "<<obj_fun<<catch_like<<endl;
   obj_fun += recr_like*recrdev_lambda(k_phase);
//   cout<<" obj_fun recr "<<obj_fun<<endl;
   obj_fun += parm_like*parm_prior_lambda(k_phase);
//   cout<<" obj_fun parm "<<obj_fun<<endl;
   obj_fun += parm_dev_like*parm_dev_lambda(k_phase);
//   cout<<" obj_fun parmdev "<<obj_fun<<endl;
   obj_fun += F_ballpark_like * F_ballpark_lambda(k_phase);
//   cout<<" obj_fun F_ballpark "<<obj_fun<<endl;
   obj_fun += CrashPen_lambda(k_phase)*CrashPen;
//   cout<<" obj_fun crash "<<obj_fun<<endl;
   obj_fun += square(dummy_datum-dummy_parm);
//   cout<<" obj_fun dummy "<<obj_fun<<endl;
   obj_fun += Fcast_recr_like;  //  lambda already factored in
//   cout<<" obj_fun forerecr "<<obj_fun<<endl;
   obj_fun += SoftBoundPen;
//   cout<<" obj_fun soft "<<obj_fun<<endl;
   if(SzFreq_Nmeth>0)  obj_fun += SzFreq_like*column(SzFreq_lambda,k_phase);
//   cout<<" obj_fun sizefreq "<<obj_fun<<endl;
   if(Do_Morphcomp>0) obj_fun += Morphcomp_lambda(k_phase)*Morphcomp_like;
   if(Do_TG>0 && Nfleet>1) obj_fun += TG_like1*column(TG_lambda1,k_phase);
   if(Do_TG>0) obj_fun += TG_like2*column(TG_lambda2,k_phase);
//   cout<<" obj_fun final "<<obj_fun<<endl;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 26 Process_STDquant */
FUNCTION void Process_STDquant()
  {
      if(rundetail>0 && mceval_counter==0) cout<<" Process STD quant "<<endl;
      for(y=styr-2; y<=YrMax;y++)
      {
        if(STD_Yr_Reverse(y)>0)
        {
          SPB_std(STD_Yr_Reverse(y))=SPB_yr(y);
          recr_std(STD_Yr_Reverse(y))=exp_rec(y,4);
        }
        if(STD_Yr_Reverse_Dep(y)>0) {depletion(STD_Yr_Reverse_Dep(y))=SPB_yr(y);}
      }

      if(rundetail>0 && mceval_counter==0) cout<<" STD OK "<<endl;

      switch(depletion_basis)
      {
        case 0:
        {
          depletion/=SPB_virgin;
          break;
        }
        case 1:
        {
          depletion/= (depletion_level*SPB_virgin);
          break;
        }
        case 2:
        {
          depletion/= (depletion_level*Bmsy);
          break;
        }
        case 3:
        {
          depletion/= (depletion_level*SPB_yr(styr));
          break;
        }
      }
      if(rundetail>0 && mceval_counter==0) cout<<" depletion OK "<<endl;

      switch (SPR_reporting)
      {
        case 0:      // keep as raw value
        {
          break;
        }
        case 1:  // compare to SPR
        {
//          SPR_std = (1.-SPR_std)/(1.-SPR_actual);
          SPR_std = (1.-SPR_std)/(1.-SPR_target);
          break;
        }
        case 2:  // compare to SPR_MSY
        {
          SPR_std = (1.-SPR_std)/(1.-MSY_SPR);
          break;
        }
        case 3:  // compare to SPR_Btarget
        {
          SPR_std = (1.-SPR_std)/(1.-SPR_Btgt);
          break;
        }
        case 4:
        {
          SPR_std = 1.-SPR_std;
          break;
        }
      }
      if(rundetail>0 && mceval_counter==0) cout<<" SPR OK "<<endl;

//  init_int Do_Forecast   //  0=none; 1=F(SPR); 2=F(MSY) 3=F(Btgt); 4=F(endyr); 5=Ave F (enter yrs); 6=read Fmult
//  Use the selected F method for the forecast as the denominator for the F_std ratio
      switch (F_report_basis)
      {
        case 0:      // keep as raw value
        {
          break;
        }
        case 1:  // compare to SPR
        {
          F_std /= Mgmt_quant(10);
          break;
        }
        case 2:  // compare to SPR_MSY
        {
          F_std /= Mgmt_quant(14);
          break;
        }
        case 3:  // compare to SPR_Btarget
        {
          F_std /= Mgmt_quant(7);
          break;
        }
      }

//  ProgLabel_7.8  get extra std quantities
    if(Selex_Std_Cnt>0)
    {
      for(i=1;i<=Selex_Std_Cnt;i++)
      {
        j=Selex_Std_Pick(i);
        if(Selex_Std_AL==1)
        {
          Extra_Std(i)=sel_l(Selex_Std_Year,Do_Selex_Std,1,j);
          if(gender==2) Extra_Std(i+Selex_Std_Cnt)=sel_l(Selex_Std_Year,Do_Selex_Std,2,j);
        }
        else
        {
          Extra_Std(i)=sel_a(Selex_Std_Year,Do_Selex_Std,1,j);
          if(gender==2) Extra_Std(i+Selex_Std_Cnt)=sel_a(Selex_Std_Year,Do_Selex_Std,2,j);
        }
      }
    }

    if(Growth_Std_Cnt>0)
    {
      for(i=1;i<=Growth_Std_Cnt;i++)
      {
        j=Growth_Std_Pick(i);  // selected age
        k=GP_finder(Do_Growth_Std,1);  // selected GP and gender  gp3
        Extra_Std(gender*Selex_Std_Cnt+i)=Ave_Size(t,1,k,j);
        if(gender==2)
        {
          k=GP_finder(Do_Growth_Std,2);  // selected GP and gender  gp3
          Extra_Std(gender*Selex_Std_Cnt+Growth_Std_Cnt+i)=Ave_Size(t,1,k,j);
        }
      }
    }

    if(NatAge_Std_Cnt>0)
    {
      if(Do_NatAge_Std<0)  // sum all areas
      {
        p1=1; p2=pop;
      }
      else  // selected area
      {
        p1=Do_NatAge_Std; p2=Do_NatAge_Std;
      }
      y=NatAge_Std_Year;
      t=styr+(y-styr)*nseas;  // first season of selected year
      for(i=1;i<=NatAge_Std_Cnt;i++)
      {
        a=NatAge_Std_Pick(i);  // selected age
        temp=0.;
        for(p=p1;p<=p2;p++)
        {
          for(g=1;g<=gmorph;g++)
          if(sx(g)==1 && use_morph(g)>0)
          {
            temp+=natage(t,p,g,a);  //  note, uses season 1 only
          }
        }
        Extra_Std(gender*(Selex_Std_Cnt+Growth_Std_Cnt)+i)=temp;
        if(gender==2)
        {
          temp=0.;
          for(p=p1;p<=p2;p++)
          {
            for(g=1;g<=gmorph;g++)
            if(sx(g)==2 && use_morph(g)>0)
            {
              temp+=natage(t,p,g,a);  //  note, uses season 1 only
            }
          }
          Extra_Std(gender*(Selex_Std_Cnt+Growth_Std_Cnt)+NatAge_Std_Cnt+i)=temp;
        }
      }
    }
    if(Extra_Std_N==1)
    {
      Extra_Std(1)=SPB_virgin;
    }
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 27 Check_Parm */
FUNCTION dvariable Check_Parm(const double& Pmin, const double& Pmax, const double& jitter, const prevariable& Pval)
  {
    dvariable NewVal;
    dvariable temp;
    NewVal=Pval;
    if(Pmin>Pmax)
    {N_warn++; cout<<" EXIT - see warning "<<endl; warning<<" parameter min > parameter max "<<Pmin<<" > "<<Pmax<<endl; cout<<" fatal error, see warning"<<endl; exit(1);}
    if(Pmin==Pmax)
    {N_warn++; warning<<" parameter min is same as parameter max, error condition for beta prior "<<Pmin<<" = "<<Pmax<<endl; NewVal=Pmin;}
    if(Pval<Pmin)
    {N_warn++; warning<<" parameter init value is less than parameter min "<<Pval<<" < "<<Pmin<<endl; NewVal=Pmin;}
    if(Pval>Pmax)
    {N_warn++; warning<<" parameter init value is greater than parameter max "<<Pval<<" > "<<Pmax<<endl; NewVal=Pmax;}
    if(jitter>0.0)
    {
      temp=log((Pmax-Pmin+0.0000002)/(NewVal-Pmin+0.0000001)-1.)/(-2.);   // transform the parameter
      temp += randn(radm) * jitter;
      NewVal=Pmin+(Pmax-Pmin)/(1.+mfexp(-2.*temp));
      if(Pmin==-99 || Pmax==99)
      {N_warn++; warning<<" use of jitter not advised unless parameter min & max are in reasonable parameter range "<<Pmin<<" "<<Pmax<<endl;}
    }
    return NewVal;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 28 Report_Parm */
FUNCTION void Report_Parm(const int NParm, const int AC, const int Activ, const prevariable& Pval, const double& Pmin, const double& Pmax, const double& RD, const double& PR, const int PR_T, const double& CV, const int PH, const prevariable& Like)
  {
    dvar_vector parm_val(1,14);
    dvar_vector prior_val(1,14);
    int i;
    dvariable parmvar;
    parmvar=0.0;
    SS2out<<NParm<<" "<<ParmLabel(NParm)<<" "<<Pval;
    if(Activ>0)
    {
      parmvar=CoVar(AC,1);
      SS2out<<" "<<AC<<" "<<PH<<" "<<Pmin<<" "<<Pmax<<" "<<RD;
      if (Pval==RD)
      {
        SS2out<<" NO_MOVE ";
      }
      else
      {
        temp=(Pval-Pmin)/(Pmax-Pmin);
        if(temp==0.0 || temp==1.0)
          {SS2out<<" BOUND "; Nparm_on_bound++;}
        else if(temp<0.01)
          {SS2out<<" LO "; Nparm_on_bound++;}
        else if(temp>=0.99)
          {SS2out<<" HI "; Nparm_on_bound++;}
        else
          {SS2out<<" OK ";}
      }
      SS2out<<" "<<parmvar;
    }
    else
    {
      SS2out<<" _ "<<PH<<" "<<Pmin<<" "<<Pmax<<" "<<RD<<" NA _ ";
    }
    if(PR_T>=0)
    {
      switch (PR_T)
      {
        case 0:
        {SS2out<<" Normal "; break;}
        case 1:
        {SS2out<<" Sym_Beta "; break;}
        case 2:
        {SS2out<<" Full_Beta "; break;}
        case 3:
        {SS2out<<" Log_Norm "; break;}
        case 4:
        {SS2out<<" Log_Norm_w/biasadj "; break;}
        case 5:
        {SS2out<<" Gamma "; break;}
      }
      SS2out<<" "<<PR<<" "<<CV<<" "<<Like<<" ";
      i=1;
      parm_val(i)=Pval;
      prior_val(i)=Get_Prior(PR_T, Pmin, Pmax, PR, CV, Pval);
      i=2;
      temp=Pval-1.96*parmvar;
      if(temp<Pmin) temp=Pmin;
      parm_val(i)=temp;
      prior_val(i)=Get_Prior(PR_T, Pmin, Pmax, PR, CV, temp);

      i=3;
      temp=Pval+1.96*parmvar;
      if(temp>Pmax) temp=Pmax;
      parm_val(i)=temp;
      prior_val(i)=Get_Prior(PR_T, Pmin, Pmax, PR, CV, temp);

      i=4;
      temp=Pmin+0.01*(Pmax-Pmin);
      parm_val(i)=temp;
      prior_val(i)=Get_Prior(PR_T, Pmin, Pmax, PR, CV, temp);
      i=14;
      temp=Pmax-0.01*(Pmax-Pmin);
      parm_val(i)=temp;
      prior_val(i)=Get_Prior(PR_T, Pmin, Pmax, PR, CV, temp);

      for(int i=5;i<=13;i++)
      {
        temp=Pmin+float(i-4)/10.0*(Pmax-Pmin);
        parm_val(i)=temp;
        prior_val(i)=Get_Prior(PR_T, Pmin, Pmax, PR, CV, temp);
      }
      SS2out<<parm_val<<" "<<prior_val;
    }
    else
    {SS2out<<" No_prior ";}
    SS2out<<endl;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 29 Get_Prior */
FUNCTION dvariable Get_Prior(const int T, const double& Pmin, const double& Pmax, const double& Pr, const double& Psd, const prevariable& Pval)
  {
    dvariable Prior_Like;
    dvariable mu;
    dvariable tau;
    dvariable Pconst=0.0001;
    dvariable Bprior;
    dvariable Aprior;
    switch (T)
    {
      case -1:  // none
      {
        Prior_Like = 0.;
        break;
      }
      case 0: // normal
      {
        Prior_Like = 0.5*square((Pval-Pr)/Psd);
        break;
      }
      case 1:  // symmetric beta    value of Psd must be >0.0
      {
        mu=-(Psd*(log( (Pmax+Pmin)*0.5- Pmin)))- (Psd*(log(0.5)));
        Prior_Like = -(mu+ (Psd*(log(Pval-Pmin+Pconst)))+(Psd*(log(1.-((Pval-Pmin-Pconst)/(Pmax-Pmin))))));
        break;
      }
      case 2:  // CASAL's Beta;  check to be sure that Aprior and Bprior are OK before running SS2!
      {
        mu=(Pr-Pmin) / (Pmax-Pmin);  // CASAL's v
        tau=(Pr-Pmin)*(Pmax-Pr)/square(Psd)-1.0;
        Bprior=tau*mu;  Aprior=tau*(1.0-mu);  // CASAL's m and n
        if(Bprior<=1.0 || Aprior <=1.0) {warning<<" bad Beta prior "<<Pval<<" "<<Pr<<endl;N_warn++;}
        Prior_Like =  (1.0-Bprior)*log(Pconst+Pval-Pmin) + (1.0-Aprior)*log(Pconst+Pmax-Pval)
              -(1.0-Bprior)*log(Pconst+Pr-Pmin) - (1.0-Aprior)*log(Pconst+Pmax-Pr);
        break;
      }
      case 3: // lognormal without bias correction
      {
        if(Pmin>0.0) {Prior_Like = 0.5*square((log(Pval)-Pr)/Psd);}
          else {N_warn++; warning<<" cannot do prior in log space for parm with min <=0.0" << endl;}
        break;
      }
      case 4: //lognormal with bias correction (from Larry Jacobson)
      {
        if(Pmin>0.0) Prior_Like=0.5*square((log(Pval)-Pr+0.5*square(Psd))/Psd);
          else {N_warn++; warning<<" cannot do prior in log space for parm with min <=0.0" << endl;}
        break;
      }
      case 5: //gamma  (from Larry Jacobson)
      {
        double warnif=1e-15;
        if(Pmin<0.0) {N_warn++; warning<<"Lower bound for gamma prior must be >=0.  Suggestion " << warnif*10.0 <<endl;}
        else 
        {
//Gamma is defined over [0,+inf) but x=zero causes trouble for some mean/variance combos.
          if(Pval < warnif) {N_warn++; warning<<"Pval too close to zero in gamma prior - can not guarantee reliable calculations.  Suggest rescaling data (e.g. * 1000)? "<<endl;}
          else
          {
            dvariable scale=square(Psd)/Pr;  // gamma parameters by method of moments
            dvariable shape=Pr/scale;
            Prior_Like= -shape*log(scale)-gammln(shape)+(shape-1.0)*log(Pval)-Pval/scale;
          }
        }
        break;
      }
    }
    return Prior_Like;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 30 Do_Equil_Calc */
FUNCTION void Do_Equil_Calc()
  {
  int t_base;
  int bio_t_base;
  int bio_t;
  dvariable N_mid;
  dvariable N_beg;
  dvariable Fishery_Survival;
  dvariable crashtemp;
  dvariable crashtemp1;
  dvar_matrix Survivors(1,pop,1,gmorph);
  dvar_matrix Survivors2(1,pop,1,gmorph);

   t_base=styr+(eq_yr-styr)*nseas-1;
   bio_t_base=styr+(bio_yr-styr)*nseas-1;
   GenTime.initialize(); Equ_penalty.initialize();
   cumF.initialize(); maxF.initialize();
   SPB_equil_pop_gp.initialize();
   if(Hermaphro_Option>0) MaleSPB_equil_pop_gp.initialize();
   equ_mat_bio=0.0;
   equ_mat_num=0.0;
   equ_catch_fleet.initialize();
   equ_numbers.initialize();
   equ_catage.initialize();
   equ_F_std=0.0;
   totbio=0.0;
   smrybio=0.0;
   smryage=0.0;
   smrynum=0.0;
     for (a=0;a<=3*nages;a++)     // go to 3x nages to approximate the infinite tail, then add the infinite tail
     {
       if(a<=nages) {a1=a;} else {a1=nages;}    // because selex and biology max out at nages

       for (s=1;s<=nseas;s++)
       {
         t=t_base+s;
         bio_t=bio_t_base+s;

         for (g=1;g<=gmorph;g++)
         if(use_morph(g)>0)
         {
           gg=sx(g);    // gender

           if(a==0 && s==Bseas(g))
           {
             if(y==styr-3)  // so in benchmark
             {
               for(p=1;p<=pop;p++)
               {equ_numbers(s,p,g,0) = equ_Recr*recr_dist_benchmark(GP(g),s,p)*submorphdist(GP2(g));}   // get the age 0 recruits for season=recr_seas
             }
             else
             {
               for(p=1;p<=pop;p++)
               {equ_numbers(s,p,g,0) = equ_Recr*recr_dist_calc(GP(g),s,p)*submorphdist(GP2(g));}   // get the age 0 recruits for season=recr_seas
             }
           }

           for (p=1;p<=pop;p++)
           {
             N_beg=equ_numbers(s,p,g,a);
             if(F_Method==1)   // Pope's approx
             {
                 N_mid = N_beg*surv1(s,GP3(g),a1);     // numbers at middle of season
                 Nsurvive=N_mid;                            // initial number of fishery survivors
                 if(Fishon==1)
                 {                       //  remove catch this round
                   // check to see if total harves would exceed max_harvest_rate
                   crashtemp=0.;  harvest_rate=1.0;
                   for (f=1;f<=Nfleet;f++)
                   if (fleet_area(f)==p && Hrate(f,t)>0.)
                   {
                     crashtemp+=Hrate(f,t)*deadfish(s,g,f,a1);
                   }

                   if(crashtemp>0.20)                  // only worry about this if the exploit rate is at all high
                   {
                     join1=1./(1.+mfexp(40.0*(crashtemp-max_harvest_rate)));  // steep joiner logistic curve at limit
                     upselex=1./(1.+mfexp(Equ_F_joiner*(crashtemp-0.2)));          //  value of a shallow logistic curve that goes through the limit
                     harvest_rate = join1 + (1.-join1)*upselex/(crashtemp);      // ratio by which all Hrates will be adjusted
                   }

                   for (f=1;f<=Nfleet;f++)
                   if (fleet_area(f)==p && Hrate(f,t)>0.)
                   {
                     temp=N_mid*Hrate(f,t)*harvest_rate;     // numbers that would be caught if fully selected
                     Nsurvive-=temp*deadfish(s,g,f,a1);       //  survival from fishery kill
                     equ_catch_fleet(f,2) += temp*deadfish_B(s,g,f,a1);
                     equ_catch_fleet(f,5) += temp*deadfish(s,g,f,a1);
                     equ_catch_fleet(f,3) += temp*sel_al_2(s,g,f,a1);      // retained fishery kill in biomass

                       equ_catch_fleet(f,1)+=temp*sel_al_1(s,g,f,a1);      //  total fishery encounter in biomass
                       equ_catch_fleet(f,4)+=temp*sel_al_3(s,g,f,a1);    // total fishery encounter in numbers
                       equ_catch_fleet(f,6)+=temp*sel_al_4(s,g,f,a1);      // retained fishery kill in numbers
                       equ_catage(s,f,g,a1)+=temp*deadfish(s,g,f,a1);      //  dead catch numbers per recruit  (later accumulate N in a1)
                   }
                 }   // end removing catch

                 Nsurvive *= surv1(s,GP3(g),a1);  // decay to end of season

                 if(a<=a1)
                 {
                   equ_Z(s,p,g,a1) = -(log((Nsurvive+1.0e-13)/(N_beg+1.0e-10)))/seasdur(s);
                   Fishery_Survival = equ_Z(s,p,g,a1)-natM(s,GP3(g),a1);
                   if(a>=Smry_Age)
                   {
                     cumF(g)+=Fishery_Survival*seasdur(s);
                     if(Fishery_Survival>maxF(g)) maxF(g)=Fishery_Survival;
                   }
                 }

             }   // end Pope's approx

             else          // Continuous F for method 2 or 3
             {
               equ_Z(s,p,g,a1)=natM(s,GP3(g),a1);
                 if(Fishon==1)
                 {
                   if(a1<=nages)
                   {
                     for (f=1;f<=Nfleet;f++)       //loop over fishing fleets to get Z
                     if (fleet_area(f)==p && Hrate(f,t)>0.0)
                     {
                       equ_Z(s,p,g,a1)+=deadfish(s,g,f,a1)*Hrate(f,t);
                     }
                     if(save_for_report>0)
                     {
                       temp=equ_Z(s,p,g,a1)-natM(s,GP3(g),a1);
                       if(a>=Smry_Age && a<=nages) cumF(g)+=temp*seasdur(s);
                       if(temp>maxF(g)) maxF(g)=temp;
                     }
                   }
                 }
                 Nsurvive=N_beg*mfexp(-seasdur(s)*equ_Z(s,p,g,a1));
             }  //  end F method
             Survivors(p,g)=Nsurvive;
           }  // end pop
         }  // end morph

         if(Hermaphro_Option>0)
         {
           if(Hermaphro_seas==-1 || Hermaphro_seas==s)
           {
             for(p=1;p<=pop;p++)
             {
               k=gmorph/2;
               for(g=1;g<=k;g++)
               if(use_morph(g)>0)
               {
                 Survivors(p,g+k) += Survivors(p,g)*Hermaphro_val(GP4(g),a1); // increment males
                 Survivors(p,g) *= (1.-Hermaphro_val(GP4(g),a1)); // decrement females
               }
             }
           }
         }
          if(do_migration>0)  // movement between areas in equil calcs
          {
            Survivors2.initialize();
            for(g=1;g<=gmorph;g++)
            if(use_morph(g)>0)
            {
              for(p=1;p<=pop;p++)
              for(p2=1;p2<=pop;p2++)
              {
                k=move_pattern(s,GP4(g),p,p2);
                if(k>0) Survivors2(p2,g) += Survivors(p,g)*migrrate(bio_yr,k,a1);
              }  // end destination pop
            }
            Survivors=Survivors2;
          }  // end do migration

          for(g=1;g<=gmorph;g++)
          if(use_morph(g)>0)
          {
            for(p=1;p<=pop;p++)
            {
              if(s==nseas)  // into next age at season 1
              {
                if(a==3*nages)
                {
                  // end of the cohort
                }
                else if(a==(3*nages-1))           // do infinite tail; note that it uses Z from nseas as if it applies annually
                {
                  if(F_Method==1)
                  {
  //                  equ_numbers(1,p,g,a+1) = Survivors(p,g);
  //                  equ_numbers(1,p,g,a+1) = Survivors(p,g)/(1.-exp(-natM(s,GP3(g),nages)));
                    equ_numbers(1,p,g,a+1) = Survivors(p,g)/(1.-exp(-equ_Z(nseas,p,g,nages)));
                  }
                  else
                  {
                    equ_numbers(1,p,g,a+1) = Survivors(p,g)/(1.-exp(-equ_Z(nseas,p,g,nages)));
                  }
                }
                else
                {
                  equ_numbers(1,p,g,a+1) = Survivors(p,g);
                }
              }
              else
              {
                equ_numbers(s+1,p,g,a) = Survivors(p,g);  // same age, next season
              }
            }
          }
        }  // end season
      }  // end age

// now calc contribution to catch and ssb
       for(g=1;g<=gmorph;g++)
       if(use_morph(g)>0)
       {
         gg=sx(g);
         for(s=1;s<=nseas;s++)
         for(p=1;p<=pop;p++)
         {
           t=t_base+s;
           bio_t=bio_t_base+s;
           equ_numbers(s,p,g,nages)+=sum(equ_numbers(s,p,g)(nages+1,3*nages));
           if(Fishon==1)
           {
             if(F_Method>=2)
             {
               Zrate2(p,g)=elem_div( (1.-mfexp(-seasdur(s)*equ_Z(s,p,g))), equ_Z(s,p,g));
               if(s<Bseas(g)) Zrate2(p,g,0)=0.0;
               for (f=1;f<=Nfleet;f++)
               if (fleet_area(f)==p && Hrate(f,t)>0.0)
               {
                 equ_catch_fleet(f,2)+=Hrate(f,t)*elem_prod(equ_numbers(s,p,g)(0,nages),deadfish_B(s,g,f))*Zrate2(p,g);      // dead catch bio
                 equ_catch_fleet(f,5)+=Hrate(f,t)*elem_prod(equ_numbers(s,p,g)(0,nages),deadfish(s,g,f))*Zrate2(p,g);      // deadfish catch numbers
                 equ_catch_fleet(f,3)+=Hrate(f,t)*elem_prod(equ_numbers(s,p,g)(0,nages),sel_al_2(s,g,f))*Zrate2(p,g);      // retained catch bio
                   equ_catage(s,f,g)=elem_prod(elem_prod(equ_numbers(s,p,g)(0,nages),deadfish(s,g,f)) , Zrate2(p,g));
                   equ_catch_fleet(f,1)+=Hrate(f,t)*elem_prod(equ_numbers(s,p,g)(0,nages),sel_al_1(s,g,f))*Zrate2(p,g);      // encountered catch bio
                   equ_catch_fleet(f,4)+=Hrate(f,t)*elem_prod(equ_numbers(s,p,g)(0,nages),sel_al_3(s,g,f))*Zrate2(p,g);      // encountered catch bio
                   equ_catch_fleet(f,6)+=Hrate(f,t)*elem_prod(equ_numbers(s,p,g)(0,nages),sel_al_4(s,g,f))*Zrate2(p,g);      // retained catch numbers
               }
             }
             else  // F_method=1
             {
               // already done in the age loop
             }
           }

           if(s==1)
           {
             totbio += equ_numbers(s,p,g)(0,nages)*Wt_Age_beg(s,g)(0,nages);
             smrybio += equ_numbers(s,p,g)(Smry_Age,nages)*Wt_Age_beg(s,g)(Smry_Age,nages);
             smrynum += sum(equ_numbers(s,p,g)(Smry_Age,nages));
             smryage += equ_numbers(s,p,g)(Smry_Age,nages) * r_ages(Smry_Age,nages);
           }
           if(s==spawn_seas)
           {
            //  need to adjust for spawn_time_seas
             if(gg==1)  // compute equilibrium spawning biomass for females
             {
               SPB_equil_pop_gp(p,GP4(g))+=elem_prod(equ_numbers(s,p,g)(0,nages),mfexp(-spawn_time_seas*equ_Z(s,p,g)(0,nages)))*fec(g);
               equ_mat_bio+=elem_prod(equ_numbers(s,p,g)(0,nages),mfexp(-spawn_time_seas*equ_Z(s,p,g)(0,nages)))*make_mature_bio(GP4(g));
               equ_mat_num+=elem_prod(equ_numbers(s,p,g)(0,nages),mfexp(-spawn_time_seas*equ_Z(s,p,g)(0,nages)))*make_mature_numbers(GP4(g));
               GenTime+=equ_numbers(s,p,g)(0,nages)*elem_prod(fec(g),r_ages);
             }
             else if(Hermaphro_Option>0 && gg==2)
             {
               MaleSPB_equil_pop_gp(p,GP4(g))+=elem_prod(equ_numbers(s,p,g)(0,nages),mfexp(-spawn_time_seas*equ_Z(s,p,g)(0,nages)))*Wt_Age_beg(s,g)(0,nages);
               equ_mat_bio+=elem_prod(equ_numbers(s,p,g)(0,nages),mfexp(-spawn_time_seas*equ_Z(s,p,g)(0,nages)))*Wt_Age_beg(s,g);
               equ_mat_num+=sum(elem_prod(equ_numbers(s,p,g)(0,nages),mfexp(-spawn_time_seas*equ_Z(s,p,g)(0,nages))));
             }
           }
         }
       }

     YPR_dead = sum(column(equ_catch_fleet,2));           // dead yield per recruit
     YPR_N_dead = sum(column(equ_catch_fleet,5));           // dead numbers per recruit
     YPR_enc = sum(column(equ_catch_fleet,1));               //  encountered yield per recruit
     YPR_ret = sum(column(equ_catch_fleet,3));           // retained yield per recruit

   if(Fishon==1)
   {
     if(F_reporting<=1)
     {
       equ_F_std=YPR_dead/smrybio;
     }
     else if(F_reporting==2)
     {
       equ_F_std=YPR_N_dead/smrynum;
     }
     else if(F_reporting==3)
     {
       if(F_Method==1)
       {
         for(s=1;s<=nseas;s++)
         {
           t=t_base+s;
           for (f=1;f<=Nfleet;f++)
           {
             equ_F_std+=Hrate(f,t);
           }
         }
       }
       else
       {
         for(s=1;s<=nseas;s++)
         {
           t=t_base+s;
           for (f=1;f<=Nfleet;f++)
           {
             equ_F_std+=Hrate(f,t)*seasdur(s);
           }
         }
       }
     }
     else if(F_reporting==4)
     {
       temp1=0.0;
       temp2=0.0;
       for(g=1;g<=gmorph;g++)
       if(use_morph(g)>0)
       {
         for(p=1;p<=pop;p++)
         {
           for(a=F_reporting_ages(1);a<=F_reporting_ages(2);a++)   //  should not let a go higher than nages-2 because of accumulator
           {
             if(nseas==1)
             {
               temp1+=equ_numbers(1,p,g,a+1);
               temp2+=equ_numbers(1,p,g,a)*mfexp(-seasdur(1)*natM(1,GP3(g),a));
             }
             else
             {
               temp1+=equ_numbers(1,p,g,a+1);
               temp3=equ_numbers(1,p,g,a);  //  numbers at begin of year
               for(int kkk=1;kkk<=nseas;kkk++) {temp3*=mfexp(-seasdur(kkk)*natM(kkk,GP3(g),a));}
               temp2+=temp3;
             }
           }
         }
       }
       equ_F_std = log(temp2)-log(temp1);
     }
   }
   SPB_equil=sum(SPB_equil_pop_gp);
   GenTime/=SPB_equil;
   smryage /= smrynum;
   cumF/=(r_ages(nages)-r_ages(Smry_Age)+1.);
   if(Hermaphro_maleSPB==1) SPB_equil+=sum(MaleSPB_equil_pop_gp);
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 31 Make_AgeLength_Key */
FUNCTION void Make_AgeLength_Key()
  {
  dvar_vector templen(2,nlength);
  dvar_vector AL(1,nlength1);
  dvar_vector use_Ave_Size_W(0,nages);
  AL(1)=0.0; AL(nlength1)=1.0;
  int gstart;
  int gpattern;
  int sex;
   smid=s+nseas;
   ALK_yr=y;
   gp=0;
//   echoinput<<y<<" "<<s<<" "<<g<<" sdsize_within_start "<<Sd_Size_within(s,gstart)<<endl;
//   echoinput<<y<<" "<<s<<" "<<g<<" sdsize_within_mid "<<Sd_Size_within(smid,gstart)<<endl;
    for(sex=1;sex<=gender;sex++)
    for(gpattern=1;gpattern<=N_GP;gpattern++)
    {
      gp=gp+1;
      gstart=(gp-1)*nseas*N_platoon+int(N_platoon/2)+1-N_platoon;
      for(birthseas=1;birthseas<=nseas;birthseas++)
      {
        gstart+=N_platoon;
        if(recr_dist_pattern_2(gpattern,birthseas)>0)
        {
          for(gp2=1;gp2<=N_platoon;gp2++)      // loop the platoons
          {
            g=gstart+ishadow(gp2);
            use_Ave_Size_W=Ave_Size(t,0,gstart);
            if(N_platoon>1) use_Ave_Size_W += shadow(gp2)*Sd_Size_between(s,gstart);
            dvar_matrix ALK_w_beg=ALK(s,g);        //  shallow copy
            gg=sx(g);
            
            if(Grow_logN==0)
            {
              for (a=0;a<=nages;a++)
              {
                templen=(len_bins(2,nlength)-use_Ave_Size_W(a))/(Sd_Size_within(s,gstart,a));
                for (z=2;z<=nlength;z++) {AL(z)= cumd_norm(templen(z));}
                ALK_w_beg(a)=first_difference(AL);
              }   // end age loop
            }
            else    //  lognormal size-at-age
            {
              for (a=0;a<=nages;a++)
              {
                temp=log(use_Ave_Size_W(a))-0.5*Sd_Size_within(s,gstart,a)*Sd_Size_within(s,gstart,a);
//  line below erroneously used sd_size_within(smid) through version 3.24S
                templen=(log(len_bins(2,nlength))-temp)/(Sd_Size_within(s,gstart,a));
                for (z=2;z<=nlength;z++) {AL(z)= cumd_norm(templen(z));}
              ALK_w_beg(a)=first_difference(AL);
              }   // end age loop
            }

            if(s==spawn_seas && gg==1)  // get fecundity at age from ALK and length-spec biology
            {
              switch(Maturity_Option)
              {
                case 1:  //  Maturity_Option=1  length logistic
                {
                  fec(g) = elem_prod(ALK_w_beg*mat_fec_len(gpattern),mat_age(gpattern));  //  reproductive output at age
                  make_mature_numbers(g)=elem_prod(ALK_w_beg*mat_len(gpattern),mat_age(gpattern));  //  mature numbers at age
                  make_mature_bio(g)=elem_prod(ALK_w_beg*elem_prod(mat_len(gpattern),wt_len(s,gg)),mat_age(gpattern));  //  mature biomass at age
                  break;
                }
                case 2:  //  Maturity_Option=2  age logistic
                {
                  fec(g) = elem_prod(ALK_w_beg*mat_fec_len(gpattern),mat_age(gpattern));  //  reproductive output at age
                  make_mature_numbers(g)=elem_prod(ALK_w_beg*mat_len(gpattern),mat_age(gpattern));  //  mature numbers at age
                  make_mature_bio(g)=elem_prod(ALK_w_beg*elem_prod(mat_len(gpattern),wt_len(s,gg)),mat_age(gpattern));  //  mature biomass at age
                  break;
                }
                case 3:  //  Maturity_Option=3  read age-maturity
                {
                  fec(g) = elem_prod(ALK_w_beg*mat_fec_len(gpattern),mat_age(gpattern));  //  reproductive output at age
                  make_mature_numbers(g)=elem_prod(ALK_w_beg*mat_len(gpattern),mat_age(gpattern));  //  mature numbers at age (Age_Maturity already copied to mat_age)
                  make_mature_bio(g)=elem_prod(ALK_w_beg*elem_prod(mat_len(gpattern),wt_len(s,gg)),mat_age(gpattern));  //  mature biomass at age
                  break;
                }
                case 4:  //  Maturity_Option=4   read age-fecundity, so no age-maturity
                {
                  fec(g)=Age_Maturity(gpattern);
                  make_mature_numbers(g)=fec(g);  //  not defined
                  make_mature_bio(g)=fec(g);   //  not defined
                  break;
                }
                case 5:  //  Maturity_Option=5   read age-fecundity from wtatage.ss
                {
                  fec(g)=WTage_emp(t,GP3(g),-2);
                  make_mature_numbers(g)=fec(g);  //  not defined
                  make_mature_bio(g)=fec(g);   //  not defined
                  break;
                }
                case 6:  //  Maturity_Option=6   read length-maturity
                {
                  fec(g) = elem_prod(ALK_w_beg*mat_fec_len(gpattern),mat_age(gpattern));  //  reproductive output at age
                  make_mature_numbers(g)=elem_prod(ALK_w_beg*mat_len(gpattern),mat_age(gpattern));  //  mature numbers at age (Length_Maturity already copied to mat_len)
                  make_mature_bio(g)=elem_prod(ALK_w_beg*elem_prod(mat_len(gpattern),wt_len(s,gg)),mat_age(gpattern));  //  mature biomass at age
                  break;
                }
              }

                save_sel_fec(t,g,0)= fec(g);   //  save sel_al_3 and save fecundity for output
                if(y==endyr) save_sel_fec(t+nseas,g,0)=fec(g);
                if(save_for_report==2 && ishadow(GP2(g))==0) bodywtout<<-y<<" "<<s<<" "<<gg<<" "<<GP4(g)<<" "<<Bseas(g)<<" "<<-2<<" "<<fec(g)<<endl;
            }

            if(WTage_rd==0)
            {Wt_Age_beg(s,g)=(ALK_w_beg*wt_len(s,gg));}   // wt-at-age at beginning of period
            else
            {Wt_Age_beg(s,g)=WTage_emp(t,GP3(g),0);}

            Save_Wt_Age(t,g)=Wt_Age_beg(s,g);
            if(save_for_report==2 && ishadow(GP2(g))==0) bodywtout<<-y<<" "<<s<<" "<<gg<<" "<<GP4(g)<<" "<<Bseas(g)<<" "<<0<<" "<<Wt_Age_beg(s,g)<<endl;

//  now do middle of season
            dvar_matrix ALK_w=ALK(smid,g);        //  shallow copy
            use_Ave_Size_W=extract_row(Ave_Size(t,1),gstart);
            if(N_platoon>1) use_Ave_Size_W += shadow(gp2)*Sd_Size_between(smid,gstart);
            if(Grow_logN==0)
            {
              for (a=0; a<=nages;a++)
              {
                templen=(len_bins(2,nlength)-use_Ave_Size_W(a))/(Sd_Size_within(smid,gstart,a));
                for (z=2;z<=nlength;z++) {AL(z)=cumd_norm(templen(z));}
                     //  note that AL(1) = 0.0 and AL(nlength+1)=1.0
                ALK_w(a)=first_difference(AL);
                if(WTage_rd==0)
                {Wt_Age_mid(s,g)=ALK_w*wt_len(s,gg);}  // use for fisheries with no size selectivity
                else
                {Wt_Age_mid(s,g)=WTage_emp(t,GP3(g),-1);}
              }   // end age loop
            }
            else  //  lognormal size-at-age
            {
              for (a=0; a<=nages;a++)
              {
                temp=log(use_Ave_Size_W(a))-0.5*Sd_Size_within(smid,gstart,a)*Sd_Size_within(smid,gstart,a);
                templen=(log(len_bins(2,nlength))-temp)/(Sd_Size_within(smid,gstart,a));
                for (z=2;z<=nlength;z++) {AL(z)=cumd_norm(templen(z));}
                     //  note that AL(1) = 0.0 and AL(nlength+1)=1.0
                ALK_w(a)=first_difference(AL);
                if(WTage_rd==0)
                {Wt_Age_mid(s,g)=ALK_w*wt_len(s,gg);}  // use for fisheries with no size selectivity
                else
                {Wt_Age_mid(s,g)=WTage_emp(t,GP3(g),-1);}
              }   // end age loop
            }

           if(save_for_report==2 && ishadow(GP2(g))==0) bodywtout<<-y<<" "<<s<<" "<<gg<<" "<<GP4(g)<<" "<<Bseas(g)<<" "<<-1<<" "<<Wt_Age_mid(s,g)<<endl;
          }  // end platoon loop
        }
      }   // end birthseason loop

    }  // end growth pattern&gender loop
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 32 Make_FishSelex */
FUNCTION void Make_FishSelex()
  {
    smid=s+nseas;

    dvar_matrix ALK_w=ALK(smid,g);        //  shallow copy
    int yf;
    int tz;

    gg=sx(g);
    if(y>endyr) {yz=endyr; } else {yz=y;}
    if(y>endyr+1) {yf=endyr+1;} else {yf=y;}    //  yf stores in endyr+1 the average selex from a range of years
    tz=styr+(y-styr)*nseas+s-1;  // can use y, not yf, because wtage_emp values are read in and can extend into forecast
    for(f=1;f<=Nfleet;f++)
    {
//      if(time_vary_sel(yz,f+Ntypes)>0 || time_vary_sel(yz,f)>0 || time_vary_MG(yz,2)>0 || time_vary_MG(yz,3)>0)
      if(time_vary_makefishsel(yf,f)>0 || save_for_report>0)
      {
        makefishsel_yr = yf;
        if (WTage_rd==1)
        {
          sel_al_1(s,g,f)=elem_prod(sel_a(yf,f,gg),WTage_emp(tz,GP3(g),f));   // Wt_Age_mid has been set to WTage_emp if necessary already
          sel_al_3(s,g,f)=sel_a(yf,f,gg);
        }
        else if(seltype(f,1)==0)  // no size_selectivity
        {
          sel_al_1(s,g,f)=elem_prod(sel_a(yf,f,gg),Wt_Age_mid(s,g));   // Wt_Age_mid has been set to WTage_emp if necessary already
          sel_al_3(s,g,f)=sel_a(yf,f,gg);
        }
        else
        {
          sel_al_1(s,g,f)=elem_prod(sel_a(yf,f,gg),(ALK_w * elem_prod( sel_l(yf,f,gg),wt_len(s,gg) )) );
          sel_al_3(s,g,f)=elem_prod(sel_a(yf,f,gg),(ALK_w * sel_l(yf,f,gg) ));
        }

        if(mceval_phase() || save_for_report>0)
        {
          if(WTage_rd==0)
          {fish_body_wt(tz,g,f)=elem_div( ALK_w*elem_prod(sel_l(yf,f,gg),wt_len(s,gg)) , ALK_w*sel_l(yf,f,gg));}
          else
          {fish_body_wt(tz,g,f)=WTage_emp(tz,GP3(g),f);}
          if(save_for_report==2 && ishadow(GP2(g))==0) bodywtout<<-y<<" "<<s<<" "<<gg<<" "<<GP4(g)<<" "<<Bseas(g)
          <<" "<<f<<" "<<fish_body_wt(tz,g,f)<<endl;
        }

        if(seltype(f,2)!=0)
        {
          sel_al_4(s,g,f)=elem_prod(sel_a(yf,f,gg),(ALK_w * sel_l_r(yf,f,gg) ));
          if(WTage_rd==0)
            {sel_al_2(s,g,f)=elem_prod(sel_a(yf,f,gg),(ALK_w * elem_prod( sel_l_r(yf,f,gg),wt_len(s,gg) )) );}
          else
            {sel_al_2(s,g,f)=elem_prod(sel_al_4(s,g,f),WTage_emp(tz,GP3(g),f));}  //   revised 3.24Q
        }
        else
        {
            sel_al_2(s,g,f)=sel_al_1(s,g,f);
            sel_al_4(s,g,f)=sel_al_3(s,g,f);
        }

        if(seltype(f,2)>=2)
        {
          deadfish(s,g,f)=elem_prod(sel_a(yf,f,gg),(ALK_w * discmort2(yf,f,gg)));  //  selected dead by numbers
          if(WTage_rd==0)
          {deadfish_B(s,g,f)=elem_prod(sel_a(yf,f,gg),(ALK_w * elem_prod(discmort2(yf,f,gg),wt_len(s,gg)) ));} // selected dead by weight
          else
          {deadfish_B(s,g,f)=sel_al_2(s,g,f);} // not quite correct, for now set equal to selected wt of retained fish without adjusting for discmort
        }
        else
        {
            deadfish_B(s,g,f)=sel_al_1(s,g,f);
            deadfish(s,g,f)=sel_al_3(s,g,f);
        }

      }  // end need to do it
      save_sel_fec(t,g,f)= value(sel_al_3(s,g,f));  //  save sel_al_3 in save_fecundity array for output

    }  // end fleet loop for mortality, retention
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 33 get_posteriors (MCMC eval) */
FUNCTION void get_posteriors()
  {
  if(rundetail>1) cout<<" mceval counter: "<<mceval_counter<<endl;
  if(rundetail==0 & double(mceval_counter)/200.==double(mceval_counter/200.)) cout<<" mceval counter: "<<mceval_counter<<endl;

  if(mceval_header==0 && mceval_phase())    // first pass through the mceval phase
  {
    // delete any old mcmc output files
    // will generate a warning if no files exist
    // but will play through just fine
    system("del rebuild.sso");
    system("del posteriors.sso");
    system("del derived_posteriors.sso");
    system("del posterior_vectors.sso");
    if(rundetail>0) cout<<" did system commands "<<endl;
  };
  // define the mcmc output files;
  ofstream rebuilder("rebuild.sso",ios::app);
  ofstream posts("posteriors.sso",ios::app);
  ofstream der_posts("derived_posteriors.sso",ios::app);
  ofstream post_vecs("posterior_vectors.sso",ios::app);

  // define the mcmc output files;
//  ofstream rebuilder("rebuild.sso",ios::app);
//  ofstream posts("posteriors.sso",ios::app);
//  ofstream der_posts("derived_posteriors.sso",ios::app);
//  ofstream post_vecs("posterior_vectors.sso",ios::app);

  if(mceval_header==0)    // first pass through the mceval phase
  {
    mceval_header=1;
    // produce the appropriate headers for the posteriors.rep
    // and derived_posteriors.rep files
    // parameters.rep matches "PARAMETERS" section in Report.SSO file
    if(rundetail>0) cout<<" write mcmc headers "<<endl;
    posts<<"Iter Objective_function ";
    for(i=1;i<=active_count;i++) {posts<<" "<<ParmLabel(active_parm(i));}
    posts << endl;

    // derived quantities
    // derived_parameters.rep matches "DERIVED_PARAMETERS" section in Report.SSO file
    NP = ParCount;
    der_posts<<"Iter Objective_function ";
    for(j=1;j<=N_STD_Yr;j++) // spawning biomass
    {
      NP++;  der_posts<<ParmLabel(NP)<<" ";
    }
    for(j=1;j<=N_STD_Yr;j++) // recruitment
    {
      NP++;  der_posts<<ParmLabel(NP)<<" ";
    }
    for(j=1;j<=N_STD_Yr_Ofish;j++) // SPRratio
    {
      NP++;  der_posts<<ParmLabel(NP)<<" ";
    }
    for(j=1;j<=N_STD_Yr_F;j++) // F
    {
      NP++;  der_posts<<ParmLabel(NP)<<" ";
    }
    for(j=1;j<=N_STD_Yr_Dep;j++) // depletion (Bratio)
    {
      NP++;  der_posts<<ParmLabel(NP)<<" ";
    }
    for(j=1;j<=N_STD_Mgmt_Quant;j++) // Management quantities
    {
      NP++;  der_posts<<ParmLabel(NP)<<" ";
    }
    for(j=1;j<=Extra_Std_N;j++)
    {
      NP++;  der_posts<<ParmLabel(NP)<<" ";
    }
    der_posts << endl;

    if(depletion_basis!=2) post_vecs<<"depletion_basis_is_not_=2;_so_info_below_is_not_B/Bmsy"<<endl;
    if(F_report_basis!=2) post_vecs<<"F_report_basis_is_not_=2;_so_info_below_is_not_F/Fmsy"<<endl;
    post_vecs<<"Endyr+1= "<<endyr+1<<endl;
    post_vecs<<"run mceval objfun Numbers Area Sex Ages:"<<age_vector<<endl;
    post_vecs<<"run mceval objfun F_yr ";
    for(y=styr-1;y<=YrMax; y++)
    {
      if(STD_Yr_Reverse_F(y)>0) post_vecs<<y<<" ";
    }
    post_vecs<<endl;
    post_vecs<<"run mceval objfun B_yr ";
    for(y=styr-1;y<=YrMax; y++)
    {
      if(STD_Yr_Reverse_Dep(y)>0) post_vecs<<y<<" ";
    }
    post_vecs<<endl;
    
  };  //  end writing headers for mceval_counter==1


  // produce standard output of all estimated parameters
  posts<<mceval_counter<<" "<<obj_fun<<" ";

  for(j=1;j<=N_MGparm2;j++)
  {
    if(active(MGparm(j))) posts<<MGparm(j)<<" ";
  }
  if(N_MGparm_dev>0)
  {
    for(i=1;i<=N_MGparm_dev;i++)
    for(j=MGparm_dev_minyr(i);j<=MGparm_dev_maxyr(i);j++)
    {
      if(active(MGparm_dev)) posts<<MGparm_dev(i,j)<<" ";
    }
  }
  for(i=1;i<=N_SRparm2;i++)
  {
    if(active(SR_parm(i))) posts<<SR_parm(i)<<" ";
  }

  if(recdev_cycle>0)
  {
    for(i=1;i<=recdev_cycle;i++)
    {
      if(active(recdev_cycle_parm(i))) posts<<recdev_cycle_parm(i)<<" ";
    }

  }
  if(recdev_do_early>0)
  {
    for(i=recdev_early_start;i<=recdev_early_end;i++)
    {
      if( active(recdev_early) ) posts<<recdev(i)<<" ";
    }
  }
  if(do_recdev>0)
  {
    for(i=recdev_start;i<=recdev_end;i++)
    {
      if( active(recdev1)||active(recdev2) ) posts<<recdev(i)<<" ";
    }
  }
  if(Do_Forecast>0)
  {
    for(i=recdev_end+1;i<=YrMax;i++)
    {
      if(active(Fcast_recruitments)) posts<<Fcast_recruitments(i)<<" ";
    }
    for(i=endyr+1;i<=YrMax;i++)
    {
      if(active(Fcast_impl_error)) posts<<Fcast_impl_error(i)<<" ";
    }

  }
  for(i=1;i<=Nfleet;i++)
  {
    if(active(init_F(i))) posts<<init_F(i)<<" ";
  }
  if(F_Method==2)
  {
    for(i=1;i<=Nfleet*(TimeMax-styr+1);i++)
    {
      if(active(F_rate(i))) posts<<F_rate(i)<<" ";
    }
  }
  for(i=1;i<=Q_Npar;i++)
  {
    if(active(Q_parm(i)))posts<<Q_parm(i)<<" ";
  }
  for (j=1;j<=N_selparm2;j++)
  {
    if(active(selparm(j))) posts<<selparm(j)<<" ";
  }
  for(i=1;i<=N_selparm_dev;i++)
  for(j=selparm_dev_minyr(i);j<=selparm_dev_maxyr(i);j++)
  {
    if(active(selparm_dev)) posts<<selparm_dev(i,j)<<" ";
  }
  if(Do_TG>0)
  {
    k=3*N_TG+2*Nfleet;
    for (j=1;j<=k;j++)
    {
      if(active(TG_parm(j))) posts<<TG_parm(j)<<" ";
    }
  }
  posts << endl;

  // derived quantities
  der_posts<<mceval_counter<<" "<<obj_fun<<" ";
  for(j=1;j<=N_STD_Yr;j++) // spawning biomass
  {
    der_posts<<SPB_std(j)<<" ";
  }
  for(j=1;j<=N_STD_Yr;j++) // recruitment
  {
    der_posts<<recr_std(j)<<" ";
  }
  for(j=1;j<=N_STD_Yr_Ofish;j++) // SPRratio
  {
    der_posts<<SPR_std(j)<<" ";
  }
  for(j=1;j<=N_STD_Yr_F;j++) // F
  {
    der_posts<<F_std(j)<<" ";
  }
  for(j=1;j<=N_STD_Yr_Dep;j++) // depletion (Bratio)
  {
    der_posts<<depletion(j)<<" ";
  }
  for(j=1;j<=N_STD_Mgmt_Quant;j++) // Management quantities
  {
    der_posts<<Mgmt_quant(j)<<" ";
  }
  for(j=1;j<=Extra_Std_N;j++)
  {
    der_posts<<Extra_Std(j)<<" ";
  }

  der_posts << endl;

  if(Do_Rebuilder>0) write_rebuilder_output();

  // derived vectors quantities
  t=styr+(endyr+1-styr)*nseas;
  for(p=1;p<=pop;p++)
  {
    for (gg=1;gg<=gender;gg++)
    {
      tempvec_a.initialize();
      for(g=1;g<=gmorph;g++)
      if(sx(g)==gg && use_morph(g)>0)
      {
        tempvec_a+=natage(t,p,g);
      }
      post_vecs<<runnumber<<" "<<mceval_counter<<" "<<obj_fun<<" N_at_Age "<<" "<<p<<" "<<gg<<" "<<tempvec_a<<endl;
    }
  }
  post_vecs<<runnumber<<" "<<mceval_counter<<" "<<obj_fun<<" F/Fmsy "<<F_std<<endl;
  post_vecs<<runnumber<<" "<<mceval_counter<<" "<<obj_fun<<" B/Bmsy "<<depletion<<endl;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 34 Get_Benchmarks;  Find Fspr, MSY */
FUNCTION void Get_Benchmarks()
  {
  int jj;  int Nloops;
  int bio_t;
  int bio_t_base;
  dvariable last_F1;  dvariable Closer;
  dvariable Vbio1_unfished;
  dvariable Vbio_MSY;
  dvariable Vbio1_MSY;
  dvariable SPR_at_target;
  dvariable junk; dvariable Nmid_c;

  dvariable df;
  dvariable BestYield;
  dvariable BestF1;
  dvar_vector F1(1,3);
  dvariable FF;
  dvar_vector yld1(1,3);
  dvariable dyld;
  dvariable dyldp;
  dvariable Fmax;
  dvariable bestF1;
  dvariable bestF2;

  show_MSY=0;
  if(mceval_phase()==0) {show_MSY=1;}
  if(show_MSY==1)
  {
    report5<<version_info_short<<endl;
    report5<<version_info<<endl<<ctime(&start);
  }
      maxpossF.initialize();
      for(g=1;g<=gmorph;g++)
        for(s=1;s<=nseas;s++)
        {
          tempvec_a.initialize();
          for(f=1;f<=Nfleet;f++) {tempvec_a+=Bmark_RelF_Use(s,f)*deadfish(s,g,f);}
          temp=max(tempvec_a);
          if(temp>maxpossF) maxpossF=temp;
        }
        maxpossF =max_harvest_rate/maxpossF;    //  applies to any F_method
        report5<<"Calculated_Max_Allowable_F "<<maxpossF<<endl<<"Bmark_relF(by_fleet_&seas)"<<endl<<Bmark_RelF_Use<<endl<<"#"<<endl;

    y=styr-3;  //  the average biology from specified benchmark years is stored here
    yz=y;
    bio_yr=y;
    eq_yr=y;
    t_base=y+(y-styr)*nseas-1;
    bio_t_base=styr+(bio_yr-styr)*nseas-1;

    for(s=1;s<=nseas;s++)
    {
      t = styr-3*nseas+s-1;
      Make_AgeLength_Key();
      Save_Wt_Age(t)=Wt_Age_beg(s);
    }

    for (s=1;s<=nseas;s++)
    for (g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
    {
      smid=s+nseas;
      Make_FishSelex();
    }

    if(SR_fxn==6 || SR_fxn==3)
    {
      alpha = 4.0 * SR_parm(2)*Recr_virgin / (5.*SR_parm(2)-1.);
      beta = (SPB_virgin*(1.-SR_parm(2))) / (5.*SR_parm(2)-1.);
    }

//  the spawner-recruitment function has Bzero based on virgin biology, not benchmark biology
//  need to deal with possibility that with time-varying biology, the SPB_virgin calculated from virgin conditions will differ from the SPB_virgin used for benchmark conditions

    //  NEED Replacement code to calc average recr_dist from the natage
//    get_MGsetup();    // in case recr_dist parameters have changed
//    get_recr_distribution();

// find Fspr             ProgLabel_710
    if(show_MSY==1)
    {
    report5<<"& & & & & Fspr search"<<endl;
    report5<<"Iter Fmult F_std SPR";
    for(f=1;f<=Nfleet;f++) {report5<<" catch:_"<<f;}
    for(p=1;p<=pop;p++)
    for(gp=1;gp<=N_GP;gp++)
    {report5<<" SSB_Area:"<<p<<"_GP:"<<gp;}
    report5<<endl;
    }
    Fmult=0.; Nloops=18; Closer=1.;
    F1(1)=log(1.0e-3); last_calc=0.; Fchange=-4.0;

    equ_Recr=1.0;

    Fishon=0;
    Do_Equil_Calc();
    SPR_unf=SPB_equil;  //  this corresponds to the biology for benchmark average years, not the virgin SPB_virgin
    Vbio1_unfished=smrybio;       // gets value from equil_calc
          if(show_MSY==1)
          {
          report5<<0<<" "<<0.00<<" "<<1.0;
          for(f=1;f<=Nfleet;f++) report5<<" "<<0.0;
          for(p=1;p<=pop;p++)
          for(gp=1;gp<=N_GP;gp++)
          {report5<<" "<<SPB_equil_pop_gp(p,gp);}
          report5<<endl;
          }

    df=1.e-5;
    Fishon=1;
    for(j=1;j<=Nloops;j++)   // loop find Fspr
    {
      if(fabs(Fchange)<=0.25)
        {
          jj=3;
          F1(2) = F1(1) + df*.5;
          F1(3) = F1(2) - df;
        }
      else
        {jj=1;}

      for(int ii=jj;ii>=1;ii--)
        {
          Fmult=mfexp(F1(ii));

          for(f=1;f<=Nfleet;f++)
          for(s=1;s<=nseas;s++)
            {t=bio_t_base+s; Hrate(f,t)=Fmult*Bmark_RelF_Use(s,f);}

          Fishon=1;
          Do_Equil_Calc();
          yld1(ii)=SPB_equil/SPR_unf;
        }
        SPR_actual=yld1(1);

          if(jj==3)
            {
            Closer*=0.5;
              dyld=(yld1(2) - yld1(3))/df;                      // First derivative (to find the root of this)
              if(dyld!=0.)
                {last_F1=F1(1); F1(1) += (SPR_target-SPR_actual)/dyld;
                 F1(1)=(1.-Closer)*F1(1)+Closer*last_F1;
                }        // averages with last good value to keep from changing too fast
              else
                {F1(1)=(F1(1)+last_F1)*0.5;}    // go halfway back towards previous value
            }
          else
            {
              if((last_calc-SPR_target)*(SPR_actual-SPR_target)<0.0) {Fchange*=-0.5;}   // changed sign, so reverse search direction
              F1(1)+=Fchange;  last_calc=SPR_actual;
            }

          if(show_MSY==1)
          {
            report5<<j<<" "<<Fmult<<" "<<equ_F_std<<" "<<SPR_actual<<" "<<column(equ_catch_fleet,2)*1.;
            for(p=1;p<=pop;p++)
            for(gp=1;gp<=N_GP;gp++)
            {report5<<" "<<SPB_equil_pop_gp(p,gp);}
            report5<<endl;
          }
    }   // end search loop

    if(fabs(SPR_actual-SPR_target)>=0.001)
    {N_warn++; warning<<" warning: poor convergence in Fspr search "<<SPR_target<<" "<<SPR_actual<<endl;}
    if(show_MSY==1)
    {
      report5<<"fleet encB deadB retB encN deadN retN): "<<endl;
      for(f=1;f<=Nfleet;f++) report5<<f<<" "<<equ_catch_fleet(f)<<endl;
    }

    SPR_temp=SPR_actual*SPR_unf;
    Get_EquilCalc = Equil_Spawn_Recr_Fxn();   // call  function

    Bspr=Get_EquilCalc(1);
    Bspr_rec=Get_EquilCalc(2);
    YPR_tgt_enc  = YPR_enc;         //  total encountered yield per recruit
    YPR_tgt_dead = YPR_dead;           // total dead yield per recruit
    YPR_tgt_N_dead = YPR_N_dead;
    YPR_tgt_ret = YPR_ret;  SPR_Fmult=Fmult;
    if(rundetail>0 && mceval_counter==0) cout<<" got Fspr "<<SPR_Fmult<<" "<<SPR_actual<<endl;
    YPR_spr=YPR_tgt_dead; Vbio_spr=totbio; Vbio1_spr=smrybio;
    Mgmt_quant(10)=equ_F_std;
    Mgmt_quant(9)=Get_EquilCalc(1);
    Mgmt_quant(11)=YPR_dead*Get_EquilCalc(2);

    SPR_at_target=SPR_actual;
//   end finding Fspr


// ******************************************************
//  find F giving Btarget      ProgLabel_720
    if(show_MSY==1)
    {
      report5<<"+ + + + + + + + + findBtarget"<<endl<<"Iter Fmult F_std SPR Catch SSB Recruits SSB/Bzero ";
      for(f=1;f<=Nfleet;f++) {report5<<" catch:_"<<f;}
      for(p=1;p<=pop;p++)
      for(gp=1;gp<=N_GP;gp++)
      {report5<<" SSB_Area:"<<p<<"_GP:"<<gp;}
      report5<<endl;
    }

    F1(1)=log(1.0e-3); last_calc=0.; Fchange=-4.0; df=1.e-5; Closer=1.;
    dvariable Closer2;
    if(SR_fxn==5) {Closer2=0.001; Nloops=40;} else {Closer2=0.10; Nloops=28;}

    Btgttgt=BTGT_target*SPB_virgin;   //  this is relative to virgin, not to the average biology from benchmark years
    for(j=0;j<=Nloops;j++)   // loop find Btarget
      {
      if(fabs(Fchange)<=Closer2)
        {
        jj=3;
        F1(2) = F1(1) + df*.5;
        F1(3) = F1(2) - df;
        }
      else
        {jj=1;}
      for(int ii=jj;ii>=1;ii--)
      {
        if(j==0) {Fmult=0.0;} else {Fmult=mfexp(F1(ii));}
        for(f=1;f<=Nfleet;f++)
        for(s=1;s<=nseas;s++)
        {
          t=bio_t_base+s; Hrate(f,t)=Fmult*Bmark_RelF_Use(s,f);
        }
        Do_Equil_Calc();
        SPR_Btgt = SPB_equil/SPR_unf;   //  here for SPR it uses benchmark's SPB_virgin for consistency
        SPR_temp=SPB_equil;
        Get_EquilCalc = Equil_Spawn_Recr_Fxn();   // call  function
        yld1(ii)=Get_EquilCalc(1);
      }

      Btgt=Get_EquilCalc(1);  //  so uses benchmark average years

      if(jj==3)
        {
        Closer *=0.5;
        dyld=(yld1(2) - yld1(3))/df;                      // First derivative
        if(dyld!=0.)
          {last_F1=F1(1); F1(1) -= (Btgt-Btgttgt)/dyld;
           F1(1)=(1.-Closer)*F1(1)+(Closer)*last_F1;
          }        // weighted average with last good value to keep from changing too fast
        else
          {F1(1)=(F1(1)+last_F1)*0.5;}    // go halfway back towards previous value
        }
      else
        {
          temp=(last_calc-Btgttgt)*(Btgt-Btgttgt)/(sfabs(last_calc-Btgttgt)*sfabs(Btgt-Btgttgt));  // values of -1 or 1
          temp1=temp-1.;  // values of -2 or 0
          Fchange*=exp(temp1/4.)*temp;
          F1(1)+=Fchange;  last_calc=Btgt;
        }

      if(show_MSY==1)
      {
        report5<<j<<" "<<Fmult<<" "<<equ_F_std<<" "<<Btgt/Get_EquilCalc(2)<<" "<<YPR_dead*Get_EquilCalc(2)<<" "<<Btgt<<" "<<Get_EquilCalc(2)
        <<" "<<Btgt/SPB_virgin<<" "<<column(equ_catch_fleet,2)*Get_EquilCalc(2);
        for(p=1;p<=pop;p++)
        for(gp=1;gp<=N_GP;gp++)
        {report5<<" "<<SPB_equil_pop_gp(p,gp)*Get_EquilCalc(2);}
        report5<<endl;
      }

      }   // end search loop

    Btgt_Rec=Get_EquilCalc(2);
    if(fabs(log(Btgt/Btgttgt))>=0.001)
    {N_warn++; warning<<" warning: poor convergence in Btarget search "<<Btgttgt<<" "<<Btgt<<endl;}
    if(show_MSY==1)
    {
      report5<<"fleet encB deadB retB encN deadN retN): "<<endl;
      for(f=1;f<=Nfleet;f++) report5<<f<<" "<<Btgt_Rec*equ_catch_fleet(f)<<endl;
    }
    Btgt_Fmult=Fmult;
    if(rundetail>0 && mceval_counter==0) cout<<" got_Btgt "<<Btgt_Fmult<<" "<<Btgt/SPB_virgin<<endl;
    YPR_Btgt_enc  = YPR_enc;         //  total encountered yield per recruit
    YPR_Btgt_dead = YPR_dead;           // total dead yield per recruit
    YPR_Btgt_N_dead = YPR_N_dead;           // total dead yield per recruit
    YPR_Btgt_ret = YPR_ret;
    Vbio_Btgt=totbio; Vbio1_Btgt=smrybio;
    Mgmt_quant(7)=equ_F_std;
    Mgmt_quant(5)=Btgt;
    Mgmt_quant(6)=SPR_Btgt;
    Mgmt_quant(8)=YPR_dead*Btgt_Rec;

//  end finding F for Btarget


// ******************************************************
//  start finding Fmsy     ProgLabel_730
//  consider using maxpossF here, instead of calculating a new Fmax

    if(Do_MSY==0)
      {
       Fmax=1.; MSY=-1; Bmsy=-1; Recr_msy=-1; MSY_SPR=-1; Yield=-1; totbio=1; smrybio=1.; MSY_Fmult=-1.;   //  use these values if MSY is not calculated
       Mgmt_quant(1)=SPB_virgin;  // this may be redundant
       if(show_MSY==1) report5<<"MSY_not_calculated;_ignore_values"<<endl;
      }
    else
    {
      if(F_Method>=2) {Fmax=maxpossF/sum(Bmark_RelF_Use);}
      switch(Do_MSY)
        {
        case 1:  // set Fmsy=Fspr
          {Fmult=SPR_Fmult;
           if(F_Method==1) {Fmax=SPR_Fmult*1.1;}
           F1(1)=-log(Fmax/SPR_Fmult-1.); last_calc=0.; Fchange=1.0; Closer=1.; Nloops=0;
           break;}
        case 2:  // calc Fmsy
          {last_calc=0.; Fchange=0.51; Closer=1.0;
           if(SR_fxn==5) {Nloops=40;} else {Nloops=19;}
          if(F_Method==1) {Fmax=(Btgt_Fmult+SPR_Fmult)*0.5*SR_parm(2)/0.05;}    //  previously /0.18
           F1(1)=-log(Fmax/Btgt_Fmult-1.);
          break;}
        case 3:  // set Fmsy=Fbtgt
          {Fmult=Btgt_Fmult;
           if(F_Method==1) {Fmax=Btgt_Fmult*1.1;}
            F1(1)=-log(Fmax/Btgt_Fmult-1.); last_calc=0.; Fchange=1.0; Closer=1.0; Nloops=0;
          break;}
        case 4:   //  set fmult for Fmsy to 1
          {Fmult=1; Fmax=1.1; F1(1)=-log(Fmax/Fmult-1.); last_calc=0.; Fchange=1.0; Closer=1.0; Nloops=0;
          break;}
        }

      if(show_MSY==1)
      {
        report5<<"+ + + + + + + + + findFmsy"<<endl<<"Iter Fmult F_std SPR Catch SSB Recruits SSB/Bzero Gradient Curvature";
        for(f=1;f<=Nfleet;f++) {report5<<" catch:_"<<f;}
        for(p=1;p<=pop;p++)
        for(gp=1;gp<=N_GP;gp++)
        {report5<<" Area:"<<p<<"_GP:"<<gp;}
        report5<<endl;
      }

        bestF1.initialize(); bestF2.initialize();

      df=0.050;
      jj=3;
      Fishon=1;
      for(j=0;j<=Nloops;j++)   // loop to find Fmsy
        {
         df*=.95;
        Closer*=0.8;
          F1(2) = F1(1) + df*.5;
          F1(3) = F1(2) - df;
        for(int ii=jj;ii>=1;ii--)
          {
          Fmult=Fmax/(1+mfexp(-F1(ii)));
          for(f=1;f<=Nfleet;f++)
          for(s=1;s<=nseas;s++)
            {t=bio_t_base+s; Hrate(f,t)=Fmult*Bmark_RelF_Use(s,f);}

          Do_Equil_Calc();
          MSY_SPR = SPB_equil/SPR_unf;
          SPR_temp=SPB_equil;
          Get_EquilCalc = Equil_Spawn_Recr_Fxn();   // call  function
          Bmsy=Get_EquilCalc(1);
          Recr_msy=Get_EquilCalc(2);
          yld1(ii)=YPR_dead*Recr_msy;   //  *mfexp(-Equ_penalty);
          Yield=YPR_dead*Recr_msy;
          bestF1+=F1(ii)*(pow(mfexp(Yield/1.0e08),5)-1.);
          bestF2+=pow(mfexp(Yield/1.0e08),5)-1.;
          }   //  end gradient calc

        dyld   = (yld1(2) - yld1(3))/df;                      // First derivative (to find the root of this)
        temp  = (yld1(2) + yld1(3) - 2.*yld1(1))/(.25*df*df);   // Second derivative (for Newton Raphson)
        dyldp = -sqrt(temp*temp+1.);   //  add 1 to keep curvature reasonably large
        last_F1=F1(1);
        temp = F1(1)-dyld*(1.-Closer)/(dyldp);
        if(show_MSY==1)
        {
          report5<<j<<" "<<Fmult<<" "<<equ_F_std<<" "<<MSY_SPR<<" "<<yld1(1)<<" "<<Bmsy<<" "<<Recr_msy<<" "<<Bmsy/SPB_virgin<<" "
          <<dyld <<" "<<dyldp<<" "<<value(column(equ_catch_fleet,2)*Recr_msy);
          for(p=1;p<=pop;p++)
          for(gp=1;gp<=N_GP;gp++)
          {report5<<" "<<SPB_equil_pop_gp(p,gp)*Recr_msy;}
          report5<<endl;
        }

        if(j<=9)
          {F1(1)=(1.-Closer)*temp+Closer*(bestF1/bestF2);}        // averages with best value to keep from changing too fast
        else
          {F1(1)=temp;}
        }   // end search loop
    if(fabs(dyld/dyldp)>=0.001 && Do_MSY==2)
    {N_warn++; warning<<" warning: poor convergence in Fmsy, final dy/dy2= "<<dyld/dyldp<<endl;}

      YPR_msy_enc = YPR_enc;
      YPR_msy_dead = YPR_dead;           // total dead yieldt
      YPR_msy_N_dead = YPR_N_dead;           // total dead yield
      YPR_msy_ret = YPR_ret;           // total retained yield
      MSY=Yield;
      MSY_Fmult=Fmult;
      Mgmt_quant(15)=Yield;
      Mgmt_quant(12)=Bmsy;
      Mgmt_quant(13)=MSY_SPR;
      Mgmt_quant(14)=equ_F_std;
      Mgmt_quant(16)=YPR_ret*Recr_msy;
      Vbio1_MSY=smrybio;
      Vbio_MSY=totbio;

      if(show_MSY==1)
      {
    report5<<"fleet encB deadB retB encN deadN retN): "<<endl;
    for(f=1;f<=Nfleet;f++) report5<<f<<" "<<Recr_msy*equ_catch_fleet(f)<<endl;
    report5<<"Equil_N_at_age_at_MSY_each"<<endl<<"Seas Area GP Sex subM"<<age_vector<<endl;
     for(s=1;s<=nseas;s++)
     for(p=1;p<=pop;p++)
     for(g=1;g<=gmorph;g++)
     {if(use_morph(g)>0) report5<<s<<" "<<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<GP2(g)<<" "<<Recr_msy*equ_numbers(s,p,g)(0,nages)<<endl;}

    report5<<"Equil_N_at_age_at_MSY_sum"<<endl<<"GP Sex N/Z"<<age_vector<<endl;
    for(gg=1;gg<=gender;gg++)
    for(gp=1;gp<=N_GP;gp++)
    {
      tempvec_a.initialize();
      for (p=1;p<=pop;p++)
      for (g=1;g<=gmorph;g++)
      if(use_morph(g)>0)
      {
        if(GP4(g)==gp && sx(g)==gg) tempvec_a+= value(Recr_msy*equ_numbers(1,p,g)(0,nages));
      }
      if(nseas>1)
      {
        tempvec_a(0)=0.;
        for(s=1;s<=nseas;s++)
        for (p=1;p<=pop;p++)
        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0 && Bseas(g)==s)
        {
          if(GP4(g)==gp && sx(g)==gg) tempvec_a(0) += value(Recr_msy*equ_numbers(1,p,g,0));
        }
      }
      report5 <<gp<<" "<<gg<<" N "<<tempvec_a<<endl;
      report5 <<gp<<" "<<gg<<" Z ";
      for(a=0;a<=nages-2;a++)
      {report5<<-log(tempvec_a(a+1)/tempvec_a(a))<<" ";}
      report5<<" NA NA"<<endl;
    }

     Fishon=0;
     Do_Equil_Calc();
    report5<<"Equil_N_at_age_M_only_Recr_MSY"<<endl<<"Seas Area GP Sex subM"<<age_vector<<endl;
     for(s=1;s<=nseas;s++)
     for(p=1;p<=pop;p++)
     for(g=1;g<=gmorph;g++)
     {if(use_morph(g)>0) report5<<s<<" "<<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<GP2(g)<<" "<<Recr_msy*equ_numbers(s,p,g)(0,nages)<<endl;}

    report5<<"Equil_N_at_age_M_only_sum"<<endl<<"GP Sex N/Z "<<age_vector<<endl;
    for(gg=1;gg<=gender;gg++)
    for(gp=1;gp<=N_GP;gp++)
    {
      tempvec_a.initialize();
      for (p=1;p<=pop;p++)
      for (g=1;g<=gmorph;g++)
      if(use_morph(g)>0)
      {
        if(GP4(g)==gp && sx(g)==gg) tempvec_a+= value(Recr_msy*equ_numbers(1,p,g)(0,nages));
      }
      if(nseas>1)
      {
        tempvec_a(0)=0.;
        for(s=1;s<=nseas;s++)
        for (p=1;p<=pop;p++)
        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0 && Bseas(g)==s)
        {
          if(GP4(g)==gp && sx(g)==gg) tempvec_a(0) += value(Recr_msy*equ_numbers(1,p,g,0));
        }
      }
      report5 <<gp<<" "<<gg<<" N "<<tempvec_a<<endl;
      report5 <<gp<<" "<<gg<<" Z ";
      for(a=0;a<=nages-2;a++)
      {report5<<-log(tempvec_a(a+1)/tempvec_a(a))<<" ";
        }
      report5<<" NA NA"<<endl;
    }

     Fishon=1;

    if(Fmult*3.0 <= SPR_Fmult) {N_warn++; warning<<" Fmsy is <1/3 of Fspr are you sure?  check for convergence "<<endl;}
    if(Fmult/3.0 >= SPR_Fmult) {N_warn++; warning<<" Fmsy is >3x of Fspr are you sure?  check for convergence "<<endl;}
    if(Fmult/0.98 >= Fmax) {N_warn++; warning<<" Fmsy is close to max allowed; check for convergence "<<endl;}
      }
    }

    if(rundetail>0 && mceval_counter==0) cout<<" got Fmsy "<<MSY_Fmult<<" "<<MSY<<endl;

// ***************** show management report   ProgLabel_740
    if(show_MSY==1)
      {
    report5<<"+ + + + +"<<endl<<"Management_report"<<endl;
    report5<<"Steepness_Recr_SPB_virgin "<<SR_parm(2)<<" "<<Recr_virgin<<" "<<SPB_virgin<<endl;
    report5<<"+"<<endl<<"Element Value Bio/Recr Bio/R0 Numbers N/R0 (B_in_mT;_N_in_thousands)"<<endl;
    report5<<"Recr_unfished "<<Recr_virgin<<" -- -- "<<endl;
    report5<<"SPR_unfished "<<SPR_unf*Recr_virgin<<" "<<SPR_unf<<endl;
    report5<<"BIO_Smry_unfished "<<Vbio1_unfished*Recr_virgin<<" "<<Vbio1_unfished<<" "<<Vbio1_unfished<<endl<<"+ + + + +"<<endl;

    report5<<"SPR_target "<<SPR_target<<endl;
    report5<<"SPR_calc "<<SPR_actual<<endl;
    report5<<"Fmult "<<SPR_Fmult<<endl;
    report5<<"F_std "<<Mgmt_quant(10)<<endl;
    report5<<"Exploit(Y/Bsmry) "<<YPR_spr/Vbio1_spr<<endl;
    report5<<"Recruits@Fspr "<<Bspr_rec<<" -- -- "<<Bspr_rec<<" "<<Bspr_rec/Recr_virgin<<" "<<endl;
    report5<<"SPBio "<<SPR_at_target*Bspr_rec*SPR_unf<<" "<<SPR_at_target*SPR_unf<<" -- "<<endl;
    report5<<"YPR_encountered "<<YPR_tgt_enc*Bspr_rec<<" "<<YPR_tgt_enc<<" -- "<<endl;
    report5<<"YPR_dead "<<YPR_tgt_dead*Bspr_rec<<" "<<YPR_tgt_dead<<" -- "<<" "<<YPR_tgt_N_dead*Bspr_rec<<endl;
    report5<<"YPR_retain "<<YPR_tgt_ret*Bspr_rec<<" "<<YPR_tgt_ret<<" -- "<<endl;
    report5<<"Biomass_Smry "<<Vbio1_spr*Bspr_rec<<" "<<Vbio1_spr<<" -- "<<endl<<"+ + + + +"<<endl;

    report5<<"Btarget  "<<BTGT_target<<endl;
    report5<<"Btgt_calc_rel_SPB_virgin "<<Btgt/SPB_virgin<<endl;
    report5<<"Btgt_calc_rel_SPB_benchmark "<<Btgt/(SPR_unf*Recr_virgin)<<endl;
    report5<<"SPR_for_Btgt "<<SPR_Btgt<<endl;
    report5<<"Fmult "<<Btgt_Fmult<<endl;
    report5<<"F_std "<<Mgmt_quant(7)<<endl;
    report5<<"Exploit(Y/Bsmry) "<<YPR_Btgt_dead/Vbio1_Btgt<<endl;
    report5<<"Recruits@Btgt "<<Btgt_Rec<<" -- -- "<<Btgt_Rec<<" "<<Btgt_Rec/Recr_virgin<<endl;
    report5<<"SPBio "<<Btgt<<" "<<Btgt/Btgt_Rec<<" -- "<<endl;
    report5<<"YPR_encountered "<<YPR_Btgt_enc*Btgt_Rec<<" "<<YPR_Btgt_enc<<" -- "<<endl;
    report5<<"YPR_dead "<<YPR_Btgt_dead*Btgt_Rec<<" "<<YPR_Btgt_dead<<" -- "<<YPR_Btgt_N_dead*Btgt_Rec<<endl;
    report5<<"YPR_retain "<<YPR_Btgt_ret*Btgt_Rec<<" "<<YPR_Btgt_ret<<" -- "<<endl;
    report5<<"Biomass_Smry "<<Vbio1_Btgt*Btgt_Rec<<" "<<Vbio1_Btgt<<" -- "<<endl<<"+ + + + +"<<endl;

        switch(Do_MSY)
          {
          case 1:  // set Fmsy=Fspr
            {report5<<"set_Fmsy=Fspr"<<endl;
            break;}
          case 2:  // calc Fmsy
            {report5<<"calculate_FMSY"<<endl;
            break;}
          case 3:  // set Fmsy=Fbtgt
            {report5<<"set_Fmsy=Fbtgt"<<endl;
            break;}
          case 4:   //  set fmult for Fmsy to 1
            {report5<<"set_Fmsy_using_Fmult=1.0"<<endl;
            break;}
          }
    report5<<"SPR "<<MSY_SPR<<endl;
    report5<<"Fmult "<<MSY_Fmult<<endl;
    report5<<"F_std "<<Mgmt_quant(14)<<endl;
    report5<<"Exploit(Y/Bsmry) "<<MSY/(Vbio1_MSY*Recr_msy)<<endl;
    report5<<"Recruits@MSY "<<Recr_msy<<" -- -- "<<Recr_msy<<" "<<Recr_msy/Recr_virgin<<endl;
    report5<<"SPBio "<<Bmsy<<" "<<Bmsy/Recr_msy<<" -- "<<endl;
    report5<<"SPBmsy/SPBzero(using_SPB_virgin) "<<Bmsy/SPB_virgin<<" -- --"<<endl;  // new version
    report5<<"SPBmsy/SPBzero(using_BenchmarkYr_biology) "<<Bmsy/(Recr_virgin*SPR_unf)<<" -- --"<<endl;
    report5<<"MSY_for_optimize "<<MSY<<" "<<MSY/Recr_msy<<" -- "<<endl;
    report5<<"MSY_encountered "<<YPR_msy_enc*Recr_msy<<" "<<YPR_msy_enc<<" -- "<<endl;
    report5<<"MSY_dead "<<YPR_msy_dead*Recr_msy<<" "<<YPR_msy_dead<<" -- "<<YPR_msy_N_dead*Recr_msy<<endl;
    report5<<"MSY_retain "<<YPR_msy_ret*Recr_msy<<" "<<YPR_msy_ret<<" -- "<<endl;
    report5<<"Biomass_Smry "<<Vbio1_MSY*Recr_msy<<" "<<Vbio1_MSY<<" -- "<<endl<<"+"<<endl;
    report5<<"Summary_age: "<<Smry_Age<<endl<<"#"<<endl;
    report5<<"#_note when there is time-varying biology"<<endl;
    report5<<"#_virgin outputs use biology at begin of time series"<<endl;
    report5<<"#_SPR outputs use biology averaged over years: "<<Bmark_Yr(1)<<" "<<Bmark_Yr(2)<<endl;
    report5<<"#_Btgt uses Bmark years biology to search for fraction of SPB_virgin, so take care in interpretation "<<endl;
    report5<<"#_MSY and Bmsy use Bmark years biology"<<endl<<"#"<<endl;
    if(F_Method==1)
    {
      report5<<"F_reported_below_is_Pope's_midseason_exploitation_rate=MSY_Fmult*Alloc"<<endl;
      report5<<"seas seas_dur "; for (f=1;f<=Nfleet;f++) {report5<<" fleet:"<<f;}
      report5<<endl;
      for(s=1;s<=nseas;s++) {report5<<s<<" "<<seasdur(s)<<" "<<MSY_Fmult*Bmark_RelF_Use(s)<<endl;}
    }
    else
    {
      report5<<"F_reported_here_is_Seasonal_apicalF=MSY_Fmult*Alloc*seas_dur_(can_be>F_std_because_of_selex)"<<endl;
      report5<<"seas seas_dur "; for (f=1;f<=Nfleet;f++) {report5<<" fleet:"<<f;}
      report5<<endl;
      for(s=1;s<=nseas;s++) {report5<<s<<" "<<seasdur(s)<<" "<<MSY_Fmult*Bmark_RelF_Use(s)*seasdur(s)<<endl;}
    }
    report5<<"#"<<endl;
      }
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 35 Get_Forecast */
FUNCTION void Get_Forecast()
  {
  if(rundetail>0 && mceval_counter==0) cout<<" Do Forecast "<<YrMax<<endl;
  t_base=styr+(endyr-styr)*nseas-1;
  int Do_4010;
  int bio_t;
  int adv_age;
  dvariable OFL_catch;
  dvariable Fcast_Crash;
  dvariable totcatch;
  dvar_matrix catage_w(1,gmorph,0,nages);
  dvar_vector tempcatch(1,Nfleet);
  dvar_vector ABC_buffer(endyr+1,YrMax);
  imatrix Do_F_tune(t_base,TimeMax_Fcast_std,1,Nfleet);  //  flag for doing F from catch
  dvar_matrix Fcast_Catch_Store(t_base,TimeMax_Fcast_std,1,Nfleet);
  dvar_vector Fcast_Catch_Current_Alloc(1,Nfleet);
  dvar_vector Fcast_Catch_Current_Allocation(1,Fcast_Catch_Allocation_Groups);
  dvar_vector Fcast_Catch_ByArea(1,pop);

    dvar_vector  H_temp(1,Nfleet);
    dvar_vector  C_temp(1,Nfleet);
    dvar_vector  H_old(1,Nfleet);
    dvar_vector  C_old(1,Nfleet);

  int Tune_F;
  int Tune_F_loops;

  int ABC_Loop_start;
  int ABC_Loop_end;

  Do_F_tune.initialize();

   report5<<"Forecast_Option:_"<<Do_Forecast<<endl;
   switch(Do_Forecast)
   {
     case 1:
       {Fcast_Fmult=SPR_Fmult; if(show_MSY==1) report5<<"Forecast_using_Fspr"<<endl; break;}
     case 2:
       {Fcast_Fmult=MSY_Fmult; if(show_MSY==1) report5<<"Forecast_using_Fmsy"<<endl; break;}
     case 3:
       {Fcast_Fmult=Btgt_Fmult; if(show_MSY==1) report5<<"Forecast_using_F(Btarget)"<<endl; break;}
     case 4:
     {
       Fcast_Fmult=0.0;
       for(y=Fcast_RelF_yr1;y<=Fcast_RelF_yr2;y++)
       for(f=1;f<=Nfleet;f++)
       for(s=1;s<=nseas;s++)
       {
         t=styr+(y-styr)*nseas+s-1;
         Fcast_Fmult+=Hrate(f,t);
       }
       Fcast_Fmult/=float(Fcast_RelF_yr2-Fcast_RelF_yr1+1);
       Fcurr_Fmult=Fcast_Fmult;
       if(show_MSY==1) report5<<"Forecast_using_ave_F_from:_"<<Fcast_RelF_yr1<<"_"<<Fcast_RelF_yr2<<endl; break;
     }
     case 5:
     {Fcast_Fmult=Fcast_Flevel; if(show_MSY==1) report5<<"Forecast_using_input_F "<<endl; break;}
   }

  if(show_MSY==1)  //  write more headers
  {
    report5<<"Annual_Forecast_Fmult: "<<Fcast_Fmult<<endl;
    report5<<"Fmultiplier_during_selected_relF_years_was: "<<Fcurr_Fmult<<endl;
    report5<<"Selectivity_averaged_over_yrs:_"<<Fcast_Sel_yr1<<"_to_"<<Fcast_Sel_yr2<<endl;
    report5<<"Cap_totalcatch_by_fleet "<<endl<<Fcast_MaxFleetCatch<<endl;
    report5<<"Cap_totalcatch_by_area "<<endl<<Max_Fcast_Catch<<endl;
    report5<<"Assign_fleets_to_allocation_groups_(0_means_not_in_a_group) "<<endl<<Allocation_Fleet_Assignments<<endl;
    report5<<"Calculated_number_of_allocation_groups "<<Fcast_Catch_Allocation_Groups<<endl;
    report5<<"Allocation_among_groups "<<Fcast_Catch_Allocation<<endl;
    if(Fcast_Catch_Basis==2)
    {report5<<"2_:_Caps_&_Alloc_use_dead_catchbio"<<endl;}
    else if(Fcast_Catch_Basis==3)
    {report5<<"3_:_Caps_&_Alloc_use_retained_catchbio"<<endl;}
    else if(Fcast_Catch_Basis==5)
    {report5<<"5_:_Caps_&_Alloc_use_dead_catchnum"<<endl;}
    else if(Fcast_Catch_Basis==6)
    {report5<<"6_:_Caps_&_Alloc_use_retained_catchnum"<<endl;}
    if(N_Fcast_Input_Catches>0)
    {
      report5<<"-1_:_Input_fixed_catches_or_F_with_fleet/time_specific_values"<<endl;
    }
    report5<<"#_Relative_F_among_fleets"<<endl;
    if(Fcast_RelF_Basis==1)
    {
      report5<<"based_on_years:_"<<Fcast_RelF_yr1<<"_to_"<<Fcast_RelF_yr2<<endl;
    }
    else
    {
      report5<<"read_from_input_file"<<endl;
    }
    if(F_Method==1)
    {
      report5<<"Pope's_midseason_exploitation_rate=Fmult*Alloc"<<endl;
      report5<<"seas seas_dur "; for (f=1;f<=Nfleet;f++) {report5<<" fleet:"<<f;}
      report5<<endl;
      for(s=1;s<=nseas;s++) {report5<<s<<" "<<seasdur(s)<<" "<<Fcast_Fmult*Fcast_RelF_Use(s)<<endl;}
    }
    else
    {
      report5<<"Seasonal_apicalF=Fmult*Alloc*seas_dur_(can_be>F_std_because_of_selex)"<<endl;
      report5<<"seas seas_dur "; for (f=1;f<=Nfleet;f++) {report5<<" fleet:"<<f;}
      report5<<endl;
      for(s=1;s<=nseas;s++) {report5<<s<<" "<<seasdur(s)<<" "<<Fcast_Fmult*Fcast_RelF_Use(s)*seasdur(s)<<endl;}
    }
    report5<<"#"<<endl;
    report5<<"N_forecast_yrs: "<<N_Fcast_Yrs<<endl;
    report5<<"OY_Control_Rule "<<" Inflection: "<<H4010_top<<" Intercept: "<<H4010_bot<<" Scale: "<<H4010_scale;
    if(HarvestPolicy==1) {report5<<" adjust_catch_below_Inflection(west_coast)"<<endl;} else {report5<<" adjust_F_below_Inflection"<<endl;}
    report5<<"#"<<endl;
  }

  for(int Fcast_Loop1=1; Fcast_Loop1<=Fcast_Loop_Control(1);Fcast_Loop1++)  //   for different forecast conditions
  {
    switch(Fcast_Loop1)  //  select which ABC_loops to use
    {
      case 1:  // do OFL only
      {
        ABC_Loop_start=1;
        ABC_Loop_end=1;
        if(show_MSY==1) report5<<"FORECAST:_With_Constant_F=Fofl;_No_Input_Catches_or_Adjustments;_Equil_Recr;_No_inpl_error"<<endl;
        break;
      }
      case 2:  //  for each year:  do 3 calculations:  (1) OFL, (2) calc ABC and apply caps and allocations, (3) get F from catch _impl
      {
        ABC_Loop_start=1;
        ABC_Loop_end=3;
        if(show_MSY==1) report5<<"FORECAST:_With_F=Fabc;_With_Input_Catches_and_Catch_Adjustments;_Equil_Recr;_No_inpl_error"<<endl;
        break;
      }
      case 3:  //  just need to get F from stored adjusted catch (after modifying stored catch by implementation error).
      {
        ABC_Loop_start=3;
        ABC_Loop_end=3;
        if(show_MSY==1) report5<<"FORECAST:_With_F_to_match_adjusted_catch;_With_Input_Catches_and_Catch_Adjustments;_Stochastic_Recr;_With_inpl_error"<<endl;
        break;
      }
    }
    if(show_MSY==1)
    {
    report5<<"pop year ABC_Loop season Ctrl_Rule bio-all bio-Smry SpawnBio Depletion recruit-0 ";
    for (f=1;f<=Nfleet;f++) {report5<<" sel(B):_"<<f<<" dead(B):_"<<f<<" retain(B):_"<<f<<
    " sel(N):_"<<f<<" dead(N):_"<<f<<" retain(N):_"<<f<<" F:_"<<f<<" R/C";}
    report5<<" Catch_Cap Total_Catch F_Std"<<endl;
    }

    //  note that spawnbio and Recruits need to retain their value from calculation in endyr,
    //  so can be used to distribute recruitment in year endyr+1 if recruitment distribution occurs before spawning season
    //  would be better to back up to last mainrecrdev and start with begin of forecast
    SPB_current=SPB_yr(endyr);
    Recruits=exp_rec(endyr,4);

    for(y=endyr+1;y<=YrMax;y++)
    {
      t_base=styr+(y-styr)*nseas-1;

      Smry_Table(y).initialize();

      if(Fcast_Loop1==3 && Do_Impl_Error>0)  //  apply implementation error, which is a random variable, so adds variance to forecast
                                             //  in future, could do this a fleet-specific implementation error
      {
        for(s=1;s<=nseas;s++)
        {
          t=t_base+s;
          for(f=1;f<=Nfleet;f++)
          {Fcast_Catch_Store(t,f)*=mfexp(Fcast_impl_error(y));}  //  should this be bias adjusted?
        }
      }

//  do biology for this year
      yz=endyr+1;  //  biology year for parameters
      if(time_vary_MG(endyr+1,2)>0 || save_for_report>0)  //  so uses endyr+1 timevary setting for duration of forecast
      {
        get_MGsetup();
        get_growth2();
      }

      // ABC_loop:  1=get OFL; 2=get_ABC, use input catches; 3=recalc with caps and allocations
      for(int ABC_Loop=ABC_Loop_start; ABC_Loop<=ABC_Loop_end;ABC_Loop++)
      {
        totcatch=0.;
        if(ABC_Loop==1) Mgmt_quant(Fcast_catch_start+N_Fcast_Yrs+y-endyr)=0.0;   // for OFL
        Mgmt_quant(Fcast_catch_start+y-endyr)=0.0;  //  for ABC
        if(Do_Retain==1) Mgmt_quant(Fcast_catch_start+2*N_Fcast_Yrs+y-endyr)=0.0;  // for retained ABC
        if(STD_Yr_Reverse_F(y)>0) F_std(STD_Yr_Reverse_F(y))=0.0;

        for(s=1;s<=nseas;s++)
        {
          t = t_base+s;
          smid=s+nseas;
          bio_t=styr+(endyr-styr)*nseas+s-1;

          if(ABC_Loop==ABC_Loop_start)  // do seasonal ALK and fishery selex
          {
            if(time_vary_MG(endyr+1,2)>0 || time_vary_MG(endyr+1,3)>0 || WTage_rd>0)
            {
              Make_AgeLength_Key();        // indexed by season and gmorph
            }
            else
            {
              Ave_Size(t)=Ave_Size(t-nseas);
              Save_Wt_Age(t)=Wt_Age_beg(s);
            }

            for (g=1;g<=gmorph;g++)
            if(use_morph(g)>0)
            {
              Make_FishSelex();   // calcs fishery selex by current season, all fleets, current gmorph
            }
          }  //  end of seasonal biology

          if(s==nseas) {adv_age=1;} else {adv_age=0;}   //      advance age or not when doing survivorship

          if(s==spawn_seas)    //  get spawnbio in a forecast year
          {
//  NEED to adjust for timing within season!!!
            SPB_pop_gp(y).initialize();
            SPB_B_yr(y).initialize();
            SPB_N_yr(y).initialize();
            for (p=1;p<=pop;p++)
            {
              for (g=1;g<=gmorph;g++)
              if(sx(g)==1 && use_morph(g)>0)     //  female
              {
                SPB_pop_gp(y,p,GP4(g))+=fec(g)*natage(t,p,g);
                SPB_B_yr(y) += make_mature_bio(GP4(g))*natage(t,p,g);
                SPB_N_yr(y) += make_mature_numbers(GP4(g))*natage(t,p,g);
              }
            }
            SPB_current=sum(SPB_pop_gp(y));
            SPB_yr(y)=SPB_current;

            if(Hermaphro_Option>0)  // get male biomass
            {
              temp=0.0;
              MaleSPB(y).initialize();
              for (p=1;p<=pop;p++)
              {
                for (g=1;g<=gmorph;g++)
                if(sx(g)==2 && use_morph(g)>0)     //  male; all assumed to be mature
                {
                  MaleSPB(y,p,GP4(g)) += Save_Wt_Age(t,g)*natage(t,p,g);   // accumulates SSB by area and by growthpattern
                  temp += sum(natage(t,p,g));
                }
              }
              if(Hermaphro_maleSPB==1)  // add MaleSPB to female SSB
              {
                SPB_current+=sum(MaleSPB(y));
                SPB_yr(y)=SPB_current;
                SPB_B_yr(y) += sum(MaleSPB(y));
                SPB_N_yr(y) += temp;
              }
            }
            Recruits=Spawn_Recr(SPB_current);    //  recruitment with deviations
            if(Fcast_Loop1<Fcast_Loop_Control(2))    //  use expected recruitment  this should include environ effect - CHECK THIS
            {
              Recruits=exp_rec(y,2);
              exp_rec(y,4)=exp_rec(y,2);  // within the spawn_recr function this has value with recrdev, so need to reset here
            }
          }  //  end of spawner-recruitment calculations

          if(ABC_Loop==1)  //  doing OFL this loop
          {
            ABC_buffer(y)=1.0;
          }
          else if(ABC_Loop==2 && s==1)  // Calc the buffer in season 1, will use last year's spawnbio if multiseas and spawnseas !=1
          {
            join1=1./(1.+mfexp(10.*(SPB_current-H4010_bot*SPB_virgin)));
            join2=1./(1.+mfexp(10.*(SPB_current-H4010_top*SPB_virgin)));

            if(HarvestPolicy==1)  // west coast
            {
              ABC_buffer(y) = H4010_scale*
              (
              (0.0001*SPB_current/(H4010_bot*SPB_virgin) ) *(join1)   // low
              +(0.0001+(1.0-0.0001)*(H4010_top*SPB_virgin/SPB_current)*(SPB_current-H4010_bot*SPB_virgin)/(H4010_top*SPB_virgin-H4010_bot*SPB_virgin)) * (1.0-join1) // curve
              )
              *(join2)   // scale combo
              +
              (H4010_scale) * (1.0-join2);    // scale right side
            }
            else if(HarvestPolicy==2)  // Alaska
            {
              ABC_buffer(y) = H4010_scale*
              (
              (0.0001*SPB_current/(H4010_bot*SPB_virgin) ) *(join1)   // low
              +(0.0001+(1.0-0.0001)*(SPB_current-H4010_bot*SPB_virgin)/(H4010_top*SPB_virgin-H4010_bot*SPB_virgin)) * (1.0-join1)   // curve
              )
              *(join2)   // scale combo
              +
              (H4010_scale) * (1.0-join2);    // scale right side
            }
            else
            {
              ABC_buffer(y)=H4010_scale;
            }
          }  // end calc of ABC buffer
          else
          {  //  ABC buffer remains at previously calculated value
          }
          for(p=1;p<=pop;p++)  //  loop areas (sub-pops)
          {
            totbio.initialize();smrybio.initialize(); smrynum.initialize();
            for (g=1;g<=gmorph;g++)
            if(use_morph(g)>0)
            {
              gg=sx(g);
//              if(s==Bseas(g)) {natage(t,p,g,0) = Recruits*recr_dist_calc(GP(g),s,p)*submorphdist(GP2(g));} //get recruits from total recruits and RecrDist
              if(s==Bseas(g)) {natage(t,p,g,0) = Recruits*recr_dist_forecast(GP(g),s,p)*submorphdist(GP2(g));} //get recruits from total recruits and RecrDist
              totbio+=natage(t,p,g)*Wt_Age_beg(s,g);
              temp=natage(t,p,g)(Smry_Age,nages)*Wt_Age_beg(s,g)(Smry_Age,nages);
              smrybio+=temp;
              smrynum+=sum(natage(t,p,g)(Smry_Age,nages));
              if(save_for_report==1)
              {
                Save_PopLen(t,p,g)=0.0;
                Save_PopLen(t,p+pop,g)=0.0;  // later put midseason here
                Save_PopWt(t,p,g)=0.0;
                Save_PopWt(t,p+pop,g)=0.0;  // later put midseason here
                Save_PopAge(t,p,g)=value(natage(t,p,g));
                for (a=0;a<=nages;a++)
                {
                  Save_PopLen(t,p,g)+=value(natage(t,p,g,a))*value(ALK(s,g,a));
                  Save_PopWt(t,p,g)+= value(natage(t,p,g,a))*value(elem_prod(ALK(s,g,a),wt_len(s,sx(g))));
                } // close age loop
              }
            }
//            report5<<y<<" start_num "<<natage(t,p,1)(0,8)<<endl;
            Tune_F_loops=1;
            for(f=1;f<=Nfleet;f++)
            {
              switch (ABC_Loop)
              {
                case 1:
                {
                  Hrate(f,t)=Fcast_Fmult*Fcast_RelF_Use(s,f);
                  break;  // no action, keep Hrate
                }
                case 2:
                {
                  Hrate(f,t)=ABC_buffer(y)*Fcast_Fmult*Fcast_RelF_Use(s,f);
                  if(Fcast_InputCatch(t,f,1)>-1.0)  //  have an input
                  {
                    if(Fcast_InputCatch(t,f,2)<=3)  //  input is catch
                      {
                        if(Fcast_InputCatch(t,f,1)==0.0)
                        {
                          Hrate(f,t)=0.0;
                          Do_F_tune(t,f)=0;
                        }
                        else
                        {
                          Tune_F_loops=8;
                          Do_F_tune(t,f)=1;
                        }
                      }
                    else
                      {Hrate(f,t)=Fcast_InputCatch(t,f,1);}  // input is as Hrate (F), but do not need tuning
                  }
                  break;
                }
                case 3:  //  always get F to match catch when in ABC_Loop==3
                {
                  Tune_F_loops=8;
                  Do_F_tune(t,f)=1;
                  break;
                }
              }
            }
            if(F_Method==1)  //  Pope's
            {
              for (g=1;g<=gmorph;g++)
              if(use_morph(g)>0)
              {
                Nmid(g) = elem_prod(natage(t,p,g),surv1(s,GP3(g)));
              }

              for(Tune_F=1;Tune_F<=Tune_F_loops;Tune_F++)
              {
                for (f=1;f<=Nfleet;f++)   // get calculated catch
                if (fleet_area(f)==p && Fcast_RelF_Use(s,f)>0.0)
                {
                  temp=0.0;
                  if(Do_F_tune(t,f)==1)  // have an input catch, so get expected catch from F and Z
                  {
                    if(ABC_Loop==2)  //  tune to input catch
                    {
                      for (g=1;g<=gmorph;g++)
                      if(use_morph(g)>0)
                      {
                        if(catchunits(f)==1)  //  catch in weight
                        {
                          if(Fcast_InputCatch(t,f,2)==2)
                            {temp+=Nmid(g)*deadfish_B(s,g,f);}      // dead catch bio
                          else if(Fcast_InputCatch(t,f,2)==3)
                            {temp+=Nmid(g)*sel_al_2(s,g,f);}      // retained catch bio
                        }
                        else   //  catch in numbers
                        {
                          if(Fcast_InputCatch(t,f,2)==2)
                          {temp+=Nmid(g)*deadfish(s,g,f);}      // deadfish catch numbers
                          else if(Fcast_InputCatch(t,f,2)==3)
                          {temp+=Nmid(g)*sel_al_4(s,g,f);}      // retained catch numbers
                        }
                      }  //close gmorph loop
                      temp=max_harvest_rate-Fcast_InputCatch(t,f,1)/(temp+NilNumbers);
                      Hrate(f,t)=max_harvest_rate-posfun(temp,0.0001,Fcast_Crash);
                    }
                    else  //  tune to adjusted catch calculated from ABC_Loop=2 (note different basis for catch)
                    {
                      for (g=1;g<=gmorph;g++)
                      if(use_morph(g)>0)
                      {
                        if(Fcast_Catch_Basis==2)
                        {temp+=Nmid(g)*deadfish_B(s,g,f);}      // dead catch bio
                        else if(Fcast_Catch_Basis==3)
                        {temp+=Nmid(g)*sel_al_2(s,g,f);}      // retained catch bio
                        else if(Fcast_Catch_Basis==5)
                        {temp+=Nmid(g)*deadfish(s,g,f);}      // deadfish catch numbers
                        else if(Fcast_Catch_Basis==6)
                        {temp+=Nmid(g)*sel_al_4(s,g,f);}      // retained catch numbers
                      }  //close gmorph loop
                      temp=max_harvest_rate-Fcast_Catch_Store(t,f)/(temp+NilNumbers);
                      Hrate(f,t)=max_harvest_rate-posfun(temp,0.0001,Fcast_Crash);
                    }
                  }  // end have fixed catch to be matched
                }  // end fishery loop
              }  //  end finding the Hrates
              
//  now get catch details and survivorship
              Nsurv=Nmid;  //  initialize the number of survivors
              for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       ProgLabel_105
              if (fleet_area(f)==p)
              {
                catch_fleet(t,f).initialize();
                temp=Hrate(f,t);
                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
//                  Nmid(g) = elem_prod(natage(t,p,g),surv1(s,GP3(g)));
                  catch_fleet(t,f,1)+=Nmid(g)*sel_al_1(s,g,f);      // encountered catch bio
                  catch_fleet(t,f,2)+=Nmid(g)*deadfish_B(s,g,f);      // dead catch bio
                  catch_fleet(t,f,3)+=Nmid(g)*sel_al_2(s,g,f);      // retained catch bio
                  catch_fleet(t,f,4)+=Nmid(g)*sel_al_3(s,g,f);      // encountered catch numbers
                  catch_fleet(t,f,5)+=Nmid(g)*deadfish(s,g,f);      // deadfish catch numbers
                  catch_fleet(t,f,6)+=Nmid(g)*sel_al_4(s,g,f);      // retained catch numbers
                  catage_w(g)= temp*elem_prod(Nmid(g),deadfish(s,g,f));
                  Nsurv(g)-=catage_w(g);
                }  //close gmorph loop
                catch_fleet(t,f)*=temp;
              }  // close fishery
              for (g=1;g<=gmorph;g++)
              if(use_morph(g)>0)
              {
                if(s<nseas) natage(t+1,p,g,0) = Nsurv(g,0)*surv1(s,GP3(g),0);  // advance age zero within year
                for (a=1;a<nages;a++) {natage(t+1,p,g,a) = Nsurv(g,a-adv_age)*surv1(s,GP3(g),a-adv_age);}
                natage(t+1,p,g,nages) = Nsurv(g,nages)*surv1(s,GP3(g),nages);   // plus group
                if(s==nseas) natage(t+1,p,g,nages) += Nsurv(g,nages-1)*surv1(s,GP3(g),nages-1);
                if(save_for_report==1)
                {
                  j=p+pop;
                  for (a=0;a<=nages;a++)
                  {
                    Save_PopLen(t,j,g)+=value(0.5*(Nmid(g,a)+Nsurv(g,a)))*value(ALK(smid,g,a));
                    Save_PopWt(t,j,g)+= value(0.5*(Nmid(g,a)+Nsurv(g,a)))*value(elem_prod(ALK(smid,g,a),wt_len(s,sx(g))));
                    Save_PopAge(t,j,g,a)=value(0.5*(Nmid(g,a)+Nsurv(g,a)));
                  } // close age loop
                }
              }
//              report5<<"NatAge"<<natage(t,1,1)(0,10)<<endl;
//              report5<<"surv1"<<surv1(1,GP3(1))(0,10)<<endl;
//              report5<<"Nmid"<<Nmid(1)(0,10)<<endl;
//              report5<<"deadfish"<<deadfish(1,1,10)(0,10)<<endl;
//              report5<<"deadfishB"<<deadfish_B(1,1,10)(0,10)<<endl;
//              report5<<"catchfleet ";
//              for(f=1;f<=Nfleet;f++) report5<<catch_fleet(t,f,2)<<" ";
//              report5<<endl;
//              report5<<"Nsurv"<<Nsurv(1)(0,10)<<endl;
//              report5<<"next_N"<<natage(t+1,1,1)(0,10)<<endl;
           }  //  end Fmethod=1 pope

            else  //  continuous F
            {
              for(Tune_F=1;Tune_F<=Tune_F_loops;Tune_F++)  //  tune F to match catch
              {
                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
                  Z_rate(t,p,g)=natM(s,GP3(g));
                  for (f=1;f<=Nfleet;f++)       //loop over fishing fleets to get Z       ProgLabel_105
                  if (fleet_area(f)==p && Fcast_RelF_Use(s,f)>0.0)
                  {
                    Z_rate(t,p,g)+=deadfish(s,g,f)*Hrate(f,t);
                  }
                  Zrate2(p,g)=elem_div( (1.-mfexp(-seasdur(s)*Z_rate(t,p,g))), Z_rate(t,p,g));
                }  //  end morph

                for (f=1;f<=Nfleet;f++)   // get calculated catch
                if (fleet_area(f)==p && Fcast_RelF_Use(s,f)>0.0)
                {
                  temp=0.0;
                  if(Do_F_tune(t,f)==1)  // have an input catch, so get expected catch from F and Z
                  {
                    if(ABC_Loop==2)  //  tune to input catch
                    {
                      for (g=1;g<=gmorph;g++)
                      if(use_morph(g)>0)
                      {
                        if(catchunits(f)==1)  //  catch in weight
                        {
                          if(Fcast_InputCatch(t,f,2)==2)
                            {temp+=elem_prod(natage(t,p,g),deadfish_B(s,g,f))*Zrate2(p,g);}      // dead catch bio
                          else if(Fcast_InputCatch(t,f,2)==3)
                            {temp+=elem_prod(natage(t,p,g),sel_al_2(s,g,f))*Zrate2(p,g);}      // retained catch bio
                        }
                        else   //  catch in numbers
                        {
                          if(Fcast_InputCatch(t,f,2)==2)
                          {temp+=elem_prod(natage(t,p,g),deadfish(s,g,f))*Zrate2(p,g);}      // deadfish catch numbers
                          else if(Fcast_InputCatch(t,f,2)==3)
                          {temp+=elem_prod(natage(t,p,g),sel_al_4(s,g,f))*Zrate2(p,g);}      // retained catch numbers
                        }
                      }  //close gmorph loop
                      temp*=Hrate(f,t);
                      H_temp(f)=Hrate(f,t);
                      C_temp(f)=temp;
                      if(Tune_F<3)
                      {
                        C_old(f)=C_temp(f);
                        H_old(f)=H_temp(f);
                        Hrate(f,t)*=(Fcast_InputCatch(t,f,1)+1.0)/(temp+1.0);  //  apply adjustment
                      }
                      else
                      {
                        Hrate(f,t)=H_old(f)+(H_temp(f)-H_old(f))/(C_temp(f)-C_old(f)+1.0e-6) * (Fcast_InputCatch(t,f,1)-C_old(f));
                        C_old(f)=C_temp(f);
                        H_old(f)=H_temp(f);
                      }
                    }
                    else  //  tune to adjusted catch calculated in ABC_Loop=2 (note different basis for catch)
                    {
                      for (g=1;g<=gmorph;g++)
                      if(use_morph(g)>0)
                      {
                        if(Fcast_Catch_Basis==2)
                        {temp+=elem_prod(natage(t,p,g),deadfish_B(s,g,f))*Zrate2(p,g);}      // dead catch bio
                        else if(Fcast_Catch_Basis==3)
                        {temp+=elem_prod(natage(t,p,g),sel_al_2(s,g,f))*Zrate2(p,g);}      // retained catch bio
                        else if(Fcast_Catch_Basis==5)
                        {temp+=elem_prod(natage(t,p,g),deadfish(s,g,f))*Zrate2(p,g);}      // deadfish catch numbers
                        else if(Fcast_Catch_Basis==6)
                        {temp+=elem_prod(natage(t,p,g),sel_al_4(s,g,f))*Zrate2(p,g);}      // retained catch numbers
                      }  //close gmorph loop
                      temp*=Hrate(f,t);
//                      Hrate(f,t)*=(Fcast_Catch_Store(t,f)+1.0)/(temp+1.0);  //  apply adjustment
                      H_temp(f)=Hrate(f,t);
                      C_temp(f)=temp;
                      if(Tune_F<3)
                      {
                        C_old(f)=C_temp(f);
                        H_old(f)=H_temp(f);
                        Hrate(f,t)*=(Fcast_Catch_Store(t,f)+1.0)/(temp+1.0);  //  apply adjustment
//                        report5<<" new "<<Hrate(f,t)<<endl;
                      }
                      else
                      {
                        if(Tune_F<7)
                        {Hrate(f,t)=(H_old(f)+(H_temp(f)-H_old(f))/(C_temp(f)-C_old(f)+1.0e-6) * (Fcast_Catch_Store(t,f)-C_old(f)));}
                        else if(Tune_F==7)
                        {Hrate(f,t)=(H_old(f)+(H_temp(f)-H_old(f))/(C_temp(f)-C_old(f)+1.0e-6) * (Fcast_Catch_Store(t,f)-C_old(f)));}
//                        report5<<" new "<<Hrate(f,t)<<endl;
                        C_old(f)=C_temp(f);
                        H_old(f)=H_temp(f);
                      }

                    }
//                    report5<<t<<" "<<Tune_F<<" "<<f<<" "<<Fcast_InputCatch(t,f)<<" "<<Fcast_Catch_Store(t,f)<<" "<<temp<<" new Hrate "<<Hrate(f,t)<<endl;
                  }  // end have fixed catch to be matched
                }  // end fishery loop
//                    if(y==endyr+2) report5<<"Tune "<<Fcast_Loop1<<" "<<ABC_Loop<<" "<<Tune_F<<" "<<Fcast_Catch_Store(t,2)<<" "<<temp<<" "<<Hrate(2,t)<<endl;
              }  //  done tuning F

              for (f=1;f<=Nfleet;f++)       //loop over fishing fleets       ProgLabel_105
              if (fleet_area(f)==p)
              {
                catch_fleet(t,f).initialize();
                for (g=1;g<=gmorph;g++)
                if(use_morph(g)>0)
                {
                  tempvec_a=Hrate(f,t)*Zrate2(p,g);
                  catch_fleet(t,f,1)+=tempvec_a*elem_prod(natage(t,p,g),sel_al_1(s,g,f));      // encountered catch bio
                  catch_fleet(t,f,2)+=tempvec_a*elem_prod(natage(t,p,g),deadfish_B(s,g,f));      // dead catch bio
                  catch_fleet(t,f,3)+=tempvec_a*elem_prod(natage(t,p,g),sel_al_2(s,g,f));      // retained catch bio
                  catch_fleet(t,f,4)+=tempvec_a*elem_prod(natage(t,p,g),sel_al_3(s,g,f));      // encountered catch numbers
                  catch_fleet(t,f,5)+=tempvec_a*elem_prod(natage(t,p,g),deadfish(s,g,f));      // deadfish catch numbers
                  catch_fleet(t,f,6)+=tempvec_a*elem_prod(natage(t,p,g),sel_al_4(s,g,f));      // retained catch numbers
                }  //close gmorph loop
              }  // close fishery
              for (g=1;g<=gmorph;g++)
              if(use_morph(g)>0)
              {
                if(s<nseas) natage(t+1,p,g,0) = natage(t,p,g,0)*mfexp(-Z_rate(t,p,g,0)*seasdur(s));  // advance age zero within year
                for (a=1;a<nages;a++) {natage(t+1,p,g,a) = natage(t,p,g,a-adv_age)*mfexp(-Z_rate(t,p,g,a-adv_age)*seasdur(s));}
                natage(t+1,p,g,nages) = natage(t,p,g,nages)*mfexp(-Z_rate(t,p,g,nages)*seasdur(s));   // plus group
//                if(g==1) report5<<y<<" survivors "<<natage(t+1,p,g)(0,8)<<endl;
                if(s==nseas) natage(t+1,p,g,nages) += natage(t,p,g,nages-1)*mfexp(-Z_rate(t,p,g,nages-1)*seasdur(s));
                if(save_for_report==1)
                {
                  j=p+pop;
                  for (a=0;a<=nages;a++)
                  {
                    Save_PopLen(t,j,g)+=value(natage(t,p,g,a)*mfexp(-Z_rate(t,p,g,a)*0.5*seasdur(s)))*value(ALK(smid,g,a));
                    Save_PopWt(t,j,g)+=value(natage(t,p,g,a)*mfexp(-Z_rate(t,p,g,a)*0.5*seasdur(s)))*value(elem_prod(ALK(smid,g,a),wt_len(s,sx(g))));
                    Save_PopAge(t,j,g,a)=value(natage(t,p,g,a)*mfexp(-Z_rate(t,p,g,a)*0.5*seasdur(s)));
                  } // close age loop
                }
              }  // end morph loop
            }  // end continuous F

//  put check for have_data here

            if(Hermaphro_Option>0)  //hermaphroditism
            {
              if(Hermaphro_seas==-1 || Hermaphro_seas==s)
              {
                k=gmorph/2;
                for (g=1;g<=k;g++)  //  loop females
                if(use_morph(g)>0)
                {
                  for (a=1;a<nages;a++)
                  {
                    natage(t+1,p,g+k,a) += natage(t+1,p,g,a)*Hermaphro_val(GP4(g),a-1); // increment males
                    natage(t+1,p,g,a) *= (1.-Hermaphro_val(GP4(g),a-1)); // decrement females
                  }
                }
              }
            }

           if(show_MSY==1)
           {
            report5<<p<<" "<<y<<" "<<ABC_Loop<<" "<<s<<" "<<ABC_buffer(y)<<" "<<totbio<<" "<<smrybio<<" ";
            if(s==spawn_seas)
            {report5<<SPB_current<<" "<<SPB_current/SPB_virgin<<" "<<Recruits;}
            else
            {report5<<0<<" "<<0<<" "<<0;}
            for (f=1;f<=Nfleet;f++)
            {
            if(fleet_area(f)==p)
            {
              if(F_Method==1)
              {report5<<" "<<catch_fleet(t,f)(1,6)<<" "<<Hrate(f,t);}
              else
              {report5<<" "<<catch_fleet(t,f)(1,6)<<" "<<Hrate(f,t)*seasdur(s);}
            }
            else
            {report5<<" - - - - - - - ";}
            if(Fcast_InputCatch(t,f,1)<0.0) {report5<<" R ";} else {report5<<" C ";}
            }
            if(s==nseas&&Max_Fcast_Catch(p)>0.) {report5<<" "<<Max_Fcast_Catch(p);} else {report5<<" NA ";} //  a max catch has been set for this area
           }
           if(p<pop && show_MSY==1) report5<<endl;

           if(s==1&&Fcast_Loop1==Fcast_Loop_Control(1))
           {
             Smry_Table(y,1)+=totbio;
             Smry_Table(y,2)+=smrybio;  // in forecast
             Smry_Table(y,3)+=smrynum;   //sums to accumulate across submorphs and birthseasons
           }
          }  //  end loop of areas
          if(do_migration>0)  // movement between areas in forecast
          {
            natage_temp=natage(t+1);
            natage(t+1).initialize();
            for (p=1;p<=pop;p++)  //   source population
            for (p2=1;p2<=pop;p2++)  //  destination population
            for (g=1;g<=gmorph;g++)
            if(use_morph(g)>0)
              {
            	 k=move_pattern(s,GP4(g),p,p2);
              if(k>0) natage(t+1,p2,g) += elem_prod(natage_temp(p,g),migrrate(bio_yr,k));
              }
          }

          if( (save_for_report>0) || ((sd_phase() || mceval_phase()) && (initial_params::mc_phase==0)) )
          {
            if(Fcast_Loop1==2 && ABC_Loop==1)  // get variance in OFL
            {
              for(f=1;f<=Nfleet;f++) {Mgmt_quant(Fcast_catch_start+N_Fcast_Yrs+y-endyr)+=catch_fleet(t,f,2);}
            }

            if(Fcast_Loop1==Fcast_Loop_Control(1) && ABC_Loop==3)  //  in final loop, so do variance quantities
            {
              if(STD_Yr_Reverse_F(y)>0)
              {
                if(F_reporting<=1)
                {
                  for(f=1;f<=Nfleet;f++) F_std(STD_Yr_Reverse_F(y))+=catch_fleet(t,f,2);   // add up dead catch biomass
                  if(s==nseas) F_std(STD_Yr_Reverse_F(y))/=Smry_Table(y,2);
                }
                else if(F_reporting==2)
                {
                  for(f=1;f<=Nfleet;f++) F_std(STD_Yr_Reverse_F(y))+=catch_fleet(t,f,5);   // add up dead catch numbers
                  if(s==nseas) F_std(STD_Yr_Reverse_F(y))/=Smry_Table(y,3);
                }
                else if(F_reporting==3)
                {
                  if(F_Method==1)
                  {
                    for(f=1;f<=Nfleet;f++)
                    {
                      F_std(STD_Yr_Reverse_F(y))+=Hrate(f,t);
                    }
                  }
                  else
                  {
                    for(f=1;f<=Nfleet;f++)
                    {
                      F_std(STD_Yr_Reverse_F(y))+=Hrate(f,t)*seasdur(s);
                    }
                  }
                }
                else if(F_reporting==4 && s==nseas)
                {
        //  sum across p and g the number of survivors to end of the year
        //  also project from the initial numbers and M, the number of survivors without F
        //  then F = ln(n+1/n)(M+F) - ln(n+1/n)(M only), but ln(n) cancels out, so only need the ln of the ratio of the two ending quantities
                  temp1=0.0;
                  temp2=0.0;
                  for(g=1;g<=gmorph;g++)
                  if(use_morph(g)>0)
                  {
                    for(p=1;p<=pop;p++)
                    {
                      for(a=F_reporting_ages(1);a<=F_reporting_ages(2);a++)   //  should not let a go higher than nages-2 because of accumulator
                      {
                        if(nseas==1)
                        {
                          temp1+=natage(t+1,p,g,a+1);
                          temp2+=natage(t,p,g,a)*mfexp(-seasdur(s)*natM(s,GP3(g),a));
                        }
                        else
                        {
                          temp1+=natage(t+1,p,g,a+1);
                          temp3=natage(t-nseas+1,p,g,a);  //  numbers at begin of year
                          for(j=1;j<=nseas;j++) {temp3*=mfexp(-seasdur(j)*natM(j,GP3(g),a));}
                          temp2+=temp3;
                        }
                      }
                    }
                  }
                  F_std(STD_Yr_Reverse_F(y)) = log(temp2)-log(temp1);
                }
              }

              for(f=1;f<=Nfleet;f++)
              {
                Mgmt_quant(Fcast_catch_start+y-endyr)+=catch_fleet(t,f,2);
                if(Do_Retain==1) Mgmt_quant(Fcast_catch_start+2*N_Fcast_Yrs+y-endyr)+=catch_fleet(t,f,3);
              }
            }
          }

          //  store catches to allow calc of adjusted F to match this catch when doing ABC_loop=3, and then when doing Fcast_loop1=3
          for(f=1;f<=Nfleet;f++)
          {
            Fcast_Catch_Store(t,f)=catch_fleet(t,f,Fcast_Catch_Basis);
            totcatch+=Fcast_Catch_Store(t,f);
          }
//        report5<<Fcast_Loop1<<" "<<y<<" "<<s<<" "<<ABC_Loop<<" "<<SPB_current<<" "<<Recruits<<" "<<ABC_buffer(y)<<" "<<Fcast_Fmult<<" "<<Hrate(1,t);
//        report5<<" "<<Hrate(2,t)<<" "<<Hrate(3,t)<<" catch "<<Fcast_Catch_Store(t)<<" ";
          if(show_MSY==1)
          {
            if(s==nseas) {report5<<" "<<totcatch<<" ";} else {report5<<" NA ";}
            if(s==nseas && STD_Yr_Reverse_F(y)>0) {report5<<F_std(STD_Yr_Reverse_F(y));} else {report5<<" NA ";}
            report5<<endl;
          }
//          report5<<y<<" end_seas "<<natage(t+1,1,1)(0,8)<<endl;
        }  //  end loop of seasons

        if(ABC_Loop==2)
        {
          // calculate annual catch for each fleet
          Fcast_Catch_Current_Alloc.initialize();
          for(f=1;f<=Nfleet;f++)
          for(s=1;s<=nseas;s++)
          {
            t=t_base+s;
            Fcast_Catch_Current_Alloc(f)+=catch_fleet(t,f,Fcast_Catch_Basis); //  accumulate annual catch according to catch basis (2=deadbio, 3=ret bio, 5=dead num, 6=ret num)
          }

          if(Fcast_Do_Fleet_Cap>0 && y>=Fcast_Cap_FirstYear)
          {
            for(f=1;f<=Nfleet;f++)   //  adjust ABC catch to fleet caps
            {
              if(Fcast_MaxFleetCatch(f)>0.)
              {
                temp = Fcast_Catch_Current_Alloc(f)/Fcast_MaxFleetCatch(f);
                join1=1./(1.+mfexp(1000.*(temp-1.0)));  // steep logistic joiner at adjustment of 1.0
                temp1=join1*1.0 + (1.-join1)*temp;
                Fcast_Catch_Current_Alloc(f)/=temp1;
                for(s=1;s<=nseas;s++)
                {Fcast_Catch_Store(t_base+s,f)/=temp1;}
              }
            }
//            report5<<Tune_F<<" tune_fleets "<<Fcast_Catch_Current_Alloc<<endl;
          }

          if(Fcast_Do_Area_Cap>0  && y>=Fcast_Cap_FirstYear)  // scale down if Totcatch exceeds Max_Fcast_Catch (in this area)
          {
            if(pop==1)  // one area
            {
              Fcast_Catch_ByArea(1)=sum(Fcast_Catch_Current_Alloc(1,Nfleet));
            }
            else
            {
              Fcast_Catch_ByArea=0.0;
              for(f=1;f<=Nfleet;f++)
              {
                Fcast_Catch_ByArea(fleet_area(f))+=Fcast_Catch_Current_Alloc(f);
              }
            }
            for(p=1;p<=pop;p++)
            if(Max_Fcast_Catch(p)>0.0)
            {
              temp = Fcast_Catch_ByArea(p)/Max_Fcast_Catch(p);
              join1=1./(1.+mfexp(1000.*(temp-1.0)));  // steep logistic joiner at adjustment of 1.0
              temp1=join1*1.0 + (1.-join1)*temp;
              for(f=1;f<=Nfleet;f++)
              if (fleet_area(f)==p)
              {
                Fcast_Catch_Current_Alloc(f)/=temp1;  // adjusts total for the year
                for(s=1;s<=nseas;s++)
                {Fcast_Catch_Store(t_base+s,f)/=temp1;}
              }
            }
//            report5<<Tune_F<<" tune_area"<<Fcast_Catch_Current_Alloc<<endl;
          }

          if(Fcast_Catch_Allocation_Groups>0  && y>=Fcast_Cap_FirstYear)  // adjust to get a specific fleet allocation
          {
            Fcast_Catch_Current_Allocation=0.0;
            for(g=1;g<=Fcast_Catch_Allocation_Groups;g++)
            for(f=1;f<=Nfleet;f++)
            if (Allocation_Fleet_Assignments(f)==g)
            {
               Fcast_Catch_Current_Allocation(g)+=Fcast_Catch_Current_Alloc(f);
            }
            temp=sum(Fcast_Catch_Current_Allocation);  // total catch for all fleets that are part of the allocation scheme
            temp1=sum(Fcast_Catch_Allocation);  // total of all allocation fractions for all fleets that are part of the allocation scheme

            for(g=1;g<=Fcast_Catch_Allocation_Groups;g++)
            {
              temp2=(Fcast_Catch_Allocation(g)/temp1) / (Fcast_Catch_Current_Allocation(g)/temp);
              for(f=1;f<=Nfleet;f++)
              if (Allocation_Fleet_Assignments(f)==g)
              {
                Fcast_Catch_Current_Alloc(f)*=temp2;
                for(s=1;s<=nseas;s++)
                {
                  Fcast_Catch_Store(t_base+s,f)*=temp2;
                }
              }
            }
//            report5<<Tune_F<<" tune_alloc"<<Fcast_Catch_Current_Alloc<<endl;
          }  //  end allocation among groups
        }
//          report5<<y<<" end_ABCloop "<<natage(t+1,1,1)(0,8)<<endl;
      }  //  end ABC_Loop

      if( (Fcast_Loop1==Fcast_Loop_Control(1) && (save_for_report>0)) || ((sd_phase() || mceval_phase()) && (initial_params::mc_phase==0)) )
      {
        eq_yr=y; equ_Recr=Recr_virgin; bio_yr=endyr;
        Fishon=0;
        Do_Equil_Calc();                      //  call function to do equilibrium calculation
        SPR_unf=SPB_equil;
        Smry_Table(y,11)=SPR_unf;
        Smry_Table(y,13)=GenTime;
        Fishon=1;
        Do_Equil_Calc();                      //  call function to do equilibrium calculation
        SPR_trial=SPB_equil;
        if(STD_Yr_Reverse_Ofish(y)>0) SPR_std(STD_Yr_Reverse_Ofish(y))=SPR_trial/SPR_unf;
        Smry_Table(y,9)=totbio;
        Smry_Table(y,10)=smrybio;
        Smry_Table(y,12)=SPR_trial;
        Smry_Table(y,14)=YPR_dead;
        for(g=1;g<=gmorph;g++)
        {
          Smry_Table(y,20+g)=(cumF(g));
          Smry_Table(y,20+gmorph+g)=(maxF(g));
        }
      }
    }  //  end year loop
  }  //  end Fcast_Loop1  for the different stages of the forecast

  }
//  end forecast function

//*********************************************************************
 /*  ProgLabel_FUNCTION 36 write_summaryoutput */
FUNCTION void write_summaryoutput()
  {
  random_number_generator radm(long(time(&finish)));

  time(&finish);
  elapsed_time = difftime(finish,start);
  report2<<runnumber<<" -logL: "<<obj_fun<<" Spbio(Vir_Start_End): "<<SPB_yr(styr-2)<<" "<<SPB_yr(styr)<<" "<<SPB_yr(endyr)<<endl;
  report2<<runnumber<<" Files: "<<datfilename<<" "<<ctlfilename;
  if(readparfile>=1) report2<<" Start_from_SS3.PAR";
  report2<<endl<<runnumber<<" N_iter: "<<niter<<" runtime(sec): "<<elapsed_time<<" starttime: "<<ctime(&start);
  report2<<runnumber<<" "<<version_info<<endl;
  report2<<runnumber<<" F_Method: "<<F_Method<<" Retro_YR: "<<retro_yr<<" Forecast_Type: "<<Do_Forecast<<" MSY_Type: "<<Do_MSY<<endl;
  if(N_SC>0)
  {
    for(j=1;j<=N_SC;j++) report2<<runnumber<<" Comment S_"<<j<<" "<<Starter_Comments(j)<<endl;
  }
  if(N_DC>0)
  {
    for(j=1;j<=N_DC;j++) report2<<runnumber<<" Comment D_"<<j<<" "<<Data_Comments(j)<<endl;
  }
  if(N_CC>0)
  {
    for(j=1;j<=N_CC;j++) report2<<runnumber<<" Comment C_"<<j<<" "<<Control_Comments(j)<<endl;
  }
  if(N_FC>0)
  {
    for(j=1;j<=N_FC;j++) report2<<runnumber<<" Comment F_"<<j<<" "<<Forecast_Comments(j)<<endl;
  }
  k=current_phase();
  if(k>max_lambda_phase) k=max_lambda_phase;
  report2<<runnumber<<" Like_Emph Total 1 "<<endl<<runnumber<<" Like_Value Total "<<obj_fun<<endl;
  if(nobs_cr>0) report2<<runnumber<<" Like_Emph Indices All "<<column(surv_lambda,k)<<endl<<
  runnumber<<" Like_Value Indices "<<surv_like*column(surv_lambda,k)<<" " <<surv_like<<endl;
  if(nobs_disc>0) report2<<runnumber<<" Like_Emph Discard All "<<column(disc_lambda,k)<<endl<<
  runnumber<<" Like_Value Discard "<<disc_like*column(disc_lambda,k)<<" " <<disc_like<<endl;
  if(nobs_mnwt>0) report2<<runnumber<<" Like_Emph MeanBodyWt All "<<column(mnwt_lambda,k)<<endl<<
  runnumber<<" Like_Value MeanBodyWt "<<mnwt_like*column(mnwt_lambda,k)<<" " <<mnwt_like<<endl;
  if(Nobs_l_tot>0) report2<<runnumber<<" Like_Emph LenComp All "<<column(length_lambda,k)<<endl<<
  runnumber<<" Like_Value LenComp "<<length_like*column(length_lambda,k)<<" " <<length_like<<endl;
  if(Nobs_a_tot>0) report2<<runnumber<<" Like_Emph AgeComp All "<<column(age_lambda,k)<<endl<<
  runnumber<<" Like_Value AgeComp "<<age_like*column(age_lambda,k)<<" " <<age_like<<endl;
  if(nobs_ms_tot>0) report2<<runnumber<<" Like_Emph MeanLAA All "<<column(sizeage_lambda,k)<<endl<<
  runnumber<<" Like_Value MeanLAA "<<sizeage_like*column(sizeage_lambda,k)<<" " <<sizeage_like<<endl;
  if(F_Method>1) report2<<runnumber<<" Like_Emph Catch All "<<column(catch_lambda,k)<<endl<<
  runnumber<<" Like_Value Catch "<<catch_like*column(catch_lambda,k)<<" " <<catch_like<<endl;
  if(SzFreq_Nmeth>0) report2<<runnumber<<" Like_Emph WeightFreq All "<<column(SzFreq_lambda,k)<<endl<<
  runnumber<<" Like_Value WeightFreq "<<SzFreq_like*column(SzFreq_lambda,k)<<" " <<SzFreq_like<<endl;
  if(Do_Morphcomp>0) report2<<runnumber<<" Like_Emph Morphcomp All "<<Morphcomp_lambda(k)<<endl<<
  runnumber<<" Like_Value Morphcomp "<<Morphcomp_like*Morphcomp_lambda(k)<<" " <<Morphcomp_like<<endl;
  if(Do_TG>0) report2<<runnumber<<" Like_Emph Tag_comp All "<<column(TG_lambda1,k)<<endl<<
  runnumber<<" Like_Value Tag_comp "<<TG_like1*column(TG_lambda1,k)<<" " <<TG_like1<<endl;
  if(Do_TG>0) report2<<runnumber<<" Like_Emph Tag_negbin All "<<column(TG_lambda2,k)<<endl<<
  runnumber<<" Like_Value Tag_negbin "<<TG_like2*column(TG_lambda2,k)<<" " <<TG_like2<<endl;

  report2<<runnumber<<" Like_Comp Equ_Catch Recruits Fcast_Recr Biasadj Priors ParmDevs CrashPen"<<endl;
  report2<<runnumber<<" Like_Emph "<<init_equ_lambda(k)<<" "<<recrdev_lambda(k)<<" " <<Fcast_recr_lambda<<" "
         <<parm_prior_lambda(k)<<" " <<parm_dev_lambda(k)<<" " <<CrashPen_lambda(k)<<endl;
  report2<<runnumber<<" Like_Value*Emph "<<equ_catch_like*init_equ_lambda(k)<<" "<<recr_like*recrdev_lambda(k)<<" "
         <<Fcast_recr_like<<" "<<parm_like*parm_prior_lambda(k)<<" "<<parm_dev_like*parm_dev_lambda(k)<<" "<<CrashPen*CrashPen_lambda(k)<<endl;

  report2 <<runnumber<<" TimeSeries Year Vir Equ "<<years<<" ";
  k=YrMax;
  if(k==endyr) k=endyr+1;
  for(y=endyr+1;y<=k;y++) {report2<<y<<"F ";}
  report2 <<endl;
  report2 <<runnumber<<" Timeseries Spbio "<<column(Smry_Table,7)<<endl;
  report2 <<runnumber<<" Timeseries Recruit "<<column(Smry_Table,8)<<endl;
  report2 <<runnumber<<" Timeseries TotBio "<<column(Smry_Table,1)<<endl;
  report2 <<runnumber<<" Timeseries SmryBio-"<<Smry_Age<<" "<<column(Smry_Table,2)<<endl;
  report2 <<runnumber<<" Timeseries TotCatch "<<column(Smry_Table,4)<<endl;
  report2 <<runnumber<<" Timeseries RetCatch "<<column(Smry_Table,5)<<endl;
  if(Do_Benchmark>0) report2<<runnumber<<" Mgmt_Quant "<<Mgmt_quant(1,6+Do_Retain)<<endl;

  report2<<runnumber<<" Parm Labels ";
  for(i=1;i<=ParCount;i++) {report2<<" "<<ParmLabel(i);}
  report2<<endl;
  report2<<runnumber<<" Parm Values ";
  report2<<" "<<MGparm<<" ";
  if(N_MGparm_dev>0) report2<<MGparm_dev<<" ";
  report2<<SR_parm<<" ";
  if(recdev_cycle>0) report2<<recdev_cycle_parm<<" ";
  if(recdev_do_early>0) report2<<recdev_early<<" ";
  if(do_recdev==1) {report2<<recdev1<<" ";}
  if(do_recdev==2) {report2<<recdev2<<" ";}
  if(Do_Forecast>0) report2<<Fcast_recruitments<<" "<<Fcast_impl_error<<" ";
  report2<<init_F<<" ";
  if(F_Method==2) report2<<" "<<F_rate;
  if(Q_Npar>0) report2<<Q_parm<<" ";
  report2<<selparm<<" ";
  if(N_selparm_dev>0) report2<<selparm_dev<<" ";
  if(Do_TG>0) report2<<TG_parm<<" ";
  report2<<endl;

  NP=0;   // count of number of parameters
  report2<<runnumber<<" MG_parm ";
  for(j=1;j<=N_MGparm2;j++)
  {NP++; report2<<" "<<ParmLabel(NP);}
  report2<<endl<<runnumber<<" MG_parm "<<MGparm<<endl;

  if(N_MGparm_dev>0)
  {
    report2<<runnumber<<" MG_parm_dev ";
    for(i=1;i<=N_MGparm_dev;i++)
    for(j=MGparm_dev_minyr(i);j<=MGparm_dev_maxyr(i);j++)
    {NP++; report2<<" "<<ParmLabel(NP);}
    report2<<endl<<runnumber<<" MG_parm_dev "<<MGparm_dev<<endl;
  }

    report2<<runnumber<<" SR_parm ";
    for(i=1;i<=N_SRparm2+recdev_cycle;i++)
    {NP++; report2<<" "<<ParmLabel(NP);}
    report2<<endl<<runnumber<<" SR_parm "<<SR_parm<<" ";
    if(recdev_cycle>0) report2<<recdev_cycle_parm;
    report2<<endl;

    if(recdev_do_early>0)
    {
      report2<<runnumber<<" Recr_early ";
      for(i=recdev_early_start;i<=recdev_early_end;i++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Recr_early ";
      for(i=recdev_early_start;i<=recdev_early_end;i++)  report2<<" "<<recdev(i);
      report2<<endl;
    }

    if(do_recdev>0)
    {
      report2<<runnumber<<" Recr_main ";
      for(i=recdev_start;i<=recdev_end;i++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Recr_main ";
      for(i=recdev_start;i<=recdev_end;i++) report2<<" "<<recdev(i);
      report2<<endl;
    }

    if(Do_Forecast>0)
    {
      report2<<runnumber<<" Recr_fore ";
      for(i=recdev_end+1;i<=YrMax;i++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Recr_fore ";
      for(i=recdev_end+1;i<=YrMax;i++) report2<<" "<<recdev(i);
      report2<<endl;
      report2<<runnumber<<" Impl_err ";
      for(i=endyr+1;i<=YrMax;i++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Impl_err ";
      for(i=endyr+1;i<=YrMax;i++) report2<<" "<<Fcast_impl_error(i);
      report2<<endl;
    }

    report2<<runnumber<<" init_F ";
    for(i=1;i<=Nfleet;i++) {NP++; report2<<" "<<ParmLabel(NP);}
    report2<<endl<<runnumber<<" init_F ";
    for(i=1;i<=Nfleet;i++) report2<<" "<<init_F(i);
    report2<<endl;

    if(F_Method==2)
    {
      report2<<runnumber<<" F_rate ";
      for(i=1;i<=Nfleet*(TimeMax-styr+1);i++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" F_rate ";
      for(i=1;i<=Nfleet*(TimeMax-styr+1);i++) report2<<" "<<F_rate(i);
      report2<<endl;
    }

    if(Q_Npar>0)
    {
      report2<<runnumber<<" Q_parm ";
      for(i=1;i<=Q_Npar;i++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Q_parm ";
      for(i=1;i<=Q_Npar;i++) report2<<" "<<Q_parm(i);
      report2<<endl;
    }

    if(N_selparm2>0)
    {
      report2<<runnumber<<" Sel_parm ";
      for(i=1;i<=N_selparm2;i++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Sel_parm "<<selparm<<endl;
    }

    if(N_selparm_dev>0)
    {
      report2<<runnumber<<" Sel_parm_dev ";
      for(i=1;i<=N_selparm_dev;i++)
      for(j=selparm_dev_minyr(i);j<=selparm_dev_maxyr(i);j++)
      {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Sel_parm_dev "<<selparm_dev<<endl;
    }

    if(Do_TG>0)
    {
      report2<<runnumber<<" Tag_parm ";
      for (f=1;f<=3*N_TG+2*Nfleet;f++) {NP++; report2<<" "<<ParmLabel(NP);}
      report2<<endl<<runnumber<<" Tag_parm "<<TG_parm<<endl;
    }


  if(Do_CumReport==2)
  {
    if(nobs_cr>0)
    for(f=1;f<=Ntypes;f++)
    if(nyr_cr(f)>0)
    {
     report2 <<runnumber<<" Index:"<<f<<" Year "<<yr_cr_y(f)<<endl;
 //    report2 <<runnumber<<" Index:"<<f<<" VBIO "<<vul_bio(f)<<endl;
     report2 <<runnumber<<" Index:"<<f<<" OBS "<<obs_cr(f)<<endl;
     if(cr_errtype(f)>=0)  // lognormal or lognormal T_dist
     {report2 <<runnumber<<" Index:"<<f<<" EXP "<<elem_prod(mfexp(log_q_cr(f)),mfexp(exp_cr(f)))<<endl;}
     else  // normal error
     {report2 <<runnumber<<" Index:"<<f<<" EXP "<<elem_prod(q_cr(f),exp_cr(f))<<endl;}
    }

    for(f=1;f<=Ntypes;f++)
    if(Nobs_l(f)>0)
    {
     report2 <<runnumber<<" Len:"<<f<<" YR "<<yr_l_y(f)<<endl;
     report2 <<runnumber<<" Len:"<<f<<" effN "<<neff_l(f)<<endl;
    }
    for(f=1;f<=Ntypes;f++)
    if(Nobs_a(f)>0)
    {
     report2 <<runnumber<<" Age:"<<f<<" YR "<<yr_a_y(f)<<endl;
     report2 <<runnumber<<" Age:"<<f<<" effN "<<neff_a(f)<<endl;
    }
    }
    report2<<endl;
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 37 write_rebuilder_output */
FUNCTION void write_rebuilder_output()
  {
    if(rundetail>0 && mceval_counter==0) cout<<" produce output for rebuilding package"<<endl;
    ofstream rebuilder("rebuild.sso",ios::app);
    ofstream rebuild_dat("rebuild.dat");

    if(mceval_counter==0) // writing to rebuild.dat
    {
      rebuild_dat<<"#Title, #runnumber: "<<runnumber<<" "<<datfilename<<" "<<ctlfilename<<
         " "<<obj_fun<<" "<<SPB_yr(styr-2)<<" "<<SPB_yr(endyr+1)<<" StartTime: "<<ctime(&start);
      rebuild_dat<<"SSv3_default_rebuild.dat"<<endl;
      rebuild_dat<<"# Number of sexes"<<endl<<gender<<endl;
      rebuild_dat<<"# Age range to consider (minimum age; maximum age)"<<endl<<0<<" "<<nages<<endl;
      rebuild_dat<<"# Number of fleets"<<endl<<Nfleet<<endl;
      rebuild_dat<<"# First year of projection (Yinit)"<<endl<<Rebuild_Yinit<<endl;
      rebuild_dat<<"# First Year of rebuilding period (Ydecl)"<<endl<<Rebuild_Ydecl<<endl;
      rebuild_dat<<"# Number of simulations"<<endl<<1000<<endl;
      rebuild_dat<<"# Maximum number of years"<<endl<<500<<endl;
      rebuild_dat<<"# Conduct projections with multiple starting values (0=No;else yes)"<<endl<<0<<endl;
      rebuild_dat<<"# Number of parameter vectors"<<endl<<1000<<endl;
      rebuild_dat<<"# Is the maximum age a plus-group (1=Yes;2=No)"<<endl<<1<<endl;
      rebuild_dat<<"# Generate future recruitments using historical recruitments (1)  historical recruits/spawner (2)  or a stock-recruitment (3)"<<endl<<3<<endl;
      rebuild_dat<<"# Constant fishing mortality (1) or constant Catch (2) projections"<<endl<<1<<endl;
      rebuild_dat<<"# Fishing mortality based on SPR (1) or actual rate (2)"<<endl<<1<<endl;
      rebuild_dat<<"# Pre-specify the year of recovery (or -1) to ignore"<<endl<<-1<<endl;
      rebuild_dat<<"# Fecundity-at-age"<<endl;
    }

    // stuff written to both rebuild.dat and rebuild.SSO
//a. "blank line" with run info
    if( mceval_phase() ) rebuilder<<"# mceval phase, cnt="<<mceval_counter<<", StartTime: "<<ctime(&start);
                    else rebuilder<<"# in maximization mode, StartTime: "<<ctime(&start);

    if(mceval_counter==0) rebuild_dat <<"#"<< age_vector << " #runnumber: "<<runnumber<<" "<<datfilename<<" "<<ctlfilename<<" "<<obj_fun<<
      " "<<SPB_yr(styr-2)<<" "<<SPB_yr(endyr+1)<<endl;

//b.  fecundity-at-age
    t = styr+(Rebuild_Yinit-styr)*nseas;
    dvar_vector tempvec2(0,nages);
    dvar_vector tempvec3(0,nages);
    tempvec_a.initialize(); tempvec2.initialize();
    for(p=1;p<=pop;p++)
    for(g=1;g<=gmorph;g++)
    if(sx(g)==1)
    {
//  NEED to adjust for spawning timing within season
      tempvec_a+=elem_prod(fec(g),natage(t+spawn_seas-1,p,g));
      tempvec2+=natage(t+spawn_seas-1,p,g);
    }
    tempvec_a= elem_div(tempvec_a,tempvec2);
    rebuilder<<tempvec_a<<" #female fecundity; weighted by N in year Y_init across morphs and areas"<<endl;
    if(mceval_counter==0)
    {
      rebuild_dat<<tempvec_a<<" #female fecundity; weighted by N in year Y_init across morphs and areas"<<endl;
      rebuild_dat<<"# Age specific selectivity and weight adjusted for discard and discard mortality"<<endl;
    }

//c.  Weight-at-age and selectivity-at-age (ordered by sex and fleet).
    // use the deadfish vectors that account for discard and for mortality of discards
    // average across morphs and areas using N_at_Age in year Yinit and across seasons using Fcast_RelF
    for (gg=1;gg<=gender;gg++)
    for (f=1;f<=Nfleet;f++)
    {
      tempvec_a.initialize();tempvec2.initialize();tempvec3.initialize(); j=0;
      for(s=1;s<=nseas;s++)
      for(p=1;p<=pop;p++)
      if (fleet_area(f)==p && Fcast_RelF_Use(s,f)>0.0)   // active fishery in this area in endyr
      {
        j=1;
        for(g=1;g<=gmorph;g++)
        if(sx(g)==gg)
        {
          tempvec_a+=elem_prod(fish_body_wt(t,g,f),natage(t+s-1,p,g)*Fcast_RelF_Use(s,f));  // body wt
          tempvec2+=elem_prod(sel_al_3(s,g,f),natage(t+s-1,p,g)*Fcast_RelF_Use(s,f));  //no wt
          tempvec3+=natage(t+s-1,p,g)*Fcast_RelF_Use(s,f);
        }
      }
      if(j==1)
      {
        tempvec_a= elem_div(tempvec_a,tempvec3);
        tempvec2 = elem_div(tempvec2,tempvec3);
        rebuilder <<tempvec_a<< " #bodywt for gender,fleet: "<<gg<<" / "<<f<< endl;
        rebuilder <<tempvec2<< " #selex for gender,fleet: "<<gg<<" / "<<f<< endl;
        if(mceval_counter==0)
        {
          rebuild_dat << " #wt and selex for gender,fleet: "<<gg<<" "<<f<< endl;
          rebuild_dat <<tempvec_a<<endl<<tempvec2<< endl;
        }
      }
    }

//d.  Natural mortality and numbers-at-age for year Yinit  (females then males).
    if(mceval_counter==0) rebuild_dat<<"# M and current age-structure in year Yinit: "<<Rebuild_Yinit<<endl;

    for (gg=1;gg<=gender;gg++)
    {
      tempvec_a.initialize();tempvec2.initialize();tempvec3.initialize();
      for(p=1;p<=pop;p++)
      {
        for(g=1;g<=gmorph;g++)
        if(sx(g)==gg && use_morph(g)>0)
        {
          tempvec_a+=elem_prod(natM(1,GP3(g)),natage(t,p,g));  tempvec2+=natage(t,p,g);  //  note, uses season 1 only
        }
      }
      tempvec_a=elem_div(tempvec_a,tempvec2);
      rebuilder<<tempvec_a<<" #mean M for year Yinit: "<<Rebuild_Yinit<<" sex: "<<gg<<endl;
      rebuilder<<tempvec2<<" #numbers for year Yinit: "<<Rebuild_Yinit<<" sex: "<<gg<<endl;
      if(mceval_counter==0) rebuild_dat<<" # gender = "<<gg<<endl<< tempvec_a<<endl<<tempvec2<< endl;
    }

//e.  Numbers-at-age for year Ydecl  (females then males).
    t = styr+(Rebuild_Ydecl-styr)*nseas;
    if(mceval_counter==0) rebuild_dat<<"# Age-structure at Ydeclare= "<<Rebuild_Ydecl<<endl;
    for (gg=1;gg<=gender;gg++)
    {
      tempvec_a.initialize();tempvec2.initialize();tempvec3.initialize();
      for(p=1;p<=pop;p++)
      {
        for(g=1;g<=gmorph;g++)
        if(sx(g)==gg && use_morph(g)>0)
        {
          tempvec2+=natage(t,p,g);
        }
      }
      rebuilder <<tempvec2<<" #numbers for year Ydeclare: "<<Rebuild_Ydecl<<" sex: "<<gg<< endl;
      if(mceval_counter==0) rebuild_dat<<tempvec2<< endl;
    }

    k=endyr;
    if(Rebuild_Yinit>k) k=Rebuild_Yinit;

//f. "blank line" used for header for following lines
    rebuilder <<"#R0 "<<years<< " #years"<<endl;

//g. recruitment
    rebuilder << exp_rec(styr-2,4)<<" ";
    for(y=styr;y<=k;y++) {rebuilder<<exp_rec(k,4)<<" "; }
    rebuilder << " #Recruits"<<endl;

//h. spawnbio
    rebuilder << SPB_yr(styr-2)<<" "<<SPB_yr(styr,k) <<" #SpawnBio"<<endl;

//i. steepness; SigmaR; rho
    rebuilder << SR_parm(2) <<" "<< sigmaR <<" "<< SR_parm(N_SRparm2) <<" # spawn-recr steepness, sigmaR, autocorr"<< endl;

    if(mceval_counter==0)
    {
      rebuild_dat<<"# Year for Tmin Age-structure (set to Ydecl by SS)"<<endl<<Rebuild_Ydecl<<endl;

      rebuild_dat<<"#  recruitment and biomass"<<endl<<"# Number of historical assessment years"<<endl<<k-styr+2<<endl;
      rebuild_dat<<"# Historical data"<<endl<<"# year recruitment spawner in B0 in R project in R/S project"<<endl;
      rebuild_dat<<styr-1<<" "<<years;
      if(Rebuild_Yinit>endyr) rebuild_dat<<" "<<Rebuild_Yinit;
      rebuild_dat<<" #years (with first value representing R0)"<<endl;
      rebuild_dat<<exp_rec(styr-2,4)<<" ";
      for(y=styr;y<=k;y++) {rebuild_dat<<exp_rec(y,4)<<" ";}
      rebuild_dat<<" #recruits; first value is R0 (virgin)"<< endl;
      rebuild_dat<<SPB_yr(styr-2)<<" "<<SPB_yr(styr,k) <<" #spbio; first value is SPB_virgin (virgin)"<< endl;
      rebuild_dat<<1<<" ";
      for(y=styr;y<=k;y++) rebuild_dat<<0<<" ";
      rebuild_dat<<" # in Bzero"<<endl;
      rebuild_dat<<0<<" ";
      for(y=styr;y<=k-3;y++) rebuild_dat<<1<<" ";
      rebuild_dat<<" 0 0 0 # in R project"<<endl;
      rebuild_dat<<0<<" ";
      for(y=styr;y<=k-3;y++) rebuild_dat<<1<<" ";
      rebuild_dat<<" 0 0 0 # in R/S project"<<endl;
      rebuild_dat<<"# Number of years with pre-specified catches"<<endl<<0<<endl;
      rebuild_dat<<"# catches for years with pre-specified catches go next"<<endl;
//      rebuild_dat<<"# Number of future recruitments to override"<<endl<<0<<endl;
//      rebuild_dat<<"# Process for overiding (-1 for average otherwise index in data list)"<<endl;

      rebuild_dat<<"# Number of future recruitments to override"<<endl;
      rebuild_dat<<Rebuild_Yinit - Rebuild_Ydecl<<endl;
      rebuild_dat<<"# Process for overiding (-1 for average otherwise index in data list)"<<endl;
      if(Rebuild_Yinit>=Rebuild_Ydecl+1)
      {
        for(y=Rebuild_Ydecl+1;y<=Rebuild_Yinit;y++) rebuild_dat<<y<<" "<<1<<" "<<y<<endl;
      }

      rebuild_dat<<"# Which probability to product detailed results for (1=0.5; 2=0.6; etc.)"<<endl<<3<<endl;
      rebuild_dat<<"# Steepness sigma-R Auto-correlation"<<endl<<SR_parm(2) <<" "<< sigmaR <<" "<< 0 << endl;
      rebuild_dat<<"# Target SPR rate (FMSY Proxy); manually change to SPR_MSY if not using SPR_target"<<endl<<SPR_target<<endl;
      rebuild_dat<<"# Discount rate (for cumulative catch)"<<endl<<0.1<<endl;
      rebuild_dat<<"# Truncate the series when 0.4B0 is reached (1=Yes)"<<endl<<0<<endl;
      rebuild_dat<<"# Set F to FMSY once 0.4B0 is reached (1=Yes)"<<endl<<0<<endl;
      rebuild_dat<<"# Maximum possible F for projection (-1 to set to FMSY)"<<endl<<-1<<endl;
      rebuild_dat<<"# Defintion of recovery (1=now only;2=now or before)"<<endl<<2<<endl;
      rebuild_dat<<"# Projection type"<<endl<<4<<endl;
      rebuild_dat<<"# Definition of the 40-10 rule"<<endl<<10<<" "<<40<<endl;
      rebuild_dat<<"# Calculate coefficients of variation (1=Yes)"<<endl<<0<<endl;
      rebuild_dat<<"# Number of replicates to use"<<endl<<10<<endl;
      rebuild_dat<<"# Random number seed"<<endl<<-99004<<endl;
      rebuild_dat<<"# File with multiple parameter vectors "<<endl<<"rebuild.SSO"<<endl;
      rebuild_dat<<"# User-specific projection (1=Yes); Output replaced (1->9)"<<endl<<"0  5"<<endl;
      rebuild_dat<<"# Catches and Fs (Year; 1/2/3 (F or C or SPR); value); Final row is -1"<<endl;
      rebuild_dat<< k <<" 1 1"<<endl<<"-1 -1 -1"<<endl;
      rebuild_dat<<"# Fixed catch project (1=Yes); Output replaced (1->9); Approach (-1=Read in else 1-9)"<<endl;
      rebuild_dat<<"0 2 -1"<<endl;
      tempvec_a(1,Nfleet)=colsum(Fcast_RelF_Use);
      rebuild_dat<<"# Split of Fs"<<endl;
      rebuild_dat<<Rebuild_Yinit<<" "<<tempvec_a(1,Nfleet)<<endl;
      rebuild_dat<<"-1 ";
      for(f=1;f<=Nfleet;f++) rebuild_dat<<" 1";
      rebuild_dat<<endl;
      rebuild_dat<<"# Yrs to define T_target for projection type 4 (a.k.a. 5 pre-specified inputs)"<<endl;
      rebuild_dat<<endyr+6<<" "<<endyr+7<<" "<<endyr+8<<" "<<endyr+9<<" "<<endyr+10<<" "<<endl;
      rebuild_dat<<"# Year for probability of recovery"<<endl;
      rebuild_dat<<endyr+10<<" "<<endyr+11<<" "<<endyr+12<<" "<<endyr+13<<" "<<endyr+14<<" "<<endyr+15<<" "<<endyr+16<<" "<<endyr+17<<endl;
      rebuild_dat<<"# Time varying weight-at-age (1=Yes;0=No)"<<endl<<0<<endl;
      rebuild_dat<<"# File with time series of weight-at-age data"<<endl<<"none"<<endl;
      rebuild_dat<<"# Use bisection (0) or linear interpolation (1)"<<endl<<1<<endl;
      rebuild_dat<<"# Target Depletion"<<endl<<0.4<<endl;
      rebuild_dat<<"# CV of implementation error"<<endl<<0<<endl;
    }
//  end output of rebuilding quantities
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 38 write_nudata */
FUNCTION void write_nudata()
  {
//  code for multinomial distribution developed by Ian Stewart, Oct 2005
  random_number_generator radm(long(time(&start)));

  dvector temp_mult(1,1000);
  dvector temp_probs(1,nlen_bin2);
  int compindex;
  dvector temp_probs2(1,n_abins2);
  int Nudat;
//  create bootstrap data files; except first file just replicates the input and second is the estimate without error
  for(i=1;i<=1234;i++) temp = randn(radm);
  cout << " N_nudata: " << N_nudata << endl;
  ofstream report1("data.ss_new");
  report1<<version_info_short<<endl;
  report1<<"#_"<<version_info<<endl<<"#_Start_time: "<<ctime(&start);
  report1  << "#_Number_of_datafiles: " << N_nudata << endl;
  for(Nudat=1;Nudat<=N_nudata;Nudat++)
  {
  if(Nudat==1)
  {
    report1<<Data_Comments<<endl;
    report1 << "#_observed data: "<< endl;
  }
  else if(Nudat==2)
  {report1 << "#_expected values with no error added " << endl;}
  else
  {report1 << "#_bootstrap file: " << Nudat-2 << endl;}

  report1 << styr << " #_styr"<<endl;
  report1 << endyr <<" #_endyr"<< endl;
  report1 << nseas <<" #_nseas"<< endl;
  report1 << 12.*seasdur<<" #_months/season"<< endl;

  report1 << spawn_seas <<" #_spawn_seas"<< endl;
  report1 << Nfleet<<" #_Nfleet"<< endl;
  report1 << Nsurvey<<" #_Nsurveys"<< endl;
  report1 << pop<<" #_N_areas"<<endl;
  report1 << fleetnameread<<endl;
  report1 << surveytime<<" #_surveytiming_in_season"<<endl;
  report1 << fleet_area<<" #_area_assignments_for_each_fishery_and_survey"<<endl;
  report1 << catchunits<<" #_units of catch:  1=bio; 2=num"<<endl;
  report1 << catch_se_rd<<" #_se of log(catch) only used for init_eq_catch and for Fmethod 2 and 3; use -1 for discard only fleets"<<endl;
  report1 << gender<<" #_Ngenders"<< endl;
  report1 << nages<<" #_Nages"<< endl;

  if(Nudat==1)  // report back the input data
  {

  report1 << obs_equ_catch<<" #_init_equil_catch_for_each_fishery"<<endl;
  report1 << N_ReadCatch<<" #_N_lines_of_catch_to_read"<<endl;
  report1 << "#_catch_biomass(mtons):_columns_are_fisheries,year,season"<<endl;
  report1 << catch_bioT<< endl << endl;
  report1<<"#"<<endl;

  report1 << nobs_cr_rd <<" #_N_cpue_and_surveyabundance_observations"<< endl;
  report1<<"#_Units:  0=numbers; 1=biomass; 2=F"<<endl;
  report1<<"#_Errtype:  -1=normal; 0=lognormal; >0=T"<<endl;
  report1<<"#_Fleet Units Errtype"<<endl;
  for(f=1;f<=Ntypes;f++) report1<<f<<" "<<cr_units(f)<<" "<<cr_errtype(f)<<" # "<<fleetname(f)<<endl;
  report1 << "#_year seas index obs err"<<endl;
  if(nobs_cr_rd>0)
  for (i=1;i<=nobs_cr_rd;i++)
  {
    f=abs(indexdata(i,3));
    report1 << indexdata(i)<<" # "<<fleetname(f)<<endl;
  }


  report1<<endl<<Ndisc_fleets<<" #_N_fleets_with_discard"<<endl;
  if(Ndisc_fleets>0)
  {
    report1<<"#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)"<< endl;
    report1<<"#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal"<<endl;
    report1<<"#_Fleet units errtype"<<endl;
    for(f=1;f<=Ntypes;f++)
    if(disc_units(f)>0) report1<<f<<" "<<disc_units(f)<<" "<<disc_errtype(f)<<" # "<<fleetname(f)<<endl;
    report1<<nobs_disc_rd<<" #_N_discard_obs"<< endl;
    report1<<"#_year seas fleet obs err"<<endl;
    for (i=1;i<=nobs_disc_rd;i++)
    {
      f=abs(discdata(i,3));
      report1 << discdata(i)<<" # "<<fleetname(f)<<endl;
    }
  }
  else
  {
    report1<<"#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)"<< endl;
    report1<<"#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal"<<endl;
    report1<<"#Fleet Disc_units err_type"<<endl;
    report1<<"0 #N discard obs"<<endl;
    report1<<"#_year seas index obs err"<<endl;
  }

  report1 <<"#"<<endl<< nobs_mnwt_rd <<" #_N_meanbodywt_obs"<< endl;
  report1<<DF_bodywt<<" #_DF_for_meanbodywt_T-distribution_like"<<endl;
  if(nobs_mnwt_rd>0)
   {
    report1 << "#_year seas index part obs cv"<<endl;
    for (i=1;i<=nobs_mnwt_rd;i++)
    {
      f=abs(mnwtdata1(i,3));
      report1 << mnwtdata1(i)<<" # "<<fleetname(f)<<endl;
    }
   }

  report1<<endl<<LenBin_option<<" # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector"<<endl;
  if(LenBin_option==1)
  {report1<<"# no additional input for option 1"<<endl;
    report1<<"# read binwidth, minsize, lastbin size for option 2"<<endl;
    report1<<"# read N poplen bins, then vector of bin lower boundaries, for option 3"<<endl;}
  else if(LenBin_option==2)
  {
    report1<<binwidth2<<" # binwidth for population size comp "<<endl;
    report1<<minLread<<" # minimum size in the population (lower edge of first bin and size at age 0.00) "<<endl;
    report1<<maxLread<<" # maximum size in the population (lower edge of last bin) "<<endl;
  }
  else
  {
    report1<<nlength<<" # number of population size bins "<<endl;
    report1<<len_bins<<endl;
  }

  report1<<endl<<min_tail<<" #_comp_tail_compression"<<endl;
  report1<<min_comp<<" #_add_to_comp"<<endl;
  report1<<CombGender_l<<" #_combine males into females at or below this bin number"<<endl;
  report1<<nlen_bin<<" #_N_LengthBins"<<endl<<len_bins_dat<<endl;

  report1<<nobsl_rd<<" #_N_Length_obs"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Nsamp datavector(female-male)"<<endl;
  if(nobsl_rd>0) report1<<lendata<<endl;

   report1 <<endl<<n_abins<<" #_N_age_bins"<<endl;
  if(n_abins>0) report1<<age_bins1<<endl;
  report1 << N_ageerr <<" #_N_ageerror_definitions"<< endl;
  if(N_ageerr>0) report1 << age_err_rd << endl;

  report1 <<nobsa_rd<<" #_N_Agecomp_obs"<<endl;
  report1<<Lbin_method<<" #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths"<<endl;
  report1<<CombGender_a<<" #_combine males into females at or below this bin number"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)"<<endl;
   if(nobsa_rd>0) report1<<agedata<<endl;

  report1 <<endl<<nobs_ms_rd<<" #_N_MeanSize-at-Age_obs"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Ageerr Ignore datavector(female-male)"<<endl;
  report1<<"#                                          samplesize(female-male)"<<endl;
   if(nobs_ms_rd>0) report1<<sizeagedata<<endl;

    report1 << endl << N_envvar<<" #_N_environ_variables"<<endl<<N_envdata<<" #_N_environ_obs"<<endl;
    if(N_envdata>0) report1<<env_temp<<endl;

  report1<<SzFreq_Nmeth<<" # N sizefreq methods to read "<<endl;
  if(SzFreq_Nmeth>0)
  {
    report1<<SzFreq_Nbins<<" #Sizefreq N bins per method"<<endl;
    report1<<SzFreq_units<<" #Sizetfreq units(bio/num) per method"<<endl;
    report1<<SzFreq_scale<<" #Sizefreq scale(kg/lbs/cm/inches) per method"<<endl;
    report1<<SzFreq_mincomp<<" #Sizefreq mincomp per method "<<endl;
    report1<<SzFreq_nobs<<" #Sizefreq N obs per method"<<endl;
    report1<<"#_Sizefreq bins "<<endl;
    for(i=1;i<=SzFreq_Nmeth;i++) {report1<<SzFreq_Omit_Small(i)*SzFreq_bins1(i,1)<<SzFreq_bins1(i)(2,SzFreq_Nbins(i))<<endl;}
    report1<<"#_Year season Fleet Gender Partition SampleSize <data> "<<endl<<SzFreq_obs1<<endl;
  }

  // tagging data (as observed)
  report1<<endl<<Do_TG<<" # do tags (0/1)"<<endl;
  if(Do_TG>0)
  {
    // info on dimensions of tagging data
    report1<<N_TG<<" # N tag groups"<<endl;
    report1<<N_TG_recap<<" # N recap events"<<endl;
    report1<<TG_mixperiod<<" # mixing latency period: N periods to delay before comparing observed to expected recoveries (0 = release period)"<<endl;
    report1<<TG_maxperiods<<" # max periods (seasons) to track recoveries, after which tags enter accumulator"<<endl;

    // tag releases
    report1<<"# Release data for each tag group.  Tags are considered to be released at the beginning of a season (period)"<<endl;
    report1<<"#<TG> area yr season <tfill> gender age Nrelease  (note that the TG and tfill values are placeholders and are replaced by program generated values)"<<endl;
    report1<<TG_release<<endl;

    // tag recaptures
    report1<<"#_TAG  Year Season Fleet Nrecap"<<endl;
    for(j=1;j<=N_TG_recap;j++)
    {
      // fill in first 4 columns:
      for(k=1;k<=5;k++) report1<<TG_recap_data(j,k)<<" ";
      report1<<endl;
    }
  }
  // end tagging data section #1 (observed data)
  
    report1<<endl<<0<<" # no morphcomp data "<<endl;
    report1 << endl << "999" << endl << endl;

  }

  else if(Nudat==2)  // report expected value with no added error
  {

  report1 << est_equ_catch<<" #_init_equil_catch_for_each_fishery"<<endl;
  report1 << (endyr-styr+1)*nseas<<" #_N_lines_of_catch_to_read"<<endl;
  report1 << "#_catch_biomass(mtons):_columns_are_fisheries,year,season"<<endl;
  for(y=styr; y<=endyr; y++)
  for(s=1; s<=nseas;s++)
  {
    t=styr+(y-styr)*nseas+s-1;
    for(f=1;f<=Nfleet;f++)
    {
      if(catch_se_rd(f)<0.0 && catch_ret_obs(f,t)>0.0)
        {
          report1<<" 0.1 ";  //  for bycatch only fleet
        }
      else if(catchunits(f)==1)
      {report1<<catch_fleet(t,f,3)<<" ";}
      else
      {report1<<catch_fleet(t,f,6)<<" ";}
    }
    report1<<y<<" "<<s<<endl;
  }
  report1<<"#"<<endl;

    report1 << nobs_cr <<" #_N_cpue_and_surveyabundance_observations"<< endl;
    report1<<"#_Units:  0=numbers; 1=biomass; 2=F"<<endl;
    report1<<"#_Errtype:  -1=normal; 0=lognormal; >0=T"<<endl;
    report1<<"#_Fleet Units Errtype"<<endl;
    for(f=1;f<=Ntypes;f++) report1<<f<<" "<<cr_units(f)<<" "<<cr_errtype(f)<<endl;
    report1 << "#_year seas index obs err"<<endl;
    if(nobs_cr>0)
    for(f=1;f<=Ntypes;f++)
    for (i=1;i<=nyr_cr(f);i++)
    {
      report1 << yr_cr_y(f,i)<<" "<<yr_cr_s(f,i)<<" "<<f*yr_cr_use(f,i)<<" ";
      if(yr_cr_use(f,i)>0)
      {
        if(cr_errtype(f)>=0)  // lognormal
        {
          report1 << mfexp(exp_cr(f,i)+log_q_cr(f,i));
        }
        else if(cr_errtype(f)==-1)  // normal
        {
          report1<<exp_cr(f,i)*q_cr(f,i);
        }
      }
      else
      {
        report1 << obs_cr(f,i);
      }
      report1 <<" "<<se_cr_obs(f,i)<<" #_orig_obs: "<<obs_cr(f,i)<<" "<<fleetname(f)<<endl;
    }

  report1<<endl<<Ndisc_fleets<<" #_N_fleets_with_discard"<<endl;
  if(Ndisc_fleets>0)
  {
    report1<<"#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)"<< endl;
    report1<<"#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal"<<endl;
    report1<<"#_Fleet units errtype"<<endl;
    for(f=1;f<=Ntypes;f++)
    if(disc_units(f)>0) report1<<f<<" "<<disc_units(f)<<" "<<disc_errtype(f)<<" # "<<fleetname(f)<<endl;
    report1<<nobs_disc<<" #_N_discard_obs"<< endl;
    report1<<"#_year seas fleet obs err"<<endl;
    for(f=1;f<=Ntypes;f++)
    if(nyr_disc(f)>0)
    for (i=1;i<=nyr_disc(f);i++)
    {
      s=abs(yr_disc_s(f,i));
      report1 << yr_disc_y(f,i)<<" "<<yr_disc_s(f,i)<<" "<<f*yr_disc_use(f,i)<<" ";
      if(yr_disc_use(f,i) >= 0.0 )
        {report1 << exp_disc(f,i);}
      else
      {report1 << obs_disc(f,i);}
      report1 << " "<< cv_disc(f,i)<<" #_orig_obs: "<<obs_disc(f,i)<<endl;
    }
  }
  else
  {
    report1<<"#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)"<< endl;
    report1<<"#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal"<<endl;
    report1<<"#Fleet Disc_units err_type"<<endl;
    report1<<"0 #N discard obs"<<endl;
    report1<<"#_year seas index obs err"<<endl;
  }

  report1 <<"#"<<endl<< nobs_mnwt <<" #_N_meanbodywt_obs"<< endl;
  report1<<DF_bodywt<<" #_DF_for_meanbodywt_T-distribution_like"<<endl;

  if(nobs_mnwt>0)
   {report1 << "#_year seas index part obs cv"<<endl;
   for(i=1;i<=nobs_mnwt;i++)
    {
     report1 << mnwtdata(1,i)<<" "<<mnwtdata(2,i)<<" "<<mnwtdata(3,i)<<" "<<mnwtdata(4,i)<<" "<<
     exp_mnwt(i)<<" "<<mnwtdata(6,i)<<" #_orig_obs: "<<mnwtdata(5,i)<< endl;
    }
   }

  report1<<endl<<LenBin_option<<" # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector"<<endl;
  if(LenBin_option==1)
  {report1<<"# no additional input for option 1"<<endl;
    report1<<"# read binwidth, minsize, lastbin size for option 2"<<endl;
    report1<<"# read N poplen bins, then vector of bin lower boundaries, for option 3"<<endl;}
  else if(LenBin_option==2)
  {
    report1<<binwidth2<<" # binwidth for population size comp "<<endl;
    report1<<minLread<<" # minimum size in the population (lower edge of first bin and size at age 0.00) "<<endl;
    report1<<maxLread<<" # maximum size in the population (lower edge of last bin) "<<endl;
  }
  else
  {
    report1<<nlength<<" # number of population size bins "<<endl;
    report1<<len_bins<<endl;
  }

  report1<<endl<<min_tail<<" #_comp_tail_compression"<<endl;
  report1<<min_comp<<" #_add_to_comp"<<endl;
  report1<<CombGender_l<<" #_combine males into females at or below this bin number"<<endl;
  report1<<nlen_bin<<" #_N_LengthBins"<<endl<<len_bins_dat<<endl;
  report1<<sum(Nobs_l)<<" #_N_Length_obs"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Nsamp datavector(female-male)"<<endl;
   for(f=1;f<=Ntypes;f++)
    {
    if(Nobs_l(f)>0)
    {
     for(i=1;i<=Nobs_l(f);i++)
     {
      if(nsamp_l(f,i)>0 && header_l(f,i,3)) // do only if this was a real observation
      {
       k=1000;  if(nsamp_l(f,i)<k) k=nsamp_l(f,i);
       exp_l_temp_dat = nsamp_l(f,i)*value(exp_l(f,i)/sum(exp_l(f,i)));
      }
      else
      {exp_l_temp_dat = obs_l(f,i);}
     report1 << header_l(f,i)(1,3)<<" "<<gen_l(f,i)<<" "<<mkt_l(f,i)<<" "<<nsamp_l(f,i)<<" "<<exp_l_temp_dat<<endl;
    }}}

   report1 <<endl<<n_abins<<" #_N_age_bins"<<endl;
  if(n_abins>0) report1<<age_bins1<<endl;
  report1 << N_ageerr <<" #_N_ageerror_definitions"<< endl;
  if(N_ageerr>0) report1 << age_err_rd << endl;

  report1 <<Nobs_a_tot<<" #_N_Agecomp_obs"<<endl;
  report1<<Lbin_method<<" #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths"<<endl;
  report1<<CombGender_a<<" #_combine males into females at or below this bin number"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)"<<endl;
   if(Nobs_a_tot>0)
   for(f=1;f<=Ntypes;f++)
    {
    if(Nobs_a(f)>=1)
    {
     for(i=1;i<=Nobs_a(f);i++)
     {
     if(nsamp_a(f,i)>0 && header_a(f,i,3)>0) // if real observation
     {
      k=1000;  if(nsamp_a(f,i)<k) k=nsamp_a(f,i);  // note that nsamp is adjusted by var_adjust, so var_adjust
                                                   // should be reset to 1.0 in control files that read the nudata.dat files
      exp_a_temp = nsamp_a(f,i)*value(exp_a(f,i)/sum(exp_a(f,i)));
     }
     else
     {exp_a_temp = obs_a(f,i);}
    report1 << header_a(f,i)(1,8)<<" "<<nsamp_a(f,i)<<" "<<exp_a_temp<<endl;
    }
    }
    }
  report1 <<endl<<nobs_ms_tot<<" #_N_MeanSize-at-Age_obs"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Ageerr Ignore datavector(female-male)"<<endl;
  report1<<"#                                          samplesize(female-male)"<<endl;
   if(nobs_ms_tot>0)
   for(f=1;f<=Ntypes;f++)
    {
    if(Nobs_ms(f)>0)
    {
     for(i=1;i<=Nobs_ms(f);i++)
     {
     report1 << header_ms(f,i)(1,7);
     for(a=1;a<=n_abins2;a++)
     {
     report1 << " " ;
         if(obs_ms_n(f,i,a)>0)
          {temp=exp_ms(f,i,a);
          if(temp<=0.) {temp=0.0001;}
          report1 << temp;}
         else
             {report1 << obs_ms(f,i,a) ;}
         }

     report1 << endl<< elem_prod(obs_ms_n(f,i),obs_ms_n(f,i)) << endl;
      }
      }
      }
    report1 << endl << N_envvar<<" #_N_environ_variables"<<endl<<N_envdata<<" #_N_environ_obs"<<endl;
    if(N_envdata>0) report1<<env_temp<<endl;

  report1<<SzFreq_Nmeth<<" # N sizefreq methods to read "<<endl;
  if(SzFreq_Nmeth>0)
  {
    report1<<SzFreq_Nbins<<" #Sizefreq N bins per method"<<endl;
    report1<<SzFreq_units<<" #Sizetfreq units(bio/num) per method"<<endl;
    report1<<SzFreq_scale<<" #Sizefreq scale(kg/lbs/cm/inches) per method"<<endl;
    report1<<SzFreq_mincomp<<" #Sizefreq mincomp per method "<<endl;
    report1<<SzFreq_nobs<<" #Sizefreq N obs per method"<<endl;
    report1<<"#_Sizefreq bins "<<endl;
    for(i=1;i<=SzFreq_Nmeth;i++) {report1<<SzFreq_Omit_Small(i)*SzFreq_bins1(i,1)<<SzFreq_bins1(i)(2,SzFreq_Nbins(i))<<endl;}
    report1<<"#_Year season Fleet Gender Partition SampleSize <data> "<<endl;
    for(iobs=1;iobs<=SzFreq_totobs;iobs++)
    {
      if(SzFreq_obs_hdr(iobs,3)>0)  // flag for date range in bounds
      {
        report1<<SzFreq_obs1(iobs)(1,7)<<" "<<SzFreq_exp(iobs)<<endl;
      }
      else
      {
        report1<<SzFreq_obs1(iobs)<<endl;
      }
    }
  }

  // begin tagging data section #2 (expected values)
  report1<<endl<<Do_TG<<" # do tags (0/1)"<<endl;
  if(Do_TG>0)
  {
    // info on dimensions of tagging data
    report1<<N_TG<<" # N tag groups"<<endl;
    report1<<N_TG_recap<<" # N recap events"<<endl;
    report1<<TG_mixperiod<<" # mixing latency period: N periods to delay before comparing observed to expected recoveries (0 = release period)"<<endl;
    report1<<TG_maxperiods<<" # max periods (seasons) to track recoveries, after which tags enter accumulator"<<endl;

    // tag releases
    report1<<"# Release data for each tag group.  Tags are considered to be released at the beginning of a season (period)"<<endl;
    report1<<"#<TG> area yr season <tfill> gender age Nrelease  (note that the TG and tfill values are placeholders and are replaced by program generated values)"<<endl;
    report1<<TG_release<<endl;

    // tag recaptures
    report1<<"#_Note: Expected values for tag recaptures are reported only for the same combinations of"<<endl; 
    report1<<"#       group, year, area, and fleet that had observed recaptures. "<<endl;
    report1<<"#_TAG  Year Season Fleet Nrecap"<<endl;
    for(j=1;j<=N_TG_recap;j++)
    {
      // fill in first 4 columns:
      for(k=1;k<=4;k++) report1<<TG_recap_data(j,k)<<" ";
      // fill in 5th column with bootstrap values
      TG=TG_recap_data(j,1);
      overdisp=TG_parm(2*N_TG+TG);
      t=styr+int((TG_recap_data(j,2)-styr)*nseas+TG_recap_data(j,3)-1) - TG_release(TG,5); // find elapsed time in terms of number of seasons
      if(t>TG_maxperiods) t=TG_maxperiods;
      report1<<value(TG_recap_exp(TG,t,0))<<" #_overdisp: "<<value(overdisp)<<endl;
    }
  }
  // end tagging data section #2 (expected values)

    report1<<endl<<0<<" # no morphcomp data "<<endl;
    report1 << endl << "999" << endl << endl;

  }

  else  //  create bootstrap data
  {

    for(f=1;f<=Nfleet;f++)
    {
      if(catchunits(f)==1)
      {report1<<est_equ_catch(f)*mfexp(randn(radm)*catch_se(styr-1,f) - 0.5*catch_se(styr-1,f)*catch_se(styr-1,f))<<" ";}
      else
      {report1<<est_equ_catch(f)*mfexp(randn(radm)*catch_se(styr-1,f) - 0.5*catch_se(styr-1,f)*catch_se(styr-1,f))<<" ";}
    }
  report1 <<" #_init_equil_catch_for_each_fishery"<<endl;
  report1 <<(endyr-styr+1)*nseas<<" #_N_lines_of_catch_to_read"<<endl;
  report1 << "#_catch_biomass(mtons):_columns_are_fisheries,year,season"<<endl;
  for(y=styr; y<=endyr; y++)
  for(s=1; s<=nseas;s++)
  {
    t=styr+(y-styr)*nseas+s-1;
    for(f=1;f<=Nfleet;f++)
    {
      if(catch_se_rd(f)<0.0 && catch_ret_obs(f,t)>0.0)
        {
          report1<<" 0.1 ";  //  for bycatch only fleet
        }
      else if(catchunits(f)==1)
      {report1<<catch_fleet(t,f,3)*mfexp(randn(radm)*catch_se(t,f) - 0.5*catch_se(t,f)*catch_se(t,f))<<" ";}
      else
      {report1<<catch_fleet(t,f,6)*mfexp(randn(radm)*catch_se(t,f) - 0.5*catch_se(t,f)*catch_se(t,f))<<" ";}
    }
    report1<<y<<" "<<s<<endl;
  }
  report1<<"#"<<endl;

  report1 << nobs_cr <<" #_N_cpue_and_surveyabundance_observations"<< endl;
  report1<<"#_Units:  0=numbers; 1=biomass; 2=F"<<endl;
  report1<<"#_Errtype:  -1=normal; 0=lognormal; >0=T"<<endl;
  report1<<"#_Fleet Units Errtype"<<endl;
  for(f=1;f<=Ntypes;f++) report1<<f<<" "<<cr_units(f)<<" "<<cr_errtype(f)<<endl;
  report1 << "#_year seas index obs err"<<endl;
  if(nobs_cr>0)
  for(f=1;f<=Ntypes;f++)
  for (i=1;i<=nyr_cr(f);i++)
  {
    report1 << yr_cr_y(f,i)<<" "<<yr_cr_s(f,i)<<" "<<f*yr_cr_use(f,i)<<" ";
    if(yr_cr_use(f,i) > 0)
    {
      if(cr_errtype(f)==-1)  // normal error
      {
        report1<<exp_cr(f,i)*q_cr(f,i)+randn(radm)*se_cr_use(f,i);    //  uses se_cr_use, not se_cr_obs to include both effect of input var_adjust and extra_sd
      }
      if(cr_errtype(f)==0)  // lognormal
      {
         report1 << mfexp(exp_cr(f,i)+log_q_cr(f,i)+ randn(radm)*se_cr_use(f,i) );    //  uses se_cr_use, not se_cr_obs to include both effect of input var_adjust and extra_sd
      }
      else if(cr_errtype(f)>0)   // lognormal T_dist
      {
        temp = sqrt( (cr_errtype(f)+1.)/cr_errtype(f));  // where df=cr_errtype(f)
        report1 << mfexp(exp_cr(f,i)+log_q_cr(f,i)+ randn(radm)*se_cr_use(f,i)*temp );    //  adjusts the sd by the df sample size
      }
    }
    else
    {
      report1 << obs_cr(f,i);
    }
    report1 <<" "<<se_cr_obs(f,i)<<" #_orig_obs: "<<obs_cr(f,i)<<" "<<fleetname(f)<<endl;
  }

  report1<<endl<<Ndisc_fleets<<" #_N_fleets_with_discard"<<endl;
  if(Ndisc_fleets>0)
  {
    report1<<"#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)"<< endl;
    report1<<"#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal"<<endl;
    report1<<"#_Fleet units errtype"<<endl;
    for(f=1;f<=Ntypes;f++)
    if(disc_units(f)>0) report1<<f<<" "<<disc_units(f)<<" "<<disc_errtype(f)<<" # "<<fleetname(f)<<endl;
    report1<<nobs_disc<<" #_N_discard_obs"<< endl;
    report1<<"#_year seas fleet obs err"<<endl;
    for(f=1;f<=Ntypes;f++)
    if(nyr_disc(f)>0)
    for (i=1;i<=nyr_disc(f);i++)
    {
      s=abs(yr_disc_s(f,i));
      report1 << yr_disc_y(f,i)<<" "<<yr_disc_s(f,i)<<" "<<f*yr_disc_use(f,i)<<" ";
      if(yr_disc_use(f,i) >= 0.0 )
      {
        if(disc_errtype(f)>=1)
        {temp=exp_disc(f,i) + randn(radm)*sd_disc(f,i)*sqrt((disc_errtype(f)+1.)/disc_errtype(f)) * exp_disc(f,i); if(temp<0.001) temp=0.001;}
        else if(disc_errtype(f)==0)
        {temp=exp_disc(f,i) + randn(radm)*sd_disc(f,i); if(temp<0.001) temp=0.001; }
        else if(disc_errtype(f)==-1)
        {temp=exp_disc(f,i) + randn(radm)*sd_disc(f,i); if(temp<0.001) temp=0.001; }
        else if(disc_errtype(f)==-2)
        {temp=exp_disc(f,i) * mfexp(randn(radm)*sd_disc(f,i));}
      }
      else
      {temp=obs_disc(f,i);}
      report1 <<" "<<temp<< " "<< cv_disc(f,i)<<" # "<<obs_disc(f,i)<<" "<<sd_disc(f,i)<<endl;
    }
  }
  else
  {
    report1<<"#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)"<< endl;
    report1<<"#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal"<<endl;
    report1<<"#Fleet Disc_units err_type"<<endl;
    report1<<"0 #N discard obs"<<endl;
    report1<<"#_year seas index obs err"<<endl;
  }

  report1 <<"#"<<endl<< nobs_mnwt <<" #_N_meanbodywt_obs"<< endl;
  report1<<DF_bodywt<<" #_DF_for_meanbodywt_T-distribution_like"<<endl;
  if(nobs_mnwt>0)
  {
    report1 << "#_year seas index part obs cv"<<endl;
    for(i=1;i<=nobs_mnwt;i++)
    {
      if(mnwtdata(3,i)>0 && mnwtdata(5,i)>0.)
      {
        temp=exp_mnwt(i)+randn(radm)*mnwtdata(6,i)*sqrt((DF_bodywt+1.)/DF_bodywt) *exp_mnwt(i);
        if(temp<=0.0) {temp=0.0001;}
      }
      else
      {
        temp=mnwtdata(5,i);
      }
      report1 << mnwtdata(1,i)<<" "<<mnwtdata(2,i)<<" "<<mnwtdata(3,i)<<" "<<mnwtdata(4,i)<<" "<<
      temp<<" "<<mnwtdata(6,i)<<" #_orig_obs: "<<mnwtdata(5,i)<< endl;
    }
  }

  report1<<endl<<LenBin_option<<" # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector"<<endl;
  if(LenBin_option==1)
  {report1<<"# no additional input for option 1"<<endl;
    report1<<"# read binwidth, minsize, lastbin size for option 2"<<endl;
    report1<<"# read N poplen bins, then vector of bin lower boundaries, for option 3"<<endl;}
  else if(LenBin_option==2)
  {
    report1<<binwidth2<<" # binwidth for population size comp "<<endl;
    report1<<minLread<<" # minimum size in the population (lower edge of first bin and size at age 0.00) "<<endl;
    report1<<maxLread<<" # maximum size in the population (lower edge of last bin) "<<endl;
  }
  else
  {
    report1<<nlength<<" # number of population size bins "<<endl;
    report1<<len_bins<<endl;
  }

  report1<<endl<<min_tail<<" #_comp_tail_compression"<<endl;
  report1<<min_comp<<" #_add_to_comp"<<endl;
  report1<<CombGender_l<<" #_combine males into females at or below this bin number"<<endl;
  report1<<nlen_bin<<" #_N_LengthBins"<<endl<<len_bins_dat<<endl;
  report1<<sum(Nobs_l)<<" #_N_Length_obs"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Nsamp datavector(female-male)"<<endl;
   for(f=1;f<=Ntypes;f++)
    {
    if(Nobs_l(f)>0)
    {
     for(i=1;i<=Nobs_l(f);i++)
     {
      if(nsamp_l(f,i)>0 && header_l(f,i,3)>0) // do only if this was a real observation
      {
       k=1000;  if(nsamp_l(f,i)<k) k=nsamp_l(f,i);
       exp_l_temp_dat.initialize();
       temp_probs = value(exp_l(f,i));
       temp_mult.fill_multinomial(radm,temp_probs);  // create multinomial draws with prob = expected values
       for (compindex=1; compindex<=k; compindex++) // cumulate the multinomial draws by index in the new data
       {exp_l_temp_dat(temp_mult(compindex)) += 1.0;}
      }
      else
      {exp_l_temp_dat = obs_l(f,i);}
     report1 << header_l(f,i)(1,3)<<" "<<gen_l(f,i)<<" "<<mkt_l(f,i)<<" "<<nsamp_l(f,i)<<" "<<exp_l_temp_dat<<endl;
    }}}

   report1 <<endl<<n_abins<<" #_N_age_bins"<<endl;
  if(n_abins>0) report1<<age_bins1<<endl;
  report1 << N_ageerr <<" #_N_ageerror_definitions"<< endl;
  if(N_ageerr>0) report1 << age_err_rd << endl;

  report1 <<Nobs_a_tot<<" #_N_Agecomp_obs"<<endl;
  report1<<Lbin_method<<" #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths"<<endl;
  report1<<CombGender_a<<" #_combine males into females at or below this bin number"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)"<<endl;
   if(Nobs_a_tot>0)
   for(f=1;f<=Ntypes;f++)
    {
    if(Nobs_a(f)>=1)
    {
     for(i=1;i<=Nobs_a(f);i++)
     {
     if(nsamp_a(f,i)>0 && header_a(f,i,3)>0) // if real observation
     {
      k=1000;  if(nsamp_a(f,i)<k) k=nsamp_a(f,i);  // note that nsamp is adjusted by var_adjust, so var_adjust
                                                   // should be reset to 1.0 in control files that read the nudata.dat files
      exp_a_temp = 0.0;
      temp_probs2 = value(exp_a(f,i));
      temp_mult.fill_multinomial(radm,temp_probs2);
      for (compindex=1; compindex<=k; compindex++) // cumulate the multinomial draws by index in the new data
      {exp_a_temp(temp_mult(compindex)) += 1.0;}
     }
     else
     {exp_a_temp = obs_a(f,i);}
    report1 << header_a(f,i)(1,8)<<" "<<nsamp_a(f,i)<<" "<<exp_a_temp<<endl;
    }
    }
    }
  report1 <<endl<<nobs_ms_tot<<" #_N_MeanSize-at-Age_obs"<<endl;
  report1<<"#Yr Seas Flt/Svy Gender Part Ageerr Ignore datavector(female-male)"<<endl;
  report1<<"#                                          samplesize(female-male)"<<endl;
   if(nobs_ms_tot>0)
   for(f=1;f<=Ntypes;f++)
    {
    if(Nobs_ms(f)>0)
    {
     for(i=1;i<=Nobs_ms(f);i++)
     {
     report1 << header_ms(f,i)(1,7);
     for(a=1;a<=n_abins2;a++)
     {
     report1 << " " ;
         if(obs_ms_n(f,i,a)>0)
          {temp=exp_ms(f,i,a)+randn(radm)*exp_ms_sq(f,i,a)/obs_ms_n(f,i,a);
          if(temp<=0.) {temp=0.0001;}
          report1 << temp;}
         else
             {report1 << exp_ms(f,i,a) ;}
         }

     report1 << endl<< elem_prod(obs_ms_n(f,i),obs_ms_n(f,i)) << endl;
      }
      }
      }
    report1 << endl << N_envvar<<" #_N_environ_variables"<<endl<<N_envdata<<" #_N_environ_obs"<<endl;
    if(N_envdata>0) report1<<env_temp<<endl;

  report1<<SzFreq_Nmeth<<" # N sizefreq methods to read "<<endl;
  if(SzFreq_Nmeth>0)
  {
    report1<<SzFreq_Nbins<<" #Sizefreq N bins per method"<<endl;
    report1<<SzFreq_units<<" #Sizetfreq units(bio/num) per method"<<endl;
    report1<<SzFreq_scale<<" #Sizefreq scale(kg/lbs/cm/inches) per method"<<endl;
    report1<<SzFreq_mincomp<<" #Sizefreq mincomp per method "<<endl;
    report1<<SzFreq_nobs<<" #Sizefreq N obs per method"<<endl;
    report1<<"#_Sizefreq bins "<<endl;
    for(i=1;i<=SzFreq_Nmeth;i++) {report1<<SzFreq_Omit_Small(i)*SzFreq_bins1(i,1)<<SzFreq_bins1(i)(2,SzFreq_Nbins(i))<<endl;}
    report1<<"#_Year season Fleet Gender Partition SampleSize <data> "<<endl;
    j=2*max(SzFreq_Nbins);
    dvector temp_probs3(1,j);
    dvector SzFreq_newdat(1,j);
    for(iobs=1;iobs<=SzFreq_totobs;iobs++)
    {
      if(SzFreq_obs_hdr(iobs,3)>0)  // flag for date range in bounds and used
      {
       j=1000;  if(SzFreq_obs1(iobs,7)<j) j=SzFreq_obs1(iobs,7);
       SzFreq_newdat.initialize();
       temp_probs3(1,SzFreq_Setup2(iobs)) = value(SzFreq_exp(iobs));
       temp_mult.fill_multinomial(radm,temp_probs3(1,SzFreq_Setup2(iobs)));  // create multinomial draws with prob = expected values
       for (compindex=1; compindex<=j; compindex++) // cumulate the multinomial draws by index in the new data
       {SzFreq_newdat(temp_mult(compindex)) += 1.0;}

        report1<<SzFreq_obs1(iobs)(1,7)<<" "<<SzFreq_newdat(1,SzFreq_Setup2(iobs))<<endl;
      }
      else
      {
        report1<<SzFreq_obs1(iobs)<<endl;
      }
    }
  }
  // begin tagging data section #3 (bootstrap data)
  report1<<endl<<Do_TG<<" # do tags (0/1)"<<endl;
  if(Do_TG>0)
  {
    dvector temp_negbin(1,1);

    // info on dimensions of tagging data
    report1<<N_TG<<" # N tag groups"<<endl;
    report1<<N_TG_recap<<" # N recap events"<<endl;
    report1<<TG_mixperiod<<" # mixing latency period: N periods to delay before comparing observed to expected recoveries (0 = release period)"<<endl;
    report1<<TG_maxperiods<<" # max periods (seasons) to track recoveries, after which tags enter accumulator"<<endl;

    // tag releases
    report1<<"# Release data for each tag group.  Tags are considered to be released at the beginning of a season (period)"<<endl;
    report1<<"#<TG> area yr season <tfill> gender age Nrelease  (note that the TG and tfill values are placeholders and are replaced by program generated values)"<<endl;
    report1<<TG_release<<endl;

    // tag recaptures
    report1<<"#_Note: Bootstrap values for tag recaptures are produced only for the same combinations of"<<endl; 
    report1<<"#       group, year, area, and fleet that had observed recaptures. "<<endl;
    report1<<"#_TAG  Year Season Fleet Nrecap"<<endl;
    for(j=1;j<=N_TG_recap;j++)
    {
      // fill in first 4 columns:
      for(k=1;k<=4;k++) report1<<TG_recap_data(j,k)<<" ";
      // fill in 5th column with bootstrap values
      temp_negbin.initialize();
      TG=TG_recap_data(j,1);
      overdisp=TG_parm(2*N_TG+TG);
      t=styr+int((TG_recap_data(j,2)-styr)*nseas+TG_recap_data(j,3)-1) - TG_release(TG,5); // find elapsed time in terms of number of seasons
      if(t>TG_maxperiods) t=TG_maxperiods;
      // some robustification of expected recaps might be needed
      // for cases where the TG_recap_exp = 0
      temp_negbin.fill_randnegbinomial(value(TG_recap_exp(TG,t,0)), value(overdisp), radm);
      report1<<temp_negbin<<" #_orig_obs: "<<TG_recap_data(j,5);
      report1<<" #_exp: "<<value(TG_recap_exp(TG,t,0))<<" #_overdisp: "<<value(overdisp)<<endl;
    }
  }
  // end tagging data section #3 (bootstrap data)

    report1<<endl<<0<<" # no morphcomp data "<<endl;
    report1 << endl << "999" << endl << endl;
  }
  }

  report1 << "ENDDATA" << endl;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 39 write_nucontrol  write new control file */
FUNCTION void write_nucontrol()
  {
  cout<<" Write new starter file "<<endl;
  ofstream NuStart("starter.ss_new");
  NuStart<<version_info_short<<endl;
  if(N_SC>0) NuStart<<Starter_Comments<<endl;
  NuStart<<datfilename<<endl<<ctlfilename<<endl;
  NuStart<<readparfile<<" # 0=use init values in control file; 1=use ss3.par"<<endl;
  NuStart<<rundetail<<" # run display detail (0,1,2)"<<endl;
  NuStart<<reportdetail<<" # detailed age-structured reports in REPORT.SSO (0,1) "<<endl;
  NuStart<<docheckup<<" # write detailed info from first call to echoinput.sso (0,1) "<<endl;
  NuStart<<Do_ParmTrace<<" # write parm values to ParmTrace.sso (0=no,1=good,active; 2=good,all; 3=every_iter,all_parms; 4=every,active)"<<endl;
  NuStart<<Do_CumReport<<" # write to cumreport.sso (0=no,1=like&timeseries; 2=add survey fits)"<<endl;
  NuStart<<Do_all_priors<<" # Include prior_like for non-estimated parameters (0,1) "<<endl;
  NuStart<<SoftBound<<" # Use Soft Boundaries to aid convergence (0,1) (recommended)"<<endl;
  NuStart<<N_nudata<<" # Number of datafiles to produce: 1st is input, 2nd is estimates, 3rd and higher are bootstrap"<<endl;
  NuStart<<Turn_off_phase<<" # Turn off estimation for parameters entering after this phase"<<endl;
  NuStart<<burn_intvl<<" # MCeval burn interval"<<endl;
  NuStart<<thin_intvl<<" # MCeval thin interval"<<endl;
  NuStart<<jitter<<" # jitter initial parm value by this fraction"<<endl;
  NuStart<<STD_Yr_min<<" # min yr for sdreport outputs (-1 for styr)"<<endl;
  NuStart<<STD_Yr_max<<" # max yr for sdreport outputs (-1 for endyr; -2 for endyr+Nforecastyrs"<<endl;
  NuStart<<N_STD_Yr_RD<<" # N individual STD years "<<endl;
  NuStart<<"#vector of year values "<<endl<<STD_Yr_RD<<endl;

  NuStart<<final_conv<<" # final convergence criteria (e.g. 1.0e-04) "<<endl;
  NuStart<<retro_yr-endyr<<" # retrospective year relative to end year (e.g. -4)"<<endl;
  NuStart<<Smry_Age<<" # min age for calc of summary biomass"<<endl;
  NuStart<<depletion_basis<<" # Depletion basis:  denom is: 0=skip; 1=rel X*B0; 2=rel X*Bmsy; 3=rel X*B_styr"<<endl;
  NuStart<<depletion_level<<" # Fraction (X) for Depletion denominator (e.g. 0.4)"<<endl;
  NuStart<<SPR_reporting<<" # SPR_report_basis:  0=skip; 1=(1-SPR)/(1-SPR_tgt); 2=(1-SPR)/(1-SPR_MSY); 3=(1-SPR)/(1-SPR_Btarget); 4=rawSPR"<<endl;
  NuStart<<F_reporting<<" # F_report_units: 0=skip; 1=exploitation(Bio); 2=exploitation(Num); 3=sum(Frates); 4=true F for range of ages"<<endl;
  if(F_reporting==4)
  {NuStart<<F_reporting_ages<<" #_min and max age over which average F will be calculated"<<endl;}
  else
  {NuStart<<"#COND 10 15 #_min and max age over which average F will be calculated with F_reporting=4"<<endl;}
  NuStart<<F_report_basis<<" # F_report_basis: 0=raw; 1=F/Fspr; 2=F/Fmsy ; 3=F/Fbtgt"<<endl;
  NuStart<<999<<" # check value for end of file"<<endl;

  cout<<" Write new forecast file "<<endl;
  ofstream NuFore("forecast.ss_new");
  NuFore<<version_info_short<<endl;
  if(N_FC>0) NuFore<<Forecast_Comments<<endl;
  NuFore<<"# for all year entries except rebuilder; enter either: actual year, -999 for styr, 0 for endyr, neg number for rel. endyr"<<endl;
  NuFore<<Do_Benchmark<<" # Benchmarks: 0=skip; 1=calc F_spr,F_btgt,F_msy "<<endl;
  NuFore<<Do_MSY<<" # MSY: 1= set to F(SPR); 2=calc F(MSY); 3=set to F(Btgt); 4=set to F(endyr) "<<endl;
  NuFore<<SPR_target<<" # SPR target (e.g. 0.40)"<<endl;
  NuFore<<BTGT_target<<" # Biomass target (e.g. 0.40)"<<endl;
  NuFore<<"#_Bmark_years: beg_bio, end_bio, beg_selex, end_selex, beg_relF, end_relF (enter actual year, or values of 0 or -integer to be rel. endyr)"<<endl<<Bmark_Yr_rd<<endl;
  NuFore<<"# "<<Bmark_Yr<<" # after processing "<<endl;
  NuFore<<Bmark_RelF_Basis<<" #Bmark_relF_Basis: 1 = use year range; 2 = set relF same as forecast below"<<endl;
  NuFore<<"#"<<endl<<Do_Forecast<<" # Forecast: 0=none; 1=F(SPR); 2=F(MSY) 3=F(Btgt); 4=Ave F (uses first-last relF yrs); 5=input annual F scalar"<<endl;
  NuFore<<N_Fcast_Yrs<<" # N forecast years "<<endl;
  NuFore<<Fcast_Flevel<<" # F scalar (only used for Do_Forecast==5)"<<endl;
  NuFore<<"#_Fcast_years:  beg_selex, end_selex, beg_relF, end_relF  (enter actual year, or values of 0 or -integer to be rel. endyr)"<<endl<<Fcast_Input(3,6)<<endl;
  NuFore<<"# "<<Fcast_yr<<" # after processing "<<endl;
  NuFore<<HarvestPolicy<<" # Control rule method (1=catch=f(SSB) west coast; 2=F=f(SSB) ) "<<endl;
  NuFore<<H4010_top<<" # Control rule Biomass level for constant F (as frac of Bzero, e.g. 0.40); (Must be > the no F level below) "<<endl;
  NuFore<<H4010_bot<<" # Control rule Biomass level for no F (as frac of Bzero, e.g. 0.10) "<<endl;
  NuFore<<H4010_scale<<" # Control rule target as fraction of Flimit (e.g. 0.75) "<<endl;

  NuFore<<Fcast_Loop_Control(1)<<" #_N forecast loops (1=OFL only; 2=ABC; 3=get F from forecast ABC catch with allocations applied)"<<endl;
  NuFore<<Fcast_Loop_Control(2)<<" #_First forecast loop with stochastic recruitment"<<endl;
  NuFore<<Fcast_Loop_Control(3)<<" #_Forecast loop control #3 (reserved for future bells&whistles) "<<endl;
  NuFore<<Fcast_Loop_Control(4)<<" #_Forecast loop control #4 (reserved for future bells&whistles) "<<endl;
  NuFore<<Fcast_Loop_Control(5)<<" #_Forecast loop control #5 (reserved for future bells&whistles) "<<endl;
  NuFore<<Fcast_Cap_FirstYear<<"  #FirstYear for caps and allocations (should be after years with fixed inputs) "<<endl;

  NuFore<<Impl_Error_Std<<" # stddev of log(realized catch/target catch) in forecast (set value>0.0 to cause active impl_error)"<<endl;

  NuFore<<Do_Rebuilder<<" # Do West Coast gfish rebuilder output (0/1) "<<endl;
  NuFore<<Rebuild_Ydecl<<" # Rebuilder:  first year catch could have been set to zero (Ydecl)(-1 to set to 1999)"<<endl;
  NuFore<<Rebuild_Yinit<<" # Rebuilder:  year for current age structure (Yinit) (-1 to set to endyear+1)"<<endl;

  NuFore<<Fcast_RelF_Basis<<" # fleet relative F:  1=use first-last alloc year; 2=read seas(row) x fleet(col) below"<<endl;
  NuFore<<"# Note that fleet allocation is used directly as average F if Do_Forecast=4 "<<endl;

  NuFore<<Fcast_Catch_Basis<<" # basis for fcast catch tuning and for fcast catch caps and allocation  (2=deadbio; 3=retainbio; 5=deadnum; 6=retainnum)"<<endl;

    NuFore<<"# Conditional input if relative F choice = 2"<<endl;
    NuFore<<"# Fleet relative F:  rows are seasons, columns are fleets"<<endl;
    if(Fcast_RelF_Basis==1)
    {
      NuFore<<"#_Fleet: ";
      for(f=1;f<=Nfleet;f++) NuFore<<" "<<fleetname(f);
      NuFore<<endl;
      for(s=1;s<=nseas;s++)
      {
        NuFore<<"# "<<Fcast_RelF_Use(s)<<endl;
      }
    }
    else
    {
      NuFore<<"#_Fleet: ";
      for(f=1;f<=Nfleet;f++) NuFore<<" "<<fleetname(f);
      NuFore<<endl;
      for(s=1;s<=nseas;s++)
      {
        NuFore<<Fcast_RelF_Use(s)<<endl;
      }
    }

  NuFore<<"# max totalcatch by fleet (-1 to have no max) must enter value for each fleet"<<endl;
  NuFore<<Fcast_MaxFleetCatch<<endl;
  NuFore<<"# max totalcatch by area (-1 to have no max); must enter value for each fleet "<<endl;
  NuFore<<Max_Fcast_Catch<<endl;
  NuFore<<"# fleet assignment to allocation group (enter group ID# for each fleet, 0 for not included in an alloc group)"<<endl;
  NuFore<<Allocation_Fleet_Assignments<<endl;
  NuFore<<"#_Conditional on >1 allocation group"<<endl<<"# allocation fraction for each of: "<<Fcast_Catch_Allocation_Groups<<" allocation groups"<<endl;
  if(Fcast_Catch_Allocation_Groups>0) {NuFore<<Fcast_Catch_Allocation<<endl;} else {NuFore<<"# no allocation groups"<<endl;}

  NuFore<<N_Fcast_Input_Catches<<" # Number of forecast catch levels to input (else calc catch from forecast F) "<<endl;
  NuFore<<-1<<" # code means to read fleet/time specific basis (2=dead catch; 3=retained catch; 99=F)  as below (units are from fleetunits; note new codes in SSV3.20)"<<endl;

  NuFore<<"# Input fixed catch values"<<endl;
  NuFore<<"#Year Seas Fleet Catch(or_F) Basis"<<endl;
  if(N_Fcast_Input_Catches>0)
  {
    for(j=1;j<=N_Fcast_Input_Catches;j++)
    {
      y=Fcast_InputCatch_rd(j,1); s=Fcast_InputCatch_rd(j,2); f=Fcast_InputCatch_rd(j,3);
      t=styr+(y-styr)*nseas +s-1;
      NuFore<<Fcast_InputCatch_rd(j)<<" "<<Fcast_InputCatch(t,f,2)<<endl;
    }
  }
  NuFore<<"#"<<endl<<999<<" # verify end of input "<<endl;

//**********************************************************
  cout<<" Write new control file "<<endl;
  ofstream report4("control.ss_new");
  report4<<version_info_short<<endl;
  if(N_CC>0) report4<<Control_Comments<<endl;
  report4 << "#_data_and_control_files: "<<datfilename<<" // "<<ctlfilename<<endl;
  report4<<"#_"<<version_info<<endl;
  report4 << N_GP << "  #_N_Growth_Patterns"<<endl;
  report4 << N_platoon << " #_N_Morphs_Within_GrowthPattern "<<endl;
  if(N_platoon==1) report4<<"#_Cond ";
  report4<<sd_ratio<<" #_Morph_between/within_stdev_ratio (no read if N_morphs=1)"<<endl;
  if(N_platoon==1) report4<<"#_Cond ";
  report4<<submorphdist(1,N_platoon)<<" #vector_Morphdist_(-1_in_first_val_gives_normal_approx)"<<endl;
  report4<<"#"<<endl;
  if(N_GP*nseas*pop==1)
  {
    report4<<"#_Cond 0  #  N recruitment designs goes here if N_GP*nseas*area>1"<<endl<<
                   "#_Cond 0  #  placeholder for recruitment interaction request"<<endl<<
                   "#_Cond 1 1 1  # example recruitment design element for GP=1, seas=1, area=1"<<endl;
  }
  else
  {
    report4<<recr_dist_read<<" #  number of recruitment assignments (overrides GP*area*seas parameter values) "<<endl<<
             recr_dist_inx<< " # recruitment interaction requested"<<endl<<
             "#GP seas area for each recruitment assignment"<<endl<<recr_dist_pattern_1<<endl;
  }
  report4<<"#"<<endl;
  if(pop==1)
  {report4<<"#_Cond 0 # N_movement_definitions goes here if N_areas > 1"<<endl
    <<"#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0"<<endl
    <<"#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10"<<endl;}
  else
  {
    report4<<endl<<do_migration<<" #_N_movement_definitions"<<endl;
    if(do_migration>0)
    {
      report4<<migr_firstage<<" # first age that moves (real age at begin of season, not integer)"<<endl
      <<"# seas,GP,source_area,dest_area,minage,maxage"<<endl<<move_def<<endl;
    }
    else
    {
    report4<<"#_Cond 1.0 # first age that moves (real age at begin of season, not integer) if do_migration>0"<<endl
    <<"#_Cond 1 1 1 2 4 10 # example move definition for seas=1, GP=1, source=1 dest=2, age1=4, age2=10"<<endl;
    }
  }
  report4<<"#"<<endl;
  report4<<N_Block_Designs<<" #_Nblock_Patterns"<<endl;
  if(N_Block_Designs>0)
  {report4<<Nblk<<" #_blocks_per_pattern "<<endl<<"# begin and end years of blocks"<<endl<<Block_Design<<endl;}
  else
  {report4<<"#_Cond "<<0<<" #_blocks_per_pattern "<<endl<<"# begin and end years of blocks"<<endl;}
  report4<<"#"<<endl;

  report4<<fracfemale<<" #_fracfemale "<<endl;
  report4<<natM_type<<" #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate"<<endl;
    if(natM_type==1)
    {report4<<N_natMparms<<" #_N_breakpoints"<<endl<<NatM_break<<" # age(real) at M breakpoints"<<endl;}
    else if(natM_type==2)
    {report4<<natM_amin<<" #_reference age for Lorenzen M; read 1P per morph"<<endl;}
    else if(natM_type>=3)
    {report4<<" #_Age_natmort_by gender x growthpattern"<<endl<<Age_NatMort<<endl;}
    else
    {report4<<"  #_no additional input for selected M option; read 1P per morph"<<endl;}

    report4<<Grow_type<<" # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_speciific_K; 4=not implemented"<<endl;
    if(Grow_type<=3)
    {report4<<AFIX<<" #_Growth_Age_for_L1"<<endl<<AFIX2<<" #_Growth_Age_for_L2 (999 to use as Linf)"<<endl;}
    else
    {report4<<" #_No Growth_Age_for_L1 and L2"<<endl;}
    if(Grow_type==3)
      {report4<<Age_K_count<<" # number of K multipliers to read"<<endl<<Age_K_points<<" # ages for K multiplier"<<endl;}

    report4<<SD_add_to_LAA<<" #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)"<<endl;   // constant added to SD length-at-age (set to 0.1 for compatibility with SS2 V1.x
    report4<<CV_depvar<<" #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)"<<endl;
    report4<<Maturity_Option<<" #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity by GP; 4=read age-fecundity by GP; 5=read fec and wt from wtatage.ss; 6=read length-maturity by GP"<<endl;
    if(Maturity_Option==3)
    {report4<<"#_Age_Maturity by growth pattern for females"<<endl<<Age_Maturity<<endl;}
    else if(Maturity_Option==4)
    {report4<<"#_Age_Fecundity by growth pattern for females"<<endl<<Age_Maturity<<endl;}
    else if(Maturity_Option==6)
    {report4<<"#_Length_Maturity by growth pattern for females"<<endl<<Length_Maturity<<endl;}
    else
    {report4<<"#_placeholder for empirical age- or length- maturity by growth pattern (female only)"<<endl;}
    report4<<First_Mature_Age<<" #_First_Mature_Age"<<endl;

    report4<<Fecund_Option<<" #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W"<<endl;
    report4<<Hermaphro_Option<<" #_hermaphroditism option:  0=none; 1=age-specific fxn"<<endl;
    if (Hermaphro_Option>0) report4<<Hermaphro_seas<<" # Hermaphro_season "<<endl<<Hermaphro_maleSPB<<" # Hermaphro_maleSPB "<<endl;
    report4<<MGparm_def<<" #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)"<<endl;
    report4<<MG_adjust_method<<
    " #_env/block/dev_adjust_method (1=standard; 2=logistic transform keeps in base parm bounds; 3=standard w/ no bound check)"<<endl;
  report4<<"#"<<endl;
  report4<<"#_growth_parms"<<endl;
  report4<<"#_LO HI INIT PRIOR PR_type SD PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn"<<endl;
  NP=0;
  for (f=1;f<=N_MGparm;f++)
  {
    NP++;
    MGparm_1(f,3)=value(MGparm(f));
    report4<<MGparm_1(f)<<" # "<<ParmLabel(NP)<<endl;
  }
  report4<<"#"<<endl;
  j=N_MGparm;
  if(N_MGparm_env>0)
  {
    report4<<1<<" #_custom_MG-env_setup (0/1)"<<endl;
    for (f=1;f<=N_MGparm_env;f++)
    {j++; NP++; if(customMGenvsetup==0) {k=1;} else {k=f;}  // use read value of custom here
    MGparm_env_1(k,3)=value(MGparm(j)); report4<<MGparm_env_1(k)<<" # "<<ParmLabel(NP)<<endl;}
   }
  else
  {
    report4<<"#_Cond 0  #custom_MG-env_setup (0/1)"<<endl;
    report4<<"#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no MG-environ parameters"<<endl;
  }
  report4<<"#"<<endl;
  if(N_MGparm_blk>0)
  {
    report4<<1<<" #_custom_MG-block_setup (0/1)"<<endl;
    for (f=1;f<=N_MGparm_blk;f++)
    {j++; NP++; if(customblocksetup_MG==0) {k=1;} else {k=f;}
    MGparm_blk_1(k,3)=value(MGparm(j));report4<<MGparm_blk_1(k)<<" # "<<ParmLabel(NP)<<endl;}
  }
  else
  {
    report4<<"#_Cond 0  #custom_MG-block_setup (0/1)"<<endl;
    report4<<"#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no MG-block parameters"<<endl;
  }
   if(N_MGparm_trend>0)
   {
     report4<<"#_MGtrend_&_cycle_parms "<<endl;
     for (f=1;f<=N_MGparm_trend2;f++)
     {
       j++;  NP++;
       MGparm_trend_1(f,3)=value(MGparm(j)); report4<<MGparm_trend_1(f)<<" # "<<ParmLabel(NP)<<endl;
     }
   }
  else
  {
    report4<<"#_Cond No MG parm trends "<<endl;
  }


  report4<<"#"<<endl;
  report4<<"#_seasonal_effects_on_biology_parms"<<endl<<MGparm_seas_effects<<" #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K"<<endl;
  if(MGparm_doseas>0)
  {
    for(f=1;f<=N_MGparm_seas;f++)
    {
      NP++; j++;  MGparm_seas_1(f,3)=value(MGparm(j));
      report4<<MGparm_seas_1(f)<<" # "<<ParmLabel(NP)<<endl;
    }
  }
  else
  {
    report4<<"#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters"<<endl;
  }
  report4<<"#"<<endl;
  if(N_MGparm_dev>0)
  {
   for(i=1;i<=N_MGparm_dev;i++)
   for(j=MGparm_dev_minyr(i);j<=MGparm_dev_maxyr(i);j++)
   {
     NP++;
    report4<<"#DisplayOnly "<<MGparm_dev(i,j)<<" # "<<ParmLabel(NP)<<endl;
   }
    report4<<MGparm_dev_PH<<" #_MGparm_Dev_Phase"<<endl;
  }
  else
  {
    report4<<"#_Cond -4 #_MGparm_Dev_Phase"<<endl;
  }

  report4<<"#"<<endl;
   report4<<"#_Spawner-Recruitment"<<endl<<SR_fxn<<" #_SR_function: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepard_3Parm"<<endl;
   report4<<"#_LO HI INIT PRIOR PR_type SD PHASE"<<endl;
   for(f=1;f<=N_SRparm2;f++)
   { NP++;
     SR_parm_1(f,3)=value(SR_parm(f));
     report4<<SR_parm_1(f)<<" # "<<ParmLabel(NP)<<endl;
   }
   report4<<SR_env_link<<" #_SR_env_link"<<endl;
   report4<<SR_env_target_RD<<" #_SR_env_target_0=none;1=devs;_2=R0;_3=steepness"<<endl;
   report4<<do_recdev<<" #do_recdev:  0=none; 1=devvector; 2=simple deviations"<<endl;
   report4<<recdev_start<<" # first year of main recr_devs; early devs can preceed this era"<<endl;
   report4<<recdev_end<<" # last year of main recr_devs; forecast devs start in following year"<<endl;
   report4<<recdev_PH<<" #_recdev phase "<<endl;
   report4<<recdev_adv<<" # (0/1) to read 13 advanced options"<<endl;
   if(recdev_adv==0) {onenum="#_Cond ";} else {onenum=" ";}
   report4<<onenum<<recdev_early_start_rd<<" #_recdev_early_start (0=none; neg value makes relative to recdev_start)"<<endl;
   report4<<onenum<<recdev_early_PH<<" #_recdev_early_phase"<<endl;
   report4<<onenum<<Fcast_recr_PH<<" #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)"<<endl;
   report4<<onenum<<Fcast_recr_lambda<<" #_lambda for Fcast_recr_like occurring before endyr+1"<<endl;
   report4<<onenum<<recdev_adj(1)<<" #_last_early_yr_nobias_adj_in_MPD"<<endl;
   report4<<onenum<<recdev_adj(2)<<" #_first_yr_fullbias_adj_in_MPD"<<endl;
   report4<<onenum<<recdev_adj(3)<<" #_last_yr_fullbias_adj_in_MPD"<<endl;
   report4<<onenum<<recdev_adj(4)<<" #_first_recent_yr_nobias_adj_in_MPD"<<endl;
   report4<<onenum<<recdev_adj(5)<<" #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)"<<endl;
   report4<<onenum<<recdev_cycle<<" #_period of cycles in recruitment (N parms read below)"<<endl;
   report4<<onenum<<recdev_LO<<" #min rec_dev"<<endl;
   report4<<onenum<<recdev_HI<<" #max rec_dev"<<endl;
   report4<<onenum<<recdev_read<<" #_read_recdevs"<<endl;
   report4<<"#_end of advanced SR options"<<endl;
   report4<<"#"<<endl;
   if(recdev_cycle>0)
   {
     for(y=1;y<=recdev_cycle;y++)
     {
       NP++;
       recdev_cycle_parm_RD(y,3)=value(recdev_cycle_parm(y));
       report4<<recdev_cycle_parm_RD(y)<<" # "<<ParmLabel(NP)<<endl;
     }
   }
   else
   {
     report4<<"#_placeholder for full parameter lines for recruitment cycles"<<endl;
   }
    if(recdev_read>0)
    {
      report4<<"# Specified recr devs to read"<<endl;
      report4<<"#_Yr Input_value # Final_value"<<endl;
      for(j=1;j<=recdev_read;j++)
      {
        y=recdev_input(j,1);
        report4<<recdev_input(j)<<" # ";
        if(y>=recdev_first)
        {report4<<recdev(y)<<endl;}
        else
        {report4<<" not used "<<endl;}
      }
    }
    else
    {
      report4<<"# read specified recr devs"<<endl;
      report4<<"#_Yr Input_value"<<endl;
    }
    report4<<"#"<<endl;
    report4<<"# all recruitment deviations"<<endl;

  if(recdev_do_early>0)
  {
    for(y=recdev_early_start;y<=recdev_early_end;y++)
    {
      NP++;  report4<<"#DisplayOnly "<<recdev(y)<<" # "<<ParmLabel(NP)<<endl;
    }
  }

    if(do_recdev>0)
    {
      for(y=recdev_start;y<=recdev_end;y++)
      {
        NP++;  report4<<"#DisplayOnly "<<recdev(y)<<" # "<<ParmLabel(NP)<<endl;
      }
    }

    if(Do_Forecast>0)
    {
      for(y=recdev_end+1;y<=YrMax;y++)
      {
        NP++;  report4<<"#DisplayOnly "<<recdev(y)<<" # "<<ParmLabel(NP)<<endl;
      }
      for(y=endyr+1;y<=YrMax;y++)
      {
        NP++;  report4<<"#DisplayOnly "<<Fcast_impl_error(y)<<" # "<<ParmLabel(NP)<<endl;
      }
    }
  report4<<"#"<<endl;
  report4<<"#Fishing Mortality info "<<endl<<F_ballpark<<" # F ballpark for annual F (=Z-M) for specified year"<<endl;
  report4<<F_ballpark_yr<<" # F ballpark year (neg value to disable)"<<endl;
  report4<<F_Method<<" # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)"<<endl;
  report4<<max_harvest_rate<<" # max F or harvest rate, depends on F_Method"<<endl;
  report4<<"# no additional F input needed for Fmethod 1"<<endl;
  report4<<"# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read"<<endl;
  report4<<"# if Fmethod=3; read N iterations for tuning for Fmethod 3"<<endl;
 if(F_Method==2)
  {
    report4<<F_setup<<" # overall start F value; overall phase; N detailed inputs to read"<<endl;
    report4<<"#Fleet Year Seas F_value se phase (for detailed setup of F_Method=2)"<<endl<<F_setup2<<endl;
  }
  else if(F_Method==3)
  {report4<<F_Tune<<"  # N iterations for tuning F in hybrid method (recommend 3 to 7)"<<endl;}

   report4<<"#"<<endl;
   report4<<"#_initial_F_parms"<<endl;
   report4<<"#_LO HI INIT PRIOR PR_type SD PHASE"<<endl;
   for(f=1;f<=Nfleet;f++)
   {
    NP++;
    init_F_parm_1(f,3)=value(init_F(f));
    report4<<init_F_parm_1(f)<<" # "<<ParmLabel(NP)<<endl;
   }

   if(F_Method==2)
   {
    report4<<endl<<"# F rates for Fmethod=2"<<endl;
    k=0;
    for(f=1;f<=Nfleet;f++)
    for(y=styr;y<=endyr;y++)
    for(s=1;s<=nseas;s++)
    {
      NP++; k++;
      report4<<"# "<<F_rate(k)<< " "<<ParmLabel(NP)<<endl;
    }
   }
   report4<<"#"<<endl;
   report4<<"#_Q_setup"<<endl<<
   " # Q_type options:  <0=mirror, 0=float_nobiasadj, 1=float_biasadj, 2=parm_nobiasadj, 3=parm_w_random_dev, 4=parm_w_randwalk, 5=mean_unbiased_float_assign_to_parm"<<endl;
   report4<<"#_for_env-var:_enter_index_of_the_env-var_to_be_linked"<<endl;
   report4<<"#_Den-dep  env-var  extra_se  Q_type"<<endl;
   for(f=1;f<=Ntypes;f++)
   report4<<Q_setup(f)<<" # "<<f<<" "<<fleetname(f)<<endl;
   report4<<"#"<<endl;
   if(ask_detail>0)  // report q_parm_detail
   {
    report4<<1<<" #_0=read one parm for each fleet with random q; 1=read a parm for each year of index"<<endl;
   }
  else
  {
    report4<<"#_Cond 0 #_If q has random component, then 0=read one parm for each fleet with random q; 1=read a parm for each year of index"<<endl;
  }

   report4<<"#_Q_parms(if_any);Qunits_are_ln(q)"<<endl;
   if(Q_Npar>0)
   {
    report4<<"# LO HI INIT PRIOR PR_type SD PHASE"<<endl;
    for(f=1;f<=Q_Npar;f++)
    {
      NP++;
      Q_parm_1(f,3)=value(Q_parm(f));
      report4<<Q_parm_1(f)<<" # "<<ParmLabel(NP)<<endl;
    }
   }
   report4<<"#"<<endl;
   report4<<"#_size_selex_types"<<endl;
   report4<<"#discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead"<<endl;
   report4<<"#_Pattern Discard Male Special"<<endl;
   for(f=1;f<=Ntypes;f++) report4<<seltype(f)<<" # "<<f<<" "<<fleetname(f)<<endl;
   report4<<"#"<<endl;
   report4<<"#_age_selex_types"<<endl;
   report4<<"#_Pattern ___ Male Special"<<endl;
   for(f=1;f<=Ntypes;f++) report4<<seltype(f+Ntypes)<<" # "<<f<<" "<<fleetname(f)<<endl;
   report4<<"#_LO HI INIT PRIOR PR_type SD PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn"<<endl;

  for (f=1;f<=N_selparm;f++)
  {
    NP++;
    selparm_1(f,3)=value(selparm(f));
    report4<<selparm_1(f)<<" # "<<ParmLabel(NP)<<endl;
  }

  j=N_selparm;

   if(N_selparm_env>0)
   {
     report4<<1<<" #_custom_sel-env_setup (0/1) "<<endl;
     for (f=1;f<=N_selparm_env;f++)
     {
       j++;  NP++;
       if(customenvsetup==0) {k=1;} else {k=f;}  // use read value of custom here
       selparm_env_1(k,3)=value(selparm(j)); report4<<selparm_env_1(k)<<" # "<<ParmLabel(NP)<<endl;
     }
   }
  else
  {
    report4<<"#_Cond 0 #_custom_sel-env_setup (0/1) "<<endl;
    report4<<"#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no enviro fxns"<<endl;
  }

   if(N_selparm_blk>0)
   {
     report4<<1<<" #_custom_sel-blk_setup (0/1) "<<endl;
     for (f=1;f<=N_selparm_blk;f++)
     {
       j++;  NP++;
       if(customblocksetup==0) {k=1;} else {k=f;}  // use read value of custom here
       selparm_blk_1(k,3)=value(selparm(j)); report4<<selparm_blk_1(k)<<" # "<<ParmLabel(NP)<<endl;
     }
   }
  else
  {
    report4<<"#_Cond 0 #_custom_sel-blk_setup (0/1) "<<endl;
    report4<<"#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no block usage"<<endl;
  }

   if(N_selparm_trend>0)
   {
     report4<<"#_seltrend_parms "<<endl;
     for (f=1;f<=N_selparm_trend2;f++)
     {
       j++;  NP++;
       selparm_trend_1(f,3)=value(selparm(j)); report4<<selparm_trend_1(f)<<" # "<<ParmLabel(NP)<<endl;
     }
   }
  else
  {
    report4<<"#_Cond No selex parm trends "<<endl;
  }


  if(N_selparm_dev>0)
  {
    for(i=1;i<=N_selparm_dev;i++)
    for(j=selparm_dev_minyr(i);j<=selparm_dev_maxyr(i);j++)
    {
      NP++; report4<<"# "<<selparm_dev(i,j)<<" # "<<ParmLabel(NP)<<endl;
    }
    report4<<selparm_dev_PH<<" #_selparmdev-phase"<<endl;
  }
  else
  {
    report4<<"#_Cond -4 # placeholder for selparm_Dev_Phase"<<endl;
  }

  if(N_selparm_env+N_selparm_blk+N_selparm_dev == 0) report4<<"#_Cond ";
   report4<<selparm_adjust_method<<
   " #_env/block/dev_adjust_method (1=standard; 2=logistic trans to keep in base parm bounds; 3=standard w/ no bound check)"<<endl;

  report4<<"#"<<endl<<"# Tag loss and Tag reporting parameters go next"<<endl;
  if(Do_TG>0)
  {
    report4<<1<<" # TG_custom:  0=no read; 1=read"<<endl;
    for (f=1;f<=3*N_TG+2*Nfleet;f++)
    {
      NP++;
      report4<<TG_parm2(f)<<" # "<<ParmLabel(NP)<<endl;
    }
  }
  else
  {
    report4<<"0  # TG_custom:  0=no read; 1=read if tags exist"<<endl
    <<"#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters"<<endl;;
  }

  report4<<"#"<<endl<<Do_Var_adjust<<" #_Variance_adjustments_to_input_values"<<endl;
  report4 <<"#_fleet: ";  for(f=1;f<=Ntypes;f++) {report4<<f<<" ";}
  report4 <<endl;
  if(Do_Var_adjust==0)
  {onenum="#_Cond ";}
  else
  {onenum=" ";}

  report4<<onenum<<var_adjust(1)<<" #_add_to_survey_CV"<<endl;
  report4<<onenum<<var_adjust(2)<<" #_add_to_discard_stddev"<<endl;
  report4<<onenum<<var_adjust(3)<<" #_add_to_bodywt_CV"<<endl;
  report4<<onenum<<var_adjust(4)<<" #_mult_by_lencomp_N"<<endl;
  report4<<onenum<<var_adjust(5)<<" #_mult_by_agecomp_N"<<endl;
  report4<<onenum<<var_adjust(6)<<" #_mult_by_size-at-age_N"<<endl;

  report4<<"#"<<endl<<max_lambda_phase<<" #_maxlambdaphase"<<endl;
  report4<<sd_offset<<" #_sd_offset"<<endl;

  report4<<"#"<<endl<<N_lambda_changes<<" # number of changes to make to default Lambdas (default value is 1.0)"<<endl;
  report4<<"# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; "<<
   endl<<"# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark"<<
   endl<<"#like_comp fleet/survey  phase  value  sizefreq_method"<<endl;

  if(N_lambda_changes>0) report4<<Lambda_changes<<endl;

  report4<<"#"<<endl<<"# lambdas (for info only; columns are phases)"<<endl;
  if(nobs_cr>0) {for(f=1;f<=Ntypes;f++) report4<<"# "<<surv_lambda(f)<<" #_CPUE/survey:_"<<f<<endl;}
  if(nobs_disc>0) {for(f=1;f<=Ntypes;f++) report4<<"# "<< disc_lambda(f)<<" #_discard:_"<<f<<endl;}
  if(nobs_mnwt>0) {for(f=1;f<=Ntypes;f++) report4<<"# "<< mnwt_lambda(f)<<" #_meanbodywt:"<<f<<endl;}
  if(Nobs_l_tot>0) {for(f=1;f<=Ntypes;f++) report4<<"# "<< length_lambda(f)<<" #_lencomp:_"<<f<<endl;}
  if(Nobs_a_tot>0) {for(f=1;f<=Ntypes;f++) report4<<"# "<< age_lambda(f)<<" #_agecomp:_"<<f<<endl;}
  if(SzFreq_Nmeth>0) for(f=1;f<=SzFreq_N_Like;f++) report4<<"# "<<SzFreq_lambda(f)<<" #_sizefreq:_"<<f<<endl;
  if(nobs_ms_tot>0) {for(f=1;f<=Ntypes;f++) report4<<"# "<< sizeage_lambda(f)<<" #_size-age:_"<<f<<endl;}
  report4<<"# "<< init_equ_lambda<<" #_init_equ_catch"<<endl;
  report4<<"# "<< recrdev_lambda<<" #_recruitments"<<endl;
  report4<<"# "<< parm_prior_lambda<<" #_parameter-priors"<<endl;
  report4<<"# "<< parm_dev_lambda<<" #_parameter-dev-vectors"<<endl;
  report4<<"# "<< CrashPen_lambda<<" #_crashPenLambda"<<endl;
  if(Do_Morphcomp>0) report4<<"# "<< Morphcomp_lambda<<" #_Morphcomplambda"<<endl;
  if(Do_TG>0)
  {
  for(TG=1;TG<=N_TG;TG++) report4<<"# "<<TG_lambda1(TG)<<" #_TG-comp_group:_"<<TG<<endl;
  for(TG=1;TG<=N_TG;TG++) report4<<"# "<<TG_lambda2(TG)<<" #_TG-negbin_group:_"<<TG<<endl;
  }
  report4<<"# "<<F_ballpark_lambda<<" # F_ballpark_lambda"<<endl;
  report4<<Do_More_Std<<" # (0/1) read specs for more stddev reporting "<<endl;
  if(Do_More_Std>0)
  {
    report4<<More_Std_Input<<" # selex type, len/age, year, N selex bins, Growth pattern, N growth ages, NatAge_area(-1 for all), NatAge_yr, N Natages"<<endl;
    if(More_Std_Input(4)>0) report4<<Selex_Std_Pick<<" # vector with selex std bin picks (-1 in first bin to self-generate)"<<endl;
    if(More_Std_Input(6)>0) report4<<Growth_Std_Pick<<" # vector with growth std bin picks (-1 in first bin to self-generate)"<<endl;
    if(More_Std_Input(9)>0) report4<<NatAge_Std_Pick<<" # vector with NatAge std bin picks (-1 in first bin to self-generate)"<<endl;
  }
  else
  {
    report4<<" # 0 1 -1 5 1 5 1 -1 5 # placeholder for selex type, len/age, year, N selex bins, Growth pattern, N growth ages, NatAge_area(-1 for all), NatAge_yr, N Natages"<<endl;
    report4<<" # placeholder for vector of selex bins to be reported"<<endl;
    report4<<" # placeholder for vector of growth ages to be reported"<<endl;
    report4<<" # placeholder for vector of NatAges ages to be reported"<<endl;
  }
  report4<<fim<<endl<<endl; // end of file indicator
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 40 write_bigoutput */
FUNCTION void write_bigoutput()
  {
  SS2out.open ("Report.sso");   // this file was created in globals so accessible to the report_parm function
  ofstream SS_compout("CompReport.sso");
//  ofstream SS_compout2("CompReport2.sso");
  ofstream SIS_table("SIS_table.sso");

  k=current_phase();
  if(k>max_lambda_phase) k=max_lambda_phase;
  SS2out<<version_info_short<<endl;
  SS2out<<version_info<<endl<<endl;
  time(&finish);
  SS_compout<<version_info_short<<endl;
  SS_compout<<version_info<<endl<<"StartTime: "<<ctime(&start);
  SIS_table<<version_info_short<<endl;
  SIS_table<<version_info<<endl<<"StartTime: "<<ctime(&start);
  SIS_table<<endl<<"Data_File: "<<datfilename<<endl;
  SIS_table<<"Control_File: "<<ctlfilename<<endl;

  cout<<" writing big output now "<<endl;
  SS2out<<"StartTime: "<<ctime(&start);
  SS2out<<"EndTime: "<<ctime(&finish);
  elapsed_time = difftime(finish,start);
  hour = long(elapsed_time)/3600;
  minute = long(elapsed_time)%3600/60;
  second = (long(elapsed_time)%3600)%60;
  SS2out<<"This run took: "<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds."<<endl;
  SS2out<<"Data_File: "<<datfilename<<endl;
  SS2out<<"Control_File: "<<ctlfilename<<endl;
  if(readparfile>=1) SS2out<<"Start_parm_values_from_SS3.PAR"<<endl;
  SS2out<<endl<<"Convergence_Level: "<<objective_function_value::pobjfun->gmax<<" is_final_gradient"<<endl<<"Hessian: "<<endl;
  if(N_SC>0)
  {
    SS2out<<"#"<<endl<<"Starter_Comments"<<endl<<Starter_Comments<<endl;
  }
  if(N_DC>0)
  {
    SS2out<<"#"<<endl<<"Data_Comments"<<endl<<Data_Comments<<endl;
  }
  if(N_CC>0)
  {
    SS2out<<"#"<<endl<<"Control_Comments"<<endl<<Control_Comments<<endl;
  }
  if(N_FC>0)
  {
    SS2out<<"#"<<endl<<"Forecast_Comments"<<endl<<Forecast_Comments<<endl;
  }

  SS2out<<endl<<"KeyWords"<<endl;
  SS2out<<"X DEFINITIONS"<<endl;
  SS2out<<"X LIKELIHOOD"<<endl;
  SS2out<<"X Input_Variance_Adjustment"<<endl;
  SS2out<<"X PARAMETERS"<<endl;
  SS2out<<"X DERIVED_QUANTITIES"<<endl;
//  SS2out<<"X MGParm_Block_Assignments"<<endl;
//  SS2out<<"X Selex_Block_Assignments"<<endl;
  SS2out<<"X MGparm_By_Year_after_adjustments"<<endl;
  SS2out<<"X selparm(Size)_By_Year_after_adjustments"<<endl;
  SS2out<<"X selparm(Age)_By_Year_after_adjustments"<<endl;
  SS2out<<"X RECRUITMENT_DIST"<<endl;
  SS2out<<"X MORPH_INDEXING"<<endl;
  SS2out<<"X SIZEFREQ_TRANSLATION"<<endl;
  SS2out<<"X MOVEMENT"<<endl;
  SS2out<<"X EXPLOITATION"<<endl;
  SS2out<<"X TIME_SERIES"<<endl;
  SS2out<<"X SPR_series"<<endl;
  SS2out<<"X Kobe_Plot"<<endl;
  SS2out<<"X SPAWN_RECRUIT"<<endl;
  SS2out<<"X Spawning_Biomass_Report_1"<<endl;
  SS2out<<"X NUMBERS_AT_AGE_Annual_1"<<endl;
  SS2out<<"X Z_AT_AGE_Annual_1"<<endl;
  SS2out<<"X Spawning_Biomass_Report_2"<<endl;
  SS2out<<"X NUMBERS_AT_AGE_Annual_2"<<endl;
  SS2out<<"X Z_AT_AGE_Annual_2"<<endl;
  SS2out<<"X INDEX_1"<<endl;
  SS2out<<"X INDEX_2"<<endl;
  SS2out<<"X INDEX_3"<<endl;
  SS2out<<"X CATCH"<<endl;
  SS2out<<"X DISCARD_SPECIFICATION"<<endl;
  SS2out<<"X DISCARD_OUTPUT"<<endl;
  SS2out<<"X DISCARD_MORT"<<endl;
  SS2out<<"X MEAN_BODY_WT"<<endl;
  SS2out<<"X FIT_LEN_COMPS"<<endl;
  SS2out<<"X FIT_AGE_COMPS"<<endl;
  SS2out<<"X FIT_SIZE_COMPS"<<endl;
  SS2out<<"X OVERALL_COMPS"<<endl;
  SS2out<<"X LEN_SELEX"<<endl;
  SS2out<<"X RETENTION"<<endl;
  SS2out<<"X KEEPERS"<<endl;
  SS2out<<"X DEADFISH"<<endl;
  SS2out<<"X AGE_SELEX"<<endl;
  SS2out<<"X ENVIRONMENTAL_DATA"<<endl;
  SS2out<<"X TAG_Recapture"<<endl;
  SS2out<<"X NUMBERS_AT_AGE"<<endl;
  SS2out<<"X NUMBERS_AT_LENGTH"<<endl;
  SS2out<<"X CATCH_AT_AGE"<<endl;
  SS2out<<"X BIOLOGY"<<endl;
  SS2out<<"X SPR/YPR_PROFILE"<<endl;
  SS2out<<"X Dynamic_Bzero "<<endl;

  SS2out<<endl<<"DEFINITIONS"<<endl;
  SS2out<<"N_seasons "<<nseas<<endl;
  SS2out<<"Sum_of_months_on_read_was:_ "<<sumseas<<" rescaled_to_sum_to: "<<sum(seasdur)<<endl;
  SS2out<<"Season_Durations "<<seasdur<<endl;
  SS2out<<"fleet_ID#: ";
  for(f=1;f<=Ntypes;f++) SS2out<<" "<<f;
  SS2out<<endl<<"fleet_names: ";
  for(f=1;f<=Ntypes;f++) SS2out<<" "<<fleetname(f);
  SS2out<<endl;
  SS2out<<"fleet_area: "<<fleet_area<<endl;
  SS2out<<"Catch_units(1=bio;2=num): "<<catchunits<<endl;
  SS2out<<"Catch_error(-1_for_bycatch_only): "<<catch_se_rd<<endl;
  SS2out<<"Survey_units(0=num;1=bio;2=F): "<<cr_units<<endl;
  SS2out<<"Survey_error(-1=normal;0=lognorm;>0=df_T): "<<cr_errtype<<endl;

  SS2out<<endl<<"LIKELIHOOD "<<obj_fun<<endl;                         //ProgLabel_310
  SS2out<<"Component logL*Lambda Lambda"<<endl;
  SS2out<<"TOTAL "<<obj_fun<<endl;
  if(F_Method>1) SS2out <<"Catch "<<catch_like*column(catch_lambda,k)<<endl;
  SS2out <<"Equil_catch "<<equ_catch_like*init_equ_lambda(k)<<" "<<init_equ_lambda(k)<<endl;
  if(nobs_cr>0) SS2out <<"Survey "<<surv_like*column(surv_lambda,k)<<endl;
  if(nobs_disc>0) SS2out <<"Discard "<<disc_like*column(disc_lambda,k)<<endl;
  if(nobs_mnwt>0) SS2out <<"Mean_body_wt "<<mnwt_like*column(mnwt_lambda,k)<<endl;
  if(Nobs_l_tot>0) SS2out <<"Length_comp "<<length_like*column(length_lambda,k)<<endl;
  if(Nobs_a_tot>0) SS2out <<"Age_comp "<<age_like*column(age_lambda,k)<<endl;
  if(nobs_ms_tot>0) SS2out <<"Size_at_age "<<sizeage_like*column(sizeage_lambda,k)<<endl;
  if(SzFreq_Nmeth>0) SS2out <<"SizeFreq "<<SzFreq_like*column(SzFreq_lambda,k)<<endl;
  if(Do_Morphcomp>0) SS2out <<"Morphcomp "<<Morphcomp_lambda(k)*Morphcomp_like<<" "<<Morphcomp_lambda(k)<<endl;
  if(Do_TG>0) SS2out <<"Tag_comp "<<TG_like1*column(TG_lambda1,k)<<endl;
  if(Do_TG>0) SS2out <<"Tag_negbin "<<TG_like2*column(TG_lambda2,k)<<endl;
  SS2out <<"Recruitment "<<recr_like*recrdev_lambda(k)<<" "<<recrdev_lambda(k)<<endl;
  SS2out <<"Forecast_Recruitment "<<Fcast_recr_like<<" "<<Fcast_recr_lambda<<endl;
  SS2out <<"Parm_priors "<<parm_like*parm_prior_lambda(k)<<" "<<parm_prior_lambda(k)<<endl;
  if(SoftBound>0) SS2out <<"Parm_softbounds "<<SoftBoundPen<<" "<<" NA "<<endl;
  SS2out <<"Parm_devs "<<parm_dev_like*parm_dev_lambda(k)<<" "<<parm_dev_lambda(k)<<endl;
  if(F_ballpark_yr>0) SS2out <<"F_Ballpark "<<F_ballpark_lambda(k)*F_ballpark_like<<" "<<F_ballpark_lambda(k)<<"  ##:est&obs: "<<annual_F(F_ballpark_yr,2)<<" "<<F_ballpark<<endl;
  SS2out <<"Crash_Pen "<<CrashPen_lambda(k)*CrashPen<<" "<<CrashPen_lambda(k)<<endl;

  SS2out<<"_"<<endl<<"Fleet:  ALL ";
  for(f=1;f<=Ntypes;f++) SS2out<<f<<" ";
  SS2out<<endl;
  if(F_Method>1) SS2out<<"Catch_lambda: _ "<<column(catch_lambda,k)<<endl<<"Catch_like: "<<catch_like*column(catch_lambda,k) <<" "<<catch_like<<endl;
  if(nobs_cr>0) SS2out<<"Surv_lambda: _ "<<column(surv_lambda,k)<<endl<<"Surv_like: "<<surv_like*column(surv_lambda,k)<<" "<<surv_like<<endl;
  if(nobs_disc>0) SS2out<<"Disc_lambda: _ "<<column(disc_lambda,k)<<endl<<"Disc_like: "<<disc_like*column(disc_lambda,k)<<" "<<disc_like<<endl;
  if(nobs_mnwt>0) SS2out<<"mnwt_lambda: _ "<<column(mnwt_lambda,k)<<endl<<"mnwt_like: "<<mnwt_like*column(mnwt_lambda,k)<<" "<<mnwt_like<<endl;
  if(Nobs_l_tot>0) SS2out<<"Length_lambda: _ "<<column(length_lambda,k)<<endl<<"Length_like: "<<length_like*column(length_lambda,k)<<" "<<length_like<<endl;
  if(Nobs_a_tot>0) SS2out<<"Age_lambda: _ "<<column(age_lambda,k)<<endl<<"Age_like: "<<age_like*column(age_lambda,k)<<" "<<age_like<<endl;
  if(nobs_ms_tot>0) SS2out<<"Sizeatage_lambda: _ "<<column(sizeage_lambda,k)<<endl<<"sizeatage_like: "<<sizeage_like*column(sizeage_lambda,k)<<" "<<sizeage_like<<endl;

  if(SzFreq_Nmeth>0)
  {
    for(j=1;j<=SzFreq_Nmeth;j++)
    {
      SS2out<<"SizeFreq_lambda:_"<<j<<"; ";
    if(j==1) {SS2out<<"_ ";} else {SS2out<<"_ ";}
      for(f=1;f<=Ntypes;f++)
      {
        if(SzFreq_LikeComponent(f,j)>0) {SS2out<<SzFreq_lambda(SzFreq_LikeComponent(f,j),k)<<" ";}
        else {SS2out<<" NA ";}
      }
      SS2out<<endl;
      SS2out<<"SizeFreq_like:_"<<j<<"; ";
    if(j==1) {SS2out<<SzFreq_like*column(SzFreq_lambda,k)<<" ";} else {SS2out<<"_ ";}
      for(f=1;f<=Ntypes;f++)
      {
        if(SzFreq_LikeComponent(f,j)>0) {SS2out<<SzFreq_like(SzFreq_LikeComponent(f,j))<<" ";}
        else {SS2out<<" NA ";}
      }
      SS2out<<endl;
    }
    SS2out<<SzFreq_like<<endl<<SzFreq_like_base<<endl;
  }

  if(Do_TG>0)
  {
    SS2out<<endl<<"_"<<endl<<"Tag_Group:  ALL ";
    for(f=1;f<=N_TG;f++) SS2out<<f<<" ";
    SS2out<<endl;
    SS2out<<"Tag_comp_Lambda _ "<<column(TG_lambda1,k)<<endl<<
    "Tag_comp_Like "<<TG_like1*column(TG_lambda1,k)<<" "<<TG_like1<<endl;
    SS2out<<"Tag_negbin_Lambda _ "<<column(TG_lambda2,k)<<endl<<
    "Tag_negbin_Like "<<TG_like2*column(TG_lambda2,k)<<" "<<TG_like2<<endl;
  }
  SS2out<<endl;

  SS2out<<"Input_Variance_Adjustment"<<endl<<"Fleet ";
  for(i=1;i<=Ntypes;i++) {SS2out<<" "<<i;}
  SS2out<<endl;
  SS2out << "Index_extra_CV "<<var_adjust(1)<<endl;
  SS2out << "Discard_extra_CV "<<var_adjust(2)<<endl;
  SS2out << "MeanBodyWt_extra_CV "<<var_adjust(3)<<endl;
  SS2out << "effN_mult_Lencomp "<<var_adjust(4)<<endl;
  SS2out << "effN_mult_Agecomp "<<var_adjust(5)<<endl;
  SS2out << "effN_mult_Len_at_age "<<var_adjust(6)<<endl;

  SS2out<<"MG_parms"<<"Using_offset_approach_#:_"<<MGparm_def<<"  (1=none, 2= M, G, CV_G as offset from female_GP1, 3=like SS2 V1.x)"<<endl;

//  SS2out<<endl<<"PARAMETERS"<<endl<<"Num Label Value Active_Cnt Phase Min Max Init Prior PR_type Pr_SD Prior_Like Parm_StDev Status Pr_atMin Pr_atMax"<<endl;
  SS2out<<endl<<"PARAMETERS"<<endl<<"Num Label Value Active_Cnt  Phase Min Max Init  Status  Parm_StDev PR_type Prior Pr_SD Prior_Like Value_again Value-1.96*SD Value+1.96*SD V_1%  V_10% V_20% V_30% V_40% V_50% V_60% V_70% V_80% V_90% V_99% P_val P_lowCI P_hiCI  P_1%  P_10% P_20% P_30% P_40% P_50% P_60% P_70% P_80% P_90% P_99%"<<endl;

  NP=0;   // count of number of parameters
  active_count=0;
  Nparm_on_bound=0;
  int Activ;
  for(j=1;j<=N_MGparm2;j++)
  {
    NP++;
    Activ=0;
    if(active(MGparm(j)))
    {
      active_count++;
      Activ=1;
    }
    Report_Parm(NP, active_count, Activ, MGparm(j), MGparm_LO(j), MGparm_HI(j), MGparm_RD(j), MGparm_PR(j), MGparm_PRtype(j), MGparm_CV(j), MGparm_PH(j), MGparm_Like(j));
  }

  if(N_MGparm_dev>0)
  {
    for(i=1;i<=N_MGparm_dev;i++)
    for(j=MGparm_dev_minyr(i);j<=MGparm_dev_maxyr(i);j++)
    {
      NP++;  SS2out<<NP<<" "<<ParmLabel(NP)<<" "<<MGparm_dev(i,j);
      if(active(MGparm_dev))
      {
        active_count++;
        SS2out<<" "<<active_count<<" "<<MGparm_dev_PH<<" _ _ _ act "<<CoVar(active_count,1);
      }
      else
      {SS2out<<" _ _ _ _ _ NA _ ";}
      SS2out<<" dev "<<endl;
    }
  }

  for(j=1;j<=N_SRparm2;j++)
  {
    NP++;
    Activ=0;
    if(active(SR_parm(j)))
    {
      active_count++;
      Activ=1;
    }
    Report_Parm(NP, active_count, Activ, SR_parm(j), SR_parm_1(j,1), SR_parm_1(j,2), SR_parm_1(j,3), SR_parm_1(j,4), SR_parm_1(j,5), SR_parm_1(j,6), SR_parm_1(j,7), SR_parm_Like(j));
  }

  if(recdev_cycle>0)
  {
    for(j=1;j<=recdev_cycle;j++)
    {
      NP++;
      Activ=0;
      if(active(recdev_cycle_parm(j)))
      {
        active_count++;
        Activ=1;
      }
      Report_Parm(NP, active_count, Activ, recdev_cycle_parm(j), recdev_cycle_parm_RD(j,1), recdev_cycle_parm_RD(j,2), recdev_cycle_parm_RD(j,3), recdev_cycle_parm_RD(j,4), recdev_cycle_parm_RD(j,5), recdev_cycle_parm_RD(j,6), recdev_cycle_parm_RD(j,7), recdev_cycle_Like(j));
    }
  }

    if(recdev_do_early>0)
      {
        for(i=recdev_early_start;i<=recdev_early_end;i++)
        {NP++;  SS2out<<NP<<" "<<ParmLabel(NP)<<" "<<recdev(i);
        if( active(recdev_early) )
        {
          active_count++;
          SS2out<<" "<<active_count<<" _ _ _ _ act "<<CoVar(active_count,1);
        }
        else
          {
            SS2out<<" _ _ _ _ _ NA _ ";
          }
        SS2out <<" dev "<<endl;
          }
      }

    if(do_recdev>0)
      {
        for(i=recdev_start;i<=recdev_end;i++)
        {NP++;  SS2out<<NP<<" "<<ParmLabel(NP)<<" "<<recdev(i);
        if( active(recdev1)||active(recdev2) )
        {
          active_count++;
          SS2out<<" "<<active_count<<" _ _ _ _ act "<<CoVar(active_count,1);
        }
        else
          {
            SS2out<<" _ _ _ _ _ NA _ ";
          }
        SS2out <<" dev "<<endl;
          }
      }

    if(Do_Forecast>0)
    {
      for(i=recdev_end+1;i<=YrMax;i++)
      {
        NP++; SS2out<<NP<<" "<<ParmLabel(NP)<<" "<<Fcast_recruitments(i);
        if(active(Fcast_recruitments))
        {active_count++;           SS2out<<" "<<active_count<<" _ _ _ _ act "<<CoVar(active_count,1);}
        else
        {SS2out<<"  _ _ _ _ _ NA _ ";}
        SS2out <<" dev "<<endl;
      }
    }

      if(Do_Forecast>0)
      {
        for(i=endyr+1;i<=YrMax;i++)
        {
          NP++; SS2out<<NP<<" "<<ParmLabel(NP)<<" "<<Fcast_impl_error(i);
          if(active(Fcast_impl_error))
          {active_count++;           SS2out<<" "<<active_count<<" _ _ _ _ act "<<CoVar(active_count,1);}
          else
          {SS2out<<"  _ _ _ _ _ NA _ ";}
        SS2out <<" dev "<<endl;
        }
      }

  for(j=1;j<=Nfleet;j++)
  {
    NP++;
    Activ=0;
    if(active(init_F(j)))
    {
      active_count++;
      Activ=1;
    }
    Report_Parm(NP, active_count, Activ, init_F(j), init_F_LO(j), init_F_HI(j), init_F_RD(j), init_F_PR(j), init_F_PRtype(j), init_F_CV(j), init_F_PH(j), init_F_Like(j));
  }

    if(F_Method==2)
    {
      for(i=1;i<=Nfleet*(TimeMax-styr+1);i++)
      {
        NP++;  SS2out<<NP<<" "<<ParmLabel(NP)<<" "<<F_rate(i);
        if(active(F_rate(i)))
        {
          active_count++;
          SS2out<<" "<<active_count<<" "<<have_catch(i)<<" 0.0  8.0  _ act "<<CoVar(active_count,1);
        }
        else
        {SS2out<<" _ _ _ _ _ NA _ ";}
        SS2out <<" F "<<endl;
      }
    }

  for(j=1;j<=Q_Npar;j++)
  {
    NP++;
    Activ=0;
    if(active(Q_parm(j)))
    {
      active_count++;
      Activ=1;
    }
    Report_Parm(NP, active_count, Activ, Q_parm(j), Q_parm_1(j,1), Q_parm_1(j,2), Q_parm_1(j,3), Q_parm_1(j,4), Q_parm_1(j,5), Q_parm_1(j,6), Q_parm_1(j,7), Q_parm_Like(j));
  }

  for(j=1;j<=N_selparm2;j++)
  {
    NP++;
    Activ=0;
    if(active(selparm(j)))
    {
      active_count++;
      Activ=1;
    }
    Report_Parm(NP, active_count, Activ, selparm(j), selparm_LO(j), selparm_HI(j), selparm_RD(j), selparm_PR(j), selparm_PRtype(j), selparm_CV(j), selparm_PH(j), selparm_Like(j));
  }

    for(i=1;i<=N_selparm_dev;i++)
    for(j=selparm_dev_minyr(i);j<=selparm_dev_maxyr(i);j++)
      {
        NP++;  SS2out<<NP<<" "<<ParmLabel(NP)<<" "<<selparm_dev(i,j);
        if(active(selparm_dev))
          {active_count++; SS2out<<" "<<active_count<<" _ _ _ _ act "<<CoVar(active_count,1);}
          else
          {SS2out<<" _ _ _ _ _ NA _ ";}
        SS2out <<" dev "<<endl;
    }

  if(Do_TG>0)
  {
     k=3*N_TG+2*Nfleet;
    for(j=1;j<=k;j++)
    {
      NP++;
      Activ=0;
      if(active(TG_parm(j)))
      {
        active_count++;
        Activ=1;
      }
      Report_Parm(NP, active_count, Activ, TG_parm(j), TG_parm_LO(j), TG_parm_HI(j), TG_parm2(j,3), TG_parm2(j,4), TG_parm2(j,5), TG_parm2(j,6), TG_parm_PH(j), TG_parm_Like(j));
    }
  }
  
  SS2out<<endl<<"Number_of_active_parameters_on_or_near_bounds: "<<Nparm_on_bound<<endl;
  SS2out<<"Active_count "<<active_count<<endl<<endl;
  SS2out<<endl<<"DERIVED_QUANTITIES"<<endl;
  SS2out<<"SPR_ratio_basis: "<<SPR_report_label<<endl;
  SS2out<<"F_report_basis: "<<F_report_label<<endl;
  SS2out<<"B_ratio_denominator: "<<depletion_basis_label<<endl;

  SS2out<<" LABEL Value  StdDev (Val-1.0)/Stddev  CumNorm"<<endl;
  for(j=1;j<=N_STD_Yr;j++)
  {
    NP++;  SS2out<<" "<<ParmLabel(NP)<<" "<<SPB_std(j);
    active_count++;
    SS2out<<" "<<CoVar(active_count,1)<<endl;
  }
  for(j=1;j<=N_STD_Yr;j++)
  {
    NP++;  SS2out<<" "<<ParmLabel(NP)<<" "<<recr_std(j);
    active_count++;
    SS2out<<" "<<CoVar(active_count,1)<<endl;
  }
  for(j=1;j<=N_STD_Yr_Ofish;j++)
  {
    NP++;  SS2out<<" "<<ParmLabel(NP)<<" "<<SPR_std(j);
    active_count++;
    SS2out<<" "<<CoVar(active_count,1);
    if( CoVar(active_count,1)>0.0)
    {
      temp=value((SPR_std(j)-1.0)/CoVar(active_count,1));
      SS2out<<" "<<temp<<" "<<cumd_norm(temp);

    }
    SS2out<<endl;
  }
  ofstream post_vecs("posterior_vectors.sso",ios::app);
  post_vecs<<runnumber<<" 0 "<<obj_fun<<" F/Fmsy_stdev ";
  for(j=1;j<=N_STD_Yr_F;j++)
  {
    NP++;  SS2out<<" "<<ParmLabel(NP)<<" "<<F_std(j);
    active_count++;
    SS2out<<" "<<CoVar(active_count,1);
    post_vecs<<CoVar(active_count,1)<<" ";
    if( CoVar(active_count,1)>0.0)
    {
      temp=value((F_std(j)-1.0)/CoVar(active_count,1));
      SS2out<<" "<<temp<<" "<<cumd_norm(temp);
    }
    SS2out<<endl;
  }
  post_vecs<<endl;
  post_vecs<<runnumber<<" 0 "<<obj_fun<<" B/Bmsy_stdev ";
 
  for(j=1;j<=N_STD_Yr_Dep;j++)
  {
    NP++;  SS2out<<" "<<ParmLabel(NP)<<" "<<depletion(j);
    active_count++;
    SS2out<<" "<<CoVar(active_count,1);
    post_vecs<<CoVar(active_count,1)<<" ";
    if( CoVar(active_count,1)>0.0)
    {
      temp=value((depletion(j)-1.0)/CoVar(active_count,1));
      SS2out<<" "<<temp<<" "<<cumd_norm(temp);
    }
    SS2out<<endl;
  }
  post_vecs<<endl;
  for(j=1;j<=N_STD_Mgmt_Quant;j++)
  {
    NP++; active_count++;
    SS2out<<" "<<ParmLabel(NP)<<" "<<Mgmt_quant(j);
    
    SS2out<<" "<<CoVar(active_count,1)<<endl;
  }

  for(j=1;j<=Extra_Std_N;j++)
  {
    NP++;      active_count++;
    SS2out<<" "<<ParmLabel(NP)<<" "<<Extra_Std(j);
    SS2out<<" "<<CoVar(active_count,1)<<endl;
  }

 /*
    SS2out<<endl<<"MGParm_Block_Assignments"<<endl;
    if(N_MGparm_blk>0)
    {
      SS2out<<"Base_parm# ";
      for(y=styr;y<=endyr;y++)
      {SS2out<<" "<<y;}
      SS2out<<endl;
      for(j=1;j<=N_MGparm;j++)
      {
        if(MGparm_1(j,13)>0) SS2out<<j<<" "<<Block_Defs_MG(j)<<endl;
      }
    }

    SS2out<<endl<<"Selex_Block_Assignments"<<endl;
    if(N_selparm_blk>0)
    {
      SS2out<<"Base_parm# ";
      for(y=styr;y<=endyr;y++) {SS2out<<" "<<y;}
      SS2out<<endl;
      for(j=1;j<=N_selparm;j++)
      {
        if(selparm_1(j,13)>0)
        SS2out<<j<<" "<<Block_Defs_Sel(j)<<endl;
      }
    }
 */
   if(reportdetail>0) {k1=endyr;} else {k1=styr;}
   SS2out<<endl<<"MGparm_By_Year_after_adjustments"<<endl<<"Year ";
   for(i=1;i<=N_MGparm;i++) SS2out<<" "<<ParmLabel(i);
   SS2out<<endl;
   for(y=styr;y<=k1;y++)
     SS2out<<y<<" "<<mgp_save(y)<<endl;

   SS2out<<endl<<"selparm(Size)_By_Year_after_adjustments"<<endl<<"Fleet/Svy Year"<<endl;
   for(f=1;f<=Ntypes;f++)
   for(y=styr;y<=k1;y++)
     {
     k=N_selparmvec(f);
     if(k>0) SS2out<<f<<" "<<y<<" "<<save_sp_len(y,f)(1,k)<<endl;
     }

   SS2out<<endl<<"selparm(Age)_By_Year_after_adjustments"<<endl<<"Fleet/Svy Year"<<endl;
   for(f=Ntypes+1;f<=2*Ntypes;f++)
   for(y=styr;y<=k1;y++)
     {
     k=N_selparmvec(f);
     if(k>0) SS2out<<f-Ntypes<<" "<<y<<" "<<save_sp_len(y,f)(1,k)<<endl;
     }

   SS2out<<endl<<"RECRUITMENT_DIST_EndYear: "<<endyr<<endl<<"G_pattern Gender Seas Area Value Used?"<<endl;
   g=0;
   for(gg=1;gg<=gender;gg++)
   for(gp=1;gp<=N_GP;gp++)
   {
    g++;
    for(s=1;s<=nseas;s++)
    for(p=1;p<=pop;p++)
    {
       SS2out<<gp<<" "<<gg<<" "<<s<<" "<<p<<" "<<recr_dist_calc(g,s,p)<<" "<<recr_dist_pattern(gp,s,p)<<endl;
    }
  }

   SS2out<<endl<<"RECRUITMENT_DIST_BENCHMARK_from_years: "<<Bmark_Yr(1)<<" "<<Bmark_Yr(2)<<endl<<"G_pattern Gender Seas Area Value Used?"<<endl;
   g=0;
   for(gg=1;gg<=gender;gg++)
   for(gp=1;gp<=N_GP;gp++)
   {
    g++;
    for(s=1;s<=nseas;s++)
    for(p=1;p<=pop;p++)
    {
       SS2out<<gp<<" "<<gg<<" "<<s<<" "<<p<<" "<<recr_dist_benchmark(g,s,p)<<" "<<recr_dist_pattern(gp,s,p)<<endl;
    }
  }
   SS2out<<endl<<"RECRUITMENT_DIST_FORECAST_from_years: "<<Fcast_yr(1)<<" "<<Fcast_yr(2)<<endl<<"G_pattern Gender Seas Area Value Used?"<<endl;
   g=0;
   for(gg=1;gg<=gender;gg++)
   for(gp=1;gp<=N_GP;gp++)
   {
    g++;
    for(s=1;s<=nseas;s++)
    for(p=1;p<=pop;p++)
    {
       SS2out<<gp<<" "<<gg<<" "<<s<<" "<<p<<" "<<recr_dist_forecast(g,s,p)<<" "<<recr_dist_pattern(gp,s,p)<<endl;
    }
  }

   SS2out<<endl<<"MORPH_INDEXING"<<endl;
   SS2out<<"Index Gpattern Gender Bseas Sub_Morph Sub_Morph_Dist Gpattern_Gender Gpattern_Gender_Bseas BirthAge_Rel_Jan1"<<endl;
   for(g=1; g<=gmorph; g++)
   {
     SS2out<<g<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<submorphdist(GP2(g))<<" "<<GP(g)<<" "<<GP3(g)<<" "<<azero_G(g)<<endl;
   }

//  3darray SzFreqTrans(1,SzFreq_Nmeth*nseas,1,nlength2,1,SzFreq_Nbins_seas_g);
   if(SzFreq_Nmeth>0)
   {
     SS2out<<endl<<"SIZEFREQ_TRANSLATION "<<SzFreq_scale<<endl;
     for(SzFreqMethod=1;SzFreqMethod<=SzFreq_Nmeth;SzFreqMethod++)
     {
       SS2out<<SzFreqMethod<<" seas len mid-len ";
           if(SzFreq_scale(SzFreqMethod)==1)
           {SS2out<<" mid-kg ";}
           else if(SzFreq_scale(SzFreqMethod)==2)
           {SS2out<<" mid-lbs ";}
           else if(SzFreq_scale(SzFreqMethod)==3)
           {SS2out<<" mid-cm ";}
           else
           {SS2out<<" mid-inch ";}
       SS2out<<SzFreq_bins1(SzFreqMethod);
       if(gender==2) SS2out<<SzFreq_bins1(SzFreqMethod);
       SS2out<<endl<<SzFreqMethod<<" seas len mid-len metric "<<SzFreq_bins(SzFreqMethod)<<endl;;
       for(s=1;s<=nseas;s++)
       {
         SzFreqMethod_seas=nseas*(SzFreqMethod-1)+s;     // index that combines sizefreqmethod and season and used in SzFreqTrans
         for(z=1;z<=nlength2;z++)
         {
           SS2out<<SzFreqMethod<<" "<<s<<" "<<len_bins2(z)<<" "<<len_bins_m2(z)<<" ";
           if(SzFreq_scale(SzFreqMethod)==1)
           {SS2out<<wt_len2(s,z)<<" ";}
           else if(SzFreq_scale(SzFreqMethod)==2)
           {SS2out<<wt_len2(s,z)/0.4536<<" ";}
           else if(SzFreq_scale(SzFreqMethod)==3)
           {SS2out<<len_bins_m2(z)<<" ";}
           else
           {SS2out<<len_bins_m2(z)/2.54<<" ";}
           for(j=1;j<=gender*SzFreq_Nbins(SzFreqMethod);j++)
           {
             SS2out<<SzFreqTrans(SzFreqMethod_seas,z,j)<<" ";
             if(SzFreqTrans(SzFreqMethod_seas,z,j)<0.0)
             {N_warn++; warning<<" CRITICAL ERROR:  Bin widths narrower than pop len bins caused negative allocation in sizefreq method:"<<
            " method, season, size, bin: "<<SzFreqMethod<<" "<<s<<" "<<len_bins2(z)<<" "<<j<<endl; exit(1);}
           }
           SS2out<<endl;
         }
       }
     }
   }

   SS2out<<"#"<<endl<<"MOVEMENT in endyear"<<endl<<" Seas Gpattern Source_area Dest_area minage maxage "<<age_vector<<endl;
   for(k=1;k<=do_migr2;k++)
   {
     SS2out<<move_def2(k)<<" "<<migrrate(endyr,k)<<endl;
   }

   SS2out<<endl<<"EXPLOITATION"<<endl<<"F_Method: "<<F_Method;
   if(F_Method==1) {SS2out<<"  Pope's_approx ";} else {SS2out<<"  Continuous_F;_(NOTE:_F_std_adjusts_for_seasdur_but_each_fleet_F_is_annual)";}
   SS2out<<endl<<"F_std_units: "<<F_reporting<<F_report_label<<endl<<"_ _ _ ";
   for(f=1;f<=Nfleet;f++) if(catchunits(f)==1) {SS2out<<" Bio ";} else {SS2out<<" Num ";}
   SS2out<<endl<<"_ _ _ ";
   for(f=1;f<=Nfleet;f++) SS2out<<" "<<f;
   SS2out<<endl<<"Yr Seas F_std";
   for(f=1;f<=Nfleet;f++) SS2out<<" "<<fleetname(f);
   SS2out<<endl<<"init_yr _ _ "<<init_F<<endl;
   for(y=styr;y<=YrMax;y++)
   for(s=1;s<=nseas;s++)
   {
     t=styr+(y-styr)*nseas+s-1;
     SS2out<<y<<" "<<s<<" ";
     if(s==1 && y>=styr && STD_Yr_Reverse_F(y)>0 ) {SS2out<<F_std(STD_Yr_Reverse_F(y));} else {SS2out<<" _ ";}
     SS2out<<" "<<column(Hrate,t)<<endl;
   }

  SS2out<<endl<<"CATCH "<<endl<<"Fleet Name Yr Seas Yr.S Obs Exp se F  Like sel_bio kill_bio ret_bio sel_num kill_num ret_num"<<endl;
  for(f=1;f<=Nfleet;f++)
  for (y=styr;y<=endyr;y++)
  for (s=1;s<=nseas;s++)
  {
    t = styr+(y-styr)*nseas+s-1;
    if(catchunits(f)==1)
    {gg=3;}  //  biomass
    else
    {gg=6;}  //  numbers
    temp = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
    SS2out<<f<<" "<<fleetname(f)<<" "<<y<<" "<<s<<" "<<temp<<" "<<catch_ret_obs(f,t)<<" "<<catch_fleet(t,f,gg)<<" "<<catch_se(t,f)<<" "<<Hrate(f,t)<<" ";
    if(catch_se_rd(f)>0.0)
      {SS2out<<0.5*square( (log(1.1*catch_ret_obs(f,t)) -log(catch_fleet(t,f,gg)+0.1*catch_ret_obs(f,t))) / catch_se(t,f));}
      else
      {SS2out<<"BYCATCH";}
      SS2out<<catch_fleet(t,f)<<endl;
  }

   int bio_t;

   SS2out<<endl<<"TIME_SERIES    BioSmry_age:_"<<Smry_Age;   // ProgLabel_320
   if(F_Method==1) {SS2out<<"  Pope's_approx"<<endl;} else {SS2out<<"  Continuous_F"<<endl;}
  SS2out<<"Area Yr Era Seas Bio_all Bio_smry SpawnBio Recruit_0 ";
  for(gp=1;gp<=N_GP;gp++) SS2out<<" Spbio_GP:"<<gp;
  if(Hermaphro_Option>0)
  {
    for(gp=1;gp<=N_GP;gp++) SS2out<<" MaleSpbio_GP:"<<gp;
  }
  for(gg=1;gg<=gender;gg++)
  for(gp=1;gp<=N_GP;gp++)
  {SS2out<<" SmryBio_SX:"<<gg<<"_GP:"<<gp;}
  for(gg=1;gg<=gender;gg++)
  for(gp=1;gp<=N_GP;gp++)
  {SS2out<<" SmryNum_SX:"<<gg<<"_GP:"<<gp;}
  dvector Bio_Comp(1,N_GP*gender);
  dvector Num_Comp(1,N_GP*gender);
  for (f=1;f<=Nfleet;f++)
  {
    SS2out<<" sel(B):_"<<f<<" dead(B):_"<<f<<" retain(B):_"<<f<<
    " sel(N):_"<<f<<" dead(N):_"<<f<<" retain(N):_"<<f<<
    " obs_cat:_"<<f;
     if(F_Method==1) {SS2out<<" Hrate:_"<<f;} else {SS2out<<" F:_"<<f;}
  }

  SS2out<<" SPB_vir_LH"<<endl;

  for (p=1;p<=pop;p++)
  {
   for (y=styr-2;y<=YrMax;y++)
   {
    if(y<=endyr && p==1) {Smry_Table(y)(1,8).initialize();  Smry_Table(y)(15,17).initialize();}
    for(s=1;s<=nseas;s++)
    {
    t = styr+(y-styr)*nseas+s-1;
    bio_t=t;
    if(y<=styr) {bio_t=styr-1+s;}
    Bio_Comp.initialize();
    Num_Comp.initialize();
    totbio.initialize(); smrybio.initialize(); smrynum.initialize(); SPB_vir_LH.initialize(); smryage.initialize();
    Recr(p,y)=0;
    for(g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
    {
     if(s==Bseas(g)) Recr(p,y)+=natage(t,p,g,0);
     gg=sx(g);
     temp=natage(t,p,g)(Smry_Age,nages)*Save_Wt_Age(bio_t,g)(Smry_Age,nages);
     Bio_Comp(GP(g))+=value(temp);   //sums to accumulate across submorphs and birthseasons
     Num_Comp(GP(g))+=value(sum(natage(t,p,g)(Smry_Age,nages)));   //sums to accumulate across submorphs and birthseasons
     totbio+= natage(t,p,g)*Save_Wt_Age(bio_t,g);
     smrybio+= temp;
     smrynum+=sum(natage(t,p,g)(Smry_Age,nages));
     smryage+=natage(t,p,g)(Smry_Age,nages)*r_ages(Smry_Age,nages);
     SPB_vir_LH += natage(t,p,g)*virg_fec(g);
     if(y<=endyr)
     {
       for (f=1;f<=Nfleet;f++)
       {
         if(fleet_area(f)==p&&y>=styr-1)
         {
           Smry_Table(y,16)+=sum(catage(t,f,g));
           Smry_Table(y,17)+=catage(t,f,g)*r_ages;
         }
       }
     }


    } //close gmorph loop
    SS2out<<p<<" "<<y;
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}

    SS2out<<s<<" "<<totbio<<" "<<smrybio<<" ";
    if(s==spawn_seas)
    {
      temp=sum(SPB_pop_gp(y,p));
      if(Hermaphro_maleSPB==1) temp+=sum(MaleSPB(y,p));
      SS2out<<temp;
    }
    else
    {SS2out<<" _ ";}
    SS2out<<" "<<Recr(p,y)<<" ";
    if(s==spawn_seas)
    {
      SS2out<<SPB_pop_gp(y,p);
      if(Hermaphro_Option>0) SS2out<<MaleSPB(y,p);
    }
    else
    {
    for(gp=1;gp<=N_GP;gp++) {SS2out<<" _ ";}

    }
    SS2out<<" "<<Bio_Comp<<" "<<Num_Comp;
    if(s==1 && y<=endyr) {Smry_Table(y,1)+=totbio; Smry_Table(y,2)+=smrybio;  Smry_Table(y,3)+=smrynum; Smry_Table(y,15)+=smryage;}  // already calculated for the forecast years
    Smry_Table(y,7)=SPB_yr(y);
    Smry_Table(y,8)=exp_rec(y,4);
    for (f=1;f<=Nfleet;f++)
    {
      if(fleet_area(f)==p&&y>=styr-1)
      {
        SS2out<<" "<<catch_fleet(t,f)<<" ";
        if(y<=endyr) {SS2out<<catch_ret_obs(f,t)<<" "<<Hrate(f,t);} else {SS2out<<" _ "<<Hrate(f,t);}
        if(y<=endyr) {Smry_Table(y,4)+=catch_fleet(t,f,1); Smry_Table(y,5)+=catch_fleet(t,f,2); Smry_Table(y,6)+=catch_fleet(t,f,3);}
        
      }
      else
      {SS2out<<" 0 0 0 0 0 0 0 0 ";}
    }
    if(s==spawn_seas)
        {SS2out<<" "<<SPB_vir_LH<<endl;}
    else
      {SS2out<<" _"<<endl;}
    }
   }
  }

    // start SPR time series                                  ProgLabel_0322
   SS2out<<endl<<"SPR_series_uses_R0= "<<Recr_virgin<<endl<<"###note_Y/R_unit_is_Dead_Biomass"<<endl;
   SS2out<<"Depletion_method: "<<depletion_basis<<" # "<<depletion_basis_label<<endl;
   SS2out<<"F_std_method: "<<F_reporting<<" # "<<F_report_label<<endl;
   SS2out<<"SPR_std_method: "<<SPR_reporting<<" # "<<SPR_report_label<<endl;
   // note  GENTIME is mean age of spawners weighted by fec(a)
   SS2out<<"Year Bio_all Bio_Smry SPBzero SPBfished SPBfished/R SPR SPR_std Y/R GenTime Deplete F_std"<<
   " Actual: Bio_all Bio_Smry Num_Smry MnAge_Smry Enc_Catch Dead_Catch Retain_Catch MnAge_Catch SPB Recruits Tot_Exploit"<<
   " More_F(by_Morph): ";
   for(g=1;g<=gmorph;g++) {SS2out<<" aveF_"<<g;}
   for(g=1;g<=gmorph;g++) {SS2out<<" maxF_"<<g;}
   SS2out<<" Enc_Catch_B Dead_Catch_B Retain_Catch_B  Enc_Catch_N Dead_Catch_N Retain_Catch_N sum_Apical_F F=Z-M";
   SS2out<<endl;

   for (y=styr;y<=YrMax;y++)
   {
     SS2out<<y<<" "<<Smry_Table(y)(9,12)<<" "<<(Smry_Table(y,12)/Recr_virgin)<<" "<<Smry_Table(y,12)/Smry_Table(y,11)<<" ";
     if(STD_Yr_Reverse_Ofish(y)>0) {SS2out<<SPR_std(STD_Yr_Reverse_Ofish(y))<<" ";} else {SS2out<<" _ ";}
     SS2out<<(Smry_Table(y,14)/Recr_virgin)<<" "<<Smry_Table(y,13)<<" ";
     if(STD_Yr_Reverse_Dep(y)>0) {SS2out<<depletion(STD_Yr_Reverse_Dep(y));} else {SS2out<<" _ ";}
     if(y>=styr && STD_Yr_Reverse_F(y)>0 ) {SS2out<<" "<<F_std(STD_Yr_Reverse_F(y));} else {SS2out<<" _ ";}
     SS2out<<" & "<<Smry_Table(y)(1,3)<<" "<<Smry_Table(y,15)/Smry_Table(y,3)<<" "<<Smry_Table(y)(4,6)<<" "<<Smry_Table(y,17)/(Smry_Table(y,16)+1.0e-06)<<" "<<Smry_Table(y)(7,8)<<" "<<Smry_Table(y,5)/Smry_Table(y,2);
     SS2out<<" & "<<Smry_Table(y)(21,20+gmorph)<<" "<<Smry_Table(y)(21+gmorph,20+2*gmorph)<<" "<<annual_catch(y)<<" "<<annual_F(y)<<endl;
   } // end year loop

  SIS_table<<"to_come_#"  <<" Best_F"<<endl;
  SIS_table<<"to_come_#"  <<"  F_units"<<endl;
  SIS_table<<endyr<<" F_year"<<endl;
  SIS_table<<"to_come_#"  <<" F_basis"<<endl;
  SIS_table<<"to_come_#"  <<" F_limit"<<endl;
  SIS_table<<"to_come_#"  <<"  F_limit_basis"<<endl;
  SIS_table<<"to_come_#"  <<" F_msy"<<endl;
  SIS_table<<"to_come_#"  <<" F_msy_basis"<<endl;
  SIS_table<<"to_come_#"  <<" F_target"<<endl;
  SIS_table<<"to_come_#"  <<" F_target_basis"<<endl;
  SIS_table<<"to_come_#"  <<" F/F_limit"<<endl;
  SIS_table<<"to_come_#"  <<" F/F_msy"<<endl;
  SIS_table<<"to_come_#"  <<" F/F_target"<<endl;
  SIS_table<<"to_come_#"  <<" Best_Bio"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_units"<<endl;
  SIS_table<<endyr<<" Bio_year"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_basis"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_limit"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_limit_basis"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_MSY"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_MSY_basis"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_MSY"<<endl;
  SIS_table<<"to_come_#"  <<" Bio_above_below"<<endl;
  SIS_table<<"to_come_#"  <<" Bio/Bio_limit"<<endl;
  SIS_table<<"to_come_#"  <<" Bio/Bio_MSY"<<endl<<endl;

 /*
  SIS_table<<"Category Year Abundance Abundance Recruitment Spawners Fmort Fmort Catch Catch"<<endl;
  SIS_table<<"Primary _ N Y Y Y Y Y N N"<<endl;
  SIS_table<<"Type _ Biomass Biomass Age Female_Mature Exploitation ";
  SIS_table<<SPR_report_label<<" Catch Catch"<<endl;
  SIS_table<<"Source _ Model Model Model Model Model Model Model Model"<<endl;
  SIS_table<<"Basis _ Biomass Biomass Numbers Eggs Rate Rate Biomass Biomass"<<endl;
  SIS_table<<"Range _ All Age_"<<Smry_Age<<"+ Age_0 Mature Rate Rate All Retained"<<endl;
  SIS_table<<"Statistic _ Mean Mean Mean Mean Mean Mean Mean Mean"<<endl;
  SIS_table<<"Scalar _ 1 1 1000 1 1 1 1 1"<<endl;

   for (y=styr;y<=YrMax;y++)
   {
     SIS_table<<"_ "<<y<<" "<<Smry_Table(y,1)<<" "<<Smry_Table(y,2)<<" "<<Smry_Table(y,8)<<" "<<Smry_Table(y,7)<<" "<<Smry_Table(y,5)/Smry_Table(y,2)<<" ";
     if(STD_Yr_Reverse_Ofish(y)>0) {SIS_table<<SPR_std(STD_Yr_Reverse_Ofish(y))<<" ";} else {SIS_table<<" _ ";}
     SIS_table<<Smry_Table(y,5)<<" "<<Smry_Table(y,6)<<endl;
   } // end year loop
 */

  SIS_table<<"Category Year Abundance Abundance Recruitment Spawners Catch Catch Catch Catch Catch Catch Fmort Fmort Fmort Fmort Fmort"<<endl;
  SIS_table<<"Primary _ N Y Y Y N Y N N N N N N N Y Y"<<endl;
  SIS_table<<"Type _ Biomass Biomass Age Female_Mature Sel_Bio Kill_Bio Retain_Bio Sel_Numbers Kill_Numbers Retain_Numbers Exploitation SPR_std F_std Sum_Fleet_Apical_Fs F=Z-M"<<endl;
  SIS_table<<"Source _ Model Model Model Model Model Model Model Model Model Model Model Model Model Model Model"<<endl;
  SIS_table<<"Basis _ Biomass Biomass Numbers Eggs Biomass Biomass Biomass Numbers Numbers Numbers Dead_Catch_Bio/Summary_Bio "<<SPR_report_label<<" "<<F_report_label<<"  Sum_Fleet_Apical_Fs F=Z-M"<<endl;
  SIS_table<<"Range _ All Age_"<<Smry_Age<<"+ Age_0 Mature Exploitable_all Exploitable_dead Exploitable_retained  Exploitable_all Exploitable_dead Exploitable_retained ";
  sprintf(onenum, "%d", int(F_reporting_ages(1)));
  anystring=onenum;
  sprintf(onenum, "%d", int(F_reporting_ages(2)));
  anystring+="_"+onenum;
  SIS_table<<" Exploitable_dead Exploitable_all Exploitable_all Exploitable_all "<<anystring<<endl;
  SIS_table<<"Statistic _ Mean Mean Mean Mean Mean Mean Mean Mean Mean Mean Mean Mean Mean Mean Mean"<<endl;
  SIS_table<<"Scalar _ 1 1 1000 1 1 1 1 1000 1000 1000 1 1 1 1 1"<<endl;

   for (y=styr;y<=YrMax;y++)
   {
     SIS_table<<"_ "<<y<<" "<<Smry_Table(y,1)<<" "<<Smry_Table(y,2)<<" "<<Smry_Table(y,8)<<" "<<Smry_Table(y,7)<<" "<<annual_catch(y)<<" "<<annual_catch(y,2)/Smry_Table(y,2)<<" ";
     if(STD_Yr_Reverse_Ofish(y)>0) {SIS_table<<SPR_std(STD_Yr_Reverse_Ofish(y))<<" ";} else {SIS_table<<" _ ";}
     if(STD_Yr_Reverse_F(y)>0) {SIS_table<<F_std(STD_Yr_Reverse_F(y))<<" ";} else {SIS_table<<" _ ";}
     SIS_table<<annual_F(y)<<endl;
   } // end year loop

// end SPR time series
  SS2out<<endl<<"NOTE:_GENTIME_is_fecundity_weighted_mean_age"<<endl<<"NOTE:_MnAgeSmry_is_numbers_weighted_meanage_at_and_above_smryage(not_accounting_for_birthseason_offsets)"<<endl;

  SS2out<<endl<<"Kobe_Plot"<<endl;
  if(F_report_basis!=2) SS2out<<"F_report_basis_is_not_=2;_so_info_below_is_not_F/Fmsy"<<endl;
  SS2out<<"MSY_basis:_";
  switch(Do_MSY)
    {
    case 1:  // set Fmsy=Fspr
      {SS2out<<"set_Fmsy=Fspr"<<endl;
      break;}
    case 2:  // calc Fmsy
      {SS2out<<"calculate_FMSY"<<endl;
      break;}
    case 3:  // set Fmsy=Fbtgt
      {SS2out<<"set_Fmsy=Fbtgt"<<endl;
      break;}
    case 4:   //  set fmult for Fmsy to 1
      {SS2out<<"set_Fmsy_using_Fmult=1.0"<<endl;
      break;}
    }
  SS2out<<"Year  B/Bmsy  F/Fmsy"<<endl;
   for (y=styr;y<=YrMax;y++)
   {
    SS2out<<y<<" "<<SPB_yr(y)/Bmsy<<" ";
     if(y>=styr && STD_Yr_Reverse_F(y)>0 ) {SS2out<<" "<<F_std(STD_Yr_Reverse_F(y));} else {SS2out<<" _ ";}
     SS2out<<endl;
   }

// ******************************************************************************
  k=Ntypes;
  if(k<4) k=4;
  dvector rmse(1,k);    //  used in the SpBio, Index, Lencomp and Agecomp reports
  dvector Hrmse(1,k);
  dvector Rrmse(1,k);
  dvector n_rmse(1,k);
  dvector mean_CV(1,k);
  dvector mean_CV2(1,k);
//                                                            ProgLabel_330
  rmse = 0.0;  n_rmse = 0.0;

   for (y=recdev_first;y<=recdev_end;y++)
   {
     temp1=recdev(y);
     if(y<recdev_start)  // so in early period
     {
       rmse(3)+=value(square(temp1)); n_rmse(3)+=1.; rmse(4)+=value(biasadj(y));
     }
     else
     {
       rmse(1)+=value(square(temp1)); n_rmse(1)+=1.; rmse(2)+=value(biasadj(y));
     }

   }
   if(n_rmse(1)>0. && rmse(1)>0.) rmse(1) = sqrt(rmse(1)/n_rmse(1));  // rmse during main period
   if(n_rmse(1)>0.) rmse(2) = rmse(2)/n_rmse(1);  // mean biasadj during main period
   if(n_rmse(3)>0. && rmse(3)>0.) rmse(3) = sqrt(rmse(3)/n_rmse(3));  //rmse during early period
   if(n_rmse(3)>0.) rmse(4) = rmse(4)/n_rmse(3);  // mean biasadj during early period

    dvariable Shepard_c;
    dvariable Shepard_c2;
    dvariable Hupper;
    if(SR_fxn==8)
    {
      Shepard_c=SR_parm(3);
      Shepard_c2=pow(0.2,Shepard_c);
      Hupper=1.0/(5.0*Shepard_c2);
      temp=0.2+(SR_parm(2)-0.2)/(0.8)*(Hupper-0.2);
    }
  SS2out<<endl<<"SPAWN_RECRUIT Function: "<<SR_fxn<<" _ _ _ _ _ _"<<endl<<
  SR_parm(1)<<" Ln(R0) "<<mfexp(SR_parm(1))<<endl<<
  SR_parm(2)<<" steep"<<endl<<
  Bmsy/SPB_virgin<<" Bmsy/Bzero ";
  if(SR_fxn==8) SS2out<<Shepard_c<<" Shepard_c "<<Hupper<<" steepness_limit "<<temp<<" Adjusted_steepness";
  SS2out<<endl;
  SS2out<<sigmaR<<" sigmaR"<<endl;
  SS2out<<SR_parm(N_SRparm2-2)<<" env_link_";
  if(SR_env_link>0)
    {
    SS2out<<"to_envvar:_"<<SR_env_link<<"_with_affect_on:";
    if(SR_env_target==1)
      {SS2out<<"_Annual_devs";}
    else if(SR_env_target==2)
      {SS2out<<"_Rzero";}
    else if(SR_env_target==3)
      {SS2out<<"_Steepness";}
    }
  SS2out<<endl<<SR_parm(N_SRparm2-1)<<" init_eq "<<mfexp(SR_parm(1)+SR_parm(N_SRparm2-1))<<endl<<
  recdev_start<<" "<<recdev_end<<" main_recdev:start_end"<<endl<<
  recdev_adj(1)<<" "<<recdev_adj(2,5)<<" breakpoints_for_bias_adjustment_ramp "<<endl;

   temp=sigmaR*sigmaR;  //  sigmaR^2
   SS2out<<"ERA    N    RMSE  RMSE^2/sigmaR^2  mean_BiasAdj"<<endl;
   SS2out<<"main  "<<n_rmse(1)<<" "<<rmse(1)<<" "<<square(rmse(1))/temp<<" "<<rmse(2);
   if(rmse(1)<0.5*sigmaR && rmse(2)>(0.01+2.0*square(rmse(1))/temp))
   {N_warn++; warning<<" Main recdev biasadj is >2 times ratio of rmse to sigmaR"<<endl; SS2out<<" # Main_recdev_biasadj_is_>2_times_ratio_of_rmse_to_sigmaR"<<endl;}
   else
   {SS2out<<endl;}
   SS2out<<"early "<<n_rmse(3)<<" "<<rmse(3)<<" "<<square(rmse(3))/temp<<" "<<rmse(4);
   if(rmse(3)<0.5*sigmaR && rmse(4)>(0.01+2.0*square(rmse(3))/temp))
   {N_warn++; warning<<" Early recdev biasadj is >2 times ratio of rmse to sigmaR"<<endl; SS2out<<" # Early_recdev_biasadj_is_>2_times_ratio_of_rmse_to_sigmaR"<<endl;}
   else
   {SS2out<<endl;}

  SS2out<<"year spawn_bio exp_recr with_env adjusted pred_recr dev biasadj era mature_bio mature_num"<<endl;
  SS2out<<"S/Rcurve "<<SPB_virgin<<" "<<Recr_virgin<<endl;
  y=styr-2;
  SS2out<<"Virg "<<SPB_yr(y)<<" "<<exp_rec(y)<<" - "<<0.0<<" Virg "<<SPB_B_yr(y)<<" "<<SPB_N_yr(y)<<endl;
  y=styr-1;
  SS2out<<"Init "<<SPB_yr(y)<<" "<<exp_rec(y)<<" - "<<0.0<<" Init "<<SPB_B_yr(y)<<" "<<SPB_N_yr(y)<<endl;

  if(recdev_first<styr)
  {
    for (y=recdev_first;y<=styr-1;y++)
    {
      SS2out<<y<<" "<<SPB_yr(styr-1)<<" "<<exp_rec(styr-1,1)<<" "<<exp_rec(styr-1,2)<<" "<<exp_rec(styr-1,3)*mfexp(-biasadj(y)*half_sigmaRsq)<<" "<<
      exp_rec(styr-1,3)*mfexp(recdev(y)-biasadj(y)*half_sigmaRsq)<<" "
      <<recdev(y)<<" "<<biasadj(y)<<" Init_age "<<SPB_B_yr(styr-1)<<" "<<SPB_N_yr(styr-1)<<endl;
    }
  }
   for (y=styr;y<=YrMax;y++)
   {
     SS2out<<y<<" "<<SPB_yr(y)<<" "<<exp_rec(y)<<" ";
     if(recdev_do_early>0 && y>=recdev_early_start && y<=recdev_early_end)
       {SS2out<<recdev(y)<<" "<<biasadj(y)<<" Early "<<SPB_B_yr(y)<<" "<<SPB_N_yr(y);}
     else if(y>=recdev_start && y<=recdev_end)
       {SS2out<<recdev(y)<<" "<<biasadj(y)<<" Main "<<SPB_B_yr(y)<<" "<<SPB_N_yr(y);}
     else if(Do_Forecast>0 && y>recdev_end)
     {
        SS2out<<Fcast_recruitments(y)<<" "<<biasadj(y);
        if(y<=endyr)
        {SS2out<<" Late "<<SPB_B_yr(y)<<" "<<SPB_N_yr(y);}
        else
        {SS2out<<" Forecast "<<SPB_B_yr(y)<<" "<<SPB_N_yr(y);}
      }
     else
       {SS2out<<" _ _ Fixed";}
//       SS2out<<" "<<recdev_cycle_parm(gg);
     SS2out<<endl;
   }

// ******************************************************************************
//                                             ProgLabel_340

  SS2out <<endl<< "INDEX_2" << endl;
  rmse = 0.0;  n_rmse = 0.0; mean_CV=0.0; mean_CV2=0.0;
  SS2out<<"Fleet Name Yr Seas Yr.S Vuln_bio Obs Exp Calc_Q Eff_Q SE Dev Like Like+log(s) Supr_Per Use"<<endl;
  if(nobs_cr>0)
  {
    for(f=1;f<=Ntypes;f++)
    {
      in_superperiod=0;
      for (i=1;i<=nyr_cr(f);i++)
      {
          y=abs(yr_cr_y(f,i));
          s=abs(yr_cr_s(f,i));
          temp = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
          SS2out<<f<<" "<<fleetname(f)<<" "<<yr_cr_y(f,i)<<" "<<abs(yr_cr_s(f,i))<<" "<<temp<<" "<<vul_bio(f,i)<<" "<<obs_cr(f,i)<<" ";
          if(cr_errtype(f)>=0)  // lognormal
          {
            temp = exp_cr(f,i)+log_q_cr(f,i);
            SS2out<<mfexp(temp)<<" "<<mfexp(log_q_cr(f,i))<<" "<<mfexp(temp)/vul_bio(f,i)<<" "<<se_cr_use(f,i);
            if(yr_cr_use(f,i) > 0)
            {
              SS2out<<" "<<Ln_obs_cr(f,i)-temp<<" ";
              if(cr_errtype(f)==0)
              {
                SS2out<<0.5*square( ( Ln_obs_cr(f,i)-temp ) / se_cr_use(f,i))<<" "
                <<0.5*square( ( Ln_obs_cr(f,i)-temp ) / se_cr_use(f,i))+log(se_cr_use(f,i));
              }
              else if (cr_errtype(f)>0)  // student's T
              {
                SS2out<<((cr_errtype(f)+1.)/2.)*log((1.+square((Ln_obs_cr(f,i)-temp ))/(cr_errtype(f)*square(se_cr_use(f,i))) ))<<" "
                <<((cr_errtype(f)+1.)/2.)*log((1.+square((Ln_obs_cr(f,i)-temp ))/(cr_errtype(f)*square(se_cr_use(f,i))) ))+log(se_cr_use(f,i));
              }
              rmse(f)+=value(square(Ln_obs_cr(f,i)-temp)); n_rmse(f)+=1.;
              mean_CV(f)+=se_cr_obs(f,i); mean_CV2(f)+=value(se_cr_use(f,i));
            }
            else
            {
              SS2out<<" _ _ _ ";
            }
          }
          else  // normal
          {
            temp = exp_cr(f,i)*q_cr(f,i);
            SS2out<<temp<<" "<<mfexp(log_q_cr(f,i))<<" "<<temp/vul_bio(f,i)<<" "<<se_cr_use(f,i);
            if(yr_cr_use(f,i)>0)
            {
              SS2out<<" "<<Ln_obs_cr(f,i)-temp<<" ";
              SS2out<<0.5*square( ( Ln_obs_cr(f,i)-temp ) / se_cr_use(f,i))<<" "
              <<0.5*square( ( Ln_obs_cr(f,i)-temp ) / se_cr_use(f,i))+log(se_cr_use(f,i));
              rmse(f)+=value(square(Ln_obs_cr(f,i)-temp)); n_rmse(f)+=1.;
              mean_CV(f)+=se_cr_obs(f,i); mean_CV2(f)+=value(se_cr_use(f,i));
            }
            else
            {
              SS2out<<" _ _ _ ";
            }
          }
          if(yr_cr_s(f,i)<0 &&in_superperiod==0)
          {in_superperiod=1; SS2out<<" beg_supr_per ";}
          else if(yr_cr_s(f,i)<0 &&in_superperiod==1)
          {in_superperiod=0; SS2out<<" end_supr_per ";}
          else if(in_superperiod==1)
          {SS2out<<" in_suprper ";}
          else{SS2out<<" _ ";}
          SS2out<<yr_cr_use(f,i);
          SS2out<<endl;
      }
      if(n_rmse(f)>0) {rmse(f) = sqrt((rmse(f)+1.0e-9)/n_rmse(f)); mean_CV(f) /= n_rmse(f); mean_CV2(f) /= n_rmse(f);}
    }
  }

  SS2out <<endl<< "INDEX_1" << endl;
  SS2out <<"Fleet Do_Power Power Do_Env_var Env_Link Do_ExtraVar Qtype  Q Num=0/Bio=1 Err_type"<<
    " N Npos r.m.s.e. mean_input_SE Input+VarAdj Input+VarAdj+extra VarAdj New_VarAdj pen_mean_Qdev rmse_Qdev"<<endl;
  for(f=1;f<=Ntypes;f++)
    {
    SS2out<<f<<" "<<Q_setup(f,1)<<" ";
    if(Q_setup(f,1)>0)
      {SS2out<<Q_parm(Q_setup(f,1))<<" ";}
    else
      {SS2out<<" 1.0 ";}
    SS2out<<Q_setup(f,2)<<" ";
    if(Q_setup(f,2)!=0)
      {SS2out<<Q_parm(Q_setup_parms(f,2));}
    else
      {SS2out<<" 0.0";}
    SS2out<<Q_setup(f,3)<<" ";
    if(Q_setup(f,3)>0)
      {SS2out<<Q_parm(Q_setup_parms(f,3));}
    else
      {SS2out<<" 0.0";}
    SS2out<<" "<<Q_setup(f,4)<<" ";
    if(nyr_cr(f)>0)
      {SS2out<<q_cr(f,1);}
    else
      {SS2out<<" _";}
    SS2out<<" "<<cr_units(f)<<" "<<cr_errtype(f)<<" "<<nyr_cr(f)<<" "<<n_rmse(f)<<" "<<rmse(f)<<" "
      <<mean_CV(f)-var_adjust(1,f)<<" "<<mean_CV(f)<<" "<<mean_CV2(f)<<" "<<var_adjust(1,f)
      <<" "<<var_adjust(1,f)+rmse(f)-mean_CV(f)
      <<" "<<Q_dev_like(f,1)<<" "<<Q_dev_like(f,2)<<" "<<fleetname(f)<<endl;
    }
    SS2out<<"rmse_Qdev_not_in_logL"<<endl<<"pen_mean_Qdev_not_in_logL_in_randwalk_approach"<<endl;

  SS2out <<"#"<<endl<< "INDEX_3"<<endl<<"Fleet  Q_parm_assignments"<<endl;
  for(f=1;f<=Ntypes;f++)
    {SS2out<<f<<" "<<Q_setup_parms(f,1)<<" _ "<<Q_setup_parms(f,2)<<" _ "<<Q_setup_parms(f)(3,4)<<" "<<fleetname(f)<<endl;}

  SS2out<<endl<<"DISCARD_SPECIFICATION ";
  SS2out<<"Discard_units_options"<<endl;
   SS2out << "1:  discard_in_biomass(mt)_or_numbers(1000s)_to_match_catchunits_of_fleet"<<endl;
   SS2out << "2:  discard_as_fraction_of_total_catch(based_on_bio_or_num_depending_on_fleet_catchunits)"<<endl;
   SS2out << "3:  discard_as_numbers(1000s)_regardless_of_fleet_catchunits"<<endl;
  SS2out<<"Discard_errtype_options"<<endl;
  SS2out << ">1:  log(L)_based_on_T_distribution_with_specified_DF"<< endl;
  SS2out << "0:  log(L)_based_on_normal_with_Std_in_as_CV"<< endl;
  SS2out << "-1:  log(L)_based_on_normal_with_Std_in_as_stddev"<< endl;
  SS2out << "-2:  log(L)_based_on_lognormal_with_Std_in_as_stddev_in_logspace"<< endl;

  SS2out<<"#_Fleet units errtype"<<endl;
  if(Ndisc_fleets>0)
  {
    for(f=1;f<=Ntypes;f++)
    if(disc_units(f)>0) SS2out<<f<<" "<<disc_units(f)<<" "<<disc_errtype(f)<<" # "<<fleetname(f)<<endl;
  }

  SS2out<<"#"<<endl<<"DISCARD_OUTPUT "<<endl;
  SS2out<<"Fleet Name Yr Seas Yr.S Obs Exp Std_in Std_use Dev Like Like+log(s) SuprPer Use Obs_cat Exp_cat F_rate"<<endl;
  if(nobs_disc>0)
  for(f=1;f<=Nfleet;f++)
  for (y=styr;y<=endyr;y++)
  for (s=1;s<=nseas;s++)
  {
    t = styr+(y-styr)*nseas+s-1;
        if(catchunits(f)==1)
        {gg=3;}  //  biomass
        else
        {gg=6;}  //  numbers
    if(yr_disc2(f,t)>0)
    {
      i=yr_disc2(f,t);
      temp = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
      SS2out<<f<<" "<<fleetname(f)<<" "<<y<<" "<<s<<" "<<temp<<" "<<obs_disc(f,i)<<" "
      <<exp_disc(f,i)<<" "<<" "<<cv_disc(f,i)<<" "<<sd_disc(f,i);

      if(yr_disc_use(f,i)>=0.0)
      {
        if(disc_errtype(f)>=1)  // T -distribution
        {
          temp=0.5*(disc_errtype(f)+1.)*log((1.+square(obs_disc(f,i)-exp_disc(f,i))/(disc_errtype(f)*square(sd_disc(f,i))) ));
          SS2out<<" "<<obs_disc(f,i)-exp_disc(f,i)<<" "<<temp<<" "<<temp + sd_offset*log(sd_disc(f,i));
        }
        else if (disc_errtype(f)==0)  // normal error, with input CV
        {
          temp=0.5*square( (obs_disc(f,i)-exp_disc(f,i)) / sd_disc(f,i));
          SS2out<<" "<<obs_disc(f,i)-exp_disc(f,i)<<" "<<temp<<" "<<temp + sd_offset*log(sd_disc(f,i));
        }
        else if (disc_errtype(f)==-1)  // normal error with input se
        {
          temp=0.5*square( (obs_disc(f,i)-exp_disc(f,i)) / sd_disc(f,i));
          SS2out<<" "<<obs_disc(f,i)-exp_disc(f,i)<<" "<<temp<<" "<<temp + sd_offset*log(sd_disc(f,i));
        }
        else  // lognormal  where input cv_disc must contain se in log space
        {
          temp=0.5*square( log(obs_disc(f,i)/exp_disc(f,i)) / sd_disc(f,i));
          SS2out<<" "<<log(obs_disc(f,i)/exp_disc(f,i))<<" "<<temp<<" "<<temp + sd_offset*log(sd_disc(f,i));
        }
      }
      else
      {
        SS2out<<"  _  _  _  ";
      }
      if(yr_disc_s(f,i)<0 &&in_superperiod==0)
      {in_superperiod=1; SS2out<<" beg_suprper ";}
      else if(yr_disc_s(f,i)<0 &&in_superperiod==1)
      {in_superperiod=0; SS2out<<" end_suprper ";}
      else if(in_superperiod==1)
      {SS2out<<" in_suprper ";}
      else{SS2out<<" _ ";}
      SS2out<<yr_disc_use(f,i);
      SS2out<<" "<<catch_ret_obs(f,t)<<" "<<catch_fleet(t,f,gg)<<" "<<Hrate(f,t);
      SS2out<<endl;
    }
  }

  SS2out <<endl<< "MEAN_BODY_WT_OUTPUT"<<endl;
  if(nobs_mnwt>0) SS2out<<"log(L)_based_on_T_distribution_with_DF=_"<<DF_bodywt<< endl;
  SS2out<<"Fleet Name Yr Seas Yr.S Mkt Obs Exp CV Dev NeglogL Neg(logL+log(s)) Use"<<endl;
  if(nobs_mnwt>0)
  for(i=1;i<=nobs_mnwt;i++)
  {
    y=abs(mnwtdata(1,i));
    s=abs(mnwtdata(2,i));
    f=abs(mnwtdata(3,i));
    temp = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
    SS2out << mnwtdata(3,i)<<" "<<fleetname(f)<<" "<<y<<" "<<s<<" "<<temp<<" "
    <<mnwtdata(4,i)<<" "<<mnwtdata(5,i)<<" "<<exp_mnwt(i)<<" "<<mnwtdata(6,i);
    if(mnwtdata(3,i)>0.)
    {
      SS2out<<" "<<mnwtdata(5,i)-exp_mnwt(i)<<" "<<
       0.5*(DF_bodywt+1.)*log(1.+square(mnwtdata(5,i)-exp_mnwt(i))/mnwtdata(8,i))<<" "<<
       0.5*(DF_bodywt+1.)*log(1.+square(mnwtdata(5,i)-exp_mnwt(i))/mnwtdata(8,i))+ mnwtdata(9,i)<<" "<<1;
    }
    else
    {
      SS2out<<" NA NA NA -1";
    }
    SS2out<<endl;
  }

  SS2out <<endl<< "FIT_LEN_COMPS" << endl;                     // ProgLabel_350
  SS2out<<"Fleet Yr Seas Gender Mkt Nsamp effN Like Period Use"<<endl;
  rmse = 0.0;  n_rmse = 0.0; mean_CV=0.0;  Hrmse=0.0; Rrmse=0.0; neff_l.initialize();
  in_superperiod=0;
   for(f=1;f<=Ntypes;f++)
   for(i=1;i<=Nobs_l(f);i++)
   {
     temp=0.00;
     neff_l(f,i)  = exp_l(f,i)*(1-exp_l(f,i))+1.0e-06;     // constant added for stability
     neff_l(f,i) /= (obs_l(f,i)-exp_l(f,i))*(obs_l(f,i)-exp_l(f,i))+1.0e-06;
     if(gen_l(f,i) !=2)
     {
       temp -= obs_l(f,i)(tails_l(f,i,1),tails_l(f,i,2)) * log(exp_l(f,i)(tails_l(f,i,1),tails_l(f,i,2)));
       temp += obs_l(f,i)(tails_l(f,i,1),tails_l(f,i,2)) * log(obs_l(f,i)(tails_l(f,i,1),tails_l(f,i,2)));
     }
     if(gen_l(f,i) >=2 && gender==2)
     {
       temp -= obs_l(f,i)(tails_l(f,i,3),tails_l(f,i,4)) * log(exp_l(f,i)(tails_l(f,i,3),tails_l(f,i,4)));
       temp += obs_l(f,i)(tails_l(f,i,3),tails_l(f,i,4)) * log(obs_l(f,i)(tails_l(f,i,3),tails_l(f,i,4)));
     }
     if(nsamp_l(f,i)>0 && header_l(f,i,3)>0)
     {
       n_rmse(f)+=1.;
       rmse(f)+=value(neff_l(f,i));
       mean_CV(f)+=nsamp_l(f,i);
       Hrmse(f)+=value(1./neff_l(f,i));
       Rrmse(f)+=value(neff_l(f,i)/nsamp_l(f,i));
     }

     SS2out<<f<<" "<<yr_l_y(f,i)<<" "<<abs(header_l(f,i,2))<<" "<<gen_l(f,i)<<" "<<mkt_l(f,i)<<" "<<nsamp_l(f,i)<<" "<<neff_l(f,i)<<" "<<
     temp*sfabs(nsamp_l(f,i))<<" ";
     if(header_l(f,i,2)<0 && in_superperiod==0)
      {SS2out<<" start "; in_superperiod=1;}
      else if (header_l(f,i,2)<0 && in_superperiod==1)
      {SS2out<<" end "; in_superperiod=0;}
      else if (in_superperiod==1)
      {SS2out<<" in ";}
      else
      {SS2out<<" _ ";}
      if(header_l(f,i,3)<0 || nsamp_l(f,i)<0)
      {SS2out<<" skip "<<endl;}
      else
      {SS2out<<" _ "<<endl;}
    }

   SS2out<<endl<<"Fleet N Npos mean_effN mean(inputN*Adj) HarMean(effN) Mean(effN/inputN) MeaneffN/MeaninputN Var_Adj"<<endl;
   for(f=1;f<=Ntypes;f++)
   {
    if(n_rmse(f)>0) {rmse(f)/=n_rmse(f); mean_CV(f)/=n_rmse(f); Hrmse(f)=n_rmse(f)/Hrmse(f); Rrmse(f)/=n_rmse(f); }
    SS2out<<f;
    if(Nobs_l(f)>0)
    {SS2out<<" "<<Nobs_l(f)<<" "<<n_rmse(f)<<" "<<rmse(f)<<" "<<mean_CV(f)<<" "<<Hrmse(f)<<" "<<Rrmse(f)<<" "<<rmse(f)/mean_CV(f)
    <<" "<<var_adjust(4,f)<<" "<<fleetname(f)<<endl;}
    else
    {SS2out<<" _ _ _ _ _ _ _ _ "<<endl;}
   }

  SS2out <<endl<< "FIT_AGE_COMPS" << endl;                        // ProgLabel_360
  SS2out<<"Fleet Yr Seas Gender Mkt Ageerr Lbin_lo Lbin_hi Nsamp effN Like 5% Mean 95% Super Use"<<endl;
  rmse = 0.0;  n_rmse = 0.0; mean_CV=0.0;  Hrmse=0.0; Rrmse=0.0;
   if(Nobs_a_tot>0)
   for(f=1;f<=Ntypes;f++)
   for(i=1;i<=Nobs_a(f);i++)
     {
      neff_a(f,i)  = exp_a(f,i)*(1-exp_a(f,i))+1.0e-06;     // constant added for stability
      neff_a(f,i) /= (obs_a(f,i)-exp_a(f,i))*(obs_a(f,i)-exp_a(f,i))+1.0e-06;
      temp=0.00;
     if(gen_a(f,i) !=2)
      {temp -= obs_a(f,i)(tails_a(f,i,1),tails_a(f,i,2)) * log(exp_a(f,i)(tails_a(f,i,1),tails_a(f,i,2)));
       temp += obs_a(f,i)(tails_a(f,i,1),tails_a(f,i,2)) * log(obs_a(f,i)(tails_a(f,i,1),tails_a(f,i,2)));
      }
     if(gen_a(f,i) >=2 && gender==2)
      {temp -= obs_a(f,i)(tails_a(f,i,3),tails_a(f,i,4)) * log(exp_a(f,i)(tails_a(f,i,3),tails_a(f,i,4)));
       temp += obs_a(f,i)(tails_a(f,i,3),tails_a(f,i,4)) * log(obs_a(f,i)(tails_a(f,i,3),tails_a(f,i,4)));
      }
     if(nsamp_a(f,i)>0 && header_a(f,i,3)>0)
       {
        n_rmse(f)+=1.;
        rmse(f)+=value(neff_a(f,i));
        mean_CV(f)+=nsamp_a(f,i);
        Hrmse(f)+=value(1./neff_a(f,i));
        Rrmse(f)+=value(neff_a(f,i)/nsamp_a(f,i));
       }
      SS2out<<f<<" "<<yr_a_y(f,i)<<" "<<abs(header_a(f,i,2))<<" "<<gen_a(f,i)<<" "<<mkt_a(f,i)<<" "<<ageerr_type_a(f,i)<<" "<<Lbin_lo(f,i)<<" "<<Lbin_hi(f,i)<<" "<<nsamp_a(f,i)<<" "<<neff_a(f,i)<<" "<<
      temp*sfabs(nsamp_a(f,i));
      SS2out<<exp_meanage(f,i)<<" ";
     if(header_a(f,i,2)<0 && in_superperiod==0)
      {SS2out<<" start "; in_superperiod=1;}
      else if (header_a(f,i,2)<0 && in_superperiod==1)
      {SS2out<<" end "; in_superperiod=0;}
      else if (in_superperiod==1)
      {SS2out<<" in ";}
      else
      {SS2out<<" _ ";}
      if(header_a(f,i,3)<0 || nsamp_a(f,i)<0)
      {SS2out<<" skip "<<endl;}
      else
      {SS2out<<" _ "<<endl;}
      }
   SS2out<<endl<<"Fleet N Npos mean_effN mean(inputN*Adj) HarMean(effN) Mean(effN/inputN) MeaneffN/MeaninputN Var_Adj"<<endl;
   for(f=1;f<=Ntypes;f++)
   {
    if(n_rmse(f)>0) {rmse(f)/=n_rmse(f); mean_CV(f)/=n_rmse(f); Hrmse(f)=n_rmse(f)/Hrmse(f); Rrmse(f)/=n_rmse(f); }
    SS2out<<f;
    if(Nobs_a(f)>0)
    {SS2out<<" "<<Nobs_a(f)<<" "<<n_rmse(f)<<" "<<rmse(f)<<" "<<mean_CV(f)<<" "<<Hrmse(f)<<" "<<Rrmse(f)<<" "<<rmse(f)/mean_CV(f)
    <<" "<<var_adjust(5,f)<<" "<<fleetname(f)<<endl;}
    else
    {SS2out<<" _ _ _ _ _ _ _ _ "<<endl;}
   }

  SS2out <<endl<< "FIT_SIZE_COMPS" << endl;                     // ProgLabel_350
  rmse = 0.0;  n_rmse = 0.0; mean_CV=0.0;  Hrmse=0.0; Rrmse=0.0;
    if(SzFreq_Nmeth>0)       //  have some sizefreq data
    {
      SzFreq_effN.initialize();
      SzFreq_eachlike.initialize();
      SS2out<<"Fleet Yr Seas Method Gender Mkt Nsamp effN Like"<<endl;
      for(iobs=1;iobs<=SzFreq_totobs;iobs++)
      {
        y=SzFreq_obs_hdr(iobs,1);
        s=SzFreq_obs_hdr(iobs,2);
        f=abs(SzFreq_obs_hdr(iobs,3));
        gg=SzFreq_obs_hdr(iobs,4);  // gender
        k=SzFreq_obs_hdr(iobs,6);
        if(SzFreq_obs_hdr(iobs,3)>0)  // flag for date range in bounds
        {
          p=SzFreq_obs_hdr(iobs,5);  // partition
          z1=SzFreq_obs_hdr(iobs,7);
          z2=SzFreq_obs_hdr(iobs,8);
          temp=0.0;
          temp1=0.0;
          for(z=z1;z<=z2;z++)
          {
            SzFreq_effN(iobs)+= value(SzFreq_exp(iobs,z)*(1.0-SzFreq_exp(iobs,z)));
            temp += square(SzFreq_obs(iobs,z)-SzFreq_exp(iobs,z));
            temp1 += SzFreq_obs(iobs,z)*log(SzFreq_obs(iobs,z))-SzFreq_obs(iobs,z)*log(SzFreq_exp(iobs,z));
          }
          SzFreq_effN(iobs) =(SzFreq_effN(iobs)+1.0e-06)/value((temp+1.0e-06));
          temp1*=SzFreq_sampleN(iobs);
          SzFreq_eachlike(iobs)=value(temp1);
          if(SzFreq_obs_hdr(iobs,3)>0)
            {
            n_rmse(f)+=1.;
            rmse(f)+=SzFreq_effN(iobs);
            mean_CV(f)+=SzFreq_sampleN(iobs);
            Hrmse(f)+=1./SzFreq_effN(iobs);
            Rrmse(f)+=SzFreq_effN(iobs)/SzFreq_sampleN(iobs);
          }
          SS2out<<SzFreq_obs_hdr(iobs,3)<<" "<<y<<" "<<s<<" "<<k<<" "<<gg<<" "<<p<<" "<<SzFreq_sampleN(iobs)<<" "<<SzFreq_effN(iobs)<<" "<<SzFreq_eachlike(iobs)<<endl;
        }
      }
      SS2out<<endl<<"Fleet N Npos mean_effN mean(inputN*Adj) HarMean(effN) Mean(effN/inputN) MeaneffN/MeaninputN Var_Adj"<<endl;
      for(f=1;f<=Ntypes;f++)
      {
        if(n_rmse(f)>0)
        {
          rmse(f)/=n_rmse(f); mean_CV(f)/=n_rmse(f); Hrmse(f)=n_rmse(f)/Hrmse(f); Rrmse(f)/=n_rmse(f);
          SS2out<<f<<" "<<"NA"<<" "<<n_rmse(f)<<" "<<rmse(f)<<" "<<mean_CV(f)<<" "<<Hrmse(f)<<" "<<Rrmse(f)<<" "<<rmse(f)/mean_CV(f)
          <<" "<<var_adjust(4,f)<<" "<<fleetname(f)<<endl;
        }
      }
    }
    else
    {SS2out<<"#_none"<<endl;}

  SS2out <<"# "<<endl<<"OVERALL_COMPS"<<endl;
  SS2out<<"Fleet N_obs len_bins "<<len_bins_dat<<endl;
  for(f=1;f<=Ntypes;f++)
  {
    if(Nobs_l(f)>0)
    {
      SS2out<<f<<" "<<Nobs_l(f)<<" freq "<<obs_l_all(1,f)<<endl;
      SS2out<<f<<" "<<Nobs_l(f)<<" cum  "<<obs_l_all(2,f)<<endl;
    }
  }

  SS2out<<"Fleet N_obs age_bins ";
  if(n_abins>1)
  {
    SS2out<<age_bins(1,n_abins)<<endl;
    for(f=1;f<=Ntypes;f++)
    {
      if(Nobs_a(f)>0)
      {
        SS2out<<f<<" "<<Nobs_a(f)<<" freq "<<obs_a_all(1,f)<<endl;
        SS2out<<f<<" "<<Nobs_a(f)<<" cum  "<<obs_a_all(2,f)<<endl;
      }
    }
  }
  else
  {SS2out<<"No_age_bins_defined"<<endl;}

  SS2out <<"# "<<endl<<"LEN_SELEX"<<endl;
  SS2out << "Lsel_is_length_selectivity" << endl;     // ProgLabel_370
  SS2out << "RET_is_retention" << endl;            // ProgLabel_390
  SS2out << "MORT_is_discard_mortality" << endl;            // ProgLabel_390
  SS2out << "KEEP_is_sel*retain" << endl;     // ProgLabel_370
  SS2out << "DEAD_is_sel*(retain+(1-retain)*discmort)";     // ProgLabel_370
  SS2out<<"; Year_styr-3_("<<styr-3<<")_stores_average_used_for_benchmark"<<endl; 
  SS2out<<"Factor Fleet year gender label "<<len_bins_m<<endl;
  for (f=1;f<=Ntypes;f++)
  {
    if(f<=Nfleet) {k=styr-3; j=endyr+1;} else {k=styr; j=endyr;}
    for (y=k;y<=j;y++)
    for (gg=1;gg<=gender;gg++)
    if(y==styr-3 || y==endyr || (y>=styr && (time_vary_sel(y,f)>0 || time_vary_sel(y+1,f)>0)))
    {
      SS2out<<"Lsel "<<f<<" "<<y<<" "<<gg<<" "<<y<<"_"<<f<<"_Lsel";
      for(z=1;z<=nlength;z++) {SS2out<<" "<<sel_l(y,f,gg,z);}
      SS2out<<endl;
    }
  }

  for (f=1;f<=Nfleet;f++)
  for (y=styr-3;y<=endyr+1;y++)
  for (gg=1;gg<=gender;gg++)
  if(y==styr-3 || y==endyr || (y>=styr && (time_vary_sel(y,f)>0 || time_vary_sel(y+1,f)>0)))
  {
    if(y>=styr && y<=endyr)
    {
      SS2out<<"Ret "<<f<<" "<<y<<" "<<gg<<" "<<y<<"_"<<f<<"_Ret";
      if(gg==1) {for(z=1;z<=nlength;z++) {SS2out<<" "<<retain(y,f,z);}}
      else
      {for(z=nlength1;z<=nlength2;z++) {SS2out<<" "<<retain(y,f,z);}}
      SS2out<<endl;
      SS2out<<"Mort "<<f<<" "<<y<<" "<<gg<<" "<<y<<"_"<<f<<"_Mort";
      if(gg==1) {for(z=1;z<=nlength;z++) {SS2out<<" "<<discmort(y,f,z);}}
      else
      {for(z=nlength1;z<=nlength2;z++) {SS2out<<" "<<discmort(y,f,z);}}
      SS2out<<endl;
    }
    SS2out<<"Keep "<<f<<" "<<y<<" "<<gg<<" "<<y<<"_"<<f<<"_Keep";
    for(z=1;z<=nlength;z++) {SS2out<<" "<<sel_l_r(y,f,gg,z);}
    SS2out<<endl;
    SS2out<<"Dead "<<f<<" "<<y<<" "<<gg<<" "<<y<<"_"<<f<<"_Dead";
    for(z=1;z<=nlength;z++) {SS2out<<" "<<discmort2(y,f,gg,z);}
    SS2out<<endl;
  }

  SS2out <<endl<< "AGE_SELEX" << endl;
  SS2out<<"Asel_is_age_selectivity_alone"<<endl;
  SS2out<<"Asel2_is_sizesel*size_at_age(ALK)"<<endl;
  SS2out<<"COMBINED_ALK*selL*selA*wtlen*ret*discmort_in_makefishsel_yr: "<<makefishsel_yr<<" With_MeanSel_From: "<<Fcast_Sel_yr1<<" - "<<Fcast_Sel_yr2;     // ProgLabel_380
  SS2out<<"; Year_styr-3_("<<styr-3<<")_stores_average_used_for_benchmark"<<endl; 

  SS2out<<"factor fleet year seas gender morph label ";
  for(a=0;a<=nages;a++) {SS2out<<" "<<a;}
  SS2out<<endl;
  for (f=1;f<=Ntypes;f++)
  {
    if(f<=Nfleet) {k=styr-3; j=endyr+1;} else {k=styr; j=endyr;}
    for (y=k;y<=j;y++)
    for(gg=1;gg<=gender;gg++)
    if(y==styr-3 || y==endyr || (y>=styr && (time_vary_sel(y,f+Ntypes)>0 || time_vary_sel(y+1,f+Ntypes)>0)))
    {
      SS2out<<"Asel "<<f<<" "<<y<<" 1 "<<gg<<" 1 "<<y<<"_"<<f<<"Asel";
      for(a=0;a<=nages;a++) {SS2out<<" "<<sel_a(y,f,gg,a);}
      SS2out<<endl;
    }
  }

  if(reportdetail>0)
  {
    if(Do_Forecast>0)
    {k=endyr+1;}
    else
    {k=endyr;}
    for (y=styr-3;y<=k;y++)
    for (s=1;s<=nseas;s++)
    {
      t=styr+(y-styr)*nseas+s-1;
      for(g=1;g<=gmorph;g++)
      if(use_morph(g)>0 && (y==styr-3 || y>=styr))
      {
        if(s==spawn_seas && (sx(g)==1 || Hermaphro_Option>0) ) SS2out<<"Fecund "<<" NA "<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<"Fecund"<<save_sel_fec(t,g,0)<<endl;
        for(f=1;f<=Nfleet;f++)
        {
          SS2out<<"Asel2 "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_Asel2"<<save_sel_fec(t,g,f)<<endl;
          SS2out<<"bodywt "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_bodywt"<<fish_body_wt(t,g,f)<<endl;
        }
      }
    }
      y=makefishsel_yr;
      for(f=1;f<=Nfleet;f++)
       for(g=1;g<=gmorph;g++)
       if(use_morph(g)>0)
       for(s=1;s<=nseas;s++)
        {
        SS2out<<"sel*wt "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_sel*wt"<<sel_al_1(s,g,f)<<endl;
        SS2out<<"sel*ret*wt "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_sel*ret*wt"<<sel_al_2(s,g,f)<<endl;
        SS2out<<"sel_nums "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_sel_nums"<<sel_al_3(s,g,f)<<endl;
        SS2out<<"sel*ret_nums "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_sel*ret_nums"<<sel_al_4(s,g,f)<<endl;
        SS2out<<"dead_nums "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_dead_nums"<<deadfish(s,g,f)<<endl;
        SS2out<<"dead*wt "<<f<<" "<<y<<" "<<s<<" "<<sx(g)<<" "<<g<<" "<<y<<"_"<<f<<"_dead*wt"<<deadfish_B(s,g,f)<<endl;
        }
    }

   SS2out << endl<<"ENVIRONMENTAL_DATA Begins_in_startyr-1"<<endl;         // ProgLabel_397
   if(N_envdata>=1)
   {
   SS2out<<"Year "; for(i=-2;i<=N_envvar;i++) SS2out<<" "<<i;
   SS2out<<endl;
    for (y=styr-1;y<=YrMax;y++)
    {
     SS2out<<y<<" "<<env_data(y)<<endl;
    }
    SS2out<<endl;
   }

  SS2out<<endl<<"TAG_Recapture"<<endl;
  SS2out<<TG_mixperiod<<" First period to use recaptures in likelihood"<<endl;
  SS2out<<TG_maxperiods<<" Accumulation period"<<endl;
  if(Do_TG>0)
  {
    SS2out<<" Tag_release_info"<<endl;
    SS2out<<"TAG Area Yr Seas Time Gender Age Nrelease Init_Loss Chron_Loss"<<endl;;
    for(TG=1;TG<=N_TG;TG++)
    {
      SS2out<<TG<<" "<<TG_release(TG)(2,8)<<" "<<TG_save(TG)(1,2)<<endl;
    }
    SS2out<<"Tags_Alive ";
    k=max(TG_endtime);
    for(t=0;t<=k;t++) SS2out<<t<<" ";
    SS2out<<endl;
    for(TG=1;TG<=N_TG;TG++)
    {
      SS2out<<TG<<" "<<TG_save(TG)(3,3+TG_endtime(TG))<<endl;
    }
    SS2out<<"Total_recaptures ";
    for(t=0;t<=k;t++) SS2out<<t<<" ";
    SS2out<<endl;
    for(TG=1;TG<=N_TG;TG++)
    {
      SS2out<<TG<<" ";
      for(TG_t=0;TG_t<=TG_endtime(TG);TG_t++) SS2out<<TG_recap_exp(TG,TG_t,0)<<" ";
      SS2out<<endl;
    }

    SS2out<<endl<<"Reporting_Rates_by_Fishery"<<endl<<"Fleet Init_Reporting Report_Decay"<<endl;
    for(f=1;f<=Nfleet;f++) SS2out<<f<<" "<<TG_report(f)<<" "<<TG_rep_decay(f)<<endl;
    SS2out<<"See_composition_data_output_for_tag_recapture_details"<<endl;
  }


// ************************                     ProgLabel_400
  if(reportdetail>0)
  {
    SS2out << endl << "NUMBERS_AT_AGE" << endl;       // ProgLabel_410
    SS2out << "Area Bio_Pattern Gender BirthSeas SubMorph Morph Yr Seas Time Beg/Mid Era"<<age_vector <<endl;
    for (p=1;p<=pop;p++)
    for (g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
      {
      for (y=styr-2;y<=YrMax;y++)
      for (s=1;s<=nseas;s++)
       {
       t = styr+(y-styr)*nseas+s-1;
       temp=double(y)+double(s-1.)/nseas;
       SS2out <<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<g<<" "<<y<<" "<<s<<" "<<temp<<" B";
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}
       SS2out<<natage(t,p,g)<<endl;
       temp=double(y)+double(s-0.5)/nseas;
       SS2out <<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<g<<" "<<y<<" "<<s<<" "<<temp<<" M";
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}
       SS2out<<Save_PopAge(t,p+pop,g)<<endl;
       }
      }

    SS2out << endl << "NUMBERS_AT_LENGTH" << endl;
    SS2out << "Area Bio_Pattern Gender BirthSeas SubMorph Morph Yr Seas Time Beg/Mid Era "<<len_bins <<endl;
    for (p=1;p<=pop;p++)
    for (g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
      {
      for (y=styr;y<=endyr;y++)
      for (s=1;s<=nseas;s++)
       {
       t = styr+(y-styr)*nseas+s-1;
       temp=double(y)+double(s-1.)/nseas;
       SS2out <<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<g<<" "<<y<<" "<<s<<" "<<temp<<" B ";
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}
       SS2out<< Save_PopLen(t,p,g) << endl;
       temp=double(y)+double(s-0.5)/nseas;
       SS2out <<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<g<<" "<<y<<" "<<s<<" "<<temp<<" M ";
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}
       SS2out<< Save_PopLen(t,p+pop,g) << endl;
       }
      }

    SS2out << endl << "BIOMASS_AT_LENGTH" << endl;
    SS2out << "Area Bio_Pattern Gender BirthSeas SubMorph Morph Yr Seas Time Beg/Mid Era "<<len_bins <<endl;
    for (p=1;p<=pop;p++)
    for (g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
      {
      for (y=styr;y<=endyr;y++)
      for (s=1;s<=nseas;s++)
       {
       t = styr+(y-styr)*nseas+s-1;
       temp=double(y)+double(s-1.)/nseas;
       SS2out <<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<g<<" "<<y<<" "<<s<<" "<<temp<<" B ";
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}
       SS2out<< Save_PopWt(t,p,g) << endl;
       temp=double(y)+double(s-0.5)/nseas;
       SS2out <<p<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<g<<" "<<y<<" "<<s<<" "<<temp<<" M ";
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}
       SS2out<< Save_PopWt(t,p+pop,g) << endl;
       }
      }

     SS2out <<endl<< "CATCH_AT_AGE" << endl;              // ProgLabel_420
     SS2out << "Area Fleet Gender  XX XX Morph Yr Seas XX Era"<<age_vector <<endl;
     for (f=1;f<=Nfleet;f++)
     for (g=1;g<=gmorph;g++)
     {
     if(use_morph(g)>0)
     {
       for (y=styr-1;y<=endyr;y++)
       for (s=1;s<=nseas;s++)
       {
         t = styr+(y-styr)*nseas+s-1;
         SS2out <<fleet_area(f)<<" "<<f<<" "<<sx(g)<<" XX XX "<<g<<" "<<y<<" "<<s;
         if(y==styr-1)
           {SS2out<<" XX INIT ";}
         else if (y<=endyr)
           {SS2out<<" XX TIME ";}
         else
           {SS2out<<" XX FORE ";}
         SS2out<<catage(t,f,g)<< endl;
       }
     }
    }
  }

  SS2out <<endl<< "BIOLOGY "<<sum(use_morph)<<" "<<nlength<<" "<<nages<<" "<<nseas<<" N_Used_morphs;_lengths;_ages;_season;_by_season_in_endyr" << endl;
   SS2out<<"GP Bin Low Mean_Size Wt_len_F Mat_len Spawn Wt_len_M Fecundity"<<endl;
   for(gp=1;gp<=N_GP;gp++)
   for (z=1;z<=nlength;z++)
     {SS2out<<gp<<" "<<z<<" "<<len_bins(z)<<" "<<len_bins_m(z)<<" "<<wt_len(1,1,z)<<" "<<mat_len(gp,z)<<" "<<mat_fec_len(gp,z);
      if(gender==2) {SS2out<<" "<<wt_len(1,2,z);}
      SS2out<<" "<<fec_len(gp,z)<<endl;
     }

    SS2out<<endl<<"Natural_Mortality Method:_"<<natM_type<<endl<<"Bio_Pattern Gender BirthSeas Seas "<<age_vector<<endl;
      g=0;
      for(gg=1;gg<=gender;gg++)
      for(gp=1;gp<=N_GP;gp++)
      for(birthseas=1;birthseas<=nseas;birthseas++)
      {
        g++;
        if(use_morph(g)>0)
        {for(s=1;s<=nseas;s++) SS2out<<gp<<" "<<gg<<" "<<birthseas<<" "<<s<<" "<<natM(s,g)<<endl;}
      }

    if(Grow_type==3)  //  age-specific K
    {
    SS2out<<endl<<"Age_Specific_K"<<endl<<"Bio_Pattern Gender "<<age_vector<<endl;
      g=0;
      for(gg=1;gg<=gender;gg++)
      for(gp=1;gp<=N_GP;gp++)
      {
        g++;
        SS2out<<gp<<" "<<gg<<" "<<-VBK(g)<<endl;
      }
    }

   SS2out<<endl<<"Growth_Parameters"<<endl<<" Count Yr Gender Platoon A1 A2 L_a_A1 L_a_A2 K A_a_L0 Linf CVmin CVmax natM_amin natM_max M_age0 M_nages"
   <<" WtLen1 WtLen2 Mat1 Mat2 Fec1 Fec2"<<endl;
   for(g=1;g<=save_gparm_print;g++) {SS2out<<save_G_parm(g)(1,2)<<" "<<sx(save_G_parm(g,3))<<" "<<save_G_parm(g)(3,22)<<endl;}

   if(MGparm_doseas>0)
    {
   SS2out<<endl<<"Seas_Effects"<<endl<<"Seas F_wtlen1 F_wtlen2 F_mat1 F_mat2 F_fec1 F_fec2 M_wtlen1 M_wtlen2 L_a_A1 VBK"<<endl;
      for(s=1;s<=nseas;s++)
      {
        SS2out<<s<<" "<<save_seas_parm(s)<<endl;
      }
    }
    dvariable Herma_Cum;

//    restore_AgeLength_Key to endyr, otherwise it will have ALK from end of forecast
      if(time_vary_MG(endyr,2)>0 || time_vary_MG(endyr,3)>0 || WTage_rd>0)
      {
        y=endyr;
        t_base=styr+(y-styr)*nseas-1;
        for (s=1;s<=nseas;s++)
        {
          t = t_base+s;
          smid=s+nseas;
          bio_t=styr+(endyr-styr)*nseas+s-1;
          {
            Make_AgeLength_Key();        // indexed by season and gmorph
          }
        }
      }

   SS2out<<endl<<"Mat*Fecund_is_Fraction_mature_times_fecundity_at_age_used_to_create_Spawning_Biomass;_Mat_Bio_and_Mat_Num_omit_fecundity"<<endl;
   SS2out<<endl<<"Biology_at_age_in_endyr_with_";
   switch(CV_depvar)
   {
   case 0:
   {SS2out<<"CV=f(LAA)"; break;}
   case 1:
   {SS2out<<"CV=F(A)"; break;}
   case 2:
   {SS2out<<"SD=F(LAA)"; break;}
   case 3:
   {SS2out<<"SD=F(A)"; break;}
   case 4:
   {SS2out<<"logSD=f(A)"; break;}
  }
    
   SS2out<<endl;
   SS2out<<"Seas Morph Bio_Pattern Gender BirthSeas SubMorph Age Age_Beg Age_Mid M Len_Beg Len_Mid SD_Beg SD_Mid Wt_Beg Wt_Mid Len_Mat Age_Mat Mat*Fecund Mat_Bio Mat_Numbers";
   if(Hermaphro_Option>0) SS2out<<" Herma_Trans Herma_Cum ";
   for(f=1;f<=Ntypes;f++) SS2out<<" Len:_"<<f<<" SelWt:_"<<f<<" RetWt:_"<<f;
   SS2out<<endl;
   for (s=1;s<=nseas;s++)
   {
      t = styr+(endyr-styr)*nseas+s-1;
      smid=s+nseas;
     for(g=1;g<=gmorph;g++)
     if(use_morph(g)>0)
     {
     Herma_Cum=fracfemale;
     for(a=0;a<=nages;a++)
     {
      SS2out<<s<<" "<<g<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<a<<" "<<a+azero_seas(s)-azero_G(g)<<" "<<a+azero_seas(s)-azero_G(g)+seasdur_half(s);
      SS2out<<" "<<natM(s,GP3(g),a)<<" "<<Ave_Size(t,0,g,a)<<" "<<Ave_Size(t,1,g,a)<<" "
        <<Sd_Size_within(s,g,a)<<" "<<Sd_Size_within(s+nseas,g,a)<<" "
      <<Wt_Age_beg(s,g,a)<<" "<<Wt_Age_mid(s,g,a)<<" "<<ALK(s,g,a)*mat_len(GP4(g))<<" ";
      if(Maturity_Option<=2)
        {SS2out<<mat_age(GP4(g),a);}
      else if(sx(g)==1 && Maturity_Option<5)
        {SS2out<<Age_Maturity(GP4(g),a);}
      else
        {SS2out<<-1.;}
      SS2out<<" "<<fec(g,a)<<" "<<make_mature_bio(g,a)<<" "<<make_mature_numbers(g,a);
      if(Hermaphro_Option>0)
      {
        if(a>1) Herma_Cum*=(1.0-Hermaphro_val(GP4(g),a-1));
        SS2out<<" "<<Hermaphro_val(GP4(g),a)<<" "<<Herma_Cum;
      }
      if(WTage_rd==0)
      {
        for(f=1;f<=Ntypes;f++) SS2out<<
        " "<<ALK(smid,g,a)*elem_prod(sel_l(endyr,f,sx(g)),len_bins_m)/(ALK(smid,g,a)*sel_l(endyr,f,sx(g)))<<
        " "<<ALK(smid,g,a)*elem_prod(sel_l(endyr,f,sx(g)),wt_len(s,sx(g)))/(ALK(smid,g,a)*sel_l(endyr,f,sx(g)))<<
        " "<<ALK(smid,g,a)*elem_prod(sel_l_r(endyr,f,sx(g)),wt_len(s,sx(g)))/(ALK(smid,g,a)*sel_l_r(endyr,f,sx(g)));
      }
      else
      {
        for(f=1;f<=Ntypes;f++) SS2out<<
        " "<<ALK(smid,g,a)*elem_prod(sel_l(endyr,f,sx(g)),len_bins_m)/(ALK(smid,g,a)*sel_l(endyr,f,sx(g)))<<
        " "<<WTage_emp(t,GP3(g),f,a)<<" "<<WTage_emp(t,GP3(g),f,a);
      }
      SS2out<<endl;
      }}}


  SS2out <<endl<< "MEAN_BODY_WT(begin)";
  if(WTage_rd>0) SS2out<<" as read from wtatage.ss";
  SS2out<<" #NOTE_yr=_"<<styr-3<<"_stores_values_for_benchmark"<<endl;
  SS2out <<"Morph Yr Seas"<<age_vector<<endl;
  if(reportdetail>0)
  {

    for(g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
    {
    for(y=styr-3;y<=YrMax;y++)
    {
      yz=y;   if(yz>endyr+2) yz=endyr+2;
    if(y==styr-3 || y==styr || time_vary_MG(yz,2)>0 || time_vary_MG(yz,3)>0 || WTage_rd>0)  // if growth or wtlen parms have changed
    for(s=1;s<=nseas;s++)
     {
      t = styr+(y-styr)*nseas+s-1;
//       SS2out<<g<<" "<<GP4(g)<<" "<<sx(g)<<" "<<Bseas(g)<<" "<<GP2(g)<<" "<<y<<" "<<s<<" "<<Save_Wt_Age(t,g)<<endl;
       SS2out<<g<<" "<<y<<" "<<s<<" "<<Save_Wt_Age(t,g)<<endl;
     }
    }
  }
  }

  SS2out <<endl<< "MEAN_SIZE_TIMESERIES" << endl;           // ProgLabel_450
  SS2out <<"Morph Yr Seas Beg/Mid"<<age_vector<<endl;
  if(reportdetail>0)
  {
    for(g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
    {
      for(y=styr-3;y<=YrMax;y++)
      {
        yz=y;   if(yz>endyr+2) yz=endyr+2;
        if(y==styr-3 || y==styr ||  time_vary_MG(yz,2)>0)
        {
          for(s=1;s<=nseas;s++)
          {
            t = styr+(y-styr)*nseas+s-1;
            for(i=0;i<=1;i++)
            {
               SS2out<<g<<" "<<y<<" "<<s<<" "<<i<<" "<< Ave_Size(t,i,g)<<endl;
            }
           }
        }
      }
    }
    s=1;
    for (i=1;i<=gender;i++)
    {
      SS2out<<endl<<"mean_size_Jan_1_for_gender: "<<i<<" NOTE:_combines_all_birthseasons_areas_GP_and_submorphs"<<endl;
      SS2out <<"Gender Yr Seas Beg "<<age_vector<<endl;
      for(y=styr;y<=YrMax;y++)
      {
        yz=y;   if(yz>endyr+2) yz=endyr+2;
        if(y<=styr || time_vary_MG(yz,2)>0 || N_platoon>1)
        {
          t = styr+(y-styr)*nseas+s-1;
          SS2out<<i<<" "<<y<<" "<<s<<" "<<0;
          for (a=0;a<=nages;a++)
          {
            temp=0.0;
            temp1=0.0;
            for(g=1;g<=gmorph;g++)
            {
              if(sx(g)==i && use_morph(g)>0)
              {
                for(p=1;p<=pop;p++)
                {
                  temp+=natage(t,p,g,a);
                  temp1+=Ave_Size(t,0,g,a)*natage(t,p,g,a);
                }  // end loop of areas
              }  //  end need to use this gender/platoon
            }  //  end loop of all platoons
            if(temp>0.0) {SS2out <<" "<< temp1/temp;} else {SS2out<<" __";}
          }  //  end loop of ages
          SS2out<<endl;
        }  // end need to report this year
      }  // end year loop
    }  // end gender loop
  }  //   end do report detail

  if(reportdetail>0)
  {
  SS2out <<endl<< "AGE_LENGTH_KEY"<<" #midseason";
  if(Grow_logN==1) SS2out<<" #Lognormal ";
  SS2out<<endl;               // ProgLabel_460
  SS2out<<" sdratio "<<sd_ratio<<endl;
  SS2out<<" sdwithin "<<sd_within_platoon<<endl;
  SS2out<<" sdbetween "<<sd_between_platoon<<endl;
   for (s=1;s<=nseas;s++)
   for (int subseas=1;subseas<=2;subseas++)
   for (g=1;g<=gmorph;g++)
   if(use_morph(g)>0)
   {
    t = styr+(endyr-styr)*nseas+s-1;
    smid=s+nseas;
    SS2out <<endl<<" Seas: "<<s<<" Sub_Seas: "<<subseas<<"   Morph: "<<g<<endl;
//    SS2out <<endl<<" Seas: "<<s<<"   Morph: "<<g<<endl;
    SS2out <<"Age:";
    for (a=0;a<=nages;a++) SS2out << " "<<a;
    SS2out<<endl;
    for (z=nlength;z>=1;z--)
     {
      SS2out << len_bins2(z) << " ";
      if(subseas==1)
      {
        for (a=0;a<=nages;a++)  SS2out << ALK(s,g,a,z) << " " ;
      }
      else
      {
        for (a=0;a<=nages;a++)  SS2out << ALK(smid,g,a,z) << " " ;
      }
      SS2out<<endl;
     }
      SS2out<<"mean " << Ave_Size(t,subseas-1,g) << endl;
      if(subseas==1)
      {
        SS2out<<"sdsize " << Sd_Size_within(s,g) << endl;
      }
      else
      {
        SS2out<<"sdsize " << Sd_Size_within(smid,g) << endl;
      }
   }
  }

  if(reportdetail>0)
  {
    SS2out <<endl<< "AGE_AGE_KEY"<<endl;              // ProgLabel_470
    if(N_ageerr>0)
    {
      for (k=1;k<=N_ageerr;k++)
      {
      SS2out << "KEY: "<<k<<endl<< "mean " << age_err(k,1) << endl<< "SD " << age_err(k,2) << endl;
      for (b=n_abins;b>=1;b--)
       {
        SS2out << age_bins(b) << " ";
        for (a=0;a<=nages;a++)
          SS2out << age_age(k,b,a) << " " ;
        SS2out<<endl;
         }
       if(gender==2)
       {
       L2=n_abins;
       A2=nages+1;
      for (b=n_abins;b>=1;b--)
       {
        SS2out << age_bins(b) << " ";
        for (a=0;a<=nages;a++)
          SS2out << age_age(k,b+L2,a+A2) << " " ;
        SS2out<<endl;
         }
       }
      }
    }
    else
    {
      SS2out<<"no_age_error_key_used"<<endl;
    }
  }

//  ********  report the composition database
  SS_compout<<endl<<"Size_Bins_pop;_(Pop_len_mid_used_for_calc_of_selex_and_bio_quantities)"<<endl;
  SS_compout<<"Pop_Bin: ";
  for(j=1;j<=nlength;j++) SS_compout<<" "<<j;
  SS_compout<<endl<<"Length: "<<len_bins<<endl;
  SS_compout<<"Len_mid: "<<len_bins_m<<endl;
  SS_compout<<"Size_Bins_dat;_(Data_len_mid_for_reporting_only)"<<endl;
  SS_compout<<"Data_Bin: ";
  for(j=1;j<=nlen_bin;j++) SS_compout<<" "<<j;
  SS_compout<<endl<<"Length: "<<len_bins_dat<<endl;
  SS_compout<<"Len_mid: "<<len_bins_dat_m<<endl;

  SS_compout<<"Combine_males_with_females_thru_sizedata_bin: "<<CombGender_l;
  if(CombGender_l>0) SS_compout<<"   Size: "<<len_bins_dat(CombGender_l);
  SS_compout<<endl<<"Combine_males_with_females_thru_agedata_bin: "<<CombGender_a;
  if(CombGender_a>0) SS_compout<<"   Age: "<<age_bins(CombGender_a);
  SS_compout<<endl<<"Method_for_Lbin_definition_for_agecomp_data: "<<Lbin_method<<endl;

  SS_compout<<"For_Sizefreq,_Lbin_Lo_is_units(bio_or_numbers);_Lbin_hi_is_scale(kg,_lb,_cm,_in),_Ageerr_is_method"<<endl;
  SS_compout<<"Yr.S_is_midpoint_of_season_as_fraction_of_year"<<endl;
  SS_compout<<"For_Tag_output,_Rep_contains_Tag_Group,_Bin_is_fleet_for_TAG1_and_Bin_is_Year.Seas_for_TAG2"<<endl;
  SS_compout<<"Column_Super?_indicates_super-periods;_column_used_indicates_inclusion_in_logL"<<endl;

  SS_compout <<endl<< "Composition_Database" << endl;           // ProgLabel_480
  SS_compout<<"Yr Seas Yr.S Fleet Rep Pick_gender Kind Part Ageerr Gender Lbin_lo Lbin_hi Bin Obs Exp Pearson N effN Like Cum_obs Cum_exp SuprPer Used?"<<endl;
  int lasttime;
  int lastfleet;
  int repli;
  int N_out;
  N_out=0;
  lasttime=0;
  lastfleet=0;

  for(f=1;f<=Ntypes;f++)
  {
    in_superperiod=0;

 /*
    if(Nobs_l(f)>0)
    {
      for(z=1;z<=nlen_bin;z++)
      {
          SS_compout<<"all all all "<<f<<" 1 0 LEN "<<0<<" 0 "<<1<<" "<<
          1<<" "<<1<<" "<<len_bins_dat2(z)<<" "<<obs_l_all(1,f,z)<<" 0 0 0 0 0 0 0 _ _"<<endl;
      }
      SS_compout<<"all all all "<<f<<" 1 0 LEN "<<0<<" 0 "<<1<<" "<<1<<" "<<1<<endl;
      if(gender=2)
      {
        for(j=3;j<=4;j++)
        {
          for(z=1;z<=nlen_bin;z++)
          {
            SS_compout<<"all all all "<<f<<" 1 3 LEN "<<0<<" 0 "<<j-2<<" "<<
            1<<" "<<1<<" "<<len_bins_dat2(z)<<" "<<obs_l_all(1,f,z)<<" 0 0 0 0 0 0 0 _ _"<<endl;
          }
          SS_compout<<"all all all "<<f<<" 1 3 LEN "<<0<<" 0 "<<j-2<<" "<<1<<" "<<1<<endl;
        }
      }
    }
 */    
    for(i=1;i<=Nobs_l(f);i++)                          // loop obs in this type/time
    {
      N_out++;
      t=yr_l_t(f,i);
      y=yr_l_y(f,i);
      s=abs(header_l(f,i,2));
      temp2=0.0;
      temp1=0.0;
      repli=1;
      if(header_l(f,i,2)<0 && in_superperiod==0)
      {in_superperiod=1; anystring="Sup";}
      else if (header_l(f,i,2)<0 && in_superperiod>0)
      {anystring="Sup"; in_superperiod=0;}
      else if (in_superperiod>0)
      {in_superperiod++; anystring="Sup";}
      else
      {anystring="_";}
      if(header_l(f,i,3)<0 || nsamp_l(f,i)<0)
      {anystring+=" skip";}
      else
      {anystring+=" _";}
      if(gen_l(f,i)!=2)
      {
        s_off=1;
        for(z=tails_l(f,i,1);z<=tails_l(f,i,2);z++)
        {
            SS_compout<<y<<" "<<s<<" "<<header_l(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_l(f,i)<<" LEN "<<mkt_l(f,i)<<" 0 "<<s_off<<" "<<
            1<<" "<<1<<" "<<len_bins_dat2(z)<<" "<<obs_l(f,i,z)<<" "<<exp_l(f,i,z)<<" ";
            temp2+=obs_l(f,i,z);
            temp1+=exp_l(f,i,z);
            if(nsamp_l(f,i)>0 && header_l(f,i,3)>0)
            {
              if(exp_l(f,i,z)!=0.0 && exp_l(f,i,z)!=1.0)
              {SS_compout<<value((obs_l(f,i,z)-exp_l(f,i,z))/sqrt( exp_l(f,i,z) * (1.-exp_l(f,i,z)) / sfabs(nsamp_l(f,i))));}
              else
              {SS_compout<<" NA ";}
              SS_compout<<" "<<nsamp_l(f,i)<<" "<<neff_l(f,i)<<" ";
              if(obs_l(f,i,z)!=0.0 && exp_l(f,i,z)!=0.0)
              {SS_compout<<" "<<value(obs_l(f,i,z)*log(obs_l(f,i,z)/exp_l(f,i,z))*nsamp_l(f,i));}
              else
              {SS_compout<<" NA ";}
            }
            else
            {SS_compout<<" NA NA NA NA ";}
         SS_compout<<" "<<temp2<<" "<<temp1<<" "<<anystring<<endl;
        }

        SS_compout<<y<<" "<<s<<" "<<header_l(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_l(f,i)<<" LEN "
        <<mkt_l(f,i)<<" 0 "<<s_off<<" "<<1<<" "<<1<<endl;
      }
      if(gen_l(f,i)>=2 && gender==2)  // do males
      {
        s_off=2;
        for(z=tails_l(f,i,3);z<=tails_l(f,i,4);z++)
        {
           SS_compout<<y<<" "<<s<<" "<<header_l(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_l(f,i)<<" LEN "<<mkt_l(f,i)<<" 0 "<<s_off<<" "<<
           1<<" "<<nlength<<" "<<len_bins_dat2(z)<<" "<<obs_l(f,i,z)<<" "<<exp_l(f,i,z)<<" ";
           temp2+=obs_l(f,i,z);
           temp1+=exp_l(f,i,z);
           if(nsamp_l(f,i)>0 && header_l(f,i,3)>0)
           {
              if(exp_l(f,i,z)!=0.0 && exp_l(f,i,z)!=1.0)
          {SS_compout<<value((obs_l(f,i,z)-exp_l(f,i,z))/sqrt( exp_l(f,i,z) * (1-exp_l(f,i,z)) / sfabs(nsamp_l(f,i))));}
          else
          {SS_compout<<" NA ";}
          SS_compout<<" "<<nsamp_l(f,i)<<" "<<neff_l(f,i)<<" ";
          if(obs_l(f,i,z)!=0.0 && exp_l(f,i,z)!=0.0)
          {SS_compout<<" "<<value(obs_l(f,i,z)*log(obs_l(f,i,z)/exp_l(f,i,z))*nsamp_l(f,i));}
          else
          {SS_compout<<" NA ";}
          }
          else
          {SS_compout<<" NA NA NA NA ";}
         SS_compout<<" "<<temp2<<" "<<temp1<<" "<<anystring<<endl;
        }
        SS_compout<<y<<" "<<s<<" "<<header_l(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_l(f,i)<<" LEN "
        <<mkt_l(f,i)<<" 0 "<<s_off<<" "<<1<<" "<<1<<endl;
      }
    }

      {
       for(i=1;i<=Nobs_a(f);i++)                          // loop all obs in this type
       {
        N_out++;
        temp2=0.0;
        temp1=0.0;
      t=yr_a_t(f,i);
      y=yr_a_y(f,i);
      s=abs(header_a(f,i,2));
      repli=1;
      if(header_a(f,i,2)<0 && in_superperiod==0)
      {in_superperiod=1; anystring="Sup";}
      else if (header_a(f,i,2)<0 && in_superperiod>0)
      {anystring="Sup"; in_superperiod=0;}
      else if (in_superperiod>0)
      {in_superperiod++; anystring="Sup";}
      else
      {anystring="_";}
      if(header_a(f,i,3)<0)
      {anystring+=" skip";}
      else
      {anystring+=" _";}

        if(gen_a(f,i)!=2)
         {s_off=1;
         for(z=tails_a(f,i,1);z<=tails_a(f,i,2);z++)
          {SS_compout<<y<<" "<<s<<" "<<header_a(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_a(f,i)<<" AGE "<<mkt_a(f,i)<<" "<<ageerr_type_a(f,i)
         <<" "<<s_off<<" "<<len_bins(Lbin_lo(f,i))<<" "<<len_bins(Lbin_hi(f,i))<<" "<<age_bins(z)<<" "<<obs_a(f,i,z)<<" " <<exp_a(f,i,z)<<" ";
           temp2+=obs_a(f,i,z);
           temp1+=exp_a(f,i,z);
          if(header_a(f,i,3)>0)
          {
            if(exp_a(f,i,z)!=0.0 && exp_a(f,i,z)!=1.0)
            {SS_compout<<value((obs_a(f,i,z)-exp_a(f,i,z))/sqrt( exp_a(f,i,z) * (1.0-exp_a(f,i,z)) / sfabs(nsamp_a(f,i))));}
            else
            {SS_compout<<" NA ";}
            SS_compout<<" "<<nsamp_a(f,i)<<" "<<neff_a(f,i)<<" ";
            if(obs_a(f,i,z)!=0.0 && exp_a(f,i,z)!=0.0)
            {SS_compout<<" "<<value(obs_a(f,i,z)*log(obs_a(f,i,z)/exp_a(f,i,z))*nsamp_a(f,i));}
            else
            {SS_compout<<" NA ";}
          }
          else
          {SS_compout<<" NA NA NA NA ";}
         SS_compout<<" "<<temp2<<" "<<temp1<<" "<<anystring<<endl;
        }

        SS_compout<<y<<" "<<s<<" "<<header_a(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_a(f,i)<<" AGE "
         <<mkt_a(f,i)<<" "<<ageerr_type_a(f,i)<<" "<<s_off<<" "<<1<<" "<<nlength<<endl;}

        if(gen_a(f,i)>=2 && gender==2)  // do males
         {s_off=2;
         for(z=tails_a(f,i,3);z<=tails_a(f,i,4);z++)
          {SS_compout<<y<<" "<<s<<" "<<header_a(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_a(f,i)<<" AGE "<<mkt_a(f,i)<<" "<<ageerr_type_a(f,i)<<" "<<s_off
         <<" "<<len_bins(Lbin_lo(f,i))<<" "<<len_bins(Lbin_hi(f,i))<<" "<<age_bins(z)<<" "<<obs_a(f,i,z)<<" "<<exp_a(f,i,z)<<" ";
           temp2+=obs_a(f,i,z);
           temp1+=exp_a(f,i,z);
          if(header_a(f,i,3)>0)
          {
            if(exp_a(f,i,z)!=0.0 && exp_a(f,i,z)!=1.0)
            {SS_compout<<value((obs_a(f,i,z)-exp_a(f,i,z))/sqrt( exp_a(f,i,z) * (1.-exp_a(f,i,z)) / sfabs(nsamp_a(f,i))));}
            else
            {SS_compout<<" NA ";}
            SS_compout<<" "<<nsamp_a(f,i)<<" "<<neff_a(f,i)<<" ";
            if(obs_a(f,i,z)!=0.0 && exp_a(f,i,z)!=0.0)
            {SS_compout<<" "<<value(obs_a(f,i,z)*log(obs_a(f,i,z)/exp_a(f,i,z))*nsamp_a(f,i));}
            else
            {SS_compout<<" NA ";}
          }
          else
          {
            SS_compout<<" NA NA NA NA ";}
            SS_compout<<" "<<temp2<<" "<<temp1<<" "<<anystring<<endl;
          }
        SS_compout<<y<<" "<<s<<" "<<header_a(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_a(f,i)<<" AGE "
         <<mkt_a(f,i)<<" "<<ageerr_type_a(f,i)<<" "<<s_off<<" "<<1<<" "<<nlength<<endl;}
       }
      }  //end have agecomp data

      {
       for(i=1;i<=Nobs_ms(f);i++)
       {
        N_out++;
      t=yr_ms_t(f,i);
      y=yr_ms_y(f,i);
      s=abs(header_ms(f,i,2));
      repli=1;
      if(header_ms(f,i,2)<0 && in_superperiod==0)
      {in_superperiod=1; anystring="Sup";}
      else if (header_ms(f,i,2)<0 && in_superperiod>0)
      {anystring="Sup"; in_superperiod=0;}
      else if (in_superperiod>0)
      {in_superperiod++; anystring="Sup";}
      else
      {anystring="_";}
      if(use_ms(f,i)<0)
      {anystring+=" skip";}
      else
      {anystring+=" _";}

       for(z=1;z<=n_abins2;z++)
       {
        if(z<=n_abins) s_off=1; else s_off=2;
        t1=obs_ms_n(f,i,z);
        //  whre:  obs_ms_n(f,i,z)=sqrt(var_adjust(6,f)*obs_ms_n(f,i,z));
        if(ageerr_type_ms(f,i)>0) {anystring2=" L@A ";} else {anystring2=" W@A ";}
        if(t1>0.) t1=square(t1);
        SS_compout<<y<<" "<<s<<" "<<header_ms(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_ms(f,i)<<anystring2<<mkt_ms(f,i)<<" "<<
         ageerr_type_ms(f,i)<<" "<<s_off<<" "<<exp_ms_sq(f,i,z)<<" "<<nlen_bin<<" "<<age_bins(z)<<" "<<
         obs_ms(f,i,z)<<" "<<exp_ms(f,i,z)<<" ";
        if(obs_ms(f,i,z)>0. && t1>0. && use_ms(f,i)>0)
        {
          SS_compout<<(obs_ms(f,i,z) -exp_ms(f,i,z)) / (exp_ms_sq(f,i,z)/obs_ms_n(f,i,z))<<" ";  // Pearson
          SS_compout<<t1<<" ";  // sample size
          SS_compout<<square(1.0/((obs_ms(f,i,z) -exp_ms(f,i,z)) / exp_ms_sq(f,i,z)))<<" "; // effectiove sample size
          SS_compout<<0.5*square((obs_ms(f,i,z) -exp_ms(f,i,z)) / (exp_ms_sq(f,i,z)/obs_ms_n(f,i,z))) + sd_offset*log(exp_ms_sq(f,i,z)/obs_ms_n(f,i,z)); //  -logL
        }
        else
        {
          SS_compout<<" NA "<<t1<<" NA NA ";
        }
        SS_compout<<" NA NA "<<anystring<<endl;
        if(z==n_abins || z==n_abins2) SS_compout<<y<<" "<<s<<" "<<header_ms(f,i,0)<<" "<<f<<" "<<repli<<" "<<gen_ms(f,i)<<
         anystring2<<mkt_ms(f,i)<<" "<<ageerr_type_ms(f,i)<<" "<<s_off<<" "<<1<<" "<<nlen_bin<<endl;
       }
      }  //end have data
      }
  }  // end fleet

    if(SzFreq_Nmeth>0)       //  have some sizefreq data
    {
      in_superperiod=0;
      for(iobs=1;iobs<=SzFreq_totobs;iobs++)
      {
        y=SzFreq_obs_hdr(iobs,1);
        if(y>=styr && y<=retro_yr)  // flag for obs that are used
        {
          N_out++;
          temp2=0.0;
          temp1=0.0;
          s=abs(SzFreq_obs_hdr(iobs,2));
//          temp=float(y)+float(abs(s)-1.)/float(nseas);
          temp = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
          f=abs(SzFreq_obs_hdr(iobs,3));
          gg=SzFreq_obs_hdr(iobs,4);  // gender
          k=SzFreq_obs_hdr(iobs,6);
          if(SzFreq_obs_hdr(iobs,2)<0 && in_superperiod==0)
          {in_superperiod=1; anystring="Sup";}
          else if (SzFreq_obs_hdr(iobs,2)<0 && in_superperiod>0)
          {anystring="Sup"; in_superperiod=0;}
          else if (in_superperiod>0)
          {in_superperiod++; anystring="Sup";}
          else
          {anystring="_";}
          if(SzFreq_obs_hdr(iobs,3)<0)
          {anystring+=" skip";}
          else
          {anystring+=" _";}
          p=SzFreq_obs_hdr(iobs,5);  // partition
          z1=SzFreq_obs_hdr(iobs,7);
          z2=SzFreq_obs_hdr(iobs,8);
          temp2=0.0;
          temp1=0.0;
          for(z=z1;z<=z2;z++)
          {
            s_off=1;
            SS_compout<<y<<" "<<s<<" "<<temp<<" "<<f<<" "<<1<<" "<<gg<<" SIZE "<<p<<" "<<k;
            if(z>SzFreq_Nbins(k)) s_off=2;
            SS_compout<<" "<<s_off<<" "<<SzFreq_units(k)<<" "<<SzFreq_scale(k)<<" ";
            if(s_off==1) {SS_compout<<SzFreq_bins1(k,z);} else {SS_compout<<SzFreq_bins1(k,z-SzFreq_Nbins(k));}
            SS_compout<<" "<<SzFreq_obs(iobs,z)<<" " <<SzFreq_exp(iobs,z)<<" ";
            temp2+=SzFreq_obs(iobs,z);
            temp1+=SzFreq_exp(iobs,z);
            if(SzFreq_obs_hdr(iobs,3)>0)
            {
              if(SzFreq_exp(iobs,z)!=0.0 && SzFreq_exp(iobs,z)!=1.0)
                {SS_compout<<(SzFreq_obs(iobs,z)-SzFreq_exp(iobs,z))/sqrt( SzFreq_exp(iobs,z) * (1.-SzFreq_exp(iobs,z)) / SzFreq_sampleN(iobs));}
              else
                {SS_compout<<" NA ";}
              SS_compout<<" "<<SzFreq_sampleN(iobs)<<" "<<SzFreq_effN(iobs)<<" ";
              if(SzFreq_obs(iobs,z)!=0.0 && SzFreq_exp(iobs,z)!=0.0)
                {SS_compout<<" "<<SzFreq_obs(iobs,z)*log(SzFreq_obs(iobs,z)/SzFreq_exp(iobs,z))*SzFreq_sampleN(iobs);}
              else
                {SS_compout<<" NA ";}
            }
            else
            {SS_compout<<" NA NA NA NA ";}
            SS_compout<<" "<<temp2<<" "<<temp1<<" "<<anystring<<endl;
            if(z==z2 || z==SzFreq_Nbins(k))
            SS_compout<<y<<" "<<s<<" "<<temp<<" "<<f<<" "<<repli<<" "<<gg<<" SIZE "<<p<<" "<<k<<" "<<s_off<<" "<<1<<" "<<2<<endl;
          }
        }
      }
    }
    if(Do_Morphcomp>0)
    {
      for(iobs=1;iobs<=Morphcomp_nobs;iobs++)
      {
        N_out++;
        y=Morphcomp_obs(iobs,1); s=Morphcomp_obs(iobs,2);
        temp1=s-1.;
        temp2=y;
        temp = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
//        temp=temp2+temp1/nseas;
        f=Morphcomp_obs(iobs,3);
        k=5+Morphcomp_nmorph;
        for(z=6;z<=k;z++)
        {
          SS_compout<<y<<" "<<s<<" "<<temp<<" "<<1<<" "<<1<<" GP% "<<0<<" "<<Morphcomp_obs(iobs,5);
         SS_compout<<" "<<0<<" "<<0<<" "<<0<<" "<<z-5<<" "<<Morphcomp_obs(iobs,z)<<" " <<Morphcomp_exp(iobs,z)<<" "<<endl;
        }
      }
    }


    if(Do_TG>0)
    {
      for(TG=1;TG<=N_TG;TG++)
      {
        y=TG_release(TG,3); s=TG_release(TG,4);
        for(TG_t=0;TG_t<=TG_endtime(TG);TG_t++)
        {
          N_out++;
          t = styr+(y-styr)*nseas+s-1;
          temp1=s-1.;
//          temp=float(y)+temp1/float(nseas);
          temp = float(y)+0.01*int(100.*(azero_seas(s)+seasdur_half(s)));
          SS_compout<<y<<" "<<s<<" "<<temp<<" "<<0<<" "<<TG<<" "<<TG_release(TG,6)<<" TAG2 NA NA NA NA NA "<<
          temp<<" "<<TG_recap_obs(TG,TG_t,0)<<" "<<TG_recap_exp(TG,TG_t,0)<<" NA NA NA NA NA NA NA ";
          if(TG_t>=TG_mixperiod) {SS_compout<<"_"<<endl;} else {SS_compout<<" skip"<<endl;}
          if(Nfleet>1)
          for(f=1;f<=Nfleet;f++)
          {
            SS_compout<<y<<" "<<s<<" "<<temp<<" "<<f<<" "<<TG<<" "<<TG_release(TG,6)<<" TAG1 NA NA NA NA NA "<<
            f<<" "<<TG_recap_obs(TG,TG_t,f)<<" "<<TG_recap_exp(TG,TG_t,f)<<" NA "<<TG_recap_obs(TG,TG_t,0)
            <<" NA NA NA NA NA ";
          if(TG_t>=TG_mixperiod) {SS_compout<<"_"<<endl;} else {SS_compout<<" skip"<<endl;}
          }
          s++; if(s>nseas) {s=1; y++;}
        }
      }
    }

  if(N_out==0) SS_compout<<styr<<" -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1"<<endl;
  SS_compout<<"-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1"<<endl<<"End_comp_data"<<endl;
  SS_compout<<"end "<<endl;

  SS2out <<endl<< "SELEX_database" << endl;
  SS2out<<"Fleet Yr Kind Gender Bin Selex"<<endl;

  if(reportdetail<=0)
  {
    SS2out<<"1 1990 L 1 30 .5"<<endl;
  }
  else
  {
  for (f=1;f<=Ntypes;f++)
  for (y=styr-3;y<=endyr;y++)
  {
   if(y==styr-3 || y==endyr || (time_vary_sel(y,f)>0 || time_vary_sel(y+1,f)>0))
   {
    for (gg=1;gg<=gender;gg++)
    {
     for(z=1;z<=nlength;z++) {SS2out<<f<<" "<<y<<" L "<<gg<<" "<<len_bins(z)<<" "<<sel_l(y,f,gg,z)<<endl;}
     if(seltype(f,2)!=0)
     {
      if(gg==1) {for(z=1;z<=nlength;z++) {SS2out<<f<<" "<<y<<" D "<<gg<<" "<<len_bins(z)<<" "<<retain(y,f,z)<<endl;}}
      else
      {for(z=1;z<=nlength;z++) {SS2out<<f<<" "<<y<<" D "<<gg<<" "<<len_bins(z)<<" "<<retain(y,f,z+nlength)<<endl;}}
     }
     if(seltype(f,2)==2)
     {
      if(gg==1) {for(z=1;z<=nlength;z++) {SS2out<<f<<" "<<y<<" DM "<<gg<<" "<<len_bins(z)<<" "<<discmort(y,f,z)<<endl;}}
      else
      {for(z=1;z<=nlength;z++) {SS2out<<f<<" "<<y<<" DM "<<gg<<" "<<len_bins(z)<<" "<<discmort(y,f,z+nlength)<<endl;}}
     }
    }
   }
   if(time_vary_sel(y,f+Ntypes)>0)
   {
    for (gg=1;gg<=gender;gg++)
    for(a=0;a<=nages;a++) {SS2out<<f<<" "<<y<<" A "<<gg<<" "<<a<<" "<<sel_a(y,f,gg,a)<<endl;}
   }
  }
  }  // end do report detail
  SS2out<<" end selex output "<<endl;

// ******************************************************
//  Do Btarget profile
  if(Do_Benchmark>0)
  {
        SS2out<<endl<<"SPR/YPR_Profile "<<endl<<"Iter Fmult F_std SPR Catch SSB Recruits SSB/Bzero YPR";
        for(f=1;f<=Nfleet;f++) {SS2out<<" catch:_"<<f;}
        for(p=1;p<=pop;p++)
        for(gp=1;gp<=N_GP;gp++)
        {SS2out<<" Area:"<<p<<"_GP:"<<gp;}
        SS2out<<endl;
        int SPRloop;
        int bio_t_base;
        dvariable Fmult2=maxpossF;
        dvariable Fcrash=Fmult2;
        dvariable Fmultchanger0=value(Fmult2/39.);
        dvariable Fmultchanger1;
        dvariable Fmultchanger2;
        dvariable Btgt_prof;
        dvariable Btgt_prof_rec;
        dvariable SPR_last;
        dvariable SPR_trial;

    y=styr-3;
    yz=y;
    bio_yr=y;
    eq_yr=y;
    t_base=y+(y-styr)*nseas-1;
    bio_t_base=styr+(bio_yr-styr)*nseas-1;

    for(s=1;s<=nseas;s++)
    {
      t = styr-3*nseas+s-1;
      Make_AgeLength_Key();
      Save_Wt_Age(t)=Wt_Age_beg(s);
    }

    for (s=1;s<=nseas;s++)
    for (g=1;g<=gmorph;g++)
    if(use_morph(g)>0)
    {
      smid=s+nseas;
      Make_FishSelex();
    }
    equ_Recr=1.0;
    Fishon=0;
    Do_Equil_Calc();
    SPR_unf=SPB_equil;

        for(int SPRloop1=0; SPRloop1<=2; SPRloop1++)
        {
          Fmultchanger1=value(pow(0.0001/Fcrash,0.025));
          Fmultchanger2=value(Fcrash/39.);
          if(SPRloop1==1)  Fmult2=Fcrash;
          for(SPRloop=1; SPRloop<=40; SPRloop++)
          {
            for(f=1;f<=Nfleet;f++)
            for(s=1;s<=nseas;s++)
            {
              t=bio_t_base+s;
              Hrate(f,t)=Fmult2*Bmark_RelF_Use(s,f);
            }

            Fishon=1;
            Do_Equil_Calc();
            SPR_temp=SPB_equil;
            Get_EquilCalc = Equil_Spawn_Recr_Fxn();   // call  function
            Btgt_prof=Get_EquilCalc(1);
            Btgt_prof_rec=Get_EquilCalc(2);
            if(SPRloop1==0)
            {
              if(Btgt_prof<0.001 && Btgt_prof_rec<0.001)
              {Fcrash=Fmult2;}
            }
            SS2out<<SPRloop<<" "<<Fmult2<<" "<<equ_F_std<<" "<<value(SPB_equil/SPR_unf)<<" "
            <<value(YPR_dead*Btgt_prof_rec)<<" "<<Btgt_prof<<" "<<Btgt_prof_rec<<" "<<value(Btgt_prof/SPB_virgin)
            <<" "<<value(YPR_dead)<<" "<<value(column(equ_catch_fleet,2)*Btgt_prof_rec);
            for(p=1;p<=pop;p++)
            for(gp=1;gp<=N_GP;gp++)
            {SS2out<<" "<<SPB_equil_pop_gp(p,gp)*Btgt_prof_rec;}
            SS2out<<endl;
            if(SPRloop1==0)
              {Fmult2-=Fmultchanger0;}
            else if(SPRloop1==1)
              {Fmult2*=Fmultchanger1;}
            else if(SPRloop1==2)
              {Fmult2+=Fmultchanger2;}
          }
        }  // end Fmult profile

        SPR_trial=value(SPB_equil/SPR_unf);
        SPR_last=SPR_trial*2.;
        SPRloop=0;
        dvariable YPR_last;
        YPR_last=-1.;
        while (SPR_trial>0.001 && SPR_last>1.00001*SPR_trial && YPR_last<YPR_dead && SPRloop<1000)
        {
          if(F_Method>1)
          {Fmult2*=1.05;}
          else
          {Fmult2=Fmult2+(1.0-Fmult2)*0.05;}
          SPR_last=SPR_trial;
          YPR_last=YPR_dead;
          for(f=1;f<=Nfleet;f++)
          for(s=1;s<=nseas;s++)
          {
            t=bio_t_base+s;
            Hrate(f,t)=Fmult2*Bmark_RelF_Use(s,f);
          }
          SPRloop+=1;
          Fishon=1;
          Do_Equil_Calc();
          SPR_temp=SPB_equil;
          Get_EquilCalc = Equil_Spawn_Recr_Fxn();   // call  function
          Btgt_prof=Get_EquilCalc(1);
          Btgt_prof_rec=Get_EquilCalc(2);
          SPR_trial=value(SPB_equil/SPR_unf);
          SS2out<<SPRloop<<" "<<Fmult2<<" "<<equ_F_std<<" "<<SPR_trial<<" "
          <<value(YPR_dead*Btgt_prof_rec)<<" "<<Btgt_prof<<" "<<Btgt_prof_rec<<" "<<value(Btgt_prof/SPB_virgin)
          <<" "<<value(YPR_dead)<<" "<<value(column(equ_catch_fleet,2)*Btgt_prof_rec);
          for(p=1;p<=pop;p++)
          for(gp=1;gp<=N_GP;gp++)
          {SS2out<<" "<<SPB_equil_pop_gp(p,gp)*Btgt_prof_rec;}
          SS2out<<endl;
        }
        // end Btarget profile
    }
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 41 write_Bzero_output */
FUNCTION void write_Bzero_output()
//  output annual time series for beginning of year and summing across areas for each GP and gender
  {
    for(fishery_on_off=1;fishery_on_off>=0;fishery_on_off--)
    {
    SS2out<<endl<<"Dynamic_Bzero"<<endl;
    SS2out<<"Spawning_Biomass_Report";
    if(fishery_on_off==0) {SS2out<<"_1 No_fishery_for_Z=M_and_dynamic_Bzero";} else {SS2out<<"_2 With_fishery";}
    SS2out<<endl<<"Yr Area: ";
    for(p=1;p<=pop;p++)
    for(gp=1;gp<=N_GP;gp++)
    {SS2out<<p<<" ";}
    SS2out<<endl<<"xxxx GP: ";
    for(p=1;p<=pop;p++)
    for(gp=1;gp<=N_GP;gp++)
    {SS2out<<gp<<" ";}
    SS2out<<endl;

      save_gparm=0;
      if(fishery_on_off==0)
      {
        get_initial_conditions();
        get_time_series();
        if(Do_Forecast>0)
        {
          show_MSY=0;
          Get_Forecast();
        }
      }

    for(y=styr-2;y<=YrMax;y++)
    {
      SS2out<<y;
       if(y==styr-2)
         {SS2out<<" VIRG ";}
       else if (y==styr-1)
         {SS2out<<" INIT ";}
       else if (y<=endyr)
         {SS2out<<" TIME ";}
       else
         {SS2out<<" FORE ";}
      for(p=1;p<=pop;p++)
      for(gp=1;gp<=N_GP;gp++)
      {SS2out<<" "<<SPB_pop_gp(y,p,gp);}
      SS2out<<endl;
    }

    SS2out << endl << "NUMBERS_AT_AGE_Annual";
    if(fishery_on_off==0) {SS2out<<"_1 No_fishery_for_Z=M_and_dynamic_Bzero";} else {SS2out<<"_2 With_fishery";}
    SS2out << endl;
    SS2out << "Bio_Pattern Gender Year "<<age_vector <<endl;
    dvector tempvec2(1,nages);
    for(gg=1;gg<=gender;gg++)
    for(gp=1;gp<=N_GP;gp++)
    for(y=styr;y<=YrMax;y++)
    {
      tempvec_a.initialize();
      t = styr+(y-styr)*nseas;  // first season only
      for (p=1;p<=pop;p++)
      for (g=1;g<=gmorph;g++)
      if(use_morph(g)>0)
      {
        if(GP4(g)==gp && sx(g)==gg) tempvec_a+= value(natage(t,p,g));
      }
      if(nseas>1)
      {
        tempvec_a(0)=0.;
        for(s=1;s<=nseas;s++)
        for (p=1;p<=pop;p++)
        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0 && Bseas(g)==s)
        {
          if(GP4(g)==gp && sx(g)==gg) tempvec_a(0) += value(natage(t,p,g,0));
        }
      }
      SS2out <<gp<<" "<<gg<<" "<<y<<" "<<tempvec_a<<endl;
    }

    SS2out << endl << "Z_AT_AGE_Annual";
    if(fishery_on_off==0) {SS2out<<"_1 No_fishery_for_Z=M_and_dynamic_Bzero";} else {SS2out<<"_2 With_fishery";}
    if(Hermaphro_Option>0) SS2out<<"_hermaphrodites_combined_gender_output";
    SS2out << endl;
    SS2out << "Bio_Pattern Gender Year "<<age_vector <<endl;
    if(Hermaphro_Option>0)
    {k=1;}
    else
    {k=gender;}
    for(gg=1;gg<=k;gg++)
    for(gp=1;gp<=N_GP;gp++)
    for(y=styr;y<=YrMax;y++)
    {
      tempvec_a.initialize();
      t = styr+(y-styr)*nseas;  // first season only
      for (p=1;p<=pop;p++)
      for (g=1;g<=gmorph;g++)
      if(use_morph(g)>0)
      {
        if(GP4(g)==gp && (sx(g)==gg || Hermaphro_Option>0)) tempvec_a+= value(natage(t,p,g));
      }
      if(nseas>1)
      {
        tempvec_a(0)=0.;
        for(s=1;s<=nseas;s++)
        for (p=1;p<=pop;p++)
        for (g=1;g<=gmorph;g++)
        if(use_morph(g)>0 && Bseas(g)==s)
        {
          if(GP4(g)==gp && (sx(g)==gg || Hermaphro_Option>0)) tempvec_a(0) += value(natage(t,p,g,0));
        }
      }
      if(y>styr)
      {
      SS2out <<gp<<" "<<gg<<" "<<y-1<<" "<<log(elem_div(tempvec2(1,nages),tempvec_a(1,nages)))<<" _ "<<endl;
      }
      for(a=0;a<=nages-1;a++) tempvec2(a+1)=value(tempvec_a(a));
      tempvec2(nages)+=value(tempvec_a(nages));
    }
    }
    SS2out<<" Note:  Z calculated as -ln(Nt+1 / Nt)"<<endl;
    SS2out<<" Note:  Z calculation for maxage-1 includes numbers at maxage, so is approximate"<<endl;
    if(nseas>1) SS2out<<" Age zero fish summed across birthseasons, but Z calc is as if all born in season 1"<<endl;
    fishery_on_off=1;
  }
  
//*********************************************************************
 /*  ProgLabel_FUNCTION 42 Join_Fxn */
FUNCTION dvariable Join_Fxn(const prevariable& MinPoss, const prevariable& MaxPoss, const prevariable& Inflec, const prevariable& Xvar, const prevariable& Y1, const prevariable& Y2)
  {
  dvariable Yresult;
  dvariable join;
  join=1.000/(1.000+mfexp(1000.0*(Xvar-Inflec)/(MaxPoss-MinPoss)));  //  steep joiner at the inflection
  Yresult=Y1*(join)+Y2*(1.000-join);
  return Yresult;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 43 Spawn_Recr */
FUNCTION dvariable Spawn_Recr(const prevariable& SPB_current)
  {
    dvariable NewRecruits;
    dvariable SPB_BH1;
    dvariable Recr_virgin_adj;
    dvariable SPB_virgin_adj;
    dvariable steepness;
    dvariable SPB_curr_adj;
    dvariable join;
    dvariable SRZ_0;
    dvariable SRZ_max;
    dvariable SRZ_surv;

//  ProgLabel_43.1  add 0.1 to input spawning biomass value to make calculation more rebust
    SPB_curr_adj = SPB_current + 0.100;   // robust
    Recr_virgin_adj=Recr_virgin;  SPB_virgin_adj=SPB_virgin;
    steepness=SR_parm(2);
    
//  ProgLabel_43.2  adjust for environmental effects on S-R parameters: Rzero or steepness
    if(SR_env_target==2)
    {
      Recr_virgin_adj*=mfexp(SR_parm(N_SRparm2-2)* env_data(y,SR_env_link));
      SPB_virgin_adj*=mfexp(SR_parm(N_SRparm2-2)* env_data(y,SR_env_link));
    }
    else if(SR_env_target==3)
    {
      temp=log((SRvec_HI(2)-SRvec_LO(2)+0.0000002)/(SR_parm(2)-SRvec_LO(2)+0.0000001)-1.)/(-2.);
      temp+=SR_parm(N_SRparm2-2)* env_data(y,SR_env_link);
      steepness=SRvec_LO(2)+(SRvec_HI(2)-SRvec_LO(2))/(1.+mfexp(-2.*temp));
    }

// functions below use Recr_virgin_adj,SPB_virgin_adj which could be adjusted to differ from R0,SPB_virgin
//  ProgLabel_43.3  calculate expected recruitment from the input spawning biomass and the SR curve
    switch(SR_fxn)
    {
      case 1: // previous placement for B-H constrained
      {
        N_warn++; cout<<"Critical error:  see warning"<<endl; warning<<"B-H constrained curve is now Spawn-Recr option #6"<<endl; exit(1);
        break;
      }
//  ProgLabel_43.3.2  Ricker
      case 2:  // ricker
      {
        NewRecruits = Recr_virgin_adj*SPB_curr_adj/SPB_virgin_adj * mfexp(steepness*(1.-SPB_curr_adj/SPB_virgin_adj));
        break;
      }
//  ProgLabel_43.3.3  Beverton-Holt
      case 3: // Beverton-Holt
      {
        NewRecruits =  (4.*steepness*Recr_virgin_adj*SPB_curr_adj) / (SPB_virgin_adj*(1.-steepness)+(5.*steepness-1.)*SPB_curr_adj);
        break;
      }
//  ProgLabel_43.3.4  constant expected recruitment
      case 4:  // none
      {
        NewRecruits=Recr_virgin_adj;
        break;
      }
//  ProgLabel_43.3.5  Hockey stick
      case 5:  // hockey stick  where "steepness" is now the fraction of B0 below which recruitment declines linearly
               //  the 3rd parameter allows for a minimum recruitment level
      {
        temp=SR_parm(3)*Recr_virgin_adj + SPB_curr_adj/(steepness*SPB_virgin_adj)*(Recr_virgin_adj-SR_parm(3)*Recr_virgin_adj);  //  linear decrease below steepness*SPB_virgin_adj
        NewRecruits=Join_Fxn(0.0*SPB_virgin_adj,SPB_virgin_adj,steepness*SPB_virgin_adj, SPB_curr_adj, temp, Recr_virgin_adj);
        break;
      }

//  ProgLabel_43.3.6  Beverton-Holt, with constraint to have constant R about Bzero
      case 6: //Beverton-Holt constrained
      {
        if(SPB_curr_adj>SPB_virgin_adj) {SPB_BH1=SPB_virgin_adj;} else {SPB_BH1=SPB_curr_adj;}
        NewRecruits=(4.*steepness*Recr_virgin_adj*SPB_BH1) / (SPB_virgin_adj*(1.-steepness)+(5.*steepness-1.)*SPB_BH1);
        break;
      }

//  ProgLabel_43.3.7  survival based
      case 7:  // survival based, so constrained such that recruits cannot exceed fecundity
      {
        // PPR_0=SPB_virgin_adj/Recr_virgin_adj;  //  pups per recruit at virgin
        // Surv_0=1./PPR_0;   //  recruits per pup at virgin
        // Pups_0=SPB_virgin_adj;  //  total population fecundity is the number of pups produced
        // Sfrac=SR_parm(2);
        SRZ_0=log(1.0/(SPB_virgin_adj/Recr_virgin_adj));
        SRZ_max=SRZ_0+SR_parm(2)*(0.0-SRZ_0);
        SRZ_surv=mfexp((1.-pow((SPB_curr_adj/SPB_virgin_adj),SR_parm(3)) )*(SRZ_max-SRZ_0)+SRZ_0);  //  survival
//        warning<<"survivalcalcs "<<SRZ_0<<" "<<SRZ_max<<" "<<SRZ_surv<<" "<<SPB_curr_adj<<" "<<SPB_virgin_adj<<" ";
        NewRecruits=SPB_curr_adj*SRZ_surv;
        exp_rec(y,1)=NewRecruits;   // expected arithmetic mean recruitment
//  ProgLabel_43.3.7.1  Do variation in recruitment by adjusting survival
        if(SR_env_target==1) SRZ_surv*=mfexp(SR_parm(N_SRparm2-2)* env_data(y,SR_env_link));   // environ effect on survival
        if(recdev_cycle>0)
        {
          gg=y - (styr+(int((y-styr)/recdev_cycle))*recdev_cycle)+1;
          SRZ_surv*=mfexp(recdev_cycle_parm(gg));
        }
        exp_rec(y,2)=SPB_curr_adj*SRZ_surv;
        SRZ_surv*=mfexp(-biasadj(y)*half_sigmaRsq);     // bias adjustment
        exp_rec(y,3)=SPB_curr_adj*SRZ_surv;
        if(y <=recdev_end)
        {
          if(recdev_doit(y)>0) SRZ_surv*=mfexp(recdev(y));  //  recruitment deviation
        }
        else if(Do_Forecast>0)
        {
          SRZ_surv *= mfexp(Fcast_recruitments(y));
        }
        join=1./(1.+mfexp(100*(SRZ_surv-1.)));
        SRZ_surv=SRZ_surv*join + (1.-join)*1.0;
//        warning<<" adjust "<<SRZ_surv<<endl;
        NewRecruits=SPB_curr_adj*SRZ_surv;
        exp_rec(y,4) = NewRecruits;
        break;
      }

//  SS_Label_43.3.8  Shepard
      case 8:  // Shepard 3-parameter SRR.  per Punt document at PFMC
      {
        dvariable Shepard_c;
        dvariable Shepard_c2;
        dvariable Hupper;
        dvariable steep2;
        Shepard_c=SR_parm(3);
        Shepard_c2=pow(0.2,Shepard_c);
        Hupper=1.0/(5.0*Shepard_c2);
        steep2=0.2+(steepness-0.2)/(0.8)*(Hupper-0.2);
        temp=(SPB_curr_adj)/(SPB_virgin_adj);
        NewRecruits =  (5.*steep2*Recr_virgin_adj*(1.-Shepard_c2)*temp) /
        (1.0 - 5.0*steep2*Shepard_c2 + (5.*steep2-1.)*pow(temp,Shepard_c));
        break;
      }
    }

    if(SR_fxn!=7)
    {
//  ProgLabel_43.4  For non-survival based SRR, get recruitment deviations by adjusting recruitment itself
      exp_rec(y,1)=NewRecruits;   // expected arithmetic mean recruitment

      if(SR_env_target==1) NewRecruits*=mfexp(SR_parm(N_SRparm2-2)* env_data(y,SR_env_link));   // environ effect on annual recruitment
      if(recdev_cycle>0)
      {
        gg=y - (styr+(int((y-styr)/recdev_cycle))*recdev_cycle)+1;
        NewRecruits*=mfexp(recdev_cycle_parm(gg));
      }
      exp_rec(y,2)=NewRecruits;
      if(SR_fxn!=4) NewRecruits*=mfexp(-biasadj(y)*half_sigmaRsq);     // bias adjustment
      exp_rec(y,3)=NewRecruits;
      if(y <=recdev_end)
      {
        if(recdev_doit(y)>0) NewRecruits*=mfexp(recdev(y));  //  recruitment deviation
      }
      else if(Do_Forecast>0)
      {
        NewRecruits *= mfexp(Fcast_recruitments(y));
      }
      exp_rec(y,4) = NewRecruits;
    }
    return NewRecruits;
  }

//*********************************************************************
 /*  ProgLabel_FUNCTION 44 Equil_Spawn_Recr_Fxn */
FUNCTION dvar_vector Equil_Spawn_Recr_Fxn()
  {
    dvar_vector EquilCalc(1,2);
    dvariable B_equil;
    dvariable R_equil;
    dvariable temp;
    dvariable join;
    dvariable steepness;
    dvariable SRZ_0;
    dvariable SRZ_max;
    dvariable SRZ_surv;
    steepness=SR_parm(2);
//  ProgLabel_44.1  calc equilibrium SpawnBio and Recruitment from input SPR_temp, which is spawning biomass per recruit at some given F level
    switch(SR_fxn)
    {
      case 1: // previous placement for B-H constrained
      {
        N_warn++; cout<<"Critical error:  see warning"<<endl; warning<<"B-H constrained curve is now Spawn-Recr option #6"<<endl; exit(1);
        break;
      }
//  ProgLabel_44.1.1  Beverton-Holt with flattop beyond Bzero
      case 6: //Beverton-Holt
      {
      alpha = 4.0 * SR_parm(2)*Recr_virgin / (5.*SR_parm(2)-1.);
      beta = (SPB_virgin*(1.-SR_parm(2))) / (5.*SR_parm(2)-1.);
        B_equil=alpha * SPR_temp - beta;
        B_equil=posfun(B_equil,0.0001,temp);
        R_equil=(4.*SR_parm(2)*Recr_virgin*B_equil) / (SPB_virgin*(1.-SR_parm(2))+(5.*SR_parm(2)-1.)*B_equil);
        break;
      }
//  ProgLabel_44.1.2  Ricker
      case 2: // Ricker
      {
        B_equil=SPB_virgin*(1.+(log(Recr_virgin/SPB_virgin)+log(SPR_temp))/SR_parm(2));
        R_equil=Recr_virgin*B_equil/SPB_virgin * mfexp(SR_parm(2)*(1.-B_equil/SPB_virgin));
        break;
      }
//  ProgLabel_44.1.3  Beverton-Holt
      case 3:  // same as case 6
      {
      alpha = 4.0 * SR_parm(2)*Recr_virgin / (5.*SR_parm(2)-1.);
      beta = (SPB_virgin*(1.-SR_parm(2))) / (5.*SR_parm(2)-1.);
        B_equil=alpha * SPR_temp - beta;
        B_equil=posfun(B_equil,0.0001,temp);
        R_equil=(4.*SR_parm(2)*Recr_virgin*B_equil) / (SPB_virgin*(1.-SR_parm(2))+(5.*SR_parm(2)-1.)*B_equil); //Beverton-Holt
        break;
      }

//  ProgLabel_44.1.4  constant recruitment
      case 4: // constant; no bias correction
      {
        B_equil=SPR_temp*Recr_virgin;  R_equil=Recr_virgin;
        break;
      }
//  ProgLabel_44.1.5  Hockey Stick
      case 5: // hockey stick
      {
        alpha=SR_parm(3)*Recr_virgin;  // min recruitment level
//        temp=SPB_virgin/R0*steepness;  // spawners per recruit at inflection
        beta=(Recr_virgin-alpha)/(steepness*SPB_virgin);   //  slope of recruitment on spawners below the inflection
        B_equil=Join_Fxn(0.0*SPB_virgin/Recr_virgin, SPB_virgin/Recr_virgin, SPB_virgin/Recr_virgin*steepness, SPR_temp, alpha/((1./SPR_temp)-beta), SPR_temp*Recr_virgin);
        R_equil=Join_Fxn(0.0*SPB_virgin, SPB_virgin, SPB_virgin*steepness, B_equil, alpha+beta*B_equil, Recr_virgin);
        break;
      }
  //  SS_Label_Info_44.1.6  #3 parameter survival based
      case 7:  // survival
      {
        SRZ_0=log(1.0/(SPB_virgin/Recr_virgin));
        SRZ_max=SRZ_0+SR_parm(2)*(0.0-SRZ_0);
        B_equil = SPB_virgin * (1. - (log(1./SPR_temp) - SRZ_0)/pow((SRZ_max - SRZ_0),(1./SR_parm(3)) ));
        SRZ_surv=mfexp((1.-pow((B_equil/SPB_virgin),SR_parm(3)) )*(SRZ_max-SRZ_0)+SRZ_0);  //  survival
        R_equil=B_equil*SRZ_surv;
        break;
      }
//  SS_Label_44.1.8  3 parameter Shepard
      case 8:  // Shepard
      {
        dvariable Shepard_c;
        dvariable Shepard_c2;
        dvariable Shep_top;
        dvariable Shep_bot;
        dvariable Hupper;
        dvariable steep2;
//  Andre's FORTRAN
//        TOP = 5*Steep*(1-0.2**POWER)*SPR/SPRF0-(1-5*Steep*0.2**POWER)
//      BOT = (5*Steep-1)
//       REC = (TOP/BOT)**(1.0/POWER)*SPRF0/SPR
// Power = exp(logC);
// Hupper = 1.0/(5.0 * pow(0.2,Power));
        Shepard_c=SR_parm(3);
        Shepard_c2=pow(0.2,Shepard_c);
        Hupper=1.0/(5.0*Shepard_c2);
        steep2=0.2+(steepness-0.2)/(0.8)*(Hupper-0.2);
        Shep_top=5.0*steep2*(1.0-Shepard_c2)*(SPR_temp*Recr_virgin)/SPB_virgin-(1.0-5.0*steep2*Shepard_c2);
        Shep_bot=5.0*steep2-1.0;
        dvariable Shep_top2;
        Shep_top2=posfun(Shep_top,0.001,temp);  
        R_equil=(SPB_virgin/SPR_temp) * pow((Shep_top2/Shep_bot),(1.0/Shepard_c));
        B_equil=R_equil*SPR_temp;
        break;
      }


    }
    EquilCalc(1)=B_equil;
    EquilCalc(2)=R_equil;
    return EquilCalc;
  }

FUNCTION void get_age_age(const int Keynum, const int AgeKey_StartAge, const int AgeKey_Linear1, const int AgeKey_Linear2)
  //*********************************************************************
  //  SS_Label_Function_45 #get_age_age 
  {
    //  FUTURE: calculate adjustment to oldest age based on continued ageing of old fish
    age_age(Keynum).initialize();
    dvariable age;
    dvar_vector age_err_parm(1,7);
    dvariable temp;

  if(Keynum==Use_AgeKeyZero)
  {
  //  SS_Label_Info_45.1 set age_err_parm to mgp_adj, so can be time-varying according to MGparm options
    for(a=1;a<=7;a++)
    {age_err_parm(a)=mgp_adj(AgeKeyParm-1+a);}
      age_err(Use_AgeKeyZero,1)(0,AgeKey_StartAge)=r_ages(0,AgeKey_StartAge)+0.5;
      age_err(Use_AgeKeyZero,2)(0,AgeKey_StartAge)=age_err_parm(5)*(r_ages(0,AgeKey_StartAge)+0.5)/(age_err_parm(1)+0.5);
  //  SS_Label_Info_45.3 calc ageing bias
      if(AgeKey_Linear1==0)
      {
        age_err(Use_AgeKeyZero,1)(AgeKey_StartAge,nages)=0.5 + r_ages(AgeKey_StartAge,nages) + age_err_parm(2)+(age_err_parm(3)-age_err_parm(2))*
        (1.0-mfexp(-age_err_parm(4)*(r_ages(AgeKey_StartAge,nages)-age_err_parm(1)))) / (1.0-mfexp(-age_err_parm(4)*(r_ages(nages)-age_err_parm(1))));
      }
      else
      {
        age_err(Use_AgeKeyZero,1)(AgeKey_StartAge,nages)=0.5 + r_ages(AgeKey_StartAge,nages) + age_err_parm(2)+(age_err_parm(3)-age_err_parm(2))*
        (r_ages(AgeKey_StartAge,nages)-age_err_parm(1))/(r_ages(nages)-age_err_parm(1));
      }
  //  SS_Label_Info_45.4 calc ageing variance
      if(AgeKey_Linear2==0)
      {
        age_err(Use_AgeKeyZero,2)(AgeKey_StartAge,nages)=age_err_parm(5)+(age_err_parm(6)-age_err_parm(5))*
        (1.0-mfexp(-age_err_parm(7)*(r_ages(AgeKey_StartAge,nages)-age_err_parm(1)))) / (1.0-mfexp(-age_err_parm(7)*(r_ages(nages)-age_err_parm(1))));
      }
      else
      {
        age_err(Use_AgeKeyZero,2)(AgeKey_StartAge,nages)=age_err_parm(5)+(age_err_parm(6)-age_err_parm(5))*
        (r_ages(AgeKey_StartAge,nages)-age_err_parm(1))/(r_ages(nages)-age_err_parm(1));
      }
  }

  /*  SS_Label_Info_45.5 calc distribution of age' for each age */
   for (a=0; a<=nages;a++)
    {
     if(age_err(Keynum,1,a)<=-1)
       {age_err(Keynum,1,a)=r_ages(a)+0.5;}
     age=age_err(Keynum,1,a);

     for (b=2;b<=n_abins;b++)     //  so the lower tail is accumulated into the first age' bin
     age_age(Keynum,b,a)= cumd_norm((age_bins(b)-age)/age_err(Keynum,2,a));

     for (b=1;b<=n_abins-1;b++)
       age_age(Keynum,b,a) = age_age(Keynum,b+1,a)-age_age(Keynum,b,a);

     age_age(Keynum,n_abins,a) = 1.-age_age(Keynum,n_abins,a) ;     // so remainder is accumulated into the last age' bin

    }

     if(gender == 2)                     //  copy ageing error matrix into male location also
     {
      L2=n_abins;
      A2=nages+1;
      for(b=1;b<=n_abins;b++)
      for(a=0;a<=nages;a++)
       {age_age(Keynum,b+L2,a+A2)=age_age(Keynum,b,a);}
     }
  }
