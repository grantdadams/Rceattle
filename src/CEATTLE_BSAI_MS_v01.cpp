#include <TMB.hpp>
// ------------------------------------------------------------------------- //
//                 CEATTLE version 3.1.2                                     //
//                  Template Model Builder                                   //
//               Multispecies Statistical Model                              //
//          Bioenergetic-based Assessment for Understanding                  //
//              Biomass Linkages To The Environment                          //
//                  in the Bering Sea                                        //
//                   Oct. 2017                                               //
//                                                                           //
// AUTHORS:   Kirstin Holsman, Jim Ianelli                                   //
//            Modified by Grant Adams                                        //
// CITATIONS:                                                                //
// 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
// 0. LOAD DEPENDENCIES                                                      //
// ------------------------------------------------------------------------- //
#include "../src/include/functions.hpp"


template<class Type>
Type objective_function<Type>::operator() (){
  // ------------------------------------------------------------------------- //
  // 1. MODEL CONFIGURATION                                                    //
  // ------------------------------------------------------------------------- //
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(debug);            // Logical vector to debug or not
  DATA_INTEGER(msmMode);
  //    0 = run in single species mode
  //    1 = run in MSM mode

  DATA_IVECTOR(logist_sel_phase); // Selectivity for BT survey
  //    0 = fit to data
  //    1 = logistic

  // 1.2. Temporal dimensions
  DATA_INTEGER( nyrs );           // Number of estimation years
  DATA_INTEGER( styr );           // Start year

  // 1.3. Number of species
  DATA_INTEGER( nspp );           // Number of species (prey)


  // 1.4. MODEL OBJECTS
  // 1.4.1. LOOPING INDICES -- k = observation, i = species/prey, j = age/prey age (yr), y = year, p = predator, a = predator age (yr)
  int  i, j, y, k, p, a;
  int fsh_yr_ind;

  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                           //
  // ------------------------------------------------------------------------- //

  // 2.1. FIXED VALUES
  int tau = 200; // Fishery age composition sample size
  Type MNConst = 0.001;           // Constant additive for logistic functions
  int nselages = 8;               // Number of selectivity ages
  Type curv_pen_fsh = 12.5;       // Fishery selectivity penalty
  Type sigma_catch = 0.05;        // SD of catch


  // 2.2. DIMENSIONS OF DATA
  int endyr = nyrs + styr;        // End year

  // -- 2.2.2. Species attributes
  DATA_IVECTOR( nages );          // Number of species (prey) ages
  int max_age = imax(nages);      // Integer of maximum nages to make the arrays.
  // DATA_INTEGER( n_pred);       // Number of predator species
  // DATA_IVECTOR( n_age_pred);   // Number of predator ages


  // 2.3. DATA INPUTS (i.e. assign data to objects)
  // -- 2.3.1 Fishery Components
  DATA_IVECTOR( nyrs_tc_biom );   // Number of years with total observed catch; n = [nspp]
  DATA_IMATRIX( yrs_tc_biom );    // Years with total observed catch; n = [nspp, nyrs_tc_biom]
  DATA_MATRIX( tcb_obs );         // Observed total yield (kg); n = [nspp, nyrs_tc_biom]

  DATA_IVECTOR( nyrs_fsh_comp );  // Number of years in the fishery sp_age composition data; n = [nspp]
  DATA_IMATRIX( yrs_fsh_comp );   // Years for the fishery sp_age composition data; n = [nspp, nyrs_fsh_comp]
  DATA_IVECTOR( fsh_age_type );   // Which method of calculating fishery age hat (2 = ATF); n = [nspp]
  DATA_IVECTOR( fsh_age_bins );   // Bins for fishery age composition data; n = [nspp]
  DATA_ARRAY( obs_catch );        // Observed fishery catch-at-age or catch-at-length; n = [nspp, fsh_age_bins, nyrs_fsh_comp]

  // -- 2.3.2 BT Survey Components
  DATA_IVECTOR( nyrs_srv_biom );  // Number of years of survey biomass data; n = [nspp]
  DATA_IMATRIX( yrs_srv_biom );   // Years of survey biomass data; n = [nspp, nyrs_srv_biom]
  DATA_MATRIX( srv_biom );        // Observed BT survey biomass (kg); n = [nspp, nyrs]
  DATA_MATRIX( srv_biom_se );     // Observed annual biomass error (SE); n = [nspp, nyrs_srv_biom]

  // -- 2.3.3. BT Survey age components
  DATA_IVECTOR( nyrs_srv_age );   // Number of years of survey age/length composition; n = [1, nspp]
  DATA_IMATRIX( yrs_srv_age );    // Years for the survey age/length composition data; n = [nspp, nyrs_srv_age]
  DATA_IVECTOR( srv_age_type );   // Type of compisition (1 = age; 2 = length); n = [1, nspp]
  DATA_IVECTOR( srv_age_bins );   // Number of size binds for the age/length comps; n = [nspp]
  DATA_MATRIX( srv_age_n );       // Sample size for the multinomial; n = [nspp, nyrs_srv_age]
  DATA_ARRAY( srv_age_obs ) ;     // Observed BT age comp; n = [nspp, nages, nyrs]
  DATA_MATRIX( srv_age_sizes );   // Observed size composition
  DATA_ARRAY( age_trans_matrix);  // observed sp_age/size compositions; n = [nspp, nages, srv_age_bins]

  // -- 2.3.4 EIT Survey components
  DATA_INTEGER( n_eit);           // Number of years with EIT data; n = [1]
  DATA_IVECTOR( yrs_eit);         // Years for available EIT data; n = [n_eit]
  DATA_VECTOR( eit_age_n);        // Number  of  EIT Hauls for multinomial; n = [yrs_eit]
  DATA_MATRIX( obs_eit_age);      // Observed EIT catch-at-age; n = [n_eit, nyrs]
  DATA_VECTOR( obs_eit);          // Observed EIT survey biomass (kg); n = [1, nyrs]
  DATA_MATRIX( eit_sel);          // Observed EIT survey selectivity; n = [eit_age, nyrs_eit_sel]

  // -- 2.3.5. Weight-at-age
  DATA_IVECTOR( nyrs_wt_at_age ); // Number of years of weight at age data; n = [nspp]
  DATA_IMATRIX( yrs_wt_at_age );  // Years of weight-at-age; data n = [nspp, nyrs_wt_at_age]
  DATA_ARRAY( wt );               // Weight-at-age by year; n = [nyrs_wt_at_age, nages, nspp]

  // -- 2.3.6. Q Diet data
  DATA_INTEGER( npred3 );         // Number of predators; n = [1]
  DATA_IVECTOR( nprey3 );         // Number of prey species for each predator; n = [1, nspp]
  DATA_IVECTOR( n_stomyrs );      // number years of stomach data; n = [1, nspp]
  DATA_IMATRIX( stomyrs );        // years of stomach data; n = [nspp, n_stomyrs]
  DATA_ARRAY( mnQ );              // meanwt of each prey spp in the stomach of each predator of sp_age a; n = [npred3,n_stomyrs,max_age,nspp+1]
  DATA_ARRAY( Qse );              // SE wt of each prey spp in the stomach of each predator of sp_age a; n = [npred3,n_stomyrs,max_age,nspp+1]


  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( other_food );      // Biomass of other prey (kg); n = [1, nspp]
  DATA_IVECTOR( useWt );          // Assign relative proportion of prey in the diet according to relative biomass in the system.,otherwise the model with use relative proportion by number; n = [1, nspp]
  DATA_IVECTOR( C_model );        // f == 1, the use Cmax*fT*P; n = [1, nspp]
  DATA_VECTOR( Pvalue );          // This scales the pvalue used if C_model ==1 , proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge; n = [1, nspp]
  DATA_IVECTOR( Ceq );            // Ceq: which Comsumption equation to use; n = [1, nspp]
  DATA_VECTOR( CA );              // Wt specific intercept of Cmax=CA*W^CB; n = [1, nspp]
  DATA_VECTOR( CB );              // Wt specific slope of Cmax=CA*W^CB; n = [1, nspp]
  DATA_VECTOR( Qc );              // used in fT, QC value; n = [1, nspp]
  DATA_VECTOR( Tco );             // used in fT, thermal optimum; n = [1, nspp]
  DATA_VECTOR( Tcm );             // used in fT, thermal max; n = [1, nspp]
  DATA_VECTOR( Tcl );             // used in fT eq 3, limit; n = [1, nspp]
  DATA_VECTOR( CK1 );             // used in fT eq 3, limit where C is .98 max (ascending); n = [1, nspp]
  DATA_VECTOR( CK4 );             // used in fT eq 3, temp where C is .98 max (descending); n = [1, nspp]
  DATA_MATRIX( S_a );             // S_a, S_b, S_b2, S_b3, S_b4, S_b5: a,L,L^2,L^3,L^4,L^5 (rows)coef for mean S=a+b*L+b2*L*L, whith a cap at 80cm for each pred spp(cols); n = [6, nspp]
  DATA_MATRIX( aLim );            // #aLim and bLim : upper limit of prey size; n = [2, nspp]


  // -- 2.4.2. von Bertalannfy growth function (VBGF)
  DATA_VECTOR( t0 );              // t0 parameter of the temp specific VonB for wt; n = [nspp, mf_type]
  DATA_VECTOR( log_mean_d );      // log mean d parameter of the temp specific von B for wt; n = [nspp, mf_type]
  DATA_VECTOR( logK );            // log k parameter of the temp specific VonB for wt; n = [nspp, mf_type]
  DATA_VECTOR( logH );            // log H parameter of the temp specific VonB for wt; n = [nspp, mf_type]
  DATA_VECTOR( Tcoef );           // T coefficent of the linear d equations of the temp specific VonB for wt; n = [nspp, mf_type]
  DATA_VECTOR( Pcoef );           // P-value coefficent of the linear d equations of the temp specific VonB for wt; n = [nspp, mf_type]


  // -- 2.4.3. Others
  DATA_MATRIX( M1_base );         // Residual natural mortality; n = [nspp, nages]
  DATA_IVECTOR( mf_type );        // Sex specific mort and weight at age? : 1 = same for both, 2 = seperate wt at sp_age for each sex
  DATA_MATRIX( propMorF );        // Proportion-at-age of females of population; n = [nspp, nages]
  DATA_MATRIX( pmature );         // Proportion of mature females at age; [nspp, nages]


  // -- 2.4.4. Weight-at-length parameters
  DATA_MATRIX( aLW );             // LW a&b regression coefs for W=a*L^b; n = [2, nspp]

  // 2.5. DIET DATA
  // -- 2.5.1. Diet proportion parameters
  DATA_MATRIX( maxK );            // Matrix of maximum diet proportions for each predd,preyy combo (used for broken stick expon); n = [nspp, nspp]
  DATA_VECTOR( klim );            // Max  pred  length  at  which K should  be  calculated  (beyond this  is  it  just  equal to  K); n = [1, nspp]
  DATA_MATRIX( prefa );           // prefa of the functional response prefa(pred,prey); n = [nspp, nspp]
  DATA_MATRIX( prefb );           // prefb of the functional response prefa(pred,prey); n = [nspp, nspp]
  DATA_MATRIX( l1 );              // l1 of the switch function l1(pred,prey); n = [nspp, nspp]
  DATA_MATRIX( l2 );              // l2 of the switch function l2(pred,prey); n = [nspp, nspp]
  DATA_MATRIX( LA_Lengths );      // length to sp_age coversion matrix (rounded to nearest integer); n = [nspp, nages]
  DATA_VECTOR( fday );            // number of foraging days for each predator; n = [1, nspp] #FIXME - assuming this is the same as fdays
  DATA_IVECTOR( nlengths );       // number of Lengths for the matrix to convert ages to lengths; n = [1, nspp]
  int maxL = imax( nlengths );    // Maximum number of lengths for the matrix to convert between ages and lengths; n = [1]
  DATA_MATRIX( lengths );         // Lengths for the matrix to convert ages to lenghths; n = [nspp, nlengths]
  // DATA_ARRAY( A2L_matrix );       // A2L_matrix : Matrix to convert ages to lenghths; n = [nspp, nages, nlengths]
  DATA_ARRAY( K );                // K(1) is nprey by maxL matrix of stomach proportions of predator 1; n = [nspp, nspp, nlengths]
  DATA_ARRAY( KAge );             // K(1) is nprey by maxL matrix of stomach proportions of predator 1; n = [nspp, nspp, nages]
  DATA_MATRIX( PAge );            // n = [nspp, nages]
  DATA_ARRAY( Pyrs );             // n = [nspp, nyrs+1, nages]: #FIXME - Assuming this is the same as Pby_yr?
  DATA_ARRAY( Uobs );             // pred, prey, predL, preyL U matrix (mean number of prey in each pred); n = [nspp, nspp, maxL, maxL]
  DATA_ARRAY( UobsWt );           // pred, prey, predL, preyL U matrix (mean wt_hat of prey in each pred); n = [nspp, nspp, maxL, maxL] #FIXME - Changed name in stomach2017.dat
  DATA_ARRAY( UobsAge );          // pred, prey, predA, preyA U matrix (mean number of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_ARRAY( UobsWtAge );        // pred, prey, predA, preyA U matrix (mean wt_hat of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_MATRIX( Mn_LatAge );       // Mean length-at-age; n = [nspp, nages]


  // -- 2.5.2. F Profile data
  /*
   DATA_INTEGER( n_f );            // Number of F vectors
   DATA_VECTOR( Frates );          // Fishing mortality vector; n = [1, n_f]
   DATA_INTEGER( np );             // Number of
   DATA_MATRIX( Fprofiles );       // F profiles; n = [np, nspp]
   */


  // 2.6. ENVIRONMENTAL DATA

  DATA_INTEGER( nTyrs );          // Number of temperature years; n = [1] #FIXME - changed the name of this in retro_data2017_asssmnt.dat
  DATA_IVECTOR( Tyrs );           // Years of hindcast data; n = [1, nTyrs] #FIXME - changed the name of this in retro_data2017_asssmnt.dat
  DATA_VECTOR( BTempC_retro );    // Vector of bottom temperature; n = [1,  nTyrs ]
  // DATA_INTEGER( ncov );        // Number of environmental covariates


  //
  // 2.5. DERIVED QUANTITIES # Calculate these in the model
  // DATA_MATRIX(d);               // VBGF allometric slope of consumption (d); n = [nspp, nyrs]
  // DATA_MATRIX(Winf);            // VBGF max asymptoptic weight; n = [nspp, nyrs]


  // ------------------------------------------------------------------------- //
  // 2.8. Debugging with data inputs                                           //
  // ------------------------------------------------------------------------- //
  if(debug == 1){
    // -- 2.8.2.1 Check to make sure the first year of survey data are not before start year
    for(int i = 0; i<nspp; i++){
      if(yrs_tc_biom(i,0) < styr){
        std::cerr<<"First year of total catch biomass of species "<< i + 1 << " is before specified start year"<<std::endl;
        return(0);
      }
    }

    // -- 2.8.2.2 Check to make sure the first year of survey data are not before start year
    for(int i = 0; i<nspp; i++){
      if(yrs_srv_biom(i,0) < styr){
        std::cerr<<"First year of survey biomass of species "<< i + 1 << " is before specified start year"<<std::endl;
        return(0);
      }
    }

    // -- 2.8.2.3. Check to make sure the years of survey data biomass and age are the same
    for(int i = 0; i<nspp; i++){
      if(nyrs_srv_biom(i) != nyrs_srv_age(i)){
        std::cerr<<"Nyrs of survey biomass and age-comp of species "<< i + 1 << " do not match"<<std::endl;
        return(0);
      }
    }

    if(yrs_eit(0) < styr){
      std::cerr<<"First year of EIT survey biomass is before specified start year"<<std::endl;
      return(0);
    }
  }


  // ------------------------------------------------------------------------- //
  // 3. INITIAL CALCULATIONS                                                   //
  // ------------------------------------------------------------------------- //

  int max_bin = imax(srv_age_bins); // Integer of maximum number of length/age bins.

  // 3.1. FISHERY DATA
  // 3.1.1. Define objects
  array<Type>  fsh_age_obs(nyrs, max_bin, nspp);     // Observed fishery age comp; n = [nyrs_fsh_comp, fsh_age_bins, nspp]
  matrix<Type> tc_obs(nspp, nyrs); tc_obs.setZero(); // Set total catch to 0 to initialize // Observed total catch (n); n = [nspp, nyrs] NOTE: This may not be necessary if loading data from tmp


  // 3.1.2. Fishery catch-at-age to age-comp
  for(i=0; i < nspp; i++){
    for(y = 0; y < nyrs_fsh_comp(i); y++){
      for(j=0; j < fsh_age_bins(i); j++){
        tc_obs(i, y) += obs_catch(y, j, i);
      }
      for(j=0; j < fsh_age_bins(i); j++){
        fsh_age_obs(y, j, i) = obs_catch(y, j, i)/(tc_obs(i, y) + 0.01); // Calculate age comp
      }
    }
  }

  // 3.2. SURVEY DATA
  // 3.2.1. Define objects
  matrix<Type> srv_biom_lse(nspp, nyrs); // Observed annual biomass CV; n = [nspp, nyrs_srv_biom]
  srv_biom_lse = srv_biom_se.array()/ srv_biom.array();          // CV estimation
  srv_biom_lse = pow( log( ( pow( srv_biom_lse.array(), Type(2) ).array() + 1).array()).array(), Type(0.5));
  matrix<Type>  eit_age_comp = obs_eit_age;          // Eit age comp; n = [n_eit, srv_age_bins(0)]

  // 3.2.2. EIT catch-at-age to age-comp
  eit_age_comp.setZero();
  for(y = 0; y < n_eit; y++){
    for(j = 0; j < srv_age_bins(0); j++){
      eit_age_comp(y, j) = obs_eit_age(y, j)/ obs_eit_age.row(y).sum(); // Convert from catch-at-age to age comp
    }
  }


  // 3.3. DIET DATA
  vector<int> prey_nlengths( nspp );
  vector<int> prey_nages( nspp );

  // FIXME - this seems redundant
  for(p = 0; p < nspp; p++){
    prey_nlengths( p ) = nlengths( p );
    prey_nages( p ) = nages( p );
  }


  // 3.3 POPULATION DYNAMICS
  matrix<Type> M1( nspp, max_age); M1.setZero();
  M1 = M1_base.array() + Type(0.0001);

  for ( i = 0; i < nspp ; i++){
    for( j = 0 ; j < nages(i); j++ ){
      pmature( i, j ) = pmature( i, j ) * propMorF(i + (mf_type(i) - 1), j);
    }
  }


  // ------------------------------------------------------------------------- //
  // 4. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //

  PARAMETER(dummy);                    // Variable to test derived quantities given input parameters; n = [1]

  // 4.1. PARAMETERS (assign parameters to objects)
  // -- 4.1.1 Recruitment parameters
  PARAMETER_VECTOR( ln_mn_rec );       // Mean recruitment; n = [1, nspp]
  PARAMETER_MATRIX( rec_dev );         // Annual recruitment deviation; n = [nspp, nyrs]
  // PARAMETER(sigma_rec);             // Standard deviation of recruitment variation # NOTE: Have this estimated if using random effects.
  // -- 4.1.2. Abundance parameters
  PARAMETER_MATRIX( init_dev );             // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this
  // -- 4.1.3. fishing mortality parameters
  PARAMETER_VECTOR( ln_mean_F );      // Log mean fishing mortality; n = [1, nspp]
  PARAMETER_MATRIX( F_dev );          // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change
  // PARAMETER_VECTOR(sigma_F); // SD of fishing mortality deviations; n = [1, nspp] # NOTE: Have this estimated if using random effects.
  // -- 4.1.4. Selectivity parameters
  PARAMETER_MATRIX( srv_sel_coff );    // Survey selectivity parameters; n = [nspp, nselages]
  PARAMETER_MATRIX( fsh_sel_coff );    // Fishery age selectivity coef; n = [nspp, nselages]
  PARAMETER_VECTOR( srv_sel_slp );     // Survey selectivity paramaters for logistic; n = [1, nspp]
  PARAMETER_VECTOR( srv_sel_inf );     // Survey selectivity paramaters for logistic; n = [1, nspp]
  PARAMETER( log_eit_q );              // EIT Catchability; n = [1]
  PARAMETER_VECTOR( log_srv_q );       // BT Survey catchability; n = [1, nspp]


  // 4.2. DERIVED QUANTITIES
  // -- 4.2.1. Fishery observations
  matrix<Type>  tc_biom_hat(nspp, nyrs);              // Estimated total yield (kg); n = [nspp, nyrs]
  array<Type>   catch_hat(nyrs, max_age, nspp);       // Estimated catch-at-age (n); n = [nspp, nages, nyrs]
  array<Type>   fsh_age_hat(nyrs, max_bin, nspp);     // Estimated fishery age comp; n = [nspp, nages, nyrs]
  matrix<Type>  tc_hat(nspp, nyrs);                   // Estimated total catch (n); n = [nspp, nyrs]
  array<Type>   F(nyrs, max_age, nspp);               // Estimated fishing mortality; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_sel(nspp, max_age);               // Log estimated fishing selectivity
  vector<Type>  avgsel_fsh(nspp);                     // Average fishery selectivity

  // -- 4.2.2. BT Survey components
  matrix<Type>  srv_bio_hat(nspp, nyrs);              // Estimated BT survey biomass (kg); n = [nspp, nyrs]
  array<Type>   srv_age_hat(nyrs, max_bin, nspp);     // Estimated BT age comp; n = [nspp, nages, nyrs]
  matrix<Type>  srv_hat(nspp, nyrs);                  // Estimated BT survey total abundance (n); n = [nspp, nyrs]
  matrix<Type>  srv_sel(nspp, max_age);               // Estimated survey selectivity at age; n = [nspp, nyrs]
  vector<Type>  avgsel_srv(nspp);                     // Average survey selectivity

  // -- 4.2.3. EIT Survey Components
  matrix<Type>  eit_age_hat(nyrs, srv_age_bins(0));   // Estimated EIT catch-at-age ; n = [nyrs, 12 ages]
  matrix<Type>  eit_age_comp_hat(nyrs, srv_age_bins(0));// Estimated EIT age comp ; n = [nyrs, 12 ages]
  vector<Type>  eit_hat(nyrs);                        // Estimated EIT survey biomass (kg); n = [nyrs]

  // -- 4.2.4. Estimated population parameters
  matrix<Type>  R(nspp, nyrs);                        // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type>   NByage(nyrs, max_age, nspp);          // Numbers at age; n = [nspp, nages, nyrs]
  array<Type>   AvgN(nyrs, max_age, nspp);            // Average numbers-at-age; n = [nspp, nages, nyrs]
  array<Type>   biomassByage(nyrs, max_age, nspp);    // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type>  biomass(nspp, nyrs);                  // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type>  biomassSSB(nspp, nyrs);               // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type>   biomassSSBByage(nyrs, max_age, nspp); // Spawning biomass at age (kg); n = [nspp, nages, nyrs]

  // -- 4.3. Parameter transformations
  Type eit_q = exp(log_eit_q);                        // EIT Catchability
  array<Type>   Zed(nyrs, max_age, nspp);             // Total mortality at age; n = [nspp, nages, nyrs]
  array<Type>   S(nyrs, max_age, nspp);               // Survival at age; n = [nspp, nages, nyrs]

  // -- 4.4. Ration components
  array<Type> ConsumAge( nyrs, max_age, nspp );       // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
  array<Type> Consum_livingAge( nyrs, max_age, nspp );// Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
  array<Type> S2Age( nyrs, max_age, nspp );           // pre-allocate mean stomach weight as a function of sp_age
  array<Type> LbyAge( nyrs, max_age, nspp );          // Length by age from LW regression
  matrix<Type> mnWt_obs( nspp, max_age );             // Mean observed weight at age (across years); n = [nspp, nages]
  matrix<Type> fT( nspp, nTyrs );                     // Pre-allocation of temperature function of consumption; n = [nspp, nTyrs]
  vector<Type> TempC( nTyrs );                        // Bottom temperature; n = [1, nTyrs]
  array<Type> ration2Age( nyrs, max_age, nspp );      // Annual ration at age (kg/yr); n = [nyrs, nages, nspp]

  // -- 4.5. Suitability components
  array<Type> suit_main(nspp, nspp, max_age, max_age);// Suitability; n = [nspp, nspp, nages, nages]
  matrix<Type> suit_other(nspp, max_age);             // Suitability not accounted for by the included prey; n = [nspp, nages]
  array<Type> stom_div_bio2(nspp, nspp, max_age, max_age, nyrs);// Stomach proportion over biomass; U/ (W * N); n = [nspp, nspp, nages, nages, nyrs]
  stom_div_bio2.setZero();
  array<Type> stomKir(nspp, nspp, max_age, max_age, nyrs);// Stomach proportion U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type> avail_food(nyrs, max_age, nspp);            // Available food to predator; n = [nyrs, nages, nspp]
  array<Type> of_stomKir(nyrs, max_age, nspp);            // Other food stomach content; n = [nyrs, nages, nspp] # FIXME - what is this?
  array<Type> M2(nyrs, max_age, nspp);                    // Predation mortality at age; n = [nyrs, nages, nspp]
  array<Type> B_eaten(nyrs, max_age, nspp);               // Biomass of prey eaten via predation; n = [nyrs, nages, nspp]

  // Calculate temperature to use
  TempC = BTempC_retro.sum()/nTyrs; // Fill with average bottom temperature

  int yr_ind = 0;
  for(y=0; y < nTyrs; y++){
    yr_ind = Tyrs( y ) - styr;
    if((yr_ind >= 0) & (yr_ind < nyrs)){
      TempC(yr_ind) = BTempC_retro( y );
    }
  }


  // Calculate length-at-age
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        LbyAge( y, j, i) = ( pow( ( 1 / aLW(0, i) ), (1 / aLW(1, i) ) ) ) * pow( ( wt(y, j, i)*1000),(1 / aLW(1, i))); // W = a L ^ b is the same as (W/a)^(1/b)
      }
    }
  }

  // Calculate stomach weight by sp age
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        S2Age(y, j, i) = S_a(0, i) + (S_a(1, i) * LbyAge(y, j, i)) + (S_a(2, i) * pow(LbyAge(y, j, i), 2)) + (S_a(3, i) * pow(LbyAge(y, j, i), 3)) + (S_a(4, i) * pow(LbyAge(y, j, i), 4)) + (S_a(5, i) * pow(LbyAge(y, j, i), 5));
        if(LbyAge(y, j, i)>80){
          S2Age(y, j, i) = S_a(0, i) + (S_a(1, i) * 80) + (S_a(2, i) * pow(80, 2)) + (S_a(3, i) * pow(80, 3)) + (S_a(4, i) * pow(80, 4)) + (S_a(5, i) * pow(80, 5)); // set everything above 80 to 80
        }
      }
    }
  }


  // Calculate temperature function of consumption
  fT.setZero();
  Type Yc = 0;
  Type Zc = 0;
  Type Vc = 0;
  Type Xc = 0;
  Type G2 = 0;
  Type L2 = 0;
  Type G1 = 0;
  Type L1 = 0;
  Type Ka = 0;
  Type Kb = 0;
  for(i=0; i < nspp; i++){
    for(y=0; y < nTyrs; y++){
      if( Ceq(i) == 1){
        fT(i, y) = exp(Qc(i) * TempC(y));
      }
      if( Ceq(i) == 2){
        Yc = log( Qc(i) ) * (Tcm(i) - Tco(i) + 2);
        Zc = log( Qc(i) ) * (Tcm(i) - Tco(i));
        Vc = (Tcm(i) - TempC(y)) / (Tcm(i) - Tco(i));
        Xc = pow(Zc, 2) * pow((1 + pow((1 + 40 / Yc), 0.5)), 2) / 400;
        fT(i, y) = pow(Vc, Xc) * exp(Xc * (1 - Vc));
      }
      if (Ceq(i) == 3){
        G2 = (1 / (Tcl(i) - Tcm(i))) * log((0.98*(1 - CK4(i)))/(CK4(i) * 0.02));
        L2 = exp(G2 * (Tcl( i ) - TempC( y )));
        Kb = (CK4(i) * L2)/(1 + CK4(i) * (L2-1));
        G1 = (1 / (Tco(i) - Qc(i))) * log((0.98 * (1 - CK1(i)))/(CK1(i) * 0.02));
        L1 = exp(G1 * (TempC(y) - Qc(i)));
        Ka = (CK1(i) * L1) / (1 + CK1(i) * (L1 - 1));
        fT(i, y) = Ka * Kb;
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 5. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0

  // 8.1. ESTIMATE FISHERY SELECTIVITY
  fsh_sel.setZero();
  avgsel_fsh.setZero();
  for (i = 0; i <nspp; i++){
    for(j=0; j < nselages; j++){
      fsh_sel(i, j) = fsh_sel_coff(i, j);
      avgsel_fsh(i) +=  exp(fsh_sel_coff(i, j));
    }
    // 8.1.3 Average selectivity up to nselages
    avgsel_fsh(i) = log(avgsel_fsh(i) / nselages);

    // 8.1.4. Plus group selectivity
    for(j = nselages; j < nages(i); j++){
      fsh_sel(i, j) = fsh_sel(i, nselages - 1);
    }

    // 8.1.5. Average selectivity across all ages
    Type avgsel_tmp = 0; // Temporary object for average selectivity across all ages
    for(j = 0; j < nages(i); j++){
      avgsel_tmp += exp(fsh_sel(i, j));
    }
    avgsel_tmp = log(avgsel_tmp / nages(i));

    // 8.1.6. Standardize selectivity
    for(j = 0; j < nages(i); j++){
      fsh_sel(i, j) -= avgsel_tmp;
      fsh_sel(i, j) = exp(fsh_sel(i, j));
    }
  }


  // 8.2. ESTIMATE FISHING MORTALITY
  for(i=0; i<nspp; i++){
    for(y=0; y<nyrs; y++){
      for(j=0; j<nages(i); j++){
        F(y, j, i) = fsh_sel(i, j) * exp(ln_mean_F(i) + F_dev(i, y));
      }
    }
  }


  // 8.4. Estimate total mortality at age NOTE: May need to go above population dynamics
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){
        Zed(y, j, i) = M1(i, j) + F(y, j, i) + M2(y, j, i);

        // Add predation mortality
        if(msmMode > 0){
          Zed(y, j, i)  += M2(y, j, i);
        }

        S(y, j, i) = exp(-Zed(y, j, i));
      }
    }
  }


  // 5.1. ESTIMATE RECRUITMENT T1.1
  for(i=0; i< nspp; i++){
    for(y=0; y < nyrs; y++){
      R(i,y) = exp(ln_mn_rec(i) + rec_dev(i,y));
      NByage(y, 0, i) = R(i, y);
    }
    //NByage(nyrs, 0, i) = R(i, nyrs-1);
  }


  // 5.2. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      if((j > 0) & (j < nages(i) - 1)){
        NByage(0, j, i) = exp(ln_mn_rec(i) - (j)*M1(i,j) + init_dev(i, j-1));
      }
      // -- 5.2.2. Where y = 1 and j > Ai.
      if(j == (nages(i)-1)){
        NByage(0, j, i) = exp(ln_mn_rec(i) - (j)*M1(i,j) + init_dev(i, j-1))/ (1-exp(-M1(i, nages(i)-1))); // NOTE: This solves for the geometric series
      }
    }
  }


  // 5.3. ESTIMATE NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
  biomass.setZero();  // Initialize Biomass
  biomassSSB.setZero(); // Initialize SSB
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=1; y < nyrs; y++){
        // -- 5.3.1.  Where 1 <= j < Ai
        if(j < (nages(i)-1)){
          NByage(y, j+1, i) = NByage(y-1, j, i) * S(y-1, j, i);
        }

        // -- 5.3.2. Plus group where j > Ai. NOTE: This is not the same as T1.3 because I used j = A_i rather than j > A_i.
        if(j == (nages(i)-1)){
          NByage(y, nages(i)-1, i) = NByage(y-1, nages(i)-2, i) * S(y-1, nages(i)-2, i) + NByage(y-1, nages(i)-1, i) * S(y-1, nages(i)-1, i);
        }
      }

      // -- 5.3.3. Estimate Biomass and SSB
      for(y=0; y < nyrs; y++){
        // NOTE: The "-1" is because the NByage has the initial population numbers
        biomassByage(y, j, i) = NByage(y, j, i) * wt(y, j, i); // 5.5.
        biomassSSBByage(y, j, i) = biomassByage(y, j, i) * pmature(i, j); // 5.6.

        biomass(i, y) += biomassByage(y, j, i);
        biomassSSB(i, y) += biomassSSBByage(y, j, i);
      }
    }
  }


  // 5.4. ESTIMATE AVERAGE NUMBERS AT AGE
  for(i = 0; i < nspp; i++){
    for(j = 0; j < nages(i); j++){
      for(y = 0; y < nyrs; y++){
        AvgN(y, j, i) = NByage(y, j, i) * (1 - S(y, j, i)) / Zed(y, j, i);
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 6. PREDATION MORTALITY EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE -- LOOPING INDICES -- k = observation, i = species/prey, j = age/prey age, y = year, p = predator, a = predator age

  // 6.1. Calculate historic ration
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){ // Caclulate ration for each species
        ConsumAge(y, j, i) = Type(24) * Type(0.0134) * exp( Type(0.0115) * TempC( y )) * Type(91.25) * S2Age(y, j, i) * wt(y, j, i); // Calculate consumption for predator-at-age; units = kg/predator
        Consum_livingAge(y, j, i) = ConsumAge(y, j, i); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        if(C_model( i )==1){
          ConsumAge(y, j, i) = CA(i) * pow(wt(y, j, i) * Type(1000), CB( i )) * fT(i, y) * fday( i ) * wt(y, j, i) * 1000;//g/pred.yr
          ConsumAge(y, j, i) = ConsumAge(y, j, i) * Pvalue(i) * Pyrs(y, j, i); //
        }
        ration2Age(y, j, i) = ConsumAge(y, j, i) / 1000; //annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
      }
    }
  }

  // 6.2. Calculate stomach content
  for (y=0; y<nyrs; y++){           // Year loop
    for(p=0; p < nspp; p++){                  // Predator species loop
      for (a=0; a < nages(p); a++){           // Predator age loop
        for (i=0; i < nspp; i++){             // Prey species loop
          for (j=0; j < nages(i); j++){       // Prey age loop
            stomKir(p, i, a, j, y) = UobsAge(p , i , a, j); //  FIXEME - I think this is equivalent to stomKir(yr,pred,pred_age,prey)=UobsAge(pred,prey,pred_age);
          }
        }
      }
    }
  }


  // 6.3. Calculate other food stomach content
  for (y=0; y<nyrs; y++){           // Year loop
    for(p=0; p < nspp; p++){                  // Predator species loop
      for (a=0; a < nages(p); a++){           // Predator age loop
        of_stomKir(y , a, p) = Type( 1 );         // Initialize other suitability
        for (i=0; i < nspp; i++){             // Prey species loop
          for (j=0; j < nages(i); j++){       // Prey age loop
            of_stomKir(y , a, p) -= stomKir(p, i, a, j, y);
          }
        }
      }
    }
  }
  of_stomKir /= other_food(0);


  if(msmMode > 0){
    // 6.4. Calculate stomach proportion over biomass; U/ (W * N)
    Type suit_tmp;                            //  Temporary storage variable
    array<Type> suma_suit(nyrs, max_age, nspp); suma_suit.setZero();
    for(y=0; y<nyrs; y++){                    // Year loop
      for(p=0; p < nspp; p++){                // Predator species loop
        for(a=0; a < nages(p); a++){          // Predator age loop
          for (i=0; i < nspp; i++){           // Prey species loop
            for (j=0; j < nages(i); j++){     // Prey age loop
              suit_tmp = stomKir(p, i, a, j, y) / (AvgN(y, j, i));
              if (wt(y, j, i)!= 0){
                stom_div_bio2(p, i, a, j, y) = suit_tmp / wt(y, j, i);
                suma_suit(y, a, p ) += stom_div_bio2(p, i, a, j, y); // Calculate sum of stom_div_bio2 across prey and  prey age for each predator, predator age, and year
              }
            }
          }
        }
      }
    }


    // 6.5. Calculate suitability
    suit_main.setZero();
    for(p=0; p < nspp; p++){                  // Predator species loop
      for (a=0; a < nages(p); a++){           // Predator age loop
        suit_other(p, a) = Type( 1 );         // Initialize other suitability
        for (i=0; i < nspp; i++){             // Prey species loop
          for (j=0; j < nages(i); j++){       // Prey age loop
            for (y=0; y<nyrs; y++){           // Year loop
              suit_main(p, i, a, j) += stom_div_bio2(p, i, a, j, y)/(suma_suit(y, a, p ) + of_stomKir(y , a, p));
            }                                 // End year loop
            suit_main(p, i, a, j) /= nyrs;
            suit_other(p, a) -= suit_main(p, i, a, j); // Subtract observed suitability from entire suitability (i.e. 1)
          }
        }
      }
    }


    // 6.6. Calculate available food
    Type tmp_othersuit = 0  ;
    avail_food.setZero();
    for(p=0; p < nspp; p++){                    // Predator species loop
      for (a=0; a < nages(p); a++){             // Predator age loop
        for(y=0; y < nyrs; y++){                // Year loop
          tmp_othersuit=0.;
          for (i=0; i < nspp; i++){             // Prey species loop
            for (j=0; j < nages(i); j++){       // Prey age loop
              avail_food(y, a, p) += suit_main(p, i, a, j) * AvgN(y, j, i) * wt(y, j, i); // FIXME - include overlap indices: FIXME - mn_wt_stom?
              tmp_othersuit += suit_main(p, i, a, j); // FIXME - include overlap indices
            }
          }
          avail_food(y, a, p) += other_food(p)*(Type(1) - (tmp_othersuit)); // FIXME - double check this is in the right loop
        }
      }
    }


    // 6.7. Calculate predation mortality
    M2.setZero();
    B_eaten.setZero();
    for (i=0; i < nspp; i++){                   // Prey species loop
      for (j=0; j < nages(i); j++){             // Prey age loop
        for(y=0; y < nyrs; y++){                // Year loop
          for (p=0; p <nspp; p++){              // Predator species loop
            for (a=0; a < nages(p); a++){       // Predator age loop
                M2(y, j, i) += (AvgN(y, a, p) * ration2Age(y, a, p) * suit_main(p , i , a, j)) / avail_food(y, a, p); // #FIXME - include indices of overlap
                B_eaten(y, j, i) += AvgN(y, a, p) * ration2Age(y, a, p) * suit_main(p , i , a, j);
            }
          }
        }
      }
    }
  }



  // ------------------------------------------------------------------------- //
  // 7. SURVEY COMPONENTS EQUATIONS                                            //
  // ------------------------------------------------------------------------- //
  // 7.1. Survey selectivity
  srv_sel.setZero();
  avgsel_srv.setZero();
  for (i=0; i<nspp; i++){
    // 7.1.1. Logisitic selectivity
    if(logist_sel_phase(i) > 0){
      for (j=0; j < nages(i); j++)
        srv_sel(i, j) = 1/ (1 + exp( -srv_sel_slp(i) * ((j+1) - srv_sel_inf(i))));
    }

    // 7.1.2. Selectivity fit to age ranges. NOTE: This can likely be improved
    if(logist_sel_phase(i) < 0){
      for(j=0; j < nselages; j++){
        srv_sel(i, j) = srv_sel_coff(i, j);
        avgsel_srv(i) +=  exp(srv_sel_coff(i, j));
      }
      // 7.1.3 Average selectivity up to nselages
      avgsel_srv(i) = log(avgsel_srv(i) / nselages);

      // 7.1.4. Plus group selectivity
      for(j = nselages; j < nages(i); j++){
        srv_sel(i, j) = srv_sel(i, nselages - 1);
      }

      // 7.1.5. Average selectivity across all ages
      Type avgsel_tmp = 0; // Temporary object for average selectivity across all ages
      for(j = 0; j < nages(i); j++){
        avgsel_tmp += exp(srv_sel(i, j));
      }
      avgsel_tmp = log(avgsel_tmp / nages(i));

      // 7.1.6. Standardize selectivity
      for(j = 0; j < nages(i); j++){
        srv_sel(i, j) -=  avgsel_tmp;
        srv_sel(i, j) = exp(srv_sel(i, j));
      }
    }
  }

  // 7.2 EIT Components
  // -- 7.2.1 EIT Survey Biomass
  int eit_yr_ind;
  eit_age_hat.setZero();
  eit_hat.setZero();
  for(j=0; j < nages(0); j++){
    for(y=0; y < n_eit; y++){
      eit_yr_ind = yrs_eit(y) - styr;
      eit_age_hat(y, j) = NByage(eit_yr_ind, j, 0) * eit_sel(eit_yr_ind, j) * eit_q; // Remove the mid-year trawl?
      eit_hat(y) += eit_age_hat(y, j) * wt(eit_yr_ind, j, 0);  //

    }
  }


  // -- 7.2.2 EIT Survey Age Composition
  eit_age_comp_hat.setZero();
  for(j=0; j < nages(0); j++){
    for (y=0; y < n_eit; y++){
      eit_age_comp_hat(y, j) = eit_age_hat(y, j) / eit_age_hat.row(y).sum(); // Divide numbers at age by total numbers for each year
    }
  }


  // 7.3 BT Components
  // -- 7.3.1 BT Survey Biomass
  srv_bio_hat.setZero();
  int srv_yr_ind;
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs_srv_biom(i); y++){
        srv_yr_ind = yrs_srv_biom(i, y) - styr; // Temporary index for years of data
        srv_bio_hat(i, y) += NByage(srv_yr_ind, j, i) * exp(-0.5 * Zed(srv_yr_ind, j, i)) * srv_sel(i, j) * exp(log_srv_q(i)) * wt(srv_yr_ind, j, i);  //
      }
    }
  }

  // -- 7.4.2 BT Survey Age Composition: NOTE: will have to alter if age comp data are not the same length as survey biomass data
  srv_age_hat.setZero();
  srv_hat.setZero();
  vector<Type> srv_age_tmp( max_age ); // Temporary vector of survey-catch-at-age for matrix multiplication
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs_srv_age(i); y++){
        srv_yr_ind = yrs_srv_age(i, y) - styr; // Temporary index for years of data
        srv_age_hat(y, j, i) = NByage(srv_yr_ind, j, i) * srv_sel(i, j) * exp(log_srv_q(i));
        srv_hat(i, y) += srv_age_hat(y, j, i);   // Total numbers
      }
    }
  }

  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_srv_age(i); y++){
      // 7.4.2.1 -- BT Survey catch-at-age
      if(srv_age_type(i)==1){
        for(j=0; j < nages(i); j++){
          srv_age_hat(y, j, i) = srv_age_hat(y, j, i) / srv_hat(i, y);
        }
      }
      // 7.4.2.2 -- Convert from catch-at-age to catch-at-length: NOTE: There has got to be a better way
      if(srv_age_type(i)!=1){
        for(j=0; j < nages(i); j++){
          srv_age_tmp(j) = srv_age_hat(y, j, i);
        }

        matrix<Type> ALK = trim_matrix( matrix_from_array(age_trans_matrix, i), nages(i), srv_age_bins(i) );
        vector<Type> srv_age_tmp_trimmed = trim_vector(srv_age_tmp, nages(i) );
        vector<Type> srv_len_tmp = vec_mat_prod( srv_age_tmp_trimmed , ALK ); // Multiply the ALK for species i against the survey catch-at-age for year y

        for(j=0; j < srv_age_bins(i); j++){
          srv_age_hat(y, j, i) = srv_len_tmp(j) / srv_hat(i, y) ; // * age_trans_matrix.col().col(i)) / srv_hat(i, y); // # NOTE: Double check the matrix algebra here
        }
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 8. FISHERY COMPONENTS EQUATIONS                                           //
  // ------------------------------------------------------------------------- //
  // 8.5. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
  tc_biom_hat.setZero(); // Initialize tc_biom_hat
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs_tc_biom(i); y++){
        fsh_yr_ind = yrs_tc_biom(i, y) - styr; // Temporary index for years of data
        tc_biom_hat(i, y) += F(fsh_yr_ind, j, i)/Zed(fsh_yr_ind, j, i) * (1 - exp(-Zed(fsh_yr_ind, j, i))) * NByage(fsh_yr_ind, j, i) * wt(y, j, i); // 5.5.
      }
    }
  }


  // 8.6. ESTIMATE FISHERY AGE COMPOSITION
  tc_hat.setZero();
  // 8.6.1. Get catch-at-age
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs_tc_biom(i); y++){
        fsh_yr_ind = yrs_tc_biom(i, y) - styr; // Temporary index for years of data
        catch_hat(y, j, i) = F(fsh_yr_ind, j, i)/Zed(fsh_yr_ind, j, i) * (1 - exp(-Zed(fsh_yr_ind, j, i))) * NByage(fsh_yr_ind, j, i); // 5.4.
        tc_hat(i, y) += catch_hat(y, j, i); // Estimate catch in numbers
      }
    }
  }

  // 8.6.2 Convert catch-at-age to age-comp
  vector<Type> fsh_age_tmp( max_age ); // Temporary vector of survey-catch-at-age for matrix multiplication
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs; y++){
      /// 8.7.2.1 -- Estimate age composition of the fishery
      if(fsh_age_type(i)==1){
        for(j=0; j < nages(i); j++){
          fsh_age_hat(y, j, i) = catch_hat(y, j, i) / tc_hat(i, y);
        }
      }

      // 8.7.2.1 -- Convert from catch-at-age to catch-at-length: NOTE: There has got to be a better way
      if(fsh_age_type(i)!=1){
        for(j=0; j < nages(i); j++){
          fsh_age_tmp(j) = catch_hat(y, j, i);
        }

        matrix<Type> ALK = trim_matrix( matrix_from_array(age_trans_matrix, i), nages(i), fsh_age_bins(i) );
        vector<Type> fsh_age_tmp_trimmed = trim_vector(fsh_age_tmp, nages(i) );
        vector<Type> fsh_len_tmp = vec_mat_prod( fsh_age_tmp_trimmed , ALK ); // Multiply the ALK for species i against the survey catch-at-age for year y

        for(j=0; j < fsh_age_bins(i); j++){
          fsh_age_hat(y, j, i) = fsh_len_tmp(j) / tc_hat(i, y) ; // * age_trans_matrix.col().col(i)) / srv_hat(i, y); // # NOTE: Double check the matrix algebra here
        }
      }
    }
  }


  // ------------------------------------------------------------------------- //
  // 9. LIKELIHOOD EQUATIONS                                                   //
  // ------------------------------------------------------------------------- //
  // 9.0. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(13,nspp); // matrix of negative log-likelihood components
  jnll_comp.setZero();
  Type jnll = 0;

  // -- Data components
  // Slot 0 -- BT survey biomass -- NFMS annual BT survey
  // Slot 1 -- BT survey age composition -- NFMS annual BT survey
  // Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
  // Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
  // Slot 4 -- Total catch -- Fishery observer data
  // Slot 5 -- Fishery age composition -- Fishery observer data
  // -- Likelihood penalties
  // Slot 6 -- Fishery selectivity
  // Slot 7 -- Fishery selectivity normalization
  // Slot 8 -- Survey selectivity
  // Slot 9 -- Survey selectivity normalization
  // -- Priors
  // Slot 10 -- Tau -- Annual recruitment deviation
  // Slot 11 -- init_dev -- Initial abundance-at-age
  // Slot 12 -- Epsilon -- Annual fishing mortality deviation


  // 9.1. OFFSETS AND PENALTIES
  // 9.1.1 -- Set up offset objects
  vector<Type> offset_srv(nspp); offset_srv.setZero(); // Offset for multinomial likelihood
  vector<Type> offset_fsh(nspp); offset_srv.setZero(); // Offset for multinomial likelihood
  Type offset_eit = 0; // Offset for multinomial likelihood


  // 9.1.1. -- Fishery age comp offsets
  for(i=0; i < nspp; i++){
    for(y = 0; y < nyrs_fsh_comp(i); y++){
      for(j=0; j < fsh_age_bins(i); j++){
        offset_fsh( i ) -= tau * (fsh_age_obs(y, j, i) + MNConst) * log(fsh_age_obs(y, j, i) + MNConst);
      }
    }
  }

  // 9.1.2. -- Survey age comp offsets
  for(i=0; i < nspp; i++){
    matrix<Type> srv_age_tmp = trim_matrix(matrix_from_array(srv_age_obs, i), nyrs_srv_age(i), srv_age_bins(i));
    for (y = 0; y < nyrs_srv_age(i); y++){
      for(j = 0; j < srv_age_bins(i); j++){
        srv_age_obs(y, j, i) = srv_age_obs(y, j, i) / srv_age_tmp.row(y).sum() ; // Convert numbers-at-age to age comp
        offset_srv(i) -= srv_age_n(i, y)*(srv_age_obs(y, j, i) + MNConst) * log(srv_age_obs(y, j, i) + MNConst ) ;
      }
    }
  }

  // 9.1.3. -- Offsets for acoustic survey
  for(y = 0; y < n_eit; y++){
    for(j = 0; j < nages(0); j++){
      offset_eit -= eit_age_n(y) * (eit_age_comp(y, j) + MNConst) * log(eit_age_comp(y, j) + MNConst );
    }
  }


  // 9.2. FIT OBJECTIVE FUNCTION
  // Slot 0 -- BT survey biomass -- NFMS annual BT survey
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_srv_biom(i); y++){
      // srv_yr_ind = yrs_srv_biom(i, y) - styr;
      jnll_comp(0, i) += pow(log(srv_biom(i, y)) - log(srv_bio_hat(i, y)), 2) / (2 * pow(srv_biom_lse(i, y), 2)); // NOTE: This is not quite the lognormal and biohat will be the median.
    }
  }


  // Slot 1 -- BT survey age composition -- NFMS annual BT survey
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_srv_age(i); y++){
      for(j=0; j < srv_age_bins(i); j++){
        // srv_yr_ind = yrs_srv_age(i, y) - styr;
        jnll_comp(1, i) -= srv_age_n(i, y) * (srv_age_obs(y, j, i) + MNConst) * log(srv_age_hat(y, j, i) + MNConst); // Should srv_age_obs  be in log space?
      }
    }
    jnll_comp(1, i) -= offset_srv(i);
  }


  // Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
  for (y=0; y < n_eit; y++){
    jnll_comp(2, 0) += 12.5 * pow(log(obs_eit(y)) - log(eit_hat(y) + 1.e-04), 2);
  }


  // Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
  for (y=0; y < n_eit; y++){
    for(j=0; j < nages(0); j++){
      // eit_yr_ind = yrs_eit(y) - styr;
      jnll_comp(3, 0) -= eit_age_n(y) *  (eit_age_comp(y, j) + MNConst) * log(eit_age_comp_hat(y, j) + MNConst);
    }
  }
  jnll_comp(3, 0) -= offset_eit;


  // Slot 4 -- Total catch -- Fishery observer data
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_tc_biom(i); y++){
      fsh_yr_ind = yrs_tc_biom(i, y) - styr;
      jnll_comp(4,i) += pow((log(tcb_obs(i, y) + Type(1.e-4)) - log(tc_biom_hat(i, y) + Type(1.e-4))), 2) / (2 * pow(sigma_catch, 2)); // T.4.5
    }
  }


  // Slot 5 -- Fishery age composition -- Fishery observer data
  for(i=0; i < nspp; i++){
    for (y=0; y < nyrs_fsh_comp(i); y++){
      for(j=0; j < fsh_age_bins(i); j++){
        fsh_yr_ind = yrs_fsh_comp(i, y) - styr; // Temporary index for years of data
        jnll_comp(5, i) -= tau * (fsh_age_obs(y, j, i) + MNConst) * log(fsh_age_hat(fsh_yr_ind, j, i) + MNConst );
      }
      jnll_comp(5, i) -= offset_fsh(i);
    }
  }


  // Slot 6 -- Fishery selectivity
  for(i=0; i < nspp; i++){
    for(j=0; j < (nages(i) - 1); j++){
      if( fsh_sel(i, j) > fsh_sel( i, j + 1 )){
        jnll_comp(6, i) += 20 * pow( log( fsh_sel(i, j) / fsh_sel(i, j + 1 ) ), 2);
      }
    }

    // Extract only the selectivities we want
    vector<Type> sel_tmp(nages(i));
    for(j=0; j < nages(i); j++){
      sel_tmp(j) = log(fsh_sel(i, j));
    }

    for(j=0; j < nages(i) - 2; j++){
      sel_tmp(j) = first_difference( first_difference( sel_tmp ) )(j);
      jnll_comp(6, i) += curv_pen_fsh * pow( sel_tmp(j) , 2); // FIX

    }
  }


  // Slot 6 -- Add fishery selectivity normalization
  for(i=0; i < nspp; i++){
    jnll_comp(7, i) += 50 * pow(avgsel_fsh(i), 2);
  }



  // Slot 7 -- Survey selectivity
  for(i=0; i < nspp; i++){
    if (logist_sel_phase(i) < 0){
      // Extract only the selectivities we want
      vector<Type> sel_tmp(nages(i));
      for(j=0; j < nages(i); j++){
        sel_tmp(j) = log(srv_sel(i, j));
      }
      for(j=0; j < nages(i) - 2; j++){
        sel_tmp(j) = first_difference( first_difference( sel_tmp ) )(j);
        jnll_comp(8, i) += curv_pen_fsh * pow( sel_tmp(j) , 2); // FIX
      }
    }
  }


  // Slot 7 -- Add survey selectivity normalization
  for(i=0; i < nspp; i++){
    jnll_comp(9, i) += 50 * pow(avgsel_srv(i), 2);
  }


  // Slots 10-12 -- PRIORS: PUT RANDOM EFFECTS SWITCH HERE
  for(i=0; i < nspp; i++){
    // Slot 10 -- init_dev -- Initial abundance-at-age
    for(j=1; j < nages(i); j++){
      jnll_comp(10, i) += pow( init_dev(i,j-1), 2);
    }

    // Slot 11 -- Tau -- Annual recruitment deviation
    // Slot 12 -- Epsilon -- Annual fishing mortality deviation
    for (y=0; y < nyrs; y++){
      jnll_comp(11, i) += pow( rec_dev(i,y), 2);     // Recruitment deviation using penalized likelihood.
      jnll_comp(12, i) += pow( F_dev(i,y), 2);       // Fishing mortality deviation using penalized likelihood.
    }
  }

  // ------------------------------------------------------------------------- //
  // 10. SIMULATION SECTION                                                    //
  // ------------------------------------------------------------------------- //


  // ------------------------------------------------------------------------- //
  // 11. REPORT SECTION                                                        //
  // ------------------------------------------------------------------------- //

  // 11.1. Population components
  REPORT( Zed );
  REPORT( NByage );
  REPORT( AvgN );
  REPORT( S );
  REPORT( biomassByage );
  REPORT( biomassSSBByage );
  REPORT( biomass );
  REPORT( biomassSSB );
  REPORT( pmature );
  REPORT( R );
  REPORT( M1 );


  // -- 11.2. Survey components
  REPORT( srv_age_obs );
  REPORT( srv_bio_hat );
  REPORT( srv_hat );
  REPORT( srv_age_hat );
  REPORT( eit_hat );
  REPORT( eit_age_hat );
  REPORT( eit_age_comp_hat )
    REPORT( obs_eit_age );
  REPORT( eit_age_comp );
  REPORT( avgsel_srv );
  REPORT( srv_sel );


  // -- 11.3. Fishery components
  REPORT( F );
  REPORT( F_dev );
  REPORT( fsh_sel );
  REPORT( avgsel_fsh );
  REPORT( tc_biom_hat );
  REPORT( catch_hat );
  REPORT( tc_hat );
  REPORT( fsh_age_hat );
  REPORT( fsh_age_obs );


  // -- 11.4. Likelihood components
  REPORT( jnll_comp );
  REPORT( offset_srv );
  REPORT( offset_fsh );
  REPORT( offset_eit );


  // -- 11.5. Ration components
  REPORT( ConsumAge );
  REPORT( Consum_livingAge );
  REPORT( S2Age );
  REPORT( LbyAge );
  REPORT( mnWt_obs );
  REPORT( fT );
  REPORT( TempC );
  REPORT( ration2Age );


  // 11.6. Suitability components
  REPORT( suit_main );
  REPORT( suit_other );
  REPORT( stom_div_bio2 );
  REPORT( stomKir );
  REPORT( avail_food );
  REPORT( of_stomKir );
  REPORT( M2 );
  REPORT( B_eaten );


  // ------------------------------------------------------------------------- //
  // END MODEL                                                                 //
  // ------------------------------------------------------------------------- //
  if(debug == 0){
    jnll = jnll_comp.sum();
  }
  if(debug == 1){
    jnll = dummy * dummy;
  }

  REPORT( jnll );


  return jnll;
}
