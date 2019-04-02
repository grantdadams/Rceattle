#include <TMB.hpp>
// ------------------------------------------------------------------------- //
//                 CEATTLE version 3.1.2                                     //
//                  Template Model Builder                                   //
//               Multispecies Statistical Model                              //
//          Bioenergetic-based Assessment for Understanding                  //
//              Biomass Linkages To The Environment                          //
//                  in the Bering Sea                                        //
//                                                                           //
// AUTHORS:   Kirstin Holsman, Jim Ianelli                                   //
//            Modified by Grant Adams                                        //
// CITATIONS:                                                                //
// 1. Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2015). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep-Sea Research Part II: Topical Studies in Oceanography, 134, 360–378. https://doi.org/10.1016/j.dsr2.2015.08.001
// ------------------------------------------------------------------------- //
//
//  INDEX:
//  0. Load dependencies
//  1. Model configuration
//  2. Model inputs
//  3. Model parameters
//  4. Derived quantities
//  5. Initial calculations
//  6. Population dynamics
//  7. Ration equations
//  8. Predation mortality equations
//  -- 8.1. Holsman et al Predation mortality
//  -- 9.2. Kinzey and Punt predation mortality
//  9. Survey components
//  10. Fishery components
//  11. Likelihood components
//  12. Report section
//  13. End
//
//
// ------------------------------------------------------------------------- //
// 0. LOAD DEPENDENCIES                                                      //
// ------------------------------------------------------------------------- //
#include "functions.hpp"


template<class Type>
Type objective_function<Type>::operator() () {
  // ------------------------------------------------------------------------- //
  // 1. MODEL CONFIGURATION                                                    //
  // ------------------------------------------------------------------------- //
  // 1.1. CONFIGURE MODEL (this section sets up the switches)
  DATA_INTEGER(debug);            // Logical to debug or not
  DATA_INTEGER(msmMode);
  //    0 = run in single species mode
  //    1 = run in MSM mode  Holsman et al (2015) MSVPA based
  // DATA_INTEGER(est_diet);         // Include diet data in the likelihood
  DATA_INTEGER(suitMode);         // Estimate suitability
  DATA_INTEGER(avgnMode);         // N used for predation function
  //    0 = AvgN
  //    1 = N*exp(-Z / 2))
  //    2 = N

  DATA_INTEGER( random_rec );     // Logical of whether to treate recruitment deviations as random effects
  DATA_INTEGER( niter );          // Number of loops for MSM mode

  DATA_IVECTOR( srv_sel_type ); // Selectivity type for BT survey
  //    0 = fit to data
  //    1 = logistic

  // 1.2. Temporal dimensions
  DATA_INTEGER( nyrs );           // Number of estimation years
  DATA_INTEGER( styr );           // Start year
  // DATA_INTEGER( proj_yr );         // Year to project with no fishing

  // 1.3. Number of species
  DATA_INTEGER( nspp );           // Number of species (prey)


  // 1.4. MODEL OBJECTS
  // 1.4.1. LOOPING INDICES -- k = observation, sp = species, age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length
  int sp, age, ln, ksp, k_age, k_ln, yr, rsp, r_age, r_ln; // k
  int fsh_yr, srv_yr, eit_yr;         // Year indices
  int ncnt;    // Pointers
  if (msmMode == 0) { niter = 1; }    // Number of iterations for SS mode

  // ------------------------------------------------------------------------- //
  // 2. MODEL INPUTS                                                           //
  // ------------------------------------------------------------------------- //

  // 2.1. FIXED VALUES
  int tau = 200;                    // Fishery age composition sample size
  Type MNConst = 0.001;             // Constant additive for logistic functions
  Type curv_pen_fsh = 12.5;         // Fishery selectivity penalty
  Type sigma_catch = 0.05;          // SD of catch
  Type sd_ration = 0.05;            // SD of ration likelihood


  // 2.2. DIMENSIONS
  // int nyrs_proj = proj_yr - styr + 1;        // End year

  // -- 2.2.2. Species attributes
  DATA_IVECTOR( nages );          // Number of species (prey) ages; n = [1, nspp]
  DATA_INTEGER(stom_tau);         // Stomach sample size
  DATA_IVECTOR( nselages );       // Number of ages to estimate selectivity; n = [1, nspp]
  int max_age = imax(nages);      // Integer of maximum nages to make the arrays; n = [1]

  // 2.3. DATA INPUTS (i.e. assign data to objects)
  // -- 2.3.1 Fishery Components
  DATA_IVECTOR( nyrs_tc_biom );   // Number of years with total observed catch; n = [nspp]
  DATA_IMATRIX( yrs_tc_biom );    // Years with total observed catch; n = [nspp, nyrs_tc_biom]
  DATA_MATRIX( tcb_obs );         // Observed total yield (kg); n = [nspp, nyrs_tc_biom]

  DATA_IVECTOR( nyrs_fsh_comp );  // Number of years in the fishery sp_age composition data; n = [nspp]
  DATA_IMATRIX( yrs_fsh_comp );   // Years for the fishery sp_age composition data; n = [nspp, nyrs_fsh_comp]
  DATA_IVECTOR( fsh_age_type );   // Composition type for fishery catch (1 = age based, 2 = length based); n = [nspp]
  DATA_IVECTOR( fsh_age_bins );   // Number of bins for fishery age/length composition data; n = [nspp]
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
  // DATA_MATRIX( srv_age_sizes );   // Observed size composition: n = [nspp, srv_age_bins]; NOTUSED
  DATA_ARRAY( age_trans_matrix);  // observed sp_age/size compositions; n = [nspp, nages, srv_age_bins]

  // -- 2.3.4 EIT Survey components
  DATA_INTEGER( n_eit);           // Number of years with EIT data; n = [1]
  DATA_IVECTOR( yrs_eit);         // Years for available EIT data; n = [n_eit]
  DATA_VECTOR( eit_age_n);        // Number  of  EIT Hauls for multinomial; n = [yrs_eit]
  DATA_MATRIX( obs_eit_age);      // Observed EIT catch-at-age; n = [n_eit, nyrs]
  DATA_VECTOR( obs_eit);          // Observed EIT survey biomass (kg); n = [1, nyrs]
  DATA_MATRIX( eit_sel);          // Observed EIT survey selectivity; n = [eit_age, nyrs_eit_sel]

  // -- 2.3.5. Weight-at-age
  // DATA_IVECTOR( nyrs_wt_at_age ); // Number of years of weight at age data; n = [nspp]: NOTUSED
  // DATA_IMATRIX( yrs_wt_at_age );  // Years of weight-at-age; data n = [nspp, nyrs_wt_at_age]: NOTUSED
  DATA_ARRAY( wt );               // Weight-at-age by year; n = [nyrs, nages, nspp]: FIXME: Change nyrs to nyrs_wt_at_age if data don't span entire bit

  // 2.3.6. Diet data
  // DATA_MATRIX( maxK );            // Matrix of maximum diet proportions for each predd,preyy combo (used for broken stick expon); n = [nspp, nspp]: NOTUSED
  // DATA_VECTOR( klim );            // Max  pred  length  at  which K should  be  calculated  (beyond this  is  it  just  equal to  K); n = [1, nspp]: NOTUSED
  // DATA_MATRIX( prefa );           // prefa of the functional response prefa(pred,prey); n = [nspp, nspp]: NOTUSED
  // DATA_MATRIX( prefb );           // prefb of the functional response prefa(pred,prey); n = [nspp, nspp]: NOTUSED
  // DATA_MATRIX( l1 );              // l1 of the switch function l1(pred,prey); n = [nspp, nspp]: NOTUSED
  // DATA_MATRIX( l2 );              // l2 of the switch function l2(pred,prey); n = [nspp, nspp]: NOTUSED
  // DATA_MATRIX( LA_Lengths );      // length to sp_age coversion matrix (rounded to nearest integer); n = [nspp, nages]: NOTUSED
  DATA_VECTOR( fday );            // number of foraging days for each predator; n = [1, nspp] #FIXME - assuming this is the same as fdays
  // DATA_IVECTOR( nlengths );       // number of Lengths for the matrix to convert ages to lengths; n = [1, nspp]: NOTUSED
  // int maxL = imax( nlengths );    // Maximum number of lengths for the matrix to convert between ages and lengths; n = [1]: NOTUSED
  // DATA_MATRIX( lengths );         // Lengths for the matrix to convert ages to lenghths; n = [nspp, nlengths]: NOTUSED
  // DATA_ARRAY( A2L_matrix );       // Age to length matrix : Matrix to convert ages to lenghths; n = [nspp, nages, nlengths]: NOTUSED
  // DATA_ARRAY( K );                // K(1) is nprey by maxL matrix of stomach proportions of predator 1; n = [nspp, nspp, nlengths]: NOTUSED
  // DATA_ARRAY( KAge );             // K(1) is nprey by maxL matrix of stomach proportions of predator 1; n = [nspp, nspp, nages]: NOTUSED
  // DATA_MATRIX( PAge );            // n = [nspp, nages]: NOTUSED
  DATA_ARRAY( Pyrs );             // n = [nspp, nyrs+1, nages]: #FIXME - Assuming this is the same as Pby_yr?
  DATA_ARRAY( Uobs );             // pred, prey, predL, preyL U matrix (mean number of prey in each pred); n = [nspp, nspp, maxL, maxL]
  DATA_ARRAY( UobsWt );           // pred, prey, predL, preyL U matrix (mean wt_hat of prey in each pred); n = [nspp, nspp, maxL, maxL] #FIXME - Changed name in stomach2017.dat
  DATA_ARRAY( UobsAge );          // pred, prey, predA, preyA U matrix (mean number of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_ARRAY( UobsWtAge );        // pred, prey, predA, preyA U matrix (mean wt_hat of prey in each pred age); n = [nspp, nspp, max_age, max_age]
  DATA_MATRIX( Mn_LatAge );       // Mean length-at-age; n = [nspp, nages], ALSO: mean_laa in Kinzey

  // -- 2.3.7. Q Diet data
  // DATA_INTEGER( npred3 );      // Number of predators; n = [1]: NOTUSED
  // DATA_IVECTOR( nprey3 );      // Number of prey species for each predator; n = [1, nspp]: NOTUSED
  // DATA_IVECTOR( n_stomyrs );      // number years of stomach data; n = [1, nspp]: NOTUSED
  // DATA_IMATRIX( stomyrs );        // years of stomach data; n = [nspp, n_stomyrs]: NOTUSED
  // DATA_ARRAY( mnQ );              // meanwt of each prey spp in the stomach of each predator of sp_age a; n = [npred3,n_stomyrs,max_age,nspp+1]: NOTUSED
  // DATA_ARRAY( Qse );              // SE wt of each prey spp in the stomach of each predator of sp_age a; n = [npred3,n_stomyrs,max_age,nspp+1]: NOTUSED

  // 2.3.8. Environmental data
  DATA_INTEGER( nTyrs );          // Number of temperature years; n = [1] #FIXME - changed the name of this in retro_data2017_asssmnt.dat
  DATA_IVECTOR( Tyrs );           // Years of hindcast data; n = [1, nTyrs] #FIXME - changed the name of this in retro_data2017_asssmnt.dat
  DATA_VECTOR( BTempC_retro );    // Vector of bottom temperature; n = [1,  nTyrs ]
  // DATA_INTEGER( ncov );        // Number of environmental covariates

  // 2.4. INPUT PARAMETERS
  // -- 2.4.1. Bioenergetics parameters (BP)
  DATA_VECTOR( other_food );      // Biomass of other prey (kg); n = [1, nspp]
  // DATA_IVECTOR( useWt );          // Assign relative proportion of prey in the diet according to relative biomass in the system.,otherwise the model with use relative proportion by number; n = [1, nspp]: NOTUSED
  DATA_IVECTOR( C_model );        // f == 1, the use Cmax*fT*P; n = [1, nspp]
  DATA_VECTOR( Pvalue );          // This scales the pvalue used if C_model ==1 , proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge; n = [1, nspp]
  DATA_IVECTOR( Ceq );            // Ceq: which Comsumption equation to use; n = [1, nspp]; Currently all sp = 1
  DATA_VECTOR( CA );              // Wt specific intercept of Cmax=CA*W^CB; n = [1, nspp]
  DATA_VECTOR( CB );              // Wt specific slope of Cmax=CA*W^CB; n = [1, nspp]
  DATA_VECTOR( Qc );              // used in fT, QC value; n = [1, nspp]
  DATA_VECTOR( Tco );             // used in fT, thermal optimum; n = [1, nspp]
  DATA_VECTOR( Tcm );             // used in fT, thermal max; n = [1, nspp]
  DATA_VECTOR( Tcl );             // used in fT eq 3, limit; n = [1, nspp]
  DATA_VECTOR( CK1 );             // used in fT eq 3, limit where C is .98 max (ascending); n = [1, nspp]
  DATA_VECTOR( CK4 );             // used in fT eq 3, temp where C is .98 max (descending); n = [1, nspp]
  DATA_MATRIX( S_a );             // S_a, S_b, S_b2, S_b3, S_b4, S_b5: a,L,L^2,L^3,L^4,L^5 (rows)coef for mean S=a+b*L+b2*L*L, whith a cap at 80cm for each pred spp(cols); n = [6, nspp]
  // DATA_MATRIX( aLim );            // #aLim and bLim : upper limit of prey size; n = [2, nspp]: NOTUSED

  // -- 2.4.2. von Bertalannfy growth function (VBGF): This is used to calculate future weight-at-age: NOT YET IMPLEMENTED
  // DATA_VECTOR( t0 );              // t0 parameter of the temp specific VonB for wt; n = [nspp, mf_type]
  // DATA_VECTOR( log_mean_d );      // log mean d parameter of the temp specific von B for wt; n = [nspp, mf_type]
  // DATA_VECTOR( logK );            // log k parameter of the temp specific VonB for wt; n = [nspp, mf_type]
  // DATA_VECTOR( logH );            // log H parameter of the temp specific VonB for wt; n = [nspp, mf_type]
  // DATA_VECTOR( Tcoef );           // T coefficent of the linear d equations of the temp specific VonB for wt; n = [nspp, mf_type]
  // DATA_VECTOR( Pcoef );           // P-value coefficent of the linear d equations of the temp specific VonB for wt; n = [nspp, mf_type]

  // -- 2.4.3. Weight-at-length parameters
  DATA_MATRIX( aLW );             // LW a&b regression coefs for W=a*L^b; n = [2, nspp]

  // -- 2.4.4. Others
  DATA_MATRIX( M1_base );         // Residual natural mortality; n = [nspp, nages]
  DATA_IVECTOR( mf_type );        // Sex specific mort and weight at age? : 1 = same for both, 2 = seperate wt at sp_age for each sex
  DATA_MATRIX( propMorF );        // Proportion-at-age of females of population; n = [nspp, nages]
  DATA_MATRIX( pmature );         // Proportion of mature females at age; [nspp, nages]

  // -- 2.4.5. F Profile data: NOTUSED
  /*
   DATA_INTEGER( n_f );            // Number of F vectors
   DATA_VECTOR( Frates );          // Fishing mortality vector; n = [1, n_f]
   DATA_INTEGER( np );             // Number of
   DATA_MATRIX( Fprofiles );       // F profiles; n = [np, nspp]
   */



  // ------------------------------------------------------------------------- //
  // 2.7. Debugging with data inputs                                           //
  // ------------------------------------------------------------------------- //
  if (debug == 1) {
    // -- 2.8.2.1 Check to make sure the first year of survey data are not before start year
    for (sp = 0; sp < nspp; sp++) {
      if (yrs_tc_biom(sp, 0) < styr) {
        std::cerr << "First year of total catch biomass of species " << sp + 1 << " is before specified start year" << std::endl;
        return (0);
      }
    }

    // -- 2.8.2.2 Check to make sure the first year of survey data are not before start year
    for (sp = 0; sp < nspp; sp++) {
      if (yrs_srv_biom(sp, 0) < styr) {
        std::cerr << "First year of survey biomass of species " << sp + 1 << " is before specified start year" << std::endl;
        return (0);
      }
    }

    // -- 2.8.2.3. Check to make sure the years of survey data biomass and age are the same
    for (sp = 0; sp < nspp; sp++) {
      if (nyrs_srv_biom(sp) != nyrs_srv_age(sp)) {
        std::cerr << "Nyrs of survey biomass and age-comp of species " << sp + 1 << " do not match" << std::endl;
        return (0);
      }
    }

    if (yrs_eit(0) < styr) {
      std::cerr << "First year of EIT survey biomass is before specified start year" << std::endl;
      return (0);
    }
  }


  // ------------------------------------------------------------------------- //
  // 3. PARAMETER SECTION                                                      //
  // ------------------------------------------------------------------------- //
  PARAMETER(dummy);                    // Variable to test derived quantities given input parameters; n = [1]

  // -- 3.1. Recruitment parameters
  PARAMETER_VECTOR( ln_mn_rec );       // Mean recruitment; n = [1, nspp]
  PARAMETER_VECTOR( ln_rec_sigma );    // Standard deviation of recruitment deviations; n = [1, nspp]
  PARAMETER_MATRIX( rec_dev );         // Annual recruitment deviation; n = [nspp, nyrs]
  // PARAMETER(sigma_rec);             // Standard deviation of recruitment variation # NOTE: Have this estimated if using random effects.

  // -- 3.2. Abundance parameters
  PARAMETER_MATRIX( init_dev );             // Initial abundance-at-age; n = [nspp, nages] # NOTE: Need to figure out how to best vectorize this

  // -- 3.3. fishing mortality parameters
  PARAMETER_VECTOR( ln_mean_F );      // Log mean fishing mortality; n = [1, nspp]
  PARAMETER_MATRIX( F_dev );          // Annual fishing mortality deviations; n = [nspp, nyrs] # NOTE: The size of this will likely change
  // PARAMETER_VECTOR(sigma_F); // SD of fishing mortality deviations; n = [1, nspp] # NOTE: Have this estimated if using random effects.

  // -- 3.4. Selectivity parameters
  PARAMETER_MATRIX( srv_sel_coff );    // Survey selectivity parameters; n = [nspp, nselages]
  PARAMETER_MATRIX( fsh_sel_coff );    // Fishery age selectivity coef; n = [nspp, nselages]
  PARAMETER_MATRIX( srv_sel_slp );     // Survey selectivity paramaters for logistic; n = [2, nspp]
  PARAMETER_MATRIX( srv_sel_inf );     // Survey selectivity paramaters for logistic; n = [2, nspp]
  PARAMETER( log_eit_q );              // EIT Catchability; n = [1]
  PARAMETER_VECTOR( log_srv_q );       // BT Survey catchability; n = [1, nspp]

  // FIXME: Create maps
  // -- 3.5. Kinzery predation function parameters
  PARAMETER_MATRIX(logH_1);            // Predation functional form; n = [nspp, nspp2]; // FIXME: make matrix; nspp2 = nspp + 1
  PARAMETER_VECTOR(logH_1a);           // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix
  PARAMETER_VECTOR(logH_1b);           // Age adjustment to H_1; n = [1, nspp]; // FIXME: make matrix

  PARAMETER_MATRIX(logH_2);            // Predation functional form; n = [nspp, nspp]
  PARAMETER_MATRIX(logH_3);            // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
  PARAMETER_MATRIX(H_4);               // Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;

  // 3.6 Gamma selectivity parameters
  PARAMETER_VECTOR( log_gam_a );       // Log predator selectivity; n = [1,nspp]; FIXME: bounds = 1.0e-10 and 19.9
  PARAMETER_VECTOR( log_gam_b );       // Log predator selectivity; n = [1,nspp]; FIXME: bounds = -5.2 and 10

  // 3.7. Preference
  PARAMETER_MATRIX( phi );             // Species preference coefficient; n = [nspp, nspp]


  // ------------------------------------------------------------------------- //
  // 4. DERIVED QUANTITIES SECTION                                             //
  // ------------------------------------------------------------------------- //

  // 4.1. Derived indices
  int max_bin = imax(srv_age_bins);                                              // Integer of maximum number of length/age bins.

  // -- 4.2. Estimated population parameters
  vector<Type>  mn_rec = exp(ln_mn_rec);                                        // Mean recruitment; n = [1, nspp]
  array<Type>   AvgN(nspp, max_age, nyrs); AvgN.setZero();                      // Average numbers-at-age; n = [nspp, nages, nyrs]
  array<Type>   biomassByage(nspp, max_age, nyrs); biomassByage.setZero();      // Estimated biomass-at-age (kg); n = [nspp, nages, nyrs]
  matrix<Type>  biomass(nspp, nyrs); biomass.setZero();                         // Estimated biomass (kg); n = [nspp, nyrs]
  matrix<Type>  biomassSSB(nspp, nyrs); biomassSSB.setZero();                   // Estimated spawning stock biomass (kg); n = [nspp, nyrs]
  array<Type>   biomassSSBByage(nspp, max_age, nyrs); biomassSSBByage.setZero();// Spawning biomass at age (kg); n = [nspp, nages, nyrs]
  array<Type>   M(nspp, max_age, nyrs); M.setZero();                            // Total natural mortality at age; n = [nyrs, nages, nspp]
  matrix<Type>  M1(nspp, max_age); M1.setZero();                                // Base natural mortality; n = [nspp, nages]
  array<Type>   M2(nspp, max_age, nyrs); M2.setZero();                          // Predation mortality at age; n = [nyrs, nages, nspp]
  array<Type>   NByage(nspp, max_age, nyrs); NByage.setZero();                  // Numbers at age; n = [nspp, nages, nyrs]
  matrix<Type>  R(nspp, nyrs); R.setZero();                                     // Estimated recruitment (n); n = [nspp, nyrs]
  array<Type>   S(nspp, max_age, nyrs); S.setZero();                            // Survival at age; n = [nspp, nages, nyrs]
  array<Type>   Zed(nspp, max_age, nyrs); Zed.setZero();                        // Total mortality at age; n = [nspp, nages, nyrs]
  vector<Type>  r_sigma(nspp); r_sigma.setZero();                               // Standard deviation of recruitment variation

  // -- 4.3. Fishery observations
  vector<Type>  fsh_age_tmp( max_age );                                         // Temporary vector of survey-catch-at-age for matrix multiplication
  vector<Type>  avgsel_fsh(nspp); avgsel_fsh.setZero();                         // Average fishery selectivity
  array<Type>   catch_hat(nspp, max_age, nyrs); catch_hat.setZero();            // Estimated catch-at-age (n); n = [nspp, nages, nyrs]
  array<Type>   F(nspp, max_age, nyrs); F.setZero();                            // Estimated fishing mortality; n = [nspp, nages, nyrs]
  array<Type>   fsh_age_obs(nspp, max_bin, nyrs); fsh_age_obs.setZero();        // Observed fishery age comp; n = [nyrs_fsh_comp, fsh_age_bins, nspp]
  array<Type>   fsh_age_hat(nspp, max_bin, nyrs); fsh_age_hat.setZero();        // Estimated fishery age comp; n = [nspp, nages, nyrs]
  matrix<Type>  fsh_sel(nspp, max_age); fsh_sel.setZero();                      // Log estimated fishing selectivity
  matrix<Type>  tc_biom_hat(nspp, nyrs); tc_biom_hat.setZero();                 // Estimated total yield (kg); n = [nspp, nyrs]
  matrix<Type>  tc_hat(nspp, nyrs); tc_hat.setZero();                           // Estimated total catch (n); n = [nspp, nyrs]
  matrix<Type>  tc_obs(nspp, nyrs); tc_obs.setZero();                           // Set total catch to 0 to initialize // Observed total catch (n); n = [nspp, nyrs] NOTE: This may not be necessary if loading data from tmp

  // -- 4.4. BT Survey components
  Type avgsel_tmp = 0;                                                          // Temporary object for average selectivity across all ages
  vector<Type>  srv_age_tmp( max_age ); srv_age_tmp.setZero();                  // Temporary vector of survey-catch-at-age for matrix multiplication
  vector<Type>  avgsel_srv(nspp); avgsel_srv.setZero();                         // Average survey selectivity; n = [1, nspp]
  array<Type>   srv_age_hat(nspp, max_bin, nyrs); srv_age_hat.setZero();        // Estimated BT age comp; n = [nspp, nages, nyrs]
  matrix<Type>  srv_bio_hat(nspp, nyrs); srv_bio_hat.setZero();                 // Estimated BT survey biomass (kg); n = [nspp, nyrs]
  matrix<Type>  srv_hat(nspp, nyrs); srv_hat.setZero();                         // Estimated BT survey total abundance (n); n = [nspp, nyrs]
  matrix<Type>  srv_sel(nspp, max_age); srv_sel.setZero();                      // Estimated survey selectivity at age; n = [nspp, nyrs]
  matrix<Type>  srv_biom_lse(nspp, nyrs); srv_biom_lse.setZero();               // Observed annual biomass CV; n = [nspp, nyrs_srv_biom]

  // -- 4.5. EIT Survey Components
  matrix<Type>  eit_age_comp_hat(srv_age_bins(0), nyrs); eit_age_comp_hat.setZero(); // Estimated EIT age comp ; n = [nyrs, 12 ages]
  matrix<Type>  eit_age_comp = obs_eit_age; eit_age_comp.setZero();              // Eit age comp; n = [n_eit, srv_age_bins(0)]
  matrix<Type>  eit_age_hat(srv_age_bins(0), nyrs); eit_age_hat.setZero();       // Estimated EIT catch-at-age ; n = [nyrs, 12 ages]
  vector<Type>  eit_hat(nyrs); eit_hat.setZero();                                // Estimated EIT survey biomass (kg); n = [nyrs]
  Type eit_q = exp(log_eit_q);                                                   // EIT Catchability

  // -- 4.6. Ration components
  array<Type>   ConsumAge( nspp, max_age, nyrs ); ConsumAge.setZero();            // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
  array<Type>   Consum_livingAge( nspp, max_age, nyrs ); Consum_livingAge.setZero(); // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
  matrix<Type>  fT( nspp, nTyrs ); fT.setZero();                                  // Pre-allocation of temperature function of consumption; n = [nspp, nTyrs]
  array<Type>   LbyAge( nspp, max_age, nyrs ); LbyAge.setZero();                  // Length by age from LW regression
  matrix<Type>  mnWt_obs( nspp, max_age ); mnWt_obs.setZero();                    // Mean observed weight at age (across years); n = [nspp, nages]
  array<Type>   ration2Age( nspp, max_age, nyrs ); ration2Age.setZero();          // Annual ration at age (kg/yr); n = [nyrs, nages, nspp]
  array<Type>   S2Age( nspp, max_age, nyrs ); S2Age.setZero();                    // pre-allocate mean stomach weight as a function of sp_age
  vector<Type>  TempC( nTyrs ); TempC.setZero();                                  // Bottom temperature; n = [1, nTyrs]

  // -- 4.7. Suitability components
  Type tmp_othersuit = 0  ;
  Type suit_tmp = 0;                                                              //  Temporary storage variable
  array<Type>   avail_food(nspp, max_age, nyrs); avail_food.setZero();            // Available food to predator; n = [nyrs, nages, nspp]
  array<Type>   B_eaten(nspp, max_age, nyrs); B_eaten.setZero();                  // Biomass of prey eaten via predation; n = [nyrs, nages, nspp]
  array<Type>   of_stomKir(nspp, max_age, nyrs); of_stomKir.setZero();            // Other food stomach content; n = [nyrs, nages, nspp]
  array<Type>   stom_div_bio2(nspp, nspp, max_age, max_age, nyrs); stom_div_bio2.setZero();// Stomach proportion over biomass; U/ (W * N) ; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   stomKir(nspp, nspp, max_age, max_age, nyrs); stomKir.setZero();       // Stomach proportion by numbers U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   stomKirWt(nspp, nspp, max_age, max_age, nyrs); stomKirWt.setZero();   // Stomach proportion by weight U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   diet_w_dat(nspp, nspp, max_bin, max_bin, nyrs); diet_w_dat.setZero(); // Observed stomach contents by weight of prey length j in predator length l
  array<Type>   diet_w_sum(nspp, nspp, max_bin, nyrs); diet_w_sum.setZero();          // Observed stomach contentes by weight of prey in predator length j
  array<Type>   suit_main(nspp, nspp, max_age, max_age, nyrs); suit_main.setZero();         // Suitability/gamma selectivity of predator age u on prey age a; n = [nspp, nspp, nages, nages]
  matrix<Type>  suit_other(nspp, max_age); suit_other.setZero();                      // Suitability not accounted for by the included prey; n = [nspp, nages]
  array<Type>   suma_suit(nspp, max_age, nyrs); suma_suit.setZero();                  // Sum of suitabilities; n = [nyrs, nages, nspp]
  array<Type>   UobsWtAge_hat(nspp, nspp, max_age, max_age, nyrs); UobsWtAge_hat.setZero();   // Estimated stomach proportion by weight U; n = [nspp, nspp, nages, nages, nyrs]
  array<Type>   mn_UobsWtAge_hat(nspp, nspp, max_age, max_age); mn_UobsWtAge_hat.setZero();   // Average estimated stomach proportion by weight U; n = [nspp, nspp, nages, nages]

  // -- 4.8. Kinzey selectivity
  vector<Type> gam_a = exp(log_gam_a); // Predator selectivity
  vector<Type> gam_b = exp(log_gam_b); // Predator selectivity

  // -- 4.9. Kinzey Functional response bits
  matrix<Type> H_1(nspp, nspp + 1); H_1 = exp(logH_1.array());                                     // FIXME: make matrix
  vector<Type> H_1a(nspp); H_1a = exp(logH_1a);                                      // FIXME: make matrix
  vector<Type> H_1b(nspp); H_1b = exp(logH_1b);                                      // FIXME: make matrix
  matrix<Type> H_2(nspp, nspp); H_2 = exp(logH_2.array());
  matrix<Type> H_3(nspp, nspp); H_3 = exp(logH_3.array());

  array<Type>  N_pred_yrs(nspp, max_age, nyrs); N_pred_yrs.setZero();                // Effective numbers of predators for each age of prey
  array<Type>  N_prey_yrs(nspp, max_age, nyrs); N_prey_yrs.setZero();                // Effective numbers of prey for each age of predator
  matrix<Type> N_pred_eq(nspp, max_age); N_pred_eq.setZero();                        // Effective numbers of predators for each age of prey (styr_pred)
  matrix<Type> N_prey_eq(nspp, max_age); N_prey_eq.setZero();                        // Effective numbers of prey for each age of predator

  array<Type>  pred_resp(nspp, nspp + 1, max_age, max_age, nyrs); pred_resp.setZero();     // Predator functional response
  array<Type>  Pred_r(nspp, max_age, nyrs); Pred_r.setZero();                        // save Pred_ratio values
  array<Type>  Prey_r(nspp, max_age, nyrs); Prey_r.setZero();                        // save Prey_ratio values

  array<Type>  Vmort_ua(nspp, nspp, max_age, max_age, nyrs); Vmort_ua.setZero();     // Predation mortality on prey age a by single predator age u
  array<Type>  Pmort_ua(nspp, nspp, max_age, nyrs); Pmort_ua.setZero();              // Predation mortality on prey age a by all predators age u

  array<Type>  eaten_la(nspp, nspp, max_bin, max_age, nyrs); eaten_la.setZero();     // Number of prey of age a eaten by predator length l
  array<Type>  eaten_ua(nspp, nspp, max_age, max_age, nyrs); eaten_ua.setZero();     // Number of prey of age a eaten by predator age u

  array<Type>  Q_mass_l(nspp, nspp + 1, max_bin, nyrs); Q_mass_l.setZero();                // Mass of each prey sp consumed by predator at length // FIXME: make into 4D array
  array<Type>  Q_mass_u(nspp, nspp + 1, max_bin, nyrs); Q_mass_u.setZero();                // Mass of each prey sp consumed by predator at age // FIXME: make into 4D array
  matrix<Type> Q_other_u(nspp, max_age); Q_other_u.setZero();                        // Mass of other prey consumed by predator at age
  array<Type>  Q_hat(nspp, nspp + 1, max_bin, nyrs); Q_hat.setZero();                      // Fraction for each prey type of total mass eaten by predator length
  array<Type>  T_hat(nspp, nspp, max_bin, max_bin); T_hat.setZero();                 // Fraction of prey of length m in predator of length l

  array<Type> omega_hat(nspp, max_age, nyrs); omega_hat.setZero();                   // Daily ration by predator age each year
  matrix<Type> omega_hat_ave(nspp, max_age); omega_hat_ave.setZero();                // Daily ration by predator age averaged over years


  // ------------------------------------------------------------------------- //
  // 5. INITIAL CALCULATIONS                                                   //
  // ------------------------------------------------------------------------- //
  // 5.1. Fishery catch-at-age to age-comp
  tc_obs.setZero();
  for (sp = 0; sp < nspp; sp++) {
    for (yr = 0; yr < nyrs_fsh_comp(sp); yr++) {
      for (ln = 0; ln < fsh_age_bins(sp); ln++) {
        tc_obs(sp, yr) += obs_catch(yr, ln, sp);
      }
      for (ln = 0; ln < fsh_age_bins(sp); ln++) {
        fsh_age_obs(sp, ln, yr) = obs_catch(yr, ln, sp) / (tc_obs(sp, yr) + 0.01); // Calculate age comp
      }
    }
  }

  // 5.2. BT survey CV estimation
  srv_biom_lse = srv_biom_se.array() / srv_biom.array();
  srv_biom_lse = pow( log( ( pow( srv_biom_lse.array(), Type(2) ).array() + 1).array()).array(), Type(0.5));


  // 5.3. EIT catch-at-age to age-comp
  for (yr = 0; yr < n_eit; yr++) {
    for (ln = 0; ln < srv_age_bins(0); ln++) {
      eit_age_comp(yr, ln) = obs_eit_age(yr, ln) / obs_eit_age.row(yr).sum(); // Convert from catch-at-age to age comp
    }
  }


  // 5.4. POPULATION DYNAMICS
  M1 = M1_base.array() + Type(0.0001);
  for ( sp = 0; sp < nspp ; sp++) {
    for ( age = 0 ; age < nages(sp); age++ ) {
      pmature( sp, age ) = pmature( sp, age ) * propMorF(sp + (mf_type(sp) - 1), age);
    }
  }

  // 5.5. Calculate temperature to use
  TempC = BTempC_retro.sum() / nTyrs; // Fill with average bottom temperature

  int yr_ind = 0;
  for (yr = 0; yr < nTyrs; yr++) {
    yr_ind = Tyrs( yr ) - styr;
    if ((yr_ind >= 0) & (yr_ind < nyrs)) {
      TempC(yr_ind) = BTempC_retro( yr );
    }
  }

  // 5.6. Calculate length-at-age
  for (sp = 0; sp < nspp; sp++) {
    for (age = 0; age < nages(sp); age++) {
      for (yr = 0; yr < nyrs; yr++) {
        LbyAge( sp, age, yr) = ( pow( ( 1 / aLW(0, sp) ), (1 / aLW(1, sp) ) ) ) * pow( ( wt(yr, age, sp) * 1000), (1 / aLW(1, sp))); // W = a L ^ b is the same as (W/a)^(1/b)
      }
    }
  }

  r_sigma = exp(ln_rec_sigma); // Convert log sd to natural scale

  // Dimensions of diet matrix
  vector<int> a1_dim = UobsAge.dim; // dimension of a1
  vector<int> a2_dim = UobsWtAge.dim; // dimension of a1


  // ------------------------------------------------------------------------- //
  // 6. POPULATION DYNAMICS EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // NOTE: Remember indexing starts at 0
  // Start iterations
  for (int iter = 0; iter < niter; iter++) {

    // 6.1. ESTIMATE FISHERY SELECTIVITY
    avgsel_fsh.setZero();
    fsh_sel.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nselages(sp); age++) {
        fsh_sel(sp, age) = fsh_sel_coff(sp, age);
        avgsel_fsh(sp) +=  exp(fsh_sel_coff(sp, age));
      }
      // 6.1.1. Average selectivity up to nselages
      avgsel_fsh(sp) = log(avgsel_fsh(sp) / nselages(sp));

      // 6.1.2. Plus group selectivity
      for (age = nselages(sp); age < nages(sp); age++) {
        fsh_sel(sp, age) = fsh_sel(sp, nselages(sp) - 1);
      }

      // 6.1.3. Average selectivity across all ages
      avgsel_tmp = 0; // Temporary object for average selectivity across all ages
      for (age = 0; age < nages(sp); age++) {
        avgsel_tmp += exp(fsh_sel(sp, age));
      }
      avgsel_tmp = log(avgsel_tmp / nages(sp));

      // 6.1.4. Standardize selectivity
      for (age = 0; age < nages(sp); age++) {
        fsh_sel(sp, age) -= avgsel_tmp;
        fsh_sel(sp, age) = exp(fsh_sel(sp, age));
      }

      vector<Type> sel_sp = fsh_sel.row(sp); // Can't max a matrix....
      fsh_sel.row(sp) /= max(sel_sp); // Standardize so max sel = 1 for each species

    }


    // 6.2. ESTIMATE FISHING MORTALITY
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs; yr++) {
        for (age = 0; age < nages(sp); age++) {
          F(sp, age, yr) = fsh_sel(sp, age) * exp(ln_mean_F(sp) + F_dev(sp, yr));
        }
      }
    }


    // 6.4. Estimate total mortality at age NOTE: May need to go above population dynamics
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          M(sp, age, yr) = M1(sp, age) + M2(sp, age, yr);
          Zed(sp, age, yr) = M1(sp, age) + F(sp, age, yr) + M2(sp, age, yr);
          S(sp, age, yr) = exp(-Zed(sp, age, yr));
        }
      }
    }


    // 6.1. ESTIMATE RECRUITMENT T1.1
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs; yr++) {
        R(sp, yr) = exp(ln_mn_rec(sp) + rec_dev(sp, yr));
        NByage(sp, 0, yr) = R(sp, yr);
      }
      //NByage(nyrs, 0, sp) = R(sp, nyrs-1);
    }


    // 6.2. ESTIMATE INITIAL ABUNDANCE AT AGE AND YEAR-1: T1.2
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        if ((age > 0) & (age < nages(sp) - 1)) {
          Type mort_sum = M1.row(sp).segment(0, age).sum(); // Sum M1 until age - 1
          NByage(sp, age, 0) = exp(ln_mn_rec(sp) - mort_sum + init_dev(sp, age - 1));
        }
        // -- 6.2.2. Where yr = 1 and age > Ai.
        if (age == (nages(sp) - 1)) {
          Type mort_sum = M1.row(sp).segment(0, age).sum(); // Sum M1 until age - 1
          NByage(sp, age, 0) = exp(ln_mn_rec(sp) - mort_sum + init_dev(sp, age - 1)) / (1 - exp(-M1(sp, nages(sp) - 1))); // NOTE: This solves for the geometric series
        }
      }
    }


    // 6.3. ESTIMATE NUMBERS AT AGE, BIOMASS-AT-AGE (kg), and SSB-AT-AGE (kg)
    biomass.setZero();
    biomassSSB.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 1; yr < nyrs; yr++) {
          // -- 6.3.1.  Where 1 <= age < Ai
          if (age < (nages(sp) - 1)) {
            NByage(sp, age + 1, yr) = NByage(sp, age, yr - 1) * S(sp, age, yr - 1);
          }

          // -- 6.3.2. Plus group where age > Ai. NOTE: This is not the same as T1.3 because sp used age = A_i rather than age > A_i.
          if (age == (nages(sp) - 1)) {
            NByage(sp, nages(sp) - 1, yr) = NByage(sp, nages(sp) - 2, yr - 1) * S(sp, nages(sp) - 2, yr - 1) + NByage(sp, nages(sp) - 1, yr - 1) * S(sp, nages(sp) - 1, yr - 1);
          }
        }

        // -- 6.3.3. Estimate Biomass and SSB
        for (yr = 0; yr < nyrs; yr++) {
          biomassByage(sp, age, yr) = NByage(sp, age, yr) * wt(yr, age, sp); // 6.5.
          biomassSSBByage(sp, age, yr) = biomassByage(sp, age, yr) * pmature(sp, age); // 6.6.

          biomass(sp, yr) += biomassByage(sp, age, yr);
          biomassSSB(sp, yr) += biomassSSBByage(sp, age, yr);
        }
      }
    }


    // 6.4. ESTIMATE AVERAGE NUMBERS AT AGE
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          if (avgnMode == 0) {
            AvgN(sp, age, yr) = NByage(sp, age, yr) * (1 - S(sp, age, yr)) / Zed(sp, age, yr); // MSVPA approach
          }
          if (avgnMode == 1) {
            AvgN(sp, age, yr) = NByage(sp, age, yr) * exp(- Zed(sp, age, yr) / 2); // Kinzey and Punt (2009) approximation
          }
          if (avgnMode == 2) {
            AvgN(sp, age, yr) = NByage(sp, age, yr); // Van Kirk et al (2010) approximation
          }
          // FIXME: put in break here
        }
      }
    }


    // ------------------------------------------------------------------------- //
    // 7. RATION EQUATIONS                                                       //
    // ------------------------------------------------------------------------- //
    // NOTE -- LOOPING INDICES -- sp = species, age = age, ln = length, ksp = prey, k_age = prey age (yr), k_ln = prey length, yr = year, rsp = predator, r_age = predator age (yr), r_ln = predator length

    // 7.1. Calculate stomach weight by sp age
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          S2Age(sp, age, yr) = S_a(0, sp) + (S_a(1, sp) * LbyAge(sp, age, yr)) + (S_a(2, sp) * pow(LbyAge(sp, age, yr), 2))
          + (S_a(3, sp) * pow(LbyAge(sp, age, yr), 3)) + (S_a(4, sp) * pow(LbyAge(sp, age, yr), 4)) + (S_a(5, sp) * pow(LbyAge(sp, age, yr), 5));

          if (LbyAge(sp, age, yr) > 80) {
            S2Age(sp, age, yr) = S_a(0, sp) + (S_a(1, sp) * 80) + (S_a(2, sp) * pow(80, 2)) + (S_a(3, sp) * pow(80, 3))
            + (S_a(4, sp) * pow(80, 4)) + (S_a(5, sp) * pow(80, 5)); // set everything above 80 to 80
          }
        }
      }
    }


    // 7.2. Calculate temperature function of consumption
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
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nTyrs; yr++) {
        if ( Ceq(sp) == 1) {
          fT(sp, yr) = exp(Qc(sp) * TempC(yr));
        }
        if ( Ceq(sp) == 2) {
          Yc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp) + 2);
          Zc = log( Qc(sp) ) * (Tcm(sp) - Tco(sp));
          Vc = (Tcm(sp) - TempC(yr)) / (Tcm(sp) - Tco(sp));
          Xc = pow(Zc, 2) * pow((1 + pow((1 + 40 / Yc), 0.5)), 2) / 400;
          fT(sp, yr) = pow(Vc, Xc) * exp(Xc * (1 - Vc));
        }
        if (Ceq(sp) == 3) {
          G2 = (1 / (Tcl(sp) - Tcm(sp))) * log((0.98 * (1 - CK4(sp))) / (CK4(sp) * 0.02));
          L2 = exp(G2 * (Tcl( sp ) - TempC( yr )));
          Kb = (CK4(sp) * L2) / (1 + CK4(sp) * (L2 - 1));
          G1 = (1 / (Tco(sp) - Qc(sp))) * log((0.98 * (1 - CK1(sp))) / (CK1(sp) * 0.02));
          L1 = exp(G1 * (TempC(yr) - Qc(sp)));
          Ka = (CK1(sp) * L1) / (1 + CK1(sp) * (L1 - 1));
          fT(sp, yr) = Ka * Kb;
        }
      }
    }


    // 7.3. Calculate historic ration
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs; yr++) {
          ConsumAge(sp, age, yr) = Type(24) * Type(0.0134) * exp( Type(0.0115) * TempC( yr )) * Type(91.25) * S2Age(sp, age, yr) * wt(yr, age, sp); // Calculate consumption for predator-at-age; units = kg/predator
          Consum_livingAge(sp, age, yr) = ConsumAge(sp, age, yr);       // Specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year

          if (C_model( sp ) == 1) {
            ConsumAge(sp, age, yr) = CA(sp) * pow(wt(yr, age, sp) * Type(1000), CB( sp )) * fT(sp, yr) * fday( sp ) * wt(yr, age, sp) * 1000;//g/pred.yr
            ConsumAge(sp, age, yr) = ConsumAge(sp, age, yr) * Pvalue(sp) * Pyrs(yr, age, sp); //
          }

          ration2Age(sp, age, yr) = ConsumAge(sp, age, yr) / 1000;      // Annual ration kg/yr //aLW(predd)*pow(lengths(predd,age),bLW(predd));//mnwt_bin(predd,age);
        }
      }
    }


    // 7.4. Calculate stomach content
    diet_w_sum.setZero();
    for (yr = 0; yr < nyrs; yr++) {                             // Year loop
      for (rsp = 0; rsp < nspp; rsp++) {                        // Predator species loop
        for (ksp = 0; ksp < nspp; ksp++) {                      // Prey species loop
          for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
            for (k_age = 0; k_age < nages(ksp); k_age++) {      // Prey age loop
              // Change to 5D array if 4D
              if(a1_dim.size() == 4){
                stomKir(rsp, ksp, r_age, k_age, yr) = UobsAge(rsp , ksp , r_age, k_age);
              }
              if(a1_dim.size() == 5){
                stomKir(rsp, ksp, r_age, k_age, yr) = UobsAge(rsp , ksp , r_age, k_age, yr);
              }

              if(a2_dim.size() == 4){
                stomKirWt(rsp, ksp, r_age, k_age, yr) = UobsWtAge(rsp, ksp, r_age, k_age);
              }
              if(a2_dim.size() == 5){
                stomKirWt(rsp, ksp, r_age, k_age, yr) = UobsWtAge(rsp, ksp, r_age, k_age, yr);
              }
            }
          }
          // By length
          for (r_ln = 0; r_ln < srv_age_bins(rsp); r_ln++) {      // Predator length loop
            for (k_ln = 0; k_ln < srv_age_bins(ksp); k_ln++) {    // Prey length loop
              diet_w_dat(rsp, ksp, r_ln, k_ln, yr) = UobsWt(rsp , ksp , r_ln, k_ln);
              diet_w_sum(rsp, ksp, r_ln, yr) += UobsWt(rsp , ksp , r_ln, k_ln);
            }
          }
        }
      }
    }


    // 7.5. Calculate other food stomach content
    of_stomKir.setZero();
    for (yr = 0; yr < nyrs; yr++) {                           // Year loop
      for (rsp = 0; rsp < nspp; rsp++) {                      // Predator species loop
        for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
          of_stomKir(rsp, r_age, yr) = Type( 1 );             // Initialize other suitability
          for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
            for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
              of_stomKir(rsp, r_age, yr) -= stomKir(rsp, ksp, r_age, k_age, yr);
            }
          }
        }
      }
    }
    of_stomKir /= other_food(0);


    // START PREDATION
    if (msmMode > 0) {

      // ------------------------------------------------------------------------- //
      // 8.1. SUITABILITY EQUATIONS                                                //
      // ------------------------------------------------------------------------- //
      // 8.1.1. Holsman and MSVPA based suitability // FIXME - not flexible for interannual variation
      if (suitMode == 0) {
        // 8.1.1.1. Calculate stomach proportion over biomass; U/ (W * N)
        suma_suit.setZero();
        for (yr = 0; yr < nyrs; yr++) {                         // Year loop
          for (ksp = 0; ksp < nspp; ksp++) {                    // Prey species loop
            for (rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop
              for (r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {  // Prey age loop
                  suit_tmp = stomKir(rsp, ksp, r_age, k_age, yr) / (AvgN(ksp, k_age, yr));
                  if (wt(yr, k_age, ksp) != 0) {
                    stom_div_bio2(rsp, ksp, r_age, k_age, yr) = suit_tmp / wt(yr, k_age, ksp);
                    suma_suit(rsp, r_age, yr ) += stom_div_bio2(rsp, ksp, r_age, k_age, yr); // Calculate sum of stom_div_bio2 across prey and  prey age for each predator, predator age, and year
                  }
                }
              }
            }
          }
        }


        // 8.1.1.2. Calculate suitability
        suit_main.setZero();
        suit_other.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                      // Predator species loop
          for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop
            suit_other(rsp, r_age) = Type( 1 );                 // Initialize other suitability
            for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
              for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, 0) += stom_div_bio2(rsp, ksp, r_age, k_age, yr) / ( suma_suit(rsp, r_age, yr ) + of_stomKir(rsp, r_age, yr) );
                }                                 // End year loop
                suit_main(rsp, ksp, r_age, k_age, 0) /= nyrs;
                suit_other(rsp, r_age) -= suit_main(rsp, ksp, r_age, k_age, 0); // Subtract observed suitability from entire suitability (ksp.e. 1)

                // Fill in years
                for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, yr) = suit_main(rsp, ksp, r_age, k_age, 0);
                }
              }
            }
          }
        }
      } // End Holsman/MSVPA suitability


      // 8.1.2. Length-based GAMMA suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 1){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
              ncnt = 0;
              gsum = 1.0e-10;                                           // Initialize
              for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                // if prey are smaller than predator:
                if (Mn_LatAge(rsp, r_age) > Mn_LatAge(ksp, k_age)) {
                  x_l_ratio = log(Mn_LatAge(rsp, r_age) / Mn_LatAge(ksp, k_age));
                  suit_main(rsp , ksp, r_age, k_age, 0) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                    (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                  ncnt += 1;
                  gsum += exp( suit_main(rsp , ksp, r_age, k_age, 0) );
                }
                else
                  suit_main(rsp , ksp, r_age, k_age, 0) = 0;
              }
              for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                // if prey are smaller than predator:
                if (Mn_LatAge(rsp, r_age) > Mn_LatAge(ksp, k_age)) {
                  suit_main(rsp , ksp, r_age, k_age, 0) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_age, k_age, 0) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                }

                // Fill in years
                for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, yr) = suit_main(rsp, ksp, r_age, k_age, 0);
                }
              }
            }
          }
        }
      } // End GAMMA selectivity


      // 8.1.2. Time-varying length-based GAMMA suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 2){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
              for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                ncnt = 0;
                gsum = 1.0e-10;                                           // Initialize
                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if ( LbyAge( rsp, r_age, yr)  > LbyAge( ksp, k_age, yr)) {
                    x_l_ratio = log(LbyAge( rsp, r_age, yr) / LbyAge( ksp, k_age, yr)); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                      (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                    ncnt += 1;
                    gsum += exp( suit_main(rsp , ksp, r_age, k_age, yr) );
                  }
                  else
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                }
                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if ( (aLW(0, rsp) * pow( wt(yr, r_age, rsp), aLW(1, rsp)))  > (aLW(0, ksp) * pow( wt(yr, k_age, ksp), aLW(1, ksp)))) {
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_age, k_age, yr) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                  }
                }
              }
            }
          }
        }
      } // End GAMMA selectivity

      // 8.1.3. Time-varying weight-based GAMMA suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 3){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type LenOpt = 0;          // Value of x_l_ratio where selectivity = 1
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                              // Pred loop
          LenOpt = Type(1.0e-10) + (gam_a(rsp) - 1) * gam_b(rsp);       // Eq. 18 Kinzey and Punt 2009
          for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
            for (r_age = 1; r_age < nages(rsp); r_age++) {              // Pred age // FIXME: start at 1?
              for (yr = 0; yr < nyrs; yr++) {                 // Year loop
                ncnt = 0;
                gsum = 1.0e-10;                                           // Initialize
                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if (wt(yr, r_age, rsp) >wt(yr, k_age, ksp)) {
                    x_l_ratio = log(wt(yr, r_age, rsp) / wt(yr, k_age, ksp)); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) +  (Type(1.0e-10) + gam_a( rsp ) - 1) * log(x_l_ratio / LenOpt + Type(1.0e-10)) -
                      (1.0e-10 + x_l_ratio - LenOpt) / gam_b(rsp);
                    ncnt += 1;
                    gsum += exp( suit_main(rsp , ksp, r_age, k_age, yr) );
                  }
                  else
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                }
                for (k_age = 0; k_age < nages(ksp); k_age++) {            // Prey age
                  // if prey are smaller than predator:
                  if (wt(yr, r_age, rsp) >wt(yr, k_age, ksp)) {
                    suit_main(rsp , ksp, r_age, k_age, yr) = Type(1.0e-10) + exp(suit_main(rsp , ksp, r_age, k_age, yr) - log(Type(1.0e-10) + gsum / Type(ncnt))); // NOT sure what this is for...
                  }
                }
              }
            }
          }
        }
      } // End GAMMA selectivity

      //-------------------------------------------------------------------------------------------------------------
      // 8.1.4. Length-based lognormal suitability // FIXME - not flexible for interannual variation in length-at-age
      // -- Turned off if not estimating selectivity
      if(suitMode == 4){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
          for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
            gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
            for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
              for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                // if prey are smaller than predator:
                if (Mn_LatAge(rsp, r_age) > Mn_LatAge(ksp, k_age)) {
                  x_l_ratio = log(Mn_LatAge(rsp, r_age) / Mn_LatAge(ksp, k_age)); // Log ratio of lengths
                  suit_main(rsp , ksp, r_age, k_age, 0) = exp(phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                  gsum += suit_main(rsp , ksp, r_age, k_age, 0);
                }
                else
                  suit_main(rsp , ksp, r_age, k_age, 0) = 0;
              }
            }
            for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
              for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                suit_main(rsp , ksp, r_age, k_age, 0) /= gsum;                // Scale, so it sums to 1.

                // Fill in years
                for (yr = 1; yr < nyrs; yr++) {                 // Year loop
                  suit_main(rsp, ksp, r_age, k_age, yr) = suit_main(rsp, ksp, r_age, k_age, 0);
                }
              }
            }
          }
        }
      } // End lognormal selectivity
      // 8.1.5. Time varying length-based lognormal suitability
      if(suitMode == 5){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
          for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
            for(yr = 0; yr < nyrs; yr++){
              gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  // if prey are smaller than predator:
                  if ( LbyAge( rsp, r_age, yr)  > LbyAge( ksp, k_age, yr)) {
                    x_l_ratio = log(LbyAge( rsp, r_age, yr) / LbyAge( ksp, k_age, yr)); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = exp(phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                    gsum += suit_main(rsp , ksp, r_age, k_age, yr);
                  }
                  else{
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                  }
                }
              }
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  suit_main(rsp , ksp, r_age, k_age, yr) /= gsum;                // Scale, so it sums to 1.
                }
              }
            }
          }
        }
      } // End lognormal selectivity


      // 8.1.6. Time varying weight-based lognormal suitability
      if(suitMode == 6){
        Type x_l_ratio = 0;       // Log(mean(predLen@age)/mean(preyLen@age))
        Type gsum = 0;

        suit_main.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                                // Pred loop
          for (r_age = 1; r_age < nages(rsp); r_age++) {                  // Pred age // FIXME: start at 1?
            for(yr = 0; yr < nyrs; yr++){
              gsum = 1;                                                     // Sum of suitability for each predator-at-age. Initialize at 1 because other biomass is assumed 1
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  // if prey are smaller than predator:
                  if (wt(yr, r_age, rsp) > wt(yr, k_age, ksp)) {
                    x_l_ratio = log(wt(yr, r_age, rsp) / wt(yr, k_age, ksp)); // Log ratio of lengths
                    suit_main(rsp , ksp, r_age, k_age, yr) = exp(phi(rsp, ksp)) * exp(-1/ (2*square(gam_a(rsp))) * square(x_l_ratio - gam_b(rsp)) );
                    gsum += suit_main(rsp , ksp, r_age, k_age, yr);
                  }
                  else{
                    suit_main(rsp , ksp, r_age, k_age, yr) = 0;
                  }
                }
              }
              for (ksp = 0; ksp < nspp; ksp++) {                            // Prey loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {              // Prey age
                  suit_main(rsp , ksp, r_age, k_age, yr) /= gsum;                // Scale, so it sums to 1.
                }
              }
            }
          }
        }
      } // End lognormal selectivity


      // ------------------------------------------------------------------------- //
      // 8. PREDATION MORTALITY EQUATIONS                                          //
      // ------------------------------------------------------------------------- //
      // -- 8.1. HOLSMAN PREDATION MORTALITY
      if (msmMode == 1) {
        // 8.1.3. Calculate available food
        avail_food.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {                        // Predator species loop
          for (r_age = 0; r_age < nages(rsp); r_age++) {          // Predator age loop
            for (yr = 0; yr < nyrs; yr++) {                       // Year loop
              tmp_othersuit = 0.;
              for (ksp = 0; ksp < nspp; ksp++) {                  // Prey species loop
                for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop
                  avail_food(rsp, r_age, yr) += suit_main(rsp, ksp, r_age, k_age, yr) * AvgN(ksp, k_age, yr) * wt(yr, k_age, ksp); // FIXME - include overlap indices: FIXME - mn_wt_stom?
                  tmp_othersuit += suit_main(rsp, ksp, r_age, k_age, yr); // FIXME - include overlap indices
                }
              }
              if(suitMode == 0){ // Holsman
                avail_food(rsp, r_age, yr) += other_food(rsp) * (Type(1) - (tmp_othersuit));
              }
              if(suitMode > 0){ // Log-normal and gamma
                avail_food(rsp, r_age, yr) += other_food(rsp) * Type(1);
              }
            }
          }
        }


        // 8.1.3. Calculate predation mortality
        M2.setZero();
        B_eaten.setZero();
        mn_UobsWtAge_hat.setZero();
        for (ksp = 0; ksp < nspp; ksp++) {                        // Prey species loop
          for (k_age = 0; k_age < nages(ksp); k_age++) {          // Prey age loop
            for (yr = 0; yr < nyrs; yr++) {                       // Year loop
              for (rsp = 0; rsp < nspp; rsp++) {                  // Predator species loop
                for (r_age = 0; r_age < nages(rsp); r_age++) {    // Predator age loop
                  M2(ksp, k_age, yr) += (AvgN(rsp, r_age, yr) * ration2Age(rsp, r_age, yr) * suit_main(rsp , ksp , r_age, k_age, yr)) / avail_food(rsp, r_age, yr); // #FIXME - include indices of overlap
                  B_eaten(ksp, k_age, yr) += AvgN(rsp, r_age, yr) * ration2Age(rsp, r_age, yr) * suit_main(rsp , ksp , r_age, k_age, yr);

                  // Estimated stomach proportion
                  UobsWtAge_hat(rsp, ksp, r_age, k_age, yr) = (AvgN(ksp, k_age, yr) * suit_main(rsp , ksp , r_age, k_age, yr) * wt(yr, k_age, ksp)) / avail_food(rsp, r_age, yr);
                  mn_UobsWtAge_hat(rsp, ksp, r_age, k_age) += UobsWtAge_hat(rsp, ksp, r_age, k_age, yr)/nyrs;
                }
              }
            }
          }
        }
      } // End 8..1. Holsman predation mortality


      // 8.2. KINZEY PREDATION EQUATIONS
      if (msmMode > 1) {
        // 8.2.2. Populate other food
        for (rsp = 0; rsp < nspp; rsp++) {
          for (r_age = 0; r_age < nages(rsp); r_age++) {
            Q_other_u(rsp, r_age) = other_food(rsp); // FIXME: may want to estimate other food exp(Q_other_est(rsp,r_age));
          }
        }


        // 8.2.3. Initialize counters
        Type Pred_ratio = 0;          // Predator ratio
        Type Prey_ratio = 0;          // Prey ratio
        Type pred_effect = 0;         // pred_resp * N(r,stage,yr)
        Type NS_Z = 0;                // N(k,yr,a) * survival/Z = 0;
        Type Tmort = 0;               // Mortality on other
        Type Q_ksum_l = 0;            // Diet sum
        Type Term = 0;                // Linear adjustment for predation


        // 8.2.4. Calculate equilibrium N predators and prey in styr_pred for each species X age: FIXME: May want to have this be the final year of a projection!
        N_pred_eq.setZero();
        N_prey_eq.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {
          for (ksp = 0; ksp < nspp; ksp++) {
            for (r_age = 0; r_age < nages(rsp); r_age++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                N_pred_eq(rsp, r_age) += NByage(rsp, r_age, 0) * suit_main(rsp , ksp, r_age, k_age, yr); // Denominator of Eq. 17 Kinzey and Punt (2009) 1st year
              }
            }
            for (k_age = 0; k_age < nages(ksp); k_age++) {
              for (r_age = 0; r_age < nages(rsp); r_age++) {
                N_prey_eq(ksp, k_age) += NByage(ksp, k_age, 0) * suit_main(rsp , ksp, r_age, k_age, yr); // Denominator of Eq. 16 Kinzey and Punt (2009) 1st year
              }
            }
          }
        }


        // 8.2.5. Calculate available prey and predator for each year
        N_pred_yrs.setZero();
        N_prey_yrs.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp < nspp; rsp++) {
            for (ksp = 0; ksp < nspp; ksp++) {
              for (r_age = 0; r_age < nages(rsp); r_age++) {
                for (k_age = 0; k_age < nages(ksp); k_age++) {
                  N_pred_yrs(rsp, r_age, yr) += NByage(rsp, r_age, yr) * suit_main(rsp , ksp, r_age, k_age, yr); // Numerator of Eq. 17 Kinzey and Punt (2009) 1st year // FIXME: Use averageN?
                }
              }
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                for (r_age = 0; r_age < nages(rsp); r_age++) {
                  N_prey_yrs(ksp, k_age, yr) += NByage(ksp, k_age, yr) * suit_main(rsp , ksp, r_age, k_age, yr); // Numerator of Eq. 16 Kinzey and Punt (2009) 1st year
                }
              }
            }
          }
        }


        // 8.2.6. Calculate predator functional response (Table 1 Kinzey and Punt (2009))
        for (rsp = 0; rsp < nspp; rsp++) {                          // Predator loop
          for (ksp = 0; ksp < (nspp + 1); ksp++) {                  // Prey loop
            for (yr = 0; yr < nyrs; yr++) {                         // Year loop
              for (r_age = 0; r_age < nages(rsp); r_age++) {        // Predator age loop

                Term = 1.0e-10 + H_1(rsp, ksp) * (Type(1) + H_1a(rsp) * H_1b(rsp) / (Type(r_age) + H_1b(rsp) + Type(1.0e-10))); // Eq. 15 Kinzey and Punt (2009)

                // Observed species
                if(ksp < nspp){
                  for (k_age = 0; k_age < nages(ksp); k_age++) {    // Prey age loop

                    // Predator-prey ratios
                    Pred_ratio = (N_pred_yrs(rsp, r_age, yr) + Type(1.0e-10)) / (N_pred_eq(rsp, r_age) + Type(1.0e-10)); // Eq. 17 Kinzey and Punt (2009): Predator biomass relative to equilibrium
                    Prey_ratio = (N_prey_yrs(ksp, k_age, yr) + Type(1.0e-10)) / (N_prey_eq(ksp, k_age) + Type(1.0e-10)); // Eq. 16 Prey Kinzey and Punt (2009): biomass relative to equilibrium
                    Pred_r(rsp, r_age, yr) = Pred_ratio;
                    Prey_r(ksp, k_age, yr) = Prey_ratio;

                    switch (msmMode) {
                    case 2: // Holling Type I (linear)
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term;
                      break;
                    case 3: // Holling Type II
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                        ( 1 + H_2(rsp, ksp) * Prey_ratio + 1.0e-10);
                      break;
                    case 4: // Holling Type III
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 +
                        Term * (1 + H_2(rsp, ksp)) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) /
                          (1 + H_2(rsp, ksp) * pow((Prey_ratio + 1.0e-10), H_4(rsp, ksp)) + 1.0e-10 );
                      break;
                    case 5: // predator interference
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                        ( 1 + H_2(rsp, ksp) * Prey_ratio + H_3(rsp, ksp) * (Pred_ratio - 1) + 1.0e-10);
                      break;
                    case 6: // predator preemption
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (1 + H_2(rsp, ksp) + Type(1.0e-10)) /
                        ( (1 + H_2(rsp, ksp) * Prey_ratio) * (1 + H_3(rsp, ksp) * (Pred_ratio - 1)) + Type(1.0e-10));
                      break;
                    case 7: // Hassell-Varley
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term * (2 + H_2(rsp, ksp) + 1.0e-10) /
                        (1.0 + H_2(rsp, ksp) * Prey_ratio + pow((Prey_ratio + Type(1.0e-10)), H_4(rsp, ksp)) + 1.0e-10 );
                      break;
                    case 8: // Ecosim
                      pred_resp(rsp, ksp, r_age, k_age, yr) = 1.0e-10 + Term /
                        (1 + H_3(rsp, ksp) * (Pred_ratio - 1 + 1.0e-10));
                      break;
                    default:
                      error("Invalid 'msmMode'");
                    }
                  }
                }
                // "other" is linear
                else{
                  pred_resp(rsp, ksp, r_age, 0, yr) = 1.0e-10 + Term;
                }
              } // end of r_ages, k_ages loop
            }   // =========================
          }
        }


        // 8.2.7. Mortality as a function of predator age: Eq. 7 Kinzey and Punt (2009)
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for (ksp = 0; ksp  <  nspp; ksp++) {
              for (r_age = 0; r_age < nages(rsp); r_age++) {
                for (k_age = 0; k_age < nages(ksp); k_age++) {
                  pred_effect = pred_resp(rsp, ksp, r_age, k_age, yr) * suit_main(rsp , ksp, r_age, k_age, yr);
                  Vmort_ua(rsp, ksp, r_age, k_age, yr) = pred_effect * NByage(rsp, r_age, yr);
                }
              }
            }
          }
        }


        // 8.2.8 Accumulate predation mortality: Eq. 6 Kinzey and Punt (2009)
        Pmort_ua.setZero();
        M2.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp < nspp; rsp++) {
            for (ksp = 0; ksp < nspp; ksp++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                for (r_age = 0; r_age < nages(rsp); r_age++) {
                  Pmort_ua(rsp, ksp, k_age, yr) += Vmort_ua(rsp, ksp, r_age, k_age, yr);
                }
                M2(ksp, k_age, yr) += Pmort_ua(rsp, ksp, k_age, yr);
              }
            }
          }
        }


        // 8.2.9. Numbers eaten (of modeled prey species); Equations 8 and 9 from Kinzey and Punt 2009
        eaten_la.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (ksp = 0; ksp  <  nspp; ksp++) {
            for (k_age = 0; k_age < nages(ksp); k_age++) {
              // Relative number
              NS_Z = NByage(ksp, k_age, yr) * (1 - exp(-Zed(ksp, k_age, yr) )) / Zed(ksp, k_age, yr);                                  // Eq. 8 Kinzey and Punt (2009) Baranov

              for (rsp = 0; rsp  <  nspp; rsp++) {
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  eaten_ua(rsp, ksp, r_age, k_age, yr) = Vmort_ua(rsp, ksp, r_age, k_age, yr) * NS_Z;                                  // Eq. 8 Kinzey and Punt (2009)

                  for (r_ln = 0; r_ln  <  srv_age_bins(rsp); r_ln++) {
                    eaten_la(rsp, ksp, r_ln, k_age, yr) += eaten_ua(rsp, ksp, r_age, k_age, yr)  * age_trans_matrix(r_age, r_ln, rsp); // Eq. 9 Kinzey and Punt (2009)
                  }
                }
              }
            }
          }
        }


        // 8.2.10. Mass eaten (including "other")
        Q_mass_l.setZero();
        Q_mass_u.setZero();
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  < nspp; rsp++) {
            for (ksp = 0; ksp  <  nspp + 1; ksp++) {

              // Species included
              if (ksp < nspp) {
                // Results by length: Eqn 10 kinda Kinzey and Punt (2009)
                for (r_ln = 0; r_ln  <  srv_age_bins(rsp); r_ln++) {
                  for (k_age = 0; k_age < nages(ksp); k_age++) {
                    Q_mass_l(rsp, ksp, r_ln, yr) += eaten_la(rsp, ksp, r_ln, k_age, yr) * wt(yr, k_age, ksp);
                  }
                }
                // Results by k_age: Eq. 11 Kinzey and Punt (2009)
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  for (k_age = 0; k_age < nages(ksp); k_age++) {
                    Q_mass_u(rsp, ksp, r_age, yr) += eaten_ua(rsp, ksp, r_age, k_age, yr) * wt(yr, k_age, ksp);
                  }
                }
              }

              // Other food
              else {
                for (r_age = 0; r_age  <  nages(rsp); r_age++) {
                  pred_effect = pred_resp(rsp, ksp, r_age, 0, yr); // Other food k_age = 0 because 1 age is used
                  Tmort = pred_effect * NByage(rsp, r_age, yr);
                  Q_mass_u(rsp, ksp, r_age, yr)  = Q_other_u(rsp, r_age) * (Type(1) - exp(-Tmort));                     // Eq. 11 Kinzey and Punt (2009)
                }
                for (r_ln = 0; r_ln < srv_age_bins(rsp); r_ln++) {
                  for (r_age = 0; r_age < nages(rsp); r_age++) {
                    Q_mass_l(rsp, ksp, r_ln, yr) += Q_mass_u(rsp, ksp, r_age, yr) * age_trans_matrix(r_age, r_ln, rsp); // Eq. 10 Kinzey and Punt (2009)
                  }
                }
              }
            }
          }
        }


        // 8.2.11. Total up the consumption by each predator and normalize (Eqn 15)
        for (yr = 0; yr < nyrs; yr++) {
          for (rsp = 0; rsp  <  nspp; rsp++) {
            for (r_ln = 0; r_ln < srv_age_bins(rsp); r_ln++) {
              Q_ksum_l = 0;
              for (ksp = 0; ksp < (nspp + 1); ksp++) {
                Q_ksum_l += Q_mass_l(rsp, ksp, r_ln, yr) + 1.0e-10; // Table 3 - Diet weight proportions Kinzey and Punt (2009)
              }
              for (ksp = 0; ksp < (nspp + 1); ksp++) {
                Q_hat(rsp, ksp, r_ln, yr) = (1.0e-10 + Q_mass_l(rsp, ksp, r_ln, yr) / Q_ksum_l); // changed paranthesis -dhk apr 28 09
              }
            }
          }
        }


        // 8.2.12. Predict ration Eq. 13 Kinzey and Punt (2009)
        Type n_avg, numer, denom;

        omega_hat_ave.setZero();
        omega_hat.setZero();
        for (rsp = 0; rsp < nspp; rsp++) {
          for (r_age = 0; r_age < nages(rsp); r_age++) {
            // Calculate year-specific values
            numer = 0;
            denom = 0;
            for (yr = 0; yr < nyrs; yr++) {
              // Average abundance
              n_avg = 1.0e-10 + NByage(rsp, r_age, yr) * sqrt(S(rsp, r_age, yr));  // FIXME: Add average N added 1.0e-10 -dhk June 24 08. not in NRM tpl -dhk apr 28 09
              // find total consumption by this r_age-class
              for (ksp = 0; ksp < (nspp + 1); ksp++) {
                omega_hat(rsp, r_age, yr) += Q_mass_u(rsp, ksp, r_age, yr);
              }
              numer += omega_hat(rsp, r_age, yr) / + Type(1.0e-10); // FIXME: Divide by 365 to make into daily ration
              denom += n_avg;
              // normalize
              omega_hat(rsp, r_age, yr) /= (n_avg); // FIXME: Divide by 365 to make into daily ration
            }
            omega_hat_ave(rsp, r_age) = numer / denom;
          }
        }


      } // End 8.2. Kinzey predation
    } // End 8. Predation mortality


    // ------------------------------------------------------------------------- //
    // 9. SURVEY COMPONENTS EQUATIONS                                            //
    // ------------------------------------------------------------------------- //
    // 9.1. Survey selectivity
    avgsel_srv.setZero();
    srv_sel.setZero();
    for (sp = 0; sp < nspp; sp++) {
      // 9.1.1. Logisitic selectivity
      if (srv_sel_type(sp) == 0) {
        for (age = 0; age < nages(sp); age++){
          srv_sel(sp, age) = 1 / (1 + exp( -srv_sel_slp(0, sp) * ((age + 1) - srv_sel_inf(0, sp))));
        }
      }

      // 9.1.2. Non-parametric selectivity fit to age ranges. NOTE: This can likely be improved
      if (srv_sel_type(sp) == 1) {
        for (age = 0; age < nselages(sp); age++) {
          srv_sel(sp, age) = srv_sel_coff(sp, age);
          avgsel_srv(sp) +=  exp(srv_sel_coff(sp, age));
        }
        // 9.1.3 Average selectivity up to nselages(sp)
        avgsel_srv(sp) = log(avgsel_srv(sp) / nselages(sp));

        // 9.1.4. Plus group selectivity
        for (age = nselages(sp); age < nages(sp); age++) {
          srv_sel(sp, age) = srv_sel(sp, nselages(sp) - 1);
        }

        // 9.1.5. Average selectivity across all ages
        avgsel_tmp = 0; // set to zero
        for (age = 0; age < nages(sp); age++) {
          avgsel_tmp += exp(srv_sel(sp, age));
        }
        avgsel_tmp = log(avgsel_tmp / nages(sp));

        // 9.1.6. Standardize selectivity
        for (age = 0; age < nages(sp); age++) {
          srv_sel(sp, age) -=  avgsel_tmp;
          srv_sel(sp, age) = exp(srv_sel(sp, age));
        }


      }

      // 9.1.2. Double logistic (Dorn and Methot 1990)
      if (srv_sel_type(sp) == 2) {
        for (age = 0; age < nages(sp); age++){
          srv_sel(sp, age) = (1 / (1 + exp( -srv_sel_slp(0, sp) * ((age + 1) - srv_sel_inf(0, sp))))) * // Upper slop
            (1 - (1 / (1 + exp( -srv_sel_slp(1, sp) * ((age + 1) - srv_sel_inf(1, sp))))));  // Downward slope;
        }
      }
      vector<Type> sel_sp = srv_sel.row(sp); // Can't max a matrix....
      srv_sel.row(sp) /= max(sel_sp); // Standardize so max sel = 1 for each species
    } // End species loop

    // 9.2 EIT Components
    // -- 9.2.1 EIT Survey Biomass
    eit_hat.setZero();
    for (age = 0; age < nages(0); age++) {
      for (yr = 0; yr < n_eit; yr++) {
        eit_yr = yrs_eit(yr) - styr;
        eit_age_hat(age, yr) = NByage(0, age, eit_yr) * eit_sel(eit_yr, age) * eit_q; // Remove the mid-year trawl?
        eit_hat(yr) += eit_age_hat(age, yr) * wt(eit_yr, age, 0);  //
      }
    }


    // -- 9.2.2 EIT Survey Age Composition
    for (age = 0; age < nages(0); age++) {
      for (yr = 0; yr < n_eit; yr++) {
        eit_age_comp_hat(age, yr) = eit_age_hat(age, yr) / eit_age_hat.col(yr).sum(); // Divide numbers at age by total numbers for each year
      }
    }


    // 9.3 BT Components
    // -- 9.3.1 BT Survey Biomass
    srv_bio_hat.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs_srv_biom(sp); yr++) {
          srv_yr = yrs_srv_biom(sp, yr) - styr; // Temporary index for years of data
          srv_bio_hat(sp, yr) += NByage(sp, age, srv_yr) * exp(-0.5 * Zed(sp, age, srv_yr)) * srv_sel(sp, age) * exp(log_srv_q(sp)) * wt(srv_yr, age, sp);  //
        }
      }
    }

    // -- 9.4.2 BT Survey Age Composition: NOTE: will have to alter if age comp data are not the same length as survey biomass data
    srv_hat.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs_srv_age(sp); yr++) {
          srv_yr = yrs_srv_age(sp, yr) - styr; // Temporary index for years of data
          srv_age_hat(sp, age, yr) = NByage(sp, age, srv_yr) * srv_sel(sp, age) * exp(log_srv_q(sp));
          srv_hat(sp, yr) += srv_age_hat(sp, age, yr);   // Total numbers
        }
      }
    }


    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs_srv_age(sp); yr++) {
        // 9.4.2.1 -- BT Survey catch-at-age
        if (srv_age_type(sp) == 1) {
          for (age = 0; age < nages(sp); age++) {
            srv_age_hat(sp, age, yr) = srv_age_hat(sp, age, yr) / srv_hat(sp, yr);
          }
        }
        // 9.4.2.2 -- Convert from catch-at-age to catch-at-length: NOTE: There has got to be a better way
        if (srv_age_type(sp) != 1) {
          for (age = 0; age < nages(sp); age++) {
            srv_age_tmp(age) = srv_age_hat(sp, age, yr);
          }

          matrix<Type> ALK = trim_matrix( matrix_from_array(age_trans_matrix, sp), nages(sp), srv_age_bins(sp) );
          vector<Type> srv_age_tmp_trimmed = srv_age_tmp.segment(0, nages(sp));
          vector<Type> srv_len_tmp = vec_mat_prod( srv_age_tmp_trimmed , ALK ); // Multiply the ALK for species sp against the survey catch-at-age for year yr

          for (ln = 0; ln < srv_age_bins(sp); ln++) {
            srv_age_hat(sp, ln, yr) = srv_len_tmp(ln) / srv_len_tmp.sum() ; // * age_trans_matrix.col().col(sp)) / srv_hat(sp, yr); // # NOTE: Double check the matrix algebra here
          }
        }
      }
    }


    // ------------------------------------------------------------------------- //
    // 10. FISHERY COMPONENTS EQUATIONS                                           //
    // ------------------------------------------------------------------------- //
    // 10.5. ESTIMATE CATCH-AT-AGE and TOTAL YIELD (kg)
    tc_biom_hat.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs_tc_biom(sp); yr++) {
          fsh_yr = yrs_tc_biom(sp, yr) - styr; // Temporary index for years of data
          tc_biom_hat(sp, yr) += F(sp, age, fsh_yr) / Zed(sp, age, fsh_yr) * (1 - exp(-Zed(sp, age, fsh_yr))) * NByage(sp, age, fsh_yr) * wt(yr, age, sp); // 5.5.
        }
      }
    }


    // 10.6. ESTIMATE FISHERY AGE COMPOSITION
    // 10.6.1. Get catch-at-age
    tc_hat.setZero();
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        for (yr = 0; yr < nyrs_tc_biom(sp); yr++) {
          fsh_yr = yrs_tc_biom(sp, yr) - styr; // Temporary index for years of data
          catch_hat(sp, age, yr) = F(sp, age, fsh_yr) / Zed(sp, age, fsh_yr) * (1 - exp(-Zed(sp, age, fsh_yr))) * NByage(sp, age, fsh_yr); // 5.4.
          tc_hat(sp, yr) += catch_hat(sp, age, yr); // Estimate catch in numbers
        }
      }
    }

    // 10.6.2 Convert catch-at-age to age-comp
    for (sp = 0; sp < nspp; sp++) {
      for (yr = 0; yr < nyrs; yr++) {
        /// 10.7.2.1 -- Estimate age composition of the fishery
        if (fsh_age_type(sp) == 1) {
          for (age = 0; age < nages(sp); age++) {
            fsh_age_hat(sp, age, yr) = catch_hat(sp, age, yr) / tc_hat(sp, yr);
          }
        }

        // 10.7.2.1 -- Convert from catch-at-age to catch-at-length: NOTE: There has got to be a better way
        if (fsh_age_type(sp) != 1) {
          for (age = 0; age < nages(sp); age++) {
            fsh_age_tmp(age) = catch_hat(sp, age, yr);
          }

          matrix<Type> ALK = trim_matrix( matrix_from_array(age_trans_matrix, sp), nages(sp), fsh_age_bins(sp) );
          vector<Type> fsh_age_tmp_trimmed = fsh_age_tmp.segment(0, nages(sp) );
          vector<Type> fsh_len_tmp = vec_mat_prod( fsh_age_tmp_trimmed , ALK ); // Multiply the ALK for species sp against the survey catch-at-age for year yr

          for (ln = 0; ln < fsh_age_bins(sp); ln++) {
            fsh_age_hat(sp, ln, yr) = fsh_len_tmp(ln) / fsh_len_tmp.sum() ; // FIXME: because the age_trans_matrix rows do not sum to 1 for PCod, we are underestimating length_comp
          }
        }
      }
    }
    // End iterations
  }

  // ------------------------------------------------------------------------- //
  // 11. LIKELIHOOD EQUATIONS                                                   //
  // ------------------------------------------------------------------------- //
  // 11.0. OBJECTIVE FUNCTION
  matrix<Type> jnll_comp(17, nspp); jnll_comp.setZero();  // matrix of negative log-likelihood components

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
  // -- M2 likelihood components
  // Slot 13 -- Ration likelihood
  // Slot 14 -- Ration penalties
  // Slot 15 -- Diet weight likelihood
  // Slot 16 -- Stomach content of prey length ln in predator length a likelihood


  // 11.1. OFFSETS AND PENALTIES
  // 11.1.1 -- Set up offset objects
  vector<Type> offset_srv(nspp); offset_srv.setZero(); // Offset for multinomial likelihood
  vector<Type> offset_fsh(nspp); offset_fsh.setZero(); // Offset for multinomial likelihood
  vector<Type> offset_diet_w(nspp); offset_diet_w.setZero(); // Offset for total stomach content weight likelihood
  vector<Type> offset_diet_l(nspp); offset_diet_l.setZero(); // Offset for total stomach content of prey length ln in predator length a
  vector<Type> offset_uobsagewt(nspp); offset_uobsagewt.setZero(); // Offset for stomach proportion by weight likelihood
  Type offset_eit = 0; // Offset for multinomial likelihood


  // 11.1.1. -- Fishery age comp offsets
  for (sp = 0; sp < nspp; sp++) {
    for (yr = 0; yr < nyrs_fsh_comp(sp); yr++) {
      for (ln = 0; ln < fsh_age_bins(sp); ln++) {
        offset_fsh( sp ) -= tau * (fsh_age_obs(sp, ln, yr) + Type(0.001)) * log(fsh_age_obs(sp, ln, yr) + Type(0.001));
      }
    }
  }

  // 11.1.2. -- Survey age comp offsets
  for (sp = 0; sp < nspp; sp++) {
    matrix<Type> srv_age_tmp = trim_matrix(matrix_from_array(srv_age_obs, sp), nyrs_srv_age(sp), srv_age_bins(sp));
    for (yr = 0; yr < nyrs_srv_age(sp); yr++) {
      for (ln = 0; ln < srv_age_bins(sp); ln++) {
        srv_age_obs(yr, ln, sp) = srv_age_obs(yr, ln, sp) / srv_age_tmp.row(yr).sum() ; // Convert numbers-at-age to age comp
        offset_srv(sp) -= srv_age_n(sp, yr) * (srv_age_obs(yr, ln, sp) + MNConst) * log(srv_age_obs(yr, ln, sp) + MNConst ) ;
      }
    }
  }

  // 11.1.3. -- Offsets for acoustic survey
  for (yr = 0; yr < n_eit; yr++) {
    for (age = 0; age < nages(0); age++) {
      offset_eit -= eit_age_n(yr) * (eit_age_comp(yr, age) + MNConst) * log(eit_age_comp(yr, age) + MNConst );
    }
  }

  // 11.1.4. -- Offsets for diet (weights)
  for (rsp = 0; rsp < nspp; rsp++) {
    for (yr = 0; yr < nyrs; yr++) {
      for (r_ln = 0; r_ln < srv_age_bins(rsp); r_ln++) {
        // if (stoms_w_N(rsp,r_ln,yr) > 0){ Sample size
        for (ksp = 0; ksp < nspp; ksp++) { // FIXME: no offset for "other food"
          if ( diet_w_sum(rsp, ksp, r_ln, yr) > 0) {
            offset_diet_w(rsp) -= stom_tau * (diet_w_sum(rsp, ksp, r_ln, yr) + MNConst) * log(diet_w_sum(rsp, ksp, r_ln, yr) + MNConst); // FIXME add actual sample size
          }
        }
      }
    }
  }

  // 11.1.5. -- Offsets for diet (lengths) FIXME: add real stomach sample size
  for (rsp = 0; rsp < nspp; rsp++) {
    for (ksp = 0; ksp < nspp; ksp++) {
      for (r_ln = 0; r_ln < srv_age_bins(rsp); r_ln++) {
        for (k_ln = 0; k_ln < srv_age_bins(ksp); k_ln++) {
          if (Uobs(rsp, ksp, r_ln, k_ln) > 0) {
            offset_diet_l(rsp) -= stom_tau * Uobs(rsp, ksp, r_ln, k_ln) * log(Uobs(rsp, ksp, r_ln, k_ln) + 1.0e-10);
          }
        }
      }
    }
  }

  // 11.1.6. -- Offsets for diet proportion by weight
  // If 5D array
  if(a2_dim.size() == 5){
    for (rsp = 0; rsp < nspp; rsp++) {
      for (yr = 0; yr < nyrs; yr++) {
        for (r_age = 0; r_age < nages(rsp); r_age++) {
          for (ksp = 0; ksp < nspp; ksp++) {
            for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
              //if (UobsWtAge(rsp, ksp, r_age, k_age, yr) > 0) { // (rsp, ksp, a, r_ln, yr)
                offset_uobsagewt(rsp) -= stom_tau * (UobsWtAge(rsp, ksp, r_age, k_age, yr) + MNConst) * log(UobsWtAge(rsp, ksp, r_age, k_age, yr) + MNConst);
              //}
            }
          }
        }
      }
    }
  }

  // If 4D array
  if(a2_dim.size() == 4){
    for (rsp = 0; rsp < nspp; rsp++) {
      for (r_age = 0; r_age < nages(rsp); r_age++) {
        for (ksp = 0; ksp < nspp; ksp++) {
          for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
            //if (UobsWtAge(rsp, ksp, r_age, k_age) > 0) { // (rsp, ksp, a, r_ln, yr)
              offset_uobsagewt(rsp) -= stom_tau * (UobsWtAge(rsp, ksp, r_age, k_age) + MNConst) * log(UobsWtAge(rsp, ksp, r_age, k_age) + MNConst);
            //}
          }
        }
      }
    }
  }




  // 11.2. FIT OBJECTIVE FUNCTION
  // Slot 0 -- BT survey biomass -- NFMS annual BT survey
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(0, sp) = 0; // FIXME: Likeliy redundant
    for (yr = 0; yr < nyrs_srv_biom(sp); yr++) {
      // srv_yr = yrs_srv_biom(sp, yr) - styr;
      jnll_comp(0, sp) += pow(log(srv_biom(sp, yr)) - log(srv_bio_hat(sp, yr)), 2) / (2 * pow(srv_biom_lse(sp, yr), 2)); // NOTE: This is not quite the lognormal and biohat will be the median.
    }
  }


  // Slot 1 -- BT survey age composition -- NFMS annual BT survey
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(1, sp) = 0; // FIXME: Likeliy redundant
    for (yr = 0; yr < nyrs_srv_age(sp); yr++) {
      for (ln = 0; ln < srv_age_bins(sp); ln++) {
        // srv_yr = yrs_srv_age(sp, yr) - styr;
        jnll_comp(1, sp) -= srv_age_n(sp, yr) * (srv_age_obs(yr, ln, sp) + MNConst) * log(srv_age_hat(sp, ln, yr) + MNConst);
      }
    }
    jnll_comp(1, sp) -= offset_srv(sp);
  }


  // Slot 2 -- EIT survey biomass -- Pollock acoustic trawl survey
  jnll_comp(2, 0) = 0; // FIXME: Likeliy redundant
  for (yr = 0; yr < n_eit; yr++) {
    jnll_comp(2, 0) += 12.5 * pow(log(obs_eit(yr)) - log(eit_hat(yr) + 1.e-04), 2);
  }


  // Slot 3 -- EIT age composition -- Pollock acoustic trawl survey
  jnll_comp(3, 0) = 0; // FIXME: Likeliy redundant
  for (yr = 0; yr < n_eit; yr++) {
    for (age = 0; age < nages(0); age++) {
      // eit_yr = yrs_eit(yr) - styr;
      jnll_comp(3, 0) -= eit_age_n(yr) *  (eit_age_comp(yr, age) + MNConst) * log(eit_age_comp_hat(age, yr) + MNConst);
    }
  }
  jnll_comp(3, 0) -= offset_eit;


  // Slot 4 -- Total catch -- Fishery observer data
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(4, sp) = 0; // FIXME: Likeliy redundant
    for (yr = 0; yr < nyrs_tc_biom(sp); yr++) {
      fsh_yr = yrs_tc_biom(sp, yr) - styr;
      jnll_comp(4, sp) += pow((log(tcb_obs(sp, yr) + Type(1.e-4)) - log(tc_biom_hat(sp, yr) + Type(1.e-4))), 2) / (2 * pow(sigma_catch, 2)); // T.4.5
    }
  }


  // Slot 5 -- Fishery age composition -- Fishery observer data
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(5, sp) = 0; // FIXME: Likeliy redundant
    for (yr = 0; yr < nyrs_fsh_comp(sp); yr++) {
      for (ln = 0; ln < fsh_age_bins(sp); ln++) {
        fsh_yr = yrs_fsh_comp(sp, yr) - styr; // Temporary index for years of data
        jnll_comp(5, sp) -= tau * (fsh_age_obs(sp, ln, yr) + MNConst) * log(fsh_age_hat(sp, ln, fsh_yr) + MNConst );
      }
    }
    jnll_comp(5, sp) -= offset_fsh(sp);
  }


  // Slot 6 -- Fishery selectivity
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(6, sp) = 0; // FIXME: Likeliy redundant
    for (age = 0; age < (nages(sp) - 1); age++) {
      if ( fsh_sel(sp, age) > fsh_sel( sp, age + 1 )) {
        jnll_comp(6, sp) += 20 * pow( log( fsh_sel(sp, age) / fsh_sel(sp, age + 1 ) ), 2);
      }
    }

    // Extract only the selectivities we want
    vector<Type> sel_tmp(nages(sp)); sel_tmp.setZero();
    for (age = 0; age < nages(sp); age++) {
      sel_tmp(age) = log(fsh_sel(sp, age));
    }

    for (age = 0; age < nages(sp) - 2; age++) {
      sel_tmp(age) = first_difference( first_difference( sel_tmp ) )(age);
      jnll_comp(6, sp) += curv_pen_fsh * pow( sel_tmp(age) , 2); // FIX

    }
  }


  // Slot 6 -- Add fishery selectivity normalization
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(7, sp) = 0; // FIXME: Likeliy redundant
    jnll_comp(7, sp) += 50 * pow(avgsel_fsh(sp), 2);
  }



  // Slot 7 -- Survey selectivity
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(8, sp) = 0; // FIXME: Likeliy redundant
    if (srv_sel_type(sp) == 1) {
      // Extract only the selectivities we want
      vector<Type> sel_tmp(nages(sp));
      for (age = 0; age < nages(sp); age++) {
        sel_tmp(age) = log(srv_sel(sp, age));
      }
      for (age = 0; age < nages(sp) - 2; age++) {
        sel_tmp(age) = first_difference( first_difference( sel_tmp ) )(age);
        jnll_comp(8, sp) += curv_pen_fsh * pow( sel_tmp(age) , 2); // FIX
      }
    }
  }


  // Slot 7 -- Add survey selectivity normalization
  for (sp = 0; sp < nspp; sp++) {
    jnll_comp(9, sp) = 0; // FIXME: Likeliy redundant
    jnll_comp(9, sp) += 50 * pow(avgsel_srv(sp), 2);
  }


  // Slots 10-12 -- PRIORS: PUT RANDOM EFFECTS SWITCH HERE
  for (sp = 0; sp < nspp; sp++) {
    // Slot 10 -- init_dev -- Initial abundance-at-age
    for (age = 1; age < nages(sp); age++) {

      if (random_rec == 0) {
        jnll_comp(10, sp) += pow( init_dev(sp, age - 1) - Type(0.25), 2);
      }

      if (random_rec == 1) {
        jnll_comp(10, sp) -= dnorm( init_dev(sp, age - 1) - square(r_sigma(sp)) / 2, Type(0.0), r_sigma(sp), true);
      }

    }

    for (yr = 0; yr < nyrs; yr++) {

      // Slot 11 -- Tau -- Annual recruitment deviation
      if (random_rec == 0) {
        jnll_comp(11, sp) += pow( rec_dev(sp, yr) - Type(0.25), 2);    // Recruitment deviation using penalized likelihood. ADDED lognormal bias correction
      }
      if (random_rec == 1) {
        jnll_comp(11, sp) -= dnorm( rec_dev(sp, yr) - square(r_sigma(sp)) / 2, Type(0.0), r_sigma(sp), true);    // Recruitment deviation using random effects.
      }

      // Slot 12 -- Epsilon -- Annual fishing mortality deviation
      jnll_comp(12, sp) += pow( F_dev(sp, yr), 2);      // Fishing mortality deviation using penalized likelihood.
    }
  }


  // 11.3. Diet likelihood components from MSVPA
  if ((msmMode == 1) & (suitMode > 0)) {
    // Slot 14 -- Diet weight likelihood
    // If 5D array
    if(a2_dim.size() == 5){
      for (rsp = 0; rsp < nspp; rsp++) {
        jnll_comp(15, rsp) = 0;
        for (yr = 0; yr < nyrs; yr++) {
          for (r_age = 0; r_age < nages(rsp); r_age++) {
            for (ksp = 0; ksp < nspp; ksp++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
                //if (UobsWtAge(rsp, ksp, r_age, k_age, yr) > 0) { // (rsp, ksp, a, r_ln, yr)
                  jnll_comp(15, rsp) -= stom_tau * (UobsWtAge(rsp, ksp, r_age, k_age, yr) + MNConst) * (log(UobsWtAge_hat(rsp, ksp, r_age, k_age, yr) + MNConst));
               // }
              }
            }
          }
        }
        jnll_comp(15, rsp) -= offset_uobsagewt(rsp);
      }
    }

    // If 4D array
    if(a2_dim.size() == 4){
      for (rsp = 0; rsp < nspp; rsp++) {
                jnll_comp(15, rsp) = 0;
        for (r_age = 0; r_age < nages(rsp); r_age++) {
          for (ksp = 0; ksp < nspp; ksp++) {
            for (k_age = 0; k_age < nages(ksp); k_age++) {                  // FIME: need to add in other food
              //if (UobsWtAge(rsp, ksp, r_age, k_age) > 0) { // (rsp, ksp, a, r_ln, yr)
                jnll_comp(15, rsp) -= stom_tau * (UobsWtAge(rsp, ksp, r_age, k_age) + MNConst) * (log(mn_UobsWtAge_hat(rsp, ksp, r_age, k_age) + MNConst));
              //}
            }
          }
        }
        jnll_comp(15, rsp) -= offset_uobsagewt(rsp);
      }
    }
  } // End diet proportion by weight component

  // 11.4. Diet likelihood components from Kinzey and Punt
  if ((msmMode > 1)) {
    // Slot 13 -- Ration likelihood
    for (yr = 0; yr < nyrs; yr++) {
      for (sp = 0; sp < nspp; sp++) {
        for (age = 1; age < nages(sp); age++) { // don't include age zero in likelihood
          jnll_comp(13, sp) += 0.5 * pow( log( omega_hat(sp, age, yr) + 1.0e-10) -
            log(ration2Age(sp, age, yr)), 2) / (sd_ration * sd_ration); // FIXME: add year indices for ration
        }
      }
    }


    // Slot 14 -- Ration penalalties FIXME: not sure if necessary: talk to Andre
    for (sp = 0; sp < nspp; sp++) {
      for (age = 0; age < nages(sp); age++) {
        Type mean_ohat = 0;
        for (yr = 0; yr < nyrs; yr++) {
          mean_ohat += omega_hat(sp, age, yr) / nyrs;
        }
        for (yr = 0; yr < nyrs; yr++) {
          jnll_comp(14, sp) += 20 *  pow(omega_hat(sp, age, yr) - mean_ohat, 2);
        }
      }
    }


    // Slot 14 -- Diet weight likelihood
    for (rsp = 0; rsp < nspp; rsp++) {
      for (yr = 0; yr < nyrs; yr++) {
        for (r_ln = 0; r_ln < srv_age_bins(rsp); r_ln++) {
          for (ksp = 0; ksp < nspp; ksp++) {                  // FIME: need to add in other food
            if (diet_w_sum(rsp, ksp, r_ln, yr) > 0) { // (rsp, ksp, a, r_ln, yr)
              jnll_comp(15, rsp) -= stom_tau * diet_w_sum(rsp, ksp, r_ln, yr) * log(Q_hat(rsp, ksp, r_ln, yr) + 1.0e-10);
            }
          }
        }
      }
      jnll_comp(15, rsp) -= offset_diet_w(rsp);
    }


    Type Denom = 0;

    // Calculate the predicted fraction by length-class (Eqn 17)
    T_hat.setZero();
    for (rsp = 0; rsp < nspp; rsp++) {
      for (ksp = 0; ksp < nspp; ksp++) {

        vector<Type> eaten_lmy(srv_age_bins(ksp)); eaten_lmy.setZero(); // no. of prey@length eaten by a predator length during iyr

        for (r_ln = 0; r_ln < srv_age_bins(rsp); r_ln++) { //
          // This is Equation 17
          for (yr = 0; yr < nyrs; yr++) { // FIXME: loop by stomach year
            // int stom_yr = yrs_stomlns(rsp, ksp, stm_yr); // FIXME: index by stomach year
            for (k_ln = 0; k_ln < srv_age_bins(ksp); k_ln ++) {
              for (k_age = 0; k_age < nages(ksp); k_age++) {
                eaten_lmy(k_ln) += eaten_la(rsp, ksp, r_ln, k_age, yr) * age_trans_matrix(k_age, k_ln, ksp);
              }
              T_hat(rsp, ksp, r_ln, k_ln) += stom_tau * eaten_lmy(k_ln); // FIXME: add actual stomach content sample size
            }
          }


          Denom = 0;
          for (k_ln = 0; k_ln < srv_age_bins(ksp); k_ln++) {
            Denom += T_hat(rsp, ksp, r_ln, k_ln);
          }

          // Renormalize the eaten vector
          for (k_ln = 0; k_ln < srv_age_bins(ksp); k_ln++) {
            T_hat(rsp, ksp, r_ln, k_ln) /= Denom;

            // Likelihood of diet length         / This is equation 16
            if (Uobs(rsp, ksp, r_ln, k_ln) > 0) {
              jnll_comp(16, rsp) -= stom_tau * Uobs(rsp, ksp, r_ln, k_ln) * log(T_hat(rsp, ksp, r_ln, k_ln)  + 1.0e-10);
            }
          }
        }
      }
      jnll_comp(16, rsp) -= offset_diet_l(rsp);
    }
  } // End if statement for diet likelihood


  // ------------------------------------------------------------------------- //
  // 12. REPORT SECTION                                                        //
  // ------------------------------------------------------------------------- //

  // 12.1. Population components
  REPORT( mn_rec );
  REPORT( Zed );
  REPORT( NByage );
  REPORT( AvgN );
  REPORT( S );
  REPORT( biomassByage );
  REPORT( biomassSSBByage );
  REPORT( biomass );
  REPORT( biomassSSB );
  REPORT( pmature );
  REPORT(r_sigma);
  REPORT( R );
  REPORT( M1 );
  REPORT( M );

  ADREPORT( Zed );
  ADREPORT( biomass );
  ADREPORT( biomassSSB );
  ADREPORT( r_sigma );
  ADREPORT( R );
  ADREPORT( M1);

  // -- 12.2. Survey components
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


  // -- 12.3. Fishery components
  REPORT( F );
  REPORT( F_dev );
  REPORT( fsh_sel );
  REPORT( avgsel_fsh );
  REPORT( tc_biom_hat );
  REPORT( catch_hat );
  REPORT( tc_hat );
  REPORT( fsh_age_hat );
  REPORT( fsh_age_obs );
  ADREPORT( F );


  // -- 12.4. Likelihood components
  REPORT( jnll_comp );
  REPORT( offset_srv );
  REPORT( offset_fsh );
  REPORT( offset_eit );
  REPORT( offset_diet_w );
  REPORT( offset_diet_l );
  REPORT( offset_uobsagewt );


  // -- 12.5. Ration components
  REPORT( ConsumAge );
  REPORT( Consum_livingAge );
  REPORT( S2Age );
  REPORT( LbyAge );
  REPORT( mnWt_obs );
  REPORT( fT );
  REPORT( TempC );
  REPORT( ration2Age );


  // 12.6. Suitability components
  REPORT( suma_suit );
  REPORT( suit_main );
  REPORT( suit_other );
  REPORT( stom_div_bio2 );
  REPORT( stomKir );
  REPORT( avail_food );
  REPORT( of_stomKir );
  REPORT( M2 );
  REPORT( B_eaten );
  REPORT( UobsWtAge_hat );
  REPORT( mn_UobsWtAge_hat );

  // -- 12.7. Kinzey predation functions
  REPORT( H_1 );
  REPORT( H_1a );
  REPORT( H_1b );
  REPORT( H_2 );
  REPORT( H_3 );
  REPORT( gam_a );
  REPORT( gam_b );

  REPORT( N_pred_yrs );
  REPORT( N_prey_yrs );
  REPORT( N_pred_eq );
  REPORT( N_prey_eq );

  REPORT( pred_resp );
  REPORT( Pred_r );
  REPORT( Prey_r );
  REPORT( Vmort_ua );
  REPORT( eaten_la );
  REPORT( eaten_ua );

  REPORT( Q_mass_l );
  REPORT( Q_mass_u );
  REPORT( Q_other_u );
  REPORT( Q_hat );
  REPORT( T_hat );

  REPORT( omega_hat );
  REPORT( omega_hat_ave );

  // ------------------------------------------------------------------------- //
  // 13. END MODEL                                                             //
  // ------------------------------------------------------------------------- //
  Type jnll = 0;

  // Estimation mode
  if (debug == 0) {
    jnll = jnll_comp.sum();
  }

  // Debug mode
  if (debug > 0) {
    jnll = dummy * dummy;
  }

  REPORT( jnll );
  return jnll;
}