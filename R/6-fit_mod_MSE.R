#' This functions runs CEATTLE for MSEs (estimation is trimmed down)
#' @description This function estimates population parameters of CEATTLE using maximum likelihood in TMB.
#'
#' @param data_list a data_list created from BSAI CEATTLE dat files \code{\link{build_dat}}, prebuilt data_list \code{\link{BS2017SS}}, or read in using \code{\link{read_excel}}.
#' @param inits (Optional) Character vector of named initial values from previous parameter estimates from Rceattle model. If NULL, will use 0 for starting parameters. Can also construct using \code{\link{build_params}}
#' @param map (Optional) A map object from \code{\link{build_map}}.
#' @param estimateMode 0 = Fit the hindcast model and projection with HCR specified via \code{HCR}. 1 = Fit the hindcast model only (no projection). 2 = Run the projection only with HCR specified via \code{HCR} given the initial parameters in \code{inits}.  3 = debug mode 1: runs the model through MakeADFun, but not nlminb, 4 = runs the model through MakeADFun and nlminb (will all parameters mapped out).
#' @param niter Number of iterations for multispecies model
#' @param msmMode The predation mortality functions to used. Defaults to no predation mortality used.
#' @param initMode how the population is initialized. 0 = initial age-structure estimated as free parameters; 1 = equilibrium age-structure estimated out from R0,  mortality (M1), and initial population deviates; 2 = non-equilibrium age-structure estimated out from initial fishing mortality (Finit), R0,  mortality (M1), and initial population deviates.
#' @param suitMode Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = length based gamma suitability, 2 = weight based gamma suitability, 3 = length based lognormal selectivity, 4 = time-varying length based lognormal selectivity.
#' @param suit_meanyr Integer. The last year used to calculate mean suitability, starting at \code{styr}. Defaults to $endyr$ in $data_list$. Used for MSE runs where suitability is held at the value estimated from the years used to condition the OM, but F is estimated for years beyond those used to condition the OM to account for projected catch.
#' @param use_gradient use the gradient to phase (default = TRUE).
#' @param rel_tol The relative tolerance for discontinuous likelihood warnings. Set to 1. This evaluates the difference between the TMB object likelihood and the nlminb likelihood.
#' @param control A list of control parameters. For details see \code{?nlminb}
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param verbose 0 = Silent, 1 = print updates of model fit, 2 = print updates of model fit and TMB estimation progress.
#'
#' @return A list of class "Rceattle" including:
#'
#' \itemize{
#'  \item{data_list: List of data inputs}
#'  \item{initial_params: List of starting parameters}
#'  \item{map: List of map used in TMB}
#'  \item{obj: TMB model object}
#'  \item{opt: Optimized model object from `nlimb`}
#'  \item{estimated_params: List of estimated parameters}
#'  \item{quantities: Derived quantities from CEATTLE}
#'  \item{run_time: Model run time}
#'  }
#' @export
fit_mod_mse <-
  function(
    data_list = NULL,
    inits = NULL,
    map = NULL,
    estimateMode = 0,
    niter = 3,
    msmMode = 0,
    initMode = 1,
    suitMode = 0,
    suit_meanyr = NULL,
    use_gradient = TRUE,
    rel_tol = 1,
    control = list(eval.max = 1e+09,
                   iter.max = 1e+09, trace = 0),
    loopnum = 5,
    verbose = 1,
    newtonsteps = 0){

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Debugging section ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
#     use_gradient = TRUE;
#     rel_tol = 1;
#     control = list(eval.max = 1e+09,
#                    iter.max = 1e+09, trace = 0);
#     loopnum = 1;
#     verbose = 1;
#     newtonsteps = 0
#     data_list = em_use$data_list;
#     inits = em_use$estimated_params;
#     map =  NULL;
#     estimateMode = 0; # Run hindcast and projection; otherwise debug
#     niter = em_use$data_list$niter;
#     msmMode = em_use$data_list$msmMode;
#     suitMode = em_use$data_list$suitMode;
#     suit_meanyr = em_use$data_list$suit_meanyr;
#     initMode = em_use$data_list$initMode;
#     loopnum = loopnum;
#     verbose = 2


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 1 - Load data ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    start_time <- Sys.time()
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 2: Load/build parameters ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    start_par <- inits

    # Set Fdev for years with 0 catch to very low number
    fsh_biom_sub <- data_list$fsh_biom %>%
      filter(Year <= data_list$endyr)
    fsh_ind <- fsh_biom_sub$Fleet_code[which(fsh_biom_sub$Catch == 0)]
    yr_ind <- fsh_biom_sub$Year[which(fsh_biom_sub$Catch == 0)] - data_list$styr + 1
    for(i in 1:length(fsh_ind)){
      start_par$F_dev[fsh_ind[i], yr_ind[i]] <- -999
    }
    rm(fsh_biom_sub)

    # Start Flimit at 0 so forecast is 0
    start_par$ln_Flimit[] <- -10
    start_par$ln_Ftarget[] <- -10



    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 3: Load/build map ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(map)) {
      map <- suppressWarnings(build_map(data_list, start_par, debug = estimateMode == 4, random_rec = FALSE, random_sel = FALSE, MSE = TRUE))
    } else{
      map <- map
    }
    if(verbose > 0) {message("Step 1: Map build complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 4: Reorganize data ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    TMBfilename <- "ceattle_v01_10"

    Rceattle:::data_check(data_list)

    data_list_reorganized <- Rceattle::rearrange_dat(data_list)
    data_list_reorganized = c(list(model = TMBfilename), data_list_reorganized)

    # Include BRPs in likelihood of hindcast
    # -- but set Ftarget and Flimit to 0
    data_list_reorganized$forecast <- TRUE
    data_list_reorganized$MSSB0[params_on] <- 0
    data_list_reorganized$MSB0[params_on] <- 0
    if(verbose > 0) {message("Step 2: Data rearranged complete")}


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 5: Fit hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if(estimateMode != 2){ # dont build if projection and estimating HCR parameters
      if(sum(as.numeric(unlist(map$mapFactor)), na.rm = TRUE) == 0){stop("Map of length 0: all NAs")}
      obj = TMB::MakeADFun(
        data_list_reorganized,
        parameters = start_par,
        DLL = TMBfilename,
        map = map$mapFactor,
        silent = verbose != 2
      )
    }

    # -- Save objects
    mod_objects <-
      list(
        TMBfilename = TMBfilename,
        initial_params = start_par,
        map = map
      )

    if(verbose > 0) {message(paste0("Step 3: final build complete. Optimizing."))}


    # -- Optimize hindcast
    fn1 <- function(Fparams){
      par <- obj$env$last.par.best
      notFind <- which(!names(par) %in% c("ln_Flimit", "ln_Ftarget"))
      par[notFind] <- Fparams
      obj$fn(par)
    }

    gr1 <- function(Fparams){
      par <- obj$env$last.par.best
      notFind <- which(!names(par) %in% c("ln_Flimit", "ln_Ftarget"))
      par[notFind] <- Fparams
      as.numeric(obj$gr(par)[notFind])
    }

    par <- obj$env$last.par.best
    notFind <- which(!names(par) %in% c("ln_Flimit", "ln_Ftarget"))


    if(estimateMode %in% c(0,1,4)){
      opt = Rceattle::fit_tmb_mse(
        fn=fn1,
        gr=gr1,
        startpar=par[notFind],
        loopnum = loopnum,
        control = control,
        quiet = verbose < 2
      )
      if(verbose > 0) {message("Step 5: Final optimization complete")
      }
    }

    # -- Get MLEs
    if (estimateMode > 1) { # Debugging and projection only: use initial parameters
      last_par <- start_par
    } else{
      last_par = try(obj$env$parList(obj$env$last.par.best)) # FIXME: maybe add obj$env$last.par.best inside?
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 6: Run HCR projections ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    fn2 <- function(Fparams){
      par <- obj$env$last.par.best
      Flimitind <- which(names(par) == "ln_Flimit")[params_on]
      Flimitind <- Flimitind[which(!is.na(Flimitind))]

      Ftargetind <- which(names(par) == "ln_Ftarget")[params_on]
      Ftargetind <- Ftargetind[which(!is.na(Ftargetind))]

      par[c(Flimitind, Ftargetind)] <- Fparams
      obj$fn(par)
    }

    gr2 <- function(Fparams){
      par <- obj$env$last.par.best
      Flimitind <- which(names(par) == "ln_Flimit")[params_on]
      Flimitind <- Flimitind[which(!is.na(Flimitind))]

      Ftargetind <- which(names(par) == "ln_Ftarget")[params_on]
      Ftargetind <- Ftargetind[which(!is.na(Ftargetind))]

      par[c(Flimitind, Ftargetind)] <- Fparams
      as.numeric(obj$gr(par)[c(Flimitind, Ftargetind)])
    }

    if(estimateMode %in% c(0,2,4)){
      if(!data_list$HCR %in% c(0, 2)){ # - All HCRs except no F and fixed F

        # * Single species mode ----
        if(msmMode == 0){

          # -- Parameters to optimize
          params_on <- 1:data_list$nspp
          par <- obj$env$last.par.best
          Flimitind <- which(names(par) == "ln_Flimit")[params_on]
          Flimitind <- Flimitind[which(!is.na(Flimitind))]

          Ftargetind <- which(names(par) == "ln_Ftarget")[params_on]
          Ftargetind <- Ftargetind[which(!is.na(Ftargetind))]

          # -- Start from previous MLEs
          Fstartpar <- suppressWarnings(
            c(inits$ln_Flimit[Flimitind - min(Flimitind)+1],
              inits$ln_Ftarget[Ftargetind - min(Ftargetind)+1]
            )
          )
          names(Fstartpar) <- names(par[c(Flimitind, Ftargetind)])

          # -- Optimize
          # opt <- optim(par = par[c(Flimitind, Ftargetind)], fn = fn2, gr = gr2, control = list(maxit = 1e9))
          opt = Rceattle::fit_tmb_mse(
            fn=fn2,
            gr=gr2,
            startpar=Fstartpar,
            loopnum = loopnum,
            control = control,
            quiet = verbose < 2
          )
        }

        # * Multi-species mode ----
        if(msmMode > 0){

          # -- Get quantities
          if(estimateMode == 2){ # Build obj if we havent done so already
            obj = TMB::MakeADFun(
              data_list_reorganized,
              parameters = last_par,
              DLL = TMBfilename,
              map = map$mapFactor,
              silent = verbose != 2
            )
          }

          # Loop across species orders
          for(HCRiter in 1:max(data_list$HCRorder)){

            # -- Species to optimize
            params_on <- c(1:data_list$nspp)[which(data_list$HCRorder == HCRiter)]

            # -- Get SB0: SSB when model is projected forward under no fishing
            quantities <- obj$report(obj$env$last.par)
            SB0 <- quantities$biomassSSB[, ncol(quantities$biomassSSB)]
            B0 <- quantities$biomass[, ncol(quantities$biomass)]
            data_list_reorganized$MSSB0[params_on] <- SB0[params_on]
            data_list_reorganized$MSB0[params_on] <- B0[params_on]
            obj$env$data$MSSB0 <- data_list_reorganized$MSSB0
            obj$env$data$MSB0 <- data_list_reorganized$MSB0

            # -- Parameters to optimize
            par <- obj$env$last.par
            Flimitind <- which(names(par) == "ln_Flimit")
            Flimitind <- Flimitind[which(!is.na(Flimitind))]

            Ftargetind <- which(names(par) == "ln_Ftarget")
            Ftargetind <- Ftargetind[which(!is.na(Ftargetind))]

            # -- Start from previous MLEs
            Fstartpar <- suppressWarnings(
              c(inits$ln_Flimit[Flimitind - min(Flimitind)+1][params_on],
                inits$ln_Ftarget[Ftargetind - min(Ftargetind)+1][params_on]
              )
            )
            Fstartpar <- Fstartpar[which(!is.na(Fstartpar))]
            Fparnames <- names(par[c(Flimitind[params_on],
                                            Ftargetind[params_on])])
            names(Fstartpar) <- Fparnames[which(!is.na(Fparnames))]

            # -- Optimize
            opt <- optim(par = Fstartpar, fn = fn2, control = list(maxit = 1e9))
            # opt = Rceattle::fit_tmb_mse(
            #   fn=fn2,
            #   gr=gr2,
            #   startpar=Fstartpar,
            #   loopnum = loopnum,
            #   control = control,
            #   quiet = verbose < 2
            # )
          }
        }

        if(verbose > 0) {message("Step 6: Projections complete")}

        # -- Update MLEs
        if (estimateMode > 2) { # Debugging, give initial parameters
          last_par <- start_par
        }else{
          last_par = try(obj$env$parList())
        }
      } # End estimable BRP/HCR projections
    } # End projection


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # STEP 7: Save output ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # - Save estimated parameters
    mod_objects$estimated_params <- last_par
    mod_objects$obj = obj

    # - Get quantities
    quantities <- obj$report(obj$env$last.par) # Last par rather than last par best

    # -- Rename jnll
    colnames(quantities$jnll_comp) <- paste0("Sp/Srv/Fsh_", 1:ncol(quantities$jnll_comp))
    rownames(quantities$jnll_comp) <- c(
      "Survey biomass",
      "Total catch",
      "Age/length composition data",
      "Sex ratio",
      "Non-parametric selectivity",
      "Selectivity deviates",
      "Selectivity normalization",
      "Catchability prior",
      "Catchability deviates",
      "Stock-recruit prior",
      "Recruitment deviates",
      "Initial abundance deviates",
      "Fishing mortality deviates",
      "SPR Calculation",
      "Zero n-at-age penalty",
      "M prior",
      "Ration",
      "Ration penalties",
      "Stomach content proportion by weight"
    )


    colnames(quantities$biomassSSB) <- data_list$styr:data_list$projyr
    colnames(quantities$R) <- data_list$styr:data_list$projyr

    rownames(quantities$biomassSSB) <- data_list$spnames
    rownames(quantities$R) <- data_list$spnames

    # -- Save derived quantities
    mod_objects$quantities <- quantities
    mod_objects$data_list <- data_list
    mod_objects$opt = opt
    mod_objects$run_time = ((Sys.time() - start_time))
    class(mod_objects) <- "Rceattle"

    # suppressWarnings(try(dyn.unload(TMB::dynlib(paste0(cpp_file)))))
    return(mod_objects)

    # Free up memory
    TMB::FreeADFun(obj)
  }

