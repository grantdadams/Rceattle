##' @title Specify the harvest control rule (HCR) used for Rceattle
##'
##' @description Defines the harvest control rule and its associated reference points.
##'
##' @param HCR Harvest control rule to use. Accepts an integer or equivalent string alias. Default = 0. See Details for the full list of available options.
##' @param DynamicHCR TRUE/FALSE. Whether to use static or dynamic reference points (default = FALSE).
##' @param Ftarget Target fishing mortality rate (yr^-1) (SPR or depletion based) or input F for projections. For example, if Ftarget is spr F40%, enter 0.40.
##' @param Flimit Limit fishing mortality rate (yr^-1) (SPR or depletion based). For example, if Flimit is spr F35%, enter 0.35.
##' @param Ptarget Target spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0 (accounts for recruitment).
##' @param Plimit Limit spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0 (accounts for recruitment).
##' @param Alpha Parameter used in NPFMC Tier 3 HCR.
##' @param Pstar Quantile used for uncertainty buffer given \code{Flimit} and \code{Sigma}.
##' @param Sigma Standard deviation used for normally distributed uncertainty buffer given \code{Flimit} and \code{Pstar}.
##' @param Fmult Multiplier on the target F (default = 1). Used to scale the target fishing mortality following NEFSC convention.
##' @param HCRorder For multi-species models, the order in which to project fishing (e.g., predators first, then prey).
##'
##' @details
##' **Harvest control rule formulations currently implemented in Rceattle:**
##'
##' \code{hcr = 0} or \code{"NoFishing"}: No catch. Estimate the hindcast.
##'
##' \code{hcr = 1} or \code{"CMSY"}: CMSY. Maximize catch across all species simultaneously. CMSY can be constrained such that depletion does not fall below \code{Plimit}.
##'
##' \code{hcr = 2} or \code{"ConstantF"}: Constant input F set at \code{Ftarget} for each species (vector or single F). SPR (single-species only) based Flimit is specified via \code{Flimit}.
##'
##' \code{hcr = 3} or \code{"ConstantFSSB"}: F that achieves \code{Ftarget}% of SSB0 in the end of the projection.
##'
##' \code{hcr = 4} or \code{"ConstantFSPR"}: Constant Fspr set at \code{Ftarget} for each species. Can be multiplied by \code{Fmult} following NEFSC.
##'
##' \code{hcr = 5} or \code{"NPFMC"}: The North Pacific Fishery Management Council (NPFMC) Tier 3 spawner-per-recruit-based harvest control rule:
##' 	Stock status: \eqn{SB > SB at Ftarget}
##' 	\deqn{Fofl = Flimit}
##' 	\deqn{Fuse = Ftarget}
##' 	Stock status: \eqn{Alpha < SB / SB at Ftarget \le 1}
##' 	\deqn{Fofl = Flimit * (SB / SB_{Ftarget} - Alpha) / (1 - Alpha)}
##' 	\deqn{Fuse = Ftarget * (SB / SB_{Ftarget} - Alpha) / (1 - Alpha)}
##' 	Stock status: \eqn{SB / SB at Ftarget \le Alpha} or \eqn{SB < Plimit * SB0}
##' 	\deqn{Fofl = 0}
##' 	\deqn{Fuse = 0}
##'
##' \code{hcr = 6} or \code{"PFMC"}: An HCR based on the Pacific Fishery Management Council (PFMC) category 1 40-10 annual catch limit (ABC) harvest control rule assuming Fofl is normally distributed with a standard deviation (sigma) = 0.5 and an uncertainty quantile buffer (P*) of 0.45 (PFMC 2020). The model uses Fspr if single-species or F that achieves X% of SSB0 for multi-species. The uncertainty buffer is the normal quantile function \code{qnorm(Pstar, mean, Sigma)}; note \code{qnorm(Pstar, Flimit, Sigma) = Flimit + qnorm(Pstar, 0, Sigma)} -- an F below \code{Flimit} when \eqn{Pstar < 0.5}. The taper runs between \eqn{SB0 * Plimit} and \eqn{SB0 * Ptarget}, so the 40-10 shape requires \code{Ptarget = 0.40} and \code{Plimit = 0.10}; the default \code{Plimit = 0} gives a 40-0 rule. Target biological reference points are:
##' 	Stock status: \eqn{SB > SB0 * Ptarget}
##' 	\deqn{Fofl = Flimit}
##' 	\deqn{Fuse = qnorm(Pstar, Flimit, Sigma)}
##' 	Stock status: \eqn{SB0 * Plimit < SB \le SB0 * Ptarget}
##' 	\deqn{Fofl = Flimit}
##' 	\deqn{Fuse = qnorm(Pstar, Flimit, Sigma) * \frac{SB0 * Ptarget * (SB - SB0 * Plimit)}{SB * SB0 * (Ptarget - Plimit)}}
##' 	Stock status: \eqn{SB < SB0 * Plimit}
##' 	\deqn{Fofl = 0}
##' 	\deqn{Fuse = 0}
##'
##' \code{hcr = 7} or \code{"SESSF"}: An HCR based on the The Southern and Eastern Scalefish and Shark Fishery (SESSF) spawner-per-recruit-based Tier 1 harvest control rule where F_Limit=F_(20%), B_Limit=SB_20, F_Target (AFMA 2017) calculated as follows:
##' 	Stock status: \eqn{SB > SB0 * Ptarget}
##' 	\deqn{Fofl = Flimit}
##' 	\deqn{Fuse = Ftarget}
##' 	Stock status: \eqn{SB0 * Ptarget > SB > SB0 * Plimit}
##' 	\deqn{Fofl = Flimit * (SB / (SB0 * Plimit) - 1)}
##' 	\deqn{Fuse = Ftarget * (SB / (SB0 * Plimit) - 1)}
##' 	Stock status: \eqn{SB < SB0 * Plimit}
##' 	\deqn{Fofl = 0}
##' 	\deqn{Fuse = 0}
##'
##' NOTE: only HCRs 0, 1, 2, 3, and 6 will work in multi-species mode.
##'
##' @return A \code{list} containing the harvest control rule and associated biological reference points.
##' @export
##'
build_hcr <- function(HCR = 0, DynamicHCR = FALSE, Ftarget = 0.40, Flimit = 0.35, Ptarget = 0.4, Plimit = 0.0, Alpha = 0.05, Pstar = 0.45, Sigma = 0.5, Fmult = 1, HCRorder = 1) {
  if(0 %in% Alpha & HCR %in% c(5, "NPFMC")){stop(paste0("Alpha = 0 for NPFMC Tier 3 HCR, divide by zero error"))}
  list(HCR = HCR, DynamicHCR = DynamicHCR, Ftarget = Ftarget, Flimit = Flimit, Ptarget = Ptarget, Plimit = Plimit, Alpha = Alpha, Pstar = Pstar, Sigma = Sigma, Fmult = Fmult, HCRorder = HCRorder)
}


#' Function to construct the TMB map argument for CEATTLE for projecting under alternative harvest control rules
#'
#' @description Reads a data list and map to update the map argument based on the HCR specified in \code{\link{build_hcr}}
#'
#' @param data_list an Rceattle data_list
#' @param map a map object created from \code{\link{build_map}}.
#' @param debug logical. If TRUE, turns off all parameters for debugging (default = FALSE).
#' @param all_params_on logical. If TRUE, leaves all hindcast parameters turned on (default = FALSE).
#' @param HCRiter for multi-species models, the order in which to project fishing (e.g. predators first, then prey)
#'
#' @return a list of map arguments for each parameter
#' @export
build_hcr_map <- function(data_list, map, debug = FALSE, all_params_on = FALSE, HCRiter = 1){

  # Turn off all population/fleet parameters ---- and turn on Fspr parameters
  params_on <- 1:data_list$nspp
  if(!all_params_on){
    map$mapList = sapply(map$mapList, function(x) replace(x, values = rep(NA, length(x))))
    yrs_proj = data_list$endyr:data_list$projyr - data_list$styr
    params_on <- c(1:data_list$nspp)[which(data_list$HCRorder <= HCRiter)]
  }

  # Turn on Fspr parameters depending on HCR ----
  # -- HCR = 0: No catch - Params off
  # -- HCR = 1: Constant catch - Params off
  # -- HCR = 2: Constant input F - Params off
  # -- HCR = 3: F that achieves X% of SSB0 in the end of the projection - Ftarget on
  # -- HCR = 4: Constant target Fspr - Ftarget on
  # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
  # -- HCR = 6: PFMC Cat 1 - Flimit on
  # -- HCR = 7: SESSF Tier 1 - Flimit and Ftarget on
  # --- Dynamic BRPS - 1 value per species and year
  if(!debug){
    if(data_list$HCR == "CMSY"){ # CMSY
      map$mapList$log_Ftarget[params_on] <- params_on
    }

    if(data_list$HCR == "ConstantF"){ # Fixed F - still have Flimit for single-species
      map$mapList$log_Flimit[params_on] <- params_on
    }
    if(data_list$HCR == "ConstantFSSB"){
      map$mapList$log_Ftarget[params_on] <- params_on
    }
    if(data_list$HCR == "ConstantFSPR"){
      map$mapList$log_Ftarget[params_on] <- params_on
      map$mapList$log_Flimit[params_on] <- params_on
    }
    if(data_list$HCR %in% c("NPFMC", "SESSF")){
      map$mapList$log_Ftarget[params_on] <- params_on
      map$mapList$log_Flimit[params_on] <- params_on
    }
    if(data_list$HCR == "PFMC"){
      map$mapList$log_Flimit[params_on] <- params_on
    }


    # Turn off SPR parameters for special cases ----
    # -- Turn off SPR parameters for species with no fishing (sum(Proj_F_proportion) == 0)
    # -- Turn off SPR parameters for species with fixed Nbyage
    for(sp in 1:data_list$nspp){

      # Check proj F if proj F prop is all 0
      prop_check <- data_list$fleet_control$Proj_F_proportion[which(data_list$fleet_control$Species == sp & data_list$fleet_control$Fleet_type == "Fishery")]
      # Turn off future F only when *every* fishery for this species is inactive.
      if(sum(prop_check, na.rm = TRUE) == 0){
        message(paste("F_prop for species",sp,"sums to 0"))
        map$mapList$log_Ftarget[sp] <- NA
        map$mapList$log_Flimit[sp] <- NA
      }

      # Fixed n-at-age: Turn off parameters
      if(data_list$estDynamics[sp] > 0){
        map$mapList$log_Ftarget[sp] <- NA
        map$mapList$log_Flimit[sp] <- NA
      }
    }
  }
  if(debug){
    map$mapList$dummy = 1
  }

  # Convert to factor ----
  map$mapFactor <- sapply(map$mapList, factor)
  return(map)
}
