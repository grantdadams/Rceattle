#' Function to construct the TMB map argument for CEATTLE for projecting under alternative harvest control rules
#'
#' @description Reads a data list and map to update the map argument based on the HCR specified in \code{\link{build_hcr}}
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param map a map object created from \code{\link{build_map}}.
#' @param TRUE/FALSE map out all parameters
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
  # -- HCR = 3: F that acheives X% of SSB0 in the end of the projection - Ftarget on
  # -- HCR = 4: Constant target Fspr - Ftarget on
  # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
  # -- HCR = 6: PFMC Cat 1 - Flimit on
  # -- HCR = 7: SESSF Tier 1 - Flimit and Ftarget on
  # --- Dynamic BRPS - 1 value per species and year
  if(!debug){
    if(data_list$HCR %in% c(1)){ # CMSY
      map$mapList$ln_Ftarget[params_on] <- params_on
    }

    if(data_list$HCR %in% c(2)){ # Fixed F - still have Flimit for single-species
      map$mapList$ln_Flimit[params_on] <- params_on
    }
    if(data_list$HCR %in% c(3)){
      map$mapList$ln_Ftarget[params_on] <- params_on
    }
    if(data_list$HCR %in% c(4)){
      map$mapList$ln_Ftarget[params_on] <- params_on
      map$mapList$ln_Flimit[params_on] <- params_on
    }
    if(data_list$HCR %in% c(5,7)){
      map$mapList$ln_Ftarget[params_on] <- params_on
      map$mapList$ln_Flimit[params_on] <- params_on
    }
    if(data_list$HCR == 6){
      map$mapList$ln_Flimit[params_on] <- params_on
    }


    # Turn off SPR parameters for special cases ----
    # -- Turn off SPR parameters for species with no fishing (sum(proj_F_prop) == 0)
    # -- Turn off SPR parameters for species with fixed Nbyage
    for(sp in 1:data_list$nspp){

      # Check proj F if proj F prop is all 0
      prop_check <- data_list$fleet_control$proj_F_prop[which(data_list$fleet_control$Species == sp & data_list$fleet_control$Fleet_type == 1)]
      if(sum(as.numeric(prop_check == 0)) != 0){ # If all fisheries for a species have no F in F_prop, turn off future F
        print(paste("F_prop for species",sp,"sums to 0"))
        map$mapList$ln_Ftarget[sp] <- NA
        map$mapList$ln_Flimit[sp] <- NA
      }

      # Fixed n-at-age: Turn off parameters
      if(data_list$estDynamics[sp] > 0){
        map$mapList$ln_Ftarget[sp] <- NA
        map$mapList$ln_Flimit[sp] <- NA
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

