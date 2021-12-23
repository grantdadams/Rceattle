#' Function to construct the TMB map argument for CEATTLE for projecting under alternative harvest control rules
#'
#' @description Reads a data list and map to update the map argument based on the HCR specified in \code{\link{build_hcr}}
#'
#' @param data_list a data_list created from \code{\link{build_dat}}.
#' @param map a map object created from \code{\link{build_map}}.
#'
#' @return a list of map arguments for each parameter
#' @export
build_hcr_map <- function(data_list, map){

  # Step 1 - Turn off all population/fleet parameters and turn on Fspr parameters
  map$mapList = sapply(map$mapList, function(x) replace(x, values = rep(NA, length(x))))

  # Step 2 - Turn on Fspr parameters depending on HCR
  # -- HCR = 0: No catch - Params off
  # -- HCR = 1: Constant F - Params off
  # -- HCR = 2: Constant catch - Params off
  # -- HCR = 3: NPFMC Tier 3 - Flimit and Ftarget on
  # -- HCR = 4: PFMC Cat 1 - Flimit on
  # -- HCR = 5: SESSF Tier 1 - Flimit and Ftarget on
  # --- Dynamic BRPS - 1 value per species and year
  if(data_list$DynamicHCR){
    if(data_list$HCR %in% c(3,5)){
      map$mapList$ln_Ftarget <- replace(map$mapList$ln_Ftarget, values = c(1:length(map$mapList$ln_Ftarget)))
      map$mapList$ln_Flimit <- replace(map$mapList$ln_Flimit, values = c(1:length(map$mapList$ln_Flimit)))
    }
    if(data_list$HCR == 4){
      map$mapList$ln_Flimit <- replace(map$mapList$ln_Flimit, values = c(1:length(map$mapList$ln_Flimit)))
    }
  }
  # --- Static BRPS - 1 value per species
  if(data_list$DynamicHCR == FALSE){
    if(data_list$HCR %in% c(3,5)){
      map$mapList$ln_Ftarget <- replace(map$mapList$ln_Ftarget, values = rep(1:data_list$nspp, ncol(length(map$mapList$ln_Ftarget))))
      map$mapList$ln_Flimit <- replace(map$mapList$ln_Flimit, values = rep(1:data_list$nspp, ncol(length(map$mapList$ln_Flimit))))
    }
    if(data_list$HCR == 4){
      map$mapList$ln_Flimit <- replace(map$mapList$ln_Flimit, values = rep(1:data_list$nspp, ncol(length(map$mapList$ln_Flimit))))
    }
  }

  # Step 3 - Turn off SPR parameters for special cases
  # -- Turn off SPR parameters for species with no fishing (sum(proj_F_prop) == 0)
  # -- Turn off SPR parameters for species with fixed Nbyage
  for(sp in 1:data_list$nspp){

    # Check proj F if proj F prop is all 0
    prop_check <- data_list$fleet_control$proj_F_prop[which(data_list$fleet_control$Species == sp & data_list$fleet_control$Fleet_type == 1)]
    if(sum(as.numeric(prop_check == 0)) != 0){ # If all fisheries for a species have no F in F_prop, turn off future F
      map$mapList$ln_Ftarget[sp,] <- NA
      map$mapList$ln_Flimit[sp,] <- NA
    }

    # Fixed n-at-age: Turn off parameters
    if(data_list$estDynamics[sp] > 0){
      map$mapList$ln_Ftarget[sp,] <- NA
      map$mapList$ln_Flimit[sp,] <- NA
    }
  }

  # STEP 5 -- Convert to factor
  map$mapFactor <- sapply(map$mapList, factor)
  return(map)
}
