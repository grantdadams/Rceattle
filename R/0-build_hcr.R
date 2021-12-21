##' Specify the harvest control rule (HCR) used for Rceattle
##'
##' @param HCR Function to be used for harvest control rule (see below). Default = 0
##' @param DynamicHCR TRUE/FALSE. Wether to use static or dynamic reference points (default = FALSE)
##' @param FXSPRtarget Percentage of spawning stock biomass-per-recruit in the absence of fishing to define the target fishing mortality rate (yr^-1). For example, if FXSPRlimit is F40%, enter 0.40.
##' @param FXSPRlimit Percentage of spawning stock biomass-per-recruit in the absence of fishing to define the limit fishing mortality rate (yr^-1). For example, if FXSPRlimit is F35%, enter 0.35.
##' @param Ctarget Target catch (mt)
##' @param Ptarget Target spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0 (accounts for recruitment)
##' @param Plimit Limit spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0 (accounts for recruitment)
##' @param Alpha Parameter used in NPFMC Tier 3 HCR
##' @param Pstar Quantile used for uncertainty buffer given \code{FXSPRlimit} and \code{Sigma}
##' @param Sigma Standard deviation used for normally distributed uncertainty buffer given \code{FXSPRlimit} and \code{Pstar}
##' @details
##'
##' **Harvest control rule formulations currently implemented in Rceattle:**
##'
##' \code{hcr = 0} No catch. Estimate the hindcast.
##'
##' \code{hcr = 1} Constant F set at \code{FXSPRtarget} for each species (vector or single F)
##'
##' \code{hcr = 2} Constant catch set at \code{Ctarget}
##'
##' \code{hcr = 3} The North Pacific Fishery Management Council (NPFMC) Tier 3 spawner-per-recruit-based harvest control rule:
##' 	Stock status: \eqn{ SB > SB at FXSPRtarget)}
##' 	\deqn{Fofl = FXSPRlimit }
##' 	\deqn{Fuse = FXSPRtarget}
##' 	Stock status: \eqn{Alpha < SB / SB at FXSPRtarget at  ≤1}
##' 	\deqn{Fofl = FXSPRlimit * (SB / Ptarget - Alpha)/(1-Alpha)}
##' 	\deqn{Fuse = FXSPRtarget * (SB/ Ptarget - Alpha)/(1-Alpha)}
##' 	Stock status: \eqn{SB/SB at FXSPRtarget ≤ Alpha} or \eqn{SB < Plimit*SB0}
##' 	\deqn{Fofl=0}
##' 	\deqn{Fuse=0}
##'
##' \code{hcr = 4} An HCR based on the The Pacific Fishery Management Council (PFMC) category 1 40-10 annual catch limit (ABC) harvest control rule assuming Fofl is normally distributed with a standard deviation (sigma) = 0.5 and an uncertainty quantile buffer (P∗) of 0.45 (PFMC 2020). Target biological reference points are calculated based on the normal cumulative distribution function (Φ(quantile,mean,standard deviation)), P∗, and sigma as follows:
##' 	Stock status: \eqn{ SB > SB0*Ptarget)}
##' 	\deqn{Fofl = FXSPRlimit}
##' 	\deqn{Fuse = Φ(Pstar, FXSPRlimit, Sigma)}
##' 	Stock status: \eqn{SB0*Plimit < SB ≤ SB0*Ptarget}
##' 	\deqn{Fofl = FXSPRlimit}
##' 	\deqn{Fuse = Φ(Pstar, FXSPRlimit, Sigma)*{SB0*Ptarget(SB - SB0*Plimit)}{SB(SB0*Ptarget - SB0*Plimit)}}
##' 	Stock status: \eqn{SB < SB0*Plimit}
##' 	\deqn{Fofl=0}
##' 	\deqn{Fuse=0}
##'
##' \code{hcr = 5} An HCR based on the The Southern and Eastern Scalefish and Shark Fishery (SESSF) spawner-per-recruit-based Tier 1 harvest control rule where F_Limit=F_(20%), B_Limit=〖SB〗_20, F_Target (AFMA 2017) calculated as follows:
##' 	Stock status: \eqn{ SB > SB0*Ptarget)}
##' 	\deqn{Fofl = FXSPRlimit}
##' 	\deqn{Fuse = FXSPRtarget}
##' 	Stock status: \eqn{Ptarget > SB > SB0*Plimit}
##' 	\deqn{Fofl = FXSPRlimit * (SB / SB0*Ptarget - 1)}
##' 	\deqn{Fuse = FXSPRtarget * (SB/ SB0*Ptarget - 1)}
##' 	Stock status: \eqn{SB < SB0*Plimit}
##' 	\deqn{Fofl=0}
##' 	\deqn{Fuse=0}
##'
##' @return A \code{list} containing the harvest control rule and associated biological reference points
##' @export
##'
build_hcr <- function(HCR = 0, DynamicHCR = FALSE, FXSPRtarget = 0.40, FXSPRlimit = 0.35, Ptarget = 0.4, Plimit = 0.2, Alpha = 0.05, Pstar = NULL, Sigma = NULL) {
  if(Alpha == 0 & HCR == 3){stop(paste0("Alpha = 0 for Tier 3 HCR, divide by zero error"))}
  list(HCR = HCR, DynamicHCR = DynamicHCR, FXSPRtarget = FXSPRtarget, FXSPRlimit = FXSPRlimit, Ptarget = Ptarget, Plimit = Plimit, Alpha = Alpha, Pstar = Pstar, Sigma = Sigma)
}
