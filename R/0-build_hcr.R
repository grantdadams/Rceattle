##' Specify the harvest control rule (HCR) used for Rceattle
##'
##' @param HCR Function to be used for harvest control rule (see below). Default = 0
##' @param DynamicHCR TRUE/FALSE. Wether to use static or dynamic reference points (default = FALSE)
##' @param FsprTarget Target Fspr or input F for projections. Percentage of spawning stock biomass-per-recruit in the absence of fishing to define the target fishing mortality rate (yr^-1). For example, if FsprLimit is F40%, enter 0.40.
##' @param FsprLimit Limit Fspr. Percentage of spawning stock biomass-per-recruit in the absence of fishing to define the limit fishing mortality rate (yr^-1). For example, if FsprLimit is F35%, enter 0.35.
##' @param Ctarget Target catch (mt)
##' @param Ptarget Target spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0 (accounts for recruitment)
##' @param Plimit Limit spawning-stock biomass as a percentage of static or dynamic spawning-stock-biomass at F = 0 (accounts for recruitment)
##' @param Alpha Parameter used in NPFMC Tier 3 HCR
##' @param Pstar Quantile used for uncertainty buffer given \code{FsprLimit} and \code{Sigma}
##' @param Sigma Standard deviation used for normally distributed uncertainty buffer given \code{FsprLimit} and \code{Pstar}
##' @details
##'
##' **Harvest control rule formulations currently implemented in Rceattle:**
##'
##' \code{hcr = 0} No catch. Estimate the hindcast.
##'
##' \code{hcr = 1} Constant catch set at \code{Ctarget}
##'
##' \code{hcr = 2} Constant input F set at \code{FsprTarget} for each species (vector or single F)
##'
##'
##' \code{hcr = 3} F that achieves \code{FsprTarget}% of SSB0 in the end of the projection
##'
##' \code{hcr = 4} Constant Fspr set at \code{FsprTarget} for each species. Can be multiplied by \code{Fmult} following NEFSC.
##'
##' \code{hcr = 5} The North Pacific Fishery Management Council (NPFMC) Tier 3 spawner-per-recruit-based harvest control rule:
##' 	Stock status: \eqn{ SB > SB at FsprTarget)}
##' 	\deqn{Fofl = FsprLimit }
##' 	\deqn{Fuse = FsprTarget}
##' 	Stock status: \eqn{Alpha < SB / SB at FsprTarget at  ≤1}
##' 	\deqn{Fofl = FsprLimit * (SB / Ptarget - Alpha)/(1-Alpha)}
##' 	\deqn{Fuse = FsprTarget * (SB/ Ptarget - Alpha)/(1-Alpha)}
##' 	Stock status: \eqn{SB/SB at FsprTarget ≤ Alpha} or \eqn{SB < Plimit*SB0}
##' 	\deqn{Fofl=0}
##' 	\deqn{Fuse=0}
##'
##' \code{hcr = 6} An HCR based on the The Pacific Fishery Management Council (PFMC) category 1 40-10 annual catch limit (ABC) harvest control rule assuming Fofl is normally distributed with a standard deviation (sigma) = 0.5 and an uncertainty quantile buffer (P∗) of 0.45 (PFMC 2020). Target biological reference points are calculated based on the normal cumulative distribution function (Φ(quantile,mean,standard deviation)), P∗, and sigma as follows:
##' 	Stock status: \eqn{ SB > SB0*Ptarget)}
##' 	\deqn{Fofl = FsprLimit}
##' 	\deqn{Fuse = Φ(Pstar, FsprLimit, Sigma)}
##' 	Stock status: \eqn{SB0*Plimit < SB ≤ SB0*Ptarget}
##' 	\deqn{Fofl = FsprLimit}
##' 	\deqn{Fuse = Φ(Pstar, FsprLimit, Sigma)*{SB0*Ptarget(SB - SB0*Plimit)}{SB(SB0*Ptarget - SB0*Plimit)}}
##' 	Stock status: \eqn{SB < SB0*Plimit}
##' 	\deqn{Fofl=0}
##' 	\deqn{Fuse=0}
##'
##' \code{hcr = 7} An HCR based on the The Southern and Eastern Scalefish and Shark Fishery (SESSF) spawner-per-recruit-based Tier 1 harvest control rule where F_Limit=F_(20%), B_Limit=〖SB〗_20, F_Target (AFMA 2017) calculated as follows:
##' 	Stock status: \eqn{ SB > SB0*Ptarget)}
##' 	\deqn{Fofl = FsprLimit}
##' 	\deqn{Fuse = FsprTarget}
##' 	Stock status: \eqn{Ptarget > SB > SB0*Plimit}
##' 	\deqn{Fofl = FsprLimit * (SB / SB0*Ptarget - 1)}
##' 	\deqn{Fuse = FsprTarget * (SB/ SB0*Ptarget - 1)}
##' 	Stock status: \eqn{SB < SB0*Plimit}
##' 	\deqn{Fofl=0}
##' 	\deqn{Fuse=0}
##'
##' @return A \code{list} containing the harvest control rule and associated biological reference points
##' @export
##'
build_hcr <- function(HCR = 0, DynamicHCR = FALSE, FsprTarget = 0.40, FsprLimit = 0.35, Ptarget = 0.4, Plimit = 0.2, Alpha = 0.05, Pstar = 0.45, Sigma = 0.5, Fmult = 1) {
  if(0 %in% Alpha & HCR == 5){stop(paste0("Alpha = 0 for Tier 3 HCR, divide by zero error"))}
  list(HCR = HCR, DynamicHCR = DynamicHCR, FsprTarget = FsprTarget, FsprLimit = FsprLimit, Ptarget = Ptarget, Plimit = Plimit, Alpha = Alpha, Pstar = Pstar, Sigma = Sigma, Fmult = Fmult)
}
