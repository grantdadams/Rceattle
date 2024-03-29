% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11b-mse_summary.R
\name{mse_summary}
\alias{mse_summary}
\title{Management strategy evaluation performance metric summary}
\usage{
mse_summary(mse)
}
\arguments{
\item{mse}{MSE runs from \code{\link{mse_run}} or \code{\link{load_mse}}}
}
\value{
Alist of two data.frames with MSE summary statistics of performance metrics including:
data.frame 1
\enumerate{
\item Average annual catch across projection years and simulations per fleet and across fleets
\item Average interannual variation in catch (IAV) across projection years (n) per fleet and across fleets
\item \% of years in which the fishery is closed across simulations (s)
}

data.frame 2
4.	Average relative mean squared error in estimate of spawning biomass in the terminal year across simulations
5. \% of years in which the population is perceived as undergoing overfishing as determined from F_Limit across simulations via \code{\link{build_hcr}} in the EM
6.	\% of years in which the population is perceived to be overfished  as determined from B_Limit across simulations via \code{\link{build_hcr}} in the EM
7. \% of years in which the population is undergoing overfishing as determined from the “true” F_Limit across simulations via \code{\link{build_hcr}} in the OM
8. \% of years in which the population is overfished as determined from the “true” B_Limit across simulations via \code{\link{build_hcr}} in the OM
9.	Average ratio of spawning biomass over B_target in the terminal year across simulations in the OM
}
\description{
Management strategy evaluation performance metric summary
}
