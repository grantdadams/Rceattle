#' Data inputs for single species CEATTLE of the Bering Sea from 1979 to 2017
#'
#' A dataset containing the inputs used for CEATTLE
#'
#' MODEL CONFIGURATION SECTION
#' @format A list of different integers and vectors, specifying model switches:
#' \describe{
#'   \item{debug}{Integer = Logical to debug (1) or not (0)}
#'   \item{msmMode}{Integer: Single species (0) or multi-species (1) mode}
#'   \item{est_diet}{Integer: include diet data in likelihood. 0 = no, 1 = yes}
#'   \item{suitMode}{Integer: Mode for suitability/functional calculation. 0 = empirical based on diet data (Holsman et al. 2015), 1 = gamma selectivity from Kinzey and Punt (2009)}
#'   \item{avgnMode}{Integer: Average numbers-at-age to be used in predation function. 0 = average numbers-at-age, 1 = \eqn{= N * e^{-Z/2} }, 2 = \eqn{N}.
#'   \item{random_rec}{Integer: Switch to estimate recruitment deviations as random effects. 0 = no, 1 = yes}
#'   \item{niter}{Integer: Number of loops around population/predation equations}
#'   \item{nspp}{Integer: number of species in model}
#'   \item{styr}{Integer: Start year}
#'   \item{nyrs}{Integer: Number of estimation years}
#'   \item{logist_sel_phase}{Vector: selectivity type. 0 = logistic, 1 = non-parametric selectivity estimated for each age; length = nspp:}
#'   \item{nselages}{Vector: if non-parametric selectivity, the number of ages to estimate selectivity; length = nspp:}
#' }
#'
#' #' DATA INPUTS SECTION
#' @format A list of different integers, vectors, matrices, and arrays. NOTE: Matrices and arrays will be filled with N
#' \describe{
#'   \item{nages}{Vector: Number of species ages; length = nspp}
#'   \item{nyrs_tc_biom}{Vector: Number of years with total observed catch; length = nspp}
#'   \item{yrs_tc_biom}{Matrix: Years with total observed catch; dim = [nspp, nyrs_tc_biom]}
#'   \item{tcb_obs}{Matrix: Observed total yield (kg); dim = [nspp, nyrs_tc_biom]}
#'   \item{}{}
#'   \item{}{}
#'   \item{}{}
#'   \item{}{}
#'
#'   ...
#' }
"BS2017SS"
