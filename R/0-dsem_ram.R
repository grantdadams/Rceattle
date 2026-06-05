# -----------------------------------------------------------------------------
# Vendored DSEM sem -> TMB-inputs pipeline
# -----------------------------------------------------------------------------
# These functions are ported (lightly adapted) from the `dsem` package
#   James-Thorson-NOAA/dsem, version 2.0.1, commit 81d3d817cc53284018ff6c3f620bd6d8bbe6ab6c
#   Thorson, J.T., Andrews, A.G., Essington, T., & Large, S. (2024).
#   Dynamic structural equation models. Methods Ecol. Evol. 15(4): 744-755.
#   dsem is distributed under GPL-3; portions of make_dsem_ram / parse_path /
#   classify_variables originate from package `sem` (GPL >= 2, used with
#   permission from John Fox).
#
# Why vendored: previously Rceattle called `dsem::dsem(run_model = FALSE)` at
# fit time purely to harvest `$tmb_inputs` / `$sem_full` internals and splice
# them into CEATTLE. That coupled Rceattle to undocumented dsem internals and
# left the vendored C++ (src/TMB/dsem.hpp) able to silently disagree with
# whatever dsem version happened to be installed. Owning this pipeline gives a
# single source of truth that is matched to src/TMB/dsem.hpp.
#
# The functions below are internal (not exported). `build_dsem_inputs()` is the
# entry point used by build_dsem_objects(); it returns the same `sem_full` /
# `tmb_inputs` contract that `dsem::dsem(run_model = FALSE)` did.
# -----------------------------------------------------------------------------

#' Parse a single SEM path specification
#'
#' Ported from \code{dsem::parse_path} (originally \code{sem::parse.path}).
#' @param path text to parse
#' @return tagged list with \code{first}, \code{second}, \code{direction}
#' @keywords internal
#' @noRd
parse_path <- function( path ){
  path.1 <- gsub("-", "", gsub(" ", "", path))
  direction <- if(regexpr("<>", path.1) > 0){
    2
  }else if(regexpr("<", path.1) > 0){
    -1
  }else if(regexpr(">", path.1) > 0){
    1
  }else{
    stop(paste("ill-formed path:", path))
  }
  path.1 <- strsplit(path.1, "[<>]")[[1]]
  list(first = path.1[1], second = path.1[length(path.1)], direction = direction)
}

#' Classify SEM variables as exogenous or endogenous
#'
#' Ported from \code{dsem::classify_variables} (originally \code{sem:::classifyVariables}).
#' @param model SEM model matrix/data.frame (first column = paths)
#' @return tagged list with \code{exogenous} and \code{endogenous}
#' @keywords internal
#' @noRd
classify_variables <- function( model ){
  variables <- logical(0)
  for (paths in model[, 1]) {
    vars <- gsub(pattern = " ", replacement = "", x = paths)
    vars <- sub("-*>", "->", sub("<-*", "<-", vars))
    if (grepl("<->", vars)) {
      vars <- strsplit(vars, "<->")[[1]]
      if (is.na(variables[vars[1]])) variables[vars[1]] <- FALSE
      if (is.na(variables[vars[2]])) variables[vars[2]] <- FALSE
    } else if (grepl("->", vars)) {
      vars <- strsplit(vars, "->")[[1]]
      if (is.na(variables[vars[1]])) variables[vars[1]] <- FALSE
      variables[vars[2]] <- TRUE
    } else if (grepl("<-", vars)) {
      vars <- strsplit(vars, "<-")[[1]]
      if (is.na(variables[vars[2]])) variables[vars[2]] <- FALSE
      variables[vars[1]] <- TRUE
    } else stop("incorrectly specified model")
  }
  list(endogenous = names(variables[variables]),
       exogenous  = names(variables[!variables]))
}

#' Make a RAM (Reticular Action Model) from SEM arrow-and-lag notation
#'
#' Ported verbatim from \code{dsem::make_dsem_ram} (v2.0.1). Returns a list of
#' class \code{dsem_ram} with elements \code{model}, \code{ram}, \code{variables},
#' \code{times}. The \code{ram} data.frame has columns
#' \code{heads, to, from, parameter, start, to_t, to_j}. NA parameter values are
#' encoded as \code{-1} (sentinel) so the matrix passes cleanly to
#' \code{DATA_IMATRIX} in TMB.
#'
#' @param sem SEM specification text (arrow-and-lag notation)
#' @param times numeric vector of times in order
#' @param variables character vector of variable names
#' @param covs variables for which to estimate a standard deviation
#' @param quiet suppress messages
#' @param remove_na retained for signature compatibility (unused; see note)
#' @keywords internal
#' @noRd
make_dsem_ram <- function( sem,
                           times,
                           variables,
                           covs = variables,
                           quiet = FALSE,
                           remove_na = TRUE ){

  ####### Error checks
  if( !is.numeric(times) ) stop("`times` must be numeric in `make_dsem_ram`")

  ####### Define local functions
  match_row = function( df, x ) which( df[1]==x[1] & df[2]==x[2] )

  add.variances <- function() {
    variables <- need.variance()
    nvars <- length(variables)
    if (nvars == 0) return(model)
    message("NOTE: adding ", nvars, " variances to the model")
    paths <- character(nvars)
    par.names <- character(nvars)
    for (i in 1:nvars) {
      paths[i] <- paste(variables[i], "<->", variables[i])
      par.names[i] <- paste("V[", variables[i], "]", sep = "")
    }
    model.2 <- cbind(
      'path'  = c(model[, 1], paths),
      'lag'   = c(model[, 2], rep(0, nvars)),
      'name'  = c(model[, 3], par.names),
      'start' = c(model[, 4], rep(NA, length(paths))) )
    model.2
  }
  need.variance <- function() {
    all.vars <- classify_variables(model)
    exo.vars <- all.vars$exogenous
    end.vars <- all.vars$endogenous
    variables <- logical(0)
    for (i in seq_len(nrow(model))) {
      paths = model[i, 1]
      lag = model[i, 2]
      vars <- gsub(pattern = " ", replacement = "", x = paths)
      vars <- sub("-*>", "->", sub("<-*", "<-", vars))
      vars <- sub("<->|<-", "->", vars)
      vars <- strsplit(vars, "->")[[1]]
      if ((vars[1] != vars[2]) | (lag != 0)) {
        for (a.variable in vars) {
          if (is.na(variables[a.variable])) variables[a.variable] <- TRUE
        }
      } else {
        variables[vars[1]] <- FALSE
      }
    }
    if (!exog.variances && length(exo.vars) > 0) variables[exo.vars] <- FALSE
    if (!endog.variances && length(end.vars) > 0) variables[end.vars] <- FALSE
    names(variables)[variables]
  }

  ####### Step 1 -- parse to data.frame
  model = scan(
    text = sem,
    what = list(path = "", lag = 1, par = "", start = 1, dump = ""),
    sep = ",",
    strip.white = TRUE,
    comment.char = "#",
    fill = TRUE,
    quiet = quiet
  )
  model$path <- gsub("\\t", " ", model$path)
  model$par[model$par == ""] <- NA
  model <- data.frame( "path" = model$path, "lag" = model$lag,
                       "name" = model$par, "start" = model$start)

  # Adding a SD automatically
  if( !is.null(covs) ){
    for (cov in covs) {
      vars <- strsplit(cov, "[ ,]+")[[1]]
      nvar <- length(vars)
      for (i in 1:nvar) {
        for (j in i:nvar) {
          p1 = paste(vars[i], "<->", vars[j])
          p2 = if (i==j) paste("V[", vars[i], "]", sep = "") else paste("C[", vars[i], ",", vars[j], "]", sep = "")
          p3 = NA
          row <- data.frame("path"=p1, "lag"=0, "name"=p2, "start"=p3)
          if( isTRUE(any((row[1] %in% model[,1]) & (row[2] %in% model[,2]))) ){
            next
          }else{
            model <- rbind(model, row, deparse.level = 0, make.row.names = FALSE)
          }
        }
      }
    }
  }

  exog.variances = endog.variances = TRUE
  model = add.variances()

  ####### Step 2 -- Make RAM
  Q_names = expand.grid( times, variables )
  ram = NULL  # heads, to, from, parameter

  # Deal with fixed values
  par.names = model[, 3]
  # Exclude varying slopes
  par.names = ifelse( par.names %in% variables, NA, par.names )
  pars = stats::na.omit(unique(par.names))

  if( length(pars)==0 ){
    par.nos = rep(0, nrow(model))
  }else{
    par.nos = apply(outer(pars, par.names, "=="), 2, which)
    par.nos = unlist(sapply( par.nos, FUN=function(x) ifelse(length(x)==0, 0, x) ))
  }
  model = cbind( model, "parameter"=par.nos )
  startvalues = model[,4]

  # Add incidence to model
  model = cbind( model, first=NA, second=NA, direction=NA )
  for( i in seq_len(nrow(model)) ){
    path = parse_path(model[i,1])
    model[i,c('first','second','direction')] = unlist( path[c('first','second','direction')] )
  }

  # ADD LOGIC FOR PATHS (varying slopes)
  which_rows = model[,'name'] %in% variables
  if( any(which_rows) ){
    model[which_rows,'direction'] = 0
  }

  # Loop through paths
  B_kk = G_kk = P_kk = Matrix::drop0(Matrix::sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*length(times),2) ))
  for( i in seq_len(nrow(model)) ){
    # Time-lag matrix ... transpose if negative lag
    lag = as.numeric(model[i,2])
    L_tt = Matrix::sparseMatrix( i = seq( abs(lag)+1, length(times)),
                                 j = seq(1, length(times)-abs(lag)),
                                 x = 1,
                                 dims = rep(length(times),2) )
    if(lag<0) L_tt = Matrix::t(L_tt)

    # Interaction matrix
    P_jj = Matrix::sparseMatrix( i = match(model[i,'second'],variables),
                                 j = match(model[i,'first'],variables),
                                 x = 1,
                                 dims = rep(length(variables),2) )

    # Combine them
    tmp_kk = Matrix::kronecker(P_jj, L_tt)
    if( abs(as.numeric(model[i,'direction'])) == 0 ){
      B_kk = B_kk + tmp_kk * match(model[i,'name'],variables)
    }
    if( abs(as.numeric(model[i,'direction'])) == 1 ){
      P_kk = P_kk + tmp_kk * i
    }
    if( abs(as.numeric(model[i,'direction'])) == 2 ){
      G_kk = G_kk + tmp_kk * i
    }
  }

  # Convert to triplet.
  # NA / NA_integer_ triggers SAN errors when passed to DATA_IMATRIX, so use -1.
  f = function( x, first_column = 1){
    triplet = Matrix::mat2triplet(x)
    if( length(triplet$x)>0 ){
      out = data.frame( first_column, triplet$i, triplet$j, triplet$x, -1, -1 )
    }else{
      out = data.frame( numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0) )
    }
    names(out) = seq_len(ncol(out))
    return(out)
  }
  # Convert to triplet for spatially/temporally varying slope
  f2 = function( x ){
    triplet = Matrix::mat2triplet(x)
    if( length(triplet$x)>0 ){
      t_k = rep( seq_along(times), length(variables) )[triplet$i]
      out = data.frame( 0, triplet$i, triplet$j, -1, t_k, triplet$x )
    }else{
      out = data.frame( numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0) )
    }
    names(out) = seq_len(ncol(out))
    return(out)
  }

  # Combine:  heads(type), to(row), from(col), parameter(beta_z index),
  #           startvalue, to_t (varying path row of tsdata), to_j (varying path col)
  ram = rbind(
    f(P_kk, 1),
    f(G_kk, 2),
    f2(B_kk)
  )
  tmp = ifelse( ram[,4] < 0, -NA, ram[,4] )
  ram = data.frame(
    ram[,1:3, drop=FALSE],
    as.numeric(par.nos)[tmp],
    as.numeric(startvalues)[tmp],
    ram[,5:6, drop=FALSE]
  )
  # swap out NAs to pass to DATA_IMATRIX
  ram[,4] = ifelse( is.na(ram[,4]), -1, ram[,4] )
  colnames(ram) = c( "heads", "to", "from", "parameter", "start", "to_t", "to_j")

  out = list( "model"=model,
              "ram"=ram,
              "variables" = variables,
              "times" = times )
  class(out) = "dsem_ram"
  return(out)
}

#' Build TMB inputs (data / parameters / random / map) for the DSEM module
#'
#' Ported from the input-assembly portion of \code{dsem::dsem()} (v2.0.1, the
#' \code{build_model = FALSE} path). Replaces the previous
#' \code{dsem::dsem(run_model = FALSE)} call inside \code{build_dsem_objects()}.
#' Returns the same contract Rceattle consumed before:
#' \code{list(sem_full = <model df>, ram = <ram df>, tmb_inputs = list(data,
#' parameters, random, map))}.
#'
#' Defaults mirror dsem 2.0.1 \code{dsem_control()} together with the controls
#' Rceattle always set (\code{use_REML = FALSE}): gmrf_parameterization
#' \code{"gmrf_project"} (\code{options[1]==3}), conditional variance
#' (\code{options[2]==0}), \code{stabilize_Q = FALSE} (\code{options[3]==0}),
#' and \code{Random = "x_tj"}.
#'
#' @param sem SEM specification text
#' @param tsdata a \code{ts} object (time x variable) of DSEM data
#' @param family character vector (one per column of tsdata): fixed/normal/bernoulli/poisson/gamma
#' @param gmrf_parameterization one of "gmrf_project","full","project","mvn_project"
#' @param constant_variance one of "conditional","marginal","diagonal"
#' @param stabilize_Q add 1e-10 to precision diagonal (logical)
#' @param estimate_delta0 estimate initial-condition offsets (logical)
#' @param estimate_mu character vector of variables for which to estimate a mean (NULL = auto)
#' @param covs variables for which to add a SD term
#' @param use_REML treat mu_j as random in addition to x_tj (logical)
#' @param project_k logical vector of variables to project (NULL = auto: zero-variance vars)
#' @param quiet suppress messages
#' @keywords internal
#' @noRd
build_dsem_inputs <- function( sem,
                               tsdata,
                               family = rep("fixed", ncol(tsdata)),
                               gmrf_parameterization = "gmrf_project",
                               constant_variance = "conditional",
                               stabilize_Q = FALSE,
                               estimate_delta0 = FALSE,
                               estimate_mu = NULL,
                               covs = colnames(tsdata),
                               use_REML = FALSE,
                               project_k = NULL,
                               quiet = TRUE ){

  if( isFALSE(inherits(tsdata, "ts")) ) stop("`tsdata` must be a `ts` object")

  # Build RAM
  out = make_dsem_ram(
    sem,
    times = as.numeric(stats::time(tsdata)),
    variables = colnames(tsdata),
    covs = covs,
    quiet = quiet
  )
  ram = out$ram

  # Error checks (mirror dsem::dsem)
  if( isTRUE(any((out$model[,'direction']==2) & (out$model[,2]!=0))) ){
    stop("All two-headed arrows should have lag=0")
  }
  if( isFALSE(all(c(out$model[,'first'],out$model[,'second']) %in% colnames(tsdata))) ){
    stop("Some variable in `sem` is not in `tsdata`")
  }
  if( isFALSE(ncol(tsdata) == length(unique(colnames(tsdata)))) ){
    stop("Please check `colnames(tsdata)` to confirm that all variables (columns) have a unique name")
  }

  # options vector (see src/TMB/dsem.hpp)
  #  [1] parameterization: full=0, project=1, mvn_project=2, gmrf_project=3
  #  [2] constant_variance: conditional=0, marginal=1, diagonal=2
  #  [3] stabilize_Q: 0/1
  options = c(
    switch(gmrf_parameterization, "full"=0, "project"=1, "mvn_project"=2, "gmrf_project"=3, NA),
    switch(constant_variance, "conditional"=0, "marginal"=1, "diagonal"=2),
    ifelse( isTRUE(stabilize_Q), 1, 0 )
  )

  # Define variables to project to `project_k` unless provided by user
  if( is.null(project_k) ){
    project_k = sapply(
      seq_len(prod(dim(tsdata))),
      FUN = function(k){
        tmp = subset( ram, (ram$heads==2) & (ram$to==k) & (ram$from==k) )
        all( (tmp$parameter==0) & (tmp$start==0) )
      }
    )
  }else{
    project_k = as.vector(project_k)
  }

  # Identify variables that serve as moderators
  moderator_tj = array( FALSE, dim = dim(tsdata) )
  tmp = ram[ which(ram[,1] %in% c(0)), 6:7]
  moderator_tj[ as.matrix(tmp) ] = 1
  moderator_k = as.vector(moderator_tj)

  # Error checks tied to parameterization
  if( any(subset(ram, ram$heads==2 & !is.na(ram$start))$start==0) & (options[1]==0) ){
    stop("Cannot use exogenous variance of zero using gmrf_parameterization=`full`")
  }
  if( any(project_k & moderator_k) ){
    if( gmrf_parameterization %in% c("mvn_project", "gmrf_project") ){
      stop("Cannot use gmrf_parameterization == `project` when using moderated variables")
    }
  }

  # Data
  Data = list(
    "options" = options,
    "RAM" = as.matrix(ram[,-5]),          # heads,to,from,parameter,to_t,to_j (drops 'start')
    "RAMstart" = as.numeric(ram[,5]),
    "familycode_j" = sapply(family, FUN=switch, "fixed"=0, "normal"=1, "bernoulli"=2, "poisson"=3, "gamma"=4 ),
    "y_tj" = tsdata
  )
  # obs_idx / unobs_idx are only used by the C++ for options 2 and 3, but
  # ceattle_v01_11.cpp always declares them via DATA_IVECTOR, so always provide
  # them (harmless for options 0/1). Forced integer for DATA_IVECTOR.
  Data$obs_idx   = as.integer(which( !project_k ) - 1L)   # CPP (0-based) indexing
  Data$unobs_idx = as.integer(which(  project_k ) - 1L)

  # Parameters
  Params = list(
    beta_z    = rep(0, max(ram[,4], na.rm=TRUE)),
    lnsigma_j = rep(0, ncol(tsdata)),
    mu_j      = rep(0, ncol(tsdata)),
    delta0_j  = rep(0, ncol(tsdata)),
    x_tj      = ifelse( is.na(tsdata), 0, tsdata )
  )
  if( estimate_delta0 == FALSE ){
    Params$delta0_j = numeric(0)
  }
  # Scale starting values: higher for two-headed than one-headed arrows
  which_nonzero = which(ram[,4] > 0)
  beta_type = tapply( ram[which_nonzero,1], INDEX=ram[which_nonzero,4], max)
  Params$beta_z = ifelse(beta_type==1, 0.01, 1)
  # Override starting values if supplied in the sem
  start_z = tapply( as.numeric(ram[which_nonzero,5]), INDEX=ram[which_nonzero,4], mean )
  Params$beta_z = ifelse( is.na(start_z), Params$beta_z, start_z)

  # Process estimate_mu
  if( is.null(estimate_mu) ){
    estimate_mu = stats::na.omit( ifelse(colSums(!is.na(tsdata))==0, NA, colnames(tsdata)) )
  }

  # Map
  Map = list()
  Map$x_tj = ifelse( is.na(as.vector(tsdata)) | (Data$familycode_j[col(tsdata)] %in% c(1,2,3,4)),
                     seq_len(prod(dim(tsdata))), NA )
  Map$lnsigma_j = factor( ifelse(Data$familycode_j %in% c(0,2,3), NA, seq_along(Params$lnsigma_j)) )
  Map$mu_j = factor( ifelse(colnames(tsdata) %in% estimate_mu, seq_len(ncol(tsdata)), NA) )
  if( (options[1] == 3) & (sum(project_k) > 0) ){
    Map$x_tj = ifelse( seq_along(Map$x_tj) %in% (Data$unobs_idx + 1), NA, Map$x_tj )
  }
  Map$x_tj = factor(Map$x_tj)

  # Random
  Random = if( isTRUE(use_REML) ) c("x_tj","mu_j") else "x_tj"

  list(
    sem_full   = out$model,
    ram        = ram,
    tmb_inputs = list("data"=Data, "parameters"=Params, "random"=Random, "map"=Map)
  )
}
