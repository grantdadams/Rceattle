library(RTMB)

## Consider a case where the same objective function has to be evaluated for many different identically shaped datasets, but the parameter vector changes
# Adapted from
# https://github.com/kaskr/RTMB/blob/bc14acc8d45bef0f9087f7cd9ade7f30be2533c2/tmb_examples/dataeval.R#L5


## Linear regressions with year specific terms and 10 replicates, but peeled 10 years
nyrs <- 20
ndat <- 10
nrep <- 10

dfs <- list()

for(i in 1:nrep){
  dfs[[i]] <- data.frame(year=rep(1:nyrs, each = ndat))
  dfs[[i]]$y <- dfs[[i]]$year + rnorm(ndat * nyrs) # Add noise
  dfs[[i]]$use <- ifelse(dfs[[i]]$year <= (nyrs - i + 1), 1, 0) # if on or before peel year, use
}

model_mats <- lapply(dfs, function(x) as.matrix(model.matrix(y ~ as.factor(year), data = x))) # Convert to factor variable and create model matrix


## Setup a parameter list for the first data
plist <- as.relistable(list(b=rep(0, nyrs), sd=1, i=1))

## Objective function for one chunk of data
F <- MakeTape(function(x) {
  ## Get parameter list
  parms <- relist(x)

  ## Fetch data from R (x and y)
  x <- DataEval( function(i) model_mats[[i]], parms$i)
  y <- DataEval( function(i) dfs[[i]]$y, parms$i)
  use <- DataEval( function(i) dfs[[i]]$use, parms$i) # Use data or no (e.g. peeled)

  ## NLL
  nll <- 0
  for(j in 1:length(y)){
    if(use[j] == 1){ # Throws error
      nll = nll - dnorm(y[j], x[j,] %*% parms$b, parms$sd, log=TRUE)
    }
  }
  b_off <- 1:length(parms$b) > (nyrs - i + 1)
  nll = nll + parms$b[b_off]^2 # Shrink to zero

  ## Add parameter penalty
  return(nll)
}, unlist(plist))

F$simplify() ## Simplify the operation sequence if possible
F <- F$atomic() ## Turn into atomic function

## Aggregate tape by summing over 'i'
obj <- MakeADFun(function(x) {
  s <- 0
  for (i in 1:10) s <- s + F(c(x,i))
  s
}, unlist(plist)[1:(nyrs + 1)])

opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
sdreport(obj)
