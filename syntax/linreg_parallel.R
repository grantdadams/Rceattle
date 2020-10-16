# https://github.com/kaskr/adcomp/blob/master/TMB/inst/examples/linreg_parallel.R
library(TMB)
compile("linreg_parallel.cpp")
dyn.load(dynlib("linreg_parallel"))

## Simulate data
set.seed(123)
x <- seq(0, 10, length = 50001)
data <- list(Y = rnorm(length(x)) + x, x = x)
parameters <- list(a=0, b=0, logSigma=0)

## Fit model
obj <- MakeADFun(data, parameters, DLL="linreg_parallel")
opt <- nlminb(obj$par, obj$fn, obj$gr)

opt$objective - obj$report(obj$env$last.par.best)$nll
opt$objective - obj$fn(obj$env$last.par.best)
