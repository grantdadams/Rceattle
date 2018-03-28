library(TMB)


setwd("src")
compile("lin_reg.cpp")
dyn.load(dynlib("lin_reg"))
setwd("../")



set.seed(123)
nobs <- 120

x <- seq(0, 10, length=nobs)
x2 <- seq(10, 20, length=nobs)
Y=rnorm(length(x), 0, 2) + 4 * x
plot(Y,x)

data <- list(Y=Y, x=array(unlist(list(cbind(x, x2), cbind(x2, x))), dim = c(nobs, 2, 2)))
parameters <- list(a=0, b=0, logSigma=0)
obj <- MakeADFun(data, parameters, DLL="lin_reg")
opt <- do.call("optim", obj)
rep <- sdreport(obj)
srep <- summary(sdreport(obj))
srep
