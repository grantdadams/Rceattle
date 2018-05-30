library(TMB)


setwd("syntax")
compile("lin_reg.cpp")
dyn.load(dynlib("lin_reg"))
setwd("../")



set.seed(123)
nobs <- 120

x <- seq(0, 10, length=nobs)
x2 <- seq(10, 20, length=nobs)
Y=rnorm(length(x), 0, 2) + 4 * x
plot(Y,x)
a1 = array(unlist(list(cbind(x, x2), cbind(x2, x))), dim = c(nobs, 2, 2))

data <- list(v1 = Y, m1 = cbind(x, x2), iv1 = round(Y, 0), a1 = a1)
parameters <- list(a= c(0, 0), logSigma=0)
obj <- MakeADFun(data, parameters, DLL="lin_reg")
opt <- do.call("optim", obj)
rep <- sdreport(obj)
srep <- summary(sdreport(obj))
srep
