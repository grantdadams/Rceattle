
# Test customized vector, matrix and array operations in TMB.

library(TMB)


setwd("syntax")
compile("vector_matrix_arrays.cpp")
dyn.load(dynlib("vector_matrix_arrays"))
setwd("../")

set.seed(123)
nobs <- 120

# Create  vector, matrix, and array
x <- seq(0, 10, length=nobs)
x2 <- seq(10, 20, length=nobs)
Y=rnorm(length(x), 0, 2) + 4 * x
plot(Y,x)
a1 = array(unlist(list(cbind(x, x2), cbind(x2, x))), dim = c(nobs, 2, 2))
data <- list(v1 = Y, m1 = cbind(x, x2), iv1 = round(Y, 0), a1 = a1)

# Test it
obj <- MakeADFun(data, parameters = list(), type="Fun",
                 checkParameterOrder=FALSE, DLL="vector_matrix_arrays")

