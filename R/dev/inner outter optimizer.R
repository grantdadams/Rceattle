# testing doing it manually
obj2 <- RTMB::MakeADFun(func=obj$env$data,
                        parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
f <- obj2$fn
g <- obj2$gr

# arbitrary point in space
theta <- c(pars$mu, pars$logtau)
u <- rep(0,8)
rand <- 1:8
finner <- function(u,theta) f(c(u,theta))
ginner <- function(u,theta) g(c(u,theta))[rand]
finner(u, theta=theta)
ginner(u, theta=theta)
uhat <- function(theta) {
  opt <- nlminb(start = u, objective = finner, gradient = ginner, theta=theta)
  if(opt$convergence!=0) warning("possibly failed convergence at ", paste(theta, collapse=', '))
  opt$par
}
fouter <- function(theta) {
  # first find uhat given theta
  utmp <- uhat(theta)
  gtmp <- function(u) ginner(u, theta)
  H <- numDeriv::jacobian(func = gtmp, x = utmp)
  -log(sqrt(2*pi)^8)+(1/2)*log(det(H)) +finner(utmp,theta)
}
fouter(theta)
obj$fn(theta)
gouter <- function(theta) {
  numDeriv::grad(fouter, x=theta)
}
gouter(theta)
obj$gr(theta)
