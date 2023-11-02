xvec <- seq(-0.002,0.002,length=601)

f <- function(x,eps=0.001) {
  eps*(1/(1-(x-eps)/eps + (x-eps)^2/eps^2 - (x-eps)^3/eps^3
          + (x-eps)^4/eps^4 - (x-eps)^5/eps^5))
}

f2 <- function(x,eps=0.001) {
  if(x < eps){
    xp <- -(x/eps-1)
    print(paste0("Penalty = ", (0.01)*(x-eps)^2))
    return(eps*(1/(1+xp+xp^2 +xp^3+xp^4+xp^5)))
  } else{
    return(x)
  }
}

all.equal(f(xvec),f2(xvec))

f_an <- Vectorize(
  function(x, eps=0.001, n=2){
    if(x<eps){
      eps / sum( (-1)^(0:n)*((x-eps)^(0:n))/(eps^(0:n)))
    }else{
      x
    }
  }
)
