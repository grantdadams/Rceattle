

mfexp <- function( val ){
  res <- c()
  for(i in 1:length(val)){
    if( val[i] > 60){
      res[i] <- exp( 60 ) * (1 + 2 * (val[i] - 60)) / (1 + val[i] - 60)
    }
    if( val[i] < -60){
      res[i] <- exp( 60 ) * (1 - val[i] - 60) / (1 + 2 * (-val[i] - 60))
    } else{
      res[i] <- exp( val[i] )
    }
  }
  return( res )
}
