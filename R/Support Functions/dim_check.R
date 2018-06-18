# Function to return turn matrices to vectors if a vector
dim_check <- function( dat ){
  # Matrix to vector
  dim_d <- dim(dat)
  if(class(dat) == "matrix" & 1 %in% dim_d & length(dim_d) == 2){
    dat <- as.vector(dat)
  }

  # List to matrix
  if( class(dat) == "list"){
    max_nrow <- max(sapply(dat, nrow))
    max_ncol <- max(sapply(dat, ncol))
    n_s <- length( dat )

    if( max_nrow == 1){
      mat <- matrix(NA, nrow = n_s, ncol = max_ncol)

      for(i in 1:n_s){
        for(c in 1:max_ncol){
          mat[i, c] <- dat[[i]][1, c]
        }
      }

      dat <- mat

    }
  }

  return(dat)
}


# Function to remove all columns to the right of a NA column including the NA column
remove_na_col <- function( dat ){

  if(class(dat) == "matrix"){
    na_col <- which(apply (is.na(dat), 2, all) == 1)
    if(length(na_col) != 0){
      dat <- dat[, 1:ncol(dat) < na_col  ]
    }
  }

  if(class(dat) == "list"){
    for(i in 1:length(dat)){
      na_col <- which(apply (is.na(dat[[i]]), 2, all) == 1)
      if(length(na_col) != 0){
        dat[[i]] <- dat[[i]][, 1:ncol(dat[[i]]) < na_col  ]
      }
    }
  }

  return(dat)
}


# Function to convert list to array
list_to_array <- function( dat ){

  if(class( dat ) == "list"){
    max_nrow <- max(sapply(dat, nrow))
    max_ncol <- max(sapply(dat, ncol))
    n_s <- length( dat )

    dat_array <- array(NA, dim = c(max_nrow, max_ncol, n_s))

    for(i in 1:n_s){
      for(x in 1:nrow(dat[[i]])){
        for(y in 1:ncol(dat[[i]])){
          dat_array[x, y, i] <- dat[[i]][x, y]
        }
      }
    }
    dat <- dat_array
  }

  return( dat )
}
