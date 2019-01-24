# Function to return turn matrices to vectors if a vector
dim_check <- function( dat ){
  # Matrix to vector
  dim_d <- dim(dat)
  if(class(dat) == "matrix" & 1 %in% dim_d & length(dim_d) == 2){
    dat <- as.vector(dat)
  }

  # List to matrix
  if( class(dat) == "list"){
    if(class(dat[[1]]) != "list"){
      max_nrow <- max(sapply(dat, nrow))
      max_ncol <- max(sapply(dat, ncol))
      n_s <- length( dat )

      if( max_nrow == 1){
        mat <- matrix(NA, nrow = n_s, ncol = max_ncol)

        for(i in 1:n_s){
          for(c in 1:ncol(dat[[i]])){
            mat[i, c] <- dat[[i]][1, c]
          }
        }

        dat <- mat
      }
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
      if(class(dat[[i]]) == "list"){
        for(j in 1:length(dat[[i]])){
          na_col <- which(apply (is.na(dat[[i]][[j]]), 2, all) == 1)
          if(length(na_col) != 0){
            dat[[i]][[j]] <- dat[[i]][[j]][, 1:ncol(dat[[i]][[j]]) < na_col  ]
          }
        }

      } else {
        na_col <- which(apply (is.na(dat[[i]]), 2, all) == 1)
        if(length(na_col) != 0){
          dat[[i]] <- dat[[i]][, 1:ncol(dat[[i]]) < na_col  ]
        }
      }
    }
  }

  return(dat)
}


# Function to convert list to array
list_to_array <- function( dat ){

  if(class( dat ) == "list"){
    if(class(dat[[1]]) == "list"){ # 4D - array

      n_s <- length( dat )
      n_t <- c()
      max_nrow <- c()
      max_ncol <- c()

      for(i in 1:n_s){
        n_t[i] <- length(dat[[i]])
        max_nrow[i] <- max(sapply(dat[[i]], nrow))
        max_ncol[i] <- max(sapply(dat[[i]], ncol))
      }

      dat_array <- array(NA, dim = c(n_s, max(n_t),max(max_nrow), max(max_ncol)))

      for(i in 1:n_s){
        for(j in 1:length(dat[[i]])){
          for(x in 1:nrow(dat[[i]][[j]])){
            for(y in 1:ncol(dat[[i]][[j]])){
              dat_array[i, j, x, y] <- dat[[i]][[j]][x, y]
            }
          }
        }
      }
      dat <- dat_array

    } else{ # 3D - array
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
  }

  return( dat )
}

# Function to condense dataframe if load input empty cells into a dataframe
condense <- function(df){
  nrow_df <- nrow(df)
  # Remove #
  df <- gsub("#", "", df)
  # Condense and remove NAs
  new_df <- df[1,which(!is.na(df[1,]))]
  for(i in 2:nrow_df){
    new_df <- rbind(new_df, df[i,which(!is.na(df[i,]))])
  }
  # Make numeric
  new_df <- matrix(as.numeric(as.character(new_df)), nrow = nrow(new_df), ncol = ncol(new_df))

  return(new_df)
}


# Safe exponent
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


