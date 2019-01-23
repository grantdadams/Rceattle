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
