#' Generate Length-at-Age Transition Matrix
#'
#' This function calculates a probability transition matrix that defines the
#' probability of a fish of a given age belonging to specific length bins.
#' It supports Von Bertalanffy and Richards growth models and includes
#' a Stock Synthesis (SS) style plus-group correction.
#'
#' @param fracyr Numeric. Fraction of the year (0 = Jan 1st).
#' @param nsex_sp Integer. Number of sexes for the species.
#' @param nages_sp Integer. Number of age classes.
#' @param nlengths_sp Integer. Number of length bins.
#' @param max_nlengths Integer. Global maximum for length bins (for array sizing).
#' @param nyrs Integer. Number of years in the simulation.
#' @param lengths_sp Vector. Boundaries of the length bins.
#' @param minage_sp Numeric. The reference age (L1) for growth estimation.
#' @param maxage_sp Numeric. The age at which growth enters the asymptotic phase.
#' @param growth_params_sp Array. Dimensions (sex, yr, 4).
#'   Params: [,,1]: K, [,,2]: L1, [,,3]: Linf, [,,4]: Richards m.
#' @param growth_ln_sd_sp Array. Dimensions (sex, 2).
#'   Log-SD of length: [s, 1] is SD at minage, [s, 2] is SD at maxage.
#' @param growth_model_sp Integer. 1 = Von Bertalanffy, 2 = Richards.
#'
#' @return A 4D array of probabilities with dimensions (sex, age, length, year).
#' @export
get_growth_matrix_r <- function(fracyr, nsex_sp, nages_sp, nlengths_sp, max_nlengths, nyrs,
                                lengths_sp, minage_sp, maxage_sp,
                                growth_params_sp, growth_ln_sd_sp, growth_model_sp) {

  # Define names for the dimensions
  dim_names <- list(
    sex    = paste0("Sex_", 1:nsex_sp),
    age    = paste0("Age_", 1:nages_sp),
    length = paste0("Len_", lengths_sp),
    year   = paste0("Year_", 1:nyrs)
  )

  # Initialize Output: (sex, age, ln, yr)
  growth_matrix <- array(0, dim = c(nsex_sp, nages_sp, nlengths_sp, nyrs),
                         dimnames = dim_names)
  length_at_age <- array(0, dim = c(nsex_sp, nages_sp, nyrs),
                         dimnames =   list(
                           sex    = paste0("Sex_", 1:nsex_sp),
                           age    = paste0("Age_", 0:(nages_sp - 1)),
                           year   = paste0("Year_", 1:nyrs)
                         ))
  length_sd     <- array(0, dim = c(nsex_sp, nages_sp, nyrs))

  l_min <- lengths_sp[1]
  l_max <- lengths_sp[nlengths_sp]

  for(s in 1:nsex_sp) {
    for(y in 1:nyrs) {
      # --- 1. Calculate Mean Length at Age ---
      # Params: 1:K, 2:L1, 3:Linf, 4:Richards_m
      k    <- growth_params_sp[s, y, 1]
      l1   <- growth_params_sp[s, y, 2]
      linf <- growth_params_sp[s, y, 3]

      b_len <- (l1 - l_min) / minage_sp

      for(a in 1:nages_sp) {
        current_age <- a + fracyr

        if (growth_model_sp == 1) { # VB
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            if(y == 1){
              length_at_age[s, a, y] <- linf + (l1 - linf) * (exp(-k * (current_age - minage_sp)))
            } else{
              if(a == nages_sp){ # linear growth + growth equation
                last_linear = l_min + b_len * minage_sp # last age (cont) with linear growth

                length_at_age[s, a, y] = last_linear + (last_linear - linf) * (exp(-k * (current_age - minage_sp)) - 1.0)
              }else{
                length_at_age[s, a, y] <- length_at_age[s, a-1, y-1] + (length_at_age[s, a-1, y-1] - growth_params_sp[s, y-1, 3]) * (exp(-growth_params_sp[s, y-1, 1]) - 1)
              }
            }
          }
        } else if (growth_model_sp == 2) { # Richards
          m <- growth_params_sp[s, y, 4]
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            if(y == 1){
              length_at_age[s, a, y] <- (linf^m + (l1^m - linf^m) * (exp(-k * (current_age - minage_sp))))^(1/m)
            } else{
              if(a == nages_sp){ # linear growth + growth equation
                last_linear = l_min + b_len * minage_sp # last age (cont) with linear growth

                length_at_age[s, a, y] = (last_linear^m + (last_linear^m - linf^m) * (exp(-k * (current_age - minage_sp)) - 1.0))^(1/m)
              }else{
                lagk <- growth_params_sp[s, y-1, 2]
                lagm <- growth_params_sp[s, y-1, 4]
                laglinf <-  growth_params_sp[s, y-1, 3]
                length_at_age[s, a, y] <- (length_at_age[s, a-1, y-1]^lagm + (length_at_age[s, a-1, y-1]^lagm - laglinf^lagm) * (exp(-lagk) - 1))^1/lagm
              }
            }
          }
        }

        # --- 2. Plus Group Correction (SS Style) ---
        if(a == nages_sp) {
          diff <- growth_params_sp[s, y, 3] - length_at_age[s, a, y]
          ages <- 0:(nages_sp-1)
          weight_a <- exp(-0.2 * ages)
          vals <- length_at_age[s, a, y] + (ages / nages_sp) * diff
          length_at_age[s, a, y] <- sum(vals * weight_a) / sum(weight_a)
        }

        # --- 3. SD Calculation ---
        sd1 <- exp(growth_ln_sd_sp[s, 1])
        sda <- exp(growth_ln_sd_sp[s, 2])

        if(current_age < minage_sp) {
          length_sd[s, a, y] <- sd1
        } else if(a == nages_sp) {
          length_sd[s, a, y] <- sda
        } else {
          slope <- (sda - sd1) / (linf - l1) # Match C++ interpolation
          length_sd[s, a, y] <- sd1 + slope * (length_at_age[s, a, y] - l1)
        }

        # --- 4. Matrix Distribution ---
        for(l in 1:nlengths_sp) {
          if(l == 1) {
            fac1 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- pnorm(fac1)
          } else if(l == nlengths_sp) {
            fac1 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- 1 - pnorm(fac1)
          } else {
            fac1 <- (lengths_sp[l+1] - length_at_age[s, a, y]) / length_sd[s, a, y]
            fac2 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- pnorm(fac1) - pnorm(fac2)
          }
        }
      }
    }
  }

  return(list(length_at_age = length_at_age, growth_matrix = growth_matrix))
}


get_growth_matrix_test <- function(fracyr, nsex_sp, nages_sp, nlengths_sp, max_nlengths, nyrs,
                                   lengths_sp, minage_sp, maxage_sp,
                                   growth_params_sp, growth_ln_sd_sp, growth_model_sp) {

  # Define names for the dimensions
  dim_names <- list(
    sex    = paste0("Sex_", 1:nsex_sp),
    age    = paste0("Age_", 1:nages_sp),
    length = paste0("Len_", lengths_sp),
    year   = paste0("Year_", 1:nyrs)
  )

  # Initialize Output: (sex, age, ln, yr)
  growth_matrix <- array(0, dim = c(nsex_sp, nages_sp, nlengths_sp, nyrs),
                         dimnames = dim_names)
  length_at_age <- array(0, dim = c(nsex_sp, nages_sp, nyrs),
                         dimnames =   list(
                           sex    = paste0("Sex_", 1:nsex_sp),
                           age    = paste0("Age_", 0:(nages_sp - 1)),
                           year   = paste0("Year_", 1:nyrs)
                         ))
  length_sd     <- array(0, dim = c(nsex_sp, nages_sp, nyrs))

  l_min <- lengths_sp[1]
  l_max <- lengths_sp[nlengths_sp]

  for(s in 1:nsex_sp) {
    for(y in 1:nyrs) {
      # --- 1. Calculate Mean Length at Age ---
      # Params: 1:K, 2:L1, 3:Linf, 4:Richards_m
      k    <- growth_params_sp[s, y, 1]
      l1   <- growth_params_sp[s, y, 2]
      linf <- growth_params_sp[s, y, 3]

      b_len <- (l1 - l_min) / minage_sp

      for(a in 1:nages_sp) {
        current_age <- a + fracyr

        if (growth_model_sp == 1) { # VB
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            length_at_age[s, a, y] <- linf + (l1 - linf) * (exp(-k * (current_age - minage_sp)))
          }
        } else if (growth_model_sp == 2) { # Richards
          m <- growth_params_sp[s, y, 4]
          if(current_age <= minage_sp) {
            length_at_age[s, a, y] <- l_min + b_len * current_age
          } else {
            length_at_age[s, a, y] <- (linf^m + (l1^m - linf^m) * (exp(-k * (current_age - minage_sp))))^(1/m)
          }
        }

        # --- 2. Plus Group Correction (SS Style) ---
        if(a == nages_sp) {
          diff <- growth_params_sp[s, y, 3] - length_at_age[s, a, y]
          ages <- 0:(nages_sp-1)
          weight_a <- exp(-0.2 * ages)
          vals <- length_at_age[s, a, y] + (ages / nages_sp) * diff
          length_at_age[s, a, y] <- sum(vals * weight_a) / sum(weight_a)
        }

        # --- 3. SD Calculation ---
        sd1 <- exp(growth_ln_sd_sp[s, 1])
        sda <- exp(growth_ln_sd_sp[s, 2])

        if(current_age < minage_sp) {
          length_sd[s, a, y] <- sd1
        } else if(a == nages_sp) {
          length_sd[s, a, y] <- sda
        } else {
          slope <- (sda - sd1) / (linf - l1) # Match C++ interpolation
          length_sd[s, a, y] <- sd1 + slope * (length_at_age[s, a, y] - l1)
        }

        # --- 4. Matrix Distribution ---
        for(l in 1:nlengths_sp) {
          if(l == 1) {
            fac1 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- pnorm(fac1)
          } else if(l == nlengths_sp) {
            fac1 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- 1 - pnorm(fac1)
          } else {
            fac1 <- (lengths_sp[l+1] - length_at_age[s, a, y]) / length_sd[s, a, y]
            fac2 <- (lengths_sp[l] - length_at_age[s, a, y]) / length_sd[s, a, y]
            growth_matrix[s, a, l, y] <- pnorm(fac1) - pnorm(fac2)
          }
        }
      }
    }
  }

  return(list(length_at_age = length_at_age, growth_matrix = growth_matrix))
}



#' Calculate Predicted Weight-at-Age
#'
#' Converts a growth matrix (length-at-age probabilities) into mean weight-at-age
#' using a length-weight relationship (W = a * L^b).
#'
#' @param nsex_sp Integer. Number of sexes.
#' @param nages_sp Integer. Number of age classes.
#' @param nlengths_sp Integer. Number of length bins.
#' @param nyrs Integer. Number of years.
#' @param lengths_sp Vector. Boundaries of the length bins.
#' @param growth_matrix Array. 4D array (sex, age, length, year) from get_growth_matrix_r.
#' @param lw_params Array. Dimensions (sex, yr, 2).
#'   Params: [s, y, 1] is alpha (a), [s, y, 2] is beta (b).
#'
#' @details The function calculates midpoints for length bins to avoid bias.
#' For the first bin, it assumes the width is equal to the second bin's width.
#' The final weight-at-age is the expected value across all length bins for that age.
#'
#' @return A 3D array of mean weights with dimensions (sex, age, year).
#' @export
get_weight_at_age_r <- function(nsex_sp, nages_sp, nlengths_sp, nyrs,
                                lengths_sp, length_at_age, growth_matrix, lw_params) {
  # Define names for the dimensions
  dim_names <- list(
    sex  = paste0("Sex_", 1:nsex_sp),
    age  = paste0("Age_", 0:(nages_sp - 1)),
    year = paste0("Year_", 1:nyrs)
  )

  # Output: (sex, age, yr)
  waa <- array(0, dim = c(nsex_sp, nages_sp, nyrs),
               dimnames = dim_names)


  for(s in 1:nsex_sp) {
    for(y in 1:nyrs) {
      alpha <- lw_params[s, y, 1]
      beta  <- lw_params[s, y, 2]

      # Weight at length for all bins
      wal <- alpha * lengths_sp ^ beta

      for(a in 1:nages_sp) {
        # Matrix multiply: Prob(length | age) * Weight(length)
        waa[s, a, y] <- sum(growth_matrix[s, a, , y] * wal)
      }
    }
  }
  return(waa)
}




# 1. Setup Mock Data
nyrs <- 5
nages <- 10
nlengths <- 20
lengths <- seq(10, 100, length.out = nlengths)

# Growth params: (sex, yr, param) -> params: K, L1, Linf, m
gp <- array(rep(c(0.2, 20, 90, 1.0), each = nyrs), dim=c(1, nyrs, 4))
gp[,,3] <- runif(nyrs, 60, 100)
gp[,,1] <- runif(nyrs, 0.1, 0.3)
gsd <- array(log(c(2, 5)), dim=c(1, 2))

# 2. Run Growth Matrix
gm <- get_growth_matrix_r(fracyr=0,
                          nsex_sp=1,
                          nages_sp=nages,
                          nlengths_sp=nlengths,
                          max_nlengths=nlengths,
                          nyrs=nyrs,
                          lengths_sp=lengths,
                          minage_sp=1,
                          maxage_sp=10,
                          growth_params_sp=gp,
                          growth_ln_sd_sp=gsd,
                          growth_model_sp=1)


gmtest <- get_growth_matrix_test(fracyr=0,
                          nsex_sp=1,
                          nages_sp=nages,
                          nlengths_sp=nlengths,
                          max_nlengths=nlengths,
                          nyrs=nyrs,
                          lengths_sp=lengths,
                          minage_sp=1,
                          maxage_sp=10,
                          growth_params_sp=gp,
                          growth_ln_sd_sp=gsd,
                          growth_model_sp=1)

# 3. Visualize a specific age/year distribution
par(mfrow = c(1,2))
plot(lengths, gm$growth_matrix[1, 5, , 2], type="h", main="Length distribution for Age 5")
plot(lengths, gmtest$growth_matrix[1, 5, , 2], type="h", main="Length distribution for Age 5")

# 4. Calculate Weight-at-Age
lw_p <- array(c(0.00001, 3.0), dim=c(1, nyrs, 2)) # a=0.00001, b=3.0
waa <- get_weight_at_age_r(nsex_sp = 1, nages_sp = nages, nlengths_sp = nlengths, nyrs = nyrs,
                           lengths_sp = lengths, growth_matrix = gm$growth_matrix, lw_params = lw_p)

print(waa[1, , 1]) # Weights for all ages in year 1
