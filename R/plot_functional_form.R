#' Plot functional form
#'
#' @describtion Function to plot the functional form estimated or specified by \code{\link{Rceattle}}
#'
#' @param params Parameter list object from \code{\link{build_params}} or \code{\link{Rceattle}}
#' @param pred Predator index
#' @param pred_age Predator age
#' @param prey Prey index
#' @export
plot_form <- function( params = NULL, pred = 1, pred_age = 1, prey = 1, msmMode = 3){

  # Get indices
  rsp = pred
  ksp = prey

  # Get parameter values
  H_1 <- exp(params$logH_1)
  H_1a <- exp(params$logH_1a)
  H_1b <- exp(params$logH_1b)

  H_2 <- exp(params$logH_2)
  H_3 <- exp(params$logH_3)

  H_4 <- params$H_4

  # Set up ratios
  Pred_r <- seq(from = 0.001, to = 5, length.out = 100) # Pred biomass relative to equilibrium
  Prey_r <- seq(from = 0.001, to = 5, length.out = 100) # Prey biomass relative to equilibrium

  # Calculate functional form
  Term = H_1[rsp, ksp] * (1 + H_1a[rsp] * H_1b[rsp] / (pred_age + H_1b[rsp]))
  response <- matrix(NA, ncol = length(Prey_r), nrow = length(Pred_r))
  rownames(response) <- Pred_r
  colnames(response) <- Prey_r

  for(i in 1:length(Pred_r)){
    for(j in 1:length(Prey_r)){
      response[i, j] <- Prey_r[j] * switch (
        as.character(msmMode),
        "2" = { # Holling Type I (linear)
          Term},
        "3" = { #  Holling Type II
          Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] )
        },
        "4" = { #  Holling Type III
          Term * (1 + H_2[rsp, ksp]) * ((Prey_r[j] ) ^ H_4[rsp, ksp]) / (1 + H_2[rsp, ksp] * ((Prey_r[j] ) ^ H_4[rsp, ksp])  )},
        "5" = { #  predator interference
          Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] + H_3[rsp, ksp] * (Pred_r[i] - 1) )
        },
        "6" = { # predator preemption
          Term * (1 + H_2[rsp, ksp] ) / ( (1 + H_2[rsp, ksp] * Prey_r[j]) * (1 + H_3[rsp, ksp] * (Pred_r[i] - 1)) )
        },
        "7" = { # Hassell-Varley
          Term * (2 + H_2[rsp, ksp] ) / (1 + H_2[rsp, ksp] * Prey_r[j] + ((Prey_r[j] ) ^ H_4[rsp, ksp]))
        },
        "8" = { #  Ecosim
          Term / (1 + H_3[rsp, ksp] * (Pred_r[i] - 1 ))},
        {
          print("msmMode not implemented")
        }
      )
    }
  }

  response_reshape <- reshape2::melt(response, id.vars = c("Pred_r", "Prey_r"), measure.vars = "Response")
  colnames(response_reshape) <- c("Pred_r", "Prey_r", "Response")

  if(msmMode %in% c(2, 3, 4, 7)){
plot(x = response_reshape$Prey_r, y = response_reshape$Response, xlab = "Prey ratio", ylab = "Functional response", type = "l")
  }

  if(msmMode %in% c(8)){
    plot(x = response_reshape$Pred_r, y = response_reshape$Response, xlab = "Pred ratio", ylab = "Functional response", type = "l")
  }

  if(msmMode %in% c(5, 6)){
    filled.contour(response, xlab = "Prey ratio", ylab = "Predator ratio", col = oce.colorsDensity(30))
  }
}
