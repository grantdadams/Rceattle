

#' Write data file
#'
#' @param data_list Rceattle data_list object
#' @param file Filname to be used. Must end with ".xlsx"
#'
#' @return
#' @export
#'
#' @examples
#'
#' library(Rceattle)
#' data(BS2017SS)
#' write_excel(data_list = BS2017SS, file = "BS2017SS.xlsx")
write_excel <- function( data_list, file = "Rceattle_data.xlsx" ){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  ## setup a workbook
  data_names <- names(data_list)
  ceattle_wb <- openxlsx::createWorkbook()
  names_used <- c()

  # control
  control <- matrix(NA, ncol = data_list$nspp, nrow = 9)
  control[1,1] <- data_list$nspp
  control[2,1] <- data_list$styr
  control[3,1] <- data_list$endyr
  control[4,] <- data_list$nages
  control[5,] <- data_list$nlengths
  control[6,] <- data_list$pop_wt_index
  control[7,] <- data_list$pop_alk_index
  control[8,] <- data_list$other_food
  control[9,] <- data_list$stom_tau
  control <- as.data.frame(control)
  colnames(control) <- paste0("Species_", 1:ncol(control))
  rownames(control) <- c("nspp", "styr", "endyr", "nages", "nlengths", "pop_wt_index", "pop_alk_index", "other_food", "stom_tau")
  names_used <- c(names_used, rownames(control))

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "control", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 1, x = control)


  # srv and fsh bits
  srv_bits <- c("srv_control", "srv_biom", "srv_emp_sel", "srv_comp", "fsh_control", "fsh_biom", "fsh_emp_sel", "fsh_comp")
  for(i in 1:length(srv_bits)){
    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = srv_bits[i], gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = i+1, x = data_list[[srv_bits[i]]])
  }
  names_used <- c(names_used, srv_bits)


  # 3D arrays by age


  # age_trans_matrix
  index_species <- data.frame(
    ALK = c(data_list$srv_control$ALK_index, data_list$fsh_control$ALK_index),
    Sp = c(data_list$srv_control$Species, data_list$fsh_control$Species))
  index_species <- index_species[!duplicated(index_species[,c('ALK','Sp')]),]
  index_species <- index_species[order(index_species$ALK),]
  nages <- data.frame(Sp = 1:data_list$nspp, Nages = data_list$nages, Nlengths = data_list$nlengths)
  index_species <- merge(index_species, nages, by = "Sp", all = TRUE )

  alk_dat <- matrix(NA, ncol = max(index_species$Nlengths) + 3, nrow = sum(index_species$Nages))
  ages_done <- 0
  for(i in 1:nrow(index_species)){
    for(age in 1:index_species$Nages[i]){
      ages_done <- ages_done + 1
      alk_dat[ages_done,1] <- index_species$ALK[i] # ALK index
      alk_dat[ages_done,2] <- index_species$Sp[i]  # Species index
      alk_dat[ages_done,3] <- age  # Age index
      alk_dat[ages_done,4:ncol(alk_dat)] <- data_list$age_trans_matrix[age,,index_species$ALK[i]]  # Age index
    }
  }
  colnames(alk_dat) <- c("ALK_index", "Species", "Age", paste0("Length_", 1:(ncol(alk_dat) - 3)))
  alk_dat <- as.data.frame(alk_dat)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "age_trans_matrix", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 10, x = alk_dat)

  names_used <- c(names_used, "age_trans_matrix")


  # age_error
  age_error <- matrix(NA, ncol = max(index_species$Nages) + 2, nrow = sum(index_species$Nages))
  ages_done <- 0
  for(sp in 1:data_list$nspp){
    for(age in 1:data_list$nages[sp]){
      ages_done <- ages_done + 1
      age_error[ages_done,1] <- sp # Species index
      age_error[ages_done,2] <- age # Species index
      age_error[ages_done,((1:data_list$nages[sp]) + 2)] <- data_list$age_error[sp, age, 1:data_list$nages[sp]] # Species index
    }
  }
  colnames(age_error) <- c("Species", "True_age", paste0("Obs_age", 1:max(data_list$nages)))
  age_error <- as.data.frame(age_error)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "age_error", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 11, x = age_error)

  names_used <- c(names_used, "age_error")


  # wt
  index_species <- data.frame(
    wt = c(data_list$srv_control$Weight_index, data_list$fsh_control$Weight_index),
    Sp = c(data_list$srv_control$Species, data_list$fsh_control$Species))
  index_species <- index_species[!duplicated(index_species[,c('wt','Sp')]),]
  index_species <- index_species[order(index_species$wt),]

  wt <- matrix(NA, ncol = max(data_list$nages) + 4, nrow = nrow(index_species) * length(data_list$styr:data_list$endyr) )
  yrs_done <- 1
  for(wt_index in 1:nrow(index_species)){
    nyears <-  length(data_list$styr:data_list$endyr)
    wt[ yrs_done : (yrs_done + nyears - 1), 1] <- paste0("Index_",wt_index)
    wt[ yrs_done : (yrs_done + nyears - 1), 2] <- wt_index
    wt[ yrs_done : (yrs_done + nyears - 1), 3] <- index_species$Sp[wt_index]
    wt[ yrs_done : (yrs_done + nyears - 1), 4] <- data_list$styr:data_list$endyr
    wt[ yrs_done : (yrs_done + nyears - 1), 5:(4 + data_list$nages[index_species$Sp[wt_index]])] <- data_list$wt[,1:(data_list$nages[index_species$Sp[wt_index]]),wt_index]
    yrs_done <- yrs_done + length(data_list$styr:data_list$endyr)
  }

  colnames(wt) <- c("Wt_name", "Wt_index", "Species", "Year", paste0("Age", 1:max(data_list$nages)))
  wt <- as.data.frame(wt)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "wt", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 12, x = wt)

  names_used <- c(names_used, "wt")


  # pmature
  pmature <- data_list$pmature
  colnames(pmature) <- paste0("Age", 1:max(data_list$nages))
  rownames(pmature) <- paste0("Species_", 1:data_list$nspp)
  pmature <- as.data.frame(pmature)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "pmature", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 13, x = pmature)

  names_used <- c(names_used, "pmature")


  # propF
  propF <- data_list$propF
  colnames(propF) <- paste0("Age", 1:max(data_list$nages))
  rownames(propF) <- paste0("Species_", 1:data_list$nspp)
  propF <- as.data.frame(propF)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "propF", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 14, x = propF)

  names_used <- c(names_used, "propF")


  # M1_base
  M1_base <- data_list$M1_base
  colnames(M1_base) <- paste0("Age", 1:max(data_list$nages))
  rownames(M1_base) <- paste0("Species_", 1:data_list$nspp)
  M1_base <- as.data.frame(M1_base)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "M1_base", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 15, x = M1_base)

  names_used <- c(names_used, "M1_base")


  # Mn_LatAge
  Mn_LatAge <- data_list$Mn_LatAge
  colnames(Mn_LatAge) <- paste0("Age", 1:max(data_list$nages))
  rownames(Mn_LatAge) <- paste0("Species_", 1:data_list$nspp)
  Mn_LatAge <- as.data.frame(Mn_LatAge)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "Mn_LatAge", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 16, x = Mn_LatAge)

  names_used <- c(names_used, "Mn_LatAge")

  # aLW
  aLW <- data_list$aLW
  rownames(aLW) <- c("alpha", "beta")
  colnames(aLW) <- paste0("Species_", 1:data_list$nspp)
  aLW <- as.data.frame(aLW)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "aLW", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 17, x = aLW)

  names_used <- c(names_used, "aLW")


  # bioenergetics_control
  bioenergetics_control <- matrix(NA, ncol = data_list$nspp, nrow = 11)
  bioenergetics_control[1,] <- data_list$Ceq
  bioenergetics_control[2,] <- data_list$Pvalue
  bioenergetics_control[3,] <- data_list$fday
  bioenergetics_control[4,] <- data_list$CA
  bioenergetics_control[5,] <- data_list$CB
  bioenergetics_control[6,] <- data_list$Qc
  bioenergetics_control[7,] <- data_list$Tco
  bioenergetics_control[8,] <- data_list$Tcm
  bioenergetics_control[9,] <- data_list$Tcl
  bioenergetics_control[10,] <- data_list$CK1
  bioenergetics_control[11,] <- data_list$CK4

  bioenergetics_control <- as.data.frame(bioenergetics_control)
  colnames(bioenergetics_control) <- paste0("Species_", 1:ncol(bioenergetics_control))
  rownames(bioenergetics_control) <- c("Ceq", "Pvalue", "fday", "CA", "CB", "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4")
  names_used <- c(names_used, rownames(bioenergetics_control))

  bioenergetics_control <- as.data.frame(bioenergetics_control)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "bioenergetics_control", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 18, x = bioenergetics_control)


  data_names[data_names %!in% names_used]


  # Temperature stuff
  Temp_data <- data.frame(Tyrs = data_list$Tyrs, BTempC = data_list$BTempC)
  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "Temp_data", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 19, x = Temp_data)
  names_used <- c(names_used, c("Tyrs", "BTempC"))



  # Diet information
  # Pyrs
  Pyrs <- matrix(NA, ncol = max(data_list$nages) + 2, nrow = data_list$nspp * length(data_list$styr:data_list$endyr) )
  yrs_done <- 1
  for(sp in 1:data_list$nspp){
    nyears <-  length(data_list$styr:data_list$endyr)
    Pyrs[ yrs_done : (yrs_done + nyears - 1), 1] <- sp
    Pyrs[ yrs_done : (yrs_done + nyears - 1), 2] <- data_list$styr:data_list$endyr
    Pyrs[ yrs_done : (yrs_done + nyears - 1), 3:(2 + data_list$nages[sp])] <- data_list$Pyrs[1:nyears,1:(data_list$nages[sp]), sp]
    yrs_done <- yrs_done + nyears
  }

  colnames(Pyrs) <- c("Species", "Year", paste0("Age", 1:max(data_list$nages)))
  Pyrs <- as.data.frame(Pyrs)

  openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "Pyrs", gridLines = FALSE)
  openxlsx::writeDataTable(wb = ceattle_wb, sheet = 20, x = Pyrs)

  names_used <- c(names_used, "Pyrs")


  # Diet UobsAge
  if(length(dim(data_list$UobsAge)) == 4){
    UobsAge <- matrix(NA, ncol = 5, nrow = length(data_list$UobsAge) )
    dims <- dim(data_list$UobsAge)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nages[pred]){
          for(prey_a in 1:data_list$nages[prey]){
            UobsAge[ind, 1] <- pred
            UobsAge[ind, 2] <- prey
            UobsAge[ind, 3] <- pred_a
            UobsAge[ind, 4] <- prey_a
            UobsAge[ind, 5] <- data_list$UobsAge[pred, prey, pred_a, prey_a]
            ind = ind + 1
          }
        }
      }
    }

    colnames(UobsAge) <- c("Pred", "Prey", "Pred_age", "Prey_age", "Stomach_proportion_by_number")
    UobsAge <- as.data.frame(UobsAge)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "UobsAge", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 21, x = UobsAge)

    names_used <- c(names_used, "UobsAge")

  }

  if(length(dim(data_list$UobsAge)) == 5){
    UobsAge <- matrix(NA, ncol = 6, nrow = length(data_list$UobsAge) )
    dims <- dim(data_list$UobsAge)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nages[pred]){
          for(prey_a in 1:data_list$nages[prey]){
            for(yr in 1:dims[5]){
              UobsAge[ind, 1] <- pred
              UobsAge[ind, 2] <- prey
              UobsAge[ind, 3] <- pred_a
              UobsAge[ind, 4] <- prey_a
              UobsAge[ind, 5] <- yr
              UobsAge[ind, 6] <- data_list$UobsAge[pred, prey, pred_a, prey_a, yr]
              ind = ind + 1
            }
          }
        }
      }
    }

    colnames(UobsAge) <- c("Pred", "Prey", "Pred_age", "Prey_age", "Year", "Stomach_proportion_by_number")
    UobsAge <- as.data.frame(UobsAge)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "UobsAge", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 21, x = UobsAge)

    names_used <- c(names_used, "UobsAge")
  }



  # Diet UobsWtAge
  if(length(dim(data_list$UobsWtAge)) == 4){
    UobsWtAge <- matrix(NA, ncol = 5, nrow = length(data_list$UobsWtAge) )
    dims <- dim(data_list$UobsWtAge)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nages[pred]){
          for(prey_a in 1:data_list$nages[prey]){
            UobsWtAge[ind, 1] <- pred
            UobsWtAge[ind, 2] <- prey
            UobsWtAge[ind, 3] <- pred_a
            UobsWtAge[ind, 4] <- prey_a
            UobsWtAge[ind, 5] <- data_list$UobsWtAge[pred, prey, pred_a, prey_a]
            ind = ind + 1
          }
        }
      }
    }

    colnames(UobsWtAge) <- c("Pred", "Prey", "Pred_age", "Prey_age", "Stomach_proportion_by_weight")
    UobsWtAge <- as.data.frame(UobsWtAge)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "UobsWtAge", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 22, x = UobsWtAge)

    names_used <- c(names_used, "UobsWtAge")

  }


  if(length(dim(data_list$UobsWtAge)) == 5){
    UobsWtAge <- matrix(NA, ncol = 6, nrow = length(data_list$UobsWtAge) )
    dims <- dim(data_list$UobsWtAge)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nages[pred]){
          for(prey_a in 1:data_list$nages[prey]){
            for(yr in 1:dims[5]){
              UobsWtAge[ind, 1] <- pred
              UobsWtAge[ind, 2] <- prey
              UobsWtAge[ind, 3] <- pred_a
              UobsWtAge[ind, 4] <- prey_a
              UobsWtAge[ind, 5] <- yr
              UobsWtAge[ind, 6] <- data_list$UobsWtAge[pred, prey, pred_a, prey_a, yr]
              ind = ind + 1
            }
          }
        }
      }
    }

    colnames(UobsWtAge) <- c("Pred", "Prey", "Pred_age", "Prey_age", "Year", "Stomach_proportion_by_weight")
    UobsWtAge <- as.data.frame(UobsWtAge)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "UobsWtAge", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 22, x = UobsWtAge)

    names_used <- c(names_used, "UobsWtAge")
  }



  # Diet Uobs
  if(length(dim(data_list$Uobs)) == 4){
    Uobs <- matrix(NA, ncol = 5, nrow = length(data_list$Uobs) )
    dims <- dim(data_list$Uobs)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nlengths[pred]){
          for(prey_a in 1:data_list$nlengths[prey]){
            Uobs[ind, 1] <- pred
            Uobs[ind, 2] <- prey
            Uobs[ind, 3] <- pred_a
            Uobs[ind, 4] <- prey_a
            Uobs[ind, 5] <- data_list$Uobs[pred, prey, pred_a, prey_a]
            ind = ind + 1
          }
        }
      }
    }

    colnames(Uobs) <- c("Pred", "Prey", "Pred_length", "Prey_length", "Stomach_proportion_by_number")
    Uobs <- as.data.frame(Uobs)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "Uobs", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 23, x = Uobs)

    names_used <- c(names_used, "Uobs")

  }

  if(length(dim(data_list$Uobs)) == 5){
    Uobs <- matrix(NA, ncol = 6, nrow = length(data_list$Uobs) )
    dims <- dim(data_list$Uobs)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nlengths[pred]){
          for(prey_a in 1:data_list$nlengths[prey]){
            for(yr in 1:dims[5]){
              Uobs[ind, 1] <- pred
              Uobs[ind, 2] <- prey
              Uobs[ind, 3] <- pred_a
              Uobs[ind, 4] <- prey_a
              Uobs[ind, 5] <- yr
              Uobs[ind, 6] <- data_list$Uobs[pred, prey, pred_a, prey_a, yr]
              ind = ind + 1
            }
          }
        }
      }
    }

    colnames(Uobs) <- c("Pred", "Prey", "Pred_length", "Prey_length", "Year", "Stomach_proportion_by_number")
    Uobs <- as.data.frame(Uobs)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "Uobs", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 23, x = Uobs)

    names_used <- c(names_used, "Uobs")
  }



  # Diet UobsWt
  if(length(dim(data_list$UobsWt)) == 4){
    UobsWt <- matrix(NA, ncol = 5, nrow = length(data_list$UobsWt) )
    dims <- dim(data_list$UobsWt)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nlengths[pred]){
          for(prey_a in 1:data_list$nlengths[prey]){
            UobsWt[ind, 1] <- pred
            UobsWt[ind, 2] <- prey
            UobsWt[ind, 3] <- pred_a
            UobsWt[ind, 4] <- prey_a
            UobsWt[ind, 5] <- data_list$UobsWt[pred, prey, pred_a, prey_a]
            ind = ind + 1
          }
        }
      }
    }

    colnames(UobsWt) <- c("Pred", "Prey", "Pred_length", "Prey_length", "Stomach_proportion_by_weight")
    UobsWt <- as.data.frame(UobsWt)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "UobsWt", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 24, x = UobsWt)

    names_used <- c(names_used, "UobsWt")

  }


  if(length(dim(data_list$UobsWt)) == 5){
    UobsWt <- matrix(NA, ncol = 6, nrow = length(data_list$UobsWt) )
    dims <- dim(data_list$UobsWt)

    ind <- 1
    for(pred in 1:dims[1]){
      for(prey in 1:dims[2]){
        for(pred_a in 1:data_list$nlengths[pred]){
          for(prey_a in 1:data_list$nlengths[prey]){
            for(yr in 1:dims[5]){
              UobsWt[ind, 1] <- pred
              UobsWt[ind, 2] <- prey
              UobsWt[ind, 3] <- pred_a
              UobsWt[ind, 4] <- prey_a
              UobsWt[ind, 5] <- yr
              UobsWt[ind, 6] <- data_list$UobsWt[pred, prey, pred_a, prey_a, yr]
              ind = ind + 1
            }
          }
        }
      }
    }

    colnames(UobsWt) <- c("Pred", "Prey", "Pred_length", "Prey_length", "Year", "Stomach_proportion_by_weight")
    UobsWt <- as.data.frame(UobsWt)

    openxlsx::addWorksheet(wb = ceattle_wb, sheetName = "UobsWt", gridLines = FALSE)
    openxlsx::writeDataTable(wb = ceattle_wb, sheet = 24, x = UobsWt)

    names_used <- c(names_used, "UobsWt")
  }




  data_names[data_names %!in% names_used]


  # write the data
  openxlsx::saveWorkbook(ceattle_wb, file, overwrite = TRUE)
}
