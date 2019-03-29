

write_excel <- function( data_list, dir = NULL, file_name = "Rceattle_data" ){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  ## setup a workbook
  data_names <- names(data_list)
  ceattle_wb <- openxlsx::createWorkbook()
  names_used <- c()

  # control
  control <- matrix(NA, ncol = data_list$nspp, nrow = 8)
  control[1,1] <- data_list$nspp
  control[2,1] <- data_list$styr
  control[3,1] <- data_list$endyr
  control[4,] <- data_list$nages
  control[5,] <- data_list$nlengths
  control[6,] <- data_list$pop_wt_index
  control[7,] <- data_list$pop_alk_index
  control[8,] <- data_list$other_food
  control <- as.data.frame(control)
  colnames(control) <- paste0("Species_", 1:ncol(control))
  rownames(control) <- c("nspp", "styr", "endyr", "nages", "nlengths", "pop_wt_index", "pop_alk_index", "other_food")
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
  data_names[data_names %!in% names_used]

  # age_trans_matrix
  index_species <- data.frame(
    ALK = c(data_list$srv_control$ALK_index, data_list$fsh_control$ALK_index),
    Sp = c(data_list$srv_control$Species, data_list$fsh_control$Species))
  index_species <- index_species[!duplicated(index_species[,c('ALK','Sp')]),]
  index_species <- index_species[order(index_species$Sp),]
  nages <- data.frame(Sp = 1:data_list$nspp, Nages = data_list$nages, Nlengths = data_list$nlengths)
  index_species <- merge(index_species, nages, by = "Sp", all = TRUE )

  alk_dat <- matrix(NA, ncol = max(index_species$nlengths) + 3, nrow = sum(index_species$Nages))
  for(i in 1:nrow(index_species)){


  }

  arrays_3d <- c("age_error", "wt")

  # write the data
  substr_dir <- substr(dir, nchar(dir), nchar(dir))
  if(substr_dir == "/"){
    file_name <- paste0(dir, file_name, ".xlxs")
  }
  if(substr_dir != "/"){
    file_name <- paste0(dir, "/", file_name, ".xlxs")
  }
  saveWorkbook(ceattle_wb, file_name, overwrite = TRUE)
}
