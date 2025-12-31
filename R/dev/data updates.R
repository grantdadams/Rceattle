data("NorthernRockfish2022")
data_list <- NorthernRockfish2022
data_list$fleet_control <- data_list$fleet_control %>%
  dplyr::mutate(Age_max_selected_upper = NA) %>%
  dplyr::rename(N_sel_bins = Nselages,
                Sel_norm_bin1 = Age_max_selected,
                Sel_norm_bin2 = Age_max_selected_upper) %>%
  mutate(Month = 0)
data_list$alpha_wt_len <- 0.0001
data_list$beta_wt_len <- 3

data_list$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = 7))
colnames(data_list$caal_data) <- c("Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", "CAAL_1")

NorthernRockfish2022 <- data_list

usethis::use_data(NorthernRockfish2022, overwrite = TRUE)
