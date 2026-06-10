# Fast CAAL obsvec-construction check (no model fit).
suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))
source(file.path("tests", "testthat", "helpers.R"))

dat <- make_test_data(nyrs = 8, nages = 5, seed = 123, growth = "vonBertalanffy")
dl  <- Rceattle::rearrange_data(dat)

cat("caal_obs rows      :", nrow(dl$caal_obs), "\n")
cat("caal obs_ctl rows  :", sum(dl$obs_ctl$type == "caal"), "\n")
cat("length col present :", "length" %in% names(dl$obs_ctl), "\n")

r1     <- which(dl$caal_obsvec_idx >= 0)[1]
sp     <- dl$caal_ctl[r1, 2]
n_caal <- dl$nages[sp]
start  <- dl$caal_obsvec_idx[r1]
expected <- (as.numeric(dl$caal_obs[r1, seq_len(n_caal)]) + 0.00001) * dl$caal_n[r1, 1]
cat("counts match log   :", isTRUE(all.equal(dl$obsvec[start + seq_len(n_caal)], expected)), "\n")

caal_rows <- dl$obs_ctl[dl$obs_ctl$type == "caal", ]
cat("one last-bin/group :", sum(caal_rows$is_last_bin) == length(unique(caal_rows$group_id)), "\n")
cat("length populated   :", all(!is.na(caal_rows$length)), "\n")
cat("age_or_length set  :", all(!is.na(caal_rows$age_or_length)), "\n")
