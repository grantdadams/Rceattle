

#' Write data file
#'
#' @param data_list Rceattle data_list object
#' @param file Filname to be used. Must end with '.xlsx'
#'
#' @export
#'
#' @examples
#'
#' library(Rceattle)
#' data(BS2017SS)
#' write_data(data_list = BS2017SS, file = 'BS2017SS.xlsx')
write_data <- function(data_list, file = "Rceattle_data.xlsx") {
    "%!in%" <- function(x, y) !(x %in% y)
    ## setup a workbook
    data_names <- names(data_list)
    names_used <- c()

    xcel_list <- list()
    # Metadata
    meta_filename <- system.file("extdata", "meta_data_names.xlsx", package = "Rceattle")
    meta_data <- suppressMessages((suppressWarnings(as.data.frame(readxl::read_xlsx(meta_filename)))))

    xcel_list$meta_data <- meta_data


    # control
    control <- matrix(NA, ncol = data_list$nspp, nrow = 17)
    control[1, 1] <- data_list$nspp
    control[2, 1] <- data_list$styr
    control[3, 1] <- data_list$endyr
    control[4, 1] <- data_list$projyr
    control[5, ] <- data_list$nsex
    control[6, ] <- data_list$spawn_month
    control[7, ] <- data_list$R_sexr
    control[8, ] <- data_list$nages
    control[9, ] <- data_list$minage
    control[10, ] <- data_list$nlengths
    control[11, ] <- data_list$pop_wt_index
    control[12, ] <- data_list$ssb_wt_index
    control[13, ] <- data_list$pop_alk_index
    control[14, ] <- data_list$sigma_rec_prior
    control[15, ] <- data_list$other_food
    control[16, ] <- data_list$estDynamics
    control[17, ] <- data_list$proj_F
    control <- as.data.frame(control)
    control <- cbind(c("nspp", "styr", "endyr", "projyr", "nsex", "spawn_month", "R_sexr", "nages", "minage", "nlengths", "pop_wt_index", "ssb_wt_index","pop_alk_index", "sigma_rec_prior",
                       "other_food", "estDynamics", "proj_F"), control)
    colnames(control) <- c("Object", data_list$spnames)
    names_used <- c(names_used, as.character(control$Object))

    xcel_list$control <- control


    # srv and fsh bits
    srv_bits <- c("fleet_control", "srv_biom", "fsh_biom", "comp_data",  "emp_sel", "NByageFixed", "age_trans_matrix")
    for (i in 1:length(srv_bits)) {
        xcel_list[[srv_bits[i]]] <- data_list[[srv_bits[i]]]
    }
    names_used <- c(names_used, srv_bits)

    # age_error
    index_species <- data.frame(ALK = c(data_list$fleet_control$ALK_index), Sp = c(data_list$fleet_control$Species))
    index_species <- index_species[!duplicated(index_species[, c("ALK", "Sp")]), ]
    index_species <- index_species[order(index_species$ALK), ]
    nages <- data.frame(Sp = 1:data_list$nspp, Nages = data_list$nages, Nlengths = data_list$nlengths)
    index_species <- merge(index_species, nages, by = "Sp", all = TRUE)
    age_error <- matrix(NA, ncol = max(index_species$Nages) + 2, nrow = sum(index_species$Nages))
    ages_done <- 0
    for (sp in 1:data_list$nspp) {
        for (age in 1:data_list$nages[sp]) {
            ages_done <- ages_done + 1
            sp_age <- data_list$minage[sp] + age - 1
            age_error[ages_done, 1] <- sp  # Species index
            age_error[ages_done, 2] <- sp_age  # Species index
            age_error[ages_done, ((1:data_list$nages[sp]) + 2)] <- data_list$age_error[sp, age, 1:data_list$nages[sp]]  # Species index
        }
    }
    colnames(age_error) <- c("Species", "True_age", paste0("Obs_age", 1:max(data_list$nages)))
    age_error <- as.data.frame(age_error)

    xcel_list$age_error <- age_error

    names_used <- c(names_used, "age_error")


    # wt
    xcel_list$wt <- data_list$wt

    names_used <- c(names_used, "wt")


    # pmature
    pmature <- data_list$pmature
    colnames(pmature) <- paste0("Age", 1:max(data_list$nages))
    rownames(pmature) <- paste0("Species_", 1:data_list$nspp)
    pmature <- as.data.frame(pmature)

    xcel_list$pmature <- pmature

    names_used <- c(names_used, "pmature")


    # propF
    propF <- data_list$propF
    colnames(propF) <- paste0("Age", 1:max(data_list$nages))
    rownames(propF) <- paste0("Species_", 1:data_list$nspp)
    propF <- as.data.frame(propF)

    xcel_list$propF <- propF

    names_used <- c(names_used, "propF")


    # M1_base
    M1_base <- data_list$M1_base
    M1_base <- as.data.frame(M1_base)

    xcel_list$M1_base <- M1_base

    names_used <- c(names_used, "M1_base")


    # Mn_LatAge
    Mn_LatAge <- data_list$Mn_LatAge
    colnames(Mn_LatAge) <- c("Species", "Sex", paste0("Age", 1:max(data_list$nages)))
    Mn_LatAge <- as.data.frame(Mn_LatAge)

    xcel_list$Mn_LatAge <- Mn_LatAge

    names_used <- c(names_used, "Mn_LatAge")

    # aLW
    aLW <- data_list$aLW
    rownames(aLW) <- c("alpha", "beta")
    colnames(aLW) <- paste0("Species_", 1:data_list$nspp)
    aLW <- as.data.frame(aLW)

    xcel_list$aLW <- aLW

    names_used <- c(names_used, "aLW")


    # bioenergetics_control
    bioenergetics_control <- matrix(NA, ncol = data_list$nspp, nrow = 11)
    bioenergetics_control[1, ] <- data_list$Ceq
    bioenergetics_control[2, ] <- data_list$Pvalue
    bioenergetics_control[3, ] <- data_list$fday
    bioenergetics_control[4, ] <- data_list$CA
    bioenergetics_control[5, ] <- data_list$CB
    bioenergetics_control[6, ] <- data_list$Qc
    bioenergetics_control[7, ] <- data_list$Tco
    bioenergetics_control[8, ] <- data_list$Tcm
    bioenergetics_control[9, ] <- data_list$Tcl
    bioenergetics_control[10, ] <- data_list$CK1
    bioenergetics_control[11, ] <- data_list$CK4

    bioenergetics_control <- as.data.frame(bioenergetics_control)

    bioenergetics_control <- cbind(c("Ceq", "Pvalue", "fday", "CA", "CB", "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4"), bioenergetics_control)
    colnames(bioenergetics_control) <- c("Object", data_list$spnames)
    names_used <- c(names_used, as.character(bioenergetics_control$Object))

    bioenergetics_control <- as.data.frame(bioenergetics_control)

    xcel_list$bioenergetics_control <- bioenergetics_control


    data_names[data_names %!in% names_used]


    # Temperature stuff
    Temp_data <- data.frame(Tyrs = data_list$Tyrs, BTempC = data_list$BTempC)
    xcel_list$Temp_data <- Temp_data
    names_used <- c(names_used, c("Tyrs", "BTempC"))



    # Diet information Pyrs
    xcel_list$Pyrs <- as.data.frame(data_list$Pyrs)
    names_used <- c(names_used, "Pyrs")


    # UobsAge
    xcel_list$UobsAge <- as.data.frame(data_list$UobsAge)
    names_used <- c(names_used, "UobsAge")

    # UobsAge
    xcel_list$UobsWtAge <- as.data.frame(data_list$UobsWtAge)
    names_used <- c(names_used, "UobsWtAge")


    # Diet Uobs
    if (length(dim(data_list$Uobs)) == 4) {
        Uobs <- matrix(NA, ncol = 5, nrow = length(data_list$Uobs))
        dims <- dim(data_list$Uobs)

        ind <- 1
        for (pred in 1:dims[1]) {
            for (prey in 1:dims[2]) {
                for (pred_a in 1:data_list$nlengths[pred]) {
                    for (prey_a in 1:data_list$nlengths[prey]) {
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

        xcel_list$Uobs <- Uobs

        names_used <- c(names_used, "Uobs")

    }

    if (length(dim(data_list$Uobs)) == 5) {
        Uobs <- matrix(NA, ncol = 6, nrow = length(data_list$Uobs))
        dims <- dim(data_list$Uobs)

        ind <- 1
        for (pred in 1:dims[1]) {
            for (prey in 1:dims[2]) {
                for (pred_a in 1:data_list$nlengths[pred]) {
                    for (prey_a in 1:data_list$nlengths[prey]) {
                        for (yr in 1:dims[5]) {
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

        xcel_list$Uobs <- Uobs

        names_used <- c(names_used, "Uobs")
    }



    # Diet UobsWt
    if (length(dim(data_list$UobsWt)) == 4) {
        UobsWt <- matrix(NA, ncol = 5, nrow = length(data_list$UobsWt))
        dims <- dim(data_list$UobsWt)

        ind <- 1
        for (pred in 1:dims[1]) {
            for (prey in 1:dims[2]) {
                for (pred_a in 1:data_list$nlengths[pred]) {
                    for (prey_a in 1:data_list$nlengths[prey]) {
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

        xcel_list$UobsWt <- UobsWt

        names_used <- c(names_used, "UobsWt")

    }


    if (length(dim(data_list$UobsWt)) == 5) {
        UobsWt <- matrix(NA, ncol = 6, nrow = length(data_list$UobsWt))
        dims <- dim(data_list$UobsWt)

        ind <- 1
        for (pred in 1:dims[1]) {
            for (prey in 1:dims[2]) {
                for (pred_a in 1:data_list$nlengths[pred]) {
                    for (prey_a in 1:data_list$nlengths[prey]) {
                        for (yr in 1:dims[5]) {
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

        xcel_list$UobsWt <- UobsWt

        names_used <- c(names_used, "UobsWt")
    }




    data_names[data_names %!in% names_used]


    # write the data
    writexl::write_xlsx(xcel_list, file)
}







#' Read a CEATTLE excel data file
#'
#' @param file Filname to be used. Must end with '.xlsx'
#'
#' @export
#'
#' @examples
#'
#' library(Rceattle)
#' data(BS2017SS)
#' write_data(data_list = BS2017SS, file = 'BS2017SS.xlsx')
#' data_list <- read_data(file = 'BS2017SS.xlsx')
read_data <- function(file = "Rceattle_data.xlsx") {
    "%!in%" <- function(x, y) !(x %in% y)
    ## setup a list object

    data_list <- list()


    sheet1 <- readxl::read_xlsx(file, sheet = "control")
    sheet1 <- as.data.frame(sheet1)

    # control
    data_list$nspp <- as.numeric(sheet1[1, 2])
    data_list$styr <- as.numeric(sheet1[2, 2])
    data_list$endyr <- as.numeric(sheet1[3, 2])
    data_list$projyr <- as.numeric(sheet1[4, 2])
    data_list$spnames <- as.character(colnames(sheet1)[2:(data_list$nspp + 1)])
    data_list$nsex <- as.numeric(sheet1[5, 2:(data_list$nspp + 1)])
    data_list$spawn_month <- as.numeric(sheet1[6, 2:(data_list$nspp + 1)])
    data_list$R_sexr <- as.numeric(sheet1[7, 2:(data_list$nspp + 1)])
    data_list$nages <- as.numeric(as.character(sheet1[8, 2:(data_list$nspp + 1)]))
    data_list$minage <- as.numeric(as.character(sheet1[9, 2:(data_list$nspp + 1)]))
    data_list$nlengths <- as.numeric(as.character(sheet1[10, 2:(data_list$nspp + 1)]))
    data_list$pop_wt_index <- as.numeric(sheet1[11, 2:(data_list$nspp + 1)])
    data_list$ssb_wt_index <- as.numeric(sheet1[12, 2:(data_list$nspp + 1)])
    data_list$pop_alk_index <- as.numeric(sheet1[13, 2:(data_list$nspp + 1)])
    data_list$sigma_rec_prior <- as.numeric(sheet1[14, 2:(data_list$nspp + 1)])
    data_list$other_food <- as.numeric(sheet1[15, 2:(data_list$nspp + 1)])
    data_list$estDynamics <- as.numeric(sheet1[16, 2:(data_list$nspp + 1)])
    data_list$proj_F <- as.numeric(sheet1[17, 2:(data_list$nspp + 1)])

    nyrs <- data_list$endyr - data_list$styr + 1


    # srv and fsh bits
    srv_bits <- c("fleet_control", "srv_biom", "fsh_biom" , "comp_data", "emp_sel", "NByageFixed")
    for (i in 1:length(srv_bits)) {
        sheet <- as.data.frame(readxl::read_xlsx(file, sheet = srv_bits[i]))
        sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]

        data_list[[srv_bits[i]]] <- sheet
    }

    data_list$fleet_control$Nselages <- suppressWarnings(as.numeric(data_list$fleet_control$Nselages))
    data_list$fleet_control$Nselages <- suppressWarnings(as.numeric(data_list$fleet_control$Nselages))

    # age_trans_matrix
    age_trans_matrix <- as.data.frame(readxl::read_xlsx(file, sheet = "age_trans_matrix"))
    data_list$age_trans_matrix <- age_trans_matrix



    # age_error
    age_error <- as.data.frame(readxl::read_xlsx(file, sheet = "age_error"))
    arm <- array(0, dim = c(data_list$nspp, max(data_list$nages, na.rm = T), max(data_list$nages, na.rm = T)))


    for (i in 1:nrow(age_error)) {
        sp <- as.numeric(as.character(age_error$Species[i]))
        true_age <- as.numeric(as.character(age_error$True_age[i])) - data_list$minage[sp] + 1

        if (true_age > data_list$nages[sp]) {
            message(paste0("Error: number of true ages specified in age_error for species: ", sp))
            message(paste0("is greater than the number of ages specified in the control"))
            message(paste0("Please remove or change nages in control"))
            stop()
        }

        arm[sp, true_age, 1:data_list$nages[sp]] <- as.numeric(as.character(age_error[i, (1:data_list$nages[sp]) + 2]))

        # Normalize
        arm[sp, true_age, 1:data_list$nages[sp]] <- arm[sp, true_age, 1:data_list$nages[sp]] / sum(arm[sp, true_age, 1:data_list$nages[sp]], na.rm = TRUE)
    }
    data_list$age_error <- arm



    # wt
    wt <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "wt")))
    data_list$wt <- wt


    # pmature
    pmature <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "pmature")))
    data_list$pmature <- pmature


    # propF
    propF <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "propF")))
    data_list$propF <- propF


    # M1_base
    M1_base <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "M1_base")))
    data_list$M1_base <- M1_base


    # Mn_LatAge
    Mn_LatAge <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "Mn_LatAge")))
    data_list$Mn_LatAge <- Mn_LatAge


    # aLW
    aLW <- as.data.frame(readxl::read_xlsx(file, sheet = "aLW"))
    data_list$aLW <- aLW


    # bioenergetics_control
    bioenergetics_control <- as.data.frame(readxl::read_xlsx(file, sheet = "bioenergetics_control"))

    for (i in 1:nrow(bioenergetics_control)) {
        data_list[[bioenergetics_control$Object[i]]] <- suppressWarnings(as.numeric(as.character(bioenergetics_control[i, ((1:data_list$nspp) +
                                                                                                                               1)])))
    }



    # Temperature stuff
    Temp_data <- as.data.frame(readxl::read_xlsx(file, sheet = "Temp_data"))
    data_list$Tyrs <- Temp_data$Tyrs
    data_list$BTempC <- Temp_data$BTempC




    # Diet information Pyrs
    pyrs_matrix <- as.data.frame(readxl::read_xlsx(file, sheet = "Pyrs"))
    data_list$Pyrs <- pyrs_matrix

    # Diet UobsAge
    UobsAge <- as.data.frame(readxl::read_xlsx(file, sheet = "UobsAge"))
    data_list$UobsAge <- UobsAge

    # Diet UobsWtAge
    UobsWtAge <- as.data.frame(readxl::read_xlsx(file, sheet = "UobsWtAge"))
    data_list$UobsWtAge <- UobsWtAge




    # Diet Uobs
    Uobs_matrix <- as.data.frame(readxl::read_xlsx(file, sheet = "Uobs"))

    # without year
    if (ncol(Uobs_matrix) == 5) {
        Uobs <- array(0, dim = c(data_list$nspp, data_list$nspp, max(data_list$nlengths), max(data_list$nlengths)))
        dims <- dim(data_list$Uobs)

        for (i in 1:nrow(Uobs_matrix)) {
            pred <- as.numeric(as.character(Uobs_matrix$Pred[i]))
            prey <- as.numeric(as.character(Uobs_matrix$Prey[i]))
            pred_length <- as.numeric(as.character(Uobs_matrix$Pred_length[i]))
            prey_length <- as.numeric(as.character(Uobs_matrix$Prey_length[i]))
            stom <- as.numeric(as.character(Uobs_matrix$Stomach_proportion_by_number[i]))

            Uobs[pred, prey, pred_length, prey_length] <- stom
        }

        data_list$Uobs <- Uobs
    }

    # with year
    if (ncol(Uobs_matrix) == 6) {
        Uobs <- array(0, dim = c(data_list$nspp, data_list$nspp, max(data_list$nlengths), max(data_list$nlengths)))
        dims <- dim(data_list$Uobs)

        for (i in 1:nrow(Uobs_matrix)) {
            pred <- as.numeric(as.character(Uobs_matrix$Pred[i]))
            prey <- as.numeric(as.character(Uobs_matrix$Prey[i]))
            pred_length <- as.numeric(as.character(Uobs_matrix$Pred_length[i]))
            prey_length <- as.numeric(as.character(Uobs_matrix$Prey_length[i]))
            year <- as.numeric(as.character(Uobs_matrix$Year[i]))
            stom <- as.numeric(as.character(Uobs_matrix$Stomach_proportion_by_number[i]))

            Uobs[pred, prey, pred_length, prey_length, year] <- stom
        }

        data_list$Uobs <- Uobs
    }



    # Diet UobsWt
    UobsWt_matrix <- as.data.frame(readxl::read_xlsx(file, sheet = "UobsWt"))

    # without year
    if (ncol(UobsWt_matrix) == 5) {
        UobsWt <- array(0, dim = c(data_list$nspp, data_list$nspp, max(data_list$nlengths), max(data_list$nlengths)))
        dims <- dim(data_list$UobsWt)

        for (i in 1:nrow(UobsWt_matrix)) {
            pred <- as.numeric(as.character(UobsWt_matrix$Pred[i]))
            prey <- as.numeric(as.character(UobsWt_matrix$Prey[i]))
            pred_length <- as.numeric(as.character(UobsWt_matrix$Pred_length[i]))
            prey_length <- as.numeric(as.character(UobsWt_matrix$Prey_length[i]))
            stom <- as.numeric(as.character(UobsWt_matrix$Stomach_proportion_by_weight[i]))

            UobsWt[pred, prey, pred_length, prey_length] <- stom
        }

        data_list$UobsWt <- UobsWt
    }

    # with year
    if (ncol(UobsWt_matrix) == 6) {
        UobsWt <- array(0, dim = c(data_list$nspp, data_list$nspp, max(data_list$nlengths), max(data_list$nlengths)))
        dims <- dim(data_list$UobsWt)

        for (i in 1:nrow(UobsWt_matrix)) {
            pred <- as.numeric(as.character(UobsWt_matrix$Pred[i]))
            prey <- as.numeric(as.character(UobsWt_matrix$Prey[i]))
            pred_length <- as.numeric(as.character(UobsWt_matrix$Pred_length[i]))
            prey_length <- as.numeric(as.character(UobsWt_matrix$Prey_length[i]))
            year <- as.numeric(as.character(UobsWt_matrix$Year[i]))
            stom <- as.numeric(as.character(UobsWt_matrix$Stomach_proportion_by_weight[i]))

            UobsWt[pred, prey, pred_length, prey_length, year] <- stom
        }

        data_list$UobsWt <- UobsWt
    }


    # write the data
    return(data_list)
}
