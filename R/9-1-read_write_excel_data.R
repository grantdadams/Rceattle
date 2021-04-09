

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
    control <- matrix(NA, ncol = data_list$nspp, nrow = 19)
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
    control[13, ] <- data_list$pop_age_transition_index
    control[14, ] <- data_list$sigma_rec_prior
    control[15, ] <- data_list$other_food
    control[16, ] <- data_list$estDynamics
    control[17, ] <- data_list$proj_F
    control[18, ] <- data_list$est_sex_ratio
    control[19, ] <- data_list$sex_ratio_sigma
    control <- as.data.frame(control)
    control <- cbind(c("nspp", "styr", "endyr", "projyr", "nsex", "spawn_month", "R_sexr", "nages", "minage", "nlengths", "pop_wt_index", "ssb_wt_index","pop_age_transition_index", "sigma_rec_prior",
                       "other_food", "estDynamics", "proj_F", "est_sex_ratio", "sex_ratio_sigma"), control)
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
    xcel_list$age_error <- as.data.frame(data_list$age_error)
    names_used <- c(names_used, "age_error")


    # wt
    xcel_list$wt <- data_list$wt

    names_used <- c(names_used, "wt")


    # pmature
    colnames(data_list$pmature) <- c("Species", paste0("Age", 1:max(data_list$nages)))
    xcel_list$pmature <- as.data.frame(data_list$pmature)
    names_used <- c(names_used, "pmature")


    # sex_ratio
    colnames(data_list$sex_ratio) <- c("Species", paste0("Age", 1:max(data_list$nages)))
    xcel_list$sex_ratio <- as.data.frame(data_list$sex_ratio)
    names_used <- c(names_used, "sex_ratio")


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
    aLW <- as.data.frame(data_list$aLW)
    xcel_list$aLW <- aLW

    names_used <- c(names_used, "aLW")


    # bioenergetics_control
    bioenergetics_control <- matrix(NA, ncol = data_list$nspp, nrow = 12)
    bioenergetics_control[1, ] <- data_list$Ceq
    bioenergetics_control[2, ] <- data_list$Cindex
    bioenergetics_control[3, ] <- data_list$Pvalue
    bioenergetics_control[4, ] <- data_list$fday
    bioenergetics_control[5, ] <- data_list$CA
    bioenergetics_control[6, ] <- data_list$CB
    bioenergetics_control[7, ] <- data_list$Qc
    bioenergetics_control[8, ] <- data_list$Tco
    bioenergetics_control[9, ] <- data_list$Tcm
    bioenergetics_control[10, ] <- data_list$Tcl
    bioenergetics_control[11, ] <- data_list$CK1
    bioenergetics_control[12, ] <- data_list$CK4

    bioenergetics_control <- as.data.frame(bioenergetics_control)

    bioenergetics_control <- cbind(c("Ceq", "Cindex","Pvalue", "fday", "CA", "CB", "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4"), bioenergetics_control)
    colnames(bioenergetics_control) <- c("Object", data_list$spnames)
    names_used <- c(names_used, as.character(bioenergetics_control$Object))

    bioenergetics_control <- as.data.frame(bioenergetics_control)

    xcel_list$bioenergetics_control <- bioenergetics_control


    data_names[data_names %!in% names_used]


    # Temperature stuff
    xcel_list$env_data <- data_list$env_data
    names_used <- c(names_used, c("env_data"))



    # Diet information Pyrs
    xcel_list$Pyrs <- as.data.frame(data_list$Pyrs)
    names_used <- c(names_used, "Pyrs")


    # UobsAge
    xcel_list$UobsAge <- as.data.frame(data_list$UobsAge)
    names_used <- c(names_used, "UobsAge")

    # UobsAge
    xcel_list$UobsWtAge <- as.data.frame(data_list$UobsWtAge)
    names_used <- c(names_used, "UobsWtAge")


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
    data_list$pop_age_transition_index <- as.numeric(sheet1[13, 2:(data_list$nspp + 1)])
    data_list$sigma_rec_prior <- as.numeric(sheet1[14, 2:(data_list$nspp + 1)])
    data_list$other_food <- as.numeric(sheet1[15, 2:(data_list$nspp + 1)])
    data_list$estDynamics <- as.numeric(sheet1[16, 2:(data_list$nspp + 1)])
    data_list$proj_F <- as.numeric(sheet1[17, 2:(data_list$nspp + 1)])
    data_list$est_sex_ratio <- as.numeric(sheet1[18, 2:(data_list$nspp + 1)])
    data_list$sex_ratio_sigma <- as.numeric(sheet1[19, 2:(data_list$nspp + 1)])

    nyrs <- data_list$endyr - data_list$styr + 1


    # srv and fsh bits
    srv_bits <- c("fleet_control", "srv_biom", "fsh_biom" , "comp_data", "emp_sel", "NByageFixed")
    for (i in 1:length(srv_bits)) {
        sheet <- as.data.frame(readxl::read_xlsx(file, sheet = srv_bits[i]))
        sheet <- sheet[rowSums(is.na(sheet)) != ncol(sheet), ]

        data_list[[srv_bits[i]]] <- sheet
    }

    data_list$fleet_control$Nselages <- suppressWarnings(as.numeric(data_list$fleet_control$Nselages))

    # age_trans_matrix
    age_trans_matrix <- as.data.frame(readxl::read_xlsx(file, sheet = "age_trans_matrix"))
    data_list$age_trans_matrix <- age_trans_matrix



    # age_error
    age_error <- as.data.frame(readxl::read_xlsx(file, sheet = "age_error"))
    age_error <- age_error[rowSums(is.na(age_error)) != ncol(age_error), ] # Remove rows with all NA's
    data_list$age_error <- age_error


    # wt
    wt <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "wt")))
    data_list$wt <- wt


    # pmature
    data_list$pmature <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "pmature")))


    # sex_ratio
    data_list$sex_ratio <- suppressMessages(as.data.frame(readxl::read_xlsx(file, sheet = "sex_ratio")))


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
    env_data <- as.data.frame(readxl::read_xlsx(file, sheet = "env_data"))
    data_list$env_data <- env_data




    # Diet information Pyrs
    pyrs_matrix <- as.data.frame(readxl::read_xlsx(file, sheet = "Pyrs"))
    data_list$Pyrs <- pyrs_matrix

    # Diet UobsAge
    UobsAge <- as.data.frame(readxl::read_xlsx(file, sheet = "UobsAge"))
    data_list$UobsAge <- UobsAge

    # Diet UobsWtAge
    UobsWtAge <- as.data.frame(readxl::read_xlsx(file, sheet = "UobsWtAge"))
    data_list$UobsWtAge <- UobsWtAge


    # write the data
    return(data_list)
}
