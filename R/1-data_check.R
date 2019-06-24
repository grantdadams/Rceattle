data_check <- function(data_list) {

    if (length(data_list$BTempC) != length(data_list$Tyrs)) {
        stop("Length of temperature series is not the same as the number of temperature years")
    }

    if (length(data_list$M1_base) == 1) {
        stop("M1 is a single value, please make it age/species specific")
    }

    if (sum(data_list$other_food < 0) > 0) {
        stop("Other food for one species is negative")
    }
}
