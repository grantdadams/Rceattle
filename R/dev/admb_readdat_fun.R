# Kirstin Holsman @kirstin.holsman@noaa.gov
readdat <- function(fn, nm, nspp) {
    # fn is the file name nm is the object name
    ifile <- scan(fn, what = "character", flush = T, blank.lines.skip = T, quiet = T)
    idx <- sapply(as.double(ifile), is.na)
    idy <- grep("#", ifile)
    datnum <- which(idx == FALSE)  # Which rows are numeric
    labnum <- which(idx == TRUE)  # Which rows are labels
    vnam <- ifile[idy]  #list names
    
    dat_tmp <- scan(file = fn, what = character(), sep = "\n", quiet = TRUE)  # Get values from each line
    dat_tmp <- gsub(" ", "", dat_tmp)  # Remove spacing
    dat_tmp <- strsplit(dat_tmp, ":")  # Substring at colon
    dat_tmp <- sapply(dat_tmp, `[[`, 1)  # Take first element.. usually name
    nmr <- grep(paste0("\\b", "#", nm, "\\b"), dat_tmp)  # Where is the element
    
    up <- nmr:length(ifile)
    skipp <- which(is.na(as.numeric(ifile[up])) == TRUE)  # Which rows are labels
    keep <- which(is.na(as.numeric(ifile[up])) == FALSE)  # Which rows are numeric
    
    # Max number of columns
    nc_vec <- c()
    
    if (nm %in% c("Qse", "mnQ")) {
        # Get data from 4D arrays
        nmr <- grep(paste0("\\b", "#", nm, "\\b"), dat_tmp) + 1  # Where is the element
        
        up <- nmr:length(ifile)
        skipp <- which(is.na(as.numeric(ifile[up])) == TRUE)  # Which rows are labels
        keep <- which(is.na(as.numeric(ifile[up])) == FALSE)  # Which rows are numeric
        
        dat_list <- list()
        
        for (sp in 1:nspp) {
            dat_list[[sp]] <- list()
            # Get number of years for each species
            st.r.sp <- up[skipp[which(diff(skipp) == 1)[sp]]] + 1  # Find first sp break point
            stp.r.sp <- up[keep[which(diff(keep) >= 3)[sp]]]  # Find last sp break point
            if (is.na(stp.r.sp)) {
                stp.r.sp <- up[skipp[length(skipp)]] - 1
            }
            
            # Number of years
            up.sp <- st.r.sp:stp.r.sp
            skipp.sp <- which(is.na(as.numeric(ifile[up.sp])) == TRUE)  # Which rows are labels
            nyr.sp <- length(skipp.sp)
            
            # Lines to be used
            up.sp <- st.r.sp:length(ifile)
            skipp.sp <- which(is.na(as.numeric(ifile[up.sp])) == TRUE)  # Which rows are labels
            keep.sp <- which(is.na(as.numeric(ifile[up.sp])) == FALSE)  # Which rows are numeric
            
            for (yr in 1:nyr.sp) {
                st.r <- up.sp[skipp.sp[which(diff(skipp.sp) != 1)[yr]]] + 1
                stp.r <- up.sp[keep.sp[which(diff(keep.sp) != 1)[yr]]]  # Find sp break point
                rr <- st.r:stp.r
                
                for (r in 1:length(rr)) {
                  nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                    quiet = T, sep = "")))
                  nc_vec <- c(nc_vec, nc_temp)
                }
                max_nc <- max(nc_vec)
                
                # Get values
                ans <- matrix(NA, length(rr), max_nc)
                for (r in 1:length(rr)) {
                  nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                    quiet = T, sep = "")))
                  for (c in 1:nc_temp) {
                    ans[r, c] <- as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                      quiet = T, sep = ""))[c]
                  }
                }
                dat_list[[sp]][[yr]] <- ans
            }
        }
        
        ans <- dat_list
        
    } else if (nm %in% c("Uobs", "UobsWt", "UobsAge", "UobsWtAge")) {
        dat_list <- list()
        
        ind <- 1
        for (sp in 1:nspp) {
            dat_list[[sp]] <- list()
            
            for (prey in 1:nspp) {
                st.r <- up[skipp[which(diff(skipp) != 1)[ind]]] + 1
                stp.r <- up[keep[which(diff(keep) != 1)[ind]]]  # Find sp break point
                rr <- st.r:stp.r
                
                for (r in 1:length(rr)) {
                  nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                    quiet = T, sep = "")))
                  nc_vec <- c(nc_vec, nc_temp)
                }
                max_nc <- max(nc_vec)
                
                # Get values
                ans <- matrix(NA, length(rr), max_nc)
                for (r in 1:length(rr)) {
                  nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                    quiet = T, sep = "")))
                  for (c in 1:nc_temp) {
                    ans[r, c] <- as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                      quiet = T, sep = ""))[c]
                  }
                }
                dat_list[[sp]][[prey]] <- ans
                ind <- ind + 1
            }
        }
        
        ans <- dat_list
        
        
    } else if (nm %in% c("obs_catch", "wt", "srv_age_obs", "age_trans_matrix", "tcb_obs", "srv_age_sizes", "K", "KAge", "PAge", 
        "Pyrs", "L2")) {
        # Get data from 3D arrays
        
        dat_list <- list()
        
        for (sp in 1:nspp) {
            st.r <- up[skipp[which(diff(skipp) != 1)[sp]]] + 1
            stp.r <- up[keep[which(diff(keep) != 1)[sp]]]  # Find sp break point
            ifile[st.r:stp.r]
            rr <- st.r:stp.r
            
            for (r in 1:length(rr)) {
                nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                  quiet = T, sep = "")))
                nc_vec <- c(nc_vec, nc_temp)
            }
            max_nc <- max(nc_vec)
            
            # Get values
            ans <- matrix(NA, length(rr), max_nc)
            for (r in 1:length(rr)) {
                nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                  quiet = T, sep = "")))
                for (c in 1:nc_temp) {
                  ans[r, c] <- as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                    quiet = T, sep = ""))[c]
                }
            }
            
            dat_list[[sp]] <- ans
            
        }
        
        ans <- dat_list
        
    } else {
        # Get data from matrices and vectors
        
        st.r <- up[keep[1]]
        stp.r <- up[skipp[skipp > keep[1]][1]] - 1  # index the nspp skipp after the first keep
        ifile[st.r:stp.r]
        rr <- st.r:stp.r
        
        for (r in 1:length(rr)) {
            nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                quiet = T, sep = "")))
            nc_vec <- c(nc_vec, nc_temp)
        }
        max_nc <- max(nc_vec)
        
        # Get values
        ans <- matrix(NA, length(rr), max_nc)
        for (r in 1:length(rr)) {
            nc_temp <- length(as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, 
                quiet = T, sep = "")))
            for (c in 1:nc_temp) {
                ans[r, c] <- as.numeric(scan(fn, what = "", flush = F, blank.lines.skip = F, skip = rr[r] - 1, nlines = 1, quiet = T, 
                  sep = ""))[c]
            }
        }
    }
    
    return(ans)
}
