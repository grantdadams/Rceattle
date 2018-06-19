readdat<-function(fn, nm , nspp){
  # fn is the file name
  #nm is the object name
  ifile <- scan(fn, what="character",flush=T,blank.lines.skip=T, quiet=T)
  iflex<-grep("#",ifile)
  #iflex<-which(is.na(ifile))
  idx <- sapply(as.double(ifile), is.na)
  #idy<-which(idx)
  idy<-grep("#",ifile)
  datnum<-which(idx==FALSE)
  labnum<-which(idx==TRUE)
  vnam <- ifile[idy] #list names

  dat_tmp <- scan(file=fn,what=character(),sep="\n",quiet=TRUE) # Get values from each line
  dat_tmp <- gsub(" ", "", dat_tmp) # Remove spacing
  dat_tmp <- strsplit(dat_tmp, ":") # Substring at colon
  dat_tmp <- sapply(dat_tmp, `[[`, 1) # Take first element.. usually name
  nmr<-grep(paste0("\\b", "#", nm, "\\b"),dat_tmp) # Where is the element

  up<-nmr:length(ifile)
  skipp<-which(is.na(as.numeric(ifile[up]))==TRUE)
  keep<-which(is.na(as.numeric(ifile[up]))==FALSE)

  # Max number of columns
  nc_vec <- c()

  if(nm %in% c("obs_catch", "wt", "srv_age_obs", "age_trans_matrix", "tcb_obs", "srv_age_sizes")){ # Get data from arrays

    dat_list <- list()

    for(sp in 1:nspp){
      st.r<-up[skipp[which(diff(skipp) != 1)[sp]]] + 1
      stp.r<-up[keep[which(diff(keep) != 1)[sp]]] # Find sp break point
      ifile[st.r:stp.r]
      rr<-st.r:stp.r

      for(r in 1:length(rr)){
        nc_temp <- length(as.numeric(scan(fn,what="",flush=F,blank.lines.skip=F,skip=rr[r]-1,nlines=1, quiet=T,sep="")))
        nc_vec <- c(nc_vec, nc_temp)
      }
      max_nc <- max(nc_vec)

      # Get values
      ans<-matrix(NA,length(rr),max_nc)
      for(r in 1:length(rr)){
        nc_temp <- length(as.numeric(scan(fn,what="",flush=F,blank.lines.skip=F,skip=rr[r]-1,nlines=1, quiet=T,sep="")))
        for(c in 1:nc_temp){
          ans[r,c] <- as.numeric(scan(fn,what="",flush=F,blank.lines.skip=F,skip=rr[r]-1,nlines=1, quiet=T,sep=""))[c]
        }
      }

      dat_list[[sp]] <- ans

    }

    ans <- dat_list

  } else { # Get data from matrices and vectors

    st.r<-up[keep[1]]
    stp.r<-up[skipp[skipp>keep[1]][1]]-1 # index the nspp skipp after the first keep
    ifile[st.r:stp.r]
    rr<-st.r:stp.r

    for(r in 1:length(rr)){
      nc_temp <- length(as.numeric(scan(fn,what="",flush=F,blank.lines.skip=F,skip=rr[r]-1,nlines=1, quiet=T,sep="")))
      nc_vec <- c(nc_vec, nc_temp)
    }
    max_nc <- max(nc_vec)

    # Get values
    ans<-matrix(NA,length(rr),max_nc)
    for(r in 1:length(rr)){
      nc_temp <- length(as.numeric(scan(fn,what="",flush=F,blank.lines.skip=F,skip=rr[r]-1,nlines=1, quiet=T,sep="")))
      for(c in 1:nc_temp){
        ans[r,c] <- as.numeric(scan(fn,what="",flush=F,blank.lines.skip=F,skip=rr[r]-1,nlines=1, quiet=T,sep=""))[c]
      }
    }
  }

  return(ans)
}
