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
  nmr<-grep(paste0("\\b", "#", nm, "\\b"),ifile)
  up<-nmr:length(ifile)
  skipp<-which(is.na(as.numeric(ifile[up]))==TRUE)
  keep<-which(is.na(as.numeric(ifile[up]))==FALSE)



  if(nm %in% c("obs_catch", "wt", "srv_age_obs", "age_trans_matrix")){
    st.r<-up[keep[1]]
    stp.r<-up[skipp[skipp>keep[1]][nspp]]-1
    ifile[st.r:stp.r]
    rr<-st.r:stp.r

    # Max number of columns
    nc_vec <- c()
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



  } else {

  }

  return(ans)
}
