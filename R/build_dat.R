build_dat <- function(ctlFilename = "asmnt2017_0", TMBfilename = "CEATTLE_BSAI_v01", dat_dir = "data/dat files/"){


  # Step 1 -- Extract data names used in TMB
  fn<-file(paste("src/", TMBfilename,".cpp",sep=""))

  cpp_file <- readLines(fn)
  skipp <- grep('MODEL INPUTS', cpp_file) # Line of data files
  nrow <- grep('PARAMETER SECTION', cpp_file) # Last line of data files
  cpp_file<- scan(fn, what="character", skip=skipp,nlines=(nrow-skipp),flush=T,blank.lines.skip=F, quiet=T)
  nt<-length(cpp_file)
  cpp_tmp <- c()
  data_lines <- grep('DATA_', cpp_file)

  for(i in 1:length(data_lines)){
    cpp_tmp[i]<- paste(scan(fn, skip=data_lines[i]+skipp-1,flush=F, sep="\t",nlines=1,quiet=TRUE, what="character",blank.lines.skip=TRUE),sep="",collapse=" ")
  }
  tt <- strsplit(cpp_tmp,split=c(" ")) # find all the text lines

  dat_names <- c() # Character string of variables used in model
  for(i in 1:length(data_lines)){
    dat_line <- grep('DATA_', tt[i][[1]])
    dat_call <- paste(tt[i][[1]][ dat_line: (dat_line + 2)], collapse="")
    dat_names[i] <- sub("\\).*", "", sub(".*\\(", "", dat_call))
  }


  # Step 2 -- Find location of data in dat files
  fn <- file(paste(dat_dir, ctlFilename,".ctl",sep=""))
  ctl_file <- readLines(fn)
  skipp <- grep('START filenames', ctl_file) # Line of data files
  nrow <- grep('END filenames', ctl_file) # Last line of data files
  ctl_file <- ctl_file[c(skipp:nrow)]
  dat_files <- ctl_file[-grep('#', ctl_file)]
  dat_files <- dat_files[grep('.dat', dat_files)]
  dat_files <- strsplit(dat_files, "/")
  dat_files <- sapply(dat_files, function(x) tail(x, n=1))

  dat_loc <- data.frame(dat_name = dat_names, datfile = rep(NA, length(dat_names)))

  for(i in 1:length(dat_files)){
    fn <- paste(dat_dir, dat_files[i],sep="")
    if(file.exists(fn)){
      dat_tmp <- scan(fn, what="character",flush=T,blank.lines.skip=F, quiet=T)
      for(j in 1:length(dat_names)){

        if(length(grep(paste0("\\b", "#", dat_names[j], "\\b"), dat_tmp)) > 0){
          dat_loc$datfile[j] <- dat_files[i]
        }
      }
    }
    if(!file.exists(fn)){
      print(paste("File", dat_files[i], "is not in the directory"))
    }
  }

  dat_loc <- dat_loc[complete.cases(dat_loc),]

  dat_list <- list()

  for(i in 1:nrow(dat_loc)){
    dat_list[[i]] <- readdat(paste(dat_dir, dat_loc[i,2],sep=""), as.character(dat_loc[i, 1]), nspp = 3)
    names(dat_list)[i] = as.character(dat_loc[i, 1])
  }


  return(dat_list)
}

