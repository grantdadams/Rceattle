#######################################
## MSMt Functions
## K. Holsman 2014
#######################################

ADMB2R<-function(ADMBfilename="tmsm",filepath="",use.headers=FALSE,datfiles="FALSE",ctlFile="TRUE",ctlFilename="tmsm.ctl"){
	# read in all .dat files used in the ADMB file named ADMBfilename
	# optional to set filepath different than working directory; include backslash at end e.g."~/Documents/Test/"
	#ADMBfilename is the name of the admb file, e.g., "tmsm"
	#filepath is the path to the folder that contains the ADMB file
	# use.headers if set to "TRUE" or "T" then it will read in data according to object names using the function readdat (must be loaded)

	fn<-file(paste(filepath,ADMBfilename,".tpl",sep=""))

	tpl.file<- scan(fn, what="character",flush=T,blank.lines.skip=F, quiet=T)
	skipp<-which(tpl.file=="DATA_SECTION")
	nrow<-which(tpl.file=="PARAMETER_SECTION")
	tpl.file<- scan(fn, what="character", skip=skipp,nlines=(nrow-skipp),flush=T,blank.lines.skip=F, quiet=T)
	nt<-length(tpl.file)
	tpl.tmp<-tpl.file
	for(i in 1:nt)
		if(tpl.file[i]!="")
			tpl.tmp[i]<- paste(scan(fn, skip=i+skipp-1,flush=F, sep="\t",nlines=1,quiet=TRUE, what="character",blank.lines.skip=TRUE),sep="",collapse=" ")
	tt<-strsplit(tpl.tmp,split=c(" ")) # find all the text lines
	tmp<-rep(0,nt) # data lines
	tmp2<-rep(0,nt) # change .dat file lines
	tmp3<-rep(0,nt) # read in data from *(ad_comm::global_datafile)
	tmp4<-rep(0,nt) # read in data from *(ad_comm::global_datafile)
	LocalCalcs<-which(tpl.file=="LOCAL_CALCS") 	# lines with Local calcs
	EndCalcs<-which(tpl.file=="END_CALCS") 		# lines with End calcs


	ai<-0  # list of file names
	# first find the data rows in the tpl file that come from .dat files
	ttCtl<-""
	ctl.file<-""
	if(ctlFile=="TRUE")
		ctl.file<- scan(paste(filepath,ctlFilename,sep=""), what="character",flush=T,blank.lines.skip=F, quiet=T)

	for (i in 1:nt){
		if(tpl.file[i]!=""){
			tt[[i]]<-tt[[i]][tt[[i]]!=""]
			if(tt[[i]][1]=="*(ad_comm::global_datafile)"){
				tmp3[i]<-1
			}else{
				tmp1<-strsplit(tt[[i]][1],split=c(">>"))
				if(tmp1[[1]][1]=="*(ad_comm::global_datafile)"){
					tmp4[i]<-1
						ttCtl<-c(ttCtl,strsplit(tmp1[[1]][2],split=";")[[1]][1])
				}
			}
			ttt<-strsplit(tt[[i]],split=c("_"))
			if(any(ttt[[1]]=="init"))
				tmp[i]<-1
			if(length(ttt)>1){
				if(any(ttt[[2]]=="comm::change")){
					t4<-strsplit(tt[[i]][2],split="\\(")
					t4<-strsplit(t4[[1]],split="\\)")
					tmp2[i]<-1
					ai<-ai+1
					if(ctlFile=="TRUE"){
						if(ai==1)
							tplA<-ctl.file[ai+1]
						if(ai>1)
							tplA[ai]<-ctl.file[ai+1]
					}else{
						if(ai==1)
							tplA<-paste(ADMBfilename,".dat",sep="")
						if(ai>1)
							tplA[ai]<-t4[[2]][1]
					}
				}
			}
		}
	}
	ttCtl<-ttCtl[-1]

	if(ctlFile=="TRUE")
		ctl.file2<-data.frame(tplFile=ttCtl,ctlFile=ctl.file)


	datfilesN<-rep(0,length(tplA));names(datfilesN)<-tplA
	datnmdim<-list()
	datanums<-which(tmp==1|tmp3==1)
	datanums2<-which(tmp3==1)  # which lines have change dat file
	datfilenums<-which(tmp2==1)
	nd<-length(datanums)
	datfileN<-1:length(datfilenums)
	dat<-list()
	dattype<-rep("",nd)
	notes<-rep("",nd)
	notes<-list()
	nmt<-rep(0,length(datfilenums))
	SVML<-rep(0,nd) # scalar(0), vector(1), matrix(2), list(3)

	## read in long data as a continuous string contained in the list ifile2
	ifile2<-list()
	for(i in 1:length(datfilenums)){
		ifile <- scan(paste(filepath,tplA[i],sep=""), what="character",flush=T,blank.lines.skip=F, quiet=T) # read in each line of the .dat file
		iflex<-which(is.na(ifile))
		idx <- sapply(as.double(ifile), is.na)
		idy<-which(idx) # which lines are text
		tmp1<--999
		for(zz in 2:length(idy)){
			ir<-idy[zz-1]
			irr<-idy[zz]
			if((irr-ir)>1){
				for(xx in 1:(irr-ir-1)){
					tmp2<-scan(paste(filepath,tplA[i],sep=""),what="",blank.lines.skip=F,skip=ir+xx-1,nlines=1, quiet=T)
					if(any(tmp2=="#"))
						tmp2<-tmp2[1:(which(tmp2=="#")-1)]  #remove values after an in-line #
					tmp2<-as.numeric(na.omit(as.double(tmp2)))
					tmp1<-c(tmp1,tmp2)
				}
			}
			if(zz==length(idy)){
				ir<-idy[zz]
				tmp2<-scan(paste(filepath,tplA[i],sep=""),what="",blank.lines.skip=F,skip=ir,nlines=1, quiet=T)
				if(any(tmp2=="#"))
					tmp2<-tmp2[1:(which(tmp2=="#")-1)]  #remove values after an in-line #
				tmp2<-as.numeric(na.omit(as.double(tmp2)))
				tmp1<-c(tmp1,tmp2)
			}
		}
		ifile2[[i]]<-tmp1[-1]
	}
	for(d in 1:nd){
		if(d>1){
			if(names(dat)[d-1]=="test_read")
				if(dat$test_read!=12345)
					stop("Problem with ",tplA[which(datfilenums<=datanums[d])][1],", test_read value is incorrect -->",dat$test_read)
			if(names(dat)[d-1]=="test_diet")
				if(dat$test_diet!=12345)
					stop("Problem with ",tplA[which(datfilenums<=datanums[d])][1],", test_diet value is incorrect -->",dat$test_diet)

			if(names(dat)[d-1]=="nages")
				max_age=max(dat$nages)
			if(names(dat)[d-1]=="nlengths")
				max_length=max(dat$nlengths)
		}

		nrow<-0;ncol<-0
		dn<-which(datfilenums<=datanums[d])
		dn<-dn[length(dn)]
		datfilesN[dn]<-d
		datfile<-tplA[dn]
		fileN<-datfileN[dn]
		nmt.use<-nmt[dn]
		fn<-paste(filepath,datfile,sep="")
		#name, type, dims, data, notes
		t1<-strsplit(tpl.tmp[datanums[d]]," ")
		t1<-t1[[1]][t1[[1]]!=""]

		if(tmp3[datanums[d]]==0){
			nm1<-t1[2] # data object name
			t2<-strsplit(nm1,"\\(")[[1]]
			nm<-strsplit(t2[1],";")[[1]]
			datnmdim[[nm]]<-nm1
		}
		if(tmp3[datanums[d]]==1){
			nmtmp1<-rep(0,nt)
			skippp<-0
			nm1<-t1[3] # find the name of the object
			t2<-strsplit(nm1,"\\(")[[1]]
			nm<-strsplit(t2[1],";")[[1]]
			# now find the dimensions in the preceeding text strings
			for (i in 1:nt){
				if(any(LocalCalcs==i))
				stop()

				if(skippp==0){
					if(tpl.file[i]!="//"){
						if(tpl.file[i]!=""){
							tttmp<-tt[[i]][tt[[i]]!=""]
							ntt<-length(tttmp)
							for(iii in 1:ntt){
								tnt<-strsplit(tttmp[iii],"\\(")[[1]]
								if(any(tnt==nm)){
									nmtmp1[i]<-1
									datnmdim[[nm]]<-tttmp[iii]
									t2<-strsplit(tttmp[iii],"\\(")[[1]]
									skippp<-1
								}
							}
						}
					}
				}
			}
		}
		dm<-"1"
		if(length(t2)>1){
			dm<-t2[2]
			dm<-strsplit(dm,"\\)")[[1]][1]
			dm<-strsplit(dm,",")[[1]]
		}
		dat[nm]<-""
		if(any(t1=="//"))
			notes[[nm]]<-paste(t1[(which(t1=="//")+1):length(t1)],collapse=" ",sep="")
		nmt.use<-nmt[dn] # last place to start from
		if(use.headers=="TRUE"|use.headers=="T"){
			lngdat<-readdat(fn,nm)
		}else{
			ir <-nmt.use+1
			lngdat<-ifile2[[fileN]][ir:length(ifile2[[fileN]])]
		}
		#now read in the corresponding .dat file
		if(length(dm[dm=="1"])==1){
			ncol<-eval(parse(text=dm[length(dm)]))
			if(length(lngdat)>1)
				SVML[d]<-1 # is a vector, else is a scalar
			dum<-lngdat[1:ncol]
			eval(parse(text=paste(nm,"<-dum",sep="")))
			dat[[nm]]<-dum
			nmtadd<-length(dum)
		}else{
			dm2<-seq(2,length(dm),2)
			ncol<-eval(parse(text=dm[length(dm)]))
			nrow<-eval(parse(text=dm[length(dm)-2]))
			ragged<-0  # 0 means it is not a ragged array
			if(any(ncol!=ncol[1]))
				ragged<-1
			if(length(dm2)<=2){
				#then it is a matrix
				if(ragged==0){
					# then this is a matrix
					SVML[d]<-2 # is a matrix
					nn<-nrow*ncol[1]
					if(length(nn)>1)
						stop(paste("problem with nn: ",nn))
					dat[[nm]]<-matrix(lngdat[1:nn],nrow,ncol,byrow=T)
					nmtadd<-sum(nrow*ncol[1])
				}else{
					SVML[d]<-3  # is a list
					dum<-list()
					# if they have the same number of rows
					if(length(nrow)==1)
						nrow<-nrow/length(ncol)
					RowCol<-rbind(nrow,ncol)
					nn<-c(0,cumsum((RowCol[1,]*RowCol[2,])))
					nmtadd<-nn[length(nn)]
					for(cl in 1:length(ncol)){
						j<-cl+1
						tmp<-lngdat[(nn[j-1]+1):(nn[j])]
						dum[[cl]]<-matrix(tmp,RowCol[1,cl],RowCol[2,cl],byrow=T)
					}
					dat[[nm]]<-dum
				}
			}
			if(length(dm2)==3){
				# 3 d array
				SVML[d]<-3  # is a list
				dum<-list()
				dms<-list();dmsn<-rep(0,3)
				for(iii in 1:3){
					dms[iii]<-eval(parse(text=dm[dm2][iii]))
					dmsn[iii]<-length(dms[iii])
				}
				if(ragged==0&all(dmsn==1)){
					RowCol<-rbind(rep(nrow,dms[[1]][1]),rep(ncol,dms[[1]][1]))
					nn<-c(0,cumsum((RowCol[1,]*RowCol[2,])))
					nmtadd<-nn[length(nn)]
					for(cl in 1:dms[[1]][1]){
						j<-cl+1
						tmp<-lngdat[(nn[j-1]+1):(nn[j])]
						dum[[cl]]<-matrix(tmp,RowCol[1,cl],RowCol[2,cl],byrow=T)
					}
				}else{
					RowCol<-rbind(nrow,ncol)
					nn<-c(0,cumsum((RowCol[1,]*RowCol[2,])))
					nmtadd<-nn[length(nn)]
					for(cl in 1:length(ncol)){
						j<-cl+1
						tmp<-lngdat[(nn[j-1]+1):(nn[j])]
						dum[[cl]]<-matrix(tmp,RowCol[1,cl],RowCol[2,cl],byrow=T)
					}
				}
				dat[[nm]]<-dum
			}
			if(length(dm2)>3){
				ndm<-length(dm2)
				# multi d array
				SVML[d]<-ndm  # is a list
				dum<-list()
				RowCol<-rbind(nrow,ncol)  # dims of the last matrix
				dms<-list();dmsn<-rep(0,ndm)
				dms.prod<-1
				for(dd in 1:ndm){
					dms[[dd]]<-eval(parse(text=dm[dm2][dd]))
					dmsn[dd]<-length(dms[[dd]])
				}
				for(dd in ndm:1)
					dms.prod<-dms.prod*sum(dms[[dd]])
				ragged<-0; warning(paste(nm, ": double check", ndm,"d array: code may not read this in correctly."))
				if(any(dmsn>1))
					ragged<-1
				if(tmp3[datanums[d]]==1){
					tpl1<-datanums[d]
					tmp1<-tpl.file[(tpl1-(ndm*4)):tpl1]
					tmp1n<-rep(0,length(tmp1))
					tpl2n<-tmp1n;tpl2n[length(tmp1)]<-tpl1
					tmptt<-1
					for(ii in length(tmp1):1){
						ttt<-strsplit(tmp1[ii],split="\\(")
						if(length(ttt[[1]])>1)
							if(ttt[[1]][1]=="for"){
								tmp1n[ii]<-1
								tpl2n[ii]<-tpl1-tmptt
								tpl2<-tpl1-tmptt		# this is where to start reading in the values
								tmptt<-1+tmptt
							}
					}
					# then read in the data
					fn2<-paste(filepath,ADMBfilename,".tpl",sep="")
					tmpt<-"";tmt<-""
					for(ii in (ndm):1){
						tmptt<-paste(scan(fn2, what="character", skip=skipp+tpl1-1-ii,nlines=1,flush=F,blank.lines.skip=F, quiet=T),collapse=" ")
						ttt<-strsplit(tmptt,split="=")[[1]]
						ttt<-strsplit(ttt,split=";")
						ob1<-strsplit(ttt[[3]][2],split="\\++")[[1]][1]
						tt1<-ttt[[2]][1]
						tt2<-strsplit(ttt[[3]][1],split="\\(")[[1]]
						tmt<-paste(tmt,"[[",ob1,"]]",sep="")
						if(ii!=1){
							if(length(tt2)>1){
								tt3<-strsplit(tt2[2],split="\\)")[[1]]
								tmpt<-paste(tmpt,"for(",ob1," in ",tt1,":",tt2[1],"[",tt3,"]","){;dum",tmt,"<-list();",sep="")
							}else{tmpt<-paste(tmpt,"for(",ob1," in ",tt1,":",tt2[1],"){;dum",tmt,"<-list();",sep="")}
						}else{
							y<-1
							tmptt<-paste(scan(fn2, what="character", skip=skipp+tpl1-1,nlines=1,flush=F,blank.lines.skip=F, quiet=T),collapse=" ")
							ttt<-strsplit(tmptt,split=">>")[[1]][2]
							ttt<-strsplit(ttt,split="\\(")[[1]][2]
							ttt<-strsplit(ttt,split="\\)")[[1]][1]
							ttt<-strsplit(ttt,split=",")[[1]]
							ttt<-ttt[1:(length(ttt)-1)]
							ttt<-paste("[[",ttt,"]]",sep="",collapse="")
							if(length(tt2)>1){
								tt3<-strsplit(tt2[2],split="\\)")[[1]]
								bb<-eval(parse(text=tt2[1]))
								if(length(dms[[ndm]])!=length(bb))
									dms[[ndm]]<-rep(dms[[ndm]],length(bb))
								tmpt<-paste(tmpt,"dum",ttt,"<-rep(0,dms[[ndm]][",tt3,"]);dum",ttt,"[",tt1,":",tt2[1],"[",tt3,"]]<-lngdat[y:(y-1+",tt2[1],"[",tt3,"])];y<-",tt2[1],"[",tt3,"]","+y",paste(rep("}",ndm-1),collapse="",sep=""),sep="")
							}else{
								tmpt<-paste(tmpt,"for(",ob1," in ",tt1,":",tt2[1],");",sep="")
								tmpt<-paste(tmpt,"dum",ttt,"<-rep(0,dms[[ndm]]);dum",ttt,"[1:",tt2[1],"]<-lngdat[y:(y-1+",tt2[1],")];y<-",tt2[1],"+y",paste(rep("}",ndm-1),collapse="",sep=""),sep="")
							}
						}
					}
					dum<-list()
					print(dms)
					eval(parse(text=tmpt))
					nmtadd<-y-1
					dat[[nm]]<-dum
					#read in from tpl.file
				}else{
					if(ragged==0){
						dms1<-dmsn
						for(cc in 1:length(dms1))
							dms1[cc]<-dms[[cc]]
						nn<-prod(dms1)
						dum<-array(lngdat[1:nn],dms1)  ### maybe not working yet
						nmtadd<-nn[length(nn)]
						dat[[nm]]<-dum
					}
					if(ragged==1)
						stop(paste(nm, ": can't read", ndm,"d ragged array directly into r: code may not read correctly from here out."))
				}
			}
		}
		nmt[dn]<-nmt[dn]+nmtadd
	}# end for each d
	print(paste("scanned data from ",paste(tplA,collapse=", ")))
	if(datfiles=="TRUE"){	return(datfilesN)}else{   return(dat)}
}
dat_names<-function(fn,nm){
	ifile <- scan(fn, what="character",flush=T,blank.lines.skip=T, quiet=T)
	iflex<-grep("#",ifile)
	return(ifile[iflex])
}
readdat<-function(fn,nm){
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
	nmr<-grep(nm,ifile)
	up<-nmr:length(ifile)
	skipp<-which(is.na(as.numeric(ifile[up]))==TRUE)
	keep<-which(is.na(as.numeric(ifile[up]))==FALSE)
	st.r<-up[keep[1]]
	stp.r<-up[skipp[skipp>keep[1]][1]]-1
	ifile[st.r:stp.r]
	rr<-st.r:stp.r
	nc<-length(as.numeric(scan(fn,what="",flush=F,blank.lines.skip=T,skip=rr[1]-1,nlines=1, quiet=T)))
	ans<-matrix(NA,length(rr),nc)
	for(r in 1:length(rr)){
		ans[r,]<-as.numeric(scan(fn,what="",flush=F,blank.lines.skip=F,skip=rr[r]-1,nlines=1, quiet=T,sep=""))
	}
	return(ans)
}
readdatOLD<-function(fn,nm){
	# fn is the file name
	#nm is the object name
	ifile <- scan(fn, what="character",flush=T,blank.lines.skip=T, quiet=T)
	iflex<-which(is.na(ifile))
	idx <- sapply(as.double(ifile), is.na)
	idy<-which(idx)
	datnum<-which(idx==FALSE)
	labnum<-which(idx==TRUE)
	vnam <- ifile[idx] #list names
	# remove vnam objects that are simply commented out numbers
	tmp<-rep(0,length(vnam))
	tt<-strsplit(vnam,split="#")
	for(i in 1:length(tmp))
		if(is.na(as.numeric(tt[[i]][2])))
			tmp[i]<-1
	vnam2<-vnam[tmp==1]
	tt<-strsplit(vnam2,split="#")
	tmp<-rep(0,length(vnam2))
	for(i in 1:length(tmp))
		if(length(tt[[i]])>1)
			tmp[i]<-1
	vnam2<-vnam2[tmp==1]
	labnum<-match(vnam2,ifile)
	ifilet<-strsplit(ifile,split="#")
	vnamt<-strsplit(vnam2,split="#")
	for(i in 1:length(vnam2))
		vnam2[i]<-vnamt[[i]][2]
	for(i in 1:length(ifile))
		ifile[i]<-ifilet[[i]][length(ifilet[[i]])]
	vnam2<-na.omit(vnam2)
	nv <- length(vnam2) #number of objects
	A <- list()
	ir <- 0
	vnam<-vnam2
	ii<-which(vnam==nm|vnam==paste(nm,":",sep="")|paste(vnam,";",sep="")==nm)
	if(length(ii)==0)
			stop (paste(nm," >> name of object in the ",ADMBfilename,".dat file does not match ",ADMBfilename,".tpl file",sep=""))
	ir <- match(vnam[ii], ifile) # find the matching name in the ifile set
	if (ii!=nv) {
		irr <- match(vnam[ii+1], ifile)
	} else {
		irr <- length(ifile)+1 #next row
	}
	ans<--999
	which(is.na(as.numeric(ifile[ir:irr]))==FALSE)
	irn<-ir+which(is.na(as.numeric(ifile[ir:irr]))==FALSE)-1
	for(i in 1:length(irn)){
		tt<-as.double(scan(fn, skip=irn[i]-1, nlines=1,quiet=TRUE, what=""))
		ans<-c(ans,as.numeric(na.omit(tt)))
	}
	ans<-ans[-1]

	return(ans)
}

#MSMdattoRlist(fn=file.path(path,"main/New_Tmsmv2.dat"))
MSMdattoRlist<-function(fn,nspp=3,skipp=2){
	ifile <- scan(fn, what="character", skip=skipp,flush=T,blank.lines.skip=FALSE, quiet=T)
	iflex<-which(is.na(ifile))
	idx <- sapply(as.double(ifile), is.na)
	datnum<-which(idx==FALSE)
	vnam <- ifile[idx] #list names
	# remove vnam objects that are simply commented out numbers
	tmp<-rep(0,length(vnam))
	tt<-strsplit(vnam,split="#")
	for(i in 1:length(tmp))
		if(is.na(as.numeric(tt[[i]][2])))
			tmp[i]<-1
	vnam2<-vnam[tmp==1]
	tt<-strsplit(vnam2,split="#")
	tmp<-rep(0,length(vnam2))
	for(i in 1:length(tmp))
		if(length(tt[[i]])>1)
			tmp[i]<-1
	vnam2<-vnam2[tmp==1]
	labnum<-match(vnam2,ifile)
	ifilet<-strsplit(ifile,split="#")
	vnamt<-strsplit(vnam2,split="#")


	for(i in 1:length(vnam2)){
		vnam2[i]<-strsplit(vnamt[[i]][2],split=":")[[1]][1]
	}
	for(i in 1:length(ifile))
		ifile[i]<-ifilet[[i]][length(ifilet[[i]])]
	vnam2<-na.omit(vnam2)
	nv <- length(vnam2) #number of objects
	A <- list()
	ir <- 0
	vnam<-vnam2

	for (i in 1:(nv)){
		ir <- match(vnam[i], ifile) # find the matching name in the ifile set
		if (i!=nv) {
			irr <- match(vnam[i+1], ifile)
		} else {
			irr <- length(ifile)+1 #next row
		}
		dum <- NA
		if (is.na(irr))
			dum <- 000
		if (irr-ir==2)
			dum <- as.double(scan(fn, skip=ir+skipp, nlines=1,quiet=TRUE, what=""))
		if (irr-ir>2){
			nall<-(irr-ir)
			tmpl<-rep(0,nall)
			nrow<-nall/nspp
			for(j in 1:length(tmpl))
				tmpl[j] <- length(as.double(scan(fn, skip=skipp+ir-1+j, nlines=1,quiet=TRUE, what="")))
			if(all(tmpl==tmpl[1])){
				dum<-as.matrix(read.table(fn, skip=ir,nrow=nall, fill=TRUE))
			}else{
				dum<-list()
				for(s in 1:nspp)
					dum[[s]]<-as.matrix(read.table(fn, skip=skipp+ir+(nrow*s-nrow),nrow=nrow, fill=TRUE))
			}
		}
		# Logical test to ensure dealing with numbers
		#if (is.numeric(dum))
			A[[vnam2[i]]] <- dum
	}
	return(A)
}
stdtoRlist<-function(tn="tmsm.std",k=3,nages=c(12,12,25)){
	# first assign names

	proj1<-(read.csv(tn,sep = "\t",header=FALSE,fill=TRUE))
	vnam<-strsplit(as.character(proj1[1,]),split=" ")
	vnam<-vnam[[1]]
	idx<-which(vnam=="")
	vnam<-vnam[-idx]
	vnam<-vnam[-length(vnam)]
	A<-(read.csv(tn,sep = "\t",header=TRUE,fill=TRUE))
	AA<-data.frame(matrix(0,length(A[,1]),length(vnam)))

	for(i in 1:length(A[,1]))
	{
		tmp<-strsplit(as.character(A[i,]),split=" ")
		tmp<-tmp[[1]]
		idx<-which(tmp=="")
		tmp<-tmp[-idx]
		tmp2<-as.numeric(tmp)
		idxs<-which(is.na(tmp2))
		txts<-tmp[idxs]
		AA[i,]<-tmp2
		AA[i,idxs]<-as.character(tmp[idxs])

	}
	names(AA)<-vnam
	return(AA)
}

reptoRlist <- function(fn) {
	ifile <- scan(fn, what="character", flush=TRUE,blank.lines.skip=FALSE, quiet=TRUE)
	idx <- sapply(as.double(ifile), is.na)
	vnam <- ifile[idx] #list names
	nv <- length(vnam) #number of objects
	A <- list()
	ir <- 0
	for (i in 1:nv) {
		ir <- match(vnam[i], ifile)
		if (i!=nv) {
			irr <- match(vnam[i+1], ifile)
		} else {
			irr <- length(ifile)+1 #next row
		}
		dum <- NA
		if (irr-ir==2) {
			dum <- as.double(scan(fn, skip=ir, nlines=1,quiet=TRUE, what=""))
		}
		if (irr-ir>2) {
			dum <- as.matrix(read.table(fn, skip=ir,nrow=irr-ir-1, fill=TRUE))
		}
		# Logical test to ensure dealing with numbers
		if (is.numeric(dum)) {
			A[[vnam[i]]] <- dum
		}
	}
	return(A)
}


cruletolist<-function(fn){
	ifile <- scan(fn, what="character", flush=TRUE,blank.lines.skip=FALSE, quiet=TRUE)
	idx <- sapply(as.double(ifile), is.na)
	idx_n<-which(idx)
	ifile_n <- as.double(scan(fn, what="numeric", sep=" ",blank.lines.skip=TRUE))

	ifile_c <- as.character(scan(fn, what="character", sep=" ",blank.lines.skip=TRUE))
	nrows<-length(ifile)
	ncols<-length(ifile_n)/nrows
	colnames<-scan(fn, what="character", sep=" ",blank.lines.skip=TRUE)[1:ncols]

	ifile_n2<-data.frame(matrix(data=ifile_n,ncol=ncols,nrow=nrows,byrow=TRUE))[-idx_n,]
	names(ifile_n2)<-colnames
	ifile_c2<-data.frame(matrix(data=ifile_c,ncol=ncols,nrow=nrows,byrow=TRUE))[-idx_n,]
	bn<-which(ifile_c2[1,]=="biomass")
	cn<-which(ifile_c2[1,]=="catch")
	blanks<-which(ifile_c2[1,]=="")
	ifile_n2[,bn]<-"biomass";ifile_n2[,cn]<-"catch"
	colnames[blanks]<-paste("blank",blanks,sep="")
	names(ifile_n2)<-colnames
	return(ifile_n2)
}


read_estrep2<-function(fn,val="M2_fut_1",printt=FALSE){
	ifile <- scan(fn, what = "character", flush = T, blank.lines.skip = T, quiet = T)
    tmp<-as.numeric(ifile)
    idx<-which(is.na(tmp))
    idy<-idx
   # idy <- which(idx2)
    datnum <- which(idx == FALSE)
    labnum <- which(idx == TRUE)
    vnam <- ifile[idx]
    special_vnam<-which(vnam=="ntemp_scen")
    tmp <- rep(0, length(vnam))
    tt <- strsplit(vnam, split = "#")
    tmp[(is.na(as.numeric(unlist(tt))))]<-1
    vnam2 <- vnam[tmp == 1]
    labnum <- match(vnam2, ifile)
    ifilet <- strsplit(ifile, split = "#")
    vnamt <- vnam2
    #for (i in 1:length(ifile))
    #	ifile[i] <- ifilet[[i]][length(ifilet[[i]])]
    ifile<-unlist(ifilet)
    vnam2 <- na.omit(vnam2)
    nv <- length(vnam2)
    A <- list()
    ir <- 0
    vnam <- vnam2
    vv<-which(vnam==val)
    all_dat<-list()
    for(ii in vv)
    {
    	ir<-which(ifile==vnam[ii])
    	if (ii != nv)
    	{
     	   irr <- which(ifile==vnam[ii + 1])
    	}else {
        	irr <- length(ifile) + 1
    	}
    	irn <- ir + which(is.na(as.numeric(ifile[ir:irr])) == FALSE) - 1
    	nr<-length(irn)
    	nc<-0
    	for (r in 1:nr){
    		tt <- as.double(scan(fn, skip = irn[r] - 1, nlines = 1, quiet = TRUE, what = ""))
    		nc<-max(nc,length(tt))
    	}
    	ans<-matrix(-999,nr,nc)
    	for (r in 1:nr)
    	{
    		 tt<- as.double(scan(fn, skip = irn[r] - 1, nlines = 1, quiet = TRUE, what = ""))
    		 if(length(tt)==length(ans[r,])){
    		 	ans[r,]<-tt
    		 }else{
    		 	ans[r,1:length(tt)]<-tt
    		 }

    	}

        eval(parse(text=paste("all_dat$",vnam[ii],"<-ans",sep="")))
        if(printt==TRUE)
        	print(paste(round(ii/nv,2)*100,"% complete :",vnam[ii]))
	}
    return(all_dat)
}
read_estrep<-function (fn,printt=TRUE)
{
    ifile <- scan(fn, what = "character", flush = T, blank.lines.skip = T, quiet = T)
    tmp<-as.numeric(ifile)
    idx<-which(is.na(tmp))

    idy<-idx
   # idy <- which(idx2)
    datnum <- which(idx == FALSE)
    labnum <- which(idx == TRUE)
    vnam <- ifile[idx]
    special_vnam<-which(vnam=="ntemp_scen")
    tmp <- rep(0, length(vnam))
    tt <- strsplit(vnam, split = "#")
    tmp[(is.na(as.numeric(unlist(tt))))]<-1
    vnam2 <- vnam[tmp == 1]
    labnum <- match(vnam2, ifile)
    ifilet <- strsplit(ifile, split = "#")
    vnamt <- vnam2
    #for (i in 1:length(ifile))
    #	ifile[i] <- ifilet[[i]][length(ifilet[[i]])]
    ifile<-unlist(ifilet)
    vnam2 <- na.omit(vnam2)
    nv <- length(vnam2)
    A <- list()
    ir <- 0
    vnam <- vnam2
    all_dat<-list()

    for(ii in 1:length(vnam))
    {
    	ir<-which(ifile==vnam[ii])
    	if (ii != nv)
    	{
     	   irr <- which(ifile==vnam[ii + 1])
    	}else {
        	irr <- length(ifile) + 1
    	}
    	irn <- ir + which(is.na(as.numeric(ifile[ir:irr])) == FALSE) - 1
    	nr<-length(irn)
    	nc<-0
    	for (r in 1:nr){
    		tt <- as.double(scan(fn, skip = irn[r] - 1, nlines = 1, quiet = TRUE, what = ""))
    		nc<-max(nc,length(tt))
    	}
    	ans<-matrix(-999,nr,nc)
    	for (r in 1:nr)
    	{
    		 tt<- as.double(scan(fn, skip = irn[r] - 1, nlines = 1, quiet = TRUE, what = ""))
    		 if(length(tt)==length(ans[r,])){
    		 	ans[r,]<-tt
    		 }else{
    		 	ans[r,1:length(tt)]<-tt
    		 }

    	}

        eval(parse(text=paste("all_dat$",vnam[ii],"<-ans",sep="")))
        if(printt==TRUE)
        	print(paste(round(ii/nv,2)*100,"% complete :",vnam[ii]))
	}

    return(all_dat)
}
read.rec<-function(fn,nspp=3){

	tmpt<-read.csv(fn,sep="")
	tmpr<-tmpt[tmpt[,2]=="rec_dev",3]
	ny<-nyrs<-length(tmpr)/3
	ln_mn_rec<-tmpt[tmpt[,2]=="ln_mn_rec",3]
	ln_mn_rec.se<-tmpt[tmpt[,2]=="ln_mn_rec",4]/sqrt(ny)
	rec_dev<-matrix(0,nspp,ny)
	rec_dev.se<-matrix(0,nspp,ny)
	for(s in 1:nspp){
		end<-ny*s
		st<-end-ny+1
		rec_dev.se[s,]<-tmpt[tmpt[,2]=="rec_dev",4][st:end]/sqrt(ny)
		rec_dev[s,]<-tmpt[tmpt[,2]=="rec_dev",3][st:end]
	}
	logR_obs<-matrix(0,nspp,nyrs)
	logR_obs.plus<-matrix(0,nspp,nyrs)
	logR_obs.minus<-matrix(0,nspp,nyrs)
	logR_obs.5.plus<-matrix(0,nspp,nyrs)
	logR_obs.5.minus<-matrix(0,nspp,nyrs)

	for (s in 1:nspp){
	 logR_obs[s,] = (ln_mn_rec[s] + rec_dev[s,] )
	 logR_obs.plus[s,] = (ln_mn_rec[s]+1.96*ln_mn_rec.se[s] + rec_dev[s,]+1.96*rec_dev.se[s,] )
	 logR_obs.minus[s,] = (ln_mn_rec[s]-1.96*ln_mn_rec.se[s] + rec_dev[s,]-1.96*rec_dev.se[s,] )
	 logR_obs.5.plus[s,] = (ln_mn_rec[s]+1*ln_mn_rec.se[s] + rec_dev[s,]+1*rec_dev.se[s,] )
	 logR_obs.5.minus[s,] = (ln_mn_rec[s]-1*ln_mn_rec.se[s] + rec_dev[s,]-1*rec_dev.se[s,] )
	}
	return(list(ln_mn_rec=ln_mn_rec,ln_mn_rec.se=ln_mn_rec.se,rec_dev=rec_dev,rec_dev.se=rec_dev.se,logR_obs=logR_obs,
		logR_obs.plus=logR_obs.plus,logR_obs.minus=logR_obs.minus,
		logR_obs.5.plus=logR_obs.5.plus,logR_obs.5.minus=logR_obs.5.minus))
}
read.Frate<-function(fn,nspp=3){

	tmpt<-read.csv(fn,sep="")
	tmpr<-tmpt[tmpt[,2]=="F_dev",3]
	ny<-nyrs<-length(tmpr)/3
	ln_mean_F<-tmpt[tmpt[,2]=="ln_mean_F",3]
	ln_mean_F.se<-tmpt[tmpt[,2]=="ln_mean_F",4]/sqrt(ny)
	F_dev<-matrix(0,nspp,ny)
	F_dev.se<-matrix(0,nspp,ny)
	for(s in 1:nspp){
		end<-ny*s
		st<-end-ny+1
		F_dev.se[s,]<-tmpt[tmpt[,2]=="F_dev",4][st:end]/sqrt(ny)
		F_dev[s,]<-tmpt[tmpt[,2]=="F_dev",3][st:end]
	}
	logF_obs<-matrix(0,nspp,nyrs)
	logF_obs.plus<-matrix(0,nspp,nyrs)
	logF_obs.minus<-matrix(0,nspp,nyrs)
	logF_obs.5.plus<-matrix(0,nspp,nyrs)
	logF_obs.5.minus<-matrix(0,nspp,nyrs)

	for (s in 1:nspp){
	 logF_obs[s,] = (ln_mean_F[s] + F_dev[s,] )
	 logF_obs.plus[s,] = (ln_mean_F[s]+1.96*ln_mean_F.se[s] + F_dev[s,]+1.96*F_dev.se[s,] )
	 logF_obs.minus[s,] = (ln_mean_F[s]-1.96*ln_mean_F.se[s] + F_dev[s,]-1.96*F_dev.se[s,] )
	 logF_obs.5.plus[s,] = (ln_mean_F[s]+1*ln_mean_F.se[s] + F_dev[s,]+1*F_dev.se[s,] )
	 logF_obs.5.minus[s,] = (ln_mean_F[s]-1*ln_mean_F.se[s] + F_dev[s,]-1*F_dev.se[s,] )
	}
	return(list(ln_mean_F=ln_mean_F,ln_mean_F.se=ln_mean_F.se,F_dev=F_dev,F_dev.se=F_dev.se,logF_obs=logF_obs,
		logF_obs.plus=logF_obs.plus,logF_obs.minus=logF_obs.minus,
		logF_obs.5.plus=logF_obs.5.plus,logF_obs.5.minus=logF_obs.5.minus))
}
