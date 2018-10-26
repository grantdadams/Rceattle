
############################################################################################
## R code to plot CEATTLE outputs
## Kirstin Holsman
## Feb 2015
##
############################################################################################

	# rm(list=ls())	# clear workspaces
	graphics.off()	# close figures
	# rm(list=ls()); setwd("/Users/kholsman/GitHub/CEATTLE/runs/assmnt_2017_0")
	#;source("/Users/kholsman/GitHub/CEATTLE/src/ceattle-master/Scripts/R_code/PLOT_CEATTLE_EST.R")
	collrange1<-colorRampPalette(colors()[c(310,24)]) # grays
	collrange1<-colorRampPalette(colors()[c(310,24)]) # grays
	col_msm<-col2<-colorRampPalette(colors()[c(71,73)])
	col_ss<-col0<-col1<-colorRampPalette(c("orange","red"))
	col2<-col_msm
    col0<-col1<-col_ss
    col_dat<-colorRampPalette(colors()[c(491,73)])

############################################################################################
## Set up
############################################################################################
	update.est<-1 		# update figures (0=no, 1=yes)
	update.proj<-1 		# update figures (0=no, 1=yes)

	# update.figs<-1 		# update figures (0=no, 1=yes)
	update.data<-1 		# update data (0=no, 1=yes)
	save.plotdata<-1    # data dump the results of this script 
	#________________________________________
	# Read in the .config file
	#________________________________________

		tt<-read.csv("tmpfile.txt",header=FALSE,colClasses="character")[[1]]
		nt<-length(grep("#",tt))
		tt_nm<-unlist(strsplit(tt[grep("#",tt)],split="#"))[seq(2,nt*2,2)]
		tt<-read.csv("tmpfile.txt",comment="#",header=FALSE,colClasses="character")[[1]]
		for(itt in 1:nt)
			eval(parse(text=paste(tt_nm[itt],"<-'",tt[itt],"'",sep="")))
		if(testread!=12345){message(paste("Error with reading tmpfile.txt  line 29 of MAKE_recruitment_files.R in ",getwd()))}
		model_path<-file.path(DIR_main,model_path)
		rm(tt_nm)
		rm(nt)
		rm(tt)
	#________________________________________
	#________________________________________

	Type1<-mode
	path<-file.path(DIR_main,output_path)
	fn.use<-paste("results/",filename,"_rs.dat",sep="")

	output_root<-strsplit(output_path,filename)[[1]]
	# filename_root<-strsplit(output_path,"0")[[1]]
	filename_root<-substring(output_path,1,nchar(output_path)-2)
	out_0<-paste(filename_root,"_0",sep="")
	out_2<-paste(filename_root,"_2",sep="")

	
	output_path<-file.path(DIR_main,output_path)
	plot_file<-plot_fileJPG<-file.path(output_path,"R_figures")
	data_file<-file.path(model_path,"Scripts/Rdata_CSV")
	
	# assmt_tmp<-read.csv(file.path("/Users/kkari/Dropbox/Manuscripts/tmsm_manuscript/r_code_tmsm/ASSMT_2012_dat.csv"),sep=",")
	# assmt<-data.frame(t(assmt_tmp[,2:5]))
	# colnames(assmt)<-assmt_tmp[,1]; rm(assmt_tmp)
	source(file.path(model_path,"Scripts/R_code/PLOT_CEATTLE_FUN.R"))	


############################################################################################
## Functions
############################################################################################
	length.na <-function(x){if (any(is.na(x)==FALSE)){length(x[is.na(x)==FALSE])}else{NA}}
	mean.na<-function(x){mean(x,na.rm=TRUE)}
	sd.na<-function(x){sd(x,na.rm=TRUE)}
	se.na<-function(x){length.na(x)/sqrt(sd(x,na.rm=TRUE))}
	sum.na<-function(x){sum(x,na.rm=TRUE)}

	max.na<-function(x){max(x,na.rm=TRUE)}

	min.na<-function(x){min(x,na.rm=TRUE)}

	first.na<-function(x)
	{
		  if (any(is.na(x)==FALSE))
		  {
		    x[is.na(x)==FALSE][1]
		   }else{NA}
	}
	plusSE.na<-function(x){mean.na(x)+1.95*se.na(x)}
	minusSE.na<-function(x){mean.na(x)-1.95*se.na(x)}
	readdat<-function(fn,nm){
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
	makeTransparent<-function(someColor, alpha=100)
	{
	  newColor<-col2rgb(someColor)
	  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
	    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
	}



############################################################################################
## LOAD DATA
############################################################################################

	if (file.exists(plot_file)){} else {dir.create(plot_file)}
	if (file.exists(plot_fileJPG)){} else {dir.create(plot_fileJPG)}
	#readdat(fn=file.path(DIR_main,"CEATTLE-master/Data/dat_input_files/diet.dat"),nm="mn_wt_stom")

	load(file.path(DIR_main,out_0,"results/CEATTLE_results.Rdata"));ceattle_0<-tmp
	load(file.path(DIR_main,out_2,"results/CEATTLE_results.Rdata"));ceattle_2<-tmp
		years<-ceattle_2$styr[1,1]+1:ceattle_2$nyrs[1,1] -1

	# load(file.path(output_path,"ceattle_2/results/CEATTLE_results.Rdata"));ceattle_2<-tmp

	nspp<-as.numeric(ceattle_0$nspp)
	nyrs<-as.numeric(ceattle_0$nyrs)
	styr<-1979
	ceattle_0<-assign.rec(target=ceattle_0,fn=file.path(DIR_main,out_0,"results/ceattle_est.std"))
	ceattle_2<-assign.rec(target=ceattle_2,fn=file.path(DIR_main,out_2,"results/ceattle_est.std"))
	if(file.exists(plot_file)){}else{dir.create(file.path(plot_file))}
	collrange1<-colorRampPalette(colors()[c(121,128)]) # blues
	collrange2_msm<-colorRampPalette(colors()[c(92,554)]) # orange
	collrange1<-colorRampPalette(colors()[c(121,173)]) # blues
	ylimm_all<-rbind(c(0,2e7),c(0,3e6),c(0,2e6))

	graphics.off()

############################################################################################
##  Estimation mode - plots
############################################################################################
	if(update.est==1)
		update.figs<-1 		# update figures (0=no, 1=yes)


	############################################################################################
	##  Recruitment
	############################################################################################
		#plot_rec(logg=FALSE,ylimm=c(.3e12,2e9,.8e9),coll=c(colors()[280],"black"),asmt=F,asmt_dat="")
		#plot_rec(logg=FALSE,ylimm=c(.3e12,2e9,.8e9),coll=c(col0(1),col2(1)),asmt=F,asmt_dat="")
		plot_rec(logg=FALSE,ylimm=c(1.6e11,2.1e9,.8e9),coll=c(col0(1),col2(1)),asmt=F,asmt_dat="")
		if(update.figs==1){
			quartz.save(file=file.path(plot_fileJPG,"Recruitment.pdf"),type="pdf",dpi=500)
		}
		plot_rec(logg=FALSE,ylimm=c(1.6e11,2.1e9,.8e9),coll=c(col0(1),col2(1)),asmt=F,asmt_dat="")
		if(update.figs==1){
			quartz.save(file=file.path(plot_fileJPG,"Recruitmentv2.pdf"),type="pdf",dpi=500)
		}


	############################################################################################
	##  Biomass
	############################################################################################
		
		# plot_srvB(coll=c(colors()[330],"black","black"),ltyy=c(2,1))
		plot_srvB(coll=c(ss=col_ss(1),msm=col_msm(1),dat=col_dat(2)[1]),ltyy=c(2,1),pchh=c(-1,-1,19),txt=c("SSM","MSM","data"))
		plot_srvB2<-function(ltyy=c(1,2,1),coll=collrange1(2),lwdd=c(2,2,1)){
			quartz(height=6,width=6)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
			layout(c(1,2,3))
			layout.show(3)
			xatt<-seq(1960,2050,10)
			ylimm<-c(1e7,2e6,1e6)
			spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")

			for(s in 1:nspp)
			{
				eval(parse(text=paste("dat1<-ceattle_0$srv_bio_",s,sep="")))
				eval(parse(text=paste("dat2<-ceattle_2$srv_bio_",s,sep="")))
				if(s==1){
					eval(parse(text=paste("dat1<-ceattle_0$srv_bio_",s,sep="")))
					eval(parse(text=paste("dat2<-ceattle_2$srv_bio_",s,sep="")))
				}
				eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_srv_biom_",s,sep="")))
				eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_srv_biom_",s,sep="")))

				eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
				eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

				eval(parse(text=paste("dat1.hat<-ceattle_0$srv_bio_hat_",s,sep="")))
				eval(parse(text=paste("dat2.hat<-ceattle_2$srv_bio_hat_",s,sep="")))

				plot(dat2.yr[1,],dat2[1,],type="p",axes=FALSE,col=coll[3],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
				yatt<-pretty(seq(0,ylimm[s],1e4),4)
				axis(2,at=yatt,lab=yatt/1e6,las=2)
				cc<-.1/1000
				for(i in 1:length(dat2.yr[1,]))
				{
					x1<-c((1-cc)*dat2.yr[1,i],(1+cc)*dat2.yr[1,i])
					y1<-dat2[1,i]+1.96*dat2.se[1,i]
					y2<-dat2[1,i]-1.96*dat2.se[1,i]
					lines(rep(dat2.yr[1,i],2),c(y1,y2),col=coll[3])
					lines(x1,rep(y1,2),col=coll[3])
					lines(x1,rep(y2,2),col=coll[3])
				}

				lines(dat2.yr[1,],dat2.hat[1,],col=coll[2],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
				lines(dat1.yr[1,],dat1.hat[1,],col=coll[1],lwd=lwdd[1],lty=ltyy[1],cex=1.2)
				
				if(s==2)
					mtext("Survey biomass (million t)",side=2,line=3,outer=FALSE,font=2)
				text(1981,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
				if(s<3)
					axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
				if(s==3)
					axis(1,at=xatt,font=1)
				if(s==2)
					legend("topright",c("data","single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,pch=c(19,0,0),box.lty=0,cex=1.2)
			}	
			mtext("Year",side=1,line=0,outer=TRUE,font=2)	
	    }
		if(update.figs==1){
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"survey_biomass.eps")))
			quartz.save(file=file.path(plot_fileJPG,"survey_biomass.pdf"),type="pdf",dpi=500)
		}
			#plot_allBiomass(coll=c(colors()[330],"black"),ltyy=c(2,1),ylimm=c(2.6e7,2.2e6,1e6))
		if(update.figs==1){
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"plot_allBiomass.eps")))
			# plot_allBiomass(coll=c(colors()[330],"black"),ltyy=c(2,1),ylimm=c(2.6e7,2.2e6,1e6))
			plot_allBiomass(coll=c(ss=col_ss(1),msm=col_msm(1),dat=col_dat(2)[1]),coll2=rep(col_dat(1),3),ylimm=c(2.6e7,2.3e6,1e6))
			quartz.save(file=file.path(plot_fileJPG,"plot_allBiomass.pdf"),type="pdf",dpi=500)
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"plot_allBiomass.jpg"),res=500,units="in",height=7,width=5);dev.off()
		}
		col1<-colorRampPalette(colors()[c(71,73)])
		col2<-colorRampPalette(c("red","orange"))
		
		plot_catchB(coll=c(colors()[330],"black"),ltyy=c(2,1),ylimm=c(2e6,5e5,4e4))
		if(update.figs==1){		
			plot_catchB(coll=c(colors()[330],"black"),ltyy=c(2,1),ylimm=c(2e6,5e5,4e4))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"catch_biomass.jpg"),res=500,units="in",height=6,width=4);dev.off()
			quartz.save(file=file.path(plot_fileJPG,"catch_biomass.pdf"),type="pdf",dpi=500)
		}

	############################################################################################
	##  Mortality & predation
	############################################################################################
		plot_M2(coll=c(colors()[330],"black"),ltyy=c(2,1))
		if(update.figs==1){
			quartz.save(file=file.path(plot_fileJPG,"Mortality.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"Mortality.eps")))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"Mortality.jpg"),res=500,units="in",height=6,width=4);dev.off()
		}

		
		# plot_M2_all(coll=c(colors()[330],"black"),ltyy=c(2,1))
		plot_M2_all(coll_SS=col_msm(2)[1],coll_MSM=col_dat,ltyy=c(2,1),ylimm=c(2,2,2))
		
		if(update.figs==1){
			quartz.save(file=file.path(plot_fileJPG,"Mortalityall.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"Mortality.eps")))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"Mortality.jpg"),res=500,units="in",height=6,width=4);dev.off()
		}

		plot_Mort_q(coll=c(colors()[330],"black"),ltyy=c(2,1))
		if(update.figs==1)
			quartz.save(file=file.path(plot_fileJPG,"Mortality_q.pdf"),type="pdf",dpi=500)


		plot_Mort(coll=c(colors()[330],"black"),ltyy=c(2,1))
		mtext("Age 1 total mortality (M1+M2)",side=2,line=1,outer=TRUE,font=2)	
		mtext("Year",side=1,line=0,outer=TRUE  ,font=2)

		if(update.figs==1){
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"Mortality2.eps")))
			quartz.save(file=file.path(plot_fileJPG,"Mortality2.pdf"),type="pdf",dpi=500)

		}

		plot_F(coll=c(colors()[330],"black"),ltyy=c(2,1),ylim=c(1,1,1))
		if(update.figs==1){
			quartz.save(file=file.path(plot_fileJPG,"plot_F.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"plot_F.eps")))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"plot_F.jpg"),res=500,units="in",height=4,width=6);dev.off()
		}

		
		calc_eaten<-function(dat1=ceattle_2,pred=2,prey=1){
			eval(parse(text=paste("ration2Age<-dat1$ration2_",pred,sep="")))
			eval(parse(text=paste("AvgN<-dat1$AvgN_",pred,sep="")))
			eval(parse(text=paste("nages_pred<-dat1$nages_",pred,sep="")))
			eval(parse(text=paste("nages_prey<-dat1$nages_",prey,sep="")))
			#suit_main_"    <<sp<<"_"<<pred_age<
			eval(parse(text=paste("overlap<-dat1$overlap_",pred,sep="")))
			eval(parse(text=paste("avail_food<-dat1$avail_food_",pred,sep="")))

			#overlap=matrix(1,nyrs,nspp)
			nyrs<-ceattle_2$nyrs
			Pred_demand<-matrix(0,nyrs,nages_pred)
			eaten<-matrix(0,nyrs,nages_prey)
			#suit_main<-rep(0,nages_prey); suit_main[1:3]=seq(1,.02,-.4); suit_main=suit_main/sum(suit_main)
			Etmp<-matrix(0,nyrs,nages_prey)
			Mtmp<-matrix(0,nyrs,nages_prey)
			#Mtmp += AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/(avail_food(pred,i,pred_age));
		    
			for( i in 1:nyrs){

				for(pred_age in 1:nages_pred){
					eval(parse(text=paste("suit_main<-dat1$suit_main_",pred,"_",pred_age,sep="")))
			        Pred_demand[i,pred_age]=AvgN[i,pred_age]*overlap[i,prey] * ration2Age[i,pred_age]*1000
			        eaten[i,]=(Pred_demand[i,pred_age]*suit_main[prey,1:nages_prey])
			        
			        Etmp[i,]=Etmp[i,]+(Pred_demand[i,pred_age]*suit_main[prey,1:nages_prey])
			        Mtmp[i,] = Mtmp[i,]+(Pred_demand[i,pred_age]*suit_main[prey,1:nages_prey])/(avail_food[i,pred_age]*1000)
			        #if(i==1&pred_age==1){
			        #	plot(rep(avail_food[i,pred_age],nages_prey)/1e9,eaten[i,]/1e9,xlim=c(0,1.1*max(avail_food/1e9)),ylim=c(0,.8*max(avail_food/1e9)))
			        #}else{
			        #	points(rep(avail_food[i,pred_age],nages_prey)/1e9,eaten[i,]/1e9)
			        #}	
			    }  # end pred sp_age loop
		    } 
		    return(list(Pred_demand=Pred_demand,Etmp=Etmp,Mtmp=Mtmp,avail_food=avail_food,eaten=eaten))
		}	

	############################################################################################
	##  Selectivity plots
	############################################################################################
		plot_sel(coll=rep("black",2))
		if(update.figs==1){
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"selectivity.eps")))
			quartz.save(file=file.path(plot_fileJPG,"selectivity.pdf"),type="pdf",dpi=500)
		}

		plot_srvAge(coll=c(colors()[330],"black"),ltyy=c(2,1))
		if(update.figs==1){
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"surv_age_comp.eps")))
			quartz.save(file=file.path(plot_fileJPG,"surv_age_comp.pdf"),type="pdf",dpi=500)
		}
		#plot_srvAge_byYr(coll=c(colors()[330],"orange","black"))
		for(sp in 1:3){
		  plot_srvAge_byYr2(s=sp)
		  if(update.figs==1)
		    quartz.save(file=file.path(plot_fileJPG,paste0("surv_age_comp_by_yr_",sp,".pdf")),type="pdf",dpi=500)
		  
		}
	

		#plot_srvAge_byYr2(coll=c(colors()[330],"orange","black"))
	############################################################################################
	##  Fishery/catch plots
	############################################################################################
		plot_fshAge(coll=c(colors()[330],"black"),ltyy=c(2,1))
		if(update.figs==1){
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"fsh_age_comp.eps")))
		 	quartz.save(file=file.path(plot_fileJPG,"fsh_age_comp.pdf"),type="pdf",dpi=500)
		}
	#	plot_fshAge_byYr(coll=c(colors()[330],"orange","black"))
	#	if(update.figs==1)
			quartz.save(file=file.path(plot_fileJPG,"fsh_age_comp_by_yr.pdf"),type="pdf",dpi=500)
		#plot_srvAge_byYr(coll=c(colors()[330],"orange","black"))
		for(sp in 1){
		  plot_fshAge_byYr2(s=sp)
		  if(update.figs==1)
		    quartz.save(file=file.path(plot_fileJPG,paste0("fsh_age_comp_by_yr_",sp,".pdf")),type="pdf",dpi=500)
		  
		}
	############################################################################################
	##  calculate amount eaten
	############################################################################################
		calc_eaten<-function(dat1=ceattle_2,pred=2,prey=1){
			eval(parse(text=paste("ration2Age<-dat1$ration2_",pred,sep="")))
			eval(parse(text=paste("AvgN<-dat1$AvgN_",pred,sep="")))
			eval(parse(text=paste("nages_pred<-dat1$nages_",pred,sep="")))
			eval(parse(text=paste("nages_prey<-dat1$nages_",prey,sep="")))
			#suit_main_"    <<sp<<"_"<<pred_age<
			eval(parse(text=paste("overlap<-dat1$overlap_",pred,sep="")))
			eval(parse(text=paste("avail_food<-dat1$avail_food_",pred,sep="")))

			#overlap=matrix(1,nyrs,nspp)
			nyrs<-dat1$nyrs
			Pred_demand<-matrix(0,nyrs,nages_pred)
			eaten<-matrix(0,nyrs,nages_prey)
			#suit_main<-rep(0,nages_prey); suit_main[1:3]=seq(1,.02,-.4); suit_main=suit_main/sum(suit_main)
			Etmp<-matrix(0,nyrs,nages_prey)
			Mtmp<-matrix(0,nyrs,nages_prey)
			#Mtmp += AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/(avail_food(pred,i,pred_age));
		    
			for( i in 1:nyrs){

				for(pred_age in 1:nages_pred){
					eval(parse(text=paste("suit_main<-dat1$suit_main_",pred,"_",pred_age,sep="")))
			        Pred_demand[i,pred_age]=AvgN[i,pred_age]*overlap[i,prey] * ration2Age[i,pred_age]*1000
			        eaten[i,]=(Pred_demand[i,pred_age]*suit_main[prey,1:nages_prey])
			        
			        Etmp[i,]=Etmp[i,]+(Pred_demand[i,pred_age]*suit_main[prey,1:nages_prey])
			        Mtmp[i,] = Mtmp[i,]+(Pred_demand[i,pred_age]*suit_main[prey,1:nages_prey])/(avail_food[i,pred_age]*1000)
			        #if(i==1&pred_age==1){
			        #	plot(rep(avail_food[i,pred_age],nages_prey)/1e9,eaten[i,]/1e9,xlim=c(0,1.1*max(avail_food/1e9)),ylim=c(0,.8*max(avail_food/1e9)))
			        #}else{
			        #	points(rep(avail_food[i,pred_age],nages_prey)/1e9,eaten[i,]/1e9)
			        #}	
			    }  # end pred sp_age loop
		    } 
		    return(list(Pred_demand=Pred_demand,Etmp=Etmp,Mtmp=Mtmp,avail_food=avail_food,eaten=eaten))
		}
		pred<-2
		prey<-3
		dat1<-ceattle_2
		nages<-unlist(dat1[grep("nages_",names(dat1))])
		#suit_1<-dat1$suit_main_1_1[prey,]
		for(pred_age in 1:nages[pred])
		{
			if(pred_age==1)
			eval(parse(text=paste("suit_1<-rbind(dat1$suit_main_",pred,"_",pred_age,"[",prey,",])",sep="")))
			else
			eval(parse(text=paste("suit_1<-rbind(suit_1,dat1$suit_main_",pred,"_",pred_age,"[",prey,",])",sep="")))
		}
		tmp1_1<-calc_eaten(dat1=ceattle_2,pred=1,prey=1)
		tmp2_1<-calc_eaten(dat1=ceattle_2,pred=2,prey=1)
		tmp3_1<-calc_eaten(dat1=ceattle_2,pred=3,prey=1)

		tmp1_2<-calc_eaten(dat1=ceattle_2,pred=1,prey=2)
		tmp2_2<-calc_eaten(dat1=ceattle_2,pred=2,prey=2)
		tmp3_2<-calc_eaten(dat1=ceattle_2,pred=3,prey=2)

		tmp1_3<-calc_eaten(dat1=ceattle_2,pred=1,prey=3)
		tmp2_3<-calc_eaten(dat1=ceattle_2,pred=2,prey=3)
		tmp3_3<-calc_eaten(dat1=ceattle_2,pred=3,prey=3)


		# plot(tmp2_3$eaten,apply(tmp2_3$avail_food,1,sum))

		alldemand_1<-apply(tmp1_1$Pred_demand,1,sum)+apply(tmp2_1$Pred_demand,1,sum)+apply(tmp3_1$Pred_demand,1,sum)
		alldemand_2<-apply(tmp1_2$Pred_demand,1,sum)+apply(tmp2_2$Pred_demand,1,sum)+apply(tmp3_2$Pred_demand,1,sum)
		alldemand_3<-apply(tmp1_3$Pred_demand,1,sum)+apply(tmp2_3$Pred_demand,1,sum)+apply(tmp3_3$Pred_demand,1,sum)

		plot(alldemand_1/1e9,ylim=c(0,26),type="l")
		par(mfrow=c(3,1))

		prey1_eaten<-apply(tmp1_1$Etmp,1,sum)+apply(tmp1_1$Etmp,1,sum)+apply(tmp1_1$Etmp,1,sum)
		prey1_M2<-apply(tmp1_1$Mtmp,1,sum)+apply(tmp2_1$Mtmp,1,sum)+apply(tmp3_1$Mtmp,1,sum)

		#plot(prey1_eaten/1e6,ylim=c(0,26),type="l")
		plot(ceattle_2$PLK_consumed[1,]/1e6,ylim=c(0,26),type="l")
		lines(ceattle_2$Biomass_1[1,]/1e6,type="l",lty=2)
		lines(apply(tmp1_1$Etmp,1,sum)/1e6,type="l",lty=1,col="red")
		lines(apply(tmp2_1$Etmp,1,sum)/1e6,type="l",lty=2,col="red")
		lines(apply(tmp3_1$Etmp,1,sum)/1e6,type="l",lty=3,col="red")
		# dat1$avail_food_1


		prey2_eaten<-apply(tmp1_2$Etmp,1,sum)+apply(tmp2_2$Etmp,1,sum)+apply(tmp3_2$Etmp,1,sum)
		prey2_M2<-apply(tmp1_2$Mtmp,1,sum)+apply(tmp2_2$Mtmp,1,sum)+apply(tmp3_2$Mtmp,1,sum)
		#plot(prey2_eaten/1e6,ylim=c(0,2),type="l")
		plot(ceattle_2$PCOD_consumed[1,]/1e6,ylim=c(0,2),type="l")

		lines(ceattle_2$Biomass_2[1,]/1e6,type="l",lty=2)
		lines(apply(tmp1_2$Etmp,1,sum)/1e6,type="l",lty=1,col="red")
		lines(apply(tmp2_2$Etmp,1,sum)/1e6,type="l",lty=2,col="red")
		lines(apply(tmp3_2$Etmp,1,sum)/1e6,type="l",lty=3,col="red")

		prey3_eaten<-apply(tmp1_3$Etmp,1,sum)+apply(tmp2_3$Etmp,1,sum)+apply(tmp3_3$Etmp,1,sum)
		prey3_M2<-apply(tmp1_3$Mtmp,1,sum)+apply(tmp2_3$Mtmp,1,sum)+apply(tmp3_3$Mtmp,1,sum)
		#plot(prey3_eaten/1e6,ylim=c(0,2),type="l")
		plot(ceattle_2$ATF_consumed[1,]/1e6,ylim=c(0,2),type="l")
		lines(ceattle_2$Biomass_3[1,]/1e6,type="l",lty=2)
		lines(apply(tmp1_3$Etmp,1,sum)/1e6,type="l",lty=1,col="red")
		lines(apply(tmp2_3$Etmp,1,sum)/1e6,type="l",lty=2,col="red")
		lines(apply(tmp3_3$Etmp,1,sum)/1e6,type="l",lty=3,col="red")
		prey1_eaten<-ceattle_2$PLK_consumed[1,]
		prey2_eaten<-ceattle_2$PCOD_consumed[1,]
		prey3_eaten<-ceattle_2$ATF_consumed[1,]

		plot(100*prey1_eaten/ceattle_2$Biomass_1[1,],type="l",ylim=c(0,105))
		plot(100*prey2_eaten/ceattle_2$Biomass_2[1,],type="l",ylim=c(0,105))
		plot(100*prey3_eaten/ceattle_2$Biomass_3[1,],type="l",ylim=c(0,305))

		plot(prey1_M2,type="l",ylim=c(0,4),lty=1);lines(prey2_M2,lty=2);lines(prey3_M2,lty=3)
		pred<-1
		eval(parse(text=paste("nages_pred<-ceattle_2$nages_",pred,sep="")))

		suit_main_tmp<-matrix(0,nages_pred,21)
		for(pred_age in 1:nages_pred)
		{
				eval(parse(text=paste("suit_main<-ceattle_2$suit_main_",pred,"_",pred_age,sep="")))
				suit_main_tmp[pred_age,]<-colSums(suit_main)
		}

		for(prey_sp in 1:nspp){
				M2_est<-list()
				E_est<-list()
				Pred_demand_all_MSM<-list()
				Pred_demand_all_SS<-list()
				for(pred_sp in 1:nspp){
					M2_est[[pred_sp]]<-calc_eaten(pred=pred_sp,prey=prey_sp)$Mtmp
					E_est[[pred_sp]]<-calc_eaten(pred=pred_sp,prey=prey_sp)$Etmp
					Pred_demand_all_MSM[[pred_sp]]<-calc_eaten(dat1=ceattle_2,pred=pred_sp,prey=prey_sp)$Pred_demand
					Pred_demand_all_SS[[pred_sp]]<-calc_eaten(dat1=ceattle_0,pred=pred_sp,prey=prey_sp)$Pred_demand
				}
				eval(parse(text=paste("M2_est_",prey_sp,"<-M2_est",sep="")))
				eval(parse(text=paste("E_est_",prey_sp,"<-E_est",sep="")))
				eval(parse(text=paste("Pred_demand_all_MSM_",prey_sp,"<-Pred_demand_all_MSM",sep="")))
				eval(parse(text=paste("Pred_demand_all_SS_",prey_sp,"<-Pred_demand_all_SS",sep="")))
				eval(parse(text=paste("E_est_all_",prey_sp,"<-E_est[[1]]+E_est[[2]]+E_est[[3]]",sep="")))
				eval(parse(text=paste("M2_est_all_",prey_sp,"<-M2_est[[1]]+M2_est[[2]]+M2_est[[3]]",sep="")))
				
				eval(parse(text=paste("Pred_demand_yr_MSM_",prey_sp,"<-rowSums(Pred_demand_all_MSM[[1]])+rowSums(Pred_demand_all_MSM[[2]])+rowSums(Pred_demand_all_MSM[[3]])",sep="")))
				eval(parse(text=paste("Pred_demand_yr_SS_",prey_sp,"<-rowSums(Pred_demand_all_SS[[1]])+rowSums(Pred_demand_all_SS[[2]])+rowSums(Pred_demand_all_SS[[3]])",sep="")))
		}	

		plot_C<-function(dat=ceattle_2,ylimm_all=rbind(c(0,1000),c(0,20000),c(0,10000)) ) {
			quartz(height=6,width=4)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(c(1,2,3))
			layout.show(3)
			par(mfrow=c(3,1))

			plot(ceattle_2$Consum_1[1,],type="l",ylim=ylimm_all[1,])
			mn<-apply(ceattle_2$Consum_1,2,mean.na)
			sd<-apply(ceattle_2$Consum_1,2,sd.na)
			se<-apply(ceattle_2$Consum_1,2,se.na)
			nn<-length(mn)
			xx<-c(1:nn,nn:1);
			yy<-c(mn+1.95*se,(mn-1.95*se)[nn:1])
			yy2<-c(mn+1.95*sd,(mn-1.95*sd)[nn:1])
			polygon(xx,yy2,border=FALSE,col=colors()[330])
			for(i in 1:34) lines(ceattle_2$Consum_1[i,],col=colors()[320])
			polygon(xx,yy,border=FALSE,col=colors()[310])
			lines(1:nn,mn,lwd=2)
			

			plot(ceattle_2$Consum_2[1,],type="l",ylim=ylimm_all[2,])
			mn<-apply(ceattle_2$Consum_2,2,mean.na)
			sd<-apply(ceattle_2$Consum_2,2,sd.na)
			se<-apply(ceattle_2$Consum_2,2,se.na)
			nn<-length(mn)
			xx<-c(1:nn,nn:1);
			yy<-c(mn+1.95*se,(mn-1.95*se)[nn:1])
			yy2<-c(mn+1.95*sd,(mn-1.95*sd)[nn:1])
			polygon(xx,yy2,border=FALSE,col=colors()[330])
				for(i in 1:34) lines(ceattle_2$Consum_2[i,],col=colors()[320])
			polygon(xx,yy,border=FALSE,col=colors()[310])
			lines(1:nn,mn,lwd=2)

			plot(ceattle_2$Consum_3[1,],type="l",ylim=ylimm_all[3,])
			mn<-apply(ceattle_2$Consum_3,2,mean.na)
			sd<-apply(ceattle_2$Consum_3,2,sd.na)
			se<-apply(ceattle_2$Consum_3,2,se.na)
			nn<-length(mn)
			xx<-c(1:nn,nn:1);
			yy<-c(mn+1.95*se,(mn-1.95*se)[nn:1])
			yy2<-c(mn+1.95*sd,(mn-1.95*sd)[nn:1])
			polygon(xx,yy2,border=FALSE,col=colors()[330])
				for(i in 1:34) lines(ceattle_2$Consum_3[i,],col=colors()[320])
			polygon(xx,yy,border=FALSE,col=colors()[310])
			lines(1:nn,mn,lwd=2)
		}
		plot_C()
		if(update.figs==1) quartz.save(file=file.path(plot_fileJPG,"plot_C.pdf"),type="pdf",dpi=500) #(dev.copy2eps(device = quartz, file =file.path(plot_file,"plot_C.eps")))

		plot_suit<-function(maxA=5,dat1){
			quartz(height=6,width=4)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(c(1,2,3))
			layout.show(3)
				spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")

			for(prey in 1:3)
			{
				eval(parse(text=paste("nages_prey<-dat1$nages_",prey,sep="")))
				for(pred in 1:3)
				{
					eval(parse(text=paste("nages_pred<-dat1$nages_",pred,sep="")))
					suit_main_tmp<-matrix(0,nages_pred,nages_prey)
					for(pred_age in 1:nages_pred)
					{
						eval(parse(text=paste("suit_main<-dat1$suit_main_",pred,"_",pred_age,sep="")))
						suit_main_tmp[pred_age,]<-suit_main[prey,1:nages_prey]
					}	
					eval(parse(text=paste("suit_main_",pred,"<-suit_main_tmp",sep="")))
				}
				#plot(rowSums(suit_main_1),type="l",ylim=c(0,1),xlim=c(0,21),lwd=2)
				#lines(rowSums(suit_main_2),lty=2,lwd=2)
				#lines(rowSums(suit_main_3),lty=1,lwd=2,col=colors()[330])

				eval(parse(text=paste("suit_main_1_",prey,"<-suit_main_1",sep="")))
				eval(parse(text=paste("suit_main_2_",prey,"<-suit_main_2",sep="")))
				eval(parse(text=paste("suit_main_3_",prey,"<-suit_main_3",sep="")))
			}
			for(pred in 1:nspp){
				eval(parse(text=paste("nages_pred<-dat1$nages_",pred,sep="")))
				eval(parse(text=paste("suit_1<-suit_main_",pred,"_1",sep="")))
				eval(parse(text=paste("suit_2<-suit_main_",pred,"_2",sep="")))
				eval(parse(text=paste("suit_3<-suit_main_",pred,"_3",sep="")))
				eval(parse(text=paste("suit_other<-ceattle_2$suit_other_",pred,sep="")))
				if(maxA==1)
					plot((suit_1[,1]),type="l",ylim=c(0,1.1),xlim=c(1,nages_pred),lwd=2,axes=FALSE)
				if(maxA==0)
					plot(rowSums(suit_1),type="l",ylim=c(0,1.1),xlim=c(1,nages_pred),lwd=2,axes=FALSE)
				if(maxA!=0&maxA!=1)
					plot(rowSums(suit_1[,1:maxA]),type="l",ylim=c(0,1.1),xlim=c(1,nages_pred),lwd=2,axes=FALSE)
				
				axis(1)
				axis(1, at=c(-10,30))
				axis(2,las=2,at=seq(0,1,.25))
				axis(2, at=c(-10,10))
				lines(suit_other[1,],lty=2,lwd=2,col=colors()[330])
				if(maxA==1){
					lines((suit_3[,1]),lty=1,lwd=2,col=colors()[330])
					lines((suit_2[,1]),lty=2,lwd=2)
				}
					
				if(maxA==0){
					lines(rowSums(suit_3),lty=1,lwd=2,col=colors()[330])
					lines(rowSums(suit_2),lty=2,lwd=2)
				}
					
				if(maxA!=0&maxA!=1){
					lines(rowSums(suit_3[,1:maxA]),lty=1,lwd=2,col=colors()[330])
					lines(rowSums(suit_2[,1:maxA]),lty=2,lwd=2)
				}

				if(pred==2){
					legend(5,1.1,
						c("Pollock","P. cod","Arrowtooth","Other"),
						col=colors()[c(24,24,330,330)],
						lty=c(1,2,1,2),lwd=2,box.lty=0)
				}
				if(pred==3){text(.25,1.1,spnames[pred],pos=4,font=2,cex=1.2)}else{
				text(0.5,1.08,spnames[pred],pos=4,font=2,cex=1.2)}
			}
			if(maxA==0)
				mtext(side=2,"Total suitability",outer=TRUE,font=2,line=2)
			if(maxA!=0)
				mtext(side=2,paste("Total suitability (prey ages 1-",maxA,")",sep=""),outer=TRUE,font=2,line=2)
			mtext(side=1,"Predator age",outer=TRUE,font=2,line=.5)
		}
		plot_suit(maxA=0,dat1=ceattle_2)
		if(update.figs==1)quartz.save(file=file.path(plot_fileJPG,"Suitability.pdf"),type="pdf",dpi=500) #(dev.copy2eps(device = quartz, file =file.path(plot_file,"Suitability.eps")))
		collrange1<-colorRampPalette(colors()[c(24,330)]);collrange<-collrange1

		plot_suit2<-function(ltyy=c(2,1),coll=rep(collrange1(2)[1],2),lwdd=c(2,2),maxA=0,datt1=ceattle_2){
			quartz(height=5,width=6)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
			layout(rbind(c(1,4,7),c(2,5,8),c(3,6,9)))
			layout.show(9)
			yatt<-seq(0,2,.5)
			xatt<-seq(0,20,4)
			#spnames<-c("W. pollock","P. cod","Arrowtooth")
			lab_txt<-rbind(
				c("a)","b)","c)"),
				c("d)","e)","f)"),
				c("g)","h)","i)"))
			spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")
			for(prey in 1:3)
			{
				eval(parse(text=paste("nages_prey<-datt1$nages_",prey,sep="")))
				for(pred in 1:3)
				{
					eval(parse(text=paste("nages_pred<-datt1$nages_",pred,sep="")))
					suit_main_tmp<-matrix(0,nages_pred,nages_prey)
					for(pred_age in 1:nages_pred)
					{
						eval(parse(text=paste("suit_main<-datt1$suit_main_",pred,"_",pred_age,sep="")))
						suit_main_tmp[pred_age,]<-suit_main[prey,1:nages_prey]
					}	
					eval(parse(text=paste("suit_main_",pred,"<-suit_main_tmp",sep="")))
				}
				eval(parse(text=paste("suit_main_1_",prey,"<-suit_main_1",sep="")))
				eval(parse(text=paste("suit_main_2_",prey,"<-suit_main_2",sep="")))
				eval(parse(text=paste("suit_main_3_",prey,"<-suit_main_3",sep="")))
			}
			for(s in 1:nspp)
			{
				eval(parse(text=paste("dat1<-ceattle_0$fsh_sel_",s,sep="")))
				eval(parse(text=paste("dat2<-ceattle_2$fsh_sel_",s,sep="")))

				plot(dat1[1,],type="l",axes=FALSE,col=coll[1],lwd=lwdd[1],lty=ltyy[1])
				lines(dat2[1,],type="l",col=coll[2],lwd=lwdd[2],lty=ltyy[2])
				mtext(lab_txt[1,s],side=3,line=-1,adj=.05,outer=FALSE,font=2,cex=.9)
				if(s==1){
					mtext("Fishery selectivity",side=2,line=3,outer=FALSE,font=2)
					axis(2,at=yatt,font=1,las=2)
				}else{
					axis(2,at=yatt,lab=rep("",length(yatt)),font=1,las=2)
				}
				mtext(spnames[s],side=3,outer=FALSE,font=2)
				axis(1,at=xatt,lab=rep("",length(xatt)),font=1);
				axis(1,at=c(-1,100));axis(2,at=c(-1,10))
						eval(parse(text=paste("dat1<-ceattle_0$srv_sel_",s,sep="")))
				eval(parse(text=paste("dat2<-ceattle_2$srv_sel_",s,sep="")))

				plot(dat1[1,],type="l",axes=FALSE,col=coll[1],lwd=lwdd[1],lty=ltyy[1])
				lines(dat2[1,],type="l",col=coll[2],lwd=lwdd[2],lty=ltyy[2])
				mtext(lab_txt[2,s],side=3,line=-1,adj=.05,outer=FALSE,font=2,cex=.9)
				if(s==1){
					mtext("Survey selectivity",side=2,line=3,outer=FALSE,font=2)
					axis(2,at=yatt,font=1,las=2)
				}else{
					axis(2,at=yatt,lab=rep("",length(yatt)),font=1,las=2)
				}
				axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
				axis(1,at=c(-1,100));axis(2,at=c(-1,10))
				if(s==3)
					legend(3,.9,c("single species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)

				#for(pred in 1:nspp){
				pred<-s
				lwdd2<-c(2,2,2,2)
				cols2<-c(collrange1(2)[1:2],collrange1(2)[1:2])
				ltyy2<-c(1,1,3,3)
				eval(parse(text=paste("nages_pred<-datt1$nages_",pred,sep="")))
				eval(parse(text=paste("suit_1<-suit_main_",pred,"_1",sep="")))
				eval(parse(text=paste("suit_2<-suit_main_",pred,"_2",sep="")))
				eval(parse(text=paste("suit_3<-suit_main_",pred,"_3",sep="")))
				eval(parse(text=paste("suit_other<-ceattle_2$suit_other_",pred,sep="")))
				if(maxA==1)
					plot((suit_1[,1]),type="l",ylim=c(0,1.1),xlim=c(1,nages_pred),lty=ltyy2[1],lwd=lwdd2[1],col=cols2[1],axes=FALSE)
				if(maxA==0)
					plot(rowSums(suit_1),type="l",ylim=c(0,1.1),xlim=c(1,nages_pred),lty=ltyy2[1],lwd=lwdd2[1],col=cols2[1],axes=FALSE)
				if(maxA!=0&maxA!=1)
					plot(rowSums(suit_1[,1:maxA]),type="l",ylim=c(0,1.1),xlim=c(1,nages_pred),lty=ltyy2[1],lwd=lwdd2[1],col=cols2[1],axes=FALSE)
				mtext(lab_txt[3,s],side=3,line=-1,adj=.05,outer=FALSE,font=2,cex=.9)
				
				lines(suit_other[1,],lty=ltyy2[4],lwd=lwdd2[4],col=cols2[4])
				if(maxA==1){
					lines((suit_3[,1]),lty=ltyy2[3],lwd=lwdd2[3],col=cols2[3])
					lines((suit_2[,1]),lty=ltyy2[2],lwd=lwdd2[2],col=cols2[2])
				}
					
				if(maxA==0){
					lines(rowSums(suit_3),lty=ltyy2[3],lwd=lwdd2[3],col=cols2[3])
					lines(rowSums(suit_2),lty=ltyy2[2],lwd=lwdd2[2],col=cols2[2])
				}
					
				if(maxA!=0&maxA!=1){
					lines(rowSums(suit_3[,1:maxA]),lty=ltyy2[3],lwd=lwdd2[3],col=cols2[3])
					lines(rowSums(suit_2[,1:maxA]),lty=ltyy2[2],lwd=lwdd2[2],col=cols2[2])
				}

				if(pred==2){
					legend(2,1.2,
						c("Pollock","P. cod","Arrowtooth","Other"),
						col=cols2,
						lty=ltyy2,lwd=lwdd2,box.lty=0)
					#lty=c(1,2,1,2),lwd=2,box.lty=0)
				}
				#if(pred==3){text(.25,1.1,spnames[pred],pos=4,font=2,cex=1.2)}else{
				#text(0.5,1.08,spnames[pred],pos=4,font=2,cex=1.2)}
				#}
				#axis(1)
				#axis(1, at=c(-10,30))
				#axis(2,las=2,at=seq(0,1,.5))
				#axis(2, at=c(-10,10))

				
				axis(1,at=xatt,font=1)
				axis(1,at=c(-1,100));axis(2,at=c(-1,10))

				axis(2,at=yatt,font=1,las=2)
				if(maxA==0)
				{
					if(s==1){
						mtext("Total suitability",side=2,line=3,outer=FALSE,font=2)
						axis(2,at=yatt,font=1,las=2)
					}else{
						axis(2,at=yatt,lab=rep("",length(yatt)),font=1,las=2)
					}
				}	
				if(maxA!=0)
				{
					if(s==1){
						mtext(paste("Total suitability (prey ages 1-",maxA,")",sep=""),side=2,line=3,outer=FALSE,font=2)
						axis(2,at=yatt,font=1,las=2)
					}else{
						axis(2,at=yatt,lab=rep("",length(yatt)),font=1,las=2)
					}
				}	
			}# end s
			mtext("Age",side=1,line=0,outer=TRUE,font=2)	
		}
		plot_suit2()
		if(update.figs==1){
			quartz.save(file=file.path(plot_fileJPG,"plot_suit2.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"plot_suit2.eps")))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"plot_suit2.jpg"),res=500,units="in",height=5,width=6);dev.off()
		}

		plot_eaten<-function(units=1e6){
			quartz(height=6,width=4)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(c(1,2))
			layout.show(2)
			
			xatt<-seq(1960,2050,10)
			ylimm<-c(4,4,4)
			spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
			years<-styr+1:nyrs
			collrange1_msm<-colorRampPalette(colors()[c(310,24,1,340)]) # grays
			all_avail<-rowSums(ceattle_2$avail_food_1)+rowSums(ceattle_2$avail_food_2)+rowSums(ceattle_2$avail_food_3)
			#plot(all_avail,type="l",ylim=c(0,max(all_avail)))
			plot(years,Pred_demand_yr_MSM_1,type="l",ylim=c(0,max(Pred_demand_yr_MSM_1)),axes=FALSE)
			yatt<-pretty(seq(0,max(Pred_demand_yr_MSM_1),1e5),4)
			axis(2,at=yatt,lab=yatt/units,las=2)
			axis(2,at=c(-1e10,1e10));axis(1);axis(1,at=c(-1000,5000))

			other_prey<-Pred_demand_yr_MSM_1-(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3))
			ny<-length(years)
			xx<-c(years,years[ny:1])

			yy<-c(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3)+other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[1],lty=0)

			yy<-c(rowSums(E_est_all_2)+rowSums(E_est_all_3)+other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[2],lty=0)

			yy<-c(rowSums(E_est_all_2)+other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[3],lty=0)

			yy<-c(other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[4],lty=0)
			lines(years,Pred_demand_yr_MSM_1,lwd=1)
			legend(2000,max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:4],c("Walleye pollock","Pacific cod","Arrowtooth flounder","other prey"),box.lty=0)
			mtext(side=2,"Prey consumed (million t)",outer=FALSE,line=2,cex=.9,font=2)
			text(1980,max(Pred_demand_yr_MSM_1),"a)",font=2)
			
			allM2<-(M2_est_1[[1]][,1])+(M2_est_1[[2]][,1])+(M2_est_1[[3]][,1])
			plot(years,allM2,type="l",ylim=c(0,max(allM2)),axes=FALSE)
			yatt<-pretty(seq(0,max(allM2),.1),4)
			axis(2,at=yatt,lab=yatt,las=2)
			axis(2,at=c(-1e10,1e10));axis(1);axis(1,at=c(-1000,5000))

			ny<-length(years)
			xx<-c(years,years[ny:1])
			yy<-c(allM2,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[1],lty=0)
			yy<-c((M2_est_1[[2]][,1])+(M2_est_1[[3]][,1]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[2],lty=0)
			yy<-c((M2_est_1[[3]][,1]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[3],lty=1)
			# lines(years,allE,lwd=1)
			text(1980,max(allM2),"b)",font=2)
			
			legend(2000,max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:4],c("Walleye pollock","Pacific cod","Arrowtooth flounder","other prey"),box.lty=0)
			mtext(side=2,"Walleye pollock M2 (age 1)",outer=FALSE,line=2,cex=.9,font=2)
		}
		graphics.off()
		plot_eaten()
		if(update.figs==1){
			quartz.save(file=file.path(plot_fileJPG,"Eaten.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"Eaten.eps")))
		}

		plot_eaten2<-function(units=1e6, new=T){
			if(new==T) quartz(height=6,width=4)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(c(1,2,3))
			layout.show(3)
			
			xatt<-seq(1960,2050,10)
			ylimm<-c(4,4,4)
			spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
			years<-styr+1:ceattle_2$nyrs
			collrange1_msm<-colorRampPalette(colors()[c(310,24,1,340)]) # grays
			all_avail<-rowSums(ceattle_2$avail_food_1)+rowSums(ceattle_2$avail_food_2)+rowSums(ceattle_2$avail_food_3)
			#plot(all_avail,type="l",ylim=c(0,max(all_avail)))
			plot(years,Pred_demand_yr_MSM_1,type="l",ylim=c(0,max(Pred_demand_yr_MSM_1)),axes=FALSE)
			yatt<-pretty(seq(0,max(Pred_demand_yr_MSM_1),1e5),4)
			axis(2,at=yatt,lab=yatt/units,las=2)
			axis(2,at=c(-1e10,1e10));axis(1);axis(1,at=c(-1000,5000))
			other_prey<-Pred_demand_yr_MSM_1-(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3))
			ny<-length(years)
			xx<-c(years,years[ny:1])
			yy<-c(Pred_demand_yr_MSM_1,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[1],lty=0)
			yy<-c(rowSums(Pred_demand_all_MSM_1[[2]])+rowSums(Pred_demand_all_MSM_1[[3]]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[2],lty=0)
			yy<-c(rowSums(Pred_demand_all_MSM_1[[3]]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[3],lty=1)
			lines(years,Pred_demand_yr_MSM_1,lwd=1)
			#legend(2000,max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:3],c("Walleye pollock","Pacific cod","Arrowtooth flounder"),box.lty=0)
			mtext(side=2,"Predator ration ",outer=FALSE,line=3,cex=.9,font=2)
			mtext(side=2,"(million t)",outer=FALSE,line=1.7,cex=.9,font=2)
			text(1980,max(Pred_demand_yr_MSM_1),"a)",font=2,cex=1.2)
			
			all_avail<-rowSums(ceattle_2$avail_food_1)+rowSums(ceattle_2$avail_food_2)+rowSums(ceattle_2$avail_food_3)
			#plot(all_avail,type="l",ylim=c(0,max(all_avail)))
			plot(years,Pred_demand_yr_MSM_1,type="l",ylim=c(0,max(Pred_demand_yr_MSM_1)),axes=FALSE)
			yatt<-pretty(seq(0,max(Pred_demand_yr_MSM_1),1e5),4)
			axis(2,at=yatt,lab=yatt/units,las=2)
			axis(2,at=c(-1e10,1e10));axis(1);axis(1,at=c(-1000,5000))
			other_prey<-Pred_demand_yr_MSM_1-(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3))
			ny<-length(years)
			xx<-c(years,years[ny:1])
			yy<-c(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3)+other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[1],lty=0)
			yy<-c(rowSums(E_est_all_2)+rowSums(E_est_all_3)+other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[2],lty=0)
			yy<-c(rowSums(E_est_all_2)+other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[3],lty=0)
			yy<-c(other_prey,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[4],lty=0)
			lines(years,Pred_demand_yr_MSM_1,lwd=1)
			legend(2000,max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:4],c("Walleye pollock","Pacific cod","Arrowtooth flounder","other prey"),box.lty=0)
			mtext(side=2,"Prey consumed",outer=FALSE,line=3,cex=.9,font=2)
			mtext(side=2,"(million t)",outer=FALSE,line=1.7,cex=.9,font=2)

			text(1980,max(Pred_demand_yr_MSM_1),"b)",font=2,cex=1.2)

			allM2<-(M2_est_1[[1]][,1])+(M2_est_1[[2]][,1])+(M2_est_1[[3]][,1])
			#allM2<-rowSums(ceattle_2$pred_E_1)
			plot(years,allM2,type="l",ylim=c(0,max(allM2)),axes=FALSE)
			yatt<-pretty(seq(0,max(allM2),.1),4)
			axis(2,at=yatt,lab=yatt,las=2)
			axis(2,at=c(-1e10,1e10));axis(1);axis(1,at=c(-1000,5000))

			ny<-length(years)
			xx<-c(years,years[ny:1])
			yy<-c(allM2,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[1],lty=0)
			yy<-c((M2_est_1[[2]][,1])+(M2_est_1[[3]][,1]),rep(0,ny))
			#yy<-c(rowSums(ceattle_2$pred_E_1[,2:3]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[2],lty=0)
			yy<-c((M2_est_1[[3]][,1]),rep(0,ny))
			#yy<-c((ceattle_2$pred_E_1[,3]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[3],lty=1)
			lines(years,allM2,lwd=1)
			text(1980,max(allM2),"c)",font=2,cex=1.2)
			
			legend(2000,max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:4],c("Walleye pollock","Pacific cod","Arrowtooth flounder","other prey"),box.lty=0)
			mtext(side=2,"Pollock M2",outer=FALSE,line=3,cex=.9,font=2)
			mtext(side=2,"(age 1)",outer=FALSE,line=1.7,cex=.9,font=2)
		}
		plot_eaten2(units=1e9)

		plot_eaten3<-function(units=1e9,new=T){
			#units of prey eatend are in Kg - or 1e-3 tons!

			if(new==T) quartz(height=6,width=4)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,4,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(c(1,2,3))
			layout.show(3)
			
			xatt<-seq(1960,2050,10)
			ylimm<-c(4,4,4)
			spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
			years<-styr+1:ceattle_2$nyrs
			collrange1_msm<-colorRampPalette(colors()[c(310,24,1,340)]) # grays
			all_avail<-rowSums(ceattle_2$avail_food_1)+rowSums(ceattle_2$avail_food_2)+rowSums(ceattle_2$avail_food_3)*1000
			#plot(all_avail,type="l",ylim=c(0,max(all_avail)))
			plot(years,Pred_demand_yr_MSM_1,type="l",ylim=c(0,max(Pred_demand_yr_MSM_1)),axes=FALSE)
			yatt<-pretty(seq(0,max(Pred_demand_yr_MSM_1),1e5),4)
			axis(2,at=yatt,lab=yatt/units,las=2)
			axis(2,at=c(-1e20,1e20));axis(1);axis(1,at=c(-1000,5000))
			other_prey<-Pred_demand_yr_MSM_1-(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3))
			ny<-length(years)
			xx<-c(years,years[ny:1])
			yy<-c(Pred_demand_yr_MSM_1,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[1],lty=0)
			yy<-c(rowSums(Pred_demand_all_MSM_1[[2]])+rowSums(Pred_demand_all_MSM_1[[3]]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[2],lty=0)
			yy<-c(rowSums(Pred_demand_all_MSM_1[[3]]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[3],lty=1)
			lines(years,Pred_demand_yr_MSM_1,lwd=1)
			#legend(2000,max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:3],c("Walleye pollock","Pacific cod","Arrowtooth flounder"),box.lty=0)
			mtext(side=2,"Predator ration ",outer=FALSE,line=3.8,cex=.9,font=2)
			mtext(side=2,paste0("( x",units*1e-3, " t)"),outer=FALSE,line=2.5,cex=.9,font=2)
				legend(2000,1.1*max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:4],border=c(collrange1_msm(4)[1:2],colors()[300],collrange1_msm(4)[4]),c("Walleye pollock","Pacific cod","Arrowtooth flounder","other prey"),box.lty=0)
			text(1980,max(Pred_demand_yr_MSM_1),"a)",font=2,cex=1.2)
			#plot(all_avail,type="l",ylim=c(0,max(all_avail)))
			plot(years,log(Pred_demand_yr_MSM_1),type="l",ylim=log(c(1000,1.1*max(Pred_demand_yr_MSM_1))),axes=FALSE,col="white")
			other_prey<-Pred_demand_yr_MSM_1-(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3))
			totalC_dat<-((ceattle_2$prey_consumed_1+ceattle_2$prey_consumed_2+ceattle_2$prey_consumed_3)[1,])*1000
			totalC<-(rowSums(E_est_all_1)+rowSums(E_est_all_2)+rowSums(E_est_all_3))
			
			other_prey<-Pred_demand_yr_MSM_1-totalC_dat
			ny<-length(years)
			xx<-c(years,years[ny:1])
			yy<-c((ceattle_2$prey_consumed_1+ceattle_2$prey_consumed_2+ceattle_2$prey_consumed_3)[1,]*1000+other_prey,rep(1e-15,ny))
			polygon(xx,log(yy),col=collrange1_msm(4)[4],lty=0)
			yy<-c((ceattle_2$prey_consumed_1+ceattle_2$prey_consumed_2+ceattle_2$prey_consumed_3)[1,]*1000,rep(1e-15,ny))
			polygon(xx,log(yy),col=collrange1_msm(4)[1],lty=0)
			
			yy<-c((ceattle_2$prey_consumed_2+ceattle_2$prey_consumed_3)[1,]*1000,rep(1e-15,ny))
			polygon(xx,log(yy),col=collrange1_msm(4)[2],lty=0)
			yy<-c(ceattle_2$prey_consumed_3[1,],rep(1e-15,ny))
			polygon(xx,log(yy),col=collrange1_msm(4)[3],lty=0)
			lines(years,Pred_demand_yr_MSM_1,lwd=1)
			#legend(2000,1.1*max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:4],c("Walleye pollock","Pacific cod","Arrowtooth flounder","other prey"),box.lty=0)
			
			mtext(side=2,"Prey consumed",outer=FALSE,line=3.8,cex=.9,font=2)
			mtext(side=2,paste0("( x",(units/1000)*1e-3," t)"),outer=FALSE,line=2.5,cex=.9,font=2)
			yatt<-c(1,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e12)
			axis(2,at=log(yatt),lab=yatt*1000/units,las=2)
			axis(2,at=c(-1e20,1e20));axis(1);axis(1,at=c(-1000,5000))
			text(1980,log(max(Pred_demand_yr_MSM_1)),"b)",font=2,cex=1.2)
			#text(1995,log(),"Pollock")

			allM2<-(M2_est_1[[1]][,1])+(M2_est_1[[2]][,1])+(M2_est_1[[3]][,1])
			#allM2<-rowSums(ceattle_2$pred_E_1)
			plot(years,allM2,type="l",ylim=c(0,max(allM2)),axes=FALSE)
			yatt<-pretty(seq(0,max(allM2),.1),4)
			axis(2,at=yatt,lab=yatt,las=2)
			axis(2,at=c(-1e20,1e20));axis(1);axis(1,at=c(-1000,5000))

			ny<-length(years)
			xx<-c(years,years[ny:1])
			yy<-c(allM2,rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[1],lty=0)
			yy<-c((M2_est_1[[2]][,1])+(M2_est_1[[3]][,1]),rep(0,ny))
			#yy<-c(rowSums(ceattle_2$pred_E_1[,2:3]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[2],lty=0)
			yy<-c((M2_est_1[[3]][,1]),rep(0,ny))
			#yy<-c((ceattle_2$pred_E_1[,3]),rep(0,ny))
			polygon(xx,yy,col=collrange1_msm(4)[3],lty=1)
			lines(years,allM2,lwd=1)
			text(1980,max(allM2),"c)",font=2,cex=1.2)
			
			#legend(2000,1.1*max(Pred_demand_yr_MSM_1),cex=.9,fill=collrange1_msm(4)[1:4],c("Walleye pollock","Pacific cod","Arrowtooth flounder","other prey"),box.lty=0)
			mtext(side=2,"Pollock M2",outer=FALSE,line=3,cex=.9,font=2)
			mtext(side=2,"(age 1)",outer=FALSE,line=1.7,cex=.9,font=2)
		}
	
		plot_eaten3(units=1e9)
		# pred_demand_1
		#lines(years,ceattle_2$M2_1[,1],col="red")

		if(update.figs==1){
			  quartz.save(file=file.path(plot_fileJPG,"Eaten2.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"Eaten2.eps")))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"Eaten2.jpg"),res=500,units="in",height=6,width=4);dev.off()
		}


		plot_propM<-function(){
			quartz(height=4,width=6)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(1)
			layout.show(1)
			mnpropM<-matrix(0,3,3)
			rownames(mnpropM)<-paste("prey",c("plk","pcod","atf"),sep="_")
			colnames(mnpropM)<-paste("pred",c("plk","pcod","atf"),sep="_")
			
			allM2<-(M2_est_2[[1]][,1])+(M2_est_2[[2]][,1])+(M2_est_2[[3]][,1])
			plot(years,(M2_est_2[[1]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,axes=FALSE,ylab="proportion of total M2",xlab="Year")
			axis(1);axis(2)
			lines(years,(M2_est_2[[2]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=2)
			lines(years,(M2_est_2[[3]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=1, col="gray")
			mnpropM[2,1]<-mean((M2_est_2[[1]][,1])/allM2)
			mnpropM[2,2]<-mean((M2_est_2[[2]][,1])/allM2)
			mnpropM[2,3]<-mean((M2_est_2[[3]][,1])/allM2)

			allM2<-(M2_est_3[[1]][,1])+(M2_est_3[[2]][,1])+(M2_est_3[[3]][,1])
			plot(years,(M2_est_3[[1]][,1])/allM2,ylim=c(0,1),type="l",lwd=2)
			lines(years,(M2_est_3[[2]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=2)
			lines(years,(M2_est_3[[3]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=1, col="gray")
			mnpropM[3,1]<-mean((M2_est_3[[1]][,1])/allM2)
			mnpropM[3,2]<-mean((M2_est_3[[2]][,1])/allM2)
			mnpropM[3,3]<-mean((M2_est_3[[3]][,1])/allM2)


			allM2<-(M2_est_1[[1]][,1])+(M2_est_1[[2]][,1])+(M2_est_1[[3]][,1])
			plot(years,(M2_est_1[[1]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,axes=FALSE,ylab="proportion of total M2",xlab="Year")
			lines(years,(M2_est_1[[2]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=2)
			lines(years,(M2_est_1[[3]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=3)
			mnpropM[1,1]<-mean((M2_est_1[[1]][,1])/allM2)
			mnpropM[1,2]<-mean((M2_est_1[[2]][,1])/allM2)
			mnpropM[1,3]<-mean((M2_est_1[[3]][,1])/allM2)
			axis(1);axis(1,at=c(1900,2050))
			axis(2,las=2);axis(2,at=c(-2,2))
			mtext(side=2,"Proportion of total M2",font=2,line=2)
			mtext(side=1,"Year",font=2,line=2)
			legend(1996,1,c("Walleye pollock","Pacific cod","Arrowtooth flounder"),box.lty=0,lty=1:3,lwd=2)
			return(mnpropM)
		}
		
		plot_propM()
		tmpp<-plot_propM()
		if(update.figs==1){
			 quartz.save(file=file.path(plot_fileJPG,"plot_propM.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"plot_propM.eps")))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"plot_propM.jpg"),res=500,units="in",height=4,width=6);dev.off()
		}
		plot_propM2<-function(){
			quartz(height=4,width=6)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(1)
			layout.show(1)
			par(mfrow=c(3,1))
			mnpropM<-matrix(0,3,3)
			rownames(mnpropM)<-paste("prey",c("plk","pcod","atf"),sep="_")
			colnames(mnpropM)<-paste("pred",c("plk","pcod","atf"),sep="_")
			for(prey in 1:3){
				mest<-eval(parse(text=paste("M2_est_",prey,sep="")))
				allM2<-(mest[[1]][,1])+(mest[[2]][,1])+(mest[[3]][,1])
				plot(years,(mest[[1]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,axes=FALSE,ylab="proportion of total M2",xlab="Year")
				lines(years,(mest[[2]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=2)
				lines(years,(mest[[3]][,1])/allM2,ylim=c(0,1),type="l",lwd=2,lty=3)
				mnpropM[prey,1]<-mean((mest[[1]][,1])/allM2)
				mnpropM[prey,2]<-mean((mest[[2]][,1])/allM2)
				mnpropM[prey,3]<-mean((mest[[3]][,1])/allM2)
				axis(2,las=2);axis(2,at=c(-2,2))
				if(prey==1) legend(1996,1,c("Walleye pollock","Pacific cod","Arrowtooth flounder"),box.lty=0,lty=1:3,lwd=2)
			}
			axis(1);axis(1,at=c(1900,2050))
			
			mtext(side=2,"Proportion of total M2",font=2,line=2)
			mtext(side=1,"Year",font=2,line=2)
			
			return(mnpropM)
		}
		plot_propM2()
		tmpp2<-plot_propM2()
		if(update.figs==1){
			 quartz.save(file=file.path(plot_fileJPG,"plot_propM2.pdf"),type="pdf",dpi=500)
			#(dev.copy2eps(device = quartz, file =file.path(plot_file,"plot_propM.eps")))
			#dev.copy(jpeg,filename=file.path(plot_fileJPG,"plot_propM.jpg"),res=500,units="in",height=4,width=6);dev.off()
		}




	############################################################################################
	##  Liklihood plots
	############################################################################################
graphics.off()
dev.new()
fun_BbyAge<-function(dat=ceattle_2,dist=-1.5,div=1e6){
	yrs<-dat$styr+1:dat$nyrs-1
	rev(yrs)[1]
	byA<-dat$biomassByage_1/div
	nages<-dim(byA)[2]
	nyrs<-dat$nyrs[[1]]
	ydat<-(1:nyrs)*dist


plot(x,y)
lines(predict(lo), col='red', lwd=2)

	y<-nyrs
	x<-1:nages
	y<-byA[y,]+ydat[y]
	lo <- smooth.spline(x, y, spar=0.1,all.knots=T)
	ysmooth<-predict(lo)
	plot(y)
	lines(ysmooth)

	barplot(x,y,ylim=c(min(ydat),2),border="white",col="black",space=0)

	polygon(c(x,c(rev(nages)[1],1)),c(byA[y,]+ydat[y],ydat[y],ydat[y]),col="black",border="white")

	for(y in 1:nyrs){
		polygon(c(1:nages,c(rev(nages)[1],1)),c(byA[y,]+ydat[y],ydat[y],ydat[y]),col="black",border="white")

		# polygon(c(1:nages,c(rev(nages)[1],1)),c(byA[y,]+ydat[y],0,0),col="black",border="white")
	}
	barplot(x,y,border="white",col="black",space=0)

	for(y in 1:nyrs){
		x<-1:nages
		y<-byA[y,]+ydat[y]
		barplot(x,y,ylim=c(min(ydat),2),border="white",col="black",space=0,add=T)

	}



}
if(save.plotdata) save.image(file=file.path(plot_file,"plotdata.Rdata"))
