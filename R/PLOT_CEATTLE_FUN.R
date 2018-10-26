### PLOT_CEATTLE_FUN.R
crain<-col<-colorRampPalette(colors()[c(92,71,257,461)])
harmonic<-function(x){
	return(length(x)/sum(1/x))	
}
plot_Temps<-function(coll=collrange1(2)){
		quartz(height=3,width=6)
		par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
		par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
		par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
		layout.show(1)
		cpoly<-colorRampPalette(c(coll[1],"white"))(10)
		yatt<-seq(0,2,.25)
		xatt<-seq(0,20,4)
		#spnames<-c("W. pollock","P. cod","Arrowtooth")
		spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")

		plot(as.numeric(names(BT.dat)),BT.dat, ylim=c(-2,6),axes=FALSE,xlim=c(1979,2050),type="l",lwd=2,col=coll[2])
		nn<-length(BT.dat)
		xx1<-c(as.numeric(names(BT.dat)),as.numeric(names(BT.dat))[nn:1])
		yy1<-c(BT.dat+1.95*BT.dat.se,(BT.dat-1.95*BT.dat.se)[nn:1])
		yy2<-c(BT.dat+BT.dat.sd,(BT.dat-1*BT.dat.sd)[nn:1])
		polygon(xx1,yy2,col=cpoly[8],border=FALSE)
		#polygon(xx1,yy1,col=cpoly[1],border=FALSE)
		
		nn2<-which(BT_fut$year==as.numeric(names(BT.dat))[nn])
		nn<-length(BT_fut$year)
		xx1<-c(BT_fut$year[nn2:nn],BT_fut$year[nn:nn2])
		yy1<-c((BT_fut$AVG+1.95*BT_fut$SD/sqrt(83))[nn2:nn],(BT_fut$AVG-1.95*BT_fut$SD/sqrt(83))[nn:nn2])
		yy2<-c((BT_fut$AVG+1*BT_fut$SD)[nn2:nn],(BT_fut$AVG-1*BT_fut$SD)[nn:nn2])
		polygon(xx1,yy2,col=cpoly[8],border=FALSE)
		#polygon(xx1,yy1,col=cpoly[1],border=FALSE)

		lines(BT_fut$year[nn2:nn],BT_fut$AVG[nn2:nn],col=coll[2],lwd=2,lty=2)
		lines(as.numeric(names(BT.dat)),BT.dat,col=coll[2],lwd=2,lty=1)
			axis(1,at=seq(1970,2052,10))
			axis(2,at=seq(-2,6,2),las=2)
		mtext("Year",side=1,line=0,outer=TRUE,font=2)	
		mtext(expression(bold(Temperature~(symbol(degree)~C))),side=2,line=2,outer=FALSE,font=2)	
		abline(v=BT_fut$year[nn2],lty=2)
	}
# plot_Temps()
read_estrep<-function (fn,printt=TRUE) {
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


read.rec<-function(fn){
	rec_dev<-matrix(0,nspp,nyrs)
	rec_dev.se<-matrix(0,nspp,nyrs)
	tmpt<-read.csv(fn,sep="")
	tmpr<-tmpt[tmpt[,2]=="rec_dev",3]
	ny<-length(tmpr)/3
	ln_mn_rec<-tmpt[tmpt[,2]=="ln_mn_rec",3]
	ln_mn_rec.se<-tmpt[tmpt[,2]=="ln_mn_rec",4]/sqrt(ny)
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
	return(list(rec_dev=rec_dev,rec_dev.se=rec_dev.se,logR_obs=logR_obs,
		logR_obs.plus=logR_obs.plus,logR_obs.minus=logR_obs.minus,
		logR_obs.5.plus=logR_obs.5.plus,logR_obs.5.minus=logR_obs.5.minus))
}
#tmp<-read.rec(fn="/Users/kkari/Dropbox/msm-master/ceattle_0/results/tmsm_est.std")
assign.rec<-function(target=ceattle_0,fn=file.path(root,"ceattle_0/results/tmsm_est.std")){
	tmp<-read.rec(fn)
	target$rec_dev<-tmp$rec_dev
	target$rec_dev.se<-tmp$rec_dev.se
	target$logR_obs<-tmp$logR_obs
	target$logR_obs.plus<-tmp$logR_obs.plus
	target$logR_obs.minus<-tmp$logR_obs.minus
	target$logR_obs.5.plus<-tmp$logR_obs.5.plus
	target$logR_obs.5.minus<-tmp$logR_obs.5.minus

	target$nages_1<-length(target$AvgN_1[1,])

	target$nages_2<-length(target$AvgN_2[1,])

	target$nages_3<-length(target$AvgN_3[1,])
	return(target)
}
plot_sel<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2)){
		quartz(height=4.5,width=6)
		par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
		par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
		par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
		layout(rbind(c(1,3,5),c(2,4,6)))
		layout.show(6)
		yatt<-seq(0,2,.25)
		xatt<-seq(0,20,4)
		#spnames<-c("W. pollock","P. cod","Arrowtooth")
		spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")
		for(s in 1:nspp)
		{
			eval(parse(text=paste("dat1<-ceattle_0$fsh_sel_",s,sep="")))
			eval(parse(text=paste("dat2<-ceattle_2$fsh_sel_",s,sep="")))

			plot(dat1[1,],type="l",axes=FALSE,col=coll[1],lwd=lwdd[1],lty=ltyy[1])
			lines(dat2[1,],type="l",col=coll[2],lwd=lwdd[2],lty=ltyy[2])
			
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
			if(s==1){
				mtext("Survey selectivity",side=2,line=3,outer=FALSE,font=2)
				axis(2,at=yatt,font=1,las=2)
			}else{
				axis(2,at=yatt,lab=rep("",length(yatt)),font=1,las=2)
			}
			axis(1,at=xatt,font=1)
			axis(1,at=c(-1,100));axis(2,at=c(-1,10))
			if(s==3)
				legend(3,.9,c("single species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)

			
		}	
		mtext("Age",side=1,line=0,outer=TRUE,font=2)	
}

plot_allBiomass<-function( coll=collrange1(3),coll2=c("black","black","black"),
							lwdd=c(2,2,1),lwdd2=c(2,1,2),ltyy=c(1,1,1),ltyy2=c(1,1,2),
							txt=c("single-species", "multi-species","data"),
							txt2=c("Total biomass", "Survey biomass", "Female spwaning biomass"),
							pchh=c(-1,-1,19),pchh2=c(-1,-1,-1),ylimm=c(2.6e7,2e6,1e6)){
		quartz(height=7,width=5)
		par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
		par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
		par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
		layout(c(1,2,3))
		layout.show(3)
		xatt<-seq(1960,2050,10)
		
		spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
		allyr<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
		for(s in 1:nspp)
		{
			eval(parse(text=paste("dat1<-ceattle_0$srv_bio_",s,sep="")))
			eval(parse(text=paste("dat2<-ceattle_2$srv_bio_",s,sep="")))
			eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_srv_biom_",s,sep="")))
			eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_srv_biom_",s,sep="")))

			eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
			eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

			eval(parse(text=paste("dat1.hat<-ceattle_0$srv_bio_hat_",s,sep="")))
			eval(parse(text=paste("dat2.hat<-ceattle_2$srv_bio_hat_",s,sep="")))

			eval(parse(text=paste("dat1_B<-ceattle_0$Biomass_",s,sep="")))
			eval(parse(text=paste("dat2_B<-ceattle_2$Biomass_",s,sep="")))
			eval(parse(text=paste("dat1_SSB<-ceattle_0$BiomassSSB_",s,sep="")))
			eval(parse(text=paste("dat2_SSB<-ceattle_2$BiomassSSB_",s,sep="")))

			plot(allyr,as.numeric(dat2_B[1,]),type="l",axes=FALSE,col=coll[3],lwd=2,lty=1,pch=pchh[3],cex=1.2,ylim=c(0,ylimm[s]))
			lines(allyr,as.numeric(dat1_B[1,]),col=coll[1],lwd=2,lty=1)
			lines(allyr,dat2_B[1,],col=coll[2],lwd=2,lty=1)

			lines(allyr,dat1_SSB[1,],col=coll[1],lwd=2,lty=2)
			lines(allyr,dat2_SSB[1,],col=coll[2],lwd=2,lty=2)
			lines(dat1.yr[1,],dat1.hat[1,],col=coll[1],lwd=1,lty=1,cex=1.2)
			lines(dat2.yr[1,],dat2.hat[1,],lwd=1,lty=1,cex=1.2,col=coll[2])

			points(dat2.yr[1,],dat2[1,],type="p",col=coll[3],lwd=1,lty=ltyy[3],pch=pchh[3],cex=1.2,ylim=c(0,ylimm[s]))
			
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

			if(s==2)
				mtext("Biomass (million t)",side=2,line=3,outer=FALSE,font=2)
			text(ceattle_0$styr[1,1],ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
			if(s<3)
				axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
			if(s==3)
				axis(1,at=xatt,font=1)
			if(s==2)
				legend("topright",txt,col=coll,lwd=lwdd,lty=ltyy,pch=pchh,box.lty=0,cex=1.2)
			if(s==1)
				legend("topright",txt2,col=coll2,lwd=lwdd2,lty=ltyy2,pch=pchh2,box.lty=0,cex=1.2)
		}	
		mtext("Year",side=1,line=0,outer=TRUE,font=2)	
}

plot_srvB<-function(ltyy=c(1,2,1),coll=collrange1(3),lwdd=c(2,2,1),txt=c("single-species", "multi-species","data"),pchh=c(-1,-1,19)){
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
			eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_srv_biom_",s,sep="")))
			eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_srv_biom_",s,sep="")))

			eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
			eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

			eval(parse(text=paste("dat1.hat<-ceattle_0$srv_bio_hat_",s,sep="")))
			eval(parse(text=paste("dat2.hat<-ceattle_2$srv_bio_hat_",s,sep="")))

			plot(dat2.yr[1,],dat2[1,],type="p",axes=FALSE,col=coll[3],lwd=1,lty=ltyy[2],pch=pchh[3],cex=1.2,ylim=c(0,ylimm[s]))
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

			lines(dat2.yr[1,],dat2.hat[1,],col=coll[2],lwd=lwdd[2],lty=ltyy[2],pch=pchh[2],cex=1.2)
			lines(dat1.yr[1,],dat1.hat[1,],col=coll[1],lwd=lwdd[1],lty=ltyy[1],pch=pchh[2],cex=1.2)
			
			if(s==2)
				mtext("Survey biomass (million t)",side=2,line=3,outer=FALSE,font=2)
			text(1981,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
			if(s<3)
				axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
			if(s==3)
				axis(1,at=xatt,font=1)
			if(s==3)
				legend("topright",txt,col=coll,lwd=lwdd,lty=ltyy,pch=pchh,box.lty=0,cex=1.2)
		}	
		mtext("Year",side=1,line=0,outer=TRUE,font=2)	
}

plot_totalBiomass_VarVSnoT<-function(ltyy=c(1,2),coll=col1(2),coll2=col2(2),lwdd=c(2,2),ylimm=c(2.6e7,2e6,1e6)){
	quartz(height=6,width=6)
	par(mar=c(2,2,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
	layout(cbind(c(1,3,5),c(2,4,6)))#,c(7,8,9)))
	layout.show(6)
	xatt<-seq(1960,2050,10)
	
	spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")
	allyr<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$srv_bio_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$srv_bio_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat1.hat<-ceattle_0$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("dat1_B<-ceattle_0$Biomass_",s,sep="")))
		eval(parse(text=paste("dat2_B<-ceattle_2$Biomass_",s,sep="")))
		eval(parse(text=paste("dat1_SSB<-ceattle_0$BiomassSSB_",s,sep="")))
		eval(parse(text=paste("dat2_SSB<-ceattle_2$BiomassSSB_",s,sep="")))


		eval(parse(text=paste("noTdat1<-ceattle_0noT$srv_bio_",s,sep="")))
		eval(parse(text=paste("noTdat2<-ceattle_2noT$srv_bio_",s,sep="")))
		eval(parse(text=paste("noTdat1.yr<-ceattle_0noT$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("noTdat2.yr<-ceattle_2noT$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("noTdat1.se<-ceattle_0noT$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("noTdat2.se<-ceattle_2noT$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("noTdat1.hat<-ceattle_0noT$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("noTdat2.hat<-ceattle_2noT$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("noTdat1_B<-ceattle_0noT$Biomass_",s,sep="")))
		eval(parse(text=paste("noTdat2_B<-ceattle_2noT$Biomass_",s,sep="")))
		eval(parse(text=paste("noTdat1_SSB<-ceattle_0noT$BiomassSSB_",s,sep="")))
		eval(parse(text=paste("noTdat2_SSB<-ceattle_2noT$BiomassSSB_",s,sep="")))

		plot(allyr,dat2_B[1,],type="l",axes=FALSE,col=coll[2],lwd=2,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		lines(allyr,noTdat1_B[1,],col=coll2[1],lwd=1,lty=1)
		lines(allyr,noTdat2_B[1,],col=coll2[2],lwd=1,lty=1)
		lines(allyr,dat1_B[1,],col=coll[1],lwd=2,lty=1)
		lines(allyr,dat2_B[1,],col=coll[2],lwd=2,lty=1)
		
		# plot(allyr,dat2_SSB[1,],type="l",axes=FALSE,col=coll[2],lwd=2,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		# lines(allyr,noTdat1_SSB[1,],col=coll2[1],lwd=1,lty=1)
		# lines(allyr,noTdat3_SSB[1,],col=coll2[2],lwd=1,lty=1)
		# lines(allyr,dat1_SSB[1,],col=coll[1],lwd=2,lty=1)
		# lines(allyr,dat3_SSB[1,],col=coll[2],lwd=2,lty=1)
		
		#lines(allyr,dat1_SSB[1,],type="l",axes=FALSE,col=coll[2],lwd=1,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]/2))
		#lines(allyr,noTdat1_SSB[1,],col=coll2[1],lwd=1,lty=2)
		#lines(allyr,noTdat2_SSB[1,],col=coll2[2],lwd=1,lty=2)
		#lines(allyr,dat1_SSB[1,],col=coll[1],lwd=1,lty=1)
		#lines(allyr,dat2_SSB[1,],col=coll[2],lwd=1,lty=1)
		
		# plot(dat1.yr[1,],dat1.hat[1,],type="l",axes=FALSE,col=coll[2],lwd=1,lty=2,pch=19,cex=1.2,ylim=c(0,ylimm[s]/2))
		# lines(dat1.yr[1,],noTdat1.hat[1,],col=coll2[1],lwd=1,lty=3)
		# lines(dat2.yr[1,],noTdat2.hat[1,],col=coll2[2],lwd=1,lty=3)
		# lines(dat1.yr[1,],dat1.hat[1,],col=coll[1],lwd=1,lty=2)
		# lines(dat2.yr[1,],dat2.hat[1,],col=coll[2],lwd=1,lty=2)
		
		#cc<-.1/1000
		#for(i in 1:length(dat2.yr[1,]))
		#{
		#	x1<-c((1-cc)*dat2.yr[1,i],(1+cc)*dat2.yr[1,i])
		#	y1<-dat2[1,i]+1.96*dat2.se[1,i]
		#	y2<-dat2[1,i]-1.96*dat2.se[1,i]
		#	lines(rep(dat2.yr[1,i],2),c(y1,y2))
		#	lines(x1,rep(y1,2))
		#	lines(x1,rep(y2,2))
		#}
		yatt<-pretty(seq(0,ylimm[s],1e4),4)
		axis(2,at=yatt,lab=yatt/1e6,las=2)
		if(s==1)
			mtext("Biomass (million t)",side=3,line=0,outer=FALSE,font=2)
		#text(ceattle_0$styr,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		mtext(spnames[s],side=2,line=3,outer=FALSE,font=2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==1)
			legend(2000,ylimm[s],c("SS with T", "MS with T","SS, no T","MS, no T"),col=c(coll,coll2),lwd=c(2,2,1,1),lty=c(1,1,1,1),box.lty=0,cex=1)

			ylimm2<-c(13,27,13)
		plot(allyr,((dat1_B[1,]-noTdat1_B[1,])/((noTdat1_B[1,])))*100,type="l",col=coll[1],lwd=2,lty=1,ylim=ylimm2[s]*c(-1,1),axes=FALSE)
		abline(h=0)
		lines(allyr,((dat2_B[1,]-noTdat2_B[1,])/((noTdat2_B[1,])))*100,type="l",col=coll[2],lwd=2,lty=1)
		
		yatt<-seq(-100,100,5)
		axis(2,at=yatt,lab=yatt,las=2)
		if(s==1)
			mtext("% Difference",side=3,line=0,outer=FALSE,font=2)
		text(ceattle_0$styr[1,1],ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)		
	}	
	mtext("Year",side=1,line=0,outer=TRUE,font=2)	
}

plot_SSBiomass_VarVSnoT<-function(ylimm2=c(13,50,13),ltyy=c(1,2),coll=col1(2),coll2=col2(2),lwdd=c(2,2),ylimm=c(.8e7,1.2e6,.5e6)){
	quartz(height=6,width=6)
	par(mar=c(2,2,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
	layout(cbind(c(1,3,5),c(2,4,6)))#,c(7,8,9)))
	layout.show(6)
	xatt<-seq(1960,2050,10)
	
	spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")
	allyr<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$srv_bio_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$srv_bio_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat1.hat<-ceattle_0$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("dat1_B<-ceattle_0$Biomass_",s,sep="")))
		eval(parse(text=paste("dat2_B<-ceattle_2$Biomass_",s,sep="")))
		eval(parse(text=paste("dat1_SSB<-ceattle_0$BiomassSSB_",s,sep="")))
		eval(parse(text=paste("dat2_SSB<-ceattle_2$BiomassSSB_",s,sep="")))


		eval(parse(text=paste("noTdat1<-ceattle_0noT$srv_bio_",s,sep="")))
		eval(parse(text=paste("noTdat2<-ceattle_2noT$srv_bio_",s,sep="")))
		eval(parse(text=paste("noTdat1.yr<-ceattle_0noT$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("noTdat2.yr<-ceattle_2noT$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("noTdat1.se<-ceattle_0noT$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("noTdat2.se<-ceattle_2noT$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("noTdat1.hat<-ceattle_0noT$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("noTdat2.hat<-ceattle_2noT$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("noTdat1_B<-ceattle_0noT$Biomass_",s,sep="")))
		eval(parse(text=paste("noTdat2_B<-ceattle_2noT$Biomass_",s,sep="")))
		eval(parse(text=paste("noTdat1_SSB<-ceattle_0noT$BiomassSSB_",s,sep="")))
		eval(parse(text=paste("noTdat2_SSB<-ceattle_2noT$BiomassSSB_",s,sep="")))

		# plot(allyr,dat2_B[1,],type="l",axes=FALSE,col=coll[2],lwd=2,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		# lines(allyr,noTdat1_B[1,],col=coll2[1],lwd=1,lty=1)
		# lines(allyr,noTdat2_B[1,],col=coll2[2],lwd=1,lty=1)
		# lines(allyr,dat1_B[1,],col=coll[1],lwd=2,lty=1)
		# lines(allyr,dat2_B[1,],col=coll[2],lwd=2,lty=1)
		
		plot(allyr,dat2_SSB[1,],type="l",axes=FALSE,col=coll[2],lwd=2,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		lines(allyr,noTdat1_SSB[1,],col=coll2[1],lwd=1,lty=1)
		lines(allyr,noTdat2_SSB[1,],col=coll2[2],lwd=1,lty=1)
		lines(allyr,dat1_SSB[1,],col=coll[1],lwd=2,lty=1)
		lines(allyr,dat2_SSB[1,],col=coll[2],lwd=2,lty=1)
		
		#lines(allyr,dat1_SSB[1,],type="l",axes=FALSE,col=coll[2],lwd=1,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]/2))
		#lines(allyr,noTdat1_SSB[1,],col=coll2[1],lwd=1,lty=2)
		#lines(allyr,noTdat2_SSB[1,],col=coll2[2],lwd=1,lty=2)
		#lines(allyr,dat1_SSB[1,],col=coll[1],lwd=1,lty=1)
		#lines(allyr,dat2_SSB[1,],col=coll[2],lwd=1,lty=1)
		
		# plot(dat1.yr[1,],dat1.hat[1,],type="l",axes=FALSE,col=coll[2],lwd=1,lty=2,pch=19,cex=1.2,ylim=c(0,ylimm[s]/2))
		# lines(dat1.yr[1,],noTdat1.hat[1,],col=coll2[1],lwd=1,lty=3)
		# lines(dat2.yr[1,],noTdat2.hat[1,],col=coll2[2],lwd=1,lty=3)
		# lines(dat1.yr[1,],dat1.hat[1,],col=coll[1],lwd=1,lty=2)
		# lines(dat2.yr[1,],dat2.hat[1,],col=coll[2],lwd=1,lty=2)
		
		#cc<-.1/1000
		#for(i in 1:length(dat2.yr[1,]))
		#{
		#	x1<-c((1-cc)*dat2.yr[1,i],(1+cc)*dat2.yr[1,i])
		#	y1<-dat2[1,i]+1.96*dat2.se[1,i]
		#	y2<-dat2[1,i]-1.96*dat2.se[1,i]
		#	lines(rep(dat2.yr[1,i],2),c(y1,y2))
		#	lines(x1,rep(y1,2))
		#	lines(x1,rep(y2,2))
		#}
		yatt<-pretty(seq(0,ylimm[s],1e4),4)
		axis(2,at=yatt,lab=yatt/1e6,las=2)
		if(s==1)
			mtext("Spawning biomass (million t)",side=3,line=0,outer=FALSE,font=2)
		#text(ceattle_0$styr,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		mtext(spnames[s],side=2,line=3,outer=FALSE,font=2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==1)
			legend(2000,ylimm[s],c("SS with T", "MS with T","SS, no T","MS, no T"),col=c(coll,coll2),lwd=c(2,2,1,1),lty=c(1,1,1,1),box.lty=0,cex=1)

			
		# plot(allyr,((dat1_B[1,]-noTdat1_B[1,])/((noTdat1_B[1,])))*100,type="l",col=coll[1],lwd=2,lty=1,ylim=ylimm2[s]*c(-1,1),axes=FALSE)
		# abline(h=0)
		# lines(allyr,((dat2_B[1,]-noTdat2_B[1,])/((noTdat2_B[1,])))*100,type="l",col=coll[2],lwd=2,lty=1)
		plot(allyr,((dat1_SSB[1,]-noTdat1_SSB[1,])/((noTdat1_SSB[1,])))*100,type="l",col=coll[1],lwd=2,lty=1,ylim=ylimm2[s]*c(-1,1),axes=FALSE)
		abline(h=0)
		lines(allyr,((dat2_SSB[1,]-noTdat2_SSB[1,])/((noTdat2_SSB[1,])))*100,type="l",col=coll[2],lwd=2,lty=1)
		

		cat("sp ",s," SS mean =",mean(abs(round((((dat1_SSB[1,]-noTdat1_SSB[1,])/((noTdat1_SSB[1,])))*100),2))),"\n")
		cat("sp ",s," SS range =",
		min(round((((dat1_SSB[1,]-noTdat1_SSB[1,])/((noTdat1_SSB[1,])))*100),2)),"-",
		max(round((((dat1_SSB[1,]-noTdat1_SSB[1,])/((noTdat1_SSB[1,])))*100),2)),"\n")
		
		cat("sp ",s," MSM mean =",mean(abs(round((((dat2_SSB[1,]-noTdat2_SSB[1,])/((noTdat2_SSB[1,])))*100),2))),"\n")
		cat("sp ",s," MSM range =",
		min(round((((dat1_SSB[1,]-noTdat2_SSB[1,])/((noTdat2_SSB[1,])))*100),2)),"-",
		max(round((((dat1_SSB[1,]-noTdat2_SSB[1,])/((noTdat2_SSB[1,])))*100),2)),"\n")

		yatt<-pretty(c(-1,1)*ylimm2[s],5)
		axis(2,at=yatt,lab=yatt,las=2)
		if(s==1)
			mtext("% Difference",side=3,line=0,outer=FALSE,font=2)
		text(ceattle_0$styr[1,1],ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)		
	}	
	mtext("Year",side=1,line=0,outer=TRUE,font=2)	
}

plot_srvyBiomass_VarVSnoT<-function(ltyy=c(1,2),coll=col1(2),coll2=col2(2),lwdd=c(2,2),ylimm=c(1e7,1.2e6,1e6)){
	quartz(height=6,width=6)
	par(mar=c(2,2,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
	layout(cbind(c(1,3,5),c(2,4,6)))#,c(7,8,9)))
	layout.show(6)
	xatt<-seq(1960,2050,10)
	
	spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")
	allyr<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
	s<-1
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$srv_bio_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$srv_bio_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat1.hat<-ceattle_0$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("dat1_B<-ceattle_0$Biomass_",s,sep="")))
		eval(parse(text=paste("dat2_B<-ceattle_2$Biomass_",s,sep="")))
		eval(parse(text=paste("dat1_SSB<-ceattle_0$BiomassSSB_",s,sep="")))
		eval(parse(text=paste("dat2_SSB<-ceattle_2$BiomassSSB_",s,sep="")))


		eval(parse(text=paste("noTdat1<-ceattle_0noT$srv_bio_",s,sep="")))
		eval(parse(text=paste("noTdat2<-ceattle_2noT$srv_bio_",s,sep="")))
		eval(parse(text=paste("noTdat1.yr<-ceattle_0noT$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("noTdat2.yr<-ceattle_2noT$yrs_srv_biom_",s,sep="")))
		eval(parse(text=paste("noTdat1.se<-ceattle_0noT$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("noTdat2.se<-ceattle_2noT$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("noTdat1.hat<-ceattle_0noT$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("noTdat2.hat<-ceattle_2noT$srv_bio_hat_",s,sep="")))
		eval(parse(text=paste("noTdat1_B<-ceattle_0noT$Biomass_",s,sep="")))
		eval(parse(text=paste("noTdat2_B<-ceattle_2noT$Biomass_",s,sep="")))
		eval(parse(text=paste("noTdat1_SSB<-ceattle_0noT$BiomassSSB_",s,sep="")))
		eval(parse(text=paste("noTdat2_SSB<-ceattle_2noT$BiomassSSB_",s,sep="")))

		# plot(allyr,dat2_B[1,],type="l",axes=FALSE,col=coll[2],lwd=2,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		# lines(allyr,noTdat1_B[1,],col=coll2[1],lwd=1,lty=1)
		# lines(allyr,noTdat2_B[1,],col=coll2[2],lwd=1,lty=1)
		# lines(allyr,dat1_B[1,],col=coll[1],lwd=2,lty=1)
		# lines(allyr,dat2_B[1,],col=coll[2],lwd=2,lty=1)
		
		# plot(allyr,dat2_SSB[1,],type="l",axes=FALSE,col=coll[2],lwd=2,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		# lines(allyr,noTdat1_SSB[1,],col=coll2[1],lwd=1,lty=1)
		# lines(allyr,noTdat2_SSB[1,],col=coll2[2],lwd=1,lty=1)
		# lines(allyr,dat1_SSB[1,],col=coll[1],lwd=2,lty=1)
		# lines(allyr,dat2_SSB[1,],col=coll[2],lwd=2,lty=1)
		
		
		plot(dat2.yr[1,],dat2.hat[1,],type="l",axes=FALSE,col=coll[2],lwd=2,lty=1,pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		points(dat2.yr[1,],dat2[1,],type="p",axes=FALSE,col=colors()[320],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		cc<-.1/1000
		for(i in 1:length(dat2.yr[1,]))
		{
			x1<-c((1-cc)*dat2.yr[1,i],(1+cc)*dat2.yr[1,i])
			y1<-dat2[1,i]+1.96*dat2.se[1,i]
			y2<-dat2[1,i]-1.96*dat2.se[1,i]
			lines(rep(dat2.yr[1,i],2),c(y1,y2),col=colors()[320])
			lines(x1,rep(y1,2),col=colors()[320])
			lines(x1,rep(y2,2),col=colors()[320])
		}
		lines(dat1.yr[1,],noTdat1.hat[1,],col=coll2[1],lwd=1,lty=1)
		lines(dat2.yr[1,],noTdat2.hat[1,],col=coll2[2],lwd=1,lty=1)
		lines(dat1.yr[1,],dat1.hat[1,],col=coll[1],lwd=2,lty=1)
		lines(dat2.yr[1,],dat2.hat[1,],col=coll[2],lwd=2,lty=1)
		
		# plot(dat1.yr[1,],dat1.hat[1,],type="l",axes=FALSE,col=coll[2],lwd=1,lty=2,pch=19,cex=1.2,ylim=c(0,ylimm[s]/2))
		# lines(dat1.yr[1,],noTdat1.hat[1,],col=coll2[1],lwd=1,lty=3)
		# lines(dat2.yr[1,],noTdat2.hat[1,],col=coll2[2],lwd=1,lty=3)
		# lines(dat1.yr[1,],dat1.hat[1,],col=coll[1],lwd=1,lty=2)
		# lines(dat2.yr[1,],dat2.hat[1,],col=coll[2],lwd=1,lty=2)
		
		#cc<-.1/1000
		#for(i in 1:length(dat2.yr[1,]))
		#{
		#	x1<-c((1-cc)*dat2.yr[1,i],(1+cc)*dat2.yr[1,i])
		#	y1<-dat2[1,i]+1.96*dat2.se[1,i]
		#	y2<-dat2[1,i]-1.96*dat2.se[1,i]
		#	lines(rep(dat2.yr[1,i],2),c(y1,y2))
		#	lines(x1,rep(y1,2))
		#	lines(x1,rep(y2,2))
		#}
		yatt<-pretty(seq(0,ylimm[s],1e4),4)
		axis(2,at=yatt,lab=yatt/1e6,las=2)
		if(s==1)
			mtext("Survey biomass (million t)",side=3,line=0,outer=FALSE,font=2)
		#text(ceattle_0$styr,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		mtext(spnames[s],side=2,line=3,outer=FALSE,font=2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==1)
			legend(2000,ylimm[s],c("SS with T", "MS with T","SS, no T","MS, no T"),col=c(coll,coll2),lwd=c(2,2,1,1),lty=c(1,1,1,1),box.lty=0,cex=1)

			ylimm2<-c(13,27,13)
		# plot(allyr,((dat1_B[1,]-noTdat1_B[1,])/((noTdat1_B[1,])))*100,type="l",col=coll[1],lwd=2,lty=1,ylim=ylimm2[s]*c(-1,1),axes=FALSE)
		# abline(h=0)
		# lines(allyr,((dat2_B[1,]-noTdat2_B[1,])/((noTdat2_B[1,])))*100,type="l",col=coll[2],lwd=2,lty=1)
		# plot(allyr,((dat1_SSB[1,]-noTdat1_SSB[1,])/((noTdat1_SSB[1,])))*100,type="l",col=coll[1],lwd=2,lty=1,ylim=ylimm2[s]*c(-1,1),axes=FALSE)
		# abline(h=0)
		# lines(allyr,((dat2_SSB[1,]-noTdat2_SSB[1,])/((noTdat2_SSB[1,])))*100,type="l",col=coll[2],lwd=2,lty=1)
		
		plot(dat1.yr[1,],((dat1.hat[1,]-noTdat1.hat[1,])/((noTdat1.hat[1,])))*100,type="l",col=coll[1],lwd=2,lty=1,ylim=ylimm2[s]*c(-1,1),axes=FALSE)
		abline(h=0)
		lines(dat2.yr[1,],((dat2.hat[1,]-noTdat2.hat[1,])/((noTdat2.hat[1,])))*100,type="l",col=coll[2],lwd=2,lty=1)
	
		cat("sp ",s," SS mean =",mean(abs(round((((dat1.hat[1,]-noTdat1.hat[1,])/((noTdat1.hat[1,])))*100),2))),"\n")
		cat("sp ",s," SS range =",
		min(round((((dat1.hat[1,]-noTdat1.hat[1,])/((noTdat1.hat[1,])))*100),2)),"-",
		max(round((((dat1.hat[1,]-noTdat1.hat[1,])/((noTdat1.hat[1,])))*100),2)),"\n")
		
		cat("sp ",s," MSM mean =",mean(abs(round((((dat2.hat[1,]-noTdat2.hat[1,])/((noTdat2.hat[1,])))*100),2))),"\n")
		cat("sp ",s," MSM range =",
		min(round((((dat2.hat[1,]-noTdat2.hat[1,])/((noTdat2.hat[1,])))*100),2)),"-",
		max(round((((dat2.hat[1,]-noTdat2.hat[1,])/((noTdat2.hat[1,])))*100),2)),"\n")

	
		yatt<-seq(-100,100,5)
		axis(2,at=yatt,lab=yatt,las=2)
		if(s==1)
			mtext("% Difference",side=3,line=0,outer=FALSE,font=2)
		text(ceattle_0$styr[1,1],ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)		
	}	
	mtext("Year",side=1,line=0,outer=TRUE,font=2)	
}

plot_catchB<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2),ylimm=c(2e6,5.5e5,3e4)){
	quartz(height=6,width=4)
	par(mar=c(2,1,1,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	par(mfrow=c(nspp,1))
	xatt<-seq(1960,2050,10)
	
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
	years<-1978+1:ceattle_0$nyrs
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$tc_biom_obs_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$tc_biom_obs_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_fs_h_comp",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_fs_h_comp",s,sep="")))

		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

		eval(parse(text=paste("dat1.hat<-ceattle_0$obs_catch_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$obs_catch_hat_",s,sep="")))

		plot(years,dat2[1,],type="p",axes=FALSE,col=coll[2],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		yatt<-pretty(seq(0,ylimm[s],1e4),4)
		axis(2,at=yatt,lab=yatt/1e6,las=2)
		cc<-.1/1000
		
		lines(years,dat1.hat[1,],col=coll[1],lwd=lwdd[1],lty=ltyy[1],cex=1.2)
		lines(years,dat2.hat[1,],col=coll[2],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
		
		if(s==2)
			mtext("Total catch (million t)",side=2,line=3,outer=FALSE,font=2)
		text(1978,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==2)
			legend(1997,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
	}	
	mtext("Year",side=1,line=0,outer=TRUE,font=2)	
}

plot_srvAge<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2)){
	quartz(height=6,width=4)
	par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,2,1))# outer margins of graph: (bottom,left, top, right)
	layout(c(1,2,3))
	layout.show(3)
	xatt<-seq(0,25,5)
	ylimm<-c(.2,.2,.2)
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")

	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$srv_age_obs_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$srv_age_obs_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$srv_age_yrs_",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$srv_age_yrs_",s,sep="")))

		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

		eval(parse(text=paste("dat1.hat<-ceattle_0$srv_age_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$srv_age_hat_",s,sep="")))

		plot(apply(dat2,2,mean.na),type="p",axes=FALSE,col=coll[2],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		axis(2,las=2)
		yatt<-pretty(seq(0,ylimm[s],1e4),4)
		#axis(2,at=yatt,lab=yatt/1e6,las=2)

		lines(apply(dat1.hat,2,mean.na),col=coll[1],lwd=lwdd[1],lty=ltyy[1],cex=1.2)
		lines(apply(dat2.hat,2,mean.na),col=coll[2],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
		
		if(s==2)
			mtext("Survey age or length composition",side=2,line=3,outer=FALSE,font=2)
		text(1,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		
		axis(1,at=pretty(-10:length(dat1[1,]),10),font=1)
		if(s==1)
			legend(7,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
		if(ceattle_0$srv_age_type[s]==1)
			mtext("Age",side=1,line=1.5,outer=FALSE,font=2)
		if(ceattle_0$srv_age_type[s]==2)
			mtext("Length",side=1,line=1.5,outer=FALSE,font=2)	
	}			
}
plot_srvAge_byYr<-function(ltyy=c(1,2),coll=collrange1(3),lwdd=c(2,2)){
	quartz(height=6,width=8)
	par(mar=c(0,0,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(.5,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	par(mfrow=c(ceattle_0$nyrs[1],nspp))
	layout(cbind(
		c(1:ceattle_0$nyrs[1]),
		c(1:ceattle_0$nyrs[1]+ceattle_0$nyrs[1]),
		c(1:ceattle_0$nyrs[1]+2*ceattle_0$nyrs[1])
		))
	# layout.show(3)
	xatt<-seq(0,25,5)
	ylimm<-c(.4,.4,.4)
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
	

	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$srv_age_obs_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$srv_age_obs_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$srv_age_yrs_",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$srv_age_yrs_",s,sep="")))

		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

		eval(parse(text=paste("dat1.hat<-ceattle_0$srv_age_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$srv_age_hat_",s,sep="")))
		years<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
		cc<-match(dat1.yr,years)
		nn2<-setdiff(years,dat1.yr)
		for(y in 1:length(cc)) {
			maxx<-1.01*max(c(dat1[y,],dat1.hat[y,],dat2.hat[y,]),na.rm=T)
			plot(dat1[y,],pch=16,axes=F,type="b",col=coll[1],ylim=c(0,maxx))
			lines(dat1.hat[y,],col=coll[2])
			lines(dat2.hat[y,],col=coll[3])
			if(s==1) mtext(side=2,dat1.yr[y],cex=.8,line=-.5,las=2)
		}
		if(length(nn2)>0) for(y in nn2) plot(1,1,col=F,axes=F)
		# plot(apply(dat2,2,mean.na),type="p",axes=FALSE,col=coll[2],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		# axis(2,las=2)
		# yatt<-pretty(seq(0,ylimm[s],1e4),4)
		#axis(2,at=yatt,lab=yatt/1e6,las=2)
	

		if(ceattle_0$fsh_age_type[s]==1)
			mtext("Age",side=1,line=0,outer=FALSE,font=2)
		if(ceattle_0$fsh_age_type[s]==2)
			mtext("Length",side=1,line=0,outer=FALSE,font=2)
		if(s==2)
			mtext("Fishery age or length composition",side=4,line=1,outer=FALSE,font=2)
	}
		# text(1,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		# axis(1,at=pretty(-10:length(dat1[1,]),10),font=1)

		# if(s==1)
		# 	legend(7,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
}

plot_srvAge_byYr2<-function(ltyy=c(1,2),coll=collrange1(3),lwdd=c(2,2),s=1){
	quartz(height=6,width=8)
	par(mar=c(0,0,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(.5,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	par(mfrow=c(ceattle_0$nyrs[1],nspp))
	nyrs<-ceattle_0$nyrs[1]
	rowss<-ceiling(nyrs/3)

	layout(cbind(seq(1,nyrs,3),
	              seq(1,nyrs,3)+1,
	              seq(1,nyrs,3)+2))
	  
	layout(cbind(
	  (1:rowss),
	  (rowss+1):(2*rowss),
	  ((rowss*2)+1):(3*rowss)))
	layout.show(nyrs)
	xatt<-seq(0,25,5)
	ylimm<-c(.4,.4,.4)
	spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")
	

	
		eval(parse(text=paste("dat1<-ceattle_0$srv_age_obs_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$srv_age_obs_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$srv_age_yrs_",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$srv_age_yrs_",s,sep="")))

		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

		eval(parse(text=paste("dat1.hat<-ceattle_0$srv_age_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$srv_age_hat_",s,sep="")))
		
		years<-1:ceattle_0$nyrs[1,1]+ceattle_0$styr[1,1]-1
		cc<-match(dat1.yr,years)
		nn2<-setdiff(years,dat1.yr)
		nages<-length(dat1[1,])
		coll1<-crain(nages)
		maxx<-1.01*max(c(dat1,dat1.hat,dat2.hat),na.rm=T)
		xx<-seq(0,nages+1,(nages+2)/nages)*((nages+2)/nages)
		xx<-((((nages+3)/nages))*(1:nages))-((((nages+3)/nages))/2)
		for(y in 1:length(years)) {
			if(years[y]%in%dat1.yr[1,]){
			  yy<-which(dat1.yr[1,]%in%years[y])
			  
			  barplot(dat1[yy,],col=coll1,ylim=c(0,maxx),las=2,axes=F)
			 
			  lines(xx,dat2.hat[yy,],col="gray")
			  points(xx,dat1.hat[yy,],pch=16,type="b")
			}
		  if(years[y]%in%dat1.yr[1,]==FALSE) barplot(dat1[1,],col=FALSE,ylim=c(0,maxx),las=2,border=FALSE,axes=F)
			mtext(side=3,years[y],cex=.8,line=-2,las=1,adj=.9)
			coll1<-c(coll1[nages],coll1[-nages])
			
		}
	
		# plot(apply(dat2,2,mean.na),type="p",axes=FALSE,col=coll[2],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		# axis(2,las=2)
		# yatt<-pretty(seq(0,ylimm[s],1e4),4)
		#axis(2,at=yatt,lab=yatt/1e6,las=2)
	
    mtext(side=3,outer=T,spnames[s],adj=0,line=-1)
		# if(ceattle_0$fsh_age_type[s]==1)
		# 	mtext("Age",side=1,line=0,outer=FALSE,font=2)
		# if(ceattle_0$fsh_age_type[s]==2)
		# 	mtext("Length",side=1,line=0,outer=FALSE,font=2)
		# if(s==2)
		# 	mtext("Fishery age or length composition",side=4,line=1,outer=FALSE,font=2)
		# 
		# text(1,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		# axis(1,at=pretty(-10:length(dat1[1,]),10),font=1)

		# if(s==1)
		# 	legend(7,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
}


plot_fshAge<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2)){
	quartz(height=6,width=4)
	par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	layout(c(1,2,3))
	layout.show(3)
	xatt<-seq(0,25,5)
	ylimm<-c(.4,.4,.4)
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")

	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$fsh_age_obs_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$fsh_age_obs_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_fs_h_comp",s,sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_fs_h_comp",s,sep="")))

		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

		eval(parse(text=paste("dat1.hat<-ceattle_0$fsh_age_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$fsh_age_hat_",s,sep="")))

		plot(apply(dat2,2,mean.na),type="p",axes=FALSE,col=coll[2],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		axis(2,las=2)
		yatt<-pretty(seq(0,ylimm[s],1e4),4)
		#axis(2,at=yatt,lab=yatt/1e6,las=2)

		lines(apply(dat1.hat,2,mean.na),col=coll[1],lwd=lwdd[1],lty=ltyy[1],cex=1.2)
		lines(apply(dat2.hat,2,mean.na),col=coll[2],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
		if(ceattle_0$fsh_age_type[s]==1)
			mtext("Age",side=1,line=1.5,outer=FALSE,font=2)
		if(ceattle_0$fsh_age_type[s]==2)
			mtext("Length",side=1,line=1.5,outer=FALSE,font=2)
		if(s==2)
			mtext("Fishery age or length composition",side=2,line=3,outer=FALSE,font=2)
		text(1,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		axis(1,at=pretty(-10:length(dat1[1,]),10),font=1)

		if(s==1)
			legend(7,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
	}	
}
plot_fshAge_byYr<-function(ltyy=c(1,2),coll=collrange1(3),lwdd=c(2,2)){
	quartz(height=6,width=8)
	par(mar=c(0,0,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(.5,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	par(mfrow=c(ceattle_0$nyrs[1],nspp))
	layout(cbind(
		c(1:ceattle_0$nyrs[1]),
		c(1:ceattle_0$nyrs[1]+ceattle_0$nyrs[1]),
		c(1:ceattle_0$nyrs[1]+2*ceattle_0$nyrs[1])
		))
	# layout.show(3)
	xatt<-seq(0,25,5)
	ylimm<-c(.4,.4,.4)
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
	

	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$fsh_age_obs_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$fsh_age_obs_",s,sep="")))
		eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_fsh_comp_",s,"[1,]",sep="")))
		eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_fsh_comp_",s,"[1,]",sep="")))

		eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
		eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))

		eval(parse(text=paste("dat1.hat<-ceattle_0$fsh_age_hat_",s,sep="")))
		eval(parse(text=paste("dat2.hat<-ceattle_2$fsh_age_hat_",s,sep="")))
		years<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
		cc<-match(dat1.yr,years)
		nn2<-setdiff(years,dat1.yr)
		for(y in 1:length(cc)) {
			maxx<-1.01*max(c(dat1[y,],dat1.hat[y,],dat2.hat[y,]),na.rm=T)
			plot(dat1[y,],pch=16,axes=F,type="b",col=coll[1],ylim=c(0,maxx))
			lines(dat1.hat[cc,][y,],col=coll[2])
			lines(dat2.hat[cc,][y,],col=coll[3])
			if(s==1) mtext(side=2,dat1.yr[y],cex=.8,line=-.5,las=2)
		}
		if(length(nn2)>0) for(y in nn2) plot(1,1,col=F,axes=F)
		# plot(apply(dat2,2,mean.na),type="p",axes=FALSE,col=coll[2],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		# axis(2,las=2)
		# yatt<-pretty(seq(0,ylimm[s],1e4),4)
		#axis(2,at=yatt,lab=yatt/1e6,las=2)
	

		if(ceattle_0$fsh_age_type[s]==1)
			mtext("Age",side=1,line=0,outer=FALSE,font=2)
		if(ceattle_0$fsh_age_type[s]==2)
			mtext("Length",side=1,line=0,outer=FALSE,font=2)
		if(s==2)
			mtext("Fishery age or length composition",side=3,line=1,outer=FALSE,font=2)
	}

}


plot_fshAge_byYr2<-function(ltyy=c(1,2),coll=collrange1(3),lwdd=c(2,2),s=1){
  quartz(height=6,width=8)
  par(mar=c(0,0,0,0)) # margins of graph: (bottom,left, top, right)
  par(mgp=c(.5,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
  par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
  par(mfrow=c(ceattle_0$nyrs[1],nspp))
  nyrs<-ceattle_0$nyrs[1]
  rowss<-ceiling(nyrs/3)
  
  layout(cbind(seq(1,nyrs,3),
               seq(1,nyrs,3)+1,
               seq(1,nyrs,3)+2))
  
  layout(cbind(
    (1:rowss),
    (rowss+1):(2*rowss),
    ((rowss*2)+1):(3*rowss)))
  layout.show(nyrs)
  xatt<-seq(0,25,5)
  ylimm<-c(.4,.4,.4)
  spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")
  
  
  
  eval(parse(text=paste("dat1<-ceattle_0$fsh_age_obs_",s,sep="")))
  eval(parse(text=paste("dat2<-ceattle_2$fsh_age_obs_",s,sep="")))
  eval(parse(text=paste("dat1.yr<-ceattle_0$yrs_fsh_comp_",s,"[1,]",sep="")))
  eval(parse(text=paste("dat2.yr<-ceattle_2$yrs_fsh_comp_",s,"[1,]",sep="")))
  
  eval(parse(text=paste("dat1.se<-ceattle_0$srv_bio_se_",s,sep="")))
  eval(parse(text=paste("dat2.se<-ceattle_2$srv_bio_se_",s,sep="")))
  
  eval(parse(text=paste("dat1.hat<-ceattle_0$fsh_age_hat_",s,sep="")))
  eval(parse(text=paste("dat2.hat<-ceattle_2$fsh_age_hat_",s,sep="")))
  
  years<-1:ceattle_0$nyrs[1,1]+ceattle_0$styr[1,1]-1
  cc<-match(dat1.yr,years)
  nn2<-setdiff(years,dat1.yr)
  nages<-length(dat1[1,])
  coll1<-crain(nages)
  maxx<-1.01*max(c(dat1,dat1.hat,dat2.hat),na.rm=T)
  xx<-seq(0,nages+1,(nages+2)/nages)*((nages+2)/nages)
  xx<-((((nages+3)/nages))*(1:nages))-((((nages+3)/nages))/2)
  for(y in 1:length(years)) {
    if(years[y]%in%dat1.yr){
      yy<-which(dat1.yr%in%years[y])
      
      barplot(dat1[yy,],col=coll1,ylim=c(0,maxx),las=2,axes=F)
      
      lines(xx,dat2.hat[yy,],col="gray")
      points(xx,dat1.hat[yy,],pch=16,type="b")
    }
    if(years[y]%in%dat1.yr==FALSE) barplot(dat1[1,],col=FALSE,ylim=c(0,maxx),las=2,border=FALSE,axes=F)
    mtext(side=3,years[y],cex=.8,line=-2,las=1,adj=.9)
    coll1<-c(coll1[nages],coll1[-nages])
    
  }
  
  # plot(apply(dat2,2,mean.na),type="p",axes=FALSE,col=coll[2],lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
  # axis(2,las=2)
  # yatt<-pretty(seq(0,ylimm[s],1e4),4)
  #axis(2,at=yatt,lab=yatt/1e6,las=2)
  
  mtext(side=3,outer=T,spnames[s],adj=0,line=-1)
  # if(ceattle_0$fsh_age_type[s]==1)
  # 	mtext("Age",side=1,line=0,outer=FALSE,font=2)
  # if(ceattle_0$fsh_age_type[s]==2)
  # 	mtext("Length",side=1,line=0,outer=FALSE,font=2)
  # if(s==2)
  # 	mtext("Fishery age or length composition",side=4,line=1,outer=FALSE,font=2)
  # 
  # text(1,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
  # axis(1,at=pretty(-10:length(dat1[1,]),10),font=1)
  
  # if(s==1)
  # 	legend(7,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
}


plot_rec<-function(ltyy=c(2,1), asmt_dat=asmt_dat1,asmt=TRUE,coll=collrange1(2),lwdd=c(2,2),err_bar=FALSE,logg=TRUE,ylimm=c(1e12,20e9,10e9)){
	quartz(height=6,width=4)
	par(mar=c(2,1,1,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	par(mfrow=c(nspp,1))

	#spnames<-c("Walley pollock","Pacific cod","Arrowtooth flounder")
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")

	cmult<-c(1.5,1.2,1.5) 
	collrange1<-colorRampPalette(colors()[c(258,54,132)]) # blues
	years<-1:nyrs+ceattle_0$styr[1,1]-1
	if(logg==TRUE)
	ylimm<-log(ylimm)
	layout.show(nspp)
	years<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
	subplot<-function(model,num=1)
	{
		if(err_bar==TRUE){
			for(i in 1:length(Zoop_years))
			{
				x1<-c((1-cc)*years[i],(1+cc)*Zoop_years[i])
				y1<-exp(model$logR_obs.plus[sp,i])
				y2<-exp(model$logR_obs.minus[sp,i])
				if(logg==TRUE){
					y1<-log(y1)
					y2<-log(y1)
				}

				lines(rep(years[i],2),c(y1,y2))
				lines(x1,rep(y1,2))
				lines(x1,rep(y2,2))
			}
		}
		
		# xx<-c(years[-1],years[-1][nn:1])
		
		xx<-c(years,rev(years))
		yy<-c(exp(model$logR_obs.plus[sp,]),rev(exp(model$logR_obs.minus[sp,])))*1000 # N is in thousands of fish
		col2<-colorRampPalette(c(coll[num],colors()[350])) # blues
		if(logg==TRUE){
			polygon(xx,log(yy),col=col2(10)[9],border=FALSE)
		}else{
			polygon(xx,yy,col=col2(10)[9],border=FALSE)
		}	
		# xx<-c(years[-1],years[-1][nn:1])
		
		yy<-c(exp(model$logR_obs.5.plus[sp,]),rev(exp(model$logR_obs.5.minus[sp,])))*1000 # N is in thousands of fish
		
		if(logg==TRUE){
			polygon(xx,log(yy),col=col2(10)[7],border=FALSE)
			lines(years,log(exp(model$logR_obs[sp,])*1000),pch=16,col=coll[num],lwd=lwdd[num],lty=ltyy[num])
			#points(log(obs)~years,pch=16,type="b",lwd=1,cex=1.2)
		
		}else{
			polygon(xx,yy,col=col2(10)[7],border=FALSE)
			lines(years,exp(model$logR_obs[sp,])*1000,pch=16,col=coll[num],lwd=lwdd[num],lty=ltyy[num])
			#points(obs~years,pch=16,type="b",lwd=1,cex=1.2)
		
		}	
		
		#lines(Zoop_years[-1],exp(mn_results$logR_hat[sp,-1]),pch=16,col=collrange1(3)[2],lwd=2,lty=2)

	}	

	for(sp in 1:nspp){
		obs<-exp(ceattle_0$logR_obs[sp,])*1000 # N is in thousands of fish
		nn<-length(years[-1])
		cc<-.1/1000
		if(logg==TRUE){
			plot(log(obs)~years,pch=16,type="p",ylim=c(1,ylimm[sp]),main="",axes=FALSE,lwd=1,cex=1.2,col="white")
		}else{
			plot((obs)~years,pch=16,type="p",ylim=c(0,ylimm[sp]),main="",axes=FALSE,lwd=1,cex=1.2,col="white")
		}
		subplot(model=ceattle_0,num=1)
		subplot(model=ceattle_2,num=2)
		
		axis(1)
		if(logg==TRUE){
			att<-(c(100,1e4,1e6,pretty(seq(0,exp(ylimm[sp]),1e6),4)))
			att<-c(100,1e4,1e6,1e8,1e10,1e12)
			axis(2,at=log(att),lab=(att)/1,las=2)
			axis(2,at=c(-1e15,1e15),lab=c("",""))
			mtext("Log Recruitment ",side=2,line=2,font=1,outer=TRUE)
		}else{
			att<-(pretty(seq(0,ylimm[sp],1e6),4))
			axis(2,at=att,lab=att/1e9,las=2)
			axis(2,at=c(-1e15,1e15),lab=c("",""))
			if(asmt==T)
				lines(asmt_dat[[sp]][,1],asmt_dat[[sp]][,2],type="l",col="red",lwd=2)
		
			mtext("Recruitment (billions) ",side=2,line=1,font=2,outer=TRUE)
		}	
		
		axis(1,at=c(1900,3000),lab=c("",""))
		text(1979,ylimm[sp]*.95,spnames[sp],pos=4,font=2,cex=1.2)
		#mtext(spnames[sp],side=3,line=-1,font=2)
		if(sp==1)
			legend("topright",c("single-species","multi-species"),lty=ltyy,cex=1.2,col=coll,lwd=lwdd,box.lty=0)
	}
	
	mtext("Year ",side=1,line=1,font=2,outer=TRUE)
}

plot_recTnvnoT<-function(ltyy=c(2,1),coll=col1(2),coll2=col2(2),lwdd=c(2,2),err_bar=FALSE,logg=TRUE,ylimm=c(1e12,20e9,10e9)){
	quartz(height=6,width=6)
	par(mar=c(2,1,1,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	layout(cbind(c(1,3,5),c(2,4,6)))
	#spnames<-c("Walley pollock","Pacific cod","Arrowtooth flounder")
	spnames<-c("Walleye pollock","Pacific cod","Arrowtooth flounder")

	cmult<-c(1.5,1.2,1.5) 
	collrange1<-colorRampPalette(colors()[c(258,54,132)]) # blues
	years<-1:ceattle_0$nyrs+1979-1
	if(logg==TRUE)
	ylimm<-log(ylimm)
	layout.show(6)
	subplot<-function(model,num=1,colluse=coll)
	{
		if(err_bar==TRUE){
			for(i in 1:length(Zoop_years))
			{
				x1<-c((1-cc)*years[i],(1+cc)*Zoop_years[i])
				y1<-exp(model$logR_obs.plus[sp,i])
				y2<-exp(model$logR_obs.minus[sp,i])
				if(logg==TRUE){
					y1<-log(y1)
					y2<-log(y1)
				}

				lines(rep(years[i],2),c(y1,y2))
				lines(x1,rep(y1,2))
				lines(x1,rep(y2,2))
			}
		}
		
		xx<-c(years[-1],years[-1][nn:1])
		yy<-c(exp(model$logR_obs.plus[sp,-1]),exp(model$logR_obs.minus[sp,-1])[nn:1])*1000 # N is in thousands of fish
		col2<-colorRampPalette(c(colluse[num],colors()[350])) # blues
		if(logg==TRUE){
			polygon(xx,log(yy),col=col2(10)[9],border=FALSE)
		}else{
			polygon(xx,yy,col=col2(10)[9],border=FALSE)
		}	
		xx<-c(years[-1],years[-1][nn:1])
		yy<-c(exp(model$logR_obs.5.plus[sp,-1]),exp(model$logR_obs.5.minus[sp,-1])[nn:1])*1000 # N is in thousands of fish
		
		if(logg==TRUE){
			polygon(xx,log(yy),col=col2(10)[7],border=FALSE)
			lines(years[-1],log(exp(model$logR_obs[sp,-1])*1000),pch=16,col=colluse[num],lwd=lwdd[num],lty=ltyy[num])
			#points(log(obs)~years,pch=16,type="b",lwd=1,cex=1.2)
		
		}else{
			polygon(xx,yy,col=col2(10)[7],border=FALSE)
			lines(years[-1],exp(model$logR_obs[sp,-1])*1000,pch=16,col=colluse[num],lwd=lwdd[num],lty=ltyy[num])
			#points(obs~years,pch=16,type="b",lwd=1,cex=1.2)
		
		}	
		#lines(Zoop_years[-1],exp(mn_results$logR_hat[sp,-1]),pch=16,col=colluserange1(3)[2],lwd=2,lty=2)

	}	
	subdiff<-function(model,modelnoT,num=1,colluse=coll,ylimm=c(-5,5),new=TRUE)
	{
		if(new==TRUE){
			plot(years,100*((model$logR_obs[sp,]-modelnoT$logR_obs[sp,])/modelnoT$logR_obs[sp,]),ylim=ylimm,type="l",col=colluse[num],lwd=2,axes=FALSE)
		}else{
			lines(years,100*((model$logR_obs[sp,]-modelnoT$logR_obs[sp,])/modelnoT$logR_obs[sp,]),ylim=ylimm,type="l",col=colluse[num],lwd=2)

		}
	}
		
	
	ylimma<-c(2,2,2)
	sp<-1
	for(sp in 1:nspp){
		obs<-exp(ceattle_0$logR_obs[sp,])*1000 # N is in thousands of fish
		nn<-length(years[-1])
		cc<-.1/1000
		if(logg==TRUE){
			plot(log(obs)~years,pch=16,type="p",ylim=c(1,ylimm[sp]),main="",axes=FALSE,lwd=1,cex=1.2,col="white")
		}else{
			plot((obs)~years,pch=16,type="p",ylim=c(0,ylimm[sp]),main="",axes=FALSE,lwd=1,cex=1.2,col="white")
		}
		
		subplot(model=ceattle_0noT,num=1,colluse=coll2)
		subplot(model=ceattle_2noT,num=2,colluse=coll2)
		subplot(model=ceattle_0,num=1)
		subplot(model=ceattle_2,num=2)
		axis(1)
		if(logg==TRUE){
			att<-(c(100,1e4,1e6,pretty(seq(0,exp(ylimm[sp]),1e6),4)))
			att<-c(100,1e4,1e6,1e8,1e10,1e12)
			axis(2,at=log(att),lab=(att)/1,las=2)
			axis(2,at=c(-1e15,1e15),lab=c("",""))
			if(sp==1) mtext("Log Recruitment ",side=3,line=-1,font=2,outer=FALSE)
		}else{
			att<-(pretty(seq(0,ylimm[sp],1e6),4))
			axis(2,at=att,lab=att/1e9,las=2)
			axis(2,at=c(-1e15,1e15),lab=c("",""))
			if(sp==1) mtext("Recruitment (billions) ",side=3,line=-1,font=2,outer=FALSE)
		}	
		mtext(spnames[sp],side=2,line=2,font=2,outer=FALSE,cex=.8)
		axis(1,at=c(1900,3000),lab=c("",""))
		#text(1979,ylimm[sp]*.95,spnames[sp],pos=4,font=2,cex=1.2)
		#mtext(spnames[sp],side=3,line=-1,font=2)
		if(sp==1)
			legend(2000,ylimm[sp],c("SS with T", "MS with T","SS, no T","MS, no T"),lty=c(ltyy,ltyy),cex=1,col=c(coll,coll2),lwd=2,box.lty=0)
		#	legend(1997,.9*ylimm[sp],c("SS with T", "MS with T","SS, no T","MS, no T"),col=c(coll,coll2),lwd=c(2,2,1,1),lty=c(1,1,1,1),box.lty=0,cex=1)

		subdiff(model=ceattle_0,modelnoT=ceattle_0noT,num=1, ylimm=c(-1,1)*ylimma[sp])
		abline(h=0)
		subdiff(model=ceattle_2,modelnoT=ceattle_2noT,num=2,new=F)

		axis(1)
		att<-pretty(seq(-ylimma[sp],ylimma[sp],.05),4)
			axis(2,at=(att),lab=(att),las=2)
			axis(2,at=c(-1e15,1e15),lab=c("",""))
		
		if(sp==1) mtext("% Difference",side=3,line=-1,font=2,outer=FALSE)
		
		axis(1,at=c(1900,3000),lab=c("",""))
		text(1979,ylimm[sp]*.95,spnames[sp],pos=4,font=2,cex=1.2)
		#mtext(spnames[sp],side=3,line=-1,font=2)
		if(sp==1)
			legend(1997,.9*ylimm[sp],c("single-species","multi-species"),lty=ltyy,cex=1.2,col=colluse,lwd=lwdd,box.lty=0)
	


	}
	
	mtext("Year ",side=1,line=1,font=2,outer=TRUE)
}

plot_M2<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2),ylimm=c(2,2,2)){
	quartz(height=6,width=4)
	par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	layout(c(1,2,3))
	layout.show(3)
	
	xatt<-seq(1960,2050,10)
	
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
	years<-1:nyrs+styr-1
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$M2_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$M2_",s,sep="")))
	
		eval(parse(text=paste("dat1_1<-ceattle_0$M1_",s,sep="")))
		eval(parse(text=paste("dat2_1<-ceattle_2$M1_",s,sep="")))

		plot(years,dat2[,1]+dat2_1[1,1],type="p",axes=FALSE,col="white",lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		yatt<-pretty(seq(0,ylimm[s],.1),4)
		axis(2,at=yatt,lab=yatt,las=2)

		lines(years,dat1[,1]+dat1_1[1,1],col=coll[1],lwd=lwdd[1],lty=ltyy[1],cex=1.2)

		
		lines(years,dat2[,1]+dat2_1[1,1],col=coll[2],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
		#lines(years,dat2[,2]+dat2_1[1,2],col=coll[2],lwd=1,lty=ltyy[2],cex=1.2)
		
		text(1979,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==2)
			legend(1997,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
	}
	mtext("Age 1 natural mortality (M1+M2)",side=2,line=1,outer=TRUE,font=2)	
	mtext("Year",side=1,line=0,outer=TRUE  ,font=2)
}

plot_M2withTvsnoT<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2)){
	quartz(height=6,width=4)
	par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	layout(c(1,2,3))
	layout.show(3)
	
	xatt<-seq(1960,2050,10)
	ylimm<-c(3,3,3)
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
	years<-1:nyrs+styr-1
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$M2_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$M2_",s,sep="")))
	
		eval(parse(text=paste("dat1_1<-ceattle_0$M1_",s,sep="")))
		eval(parse(text=paste("dat2_1<-ceattle_2$M1_",s,sep="")))

		eval(parse(text=paste("dat1noT<-ceattle_0noT$M2_",s,sep="")))
		eval(parse(text=paste("dat2noT<-ceattle_2noT$M2_",s,sep="")))
	
		eval(parse(text=paste("dat1_1noT<-ceattle_0noT$M1_",s,sep="")))
		eval(parse(text=paste("dat2_1noT<-ceattle_2noT$M1_",s,sep="")))

		plot(years,dat2[,1]+dat2_1[1,1],type="p",axes=FALSE,col="white",lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		yatt<-pretty(seq(0,ylimm[s],.1),4)
		axis(2,at=yatt,lab=yatt,las=2)

		lines(years,dat1noT[,1]+dat1_1noT[1,1],col=coll[1],lwd=lwdd[1],lty=ltyy[1],cex=1.2)
		lines(years,dat1[,1]+dat1_1[1,1],col=coll[2],lwd=lwdd[1],lty=ltyy[1],cex=1.2)

		lines(years,dat2noT[,1]+dat2_1noT[1,1],col=coll[1],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
		lines(years,dat2[,1]+dat2_1[1,1],col=coll[2],lwd=lwdd[2],lty=ltyy[2],cex=1.2)			
		#lines(years,dat2[,2]+dat2_1[1,2],col=coll[2],lwd=1,lty=ltyy[2],cex=1.2)
		
		text(1979,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==2)
			legend(1993,.95*ylimm[s],c("single-species (i.e., M1)", "multi-species, with T","multi-species, without T"),col=coll[c(2,2,1)],lwd=lwdd[c(2,2,2)],lty=ltyy[c(1,2,2)],box.lty=0,cex=1.2)
	}
	mtext("Age 1 natural mortality (M1+M2)",side=2,line=1,outer=TRUE,font=2)	
	mtext("Year",side=1,line=0,outer=TRUE  ,font=2)
}

plot_M2_all<-function(ltyy=c(1,1),coll_SS="black", coll_MSM=collrange1,lwdd=c(2,2),ylimm=c(2,1,1)){
			quartz(height=6,width=4)
			par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
			par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
			par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
			layout(c(1,2,3))
			layout.show(3)
			
			xatt<-seq(1960,2050,10)
			spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")

			for(s in 1:nspp)
			{
				eval(parse(text=paste("dat1<-ceattle_0$M2_",s,sep="")))
				eval(parse(text=paste("dat2<-ceattle_2$M2_",s,sep="")))
			
				eval(parse(text=paste("dat1_1<-ceattle_0$M1_",s,sep="")))
				eval(parse(text=paste("dat2_1<-ceattle_2$M1_",s,sep="")))
				nyrs<-ceattle_0$nyrs[1,1]
				years<-1:nyrs+ceattle_0$styr[1,1]-1

				plot(years,dat2[,1]+dat2_1[1,1],type="p",axes=FALSE,col="white",lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
				yatt<-pretty(seq(0,ylimm[s],.1),4)
				axis(2,at=yatt,lab=yatt,las=2)
				nages<-dim(dat2)[2]
				
				lines(years,dat2[,1]+dat2_1[1,1],col=coll_MSM(nages)[1],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
				
				for(a in 1:nages){
					lines(years,dat2[,a]+dat2_1[1,1],col=coll_MSM(nages)[a],lwd=1,lty=ltyy[2],cex=1.2)
				}
				lines(years,dat1[,1]+dat1_1[1,1],col=coll_SS,lwd=lwdd[1],lty=ltyy[1],cex=1.2)
				#lines(years,dat2[,2]+dat2_1[1,2],col=coll[2],lwd=1,lty=ltyy[2],cex=1.2)
				
				text(1979,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
				if(s<3)
					axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
				if(s==3)
					axis(1,at=xatt,font=1)
				if(s==2)
					legend("topright",c("single-species", "multi-species"),col=c(coll_SS,coll_MSM(nages)[1]),lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)
			}
			mtext("Natural mortality (M1+M2)",side=2,line=1,outer=TRUE,font=2)	
			mtext("Year",side=1,line=0,outer=TRUE  ,font=2)
}

plot_Mort_q<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2)){
	quartz(height=6,width=4)
	par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	layout(c(1,2,3))
	layout.show(3)
	
	xatt<-seq(1960,2050,10)
	ylimm<-c(2,1,1)
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
	years<-1:nyrs+styr-1
	col1<-colorRampPalette(colors()[c(71,73)])
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$M2_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$M2_",s,sep="")))
	
		eval(parse(text=paste("dat1_1<-ceattle_0$M1_",s,sep="")))
		eval(parse(text=paste("dat2_1<-ceattle_2$M1_",s,sep="")))

		eval(parse(text=paste("dat1_F<-ceattle_0$F_",s,sep="")))
		eval(parse(text=paste("dat2_F<-ceattle_2$F_",s,sep="")))

		apply(dat1_F,2,mean)
		nage<-dim(dat1)[2]

		plot(1:nage,apply(dat2,2,mean)+apply(dat2_F,2,mean)+dat2_1[1,],type="p",axes=FALSE,col="white",
			lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		yatt<-pretty(seq(0,ylimm[s],.1),4)
		axis(2,at=yatt,lab=yatt,las=2)
		xatt<-seq(0,30,5)
		
		#lines(1:nage,dat2_1[1,],lwd=1,lty=1,cex=1.2)
		qq<-matrix(0,nage,5);colnames(qq)<-names(quantile(dat2_F[,1]))
		quantF<-quantM2<-qq
		for(a in 1:nage){
			quantF[a,]<-quantile(dat2_F[,a])
			quantM2[a,]<-quantile(dat2[,a])
		}
		polygon(c(1:nage,nage:1),c(quantF[,1],rev(quantF[,5])),col=makeTransparent(col1(3)[1],100),border=FALSE)
		polygon(c(1:nage,nage:1),c(quantF[,2],rev(quantF[,4])),col=makeTransparent(col1(3)[1],150),border=FALSE)
		lines(1:nage,quantF[,3],lwd=2,col=makeTransparent(col1(3)[2],200))
		#lines(1:nage,apply(dat2_F,2,mean),lwd=2,col=makeTransparent(col1(3)[2],200))

		polygon(c(1:nage,nage:1),c(quantM2[,1],rev(quantM2[,5])),col=makeTransparent(col1(3)[3],100),border=FALSE)
		polygon(c(1:nage,nage:1),c(quantM2[,2],rev(quantM2[,4])),col=makeTransparent(col1(3)[3],150),border=FALSE)
		lines(1:nage,quantM2[,3],lwd=2,col=makeTransparent(col1(3)[3],200))

		lines(1:nage,dat2_1,col=col1(3)[2],lwd=lwdd[2],lty=2,cex=1.2)
		
		#for (y in 1:34)
		#{
		#	lines(1:nage,dat2_F[y,],lwd=2,lty=2,cex=1.2)
		#	lines(1:nage,dat2[y,],lwd=2,lty=1,cex=1.2)
		#}
		
		#lines(years,dat2[,2]+dat2_1[1,2],col=coll[2],lwd=1,lty=ltyy[2],cex=1.2)
		
		text(1979,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==1){
			legend(8,2,c("F", "M2", "M1"),col=makeTransparent(col1(3),100),lwd=c(10,10,0),lty=c(1,1,2),box.lty=0,cex=1.2)
			legend(8,2,c("F", "M2", "M1"),col=makeTransparent(col1(3),150),lwd=c(6,6,0),lty=c(1,1,2),box.lty=0,cex=1.2)
			legend(8,2,c("F", "M2", "M1"),col=col1(3)[c(2,3,2)],lwd=c(2,2,2),lty=c(1,1,2),box.lty=0,cex=1.2)
		}
	}	
	mtext("Mortality",side=2,line=1,outer=TRUE,font=2)	
	mtext("Age",side=1,line=0,outer=TRUE  ,font=2)
}

plot_Mort<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2)){
	quartz(height=6,width=4)
	par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	layout(c(1,2,3))
	layout.show(3)
	
	xatt<-seq(1960,2050,10)
	ylimm<-c(2,1,1)
	spnames<-c("a) Walleye pollock","b) Pacific cod","c) Arrowtooth flounder")
	years<-1:nyrs+styr-1
	for(s in 1:nspp)
	{
		eval(parse(text=paste("dat1<-ceattle_0$M2_",s,sep="")))
		eval(parse(text=paste("dat2<-ceattle_2$M2_",s,sep="")))
	
		eval(parse(text=paste("dat1_1<-ceattle_0$M1_",s,sep="")))
		eval(parse(text=paste("dat2_1<-ceattle_2$M1_",s,sep="")))

		eval(parse(text=paste("dat1_F<-ceattle_0$F_",s,sep="")))
		eval(parse(text=paste("dat2_F<-ceattle_2$F_",s,sep="")))

		apply(dat1_F,2,mean)
		nage<-dim(dat1)[2]

		plot(1:nage,apply(dat2,2,mean)+apply(dat2_F,2,mean)+dat2_1[1,],type="p",axes=FALSE,col="white",
			lwd=1,lty=ltyy[2],pch=19,cex=1.2,ylim=c(0,ylimm[s]))
		yatt<-pretty(seq(0,ylimm[s],.1),4)
		axis(2,at=yatt,lab=yatt,las=2)
		axis(1,at=seq(0,25,5))
		lines(1:nage,dat2_1[1,],lwd=1,lty=1,cex=1.2)
		lines(1:nage,apply(dat2_F,2,mean),lwd=2,lty=2,cex=1.2)
		lines(1:nage,apply(dat2,2,mean),lwd=2,lty=1,cex=1.2)

		lines(1:nage,apply(dat2_F,2,mean)+1.95*apply(dat2_F,2,se.na),lwd=2,lty=2,cex=1.2)
		lines(1:nage,apply(dat2_F,2,mean)-1.95*apply(dat2_F,2,se.na),lwd=2,lty=2,cex=1.2)

		lines(1:nage,apply(dat2,2,mean)+1.95*apply(dat2,2,se.na),lwd=1,lty=1,cex=1.2)
		lines(1:nage,apply(dat2,2,mean)-1.95*apply(dat2,2,se.na),lwd=1,lty=1,cex=1.2)
		
		#for (y in 1:34)
		#{
		#	lines(1:nage,dat2_F[y,],lwd=2,lty=2,cex=1.2)
		#	lines(1:nage,dat2[y,],lwd=2,lty=1,cex=1.2)
		#}
		
	}
		lines(years,dat2[,1]+dat2_1[1,1],col=coll[2],lwd=lwdd[2],lty=ltyy[2],cex=1.2)
		#lines(years,dat2[,2]+dat2_1[1,2],col=coll[2],lwd=1,lty=ltyy[2],cex=1.2)
		
		text(1979,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
		if(s<3)
			axis(1,at=xatt,lab=rep("",length(xatt)),font=1)
		if(s==3)
			axis(1,at=xatt,font=1)
		if(s==2)
			legend(1997,.95*ylimm[s],c("single-species", "multi-species"),col=coll,lwd=lwdd,lty=ltyy,box.lty=0,cex=1.2)

}


plot_F<-function(ltyy=c(1,2),coll=collrange1(2),lwdd=c(2,2),ylimm=c(.5,3,3)){
	quartz(height=4,width=6)
	par(mar=c(2,1,0,0)) # margins of graph: (bottom,left, top, right)
	par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
	par(oma=c(2,3.5,1,1))# outer margins of graph: (bottom,left, top, right)
	layout(1)
	layout.show(1)
	xatt<-seq(1960,2050,10)
	
	spnames<-c("Pollock","P. cod","Arrowtooth")
	allyr<-ceattle_0$styr[1,1]+1:ceattle_0$nyrs[1,1]-1
	since_yr<-c(1982,1982,1982)
	slopes<-matrix(0,2,3)
	rownames(slopes)<-c("ss","msm")
	colnames(slopes)<-spnames
	for(s in 1:nspp)
	{
		eval(parse(text=paste("Fage<-ceattle_0$F_",s,sep="")))
		eval(parse(text=paste("Fsel<-ceattle_0$fsh_sel_",s,sep="")))
		dat1<-rep(0,dim(Fage)[1])
		nages<-dim(Fage)[2]
		for(y in 1:length(dat1))
		 	dat1[y]<-(Fage[y,]/Fsel[1,])[nages]

		eval(parse(text=paste("Fage<-ceattle_2$F_",s,sep="")))
		eval(parse(text=paste("Fsel<-ceattle_2$fsh_sel_",s,sep="")))
		dat2<-rep(0,dim(Fage)[1])
		nages<-dim(Fage)[2]
		for(y in 1:length(dat1))
		 	dat2[y]<-(Fage[y,]/Fsel[1,])[nages] 	
		if(s==1)
			plot(allyr,dat2,type="l",axes=FALSE,col=coll[2],lwd=2,lty=s,cex=1.2,ylim=c(0,ylimm[s]))
		lines(allyr,dat1,col=coll[1],lwd=2,lty=s)
		lines(allyr,dat2,col=coll[2],lwd=2,lty=s)
		
		yatt<-pretty(seq(0,ylimm[s],.1),4)
		axis(2,at=yatt,lab=yatt,las=2)
		nmn<-which((allyr)==since_yr[s])
		nmn2<-length(allyr)
		nn<-0:(length(nmn:nmn2)-1)
		m0<-lm(Frate~year,data=data.frame(Frate=(dat1[nmn:nmn2]),year=nn))
		m2<-lm(Frate~year,data=data.frame(Frate=(dat2[nmn:nmn2]),year=nn))
		slopes[1,s]<-100*coef(m0)[2]/abs(dat1[nmn])
		slopes[2,s]<-100*coef(m2)[2]/abs(dat2[nmn])
		slopes[1,s]<-(100/length(nmn:nmn2))*(dat1[nmn2]-dat1[nmn])/abs(dat1[nmn])
		slopes[2,s]<-(100/length(nmn:nmn2))*(dat2[nmn2]-dat2[nmn])/abs(dat2[nmn])


	}	
	mtext("Fishing mortality (F)",side=2,line=2,outer=FALSE,font=2)
	#text(1981,ylimm[s]*.95,spnames[s],pos=4,font=2,cex=1.2)
	axis(1,at=xatt,font=1)
	legend(1979,1.08*ylimm[1],c("single-species", "multi-species"),col=coll,lwd=2,lty=1,box.lty=0,cex=1)
	legend(1994,1.08*ylimm[1],spnames,col="black",lwd=2,lty=1:3,box.lty=0,cex=1)

	mtext("Year",side=1,line=0,outer=TRUE,font=2)	
	return(slopes)
}

