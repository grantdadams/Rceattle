rm(list=ls())
library(TMB)
setwd("/Users/kkari/GitHub/TMB_models/ceattle_TMB/")

#--------------------------------------------------
# load data and compile model
#--------------------------------------------------
	load("data/bsai_ceattle.RData") # gets "dat"
	setwd("src")
		compile("fsa2.cpp")
		dyn.load(dynlib("fsa2"))
	setwd("../")

#--------------------------------------------------
# first set up switches
#--------------------------------------------------
	update.est<-TRUE		# TRUE: update the estimation? FLASE = run projections
	# multispp<-TRUE   	# run the model with spp interactions?
	filename<-"test"
	nspp<-3 			# how many species in the model?
	nyrs<-45

# 2. fsh sel mode
# 3. srvy sel mode
# 4. rec_mode mode
# 5. estimate M or calc
	# nrep<-4 			# number of times to repeat estimates
#--------------------------------------------------
# set up controls - could read this in from ctl. file
#--------------------------------------------------

#--------------------------------------------------
# preallocate arrays & lists
#--------------------------------------------------
	ration<-totRat<-matrix(0,nspp,nyrs)
	parms<-list()
	M<-M2<-list()
	for(s in 1:nspp)
		M[[s]]<-dat[[1]]$M
	for(s in 1:nspp)
		M2[[s]]<-dat[[1]]$M*0
	pref<-matrix(0,nspp,nspp);colnames(pref)<-paste0("prey ",1:nspp);rownames(pref)<-paste0("pred ",1:nspp);
	diag(pref)<-c(1,.5,.2)  # cannibalism
	pref[1,2:3]<-c(0.1,.01)  # pred 1
	pref[2,c(1,3)]<-c(1,.01)  # pred 2
	pref[3,c(1,2)]<-c(1,.01)  # pred 3

	calc_M2<-function(sp=1,rat=ration,p=pref[,1],MIn=M[[1]]){

		M2<-MIn*0
		totRat<-apply(rat*p,2,sum,na.rm=T)

		M2[1,]<-exp(as.numeric(scale(totRat))*.1)*log(sp/3)*-1
		M2[2,]<-exp(as.numeric(scale(totRat))*.1)*.1*log(sp/3)*-1
		return(M2)
	}


#--------------------------------------------------
# Fit the model in estimation mode:
#--------------------------------------------------
	runModel<-function(multispp=TRUE,nrep=10){
		setwd("src")
		obj<-opt<-rep<-srep<-list()
		for(s in 1:nspp){
			parameters <- list(
				estpar=0,
			  logN1Y=rep(0,nrow(dat[[s]]$catchNo)),
			  logN1A=rep(0,ncol(dat[[s]]$catchNo)-1),
			  logFY=rep(0,ncol(dat[[s]]$catchNo)),
			  logFA=rep(0,nrow(dat[[s]]$catchNo)),
			  logVarLogCatch=c(0,0),
			  logQ=rep(0,nrow(dat[[s]]$Q1)),
			  logVarLogSurvey=0
			)
			obj[[s]] <- MakeADFun(dat[[s]],parameters,DLL="fsa2", map=list(estpar=factor(NA),logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)
			opt[[s]] <- nlminb(obj[[s]]$par, obj[[s]]$fn, obj[[s]]$gr, control=list(iter.max=1000,eval.max=1000))
			rep[[s]] <- sdreport(obj[[s]])
			srep[[s]]<-summary(sdreport(obj[[s]]))
			ssb<-srep[[s]][rownames(srep[[s]])=="ssb",]
			ration[s,]<-srep[[s]][rownames(srep[[s]])=="ration",1]*log(s/.8)

		}



		for(l in 1:nrep){
			for(s in 1:nspp){
				dat[[s]]$M<-M[[s]]
				if(multispp){
					M2[[s]]<-calc_M2(sp=s,rat=ration,p=pref[,s],MIn=M[[s]])
					dat[[s]]$M<-M[[s]]+M2[[s]]
				}
			}

			# now refit the model with M2
			for(s in 1:nspp){
				parameters <- list(
					estpar=0,
					  logN1Y=rep(0,nrow(dat[[s]]$catchNo)),
					  logN1A=rep(0,ncol(dat[[s]]$catchNo)-1),
					  logFY=rep(0,ncol(dat[[s]]$catchNo)),
					  logFA=rep(0,nrow(dat[[s]]$catchNo)),
					  logVarLogCatch=c(0,0),
					  logQ=rep(0,nrow(dat[[s]]$Q1)),
					  logVarLogSurvey=0
				)
				obj[[s]] <- MakeADFun(dat[[s]],parameters,DLL="fsa2", map=list(estpar=factor(NA),logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)
				opt[[s]] <- nlminb(obj[[s]]$par, obj[[s]]$fn, obj[[s]]$gr, control=list(iter.max=1000,eval.max=1000))
				rep[[s]] <- sdreport(obj[[s]])
				srep[[s]]<-summary(sdreport(obj[[s]]))
				ssb<-srep[[s]][rownames(srep[[s]])=="ssb",]
				# ration<-srep[rownames(srep)=="ration",]
				ration[s,]<-.4*ssb[,1]*log(s/.8)
				# parms[[s]]<-list(
				#   logN1Y=as.numeric(srep[[s]][rownames(srep[[s]])=="logN1Y",1]),
				#   logN1A=as.numeric(srep[[s]][rownames(srep[[s]])=="logN1A",1]),
				#   logFY=as.numeric(srep[[s]][rownames(srep[[s]])=="logFY",1]),
				#   logFA=as.numeric(srep[[s]][rownames(srep[[s]])=="logFA",1]),
				#   logVarLogCatch=as.numeric(srep[[s]][rownames(srep[[s]])=="logVarLogCatch",1]),
				#   logQ=as.numeric(srep[[s]][rownames(srep[[s]])=="logQ",1]),
				#   logVarLogSurvey=as.numeric(srep[[s]][rownames(srep[[s]])=="logVarLogSurvey",1])
				# )
			}

		}
		setwd("../")
		return(list(obj=obj,opt=opt,rep=rep,srep=srep,ration=ration,parms=parms,M2=M2))

	}
	if(update.est){
		msm<-runModel(multispp=TRUE)
		ss<-runModel(multispp=FALSE)
		est<-list(msm=msm,ss=ss)
		if(file.exists(file.path("runs",filename))){
			message(paste0(
				"__________________________________________\n",
				'"',filename,'"'," already exists; ",'"',filename,'"'," will be overwritten if you continue\n",
				"__________________________________________\n"))
		}else{
			dir.create(file.path("runs",filename))
		}

		save(est,file=file.path("runs",filename,"est.Rdata"))
	}else{
		load(file=file.path("runs",filename,"est.Rdata"))
		run_proj()
	}
#--------------------------------------------------
# Plot results
#--------------------------------------------------
	source("scripts/plot_est.R")

#--------------------------------------------------
# Project the model as a MSE
#--------------------------------------------------
# all models:
	# fit recruitment model to est:

	# set ABC using harvest control rule:

	# determine ABC and catch:
#--------------------------------------------------
# now specific projections:
#--------------------------------------------------

proj_filename<-"test1"

mcmc<-1;multispp<-TRUE;estIn<-est
proj_mse<-function(MSE=TRUE,mcmc=100,Recmode=0,Hmode=0, SurveyMode=0, nyrs_fut=50,multispp=TRUE,nrep=10,estIn=est){
		setwd("src")
		rand1<-FALSE
		if(mcmc>1) rand1<-TRUE
		if(multispp==TRUE) srep<-estIn$msm$srep
		if(multispp==FALSE) srep<-estIn$ss$srep
		get_val<-function(val="logN1Y",rand=TRUE,sp=1){
			mn<-as.numeric(srep[[sp]][rownames(srep[[sp]])==val,1])
			sd<-as.numeric(srep[[sp]][rownames(srep[[sp]])==val,2])
			out<-mn
			if(rand==TRUE){
				for(i in 1:length(mn)) out[i]<-rnorm(1,mn[i],sd[i])
			}
			return(out)

		}
		parameters <- list(
				estpar=0,
			  logN1Y=rep(0,nrow(dat[[s]]$catchNo)),
			  logN1A=rep(0,ncol(dat[[s]]$catchNo)-1),
			  logFY=rep(0,ncol(dat[[s]]$catchNo)),
			  logFA=rep(0,nrow(dat[[s]]$catchNo)),
			  logVarLogCatch=c(0,0),
			  logQ=rep(0,nrow(dat[[s]]$Q1)),
			  logVarLogSurvey=0
			)

		parms<-list()
		assign<-function(obj="estpar",parIn=parms[[s]],datIn=0){
			if(length(parIn[[obj]])==length(datIn)){
				return(datIn)
			}else{
				message("object and data lengths don't match")
				return(NA)
			}
		}
		for(s in 1:nspp){
			parms[[s]]<-tmppar<-parameters
			parms[[s]]$estpar<-assign(parIn=parms[[s]],obj="estpar",datIn=0)
			parms[[s]]$logN1Y<-assign(parIn=parms[[s]],obj="logN1Y",datIn=get_val(val="logN1Y",rand=rand1,sp=s))
			parms[[s]]$logN1A=assign(parIn=parms[[s]],obj="logN1A",datIn=get_val(val="logN1A",rand=rand1,sp=s))
			parms[[s]]$logFY<-assign(parIn=parms[[s]],obj="logFY",datIn=get_val(val="logFY",rand=rand1,sp=s))
			parms[[s]]$logFA<-assign(parIn=parms[[s]],obj="logFA",datIn=c(get_val(val="logFA",rand=rand1,sp=s),0,0,0))
			parms[[s]]$logVarLogCatch<-assign(parIn=parms[[s]],obj="logVarLogCatch",datIn=get_val(val="logVarLogCatch",rand=rand1,sp=s) )
			parms[[s]]$logQ<-assign(parIn=parms[[s]],obj="logQ",datIn=get_val(val="logQ",rand=rand1,sp=s))
			parms[[s]]$logVarLogSurvey<-assign(parIn=parms[[s]],obj="logVarLogSurvey",datIn=get_val(val="logVarLogSurvey",rand=rand1,sp=s))

		}


			ny<-ncol(dat[[s]]$catchNo)
			nages<-nrow(dat[[s]]$catchNo)


			map <- list(
				estpar=factor(1),
			  logN1Y=factor(rep(NA,nages)),
			  logN1A=factor(rep(NA,ny-1)),
			  logFY=factor(rep(NA,ny)),
			  logFA=factor(rep(NA,nages)),
			  logVarLogCatch=factor(c(NA,NA)),
			  logQ=factor(rep(NA,nrow(dat[[s]]$Q1))),
			  logVarLogSurvey=factor(NA)
			)
			datp<-dat[[s]]
			datp$est<-0

			dat1<-dat[[s]];dat1$est<-1
			# parms[[s]]<-parameters
			srep_proj<-obj_proj<-opt_proj<-rep_proj<-list()
			obj_proj[[s]] <- MakeADFun(datp,parms[[s]],DLL="fsa2", map=map, silent=TRUE)
			opt_proj[[s]] <- nlminb(obj_proj[[s]]$par, obj_proj[[s]]$fn, obj_proj[[s]]$gr, control=list(iter.max=1000,eval.max=1000))
			rep_proj[[s]] <- sdreport(obj_proj[[s]])
			srep_proj[[s]]<-summary(sdreport(obj_proj[[s]]),ignore.parm.uncertainty = TRUE,par.fixed=T)
			ssb<-srep_proj[[s]][rownames(srep_proj[[s]])=="ssb",]
			ssb2<-obj_proj[[s]]$report()$ssb
			plot(ssb[,1],type="l",col="red")
			lines(ssb2)


# STOPPING HERE FOR NOW - need to now make data full length of the projection, etc.
## Need to make catch years, survey years different lengths
## Need to add flex to add multiple surveys
## Need to add CA, CB FT() function for the ration calcs - do there or in R?
## Need to add temp specific weight at AGE
## Need to add M est




		}



		for(l in 1:nrep){
			for(s in 1:nspp){
				dat[[s]]$M<-M[[s]]
				if(multispp){
					M2[[s]]<-calc_M2(sp=s,rat=ration,p=pref[,s],MIn=M[[s]])
					dat[[s]]$M<-M[[s]]+M2[[s]]
				}
			}

			# now refit the model with M2
			for(s in 1:nspp){
				parameters <- list(
					  logN1Y=rep(0,nrow(dat[[s]]$catchNo)),
					  logN1A=rep(0,ncol(dat[[s]]$catchNo)-1),
					  logFY=rep(0,ncol(dat[[s]]$catchNo)),
					  logFA=rep(0,nrow(dat[[s]]$catchNo)),
					  logVarLogCatch=c(0,0),
					  logQ=rep(0,nrow(dat[[s]]$Q1)),
					  logVarLogSurvey=0
				)
				obj[[s]] <- MakeADFun(dat[[s]],parameters,DLL="fsa2", map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)
				opt[[s]] <- nlminb(obj[[s]]$par, obj[[s]]$fn, obj[[s]]$gr, control=list(iter.max=1000,eval.max=1000))
				rep[[s]] <- sdreport(obj[[s]])
				srep[[s]]<-summary(sdreport(obj[[s]]))
				ssb<-srep[[s]][rownames(srep[[s]])=="ssb",]
				# ration<-srep[rownames(srep)=="ration",]
				ration[s,]<-.4*ssb[,1]*log(s/.8)


		}
		setwd("../")
		return(list(obj=obj,opt=opt,rep=rep,srep=srep,ration=ration,parms=parms,M2=M2))



	# MSE : if set to TRUE, the model will be refit for each year of the projection - if set to false, projection is deterministic
	# mcmc = if set to 1, no random draws, if set to >1 random draws to be pulled from parameter estimates

	# generate survey index from availble biomass
	# project rectruits
	# update the data with an additional year
	# first run the model forward a year by re-estimating the


}





