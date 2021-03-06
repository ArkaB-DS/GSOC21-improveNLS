options(warn=-1) # remove this later!
		
## get system info
sy<-Sys.info()
cpu<-benchmarkme::get_cpu()
ramtry<-as.numeric( benchmarkme::get_ram() )
if (!anyNA(ramtry)) { 
     ram <- ramtry/(1024^3) 
	} else { 
  	ram<-paste(sum(as.numeric(system('wmic MemoryChip get Capacity',intern=TRUE)),na.rm=T)/1024^3,"GB")
	}
machid<-paste(sy["nodename"],":",sy["user"],"-",sy["sysname"],"-",sy["release"],
              "|",cpu$model_name,"|",ram," RAM", sep='')

## SETTING DIRECTORY
#scr_dir <- dirname(sys.frame(1)$ofile)
#setwd(scr_dir)
setwd("C:\\Users\\ARKAJYOTI\\Desktop\\runner")

## create final spreadsheet
spreadsheet <- data.frame(DateTime="",FileName="",Solver="",
	Algorithm="",Control="",Residuals="",Deviance="",Gradient="",
	Parameters="",Rmat="",Convergence="",Passed="", Tags="",MachID="")

## create error spreadsheet
spreadsheet_error <- data.frame(DateTime="",FileName="",Solver="",WhatFails="",
	Algorithm="",Control="",Message="",MachID="")

## see if files exist
if (!file.exists("nlsDatabase.csv")){
	write.table(spreadsheet,file='nlsDatabase.csv',append=TRUE,
		sep=",",col.names=TRUE,row.names=FALSE)
}

if (!file.exists("nlsErrorLog.csv")){
	write.table(spreadsheet_error,file='nlsErrorLog.csv',append=TRUE,
		sep=",",col.names=TRUE,row.names=FALSE)
}


NLSproblems <- read.table("problems.csv",header=T)
NLSmethods <- read.table("methods.csv",header=T,sep=",")

problemNumber <- 1
errorNumber <- 1

for(i in 1:nrow(NLSproblems)){
	source(paste("test_files\\",NLSproblems$Name[i],sep=""))
        cat("Sourced problem number",i,"\n")
	for(j in 1:nrow(NLSmethods)){
        cat("Using method number",j,"\n")
		errorNLSFlag <- 0
		errorOtherFlag<-0
		#if(NLSmethods[j,1]=="nlsr::nlxb"){
		#	NLSrunline <- "(formula=NLSformula, data=NLSdata, start=NLSstart,
      	#			    control=List)"
		#}else{###ERROR:Error: $ operator is invalid for atomic vectors
		NLSrunline <- "(formula=NLSformula, data=NLSdata, start=NLSstart,
					lower=NLSlower,upper=NLSupper,
					weights=NLSweights, 
					subset=NLSsubset,
					#na.action=
					algorithm=NLSmethods$algorithm[j],
					control=eval(parse(text=NLSmethods$control[j])))"
		## This is the runline with default algorithm when comparing "marquardt" of nlsj
		NLSrunline0 <- "(formula=NLSformula, data=NLSdata, start=NLSstart,
					lower=NLSlower,upper=NLSupper,
					weights=NLSweights, 
					subset=NLSsubset,
					#na.action=
					control=eval(parse(text=NLSmethods$control[j])))"
			#}
		if (NLSmethods$algorithm[j]=="marquardt"){
			checker.nls<-try(output_nls <- eval(parse(text=paste("nls",NLSrunline0))),silent=TRUE) 
			}else{
			checker.nls<-try(output_nls <- eval(parse(text=paste("nls",NLSrunline))),silent=TRUE) 
			}
		if (inherits(checker.nls,"try-error")) {
			errorNLSFlag = errorNLSFlag + 1
			}
		checker.other<-try(output <- eval(parse(text=paste(NLSmethods[j,1],NLSrunline))),silent=TRUE)
		if (inherits(checker.other,"try-error")) {
			errorOtherFlag = errorOtherFlag + 1
			}
		if (errorNLSFlag==1 | errorOtherFlag==1 ){
			errorNumber <-errorNumber + 1
			spreadsheet_error[errorNumber,3]<-NLSmethods[j,1]
			spreadsheet_error[errorNumber,5]<-NLSmethods[j,2]
			spreadsheet_error[errorNumber,6]<-NLSmethods[j,3]
			spreadsheet_error[errorNumber,8]<-machid
			spreadsheet_error[errorNumber,2]<-NLSproblems$Name[i]
			spreadsheet_error[errorNumber,1]<-format(Sys.time(), "%Y-%m-%d %H:%M")
			
			if(errorNLSFlag==1 & errorOtherFlag==0){
				spreadsheet_error[errorNumber,4]<-"NLS fails."
				spreadsheet_error[errorNumber,7]<-attr(checker.nls,"condition")$message

				}else if (errorNLSFlag==0 & errorOtherFlag==1){
				spreadsheet_error[errorNumber,4]<-paste(NLSmethods$solver[j]," fails.",sep="")
				spreadsheet_error[errorNumber,7]<-attr(checker.other,"condition")$message 				
				} else if (errorNLSFlag==1 & errorOtherFlag==1){
				spreadsheet_error[errorNumber,7]<-paste(attr(checker.nls,"condition")$message,
										attr(checker.other,"condition")$message,
										sep="|")
				spreadsheet_error[errorNumber,4]<- "Both fail."
				}
			next
			}	
		#### TESTING 
		
		## SETTING TOLERANCE
		epstol <- sqrt(.Machine$double.eps*100) # Can replace 100 with nls.control()$offset
	      
		## residuals
		checker.resid<-try(Residuals<-all.equal(as.vector(resid(output_nls)),
		    as.vector(resid(output)),
		    tolerance=epstol*(max(abs(c(as.vector(resid(output_nls)),
					as.vector(resid(output)))
					)) + epstol)))
		if(inherits(checker.resid,"try-error")){
			Residuals <- attr(checker.resid,"condition")$message
			}
		if(grepl("Mean relative difference: ",Residuals)){
			Residuals <- gsub("^.*: *","",Residuals)
			}

		#	## fitted
		#	all.equal(as.vector(fitted(output_nls)),
		#		    as.vector(fitted(output)),
		#		    tolerance=epstol*(max(abs(c(as.vector(fitted(output_nls)),
		#					as.vector(fitted(output)))
		#					)) + epstol))
		## deviance
		checker.dev<-try(Deviance<-all.equal(deviance(output_nls),
			 deviance(output)))
		if(inherits(checker.dev,"try-error")){
			Deviance <- attr(checker.dev,"condition")$message
			}
		if(grepl("Mean relative difference: ",Deviance)){
			Deviance <- gsub("^.*: *","",Deviance)
			}
		## gradient
		checker.grad<-try(Gradient<-all.equal( output_nls$m$gradient(),
			 	output$m$gradient()))
		if(inherits(checker.grad,"try-error")){
			Gradient <- attr(checker.grad,"condition")$message
			}
		if(grepl("Mean relative difference: ",Gradient)){
			Gradient <- gsub("Mean relative difference: ","",Gradient)
			}
		## getPars # difference between getAllPars and getPars?
		checker.pars<-try(Parameters<-all.equal( output_nls$m$getPars(),
			  output$m$getPars()))
		if(inherits(checker.pars,"try-error")){
			Parameters <- attr(checker.pars,"condition")$message
			}
		if(grepl("Mean relative difference: ",Parameters)){
			Parameters <- gsub("^.*: *","",Parameters)
			}
		## Rmat
		checker.rmat<-try(Rmat<-all.equal( as.numeric(output_nls$m$Rmat()), #!!!!NOTE THIS
			  as.numeric(output$m$Rmat())))
		if(inherits(checker.rmat,"try-error")){
			Rmat <- attr(checker.rmat,"condition")$message
			}
		if(grepl("Mean relative difference: ",Rmat)){
			Rmat <- gsub("^.*: *","",Rmat)
			}
		#	## predict
		#	all.equal( output_nls$m$predict(),
		#			  output$m$predict())

		## testing convInfo # FAILED
		checker.conv<-try(Convergence<-all.equal(as.numeric(output_nls$convInfo$isConv),
			 as.numeric(ifelse(NLSmethods[j,1]=="minpack.lm::nlsLM",
						output$convInfo$isConv,output$convInfo))))
		if(inherits(checker.conv,"try-error")){
			Convergence <- attr(checker.conv,"condition")$message
			}
		if(grepl("Mean relative difference: ",Convergence)){
			Convergence <- gsub("^.*: *","",Convergence)
			}
		## write in spreadsheet
		spreadsheet[problemNumber,3] <- NLSmethods[j,1]
		spreadsheet[problemNumber,4] <- NLSmethods[j,2]
		spreadsheet[problemNumber,5] <- NLSmethods[j,3]
		
		spreadsheet[problemNumber,14] <- machid
		spreadsheet[problemNumber,2] <- NLSproblems$Name[i]
		spreadsheet[problemNumber,6] <- Residuals
		spreadsheet[problemNumber,7] <- Deviance
		spreadsheet[problemNumber,8] <-  Gradient
		spreadsheet[problemNumber,9] <- Parameters
		spreadsheet[problemNumber,10] <- Rmat
		spreadsheet[problemNumber,11] <- Convergence
		spreadsheet[problemNumber,12] <- ifelse(isTRUE(all.equal(as.numeric(c(Residuals,Deviance,
									Gradient,
									Parameters,Rmat,Convergence)),rep(1,6))),"Passed",
							   ifelse(isTRUE(all.equal(as.numeric(c(Residuals,Deviance,
									Gradient,
									Parameters,Rmat,Convergence)),rep(0,6))),"Failed",
								"Indeterminate"))
		spreadsheet[problemNumber,1] <- format(Sys.time(), "%Y-%m-%d %H:%M")		
		(problemNumber <- problemNumber +1)
		cat("j is>>>>>>",j,"\n\n")
	}
	cat("Successful problem-->",i,"\n")
}

write.table(spreadsheet,file='nlsDatabase.csv',append=TRUE,
		sep=",",col.names=FALSE,row.names=FALSE)
write.table(spreadsheet_error,file='nlsErrorLog.csv',append=TRUE,
		sep=",",col.names=FALSE,row.names=FALSE)

rm(Convergence,cpu,Deviance,epstol,Gradient,i,j,machid,
	NLSdata,NLSformula,NLSlower,NLSmethods,NLSproblems,NLSrunline,
	NLSstart,NLSsubset,NLSupper,NLSweights,output,output_nls,
	Parameters,problemNumber,ram,Residuals,Rmat,spreadsheet,
	sy,scr_dir,spreadsheet_error,ramtry,NLSrunline,errorOtherFlag,
	errorNumber,errorNLSFlag,checker.mat,checker.resid,checker.pars,
	checker.rmat.checker.resid,checker.pars,checker.other,checker.nls,
	checker.grad,checker.dev,checker.conv,checker.rmat,NLSrunline0)