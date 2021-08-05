options(warn=-1) # remove this later!

## check function for error log
check<-function(){
	if (inherits(checker,"try-error")){
		spreadsheet_error[errorNumber,4]<<-NLSmethods[j,1]
		spreadsheet_error[errorNumber,5]<<-NLSmethods[j,2]
		spreadsheet_error[errorNumber,6]<<-NLSmethods[j,3]
		spreadsheet_error[errorNumber,2]<<-machid
		spreadsheet_error[errorNumber,3]<<-NLSproblems$Name[i]
		spreadsheet_error[errorNumber,1]<<-format(Sys.time(), "%Y-%m-%d %H:%M")
		spreadsheet_error[errorNumber,7]<<-attr(checker,"condition")$message
		TRUE		
	}else	{ FALSE}
}
		
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
setwd("C:\\Users\\ARKAJYOTI\\Desktop\\runner_v2")

## create final spreadsheet
spreadsheet <- data.frame(DateTime="",MachID="",FileName="",Solver="",
	Algorithm="",Control="",Residuals="",Deviance="",Gradient="",
	Parameters="",Rmat="",Convergence="",Passed="", Tags="" 
	)

## create error spreadsheet
spreadsheet_error <- data.frame(DateTime="",MachID="",FileName="",Solver="",
	Algorithm="",Control="",Message="")

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
	for(j in 1:nrow(NLSmethods)){
		#if(NLSmethods[j,1]=="nlsr::nlxb"){
		#	NLSrunline <- "(formula=NLSformula, data=NLSdata, start=NLSstart,
      	#			    control=List)"
		#}else{###ERROR:Error: $ operator is invalid for atomic vectors
		NLSrunline <- "(formula=NLSformula, data=NLSdata, start=NLSstart,
					algorithm=NLSmethods$algorithm[j],
					control=eval(parse(text=NLSmethods$control[j])))"
		#NLSrunline0 <- "(formula=NLSformula, data=NLSdata, start=NLSstart,
		#			algorithm=NLSmethods$algorithm[j],
		#			control=eval(parse(text=NLSmethods$control[j])))"
			#}
		checker<-try(output_nls <- eval(parse(text=paste("nls",NLSrunline))),silent=TRUE) 
		if (check()) {
			errorNumber <-errorNumber + 1
			next
			}
		checker<-try(output <- eval(parse(text=paste(NLSmethods[j,1],NLSrunline))),silent=TRUE)
		if (check()) {
			errorNumber <-errorNumber + 1
			next
			}

		#### TESTING 
		
		## SETTING TOLERANCE
		epstol <- sqrt(.Machine$double.eps*100) # Can replace 100 with nls.control()$offset
	      
		## residuals
		Residuals<-all.equal(as.vector(resid(output_nls)),
		    as.vector(resid(output)),
		    tolerance=epstol*(max(abs(c(as.vector(resid(output_nls)),
					as.vector(resid(output)))
					)) + epstol))
		#	## fitted
		#	all.equal(as.vector(fitted(output_nls)),
		#		    as.vector(fitted(output)),
		#		    tolerance=epstol*(max(abs(c(as.vector(fitted(output_nls)),
		#					as.vector(fitted(output)))
		#					)) + epstol))
		## deviance
		Deviance<-all.equal(deviance(output_nls),
			 deviance(output))
		## gradient
		Gradient<-all.equal( output_nls$m$gradient(),
			 	output$m$gradient())
		## getPars # difference between getAllPars and getPars?
		Parameters<-all.equal( output_nls$m$getPars(),
			  output$m$getPars())
		## Rmat
		Rmat<-all.equal( as.numeric(output_nls$m$Rmat()), #!!!!NOTE THIS
			  as.numeric(output$m$Rmat()))
		#	## predict
		#	all.equal( output_nls$m$predict(),
		#			  output$m$predict())

		## testing convInfo # FAILED
		Convergence<-all.equal(as.numeric(output_nls$convInfo$isConv),
			 as.numeric(ifelse(NLSmethods[j,1]=="minpack.lm::nlsLM",
						output$convInfo$isConv,output$convInfo)))
		
		## write in spreadsheet
		spreadsheet[problemNumber,4] <- NLSmethods[j,1]
		spreadsheet[problemNumber,5] <- NLSmethods[j,2]
		spreadsheet[problemNumber,6] <- NLSmethods[j,3]
		
		spreadsheet[problemNumber,2] <- machid
		spreadsheet[problemNumber,3] <- NLSproblems$Name[i]
		spreadsheet[problemNumber,7] <- Residuals
		spreadsheet[problemNumber,8] <- Deviance
		spreadsheet[problemNumber,9] <-  Gradient
		spreadsheet[problemNumber,10] <- Parameters
		spreadsheet[problemNumber,11] <- Rmat
		spreadsheet[problemNumber,12] <- Convergence
		spreadsheet[problemNumber,13] <- ifelse(isTRUE(all.equal(as.numeric(c(Residuals,Deviance,
									Gradient,
									Parameters,Rmat,Convergence)),rep(1,6))),"Passed",
							   ifelse(isTRUE(all.equal(as.numeric(c(Residuals,Deviance,
									#radient,
									Parameters,Rmat,Convergence)),rep(0,6))),"Failed",
								"Indeterminate"))
		spreadsheet[problemNumber,1] <- format(Sys.time(), "%Y-%m-%d %H:%M")		
		(problemNumber <- problemNumber +1)
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
	Parameters,problemNumber,ram,Residuals,Rmat,spreadsheet,sy,scr_dir)