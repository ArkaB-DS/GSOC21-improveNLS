options(warn=-1) # remove this later!
library(benchmarkme)
sy<-Sys.info()
cpu<-get_cpu()

ramtry<-as.numeric( get_ram() )

if (ramtry > 0) { 
     ram <- ramtry/(1024^3) 
} else { 
  ram<-paste(sum(as.numeric(system('wmic MemoryChip get Capacity',intern=TRUE)),na.rm=T)/1024^3,"GB")
}

machid<-paste(sy["nodename"],":",sy["user"],"-",sy["sysname"],"-",sy["release"],
              "|",cpu$model_name,"|",ram," RAM", sep='')

# SETTING DIRECTORY ## ?? DID NOT WORK FOR JN
##scr_dir <- dirname(sys.frame(1)$ofile)
## scr_dir <- dirname(sys.frame(0)$ofile)
scr_dir <- toupper(basename(normalizePath("./")))
if (scr_dir != "IMPROVENLS") {
     msg <- paste("Do you want to run tests in the directory ",scr_dir,"?")
     tmp <- readline(msg)
     if (substr(toupper(tmp),1,1) != "Y") stop("Wrong directory!")
}
## setwd(scr_dir) # Go to directory and run.
#setwd("C:\\Users\\ARKAJYOTI\\Desktop\\runner_v2")


# create final spreadsheet
spreadsheet <- data.frame(DateTime="",MachID="",FileName="",Solver="",
	Algorithm="",Control="",Residuals="",Deviance="",Gradient="",
	Parameters="",Rmat="",Convergence="",Passed="", Tags="" 
	)

if (!file.exists("nlsDatabase.csv")){
	write.table(spreadsheet,file='nlsDatabase.csv',append=TRUE,
		sep=",",col.names=TRUE,row.names=FALSE)
}


NLSproblems <- read.table("problems.csv",header=T)
NLSmethods <- read.table("methods.csv",header=T,sep=",")

problemNumber <- 1
for(i in 1:nrow(NLSproblems)){
	source(paste("test_files/",NLSproblems$Name[i],sep=""))
	for(j in 1:nrow(NLSmethods)){
		NLSrunline <- "(formula=NLSformula, data=NLSdata, start=NLSstart,
					algorithm=NLSmethods$algorithm[j])"#,
					#control=list(NLSmethods$control[j]))"
		output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
		output <- eval(parse(text=paste(NLSmethods[j,1],NLSrunline))) 
		#### TESTING 
		# SETTING TOLERANCE
		epstol <- sqrt(.Machine$double.eps*100) # Can replace 100 with nls.control()$offset
	      # residuals
		Residuals<-all.equal(as.vector(resid(output_nls)),
		    as.vector(resid(output)),
		    tolerance=epstol*(max(abs(c(as.vector(resid(output_nls)),
					as.vector(resid(output)))
					)) + epstol))
		#	# fitted
		#	all.equal(as.vector(fitted(output_nls)),
		#		    as.vector(fitted(output)),
		#		    tolerance=epstol*(max(abs(c(as.vector(fitted(output_nls)),
		#					as.vector(fitted(output)))
		#					)) + epstol))
		# deviance
		Deviance<-all.equal(deviance(output_nls),
			 deviance(output))
		# gradient
		Gradient<-all.equal( output_nls$m$gradient(),
			  attr(output$m$resid(),"gradient"))
		# getPars # difference between getAllPars and getPars?
		Parameters<-all.equal( output_nls$m$getPars(),
			  output$m$getPars())
		# Rmat
		Rmat<-all.equal( as.numeric(output_nls$m$Rmat()), #!!!!NOTE THIS
			  as.numeric(output$m$Rmat()))
		#	# predict
		#	all.equal( output_nls$m$predict(),
		#			  output$m$predict())

		# testing convInfo # FAILED
		Convergence<-all.equal(as.numeric(output_nls$convInfo$isConv),
			 as.numeric(output$convInfo))
		## write in spreadsheet
		spreadsheet[problemNumber,4] <- NLSmethods[j,1]
		spreadsheet[problemNumber,5] <- NLSmethods[j,2]
		spreadsheet[problemNumber,6] <- NLSmethods[j,3]
		
		spreadsheet[problemNumber,2] <- machid
		spreadsheet[problemNumber,3] <- NLSproblems$Name[i]
		spreadsheet[problemNumber,7] <- Residuals
		spreadsheet[problemNumber,8] <- Deviance
		#spreadsheet[problemNumber,9] <- Gradient
		spreadsheet[problemNumber,10] <- Parameters
		spreadsheet[problemNumber,11] <- Rmat
		spreadsheet[problemNumber,12] <- Convergence
		spreadsheet[problemNumber,13] <- ifelse(isTRUE(all.equal(c(Residuals,Deviance,
									Gradient,Parameters,Rmat,Convergence),rep("TRUE",6))),"Passed",
							   ifelse(isTRUE(all.equal(c(Residuals,Deviance,
									Gradient,Parameters,Rmat,Convergence),rep("FALSE",6))),"Failed",
								"Indeterminate"))
		spreadsheet[problemNumber,1] <- format(Sys.time(), "%Y-%m-%d %H:%M")		
		(problemNumber <- problemNumber +1)
	}
	cat("Successful problem-->",i,"\n")
}
write.table(spreadsheet,file='nlsDatabase.csv',append=TRUE,
		sep=",",col.names=FALSE,row.names=FALSE)
		
rm(Convergence,cpu,Deviance,epstol,Gradient,i,j,machid,
	NLSdata,NLSformula,NLSlower,NLSmethods,NLSproblems,NLSrunline,
	NLSstart,NLSsubset,NLSupper,NLSweights,output,output_nls,
	Parameters,problemNumber,ram,Residuals,Rmat,spreadsheet,sy,scr_dir)
# ?? do we need to remove everything?
