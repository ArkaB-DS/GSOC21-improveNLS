options(warn=-1) # remove this later!

sy<-Sys.info()
cpu<-benchmarkme::get_cpu()
ram<-benchmarkme::get_ram()
machid<-paste(sy["nodename"],":",sy["user"],"-",sy["sysname"],"-",sy["release"],
              "|",cpu$model_name,"|",ram,"bytesRAM", sep='')


#set working directory
setwd("C:\\Users\\ARKAJYOTI\\Desktop\\runner")

# give location of your test files
testfilenames <- list.files("C:\\Users\\ARKAJYOTI\\Desktop\\runner\\test_files")

# create final spreadsheet
(spreadsheet <- data.frame(DateTime="",MachID="",FileName="",
Tags="",Passed="" ))

# uncomment below code only once to create the spreadsheet
# write.table(spreadsheet,file='nlsDatabase.csv',append=TRUE,
#	sep=",",col.names=TRUE,row.names=FALSE) 

j=1
for (i in testfilenames){
	passed <- 'F'
	location <- paste("C:\\Users\\ARKAJYOTI\\Desktop\\runner\\test_files",
				"\\",i,sep="")
	source(location)
	cat("Files run: ",j,"\\",length(testfilenames),"\n")
	j=j+1
	print("dataframe matrix--->")
	print(NLSmatrix)
	if(sum(NLSmatrix==TRUE,na.rm=TRUE)==(nrow(NLSmatrix)*ncol(NLSmatrix))-5){
		passed <- 'T'
	}
	spreadsheet[,"MachID"] <- machid
	spreadsheet[,"FileName"] <- i
	spreadsheet[,"Passed"] <- passed	
	spreadsheet[,"DateTime"] <- format(Sys.time(), "%Y-%m-%d %H:%M")
	write.table(spreadsheet,file='nlsDatabase.csv',append=TRUE,
		sep=",",col.names=FALSE,row.names=FALSE)
}		


rm(testfilenames,location,j,i,cpu,machid,NLSmatrix,passed,ram,spreadsheet,sy)
ls()

