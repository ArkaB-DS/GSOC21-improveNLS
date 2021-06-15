sy<-Sys.info()
tsstr<- format(Sys.time(), "%Y%m%d%H%M")
cpu<-benchmarkme::get_cpu()
ram<-benchmarkme::get_ram()
machid<-paste(sy["nodename"],":",sy["user"],"-",sy["sysname"],"-",sy["release"],
              "|",cpu$model_name,"|",ram,"bytesRAM", sep='')
cat(machid,"\n")
