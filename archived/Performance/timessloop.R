# time ss loops
# JN 2021-7-14
rm(list=ls())
fn <- readline("filename to record sumsquares times=")
sink(fn, split=TRUE) # you will still see on screen, but saved in fn
sy<-Sys.info()
tsstr<- format(Sys.time(), "%Y%m%d%H%M")
cpu<-benchmarkme::get_cpu()
ram<-benchmarkme::get_ram()
machid<-paste(sy["nodename"],":",sy["user"],"-",sy["sysname"],"-",sy["release"],
              "|",cpu$model_name,"|",ram,"bytesRAM", sep='')
cat(machid,"\n")
nn <- c(10, 100, 1000, 10000)
lnn<-length(nn)
meth <- c("loopx2", "loopxx", "sumx2", "sumxx", "crossprod")
lm <- length(meth)
tmean<-matrix(NA, nrow=lnn, ncol=lm)
tmax<-matrix(NA, nrow=lnn, ncol=lm)
tmin<-matrix(NA, nrow=lnn, ncol=lm)
tsd <-matrix(NA, nrow=lnn, ncol=lm)
library(microbenchmark)
set.seed(12345)
for (n in 1:lnn){
   vv <- runif(nn[n])
   nval <- nn[n]
   sslx2<-0.0
   time<-microbenchmark(for (j in 1:nval) {sslx2 <- sslx2+vv[j]^2}, unit="ns")$time
   tmean[n, 1]<-mean(time)
   tsd[n, 1]<-sd(time)
   tmin[n, 1]<-min(time)
   tmax[n, 1]<-max(time)
   ## Note how we force units to be consistent at nanoseconds
   sslxx <- 0.0 
   time<-microbenchmark(for (j in 1:nval) {sslxx <- sslxx+vv[j]*vv[j]}, unit="ns")$time
   tmean[n, 2]<-mean(time)
   tsd[n, 2]<-sd(time)
   tmin[n, 2]<-min(time)
   tmax[n, 2]<-max(time)
   time<-microbenchmark(sssx2<-sum(vv^2), unit="ns")$time
   tmean[n, 3]<-mean(time)
   tsd[n, 3]<-sd(time)
   tmin[n, 3]<-min(time)
   tmax[n, 3]<-max(time)
   time<-microbenchmark(sssxx<-sum(vv*vv), unit="ns")$time
   tmean[n, 4]<-mean(time)
   tsd[n, 4]<-sd(time)
   tmin[n, 4]<-min(time)
   tmax[n, 4]<-max(time)
   time<-microbenchmark(sscp<-as.numeric(crossprod(vv)), unit="ns")$time
   tmean[n, 5]<-mean(time)
   tsd[n, 5]<-sd(time)
   tmin[n, 5]<-min(time)
   tmax[n, 5]<-max(time)
}
colnames(tmean)<- meth
rownames(tmean)<-as.character(nn)
tmean
colnames(tsd)<- meth
rownames(tsd)<-as.character(nn)
tsd
colnames(tmin)<- meth
rownames(tmin)<-as.character(nn)
tmin
colnames(tmax)<- meth
rownames(tmax)<-as.character(nn)
tmax
sink()