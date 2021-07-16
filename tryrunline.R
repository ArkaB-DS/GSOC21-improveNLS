xdata = c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata = c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)
Croucherdata<-data.frame(xdata, ydata)
Croucherform <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
# ?? Do we want an array of starts? Some problems have > 1
Croucherstart<-list(p1=1,p2=0.2)
Crouchersubset<-1:8 # Possibly have subsets -- multiple tests.
# ?? do we want separate files for each test. More storage, but then each test
# has a separate ID and can be run independently
## Example use 
# library(nlsj)
# fitj = nlsjx(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), subset=1:8, trace=TRUE)
# summarise
# summary(fitj)
callstring<-"(Croucherform, start=Croucherstart, data=Croucherdata)"
#?? can we embed this in different calls e.g., 

progs <- c("nlsr::nlxb", "nls", "minpack.lm::nlsLM")

for (j in (1:length(progs))){
  runline <- paste("result<-",progs[j],callstring)
  cat("about to run:",runline,"\n")
  eval(parse(text=runline))
  print(result) #?? need to process!
  tmp <- readline("next?")
}
