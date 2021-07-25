## lower and upper in algorithm="port" ?? Won't use port for now 210716, but put in bounds
rm(list=ls())
set.seed(123)
x <- runif(200)
a <- b <- 1; c <- -0.1
y <- a+b*x+c*x^2+rnorm(200, sd=0.05)
## plot(x,y)
## curve(a+b*x+c*x^2, add = TRUE)
df <- data.frame(y=y, x=x)
form<-y ~ a+b*x+c*I(x^2)
# vnames<-all.vars(form)
# cat("vnames:"); print(vnames)
# dnames<-ls(df)
# print(dnames)
getlen <- function(lnames) {
  #   print(lnames)
  nn<-length(lnames)
  #    cat("nn=",nn,"\n")
  ll<-rep(NA,nn)
  for (i in 1:nn) ll[i]=length(get(lnames[i]))
  ll
}
lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)),
                                 error = function(e) -1L)

# rm(dnames)
# lnames<-getlen(vnames)
# print(lnames)
# dnames<-vnames[which(lnames>1)]
# print(dnames)
# pnames<-vnames[which(lnames==1)]
# print(pnames)

cat("with dataframe\n")
withdf<-nlsj(form, start = c(a=1, b=1, c=0.1), data=df)
summary(withdf)
cat("no explicit data\n")
nodf<-nlsj(form, start = c(a=1, b=1, c=0.1))
summary(nodf)
