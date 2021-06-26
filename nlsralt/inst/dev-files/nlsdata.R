## @knitr nlsdata.R
# try different ways of supplying data to R nls stuff
ydata <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558,
50.156, 62.948, 75.995, 91.972)
ttdata <- seq_along(ydata) # for testing

mydata <- data.frame(y = ydata, tt = ttdata)

hobsc <- y ~ 100*b1/(1 + 10*b2 * exp(-0.1 * b3 * tt))

ste <- c(b1 = 2, b2 = 5, b3 = 3)

# let's try finding the variables

findmainenv <- function(formula, prm) {
   vn <- all.vars(formula)
   pnames <- names(prm)
   ppos <- match(pnames, vn)
   datvar <- vn[-ppos]
   cat("Data variables:")
   print(datvar)
   cat("Are the variables present in the current working environment?\n")
   for (i in seq_along(datvar)){
       cat(datvar[[i]]," : present=",exists(datvar[[i]]),"\n")
   }
}


findmainenv(hobsc, ste)
y <- ydata
tt <- ttdata
findmainenv(hobsc, ste)
rm(y)
rm(tt)

# ===============================


# let's try finding the variables in dotargs

finddotargs <- function(formula, prm, ...) {
   dots <- list(...)
   cat("dots:")
   print(dots)
   cat("names in dots:")
   dtn <- names(dots)
   print(dtn)
   vn <- all.vars(formula)
   pnames <- names(prm)
   
   ppos <- match(pnames, vn)
   datvar <- vn[-ppos]
   cat("Data variables:")
   print(datvar)
   cat("Are the variables present in the dot args?\n")
   for (i in seq_along(datvar)){
       dname <- datvar[[i]]
       cat(dname," : present=",(dname %in% dtn),"\n")
   }
}

finddotargs(hobsc, ste, y=ydata, tt=ttdata)
# ===============================

y <- ydata
tt <- ttdata
tryq <- try(nlsquiet <- nls(formula=hobsc, start=ste))
if (class(tryq) != "try-error") {print(nlsquiet)} else {cat("try-error\n")}
#- OK
rm(y)
rm(tt)


tdots<-try(nlsdots <- nls(formula=hobsc, start=ste, y=ydata, tt=ttdata))
if ( class(tdots) != "try-error") {print(nlsdots)} else {cat("try-error\n")}
#- Fails
tframe <- try(nlsframe <- nls(formula=hobsc, start=ste, data=mydata))
if (class(tframe) != "try-error") {print(nlsframe)} else {cat("try-error\n")}
#- OK

library(nlsr)
y <- ydata
tt <- ttdata
tquiet <- try(nlsrquiet <- nlxb(formula=hobsc, start=ste))
if ( class(tquiet) != "try-error") {print(nlsrquiet)}
#- OK
rm(y)
rm(tt)
## this will fail
## test <- try(nlsrdots <- nlxb(formula=hobsc, start=ste, y=ydata, tt=ttdata))
## but ...
mydata<-data.frame(y=ydata, tt=ttdata)
test <- try(nlsrdots <- nlxb(formula=hobsc, start=ste, data=mydata))
  if (class(test) != "try-error") { print(nlsrdots) } else {cat("Try error\n") }
#- Note -- does NOT work -- do we need to specify the present env. in nlfb for y, tt??
test2 <- try(nlsframe <- nls(formula=hobsc, start=ste, data=mydata))
if (class(test) != "try-error") {print(nlsframe) } else {cat("Try error\n") }
#- OK


library(minpack.lm)
y <- ydata
tt <- ttdata
nlsLMquiet <- nlsLM(formula=hobsc, start=ste)
print(nlsLMquiet)
#- OK
rm(y)
rm(tt)
## Dotargs
tdots <- try(nlsLMdots <- nlsLM(formula=hobsc, start=ste, data=mydata))
if (class(tdots) != "try-error") { print(nlsLMdots) } else {cat("try-error\n") }
#-  Note -- does NOT work
## dataframe
tframe <- try(nlsLMframe <- nlsLM(formula=hobsc, start=ste, data=mydata) )
if (class(tdots) != "try-error") {print(nlsLMframe)} else {cat("try-error\n") }
#- does not work

## detach("package:nlsr", unload=TRUE)
## Uses nlmrt here for comparison
## library(nlmrt)
## txq <- try( nlxbquiet <- nlxb(formula=hobsc, start=ste))
## if (class(txq) != "try-error") {print(nlxbquiet)} else { cat("try-error\n")}
#- Note -- does NOT work
## txdots <- try( nlxbdots <- nlxb(formula=hobsc, start=ste, y=y, tt=tt) )
## if (class(txdots) != "try-error") {print(nlxbdots)} else {cat("try-error\n")}
#- Note -- does NOT work
## dataframe
## nlxbframe <- nlxb(formula=hobsc, start=ste, data=mydata)
## print(nlxbframe)
#- OK

