gnjn <- function(start, resfn, jacfn = NULL, trace = FALSE, 
		data=NULL, control=list(), ...){
# simplified Gauss Newton
   klimit =4000 # change as needed
   offset = 1e6 # for no change in parms
   stepred <- 0.5 # start with this as per nls()
   par <- start
   cat("starting parameters:")
   print(par)
   res <- resfn(par, data, ...)
   ssbest <- as.numeric(crossprod(res))
   cat("initial ss=",ssbest,"\n")
   par0 <- par
   kres <- 1
   kjac <- 0
   keepon <- TRUE
   while (keepon) {
      cat("kjac=",kjac,"  kres=",kres,"  SSbest now ", ssbest,"\n")
      print(par)
      JJ <- jacfn(par, data, ...)
      kjac <- kjac + 1
      QJ <- qr(JJ)
      delta <- qr.coef(QJ, -res)
      ss <- ssbest + offset*offset # force evaluation
      step <- 1.
      while (ss > ssbest) {
        par <- par0+delta * step
        if (max(par + offset) != max(par0 + offset)){
           res <- resfn(par, data, ...)
           ss <- as.numeric(crossprod(res))
           kres <- kres + 1
           cat("step =", step,"  ss=",ss,"\n")
##           tmp <- readline("continue")
           if (ss > ssbest) {
              step <- step * stepred
           } else {
              par0 <- par
              ssbest <- ss
           }
        } # end inner loop
        if (kjac >= klimit)  { 
            keepon = FALSE
            cat("artificial stop at kjac=",klimit," -- we only want to check output\n") 
            break
        }
      } else { keepon <- FALSE # no change in parameters, no lower ss }
   } # end main iteration (keepon still TRUE)
} # seems to need this

} # end gnjne
