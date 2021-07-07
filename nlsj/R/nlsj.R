nlsj <- function (formula, data = parent.frame(), start, control = nlsj.control(),
            algorithm = "default", weights=NULL, subset, trace = FALSE,
            lower = -Inf, upper = Inf, ...) {
# ?? left out -- FIXME??    na.action, model = FALSE, (masked from nlxb)
# ?? at this stage ONLY treat "default", but will add bounds
# ?? data in .GlobalEnv -- should be OK
##?? Should a lot of this material be in nlsjModel() to build tools for problem??

   if (is.null(algorithm)) algorithm<-"default"
   algchoice <- c("default", "port", "plinear", "marquardt")
   alg <- which(algorithm %in% algchoice)
   switch(alg,
      "default" = if (trace) cat("nlsj: Using default algorithm\n"),
      # rest
      { msg <- paste("Algorithm choice '",alg,"' not yet implemented or does not exist")
        stop(msg) }
   ) # end switch choices

   dnames <- all.vars(formula)[which(all.vars(formula) %in% ls(data))]
   if (length(dnames) < 1) stop("No data found")
   vnames <- all.vars(formula) # all names in the formula
   pnames <- vnames[ - which(vnames %in% dnames)] # the "non-data" names in the formula
   npar <- length(pnames)
   # cat("npar=",npar,"\n") # Proves this works (??put in nlsjModel()?)
   mraw <- length(eval(as.name(dnames[1])))
   # cat("mraw=",mraw,"\n") # Proves this works (raw length)
   if(is.null(weights)) {
       weights<-rep(1.0, mraw) # set all weights to 1.0
   } 
   else if (any(weights < 0.0)) stop("weights must be non-negative")

   ## Process formula to get names of parameters -- differs from nlsr!!
   stopifnot(inherits(formula, "formula") )
   if (length(formula) == 2) {
       residexpr <- formula[[2]]
   } else if (length(formula) == 3) {
       residexpr <- call("-", formula[[3]], formula[[2]])
   } else stop("Unrecognized formula")


   if( ! missing(subset) && ! is.null(subset) ){ 
     # we need to subset, which we do via the weights
     if (! all(is.integer(subset))) stop("subset must have integer entries")
     if ( any(subset < 1) || any(subset > mraw) ) stop("subset entry out of range")
     #?? need to test these possibilities
     weights[- subset]<-0.0 # NOTE: the minus gets the values NOT in the dataset
   }
   cat("nlsj: weights:")
   print(weights)  

# from nls.R ?? 
#    QR.rhs <- qr(.swts * rhs)
#    lin <- qr.coef(QR.rhs, .swts * lhs)
#    resid <- qr.resid(QR.rhs, .swts * lhs)

##??210702 -- sort out missing start. Is is.null() or is.missing() better?
    if (is.null(start)) { # start not specified
       warning("start vector not specified for nlsj")
       start<-0.9+seq(npar)*0.1 # WARNING: very crude??
       names(start)<-pnames # and make sure these are named?? necessary??
       ## ??? put in ways to get at selfstart models 
    }
    else { # we have a start vector
       snames<-names(start) # names in the start vector
       if ((length(snames) != length(pnames)) || (! all.equal(snames, pnames))) {
           stop("Start names differ in number or name from formula parameter names")
       }
       start <- as.numeric(start) # ensure we convert (e.g., if matrix)
       names(start) <- snames ## as.numeric strips names, so this is needed ??
    }

    #### Build the "model" object ####
    m <- nlsjModel(formula, data, start, weights, lower=lower, upper=upper, control=control)

    # ?? Are we ready to solve?
    njac <- 0 # number of jacobians so far
    nres <- 1
    ## ?? need current prm available, prm_old??
    resraw <- m$rjfun(start) # 
    ssmin <- m$deviance() # get the sum of squares (this is weighted)
    while (! (confInfo <- m$conv()) ) { # main loop ##?? conv NOT correct for rel. offset
       # Here we have to choose method based on controls
       # Inner loop over either line search or Marquardt stabilization
       # returns delta, break from "while" if prm unchanged by delta
       # resid is recalculated, as is new deviance, but NOT jacobian
       ##?? How to track if jacobian needs updating?

       J <- m$jacobian() # this calls derivatives of different types
       delta <- m$incr()
       print(delta)
       prm <- m$getPars()
       fac <- 1.0
       ssnew<-ssmin
       while ((ssnew >= ssmin)  && (fac > control$minFactor)) {
           newp <- prm + fac * delta
           fac <- 0.5 * fac # ?? this is fixed in nls(), but we could alter
           eq <- m$parset(newp)
           cat("newprm:"); print(m$getPars())
           if (! eq) {
              # ?? trying to get NEW values here. Does not want to re-evaluate when running!
              # ?? Is it using p1 and p2 rather than prm??
              ssnew <- m$deviance()
              if (trace) cat("fac=",fac,"   ss=",ssnew,"\n")
              if ( ssnew < ssmin) break
           }
           else {
              cat("Parameters unchanged\n")
              break
           }
       } # end inner while
       cat("after while loop over fac, ssnew=",ssnew," fac=",fac,"\n")
       if (is.na(ssnew)) stop("failed to find lower ss -- parameters unchanged") #?? fix
       tmp <- readline("next")
#       if (trace) cat("Here report progress\n")
    } # end outer while


    ## names(prm) <- pnames # Make sure names re-attached. ??Is this needed??
    result <- list(m=m, convInfo=convInfo, control=ctrl)
    ##?? Add call -- need to set up properly somehow. Do we need model.frame?

    class(result) <- "nlsj" ## CAUSES ERRORS ?? Does it?? 190821
    result
}

