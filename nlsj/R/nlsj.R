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

    cat("nlsj env:"); print(ls())

    #### Build the "model" object ####
    m <- nlsjModel(formula, data, start, weights, lower=lower, upper=upper, control=control)
    cat("nlsj env after m:"); print(ls())
    nlenv <- m$getEnv() # we can now store things in nlenv
    cat("contents of nlenv:"); print(ls(nlenv))
    tmp <- readline("cont.")
    # ?? Are we ready to solve?
    nlenv$njac <- 0 # number of jacobians so far
    nlenv$nres <- 1
    swts <- nlenv$swts # ?? do we need to do this sort of assignment (i.e, can we aling envs?)
    resraw <- m$rjfun(start) # includes Jacobian
    J <- swts * attr(resraw,"gradient") 
    nlenv$QRJ <-qr(J)
    ssmin <- m$deviance() # get the sum of squares (this is weighted)
    cat("Before iteration, deviance=",ssmin,"\n")
    prm <- m$getPars()
    #?? At this stage, have res and J
    # Do we want to have a QR or SVD available to compute conv. 
    while (! m$conv() ) { #?? conv() consolidates several options (as available) 
       # Here we have to choose method based on controls
       # Inner loop over either line search or Marquardt stabilization
       # returns delta, break from "while" if prm unchanged by delta
       # resid is recalculated, as is new deviance, but NOT jacobian
       ##?? How to track if jacobian needs updating?

       J <- swts * attr(resraw,"gradient") 
       wres <- swts * resraw # get wts residuals
       cat("wres:"); print(wres)
       nlenv$QRJ <- qr(J)
       cat("QRJ:"); print(str(nlenv$QRJ))
       qrDim <- min(dim(nlenv$QRJ$qr))
       if (nlenv$QRJ$rank < qrDim) stop("Singular jacobian") # for now don't continue

       delta <- qr.coef(nlenv$QRJ, -wres) # LS solve of J delta ~= -wres
       cat("delta:")
       print(delta)
       fac <- 1.0
       ssnew<-ssmin
       while ((ssnew >= ssmin)  && (fac > control$minFactor)) {
           newp <- prm + fac * delta
           fac <- 0.5 * fac # ?? this is fixed in nls(), but we could alter
           cat("newp:"); print(as.numeric(newp))
           eq <- all( (prm+control$offset) == (newp+control$offset) )
           if (! eq ) {
              # ?? trying to get NEW values here. Does not want to re-evaluate when running!
              # ?? Is it using p1 and p2 rather than prm??
              newwres <- swts*m$resfun(newp)
              ssnew <- as.numeric(crossprod(newwres))
              if (trace) cat("fac=",fac,"   ss=",ssnew,"\n")
              if ( ssnew < ssmin) break
           }
           else {
              cat("Parameters unchanged\n")
              break
           }
       } # end inner while
       cat("after inner while loop over fac, ssnew=",ssnew," fac=",fac,"\n")
       if (is.na(ssnew)) stop("NA ss -- parameters unchanged") #?? fix
       if (ssnew < ssmin) {
	 prm <- newp
       	 ssmin <- ssnew
         resraw <- m$rjfun(prm) # ?? new res and gradient -- does extra res??
         tmp <- readline("next")
        }
        else tconv <- TRUE
#       if (trace) cat("Here report progress\n")
    } # end outer while
    ## names(prm) <- pnames # Make sure names re-attached. ??Is this needed??
    result <- list(m=m, convInfo=convInfo, control=ctrl)
    ##?? Add call -- need to set up properly somehow. Do we need model.frame?

    class(result) <- "nlsj" ## CAUSES ERRORS ?? Does it?? 190821
    result
}
