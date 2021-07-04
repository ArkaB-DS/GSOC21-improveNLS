nlsj <- function (formula, data = parent.frame(), start, control = nlsj.control(),
            algorithm = "default", weights=NULL, trace = FALSE,
            lower = -Inf, upper = Inf, ...) {
# ?? left out -- FIXME??    subset,  na.action, model = FALSE, (masked from nlxb)
# ?? at this stage ONLY treat "default", but will add bounds
# ?? MAY add analytic derivatives


if (is.null(algorithm)) algorithm<-"default"
algchoice <- c("default", "port", "plinear", "marquardt")
alg <- which(algorithm %in% algchoice)
switch(alg,
   "default" = cat("nlsj: Using default algorithm\n"),
   # rest
    { msg <- paste("Algorithm choice '",alg,"' not yet implemented or does not exist")
     stop(msg) }
)

## First process formula to get names of parameters -- differs from nlsr!!
    stopifnot(inherits(formula, "formula"))
    if (length(formula) == 2) {
        residexpr <- formula[[2]]
    } else if (length(formula) == 3) {
        residexpr <- call("-", formula[[3]], formula[[2]])
    } else stop("Unrecognized formula")

    if (is.null(data)) {
	data <- environment(formula) # this will handle variables in the parent frame
        }
    else if (is.list(data))
	data <- list2env(data, parent = environment(formula))
    else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")

# Now check that data is of correct structure
    if (is.null(data)) {
	data <- environment(formula) # this will handle variables in the parent frame
        }
    else if (is.list(data))
	data <- list2env(data, parent = environment(formula))
    else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")

dnames <- ls(data)
vnames <- all.vars(formula) # all names in the formula
pnames <- vnames[ - which(vnames %in% dnames)] # the "non-data" names in the formula
npar <- length(pnames)

# from nls.R ?? 
#    QR.rhs <- qr(.swts * rhs)
#    lin <- qr.coef(QR.rhs, .swts * lhs)
#    resid <- qr.resid(QR.rhs, .swts * lhs)

##??210702 -- sort out missing start. Is is.null() or is.missing() better?
    if (is.null(start)) { # start not specified
       warning("start vector not specified for nlsj")
       start<-0.9+seq(npar)*0.1 # WARNING: very crude??
       ## ??? put in ways to get at selfstart models 
    }
    else { # we have a start vector
       snames<-names(start) # names in the start vector
       if ((length(snames) != length(pnames)) || (! all.equal(snames, pnames))) {
           stop("Start names differ in number or name from formula parameter names")
       }
    }

# ensure params in vector
    start <- as.numeric(start) # ensure we convert (e.g., if matrix)
    names(start) <- pnames ## as.numeric strips names, so this is needed??
    # bounds
    if (length(lower) == 1) 
        lower <- rep(lower, npar) # expand to full dimension
    if (length(upper) == 1) 
        upper <- rep(upper, npar)
# more tests on bounds
    if (length(lower) != npar) 
        stop("Wrong length: lower")
    if (length(upper) != npar) 
        stop("Wrong length: upper")
    if (any(start < lower) || any(start > upper)) 
        stop("Infeasible start")
    if (trace) {
        cat("formula: ")
        print(formula)
        cat("lower:")
        print(lower)
        cat("upper:")
        print(upper)
    }

    bdmsk <- rep(1, npar)  # set all params free for now
    maskidx <- which(lower==upper)
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters ??? do we want to define U and L parms

    if (is.null(dnames)) stop("No data") # ?? fixup for "no data" problems using resid length
    else mres <- length(get(dnames[1], envir=data))
    cat("mres=",mres,"\n") 

    if (is.null(weights)) weights<-rep(1,mres)

    m <- nlsjModel(formula, data, start, weights, upper=NULL, lower=NULL, control=control)
          # ?? can we return / use nlenv??
    # ?? Are we ready to solve?
    njac <- 0 # number of jacobians so far
    nres <- 1
    ## ?? need current prm available, prm_old??
    resbest <- m$resid() # ?? do we need to specify environment nlenv?
       # ?? may or may not have attr "gradient"
    ssmin <- m$deviance() # get the sum of squares
    while (! (confInfo <- m$conv())  ) { # main loop
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
       while ((m$deviance() >= ssmin)  && (fac > control$minFactor)) {
           newp <- prm + fac * delta
           fac <- 0.5 * fac # ?? this is fixed in nls(), but we could alter
           ssnew <- NA
           eq <- m$parset(newp)
           cat("newprm:")
           print(m$getPars())
            if (! eq) {
              # ?? trying to get NEW values here. Does not want to re-evaluate when running!
              ssnew <- eval(m$deviance, m$getEnv())()
              if (trace) cat("fac=",fac,"   ss=",ssnew,"\n")
              if ( ssnew < ssmin) break
           }
           else {
              cat("Parameters unchanged\n")
              break
           }
       }
       cat("after while loop over fac, ssnew=",ssnew,"\n")
       if (is.na(ssnew)) stop("failed to find lower ss -- parameters unchanged") #?? fix
       tmp <- readline("next")
#       if (trace) cat("Here report progress\n")
    }


    ## names(prm) <- pnames # Make sure names re-attached. ??Is this needed??
    result <- list(m=m, convInfo=convInfo, control=ctrl)

    class(result) <- "nlsj" ## CAUSES ERRORS ?? Does it?? 190821
    result
}

