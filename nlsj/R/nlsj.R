nlsj <- function (formula, data = parent.frame(), start, control = nlsj.control(),
            algorithm = "default", weights, trace = FALSE,
            lower = -Inf, upper = Inf, ...) {
# ?? left out -- FIXME??    subset,  na.action, model = FALSE, (masked from nlxb)
# ?? at this stage ONLY treat "default", but will add bounds
# ?? MAY add analytic derivatives

if (is.null(algorithm)) algorithm<-"default"
algchoice <- c("default", "port", "plinear", "marquardt")
alg <- which(algorithm %in% algchoice)
switch(alg,
   "default" = cat("Continuing with default algorithm\n"),
   # rest
    { msg <- paste("Algorithm choice '",alg,"' not yet implemented or does not exist")
     stop(msg) }
)


## First process formula to get names of parameters -- differs from nlsr!!
    stopifnot(inherits(form, "formula"))
    if (length(form) == 2) {
        residexpr <- form[[2]]
    } else if (length(form) == 3) {
        residexpr <- call("-", form[[3]], form[[2]])
    } else stop("Unrecognized formula")

    if (is.null(data)) {
	data <- environment(form) # this will handle variables in the parent frame
        }
    else if (is.list(data))
	data <- list2env(data, parent = environment(form))
    else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")

# Now check that data is of correct structure
    if (is.null(data)) {
	data <- environment(form) # this will handle variables in the parent frame
        }
    else if (is.list(data))
	data <- list2env(data, parent = environment(form))
    else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")

dnames <- colnames(data)
pnames <- all.vars(formula) # all names in the formula
pnames <- pnames[! (pnames %in% dnames)] # the "non-data" names in the formula
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
           cat("Start names:")
           print(snames)
           cat("Formula parameter names:")
           print(pnames)
           stop("Start names differ in number or name from formula parameter names")
       }
    }

# ensure params in vector
    pnames <- names(start)
    start <- as.numeric(start) # ensure we convert (e.g., if matrix)
    names(start) <- pnames ## as.numeric strips names, so this is needed
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


    # controls
    ctrl <- list(watch = FALSE, phi = 1, lamda = 1e-04, offset = 100, 
        laminc = 10, lamdec = 4, femax = 10000, jemax = 5000, rofftest = TRUE, 
        smallsstest = TRUE)
    ncontrol <- names(control)
    nctrl <- names(ctrl)
    for (onename in ncontrol) {
        if (!(onename %in% nctrl)) {
            if (trace) cat("control ", onename, " is not in default set\n")
            stop(onename," is not a control for nlxb")
        }
        ctrl[onename] <- control[onename]
    }
    if (trace) print(ctrl)
#??    phiroot <- sqrt(ctrl$phi) -- needed later


    bdmsk <- rep(1, npar)  # set all params free for now
    maskidx <- which(lower==upper)
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters ??? do we want to define U and L parms

    m <- nlsjModel(formula, data, start, weights, upper=NULL, lower=NULL, control=ctrl)
          # ?? can we return / use nlenv??
    # ?? Are we ready to solve?
    njac <- 0 # number of jacobians so far
    nres <- 1
    ## ?? need current prm available, prm_old??
    resbest <- m$resid() # ?? do we need to specify environment nlenv?
       # ?? may or may not have attr "gradient"
    ssmin <- m$deviance() # get the sum of squares
    while (! m$conv() ) { # main loop
       J <- m$jacobian() # this calls derivatives of different types
       # Here we have to choose method based on controls
       # Inner loop over either line search or Marquardt stabilization
       # returns delta, break from "while" if prm unchanged by delta
       # resid is recalculated, as is new deviance, but NOT jacobian
       ##?? How to track if jacobian needs updating?
       if (trace) cat("Here report progress\n")
    }


    ## names(prm) <- pnames # Make sure names re-attached. ??Is this needed??
    result <- list(m=m, convInfo=convInfo, control=ctrl)

    class(result) <- "nlsj" ## CAUSES ERRORS ?? Does it?? 190821
    result
}

