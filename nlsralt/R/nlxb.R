nlxb <- function(formula, start, trace = FALSE, data=NULL, lower = -Inf,
                 upper = Inf, masked = NULL, weights=NULL, control=list()) {
    # A simplified and hopefully robust alternative to finding
    # the nonlinear least squares minimizer that causes
    # 'formula' to give a minimal residual sum of squares.
    # 
    # nlxb is particularly intended to allow for the
    # resolution of very ill-conditioned or else near
    # zero-residual problems for which the regular nls()
    # function is ill-suited. 
    # 
    # J C Nash 2014-7-16   nashjc _at_ uottawa.ca
    # 
    # formula looks like 'y~b1/(1+b2*exp(-b3*T))' start MUST be
    # a vector where all the elements are named: e.g.,
    # start=c(b1=200, b2=50, b3=0.3) trace -- TRUE for console
    # output data is a data frame containing data for variables
    # used in the formula that are NOT the parameters. This
    # data may also be defined in the parent frame i.e.,
    # 'global' to this function lower is a vector of lower
    # bounds upper is a vector of upper bounds masked is a
    # character vector of names of parameters that are fixed.
    # control is a list of control parameters. These are: ...
    # 
    # This variant uses a qr solution without forming the sum
    # of squares and cross products t(J)%*%J
    # 
# ?? and put in the weights
#    ######### get data from data frame if exists
#    ######### print(str(data))
#    if (!is.null(data)) {
#        for (dfn in names(data)) {
#            cmd <- paste(dfn, "<-data$", dfn, "")
#            eval(parse(text = cmd))
#        }
#    } else stop("'data' must be a list or an environment")
# ensure params in vector
    pnames <- names(start)
    start <- as.numeric(start) # ensure we convert (e.g., if matrix)
    names(start) <- pnames ## as.numeric strips names, so this is needed
    # bounds
    npar <- length(start)  # number of parameters
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
##      maxlamda <- 1e+60) ## dropped 130709 ??why?
##    epstol <- (.Machine$double.eps) * ctrl$offset # ??161018 - not used elsewhere
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
    phiroot <- sqrt(ctrl$phi)
 # Note spelling of lamda -- a throwback to Ag Can 1974 and way to see if folk are copying code.
 # First get all the variable names:
 #    vn <- all.vars(parse(text = formula))
 # ??? need to fix -- why?, what is wrong
    vn <- all.vars(formula)
    # Then see which ones are parameters (get their positions
    # in the set xx
    pnum <- start  # may simplify later??
    pnames <- names(pnum)
    bdmsk <- rep(1, npar)  # set all params free for now
    # ?? put in  lower==upper mask defn
    maskidx <- union(which(lower==upper), which(pnames %in% masked)) # use both methods for masks  
          # NOTE: %in% not == or order gives trouble
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters
    if (trace) { # diagnostic printout
        cat("Finished masks check\n")
        parpos <- match(pnames, vn) # ?? check this is right??
        datvar <- vn[-parpos]  # NOT the parameters
        cat("datvar:")
        print(datvar)
        for (i in 1:length(datvar)) {
            dvn <- datvar[[i]]
            cat("Data variable ", dvn, ":")
            if (is.null(data)) { 
                print(eval(parse(text = dvn)))
            } else {
                print(with(data, eval(parse(text = dvn)))) 
            }
        }
    }
    trjfn<-model2rjfun(formula, pnum, data=data) 
    if (trace) {
       cat("trjfn:\n")
       print(trjfn)
    }
    ## Call the nlfb function here
## ?? problem is getting the data into the tresfn and tjacfn?? How?
## which gets data into the functions
    resfb <- nlfb(start=pnum, resfn=trjfn, jacfn=trjfn, trace=trace, 
            data=data, lower=lower, upper=upper, maskidx=maskidx, 
            weights=weights, control=ctrl)
##	    control=ctrl, ...)
    resfb$formula <- formula # 190805 to add formula
# ?? should there be any ... arguments
    pnum <- as.vector(resfb$coefficients)
    names(pnum) <- pnames # Make sure names re-attached. ??Is this needed??
##    resfb$coefficients <- pnum ## commented 190821
    result <- resfb
##    attr(result, "pkgname") <- "nlsr"
    class(result) <- "nlsr" ## CAUSES ERRORS ?? Does it?? 190821
    result
}
