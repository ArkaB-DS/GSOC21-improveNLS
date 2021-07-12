nlsjModel <- function(form, data, start, wts=NULL, lower=-Inf, upper=Inf, control)
{
##?? Note that there is a lot of duplication with nlsj() if we are not VERY careful
##?? Many items put in nlenv$... Can possibly simplify.

# This function is designed to set up the information needed to estimate the model
# which it defines and for which it provides the tools
# Needed:
# (formula)		formula: form
# (char vec)		parameters names
# (char vec)		data variable names
# (df)			data frame: data
# (num vec)		weights: wts
# (named num vec)	??not yet -- bounds upper and lower
# (named num vec)	starting parameters (or ...)
# (function)		residual function --> resid vector
# (function)		residual plus jacobian ("gradient) --> residj
# (num vec)		resid
#?? (function)		jacobian function (puts matrix J in attribute "gradient" of resid vector)
# (num matrix)		J -- ?? here or in nlsj()
# (list)		controls: control 
# (function)		convtest function --> returns a logical TRUE when we can TERMINATE,
#			Do we want attributes e.g., message, other information?
# (num)			npar (number of parameters) ?? not critical but nice
# (num)			mres (number of residuals) ?? not critical but nice
# (function)		qsol - set up the solution and produce matrix decomp.
#                       ?? this may have several forms eventually
# (function)		incr -- increment the parameters. May include bounds and masks??
# (environment)         nlenv -- do we need this or not?? If so, we need it better documented.


# Example of nlenv for Croucher: ?? make sure to check doc
# "mres"  
# "njac"  
# "npar"  
# "nres"  
# ?? DO NOT WANT THESE IN THIS FORM "p1"    "p2"    
# "prm"   
# "swts"  "wts"   
# "xdata" "ydata"

# nlsModel() gives m as follows. XX means we have it (or XX Name for replacement):
#  "conv" 
#  "deviance"  
# "fitted"   
#  "formula"  
#  "getAllPars" 
#   "getEnv"   
#  "getPars"   
#  "gradient" 
#  "incr"  -- may not work because we no longer just use Gauss Newton ??
##     to compute the delta (?? with or without fac)
#  "lhs"  
#   "predict"  
#  "resid"  -- ?? would we be better with residu and swts separately
#   "Rmat"  -- ?? only makes sense with QR 
#    "setPars"  -- ?? simpler setPars. If pars changed, null resid attr "gradient"
#  "setVarying" 
#  "trace"  -- ?? this is a FUNCTION in nls(). Is it needed?
##?? Note getRHS is NOT exposed, but is used
## ?? Don't yet handle variable in dot args. But no dot args here.
print(control)
# First segment of this function is to check the specification is valid 
# beyond checks already done in nlsj before call to nlsjModel()

    if (is.null(control$derivmeth)) control$derivmeth="default" # for safety

    stopifnot(inherits(form, "formula"))

    if (is.null(data)) {
	data <- environment(form) # this will handle variables in the parent frame
    }
    else if (is.list(data)){
	    data <- list2env(data, parent = environment(form))
        }
        else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")

    if (is.null(wts)) stop("Weights MUST be declared for nlsjModel")

    cat("nlsjModel: ls(data) =")
    print(ls(data))

    pnames<-names(start)
    resjac <- NULL # to start with
    resvec <- NULL 
    nlenv <- new.env(hash = TRUE, parent = environment(form))
    cat("nlenv created. ls(nlenv):")
    ls(nlenv)
    cat("\n")

##?? Duplication from nlsj, but needed there for subset
    dnames <- all.vars(form)[which(all.vars(form) %in% ls(data))]
    if (length(dnames) < 1) stop("No data found")
    vnames <- all.vars(form) # all names in the formula
    pnames <- vnames[ - which(vnames %in% dnames)] # the "non-data" names in the formula
    npar <- length(pnames)

    nlnames <- ls(nlenv)
    cat("nlnames:"); print(nlnames)
    for (v in dnames){
       if (v %in% nlnames) warning(sprintf("Variable %s already in nlenv!", v))
       nlenv[[v]] = data[[v]]
    }
    cat("nlnames again:")
    print(ls(nlenv))
    nlenv <- list2env(as.list(start),nlenv) # make sure we add parameters. Note "as.list"
    # Above line needed so formulas can be evaluated
    nlenv$epstol4 <- (.Machine$double.eps * control$offset)^4 # used for smallsstest
    nlenv$prm <- start # start MUST be defined at this point
    npar <- length(start)
    nlenv$npar <- npar
    nlenv$QRJ <- NA # set NA so we can test
    cat("nlnames again2:")
    print(ls(nlenv))

    # bounds
    if (length(lower) == 1) lower <- rep(lower, npar) # expand to full dimension
    if (length(upper) == 1) upper <- rep(upper, npar)
    # more checks on bounds
    if (length(lower) != npar) stop("Wrong length: lower")
    if (length(upper) != npar) stop("Wrong length: upper")
    if (any(start < lower) || any(start > upper)) 
        stop("Infeasible start")

    bdmsk <- rep(1, npar)  # set all params free for now
    maskidx <- which(lower==upper)
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters ??? do we want to define U and L parms
    nlenv$bdmsk <- bdmsk
    njac <- 0
    nres <- 0
    nlenv$njac <- njac # Count of jacobian evaluations ("iterations")
    nlenv$nres <- nres # Count of residual evaluations

    if (is.null(wts)) wts <- rep(1, mres) # ?? more complicated if subsetting
    mres<-length(wts[which(wts > 0.0)])
    nlenv$mres <- mres
    nlenv$wts <- wts # 
    swts<-sqrt(wts)
    nlenv$swts <- swts ## ?? can we save a line or two by only using swts
    nlenv$control<-control

# By the time we get here, we assume nlsj has checked data and parameter names.
# If nlsjModel called otherwise, be it on head of caller!

# ?? Can we simplify this -- it seems to set up the parameters
    getPars <- function() unlist(get("prm", nlenv)) # Is this sufficient?? Simplest form??

    # oneSidedFormula ?? Should we be more explicit?
    if (length(form) == 2) {
        residexpr <- form[[2]] ##?? Need to make sure this works -- may need to be call
        lhs <- NULL
        rhs <- eval(form[[2L]], envir = nlenv)
    } else if (length(form) == 3) {
           residexpr <- call("-", form[[3]], form[[2]])
           lhs <- eval(form[[2L]], envir = nlenv)
           # set the "promise" for the lhs of the model
           rhs <- eval(form[[3L]], envir = nlenv)
           lnames<-all.vars(form[[2L]]) # Check that lhs is a single variable from data
#           print(lnames)
           ldname <- which(dnames %in% lnames)
#           str(ldname)
#           print(ldname)
           if (length(lnames) != 1L) {
#              cat("lhs has names:")
#              print(lnames)
              warning("lhs has either no named variable or more than one")
           }
           else { if (control$trace) cat("lhs has just the variable ",lnames,"\n")}
       } 
       else stop("Unrecognized formula")

    resfun <- function(prm) { # only computes the residuals (unweighted)
        #?? What about subsetting??
        nres <- nres + 1
        if (is.null(names(prm))) names(prm) <- names(start)
	  localdata <- list2env(as.list(prm), parent = data)
	  eval(residexpr, envir = localdata) 
    }

    if (all(nlsderivchk(residexpr, names(start)))) { # all derivs can be computed
       rjexpr <- deriv(residexpr, names(start)) ##?? will fail on some functions
    } 
    else rjexpr <- NULL

    if (is.null(rjexpr) && (control$derivmeth == "default")) 
        control$derivmeth <- nlsjcontrol()$altderivmeth

    rjfun <- function(prm) { # Computes residuals and jacobian
        if (is.null(names(prm))) names(prm) <- names(start)
        localdata <- list2env(as.list(prm), parent = data)
        nres <- nres + 1
        if (control$derivmeth == "numericDeriv"){ # use numerical derivatives
           val <- numericDeriv(residexpr, names(prm), rho=localdata)
        }
        else if(control$derivmeth == "default"){ # use analytic
	  val <- eval(rjexpr, envir = localdata)
        } 
        njac <- njac + 1
        val
    }

    nlenv$resid <- swts * resfun(nlenv$prm) # ?? is this the best way to do this
    nlenv$dev <- sum(nlenv$resid^2)
    
#    JJ <- attr(rjfun(start),"gradient")  # ?? would like to save this in nlenv for later use
    # ?? for now use QR -- may have other methods later, e.g., svd
#    QR <- qr(swts * JJ) # ensure we get the jacobian JJ
#    qrDim <- min(dim(QR$qr))
    # ?? Does it make sense to do this here -- probably NOT -- put in nlsj()
#    if(QR$rank < qrDim)
#        stop("singular jacobian matrix at initial parameter estimates")

    parset <- function(newPars){ # nlsj new function
       if (! all.equal(names(newPars), pnames)) stop("newPars has wrong names!")
       newPars<-as.numeric(newPars)
       oprm<-as.numeric(nlenv$prm)
       eq <- all( (oprm+control$offset) == (newPars+control$offset) )
       names(newPars)<-pnames # as.numeric strips names
       if (! eq) {
          nlenv$prm <- newPars
       }
       attr(eq,"prm")<-newPars # ?? do we want this
       eq
    } # return eq -- Are parameters the same?

# All the following are defined at START of problem
    cjmsg <- paste("Max. jacobian evaluations (",control$maxiter,") exceeded")
    crmsg <- paste("Max. residual evaluations (",control$maxres,") exceeded")
    csmsg <- paste("Small wtd sumsquares (deviance) <= ",nlenv$epstol4)
    # above could depend on the problem
    comsg <- paste("Relative offset less than ",control$tol)
    cnmsg <- "No parameters for this problem"
    cxmsg <- "Not yet defined"
    cvmsg <- c(cjmsg, crmsg, csmsg, comsg, cnmsg, cxmsg)

    convCrit <- function() { # defaults
#        cat("convCrit: scaleOffset=",control$scaleOffset,"\n")
##?? can put trace in here
        if (control$trace) {
           cat("convCrit: counts are ",nlenv$njac," ",nlenv$nres,"\n")
           #?? other stuff
        }
        cval <- FALSE # Initially NOT converged
        # max jacs
        cj <- (nlenv$njac > control$maxiter)
        # max fns
        cr <- (nlenv$nres > control$resmax)
        # small ss
        cat("nlenv$dev is ",nlenv$dev,"\n")
        cs <- (control$smallsstest && (nlenv$dev <= nlenv$epstol4))
        # roffset < tolerance for relative offset test
        co <- FALSE # ?? for the moment
        ctol <- NA
#        if (! is.na(nlenv$QRJ)) { # QR available for relative offset test
         ## Make it always available via nlsj()
            scoff <- control$scaleOffset
#            cat("scoff now1=",scoff,"\n")
            if (scoff) scoff <- (mres - npar) * scoff^2 # adjust for problem 
#            cat("scoff now2=",scoff,"\n")
	    rr <- qr.qty(nlenv$QRJ, c(nlenv$resid))
            cat("rr:"); print(rr)
            print(rr[1L:npar])
            print(rr[-(1L:npar)])
            ctol <- sqrt( sum(rr[1L:npar]^2) / (scoff + sum(rr[-(1L:npar)]^2)))
            cat("ctol =",ctol,"\n")
            ##            projected resids ss             deviance
            co <- (ctol <= control$tol) # compare relative offset criterion
#        } else cat("No QRJ\n")
        tmp <- readline("continue")
        # other things??
        cn <- (npar < 1) # no parameters
        cx <- FALSE # Anything else to be added
        if(npar == 0) { # define to avoid exceptions
            cval <- TRUE
            ctol <- NA
        }
        cvec <- c(cj, cr, cs, co, cn, cx)
        cat("cvec:",cvec,"\n")
        cmsg <- "Termination msg: "
        for (i in 1:length(cvec)){
            if (cvec[i]) cmsg <- paste(cmsg,cvmsg[i],"&&")
        }
        attr(cval, "cmsg") <- cmsg
        attr(cval, "ctol") <- ctol
        attr(cval, "nres") <- nlenv$nres
        attr(cval, "njac") <- nlenv$njac
        cval
    } # end convCrit()

#--->
#?? needed?
##    on.exit(remove(i, data, parLength, start, temp, m, gr,
##                   marg, dimGrad, qrDim, gradSetArgs))
## must use weighted resid for use with "port" algorithm.

##?? Those functions not yet working marked with #?? Others seem OK
    m <-
	list(resfun = function(prm) resfun(prm), # ??
             resid = function() resid, # ?? weighted
             rjfun = function(prm) rjfun(prm), # ??
	     fitted = function() rhs, # OK
	     formula = function() form, #OK
	     deviance = function() nlenv$dev, #OK weighted
	     lhs = function() lhs, #OK
##	     jacobian = function() jacobian, #??
	     conv = function() convCrit(), # possibly OK
##	     incr = function() qr.coef(QR, resid),  #?? 
# ?? Need to modify this for different iteration approaches besides Gauss Newton
	     getPars = function() getPars(), #OK
	     getEnv = function() nlenv, #OK
             parset = function(newPars) parset(newPars), #??
	     predict = function(newdata = list(), qr = FALSE)
                 eval(form[[3L]], as.list(newdata), nlenv) #OK
	     )
    class(m) <- "nlsModel"
#    cat("Contents of 'nlenv':")
#    ls.str(nlenv)
    m
}
