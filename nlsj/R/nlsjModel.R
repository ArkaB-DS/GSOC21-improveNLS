nlsjModel <- function(form, data, start, wts=NULL, upper=NULL, lower=NULL, control)
{
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


# Example of nlenv for Croucher:
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
#  "trace"    
##?? Note getRHS is NOT exposed, but is used
## ?? Don't yet handle variable in dot args. But no dot args here.

# First segment of this function is to check the specification is valid 
# beyond checks already done in nlsj before call to nlsjModel()

    stopifnot(inherits(form, "formula"))
    cat("at start:")
    print(str(data))
    if (is.null(data)) {
        cat("in is.null(data)\n")
	data <- environment(form) # this will handle variables in the parent frame
    }
    else if (is.list(data)){
            cat("in is.list(data)\n")
            print(ls(data))
	    data <- list2env(data, parent = environment(form))
            cat("after list2env:"); print(ls(data))
        }
        else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")

    cat("str(data):")
    print(str(data))
    print(ls(data))

    pnames<-names(start)
    resjac <- NULL # to start with
    resvec <- NULL 
    nlenv <- new.env(hash = TRUE, parent = environment(form))
    cat("nlenv created. ls(nlenv):")
    ls(nlenv)
    dnames<-names(data) # Possibly could pull in from nlsj ?? what if data an environment?
    cat("dnames:"); print(dnames)
    nlnames <- deparse(substitute(nlenv)) 
    cat("nlnames:"); print(nlnames)
    
    for (v in dnames){
       if (v %in% nlnames) warning(sprintf("Variable %s already in nlenv!", v))
       nlenv[[v]] = data[[v]]
    }
    nlenv$prm <- start # start MUST be defined at this point
    npar <- length(start)
    nlenv$npar <- npar
    nlenv$njac <- 0 # Count of jacobian evaluations
    nlenv$nres <- 0 # Count of residual evaluations
# ?? number of residuals altered by "subset" which is NOT yet functional !!
# ?? Do we want to copy the data to a working array, or subset on the fly?
    mres <- dim(data)[1] # Get the number of residuals (no subset)
    #?? maybe mres<-length(subset)
    nlenv$mres <- mres
    if (is.null(wts)) wts <- rep(1, mres) # ?? more complicated if subsetting
    nlenv$wts <- wts # 
    swts<-sqrt(wts)
    nlenv$swts <- swts ## ?? can we save a line or two

# By the time we get here, we assume nlsj has checked data and parameter names.
# If nlsjModel called otherwise, be it on head of caller!

# ?? Can we simplify this -- it seems to set up the parameters
    getPars <- function() unlist(get("prm", nlenv)) # Is this sufficient?? Simplest form??

    nlenv <- list2env(start,nlenv) # make sure we add parameters
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
           print(lnames)
           ldname <- which(dnames %in% lnames)
           str(ldname)
           print(ldname)
           if (length(lnames) != 1L) { # ?? why do we get to this bit -- 
              cat("lhs has names:")
              print(lnames)
              stop("lhs has either no named variable or more than one")
           }
           else { cat("lhs has just the variable ",lnames,"\n")}
       } 
       else stop("Unrecognized formula")

    resfun <- function(prm) { # only computes the residuals (unweighted)
        #?? What about subsetting??
        if (is.null(names(prm))) 
	    names(prm) <- names(start)
	  localdata <- list2env(as.list(prm), parent = data)
	  eval(residexpr, envir = localdata) 
    }

    if (all(nlsderivchk(residexpr, names(start)))) { # all derivs can be computed
       rjexpr <- deriv(residexpr, names(start)) ##?? will fail on some functions
    } 
    else rjexpr <- NULL

    rjfun <- function(prm) {
        localdata <- list2env(as.list(prm), parent = data)
        if (is.null(names(prm))) names(prm) <- names(start)
        if (is.null(rjexpr)){ # use numerical derivatives
           val <- numericDeriv(residexpr, names(prm), rho=localdata)
        }
        else
	  val <- eval(rjexpr, envir = localdata) 
        val
    }

    addjac <- function(){# ?? assume resjac in nlenv
        if(is.null(resjac)) resjac<-residu()
        if(is.null(attr(resjac, "gradient"))) { # if not defined, get it
           #?? for the moment, just numericDeriv()
           resjac <- numericDeriv(residexpr, pnames, rho = nlenv) # unweighted
        }
        nlenv$njac <- nlenv$njac+1
        resjac # invisible to avoid double display in calling space
    }

    residu <- function() { # ?? no subset here yet
        nlenv$nres <- nlenv$nres+1
        (lhs - rhs) # this should be fine for returning unweighted resids
    }
    resid <- swts * residu()
    dev <- sum(resid^2)
    
    jacobian <- function() { 
        ## cat("In jacobian, resjac null?:", is.null(attr(resjac,"gradient")),"\n")
        if (is.null(attr(resjac,"gradient"))) { resjac <- addjac()}
        JJ <- attr(resjac,"gradient")
        JJ # invisible to avoid double display in calling space
    }
    JJ <- jacobian() # ?? would like to save this in nlenv for later use
    # ?? for now use QR -- may have other methods later, e.g., svd
    QR <- qr(swts * JJ) # ensure we get the jacobian JJ
    qrDim <- min(dim(QR$qr))
    # ?? Does it make sense to do this here
    if(QR$rank < qrDim)
        stop("singular jacobian matrix at initial parameter estimates")

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

    convCrit <- function() { # defaults
        cval <- FALSE # Initially NOT converged
        cmsg <- "Termination msg: "
        if(npar == 0) {
            cval <- TRUE
            cmsg <- paste(cmsg, "No parameters for this problem")
            ctol <- NA
        }
        else {
          rr <- qr.qty(QR, c(resid)) # rotated residual vector
          scoff <- control$scaleOffset
          if(scoff) scoff <- (length(resid) - npar) * scoff^2 # adjust for problem 
          ctol <- sqrt( sum(rr[1L:npar]^2) / (scoff + sum(rr[-(1L:npar)]^2)))
          cval <- (ctol <= control$tol) # compare relative offset criterion
          if (scoff > 0.0) 
               cmsg <- paste("Check relative offset criterion, scaleOffset = ",scoff)
          else cmsg <- "Check relative offset criterion - default"
          if (nlenv$njac > control$maxiter) cmsg<-paste(cmsg,"\n","Too many jacobians")
	}
        attr(cval, "cmsg") <- cmsg
        attr(cval, "ctol") <- ctol
        attr(cval, "nres") <- nlenv$nres
        attr(cval, "njac") <- nlenv$njac
        cval
    } # end convCrit()

#    ctest <- convCrit()
#    cat("ctest:")
#    print(ctest)

#--->
#?? needed?
##    on.exit(remove(i, data, parLength, start, temp, m, gr,
##                   marg, dimGrad, qrDim, gradSetArgs))
## must use weighted resid for use with "port" algorithm.

##?? Those functions not yet working marked with #?? Others seem OK
    m <-
	list(resid = function() resid, # OK
             residu = function() residu(), # ??
	     fitted = function() rhs, # OK
	     formula = function() form, #OK
	     deviance = function() dev, #OK
	     lhs = function() lhs, #OK
             addjac = function() addjac(), #OK with "invisible"??
	     jacobian = function() jacobian(), #??
	     conv = function() convCrit(), # possibly OK
	     incr = function() qr.coef(QR, resid),  #OK
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
