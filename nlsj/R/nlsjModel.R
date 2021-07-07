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
#  "trace"  -- ?? this is a FUNCTION in nls(). Is it needed?
##?? Note getRHS is NOT exposed, but is used
## ?? Don't yet handle variable in dot args. But no dot args here.

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

    pnames<-names(start)
    resjac <- NULL # to start with
    resvec <- NULL 
    nlenv <- new.env(hash = TRUE, parent = environment(form))
    cat("nlenv created. ls(nlenv):")
    ls(nlenv)
    dnames<-names(data) # Note NOT colnames, as we now have environment
#    cat("dnames:"); print(dnames)
    nlnames <- deparse(substitute(nlenv)) 
#    cat("nlnames:"); print(nlnames)
    
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
    nlenv$control<-control
    nlenv <- list2env(start,nlenv) # make sure we add parameters
    nlenv$prm <- start # ensure we have some parameters

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
        if (control$derivmeth == "numericDeriv"){ # use numerical derivatives
           val <- numericDeriv(residexpr, names(prm), rho=localdata)
        }
        else if(control$derivmeth == "default"){ # use analytic
	  val <- eval(rjexpr, envir = localdata)
        } 
        val
    }

    resid <- swts * resfun(nlenv$prm) # ?? is this the best way to do this
    dev <- sum(resid^2)
    
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

    convCrit <- function() { # defaults
        cval <- FALSE # Initially NOT converged
        cmsg <- "Termination msg: "
        if(npar == 0) {
            cval <- TRUE
            cmsg <- paste(cmsg, "No parameters for this problem")
            ctol <- NA
        }
        else { # Do not have qr results available -- need to recompute!!
          #?? Need to compute what we can -- may be tricky with different criteria
          # e.g., small SS test.
#??          rr <- qr.qty(QR, c(resid)) # rotated residual vector
          scoff <- control$scaleOffset
          if(scoff) scoff <- (length(resid) - npar) * scoff^2 # adjust for problem 
##          ctol <- sqrt( sum(rr[1L:npar]^2) / (scoff + sum(rr[-(1L:npar)]^2)))
          ctol <- 1e-4 #?? temporary value to get function available for checking
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
	     deviance = function() dev, #OK weighted
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
