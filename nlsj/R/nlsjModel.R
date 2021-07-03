nlsjModel <- function(form, data, start, wts, upper=NULL, lower=NULL, control)
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
# (num vec)		resid
# (function)		jacobian function (puts matrix J in attribute "gradient" of resid vector)
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
#  "resid"   
#   "Rmat"   
#    "setPars"  -- ?? simpler setPars. If pars changed, null resid attr "gradient"
#  "setVarying" 
#  "trace"    
##?? Note getRHS is NOT exposed, but is used

# First segment of this function is to check the specification is valid 
# beyond checks already done in nlsj before call to nlsjModel()

    stopifnot(inherits(form, "formula"))

    nlenv <- new.env(hash = TRUE, parent = environment(form))
    dnames<-colnames(data) # Possibly could pull in from nlsj
    nlnames <- deparse(substitute(nlenv)) 
    for (v in dnames){
       if (v %in% nlnames) warning(sprintf("Variable %s already in nlenv!", v))
       nlenv[[v]] = data[[v]]
    }
    nlenv$prm <- start # start MUST be defined at this point
    npar <- length(start)
    nlenv$npar <- npar
# ?? number of residuals altered by "subset" which is NOT yet functional
# ?? Do we want to copy the data to a working array, or subset on the fly?
    mres <- dim(data)[1] # Get the number of residuals (no subset)
    nlenv$mres <- mres
    if (is.null(wts)) wts <- rep(1, mres) # ?? more complicated if subsetting
    nlenv$wts <- wts # 
    swts<-sqrt(wts)
    nlenv$swts <- swts ## ?? can we save a line or two

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
           lnames<-all.vars(lhs) # Check that lhs is a single variable from data
           ldname <- which(lnames %in% dnames)
           if (length(ldname) != 1L) {
              cat("lhs has names:")
              print(ldname)
              stop("lhs has either no named variable or more than one")
           }
       } 
       else stop("Unrecognized formula")

    residu <- (lhs - rhs) # this should be fine for returning unweighted resids
    resid <- swts * residu
    dev <- sum(resid^2) 

    addjac <- function(){# ?? assume resjac in nlenv
        if(is.null(resjac)) resjac<-residu()
        if(is.null(attr(resjac, "gradient"))) { # if not defined, get it
           #?? for the moment, just numericDeriv()
           resjac <- numericDeriv(residexpr, pnames, rho = nlenv) # unweighted
        }
        resjac
    }
    jacobian <- function() { 
        if (is.null(attr(resjac,"gradient") { resjac <- addjac()}
        JJ <- attr(resjac,"gradient")
        JJ
    }
    JJ <- jacobian() # ?? would like to save this in nlenv for later use
    # ?? for now use QR -- may have other methods later, e.g., svd
    QR <- qr(swts * JJ) # ensure we get the jacobian JJ

    qrDim <- min(dim(QR$qr))
    if(QR$rank < qrDim)
        stop("singular gradient matrix at initial parameter estimates")

    getPars.varying <- function() unlist(mget(names(ind), nlenv))[useParams]
    setPars.noVarying <- function(newPars)
    {
        internalPars <<- newPars # envir = thisnlenv
        for(i in names(ind))
            nlenv[[i]] <- unname(newPars[ ind[[i]] ])
    }
    setPars.varying <- function(newPars)
    {
        internalPars[useParams] <<- newPars
        for(i in names(ind))
            nlenv[[i]] <- unname(internalPars[ ind[[i]] ])
    }
    setPars <- setPars.noVarying

    if(scaleOffset) scaleOffset <- (length(resid)-npar) * scaleOffset^2
    convCrit <- function() {
        if(npar == 0) return(0)
        rr <- qr.qty(QR, c(resid)) # rotated residual vector
        sqrt( sum(rr[1L:npar]^2) / (scaleOffset + sum(rr[-(1L:npar)]^2)))
    }

    on.exit(remove(i, data, parLength, start, temp, m, gr,
                   marg, dimGrad, qrDim, gradSetArgs))
    ## must use weighted resid for use with "port" algorithm.
    m <-
	list(resid = function() resid,
	     fitted = function() rhs,
	     formula = function() form,
	     deviance = function() dev,
	     lhs = function() lhs,
	     gradient = function() .swts * attr(rhs, "gradient"),
	     conv = function() convCrit(),
	     incr = function() qr.coef(QR, resid),
             ##?? Why all the global (scoping) assignments? To put in .GlobalEnv??
	     setVarying = function(vary = rep_len(TRUE, np)) {
                 np <- length(useParams)
		 useParams <<- useP <-
                     if(is.character(vary)) {
                         temp <- logical(np)
                         temp[unlist(ind[vary])] <- TRUE
                         temp
                     } else if(is.logical(vary) && length(vary) != np)
                         stop("setVarying : 'vary' length must match length of parameters")
                     else
                         vary # envir = thisnlenv
		 gradCall[[length(gradCall) - 1L]] <<- useP
		 if(all(useP)) {
		     setPars <<- setPars.noVarying
		     getPars <<- getPars.noVarying
		     getRHS  <<-  getRHS.noVarying
		     npar    <<- length(useP)
		 } else {
		     setPars <<- setPars.varying
		     getPars <<- getPars.varying
		     getRHS  <<-  getRHS.varying
		     npar    <<- sum(useP)
		 }
	     },
	     setPars = function(newPars) {
		 setPars(newPars)
		 resid <<- .swts * (lhs - (rhs <<- getRHS())) # envir = thisnlenv {2 x}
		 dev   <<- sum(resid^2) # envir = thisnlenv
		 if(length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
		 QR <<- qr(.swts * gr) # envir = thisnlenv
		 (QR$rank < min(dim(QR$qr))) # to catch the singular gradient matrix
	     },
             setparj = function(newPars) {
		#   sP2(newPars)
		#   rhs <- m$getRHS()
		#   lhs <- m$lhs
		   resid <- wts * (lhs() - rhs())
		   dev   <- sum(resid^2) # envir = thisnlenv
		   gr <- NULL # to define in case no attr.
		#   if(length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
		   if(length(gr <- attr(resid, "gradient")) == 1L) gr <- c(gr)
		   QR <- qr(wts * gr) # envir = thisnlenv
		   singular <- (QR$rank < min(dim(QR$qr))) # to catch the singular gradient matrix
		   spobj <- list(singular=singular, QR=QR, gr=gr, dev=dev, resid=resid)
	      },
	     getPars = function() getPars(),
	     getAllPars = function() getPars(),
	     getEnv = function() nlenv,
	     trace = function() {
		 d <- getOption("digits")
		 cat(sprintf("%-*s (%.2e): par = (%s)\n", d+4L+2L*(scaleOffset > 0),
			     formatC(dev, digits=d, flag="#"),
			     convCrit(),
			     paste(vapply(getPars(), format, ""), collapse=" ")))
	     },
	     Rmat = function() qr.R(QR),
	     predict = function(newdata = list(), qr = FALSE)
                 eval(form[[3L]], as.list(newdata), nlenv),
             getRHS = function() getRHS # JN added 20210630 temporarily
	     )
    class(m) <- "nlsModel"
    cat("Contents of 'nlenv':")
    ls.str(nlenv)
    m
}
