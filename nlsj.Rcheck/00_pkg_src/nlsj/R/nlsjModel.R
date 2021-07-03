nlsjModel <- function(form, data, start, wts, upper=NULL, lower=NULL)
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
# (num matrix)		J
# (list)		controls: ctrl 
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
#  "incr" 
#  "lhs"  
#   "predict"  
#  "resid"   
#   "Rmat"   
#    "setPars"   
#  "setVarying" 
#  "trace"    

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

# By the time we get here, we assume nlsj has checked data and parameter names.
# If nlsjModel called otherwise, be it on head of caller!

# ?? Can we simplify this -- it seems to set up the param
    getPars <- function() unlist(get("prm", nlenv)) # Is this sufficient?? Simplest form??
    # oneSidedFormula ?? Should we be more explicit?
    if (length(form) == 2) {
        residexpr <- form[[2]]
        lhs <- NULL
        rhs <- eval(form[[3L]], envir = nlenv)
    } else if (length(form) == 3) {
        residexpr <- call("-", form[[3]], form[[2]])
        lhs <- eval(form[[2L]], envir = nlenv)
        # set the "promise" for the lhs of the model
        rhs <- eval(form[[3L]], envir = nlenv)
    } else stop("Unrecognized formula")


    if(!is.null(upper)) upper <- rep_len(upper, parLength)
    # Expand upper and lower to full vectors
    useParams <- rep_len(TRUE, parLength)
    # JN Seems to be a logical vector to indicate free parameters.
    lhs <- eval(form[[2L]], envir = nlenv)
    # set the "promise" for the lhs of the model
    lnames<-all.vars(lhs)
    if (any(lnames %in% pnames)) stop("lhs has parameters -- reformulate!") ##?? fix?
    rhs <- eval(form[[3L]], envir = nlenv)
    # similarly for the rhs
    #?? WHY THE DOT -- does this not hide the weights??
    ## non-zero is TRUE; note that missing() is a VERY special function 
# .swts are square roots of wts. But why make them hidden.
    .swts <- if(!missing(wts) && length(wts))
        sqrt(wts) else rep_len(1, length(rhs))
    ##JN: the weights, which are put into the nlenv
#    cat("nlenv$.swts:")
#    print(nlenv$.swts)
#    readline("cont.")

    resid <- .swts * (lhs - rhs)
    cat("resid:")
    print(resid)

    dev <- sum(resid^2)
    cat("dev=",dev,"\n")
    if(is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function() {
            if(is.null(upper)) # always for "default"
                numericDeriv(form[[3L]], names(ind), nlenv, central = nDcentral)
            else # possibly with "port"
                numericDeriv(form[[3L]], names(ind), nlenv,
                             dir = ## ifelse(internalPars < upper, 1, -1)
                                 -1 + 2*(internalPars < upper), central = nDcentral)
        }
## ?? above changes from forward to backward diff if on upper bound. 
## assumes that we step away from lower bound
        getRHS <- getRHS.noVarying
        rhs <- getRHS()
        cat("when no gradient, rhs:")
        print(rhs)
    } else {
        getRHS.noVarying <- function() eval(form[[3L]], envir = nlenv)
        getRHS <- getRHS.noVarying
        cat("Have gradient, getRHS:")
        print(getRHS)
    }
    dimGrad <- dim(attr(rhs, "gradient"))
    marg <- length(dimGrad)
    cat("dimGrad, marg:",dimGrad,marg,"\n")
    ##JN marg is number of dimensions of "gradient" matrix. Should be 2??
    if(marg > 0L) {
        gradSetArgs <- vector("list", marg + 1L)
        for(i in 2L:marg)
            gradSetArgs[[i]] <- rep_len(TRUE, dimGrad[i-1L])
        useParams <- rep_len(TRUE, dimGrad[marg])
    } else {
        gradSetArgs <- vector("list", 2L)
        useParams <- rep_len(TRUE, length(attr(rhs, "gradient")))
    }
    cat("gradSetArgs:")
    print(gradSetArgs)
    # JN: What is purpose of gradSetArgs??

    npar <- length(useParams)
    cat("npar=",npar,"\n")
    gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
    #JN: Note that this assigns a FORMULA object
    #JN??: we need to be VERY careful to get the right pieces of the model!
    cat("new gradSetArgs[[1L]]:")
    print(gradSetArgs[[1L]])
    #JN: This appears to be the EXPRESSION for the jacobian matrix.

    gradCall <-
        switch(length(gradSetArgs) - 1L,
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], drop = FALSE),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    drop = FALSE),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    gradSetArgs[[3L]], drop = FALSE),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    gradSetArgs[[3L]], gradSetArgs[[4L]], drop = FALSE))
    getRHS.varying <- function()
    {
        ans <- getRHS.noVarying()
        attr(ans, "gradient") <- eval(gradCall)
        ans
    }
    #JN??: next line seems to be working with 1 dimensional model.
    if(length(gr <- attr(rhs, "gradient")) == 1L)
		    attr(rhs, "gradient") <- gr <- as.vector(gr)
    QR <- qr(.swts * gr)
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
