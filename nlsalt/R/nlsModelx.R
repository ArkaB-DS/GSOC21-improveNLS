## This file to generate a lot of output so we can see where information
##  is generated. JN 2021-6-25
nlsModelx <- function(form, data, start, wts, upper=NULL, scaleOffset = 0, nDcentral = FALSE)
{
    ## thisEnv <- environment() # shared by all functions in the 'm' list; variable no longer needed
    env <- new.env(hash = TRUE, parent = environment(form))
    cat("env content (at this stage may be empty):")
    print(ls(env))
    for(i in names(data)) env[[i]] <- data[[i]]
    ind <- as.list(start)
    parLength <- 0L
    for(i in names(ind) ) {
        temp <- start[[i]]
        storage.mode(temp) <- "double"
        env[[i]] <- temp
        ind[[i]] <- parLength + seq_along(temp)
        parLength <- parLength + length(temp)
    }
    cat("After looping over names in data -- env:")
    print(ls(env))

    getPars.noVarying <- function() unlist(mget(names(ind), env))
    getPars <- getPars.noVarying
    cat("getPars:")
    print(getPars)
    cat("Running nlsModel -- just set getPars\n")
    print(str(getPars))
    internalPars <- getPars()
    cat("internalPars:")
    print(internalPars)

    if(!is.null(upper)) upper <- rep_len(upper, parLength)
    # Expand upper to a full vector
    useParams <- rep_len(TRUE, parLength)
    # JN Seems to be a logical vector to indicate free parameters.
    lhs <- eval(form[[2L]], envir = env)
    # set the "promise" for the lhs of the model
    rhs <- eval(form[[3L]], envir = env)
    # similarly for the rhs
    .swts <- if(!missing(wts) && length(wts))
        sqrt(wts) else rep_len(1, length(rhs))
    ##JN: the weights, which are put into the env
    cat(".swts:")
    print(.swts)
    env$.swts <- .swts

    resid <- .swts * (lhs - rhs)
    cat("resid:")
    print(resid)

    dev <- sum(resid^2)
    cat("dev=",dev,"\n")
    if(is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function() {
            if(is.null(upper)) # always for "default"
                numericDeriv(form[[3L]], names(ind), env, central = nDcentral)
            else # possibly with "port"
                numericDeriv(form[[3L]], names(ind), env,
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
        getRHS.noVarying <- function() eval(form[[3L]], envir = env)
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

    getPars.varying <- function() unlist(mget(names(ind), env))[useParams]
    setPars.noVarying <- function(newPars)
    {
        internalPars <<- newPars # envir = thisEnv
        for(i in names(ind))
            env[[i]] <- unname(newPars[ ind[[i]] ])
    }
    setPars.varying <- function(newPars)
    {
        internalPars[useParams] <<- newPars
        for(i in names(ind))
            env[[i]] <- unname(internalPars[ ind[[i]] ])
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
                         vary # envir = thisEnv
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
		 resid <<- .swts * (lhs - (rhs <<- getRHS())) # envir = thisEnv {2 x}
		 dev   <<- sum(resid^2) # envir = thisEnv
		 if(length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
		 QR <<- qr(.swts * gr) # envir = thisEnv
		 (QR$rank < min(dim(QR$qr))) # to catch the singular gradient matrix
	     },
	     getPars = function() getPars(),
	     getAllPars = function() getPars(),
	     getEnv = function() env,
	     trace = function() {
		 d <- getOption("digits")
		 cat(sprintf("%-*s (%.2e): par = (%s)\n", d+4L+2L*(scaleOffset > 0),
			     formatC(dev, digits=d, flag="#"),
			     convCrit(),
			     paste(vapply(getPars(), format, ""), collapse=" ")))
	     },
	     Rmat = function() qr.R(QR),
	     predict = function(newdata = list(), qr = FALSE)
                 eval(form[[3L]], as.list(newdata), env)
	     )
    class(m) <- "nlsModel"
    m
}

