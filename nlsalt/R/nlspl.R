#  File src/library/stats/R/nlspl.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 2000-2020 The R Core Team
#  Copyright (C) 1999-1999 Saikat DebRoy, Douglas M. Bates, Jose C. Pinheiro
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

###
###            plinear Model for Nonlinear least squares for R
###

nlsModel.plinear <- function(form, data, start, wts, scaleOffset = 0, nDcentral = FALSE)
{
    ## thisEnv <- environment() # shared by all functions in the 'm' list; variable no longer needed
    env <- new.env(hash = TRUE, parent=environment(form))
    for(i in names(data)) env[[i]] <- data[[i]]
    ind <- as.list(start)
    p2 <- 0L #{non-linear parameters}
    for(i in names(ind)) {
        temp <- start[[i]]
        storage.mode(temp) <- "double"
        env[[i]] <- temp
        ind[[i]] <- p2 + seq_along(temp)
        p2 <- p2 + length(temp)
    }
    lhs <- eval(form[[2L]], envir = env); storage.mode(lhs) <- "double"
    rhs <- eval(form[[3L]], envir = env); storage.mode(rhs) <- "double"
    .swts <- if(!missing(wts) && length(wts))
        sqrt(wts) else 1 # more efficient than  rep_len(1, NROW(rhs))
    env$.swts <- .swts
    p1 <- NCOL(rhs) #{linear par.}
    p <- p1 + p2    # total #{param.}
    n <- length(lhs)
    fac <- (n -  p)/p
    cc <- QR.B <- NA
    useParams <- rep_len(TRUE, p2)
    if(is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function()
            numericDeriv(form[[3L]], names(ind), env, central = nDcentral)
        getRHS <- getRHS.noVarying
        rhs <- getRHS()
    } else {
        getRHS.noVarying <- function() eval(form[[3L]], envir = env)
        getRHS <- getRHS.noVarying
    }
    dimGrad <- dim(attr(rhs, "gradient"))
    marg <- length(dimGrad)
    if(marg > 0) {
        if(marg < 2L) stop("invalid  'attr(rhs, \"gradient\")'")
        gradSetArgs <- vector("list", marg + 1L)
        for(i in 2:marg)
            gradSetArgs[[i]] <- rep_len(TRUE, dimGrad[i-1L])
        useParams <- rep_len(TRUE, dimGrad[marg])
    } else {
        gradSetArgs <- vector("list", 2L)
        useParams <- rep_len(TRUE, length(attr(rhs, "gradient")))
    }
    gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
    gradCall <-
        switch(length(gradSetArgs) - 1L,
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]]),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]]),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]], gradSetArgs[[3L]]),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]], gradSetArgs[[3L]], gradSetArgs[[4L]]))
    getRHS.varying <- function()
    {
        ans <- getRHS.noVarying()
        attr(ans, "gradient") <- eval(gradCall)
        ans
    }
    QR.rhs <- qr(.swts * rhs)
    lin <- qr.coef(QR.rhs, .swts * lhs)
    resid <- qr.resid(QR.rhs, .swts * lhs)
    topzero <- double(p1)
    dev <- sum(resid^2)
    if(marg <= 1) {
        ddot <- function(A, b) A %*% b
        dtdot <- function(A, b) t(A) %*% b
    } else if(marg == 2) {
        if(p1 == 1) {
            ddot <- function(A, b) as.matrix(A*b)
            dtdot <- function(A, b) t(b) %*% A
        } else if(p2 == 1) {
            ddot  <- function(A, b)   A  %*% b
            dtdot <- function(A, b) t(A) %*% b
        }
    } else {
        ddot  <- function(A, b) apply(A, MARGIN =   3L    , FUN = "%*%", b)
        dtdot <- function(A, b) apply(A, MARGIN = c(2L,3L), FUN = "%*%", b)
    }

    getPars.noVarying <- function() unlist(mget(names(ind), env))
    getPars.varying   <- function() unlist(mget(names(ind), env))[useParams]
    getPars <- getPars.noVarying

    internalPars <- getPars()
    setPars.noVarying <- function(newPars)
    {
        internalPars <<- newPars # into thisEnv
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
    getPred <-
        if(is.matrix(rhs)) function(X) as.numeric(X %*% lin)
        else function(X) X * lin

    if(scaleOffset) scaleOffset <- (n-p1) * scaleOffset^2
    convCrit <- function() {
        cc <<- c(topzero, qr.qty(QR.rhs, .swts * lhs)[ -(1L:p1)]) # envir = thisEnv
        rr <- qr.qy(QR.rhs, cc)
        B <- qr.qty(QR.rhs, .swts * ddot(attr(rhs, "gradient"), lin))
        B[1L:p1, ] <- dtdot(.swts *      attr(rhs, "gradient"), rr)
        R <- t( qr.R(QR.rhs)[1L:p1, ] )
        if(p1 == 1)
            B[1L,] <- B[1L,]/c(R)
        else
            B[1L:p1, ] <- forwardsolve(R, B[1L:p1, ])
        QR.B <<- qr(B) ## envir = thisEnv
        rr <- qr.qty(QR.B, cc)
        sqrt( fac*sum(rr[1L:p1]^2) / (scaleOffset + sum(rr[-(1L:p1)]^2)) )
    }

    m <-
        list(resid = function() resid,
             fitted = function() getPred(rhs),
             formula = function() form,
             deviance = function() dev,
             lhs = function() lhs,
             gradient = function() attr(rhs, "gradient"),
	     conv = function() convCrit(),
             incr = function() qr.solve(QR.B, cc),
             setVarying = function(vary = rep_len(TRUE, np)) {
                 np <- length(useParams)
                 useParams <<-
                     if(is.character(vary)) {
                         temp <- logical(np)
                         temp[unlist(ind[vary])] <- TRUE
                         temp
                     } else if(is.logical(vary) && length(vary) != np)
                        stop("setVarying : 'vary' length must match length of parameters")
                     else
                         vary # envir = thisEnv
                 gradCall[[length(gradCall)]] <<- useParams
                 if(all(useParams)) {
		     setPars <<- setPars.noVarying
		     getPars <<- getPars.noVarying
		     getRHS  <<-  getRHS.noVarying
                 } else {
		     setPars <<- setPars.varying
		     getPars <<- getPars.varying
		     getRHS  <<-  getRHS.varying
                 }
             },
             setPars = function(newPars) {
                 setPars(newPars)
                 QR.rhs <<- qr(.swts * (rhs <<- getRHS())) # envir = thisEnv
                 resid <<- qr.resid(QR.rhs, .swts * lhs) # envir = thisEnv
                 dev <<- sum(resid^2) # envir = thisEnv
                 if(QR.rhs$rank < p1) { # singular gradient
                     TRUE
                 } else {
                     lin <<- qr.coef(QR.rhs, .swts * lhs) # envir = thisEnv
                     FALSE
                 }
             },
             getPars = function() getPars(),
             getAllPars = function() c( getPars(), c( .lin = lin ) ),
	     getEnv = function() env,
	     trace = function() {
		 d <- getOption("digits")
		 cat(sprintf("%-*s (%.2e): par = (%s)\n", d+4L+2L*(scaleOffset > 0),
			     formatC(dev, digits=d, flag="#"),
			     convCrit(),
			     paste(vapply(c(getPars(), lin), format, ""), collapse=" ")))
	     },
             Rmat = function()
                 qr.R(qr(.swts * cbind(ddot(attr(rhs, "gradient"), lin), rhs))),
             predict = function(newdata = list(), qr = FALSE)
                 getPred(eval(form[[3L]], as.list(newdata), env))
             )
    class(m) <- c("nlsModel.plinear", "nlsModel")
    m$conv()
    on.exit(remove(data, i, m, marg, dimGrad, n, p, p2, start, temp, gradSetArgs) )
    m
}

