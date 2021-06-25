nlsModel <- function (form, data, start, wts, upper = NULL) 
{
    thisEnv <- environment()
    env <- new.env(hash = TRUE, parent = environment(form))
    for (i in names(data)) assign(i, data[[i]], envir = env)
    ind <- as.list(start)
    parLength <- 0
    for (i in names(ind)) {
        temp <- start[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp, envir = env)
        ind[[i]] <- parLength + seq_along(start[[i]])
        parLength <- parLength + length(start[[i]])
    }
    getPars.noVarying <- function() unlist(setNames(lapply(names(ind), 
        get, envir = env), names(ind)))
    getPars <- getPars.noVarying
    internalPars <- getPars()
    if (!is.null(upper)) upper <- rep_len(upper, parLength)
    useParams <- rep(TRUE, parLength)
    lhs <- eval(form[[2L]], envir = env)
    rhs <- eval(form[[3L]], envir = env)
    .swts <- if (!missing(wts) && length(wts)) sqrt(wts) else rep_len(1, length(rhs))
    assign(".swts", .swts, envir = env)
    resid <- .swts * (lhs - rhs)
    dev <- sum(resid^2)
    if (is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function() {
            if (is.null(upper)) 
                numericDeriv(form[[3L]], names(ind), env)
            else numericDeriv(form[[3L]], names(ind), env, ifelse(internalPars < 
                upper, 1, -1))
        }
        getRHS <- getRHS.noVarying
        rhs <- getRHS()
    }
    else {
        getRHS.noVarying <- function() eval(form[[3L]], envir = env)
        getRHS <- getRHS.noVarying
    }
    dimGrad <- dim(attr(rhs, "gradient"))
    marg <- length(dimGrad)
    if (marg > 0L) {
        gradSetArgs <- vector("list", marg + 1L)
        for (i in 2L:marg) gradSetArgs[[i]] <- rep(TRUE, dimGrad[i - 
            1])
        useParams <- rep(TRUE, dimGrad[marg])
    }
    else {
        gradSetArgs <- vector("list", 2L)
        useParams <- rep(TRUE, length(attr(rhs, "gradient")))
    }
    npar <- length(useParams)
    gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
    gradCall <- switch(length(gradSetArgs) - 1L, call("[", gradSetArgs[[1L]], 
        gradSetArgs[[2L]], drop = FALSE), call("[", gradSetArgs[[1L]], 
        gradSetArgs[[2L]], gradSetArgs[[2L]], drop = FALSE), 
        call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]], 
            gradSetArgs[[3L]], drop = FALSE), call("[", gradSetArgs[[1L]], 
            gradSetArgs[[2L]], gradSetArgs[[2L]], gradSetArgs[[3L]], 
            gradSetArgs[[4L]], drop = FALSE))
    getRHS.varying <- function() {
        ans <- getRHS.noVarying()
        attr(ans, "gradient") <- eval(gradCall)
        ans
    }
    QR <- qr(.swts * attr(rhs, "gradient"))
    qrDim <- min(dim(QR$qr))
    if (QR$rank < qrDim) stop("singular gradient matrix at initial parameter estimates")
    
    getPars.varying <- function() unlist(setNames(lapply(names(ind), 
        get, envir = env), names(ind)))[useParams]
    setPars.noVarying <- function(newPars) {
        assign("internalPars", newPars, envir = thisEnv)
        for (i in names(ind)) assign(i, unname(newPars[ind[[i]]]), 
            envir = env)
    }
    setPars.varying <- function(newPars) {
        internalPars[useParams] <- newPars
        for (i in names(ind)) assign(i, unname(internalPars[ind[[i]]]), 
            envir = env)
    }
    setPars <- setPars.noVarying
    on.exit(remove(i, data, parLength, start, temp, m))
    m <- list(resid = function() resid, fitted = function() rhs, 
        formula = function() form, deviance = function() dev, 
        lhs = function() lhs, gradient = function() .swts * attr(rhs, 
            "gradient"), conv = function() {
            if (npar == 0) return(0)
            rr <- qr.qty(QR, resid)
            sqrt(sum(rr[1L:npar]^2)/sum(rr[-(1L:npar)]^2))
        }, incr = function() qr.coef(QR, resid), setVarying = function(vary = rep(TRUE, 
            length(useParams))) {
            assign("useParams", if (is.character(vary)) {
                temp <- logical(length(useParams))
                temp[unlist(ind[vary])] <- TRUE
                temp
            } else if (is.logical(vary) && length(vary) != length(useParams)) stop("setVarying : 'vary' length must match length of parameters") else {
                vary
            }, envir = thisEnv)
            gradCall[[length(gradCall) - 1L]] <<- useParams
            if (all(useParams)) {
                assign("setPars", setPars.noVarying, envir = thisEnv)
                assign("getPars", getPars.noVarying, envir = thisEnv)
                assign("getRHS", getRHS.noVarying, envir = thisEnv)
                assign("npar", length(useParams), envir = thisEnv)
            } else {
                assign("setPars", setPars.varying, envir = thisEnv)
                assign("getPars", getPars.varying, envir = thisEnv)
                assign("getRHS", getRHS.varying, envir = thisEnv)
                assign("npar", length(seq_along(useParams)[useParams]), 
                  envir = thisEnv)
            }
        }, setPars = function(newPars) {
            setPars(newPars)
            assign("resid", .swts * (lhs - assign("rhs", getRHS(), 
                envir = thisEnv)), envir = thisEnv)
            assign("dev", sum(resid^2), envir = thisEnv)
            assign("QR", qr(.swts * attr(rhs, "gradient")), envir = thisEnv)
            return(QR$rank < min(dim(QR$qr)))
        }, getPars = function() getPars(), getAllPars = function() getPars(), 
        getEnv = function() env, trace = function() {
            cat(format(dev), ": ", format(getPars()))
            cat("\n")
        }, Rmat = function() qr.R(QR), predict = function(newdata = list(), 
            qr = FALSE) eval(form[[3L]], as.list(newdata), env))
    class(m) <- "nlsModel"
    m
}
