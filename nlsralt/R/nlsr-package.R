# nlsr-package.R -- print and summary methods

summary.nlsr <- function(object, ...) {
    smalltol <- .Machine$double.eps * 1000
    options(digits = 5) # 7 is default
    resname <- deparse(substitute(object))
    JJ <- object$jacobian
    res <- object$resid
    coeff <- object$coefficients
    pnames<-names(coeff)
##    param <- coef(object)
    npar <- length(coeff)
    w <- object$weights
    nobs <- if (!is.null(w)) sum(w > 0) else length(res)
    rdf <- nobs - npar
    lo <- object$lower
    if (is.null(lo)) lo <- rep( -Inf, npar)
    up <- object$upper
    if (is.null(up)) up <- rep( Inf, npar)
    mi <- object$maskidx
    mt <- rep(" ",npar) # start with all "unmasked"
    mt[mi] <- "M" # Put in the masks
    bdmsk <- rep(1, npar) # bounds and masks indicator ?? should it be 1L
    bdmsk[mi] <- 0 # masked
    ct <- rep(" ",npar) # start with all "free"
    for (i in seq_along(coeff)){
       if (lo[[i]] - coeff[[i]] > 0) {
          ct[[i]] <- "-" # lower bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -3
       } else { 
          if (coeff[[i]] - lo[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
             ct[[i]] <- "L" # "at" lower bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -3 # leave mask indication intact
          }
       }
       if (coeff[[i]] - up[[i]] > 0) {
          ct[[i]] <- "+" # upper bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -1
       } else { 
          if (up[[i]] - coeff[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
             ct[[i]] <- "U" # "at" upper bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -1 # leave mask indication intact
          }
       }
    }
    ss <- object$ssquares
    rdf <- nobs - npar
    if (rdf <= 0) {
          if (rdf < 0) { stop(paste("Inadmissible degrees of freedom =",rdf,sep='')) }
          else { resvar <- Inf }
    } else {
       resvar <- ss/(rdf)
    }
    dec <- svd(JJ)
    U <- dec$u
    V <- dec$v
    Sd <- dec$d
    if (min(Sd) <= smalltol * max(Sd)) { # singular
       SEs <- rep(NA, npar) # ?? Inf or NA
       XtXinv <- matrix(NA, nrow=npar, ncol=npar)
    } else {
       Sinv <- 1/Sd
##       Sinv[which(bdmsk != 1)] <- 0 # ?? 140714 maybe don't want this
       if (npar > 1) {
           VS <- crossprod(t(V), diag(Sinv))
       } else {
           VS <- V/Sinv
       }
       XtXinv <- crossprod(t(VS))
       SEs <- sqrt(diag(XtXinv) * resvar)
    }
##    cat("CHECK XtXinv:")
##    print(XtXinv)
    dimnames(XtXinv) <- list(pnames, pnames)
    gr <- crossprod(JJ, res)
    if (any(is.na(SEs))) {
        tstat<-rep(NA, npar)
    } else {
        tstat <- coeff/SEs
    }
    tval <- coeff/SEs
    param <- cbind(coeff, SEs, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <-
        list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

# Note: We don't return formula because we may be doing nlfb summary 
#   i.e., resfn and jacfn approach  ?? But could we??
    ans <- list(residuals = res, sigma = sqrt(resvar),  
                df = c(npar, rdf), cov.unscaled = XtXinv,
                param = param, resname=resname, ssquares=ss, nobs=nobs, 
                ct=ct, mt=mt, Sd=Sd, gr=gr, jeval=object$jeval,feval=object$feval)
                ## parameters = param)# never documented, for back-compatibility
    attr(ans,"pkgname") <- "nlsr"
    class(ans) <- "summary.nlsr"
    ans
}

# coef() function
coef.nlsr <- function(object, ...) {
       out <- object$coefficients
       attr(out,"pkgname")<-"nlsr"
##       invisible(out)
       out # JN 170109
}

print.nlsr <- function(x, ...) {
    xx<-summary(x)
    with(xx, { 
	cat("nlsr object:",resname,"\n")
	pname<-dimnames(param)[[1]] # param is augmented coefficients with SEs and tstats
	npar <- dim(param)[1] # ?? previously length(coeff) 
        cat("residual sumsquares = ",ssquares," on ",nobs,"observations\n")
        cat("    after ",jeval,"   Jacobian and ",feval,"function evaluations\n")
        cat("  name     ","      coeff    ","     SE   ","   tstat  ",
             "   pval  ","   gradient  "," JSingval  ","\n")
        SEs <- param[,2]
        tstat <- param[,3]
        pval <- param[,4]
        for (i in seq_along(param[,1])){
            tmpname<-pname[i]
            if (is.null(tmpname)) {tmpname <- paste("p_",i,sep='')}
            cat(format(tmpname, width=10)," ")
            cat(format(param[[i]], digits=6, width=12))
            cat(ct[[i]],mt[[i]]," ")
            cat(format(SEs[[i]], digits=4, width=9)," ")
            cat(format(tstat[[i]], digits=4, width=9)," ")
            cat(format(pval[[i]], digits=4, width=9)," ")
            cat(format(gr[[i]], digits=4, width=10)," ")
            cat(format(Sd[[i]], digits=4, width=10)," ")
            cat("\n")
        }
    }) # remember to close with()
  invisible(x)
}


## ?? can we do a generic resids??
res <- function(object){
  resids <- object$resid
  resids # so the function prints
}

predict.nlsr <- function(object=NULL, newdata=list(), ...) { 
#  This ONLY works if we have used nlxb. Do we want to add class 'nlxb'??
    if( is.null(object) ) stop("predict.nlsr REQUIRES an nlsr solution object")
    form <- object$formula
    if (is.null(form)) stop("nlsr.predict works only if formula is defined")
# ?? give more output
#
#  we assume a formula of style y~something, and use the something
#  In some ways need to check this more carefully
    env4coefs <- list2env(as.list(object$coefficients))
    preds <- eval(form[[3]], as.list(newdata), env4coefs)
    class(preds)<- "predict.nlsr"
    attr(preds,"pkgname") <- "nlsr"
    preds
}
