#  File src/library/stats/R/nls.R
#     nls-noc.R -- replacement for nls.R
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
###            Nonlinear least squares for R
###
## port_cpos, port_msg() , ... are in  ==> ./nlminb.R

nlsx <-
  function (formula, data = parent.frame(), start, control = nls.control(),
            algorithm = c("default", "plinear", "port"), trace = FALSE,
            subset, weights, na.action, model = FALSE,
            lower = -Inf, upper = Inf, ...)
{
    ## canonicalize the arguments
    formula <- as.formula(formula)
    algorithm <- match.arg(algorithm)

    if(!is.list(data) && !is.environment(data))
        stop("'data' must be a list or an environment")

    mf <- cl <- match.call()		# for creating the model frame
    varNames <- all.vars(formula) # parameter and variable names from formula
    ## adjust a one-sided model formula by using 0 as the response
    if (length(formula) == 2L) {
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }
    ## for prediction we will need to know those which are in RHS
    form2 <- formula; form2[[2L]] <- 0
    varNamesRHS <- all.vars(form2)
    mWeights <- missing(weights)

    ## get names of the parameters from the starting values or selfStart model
    pnames <-
	if (missing(start)) {
	    if(!is.null(attr(data, "parameters"))) {
		names(attr(data, "parameters"))
	    } else { ## try selfStart - like object
		cll <- formula[[length(formula)]]
		if(is.symbol(cll)) { ## replace  y ~ S   by   y ~ S + 0 :
		    ## formula[[length(formula)]] <-
		    cll <- substitute(S + 0, list(S = cll))
		}
		fn <- as.character(cll[[1L]])
		if(is.null(func <- tryCatch(get(fn), error=function(e)NULL)))
		    func <- get(fn, envir=parent.frame()) ## trying "above"
		if(!is.null(pn <- attr(func, "pnames")))
		    as.character(as.list(match.call(func, call = cll))[-1L][pn])
	    }
	} else
	    names(start)

    env <- environment(formula) %||% parent.frame()

    ## Heuristics for determining which names in formula represent actual
    ## variables :

    ## If it is a parameter it is not a variable (nothing to guess here :-)
    if(length(pnames))
        varNames <- varNames[is.na(match(varNames, pnames))]

    ## This aux.function needs to be as complicated because
    ## exists(var, data) does not work (with lists or dataframes):
    lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)),
				     error = function(e) -1L)
    if(length(varNames)) {
        n <- vapply(varNames, lenVar, 0)
        if(any(not.there <- n == -1L)) {
            nnn <- names(n[not.there])
            if(missing(start)) {
                if(algorithm == "plinear")
                    ## TODO: only specify values for the non-lin. parameters
                    stop("no starting values specified")
                ## Provide some starting values instead of erroring out later;
                ## '1' seems slightly better than 0 (which is often invalid):
                warning("No starting values specified for some parameters.\n",
                        "Initializing ", paste(sQuote(nnn), collapse=", "),
                        " to '1.'.\n",
                        "Consider specifying 'start' or using a selfStart model", domain = NA)
		start <- setNames(as.list(rep_len(1., length(nnn))), nnn)
                varNames <- varNames[i <- is.na(match(varNames, nnn))]
                n <- n[i]
            }
            else                        # has 'start' but forgot some
                stop(gettextf("parameters without starting value in 'data': %s",
                              paste(nnn, collapse=", ")), domain = NA)
        }
    }
    else { ## length(varNames) == 0
	if(length(pnames) && any((np <- sapply(pnames, lenVar)) == -1)) {
            ## Can fit a model with pnames even if no varNames
            message(sprintf(ngettext(sum(np == -1),
                                     "fitting parameter %s without any variables",
                                     "fitting parameters %s without any variables"),
                            paste(sQuote(pnames[np == -1]), collapse=", ")),
                    domain = NA)
            n <- integer()
        }
	else
	    stop("no parameters to fit")
    }

    ## If its length is a multiple of the response or LHS of the formula,
    ## then it is probably a variable.
    ## This may fail (e.g. when LHS contains parameters):
    respLength <- length(eval(formula[[2L]], data, env))

    if(length(n) > 0L) {
	varIndex <- n %% respLength == 0
	if(is.list(data) && diff(range(n[names(n) %in% names(data)])) > 0) {
	    ## 'data' is a list that can not be coerced to a data.frame
            ## (not using varNames, varIndex at all - inconsistency FIXME?)
	    mf <- data
            if(!missing(subset))
                warning("argument 'subset' will be ignored")
            if(!missing(na.action))
                warning("argument 'na.action' will be ignored")
	    if(missing(start))
		start <- getInitial(formula, data=mf, control=control, trace=trace)
	    startEnv <- new.env(hash = FALSE, parent = environment(formula)) # small
	    for (i in names(start))
		startEnv[[i]] <- start[[i]]
	    rhs <- eval(formula[[3L]], data, startEnv)
	    n <- NROW(rhs)
            ## mimic what model.frame.default does
	    wts <- if (mWeights) rep_len(1, n)
		   else eval(substitute(weights), data, environment(formula))
	}
        else {
	    vNms <- varNames[varIndex]
	    if(any(nEQ <- vNms != make.names(vNms))) vNms[nEQ] <- paste0("`", vNms[nEQ], "`")
            mf$formula <-  # replace by one-sided linear model formula
		as.formula(paste("~", paste(vNms, collapse = "+")),
                           env = environment(formula))
            mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
            mf$lower <- mf$upper <- NULL
            ## need stats:: for non-standard evaluation
            mf[[1L]] <- quote(stats::model.frame)
            mf <- eval.parent(mf)
            n <- nrow(mf)
            mf <- as.list(mf)
            wts <- if (!mWeights) model.weights(mf) else rep_len(1, n)
        }
        if (any(wts < 0 | is.na(wts)))
            stop("missing or negative weights not allowed")
    }
    else {
        ## length(n) == 0 : Some problems might have no official varNames
        ##                  but still parameters to fit
        varIndex <- logical()
        mf <- list(0)
        wts <- numeric()
    }

    ## set up iteration
    if(missing(start))
        start <- getInitial(formula, data=mf, control=control, trace=trace)
    for(var in varNames[!varIndex])
        mf[[var]] <- eval(as.name(var), data, env)
    varNamesRHS <- varNamesRHS[ varNamesRHS %in% varNames[varIndex] ]

    ## requires 'control' to not contain extra entries (not fulfilled for several CRAN packages)
    ## ctrl <- do.call(nls.control, as.list(if(!missing(control)) control))
    ## Less nice, but more tolerant (to "garbage" which is also put into 'ctrl'):
    ctrl <- nls.control()
    if(!missing(control)) {
	control <- as.list(control)
	ctrl[names(control)] <- control
    }
    scOff  <- ctrl$scaleOffset
    nDcntr <- ctrl$nDcentral
    m <- switch(algorithm,
	plinear = nlsModel.plinear(formula, mf, start, wts, scaleOffset=scOff, nDcentral=nDcntr),
	port = nlsModel (formula, mf, start, wts, upper, scaleOffset=scOff, nDcentral=nDcntr),
        default = nlsModel (formula, mf, start, wts, scaleOffset=scOff, nDcentral=nDcntr))
    cat("Right after setup -- m:\n")
    print(str(m))
    ## Iterate
    if (algorithm != "port") { ## i.e. "default" or  "plinear" :
	if (!identical(lower, -Inf) || !identical(upper, +Inf)) {
	    warning('upper and lower bounds ignored unless algorithm = "port"')
	    cl$lower <- NULL # see PR#15960 -- confint() would use these regardless of algorithm
	    cl$upper <- NULL
	}
#        convInfo <- nlsiter(m, ctrl, trace) -- we will build convInfo
##------ Beginning of nls_iter portion	
# nls_iter(SEXP m, SEXP control, SEXP doTraceArg)
#  nls-noc:      m      ctrl        trace
# {                                 doTrace
# ?? Some of these tests might be better earlier
      if (! is.list(control) ) stop("'control' must be a list")
      if (! is.list(m) ) stop("'m' must be a list")
      if ( is.null(ctrl$maxiter) || (! is.numeric(ctrl$maxiter)) ) stop("Missing ctrl$maxiter")
      if ( is.null(ctrl$tol) || ! is.numeric(ctrl$tol) ) stop("Missing ctrl$tol")
#     
#     conv = getListElement(control, tmp, "minFactor");
      if( is.null(ctrl$minFactor) || ! is.numeric(ctrl$minFactor) ) stop("Missing ctrl$minFactor")

#     conv = getListElement(control, tmp, "warnOnly");
      if( is.null(ctrl$warnOnly) || ! is.logical(ctrl$warnOnly))
            stop("Missing ctrl$warnOnly")
#     
#     conv = getListElement(control, tmp, "printEval");
      if( is.null(ctrl$printEval) || ! is.logical(ctrl$printEval)) stop("Missing ctrl$printEval")
      cat("m:\n")
      print(str(m))
#     
#     // now get parts from 'm'  ---------------------------------
#      tmp = getAttrib(m, R_NamesSymbol);
      tmpa <- attributes(m)
      conv<-m$conv
      if (is.null(conv) || ! is.function(conv)) stop("m$conv missing")
      incr<-m$incr
      if (is.null(incr) || ! is.function(incr)) stop("m$incr missing")
      deviance<-m$deviance
      if (is.null(deviance) || ! is.function(deviance)) stop("m$deviance missing")
      tracefn<-m$trace # BE CAREFUL WITH SAME NAMES !!
      if (is.null(tracefn) || ! is.function(tracefn)) stop("m$trace missing")
      setPars<-m$setPars
      if (is.null(setPars) || ! is.function(setPars)) stop("m$setPars missing")
      getPars<-m$getPars
      if (is.null(getPars) || ! is.function(getPars)) stop("m$getPars missing")
      pars <- eval(getPars(), .GlobalEnv)
#      cat("pars=")
#      print(pars)
      nPars <- length(pars)
      dev <- eval(deviance(), .GlobalEnv)
      cat("dev=",dev,"\n")
      if (trace) {
         cat("result of eval of tracefn:\n")
         eval(tracefn, .GlobalEnv)
      }
      fac <- 1.0;
      hasConverged <- FALSE;
#     SEXP newPars = PROTECT(allocVector(REALSXP, nPars));
      newPars <- rep(NA,nPars)
#     int evaltotCnt = 1;
      evaltotCnt <- 1
      convNew <- -1.0 # double convNew = -1. /* -Wall */; ?? what is -Wall about?
      for (i in 1:ctrl$maxiter){ # top of main iteration -- are there better ways
#       Test for termination
         cat("str(conv):")
         print(str(conv))
         cat("conv:")
         print(conv)
         conv<-as.call(conv)
         cat("str(as.call(conv)):")
         print(str(conv))
         convNew <-eval(conv, .GlobalEnv) # ?? do we need as.call
         cat("str(convNew):")
         print(str(convNew))
         if (convNew <= ctrl$tol) {
            hasConverged <- TRUE
            break;
         }
         # incr computes the delta for Gauss-Newton. Bit incr calls QR, which uses 
         # setPars to create the QR matrix.
         newIncr <- eval(as.call(incr), .GlobalEnv)
	 # BUT PROTECT(incr = lang1(incr)); in C code converts incr to a call with 0 args
         cat("newIncr:")
         print(str(newIncr))
#         double
#         *par   = REAL(pars),
#         *npar  = REAL(newPars),
#         *nIncr = REAL(newIncr);
#         int evalCnt = -1;
         if(ctrl$printEval) evalCnt <- 1
#         
#         while(fac >= minFac) { // 1-dim "line search"
             if(ctrl$printEval) {
#                 Rprintf("  It. %3d, fac= %11.6g, eval (no.,total): (%2d,%3d):",
#                         i+1, fac, evalCnt, evaltotCnt);
#                 evalCnt++;
#                 evaltotCnt++;
             }
#             for(int j = 0; j < nPars; j++)
#                 npar[j] = par[j] + fac * nIncr[j];
#             
#             PROTECT(tmp = lang2(setPars, newPars));
#             if (asLogical(eval(tmp, R_GlobalEnv))) { /* singular gradient */
#                     UNPROTECT(11);
#                 
#                 NON_CONV_FINIS(1, _("singular gradient"));
#             }
#             UNPROTECT(1);
#             
#             double newDev = asReal(eval(deviance, R_GlobalEnv));
#             if(ctrl$printEval)
#                 Rprintf(" new dev = %g\n", newDev);
#             if(newDev <= dev) {
#                 dev = newDev;
#                 tmp = newPars;
#                 newPars = pars;
#                 pars = tmp;
#                 fac = MIN(2*fac, 1);
#                 break;
#             } // else
#                 fac /= 2.;
#         }
#         UNPROTECT(1);
#         if(doTrace) eval(trace, R_GlobalEnv);
#         if( fac < minFac ) {
#             UNPROTECT(9);
#             NON_CONV_FINIS_2(2,
#                              _("step factor %g reduced below 'minFactor' of %g"),
#                              fac, minFac);
#         }
#     }
#     
#     UNPROTECT(9);
#     if(!hasConverged) {
#         NON_CONV_FINIS_1(3,
#                          _("number of iterations exceeded maximum of %d"),
#                          maxIter);
#     }
#     else
#         return CONV_INFO_MSG(_("converged"), 0);
    } # end of main loop for default iteration
# #undef CONV_INFO_MSG
# #undef NON_CONV_FINIS
# #undef NON_CONV_FINIS_1
# #undef NON_CONV_FINIS_2
##------ End of nls_iter portion
        convInfo <- NULL # temporary assignment	
	nls.out <- list(m = m, convInfo = convInfo,
			data = substitute(data), call = cl)
    } # end of default and plinear options
    else { ## "port" i.e., PORT algorithm
	pfit <- nls_port_fit(m, start, lower, upper, control, trace,
			     give.v=TRUE)
        iv <- pfit[["iv"]]
	msg.nls <- port_msg(iv[1L])
	conv <- (iv[1L] %in% 3:6)
	if (!conv) {
	    msg <- paste("Convergence failure:", msg.nls)
	    if(ctrl$warnOnly) warning(msg) else stop(msg)
	}
	v. <- port_get_named_v(pfit[["v"]])
	## return a 'convInfo' list compatible to the non-PORT case:
	cInfo <- list(isConv = conv,
		      finIter = iv[31L], # 31: NITER
		      finTol  =	 v.[["NREDUC"]],
		      nEval = c("function" = iv[6L], "gradient" = iv[30L]),
		      stopCode = iv[1L],
		      stopMessage = msg.nls)
        ## we need these (evaluated) for profiling
	cl$lower <- lower
	cl$upper <- upper
	nls.out <- list(m = m, data = substitute(data),
                        call = cl, convInfo = cInfo,
                        ## UGLY: this is really a logical for  *NON*convergence:
                        ## deprecate these two, as they are now part of convInfo
			convergence = as.integer(!conv),
			message = msg.nls)
    }

    ## we need these (evaluated) for profiling
    nls.out$call$algorithm <- algorithm
    nls.out$call$control <- ctrl
    nls.out$call$trace <- trace

    nls.out$na.action <- attr(mf, "na.action")
    nls.out$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
    if(model) nls.out$model <- mf
    if(!mWeights) nls.out$weights <- wts
    nls.out$control <- control
    class(nls.out) <- "nls"
    nls.out
}
