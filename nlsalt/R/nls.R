#  File src/library/stats/R/nls.R
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

nls <-
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
		plinear = nlsModel.plinear(formula, mf, start, wts,        scaleOffset=scOff, nDcentral=nDcntr),
		port    = nlsModel        (formula, mf, start, wts, upper, scaleOffset=scOff, nDcentral=nDcntr),
                default = nlsModel        (formula, mf, start, wts,        scaleOffset=scOff, nDcentral=nDcntr))

    ## Iterate
    if (algorithm != "port") { ## i.e. "default" or  "plinear" :
	if (!identical(lower, -Inf) || !identical(upper, +Inf)) {
	    warning('upper and lower bounds ignored unless algorithm = "port"')
	    cl$lower <- NULL # see PR#15960 -- confint() would use these regardless of algorithm
	    cl$upper <- NULL
	}
        convInfo <- .Call(C_nls_iter, m, ctrl, trace)
##------ Beginning of nls_iter portion	
# /*
#     *  call to nls_iter from R --- .Call("nls_iter", m, control, doTrace)
# *  where m and control are nlsModel and nlsControl objects
# *             doTrace is a logical value.
# *  m is modified; the return value is a "convergence-information" list.
# */
#     SEXP
# nls_iter(SEXP m, SEXP control, SEXP doTraceArg)
# {
#     int doTrace = asLogical(doTraceArg);
#     
#     if(!isNewList(control))
#         error(_("'control' must be a list"));
#     if(!isNewList(m))
#         error(_("'m' must be a list"));
#     
#     SEXP tmp, conv;
#     PROTECT(tmp = getAttrib(control, R_NamesSymbol));
#     
#     conv = getListElement(control, tmp, "maxiter");
#     if(conv == NULL || !isNumeric(conv))
#         error(_("'%s' absent"), "control$maxiter");
#     int maxIter = asInteger(conv);
#     
#     conv = getListElement(control, tmp, "tol");
#     if(conv == NULL || !isNumeric(conv))
#         error(_("'%s' absent"), "control$tol");
#     double tolerance = asReal(conv);
#     
#     conv = getListElement(control, tmp, "minFactor");
#     if(conv == NULL || !isNumeric(conv))
#         error(_("'%s' absent"), "control$minFactor");
#     double minFac = asReal(conv);
#     
#     conv = getListElement(control, tmp, "warnOnly");
#     if(conv == NULL || !isLogical(conv))
#         error(_("'%s' absent"), "control$warnOnly");
#     int warnOnly = asLogical(conv);
#     
#     conv = getListElement(control, tmp, "printEval");
#     if(conv == NULL || !isLogical(conv))
#         error(_("'%s' absent"), "control$printEval");
#     Rboolean printEval = asLogical(conv);
#     
#     // now get parts from 'm'  ---------------------------------
#         tmp = getAttrib(m, R_NamesSymbol);
#     
#     conv = getListElement(m, tmp, "conv");
#     if(conv == NULL || !isFunction(conv))
#         error(_("'%s' absent"), "m$conv()");
#     PROTECT(conv = lang1(conv));
#     
#     SEXP incr = getListElement(m, tmp, "incr");
#     if(incr == NULL || !isFunction(incr))
#         error(_("'%s' absent"), "m$incr()");
#     PROTECT(incr = lang1(incr));
#     
#     SEXP deviance = getListElement(m, tmp, "deviance");
#     if(deviance == NULL || !isFunction(deviance))
#         error(_("'%s' absent"), "m$deviance()");
#     PROTECT(deviance = lang1(deviance));
#     
#     SEXP trace = getListElement(m, tmp, "trace");
#     if(trace == NULL || !isFunction(trace))
#         error(_("'%s' absent"), "m$trace()");
#     PROTECT(trace = lang1(trace));
#     
#     SEXP setPars = getListElement(m, tmp, "setPars");
#     if(setPars == NULL || !isFunction(setPars))
#         error(_("'%s' absent"), "m$setPars()");
#     PROTECT(setPars);
#     
#     SEXP getPars = getListElement(m, tmp, "getPars");
#     if(getPars == NULL || !isFunction(getPars))
#         error(_("'%s' absent"), "m$getPars()");
#     PROTECT(getPars = lang1(getPars));
#     
#     SEXP pars = PROTECT(eval(getPars, R_GlobalEnv));
#     int nPars = LENGTH(pars);
#     
#     double dev = asReal(eval(deviance, R_GlobalEnv));
#     if(doTrace) eval(trace,R_GlobalEnv);
#     
#     double fac = 1.0;
#     Rboolean hasConverged = FALSE;
#     SEXP newPars = PROTECT(allocVector(REALSXP, nPars));
#     int evaltotCnt = 1;
#     double convNew = -1. /* -Wall */;
#     int i;
#     #define CONV_INFO_MSG(_STR_, _I_)				\
#     ConvInfoMsg(_STR_, i, _I_, fac, minFac, maxIter, convNew)
#     
#     #define NON_CONV_FINIS(_ID_, _MSG_)		\
#     if(warnOnly) {				\
#         warning(_MSG_);				\
#         return CONV_INFO_MSG(_MSG_, _ID_);      \
#     }						\
#     else					\
#     error(_MSG_);
#     
#     #define NON_CONV_FINIS_1(_ID_, _MSG_, _A1_)	\
#     if(warnOnly) {				\
#         char msgbuf[1000];			\
#         warning(_MSG_, _A1_);			\
#         snprintf(msgbuf, 1000, _MSG_, _A1_);	\
#         return CONV_INFO_MSG(msgbuf, _ID_);	\
#     }						\
#     else					\
#     error(_MSG_, _A1_);
#     
#     #define NON_CONV_FINIS_2(_ID_, _MSG_, _A1_, _A2_)	\
#     if(warnOnly) {					\
#         char msgbuf[1000];				\
#         warning(_MSG_, _A1_, _A2_);			\
#         snprintf(msgbuf, 1000, _MSG_, _A1_, _A2_);	\
#         return CONV_INFO_MSG(msgbuf, _ID_);		\
#     }							\
#     else						\
#     error(_MSG_, _A1_, _A2_);
#     
#     for (i = 0; i < maxIter; i++) { // ---------------------------------------------
#             
#             if((convNew = asReal(eval(conv, R_GlobalEnv))) <= tolerance) {
#                 hasConverged = TRUE;
#                 break;
#             }
#         
#         SEXP newIncr = PROTECT(eval(incr, R_GlobalEnv));
#         double
#         *par   = REAL(pars),
#         *npar  = REAL(newPars),
#         *nIncr = REAL(newIncr);
#         int evalCnt = -1;
#         if(printEval)
#             evalCnt = 1;
#         
#         while(fac >= minFac) { // 1-dim "line search"
#             if(printEval) {
#                 Rprintf("  It. %3d, fac= %11.6g, eval (no.,total): (%2d,%3d):",
#                         i+1, fac, evalCnt, evaltotCnt);
#                 evalCnt++;
#                 evaltotCnt++;
#             }
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
#             if(printEval)
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
# }
# #undef CONV_INFO_MSG
# #undef NON_CONV_FINIS
# #undef NON_CONV_FINIS_1
# #undef NON_CONV_FINIS_2
##------ End of nls_iter portion	
	nls.out <- list(m = m, convInfo = convInfo,
			data = substitute(data), call = cl)
	}
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
    nls.out$dataClasses <-
        attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
    if(model)
	nls.out$model <- mf
    if(!mWeights)
	nls.out$weights <- wts
    nls.out$control <- control
    class(nls.out) <- "nls"
    nls.out
}
