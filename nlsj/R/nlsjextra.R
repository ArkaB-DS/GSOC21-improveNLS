coef.nlsj <- function(object, ...) object$m$getPars() # not getAllPars since no plinear yet

summary.nlsj <-
    function (object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    r <- as.vector(object$m$resid()) # These are weighted residuals.
    w <- object$weights
    n <- if (!is.null(w)) sum(w > 0) else length(r)
    param <- coef(object)
    pnames <- names(param)
    p <- length(param)
    rdf <- n - p
    resvar <- if(rdf <= 0) NaN else deviance(object)/rdf
    XtXinv <- chol2inv(object$m$Rmat())
    dimnames(XtXinv) <- list(pnames, pnames)
    se <- sqrt(diag(XtXinv) * resvar)
    tval <- param/se
    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <-
        list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans <- list(formula = formula(object), residuals = r, sigma = sqrt(resvar),
                df = c(p, rdf), cov.unscaled = XtXinv,
                call = object$call,
                convInfo = object$convInfo,
                control = object$control,
                na.action = object$na.action,
                coefficients = param,
                parameters = param)# never documented, for back-compatibility
    if(correlation && rdf > 0) {
        ans$correlation <- (XtXinv * resvar)/outer(se, se)
        ans$symbolic.cor <- symbolic.cor
    }
    ## if(identical(object$call$algorithm, "port"))
    ##     ans$message <- object$message
    class(ans) <- "summary.nlsj"
    ans
}

print.nlsj <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("Nonlinear regression model\n")
    cat("  model: ", deparse(formula(x)), "\n", sep = "")
    cat("   data: ", deparse(x$data), "\n", sep = "")
#    print(x$m$getAllPars(), digits = digits, ...)
#?? getAllPars includes plinear ones??
     print(x$m$getPars(), digits = digits, ...)
    cat(" ", if(!is.null(x$weights) && diff(range(x$weights))) "weighted ",
	"residual sum-of-squares: ", format(x$m$deviance(), digits = digits),
	"\n", sep = "")
#    .p.nlsj.convInfo(x, digits = digits)
    print(x$convInfo)
    invisible(x)
}

print.summary.nlsj <-
  function (x, digits = max(3L, getOption("digits") - 3L),
            symbolic.cor = x$symbolic.cor,
            signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nFormula: ",
	paste(deparse(x$formula), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    df <- x$df
    rdf <- df[2L]
    cat("\nParameters:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                 ...)
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
    cat("\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Parameter Estimates:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {
		print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2), nsmall = 2L, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
        }
    }

#    .p.nlsj.convInfo(x, digits = digits)
    print(x$convInfo)

    if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("\n")
    invisible(x)
}

weights.nlsj <- function(object, ...) object$weights

predict.nlsj <-
  function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
           interval = c("none", "confidence", "prediction"), level = 0.95,
           ...)
{
    if (missing(newdata)) return(as.vector(fitted(object)))
    if(!is.null(cl <- object$dataClasses)) .checkMFClasses(cl, newdata)
    object$m$predict(newdata)
}

fitted.nlsj <- function(object, ...)
{
    val <- as.vector(object$m$fitted())
    if(!is.null(object$na.action)) val <- napredict(object$na.action, val)
    lab <- "Fitted values"
    if (!is.null(aux <- attr(object, "units")$y)) lab <- paste(lab, aux)
    attr(val, "label") <- lab
    val
}

formula.nlsj <- function(x, ...) x$m$formula()

residuals.nlsj <- function(object, type = c("response", "pearson"), ...)
{ # 210723 -- ?? error. Need to verify!!
    type <- match.arg(type)
    if (type == "pearson") {
        val <- as.vector(object$m$resid())
        std <- sqrt(sum(val^2)/(length(val) - length(coef(object))))
        val <- val/std
        if(!is.null(object$na.action)) val <- naresid(object$na.action, val)
        attr(val, "label") <- "Standardized residuals"
    } else {
##??WRONG, but why?        val <- as.vector(object$m$lhs() - object$m$fitted())
        val <- as.numeric(object$m$lhs()-object$m$fitted())
        if(!is.null(object$na.action))
            val <- naresid(object$na.action, val)
        lab <- "Residuals"
        if (!is.null(aux <- attr(object, "units")$y)) lab <- paste(lab, aux)
        attr(val, "label") <- lab
    }
    val
}

logLik.nlsj <- function(object, REML = FALSE, ...)
{
    if (REML)
        stop("cannot calculate REML log-likelihood for \"nls\" objects")
    res <- object$m$resid() # These are weighted residuals.
    N <- length(res)
    w <- object$weights %||% rep_len(1, N)
    ## Note the trick for zero weights
    zw <- w == 0
    N <- sum(!zw)
    val <-  -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(res^2)))/2
    ## the formula here corresponds to estimating sigma^2.
    attr(val, "df") <- 1L + length(coef(object))
    attr(val, "nobs") <- attr(val, "nall") <- N
    class(val) <- "logLik"
    val
}

df.residual.nlsj <- function(object, ...)
{
    w <- object$weights
    n <- if(!is.null(w)) sum(w != 0) else length(object$m$resid())
    n - length(coef(object))
}

deviance.nlsj <- function(object, ...) object$m$deviance()

vcov.nlsj <- function(object, ...)
{
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}


anova.nlsj <- function(object, ...)
{
    if(length(list(object, ...)) > 1L) return(anovalist.nlsj(object, ...))
    stop("anova is only defined for sequences of \"nls\" objects")
}

anovalist.nlsj <- function (object, ..., test = NULL)
{
    objects <- list(object, ...)
    responses <- as.character(lapply(objects,
				     function(x) formula(x)[[2L]]))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
	objects <- objects[sameresp]
        warning(gettextf("models with response %s removed because response differs from model 1",
                         sQuote(deparse(responses[!sameresp]))),
                domain = NA)
    }
    ## calculate the number of models
    nmodels <- length(objects)
    if (nmodels == 1L)
        stop("'anova' is only defined for sequences of \"nls\" objects")

    models <- as.character(lapply(objects, function(x) formula(x)))

    ## extract statistics
    df.r <- unlist(lapply(objects, df.residual))
    ss.r <- unlist(lapply(objects, deviance))
    df <- c(NA, -diff(df.r))
    ss <- c(NA, -diff(ss.r))
    ms <- ss/df
    f <- p <- rep_len(NA_real_, nmodels)
    for(i in 2:nmodels) {
	if(df[i] > 0) {
	    f[i] <- ms[i]/(ss.r[i]/df.r[i])
	    p[i] <- pf(f[i], df[i], df.r[i], lower.tail = FALSE)
	}
	else if(df[i] < 0) {
	    f[i] <- ms[i]/(ss.r[i-1]/df.r[i-1])
	    p[i] <- pf(f[i], -df[i], df.r[i-1], lower.tail = FALSE)
	}
	else {                          # df[i] == 0
            ss[i] <- 0
	}
    }
    table <- data.frame(df.r,ss.r,df,ss,f,p)
    dimnames(table) <- list(1L:nmodels, c("Res.Df", "Res.Sum Sq", "Df",
					 "Sum Sq", "F value", "Pr(>F)"))
    ## construct table and title
    title <- "Analysis of Variance Table\n"
    topnote <- paste0("Model ", format(1L:nmodels), ": ", models,
                      collapse = "\n")

    ## calculate test statistic if needed
    structure(table, heading = c(title, topnote),
	      class = c("anova", "data.frame")) # was "tabular"
}

###
### asOneSidedFormula is extracted from the NLME-3.1 library for S
###

# asOneSidedFormula <-
#   ## Converts an expression or a name or a character string
#   ## to a one-sided formula
#   function(object)
# {
#     if ((mode(object) == "call") && (object[[1L]] == "~")) {
#         object <- eval(object)
#     }
#     if (inherits(object, "formula")) {
#         if (length(object) != 2L) {
#             stop(gettextf("formula '%s' must be of the form '~expr'",
#                           deparse(as.vector(object))), domain = NA)
#         }
#         return(object)
#     }
#     do.call("~",
#             list(switch(mode(object),
#                         name = ,
#                         numeric = ,
#                         call = object,
#                         character = as.name(object),
#                         expression = object[[1L]],
#                         stop(gettextf("'%s' cannot be of mode '%s'",
#                                       substitute(object), mode(object)),
#                              domain = NA)
#                         ))
#             )
# }
# 
# ## "FIXME": move to 'base' and make .Internal or even .Primitive
# setNames <- function(object = nm, nm)
# {
#     names(object) <- nm
#     object
# }
# nobs.nlsj <- function(object, ...){
#     if (is.null(w <- object$weights)) length(object$m$resid()) else sum(w != 0)
# }