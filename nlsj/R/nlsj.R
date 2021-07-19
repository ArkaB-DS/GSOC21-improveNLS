nlsj <- function (formula, data = parent.frame(), start, control = nlsj.control(),
            algorithm = "default", weights=NULL, subset, trace = FALSE,
            na.action, model=FALSE, lower = -Inf, upper = Inf, ...) {
# ?? left out -- FIXME??    na.action, model = FALSE, (masked from nlxb)
# ?? at this stage ONLY treat "default", but will add bounds
# ?? data in .GlobalEnv -- should be OK
##?? Should a lot of this material be in nlsjModel() to build tools for problem??
# DOES NOT CALL nlsjModel, but does everything here

# Controls
#   cat("control:"); print(control)
   if (is.null(control$derivmeth)) control$derivmeth="default" # for safety
   epstol <- (.Machine$double.eps * control$offset) 
   epstol4 <- epstol^4 # used for smallsstest
   ##?? may want these in nlsj.control
   if (control$derivmeth == "numericDeriv") warning("Forcing numericDeriv")

# Algorithm
   if (is.null(algorithm)) algorithm<-"default"
   algchoice <- c("default", "port", "plinear", "marquardt")
   alg <- which(algorithm %in% algchoice)
   switch(alg,
      "default" = if (trace) cat("nlsj: Using default algorithm\n"),
      # rest
      { msg <- paste("Algorithm choice '",alg,"' not yet implemented or does not exist")
        stop(msg) }
   ) # end switch choices

   
# Data
    stopifnot(inherits(formula, "formula"))
    if (missing(data)) {  ## rather than   #  if (is.null(data)) {
        warning("Data is not declared explicitly. Caution!")
	data <- environment(formula) # this will handle variables in the parent frame
    }
    else if (is.list(data)){
	    data <- list2env(data, parent = environment(formula))
        }
        else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")
   cat("ls(data):"); print(ls(data))
   vnames <- all.vars(formula) # all names in the formula
   dnames <- vnames[which(vnames %in% ls(data))]
   if (length(dnames) < 1) stop("No data found")
# ?? This fails when we have the parameters in variables as well. So we need
# ?? to figure out which are the true "variables" of the problem and which are
# ?? parameters.
   pnames <- vnames[ - which(vnames %in% dnames)] # the "non-data" names in the formula
#   ll <- vapply(vnames, function(xx) {length(eval(parse(text=xx)))}, numeric(1) )
#   pnames <- vnames[which(ll == 1)] ?? doesn't seem to be working
   npar <- length(pnames)
   cat("vnames:"); print(vnames)
   cat("npar=",npar," pnames:"); print(pnames)

# Start vector
   if (is.null(start)) { # start not specified
     warning("start vector not specified for nlsj")
     start<-0.9+seq(npar)*0.1 # WARNING: very crude??
     names(start)<-pnames # and make sure these are named?? necessary??
     ## ??? put in ways to get at selfstart models 
   }
   else { # we have a start vector
     snames<-names(start) # names in the start vector
#     if ((length(snames) != length(pnames)) || (! all.equal(snames, pnames))) {
#       cat("snames:"); print(snames)
#       cat("pnames:"); print(pnames)
#       stop("Start names differ in number or name from formula parameter names")
#     }
     start <- as.numeric(start) # ensure we convert (e.g., if matrix)
     names(start) <- snames ## as.numeric strips names, so this is needed ??
   }
   prm <- start # start MUST be defined at this point
   localdata <- list2env(as.list(prm), parent = data)

   
   # ?? How nls() gets pnames
   ## get names of the parameters from the starting values or selfStart model
   # pnames <-
   #   if (missing(start)) {
   #     if(!is.null(attr(data, "parameters"))) {
   #       names(attr(data, "parameters"))
   #     } else { ## try selfStart - like object
   #       cll <- formula[[length(formula)]]
   #       if(is.symbol(cll)) { ## replace  y ~ S   by   y ~ S + 0 :
   #         ## formula[[length(formula)]] <-
   #         cll <- substitute(S + 0, list(S = cll))
   #       }
   #       fn <- as.character(cll[[1L]])
   #       if(is.null(func <- tryCatch(get(fn), error=function(e)NULL)))
   #         func <- get(fn, envir=parent.frame()) ## trying "above"
   #       if(!is.null(pn <- attr(func, "pnames")))
   #         as.character(as.list(match.call(func, call = cll))[-1L][pn])
   #     }
   #   } else
   #     names(start)
   # 
   



# Weights
   mraw <- length(eval(as.name(dnames[1])))
   if(is.null(weights)) {
       weights<-rep(1.0, mraw) # set all weights to 1.0
   } 
   else if (any(weights < 0.0)) stop("weights must be non-negative")

# Subsetting -- nls() uses model frame
   if( ! missing(subset) && ! is.null(subset) ){
     # we need to subset, which we do via the weights
     if (! all(is.integer(subset))) stop("subset must have integer entries")
     if ( any(subset < 1) || any(subset > mraw) ) stop("subset entry out of range")
     #?? need to test these possibilities
     weights[- subset] <- 0.0 # NOTE: the minus gets the values NOT in the dataset
   }
   mres<-length(weights[which(weights > 0.0)]) # number of residuals (observations)
   swts<-sqrt(weights)

# Bounds and masks (fixed parameters)
   if (length(lower) == 1) lower <- rep(lower, npar) # expand to full dimension
   if (length(upper) == 1) upper <- rep(upper, npar)
   # more checks on bounds
   if (length(lower) != npar) stop("Wrong length: lower")
   if (length(upper) != npar) stop("Wrong length: upper")
   if (any(start < lower) || any(start > upper)) 
        stop("Infeasible start")

   bdmsk <- rep(1, npar)  # set all params free for now
   maskidx <- which(lower==upper)
   if (length(maskidx) > 0 && trace) {
       cat("The following parameters are masked:")
       print(pnames[maskidx])
   }
   bdmsk[maskidx] <- 0  # fixed parameters ??? do we want to define U and L parms
   bdmsk <- bdmsk

# Formula processing
   # oneSidedFormula ?? Should we be more explicit?
   if (length(formula) == 2) {
        residexpr <- formula[[2L]] ##?? Need to make sure this works 
        ## ?? -- may need to be a call
        lhs <- NULL
        rhs <- eval(formula[[2L]], envir=localdata)
   } else if (length(formula) == 3) {
         ##?? WARNING: seems to disagree with nls()
             residexpr <- call("-", formula[[3]], formula[[2]])
             lhs <- eval(formula[[2L]], envir=localdata)
             # set the "promise" for the lhs of the model
             rhs <- eval(formula[[3L]], envir=localdata)
             lnames<-all.vars(formula[[2L]]) 
             # Check that lhs is a single variable from data ??
             ldname <- which(dnames %in% lnames)
             if (length(lnames) != 1L) {
                warning("lhs has either no named variable or more than one")
             }
             else { if (control$trace) cat("lhs has just the variable ",lnames,"\n")}
           } 
           else stop("Unrecognized formula")
   if (control$derivmeth == "numericDeriv") {
         rjexpr <- residexpr # unchanged -- numeric deriv put into rjfun
   } else
        if (all(nlsderivchk(residexpr, names(start)))) { # all derivs can be computed
           rjexpr <- deriv(residexpr, names(start)) ##?? could fail on some functions
        }
        else  rjexpr <- NULL
   if (is.null(rjexpr) && (control$derivmeth == "default")) {
        warning("Changing to alternative derivative method")
        control$derivmeth <- nlsj.control()$altderivmeth
   }

# Define functions for residual and (residual + jacobian)
   resfun <- function(prm) { # only computes the residuals (unweighted)
      if (is.null(names(prm))) names(prm) <- names(start)
      localdata <- list2env(as.list(prm), parent = data)
      eval(residexpr, envir = localdata) # ?? needed? or is eval(residexpr) enough?
   }

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

##   resid <- function(prm) {- swts * rjfun(prm)} # used to return in m ?? CHECK SIGN
   ## define for m$resid in return section

##   deviance <- function(prm) { sum(resid(prm)^2)}
   # ?? nls doesn't have prm in this function

## tracefn -- called trace in nls(), but trace is also logical argument in call
   tracefn = function() {
      d <- getOption("digits")
  # Note that we assume convInfo is available
#      cat(sprintf("%-*s (%.2e): par = (%s)\n", d+4L+2L*(control$scaleOffset > 0),
# ??     format(dev, digits=d, flag="#"),
#      format(ssmin, digits=d, flag="#"),
#      convInfo$ctol,
#      paste(vapply(getPars(), format, ""), collapse=" ")))
    cat(ssmin,":(")
    for (ii in 1:npar) cat(prm[ii]," ")
    cat(")  rofftest=",attr(convInfo,"ctol"),"\n")
   } # end tracefn()


# All the following are defined at START of problem
   cjmsg <- paste("Max. jacobian evaluations (",control$maxiter,") exceeded")
   crmsg <- paste("Max. residual evaluations (",control$maxres,") exceeded")
   csmsg <- paste("Small wtd sumsquares (deviance) <= ",epstol4)
   # above could depend on the problem
   comsg <- paste("Relative offset less than ",control$tol)
   cnmsg <- "No parameters for this problem"
   cxmsg <- "Not yet defined"
   cvmsg <- c(cjmsg, crmsg, csmsg, comsg, cnmsg, cxmsg)

   convCrit <- function() { # returns TRUE if should terminate algorithm
        # put the trace function here and don't include in m
        cval <- FALSE # Initially NOT converged
        # max jacs
        cj <- (njac > control$maxiter)
        # max fns
        cr <- (nres > control$resmax)
        # small ss
        cs <- (control$smallsstest && (ssmin <= epstol4))
        # roffset < tolerance for relative offset test
        co <- FALSE # ?? for the moment
        ctol <- NA
        ## Make it always available via nlsj()
        scoff <- control$scaleOffset
        if (scoff) scoff <- (mres - npar) * scoff^2 # adjust for problem 
        ## at this point, QRJ and wres should be defined
        if (haveQRJ) {
           rr <- qr.qty(QRJ, wres)
           print(rr[1L:npar])
           print(rr[-(1L:npar)])
           ctol <- sqrt( sum(rr[1L:npar]^2) / (scoff + sum(rr[-(1L:npar)]^2)))
           ##            projected resids ss             deviance
        } else ctol <- .Machine$double.xmax ## We don't have QRJ defined, so big value
        co <- (ctol <= control$tol) # compare relative offset criterion
        # other things??
        cn <- (npar < 1) # no parameters
        cx <- FALSE # Anything else to be added
        if(npar == 0) { # define to avoid exceptions
            cval <- TRUE
            ctol <- NA
        }
        cvec <- c(cj, cr, cs, co, cn, cx)
        cmsg <- "Termination msg: "
        for (i in 1:length(cvec)){
            if (cvec[i]) cmsg <- paste(cmsg,cvmsg[i],"&&")
        }
        cval <- any(cvec)
        attr(cval, "cmsg") <- cmsg
        attr(cval, "ctol") <- ctol
        attr(cval, "nres") <- nres
        attr(cval, "njac") <- njac
        cval
   } # end convCrit()

<<<<<<< HEAD
# Initialization of iteration
# counts of evaluations
   xcmsg <- NULL
   haveQRJ <- FALSE # QR of Jacobian NOT available
   keepgoing <- TRUE # use to control the main loop
   haveJ <- FALSE # just in case -- inform program J NOT available
   nres <- 0 # Count of residual evaluations
   njac <- 0 # Count of jacobian evaluations ("iterations")
   while (keepgoing) { # Top main loop 
      if (! haveJ) resraw <- rjfun(start) # includes Jacobian
      njac <- njac + 1 
      nres <- nres + 1 # ?? be nice to ONLY compute Jacobian
      wres <- swts * resraw 
      J <- swts * attr(resraw,"gradient")  
      # NOTE: multiplication for wres does NOT chg attribute
      haveJ <- TRUE # to record whether a current Jacobian is available
      QRJ <-qr(J)
      haveQRJ <- TRUE # ?? may not need all these later, but for now ...
      qrDim <- min(dim(QRJ$qr))
      if ((algorithm == "default") && (QRJ$rank < qrDim)) stop("Singular jacobian") 
      # ?? for now don't continue with Gauss-Newton, since delta can't be computed
      ssnew <- sum(wres^2) # get the sum of squares (this is weighted)
      if (nres == 1) { # Set the starting sumsquares. Do we want to keep this??
         ssmin <- ssnew
         ssnew <- .Machine$double.xmax # To ensure we continue
      }
      convInfo <- convCrit() # This is essentially the convergence test
      if (convInfo) {
         keepgoing <- FALSE
         break # to escape the main loop
      }
      if (trace) tracefn() # printout of tracking information
      # ?? "default" algorithm
      delta <- qr.coef(QRJ, -wres) # LS solve of J delta ~= -wres
      cat("delta:"); print(as.numeric(delta))
      fac <- 1.0
      ssnew<-ssmin # initialize so we do one loop at least
      eq <- FALSE # In case it is needed below to check parameters changed
      while ((ssnew >= ssmin)  && (fac > control$minFactor)) {
         newp <- prm + fac * delta
         fac <- 0.5 * fac # ?? this is fixed in nls(), but we could alter
         #  cat("newp:"); print(as.numeric(newp))
         eq <- all( (prm+control$offset) == (newp+control$offset) )
         # We check if the parameters have been changed (eq TRUE) 
         if (! eq ) { # parameters have changed, we can try to find new sumsquares
             # ?? trying to get NEW values here. 
             # ?? Does not want to re-evaluate when running!
             newwres <- swts*resfun(newp) # ?? assume is evaluated
             ## ?? should we introduce a non-compute flag?
             nres <- nres + 1 # Only residual evaluated
             ssnew <- sum(newwres^2) 
             if (trace) cat("fac=",fac,"   ssnew=",ssnew,"\n")
             # ?? be nice to have multi-level trace
             if ( ssnew < ssmin) break
         }
         else {
=======

# Initialization
   njac <- 1 # number of jacobians so far
   nres <- 1
   resraw <- rjfun(start) # includes Jacobian
   wres <- swts * resraw 
   J <- swts * attr(resraw,"gradient")  # multiplication for wres does NOT chg attribute
   QRJ <-qr(J)
   qrDim <- min(dim(QRJ$qr))
   if (QRJ$rank < qrDim) stop("Singular initial jacobian") # for now don't continue
   ssmin <- sum(wres^2) # get the sum of squares (this is weighted)
   cat("Before iteration, deviance=",ssmin,"\n")
   while (! (convInfo <- convCrit()) ) { # Top main loop -- save convInfo at same time
       delta <- qr.coef(QRJ, -wres) # LS solve of J delta ~= -wres
 #      cat("delta:"); print(delta)
       fac <- 1.0
       ssnew<-ssmin # initialize so we do one loop at least
       while ((ssnew >= ssmin)  && (fac > control$minFactor)) {
# with bounds and masks will modify this

           newp <- prm + fac * delta




           fac <- 0.5 * fac # ?? this is fixed in nls(), but we could alter
#           cat("newp:"); print(as.numeric(newp))
           eq <- all( (prm+control$offset) == (newp+control$offset) )
           if (! eq ) {
              # ?? trying to get NEW values here. Does not want to re-evaluate when running!
              # ?? Is it using p1 and p2 rather than prm??
              newwres <- swts*resfun(newp)
              nres <- nres + 1 # Only residual evaluated
              ssnew <- sum(newwres^2) 
              if (trace) cat("fac=",fac,"   ssnew=",ssnew,"\n")
              if ( ssnew < ssmin) break
           }
           else {
>>>>>>> e0b838166b1aca7c461bedf6517bc6c0081e185b
              cat("Parameters unchanged\n")
              break
         }
      } # end inner while
      if (is.na(ssnew)) stop("ss is NA -- parameters unchanged") #?? fix
      if (trace) cat("after inner while loop over fac, ssnew=",ssnew," fac=",fac,"\n")
      if (eq || (ssnew >= ssmin)) { # no progress in "default". Done!
         keepgoing<-FALSE # No progress achieved ?? message or indicator?
         xcmsg <- "Default backtrack search failed"
      } else {         
	       prm <- newp
       	 ssmin <- ssnew
         resraw <- rjfun(prm) # ?? new res and gradient -- does extra res??
         haveJ <- TRUE
         njac <- njac + 1
         nres <- nres + 1 # Both residual and jacobian are evaluated
         if (control$watch) tmp <- readline("next iteration")
       }
#      if (trace) cat("Here report progress\n")
       fac <- 2.0 * fac # to reset for next iteration
    } # end outer while
    ## names(prm) <- pnames # Make sure names re-attached. ??Is this needed??
    if (! is.null(xcmsg)) cmsg <- paste(convInfo$cmsg,"&&",xcmsg) # include extra info
    m <- list(resfun = function(prm) resfun(prm), # ??
             resid = function() {- wres}, # ?? weighted. NOTE SIGN?? ??callable?
             rjfun = function(prm) rjfun(prm), # ??
	     fitted = function() rhs, # OK
	     formula = function() formula, #OK
	     deviance = function() ssmin, # ?? Probably wrong -- needs to be callable
	     lhs = function() lhs, #OK
	     conv = function() convCrit(), # possibly OK
	     getPars = function() {prm},
	     Rmat = function() qr.R(QRJ), # ?? we need environment 
	     predict = function(newdata = list(), qr = FALSE)
                 eval(formula[[3L]], as.list(newdata)) 
                 # ?? do we need to specify environment
	   )
    class(m) <- "nlsModel"
    result <- list(m=m, convInfo=convInfo, control=control)
    ##?? Add call -- need to set up properly somehow. Do we need model.frame?
    class(result) <- "nlsj" ## CAUSES ERRORS ?? Does it?? 190821
    result
}
