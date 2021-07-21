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

#   bdmsk <- rep(1, npar)  # set all params free for now -- done each iteration anyway
   maskidx <- which(lower==upper) # ?? may want to put in tolerance??
   if (length(maskidx) > 0 && trace) {
       cat("The following parameters are masked:")
       print(pnames[maskidx])
   }
##   bdmsk[maskidx] <- 0  # fixed parameters ??? do we want to define U and L parms
  
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

# Initialization of iteration
# counts of evaluations
   xcmsg <- NULL
   haveQRJ <- FALSE # QR of Jacobian NOT available
   keepgoing <- TRUE # use to control the main loop
   haveJ <- FALSE # just in case -- inform program J NOT available
   nres <- 0 # Count of residual evaluations
   njac <- 0 # Count of jacobian evaluations ("iterations")
   while (keepgoing) { # Top main loop 
      if (! haveJ) { # need to get new Jacobian and reset bounds constraints
         resraw <- rjfun(start) # includes Jacobian
         njac <- njac + 1 
         nres <- nres + 1 # ?? be nice to ONLY compute Jacobian
      } # Probably ALWAYS recompute J here??
      wres <- swts * resraw 
      J <- swts * attr(resraw,"gradient")  
      # NOTE: multiplication for wres does NOT chg attribute
      haveJ <- TRUE # to record whether a current Jacobian is available
      bdmsk<-rep(1,npar)
      bdmsk[maskidx]<-0 # masked
      bdmsk[which(prm-lower<epstol*(abs(lower)+epstol))]<- -3 # at lower bounds
      bdmsk[which(upper-prm<epstol*(abs(upper)+epstol))]<- -1 # at upper bounds
      # Here we use a tolerance to see if we are close to bounds
      QRJ <-qr(J)
      haveQRJ <- TRUE # ?? may not need all these later, but for now ...
      qrDim <- min(dim(QRJ$qr))
      if ((algorithm == "default") && (QRJ$rank < qrDim)) stop("Singular jacobian") 
      # ?? for now don't continue with Gauss-Newton, since delta can't be computed
      ssnew <- sum(wres^2) # get the sum of squares (this is weighted)
      if (nres == 1) { # Set the starting sumsquares on first iteration
         ssmin <- ssnew # the best ss (prm are parameters)
         ssnew <- .Machine$double.xmax # To ensure we continue
      }
      convInfo <- convCrit() # This is essentially the convergence test
      if (trace) tracefn() # printout of tracking information
      if (! convInfo) {??wrong
         keepgoing <- FALSE
         break # to escape the main loop
      }
      # "default" algorithm -- added "try"
      delta <- try(qr.coef(QRJ, -wres)) # LS solve of J delta ~= -wres
      if (inherits(delta,"try-error")) stop("Cannot solve Gauss-Newton equations")
      cat("delta:"); print(as.numeric(delta))
      # Do we need the gradient projection??
      gjty <- t(J) %*% wres
      cat("bdmsk:"); print(bdmsk)
      for (i in 1:npar){ #?? Improve using subsets
         bmi<-bdmsk[i]
         if (bmi==0) {
             gjty[i]<-0 # masked
             J[,i]<-0
          }
# Following code works by using sign = 2+bdmsk = 1 for upper bd, -1 for lower bd
# If gradient is negative (downhill) at LB, sign*gradient is > 0. OK if delta > 0
# If gradient is positive (uphill) at UB, we can step backwards, 
#    sign * gradient > 0, and we can proceed if delta < 0 
          if (bmi<0) {
             if((2+bmi)*gjty[i] > 0) { # free parameter
                bdmsk[i]<-1
                if (control$watch) cat("freeing parameter ",i," now at ",prm[i],"\n")
             } else {
                gjty[i]<-0 # active bound
                J[,i]<-0
                if (control$watch) cat("active bound ",i," at ",prm[i],"\n") 
             }
          } # bmi
      } # end for loop on i
      gproj <- crossprod(delta, gjty)
      if (is.na(gproj) || (gproj >= 0) ) {
         if (trace) cat("Uphill step direction") # should NOT be possible, but at end??
         keepgoing<-FALSE #?? may want cleaner exit
         xcmsg <- "Uphill search direction"
         break
      }
      gangle <- gproj/sqrt(sum(gjty^2) * sum(delta^2))
      gangle <- 180 * acos(sign(gangle)*min(1, abs(gangle)))/pi
      if (control$watch) cat("gradient projection = ",gproj,
                      " g-delta-angle=",gangle,"\n")
      step<-rep(1,npar)  # Check how far to bounds
      cat("chk bds, pmr:");print(as.numeric(prm))
      cat("       upper:");print(as.numeric(upper))
      cat("       lower:");print(as.numeric(lower))
      for (i in 1:npar){
          bd<-bdmsk[i]
          da<-delta[i]
          if (control$watch) cat(i," bdmsk=",bd,"  delta=",da,"\n")
          if (bd==0 || ((bd==-3) && (da<0)) ||((bd==-1) && (da>0))) {
             delta[i]<-0 # Cases where bounds or masks active
          } else {
             if ((da > 0) && (is.finite(upper[i]))) step[i]<-(upper[i]-prm[i])/da
             if ((da < 0) && (is.finite(lower[i]))) step[i]<-(lower[i]-prm[i])/da
             # positive steps in both cases
          }
      } # end loop (?? can we make more efficient?
      fac <- min(1.0, step[which(delta!=0)]) # stepsize control
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
              cat("Parameters unchanged\n")
              xcmsg <- "Parameters unchanged"
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
    cat("outside outer while loop\n")
    if (! is.null(xcmsg)) cmsg <- paste(attr(convInfo,"cmsg"),"&&",xcmsg) # include extra info
    cat("build m\n")
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
