nlsgn <- function (formula, data = parent.frame(), start, control = nlsj.control(),
            weights=NULL, subset, trace = FALSE) {
# ?? left out -- FIXME??    na.action, model = FALSE, (masked from nlxb)
# Controls
   if (is.null(control$derivmeth)) control$derivmeth="default" # for safety
   epstol <- (.Machine$double.eps * control$offset) 
   epsh <- sqrt(epstol)
   epstol4 <- epstol^4 # used for smallsstest
   if (control$derivmeth == "numericDeriv") warning("Forcing numericDeriv")

getlen <- function(lnames) {
      nn<-length(lnames)
      ll<-rep(NA,nn)
      for (i in 1:nn) ll[i]=length(get(lnames[i]))
      ll
}
   
# Data
   stopifnot(inherits(formula, "formula"))
   vnames <- all.vars(formula) # all names in the formula
   if (missing(data)) {  
      warning("Data is not declared explicitly. Caution!")
      data <- environment(formula) # this will handle variables in the parent frame??
      dnames <- vnames[which(vnames %in% ls(data))]
      ldata<-getlen(vnames) # lengths of data
      dnames <- vnames[which(ldata > 1)]
      pnames <- vnames[which(ldata == 1)] # could 
   } else if (is.list(data)) {
          data <- list2env(data, parent = environment(formula))
          dnames <- vnames[which(vnames %in% ls(data))]
          if (length(dnames) < 1) stop("No data found")
          pnames <- vnames[ - which(vnames %in% dnames)] 
          # the "non-data" names in the formula
        } else if (!is.environment(data)) stop("'data' must be a dataframe, list, or environment")

   npar <- length(pnames)
   if (npar < 1) stop("No parameters!")

# Start vector
   if (is.null(start)) { # start not specified
     warning("start vector not specified for nlsj")
     start<-0.9+seq(npar)*0.1 # WARNING: very crude??
     names(start)<-pnames # and make sure these are named?? necessary??
     ## ??? put in ways to get at selfstart models 
   } else { # we have a start vector
     snames<-names(start) # names in the start vector
     start <- as.numeric(start) # ensure we convert (e.g., if matrix)
     names(start) <- snames ## as.numeric strips names, so this is needed ??
   }
   prm <- start # start MUST be defined at this point
   localdata <- list2env(as.list(prm), parent = data)

# Weights
   mraw <- length(eval(as.name(dnames[1]), envir=data))
   if(is.null(weights)) {
       weights<-rep(1.0, mraw) # set all weights to 1.0
   } else if (any(weights < 0.0)) stop("weights must be non-negative")

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

# Formula processing ?? Do we need all these? Can we simplify?f
   # oneSidedFormula ?? Should we be more explicit?
   if (length(formula) == 2) {
        residexpr <- formula[[2L]] ##?? Need to make sure this works 
        ## ?? -- may need to be a call
        lhsexpr <- NULL
        rhsexpr <- eval(formula[[2L]], envir=localdata)
   } else if (length(formula) == 3) {
         ##?? WARNING: seems to disagree with nls()
             residexpr <- call("-", formula[[3]], formula[[2]])
             lhsexpr <- formula[[2L]]
             
             # set the "promise" for the lhs of the model
             rhsexpr <- formula[[3L]]
             lnames<-all.vars(formula[[2L]]) 
             # Check that lhs is a single variable from data ??
             ldname <- which(dnames %in% lnames)
             if (length(lnames) != 1L) {
                warning("lhs has either no named variable or more than one")
             }
             else { if (trace) cat("lhs has just the variable ",lnames,"\n")}
           } else stop("Unrecognized formula")
   if (control$derivmeth == "numericDeriv") {
         rjexpr <- residexpr # unchanged -- numeric deriv put into rjfun
   } else
        if (all(nlsderivchk(residexpr, names(start)))) { # all derivs can be computed
           rjexpr <- deriv(residexpr, names(start)) ##?? could fail on some functions
        } else  rjexpr <- NULL
   if (is.null(rjexpr) && (control$derivmeth == "default")) {
        warning("Changing to alternative derivative method")
        control$derivmeth <- nlsj.control()$altderivmeth
   }
   rhs <- eval(rhsexpr, envir=localdata)

   lhs <- eval(lhsexpr, envir=localdata)

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

## tracefn -- called trace in nls(), but trace is also logical argument in call
   tracefn = function() {
      d <- getOption("digits")
  # Note that we assume convInfo is available
    cat(" ",nres,"/",njac," ")
    cat(ssmin,":(") # ?? deviance??
    for (ii in 1:npar) cat(prm[ii]," ")
    cat(")  rofftest=",attr(convInfo,"ctol"),"\n")
   } # end tracefn()


# All the following are defined at START of problem
   cjmsg <- paste("Max. jacobian evaluations (",control$maxiter,") exceeded")
   crmsg <- paste("Max. residual evaluations (",control$maxres,") exceeded")
   csmsg <- paste("Small wtd sumsquares (deviance) <= ",epstol4)
   # above could depend on the problem
   comsg <- paste("Relative offset less than ",control$tol)
   cfmsg <- paste("Step factor less than ",control$minFactor)
   cnmsg <- "No parameters for this problem"
   cxmsg <- "Not yet defined"
   cvmsg <- c(cjmsg, crmsg, csmsg, comsg, cfmsg, cnmsg, cxmsg)

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
        if (scoff) scoff <- (mres - npar) * scoff^2 # adjust for problem ??
        ## at this point, QRJ and wresy should be defined
        if (haveQRJ) {
           rr <- qr.qty(QRJ, -wresb) # ?? + or -?
           ctol <- sqrt( sum(rr[1L:npar]^2) / (scoff + sum(rr[-(1L:npar)]^2)))
           ##            projected resids ss             deviance
        } else ctol <- .Machine$double.xmax ## We don't have QRJ defined, so big value
        co <- (ctol <= control$tol) # compare relative offset criterion
        cf <- (fac <= control$minFactor)
        
        # other things??
        cn <- (npar < 1) # no parameters
        cx <- FALSE # Anything else to be added
        if(npar == 0) { # define to avoid exceptions
            cval <- TRUE
            ctol <- NA
        }
        cvec <- c(cj, cr, cs, co, cf, cn, cx)
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
#   Algorithm controls
   slam <- 0.0 #?? to avoid trouble, cleanup later
# counts of evaluations
   xcmsg <- NULL
   haveQRJ <- FALSE # QR of Jacobian NOT available
   keepgoing <- TRUE # use to control the main loop
   haveJ <- FALSE # J NOT pulled from residual object
   nres <- 1 # Count of residual evaluations
   njac <- 1 # Count of jacobian evaluations ("iterations")
   resb <- rjfun(start) # "best" so farl
   wresb <- swts* resb # as.numeric takes twice as long?!
   # NOTE: multiplication for wresb does NOT chg attribute
   prm <- start

   ssnew <- sum(wresb^2) # get the sum of squares (this is weighted)
   ssmin <- ssnew # the best ss (prm are parameters)
   fac <- 1.0 # to ensure initially defined
   while (keepgoing) { # Top main loop 
      if (! haveJ) { # need to get new Jacobian and reset bounds constraints
         J <- swts * attr(resb,"gradient")
         njac <- njac + 1
         haveJ <- TRUE # to record whether a current Jacobian is available
         # Need to check bounds ONLY when new jacobian (and new gradient)
         J0 <- J # save raw Jacobian
         wresy <- wresb # needed for default
      }
      gjty <- t(J0) %*% wresb   # Need gradient projection for bounds
      #?? Should this be wresb or wresy and J or J0. Probably J0 (original)
      if (! haveQRJ) QRJ <-qr(J)
      haveQRJ <- TRUE # Since we will have QRJ for sure here
      qrDim <- min(dim(QRJ$qr))
      if (QRJ$rank < qrDim) {
            if (trace) print(J)
            stop("Singular jacobian")
      }
      convInfo <- convCrit() # This is the convergence test
      if (control$watch) print(convInfo)
      if (trace) tracefn() # printout of tracking information
      if (convInfo) { # Stop if we have converged or must terminate
         if (control$watch) cat("convInfo TRUE -- stop iteration\n")
         keepgoing <- FALSE
         break # to escape the main loop
      }
      delta <- try(qr.coef(QRJ, -wresy)) # LS solve of J delta ~= -wres
      deltaOK <- TRUE
      if (inherits(delta,"try-error")) {
          stop("Cannot solve Gauss-Newton equations")
          deltaOK<-FALSE # cannot get here!!
      }
      if (deltaOK) {
         gproj <- as.numeric(crossprod(delta, gjty)) # ?? do we need as.numeric
         if (is.na(gproj) || (gproj >= 0) ) {
           if (trace) cat("Uphill or NA step direction")
           # should NOT be possible, except possibly when converged??
#           keepgoing<-FALSE #?? may want cleaner exit
           xcmsg <- "Uphill search direction"
#           tmp <- readline("uphill search")
          ## ?? break
         }
         eq <- FALSE # In case it is needed below to check parameters changed
         if (control$watch) {cat("stepsize=",fac," delta:"); print(as.numeric(delta))}
         while ((ssnew >= ssmin)  && (fac > control$minFactor)) {
           newp <- prm + fac * delta
           eq <- all( (prm+control$offset) == (newp+control$offset) )
           # We check if the parameters have been changed (eq TRUE)
           if (! eq ) { # parameters have changed, we can try to find new sumsquares
             res <- rjfun(newp) # ?? assume is evaluated
             nres <- nres + 1 # but count just the residual use
             if (any(is.na(res))) {
                   warning("NA in residuals")
                   ssnew <- .Machine$double.max # ensure BIG
             }
             else { # Do NOT recompute wresb -- stays until new J used
                ssnew <- sum((swts * res)^2)
             }
             if (control$watch) {
                 cat("fac=",fac,"   ssnew=",ssnew,"  ",(ssnew<ssmin)," ");
                 print(as.numeric(newp))
             } # ?? be nice to have multi-level trace
             if (ssnew <= ssmin) { # nls() uses <=. JN prefers <
                fac <- min(1.0, 2.0*fac)
                if (trace) cat("<")
                break # finished inner loop in all algs
             } else { fac <- 0.5*fac }
             # assume backtrack continues with new fac
           } # end ! eq (parameters changed)
           else { ## eq TRUE
              if (trace) cat("Parameters unchanged\n")
              xcmsg <- "Parameters unchanged"
              keepgoing <- FALSE # all methods??
#              tmp <- readline("Unchanged")
              break
           }
         } # end inner while
         if (keepgoing && (ssnew < ssmin)) {# found better point
             if (trace) cat("Backtrack: ssnew=",ssnew," fac=",fac,"\n")
             prm <- newp
             ssmin <- ssnew
             resb<-res
             wresb <- swts*resb # ??can we rationalize somehow
             haveJ <- FALSE # need new Jacobian
             haveQRJ <- FALSE #!! must reset
             if (control$watch) tmp <- readline("next iteration")
         } # don't need to use else
      } else warning("delta NOT computed OK")
      # at this point, we have either made progress (new ssmin)
      # or we are going round again
   } # end outer while
#   tmp <- readline("end of while loop in nlsj")
   if (! is.null(xcmsg)) cmsg <- paste(attr(convInfo,"cmsg"),"&&",xcmsg) # include extra info
# Ensure reset of values of parameters.  Needed!
   localdata <- list2env(as.list(prm), parent = data)
   rhs <- eval(rhsexpr, envir=localdata)
   lhs <- eval(lhsexpr, envir=localdata)

   m <- list(resfun = function(prm) as.numeric(rjfun(prm)), # ??
             resid = function() {- wresb}, #  weighted. NOTE SIGN?? ??callable?
             rjfun = function(prm) rjfun(prm), # ??
             jacobian = function() attr(resb,"gradient"),
             gradient = function() attr(resb,"gradient"),
	     fitted = function() rhs, # working??  UNWEIGHTED
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
    result <- list(m=m, convInfo=convInfo, weights=weights, control=control)
    ##?? Add call -- need to set up properly somehow. Do we need model.frame?
    class(result) <- "nlsj" ## CAUSES ERRORS ?? Does it?? 190821
    result
}
