nlsj <- function (formula, data = parent.frame(), start, control = nlsj.control(),
            algorithm = "default", weights=NULL, subset=NULL, trace = FALSE,
            na.action, model=FALSE, lower = -Inf, upper = Inf, ...) {
# ?? left out -- FIXME??    na.action, model = FALSE, (masked from nlxb)
# ?? at this stage ONLY treat "default", but will add bounds
# ?? data in .GlobalEnv -- should be OK
##?? Should a lot of this material be in nlsjModel() to build tools for problem??
# DOES NOT CALL nlsjModel, but does everything here

##?? ... args may NOT be well-defined for this function. CAUTION  
  
# Controls
   if (is.null(control$derivmeth)) control$derivmeth="default" # for safety
   epstol <- (.Machine$double.eps * control$offset) 
   epsh <- sqrt(epstol)
   epstol4 <- epstol^4 # used for smallsstest
   ##?? may want these in nlsj.control
   if (control$derivmeth == "numericDeriv") warning("Forcing numericDeriv")
#   cat("control:"); print(control)
#   tmp <- readline("cont.")
#   control$watch<-TRUE

# Algorithm
   if (is.null(algorithm)) algorithm<-"default"
   algchoice <- c("default", "port", "plinear", "marquardt")
   alg <- which(algorithm %in% algchoice)
   switch(alg,
      "default" = if (trace) cat("nlsj: Using default algorithm\n"),
      "marquardt" = if (trace) cat("nlsj: Using marquardt algorithm\n"),
      # all other choices
      { msg <- paste("Algorithm choice '",alg,
                     "' not yet implemented or does not exist")
        stop(msg) }
   ) # end switch choices

getlen <- function(lnames) {
      nn<-length(lnames)
      ll<-rep(NA,nn)
      for (i in 1:nn) ll[i]=length(get(lnames[i]))
      ll
}
   
# Data
   stopifnot(inherits(formula, "formula"))
   vnames <- all.vars(formula) # all names in the formula
   if (missing(data)) {  ## rather than   #  if (is.null(data)) {
      warning("Data is not declared explicitly. Caution!")
      data <- environment(formula) # this will handle variables in the parent frame??
      dnames <- vnames[which(vnames %in% ls(data))]
      ldata<-getlen(vnames) # lengths of data
      dnames <- vnames[which(ldata > 1)]
      pnames <- vnames[which(ldata == 1)] # could 
   }
   else if (is.list(data)){
          data <- list2env(data, parent = environment(formula))
          dnames <- vnames[which(vnames %in% ls(data))]
          if (length(dnames) < 1) stop("No data found")
          pnames <- vnames[ - which(vnames %in% dnames)] 
          # the "non-data" names in the formula
        }
        else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")

   npar <- length(pnames)
   if (npar < 1) stop("No parameters!")

# Start vector
   if (is.null(start)) { # start not specified
     warning("start vector not specified for nlsj")
     start<-0.9+seq(npar)*0.1 # WARNING: very crude??
     names(start)<-pnames # and make sure these are named?? necessary??
     ## ??? put in ways to get at selfstart models 
   }
   else { # we have a start vector
     snames<-names(start) # names in the start vector
     start <- as.numeric(start) # ensure we convert (e.g., if matrix)
     names(start) <- snames ## as.numeric strips names, so this is needed ??
   }
   prm <- start # start MUST be defined at this point
   localdata <- list2env(as.list(prm), parent = data)

   #  model frame
   # if (model) {
   #   cat("model frame args:\n")
   #   args<-list(formula=formula, data=data, subset=subset, ...)
   #   # We include ... for completeness, but may be undefined or wrong
   #   print(str(args))
   #   tmp<-readline("continue.")
   #   mf <- do.call(stats::model.frame, args)
   #   print(str(mf))
   #   tmp<-readline("more.")
   # } else 
        mf <- NULL
   
   # Weights
   mraw <- length(eval(as.name(dnames[1]), envir=data))
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
     weights[- subset] <- 0.0 # NOTE: the minus sgets the values NOT in the dataset
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

   maskidx <- which(lower==upper) # ?? may want to put in tolerance??
   if (length(maskidx) > 0 && trace) {
       cat("The following parameters are masked:")
       print(pnames[maskidx])
   } ## NOT NEEDED else {cat("maskidx:");print(maskidx)}
     
  
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
   rhs <- eval(rhsexpr, envir=localdata)

   lhs <- eval(lhsexpr, envir=localdata)

# Define functions for residual and (residual + jacobian)
# Currently do not use resfun -- just rjfun
#   resfun <- function(prm) { # only computes the residuals (unweighted)
#      if (is.null(names(prm))) names(prm) <- names(start)
#      localdata <- list2env(as.list(prm), parent = data)
#      eval(residexpr, envir = localdata) # ?? needed? or is eval(residexpr) enough?
#   }

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
        if (scoff) scoff <- (mres - npar) * scoff^2 # adjust for problem 
        ## at this point, QRJ and wresy should be defined
        if (haveQRJ) {
           rr <- qr.qty(QRJ, - wresy) # ?? + or -?
        #   print(rr[1L:npar])
        #   print(rr[-(1L:npar)])
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
   marqalg<-FALSE # Is there a better way??
   slam <- 0.0 #?? to avoid trouble, cleanup later
   if (algorithm=="default") defalg<-TRUE else defalg<-FALSE
   if (algorithm=="marquardt") {
      marqalg<-TRUE
      slam <- sqrt(control$lamda)
#      cat("initial slam=",slam,"\n")
#??NOT needed as is included later
      if (slam <= epsh) slam <- epsh*10. # ensure NOT too small
      slinc <- sqrt(control$laminc)
      sldec <- sqrt(control$lamdec)
   }
# counts of evaluations
   xcmsg <- NULL
   haveQRJ <- FALSE # QR of Jacobian NOT available
   keepgoing <- TRUE # use to control the main loop
   haveJ <- FALSE # J NOT pulled from residual object
   nres <- 1 # Count of residual evaluations
   njac <- 1 # Count of jacobian evaluations ("iterations")
   resb <- rjfun(start) # "best" so far
   wresb <- swts* resb # as.numeric takes twice as long?!
   # NOTE: multiplication for wresb does NOT chg attribute
   prm <- start
#   cat("npar=",npar," pnames:"); print(pnames)

   ssnew <- sum(wresb^2) # get the sum of squares (this is weighted)
   ssmin <- ssnew # the best ss (prm are parameters)
   fac <- 1.0 # to ensure initially defined
#   cat("Top - slam=",slam," ssmin=",ssmin," at "); print(as.numeric(prm))
   # ?? Do we want to record ss0, res0 ?
#   cat("npar=",npar,"\n")
   bdmsk<-rep(1,npar)
#   print(bdmsk)
   while (keepgoing) { # Top main loop 
      if (! haveJ) { # need to get new Jacobian and reset bounds constraints
         J <- swts * attr(resb,"gradient")
         njac <- njac + 1
         haveJ <- TRUE # to record whether a current Jacobian is available
         # Need to check bounds ONLY when new jacobian (and new gradient)
         if (length(bdmsk) > 0) bdmsk[maskidx]<-0 # masked
         bdmsk[which(prm-lower<epstol*(abs(lower)+epstol))]<- -3 # at lower bounds
         bdmsk[which(upper-prm<epstol*(abs(upper)+epstol))]<- -1 # at upper bounds
         # Here we use a tolerance to see if we are close to bounds
         J0 <- J # save raw Jacobian
         wresy <- wresb # needed for default
      }
      gjty <- t(J0) %*% wresb   # Need gradient projection for bounds
      #?? Should this be wresb or wresy and J or J0. Probably J0 (original)
#      print(bdmsk)
      for (i in 1:npar){ # Tried subsets but slower
        bmi<-bdmsk[i]
#        cat("i=",i," bmi=",bmi,"\n")
        if (bmi==0) {
           gjty[i]<-0 # masked
           J0[,i]<-0
        }
# Following code works by using sign = 2+bdmsk = 1 for upper bd, -1 for lower bd
# If gradient is negative (downhill) at LB, sign*gradient is > 0. OK if delta > 0
# If gradient is positive (uphill) at UB, we can step backwards, 
#   sign * gradient > 0, and we can proceed if delta < 0 
        if (bmi<0) {
          if((2+bmi)*gjty[i] > 0) { # free parameter
             bdmsk[i]<-1
             if (control$watch) cat("freeing parameter ",i," now at ",prm[i],"\n")
          } else {
             gjty[i]<-0 # active bound
             J0[,i]<-0
             if (control$watch) cat("active bound ",i," at ",prm[i],"\n") 
          }
        } # bmi < 0. bmi > 0 is free parameter, so no fuss
      } # end for loop on i
      if (marqalg) { # This is whether or not new Jacobian
	 if (slam <= epsh) slam <- epsh*10. # ensure NOT too small
         J <- rbind(J0, diag(slam, npar)) # augment the Jacobian
         wresy <- c(wresb, rep(0,npar)) # and the residuals (y)
         # ?? This does NOT use any scaling. To be considered??
         # Could mean(abs(diag(QRJ$qr))) be used to scale epsh as slam?
         haveQRJ <- FALSE # since J has changed
         cat("Adjusted J0 to J, haveQRJ=",haveQRJ,"\n")
      }
#      cat("haveQRJ:",haveQRJ," J:"); print(J)
#      cat("wresb:"); print(as.numeric(wresb))
      if (! haveQRJ) QRJ <-qr(J)
      haveQRJ <- TRUE # Since we will have QRJ for sure here
      singJ <- FALSE # singular Jacobian? ?? may not keep this indicator
      qrDim <- min(dim(QRJ$qr))
      if (QRJ$rank < qrDim) singJ <- TRUE
      if (trace && singJ) {cat("Jacobian when rank<dim:\n"); print(J)}
      if (defalg && singJ) stop("Singular jacobian")
      # ?? Don't continue with Gauss-Newton; delta can't be computed
      # ?? marquardt is using J (augmented Jacobian) in convergence test. This may
      # ?? not be exactly what we want. Leave it there for now.
      convInfo <- convCrit() # This is the convergence test
      if (control$watch) print(convInfo)
      if (marqalg && trace) cat("slam =",slam," ")
      if (trace) tracefn() # printout of tracking information <<<<<<<<<<<< printout
      if (convInfo) { # Stop if we have converged or must terminate
         if (control$watch) cat("convInfo TRUE -- stop iteration\n")
         keepgoing <- FALSE
         break # to escape the main loop
      }
      if (! singJ) {
        delta <- try(qr.coef(QRJ, -wresy)) # LS solve of J delta ~= -wres
        deltaOK <- TRUE
        if (inherits(delta,"try-error")) {
          if (! marqalg) stop("Cannot solve Gauss-Newton equations")
          # ?? again consider a switch to marquardt
          deltaOK <- FALSE # Tells program we want to increase lambda / slam
        }
      }
      else deltaOK <- FALSE # delta NOT OK if singJ
      if (deltaOK) {
         gproj <- as.numeric(crossprod(delta, gjty)) # ?? do we need as.numeric
         if (is.na(gproj) || (gproj >= 0) ) {
           if (trace) cat("Uphill or NA step direction\n") 
           # should NOT be possible, except possibly when converged??
           # keepgoing<-FALSE #?? may want cleaner exit
           # ??? do we want to increase slam if marquardt?
           xcmsg <- "Uphill search direction"
##           tmp <- readline("uphill search")
##           break #??  See what happens if we just continue
         }
         gangle <- gproj/sqrt(sum(gjty^2) * sum(delta^2))
         gangle <- 180 * acos(sign(gangle)*min(1, abs(gangle)))/pi
         if (control$watch) cat("gradient projection = ",gproj,
                      " g-delta-angle=",gangle,"\n")
         step<-rep(1,npar)  # Check how far to bounds
         print(delta)
         for (i in 1:npar){
            bd<-bdmsk[i]
            da<-delta[i]
#            if (control$watch) 
            cat(i," bdmsk=",bd,"  delta=",da,"\n")
            if (bd==0 || ((bd==-3) && (da<0)) ||((bd==-1) && (da>0))) {
               delta[i]<-0 # Cases where bounds or masks active
            } else {
               if ((da > 0) && (is.finite(upper[i]))) step[i]<-(upper[i]-prm[i])/da
               if ((da < 0) && (is.finite(lower[i]))) step[i]<-(lower[i]-prm[i])/da
               # positive steps in both cases
            }
         } # end loop
         eq <- FALSE # In case it is needed below to check parameters changed
         cat("fac before =",fac)
         if (defalg) fac <- min(fac, step[which(delta!=0)]) # stepsize control
         else fac <- min(1.0, step[which(delta!=0)]) # stepsize control
         cat("  after=",fac,"\n")
         if (control$watch) {cat("stepsize=",fac," delta:"); print(as.numeric(delta))}
         while ((ssnew >= ssmin)  && (fac > control$minFactor)) {
           cat("fac=",fac," ssnew=")
           newp <- prm + fac * delta
           #  cat("newp:"); print(as.numeric(newp))
           eq <- all( (prm+control$offset) == (newp+control$offset) )
           # We check if the parameters have been changed (eq TRUE) 
           if (! eq ) { # parameters have changed, we can try to find new sumsquares
             res <- rjfun(newp) # ?? assume is evaluated
             ## ?? Currently only use rjfun -- then don't need to call again for J
             nres <- nres + 1 # but count just the residual use
             ## Put at least one NA in res if "not computable", but this may not be
             ## something the user can control in detail
             if (any(is.na(res))) {
                if (marqalg) { # note: could have both defalg and marqalg true if we
                ## allow switch when singular gradient ??
                   slam <- slam*slinc
                   warning("NA in residuals")
                   break # out of inner loop
                }
                ssnew <- .Machine$double.max # ensure BIG
             } 
             else { # Do NOT recompute wresb -- stays until new J used
                ssnew <- sum((swts * res)^2) 
             }   
             cat(ssnew,"\n")
             if (control$watch) { 
                 cat("fac=",fac,"   ssnew=",ssnew,"  ",(ssnew<ssmin)," "); 
                 print(as.numeric(newp))
             } # ?? be nice to have multi-level trace
             if (ssnew < ssmin) {
                if (trace) cat("<")
                fac <- min(2.0*fac, 1.0) # reset for next iteration to match nls()
                break # finished inner loop in all algs
             }
             else if (marqalg) {
                 slam <- slam*slinc
                 break # out of inner loop
             }
             # not marqalg, assume backtrack continues with new fac
           } # end ! eq (parameters changed)
           else { ## eq TRUE
              if (trace) cat("Parameters unchanged\n")
              xcmsg <- "Parameters unchanged"
              keepgoing <- FALSE # all methods??
              tmp <- readline("Unchanged")
              break
           }
           fac <- 0.5 * fac # ?? this is fixed in nls(), but we could alter
         } # end inner while
         if (keepgoing && (ssnew < ssmin)) {# found better point
             if (defalg && control$watch) 
                  cat("Backtrack: ssnew=",ssnew," fac=",fac,"\n")
             prm <- newp
             ssmin <- ssnew
             resb<-res
             wresb <- swts*resb # ??can we rationalize somehow
             haveJ <- FALSE # need new Jacobian
             haveQRJ <- FALSE 
             if (marqalg) slam <- sldec * slam # reduce lamda in Marquardt
             if (control$watch) tmp <- readline("next iteration")
         } # don't need to use else
      } # deltaOK
      else { # delta NOT OK -- report
         # warning("delta NOT computed OK")
	 cat("delta NOT computed OK, haveJ=",haveJ,"\n")
         slam <- slinc * slam # will have stopped if not marqalg
         tmp <- readline("cont.")
      }       
      # at this point, we have either made progress (new ssmin)
      # or we are going round again
   } # end outer while
   if (! is.null(xcmsg)) cmsg <- paste(attr(convInfo,"cmsg"),"&&",xcmsg) # include extra info
#   cat("build m\n")
# Ensure reset of values of parameters.  Needed!
   localdata <- list2env(as.list(prm), parent = data)
   rhs <- eval(rhsexpr, envir=localdata)
   lhs <- eval(lhsexpr, envir=localdata)

   m <- list(resfun = function(prm) as.numeric(rjfun(prm)), # to get just resdiduals
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
    # ?? as of 20210802, missing na.action, call, convergence, message, dataClasses
    result <- list(m=m, convInfo=convInfo, weights=weights, model=mf, control=control)
    ##?? Add call -- need to set up properly somehow. Do we need model.frame?
    class(result) <- "nlsj" ## CAUSES ERRORS ?? Does it?? 190821
    result
}
