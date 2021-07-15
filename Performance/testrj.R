# test rjfun values with numeric/analytic derivs: testrj.R
rm(list=ls())
library(nlsj)

# A function to evaluate f in environment e
# with_env <- function(f, e=parent.frame()) {
#   stopifnot(is.function(f))
#   environment(f) <- e
#   f
# }
# construct the data vectors using c()
xdata = c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata = c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)
Cform <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
Ccall<-call("-",Cform[[3]], Cform[[2]])
Cstart<-list(p1=1,p2=0.2)
Cdata<-data.frame(xdata, ydata)
Ctheta<-c("p1","p2")
mdata<-length(xdata)
Cwts <- rep(0.25, mdata)
Csubset<-1:8

formula <- Cform
start <- Cstart
data <- Cdata
subset <- NULL # to start with
weights<-NULL
trace <- TRUE
control <- nls.control()
control$trace <- trace
# 
# modj <- nlsjModel(form=Cform, data=Cdata, start=Cstart, wts=Cwts, control=nlsj.control())

## nlsjx <- function (formula, data = parent.frame(), start, control = nlsj.control(),
##                   algorithm = "default", weights=NULL, subset, trace = FALSE,
##                   lower = -Inf, upper = Inf, ...) {
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
  
  # Data
  stopifnot(inherits(formula, "formula"))
  if (is.null(data) ) {
    data <- environment(formula) # this will handle variables in the parent frame
  } else {if (is.list(data)){
    data <- list2env(data, parent = environment(formula))
  }
  else {if (!is.environment(data))
    stop("'data' must be a dataframe, list, or environment")
  }
  }
  dnames <- all.vars(formula)[which(all.vars(formula) %in% ls(data))]
  if (length(dnames) < 1) stop("No data found")
  vnames <- all.vars(formula) # all names in the formula
  pnames <- vnames[ - which(vnames %in% dnames)] # the "non-data" names in the formula
  npar <- length(pnames)
  
  # Start
  if (is.null(start)) { # start not specified
    warning("start vector not specified for nlsj")
    start<-0.9+seq(npar)*0.1 # WARNING: very crude??
    names(start)<-pnames # and make sure these are named?? necessary??
    ## ??? put in ways to get at selfstart models 
  } else { # we have a start vector
    snames<-names(start) # names in the start vector
    if ((length(snames) != length(pnames)) || (! all.equal(snames, pnames))) {
      stop("Start names differ in number or name from formula parameter names")
    }
    start <- as.numeric(start) # ensure we convert (e.g., if matrix)
    names(start) <- snames ## as.numeric strips names, so this is needed ??
  }
  prm <- start # start MUST be defined at this point
  localdata <- list2env(as.list(prm), parent = data)
  
  # Weights
  mraw <- length(eval(as.name(dnames[1])))
  if(is.null(weights)) {
    weights<-rep(1.0, mraw) # set all weights to 1.0
  } else if (any(weights < 0.0)) stop("weights must be non-negative")
  
  # Subsetting
  if( ! missing(subset) && ! is.null(subset) ){ 
    # we need to subset, which we do via the weights
    if (! all(is.integer(subset))) stop("subset must have integer entries")
    if ( any(subset < 1) || any(subset > mraw) ) stop("subset entry out of range")
    #?? need to test these possibilities
    weights[- subset] <- 0.0 # NOTE: the minus gets the values NOT in the dataset
  }
  mres<-length(weights[which(weights > 0.0)]) # number of residuals (observations)
  swts<-sqrt(weights)
  
  # Formula processing
  # oneSidedFormula ?? Should we be more explicit?
  if (length(formula) == 2) {
    residexpr <- formula[[2L]] ##?? Need to make sure this works -- may need to be call
    lhs <- NULL
    rhs <- eval(formula[[2L]], envir=localdata)
  } else if (length(formula) == 3) {
    ##?? WARNING: seems to disagree with nls()
    residexpr <- call("-", formula[[3]], formula[[2]])
    lhs <- eval(formula[[2L]], envir=localdata)
    # set the "promise" for the lhs of the model
    rhs <- eval(formula[[3L]], envir=localdata)
    lnames<-all.vars(formula[[2L]]) # Check that lhs is a single variable from data
    ldname <- which(dnames %in% lnames)
    if (length(lnames) != 1L) {
      warning("lhs has either no named variable or more than one")
    }
    else { if (control$trace) cat("lhs has just the variable ",lnames,"\n")}
  } else stop("Unrecognized formula")
  if (control$derivmeth == "numericDeriv") {
    rjexpr <- residexpr # unchanged -- numeric deriv put into rjfun
  } else
    if (all(nlsderivchk(residexpr, names(start)))) { # all derivs can be computed
      rjexpr <- deriv(residexpr, names(start)) ##?? could fail on some functions
    } else  rjexpr <- NULL
  if (is.null(rjexpr) && (control$derivmeth == "default")) {
    warning("Changing to alternative derivative method")
    control$derivmeth <- nlsjcontrol()$altderivmeth
  }
  
  # Define functions
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
  
ures <- resfun(prm) 
ures
wres <- swts*ures
wres
cat("derivmeth=", control$derivmeth,"\n")
ujres <- rjfun(prm)
ujres
control$derivmeth<-"numericDeriv"
cat("derivmeth=", control$derivmeth,"\n")
ujresn <- rjfun(prm)
ujresn
max(abs(ujres-ujresn))
J<-attr(ujres,"gradient")
Jn<-attr(ujresn,"gradient")
J
Jn
max(abs(J-Jn))
