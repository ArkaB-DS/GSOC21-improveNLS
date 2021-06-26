rjfundoc <- function(fun, savefile=NULL) {
  efun <- environment(fun) # get the environment with data, expression, etc
  avn <- all.vars(efun$modelformula) # vars and parameters
  pnames <- names(efun$pvec) # assumes that vector is named (normal)
# ?? what do we do when it is not -- need to name it p1, p2, etc.
# DJM:  the model2rjfun* functions do this; if a user creates fun some other way,
# they'd better do it too!
  iprm <- match(pnames, avn)
  if (length(iprm))
    notprm <- avn[-iprm]
  else
    notprm <- avn
  funname <- deparse(match.call()[[2]])
  data <- efun$data
  # data is an environment, but not all variables are in the top level, so use get() to find them
  modeldata <- mget(notprm, envir=data, inherits=TRUE, ifnotfound=list(function(n) NULL))
  notdata <- setdiff(notprm, names(modeldata))
  resids <- fun(efun$pvec)
  n <- length(resids)
  islengthn <- sapply(modeldata, function(col) length(col) == n)
  result <- structure(list(funname=funname, modelformula=efun$modelformula,
                           modelexpr=modelexpr(fun), n=n,
                           pvec=efun$pvec, data=as.data.frame(modeldata[islengthn]), 
                           extradata=modeldata[!islengthn], unknown = notdata),
                      class = "rjfundoc")
  if (!is.null(savefile)) {
    sink(savefile)
    print(result)
    sink()
  }
  result
}
  
print.rjfundoc <- function(x, ...) {
  cat("FUNCTION", x$funname, "\n")
  cat("Formula:\t")
  print(x$modelformula)
  cat("Code:\t\t")
  print(x$modelexpr)
  cat("Parameters:\t", paste(names(x$pvec), collapse=", "), "\n")
  cat("Data:\t\t", paste(names(x$data), collapse=", "), "\n")
  if (length(x$extradata))
    cat("Extra:\t\t", paste(names(x$extradata), collapse=", "), "\n")
  if (length(x$unknown)) 
    cat("Unknown symbols:\t", paste(x$unknown, collapse=", "), "\n")
  cat("\nVALUES\n")
  cat("Observations:\t", x$n, "\n")
  cat("Parameters:\n")
  print(x$pvec)
  if (length(x$data)) {
    cat("Data (length ", x$n, "):\n", sep="")
    print(x$data)
  }
  if (length(x$extradata)) {
    cat("Extra:\n")
    print(x$extradata)
  }
  x
}
