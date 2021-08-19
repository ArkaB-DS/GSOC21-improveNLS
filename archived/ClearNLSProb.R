# ClearNLSProb.R
# J C Nash 2021-7-22
# Clear objects created by NLS test problems
rmNLS <- function(envir=.GlobalEnv) { # Do we need to specify .GlobalEnv or other environment?
  g <- .GlobalEnv
  removed <- NULL
  torm <- c("NLSformula", "NLStestdata", "NLSstart", "NLSupper", "NLSlower", "NLSweights", "NLSsubset")
  nrm <- length(torm)
  print(ls(g))
# cat("Is g locked:", environmentIsLocked(g),"\n")
  cat("Search list element 1 is ", search()[1],"\n")
  for (i in 1:nrm){
     ob <- torm[i]    
     if (exists(ob, envir=g, inherits=TRUE)) {
       result<-rm(list=ob, envir=g) # NOTE: "list=" is critical
       removed <- c(removed, ob)
       cat("removed ",ob," from .GlobalEnv\n")
     }
  }
  removed
}
