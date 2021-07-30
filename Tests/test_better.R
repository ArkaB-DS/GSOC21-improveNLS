test_better <- function( objectnew, objectold){
## 2021-7-30 Attempt to describe philosophy when new method does not fail
##  example: singular gradient
# We will want to use skip() in some of the cases
   pass <- FALSE
   if (expect_error(objectold)){ # this is as we expect
          cat("Old method produces error\n")
   } else { stop("Old method somehow managed to avoid error") 
          # return(FALSE)}
   if (expect_error(objectnew) { # this is not as we want
          stop("New method still failing")
   } else { cat("New method works when old one does not\n") }
   pass <- TRUE
   pass
}