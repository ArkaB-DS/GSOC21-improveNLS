#  File src/library/stats/R/nlsnd.R
#  Part of the modified R package, https://www.R-project.org
#
#  Copyright (C) 2000-2020 The R Core Team
#  Copyright (C) 1999-1999 Saikat DebRoy, Douglas M. Bates, Jose C. Pinheiro
#  J C Nash 2021
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
###            numeric Jacobian for Nonlinear least squares for R
###

numericDeriv <- function(expr, theta, rho = parent.frame(), dir = 1,
                 eps = .Machine$double.eps ^ (1/if(central) 3 else 2), central = FALSE)
## Note: this expr must be set up as a call to work properly according to JN??
## ?? we set eps conditional on central. But central set AFTER eps. Is this OK.
{   ndtrace<-FALSE
#    if(ndtrace) cat("numericDeriv-Alt\n")
    dir <- rep_len(dir, length(theta))
    stopifnot(is.finite(eps), eps > 0)
    rho1 <- new.env(FALSE, rho, 0) # 0 will be ignored since hash is FALSE
    if (!is.character(theta) ) {stop("'theta' should be of type character")}
    if (is.null(rho)) {
            stop("use of NULL environment is defunct")
    } else {
          if(! is.environment(rho)) {stop("'rho' should be an environment")}
    }
    if( ! ((length(dir) == length(theta) ) & (is.numeric(dir) ) ) )
              {stop("'dir' is not a numeric vector of the correct length") }
    if(is.na(central)) { stop("'central' is NA, but must be TRUE or FALSE") }
    ##?? should check if we have the residuals already!
    res0 <- eval(expr, rho) # the base residuals. ?? C has a check for REAL ANS=res0
    if (any(is.infinite(res0)) ) {stop("residuals cannot be evaluated at base point")}
    ##  CHECK_FN_VAL(res, ans);  ?? how to do this. Is it necessary?
    nt <- length(theta) # number of parameters
    mr <- length(res0) # number of residuals
    JJ <- matrix(NA, nrow=mr, ncol=nt) # Initialize the Jacobian
    for (j in 1:nt){
       origPar<-get(theta[j],rho)
       xx <- abs(origPar)
       delta <- if (xx == 0.0) {eps} else { xx*eps }
       ## JN: I prefer eps*(xx + eps)  which is simpler ?? Should we suggest / use a control switch
       prmx<-origPar+delta*dir[j]
       assign(theta[j],prmx,rho)
       res1 <- eval(expr, rho) # new residuals (forward step)
       if (central) { # compute backward step resids for central diff
          prmb <- origPar - dir[j]*delta
          assign(theta[j], prmb, envir=rho) # may be able to make more efficient later??
          resb <- eval(expr, rho)
          JJ[, j] <- dir[j]*(res1-resb)/(2*delta) # vectorized
       } else { ## forward diff
          JJ[,j] <- dir[j]*(res1-res0)/delta
       }  # end forward diff
       assign(theta[j],origPar,rho) # restore the parameter value !! IMPORTANT
    } # end loop over the parameters
    attr(res0, "gradient") <- JJ
    if (ndtrace){
       cat("par:")
       for (j in 1:nt){ cat(get(theta[j],rho)," ") }
       cat("\n")
       print(res0)
    }
    return(res0)
}


