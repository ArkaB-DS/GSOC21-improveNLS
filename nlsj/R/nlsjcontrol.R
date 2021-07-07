#  File src/library/stats/R/nlscontrol.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 2000-2020 The R Core Team
#  Copyright (C) 1999-1999 Saikat DebRoy, Douglas M. Bates, Jose C. Pinheiro
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
###     nlscontrol for Nonlinear least squares for R
###


nlsj.control <- function(maxiter = 500, tol = 0.00001, minFactor = 1/1024,
			printEval = FALSE, warnOnly = FALSE, scaleOffset = 0,
                        nDcentral = FALSE, watch = FALSE, phi = 1, lamda = 0, 
			offset = 100, laminc = 10, lamdec = 0.4, resmax = 10000, 
			rofftest = TRUE, smallsstest = TRUE,
			derivmeth="numericDeriv", altderivmeth="numericDeriv",
			trace=FALSE){
# Note that trace has been added to nlsj.control(). It may appear also in the 
# main call of nlsj(), so we need to ensure these do NOT conflict.??
    stopifnot(is.numeric(tol), length(tol) == 1L, tol > 0,
              is.numeric(minFactor),   length(minFactor) == 1L,
              is.numeric(scaleOffset), length(scaleOffset) == 1L,
              is.logical(nDcentral), length(nDcentral) == 1L, !is.na(nDcentral),
              (derivmeth %in% c("default", "numericDeriv", "numDer0", "numDerC")),
              altderivmeth=="numericDeriv")
# derivmeth:
#  default -- try to use analytic derivs, and if that fails, then numericDeriv
#  numericDeriv -- R version of original base R numericDeriv from nls.R
#  numDer0 -- default jacobian from numDeriv package -- ?? NOT yet active
#  numDerC -- complex derivative option of numDeriv package -- ?? NOT yet active

    list(maxiter = maxiter, tol = tol, minFactor = minFactor,
	 printEval = printEval, warnOnly = warnOnly,
         scaleOffset = scaleOffset, nDcentral = nDcentral,
	 watch = watch, phi = phi, lamda = lamda, offset = offset, 
         laminc = laminc, lamdec = lamdec, resmax = resmax, 
         rofftest = rofftest, smallsstest = smallsstest, derivmeth=derivmeth,
         altderivmeth=altderivmeth, trace=trace)
}
