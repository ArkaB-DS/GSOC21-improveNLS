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


nlsj.control <- function(maxiter = 50, tol = 0.00001, minFactor = 1/1024,
			printEval = FALSE, warnOnly = FALSE, scaleOffset = 0,
                        nDcentral = FALSE) {
    stopifnot(is.numeric(tol), length(tol) == 1L, tol > 0,
              is.numeric(minFactor),   length(minFactor) == 1L,
              is.numeric(scaleOffset), length(scaleOffset) == 1L,
              is.logical(nDcentral), length(nDcentral) == 1L, !is.na(nDcentral))
    list(maxiter = maxiter, tol = tol, minFactor = minFactor,
	 printEval = printEval, warnOnly = warnOnly,
         scaleOffset = scaleOffset, nDcentral = nDcentral)
}
