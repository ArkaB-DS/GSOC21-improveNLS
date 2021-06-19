#  File src/library/stats/R/nlsportfit.R
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
###            Nonlinear least squares for R
###
## port_cpos, port_msg() , ... are in  ==> ./nlminb.R

nls_port_fit <- function(m, start, lower, upper, control, trace, give.v=FALSE)
{
    ## Establish the working vectors and check and set options
    p <- length(par <- as.double(unlist(start)))
    iv <- integer(4L*p + 82L)
    v <- double(105L + (p * (2L * p + 20L)))
    .Call(C_port_ivset, 1, iv, v)
    if (length(control)) {
	if (!is.list(control) || is.null(nms <- names(control)))
	    stop("'control' argument must be a named list")
	## remove components that do not apply here
	for(noN in intersect(nms, c("tol", "minFactor", "warnOnly", "printEval",
                                    "scaleOffset", "nDcentral")))
	    control[[noN]] <- NULL
	nms <- names(control)
	pos <- pmatch(nms, names(port_cpos))
	if (any(nap <- is.na(pos))) {
            warning(sprintf(ngettext(length(nap),
                                     "unrecognized control element named %s ignored",
                                     "unrecognized control elements named %s ignored"),
                            paste(nms[nap], collapse = ", ")),
                    domain = NA)
	    pos <- pos[!nap]
	    control <- control[!nap]
	}
	ivpars <- pos <= 4 ; vpars <- !ivpars
	if (any(ivpars))
	    iv[port_cpos[pos[ivpars]]] <- as.integer(unlist(control[ivpars]))
	if (any(vpars))
	    v [port_cpos[pos[ vpars]]] <- as.double(unlist(control[vpars]))
    }
    if (trace)
        iv[port_cpos[["trace"]]] <- 1L
    scale <- 1
    if (any(lower != -Inf) || any(upper != Inf)) {
        low <- rep_len(as.double(lower), length(par))
        upp <- rep_len(as.double(upper), length(par))
        if(any(unlist(start) < low) ||any( unlist(start) > upp)) {
            iv[1L] <- 300
	    return(if(give.v) list(iv = iv, v = v[seq_len(18L)]) else iv)
        }
    } else
    	low <- upp <- numeric()

    if(p > 0) {
        ## driver routine port_nlsb() in ../src/port.c -- modifies m & iv
        .Call(C_port_nlsb, m,
              d = rep_len(as.double(scale), length(par)),
              df = m$gradient(), iv, v, low, upp)
    } else iv[1L] <- 6

    if(give.v)## also want v[] e.g., for attained precision
        ## v[1:18] --> ../src/portsrc.f
        list(iv = iv, v = v[seq_len(18L)]) else iv
}
