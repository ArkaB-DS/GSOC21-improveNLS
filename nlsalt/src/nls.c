/*
 *  Routines used in calculating least squares solutions in a
 *  nonlinear model in nls library for R.
 *
 *  Copyright 2005-2020 The R Core Team
 *  Copyright 2006      The R Foundation
 *  Copyright 1999-2001 Douglas M. Bates
 *                      Saikat DebRoy
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include "nls.h"

#include "internals.h"

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

/*
 * get the list element named str. names is the name attribute of list
 */

static SEXP
getListElement(SEXP list, SEXP names, const char *str)
{
    SEXP elmt = (SEXP) NULL;
    const char *tempChar;
    int i;

    for (i = 0; i < length(list); i++) {
	tempChar = CHAR(STRING_ELT(names, i)); /* ASCII only */
	if( strcmp(tempChar,str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    }
    return elmt;
}

/*
 * put some convergence-related information into list
 */
static SEXP
ConvInfoMsg(char* msg, int iter, int whystop, double fac,
	    double minFac, int maxIter, double convNew)
{
    const char *nms[] = {"isConv", "finIter", "finTol",
			 "stopCode", "stopMessage",  ""};
    SEXP ans;
    PROTECT(ans = mkNamed(VECSXP, nms));

    SET_VECTOR_ELT(ans, 0, ScalarLogical(whystop == 0)); /* isConv */
    SET_VECTOR_ELT(ans, 1, ScalarInteger(iter));	 /* finIter */
    SET_VECTOR_ELT(ans, 2, ScalarReal   (convNew));	 /* finTol */
    SET_VECTOR_ELT(ans, 3, ScalarInteger(whystop));      /* stopCode */
    SET_VECTOR_ELT(ans, 4, mkString(msg));               /* stopMessage */

    UNPROTECT(1);
    return ans;
}


/*
 *  call to nls_iter from R --- .Call("nls_iter", m, control, doTrace)
 *  where m and control are nlsModel and nlsControl objects
 *             doTrace is a logical value.
 *  m is modified; the return value is a "convergence-information" list.
 */
SEXP
nls_iter(SEXP m, SEXP control, SEXP doTraceArg)
{
    int doTrace = asLogical(doTraceArg);

    if(!isNewList(control))
	error(_("'control' must be a list"));
    if(!isNewList(m))
	error(_("'m' must be a list"));

    SEXP tmp, conv;
    PROTECT(tmp = getAttrib(control, R_NamesSymbol));

    conv = getListElement(control, tmp, "maxiter");
    if(conv == NULL || !isNumeric(conv))
	error(_("'%s' absent"), "control$maxiter");
    int maxIter = asInteger(conv);

    conv = getListElement(control, tmp, "tol");
    if(conv == NULL || !isNumeric(conv))
	error(_("'%s' absent"), "control$tol");
    double tolerance = asReal(conv);

    conv = getListElement(control, tmp, "minFactor");
    if(conv == NULL || !isNumeric(conv))
	error(_("'%s' absent"), "control$minFactor");
    double minFac = asReal(conv);

    conv = getListElement(control, tmp, "warnOnly");
    if(conv == NULL || !isLogical(conv))
	error(_("'%s' absent"), "control$warnOnly");
    int warnOnly = asLogical(conv);

    conv = getListElement(control, tmp, "printEval");
    if(conv == NULL || !isLogical(conv))
	error(_("'%s' absent"), "control$printEval");
    Rboolean printEval = asLogical(conv);

    // now get parts from 'm'  ---------------------------------
    tmp = getAttrib(m, R_NamesSymbol);

    conv = getListElement(m, tmp, "conv");
    if(conv == NULL || !isFunction(conv))
	error(_("'%s' absent"), "m$conv()");
    PROTECT(conv = lang1(conv));

    SEXP incr = getListElement(m, tmp, "incr");
    if(incr == NULL || !isFunction(incr))
	error(_("'%s' absent"), "m$incr()");
    PROTECT(incr = lang1(incr));

    SEXP deviance = getListElement(m, tmp, "deviance");
    if(deviance == NULL || !isFunction(deviance))
	error(_("'%s' absent"), "m$deviance()");
    PROTECT(deviance = lang1(deviance));

    SEXP trace = getListElement(m, tmp, "trace");
    if(trace == NULL || !isFunction(trace))
	error(_("'%s' absent"), "m$trace()");
    PROTECT(trace = lang1(trace));

    SEXP setPars = getListElement(m, tmp, "setPars");
    if(setPars == NULL || !isFunction(setPars))
	error(_("'%s' absent"), "m$setPars()");
    PROTECT(setPars);

    SEXP getPars = getListElement(m, tmp, "getPars");
    if(getPars == NULL || !isFunction(getPars))
	error(_("'%s' absent"), "m$getPars()");
    PROTECT(getPars = lang1(getPars));

    SEXP pars = PROTECT(eval(getPars, R_GlobalEnv));
    int nPars = LENGTH(pars);

    double dev = asReal(eval(deviance, R_GlobalEnv));
    if(doTrace) eval(trace,R_GlobalEnv);

    double fac = 1.0;
    Rboolean hasConverged = FALSE;
    SEXP newPars = PROTECT(allocVector(REALSXP, nPars));
    int evaltotCnt = 1;
    double convNew = -1. /* -Wall */;
    int i;
#define CONV_INFO_MSG(_STR_, _I_)				\
	ConvInfoMsg(_STR_, i, _I_, fac, minFac, maxIter, convNew)

#define NON_CONV_FINIS(_ID_, _MSG_)		\
    if(warnOnly) {				\
	warning(_MSG_);				\
	return CONV_INFO_MSG(_MSG_, _ID_);      \
    }						\
    else					\
	error(_MSG_);

#define NON_CONV_FINIS_1(_ID_, _MSG_, _A1_)	\
    if(warnOnly) {				\
	char msgbuf[1000];			\
	warning(_MSG_, _A1_);			\
	snprintf(msgbuf, 1000, _MSG_, _A1_);	\
	return CONV_INFO_MSG(msgbuf, _ID_);	\
    }						\
    else					\
	error(_MSG_, _A1_);

#define NON_CONV_FINIS_2(_ID_, _MSG_, _A1_, _A2_)	\
    if(warnOnly) {					\
	char msgbuf[1000];				\
	warning(_MSG_, _A1_, _A2_);			\
	snprintf(msgbuf, 1000, _MSG_, _A1_, _A2_);	\
	return CONV_INFO_MSG(msgbuf, _ID_);		\
    }							\
    else						\
	error(_MSG_, _A1_, _A2_);

    for (i = 0; i < maxIter; i++) { // ---------------------------------------------

	if((convNew = asReal(eval(conv, R_GlobalEnv))) <= tolerance) {
	    hasConverged = TRUE;
	    break;
	}

	SEXP newIncr = PROTECT(eval(incr, R_GlobalEnv));
	double
	    *par   = REAL(pars),
	    *npar  = REAL(newPars),
	    *nIncr = REAL(newIncr);
	int evalCnt = -1;
	if(printEval)
	    evalCnt = 1;

	while(fac >= minFac) { // 1-dim "line search"
	    if(printEval) {
		Rprintf("  It. %3d, fac= %11.6g, eval (no.,total): (%2d,%3d):",
			i+1, fac, evalCnt, evaltotCnt);
		evalCnt++;
		evaltotCnt++;
	    }
	    for(int j = 0; j < nPars; j++)
		npar[j] = par[j] + fac * nIncr[j];

	    PROTECT(tmp = lang2(setPars, newPars));
	    if (asLogical(eval(tmp, R_GlobalEnv))) { /* singular gradient */
		UNPROTECT(11);

		NON_CONV_FINIS(1, _("singular gradient"));
	    }
	    UNPROTECT(1);

	    double newDev = asReal(eval(deviance, R_GlobalEnv));
	    if(printEval)
		Rprintf(" new dev = %g\n", newDev);
	    if(newDev <= dev) {
		dev = newDev;
		tmp = newPars;
		newPars = pars;
		pars = tmp;
		fac = MIN(2*fac, 1);
		break;
	    } // else
	    fac /= 2.;
	}
	UNPROTECT(1);
	if(doTrace) eval(trace, R_GlobalEnv);
	if( fac < minFac ) {
	    UNPROTECT(9);
	    NON_CONV_FINIS_2(2,
			     _("step factor %g reduced below 'minFactor' of %g"),
			     fac, minFac);
	}
    }

    UNPROTECT(9);
    if(!hasConverged) {
	NON_CONV_FINIS_1(3,
			 _("number of iterations exceeded maximum of %d"),
			 maxIter);
    }
    else
	return CONV_INFO_MSG(_("converged"), 0);
}
#undef CONV_INFO_MSG
#undef NON_CONV_FINIS
#undef NON_CONV_FINIS_1
#undef NON_CONV_FINIS_2
