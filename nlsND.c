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
 *  call to numeric_deriv from R -
 *  .Call("numeric_deriv", expr, theta, rho, dir = 1., eps = .Machine$double.eps, central=FALSE)
 *  Returns: ans
 */
SEXP
numeric_deriv(SEXP expr, SEXP theta, SEXP rho, SEXP dir, SEXP eps_, SEXP centr,
              SEXP rho1)
{
    if(!isString(theta))
	error(_("'theta' should be of type character"));
    if (isNull(rho)) {
	error(_("use of NULL environment is defunct"));
	rho = R_BaseEnv;
    } else
	if(!isEnvironment(rho))
	    error(_("'rho' should be an environment"));
    int nprot = 3;
    if(TYPEOF(dir) != REALSXP) {
	PROTECT(dir = coerceVector(dir, REALSXP)); nprot++;
    }
    if(LENGTH(dir) != LENGTH(theta))
	error(_("'dir' is not a numeric vector of the correct length"));
    Rboolean central = asLogical(centr);
    if(central == NA_LOGICAL)
	error(_("'central' is NA, but must be TRUE or FALSE"));
//    SEXP rho1 = PROTECT(R_NewEnv(rho, FALSE, 0));
//    nprot++;
    SEXP
	pars = PROTECT(allocVector(VECSXP, LENGTH(theta))),
        ans  = PROTECT(duplicate(eval(expr, rho1)));
    double *rDir = REAL(dir),  *res = NULL; // -Wall
#define CHECK_FN_VAL(_r_, _ANS_) do {					\
    if(!isReal(_ANS_)) {						\
	SEXP temp = coerceVector(_ANS_, REALSXP);			\
	UNPROTECT(1);/*: _ANS_ *must* have been the last PROTECT() ! */ \
	PROTECT(_ANS_ = temp);						\
    }									\
    _r_ = REAL(_ANS_);							\
    for(int i = 0; i < LENGTH(_ANS_); i++) {				\
	if (!R_FINITE(_r_[i]))						\
	    error(_("Missing value or an infinity produced when evaluating the model")); \
    }									\
} while(0)

    CHECK_FN_VAL(res, ans);

    const void *vmax = vmaxget();
    int lengthTheta = 0;
    for(int i = 0; i < LENGTH(theta); i++) {
	const char *name = translateChar(STRING_ELT(theta, i));
	SEXP s_name = install(name);
	SEXP temp = findVar(s_name, rho1);
	if(isInteger(temp))
	    error(_("variable '%s' is integer, not numeric"), name);
	if(!isReal(temp))
	    error(_("variable '%s' is not numeric"), name);
	// We'll be modifying the variable, so need to make a copy PR#15849
	defineVar(s_name, temp = duplicate(temp), rho1);
	MARK_NOT_MUTABLE(temp);
	SET_VECTOR_ELT(pars, i, temp);
	lengthTheta += LENGTH(VECTOR_ELT(pars, i));
    }
    vmaxset(vmax);
    SEXP gradient = PROTECT(allocMatrix(REALSXP, LENGTH(ans), lengthTheta));
    double *grad = REAL(gradient);
    double eps = asReal(eps_); // was hardcoded sqrt(DOUBLE_EPS) { ~= 1.49e-08, typically}
    for(int start = 0, i = 0; i < LENGTH(theta); i++) {
	double *pars_i = REAL(VECTOR_ELT(pars, i));
	for(int j = 0; j < LENGTH(VECTOR_ELT(pars, i)); j++, start += LENGTH(ans)) {
	    double
		origPar = pars_i[j],
		xx = fabs(origPar),
		delta = (xx == 0) ? eps : xx*eps;
	    pars_i[j] += rDir[i] * delta;
	    SEXP ans_del = PROTECT(eval(expr, rho1));
	    double *rDel = NULL;
	    CHECK_FN_VAL(rDel, ans_del);
	    if(central) {
		pars_i[j] = origPar - rDir[i] * delta;
		SEXP ans_de2 = PROTECT(eval(expr, rho1));
		double *rD2 = NULL;
		CHECK_FN_VAL(rD2, ans_de2);
		for(int k = 0; k < LENGTH(ans); k++) {
		    grad[start + k] = rDir[i] * (rDel[k] - rD2[k])/(2 * delta);
		}
	    } else { // forward difference  (previously hardwired):
		for(int k = 0; k < LENGTH(ans); k++) {
		    grad[start + k] = rDir[i] * (rDel[k] - res[k])/delta;
		}
	    }
	    UNPROTECT(central ? 2 : 1); // ansDel & possibly ans
	    pars_i[j] = origPar;
	}
    }
    setAttrib(ans, install("gradient"), gradient);
    UNPROTECT(nprot);
    return ans;
}
