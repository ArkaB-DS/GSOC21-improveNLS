## ReplaceNLS.R -- code to source() to replace functions in base-R
numericDeriv <- function(expr, theta, rho = parent.frame(), dir = 1,
                         eps = .Machine$double.eps ^ (1/if(central) 3 else 2), central = FALSE)
{
    cat("Replacement numericDeriv\n")  # At the moment simply prints this
    dir <- rep_len(dir, length(theta))
    stopifnot(is.finite(eps), eps > 0)
    val <- .Call(C_numeric_deriv, expr, theta, rho, dir, eps, central) ## ../src/nls.c
    # here
    if (!is.null(d <- dim(val))) {
        if(d[length(d)] == 1L)
            d <- d[-length(d)]
        if(length(d) > 1L)
            dim(attr(val, "gradient")) <- c(d, dim(attr(val, "gradient"))[-1L])
    }
    val
}