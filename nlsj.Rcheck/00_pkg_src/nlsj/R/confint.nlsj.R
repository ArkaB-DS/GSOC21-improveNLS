confint.nlsj <- function(object, parm, level = 0.95, ...)
{
    if(!requireNamespace("MASS", quietly = TRUE))
        stop("package 'MASS' must be installed")
    confint.nlsj <- get("confint.nls", asNamespace("MASS"), inherits = FALSE)
    confint.nlsj(object, parm, level, ...)
}

