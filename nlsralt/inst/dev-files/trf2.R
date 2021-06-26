trf2 <- function(par, data) {
    tt <- data[,"tt"]
    res <- par["aa"]*exp(-par["bb"]*tt) + par["cc"] - y2
}