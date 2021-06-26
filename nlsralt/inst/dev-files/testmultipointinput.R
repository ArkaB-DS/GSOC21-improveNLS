library(nlsr)
ff <- ~ (x^2 + 4 * y^2)
myder <- fnDeriv(ff, c("x","y"), verbose=TRUE, hessian=TRUE)
myder
myder(1,2)
x <- 1:3
y <- c(2, 4, 6)
myder(x, y)

