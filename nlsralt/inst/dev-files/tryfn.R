tryfn <- function (Input, hi.y, HillSlope, logEC50) 
{
  .expr1 <- -HillSlope
  .expr2 <- log(Input)
  .expr3 <- .expr2 - logEC50
  .expr4 <- .expr1 * .expr3
  .expr5 <- exp(.expr4)
  .expr6 <- 1 + .expr5
  .expr7 <- .expr6^2
  .expr8 <- -.expr3
  .expr9 <- .expr5 * .expr8
  .expr10 <- .expr5 * HillSlope
  .expr11 <- hi.y * .expr9
  .expr12 <- 2 * .expr6
  .expr13 <- .expr7^2
  .expr14 <- hi.y * .expr10
  .expr15 <- .expr12 * .expr10
  .value <- hi.y/(1 + exp(.expr1 * (.expr3)))
  .grad <- array(0, c(length(.value), 3L), list(NULL, c("hi.y", 
                                                        "HillSlope", "logEC50")))
  .grad[, "hi.y"] <- 1/.expr6
  .grad[, "HillSlope"] <- -(.expr11/.expr7)
  .grad[, "logEC50"] <- -(.expr14/.expr7)
  .hessian <- array(0, c(length(.value), 3L, 3L), list(NULL,
     c("hi.y", "HillSlope", "logEC50"), c("hi.y", "HillSlope","logEC50")))
  .hessian[, "hi.y", "hi.y"] <- 0
  .hessian[, "hi.y", "HillSlope"] <- .hessian[, "HillSlope", 
                                              "hi.y"] <- -(.expr9/.expr7)
  .hessian[, "hi.y", "logEC50"] <- .hessian[, "logEC50", "hi.y"] <- -(.expr10/.expr7)
  .hessian[, "HillSlope", "HillSlope"] <- -(hi.y * (.expr9 * 
                                                      .expr8)/.expr7 - .expr11 * (.expr12 * .expr9)/.expr13)
  .hessian[, "HillSlope", "logEC50"] <- .hessian[, "logEC50", 
                                                 "HillSlope"] <- -(hi.y * (.expr5 + .expr10 * .expr8)/.expr7 - 
                                                                     .expr11 * .expr15/.expr13)
  .hessian[, "logEC50", "logEC50"] <- -(hi.y * (.expr10 * HillSlope)/.expr7 - 
                                          .expr14 * .expr15/.expr13)
  attr(.value, "hessian") <- .hessian
  attr(.value, "gradient") <- .grad
  .value
}
