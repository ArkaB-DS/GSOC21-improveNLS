## weighted nlsj fit
rm(list=ls())
set.seed(123)
y <- x <- 1:10
yeps <- y + rnorm(length(y), sd = 0.01)
wts <- rep(c(1, 2), length = 10); wts[5] <- 0
fit0 <- lm(yeps ~ x, weights = wts)
## DOES NOT WORK: fit0x <- lm(yeps ~ a+b*x, weights=wts)
## IGNORE_RDIFF_BEGIN
summary(fit0, cor = TRUE)
cf0 <- coef(summary(fit0))[, 1:2]
fit <- nlsj(yeps ~ a + b*x, start = list(a = 0.12345, b = 0.54321),  weights = wts, trace = TRUE)
summary(fit, cor = TRUE)
df <- data.frame(x=x, yeps=yeps)
## IGNORE WHEN RUNNING R CMD check since needs external pkg
## fit2 <- nlsr::nlxb(yeps ~ a + b*x, start = list(a = 0.12345, b = 0.54321), data=df,  weights = wts, trace = TRUE)
fitn <- nls(yeps ~ a + b*x, start = list(a = 0.12345, b = 0.54321), weights = wts, trace = TRUE)
summary(fitn)
## IGNORE_RDIFF_END
stopifnot(all.equal(as.numeric(residuals(fit)), as.numeric(residuals(fit0)), tolerance = 1e-5,
                    check.attributes = FALSE))
stopifnot(all.equal(residuals(fit), residuals(fit0), tolerance = 1e-5,
                    check.attributes = FALSE))
all.equal(residuals(fit), residuals(fit0), tolerance = 1e-5, check.attributes = FALSE)
(as.numeric(fitn$m$resid())/as.numeric(resid(fitn)))^2-wts
fnlsj <- fit$m$fitted()
fnlsj
stopifnot(df.residual(fit) == df.residual(fit0))
stopifnot(all.equal(logLik(fit), logLik(fit0), tolerance = 1e-8))
cf1 <- coef(summary(fit))[, 1:2]
