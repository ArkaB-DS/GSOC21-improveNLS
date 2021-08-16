# NLSProbName: hobbs-base
# NLSProbDescription: {The Hobbs weed infestation problem to estimate a 
#    3-parameter logistic S curve in its unscaled form from a reasonably
#    easy starting point of (200, 50, 0.3)
#  }
## Arkajyoti: Do we want to use names like these so we can find things
## with programs? There may be useful ideas like ROxygen, but I haven't
## used them. 

SS3logis <- # selfStart(~ Asym/(1 + exp((xmid - input)/scal)) 
	     ##  ==> b1/(1+b2*exp(-b3*input))
    selfStart(
        function(input, b1, b2, b3)
        {
               Asym <- b1
              # b2 = exp(xmid/scal)  --> b2 = exp(xmid*b3) -> log(b2) = xmid*b3
               xmid <- log(b2)/b3
              # b3 = 1/scal  ---> 
               scal <- 1/b3
              .expr1 <- xmid - input
              .expr3 <- exp(.e2 <- .expr1/scal)
              .expr4 <- 1 + .expr3
              .value <- Asym/.expr4
              .actualArgs <- as.list(match.call()[c("Asym", "xmid", "scal")])
              if(all(vapply(.actualArgs, is.name, NA)))
              {
		  .expr10 <- .expr4^2
                  .grad <- array(0, c(length(.value), 3L), list(NULL, c("Asym", "xmid", "scal")))
                  .grad[, "Asym"] <- 1/.expr4
		  .grad[, "xmid"] <- - (xm <- Asym * .expr3/scal/.expr10)
		  .grad[, "scal"] <- xm * .e2
                  dimnames(.grad) <- list(NULL, .actualArgs)
                  attr(.value, "gradient") <- .grad
              }
              .value
        },
        initial = function(mCall, data, LHS, ...) {
              xy <- sortedXyData(mCall[["input"]], LHS, data)
              if(nrow(xy) < 4) {
                  stop("too few distinct input values to fit a logistic model")
              }
              z <- xy[["y"]]
              ## transform to proportion, i.e. in (0,1) :
              rng <- range(z); dz <- diff(rng)
              z <- (z - rng[1L] + 0.05 * dz)/(1.1 * dz)
              xy[["z"]] <- log(z/(1 - z))		# logit transformation
              aux <- coef(lm(x ~ z, xy))
              ## ?? THIS IS CIRCULAR!! Depends on nls() plinear
              pars <- coef(nls(y ~ 1/(1 + exp((xmid - x)/scal)),
                               data = xy,
                               start = list(xmid = aux[[1L]], scal = aux[[2L]]),
                               algorithm = "plinear", ...))
              setNames(pars [c(".lin", "xmid", "scal")],
                       mCall[c("b1", "b2", "b3")])
        },
        parameters = c("b1", "b2", "b3"))


# Use the Hobbs Weed data
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12

NLSformula0 <- y ~ b1/(1+b2*exp(-b3*tt))
NLSformula <- y ~ SS3logis(tt, b1, b2, b3)
NLStestdata <- data.frame(y=weed, tt=tt) # should we use standard name?
## NLSstart <- c(b1=200, b2=50, b3=0.3) # an easy starting vector (named!)
NLSstart <- getInitial(NLSformula, NLStestdata)
cat("SS3logis selfStart parameters:")
print(as.numeric(NLSstart))
## Unscaled model
NLSlower <- NULL
NLSupper <- NULL
# We want to get the initial vector then run with ORIGINAL formula as a test
NLSrunline <- "(formula=NLSformula0, data=NLStestdata, start=NLSstart)"
