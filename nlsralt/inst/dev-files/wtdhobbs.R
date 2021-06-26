traceval  <-  TRUE  # traceval set TRUE to debug or give full history

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)

eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
weeddata1  <-  data.frame(y=ydat, tt=tdat)

anlxb0  <-  try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata1))
print(anlxb0)

start2 <- coef(anlxb0)
wts <- 1/sqrt(ydat)
anlsw <- nls(eunsc, start=start2, trace=traceval, data=weeddata1, weights=wts)
anlsw
anlsw2 <- nls(eunsc, start=start2, trace=traceval, data=weeddata1, weights=wts^2)
anlsw2 # shows that nls() wants to minimize sum(wtsx*resid*resid)
anlxbw <- nlxb(eunsc, start=start2, trace=traceval, data=weeddata1, weights=wts)
anlxbw # but nlxb minimizes sum(wts*wts*resid*resid)
