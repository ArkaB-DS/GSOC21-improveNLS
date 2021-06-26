require(Deriv)
#- First must provide some data
x <- 3
y <- 4

f <- function(x) x^2
Deriv(f)
#- Should see
#- function (x)
#- 2 * x
#- Now save the derivative
f1 <- Deriv(f)
f1 #- print it
f2 <- Deriv(f1) #- and take second derivative
f2 #- print it

f <- function(x, y) sin(x) * cos(y)
f_ <- Deriv(f)
f_ #- print it
#- Should see
#- function (x, y)
#- c(x = cos(x) * cos(y), y = -(sin(x) * sin(y)))
f_(3, 4)
#- Should see
#-              x         y
#- [1,] 0.6471023 0.1068000

f2 <- Deriv(~ f(x, y^2), "y") #- This has a tilde to render the 1st argument as a formula object
#- Also we are substituting in y^2 for y
f2 #- print it
#- -(2 * (y * sin(x) * sin(y^2)))
mode(f2) #- check what type of object it is
arg1 <- ~ f(x,y^2)
mode(arg1) #- check the type
f2a <- Deriv(arg1, "y")
f2a #- and print to see if same as before
#- try evaluation of f using current x and y
x
y
f(x,y^2)
eval(f2a) #- We need x and y defined to do this.

f3 <- Deriv(quote(f(x, y^2)), c("x", "y"), cache.exp=FALSE) #- check cache.exp operation
#- Note that we need to quote or will get evaluation at current x, y values (if they exist)
f3 #- print it
#- c(x = cos(x) * cos(y^2), y = -(2 * (y * sin(x) * sin(y^2))))
f3c <- Deriv(quote(f(x, y^2)), c("x", "y"), cache.exp=TRUE) #- check cache.exp operation
f3c #- print it
#- Now want to evaluate the results
eval(f3c)
#- Should see
#- x         y 
#- 0.9480757 0.3250313 
eval(f3) #- check this also
#- or we can create functions
f3cf <- function(x, y){eval(f3c)}
f3cf(x=1, y=2)
#-          x          y 
#- -0.3531652  2.5473094 
f3f <- function(x,y){eval(f3)}
f3f(x=3, y=4)
#-         x         y 
#- 0.9480757 0.3250313 

#- try an expression
Deriv(expression(sin(x^2) * y), "x")
#- should see
#- expression(2 * (x * y * cos(x^2)))

#- quoted string
Deriv("sin(x^2) * y", "x") # differentiate only by x
#- Should see
#- "2 * (x * y * cos(x^2))"

Deriv("sin(x^2) * y", cache.exp=FALSE) #- differentiate by all variables (here by x and y)
#- Note that default is to differentiate by all variables.
#- Should see
#- "c(x = 2 * (x * y * cos(x^2)), y = sin(x^2))"

#- Compound function example (here abs(x) smoothed near 0)
#- Note that this introduces the possibilty of `if` statements in the code
#- BUT (JN) seems to give back quoted string, so we must parse.
fc <- function(x, h=0.1) if (abs(x) < h) 0.5*h*(x/h)**2 else abs(x)-0.5*h
efc1 <- Deriv("fc(x)", "x", cache.exp=FALSE) 
#- "if (abs(x) < h) x/h else sign(x)"
#- A few checks on the results
efc1

fc1 <- function(x,h=0.1){ eval(parse(text=efc1)) }
fc1
## h=0.1
fc1(1)
fc1(0.001)
fc1(-0.001)
fc1(-10)
fc1(0.001, 1)

#- Example of a first argument that cannot be evaluated in the current environment:
try(suppressWarnings(rm("xx", "yy"))) #- Make sure there are no objects xx or yy
Deriv(~xx^2+yy^2)
#- Should show
#- c(xx = 2 * xx, yy = 2 * yy)
#- ?? What is the meaning / purpose of this construct?

#- ?? Is following really AD?  
#- Automatic differentiation (AD), note intermediate variable 'd' assignment
Deriv(~{d <- ((x-m)/s)^2; exp(-0.5*d)}, "x")
# Note that the result we see does NOT match what follows in the example(Deriv) (JN ??)
#{
#   d <- ((x - m)/s)^2
#   .d_x <- 2 * ((x - m)/s^2)
#   -(0.5 * (.d_x * exp(-(0.5 * d))))
#}
#- For some reason the intermediate variable d is NOT included.??


#- Custom derivative rule. Note that this needs explanations??
myfun <- function(x, y=TRUE) NULL #- do something useful
dmyfun <- function(x, y=TRUE) NULL #- myfun derivative by x.
drule[["myfun"]] <- alist(x=dmyfun(x, y), y=NULL) #- y is just a logical
Deriv(myfun(z^2, FALSE), "z")
# 2 * (z * dmyfun(z^2, FALSE))

#- Differentiation by list components
theta <- list(m=0.1, sd=2.) #- Why do we set values??
x <- names(theta) #- and why these particular names??
names(x)=rep("theta", length(theta))
Deriv(~exp(-(x-theta$m)**2/(2*theta$sd)), x, cache.exp=FALSE)
#- Should show the following (but why??)
#- c(theta_m = exp(-((x - theta$m)^2/(2 * theta$sd))) *
#-  (x - theta$m)/theta$sd, theta_sd = 2 * (exp(-((x - theta$m)^2/
#-  (2 * theta$sd))) * (x - theta$m)^2/(2 * theta$sd)^2))
lderiv <-  Deriv(~exp(-(x-theta$m)**2/(2*theta$sd)), x, cache.exp=FALSE)
fld <- function(x){ eval(lderiv)} #- put this in a function
fld(2) #- and evaluate at a value
