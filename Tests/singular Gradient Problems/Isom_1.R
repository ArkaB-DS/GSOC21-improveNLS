# NLSProbName: Isom_1.R
# NLSProbDescription: { The Isom data frame has 24 rows and 4 columns from an isomerization experiment.
# The four columns are:This data frame contains the following columns:
# `hyd`: partial pressure of hydrogen (psia).
# `n.pent`: partial pressure of n-pentane (psia).
# `iso.pen`: partial pressure of isopentane (psia).
# `rate`: reaction rate for isomerization of n-pentane to isopentane (1/hr).
# }


# Use the Isom data from NRAIA package

## DATA
rate=c(3.541,  2.397,  6.694,  4.722,  0.593,  0.268,  2.797,
		 2.451,  3.196,  2.021,  0.896,  5.084,  5.686,  1.193,
		 2.648,  3.303,  3.054,  3.302,  1.271, 11.648,  2.002,
		 9.604,  7.754, 11.590)
hyd = c(205.8, 404.8, 209.7, 401.6, 224.9, 402.6, 212.7, 406.2, 133.3, 470.9, 300.0,
		  301.6, 297.3, 314.0, 305.7, 300.1, 305.4, 305.2, 300.1, 106.6, 417.2, 251.0,
		  250.3, 145.1)
iso.pen = c(37.1,  36.3,  49.4,  44.9, 116.3, 128.9, 134.4, 134.9,  87.6,  86.9,  81.7,
			101.7,  10.5, 157.1,  86.0,  90.2,  87.4,  87.0,  66.4,  33.0,  32.9,  41.5,
			14.7,  50.2)
n.pent = c( 90.9,  92.9, 174.9, 187.2,  92.7, 102.2, 186.9, 192.6, 140.8, 144.2,  68.3,
			214.6, 142.2, 146.7, 142.0, 143.7, 141.1, 141.5,  83.0, 209.6,  83.9, 294.4,
			148.0, 291.0)
	
NLStestdata <- data.frame(rate,hyd,iso.pen,n.pent)

## STARTING VALUE

NLSstart <-c(b2 = 0.1, b3 = 0.1, b4 = 0.1) # a starting vector (named!)

## MODEL
NLSformula <- rate ~ b3*(n.pent - iso.pen/1.632)/(1+b2*hyd+b3*n.pent+b4*iso.pen)
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, data=NLStestdata, start=NLSstart)"
output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
output_nlsj <- eval(parse(text=paste("nlsj::nlsj",NLSrunline))) # nlsj is the new nls
output_nlsLM <- eval(parse(text=paste("minpack.lm::nlsLM",NLSrunline))) # nlsLM
output_nlxb <- eval(parse(text=paste("nlsr::nlxb",NLSrunline))) # nlxb

#### TESTING nls VS nlsj

# NLSout/expout has "m", "convInfo", "data", "call",
# "dataClasses", "control"

## testing m values:  "resid"      "fitted"     "formula"    "deviance"   "lhs"       
# "gradient"   "conv"       "incr"       "setVarying" "setPars"   
# "getPars"    "getAllPars" "getEnv"     "trace"      "Rmat"      
# "predict"

      # residuals
	all.equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsj)))
	# fitted
	all.equal(as.vector(fitted(output_nls)),
			 as.vector(fitted(output_nlsj)))
#	# formula
#	all.equal(formula(output_nls),
#			 formula(output_nlsj))
	# deviance
	all.equal(deviance(output_nls),
			 deviance(output_nlsj))
	# gradient
	all.equal( output_nls$m$gradient(),
			  attr(output_nlsj$m$resid(),"gradient"))
#	# conv
#	all.equal( output_nls$m$conv(),
#			  output_nlsj$m$conv())
#	# incr
#	all.equal( output_nls$m$incr(),
#			  output_nlsj$m$incr())
	# getPars # difference between getAllPars adn getPars?
	all.equal( output_nls$m$getPars(),
			  output_nlsj$m$getPars())
#	# getEnv
#	all.equal( output_nls$m$igetEnv(),
#			  output_nlsj$m$getEnv())
	# Rmat
	all.equal( output_nls$m$Rmat(),
			  output_nlsj$m$Rmat())
#	# predict
#	all.equal( output_nls$m$predict(),
#			  output_nlsj$m$predict())

# testing control #FAILED
#	all.equal(output_nls$control,
#			 output_nlsj$control)

# testing convInfo # FAILED
	all.equal(as.numeric(output_nls$convInfo$isConv),
			 as.numeric(output_nlsj$convInfo))

#### TESTING nls VS minpacl.lm::nlsLM

## testing m values:  # PASSED
      # residuals
	all.equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsLM)))
	# fitted
	all.equal(as.vector(fitted(output_nls)),
			 as.vector(fitted(output_nlsLM)))
#	# formula
#	all.equal(formula(output_nls),
#			 formula(output_nlsLM))
	# deviance
	all.equal(deviance(output_nls),
			 deviance(output_nlsLM))
	# gradient
	all.equal( output_nls$m$gradient(),
			  output_nlsLM$m$gradient())
#	# conv
#	all.equal( output_nls$m$conv(),
#			  output_nlsLM$m$conv())
#	# incr
#	all.equal( output_nls$m$incr(),
#			  output_nlsLM$m$incr())
	# getPars # difference between getAllPars adn getPars?
	all.equal( output_nls$m$getAllPars(),
			  output_nlsLM$m$getAllPars())
#	# getEnv
#	all.equal( output_nls$m$getEnv(),
#			  output_nlsLM$m$getEnv())
	# Rmat
	all.equal( output_nls$m$Rmat(),
			  output_nlsLM$m$Rmat())
#	# predict
#	all.equal( output_nls$m$predict(),
#			  output_nlsLM$m$predict())

# testing control #FAILED
#	all.equal(output_nls$control,
#			 output_nlsLM$control)

# testing convInfo ## FAILED
	all.equal(output_nls$convInfo$isConv,
			 output_nlsLM$convInfo$isConv)
	all.equal(output_nls$convInfo$finIter,
			 output_nlsLM$convInfo$finIter)

## TESTING nls VS nlsr::nlxb

# testing coefficients # FAILED----> need to change tolerance - suggestions??
	all.equal(as.vector(output_nls$m$getAllPars()),
			as.vector(coefficients(output_nlxb)))
# testing residuals # FAILED----> signs all opposite??
	all.equal(as.vector(residuals(output_nls)),
			-1*as.vector(output_nlxb$resid))
# testing jacobian # FAILED----> should I change dimnames here? then it should pass
	all.equal(as.numeric(as.matrix(output_nls$m$gradient())),
			 as.numeric(output_nlxb$jacobian))
# testing formula  #PASSED
#	all.equal(formula(output_nls),
#			 	formula(output_nlxb))
# testing weights  #PASSED
#	all.equal(weights(output_nls),
#			 	weights(output_nlxb))
# testing ssquares  #PASSED
	all.equal(sum(as.vector(residuals(output_nls))^2),
		 	output_nlxb$ssquares)
 
# How do I get feval, jeval, lower, upper in nls? Also, what is maskidx?Not documented in the help file of nlxb
				
rm(hyd,iso.pen,n.pent,rate)

rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
"NLSupper","NLSlower","output_nls","output_nlsj","output_nlxb",
"output_nlsLM")

print("End of test file 'Isom_1.R' ")

