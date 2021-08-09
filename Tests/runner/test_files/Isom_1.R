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
	
NLSdata <- data.frame(rate,hyd,iso.pen,n.pent)

## STARTING VALUE
b2 = 0.1
b3 = 0.1
b4 = 0.1
NLSstart <-c(b2 = b2, b3 = b3, b4 = b4) # a starting vector (named!)

## MODEL
NLSformula <- rate ~ b3*(n.pent - iso.pen/1.632)/(1+b2*hyd+b3*n.pent+b4*iso.pen)
NLSlower<- c(-Inf,-Inf,-Inf)
NLSupper<- c(Inf,Inf,Inf)
NLSweights <- rep(1,length(n.pent))
NLSsubset <- 1:length(n.pent)
rm(b2,b3,b4,rate,hyd,iso.pen,n.pent)
