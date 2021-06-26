# test name retention with as.numeric
# 161018
start <- c(one=1, two=2)
str(start)
s1 <- as.numeric(start)
str(s1)
names(s1)<-names(start)
str(s1)
strt2 <- as.matrix(start)
str(strt2)
s21 <- as.numeric(strt2)
s21
str(s21)
savehistory("tasnumeric.R")
