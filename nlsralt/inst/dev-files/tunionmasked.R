prm <- c(one=1, two=2, 3=3, 4=4)
prm <- c(one=1, two=2, "3"=3, "4"=4)
prm
maskidx <- c(1)
maskidx <- c("one")
lower <- c(0,0,0,1)
upper <- c(1,1,1,1)
tmask <- which(lower==upper)
tmask
newmask <- maskidx && tmask
newmask <- maskidx & tmask
tmask
str(tmask)
str(maskidx)
pnames <- names(prm)
pnames
masked <- maskidx # fix this here -- wrong var name
masked
maskidx <- which(pnames %in% masked)
maskidx
newmask <- maskidx & tmask
newmask
# not the right type
newmask <- union(tmask, maskidx)
newmask
union(tmask, tmask)
savehistory("tunionmasked.R")
