##' tours data with covariates : age, race, baseline weight
##' 2013/06/23

dat1 <- read.table('toursbootagebbase.txt')
dat2 <- read.table('toursbootagebbase2.txt')
dat3 <- read.table('toursbootagebbase3.txt')
dat4 <- read.table('toursbootagebbase4.txt')
dat5 <- read.table('toursbootagebbase5.txt')

dat <- rbind(dat1, dat2, dat3, dat4, dat5)

datsummary <- apply(dat, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

datest <- read.table('ageracebase.txt', header = T)

ci <- matrix(c(datsummary), 10, 2*dim(datest)[2], byrow = TRUE)
colnames(ci) <- c('Int.lo', 'Int.up', 'Age.lo', 'Age.up', 'White.lo', 'White.up', 'Base.lo', 'Base.up')

total <- cbind(datest, ci)[, c(1, 5,6,2,7,8,3,9,10,4,11,12)]

library(xtable)
print(xtable(total))
