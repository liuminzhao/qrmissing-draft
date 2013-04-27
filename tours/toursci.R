dat1 <- read.table('toursboot.txt')
dat2 <- read.table('toursboot2.txt')
dat3 <- read.table('toursboot3.txt')
dat4 <- read.table('toursboot4.txt')
dat5 <- read.table('toursboot5.txt')
dat6 <- read.table('toursboot6.txt')

dat <- rbind(dat1, dat2, dat3, dat4, dat5, dat6)

datsummary <- apply(dat, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

datest <- read.table('weightresult.txt', header = T)

ci <- matrix(c(datsummary), 10, 8, byrow = TRUE)
colnames(ci) <- c('Int.lo', 'Int.up', 'Trt.O.lo', 'Trt.O.up', 'Trt.T.lo', 'Trt.T.up', 'Race.3.lo', 'Race.3.up')

total <- cbind(datest, ci)[, c(1, 5,6,2,7,8,3,9,10, 4, 11,12)]

library(xtable)
print(xtable(total))
