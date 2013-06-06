dat1 <- read.table('toursbootage.txt')
dat2 <- read.table('toursbootage2.txt')
dat3 <- read.table('toursbootage3.txt')
dat4 <- read.table('toursbootage4.txt')
dat5 <- read.table('toursbootage5.txt')

dat <- rbind(dat1, dat2, dat3, dat4, dat5)

datsummary <- apply(dat, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

datest <- read.table('ageresult.txt', header = T)

ci <- matrix(c(datsummary), 10, 6, byrow = TRUE)
colnames(ci) <- c('Int.lo', 'Int.up', 'Age.lo', 'Age.up', 'White.lo', 'White.up')

total <- cbind(datest, ci)[, c(1, 4,5,2,6,7,3,8,9)]

library(xtable)
print(xtable(total))
