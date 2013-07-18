##' Time-stamp: "liuminzhao 07/17/2013 17:10:12"
##' ##' tours data with covariates : age, race, baseline weight
##' 2013/07/14 MNAR

dat1 <- read.table('toursbootagebbasemnar.txt')
dat2 <- read.table('toursbootagebbasemnar2.txt')
dat3 <- read.table('toursbootagebbasemnar3.txt')
dat4 <- read.table('toursbootagebbasemnar4.txt')
dat5 <- read.table('toursbootagebbasemnar5.txt')

dat <- rbind(dat1, dat2, dat3, dat4, dat5)

datsummary <- apply(dat, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

datest <- read.table('ageracebasemnar.txt', header = T)

ci <- matrix(c(datsummary), 10, 2*dim(datest)[2], byrow = TRUE)
colnames(ci) <- c('Int.lo', 'Int.up', 'Age.lo', 'Age.up', 'White.lo', 'White.up', 'Base.lo', 'Base.up')

total <- cbind(datest, ci)[, c(1, 5,6,2,7,8,3,9,10,4,11,12)]

library(xtable)
print(xtable(total))

## for mod1 and mod9 only

dat1 <- read.table('toursbootagebbasemnar-19.txt')
dat2 <- read.table('toursbootagebbasemnar-19-2.txt')
dat <- rbind(dat1, dat2)
apply(dat[, 17:18], 2, sum)

index1 <- which(dat[,17] == 1)
index2 <- which(dat[,18] == 1)
coef1 <- dat[index1, c(1:4, 9:12)]
coef2 <- dat[index2, c(5:8, 13:16)]

datsummary1 <- apply(coef1, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary2 <- apply(coef2, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

datest <- read.table('ageracebasemnar.txt', header = T)[c(1, 5, 6, 10),]

ci1 <- matrix(c(datsummary1), 2, 8, byrow = TRUE)
ci2 <- matrix(c(datsummary2), 2, 8, byrow = TRUE)

ci <- rbind(ci1[1, ], ci2[1, ], ci1[2,] , ci2[2,])
colnames(ci) <- c('Int.lo', 'Int.up', 'Age.lo', 'Age.up', 'White.lo', 'White.up', 'Base.lo', 'Base.up')

total19 <- cbind(datest, ci)[, c(1, 5,6,2,7,8,3,9,10,4,11,12)]

rownames(total19) <- c('Y1 0.1', 'Y1 0.9', 'Y2 0.1', 'Y2 0.9')

print(xtable(total19))

print(total19)
