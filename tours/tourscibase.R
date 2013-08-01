##' tours data with covariates : age, race, baseline weight
##' 2013/06/23
##' 2013/07/29 using new initial and new ci

dat1 <- read.table('toursbootagebase-0729.txt')
dat2 <- read.table('toursbootagebase2.txt')
dat3 <- read.table('toursbootagebase3.txt')
dat4 <- read.table('toursbootagebase4.txt')
dat5 <- read.table('toursbootagebase5.txt')

dat <- rbind(dat1, dat2, dat3, dat4, dat5)


index1 <- which(dat1[, 41] == 1)
index3 <- which(dat1[,42] == 1)
index5 <- which(dat1[,43] == 1)
index7 <- which(dat1[,44] == 1)
index9 <- which(dat1[, 45] == 1)
coef1 <- dat1[index1, c(1:4, 21:24)]
coef3 <- dat1[index3, c(5:8, 25:28)]
coef5 <- dat1[index5, c(9:12, 29:32)]
coef7 <- dat1[index7, c(13:16, 33:36)]
coef9 <- dat1[index9, c(17:20, 37:40)]
datsummary1 <- apply(coef1, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary9 <- apply(coef9, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary3 <- apply(coef3, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary5 <- apply(coef5, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary7 <- apply(coef7, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

datest <- read.table('ageracebase.txt', header = T)

ci1 <- matrix(c(datsummary1), 2, 8, byrow = TRUE)
ci3 <- matrix(c(datsummary3), 2, 8, byrow = TRUE)
ci5 <- matrix(c(datsummary5), 2, 8, byrow = TRUE)
ci7 <- matrix(c(datsummary7), 2, 8, byrow = TRUE)
ci9 <- matrix(c(datsummary9), 2, 8, byrow = TRUE)

ci <- rbind(ci1[1, ], ci3[1, ], ci5[1,], ci7[1,], ci9[1,],
            ci1[2,] , ci3[2,], ci5[2,], ci7[2,], ci9[2,])

colnames(ci) <- c('Int.lo', 'Int.up', 'Age.lo', 'Age.up', 'White.lo', 'White.up', 'Base.lo', 'Base.up')

total <- cbind(datest, ci)[, c(1, 5,6,2,7,8,3,9,10,4,11,12)]

rownames(total) <- c('Y1 0.1', 'Y1 0.3', 'Y1 0.5', 'Y1 0.7', 'Y1 0.9',
                     'Y2 0.1', 'Y2 0.3', 'Y2 0.5', 'Y2 0.7', 'Y2 0.9')

library(xtable)
print(xtable(total))

print(total)
