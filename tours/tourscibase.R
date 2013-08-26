##' tours data with covariates : age, race, baseline weight
##' Time-stamp: <liuminzhao 08/26/2013 08:27:56>
##' 2013/06/23
##' 2013/07/29 using new initial and new ci
##' 2013/08/13 new method from  uobyqa
##' 2013/08/26 qrmissing package

dat <- read.table('toursbootageracebase-0825.txt')

index1 <- which(dat[,41] == 0)
index3 <- which(dat[,42] == 0)
index5 <- which(dat[,43] == 0)
index7 <- which(dat[,44] == 0)
index9 <- which(dat[,45] == 0)
coef1 <- dat[index1, c(1:4, 21:24)]
coef3 <- dat[index3, c(5:8, 25:28)]
coef5 <- dat[index5, c(9:12, 29:32)]
coef7 <- dat[index7, c(13:16, 33:36)]
coef9 <- dat[index9, c(17:20, 37:40)]
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

rownames(total) <- c('Y10.1', 'Y10.3', 'Y10.5', 'Y10.7', 'Y10.9',
                     'Y20.1', 'Y20.3', 'Y20.5', 'Y20.7', 'Y20.9')

## total[, -c(4:6, 9:12)] <- total[, -c(4:6, 9:12)]*100

library(xtable)
print(xtable(total))

print(total)
