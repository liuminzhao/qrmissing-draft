##' Time-stamp: "liuminzhao 07/28/2013 23:14:27"
##' ##' tours data with covariates : age, race, baseline weight
##' 2013/07/14 MNAR
##' 2013/07/27 new mnar
##' new mnar with 200 boot for mod1 and 200 for mod9
##' 2013/07/28

dat11 <- read.table('toursbootagebasemnar-0727-1.txt')
dat12 <- read.table('toursbootagebasemnar-0727-1-2.txt')
dat13 <- read.table('toursbootagebasemnar-0727-1-3.txt')
dat91 <- read.table('toursbootagebasemnar-0727-9.txt')
dat92 <- read.table('toursbootagebasemnar-0727-9-2.txt')
dat93 <- read.table('toursbootagebasemnar-0727-9-3.txt')
dat94 <- read.table('toursbootagebasemnar-0727-9-4.txt')
dat95 <- read.table('toursbootagebasemnar-0727-9-5.txt')
dat1 <- rbind(dat11, dat12, dat13)
dat9 <- rbind(dat91, dat92, dat93, dat94, dat95)

sum(dat1[,9])
sum(dat9[,9])

index1 <- which(dat1[, 9] == 1)
index9 <- which(dat9[, 9] == 1)
coef1 <- dat1[index1, 1:8]
coef9 <- dat9[index9, 1:8]
datsummary1 <- apply(coef1, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary9 <- apply(coef9, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

dat37 <- read.table('toursbootagebasemnar-0724-3-7.txt')

apply(dat37[, 25:27], 2, sum)

index3 <- which(dat37[,25] == 1)
index5 <- which(dat37[,26] == 1)
index7 <- which(dat37[,27] == 1)
coef3 <- dat37[index3, c(1:4, 13:16)]
coef5 <- dat37[index5, c(5:8, 17:20)]
coef7 <- dat37[index7, c(9:12, 21:24)]

datsummary3 <- apply(coef3, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary5 <- apply(coef5, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
datsummary7 <- apply(coef7, 2, function(x) quantile(x, probs = c(0.025, 0.975)))


datest <- read.table('ageracebasemnar.txt', header = T)

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
