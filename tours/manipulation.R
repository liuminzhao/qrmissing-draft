## Time-stamp: <liuminzhao 04/20/2013 12:04:17>
## manipulate data TOURS
## 2012/06/06 add age.center
rm(list=ls())
TOURS <- read.csv('../tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
age_center <- (age-mean(age))/sd(age)
race3 <- as.numeric(TOURS$RACE == 3)

## change <- weight2-weight1
change1 <- weight1 - weight2
change2 <- weight2 - weight3
change3 <- weight1 - weight3

n <- length(change1)
y <- matrix(0, n, 2)
y[,1] <- change1
y[,2] <- change2


X <- matrix(0, n, 3)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3

R <- 1 - as.numeric(is.na(change2))

dat <- data.frame(change1, change2, change3, age_center, race = factor(TOURS$RACE))

## PLOT
library(ggplot2)
ggplot(data = dat, aes(x = age_center, y = change1)) + geom_point()
ggplot(data = dat, aes(x = age_center, y = change2)) + geom_point()
ggplot(data = dat, aes(x = age_center, y = change3)) + geom_point()


scat1 <- ggplot(data = dat, aes(x = age_center, y = change1, color = race)) + geom_point()
scat2 <- ggplot(data = dat, aes(x = age_center, y = change2, color = race)) + geom_point()
scat3 <- ggplot(data = dat, aes(x = age_center, y = change3, color = race)) + geom_point()

box1 <- ggplot(data = dat, aes(race, y = change1)) + geom_boxplot()
box2 <- ggplot(data = dat, aes(race, y = change2)) + geom_boxplot()
box3 <- ggplot(data = dat, aes(race, y = change3)) + geom_boxplot()

sds <- grid.arrange(scat1, box1, scat2, box2, scat3, box3, nrow = 3, ncol = 2)

pdf('tours.pdf')
sds <- grid.arrange(scat1, box1, scat2, box2, scat3, box3, nrow = 3, ncol = 2)
dev.off()

## ANALYSIS
y[is.na(y[,2]),2] <- 0

source('BiMLESigma.R')
mod1 <- BiQRGradient(y, R, X, tau = 0.1, method = 'heter2')
mod3 <- BiQRGradient(y, R, X, tau = 0.3, method = 'heter2')
mod5 <- BiQRGradient(y, R, X, tau = 0.5, method = 'heter2')
mod7 <- BiQRGradient(y, R, X, tau = 0.7, method = 'heter2')
mod9 <- BiQRGradient(y, R, X, tau = 0.9, method = 'heter2')

plot(change1 ~ age_center)
abline(mod1$param[1:2], col = 'red')
abline(mod3$param[1:2], col = 'red')
abline(mod5$param[1:2], col = 'red')
abline(mod7$param[1:2], col = 'red')
abline(mod9$param[1:2], col = 'red')

modrq = rq(change1 ~ age_center + factor(race), tau = seq(0.1,0.9,len =5))

for (i in 1:5)
  abline(modrq$coef[-3, i], col = 'green')


plot(change2 ~ age_center)
abline(mod1$param[13:14], col = 'red')
abline(mod3$param[13:14], col = 'red')
abline(mod5$param[13:14], col = 'red')
abline(mod7$param[13:14], col = 'red')
abline(mod9$param[13:14], col = 'red')

modrq2 = rq(change2[R==1] ~ age_center[R==1] + factor(race)[R==1], tau = seq(0.1,0.9,len =5))

for (i in 1:5)
  abline(modrq2$coef[-3, i], col = 'green')

## RESULTS

result1 <- cbind(mod1$param[c(1:3, 13:15)], mod3$param[c(1:3, 13:15)],
                 mod5$param[c(1:3, 13:15)], mod7$param[c(1:3, 13:15)],
                 mod9$param[c(1:3, 13:15)])
      
resultrq <- rbind(modrq$coef, modrq2$coef)

library(xtable)
print(xtable(result1))
print(xtable(resultrq))
write.table(result1, '../tours/result1.txt')
write.table(resultrq, '../tours/resultrq.txt')
