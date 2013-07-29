#!/bin/Rscript
##' Time-stamp: <liuminzhao 07/29/2013 05:19:25>
##' manipulate data TOURS
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate

rm(list=ls())
TOURS <- read.csv('~/Documents/qrmissing/tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
trt <- TOURS$TREATMENT
age_center <- (age-mean(age))/sd(age)
race3 <- as.numeric(TOURS$RACE == 3)

## center weight2?
## weight2 <- scale(weight2)
## weight3 <- scale(weight3)
weight1 <- weight1/100
weight2 <- weight2/100
weight3 <- weight3/100


n <- length(age)
y <- matrix(0, n, 2)
y[,1] <- weight2
y[,2] <- weight3


X <- matrix(0, n, 4)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3
X[,4] <- weight1

R <- 1 - as.numeric(is.na(weight3))

dat <- data.frame(weight2, weight3, trt, age_center, age=TOURS$AGE, race = factor(TOURS$RACE), base = weight1)

## PLOT
## library(ggplot2)
## library(gridExtra)
## ggplot(data = dat, aes(x = age_center, y = weight2)) + geom_point()
## ggplot(data = dat, aes(x = age_center, y = weight3)) + geom_point()

## box1 <- ggplot(data = dat, aes(x = age, y = weight2*100)) + geom_point()+ ylab('Weight (Kg) at 6 months') + xlab('Age')
## box2 <- ggplot(data = dat, aes(x = age, y = weight3*100)) + geom_point()+ ylab('Weight (Kg) at 18 months') + xlab('Age')
## box3 <- ggplot(data = dat, aes(x = race, y = weight2*100)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight (Kg) at 6 months') + xlab('Race')
## box4 <- ggplot(data = dat, aes(race, y = weight3*100)) + geom_boxplot() + ylab('Weight (Kg) at 18 months') + xlab('Race') + scale_x_discrete(labels=c("Black", "White"))
## box5 <- ggplot(data = dat, aes(x = weight1*100, y = weight2*100)) + geom_point() + ylab('Weight (Kg) at 6 months') + xlab('Weight (Kg) at baseline')
## box6 <- ggplot(data = dat, aes(x = weight1*100, y = weight3*100)) + geom_point() + ylab('Weight (Kg) at 18 months') + xlab('Weight (Kg) at baseline')

## scat1 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat2 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()
## scat3 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat4 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()

## pdf('weight-age-race-base.pdf')
## sds <- grid.arrange(box1, box2, box3, box4, box5, box6, nrow = 2, ncol = 3)
## dev.off()


###############
## ANALYSIS
###############

y[is.na(y[,2]),2] <- 0

## initial values
init9 <- c(0.0198771492948432, 0.00419168481637642, -0.0286769281833563,
0.962705938424494, 0.124463473412154, 0.0082669321531808, -0.0332265872488645,
-0.101252841224114, -4.00482423046243, 0.00384361643406015, 0.242014360378546,
0.774992605720537, -4.41878993116355, -0.0893229764233047, 0.235700725875636,
0.367955520969719, 0.0365627649996702, -0.000848529569421987,
-0.0131870438707394, 1.00131102458501, 0.036, 0, 0, 0, -4.12227104972579,
0.00335677750788003, 0.245472347196181, 1.13619092057686, 0,
0, 0, 0, 1.15785690811763, 0, 0.943612829565736)

init1 <- c(0.00128500849891446, 0.00165353550525411, -0.0457700527122326,
0.852903241638371, 0.0965964850319247, 0.0106681055898985, -0.00833488917465255,
-0.0983272783788632, -3.97083966918167, 0.0224131103967352, 0.143912587571141,
0.829355087608916, -2.06231973874536, -0.197632882886333, -1.70344194695517,
-0.0974725348981226, 0.0782115557466507, -0.00109541725928936,
-0.0654920860767664, 0.755271619906803, 0.036, 0, 0, 0, -4.2102888327682,
-0.00912613468204781, 0.286319221594694, 1.19375723405397, 0,
0, 0, 0, 1.14599024095279, 0, 0.946162518243309)

init3 <- c(0.00346754533524767, 0.00340547757160402, -0.0424477368487844,
0.903304542638568, 0.232602165858655, 0.0066085533754103, -0.047828380774495,
-0.205949076863877, -3.98127225565636, 0.00826563282087521, 0.202449750934045,
0.788235833662693, -3.5365921209912, -0.143132715223207, 0.792684508025302,
-1.21289657040325, 0.0522199852142105, -0.000783994978801458,
-0.0436083360986863, 0.856314796337822, 0, 0, 0, 0, -4.33162748420352,
-0.0108090053648665, 0.327970623481113, 1.28007377111706, 0,
0, 0, 0, 1.14599514929207, 0, 0.937596514557817)

init5 <- c(0.00817699838851987, 0.00330942186819275, -0.0399131709475386,
0.922897301066523, 0.126750584163812, 0.0058020299150031, -0.0355863273100943,
-0.102391762777324, -4.01950753497282, 0.00628486521066096, 0.243707932246988,
0.792047381377232, -3.28627563586285, 0.0348857118359061, 0.355686812534977,
-0.895834755931226, 0.0489958299031745, -0.00101797992303196,
-0.0362702460564401, 0.898227827546846, 0, 0, 0, 0, -4.2666036648914,
-0.00736395189520582, 0.307366569878983, 1.23213379914368, 0,
0, 0, 0, 1.15359495226242, 0, 0.941404809334169)

init7 <- c(0.00665810647810459, 0.00390179829308831, -0.0355819879382891,
0.946502966557501, 0.151101151390684, 0.00565476504014295, -0.0395569668320582,
-0.125689555220142, -3.93441276791615, 0.00427008684110161, 0.233183683079389,
0.714398756168959, -3.53017373626221, 0.0480494948493253, 0.575154747202667,
-0.909814475356364, 0.0394362427537354, -0.000588417299503833,
-0.0268498352569201, 0.944566617750712, 0, 0, 0, 0, -4.07550686393392,
-0.00479128863405338, 0.276711851665603, 1.06053661934902, 0,
0, 0, 0, 1.15556927936215, 0, 0.9391757154398)

source('~/Documents/qrmissing/code/BiMLESigma.R')
mod1 <- BiQRGradient(y, R, X, tau = 0.05, niter = 2000,  method = 'heter2', init = init1)
mod3 <- BiQRGradient(y, R, X, tau = 0.3, niter = 2000,  method = 'heter2', init = init3)
mod5 <- BiQRGradient(y, R, X, tau = 0.5, niter = 2000,  method = 'heter2', init = init5)
mod7 <- BiQRGradient(y, R, X, tau = 0.7, niter = 2000,  method = 'heter2', init = init7)
mod9 <- BiQRGradient(y, R, X, tau = 0.95, niter = 2000,  method = 'heter2', init = init9)

coef1 <- coef(mod1)
coef3 <- coef(mod3)
coef5 <- coef(mod5)
coef7 <- coef(mod7)
coef9 <- coef(mod9)

coefw2 <- rbind(coef1[1,], coef3[1, ], coef5[1, ], coef7[1, ], coef9[1,])
coefw3 <- rbind(coef1[2,], coef3[2, ], coef5[2, ], coef7[2, ], coef9[2,])

rownames(coefw2) <- c('tau = 0.05', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.95')
rownames(coefw3) <- c('tau = 0.05', 'tau = 0.3', 'tau = 0.5', 'tau = 0.7', 'tau = 0.95')

colnames(coefw2) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')
colnames(coefw3) <- c('Intercept', 'Age(centered)', 'White', 'BaseWeight')

library(xtable)
print(xtable(coefw2))
print(xtable(coefw3))

write.table(rbind(coefw2, coefw3), 'ageracebase.txt', row.names=FALSE)
