##' 2014-04-12 An example showing that an interesting covariate efect (for example, a covariate with an effect only for upper quantiles) can be captured by the proposed model would strengthen the paper
##' from reviewer
##' borrow idea from sim-normal-mar.R

rm(list = ls())
library(qrmissing)
library(ggplot2)
set.seed(1)

###############
## PARAMETER
###############
n <- 1000
p <- 0.5
alpha <- 0

R <- rbinom(n, 1, p)
x <- runif(n, 0, 2)
y <- matrix(0, n, 2)
for (i in 1:n){
    if (R[i] == 1){
        y[i, 1] <- 2 + x[i] +(1 + alpha*x[i])*rnorm(1)
        y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + (1+alpha*x[i])*rnorm(1)
    } else {
        y[i, 1] <- -2 - x[i] +(1 + alpha*x[i])*rnorm(1)
        y[i, 2] <- 1 - x[i] - 0.5 * y[i, 1] + (1+alpha*x[i])*rnorm(1)
    }
}

X <- matrix(0, n, 2)
X[,1] <- 1
X[,2] <- x

mod1 <- QRMissingBi(y, R, X, tau = 0.1)
mod3 <- QRMissingBi(y, R, X, tau = 0.3)
mod5 <- QRMissingBi(y, R, X, tau = 0.5)
mod7 <- QRMissingBi(y, R, X, tau = 0.7)
mod9 <- QRMissingBi(y, R, X, tau = 0.9)

## plot
data <- data.frame(x = x, y1 = y[, 1], y2 = y[, 2])
coefmat <- rbind(coef(mod1), coef(mod3), coef(mod5), coef(mod7), coef(mod9))
coefmat <- cbind(coefmat, rep(seq(0.1, 0.9, 0.2), each = 2))
colnames(coefmat) <- c('intercept', 'slope', 'tau')
coefmat <- data.frame(as.matrix(coefmat))
coefmat$term <- rep(c('Y1', 'Y2'), 5)

## y1
y1plot <- ggplot(data, aes(x = x, y = y1)) + geom_point() + geom_abline(data=subset(coefmat, term == 'Y1'), aes(intercept = intercept, slope = slope, color = as.factor(tau)), show_guide = TRUE)

## y2
y2plot <- ggplot(data, aes(x = x, y = y2)) + geom_point() + geom_abline(data=subset(coefmat, term == 'Y2'), aes(intercept = intercept, slope = slope, color = as.factor(tau)), show_guide = TRUE)

png('diffquantiley1.png')
y1plot
dev.off()

png('diffquantiley2.png')
y2plot
dev.off()
