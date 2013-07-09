#!/bin/Rscript
##' Time-stamp: <liuminzhao 07/09/2013 05:07:28>
##' 2013/07/09 Bottai's algorithm
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Bottai's algorithm for bivariate data
##' @param y
##' @param R
##' @param X
##' @param tau
##' @param M
##' @return
##' @author Minzhao Liu
Bottai <- function(y, R, X, tau = 0.5, M = 5){
  require(quantreg)
  ans <- matrix(0, ncol(X), 2)
  for (i in 1:M) {
    mod1 <- rq(y[, 1] ~ X[, -1], tau = tau)
    ans[, 1] <- ans[, 1] + coef(mod1)
    index <- which(R == 0)
    for (j in 1:length(index)){
      u <- runif(1)
      mod21 <- rq(y[R==1, 2] ~ X[R==1, -1] + y[R==1, 1], tau = u)
      coef21 <- coef(mod21)
      y[index[j], 2] <- (c(X[index[j],], y[index[j],1]))%*%coef21
    }
    mod2 <- rq(y[, 2] ~ X[, -1], tau = tau)
    ans[, 2] <- ans[, 2] + coef(mod2)
  }
  ans <- ans / M
  colnames(ans) <- c('Y1', 'Y2')
  return(ans)
}
