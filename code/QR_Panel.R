##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param X [n, xdim]
##' @param y [n, 2]
##' @param R [n] indicator if it 2nd component is missing, 1 is not missing, 0 is missing
##' @param w weights
##' @param taus quantiles of interest
##' @param lambda
##' @return
##' @author Minzhao Liu
QR.Panel <- function(X, y, R, w = c(.25,.5,.25), taus=(1:3)/4, lambda = 1){
# prototype function for panel data fitting of QR models
# the matrix X is assumed to contain an intercept
# the vector s is a strata indicator assumed (so far) to be a one-way layout
# NB:
# 1.  The value of the shrinkage parameter lambda is an open research problem in
# 	the simplest homogneous settings it should be the ratio of the scale parameters
# 	of the fixed effects and the idiocyncratic errors
# 2.  On return the coefficient vector has m*p + n elements where m is the number
#	quantiles being estimated, p is the number of colums of X, and n is the
#	number of distinct values of s.  The first m*p coefficients are the
#	slope estimates, and the last n are the "fixed effects"
# 3.  Like all shrinkage (regularization) estimators, asymptotic inference is somewhat
#	problematic... so the bootstrap is the natural first resort.

    n <- dim(y)[1]
    yc <- c(t(y))
    xdim <- dim(X)[2]

    oddind <- c(seq(1, n * 2 - 1, by = 2))
    evenind <- seq(2, n * 2, by = 2)

    Xstar <- matrix(0, n * 2, xdim * 2)
    Xstar[oddind, 1:xdim] <- X
    Xstar[evenind, (xdim + 1):(xdim * 2)] <- X

    Ry1 <- rep(1, n)
    Rc <- c(rbind(Ry1, R))

    yy <- yc[Rc == 1]
    XXstar <- Xstar[Rc == 1, ]

    s <- rep(1:n, R + 1)

    y <- yy
    X <- XXstar

    require(SparseM)
    require(quantreg)
    K <- length(w)
    if(K != length(taus)) stop("length of w and taus must match")
    X <- as.matrix(X)
    p <- ncol(X)
    n <- length(levels(as.factor(s)))
    N <- length(y)
    if(N != length(s) || N != nrow(X)) stop("dimensions of y,X,s must match")
    Z <- as.matrix.csr(model.matrix(~as.factor(s)-1))
    Fidelity <- cbind(as(w,"matrix.diag.csr") %x% X,as.matrix(w) %x% Z)
    Penalty <- cbind(as.matrix.csr(0,n,K*p),lambda*as(n,"matrix.diag.csr"))
    D <- rbind(Fidelity,Penalty)
    y <- c(w %x% y,rep(0,n))
    a <- c((w*(1-taus)) %x% (t(X)%*%rep(1,N)),
           sum(w*(1-taus)) * (t(Z) %*% rep(1,N)) + lambda * rep(1/2,n))

    mod <- rq.fit.sfn(D,y,rhs=a)

    coef <- matrix(mod$coef[1:(K*p)], p, K)
    coef[xdim + 1, ] <- coef[xdim + 1, ] + coef[1, ]
    intercept <- mod$coef[-(1:(K*p))]
    colnames(coef) <- taus

    return(list(coef = coef, intercept = intercept, mod = mod, taus = taus, X = X, y = y, s = s, w = w, lambda = lambda))
}
