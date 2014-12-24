
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param X [n, p] : same for same subject
##' @param y [n, 2] : y[, 2] maybe missing
##' @param R [n] indicator if it 2nd component is missing, 1 is not missing, 0 is missing
##' @param tau quantile of interest
##' @return
##' @author Minzhao Liu
Lipsitz <- function(X, y, R, tau) {
    n <- dim(y)[1]
    y1 <- y[, 1]
    xdim <- dim(X)[2]

    mod <- glm(R ~ y1 + X[, -1], family = "binomial")

    eta <- fitted(mod)

    pidi <- R * eta + (1 - R) * (1 - eta)

    ystar <- y/pidi
    Xstar <- X
    Xstar[, -1] <- X[, -1]/pidi
    oddind <- c(seq(1, n * 2 - 1, by = 2))
    evenind <- seq(2, n * 2, by = 2)

    ## new X2 design matrix to allow different coef for each components, X2[n*2, p*2]
    Xstarc <- Xstar
    Xstarc2 <- matrix(0, n * 2, xdim * 2)
    Xstarc2[oddind, 1:xdim] <- Xstarc
    Xstarc2[evenind, (xdim + 1):(xdim * 2)] <- Xstarc
    Xstarc2[, 1] <- 1

    ## new y
    ystarc <- c(t(ystar))

    ## new R, for subset condition
    Ry1 <- rep(1, n)
    Rc <- c(rbind(Ry1, R))

    mod <- rq(ystarc ~ Xstarc2[, -1], tau = tau, subset = Rc == 1)

    ## coef
    coef <- coef(mod)
    coef[xdim + 1, ] <- coef[xdim + 1, ] + coef[1, ]

    gamma1 <- coef[1:xdim]
    gamma2 <- coef[(xdim + 1): (xdim * 2)]
    gamma2[1] <- gamma2[1] + gamma1[1]
    return(list(mod = mod, gamma1 = gamma1, gamma2 = gamma2, coef = coef))
}
