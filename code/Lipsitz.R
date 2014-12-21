
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

    mod <- glm(R ~ y1 + X[, -1], family = "binomial")

    eta <- fitted(mod)

    pidi <- R * eta + (1 - R) * (1 - eta)

    ystar <- y/pidi
    Xstar <- X
    Xstar[, -1] <- X[, -1]/pidi
    ind <- c(seq(1, n * 2 - 1, by = 2), seq(2, n * 2, by = 2))

    Xstarc <- rbind(Xstar, Xstar)
    Xstarc2[ind, ] <- Xstarc

    ystarc <- c(t(ystar))
    Ry1 <- rep(1, n)
    Rc <- c(rbind(Ry1, R))

    mod <- rq(ystarc ~ Xstarc2[, -1], tau = tau, subset = Rc == 1)
    return(mod)
}
