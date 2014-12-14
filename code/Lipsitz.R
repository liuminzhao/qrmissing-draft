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
