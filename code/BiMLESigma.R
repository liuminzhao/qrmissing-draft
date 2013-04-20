## Time-stamp: <liuminzhao 04/20/2013 18:19:28>
## WRAP UP BiMLESigma.f

dyn.load('BiMLESigma.so')
dyn.load('BiMLESigmaH1.so')
dyn.load('BiMLESigmaH2.so')
BiQRGradient <- function(y, R, X, tau=0.5, niter = 1000, sp = rep(0, 1 + 2*dim(X)[2]), method = 'homo'){
  n <- length(R)
  if (method == 'homo') {
    param <- rep(0, 14)
    param[1] = 0 
    param[2] = 0 
    param[3] = 0 
    param[4] = 0 
    param[5] = 1 
    param[6] = 1 
    param[7] = 0 
    param[8] = 0 
    param[9] = sp[1]
    param[10]= sp[2]
    param[11]= sp[3] 
    param[12]= 1 
    param[13]= sp[4]
    param[14]= 0.5 

    paramsave <- matrix(0, niter, 15)
    mod <- .Fortran("BiQRGradientf",
                    y = as.double(y),
                    R = as.integer(R),
                    x = as.double(x),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave)
                    )
  } else if (method == 'heter1') {
    param <- rep(0, 18)
    param[1] = 0 
    param[2] = 0 
    param[3] = 0 
    param[4] = 0 
    param[5] = 1 
    param[6] = 1 
    param[7] = 0 
    param[8] = 0 
    param[9] = sp[1]
    param[10]= sp[2]
    param[11]= sp[3] 
    param[12]= 1 
    param[13]= sp[4]
    param[14]= 0.5 
    param[15]= 0
    param[16]= 0
    param[17]= 0
    param[18]= sp[5]

    paramsave <- matrix(0, niter, 19)
    mod <- .Fortran("BiQRGradientH1f",
                    y = as.double(y),
                    R = as.integer(R),
                    x = as.double(x),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave)
                    )
  } else if (method == 'heter2') {
    xdim <- dim(X)[2]
    lmcoef <- solve(t(X)%*%X)%*%(t(X)%*%y)
    param <- rep(0, 8*xdim + 3)
    param[1:xdim] <- lmcoef[, 1]
    param[(4*xdim + 1):(5*xdim)] <- lmcoef[, 2]
    param[(5*xdim + 1):(6*xdim)] = sp[1:xdim]
    param[(7*xdim + 1):(8*xdim)] = sp[(xdim + 2):(2*xdim + 1)]
    param[8*xdim + 2]= sp[xdim + 1]
    param[8*xdim + 3] = 0.5

    paramsave <- matrix(0, niter, 8*xdim + 4)
    mod <- .Fortran("BiQRGradientH2f",
                    y = as.double(y),
                    R = as.integer(R),
                    x = as.double(X),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave),
                    xdim = as.integer(xdim)
                    )
  }
  return(mod)
}

Diagnose <- function(mod){
  a <- matrix(mod$paramsave, mod$niter, 19)
  colnames(a) <- c('Gamma0', 'Gamma1', 'Beta0', 'Beta1', 'Sigma1', 'Sigma0', 'Gamma02', 'Gamma12', 'beta02','beta12','beta22','sigma2','lambda','phi', 'heter1(1)', 'heter1(2)', 'heter2(1)', 'heter2(2)','NLL')
  par(mfrow = c(2, 2))
  for (i in 1:19){
    plot(ts(a[, i]), main = colnames(a)[i])
  }
}

mysummary <- function(mod){
  q <- mod$xdim
  param <- mod$param[1:(8*q)]
  coef <- matrix(param, 8, q, byrow = T)
  coef <- coef[c(1, 5, 2, 6, 3, 4, 7, 8), ]
  p <- mod$param[8*q + 3]
  beta22 <- mod$param[8*q + 1]
  h <- mod$param[8*q + 2]
  rownames(coef) <- c('Q1Coef', 'Q2Coef', 'beta1(0)', 'beta2(0)',
                      'Sigma1(1)', 'Sigma1(0)', 'Sigma2(1)', 'Sigma2(0)')
  cat("\n Coefficients: \n")
  print(coef)
  cat("\n beta22 is ", beta22, " h is ", h, " p is ", p, "\n")
}
