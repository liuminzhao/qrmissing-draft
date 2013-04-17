## Time-stamp: <liuminzhao 04/16/2013 22:56:58>
## WRAP UP BiMLESigma.f

dyn.load('BiMLESigma.so')
dyn.load('BiMLESigmaH1.so')
dyn.load('BiMLESigmaH2.so')
BiQRGradient <- function(y, R, x, tau=0.5, niter = 1000, sp = c(0,0,0,1,1), method = 'homo'){
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
    param <- rep(0, 19)
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
    param[11]= 0
    param[12]= 1 
    param[13]= sp[4]
    param[14]= 0.5 
    param[15]= 0
    param[16]= 0
    param[17]= 0
    param[18]= sp[5]
    param[19] = sp[3]

    paramsave <- matrix(0, niter, 20)
    mod <- .Fortran("BiQRGradientH2f",
                    y = as.double(y),
                    R = as.integer(R),
                    x = as.double(x),
                    tau = as.double(tau),
                    n = as.integer(n),
                    niter = as.integer(niter),
                    param = as.double(param),
                    paramsave = as.double(paramsave)
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
