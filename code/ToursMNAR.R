#!/bin/Rscript
##' Time-stamp: <liuminzhao 07/27/2013 13:44:20>
##' 2013/07/05

ToursMNAR <- function(y, R, X, tau = 0.5, niter = 1000,  method="heter2",
                      init = NULL, tol = 0.00001){
  if (!is.loaded('toursmnarh2f')) {
    if (file.exists('toursmnar.so')) {
      dyn.load('toursmnar.so')
    } else {
      if (file.exists('~/Documents/qrmissing/code/toursmnar.so')) {
        dyn.load("~/Documents/qrmissing/code/toursmnar.so")
      } else stop('no shared library found.')
    }
  }
  n <- length(R)
  xdim <- dim(X)[2]
  sp <- rep(0, 2*xdim + 1)
  sp[1] <- 3.6/100
  if (is.null(init)) {
    require(quantreg)
    lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
    lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))
    if (tau == 0.05) {
      lmcoef1 <- c(0.080726637783397, 0.00451080465161577, -0.442936708038249, 0.753033829630134)
      lmcoef2 <- c(-0.206850522306966, -0.00232449257256374, -0.0183695937996462, 0.310703921030932)
    }
    if (tau == 0.95) {
      lmcoef1 <- c(0.112940628486934,  0.00540403561172756, -0.0692974826323355, 0.927829207235262)
      lmcoef2 <- c(0.217754080961376, -0.0026205522008995, -0.00553671984010173, 0.325975371558977)
    }
    param <- rep(0, 8*xdim + 3)
    param[1:xdim] <- lmcoef1
    param[(4*xdim + 1):(5*xdim)] <- lmcoef2
    param[(5*xdim + 1):(6*xdim)] = sp[1:xdim]
    param[(7*xdim + 1):(8*xdim)] = sp[(xdim + 2):(2*xdim + 1)]
    param[8*xdim + 2]= sp[xdim + 1]
    param[8*xdim + 3] = sum(R)/dim(X)[1]
    ## print(param)
    if (tau == 0.05) {
      param <- c(-0.000572618960159102, -0.00284943073315646, 0.0387886456476203,
                 0.731309636261029, -0.110853158691638, 0.0018238760530993, -0.0243176721533815,
                 0.179362520534912, -2.33848313046036, 0.0578751990667433, -0.551172968932799,
                 -0.134787034410123, -0.864886146000788, 0.0841551044842206, -0.374143122455584,
                 -0.872461326711233, -0.118910450477788, -0.00597614317844824,
                 0.0551509758919756, 0.827440939184248, 0.036, 0, 0, 0, -2.16130450605715,
                 0.0197336784662978, -0.225014982283925, -0.434830463466941, 0,
                 0, 0, 0, 1.18819130154664, 0, 0.94)
    }
    if (tau == 0.3) {
      param <- c(0.00207164856012197, 0.00348698431059567, -0.0434329086796143,
                 0.906494665358308, 0.0347767597803538, 0.00366260485541747, -0.0222134329792425,
                 -0.0174913111911812, -3.89434210328669, 0.00938039782029561,
                 0.211631973546578, 0.689616532836046, -4.21150770757033, 0.298592858625249,
                 1.93273639967165, -1.23939635465963, 0.052067291667809, -0.000102530353471596,
                 -0.0471642463757635, 0.861947813331854, 0.036, 0, 0, 0, -4.05935053112,
                 -0.00831585400869144, 0.279920867519652, 1.04348694211665, 0,
                 0, 0, 0, 1.15338496365044, 0, 0.94187458632803)
    }
    if (tau == 0.5) {
      param <- c(0.00994876514109072, 0.00419433932669911, -0.0383801425206726,
                 0.919951346365731, 0.0423430212405825, 0.00066663806352318, -0.0239859167896611,
                 -0.0283541244706559, -3.97942376952496, 0.00479280379443348,
                 0.229934102494423, 0.763430599426436, -1.02355045554959, 0.418650693153796,
                 0.429279425277303, -3.03175556803254, 0.049696448456005, 0.000342935729748995,
                 -0.0359943125714626, 0.898718393848919, 0.036, 0, 0, 0, -4.12006548592443,
                 -0.00561299217709275, 0.277345926929968, 1.1073039021098, 0,
                 0, 0, 0, 1.16365245384651, 0, 0.940214876125023)
    }
    if (tau == 0.7){
      param <- c(0.0116084984524679, 0.00477481026268049, -0.0346469683144935,
                 0.939746359895385, 0.0253558014316562, 0.00311442782148611, -0.0221453029998786,
                 -0.00924082301211957, -3.83950410842505, 0.00328283576575807,
                 0.216349315919138, 0.635663157238318, -4.10109210403393, 0.223481297781384,
                 2.26809120056051, -1.77503535231376, 0.0372299629093819, 0.000536074356143069,
                 -0.0263191599334046, 0.946838341736514, 0.036, 0, 0, 0, -3.98426874678104,
                 -0.00116401600961091, 0.242839644419077, 0.996126674727847, 0,
                 0, 0, 0, 1.15633120323484, 0, 0.941387714428302)
    }
    if (tau == 0.95) {
      param <- c(0.00991529110851534, 0.00309037437311279, -0.0183157387476969,
                 0.985564698371604, 0.163977788803744, 0.00798997538413755, -0.0364156667881455,
                 -0.140833953471624, -3.87054946009393, -0.0129390031193441, 0.28825171824574,
                 0.624904276468576, -2.60435609380152, -0.210512471011307, 0.122738068201498,
                 -1.41323625802395, 0.0675256965779804, -0.00153082235634583,
                 -0.00231632297479363, 0.991348196699425, 0.036, 0, 0, 0, -3.41404556047056,
                 -0.0076215084093807, 0.213233822326791, 0.449468768996346, 0,
                 0, 0, 0, 1.14318287144912, 0, 0.945803282246696)
    }
  } else param <- init

  paramsave <- matrix(0, niter, 8*xdim + 4)
  mod <- .Fortran("toursmnarh2f",
                  y = as.double(y),
                  R = as.integer(R),
                  x = as.double(X),
                  tau = as.double(tau),
                  n = as.integer(n),
                  niter = as.integer(niter),
                  param = as.double(param),
                  paramsave = as.double(paramsave),
                  xdim = as.integer(xdim),
                  converge = as.logical(TRUE),
                  tol = as.double(tol)
                  )
  mod$method <- method
  class(mod) <- "BiQRGradient"
  return(mod)
}


coef.BiQRGradient <- function(mod, ...){
  q <- mod$xdim
  param <- mod$param[1:(8*q)]
  coef <- matrix(param, 8, q, byrow = T)
  coef <- coef[c(1, 5, 2, 6, 3, 4, 7, 8), ]
  p <- mod$param[8*q + 3]
  beta22 <- mod$param[8*q + 1]
  h <- mod$param[8*q + 2]
  rownames(coef) <- c('Q1Coef', 'Q2Coef', 'beta1(0)', 'beta2(0)',
                      'Sigma1(1)', 'Sigma1(0)', 'Sigma2(1)', 'Sigma2(0)')
  return(coef[c(1, 2),])
}

print.BiQRGradient <- function(mod, ...){
  cat('Coefficients: \n')
  print(coef(mod))
}

summary.BiQRGradient <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau
  param <- mod$param
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', param[8*q + 3], '\n')
  cat('Quantile: ', tau, '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
}

plot.BiQRGradient <- function(mod, ...){
  a <- matrix(mod$paramsave, mod$niter, length(mod$param) + 1)
  xdim <- mod$xdim
  colnames(a) <- c(paste('gamma1', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma1', 1:xdim, sep = ''), paste('sigma0', 1:xdim, sep = ''), paste('gamma2', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma2', 1:xdim, sep = ''), paste('sp', 1:xdim, sep = ''), 'beta22', 'h', 'p', 'nll')
  par(mfrow = c(2, 2))
  for (i in 1:dim(a)[2]){
    plot(ts(a[, i]), main = colnames(a)[i])
  }
}

Diagnose <- function(mod){
  a <- matrix(mod$paramsave, mod$niter, length(mod$param) + 1)
  xdim <- mod$xdim
  colnames(a) <- c(paste('gamma1', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma1', 1:xdim, sep = ''), paste('sigma0', 1:xdim, sep = ''), paste('gamma2', 1:xdim, sep=''), paste('beta1', 1:xdim, sep=''), paste('sigma2', 1:xdim, sep = ''), paste('sp', 1:xdim, sep = ''), 'beta22', 'h', 'p', 'nll')
  par(mfrow = c(2, 2))
  for (i in 1:dim(a)[2]){
    plot(ts(a[, i]), main = colnames(a)[i])
  }
}
