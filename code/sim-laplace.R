#!/bin/Rscript
##' simulation laplace distribution
##' 2013/06/27

rm(list = ls())
source('BiMLESigma.R')
source('sendEmail.R')
library(quantreg)
library(xtable)
library(doMC)
registerDoMC()
options(cores = 8)
set.seed(1)
