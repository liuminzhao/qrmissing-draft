dat1 <- read.table('ageracebaseMixBayes0508-1.txt')
dat3 <- read.table('ageracebaseMixBayes0508-3.txt')
dat5 <- read.table('ageracebaseMixBayes0508-5.txt')
dat7 <- read.table('ageracebaseMixBayes0508-7.txt')
dat9 <- read.table('ageracebaseMixBayes0508-9.txt')

gamma1 <- rbind(dat1[1, ], dat3[1, ], dat5[1, ], dat7[1, ], dat9[1, ])
gamma2 <- rbind(dat1[2, ], dat3[2, ], dat5[2, ], dat7[2, ], dat9[2, ])

library(xtable)

xtable(gamma1)
xtable(gamma2)
