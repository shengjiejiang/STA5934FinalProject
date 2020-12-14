library(pryr)
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(data.table)
sourceCpp("twas_perm.cpp")
source("aSPUwscore.R")
source("ACAT.R")
source("dist_support.R")
source("twas_support.R")
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))
source("aSPUwscore.R")
source("ACAT.R")
source("dist_support.R")
source("twas_support.R")

U1 <- rep(NA, 10000)
U2 <- rep(NA, 10000)
U3 <- rep(NA, 10000)
X1 <- rep(NA, 10000)
X2 <- rep(NA, 10000)
X3 <- rep(NA, 10000)
Y11 <- rep(NA, 10000)
Y12 <- rep(NA, 10000)
Y13 <- rep(NA, 10000)
Y21 <- rep(NA, 10000)
Y22 <- rep(NA, 10000)
Y23 <- rep(NA, 10000)

betaXG1 <- rep(NA, 25)
betaXG2 <- rep(NA, 25)
betaXG3 <- rep(NA, 25)
sebetaXG1 <- rep(NA, 25)
sebetaXG2 <- rep(NA, 25)
sebetaXG3 <- rep(NA, 25)
betaYG11 <- rep(NA, 25)
betaYG12 <- rep(NA, 25)
betaYG21 <- rep(NA, 25)
betaYG22 <- rep(NA, 25)
betaYG31 <- rep(NA, 25)
betaYG32 <- rep(NA, 25)

sebetaYG11 <- rep(NA, 25)
sebetaYG12 <- rep(NA, 25)
sebetaYG21 <- rep(NA, 25)
sebetaYG22 <- rep(NA, 25)
sebetaYG31 <- rep(NA, 25)
sebetaYG32 <- rep(NA, 25)

beta_median1 <- rep(NA, 5000)
beta_median2 <- rep(NA, 5000)
beta_median3 <- rep(NA, 5000)
beta_median4 <- rep(NA, 5000)
beta_median5 <- rep(NA, 5000)
beta_median6 <- rep(NA, 5000)

sebeta_median1 <- rep(NA, 5000)
sebeta_median2 <- rep(NA, 5000)
sebeta_median3 <- rep(NA, 5000)
sebeta_median4 <- rep(NA, 5000)
sebeta_median5 <- rep(NA, 5000)
sebeta_median6 <- rep(NA, 5000)

betaWM1 <- rep(NA, 5000)
betaWM2 <- rep(NA, 5000)
betaWM3 <- rep(NA, 5000)
betaWM4 <- rep(NA, 5000)
betaWM5 <- rep(NA, 5000)
betaWM6 <- rep(NA, 5000)

sebetaWM1 <- rep(NA, 5000)
sebetaWM2 <- rep(NA, 5000)
sebetaWM3 <- rep(NA, 5000)
sebetaWM4 <- rep(NA, 5000)
sebetaWM5 <- rep(NA, 5000)
sebetaWM6 <- rep(NA, 5000)


betaPWM1 <- rep(NA, 5000)
betaPWM2 <- rep(NA, 5000)
betaPWM3 <- rep(NA, 5000)
betaPWM4 <- rep(NA, 5000)
betaPWM5 <- rep(NA, 5000)
betaPWM6 <- rep(NA, 5000)


sebetaPWM1 <- rep(NA, 5000)
sebetaPWM2 <- rep(NA, 5000)
sebetaPWM3 <- rep(NA, 5000)
sebetaPWM4 <- rep(NA, 5000)
sebetaPWM5 <- rep(NA, 5000)
sebetaPWM6 <- rep(NA, 5000)

betaIVW1 <- rep(NA, 5000)
betaIVW2 <- rep(NA, 5000)
betaIVW3 <- rep(NA, 5000)
betaIVW4 <- rep(NA, 5000)
betaIVW5 <- rep(NA, 5000)
betaIVW6 <- rep(NA, 5000)

sebetaIVW1 <- rep(NA, 5000)
sebetaIVW2 <- rep(NA, 5000)
sebetaIVW3 <- rep(NA, 5000)
sebetaIVW4 <- rep(NA, 5000)
sebetaIVW5 <- rep(NA, 5000)
sebetaIVW6 <- rep(NA, 5000)

betaEGGER1 <- rep(NA, 5000)
betaEGGER2 <- rep(NA, 5000)
betaEGGER3 <- rep(NA, 5000)
betaEGGER4 <- rep(NA, 5000)
betaEGGER5 <- rep(NA, 5000)
betaEGGER6 <- rep(NA, 5000)

sebetaEGGER1 <- rep(NA, 5000)
sebetaEGGER2 <- rep(NA, 5000)
sebetaEGGER3 <- rep(NA, 5000)
sebetaEGGER4 <- rep(NA, 5000)
sebetaEGGER5 <- rep(NA, 5000)
sebetaEGGER6 <- rep(NA, 5000)




genos <- read_plink("/gpfs/research/chongwu/shengjie/project3/simulation/usedSNP_PVRL2_19")
genevar <- genos$bed[, sample(ncol(genos$bed), 25)]
gamma <- runif(25, 0.03, 0.1)
alpha <- runif(25, -0.2, 0.2)
phi <- runif(25, -0.2, 0.2)

for (k in 1:5000){
print(k)
e_u <- rnorm(10000, 0, 1)
e_x <- rnorm(10000, 0, 1)
e_y <- rnorm(10000, 0, 1)

for (i in 1:10000){
 	indd <- which(genevar[i,] != "NA")
	ind_0.1 <- sample(indd, 3)
	ind_0.2 <- c(ind_0.1, sample(indd, 2))
	ind_0.3 <- c(ind_0.2, sample(indd, 3))
	U1[i] <- phi[ind_0.1] %*% genevar[i, ind_0.1] + e_u[i]
	X1[i] <- gamma[indd] %*% genevar[i,indd] + U1[i] + e_x[i]
	U2[i] <- phi[ind_0.2] %*% genevar[i, ind_0.2] + e_u[i]
	X2[i] <- gamma[indd] %*% genevar[i,indd] + U2[i] + e_x[i]
	U3[i] <- phi[ind_0.3] %*% genevar[i, ind_0.3] + e_u[i]
	X3[i] <- gamma[indd] %*% genevar[i,indd] + U3[i] + e_x[i]

	Y11[i] <- alpha[ind_0.1] %*% genevar[i,ind_0.1] + U1[i] + e_y[i]
	Y12[i] <- alpha[ind_0.2] %*% genevar[i,ind_0.2] + U2[i] + e_y[i]
	Y13[i] <- alpha[ind_0.3] %*% genevar[i,ind_0.3] + U3[i] + e_y[i]
	Y21[i] <- alpha[ind_0.1] %*% genevar[i,ind_0.1] + 0.1 * X1[i] + U1[i] + e_y[i]
	Y22[i] <- alpha[ind_0.2] %*% genevar[i,ind_0.2] + 0.1 * X2[i] + U2[i] + e_y[i]
	Y23[i] <- alpha[ind_0.3] %*% genevar[i,ind_0.3] + 0.1 * X3[i] + U3[i] + e_y[i]
}
	
	
for ( j in 1:25){
	fit1 <- lm(X1[1:5000] ~ genevar[1:5000,j])
	test1 <- summary(fit1)
	betaXG1[j] <- test1$coefficient[2,1]
	sebetaXG1[j] <- test1$coefficient[2,2]
	
	fit2 <- lm(Y11[5001:10000] ~ genevar[5001:10000,j])
	test2 <- summary(fit2)
	betaYG11[j] <- test2$coefficient[2,1]
	sebetaYG11[j] <- test2$coefficient[2,2]

	fit3 <- lm(Y21[5001:10000] ~ genevar[5001:10000,j])
	test3 <- summary(fit3)
	betaYG12[j] <- test3$coefficient[2,1]
	sebetaYG12[j] <- test3$coefficient[2,2]

	fit4 <- lm(X2[1:5000] ~ genevar[1:5000,j])
	test4 <- summary(fit4)
	betaXG2[j] <- test4$coefficient[2,1]
	sebetaXG2[j] <- test4$coefficient[2,2]

	fit5 <- lm(Y12[5001:10000] ~ genevar[5001:10000,j])
	test5 <- summary(fit5)
	betaYG21[j] <- test5$coefficient[2,1]
	sebetaYG21[j] <- test5$coefficient[2,2]

	fit6 <- lm(Y22[5001:10000] ~ genevar[5001:10000,j])
	test6 <- summary(fit6)
	betaYG22[j] <- test6$coefficient[2,1]
	sebetaYG22[j] <- test6$coefficient[2,2]

	fit7 <- lm(X3[1:5000] ~ genevar[1:5000,j])
	test7 <- summary(fit7)
	betaXG3[j] <- test7$coefficient[2,1]
	sebetaXG3[j] <- test7$coefficient[2,2]
	
	fit8 <- lm(Y13[5001:10000] ~ genevar[5001:10000,j])
	test8 <- summary(fit8)
	betaYG31[j] <- test8$coefficient[2,1]
	sebetaYG31[j] <- test8$coefficient[2,2]
	
	fit9 <- lm(Y23[5001:10000] ~ genevar[5001:10000,j])
	test9 <- summary(fit9)
	betaYG32[j] <- test9$coefficient[2,1]
	sebetaYG32[j] <- test9$coefficient[2,2]
}

 betaIV1 = betaYG11/betaXG1 # ratio estimates
 weights1 = (sebetaYG11/betaXG1)^-2 # inverse-variance weights
 betaIVW1[k] = sum(betaYG11*betaXG1*sebetaYG11^-2)/sum(betaXG1^2*sebetaYG11^-2)


 betaIV2 = betaYG12/betaXG1 # ratio estimates
 weights2 = (sebetaYG12/betaXG1)^-2 # inverse-variance weights
 betaIVW2[k] = sum(betaYG12*betaXG1*sebetaYG12^-2)/sum(betaXG1^2*sebetaYG12^-2)


 betaIV3 = betaYG21/betaXG2 # ratio estimates
 weights3 = (sebetaYG21/betaXG2)^-2 # inverse-variance weights
 betaIVW3[k] = sum(betaYG21*betaXG2*sebetaYG21^-2)/sum(betaXG2^2*sebetaYG21^-2)


 betaIV4 = betaYG22/betaXG2 # ratio estimates
 weights4 = (sebetaYG22/betaXG2)^-2 # inverse-variance weights
 betaIVW4[k] = sum(betaYG22*betaXG2*sebetaYG22^-2)/sum(betaXG2^2*sebetaYG22^-2)


 betaIV5 = betaYG31/betaXG3 # ratio estimates
 weights5 = (sebetaYG31/betaXG3)^-2 # inverse-variance weights
 betaIVW5[k] = sum(betaYG31*betaXG3*sebetaYG31^-2)/sum(betaXG3^2*sebetaYG31^-2)


 betaIV6 = betaYG32/betaXG3 # ratio estimates
 weights6 = (sebetaYG32/betaXG3)^-2 # inverse-variance weights
 betaIVW6[k] = sum(betaYG32*betaXG3*sebetaYG32^-2)/sum(betaXG3^2*sebetaYG32^-2)

weighted.median <- function(betaIV.in, weights.in) {
betaIV.order = betaIV.in[order(betaIV.in)]
weights.order = weights.in[order(betaIV.in)]
weights.sum = cumsum(weights.order)-0.5*weights.order
weights.sum = weights.sum/sum(weights.order)
below = max(which(weights.sum<0.5))
weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
(0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
return(weighted.est) }


 beta_median1[k] <- weighted.median(betaIV1, weights1)
 beta_median2[k] <- weighted.median(betaIV2, weights2)
 beta_median3[k] <- weighted.median(betaIV3, weights3)
 beta_median4[k] <- weighted.median(betaIV4, weights4)
 beta_median5[k] <- weighted.median(betaIV5, weights5)
 beta_median6[k] <- weighted.median(betaIV6, weights6)


 weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in){
 med = NULL
 for(i in 1:1000){
 betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
 betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
 betaIV.boot = betaYG.boot/betaXG.boot
 med[i] = weighted.median(betaIV.boot, weights.in)
 }
 return(sd(med)) }

 sebeta_median1 [k] <- weighted.median.boot(betaXG1, betaYG11, sebetaXG1, sebetaYG11, weights1)
 sebeta_median2 [k] <- weighted.median.boot(betaXG1, betaYG12, sebetaXG1, sebetaYG12, weights2)
 sebeta_median3 [k] <- weighted.median.boot(betaXG2, betaYG21, sebetaXG2, sebetaYG21, weights3)
 sebeta_median4 [k] <- weighted.median.boot(betaXG2, betaYG22, sebetaXG2, sebetaYG22, weights4)
 sebeta_median5 [k] <- weighted.median.boot(betaXG3, betaYG31, sebetaXG3, sebetaYG31, weights5)
 sebeta_median6 [k] <- weighted.median.boot(betaXG3, betaYG32, sebetaXG3, sebetaYG32, weights6)

 # IVW estimate
 penalty1 = pchisq(weights1*(betaIV1-betaIVW1[k])^2, df=1, lower.tail=FALSE)
 penalty2 = pchisq(weights2*(betaIV2-betaIVW2[k])^2, df=1, lower.tail=FALSE)
 penalty3 = pchisq(weights3*(betaIV3-betaIVW3[k])^2, df=1, lower.tail=FALSE)
 penalty4 = pchisq(weights4*(betaIV4-betaIVW4[k])^2, df=1, lower.tail=FALSE)
 penalty5 = pchisq(weights5*(betaIV5-betaIVW5[k])^2, df=1, lower.tail=FALSE)
 penalty6 = pchisq(weights6*(betaIV6-betaIVW6[k])^2, df=1, lower.tail=FALSE)

 pen.weights1 = weights1*pmin(1, penalty1*20) # penalized weights
 pen.weights2 = weights2*pmin(1, penalty2*20) # penalized weights
 pen.weights3 = weights3*pmin(1, penalty3*20) # penalized weights
 pen.weights4 = weights4*pmin(1, penalty4*20) # penalized weights
 pen.weights5 = weights5*pmin(1, penalty5*20) # penalized weights
 pen.weights6 = weights6*pmin(1, penalty6*20) # penalized weights

 betaWM1[k] = weighted.median(betaIV1, weights1) # weighted median estimate
 betaWM2[k] = weighted.median(betaIV2, weights2) # weighted median estimate
 betaWM3[k] = weighted.median(betaIV3, weights3) # weighted median estimate
 betaWM4[k] = weighted.median(betaIV4, weights4) # weighted median estimate
 betaWM5[k] = weighted.median(betaIV5, weights5) # weighted median estimate
 betaWM6[k] = weighted.median(betaIV6, weights6) # weighted median estimate

 sebetaWM1[k] = weighted.median.boot(betaXG1, betaYG11, sebetaXG1, sebetaYG11, weights1)
 sebetaWM2[k] = weighted.median.boot(betaXG1, betaYG12, sebetaXG1, sebetaYG12, weights2)
 sebetaWM3[k] = weighted.median.boot(betaXG2, betaYG21, sebetaXG2, sebetaYG21, weights3)
 sebetaWM4[k] = weighted.median.boot(betaXG2, betaYG22, sebetaXG2, sebetaYG22, weights4)
 sebetaWM5[k] = weighted.median.boot(betaXG3, betaYG31, sebetaXG3, sebetaYG31, weights5)
 sebetaWM6[k] = weighted.median.boot(betaXG3, betaYG32, sebetaXG3, sebetaYG32, weights6)

 # standard error
 betaPWM1[k] = weighted.median(betaIV1, pen.weights1) # penalized weighted median estimate
 betaPWM2[k] = weighted.median(betaIV2, pen.weights2) # penalized weighted median estimate
 betaPWM3[k] = weighted.median(betaIV3, pen.weights3) # penalized weighted median estimate
 betaPWM4[k] = weighted.median(betaIV4, pen.weights4) # penalized weighted median estimate
 betaPWM5[k] = weighted.median(betaIV5, pen.weights5) # penalized weighted median estimate
 betaPWM6[k] = weighted.median(betaIV6, pen.weights6) # penalized weighted median estimate

 sebetaPWM1[k] = weighted.median.boot(betaXG1, betaYG11, sebetaXG1, sebetaYG11, pen.weights1)
 sebetaPWM2[k] = weighted.median.boot(betaXG1, betaYG12, sebetaXG1, sebetaYG12, pen.weights2)
 sebetaPWM3[k] = weighted.median.boot(betaXG2, betaYG21, sebetaXG2, sebetaYG21, pen.weights3)
 sebetaPWM4[k] = weighted.median.boot(betaXG2, betaYG22, sebetaXG2, sebetaYG22, pen.weights4)
 sebetaPWM5[k] = weighted.median.boot(betaXG3, betaYG31, sebetaXG3, sebetaYG31, pen.weights5)
 sebetaPWM6[k] = weighted.median.boot(betaXG3, betaYG32, sebetaXG3, sebetaYG32, pen.weights6)

 betaIVW1[k] = summary(lm(betaYG11~betaXG1-1, weights=sebetaYG11^-2))$coef[1,1]
 betaIVW2[k] = summary(lm(betaYG12~betaXG1-1, weights=sebetaYG12^-2))$coef[1,1]
 betaIVW3[k] = summary(lm(betaYG21~betaXG2-1, weights=sebetaYG21^-2))$coef[1,1]
 betaIVW4[k] = summary(lm(betaYG22~betaXG2-1, weights=sebetaYG22^-2))$coef[1,1]
 betaIVW5[k] = summary(lm(betaYG31~betaXG3-1, weights=sebetaYG31^-2))$coef[1,1]
 betaIVW6[k] = summary(lm(betaYG32~betaXG3-1, weights=sebetaYG32^-2))$coef[1,1]

 sebetaIVW1[k] = summary(lm(betaYG11~betaXG1-1, weights=sebetaYG11^-2))$coef[1,2]/
 min(summary(lm(betaYG11~betaXG1-1, weights=sebetaYG11^-2))$sigma, 1)
 sebetaIVW2[k] = summary(lm(betaYG12~betaXG1-1, weights=sebetaYG12^-2))$coef[1,2]/
 min(summary(lm(betaYG12~betaXG1-1, weights=sebetaYG12^-2))$sigma, 1)
 sebetaIVW3[k] = summary(lm(betaYG21~betaXG2-1, weights=sebetaYG21^-2))$coef[1,2]/
 min(summary(lm(betaYG21~betaXG2-1, weights=sebetaYG21^-2))$sigma, 1)
 sebetaIVW4[k] = summary(lm(betaYG22~betaXG2-1, weights=sebetaYG22^-2))$coef[1,2]/
 min(summary(lm(betaYG22~betaXG2-1, weights=sebetaYG22^-2))$sigma, 1)
 sebetaIVW5[k] = summary(lm(betaYG31~betaXG3-1, weights=sebetaYG31^-2))$coef[1,2]/
 min(summary(lm(betaYG31~betaXG3-1, weights=sebetaYG31^-2))$sigma, 1)
 sebetaIVW6[k] = summary(lm(betaYG32~betaXG3-1, weights=sebetaYG32^-2))$coef[1,2]/
 min(summary(lm(betaYG32~betaXG3-1, weights=sebetaYG32^-2))$sigma, 1)

 betaEGGER1[k] = summary(lm(betaYG11~betaXG1, weights=sebetaYG11^-2))$coef[2,1]
 betaEGGER2[k] = summary(lm(betaYG12~betaXG1, weights=sebetaYG12^-2))$coef[2,1]
 betaEGGER3[k] = summary(lm(betaYG21~betaXG2, weights=sebetaYG21^-2))$coef[2,1]
 betaEGGER4[k] = summary(lm(betaYG22~betaXG2, weights=sebetaYG22^-2))$coef[2,1]
 betaEGGER5[k] = summary(lm(betaYG31~betaXG3, weights=sebetaYG31^-2))$coef[2,1]
 betaEGGER6[k] = summary(lm(betaYG32~betaXG3, weights=sebetaYG32^-2))$coef[2,1]


 sebetaEGGER1[k] = summary(lm(betaYG11~betaXG1, weights=sebetaYG11^-2))$coef[2,2]
 min(summary(lm(betaYG11~betaXG1, weights=sebetaYG11^-2))$sigma, 1)
 sebetaEGGER2[k] = summary(lm(betaYG12~betaXG1, weights=sebetaYG21^-2))$coef[2,2]/
 min(summary(lm(betaYG12~betaXG1, weights=sebetaYG12^-2))$sigma, 1)
 sebetaEGGER3[k] = summary(lm(betaYG21~betaXG2, weights=sebetaYG21^-2))$coef[2,2]/
 min(summary(lm(betaYG21~betaXG2, weights=sebetaYG21^-2))$sigma, 1)
 sebetaEGGER4[k] = summary(lm(betaYG22~betaXG2, weights=sebetaYG22^-2))$coef[2,2]/
 min(summary(lm(betaYG22~betaXG2, weights=sebetaYG22^-2))$sigma, 1)
 sebetaEGGER5[k] = summary(lm(betaYG31~betaXG3, weights=sebetaYG31^-2))$coef[2,2]/
 min(summary(lm(betaYG31~betaXG3, weights=sebetaYG31^-2))$sigma, 1)
 sebetaEGGER6[k] = summary(lm(betaYG32~betaXG3, weights=sebetaYG32^-2))$coef[2,2]/
 min(summary(lm(betaYG32~betaXG3, weights=sebetaYG32^-2))$sigma, 1)
 }

saveRDS(betaWM1, "betaWM_null_S3_0.1.rds")
saveRDS(betaWM2, "betaWM_null_S3_0.2.rds")
saveRDS(betaWM3, "betaWM_null_S3_0.3.rds")
saveRDS(betaWM4, "betaWM_positive_S3_0.1.rds")
saveRDS(betaWM5, "betaWM_positive_S3_0.2.rds")
saveRDS(betaWM6, "betaWM_positive_S3_0.3.rds")
 print(mean(betaWM1))
 print(mean(betaWM2))
 print(mean(betaWM3))
 print(mean(betaWM4))
 print(mean(betaWM5))
 print(mean(betaWM6))



saveRDS(sebetaWM1, "sebetaWM_null_S3_0.1.rds")
saveRDS(sebetaWM2, "sebetaWM_null_S3_0.2.rds")
saveRDS(sebetaWM3, "sebetaWM_null_S3_0.3.rds")
saveRDS(sebetaWM4, "sebetaWM_positive_S3_0.1.rds")
saveRDS(sebetaWM5, "sebetaWM_positive_S3_0.2.rds")
saveRDS(sebetaWM6, "sebetaWM_positive_S3_0.3.rds")
 print(mean(sebetaWM1))
 print(mean(sebetaWM2))
 print(mean(sebetaWM3))
 print(mean(sebetaWM4))
 print(mean(sebetaWM5))
 print(mean(sebetaWM6))



saveRDS(betaPWM1, "betaPWM_null_S3_0.1.rds")
saveRDS(betaPWM2, "betaPWM_null_S3_0.2.rds")
saveRDS(betaPWM3, "betaPWM_null_S3_0.3.rds")
saveRDS(betaPWM4, "betaPWM_positive_S3_0.1.rds")
saveRDS(betaPWM5, "betaPWM_positive_S3_0.2.rds")
saveRDS(betaPWM6, "betaPWM_positive_S3_0.3.rds")
 print(mean(betaPWM1))
 print(mean(betaPWM2))
 print(mean(betaPWM3))
 print(mean(betaPWM4))
 print(mean(betaPWM5))
 print(mean(betaPWM6))


saveRDS(sebetaPWM1, "sebetaPWM_null_S3_0.1.rds")
saveRDS(sebetaPWM2, "sebetaPWM_null_S3_0.2.rds")
saveRDS(sebetaPWM3, "sebetaPWM_null_S3_0.3.rds")
saveRDS(sebetaPWM4, "sebetaPWM_positive_S3_0.1.rds")
saveRDS(sebetaPWM5, "sebetaPWM_positive_S3_0.2.rds")
saveRDS(sebetaPWM6, "sebetaPWM_positive_S3_0.3.rds")
 print(mean(sebetaPWM1))
 print(mean(sebetaPWM2))
 print(mean(sebetaPWM3))
 print(mean(sebetaPWM4))
 print(mean(sebetaPWM5))
 print(mean(sebetaPWM6))



saveRDS(betaIVW1, "betaIVW_null_S3_0.1.rds")
saveRDS(betaIVW2, "betaIVW_null_S3_0.2.rds")
saveRDS(betaIVW3, "betaIVW_null_S3_0.3.rds")
saveRDS(betaIVW4, "betaIVW_positive_S3_0.1.rds")
saveRDS(betaIVW5, "betaIVW_positive_S3_0.2.rds")
saveRDS(betaIVW6, "betaIVW_positive_S3_0.3.rds")
print(mean(betaIVW1))
print(mean(betaIVW2))
print(mean(betaIVW3))
print(mean(betaIVW4))
print(mean(betaIVW5))
print(mean(betaIVW6))


saveRDS(sebetaIVW1, "sebetaIVW_null_S3_0.1.rds")
saveRDS(sebetaIVW2, "sebetaIVW_null_S3_0.2.rds") 
saveRDS(sebetaIVW3, "sebetaIVW_null_S3_0.3.rds") 
saveRDS(sebetaIVW4, "sebetaIVW_positive_S3_0.1.rds") 
saveRDS(sebetaIVW5, "sebetaIVW_positive_S3_0.2.rds") 
saveRDS(sebetaIVW6, "sebetaIVW_positive_S3_0.3.rds") 
print(mean(sebetaIVW1))
print(mean(sebetaIVW2))
print(mean(sebetaIVW3))
print(mean(sebetaIVW4)) 
print(mean(sebetaIVW5))
print(mean(sebetaIVW6))

saveRDS(betaEGGER1, "betaEGGER_null_S3_0.1.rds")
saveRDS(betaEGGER2, "betaEGGER_null_S3_0.2.rds")
saveRDS(betaEGGER3, "betaEGGER_null_S3_0.3.rds")
saveRDS(betaEGGER4, "betaEGGER_positive_S3_0.1.rds")
saveRDS(betaEGGER5, "betaEGGER_positive_S3_0.2.rds")
saveRDS(betaEGGER6, "betaEGGER_positive_S3_0.3.rds")
print(mean(betaEGGER1))
print(mean(betaEGGER2))
print(mean(betaEGGER3))
print(mean(betaEGGER4))
print(mean(betaEGGER5))
print(mean(betaEGGER6))

saveRDS(sebetaEGGER1, "sebetaEGGER_null_S3_0.1.rds")
saveRDS(sebetaEGGER2, "sebetaEGGER_null_S3_0.2.rds")
saveRDS(sebetaEGGER3, "sebetaEGGER_null_S3_0.3.rds")
saveRDS(sebetaEGGER4, "sebetaEGGER_positive_S3_0.1.rds")
saveRDS(sebetaEGGER5, "sebetaEGGER_positive_S3_0.2.rds")
saveRDS(sebetaEGGER6, "sebetaEGGER_positive_S3_0.2.rds")
print(mean(sebetaEGGER1))
print(mean(sebetaEGGER2))
print(mean(sebetaEGGER3))
print(mean(sebetaEGGER4))
print(mean(sebetaEGGER5))
print(mean(sebetaEGGER6))

