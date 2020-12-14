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

U <- rep(NA, 10000)
X <- rep(NA, 10000)
Y11 <- rep(NA, 10000)
Y12 <- rep(NA, 10000)
Y13 <- rep(NA, 10000)
Y21 <- rep(NA, 10000)
Y22 <- rep(NA, 10000)
Y23 <- rep(NA, 10000)

betaXG <- rep(NA, 25)
sebetaXG <- rep(NA, 25)
betaYG <- rep(NA, 25)
sebetaYG <- rep(NA, 25)
beta_median <- rep(NA, 5000)
sebeta_median <- rep(NA, 5000)
betaWM <- rep(NA, 5000)
sebetaWM <- rep(NA, 5000)
betaPWM <- rep(NA, 5000)
sebetaPWM <- rep(NA, 5000)
betaIVW <- rep(NA, 5000)
sebetaIVW <- rep(NA, 5000)
betaEGGER <- rep(NA, 5000)
sebetaEGGER <- rep(NA, 5000)



genos <- read_plink("/gpfs/research/chongwu/shengjie/project3/simulation/usedSNP_PVRL2_19")
genevar <- genos$bed[, sample(ncol(genos$bed), 25)]


#Scenario 1 (phi=0 and alpha~Unifrom(-0.2,0.2) with null causal effect)

for (k in 1:5000){
	print(k)

	e_u <- rnorm(10000, 0, 1)
	e_x <- rnorm(10000, 0, 1)
	e_y <- rnorm(10000, 0, 1)


	gamma <- runif(25, 0.03, 0.1)
	alpha <- runif(25, -0.2, 0.2)


	for (i in 1:10000){
		
		U[i] <- e_u[i]
		indd <- which(genevar[i,] != "NA")
		X[i] <- gamma[indd] %*% genevar[i,indd] + U[i] + e_x[i]
		ind_0.1 <- sample(indd, 3) 
		ind_0.2 <- sample(indd, 5)
		ind_0.3 <- sample(indd, 8)
		Y1[i] <- alpha[ind_0.2] %*% genevar[i,ind_0.2] + U[i] + e_y[i]

	}

	for ( j in 1:25){
		
		fit1 <- lm(X[1:5000] ~ genevar[1:5000,j])
		test1 <- summary(fit1)
		betaXG[j] <- test1$coefficient[2,1]
		sebetaXG[j] <- test1$coefficient[2,2]

		fit2 <- lm(Y[5001:10000] ~ genevar[5001:10000,j])
		test2 <- summary(fit2)
		betaYG[j] <- test2$coefficient[2,1]
		sebetaYG[j] <- test2$coefficient[2,2]
	}

	#betaXG
	#sebetaXG
	#betaYG
	#sebetaYG


	betaIV = betaYG/betaXG # ratio estimates
	weights = (sebetaYG/betaXG)^-2 # inverse-variance weights
	betaIVW[k] = sum(betaYG*betaXG*sebetaYG^-2)/sum(betaXG^2*sebetaYG^-2)

	weighted.median <- function(betaIV.in, weights.in) {
		betaIV.order = betaIV.in[order(betaIV.in)]
		weights.order = weights.in[order(betaIV.in)]
		weights.sum = cumsum(weights.order)-0.5*weights.order
		weights.sum = weights.sum/sum(weights.order)
		below = max(which(weights.sum<0.5))
		weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
		(0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
		return(weighted.est) }

	beta_median[k] <- weighted.median(betaIV, weights)


	weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in){
		med = NULL
		for(i in 1:1000){
		betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
		betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
		betaIV.boot = betaYG.boot/betaXG.boot
		med[i] = weighted.median(betaIV.boot, weights.in)
		}
		return(sd(med)) }

	sebeta_median [k] <- weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, weights)


	# IVW estimate
	penalty = pchisq(weights*(betaIV-betaIVW[k])^2, df=1, lower.tail=FALSE)
	pen.weights = weights*pmin(1, penalty*20) # penalized weights
	betaWM[k] = weighted.median(betaIV, weights) # weighted median estimate
	sebetaWM[k] = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, weights)

	# standard error
	betaPWM[k] = weighted.median(betaIV, pen.weights) # penalized weighted median estimate
	sebetaPWM[k] = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, pen.weights)


	betaIVW[k] = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,1]
	sebetaIVW[k] = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,2]/
	min(summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$sigma, 1)
	betaEGGER[k] = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$coef[2,1]
	sebetaEGGER[k] = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$coef[2,2]/
	min(summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$sigma, 1)

}

saveRDS(betaWM, "betaWM_null_S1_0.1.rds")
saveRDS(sebetaWM, "sebetaWM_null_S1_0.1.rds")
saveRDS(betaPWM, "betaPWM_null_S1_0.1.rds")
saveRDS(sebetaPWM, "sebetaPWM_null_S1_0.1.rds")
saveRDS(betaIVW, "betaIVW_null_S1_0.1.rds")
saveRDS(sebetaIVW, "sebetaIVW_null_S1_0.1.rds") 
saveRDS(betaEGGER, "betaEGGER_null_S1_0.1.rds")
saveRDS(sebetaEGGER, "sebetaEGGER_null_S1_0.1.rds")

