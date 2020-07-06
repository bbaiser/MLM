#FULL VERSION

library(vegan)
install.packages("lme4")
library(lme4)

nreps <- 20

ddddd

# binary or count data
binaryflag <- T

# number of sites, species, and environmental variables
nsite <- 30;
nspp <- 20
nx <- 5

# construct environmental variables with covariance matrix V
Vx <- rbind(c(1,0,0,0,0),c(0,1,.7,0,0),c(0,.7,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1))
v <- svd(Vx)
C <- v$u %*% diag(v$d^(.5))

MLM.ranef <- NULL
MLM.fixef <- NULL
MLM.pvals <- NULL
MLM.fixef.pvals <- NULL
#
# CCA.pvals <- NULL
# RDA.pvals <- NULL
# NMDS.pvals <- NULL

for(rep in 1:nreps){
	show(rep)

	# simulate data
	X <- array(0,c(nx,nsite))
	for(j in 1:nsite){
		X[,j] <- C %*% rnorm(nx)
	}
	X <- t(X)

	for(j in 1:nx) X[,j] <- (X[,j] - mean(X[,j]))/sd(X[,j])
	cov(X)

	# b0 and b1 for each species-environmental variable combination
	b0bar <- 0
	b1bar <- c(0,0,0,.5,0)

	b0sig <- 1
	b1sig <- c(.4,.4,0,0,0)

	b0 <- rnorm(nspp,mean=b0bar,sd=b0sig)

	b1 <- array(0,c(nx,nspp))
	for(j in 1:nx){
		b1[j,] <- rnorm(nspp,mean=b1bar[j],sd=b1sig[j])
	}

	# construct data set

	# lam is the parameter from a Poisson distribution giving the probability of occurrence
	lam <- array(1,c(nsite,1)) %*% b0 + X %*% b1

	Y <- array(0,c(nsite,nspp))
	for(j in 1:nsite){
		Y[j,] <- rpois(nspp,exp(lam[j,]))
	}

	# plots Y
	#if(binaryflag == T)	image(Y>0)
	#if(binaryflag == F)	image(Y)

	# construct data.frame
	dim(Y) <- c(nspp*nsite,1)

	XX <- NULL
	for(j in 1:nx){
		XX <- cbind(XX,rep(X[,j],nspp))
	}

	sitenames <- array(1:nsite,c(nspp*nsite,1))
	sppnames <- kronecker(1:nspp,array(1,c(nsite,1)))

	d <- data.frame(site = as.factor(sitenames), spp = as.factor(sppnames), Y = Y, X = XX)

	# analyze data: lmer
	if(binaryflag == T){
		d$YY <- (d$Y > 0)
		z <- glmer(YY ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=binomial)
		# compute ranef pvalues
		z1 <- glmer(YY ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=binomial)
		z2 <- glmer(YY ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=binomial)
		z3 <- glmer(YY ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=binomial)
		z4 <- glmer(YY ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.5 | spp), data=d, family=binomial)
		z5 <- glmer(YY ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp), data=d, family=binomial)
	}

	if(binaryflag == F){
		z <- glmer(Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=poisson)
		# compute ranef pvalues
		z1 <- glmer(Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=poisson)
		z2 <- glmer(Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=poisson)
		z3 <- glmer(Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.4 | spp) + (0 + X.5 | spp), data=d, family=poisson)
		z4 <- glmer(Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.5 | spp), data=d, family=poisson)
		z5 <- glmer(Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + (1 | spp) + (0 + X.1 | spp) + (0 + X.2 | spp) + (0 + X.3 | spp) + (0 + X.4 | spp), data=d, family=poisson)
	}

	LL <- c(logLik(z),logLik(z1),logLik(z2),logLik(z3),logLik(z4),logLik(z5))

	mlm.pvals <- c(1-pchisq(LL[2] - LL[1],1), 1-pchisq(LL[3] - LL[1],1), 1-pchisq(LL[4] - LL[1],1), 1-pchisq(LL[5] - LL[1],1), 1-pchisq(LL[6] - LL[1],1))

	# accumulate results
	MLM.pvals <- rbind(MLM.pvals,mlm.pvals)
	MLM.ranef <- rbind(MLM.ranef, ranef(z))
	MLM.fixef <- rbind(MLM.fixef, fixef(z))
	MLM.fixef.pvals <- rbind(MLM.fixef.pvals, summary(z)$coef[,4])
}
################################################
# output results

MLM.ranef
MLM.fixef
MLM.pvals
# CCA.pvals
# NMDS.pvals


barplot(colMeans(MLM.pvals < 0.05),main="MLM", xlab="environ variable", ylab="prop", ylim=c(0,1))

# barplot(colMeans(CCA.pvals < 0.05),main="CCA", xlab="environ variable", ylab="prop", ylim=c(0,1))
#
# barplot(colMeans(NMDS.pvals < 0.05),main="NMDS", xlab="environ variable", ylab="prop", ylim=c(0,1))

## Exercise: add CCA and RDA into simulation for loop.
