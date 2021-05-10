#### GLMER for LPK Presence/Absence Data ####
# 1) Read in Data
# 2) Brief data exploration (see DataExploration.R for more)
# 3) Final Models for both time periods
# 4) Univariate models (LONG)

library(tidyverse)
library(lme4)
library(ggplot2)
library(ggfortify)
library(dfoptim)
library(optimx)
library(lattice)




#################### 1 ####################

#Read in data for binary comparison
#New vegetation
veg<- read_csv("data/1m_New_Wide.csv")
veg <- veg[-c(41:88, 137:184),] #remove sites without explanatory variables
veg <- veg[,-1]
veg <- veg[, colSums(veg != 0) > 0] # remove species that are not present

#Slocum vegetation
veg.slocum <- read_csv("data/1m_Old_Wide.csv")
veg.slocum <- veg.slocum[-c(41:88, 137:184),] #remove sites without variables
veg.slocum <- veg.slocum[,-1]
veg.slocum <- veg.slocum[, colSums(veg.slocum != 0) > 0] # remove species that are not present

#Variable sheet
var <- read_csv("data/Variables.csv")

#combine variable with each time series of veg data.
veg <- cbind(var, veg)
veg.slocum<- cbind(var, veg.slocum)

#Remove unnecessary variables
veg <- veg[,c(1:3, 7, 18:214)]
veg.slocum <- veg.slocum[,c(1:3, 7, 11:17, 25:217)]

#Save Names for later
veg.names <- veg[,c(12:201)]
veg.slocum.names <- veg.slocum[,c(12:204)]

###########################################################################
#Normalize predictor variables (could also just use scale())
for (i in 4:11) {
  veg[, i] <- (veg[, i] - mean(veg[, i])) / sd(veg[, i])
}

#Long format
veg<- gather(veg, key = "SPP", value = "Presence", ACALCHAM:ZAMIINTE)
veg$Presence[veg$Presence>0] <- 1

#Need to make specific data names for sites
veg$Submodule <- c("A", "B", "C", "D", "E", "F", "G", "H")
Submodule2 <- unite(veg, "Submodule2", c("Block", "Plot", "Submodule"), sep = "_", remove=TRUE)
Submodule2 <- Submodule2[,1]
Plot2 <- unite(veg, "Plot2", c("Block", "Plot"), sep = "_", remove = TRUE)
Plot2 <- Plot2[,1]
veg <- cbind(Submodule2, Plot2, veg)

###########################################################################
#Normalize predictor variables (could also just use scale())
for (i in 4:11) {
  veg.slocum[, i] <- (veg.slocum[, i] - mean(veg.slocum[, i])) / sd(veg.slocum[, i])
}

#Long format
veg.slocum<- gather(veg.slocum, key = "SPP", value = "Presence", ACALCHAM:ZAMIINTE)
veg.slocum$Presence[veg.slocum$Presence>0] <- 1

#Need to make specific data names for sites
veg.slocum$Submodule <- c("A", "B", "C", "D", "E", "F", "G", "H")
Submodule2 <- unite(veg.slocum, "Submodule2", c("Block", "Plot", "Submodule"), sep = "_", remove=TRUE)
Submodule2 <- Submodule2[,1]
Plot2 <- unite(veg.slocum, "Plot2", c("Block", "Plot"), sep = "_", remove = TRUE)
Plot2 <- Plot2[,1]
veg.slocum <- cbind(Submodule2, Plot2, veg.slocum)






#################### 2 ####################

# Data Exploration (Basic, see Data_Exploration.R for more)
table(veg$Presence) #number of pres/abs
table(veg.slocum$Presence) #number of pres/abs

#Pairs plot new (these take a long time to load)
pairs(veg[,6:13], pch=15, cex=0.8,
      lower.panel = NULL)

pairs(veg.slocum[,6:13], pch=15, cex=0.8,
      lower.panel = NULL)

#### VIF ####
# VIF < 5         Good
# 10 >= VIF >= 5  Moderate
# VIF > 10        Bad

source("code/HighstatLibV10.r")
veg.x <- veg[,6:13]
X <- corvif(veg.x)#In cov2cor(v) : diag(.) had 0 or NA entries; non-finite result is doubtful

veg.x <- veg.x[,-3] #remove Water Surface...don't need it. 
X <- corvif(veg.x)

veg.x <- veg.x[,-5] #Remove Fire Season
X <- corvif(veg.x)

veg.x <- veg.x[,-3] #Remove Water Depth
X <- corvif(veg.x)

veg.x <- veg.x[,-4] #Remove Wet Season
X <- corvif(veg.x) 

veg.x <- veg.x[,-3] #Remove Dry Season
X <- corvif(veg.x) #GOOD!





#################### 3 ####################
# Present Data
# Remove RE systematically to determine mlm pvalues
z1 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              IN_20 + (0 + IN_20|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial",
            control=glmerControl(optimizer="bobyqa"),
            data = veg)
summary(z1)

ranef(z1)$SPP$Elevation
fixef(z1)

#Without Burn Fequency
z2 <- glmer(Presence ~ (1|SPP) + BF_2000 +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(IN_20, degree = 2) + (0 + IN_20|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(z2)

#Without Elevation ###this originally contained elevation
z3 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
              poly(IN_20, degree = 2) + (0 + IN_20|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(z3)
logit.predictions<-predict(z3)

prob.predictions <- 1 / (1 + exp(-logit.predictions))

ff<-round(prob.predictions,2)
hist(ff)
ranef(z1)$SPP$Elevation
fixef(z1)


install.packages("detectseperation")
library(brglm2)

sep <- glm(Presence ~ BF_2000 +poly(Elevation, degree = 2) + IN_20 , data = veg,
                       family = binomial("logit"),
                       method = "detectseparation")
summary(sep)
diag(9,4)

Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
  poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
  IN_20 + (0 + IN_20|SPP) + 
  (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
family = "binomial",

#Without Inundation
z4 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(IN_20, degree = 2) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(z4)

LL <- c(deviance(z1), deviance(z2), deviance(z3), deviance(z4))

x <- 1 - pchisq(LL[2] - LL[1], 1)
#these are chi-square distribution functions). 
mlm.pvals <- c(
  1 - pchisq(LL[2] - LL[1], 1),
  1 - pchisq(LL[3] - LL[1], 1),
  1 - pchisq(LL[4] - LL[1], 1))

mlm.pvals

# compute residual effect of random effects environmental variables
nsite <- 176
nspp <- 190

glmer.resid <- glmer(Presence ~ (1 | SPP) + BF_2000 +
                       poly(Elevation, degree = 2) + 
                       poly(IN_20, degree = 2),
                     data = veg,
                     family = "binomial"
)

MLM.fitted <- array(fitted(z1) - fitted(glmer.resid), c(nsite, nspp))

rownames(MLM.fitted) = c(1:176)
colnames(MLM.fitted) = names(veg.names)[1:190]

# standardize over spp
MLM.fitted.standard <- MLM.fitted
for (j in 1:nspp)
  MLM.fitted.standard[, j] <-
  (MLM.fitted[, j] - mean(MLM.fitted[, j])) /
  sd(MLM.fitted[, j])

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[, 1:2]

# environmental variables (only those with random effects)
envir.vars <- cbind(veg$Elevation, veg$BF_2000, veg$IN_20)
mlm.envir <- NULL
for (j in 1:3)
  mlm.envir <-
  cbind(mlm.envir, envir.vars[, j] * mlm.fit[, 1], envir.vars[, j] * mlm.fit[, 2])

envir.points <- matrix(colMeans(mlm.envir), nrow = ncol(envir.vars), ncol = 2, byrow = T)

# plot mlm
par(mfcol = c(1, 1))
plot(-mlm.fit,
     xlab = "PC1",
     ylab = "PC2",
     type = "n")
text(-mlm.fit, label = 1:nsite, cex = .5)

arrow.coordMLM <- cbind(array(0, dim(envir.points)), -envir.points)

arrows(
  arrow.coordMLM[, 1],
  arrow.coordMLM[, 2],
  arrow.coordMLM[, 3],
  arrow.coordMLM[, 4],
  col = "black",
  length = 0.1
)

text(
  1.3 * -envir.points,
  label = c("Elevation", "BF", "IN"),
  cex = .7
)

title(main = "PCA 2010-Present, BF, Elev, Inundation")

################ Slocum Data ####################
# Remove RE systematically to determine mlm pvalues
z1.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
                     poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
                     IN_20_Old + (0 + IN_20_Old|SPP) + 
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial", 
                   control=glmerControl(optimizer="bobyqa"),
                   data = veg.slocum)
summary(z1.o)

#Without Burn Fequency
z2.o <- glmer(Presence ~ (1|SPP) + BF_1980 +
                poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
                IN_20_Old + (0 + IN_20_Old|SPP) + 
                (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
              family = "binomial", 
              data = veg.slocum)
summary(z2.o)

#Without Elevation
z3.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
                poly(Elevation, degree = 2) + 
                IN_20_Old + (0 + IN_20_Old|SPP) + 
                (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
              family = "binomial", 
              data = veg.slocum)
summary(z3.o)

#Without Inundation
z4.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
                poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
                IN_20_Old + 
                (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
              family = "binomial", 
              data = veg.slocum)
summary(z4.o)

LL <- c(deviance(z1.o), deviance(z2.o), deviance(z3.o), deviance(z4.o))

x <- 1 - pchisq(LL[2] - LL[1], 1)
#these are chi-square distribution functions). 
mlm.pvals <- c(
  1 - pchisq(LL[2] - LL[1], 1),
  1 - pchisq(LL[3] - LL[1], 1),
  1 - pchisq(LL[4] - LL[1], 1))

mlm.pvals

# compute residual effect of random effects environmental variables
nsite <- 176
nspp <- 193
glmer.resid.o <- glmer(Presence ~ (1 | SPP) + BF_1980 + 
                         poly(Elevation, degree = 2) + 
                         poly(IN_20_Old, degree = 2),
                       data = veg.slocum,
                       family = "binomial"
)

MLM.fitted <- array(fitted(z1.o) - fitted(glmer.resid.o), c(nsite, nspp))

rownames(MLM.fitted) = c(1:176)
colnames(MLM.fitted) = names(veg.slocum.names)[1:193]

# standardize over spp
MLM.fitted.standard <- MLM.fitted
for (j in 1:nspp)
  MLM.fitted.standard[, j] <-
  (MLM.fitted[, j] - mean(MLM.fitted[, j])) /
  sd(MLM.fitted[, j])

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[, 1:2]

# environmental variables (only those with random effects)
envir.vars <- cbind(veg.slocum$Elevation, veg.slocum$BF_1980, veg.slocum$IN_20_Old)
mlm.envir <- NULL
for (j in 1:3)
  mlm.envir <-
  cbind(mlm.envir, envir.vars[, j] * mlm.fit[, 1], envir.vars[, j] * mlm.fit[, 2])

envir.points <- matrix(colMeans(mlm.envir), nrow = ncol(envir.vars), ncol = 2, byrow = T)

# plot mlm
par(mfcol = c(1, 1))
plot(-mlm.fit,
     xlab = "PC1",
     ylab = "PC2",
     type = "n")
text(-mlm.fit, label = 1:nsite, cex = .5)

arrow.coordMLM <- cbind(array(0, dim(envir.points)), -envir.points)

arrows(
  arrow.coordMLM[, 1],
  arrow.coordMLM[, 2],
  arrow.coordMLM[, 3],
  arrow.coordMLM[, 4],
  col = "black",
  length = 0.1
)

text(
  1.3 * -envir.points,
  label = c("Elevation", "BF", "IN"),
  cex = .7
)

title(main = "PCA, 1990-2000, BF, Elev, Inundation")







#######################################################
#### Run univariate models ####
#######################################################
#Elevation. AIC = 18518.4
glmer.elev <- glmer(Presence ~ (1|SPP) + Elevation,
                    family = "binomial", 
                    data = veg)
summary(glmer.elev)

#Elevation + Elevation2 AIC = 18484.1
glmer.elev2 <- glmer(Presence ~ (1|SPP) + poly(Elevation, degree = 2), 
                     family = "binomial",
                     data = veg)
summary(glmer.elev2)

#Elevation + Sites + SPP RE. AIC = 16551.7
glmer.elev3 <- glmer(Presence ~ (1|SPP) + Elevation +
                     (0 + Elevation|SPP) + 
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial",
                     data = veg)
summary(glmer.elev3)

#Elevation + Elevation2 + Sites + RE. AIC = 16354.3 **BEST**
glmer.elev4 <- glmer(Presence ~ (1|SPP) + poly(Elevation, degree = 2) +
                       (0 + Elevation|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial",
                     data = veg)
summary(glmer.elev4)
###########################################################################

#WD. AIC = 18522.8
glmer.wd <- glmer(Presence ~ (1|SPP) + WD,
                    family = "binomial", 
                    data = veg)
summary(glmer.wd)

#Elevation + Elevation2 AIC = 18496.1
glmer.wd2 <- glmer(Presence ~ (1|SPP) + poly(WD, degree = 2), 
                     family = "binomial",
                     data = veg)
summary(glmer.wd2)

#WD + Sites + RE SPP. AIC = 18380.9
glmer.WD3 <- glmer(Presence ~ (1|SPP) + WD + (0 + WD|SPP) +
                    (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial",
                   data = veg)
summary(glmer.WD3)

#WD + WD2 + Sites + RE. AIC = 16445.5
glmer.wd4 <- glmer(Presence ~ (1|SPP) + poly(WD, degree = 2) +
                    (0 + WD|SPP) +
                    (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial",
                     data = veg)
summary(glmer.wd4)

###########################################################################

#BF. AIC = 18503.3
glmer.bf <- glmer(Presence ~ (1|SPP) + BF_2000,
                  family = "binomial", 
                  data = veg)
summary(glmer.bf)

#BF + BF2 AIC = 18504.1
glmer.bf2 <- glmer(Presence ~ (1|SPP) + poly(BF_2000, degree = 2), 
                   family = "binomial",
                   data = veg)
summary(glmer.bf2)

#BF + Sites + SPP RE. AIC = 18389.2
glmer.bf3 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial",
                   data = veg)
summary(glmer.bf3)

#BF + BF2 + Sites + RE. AIC = 18222.7
glmer.bf4 <- glmer(Presence ~ (1|SPP) + poly(BF_2000, degree = 2) +
                    (0 + BF_2000|SPP) +
                    (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial",
                   data = veg)
summary(glmer.bf4)

###############################################################################
#Inundation AIC = 18530.2
glmer.inun <- glmer(Presence ~ (1|SPP) + 
                    IN_20,
                  family = "binomial", 
                  data = veg)
summary(glmer.inun)

#Inundation + Inundation2, AIC = 18530.2
glmer.inun.2 <- glmer(Presence ~ (1|SPP) + poly(IN_20, degree = 2),
                    family = "binomial", 
                    data = veg)
summary(glmer.inun.2)

#Inundation + SPP + Site, AIC = 17986.0
glmer.inun.3 <- glmer(Presence ~ (1|SPP) + IN_20 + 
                      (0 + IN_20|SPP) + 
                      (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                    family = "binomial", 
                    data = veg)
summary(glmer.inun.3)

#Inundation + SPP + Site, AIC = 17986.0
glmer.inun.4 <- glmer(Presence ~ (1|SPP) +  poly(IN_20, degree = 2) +
                      (0 + IN_20|SPP) + 
                      (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                    family = "binomial", 
                    data = veg)
summary(glmer.inun.4)

###########################################
#Dry Season Water WIndow, AIC = 18537.3
glmer.dry <- glmer(Presence ~ (1|SPP) + WDD,
                    family = "binomial", 
                    data = veg)
summary(glmer.dry)

#Dry Season Water Window + 2, AIC = 18537.3
glmer.dry.2 <- glmer(Presence ~ (1|SPP) + poly(WDD, degree = 2),
                     family = "binomial", 
                     data = veg)
summary(glmer.dry.2)

#Dry Season Water Window + 2, AIC = 18537.3
glmer.dry.3 <- glmer(Presence ~ (1|SPP) + WDDy +(0 + WDD|SPP) + 
                   (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial", 
                   data = veg)
summary(glmer.dry.3)

#Dry Season Water Window + Full RE, AIC = 18537.3
glmer.dry.4 <- glmer(Presence ~ (1|SPP) + poly(WDD, degree = 2) +
                     (0 + WDD|SPP) + 
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial", 
                     data = veg)
summary(glmer.dry.4)

#########################################################
#Fire Season Water Window, AIC = 18532.9
glmer.fire <- glmer(Presence ~ (1|SPP) + WDF,
                   family = "binomial", 
                   data = veg)
summary(glmer.fire)
tab_model(glmer.fire)

#Fire Season Water WIndow + Fire2, AIC = 18505.0
glmer.fire2 <- glmer(Presence ~ (1|SPP) + poly(WDF, degree = 2),
                    (0 + WDF|SPP),
                    family = "binomial", 
                    data = veg)
summary(glmer.fire2)
tab_model(glmer.fire)

#Water window Fire + Fire2 + Site, AIC = 18381.5
glmer.fire.3 <- glmer(Presence ~ (1|SPP) + (0 + WDF|SPP) + 
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial", 
                     data = veg)
summary(glmer.fire.3)

#Water window Fire + Fire2 + Site, AIC = 16424.1
glmer.fire.4 <- glmer(Presence ~ (1|SPP) + poly(Water_Depth_Fire, degree = 2) + 
                       (0 + WDF|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial", 
                     data = veg)
summary(glmer.fire.4)

##################################################################
#Wet Season Water WIndow, AIC = 18537.3
glmer.wet <- glmer(Presence ~ (1|SPP) + WDW,
                   family = "binomial", 
                   data = veg)
summary(glmer.wet)

#Wet Season Water Window + Wet2, AIC = 18537.3
glmer.wet.2 <- glmer(Presence ~ (1|SPP) + poly(WDW, degree = 2),
                   family = "binomial", 
                   data = veg)
summary(glmer.wet.2)

#Water window Wet + Wet2 + Site, AIC = 16578.7
glmer.wet.3 <- glmer(Presence ~ (1|SPP) + WDW + (0 + WDW|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial", 
                     data = veg)
summary(glmer.wet.3)

#Water window Wet + Wet2 + Site, AIC = 16578.7
glmer.wet.4 <- glmer(Presence ~ (1|SPP) + poly(WDW, degree = 2) + (0 + WDW|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial", 
                     data = veg)
summary(glmer.wet.4)


#Water Window all, AIC = 18525.9
glmer.water.win <- glmer(Presence ~ (1|SPP) + 
                     WDD + WDF + WDW,
                   family = "binomial", 
                   data = veg)
summary(glmer.water.win)
tab_model(glmer.dry)

#### Add Variables and Build ####
# Work through the variables 
m1 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
            (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(m1)

m2 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
            poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
            (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(m2)

m3 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(IN_20, degree = 2) + (0 + IN_20|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(m3)

m4 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(WD_New, degree = 2) + (0 + WD_New|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(m4)

m5 <- glmer(Presence ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(WDF, degree = 2) + (0 + WDF|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg)
summary(m5)

m6 <- glmer(ABUND ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
                poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
                poly(WDD, degree = 2) + (0 + WDD|SPP) + 
                (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
              family = "binomial", 
              data = veg)
summary(m6)

m7 <- glmer(ABUND ~ (1|SPP) + BF_2000 + (0 + BF_2000|SPP) +
                poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
                poly(WDW, degree = 2) + (0 + WDW|SPP) + 
                (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
              family = "binomial",
              data = veg)
summary(m7)








#### Slocum Data ####

#Elevation. AIC = 20273.3
glmer.elev.o <- glmer(Presence ~ (1|SPP) + Elevation,
                    family = "binomial", 
                    data = veg.slocum)
summary(glmer.elev.o)

#Elevation + Elevation2 AIC = 20270.4
glmer.elev.o2 <- glmer(Presence ~ (1|SPP) + poly(Elevation, degree = 2), 
                     family = "binomial",
                     data = veg.slocum)
summary(glmer.elev.o2)

#Elevation + Elevation2 + Sites. AIC = 20071.4
glmer.elev.o3 <- glmer(Presence ~ (1|SPP) + Elevation +
                       (0 + Elevation|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial",
                     data = veg.slocum)
summary(glmer.elev.o3)

#Elevation + Elevation2 + Sites + RE. AIC = 17515.4 **BEST**
glmer.elev.o4 <- glmer(Presence ~ (1|SPP) + poly(Elevation, degree = 2) +
                       (0 + Elevation|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial",
                     data = veg.slocum)
summary(glmer.elev.o4)

###########################################################################
#WD. AIC = 20273.9
glmer.wd.o <- glmer(Presence ~ (1|SPP) + WD_Old,
                  family = "binomial", 
                  data = veg.slocum)
summary(glmer.wd.o)

#Elevation + Elevation2 AIC = 20269.8
glmer.wd.o2 <- glmer(Presence ~ (1|SPP) + poly(WD_Old, degree = 2), 
                   family = "binomial",
                   data = veg.slocum)
summary(glmer.wd.o2)

#WD + WD2 + Sites. AIC = 20072.6
glmer.WD.o3 <- glmer(Presence ~ (1|SPP) + WD_Old +
                     (0 + WD_Old|SPP)+
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial",
                   data = veg.slocum)
summary(glmer.WD.o3)

#WD + WD2 + Sites + RE. AIC = 17728.6
glmer.wd.o4 <- glmer(Presence ~ (1|SPP) + poly(WD_Old, degree = 2) + 
                     (0 + WD_Old|SPP)+
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial",
                   data = veg.slocum)
summary(glmer.wd.o4)

###########################################################################
#BF. AIC = 20306.2
glmer.bf.o <- glmer(Presence ~ (1|SPP) + BF_1980,
                  family = "binomial", 
                  data = veg.slocum)
summary(glmer.bf.o)

#BF + BF2 AIC = 20272.0
glmer.bf.o2 <- glmer(Presence ~ (1|SPP) + poly(BF_1980, degree = 2), 
                   family = "binomial",
                   data = veg.slocum)
summary(glmer.bf.o2)

#BF + BF2 + Sites. AIC = 20077.8
glmer.bf.o3 <- glmer(Presence ~ (1|SPP) + BF_1980 + 
                     (0 + BF_1980|SPP) +
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial",
                   data = veg.slocum)
summary(glmer.bf.o3)

#BF + BF2 + Sites + RE. AIC = 19227.4
glmer.bf.o4 <- glmer(Presence ~ (1|SPP) + poly(BF_1980, degree = 2) +
                     (0 + BF_1980|SPP) + 
                     (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                   family = "binomial",
                   data = veg.slocum)
summary(glmer.bf.o4)

###############################################################################
#Inundation AIC = 20310.5
glmer.inun.o <- glmer(Presence ~ (1|SPP) + IN_20_Old,
                    family = "binomial", 
                    data = veg.slocum)
summary(glmer.inun.o)

#Inundation + SPP + Site, AIC = 19463.7
glmer.inun.o2 <- glmer(Presence ~ (1|SPP) + poly(IN_20_Old, degree = 2),
                    family = "binomial", 
                    data = veg.slocum)
summary(glmer.inun.o2)

#Inundation + SPP + Site, AIC = 19463.7
glmer.inun.o3 <- glmer(Presence ~ (1|SPP) + IN_20_Old +  
                       (0 + IN_20_Old|SPP),
                       family = "binomial", 
                       data = veg.slocum)
summary(glmer.inun.o3)

#Inundation + SPP + Site, AIC = 19463.7
glmer.inun.o4 <- glmer(Presence ~ (1|SPP) + poly(IN_20_Old, degree = 2) + 
                       (0 + IN_20_Old|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                       family = "binomial", 
                       data = veg.slocum)
summary(glmer.inun.o4)

#####################################################
#Dry Season Water Window, AIC = 20276.9
glmer.dry.o <- glmer(Presence ~ (1|SPP) + WDD_Old,
                   family = "binomial", 
                   data = veg.slocum)
summary(glmer.dry.o)

#Dry Season Water Window + Dry2, AIC = 20276.9
glmer.dry.o2 <- glmer(Presence ~ (1|SPP) + poly(WDD_Old, degree = 2),
                     family = "binomial", 
                     data = veg.slocum)
summary(glmer.dry.o2)

#Dry Season Water Window + Full RE, AIC = 20276.9
glmer.dry.o3 <- glmer(Presence ~ (1|SPP) + WDD_Old +
                      (0 + WDD_Old|SPP) +
                      (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                      family = "binomial", 
                      data = veg.slocum)
summary(glmer.dry.o3)

#Dry Season Water Window + Dyr 2 + Full RE, AIC = 20276.9
glmer.dry.o4 <- glmer(Presence ~ (1|SPP) + poly(WDD_Old, degree = 2) +
                        (0 + WDD_Old|SPP) +
                        (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                      family = "binomial", 
                      data = veg.slocum)
summary(glmer.dry.o4)

#######################################################
#Fire Season Water Window, AIC = 20282.1
glmer.fire.o <- glmer(Presence ~ (1|SPP) + WDF_Old,
                    family = "binomial", 
                    data = veg.slocum)
summary(glmer.fire.o)

#Fire Season Water Window + Fire2, AIC = 20280.4
glmer.fire.o2 <- glmer(Presence ~ (1|SPP) + poly(WDF_Old, degree = 2),
                     family = "binomial", 
                     data = veg.slocum)
summary(glmer.fire.o2)

#Fire Season Water Window + Full RE, AIC = 20072.5
glmer.fire.o3 <- glmer(Presence ~ (1|SPP) + WDF_Old + 
                       (0 + WDF_Old|SPP) + 
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial", 
                     data = veg.slocum)
summary(glmer.fire.o3)

#Fire Season Water Window + Fire2 + Full RE, AIC = 17672.6
glmer.fire.o4 <- glmer(Presence ~ (1|SPP) + poly(WDF_Old, degree = 2) + 
                       (0 + WDF_Old|SPP) +
                       (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                     family = "binomial", 
                     data = veg.slocum)
summary(glmer.fire.o4)

##########################################################
#Wet Season Water Window, AIC = 20282.1
glmer.wet.o <- glmer(Presence ~ (1|SPP) + WDW_Old,
                      family = "binomial", 
                      data = veg.slocum)
summary(glmer.wet.o)

#Wet Season Water Window + Wet2, AIC = 20280.4
glmer.wet.o2 <- glmer(Presence ~ (1|SPP) + poly(WDW_Old, degree = 2),
                       family = "binomial", 
                       data = veg.slocum)
summary(glmer.wet.o2)

#Wet Season Water Window + Full RE, AIC = 20072.5
glmer.wet.o3 <- glmer(Presence ~ (1|SPP) + WDW_Old + 
                         (0 + WDW_Old|SPP) + 
                         (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                       family = "binomial", 
                       data = veg.slocum)
summary(glmer.wet.o3)

#Wet Season Water Window + Wet2 + Full RE, AIC = 17672.6
glmer.wet.o4 <- glmer(Presence ~ (1|SPP) + poly(WDW_Old, degree = 2) + 
                         (0 + WDW_Old|SPP) +
                         (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
                       family = "binomial", 
                       data = veg.slocum)
summary(glmer.wet.o4)

#############################################################
#Water Window all, AIC = 20246.9
glmer.water.win.o <- glmer(Presence ~ (1|SPP) + 
                           WDD_Old + 
                           WDF_Old + 
                           WDW_Old,
                         family = "binomial", 
                         data = veg.slocum)
summary(glmer.water.win.o)

##########################################
######## Combine my Variables now ########
##########################################
# Work through the variables. 
m1.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(m1.o)

m2.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(m2.o)

m3.o <- glmer(Presence~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(IN_20_Old, degree = 2) + (0 + IN_20_Old|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(m3.o)

m4.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(WD_Old, degree = 2) + (0 + WD_Old|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(m4.o)

m5.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(WDF_Old, degree = 2) + (0 + WDF_Old|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(m5.o)

m6.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
                poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
                poly(WDD_Old, degree = 2) + (0 + WDD_Old|SPP) + 
                (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
              family = "binomial", 
              data = veg.slocum)
summary(m6.o)

m7.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
                poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
                poly(WDW_Old, degree = 2) + (0 + WDW_Old|SPP) + 
                (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
              family = "binomial",
              data = veg.slocum)
summary(m7.o)

######################################################################
# Remove RE systematically to determine mlm pvalues
z1.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(IN_20_Old, degree = 2) + (0 + IN_20_Old|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(z1.o)

#Without Burn Fequency
z2.o <- glmer(Presence ~ (1|SPP) + BF_1980 +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(IN_20_Old, degree = 2) + (0 + IN_20_Old|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(z2.o)

#Without Elevation
z3.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              poly(Elevation, degree = 2) + 
              poly(IN_20_Old, degree = 2) + (0 + IN_20_Old|SPP) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(z3.o)

#Without Inundation
z4.o <- glmer(Presence ~ (1|SPP) + BF_1980 + (0 + BF_1980|SPP) +
              poly(Elevation, degree = 2) + (0 + Elevation|SPP) + 
              poly(IN_20_Old, degree = 2) + 
              (1|Block) + (1|Block:Plot2) + (1|Block:Plot2:Submodule2),
            family = "binomial", 
            data = veg.slocum)
summary(z4.o)

LL <- c(deviance(z1.o), deviance(z2.o), deviance(z3.o), deviance(z4.o))

x <- 1 - pchisq(LL[2] - LL[1], 1)
#these are chi-square distribution functions). 
mlm.pvals <- c(
  1 - pchisq(LL[2] - LL[1], 1),
  1 - pchisq(LL[3] - LL[1], 1),
  1 - pchisq(LL[4] - LL[1], 1))

mlm.pvals

# compute residual effect of random effects environmental variables
nsite <- 176
nspp <- 193
glmer.resid.o <- glmer(Presence ~ (1 | SPP) + BF_1980 + 
                       poly(Elevation, degree = 2) + 
                       poly(IN_20_Old, degree = 2),
                     data = veg.slocum,
                     family = "binomial"
)

MLM.fitted <- array(fitted(z1.o) - fitted(glmer.resid.o), c(nsite, nspp))

rownames(MLM.fitted) = c(1:176)
colnames(MLM.fitted) = names(veg.names)[1:193]

# standardize over spp
MLM.fitted.standard <- MLM.fitted
for (j in 1:nspp)
  MLM.fitted.standard[, j] <-
  (MLM.fitted[, j] - mean(MLM.fitted[, j])) /
  sd(MLM.fitted[, j])

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[, 1:2]

# environmental variables (only those with random effects)
envir.vars <- cbind(veg.slocum$Elevation, veg.slocum$BF_1980, veg.slocum$IN_20_Old)
mlm.envir <- NULL
for (j in 1:3)
  mlm.envir <-
  cbind(mlm.envir, envir.vars[, j] * mlm.fit[, 1], envir.vars[, j] * mlm.fit[, 2])

envir.points <- matrix(colMeans(mlm.envir), nrow = ncol(envir.vars), ncol = 2, byrow = T)

# plot mlm
par(mfcol = c(1, 1))
plot(-mlm.fit,
     xlab = "PC1",
     ylab = "PC2",
     type = "n")
text(-mlm.fit, label = 1:nsite, cex = .5)

arrow.coordMLM <- cbind(array(0, dim(envir.points)), -envir.points)

arrows(
  arrow.coordMLM[, 1],
  arrow.coordMLM[, 2],
  arrow.coordMLM[, 3],
  arrow.coordMLM[, 4],
  col = "black",
  length = 0.1
)

text(
  1.3 * -envir.points,
  label = c("Elevation", "BF", "IN_Old"),
  cex = .7
)

title(main = "PCA, 1990-2000, BF, Elev, Total Inundation")

####################################################################
fixef(z1)
spp.bf.1 <- ranef(z1)$SPP$BF + 0.180
spp.TI.1 <- ranef(z1)$SPP$Inundation_20 - 153.460
spp.Elev.1 <- ranef(z1)$SPP$Elevation - 31.96
names <- unique(veg.names)
z1.spp <- cbind(names, spp.bf.1, spp.TI.1, spp.Elev.1)
z1.spp <- as.data.frame(z1.spp)
write_csv(z1.spp, "Present_SPP_RE.csv")

spp.bf <- ranef(Model.2.o)$SPP$BF + 0.273
spp.WD <- ranef(Model.2.o)$SPP$WD + 44.574
spp.TI <- ranef(Model.2.o)$SPP$Inundation_2_Old -1.529
slocum.names <- veg.slocum[,15]
slocum.names<- unique(slocum.names)
Model.2.o.spp <- cbind(slocum.names, spp.bf, spp.TI, spp.WD)
Model.2.o.spp <- as.data.frame(Model.2.o.spp)
write_csv(Model.2.o.spp, "Past_SPP_RE.csv")
