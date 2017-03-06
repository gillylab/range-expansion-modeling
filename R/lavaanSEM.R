## SEM analysis with lavaan for Humboldt squid
## Stewart et al. 2014 Global Change Biology: doi: 10.1111/gcb.12502

## setup ----

library(lavaan) # install.packages('devtools'); devtools::install_github('yrosseel/lavaan')
library(tidyverse) # install.packages('tidyverse')

## load data 
load("data/mysquid.RData")
head(mysquid); summary(mysquid) # see data/README.md for column headers

## preparations for SEM ----

## curate data
mysquid <- mysquid %>%
  mutate(Hake   = log(Hake+0.1), 
         Loligo = log(Loligo+0.1),
         Myct   = log(Myct+0.1),
         UIwin  = UIwin/1000, # otherwise it is too large and causes problems.
         Strat = as.numeric(Strat),
         SqMonth = as.numeric(as.character(factor(SqMonth))), # factors to numeric
         SqYear = as.numeric(as.character(factor(SqYear))),
         ROV_ID = as.numeric(as.character(factor(ROV_ID))))

## correlation and covariance matrices
cormat <-  cor(mysquid)
sds <-  sapply(mysquid,sd)
covmat <-  cor2cov(cormat, sds, names = NULL)
cormat
covmat

## model definition ----

squid.path <- '
  NPGO ~~ NOI 
  NOI ~~ PDO
  UIwin ~ NOI + PDO + NPGO
  tcstren ~ UIwin + Distkm + NOI + PDO + NPGO
  DepthOMZ ~ UIwin + Distkm + tcstren + NOI + PDO + NPGO
  DepthOMZ ~~ tcstren 
  Hake ~ Distkm + tcstren + DepthOMZ          
  Myct ~ Distkm + tcstren + DepthOMZ 
  Squid ~ Distkm + Hake  + Myct + DepthOMZ
'


## bootstrapping ----

# currently runs using only positives (fit.boot.pos)
fit.boot.pos<-sem(squid.path, data=mysquid[mysquid$Squid>0,], fixed.x=F, test="boot", se="boot", bootstrap=1000) 
summary(fit.boot.pos, standardized = TRUE)
coef(fit.boot.pos)
parameterEstimates(fit.boot.pos)
modificationIndices(fit.boot.pos)


## zero inflation Jarrett Skype Call, Feb 26 2013
require(ggplot2)
require(pscl)
require(boot)

s1 <- zeroinfl(round(Squid) ~ Distkm + tcstren + DepthOMZ + Hake + Loligo + Myct, 
               data = mysquid)

dmod <- lm(Squid ~ (Distkm + tcstren + DepthOMZ + Hake + Loligo + Myct)* SquidDetect, 
           data = mysquid)
summary(dmod)


#######
# Satorra-Bentler Correction --does not work with my data as is (see Yves comment)
########
fit.sb <- sem(squid.path, data = mysquid[mysquid$Squid>0,], 
              estimator="mlm", fixed.x=F) 
summary(fit.sb, standardized = T, fit.measures=T)
modificationIndices(fit.sb)
modindices(fit.sb)


## detection probit by Jarrett Feb 24 2013 ----

mysquid$SquidDetect <- as.numeric(mysquid$Squid>0)

squid.pathDetect <- '
 UIwin ~ NOI
 tcstren ~ UIwin + Distkm
 DepthOMZ ~ UIwin + Distkm + tcstren
 Loligo ~ Distkm + tcstren
 Hake ~ Distkm + tcstren + DepthOMZ
 Myct ~ Distkm + tcstren + DepthOMZ
 SquidDetect ~ Distkm + tcstren + DepthOMZ + Hake + Loligo + Myct
 Squid ~ Distkm + tcstren + DepthOMZ + Hake + Loligo + Myct + SquidDetect
'
detectFull_probit <- sem(squid.pathDetect, data = mysquid, 
                         fixed.x=F,test="boot", se="boot", bootstrap=5)
summary(detectFull_probit, standardized = T, fit.measures=T)

## probit ----

# Model 1
squid.pathProbit <- '
  NOI ~ PDO
  #NOI ~~ PDO
  UIwin ~ NOI + PDO 
  tcstren ~ UIwin + Distkm 
  DepthOMZ ~ UIwin + Distkm + tcstren
  #DepthOMZ ~~ tcstren 
  Hake ~ Distkm + tcstren + DepthOMZ         
  Myct ~ Distkm + tcstren + DepthOMZ
  SquidDetect ~ Distkm + Hake + Myct + DepthOMZ
  #Squid ~ Distkm + Hake + Myct + DepthOMZ + SquidDetect
'
fitFull_probitBoot <- sem(squid.pathProbit, data = mysquid, 
                          fixed.x=F, test="boot", se="boot", bootstrap=1000)
summary(fitFull_probitBoot, standardized = T, fit.measures=T)
#modificationIndices(fitFull_probitBoot)

fitFull_probit <- sem(squid.pathProbit, data = mysquid, 
                      estimator="wlsmv", ordered="SquidDetect") 
summary(fitFull_probit, standardized = T, fit.measures=T)
modificationIndices(fitFull_probit) # may have trouble with bootstrapped 
