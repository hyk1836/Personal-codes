# 0. PREPARATION
setwd("C:/Users/Main/Desktop/업무_효연/기후보건영향평가/현행화")
data <- readRDS("deathcount_sido_PM25_short_2015-2019_we.rds")

# LOAD THE LIBRARIES
library(mediation)
library(mgcv)
library(dlnm)
library(bda)
library(metafor)
library(splines)
library(zoo)

# CREATE SIDO-SPECIFIC SUBSETS
seoul<-subset(data,sido==11)
busan<-subset(data,sido==21)
daegu<-subset(data,sido==22)
incheon<-subset(data,sido==23)
gwangju<-subset(data,sido==24)
daejeon<-subset(data,sido==25)
ulsan<-subset(data,sido==26)
sejong<-subset(data,sido==29)
gyeonggi<-subset(data,sido==31)
gangwon<-subset(data,sido==32)
chungbuk<-subset(data,sido==33)
chungnam<-subset(data,sido==34)
jeonbuk<-subset(data,sido==35)
jeonnam<-subset(data,sido==36)
gyeongbuk<-subset(data,sido==37)
gyeongnam<-subset(data,sido==38)
jeju<-subset(data,sido==39)

# 1. FIND THE MINIMUM-MORTALITY TEMPERATURE (MMT) W/ DLNM

# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

### Example of DLNM ###
dfseas <- 7 #Seasonality
spltime <- ns(seoul$death_date, df=dfseas*length(unique(seoul$year)))

# Modeling spec
varfun = "bs" #Spline function for the exposure-response relationship
varper <- c(10,75,90) #Location of Knots
degree=2 #Degree for the bs spline

lag <- 4 #Lag days
lagnk <- 2 #Df of the lag-response relationship

# Formula
formula <- TOT_nonacc ~ cb + ns(hum, df = 3) + DOW + spltime   

#Run a model
argvar <- list(fun=varfun,degree=degree,knots=quantile(seoul$tem,varper/100,na.rm=T))
cb <- crossbasis(seoul$tem,lag=lag, argvar=argvar, 
                 arglag=list(knots=logknots(lag,lagnk)))

model <-glm(formula ,family=quasipoisson,data=seoul)
fqaic(model)

pred<-crosspred(cb,model) #Risk prediction
mmt<-pred$predvar[which.min(pred$allRRfit)] #Optimal point (minimum mortality temp.)
pred_fin<-crosspred(cb,model,cen=mmt) #Centering
plot(pred_fin) #3D plot

#Overall effect (lag-cumulative exposure-response relationship)
red<-crossreduce(cb,model,cen=mmt)
plot(red,ylab="RR",xlab="Temperature")

# 2. SEPARATE THE DATA INTO overMMT AND underMMT 

seoul_over<-subset(seoul,tem >= 27)
seoul_under<-subset(seoul,tem < 27)

# 3. CONDUCT MEDIATION ANALYSIS

# Create moving averages
seoul_over$tem_ma7 <- zoo::rollapply(seoul_over$tem, width=7, FUN=mean, align="right", fill=NA)
seoul_over$tem_ma1 <- zoo::rollapply(seoul_over$tem, width=1, FUN=mean, align="right", fill=NA)

# Equation 1. the association between daily mean temperature and daily mean PM2.5 concentration
med.fit<-gam(PM25_mean~tem+s(ddd,k=7*10)+hum,family=gaussian(), data=seoul_over)
summary(med.fit)
plot(med.fit)

# Equation 2. the association between the daily mean temperature and daily number of deaths
out.fit<-gam(TOT_nonacc~PM25_mean+tem+DOW+s(ddd,k=7*10)+hum,family=poisson(), data=seoul_over)
summary(out.fit)
plot(out.fit)

# Mediation analysis
med.out<-mediate(med.fit, out.fit, treat="tem", mediator="PM25_mean", boot=T, sims=5000)
summary(med.out)
plot(med.out)

# 4. SOBEL TEST TO ACQUIRE P-VALUE

# 5. CONDUCT META-ANALYSIS TO OBTAIN POOLED ESTIMATES

# 6. ESTIMATE ATTRIBUTABLE FRACTION (AF)


