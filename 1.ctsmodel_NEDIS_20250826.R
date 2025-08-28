# 0. LOAD FILES
setwd("C:/Users/Main/Desktop/업무_효연/NEDIS/20250826")
NEDIS <-readRDS("NEDIS_timeseries_20250821.rds")

# LOAD LIBRARIES
library(xtable)
library(data.table) 
library(gnm)
library(splines)
library(dlnm) 

# REGIONS
regions <- as.character(unique(NEDIS$ptmigucd_2023))

# CREATE A LIST WITH THE REGIONAL SERIES
datalist <- lapply(regions,function(x) NEDIS[NEDIS$ptmigucd_2023==x,])
names(datalist) <- regions
m <- length(regions)
NEDIS$time <- as.numeric(factor(NEDIS$ddate))
NEDIS$dow <- factor(NEDIS$dow)

# LIST OF DATAFRAMES  (중복으로 text 처리함)
# datalist <- lapply(regions, function(region) NEDIS[NEDIS$ptmigucd_2023==region,])
# names(datalist) <- regions

# TEMPERATURE RANGES
tem_ranges <- t(sapply(datalist, function(x) range(x$tem,na.rm=T)))

####################################################################
# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

####################################################################
# MAIN MODEL ON CASE TIME SERIES DATA
memory.limit(size=100000)
group <- factor(paste(NEDIS$ptmigucd_2023, NEDIS$year, sep="-"))

# DEFINE THE STRATA 
NEDIS <- data.table(NEDIS)
NEDIS[, stratum := factor(paste(ptmigucd_2023, year, month, sep=":"))]

# RUN THE MODEL
# NB: EXCLUDE EMPTY STRATA, OTHERWISE BIAS IN gnm WITH quasipoisson
NEDIS[,  keep := sum(AID23) > 0, by = stratum]

##################################################################
dftrend <-round(as.numeric((diff(range(NEDIS$ddate))/365.25)* 7))   # df = 7
btrend <- ns(NEDIS$ddate, knots=equalknots(NEDIS$ddate, dftrend-1))

lag <- 21
lagnk <- 2
varknots <- quantile(NEDIS$tem, probs = c(.10, .75, .90), na.rm = TRUE)
lagknots <- logknots(lag, lagnk)
degree = 2 

#lag 7까지 포함하는 crossbasis -> lag 4, 14, 21 등등 수정하여 여러 모델 비교할 수 있음
cbtem <- crossbasis(NEDIS$tem, lag = lag, group = group,
                    argvar = list(fun = "bs", degree = degree, knots = varknots),
                    arglag = list(knots = lagknots))


mod4 <- gnm(AID23 ~ cbtem + btrend + ns(hum, df = 3) + dow, 
            eliminate=stratum, family = "quasipoisson", data = NEDIS, subset = keep)

# 예측/요약: 중심화값 설정 (MMT)
pred0 <- crosspred(cbtem, mod4)
mmt   <- pred0$predvar[which.min(pred0$allRRfit)]  # 최소위험기온

pred   <- crosspred(cbtem, mod4, cen = mmt, by = 1)
red  <- crossreduce(cbtem, mod4, cen = mmt, by = 1)

####################################################################
###                  결과 plot 그리기
####################################################################
## 01 main plot
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
layout(matrix(1:2,ncol=2))

plot(pred, "overall", xlab = "Temperature", ylab = "RR", col = "blue",
     main = "overall cumulative exposure-response")
legend("topleft", c("Lag 0-21(with 95%CI)"), lty=c(1,2), lwd=1.5, col=c("blue"),bty="n",inset=c(0.05,0.1),cex=0.8)
plot(pred, xlab="Temperature", zlab="RR", ltheta=240, lphi=60, cex.axis=0.8,
     main = "Temperature: exposure-lag-response", col=4)

## 02 Exposure-response at a certain lag point
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
layout(matrix(1:2,ncol=2))

red_lag<-crossreduce(cbtem,mod4,type="lag",value=14,cen=mmt)             ##Lag=0
plot(red_lag,ylab="RR",xlab="Temperature", col = "blue")
legend("topleft", c("Lag 14(with 95%CI)"), lty=c(1,2), lwd=1.5, col=c("blue"),bty="n",inset=c(0.05,0.1),cex=0.8)

red_lag<-crossreduce(cbtem,mod4,type="lag",value=21,cen=mmt)             ##Lag=2
plot(red_lag,ylab="RR",xlab="Temperature", col = "blue")
legend("topleft", c("Lag 21(with 95%CI)"), lty=c(1,2), lwd=1.5, col=c("blue"),bty="n",inset=c(0.05,0.1),cex=0.8)

## 03 Lag-response at a certain exposure point (30도, -10도)
red_exp<-crossreduce(cbtem,mod4,type="var",value=30,cen=mmt) 
plot(red_exp,ylab="RR",xlab="Lag days", col = 2)
legend("topleft", legend = expression(30 * degree * "C (with 95% CI)"), lty=c(1,2), lwd=1.5, col=2, bty="n",inset=c(0.05,0.1),cex=0.8)

red_exp<-crossreduce(cbtem,mod4,type="var",value=-10,cen=mmt) 
plot(red_exp,ylab="RR",xlab="Lag days", col = 2)
legend("topleft", legend = expression(-10 * degree * "C (with 95% CI)"), lty=c(1,2), lwd=1.5, col=2, bty="n",inset=c(0.05,0.1),cex=0.8)

####################################################################
##                      결과 숫자로 추출하기
####################################################################
## 01 기온별 RR 구하기
pred_tem <- crosspred(cbtem, mod4, cen = mmt, at = c(-10, 30))

result_tem <- matrix(NA, 2, 3)
colnames(result_tem) <- c("RRfit", "RRlow", "RRhigh")
rownames(result_tem) <- c('-10°C', '30°C')

result_tem[1, ] <- c(RRfit = pred_tem$allRRfit["-10"], RRlow = pred_tem$allRRlow["-10"], RRhigh = pred_tem$allRRhigh["-10"])
result_tem[2, ] <-c(RRfit = pred_tem$allRRfit["30"],  RRlow = pred_tem$allRRlow["30"],  RRhigh = pred_tem$allRRhigh["30"])

## 02 퍼센타일 별 RR 구하기
q <- quantile(NEDIS$tem, probs = c(.05, .95), na.rm = TRUE)
pred_q <- crosspred(cbtem, mod4, cen = mmt, at = q)

result_q <- matrix(NA, 2, 3)
colnames(result_q) <- c("RRfit", "RRlow", "RRhigh")
rownames(result_q) <- c('P5', 'P95')

result_q[1, ] <- c(RRfit = pred_q$allRRfit[1], RRlow = pred_q$allRRlow[1], RRhigh = pred_q$allRRhigh[1])
result_q[2, ] <-c(RRfit = pred_q$allRRfit[2],  RRlow = pred_q$allRRlow[2],  RRhigh = pred_q$allRRhigh[2])
