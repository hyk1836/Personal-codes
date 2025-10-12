##read data
data <- readRDS("C:/Users/cmc/Desktop/과제/1.기후변화_자가면역/NEDIS/1_data/1_rds/20251010_casecross/NEDIS_casecross_primary_final_20251010.rds")

##Reading R packages
library(dlnm)
library(splines)
library(survival)

head(data)

# Set formula
fmla<-as.formula(AID11 ~ cbtem + strata(case_id) + ns(hum, df = 3) + dow)

###DLNM function
## 7-day lag with 2 degrees of freedom
##################################################################
dftrend <-round(as.numeric((diff(range(data$ddate))/365.25)* 7))   # df = 7
btrend <- ns(data$ddate, knots=equalknots(data$ddate, dftrend-1))

lag <- 7
lagnk <- 2
varknots <- quantile(data$tem, probs = c(.10, .75, .90), na.rm = TRUE)
lagknots <- logknots(lag, lagnk)
degree = 2 

#lag 7까지 포함하는 crossbasis -> lag 4, 14, 21 등등 수정하여 여러 모델 비교할 수 있음
cbtem <- crossbasis(data$tem, lag = lag, 
                    argvar = list(fun = "bs", degree = degree, knots = varknots),
                    arglag = list(knots = lagknots))

mod4 <- clogit(fmla,data=data,na.action=na.exclude)

# 예측/요약: 중심화값 설정 (MMT)
pred0 <- crosspred(cbtem, mod4)
mmt   <- pred0$predvar[which.min(pred0$allRRfit)]  # 최소위험기온

pred   <- crosspred(cbtem, mod4, cen = mmt, by = 1)
red  <- crossreduce(cbtem, mod4, cen = mmt, by = 1)

## Removing Temperature Extremes
wbgt <- data[,8:(8+7)]
p99<-quantile(wbgt, probs = 0.999,na.rm = T)
p1<-quantile(wbgt, probs = 0.001,na.rm = T)
wbgt[wbgt > p99] <- NA
wbgt[wbgt < p1] <- NA
bound<-c(round(p1,1),round(p99,1))

##
t1<-round(quantile(wbgt,0.01,na.rm = T),1)
t2<-round(quantile(wbgt,0.025,na.rm = T),1)
t3<-round(quantile(wbgt,0.975,na.rm = T),1)
t4<-round(quantile(wbgt,0.99,na.rm = T),1)

## building up cross-bases
argvar <- list(fun="ns", knots = equalknots(wbgt,fun = "ns",df=4),
               Bound=bound)
cb <- crossbasis(wbgt,lag=lag,argvar=argvar,
                 arglag=list(knots=logknots(lag,lagnk)))
model<-clogit(fmla,data=data,na.action=na.exclude)
pred <- crosspred(cb,model,by=0.1,cumul = TRUE)
## centering value unspecified. Automatically set to 10
cen<-pred$predvar[which.min(pred$allfit)]
print(cen) ## MMT
## [1] 15
pred_S <- crosspred(cb,model,by=0.1,cumul = TRUE,cen=cen,
                    at=c(seq(bound[1],bound[2],by=0.1),t3,t4))

##Calculating OR
x1<-t(as.matrix(round(with(pred,cbind(allRRfit,allRRlow,allRRhigh)[as.character(t1),]),3)))
x2<-t(as.matrix(round(with(pred,cbind(allRRfit,allRRlow,allRRhigh)[as.character(t2),]),3)))
x3<-t(as.matrix(round(with(pred,cbind(allRRfit,allRRlow,allRRhigh)[as.character(t3),]),3)))
x4<-t(as.matrix(round(with(pred,cbind(allRRfit,allRRlow,allRRhigh)[as.character(t4),]),3)))
result <- data.frame(rbind(x1,x2,x3,x4))
rownames(result) <- c(t1,t2,t3,t4)

##Exposure-response curve
# Define the x-axis label with the degree symbol
xlab <- expression(paste("WBGT (", degree, "C)"))

# Set up plotting parameters
par(mar = c(4, 4, 3, 1), cex.axis = 0.9, mgp = c(2.5, 1, 0), cex.lab = 1)

# Ensure that pred_S is properly defined and that the plot function supports your arguments
if (exists("pred_S")) {
  plot(pred_S, "overall", col = 2, cex.axis = 1, xlab = xlab, ylab = "Odds ratio", ylim = c(0.8, 2.0), 
       main = paste0("Exposure-response curve"), cex.main = 1, ci = "area", ci.arg = list(col = grey(0.85)))
  lines(pred_S, col = 2, lwd = 2)
} else {
  warning("The object 'pred_S' does not exist.")
}

# Save the plot as a TIFF file
tiff("figure_ER.tiff", res = 300, height = 1200, width = 1200)

# Reproduce the plot for saving
par(mar = c(4, 4, 3, 1), cex.axis = 0.9, mgp = c(2.5, 1, 0), cex.lab = 1)
plot(pred_S, "overall", col = 2, cex.axis = 1, xlab = xlab, ylab = "Odds ratio", ylim = c(0.8, 2.0), 
     main = paste0("Exposure-response curve"), cex.main = 1, ci = "area", ci.arg = list(col = grey(0.85)))
lines(pred_S, col = 2, lwd = 2)

# Close the TIFF device
dev.off()
## quartz_off_screen 
##                 2
## lag pattern curve
# Define the x-axis label
xlab <- "Lag (day)"
ylab <- "Odds ratio"

# Set up plotting parameters
par(mar = c(4, 4, 3, 1), cex.axis = 0.9, mgp = c(2.5, 1, 0), cex.lab = 1)

# Ensure that pred_S and t4 are properly defined
if (exists("pred_S") && exists("t4")) {
  plot(pred_S, "slices", var = t4,
       xlim = c(0, 7), cex.main = 1,
       col = 2, lwd = 2,
       xlab = xlab,
       ylab = ylab, cex.lab = 1,
       ci = "area", ci.arg = list(col = grey(0.85)),
       main = paste0("Lag pattern"))
  lines(pred_S, col = 1, lwd = 2)
} else {
  warning("The objects 'pred_S' or 't4' do not exist.")
}

# Save the plot as a TIFF file
tiff("figure_lag.tiff", res = 300, height = 1200, width = 1200)

# Set up plotting parameters
par(mar = c(4, 4, 3, 1), cex.axis = 0.9, mgp = c(2.5, 1, 0), cex.lab = 1)

# Ensure that pred_S and t4 are properly defined
if (exists("pred_S") && exists("t4")) {
  plot(pred_S, "slices", var = t4,
       xlim = c(0, 7), cex.main = 1,
       col = 2, lwd = 2,
       xlab = xlab,
       ylab = ylab, cex.lab = 1,
       ci = "area", ci.arg = list(col = grey(0.85)),
       main = paste0("Lag pattern"))
  lines(pred_S, col = 1, lwd = 2)
} else {
  warning("The objects 'pred_S' or 't4' do not exist.")
}

# Close the TIFF device
dev.off()
## quartz_off_screen 
##                 2