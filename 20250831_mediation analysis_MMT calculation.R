# ===== 0) 패키지 로드 =====
library(dplyr)
library(splines)
library(dlnm)
library(lubridate)
library(readr)

setwd("C:/Users/Main/Desktop/업무_효연/기후보건영향평가/현행화")
data <- readRDS("deathcount_sido_PM25_short_2015-2019_we.rds")

# ===== 1) 메타데이터/리스트 =====
SIDO  <- unique(data[c("sido")])
dlist <- split(data, data$sido)[SIDO$sido]

# ===== 2) 파라미터 =====
arglagtmean <- list(fun="strata", breaks=1)
lagtmean    <- 4
dftime      <- 7
fmod        <- TOT_nonacc ~ cbtmean + ns(hum, 3) + DOW + spltime

QAIC <- function(model){
  phi    <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  -2*loglik + 2*summary(model)$df[3]*phi
}

# ===== 3) 결과 저장용 폴더 생성 =====
setwd("C:/Users/Main/Desktop/업무_효연/기후보건영향평가/20250830")
if(!dir.exists("mmt_plots")) dir.create("mmt_plots")

# ===== 4) 시도별 MMT 추출 + 플롯 저장 =====
stage1list <- vector("list", length(dlist))

for(i in seq_along(dlist)){
  df_raw <- dlist[[i]]
  
  res <- try({
    df <- df_raw %>%
      mutate(
        DOW = factor(lubridate::wday(death_date, week_start=1),
                     levels=1:7,
                     labels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")),
        spltime = ns(as.numeric(death_date),
                     df=dftime * dplyr::n_distinct(year))
      ) %>%
      filter(!is.na(tem), !is.na(TOT_nonacc), !is.na(hum))
    
    if(nrow(df) < 30){
      return(list(sido=SIDO$sido[i], mmt=NA_real_, qaic=NA_real_,
                  conv=FALSE, knots=NA, notes="insufficient_rows"))
    }
    
    knotstmean  <- stats::quantile(df$tem, c(0.10,0.75,0.90), na.rm=TRUE)
    argvartmean <- list(fun="bs", knots=knotstmean, degree=2)
    cbtmean     <- crossbasis(df$tem, lag=lagtmean,
                              argvar=argvartmean, arglag=arglagtmean)
    
    mod <- glm(fmod, data=df, family=quasipoisson)
    
    if(!isTRUE(mod$converged)){
      return(list(sido=SIDO$sido[i], mmt=NA_real_, qaic=NA_real_,
                  conv=FALSE, knots=as.numeric(knotstmean),
                  notes="glm_not_converged"))
    }
    
    # 예측 & MMT
    pred   <- crosspred(cbtmean, mod, by=0.1)               
    mmt    <- pred$predvar[which.min(pred$allRRfit)]        
    pred_cen <- crosspred(cbtmean, mod, cen=mmt, by=0.1)    
    qaic   <- QAIC(mod)
    
    # --------- 플롯 자동 저장 ---------
    png(filename = file.path("mmt_plots",
                             paste0("RRcurve_", SIDO$sido[i], ".png")),
        width = 1200, height = 900, res = 150)
    plot(pred_cen, "overall",
         xlab="Temperature (°C)",
         ylab="Relative Risk (lag-cumulative)",
         main=paste0("Lag 0–", lagtmean, " | SIDO=", SIDO$sido[i]))
    abline(v=mmt, col="red", lty=2, lwd=2)  # MMT 표시
    dev.off()
    # ----------------------------------
    
    list(
      sido  = SIDO$sido[i],
      mmt   = as.numeric(mmt),
      qaic  = qaic,
      conv  = TRUE,
      knots = as.numeric(knotstmean),
      notes = "ok"
    )
  }, silent=TRUE)
  
  if(inherits(res, "try-error")){
    res <- list(sido=SIDO$sido[i], mmt=NA_real_, qaic=NA_real_,
                conv=FALSE, knots=NA,
                notes=paste0("error: ", as.character(res)[1]))
  }
  
  stage1list[[i]] <- res
}

names(stage1list) <- SIDO$sido

# ===== 5) 결과 테이블/저장 =====

mmt_tbl <- purrr::map_dfr(stage1list, ~{
  tibble::tibble(
    sido  = .x$sido,
    mmt   = .x$mmt,
    qaic  = .x$qaic,
  )
})

saveRDS(stage1list, "stage1_mmt_list.rds")
write_csv(mmt_tbl, "stage1_mmt_table.csv")

mmt_tbl

# ===== 6) 지역별 mmt 기준 데이터 분리 =====
# 1) sido별 mmt lookup 벡터 만들기
mmt_lookup <- setNames(mmt_tbl$mmt, mmt_tbl$sido)

# 2) dlist를 over/under로 분리해서 다시 저장
for(s in names(dlist)){
  df <- dlist[[s]]
  mmt <- mmt_lookup[[s]]
  
  if(is.na(mmt)) {
    dlist[[s]] <- list(over=NULL, under=NULL, mmt=NA)
  } else {
    dlist[[s]] <- list(
      over  = subset(df, tem >= mmt),
      under = subset(df, tem <  mmt),
      mmt   = mmt
    )
  }
}
  