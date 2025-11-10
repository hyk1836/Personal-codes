## 패키지
library(dplyr)
library(zoo)
library(mgcv)
library(mediation)

## 1. 데이터 불러오기
data <- readRDS("C:/Users/Main/Desktop/업무_효연/기후보건영향평가/20251104_mediation troubleshoot/deathcount_sido_O3_short_2020-2023_we.rds")
mmt_tbl <- read.csv("C:/Users/Main/Desktop/업무_효연/기후보건영향평가/20251104_mediation troubleshoot/MMT_blup_by_sido.csv")

## 2. 시도별 MMT 붙이고 서울만 추출
data2 <- data %>%
  mutate(sido = as.integer(sido)) %>%
  left_join(mmt_tbl, by = "sido")

data_seoul <- data2 %>%
  filter(sido == 11) %>%
  arrange(ddd)   # 시계열 정렬

## 3. 기온 이동평균(0–1일) 만들기 + segment는 일단 참고용만 생성
data_seoul <- data_seoul %>%
  mutate(
    tem_ma01 = zoo::rollapply(tem, width = 2, FUN = mean,
                              align = "right", fill = NA),
    segment  = if_else(tem >= mmt, "over", "under")
  )

## 4. 처리값(treat/control) 결정
tm_range  <- range(data_seoul$tem_ma01, na.rm = TRUE)
mmt_seoul <- unique(data_seoul$mmt)[1]

treat_val   <- mmt_seoul + 1
control_val <- mmt_seoul

# 자료 범위를 벗어나면 안쪽으로 보정
if (treat_val > tm_range[2]) treat_val <- tm_range[2]
if (control_val < tm_range[1]) control_val <- tm_range[1]

## 5. 매개모형: O3 ~ 기온 + 공변량
med.fit <- gam(
  o3_8h_daily_max ~ tem_ma01 + s(ddd, k = 7 * 5) + hum + mean_tsr + DOW,
  family = gaussian(),
  data   = data_seoul,
  method = "REML",
  model  = TRUE
)
summary(med.fit)

## 6. 결과모형: 사망 ~ O3 + 기온 + 공변량
out.fit <- gam(
  TOT_nonacc ~ o3_8h_daily_max + tem_ma01 + DOW + s(ddd, k = 7 * 5) + hum + mean_tsr,
  family = poisson(),
  data   = data_seoul,
  method = "REML",
  model  = TRUE
)
summary(out.fit)

## 7. 인과 매개분석 실행
set.seed(2025)

mo <- mediate(
  model.m = med.fit,
  model.y = out.fit,
  treat   = "tem_ma01",
  mediator = "o3_8h_daily_max",
  treat.value   = treat_val,
  control.value = control_val,
  boot = TRUE,
  sims = 1000
  # sim.outcomes = TRUE 는 mgcv와 함께일 때 에러가 날 수 있어 일단 뺐습니다.
)

summary(mo)

## 1) mediator 모형으로 두 가지 오존을 예측
nd_control <- data_seoul
nd_control$tem_ma01 <- control_val
o3_control <- predict(med.fit, newdata = nd_control, type = "response")

nd_treated <- data_seoul
nd_treated$tem_ma01 <- treat_val
o3_treated <- predict(med.fit, newdata = nd_treated, type = "response")

## 2) outcome 모형으로 세 가지 세계 예측

# (a) control 세계: 기온=control, 오존=control에서 예측된 값
new0 <- data_seoul
new0$tem_ma01 <- control_val
new0$o3_8h_daily_max <- o3_control
mu0 <- predict(out.fit, newdata = new0, type = "response")

# (b) treated-but-mediator-frozen: 기온만 treat로 올리고 오존은 control 거 사용
new1m0 <- data_seoul
new1m0$tem_ma01 <- treat_val
new1m0$o3_8h_daily_max <- o3_control
mu1m0 <- predict(out.fit, newdata = new1m0, type = "response")

# (c) treated 세계: 기온=treat, 오존=treat에서 예측된 값
new1m1 <- data_seoul
new1m1$tem_ma01 <- treat_val
new1m1$o3_8h_daily_max <- o3_treated
mu1m1 <- predict(out.fit, newdata = new1m1, type = "response")

## 3) 각 시나리오의 평균
m0    <- mean(mu0,    na.rm = TRUE)
m1m0  <- mean(mu1m0,  na.rm = TRUE)
m1m1  <- mean(mu1m1,  na.rm = TRUE)

## 4) RR와 %change 계산
RR_total <- m1m1 / m0        # 전체효과 RR
RR_NDE   <- m1m0 / m0        # 직접효과 RR
RR_NIE   <- m1m1 / m1m0      # 간접효과 RR

PC_total <- (RR_total - 1) * 100
PC_NDE   <- (RR_NDE   - 1) * 100
PC_NIE   <- (RR_NIE   - 1) * 100

RR_total; PC_total
RR_NDE;   PC_NDE
RR_NIE;   PC_NIE

