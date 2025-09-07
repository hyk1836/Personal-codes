data <- readRDS("C:/Users/Main/Desktop/업무_효연/기후보건영향평가/현행화/deathcount_sido_PM25_short_2015-2019_we.rds")
mmt_tbl <- read.csv("C:/Users/Main/Desktop/업무_효연/기후보건영향평가/20250906/MMT_blup_by_sido.csv")

# 1. SEPARATE DATA

library(dplyr)

# 1) 원데이터에 시도별 mmt를 붙인다
data2 <- data %>%
  mutate(sido = as.integer(sido)) %>%
  left_join(mmt_tbl, by = "sido") %>%
  mutate(segment = if_else(tem >= mmt, "over", "under"))


# 2) 전체 over/under 한 번에 분리
data_over  <- data2 %>% filter(segment == "over")
data_under <- data2 %>% filter(segment == "under")

# 3) 시도별 리스트로 쪼개기
dlist_over  <- split(data_over,  data_over$sido)
dlist_under <- split(data_under, data_under$sido)

# (선택) 각 시도의 데이터 개수/확인
sapply(dlist_over,  nrow)
sapply(dlist_under, nrow)

# 3. CONDUCT MEDIATION ANALYSIS

# SET FUNCTION

library(mgcv)

run_models <- function(df, region_name){
  cat("\n====================================\n")
  cat("▶ Region:", region_name, "\n")
  cat("====================================\n")
  
  # Equation 1
  cat("\n--- Equation 1: Temperature → PM2.5 ---\n")
  med.fit <- gam(PM25_mean ~ tem + s(ddd, k=7*5) + hum, 
                 family = gaussian(), data = df, method = "REML")
  print(summary(med.fit))
  plot(med.fit, main = paste("Eq1:", region_name))
  
  # Equation 2
  cat("\n--- Equation 2: Temperature/PM2.5 → Deaths ---\n")
  out.fit <- gam(TOT_nonacc ~ PM25_mean + tem + DOW + s(ddd, k=7*5) + hum,
                 family = poisson(), data = df, method = "REML")
  print(summary(out.fit))
  plot(out.fit, main = paste("Eq2:", region_name))
  
  return(list(med.fit = med.fit, out.fit = out.fit))
}

# REPEAT OVER SEPARATE DATA(OVER/UNDER)

# 제외할 sido 번호
exclude_id <- "29"

# over 데이터 루프
over_names <- setdiff(names(dlist_over), exclude_id)
results_over <- lapply(over_names, function(s){
  run_models(dlist_over[[s]], s)
})
names(results_over) <- over_names

# under 데이터 루프
under_names <- setdiff(names(dlist_under), exclude_id)
results_under <- lapply(under_names, function(s){
  run_models(dlist_under[[s]], s)
})
names(results_under) <- under_names

# Mediation analysis
# results_*에 대해 일괄 수행 ----
library(mediation)

sims_n <- 100  # 필요시 조정

med_once <- function(res, sido, segment, sims = sims_n) {
  tryCatch({
    mo <- mediate(res$med.fit, res$out.fit,
                  treat = "tem", mediator = "PM25_mean",
                  boot = TRUE, sims = sims)
    
    # 필요하면 플롯 저장 (주석 해제)
    # pdf(sprintf("med_%s_%s.pdf", sido, segment)); plot(mo); dev.off()
    
    data.frame(
      sido = sido, segment = segment,
      ACME = mo$d0, ADE = mo$z0, TE = mo$tau.coef, PM = mo$n0,
      p_ACME = mo$d0.p, p_ADE = mo$z0.p, p_TE = mo$tau.p, p_PM = mo$n0.p,
      notes = "ok", row.names = NULL
    )
  }, error = function(e) {
    data.frame(
      sido = sido, segment = segment,
      ACME = NA, ADE = NA, TE = NA, PM = NA,
      p_ACME = NA, p_ADE = NA, p_TE = NA, p_PM = NA,
      notes = paste0("error: ", e$message), row.names = NULL
    )
  })
}

# over / under 각각 루프 실행
res_over_med <- do.call(rbind, lapply(names(results_over), function(s) {
  med_once(results_over[[s]], s, "over")
}))

res_under_med <- do.call(rbind, lapply(names(results_under), function(s) {
  med_once(results_under[[s]], s, "under")
}))

# 최종 테이블
med_results <- rbind(res_over_med, res_under_med)
med_results
readr::write_excel_csv(med_results, "med_results.csv")

# 4. SOBEL TEST TO ACQUIRE P-VALUE

#install.packages('bda')
library(bda)
mediation_results <- matrix(NA,  14,  1)
colnames(mediation_results) <- c("sobel p")
rownames(mediation_results)<-c ("seoulo","seoulu","busano","busanu","daeguo","daeguu","incheono","incheonu",
                                "gwangjuo","gwangjuu","daejeono","daejeonu","ulsano","ulsanu")
data_list <- list(seoulo,seoulu,busano,busanu, daeguo,daeguu,incheono,incheonu,
                  gwangjuo, gwangjuu, daejeono, daejeonu, ulsano, ulsanu)
# data_list에 있는 각 도시 데이터셋을 순회하면서 미디에이션 테스트 결과 계산
for (i in seq_along(data_list)) {
  city_data <- data_list[[i]]
  # mediation.test 실행
  mediation_result <- mediation.test(city_data$pm25, city_data$meanT_lag01, city_data$nonacc_tot)
  # 결과 리스트에 추가
  mediation_results[[i]] <- mediation_result$Sobel[2]
}
mediation_results<-data.frame(mediation_results)
mediation_results<-round(mediation_results,digits = 2)
under_mediation_results <- mediation_results[seq(2, nrow(mediation_results), 2), ]
over_mediation_results  <- mediation_results[seq(1, nrow(mediation_results), 2), ]
under_mediation_results
over_mediation_results

# 5. CONDUCT META-ANALYSIS TO OBTAIN POOLED ESTIMATES

# 6. ESTIMATE ATTRIBUTABLE FRACTION (AF)