library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)

# 1) 데이터 불러오기
data <- readRDS("C:/Users/cmc/Desktop/과제/1.기후변화_자가면역/NEDIS/1_data/1_rds/20251002_casecross/NEDIS_casecross_primary_20251002.rds")

# 2) adm_date를 Date로 변환 + 사건 ID 부여
cases <- data %>%
  mutate(
    adm_date = as.Date(adm_date),
    case_id = row_number()   # 고유 사건 식별자
  )

# 3) 동일 월·요일 컨트롤 날짜 리스트 생성
same_wday_controls <- function(x_date) {
  m_start <- floor_date(x_date, "month")
  m_end   <- ceiling_date(x_date, "month") - days(1)
  cand <- seq(m_start, m_end, by = "day")
  cand[wday(cand, week_start = 1) == wday(x_date, week_start = 1) & cand != x_date]
}

cases_ctrl <- cases %>%
  mutate(
    controls = map(adm_date, same_wday_controls),
    n_ctrl = lengths(controls)
  )

# 4) AID 변수명 자동 탐색
aid_vars <- grep("^AID[0-9]+$", names(cases), value = TRUE)

# 5) 케이스 데이터 (case = 1)
long_case <- cases_ctrl %>%
  mutate(case = 1L) %>%
  select(case_id, adm_date, ptmigucd_2023, case, all_of(aid_vars))

# 6) 컨트롤 데이터 (case = 0) + AID 0으로 세팅
long_ctrl <- cases_ctrl %>%
  select(case_id, ptmigucd_2023, controls) %>%
  unnest_longer(controls, values_to = "adm_date") %>%
  mutate(case = 0L)

# 컨트롤용 AID열 0으로 생성
aid_zero <- as_tibble(matrix(0, nrow = nrow(long_ctrl), ncol = length(aid_vars)))
names(aid_zero) <- aid_vars
long_ctrl <- bind_cols(long_ctrl, aid_zero)

# 7) 결합
long_dat <- bind_rows(long_case, long_ctrl) %>%
  mutate(
    year = year(adm_date),
    month = month(adm_date),
    wdayN = wday(adm_date, week_start = 1)
  ) %>%
  group_by(case_id, year, month, wdayN) %>%
  mutate(set_id = cur_group_id()) %>%
  ungroup()

# 8) 확인
long_dat %>%
  count(set_id, case) %>%
  tidyr::pivot_wider(names_from = case, values_from = n,
                     names_prefix = "case_", values_fill = 0) %>%
  print(n = 10)

saveRDS(long_dat, file = "C:/Users/cmc/Desktop/과제/1.기후변화_자가면역/NEDIS/1_data/1_rds/20251010_casecross/NEDIS_casecross_primary_20251010.rds")
