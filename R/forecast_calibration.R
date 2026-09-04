# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Forecast calibration: rolling-origin backtest + bias-correction + logit split-conformal + Jeffreys floor for the 2025-2028 projection.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressMessages({ library(tidyverse); library(forecast); library(lme4) })

source("R/amrsn_data_shared.R")
OUT <- "output/ml"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ORIG_FC <- "output/ml/table_forecast_resistance_2025_2028.csv"

Z80 <- qnorm(0.90); Z95 <- qnorm(0.975)
EPS <- 5e-4
lg  <- function(p) { pp <- pmax(EPS, pmin(1 - EPS, p / 100)); log(pp / (1 - pp)) }
inv <- function(x) 100 / (1 + exp(-x))
carb_drugs <- c("Imipenem", "Meropenem")
strat_of   <- function(drug) ifelse(drug %in% carb_drugs, "carbapenem", "other")
conf_q <- function(resid, level) {
  resid <- resid[is.finite(resid)]; n <- length(resid)
  if (n == 0) return(NA_real_)
  if (n < 5)  return(stats::sd(resid) * qnorm(1 - (1 - level) / 2))
  sort(abs(resid))[min(ceiling((n + 1) * level), n)]
}
jeffreys_int <- function(p_pct, n, level) {
  a <- 0.5 + p_pct / 100 * n; b <- 0.5 + (1 - p_pct / 100) * n
  c(lo = 100 * qbeta((1 - level) / 2, a, b), hi = 100 * qbeta(1 - (1 - level) / 2, a, b))
}
CHECKS <- list()
chk <- function(name, ok, detail = "") {
  CHECKS[[length(CHECKS) + 1]] <<- tibble(check = name, status = ifelse(ok, "PASS", "FAIL"), detail = detail)
  cat(sprintf("  [%s] %s %s\n", ifelse(ok, "PASS", "FAIL"), name, detail))
}
near <- function(a, b, tol) is.finite(a) && abs(a - b) <= tol

fl <- full_long %>% mutate(pair_id = paste(organism, drug, sep = " | "))
future_years <- 2025L:2028L

arima_pair <- function(yv, w, yrs, fut) {
  h <- length(fut)
  if (length(unique(yv)) == 1L || stats::sd(yv) < 1e-6)
    return(tibble(pred_arima = rep(yv[length(yv)], h), se_arima = 5 / Z95))
  ts_y <- stats::ts(yv, start = min(yrs), frequency = 1)
  fit  <- tryCatch(forecast::auto.arima(ts_y, stepwise = TRUE, approximation = TRUE), error = function(e) NULL)
  if (is.null(fit)) {
    m  <- stats::lm(yv ~ yrs, weights = w)
    pr <- predict(m, newdata = data.frame(yrs = fut), interval = "prediction", level = 0.95)
    return(tibble(pred_arima = pr[, "fit"], se_arima = (pr[, "upr"] - pr[, "lwr"]) / (2 * Z95)))
  }
  fc <- forecast::forecast(fit, h = h, level = 95)
  tibble(pred_arima = as.numeric(fc$mean),
         se_arima  = (as.numeric(fc$upper[, 1]) - as.numeric(fc$lower[, 1])) / (2 * Z95))
}

cat("\n=========================== SECTION 1: FIXED-PI FORECAST ===========================\n")
ab_full <- fl %>% group_by(organism, drug, pair_id) %>%
  group_modify(function(.x, .y) { .x <- arrange(.x, year)
    out <- arima_pair(.x$res_pct, .x$tested, .x$year, future_years); out$year <- future_years; out }) %>% ungroup()
fm_full <- tryCatch(lme4::lmer(res_pct ~ year + (year | pair_id), data = fl, weights = fl$tested,
                    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), error = function(e) NULL)
pm <- fl %>% distinct(organism, drug, pair_id); nd <- tidyr::crossing(pm, year = future_years)
nd$pred_lmer <- if (!is.null(fm_full)) pmax(0, pmin(100, as.numeric(predict(fm_full, newdata = nd, re.form = NULL)))) else {
  fg <- lm(res_pct ~ year + pair_id, data = fl, weights = tested); pmax(0, pmin(100, as.numeric(predict(fg, newdata = nd)))) }

forecast_tbl <- nd %>% left_join(ab_full, by = c("organism", "drug", "pair_id", "year")) %>%
  mutate(pred_arima = pmax(0, pmin(100, pred_arima)),
         pred_ensemble = (pred_arima + pred_lmer) / 2,
         sd_ens = sqrt(se_arima^2 + (pred_arima - pred_lmer)^2 / 4),
         lo80 = pmax(0, pmin(100, pred_ensemble - Z80 * sd_ens)), hi80 = pmax(0, pmin(100, pred_ensemble + Z80 * sd_ens)),
         lo95 = pmax(0, pmin(100, pred_ensemble - Z95 * sd_ens)), hi95 = pmax(0, pmin(100, pred_ensemble + Z95 * sd_ens)))
write_csv(forecast_tbl %>% select(organism, drug, pair_id, year, pred_lmer, pred_arima,
                                  lo80, hi80, lo95, hi95, pred_ensemble, sd_ens, se_arima),
          file.path(OUT, "table_forecast_resistance_2025_2028_fixedPI.csv"))
if (file.exists(ORIG_FC)) {
  orig <- read_csv(ORIG_FC, show_col_types = FALSE) %>% select(organism, drug, year, pe_old = pred_ensemble)
  cmp  <- forecast_tbl %>% left_join(orig, by = c("organism", "drug", "year"))
  chk("fixedPI points reproduce canonical forecast", all(abs(cmp$pred_ensemble - cmp$pe_old) < 1e-6, na.rm = TRUE),
      sprintf("max|Δ|=%.2e", max(abs(cmp$pred_ensemble - cmp$pe_old), na.rm = TRUE)))
}
chk("all ensemble points inside their fixed 95% PI",
    all(forecast_tbl$pred_ensemble >= forecast_tbl$lo95 - 1e-9 & forecast_tbl$pred_ensemble <= forecast_tbl$hi95 + 1e-9))

cat("\n=========================== SECTION 2: ROLLING-ORIGIN BACKTEST ===========================\n")
origins <- c(2021L, 2022L, 2023L); rows <- list()
for (O in origins) {
  tr  <- fl %>% filter(year <= O); fut <- (O + 1L):2024L
  ab <- tr %>% group_by(organism, drug, pair_id) %>%
    group_modify(function(.x, .y) { .x <- arrange(.x, year)
      out <- arima_pair(.x$res_pct, .x$tested, .x$year, fut); out$year <- fut; out }) %>% ungroup()
  fm <- tryCatch(lme4::lmer(res_pct ~ year + (year | pair_id), data = tr, weights = tr$tested,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), error = function(e) NULL)
  pmO <- tr %>% distinct(organism, drug, pair_id); ndO <- tidyr::crossing(pmO, year = fut)
  ndO$pred_lmer <- if (!is.null(fm)) pmax(0, pmin(100, as.numeric(predict(fm, newdata = ndO, re.form = NULL)))) else {
    fg <- lm(res_pct ~ year + pair_id, data = tr, weights = tested); pmax(0, pmin(100, as.numeric(predict(fg, newdata = ndO)))) }
  last_val <- tr %>% group_by(pair_id) %>% slice_max(year, n = 1, with_ties = FALSE) %>% transmute(pair_id, naive = res_pct) %>% ungroup()
  lin <- tr %>% group_by(organism, drug, pair_id) %>%
    group_modify(function(.x, .y) { m <- tryCatch(stats::lm(res_pct ~ year, data = .x, weights = tested), error = function(e) NULL)
      tibble(year = fut, linear = if (is.null(m)) NA_real_ else pmax(0, pmin(100, as.numeric(predict(m, newdata = data.frame(year = fut)))))) }) %>% ungroup()
  actual <- fl %>% filter(year %in% fut) %>% transmute(pair_id, year, actual = res_pct)
  df <- ab %>% left_join(ndO, by = c("organism", "drug", "pair_id", "year")) %>%
    left_join(last_val, by = "pair_id") %>% left_join(lin, by = c("organism", "drug", "pair_id", "year")) %>%
    left_join(actual, by = c("pair_id", "year")) %>%
    mutate(origin = O, horizon = year - O, pred_arima = pmax(0, pmin(100, pred_arima)),
           ensemble = (pred_arima + pred_lmer) / 2, sd_ens = sqrt(se_arima^2 + (pred_arima - pred_lmer)^2 / 4)) %>%
    filter(!is.na(actual))
  rows[[length(rows) + 1]] <- df
}
bt <- bind_rows(rows) %>% mutate(stratum = strat_of(drug),
                                 error_pp = actual - ensemble, e_logit = lg(actual) - lg(ensemble))
write_csv(bt %>% select(pair_id, organism, drug, origin, horizon, year, ensemble, naive, linear,
                        actual, error_pp, e_logit, sd_ens, stratum), file.path(OUT, "forecast_backtest_rows.csv"))
bt_sum <- function(d, lab) tibble(set = lab, n = nrow(d),
  ens_MAE = mean(abs(d$ensemble - d$actual)), naive_MAE = mean(abs(d$naive - d$actual)),
  lin_MAE = mean(abs(d$linear - d$actual), na.rm = TRUE), bias = mean(d$ensemble - d$actual),
  cov80 = mean(d$actual >= d$ensemble - Z80 * d$sd_ens & d$actual <= d$ensemble + Z80 * d$sd_ens),
  cov95 = mean(d$actual >= d$ensemble - Z95 * d$sd_ens & d$actual <= d$ensemble + Z95 * d$sd_ens))
bt_metrics <- bind_rows(bt_sum(bt, "ALL"), bt_sum(bt %>% filter(stratum == "carbapenem"), "Carbapenems"))
write_csv(bt_metrics, file.path(OUT, "backtest_metrics.csv"))
print(as.data.frame(bt_metrics %>% mutate(across(where(is.numeric), ~round(., 2)))), row.names = FALSE)
bA <- bt_metrics %>% filter(set == "ALL"); bC <- bt_metrics %>% filter(set == "Carbapenems")
chk("backtest ALL ens_MAE~4.99 & naive~4.71", near(bA$ens_MAE, 4.99, 0.2) && near(bA$naive_MAE, 4.71, 0.2),
    sprintf("ens=%.2f naive=%.2f", bA$ens_MAE, bA$naive_MAE))
chk("backtest ALL biased low (~-1.87) & PIs overconfident", near(bA$bias, -1.87, 0.3) && bA$cov95 < 0.85,
    sprintf("bias=%.2f cov80=%.2f cov95=%.2f", bA$bias, bA$cov80, bA$cov95))
chk("backtest carbapenem linear beats ensemble", bC$lin_MAE < bC$ens_MAE,
    sprintf("lin=%.2f ens=%.2f", bC$lin_MAE, bC$ens_MAE))

loo_eval <- function(scale) {
  res <- list()
  for (Ostar in origins) {
    cal <- bt %>% filter(origin != Ostar); tst <- bt %>% filter(origin == Ostar)
    key <- if (scale == "logit") "e_logit" else "error_pp"
    cp <- cal %>% group_by(stratum) %>%
      summarise(bias = mean(.data[[key]]), q80 = conf_q(.data[[key]] - mean(.data[[key]]), 0.80),
                q95 = conf_q(.data[[key]] - mean(.data[[key]]), 0.95), .groups = "drop")
    e <- tst %>% left_join(cp, by = "stratum")
    if (scale == "logit") {
      e <- e %>% mutate(center = lg(ensemble) + bias, point_new = inv(center), r = lg(actual) - center,
                        in80 = abs(r) <= q80, in95 = abs(r) <= q95, ae = abs(point_new - actual))
    } else {
      e <- e %>% mutate(point_new = pmax(0, pmin(100, ensemble + bias)),
                        in80 = abs(actual - point_new) <= q80, in95 = abs(actual - point_new) <= q95, ae = abs(actual - point_new))
    }
    e <- e %>% mutate(old_in80 = actual >= ensemble - Z80 * sd_ens & actual <= ensemble + Z80 * sd_ens,
                      old_in95 = actual >= ensemble - Z95 * sd_ens & actual <= ensemble + Z95 * sd_ens,
                      old_ae = abs(ensemble - actual))
    res[[length(res) + 1]] <- e
  }
  bind_rows(res)
}
loo_tab <- function(ev, scale) {
  f <- function(d, lab) tibble(scale = scale, set = lab, n = nrow(d),
    old_MAE = round(mean(d$old_ae), 2), new_MAE = round(mean(d$ae), 2),
    old_cov80 = round(mean(d$old_in80), 2), new_cov80 = round(mean(d$in80), 2),
    old_cov95 = round(mean(d$old_in95), 2), new_cov95 = round(mean(d$in95), 2))
  bind_rows(f(ev, "ALL pairs"), f(ev %>% filter(stratum == "carbapenem"), "Carbapenems"))
}

calibrate_forecast <- function(scale) {
  key <- if (scale == "logit") "e_logit" else "error_pp"
  cal_sh <- bt %>% group_by(stratum, horizon) %>%
    summarise(bias = mean(.data[[key]]), q80 = conf_q(.data[[key]] - mean(.data[[key]]), 0.80),
              q95 = conf_q(.data[[key]] - mean(.data[[key]]), 0.95), n = n(), .groups = "drop")
  ext <- cal_sh %>% group_by(stratum) %>%
    summarise(bias = bias[horizon == 3],
              q80 = q80[horizon == 3] * pmin(pmax(q80[horizon == 3] / q80[horizon == 2], 1), 1.6),
              q95 = q95[horizon == 3] * pmin(pmax(q95[horizon == 3] / q95[horizon == 2], 1), 1.6),
              n = 0L, .groups = "drop") %>% mutate(horizon = 4L) %>% select(stratum, horizon, bias, q80, q95, n)
  cal_all <- bind_rows(cal_sh, ext) %>% arrange(stratum, horizon)
  latest_n <- fl %>% group_by(pair_id) %>% slice_max(year, n = 1, with_ties = FALSE) %>% ungroup() %>% transmute(pair_id, n_proj = tested)
  base <- forecast_tbl %>% mutate(stratum = strat_of(drug), horizon = year - 2024L) %>%
    left_join(cal_all, by = c("stratum", "horizon")) %>% left_join(latest_n, by = "pair_id")
  if (scale == "logit") {
    out <- base %>% rowwise() %>%
      mutate(center = lg(pred_ensemble) + bias, point_cal = inv(center),
             lo80_c = inv(center - q80), hi80_c = inv(center + q80),
             lo95_c = inv(center - q95), hi95_c = inv(center + q95),
             j_lo = jeffreys_int(point_cal, n_proj, 0.95)["lo"], j_hi = jeffreys_int(point_cal, n_proj, 0.95)["hi"],
             lo95_c = min(lo95_c, j_lo), hi95_c = max(hi95_c, j_hi),
             floor_binds = (j_lo < lo95_c) | (j_hi > hi95_c)) %>% ungroup()
  } else {
    out <- base %>% rowwise() %>%
      mutate(point_cal = pmax(0, pmin(100, pred_ensemble + bias)),
             half80 = max(q80, (jeffreys_int(point_cal, n_proj, 0.80)["hi"] - jeffreys_int(point_cal, n_proj, 0.80)["lo"]) / 2),
             half95 = max(q95, (jeffreys_int(point_cal, n_proj, 0.95)["hi"] - jeffreys_int(point_cal, n_proj, 0.95)["lo"]) / 2),
             lo80_c = pmax(0, point_cal - half80), hi80_c = pmin(100, point_cal + half80),
             lo95_c = pmax(0, point_cal - half95), hi95_c = pmin(100, point_cal + half95),
             floor_binds = half95 > q95) %>% ungroup()
  }
  out %>% select(organism, drug, pair_id, year, horizon, stratum,
                 pred_ensemble_old = pred_ensemble, lo95_old = lo95, hi95_old = hi95,
                 point_cal, lo80_cal = lo80_c, hi80_cal = hi80_c, lo95_cal = lo95_c, hi95_cal = hi95_c, floor_binds)
}

cat("\n=========================== SECTION 3: TIER 1 (RAW SCALE) ===========================\n")
ev_raw <- loo_eval("raw"); loo_raw <- loo_tab(ev_raw, "raw")
print(as.data.frame(loo_raw), row.names = FALSE)
fc_raw <- calibrate_forecast("raw")
write_csv(fc_raw, file.path(OUT, "table_forecast_resistance_2025_2028_tier1cal.csv"))
write_csv(bt %>% select(pair_id, organism, drug, origin, horizon, year, ensemble, actual, error_pp, sd_ens, stratum),
          file.path(OUT, "tier1_backtest_errors.csv"))
rA <- loo_raw %>% filter(set == "ALL pairs"); rC <- loo_raw %>% filter(set == "Carbapenems")
chk("Tier1 raw LOO cov95 restored (~0.95)", near(rA$new_cov95, 0.95, 0.06), sprintf("cov95=%.2f (was %.2f)", rA$new_cov95, rA$old_cov95))
chk("Tier1 raw carbapenem MAE improved (<4.28 linear)", rC$new_MAE < 4.28, sprintf("new=%.2f", rC$new_MAE))

cat("\n=========================== SECTION 4: TIER 1 (LOGIT SCALE) ===========================\n")
ev_lgt <- loo_eval("logit"); loo_lgt <- loo_tab(ev_lgt, "logit")
print(as.data.frame(loo_lgt), row.names = FALSE)
fc_lgt <- calibrate_forecast("logit")
write_csv(fc_lgt, file.path(OUT, "table_forecast_resistance_2025_2028_tier1logit.csv"))
lA <- loo_lgt %>% filter(set == "ALL pairs"); lC <- loo_lgt %>% filter(set == "Carbapenems")
ab_max <- fc_lgt %>% filter(organism == "A. baumannii", drug %in% carb_drugs) %>% summarise(mx = max(point_cal)) %>% pull(mx)
chk("Tier1 logit LOO coverage near nominal (80/95)", near(lA$new_cov80, 0.80, 0.05) && near(lA$new_cov95, 0.95, 0.05),
    sprintf("cov80=%.2f cov95=%.2f", lA$new_cov80, lA$new_cov95))
chk("Tier1 logit carbapenem MAE improved (<3.60 raw)", lC$new_MAE < 3.60, sprintf("new=%.2f", lC$new_MAE))
chk("Tier1 logit FIXES A.baumannii saturation (max point <100)", ab_max < 99.9, sprintf("A.b max point=%.1f%%", ab_max))

poster_bg<-"#3C3C3C"; poster_text<-"#E8E8E8"; poster_grid<-"#525252"
theme_ml <- theme_minimal(base_size = 11) + theme(
  plot.background=element_rect(fill=poster_bg,colour=NA), panel.background=element_rect(fill=poster_bg,colour=NA),
  panel.grid.major=element_line(colour=poster_grid,linewidth=0.3), text=element_text(colour=poster_text),
  plot.title=element_text(colour="#FFFFFF",face="bold",hjust=0.5), plot.subtitle=element_text(hjust=0.5,colour="#AAAAAA",size=9),
  axis.text=element_text(colour=poster_text), axis.title=element_text(colour=poster_text), legend.position="none")
hist_fc <- fl %>% select(organism, drug, year, observed = res_pct)
figt <- fc_lgt %>% ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = lo80_cal, ymax = hi80_cal), alpha = 0.25, fill = "#48C9B0") +
  geom_line(aes(y = point_cal, colour = drug), linewidth = 0.8) + geom_point(aes(y = point_cal, colour = drug), size = 1.4) +
  geom_line(data = hist_fc, aes(x = year, y = observed, colour = drug), linewidth = 0.5, linetype = "dotted", alpha = 0.7) +
  facet_wrap(~organism, scales = "free_y", ncol = 2) + scale_x_continuous(breaks = seq(2017, 2028, 2)) +
  labs(title = "Tier 1 (logit) calibrated resistance forecast 2025-2028",
       subtitle = "Point = bias-corrected; shaded = 80% backtest-calibrated (conformal) PI; dotted = observed",
       x = NULL, y = "Resistance (%)") + theme_ml
ggsave(file.path(OUT, "fig_forecast_tier1logit_facets.png"), figt, width = 12, height = 10, dpi = 300)

cat("\n=========================== SECTION 5: CONSOLIDATED METRICS + SUMMARY ===========================\n")
metrics <- bind_rows(loo_raw, loo_lgt) %>% mutate(across(where(is.numeric), ~round(., 3)))
write_csv(metrics, file.path(OUT, "tier_metrics_summary.csv"))
checks_df <- bind_rows(CHECKS); write_csv(checks_df, file.path(OUT, "verification_checks.csv"))
all_pass <- all(checks_df$status == "PASS")

kp <- fc_lgt %>% filter(organism == "K. pneumoniae", drug == "Meropenem", year == 2028)
summ_md <- c(
  "# CARE-INDIA forecast calibration (self-verified)",
  "",
  sprintf("Generated by `forecast_calibration.R`. Overall gate: **%s** (%d/%d checks passed).",
          ifelse(all_pass, "PASS", "FAIL"), sum(checks_df$status == "PASS"), nrow(checks_df)),
  "",
  "## Outputs in this folder",
  "- `table_forecast_resistance_2025_2028_fixedPI.csv` — Section 1: ensemble point + mixture PI (points identical to canonical forecast).",
  "- `forecast_backtest_rows.csv`, `backtest_metrics.csv` — Section 2: rolling-origin backtest vs naive/linear.",
  "- `table_forecast_resistance_2025_2028_tier1cal.csv` — Section 3: Tier 1 raw-scale calibrated forecast.",
  "- `table_forecast_resistance_2025_2028_tier1logit.csv` — Section 4: **Tier 1 logit-scale (carry-forward)**.",
  "- `tier1_backtest_errors.csv`, `tier_metrics_summary.csv`, `verification_checks.csv`, `fig_forecast_tier1logit_facets.png`.",
  "",
  "## Leave-one-origin-out coverage & MAE (out-of-sample)",
  "| scale | set | old MAE | new MAE | old cov80 | new cov80 | old cov95 | new cov95 |",
  "|---|---|---|---|---|---|---|---|",
  sprintf("| raw | ALL | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |", rA$old_MAE, rA$new_MAE, rA$old_cov80, rA$new_cov80, rA$old_cov95, rA$new_cov95),
  sprintf("| raw | Carb | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |", rC$old_MAE, rC$new_MAE, rC$old_cov80, rC$new_cov80, rC$old_cov95, rC$new_cov95),
  sprintf("| logit | ALL | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |", lA$old_MAE, lA$new_MAE, lA$old_cov80, lA$new_cov80, lA$old_cov95, lA$new_cov95),
  sprintf("| logit | Carb | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |", lC$old_MAE, lC$new_MAE, lC$old_cov80, lC$new_cov80, lC$old_cov95, lC$new_cov95),
  "",
  sprintf("Headline K. pneumoniae meropenem 2028 (logit): **%.1f%%** (95%% PI %.1f-%.1f).", kp$point_cal, kp$lo95_cal, kp$hi95_cal),
  sprintf("A. baumannii carbapenem max calibrated point: **%.1f%%** (raw-scale pinned at 100.0 — artifact fixed).", ab_max),
  "",
  "## Verification checks",
  paste(sprintf("- [%s] %s %s", checks_df$status, checks_df$check, checks_df$detail), collapse = "\n"),
  "",
  "## Scope note",
  "This script covers the forecast calibration only (rolling-origin backtest, mixture-PI fix, Tier 1 raw + logit).",
  "Descriptive resistance prevalence and the ML classifier baselines are produced by the other deposit scripts.")
writeLines(summ_md, file.path(OUT, "RUN_SUMMARY.md"))

cat("\n================================ GATE ================================\n")
cat(sprintf("Checks passed: %d/%d  =>  %s\n", sum(checks_df$status == "PASS"), nrow(checks_df),
            ifelse(all_pass, "ALL CHECKS PASS", "SOME CHECKS FAILED")))
cat("Wrote all forecast-calibration outputs + RUN_SUMMARY.md\n")
if (!all_pass) quit(status = 1)
