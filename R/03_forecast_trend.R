# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Trend-based resistance forecast 2025-2028: per-pair weighted logit regression with prediction intervals; rolling-origin backtest against naive persistence and an ARIMA + mixed-model ensemble comparator.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages(library(tidyverse))

root <- getwd()
if (!file.exists(file.path(root, "R", "amrsn_data_shared.R"))) stop("Run from the repository root; R/amrsn_data_shared.R not found.")
source(file.path(root, "R", "amrsn_data_shared.R"))
OUT <- file.path(root, "output", "forecast")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

EPS     <- 5e-4
lgt     <- function(p) { pp <- pmax(EPS, pmin(1 - EPS, p / 100)); log(pp / (1 - pp)) }
inv     <- function(x) 100 / (1 + exp(-x))
CARB    <- c("Imipenem", "Meropenem")
FUTURE  <- 2025:2028
ORIGINS <- 2021:2023
has_ens <- requireNamespace("forecast", quietly = TRUE) && requireNamespace("lme4", quietly = TRUE)

fl <- full_long %>% mutate(pair_id = paste(organism, drug, sep = " | ")) %>% arrange(organism, drug, year)

fit_trend <- function(d) stats::lm(lgt(res_pct) ~ year, data = d, weights = tested)

predict_trend <- function(m, years, w_new) {
  nd  <- data.frame(year = years); w <- rep(w_new, length(years))
  p95 <- stats::predict(m, nd, interval = "prediction", level = 0.95, weights = w)
  p80 <- stats::predict(m, nd, interval = "prediction", level = 0.80, weights = w)
  tibble(year = years,
         point = inv(p95[, "fit"]),
         lo80 = inv(p80[, "lwr"]), hi80 = inv(p80[, "upr"]),
         lo95 = inv(p95[, "lwr"]), hi95 = inv(p95[, "upr"]))
}

ensemble_forecast <- function(tr, fut) {
  h <- length(fut)
  ab <- tr %>% group_by(organism, drug, pair_id) %>% group_modify(function(.x, .y) {
    yv <- .x$res_pct
    if (length(unique(yv)) == 1L || stats::sd(yv) < 1e-6) return(tibble(year = fut, pred_arima = rep(yv[length(yv)], h)))
    fit <- tryCatch(forecast::auto.arima(stats::ts(yv, start = min(.x$year), frequency = 1), stepwise = TRUE, approximation = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) {
      m <- stats::lm(res_pct ~ year, data = .x, weights = tested)
      return(tibble(year = fut, pred_arima = as.numeric(stats::predict(m, data.frame(year = fut)))))
    }
    tibble(year = fut, pred_arima = as.numeric(forecast::forecast(fit, h = h)$mean))
  }) %>% ungroup()
  fm <- tryCatch(suppressMessages(lme4::lmer(res_pct ~ year + (year | pair_id), data = tr, weights = tr$tested,
                            control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))), error = function(e) NULL)
  nd <- tidyr::crossing(tr %>% distinct(organism, drug, pair_id), year = fut)
  nd$pred_lmer <- if (!is.null(fm)) as.numeric(stats::predict(fm, newdata = nd, re.form = NULL)) else {
    fg <- stats::lm(res_pct ~ year + pair_id, data = tr, weights = tested); as.numeric(stats::predict(fg, newdata = nd)) }
  nd %>% left_join(ab, by = c("organism", "drug", "pair_id", "year")) %>%
    mutate(ensemble = (pmax(0, pmin(100, pred_arima)) + pmax(0, pmin(100, pred_lmer))) / 2) %>%
    select(organism, drug, pair_id, year, ensemble)
}

fc <- fl %>% group_by(organism, drug, pair_id) %>% group_modify(function(.x, .y) {
  m <- fit_trend(.x)
  predict_trend(m, FUTURE, mean(.x$tested)) %>%
    mutate(slope_logit_per_year = unname(stats::coef(m)["year"]), resid_sigma = summary(m)$sigma,
           n_years = nrow(.x), pi_denominator = round(mean(.x$tested)))
}) %>% ungroup()
write_csv(fc, file.path(OUT, "table_forecast_trend_2025_2028.csv"))

s8 <- fc %>% filter(drug %in% CARB) %>%
  transmute(organism, drug, year, point = round(point, 1),
            pi80 = sprintf("%.1f-%.1f", lo80, hi80), pi95 = sprintf("%.1f-%.1f", lo95, hi95))
write_csv(s8, file.path(OUT, "table_S8_forecast_carbapenems.csv"))

bt <- list()
for (O in ORIGINS) {
  tr <- fl %>% filter(year <= O); te <- fl %>% filter(year > O)
  fut <- sort(unique(te$year))
  trend <- tr %>% group_by(organism, drug, pair_id) %>% group_modify(function(.x, .y) {
    m <- fit_trend(.x)
    predict_trend(m, fut, mean(.x$tested)) %>% mutate(naive = .x$res_pct[which.max(.x$year)])
  }) %>% ungroup()
  ens <- if (has_ens) ensemble_forecast(tr, fut) else tibble(organism = character(), drug = character(), pair_id = character(), year = integer(), ensemble = numeric())
  bt[[length(bt) + 1]] <- trend %>%
    left_join(ens, by = c("organism", "drug", "pair_id", "year")) %>%
    inner_join(te %>% select(pair_id, year, actual = res_pct), by = c("pair_id", "year")) %>%
    mutate(origin = O, horizon = year - O)
}
bt <- bind_rows(bt) %>%
  mutate(stratum = if_else(drug %in% CARB, "carbapenem", "other"),
         err_trend = actual - point, err_naive = actual - naive, err_ensemble = actual - ensemble,
         in80 = actual >= lo80 & actual <= hi80, in95 = actual >= lo95 & actual <= hi95)
write_csv(bt, file.path(OUT, "table_forecast_backtest_rows.csv"))

summ <- function(d, lab) tibble(set = lab, n = nrow(d),
  MAE_trend = mean(abs(d$err_trend)), MAE_naive = mean(abs(d$err_naive)), MAE_ensemble = mean(abs(d$err_ensemble), na.rm = TRUE),
  bias_trend = mean(d$err_trend), coverage80 = mean(d$in80), coverage95 = mean(d$in95))
bt_summary <- bind_rows(
  summ(bt, "all pairs"),
  summ(filter(bt, stratum == "carbapenem"), "carbapenems"),
  bind_rows(lapply(sort(unique(bt$horizon)), function(h) summ(filter(bt, horizon == h), paste0("all pairs, horizon ", h)))),
  bind_rows(lapply(sort(unique(bt$horizon)), function(h) summ(filter(bt, stratum == "carbapenem", horizon == h), paste0("carbapenems, horizon ", h)))))
write_csv(bt_summary, file.path(OUT, "table_forecast_backtest_summary.csv"))

obs <- fl %>% filter(drug %in% CARB) %>% select(organism, drug, year, res_pct)
fitted_line <- fl %>% filter(drug %in% CARB) %>% group_by(organism, drug) %>% group_modify(function(.x, .y) {
  m <- fit_trend(.x); tibble(year = 2017:2024, fitted = inv(stats::predict(m, data.frame(year = 2017:2024)))) }) %>% ungroup()
fcc  <- fc %>% filter(drug %in% CARB)
cols <- c(Imipenem = "#2166AC", Meropenem = "#B2182B")
fig5 <- ggplot() +
  geom_ribbon(data = fcc, aes(x = year, ymin = lo95, ymax = hi95, fill = drug), alpha = 0.12) +
  geom_ribbon(data = fcc, aes(x = year, ymin = lo80, ymax = hi80, fill = drug), alpha = 0.28) +
  geom_line(data = fitted_line, aes(year, fitted, colour = drug), linewidth = 0.5, linetype = "dotted") +
  geom_line(data = obs, aes(year, res_pct, colour = drug), linewidth = 0.8) +
  geom_point(data = obs, aes(year, res_pct, colour = drug), size = 1.8) +
  geom_line(data = fcc, aes(year, point, colour = drug), linewidth = 0.9) +
  geom_point(data = fcc, aes(year, point, colour = drug), shape = 17, size = 2.2) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", colour = "grey55") +
  facet_wrap(~organism, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = seq(2017, 2028, 2)) +
  labs(x = NULL, y = "Resistance (%)", colour = NULL, fill = NULL,
       title = "Trend-based carbapenem resistance forecast, 2025-2028",
       subtitle = "Observed 2017-2024 (circles; dotted = fitted weighted logit trend) and projection (triangles); shading = 80% (dark) and 95% (light) prediction intervals; dashed line = last observed year") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "top", strip.text = element_text(face = "bold.italic", size = 11),
        plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 8.5))
ggsave(file.path(OUT, "fig5_forecast_trend.png"), fig5, width = 12, height = 10, dpi = 300)
ggsave(file.path(OUT, "fig5_forecast_trend.pdf"), fig5, width = 12, height = 10)

kp <- fc %>% filter(organism == "K. pneumoniae", drug == "Meropenem")
a  <- bt_summary %>% filter(set == "all pairs"); cb <- bt_summary %>% filter(set == "carbapenems")
writeLines(c(
  "# Trend-based forecast: run summary",
  sprintf("K. pneumoniae meropenem, weighted logit trend: %s", paste(sprintf("%d %.1f%% (80%% PI %.1f-%.1f; 95%% PI %.1f-%.1f)", kp$year, kp$point, kp$lo80, kp$hi80, kp$lo95, kp$hi95), collapse = "; ")),
  sprintf("Backtest, all pairs (n=%d): MAE trend %.2f, naive %.2f, ensemble %.2f pp; bias trend %.2f; PI coverage 80%%: %.2f, 95%%: %.2f", a$n, a$MAE_trend, a$MAE_naive, a$MAE_ensemble, a$bias_trend, a$coverage80, a$coverage95),
  sprintf("Backtest, carbapenems (n=%d): MAE trend %.2f, naive %.2f, ensemble %.2f pp; bias trend %.2f; PI coverage 80%%: %.2f, 95%%: %.2f", cb$n, cb$MAE_trend, cb$MAE_naive, cb$MAE_ensemble, cb$bias_trend, cb$coverage80, cb$coverage95),
  if (has_ens) "Ensemble comparator: auto.arima + lme4 random-slope model, equally weighted." else "Ensemble comparator skipped (forecast/lme4 not installed)."),
  file.path(OUT, "RUN_SUMMARY.md"))
cat(readLines(file.path(OUT, "RUN_SUMMARY.md")), sep = "\n")
