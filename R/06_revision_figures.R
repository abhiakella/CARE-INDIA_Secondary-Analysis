# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Journal-style main figures (no in-figure titles, lowercase panel tags, no panel borders; TIFF 300 dpi >= 107 mm) and
# the year-by-year antibiogram supplementary figure. Figure 1 shows all five organisms (E. coli included).
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
root <- getwd()
if (!file.exists(file.path(root, "R", "amrsn_data_shared.R"))) stop("Run from the repository root; R/amrsn_data_shared.R not found.")
source(file.path(root, "R", "amrsn_data_shared.R"))
OUT <- file.path(root, "output", "figures_revision"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

th <- theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        legend.position = "bottom", legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        strip.text = element_text(face = "bold.italic", size = 9.5), plot.tag = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 10), axis.text = element_text(size = 9))
org_cols <- c("K. pneumoniae" = "#E41A1C", "A. baumannii" = "#377EB8", "P. aeruginosa" = "#4DAF4A",
              "Enterobacter spp." = "#FF7F00", "E. coli" = "#984EA3", "Enterobacterales" = "#A65628")
org_ord <- c("E. coli", "K. pneumoniae", "A. baumannii", "P. aeruginosa", "Enterobacter spp.")
save_fig <- function(p, name, w, h) {   # w, h in mm; TIFF for submission, PNG for embedding
  ggsave(file.path(OUT, paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 300, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")),  p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".jpg")),  p, width = w, height = h, units = "mm", dpi = 300, device = "jpeg", quality = 95, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".pdf")),  p, width = w, height = h, units = "mm", device = grDevices::cairo_pdf, bg = "white")
}

# ---- Figure 1: meropenem resistance trends, five organisms; S1: imipenem; S2: faceted both carbapenems ----
trend_plot <- function(abx) {
  amrsn_ee_long %>% filter(antibiotic == abx) %>% mutate(organism = factor(organism, levels = org_ord)) %>%
    ggplot(aes(year, resistance_pct, colour = organism, shape = organism)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
    scale_colour_manual(values = org_cols, name = NULL) + scale_shape_discrete(name = NULL) +
    scale_x_continuous(breaks = 2017:2024) + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(x = "Year", y = paste(abx, "resistance (%)")) + th + theme(legend.text = element_text(face = "italic"))
}
save_fig(trend_plot("Meropenem"), "Figure1_meropenem_trends", 170, 95)
save_fig(trend_plot("Imipenem"),  "FigureS1_imipenem_trends", 170, 95)
s2 <- amrsn_ee_long %>% mutate(organism = factor(organism, levels = org_ord)) %>%
  ggplot(aes(year, resistance_pct, colour = antibiotic, linetype = antibiotic)) + geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  facet_wrap(~organism, ncol = 3) + scale_colour_manual(values = c(Imipenem = "#D95F02", Meropenem = "#1B9E77"), name = NULL) +
  scale_linetype_discrete(name = NULL) + scale_x_continuous(breaks = seq(2017, 2024, 2)) + scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Year", y = "Resistance (%)") + th
save_fig(s2, "FigureS2_carbapenems_by_organism", 170, 110)

# ---- Figure 2: 2024 cumulative antibiogram heatmap; Figure 3: therapeutic tiers (80 / 60 / 10 boundaries) ----
drug_ord <- c("Piperacillin-tazobactam", "Cefotaxime", "Ceftazidime", "Cefepime", "Ertapenem", "Imipenem", "Meropenem",
              "Amikacin", "Gentamicin", "Tobramycin", "Ciprofloxacin", "Levofloxacin", "Colistin", "Minocycline")
d24 <- full_long %>% filter(year == 2024) %>% select(organism, drug, tested, susc_pct) %>%
  tidyr::complete(organism = org_ord, drug = drug_ord) %>%                      # untested agent-organism pairs become explicit NA cells
  mutate(drug = factor(drug, levels = rev(drug_ord)), organism = factor(organism, levels = org_ord), tested_flag = !is.na(susc_pct))
NT_FILL <- "#C8C8C8"
f2 <- ggplot(d24, aes(organism, drug, fill = susc_pct)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = d24 %>% filter(tested_flag), aes(label = round(susc_pct)), size = 2.8) +
  geom_text(data = d24 %>% filter(!tested_flag), label = "NT", size = 2.4, colour = "grey35") +
  scale_fill_gradientn(colours = c("#D73027", "#FC8D59", "#FEE090", "#ABD9E9", "#4575B4"), values = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 100), name = "Susceptible (%)",
                       na.value = NT_FILL) +
  labs(caption = "Grey cells (NT): agent not tested for that organism") +
  labs(x = NULL, y = NULL) + th + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1, face = "italic"))
save_fig(f2, "Figure2_antibiogram_2024", 150, 120)

TIER_BREAKS <- c(-1, 10, 60, 80, 101); TIER_LABELS <- c("Inactive (<10%)", "Not recommended (10-59%)", "Conditional (60-79%)", "Preferred (>=80%)")
tier24 <- d24 %>% mutate(tier = cut(susc_pct, TIER_BREAKS, labels = TIER_LABELS),
                         tier = factor(ifelse(is.na(tier), "Not tested", as.character(tier)), levels = c(TIER_LABELS, "Not tested")))
write_csv(tier24 %>% filter(tested_flag) %>% select(organism, drug, tested, susc_pct, tier) %>% arrange(organism, drug), file.path(OUT, "table_tiers_2024_80_60.csv"))
f3 <- ggplot(tier24, aes(organism, drug, fill = tier)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = tier24 %>% filter(tested_flag), aes(label = round(susc_pct)), size = 2.6) +
  scale_fill_manual(values = c("#D73027", "#FC8D59", "#FEE090", "#4575B4", NT_FILL), name = "Tier", drop = FALSE) + guides(fill = guide_legend(nrow = 3, byrow = FALSE)) +
  labs(x = NULL, y = NULL) + th + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1, face = "italic"))
save_fig(f3, "Figure3_tier_map_2024", 150, 125)

# ---- Figure 4: carbapenemase genes, (a) ICMR-AMRSN 2024 panel and (b) curated Indian literature ----
gt <- read_csv(file.path(root, "data", "gene_trends.csv"), show_col_types = FALSE)
genes <- c("NDM", "OXA-23", "OXA-48", "KPC", "VIM", "IMP"); gene_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#8BC34A", "#F2C80F")
ga <- gt %>% filter(gene %in% genes, year == 2024, organism %in% c("E. coli", "K. pneumoniae", "A. baumannii", "P. aeruginosa")) %>%
  mutate(organism = factor(organism, levels = c("E. coli", "K. pneumoniae", "A. baumannii", "P. aeruginosa")), gene = factor(gene, levels = genes))
f4a <- ggplot(ga, aes(organism, prevalence_pct, fill = gene)) + geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  scale_fill_manual(values = setNames(gene_cols, genes), name = NULL, drop = FALSE) + guides(fill = guide_legend(nrow = 1)) + scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Detection frequency (%)") + th + theme(axis.text.x = element_text(face = "italic"), legend.position = "top")
gb <- read_csv(file.path(root, "data", "fig4_panel_b_curated_literature.csv"), show_col_types = FALSE) %>%
  mutate(label = factor(label, levels = rev(label)), txt = ifelse(is.na(lo), sprintf("%g%%", point), sprintf("%g%% [%g-%g]", point, lo, hi)))
f4b <- ggplot(gb, aes(point, label, fill = gene)) + geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", na.rm = TRUE) + geom_text(aes(x = pmax(point, hi, na.rm = TRUE) + 2, label = txt), hjust = 0, size = 2.8) +
  scale_fill_manual(values = setNames(gene_cols, genes), guide = "none") + scale_x_continuous(limits = c(0, 125), breaks = seq(0, 100, 20)) +
  labs(x = "Reported prevalence in India (%)", y = NULL) + th + theme(axis.text.y = element_text(face = "italic", size = 9))
f4 <- (f4a / f4b) + plot_layout(heights = c(1, 0.9)) + plot_annotation(tag_levels = "a")
save_fig(f4, "Figure4_carbapenemase_genes", 170, 170)

# ---- Figure 5: trend-based carbapenem forecast (needs output/forecast/table_forecast_trend_2025_2028.csv) ----
fcf <- file.path(root, "output", "forecast", "table_forecast_trend_2025_2028.csv")
if (file.exists(fcf)) {
  EPS <- 5e-4; lgt <- function(p) { pp <- pmax(EPS, pmin(1 - EPS, p / 100)); log(pp / (1 - pp)) }; inv <- function(x) 100 / (1 + exp(-x))
  CARB <- c("Imipenem", "Meropenem"); fc <- read_csv(fcf, show_col_types = FALSE) %>% filter(drug %in% CARB)
  obs <- full_long %>% filter(drug %in% CARB) %>% select(organism, drug, year, res_pct)
  fitted_line <- full_long %>% filter(drug %in% CARB) %>% group_by(organism, drug) %>% group_modify(function(.x, .y) {
    m <- stats::lm(lgt(res_pct) ~ year, data = .x, weights = tested); tibble(year = 2017:2024, fitted = inv(stats::predict(m, data.frame(year = 2017:2024)))) }) %>% ungroup()
  cols <- c(Imipenem = "#2166AC", Meropenem = "#B2182B")
  f5 <- ggplot() +
    geom_ribbon(data = fc, aes(year, ymin = lo95, ymax = hi95, fill = drug), alpha = 0.12) +
    geom_ribbon(data = fc, aes(year, ymin = lo80, ymax = hi80, fill = drug), alpha = 0.28) +
    geom_line(data = fitted_line, aes(year, fitted, colour = drug), linewidth = 0.45, linetype = "dotted") +
    geom_line(data = obs, aes(year, res_pct, colour = drug), linewidth = 0.7) + geom_point(data = obs, aes(year, res_pct, colour = drug), size = 1.5) +
    geom_line(data = fc, aes(year, point, colour = drug), linewidth = 0.8) + geom_point(data = fc, aes(year, point, colour = drug), shape = 17, size = 1.9) +
    geom_vline(xintercept = 2024.5, linetype = "dashed", colour = "grey55") +
    facet_wrap(~organism, scales = "free_y", ncol = 2) + scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
    scale_x_continuous(breaks = seq(2017, 2028, 2)) + labs(x = NULL, y = "Resistance (%)") + th
  save_fig(f5, "Figure5_forecast_carbapenems", 170, 150)
} else cat("Figure 5 skipped: run R/03_forecast_trend.R first\n")

# ---- Figure 6: meropenem susceptibility by infection syndrome, 2024, with 80% and 60% reference lines ----
syn <- read_csv(file.path(root, "data", "syndrome_key_antibiogram_2024.csv"), show_col_types = FALSE) %>% filter(drug == "Meropenem") %>%
  mutate(syndrome = factor(syndrome, levels = c("BSI", "UTI", "VAP"), labels = c("Bloodstream infection", "Urinary tract infection", "Ventilator-associated pneumonia")),
         organism = factor(organism, levels = c("K. pneumoniae", "E. coli", "Enterobacterales", "A. baumannii", "P. aeruginosa")))
f6 <- ggplot(syn, aes(syndrome, susc_pct, fill = organism)) + geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey40") + geom_hline(yintercept = 60, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 3.45, y = 82, label = "80%", size = 2.8, colour = "grey30", hjust = 1) + annotate("text", x = 3.45, y = 62, label = "60%", size = 2.8, colour = "grey30", hjust = 1) +
  scale_fill_manual(values = org_cols, name = NULL, labels = c(expression(italic("K. pneumoniae")), expression(italic("E. coli")), "Enterobacterales", expression(italic("A. baumannii")), expression(italic("P. aeruginosa")))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Meropenem susceptibility (%)") + th
save_fig(f6, "Figure6_syndrome_meropenem_2024", 170, 95)

# ---- Supplementary: year-by-year antibiogram heatmap, all 47 organism-drug pairs ----
hm <- full_long %>% select(organism, drug, year, susc_pct) %>% tidyr::complete(organism = org_ord, drug = drug_ord, year = 2017:2024) %>%
  mutate(drug = factor(drug, levels = rev(drug_ord)), organism = factor(organism, levels = org_ord), yr = factor(year), tested_flag = !is.na(susc_pct))
fs <- ggplot(hm, aes(yr, drug, fill = susc_pct)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(data = hm %>% filter(tested_flag), aes(label = round(susc_pct)), size = 2.1) +
  facet_wrap(~organism, ncol = 3) +
  scale_fill_gradientn(colours = c("#D73027", "#FC8D59", "#FEE090", "#ABD9E9", "#4575B4"), values = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 100), name = "Susceptible (%)",
                       na.value = NT_FILL) +
  labs(caption = "Grey cells: agent not tested for that organism") +
  scale_x_discrete(labels = function(x) paste0("'", substr(x, 3, 4))) +
  labs(x = "Year", y = NULL) + th + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7.5))
save_fig(fs, "FigureS_antibiogram_by_year", 190, 215)
cat("Revision figures written to output/figures_revision/\n")
