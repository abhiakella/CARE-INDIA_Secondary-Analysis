# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Resistance archetypes: k-medoids (PAM) clustering of organisms on their per-drug resistance-slope trajectories.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(cluster)
})
set.seed(42)

root <- getwd()
dir.create("output", recursive = TRUE, showWarnings = FALSE)
csv_in <- "data/AMRSN_GN_ESKAPEE_Susceptibility_Trends.csv"

amrsn <- read_csv(csv_in, show_col_types = FALSE) |>
  filter(grepl("^All", Specimen_Type)) |>
  filter(!is.na(Tested_N), !is.na(Susceptible_N), Tested_N >= 30) |>
  mutate(
    organism   = Pathogen,
    antibiotic = Antibiotic,
    year       = as.integer(Year),
    res_pct    = as.numeric(Resistance_pct),
    tested_n   = as.integer(Tested_N)
  ) |>
  select(organism, antibiotic, year, res_pct, tested_n)

pair_counts <- amrsn |>
  count(organism, antibiotic, name = "n_years") |>
  filter(n_years >= 5)
amrsn <- amrsn |> semi_join(pair_counts, by = c("organism", "antibiotic"))

slope_tbl <- amrsn |>
  group_by(organism, antibiotic) |>
  arrange(year, .by_group = TRUE) |>
  summarise(
    slope_ppt_per_yr = tryCatch(
      coef(lm(res_pct ~ year, weights = tested_n))[["year"]],
      error = function(e) NA_real_
    ),
    n_years = n(),
    .groups = "drop"
  )
write_csv(slope_tbl, "output/table_organism_drug_slopes.csv")

slope_mat <- slope_tbl |>
  select(organism, antibiotic, slope_ppt_per_yr) |>
  pivot_wider(names_from = antibiotic, values_from = slope_ppt_per_yr) |>
  column_to_rownames("organism") |>
  as.matrix()
drugs_complete <- colnames(slope_mat)[colSums(is.na(slope_mat)) == 0]
slope_mat_complete <- slope_mat[, drugs_complete, drop = FALSE]

cat(sprintf("Trajectory matrix: %d organisms x %d drugs.\n",
            nrow(slope_mat_complete), ncol(slope_mat_complete)))

n_org   <- nrow(slope_mat_complete)
k_range <- 2:max(2, n_org - 1)
sil_widths <- sapply(k_range, function(k) pam(slope_mat_complete, k = k, diss = FALSE)$silinfo$avg.width)
k_best <- k_range[which.max(sil_widths)]
cat(sprintf("Best k by average silhouette width: k=%d (silhouette=%.3f)\n", k_best, max(sil_widths)))

pam_final <- pam(slope_mat_complete, k = k_best, diss = FALSE)
cluster_tbl <- tibble(
  organism    = rownames(slope_mat_complete),
  pam_cluster = pam_final$clustering,
  silhouette  = pam_final$silinfo$widths[match(rownames(slope_mat_complete),
                                               rownames(pam_final$silinfo$widths)), "sil_width"]
)
write_csv(cluster_tbl, "output/table_pam_trajectory_clusters.csv")
cat("\nCluster assignments:\n"); print(cluster_tbl)

heat_long <- slope_tbl |>
  filter(antibiotic %in% drugs_complete) |>
  left_join(cluster_tbl |> select(organism, pam_cluster), by = "organism") |>
  mutate(organism = factor(organism, levels = cluster_tbl$organism[order(cluster_tbl$pam_cluster)]))

p_heat <- ggplot(heat_long, aes(x = antibiotic, y = organism, fill = slope_ppt_per_yr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.1f", slope_ppt_per_yr)), size = 3) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                       midpoint = 0, name = "Slope\n(pp/yr R)") +
  labs(title = sprintf("Resistance trajectories (pp/yr) - k=%d trajectory clusters", k_best),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave("output/fig_trajectory_slope_heatmap.png", p_heat, width = 9, height = 4.5, dpi = 200)

png("output/fig_silhouette_trajectory.png", width = 1400, height = 900, res = 200)
plot(silhouette(pam_final), main = sprintf("Silhouette - trajectory clustering (k=%d)", k_best), col = "steelblue")
dev.off()

cat("\n  Wrote output/table_organism_drug_slopes.csv, table_pam_trajectory_clusters.csv, and figures.\n")
