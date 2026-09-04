# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# R ML modules: clustering/NMF resistance archetypes and PELT changepoints
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

library(tidyverse)
library(scales)

ensure_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      tryCatch(
        install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE),
        error = function(e) invisible(NULL)
      )
    }
  }
}

ensure_pkgs(c("changepoint", "cluster"))

library(changepoint)
library(cluster)

root <- getwd()
if (!file.exists(file.path(root, "R", "amrsn_data_shared.R"))) {
  stop("Run from project root; R/amrsn_data_shared.R not found.")
}
source(file.path(root, "R", "amrsn_data_shared.R"))

dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

poster_bg <- "#3C3C3C"
poster_red <- "#8B0000"
poster_text <- "#E8E8E8"
poster_white <- "#FFFFFF"
poster_grid <- "#525252"

theme_ml <- theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = poster_bg, colour = NA),
    panel.background = element_rect(fill = poster_bg, colour = NA),
    panel.grid.major = element_line(colour = poster_grid, linewidth = 0.3),
    text = element_text(colour = poster_text),
    plot.title = element_text(colour = poster_white, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "#AAAAAA", size = 9),
    axis.text = element_text(colour = poster_text),
    axis.title = element_text(colour = poster_text),
    legend.position = "bottom",
    legend.background = element_rect(fill = poster_bg, colour = NA)
  )

theme_ml_light <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    text             = element_text(colour = "#1f1f1f"),
    plot.title       = element_text(colour = "#1f1f1f", face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = 9),
    axis.text        = element_text(colour = "#1f1f1f"),
    axis.title       = element_text(colour = "#1f1f1f"),
    legend.position  = "bottom",
    legend.background = element_rect(fill = "white", colour = NA)
  )

full_long_ml <- full_long %>%
  mutate(pair_id = paste(organism, drug, sep = " | "))

cat("\n=== MODULE 9: CLUSTERING / ARCHETYPES ===\n")

abg_wide_2024 <- full_long %>%
  filter(year == 2024) %>%
  select(organism, drug, susc_pct) %>%
  pivot_wider(names_from = drug, values_from = susc_pct)

org_prof <- abg_wide_2024 %>%
  tibble::column_to_rownames("organism")
org_prof_imp <- org_prof
for (j in seq_len(ncol(org_prof_imp))) {
  ind <- is.na(org_prof_imp[, j])
  org_prof_imp[ind, j] <- mean(org_prof_imp[, j], na.rm = TRUE)
}

d_org <- dist(org_prof_imp, method = "euclidean")
hc_org <- hclust(d_org, method = "average")
write_csv(
  tibble(order = hc_org$order, organism = rownames(org_prof_imp)[hc_org$order]),
  "output/ml/table_hclust_organism_order.csv"
)

png("output/ml/fig_dendrogram_organisms_2024.png", width = 2400, height = 1400, res = 200)
plot(hc_org, main = "Hierarchical clustering - organism antibiogram profiles (2024)",
     xlab = "", sub = "")
dev.off()

set.seed(42)
pam_k <- cluster::pam(scale(as.matrix(org_prof_imp)), k = 3)
org_clusters <- tibble(
  organism = rownames(org_prof_imp),
  pam_cluster = as.integer(pam_k$clustering)
)
write_csv(org_clusters, "output/ml/table_pam_organism_clusters.csv")

pair_features <- full_long_ml %>%
  group_by(organism, drug, pair_id) %>%
  summarise(
    mean_susc = mean(susc_pct, na.rm = TRUE),
    delta_susc = susc_pct[year == max(year)] - susc_pct[year == min(year)],
    slope_w = {
      fit <- tryCatch(
        stats::lm(susc_pct ~ year, weights = tested),
        error = function(e) NULL
      )
      if (is.null(fit)) NA_real_ else unname(stats::coef(fit)["year"])
    },
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))

pfm <- pair_features %>%
  select(mean_susc, delta_susc, slope_w) %>%
  scale()
rownames(pfm) <- pair_features$pair_id

d_pair <- dist(pfm, method = "euclidean")
hc_pair <- hclust(d_pair, method = "average")
pair_cut <- tibble(
  pair_id = rownames(pfm),
  pair_cluster = as.integer(cutree(hc_pair, k = 4L))
)

pair_clusters <- pair_features %>%
  left_join(pair_cut, by = "pair_id")
write_csv(pair_clusters, "output/ml/table_pair_clusters_trajectory.csv")

V_raw <- full_long %>%
  filter(year == 2024) %>%
  select(drug, organism, susc_pct) %>%
  mutate(susc_pct = pmax(0, susc_pct)) %>%
  pivot_wider(names_from = organism, values_from = susc_pct) %>%
  tibble::column_to_rownames("drug")
V <- as.matrix(V_raw)
V[is.na(V)] <- min(V, na.rm = TRUE) * 0.5
V <- V / 100

if (!requireNamespace("NMF", quietly = TRUE)) {
  tryCatch({
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    BiocManager::install("Biobase", ask = FALSE, update = FALSE)
    install.packages("NMF", repos = "https://cloud.r-project.org", quiet = TRUE)
  }, error = function(e) {
    message("  NMF/Biobase install failed: ", conditionMessage(e))
  })
}
nmf_ok <- requireNamespace("NMF", quietly = TRUE)
if (nmf_ok) {
  nmff <- tryCatch(
    NMF::nmf(V, rank = 3, seed = 42, .options = list(maxIter = 500L, keep.all = FALSE)),
    error = function(e) NULL
  )
  if (!is.null(nmff)) {
    W <- NMF::basis(nmff)
    H <- NMF::coef(nmff)
    archetype_drug <- as_tibble(W, rownames = "drug") %>%
      pivot_longer(-drug, names_to = "archetype", values_to = "loading")
    archetype_org <- as_tibble(t(H), rownames = "organism") %>%
      pivot_longer(-organism, names_to = "archetype", values_to = "score")
    write_csv(archetype_drug, "output/ml/table_nmf_archetype_drug_loadings.csv")
    write_csv(archetype_org, "output/ml/table_nmf_archetype_organism_scores.csv")
    cat("  NMF archetypes (rank=3) saved.\n")
  }
} else {
  cat("  NMF package not installed; skipped (install.packages('NMF')).\n")
}

cat("  Clustering tables and dendrogram written to output/ml/.\n")

cat("\n=== MODULE 10: PELT CHANGEPOINT SENSITIVITY ===\n")

eskape_organisms <- c(
  "K. pneumoniae", "A. baumannii",
  "P. aeruginosa", "Enterobacter spp."
)

cpt_rows <- list()
for (org in eskape_organisms) {
  d <- amrsn %>% filter(organism == org) %>% arrange(year)
  x_mer <- d$mer_res_pct / 100
  x_imi <- d$imi_res_pct / 100
  cp_mer <- tryCatch(
    changepoint::cpt.meanvar(x_mer, method = "PELT"),
    error = function(e) NULL
  )
  cp_imi <- tryCatch(
    changepoint::cpt.meanvar(x_imi, method = "PELT"),
    error = function(e) NULL
  )
  add_cp <- function(cp_obj, abx) {
    if (is.null(cp_obj)) {
      return(tibble(organism = org, antibiotic = abx, n_cpt = NA_integer_,
                    cpt_years = NA_character_))
    }
    ix <- changepoint::cpts(cp_obj)
    yrs <- if (length(ix)) d$year[ix] else NA_integer_
    tibble(
      organism = org,
      antibiotic = abx,
      n_cpt = length(ix),
      cpt_years = paste(unique(yrs), collapse = ";")
    )
  }
  cpt_rows[[length(cpt_rows) + 1]] <- add_cp(cp_mer, "Meropenem")
  cpt_rows[[length(cpt_rows) + 1]] <- add_cp(cp_imi, "Imipenem")
}

cpt_tbl <- bind_rows(cpt_rows)
write_csv(cpt_tbl, "output/ml/table_changepoint_pelt_carbapenem.csv")

reg_list <- list()
for (org in eskape_organisms) {
  d <- amrsn %>% filter(organism == org) %>% arrange(year)
  n <- nrow(d)
  cp <- tryCatch(
    changepoint::cpt.meanvar(d$mer_res_pct / 100, method = "PELT"),
    error = function(e) NULL
  )
  if (is.null(cp)) {
    seg <- rep(1L, n)
  } else {
    cp_i <- sort(unique(as.integer(changepoint::cpts(cp))))
    br <- sort(unique(c(0L, cp_i, n)))
    seg <- integer(n)
    lab <- 1L
    for (j in seq_len(length(br) - 1L)) {
      lo <- br[j] + 1L
      hi <- br[j + 1L]
      seg[lo:hi] <- lab
      lab <- lab + 1L
    }
  }
  reg_list[[length(reg_list) + 1]] <- tibble(
    organism = org,
    year = d$year,
    mer_res_pct = d$mer_res_pct,
    segment = seg
  )
}
regime_long <- bind_rows(reg_list)
write_csv(regime_long, "output/ml/table_meropenem_regime_map_pelt.csv")

fig_regime <- regime_long %>%
  ggplot(aes(
    x = year,
    y = mer_res_pct,
    colour = factor(segment),
    group = interaction(organism, segment, drop = TRUE)
  )) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~organism, scales = "fixed", ncol = 2) +
  scale_x_continuous(breaks = 2017:2024) +
  labs(
    title = "Meropenem resistance with PELT mean-variance segments",
    x = NULL, y = "Resistance (%)", colour = "Segment"
  ) +
  theme_ml

ggsave("output/ml/fig_changepoint_meropenem_regimes.png", fig_regime, width = 10, height = 7, dpi = 300)
ggsave("output/ml/fig_changepoint_meropenem_regimes_light.png",
       fig_regime + theme_ml_light, width = 10, height = 7, dpi = 300)

cat("  Wrote changepoint and regime-map outputs.\n")

cat("\n  amrsn_ml_enhancements.R completed - see output/ml/\n")
