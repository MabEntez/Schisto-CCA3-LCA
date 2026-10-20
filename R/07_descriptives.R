################################################################################
## 07_descriptives.R -- Table 1, Spearman's correlations, Supplementary Fig 2
## Fixes:
##   * Average KK intensity is the mean EPG over all smears. The original
##     mean(unlist(dt)) averaged EPG and G-scores together (Table 1 showed 6 and
##     82.6 EPG; the correct values are ~11.5 and ~210.5).
##   * Mayuge POC-CCA3 correlation with KK used POC-CCA day-1 columns (C1).
##   * Spearman's p-values from cor.test treat repeated reads of the same child
##     as independent; rho is therefore also reported with a 95% interval from
##     a bootstrap that resamples children.
################################################################################

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

EPG_FACTOR <- 24                # 41.7 mg template
GS <- c(G2 = 2, G2.5 = 2.5, G3 = 3, G4 = 4)
DAYS <- list(d1 = 1:3, d2 = 4, d3 = 5)    # positions within COLS$CCA / COLS$CCA3
KK_DAYS <- list(d1 = 1:2, d2 = 3:4, d3 = 5:6)

pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)
prev_by_g <- function(score) paste(sprintf("%s: %.1f", names(GS), sapply(GS, function(g) pct(score >= g))), collapse = ", ")

## ---- Table 1 ---------------------------------------------------------------------
tab1 <- list()
for (s in names(SITES)) {
  dt  <- read_site(s)
  kk  <- as.matrix(dt[, COLS$KK])
  cca <- as.matrix(dt[, COLS$CCA]); cca3 <- as.matrix(dt[, COLS$CCA3])
  ## per-child mean over all five readings (children with any missing reading
  ## are excluded)
  all_cca  <- rowMeans(cca);  all_cca3 <- rowMeans(cca3)
  tab1[[s]] <- c(
    n                              = nrow(dt),
    KK_prev_all_pct                = pct(rowSums(kk) > 0),
    KK_prev_first_pct              = pct(rowSums(kk[, 1:2]) > 0),
    KK_mean_EPG                    = round(mean(kk * EPG_FACTOR, na.rm = TRUE), 1),
    CCA_prev_all                   = prev_by_g(all_cca),
    CCA_mean_G_all                 = round(mean(cca, na.rm = TRUE), 1),
    CCA_prev_first                 = prev_by_g(cca[, 1]),
    CCA_mean_G_first               = round(mean(cca[, 1], na.rm = TRUE), 1),
    CCA3_prev_all                  = prev_by_g(all_cca3),
    CCA3_mean_G_all                = round(mean(cca3, na.rm = TRUE), 1),
    CCA3_prev_first                = prev_by_g(cca3[, 1]),
    CCA3_mean_G_first              = round(mean(cca3[, 1], na.rm = TRUE), 1)
  )
}
tab1 <- data.frame(Metric = names(tab1[[1]]), Tororo = tab1$Tororo, Mayuge = tab1$Mayuge, row.names = NULL)
write.csv(tab1, file.path(PATHS$tab, "table1_descriptives.csv"), row.names = FALSE)

## ---- Spearman's correlations -----------------------------------------------------
daily_mean <- function(m, days) do.call(c, lapply(days, function(j) rowMeans(m[, j, drop = FALSE], na.rm = TRUE)))

rho_boot <- function(x, y, id, B = 2000, seed = 2026) {
  ok <- complete.cases(x, y) & is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; id <- id[ok]
  rho <- cor(x, y, method = "spearman")
  rows <- split(seq_along(id), id)
  set.seed(seed)
  b <- replicate(B, { i <- unlist(rows[sample.int(length(rows), replace = TRUE)], use.names = FALSE)
                      suppressWarnings(cor(x[i], y[i], method = "spearman")) })
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(rho = rho, lwr = unname(quantile(b, .025, na.rm = TRUE)), upr = unname(quantile(b, .975, na.rm = TRUE)),
             n_pairs = length(x), n_children = length(rows), p_naive = ct$p.value)
}

sp <- list()
for (s in names(SITES)) {
  dt  <- read_site(s); id <- seq_len(nrow(dt))
  kk  <- as.matrix(dt[, COLS$KK]); cca <- as.matrix(dt[, COLS$CCA]); cca3 <- as.matrix(dt[, COLS$CCA3])
  kk_d <- daily_mean(kk, KK_DAYS)
  sp[[s]] <- rbind(
    cbind(Setting = s, Comparison = "POC-CCA3 vs POC-CCA (all paired reads)",
          rho_boot(as.vector(cca3), as.vector(cca), rep(id, 5))),
    cbind(Setting = s, Comparison = "POC-CCA3 vs KK (daily means)",
          rho_boot(daily_mean(cca3, DAYS), kk_d, rep(id, 3))),
    cbind(Setting = s, Comparison = "POC-CCA vs KK (daily means)",
          rho_boot(daily_mean(cca, DAYS), kk_d, rep(id, 3))))
}
sp <- do.call(rbind, sp)
write.csv(sp, file.path(PATHS$tab, "spearman.csv"), row.names = FALSE)
print(sp, digits = 3)

## ---- Supplementary Figure 2: G-score by KK intensity category ---------------------
INT_LEVELS <- c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram")
sfig2_panel <- function(s) {
  dt <- read_site(s)
  epg <- rowMeans(as.matrix(dt[, COLS$KK]) * EPG_FACTOR, na.rm = TRUE)
  cat_ <- cut(epg, c(-Inf, 0, 100, 400, Inf), labels = INT_LEVELS, include.lowest = TRUE)
  d <- rbind(data.frame(intensity = cat_, assay = "POC-CCA",  G = rowMeans(as.matrix(dt[, COLS$CCA]),  na.rm = TRUE)),
             data.frame(intensity = cat_, assay = "POC-CCA3", G = rowMeans(as.matrix(dt[, COLS$CCA3]), na.rm = TRUE)))
  d <- d[!is.na(d$intensity) & is.finite(d$G), ]
  ggplot(d, aes(intensity, G, fill = intensity)) +
    geom_boxplot(show.legend = FALSE) + facet_wrap(~assay, nrow = 1) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = grey.colors(4, start = 1, end = 0)) +
    scale_y_continuous(breaks = seq(2, 10, 2), limits = c(1, 10)) +
    labs(x = "Kato-Katz eggs per gram", y = "Mean G-score", title = s) +
    theme_minimal(base_size = 13)
}
sfig2 <- function() print(plot_grid(sfig2_panel("Tororo"), sfig2_panel("Mayuge"), ncol = 1))
png(file.path(PATHS$fig, "sfig2_gscore_by_kk.png"), 1100, 800, res = 100); sfig2(); dev.off()
tiff(file.path(PATHS$fig, "sfig2_gscore_by_kk.tif"), 3300, 2400, res = 300, compression = "lzw"); sfig2(); dev.off()
