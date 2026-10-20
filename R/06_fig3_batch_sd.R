################################################################################
## 06_fig3_batch_sd.R -- Figure 3: posterior SD of each batch, lab and field
################################################################################

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

MODEL <- if (exists("MODEL")) MODEL else "inter"
DISPLAY_BATCH <- list(Tororo = c(2, 1, 3), Mayuge = c(1, 2, 3), Lab = c(1, 2, 3))

sd_long <- function(m, pars, src) {
  do.call(rbind, lapply(1:3, function(b) data.frame(
    Batch = paste0("Batch_", DISPLAY_BATCH[[src]][b]), SD = 1 / sqrt(m[, pars[b]]))))
}

panels <- list()
lab <- lapply(c(CCA = "CCA", CCA3 = "CCA3"), function(t) draws_matrix(read_posterior("lab", t)))
for (t in c("CCA", "CCA3")) panels[[paste("Lab", t)]] <- sd_long(lab[[t]], paste0("tau", 1:3), "Lab")
for (s in c("Tororo", "Mayuge")) {
  m <- draws_matrix(read_posterior(MODEL, s))
  panels[[paste(s, "CCA")]]  <- sd_long(m, sprintf("tauCCA_batch[%d]", 1:3), s)
  panels[[paste(s, "CCA3")]] <- sd_long(m, sprintf("tauRECCCA_batch[%d]", 1:3), s)
}

FILLS <- list(CCA = c("#2b6f04", "#f1f1f1", "#0ade4c"), CCA3 = c("#003f5c", "#bc5090", "#ffa600"))
dens <- function(key) {
  t <- if (grepl("CCA3$", key)) "CCA3" else "CCA"
  where <- sub(" .*", "", key)
  ggplot(panels[[key]], aes(SD, fill = Batch, y = after_stat(scaled))) +
    geom_density(alpha = 0.5) + xlim(0, 3.5) +
    scale_fill_manual(values = FILLS[[t]], name = paste(if (t == "CCA") "CCA" else "CCA3", "Batches")) +
    labs(x = "Standard Deviation", y = "Density",
         title = paste(if (t == "CCA") "POC-CCA" else "POC-CCA3", "-", where)) +
    theme_bw()
}
row_ <- function(where, legend = NULL) plot_grid(
  dens(paste(where, "CCA")) + theme(legend.position = "none"),
  dens(paste(where, "CCA3")) + theme(legend.position = "none"),
  legend, ncol = 3, rel_widths = c(1, 1, .3))
fig3 <- function() print(plot_grid(
  row_("Lab", get_legend(dens("Lab CCA"))),
  row_("Tororo", get_legend(dens("Lab CCA3"))),
  row_("Mayuge"), nrow = 3))

dir.create(file.path(PATHS$fig, MODEL), showWarnings = FALSE)
png(file.path(PATHS$fig, MODEL, "fig3_batch_sd.png"), 1000, 700, res = 100); fig3(); dev.off()
tiff(file.path(PATHS$fig, MODEL, "fig3_batch_sd.tif"), 3000, 2100, res = 300, compression = "lzw"); fig3(); dev.off()

## Day-to-day SD (inter model only)
if (MODEL == "inter") {
  dd <- do.call(rbind, lapply(c("Tororo", "Mayuge"), function(s) {
    m <- draws_matrix(read_posterior("inter", s))
    rbind(data.frame(Setting = s, Test = "POC-CCA",  SD = 1 / sqrt(m[, "tauCCA_inter"])),
          data.frame(Setting = s, Test = "POC-CCA3", SD = 1 / sqrt(m[, "tauRECCCA_inter"])))
  }))
  p <- ggplot(dd, aes(SD, fill = Test, y = after_stat(scaled))) + geom_density(alpha = .5) +
    facet_wrap(~Setting) + scale_fill_manual(values = c("#2b6f04", "#003f5c")) +
    labs(x = "Day-to-day standard deviation", y = "Density") + theme_bw()
  ggsave(file.path(PATHS$fig, MODEL, "sfig_day_to_day_sd.png"), p, width = 8, height = 3.5, dpi = 150)
}
