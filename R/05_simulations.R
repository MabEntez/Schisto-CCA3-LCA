################################################################################
## 05_simulations.R -- predicted prevalence, test performance, ROC/AUC, WHO TPP
##                     and the 10%-prevalence scenario
################################################################################

source(here::here("R", "00_setup.R"))
source(here::here("R", "04_sim_functions.R"))
suppressPackageStartupMessages({ library(coda); library(ggplot2); library(cowplot) })

MODEL  <- if (exists("MODEL")) MODEL else "inter"
NDRAWS <- 1000      # posterior draws
NIND   <- 1000      # children per simulated population (sens/spec, AUC)
NTPP   <- 100       # children per simulated survey (TPP)
SEED   <- 100

fig_dir <- file.path(PATHS$fig, MODEL); dir.create(fig_dir, showWarnings = FALSE)
tab_dir <- file.path(PATHS$tab, MODEL); dir.create(tab_dir, showWarnings = FALSE)
TESTS  <- c(CCA = "POC-CCA", CCA3 = "POC-CCA3")

save_fig <- function(plot_fun, name, width, height) {
  png(file.path(fig_dir, paste0(name, ".png")), width = width, height = height, res = 100)
  plot_fun(); dev.off()
  tiff(file.path(fig_dir, paste0(name, ".tif")), width = width * 3, height = height * 3,
       res = 300, compression = "lzw")
  plot_fun(); dev.off()
}

q3 <- function(x) quantile(x, c(.025, .5, .975), na.rm = TRUE)

## ---- posterior draws, aligned across parameters ------------------------
SITES_ <- names(SITES)
draws <- lapply(SITES_, function(s) posterior_draws(read_posterior(MODEL, s), NDRAWS, SEED))
names(draws) <- SITES_

## ---- Figure 2: predicted prevalence ---------------------------------------------
prev_df <- do.call(rbind, lapply(SITES_, function(s) data.frame(site = s, P = 100 * draws[[s]]$P)))
prev_df$site <- factor(prev_df$site, levels = c("Tororo", "Mayuge"))
prev_summary <- aggregate(P ~ site, prev_df, function(x) round(q3(x), 1))
write.csv(do.call(data.frame, prev_summary), file.path(tab_dir, "predicted_prevalence.csv"), row.names = FALSE)

fig2 <- function() {
  panels <- lapply(levels(prev_df$site), function(s) {
    ggplot(subset(prev_df, site == s), aes(P)) +
      geom_histogram(binwidth = 1, fill = "grey45", colour = "black") +
      scale_x_continuous(labels = function(x) paste0(x, "%")) +
      labs(title = c(Tororo = "A", Mayuge = "B")[[s]], x = "Prevalence", y = "Frequency") +
      theme_minimal(base_size = 13)
  })
  print(plot_grid(plotlist = panels, ncol = 2))
}
save_fig(fig2, "fig2_predicted_prevalence", 800, 340)

## ---- Sensitivity / specificity simulations (fix S1-S4) --------------------------
set.seed(SEED)
sim <- lapply(SITES_, function(s) {
  d <- draws[[s]]
  res <- lapply(seq_len(nrow(d)), function(i) {
    pop <- simulate_population(NIND, d[i, ])
    rbind(cbind(test = "CCA",  sens_spec(pop, "CCA")),
          cbind(test = "CCA3", sens_spec(pop, "CCA3")),
          cbind(test = "KK", threshold = "any egg", sens_spec_kk(pop)))[, c("test", "days", "threshold", "sens", "spec")]
  })
  out <- do.call(rbind, Map(function(r, i) cbind(draw = i, r), res, seq_along(res)))
  cbind(site = s, out)
})
sim <- do.call(rbind, sim)
saveRDS(sim, file.path(tab_dir, "sens_spec_draws.rds"))

## Supplementary Table 3: median (95% CrI)
cs <- subset(sim, test != "KK")
st3 <- do.call(rbind, lapply(split(cs, list(cs$site, cs$test, cs$days, cs$threshold), drop = TRUE), function(x) {
  se <- q3(x$sens); sp <- q3(x$spec)
  data.frame(Setting = x$site[1], Diagnostic = TESTS[[x$test[1]]], Days = x$days[1],
             Threshold = x$threshold[1],
             Sensitivity = fmt_ci(se[2], se[1], se[3]), Specificity = fmt_ci(sp[2], sp[1], sp[3]),
             sens_median = se[2], spec_median = sp[2])
}))
st3 <- st3[order(match(st3$Setting, c("Tororo", "Mayuge")), st3$Diagnostic, st3$Days,
                 match(st3$Threshold, names(THRESHOLDS))), ]
write.csv(st3, file.path(tab_dir, "stable3_sens_spec.csv"), row.names = FALSE)

## Supplementary Table 4: KK sensitivity (specificity is 100% by construction)
kk <- subset(sim, test == "KK")
st4 <- do.call(rbind, lapply(split(kk, list(kk$site, kk$days), drop = TRUE), function(x) {
  se <- q3(x$sens)
  data.frame(Setting = x$site[1], Days = x$days[1], Sensitivity = fmt_ci(se[2], se[1], se[3]))
}))
st4 <- st4[order(match(st4$Setting, c("Tororo", "Mayuge")), st4$Days), ]
write.csv(st4, file.path(tab_dir, "stable4_kk_sens.csv"), row.names = FALSE)

## ---- Table 2: AUC per draw ------------------------------------------------------
auc <- do.call(rbind, lapply(split(cs, list(cs$site, cs$test, cs$days, cs$draw), drop = TRUE), function(x)
  data.frame(site = x$site[1], test = x$test[1], days = x$days[1], draw = x$draw[1],
             auc = auc_points(x$sens, x$spec))))
tab2 <- do.call(rbind, lapply(split(auc, list(auc$days, auc$site, auc$test), drop = TRUE), function(x) {
  q <- q3(x$auc); data.frame(Days = x$days[1], col = paste(x$site[1], TESTS[[x$test[1]]]),
                             AUC = fmt_ci(q[2], q[1], q[3]))
}))
tab2 <- reshape(tab2, idvar = "Days", timevar = "col", direction = "wide")
names(tab2) <- sub("^AUC\\.", "", names(tab2))
tab2 <- tab2[order(tab2$Days), c("Days", "Tororo POC-CCA", "Tororo POC-CCA3", "Mayuge POC-CCA", "Mayuge POC-CCA3")]
write.csv(tab2, file.path(tab_dir, "table2_auc.csv"), row.names = FALSE)

## AUC comparison.
cmp <- do.call(rbind, lapply(split(auc, list(auc$site, auc$days), drop = TRUE), function(x) {
  a1 <- x$auc[x$test == "CCA"][order(x$draw[x$test == "CCA"])]
  a3 <- x$auc[x$test == "CCA3"][order(x$draw[x$test == "CCA3"])]
  data.frame(Setting = x$site[1], Days = x$days[1],
             median_diff = median(a3 - a1),
             prob_CCA3_higher = mean(a3 > a1),
             wilcoxon_p_one_tailed = suppressWarnings(wilcox.test(a1, a3, paired = TRUE, alternative = "less")$p.value))
}))
write.csv(cmp[order(cmp$Setting, cmp$Days), ], file.path(tab_dir, "auc_comparison.csv"), row.names = FALSE)

## ---- Figure 4: ROC curves -------------------------------------------------------
## Points at G4, G3 and G2 (left to right) for 1, 2 and 3 days of sampling.
roc_q <- do.call(rbind, lapply(split(cs, list(cs$site, cs$test, cs$days, cs$threshold), drop = TRUE), function(x) {
  se <- q3(x$sens); fp <- q3(1 - x$spec)
  data.frame(site = x$site[1], test = x$test[1], days = x$days[1], threshold = x$threshold[1],
             fpr = fp[2], fpr_lo = fp[1], fpr_hi = fp[3], sens = se[2], sens_lo = se[1], sens_hi = se[3])
}))
roc_q <- subset(roc_q, threshold %in% c("G2", "G3", "G4"))
DAY_COL <- c("#f78f0b", "#59ce48", "#8628af")

fig4 <- function() {
  par(mfcol = c(2, 2), font = 2, lwd = 2, mar = c(4, 4, 2, 1), mgp = c(2.2, 0.6, 0), oma = c(0, 2, 1.5, 0))
  for (s in c("Tororo", "Mayuge")) for (t in c("CCA", "CCA3")) {
    plot(NA, xlim = c(0, 0.6), ylim = c(0.7, 1), xlab = "1 - Specificity", ylab = "Sensitivity",
         las = 1, cex.lab = 1.2)
    for (day in 1:3) {
      r <- subset(roc_q, site == s & test == t & days == day)
      r <- r[order(r$fpr), ]
      lines(r$fpr, r$sens, col = DAY_COL[day]); points(r$fpr, r$sens, pch = 19, col = DAY_COL[day])
      arrows(r$fpr, r$sens_lo, r$fpr, r$sens_hi, angle = 90, code = 3, length = .04, col = DAY_COL[day])
      arrows(r$fpr_lo, r$sens, r$fpr_hi, r$sens, angle = 90, code = 3, length = .04, col = DAY_COL[day])
      if (s == "Tororo" && t == "CCA" && day == 1) text(r$fpr, r$sens, r$threshold, pos = 4, cex = .8)
    }
    if (s == "Mayuge" && t == "CCA3") legend("bottomright", c("1 Day", "2 Days", "3 Days"),
                                             col = DAY_COL, lty = 1, bty = "n")
  }
  mtext(c("Tororo", "Mayuge"), side = 3, outer = TRUE, at = c(.27, .77), cex = 1.3)
  mtext(c("POC-CCA", "POC-CCA3"), side = 2, outer = TRUE, at = c(.77, .27), cex = 1.3)
}
save_fig(fig4, "fig4_roc", 800, 600)

## ---- Figure 5: probability of meeting the WHO TPP --------------------------------
set.seed(SEED + 1)
tpp <- do.call(rbind, lapply(SITES_, function(s) {
  d <- draws[[s]]
  do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    pop <- simulate_population(NTPP, d[i, ])
    rbind(cbind(site = s, test = "CCA",  sens_spec(pop, "CCA")),
          cbind(site = s, test = "CCA3", sens_spec(pop, "CCA3")))
  }))
}))
tpp_pct <- do.call(rbind, lapply(names(TPP), function(lvl) {
  ok <- tpp$sens >= TPP[[lvl]]["sens"] & tpp$spec >= TPP[[lvl]]["spec"]
  a <- aggregate(ok ~ site + test + days + threshold, data = cbind(tpp, ok = ok), FUN = function(z) 100 * mean(z))
  cbind(level = lvl, a)
}))
names(tpp_pct)[names(tpp_pct) == "ok"] <- "percent"
tpp_pct$setting <- factor(paste(tpp_pct$site, tpp_pct$test),
                          levels = c("Tororo CCA", "Tororo CCA3", "Mayuge CCA", "Mayuge CCA3"),
                          labels = c("A", "B", "C", "D"))
write.csv(tpp_pct, file.path(tab_dir, "fig5_tpp_percent.csv"), row.names = FALSE)

fig5 <- function() {
  dd <- subset(tpp_pct, threshold %in% c("G2", "G3", "G4"))
  dd$threshold <- factor(dd$threshold, labels = c("2 G score threshold", "3 G score threshold", "4 G score threshold"))
  dd$level <- factor(dd$level, levels = c("optimum", "minimum"),
                     labels = c("Percentage achieving\noptimum TPP", "Percentage achieving\nminimum TPP"))
  dd$days <- factor(dd$days, labels = paste(1:3, c("day sampling", "days sampling", "days sampling")))
  p <- ggplot(dd, aes(setting, percent, fill = threshold)) +
    geom_col(position = "dodge", colour = "black", width = .9) +
    facet_grid(level ~ days, switch = "y") +
    scale_fill_manual(values = c("grey0", "grey50", "grey100"), name = "Threshold") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
    labs(x = "Setting", y = NULL) +
    theme_bw(base_size = 14) + theme(strip.placement = "outside", strip.background = element_blank())
  print(p)
}
save_fig(fig5, "fig5_tpp", 1300, 500)

## ---- Table 3: setting with ~10% prevalence by a single KK -------------------------
## Uses Tororo posterior draws with the lowest prevalence for the intensity
## distribution and test parameters, sets true prevalence to 20%, simulates
## 100 children, and keeps the 1,000 simulations whose single-smear KK
## prevalence is closest to 10%.
set.seed(SEED + 2)
mT <- draws_matrix(read_posterior(MODEL, "Tororo"))
low_rows <- order(mT[, "P"])[1:100]
d10 <- posterior_draws(read_posterior(MODEL, "Tororo"), rows = sample(low_rows, 10000, replace = TRUE))
p10 <- t(sapply(seq_len(nrow(d10)), function(i) {
  pop <- simulate_population(100, d10[i, ], P = 0.2)
  c(KK = mean(pop$KK[, 1] > 0), setNames(diag_prev(pop, "CCA"), paste0("CCA_", names(THRESHOLDS))),
    setNames(diag_prev(pop, "CCA3"), paste0("CCA3_", names(THRESHOLDS))))
}))
keep <- order(abs(p10[, "KK"] - 0.10))[1:1000]
p10 <- p10[keep, ]
summ <- function(x) { q <- quantile(x, c(.025, .975)); fmt_ci(mean(x), q[1], q[2]) }
tab3 <- data.frame(G_score = names(THRESHOLDS),
                   `POC-CCA prevalence`  = sapply(names(THRESHOLDS), function(t) summ(p10[, paste0("CCA_", t)])),
                   `POC-CCA3 prevalence` = sapply(names(THRESHOLDS), function(t) summ(p10[, paste0("CCA3_", t)])),
                   `KK prevalence` = summ(p10[, "KK"]), check.names = FALSE)
write.csv(tab3, file.path(tab_dir, "table3_prev10.csv"), row.names = FALSE)

message("Simulations for the '", MODEL, "' model written to ", tab_dir, " and ", fig_dir)
