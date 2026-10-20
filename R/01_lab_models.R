################################################################################
## 01_lab_models.R -- laboratory evaluation: logistic G-score curves and
##                    batch-specific variation in CCA-spiked urine
################################################################################

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(runjags); library(coda); library(ggplot2); library(cowplot) })

REFIT <- FALSE
NG_PER_UNIT <- 0.03

LAB <- c(CCA = "spiked_POC-CCA.csv", CCA3 = "spiked_POC-CCA3.csv")
TEST_LABEL <- c(CCA = "POC-CCA", CCA3 = "POC-CCA3")

read_lab <- function(test) {
  dt <- read.csv(file.path(PATHS$lab, LAB[[test]]), header = FALSE, fileEncoding = "UTF-8-BOM")[, 1:4]
  names(dt) <- c("dose", "Batch_1", "Batch_2", "Batch_3")
  stopifnot(nrow(dt) == 36, !anyNA(dt))
  dt
}

lab_model <- "model {
  for (i in 1:n){
    x1[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(log(y[i]+ 1e-6)-multiParam[2]))),tau1)
    x2[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(log(y[i]+ 1e-6)-multiParam[2]))),tau2)
    x3[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(log(y[i]+ 1e-6)-multiParam[2]))),tau3)
  }
  tau1 ~ dgamma(0.001,0.001)
  tau2 ~ dgamma(0.001,0.001)
  tau3 ~ dgamma(0.001,0.001)
  multiParam[1] ~ dnorm(0, 0.0001)
  multiParam[2] ~ dnorm(0, 0.0001)
}"

fit_lab <- function(test) {
  dt <- read_lab(test)
  dat <- list(n = nrow(dt), x1 = dt$Batch_1 - 1, x2 = dt$Batch_2 - 1, x3 = dt$Batch_3 - 1, y = dt$dose)
  inits <- lapply(1:4, function(ch) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 100 + ch))
  res <- run.jags(lab_model, data = dat, monitor = c("tau1", "tau2", "tau3", "multiParam"),
                  inits = inits, n.chains = 4, burnin = 3000, sample = 10000, silent.jags = TRUE)
  post <- as.mcmc.list(res)
  saveRDS(post, file.path(PATHS$post, paste0("lab_", test, ".rds")))
  post
}

lab_post <- lapply(names(LAB), function(t) {
  f <- file.path(PATHS$post, paste0("lab_", t, ".rds"))
  if (REFIT || !file.exists(f)) fit_lab(t) else readRDS(f)
})
names(lab_post) <- names(LAB)

## ---- Figure 1: fitted saturation curves with lab readings -------------------------
set.seed(1)
curve_panel <- function(test, cols) {
  m  <- draws_matrix(lab_post[[test]])
  mp <- m[sample.int(nrow(m), 1000), c("multiParam[1]", "multiParam[2]")]
  x  <- seq(10, 10000, by = 10)                  # 0.3-300 ng/ml, as in the original figure
  g  <- sapply(seq_len(nrow(mp)), function(j) 1 + 9 / (1 + exp(-mp[j, 1] * (log(x) - mp[j, 2]))))
  curves <- data.frame(ng = rep(x * NG_PER_UNIT, ncol(g)), G = as.vector(g), id = rep(seq_len(ncol(g)), each = length(x)))
  mean_curve <- data.frame(ng = x * NG_PER_UNIT, G = rowMeans(g))
  obs <- read_lab(test)
  obs <- data.frame(ng = rep(obs$dose * NG_PER_UNIT, 3), G = unlist(obs[, 2:4]),
                    Batch = rep(paste("Batch", 1:3), each = nrow(obs)))
  obs <- subset(obs, ng > 0)                     # zero dose cannot sit on a log axis
  ggplot() +
    geom_line(data = curves, aes(ng, G, group = id), colour = "black", alpha = 0.1, linewidth = 1) +
    geom_line(data = mean_curve, aes(ng, G), colour = "grey", linetype = "dashed", linewidth = 1.3) +
    geom_point(data = obs, aes(ng, G, colour = Batch, shape = Batch), alpha = 0.5, size = 5) +
    scale_x_log10() + annotation_logticks(sides = "b") +
    scale_colour_manual(values = cols) + scale_shape_manual(values = c(16, 17, 18)) +
    labs(x = "Spiked CCA Concentration (ng/ml)", y = "G Score", title = TEST_LABEL[[test]],
         colour = NULL, shape = NULL) +
    theme_bw(base_size = 13)
}
fig1 <- function() print(plot_grid(
  curve_panel("CCA",  c("#46c19a", "#7f62b8", "#b84c7d")),
  curve_panel("CCA3", c("#b8533c", "#bc9b3c", "#6ca24d")), ncol = 2))

png(file.path(PATHS$fig, "fig1_lab_curves.png"), 1000, 420, res = 100); fig1(); dev.off()
tiff(file.path(PATHS$fig, "fig1_lab_curves.tif"), 3000, 1260, res = 300, compression = "lzw"); fig1(); dev.off()

## Summary of lab batch SDs
lab_sd <- do.call(rbind, lapply(names(LAB), function(t) {
  m <- draws_matrix(lab_post[[t]])
  do.call(rbind, lapply(1:3, function(b) {
    q <- quantile(1 / sqrt(m[, paste0("tau", b)]), c(.025, .5, .975))
    data.frame(test = TEST_LABEL[[t]], batch = b, sd_median = q[2], sd_lo = q[1], sd_hi = q[3])
  }))
}))
write.csv(lab_sd, file.path(PATHS$tab, "lab_batch_sd.csv"), row.names = FALSE)
message("Lab models done.")
