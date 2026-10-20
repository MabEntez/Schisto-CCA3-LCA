################################################################################
## 03_diagnostics.R -- MCMC convergence and a posterior predictive check
##
## Supplementary Note 3 states that convergence was assessed with trace plots
## and the Gelman-Rubin statistic; this script produces and saves them.
##
## The posterior predictive check compares the observed proportion of children
## at or above G2/G3/G4 on the FIRST day-1 read (batch 1) with the proportion
## predicted by the fitted model (the model's own, untruncated, observation
## distribution). Large gaps indicate the Gaussian observation model is not
## reproducing the data at that threshold.
################################################################################

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(coda) })

fits <- expand.grid(model = c("batch", "inter"), site = names(SITES), stringsAsFactors = FALSE)
fits <- rbind(fits, data.frame(model = "lab", site = c("CCA", "CCA3")))

conv <- do.call(rbind, lapply(seq_len(nrow(fits)), function(i) {
  post <- read_posterior(fits$model[i], fits$site[i])
  ## drop the dummy precision for missing tests (index 4): it is sampled from
  ## its prior and is not a model parameter
  keep <- !grepl("_batch\\[4\\]", varnames(post))
  post <- post[, keep, drop = FALSE]
  pdf(file.path(PATHS$diag, sprintf("trace_%s_%s.pdf", fits$model[i], fits$site[i])), 8, 10)
  plot(post, density = FALSE); dev.off()
  q <- summary(post)$quantiles
  data.frame(model = fits$model[i], site = fits$site[i], parameter = varnames(post),
             median = q[, "50%"], lower95 = q[, "2.5%"], upper95 = q[, "97.5%"],
             Rhat = gelman.diag(post, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1],
             ESS = effectiveSize(post), row.names = NULL)
}))
write.csv(conv, file.path(PATHS$diag, "convergence_summary.csv"), row.names = FALSE)
message(sprintf("Max R-hat = %.3f; min ESS = %.0f (%s)", max(conv$Rhat), min(conv$ESS),
                paste(conv$model, conv$site, conv$parameter)[which.min(conv$ESS)]))

## ---- posterior predictive check (requires the field data) ----------------------
if (all(file.exists(file.path(PATHS$raw, SITES)))) {
  set.seed(7)
  TH <- c(G2 = 1, G3 = 2, G4 = 3)                          # 0-9 scale
  ppc <- do.call(rbind, lapply(c("batch", "inter"), function(model) do.call(rbind, lapply(names(SITES), function(s) {
    d <- model_data(s)
    m <- draws_matrix(read_posterior(model, s)); m <- m[sample.int(nrow(m), 1000), ]
    do.call(rbind, lapply(c(CCA = "CCA", CCA3 = "REC_CCA"), function(tv) {
      obs <- d[[tv]][, 1]; obs <- obs[!is.na(obs)]
      pre <- if (tv == "CCA") "CCA" else "RECCCA"
      sim <- sapply(seq_len(nrow(m)), function(j) {
        inf <- runif(length(obs)) < m[j, "P"]
        inten <- ifelse(inf, rgamma(length(obs), m[j, "sh"], m[j, "rt"]), 0)
        mu <- ifelse(inf, 9 / (1 + exp(-m[j, paste0("multiParam", pre, "[1]")] *
                                         (inten - m[j, paste0("multiParam", pre, "[2]")]))), 0)
        sdi <- if (paste0("tau", pre, "_inter") %in% colnames(m)) 1 / sqrt(m[j, paste0("tau", pre, "_inter")]) else 0
        s <- sqrt(1 / m[j, paste0("tau", pre, "_batch[1]")] + sdi^2)
        y <- round(rnorm(length(obs), mu, s))
        sapply(TH, function(t) mean(y >= t))
      })
      data.frame(model = model, site = s, test = names(which(c(CCA = "CCA", CCA3 = "REC_CCA") == tv)),
                 threshold = names(TH),
                 observed = sapply(TH, function(t) mean(obs >= t)),
                 predicted = apply(sim, 1, median),
                 pred_lo = apply(sim, 1, quantile, .025), pred_hi = apply(sim, 1, quantile, .975))
    }))
  }))))
  ppc$within_95 <- ppc$observed >= ppc$pred_lo & ppc$observed <= ppc$pred_hi
  write.csv(ppc, file.path(PATHS$diag, "ppc_first_read.csv"), row.names = FALSE)
  print(ppc, digits = 3, row.names = FALSE)
}

