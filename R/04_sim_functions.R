################################################################################
## 04_sim_functions.R -- simulate test performance from the model posteriors
##
## Every function takes the posterior draw it needs as an argument (`d`, one row
## of the table returned by posterior_draws())
################################################################################

## G-score thresholds on the model (0-9) scale and their labels on the 1-10 scale
THRESHOLDS <- c(G2 = 1, G2.5 = 1.5, G3 = 2, G4 = 3)

## WHO diagnostic TPP for S. mansoni (sensitivity, specificity)
TPP <- list(minimum = c(sens = 0.60, spec = 0.95),
            optimum = c(sens = 0.75, spec = 0.965))

## Truncated normal via inverse CDF (vectorised over mu, s and lower)
rnormt <- function(n, lower, mu, s) {
  Fa <- pnorm(lower, mean = mu, sd = s)
  u  <- runif(n, min = Fa, max = 1)
  qnorm(u, mean = mu, sd = s)
}

## ---- posterior draws ------------------------------------------------------------
## One row per sampled MCMC iteration; all parameters from that same iteration.
posterior_draws <- function(post, ndraws = 1000, seed = 100, rows = NULL) {
  m <- draws_matrix(post)
  if (is.null(rows)) {
    set.seed(seed)
    rows <- sample.int(nrow(m), ndraws)
  }
  sd_of <- function(par) if (par %in% colnames(m)) 1 / sqrt(m[rows, par]) else rep(0, length(rows))
  data.frame(
    row       = rows,
    P         = m[rows, "P"],
    sh        = m[rows, "sh"],
    rt        = m[rows, "rt"],
    k         = m[rows, "k"],
    a_CCA     = m[rows, "multiParamCCA[1]"],
    b_CCA     = m[rows, "multiParamCCA[2]"],
    a_CCA3    = m[rows, "multiParamRECCCA[1]"],
    b_CCA3    = m[rows, "multiParamRECCCA[2]"],
    sdb1_CCA  = sd_of("tauCCA_batch[1]"),
    sdb2_CCA  = sd_of("tauCCA_batch[2]"),
    sdb3_CCA  = sd_of("tauCCA_batch[3]"),
    sdb1_CCA3 = sd_of("tauRECCCA_batch[1]"),
    sdb2_CCA3 = sd_of("tauRECCCA_batch[2]"),
    sdb3_CCA3 = sd_of("tauRECCCA_batch[3]"),
    sdi_CCA   = sd_of("tauCCA_inter"),      # 0 for the batch model
    sdi_CCA3  = sd_of("tauRECCCA_inter")
  )
}

## ---- simulate one population ---------------------------------------------------
## Returns infection status, intensity, three days of G-scores for each test
## (0-9 scale, one random batch per child per day) and six KK smear counts.
simulate_population <- function(nInd, d, P = d$P) {
  nInf <- rbinom(1, nInd, P)
  infected  <- c(rep(FALSE, nInd - nInf), rep(TRUE, nInf))
  intensity <- numeric(nInd)
  intensity[infected] <- rgamma(nInf, shape = d$sh, rate = d$rt)

  gscores <- function(test) {
    a  <- d[[paste0("a_", test)]];  b <- d[[paste0("b_", test)]]
    sdb <- c(d[[paste0("sdb1_", test)]], d[[paste0("sdb2_", test)]], d[[paste0("sdb3_", test)]])
    sdi <- d[[paste0("sdi_", test)]]
    mu    <- ifelse(infected, 9 / (1 + exp(-a * (intensity - b))), 0)
    lower <- ifelse(infected, 0, -3)
    sapply(1:3, function(day) {
      batch <- sample.int(3, nInd, replace = TRUE)
      rnormt(nInd, lower, mu, sqrt(sdb[batch]^2 + sdi^2))
    })
  }

  kk <- sapply(1:6, function(i) rnbinom(nInd, size = d$k, mu = intensity))

  list(infected = infected, CCA = gscores("CCA"), CCA3 = gscores("CCA3"), KK = kk)
}

## Score after 1, 2 or 3 days. A single day's reading is rounded to the
## nearest whole G-score (as read in the field); multi-day scores are the mean
## of the daily latent readings
day_scores <- function(g) cbind(day1 = round(g[, 1]),
                                day2 = rowMeans(g[, 1:2]),
                                day3 = rowMeans(g))

## Sensitivity/specificity for each day x threshold -> long data.frame
sens_spec <- function(pop, test) {
  s <- day_scores(pop[[test]])
  inf <- pop$infected
  do.call(rbind, lapply(1:3, function(day) do.call(rbind, lapply(names(THRESHOLDS), function(th) {
    pos <- s[, day] >= THRESHOLDS[[th]]
    data.frame(days = day, threshold = th,
               sens = if (any(inf)) mean(pos[inf]) else NA_real_,
               spec = if (any(!inf)) mean(!pos[!inf]) else NA_real_)
  }))))
}

sens_spec_kk <- function(pop) {
  inf <- pop$infected
  do.call(rbind, lapply(1:3, function(day) {
    pos <- rowSums(pop$KK[, seq_len(2 * day), drop = FALSE]) > 0     # S4
    data.frame(days = day, sens = mean(pos[inf]), spec = mean(!pos[!inf]))
  }))
}

## Trapezoidal AUC through (0,0), the four threshold points and (1,1)
auc_points <- function(sens, spec) {
  x <- c(0, 1 - spec, 1); y <- c(0, sens, 1)
  o <- order(x, y); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

## Observed (diagnostic) prevalence from a single day's test
diag_prev <- function(pop, test) {
  s <- round(pop[[test]][, 1])
  sapply(THRESHOLDS, function(t) mean(s >= t))
}
