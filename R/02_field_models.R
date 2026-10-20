################################################################################
## 02_field_models.R -- Bayesian latent class models for the field data (JAGS)
##
## Two model families are fitted to each site (Tororo, Mayuge):
##
##   "batch"  Each G-score reading ~ N(logistic(intensity), tau_batch[b]).
##            Batch-specific precision only; informative gamma(204.6, 956.9)
##            prior on every precision.
##
##   "inter"  Adds a day-level latent G-score per child, X[n,d] ~ N(logistic,
##            tau_inter), so day-to-day variation is separated from
##            batch-specific (within-day) variation.
##
## REFIT <- TRUE to refit everything.
################################################################################

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(runjags); library(rjags); library(coda) })

REFIT   <- FALSE
MCMC    <- list(burnin = 1000, sample = 5000, thin = 1, n.chains = 2)

## Fixed RNG streams per chain
.RNG.seed <- function(chain) switch(chain, "1" = 1, "2" = 2)
.RNG.name <- function(chain) switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill")

## Informative precision prior used for batch precisions in the batch model and
## for tau_inter in Tororo's inter model.
TAU_INFORMATIVE <- "dgamma(204.6194,956.8694)"
TAU_VAGUE       <- "dgamma(0.001,0.001)"

## Shared priors on the infection-intensity distribution
INTENSITY_PRIORS <- "
  P ~ dbeta(1, 1)
  k  ~ dgamma(0.001,0.001)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  multiParamCCA[1:2] ~ dmnorm.vcov(multim, multicov)
  multiParamRECCCA[1] ~ dnorm(multim[1],1/multicov[1,1])
  multiParamRECCCA[2] ~ dnorm(multim[2],1/multicov[2,2])
"

model_batch <- function() paste0("model {
  for (n in 1:N){
    status[n] ~ dbern(P)
    infect_intensity[n] ~ dgamma(sh,rt)
    for (r in 1:nKK){
      KK[n,r] ~ dnegbin(k/(infect_intensity[n]*status[n] +k),k)
    }
    for(d in 1:5){
      CCAintensity[n,d,1] <- 0
      CCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamCCA[1]*(infect_intensity[n]-multiParamCCA[2]))))
      RECCCAintensity[n,d,1] <- 0
      RECCCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamRECCCA[1]*(infect_intensity[n]-multiParamRECCCA[2]))))
      CCA[n,d] ~ dnorm(CCAintensity[n,d,status[n]+1],tauCCA_batch[CCA_batch[n,d]])
      REC_CCA[n,d] ~ dnorm(RECCCAintensity[n,d,status[n]+1],tauRECCCA_batch[REC_CCA_batch[n,d]])
    }
  }
", INTENSITY_PRIORS, "
  for (b in 1:4) {                       # b = 4 is a dummy for missing tests
    tauCCA_batch[b] ~ ", TAU_INFORMATIVE, "
    tauRECCCA_batch[b] ~ ", TAU_INFORMATIVE, "
  }
}")

model_inter <- function(tau_inter_prior) paste0("model {
  for (n in 1:N){
    status[n] ~ dbern(P)
    infect_intensity[n] ~ dgamma(sh,rt)
    for (r in 1:nKK){
      KK[n,r] ~ dnegbin(k/(infect_intensity[n]*status[n] +k),k)
    }
    for(d in 1:3){
      CCAintensity[n,d,1] <- 0
      CCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamCCA[1]*(infect_intensity[n]-multiParamCCA[2]))))
      RECCCAintensity[n,d,1] <- 0
      RECCCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamRECCCA[1]*(infect_intensity[n]-multiParamRECCCA[2]))))
      X_CCA[n,d,1] ~ dnorm(CCAintensity[n,d,1], tauCCA_inter)
      X_CCA[n,d,2] ~ dnorm(CCAintensity[n,d,2], tauCCA_inter)
      X_RECCCA[n,d,1] ~ dnorm(RECCCAintensity[n,d,1], tauRECCCA_inter)
      X_RECCCA[n,d,2] ~ dnorm(RECCCAintensity[n,d,2], tauRECCCA_inter)
    }
    for (r in 1:3){                      # day 1: three reads, one per batch
      CCA1[n,r] ~ dnorm(X_CCA[n,1,status[n]+1],tauCCA_batch[CCA_batch[n,r]])
      REC_CCA1[n,r] ~ dnorm(X_RECCCA[n,1,status[n]+1],tauRECCCA_batch[REC_CCA_batch[n,r]])
    }
    CCA2[n] ~ dnorm(X_CCA[n,2,status[n]+1],tauCCA_batch[CCA_batch[n,4]])
    CCA3[n] ~ dnorm(X_CCA[n,3,status[n]+1],tauCCA_batch[CCA_batch[n,5]])
    REC_CCA2[n] ~ dnorm(X_RECCCA[n,2,status[n]+1],tauRECCCA_batch[REC_CCA_batch[n,4]])
    REC_CCA3[n] ~ dnorm(X_RECCCA[n,3,status[n]+1],tauRECCCA_batch[REC_CCA_batch[n,5]])
  }
", INTENSITY_PRIORS, "
  tauCCA_inter ~ ", tau_inter_prior, "
  tauRECCCA_inter ~ ", tau_inter_prior, "
  for (b in 1:4) {                       # b = 4 is a dummy for missing tests
    tauCCA_batch[b] ~ ", TAU_VAGUE, "
    tauRECCCA_batch[b] ~ ", TAU_VAGUE, "
  }
}")

## What was run for each site
MODELS <- list(
  batch = list(Tororo = model_batch(), Mayuge = model_batch()),
  inter = list(Tororo = model_inter(TAU_INFORMATIVE), Mayuge = model_inter(TAU_VAGUE))
)

MONITOR <- list(
  batch = c("P", "sh", "rt", "tauCCA_batch", "tauRECCCA_batch", "k",
            "multiParamCCA", "multiParamRECCCA"),
  inter = c("P", "sh", "rt", "tauCCA_inter", "tauCCA_batch", "tauRECCCA_inter",
            "tauRECCCA_batch", "k", "multiParamCCA", "multiParamRECCCA")
)

jags_data <- function(model, site) {
  d  <- model_data(site)
  pr <- logistic_prior()
  common <- list(N = d$N, nKK = d$nKK, KK = d$KK,
                 CCA_batch = d$CCA_batch, REC_CCA_batch = d$REC_CCA_batch,
                 multim = pr$multim, multicov = pr$multicov)
  if (model == "batch") {
    c(common, list(CCA = d$CCA, REC_CCA = d$REC_CCA))
  } else {
    c(common, list(CCA1 = d$CCA[, 1:3], CCA2 = d$CCA[, 4], CCA3 = d$CCA[, 5],
                   REC_CCA1 = d$REC_CCA[, 1:3], REC_CCA2 = d$REC_CCA[, 4],
                   REC_CCA3 = d$REC_CCA[, 5]))
  }
}

fit_field_model <- function(model, site) {
  dat <- jags_data(model, site)
  inits <- function(chain) list(.RNG.seed = .RNG.seed(chain), .RNG.name = .RNG.name(chain),
                                status = rep(1, dat$N))
  res <- run.jags(MODELS[[model]][[site]], data = dat, monitor = MONITOR[[model]],
                  inits = list(inits(1), inits(2)),
                  burnin = MCMC$burnin, sample = MCMC$sample, thin = MCMC$thin,
                  n.chains = MCMC$n.chains, method = "parallel", silent.jags = TRUE)
  post <- as.mcmc.list(res)
  attr(post, "fit_info") <- list(model = model, site = site, N = dat$N,
                                 date = Sys.time(), mcmc = MCMC)
  saveRDS(post, posterior_file(model, site))      # replaces save.image()
  invisible(post)
}

if (sys.nframe() == 0L || isTRUE(getOption("cca3.run_all"))) {
  for (model in names(MODELS)) for (site in names(SITES)) {
    f <- posterior_file(model, site)
    if (REFIT || !file.exists(f)) {
      message("Fitting ", model, " model for ", site, " ...")
      fit_field_model(model, site)
    } else {
      message("Posterior exists, skipping: ", basename(f))
    }
  }
}
