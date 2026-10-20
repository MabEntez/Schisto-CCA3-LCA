################################################################################
## 00_setup.R -- shared setup for every script in the pipeline
##
## Sourced at the top of each numbered script. Defines paths (via {here}, so the
## pipeline runs from a clean checkout under Rscript or RStudio), the column
## names used from the field data, and the functions that turn a site's CSV
## into the data objects the JAGS models need.
##
## Conventions
##   * G-scores are shifted to a 0-9 scale (raw score - 1) for modelling. The
##     raw 1-10 scale is used only for descriptive prevalence (07_descriptives.R).
##     Model thresholds 1, 1.5, 2, 3 therefore correspond to G2, G2.5, G3, G4.
################################################################################

suppressPackageStartupMessages({
  library(here)
})

PATHS <- list(
  raw    = here::here("data", "raw"),
  lab    = here::here("data", "lab"),
  priors = here::here("data", "priors"),
  post   = here::here("outputs", "posteriors"),
  fig    = here::here("outputs", "figures"),
  tab    = here::here("outputs", "tables"),
  diag   = here::here("outputs", "diagnostics")
)

invisible(lapply(PATHS[c("post", "fig", "tab", "diag")], dir.create,
                 recursive = TRUE, showWarnings = FALSE))

SITES <- c(Tororo = "CORNTD_recPOCCCA_2023_Tororo.csv",
           Mayuge = "CORNTD_recPOCCCA_2023_Mayuge.csv")

## ---- field-data columns ----------------
## Day 1: three POC-CCA and three POC-CCA3 reads (one per batch); days 2-3: one.
## KK: duplicate S. mansoni slides per day.
COLS <- list(
  KK         = c("KK_Sm1A", "KK_Sm1B", "KK_Sm2A", "KK_Sm2B", "KK_Sm3A", "KK_Sm3B"),
  CCA        = c("POC_CCA_1_R1", "POC_CCA_1_R2", "POC_CCA_1_R3", "POC_CCA_2", "POC_CCA_3"),
  CCA_batch  = c("Batch_POC_CCA_urine_1R1", "Batch_POC_CCA_urine_1R2",
                 "Batch_POC_CCA_urine_1R3", "Batch_POC_CCA_urine_2",
                 "Batch_POC_CCA_urine_3"),
  CCA3       = c("recPOC_1_R1", "recPOC_1_R2", "recPOC_1_R3", "recPOC_2", "recPOC_3"),
  CCA3_batch = c("Batch_recPOC_urine_1R1", "Batch_recPOC_urine_1R2",
                 "Batch_recPOC_urine_1R3", "Batch_recPOC_urine_2",
                 "Batch_recPOC_urine_3")
)

## The Mayuge POC-CCA batch columns hold lot numbers rather than 1/2/3.
LOT_TO_BATCH <- c("220701075" = 1, "220902098" = 2, "221117133" = 3)

## Columns never carried into analysis objects (re-identification risk).
DROP_COLS <- c("Initials")

## ---- reading -------------------------------------------------------------------
read_site <- function(site) {
  f <- file.path(PATHS$raw, SITES[[site]])
  if (!file.exists(f)) {
    stop("Field data not found: ", f,
         call. = FALSE)
  }
  ## The exports carry a UTF-8 byte-order mark; without fileEncoding the first
  ## column is read as `X...ID` / `ï..ID` on some platforms.
  dt <- read.csv(f, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  if (grepl("ID$", names(dt)[1])) names(dt)[1] <- "ID"

  missing <- setdiff(c("ID", unlist(COLS)), names(dt))
  if (length(missing)) stop("Missing columns in ", f, ": ",
                            paste(missing, collapse = ", "), call. = FALSE)

  dt <- dt[!is.na(dt$ID) & nzchar(trimws(as.character(dt$ID))), ]
  dt[, intersect(DROP_COLS, names(dt))] <- NULL
  rownames(dt) <- NULL
  dt
}

recode_batch <- function(m) {
  m <- as.matrix(m)
  lot <- !is.na(m) & as.character(m) %in% names(LOT_TO_BATCH)
  m[lot] <- LOT_TO_BATCH[as.character(m[lot])]
  storage.mode(m) <- "numeric"
  if (any(!is.na(m) & !(m %in% 1:3))) stop("Unrecognised batch code", call. = FALSE)
  m[is.na(m)] <- 4                     # dummy precision for missing tests
  m
}

## Data list for the JAGS models. Values are on the 0-9 G-score scale.
model_data <- function(site) {
  dt <- read_site(site)
  KK      <- as.matrix(dt[, COLS$KK])
  CCA     <- as.matrix(dt[, COLS$CCA])  - 1
  REC_CCA <- as.matrix(dt[, COLS$CCA3]) - 1
  list(
    N             = nrow(dt),
    nKK           = ncol(KK),
    KK            = unname(KK),
    CCA           = unname(CCA),
    REC_CCA       = unname(REC_CCA),
    CCA_batch     = unname(recode_batch(dt[, COLS$CCA_batch])),
    REC_CCA_batch = unname(recode_batch(dt[, COLS$CCA3_batch]))
  )
}

## ---- prior for the logistic G-score curve -------------------------------------
## Multivariate normal fitted to the (slope, midpoint) pairs in kints.csv
## (G-score model from the earlier POC-CCA reproducibility study).
logistic_prior <- function() {
  kints <- read.csv(file.path(PATHS$priors, "kints.csv"))
  fit <- MGMM::FitGMM(data = as.matrix(kints[, c("k_GScore", "intercept_Gscore")]))
  list(multim = as.numeric(fit@Mean), multicov = unname(as.matrix(fit@Covariance)))
}

## ---- posterior helpers ---------------------------------------------------------
posterior_file <- function(model, site) file.path(PATHS$post, sprintf("%s_%s.rds", model, site))

read_posterior <- function(model, site) {
  f <- posterior_file(model, site)
  if (!file.exists(f)) stop("Posterior not found: ", f,
                            "\nRun R/02_field_models.R first.", call. = FALSE)
  readRDS(f)
}

## Stack chains and keep parameters by NAME (never by column position).
draws_matrix <- function(mcmc_list) do.call(rbind, lapply(mcmc_list, as.matrix))

fmt_ci <- function(mid, lo, hi, digits = 1) {
  sprintf(paste0("%.", digits, "f%% (%.", digits, "f–%.", digits, "f)"),
          100 * mid, 100 * lo, 100 * hi)
}
