################################################################################
## tools/import_legacy_posteriors.R -- one-off conversion of the original
## workspace images (save.image) into parameter-only posterior files.
##
## The original JAGS scripts saved their whole workspace, including the field
## data. This extracts only the mcmc.list of monitored parameters so the
## downstream pipeline can run without the field data and without the JAGS
## re-fit. Usage:
##   Rscript tools/import_legacy_posteriors.R "<folder with the .RData files>"
################################################################################
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(runjags); library(coda) })

args <- commandArgs(trailingOnly = TRUE)
legacy_dir <- if (length(args)) args[1] else "."

legacy <- c(batch_Tororo = "Results_tor_batch.RData",
            batch_Mayuge = "Results_may_batch.RData",
            inter_Tororo = "Results_inter_t.RData",
            inter_Mayuge = "Results_inter_m.RData")

for (nm in names(legacy)) {
  e <- new.env()
  load(file.path(legacy_dir, legacy[[nm]]), envir = e)
  post <- as.mcmc.list(e$Results)
  parts <- strsplit(nm, "_")[[1]]
  attr(post, "fit_info") <- list(model = parts[1], site = parts[2], N = e$N,
                                 source = legacy[[nm]],
                                 note = "Imported from original run (not re-fitted)")
  saveRDS(post, posterior_file(parts[1], parts[2]))
  message(sprintf("%-14s N = %d  iterations = %d x %d chains", nm, e$N,
                  coda::niter(post), coda::nchain(post)))
}
