################################################################################
## run_all.R -- run the full pipeline in order from a clean checkout
##
##   Rscript run_all.R
##
## Steps that need the field data (data/raw/) are skipped with a message if the
## data are absent. JAGS models are refitted only if their posterior file is
## missing (set REFIT <- TRUE inside 01/02 to force).
################################################################################
options(cca3.run_all = TRUE)
step <- function(f, ...) {
  message("\n==== ", f, " ====")
  env <- new.env()
  list2env(list(...), env)
  sys.source(here::here("R", f), envir = env)
}
have_field_data <- all(file.exists(here::here("data", "raw",
  c("CORNTD_recPOCCCA_2023_Tororo.csv", "CORNTD_recPOCCCA_2023_Mayuge.csv"))))

step("01_lab_models.R")
if (have_field_data) step("02_field_models.R") else message("No field data: using shipped posteriors in outputs/posteriors/")
step("03_diagnostics.R")
step("05_simulations.R", MODEL = "inter")
step("06_fig3_batch_sd.R", MODEL = "inter")
if (have_field_data) step("07_descriptives.R") else message("No field data: skipping 07_descriptives.R")

writeLines(capture.output(sessionInfo()), here::here("outputs", "sessionInfo.txt"))
message("\nDone.")
