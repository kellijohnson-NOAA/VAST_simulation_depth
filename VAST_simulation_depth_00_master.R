###############################################################################
###############################################################################
## Purpose:    VAST simulation for covariates
## Author:     Kelli Faye Johnson
## Contact:    kelli.johnson@noaa.gov
## Date:       2018-03-26
## Comments:   Adapted from VAST_simulation_2018-03-25.R (JTT)
###############################################################################
###############################################################################

###############
# Settings
###############
in_rep <- 100
Date <- "depthsquared"
Species <- "WCGBTS_Sebastes_crameri"

###############
# Use inputs
###############
# Find the database directory titled "stockAssessment" on a drive
ccc <- 27; RootDir <- NULL
while(length(RootDir) == 0 & ccc > -1) {
  ccc <- ccc - 1
  RootDir <- file.path(dir(paste0(letters[ccc], ":", .Platform$file.sep),
    pattern = "StockAssessment",
    ignore.case = TRUE, full.names = TRUE), "VAST_simulation_depth")
}
# Load base-case conditions and a bunch of logistical stuff
source(dir(RootDir, pattern = "_00.R", full.names = TRUE))
Sim_Settings <- list(
  "beta1_mean" = 0, "beta2_mean" = 0,
  # The slope remains zero and is not changed.
  "beta1_slope" = 0, "beta2_slope" = 0,
  "beta1_sd" = 0, "beta2_sd" = 0,
  "Nyears" = 10, "Nsamp_per_year" = 600,
  "Depth1_km" = 0, "Depth1_km2" = 0,
  "Dist_sqrtkm" = 0,
  "SigmaO1" = 0.5, "SigmaO2" = 0.5,
  "SigmaE1" = 0.5, "SigmaE2" = 0.5,
  "SigmaV1" = 0, "SigmaV2" = 0,
  "SigmaVY1" = 0, "SigmaVY2" = 0,
  "Range1" = 1000, "Range2" = 500, "SigmaM" = 1,
  "ObsModel" = c(2, 0),
  "ObsModelcondition" = c(2,1),
  "ObsModelEM" = c(2, 0),
  "nknots" = n_x,
  "strata" = strata.limits,
  "depth" = "no",
  "Species" = Species,
  "version" = Version,
  "changepar" = c(
    "SigmaO1", "SigmaO2", "SigmaE1", "SigmaE2",
    "Range1", "Range2"),
  "replicates" = 1:in_rep, "replicatesneeded" = in_rep)

devtools::build("d:/stockAssessment/VASTWestCoast")
install.packages("d:/stockAssessment/VASTWestCoast_0.1.tar.gz",
  type = "source")
library(VASTWestCoast)
# detach("package:VASTWestCoast", unload = TRUE)

Sim_Settings$depth <- "linear"
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2)
Sim_Settings$depth <- "squared"
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2)
Sim_Settings$depth <- "no"
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2)
Sim_Settings$depth <- "squared"
Sim_Settings$nknots <- 500
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2)

# Sensitivity
Sim_Settings$nknots <- n_x
Sim_Settings$changepar <- c("none")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2, conditiondir = file.path(RootDir, "02_VAST_simulation"))

Sim_Settings$changepar <- c("SigmaO1", "SigmaO2", "SigmaE1", "SigmaE2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2, conditiondir = file.path(RootDir, "02_VAST_simulation"))

Sim_Settings$changepar <- c("Range1", "Range2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2, conditiondir = file.path(RootDir, "02_VAST_simulation"))

Sim_Settings$depth <- "linear"
Sim_Settings$nknots <- n_x
Sim_Settings$changepar <- c("none")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2, conditiondir = file.path(RootDir, "01_VAST_simulation"))

Sim_Settings$changepar <- c("SigmaO1", "SigmaO2", "SigmaE1", "SigmaE2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2, conditiondir = file.path(RootDir, "01_VAST_simulation"))

Sim_Settings$changepar <- c("Range1", "Range2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = 2, conditiondir = file.path(RootDir, "01_VAST_simulation"))

# Delete old results files if they exist and read in the results
sapply(dir(RootDir, pattern = "results.csv", recursive = TRUE, full.names = TRUE),
  unlink)
results <- lapply(
  dir(RootDir, pattern = "^[[:digit:]]{2}_", full.names = TRUE),
  get_results)
results <- do.call("rbind", results)

keep <- results[
    results$hessian == TRUE &
    !is.na(results$hessian) &
    results$gradient < 0.0001, ]


#### EndOfFile
