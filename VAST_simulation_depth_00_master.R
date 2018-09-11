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
clusters <- 3

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
  # Empirical data is used to generate samples
  "Nyears" = NULL, "Nsamp_per_year" = NULL,
  "Depth1_km" = 0, "Depth1_km2" = 0,
  "Dist_sqrtkm" = 0,
  "SigmaO1" = 0.5, "SigmaO2" = 0.5,
  "SigmaE1" = 0.5, "SigmaE2" = 0.5,
  "SigmaV1" = 0, "SigmaV2" = 0,
  "SigmaVY1" = 0, "SigmaVY2" = 0,
  "Range1" = 1000, "Range2" = 500, "SigmaM" = 1,
  "ObsModelcondition" = c(2, 0),
  "ObsModelEM" = c(2, 0),
  "nknots" = n_x,
  "strata" = strata.limits,
  "depth" = "no",
  "Species" = Species,
  "version" = Version,
  "changepar" = c(
    "SigmaO1", "SigmaO2", "SigmaE1", "SigmaE2",
    "Range1", "Range2", "Depth1_km", "Depth2_km"),
  "replicates" = 1:in_rep, "replicatesneeded" = in_rep)

if ("package:VASTWestCoast" %in% search()) {
  detach("package:VASTWestCoast", unload = TRUE)
}
devtools::build("d:/stockAssessment/VASTWestCoast")
install.packages("d:/stockAssessment/VASTWestCoast_0.1.tar.gz",
  type = "source")
library(VASTWestCoast)

RootDir <- file.path(dirname(RootDir), "VAST_simulation_new")
Sim_Settings$depth <- "linear"
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters, conditiondir = file.path(RootDir, "01_VAST_simulation"))
Sim_Settings$depth <- "squared"
Sim_Settings$changepar <- c(Sim_Settings$changepar,
  "Depth1_km2", "Depth2_km2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters)
Sim_Settings$depth <- "no"
Sim_Settings$changepar <- Sim_Settings$changepar[!grepl("Depth",
  Sim_Settings$changepar)]
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters)
Sim_Settings$depth <- "squared"
Sim_Settings$changepar <- c(Sim_Settings$changepar,
  "Depth1_km", "Depth2_km", "Depth1_km2", "Depth2_km2")
Sim_Settings$nknots <- 500
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters)

# Sensitivity
Sim_Settings$nknots <- n_x
Sim_Settings$depth <- "squared"
Sim_Settings$changepar <- c(
  "Depth1_km", "Depth2_km", "Depth1_km2", "Depth2_km2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters, conditiondir = file.path(RootDir, "02_VAST_simulation"))

Sim_Settings$changepar <- c(
  "Depth1_km", "Depth2_km", "Depth1_km2", "Depth2_km2",
  "SigmaO1", "SigmaO2", "SigmaE1", "SigmaE2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters, conditiondir = file.path(RootDir, "02_VAST_simulation"))

Sim_Settings$changepar <- c(
  "Depth1_km", "Depth2_km", "Depth1_km2", "Depth2_km2",
  "Range1", "Range2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters, conditiondir = file.path(RootDir, "02_VAST_simulation"))

Sim_Settings$depth <- "linear"
Sim_Settings$nknots <- n_x
Sim_Settings$changepar <- c(
  "Depth1_km", "Depth2_km")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters, conditiondir = file.path(RootDir, "01_VAST_simulation"))

Sim_Settings$changepar <- c(
  "Depth1_km", "Depth2_km", ,
  "SigmaO1", "SigmaO2", "SigmaE1", "SigmaE2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters, conditiondir = file.path(RootDir, "01_VAST_simulation"))

Sim_Settings$changepar <- c(
  "Depth1_km", "Depth2_km", ,
  "Range1", "Range2")
VAST_simulation(maindir = RootDir, globalsettings = Sim_Settings,
  n_cluster = clusters, conditiondir = file.path(RootDir, "01_VAST_simulation"))

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
