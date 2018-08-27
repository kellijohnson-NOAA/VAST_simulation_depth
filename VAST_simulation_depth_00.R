
#### Libraries
if (!"FishData" %in% installed.packages())
  devtools::install_github("james-thorson/FishData")
if (!"INLA" %in% installed.packages())
  source("http://www.math.ntnu.no/inla/givemeINLA.R")
# SpatialDeltaGLMM is now deprecated, but I started with this package and
# so the code references a given commit to make sure that the correct
# version is installed.
if (!"SpatialDeltaGLMM" %in% installed.packages() |
  substring(packageDescription("SpatialDeltaGLMM")$Version, 1, 3) < 3.4)
  devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM",
    ref = "da3a3badfd0c50cc805fc5439149fb3370c3791e")
if (!"TMB" %in% installed.packages())
  devtools::install_github("kaskr/adcomp/TMB")
# Used to optimize the model
if (!"TMBhelper" %in% installed.packages())
  devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
if (!"VAST" %in% installed.packages())
  devtools::install_github("james-thorson/VAST")

# Load libraries
library(doParallel, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(foreach, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(TMB, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(TMBhelper, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(VAST, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)

library(dplyr, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(ggplot2, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(ggmap, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(ggsn, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(maps, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(mapdata, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(MBA, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(plyr, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)
library(tidyr, warn.conflicts = TRUE, quietly = TRUE, verbose = FALSE)

#### Directory
# Make specific folder for this run series
DateDir <- file.path(
  RootDir,
  paste0("VAST_simulation_depth_", Date),
  .Platform$file.sep)
dir.create(DateDir, showWarnings = FALSE)
DownloadDir <- file.path(RootDir, "downloads")
dir.create(DownloadDir, recursive = TRUE, showWarnings = FALSE)
# Download the data
initialDataDownload <- FishData::download_catch_rates(
  survey = "WCGBTS",
  species_set = "Sebastes crameri",
  error_tol = 0.01, localdir = paste0(DownloadDir, .Platform$file.sep))
initialDataDownloadEB <- FishData::download_catch_rates(
  survey = "EBSBTS",
  species_set = "Sebastes alutus",
  error_tol = 0.01, localdir = paste0(DownloadDir, .Platform$file.sep))

#### Options
options(scipen = 999)
VAST_simulation_depth_theme <- theme_bw() + theme(
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1))

#### Model configuration
# Number of stations
strata.limits <- data.frame("STRATA" = "All_areas")
n_x <- ifelse(grepl("^WC", Species), 250, 100)
  if (Date == "depthsquared500") n_x <- 500
Version <- "VAST_v4_0_0"

OverdispersionConfig <- c("eta1"=0, "eta2"=0)
  if (grepl("^WC", Species)) OverdispersionConfig <- c("Delta1"=1, "Delta2"=1)

# Gamma p_1
# "logit-link" == 0 vs. "log-link" == 1
ObsModel_Set <- list(
  "Conventional delta" = c(2, 0),
  "Poisson-process link" = c(2, 1))

#### Plotting
gdURL <- "http://www.stat.ubc.ca/~jenny/notOcto/STAT545A/examples/gapminder/data/gapminderCountryColors.txt"
countryColors <- read.delim(file = gdURL, as.is = 3) # protect color
