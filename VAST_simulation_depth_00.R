
#### Libraries
#### Libraries
if (!"FishData" %in% installed.packages()) devtools::install_github("james-thorson/FishData")
if (!"INLA" %in% installed.packages()) source("http://www.math.ntnu.no/inla/givemeINLA.R")
if (!"SpatialDeltaGLMM" %in% installed.packages()) devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
if (!"TMB" %in% installed.packages()) devtools::install_github("kaskr/adcomp/TMB")
if (!"TMBhelper" %in% installed.packages()) devtools::install_github("kaskr/TMB_contrib_R/TMBhelper") # Optimize
if (!"VAST" %in% installed.packages()) devtools::install_github("james-thorson/VAST")

# Load libraries
library(doParallel, warn.conflicts = TRUE, quietly = TRUE)
library(foreach, warn.conflicts = TRUE, quietly = TRUE)
library(TMB, warn.conflicts = TRUE, quietly = TRUE)
library(TMBhelper, warn.conflicts = TRUE, quietly = TRUE)
library(VAST, warn.conflicts = TRUE, quietly = TRUE)

library(ggplot2, warn.conflicts = TRUE, quietly = TRUE)
library(ggmap, warn.conflicts = TRUE, quietly = TRUE)
library(ggsn, warn.conflicts = TRUE, quietly = TRUE)
library(maps, warn.conflicts = TRUE, quietly = TRUE)
library(mapdata, warn.conflicts = TRUE, quietly = TRUE)
library(MBA, warn.conflicts = TRUE, quietly = TRUE)
library(plyr, warn.conflicts = TRUE, quietly = TRUE)

options(sciepen = 999)

# Set up the parallel processors
# Set up parallel running
if (FALSE) {
  numcores <- Sys.getenv("NUMBER_OF_PROCESSORS")
  numcores <- 3
  mode(numcores) <- "numeric"
  cl <- makeCluster(numcores)
  registerDoParallel(cl)
}
# Instructions: foreach for loop
# foreach(it = dir()) %dopar% {
#   dir(it)
# }

  DownloadDir <- file.path(RootDir, "downloads")
  dir.create(DownloadDir, recursive = TRUE, showWarnings = FALSE)
  source(dir(RootDir, pattern = "functions.R", full.names = TRUE))


  Version <- "VAST_v4_0_0"
  Method <- c("Grid", "Mesh")[2]
  grid_size_km <- 25
