###############################################################################
###############################################################################
## Purpose:    VAST simulation for covariates with several species and knots
## Author:     Kelli Faye Johnson
## Contact:    kelli.johnson@noaa.gov
## Date:       2018-09-19
## Comments:   Adapted from VAST_simulation_2018-03-25.R (JTT)
###############################################################################
###############################################################################

# Database directory
Date <- "spp"
Species <- "WCGBTS_Sebastes_crameri"
ccc <- 27; RootDir <- NULL
while(length(RootDir) == 0) {
  ccc <- ccc - 1
  RootDir <- file.path(dir(paste0(letters[ccc], ":/"),
    pattern = "StockAssessment",
    ignore.case = TRUE, full.names = TRUE), "VAST_simulation_depth")
}
source(dir(RootDir, pattern = "_00.R", full.names = TRUE))

# SpeciesList
SpeciesList <- c(
  "WCGBTS_Atheresthes_stomias", # Arrowtooth 2017
  "WCGBTS_Sebastes_alutus", # POP 2017
  "WCGBTS_Sebastes_pinniger", # Canary rockfish 2015
  "WCGBTS_Sebastes_goodei", # Chilipepper rockfish 2015
  "WCGBTS_Eopsetta_jordani", # Petrale 2015
  "WCGBTS_Anoplopoma_fimbria", # Sablefish 2015
  "WCGBTS_Sebastes_entomelas", # Widow rockfish 2015
  "WCGBTS_Sebastolobus_altivelis", # Longspine Thornyhead 2013
  "WCGBTS_Sebastolobus_alascanus", # Shortspine Thornyhead 2013
  "WCGBTS_Citharichthys_sordidus", # Pacific sanddab 2013
  "WCGBTS_Sebastes_aurora", # Aurora rockfish 2013
  "WCGBTS_Microstomus_pacificus", # Dover sole 2011
  "WCGBTS_Squalus_suckleyi", # Spiny dogfish 2011
  "WCGBTS_Sebastes_crameri", # Darkblotched 2017
  "EBSBTS_Sebastes_alutus", # EBS pop
  "EBSBTS_Atheresthes_stomias" # EBS Arrowtooth
  )

for (Speciesloop in SpeciesList) {
for (n_x in c(125, 250, 500, 750)) {
for (yesnodepth in c("no", "linear", "squared")){
  DateDir <- file.path(RootDir, "spp",
    paste0(Speciesloop, "_", n_x, "_", yesnodepth))
  dir.create(DateDir, recursive = TRUE, showWarnings = FALSE)

spp_Settings <- list(
  "Species" = Speciesloop,
  "ObsModelcondition" = c(2, 0),
  "nknots" = n_x,
  "strata" = data.frame("STRATA" = "All_areas"),
  "depth" = yesnodepth,
  "version" = Version,
  "Passcondition" = FALSE)

test <- VAST_condition(
  conditiondir = DateDir,
  settings = spp_Settings, spp = spp_Settings$Species,
  datadir = file.path(RootDir, "downloads"),
  overdisperion = NULL)
}}}

dirs <- dir(file.path(RootDir, "spp"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "Save.RData")
dirs <- dirs[!grepl("DatabaseSave", dirs)]

allaic <- do.call("rbind",
  mapply(VASTWestCoast:::getcompare, dirs, SIMPLIFY = FALSE))
rownames(allaic) <- NULL
colnames(allaic) <-
  c("region", "species", "n", "depth", "AIC", "d1", "d12", "d2", "d22")
allaic <- as.data.frame(allaic, stringsAsFactors = FALSE)
allaic$depth <- ifelse(allaic$d1 == 0 & allaic$d12 == 0, "no", "linear")
allaic$depth <- ifelse(allaic$d12 != 0 & allaic$d22 != 0,
  "quadratic", allaic$depth)
allaicwide <- reshape(allaic,
  direction = "wide",
  timevar = "depth", idvar = c("region", "species", "n"))
rownames(allaicwide) <- NULL
allaicwide$delaicL <- format(as.numeric(allaicwide$"AIC.no") -
  as.numeric(allaicwide$"AIC.linear"), digits = 2, nsmall = 2)
allaicwide$delaicQ <- format(as.numeric(allaicwide$"AIC.no") -
  as.numeric(allaicwide$"AIC.quadratic"), digits = 2, nsmall = 2)
cols <- c("species", "n", 
  grep("d.+quadratic", colnames(allaicwide), value = TRUE),
  grep("AIC", colnames(allaicwide), value = TRUE, ignore.case = TRUE))
allaicwide <- allaicwide[, cols]
write.table(allaicwide, sep = ",", row.names = FALSE,
  file = file.path(RootDir, "VAST_simulation_depth_SpeciesComparison.csv"))
