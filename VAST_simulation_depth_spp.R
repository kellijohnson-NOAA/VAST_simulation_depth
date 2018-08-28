###############################################################################
###############################################################################
## Purpose:    VAST simulation for covariates
## Author:     Kelli Faye Johnson
## Contact:    kelli.johnson@noaa.gov
## Date:       2018-03-26
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
for (yesnodepth in c(1, 3)[1]){
  DateDir <- file.path(RootDir, "spp",
    paste0(Speciesloop, "_", n_x, "_",
      ifelse(yesnodepth == 1, "nodepth", "")))
  dir.create(DateDir, recursive = TRUE, showWarnings = FALSE)

Sim_Settings <- list(
  "Species" = Speciesloop,
  "ObsModelcondition" = c(2, 1),
  "nknots" = n_x,
  "strata" = data.frame("STRATA" = "All_areas"),
  "depth" = c("no", "linear", "squared")[yesnodepth],
  "version" = Version,
  "Passcondition" = FALSE)

test <- VAST_condition(
  conditiondir = DateDir,
  settings = Sim_Settings, spp = Sim_Settings$Species,
  datadir = file.path(RootDir, "downloads"),
  overdisperion = NULL)


# indexd <- read.table(dir(dir(DateDir,
#   pattern = "delta_ob", include.dirs = TRUE, full.names = TRUE),
#   "Table", full.names = TRUE), sep = ",", header = TRUE)
# indexn <- read.table(dir(dir(DateDir,
#   pattern = "nodepth", include.dirs = TRUE, full.names = TRUE),
#   "Table", full.names = TRUE), sep = ",", header = TRUE)

# gg <- ggplot() +
# geom_ribbon(data = indexn, aes(x = Year,
#   ymin = log(Estimate_metric_tons) - 1.96*SD_log,
#   ymax = log(Estimate_metric_tons) + 1.96*SD_log),
#   fill = countryColors$color[1], alpha = 0.3) +
# geom_ribbon(data = indexd, aes(x = Year,
#   ymin = log(Estimate_metric_tons) - 1.96*SD_log,
#   ymax = log(Estimate_metric_tons) + 1.96*SD_log),
#   fill = countryColors$color[80], alpha = 0.3) +
# geom_line(data = indexd,
#   aes(x = Year, y = log(Estimate_metric_tons)),
#   col = countryColors$color[1]) +
# geom_line(data = indexn,
#   aes(x = Year, y = log(Estimate_metric_tons)),
#   col = countryColors$color[80]) +
#   ylab("ln abundance (mt)") +
#   theme_bw() + theme(
#   strip.background = element_blank(),
#   panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   scale_x_continuous(breaks = indexn$Year)

# ggsave(paste0(DateDir, basename(DateDir), "_index_log.jpeg"),
#   gg, dev = "jpeg")
}}}

dirs <- dir(file.path(RootDir, "spp"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "Save.RData")
dirs <- dirs[!grepl("DatabaseSave", dirs)]

allaic <- do.call("rbind",
  mapply(VASTWestCoast:::getcompare, dirs, SIMPLIFY = FALSE))
rownames(allaic) <- NULL
colnames(allaic)[1:5] <- c("region", "species", "n", "AIC", "depth")
allaic <- as.data.frame(allaic, stringsAsFactors = FALSE)
allaicwide <- reshape(allaic,
  direction = "wide", timevar = "depth", idvar = c("region", "species", "n"))
rownames(allaicwide) <- NULL
allaicwide$delaic <- as.numeric(allaicwide$"AIC.depth") -
  as.numeric(allaicwide$"AIC.nodepth")
  # todo: fix the next lines
allaicwide <- allaicwide[, -grep("FALSE|^AIC", colnames(allaicwide))]
colnames(allaicwide) <- gsub("TRUE|\\.", "", colnames(allaicwide))
write.table(format(allaicwide, digits = 2), sep = ",", row.names = FALSE,
  file = file.path(RootDir, "VAST_simulation_depth_SpeciesComparison.csv"))
