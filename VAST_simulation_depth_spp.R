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
ccc <- 27; RootDir <- NULL
while(length(RootDir) == 0) {
  ccc <- ccc - 1
  RootDir <- file.path(dir(paste0(letters[ccc], ":/"),
    pattern = "StockAssessment",
    ignore.case = TRUE, full.names = TRUE), "VAST_simulation_depth")
}
source(dir(RootDir, pattern = "_00.R", full.names = TRUE))
gdURL <- "http://www.stat.ubc.ca/~jenny/notOcto/STAT545A/examples/gapminder/data/gapminderCountryColors.txt"
countryColors <- read.delim(file = gdURL, as.is = 3) # protect color

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
  "WCGBTS_Sebastes_ruberrimus", # Yelloweye rockfish 2011
  "WCGBTS_Sebastes_crameri", # Darkblotched 2017
  "EBSBTS_Anoplopma_fimbria", # EBS Sablefish
  "EBSBTS_Atheresthes_stomias" # EBS Arrowtooth
  )

RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
OverdispersionConfig = c("eta1"=0, "eta2"=0)
Kmeans_Config = list(
  "randomseed" = 1,
  "nstart" = 100,
  "iter.max" = 1e3)     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
Options = c(
  SD_site_density = 0, SD_site_logdensity = 0,
  Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
  Calculate_Cov_SE = 0, Calculate_Synchrony = 0, Calculate_Coherence = 0)
# Gamma p_1
# "logit-link" == 0 vs. "log-link" == 1
ObsModel_Set = list("Conventional delta"=c(2,0), "Poisson-process link"=c(2,1) )
# 1=Presence-absence; 2=Density given presence;
# Epsilon=Spatio-temporal; Omega=Spatial
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
# SpeciesList

# Download the data
initialDataDownload <- FishData::download_catch_rates(
  survey = "WCGBTS",
  species_set = "Sebastes crameri",
  error_tol = 0.01, localdir = paste0(DownloadDir, "/"))

for (Species in SpeciesList) {
for (n_x in c(100, 250, 500, 1000)) {
	# n_x <- ifelse(grepl("^WC", Species), 250, 100)  # Number of stations
  Date <- paste0("VAST_simulation_depth_", Species, "_", n_x)
  DateDir <- file.path(RootDir, Date, "/")
  dir.create(DateDir, recursive = TRUE, showWarnings = FALSE)
  if (grepl("^WC", Species)) OverdispersionConfig = c("Delta1"=1, "Delta2"=1)

  # Decide on case-specific settings for use when calculating indices
  # strata.limits <- data.frame(
  #   "STRATA" = c("Coastwide","CA","OR","WA"),
  #   "north_border" = c(49.0, 42.0, 46.0, 49.0),
  #   "south_border" = c(32.0, 32.0, 42.0, 46.0),
  #   "shallow_border" = c(55, 55, 55, 55),
  #   "deep_border" = c(1280, 1280, 1280, 1280)
  #   )
  strata.limits <- data.frame('STRATA'="All_areas")

  # Save
  Record = list("Species"=Species, "Version"=Version,
  	"Method"=Method,
  	"grid_size_km"=grid_size_km,
  	"n_x"=n_x,
  	"FieldConfig"=FieldConfig,
  	"RhoConfig"=RhoConfig,
  	"OverdispersionConfig"=OverdispersionConfig,
  	"Kmeans_Config"=Kmeans_Config,
  	"Options"=Options,
  	"ObsModel_Set"=ObsModel_Set,
  	"strata.limits"=strata.limits)
  save(Record, file=paste0(DateDir,"Record.RData"))


###############
# Run multiple operating models
###############

repI = omI = emI = 1
observed = "_observed"
for(omdepth in c("", "_nodepth")) {

  # Create directory for OM
  OmDir = paste0(DateDir,
    "OM=", names(ObsModel_Set)[omI], omdepth, observed, "/")
  RepDir = paste0(OmDir,"rep=",repI,"/")
  EmDir = paste0(RepDir,"EM=",names(ObsModel_Set)[emI],"/TRUE/")
    dir.create(EmDir, recursive=TRUE)

    Database = FishData::download_catch_rates(
      survey = strsplit(Species, "_")[[1]][1],
      species_set = gsub("_", " ", gsub("[A-Z]{3}BTS_", "", Species)),
      # species_set = 25,
      error_tol = 0.01, localdir = paste0(DownloadDir, "/"))
    Database$Vessel <- as.factor(
      paste(Database$Vessel, Database$Year, sep = "_"))

    if(!("Vessel" %in% names(Database))) Database = cbind(Database, 'Vessel'=1 )
    Database = ThorsonUtilities::rename_columns( Database[,
      c('Sci','Wt','Year','Long','Lat', "Vessel")], newname=
      c('Sci','Catch_KG','Year','Lon','Lat', "Vessel"))
    Database = cbind(Database,
      'AreaSwept_km2' = 0.01)  # WCGBTS and all AFSC surveys are in KG/Hectare

    Database = na.omit( Database )
    # tapply( ifelse(Database$Catch_KG>0,1,0), INDEX=Database$Year, FUN=mean )

    # Save
    DatabaseSave = list("Database"=Database)
    save( DatabaseSave, file=paste0(DateDir,"DatabaseSave.RData"))
  Data_Geostat = DatabaseSave$Database

  # Get extrapolation data
  if(grepl("EBSBTS", Species)){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="eastern_bering_sea", strata.limits=strata.limits )
    Extrapolation_List$Data_Extrap = ThorsonUtilities::rename_columns( Extrapolation_List$Data_Extrap, origname="Mid_Depth", newname="Depth_km")
  }
stand01 <- function(data) {
  (data - mean(data)) / sd(data)
}
  if(grepl("WCGBTS", Species)){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="California_current", strata.limits=strata.limits )
    Extrapolation_List_keep <-
      SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(
        Region = "California_current", strata.limits = strata.limits)
  }
  # Standardize the depth data
  Extrapolation_List$Data_Extrap[,'Depth_km'] <-
    (Extrapolation_List$Data_Extrap[,'Depth_km'] -
    mean(Extrapolation_List$Data_Extrap[,'Depth_km'])) /
    sd(Extrapolation_List$Data_Extrap[,'Depth_km'])


  # Calculate spatial information for SPDE mesh, strata areas, and AR1 process
  Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateDir )
  Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

  # Generate covariate
  Depth_x = tapply( Extrapolation_List$Data_Extrap[,'Depth_km'], INDEX=Spatial_List$PolygonList$NN_Extrap$nn.idx, FUN=mean )
  X_xtp = Depth_x %o% rep(1,diff(range(Data_Geostat[,'Year']))+1) %o% 1
  X_xtp <- array(c(X_xtp, X_xtp^2), dim = c(dim(X_xtp)[1:2], 2))
  Data_raw <- list()
  Data_raw$X_xtp <- X_xtp
  # Change depth in Extrapolation_List so that it matches EM
  Extrapolation_List$Data_Extrap[,'Depth_km'] = Depth_x[Spatial_List$PolygonList$NN_Extrap$nn.idx]

  # Plot settings
  Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

    OmDir = paste0(DateDir,
      "OM=", names(ObsModel_Set)[omI], omdepth, observed, "/")
    Save <- list("Data_Geostat" = Data_Geostat, "Spatial_raw" = Spatial_List,
      "Data_raw" = Data_raw)
    Sim <- list("Data_Geostat" = Data_Geostat)
    save(Save, file = paste0(OmDir, "Save.RData"))
    save(Sim, file = paste0(OmDir, "rep=", repI, "/Sim.RData"))

    finalEM <- rerun(Date = gsub("VAST_simulation_depth_", "", Date),
        RootDir = RootDir,
        OM = paste0("OM=", names(ObsModel_Set)[omI], omdepth, observed),
        EM = paste0("OM=", names(ObsModel_Set)[omI], omdepth, observed),
        rep = 1,
        depth = ifelse(omdepth == "", TRUE, FALSE))
    }

indexd <- read.table(dir(dir(DateDir,
  pattern = "delta_ob", include.dirs = TRUE, full.names = TRUE),
  "Table", full.names = TRUE), sep = ",", header = TRUE)
indexn <- read.table(dir(dir(DateDir,
  pattern = "nodepth", include.dirs = TRUE, full.names = TRUE),
  "Table", full.names = TRUE), sep = ",", header = TRUE)

gg <- ggplot() +
geom_ribbon(data = indexn, aes(x = Year,
  ymin = log(Estimate_metric_tons) - 1.96*SD_log,
  ymax = log(Estimate_metric_tons) + 1.96*SD_log),
  fill = countryColors$color[1], alpha = 0.3) +
geom_ribbon(data = indexd, aes(x = Year,
  ymin = log(Estimate_metric_tons) - 1.96*SD_log,
  ymax = log(Estimate_metric_tons) + 1.96*SD_log),
  fill = countryColors$color[80], alpha = 0.3) +
geom_line(data = indexd,
  aes(x = Year, y = log(Estimate_metric_tons)),
  col = countryColors$color[1]) +
geom_line(data = indexn,
  aes(x = Year, y = log(Estimate_metric_tons)),
  col = countryColors$color[80]) +
  ylab("ln abundance (mt)") +
  theme_bw() + theme(
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_x_continuous(breaks = indexn$Year)

ggsave(paste0(DateDir, basename(DateDir), "_index_log.jpeg"),
  gg, dev = "jpeg")
}}

dirs <- sapply(dir(RootDir, pattern = "BTS_", full.names = TRUE),
  dir, recursive = TRUE, pattern = "Save.RData")
dirs <- dirs[!grepl("DatabaseSave", dirs)]

getcompare <- function(file) {
  load(file)
  # browser()
  aic <- AIC
  temp <- gsub("VAST_simulation_depth_|/Save.RData|EBSBTS_|WCGBTS_", "", file)
  temp <- strsplit(temp, "/")[[1]]
  species <- gsub("_", " ", temp[1])
  area <- ifelse(grepl("EBSBTS_", file), "EBS", "WC")
  type <- ifelse(grepl("nodepth", temp[2]), "FALSE", "TRUE")
  depth <- unique(unlist(Save$ParHat[c("gamma1_ctp", "gamma2_ctp")]))
  if (length(par) == 1) par <- c(0, 0)
  returnme <- data.frame(area, species, type, aic, depth[1], depth[2], stringsAsFactors = FALSE)
}
allaic <- do.call("rbind", mapply(getcompare, dirs, SIMPLIFY = FALSE))
allaicwide <- reshape(allaic, direction = "wide",
  idvar = "species", timevar = "type")
rownames(allaicwide) <- NULL
allaicwide <- allaicwide[, !grepl("FALSE|aic", colnames(allaicwide))]
allaicwide$delaic <- tapply(allaic[, 3], allaic[, 1], diff)
colnames(allaicwide) <- gsub("TRUE|\\.", "", colnames(allaicwide))
write.table(format(allaicwide, digits = 2), sep = ",", row.names = FALSE,
  file = file.path(RootDir, "VAST_simulation_depth_SpeciesComparison.csv"))
