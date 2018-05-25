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

# DateDir
Date <- Sys.Date()
  DateDir <- file.path(RootDir, paste0("VAST_simulation_depth_", Date), "/")

###############
# Settings
###############
if(file.exists(paste0(DateDir,"Record.RData"))){
  load(file = paste0(DateDir,"Record.RData"))
  attach(Record)
  message("Loaded Record")
}else{
  # Species <- "WCGBTS_Anoplopoma_fimbria"
  Species <- "WCGBTS_Sebastes_crameri"
  n_x <- ifelse(grepl("^WC", Species), 250, 100)  # Number of stations
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
  OverdispersionConfig = c("eta1"=0, "eta2"=0)
  if (grepl("^WC", Species)) OverdispersionConfig = c("Delta1"=1, "Delta2"=1)
  Kmeans_Config = list(
    "randomseed" = 1,
    "nstart" = 100,
    "iter.max" = 1e3)     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
  Options = c(
    SD_site_density = 0, SD_site_logdensity = 0,
    Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
    Calculate_Cov_SE = 0, Calculate_Synchrony = 0, Calculate_Coherence = 0)
  Sim_Settings = list(
    "beta1_mean" = 0, "beta2_mean" = 0,
    "beta1_slope" = 0, "beta2_slope" = 0,
    "beta1_sd" = 0, "beta2_sd" = 0,
    "Nyears" = 10, "Nsamp_per_year" = 600,
    "Depth_km" = 0, "Depth_km2" = 0, "Dist_sqrtkm" = 0,
    "SigmaO1" = 0, "SigmaO2" = 0,
    "SigmaE1" = 0, "SigmaE2" = 0,
    "SigmaV1" = 0, "SigmaV2" = 0,
    "SigmaVY1" = 0, "SigmaVY2" = 0,
    "Range1" = 1000, "Range2" = 500, "SigmaM" = 1,
    "ObsModel" = c(2, 0))

  # Gamma p_1
  # "logit-link" == 0 vs. "log-link" == 1
  ObsModel_Set = list( "Conventional delta"=c(2,0), "Poisson-process link"=c(2,1) )
  # 1=Presence-absence; 2=Density given presence;
  # Epsilon=Spatio-temporal; Omega=Spatial
  FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
  n_rep = 100

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
  Record = list("Species"=Species, "Version"=Version, "Method"=Method, "grid_size_km"=grid_size_km, "n_x"=n_x, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig, "Kmeans_Config"=Kmeans_Config, "Options"=Options, "ObsModel_Set"=ObsModel_Set, "n_rep"=n_rep, "strata.limits"=strata.limits, "Sim_Settings"=Sim_Settings)
  save( Record, file=paste0(DateDir,"Record.RData"))
}

###############
# Run multiple operating models
###############

repI = omI = emI = 1
omdepth = "_nodepth"; observed = "_observed"
if (Date > "2018-05-12") rawrange = TRUE

for(observed in c("_observed", "_5")) {
for(omdepth in c("", "_nodepth")) {
for(repI in 1:n_rep) {
for(omI in 1:length(ObsModel_Set[1])) {
for(emI in 1:length(ObsModel_Set[1])) {

  # Create directory for OM
  OmDir = paste0(DateDir,
    "OM=", names(ObsModel_Set)[omI], omdepth, observed, "/")
  RepDir = paste0(OmDir,"rep=",repI,"/")
  EmDir = paste0(RepDir,"EM=",names(ObsModel_Set)[emI],"/TRUE/")
    dir.create(EmDir, recursive=TRUE)

  if( file.exists(paste0(DateDir,"DatabaseSave.RData")) ){
    load( file=paste0(DateDir,"DatabaseSave.RData"))
  }else{
    # Surveys with a public API  #   # FishData::
    Database = FishData::download_catch_rates(
      survey = strsplit(Species, "_")[[1]][1],
      species_set = gsub("_", " ", gsub("[A-Z]{3}BTS_", "", Species)),
      # species_set = 25,
      error_tol = 0.01, localdir = paste0(DownloadDir, "/"))
    Database$Vessel <- as.factor(
      paste(Database$Vessel, Database$Year, sep = "_"))

    if(!("Vessel" %in% names(Database)) ) Database = cbind( Database, 'Vessel'=1 )
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
  }
  Data_Geostat = DatabaseSave$Database

  # Get extrapolation data
  if(grepl("EBSBTS", Species)){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="eastern_bering_sea", strata.limits=strata.limits )
    Extrapolation_List$Data_Extrap = ThorsonUtilities::rename_columns( Extrapolation_List$Data_Extrap, origname="Mid_Depth", newname="Depth_km")
  }
  if(grepl("WCGBTS", Species)){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="California_current", strata.limits=strata.limits )
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
  # Change depth in Extrapolation_List so that it matches EM
  Extrapolation_List$Data_Extrap[,'Depth_km'] = Depth_x[Spatial_List$PolygonList$NN_Extrap$nn.idx]

  # Plot settings
  Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

  ##############################
  # Calculate parameters for OM
  ##############################

  if( file.exists(paste0(OmDir,"Save.RData")) ){
    load( file=paste0(OmDir,"Save.RData") )
  }else{
    # Make TMB data list
    TmbDataMaster = VAST::Data_Fn(
      "Version" = Version,
      "X_xtp" = X_xtp,
      "FieldConfig" = FieldConfig,
      "OverdispersionConfig"=OverdispersionConfig,
      "RhoConfig" = RhoConfig,
      "ObsModel" = ObsModel_Set[[omI]],
      "c_i" = rep(0,nrow(Data_Geostat)),
      "b_i" = Data_Geostat[,'Catch_KG'],
      "a_i" = Data_Geostat[,'AreaSwept_km2'],
      "v_i" = as.numeric(Data_Geostat[,'Vessel'])-1,
      "s_i" = Data_Geostat[,'knot_i']-1,
      "t_i" = Data_Geostat[,'Year'],
      "a_xl" = Spatial_List$a_xl,
      "MeshList" = Spatial_List$MeshList,
      "GridList" = Spatial_List$GridList,
      "Method" = Spatial_List$Method,
      "Options" = Options)
    TmbData <- TmbDataMaster
    OmDir = paste0(DateDir,
      "OM=", names(ObsModel_Set)[omI], omdepth, observed, "/")
    # Make TMB data list
    if (omdepth != ""){
        TmbData = VAST::Data_Fn(
          "Version" = Version,
          "X_xtp" = NULL,
          "FieldConfig" = FieldConfig,
          "OverdispersionConfig"=OverdispersionConfig,
          "RhoConfig" = RhoConfig,
          "ObsModel" = ObsModel_Set[[omI]],
          "c_i" = rep(0,nrow(Data_Geostat)),
          "b_i" = Data_Geostat[,'Catch_KG'],
          "a_i" = Data_Geostat[,'AreaSwept_km2'],
          "v_i" = as.numeric(Data_Geostat[,'Vessel'])-1,
          "s_i" = Data_Geostat[,'knot_i']-1,
          "t_i" = Data_Geostat[,'Year'],
          "a_xl" = Spatial_List$a_xl,
          "MeshList" = Spatial_List$MeshList,
          "GridList" = Spatial_List$GridList,
          "Method" = Spatial_List$Method,
          "Options" = Options)
          OmDir = paste0(DateDir,
    "OM=", names(ObsModel_Set)[omI], omdepth, observed, "/")
    }
    # Make TMB object
    TmbList = VAST::Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateDir,
      "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x)

    # Run model
    Obj = TmbList[["Obj"]]
    OptOM = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=OmDir, newtonsteps=1, bias.correct = TRUE)
    Report = Obj$report()
    Save = list("Report" = Report, "ParHat" = Obj$env$parList(OptOM$par),
      "OptOM" = OptOM, "Data_raw" = TmbData, "Spatial_raw" = Spatial_List,
      "Version" = Version, "RhoConfig" = RhoConfig)
    save(Save, file=paste0(OmDir,"Save.RData"))

    # Plot index
    Index = SpatialDeltaGLMM::PlotIndex_Fn(DirName=OmDir, TmbData=TmbData,
      Sdreport=OptOM[["SD"]], Year_Set=Year_Set,
      strata_names=strata.limits[,1], use_biascorr=TRUE)
  }

  ##############################
  # Simulate new data
  ##############################

  if( file.exists(paste0(RepDir,"Sim.RData")) ){
    load( file=paste0(RepDir,"Sim.RData") )
  }else{
    # Define parameters
    Sim_Settings[["ObsModel"]] = ObsModel_Set[[omI]]
    Sim_Settings[["beta1_mean"]] = mean(Save$ParHat$beta1_ct)
    Sim_Settings[["beta1_sd"]] = sd(Save$ParHat$beta1_ct)
    Sim_Settings[["beta2_mean"]] = mean(Save$ParHat$beta2_ct)
    Sim_Settings[["beta2_sd"]] = sd(Save$ParHat$beta2_ct)
    Sim_Settings[["SigmaM"]] = exp(Save$ParHat$logSigmaM[1,1])
    Sim_Settings[["Depth1_km"]] = mean(Save$ParHat$gamma1_ctp)
    Sim_Settings[["Depth2_km"]] = mean(Save$ParHat$gamma2_ctp)

    # Define spatial and spatio-temporal variation
    # I'm using lower values than observed so that its less likely to have replicates with 0% or 100% encounter rates
    Sim_Settings[["SigmaO1"]] = 0.5 # abs(Save$ParHat$L_omega1_z)
    Sim_Settings[["SigmaO2"]] = 0.5 # abs(Save$ParHat$L_omega2_z)
    Sim_Settings[["SigmaE1"]] = 0.5 # abs(Save$ParHat$L_epsilon1_z)
    Sim_Settings[["SigmaE2"]] = 0.5 # abs(Save$ParHat$L_epsilon1_z)

    if (observed == "_observed") {
      Sim_Settings[["SigmaO1"]] = abs(Save$ParHat$L_omega1_z)
      Sim_Settings[["SigmaO2"]] = abs(Save$ParHat$L_omega2_z)
      Sim_Settings[["SigmaE1"]] = abs(Save$ParHat$L_epsilon1_z)
      Sim_Settings[["SigmaE2"]] = abs(Save$ParHat$L_epsilon1_z)
    }
    if (rawrange == TRUE) {
      Sim_Settings[["Range1"]] = Save$Report$Range_raw1
      Sim_Settings[["Range2"]] = Save$Report$Range_raw2
    }

    # Check encounter rate
    Prop_t = 0
    counter <- 0
    while(any(any(Prop_t==0) | any(Prop_t==1)) & counter < 10) {
      counter <- counter + 1
      Sim = SpatialDeltaGLMM::Geostat_Sim(Sim_Settings=Sim_Settings,
        Extrapolation_List=Extrapolation_List, Data_Geostat=Data_Geostat )
      Prop_t = tapply(Sim$Data_Geostat[,'Catch_KG'],
        INDEX=Sim$Data_Geostat[,'Year'], FUN=function(vec){mean(vec>0)} )
    }
    save(Sim, file=paste0(RepDir,"Sim.RData"))
    # Plot
    ThorsonUtilities::save_fig(file=paste0(RepDir,"Index-Sim"),
      width=5, height=5, res=200, units='in')
    plot( x=Year_Set, y=Sim$B_tl[,1]/1000, type="b",
      ylim=c(0,max(Sim$B_tl[,1]/1000)) )
    dev.off()
  }

  ##############################
  # Fit using different models
  ##############################

  # Run model if necessary
  if( !file.exists(paste0(EmDir,"parameter_estimates.RData")) ){
    for(emCov in c(TRUE, FALSE)) {
      if (!emCov) {
          EmDir <- gsub("TRUE", "FALSE", EmDir)
          dir.create(EmDir, recursive = TRUE)
      }
      finalEM <- rerun(Date = Date,
        RootDir = RootDir,
        OM = paste0("OM=", names(ObsModel_Set)[omI], omdepth, observed),
        EM = names(ObsModel_Set)[emI], rep = repI,
        depth = emCov)
      # counter <- 0
      # OptEM <- list("SD" = NULL)
      # while (is.null(OptEM[["SD"]]) & counter < 2) {
      #   TmbListEM  = Build_TMB_Fn("TmbData"=TmbDataEM, "RunDir"=DateDir,
      #     "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x)
      #   ObjEM = TmbListEM[["Obj"]]
      #   OptEM = try(TMBhelper::Optimize(obj=ObjEM,
      #     lower=TmbListEM[["Lower"]], upper=TmbListEM[["Upper"]],
      #     getsd=TRUE, savedir=EmDir, newtonsteps=1, bias.correct=TRUE),
      #     silent = TRUE)
      #   if (class(OptEM) == "try-error") OptEM <- list("SD" = NULL)
      # }
    }
  }
}}}}} # omI

keep <- get_results(gsub("/$", "", DateDir))
keep$sigmao1_factor <- ifelse(keep$om_sigmao1 == 0.5, 0.5, "Observed")
keep_re <- keep_num(re_all(data = keep))
keep_ae <- keep_num(re_all(data = keep, type = "ae"))

dim(keep);dim(keep_re);dim(keep_ae)
aggregate(hessian ~ om_type + om_depth + em_depth, data = keep, length)
aggregate(hessian ~ om_type + om_depth + em_depth, data = keep_ae, length)

## Now make running MAE calculations. Takes a while to run these calcs!!!
my.median <- function(x) {sapply(1:length(x), function(i) {median(x[1:i], na.rm=TRUE)})}
test <- plyr::ddply(subset(keep_ae, sigmao1_factor == "Observed"),
  .(om_depth, em_depth),
  .fun=summarize,
  replicate2=1:length(rep),
  MAE_range=my.median(abs(ae_range1)) - median(abs(ae_range1), na.rm = TRUE),
  MAE_endyear=my.median(abs(ae_index_13)) - median(abs(ae_index_13), na.rm = TRUE))
scalars.mares <- reshape::melt(test,
  c("om_depth", "em_depth", "replicate2"))
gg <- ggplot(data=scalars.mares,
  aes(x=replicate2, y=value, col = em_depth))+
    geom_vline(xintercept=100, col= "red") +
    ylab("Centered MAE") +xlab("replicate")
gg <- gg + geom_line(lwd=.5, alpha=1.0)+
  facet_grid(om_depth ~ variable,
  scales = "free") + VAST_simulation_depth_theme +
  scale_y_continuous(limits = c(-320, 800), expand = c(-0.10, -0.25)) +
  guides(colour = guide_legend(title = "Depth in EM"))
ggsave(file.path(DateDir,
  "VAST_simulation_depth_runningmean.jpeg"),
  gg, dev = "jpeg", height = 8, width = 8)

pdf(file.path(DateDir, "VAST_simulation_depth_plots.pdf"))
for (ii_name in grep("re_", colnames(keep_re), value = TRUE)) {
  for (iii_name in grep("re_", colnames(keep_re), value = TRUE)) {
    if (ii_name == iii_name) next
    if (grepl("index", ii_name) | grepl("index", iii_name)) next
    ggplotre(keep_re, ii_name, iii_name, type = "points",
      facety = c("om_type", "om_depth", "om_sigmao1"), gradient = TRUE)
  }
}
dev.off()


###############################################################################
#### Results
###############################################################################

###############################################################################
## Tables
# Calculate the depth range
utils::data(california_current_grid, package = "SpatialDeltaGLMM")
depthsmean <- tapply( california_current_grid[,'Depth_km'],
  INDEX=Spatial_List$PolygonList$NN_Extrap$nn.idx, FUN=mean)
range(depthsmean)*1000
range(depthsmean[unique(Data_Geostat[Data_Geostat$Catch_KG>0,]$knot_i)])*1000
mean(aggregate(Catch_KG ~ Year, data = DatabaseSave$Database,
  function(x) sum(x>0)/length(x))[, 2])

# Info about raw empirical data
mean(with(DatabaseSave$Database,
  tapply(Catch_KG, Year, function(x) sum(x>0)/length(x))))
range(with(DatabaseSave$Database, tapply(Catch_KG, Year, length)))

# Fits to empirical data
rawests <- Reduce(function(x,y) merge(x, y, all = TRUE),
  lapply(dir(dir(DateDir, full.names = TRUE, pattern = "observed"),
    recursive = FALSE, "parameter_estimates.RData", full.names = TRUE),
  getrdata))
rawests <- reshape(rawests,
  direction = "wide", idvar = "par", timevar = "om_name")
colnames(rawests) <- gsub(" delta|_observed|entional", "", colnames(rawests))
rawests[grep("L", rawests$par), ]
rawests[grep("gamma", rawests$par), ]

inpar <- int95(rawests[, 2], rawests[, 4], rawests[, 3], rawests[, 5], rawests$par)
write.table(inpar[inpar[,7]!=1, ],
  file = file.path(DateDir, "parameter_estimates_intervals.txt"))

# AIC of empirical models
rawaic <- mapply(readLines,
  dir(dir(DateDir, full.names = TRUE, pattern = "observed"),
    recursive = FALSE, "parameter_estimates.txt", full.names = TRUE),
  MoreArgs = list(n = 40L))
rawaic <- gsub("\\[1\\]| ", "", rawaic[grep("AIC", rawaic[, 1]) + 1, ])
names(rawaic) <- basename(dirname(names(rawaic)))
diff(as.numeric(rawaic))

###############################################################################
## Plots

#### Map of study area
# download the data from
# https://www.arcgis.com/home/item.html?id=f7f805eb65eb4ab787a0a3e1116ca7e5
url <- 'https://www.arcgis.com/sharing/rest/content/items/f7f805eb65eb4ab787a0a3e1116ca7e5/data'
fil <- "states_21basic.zip"
if (!file.exists(file.path(DownloadDir, fil))) {
  download.file(url, file.path(DownloadDir, fil))
}
  unzip(file.path(DownloadDir, fil), exdir = DownloadDir)

states <- rgdal::readOGR(paste0(DownloadDir, "states.shp"),
  layer = "states", stringsAsFactors = FALSE)
states <- sp::spTransform(states, sp::CRS("+proj=longlat +datum=WGS84"))
states <- ggplot2::fortify(states)

box <- make_bbox(lat = Lat, lon = Lon, data = DatabaseSave$Database)
big <- get_map(location = box, source = "stamen", maptype = "watercolor",
  zoom = 5)
gg <- ggmap(big) +
geom_map(data=states, map = states,aes(x=long,y=lat,map_id=id),
  col = "white", fill = "orangered4", alpha = .2, size = .2) +
geom_point(data = DatabaseSave$Database[DatabaseSave$Database$Catch_KG > 0, ],
  aes(x = Lon, y = Lat, cex = Catch_KG), alpha = 0.2) +
xlab("longitude") + ylab("latitude") + labs(cex="catch (kg)") +
theme(legend.position = c(.2, .1),
  legend.key = element_rect(fill=alpha('white', 0.1)),
  legend.background = element_rect(fill=alpha('blue', 0.01)))

gg <- gg + ggsn::north2(gg, .65, .91)
ggsave(file.path(DateDir,
  "VAST_simulation_depth_map_rawdata.jpeg"),
  gg, dev = "jpeg", height = 8, width = 3.8)

# mba <- MBA::mba.surf(
#   california_current_grid[, c("Lon", "Lat", "Depth_km")], 100, 100)
# dimnames(mba$xyz.est$z) <- list(mba$xyz.est$x, mba$xyz.est$y)
# df3 <- reshape::melt(mba$xyz.est$z,
#   varnames = c("Lon", "Lat"), value.name = "Depth_km")
# ggplot(data=df3, aes(Lon, Lat))+
# geom_raster(aes(fill = value), interpolate = F, hjust = 0.5, vjust = 0.5) +
# geom_contour(aes(z = value))

#### Empirical index plot
fits <- lapply(sapply(dir(DateDir, pattern = "observed", full.names = TRUE),
  dir, full.names = TRUE, pattern = "Table"),
  read.csv, header = TRUE, sep = ",")
fits[[1]]$depth <- "FALSE"
fits[[2]]$depth <- "TRUE"
fits <- do.call("rbind", fits)
fits <- with(fits,
  data.frame(fits,
    fit = log(fits$Estimate_metric_tons),
    lwr = log(fits$Estimate_metric_tons) - 1.96*SD_log,
    upr = log(fits$Estimate_metric_tons) + 1.96*SD_log))
rownames(fits) <- basename(rownames(fits))

# need raw data
gg <- ggplot(fits, aes(Year, fit, color = as.factor(depth),
  fill = as.factor(depth))) +
    geom_ribbon(data=fits,aes(ymin=lwr,ymax=upr),alpha=0.2,
      show.legend = FALSE) +
    geom_line(data=fits, lwd = 2)+
    theme_bw() +
    scale_x_continuous(limits = range(Data_Geostat$Year), expand = c(0, 0.2),
      labels = scaleFUN) +
    scale_colour_grey(name = 'depth', end = 0.5, start = 0.8) +
    scale_fill_grey(end = 0.5, start = 0.8) +
    guides(fill = FALSE) +
    theme(legend.key = element_rect(colour = NA, fill = NA),
        legend.justification = c(1,0), legend.position = c(0.85, .8)) +
    ylab("ln abundance (mt)")
ggsave(file.path(DateDir,
  "VAST_simulation_depth_rawTimeSeries.jpeg"),
  gg, dev = "jpeg")

#### Diagnostic plots
# Get region-specific settings for plots
MapDetails_List <- SpatialDeltaGLMM::MapDetails_Fn(
  Region = ifelse(grep("WC", Species), "California_current", "NA"),
  NN_Extrap = Spatial_List$PolygonList$NN_Extrap,
  Extrapolation_List = Extrapolation_List)

#Plot Anisotropy
tempdir <- dir(DateDir, pattern = "observed", full.names = TRUE)
mapply(dir, tempdir)
e1 <- new.env(parent = baseenv())
e2 <- new.env(parent = baseenv())
load(file.path(tempdir[1], "Save.RData"), envir = e1)
load(file.path(tempdir[2], "Save.RData"), envir = e2)
SpatialDeltaGLMM::PlotAniso_Fn(
  FileName = file.path(tempdir[1], "Aniso.png"),
  Report = get("Save", envir = e1)$Report,
  TmbData = TmbDataMaster)
SpatialDeltaGLMM::PlotAniso_Fn(
  FileName = file.path(tempdir[2], "Aniso.png"),
  Report = get("Save", envir = e2)$Report,
  TmbData = TmbDataMaster)
Q <- SpatialDeltaGLMM::QQ_Fn(TmbData = TmbDataMaster,
  Report = get("Save", envir = e1)$Report,
  save_dir = tempdir[1],
  FileName_PP = "Posterior_Predictive",
  FileName_Phist = "Posterior_Predictive-Histogram",
  FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")
Q <- SpatialDeltaGLMM::QQ_Fn(TmbData = TmbDataMaster,
  Report = get("Save", envir = e2)$Report,
  save_dir = tempdir[2],
  FileName_PP = "Posterior_Predictive",
  FileName_Phist = "Posterior_Predictive-Histogram",
  FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")
Enc_prob1 <- SpatialDeltaGLMM::Check_encounter_prob(
  Report=get("Save", envir = e1)$Report,
  Data_Geostat=Data_Geostat, DirName=tempdir[1])
Enc_prob2 <- SpatialDeltaGLMM::Check_encounter_prob(
  Report=get("Save", envir = e2)$Report,
  Data_Geostat=Data_Geostat, DirName=tempdir[2])

#### AIC plot
aicinfo <- aggregate(AIC ~ om_depth + om_sigmao1 + rep + om_type + em_type,
  data = keep_re, FUN = function(x) ifelse(any(is.na(x)), NA, diff(x)))
aicinfo$AIC <- unlist(aicinfo$AIC)
aicinfo$om_sigmao1 <- round(aicinfo$om_sigmao1, 2)
aicinfo$fixed <- ifelse(aicinfo$om_sigmao1 == 0.5, "0.5", "Observed")

temp <- aggregate(AIC ~ om_depth + fixed + om_type + em_type,
  data = aicinfo,
  function(x) format(sum(x<=2) / length(x), 2, nsmall = 2, trim = TRUE))
gg <- ggplot(aicinfo) +
  geom_point(aes(x = fixed, y = AIC)) +
  geom_hline(yintercept = c(-2), col = "red") +
  # geom_hline(yintercept = c(-5), col = "red") +
  # geom_hline(yintercept = c(-10), col = "red") +
  facet_grid(om_depth ~ ., scales = "free") +
  xlab(expression(sigma[omega[1]])) +
  ylab(expression(Delta~AIC)) +
  VAST_simulation_depth_theme +
  geom_text(data = temp, vjust = 1.2, hjust = -0.5,
    y = -2, aes(x = fixed, label = AIC))
ggsave(file.path(DateDir,
  "VAST_simulation_depth_em_delAIC.jpeg"),
  gg, dev = "jpeg")
rm(temp)

results <- keep_ae[grepl("Conventional", keep_ae$om_type) &
  grepl("Conventional", keep_ae$em_type), ]
results$om_sigmao1_factor <- ifelse(results$om_sigmao1 == 0.5, 0.5, "Observed")
gg <- ggplotre(results[results$em_depth == TRUE, ],
  "em_depth1_km", "em_depth2_km",
  labels = c(
    "mean depth effect for encounters",
    "mean depth effect for positive catch rates"),
  type = "points", facetx = c("om_sigmao1_factor"), facety = c("om_depth"),
  dir = NULL, scales = "fixed", print = FALSE)
gg <- gg  +
  geom_vline(aes(xintercept = om_depth1_km),
    col = "red") +
  geom_hline(aes(yintercept = om_depth2_km),
    col = "red")
gg <- addMAE(results[results$em_depth == TRUE, ],
  c("ae_depth1_km", "ae_depth2_km"),
  "om_depth + om_sigmao1_factor + em_depth", gg = gg)
ggsave(file.path(DateDir,
  "VAST_simulation_depth_em_depth1kmVSem_depth2_km.jpeg"),
  gg, dev = "jpeg")

gg <- ggplotre(results,
  "em_depth", "re_logratiodiff",
  labels = c(
    "", "error in the estimated trend of the index"),
  type = "box", facetx = c("om_sigmao1_factor"), facety = c("om_depth"),
  dir = DateDir, scales = "fixed")
ggsave(file.path(DateDir,
  "VAST_simulation_depth_em_typeVSre_logratiodiff.jpeg"),
  gg, dev = "jpeg")

gg <- ggplotre(results,
  "em_range1", "em_sigmao1",
  labels = c(
    "range the spatial and spatiotemporal fields for encounters",
    "standard deviation of the spatial field for encounters"),
  type = "points", facetx = c("om_sigmao1_factor", "em_depth"),
  facety = c("om_depth"),
  dir = NULL, scales = "fixed", print = FALSE, gradient = FALSE)
gg <- gg + geom_hline(aes(yintercept = as.numeric(om_sigmao1)), col = "red") +
  geom_vline(aes(xintercept = om_range1), col = "red")
gg <- addMAE(results, y = c("ae_range1", "ae_sigmao1"),
  "om_depth + om_sigmao1_factor + em_depth", gg = gg,
  nsmall = 2)
ggsave(file.path(DateDir,
  "VAST_simulation_depth_re_range1VSre_sigmao1.jpeg"),
  gg, dev = "jpeg")
# ggg <- ggplotre(
#   results[results$em_depth == FALSE & results$om_depth == TRUE, ],
#   "re_range1", "re_sigmao1",
#   labels = c("", ""),
#   type = "points", facetx = c("em_depth"), facety = c("om_depth"),
#   dir = NULL, scales = "fixed", print = FALSE, gradient = TRUE)
# ggg <- ggg + xlim(c(-0.9, 0.9)) + ylim(c(-0.9, 0.9)) + ylab("") + xlab("")+
#   scale_colour_gradient(guide = FALSE)
# vp <- grid::viewport(width = 0.35, height = 0.35,
#   x = 0.26, y = 0.312)
# jpeg(file.path(DateDir, "VAST_simulation_depth_re_range1VSre_sigmao1.jpeg"),
#   res = 600, units = "in", height = 8, width = 8)
#   print(gg + aes(shape = om_sigmao1))
#   theme_set(theme_bw(base_size = 8))
#   # print(ggg, vp = vp)
# dev.off()
# Can't use the gradient to tell the difference between results because
# some are unbiased but have a high gradient. Turning this off.


#### Appendix A
matched <- keep_re[keep_re$om_depth == TRUE & keep_re$om_sigmao1 == 0.5, ]
matched <- matched[ifelse(matched$om_type != matched$em_type, FALSE, TRUE), ]

# Error in log ratio of last and first years
gg <- ggplotre(matched, "re_logratiodiff", "re_sigmao1",
  labels = c(
    "error in logged ratio of the first and last year",
    "RE in standard deviation of the spatial field for the first component"),
  type = "points",
  facety = c("om_type"), facetx = c("em_depth"),
  dir = NULL, gradient = TRUE, scales = "free")
ggsave(file.path(DateDir,
  "VAST_simulation_depth_AppendixA_re_sigmao1VSre_logratiodiff.jpeg"),
  gg, dev = "jpeg")

gg <- ggplotre(matched[matched$em_depth == TRUE, ],
  "em_depth1_km", "em_depth2_km",
  labels = c(
    "mean depth effect for the first model component",
    "mean depth effect for the second model component"),
  type = "points", facetx = c(""), facety = c(""),
  dir = NULL)
gg <- gg + geom_line(data = data.frame(
    x = c(-0.25, 0.5), y = as.numeric(unique(matched$om_depth2_km)[2])),
    aes(x = x, y = y), col = 2) +
  geom_line(data = data.frame(
    x = c(1.0, 1.4), y = as.numeric(unique(matched$om_depth2_km)[1])),
    aes(x = x, y = y), col = 1) +
  geom_vline(aes(xintercept = as.numeric(om_depth1_km),
    col = as.factor(om_type))) +
  scale_colour_manual(guide = FALSE, values = 1:2)
ggsave(file.path(DateDir,
  "VAST_simulation_depth_AppendixA_em_depth1_kmVSem_depth2_km.jpeg"),
  gg, dev = "jpeg")
dev.off()


# Linear trend in log-abundance
# We next compare the trend in true abundance (fitting a log-linear model to true abundance) to the trend in estimated abundance (Fig. 2).  This again shows that the model is unbiased in each scenario.
# caption -- Error in estimates of trend for scenarios with stable, increasing, or decreasing trends in abundance, where blue line is true expected trend, and the number in the top-left is the average estimated trend for each scenario (true expected trend is 0, -0.1, and 0.1 for the 1st, 2nd, and 3rd panels, respectively)
# with(keep, plot(as.numeric(em_linear) - as.numeric(om_linear),
#   gradient))

# Check density-dependent catchability
# Finally, we calculate density-dependent catchability (whether estimated indices are hyperstable or hypersensitive to changes in the true index).  This shows a delta of approx. 1.0 for every scenario, as also shown in Thorson et al. 2015 ICESJMS.
# caption -- Test for hypersensitive or hyperstable indices (Estimate = 1.0 implies well-calibrated sensitivity)


# TODO: all the below needs to be altered.
# DF = data.frame("True"=as.vector(Results_index[,,"True","Value",]), "Est"=as.vector(Results_index[,,"Est","Value",]), "Scenario_"=dimnames(Results_index)[[1]], "RepNum"=rep(1:dim(Results_index)[2],each=dim(Results_index)[1]), "Year"=rep(1:dim(Results_index)[5],each=prod(dim(Results_index)[1:2])))
# Lm = lme4::lmer(log(Est) ~ 0 + Scenario_ + log(True):Scenario_ + (1 | RepNum:Scenario_), data=DF)
# Coef = summary(Lm)$coef
# Coef[grep("True", rownames(Coef)), ]

if (exists("cl")) stopCluster(cl)
