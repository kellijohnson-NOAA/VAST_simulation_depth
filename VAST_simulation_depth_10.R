###############################################################################
###############################################################################
#### Results
###############################################################################
###############################################################################
keep$sigmao1_factor <- ifelse(keep$om_sigmao1 == 0.5, 0.5, "Observed")
keep_re <- VASTWestCoast:::re_all(data = keep)
keep_ae <- VASTWestCoast:::re_all(data = keep, type = "ae")

# depth2 file
filesquared <- file.path(dir(RootDir, pattern = "^02", full.names = TRUE))
envsquared <- new.env()
load(dir(filesquared, pattern = "setup", full.names = TRUE), envir = envsquared)
if (get("settings", envir = envsquared)$depth != "squared") stop("The wrong file",
  " is being called with filesquared.")

###############################################################################
###############################################################################
#### Tables
#### Calculate the depth range
###############################################################################
###############################################################################
locallimits <- data.frame('STRATA'="All_areas")
localregion <- ifelse(grepl("WCGBTS", Sim_Settings$Species),
  "California_current", NA)
if (is.na(localregion)) stop("Need to set the local region because",
  " it is not the California current.")
localinfo <- get("info", envir = envsquared)
localextrapolationlist <- localinfo$Extrapolation_List
localspatiallist <- localinfo$Spatial_List
initialDataDownload$knot_i <- localspatiallist$knot_i
utils::data(california_current_grid, package = "SpatialDeltaGLMM")
depthsmean <- tapply(california_current_grid[, "Depth_km"],
  INDEX = localspatiallist$PolygonList$NN_Extrap$nn.idx, FUN = mean)
# Get region-specific settings for plots
MapDetails_List <- FishStatsUtils::make_map_info(
  Region = localregion,
  NN_Extrap = localspatiallist$PolygonList$NN_Extrap,
  Extrapolation_List = localextrapolationlist)

# Spit out information to console about raw empirical data
range(depthsmean)*1000
range(depthsmean[unique(initialDataDownload[
  initialDataDownload$Wt>0,]$knot_i)])*1000
mean(aggregate(Wt ~ Year, data = initialDataDownload,
  function(x) sum(x > 0)/length(x))[, 2])
mean(with(initialDataDownload,
  tapply(Wt, Year, function(x) sum(x>0)/length(x))))
range(with(initialDataDownload, tapply(Wt, Year, length)))

# Fits to empirical data
rawests <- Reduce(function(x,y) merge(x, y, all = TRUE),
  lapply(dir(dir(RootDir, full.names = TRUE, pattern = "^[[:digit:]]+_VAST_"),
    recursive = FALSE, "parameter_estimates.RData", full.names = TRUE),
  VASTWestCoast:::getrdata))
rawests <- reshape(rawests,
  direction = "wide", idvar = "par", timevar = "om_name")
colnames(rawests) <- gsub(" delta|_observed|entional", "", colnames(rawests))
rawests[grep("L", rawests$par), ]
rawests[grep("gamma", rawests$par), ]

# AIC of empirical models
# todo: this needs to be made more robust
newenv <- new.env()
diraic <- dir(dir(RootDir, full.names = TRUE, pattern = "[[:digit:]]+_VAST_"),
    recursive = FALSE, "parameter_estimates.RData", full.names = TRUE)
diraic <- sapply(diraic, function(x) {
  load(x, envir = newenv)
  get("parameter_estimates", envir = newenv)$AIC
})
diraic - min(diraic[-4])
rm(newenv)

###############################################################################
###############################################################################
#### Simulation table
###############################################################################
###############################################################################
write.csv(aggregate(coverage ~ nknots + EM +
  om_depth1_km + om_depth1_km2 + om_sigmao1 + om_range1,
  data = keep, mean),
  row.names = FALSE,
  file = file.path(RootDir, "VAST_simulation_depth_coverage.csv"))

###############################################################################
###############################################################################
#### Plots
###############################################################################
###############################################################################

#### Map of study area
# download the data from
# https://www.arcgis.com/home/item.html?id=f7f805eb65eb4ab787a0a3e1116ca7e5
url <- 'https://www.arcgis.com/sharing/rest/content/items/f7f805eb65eb4ab787a0a3e1116ca7e5/data'
fil <- "states_21basic.zip"
if (!file.exists(file.path(DownloadDir, fil))) {
  download.file(url, file.path(DownloadDir, fil))
}
  unzip(file.path(DownloadDir, fil), exdir = DownloadDir)

states <- rgdal::readOGR(file.path(DownloadDir, "states.shp"),
  layer = "states", stringsAsFactors = FALSE)
states <- sp::spTransform(states, sp::CRS("+proj=longlat +datum=WGS84"))
states <- ggplot2::fortify(states)

box <- make_bbox(lat = Lat, lon = Long, data = initialDataDownload)
big <- get_map(location = box,
  color = "bw",
  maptype = "terrain",
  zoom = 5)

gg <- ggplot() +
geom_point(data = initialDataDownload[initialDataDownload$Wt > 0, ],
  aes(x = Long, y = Lat, cex = Wt), alpha = 0.2) +
coord_map(xlim = c(-126, -118), ylim = c(33, 48.5)) +
geom_map(data=states, map = states, aes(x=long,y=lat,map_id=id),
  col = "black", alpha = 0.2, size = 0.2) +
VAST_simulation_depth_theme +
xlab("longitude") + ylab("latitude") + labs(cex="catch (kg)") +
theme(legend.position = c(0.2, 0.15),
  legend.key = element_rect(fill=alpha('white', 0.1)),
  legend.background = element_rect(fill=alpha('blue', 0.01)))
# gg <- gg + ggsn::north2(gg, .65, .91)
ggsave(file.path(RootDir,
  "VAST_simulation_depth_map_rawdata.jpeg"),
  gg, dev = "jpeg", height = 8, width = 3.8)
dev.off()

#### Diagnostic plots

#Plot Anisotropy
for (ii in dir(RootDir, pattern = "^[[:digit:]]{2}", full.names = TRUE)) {
  e1 <- new.env(parent = baseenv())
  load(dir(ii, pattern = "^Save.RData", full.names = TRUE), envir = e1)
  SpatialDeltaGLMM::PlotAniso_Fn(
    FileName = file.path(ii, "Aniso.png"),
    Report = get("Report", envir = e1),
    TmbData = get("TmbData", envir = e1))
  Q <- SpatialDeltaGLMM::QQ_Fn(
    TmbData = get("TmbData", envir = e1),
    Report = get("Report", envir = e1),
    save_dir = ii,
    FileName_PP = "Posterior_Predictive",
    FileName_Phist = "Posterior_Predictive-Histogram",
    FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")
  Enc_prob1 <- SpatialDeltaGLMM::Check_encounter_prob(
    Report = get("Report", envir = e1),
    Data_Geostat = data.frame("Catch_KG" = get("TmbData", envir = e1)$b_i),
    DirName = ii)
}

#### AIC plot
names <- c("om_sigmao1", "om_depth1_km", "om_depth1_km2",
  "rep", "nknots", "om_range1", "om_depthtype")
aicinfo <- suppressWarnings(reshape(keep[,c(names, "AIC", "EM")],
  direction = "wide", v.names = "AIC",
  idvar = names, timevar = "EM"))
allContrasts <- outer(
  paste0("`", colnames(aicinfo)[-c(seq_along(names))], "`"),
  paste0("`", colnames(aicinfo)[-c(seq_along(names))], "`"),
  paste, sep = " - ") %>%
as.character %>% setNames(., .) %>% as.list()
forPlotting <-
  aicinfo %>%
  mutate_(.dots = allContrasts) %>%
  select_(.dots = paste0("-`", colnames(aicinfo)[-c(seq_along(names))], "`")) %>%
  gather(Comparison, Difference, -c(om_sigmao1, om_depth1_km, om_depth1_km2,
  rep, nknots, om_range1, om_depthtype)) %>%
  separate(Comparison, c("EMa", "EMb"), " - ") %>%
  filter(EMa != EMb) %>%
  mutate_each(funs(gsub("`", "", .)), EMa, EMb)
forPlotting$EMa <- gsub("AIC.01", "none", forPlotting$EMa)
forPlotting$EMa <- gsub("AIC.02", "linear", forPlotting$EMa)
forPlotting$EMa <- gsub("AIC.03", "quadratic", forPlotting$EMa)
forPlotting$EMb <- gsub("AIC.01", "none", forPlotting$EMb)
forPlotting$EMb <- gsub("AIC.02", "linear", forPlotting$EMb)
forPlotting$EMb <- gsub("AIC.03", "quadratic", forPlotting$EMb)
gg <-
  ggplot(subset(forPlotting, om_sigmao1 != 0.5 & om_range1 != 1000 &
    nknots == 250),
    aes(x = relevel(as.factor(EMa), "none"),
      y = Difference, col = relevel(as.factor(EMb), "none"))) +
  # geom_point() +
  geom_boxplot() +
  facet_grid(. ~ om_depthtype,
    scales = "free_y") +
  xlab("EM") +
  VAST_simulation_depth_theme +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(colour = guide_legend(title = "EM"))
ggsave(file.path(RootDir,
  "VAST_simulation_depth_em_AIC.jpeg"),
  gg, dev = "jpeg")
dev.off()

###############################################################################
#### Depth effects
for (ii in unique(keep$nknots)) {
for (iiii in c("0.5", "conditioning")) {
use <- keep_ae[
  grepl("Conventional", keep_ae$om_type) &
  grepl("Conventional", keep_ae$em_type) &
  keep_ae$nknots == ii &
  keep_ae$om_range1 != 1000, ]
if (iiii == "0.5") {
  use <- subset(use, om_sigmao1 == 0.5)
} else {
  use <- subset(use, om_sigmao1 != 0.5)
}
if (dim(use)[1] == 0) next
aa <- use[use$EM != "01", ]
gg <- VASTWestCoast:::ggplotre(aa,
  "em_depth1_km", "em_depth2_km",
  labels = c(
    "linear depth effect for encounters",
    "linear depth effect for positive catch rates"),
  type = "points",
  facetx = c("ifelse(EM == '02','linear EM', 'quadratic EM')"),
    facety = c("om_depthtype"),
  dir = NULL, scales = "fixed", print = FALSE)
gg <- gg  +
  geom_vline(aes(xintercept = om_depth1_km),
    col = "red") +
  geom_hline(aes(yintercept = om_depth2_km),
    col = "red")
gg <- VASTWestCoast:::addMAE(aa,
  c("ae_depth1_km", "ae_depth2_km"),
  "om_depthtype + sigmao1_factor + em_depth + EM", gg = gg)
ggsave(file.path(RootDir, paste0("VAST_simulation_depth_depth_",
  ii, "knots_", iiii, ".jpeg")),
  gg, dev = "jpeg")

aa <- use[use$EM == "03", ]
gg <- VASTWestCoast:::ggplotre(aa,
  "em_depth1_km2", "em_depth2_km2",
  labels = c(
    "depth-squared effect for encounters",
    "depth-squared effect for positive catch rates"),
  type = "points",
  facetx = c("."), facety = c("om_depthtype"),
  dir = NULL, scales = "fixed", print = FALSE)
gg <- gg  +
  geom_vline(aes(xintercept = om_depth1_km2),
    col = "red") +
  geom_hline(aes(yintercept = om_depth2_km2),
    col = "red")
gg <- VASTWestCoast:::addMAE(aa,
  c("ae_depth1_km2", "ae_depth2_km2"),
  "om_depthtype + sigmao1_factor + em_depth + EM", gg = gg)
ggsave(file.path(RootDir, paste0(
  "VAST_simulation_depth_depth2_quadratic", "_", ii, "knots_",
  iiii, ".jpeg")),
  gg, dev = "jpeg")
}}
dev.off()

###############################################################################
use <- keep_ae[
  grepl("Conventional", keep_ae$om_type) &
  grepl("Conventional", keep_ae$em_type) &
  keep_ae$EM %in% c("01", "03") &
  keep_ae$nknots == 250 &
  keep_ae$om_range1 != 1000 &
  keep_ae$om_sigmao1 != 0.5, ]

gg <- VASTWestCoast:::ggplotre(use,
  "em_depth", "re_logratiodiff",
  labels = c(
    "", "error in the estimated trend of the index"),
  type = "box",
  # facetx = c("sigmao1_factor"),
  facety = c("om_depthtype"),
  dir = RootDir, scales = "fixed")

gg <- VASTWestCoast:::ggplotre(use,
  "em_range1", "em_sigmao1",
  labels = c(
    "range of the spatial and spatiotemporal fields for encounters",
    "standard deviation of the spatial field for encounters"),
  type = "points",
  facetx = c("em_depth"),
  facety = c("om_depthtype"),
  dir = NULL, scales = "fixed", print = FALSE, gradient = FALSE)
gg <- gg + geom_hline(aes(yintercept = as.numeric(om_sigmao1)),
  col = "red") +
  geom_vline(aes(xintercept = om_range1), col = "red")
gg <- VASTWestCoast:::addMAE(use,
  y = c("ae_range1", "ae_sigmao1"),
  "om_depthtype + sigmao1_factor + em_depth", gg = gg,
  nsmall = 2)
ggsave(file.path(RootDir,
  "VAST_simulation_depth_re_range1VSre_sigmao1.jpeg"),
  gg, dev = "jpeg")

# Plot replicate time series.
test <- reshape(use[, c("om_depthtype", "em_depth", "rep", "EM",
    paste0("em_index_", 1:13),
    paste0("om_index_", 1:13))],
  varying = list(paste0("em_index_", 1:13),
    paste0("om_index_", 1:13)),
  v.names = c("em", "om"),
  times = 1:13, direction = "long")
gg <- ggplot(data = test[as.numeric(test$rep) %in% c(1,11), ],
  aes(col = rep, group = id)) +
  geom_line(aes(time, em, lty = em_depth), lwd = 0.5) +
  geom_line(aes(time, om), lwd = 1.5) +
  facet_grid(om_depthtype ~ ., scales = "free_y") +
  VAST_simulation_depth_theme +
  ylab("index") +
  xlab("year") +
  guides(color = FALSE)
ggsave(file.path(RootDir,
  "VAST_simulation_depth_rep1and10.jpeg"),
  gg, dev = "jpeg")

# Check density-dependent catchability
# Finally, we calculate density-dependent catchability (whether estimated indices are hyperstable or hypersensitive to changes in the true index).  This shows a delta of approx. 1.0 for every scenario, as also shown in Thorson et al. 2015 ICESJMS.
# caption -- Test for hypersensitive or hyperstable indices (Estimate = 1.0 implies well-calibrated sensitivity)

  #DF = data.frame(
  #"True"=as.vector(Results_index[,,"True","Value",]),
  #"Est"=as.vector(Results_index[,,"Est","Value",]),
  #"Scenario_"=dimnames(Results_index)[[1]],
  #"RepNum"=rep(1:dim(Results_index)[2],each=dim(Results_index)[1]),
  #"Year"=rep(1:dim(Results_index)[5],each=prod(dim(Results_index)[1:2])))
# Lm = lme4::lmer(log(Est) ~ 0 + Scenario_ + log(True):Scenario_ + (1 | RepNum:Scenario_), data=DF)
  test$scenario_ <- with(test, paste(om_depthtype, EM))
  lmres <- lme4::lmer(log(em) ~ 0 + scenario_ + log(om):scenario_ +
    (1 | rep:scenario_), data = test)
sink(file.path(RootDir, "VAST_simulation_depth_hypersable.txt"))
summary(lmres)
  Coef = summary(lmres)$coef
Coef[grep(":", rownames(Coef)), ]
sink()

#### Empirical index plot
fits <- lapply(sapply(dir(RootDir,
  pattern = "0[1-3]_VAST",
  full.names = TRUE, recursive = FALSE),
  dir, full.names = TRUE, pattern = "Table"),
  read.csv, header = TRUE, sep = ",")
fits[[3]]$depth <- "no depth"
fits[[1]]$depth <- "linear depth"
fits[[2]]$depth <- "quadratic depth"
fits <- do.call("rbind", fits)
fits <- with(fits,
  data.frame(fits,
    fit = log(fits$Estimate_metric_tons),
    lwr = log(fits$Estimate_metric_tons) - 1.96*SD_log,
    upr = log(fits$Estimate_metric_tons) + 1.96*SD_log))
rownames(fits) <- basename(rownames(fits))

gg <- ggplot(fits, aes(Year, fit, color = as.factor(depth),
  fill = as.factor(depth))) +
    geom_ribbon(data = fits,
      aes(ymin = lwr, ymax = upr), alpha = 0.6,
      show.legend = FALSE) +
    geom_line(data = fits, lwd = 2)+
    theme_bw() +
    scale_x_continuous(limits = range(initialDataDownload$Year),
      expand = c(0, 0.2),
      labels = VASTWestCoast:::scaleFUN) +
    scale_colour_brewer(palette="Spectral", name = "") +
    scale_fill_brewer(palette="Spectral", name = "") +
    # guides(fill = FALSE) +
    theme(legend.key = element_rect(colour = NA, fill = NA),
        legend.justification = c(1, 0),
        legend.position = c(0.85, 0.8)) +
    ylab("ln abundance (mt)")
ggsave(file.path(RootDir,
  "VAST_simulation_depth_rawTimeSeries.jpeg"),
  gg, dev = "jpeg")
dev.off()

###############################################################################
###############################################################################
#### Median absolute error calcs
###############################################################################
###############################################################################
# Now make running MAE calculations. Takes a while to run these calcs!!!
# todo: redo code for new use
my.median <- function(x) {
  sapply(1:length(x), function(i) {median(x[1:i], na.rm = TRUE)})
}
test <- plyr::ddply(subset(keep_ae,
  sigmao1_factor == "Observed" &
  om_range1 != 1000 &
  EM %in% c("01", "03") & sim == "02"),
  .(om_depthtype, em_depth),
  .fun = summarize,
  replicate2 = 1:length(rep),
  "MAE range in observation model" = my.median(abs(ae_range1)) -
    median(abs(ae_range1), na.rm = TRUE),
  "MAE depth for catch-rate model" = my.median(abs(ae_depth2_km)) -
    median(abs(ae_depth2_km), na.rm = TRUE))
scalars.mares <- reshape::melt(test,
  c("om_depthtype", "em_depth", "replicate2"))
gg <- ggplot(data = scalars.mares,
  aes(x = replicate2, y = value, col = em_depth)) +
    ylab("Centered MAE") + xlab("replicate")
gg <- gg + geom_line(lwd = 0.5, alpha = 1.0) +
  facet_grid(variable ~ om_depthtype,
  scales = "free") + VAST_simulation_depth_theme +
  guides(colour = guide_legend(title = "Depth in EM"))
ggsave(file.path(RootDir,
  "VAST_simulation_depth_runningmean.jpeg"),
  gg, dev = "jpeg", height = 8, width = 8)
dev.off()

pdf(file.path(RootDir, "VAST_simulation_depth_plots.pdf"))
for (ii_name in grep("re_", colnames(keep_re), value = TRUE)) {
  for (iii_name in grep("re_", colnames(keep_re), value = TRUE)) {
    if (ii_name == iii_name) next
    if (grepl("index", ii_name) | grepl("index", iii_name)) next
    ggplotre(keep_re, ii_name, iii_name, type = "points",
      facety = c("om_type", "om_depth", "om_sigmao1"), gradient = TRUE)
  }
}
dev.off()
