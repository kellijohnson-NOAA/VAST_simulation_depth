VAST_simulation_depth_theme <- theme_bw() + theme(
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1))

get_results <- function(dir, truncate = NULL) {
  dir.reps <- dir(dir, pattern = "rep=", full.names = TRUE,
    recursive = TRUE, include.dirs = TRUE)
  if (is.null(truncate)) truncate <- seq_along(dir.reps)

  get_OM <- function(file) {
    if (!file.exists(file)) return(NULL)
    load(file)
    index <- Sim$B_tl[,1]/1000
    names(index) <- paste0("om_index_", seq_along(index))
    linear <- lm(log(index) ~ 1 + I(seq_along(index)))$coef
    names(linear) <- NULL
    sets <- Sim$Sim_Settings[c("Range1", "Range2", "SigmaO1",
      "SigmaO2", "SigmaE1", "SigmaE2", "Depth1_km", "Depth2_km")]
    names(sets) <- paste0("om_", tolower(names(sets)))
    om <- unlist(c(index, sets))
    return(c("om_type" = gsub(".+OM=|/.+$", "", file), om,
      "om_linear" = linear[2],
      "n_t" = length(index),
      "om_depth" = ifelse(grepl("nodepth", file), FALSE, TRUE)))
  }
  get_EM <- function(file) {
    if (!file.exists(file)) return(NULL)
    load(file)
    if (is.na(EmSave$Report$jnll)) return(NULL)
    if (any(is.na(EmSave$ParHat))) return(NULL)
    em <- unlist(c(
      setNames(EmSave$Report[c("Range_raw1", "Range_raw2")],
        c("em_range1", "em_range2")),
      setNames(EmSave$ParHat[
        c("L_omega1_z", "L_omega2_z", "L_epsilon1_z", "L_epsilon2_z")],
        c("em_sigmao1", "em_sigmao2", "em_sigmae1", "em_sigmae2")),
      setNames(
        c(mean(EmSave$ParHat[['gamma1_ctp']]),
          mean(EmSave$ParHat[['gamma2_ctp']])),
        c("em_depth1_km", "em_depth2_km"))))
    em[grep("sigma", names(em))] <- abs(em[grep("sigma", names(em))])
    em_type <- gsub(".+rep=[0-9]+/EM=|/EmSave.+$", "", file)
    if (length(dir(dirname(file),
      full.names = TRUE, pattern = "Table_for_SS3")) == 0) {
      index <- rep(NA, times = EmSave$Data_EM$n_t)
      linear <- rep(NA, 2)
    } else {
      index <- read.table(dir(dirname(file),
        full.names = TRUE, pattern = "Table_for_SS3"), sep = ",", header = TRUE)
      index <- index[, "Estimate_metric_tons"]
      linear <- lm(log(index) ~ 1 + I(seq_along(index)))$coef
      names(linear) <- NULL
    }
    AIC <- dir(dirname(file),
      pattern = "estimates.txt", full.name = TRUE)
    if (length(AIC) == 0) {
      AIC <- NA
    } else {
      AIC <- readLines(dir(dirname(file),
        pattern = "estimates.txt", full.name = TRUE))
      AIC <- as.numeric(gsub("\\[1\\] ", "", AIC[grep("AIC", AIC) + 1]))
    }
    names(index) <- paste0("em_index_", seq_along(index))
    return(c("em_type" = gsub("TRUE|FALSE|/", "", em_type),
      "em_depth" = grepl("TRUE", em_type),
      em, index, "em_linear" = linear[2],
      "AIC" = AIC,
      "gradient" = ifelse(is.null(EmSave$Sdreport), NA, max(abs(EmSave$Sdreport$gradient.fixed))),
      "hessian" = ifelse(is.null(EmSave$Sdreport), NA, EmSave$Sdreport$pdHess))) #<1e-6
  }
  info <- lapply(dir.reps[truncate], function(x) {
    omom <- get_OM(file.path(x, "Sim.RData"))
    if (is.null(omom[1])) return(NULL)
    dir.ems <- dir(x, pattern = "EmSave", recursive = TRUE,
      full.names = TRUE, include.dirs = TRUE)
    emem <- mapply(get_EM, dir.ems)
    if (is.list(emem)) {
      if (TRUE) return(NULL)
      new <- as.data.frame(emem[!sapply(emem, is.null)])
      colnames(new) <- dir.ems[!sapply(emem, is.null)]
      emem <- new
      rm(new)
    }
    emem <- t(emem)
    rownames(emem) <- NULL
    if (length(emem) == 0) return(NULL)
    stuff <- data.frame(t(omom), emem, stringsAsFactors = FALSE)
    # Make long dataframe of index
    indexdatatouse <- stuff[, grepl("index|rep", colnames(stuff))]
    indexdata <- reshape(indexdatatouse, direction = "long",
      varying = colnames(indexdatatouse), sep = "_index_")
    indexdata <- as.data.frame(apply(indexdata, 2, as.numeric))
    indexcoef <- mapply(coef, tapply(1:NROW(indexdata), indexdata$id,
      function(x) {
      lm(I(log(em)) ~ time + I(log(om)), data = indexdata[x, ])
    }))[3, ]
    stuff$delta <- indexcoef
    stuff$rep <- gsub(".+rep=", "", x)
    return(stuff)
  })

  dims <- sapply(lapply(info, dim), "[[", 2)
  bad <- which(sapply(dims, is.null))
  bad <- c(bad, which(sapply(dims,
    function(x) any(x %in%
      names(table(unlist(dims))[-length(unique(unlist(dims)))])))))
  info <- info[-bad]

  info <- do.call("rbind", info)
  info$om_type <- gsub("_[0-9]+", "", info$om_type)
  # colschange <- which(apply(info, 2,
  #   function(x) any(grepl("e-\\d+|e+\\d+", x))))
  # numbers <- apply(info, 2,
  #   function(x) stringr::str_extract(x, "\\-*\\d+\\.*\\d*"))
  # numbers[, names(colschange)] <- apply(info[, names(colschange)],
  #   2, as.numeric)

  # info <- cbind(
  #   info[, apply(numbers, 2, function(x) all(is.na(x)))],
  #   apply(numbers[, !apply(numbers, 2, function(x) all(is.na(x)))], 2,
  #     as.numeric))

  info$om_logratio <-
    log(as.numeric(info[, paste0("om_index_", info[1, "n_t"])])) -
    log(as.numeric(info[, paste0("om_index_1")]))
  info$em_logratio <-
    log(as.numeric(info[, paste0("em_index_", info[1, "n_t"])])) -
    log(as.numeric(info[, paste0("em_index_1")]))

  return(info)
}

re_all <- function(data, bind = TRUE, type = c("re", "ae")) {
  type <- match.arg(type)
  re <- function(text = "index_1", data) {
    omd <- as.numeric(data[, paste0("om_", text)])
    emd <- as.numeric(data[, paste0("em_", text)])

    return((emd - omd) / omd)
  }
  if (type == "ae") {
    re <- function(text = "index_1", data) {
      omd <- as.numeric(data[, paste0("om_", text)])
      emd <- as.numeric(data[, paste0("em_", text)])
      return(abs(emd - omd))
    }
  }
  names <- gsub("om_", "", grep("om_", colnames(data), value = TRUE))
  names <- names[!(names %in% c("type"))]
  names <- names[!(names %in% c("logratio"))]
  names <- names[!grepl("depth$|linear|logratio", names)]

  res <- mapply(re, names, MoreArgs = list(data = data))
  colnames(res) <- paste0(type, "_", colnames(res))
  res <- data.frame(res,
    "re_logratiodiff" = data$em_logratio - data$om_logratio)

  indexnames <- grep("index", names, value = TRUE)
  om <- data[, grepl("om_index_", colnames(data))]
  em <- data[, grepl("em_index_", colnames(data))]
  om <- apply(om, 2, as.numeric)
  em <- apply(em, 2, as.numeric)
  om <- t(apply(om, 1, function(x) log(x) - mean(log(x))))
  em <- t(apply(em, 1, function(x) log(x) - mean(log(x))))
  logre <- em - om
  colnames(logre) <- gsub("em", "lncentdiff", colnames(logre))

  res <- cbind(res, logre)
  res$rmse_index <- apply(logre, 1, function(x) sqrt(mean(x^2)))

  if (bind) {
    return(cbind(data, res))
  } else {return(res)}

}

ggplotre <- function(data, x, y, print = TRUE, gradient = FALSE,
  facetx = c("emname", "em_depth"), facety = "omname",
  labels = NULL, type = "box", scales = c("fixed", "free"),
  dir = NULL, lim = list(x = NULL, y = NULL)) {

  scales <- match.arg(scales)
  facetx <- facetx[!facetx %in% x]

  data$em_depth <- paste("depth in EM =", data$em_depth)

  data$em_type <- gsub("process ", "", data$em_type)
  data$om_type <- gsub("process ", "", data$om_type)
  data$em_type <- gsub("Conventional delta", "Delta-model", data$em_type)
  data$om_type <- gsub("Conventional delta", "Delta-model", data$om_type)
  data$omname <- paste("OM =", data$om_type)
  data$emname <- paste("EM =", data$em_type)

  if (x != "em_depth") data[, x] <- as.numeric(data[, x])
  data[, y] <- as.numeric(data[, y])

  gg <- ggplot(data)
  if (gradient) gg <- ggplot(data,
    aes(col = as.numeric(gradient)))

  if (type == "points"){
    gg <- gg +
      geom_point(aes_string(x = x, y = y),
      alpha = 0.5, size = 2.5) +
      scale_shape_manual(name = "", values = c(15, 19), guide = FALSE) +
      scale_color_gradient(name = "gradient")
    if (grepl("^re", y)) {
      gg <- gg + geom_vline(xintercept = 0, lty = 2, col = "red")
    }
  }
  if (type == "box") {
    gg <- gg +
      geom_boxplot(aes_string(x = x, y = y))
  }
  if (facetx[1] == "" & facety[1] == "") {
    facet <- "."
  } else {
    gg <- gg + facet_grid(as.formula(paste(paste(facety, collapse = "+"), "~",
      paste(facetx, collapse = "+"))), scales = scales)
  }

  gg <- gg + VAST_simulation_depth_theme
  if (grepl("^re", y)) {
    gg <- gg + geom_hline(yintercept = 0, lty = 2, col = "red")
  }


  if(!is.null(labels[1])) gg <- gg +
    xlab(labels[1]) + ylab(labels[2])
  if (!is.null(lim$x)) gg <- gg + xlim(lim$x)
  if (!is.null(lim$y)) gg <- gg + ylim(lim$y)

  if (print) print(gg)
  if (!is.null(dir)) {
    jpeg(
      filename = file.path(dir,
        paste0("VAST_simulation_depth_", x, "VS", y, ".jpeg")),
      res = 600, units = "in", width = 8, height = 8)
    print(gg)
    dev.off()
  }
  invisible(gg)
}

getrdata <- function(file){
  ne <- new.env()
  load(file, env = ne)

  out <- data.frame(
    "par" = names(
      get("parameter_estimates", env = ne)$SD[c("par.fixed")][[1]]),
    "val" = get("parameter_estimates", env = ne)$SD[c("par.fixed")],
    "se" = sqrt(diag(
      get("parameter_estimates", env = ne)$SD[c("cov.fixed")][[1]])))
  colnames(out)[2] <- "val"
  report <- data.frame(
    "par" = names(get("parameter_estimates", env = ne)$SD$value),
    "val" = get("parameter_estimates", env = ne)$SD$value,
    "se" = sqrt(diag(get("parameter_estimates", env = ne)$SD$cov)))
out$par <- make.unique(as.character(out$par))
report$par <- make.unique(as.character(report$par))

  # ne2 <- new.env()
  # load(file.path(dirname(file), "Save.RData"), env = ne2)
  # depth <- data.frame(
  #   "par" = c("gamma1_ctp", "gamma2_ctp"),
  #   "val" = sapply(get("Save",
  #     env = ne2)$ParHat[c("gamma1_ctp", "gamma2_ctp")], mean),
  #   "se" = NA)
  # report <- rbind(report, depth)

  out <- rbind(out, report)
  out$om_name <- basename(dirname(file))

  out[grep("L_", out$par), "val"] <- abs(out[grep("L_", out$par), "val"])
  return(out)
}

#' Calculate confidence intervals for two sets of data (x and y) and
#' determine if x is contained in y and vice-versa.
#'
#' @param names A vector of names for the rows in the resulting data.
#'
int95 <- function(x, y, x.se, y.se, names = NULL) {
  x.se <- ifelse(is.na(x.se), 0, x.se)
  y.se <- ifelse(is.na(y.se), 0, y.se)
  y <- ifelse(is.na(y), 0, y)
  low <- cbind(x, y) - 1.96 * cbind(x.se, y.se)
  hig <- cbind(x, y) + 1.96 * cbind(x.se, y.se)
  x.int <- cbind(y, low[, 1], hig[, 1])
  y.int <- cbind(x, low[, 2], hig[, 2])

  giveme <- data.frame(
    x, x.low = low[, 1], x.upp = hig[, 1],
    y, y.low = low[, 2], y.upp = hig[, 2],
    yinx = apply(x.int, 1, function(xx) findInterval(xx[1], xx[2:3][order(xx[2:3])])),
    xiny = apply(y.int, 1, function(xx) findInterval(xx[1], xx[2:3][order(xx[2:3])])))
  row.names(giveme) <- names
  return(giveme)
}

#' Calculate the median absolute error over categories and add to a plot
#'
#' @param data The data used in the \code{\link{ggplot2}} plot
#' @param y A character vector of length 2 that specifies the
#' columns to summarize
#' @param x A character value providing the right side of the
#' formula used for aggregation (e.g., \code{"depth + year"})
#' @param gg The \code{\link{ggplot2}} you want to add the MAE to
#' @param nsmall An integer value providing the number of digits you
#' want after the decimal for the MAE
#'
addMAE <- function(data, y, x, gg = NULL, nsmall = 3) {
  if (grepl("em_depth", x)) {
    data$em_depth <- paste("depth in EM =", data$em_depth)
  }
  agg1 <- aggregate(as.formula(paste(y[1], "~", x)),
    data = data,
    function(xx) format(median(xx, na.rm = TRUE), nsmall = nsmall,
      digits = 1))
  agg2 <- aggregate(as.formula(paste(y[2], "~", x)),
    data = data,
    function(xx) format(median(xx, na.rm = TRUE), nsmall = nsmall,
      digits = 1))
  aggs <- merge(agg1, agg2)
  if (is.null(gg)) return(aggs)

  gg <- gg + geom_text(data = aggs, x = Inf, y = -Inf,
  vjust = -1.5, hjust = 1.5, aes_string(label = y[1])) +
  geom_text(data = aggs, x = -Inf, y = Inf,
    vjust = 2.5, hjust = -0.5, aes_string(label = y[2]))
  return(gg)
}

scaleFUN <- function(x) sprintf("%.0f", x)

keep_num <- function(data, number = 100, gradientval = 0.0001) {
  data <- subset(data, gradient < gradientval)
  data <- data[with(data, order(om_type, om_depth, em_depth, sigmao1_factor)), ]
  data$id <- as.factor(apply(data[, cols], 1, paste0, collapse = "__"))
  data$id <- sequence(tabulate(data$id))
  data <- data[data$id <= number, ]
  return(data)
}

#' Run the Estimation model
#'
#' @param Date A character value giving the date
#' @param RootDir A file path to the base directory
#' @param OM The folder for the OM inside of the \code{RootDir}
#' @param EM The folder for the EM inside of the \code{rep} directory
#' @param rep An integer giving the replicate number
#' @param depth A logical value if depth is included in the \code{EM}
#'
rerun <- function(Date = "2018-05-12",
  RootDir, OM, EM, rep, depth) {

  DateDir <- file.path(RootDir, paste0("VAST_simulation_depth_", Date), "/")
  DownloadDir <- paste0(DateDir, "downloads/")
  OmDir <- paste0(DateDir, OM, "/")
  RepDir <- paste0(OmDir, "rep=", rep, "/")
  load(file = paste0(DateDir, "Record.RData"))
  attach(Record)
  if (grep("OM=", OmDir)) {
    EmDir <- OmDir
  } else { EmDir = paste0(RepDir,"EM=", EM, "/", depth, "/") }

  load(dir(RepDir, pattern = "Sim.R", full.names = TRUE))
  # load(dir(EmDir, pattern = "EmSave", full.names = TRUE))
  load(dir(OmDir, pattern = "Save.RData", full.names = TRUE))
  for(ii in seq_along(Save)) {
    assign(names(Save)[ii], Save[[ii]])
  }
  for(ii in seq_along(Record)) {
    assign(names(Record)[ii], Record[[ii]])
  }
  # for(ii in seq_along(EmSave)) {
  #   assign(names(EmSave)[ii], EmSave[[ii]])
  # }
  for(ii in seq_along(Sim)) {
    assign(names(Sim)[ii], Sim[[ii]])
  }
  region <- ifelse(grepl("WCGBTS", Species),
    "California_current",
    "eastern_bering_sea")
  Extrapolation_List <- SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(
    Region=region, strata.limits=strata.limits)
 Spatial_List <- Spatial_raw
   Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
   Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
  if (file.exists(paste0(EmDir, "parameter_estimates.RData"))) {
    parenv <- new.env()
    load(paste0(EmDir, "parameter_estimates.RData"), envir = parenv)
    parvals <- get("parameter_estimates", envir = parenv)$diagnostics
    badaniso <- as.character(parvals$Param[
    which.max(parvals$final_gradient)]) == "ln_H_input"
    badaniso <- ifelse(badaniso, 1, 0)
  } else {badaniso = 1}

  depthdata <- NULL
  if (depth) depthdata <- Data_raw$X_xtp
  if (!grepl("OM", EM)) {
    OverdispersionConfig <- c("eta1"=0, "eta2"=0)
    vyeffect <- rep(0, length(as.numeric(Sim$Data_Geostat[,'Vessel'])-1))
  } else {
    vyeffect <- as.numeric(Sim$Data_Geostat[,'Vessel'])-1
  }

      TmbDataEM <- VAST::Data_Fn(
        "Version"=Version,
        "X_xtp"=depthdata,
        "FieldConfig"=FieldConfig,
        "OverdispersionConfig"=OverdispersionConfig,
        "RhoConfig"=RhoConfig,
        "ObsModel"=ObsModel_Set[[which(sapply(names(ObsModel_Set), grep, EM) ==1)]],
        "c_i"=rep(0,nrow(Sim$Data_Geostat)),
        "b_i"=Sim$Data_Geostat[,'Catch_KG'],
        "a_i"=Sim$Data_Geostat[,'AreaSwept_km2'],
        "v_i"=vyeffect,
        "s_i"=Spatial_List$knot_i-1,
        "t_i"=Sim$Data_Geostat[,'Year'],
        "a_xl"=Spatial_List$a_xl,
        "MeshList"=Spatial_List$MeshList,
        "GridList"=Spatial_List$GridList,
        "Method"=Spatial_List$Method, "Options"=Options,
        "Aniso" = badaniso)
      TmbListEM  = Build_TMB_Fn(
        "TmbData"=TmbDataEM, "RunDir"=RootDir,
        "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x)
      # Run model
      ObjEM = TmbListEM[["Obj"]]
      OptEM = TMBhelper::Optimize(obj=ObjEM,
        lower=TmbListEM[["Lower"]], upper=TmbListEM[["Upper"]],
        getsd=TRUE, savedir=EmDir, newtonsteps=1, bias.correct=TRUE,
        bias.correct.control = list(sd = FALSE, split = NULL, nsplit = 1,
          vars_to_correct = "Index_cyl"))
      ReportEM = ObjEM$report()
      if (is.null(class(OptEM$par))) {
        ParHat <- NULL
      } else {ParHat = ObjEM$env$parList(OptEM$par)}
      if (!is.null(OptEM[["SD"]])){
        IndexEM = SpatialDeltaGLMM::PlotIndex_Fn(DirName=EmDir,
          TmbData=TmbDataEM, Sdreport=OptEM[["SD"]],
          Year_Set=Year_Set, strata_names=strata.limits[,1], use_biascorr=TRUE)
      } else {
        IndexEM$Index_ctl <- array(NA, dim = dim(Index$Index_ctl),
          dimnames = list(NULL, NULL, NULL, c("Estimate", "Std. Error")))
        ParHat <- list('L_omega1_z' = NA,'L_epsilon1_z' = NA,
          'L_omega2_z' = NA,'L_epsilon2_z' = NA,
          'gamma1_ctp' = NA, 'gamma2_ctp' = NA)
      }
    AIC <- OptEM$AIC
      EmSave = list("Report"=ReportEM, "ParHat"=ParHat, "Data_EM" = TmbDataEM,
        "Sdreport" = OptEM[["SD"]], "AIC" = AIC)
      save(EmSave, AIC, file=paste0(EmDir,"EmSave.RData"))
      invisible(EmSave)
}
