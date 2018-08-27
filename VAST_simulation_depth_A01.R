
#### Appendix A
appADir <- "VAST_simulation_depth_2018-03-27"
appA <- get_results(file.path(RootDir, appADir))
appA_re <- re_all(data = appA)
matched <- appA_re[appA_re$om_depth == "depth" & appA_re$om_sigmao1 == 0.5, ]
# matched <- matched[ifelse(matched$om_type != as.character(matched$em_type), FALSE, TRUE), ]
matched <- matched[matched$om_type == "Conventional delta", ]

# Error in log ratio of last and first years
gg <- ggplotre(matched, "re_logratiodiff", "re_sigmao1",
  labels = c(
    "error in logged ratio of the first and last year",
    "RE in standard deviation of the spatial field for the occupancy model"),
  type = "points",
  facety = c("em_type"), facetx = c("em_depth", "om_depth"),
  dir = NULL, gradient = FALSE, scales = "free")
ggsave(paste0(gsub(paste0("VAST_simulation_depth_", Date),
  appADir, DateDir),
  "VAST_simulation_depth_AppendixA_re_sigmao1VSre_logratiodiff.jpeg"),
  gg, dev = "jpeg")

gg <- ggplotre(matched[matched$em_depth == TRUE, ],
  "em_depth1_km", "em_depth2_km",
  labels = c(
    "depth effect for the occupancy model",
    "depth effect for the catch-rate model"),
  type = "points", facetx = c("em_type"), facety = c(""),
  dir = NULL)
gg <- gg +
  geom_vline(aes(xintercept = as.numeric(om_depth1_km)), col = "red", lty = 2) +
  geom_hline(aes(yintercept = as.numeric(om_depth2_km)), col = "red", lty = 2) +
  scale_colour_manual(guide = FALSE, values = 1:2)
ggsave(paste0(gsub(paste0("VAST_simulation_depth_", Date),
  appADir, DateDir),
  "VAST_simulation_depth_AppendixA_em_depth1_kmVSem_depth2_km.jpeg"),
  gg, dev = "jpeg")
dev.off()

