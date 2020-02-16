library(dplyr)
library(ggplot2)
library(viridis)

###plots of output

spc <- list.files("H:/SHEFS_joint_project/sdm_pipeline_output/ensemble/weighted/")

spc <- gsub("_ensemble.tif", "", spc)

sp_name <- "Zizula_hylax"

dir.create("H:/SHEFS_joint_project/sdm_pipeline_output/plots")

gg_ras <- function(ras){
  
  ras.p <- rasterToPoints(ras)
  ras.g <- data.frame(ras.p) 
#  ras.g$model <- names(ras)
  colnames(ras.g) <- c("Longitude", "Latitude", "Value")
  return(ras.g)
  
  }


ggplotter <- function(sp_name){
  
  xy <- read.csv(paste0("H:/SHEFS_joint_project/sdm_pipeline_output/points/", sp_name, ".csv"),  colClasses=c("NULL", NA, NA))
  # bc <- raster(paste0("H:/SHEFS_joint_project/sdm_pipeline_output/predictions/", sp_name, "_bioclim.tif"))
  # gm <- raster(paste0("H:/SHEFS_joint_project/sdm_pipeline_output/predictions/", sp_name, "_glm.tif"))
  # rf <- raster(paste0("H:/SHEFS_joint_project/sdm_pipeline_output/predictions/", sp_name, "_rf.tif"))
  ens <- raster(paste0("H:/SHEFS_joint_project/sdm_pipeline_output/ensemble/weighted/", sp_name, "_ensemble.tif"))
  
  gg_ens <- gg_ras(ens)

  ggplot(gg_ens, aes(y = Latitude, x = Longitude)) +
    geom_tile(aes(fill = Value)) +
    scale_fill_viridis_c(limits = c(0,1)) +
    geom_point(xy, mapping = aes(x = x, y = y), color = "white", shape = 4) +
    coord_quickmap() +
    theme_bw() +
    ggtitle(sp_name)
  
  ggsave(paste0("H:/SHEFS_joint_project/sdm_pipeline_output/plots/",sp_name, ".png"), width = 8, height = 8)
  print(sp_name)
  
  }


lapply(spc[60:length(spc)], ggplotter)

