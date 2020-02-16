.libPaths(c('/lustre/scratch/scratch/ucfafsp/SHEFS/R_Library',.libPaths()))

library(CoordinateCleaner)
library(dplyr)
library(dismo)
library(here)
library(lwgeom)
library(raster)
library(readr)
library(sf)
library(sp)
library(stringr)
library(vroom)

#wd <- getwd()
base_dir <- "/lustre/scratch/scratch/ucfafsp/SHEFS/"

source("/lustre/scratch/scratch/ucfafsp/SHEFS/SA_pipeline_functions.R")


###Data loading and prepping
###Data loading and prepping
bm <- st_read(paste0(base_dir,"WWF_Ecoregions/Species_Selection_Area_True.shp"))
ext_sp<-extent(bm)

xmin <- 4
xmax <- 59 # could be 52 if we don't want to include Mauritius 
ymin <- -35
ymax <- 0

min_occ <- 40 # minimum number of occurrence points to be included in the models

bioclim_layers <- c(1, 5, 6, 13, 14)#the number of the bioclim layers to be included as environmental variables - https://worldclim.org/bioclim

ext_occ<-extent(xmin, xmax, ymin, ymax)

###Loading and cropping the environmental data
env_layers <- loadBioclim(path = paste0(base_dir, "CHELSA/"), extension = ".tif")
env_crop <- raster::crop(env_layers, ext_occ)


##Loading in the species names
#files <- list.files(paste0(base_dir, "species_lists"), full.names = TRUE)

sp_names <- read.csv(paste0(base_dir, "/species_lists","/","job_001.csv"),stringsAsFactors = FALSE)
sp_names <- sp_names$sp_name
#sp_names <- "Scotopelia_peli"

### Downloading the gbif data
if(!dir.exists(paste0(base_dir, "points/raw"))){
  dir.create(paste0(base_dir, "/points/raw/"), recursive = TRUE)
}

spxy_out <- lapply(X = sp_names, FUN = gbifData, ext_sp = ext_sp, ext_occ = ext_occ, out_dir = paste0(base_dir, "/points/raw"), min_occ = min_occ)

#### Cleaning the coordinates with coordinate cleaner
if(!dir.exists(paste0(base_dir,"points/cleaned_raw"))){
  dir.create(paste0(base_dir, "points/cleaned_raw/"))
}

###would probably be better to clean multiple species at once because of how coordinate cleaner works  

sp_names_raw <- gsub(".csv", "", list.files(paste0(base_dir,"points/raw")))

sp_names <- sp_names[sp_names %in% sp_names_raw]

lapply(X = sp_names, FUN = cc_wrapper, in_dir = paste0(base_dir, "points/raw"), out_dir = paste0(base_dir, "points/cleaned_raw/"), min_occ = min_occ)

####


####Rarefying points
if(!dir.exists(paste0(base_dir, "points/rarefied"))){
  dir.create(paste0(base_dir,"/points/rarefied/"))
}

ref_map <- env_crop[[1]]
ref_map[!is.na(ref_map)] <- 0   #ref_map should be full of non-1 values

sp_names_cln <- gsub(".csv", "", list.files(paste0(base_dir,"points/cleaned_raw")))

sp_names <- sp_names[sp_names %in% sp_names_cln]

lapply(X = sp_names, FUN = rarefyPoints,in_dir = paste0(base_dir, "points/cleaned_raw/"), out_dir = paste0(base_dir, "points/rarefied/"), ref_map = ref_map, min_occ = min_occ)

#####

#### Extract data for presence points

if(!dir.exists(paste0(base_dir,"environmental/presence/"))){
  dir.create(paste0(base_dir,"environmental/presence/"), recursive = TRUE)
}

sp_names_rar <- gsub(".csv", "", list.files(paste0(base_dir, "points/rarefied/")))

sp_names <- sp_names[sp_names %in% sp_names_rar]

lapply(X= sp_names, FUN = ras_extract, in_dir = paste0(base_dir, "points/rarefied/"), out_dir = paste0(base_dir, "environmental/presence/"), raster_in = env_crop)

####

####Background points creation

ecoreg <-
  sf::st_read(paste0(base_dir, "WWF_Ecoregions/wwf_terr_ecos.shp")) %>% 
  sf::st_crop(.,ext_occ) %>%  ##cropping to the area of intereste=
  dplyr::select(OBJECTID, ECO_NAME) ##just selecting out the columns we're interested in

####can replace pseudoabsence with background

if(!dir.exists(paste0(base_dir, "points/pseudoabsence/"))){
  dir.create(paste0(base_dir,"points/pseudoabsence/"), recursive = TRUE)
}

lapply(X = sp_names, background_sampler, in_dir = paste0(base_dir, "points/rarefied/"), out_dir = paste0(base_dir,"points/pseudoabsence/"), dens_abs = "density", density = 250, type = "pseudoabsence", buffer = 25, polygon = ecoreg)

#### Extract environmental data for pseudoabsence points
if(!dir.exists(paste0(base_dir,"environmental/pseudoabsence/"))){
  dir.create(paste0(base_dir,"environmental/pseudoabsence/"))
}

lapply(X= sp_names, FUN = ras_extract,in_dir = paste0(base_dir,"points/pseudoabsence/"), out_dir = paste0(base_dir,"environmental/pseudoabsence/"), raster_in = env_crop)

#### Subset to species we have all data for
sp_names_env <- gsub("*.csv","",list.files(paste0(base_dir, "environmental/presence"),
                                           pattern = "*.csv"))

sp_names <- sp_names[sp_names %in% sp_names_env]

sp_names <- sp_names[sp_names %in% gsub(".csv","",list.files(paste0(base_dir, "environmental/pseudoabsence/")))]
####

####
if(!dir.exists(paste0(base_dir,"predictions/bioclim/"))){
  dir.create(paste0(base_dir,"predictions/bioclim"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir, "predictions/glm/"))){
  dir.create(paste0(base_dir,"predictions/glm"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir,"predictions/rf/"))){
  dir.create(paste0(base_dir,"predictions/rf"), recursive = TRUE)
}
####

#### Fit Bioclim Models
lapply(
  X = sp_names,
  fitBC,
  pres_dir = paste0(base_dir,"environmental/presence"),
  backg_dir = paste0(base_dir, "environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir = paste0(base_dir,"predictions/bioclim/"),
  eval_out_dir = paste0(base_dir,"evaluation/bioclim/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)
####

#### Fit GLM Models
lapply(
  X = sp_names,
  fitGLM,
  pres_dir = paste0(base_dir,"environmental/presence"),
  backg_dir = paste0(base_dir, "environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,   #need to think about which ones we want to include
  predictors = env_crop,
  pred_out_dir =  paste0(base_dir,"predictions/glm/"),
  eval_out_dir =  paste0(base_dir,"evaluation/glm/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)
####

#### Fit random forest models
lapply(
  X = sp_names,
  fitRF,
  pres_dir =  paste0(base_dir,"environmental/presence"),
  backg_dir =  paste0(base_dir,"environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir =  paste0(base_dir,"predictions/rf/"),
  eval_out_dir =  paste0(base_dir,"evaluation/rf/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)
####

#### Getting evaluation and AUCs out

eval_files <-
  list.files(
    paste0(base_dir,"evaluation/"),
    full.names = TRUE,
    recursive = TRUE,
    pattern = "*eval")


evals_out <- lapply(eval_files, get_eval, threshold = "tss")

eval_df <- do.call(rbind, evals_out)
eval_df$sp_name <- as.character(eval_df$sp_name)

####
if(!dir.exists( paste0(base_dir,"predictions/ensemble/majority_pa"))){
  dir.create(paste0(base_dir,"predictions/ensemble/majority_pa"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir,"predictions/ensemble/weighted"))){
  dir.create(paste0(base_dir,"predictions/ensemble/weighted"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir,"predictions/ensemble/mean"))){
  dir.create(paste0(base_dir,"predictions/ensemble/mean"), recursive = TRUE)
}

####

#### Build ensemble models using either a majority pa or weighted mean strategy

preds <- list.files(paste0(base_dir, "predictions"), full.names = TRUE, recursive = TRUE)

preds <- preds[!grepl("/ensemble/", preds)]

get_sp <- function(name){
  if(sum(stringr::str_detect(preds, name)) >= 1){
    return(name)
  }
}

sp_out <- lapply(sp_names, get_sp)
eval_names <- do.call(rbind, sp_out)  

lapply(X = eval_names, FUN = ensemble_model, eval_df = eval_df, preds = preds, out_dir = paste0(base_dir, "predictions/ensemble/"), method = "majority_pa")

lapply(X = eval_names, FUN = ensemble_model, eval_df = eval_df, preds = preds, out_dir = paste0(base_dir, "predictions/ensemble/"), method = "weighted")


########Plots for eyeballing


if(!dir.exists(paste0(base_dir, "/plots/ensemble/weighted"))){
  dir.create(paste0(base_dir, "/plots/ensemble/weighted"), recursive = TRUE)
}


ggplot_out <- function(sp_name, points_dir, rast_dir, out_dir){
  
  points <- read.csv(paste0(points_dir,"/" ,sp_name,".csv"))
  ras <- raster(paste0(rast_dir, "/", sp_name, "_ensemble.tif"))
  
  ras.p <-  rasterToPoints(ras)
  df <- data.frame(ras.p)
  colnames(df) <- c("Longitude", "Latitude", "Prob")
  
  ggplot(df, aes(y = Latitude, x = Longitude))+
    geom_tile(aes(fill = Prob))+
    geom_point(data = points, aes(x = x, y = y), colour = "red", shape = 4)+
    coord_equal()+
    theme_bw()+
    ggsave(paste0(out_dir,"/",sp_name, ".png"))
}

#ggplot_out(sp_name = sp_name, points_dir = paste0(base_dir, "/points/rarefied/"), rast_dir = paste0(base_dir,"/predictions/ensemble/weighted"), out_dir = paste0(base_dir, "/plots/ensemble/weighted"))

lapply(X = sp_names, points_dir = paste0(base_dir, "/points/rarefied/"), rast_dir = paste0(base_dir,"/predictions/ensemble/weighted"), out_dir = paste0(base_dir, "/plots/ensemble/weighted"))

