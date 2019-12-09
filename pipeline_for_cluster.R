library(CoordinateCleaner)
library(dplyr)
library(here)
library(lwgeom)
library(raster)
library(readr)
library(sf)

source("SA_pipeline_functions.R")


###Data loading and prepping
bm<-st_read(here::here("WWF_Ecoregions/Species_Selection_Area_True.shp"))
ext_sp<-extent(bm)

xmin <- 4
xmax <- 59 # could be 52 if we don't want to include Mauritius 
ymin <- -35
ymax <- 0

min_occ <- 40 # minimum number of occurrence points to be included in the models

bioclim_layers <- c(1, 5, 6, 13, 14)#the number of the bioclim layers to be included as environmental variables - https://worldclim.org/bioclim

ext_occ<-extent(xmin, xmax, ymin, ymax)

###Loading and cropping the environmental data
env_layers <- loadBioclim(path = here::here("CHELSA/"), extension = ".tif")
env_crop <- raster::crop(env_layers, ext_occ)


##Loading in the species names
sp_names <- gsub(".csv", "", list.files(here::here("points/raw")))

#sp_names <- "Pycnonotus_tricolor"

### Downloading the gbif data
if(!dir.exists(here::here("points/raw"))){
  dir.create(here::here("points/raw/"), recursive = TRUE)
}

spxy_out <- lapply(X = sp_names[1:100], FUN = gbifData, ext_sp = ext_sp, ext_occ = ext_occ, out_dir = here::here("points/raw"), min_occ = min_occ)
##15:48 start

#### Cleaning the coordinates with coordinate cleaner
if(!dir.exists(here::here("points/cleaned_raw"))){
  dir.create(here::here("points/cleaned_raw/"))
}

###would probably be better to clean multiple species at once because of how coordinate cleaner works  

lapply(X = sp_names, FUN = cc_wrapper, in_dir = here::here("points/raw"), out_dir = here::here("points/cleaned_raw/"), min_occ = min_occ)

####


####Rarefying points
if(!dir.exists(here::here("points/rarefied"))){
  dir.create(here::here("points/rarefied/"))
}

ref_map <- env_crop[[1]]
ref_map[!is.na(ref_map)] <- 0   #ref_map should be full of non-1 values

sp_names <- gsub(".csv", "", list.files(here::here("points/cleaned_raw")))

lapply(X = sp_names, FUN = rarefyPoints,in_dir = here::here("points/cleaned_raw/"), out_dir = here::here("points/rarefied/"), ref_map = ref_map, min_occ = min_occ)

#####

#### Extract data for presence points

if(!dir.exists(here::here("environmental/presence/"))){
  dir.create(here::here("environmental/presence/"), recursive = TRUE)
}

sp_names <- gsub(".csv", "", list.files(here::here("points/rarefied/")))

lapply(X= sp_names, FUN = ras_extract,in_dir = here::here("points/rarefied/"), out_dir = here::here("environmental/presence/"), raster_in = env_crop)

####

####Background points creation

sp_names <- gsub("*.csv","",list.files(here::here("points/rarefied"),
                                       pattern = "*.csv"))
ecoreg <-
  sf::st_read(here::here("WWF_Ecoregions/wwf_terr_ecos.shp")) %>% 
  sf::st_crop(.,ext_occ) %>%  ##cropping to the area of intereste=
  dplyr::select(OBJECTID, ECO_NAME) ##just selecting out the columns we're interested in

####can replace pseudoabsence with background

if(!dir.exists(here::here("points/pseudoabsence/"))){
  dir.create(here::here("points/pseudoabsence/"), recursive = TRUE)
}

lapply(X = sp_names, background_sampler, in_dir = here::here("points/rarefied/"), out_dir = here::here("points/pseudoabsence/"), dens_abs = "density", density = 250, type = "pseudoabsence", buffer = 25, polygon = ecoreg)

#### Extract environmental data for pseudoabsence points
if(!dir.exists(here::here("environmental/pseudoabsence/"))){
  dir.create(here::here("environmental/pseudoabsence/"))
}

lapply(X= sp_names, FUN = ras_extract,in_dir = here::here("points/pseudoabsence/"), out_dir = here::here("environmental/pseudoabsence/"), raster_in = env_crop)

#### Subset to species we have all data for
sp_names <- gsub("*.csv","",list.files(here::here("environmental/presence"),
                                       pattern = "*.csv"))
sp_names <- sp_names[sp_names %in% gsub(".csv","",list.files(here::here("environmental/pseudoabsence/")))]
####

####
if(!dir.exists(here::here("predictions/bioclim/"))){
  dir.create(here::here("predictions/bioclim"), recursive = TRUE)
}

if(!dir.exists(here::here("predictions/glm/"))){
  dir.create(here::here("predictions/glm"), recursive = TRUE)
}

if(!dir.exists(here::here("predictions/rf/"))){
  dir.create(here::here("predictions/rf"), recursive = TRUE)
}
####

#### Fit Bioclim Models
lapply(
  X = sp_names,
  fitBC,
  pres_dir = here::here("environmental/presence"),
  backg_dir = here::here("environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir = here::here("predictions/bioclim/"),
  eval_out_dir = here::here("evaluation/bioclim/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)
####

#### Fit GLM Models
lapply(
  X = sp_names,
  fitGLM,
  pres_dir = here::here("environmental/presence"),
  backg_dir = here::here("environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,   #need to think about which ones we want to include
  predictors = env_crop,
  pred_out_dir = here::here("predictions/glm/"),
  eval_out_dir = here::here("evaluation/glm/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)
####

#### Fit random forest models
lapply(
  X = sp_names,
  fitRF,
  pres_dir = here::here("environmental/presence"),
  backg_dir = here::here("environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir = here::here("predictions/rf/"),
  eval_out_dir = here::here("evaluation/rf/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)
####

#### Getting evaluation and AUCs out

eval_files <-
  list.files(
    here::here("evaluation/"),
    full.names = TRUE,
    recursive = TRUE,
    pattern = "*eval")


evals_out <- lapply(eval_files, get_eval, threshold = "tss")

eval_df <- do.call(rbind, evals_out)
eval_df$sp_name <- as.character(eval_df$sp_name)

####
if(!dir.exists(here::here("predictions/ensemble/majority_pa"))){
  dir.create(here::here("predictions/ensemble/majority_pa"), recursive = TRUE)
}

if(!dir.exists(here::here("predictions/ensemble/weighted"))){
  dir.create(here::here("predictions/ensemble/weighted"))
}

if(!dir.exists(here::here("predictions/ensemble/mean"))){
  dir.create(here::here("predictions/ensemble/mean"))
}

####

#### Build ensemble models using either a majority pa or weighted mean strategy

preds <- list.files(here::here("predictions"), full.names = TRUE, recursive = TRUE)

preds <- preds[!grepl("/ensemble/", preds)]

get_sp <- function(name){
  if(sum(stringr::str_detect(preds, name)) >= 1){
    return(name)
  }
}

sp_out <- lapply(sp_names, get_sp)
eval_names <- do.call(rbind, sp_out)  

lapply(X = eval_names, FUN = ensemble_model, eval_df = eval_df, preds = preds, method = "majority_pa")

lapply(X = eval_names, FUN = ensemble_model, eval_df = eval_df, preds = preds, method = "weighted")








