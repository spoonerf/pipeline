#source(here::here("combine/FS","pat.R"))
#devtools::install_github("https://github.com/hferg/sdmpl", auth_token = GITHUB_PAT)
#library(sdmpl)

library(CoordinateCleaner)
library(dplyr)
library(here)
library(lwgeom)
library(raster)
library(readr)
library(sf)
library(caret)
library(randomForest)

source("SA_pipeline_functions.R")

# Data Preparation

sp_name <-
  "Martes_martes"  #add species name with underscore between genus and species

#We are only interested in species occurrences in the study region - covered in the extents below

xmin <- -13
xmax <- 2
ymin <- 50
ymax <- 62

ext_occ <- extent(xmin, xmax, ymin, ymax)

min_occ <-
  20 # minimum number of occurrence points needed for the species to be included in the models

bioclim_layers <-
  c(1, 5, 6, 13, 14)#the number of the bioclim layers to be included as environmental variables - https://worldclim.org/bioclim

lapply(bioclim_layers, chelsa_bioclim_get)  #this will download the chelsa bioclim data for the layers that you have selected 

###Loading and cropping the environmental data
env_layers <-
  loadBioclim(path = ("CHELSA/"), extension = ".tif")  #add the path linking to your environment variables. CHELSA bioclim variables from here: http://chelsa-climate.org/downloads/ 
env_crop <- raster::crop(env_layers, ext_occ)


#This code creates this directory  if itdoes not already exist
if (!dir.exists(here::here("points/raw"))) {
  dir.create(here::here("points/raw/"), recursive = TRUE)
}

spxy_out <- gbifData(
  sp_name = sp_name,
  ext_occ = ext_occ,#area over which occurrence points will be downloaded
  out_dir = here::here("points/raw"),#where points will be saved
  min_occ = min_occ
)

##Coordinate Cleaner

if (!dir.exists(here::here("points/cleaned_raw"))) {
  dir.create(here::here("points/cleaned_raw/"))
}

###Cleaning the coordinates using the CoordinateCleaner package - https://cran.r-project.org/web/packages/CoordinateCleaner/CoordinateCleaner.pdf
cc_wrapper(
  sp_name = sp_name,
  in_dir = here::here("points/raw"),
  out_dir = here::here("points/cleaned_raw")
)

###Rarefy Points - so there is only one occurrence point per grid cell
if (!dir.exists(here::here("points/rarefied"))) {
  dir.create(here::here("points/rarefied/"))
}

ref_map <- env_crop[[1]]
ref_map[!is.na(ref_map)] <- 0   #ref_map should be full of non-1 values

rarefyPoints(
  sp_name = sp_name,
  in_dir = here::here("points/cleaned_raw"),
  out_dir = here::here("points/rarefied/"),
  ref_map = ref_map
)

###Extract Data for presence points
if (!dir.exists(here::here("environmental/presence/"))) {
  dir.create(here::here("environmental/presence/"), recursive = TRUE)
}

ras_extract(
  sp_name = sp_name,
  in_dir = here::here("points/rarefied"),
  out_dir = here::here("environmental/presence"),
  raster_in = env_crop
)

###Background Data - create background and pseudoabsence points - need to think about buffer and density values
if (!dir.exists(here::here("points/background/"))) {
  dir.create(here::here("points/background/"))
}

if (!dir.exists(here::here("points/pseudoabsence/"))) {
  dir.create(here::here("points/pseudoabsence/"))
}

wwf_ecoregions_get() # this should download the wwf ecoregions data and put it in the right place for you. 

ecoreg <-
  sf::st_read(here::here("WWF_Ecoregions/wwf_terr_ecos.shp")) %>%
  sf::st_crop(., ext_occ) %>%  ##cropping to the area of intereste=
  dplyr::select(OBJECTID, ECO_NAME) ##just selecting out the columns we're interested in

background_sampler(
  sp_name = sp_name,
  in_dir = here::here("points/rarefied"),
  out_dir = here::here("points/background"),
  dens_abs = "density",
  density = 100,
  type = "background",
  polygon = ecoreg
)

background_sampler(
  sp_name = sp_name,
  in_dir = here::here("points/rarefied"),
  out_dir = here::here("points/pseudoabsence"),
  dens_abs = "density",
  density = 100,
  type = "pseudoabsence",
  buffer = 100,
  polygon = ecoreg
)


###Extract environmental data for background points

if (!dir.exists(here::here("environmental/background/"))) {
  dir.create(here::here("environmental/background/"))
}

ras_extract(
  sp_name = sp_name,
  in_dir = here::here("points/background"),
  out_dir = here::here("environmental/background"),
  raster_in = env_crop
)

###Extract environmental data for pseudoabsence points

if (!dir.exists(here::here("environmental/pseudoabsence/"))) {
  dir.create(here::here("environmental/pseudoabsence/"))
}

ras_extract(
  sp_name = sp_name,
  in_dir = here::here("points/pseudoabsence"),
  out_dir = here::here("environmental/pseudoabsence/"),
  raster_in = env_crop
)


if (!dir.exists(here::here("predictions/bioclim/"))) {
  dir.create(here::here("predictions/bioclim"), recursive = TRUE)
}

if (!dir.exists(here::here("predictions/glm/"))) {
  dir.create(here::here("predictions/glm"), recursive = TRUE)
}

if (!dir.exists(here::here("predictions/rf/"))) {
  dir.create(here::here("predictions/rf"), recursive = TRUE)
}


###Fit Bioclim

fitBC(
  sp_name = sp_name,
  pres_dir = here::here("environmental/presence/"),
  backg_dir = here::here("environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir = here::here("predictions/bioclim/"),
  eval_out_dir = here::here("evaluation/bioclim/"),
  overwrite = TRUE,
  eval = TRUE
)



###Fit GLM

fitGLM(
  sp_name = sp_name,
  pres_dir = here::here("environmental/presence/"),
  backg_dir = here::here("environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir = here::here("predictions/glm/"),
  eval_out_dir = here::here("evaluation/glm/"),
  overwrite = TRUE,
  eval = TRUE
)


###Fit RF

fitRF(
  sp_name = sp_name,
  pres_dir = here::here("environmental/presence/"),
  backg_dir = here::here("environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir = here::here("predictions/rf/"),
  eval_out_dir = here::here("evaluation/rf/"),
  overwrite = TRUE,
  eval = TRUE
)



###Get Evaluations and AUCs

eval_files <-
  list.files(
    here::here("evaluation/"),
    full.names = TRUE,
    recursive = TRUE,
    pattern = paste0("*", sp_name)
  )

evals_out <- lapply(eval_files, get_eval, threshold = "tss")

eval_df <- do.call(rbind, evals_out)
eval_df$sp_name <- as.character(eval_df$sp_name)

eval_df
###Plotting the Model Results

bc_plot <- raster(paste0(here::here("predictions/bioclim/"), "/", sp_name, "_bioclim.tif"))
glm_plot <- raster(paste0(here::here("predictions/glm/"), "/", sp_name, "_glm.tif"))
rf_plot <-  raster(paste0(here::here("predictions/rf/"), "/", sp_name, "_rf.tif"))

cuts=seq(0,1,0.05)#set breaks
pal <- colorRampPalette(c("grey","forestgreen","darkgreen"))

par(mfrow = c(2,2)) #change this to c(1,1) if you don't want all the plots in a pane
plot(bc_plot, main = "Bioclim", breaks=cuts, col = pal(length(cuts)))
plot(glm_plot, main = "GLM", breaks=cuts, col = pal(length(cuts)))
plot(rf_plot, main = "RF", breaks=cuts, col = pal(length(cuts)))


### Plotting the model results with thresholds from eval_df

plot(bc_plot > eval_df$threshold[which(eval_df$model == "bioclim")], main = "Bioclim")
plot(glm_plot > eval_df$threshold[which(eval_df$model == "glm")], main = "GLM")
plot(rf_plot > eval_df$threshold[which(eval_df$model == "rf")], main = "RF")

###Build Ensemble Model

if (!dir.exists(here::here("predictions/ensemble/majority_pa"))) {
  dir.create(here::here("predictions/ensemble/majority_pa"),
             recursive = TRUE)
}

if (!dir.exists(here::here("predictions/ensemble/weighted"))) {
  dir.create(here::here("predictions/ensemble/weighted"))
}

preds <-
  list.files(
    here::here("predictions"),
    full.names = TRUE,
    recursive = TRUE,
    pattern =  paste0("*", sp_name)
  )

preds <- preds[!grepl("/ensemble/", preds)]

#maj_pa takes the thresholded outputs from the models and stacks them on top of eachother.
#For each cell if the majority of the models predict a presence then the output of maj_pa shows a predicted presence

maj_pa <- ensemble_model(
  sp_name = sp_name,
  eval_df = eval_df,
  preds = preds,
  method = "majority_pa",
  out_dir = here::here("predictions/ensemble")
)

#weighted takes the outputs from the models and applies a weighted average 
#to them, so that the best performing model (based on AUC) has a bigger influence on the output.

weighted <- ensemble_model(
  sp_name = sp_name,
  eval_df = eval_df,
  preds = preds,
  method = "weighted",
  out_dir = here::here("predictions/ensemble")
)

xy <-  read.csv(paste0(here::here("points/rarefied/"), "/", sp_name, ".csv"))


plot(maj_pa, main = "Majority Presence/Absence")
points(xy$x, xy$y)

plot(weighted, main = "Weighted Ensemble Model", breaks=cuts, col = pal(10))
points(xy$x, xy$y)

