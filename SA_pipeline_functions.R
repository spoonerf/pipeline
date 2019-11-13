



gbifData <- function(species, ext_sp, ext_occ) {
  # Include something for if there is nothing in GBIF...
  gen <- strsplit(species, " ")[[1]][1]
  sp <- strsplit(species, " ")[[1]][2]
  
  # count records.
  .count <- dismo::gbif(
    genus = gen,
    species = sp,
    ext = ext_sp,
    geo = TRUE,
    removeZeros = TRUE,
    download = FALSE
  )
  if (.count <= 200000) {
    .xx <- dismo::gbif(
      genus = gen,
      species = sp,
      ext = ext_occ,
      geo = TRUE,
      removeZeros = TRUE,
      download = TRUE
    )
    # remove NAs and duplicates
    if (is.null(.xx)) {
      output_data <- NULL
    } else {
      if (all(c("lon", "lat") %in% colnames(.xx))) {
        xx <- cbind(.xx$lon, .xx$lat)
        output_data <-
          matrix(unique(xx[complete.cases(xx),]), ncol = 2)
        output_data <- cbind(species, output_data)
        colnames(output_data) <- c("species", "x", "y")
      } else {
        output_data <- NULL
      }
    }
  } else {
    output_data <- paste0(species, " too many records")
  }
  print(paste(species, "done!"))
  return(output_data)
}


tax_out <- function(id) {
  out <- taxize::classification(id, db = "ncbi")
  sp_out <- NULL
  
  if (!is.na(out[[1]])) {
    profvis::pause(3)
    sp_out <-
      out[[1]][out[[1]]$rank %in% c("kingdom",
                                    "phylum",
                                    "class",
                                    "order",
                                    "family",
                                    "genus",
                                    "species"), ]
    
    sp_out$kingdom <- sp_out$name[sp_out$rank == "kingdom"]
    sp_out$phylum <- sp_out$name[sp_out$rank == "phylum"]
    sp_out$class <- sp_out$name[sp_out$rank == "class"]
    sp_out$order <- sp_out$name[sp_out$rank == "order"]
    sp_out$family <- sp_out$name[sp_out$rank == "family"]
    sp_out$genus <- sp_out$name[sp_out$rank == "genus"]
    sp_out$species <- sp_out$name[sp_out$rank == "species"]
    
    
    sp_out <- sp_out %>%
      dplyr::select(kingdom, phylum, class, order, family, genus, species) %>%
      distinct()
    
  }
  return(sp_out)
}



loadBioclim <- function(path, extension, extent = NULL, prjctn = NULL) {
  bio_layers <- list.files(
    path, pattern = paste0(extension, "$"), full.names = TRUE
  )
   bioclim <- raster::stack(bio_layers)
#  bioclim <- raster::stack()    #don't think this section is needed, have replaced with the line above but maybe I'm missing something
#  for (i in bio_layers) {
#    bioclim <- raster::stack(bioclim, raster[i])
#  }
  if (!is.null(extent)) {
    bioclim <- raster::crop(bioclim, extent)
  }
  if (!is.null(prjctn)) {
    if (class(prjctn) != "CRS") {
      stop("Projection must be an object of class CRS.")
    }
    bioclim <- raster::projectRaster(bioclim, crs = prjctn)
  }
  return(bioclim)
}





rarefyPoints <- function(ref_map, pnts) {
  cells <- raster::cellFromXY(ref_map, pnts)
  pres_cells <- ref_map
  pres_cells[unique(cells)] <- 1
  rarefied_presence <- raster::rasterToPoints(
    pres_cells,
    fun = function(x) {
      x == 1
    }
  )[, c(1, 2)]
  rarefied_presence <- sp::SpatialPoints(rarefied_presence)
  sp::proj4string(rarefied_presence) <- sp::proj4string(ref_map)
  return(rarefied_presence)
}



background_sampler <- function(in_file, no_pnts, out_file) {
 
  sp_name <- gsub(".csv", "", basename(in_file))
  
   sf_int <- read_csv(in_file) %>%
    dplyr::select("x", "y") %>%
    distinct() %>%
    st_as_sf(.,
             coords = c("x", "y"),
             crs = 4326) %>%
    st_intersection(., ecoreg)
  
  bkg_ecoreg <- ecoreg %>%
    filter(ECO_NAME %in% sf_int$ECO_NAME) %>%
    st_sample(., size = no_pnts, type = "random")
  
  tibb <- as_tibble(bkg_ecoreg)
  
  sep_df <- tibb %>% 
    mutate(x = unlist(map(tibb$geometry,1)),
           y = unlist(map(tibb$geometry, 2))) %>% 
    dplyr::select(x, y)
  
  df_out <- data.frame(sp_name, sep_df)
  
  write.csv(x = df_out,
            file = paste0(out_file),
            row.names = FALSE)
  
  print(basename(in_file))
  return(df_out)
}




fitBC <- function(pres_file_in, backg_file_in, predictor_names,predictors, filename, overwrite, threads = 4) {
  
  predictor_names <- stringr::str_pad(predictor_names, 2, pad = "0")
  
  predictor_names <- paste0("CHELSA_bio10_", predictor_names)
  
  curr <- raster::dropLayer(predictors, which(!names(predictors) %in% predictor_names)
  )
  
  sp_name <- gsub(".csv", "",basename(pres_file_in))
  
  pres_pts <- read.table(pres_file_in, stringsAsFactors = FALSE, colClasses = c("NULL", rep("numeric",2), rep("NULL", 19)), sep = ",", header = TRUE)
  
  backg_pts <- 
    
  cat("Fitting bioclim model...\n")
  bc <- dismo::bioclim(curr, pres_pts)
  cat("Done.\n")
  cat("...\n")
  cat("Evaluating bioclim model...\n")
 
  kfolds_p <- dismo::kfold(pa_raster$combined_presence, 4)
  kfolds_a <- dismo::kfold(pa_raster$background, 4)
    if (.Platform$OS.type == "unix") {
      cl <- parallel::makeForkCluster(4)
    } else {
      cl <- parallel::makeCluster(4)
    }
    parallel::clusterExport(
      cl, varlist = c("pa_raster", "kfolds_p", "kfolds_a", "curr"),
      envir = environment()
    )
    aucs <- parallel::clusterApply(cl, 1:4, function(x) {
      clustEvalPa(x, pa_raster, kfolds_p, kfolds_a, curr, mod = "bioclim")
    })
    parallel::stopCluster(cl)
    cat("Done.\n")
    cat("...\n")
    thresholds <- getThresholds(aucs)
  }
  
  cat("Predicting from bioclim model...\n")
  if (.Platform$OS.type == "unix") {
    cl <- parallel::makeForkCluster(threads)
  } else {
    cl <- parallel::makeCluster(threads)
    parallel::clusterExport(
      cl, varlist = c("predictors", "bc"), envir = environment()
    )
  }
  
  res <- parallel::clusterApply(
    cl, 1:length(predictors), function(x) {
      dismo::predict(predictors[[x]], bc)
    }
  )
  parallel::stopCluster(cl)
  res <- raster::stack(res)
  names(res) <- names(predictors)
  cat("Done.\n")
  cat("...\n")
  cat("Writing bioclim predictions...\n")
  out_file <- paste0(filename, "bioclim.grd")
  raster::writeRaster(res, filename = out_file, format = "raster",
                      overwrite = overwrite)
  gc()
  cat("Done.\n")
  cat("...\n")
  if (eval) {
    return(list(aucs = aucs, thresholds = thresholds))
  }
}

