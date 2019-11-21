




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
          matrix(unique(xx[complete.cases(xx), ]), ncol = 2)
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
                                    "species"),]
    
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


ras_extract <- function(in_file, out_file) {
  df <- vroom::vroom(in_file)
  xy <- SpatialPointsDataFrame(matrix(c(df$x, df$y), ncol = 2), df)
  ras_ext <- raster::extract(env_crop, xy)
  pres_ext <- data.frame(df, ras_ext)
  write.csv(x = pres_ext,
            file = paste0(out_file),
            row.names = FALSE)
  
}


loadBioclim <-
  function(path,
           extension,
           extent = NULL,
           prjctn = NULL) {
    bio_layers <- list.files(path,
                             pattern = paste0(extension, "$"),
                             full.names = TRUE)
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


####
###@dens_abs - whether you want to use sampling based on density or an absolute number e.g. 10000
###@density - what that density is in terms of per km. e.g. 1000 would be 1 point per 1000km2
###@no_pnts - number of points if using an absolute number
###@type - type of points to use, either background or pseudoabsence. Background will distribute points randomly across the polygons in which points fall. 
### Pseudo-absence will also do this but exclude points less than @buffer km from a presence point
###@buffer - the size of the buffer used in type = "pseudoabsence"
###@polygon - the polygon shapefile used as the background - here we use ecoregions but could be any shapefile with multiple polygons

background_sampler <- function(sp_name, in_dir, out_dir, dens_abs = "absolute", 
                               density = NULL, no_pnts = NULL, type = "background", buffer = NULL, polygon = NULL) 
{
  
  in_file <- list.files(in_dir, full.names = TRUE)[grepl(sp_name, list.files(in_dir))]
  
  sf_int <- read_csv(paste0(in_dir, "/", sp_name, ".csv")) %>%
    dplyr::select("x", "y") %>%
    dplyr::distinct() %>%
    sf::st_as_sf(.,
             coords = c("x", "y"),
             crs = 4326) %>%
    sf::st_intersection(., polygon)
  
  bkg_polygon <- polygon %>%
    dplyr::filter(ECO_NAME %in% sf_int$ECO_NAME)
  
  if (dens_abs == "density"){
    no_pnts <- round(as.numeric(sum(st_area(bkg_polygon)))/(1000000*density))   
  }
  
  if (type == "background"){
    points_out <- bkg_polygon %>% 
      sf::st_sample(., size = no_pnts, type = "random")  
  }
  
  if (type == "pseudoabsence"){
    
    diss_bkg_polygon <- sf::st_union(bkg_polygon)
    sf_int_trans <- st_transform(sf_int, 54030) #robinson projection
    buff_pnts <- sf::st_buffer(sf_int_trans, buffer*1000)  
    buff_pnts <- st_transform(buff_pnts, 4326)
    buff_pnts <- sf::st_union(buff_pnts)
    diff_bkg_polygon <- sf::st_difference(diss_bkg_polygon, buff_pnts)  
    points_out <- diff_bkg_polygon %>% 
      sf::st_sample(., size = no_pnts, type = "random")
    
  }
  
  tibb <- as_tibble(points_out)
  
  sep_df <- tibb %>%
    mutate(x = unlist(purrr::map(tibb$geometry, 1)),
           y = unlist(purrr::map(tibb$geometry, 2))) %>%
    dplyr::select(x, y)
  
  df_out <- data.frame(sp_name, sep_df)
  
  write.csv(x = df_out,
            file = paste0(out_dir, "/", sp_name, ".csv"),
            row.names = FALSE)
  
  print(basename(in_file))
  return(df_out)
}



clustEvalPa <-
  function(num,
           pres_pts,
           backg_pts,
           kfolds_p,
           kfolds_a,
           curr,
           mod) {
    pres_train <- pres_pts[kfolds_p != num,]
    pres_test <- pres_pts[kfolds_p == num,]
    backg_test <- backg_pts[kfolds_a == num,]
    if (mod == "bioclim") {
      .m <- dismo::bioclim(curr, pres_train)
    } else if (mod == "domain") {
      .m <- dismo::domain(curr, pres_train)
    }
    e <- dismo::evaluate(pres_test, backg_test, .m, curr)
    return(e)
  }


clustEvalSdm <- function(num, sdm_set, kfolds, model, mod) {
  train <- sdm_set[kfolds != num,]
  test_p <- sdm_set[kfolds == num & sdm_set[, "pb"] == 1,]
  test_a <- sdm_set[kfolds == num & sdm_set[, "pb"] == 0,]
  if (mod == "glm") {
    .m <- stats::glm(stats::formula(model),
                     data = train,
                     family = binomial(link = "logit"))
  } else if (mod == "rf") {
    .m <-
      suppressWarnings(randomForest::randomForest(model, data = train))
  }
  
  e <- dismo::evaluate(test_p, test_a, .m)
  e
}



getThresholds <- function(aucs) {
  thresholds <- vector(mode = "list", length = 5)
  names(thresholds) <- c("spec_sens",
                         "no_omission",
                         "prevalence",
                         "equal_sens_spec",
                         "tss")
  thresholds[[1]] <- sapply(aucs, function(x)
    dismo::threshold(x, "spec_sens"))
  thresholds[[2]] <- sapply(aucs, function(x)
    dismo::threshold(x, "no_omission"))
  thresholds[[3]] <- sapply(aucs, function(x)
    dismo::threshold(x, "prevalence"))
  thresholds[[4]] <- sapply(aucs, function(x)
    dismo::threshold(x, "equal_sens_spec"))
  thresholds[[5]] <- sapply(aucs, function(x)
    tssCalc(x))
  return(thresholds)
}


tssCalc <- function(eval) {
  res <- data.frame(threshold = eval@t,
                    tss = apply(eval@confusion, 1, function(x) {
                      cm <- t(matrix(rev(x), nrow = 2))
                      dimnames(cm) <- list(pred = c("0", "1"),
                                           obs = c("0", "1"))
                      class(cm) <- "table"
                      sens <- caret::sensitivity(cm)
                      spec <- caret::specificity(cm)
                      tss <- sens + spec - 1
                      return(tss)
                    }))
  thresh <- res$threshold[which.max(res$tss)]
  return(thresh)
}

fitBC <-
  function(sp_name,
           pres_dir,
           backg_dir,
           predictor_names,
           predictors,
           pred_out_dir,
           eval_out_dir,
           overwrite,
           threads = 4,
           eval = TRUE) {
    
    print(sp_name)
    
    predictor_names <-
      stringr::str_pad(predictor_names, 2, pad = "0")
    
    CHELSA_predictor_names <-
      paste0("CHELSA_bio10_", predictor_names)
    
    curr <-
      raster::dropLayer(predictors, which(!names(predictors) %in% CHELSA_predictor_names))
    
    pres_pts <-
      as.matrix(data.table::fread(paste0(pres_dir, "/", sp_name, ".csv"), select = c("x", "y")), ncol = 2)
    
    if (nrow(pres_pts) < 10) {
      cat("Fewer than 10 data points - cannot fit model!")
    } else {
      backg_pts <-
        as.matrix(data.table::fread(paste0(backg_dir, "/", sp_name, ".csv"), select = c("x", "y")))
      
      cat("Fitting bioclim model...\n")
      bc <- dismo::bioclim(curr, pres_pts)
      cat("Done.\n")
      cat("...\n")
      cat("Evaluating bioclim model...\n")
      
      kfolds_p <- dismo::kfold(pres_pts, 4)
      kfolds_a <- dismo::kfold(backg_pts, 4)
      
      if (.Platform$OS.type == "unix") {
        cl <- parallel::makeForkCluster(threads)
      } else {
        cl <- parallel::makeCluster(threads)
      }
      parallel::clusterExport(
        cl,
        varlist = c(
          "pres_pts",
          "backg_pts",
          "kfolds_p",
          "kfolds_a",
          "curr",
          "clustEvalPa"
        ),
        envir = environment()
      )
      aucs <- parallel::clusterApply(cl, 1:4, function(x) {
        clustEvalPa(x, pres_pts, backg_pts, kfolds_p, kfolds_a, curr, mod = "bioclim")
      })
      parallel::stopCluster(cl)
      cat("Done.\n")
      cat("...\n")
      thresholds <- getThresholds(aucs)
      
      cat("Predicting from bioclim model...\n")
      res <- dismo::predict(curr, bc)
      cat("Done.\n")
      cat("...\n")
      cat("Writing bioclim predictions...\n")
      out_file <- paste0(sp_name, "_bioclim.tif")
      
      if (!dir.exists(pred_out_dir)) {
        dir.create(pred_out_dir)
      }
      
      raster::writeRaster(
        res,
        filename = paste(pred_out_dir, out_file, sep = "/"),
        format = "GTiff",
        overwrite = overwrite
      )
      gc()
      cat("Done.\n")
      cat("...\n")
      if (!dir.exists(eval_out_dir)) {
        dir.create(eval_out_dir, recursive = TRUE)
      }
      
      if (eval) {
        evals <- list(sp_name = sp_name, model = "bioclim", aucs = aucs, thresholds = thresholds)
        #evals <- data.frame(sp_name = sp_name, model = "bioclim", aucs = unlist(aucs), thresholds = unlist(thresholds))
        save(evals, file = paste0(eval_out_dir, "/", sp_name, "_bioclim_eval.RDA"))
        #save(evals, file = paste0(eval_out_dir, "/", sp_name, "_bioclim_eval.csv"))      
      
        }
    }

  }


fitGLM <-
  function(sp_name,
           pres_dir,
           backg_dir,
           predictor_names,
           predictors,
           pred_out_dir,
           eval_out_dir,
           overwrite,
           threads = 4,
           eval = TRUE) {
    print(sp_name)
    
    predictor_names <- stringr::str_pad(predictor_names, 2, pad = "0")
    
    predictor_names <- paste0("CHELSA_bio10_", predictor_names)
    
    curr <-
      raster::dropLayer(predictors, which(!names(predictors) %in% predictor_names))
    
    model <-
      stats::formula(paste("pb ~", paste(predictor_names, collapse = "+")))
    
    pres <-
      data.frame(pb = 1,
                 data.table::fread(paste0(pres_dir, "/", sp_name, ".csv"), select = predictor_names))
    
    
    if (nrow(pres) < 10) {
      cat("Fewer than 10 data points - cannot fit model!")
      
    } else {
      backg <-
        data.frame(pb = 0,
                   data.table::fread(paste0(backg_dir, "/", sp_name, ".csv"), select = predictor_names))
      
      sdm_set <- rbind(pres, backg)
      
      cat("Fitting GLM...\n")
      m <- stats::glm(stats::formula(model),
                      data = sdm_set,
                      family = binomial(link = "logit"))
      cat("Done.\n")
      cat("...\n")
      if (eval) {
        cat("Evaluating GLM...\n")
        kfolds <- dismo::kfold(sdm_set, 4)
        if (.Platform$OS.type == "unix") {
          cl <- parallel::makeForkCluster(4)
        } else {
          cl <- parallel::makeCluster(4)
          parallel::clusterExport(
            cl,
            varlist = c("sdm_set", "kfolds", "model", "clustEvalSdm"),
            envir = environment()
          )
          parallel::clusterCall(cl, function()
            library(dismo))
        }
        aucs <- parallel::clusterApply(cl, 1:4, function(x) {
          clustEvalSdm(x, sdm_set, kfolds, model, mod = "glm")
        })
        parallel::stopCluster(cl)
        cat("Done.\n")
        cat("...\n")
        thresholds <- getThresholds(aucs)
      }
      
      cat("Predicting from GLM...\n")
      res <- dismo::predict(curr, m)
      res <-
        raster:::calc(
          res,
          fun = function(x) {
            exp(x) / (1 + exp(x))
          }
        ) #backtransforming from logit space
      cat("Done.\n")
      cat("...\n")
      cat("Writing GLM predictions...\n")
      out_file <- paste0(sp_name, "_glm.tif")
      
      if (!dir.exists(pred_out_dir)) {
        dir.create(pred_out_dir, recursive = TRUE)
      }
      
      raster::writeRaster(
        res,
        filename = paste(pred_out_dir, out_file, sep = "/"),
        format = "GTiff",
        overwrite = overwrite
      )
      gc()
      cat("Done.\n")
      cat("...\n")
      
      if (!dir.exists(eval_out_dir)) {
        dir.create(eval_out_dir)
      }
      
      if (eval) {
        evals = data.frame(sp_name = sp_name, model = "glm",aucs = aucs, thresholds = thresholds)
        save(evals, file = paste0(eval_out_dir, "/", sp_name, "_glm_eval.csv"))      
      }
    }
    
    
  }




fitRF <-   function(sp_name,
                    pres_dir,
                    backg_dir,
                    predictor_names,
                    predictors,
                    pred_out_dir,
                    eval_out_dir,
                    overwrite,
                    threads = 4,
                    eval = TRUE) {
  print(sp_name)
  
  predictor_names <- stringr::str_pad(predictor_names, 2, pad = "0")
  
  predictor_names <- paste0("CHELSA_bio10_", predictor_names)
  
  curr <-
    raster::dropLayer(predictors, which(!names(predictors) %in% predictor_names))
  
  model <-
    stats::formula(paste("pb ~", paste(predictor_names, collapse = "+")))
  
  pres <-
    data.frame(pb = 1,
               data.table::fread(paste0(pres_dir, "/", sp_name, ".csv"), select = predictor_names))
  
  
  if (nrow(pres) < 10) {
    cat("Fewer than 10 data points - cannot fit model!")
    
  } else {
    backg <-
      data.frame(pb = 0,
                 data.table::fread(paste0(backg_dir, "/", sp_name, ".csv"), select = predictor_names))
    
    sdm_set <- rbind(pres, backg)
    
    sdm_set <-
      sdm_set[complete.cases(sdm_set),]  #should remove this line and just have no NAs background data
    
    cat("Fitting random forest model...\n")
    m <-
      suppressWarnings(randomForest::randomForest(model, data = sdm_set))
    cat("Done.\n")
    cat("...\n")
    if (eval) {
      cat("Evaluating Random Forest model...\n")
      kfolds <- dismo::kfold(sdm_set, 4)
      if (.Platform$OS.type == "unix") {
        cl <- parallel::makeForkCluster(4)
      } else {
        cl <- parallel::makeCluster(4)
      }
      parallel::clusterExport(
        cl,
        varlist = c("sdm_set", "kfolds", "model", "clustEvalSdm"),
        envir = environment()
      )
      aucs <- parallel::clusterApply(cl, 1:4, function(x) {
        clustEvalSdm(x, sdm_set, kfolds, model, mod = "rf")
      })
      parallel::stopCluster(cl)
      cat("Done.\n")
      cat("...\n")
      thresholds <- getThresholds(aucs)
    }
    
    cat("Predicting from random forest...\n")
    res <- dismo::predict(curr, m)
    cat("Done.\n")
    cat("...\n")
    cat("Writing random forest predictions...\n")
    out_file <- paste0(sp_name, "_rf.tif")
    
    if (!dir.exists(pred_out_dir)) {
      dir.create(pred_out_dir, recursive = TRUE)
    }
    
    raster::writeRaster(
      res,
      filename = paste(pred_out_dir, out_file, sep = "/"),
      format = "GTiff",
      overwrite = overwrite
    )
    gc()
    cat("Done.\n")
    cat("...\n")
    
    if (!dir.exists(eval_out_dir)) {
      dir.create(eval_out_dir, recursive = TRUE)
    }
    
    if (eval) {
      
      evals = data.frame(sp_name = sp_name, model = "rf", aucs = aucs, thresholds = thresholds)
      save(evals, file = paste0(eval_out_dir, "/", sp_name, "_rf_eval.csv"))      
    }
    
  }
}
