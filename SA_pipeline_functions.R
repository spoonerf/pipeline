


gbifData <- function(species, ext) {
  # Include something for if there is nothing in GBIF...
  gen <- strsplit(species, " ")[[1]][1]
  sp <- strsplit(species, " ")[[1]][2]
  
  # count records.
  .count <- dismo::gbif(
    genus = gen,
    species = sp,
    ext = ext,
    geo = TRUE,
    removeZeros = TRUE,
    download = FALSE
  )
  if (.count <= 200000) {
    .xx <- dismo::gbif(
      genus = gen,
      species = sp,
      ext = ext,
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
        output_data <- unique(xx[complete.cases(xx), ])
        output_data <- cbind(species, output_data)
        if (length(output_data) == 3) {
          names(output_data) <- c("species", "x", "y")
        } else {
          colnames(output_data) <- c("species", "x", "y")
        }
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




background_sampler <- function(occ_file, no_pnts) {
  sf_int <- read_csv(occ_file) %>%
    dplyr::select("decimalLongitude", "decimalLatitude") %>%
    distinct() %>%
    st_as_sf(.,
             coords = c("decimalLongitude", "decimalLatitude"),
             crs = 4326) %>%
    st_intersection(., ecoreg)
  
  bkg_ecoreg <- ecoreg %>%
    filter(ECO_NAME %in% sf_int$ECO_NAME) %>%
    st_sample(., size = no_pnts, type = "random")
  
  print(basename(occ_file))
  return(bkg_ecoreg)
}
