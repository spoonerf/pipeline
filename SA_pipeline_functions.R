
gbifData <- function(species, ext) {
  # Include something for if there is nothing in GBIF...
  gen <- strsplit(species, " ")[[1]][1]
  sp <- strsplit(species, " ")[[1]][2]
  
  # count records.
  .count <- dismo::gbif(genus = gen, 
                        species = sp, 
                        ext = ext,
                        geo = TRUE,
                        removeZeros = TRUE,
                        download = FALSE) 
  if (.count <= 200000) {
    .xx <- dismo::gbif(genus = gen, 
                       species = sp, 
                       ext = ext,
                       geo = TRUE,
                       removeZeros = TRUE,
                       download = TRUE)
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

