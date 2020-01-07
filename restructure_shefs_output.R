dirs <- list.dirs("H:/SHEFS_joint_project/sdm_pipeline_output/unzipped_jobs/", recursive = FALSE)
dirs <- dirs[grepl("job", dirs)]
#jobs <- basename(dirs)
pattern <- "/home/ucfafsp/Scratch/SHEFS/"


eval_dir <- "H:/SHEFS_joint_project/sdm_pipeline_output/evaluation/"
points_dir <- "H:/SHEFS_joint_project/sdm_pipeline_output/points/"
pred_dir <- "H:/SHEFS_joint_project/sdm_pipeline_output/predictions/"
maj_dir <- "H:/SHEFS_joint_project/sdm_pipeline_output/ensemble/majority_pa"
wei_dir <- "H:/SHEFS_joint_project/sdm_pipeline_output/ensemble/weighted"

dir.create(eval_dir)
dir.create(points_dir)
dir.create(pred_dir)
dir.create(maj_dir, recursive = TRUE)
dir.create(wei_dir, recursive = TRUE)

get_files <- function(in_dir){
  
  eval_files <- list.files(paste0(in_dir, pattern, basename(in_dir), "/evaluation"), recursive = TRUE, full.names = TRUE)
  points_files <- list.files(paste0(in_dir, pattern, basename(in_dir), "/points/rarefied"), recursive = TRUE, full.names = TRUE)
  preds_files <- list.files(paste0(in_dir, pattern, basename(in_dir), "/predictions"), recursive = TRUE, full.names = TRUE)
  maj_files <- preds_files[grepl("majority_pa",  preds_files)]
  wei_files <- preds_files[grepl("weighted",  preds_files)]
  
  file.copy(eval_files, paste0(eval_dir, basename(eval_files)), overwrite = TRUE)
  file.copy(points_files, paste0(points_dir, basename(points_files)), overwrite = TRUE)
  file.copy(preds_files, paste0(pred_dir, basename(preds_files)), overwrite = TRUE)
  file.copy(maj_files, paste0(maj_dir, "/", basename(maj_files)), overwrite = TRUE)
  file.copy(wei_files, paste0(wei_dir, "/",basename(wei_files)), overwrite = TRUE)
  print(in_dir)
  
  }

##need to untar 1:20 again
lapply(dirs[158], get_files)

