lf <- list.files("H:/SHEFS_joint_project/sdm_pipeline_output/", pattern = "*.tar.gz$", full.names = TRUE)

#dir.create("H:/SHEFS_joint_project/sdm_pipeline_output/exdir")

cmd_make <- function(in_file){

  out_folder <- gsub(".tar.gz", "",in_file)
  dir.create(out_folder)
  cmd_call <- noquote(paste0("tar -xvzf ",in_file, " -C ", out_folder, " &"))  
  return(cmd_call)
}

cmd_out <- lapply(lf,cmd_make)

cmd_all <- do.call(rbind, cmd_out)

write.table(cmd_all, "cmd_call.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)


