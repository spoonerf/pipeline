######################################################################################################
### creating an R script for each job - essentially changing the lines which refer to the job name ###
######################################################################################################

files <- list.files(here::here("species_lists/"))

if(!dir.exists(here::here("r_jobs"))){
  dir.create(here::here("r_jobs"))
}
               
           
for(file in files){
  
  ptc <- readLines(here::here("r_jobs", "job_001.R"))  #Read in the template file
  job <- gsub(".csv", "", file)
  ptc[17] <- gsub("job_001", job, ptc[17])
  ptc[46] <- gsub("job_001.csv", file, ptc[46])
  file_out <- gsub(".csv", ".R", file)
  writeLines(ptc, here::here("r_jobs",file_out))  #Use this as the file name and write it to the output directory
  
}

################################################################################
### R file to delete files created during job - after it has been tarred up. ###
################################################################################

if(!dir.exists(here::here("delete_jobs"))){
  dir.create(here::here("delete_jobs"))
}


for(file in files){
  
  ptc <- readLines(here::here("delete_jobs", "job_001.R"))  #Read in the template file
  job <- gsub(".csv", "", file)
  ptc[1] <- gsub("job_001", job, ptc[1])
  ptc[3] <- gsub("job_001.csv", file, ptc[3])
  file_out <- gsub(".csv", ".R", file)
  writeLines(ptc, here::here("delete_jobs",file_out))  #Use this as the file name and write it to the output directory
  
}

########################################### 
### bash_scripts to submit all the jobs ###
###########################################

if(!dir.exists(here::here("bash_scripts"))){
  dir.create(here::here("bash_scripts"))
}

files <- gsub("\\..*","",files)

for(file in files){
  
  ptc <- readLines(here::here("bash_scripts","job_001.sh"))  #Read in the template file
  ptc[5] <- gsub("job_001", file, ptc[5])
  ptc[23] <- gsub("job_001", file, ptc[23])
  ptc[24] <- gsub("job_001", file, ptc[24])
  ptc[25] <- gsub("job_001", file, ptc[25])
  writeLines(ptc, here::here("bash_scripts",paste0(file, ".sh")), sep = "\n")  #Use this as the file name and write it to the output directory
  
}

####################################################
### Changing the EOL formatting from dos to unix ###
### Code to be run in cmd on windows             ###
####################################################

unix_convert <- "tr -d '\015' <job_001.sh > ujob_001.sh"

convert_out <- NULL

for (file in files){
  
  ptc <- "tr -d '\015' <job_001.sh > ujob_001.sh"
  ptc <- gsub("job_001", file, ptc)
  
  convert_out <- rbind(convert_out, ptc)

    }

write.table(convert_out, file = "convert_unix.txt", row.names = FALSE)

