## scenario_folder.m ###########################################################
scenario_folder <- function(scn) {
  
  # Delete and recreate folders for Figures and Results.
  unlink(    paste0("../R/Results-temp/", scn, "/"),               recursive = TRUE) # Delete the folder and its contents
  dir.create(paste0("../R/Results-temp/", scn, "/"),               recursive = TRUE) # Create the new folder
  dir.create(paste0("../R/Results-temp/", scn, "/RData/"),         recursive = TRUE) # Create the new folder
  dir.create(paste0("../R/Results-temp/", scn, "/Figures/"),       recursive = TRUE) # Create the new folder
}
