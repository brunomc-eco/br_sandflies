# ENM output area calculations
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)
library(rgdal)


# load data and set values ------------------------------------------------

# name the folder where previous final models were saved
run_name <- c("./outputs/models_hist_100km/")

# scenario and gcm names
gcm_names <- c("BCC-CSM2-MR", "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1",
               "IPSL-CM6A-LR", "MIROC6", "MIROC-ES2L", "MRI-ESM2-0")

ssp_names <- c("ssp126_2041-2060", "ssp245_2041-2060", "ssp585_2041-2060",
               "ssp126_2081-2100", "ssp245_2081-2100", "ssp585_2081-2100")


# species names
study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")


# mask for cropping results to Brazil
br_mask <- readOGR("C:/layers/vector/GADM_Continental South America/BRA_adm/BRA_adm0.shp")



# suitability differences -------------------------------------------------

for(i in 1:length(study_sp)){
  
  # load present ensemble
  present <- list.files(path = paste0(run_name, study_sp[i], "/present/ensemble"),
                        pattern = "mean",
                        full.names = TRUE) %>%
    raster() 
  
  for(g in 1:length(ssp_names)){
   
     # load future ensemble for this scenario
    future <- list.files(path = paste0(run_name, study_sp[i], "/future_consensus"),
                         pattern = paste0(ssp_names[g], "_mean"),
                         full.names = TRUE) %>%
      raster() 
    
    diff = future - present
    plot(diff, main = paste(study_sp[i], "diff", ssp_names[g]))
    
    # save diff output
    if (!file.exists(paste0(run_name, study_sp[i], "/diffs"))) {
      dir.create(paste0(run_name, study_sp[i], "/diffs"))
    }
    
    writeRaster(diff, filename = paste0(run_name, study_sp[i], "/diffs/", "diff_", 
                                        ssp_names[g], study_sp[i], ".tif"), 
                overwrite = TRUE)
  }
  
}


# differences in consensus ------------------------------------------------

for(i in 1:length(study_sp)){
  
  # load species calibration area
  m <- readOGR(dsn = paste0(run_name, study_sp[i]),
               layer = paste0("03_calibration_area_", study_sp[i]))
  
  # load present ensemble
  present <- list.files(path = paste0(run_name, study_sp[i], "/present/ensemble"),
                        pattern = "consensus",
                        full.names = TRUE) %>%
    raster() %>%
    #crop(m) %>%
    mask(m)
  
  plot(present, main = paste(study_sp[i], "present"))
  plot(m, add=TRUE)
  
  for(g in 1:length(ssp_names)){
    
    # load future ensemble for this scenario
    future <- list.files(path = paste0(run_name, study_sp[i], "/future_consensus"),
                         pattern = paste0(ssp_names[g], "_consensus"),
                         full.names = TRUE) %>%
      raster() 
    
    plot(future, main = paste(study_sp[i], "future"))
    diff = future - present
    plot(diff, main = paste(study_sp[i], "diff", ssp_names[g]))
    
    # save diff output
    if (!file.exists(paste0(run_name, study_sp[i], "/diffs"))) {
      dir.create(paste0(run_name, study_sp[i], "/diffs"))
    }
    
    writeRaster(diff, filename = paste0(run_name, study_sp[i], "/diffs/", "diff_", 
                                        ssp_names[g], study_sp[i], ".tif"), 
                overwrite = TRUE)
  }
  
}


# crop and calculate area -------------------------------------------------


for(i in 1:length(study_sp)){
  
  ## load model outputs
  
  present <- list.files(path = paste0(run_name, study_sp[i], "/present/ensemble"),
                        pattern = "consensus",
                        full.names = TRUE) %>%
    raster() %>%
    mask(br_mask)
  
  plot(present)
}




# save output tables ------------------------------------------------------


