# Projecting selected models into different GCMs and scenarios
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)
library(gbm)
library(randomForest)

# load data and set values ------------------------------------------------

# species names
study_sp <- c("L_complexa", "L_cruzi", "L_flaviscutellata", "L_intermedia",
              "L_longipalpis", "L_migonei", "L_neivai", "L_umbratilis",
              "L_wellcomei", "L_whitmani")

# name the folder where previous final models were saved
run_name <- c("./outputs/models_chelsa_100km/")

# name the folder with gcm layers
gcm_folder <- c("C:/layers/raster/chelsa/2071-2100/")

# selected variables, extent, and wc sample for resolution resampling
wc_sel_names <- read_lines("./outputs/02_chelsa_selected_variable_names.txt")
ext <- readRDS("./outputs/02_study_extent.rds")
wc_sample <- raster("C:/layers/raster/worldclim2_historical_SA_2.5/bio1_SA.tif") %>%
  crop(ext)

# scenario and gcm names
gcm_names <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")

ssp_names <- c("ssp126", "ssp370", "ssp585")

# project selected models -------------------------------------------------

start <- Sys.time()
for(g in 1:length(ssp_names)){
  
  for(h in 1:length(gcm_names)){
    
    # prepare projection layers
    message(paste("Preparing projection layers of", gcm_names[h], 
                  ssp_names[g], sep = " "))
    # load files
    pred_gcm <- list.files(path = paste0(gcm_folder, gcm_names[h], "/", ssp_names[g], "/bio"),
                           full.names = TRUE) %>%
      stack()
    
    # rename for variable selection (according to script 02)
    names(pred_gcm) <- c("CHELSA_bio1_1981.2010_V.2.1", "CHELSA_bio10_1981.2010_V.2.1",
                         "CHELSA_bio11_1981.2010_V.2.1", "CHELSA_bio12_1981.2010_V.2.1",
                         "CHELSA_bio13_1981.2010_V.2.1", "CHELSA_bio14_1981.2010_V.2.1",
                         "CHELSA_bio15_1981.2010_V.2.1", "CHELSA_bio16_1981.2010_V.2.1",
                         "CHELSA_bio17_1981.2010_V.2.1", "CHELSA_bio18_1981.2010_V.2.1",
                         "CHELSA_bio19_1981.2010_V.2.1", "CHELSA_bio2_1981.2010_V.2.1",
                         "CHELSA_bio3_1981.2010_V.2.1", "CHELSA_bio4_1981.2010_V.2.1",
                         "CHELSA_bio5_1981.2010_V.2.1", "CHELSA_bio6_1981.2010_V.2.1",
                         "CHELSA_bio7_1981.2010_V.2.1", "CHELSA_bio8_1981.2010_V.2.1",
                         "CHELSA_bio9_1981.2010_V.2.1", "CHELSA_npp_1981.2010_V.2.1")
    
    # subset, crop and resample (takes a while)
    pred_gcm <- subset(pred_gcm, wc_sel_names) %>%
      crop(ext) %>%
      resample(wc_sample, method = "bilinear") %>%
      mask(wc_sample)
    
    for(i in 1:length(study_sp)){
      
      # selected algorithms for this species
      selected_algo <- read_csv(paste0(run_name, study_sp[i], 
                                       "/present/final_models/", study_sp[i], 
                                       "_external_validation.csv")) %>%
        filter(sens_pass == 1) %>%
        pull(algo)
      
      
      for(j in 1:length(selected_algo)){
        
        # load model objects from selected algorithm
        proj_partitions <- list.files(path = paste0(run_name, study_sp[i], "/present/partitions/"), 
                                      pattern = paste0(selected_algo[j], "_model"), 
                                      full.names = TRUE)
        
        partitions <- raster::stack()
        
        message(paste("Projecting partitions of", study_sp[i], gcm_names[h], 
                      ssp_names[g], selected_algo[j], sep = " "))
        
        for(k in 1:length(proj_partitions)){
          
          load(proj_partitions[k])
          
          if (selected_algo[j] == "brt") {
            mod_proj_cont <- dismo::predict(pred_gcm, mod, n.trees = mod$n.trees, progress = 'text')
          }
          if (selected_algo[j] == "glm") {
            mod_proj_cont <- raster::predict(pred_gcm, mod, type = "response", progress = 'text')
          }
          if (selected_algo[j] %in% c("bioclim", "domain", "maxent", "mahal")) {
            mod_proj_cont <- dismo::predict(pred_gcm, mod, progress = 'text')
          }
          if (selected_algo[j] %in% c("svmk", "maxnet", "svme", "rf")) {
            mod_proj_cont <- raster::predict(pred_gcm, mod, progress = 'text')
          }
          
          # saving partitions
          partitions <- addLayer(partitions, mod_proj_cont)
          
        }
        
        final_raw_mean <- mean(partitions)
        
        selected_thres <- read_csv(paste0(run_name, study_sp[i], 
                                         "/present/final_models/", study_sp[i], 
                                         "_mean_statistics.csv"))
        
        selected_thres <- pull(selected_thres[selected_thres$algorithm == selected_algo[j], "spec_sens"])
        
        final_raw_mean_th <- final_raw_mean > selected_thres
        
        projection_folder <- paste0(run_name, study_sp[i], "/", gcm_names[h], "/", ssp_names[g], "/")
        
        if (file.exists(projection_folder) == FALSE){
          dir.create(paste0(projection_folder), recursive = TRUE, showWarnings = FALSE)
        } 
        
        writeRaster(final_raw_mean, filename = paste0(projection_folder, study_sp[i], "_", 
                                                      selected_algo[j], "_raw_mean.tif"), overwrite = TRUE)
        writeRaster(final_raw_mean_th, filename = paste0(projection_folder, study_sp[i], "_", 
                                                      selected_algo[j], "_raw_mean_th.tif"), overwrite = TRUE)
        
        
      }
      
    }
    
  }
  
}
