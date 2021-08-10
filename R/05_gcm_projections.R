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
study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")

# name the folder where previous final models were saved
run_name <- c("./outputs/models_hist_100km/")

# name the folder with gcm layers
gcm_folder <- c("C:/layers/raster/worldclim2_cmip6_2.5/")

# selected variables and extent
wc_sel_names <- read_lines("./outputs/02_selected_variable_names.txt")
ext <- readRDS("./outputs/02_study_extent.rds")

# scenario and gcm names
gcm_names <- c("BCC-CSM2-MR", "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1",
               "IPSL-CM6A-LR", "MIROC6", "MIROC-ES2L", "MRI-ESM2-0")

ssp_names <- c("ssp126_2041-2060", "ssp245_2041-2060", "ssp585_2041-2060")


# project selected models -------------------------------------------------

start <- Sys.time()
for(g in 1:length(ssp_names)){
  
  for(h in 1:length(gcm_names)){
    
    pred_gcm <- stack(paste0(gcm_folder, "wc2.1_2.5m_bioc_", gcm_names[h], "_", ssp_names[g], ".tif"))
    
    names(pred_gcm) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_3",
                         "wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6",
                         "wc2.1_2.5m_bio_7", "wc2.1_2.5m_bio_8", "wc2.1_2.5m_bio_9",
                         "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", "wc2.1_2.5m_bio_12",
                         "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_14", "wc2.1_2.5m_bio_15",
                         "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", "wc2.1_2.5m_bio_18",
                         "wc2.1_2.5m_bio_19")
    
    pred_gcm <- subset(pred_gcm, wc_sel_names) %>%
      crop(ext)
    
    for(i in 1:length(study_sp)){
      
      # selected algorithms for this species
      selected_algo <- read_csv(paste0(run_name, study_sp[i], 
                                       "/present/final_models/", study_sp[i], 
                                       "_external_validation.csv")) %>%
        filter(pass == 1) %>%
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
