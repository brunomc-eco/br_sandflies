# Model ensemble based on sensitivity
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)


# load data and set values ------------------------------------------------

# species names
study_sp <- c("L_complexa", "L_cruzi", "L_flaviscutellata", "L_intermedia",
              "L_longipalpis", "L_migonei", "L_neivai", "L_umbratilis",
              "L_wellcomei", "L_whitmani")

# name the folder where previous final models were saved
run_name <- c("./outputs/models_chelsa_100km/")

# scenario and gcm names
gcm_names <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")

ssp_names <- c("ssp126", "ssp370", "ssp585")

# set consensus level for ensemble binary
consensus_level <- 0.5 # majority rule

# load validation dataset
valid <- read_csv("./data/01_validation_std_presences.csv")

# set sensitivity threshold for keeping algos in the ensemble
sens_thres <- 0.8


# historical models -------------------------------------------------------

valid_ensembles <- list()

for(i in 1:length(study_sp)){
  
  # selected algorithms for this species
  selected_algo <- read_csv(paste0(run_name, study_sp[i], 
                                   "/present/final_models/", study_sp[i], 
                                   "_external_validation.csv")) %>%
    filter(sens_pass == 1) %>%
    pull(algo)
  
  selected_algo <- as.character(selected_algo)
  
  # weights for ensemble
  w <- read_csv(paste0(run_name, study_sp[i], 
                       "/present/final_models/", study_sp[i], 
                       "_external_validation.csv")) %>%
    filter(sens_pass == 1) %>%
    pull(sensitivity)
  
  # load final models for this species, only selected algo
  raw_mean_models <- list.files(path = paste0(run_name, study_sp[i],
                                              "/present/final_models"),
                                pattern = "raw_mean.tif",
                                full.names = TRUE) %>%
    stack()
  
  raw_mean_th_models <- list.files(path = paste0(run_name, study_sp[i],
                                                 "/present/final_models"),
                                   pattern = "raw_mean_th.tif",
                                   full.names = TRUE) %>%
    stack()
  
  algo <- read_csv(paste0(run_name, study_sp[i], 
                          "/present/final_models/", study_sp[i], 
                          "_mean_statistics.csv")) %>%
    pull(algorithm)
  
  names(raw_mean_models) <- algo
  names(raw_mean_th_models) <- algo
  
  raw_mean_models <- raster::subset(raw_mean_models, selected_algo)
  raw_mean_th_models <- raster::subset(raw_mean_th_models, selected_algo)
  
  # create ensembles
  ens <- raster::weighted.mean(raw_mean_models, w)
  ens_bin <- mean(raw_mean_th_models) > consensus_level
  
  # evaluate ensembles
  # get validation records for this species
  valid_sp <- filter(valid, species == study_sp[i])
  coordinates(valid_sp) <- c('lon', 'lat')
  
  # extract predicted values
  vals <- extract(ens_bin, valid_sp, df = TRUE) %>%
    dplyr::select(-"ID")
  
  # calculate sensitivity
  
  t <- data.frame(pres_predicted_as_pres = NA,
                  pres_predicted_as_abs = NA)
  t[1,1] <- count(filter(vals, vals[1] == 1))
  t[1,2] <- count(filter(vals, vals[1] == 0))
  
  valid_ensembles[[i]] <- as_tibble(t) %>%
    mutate(total_valid_records = pres_predicted_as_pres + pres_predicted_as_abs,
           species = study_sp[i],
           algo = "ensemble_consensus",
           sensitivity = pres_predicted_as_pres/total_valid_records,
           sens_pass = ifelse(sensitivity >= sens_thres, 1, 0)) %>%
    relocate(species, .before = pres_predicted_as_pres) %>%
    relocate(algo, .before = pres_predicted_as_pres) %>%
    relocate(total_valid_records, .before = pres_predicted_as_pres)
  
  
  # save output maps and table
  ensemble_folder <- paste0(run_name, study_sp[i], "/present/ensemble/")
  if (!file.exists(ensemble_folder)) {
    dir.create(ensemble_folder)
  }
  
  writeRaster(ens, filename = paste0(ensemble_folder, study_sp[i], "_weighted_mean_ensemble.tif"),
              overwrite = TRUE)
  writeRaster(ens_bin, filename = paste0(ensemble_folder, study_sp[i], "_consensus_ensemble.tif"),
              overwrite = TRUE)
  
  write_csv(t, file = paste0(ensemble_folder, study_sp[i], "_ensemble_validation.csv"))
 
}

valid_ensembles <- rbindlist(valid_ensembles)

write_csv(valid_ensembles, file = paste0(run_name, "06_validation_ensembles.csv"))



# climate change projections ----------------------------------------------

for(i in 1:length(study_sp)){
  
  # selected algorithms for this species
  selected_algo <- read_csv(paste0(run_name, study_sp[i], 
                                   "/present/final_models/", study_sp[i], 
                                   "_external_validation.csv")) %>%
    filter(sens_pass == 1) %>%
    pull(algo)
  
  selected_algo <- as.character(selected_algo)
  
  # weights for ensemble
  w <- read_csv(paste0(run_name, study_sp[i], 
                       "/present/final_models/", study_sp[i], 
                       "_external_validation.csv")) %>%
    filter(sens_pass == 1) %>%
    pull(sensitivity)
  
  # load final models for this species, only selected algo
  
  for(h in 1:length(gcm_names)){
    
    
    for(g in 1:length(ssp_names)){
      
      projection_folder <- paste0(run_name, study_sp[i], "/", gcm_names[h], "/", ssp_names[g], "/")
      
      raw_mean_models <- list.files(path = projection_folder, pattern = "raw_mean.tif", full.names = TRUE) %>%
        stack()
      names(raw_mean_models) <- selected_algo
      
      raw_mean_th_models <- list.files(path = projection_folder, pattern = "raw_mean_th.tif", full.names = TRUE) %>%
        stack()
      names(raw_mean_th_models) <- selected_algo
      
      message(paste("Calculating ensembles of", study_sp[i], gcm_names[h], 
                    ssp_names[g], sep = " "))
      
      # create ensembles
      ens <- raster::weighted.mean(raw_mean_models, w)
      ens_bin <- mean(raw_mean_th_models) > consensus_level
      
      # save output maps
      ensemble_folder <- paste0(projection_folder, "ensemble/")
      
      if (!file.exists(ensemble_folder)) {
        dir.create(ensemble_folder)
      }
      
      writeRaster(ens, filename = paste0(ensemble_folder, study_sp[i], "_weighted_mean_ensemble.tif"),
                  overwrite = TRUE)
      writeRaster(ens_bin, filename = paste0(ensemble_folder, study_sp[i], "_consensus_ensemble.tif"),
                  overwrite = TRUE)
      
    }
    
  }
    
}



# multi-gcm consensus -----------------------------------------------------


# rescale_layer function from modleR:

rescale_layer <- function(layers) {
  message("Standardizing models from 0 to 1")
  if (missing(layers)) {
    stop("No layers were provided. Please enter a Raster or a RasterStack")
  }
  for (i in 1:dim(layers)[3]) {
    stand <- function(x) {
      (x - min(layers[[i]][], na.rm = TRUE)) /
        (max(layers[[i]][], na.rm = TRUE) - min(layers[[i]][],
                                                na.rm = TRUE))
    }
    bb <- raster::calc(layers[[i]], stand)
    bb
    if (i == 1) {
      cc <- raster::stack(bb)
      names(cc)[i] <- names(layers)[i]
    }
    else{
      cc <- raster::stack(cc, bb)
      names(cc)[i] <- names(layers)[i]
    }
    if (i == dim(layers)[3]) {
      layers <- cc
      rm(cc)
      return(layers)
    }
  }
}


for(i in 1:length(study_sp)){
  
  for(g in 1:length(ssp_names)){
    
    ssp_consensus <- stack()
    ssp_mean <- stack()
    for(h in 1:length(gcm_names)){
      
      consensus_ensemble <- list.files(path = paste0(run_name, study_sp[i], "/", 
                                                     gcm_names[h], "/", ssp_names[g],
                                                     "/", "ensemble/"),
                                       pattern = "consensus",
                                       full.names = TRUE) %>%
        raster()
      
      ssp_consensus <- addLayer(ssp_consensus, consensus_ensemble)
      
      mean_ensemble <- list.files(path = paste0(run_name, study_sp[i], "/", 
                                                gcm_names[h], "/", ssp_names[g],
                                                "/", "ensemble/"),
                                       pattern = "mean",
                                       full.names = TRUE) %>%
        raster()
      
      ssp_mean <- addLayer(ssp_mean, mean_ensemble)
      
    }
    
    message(paste("Calculating ssp consensus of", study_sp[i], 
                  ssp_names[g], sep = " "))
    
    final_consensus <- calc(ssp_consensus, sum)
    
    message(paste("Calculating ssp means of", study_sp[i], 
                  ssp_names[g], sep = " "))
    
    final_mean <- calc(ssp_mean, mean)
    final_median <- calc(ssp_mean, median)
    final_sd <- calc(ssp_mean, sd)
    
    # saving outputs
    
    future_folder <- paste0(run_name, study_sp[i], "/future_consensus/")
    
    if (!file.exists(future_folder)) {
      dir.create(future_folder)
    }
    
    writeRaster(final_consensus, filename = paste0(future_folder, study_sp[i], "_", ssp_names[g], "_consensus.tif"),
                overwrite = TRUE)
    writeRaster(final_mean, filename = paste0(future_folder, study_sp[i], "_", ssp_names[g], "_mean.tif"),
                overwrite = TRUE)
    writeRaster(final_median, filename = paste0(future_folder, study_sp[i], "_", ssp_names[g], "_median.tif"),
                overwrite = TRUE)
    writeRaster(final_sd, filename = paste0(future_folder, study_sp[i], "_", ssp_names[g], "_sd.tif"),
                overwrite = TRUE)
    
  }
  
}

