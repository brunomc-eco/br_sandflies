# Model ensemble
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)
library(modleR)


# ensemble based on sensitivity -------------------------------------------


for(i in 1:length(study_sp)){
  
  # selected algorithms for this species
  selected_algo <- read_csv(paste0(run_name, study_sp[i], 
                                   "/present/final_models/", study_sp[i], 
                                   "_external_validation.csv")) %>%
    filter(pass == 1) %>%
    pull(algo)
  
  selected_algo <- as.character(selected_algo)
  
  # weights for ensemble
  w <- read_csv(paste0(run_name, study_sp[i], 
                       "/present/final_models/", study_sp[i], 
                       "_external_validation.csv")) %>%
    filter(pass == 1) %>%
    pull(sensitivity)
  
  #for(j in paths){
  
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
  #ens_avg <- raster::mean(raw_mean_models)
  ens_bin <- mean(raw_mean_th_models) > consensus_level
  
  # evaluate ensembles
  # get validation records for this species
  valid_sp <- filter(valid, species == study_sp[i])
  coordinates(valid_sp) <- c('lon', 'lat')
  
  # extract predicted values
  vals <- extract(ens_bin, valid_sp, df = TRUE) %>%
    dplyr::select(-"ID")
  
  # calculate sensitivity
  t <- tibble(species = study_sp[i],
              predicted_as_abs = count(filter(vals, vals[1] == 0)),
              predicted_as_pres = count(filter(vals, vals[1] == 1)),
              algo = "ensemble_consensus") %>%
    mutate(total_valid = predicted_as_abs + predicted_as_pres,
           sensitivity = predicted_as_pres/total_valid,
           specificity = predicted_as_abs/total_valid)
  
  # save output maps and table
  ensemble_folder <- paste0(run_name, study_sp[i], "/",
                            "/present/ensemble/")
  if (!file.exists(ensemble_folder)) {
    dir.create(ensemble_folder)
  }
  
  writeRaster(ens, filename = paste0(ensemble_folder, study_sp[i], "_weighted_mean_ensemble.tif"),
              overwrite = TRUE)
  writeRaster(ens_bin, filename = paste0(ensemble_folder, study_sp[i], "_consensus_ensemble.tif"),
              overwrite = TRUE)
  
  write_csv(t, file = paste0(ensemble_folder, study_sp[i], "_ensemble_validation.csv"))
  
  #}
}
