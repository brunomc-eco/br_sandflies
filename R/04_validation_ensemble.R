# External validation & model ensembles
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)
library(modleR)


# load data and set values ------------------------------------------------


# load validation dataset
valid <- read_csv("./data/01_occ_hist_100km_valid.csv")

# species names
study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")

# set sensitivity threshold for keeping algos in the ensemble
sens_thres <- 0.8

# set consensus level for ensemble binary
consensus_level <- 0.5

# name the folder where previous final models were saved
run_name <- c("./outputs/models_hist_100km/")

# calculate sensitivity/omission ------------------------------------------

for(i in 1:length(study_sp)){
  
  # load final models for this species
  final_models <- list.files(path = paste0(run_name, study_sp[i], 
                                           "/present/final_models"),
                             pattern = "raw_mean_th.tif",
                             full.names = TRUE) %>%
    stack()
  
  algo <- read_csv(paste0(run_name, study_sp[i], 
                          "/present/final_models/", study_sp[i], 
                          "_mean_statistics.csv")) %>%
    pull(algorithm)
    
  names(final_models) <- algo
  
  
  # get validation records for this species
  valid_sp <- filter(valid, species == study_sp[i])
  coordinates(valid_sp) <- c('lon', 'lat')
  
  # extract predicted values
  vals <- extract(final_models, valid_sp, df = TRUE) %>%
    dplyr::select(-"ID")
  
  # calculate sensitivity
  t <- apply(vals, 2, table) %>%
    t() %>%
    as_tibble() %>%
    rename(zeros = "0",
           ones = "1") %>%
    mutate(algo = algo,
           sensitivity = ones/(zeros+ones),
           pass = ifelse(sensitivity > sens_thres, 1, 0))
  
  write_csv(t, file = paste0(run_name, study_sp[i], 
                             "/present/final_models/", study_sp[i], 
                             "_external_validation.csv"))
}



# ensemble based on sensitivity -------------------------------------------

# projections path names (GCMs)
paths <- c("present", 
           "BCC.CSM2.MR_ssp245_2021.2040", 
           "BCC.CSM2.MR_ssp585_2021.2040", 
           "CanESM5_ssp245_2021.2040", 
           "CanESM5_ssp585_2021.2040", 
           "CNRM.CM6.1_ssp245_2021.2040", 
           "CNRM.CM6.1_ssp585_2021.2040", 
           "CNRM.ESM2.1_ssp245_2021.2040", 
           "CNRM.ESM2.1_ssp585_2021.2040", 
           "MIROC.ES2L_ssp245_2021.2040", 
           "MIROC.ES2L_ssp585_2021.2040", 
           "MIROC6_ssp245_2021.2040", 
           "MIROC6_ssp585_2021.2040")

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
  
  for(j in paths){
    
    # load final models for this species, only selected algo
    raw_mean_models <- list.files(path = paste0(run_name, study_sp[i], j,
                                                "/final_models"),
                                  pattern = "raw_mean.tif",
                                  full.names = TRUE) %>%
      stack()
    
    raw_mean_th_models <- list.files(path = paste0(run_name, study_sp[i], j,
                                                   "/final_models"),
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
    t <- apply(vals, 2, table) %>%
      t() %>%
      as_tibble() %>%
      rename(zeros = "0",
             ones = "1") %>%
      mutate(algo = c("consensus_ensemble"),
             sensitivity = ones/(zeros+ones))
    
    # save output maps and table
    ensemble_folder <- paste0(run_name, study_sp[i], "/", j, 
                              "/ensemble/")
    if (!file.exists(ensemble_folder)) {
      dir.create(ensemble_folder)
    }
    
    writeRaster(ens, filename = paste0(ensemble_folder, study_sp[i], "_weighted_mean_ensemble.tif"))
    writeRaster(ens_bin, filename = paste0(ensemble_folder, study_sp[i], "_consensus_ensemble.tif"))
    
    write_csv(t, file = paste0(ensemble_folder, study_sp[i], "_ensemble_validation.csv"))
    
  }
}
