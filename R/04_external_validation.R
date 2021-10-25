# External validation of final models
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)


# load data and set values ------------------------------------------------

# load validation datasets
#valid <- read_csv("./data/01_validation_dataset.csv")
valid <- read_csv("./data/01_validation_std_presences.csv")

# species names
study_sp <- c("L_complexa", "L_cruzi", "L_flaviscutellata", "L_intermedia",
              "L_longipalpis", "L_migonei", "L_neivai", "L_umbratilis",
              "L_wellcomei", "L_whitmani")

# set thresholds for keeping algos in the ensemble
sens_thres <- 0.8
tss_thres <- 0.5

# set consensus level for ensemble binary
consensus_level <- 0.5 # majority rule

# name the folder where previous final models were saved
run_name <- c("./outputs/models_chelsa_100km/")


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
  vals <- raster::extract(final_models, valid_sp, df = TRUE) %>%
    dplyr::select(-"ID") %>%
    mutate(pa = valid_sp$pa)
  
  
  # calculate sensitivity ##############################
 
  t <- data.frame(pres_predicted_as_pres = rep(NA, length(algo)),
                  pres_predicted_as_abs = rep(NA, length(algo)),
                  abs_predicted_as_pres = rep(NA, length(algo)),
                  abs_predicted_as_abs = rep(NA, length(algo)))
  for(j in 1:length(algo)){
    t[j, 1] <- count(filter(vals, pa == 1, vals[j] == 1))
    t[j, 2] <- count(filter(vals, pa == 1, vals[j] == 0))
    t[j, 3] <- count(filter(vals, pa == 0, vals[j] == 1))
    t[j, 4] <- count(filter(vals, pa == 0, vals[j] == 0))
  }
  
  t <- as_tibble(t) %>%
    mutate(species = rep(study_sp[i], length(algo)),
           algo = algo,
           sensitivity = pres_predicted_as_pres/(pres_predicted_as_abs+pres_predicted_as_pres),
           specificity = abs_predicted_as_abs/(abs_predicted_as_abs+abs_predicted_as_pres),
           TSS = sensitivity + specificity - 1,
           tss_pass = ifelse(TSS > tss_thres, 1, 0),
           sens_pass = ifelse(sensitivity > sens_thres, 1, 0)) %>%
    relocate(algo, .before = pres_predicted_as_pres) %>%
    relocate(species, .before = algo)
  
  write_csv(t, file = paste0(run_name, study_sp[i], 
                             "/present/final_models/", study_sp[i], 
                             "_external_validation.csv"))
}

summary_valid <- list()
for(i in 1:length(study_sp)){
  summary_valid[[i]] <- read_csv(paste0(run_name, study_sp[i], 
                                     "/present/final_models/", study_sp[i],
                                     "_external_validation.csv"))
}

summary_valid <- data.table::rbindlist(summary_valid) %>%
  mutate(total_valid_records = pres_predicted_as_pres + pres_predicted_as_abs) %>%
  relocate(total_valid_records, .before = sensitivity)

write_csv(summary_valid, file = paste0(run_name, "validation_summary.csv"))
