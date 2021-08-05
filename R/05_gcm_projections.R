# Projecting selected models to different GCMs and scenarios
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)
library(modleR)

# load data and set values ------------------------------------------------

# species names
study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")

# name the folder where previous final models were saved
run_name <- c("./outputs/models_hist_100km/")



# project selected models -------------------------------------------------



# get validation results for this species
valid <- read_csv(paste0(run_name, study_sp[i], 
                         "/present/final_models/", study_sp[i], 
                         "_external_validation.csv"))







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

# testing project later

load("./outputs/models_hist_100km/L_neivai/present/partitions/bioclim_model_L_neivai_1_5.rda")

pred_proj <- raster::stack(list.files("C:/layers/wc2_cmip6_SA_2.5/BCC.CSM2.MR_ssp245_2021.2040/", full.names = TRUE))
mod_proj_cont <- dismo::predict(pred_proj, mod)
mod_proj_bin <- mod_proj_cont > th_mod

proj_data_folder = "future_layers"
pfiles <- list.dirs(proj_data_folder, recursive = FALSE)
for (proje in pfiles) {
  v <- strsplit(proje, "/")
  name_proj <- v[[1]][length(v[[1]])]
  projection_folder <- paste0(models_dir, "/", 
                              species_name, "/", name_proj, "/partitions")
  if (file.exists(projection_folder) == FALSE) 
    dir.create(paste0(projection_folder), recursive = TRUE, 
               showWarnings = FALSE)
  pred_proj <- raster::stack(list.files(proje, 
                                        full.names = TRUE))
  pred_proj <- raster::subset(pred_proj, names(predictors))
  message(name_proj)
  message("projecting models")
  if (algorithm == "brt") {
    mod_proj_cont <- dismo::predict(pred_proj, 
                                    mod, n.trees = n.trees)
  }
  if (algorithm == "glm") {
    mod_proj_cont <- raster::predict(pred_proj, 
                                     mod, type = "response")
  }
  if (algorithm %in% c("bioclim", "domain", 
                       "maxent", "mahal")) {
    mod_proj_cont <- dismo::predict(pred_proj, 
                                    mod)
  }
  if (algorithm %in% c("svmk", "maxnet", 
                       "svme", "rf")) {
    mod_proj_cont <- raster::predict(pred_proj, 
                                     mod)
  }
  if (write_bin_cut == TRUE) {
    mod_proj_bin <- mod_proj_cont > th_mod
    mod_proj_cut <- mod_proj_bin * mod_proj_cont
  }
  if (class(mask) == "SpatialPolygonsDataFrame") {
    mod_proj_cont <- crop_model(mod_proj_cont, 
                                mask)
    mod_proj_bin <- crop_model(mod_proj_bin, 
                               mask)
    mod_proj_cut <- crop_model(mod_proj_cut, 
                               mask)
  }
  message("writing projected models raster")
  raster::writeRaster(x = mod_proj_cont, filename = paste0(projection_folder, 
                                                           "/", algorithm, "_cont_", species_name, 
                                                           "_", i, "_", g, ".tif"), 
                      overwrite = TRUE)
  if (write_bin_cut == TRUE) {
    raster::writeRaster(x = mod_proj_bin, filename = paste0(projection_folder, 
                                                            "/", algorithm, "_bin_", 
                                                            species_name, "_", i, "_", 
                                                            g, ".tif"), overwrite = TRUE)
    raster::writeRaster(x = mod_proj_cut, filename = paste0(projection_folder, 
                                                            "/", algorithm, "_cut_", 
                                                            species_name, "_", i, "_", 
                                                            g, ".tif"), overwrite = TRUE)
  }
  
  
  