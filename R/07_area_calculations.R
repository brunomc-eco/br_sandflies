# ENM output area calculations
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)
library(rgdal)



# load data and set values ------------------------------------------------

# paths to layers
layer_path <- c("C:/layers/raster/worldclim2_historical/")

# selected predictors and extent

wc_sel_names <- read_lines("./outputs/02_selected_variable_names.txt")
ext <- readRDS("./outputs/02_study_extent.rds")


# species names
study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")

# name the folder where previous final models were saved
run_name <- c("./outputs/models_hist_100km/")

# scenario and gcm names
gcm_names <- c("BCC-CSM2-MR", "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1",
               "IPSL-CM6A-LR", "MIROC6", "MIROC-ES2L", "MRI-ESM2-0")

ssp_names <- c("ssp126_2041-2060", "ssp245_2041-2060", "ssp585_2041-2060")

# mask for cropping results to Brazil
br_mask <- readOGR("C:/layers/vector/GADM_Continental South America/BRA_adm/BRA_adm0.shp")







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


