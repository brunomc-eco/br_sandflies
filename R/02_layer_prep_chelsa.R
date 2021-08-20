# Preparing predictor layers - CHELSA dataset
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(raster)
library(dplyr)
library(readr)
library(caret)
library(kuenm)



# load data ---------------------------------------------------------------

layer_path <- c("C:/layers/raster/chelsa/1981-2010/bio/")

ext <- readRDS("./outputs/02_study_extent.rds")

wc_sample <- raster("C:/layers/raster/worldclim2_historical_SA_2.5/bio1_SA.tif") %>%
  crop(ext)


# crop & resample ---------------------------------------------------------

chelsa <- list.files(layer_path, pattern = ".tif", full.names = TRUE) %>%
  stack() %>%
  dropLayer(c("CHELSA_bio8_1981.2010_V.2.1", "CHELSA_bio9_1981.2010_V.2.1", 
              "CHELSA_bio18_1981.2010_V.2.1", "CHELSA_bio19_1981.2010_V.2.1"))

chelsa_resampled <- stack()
for(i in 1:nlayers(chelsa)){
  
  a <- crop(chelsa[[i]], wc_sample) %>%
    resample(wc_sample, method = "bilinear") %>%
    mask(wc_sample)
  
  chelsa_resampled <- addLayer(chelsa_resampled, a)
  
}


# select uncorrelated variables -------------------------------------------

# load historical records
points <- read_csv("./data/01_occ_chelsa.csv") %>%
  dplyr::select(!species)

# extract predictor values in each record
vals <- raster::extract(chelsa_resampled, points)

# check which variables are correlated
exclude_vars <- findCorrelation(cor(vals, method = 'spearman'),
                                cutoff = 0.8, names = TRUE)

# selecting variables with lower correlation
chelsa_sel <- chelsa_resampled[[which(!names(chelsa_resampled) %in% exclude_vars)]]

chelsa_sel_names <- names(chelsa_sel)

# save txt with selected variable names
write_lines(chelsa_sel_names, file = "./outputs/02_chelsa_selected_variable_names.txt")
