# External validation & model ensembles
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)

# load data ---------------------------------------------------------------


## load validation datasets
occ_hist_valid_100 <- read_csv("./data/01_occ_hist_100km_valid.csv")
coordinates(occ_hist_valid_100) <- c('lon', 'lat')


## load final models

final_models <- list.files(path = "./outputs/models/L_longipalpis/present/final_models",
                           pattern = "raw_mean_th.tif",
                           full.names = TRUE) %>%
                stack()

names(final_models) <- c('bioclim', 'brt', 'glm', 'maxent', 'rf', 'svmk')



# calculate omission rates ------------------------------------------------

vals <- extract(final_models, occ_hist_valid_100[,2:3], df = TRUE)


t <- apply(vals[,2:7], 2, table) %>%
  t() %>%
  data.frame() %>%
  tibble() %>%
  mutate(omission = X0/(X0+X1))



# model ensemble ----------------------------------------------------------



# save outputs ------------------------------------------------------------


