# Preparing predictor layers
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(raster)
library(dplyr)
library(readr)
library(caret)


# define study extent -----------------------------------------------------

# based on total distribution of occs + 2º

points <- read_csv("./data/01_occ_full.csv")
coords <- points[ ,2:3]
coordinates(coords) <- c("lon", "lat")
proj4string(coords) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

extent(coords)

# final extent of study area

ext <- extent(-83, -34, -34, 13) 

saveRDS(ext, file = "./data/02_study_extent.rds")


# load predictors ---------------------------------------------------------

layer_path <- c("C:/layers/wc2_current_SA_2.5/")

wc <- list.files(layer_path, pattern = "_SA.asc", full.names = TRUE) %>%
  stack() %>%
  # exclude biovars that combine temp and prec (bio8, bio9, bio18, bio19)
  dropLayer(c("bio8_SA", "bio9_SA", "bio18_SA", "bio19_SA"))


# select uncorrelated variables -------------------------------------------

# load historical records
points <- read_csv("./data/01_occ_hist.csv") %>%
  dplyr::select(!species)

# extract predictor values in each record
vals <- raster::extract(wc, points)

# check which variables are correlated
exclude_vars <- caret::findCorrelation(cor(vals, method = 'spearman'), 
                                       cutoff = 0.8, names = TRUE)

# selecting variables with lower correlation
wc_sel <- wc[[which(!names(wc) %in% exclude_vars)]]

# save txt with selected variable names
write_lines(names(wc_sel), file = "./data/02_selected_variable_names.txt")
