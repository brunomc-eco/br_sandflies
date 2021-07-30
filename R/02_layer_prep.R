# Preparing predictor layers
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(raster)
library(dplyr)
library(readr)


layer_path <- c("D:/OneDrive/work/layers/raster/WorldClim_v2/current_SA_WorldClim2_2.5")

wc <- list.files(layer_path, pattern = "_SA.asc", full.names = TRUE) %>%
  stack() %>%
  # exclude biovars that combine temp and prec (bio8, bio9, bio18, bio19)
  dropLayer(c("bio8_SA", "bio9_SA", "bio18_SA", "bio19_SA"))


# select uncorrelated variables -------------------------------------------

# load historical records
points <- read_csv("./data/01_occ_hist.csv") %>%
  select(!species)

# extract predictor values in each record
vals <- extract(wc, points)

# check which variables are correlated
exclude_vars <- caret::findCorrelation(cor(vals2, method = 'spearman'), cutoff = 0.8, names = TRUE)

# selecting variables with lower correlation
wc_sel <- wc[[which(!names(wc) %in% exclude_vars)]]

# save txt with selected variable names
write_lines(names(wc_sel), file = "./data/02_selected_variable_names.txt")
