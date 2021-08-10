# Preparing predictor layers
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(raster)
library(dplyr)
library(readr)
library(caret)
library(kuenm)


# define study extent -----------------------------------------------------

# based on total distribution of occs + 2º

points <- read_csv("./data/01_occ_hist.csv")
coords <- points[ ,2:3]
coordinates(coords) <- c("lon", "lat")
proj4string(coords) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

extent(coords)

# final extent of study area

ext <- extent(-83, -34, -34, 13) 

saveRDS(ext, file = "./outputs/02_study_extent.rds")


# load predictors ---------------------------------------------------------

layer_path <- c("C:/layers/raster/worldclim2_historical/")

wc <- list.files(layer_path, pattern = ".tif", full.names = TRUE) %>%
  stack() %>%
  # exclude biovars that combine temp and prec (bio8, bio9, bio18, bio19)
  dropLayer(c("wc2.1_2.5m_bio_8", "wc2.1_2.5m_bio_9", 
              "wc2.1_2.5m_bio_18", "wc2.1_2.5m_bio_19"))


# select uncorrelated variables -------------------------------------------

# load historical records
points <- read_csv("./data/01_occ_hist.csv") %>%
  dplyr::select(!species)

# extract predictor values in each record
vals <- raster::extract(wc, points)

# check which variables are correlated
exclude_vars <- findCorrelation(cor(vals, method = 'spearman'),
                                cutoff = 0.8, names = TRUE)

# selecting variables with lower correlation
wc_sel <- wc[[which(!names(wc) %in% exclude_vars)]]

wc_sel_names <- names(wc_sel)

# save txt with selected variable names
write_lines(wc_sel_names, file = "./outputs/02_selected_variable_names.txt")


# MOP analysis ------------------------------------------------------------

gcm_folder <- c("C:/layers/raster/worldclim2_cmip6_2.5/")

# scenario and gcm names
gcm_names <- c("BCC-CSM2-MR", "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1",
               "IPSL-CM6A-LR", "MIROC6", "MIROC-ES2L", "MRI-ESM2-0")

ssp_names <- c("ssp126_2041-2060", "ssp245_2041-2060", "ssp585_2041-2060")

wc_sel <- wc_sel %>%
  crop(ext)


for(g in 1:length(ssp_names)){
  
  for(h in 1:length(gcm_names)){
    
    pred_gcm <- stack(paste0(gcm_folder, "wc2.1_2.5m_bioc_", gcm_names[h], "_", ssp_names[g], ".tif"))
    
    names(pred_gcm) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_3",
                         "wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6",
                         "wc2.1_2.5m_bio_7", "wc2.1_2.5m_bio_8", "wc2.1_2.5m_bio_9",
                         "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", "wc2.1_2.5m_bio_12",
                         "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_14", "wc2.1_2.5m_bio_15",
                         "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", "wc2.1_2.5m_bio_18",
                         "wc2.1_2.5m_bio_19")
    
    pred_gcm <- subset(pred_gcm, wc_sel_names) %>%
      crop(ext)
    
    mop_layer <- kuenm::kuenm_mop(M.variables = wc_sel, G.stack = pred_gcm, percent = 5,
                                  parallel = TRUE, n.cores = 3)
    
    #mop_layer <- kuenm_mmop(G.var.dir = "C:/layers/raster/worldclim2_testing_future_SA_10/",
    #                        M.var.dir = "C:/layers/raster/worldclim2_testing_prec_SA_10/",
    #                        is.swd = FALSE,
    #                        sets.var = c("prec_01_SA.tif", "prec_02_SA.tif", "prec_03_SA.tif"),
    #                        percent = 5,
    #                        out.mop = "C:/layers/raster/mop_test/")
    
  }
  
}

