# Running Ecological Niche Models
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(raster)
library(rgeos)
library(dplyr)
library(readr)
library(modleR)


# loading species occurrence
occ_hist <- read_csv("./data/01_occ_hist.csv")
occ_hist_100 <- read_csv("./data/01_occ_full_100km.csv")

study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", "L_umbratilis",
              "L_flaviscutellata", "L_cruzi", "L_intermedia", "L_complexa", "L_wellcomei")

# set to be modelled
occs <- occ_hist_100

# loading selected predictors
layer_path <- c("D:/OneDrive/work/layers/raster/WorldClim_v2/current_SA_WorldClim2_2.5")

wc_sel_names <- read_lines("./data/02_selected_variable_names.txt")

wc <- list.files(layer_path, pattern = "_SA.asc", full.names = TRUE) %>%
  stack() %>%
  subset(wc_sel_names)

# projections
crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
crs.albers <- CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs") # projected, South America Albers Equal Area Conic


# modleR 1/4: Setup data --------------------------------------------------

for(i in 1:length(study_sp)){
  
  # getting only occurrences for this species
  species_df <- occs[occs$species == study_sp[i], ]
  
  # creating model calibration area for this species
  coords <- species_df[ ,2:3]
  coordinates(coords) <- c("lon", "lat")
  proj4string(coords) <- crs.wgs84  # define original projection - wgs84
  coords <- spTransform(coords, crs.albers)  # project to Albers Equal Area
  mcp <- gConvexHull(coords) # create minimum convex polygon
  mcp_buf <- mcp %>%
    gBuffer(width = gArea(mcp)*2e-07) %>% # draw a buffer of 20% in km around of the minimum convex polygon (Barve et al. 2011)
    spTransform(crs.wgs84) %>% #convert back into wgs84
    SpatialPolygonsDataFrame(data = data.frame("val" = 1, row.names = "buffer"))
  
  # choosing the type of partition depending on the number of records
  partition_type <- ifelse(nrow(species_df) > 50, "crossvalidation", "bootstrap")
  
  #running setup_sdmdata
  setup_sdmdata(species_name = study_sp[i],
                occurrences = as.data.frame(species_df),
                predictors = wc, # set of predictors for running the models
                models_dir = "./outputs/models/", # folder to save partitions
                seed = 123, # set seed for random generation of pseudoabsences
                buffer_type = "user", # buffer type for sampling pseudoabsences
                buffer_shape = mcp_buf, # buffer created in previous step
                clean_dupl = TRUE, # remove duplicate records (done in script 01)
                clean_nas = TRUE, # remove records with na values from variables
                clean_uni = TRUE, # remove records falling at the same pixel (done in script 01)
                png_sdmdata = TRUE, # save minimaps in png
                n_back = nrow(species_df) * 10, # number of pseudoabsences
                partition_type = partition_type,
                cv_partitions = 10, # number of folds for crossvalidation
                cv_n = 1,# number of crossvalidation runs
                boot_n = 10, # number of crossvalidation runs
                boot_proportion = 0.1) # number of partitions in the crossvalidation
}




# modleR 2/4: model calibration -------------------------------------------

for(i in 1:length(study_sp)){
  
  # run selected algorithms for each partition
  do_many(species_name = study_sp[i],
          predictors = wc,
          models_dir = "./outputs/models",
          project_model = F, # project models into other sets of variables
          #proj_data_folder = "./data/env_sel/future/", # folder with projection variables
          #mask = ma_mask, # mask for projecting the models
          png_partitions = TRUE, # save minimaps in png?
          write_bin_cut = TRUE, # save binary and cut outputs?
          dismo_threshold = "spec_sens", # threshold rule for binary outputs
          equalize = TRUE, # equalize numbers of presence and pseudoabsences for random forest and brt
          bioclim = TRUE,
          glm = TRUE,
          maxent = TRUE,
          rf = TRUE,
          svmk = TRUE,
          brt = TRUE,
          svme = TRUE)
}


# modleR 3/4: final models by algo ----------------------------------------

# projections path names
paths <- c("present", "ac85bi50", "he85bi50", "mp85bi50")

for(i in 1:length(study_sp)){
  
  #combine partitions into one final model per algorithm
  final_model(species_name = study_sp[i],
              models_dir = "./outputs/models/",
              which_models = c("raw_mean", "bin_consensus"),
              proj_dir = "present", #"present" "ac85bi50" "he85bi50" "mp85bi50" change paths later
              consensus_level = 0.5, # proportion of models in the binary consensus
              png_final = TRUE,
              overwrite = TRUE)
}

# modleR 4/4: model ensemble ----------------------------------------------

# select algorithms to do the ensemble models
#algo <- c("glm", "maxent", "rf", "svme", "brt")

for(i in 1:length(study_sp)){
  
  # getting only occurrences for this species
  species_df <- occs[occs$species == study_sp[i], ]
  
  #generate ensemble models, combining final models from all algorithms
  ensemble_model(species_name = study_sp[i],
                  occurrences = species_df,
                  #algorithms = algo,
                  models_dir = "./outputs/models/",
                  performance_metric = "TSSmax",
                  proj_dir = "present",
                  which_ensemble = c("weighted_average", "consensus"),
                  which_final = c("raw_mean", "bin_consensus"),
                  #ensemble_dir = "ensemble_five_algo",
                  ensemble_dir = "ensemble",
                  consensus_level = 0.5,
                  png_ensemble = TRUE,
                  uncertainty = TRUE,
                  overwrite = TRUE)
}

