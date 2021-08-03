# Running Ecological Niche Models
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(raster)
library(rgeos)
library(rgdal)
library(dplyr)
library(readr)
library(modleR)


# loading occs and predictors ---------------------------------------------

## species occurrence
occ_hist <- read_csv("./data/01_occ_hist.csv")
occ_hist_100 <- read_csv("./data/01_occ_full_100km.csv")

study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")

# set to be modelled
occs <- occ_hist_100


## selected predictors

wc_sel_names <- read_lines("./data/02_selected_variable_names.txt")
ext <- readRDS("./data/02_study_extent.rds")

#test
#layer_path <- c("C:/layers/wc2_test_prec_SA_10/")
#wc <- list.files(layer_path, pattern = ".tif", full.names = TRUE) %>%
#  stack() %>%
#  crop(ext)

layer_path <- c("C:/layers/wc2_current_SA_2.5/")

wc <- list.files(layer_path, pattern = "_SA.asc", full.names = TRUE) %>%
  stack() %>%
  subset(wc_sel_names) %>%
  crop(ext)

# projections
crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, WGS84
crs.albers <- CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60
                  +x_0=0 +y_0=0 +ellps=aust_SA 
                  +units=m +no_defs") # projected, SA Albers Equal Area Conic


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
                select_variables = FALSE, # do not select/exclude variables
                clean_nas = TRUE, # remove records with na values from variables
                clean_uni = TRUE, # remove records falling at the same pixel
                png_sdmdata = TRUE, # save minimaps in png
                n_back = nrow(species_df) * 10, # number of pseudoabsences
                partition_type = partition_type,
                cv_partitions = 10, # number of folds for crossvalidation
                cv_n = 1,# number of crossvalidation runs
                boot_n = 10, # number of bootstrap runs
                boot_proportion = 0.1) # number of partitions in the bootstrap
}



# modleR 2/4: model calibration -------------------------------------------


for(i in 1:length(study_sp)){
  
  # run selected algorithms for each partition
  do_many(species_name = study_sp[i],
          predictors = wc,
          models_dir = "./outputs/models",
          project_model = TRUE, # project models into other sets of variables
          proj_data_folder = "C:/layers/future_test", # folder with GCM projection variables
          png_partitions = TRUE, # save minimaps in png?
          write_bin_cut = FALSE, # save binary and cut outputs?
          dismo_threshold = "spec_sens", # threshold rule for binary outputs
          equalize = TRUE, # equalize numbers of presence and pseudoabsences for random forest and brt
          bioclim = TRUE,
          glm = TRUE,
          maxent = TRUE,
          rf = TRUE,
          svmk = TRUE,
          brt = TRUE)
}


# modleR 3/4: final models by algo ----------------------------------------

# projections path names (GCMs)
paths <- c("present", "futtest1", "futtest2")

for(i in 1:length(study_sp)){
  
  for(path in paths){
    
    #combine partitions into one final model per algorithm
    final_model(species_name = study_sp[i],
                models_dir = "./outputs/models/",
                scale_models = TRUE, # convert model outputs to 0-1
                which_models = c("raw_mean", "raw_mean_th"), # mean outputs by algo and binarise them using the mean of maxTSS thresholds of each partition
                mean_th_par = "spec_sens",
                proj_dir = path,
                #consensus_level = 0.5, # proportion of models in the binary consensus
                png_final = TRUE,
                overwrite = TRUE)
    
  }
  
}
