# Running Ecological Niche Models
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(raster)
library(rgeos)
library(rgdal)
library(dplyr)
library(readr)
library(modleR)


# define settings of this run ---------------------------------------------

## species occurrence set
occ <- read_csv("./data/01_occ_hist_100km.csv")

# paths to layers
layer_path <- c("C:/layers/raster/worldclim2_historical/")

# name a folder to save outputs of this run
run_name <- c("./outputs/models_hist_100km/")


# load data ---------------------------------------------------------------

# selected predictors and extent

wc_sel_names <- read_lines("./outputs/02_selected_variable_names.txt")
ext <- readRDS("./outputs/02_study_extent.rds")

#test
#layer_path <- c("C:/layers/wc2_test_prec_SA_10/")
#wc <- list.files(layer_path, pattern = ".tif", full.names = TRUE) %>%
#  stack() %>%
#  crop(ext)

wc <- list.files(layer_path, pattern = "tif", full.names = TRUE) %>%
  stack() %>%
  subset(wc_sel_names) %>%
  crop(ext)

# projections
## geographical, WGS84
crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  

## projected, South America Albers Equal Area Conic
crs.albers <- CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60
                  +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs") 

study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")


# modleR 1/3: Setup data --------------------------------------------------

for(i in 1:length(study_sp)){
  
  # getting only occurrences for this species
  species_df <- occ[occ$species == study_sp[i], ]
  
  # creating model calibration area for this species
  coords <- species_df[ ,2:3]
  coordinates(coords) <- c("lon", "lat")
  proj4string(coords) <- crs.wgs84  # define original projection - wgs84
  coords <- spTransform(coords, crs.albers)  # project to Albers Equal Area
  mcp <- gConvexHull(coords) # create convex hull of occurrence records
  mcp_buf <- mcp %>%
    gBuffer(width = gArea(mcp)*2e-07) %>% # draw a buffer of +20% of the convex hull (Barve et al. 2011)
    spTransform(crs.wgs84) %>% #convert back into wgs84
    SpatialPolygonsDataFrame(data = data.frame("val" = 1, row.names = "buffer"))
  
  # modleR will automatically do jackknife if number of records is <= 10
  
  #running setup_sdmdata
  setup_sdmdata(species_name = study_sp[i],
                occurrences = as.data.frame(species_df),
                predictors = wc, # set of predictors for running the models
                models_dir = run_name, # folder to save partitions
                seed = 123, # set seed for random generation of pseudoabsences
                buffer_type = "user", # buffer type for sampling pseudoabsences
                buffer_shape = mcp_buf, # buffer created in previous step
                select_variables = FALSE, # do not select/exclude variables
                clean_nas = TRUE, # remove records with na values from variables
                clean_uni = TRUE, # remove records falling at the same pixel
                png_sdmdata = TRUE, # save minimaps in png
                n_back = nrow(species_df) * 10, # number of pseudoabsences
                partition_type = "crossvalidation", 
                cv_partitions = 10, # number of folds for crossvalidation
                cv_n = 1)# number of crossvalidation runs
}



# modleR 2/3: model calibration -------------------------------------------


for(i in 1:length(study_sp)){
  
  # run selected algorithms for each partition
  do_many(species_name = study_sp[i],
          predictors = wc,
          models_dir = run_name,
          project_model = FALSE, # gcm projections will be done later
          write_rda = TRUE, # will need this for the gcm projections
          png_partitions = TRUE, # save minimaps in png?
          write_bin_cut = FALSE, # save binary and cut outputs?
          dismo_threshold = "spec_sens", # threshold rule for binary outputs
          equalize = TRUE, # equalize numbers of pres and pseudoabs for random forest
          bioclim = TRUE,
          glm = TRUE,
          maxent = TRUE,
          rf = TRUE,
          svmk = TRUE)
}


# L. cruzi, L. complexa and L. wellcomei did not have enough unique records to run BRT,
# so using the following species subset to run BRT for them:

study_sp_brt <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
                  "L_umbratilis", "L_flaviscutellata", "L_intermedia")


for(i in 1:length(study_sp_brt)){
  
  # run selected algorithms for each partition
  do_many(species_name = study_sp_brt[i],
          predictors = wc,
          models_dir = run_name,
          project_model = FALSE, # gcm projections will be done later
          write_rda = TRUE, # will need this for the gcm projections
          png_partitions = TRUE, # save minimaps in png?
          write_bin_cut = FALSE, # save binary and cut outputs?
          dismo_threshold = "spec_sens", # threshold rule for binary outputs
          equalize = FALSE, # do not equalize numbers of pres and pseudoabs for brt
          brt = TRUE)
}


# modleR 3/3: final models by algo ----------------------------------------


for(i in 1:length(study_sp)){
  
  #combine partitions into one final model per algorithm
  final_model(species_name = study_sp[i],
              models_dir = run_name,
              scale_models = TRUE, # convert model outputs to 0-1
              which_models = c("raw_mean", "raw_mean_th"), # mean outputs by algo and 
              # binarise them using the mean of maxTSS thresholds of each partition
              mean_th_par = "spec_sens",
              png_final = TRUE,
              overwrite = TRUE)

}
