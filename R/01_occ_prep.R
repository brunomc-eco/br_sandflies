# Preparing occurrence data
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(spThin)
library(raster)
library(rgdal)
library(rgeos)
library(data.table)


# load data ---------------------------------------------------------------

# loading occurrence data files, historical for chelsa dataset (1981-2010)
occ_chelsa <- read_csv("./data/raw/occs_1981-2010.csv")

# select study species
study_sp <- c("L_complexa", "L_cruzi", "L_flaviscutellata", "L_intermedia",
              "L_longipalpis", "L_migonei", "L_neivai", "L_umbratilis",
              "L_wellcomei", "L_whitmani")


# shapefiles with model calibration area by species
calib_files <- list.files("./data/shp/calib", pattern = ".shp", full.names = TRUE)

# WGS84 datum CRS
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")


# functions ---------------------------------------------------------------

# thin records and build species/lat/lon file
thinning <- function(x,y){
  l <- list()
  for(i in 1:length(study_sp)){
    df <- x %>%
      filter(species %in% study_sp[i])
    thinned <- thin(df,
                    lat.col = "lat",
                    long.col = "lon",
                    spec.col = "species",
                    thin.par = y, # distance in km
                    reps = 1,
                    locs.thinned.list.return = TRUE,
                    write.files = FALSE,
                    write.log.file = FALSE)
    l[[i]] <- data.frame(species = rep(study_sp[i], nrow(thinned[[1]])),
                         lon = thinned[[1]]$Longitude,
                         lat = thinned[[1]]$Latitude)
  }
  rbindlist(l)
}

# counting records by species in different sets
count_unique <- function(x){
  x %>%
    group_by(species) %>%
    summarize(count = n())
}

# filter records for model calibration ------------------------------------

# retaining only records from study species, with unique coordinates
occ_chelsa <- occ_chelsa %>%
  filter(lon_dec != "NA",
         lat_dec != "NA",
         considered.species %in% study_sp) %>%
  rename(species = considered.species,
         lon = lon_dec,
         lat = lat_dec) %>%
  dplyr::select(species, lon, lat) %>%
  distinct()

occ_chelsa_10km <- thinning(occ_chelsa, 10)
occ_chelsa_100km <- thinning(occ_chelsa, 100)

# saving formatted datasets
write_csv(occ_chelsa, "./data/01_occ_chelsa.csv")
write_csv(occ_chelsa_10km, "./data/01_occ_chelsa_10km.csv")
write_csv(occ_chelsa_100km, "./data/01_occ_chelsa_100km.csv")


# generating validation datasets ------------------------------------------

# presence records that were discarded in the spatial thinning

occ_chelsa_100km_valid <- occ_chelsa %>%
  anti_join(occ_chelsa_100km)

# randomly selecting the same number of validation records by species
# here using 20 records because that's the smallest possible number of validation occs, for L. wellcomei

valid_spp <- list()
for(i in 1:length(study_sp)){
  valid_spp[[i]] <- filter(occ_chelsa_100km_valid, species == study_sp[i]) %>%
    slice_sample(n=20)
}
occ_valid_chelsa_std <- rbindlist(valid_spp) %>%
  mutate(pa = 1)



# thinning validation records
occ_chelsa_100km_valid <- thinning(occ_chelsa_100km_valid, 100) %>%
  mutate(pa = 1)

table(occ_chelsa_100km_valid$species)


# absence records criteria: same number as presence records,
# sampled inside model calibration area, but 
# outside a 1 arc-degree buffer from each presence record

validation <- list()
for(i in 1:length(study_sp)){
  
  # validation presences for this species
  valid_pres <- occ_chelsa_100km_valid %>%
    filter(species == study_sp[i])
  
  coordinates(valid_pres) <- c("lon", "lat")
  proj4string(valid_pres) <- wgs84
  
  # create mask for sampling absences
  calib <- readOGR(calib_files[i])
  pres_buf <- rgeos::gBuffer(valid_pres, width = 1) # 1 arc-degree buffer
  
  abs_mask <- calib - pres_buf
  proj4string(abs_mask) <- wgs84
  
  # get absence dataset and thin records
  valid_abs <- occ_chelsa_100km_valid %>%
    filter(species != study_sp[i]) %>%
    dplyr::select(-species, -pa) %>%
    mutate(species = study_sp[i]) %>%
    thinning(., 100)
  
  coordinates(valid_abs) <- c("lon", "lat")
  proj4string(valid_abs) <- wgs84
  
  # sample same number of abs as pres
  valid_abs <- rgeos::gIntersection(valid_abs, abs_mask)
  
  if(length(valid_abs) < length(valid_pres)){
    
    valid_pres <- sample_n(data.frame(valid_pres), length(valid_abs))
    
  } else {
    
    valid_abs <- sample_n(data.frame(valid_abs), length(valid_pres))
    
  }
  
  # merging final validation dataset for this species
  p <- data.frame(valid_pres) %>%
    dplyr::select(-optional)
  
  a <- data.frame(valid_abs) %>%
    remove_rownames() %>%
    rename( lon = "x", lat = "y") %>%
    mutate(species = study_sp[i],
           pa = 0)
  
  validation[[i]] <- bind_rows(p, a)

}

validation_all <- rbindlist(validation)

# saving validation datasets
write_csv(validation_all, "./data/01_validation_dataset.csv")
write_csv(occ_valid_chelsa_std, "./data/01_validation_std_presences.csv")

# counting records --------------------------------------------------------

n_records <- count_unique(occ_chelsa) %>%
  left_join(count_unique(occ_chelsa_10km), by = "species") %>%
  left_join(count_unique(occ_chelsa_100km), by = "species") %>%
  left_join(count_unique(filter(validation_all, pa == 1)), by = "species") %>%
  left_join(count_unique(filter(validation_all, pa == 0)), by = "species")
  
names(n_records) <- c("species", "n_chelsa", "n_chelsa_10km", 
                      "n_chelsa_100km", "n_validation_pres",
                      "n_validation_abs")

# saving record counts
write_csv(n_records, "./data/01_n_records.csv")
