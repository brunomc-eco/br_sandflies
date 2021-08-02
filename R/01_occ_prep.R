# Preparing occurrence data
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(ggplot2)
library(spThin)
library(data.table)


# load data ---------------------------------------------------------------


# select study species
study_sp <- c("L_neivai", "L_migonei", "L_longipalpis", "L_whitmani", 
              "L_umbratilis", "L_flaviscutellata", "L_cruzi", 
              "L_intermedia", "L_complexa", "L_wellcomei")

# function for quick visual inspection of records
visual_check <- function(x){
  ggplot() +
    borders("world", colour="white", fill="gray50") +
    geom_point(data = x, 
               aes(x = lon, y = lat, colour = species), size = 1) +
    coord_sf(xlim = c(min(x$lon) -0.5,
                      max(x$lon) +0.5), 
             ylim = c(min(x$lat) -5,
                      max(x$lat) +5))
}


# loading occurrence data files, full and historical (1970-2000)
occ_full <- read_csv("./data/raw/ocorrencias_total_limpando.csv")
occ_hist <- read_csv("./data/raw/ocorrencias_total_limpando_1970-2000.csv")


# retaining only records from study species, with unique coordinates
occ_full <- occ_full %>%
  filter(lon_dec != "NA",
         lat_dec != "NA",
         considered.species %in% study_sp) %>%
  rename(species = considered.species,
         lon = lon_dec,
         lat = lat_dec) %>%
  select(species, lon, lat) %>%
  distinct()

occ_hist <- occ_hist %>%
  filter(lon_dec != "NA",
         lat_dec != "NA",
         considered.species %in% study_sp) %>%
  rename(species = considered.species,
         lon = lon_dec,
         lat = lat_dec) %>%
  select(species, lon, lat) %>%
  distinct()


# see how it looks

visual_check(occ_full)
visual_check(occ_hist)


# saving formatted data
write_csv(occ_full, "./data/01_occ_full.csv")
write_csv(occ_hist, "./data/01_occ_hist.csv")


# thin records ------------------------------------------------------------


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


occ_full_10km <- thinning(occ_full, 10)
occ_full_100km <- thinning(occ_full, 100)
occ_hist_10km <- thinning(occ_hist, 10)
occ_hist_100km <- thinning(occ_hist, 100)


# generating validation datasets ------------------------------------------

occ_hist_valid_100 <- occ_hist %>%
  anti_join(occ_hist_100km)

write_csv(occ_hist_valid_100, "./data/01_occ_hist_100km_valid.csv")


# counting records --------------------------------------------------------


# counting records by species in different sets

count_unique <- function(x){
  x %>%
    group_by(species) %>%
    summarize(count = n())
}


n_records <- count_unique(occ_full) %>%
  left_join(count_unique(occ_hist), by = "species") %>%
  left_join(count_unique(occ_full_10km), by = "species") %>%
  left_join(count_unique(occ_full_100km), by = "species") %>%
  left_join(count_unique(occ_hist_10km), by = "species") %>%
  left_join(count_unique(occ_hist_100km), by = "species")

names(n_records) <- c("species", "n_full", "n_hist", "n_full_10km", 
                      "n_full_100km", "n_hist_10km", "n_hist_100km")



# saving outputs ----------------------------------------------------------

write_csv(occ_full_10km, "./data/01_occ_full_10km.csv")
write_csv(occ_full_100km, "./data/01_occ_full_100km.csv")
write_csv(occ_hist_10km, "./data/01_occ_hist_10km.csv")
write_csv(occ_hist_100km, "./data/01_occ_hist_100km.csv")
write_csv(n_records, "./data/01_n_records.csv")