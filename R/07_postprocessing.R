# ENM output post-processing and  area calculations
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(raster)
library(rgdal)
library(ggplot2)


# load data and set values ------------------------------------------------

# name the folder where previous final models were saved
run_name <- c("./outputs/models_chelsa_100km/")

# scenario and gcm names
gcm_names <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
ssp_names <- c("ssp126", "ssp370", "ssp585")


# species names
study_sp <- c("L_complexa", "L_cruzi", "L_flaviscutellata", "L_intermedia",
              "L_longipalpis", "L_migonei", "L_neivai", "L_umbratilis",
              "L_wellcomei", "L_whitmani")


# model calibration area by species
calib_files <- list.files("./data/shp/calib", pattern = ".shp", full.names = TRUE)

# admin1 shapefile for mapping
admin1 <- readOGR("C:/layers/vector/GADM_Continental South America/AmSulNew/AmSulNew.shp")



# historical outputs ------------------------------------------------------

# reclass matrix for final outputs - 0s to NAs to facilitate count and plot
rclmat <- matrix(c(0, NA), ncol=2, byrow = TRUE)

for(i in 1:length(study_sp)){
  
  # load historical final ensemble
  hist_mean <- raster(paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], 
                             "_weighted_mean_ensemble.tif")) %>%
    mask(readOGR(calib_files[i]))
  
  png(filename = paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], "_mini_ensemble_mean.png"),
      res = 100)
  plot(hist_mean, main = paste(study_sp[i]), "historical (1981-2010)")
  plot(admin1, add=T)
  dev.off()
  
  hist_bin <- raster(paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], 
                            "_consensus_ensemble.tif")) %>%
    mask(readOGR(calib_files[i])) %>%
    reclassify(rclmat)
  
  png(filename = paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], "_mini_ensemble_consensus.png"),
      res = 100)
  plot(hist_bin, main = paste(study_sp[i], "historical (1981-2010)"), legend = FALSE, col = "#00A600")
  legend("topright", legend = c("Unsuitable", "Suitable"), fill = rev(terrain.colors(2)))
  plot(admin1, add=T)
  dev.off()
  
  writeRaster(hist_bin, filename = paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i],
                                          "_final_consensus_ensemble.tif"),
              overwrite = TRUE)
  
 
}



# suitability differences -------------------------------------------------


for(i in 1:length(study_sp)){
  
  # load historical ensemble
  hist_mean <- raster(paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], 
                             "_weighted_mean_ensemble.tif"))
  # load future ssps
  
  ssp126 <- raster(paste0(run_name, study_sp[i], "/future_consensus/", study_sp[i], 
                       "_ssp126_mean.tif"))
  ssp370 <- raster(paste0(run_name, study_sp[i], "/future_consensus/", study_sp[i], 
                          "_ssp370_mean.tif"))
  ssp585 <- raster(paste0(run_name, study_sp[i], "/future_consensus/", study_sp[i], 
                          "_ssp585_mean.tif"))
  
  # calculate differences in pixel values
  
  diff_ssp126 = ssp126 - hist_mean
  diff_ssp370 = ssp370 - hist_mean
  diff_ssp585 = ssp585 - hist_mean
  
  # save diff outputs
  if (!file.exists(paste0(run_name, study_sp[i], "/diffs"))) {
    dir.create(paste0(run_name, study_sp[i], "/diffs"))
  }
  
  writeRaster(diff_ssp126, filename = paste0(run_name, study_sp[i], "/diffs/",
                                             study_sp[i], "_diff_ssp126.tif"), 
              overwrite = TRUE)
  
  writeRaster(diff_ssp370, filename = paste0(run_name, study_sp[i], "/diffs/",
                                             study_sp[i], "_diff_ssp370.tif"), 
              overwrite = TRUE)
  
  writeRaster(diff_ssp585, filename = paste0(run_name, study_sp[i], "/diffs/",
                                             study_sp[i], "_diff_ssp585.tif"), 
              overwrite = TRUE)
  
  # set dataframe with pixel values for plots
  
  #pixels <- tibble(id = seq(1:ncell(diff_ssp126)),
  #                 diff_ssp126_vals = getValues(diff_ssp126),
  #                 diff_ssp370_vals = getValues(diff_ssp370),
  #                 diff_ssp585_vals = getValues(diff_ssp585))
  
  pixels <- tibble(scenario = c(rep("ssp126", ncell(diff_ssp126)),
                                rep("ssp370", ncell(diff_ssp370)),
                                rep("ssp585", ncell(diff_ssp585))),
                   values = c(getValues(diff_ssp126),
                              getValues(diff_ssp370),
                              getValues(diff_ssp585)))
  
  # plot density (histograms)
  ggplot(pixels, aes(x = values, fill = scenario)) +
    geom_density(alpha=.2) +
    xlim(c(-1,1)) +
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1) +
    labs(title= paste(study_sp[i], "suitability diff"))
  # save
  ggsave(filename = paste0(run_name, study_sp[i], "/diffs/", study_sp[i], "_densityplot.tif"),
         width = 4200, height = 1800, units = "px", device='tiff', dpi=300)
}

