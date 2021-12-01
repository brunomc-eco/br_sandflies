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

# maximum distance shapefiles
max_dist <- list.files("./data/shp/max_distance/", pattern = ".shp", full.names = TRUE)

# admin0 shapefile for mini-maps
admin0 <- readOGR("C:/layers/vector/GADM/GADM_Continental South America/AmSulNew/AmSulNew.shp")


# historical outputs ------------------------------------------------------

for(i in 1:length(study_sp)){
  
  # load historical final ensemble
  hist_mean <- raster(paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], 
                             "_weighted_mean_ensemble.tif")) %>%
    mask(readOGR(calib_files[i]), updatevalue = 0)
  
  png(filename = paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], "_mini_ensemble_mean.png"),
      res = 100)
  plot(hist_mean, main = paste(study_sp[i], "historical (1981-2010)"))
  plot(admin0, add=T)
  dev.off()
  
  hist_bin <- raster(paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], 
                            "_consensus_ensemble.tif")) %>%
    mask(readOGR(calib_files[i]), updatevalue = 0) 
  
  png(filename = paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], "_mini_ensemble_consensus.png"),
      res = 100)
  plot(hist_bin, main = paste(study_sp[i], "historical (1981-2010)"), legend = FALSE, col = rev(terrain.colors(2)))
  legend("topright", legend = c("Unsuitable", "Suitable"), fill = rev(terrain.colors(2)), bty = "n")
  plot(admin0, add=T)
  dev.off()
  
  writeRaster(hist_bin, filename = paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i],
                                          "_final_consensus_ensemble.tif"),
              overwrite = TRUE)
  
 
}



# suitability diffs, continuous -------------------------------------------


for(i in 1:length(study_sp)){
  
  # load historical ensemble
  hist_mean <- raster(paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], 
                             "_weighted_mean_ensemble.tif"))
  
  # calculate diffs
  diffs <- list()
  for(g in 1:length(ssp_names)){
    ssp <- raster(paste0(run_name, study_sp[i], "/future_consensus/", study_sp[i], 
                         "_", ssp_names[g], "_mean.tif"))
    diffs[[g]] <- ssp - hist_mean
    
    if (!file.exists(paste0(run_name, study_sp[i], "/diffs"))) {
      dir.create(paste0(run_name, study_sp[i], "/diffs"))
    }
    
    # save tiff file
    writeRaster(diffs[[g]], filename = paste0(run_name, study_sp[i], "/diffs/", study_sp[i],
                                              "_diff_", ssp_names[g], ".tif"), 
                overwrite = TRUE)
  }
  
  # pixel density plots
  pixels <- tibble(scenario = c(rep("ssp126", ncell(diffs[[1]])),
                                rep("ssp370", ncell(diffs[[2]])),
                                rep("ssp585", ncell(diffs[[3]]))),
                   values = c(getValues(diffs[[1]]),
                              getValues(diffs[[2]]),
                              getValues(diffs[[3]])))
  
  ggplot(pixels, aes(x = values, fill = scenario)) +
    geom_density(alpha=.2) +
    xlim(c(-1,1)) +
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1) +
    labs(title= paste(study_sp[i], "suitability diff"))
    
  ggsave(filename = paste0(run_name, study_sp[i], "/diffs/", study_sp[i], "_diff_densityplot.tif"),
         width = 4200, height = 1800, units = "px", device='tiff', dpi=300)
  
}



# suitability diffs, binary -----------------------------------------------


for(i in 1:length(study_sp)){
  
  # load historical binary
  hist_bin <- raster(paste0(run_name, study_sp[i], "/present/ensemble/", study_sp[i], 
                            "_final_consensus_ensemble.tif"))
  
  for(g in 1:length(ssp_names)){
    ssp_bin <- raster(paste0(run_name, study_sp[i], "/future_consensus/", study_sp[i], 
                            "_", ssp_names[g], "_consensus.tif"))
    ssp_bin <- ssp_bin >= 3
    
    expansion <- ((ssp_bin - hist_bin) == 1) * 4
    contraction <- ((hist_bin - ssp_bin) == 1) * 2
    stable <- (hist_bin + ssp_bin) == 2 
    
    diff_cat = expansion + contraction + stable # values: 1 = stable, 2 = contraction, 4 = expansion
    
    writeRaster(diff_cat, filename = paste0(run_name, study_sp[i], "/diffs/", study_sp[i],
                                            "_diff_categorical_", ssp_names[g], ".tif"), 
                overwrite = TRUE)
  }
  
}


# suitability diffs, max distance -----------------------------------------

for(i in 1:length(study_sp)){
  
  for(g in 1:length(ssp_names)){
    
    # load categorical diffs
    diff_cat <- raster(paste0(run_name, study_sp[i], "/diffs/", study_sp[i],
                              "_diff_categorical_", ssp_names[g], ".tif"))
    
    # load max distance mask
    dist_mask <- readOGR(max_dist[i])
    
    # apply mask 
    diff_masked <- mask(diff_cat, dist_mask, updatevalue = 0)
    
    # save output
    writeRaster(diff_masked, 
                filename = paste0(run_name, study_sp[i], "/diffs/", 
                                  study_sp[i], "_diff_categorical_maxdist", 
                                  ssp_names[g], ".tif"))
  }
  
}
