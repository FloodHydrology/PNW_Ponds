#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: PNW Ponds Initial Analysis
#Coder: Nate Jones (cnjones7@ua.edu)
#Date: 6/27/2019
#Description: The goal of this script is to kick off exploratory analysis of hydrogeomorphic 
#             characteristics of ponds created during the 1980 Mount St. Helen 
#             erruption. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#1.0 Setup Workspace------------------------------------------------------------
#Clear memory
remove(list=ls())

#download required packages
library(foreign)
library(stars)
library(sf)
library(raster)
library(tidyverse)

#Load scripts
source("R/pond_identification.R")

#Define directories
data_dir<-"//storage.research.sesync.org/njones-data/Research Projects/PNW_Ponds/spatial_data/"
scratch_dir<-"C:\\ScratchWorkspace\\"
wbt_dir<-"C:/WBT/whitebox_tools"

#Load spatial data
dem_1m<-raster(paste0(data_dir, "raw_data/1m_DEM/msh2009dem")) #data source: https://pubs.er.usgs.gov/publication/ds904
dem_10m<-raster(paste0(data_dir, "raw_data/10m_DEM/msh2009dem10")) #data source: https://pubs.er.usgs.gov/publication/ds904
ponds<-st_read(paste0(data_dir,"raw_data/pond_polygons_shape.shp")) #Emperical Data Collection 

#2.0 Prep 10m DEM for watershed delineation-------------------------------------
#Steps: 
  #(1) Filter DEM
  #(2) Burn Ponds into Raster
  #(3) Breach depressions
  #(4) Flow direction analysis
  #(5) Flow accumulation analysis

#2.1 Filter DEM-----------------------------------------------------------------
#Export rater to scratch workspace
writeRaster(dem_10m, paste0(scratch_dir,"dem.tif"), overwrite=T)

#fill missing data
system(paste(paste(wbt_dir),
             "-r=FillMissingData",
             paste0("--wd=",scratch_dir),
             "-i='dem.tif'", 
             "-o='dem_gapfilled.tif",
             "--filter=25"))

#Fill Depressions and correct flat areas
system(paste(paste(wbt_dir), 
             "-r=FillDepressions", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_gapfilled.tif'", 
             "-o='dem_filled.tif'"))

#2.2 Burn ponds into dem--------------------------------------------------------
#Read DEM back in from scratch workspace
dem<-raster(paste0(scratch_dir,"dem_filled.tif"))

#Create DEM mask
mask<-rasterize(ponds, dem, 1)

#Create minimum raster
dem_min<-mask*500 #Make the burn a constant elevation of 500 m
dem_min[is.na(dem_min)]<-0

#Replace masked location with min raster value
dem_mask<-mask*0
dem_mask[is.na(dem_mask)]<-1
dem_burn<-dem*dem_mask+dem_min

#Cleanup workspace
remove(list=c("dem","mask","dem_min",'dem_mask'))

#2.3 Breach depressions---------------------------------------------------------
#Export rater to scratch workspace
writeRaster(dem_burn, paste0(scratch_dir,"dem.tif"), overwrite=T)

#Breach Analysis of the DEM
system(paste(paste(wbt_dir),
             "-r=BreachDepressions",
             paste0("--wd=",scratch_dir),
             "--dem=dem.tif",
             "-o=dem_breach.tif"))

#2.4 Flow Accumulation and Flow Direction Anlaysis------------------------------
#Flow Direction
system(paste(paste(wbt_dir), 
             "-r=D8Pointer", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach.tif'", 
             "-o='fdr.tif'",
             "--out_type=sca"))

#Flow Accumulation
system(paste(paste(wbt_dir), 
             "-r=DInfFlowAccumulation", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach.tif'", 
             "-o='fac.tif'",
             "--out_type=sca"))

#Read fac and fdr into R environment
fac<-raster(paste0(scratch_dir,"fac.tif"))
fdr<-raster(paste0(scratch_dir,"fdr.tif"))

#3.0 Delineate entire upstream watershe for each pond---------------------------
#Create folder to export watershed shapefile too
dir.create(paste0(data_dir, 'modified_data/watersheds'))

#Create function for delineation
fun<-function(n){

  #Define pond id
  id<-ponds$pond_id[n]
  
  #isolate pond
  pond<-ponds %>% filter(pond_id==id)
  
  #Add 10 m buffer (10 m comes from the resolutin of the raster)
  pond<-st_buffer(pond, 10)
  
  #Find max fac point within pond
  fac_pond<-crop(fac, as_Spatial(pond))
  fac_pond<-mask(fac_pond, as_Spatial(pond))
  
  #create pour point
  pp<-rasterToPoints(fac_pond) %>% 
    #convert to tibble
    as_tibble() %>%
    #select point with max fac
    filter(fac==max(fac))
  
  #Make pour point an sf shape
  pp<-st_as_sf(pp, coords = c("x", "y"), crs = st_crs(pond))
  
  #Export pour point to scratch workspace
  write_sf(pp, paste0(scratch_dir,"pp.shp"), delete_layer = T)
  
  #Delineate Watershed w/ WBT
  system(paste(paste(wbt_dir), 
               "-r=Watershed", 
               paste0("--wd=",scratch_dir),
               "--d8_pntr='fdr.tif''", 
               "--pour_pts='pp.shp'",
               "-o='watershed.tif"))
  
  #Read watershed in
  w_grd<-raster(paste0(scratch_dir,"watershed.tif"))
  
  #Convert watershed to polygon
  w_shp<- w_grd %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  
  #Add id info and write to output folder
  w_shp$pond_id<-id
  write_sf(w_shp, paste0(data_dir, "modified_data/watersheds/",id,".shp"))
}

#Apply function (note, because you are using WBT, this has to be single threaded)
lapply(fun, seq(1,nrow(ponds)))




