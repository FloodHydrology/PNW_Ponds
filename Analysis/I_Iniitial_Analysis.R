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
library(fasterize)
library(stars)
library(sf)
library(raster)
library(tidyverse)
library(parallel)

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

#3.0 Delineate entire upstream watershed for each pond--------------------------
#Create folder to export watershed shapefile too
dir.create(paste0(data_dir, 'modified_data/watersheds'))

#Create function for delineation
fun<-function(n){

  #Step 1: Delineate watershed with low res DEM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Define pond id
  id<-ponds$pond_id[n]
  
  #isolate pond
  pond<-ponds %>% filter(pond_id==id)
  
  #Add 3 m buffer (10 m comes from the resolutin of the raster)
  pond<-st_buffer(pond, 10/3)
  
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
  
  #Step 2: Re-delineate watershed with high res DEM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #Clip high resolution DEM with w_shp
  mask<-st_buffer(w_shp, 100)
  dem_hr<-crop(dem_1m, as_Spatial(mask))

  #Write raster to scratch workspace
  writeRaster(dem_hr,paste0(scratch_dir,"dem_hr.tif"), overwrite=T)
  
  #Breach Analysis of the DEM
  system(paste(paste(wbt_dir),
               "-r=BreachDepressions",
               paste0("--wd=",scratch_dir),
               "--dem=dem_hr.tif",
               "-o=dem_breach_hr.tif"))
  
  #Flow Direction
  system(paste(paste(wbt_dir), 
               "-r=D8Pointer", 
               paste0("--wd=",scratch_dir),
               "--dem='dem_breach_hr.tif'", 
               "-o='fdr_hr.tif'",
               "--out_type=sca"))
  
  #Flow Accumulation
  system(paste(paste(wbt_dir), 
               "-r=DInfFlowAccumulation", 
               paste0("--wd=",scratch_dir),
               "--dem='dem_breach_hr.tif'", 
               "-o='fac_hr.tif'",
               "--out_type=sca"))
  
  #Read fac and fdr into R environment
  fac_hr<-raster(paste0(scratch_dir,"fac_hr.tif"))
  fdr_hr<-raster(paste0(scratch_dir,"fdr_hr.tif"))
  
  #Find max fac point within pond
  fac_pond_hr<-crop(fac_hr, as_Spatial(pond))
  fac_pond_hr<-mask(fac_pond_hr, as_Spatial(pond))
  
  #create pour point
  pp_hr<-rasterToPoints(fac_pond_hr) %>% 
    #convert to tibble
    as_tibble() %>%
    #select point with max fac
    filter(fac_hr==max(fac_hr))
  
  #Make pour point an sf shape
  pp_hr<-st_as_sf(pp_hr, coords = c("x", "y"), crs = st_crs(pond))
  
  #Export pour point to scratch workspace
  write_sf(pp_hr, paste0(scratch_dir,"pp_hr.shp"), delete_layer = T)
  
  #Delineate Watershed w/ WBT
  system(paste(paste(wbt_dir), 
               "-r=Watershed", 
               paste0("--wd=",scratch_dir),
               "--d8_pntr='fdr_hr.tif''", 
               "--pour_pts='pp_hr.shp'",
               "-o='watershed_hr.tif"))
  
  #Read watershed in
  w_hr_grd<-raster(paste0(scratch_dir,"watershed_hr.tif"))
  
  #Convert watershed to polygon
  w_hr_shp<- w_hr_grd %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  
  #Add id info and write to output folder
  w_hr_shp$pond_id<-id
  write_sf(w_hr_shp, paste0(data_dir, "modified_data/watersheds/",id,".shp"), delete_layer = T)
}

#Apply function (note, because you are using WBT, this has to be single threaded)
lapply(seq(1,nrow(ponds)),fun)

#4.0 Delineate pond subshed-----------------------------------------------------
#Create folder to export subshed shapefile too
dir.create(paste0(data_dir, 'modified_data/subsheds'))

#Create function to identify subsheds
fun<-function(n){

  #Define pond id
  id<-ponds$pond_id[n]
  
  #isolate pond
  pond<-ponds %>% filter(pond_id==id)
  
  #load pond watershed
  ws_shp<-st_read(paste0(data_dir, "modified_data/watersheds/",id,".shp"))

  #Determine if there are any upstream ponds
  upstream_ponds<-ponds[ws_shp,] %>% filter(pond_id!=id)

  #If there are upstream ponds, clip them out
  if(nrow(upstream_ponds)>0){
    #Create loop to crop polygons that are completely contained by upstream watershed
    for(i in 1:nrow(upstream_ponds)){
      upstream_subshed<-st_read(paste0(data_dir, "modified_data/watersheds/",upstream_ponds$pond_id[i],".shp"))
      if(nrow(upstream_subshed[pond,])==0){
        ws_shp<-st_difference(ws_shp, upstream_subshed)
      }
    }
  }
  
  #Export subshed shape
  ws_shp$pond_id<-id
  write_sf(ws_shp, paste0(data_dir, "modified_data/subsheds/",id,".shp"), delete_layer = T)
}
  
#Apply function 
lapply(seq(1,nrow(ponds)),fun)

#5.0 Inundation analysis--------------------------------------------------------
#Create folders to export bathymetric relationships and estimated inundation too
dir.create(paste0(data_dir, 'modified_data/bathymetry_relationships'))
dir.create(paste0(data_dir, 'modified_data/inundation/'))

#Create function to inundate ponds
fun<-function(n){
  #Setup workspace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Define pond id
  id<-ponds$pond_id[n]
  
  #isolate pond
  pond<-ponds %>% filter(pond_id==id)
  
  #isolate subshed
  subshed<-st_read(paste0(data_dir,"modified_data/subsheds/",id,'.shp'))
  
  #Crop dem
  dem<-crop(dem_1m, as_Spatial(subshed))
  dem<-mask(dem, as_Spatial(subshed))
  
  #Inundate Analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create Minimum Raster
  dem_min<-dem*0+minValue(dem)
  dem_min@crs<-dem@crs
  
  #Create minimum point
  pnt_min<-cellStats(dem, min)
  pnt_min<-rasterToPoints(dem,spatial=T, fun=function(x) x ==pnt_min)
  
  #Create function to return conditional raster
  Con<-function(condition, trueValue, falseValue){
    return(condition * trueValue + (!condition)*falseValue)
  }
  
  #Create function to calcluate inundation area, volume, and spill boundary length
  inundate<-function(z){
    
    #define inundated area [area connected to pnt_min]
    area<-Con(dem>(dem_min+z),0,1)
    group<-raster::extract(x=area, y=pnt_min)
    area[area!=group]<-0
    area[area==group]<-1
    
    #Create metrics to estimate area, volume, and spill boundary length
    volume<-(((z+dem)-dem_min)*area)*res(area)[1]*res(area)[2]
    outflow<-cellStats(area*boundaries(dem_min, type="inner"), 'sum')
    
    #Estimate number of inundated cells within delineated pond boundary
    mask_inside<-fasterize(pond, dem)
    mask_outside<-mask_inside*0
    mask_outside[is.na(mask_outside)]<-1
    inside<-mask_inside*area
    outside<-mask_outside*area
    
    #Export Data
    c(z, #inundation depth
      cellStats(area, 'sum')*res(area)[1]*res(area)[2], #area (m^2)
      cellStats(volume, 'sum'), #volume (m^3)
      cellStats(inside, 'sum'), #area inside of delineated pond
      cellStats(outside, 'sum'), #area outside of the delineated pond
      outflow) #Outflow length (increments = raster cell resolution)
  }
  
  #Create function to calculate inundation area/volume
  df<-lapply(seq(0,3,0.01),inundate)
  df<-do.call(rbind, df)
  df<-data.frame(df)  
  colnames(df)<-c("z", "area","volume", 'area_inside', 'area_outside', "outflow_length")
  
  #Export to new csv file
  write_csv(df, paste0(data_dir,'modified_data/bathymetry_relationships/',id,".csv"))
  
  #Estimate max inundation based on emperical pond delineation~~~~~~~~~~~~~~~~~~
  #Select stage where the number of "outside" cells and difference in area estimates
  #(e.g., emperical vs modeled) is minimized  
  df <- df %>%
    mutate(pond_id = id,
           area_percent_diff = abs(area-as.numeric(paste(st_area(pond))))/area*100, 
           rank      = rank(area_percent_diff)) %>%
    filter(rank<=15) %>%
    arrange(rank) %>%
    filter(area_outside == min(area_outside, na.rm=T))
  
  #if nrow(df)>1, then limit to lowest error
  if(nrow(df)>1){df<-df[1,]}
  
  #Creat inundation shape and export 
  inundation<-Con(dem>(dem_min+df$z),0,1) 
  inundation[inundation==0]<-NA
  inundation<-inundation %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  st_write(inundation, paste0(data_dir, "modified_data/inundation/",id,'.shp'), delete_layer = T)

  #Export metrics to global environment
  df
}

#Apply function in parallel (note, use mclapply if using a non-PC -- much easier!)
t0<-Sys.time()
n.cores<-detectCores()
cl <- makePSOCKcluster(n.cores) #Create Clusters
clusterEvalQ(cl, {library(fasterize)
                  library(stars)
                  library(sf)
                  library(raster)
                  library(tidyverse)})  #Send clusters the required libraries
clusterExport(cl, c('fun', 'data_dir', 'ponds', 'dem_1m'), env=.GlobalEnv)  #Send Clusters function with the execute function
x<-parLapply(cl, seq(1,nrow(ponds)), fun)  #Execute function
stopCluster(cl)  #Turn clusters off
tf<-Sys.time()
tf-t0

#Gather information
output <- x %>% bind_rows(.) %>% select(pond_id, z, area, volume, area_percent_diff)

#bind to ponds
ponds<-ponds %>% 
  left_join(.,output) %>% 
  rename(pond_depth = z, 
         pond_area =  area, 
         pond_volume = volume)

#6.0 Estimate Metrics-----------------------------------------------------------
#Create function to estimate metrics for individual ponds
fun<-function(n){
  #Setup workspace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Define pond id
  id<-ponds$pond_id[n]
  
  #isolate pond
  pond<-ponds %>% filter(pond_id==id)
  
  #isolate subshed
  subshed<-st_read(paste0(data_dir,"modified_data/subsheds/",id,'.shp'))
  
  #Estimate metrics from watershed delineations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Read watershed/subshed shapefiles
  subshed<-read_sf(paste0(data_dir, 'modified_data/subsheds/',id,".shp")) %>%
    st_union(.)
  watershed<-read_sf(paste0(data_dir, 'modified_data/watersheds/',id,".shp"))
  
  #Estiamte pond metrics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subshed_area<-st_area(subshed)
  pond_specific_volume<-as.numeric(paste(pond$pond_volume/subshed_area))
  
  #Estimate total watershed metrics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  watershed_area<-st_area(watershed)
  ponds_watershed<-ponds[watershed,]
  total_pond_area<-sum(ponds_watershed$pond_area, na.rm=T)
  total_pond_volume<-sum(ponds_watershed$pond_volume, na.rm = T)
  total_pond_specific_volume<-as.numeric(paste(total_pond_volume/watershed_area))
  
  #Estimate upstream watershed metrics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  upstream_watershed_area<-watershed_area-subshed_area
  upstream_pond_area<-total_pond_area- pond$pond_area
  upstream_pond_volume<-total_pond_volume - pond$pond_volume
  upstream_pond_specific_volume<-as.numeric(paste(upstream_pond_volume/upstream_watershed_area))
    upstream_pond_specific_volume[is.na(upstream_pond_specific_volume)]<-0
  
  #Export results
  tibble(pond_id = id, 
         subshed_area, pond_specific_volume,
         upstream_watershed_area, upstream_pond_area, upstream_pond_volume, upstream_pond_specific_volume,
         watershed_area, total_pond_area, total_pond_volume, total_pond_specific_volume)
}

#apply function
output<-lapply(seq(1, nrow(ponds)), fun) %>% bind_rows(.) 

#join to pond dataset
ponds<-left_join(ponds, output)

#Export outputs
dir.create(paste0(data_dir, 'output/'))
st_write(ponds, paste0(data_dir,"output/ponds.shp"))
output<-as_tibble(ponds) %>% select(-geometry)
write_csv(output, paste0(data_dir,"output/ponds_metrics.csv"))
  