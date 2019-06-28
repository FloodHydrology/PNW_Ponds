#2.0 Clip DEM to NHDPlus "Neighborhood"-----------------------------------------
#Here, we are going to clip the DEM to relevant drainage basins. This involves 
#(1) selecting the NHDPlusCatchments (hillslopes) that intersect with the ponds, 
#(2) selecting the associated NHDStreamLine (stream), 
#(3) identifying upstream NHDPlusStreamLines, 
#(4) selecting NHDPlusCatchments that intersect those NHDPlusStreamLines, and 
#(5) clipping the DEM to those catchments. 
#Notably, there's most definitely a better way to do all of this using graph theory. 

#Download data
streams<-st_read(paste0(data_dir,"modified_data/NHDPlus_HR_FlowLines.shp")) #data source: https://www.usgs.gov/core-science-systems/ngp/national-hydrography/access-national-hydrography-products
hillslopes<-st_read(paste0(data_dir,"modified_data/NHDPlusCatchment.shp")) #data source: https://www.usgs.gov/core-science-systems/ngp/national-hydrography/access-national-hydrography-products
connectivity<-read.dbf(paste0(data_dir,"modified_data/NHDPlusFlow.dbf")) %>% as_tibble() #data source: https://www.usgs.gov/core-science-systems/ngp/national-hydrography/access-national-hydrography-products 

#Standardize CRS accross spatial datasets
p<-st_crs(ponds)
streams<-st_transform(streams, p)
hillslopes<-st_transform(hillslopes, p)
streams<-st_zm(streams)
hillslopes<-st_zm(hillslopes)


#2.1 Select relevent hillslope and streamline areas~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Select hillslopes that intersect with pond polygons
h<-hillslopes[ponds,] 

#Select streams taht intersect with hillslopes
s<-streams[h,] %>% 
  as_tibble() %>% 
  select(Permanent_) %>% 
  rename(UID = Permanent_) %>%
  mutate(UID = as.numeric(paste(UID))) %>%
  na.omit()

#2.2 Select upstream streamlines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create edge list
connectivity <- connectivity %>% 
  select(FromPermID, ToPermID) %>%
  rename(FROM = FromPermID, 
         TO = ToPermID) %>%
  mutate(FROM = as.numeric(paste(FROM)), 
         TO = as.numeric(paste(TO))) %>%
  na.omit()

#Conduct while statement to find upstream stream reaches
m<-nrow(s)
while(m>0){
  #Find upstream reaches
  temp<-s %>%
    rename(TO = UID) %>%
    left_join(.,connectivity) %>%
    select(FROM) %>%
    filter(FROM>0) %>%
    rename(UID=FROM) %>%
    distinct(.)
  
  #Bind rows
  s<-bind_rows(s,temp) %>% distinct(.)
  
  #Count number of upstream reaches.If same as last iteration, then kill the loop
  if(nrow(s)==m){
    m<-0}else{
      m<-nrow(s)
    }
  print(m)
}

#2.3 Find upstream areas~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Select upstream flowlines
streams<-streams %>%
  rename(UID = Permanent_) %>%
  mutate(UID = as.numeric(paste(UID))) %>%
  right_join(.,s)

#Intersect streams with hillslopes
hillslopes<-hillslopes[streams, ] 

#2.4 Crop DEM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create mask
mask<-st_union(hillslopes)
mask<-st_buffer(mask, dist = 1000)
mask<-as_Spatial(mask)

#Crop DEM
dem<-crop(dem, mask)
dem<-mask(dem, mask)
remove(mask)

#2.5 Clean up workspace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
writeRaster(dem, paste0(data_dir,"modified_data/dem_clipped.tif"))
remove(list=ls()[ls()!='dem' &
                   ls()!='streams' &
                   ls()!='ponds' &
                   ls()!='data_dir' &
                   ls()!='scratch_dir' &
                   ls()!='wbt_dir'])