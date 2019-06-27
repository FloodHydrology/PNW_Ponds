#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: PNW Ponds Initial Analysis
#Coder: Nate Jones (cnjones7@ua.edu)
#Date: 6/27/2019
#Description: The goal of this script is to kick off exploratory analysis of hydrogeomorphic 
#             characteristics of ponds created during the 1980 Mount St. Helen 
#             erruption. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Setup Workspace----------------------------------------------------------------
#Clear memory
remove(list=ls())

#download required packages
library(tidyverse)
library(raster)
library(sf)

#Define directories
data_dir<-

