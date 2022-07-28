# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:46:06 2022

@author: mv296


Script to prepare the predictor layers for Sumatra

# one script based on Hansen
# one script based on JRC

"""


import numpy as np

import os

from osgeo import ogr



from glob import glob

import csv
from itertools import chain
import math

import macpyver as mp

# unit conversion
km2_ha = 100
ha_km2 = 0.01
km2_Mha = 0.0001
m2_km2 = 0.000001

# resolution (180m)
res = str(180)


out_path = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res'

workaround = r"C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts"



"""
# Run this if you want to jump in the middle of the script

extent = "-3292738.004 962289.385 -1721878.004 2256489.385"
base_map = mp.tif.read('C:\\Users\\mv296\\work\\Sumatra\\data\\model_input\\repro_res\\sumatra_complete_180_m_repro.tif', 1)


forest_2000 = mp.tif.read(out_path + "\\tree_cover_repro_70_clip_" + res + "m_repro_res.tif", 1)
forest_loss = mp.tif.read(out_path + "\\Hansen_GFC-2021-v1.9_lossyear_"+ res + "m_repro_res.tif" , 1)


"""




#-----------#
# Baseshape #
#-----------#

# process the base shape
# select the regions by name
# Aceh, Sumatera Utara, Riau, Sumatera Barat, Jambi, Bengkulu, Sumatera Selatan, 
# Bangka Belitung, Lampung, Kepulauan Riau

# check how much forest in each


unit_filter = str("NAME_1 LIKE '%Aceh%' OR NAME_1 LIKE '%Sumatera%'  OR NAME_1 LIKE '%Jambi%' OR NAME_1 LIKE '%Bengkulu%' OR NAME_1 LIKE '%Bangka Belitung%' OR NAME_1 LIKE '%Lampung%'  OR NAME_1 LIKE '%Riau%' AND NAME_1 NOT LIKE '%Kepulauan Riau%'")

# i am excluding  OR NAME_1 LIKE '%Kepulauan Riau%' because it includes islands in the West of Sarawak and is just a bit messy
gadm_path = r'C:\Users\mv296\work\Indonesia\gadm\gadm36_IDN_1\gadm36_IDN_1.shp'              
sumatra_shape_path = r'C:\Users\mv296\work\Sumatra\data\gadm\sumatra_complete_shape.shp'        
    
mp.get_stdout
print("""ogr2ogr  -overwrite -a_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -s_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -sql "SELECT * FROM  gadm36_IDN_1 WHERE """+ unit_filter  + """ " """ + sumatra_shape_path + " " + gadm_path)

sumatra_basename = os.path.basename(sumatra_shape_path)[:-4]

# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + 
    sumatra_shape_path[:-4] + "_repro.shp" + " " + sumatra_shape_path)
# the file sumatra_shape_path[:-4] + "_repro.shp" will be used for regions!

mp.get_stdout("ogrinfo -so -al " +  sumatra_shape_path[:-4] + "_repro.shp" )

# extract extent:
pre_extent = "-3292738.004418 962243.639245 -1721896.360977 2256489.384612"

mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l "+ sumatra_basename +"_repro -tr " + res + " -" + res + " -te " + pre_extent +" "+  sumatra_shape_path[:-4] + "_repro.shp " +  out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")

mp.get_stdout("gdalinfo " +  out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")

extent = "-3292738.004 962289.385 -1721878.004 2256489.385"


# save the path here
base_path = out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif"

base_map = mp.tif.read(base_path, 1)

# outer shape (digitized in ArcGIS for clipping and others)

outer_base_shape = r'C:\Users\mv296\work\Sumatra\data\Sumatra_outer_shape\Sumatra_outer_shape.shp'


#--------#
# Forest #
#--------#

# primary forest
# PRIMARY FOREST LAYER
# reproject primary forest layer from MARGONO / Hansen to clip the forest with

mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent +
        r" C:\Users\mv296\work\Indonesia\data\Margono_primary_forest_change\timeseq_change00_12.tif"  + " " + out_path + "//timeseq_change00_12_" + res + "m_repro_res.tif" )  
        
# recode
primary = mp.tif.read(out_path + "\\timeseq_change00_12_" + res + "m_repro_res.tif", 1)
# only category 3 is non-forets because all the other change will be captured within the change  
# from hansen
primary = np.where((primary == 3) | (primary ==0), 0, 1)      

mp.tif.write(out_path + "\\timeseq_change00_12_" + res + "m_repro_res.tif",
                          out_path +"\\primary_forest_"+res+"m_repro_res.tif",
                       primary, nodata = -9999, option='compress=deflate')
                       

# continue checking this

# mangrove forest from mangrove atlas

# I used this layer C:\Users\mv296\work\Wallacea\deforestation_model\data\Mangroves\GMW_001_GlobalMangroveWatch\GMW_001_GlobalMangroveWatch\01_Data\GMW_2016_v2.shp
# and the 
# to clip in ArcGIS and produce the Sumatra layer
mangrove_path = r'C:\Users\mv296\work\Sumatra\data\forest\Mangroves_GMW_2016_v2_Sumatra.shp'

# reproject 
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
mangrove_path[:-4] + "_repro.shp "  + mangrove_path )

# then rasterize
mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l Mangroves_GMW_2016_v2_Sumatra_repro -tr " + res + " -" + res + " -te " + extent +" "+ mangrove_path[:-4] + "_repro.shp " +  out_path + "\\" + "GMW_2016_v2_"+ res + "_m_repro.tif")


# combine mangrove and primary forest and Hansen layer// what about David



mangrove = mp.tif.read(mangrove_path[:-4] + "_repro.shp ", 1)
primary = np.where((primary == 1) | (mangrove == 1), 1, -9999)

# forest_loss

forest_loss_path = r'C:\Users\mv296\work\Indonesia\deforestation\forest_loss_Hansen\v1.9\Hansen_GFC-2021-v1.9_lossyear.tif'

file_list = glob(r'C:\Users\mv296\work\Indonesia\deforestation\forest_loss_Hansen\v1.9\Hansen_GFC-2021-v1.9_lossyear*.tif')

files_string = " ".join(file_list)

command = workaround + "\gdal_merge.py -o "+ forest_loss_path +" -of gtiff " + files_string

mp.get_stdout(command)


mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
              "-tr " + res + " -" + res + " -te " + extent + " "+
              forest_loss_path + " " +  out_path + "\\Hansen_GFC-2021-v1.9_lossyear_"+ res + "m_repro_res.tif" )  
        





# forest cover

# forest cover
# combining all of Indonesia, just because its easier
forest_cover_path = r'C:\Users\mv296\work\Indonesia\deforestation\forest_loss_Hansen\v1.9\Hansen_GFC-2021-v1.9_treecover2000.tif'

file_list = glob(r'C:\Users\mv296\work\Indonesia\deforestation\forest_loss_Hansen\v1.9\Hansen_GFC-2021-v1.9_treecover2000_*.tif')

files_string = " ".join(file_list)

command = workaround + "\gdal_merge.py -o "+ forest_cover_path +" -of gtiff " + files_string

mp.get_stdout(command)


mp.get_stdout("""gdalwarp -r "bilinear" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
              "-tr " + res + " -" + res + " -te " + extent + " "+
              forest_cover_path + " " +  out_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif" )  
        

forest_2000 = mp.tif.read( out_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif", 1)


# then fix cover > 70% 
forest_2000_70 = np.where(forest_2000 >= 70, 1, -9999) 

np.unique(forest_2000_70)


# mp.tif.write(out_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif", 
#              out_path +"\\tree_cover_repro_70.tif",
#                        forest_2000_70,  
#                        nodata = -9999,
#                        option='compress=deflate')



forest_2000_clip = np.where(primary ==1, forest_2000_70, -9999)
 
mp.tif.write(out_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif", 
             out_path + "\\tree_cover_repro_70_clip_" + res + "m_repro_res.tif", 
             forest_2000_clip, 
             nodata = -9999,
             option='compress=deflate')


# prepare forest loss layers for model

#out_path_small_islands = "C:\Users\mv296\work\Wallacea\deforestation_model\data/wallacea_delim/wallacea_no_small_islands.shp"


forest_2000 = mp.tif.read(out_path + "\\tree_cover_repro_70_clip_" + res + "m_repro_res.tif", 1)
forest_loss = mp.tif.read(out_path + "\\Hansen_GFC-2021-v1.9_lossyear_"+ res + "m_repro_res.tif" , 1)
np.unique(forest_loss)
# test 
base_map.shape == forest_2000.shape == forest_loss.shape



 # prepare two layers:
    # one is:  0 = previous deforestation, 1= forest, -9999 = sea, non-forest vegetation
    #  one is: -9999 = not deforested or not considered, 1 = forest converted to non-forest 
#PROBLEM: IN NON PRIMARY FOREST THERE IS FOREST
forest_1 = np.where((forest_2000 == 0), -9999, forest_2000)
    # this is previous deforestation
forest_1 = np.where((forest_loss >= 1) & (forest_loss <= 16) & (forest_1 != -9999), 0, forest_1)
    # this is current forest
forest_1 = np.where((forest_loss >= 17) & (forest_2000 == 1), 1, forest_1)
forest_1 = np.where((base_map == 0), -9999, forest_1)

mp.tif.write(base_path,
                         out_path + '/forest_2016_21_'+str(res) + 'm_repro_res.tif', 
                            forest_1, nodata=-9999, option='compress=deflate')
                         
    # layer 2 (deforestation in 2016-2018)
forest_2 = np.where((base_map == 0)|(forest_2000 == 0), -9999, forest_loss)
forest_2 = np.where((forest_2 >= 17) & (forest_2000 == 1), 1, -9999)
forest_2 = np.where((base_map == 0), -9999, forest_2)

mp.tif.write(base_path,
                             out_path + '/deforestation_2016_21_'+str(res) + 'm_repro_res.tif',
                             forest_2, nodata=-9999, option='compress=deflate')
     
     
# prep

# TMF layer

# merge TMF layer, think which layers and how? 

# yearly layer or the year layer ? because that is when first deforested
# merge the following layers 
years = range(1990, 2022)
tile_ids = ["ID61_N10_E100", "ID60_N10_E90", "ID39_N0_E100", "ID38_N0_E90"]
# we need these because we only want the ones from sumatra
tmf_path = r'N:\Landcover\Indonesia\TMF_JRC_forest'

middle_path = r'\forobs\products\tmf_v1\AnnualChange'

# 0 - NA/out of info
# 1 undisturbed tropical moist forest
# 2 degraded tropical moist forest
# 3 deforested land
# 4 forest regrowth
# 5 is water, so don't wonder about this
# 6 other landcover


# I am only defining 1 and 2 as forest, not regrowth, deforested land, water or other

# for each year assemble the four tiles

tmf_out_path = r'N:\Landcover\Sumatra\processed'


for year in years:
    print(year)
    path_list = []
    for tile_id in tile_ids:
            path = tmf_path + "\\JRC_TMF_AnnualChange_v1_ASI_" + tile_id + "\\" + middle_path + "\\JRC_TMF_AnnualChange_v1_"+str(year)+"_ASI_"+tile_id+".tif"
            path_list.append(path)
            # here continue merging
            files_string = " ".join(path_list)
            command = workaround + "\gdal_merge.py -o "+ tmf_out_path + "\\JRC_TMF_AnnualChange_v1_Sumatra_"+str(year)+".tif" + " -of gtiff " + files_string
            mp.get_stdout(command)
            # reproject
            mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
                          "-tr " + res + " -" + res + " -te " + extent + " "+
                          tmf_out_path+ "\\JRC_TMF_AnnualChange_v1_Sumatra_"+str(year)+".tif "  +  tmf_out_path+ "\\JRC_TMF_AnnualChange_v1_Sumatra_"+str(year)+"_"+ res + "m_repro_res.tif" )  
        
    
# use the first year 1990, then 2017-2021  to construct layer 1 and to
# and clip with basemap 

 # prepare two layers:
    # one is:  0 = previous deforestation, 1= forest, -9999 = sea, non-forest vegetation in trainingsperiod
    #  one is: -9999 = not deforested or not considered, 1 = forest converted to non-forest in trainingsperiod
        
# deforestation including 2016 is 0, where there is forest is 1, rest NA
# so pixels that were there in first year, but not in 2017, that is 0 
# pixels that are there in 2017 are 1

forest_1990 = mp.tif.read(tmf_out_path+ "\\JRC_TMF_AnnualChange_v1_Sumatra_1990_"+ res + "m_repro_res.tif", 1)
np.unique(forest_1990)
forest_1990 = np.where((forest_1990 == 1) | (forest_1990 == 2), 1, 0)
np.unique(forest_1990)


forest_2017 = mp.tif.read(tmf_out_path+ "\\JRC_TMF_AnnualChange_v1_Sumatra_2017_"+ res + "m_repro_res.tif", 1)
np.unique(forest_2017) 
forest_2017 = np.where((forest_2017 == 1) | (forest_2017 == 2), 1, 0)
np.unique(forest_2017) 

forest_tmf_1 = np.where((forest_1990 == 1) & (forest_2017 == 0),0,-9999)

forest_tmf_1 = np.where((forest_2017 == 1),1, forest_tmf_1)

forest_tmf_1  = np.where((base_map == 0), -9999, forest_tmf_1 )


mp.tif.write(base_path,
                         out_path + '/forest_tmf_2017_21_'+str(res) + 'm_repro_res.tif', 
                            forest_tmf_1, nodata=-9999, option='compress=deflate')

# deforestation in 2017 to 2021
# so pixels that had forest in 2016 but not in 2021

forest_2016 = mp.tif.read(tmf_out_path+ "\\JRC_TMF_AnnualChange_v1_Sumatra_2016_"+ res + "m_repro_res.tif", 1)
np.unique(forest_2016)
forest_2016 = np.where((forest_2016 == 1) | (forest_2016 == 2), 1, 0)
np.unique(forest_2016)


forest_2021 = mp.tif.read(tmf_out_path+ "\\JRC_TMF_AnnualChange_v1_Sumatra_2021_"+ res + "m_repro_res.tif", 1)
np.unique(forest_2021) 
forest_2021 = np.where((forest_2021 == 1) | (forest_2021 == 2), 1, 0)
np.unique(forest_2021) 


forest_tmf_2 = np.where((forest_2016 == 1) & (forest_2021 == 0), 1, -9999)
forest_tmf_2  = np.where((base_map == 0), -9999, forest_tmf_2)
np.unique(forest_tmf_2) 


mp.tif.write(base_path,
                         out_path + '/deforestation_tmf_2017_21_'+str(res) + 'm_repro_res.tif', 
                            forest_tmf_2, nodata=-9999, option='compress=deflate')


# add predictors to a string, then apply the forest in a concerted manner
# in the end ,exporting one into layer for tmf and one for Hansen


#-------#
# Slope #
#-------#
# prepare slope

# combine all elevation layers

# merged 30m layer
# list of all files
dem_30m_path = r'N:\Elevation\Sumatra\raw'
dem_30m_out_path = r'N:\Elevation\Sumatra\processed'
file_list_30_dem = glob(dem_30m_path + "\*.hgt\*.hgt")

for file_30_dem in file_list_30_dem:
    os.system("gdal_translate " + file_30_dem + " " + file_30_dem[:-4]+ ".tif") 

file_list_30_dem = glob(dem_30m_path + "\*.hgt\*.tif", recursive = True)



# so if you just run it over all files, it will only merge some tiles and not all
# so I make a smaller number of tile merges and then merge them together again
# here I am using 20, for Wallacea I used 50

chunk_size = 20
multiplier = math.ceil(len(file_list_30_dem)/chunk_size)

for i in range(0, int(multiplier)):
      print(i)
      # remember here that range always does one less than the final number, so should be -1 + 1 for the second number
      files_string = " ".join(file_list_30_dem[(chunk_size*i):min((chunk_size*(i+1)), len(file_list_30_dem))])
      # this weird coding is only becaues the file S03E118 is the chunk_sizeth and thus goes in 0,
      # but then is overwritten in i=1
      if(i==0):
          files_string = " ".join(file_list_30_dem[(chunk_size*i):min((chunk_size*(i+1)-1), len(file_list_30_dem))])
      if(i==1):
          files_string = " ".join(file_list_30_dem[(chunk_size*i)-1:min((chunk_size*(i+1)-1), len(file_list_30_dem))])
      print(files_string)
      command = workaround + "//gdal_merge.py -o " + dem_30m_out_path +"\\dem_i" + str(i) + "_30m.tif -of gtiff " + files_string
      mp.get_stdout(command)

# merge the layers
file_list_merge_30_dem = glob(dem_30m_out_path +"\\*_30m.tif")
files_string = " ".join(file_list_merge_30_dem)

files_string = " ".join(file_list_30_dem)

command = workaround +"//gdal_merge.py -o " + dem_30m_out_path +"\\dem_30m.tif -of gtiff " + files_string
mp.get_stdout(command)

mp.get_stdout("gdalinfo "+ dem_30m_out_path +"\\dem_30m.tif")

 # warp to 180 m, pay attention to reprojection
mp.get_stdout("""gdalwarp -r "bilinear" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " " + res + " -te " + extent +" "+
       dem_30m_out_path +"\\dem_30m.tif " +  out_path + "\\Sumatra_Elevation_" + res + "m_repro_res.tif" )  
        
mp.get_stdout("gdalinfo "+ out_path + "\\Sumatra_Elevation_" + res + "m_repro_res.tif" )


# then let the slope function run

# create slope map
#https://gdal.org/programs/gdaldem.html
mp.get_stdout("gdaldem slope -compute_edges -p " +  out_path + "\\Sumatra_Elevation_" + res + "m_repro_res.tif " + 
          dem_30m_out_path + "\\Sumatra_slope_" + res + "m_repro_res.tif")

# I can't see any voids here, but if there were I could fill them up with the 90 m?!
# what I would want to do is use the void-filled 90m layer to inform the value of the voids
# for now fill all no-data with 0
slope = mp.tif.read(dem_30m_out_path + "\\Sumatra_slope_" + res + "m_repro_res.tif", 1)

np.min(slope)
# 0.0

# also setting the ocean to nodata, which is not strictly necessary but less confusing

slope = np.where((base_map == 0), -9999, slope)
np.max(slope)

mp.tif.write(dem_30m_out_path + "\\Sumatra_slope_" + res + "m_repro_res.tif", 
             out_path + "\\slope_" + res + "m_repro_res.tif",
             slope,
             nodata = -9999, 
             option='compress=deflate')
                     



#---------------#
# Fire activity # 
#---------------#

processed_path = r'N:\Fire_burnt_areas\Sumatra\Fire\processed'

# prepare fire density
#DL_FIRE_M-C61_xxxx if you requested MODIS data (M6 stands for MODIS Collection 6, combined Aqua and Terra), or 
# fire_archive_M-C61_xx = MODIS standard quality Thermal Anomalies / Fire locations 
#                         processed by the University of Maryland with a 3-month
#                         lag and distributed by FIRMS. These standard data (MCD14ML) 
#                         replace the NRT (MCD14DL) files when available.
                        
                        
                        
#DL_FIRE_SV-C2_xxxx if you requested VIIRS 375m data from S-NPP
# fire_archive_SV-C2_280143
# SUOMI VIIRS C2 (S-NPP and/or NOAA-20)
# -- fire_archive_SV-C2_xx = VIIRS 375m standard Active Fire and Thermal Anomalies product
#                         processed by the University of Maryland with a 3-month lag and
#                         distributed by FIRMS. These standard data (VNP14IMGTML) replace
#                         the NRT files (VNP14IMGTDL ) when available.


# VLIIRS
# for information on attributes: https://earthdata.nasa.gov/earth-observation-data/near-real-time/firms/viirs-i-band-active-fire-data
# requested 2000-11-01 to 2022-07-04 

VIIRS_path = r'N:\Fire_burnt_areas\Sumatra\Fire\DL_FIRE_SV-C2_280143\fire_archive_SV-C2_280143.shp'

viirs_name = os.path.basename(VIIRS_path)[:-4]


mp.get_stdout("""ogr2ogr -where "CONFIDENCE = 'n' OR CONFIDENCE = 'h'" -overwrite """ + 
              processed_path +"\\" +viirs_name + "_filtered.shp " + VIIRS_path)

#1. reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
processed_path + "\\" + viirs_name + "_filtered_repro.shp " + VIIRS_path[:-4]+ "_filtered.shp")

mp.get_stdout("ogrinfo -so -al "+ processed_path + "\\" + viirs_name + "_filtered_repro.shp ")

# rasterize (with 375m resolution)
mp.get_stdout(r'C:\Anaconda3\envs\geo_py37\Library\bin\gdal_rasterize ' + 
" -burn 1 -add -l " + viirs_name + "_filtered_repro -te " + extent + "  -tr 375 -375 " + 
    processed_path + "\\" + viirs_name + "_filtered_repro.shp " + processed_path + "\\" + viirs_name + "_filtered_repro.tif") 

# Change file path to outpath here 

# convert to yearly average
viirs = mp.tif.read(processed_path + "\\" + viirs_name + "_filtered_repro.tif", 1)



# viirs runs for almost 10 years ( 9.946612) --> yearly average value
#in R
#library(foreign)
#library(lubridate)
#viirs <- read.dbf("N:\Fire_burnt_areas\Sumatra\Fire\DL_FIRE_SV-C2_280143\fire_archive_SV-C2_280143.dbf")
#as.duration(( max(viirs$ACQ_DATE) - min(viirs$ACQ_DATE)))/dyears(1) # using package lubridate

viirs = viirs / 9.946612

mp.tif.write(processed_path + "\\" + viirs_name + "_filtered_repro.tif",
                        processed_path + "\\" + viirs_name + "_yearly_average_repro.tif",
                         viirs, 
                         dtype = 4, 
                         nodata=-9999, 
                         option='compress=deflate')

# I am resampling this with near, because I don't want interpolation of values
# change resolution
mp.get_stdout("gdalwarp -r near -overwrite "+
        "-tr " + res + " -" + res + " -te " + extent +" "+
        processed_path + "\\" + viirs_name + "_yearly_average_repro.tif " +
        processed_path + "\\" + viirs_name + "_yearly_average_" + res + "m_repro_res.tif")
     



# MODIS

MODIS_path = r'N:\Fire_burnt_areas\Sumatra\Fire\DL_FIRE_M-C61_280141\fire_archive_M-C61_280141.shp'

MODIS_name = os.path.basename(MODIS_path)[:-4]


mp.get_stdout("""ogr2ogr -where "CONFIDENCE > 50" -overwrite """ + processed_path + "\\" + MODIS_name + "_filtered.shp " + MODIS_path)
# mask out for time for which we have VIIRS (20.Januar 2012)
# to not overinflate later years by doubling the fire data
# viirs and modis in sumatra spatially complimentery in a way but ...
# using min(viirs$ACQ_DATE) =  "2012-01-20"

mp.get_stdout("""ogr2ogr -where "CAST(ACQ_DATE as date) < CAST('2012/01/20' as date)" -overwrite """ +  processed_path + "\\" + MODIS_name + "_filtered_date.shp "  +  processed_path + "\\" + MODIS_name + "_filtered.shp " )
# test with 
# modis <- read.dbf("C:\Users\mv296\work\Wallacea\deforestation_model\data/Fire/include_NusaTenggara/DL_FIRE_M6_103334/fire_archive_M6_103334_filtered_date.dbf")
# max(modis$ACQ_DATE) 
# 
# as.duration(( max(modis$ACQ_DATE) - min(modis$ACQ_DATE)))/dyears(1) # using package lubridate

# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
 processed_path + "\\" + MODIS_name + "_filtered_date_repro.shp " +  processed_path + "\\" + MODIS_name + "_filtered_date.shp ")
mp.get_stdout("ogrinfo -so -al "+ processed_path + "\\" + MODIS_name + "_filtered_date_repro.shp ")


# rasterize (with 1km resolution)
mp.get_stdout(r'C:\Anaconda3\envs\geo_py37\Library\bin\gdal_rasterize ' + " -burn 1 -add -l " + MODIS_name + "_filtered_date_repro -te " + extent + " -tr 1000 -1000 " + 
processed_path + "\\" + MODIS_name + "_filtered_date_repro.shp " + processed_path + "\\" + MODIS_name + "_filtered_date_repro.tif ") 

# MODIS runs for 11.21918 years

# convert to yearly average

modis = mp.tif.read(processed_path + "\\" + MODIS_name + "_filtered_date_repro.tif ", 1)

modis = modis / 11.21918

mp.tif.write(processed_path + "\\" + MODIS_name + "_filtered_date_repro.tif ",
             processed_path + "\\" + MODIS_name + "_filtered_date_yearly_average_repro.tif",
                         modis, 
                         nodata=-9999, 
                         dtype = 4, 
                         option='compress=deflate')

# I am resampling this with near, because I don't want interpolation of values
# change resolution
mp.get_stdout("gdalwarp -r bilinear -overwrite "+
        "-tr " + res + " -" + res + " -te " + extent +" "+ 
        processed_path + "\\" + MODIS_name + "_filtered_date_yearly_average_repro.tif " + 
        processed_path + "\\" + MODIS_name + "_filtered_date_yearly_average_repro_" + res + "m_repro_res.tif")
     

# read the two layers in, combine them and then write out
viirs_repro_res =  mp.tif.read(processed_path + "\\" + viirs_name + "_yearly_average_" + res + "m_repro_res.tif", 1)

modis_repro_res = mp.tif.read(processed_path + "\\" + MODIS_name + "_filtered_date_yearly_average_repro_" + res + "m_repro_res.tif", 1)

np.max(viirs_repro_res)
np.max(modis_repro_res)

# Add rasters together and divide by 2 ignoring -9999 nodata values (-9999 not recognised as NoData in table, so is alos divided) 
fire_density = np.where((viirs_repro_res == -9999) | (modis_repro_res == -9999), -9999,(viirs_repro_res + modis_repro_res)/2)
np.max(fire_density)


fire_density = np.where((base_map == 0), -9999, fire_density)
np.max(fire_density)
 


mp.tif.write(processed_path + "\\" + MODIS_name + "_filtered_date_yearly_average_repro_" + res + "m_repro_res.tif",
                         out_path +  "\\fire_yearly_average_" + res + "m_repro_res.tif",
                         fire_density, 
                         nodata=-9999, 
                         dtype = 4, 
                         option='compress=deflate')




#-------#
# Roads #
#-------#
## roads

# we use OSM roads
# and for now we decided to use them all, because even smaller motorbikes mean access
road_processed_path = r'N:\Roads_Transport\Sumatra\processed'

road_path = r'N:\Roads_Transport\Sumatra\gis_osm_roads_free_1.shp'
mp.get_stdout("ogrinfo -so -al "+ road_path )
road_name = os.path.basename(road_path)[:-4]


# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
 road_processed_path + "\\" + road_name + "_repro.shp " +  road_path)

mp.get_stdout("ogrinfo -so -al "+ road_processed_path + "\\" + road_name + "_repro.shp ")

mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l " + road_name + "_repro -tr 30 -30 -te " + 
              extent +" "+  road_processed_path + "\\" + road_name + "_repro.shp " + 
               road_processed_path + "\\" + road_name + "_30_m_repro.tif")

# change resolution
mp.get_stdout("gdalwarp -r max -overwrite "+
        "-tr " + res + " -" + res + " -te " + extent +" "+ 
       road_processed_path + "\\" + road_name + "_30_m_repro.tif " + 
       road_processed_path + "\\" + road_name + "_" + res + "m_repro_res.tif")
     


mp.get_stdout 
# get this to run in anaconda, because it takes forever
print("""C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_proximity.py """ + 
               road_processed_path + "\\" + road_name + "_30_m_repro.tif "  +
                out_path + "\\" + road_name + "_30_m_repro_distance.tif " +
               " -nodata -9999 -distunits GEO """)




#----------#
# Rivers
#----------#
# we use OSM rivers
# and for now we decided to use them all, because even smaller motorbikes mean access
river_processed_path = r'N:\Rivers\OSM\Sumatra\processed'

river_path = r'N:\Rivers\OSM\Sumatra\gis_osm_waterways_free_1.shp'
mp.get_stdout("ogrinfo -so -al "+ river_path )
river_name = os.path.basename(river_path)[:-4]


# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
 river_processed_path + "\\" + river_name + "_repro.shp " +  river_path)

mp.get_stdout("ogrinfo -so -al "+ river_processed_path + "\\" + river_name + "_repro.shp ")

mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l " + river_name + "_repro -tr 30 -30 -te " + 
              extent +" "+  river_processed_path + "\\" + river_name + "_repro.shp " + 
               river_processed_path + "\\" + river_name + "_30_m_repro.tif")

# change resolution
mp.get_stdout("gdalwarp -r max -overwrite "+
        "-tr " + res + " -" + res + " -te " + extent +" "+ 
       river_processed_path + "\\" + river_name + "_30_m_repro.tif " + 
       river_processed_path + "\\" + river_name + "_" + res + "m_repro_res.tif")
     


mp.get_stdout 
# get this to run in anaconda, because it takes forever
print("""C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_proximity.py """ + 
               river_processed_path + "\\" + river_name + "_30_m_repro.tif "  +
                out_path + "\\" + river_name + "_30_m_repro_distance.tif " +
               " -nodata -9999 -distunits GEO """)



#---------------------------#
# Human population pressure # 
#---------------------------#

pop_pressure_path = r"N:\Accessibility_And_PopPressure\Indonesia\pressure_and_pop_density\pressurelog10_sigma"

sigma_list = [1, 2, 5, 15, 25, 50]

for sigma in sigma_list:
    print(sigma)
    mp.get_stdout("""gdalwarp -r near -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent +" "+
        pop_pressure_path + str(sigma) +  ".tif " + 
       out_path + "\\pressurelog10_sigma" + str(sigma) + "_" + res + "m_repro_res.tif")   



#---------------#
# Accessibility #
#---------------#

access_path = r'N:\Accessibility_And_PopPressure\Indonesia\travel_cost_distance_Indonesia\IDN_TTCSM_hrs.tif'
access_name = os.path.basename(access_path)[:-4]


# change resolution
mp.get_stdout("gdalwarp  -overwrite -r bilinear "+
        "-tr " + res + " -" + res + " -te " + extent +" "+ 
      access_path  + " "+
       out_path +"\\" + access_name + "_" + res + "m_repro_res.tif")

access = mp.tif.read( out_path +"\\" + access_name + "_" + res + "m_repro_res.tif", 1)
np.min(access)
np.max(access)

access = np.where((access < 0), 10000, access) # setting the NA value to a value very high

mp.tif.write(out_path +"\\" + access_name + "_" + res + "m_repro_res.tif",
                       out_path +"\\" + access_name + "_" + res + "m_repro_res_filled.tif",
                        access, 
                         nodata=-9999, 
                         dtype = 4, 
                         option='compress=deflate')

#----------#
# Land use + TORA #
#----------#



#----------#
# PIAPS #
#----------#




#--------------------------------------#
# Smallholder vs industrial concession #
# 10 m resolution 
# Source: Descals et al 2021
# https://doi.org/10.5281/zenodo.4617748 #
# layer is 
# 0 - NA / is outside of the grid (grid_withOP.shp in same folder
# 1 - industrial closed-canopy oil palm plantations
# 2 -  smallholder closed-canopy oil palm plantations
# 3 - other land covers and/or uses that are not closed-canopy oil palm.
#----------------------------------------#

# set the extent
#in_path = r''
plantation_in_path = r'C:\Users\mv296\work\world\oil_palm_and_smallholder\High_resolution_global_industrial_and_smallholder_oil_palm_map_for_2019\oil_palm_map'


with open(r'C:\Users\mv296\work\Sumatra\deforestation_model\model_input\Sumatra_cells_plantation.csv', 'r') as read_obj: # read csv file as a list of lists
  csv_reader = csv.reader(read_obj) # pass the file object to reader() to get the reader object
  list_of_rows = list(csv_reader) # Pass reader object to list() to get a list of lists

print(list_of_rows)
# unlist
rows = list(chain(*list_of_rows))


file_list = list()
#file_list.append(4)
for i in range(0, len(list_of_rows)):
  #  print("L2_2019b_"+str(rows[i]))
    file_list.append(plantation_in_path+"\L2_2019b_"+str(rows[i])+".tif")

print(file_list)
# split the list because as a whole seemingly too large
file_list_a = file_list[0:round(len(file_list)/2)]
file_list_b = file_list[round(len(file_list)/2)+1:len(file_list)]


files_string_a = " ".join(file_list_a)
command = "C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_merge.py -o " + out_path +"\\plantations_a.tif -of gtiff " + files_string_a
mp.get_stdout(command)


files_string_b = " ".join(file_list_b)
command = "C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_merge.py -o " + out_path +"\\plantations_b.tif -of gtiff " + files_string_b
mp.get_stdout(command)

# join the two
command = "C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_merge.py -o " + out_path +"\\plantations.tif -of gtiff " + out_path +"\\plantations_a.tif " + out_path +"\\plantations_b.tif "
mp.get_stdout(command)



# then reproject 

# code 0 as NA!!!!


# and use extent of base/forest



#------------------------#
# check infrastructure
# check peat
# Check mining, 
# check main commodity (!!!!!)
# check transmigrant 
#------------------------#



# clip with forest for all layers





##################################################
##################################################

# Appendix


# I first thought to use a combo of forest now and first year of defor
# but I keep deciding that this is not a good approach because i dont know if there was 
# regrowth in th emeantime

# use a comination of undisturbed-degraded tropical moist forest + deforestation year (which is the first year a pixel was deforested)

# 
TMF_path = r'N:\Landcover\Sumatra\JRC_TMF'

file_list = glob(TMF_path + '//JRC_TMF_UndisturbedDegradedForest_v1_1982_2021_ASI_ID*.tif')

files_string = " ".join(file_list)

command = workaround + "\gdal_merge.py -o "+ TMF_path+ "JRC_TMF_UndisturbedDegradedForest_v1_1982_2021.tif" + " -of gtiff " + files_string

mp.get_stdout(command)

# reproject

mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent + " "+
         TMF_path+ "JRC_TMF_UndisturbedDegradedForest_v1_1982_2021.tif" + " " +TMF_path+ "JRC_TMF_UndisturbedDegradedForest_v1_1982_2021_" + res + "m_repro_res.tif" )  
        


file_list = glob(TMF_path + '//JRC_TMF_DeforestationYear_INT_*.tif')

files_string = " ".join(file_list)

command = workaround + "\gdal_merge.py -o "+ TMF_path+ "JRC_TMF_DeforestationYear.tif" + " -of gtiff " + files_string

mp.get_stdout(command)


# reproject

mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent + " "+
         TMF_path+ "JRC_TMF_DeforestationYear.tif" + " " +TMF_path+ "JRC_TMF_DeforestationYear_" + res + "m_repro_res.tif" )  
        
