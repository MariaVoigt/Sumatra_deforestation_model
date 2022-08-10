# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:46:06 2022

@author: mv296


Script to prepare the predictor layers for Sumatra
# use the raw layers from the NAS drive, process them and them export them into model_input/repro_res
# in the next scrpt, will work through them and prep them with
# Hansen or TMF layer

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

out_path_forest = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\forest'
out_path_base = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\base'
out_path_pred = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\predictors'


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
gadm_path = r'N:\Admin_boundaries\GADM\Asia\gadm36_IDN_shp\gadm36_IDN_1.shp'              
sumatra_shape_path = r'N:\Admin_boundaries\GADM\Sumatra\sumatra_complete_shape.shp'        
    
mp.get_stdout("""ogr2ogr  -overwrite -a_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -s_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -sql "SELECT * FROM  gadm36_IDN_1 WHERE """+ unit_filter  + """ " """ + sumatra_shape_path + " " + gadm_path)

sumatra_basename = os.path.basename(sumatra_shape_path)[:-4]

# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + 
    sumatra_shape_path[:-4] + "_repro.shp" + " " + sumatra_shape_path)
# the file sumatra_shape_path[:-4] + "_repro.shp" will be used for regions!

mp.get_stdout("ogrinfo -so -al " +  sumatra_shape_path[:-4] + "_repro.shp" )

# extract extent:
pre_extent = "-3292738.004418 962243.639245 -1721896.360977 2256489.384612"

mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l "+ sumatra_basename +"_repro -tr " + res + " -" + res + " -te " + pre_extent +" "+  sumatra_shape_path[:-4] + "_repro.shp " +  out_path_base + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")

mp.get_stdout("gdalinfo " +  out_path_base + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")

extent = "-3292738.004 962289.385 -1721878.004 2256489.385"


# save the path here
base_path =out_path_base + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif"

base_map = mp.tif.read(base_path, 1)

# outer shape (digitized in ArcGIS for clipping and others)

outer_base_shape = r'N:\Admin_boundaries\GADM\Sumatra\Sumatra_outer_shape.shp'


#--------#
# Forest #
#--------#

# primary forest
# PRIMARY FOREST LAYER
# reproject primary forest layer from MARGONO / Hansen to clip the forest with
primary_processed_path = r'N:\Landcover\Indonesia\Margono_primary_forest_change\processed'


mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent +
        r" N:\Landcover\Indonesia\Margono_primary_forest_change\timeseq_change00_12.tif"  + " " + primary_processed_path + "//timeseq_change00_12_" + res + "m_repro_res.tif" )  
        
# recode
primary = mp.tif.read(primary_processed_path + "\\timeseq_change00_12_" + res + "m_repro_res.tif", 1)
# only category 3 is non-forets because all the other change will be captured within the change  
# from hansen
primary = np.where((primary == 3) | (primary ==0), 0, 1)      

mp.tif.write(primary_processed_path + "\\timeseq_change00_12_" + res + "m_repro_res.tif",
                          primary_processed_path +"\\primary_forest_"+res+"m_repro_res.tif",
                       primary, nodata = -9999, option='compress=deflate')
                       

# continue checking this

# mangrove forest from mangrove atlas

# I used this layer C:\Users\mv296\work\Wallacea\deforestation_model\data\Mangroves\GMW_001_GlobalMangroveWatch\GMW_001_GlobalMangroveWatch\01_Data\GMW_2016_v2.shp
# and the 
# to clip in ArcGIS and produce the Sumatra layer
mangrove_path = r'C:\Users\mv296\work\Sumatra\data\forest\Mangroves_GMW_2016_v2_Sumatra.shp'
mangrove_processed_path = r'N:\Mangroves\GMW_001_GlobalMangroveWatch\processed'
# reproject 
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
mangrove_path[:-4] + "_repro.shp "  + mangrove_path )

# then rasterize
mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l Mangroves_GMW_2016_v2_Sumatra_repro -tr " + res + " -" + res + " -te " + extent +" "+ mangrove_path[:-4] + "_repro.shp " +  mangrove_processed_path + "\\" + "GMW_2016_v2_"+ res + "_m_repro.tif")


# combine mangrove and primary forest and Hansen layer// what about David



mangrove = mp.tif.read(mangrove_processed_path + "\\" + "GMW_2016_v2_"+ res + "_m_repro.tif", 1)
primary = np.where((primary == 1), 1, -9999)




# forest_loss
Hansen_forest_processed_path = r'N:\Landcover\Sumatra\Hansen_GFC-2021-v1.9\processed'

forest_loss_path = r'N:\Landcover\Indonesia\Forest_Hansen\v1.9\Hansen_GFC-2021-v1.9_lossyear.tif'

file_list = glob(r'N:\Landcover\Indonesia\Forest_Hansen\v1.9\Hansen_GFC-2021-v1.9_lossyear*.tif')

files_string = " ".join(file_list)

command = workaround + "\gdal_merge.py -o "+ forest_loss_path +" -of gtiff " + files_string

mp.get_stdout(command)


mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
              "-tr " + res + " -" + res + " -te " + extent + " "+
              forest_loss_path + " " +  Hansen_forest_processed_path + "\\Hansen_GFC-2021-v1.9_lossyear_"+ res + "m_repro_res.tif" )  
        



# forest cover

# forest cover
# combining all of Indonesia, just because its easier
forest_cover_path = r'N:\Landcover\Indonesia\Forest_Hansen\v1.9\Hansen_GFC-2021-v1.9_treecover2000.tif'

file_list = glob(r'N:\Landcover\Indonesia\Forest_Hansen\v1.9\Hansen_GFC-2021-v1.9_treecover2000_*.tif')

files_string = " ".join(file_list)

command = workaround + "\gdal_merge.py -o "+ forest_cover_path +" -of gtiff " + files_string

mp.get_stdout(command)


mp.get_stdout("""gdalwarp -r "bilinear" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
              "-tr " + res + " -" + res + " -te " + extent + " "+
              forest_cover_path + " " +  Hansen_forest_processed_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif" )  
        

forest_2000 = mp.tif.read(Hansen_forest_processed_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif", 1)


# then fix cover > 70% 
forest_2000_70 = np.where(forest_2000 >= 70, 1, -9999) 

np.unique(forest_2000_70)


# mp.tif.write(out_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif", 
#              out_path +"\\tree_cover_repro_70.tif",
#                        forest_2000_70,  
#                        nodata = -9999,
#                        option='compress=deflate')



forest_2000_clip = np.where(primary ==1, forest_2000_70, -9999)
# add mangrove forest in
forest_2000_clip_mangrove = np.where(mangrove == 1, 1, forest_2000_clip)
 
mp.tif.write(Hansen_forest_processed_path + "\\Hansen_GFC-2021-v1.9_treecover2000_"+ res + "m_repro_res.tif", 
             Hansen_forest_processed_path + "\\tree_cover_repro_70_clip_" + res + "m_repro_res.tif", 
             forest_2000_clip_mangrove, 
             nodata = -9999,
             option='compress=deflate')




# prepare forest loss layers for model

forest_2000 = mp.tif.read(Hansen_forest_processed_path + "\\tree_cover_repro_70_clip_" + res + "m_repro_res.tif", 1)
forest_loss = mp.tif.read(Hansen_forest_processed_path + "\\Hansen_GFC-2021-v1.9_lossyear_"+ res + "m_repro_res.tif" , 1)
np.unique(forest_loss)
# test 
base_map.shape == forest_2000.shape == forest_loss.shape



 # prepare two layers:
    # one is:  0 = previous deforestation, 1= forest, -9999 = sea, non-forest vegetation
    #  one is: -9999 = not deforested or not considered, 1 = forest converted to non-forest 

forest_1 = np.where((forest_2000 == 0), -9999, forest_2000)
    # this is previous deforestation
forest_1 = np.where((forest_loss >= 1) & (forest_loss <= 16) & (forest_1 != -9999), 0, forest_1)
    # this is current forest
forest_1 = np.where((forest_loss >= 17) & (forest_2000 == 1), 1, forest_1)
forest_1 = np.where((base_map == 0), -9999, forest_1)

mp.tif.write(base_path,
                         out_path_forest + '/forest_hansen_2017_21_'+str(res) + 'm_repro_res.tif', 
                            forest_1, nodata=-9999, option='compress=deflate')
                         
    # layer 2 (deforestation in 2016-2018)
forest_2 = np.where((base_map == 0)|(forest_2000 == 0), -9999, forest_loss)
forest_2 = np.where((forest_2 >= 17) & (forest_2000 == 1), 1, -9999)
forest_2 = np.where((base_map == 0), -9999, forest_2)

mp.tif.write(base_path,
                             out_path_forest + '/deforestation_hansen_2017_21_'+str(res) + 'm_repro_res.tif',
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

tmf_out_path = r'N:\Landcover\Sumatra\TMF_JRC_forest\processed'


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
                         out_path_forest  + '/forest_tmf_2017_21_'+str(res) + 'm_repro_res.tif', 
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
                         out_path_forest + '/deforestation_tmf_2017_21_'+str(res) + 'm_repro_res.tif', 
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
       dem_30m_out_path +"\\dem_30m.tif " +  out_path_pred + "\\Sumatra_elevation_" + res + "m_repro_res.tif" )  
        
mp.get_stdout("gdalinfo "+ dem_30m_out_path + "\\Sumatra_elevation_" + res + "m_repro_res.tif" )


# then let the slope function run

# create slope map
#https://gdal.org/programs/gdaldem.html
mp.get_stdout("gdaldem slope -compute_edges -p " +  dem_30m_out_path + "\\Sumatra_elevation_" + res + "m_repro_res.tif " + 
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
             out_path_pred + "\\slope_" + res + "m_repro_res.tif",
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
# modis <- read.dbf()
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
                         out_path_pred +  "\\fire_yearly_average_" + res + "m_repro_res.tif",
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
     

mp.get_stdout("""C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_proximity.py """ + 
              road_processed_path + "\\" + road_name + "_" + res + "m_repro_res.tif "  +
               out_path_pred + "\\road_" + res + "m_repro_distance.tif " +
               " -nodata -9999 -distunits GEO """)


# reproject to 180 m

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
     


mp.get_stdout("""C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_proximity.py """ + 
               river_processed_path + "\\" + river_name + "_" + res + "m_repro_res.tif "  +
               out_path_pred + "\\river_" + res + "m_repro_distance.tif " +
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
       out_path_pred + "\\pressurelog10_sigma" + str(sigma) + "_" + res + "m_repro_res.tif")   



#---------------#
# Accessibility #
#---------------#

access_path = r'N:\Accessibility_And_PopPressure\Indonesia\travel_cost_distance_Indonesia\IDN_TTCSM_hrs.tif'
access_name = os.path.basename(access_path)[:-4]


# change resolution
mp.get_stdout("gdalwarp  -overwrite -r bilinear "+
        "-tr " + res + " -" + res + " -te " + extent +" "+ 
      access_path  + " "+
       access_path[:-17] + access_name + "_" + res + "m_repro_res.tif")

access = mp.tif.read(access_path[:-17] + access_name + "_" + res + "m_repro_res.tif", 1)
np.min(access)
np.max(access)

access = np.where((access < 0), 10000, access) # setting the NA value to a value very high

mp.tif.write(access_path[:-17] + access_name + "_" + res + "m_repro_res.tif",
                     out_path_pred +"\\" + access_name + "_" + res + "m_repro_res_filled.tif",
                        access, 
                         nodata=-9999, 
                         dtype = 4, 
                         option='compress=deflate')

#----------#
# Land use #
# we are using hte 2010 layer, although that is outdated, 
# but I need to get the layer on the road
#----------#

infile = r'N:\Landuse\Indonesia\Indonesia_legal_classification'
    
# reproject
  
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + 
    infile + "\\processed\\Indonesia_legal_classification_repro.shp " + infile + "\\Indonesia_legal_classification.shp")

mp.get_stdout("ogrinfo  -so -al "+   infile + "\\processed\\Indonesia_legal_classification_repro.shp ")

"""    
we will have classes (from state of forests)
#1 - Non Forestland - Non-Protected Areas (APL) --> can be converted to agriculture, settlement etc. ,
#2 - permanent production (HP) --> clear cutting and timber plantation,
#2 - convertible production forest (HPK) --> clear cutting and ind plantation (or released to non-forest land),
#3 - limited production forest (HTP) --> logging,
#4 - protection forest (hutan lindung) -> protect buffer for water systems, flood prevention, erosion protection, etc.,
#5 - conservation forest (hutan konservasi) -> particular characteristic and main function to protect biodiversity and ecosystem

#4 and 5 will be lumped in image and model
9 - other (water?) --> -9999
-> first code 0-4 +1 because of the already existing nodata 0
"""
#os.system("ogrinfo  "+ infile[:-4] + "_repro.shp" + """ -sql "ALTER TABLE Indonesia_legal_classification_repro DROP COLUMN st_areasha" """)

mp.get_stdout("ogrinfo  "+   infile + "\\processed\\Indonesia_legal_classification_repro.shp " + """ -sql "ALTER TABLE Indonesia_legal_classification_repro DROP COLUMN LU_id" """)



source = ogr.Open( infile + "\\processed\\Indonesia_legal_classification_repro.shp" , update=True)
layer = source.GetLayer()
layer_defn = layer.GetLayerDefn()
field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]

feature_count = layer.GetFeatureCount()
feature_count = list(range(1, feature_count + 1) )

# Add a new field

field_defn = ogr.FieldDefn( "LU_id", ogr.OFTReal )
layer.CreateField(field_defn)

for i in layer:
    kh_fungsi = i.GetField("kh_fungsi_")
    if kh_fungsi == "APL":
        i.SetField( "LU_id", 1) 
     
    if kh_fungsi == "HP":
        i.SetField( "LU_id", 2)
    if kh_fungsi == "HPK":
        i.SetField( "LU_id", 2) 
    if kh_fungsi == "HPT":
        i.SetField( "LU_id", 3)  
    if kh_fungsi == "HL":
        i.SetField( "LU_id", 4)
    if kh_fungsi  == "CA" or kh_fungsi == "HSAW" or kh_fungsi == "KSPA" or kh_fungsi == "SM" or kh_fungsi == "TN" or kh_fungsi == "TAHURA" or kh_fungsi == "TNL" or kh_fungsi == "TWA" or kh_fungsi == "TWA/HW" or kh_fungsi == "TWAL" or kh_fungsi == "TB":
        i.SetField( "LU_id", 5)
    if kh_fungsi != "HL" and kh_fungsi  != "CA" and kh_fungsi != "HSAW" and kh_fungsi != "KSPA" and kh_fungsi != "SM" and kh_fungsi != "TN" and kh_fungsi != "TAHURA" and kh_fungsi != "TNL" and kh_fungsi != "TWA" and kh_fungsi != "TWA/HW" and kh_fungsi != "TWAL" and kh_fungsi != "TB" and kh_fungsi != "APL" and kh_fungsi != "HP" and kh_fungsi != "HPK" and kh_fungsi != "HPT":
        i.SetField( "LU_id", 9)   
    layer.SetFeature(i)
    
source = None       

# rasterize burning LU_id
os.system("""gdal_rasterize -a_nodata 0 -a "LU_id" -l Indonesia_legal_classification_repro -tr """ + res + " -" + res +
" -te " + extent +" "+ infile + "\\processed\\Indonesia_legal_classification_repro.shp " + 
infile + "\\processed\\Indonesia_legal_classification_" + res + "m_repro_res.tif")

# recode
lu = mp.tif.read(infile + "\\processed\\Indonesia_legal_classification_" + res + "m_repro_res.tif", 1)
# 0 becomes -9999, 0 becomes -9999, rest goes down 1
lu = np.where((lu == 0) | (lu == 9) | (lu == -9999), -9999, lu-1)


# clip landuse with forest extent (and island extent included)
# I will also lump HL and CA --> 3 and 4
np.unique(lu)
lu =  np.where(lu == 4, 3, lu)
lu =np.where(lu >= 0, lu + 1, lu)
np.unique(lu)
# this is the reference class
lu = np.where(lu == 4, 0, lu)


mp.tif.write(infile + "\\processed\\Indonesia_legal_classification_" + res + "m_repro_res.tif", 
             out_path_pred +  "\\lu_" + res + "m_repro_res.tif",
                       lu,
                       nodata = -9999, 
                       option='compress=deflate')



# final layer classes
# 0 -PA + HL 
# 1 - APL
# 2 - Production forest
# 3 - Limited production forestt 

#----------#
#PIAPS #
# social forestry started in xxx
# pray that I processed
# that indonesia wide
# the two layers are almost everywhere complimentary
# so I am lumping them for now
# no social forestry
# proposed
# implemented
# use script C:/Users/mv296/work/Sumatra/deforestation_model/src/scripts/PIAPS_2020_processing_shapefiles.R

#----------#

piaps_path = r"N:/Landuse/Indonesia/PIAPS/PIAPS_Sept2020_Kraus/PIAPS_Sept2020_short_selected/piaps_combined.shp"

# implemented - 1
# proposed/available 2

# reproject
mp.get_stdout("ogrinfo -so -al "+ piaps_path )
piaps_name = os.path.basename(piaps_path)[:-4]


# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
 piaps_path[:-4] + "_repro.shp " +  piaps_path)

mp.get_stdout("ogrinfo -so -al "+piaps_path[:-4] + "_repro.shp ")




# rasterise

mp.get_stdout(""" gdal_rasterize -a_nodata 0 -a status_id -l """ + piaps_name + "_repro " +
                    "-tr " + res + " -" + res + " -te " + extent +" " + 
                        piaps_path[:-4] + "_repro.shp " +
                         piaps_path[:-4] + "_repro.tif")


#----------#
# I thought about TORA
# but it has only been implemented 2020/21
# so not really a pattern there yet
#-------------

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

plantation_in_path = r'N:\Landcover\world\oil_palm_and_smallholder\High_resolution_global_industrial_and_smallholder_oil_palm_map_for_2019\oil_palm_map'
plantation_processed_path =  r'N:\Landcover\world\oil_palm_and_smallholder\High_resolution_global_industrial_and_smallholder_oil_palm_map_for_2019\processed'

with open(plantation_processed_path + '\\Sumatra_cells_plantation.csv', 'r') as read_obj: # read csv file as a list of lists
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
file_list_b = file_list[round(len(file_list)/2):len(file_list)+1]





files_string_a = " ".join(file_list_a)
command_a = "C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_merge.py -o " + plantation_processed_path  +"\\plantations_a.tif -of gtiff " + files_string_a
mp.get_stdout
print(command_a)


files_string_b = " ".join(file_list_b)
command_b = "C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_merge.py -o " + plantation_processed_path +"\\plantations_b.tif -of gtiff " + files_string_b
mp.get_stdout
print(command_b)


# join the two
# the n 0 is important here, otherwise it will write over a few tiles wiht 0 info
command = "C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_merge.py -n 0 -o " + plantation_processed_path +"\\plantations.tif -of gtiff " +plantation_processed_path +"\\plantations_a.tif " + plantation_processed_path +"\\plantations_b.tif "
mp.get_stdout
print(command)



# then reproject 
mp.get_stdout("gdalinfo "+ plantation_processed_path  +"\\plantations.tif")
mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent + " "+
        plantation_processed_path +"\\plantations.tif "+
        plantation_processed_path +"\\plantations_" + res + "m_repro_res.tif")

# code 0 as NA!!!!
plantation = mp.tif.read(plantation_processed_path +"\\plantations_" + res + "m_repro_res.tif", 1)

np.unique(plantation )


# recode 
plantation  = np.where((plantation  == 0)|(plantation == 3),-9999, plantation)

mp.tif.write(plantation_processed_path +"\\plantations_" + res + "m_repro_res.tif", 
            out_path_pred  +"\\plantations_" + res + "m_repro_res.tif",
                       plantation,
                       nodata = -9999,
                       option='compress=deflate')



#------------------------#
# Peat
# checked a layer from the government
# both available through menlhk and from ministry of agriculture
# but decided going with CRIS{}
#---------------------# 


peat_path = r'N:\Landcover\Indonesia\CRISP\CRISP_SEA_2015_deliverable'
peat_name = 'Per-humid_SEA_LC_2015_CRISP_Geotiff_indexed_colour.tif'

# reproject first
mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent + " "+
        peat_path  + "\\" + peat_name + " " +
        peat_path  + "\\processed\\" + peat_name[:-4] +"_" + res + "m_repro_res.tif")
        
# import


peat = mp.tif.read(peat_path  + "\\processed\\" + peat_name[:-4] + "_repro.tif", 1)

np.unique(peat)
# we are looking at 3, which is peatforest

# recode 
peat = np.where((peat == 3),1, 0)



mp.tif.write( peat_path  + "\\processed\\" + peat_name[:-4] + "_repro.tif", 
             out_path_pred +  "\\peat_" + res + "m_repro_res.tif",
                       peat,
                       nodata = -9999,
                       option='compress=deflate')


# export


#---------#
# Mining
#----------#
mining_path  = r'N:\Landuse\Indonesia\WRI\Mining\processed\WRI_IUP_tambang_repro.shp'
# this path is already reprojected and has the type info in column 'type'
mp.get_stdout("ogrinfo  "+ mining_path + """ -sql "ALTER TABLE WRI_IUP_tambang_repro DROP COLUMN type_1 " """)


#reproject
mp.get_stdout("ogrinfo -so -al "+ mining_path  )
# no need to reproject

# category 1
# Eksplorasi   // Ekplorasi //  Eskplorasi  
# Pencadangan Wilayah - regional backup
# Studi Kelayakan - feasibility study
# WIUPK:  Wiupk (Special Mining Business Licence Area (Wilayah Izin Usaha Pertambangan Khusus – “WIUPK”) means an area that is authorised to a Special Mining Business Licence holder.)
        # An IUPK will be granted after the mining company has secured a WIUPK. IUPK holders are permitted to carry out mining activities only in the WIUPK. S

#category 2
# Konstruksi
# Eksploitasi
# Operasi Prodduksi //Operasi Produiksi    
# Eksplorasi operasi Produksi 


# rasterise

mp.get_stdout(""" gdal_rasterize -a type -l WRI_IUP_tambang_repro """ +
                    "-tr " + res + " -" + res + " -te " + extent +" " + 
                        mining_path + " "+
                        out_path_pred + "\\mining_" + res + "m_repro_res.tif")

#------------------------#

# main commodity  from PODES
# podes main commodity is first processed in  
# src/scripts/prepare_podes_livelihood.R
#------------------------#

podes_in_path = r'C:\Users\mv296\work\Sumatra\data\model_input\PODES_processed/Sumatra_PODES2018_livelihood.shp'
podes_processed_path = r'C:\Users\mv296\work\Sumatra\data\model_input\PODES_processed'
# reproject
mp.get_stdout("ogrinfo -so -al "+ podes_in_path )
podes_name = os.path.basename(podes_in_path)[:-4]


# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" """ + 
podes_in_path[:-4] + "_repro.shp " +  podes_in_path)

mp.get_stdout("ogrinfo -so -al "+podes_in_path[:-4] + "_repro.shp ")


"""

Lookup
1 Rice
2 Palawija (corn, beans, sweet potatoes)
3 Horticulture (fruits, vegetables, ornamental plants, medicinal plants, etc.)
4 Rubber
5 Palm oil
6 Coffee
7 Cocoa
8 Coconut
9 Pepper
10 Cloves
11Tobacco
12 Sugar cane
13 Animal Husbandry (cattle, sheep, chickens, etc.)
14 Capture fisheries (including other biota)
15 Aquaculture (including other biota)
16 Forestry cultivation (teak, mahogany, sengon, bamboo, etc.)
17 Collection of forest products (honey, agarwood, fruits, firewood, etc.)
18 Capture of wild animals (pigs, partridges, deer, etc.)
19 Captivity of wild animals / plants (arowana, crocodile, orchid, etc.)
20 Agricultural services (hatchery, tractor rental, rattan, etc.)
21 Other

subsistence agriculture: Rice, Palawija, Coffee, Cocoa, Coconut, Pepper, Cloves, 
Tobacco and Sugar cane)

Plantation agriculture: Rubber, Palm oil, and forest culivataion as main commoditiy were 
 grouped in plantation category. 
 
Non agriculture: All other main commodities (horticulture, animal husbandry, capture fisheries,
 aquaculture, collection of forest products, capture of wild animals, 
 captivity of wild animals and agricultural services) 
 
 # then all non-forest pixels with these categories get filtered out and all
 # forest pixels get a distance to these non-forest pixels (polygonize and then 
 distance function)
 
#
to figure out whether rubber is subsistence or plantation
# we extract the majority of landuse in a desa
#
"""




mp.get_stdout(""" gdal_rasterize -a R403B -l Sumatra_PODES2018_livelihood_repro """ +
                    "-tr " + res + " -" + res + " -te " + extent +" " + 
                       podes_in_path[:-4] + "_repro.shp "+
                        podes_in_path[:-4] + "_" + res + "m_repro_res.tif")


# lookup
# see R script plot_PODES_livelihood_types

subsistence = [1, 2, 6, 7, 8, 9, 10, 11, 12]
plantations =  [4, 5, 16]
non_agri = [3, 13, 17, 18, 19, 20, 21, 14, 15] 


# rasterize the layer and then np.where the aboe classes and

livelihood = mp.tif.read(podes_in_path[:-4] + "_" + res + "m_repro_res.tif", 1)

# make layers with each type

livelihood_recode = np.where(np.isin(livelihood,subsistence), 100, livelihood )
livelihood_recode = np.where(np.isin(livelihood_recode, plantations), 200, livelihood_recode )
livelihood_recode = np.where(np.isin(livelihood_recode, non_agri), 300, livelihood_recode )


np.unique(livelihood_recode)


mp.tif.write(  podes_in_path[:-4] + "_" + res + "m_repro_res.tif", 
            podes_processed_path+"\\livelihood_recode_" + res + "m_repro_res.tif",
                       livelihood_recode,
                       nodata = 0,
                       option='compress=deflate')


# additionally make layer with different livelihood classes that are categorical


# go through subsistence, plantation, non_agri
# we want the distance of non-forest pixels with a certain livelihood
# for each forested pixel 

livelihood_class = ["subsistence", "plantation", "non_agri"]
livelihood_class_value = [100, 200, 300,]

for i in range(0,len(livelihood_class)):
    print(i)
    livelihood_layer= np.where(livelihood_recode == livelihood_class_value[i], 1, -9999)
    mp.tif.write(podes_processed_path+"\\livelihood_recode_" + res + "m_repro_res.tif",
                          podes_processed_path + "\\" + livelihood_class[i]+
                            "_" + res + "m_repro_res.tif",
                            livelihood_layer, dtype = 1, 
                            nodata = -9999, option='compress=deflate') 
    # gdal proximity.py
    mp.get_stdout("""C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_proximity.py -values 1 """ + 
    podes_processed_path + "\\" + livelihood_class[i] + "_" + res + "m_repro_res.tif " +
    out_path_pred + "\\" + livelihood_class[i] +"_" + res + "m_repro_res_distance.tif" ) 



#----------#
# Transmigrant 
#------------------------#

# WORK WITH THE LANDCOVER FROM MEF TO FILTER OUT TRANSMIGRATION
# reproject and bring it in, 
# filter out category 20122
# distance? 
#----------------------#

LC_transmigrant_path = r'N:\Landcover\Sumatra\MoF_2019'
LC_transmigrant_processed_path = r'N:\Landcover\Sumatra\MoF_2019\processed\land_cover_2019.shp'


file_list = glob(LC_transmigrant_path + "\\**\*.shp", recursive = True)

# append each file
mp.get_stdout("""ogr2ogr -f "Esri Shapefile" """ + LC_transmigrant_processed_path + " " + file_list[0])

#Then merge the following files by using:
for file_path in file_list[1:]:
    print(file_path)
    mp.get_stdout("""ogr2ogr -f "ESRI Shapefile" -update -append """ + LC_transmigrant_processed_path + " " +  file_path + " -nln land_cover_2019")



# reproject 
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + 
    LC_transmigrant_processed_path[:-4] + "_repro.shp" + " " + LC_transmigrant_processed_path)
    
mp.get_stdout("ogrinfo -so -al "+  LC_transmigrant_processed_path[:-4] + "_repro.shp")
# raster

# it can have information on transmigration in either LC15 or LC96, so rasterise both and then pu ttogether
# eh need to have codes instead


mp.get_stdout("gdal_rasterize -a_nodata 0 -a PL_19_R -l land_cover_2019_repro -tr " + res + " -" + res + " -te " + extent +" "+ LC_transmigrant_processed_path[:-4] + "_repro.shp " + LC_transmigrant_processed_path[:-4] + "_"+ res + "_m_repro.tif")



        
LC_transmigrant = mp.tif.read(LC_transmigrant_processed_path[:-4] + "_"+ res + "_m_repro.tif", 1)
np.unique(LC_transmigrant)
LC_transmigrant = np.where(LC_transmigrant == 20122, 1, -9999)
np.unique(LC_transmigrant)


# use the distance to non-forested AND forested transmigrant pixels
    
mp.tif.write(LC_transmigrant_processed_path[:-4] + "_"+ res + "_m_repro.tif",
                          LC_transmigrant_processed_path[:-4] + "_recode_"+ res + "_m_repro.tif",
                         LC_transmigrant, nodata=-9999, dtype = 0, option='compress=deflate')
                         
mp.get_stdout("C:\Anaconda3\envs\geo_py37\python.exe C:\Anaconda3\envs\geo_py37\Scripts\gdal_proximity.py -values 1 " + 
    LC_transmigrant_processed_path[:-4] + "_recode_"+ res + "_m_repro.tif " +
    out_path_pred + "//transmigrant_distance_"+ res + "_m_repro.tif") 



# THE END
# continue in script C:\Users\mv296\work\Sumatra\deforestation_model\analysis\src\scripts\prepare_layers_defor_model_sumatra_units_b.py 
