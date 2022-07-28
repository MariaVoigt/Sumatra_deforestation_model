# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:46:06 2022

@author: mv296
"""


import numpy as np
import os
from osgeo import ogr
from glob import glob
import csv
from itertools import chain

import macpyver as mp

# unit conversion
km2_ha = 100
ha_km2 = 0.01
km2_Mha = 0.0001
m2_km2 = 0.000001

# resolution (180m)
res = str(180)

out_path = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res'

#-----------#
# Baseshape #
#-----------#

# process the base shape
# select the regions by name
# Aceh, Sumatera Utara, Riau, Sumatera Barat, Jambi, Bengkulu, Sumatera Selatan, 
# Bangka Belitung, Lampung, Kepulauan Riau

# check how much forest in each

# create a wallacea_shape file
# WE ARE CUTTING AT LYDEKKERS LINE EXCLUDING ARU SEE EG https://www.starfish.ch/dive/Wallacea.html
unit_filter = str("NAME_1 LIKE '%Aceh%' OR NAME_1 LIKE '%Sumatera%'  OR NAME_1 LIKE '%Jambi%' OR NAME_1 LIKE '%Bengkulu%' OR NAME_1 LIKE '%Bangka Belitung%' OR NAME_1 LIKE '%Lampung%'  OR NAME_1 LIKE '%Riau%' AND NAME_1 NOT LIKE '%Kepulauan Riau%'")

# i am excluding  OR NAME_1 LIKE '%Kepulauan Riau%' because it includes islands in the West of Sarawak and is just a bit messy

gadm_path = r'C:\Users\mv296\work\Indonesia\gadm\gadm36_IDN_1\gadm36_IDN_1.shp'              
sumatra_shape_path = r'C:\Users\mv296\work\Sumatra\data\gadm\sumatra_complete_shape.shp'            
mp.get_stdout("""ogr2ogr  -overwrite -a_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -s_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -sql "SELECT * FROM  gadm36_IDN_1 WHERE """+ unit_filter  + """ " """ + sumatra_shape_path + " " + gadm_path)

sumatra_basename = os.path.basename(sumatra_shape_path)[:-4]

# reproject
os.system("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + 
    sumatra_shape_path[:-4] + "_repro.shp" + " " + sumatra_shape_path)
# the file sumatra_shape_path[:-4] + "_repro.shp" will be used for regions!

mp.get_stdout("ogrinfo -so -al " +  sumatra_shape_path[:-4] + "_repro.shp" )

# extract extent:
extent = "-3292738.004418 962243.639245 -1721896.360977 2256489.384612"

mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l "+ sumatra_basename +"_repro -tr " + res + " -" + res + " -te " + extent +" "+  sumatra_shape_path[:-4] + "_repro.shp " +  out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")

mp.get_stdout("gdalinfo " +  out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")

extent = "-3292738.004 962289.385 -1721878.004 2256489.385"
# save the path here
base_path = out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif"

#--------#
# Forest #
#--------#





#-------#
# Slope #
#-------#



#---------------#
# Fire activity # 
#---------------#


#---------------------------#
# Human population pressure # 
#---------------------------#

#---------------#
# Accessibility #
#---------------#


#----------#
# Land use #
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
# and use extent of base/forest



#------------------------#
# check infrastructure
# check peat
# Check mining, 
# check main commodity
# check transmigrant 
#------------------------#



# clip with forest for all layers
