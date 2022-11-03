# -*- coding: utf-8 -*-
"""
Created on Mon Dec 02 17:31:48 2019

@author: mv39zilo


"""
import numpy as np
import os
from osgeo import ogr
import shutil
from glob import glob
import macpyver as mp
res = str(180)
# unit area


# reproject area according to predicted deforestation, then I can use this to calculate
# the percent as well 
forest_layer = 'hansen'
out_path = r"C:\Users\mv296\work\Sumatra\deforestation_model\results" +"\\"+ forest_layer  +"\\"+"deforestation_quantification"




unit_filter = str("NAME_1 LIKE '%Aceh%' OR NAME_1 LIKE '%Sumatera%'  OR NAME_1 LIKE '%Jambi%' OR NAME_1 LIKE '%Bengkulu%' OR NAME_1 LIKE '%Lampung%'  OR NAME_1 LIKE '%Riau%' AND NAME_1 NOT LIKE '%Kepulauan Riau%'")

# i am excluding  OR NAME_1 LIKE '%Kepulauan Riau%' because it includes islands in the West of Sarawak and is just a bit messy
gadm_path = r'N:\Admin_boundaries\GADM\Asia\gadm41\gadm41_IDN_shp\gadm41_IDN_1.shp'              
sumatra_shape_path = r'N:\Admin_boundaries\GADM\Sumatra\sumatra_complete_shape_no_BB.shp'        
    
mp.get_stdout("""ogr2ogr  -overwrite -a_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -s_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -sql "SELECT * FROM  gadm41_IDN_1 WHERE """+ unit_filter  + """ " """ + sumatra_shape_path + " " + gadm_path)

sumatra_basename = os.path.basename(sumatra_shape_path)[:-4]

# reproject
mp.get_stdout("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + 
    sumatra_shape_path[:-4] + "_repro.shp" + " " + sumatra_shape_path)
# the file sumatra_shape_path[:-4] + "_repro.shp" will be used for regions!

mp.get_stdout("ogrinfo -so -al " +  sumatra_shape_path[:-4] + "_repro.shp" )

# extract extent:
pre_extent = "-3292737.958267 962243.651532 -1981588.802986 2256489.436391"

mp.get_stdout("""gdal_rasterize -a_nodata 0 -a "CC_1" -l """+ sumatra_basename +"_repro -tr " + res + " -" + res + " -te " + pre_extent +" "+  sumatra_shape_path[:-4] + "_repro.shp " +  out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")

mp.get_stdout("gdalinfo " +  out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif")







infile =  out_path + "\\" + sumatra_basename[:-6] + "_"+ res + "_m_repro.tif"
unit_layer = mp.tif.read(infile, 1) 

units = np.unique(unit_layer)[1:]





# make unit area lookup table

unit_lookup = open(out_path  + "\\unit_area_lookup.csv",'w')
# write header
unit_lookup.write("unit,nr_pixel \n")
for unit in units: 
    print(unit)
    area = np.sum(np.where(unit_layer == unit, 1, 0))
    unit_lookup.write(str(unit) + "," + str(area) + " \n")


print("finished")

unit_lookup.close()  

