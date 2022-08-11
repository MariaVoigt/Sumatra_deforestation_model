# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:18:44 2019

@author: mv39zilo

# prepares predictors

# first filter the forest and basemap
# then filter for provinces
# also think about border regions this time
"""

import numpy as np
import os
from osgeo import ogr
import shutil
from glob import glob
import subprocess
import math

import macpyver as mp

    
# functions
# this extracts texts after or between string bits    
def after(value, a): # from here https://www.dotnetperls.com/between-before-after-python
    # Find and validate first part.
    pos_a = value.rfind(a)
    if pos_a == -1: return ""
    # Returns chars after the found string.
    adjusted_pos_a = pos_a + len(a)
    if adjusted_pos_a >= len(value): return ""
    return value[adjusted_pos_a:]
    
def between(value, a, b):
    # Find and validate before-part.
    pos_a = value.find(a)
    if pos_a == -1: return ""
    # Find and validate after part.
    pos_b = value.rfind(b)
    if pos_b == -1: return ""
    # Return middle part.
    adjusted_pos_a = pos_a + len(a)
    if adjusted_pos_a >= pos_b: return ""
    return value[adjusted_pos_a:pos_b]
    

# this compares the extent of rasters
def compare(rast_path_1, rast_path_2):
    rast1_ext = mp.raster.tiff.get_extent(rast_path_1)
    rast2_ext = mp.raster.tiff.get_extent(rast_path_2)
    error= False
    if rast1_ext.px_size != rast2_ext.px_size:
        error = True
    if rast1_ext.columns != rast2_ext.columns:
        error = True
    if rast1_ext.left != rast2_ext.left:
        error = True
    if error:
        return('files are not matching')
    else:
        print('alles ist gut')

km2_ha = 100
ha_km2 = 0.01
km2_Mha = 0.0001
m2_km2 = 0.000001
res = str(180)
res_f = float(res)



# read all files in, check what is in the output folder
# then np where with basemap
# then np where with forest, here split hansen and jrc forest

##
in_path_forest = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\forest'
in_path_base = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\base'
in_path_pred = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\predictors'

out_path_pred = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\predictors_final'

base_map = mp.tif.read(in_path_base + "\\sumatra_complete_"+ res + "_m_repro.tif")
units_path = r'C:\Users\mv296\work\Sumatra\data\model_input\units'


forest_layer_types = ["tmf", "hansen"]
years = "2017_21"
res = str(180)

pred_list = glob(in_path_pred + "//*.tif")

"""
 we will subdivide Sumatra in a number of units that are managable at a high resolution

# there are some potential issues, mainly Bangka Belitung being relatively small 
# and some others like Sumatera Utara being too large
# first we start with admin boundaries though
# also wanted to see if I can explore the boundary regions as well


"""

sumatra_shape_path = r'N:\Admin_boundaries\GADM\Sumatra\sumatra_complete_shape_repro.shp'        


list_units = [["NAME_1 LIKE '%Aceh%'", "Aceh"],
              ["NAME_1 LIKE '%Sumatera Barat%'","W_Sumatra"],
              ["NAME_1 LIKE '%Sumatera Selatan%'","S_Sumatra"],
              ["NAME_1 LIKE '%Sumatera Utara%'","N_Sumatra"],
              ["NAME_1 LIKE '%Jambi%'", "Jambi",],
              ["NAME_1 LIKE '%Bengkulu%'", "Bengkulu"],
              ["NAME_1 LIKE '%Bangka Belitung%'", "B_Belitung"],
              ["NAME_1 LIKE '%Lampung%'", "Lampung"],
              ["NAME_1 LIKE '%Riau%' AND NAME_1 NOT LIKE '%Kepulauan Riau%'", "Riau"]]
              
             
predictor_names = [os.path.basename(x)[:-4] for x in glob(out_path_pred + "//*.tif")]   

# get forest dimensions to align the regions to 
dimension_info_forest = mp.get_stdout("gdalinfo " + in_path_forest + "\\forest_tmf_" + years + "_" + res + "m_repro_res.tif")
dimension_info_forest = " ".join(dimension_info_forest)
dimension_info_forest 
lower_left_forest = after(dimension_info_forest, "Lower Left  (")[:25] 
lower_left_forest_x = float([x.strip() for x in lower_left_forest.split(',')][0])
lower_left_forest_y = float([x.strip() for x in lower_left_forest.split(',')][1])
upper_right_forest = after(dimension_info_forest, "Upper Right (")[:24]
upper_right_forest_x = float([x.strip() for x in upper_right_forest.split(',')][0])
upper_right_forest_y = float([x.strip() for x in upper_right_forest.split(',')][1])


# Clip forest layers with the basemap

for forest_layer_type in forest_layer_types:
    # forest_layer_type = forest_layer_types[0]
    print(forest_layer_type)
    # we need to clip the forest layers with base map as well
    forest_1  = mp.tif.read(in_path_forest + "\\" + "forest_" + forest_layer_type + "_" + years + "_" + res + "m_repro_res.tif", 1)
    forest_2  = mp.tif.read(in_path_forest + "\\" + "deforestation_" + forest_layer_type + "_" + years + "_" + res + "m_repro_res.tif", 1)

    forest_1 =  np.where((base_map == 1), forest_1, -9999)
    mp.tif.write(in_path_forest + "\\" + "forest_" + forest_layer_type + "_" + years + "_" + res + "m_repro_res.tif",
                    out_path_pred  + "\\" + forest_layer_type + "\\" + "forest_" + years + "_" + res + "m_repro_res_"  + forest_layer_type +".tif",
                    forest_1, 
                    nodata = -9999, 
                    option='compress=deflate')
    
    forest_2 =  np.where((base_map == 1), forest_2, -9999)
    mp.tif.write(in_path_forest + "\\" + "deforestation_" + forest_layer_type + "_" + years + "_" + res + "m_repro_res.tif",
                    out_path_pred  + "\\" + forest_layer_type + "\\" + "deforestation_"  + years + "_" + res + "m_repro_res_"  + forest_layer_type +".tif",
                    forest_2, 
                    nodata = -9999, 
                    option='compress=deflate')

    for pred in pred_list:
     #pred = pred_list[0]
       pred_name = os.path.basename(pred)[:-4]
       print(pred_name)
       pred_layer = mp.tif.read(pred, 1)
       pred_layer = np.where(forest_1  == -9999, -9999, pred_layer)
       mp.tif.write(pred,
                    out_path_pred+ "\\" + forest_layer_type + "\\" + pred_name + "_" + forest_layer_type + ".tif",
                    pred_layer, 
                    nodata = -9999, 
                    option='compress=deflate')

# CONTINUE HERE INCLUDING FOREST IN PREDICTORS
    
#list_unit = list_units[4]
for list_unit in list_units:    
    print("starting units " + str(list_unit))
    mp.get_stdout("ogr2ogr -overwrite -where \"" + str(list_unit[0]) + """\" -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + units_path+ "//"+ str(list_unit[1]) + ".shp " + sumatra_shape_path)
   
    # Clip out the smaller units
    # to align their origin, I am extracting the corner information and add/subtract residual difference
    # from smaller extent    


    dimension_info_unit = mp.get_stdout("ogrinfo -so -al " + units_path+ "//"+ str(list_unit[1]) + ".shp")
    dimension_info_unit = " ".join(dimension_info_unit)
    
    lower_left_unit = between(dimension_info_unit, "Extent: (", ") - (" ) 
    lower_left_unit_x = float([x.strip() for x in lower_left_unit.split(',')][0])
    lower_left_unit_y = float([x.strip() for x in lower_left_unit.split(',')][1])
    upper_right_unit = between(dimension_info_unit, ") - (", ") Layer " )
    upper_right_unit_x = float([x.strip() for x in upper_right_unit.split(',')][0])
    upper_right_unit_y =float( [x.strip() for x in upper_right_unit.split(',')][1])
    
    # this finds the residual between the differences in origins and subtracts this from the lower
    # corner and adds it to the upper corner
    # to take out the offset between the layers

    #lower left corner
    residual_ll_x = ( np.abs(lower_left_forest_x - lower_left_unit_x) /  res_f - 
                 np.floor(np.abs(lower_left_forest_x - lower_left_unit_x) / res_f)) * res_f

    residual_ll_y = ( np.abs(lower_left_forest_y - lower_left_unit_y) / res_f - 
                np.floor(np.abs(lower_left_forest_y - lower_left_unit_y) / res_f)) * res_f

    lower_left_unit_x_new = lower_left_unit_x - residual_ll_x 
    lower_left_unit_y_new = lower_left_unit_y - residual_ll_y

    #upper right corner
    residual_ur_x = ( np.abs(upper_right_forest_x - upper_right_unit_x) / res_f - 
                np.floor(np.abs(upper_right_forest_x - upper_right_unit_x) / res_f)) * res_f

    residual_ur_y = ( np.abs(upper_right_forest_y - upper_right_unit_y) / res_f -
                 np.floor(np.abs(upper_right_forest_y - upper_right_unit_y) / res_f)) * res_f

    upper_right_unit_x_new = upper_right_unit_x + residual_ur_x 
    upper_right_unit_y_new = upper_right_unit_y + residual_ur_y
    
    # put it all together in the new te info
    te_info_new = str(lower_left_unit_x_new) + " " + str(lower_left_unit_y_new) + " " +  str(upper_right_unit_x_new) + " " + str(upper_right_unit_y_new)

    print(str(te_info_new))
    # clip the state, so I can clip out the pixels that fall within extent but not state

    mp.get_stdout("gdal_rasterize -a_nodata 0 -burn 1 -l " + str(list_unit[1]) + '' + "   -tr " + res + " -" + res + " -te "+ str(te_info_new) + " " + units_path + "//"  + str(list_unit[1]) + '.shp ' + units_path +"//" + str(list_unit[1]) + '.tif')
    
    ref_unit_path = (units_path +"//" + str(list_unit[1]) + '.tif')
    unit_delim = mp.tif.read(ref_unit_path, 1)
    
    for forest_layer_type in forest_layer_types:
    # forest_layer_type = forest_layer_types[0]
        print(forest_layer_type)
        out_path_units = out_path_pred+ "\\" + forest_layer_type + "\\units"
        out_path_asci = out_path_pred+ "\\" + forest_layer_type + "\\asci"

       # update pred list to include the two forest layers
        pred_list = glob( out_path_pred+ "\\" + forest_layer_type + "//*.tif")

        #clip the predictors and convert to asci
        for pred in pred_list:
            #pred = pred_list[0]
            pred_name = os.path.basename(pred)[:-4]
            print(pred_name)
            mp.get_stdout("gdalwarp -overwrite -dstNodata -9999  -tr " + res + " -" + res + 
            " -te "+ te_info_new   + " " +   
            out_path_pred+ "\\" + forest_layer_type + "\\" + pred_name + "_" + forest_layer_type + ".tif " + 
            out_path_units + '\\' + str(pred_name) +
            "_" + str(list_unit[1])+  "_" + forest_layer_type + ".tif")  
            #read in,
            pred_unit = mp.tif.read( out_path_units + '\\' + str(pred_name) +
            "_" + str(list_unit[1])+  "_" + forest_layer_type + ".tif", 1) 
            pred_unit_cut = np.where(unit_delim == 1, pred_unit, -9999)
            mp.tif.write(ref_unit_path,
                               out_path_units + '\\' + str(pred_name)+
                               "_" + str(list_unit[1])+  "_" + forest_layer_type + "_cut.tif",
                                 pred_unit_cut, nodata=-9999, dtype = 4, option='compress=deflate')                         
            # translate to asci
            mp.get_stdout("""gdal_translate -of "AAIGrid" -a_nodata -9999 """ +  out_path_units + "\\" + str(pred_name)+"_" + str(list_unit[1])+  "_" + forest_layer_type + "_cut.tif "+   out_path_asci + "\\" + str(pred_name)+ "_" + str(list_unit[1])+  "_" + forest_layer_type + ".ascii" )
       
    
