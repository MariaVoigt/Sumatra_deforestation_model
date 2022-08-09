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
    

km2_ha = 100
ha_km2 = 0.01
km2_Mha = 0.0001
m2_km2 = 0.000001
res = str(180)
res_f = float(res)



# read all files in, check what is in the output folder
# then np where with basemap
# then np where with forest, here split hansen and jrc forest

####
##  continue here!!!
##
forest_layers = ["tmf", "hansen"]


in_path_forest = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\forest'

in_path_base = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\base'


in_path_pred = r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\predictors'

pred_list = glob(in_path_pred + "//*.tif")


for forest_layer in forest_layers:
    print(forest_layer)
    for pred in pred_list:
       print(os.path.basename(pred))
       # 
       pred = mp.tif.read(pred, 1)
       pred =  np.where((base_map == 0), -9999, pred)
       pred = np.where(forest_1 == -9999, -9999, pred)
       mp.tif.write(...,
                    outpath-to-new-out-folder,
                    pred, 
                    ...)




# make the regions based on this file
# 'C:\\Users\\mv296\\work\\Sumatra\\data\\gadm\\sumatra_complete_shape_repro.shp

in_path = (r'I:\biocon\Maria_Voigt\deforestation_model\data\Gfw')

out_path = (r'I:\biocon\Maria_Voigt\deforestation_model\data\model_input\repro_res')

out_path_final = (r'I:\biocon\Maria_Voigt\deforestation_model\data\model_input\units')

out_path_asci = (r'I:\biocon\Maria_Voigt\deforestation_model\data\model_input\asci\units')
# then read it

                                             

# maybe all until here could happen in the other script?
                      
"""
 we will subdivide Wallacea in a number of units that are managable at a high resolution
 (sulawesi division based on historic admin borders before split)
- North-Sulawesi (Utara) + Gorontalo (combined until 2000)
- West (Barat) and South Sulawesi (Selatan) (split off in 2004)
- Central Sulawesi (tengah)
- Southeast Sulawesi (Tenggara)
- Maluku
-  Maluku will have to be split up in a Northern and a Southern part
    Name_2 for this division is:
        - Buru, Buru Selatan, Seram Bagian Barat, Ambon, Maluku Tengah, Seram Bagian Timur
        -Tual, Maluku Tenggara, Maluku Tenggara Barat, Maluku Barat Daya
- North Maluku (Utara)
"""
# reproject the gadm_1 layer
infile = r'I:\biocon\Maria_Voigt\deforestation_model\data\gadm\gadm36_IDN_shp\gadm36_IDN_1.shp'

os.system("""ogr2ogr -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + 
infile[:-4] + "_repro.shp" + " " + infile)


# extract extent, check if it is ok to use it
# then extract forests for this dimension
# turn into asci

# list of qualifiers


                

# should we clip everything with these dimensions or should we make a shapefile and use this to cut?
# I think a shapefile is better, becaues if we are doing it with admin units
# the shapes will be reallz odd and look weird in a publication
# in arcgis, you go in catalogue on right hand corner to folder, add new, polygon
# then reproject using ogr2ogr because would complain about projection even if selecting the right one
#then http://gis.yohman.com/up206a/how-tos/how-to-create-your-own-shapefile/

# HERE I NEED TO ADD NUSA TENGARRA
list_units = [["NAME_1 LIKE '%Sulawesi Utara%' OR NAME_1 LIKE '%Gorontalo%'", "N_Sulawesi_Gorontalo"], 
              ["NAME_1 LIKE '%Sulawesi Barat%' OR NAME_1 LIKE '%Sulawesi Selatan%'", "W_S_Sulawesi"], 
              ["NAME_1 LIKE '%Sulawesi Tenggara%'","SE_Sulawesi"],
              ["NAME_1 LIKE '%Sulawesi Tengah%'","C_Sulawesi"],
           #   ["NAME_1 LIKE '%Maluku%' AND NAME_1 NOT LIKE  '%Maluku Utara%'","Maluku"],   
              ["NAME_2 LIKE '%Buru%' OR NAME_2 LIKE '%Seram%' OR NAME_2 LIKE '%Ambon%' OR NAME_2 LIKE '%Maluku Tengah%'","C_Maluku"],   
              ["NAME_2 LIKE '%Tual%' OR NAME_2 LIKE '%Maluku Tenggara%' OR NAME_2 LIKE '%Maluku Barat Daya%'","S_Maluku"],   
              ["NAME_1 LIKE '%Maluku Utara%'","N_Maluku"]  , 
              ["NAME_1 LIKE '%Nusa Tenggara Barat%'","W_Nusa_Tenggara"],   
              ["NAME_1 LIKE '%Nusa Tenggara Timur%'","E_Nusa_Tenggara"]
              ]                  

              
predictor_names = [#"forest_2014_18", "deforestation_2014_18", "lu",
                   #"slope",
                #   "fire_yearly_average", #"op_suitability","op_suitability_recode", 
                 #  "gaez_op","gaez_coffee", "gaez_cacao","gaez_dry_rice", 
                  # "gaez_wet_rice","gaez_maize", 
                 #  "gaez_comb", "access_hrs",
                   "pop_pressurelog_sigma1", "pop_pressurelog_sigma2",  "pop_pressurelog_sigma5",
                  "pop_pressurelog_sigma15", "pop_pressurelog_sigma25", 
                   "pop_pressurelog_sigma50"#,
                #   "subsistence_distance_non_forest", "plantation_distance_non_forest", 
               #     "non_agri_distance_non_forest","fisheries_distance_non_forest",
               #    "transmigrant_distance" 
                 #  "mining"
                   ]        

in_path = "I:/biocon/Maria_Voigt/deforestation_model/data/wallacea_delim"

#list_unit = list_units[4]
for list_unit in list_units:    
    print("starting units " + str(list_unit))
    os.system("ogr2ogr -overwrite -where \"" + str(list_unit[0]) + """\" -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """ + in_path + "/units/"+ str(list_unit[1]) + ".shp " + in_path + "/wallacea_complete_shape_repro.shp")
    # Clip out the smaller units
    # to align their origin, I am extracting the corner information and add/subtract residual difference
    # from smaller extent    
    dimension_info_forest = subprocess.check_output(("gdalinfo " + out_path + "\\forest_2014_18_" + str(res) + "m_repro_res.tif "), shell=True).rstrip()
    lower_left_forest = after(dimension_info_forest, "Lower Left  (")[:25] 
    lower_left_forest_x = float([x.strip() for x in lower_left_forest.split(',')][0])
    lower_left_forest_y = float([x.strip() for x in lower_left_forest.split(',')][1])
    upper_right_forest = after(dimension_info_forest, "\r\nUpper Right ( ")[:24]
    upper_right_forest_x = float([x.strip() for x in upper_right_forest.split(',')][0])
    upper_right_forest_y = float([x.strip() for x in upper_right_forest.split(',')][1])

    dimension_info_unit = subprocess.check_output("ogrinfo -so -al " + in_path + '/units/'+str(list_unit[1]) + '.shp', shell=True).rstrip()  
    lower_left_unit = between(dimension_info_unit, "Extent: (", ") - (" ) 
    lower_left_unit_x = float([x.strip() for x in lower_left_unit.split(',')][0])
    lower_left_unit_y = float([x.strip() for x in lower_left_unit.split(',')][1])
    upper_right_unit = between(dimension_info_unit, ") - (", ")\r\nLayer" )
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
    os.system("gdal_rasterize -a_nodata 0 -burn 1 -l " + str(list_unit[1]) + '' + "   -tr " + res + " -" + res + " -te "+ str(te_info_new) + " " + in_path + "/units/"  + str(list_unit[1]) + '.shp ' + in_path + "/units/" + str(list_unit[1]) + '.tif')
    
    
    ref_unit_path = out_path_final + '/forest_2014_18_'+str(list_unit[1])+ "_" +str(res) + 'm_repro_res.tif'
    unit_delim = mp.raster.tiff.read_tif(in_path + "/units/" + str(list_unit[1]) + '.tif', 1)
    #clip the predictors
    for predictor_name in predictor_names:
        print(predictor_name)
        os.system("gdalwarp -overwrite -dstNodata -9999  -tr " + res + " -" + res + 
        " -te "+ te_info_new   + " " +  out_path + '\\' + str(predictor_name) + "_"+
        str(res) + 'm_repro_res.tif '+ out_path_final + '\\' + str(predictor_name) +
        "_" + str(list_unit[1])+ "_" +str(res) + 'm_repro_res.tif')  
        #read in,
        pred_unit = mp.raster.tiff.read_tif(out_path_final + '\\' +
        str(predictor_name) + "_" + str(list_unit[1])+ "_" +str(res) +
        'm_repro_res.tif', 1) 
        pred_unit_cut = np.where(unit_delim == 1, pred_unit, -9999)
        mp.raster.tiff.write_tif(ref_unit_path,
                             out_path_final + '\\' + str(predictor_name) + 
                             "_" + str(list_unit[1])+ "_" +str(res) + 
                             'm_cut_repro_res.tif',
                             pred_unit_cut, nodata=-9999, dtype = 4, option='compress=deflate')                         
        # translate to asci
        os.system("""gdal_translate -of "AAIGrid" -a_nodata -9999 """ + 
        out_path_final +  "\\" + str(predictor_name) + "_" +
        str(list_unit[1]) + "_"+ str(res) + 'm_cut_repro_res.tif ' + 
        out_path_asci + "\\" +  str(predictor_name) + "_" +str(list_unit[1])+
        "_" +str(res) + 'm_repro_res.ascii' )
   

