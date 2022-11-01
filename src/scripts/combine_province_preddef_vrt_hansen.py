# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 21:36:39 2019

@author: mv39zilo

join the deforestation maps from each province

"""


import numpy as np
import os

import glob

import macpyver as mp

km2_ha = 100

res = str(180)


forest_layer = 'hansen'


# list_units_tmf = rbind(c(1, "Aceh", "Aceh", "A"),
#                        c(2, "N_Sumatra", "North Sumatra", "B"),
#                        c(3, "W_Sumatra", "West Sumatra", "B"),
#                        c(4, "Riau", "Riau", "A"),
#                        c(5, "Jambi", "Jambi", "D"),
#                        c(6, "Bengkulu", "Bengkulu", "F"),
#                        c(7, "S_Sumatra", "South Sumatra", "A"),
#                        c(8, "Lampung", "Lampung",  "I"))


list_units = [[ "Aceh", "Aceh", "A"],
                         ["N_Sumatra", "North Sumatra", "G"],
                         [ "W_Sumatra", "West Sumatra", "C"],
                         [ "Riau", "Riau", "A"],
                         [ "Jambi", "Jambi", "E"],
                         [ "Bengkulu", "Bengkulu", "G"],
                         [ "S_Sumatra", "South Sumatra", "A"],
                         [ "Lampung", "Lampung",  "H"]
]

#NAS_path = r'C:\Users\mv296\SynologyDrive\My_drive\SynologyDrive\Sumatra_model_August22'
drive_path = r'E:\Sumatra_model_August22'


out_path = drive_path + '\\' + forest_layer + '\\' + "projected_deforest\provinces_combined"
unit_shape_path = r'C:\Users\mv296\work\Sumatra\data\model_input\units'
        


# Sumatra extent


#we can't use a base path anymore because we have excluded bangka belitung for now
# so I created a new one  
# getting the extent from running a vrt further down 
extent = "-3292738.005  962109.389 -1981438.005 2256489.389"

base_path = (r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\predictors_final' +
 '\\' + forest_layer +  '\\' + "forest_2017_21_" + res + "m_repro_res_" + 
 forest_layer + ".tif")

 
mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent + " "+
     base_path + " " + 
     out_path + "\\" + "forest_2017_21_" + res + "m_repro_res_" +  forest_layer + "_base.tif")
 

base_map = mp.tif.read( out_path + "\\" + "forest_2017_21_" + res + "m_repro_res_" +  forest_layer + "_base.tif", 1)


in_path = drive_path + '\\' + forest_layer 


# extend all the province tifs to sumatra

# gdal build vrt for each year and i (nodata of asci is 0)
# then convert to raster
# nodata is -9999

#dalbuildvrt -separate stacked.vrt [in vrts or rasters]
#gdal_translate stacked.vrt stacked.tif

right_shape = (7191, 7285)

for i in range(76,100):
#for i in range(0,5):    
    print("i" + str(i))
    # initialize empty sumprob
    # will already expand extend for each year, because we need this later for not cummulative years
    for year in range(0, 8): # year 0 not included is included, because we need that for TS analysis
        print("year" + str(year))
        preddef_list_i_yr = []
        # combine all units here
        for unit in range(0, len(list_units)):
            model_unit = list_units[unit][1]
            list_unit = list_units[unit][0]
            scenario = list_units[unit][2]
            print(model_unit + ", " + list_unit + ", " + scenario)
            province_path = unit_shape_path + "\\" +  list_unit + ".tif"
            in_path_unit =  "\""+ in_path + "\\" + list_unit + "\\" + "model_" + list_unit +"_" + scenario + "\\" + "preddef_i" + str(i) + "_" + str(year) + "yr.asc"+ "\" "
            # make a list of asci files
            preddef_list_i_yr = preddef_list_i_yr + [in_path_unit] # how to make this list here
            
        os.system("""gdalbuildvrt -r "near" -overwrite -srcnodata 0 -a_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def"  """  + 
                  out_path + "\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".vrt "+ " ".join(preddef_list_i_yr))
        os.system("gdal_translate "+ out_path + "\\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".vrt " + 
                  out_path + "\\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".tif") 
        
        mp.get_stdout(("gdalinfo " + out_path + "\\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".tif" ))

            # this worked, although no wthere is 1 and zero
        preddef_combined_i_yr = mp.tif.read(out_path + "\\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".tif", 1)   
    
        # this filters out past deforestation and zeros in preddef
        # but what about NAs in base map
        preddef_combined_i_yr = np.where(preddef_combined_i_yr == 0, -9999, preddef_combined_i_yr)
            # at the border between states there can be an overlap
        preddef_combined_i_yr = np.where(preddef_combined_i_yr < 0, -9999, preddef_combined_i_yr)
        preddef_combined_i_yr = np.where(preddef_combined_i_yr > 1, 1, preddef_combined_i_yr)  
        if np.shape(preddef_combined_i_yr) != right_shape:
            print("shape not right: " + str(np.shape(preddef_combined_i_yr)))
            break
    # I need the non-cummulative preddef for the future_forest_loss_5y script
        mp.tif.write( out_path + "\\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".tif" ,
                    out_path + "\\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".tif",
                     preddef_combined_i_yr,
                     nodata = -9999,
                     option='compress=deflate')

#check whether all output tifs were written
total_number_files = 100 *8
list_files = glob.glob(out_path + "\\Sumatra_"+ forest_layer +"_preddef_i*_yr[0-7].tif")
len(list_files) == total_number_files
    
#THIS IS CUMMULATIVE PREDDEF
#I already have the units with the extent of borneo


path_defor_17_21 = (r'C:\Users\mv296\work\Sumatra\data\model_input\repro_res\predictors_final' +
                    '\\' + forest_layer +  '\\' + "deforestation_2017_21_" + res + "m_repro_res_" + 
 forest_layer + ".tif")

mp.get_stdout("""gdalwarp -r "near" -overwrite -t_srs "+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_def" """+
        "-tr " + res + " -" + res + " -te " + extent + " "+
        path_defor_17_21 + " " + 
        out_path + "\\" + "deforestation_2017_21_" + res + "m_repro_res_" +  forest_layer + "_base.tif")


defor_17_21 = mp.tif.read(out_path + "\\" + "deforestation_2017_21_" + res + "m_repro_res_" +  forest_layer + "_base.tif", 1)



forest_17_21 = base_map 


defor_00_16 = np.where(forest_17_21 == 0, 1, -9999)

forest_17_21 = np.where(forest_17_21 == 0, -9999, forest_17_21)


for i in range(0,100):
#for i in range(0,5):
    print("i"+str(i))
    # initialize empty sumprob
    preddef_yr_cum = np.empty(base_map.shape)
          # will already expand extend for each year, because we need this later for not cummulative years
 
    for year in range(1, 8): # year 0 not included
        print("year" + str(year))
        #import year_i_sumatra one 
        in_path_sumatra_i_yr = out_path + "\\Sumatra_"+ forest_layer +"_preddef_i"+ str(i) + "_yr" + str(year) + ".tif"
        preddef_i_yr =  mp.tif.read(in_path_sumatra_i_yr, 1)
        if year == 1:
            preddef_yr_cum = preddef_i_yr
        else: 
            preddef_yr_cum =np.where((preddef_i_yr == 1) & (preddef_yr_cum != -9999), preddef_yr_cum + preddef_i_yr, preddef_yr_cum)
            preddef_yr_cum =np.where((preddef_i_yr == 1) & (preddef_yr_cum == -9999), 1, preddef_yr_cum)
            # only save if it is not the first year as there is nothing cummulative
            mp.tif.write(in_path_sumatra_i_yr,
                                 in_path_sumatra_i_yr[:-4 ] +  "cum.tif",
                                 preddef_yr_cum,
                                 nodata = -9999, 
                                 option='compress=deflate')
            
        # add the actual deforestation in year 0 for mapping
        # in repro_res
        preddef_yr_cum_cor = np.where(defor_17_21 == 1, 1, preddef_yr_cum)
        out_path_preddef_yr_cum_cor =  in_path_sumatra_i_yr[:-4 ] + "_defor_17_21.tif"
        mp.tif.write( in_path_sumatra_i_yr,
                out_path_preddef_yr_cum_cor,
                preddef_yr_cum_cor,
                nodata = -9999, option='compress=deflate')   
    
        
        # add the deforestation in year 2000-2013 in a different number
        preddef_yr_cum_cor_past_defor = np.where(defor_00_16 == 1, 2, preddef_yr_cum_cor)
        out_path_preddef_yr_cum_cor_past_defor = in_path_sumatra_i_yr[:-4 ] + "_defor_01_18.tif"
        mp.tif.write( in_path_sumatra_i_yr,
                        out_path_preddef_yr_cum_cor_past_defor,
                           preddef_yr_cum_cor_past_defor, nodata = -9999, option='compress=deflate')   
        
        # compute remaining forest (where there is a 1 in forest_17_21 and a 1 in 
        # predicted defor, there a 0)
        forest_yr_cum_cor = np.where(preddef_yr_cum_cor == 1, -9999, forest_17_21)
        out_path_forest_cum_cor_past_defor = out_path + "\\sumatra_forest_i"+str(i) + "_yr" + str(year) + "cum_cor.tif"
        mp.tif.write(in_path_sumatra_i_yr,
                        out_path_forest_cum_cor_past_defor,
                       forest_yr_cum_cor, nodata = -9999, option='compress=deflate')   
    
        
     
# to make sumprob I only need to add up all 
# we have years 1 to 7, year 1 is 2022+5 etc
years = list()
sumprobs = list()
start_sumprob = 2021
shape_out = forest_17_21.shape


# create base map for sumprob, to know what is na and what isnt
#N:\Admin_boundaries\GADM\Sumatra


for year in range(1, 8):
    print(year)

    if year == 1:
        start_sumprob = start_sumprob +1
    else:
        start_sumprob = start_sumprob + 5
    end_sumprob = (start_sumprob-1) + 5
    
    if year == 1:
        # here we have no sumprob, so 
        in_path_sumatra_yr_cum = out_path + "\\Sumatra_"+ forest_layer +"_preddef_i*_yr" + str(year) + ".tif"
    else:
        in_path_sumatra_yr_cum = out_path + "\\Sumatra_"+ forest_layer +"_preddef_i*_yr" + str(year) + "cum.tif"
    preddef_list_yr = glob.glob(in_path_sumatra_yr_cum)
    # make an empty tif file 
    sumprob =  np.zeros(shape_out)
    sumprob_path = out_path + "\\Sumatra_"+ forest_layer +"_yr" + str(year) + "_sumprob.tif"

        
    for i in range(0,len(preddef_list_yr)):
    #for i in range(0,10):
        print(preddef_list_yr[i])
        pred_def_i = mp.tif.read(preddef_list_yr[i], 1)
        sumprob = np.where(pred_def_i == -9999, sumprob, sumprob + pred_def_i)
    mp.tif.write(out_path + "\\" + "deforestation_2017_21_" + res + "m_repro_res_" +  forest_layer + "_base.tif",
                sumprob_path,
                sumprob, nodata = -9999, option='compress=deflate')      
        