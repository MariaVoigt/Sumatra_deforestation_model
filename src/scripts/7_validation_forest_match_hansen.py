# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 13:00:55 2020

@author: mv39zilo
"""

# testing validation 
import numpy as np
import macpyver as mp
import csv 


def mw_sum(in_arr, x_win, y_win):

    out_arr = np.zeros(in_arr.shape)

    xn,yn = in_arr.shape

    for x in range(xn):

        xmin = max([0,x - x_win])

        xmax = min([xn, x + x_win])

        for y in range(yn):

            ymin = max([0,y - y_win])

            ymax = min([yn, y + y_win])

 

            out_arr[x,y] = in_arr[xmin:xmax+1, ymin:ymax+1].sum()

    return out_arr

"""
--> we are only working with the calibration, so year == 0, but would also be good maybe to see how it 
changes over year and whether then there is more overlap between observed and predicted

focal function
"""
# NEED TO FIX THESE
PCVANT_path = r'C:\Users\mv296\work\Wallacea\deforestation_model'

PCVANT_model_path = PCVANT_path + "\\" + r'Model_wallacea\src\final_pred_set'
out_path = PCVANT_path + "\\"+ r'results\validation'
# make a table with units and model-runs
#unit = "N_Sulawesi_Gorontalo"
#model_run =  "XX"
res = 180
forest_layer = 'hansen'

list_scenarios = [[ "Aceh", "Aceh", "A"],
                         ["N_Sumatra", "North Sumatra", "G"],
                         [ "W_Sumatra", "West Sumatra", "C"],
                         [ "Riau", "Riau", "A"],
                         [ "Jambi", "Jambi", "E"],
                         [ "Bengkulu", "Bengkulu", "G"],
                         [ "S_Sumatra", "South Sumatra", "A"],
                         [ "Lampung", "Lampung",  "H"]]
  
inpath =r'M:\Sumatra_model_August22' + "\\"+forest_layer

out_path = r"C:\Users\mv296\work\Sumatra\deforestation_model\results"


# make a table
year = str(0)


validation_data = open(out_path  + "\\"+ forest_layer +"\\validation_full_map_"+forest_layer+".csv",'w')

# write header
# province, scenario, year(0 for now), i(0:99), radius(0, 1Xres, 2Xres, 10Xres), match_type(perfect, omission, comission), nr_pixels
validation_data.write("province,scenario,year,i,match_type,sum_detected,total_forested, total_forested_projected \n")# what is the total numbers of pixels?

for province in range(0, len(list_scenarios)):
    model_unit = list_scenarios[province][1]
    list_unit = list_scenarios[province][0]
    scenario = list_scenarios[province][2]
    print(model_unit + ", " + list_unit + ", " + scenario)
    in_path_scenario = inpath + "\\" + list_unit + "\\model_" + list_unit + "_" + scenario + "\\" 
    year_0_deforest_path = in_path_scenario + "\\deforestation_2017_21_"+ str(res) + "m_repro_res_"+forest_layer+"_"+ list_unit+".ascii"
    year_0_forest_path = in_path_scenario + "\\forest_2017_21_"+str(res) + "m_repro_res_"+forest_layer+"_"+ list_unit+".ascii"
    deforest = np.loadtxt(year_0_deforest_path, skiprows=6)
    forest =    np.loadtxt(year_0_forest_path, skiprows=6)
    # forest in year 2021 for each province = forest_2017_21 (0+1) - deforest_2017_21 
   
    forest_2021 = np.where(deforest == 1, 0, forest)
    forest_2021 = np.where(forest_2021 == 0, - 9999, forest )
   # just take past deforestation events away
    forest_2016 = np.where(forest == 0, -9999, forest)

    # mp.tif.write(r"C:\Users\mv296\work\Wallacea\deforestation_model\data\model_input\units" + "/forest_2017_21_W_S_Sulawesi_210m_repro_res.tif", 
    #                      r"C:\Users\mv296\work\Wallacea\deforestation_model\data\model_input\units"  + "/forest_test2.tif" , 
    #                      forest_2016,
    #                      dtype = 0, 
    #                      nodata = -9999)
    nr_forested_2021 = np.sum(np.where(forest_2021==1, 1, 0))     
    nr_forested_2016 = np.sum(np.where(forest_2016==1, 1, 0))     
    for i in range(0,100):
#    for i in range(0,11):
        print(i)
        #CHECK THIS
        pred_path = inpath + "/" + list_unit + "/" + "model_" + list_unit +"_" + scenario +"/preddef_i"+ str(i) + "_"+ year+ "yr.asc"
        print(pred_path)
        pred_def = np.loadtxt(pred_path, skiprows=6)
        # predicted forest in 2021
        pred_forest_2021 = forest_2016-pred_def
        
        pred_nr_forested_2021 = np.sum(np.where(pred_forest_2021==1, 1, 0))     
        # omission: forest and in projection no forest
        omission_nr = np.where((forest_2021 == 1) & (pred_forest_2021 < 1), 1, 0).sum()
        validation_data.write(model_unit+ "," + scenario + "," + year + "," + 
                                        str(i) + "," +"omission," + str(omission_nr) +"," + 
                                        str(nr_forested_2021) + ","  + str(pred_nr_forested_2021) +" \n")
        # comission: no forest but in projection there is forest
        comission_nr = np.where((forest_2021 == -9999) & (pred_forest_2021 > 0), 1, 0).sum()
        validation_data.write(model_unit+ "," + scenario + "," + year + "," + 
                                        str(i) + ","  +"comission," + str(comission_nr) +"," +  
                                        str(nr_forested_2021) + ","  + str(pred_nr_forested_2021) +" \n")
            
        match_nr = np.where((forest_2021 == 1) & ( pred_forest_2021 == 1), 1, 0) .sum() # this is the total number of pixels with a match per province
        validation_data.write(model_unit+ "," + scenario + "," + year + "," + 
            str(i) + "," +"match," + str(match_nr) +"," +  
            str(nr_forested_2021)  + "," + str(pred_nr_forested_2021) +" \n")
              
validation_data.close()  

    