rm(list = ls())
library(tidyverse)
library(reshape2)
library(raster)
library(rgdal)
library(foreach)
#install.packages("corrgram")
library(corrgram)
#install.packages("EnvStats")
library(EnvStats)
#install.packages('car')
library(car)

# Script to plot the rates of deforestation 
#set working directory
predictors_path <- 
in_path <- ("C:/Users/mv296/work/Sumatra/data/model_input/repro_res/predictors_final")
outpath <- ("C:/Users/mv296/work/Sumatra/deforestation_model/results")
options(scipen=999)

res = "180"
# make one table with the id, deforestation and all parameters
# ADD ADDITIONAL REGIONS HERE

#read unit and variable names with this, but ultimately manually edited them
#unit_names <- sapply(list.files("C:/Users/mv296/work/Sumatra/data/model_input/units/", "*.dbf"), tools::file_path_sans_ext)
# variable_names <- sapply(list.files("C:/Users/mv296/work/Sumatra/data/model_input/repro_res/predictors_final/hansen", "*.tif"), tools::file_path_sans_ext)
  
unit_names <-  c("Aceh",
                 "N_Sumatra",
                 "W_Sumatra",
                 "Riau",
                 "Jambi",
                 "Bengkulu",
                 "S_Sumatra",
                # "B_Belitung",
                 "Lampung")
               
forest_layers <- c("tmf", "hansen")

variable_names <- c(  "forest_2017_21",
                      "slope",
                      "fire_yearly_average",
                      "IDN_TTCSM_hrs",
                      "pressurelog10_sigma1",
                      # "pressurelog10_sigma2", 
                      # "pressurelog10_sigma5", 
                      # "pressurelog10_sigma15", 
                      # "pressurelog10_sigma25", 
                      "pressurelog10_sigma50", 
                      "river_distance",           
                      "road_distance",      
                      "subsistence_LH_distance", 
                      "plantation_LH_distance",
                      "non_agri_LH_distance", 
                      "transmigrant_distance" ,
                      "small_plantations_distance", 
                      "ind_plantations_distance",
                      "peat",    
                      "mining",
                      "piaps",
                      "lu_new_class",
                      "soc_econMPI")

# 
# unit <- unit_names[1]
forest_layer <- "tmf"
unit = "Bengkulu"
# variable_names <- variable_names[c(2:3)]
# variable <- variable_names[2]

#for (forest_layer in forest_layers) {
#  for (unit in unit_names) {
    print(paste("Start unit ", unit, "at ", Sys.time()))
    in_path_unit <- file.path("C:/Users/mv296/work/Sumatra/data/model_input/repro_res/predictors_final", forest_layer , "asci")
    deforestation <-
      as.data.frame(raster(file.path(
        in_path_unit,
        paste0("forest_2017_21_", res, "m_repro_res_", forest_layer, "_",  unit,".ascii")
      )), xy = T) %>%
      mutate(id = paste0(round(x), "_", round(y)))
    
    names(deforestation) <- c("x", "y", "forest_2017_21", "id")
    
    for (variable in variable_names[1:length(variable_names)]) {
      print(paste(
        "Started variable ",
        variable,
        "for unit ",
        unit,
        "at ",
        Sys.time()
      ))
      variable_frame <-
        as.data.frame(raster(file.path(
          in_path_unit ,
          paste0(variable, "_",  res, "m_repro_res_", forest_layer, "_", unit, ".ascii")
        )), xy = T) %>%
        mutate(id = paste0(round(x), "_", round(y)))#%>%
      #   dplyr::select(-x,-y)
      print(paste("Finished loading ", variable, "for unit ", unit, "at ",  Sys.time() ))
      names(variable_frame) <- c("x", "y", variable, "id")
      print(paste("Loaded variable ",variable, "for forest",forest_layer,   "and unit",  unit,  "at ", Sys.time() ))
      deforestation[, variable] <- variable_frame[, variable]
    }
    
    
    print(paste("Finished all variables for unit ", unit, "at ", Sys.time()))
    
    deforestation_filter <-
      deforestation[complete.cases(deforestation),]
    
    
    deforestation_filter <-
      dplyr::select(deforestation_filter, -id,-x ,-y)
    deforestation_filter_cor <-   deforestation_filter
    
    names(deforestation_filter_cor) <-variable_names
    defor_cor <- cor(deforestation_filter_cor, method = c("pearson"))
    
    
    pdf(
      file.path(outpath, paste0("corrgram_", forest_layer,"_", unit , ".pdf")),
      width = 0,
      height = 0,
      paper = "a4r"
    )
    corrgram(
      defor_cor,
      lower.panel = panel.shade,
      upper.panel = panel.cor,
      main = unit
    )
    dev.off()



    # need to exclude these because these variables are aliased, i.e. directtly dependant on each other

    pred_names <-
      names(
        dplyr::select(
          deforestation,
          -id,
          -x,
          -y,
          # -pressurelog10_sigma2, #these are collinear
          # -pressurelog10_sigma5,
          # -pressurelog10_sigma15,
          -pressurelog10_sigma50,
          -plantation_LH_distance ,
         # -transmigrant_distance,
        #-ind_plantations_distance,
          -small_plantations_distance,
           # this doesnt exist in Aceh
        )      )
    
    
 #   pred_names <- variable_names
    defor_pred <- (deforestation[, pred_names])
    
    
    
    # check for multicoliniarity
    lm_res <-
      lm(paste("forest_2017_21 ~ ", paste(pred_names[2:length(pred_names)], collapse = "+")), data =
           defor_pred )
    
    # save this somewhere
    vif_res <- as.data.frame(vif(lm_res))
    # dataframe and hten save that
    vif_res
    
    write.csv(vif_res , file.path(outpath, paste0("vif_", forest_layer, "_", unit,".csv")))
    
    
    
    print(paste("Finished unit ", unit, "at ", Sys.time()))
    
    
  
  }
  }

