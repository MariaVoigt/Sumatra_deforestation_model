rm(list = ls())
library(foreign)
library(terra)
library(dplyr)
library(tidyr)
library(rgdal)

Sumatra_provinces <- c("Sumatera Selatan", "Sumatera Barat","Aceh", "Bengkulu",
                       "Lampung", "Jambi", "Sumatera Utara", "Riau") 

# inpath

podes_path = "P:/Sumatra/MPI/2018"

# clip this with sumatra borders
borders_shape <- readOGR(dsn = "N:/PODES_Boundaries/PODES_2018_boundaries",
                         layer = "PODES2018")

unique(borders_shape@data$PROPINSI)

borders_shape_sumatra <- borders_shape %>% 
  subset(PROPINSI %in% Sumatra_provinces)

#plot(borders_shape_sumatra)
# combine the deprivation minus environment



# import mpi  to join
mpi <- read.dbf(file.path(podes_path, "Sumatra_2018_MPI.dbf"))%>% 
  dplyr::filter(!is.na(DESA_ID)) 

str(mpi)

borders_shape_sumatra_data <- borders_shape_sumatra@data %>%
  mutate(DESA_ID = as.numeric(substr(ID, 3, 12))) %>% 
  left_join(mpi, by = "DESA_ID") %>% 
  dplyr::select(PROPINSI,
                KABUPATEN,
                KECAMATAN,
                DESA,
                ID, 
                DESA_ID, 
                LS_MPI_tot, # living standard
                H_MPI_tota, # health
                I_MPI_tota, # infrastructure
                Ed_MPI_tot, # education
                S_MPI_tota) %>% # social
  # add the total counts of all dimensions minus environment
  mutate(MPInoenv = LS_MPI_tot+ # living standard
                             H_MPI_tota+ # health
                             I_MPI_tota+ # infrastructure
                             Ed_MPI_tot+ # education
                             S_MPI_tota  , # social
                          
         MPInoenv_score = MPInoenv/15) #maximum count without environment

max(borders_shape_sumatra_data$MPInoenv_score, na.rm = T)

borders_shape_sumatra@data <- borders_shape_sumatra_data

summary( borders_shape_sumatra_data$MPInoenv_score)
# plot(borders_shape_sumatra)
# there are 154 NAs

# for  MPInoenv_score I will sample randomly from the distribution of values and assign those to the 
# missing values to not keep them at NA
plot(table(borders_shape_sumatra_data$MPInoenv_score))
s_NA <- sample(borders_shape_sumatra_data$MPInoenv_score, replace = F, size = sum(is.na(borders_shape_sumatra_data$MPInoenv_score)))
plot(table(s_NA))
borders_shape_sumatra_data[is.na(borders_shape_sumatra_data$MPInoenv_score), "MPInoenv_score"] <- s_NA
summary(borders_shape_sumatra_data$MPInoenv_score)
plot(table(borders_shape_sumatra_data$MPInoenv_score))

writeOGR(borders_shape_sumatra, 
        dsn = "C:/Users/mv296/work/Sumatra/data/model_input/PODES_processed",
        layer = "MPI_2018_sumatra",
        driver = "ESRI Shapefile", 
        overwrite = T,
        verbose = F)
