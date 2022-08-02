# script to combine and interprete 2020 piaps data
rm(list = ls())

library(foreign)
library(dplyr)
library(tools)
library(rgdal)
library(raster)
library(ggplot2)

piaps_path <- "N:/Landuse/Indonesia/PIAPS/PIAPS_Sept2020_Kraus/PIAPS_Sept2020_orig"
piaps_edit_path <- "N:/Landuse/Indonesia/PIAPS/PIAPS_Sept2020_Kraus/PIAPS_Sept2020_short_selected"
file_path_sans_ext(list.files(piaps_path, "*.dbf", full.names = F))

res = 100
res2ha = 0.0001
options(scipen = -999)

#"PIAPS-Revisi3" - not assigned to a specific scheme,
# under criteria they have Blok Pemberdayaan, Gambut Bebas Izin,
# Perhutani, Proses PS, Usulan PS, kelola sosial

#---Hutan Desa ----
# these are the rights
#"HPHD" - (Hak Pengelolaan Hutan Desa - Village Forest Management Rights
# # "Usulan_HD"- proposed HD


#------- Hutan Kemasyarakatan - Community utiliyation------
# "IUPHKM" - IUPHKm (Izin Usaha Pemanfaatan Hutan Kemasyarakatan) 
#                   Community Forest Utilization Business Permit), 
# Usulan_HKm 

#------Hutan Tanaman Rakyat (plantation/utilization) ------
# "IUPHHKHTR" - Izin Usaha Pemanfaatan Hasil Hutan Kayu dalam Hutan Tanaman Rakyat 
#             Business Permit for Utilization of Timber Forest Products in Community Plantation Forests) 
# --> seems to be related to HTR: "IPHPS" - Social Forestry Forest Utilization Permit (IPHPS)
# Usulan_HTR proposed plantation
# Usulan IPHPS  - proposed utilization permit

# Others --------
# "KULINKK" # this only on Java, Bali and Sumatra
# UsulanKulinKK 
#----------

#-------
#"PENETAPAN_PENCANTUMAN_HA" # determination of inclusion hectares



# "INDIKATIF_HA" - points needs excluding !!!

# what are important things:
# PS scheme
# status
# year of implementation (although not necessarily important)
# LU status  (although might be different ones)
# try and figure out what I should retain of each

# they have the area of landuse class, but I don't think it is very accurate or
# at least I am not sure, so rather calculate that again


# join them
available <- readOGR(file.path(piaps_path, "PIAPS-Revisi3.shp")) 
                      
available_d <- available@data %>%                  
  mutate(id = paste0(objectd, "_available"),
         type = NA,
         status = "available",
         permit = NA) %>% 
dplyr::select(id, 
              type,
              status, 
              permit)

available@data <- available_d
writeOGR(available, dsn = piaps_edit_path, layer = "available", overwrite = T, 
         driver = "ESRI Shapefile")

#"HPHD" - (Hak Pengelolaan Hutan Desa - Village Forest Management Rights
hphd <- readOGR(file.path(piaps_path, "HPHD.shp" ))
hd_implemented <- hphd@data %>% 
  mutate(id = paste0(objectd, "_HPHD"),
         type = "HD",
                        status= "implemented") %>% 
  
  dplyr::select(id, 
                type,
                status, 
                permit = n_sk_hp)
hphd@data <- hd_implemented

writeOGR(hphd, dsn = piaps_edit_path, layer = "hd_implemented", overwrite = T,
         driver = "ESRI Shapefile")


# # "Usulan_HD"- proposed HD
usulan_hd <- readOGR(file.path(piaps_path, "Usulan_HD.shp"))

hd_proposal <-usulan_hd@data %>% 
  mutate(id = paste0(objectid, "_proposed_HD"),
         type = "HD",
         status= "proposed",
         permit = NA) %>% 
  dplyr::select(id, 
                type,
                status, 
                permit)

usulan_hd@data <- hd_proposal

writeOGR(usulan_hd, dsn = piaps_edit_path, layer = "hd_proposal", overwrite = T,
         driver = "ESRI Shapefile")

# # "IUPHKM" - IUPHKm (Izin Usaha Pemanfaatan Hutan Kemasyarakatan) 
#                   Community Forest Utilization Business Permit), 


iuphkm <- readOGR(file.path(piaps_path, "IUPHKM.shp"))

hkm_implemented <-  iuphkm@data %>% 
  mutate(id = paste0(objectid_1, "_HKM"),
         type = "HKM",
         status= "implemented") %>% 
  
  dplyr::select(id, 
                type,
                status, 
                permit = no_sk_iuph)

iuphkm@data <- hkm_implemented

writeOGR(iuphkm, dsn = piaps_edit_path, layer = "hkm_implemented", overwrite = T,
         driver = "ESRI Shapefile")


# Usulan_HKm 
usulan_hkm <- readOGR(file.path(piaps_path, "Usulan_HKm.shp"))
hkm_proposal <-usulan_hkm@data %>% 
  mutate(id = paste0(1:nrow(.), "_proposed_HKM"),
         type = "HKM",
         status= "proposed",
         permit = NA) %>% 
  dplyr::select(id, 
                type,
                status, 
                permit)

usulan_hkm@data <- hkm_proposal

writeOGR(usulan_hkm, dsn = piaps_edit_path, layer = "hkm_proposal",
         driver = "ESRI Shapefile", overwrite = T)


# these are few areas, but I am including these
penetapan_pencantuman_ha <- readOGR(file.path(piaps_path, "PENETAPAN_PENCANTUMAN_HA.shp"))
added_areas <- penetapan_pencantuman_ha@data %>% 
  mutate(id = paste0(1:nrow(.), "_added_ha"),
         type = "added_ha",
         status= "proposed",
         permit = NA) %>% 
  dplyr::select(id, 
                type,
                status, 
                permit)

penetapan_pencantuman_ha@data <- added_areas


writeOGR(penetapan_pencantuman_ha, dsn = piaps_edit_path, layer = "added_areas",
         overwrite = T,
         driver = "ESRI Shapefile")

# also need to write the dbfs back out... and/or edit shapefile info


piaps <- raster::bind(penetapan_pencantuman_ha, usulan_hkm, iuphkm, 
                      usulan_hd, hphd, available)

piaps_data <- piaps@data
head(piaps_data)
# I want to extract the permit year
unique(piaps_data$permit)

piaps_data$id[duplicated(piaps_data$id)]

piaps_data$perm_yr <- stringr::str_extract(piaps_data$permit, "(?<=/20)[:digit:]{2}$")
testtest <- unique(piaps_data[ , c("permit", "perm_yr")])
unique(piaps_data$perm_yr)


test <- filter(piaps_data, !is.na(permit))
test2 <- filter(test, is.na(perm_yr))


difficult_years <- filter(piaps_data, !is.na(permit) & is.na(perm_yr))

difficult_years$perm_yr <- stringr::str_extract(difficult_years$permit, "(?<=20)[:digit:]{2}")


piaps_data[!is.na(piaps_data$permit) &
             is.na(piaps_data$perm_yr), "perm_yr"] <- stringr::str_extract(
               piaps_data[!is.na(piaps_data$permit) &
                                  is.na(piaps_data$perm_yr), "permit"],
                     "(?<=20)[:digit:]{2}")

piaps_data[!is.na(piaps_data$permit) &
             (piaps_data$permit == "SK.5906/MENLHK-PSKL/PKPS/PSL.0/10/2107"), "perm_yr"] <- 17
piaps_data[!is.na(piaps_data$permit) &
             (piaps_data$permit == "SK.6706/MENLHK-PSKL/PKPS/PSL.0/12/2107"), "perm_yr"] <- 17

test <- piaps_data[!is.na(piaps_data$permit) &
                     is.na(piaps_data$perm_yr), ]
# one observation left, but for that I erally dont know the year

test3 <- filter(piaps_data, !is.na(permit))
test4 <- filter(test3, is.na(perm_yr))

unique(piaps_data$perm_yr)
piaps_data$perm_yr <-as.integer(piaps_data$perm_yr) + 2000

sort(unique(piaps_data$perm_yr))

# add a joint id for presentation in arcgis
head(piaps)

piaps_data$label <- as.factor(paste0(piaps_data$type, "_", piaps_data$status))
unique(piaps_data$label)
piaps_data$label2 <- as.numeric(piaps_data$label)
unique(piaps_data$label2)
piaps_label_lookup <- piaps_data %>% 
  group_by(label) %>% 
  summarise(n = n(),
            label2 = unique(label2))

# here we want to give an id according to status

# implemented - 1
# proposed/available 2
piaps_data$status_id <- NA
piaps_data[piaps_data$status == "implemented", "status_id"] <- 1
piaps_data[piaps_data$status %in% c("proposed", "available"), "status_id"] <- 2



summary(piaps_data) 


piaps@data <- piaps_data

writeOGR(piaps, piaps_edit_path, "piaps_combined", overwrite = T, driver = "ESRI Shapefile")
