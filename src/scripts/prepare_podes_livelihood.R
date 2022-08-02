# assemble sumatra wide podes 

rm(list = ls())
library(foreign)
library(foreach)
library(tidyverse)
library(stringr)
options(scipen = 999)
# descend into all files 
# get from the first 2018 one the livelihood for each region, then put in long format and plot

#files on NAS under podes mapped to drive P
in_path <- "P:/Sumatra/Raw data"
# gettint all files in the 2018 files
# a bit clumsy because most aer called podes_desa* the for one province only podes!

dirs <- list.dirs(in_path)
dir_list <- dirs[ grepl("2018", dirs) ]

file_list <- list.files(dir_list , pattern = "*.dbf", recursive = T)


#in_path <-file.path(PCVANT_path, "data/PODES/Sulawesi_Maluku_complete") 
#out_path <- file.path(PCVANT_path, "data/PODES/results")

boundaries <- read.dbf("N:/PODES_Boundaries/PODES_2018_boundaries/PODES2018.dbf") 

livelihood_data <- as.data.frame(foreach (dir = dir_list, .combine = rbind)%do%{
#for(dir in dir_list){
  path = list.files(dir, "*_d1.dbf", full.names = T )
  region_name = strsplit(dir, "/")[[1]][4]
  print(region_name)
  podes_data <- read.dbf(path) %>% 
  
  # pad colums if necessary
  
    mutate(ID = paste0("ID",  R101,  
                       str_pad(R102,2,"left",pad = "0"),
                       str_pad(R103,3, "left", pad = "0"),
                       str_pad(R104, 3, "left", pad = "0")
                       )) %>% 
    dplyr::select(ID, R403A, R403B)
# this is a side effect, be careful here  
    
})

boundaries <- left_join(boundaries,livelihood_data, by = "ID")

summary(boundaries)

# rowbind them all, then left join

write.dbf(boundaries, "C:/Users/mv296/work/Sumatra/data/PODES_2018_boundaries/PODES2018.dbf")



names(livelihood_data) <- c("region", "question", "answer")
str(livelihood_data)
livelihood_data$answer <- as.integer(as.character(livelihood_data$answer) )
# separate the two questions and join with lookup


livelihood_name <- c("Agriculture", "Mining and Excavation","Manufacturing industry (factories, crafts, etc.)",
                     "Wholesale / retail trade and restaurants", "Transportation, warehousing, communication",
                     "Services", "Other")

lookup_403A <- data.frame(answer= c(1:7), livelihood_name = as.factor(livelihood_name))

lookup_403A$livelihood_name <- factor(livelihood_name, levels = livelihood_name)


livelihood_data_A <- livelihood_data %>%
  filter(question == "R403A")%>%
  left_join(lookup_403A, by="answer" )

pdf(file.path(out_path , "main_source_income.pdf"), paper = "a4r", width = 0, height = 0)
ggplot()+
  geom_bar(data = livelihood_data_A,
           aes(x = as.factor(livelihood_name), fill = as.factor(region)))+
  labs( x="main source of income")+
  #facet_wrap(~region)+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x=element_text(size = 10, angle = 90, hjust=0.5,vjust=0.2),
    # axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 10, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=10),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 10, face="bold")) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

dev.off()

activity_name <- c(
  "Rice", "Palawija (corn, beans, sweet potatoes)", 
  "Horticulture (fruits, vegetables, ornamental plants, medicinal plants, etc.)",
  "Rubber","Palm oil", "Coffee","Cocoa","Coconut", "Pepper","Cloves","Tobacco","Sugar cane",
  "Animal Husbandry (cattle, sheep, chickens, etc.)", "Capture fisheries (including other biota)",
  "Aquaculture (including other biota)","Forestry cultivation (teak, mahogany, sengon, bamboo, etc.)",
  "Collection of forest products (honey, agarwood, fruits, firewood, etc.)", 
  "Capture of wild animals (pigs, partridges, deer, etc.)",
  "Captivity of wild animals / plants (arowana, crocodile, orchid, etc.)",
  "Agricultural services (hatchery, tractor rental, rattan, etc.)", "Other")



lookup_403B <- data.frame(answer= c(1:21),  activity_name = as.factor(activity_name))

lookup_403B$activity_name <- factor(activity_name, levels = activity_name)


livelihood_data_B <- livelihood_data %>%
  filter(question == "R403B")%>%
  left_join(lookup_403B, by="answer" )

#livelihood_data_B$activity_name <- as.factor(livelihood_data_B$activity_name)

#levels(livelihood_data_B$activity_name) <- activity_name

pdf(file.path(out_path, "main_commodity.pdf"), paper = "a4r", width = 0, height = 0)
ggplot()+
  geom_bar(data = livelihood_data_B,
           aes(x=activity_name, fill = as.factor(region)))+
  labs( x="main commodity produced by villages")+
  
  #facet_wrap(~region)+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x=element_text(size = 10, angle = 90, hjust=1,vjust=0.2),
    # axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 10, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=10),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 10, face="bold")) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

dev.off()
