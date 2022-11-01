# 
# script to validate the future deforestation
# by checking perfect match, omission and comission errors in focal and surrounding pixels 
# we can only do this for years for which we have observations
# 2001:2016
# first only look at the respective year, later also the years prior ?
# first 2015, 2016

rm(list=ls())
library(tidyverse)
library(foreign)
library(ggplot2)

res2km <- 180/1000

PCVANT_path <- "C:/Users/mv296/work/Sumatra/deforestation_model/results/hansen"


units_list = rbind(c("Aceh", "Aceh", "A"),
                          c("N_Sumatra", "North Sumatra", "G"),
                          c( "W_Sumatra", "West Sumatra", "C"),
                          c("Riau", "Riau", "A"),
                          c("Jambi", "Jambi", "E"),
                          c("Bengkulu", "Bengkulu", "G"),
                          c("S_Sumatra", "South Sumatra", "A"),
                          c("Lampung", "Lampung",  "H"))

units_list <- as.data.frame(units_list)
names(units_list) <- c("unit_name","unit_name_long", "model_run")



# then also load the dbf so we know the codes
unit_codes<- read.dbf(file.path(PCVANT_path, "deforestation_quantification/sumatra_complete_shape_no_BB_repro.dbf")) %>% 
  dplyr::select(NAME_1, CC_1) %>% 
  mutate(unit_name = c("Aceh",
                       "Bengkulu",
                       "Jambi", 
                       "Lampung",
                       "Riau",
                       "W_Sumatra",
                       "S_Sumatra",
                       "N_Sumatra")) %>% 
  dplyr::select(unit = CC_1, unit_name)


# i need the unit area lookup
units_area_lookup <- read.csv(file.path(PCVANT_path, "deforestation_quantification/unit_area_lookup.csv"), stringsAsFactors = F) %>% 
  mutate(unit = as.factor(unit)) %>% 
  left_join(unit_codes, by = "unit") %>% 
  left_join(units_list, by = "unit_name")



# read csv

validation_data <- read.csv(file.path(PCVANT_path, "/validation/validation_full_map_hansen.csv"), stringsAsFactors = F) %>%
  rename(unit = province)%>%
  dplyr::select(-year, -scenario) %>%
 filter(match_type == "match" | match_type != "match" ) %>%
  rename(unit_name_long = unit) %>% 
  left_join(units_list, by = "unit_name_long") %>% 
  rename(unit_short = unit_name )
# exclude the out of radius 0 for comission / omission
  # nr deforested is the total number of deforested pixels within the province
  # but we are summing over the province, so this makes sense (match is total number of pixels with match)
 


#I AM LEAVING SUM DETECTED IN PIXEL; BECAUSE THIS IS MORE IMPORTANT THAN CONVERTING TO KM2 IMO

validation_total <- validation_data %>% 
  group_by(i, unit_short) %>% 
  summarise(total = sum(sum_detected))
  
validation_match <- validation_data %>%
  filter(match_type == "match" ) %>% 
  left_join(validation_total, by = c("i", "unit_short")) %>% 
  mutate(perc_of_total = sum_detected * 100/total)  # the value in radius is the number of pixels of 180m

validation_match_observed <- validation_data %>%
  filter(match_type == "match") %>% 
  left_join(validation_total, by = c("i", "unit_short")) %>% 
  mutate(perc_of_total= sum_detected * 100/total_forested)


validation_omission <- validation_data %>%
  filter(match_type == "omission") %>% 
  left_join(validation_total, by = c("i", "unit_short")) %>% 
  mutate(
    perc_of_total= sum_detected * 100/total_forested)

validation_comission <- validation_data %>%
  filter(match_type == "comission") %>% 
  left_join(validation_total, by = c("i", "unit_short")) %>% 
  mutate(
    perc_of_total= sum_detected * 100/total_forested_projected)


validation_all <- rbind(validation_match_observed, validation_omission, validation_comission)


validation_all$match_type <- factor(validation_all$match_type,
                                    levels = c("match", "omission", "comission"))

pdf(file.path(PCVANT_path, "validation/validation_all_forest.pdf"), paper = "a4r", width = 0, height = 0)

ggplot(data =validation_all, aes(x = as.factor(match_type), y = perc_of_total, colour = as.factor(match_type))) +  # at what are we looking here
  geom_boxplot(size = 1, outlier.size = 2, notchwidth = 0.5)+
  # xlab("radius in m")+
 #ylab("% match, omission and comission") +
  scale_color_manual(values = c("black", "blue", "green"))+
  facet_wrap(~unit_name_long)+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(from=0, to=100, by=20),
                     labels = seq(from=0, to=100, by=20))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 18, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=15),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 16, face="bold"))+
  theme(plot.margin = unit(c(0.2,-0.2,0.1,0.1),"cm"))

dev.off()

ggplot(data =validation_all, aes(x = as.factor(unit_name_long), y = perc_of_total, colour = as.factor(match_type))) +  # at what are we looking here
  geom_boxplot(size = 1, outlier.size = 2, notchwidth = 0.5)+
  # xlab("radius in m")+
  #ylab("% match, omission and comission") +
  scale_color_manual(values = c("black", "blue", "green"))+
  facet_wrap(~match_type)+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(from=0, to=100, by=20),
                     labels = seq(from=0, to=100, by=20))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 18, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=15),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 16, face="bold"))+
  theme(plot.margin = unit(c(0.2,-0.2,0.1,0.1),"cm"))


validation_summary_all <- validation_all %>%
  # filter(radius == 0) %>%
  group_by(unit_name_long, match_type) %>%
  summarize(median= round(median(perc_of_total)),
            lCI = round(quantile(perc_of_total, probs = 0.025, na.rm =T),2),
            uCI = round(quantile(perc_of_total, probs = 0.975, na.rm =T),2 ))

val_sum_S <- validation_all %>%
  group_by(match_type) %>%
  summarize(median= round(median(perc_of_total)),
            lCI = round(quantile(perc_of_total, probs = 0.025, na.rm =T),2),
            uCI = round(quantile(perc_of_total, probs = 0.975, na.rm =T),2))%>%
  mutate(unit_name_long = "Sumatra") %>%
  dplyr::select(unit_name_long,  match_type, median, lCI, uCI) %>%
  bind_rows(validation_summary_all) 


val_sum_S_out <- val_sum_S %>% 
  pivot_longer(cols= c(median, lCI, uCI)) %>% 
  pivot_wider(names_from = unit_name_long, values_from = value) %>% 
  as.data.frame() 

write.csv(val_sum_S_out, file.path(PCVANT_path, "validation/validation_table_out.csv"))


val_sum_S_all <- validation_all %>%
  # filter(radius == 0) %>%
  group_by(match_type) %>%
  summarize(median= round(median(perc_of_total)),
            lCI = round(quantile(perc_of_total, probs = 0.025, na.rm =T)),
            uCI = round(quantile(perc_of_total, probs = 0.975, na.rm =T)))%>%
  mutate(unit_name_long = "Sumatra") %>%
  dplyr::select(unit_name_long, match_type, median, lCI, uCI)

paste0("The overall prevalence of ",
      "perfect matches was ", val_sum_S_all[val_sum_S_all$match_type == "match" , "median"],
       "%, ",
       " false positives (commission errors) was at ",
       val_sum_S_all[val_sum_S_all$match_type == "comission" , "median"],
        "%, and the ",
       "prevalence of false negatives (omission errors) at ",
       val_sum_S_all[val_sum_S_all$match_type == "omission" , "median"],"%.")

val_sum_S_all <- val_sum_S_all%>%
  bind_rows(val_sum_S) %>% 
  as.data.frame()
