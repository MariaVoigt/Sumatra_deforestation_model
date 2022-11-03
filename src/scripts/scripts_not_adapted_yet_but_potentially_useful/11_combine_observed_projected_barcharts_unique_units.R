# script to produce both observed and predicted deforestation in Wallacea

rm(list=ls())
#install.packages('tidyverse')
library(tidyverse)
#install.packages('foreach')
library(foreach)
#install.packages('reshape2')
library(reshape2)
#install.packages('mgcv')
library(mgcv)

options("scipen" = 100, "digits" = 10)

res2km2 = (180*180)/(1000*1000)

PCVANT_path <- "C:/Users/mv296/work/Wallacea/deforestation_model"

calc.annual.change <- function(value_start, value_end, year_start, year_end, round_to = 2){
  result <- (((value_end / value_start )^(1/(year_end - year_start))) - 1) * -100
  return(round(result, round_to))
}




#----------------#
#   OBSERVATION  #
#----------------#


# 1 year
# calculate_observed_change.py
data <- read.csv(file.path(PCVANT_path, "results/deforestation_quantification/observed_yearly_deforestation_2018.csv"), stringsAsFactors = F)
#CHECKED


# deforestation rate
forest_data <- data %>%
  dplyr::filter(type == "loss"  & unit != "Wallacea") %>%
  mutate(year = year + 2000)

forest_data[!is.na(forest_data$nr_pixel), "abs_loss_km2"] <- as.numeric(forest_data[!is.na(forest_data$nr_pixel), "nr_pixel"])*res2km2


units <- unique(forest_data$unit)

lookup_units <- data.frame(unit = units, unit_name =  c("N_Sulawesi_Gorontalo", 
                                                        "W_S_Sulawesi", "SE_Sulawesi",
                                                        "C_Sulawesi", "C_Maluku", "S_Maluku","N_Maluku", 
                                                        "E_Nusa_Tenggara","W_Nusa_Tenggara"))

# careful if I was to change the order can only change the facvor levels
lookup_units$unit_name <- factor(lookup_units$unit_name, levels = c("N_Sulawesi_Gorontalo", "C_Sulawesi",
                                                                    "W_S_Sulawesi", "SE_Sulawesi",
                                                                    "N_Maluku","C_Maluku", "S_Maluku", 
                                                                    "W_Nusa_Tenggara", "E_Nusa_Tenggara"))

forest_data_1 <- forest_data %>%
  left_join(lookup_units, by = "unit") %>%
  mutate(interval_length = 1) %>%
  filter(year != 2000)





## 5 years
# alternatively sum this up together, means we have one less script and things that can go wrong

# deforestation rate
forest_data_5 <- data %>%
  dplyr::filter(type == "loss"  & unit != "Wallacea") %>%
  mutate(year = year + 2000) %>%
  filter(year != 2000) %>%
  left_join(lookup_units, by = "unit")


# combine loss for everything that is 2001-2013
# combine loss for everything that is 2014-2018
forest_data_5[forest_data_5$year >= 2001& forest_data_5$year <= 2013, "year_interval"] <- "2001-2013"
forest_data_5[forest_data_5$year >= 2014& forest_data_5$year <= 2018, "year_interval"] <- "2014-2018"

forest_data_5 <- forest_data_5 %>%
  group_by(unit_name, year_interval )%>%
  summarize(abs_loss_km2 = sum(as.numeric(nr_pixel)*res2km2)) # check whether this is the solution

forest_data_5[forest_data_5$year_interval == "2001-2013", "interval_length"] <- length(c(2001:2013))
forest_data_5[forest_data_5$year_interval == "2001-2013", "year"] <- mean(c(2001:2013))
forest_data_5[forest_data_5$year_interval == "2014-2018", "interval_length"] <- length(c(2014:2018))
forest_data_5[forest_data_5$year_interval == "2014-2018", "year"] <- mean(c(2014:2018))


#------------------#
# Prepare baseline #
#------------------#

# here I need cover which is missing for year 0

baseline <-  data %>%
  dplyr::filter(type == "forest" & year == 0  & unit != "Wallacea")%>%
  mutate(area_km2 = as.numeric(nr_pixel) * res2km2) %>%
  dplyr::select(unit, forest_area_2000_km2 = area_km2)%>%
  left_join(lookup_units, by= "unit")

endline <- data %>%
  dplyr::filter(type == "forest" & year == 18 &  unit != "Wallacea")%>%
  mutate(area_km2 = as.numeric(nr_pixel) * res2km2) %>%
  dplyr::select(unit, forest_area_2018_km2 = area_km2)%>%
  left_join(lookup_units, by= "unit")


# calculate total observed forest loss

obs_loss <- data %>% 
  dplyr::filter(type == "forest"  &  unit == "Wallacea")
  
obs_loss_2000 <- obs_loss %>% 
  filter(year == 0) %>% 
  mutate(abs_area_km2 =  as.numeric(nr_pixel)*res2km2)

obs_loss_2018 <- obs_loss %>% 
  filter(year == 18) %>% 
  mutate(abs_area_km2 =  as.numeric(nr_pixel)*res2km2)


# obs percent loss
100-(obs_loss_2018$abs_area_km2 *100/obs_loss_2000$abs_area_km2)



# observed loss
obs_loss_2000$abs_area_km2 -obs_loss_2018$abs_area_km2



#------------#
#   FUTURE   #
#------------#
# from: future_forest_loss_5y_I_roads_each_unit.py
forest_data_5_project <- read.csv(file.path(PCVANT_path, "results/deforestation_quantification/future_area_loss_per_unit_i.csv")) %>%
  filter(type == "loss") %>%
  mutate(abs_loss_km2 = as.numeric(nr_pixel)*res2km2)

#test
max(forest_data_5_project$i)

years = unique(forest_data_5_project$year)
units <- unique(forest_data_5_project$unit)

# this is always the end-year--> why not the middle year
forest_data_5_project <- forest_data_5_project %>%
  mutate(year_2 = 2016 + year * 5, # originally 18
         interval_length = 5)

unique(forest_data_5_project$year_2)






# combine both and have a label that we can use to wrap it with (observed - projected)


forest_data_observed <- forest_data_1 %>%
  mutate(yearly_loss_km2 = abs_loss_km2 / interval_length) %>%
  dplyr::select(-unit)%>%
  left_join(baseline, by = "unit_name")%>%
  mutate(perc_loss = yearly_loss_km2*100/as.numeric(forest_area_2000_km2),
         type = "observation")%>%
  dplyr::select(unit_name, year, forest_area_2000_km2, abs_loss_km2, yearly_loss_km2, perc_loss, type ) 



forest_data_5_observed <- forest_data_5 %>%
  mutate(yearly_loss_km2 = abs_loss_km2 / interval_length)%>%
  left_join(baseline, by = "unit_name")%>%
  mutate(type = "observation", 
         perc_loss = yearly_loss_km2*100/forest_area_2000_km2) %>%
  dplyr::select(unit_name, year, forest_area_2000_km2, abs_loss_km2, yearly_loss_km2, perc_loss, type ) 

forest_data_5_project <- forest_data_5_project %>%
  mutate(yearly_loss_km2 = abs_loss_km2 / interval_length,
         unit = as.character(X..unit)) %>%
  left_join(baseline, by = "unit")%>%
  mutate(type = "projection", 
         perc_loss = yearly_loss_km2*100/forest_area_2000_km2) %>%
  dplyr::select(unit_name, year = year_2, abs_loss_km2, yearly_loss_km2, perc_loss, type)


# for the total of wallacea for output
# median annual future deforestation in percent of the 2000 baseline
#OUTPUT-OUTPUT-OUTPUT
forest_data_proj_ci_wallacea<- forest_data_5_project %>%
  filter(year != 2016) %>% 
  summarize(
    median = median(perc_loss),
    lower_ci = quantile(perc_loss, probs = 0.025),
    upper_ci = quantile(perc_loss, probs = 0.975))

# for plotting
forest_data_proj_ci <- forest_data_5_project %>%
  group_by(unit_name, year)%>%
  summarize(
    median = median(perc_loss),
    lower_ci = quantile(perc_loss, probs = 0.025),
    upper_ci = quantile(perc_loss, probs = 0.975))

# this is the overshoot

overshoot<- foreach(unit = lookup_units$unit_name, .combine = 'rbind')%do%{
  divergence <- (forest_data_proj_ci[forest_data_proj_ci$unit_name == unit &
                        forest_data_proj_ci$year ==2016, "median" ]- # changed to middle here
                        forest_data_5_observed[forest_data_5_observed$unit_name == unit &
                           forest_data_5_observed$year == 2016, "perc_loss"])

 return(c(as.character(unit), divergence) )
}

overshoot <- data.frame(unit = lookup_units$unit_name, median = unlist(overshoot[, 2]))
round(min(overshoot$median), 2)
as.character(overshoot[overshoot$median == min(overshoot$median), "unit"])
round(max(overshoot$median), 2)
as.character(overshoot[overshoot$median == max(overshoot$median), "unit"])



#forest_data_observed <- forest_data_observed[forest_data_observed$year != 2000, ]

forest_data <- rbind(forest_data_observed, forest_data_5_observed)


# test whether I can smooth out & check trends of observations with a gam




x = forest_data_observed$year
y = forest_data_observed$perc_loss
# we want a barchart that extends until 2031
# which plots first the five year observation, then the one year on top

# FIGURE PART B
  
# from: future_forest_loss_5y_I_roads_each_unit.py
  
pdf(file.path(PCVANT_path, "results/deforestation_quantification/proj_observ_forest_loss_bar_units.pdf"), paper = "a4r", width = 0, height = 0)
ggplot()+
 geom_bar(data = forest_data_proj_ci,
          aes(x = year, y = median), color = "blue", fill = "white", stat = "identity", position = "dodge")+
 geom_errorbar(data = forest_data_proj_ci,
              aes(ymin=lower_ci, ymax=upper_ci, x = year), size = 0.6, width=1.2)+
  geom_bar(data = forest_data_5_observed,
            aes(x = year, y = perc_loss), color = "red", fill = "grey", alpha = 0.4, size = 0.5, stat = "identity", position = "dodge")+
  geom_point(data = forest_data_observed,
             aes(x = year, y = perc_loss), color = "black")+
  geom_line(data = forest_data_observed,
            aes(x = year, y = perc_loss), color = "red", linetype = "dashed", size = 0.5)+
  facet_wrap(~unit_name) +
 #   scale_x_continuous(breaks = c(2001, 2007, 2016, 2021, 2026, 2031, 2036, 2041,2046,2051 ), 
#           labels = c(2001, 2011, 2016, 2021, 2026, 2031, 2036, 2041,2046,2051))+
  # scale_x_continuous(breaks = c(2001:2053),
  #                    labels = c(2001,rep("", times = 11),2013,2014,rep("", times = 3),
  #                              2018,2019,rep("", times = 3),
  #                              2023,2024, rep("", times = 3), 2028,2029,rep("", times = 3), 2033,2034,
  #                              rep("", times = 3), 2038,2039, rep("", times = 3) ,
  #                              2043,2044, rep("", times = 3),
  #                              2048,2049,rep("", times = 3),  2053))+
  scale_x_continuous(breaks = c(2001:2053),
                   labels = c(2001,rep("", times = 12),2014,rep("", times = 4),
                             2019,rep("", times = 4),
                             2024, rep("", times = 3), 2029,rep("", times = 4), 2034,
                             rep("", times = 4), 2039, rep("", times = 4) ,
                             2044, rep("", times = 4),
                             2049,rep("", times = 4),  2053))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x=element_text(size = 12, angle = 90, hjust=0.5,vjust=0.2),
    # axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 18, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 16, face="bold")) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
dev.off()


#plot for overview figure, annual for whole of Wallacea
# sum everything togther to get wallacea

forest_data_observed_annual <- data %>%
  dplyr::filter(type == "loss"  & unit == "Wallacea" &
  year != 0) %>%
  mutate(year = year + 2000,
         forest_loss_km2 = as.numeric(nr_pixel)*res2km2) 

# i need to add baseline
pdf(file.path(PCVANT_path, "results/deforestation_quantification/proj_observ_forest_loss_Wallacea_anual.pdf"), paper = "a4r", width = 0, height = 0)

ggplot()+
    geom_point(data = forest_data_observed_annual,
             aes(x = year, y = forest_loss_km2), color = "darkgreen")+
  geom_line(data = forest_data_observed_annual,
            aes(x = year, y = forest_loss_km2), color = "darkgreen", linetype = "dashed", size = 0.5)+
 scale_x_continuous(breaks = c(2001:2018),
                    labels = c(2001,rep("", times = 11),2013,
                               2014,rep("", times = 3),2018))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x=element_text(size = 12, angle = 90, hjust=0.5,vjust=0.2),
    # axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 18, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 16, face="bold")) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
dev.off()








# Numbers for manuscript AND FIGURE PART A
# median projected
range(forest_data_proj_ci$median)

# comparison 

# comparison of forest area
baseline <- mutate(baseline, year = 2000)%>%
  rename(forest_area_km2 = forest_area_2000_km2)

endline <- mutate(endline, year = 2018) %>%
  rename(forest_area_km2 = forest_area_2018_km2)

# to know where we are in the future, we need to use forest 2019 and deduce the summed loss for each of the 100
# calculate mean and then per year 

# Figure part A
# this is the upper part of the graph
# add year 33 as well


forest_data_5_project <- read.csv(file.path(PCVANT_path, "results/deforestation_quantification/future_area_loss_per_unit_i.csv")) %>%
  filter(type == "loss") %>%
  rename( unit = X..unit) %>%
  mutate( abs_loss_km2 = nr_pixel*res2km2,
          unit = as.factor(unit))

# is there something wrong with endline?

futureline_33 <- forest_data_5_project %>%
  filter(year != 0 & year < 4) %>% # exclude the first year
  # loss over all years for each iteration
  group_by(unit, i)%>%
  summarize(future_loss_sum_km2 = sum(abs_loss_km2))%>% # here you sum the loss for each year together to the grant total which is 53
  #calculate final area
  left_join(endline, by = "unit") %>%
  dplyr::select(unit, i, future_loss_sum_km2, forest_area_18_km2 = forest_area_km2) %>% 
  left_join(baseline, by = "unit")%>%
  rename(forest_area_00_km2 = forest_area_km2) %>% 
  mutate(final_area_km2 = forest_area_18_km2 - future_loss_sum_km2,
         percentage_area_left_2000 = round((final_area_km2 * 100 / forest_area_00_km2) , 2),
         percentage_area_lost_2000 = round(100-(final_area_km2 * 100 / forest_area_00_km2) , 2),
         percentage_area_loss = round(100- (final_area_km2 * 100 / forest_area_18_km2) , 2))%>%
   group_by(unit, unit_name)%>%
  summarize(
    forest_area_km2 = median(final_area_km2),
    lower_ci_area = quantile(final_area_km2, probs = 0.025),
    upper_ci_area = quantile(final_area_km2, probs = 0.975),
    median_percent_loss = median(percentage_area_loss),
    lower_ci_percent_loss = quantile(percentage_area_loss, probs = 0.025),
    upper_ci_percent_loss = quantile(percentage_area_loss, probs = 0.975),
    median_percent_lost_00 = round(median(percentage_area_lost_2000), 2),
    lower_ci_percent_lost_00 = round(quantile(percentage_area_lost_2000, probs = 0.025), 2),
    upper_ci_percent_lost_00 = round(quantile(percentage_area_lost_2000, probs = 0.975), 2),
    median_percent_left_00 = round(median(percentage_area_left_2000), 2),
    lower_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.025), 2),
    upper_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.975), 2)  )%>%
  mutate(year = 2033) %>% # this is the end of last interval
  as.data.frame()

futureline_33_ci <- futureline_33 %>%
  dplyr::select(unit_name, year, lower_ci_area, upper_ci_area) %>%
  mutate(unit = 0)

futureline_33_ci$unit <- c(1:9)

futureline_33_area <- futureline_33 %>%
  dplyr::select(unit_name, forest_area_km2, year)



futureline_53 <- forest_data_5_project %>%
  filter(year != 0) %>%
  # loss over all years for each iteration
  group_by(unit, i)%>%
  summarize(future_loss_sum_km2 = sum(abs_loss_km2)) %>%  # here you sum the loss for each year together to the grant total which is 53
   #calculate final area
  left_join(endline, by = "unit") %>%
  dplyr::select(unit, i, future_loss_sum_km2, forest_area_18_km2 = forest_area_km2) %>% 
  left_join(baseline, by = "unit")%>%
  rename(forest_area_00_km2 = forest_area_km2) %>% 
  mutate(final_area_km2 = forest_area_18_km2 - future_loss_sum_km2,
         percentage_area_left_2000 = round((final_area_km2 * 100 / forest_area_00_km2) , 2),
         percentage_area_lost_2000 = round(100-(final_area_km2 * 100 / forest_area_00_km2) , 2),
         percentage_area_loss = round(100- (final_area_km2 * 100 / forest_area_18_km2) , 2),
         yearly_loss_total = calc.annual.change(forest_area_00_km2 , final_area_km2, 2000, 2053),
         yearly_loss_observed = calc.annual.change(forest_area_00_km2 , forest_area_18_km2, 2000, 2018),
         yearly_loss_future = calc.annual.change(forest_area_18_km2, final_area_km2, 2018, 2053))%>%
  group_by(unit, unit_name)%>%
  summarize(
    forest_area_km2 = median(final_area_km2),
    lower_ci_area = quantile(final_area_km2, probs = 0.025),
    upper_ci_area = quantile(final_area_km2, probs = 0.975),
    median_percent_loss = median(percentage_area_loss),
    lower_ci_percent_loss = quantile(percentage_area_loss, probs = 0.025),
    upper_ci_percent_loss = quantile(percentage_area_loss, probs = 0.975),
    median_percent_lost_00 = round(median(percentage_area_lost_2000), 2),
    lower_ci_percent_lost_00 = round(quantile(percentage_area_lost_2000, probs = 0.025), 2),
    upper_ci_percent_lost_00 = round(quantile(percentage_area_lost_2000, probs = 0.975), 2),
    median_percent_left_00 = round(median(percentage_area_left_2000), 2),
    lower_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.025), 2),
    upper_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.975), 2),
    median_percent_yearly_loss_total = round(median(yearly_loss_total), 2),
    lower_ci_percent_yearly_loss_total = round(quantile(yearly_loss_total, probs = 0.025), 2),
    upper_ci_percent_yearly_loss_total = round(quantile(yearly_loss_total, probs = 0.975), 2),
    
    median_percent_yearly_loss_obs = round(median(yearly_loss_observed), 2),
    median_percent_yearly_loss_future = round(median(yearly_loss_future), 2),
    lower_ci_percent_yearly_loss_future = round(quantile(yearly_loss_future, probs = 0.025), 2),
    upper_ci_percent_yearly_loss_future = round(quantile(yearly_loss_future, probs = 0.975), 2)
    
  
      )%>%
  mutate(year = 2053) %>% # this is the end of last interval
  as.data.frame()

futureline_53_ci <- futureline_53 %>%
  dplyr::select(unit_name, year, lower_ci_area, upper_ci_area) %>%
  mutate(unit = 0)

futureline_53_ci$unit <- c(1:9)

futureline_53_area <- futureline_53 %>%
  dplyr::select(unit_name, forest_area_km2, year)


#here

forest_area <- rbind(baseline, endline)%>%
  dplyr::select(unit_name, forest_area_km2, year)%>%
  rbind(., futureline_33_area, futureline_53_area) %>%
  melt(id = c("unit_name", "year"))

futureline_ci <- futureline_33_ci %>%
  rbind(., futureline_53_ci)


# forest_area$value <- forest_area$value / 10000
# futureline_33_ci$lower_ci_area <- futureline_33_ci$lower_ci_area / 10000
# futureline_33_ci$upper_ci_area <- futureline_33_ci$upper_ci_area / 10000
# 
# futureline_53_ci$lower_ci_area <- futureline_53_ci$lower_ci_area / 10000
# futureline_53_ci$upper_ci_area <- futureline_53_ci$upper_ci_area / 10000

# prep for plotting
pdf(file.path(PCVANT_path, "results/deforestation_quantification/proj_observ_forest_unit_unique.pdf"), paper = "a4r", width = 0, height = 0)
ggplot()+
  # here!!!
  geom_bar(data = forest_area,
           aes(x = as.factor(year), y = value, fill = as.factor(year)), 
           # fill = alpha(rgb(0,97,0, max = 255), 0.5), 
           stat = "identity", width = 0.9)+ #, position = position_dodge(width = 0.2))+
  geom_errorbar(data = futureline_33_ci, # this I might have to tweak
                aes(ymin=lower_ci_area, ymax=upper_ci_area, x =as.factor(year)),
                size = 1, width=0.16)+
  geom_errorbar(data = futureline_53_ci, # this I might have to tweak
                aes(ymin=lower_ci_area, ymax=upper_ci_area, x = as.factor(year)),
                size = 1, width=0.16)+
  facet_wrap(~unit_name, scale = "free") +
  ylim(0, 41000)+
  
  scale_fill_manual(values=alpha(c("darkgreen", "darkgreen", "green", "green"), 0.5)) +
  
  #scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  theme_bw() +
  xlab("year") +
  ylab(as.expression(bquote("forest cover in " ~ km^2 ~ ""))) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size = 14),
 #   axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 18, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 16, face="bold")) +
  theme(plot.margin = unit(c(0.2,-0.2,0.1,0.1),"cm"))
dev.off()

# would like toa add percentage loss
pdf(file.path(PCVANT_path, "results/deforestation_quantification/proj_observ_forest_unit_unique_together.pdf"), paper = "a4r", width = 0, height = 0)
ggplot()+
  # here!!!
  geom_bar(data = forest_area,
           aes(x = as.factor(unit_name), y = value, fill = as.factor(year)), 
           # fill = alpha(rgb(0,97,0, max = 255), 0.5), 
           stat = "identity", width = 0.9, position = position_dodge(width = 1))+
  geom_label(data = futureline_53,
             aes(x = as.factor(unit_name), y = 3000, label=round(median_percent_lost_00,0)))+
   scale_fill_manual(values=alpha(c("darkgreen", "darkgreen", "green", "green"), 0.5)) +
  
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  theme_bw() +
  xlab("year") +
  ylab(as.expression(bquote("forest cover in " ~ km^2 ~ ""))) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size = 14),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title=element_text(size = 18, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 16, face="bold")) +
  theme(plot.margin = unit(c(0.2,-0.2,0.1,0.1),"cm"))
dev.off()


# add the circle graph for poster
forest_area_perc <- rbind(baseline, endline)%>%
  dplyr::select(unit_name, forest_area_km2, year)%>%
  rbind(., futureline_33_area, futureline_53_area) %>%
  pivot_wider(id_cols = unit_name, names_from = year, names_prefix = "year_", values_from = forest_area_km2) %>% 
  mutate(perc_53_loss = (year_2033-year_2053)*100/year_2000,
         perc_33_loss = (year_2018-year_2033)*100/year_2000,
         perc_18_loss = (year_2000-year_2018)*100/year_2000,
         perc_00_rem = 100-perc_53_loss - perc_33_loss - perc_18_loss) %>% 
  pivot_longer(cols = starts_with("perc_"), names_to = "year", values_to = "percent")



# would like toa add percentage loss
hsize = 2
forest_area_perc$year <- as.factor(forest_area_perc$year)
str(forest_area_perc$year)

forest_area_perc$year <- sizes <- factor(forest_area_perc$year , levels = rev(levels(forest_area_perc$year)))
colors <- c("267300ff", "269a00ff", "94cc46ff", "dced36ff")


pdf(file.path(PCVANT_path, "results/deforestation_quantification/proj_observ_forest_unit_poster.pdf"), paper = "a4r", width = 0, height = 0)

ggplot()+
  geom_col(data = forest_area_perc,
           aes(x = hsize, y =percent, fill = year))+
  xlim(c(0.2, hsize + 0.5))+
  
  coord_polar(theta = "y") +
  facet_wrap( ~unit_name)+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size = 14),
    axis.text.x = element_blank(),
    axis.title=element_text(size = 18, face = "bold"),
    #legend.title=element_text(size=18, face = "bold") ,
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    legend.key.height=unit(0.8,"cm"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 16, face="bold")) +
  theme(plot.margin = unit(c(0.2,-0.2,0.1,0.1),"cm"))
dev.off()






loss_C_SLW_18_53 <- forest_area[forest_area$year == 2018 & forest_area$unit_name == "C_Sulawesi","value" ]-
  forest_area[forest_area$year == 2053 & forest_area$unit_name == "C_Sulawesi","value" ] 

round(loss_C_SLW_18_53)

loss_C_SLW_00_53 <- forest_area[forest_area$year == 2000 & forest_area$unit_name == "C_Sulawesi","value" ]-
  forest_area[forest_area$year == 2053 & forest_area$unit_name == "C_Sulawesi","value" ] 
round(loss_C_SLW_00_53)
# how much % of the area of a state is forested in 2000
# load unit lookup, merge with unit lookup containing names
# calculate percents

unit_area <- read.csv(file.path(PCVANT_path, "results/deforestation_quantification/unit_area_lookup.csv"), stringsAsFactors = F) %>%
  mutate(unit = as.factor(unit),
         area_km2 = nr_pixel*res2km2,
         year = 2000, 
         type = "area") %>%
  left_join(lookup_units, by ="unit")  %>%
  dplyr::select(unit, forest_area_km2 = area_km2, unit_name, year, type)



forest_area <- bind_rows(baseline, endline)%>%
  bind_rows(.,futureline_33_area)%>%
  bind_rows(., futureline_53_area)%>%
  mutate(type = "forest")%>%
  bind_rows(., unit_area)%>%
  mutate(id = paste0(type, "_", year)) %>%
  dplyr::select(unit_name, id, forest_area_km2) 




# construct a proper table around this
# unit, area, forest_cover in 2000, percent forest cover in 2000 of area, forest in 2018, 
# forest in 2032, percent forest in 31 of 2018 (loss from 18 - 33/53)

# format futurelines
futureline_33 <- futureline_33 %>%
  dplyr::select(-unit_name, -year) 
names(futureline_33)[2:length(futureline_33)] <- paste0(names(futureline_33)[2:length(futureline_33)], "_2033") 

futureline_53 <- futureline_53 %>%
  dplyr::select(-unit_name, -year) 
names(futureline_53)[2:length(futureline_53)] <- paste0(names(futureline_53)[2:length(futureline_53)], "_2053") 



unit_area <- read.csv(file.path(PCVANT_path, "results/deforestation_quantification/unit_area_lookup.csv"), stringsAsFactors = F) %>%
  mutate(area_km2 = nr_pixel * res2km2,
         unit = as.factor(unit))%>%
  dplyr::select(unit, area_km2) %>%
  left_join(baseline[ , c("unit", "forest_area_km2")], by = "unit")%>%
  rename(forest_area_km2_2000 = forest_area_km2) %>%
  left_join(endline[ , c("unit", "forest_area_km2")], by = "unit")%>%
  rename(forest_area_km2_2018 = forest_area_km2)%>%
  left_join(lookup_units, by = "unit")%>%
  # how much left in 2018 in %
  mutate(percent_left_00_2018 = round(forest_area_km2_2018 * 100/forest_area_km2_2000,2)) %>% 
    # no unit in futureline
  left_join(futureline_33, by = "unit")%>%
  left_join(futureline_53, by = "unit")%>%
  # % forest cover in 2000 from total area
  mutate(perc_cover_with_forest_2000 = round(forest_area_km2_2000*100/area_km2, 2),
         perc_cover_with_forest_2018 = round(forest_area_km2_2018*100/area_km2, 2))
# % forest cover loss into future (based on all 100 values)
# here you should actually use the raw values to calculate the percent decline 

# ggplot()+
#   geom_bar(data = forest_area_km2, 
#            aes(x = unit_name, y = area, fill = as.factor(id)), stat = "identity", position = "dodge") + 
#   geom_label(data = unit_area, aes(x = unit_name, y =100000, label = round(median_percent_loss, 0)))
# 
# round(unit_area$median_percent_loss - unit_area$lower_ci_percent_loss,0)
# round(unit_area$upper_ci_percent_loss - unit_area$median_percent_loss, 0)
# 
# export in decimal separated notation
unit_area_out <- unit_area %>%
  dplyr::select(unit_name, area_km2:forest_area_km2_2000,perc_cover_with_forest_2000,forest_area_km2_2018,
                perc_cover_with_forest_2018, 
                forest_area_km2_2033:upper_ci_percent_loss_2053,
                median_percent_yearly_loss_total_2053:upper_ci_percent_yearly_loss_future_2053) 

write.csv(unit_area_out , file.path(PCVANT_path, "results/deforestation_quantification/unit_area_out.csv"), row.names = F)



# reformat to nicer format for table
for (name in names(unit_area_out)[c(2,3, 5,7:9, 19:21)]){
  print(name)
  unit_area_out[, name] <- as.character(format(round(unit_area_out[ , name], 0), big.mark = ","))
}

for (name in names(unit_area_out)[c(4, 6, 10:18, 22:24)]){
  unit_area_out[, name] <- as.character(round(unit_area_out[ , name],0))
}

for (name in names(unit_area_out)[c(25:length(unit_area_out))]){
  print(name)
  unit_area_out[, name] <- as.character(round(unit_area_out[ , name],2))
}

# remember observed doesnt have CI!!!


plot_area_out <- unit_area_out %>% 
  arrange(unit_name) %>% 
  pivot_longer(cols = area_km2:upper_ci_percent_yearly_loss_future_2053) %>% 
  pivot_wider(names_from = unit_name, values_from = value) %>% 
  mutate(year = substr(name, start =nchar(name)-3, stop = nchar(name))) %>% 
  dplyr::select(year, name,  N_Sulawesi_Gorontalo:E_Nusa_Tenggara) %>% 
  as.data.frame()

name_new <- c()
for (i in 2:length(plot_area_out$name)){
  name <- substr(plot_area_out$name[i], start= 1, 
         stop = nchar(plot_area_out$name[i])-5)
  name_new <- c(name_new, name)
}

plot_area_out$name <- c("area_km2", name_new)

plot_area_out[plot_area_out$name == "area_km2", "year"] <- "-"

  
write.csv(plot_area_out, file.path(PCVANT_path, "results/deforestation_quantification/overview_table_fig_3.csv"), row.names = F)


# text
baseline_wallacea <- sum(baseline$forest_area_km2)
endline_wallacea <- sum(endline$forest_area_km2)

# forest area in 2000, in 18 as percent of 2000, and 33 and 53  each as % of 2000

perc_left_18 <- endline_wallacea * 100/ baseline_wallacea

# take all out of unit_area
futureline_wallacea_33 <- forest_data_5_project %>%
  filter(year != 0 & year < 4) %>% # exclude the first year
group_by( i)%>%
  summarize(future_loss_sum_km2 = sum(abs_loss_km2))%>% # here you sum the loss for each year together to the grant total which is 53
  #calculate final area
  mutate(final_area_km2 = endline_wallacea - future_loss_sum_km2,
         percentage_area_left_2000 = round((final_area_km2 * 100 / baseline_wallacea) , 2)) %>% 
  summarize(
    forest_area_km2 = median(final_area_km2),
    lower_ci_area = quantile(final_area_km2, probs = 0.025),
    upper_ci_area = quantile(final_area_km2, probs = 0.975),
    median_percent_left_00 = round(median(percentage_area_left_2000), 2),
    lower_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.025), 2),
    upper_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.975), 2)  )%>%
  mutate(year = 2033) %>% # this is the end of last interval
  as.data.frame()

  
# take all out of unit_area
futureline_wallacea_53 <- forest_data_5_project %>%
  filter(year != 0) %>%
  group_by( i)%>%
  summarize(future_loss_sum_km2 = sum(abs_loss_km2))%>% # here you sum the loss for each year together to the grant total which is 53
  #calculate final area
  mutate(final_area_km2 = endline_wallacea - future_loss_sum_km2,
         percentage_area_left_2000 = round((final_area_km2 * 100 / baseline_wallacea) , 2),
         percentage_area_loss_2018 = round(100- (final_area_km2 * 100 / endline_wallacea) , 2),
         yearly_loss_total = calc.annual.change(baseline_wallacea, final_area_km2, 2000, 2053),
         yearly_loss_future = calc.annual.change(endline_wallacea, final_area_km2, 2018, 2053))%>%
  summarize(
    forest_area_km2 = median(final_area_km2),
    lower_ci_area = quantile(final_area_km2, probs = 0.025),
    upper_ci_area = quantile(final_area_km2, probs = 0.975),
    median_percent_left_00 = round(median(percentage_area_left_2000), 2),
    lower_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.025), 2),
    upper_ci_percent_left_00 = round(quantile(percentage_area_left_2000, probs = 0.975), 2) ,
    median_percent_loss_18 = round(median(percentage_area_loss_2018), 2),
    lower_ci_percent_loss_18 = round(quantile(percentage_area_loss_2018, probs = 0.025), 2),
    upper_ci_percent_loss_18 = round(quantile(percentage_area_loss_2018, probs = 0.975), 2) ,
    median_yearly_loss_total = round(median(yearly_loss_total), 2),
    lower_ci_yearly_loss_total = round(quantile(yearly_loss_total, probs = 0.025), 2),
    upper_ci_yearly_loss_total = round(quantile(yearly_loss_total, probs = 0.975), 2),
    median_yearly_loss_future = round(median(yearly_loss_future), 2),
    lower_ci_yearly_loss_future = round(quantile(yearly_loss_future, probs = 0.025), 2),
    upper_ci_yearly_loss_future = round(quantile(yearly_loss_future, probs = 0.975), 2),
    )%>%
  mutate(year = 2053) %>% # this is the end of last interval
  as.data.frame()
  
total_loss_00_18 <-  baseline_wallacea - endline_wallacea
total_loss_00_53 <-  baseline_wallacea - futureline_wallacea_53$forest_area_km2
total_loss_18_53 <-  round(endline_wallacea - futureline_wallacea_53$forest_area_km2)


# text
paste(" Until 2018 Wallacean forests had decreased by ",round(perc_left_18,0)  ,"% of its 2000 cover. Future deforestation between 2018", 
      "and 2033 was projected to further decrease forest cover to ", round(futureline_wallacea_33$median_percent_left_00,0) ,
      "%, and to ",round(futureline_wallacea_53$median_percent_left_00, 0), "% by 2053.")

paste("Total observed and projected forest loss between 2000 and 2053 amounted to ", round(total_loss_00_53, 0)," km2.", 
      "Median regional deforestation ranged from ", min(unit_area$median_percent_left_00_2053), "(",
      unit_area[unit_area$median_percent_left_00_2053 == min(unit_area$median_percent_left_00_2053), "unit_name"], ") to ", 
      max(unit_area$median_percent_left_00_2053), "(",
      unit_area[unit_area$median_percent_left_00_2053 == max(unit_area$median_percent_left_00_2053), "unit_name"], 
      ")") 

paste("Total observed and projected forest loss between 2000 and 2053 amounted to ", round(total_loss_00_53, 0)," km2.", 
      "Median regional deforestation ranged from ", round(min(unit_area$median_percent_left_00_2053)), "(",
      unit_area[unit_area$median_percent_left_00_2053 == min(unit_area$median_percent_left_00_2053), "unit_name"], ") to ", 
     round(max(unit_area$median_percent_left_00_2053)), "(",
      unit_area[unit_area$median_percent_left_00_2053 == max(unit_area$median_percent_left_00_2053), "unit_name"], 
      ")")
      

paste("Total loss between 2000 and 2018 amounted to ", round(total_loss_00_18, 0)," km2.")

# forest in 2030 as % of forest in 2000

paste("projected future deforestation in Wallacea [...] ",
      round(100-futureline_wallacea_33$median_percent_left_00),
            "% for Wallacea, compared to 2000.")




observed_annual_loss <- calc.annual.change(baseline_wallacea, endline_wallacea, 2000, 2018)
futureline_wallacea_53$median_yearly_loss_future

paste0("The observed annual losses of forest across Wallacea between 2000 and 2018 was ", observed_annual_loss,
       "%, which was projected to increase to ", futureline_wallacea_53$median_yearly_loss_future, "%")

       
max(futureline_53$median_percent_yearly_loss_future_2053)
futureline_53[futureline_53$median_percent_yearly_loss_future_2053 == max(futureline_53$median_percent_yearly_loss_future_2053), "unit" ]
#N MLK

min(futureline_53$median_percent_yearly_loss_future_2053)
futureline_53[futureline_53$median_percent_yearly_loss_future_2053 == max(futureline_53$median_percent_yearly_loss_future_2053), "unit" ]
#East NT

# borneo rate from reconciling forest conservation
forest_start_borneo <- 303525
forest_end_borneo <- 303525 - 14212 
start_year_Borneo = 2000
end_year_Borneo = 2010

annual_loss_Borneo = calc.annual.change(forest_start_borneo, forest_end_borneo , 2000, 2010)


# getting the data from gaveaus 21 paper

defor_overview <- read.csv("C:/Users/mv296/work/Indonesia/deforestation/gaveau_deforestation_table.csv", stringsAsFactors = F)

forest_Kali_2019 <- defor_overview[defor_overview$area_ha == "2019 forest area", "Kalimantan"]
loss_Kali_01_19 <- defor_overview[defor_overview$area_ha == "Forest loss 2001-2019", "Kalimantan"]
forest_Kali_2000 <- forest_Kali_2019 + loss_Kali_01_19
Kali_change <- calc.annual.change(forest_Kali_2000, forest_Kali_2019, 2000, 2019, round_to = 2)
round(Kali_change, 2)
#0.76%

forest_Sumatra_2019 <- defor_overview[defor_overview$area_ha == "2019 forest area", "Sumatra"]
loss_Sumatra_01_19 <- defor_overview[defor_overview$area_ha == "Forest loss 2001-2019", "Sumatra"]
forest_Sumatra_2000 <- forest_Sumatra_2019 + loss_Sumatra_01_19
Sumatra_change <- calc.annual.change(forest_Sumatra_2000, forest_Sumatra_2019, 2000, 2019, round_to = 2)
round(Sumatra_change,2)
# 1.52%

# of 0.39% in Wallacea not including 2019 loss
Wallacea_change <- observed_annual_loss


round(Kali_change/Wallacea_change,0)

round(Wallacea_change/Kali_change,2)

round(Sumatra_change/Wallacea_change, 0)
round(Wallacea_change/Sumatra_change, 2)
                                           
# old text
# observed forest
max_area <- max(unit_area$forest_area_km2_2018)
max_area_province <- unit_area[unit_area$forest_area_km2_2018== max_area , "unit_name"]
max_area_perc <-  unit_area[unit_area$forest_area_km2_2018== max_area , "perc_cover_with_forest_2018"]

min_area <- min(unit_area$forest_area_km2_2018)
min_area_province <- unit_area[unit_area$forest_area_km2_2018== min_area , "unit_name"]
min_area_perc <-  unit_area[unit_area$forest_area_km2_2018== min_area , "perc_cover_with_forest_2018"]

  
paste0("Total forest area at the end of the observation period (2018) varied across provinces ", 
"(Fig. 3 and Table S3) with ", max_area_province ," having the largest area (", max_area ," km2) in 2018, representing ", 
max_area_perc, "%  of the total province area, and ", min_area_province ," having the smallest area (", 
min_area ,"km2, ", min_area_perc, "% of its area).")

  # (95%CI: ", xx, "-", xx, "%)
  
# loss between 33 and 53
min_forest_perc_loss_33 <- min(unit_area$median_percent_loss_2033) 
min_forest_perc_loss_33_lci <-  unit_area[unit_area$median_percent_loss_2033== min_forest_perc_loss_33 , "lower_ci_percent_loss_2033"]
min_forest_perc_loss_33_uci <-  unit_area[unit_area$median_percent_loss_2033== min_forest_perc_loss_33 , "upper_ci_percent_loss_2033"]

min_forest_perc_loss_33_region <-  unit_area[unit_area$median_percent_loss_2033== min_forest_perc_loss_33 , "unit_name"]

max_forest_perc_loss_33 <- max(unit_area$median_percent_loss_2033) 
max_forest_perc_loss_33_lci <-  unit_area[unit_area$median_percent_loss_2033== max_forest_perc_loss_33 , "lower_ci_percent_loss_2033"]
max_forest_perc_loss_33_uci <-  unit_area[unit_area$median_percent_loss_2033== max_forest_perc_loss_33 , "upper_ci_percent_loss_2033"]

max_forest_perc_loss_33_region <-  unit_area[unit_area$median_percent_loss_2033== max_forest_perc_loss_33 , "unit_name"]



min_forest_perc_loss_53 <- min(unit_area$median_percent_loss_2053) 
min_forest_perc_loss_53_lci <-  unit_area[unit_area$median_percent_loss_2053== min_forest_perc_loss_53 , "lower_ci_percent_loss_2033"]
min_forest_perc_loss_53_uci <-  unit_area[unit_area$median_percent_loss_2053== min_forest_perc_loss_53 , "upper_ci_percent_loss_2033"]

min_forest_perc_loss_53_region <-  unit_area[unit_area$median_percent_loss_2053== min_forest_perc_loss_53 , "unit_name"]

max_forest_perc_loss_53 <- max(unit_area$median_percent_loss_2053) 
max_forest_perc_loss_53_lci <-  unit_area[unit_area$median_percent_loss_2053== max_forest_perc_loss_53 , "lower_ci_percent_loss_2053"]
max_forest_perc_loss_53_uci <-  unit_area[unit_area$median_percent_loss_2053== max_forest_perc_loss_53 , "upper_ci_percent_loss_2053"]

max_forest_perc_loss_53_region <-  unit_area[unit_area$median_percent_loss_2053== max_forest_perc_loss_53 , "unit_name"]

#CONTINUE HERE

paste0("In all provinces natural forest decreased in extent between the start ",
"of the observation period in 2000 and 2018, and was projected to decrease further until 2033, ",
"ranging from ", min_forest_perc_loss_33 , "% (CI: ", round(min_forest_perc_loss_33_lci, 2), "-",
       round(min_forest_perc_loss_33_uci,2) , ", ", min_forest_perc_loss_33_region , "), ", " to ",
       max_forest_perc_loss_33, "% (CI: ", round(max_forest_perc_loss_33_lci,2) ,"-",
      round(max_forest_perc_loss_33_uci,2) ,
       ", ",max_forest_perc_loss_33_region, ") reduction (Fig. 2b and Supporting Information).", 
       " Until 2053 forest was projected to decrease further from ", 
min_forest_perc_loss_53 , "% (CI: ", round(min_forest_perc_loss_53_lci,2), "-",
round(min_forest_perc_loss_53_uci,2) , ", ", min_forest_perc_loss_53_region , ") to ",
max_forest_perc_loss_53, "% (CI: ", round(max_forest_perc_loss_53_lci,2) ,"-",
round(max_forest_perc_loss_53_uci,2) ,
", ",max_forest_perc_loss_53_region, ")")


#forest_data_5_observed // forest_data_observed // forest_data_proj_ci

min_observed_change <- min(forest_data_observed$perc_loss)
min_observed_change_unit <- forest_data_observed[forest_data_observed$perc_loss == min_observed_change, "unit_name"]

max_observed_change <- max(forest_data_observed$perc_loss)
max_observed_change_unit <- forest_data_observed[forest_data_observed$perc_loss == max_observed_change, "unit_name"]

forest_data_E_NT <- forest_data_5_observed %>%
  filter(unit_name == "E_Nusa_Tenggara")
 
km_difference_E_NT <- forest_data_E_NT[forest_data_E_NT$year == 2007, "yearly_loss_km2"] - forest_data_E_NT[forest_data_E_NT$year == 2016, "yearly_loss_km2"] 
perc_difference_E_NT <- forest_data_E_NT[forest_data_E_NT$year == 2007, "perc_loss"] - forest_data_E_NT[forest_data_E_NT$year == 2016, "perc_loss"] 
#   
forest_data_proj_ci_out <- forest_data_proj_ci %>%
  filter(year != 2016) %>%
  as.data.frame()


perc_min_unit <- min(forest_data_proj_ci_out$median)
perc_min_unit_region <- forest_data_proj_ci_out[forest_data_proj_ci_out$median == perc_min_unit, "unit_name"]
perc_min_unit_median <- median(forest_data_proj_ci_out[forest_data_proj_ci_out$median == perc_min_unit, "median"])

perc_max_unit <- max(forest_data_proj_ci_out$median)
perc_max_unit_region <- as.character(forest_data_proj_ci_out[forest_data_proj_ci_out$median == perc_max_unit, "unit_name"])
perc_max_unit_year <- forest_data_proj_ci_out[forest_data_proj_ci_out$median == perc_max_unit, 
                                                           "year"] 

perc_max_unit_year+3

# hard coding the smallest and largest change between rates in regions
# med_perc_change
forest_data_W_NT <- forest_data_proj_ci[forest_data_proj_ci$unit_name == "W_Nusa_Tenggara" ,  ]
change_W_NT <- round(abs(min(forest_data_W_NT$median)-max(forest_data_W_NT$median)), 2)

forest_data_E_NT <- forest_data_proj_ci[forest_data_proj_ci$unit_name == "E_Nusa_Tenggara" ,  ]
change_E_NT <- round(abs(min(forest_data_E_NT$median)-max(forest_data_E_NT$median)), 2)

if (change_W_NT < change_E_NT){
  min_change_unit <- "West Nusa Tenggara"
  min_change <- change_W_NT 
}else{
  min_change_unit <- "East Nusa Tenggara"
  min_change <- change_E_NT 
}

# max change
forest_data_N_MLK <- forest_data_proj_ci[forest_data_proj_ci$unit_name == "N_Maluku" ,  ]
change_N_MLK <- round(abs(min(forest_data_N_MLK$median)-max(forest_data_N_MLK$median)), 2)

paste0("Observed annual deforestation (2000-2018) ranged between ", round(min_observed_change, 2), " (",min_observed_change_unit ,
")-", round(max_observed_change,2), "% (",max_observed_change_unit ,
"), relative to the forest at the beginning of", 
"the observation period (2000), with high inter-annual fluctuations (Fig. 3b).",
"On average, the rate increased from the first (i.e., deforestation before the calibration period) 12 years",
" to the calibration period in all provinces, with the exception of East Nusa Tenggara where the deforestation rates",
"were relatively low and decreased slightly (", round(km_difference_E_NT,2),
"km2 or ",round(perc_difference_E_NT,2) , " percent point difference).",
"Projected median annual deforestation rates for regions ranged between ",round(perc_min_unit,2) ,
"% (",perc_min_unit_region , ") and ", round(perc_max_unit,2),"% (", perc_max_unit_region ,").",
" Of all regions ",forest_data_proj_ci_out[forest_data_proj_ci_out$median == max(forest_data_proj_ci_out$median), "unit_name"]," has the steepest increase to highest levels with highest levels in the",
 " period ",perc_max_unit_year-2," to ",perc_max_unit_year+3,".",
"Over the projection period, the increase in projected median forest loss increased and then decreased",
"in all regions except in South Maluku where it continues to increase until the last step, ",
"although starting at a lower level than the other provinces, and East Nusa Tenggara, where ",
"it stays relatively low (median percent change over projected years ",round(perc_min_unit_median,2) , 
"). ")

paste0("Both modelling regions in Nusa Tenggara have overall low levels of change in deforestation levels", 
       " compared to the other regions (lowest difference in median percent change over projected years of ",
       min_change, " percent points in ", 
       min_change_unit, ".", 
       "The increase was steepest in North Maluku (",change_N_MLK  ," percent points), with highest levels in the period ",
                                        perc_max_unit_year-2," to ",perc_max_unit_year+3,".") 


# 
