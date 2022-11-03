rm(list = ls())
library(tidyverse)
library(dplyr)
library(foreach)
library("scales")   

options(scipen=999)
# Script to plot the rates of deforestation 


M_path <- "M:/Sumatra_model_August22"
forest_layers <- c("hansen", "tmf")


list_units_tmf = rbind(c(1, "Aceh", "Aceh", "A"),
                       c(2, "N_Sumatra", "North Sumatra", "B"))


list_units_hansen = rbind(c(1, "Aceh", "Aceh", "A"),
                          c(2, "N_Sumatra", "North Sumatra", "G"))



for(forest_layer in forest_layers){
  print(forest_layer)
pred_names_inpath <- file.path("C:/Users/mv296/OneDrive - University of Kent/Sumatra/deforstation_model/predictors_final", forest_layer )

out_path <- file.path("C:/Users/mv296/OneDrive - University of Kent/Sumatra/deforstation_model/results")
# loop over units

list_units <- get(paste0("list_units_", forest_layer))

list_units <- as.data.frame(list_units)

list_units <- as.data.frame(list_units)
names(list_units) <- c("unit", "unit_name", "unit_long", "model_run")


list_units$unit_name <- factor(list_units$unit_name, levels = c("Aceh",
                                                                "N_Sumatra"))

lookup_names <- read_csv( file.path(pred_names_inpath, 
                                    paste0("lookup_pred_names_", forest_layer, ".csv")))



# i will import a rooster what the predictor names mean for each model




parameters_m_s_units <- foreach(i = 1:nrow(list_units), .combine = rbind)%do%{
  print(i)
  unit_number <- list_units$unit[i]
  unit_name <- list_units$unit_name[i]
  model_run <- list_units$model_run[i]
  lookup_names_unit <- lookup_names[ , c("variables")]
  lookup_names_unit <- cbind(lookup_names_unit, lookup_names[ , model_run])
  
  names(lookup_names_unit) <- c("pred_name", "variable")
  lookup_names_unit <- lookup_names_unit[!is.na(lookup_names_unit$variable), ]

  in_path_unit <- file.path(M_path, forest_layer, unit_name, paste0("model_", unit_name,"_", model_run ))
    #Script to plot the Parameter Values
  parameters<-read.table(file.path(in_path_unit, "Parameters_Values.txt"), header=T)
  parameters_m <- parameters %>%
    dplyr::select(-Iteration) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value")%>%
    left_join(lookup_names_unit, by ="variable") 
  
  
  for (variable in lookup_names_unit$variable){
    if (sum(parameters_m[parameters_m$variable == variable, "value"]) == 0){
      parameters_m[parameters_m$variable == variable, "value"] <- NA
    }
  }
  
  
  # the measure would be the distance of the median to zero and then plot this with a 
  # a grey scale
  
  
  
  parameters_m_s <- parameters_m %>%
    group_by(variable) %>%
    summarize(value = median(value, na.rm =T),
              value_lower_ci = quantile(value, probs = 0.025, na.rm =T),
              value_upper_ci = quantile(value, probs = 0.975, na.rm =T),
              #return to unique(name)
              name = unique(pred_name))%>%
    as.data.frame()
  

  parameters_m_s$unit <- unit_number

  # this needs to be exported and combined


    return(parameters_m_s)
}

param_levels <- unique(parameters_m_s_units$name)

#???????
#parameters_m_s_units$name <-factor(parameters_m_s$name, levels = rev(param_levels))




parameters_m_s_units <- parameters_m_s_units %>%
  left_join(list_units, by="unit") %>%
  mutate(pred_names_char = as.character(name))
  
parameters_m_s_units$unit_long <- as.factor(parameters_m_s_units$unit_long)
parameters_m_s_units$unit_long <-factor(parameters_m_s_units$unit_long, levels =  c("Aceh",
                                                                                      "North Sumatra")) 




# color scheme differentiating the states and units
# and maybe also shapes
# hex digits from here: https://statisticsglobe.com/identify-default-color-palette-names-of-ggplot2-in-r
cbPalette_custom <- c("#F8766D", "#F37B59", "#ED8141", "#E08B00",
                      "#00BFC4", "#00A5FF", "#7997F", 
                      "#BF80FF", "#FC717F")# 
#
hex_codes2 <- hue_pal()(50)
cbPalette_custom <- hex_codes2[c(1,
                                #  3,5,9,
                                26
                              #  ,33,35,
                                # 38, 
                                #47
                                )]

shape_custom <- c(rep(19, times = 4), rep(15, times = 3), rep(17, times = 2)) # 16 or 19 as circle


parameters_m_s_units$pred_names_char <- factor(parameters_m_s_units$pred_names_char, levels =
                                             rev(c( "intercept",
                                                    "past deforestation",
                                                    "slope",  
                                                    "fire (yearly average)",
                                                    "access (hrs)",
                                                    "roads",
                                                    "rivers",
                                                    "population pressure",
                                                    "subsistence livelihood",
                                                    "plantation livelihood",
                                                    "non-agricultural livelihood",
                                                    "small scale plantations",
                                                    "industrial scale plantation distance",
                                                    "socio-economic deprivations",
                                                    "transmigrant population",
                                                    "peat*",
                                                      "mining (exploration)*",
                                                      "mining (production)*",
                                                      "social forest (implemented)*",
                                                      "social forest (proposed)*",
                                                      "strict protected (IUCN 1 and 2 in WDPA)*",
                                                      "other protected(IUCN further down and national)*",
                                                      "watershed protection forest*",
                                                      "production forest*",
                                                      "conversion forest*"  )))

parameters_m_s_units_filter <- parameters_m_s_units[!is.na(parameters_m_s_units$value),]



# add the number of models it appears in for each predictor
# 
# table to see predictor numbers
parameters_m_s_units_filter_overview <- parameters_m_s_units_filter %>%
  dplyr::select(pred_names_char, unit_long, value) %>%
  group_by(pred_names_char) %>%
  summarize(n = n(),
            mean_coef = round(mean(value), 4),
            min_coef = round(min(abs(value)), 4),
            max_coef = round(max(abs(value)), 4),
            spread =round( max_coef + min_coef), 4)%>%
  arrange(mean_coef)

# export table S3
table_S3 <- parameters_m_s_units %>% 
  dplyr::select(pred_names_char, value, unit_long) %>% 
  arrange(unit_long) %>% 
  mutate(value = round(value,3)) %>% 
  pivot_wider(names_from = unit_long, values_from = value) 

write.csv(table_S3,  file.path(out_path, paste0("table_coefficient_values_Aceh_N_Sumatra_",forest_layer, ".csv")), row.names = F)

variables_1 <- parameters_m_s_units_filter$pred_names_char
parameters_m_s_units_filter <- parameters_m_s_units_filter %>%
  left_join(parameters_m_s_units_filter_overview, by = "pred_names_char") %>%
  filter(!(abs(mean_coef) < 0.05 & spread < 0.1) ) %>%
  filter( pred_names_char != "intercept") 

variables_2 <- parameters_m_s_units_filter$pred_names_char

#variables_1[
 unique(variables_1)[!(unique(variables_1) %in% unique(variables_2))]
 
p1 <- ggplot(data = parameters_m_s_units ) + # 
  geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed")+
  geom_boxplot(aes(y = value, x = pred_names_char ), color = "black", fill = "gray89", width = 0.8, outlier.shape=NA)+
  geom_jitter(aes( y = value, x = pred_names_char , color = unit_long), shape = 19, size=3, width = 0.3, alpha = 0.9)+
  # scale_shape_manual(values = shape_custom)+
  scale_color_manual(values = cbPalette_custom)+
  coord_flip()+
  # labs(fill = "test") +
  theme_bw() +
  ylab("predictor coefficient") +
  xlab("predictors") +
  ggtitle(paste0("all predictors for ", forest_layer, " layer"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size = 12, face = "bold", hjust = 1),
    axis.title=element_text(size = 12, face = "bold"),
    #  axis.text.x=element_text(angle=45),
    legend.text=element_text(size=12),
    legend.key.height=unit(0.8,"cm"),
    legend.title = element_blank(),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 12, face="bold")) +
  theme(plot.margin = unit(c(0.2,-0.2,0.1,0.1),"cm")) 


pdf(file.path(out_path, paste0("Parameter_effect_size_Aceh_N_Sumatra_all_",forest_layer,".pdf")), width =0, height = 0, paper = "a4r")
print(p1)
dev.off()

plot2 <-  ggplot(data = parameters_m_s_units_filter ) + # 
  #ggplot(data = parameters_m_s_units ) + # 
  geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed")+
  geom_boxplot(aes(y = value, x = pred_names_char ), color = "black", fill = "gray89", width = 0.8, outlier.shape=NA)+
  geom_jitter(aes( y = value, x = pred_names_char , color = unit_long), shape = 19, size=3, width = 0.3, alpha = 0.9)+
  #scale_shape_manual(values = shape_custom)+
  scale_color_manual(values = cbPalette_custom)+
  coord_flip()+
  # labs(fill = "test") +
  theme_bw() +
  ylab("predictor coefficient") +
  xlab("predictors") +
  ggtitle(paste0("predictors with effect for ", forest_layer, " layer"))+
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size = 12, face = "bold", hjust = 1),
    axis.title=element_text(size = 12, face = "bold"),
    #  axis.text.x=element_text(angle=45),
    legend.text=element_text(size=12),
    legend.key.height=unit(0.8,"cm"),
    legend.title = element_blank(),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(size = 12, face="bold")) +
  theme(plot.margin = unit(c(0.2,-0.2,0.1,0.1),"cm")) 



pdf(file.path(out_path, paste0("Parameter_effect_size_Aceh_N_Sumatra_important_",forest_layer,".pdf")), width =0, height = 0, paper = "a4r")
print(plot2)
dev.off()

# 
} # close forest layer loop
#