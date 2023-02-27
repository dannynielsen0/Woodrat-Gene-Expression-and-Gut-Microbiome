###Analysis for 2019 laboratory feeding experiments

#clear working environment and set directory
rm(list=ls())
setwd("/Volumes/GoogleDrive/My Drive/dissertation/ch. 3/Feeding_trials/Trial_data")

#load some libraries

library(ggplot2)
library(data.table)
library(lme4)
library(sjPlot)
library(car)
library(emmeans)
library(ggpubr)
library(estimatr)
library(recipes)
library(car)


#read in data
diet_data <- read.csv("All_trial_data_Wwheel.csv", header=TRUE)
diet_data$Day_of_trial <- as.factor(diet_data$Day_of_trial)
diet_data$Diet_threshold <- as.factor(diet_data$Diet_threshold)

#filter data to only beginning and end points

high_low_data <- diet_data %>%
  filter(Diet_threshold == "Chow" | Diet_threshold == "Max")

#create a variable for species, diet and percentage
high_low_data$sp_diet_thresh <- paste(high_low_data$Species, high_low_data$Diet_type, 
                                      high_low_data$Diet_threshold, sep="_")




#run model with lm_robust from estimatr model

water_model <- lm_robust(Water_consumed_day~sp_diet_thresh,
                         clusters=Woodrat_id,
                         data=high_low_data)

summary(water_model)

#minutes active per day
activity_model <- lm_robust(Minutes_active~sp_diet_thresh + WR_bm + Sex,
                         clusters=Woodrat_id,
                         data=high_low_data)
summary(activity_model)

#max speed
max_speed_model <- lm_robust(Max_speed_ms~sp_diet_thresh + WR_bm + Sex,
                            clusters=Woodrat_id,
                            data=high_low_data)
summary(max_speed_model)


#food/day
food_model <- lm_robust(Food_eaten_day~sp_diet_thresh + WR_bm + Sex,
                             clusters=Woodrat_id,
                             data=high_low_data)
summary(food_model)






