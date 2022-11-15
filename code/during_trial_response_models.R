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

water_model <- lm_robust(Water_consumed_day~sp_diet_thresh + WR_bm + Sex,
                         clusters=Woodrat_id,
                         data=high_low_data)

summary(water_model)
car::Anova(water_model)

water_model <- aov(Water_consumed_day~Species * Diet_type * Diet_threshold + WR_bm +
                     Species + Sex,
                         data=high_low_data)

summary(water_model)
car::Anova(water_model)
TukeyHSD(water_model, which = c("Diet_type","Diet_threshold"))


mixed <- lmer(
  Water_consumed_day ~ Species * Diet_type * Diet_threshold + WR_bm + Sex + (1|Woodrat_id),
  data = high_low_data
)

summary(mixed, correlation = FALSE)
anova(mixed)

library(emmeans)

marginals <- emmeans(mixed, ~Species * Diet_type * Diet_threshold)
plot(marginals)

emmip(water_model, Species ~ Diet_threshold | Diet_type) +
  scale_color_brewer(type = "qual") +
  theme_classic(16)

