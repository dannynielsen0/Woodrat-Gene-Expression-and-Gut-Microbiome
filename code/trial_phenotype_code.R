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

#read in data
diet_data <- read.csv("All_trial_data_Wwheel.csv", header=TRUE)
diet_data$Day_of_trial <- as.factor(diet_data$Day_of_trial)

#melt the max dose data
max_dose_long <- melt(diet_data, id=c("Diet_threshold","Woodrat_id", "Diet_type", "Diet_perc", "Species", "Pop", "Trial", "Sex", "WR_bm"), measure="mass.adjusted_dose") #make long format
colnames(max_dose_long)[c(1,11)] <- c("Dose", "MA_dose") #rename columns
max_dose_long$Dose[max_dose_long$Dose=="Chow"] <- "0" 

#subset to just 0% and max dose
max_dose_long_sub <- subset(max_dose_long, max_dose_long$Dose == "Max" | max_dose_long$Dose=="0")
#then average over treatment groups
max_dose_long_sub_ave <- aggregate(Diet_perc~ Woodrat_id + Diet_type + Species, FUN=mean, data=max_dose_long_sub) #get average water intake for each woodrat at each dose
max_dose_long_sub_ave_bm <- aggregate(WR_bm~ Woodrat_id, FUN=mean, data=max_dose_long_sub) #get average body mass for each woodrat
max_dose_long_sub_ave$WR_bm <- max_dose_long_sub_ave_bm$WR_bm[match(max_dose_long_sub_ave$Woodrat_id, max_dose_long_sub_ave_bm$Woodrat_id)]

#make species_diet variable
max_dose_long_sub_ave$Sp_diet <- paste(max_dose_long_sub_ave$Species,max_dose_long_sub_ave$Diet_type, sep="_")
max_dose_long_sub_ave$Sp_diet <- factor(max_dose_long_sub_ave$Sp_diet, levels=c("N. lepida_PRFA", "N. bryanti_FRCA", "N. lepida_FRCA", "N. bryanti_PRFA"))

max_dose_long_sub_ave$MA_dose <-max_dose_long_sub_ave$Diet_perc/max_dose_long_sub_ave$WR_bm

max_dose_mod <- aov(MA_dose ~ Sp_diet , data=max_dose_long_sub_ave)
max_dose_mod_mm <- lmer(MA_dose ~ Species * Diet_type + Pop + (1|Woodrat_id), data=max_dose_long_sub)
drop1(max_dose_mod_mm, test="Chisq")
summary(max_dose_mod_mm)
coef(max_dose_mod_mm)

sjPlot::plot_model(max_dose_mod_mm, show.values = TRUE, show.p = TRUE)

summary(max_dose_mod)
coef(max_dose_mod)
TukeyHSD(max_dose_mod, which="Sp_diet")

max_crap_is <- emmeans(max_dose_mod, ~ Sp_diet)
max_tukey <- multcomp::cld(max_crap_is, alpha=0.05, Letters=letters)


#mass adjusted max dose plot

max_plot <- ggplot(data=max_dose_long_sub_ave, aes(x=Sp_diet , y=MA_dose)) + 
  geom_boxplot() + theme_bw() + geom_jitter() +
  theme(axis.text.x = element_text(size = 20, face = "italic")) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + 
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  #facet_wrap(~Species) + 
  stat_summary(fun="mean", geom="point", col="red") +
  theme(strip.text.x = element_text(size = 24)) + theme(strip.text = element_text(face = "italic"))+
  xlab("Species/Diet") + ylab("Mass Adjusted Max Dose")

max_plot

ggsave(plot=max_plot, "Lab-diet-trial-16S-analysis/figures/Max_dose.jpg", width = 11, height = 6, device='jpg', dpi=500)


#reformat water intake data for plotting

water_long <- melt(diet_data, id=c("Diet_threshold","Woodrat_id", "Diet_type", "Species", "Trial", "Sex", "WR_bm"), measure="Water_consumed_day") #make long format
names(water_long)[9] <- "Water_intake"
names(water_long)[1] <- "Dose"
water_long$MA_water <- water_long$Water_intake/water_long$WR_bm
water_long <- aggregate(MA_water ~ Woodrat_id + Diet_type + Species + Dose, FUN=mean, data=water_long) #get average water intake for each woodrat at each dose
water_long <- subset(water_long, water_long$Dose!= "Dose")

#rename chow to number 0
water_long$Dose[water_long$Dose=="Chow"] <- "0" 

#reorder variabels for plotting
water_long$Species <- factor(water_long$Species, levels=c("N. lepida", "N. bryanti"))
water_long$Diet_type <- factor(water_long$Diet_type, levels=c("PRFA", "FRCA"))
water_long$Sp_diet <- paste(water_long$Species, water_long$Diet_type, sep="_")
water_long$Sp_diet <- factor(water_long$Sp_diet, levels=c("N. lepida_PRFA", "N. bryanti_FRCA", "N. lepida_FRCA", "N. bryanti_PRFA"))
water_long$Sp_diet_dose <- paste(water_long$Species, water_long$Diet_type, water_long$Dose, sep="_")

#conduct statistical analysis for each dietXspecies treatment
shapiro.test(water_long$MA_water) #W = 0.57666, p-value = 1.131e-1;  not normal, justify kruskal wallis test for each speciesxdiet treatment group
bry_prfa <- subset(water_long, water_long$Sp_diet =="N. bryanti_PRFA")
bry_frca <- subset(water_long, water_long$Sp_diet =="N. bryanti_FRCA")
lep_prfa <- subset(water_long, water_long$Sp_diet =="N. lepida_PRFA")
lep_frca <- subset(water_long, water_long$Sp_diet =="N. lepida_FRCA")

#perform kruskal test for each
kruskal.test(MA_water~Dose, data=bry_prfa) #Kruskal-Wallis chi-squared = 0.39706, df = 1, p-value = 0.5286; no difference in bry
kruskal.test(MA_water~Dose, data=bry_frca) #Kruskal-Wallis chi-squared = 0.10204, df = 1, p-value = 0.7494; no difference in bry
kruskal.test(MA_water~Dose, data=lep_frca) #Kruskal-Wallis chi-squared = 3.9224, df = 1, p-value = 0.04765; df = 1, p-value = 0.7494; no difference in bry
kruskal.test(MA_water~Dose, data=lep_prfa) #Kruskal-Wallis chi-squared = 0.2, df = 1, p-value = 0.6547; df = 1, p-value = 0.7494; no difference in bry


water_plot <- ggplot(data=water_long, aes(x=Dose, y=MA_water)) + geom_boxplot() + 
  theme_bw() + facet_grid(~ Sp_diet) +
  stat_summary(fun="mean", geom="point", col="red") +
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  theme(strip.text = element_text(face = "italic")) + ylab("Mass adjusted water/day") + xlab("% Dose")

#save plot
ggsave(plot=water_plot, "Lab-diet-trial-16S-analysis/figures/MA_adjusted_water.jpg", width = 12, height = 6, device='jpg', dpi=500)




#reformat food intake data for plotting

food_long <- melt(diet_data, id=c("Diet_threshold","Woodrat_id", "Diet_type", "Species", "Trial", "Sex", "WR_bm"), measure="Food_eaten_day") #make long format
names(food_long)[9] <- "Food_intake"
names(food_long)[1] <- "Dose"
food_long$MA_food <- food_long$Food_intake/food_long$WR_bm
food_long <- aggregate(MA_food ~ Woodrat_id + Diet_type + Species + Dose, FUN=mean, data=food_long) #get average water intake for each woodrat at each dose
food_long <- subset(food_long, food_long$Dose!= "Dose")

#rename chow to number 0
food_long$Dose[food_long$Dose=="Chow"] <- "0" 

#reorder variabels for plotting
food_long$Species <- factor(food_long$Species, levels=c("N. lepida", "N. bryanti"))
food_long$Diet_type <- factor(food_long$Diet_type, levels=c("PRFA", "FRCA"))
food_long$Sp_diet <- paste(food_long$Species, food_long$Diet_type, sep="_")
food_long$Sp_diet <- factor(food_long$Sp_diet, levels=c("N. lepida_PRFA", "N. bryanti_FRCA", "N. lepida_FRCA", "N. bryanti_PRFA"))
food_long$Sp_diet_dose <- paste(food_long$Species, food_long$Diet_type, food_long$Dose, sep="_")

#conduct statistical analysis for each dietXspecies treatment
shapiro.test(food_long$MA_food) #W = 0.97403, p-value = 0.2473; normal, justify t.test test for each speciesxdiet treatment group
bry_prfa <- subset(food_long, food_long$Sp_diet =="N. bryanti_PRFA")
bry_frca <- subset(food_long, food_long$Sp_diet =="N. bryanti_FRCA")
lep_prfa <- subset(food_long, food_long$Sp_diet =="N. lepida_PRFA")
lep_frca <- subset(food_long, food_long$Sp_diet =="N. lepida_FRCA")

#perform kruskal test for each
t.test(MA_food~Dose, data=bry_prfa) #t.test; t = -1.7262, df = 13.316, p-value = 0.1074; no difference in bry
t.test(MA_food~Dose, data=bry_frca) #t.test; t = -0.40605, df = 11.204, p-value = 0.6924; no difference in bry
t.test(MA_food~Dose, data=lep_frca) #t = 1.2604, df = 11.662, p-value = 0.2322; df = 1, p-value = 0.7494; no difference in bry
t.test(MA_food~Dose, data=lep_prfa) #t = -1.7199, df = 11.764, p-value = 0.1116; no difference in bry


food_plot <- ggplot(data=food_long, aes(x=Dose, y=MA_food)) + geom_boxplot() + 
  theme_bw() + facet_grid(~ Sp_diet) +
  stat_summary(fun="mean", geom="point", col="red") +
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  theme(strip.text = element_text(face = "italic")) + ylab("Mass adjusted food/day") + xlab("% Dose")

#save plot
ggsave(plot=food_plot, "Lab-diet-trial-16S-analysis/figures/MA_adjusted_food.jpg", width = 12, height = 6, device='jpg', dpi=500)





#reformat wheel data intake data for plotting

minutes_long <- melt(diet_data, id=c("Diet_threshold","Woodrat_id", "Diet_type", "Species", "Trial", "Sex", "WR_bm"), measure="Minutes_active") #make long format
names(minutes_long)[9] <- "minutes_active"
names(minutes_long)[1] <- "Dose"
minutes_long$MA_minutes <- minutes_long$minutes_active/minutes_long$WR_bm
minutes_long <- aggregate(MA_minutes ~ Woodrat_id + Diet_type + Species + Dose, FUN=mean, data=minutes_long) #get average water intake for each woodrat at each dose
minutes_long <- subset(minutes_long, minutes_long$Dose!= "Dose")

#rename chow to number 0
minutes_long$Dose[minutes_long$Dose=="Chow"] <- "0" 

#reorder variabels for plotting
minutes_long$Species <- factor(minutes_long$Species, levels=c("N. lepida", "N. bryanti"))
minutes_long$Diet_type <- factor(minutes_long$Diet_type, levels=c("PRFA", "FRCA"))
minutes_long$Sp_diet <- paste(minutes_long$Species, minutes_long$Diet_type, sep="_")
minutes_long$Sp_diet <- factor(minutes_long$Sp_diet, levels=c("N. lepida_PRFA", "N. bryanti_FRCA", "N. lepida_FRCA", "N. bryanti_PRFA"))
minutes_long$Sp_diet_dose <- paste(minutes_long$Species, minutes_long$Diet_type, minutes_long$Dose, sep="_")

#conduct statistical analysis for each dietXspecies treatment
shapiro.test(minutes_long$MA_minutes) #W = 0.93062, p-value = 0.002576; not normal, justify kruskal.wallis test for each speciesxdiet treatment group
bry_prfa <- subset(minutes_long, minutes_long$Sp_diet =="N. bryanti_PRFA")
bry_frca <- subset(minutes_long, minutes_long$Sp_diet =="N. bryanti_FRCA")
lep_prfa <- subset(minutes_long, minutes_long$Sp_diet =="N. lepida_PRFA")
lep_frca <- subset(minutes_long, minutes_long$Sp_diet =="N. lepida_FRCA")

#perform kruskal test for each
kruskal.test(MA_minutes~Dose, data=bry_prfa) #Kruskal-Wallis chi-squared = 0, df = 1, p-value = 1; no difference in bry
kruskal.test(MA_minutes~Dose, data=bry_frca) #Kruskal-Wallis chi-squared = 0.10204, df = 1, p-value = 0.7494; no difference in bry
kruskal.test(MA_minutes~Dose, data=lep_frca) #Kruskal-Wallis chi-squared = 2.9755, df = 1, p-value = 0.08453; no difference in bry
kruskal.test(MA_minutes~Dose, data=lep_prfa) #Kruskal-Wallis chi-squared = 0.6898, df = 1, p-value = 0.4062; no difference in bry


minutes_plot <- ggplot(data=minutes_long, aes(x=Dose, y=MA_minutes)) + geom_boxplot() + 
  theme_bw() + facet_grid(~ Sp_diet) +
  stat_summary(fun="mean", geom="point", col="red") +
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  theme(strip.text = element_text(face = "italic")) + ylab("Mass adjusted minutes/day") + xlab("% Dose")

#save plot
ggsave(plot=minutes_plot, "Lab-diet-trial-16S-analysis/figures/MA_adjusted_minutes.jpg", width = 12, height = 6, device='jpg', dpi=500)


#reformat wheel data intake data for plotting

rotations_long <- melt(diet_data, id=c("Diet_threshold","Woodrat_id", "Diet_type", "Species", "Trial", "Sex", "WR_bm"), measure="rotations") #make long format
names(rotations_long)[9] <- "rotations"
names(rotations_long)[1] <- "Dose"
rotations_long$MA_rotations <- rotations_long$rotations/rotations_long$WR_bm
rotations_long <- aggregate(MA_rotations ~ Woodrat_id + Diet_type + Species + Dose, FUN=mean, data=rotations_long) #get average water intake for each woodrat at each dose
rotations_long <- subset(rotations_long, rotations_long$Dose!= "Dose")

#rename chow to number 0
rotations_long$Dose[rotations_long$Dose=="Chow"] <- "0" 

#reorder variabels for plotting
rotations_long$Species <- factor(rotations_long$Species, levels=c("N. lepida", "N. bryanti"))
rotations_long$Diet_type <- factor(rotations_long$Diet_type, levels=c("PRFA", "FRCA"))
rotations_long$Sp_diet <- paste(rotations_long$Species, rotations_long$Diet_type, sep="_")
rotations_long$Sp_diet <- factor(rotations_long$Sp_diet, levels=c("N. lepida_PRFA", "N. bryanti_FRCA", "N. lepida_FRCA", "N. bryanti_PRFA"))
rotations_long$Sp_diet_dose <- paste(rotations_long$Species, rotations_long$Diet_type, rotations_long$Dose, sep="_")

#conduct statistical analysis for each dietXspecies treatment
shapiro.test(rotations_long$MA_rotations) #W = 0.91212, p-value = 0.0004735; not normal, justify kruskal.wallis test for each speciesxdiet treatment group
bry_prfa <- subset(rotations_long, rotations_long$Sp_diet =="N. bryanti_PRFA")
bry_frca <- subset(rotations_long, rotations_long$Sp_diet =="N. bryanti_FRCA")
lep_prfa <- subset(rotations_long, rotations_long$Sp_diet =="N. lepida_PRFA")
lep_frca <- subset(rotations_long, rotations_long$Sp_diet =="N. lepida_FRCA")

#perform kruskal test for each
kruskal.test(MA_rotations~Dose, data=bry_prfa) #Kruskal-Wallis chi-squared = 0.17647, df = 1, p-value = 0.6744; no difference in bry
kruskal.test(MA_rotations~Dose, data=bry_frca) #Kruskal-Wallis chi-squared = 0.036735, df = 1, p-value = 0.848; no difference in bry
kruskal.test(MA_rotations~Dose, data=lep_frca) #Kruskal-Wallis chi-squared = 1.4735, df = 1, p-value = 0.2248; no difference in bry
kruskal.test(MA_rotations~Dose, data=lep_prfa) #Kruskal-Wallis chi-squared = 0.2, df = 1, p-value = 0.6547; no difference in bry


rotations_plot <- ggplot(data=rotations_long, aes(x=Dose, y=MA_rotations)) + geom_boxplot() + 
  theme_bw() + facet_grid(~ Sp_diet) +
  stat_summary(fun="mean", geom="point", col="red") +
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  theme(strip.text = element_text(face = "italic")) + ylab("Mass adjusted rotations/day") + xlab("% Dose")

#save plot
ggsave(plot=rotations_plot, "Lab-diet-trial-16S-analysis/figures/MA_adjusted_rotations.jpg", width = 12, height = 6, device='jpg', dpi=500)


#reformat wheel data intake data for plotting max speed

speed_long <- melt(diet_data, id=c("Diet_threshold","Woodrat_id", "Diet_type", "Species", "Trial", "Sex", "WR_bm"), measure="Max_speed_ms") #make long format
names(speed_long)[9] <- "speed"
names(speed_long)[1] <- "Dose"
speed_long$MA_speed <- speed_long$speed/speed_long$WR_bm
speed_long <- aggregate(MA_speed ~ Woodrat_id + Diet_type + Species + Dose, FUN=mean, data=speed_long) #get average water intake for each woodrat at each dose
speed_long <- subset(speed_long, speed_long$Dose!= "Dose")

#rename chow to number 0
speed_long$Dose[speed_long$Dose=="Chow"] <- "0" 

#reorder variabels for plotting
speed_long$Species <- factor(speed_long$Species, levels=c("N. lepida", "N. bryanti"))
speed_long$Diet_type <- factor(speed_long$Diet_type, levels=c("PRFA", "FRCA"))
speed_long$Sp_diet <- paste(speed_long$Species, speed_long$Diet_type, sep="_")
speed_long$Sp_diet <- factor(speed_long$Sp_diet, levels=c("N. lepida_PRFA", "N. bryanti_FRCA", "N. lepida_FRCA", "N. bryanti_PRFA"))
speed_long$Sp_diet_dose <- paste(speed_long$Species, speed_long$Diet_type, speed_long$Dose, sep="_")

#conduct statistical analysis for each dietXspecies treatment
shapiro.test(speed_long$MA_speed) #W = 0.91227, p-value = 0.0004797; not normal, justify kruskal.wallis test for each speciesxdiet treatment group
bry_prfa <- subset(speed_long, speed_long$Sp_diet =="N. bryanti_PRFA")
bry_frca <- subset(speed_long, speed_long$Sp_diet =="N. bryanti_FRCA")
lep_prfa <- subset(speed_long, speed_long$Sp_diet =="N. lepida_PRFA")
lep_frca <- subset(speed_long, speed_long$Sp_diet =="N. lepida_FRCA")

#perform kruskal test for each
kruskal.test(MA_speed~Dose, data=bry_prfa) #Kruskal-Wallis chi-squared = 0.099265, df = 1, p-value = 0.7527, df = 1, p-value = 0.6744; no difference in bry
kruskal.test(MA_speed~Dose, data=bry_frca) #Kruskal-Wallis chi-squared = 0.6898, df = 1, p-value = 0.4062, df = 1, p-value = 0.848; no difference in bry
kruskal.test(MA_speed~Dose, data=lep_frca) #Kruskal-Wallis chi-squared = 0.49388, df = 1, p-value = 0.4822, df = 1, p-value = 0.2248; no difference in bry
kruskal.test(MA_speed~Dose, data=lep_prfa) #Kruskal-Wallis chi-squared = 0.91837, df = 1, p-value = 0.3379; no difference in bry


speed_plot <- ggplot(data=speed_long, aes(x=Dose, y=MA_speed)) + geom_boxplot() + 
  theme_bw() + facet_grid(~ Sp_diet) +
  stat_summary(fun="mean", geom="point", col="red") +
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  theme(strip.text = element_text(face = "italic")) + ylab("Mass adjusted max m/s") + xlab("% Dose")

#save plot
ggsave(plot=speed_plot, "Lab-diet-trial-16S-analysis/figures/MA_adjusted_speed.jpg", width = 12, height = 6, device='jpg', dpi=500)





