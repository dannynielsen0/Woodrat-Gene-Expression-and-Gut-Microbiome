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
  geom_text(data=max_tukey,aes(x=Sp_diet, y=0.75,label=max_tukey$.group),vjust=-0.15) + xlab("Species/Diet") + ylab("Mass Adjusted Max Dose")

max_plot

ggsave(plot=max_plot, "Max_dose.jpg", width = 11, height = 6, device='jpg', dpi=500)


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


ggsave(plot=water_plot, "MA_adjusted_water.jpg", width = 12, height = 6, device='jpg', dpi=500)




#reformat food intake data for plotting

food_long <- melt(diet_data, id=c("Diet.","Woodrat_id", "Diet_type", "Species", "Trial", "Sex"), measure="Food_eaten_day") #make long format
food_long <- aggregate(value ~ Woodrat_id + Diet., FUN=mean, data=food_long) #get average water intake for each woodrat at each dose

#match back the species, sex, and trial data to this new averaged dataset
food_long$Species <- diet_data$Species[match(food_long$Woodrat_id, diet_data$Woodrat_id)]
food_long$Sex <- diet_data$Sex[match(food_long$Woodrat_id, diet_data$Woodrat_id)]
food_long$Trial <- diet_data$Trial[match(food_long$Woodrat_id, diet_data$Woodrat_id)]
food_long$Diet <- diet_data$Diet_type[match(food_long$Woodrat_id, diet_data$Woodrat_id)]

colnames(food_long)[2:3] <- c("Dose", "Food_intake") #rename columns


###Plotting food intake

food_plot <- ggplot(data=food_long, aes(x=Dose, y=Food_intake)) + geom_point() +
  facet_grid(~Species + Diet) + stat_summary(fun="median", color="red")

#run a model on food

food_model <- aov(Food_intake ~ Dose + Species + Diet + Species*Diet + Error(Sex,Trial), data=food_long)
food_model <- aov(Food_intake ~ Dose + Species + Diet + Species*Diet, data=food_long)

summary(food_model, type="III")

#mixed effects model food

mod_mm_food <- lmer(Food_intake ~ Dose + Species + Diet + Species*Diet + (1|Trial) + (1|Sex), data=food_long)
summary(mod_mm_food)
coef(mod_mm_food)

#plots and tables
#water
water_mm_mod_plot <- sjPlot::plot_model(mod_mm, show.values = TRUE, show.p = TRUE)
water_mm_mod_table <- sjPlot::tab_model(mod_mm)
save.image(water_mm_mod_table, file="water_mm_table.jpg")


#food
food_mm_mod_plot <- sjPlot::plot_model(mod_mm_food, show.values = TRUE, show.p = TRUE)
food_mm_mod_table <- sjPlot::tab_model(mod_mm_food)

diet_trial_effectsizes <- ggarrange(water_mm_mod_plot, food_mm_mod_plot, nrow=1)
ggsave(plot=last_plot(), filename = "effect_sizes.jpg", dpi=500)



###Individual tracking plots

water_plot <- ggplot(data=diet_data, aes(x=Diet., y=Water_consumed_day, color = Diet_type)) + 
  geom_line(aes(linetype=Woodrat_id)) +
  #geom_point()  +
  facet_wrap(~Species)

water_plot + ggtitle("Water intake") + theme_bw() + theme(panel.grid.minor = element_blank()) + scale_x_continuous(breaks = 1:15)


#reformat wheel data intake data for plotting

activity_long <- melt(diet_data, id=c("Diet.", "Diet_threshold","Woodrat_id", "Diet_type", "Species", "Pop", "Trial", "Sex", "WR_bm"), measure="Minutes_active") #make long format
activity_long$MA_activity <- activity_long$value/activity_long$WR_bm
activity_long <- aggregate(MA_activity ~ Woodrat_id + Diet_type + Species + Pop + Diet. + Diet_threshold, FUN=mean, data=activity_long) #get average water intake for each woodrat at each dose
activity_long <- subset(activity_long, activity_long$Diet_threshold != "Dose")


#match back the species, sex, and trial data to this new averaged dataset
activity_long$Species <- diet_data$Species[match(activity_long$Woodrat_id, diet_data$Woodrat_id)]
activity_long$Pop <- diet_data$Species[match(activity_long$Woodrat_id, diet_data$Woodrat_id)]
activity_long$Sex <- diet_data$Sex[match(activity_long$Woodrat_id, diet_data$Woodrat_id)]
activity_long$Trial <- diet_data$Trial[match(activity_long$Woodrat_id, diet_data$Woodrat_id)]
activity_long$Diet <- diet_data$Diet_type[match(activity_long$Woodrat_id, diet_data$Woodrat_id)]
activity_long$Woodrat_id <- diet_data$Woodrat_id[match(activity_long$Woodrat_id, diet_data$Woodrat_id)]
activity_long$Woodrat_id <- diet_data$Woodrat_id[match(activity_long$Woodrat_id, diet_data$Woodrat_id)]
activity_long$Species_diet_dose <- paste(activity_long$Species, activity_long$Diet_type, activity_long$Diet., sep="_")
activity_long$Species_diet <- paste(activity_long$Species, activity_long$Diet_type, sep="_")
activity_long$Species_diet_threshold <- paste(activity_long$Species, activity_long$Diet_type, activity_long$Diet_threshold, sep="_")



colnames(activity_long)[c(4,6)] <- c("Dose", "activity") #rename columns
activity_long$Dose[activity_long$Diet_threshold=="Chow"] <- "0" 

#reorder variabels for plotting
activity_long$Species <- factor(activity_long$Species, levels=c("N. lepida", "N. bryanti"))
activity_long$Diet_type <- factor(activity_long$Diet_type, levels=c("PRFA", "FRCA"))

activity_mod_mm <- glmmTMB(MA_activity ~ Species + Diet_type +  Diet_threshold + Species*Diet_type + (1|Woodrat_id), data=activity_long)
activity_mod_mm <- glmmTMB(activity ~ Species_diet_threshold + (1|Woodrat_id), data=activity_long)

summary(activity_mod_mm)
coef(activity_mod_mm)

qqnorm(log10(activity_long$activity))
activity_long$activity <- sqrt(activity_long$activity)
ggdensity(activity_long, x = "activity", fill = "lightgray") +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

library(emmeans)
activity_crap_is <- emmeans(activity_mod_mm, ~ Species_diet_threshold)
activity_tukey <- multcomp::cld(activity_crap_is, alpha=0.05, Letters=letters)
TukeyHSD(activity_mod_mm)

###Plotting activity
activity_long$Sp_diet <- paste(activity_long_sub$Species, activity_long$Diet_type, sep="_")
activity_long$Sp_diet <- factor(activity_long$Sp_diet, levels=c("N. lepida_PRFA", "N. bryanti_FRCA", "N. lepida_FRCA", "N. bryanti_PRFA"))
activity_long$Sp_diet_dose <- paste(activity_long$Species, activity_long$Diet_type, activity_long$Dose, sep="_")


activity_plot <- ggplot(data=activity_long, aes(x=Species_diet_threshold, y=activity)) + geom_boxplot() + 
  theme_bw() + #facet_grid(~ Species) +
  stat_summary(fun="mean", geom="point", col="red") +
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  theme(strip.text = element_text(face = "italic")) + ylab("rotations") + #xlab("% Dose") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  geom_text(data=activity_tukey,aes(x=Species_diet_threshold, y=2.2,label=activity_tukey$.group),vjust=-0.15)

ggsave("rotations.jpg", width = 12, height = 6, device='jpg', dpi=300)

#mixed model for activity with individual (and trial?) as random effect



stat.test.activity <- ggpubr::compare_means(
  MA_water ~ Dose, group="Sp_diet",
  data = water_long_sub, method = "kruskal",
  p.adjust.method = "none"
)
stat.test.water

my_comps <- stat.test.water %>% mutate(y.position=c(0.3, 0.3,0.55,0.35))



water_plot <- ggplot(data=water_long_sub, aes(x=Dose, y=MA_water)) + geom_boxplot() + 
  theme_bw() + facet_grid(~ Sp_diet) +
  theme(panel.grid.major = element_blank()) + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20,face="bold")) + theme(strip.text.x = element_text(size = 24)) +
  theme(strip.text = element_text(face = "italic")) + ylab("Mass adjusted water/day") + xlab("% Dose") +
  stat_pvalue_manual(my_comps, label="p.signif") + 
  stat_summary(fun="mean", geom="point", col="lightgrey")


water_plot
ggsave(plot=water_plot, "MA_adjusted_water.jpg", width = 12, height = 6, device='jpg', dpi=500)





#run a model on water
water_long$Dose <- as.factor(water_long$Dose)
water_long$Species <- as.factor(water_long$Species)
water_long$Sex <- as.factor(water_long$Sex)
water_long$Trial <- as.factor(water_long$Trial)
water_long$Diet <- as.factor(water_long$Diet)
water_long$Woodrat_id <- as.factor(water_long$Woodrat_id)
water_long$Diet_type[water_long$Dose=="0"] <- "Chow"
water_long_sub <- subset(water_long, water_long$Dose != "Dose")

water_model <- aov(MA_water ~ Species_diet_dose, data=water_long)

max_crap_is <- emmeans(water_model, ~ Species_diet_dose)
max_tukey <- multcomp::cld(max_crap_is, alpha=0.05, Letters=letters)

water_model <- aov(MA_water ~ Dose + Species * Diet_type + WR_bm , data=water_long)
Anova(water_model, type="III")

water_model <- lm(Water_intake ~ Dose * Species * Diet_type + Woodrat_id + WR_bm, data=water_long)
wilcox.test(MA_water ~ Dose, paired=TRUE, data=water_long)

summary(water_model)
coef(water_model)
TukeyHSD(water_model, which="Sp_diet_dose")


water_lin_model <- lm(Water_intake ~ Dose + Species + Diet_type + Species*Diet_type + Woodrat_id, data=water_long)
summary(water_lin_model)

water_anova <- nlme::lme(MA_water ~ Dose + Species + Diet_type + Species*Diet_type, random = ~ 1|Woodrat_id, data = water_long, method = "REML")
#repeated measures anova
summary(water_anova)

TukeyHSD(water_anova)


#mixed effects model to include trial as a random effect

qqnorm(water_long$MA_water)
shapiro.test(water_long$ma_water_sqrt)
water_long_sub$ma_water_sqrt <- sqrt(water_long_sub$MA_water)

mod_mm <- lmer(MA_water ~ Dose + Species * Diet + (1|Trial) + (1|Woodrat_id), data=water_long_sub)
mod_mm <- lmer(MA_water ~ Sp_diet_dose + (1|Trial) + (1|Woodrat_id), data=water_long_sub)

drop1(mod_mm,test = "Chisq") #likelihood ratio test

mod_mm <- lmer(Water_intake ~ Dose * Species * Diet + (1|Trial) + (1|Sex), data=water_long)

summary(mod_mm)
coef(mod_mm)




