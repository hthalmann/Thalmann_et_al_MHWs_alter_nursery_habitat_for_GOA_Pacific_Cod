## Linear Mixed Effects Model Selection and Analysis to evaluate the factors influencing juvenile Pacific Cod relative growth rates in July and August. Calculation of predicted size for August fish based on growth rates from July. 

## Date Created: 1/11/24
## R version: 4.2.2


## Load Libraries
library(rstudioapi) #load R working directory
library(tidyverse)  #data wrangling and visualization
library(patchwork) #Data Visualization
library(lme4) # Run mixed effects linear models
library(nlme) #Run mixed effects models with autocorrelation 
library(effects) #Marginal means
library(ggeffects)#Visualize marginal means https://strengejacke.github.io/ggeffects/
library(car) #ANOVA function for models
library(MuMIn) #Runs the R.squaredGLMM function
library(viridis) #Color Palette

## Paired with:
#Relative Growth with the Biological Intercept Back calculation
#Size and Condition Data
#Prey-Specific Index of Relative Importance (PSIRI) data matrices for July and August
#Stomach Fullness
#Trident Bay Daily Temperature
#Final 21 Days with Date to match date to daily increments and temperature
#Sampling Data to match Seine ID, Site ID, and Fish ID
#Seine Data for cod abundance

##Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#Sets directory to folder that this R script is saved to


#### Growth Data Wrangling ####

Growth <- read.csv("D12.BackCalculatedGrowth_ROutput.csv")

#Growth Summary
GrowthSummary <- Growth %>%
  group_by(Month,Heatwave) %>%
  reframe(mm_day_mean = mean(mm_day_BIC), mm_day_sd = sd(mm_day_BIC), mm_mm_day_mean = mean(mm_mm_day_BIC), mm_mm_day_sd = sd(mm_mm_day_BIC), mm_mm_day_se = (se = sd(mm_mm_day_BIC, na.rm = T)/sqrt(length(mm_mm_day_BIC))))


## Diet Files 
JulyNMS <- read.csv("D13.JulyNMS_ROutput.csv") %>%
  select(FISHID, NMDS1, NMDS2, NMDS3)
AugustNMS <- read.csv("D14.AugNMS_ROutput.csv") %>%
  select(FISHID, NMDS1, NMDS2, NMDS3)

AllDiet <- rbind(JulyNMS, AugustNMS)

Fullness <- read.csv("D15.StomachFullness_ROutput.csv") %>%
  select(FISHID, StomachFullness, sqrtFullness)

## Temperature Files
TridentTemps <- read.csv("D4.TridentBay_DailyTemperature.csv") %>%
  mutate(Date = ymd(Date))%>%
  select(Date, TridentTemp)

IncrementDate <- read.csv("D16.Final21Days_withDate.csv") %>%
  mutate(Date = ymd(Date))%>%
  mutate(IncrementNumber = case_when(
   Increment_21Edge == 'increment1' ~ 1,
   Increment_21Edge == 'increment2' ~ 2,
   Increment_21Edge == 'increment3' ~ 3,
   Increment_21Edge == 'increment4' ~ 4,
   Increment_21Edge == 'increment5' ~ 5,
   Increment_21Edge == 'increment6' ~ 6,
   Increment_21Edge == 'increment7' ~ 7,
   Increment_21Edge == 'increment8' ~ 8,
   Increment_21Edge == 'increment9' ~ 9,
   Increment_21Edge == 'increment10' ~ 10,
   Increment_21Edge == 'increment11' ~ 11,
   Increment_21Edge == 'increment12' ~ 12,
   Increment_21Edge == 'increment13' ~ 13,
   Increment_21Edge == 'increment14' ~ 14,
   Increment_21Edge == 'increment15' ~ 15,
   Increment_21Edge == 'increment16' ~ 16,
   Increment_21Edge == 'increment17' ~ 17,
   Increment_21Edge == 'increment18' ~ 18,
   Increment_21Edge == 'increment19' ~ 19,
   Increment_21Edge == 'increment20' ~ 20,
   Increment_21Edge == 'increment21' ~ 21
  )) 
#Note: sometimes has issues loading, data file may have converted dates to different format and need to be converted back to a year-month-day formate

Trident_byFish <- merge(IncrementDate, TridentTemps, by = "Date") %>%
  select(FISHID, Date, IncrementNumber, TridentTemp)


##Habitat and Abundance Files 
Seine <- read.csv("D2.SeineData.csv") %>%
  mutate(DateCaptured = ymd(Date)) %>%
  select(SeineID, SITEID, DateCaptured, No.PacificCod)
Sampling <- read.csv("D10.SamplingData_Master.csv") %>%
  mutate(Yearf = as.factor(Year))

Habitat <- Sampling %>%
  left_join(Seine) 

##Fish Dissection Files
all.fish <- read.csv("D6.Size_Condition_ROutput.csv") %>%
  select(FISHID, SL_mm, WholeBodyWW_g, LWResiduals_All, HSI)

## Merge Data for Growth Model
AllGrowth <- Growth %>%
  select(FISHID, IncrementNumber, Size_Initial_BIC, mm_day_BIC, mm_mm_day_BIC) %>%
  left_join(Trident_byFish, by = c('FISHID', 'IncrementNumber'))%>%
  left_join(Habitat, by = "FISHID") %>%
  left_join(all.fish, by = 'FISHID') %>%
  left_join(AllDiet, by = 'FISHID') %>%
  left_join(Fullness, by = 'FISHID') 

## Subset by July and August

JulyGrowth <- AllGrowth[AllGrowth$Month ==7,]
AugGrowth <- AllGrowth[AllGrowth$Month ==8,]

#### July Growth Model #####

##### Apply Centering Function and Scaling #####

c. <- function (x) scale(x, scale = TRUE) 

JulyGrowth_Model <- JulyGrowth %>%
  mutate(log_mm_mm_day = log(mm_mm_day_BIC))%>%
  mutate(scaledTemp = c.(TridentTemp)) %>%
  mutate(scaledSize = c.(SL_mm))%>%
  mutate(scaledNMDS1 = c.(NMDS1)) %>%
  mutate(scaledPcodAbundance = c.(No.PacificCod)) %>%
  mutate(scaledFullness = c.(sqrtFullness)) %>%
  mutate(Heatwave = as.factor(Heatwave)) %>%
  select(FISHID, IncrementNumber, log_mm_mm_day, scaledTemp, scaledSize, scaledNMDS1, scaledFullness, scaledPcodAbundance, Heatwave, Year)


#### Correlation between Variables

cor(JulyGrowth_Model$log_mm_mm_day, JulyGrowth_Model$scaledTemp, use="complete.obs")  #r = -0.33
cor(JulyGrowth_Model$log_mm_mm_day, JulyGrowth_Model$scaledSize, use="complete.obs")  #r = -0.41
cor(JulyGrowth_Model$log_mm_mm_day, JulyGrowth_Model$scaledNMDS1, use="complete.obs") # r = -0.25
cor(JulyGrowth_Model$log_mm_mm_day, JulyGrowth_Model$scaledPcodAbundance, use="complete.obs") # r = -0.017
cor(JulyGrowth_Model$log_mm_mm_day, JulyGrowth_Model$scaledFullness, use="complete.obs") # r = -0.022
cor(JulyGrowth_Model$scaledSize, JulyGrowth_Model$scaledTemp, use="complete.obs")  #r = 0.38
cor(JulyGrowth_Model$scaledPcodAbundance, JulyGrowth_Model$scaledSize, use="complete.obs") # r = -0.15
cor(JulyGrowth_Model$scaledPcodAbundance, JulyGrowth_Model$scaledTemp, use="complete.obs") # r = -0.31
cor(JulyGrowth_Model$scaledNMDS1, JulyGrowth_Model$scaledSize, use="complete.obs") # r = 0.58
cor(JulyGrowth_Model$scaledNMDS1, JulyGrowth_Model$scaledTemp, use="complete.obs") # r = 0.19
cor(JulyGrowth_Model$scaledNMDS1, JulyGrowth_Model$scaledPcodAbundance, use="complete.obs") # r = 0.26
cor(JulyGrowth_Model$scaledNMDS1, JulyGrowth_Model$scaledFullness, use="complete.obs") # r = 0.29

#### July Model Selection ####

##### Optimal Random Structure ####

#No Random Effects, all Fixed Effects and Interactions
M1 <- lm(log_mm_mm_day ~ scaledSize * scaledTemp * scaledNMDS1 * scaledFullness  * scaledPcodAbundance * Heatwave, na.action = na.exclude, JulyGrowth_Model) 
summary(M1)
anova(M1)


#Fish ID as a random Intercept
M2<- lmer(log_mm_mm_day ~ scaledSize *  scaledTemp * scaledNMDS1 * scaledFullness *  scaledPcodAbundance * Heatwave + (1|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(M2) 
anova(M2)
Anova(M2)


# Day of Life as a random slope and Fish ID as a random intercept
M3<- lmer(log_mm_mm_day ~ scaledSize *  scaledPcodAbundance * scaledTemp * scaledNMDS1 * scaledFullness * Heatwave  + (IncrementNumber|FISHID), na.action = na.exclude,JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(M3)
Anova(M3)


#model comparison 
#AIC with REML
AIC(M1, M2, M3)

anova(M2, M3) #AIC with ML
#Still Model 3

r.squaredGLMM(M2)
r.squaredGLMM(M3)

##### Optimal Fixed Structure ####

#Full Model
MA<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp * scaledNMDS1 * scaledFullness * Heatwave * scaledPcodAbundance + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MA)
Anova(MA)


#Remove Abundance
MB <- lmer(log_mm_mm_day ~  scaledSize * scaledTemp * scaledNMDS1 *  scaledFullness * Heatwave   + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MB) 
anova(MB)
Anova(MB)
AIC(MA, MB)
# Model w/out abundance is better

#Remove NMDS Axis 1
MC <- lmer(log_mm_mm_day ~  scaledSize * scaledTemp *   scaledFullness * Heatwave   + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MC) 
anova(MC)
Anova(MC)
AIC(MA, MB, MC)
#Model w/out NMS is better

#Remove stomach fullness
MD <- lmer(log_mm_mm_day ~  scaledSize * scaledTemp *  Heatwave   + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MD) 
anova(MD)
Anova(MD)
AIC(MA, MB, MC, MD)


#Size and Temp

ME <- lmer(log_mm_mm_day ~  scaledSize * scaledTemp + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(ME) 
anova(ME)
Anova(ME)
AIC(MA, MB, MC, MD, ME)


#Size and HW

MF <- lmer(log_mm_mm_day ~  scaledSize * Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MF) 
anova(MF)
Anova(MF)
AIC(MA, MB, MC, MD, ME, MF)


#HW and Temp

MG <- lmer(log_mm_mm_day ~  scaledTemp * Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MG) 
anova(MG)
Anova(MG)
AIC(MA, MB, MC, MD, ME, MF, MG)


#HW Only

MH <- lmer(log_mm_mm_day ~  Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MH) 
anova(MH)
Anova(MH)
AIC(MA, MB, MC, MD, ME, MF, MG, MH)

#Size Only

MI <- lmer(log_mm_mm_day ~  scaledSize + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MI) 
anova(MI)
Anova(MI)
AIC(MA, MB, MC, MD, ME, MF, MG, MH, MI)


#Temp Only

MJ <- lmer(log_mm_mm_day ~  scaledTemp + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MJ) 
anova(MJ)
Anova(MJ)
AIC(MA, MB, MC, MD, ME, MF, MG, MH, MI, MJ)



#Size, Temp, HW, no interactions
MK <- lmer(log_mm_mm_day ~  scaledSize + scaledTemp +  Heatwave   + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(MK) 
anova(MK)
Anova(MK)
AIC(MA, MB, MC, MD, ME, MF, MG, MH, MI, MJ, MK)


##### Optimal Error Structure ####
#Fit with lmer()
JulySize_reml<- lmer(log_mm_mm_day ~ scaledSize + scaledTemp  +  Heatwave  + (IncrementNumber|FISHID), na.action = na.exclude, JulyGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML = TRUE)
summary(JulySize_reml)
anova(JulySize_reml)
Anova(JulySize_reml)

r.squaredGLMM(JulySize_reml)

#Fit with lme()
JulySize_reml2 <- lme(log_mm_mm_day ~ scaledSize + scaledTemp + Heatwave , method = "REML",
                      JulyGrowth_Model, na.action = na.exclude, random = ~IncrementNumber|FISHID)
summary(JulySize_reml2) 
anova(JulySize_reml2)
Anova(JulySize_reml2)
r.squaredGLMM(JulySize_reml2)

##Plot Autocorrelation
ACF(JulySize_reml2,maxLag = 10) 
plot(ACF(JulySize_reml2), alpha = 0.05, resType = "normalized")

## Autocorrelation error structure
JulySize_reml_acor <- update(JulySize_reml2, 
                             correlation=corAR1(form=~IncrementNumber|FISHID))
plot(ACF(JulySize_reml_acor, resType = "normalized"), alpha = 0.05)
summary(JulySize_reml_acor)
anova(JulySize_reml_acor)
Anova(JulySize_reml_acor)
r.squaredGLMM(JulySize_reml_acor)

AIC(JulySize_reml2, JulySize_reml_acor)

#R2m = 0.21
#R2c = 0.691
#Phi is an estimate of the lag 1 autocorrelation and is telling us that the residuals are correlated at lag 1, with a correlation coefficient (recall Pearson’s) of 0.27



#### July Model Validation #####


res<-residuals(JulySize_reml_acor,type='normalized')

hist(res,main='July Normalized residuals',
     xlab="Normalized residuals")

#Check Homogeneity of Variance
plot(log(JulyGrowth_Model$scaledSize),res,main='Homogeneity of variance: Size',
     xlab='Scaled Size',
     ylab='Normalized residuals')
abline(h=0,lty=2)
plot(log(JulyGrowth_Model$scaledTemp),res,main='Homogeneity of variance: Temp',
     xlab='Scaled Temp',
     ylab='Normalized residuals')
abline(h=0,lty=2)


#Check Autocorrelation
plot(ACF(JulySize_reml_acor, resType = "normalized"), alpha = 0.05)




#Check normality
qqnorm(resid(JulySize_reml_acor, type = "normalized"))
qqline(resid(JulySize_reml_acor, type = "normalized"))


#Check Fitted Values
plot(predict(JulySize_reml_acor),res,main = "July Fitted Values",
     ylab='Normalized residuals',
     xlab='Fitted values')
abline(h=0,lty=2)


#Check Independence and spread of variance
boxplot(res~JulyGrowth_Model$Heatwave,main='July Independence',
        xlab='Heatwave',ylab='Normalized residuals')


## Overall, the Model fits all model assumptions


#### July Predicted Growth Plots  #####

JulyGrowth_Model$predictedSizeModel <- predict(JulySize_reml_acor)

#Unscale Scaled Variables
JulyGrowth_Model  <- JulyGrowth_Model %>%
  mutate(temp_unscaled = scaledTemp * attr(scaledTemp, 'scaled:scale') + attr(scaledTemp, 'scaled:center')) %>%
  mutate(size_unscaled = scaledSize * attr(scaledSize, 'scaled:scale') + attr(scaledSize, 'scaled:center'))

#July Predicted Growth against Temperature (Figure 4a)
JulyTemp <- ggplot(NULL, mapping = aes(temp_unscaled, exp(predictedSizeModel))) +
  geom_point(data = JulyGrowth_Model, aes(x = temp_unscaled, y = exp(log_mm_mm_day), fill = Heatwave), alpha = 1, size =3, shape =21, stroke = 0)+
  geom_smooth(data =JulyGrowth_Model, aes(x = temp_unscaled, y = exp(predictedSizeModel), color = Heatwave), linewidth = 1.5, method = lm, se = TRUE )+
  theme_classic(base_size =20, base_family = "") +
  scale_color_manual(values=c("#402aaa", "grey1", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea", "#888888", "#e5b2bb")) +
  xlab("Temperature (ºC)") +
  ylab("Predicted Growth mm/mm/day") +
  ggtitle("July") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(JulyTemp)

#July Predicted Growth against Size (Figure 4c)
JulySize <- ggplot(NULL, mapping = aes(size_unscaled, exp(predictedSizeModel))) +
  geom_point(data = JulyGrowth_Model, aes(x = size_unscaled, y = exp(log_mm_mm_day), fill = Heatwave), alpha = 1, size =3, shape =21, stroke = 0)+
  geom_smooth( data =JulyGrowth_Model, aes(x = size_unscaled, y = exp(predictedSizeModel), color = Heatwave), linewidth = 1.5, method = lm, se = TRUE)+
  theme_classic(base_size =20, base_family = "") +
  scale_color_manual(values=c("#402aaa", "grey1", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea", "#888888", "#e5b2bb")) +
  xlab("Size (mm)") +
  ylab("Predicted Growth mm/mm/day")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"))
print(JulySize)


#### July Effects Plots ####

mydf_JulySize_LMM <- ggemmeans(JulySize_reml_acor, terms = c("scaledSize", "scaledTemp", "Heatwave"))
mydf_JulySize_LMM

plot(mydf_JulySize_LMM)+
  theme_classic(base_size = 17, base_family = "")

July_SizeModel_byHW <- as.data.frame(mydf_JulySize_LMM) %>%
  mutate(Heatwave = facet) %>%
  select(!facet)
labels<- c('Before' = "Before" ,'Between' = "Between" , 'Heatwave' = "Heatwave")

#Calculate scaling values to include in model
attr(JulyGrowth_Model$scaledSize, 'scaled:scale') #11.52
attr(JulyGrowth_Model$scaledSize, 'scaled:center') # 49.73
attr(JulyGrowth_Model$scaledTemp, 'scaled:scale') #1.08
attr(JulyGrowth_Model$scaledTemp, 'scaled:center') # 9.04

##### Marginal Means (Supplemental Table S7) #####
July_MarginalMeans_All <- July_SizeModel_byHW  %>%
  mutate(Size = x * 11.52 + 49.73) %>%
  mutate(RelativeGrowth_Mean = exp(predicted)) %>%
  mutate(RelativeGrowth_SE = exp(std.error)) %>%
  mutate(RelativeGrowth_conf.low = exp(conf.low)) %>%
  mutate(RelativeGrowth_conf.high = exp(conf.high)) %>%
  mutate(Temp_Unscaled = as.numeric(group)) %>%
  mutate(Temp = (Temp_Unscaled - 2) * 1.08 + 9.04) %>%
  mutate(Temperature = ((as.numeric(group) - 2) *1.08) + 9.04) %>%
  select(Size, Temperature, Heatwave, RelativeGrowth_Mean, RelativeGrowth_SE, RelativeGrowth_conf.low, RelativeGrowth_conf.high)

## Trim size lines that extend beyond data

#Min and Max Sizes
JulyGrowth_Model %>% 
  group_by(Heatwave) %>%
  reframe(min_SL = min(size_unscaled), max.SL = max(size_unscaled))

July_SizeModel_byHW_TRIM_Before <- July_SizeModel_byHW %>%
  mutate(UnscaledSize = x * 11.52 + 49.73) %>%
  filter(Heatwave == "Before") %>%
  filter(UnscaledSize <= 67.01) %>%
  filter(UnscaledSize >= 31) %>%
  mutate(TempNumeric = as.numeric(group)) %>%
  mutate(UnscaledTemp = (TempNumeric - 2) * 1.08 + 9.04) %>%
  mutate(TempFactor = as.factor(UnscaledTemp))

July_SizeModel_byHW_TRIM_Between<- July_SizeModel_byHW %>%
  mutate(UnscaledSize = x * 11.52 + 49.73) %>%
  filter(Heatwave == "Between") %>%
  filter(UnscaledSize <= 96) %>%
  filter(UnscaledSize >= 32) %>%
  mutate(TempNumeric = as.numeric(group)) %>%
  mutate(UnscaledTemp = (TempNumeric - 2) * 1.08 + 9.04) %>%
  mutate(TempFactor = as.factor(UnscaledTemp))

July_SizeModel_byHW_TRIM_HW <- July_SizeModel_byHW %>%
  mutate(UnscaledSize = x * 11.52 + 49.73) %>%
  filter(Heatwave == "Heatwave") %>%
  filter(UnscaledSize <= 79) %>%
  filter(UnscaledSize >= 32) %>%
  mutate(TempNumeric = as.numeric(group)) %>%
  mutate(UnscaledTemp = (TempNumeric - 2) * 1.08 + 9.04) %>%
  mutate(TempFactor = as.factor(UnscaledTemp))

July_SizeModel_byHW_TRIM <- rbind(July_SizeModel_byHW_TRIM_Before, July_SizeModel_byHW_TRIM_Between, July_SizeModel_byHW_TRIM_HW) 

#July Marginal Means Plot (Figure 5a)
JulyEffects <- ggplot() +
  geom_point(data = JulyGrowth_Model, aes(x =size_unscaled, y = exp(log_mm_mm_day)),color = "grey30", alpha = 0.1)+
  geom_smooth(data =July_SizeModel_byHW_TRIM, aes(x = UnscaledSize, y = exp(predicted), color = TempFactor), linewidth = 1.5)+
  facet_wrap( ~ Heatwave) +
  geom_ribbon(data = July_SizeModel_byHW_TRIM, aes(x = UnscaledSize, y = exp(predicted), ymin = exp(conf.low), ymax = exp(conf.high), fill = TempFactor, group = TempFactor), alpha = 0.25, show.legend = FALSE) +
  theme_classic(base_size =20, base_family = "") +
  scale_color_viridis_d(option="inferno", labels = c("7.96", "9.04", "10.12")) +
  scale_fill_viridis_d(option="inferno", labels = c("7.96", "9.04", "10.12")) +
  coord_cartesian(ylim = c(0, 0.04)) + 
  theme(axis.text.x = element_text(size = 12, angle = 45),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle("July") +
  guides(color = guide_legend(title = "Temp (ºC)"))+
  #xlim(75, 150) +
  ylab("Predicted Growth mm/mm/day") +
  xlab("Size (mm)") 
print(JulyEffects)

#### August Growth Model ####

###### Apply Centering Function and Scaling ######

c. <- function (x) scale(x, scale = TRUE) 

AugGrowth_Model <- AugGrowth %>%
  mutate(log_mm_mm_day = log(mm_mm_day_BIC))%>%
  mutate(scaledTemp = c.(TridentTemp)) %>%
  mutate(scaledSize = c.(SL_mm))%>%
  mutate(scaledNMDS1 = c.(NMDS1)) %>%
  mutate(scaledPcodAbundance = c.(No.PacificCod)) %>%
  mutate(scaledFullness = c.(sqrtFullness)) %>%
  mutate(Heatwave = as.factor(Heatwave)) %>%
  select(FISHID, IncrementNumber, log_mm_mm_day, scaledTemp, scaledSize, scaledNMDS1, scaledFullness, scaledPcodAbundance, Heatwave, Year)

## Correlation between Variables

cor(AugGrowth_Model$log_mm_mm_day, AugGrowth_Model$scaledTemp, use="complete.obs")  #r = -0.25
cor(AugGrowth_Model$log_mm_mm_day, AugGrowth_Model$scaledSize, use="complete.obs")  #r = -0.58
cor(AugGrowth_Model$log_mm_mm_day, AugGrowth_Model$scaledNMDS1, use="complete.obs") # r = -0.41
cor(AugGrowth_Model$log_mm_mm_day, AugGrowth_Model$scaledPcodAbundance, use="complete.obs") # r = -0.09
cor(AugGrowth_Model$log_mm_mm_day, AugGrowth_Model$scaledFullness, use="complete.obs") # r = 0.069
cor(AugGrowth_Model$scaledSize, AugGrowth_Model$scaledTemp, use="complete.obs")  #r = 0.54
cor(AugGrowth_Model$scaledPcodAbundance, AugGrowth_Model$scaledSize, use="complete.obs") # r = -0.01
cor(AugGrowth_Model$scaledPcodAbundance, AugGrowth_Model$scaledTemp, use="complete.obs") # r = -0.07
cor(AugGrowth_Model$scaledNMDS1, AugGrowth_Model$scaledSize, use="complete.obs") # r = 0.69
cor(AugGrowth_Model$scaledNMDS1, AugGrowth_Model$scaledTemp, use="complete.obs") # r = 0.57
cor(AugGrowth_Model$scaledNMDS1, AugGrowth_Model$scaledPcodAbundance, use="complete.obs") # r = =0.01

#### August Model Selection ####

##### Optimal Random Structure August Model ####

#No Random Effects, all Fixed Effects and Interactions
N1 <- lm(log_mm_mm_day ~ scaledSize *  scaledTemp * scaledNMDS1 * scaledFullness * scaledPcodAbundance * Heatwave, na.action = na.exclude, AugGrowth_Model) 
summary(N1)
anova(N1)


#Fish ID as a random Intercept
N2<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp * scaledNMDS1 * scaledFullness * scaledPcodAbundance * Heatwave + (1|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(N2) 
anova(N2)
Anova(N2)

# Day of Life as a random slope and Fish ID as a random intercept
N3<- lmer(log_mm_mm_day ~ scaledSize  * scaledTemp * scaledNMDS1 * scaledFullness * scaledPcodAbundance * Heatwave + (IncrementNumber|FISHID), na.action = na.exclude,AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(N3)
Anova(N3)

AIC(N1,N2, N3)

anova(N2, N3)
#Model 3 best 

r.squaredGLMM(N2)
r.squaredGLMM(N3)

##### Optimal Fixed Structure ####

#Full Model
NA1<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp * scaledNMDS1 * scaledFullness * scaledPcodAbundance * Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NA1)
Anova(NA1)

#Remove Abundance
#NMDS Axis 1, Fullness, Size, Temp, HW
NB<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp * scaledNMDS1 * scaledFullness *  Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NB)
Anova(NB)
AIC(NA1,NB)

#Remove NMS Axis 1
#Stomach Fullness, Size, Temp, HW
NC<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp * scaledFullness *   Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NC)
Anova(NC)
AIC(NA1,NB, NC)

#Remove Fullness
#Size, Temp, HW
ND<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp *    Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(ND)
Anova(ND)
AIC(NA1,NB, NC, ND)

#Remove HW
#Size, Temp
NE<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp *  (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NE)
Anova(NE)
AIC(NA1,NB, NC, ND, NE)

#Remove Temp
#HW, Size
NF<- lmer(log_mm_mm_day ~ Heatwave * scaledSize  + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NF)
Anova(NF)
AIC(NA1,NB, NC, ND, NE, NF)

#Remove Size
#HW, Temp
NG<- lmer(log_mm_mm_day ~ Heatwave * scaledTemp  + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NG)
Anova(NG)
AIC(NA1,NB, NC, ND, NE, NF, NG)

#HW Only
NH <- lmer(log_mm_mm_day ~ Heatwave  + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NH)
Anova(NH)
AIC(NA1,NB, NC, ND, NE, NF, NG, NH)

#Size Only
NI <- lmer(log_mm_mm_day ~ scaledSize  + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NI)
Anova(NI)
AIC(NA1,NB, NC, ND, NE, NF, NG, NH, NI)

#Temp Only
NJ <- lmer(log_mm_mm_day ~ scaledTemp  + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NJ)
Anova(NJ)
AIC(NA1,NB, NC, ND, NE, NF, NG, NH, NI, NJ)

#Size, Temp, HW, no Interactions
NK <- lmer(log_mm_mm_day ~ scaledTemp +scaledSize + Heatwave + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML =FALSE) 
summary(NK)
Anova(NK)
AIC(NA1,NB, NC, ND, NE, NF, NG, NH, NI, NJ, NK)

##### Optimal Error Structure ####

#Fit with lmer()
AugSize_reml<- lmer(log_mm_mm_day ~ scaledSize * scaledTemp  *  Heatwave  + (IncrementNumber|FISHID), na.action = na.exclude, AugGrowth_Model, control = lmerControl(optimizer = "Nelder_Mead"), REML = TRUE)
summary(AugSize_reml)
anova(AugSize_reml)
Anova(AugSize_reml)
r.squaredGLMM(AugSize_reml)

#Fit with lme()
AugSize_reml2 <- lme(log_mm_mm_day ~ scaledSize * scaledTemp * Heatwave , method = "REML",
                     AugGrowth_Model, na.action = na.exclude, random = ~IncrementNumber|FISHID)

summary(AugSize_reml2) 
anova(AugSize_reml2)
Anova(AugSize_reml2)
r.squaredGLMM(AugSize_reml2)

##Plot Autocorrelation
ACF(AugSize_reml2, maxLag = 10) 
plot(ACF(AugSize_reml2), alpha = 0.05, resType = "normalized")

## Autocorrelation error structure
AugSize_reml_acor <- update(AugSize_reml2, 
                            correlation=corAR1(form=  ~IncrementNumber|FISHID))
plot(ACF(AugSize_reml_acor, resType = "normalized"), alpha = 0.05)
summary(AugSize_reml_acor)
anova(AugSize_reml_acor)
Anova(AugSize_reml_acor)
r.squaredGLMM(AugSize_reml_acor)

#R2m = 0.424
#R2c = 0.861

#Phi is an estimate of the lag 1 autocorrelation and is telling us that the residuals are correlated at lag 1, with a correlation coefficient (recall Pearson’s) of 0.34

# The Random Effects:StdDev is the standard deviation of the random terms associated with each coefficient as well as the residual.  Small StdDev values indicate that there is not a very strong difference among the individual fish. The Corr value is the correlation between the random terms of each coefficient. 

AIC(AugSize_reml2, AugSize_reml_acor)

#### August Model validation ####

res1<-residuals(AugSize_reml_acor,type='normalized')

hist(res1,main='August Normalized residuals',
     xlab="Normalized residuals")

#Check normality
qqnorm(resid(AugSize_reml_acor, type = "normalized"))
qqline(resid(AugSize_reml_acor, type = "normalized"))

#Check spread of residuals
plot(predict(AugSize_reml_acor),res1,main='August Fitted values',
     ylab='Normalized residuals',
     xlab='Fitted values')
abline(h=0,lty=2)
plot(log(AugGrowth_Model$scaledSize),res1,main='Homogeneity of Variance: Size',
     xlab='Scaled Age',
     ylab='Normalized residuals')
abline(h=0,lty=2)
plot(log(AugGrowth_Model$scaledTemp),res1,main='Homogeneity of Variance',
     xlab='Scaled Temp',
     ylab='Normalized residuals')
abline(h=0,lty=2)

#Check Independence and spread of variance
boxplot(res1~AugGrowth_Model$Heatwave,main='August Independence',
        xlab='Heatwave',ylab='Normalized residuals')

### Overall, the Model fits all model assumptions



#### August Predicted Growth Plots ####

AugGrowth_Model$predictedSizeModel <- predict(AugSize_reml_acor)

#Unscale scaled variables

AugGrowth_Model  <- AugGrowth_Model %>%
  mutate(temp_unscaled = scaledTemp * attr(scaledTemp, 'scaled:scale') + attr(scaledTemp, 'scaled:center')) %>%
  mutate(size_unscaled = scaledSize * attr(scaledSize, 'scaled:scale') + attr(scaledSize, 'scaled:center'))

#August Predicted Growth against Temp (Figure 4b)
AugTemp <- ggplot(NULL, mapping = aes(temp_unscaled, exp(predictedSizeModel))) +
  geom_point(data = AugGrowth_Model, aes(x = temp_unscaled, y = exp(log_mm_mm_day), fill = Heatwave), alpha = 1, size =3, shape =21, stroke = 0)+
  geom_smooth(data =AugGrowth_Model, aes(x = temp_unscaled, y = exp(predictedSizeModel), color = Heatwave), linewidth = 1.5, method = lm, se = TRUE )+
  theme_classic(base_size =20, base_family = "") +
  scale_color_manual(values=c("#402aaa", "grey1", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea", "#888888", "#e5b2bb")) +
  xlab("Temperature (ºC)") +
  ylab("Predicted Growth mm/mm/day") +
  ggtitle("August") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(AugTemp)


#August Predicted Growth Against Size (Figure 4d)
AugSize <- ggplot(NULL, mapping = aes(size_unscaled, exp(predictedSizeModel))) +
  geom_point(data = AugGrowth_Model, aes(x = size_unscaled, y = exp(log_mm_mm_day), fill = Heatwave), alpha = 1, size =3, shape =21, stroke = 0)+
  geom_smooth(data =AugGrowth_Model, aes(x = size_unscaled, y = exp(predictedSizeModel), color = Heatwave), linewidth = 1.5, method = lm, se = TRUE )+
  theme_classic(base_size =20, base_family = "") +
  scale_color_manual(values=c("#402aaa", "grey1", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea", "#888888", "#e5b2bb")) +
  xlab("Size (mm)") +
  ylab("Predicted Growth mm/mm/day") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"))
print(AugSize)

#### August Effects Plots ####


mydf_AugSize_LMM <- ggemmeans(AugSize_reml_acor, terms = c("scaledSize", "scaledTemp", "Heatwave"))
mydf_AugSize_LMM

Aug_SizeTempHW_byHW <- as.data.frame(mydf_AugSize_LMM) %>%
  mutate(Heatwave = facet) %>%
  select(!facet)

labels<- c('Before' = "Before" ,'Between' = "Between" , 'Heatwave' = "Heatwave")


#Calculate scaling values to include in model
attr(AugGrowth_Model$scaledSize, 'scaled:scale') #16.92
attr(AugGrowth_Model$scaledSize, 'scaled:center') # 75.60
attr(AugGrowth_Model$scaledTemp, 'scaled:scale') #1.324
attr(AugGrowth_Model$scaledTemp, 'scaled:center') # 10.27

##### Marginal Means (Table S7) #####
Aug_MarginalMeans_All <- Aug_SizeTempHW_byHW %>%
  mutate(Size = x * 16.92 + 75.6) %>%
  mutate(RelativeGrowth_Mean = exp(predicted)) %>%
  mutate(RelativeGrowth_SE = exp(std.error)) %>%
  mutate(RelativeGrowth_conf.low = exp(conf.low)) %>%
  mutate(RelativeGrowth_conf.high = exp(conf.high)) %>%
  mutate(Temperature = ((as.numeric(group) - 2) *1.324) + 10.27) %>%
  select(Size, Temperature, Heatwave, RelativeGrowth_Mean, RelativeGrowth_SE, RelativeGrowth_conf.low, RelativeGrowth_conf.high)
#write.csv(Aug_MarginalMeans_All, "AugustMarginalMeans_All.csv")

## Trim size lines that extend beyond data

#Min and Max Sizes
AugGrowth_Model %>% 
  group_by(Heatwave) %>%
  reframe(min_SL = min(size_unscaled), max.SL = max(size_unscaled))

Aug_SizeModel_byHW_TRIM_Before <- Aug_SizeTempHW_byHW %>%
  mutate(UnscaledSize = x * 16.92 + 75.6) %>%
  filter(Heatwave == "Before") %>%
  filter(UnscaledSize <= 103) %>%
  filter(UnscaledSize >= 41) %>%
  mutate(TempNumeric = as.numeric(group)) %>%
  mutate(UnscaledTemp = (TempNumeric - 2) * 1.324 + 10.27) %>%
  mutate(TempFactor = as.factor(UnscaledTemp))

Aug_SizeModel_byHW_TRIM_Between<- Aug_SizeTempHW_byHW %>%
  mutate(UnscaledSize = x * 16.92 + 75.6) %>%
  filter(Heatwave == "Between") %>%
  filter(UnscaledSize <= 135) %>%
  filter(UnscaledSize >= 50) %>%
  mutate(TempNumeric = as.numeric(group)) %>%
  mutate(UnscaledTemp = (TempNumeric - 2) * 1.324 + 10.27) %>%
  mutate(TempFactor = as.factor(UnscaledTemp))

Aug_SizeModel_byHW_TRIM_HW <- Aug_SizeTempHW_byHW %>%
  mutate(UnscaledSize = x * 16.92 + 75.6) %>%
  filter(Heatwave == "Heatwave") %>%
  filter(UnscaledSize <= 135) %>%
  filter(UnscaledSize >= 67) %>%
  mutate(TempNumeric = as.numeric(group)) %>%
  mutate(UnscaledTemp = (TempNumeric - 2) * 1.324 + 10.27) %>%
  mutate(TempFactor = as.factor(UnscaledTemp))

Aug_SizeModel_byHW_TRIM <- rbind(Aug_SizeModel_byHW_TRIM_Before, Aug_SizeModel_byHW_TRIM_Between, Aug_SizeModel_byHW_TRIM_HW)


#August Marginal Means Plot (Figure 5b)
AugEffects <- ggplot() +
  geom_point(data = AugGrowth_Model, aes(x =size_unscaled, y = exp(log_mm_mm_day)), color = "grey30", alpha = 0.1)+
  geom_smooth(data =Aug_SizeModel_byHW_TRIM, aes(x = UnscaledSize, y = exp(predicted), color = group), linewidth = 1.5)+
  facet_wrap( ~ Heatwave) +
  coord_cartesian(ylim = c(0, 0.04)) + 
  geom_ribbon(data = Aug_SizeModel_byHW_TRIM, aes(x = UnscaledSize, y = exp(predicted), ymin = exp(conf.low), ymax = exp(conf.high), fill = group, group = group), alpha = 0.25, show.legend = FALSE) +
  theme_classic(base_size =20, base_family = "") +
  scale_color_viridis_d(option="inferno", labels = c("8.95", "10.27", "11.59")) +
  scale_fill_viridis_d(option="inferno", labels = c("8.95", "10.27", "11.59")) +
  theme(axis.text.x = element_text(size = 12, angle = 45),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle("August") +
  guides(color = guide_legend(title = "Temp (ºC)"))+
  #xlim(75, 150) +
  ylab("Predicted Growth mm/mm/day") +
  xlab("Size (mm)")
print(AugEffects)

#### Manuscript Figures ####

#Predicted Growth Plots (Figure 4)
PredictedGrowth_Plots <- JulyTemp + AugTemp + JulySize + AugSize +
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(PredictedGrowth_Plots)
#dev.copy(jpeg,'Figure_4.jpg', width=14, height=10, units='in', res=300)
#dev.off()

#Marginal Means Plots (Figure 5)
GrowthEffects_Plots <- JulyEffects + AugEffects + 
  plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(GrowthEffects_Plots)
#dev.copy(jpeg,'Figure_5.jpg', width=8, height=8, units='in', res=300)
#dev.off()
