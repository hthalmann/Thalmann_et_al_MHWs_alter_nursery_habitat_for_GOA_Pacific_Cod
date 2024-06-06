#Predicted August Size Based on July Size and Growth
#Model simulations for July individuals growing out to August


library(rstudioapi) #load R working directory
library(tidyverse)  #data wrangling and visualization
library(lme4) #linear models
library(car) #ANOVA function for models
library(MuMIn) #Runs the R.squaredGLMM function
library(ggpmisc) #LM equations on plot
library(patchwork) #data visualization

##Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#Sets directory to folder that this R script is saved to

## Data Summaries ##
MeanGrowth <- read.csv("D12.BackCalculatedGrowth_ROutput.csv") %>%
  select(FISHID, Year, Month, mm_day_BIC, mm_mm_day_BIC) %>%
  group_by(FISHID) %>%
  summarise_all(list(mean = ~mean(., na.rm = T))) %>%
  mutate(mm_mm_day = mm_mm_day_BIC_mean) %>%
  mutate(mm_day = mm_day_BIC_mean) %>%
  select( -Year_mean, -Month_mean, -mm_mm_day_BIC_mean, -mm_day_BIC_mean) %>%ungroup()
Growth_Size <- read.csv("D6.Size_Condition_ROutput.csv") %>%
  select(FISHID,  SL_mm, Year, Month, Heatwave) %>%
  mutate(Year = as.factor(Year)) %>%
  left_join(MeanGrowth) 


Growth_Size_Summary <- Growth_Size %>%
  select(-Year, -FISHID) %>% 
group_by(Month,Heatwave) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

#Number of Days between sampling
Seine <- read.csv("D2.SeineData.csv") %>%
  mutate(DateCaptured = ymd(Date)) %>%
  select(SeineID, SITEID, DateCaptured, No.PacificCod)

Sampling <- read.csv("D10.SamplingData_Master.csv") %>%
  mutate(Year = as.factor(Year)) %>%
  left_join(Seine) 

JulySamplingSummary <- Sampling %>%
  filter(Month == '7') %>%
  group_by(Year, Heatwave) %>%
  summarise_at(vars(DateCaptured), list(MeanJulyCaptureDate = mean), na.rm = T) %>%
  mutate(MeanJulyCaptureDate_Julian = yday(MeanJulyCaptureDate)) 

MeanSamplingDate <- read.csv("D10.SamplingData_Master.csv") %>%
  mutate(Year = as.factor(Year)) %>%
  left_join(Seine) %>%
  filter(Month == 8) %>%
  group_by(Year, Heatwave) %>%
  summarise_at(vars(DateCaptured), list(MeanAugCaptureDate = mean), na.rm = T) %>%
  mutate(MeanAugCaptureDate_Julian = yday(MeanAugCaptureDate)) %>%
  left_join(JulySamplingSummary) %>%
  mutate(Days_Between_Sampling = MeanAugCaptureDate_Julian -MeanJulyCaptureDate_Julian) 

MeanSamplingDate_Diff <- MeanSamplingDate %>%
  group_by(Heatwave) %>%
  summarise_at(vars(Days_Between_Sampling), list(MeanDays_Between_Sampling = mean), na.rm = T)



####Predicted Size in August based on July sizes and growth rates####

AugustPredictedSize <- Growth_Size  %>%
  left_join(Sampling) %>%
  filter(Month == 7) %>%
  mutate(JulianDateCaptured = yday(DateCaptured)) %>%
  left_join(MeanSamplingDate) %>%
  select(FISHID, Year, Heatwave, mm_mm_day, SL_mm, JulianDateCaptured, MeanAugCaptureDate_Julian) %>% 
  drop_na(mm_mm_day) %>%
  mutate(DaysbetweenCapture = MeanAugCaptureDate_Julian - JulianDateCaptured) %>%
  mutate(mm_day = mm_mm_day * SL_mm) %>%
  mutate(mm_betweenJulyAug = mm_day * DaysbetweenCapture) %>%
  mutate(EstAugSize_fromJulyGrowth = SL_mm + mm_betweenJulyAug) %>%
  mutate(Group = "PredictedAugustSize") %>%
  ungroup() %>%
  select(FISHID, Year,  Heatwave, Group, EstAugSize_fromJulyGrowth) %>%
  mutate(SL_mm = EstAugSize_fromJulyGrowth) %>%
  select(-EstAugSize_fromJulyGrowth) 

## Actual August Size
AugActualSize <- Growth_Size %>%
  filter(Month ==8) %>%
  mutate(Group = "ActualAugustSize") %>%
  select(FISHID, Year, Heatwave, Group, SL_mm) 
  
## Actual July Size
JulyActualSize <- Growth_Size %>%
  filter(Month ==7) %>%
  mutate(Group = "ActualJulySize") %>%
  select(FISHID, Year, Heatwave, Group, SL_mm) 

Predicted_and_Actual_Size_All <- rbind(AugActualSize, JulyActualSize, AugustPredictedSize)

Predicted_and_Actual_Size_NoUnpairedYears <- Predicted_and_Actual_Size_All %>%
  filter(!Year %in% 2006) %>%
  filter(!Year %in% 2008) %>%
  filter(!Year %in% 2016) 


#Percent Increase in Mean Size between Actual August Sizes and Predicted August Sizes
PredictedActualSize_PercentIncrease_Year <- Predicted_and_Actual_Size_NoUnpairedYears %>%
  filter(Group != "ActualJulySize") %>%
  group_by(Year, Group) %>%
  reframe(SL_mean = mean(SL_mm, na.rm = T)) %>%
  pivot_wider(names_from = Group, values_from = c(SL_mean)) %>%
  mutate(MeanPercentIncrease_PredictedAugSL_ActualAugSL = ((ActualAugustSize - PredictedAugustSize)/PredictedAugustSize) * 100)  %>%
  mutate(Heatwave = case_when(
    Year == '2007'|Year == '2009' | Year == '2010' | Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) %>%
  ungroup() %>%
  group_by(Year) %>%
  reframe(Mean_byYear = mean(MeanPercentIncrease_PredictedAugSL_ActualAugSL))

PredictedActualSize_PercentIncrease_HW <- Predicted_and_Actual_Size_NoUnpairedYears %>%
  filter(Group != "ActualJulySize") %>%
  group_by(Year, Group) %>%
  reframe(SL_mean = mean(SL_mm, na.rm = T)) %>%
  pivot_wider(names_from = Group, values_from = c(SL_mean)) %>%
  mutate(MeanPercentIncrease_PredictedAugSL_ActualAugSL = ((ActualAugustSize - PredictedAugustSize)/PredictedAugustSize) * 100)  %>%
  mutate(Heatwave = case_when(
    Year == '2007'|Year == '2009' | Year == '2010' | Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) %>%
  ungroup() %>%
  group_by(Heatwave) %>%
  reframe(Mean_byHW = mean(MeanPercentIncrease_PredictedAugSL_ActualAugSL), sd_byHW = sd(MeanPercentIncrease_PredictedAugSL_ActualAugSL), se_byHW = (sd(MeanPercentIncrease_PredictedAugSL_ActualAugSL)/sqrt(length(MeanPercentIncrease_PredictedAugSL_ActualAugSL))))

#Plot Predicted August Size based on July Growth (Fig. 6)
FakeData <-  data.frame(c('BLXXX', 'BLXYY', 'BLXYZ', 'BLWXY'), c(2006, 2008, 2011, 2016), c("Before", "Before", "Before", "Heatwave"), c(1,1,1,1), c("ActualAugustSize","ActualAugustSize", "ActualAugustSize", "ActualAugustSize"),c(1,1,1,1))
colnames(FakeData) <- c("FISHID","Year", "Heatwave", "SL_mm", "Group", "GroupingFactor")
FakeData <- FakeData %>%
  mutate(Year = as.factor(Year)) %>%
  mutate(GroupingFactor = as.character(GroupingFactor)) %>%
  select(FISHID, Year, Heatwave, SL_mm, Group)

Predicted_and_Actual_Size_WithPlaceholders <- rbind(Predicted_and_Actual_Size_NoUnpairedYears, FakeData)%>%
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))) %>% filter(Group != "ActualJulySize")

new_rect_test2 <- data.frame(Year = as.integer(c(2006, 2007, 2008, 2009, 2010,2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)), xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5), Color = c("Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z2", "Z2", "Z2",  "Z3", "Z3", "Z2")) %>% 
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))) %>%
  mutate(Line = "black")

labelb <- data.frame(y = rep(140,14), x = c("2006" , "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"),labels = c("", "-15.7%", "", "-0.7%", "2.2%", "", "7.3%", "9.0%", "39.7%", "34.2%", "", "-8.4%", "-8.0%", "19.8%"), GroupingFactor = c(50,50,50,50,50,50,50,50,50,50, 50, 50, 50, 50)) %>%
  mutate(x= as.factor(x))

PredictedSize <- ggplot() +
  geom_boxplot(data = Predicted_and_Actual_Size_WithPlaceholders, aes(x = Year,  y = SL_mm, fill=Group), position = position_dodge(preserve = "single")) +
  geom_text(data = labelb, aes(x = x, y = y, label = labels), angle = 0, size = 4, fontface = "bold") +
  geom_rect(data = new_rect_test2, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),colour = "white",fill = c("#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea","#b6abea","#b6abea", "#e5b2bb", "#e5b2bb","#e5b2bb", "grey90", "grey90", "#e5b2bb"), alpha = .5, show.legend = F, inherit.aes=FALSE)  +
  geom_boxplot(data = Predicted_and_Actual_Size_WithPlaceholders, aes(x = Year,  y = SL_mm, fill=Group), position = position_dodge(preserve = "single"), color = "grey50") +
  theme_classic() +
  scale_fill_manual(labels = c( "Actual August Size", "Predicted August Size"), values=c(  "grey30","cornsilk"))  +
  coord_cartesian(ylim = c(25, 140)) +
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = (element_text(size = 14)), legend.title=element_blank()) +
  xlab("Year") + 
  ylab("Standard Length (mm)")
print(PredictedSize) 
#dev.copy(jpeg,'Figure_6.jpg', width=11, height=6, units='in', res=300)
#dev.off()

#Figure S6; Actual and Predicted August Size as a Density Plot
RealData_DensityPlot <- Predicted_and_Actual_Size_NoUnpairedYears %>%
  filter(Group != "ActualJulySize") %>%
ggplot() +
  geom_density( aes(x = SL_mm, fill=Group), alpha = 0.5) +
  facet_wrap(~Heatwave, nrow = 3) +
  theme_classic() +
  scale_fill_manual(labels = c( "Actual August Size", "Predicted August Size"), values=c(  "grey30","cornsilk"))  +
  theme(axis.text.x=element_text(size=10, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 10),  axis.title=element_text(size=10,face="bold"), legend.text = (element_text(size = 8)), legend.title=element_blank()) +
  xlab("Standard Length (mm)") +
  ggtitle("Real Data") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.text = (element_text(size = 12)),  axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), strip.text.x=element_text(size = 14)) +
  theme(legend.position = "bottom") +
  xlim(40, 150)
print(RealData_DensityPlot) 


Predicted_and_Actual_Size_NoUnpairedYears_NoJuly <- Predicted_and_Actual_Size_NoUnpairedYears %>%
  filter(Group != "ActualJulySize")

Summary_RealData <- Predicted_and_Actual_Size_NoUnpairedYears_NoJuly %>%
  select(-FISHID, -Year) %>%
  group_by(Heatwave, Group) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)
##92.66 mm mean SL for HWs

##Linear Mixed Effects Models comparing predicted August size to actual August size

P1 <- lm(SL_mm ~  Heatwave * Group, Predicted_and_Actual_Size_NoUnpairedYears_NoJuly, na.action = na.exclude) 
summary(P1)
anova(P1)

P2 <- lmer(SL_mm ~  Heatwave * Group + (1|Year), na.action = na.exclude, Predicted_and_Actual_Size_NoUnpairedYears_NoJuly, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(P2) 
Anova(P2, type = 3)

AIC(P1, P2)
r.squaredGLMM(P2)

#### Predicted Sizes if only the largest fish in the upper quantiles survive ####

JulyMHW<- Growth_Size %>%
  filter(Month ==7) %>%
  filter(Heatwave == "Heatwave") %>%
  drop_na(mm_mm_day) 

Quantiles_JulyMHW <- JulyMHW %>%
  reframe(quantile(SL_mm, na.rm = T,probs = c(0, .25, .5, .75, .80, .85, .9, .99, 1)))

q0 <- as.numeric(Quantiles_JulyMHW[1,1])
q25<- as.numeric(Quantiles_JulyMHW[2,1])
q50<- as.numeric(Quantiles_JulyMHW[3,1])
q75 <- as.numeric(Quantiles_JulyMHW[4,1])
q80 <- as.numeric(Quantiles_JulyMHW[5,1])
q85 <- as.numeric(Quantiles_JulyMHW[6,1])
q90 <- as.numeric(Quantiles_JulyMHW[7,1])
q99 <- as.numeric(Quantiles_JulyMHW[8,1])
q100 <- as.numeric(Quantiles_JulyMHW[9,1])


SizesBetween_q75_q100 <- JulyMHW %>%
  filter(SL_mm > q75 & SL_mm < q100) %>%
  mutate(Size_Diff = mm_day * 40.3) %>%
  mutate(PredictedAugSize = SL_mm + Size_Diff) %>%
  select(PredictedAugSize) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

SizesBetween_q80_q100 <- JulyMHW %>%
  filter(SL_mm > q80 & SL_mm < q100) %>%
  mutate(Size_Diff = mm_day * 40.3) %>%
  mutate(PredictedAugSize = SL_mm + Size_Diff) %>%
  select(PredictedAugSize) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

SizesBetween_q85_q100 <- JulyMHW %>%
  filter(SL_mm > q85 & SL_mm < q100) %>%
  mutate(Size_Diff = mm_day * 40.3) %>%
  mutate(PredictedAugSize = SL_mm + Size_Diff) %>%
  select(PredictedAugSize) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)
#If only the upper 15% of fish grow out to August, we get the actual August size distribution

SizesBetween_q90_q100 <- JulyMHW %>%
  filter(SL_mm > q90 & SL_mm < q100) %>%
  mutate(Size_Diff = mm_day * 40.3) %>%
  mutate(PredictedAugSize = SL_mm + Size_Diff) %>%
  select(PredictedAugSize) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

#Overshoots a bit, but still quite reasonable


####Simulated populations to attempt to evaluate uncertainty in August Predicted Size ####

#Relationship between Relative Growth and Size for July, all MHW classes

July_Growth_Size_Relationship <- Growth_Size %>%
  filter(Month == "7") %>%
  ggplot(aes(x=SL_mm, y=mm_mm_day, color=Heatwave)) + 
  geom_point(size=3) +
  theme_classic() +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  scale_color_manual(values=c("#402aaa", "grey1", "#c7576a")) +
  scale_fill_manual(values=c("#b6abea", "#888888", "#e5b2bb")) 

Heatwave <- Growth_Size %>%
  filter(Month == "7") %>%
  filter(Heatwave == "Heatwave") 
lm_Heatwave <- lm(mm_mm_day ~ SL_mm, data = Heatwave)

Between <- Growth_Size %>%
  filter(Month == "7") %>%
  filter(Heatwave == "Between") 
lm_Between <- lm(mm_mm_day ~ SL_mm, data = Between)

Before <- Growth_Size %>%
  filter(Month == "7") %>%
  filter(Heatwave == "Before") 
lm_Before <- lm(mm_mm_day ~ SL_mm, data = Before)


# Simulated population of 10,000 During MHWs
#55.08 mean size in July in HWs
#9.49 sd
#40.3 days between sampling during MHWs

set.seed(1234)
SL_mm = rnorm(mean = 55.08, sd = 9.49, n = 10000)
hist(SL_mm)
fish.id <- c(1:10000)

July.During <- as.data.frame(cbind(fish.id, SL_mm)) 
July.During.Growth <- as.data.frame(predict(lm_Heatwave, July.During, interval = "confidence")) %>%
    mutate(mm_mm_day = fit) %>%
  select(mm_mm_day) 

July.During.Population <- cbind(July.During, July.During.Growth) %>%
  mutate(mm_day = mm_mm_day * SL_mm) %>%
  mutate(mm_betweenJulyAug = mm_day * 40.3) %>%
  mutate(EstAugSize_fromJulyGrowth = SL_mm + mm_betweenJulyAug) %>%
  mutate(Heatwave = "Heatwave") %>%
  mutate(Group = "PredictedAugustSize") %>%
  select(-SL_mm) %>%
  mutate(SL_mm = EstAugSize_fromJulyGrowth) %>%
  select(fish.id, Heatwave, Group, SL_mm)

# Simulated population of 10,000 Before MHWs
#43.38 mean size 
#6.33 sd
#37.2 days between sampling

set.seed(1234)
SL_mm = rnorm(mean = 43.38, sd = 6.325, n = 10000)
hist(SL_mm)
fish.id <- c(1:10000)

July.Before<- as.data.frame(cbind(fish.id, SL_mm)) 
July.Before.Growth <- as.data.frame(predict(lm_Before, July.Before, interval = "confidence")) %>%
  mutate(mm_mm_day = fit) %>%
  select(mm_mm_day) 

July.Before.Population <- cbind(July.Before, July.Before.Growth) %>%
  mutate(mm_day = mm_mm_day * SL_mm) %>%
  mutate(mm_betweenJulyAug = mm_day * 37.2) %>%
  mutate(EstAugSize_fromJulyGrowth = SL_mm + mm_betweenJulyAug) %>%
  mutate(Heatwave = "Before") %>%
  mutate(Group = "PredictedAugustSize") %>%
  select(-SL_mm) %>%
  mutate(SL_mm = EstAugSize_fromJulyGrowth) %>%
  select(fish.id, Heatwave, Group, SL_mm)

# Simulated population of 10,000 Between MHWs
#62.67 mean size
#12.937 sd
#36.0 days between sampling during MHWs

set.seed(1234)
SL_mm = rnorm(mean = 62.67, sd = 12.937, n = 10000)
fish.id <- c(1:10000)
hist(SL_mm)

July.Between<- as.data.frame(cbind(fish.id, SL_mm)) 
July.Between.Growth <- as.data.frame(predict(lm_Between, July.Between, interval = "confidence")) %>%
  mutate(mm_mm_day = fit) %>%
  select(mm_mm_day) 

July.Between.Population <- cbind(July.Between, July.Between.Growth) %>%
  mutate(mm_day = mm_mm_day * SL_mm) %>%
  mutate(mm_betweenJulyAug = mm_day * 36.0) %>%
  mutate(EstAugSize_fromJulyGrowth = SL_mm + mm_betweenJulyAug) %>%
  mutate(Heatwave = "Between") %>%
  mutate(Group = "PredictedAugustSize") %>%
  select(-SL_mm) %>%
  mutate(SL_mm = EstAugSize_fromJulyGrowth) %>%
  select(fish.id, Heatwave, Group, SL_mm)


## Predicted August Size from Simulation

PredictedAugSize_Simulation <- rbind(July.During.Population, July.Before.Population, July.Between.Population)

### Actual August Size from Simulation
set.seed(1234)
Aug.Size.During = rnorm(mean = 92.66, sd = 11.014, n = 10000)
fish.id <- c(1:10000)

Aug.During <- as.data.frame(cbind(fish.id, Aug.Size.During)) %>%
  mutate(Heatwave = "Heatwave") %>%
  mutate(SL_mm = Aug.Size.During) %>%
  mutate(Group = "ActualAugustSize") %>%
  select(fish.id, Heatwave, Group, SL_mm)

set.seed(1234)
Aug.Size.Before = rnorm(mean = 67.12, sd = 10.618, n = 10000)
fish.id <- c(1:10000)
Aug.Before <- as.data.frame(cbind(fish.id, Aug.Size.Before)) %>%
  mutate(Heatwave = "Before") %>%
  mutate(SL_mm = Aug.Size.Before) %>%
  mutate(Group = "ActualAugustSize") %>%
  select(fish.id, Heatwave, Group, SL_mm)

set.seed(1234)
Aug.Size.Between = rnorm(mean = 83.95, sd = 17.885, n = 10000)
fish.id <- c(1:10000)
Aug.Between <- as.data.frame(cbind(fish.id, Aug.Size.Between)) %>%
  mutate(Heatwave = "Between") %>%
  mutate(SL_mm = Aug.Size.Between) %>%
  mutate(Group = "ActualAugustSize") %>%
  select(fish.id, Heatwave, Group, SL_mm)

ActualAugSize_Simulation <- rbind(Aug.During, Aug.Before, Aug.Between)

Simulation <- rbind(ActualAugSize_Simulation, PredictedAugSize_Simulation)

#Figure S6b
SimulationPlot_FullSimulatedPopulation <- ggplot() +
  geom_density(data = Simulation, aes(x = SL_mm, fill=Group), alpha = 0.5) +
  facet_wrap(~Heatwave, nrow = 3) +
  theme_classic() +
  scale_fill_manual(labels = c( "Simulation of Actual August Size", "Simulation of Predicted August Size"), values=c(  "grey30","cornsilk"))  +
  theme(axis.text.x=element_text(size=10, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 10),  axis.title=element_text(size=10,face="bold"), legend.text = (element_text(size = 8)), legend.title=element_blank()) +
  xlab("Standard Length (mm)") + 
  ggtitle("Simulated Population (n = 10,000)") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.text = (element_text(size = 12)),  axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), strip.text.x=element_text(size = 14)) +
  theme(legend.position = "bottom") +
  xlim(40, 150)
print(SimulationPlot_FullSimulatedPopulation) 

#### Now, pull the same size sample as what we pulled for our fish samples
#Use the Growth n for July and August (Table 1)
#
set.seed(1234)
rand_PredictedAugust_During <- July.During.Population[sample(nrow(July.During.Population), size=69), ]
set.seed(1234)
rand_PredictedAugust_Before <- July.Before.Population[sample(nrow(July.Before.Population), size=137), ]
set.seed(1234)
rand_PredictedAugust_Between <- July.Between.Population[sample(nrow(July.Between.Population), size=44), ]

set.seed(1234)
rand_Aug.During <- Aug.During[sample(nrow(Aug.During), size = 50), ]
set.seed(1234)
rand_Aug.Before <- Aug.Before[sample(nrow(Aug.Before), size = 144), ]
set.seed(1234)
rand_Aug.Between <- Aug.Between[sample(nrow(Aug.Between), size = 44), ]

Subset_SimulatedPopulation <- rbind(rand_PredictedAugust_During, rand_PredictedAugust_Before, rand_PredictedAugust_Between,  rand_Aug.During,  rand_Aug.Before, rand_Aug.Between  )

#Figure S6c
SimulationPlot_SubsetPopulation <- ggplot() +
  geom_density(data = Subset_SimulatedPopulation, aes(x = SL_mm, fill=Group), alpha = 0.5) +
  facet_wrap(~Heatwave, nrow = 3) +
  theme_classic() +
  scale_fill_manual(labels = c( "Simulation of Actual August Size", "Simulation of Predicted August Size"), values=c(  "grey30","cornsilk"))  +
  theme(axis.text.x=element_text(size=10, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 10),  axis.title=element_text(size=10,face="bold"), legend.text = (element_text(size = 8)), legend.title=element_blank()) +
  xlab("Standard Length (mm)") +
  ggtitle("Random Draws from Simulated Population") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.text = (element_text(size = 12)),  axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), strip.text.x=element_text(size = 14)) +
  theme(legend.position = "bottom") +
  xlim(40, 150)
print(SimulationPlot_SubsetPopulation) 


#Calculate means and confidence intervals for ALL simulated data

SimulatedSizeSummary <- Simulation %>%
  select(-fish.id) %>% 
  group_by(Group, Heatwave) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

#Statistics
S1 <- lm(SL_mm ~ Heatwave * Group, Simulation) 
summary(S1)
Anova(S1, type =3 )

#Calculate means and confidence intervals for SUBSET simulated data

Subset_SimulatedSizeSummary <- Subset_SimulatedPopulation %>%
  select(-fish.id) %>% 
  group_by(Group, Heatwave) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)


#Statistics
S2 <- lm(SL_mm ~ Heatwave * Group, Subset_SimulatedPopulation) 
summary(S2)
Anova(S2, type =3 )


#### Plot All Data Distributions: Figure S6

SizeDistributions <- RealData_DensityPlot + SimulationPlot_FullSimulatedPopulation + SimulationPlot_SubsetPopulation +  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = list(c('a', 'b', 'c'))) &
  theme(plot.tag = element_text(face = 'bold', size = 16)) 
#dev.copy(jpeg,'Figure_S6.jpg', width=20, height=8, units='in', res=300)
#dev.off()

  
### Simulate what size quantiles DURING MHWS you need to get August Sizes ####
#92.7 mm mean size for August during HWs 

JulyMHW_FullSimulation<- cbind(July.During, July.During.Growth) %>%
  mutate(mm_day = mm_mm_day * SL_mm) %>%
  mutate(mm_betweenJulyAug = mm_day * 40.3) %>%
  mutate(EstAugSize_fromJulyGrowth = SL_mm + mm_betweenJulyAug) 

Quantiles_JulyMHW_FullSimulation <- JulyMHW_FullSimulation %>%
  reframe(quantile(SL_mm, na.rm = T,probs = c(0, .25, .5, .70, .75, .80, .85, .9, .99, 1)))

q0 <- as.numeric(Quantiles_JulyMHW_FullSimulation[1,1])
q25<- as.numeric(Quantiles_JulyMHW_FullSimulation[2,1])
q50<- as.numeric(Quantiles_JulyMHW_FullSimulation[3,1])
q70 <- as.numeric(Quantiles_JulyMHW_FullSimulation[4,1])
q75 <- as.numeric(Quantiles_JulyMHW_FullSimulation[5,1])
q80 <- as.numeric(Quantiles_JulyMHW_FullSimulation[6,1])
q85 <- as.numeric(Quantiles_JulyMHW_FullSimulation[7,1])
q90 <- as.numeric(Quantiles_JulyMHW_FullSimulation[8,1])
q99 <- as.numeric(Quantiles_JulyMHW_FullSimulation[9,1])
q100 <- as.numeric(Quantiles_JulyMHW_FullSimulation[10,1])

FullSimSizesBetween_q70_q100 <- JulyMHW_FullSimulation %>%
  filter(SL_mm > q70 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)
#Top 30% of fish to get to Actual August size

FullSimSizesBetween_q75_q100 <- JulyMHW_FullSimulation %>%
  filter(SL_mm > q75 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)
#25% overshoots a little, but still looks pretty good

FullSimSizesBetween_q80_q100 <- JulyMHW_FullSimulation %>%
  filter(SL_mm > q80 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

FullSimSizesBetween_q85_q100 <- JulyMHW_FullSimulation %>%
  filter(SL_mm > q85 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)


#### Quantiles from randomly drawn subset During MHWs
set.seed(1234)
rand_JulyMHW <- JulyMHW_FullSimulation[sample(nrow(JulyMHW_FullSimulation), size=69), ]

Quantiles_JulyMHW_Subset <- rand_JulyMHW  %>%
  reframe(quantile(SL_mm, na.rm = T,probs = c(0, .25, .5, .7,  .75, .80, .85, .9, .95, .99, 1)))

q0 <- as.numeric(Quantiles_JulyMHW_Subset [1,1])
q25<- as.numeric(Quantiles_JulyMHW_Subset [2,1])
q50<- as.numeric(Quantiles_JulyMHW_Subset [3,1])
q70<- as.numeric(Quantiles_JulyMHW_Subset [4,1])
q75 <- as.numeric(Quantiles_JulyMHW_Subset [5,1])
q80 <- as.numeric(Quantiles_JulyMHW_Subset [6,1])
q85 <- as.numeric(Quantiles_JulyMHW_Subset [7,1])
q90 <- as.numeric(Quantiles_JulyMHW_Subset [8,1])
q95 <- as.numeric(Quantiles_JulyMHW_Subset [9,1])
q99 <- as.numeric(Quantiles_JulyMHW_Subset [10,1])
q100 <- as.numeric(Quantiles_JulyMHW_Subset [11,1])

SubsetSimSizesBetween_q70_q100 <- rand_JulyMHW  %>%
  filter(SL_mm > q70 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)


SubsetSimSizesBetween_q75_q100 <- rand_JulyMHW  %>%
  filter(SL_mm > q75 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)
#Top 25% is right on 


SubsetSimSizesBetween_q80_q100 <- rand_JulyMHW %>%
  filter(SL_mm > q80 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

FullSimSizesBetween_q85_q100 <- rand_JulyMHW  %>%
  filter(SL_mm > q85 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

FullSimSizesBetween_q90_q100 <- rand_JulyMHW  %>%
  filter(SL_mm > q90 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

FullSimSizesBetween_q95_q100 <- rand_JulyMHW  %>%
  filter(SL_mm > q95 & SL_mm < q100) %>%
  select(EstAugSize_fromJulyGrowth) %>%
  summarise_all(list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T), cv = ~sd(., na.rm = T)/mean(., na.rm =T), se = ~sd(., na.rm = T)/sqrt(length(.)), CI = ~1.96*sd(.,na.rm = T)/sqrt(length(.)), length = ~length(.)),na.rm = TRUE)

