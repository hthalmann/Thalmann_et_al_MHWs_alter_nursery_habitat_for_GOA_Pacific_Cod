#### Kodiak Island Juvenile Pacific Cod Abundance, Size, and Condition 
### Date Created: 1/11/2024
# R v. 4.2.2

##Load libraries 
library(rstudioapi) #Load Working Directory
library(tidyverse) #For data wrangling and visualization
library(lubridate) #For working with dates
library(patchwork) #Data visualization
library(car) #Run ANOVA function for linear models
library(lme4) #Runs linear mixed effects models
library(MuMIn) #Runs the R.squaredGLMM function
library(viridis) #Color Palette


#Paired with:
#Seine data 
#Site data
#Raw Size data, which includes fish standard length, weight, and liver weight

##Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#Sets directory to folder that this R script is saved to

#### Calculate Catch Per Unit Effort ####
Site <- read.csv("D1.SiteData.csv") %>%
  select(SITEID, Bay, Habitat)

CPUE <- read.csv("D2.SeineData.csv") %>%
  mutate(Date = ymd(Date)) %>%
  mutate(JulianDate = yday(Date)) %>%
  filter(!Year %in% 2020) %>%
  filter(!Year %in% 2021) %>%
  mutate(Heatwave = case_when(
    Year == 2006 | Year == 2007 | Year == 2008 | Year == 2009 | Year == 2010 | Year == 2011 | Year == 2012| Year == 2013 ~ "Before",
    Year == 2014 | Year == 2015 | Year == 2016 | Year == 2019~ "During",
    Year == 2017 | Year == 2018 ~ "Between"
  )) %>%
  select(SeineID, Month, Year, Heatwave,  SITEID, Date, JulianDate, TotalGadidCount, No.PacificCod) %>%
  left_join(Site, by = "SITEID") %>%
  mutate(NoPacificCod = as.numeric(No.PacificCod))%>%
  mutate(TotalGadidCount = as.numeric(TotalGadidCount))%>%
  filter(!Bay %in% "Offshore") %>%
  group_by(Year, Month) %>%
  mutate(CPUE_PCod = mean(NoPacificCod, na.rm = T)) %>%
  mutate(CPUE_TotalGadid = mean(TotalGadidCount, na.rm = T)) %>%
  mutate(AvgJulianDate = mean(JulianDate)) %>%
  filter(Month == 7 | Month == 8)


#Summary of CPUE (Table 1, Table S1)

CPUE_summary <- CPUE %>%
  group_by(Month,Year) %>%
  reframe(No.PacificCod = mean(NoPacificCod, na.rm = T), CPUE_se = (sd(NoPacificCod, na.rm = T)/sqrt(length(NoPacificCod))), AverageDay = mean(JulianDate, na.rm = T))

CPUE_summary_HW <- CPUE %>%
  group_by( Heatwave, Month) %>%
  reframe(No.PacificCod = mean(NoPacificCod, na.rm = T), CPUE_se = (sd(NoPacificCod, na.rm = T)/sqrt(length(NoPacificCod))))


CPUE_2016<- data.frame(2016, 8,NA, NA, NA)
colnames(CPUE_2016) <- c("Year","Month", "No.PacificCod", "CPUE_se", "AverageDay")
CPUE_2016

CPUE_all <- rbind(CPUE_summary, CPUE_2016) %>%
  mutate(Heatwave = case_when(
    Year == '2006'| Year == '2007'|Year == '2008'|Year == '2009' | Year == '2010' | Year == '2011'|Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year == '2016' | Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) 


## Plot Nursery CPUE by month (Figure 1a)

NurseryCPUE<- ggplot(CPUE_all, aes(x=Year, y=No.PacificCod, color = as.factor(Month), na.rm = T)) + 
  geom_errorbar(aes(ymin=No.PacificCod-CPUE_se, ymax=No.PacificCod+CPUE_se), width=.1, na.rm = T) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(labels = c( "July","August"), values=c( "grey60", "grey30"))  +
  geom_point() +
  scale_x_continuous("Year", labels = as.character(CPUE_all$Year), breaks = CPUE_all$Year) +
  ylab("CPUE") + 
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = element_text(size = 14),  legend.title=element_blank(), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
print(NurseryCPUE)



#### Size, Body Mass, and Condition  ####
Morphology <- read.csv("D3.Size_Condition_Raw.csv") %>%
  select(FISHID, Year, Month, SL_mm, WholeBodyWW_g, TOTAL_Liver_WW_g) %>%
  mutate(Heatwave = case_when(
    Year == '2006'| Year == '2007'|Year == '2008'|Year == '2009' | Year == '2010' | Year == '2011'| Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year =='2016'| Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) %>%
  filter(!Month %in% 9) %>%
  mutate(logSL = log(SL_mm)) %>%
  mutate(logWeight = log(WholeBodyWW_g)) %>%
  mutate(Year = as.factor(Year)) %>%
  mutate(Month = as.factor(Month)) %>%
  mutate(HSI = (TOTAL_Liver_WW_g/WholeBodyWW_g)*100) %>%
  mutate(HSI_numeric = as.numeric(HSI))

#Summary of Size and Body Mass (Table 1, Table S2)
Size_Mass_Summary_byYear <- Morphology %>%
  group_by(Month,Year) %>%
  reframe(SampleSize = length(SL_mm), SL_mean = mean(SL_mm), SL_se = (sd(SL_mm, na.rm = T)/sqrt(length(SL_mm))), SL_min = min(SL_mm), SL_max = max(SL_mm), BodyMass_mean = mean(WholeBodyWW_g), BodyMass_se = (sd(WholeBodyWW_g, na.rm = T)/sqrt(length(WholeBodyWW_g))), BodyMass_min = min(WholeBodyWW_g), BodyMass_max = max(WholeBodyWW_g))


Size_Mass_Summary_byHW <- Morphology %>%
  group_by(Month,Heatwave) %>%
  reframe(SampleSize = length(SL_mm), SL_mean = mean(SL_mm), SL_se = (sd(SL_mm, na.rm = T)/sqrt(length(SL_mm))), SL_min = min(SL_mm), SL_max = max(SL_mm), BodyMass_mean = mean(WholeBodyWW_g), BodyMass_se = (sd(WholeBodyWW_g, na.rm = T)/sqrt(length(WholeBodyWW_g))), BodyMass_min = min(WholeBodyWW_g), BodyMass_max = max(WholeBodyWW_g))

# Intra-annual Percent Increase in Standard Length between July and August 

StandardLength_PercentIncrease <- Morphology %>%
  mutate(MonthName = case_when(Month == 7 ~ "July", Month == 8 ~ "August")) %>%
  group_by(Year, MonthName) %>%
  reframe(SL_mean = mean(SL_mm, na.rm = T)) %>%
  pivot_wider(names_from = MonthName, values_from = c(SL_mean)) %>%
  mutate(MeanPercentIncrease_JulySL_AugustSL = ((August - July)/July) * 100)  %>%
  mutate(Heatwave = case_when(
    Year == '2007'|Year == '2009' | Year == '2010' | Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) %>%
  ungroup() %>%
  group_by(Heatwave) %>%
  reframe(Mean_byHW = mean(MeanPercentIncrease_JulySL_AugustSL), sd_byHW = sd(MeanPercentIncrease_JulySL_AugustSL), se_byHW = (sd(MeanPercentIncrease_JulySL_AugustSL, na.rm = T)/sqrt(length(MeanPercentIncrease_JulySL_AugustSL))))

# Intra-annual Percent Increase in Body Mass between July and August 

BodyMass_PercentIncrease <- Morphology %>%
  mutate(MonthName = case_when(Month == 7 ~ "July", Month == 8 ~ "August")) %>%
  group_by(Year, MonthName) %>%
  reframe(Mass_mean = mean(WholeBodyWW_g, na.rm = T)) %>%
  pivot_wider(names_from = MonthName, values_from = c(Mass_mean)) %>%
  mutate(MeanPercentIncrease_JulyMass_AugustMass = ((August - July)/July) * 100)  %>%
  mutate(Heatwave = case_when(
    Year == '2007'|Year == '2009' | Year == '2010' | Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) %>%
  ungroup() %>%
  group_by(Heatwave) %>%
  reframe(Mean_byHW = mean(MeanPercentIncrease_JulyMass_AugustMass), sd_byHW = sd(MeanPercentIncrease_JulyMass_AugustMass), se_byHW = (sd(MeanPercentIncrease_JulyMass_AugustMass, na.rm = T)/sqrt(length(MeanPercentIncrease_JulyMass_AugustMass))))

# Intra-annual Percent Increase in HSI between July and August 

HSI_PercentIncrease <- Morphology %>%
  mutate(MonthName = case_when(Month == 7 ~ "July", Month == 8 ~ "August")) %>%
  group_by(Year, MonthName) %>%
  reframe(HSI_mean = mean(HSI_numeric, na.rm = T)) %>%
  pivot_wider(names_from = MonthName, values_from = c(HSI_mean)) %>%
  mutate(MeanPercentIncrease_JulyHSI_AugustHSI = ((August - July)/July) * 100)  %>%
  mutate(Heatwave = case_when(
    Year == '2007'|Year == '2009' | Year == '2010' | Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) %>%
  ungroup() %>%
  group_by(Heatwave) %>%
  reframe(Mean_byHW = mean(MeanPercentIncrease_JulyHSI_AugustHSI), sd_byHW = sd(MeanPercentIncrease_JulyHSI_AugustHSI), se_byHW = (sd(MeanPercentIncrease_JulyHSI_AugustHSI, na.rm = T)/sqrt(length(MeanPercentIncrease_JulyHSI_AugustHSI))))

## Plot log-transformed Size and Weight Relationship (Figure S8) 
LengthWeightRelationship <- ggplot(Morphology, aes(y=logWeight, x = logSL)) +
  geom_point(aes(color = Year), alpha = 1) +
  scale_color_viridis_d(option="turbo") +
  geom_smooth(method = lm, level = 0.95, color = "grey60", linewidth = 1.5)+
  theme_classic(base_size = 17, base_family = "")+
  theme(legend.position = "right")+
  labs(y = "log (Whole Body Wet Weight (g))", x = "log (Standard Length (mm))") +
  theme(axis.text.x=element_text(size=15, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = (element_text(size = 14)), legend.title = (element_text(size = 15)))+
  guides(fill=guide_legend(title="")) 
print(LengthWeightRelationship)
#dev.copy(jpeg,'Figure_S7.jpg', width=5, height=4, units='in', res=300)
#dev.off()


## Calculate Length-Weight Condition Residuals (combine months for equation)
lm1.Weight.Size <- lm(logWeight ~ logSL, data = Morphology)
summary(lm1.Weight.Size)
resid1 <- resid(lm1.Weight.Size)
plot(resid1)


## Add LW Condition Residuals to Morphology File
Morphology$LWResiduals_All <- resid1
summary(Morphology)

#Create Size and Condition R Output File for Subsequent Analyses
#write.csv(Morphology, "Size_Condition_ROutput.csv")


##### Plot Standard Length (Figure 1b) ####

#Add a placeholder for 2011
FakeData <-  data.frame('BLXXX', 2011, "Before", 7, 1)
colnames(FakeData) <- c("FISHID","Year", "Heatwave", "Month", "SL_mm")
FakeData <- FakeData %>%
  mutate(Year = as.factor(Year))

StandardLength <- Morphology %>%
  select(FISHID, Year, Month, Heatwave, SL_mm)

ActualSize_EmptyYears <- rbind(StandardLength,FakeData) %>%
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")))

new_rect_1 <- data.frame(Year = as.integer(c(2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)), xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5), Color = c("Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z2", "Z2", "Z2", "Z3", "Z3", "Z2")) %>% 
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))) %>%
  mutate(Line = "black")

labela <- data.frame(y = rep(140,14), x = c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"),labels = c("", "22.8%", "", "68.9%", "56.1%", "", "60.0%", "54.4%", "103.0%", "96.9%", "", "33.8%", "34.9%", "83.4%"), Month = c(50,50,50,50,50,50,50,50,50,50, 50, 50, 50, 50)) %>%
  mutate(x= as.factor(x))

StandardLength_Plot <- ggplot() +
  geom_boxplot(data = ActualSize_EmptyYears, aes(x = Year,  y = SL_mm, fill= Month), position = position_dodge(preserve = "single")) +
  geom_text(data = labela, aes(x = x, y = y, label = labels), angle = 0, size = 4, fontface = "bold") +
  geom_rect(data = new_rect_1, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),colour = "white",fill = c("#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#e5b2bb", "#e5b2bb","#e5b2bb", "grey90", "grey90", "#e5b2bb"), alpha = .5, show.legend = F, inherit.aes=FALSE)  +
  geom_boxplot(data = ActualSize_EmptyYears, aes(x = Year,  y = SL_mm, fill=Month), position = position_dodge(preserve = "single"), color = "grey50") +
  theme_classic() +
  scale_fill_manual(labels = c( "July","August"), values=c( "grey60", "grey30"))  +
  coord_cartesian(ylim = c(25, 140)) +
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = (element_text(size = 14)), legend.title=element_blank(), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  xlab("Year") + 
  ylab("Standard Length (mm)") 
print(StandardLength_Plot) 


##### Plot Minimum Standard Length for July and August (Figure 1c) ####

MinimumSL_Summary <- Morphology %>%
  group_by(Month, Year, Heatwave) %>%
  reframe(meanSL = mean(SL_mm, na.rm = T), SL_sd = sd(SL_mm, na.rm = T), SL_se = (sd(SL_mm, na.rm = T)/sqrt(length(SL_mm))),  maxSL = max(SL_mm), minSL = min(SL_mm), SL_range = maxSL - minSL) 

MinimumSL_Summary_Aug <- MinimumSL_Summary %>%
  filter(Month == 8) 

MinimumSL_Summary_July <- MinimumSL_Summary %>%
  filter(Month == 7) 

#Simple Linear Regression: July Minimum SL by Heatwave
lm_JulyMinSL <- lm(minSL ~ Heatwave,data = MinimumSL_Summary_July)
summary(lm_JulyMinSL)

#Simple Linear Regression: August Minimum SL by Heatwave
lm_AugMinSL <- lm(minSL ~ Heatwave,data = MinimumSL_Summary_Aug)
summary(lm_AugMinSL)


July_minSL <- ggplot(MinimumSL_Summary_July) +
  geom_boxplot( aes(x = Heatwave, y = minSL, fill = Heatwave)) + 
  ylab("Minimum Standard Length (mm)") + 
  xlab("Heatwave Class") + 
  ggtitle("July") +
  theme_classic() +
  scale_fill_manual(values=c("#b6abea", "grey90", "#e5b2bb")) +
  annotate("text", x = 1, y=45, label = "a", fontface = "bold", size = 5) +
  annotate("text", x = 2, y=45, label = "a", fontface = "bold", size = 5) +
  annotate("text", x = 3, y=54, label = "a", fontface = "bold", size = 5) +
  scale_fill_manual(values=c("#b6abea", "grey90", "#e5b2bb")) +
  coord_cartesian(ylim = c(20, 105)) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_blank(), plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(July_minSL)

Aug_minSL <- ggplot(MinimumSL_Summary_Aug) +
  geom_boxplot( aes(x = Heatwave, y = minSL, fill = Heatwave)) + 
  ylab("Minimum Standard Length (mm)") + 
  xlab("Heatwave Class") + 
  ggtitle("August") +
  theme_classic() +
  annotate("text", x = 1, y=60, label = "b", fontface = "bold", size = 5) +
  annotate("text", x = 2, y=60, label = "b", fontface = "bold", size = 5) +
  annotate("text", x = 3, y=99, label = "c", fontface = "bold", size = 5) +
  scale_fill_manual(values=c("#b6abea", "grey90", "#e5b2bb")) +
  coord_cartesian(ylim = c(20, 105)) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_blank(), plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(Aug_minSL)

MinSL <- July_minSL + Aug_minSL +
  plot_layout(nrow = 1, guides = "collect") + 
  theme(plot.tag = element_text(face = 'bold', size = 16)) 
print(MinSL)

##### Plot HSI Condition (Figure 1d) ####

#Add a placeholder for 2011
FakeData <-  data.frame('BLXXX', 2011, "Before", 7, 22)
colnames(FakeData) <- c("FISHID","Year", "Heatwave", "Month", "HSI_numeric")
FakeData <- FakeData %>%
  mutate(Year = as.factor(Year))

HSI <- Morphology %>%
  select(FISHID, Year, Month, Heatwave, HSI_numeric)

ActualHSI_EmptyYears <- rbind(HSI,FakeData) %>%
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")))

new_rect_1 <- data.frame(Year = as.integer(c(2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)), xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5), Color = c("Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z2", "Z2", "Z2", "Z3", "Z3", "Z2")) %>% 
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))) %>%
  mutate(Line = "black")


HSI_Plot <- ggplot() +
  geom_boxplot(data = ActualHSI_EmptyYears, aes(x = Year,  y = HSI_numeric, fill= Month), position = position_dodge(preserve = "single")) +
  geom_rect(data = new_rect_1, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),colour = "white",fill = c("#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#e5b2bb", "#e5b2bb","#e5b2bb", "grey90", "grey90", "#e5b2bb"), alpha = .5, show.legend = F, inherit.aes=FALSE)  +
  geom_boxplot(data = ActualHSI_EmptyYears, aes(x = Year,  y = HSI_numeric, fill=Month), position = position_dodge(preserve = "single"), color = "grey50") +
  theme_classic() +
  scale_fill_manual(labels = c( "July","August"), values=c( "grey60", "grey30"))  +
  coord_cartesian(ylim = c(0, 10.65)) +
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = (element_text(size = 14)), legend.title=element_blank(), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  xlab("Year") + 
  ylab("Hepatosomatic Index") 
print(HSI_Plot) 


#Combined Abundance and Size and Minimum Size File (Figure 1)
Abundance_Size_MinSize <- NurseryCPUE + StandardLength_Plot  + MinSL  + HSI_Plot +
  plot_layout(nrow = 4) + 
  plot_annotation(tag_levels = list(c('a', 'b', 'c', '', 'd'))) &
  theme(plot.tag = element_text(face = 'bold', size = 16)) 
print(Abundance_Size_MinSize)
#dev.copy(jpeg,'Figure_1.jpg', width=14, height=18, units='in', res=300)
#dev.off()

##### Plot Body Mass (Figure S2a) ####

#Add a placeholder for 2011
FakeData_Mass <-  data.frame('BLXXX', 2011, "Before", 7, 30)
colnames(FakeData_Mass) <- c("FISHID","Year", "Heatwave", "Month", "WholeBodyWW_g")
FakeData_Mass <- FakeData_Mass %>%
  mutate(Year = as.factor(Year))

Mass <- Morphology %>%
  select(FISHID, Year, Month, Heatwave, WholeBodyWW_g)

ActualMass_EmptyYears <- rbind(Mass,FakeData_Mass) %>%
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")))

new_rect_1 <- data.frame(Year = as.integer(c(2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)), xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5), Color = c("Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z2", "Z2", "Z2", "Z3", "Z3", "Z2")) %>% 
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))) %>%
  mutate(Line = "black")

labelc<- data.frame(y = rep(24.3,14), x = c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"),labels = c("", "80.0%", "", "392%", "304%", "", "327%", "339%", "804%", "808%", "", "123%", "170%", "565%"), Month = c(50,50,50,50,50,50,50,50,50,50, 50, 50, 50, 50)) %>%
  mutate(x= as.factor(x))


Mass_Plot <- ggplot() +
  geom_boxplot(data = ActualMass_EmptyYears, aes(x = Year,  y = WholeBodyWW_g, fill= Month), position = position_dodge(preserve = "single")) +
    geom_text(data = labelc, aes(x = x, y = y, label = labels), angle = 0, size = 4, fontface = "bold") +
  geom_rect(data = new_rect_1, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),colour = "white",fill = c("#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#e5b2bb", "#e5b2bb","#e5b2bb", "grey90", "grey90", "#e5b2bb"), alpha = .5, show.legend = F, inherit.aes=FALSE)  +
  geom_boxplot(data = ActualMass_EmptyYears, aes(x = Year,  y = WholeBodyWW_g, fill=Month), position = position_dodge(preserve = "single"), color = "grey50") +
  theme_classic() +
  scale_fill_manual(labels = c( "July","August"), values=c( "grey60", "grey30"))  +
  coord_cartesian(ylim = c(0, 24)) +
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = (element_text(size = 14)), legend.title=element_blank(), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  xlab("Year") + 
  ylab("Body Mass (g)") 
print(Mass_Plot) 

##### Plot LW Resids (Figure S2b) ####

#Add a placeholder for 2011
FakeData_LW <-  data.frame('BLXXX', 2011, "Before", 7, 1.5)
colnames(FakeData_LW) <- c("FISHID","Year", "Heatwave", "Month", "LWResiduals_All")
FakeData_LW <- FakeData_LW %>%
  mutate(Year = as.factor(Year))

LW <- Morphology %>%
  select(FISHID, Year, Month, Heatwave, LWResiduals_All)

ActualLW_EmptyYears <- rbind(LW,FakeData_LW) %>%
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019")))

new_rect_1 <- data.frame(Year = as.integer(c(2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)), xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5), xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5), Color = c("Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z2", "Z2", "Z2", "Z3", "Z3", "Z2")) %>% 
  mutate(Year = factor(Year, levels=c("2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))) %>%
  mutate(Line = "black")


LW_Plot <- ggplot() +
  geom_boxplot(data = ActualLW_EmptyYears, aes(x = Year,  y = LWResiduals_All, fill= Month), position = position_dodge(preserve = "single")) +
  geom_rect(data = new_rect_1, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),colour = "white",fill = c("#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#b6abea", "#e5b2bb", "#e5b2bb","#e5b2bb", "grey90", "grey90", "#e5b2bb"), alpha = .5, show.legend = F, inherit.aes=FALSE)  +
  geom_boxplot(data = ActualLW_EmptyYears, aes(x = Year,  y = LWResiduals_All, fill=Month), position = position_dodge(preserve = "single"), color = "grey50") +
  theme_classic() +
  scale_fill_manual(labels = c( "July","August"), values=c( "grey60", "grey30"))  +
  coord_cartesian(ylim = c(-1, 1)) +
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = (element_text(size = 14)), legend.title=element_blank(), plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  xlab("Year") + 
  ylab("Length-Weight Condition Residuals") 
print(LW_Plot) 


#Combine to create Figure S2
Condition_Mass <- Mass_Plot + LW_Plot+ 
  plot_layout(nrow = 2, guides = "collect") + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16)) 
print(Condition_Mass)
#dev.copy(jpeg,'Figure_S2.jpg', width=12, height=10, units='in', res=300)
#dev.off()


#### Linear Mixed Models for Size and Condition (Table S3) ####

##### Abundance Model #####

#No Random Effects, all Fixed Effects and Interactions
A1 <- lm(log(No.PacificCod) ~  Heatwave * Month , na.action = na.exclude, CPUE_all) 
summary(A1)
anova(A1)


#Year as a random Intercept
A2<- lmer(log(No.PacificCod) ~  Heatwave *Month+ (1|Year), na.action = na.exclude, CPUE_all, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(A2) 
Anova(A2, type =3)

AIC(S1, S2)
r.squaredGLMM(S2)
##### Size Model #####

#No Random Effects, all Fixed Effects and Interactions
S1 <- lm(logSL ~  Heatwave * Month , na.action = na.exclude, Morphology) 
summary(S1)
anova(S1)


#Year as a random Intercept
S2<- lmer(logSL ~  Heatwave *Month+ (1|Year), na.action = na.exclude, Morphology, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(S2) 
Anova(S2, type =3)

AIC(S1, S2)
r.squaredGLMM(S2)


##### Body Mass Model #####

#No Random Effects, all Fixed Effects and Interactions
W1 <- lm(logWeight ~  Heatwave * Month , na.action = na.exclude, Morphology) 
summary(W1)
anova(W1)
#Anova(M1)

#Year as a random Intercept
W2<- lmer(logWeight ~  Heatwave *Month+ (1|Year), na.action = na.exclude, Morphology, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(W2) 
Anova(W2, type = 3)

AIC(W1, W2)
r.squaredGLMM(W2)


##### HSI Model ####
#No Random Effects, all Fixed Effects and Interactions
Morphology_NA <- na.omit(Morphology) #Remove all fish without liver weights
H1 <- lm(HSI_numeric ~  Heatwave * Month, Morphology_NA, na.action = na.exclude) 
summary(H1)
anova(H1)

#Year as a random Intercept
H2<- lmer(HSI_numeric ~  Heatwave *Month+ (1|Year), na.action = na.exclude, Morphology_NA, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(H2) 
Anova(H2, type = 3)

AIC(H1, H2)
r.squaredGLMM(H2)

#Pull out 2019 and rerun

Morphology_NA_No2019 <- Morphology_NA %>%
  filter(Year != '2019')

H3 <- lm(HSI_numeric ~  Heatwave * Month, Morphology_NA_No2019, na.action = na.exclude) 
summary(H3)
anova(H3)

H4<- lmer(HSI_numeric ~  Heatwave * Month+ (1|Year), na.action = na.exclude, Morphology_NA_No2019, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(H4) 
Anova(H4, type = 3)

AIC(H3, H4)
r.squaredGLMM(H4)


##### LW Residual Condition #####
#No Random Effects, all Fixed Effects and Interactions

L1 <- lm(LWResiduals_All ~  Heatwave * Month , na.action = na.exclude, Morphology) 
summary(L1)
anova(L1)

#Year as a random Intercept
L2<- lmer(LWResiduals_All ~  Heatwave *Month+ (1|Year), na.action = na.exclude, Morphology, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(L2) 
Anova(L2, type = 3)

AIC(L1, L2)
r.squaredGLMM(L2)


