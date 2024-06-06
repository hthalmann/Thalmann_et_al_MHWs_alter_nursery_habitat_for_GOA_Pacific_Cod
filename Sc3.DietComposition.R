#### Kodiak Island Juvenile Pacific Cod Stomach Fullness and Diet Composition
### Date Created: 1/11/2024
# R v. 4.2.2

## Load libraries 
library(rstudioapi) #Load Working Directory
library(tidyverse) #For data wrangling and visualization
library(lubridate) #For working with dates
library(patchwork) #Data visualization
library(lme4) #Runs linear models
library(car)  #Runs ANOVAs for linear models
library(MuMIn) #Runs the R.squaredGLMM function
library(PNWColors) #Color Palette


#Paired with:
#Raw stomach fullness data
#Size and Condition data
#Raw stomach contents data with counts and weights of prey species for each individual fish


##Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#Sets directory to folder that this R script is saved to


#### Calculate Stomach Fullness ####
StomachFullness <- read.csv('D5.StomachFullness_Raw.csv') 
SizeCondition <- read.csv("D6.Size_Condition_ROutput.csv") %>%
  select(FISHID, SL_mm, WholeBodyWW_g, Heatwave) 

Fullness <- StomachFullness %>%
  left_join(SizeCondition) %>%
  mutate(Month = as.factor(Month)) %>%
  mutate(Year = as.factor(Year)) %>%
  mutate(StomachFullness = ContentsWW_g/(WholeBodyWW_g - ContentsWW_g))%>%
  mutate(sqrtFullness = sqrt(StomachFullness)) #Square-root transform Stomach Fullness

#Create Stomach Fullness R Output File for Subsequent Analyses
#write.csv(Fullness, "StomachFullness_ROutput.csv")

#Determine number of empty stomachs
Empty <- Fullness %>%
  filter(StomachFullness == 0) %>%
  summarise(length(StomachFullness))
#12 empty stomachs


##### Plot Stomach Fullness (Figure S3) ####

new_rect <- data.frame(Heatwave =c("Before", "Between", "Heatwave"), xmin = c(0.5, 1.5, 2.5), xmax = c(1.5, 2.5, 3.5), Color = c("Z1", "Z2", "Z3"))

StomachFullness_SqrtTransform <- ggplot()+ 
  geom_boxplot(data = Fullness, aes(x = Heatwave, y = sqrtFullness,  fill=Month), position = position_dodge(preserve = "single")) +
  geom_rect(data = new_rect, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),colour = "white",fill = c("#b6abea", "grey90", "#e5b2bb"), alpha = .5, show.legend = F, inherit.aes=FALSE)  +
  geom_boxplot(data = Fullness, aes(x = Heatwave, y = sqrtFullness,  fill=Month), position = position_dodge(preserve = "single")) +
  theme_classic() +
  scale_fill_manual(labels = c( "July","August"), values=c( "grey60", "grey30"))  +
  ylab("Stomach Fullness") +
  xlab("Heatwave Class") +
  theme( axis.text = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = (element_text(size = 12))) 
print(StomachFullness_SqrtTransform)
#dev.copy(jpeg,'Figure_S3.jpg', width=7, height=5, units='in', res=300)
#dev.off()


##### Run Mixed Effects Models on Stomach Fullness to evaluate impact of month and heatwave status (Table S3) ####

#No Random Effects, all Fixed Effects and Interactions
S10 <- lm(sqrtFullness ~ Heatwave  * Month, na.action = na.exclude, Fullness) 
summary(S10)
Anova(S10, type =3)
r.squaredGLMM(S10)

#Year as a random Intercept
S20<- lmer(sqrtFullness ~  Heatwave  * Month +  (1|Year), na.action = na.exclude, Fullness, control = lmerControl(optimizer = "Nelder_Mead"), REML =TRUE) 
summary(S20) 
Anova(S20, type =3)

AIC(S10, S20)
r.squaredGLMM(S20)

#Percent Increase Between Stomach Fullness in July and August
FullnessMean <-Fullness %>%
  mutate(MonthName = case_when(Month == 7 ~ "July", Month == 8 ~ "August")) %>%
  group_by(Year, MonthName) %>%
  reframe(meanFullness = mean(sqrtFullness, na.rm = T)) %>%
  pivot_wider(names_from = c(MonthName), values_from = meanFullness) %>% 
  mutate(MeanPercentIncrease_July_Aug = ((August - July)/July) * 100) %>%
  ungroup() %>%
  mutate(Heatwave = case_when(
    Year == '2007'|Year == '2009' | Year == '2010' | Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  )) %>%
  group_by( Heatwave) %>%
  reframe(MeanPercentIncrease = mean(MeanPercentIncrease_July_Aug), se = (sd(MeanPercentIncrease_July_Aug)/sqrt(length(MeanPercentIncrease_July_Aug))))
FullnessMean

#### Calculate prey IRI (PSIRI) for juvenile Pacific Cod stomachs ####
#following Brown et al. 2012

StomachContents <- read.csv("D7.PreyData_Raw.csv") %>%
  select(FISHID, Year, Month,PreySpecies, PreyCount, Prey_WW_g) %>%
  mutate(FISHID = as.factor(FISHID)) %>%
  filter(!PreySpecies %in% 'EMPTY') %>%
  filter(!PreySpecies %in% 'NEMATODA') %>%
  filter(!PreySpecies %in% 'PLANTAE') %>%
  filter(!PreySpecies %in% 'PLASTIC') %>%
  filter(!PreySpecies %in% 'ROCK') %>%
  filter(!PreySpecies %in% 'UNIDENTIFIED') %>%
  filter(!PreySpecies %in% 'WOOD') %>%
  mutate(Heatwave = case_when(Year == '2006'|Year == '2007'|Year == '2008'|Year == '2009' | Year == '2010' | Year == '2011'| Year == '2012'| Year == '2013'~ "Before", Year == '2014' | Year == '2015' | Year =='2016'| Year == '2019' ~ "Heatwave", Year == "2017" | Year == "2018" ~ "Between"
  ))


#Run Separate pIRI and IRI for each month

Aug_Stomachs <- StomachContents %>%
  filter(Month == 8)

July_Stomachs <- StomachContents %>%
  filter(Month == 7)

##### August IRI and Prey-Specific IRI ####

All_PreyCounts_Fish <- Aug_Stomachs %>% group_by(FISHID) %>%  #Total counts of prey per stomach
  summarise(AllPreyNo = sum(PreyCount)) 

All_PreyWeights_Fish <- Aug_Stomachs %>% group_by(FISHID) %>%  #Total weights of prey per stomach
  summarise(AllPreyWt = sum(Prey_WW_g)) 

Aug.1 <- left_join(All_PreyCounts_Fish, All_PreyWeights_Fish, by = c("FISHID")) #Join total counts and weights

Aug.2 <- left_join(Aug_Stomachs,Aug.1, by = c("FISHID")) #Join with dataframe


# Calculate % Frequency of Occurrence

Aug_FrequencyTable = as.data.frame(table(Aug.2$PreySpecies)) %>%
  mutate(Freq = as.numeric(Freq)) %>%
  rename(PreySpecies = Var1)

Aug.3 <- left_join(Aug.2,Aug_FrequencyTable, by = c("PreySpecies"))

#Calculate IRI
Aug.IRI <- Aug.3 %>% 
  group_by(FISHID) %>%
  mutate(NC = (PreyCount/AllPreyNo)* 100) %>%
  mutate(GC = (Prey_WW_g/AllPreyWt) * 100) %>%
  mutate(FO = (Freq/length(unique(Aug_Stomachs$FISHID)))) %>%
  mutate(IRI = (NC + GC) * FO) %>%
  mutate(percent_IRI = 100 * IRI/sum(IRI))


# Calculate PSIRI

PSIRI_Spp_Aug <- Aug.3 %>% 
  group_by(FISHID, PreySpecies) %>%  
  mutate(PN = ((PreyCount/AllPreyNo)*100)/Freq) %>% #Prey-specific counts
  mutate(PW = ((Prey_WW_g/AllPreyWt)*100)/Freq) %>%  #Prey-specific weight
  ungroup() %>%
  group_by(PreySpecies,FISHID,Freq, Year, Month, Heatwave) %>% 
  summarise(PN = sum(PN),PW = sum(PW)) %>%
  mutate(FO = (Freq/length(unique(Aug_Stomachs$FISHID)))) %>% 
  mutate(PSIRI = (FO * (PN+PW))/2) %>% 
  group_by(FISHID) %>%
  mutate(percent_PSIRI =  PSIRI/sum(PSIRI)) %>%
  arrange(-PSIRI)  


Wide_Aug_PSIRI <- PSIRI_Spp_Aug %>%
  select(FISHID, Year, Month, Heatwave, PreySpecies, percent_PSIRI) %>%
  #group_by(FISHID, PreySpecies) %>%
  pivot_wider(names_from = PreySpecies, values_from = percent_PSIRI) %>%
  replace(is.na(.), 0) 


ColumnMeans <- Wide_Aug_PSIRI %>%
  ungroup() %>%
  select( -Year, -Month, -Heatwave, -FISHID) %>%
  colMeans() 

ColumnMeans  <- as.data.frame(ColumnMeans) %>%
  mutate(Percent = ColumnMeans * 100) %>%
  format(scientific = F) 
#write.csv(ColumnMeans, "AugColumnMeans_PSIRI_ALL.csv")
#Group similar taxa and add taxa < 3.5% PSIRI to other

Wide_Aug_PSIRI_Greaterthan3.5 <- Wide_Aug_PSIRI %>%
  mutate(CALANOIDA = ACARTIALONGIREMIS + CALANUS + COPEPODA + EURYTEMORA + METRIDIA + UNKNOWNSMALLCOPEPOD) %>%
  select(-ACARTIALONGIREMIS, -CALANUS, -COPEPODA, -EURYTEMORA, -METRIDIA, -UNKNOWNSMALLCOPEPOD) %>%
  mutate(GAMMARIDAE = GAMMARIDAE + CORPHIIDAE + UNKNOWNAMPHIPOD) %>%
  select(-CORPHIIDAE, - UNKNOWNAMPHIPOD) %>%
  mutate(SHRIMP = UNKNOWNPANDALID + CRANGONALASKENSIS + STRIPEDSHRIMP +  SPIRONTOCARIS ) %>%
  select(-UNKNOWNPANDALID, -CRANGONALASKENSIS, -STRIPEDSHRIMP, -SPIRONTOCARIS) %>%
  mutate(ANNELIDA = NEPHTYSSP  + SIPUNCULIDA + NEREISSP + POLYCHAETA ) %>%
  select( -NEPHTYSSP,  -SIPUNCULIDA, -NEREISSP, -POLYCHAETA) %>%
  mutate(CLADOCERA = CLADOCERA + ClADOCERA) %>%
  select(-ClADOCERA) %>%
  mutate(OTHER = ZOEA + OITHONIA + SALPS + BarnNaup + BarnCyp + OSTEICHTHYES + DECAPODA + CLUPEAPALLASII + MEGALOPE + INSECTA + MICROGADUSPROXIMUS + CIRRIPEDIA + SHRIMPLARVAE + ISOPODA+ HYPERIIDAE + MEGAYOLDIA + NEVERITA + SHRIMP + CUMACEA + PTEROPODA) %>%
  select(-ZOEA, -OITHONIA, -SALPS, -BarnNaup, -BarnCyp, -OSTEICHTHYES, -DECAPODA ,- CLUPEAPALLASII, -MEGALOPE, -INSECTA, -MICROGADUSPROXIMUS, -CIRRIPEDIA, -SHRIMPLARVAE, - ISOPODA, -HYPERIIDAE, -MEGAYOLDIA, -NEVERITA, -SHRIMP, -CUMACEA , -PTEROPODA)

#Create Data Sheet for NMS Ordinations
#write.csv(Wide_Aug_PSIRI_Greaterthan3.5, "AugPSIRI_ROutput.csv")

#Look at prey taxa percentages for August overall
ColumnMeans_Trim <- Wide_Aug_PSIRI_Greaterthan3.5 %>%
  ungroup() %>%
  select( -Year, -Month, -Heatwave, -FISHID) %>%
  colMeans() 
ColumnMeans_Trim  <- as.data.frame(ColumnMeans_Trim) %>%
  mutate(Percent = ColumnMeans_Trim * 100) %>%
  format(scientific = F) 

#Look at prey taxa percentages for August by Heatwave class
PSIRI_Aug_Trim_byHW <- Wide_Aug_PSIRI_Greaterthan3.5 %>%
  group_by(Heatwave) %>%
  summarise_all("mean") 

##### July IRI and Prey-Specific IRI ####

All_PreyCounts_July <- July_Stomachs %>% group_by(FISHID) %>%  #Total counts of prey per stomach
  summarise(AllPreyNo = sum(PreyCount)) 

All_PreyWeights_July <- July_Stomachs %>% group_by(FISHID) %>%  #Total weights of prey per stomach
  summarise(AllPreyWt = sum(Prey_WW_g)) 

July.1 <- left_join(All_PreyCounts_July, All_PreyWeights_July, by = c("FISHID")) #Join total counts and weights

July.2 <- left_join(July_Stomachs,July.1, by = c("FISHID")) #Join with dataframe


# Calculate % Frequency of Occurrence

July_FrequencyTable = as.data.frame(table(July.2$PreySpecies)) %>%
  mutate(Freq = as.numeric(Freq)) %>%
  rename(PreySpecies = Var1)

July.3 <- left_join(July.2,July_FrequencyTable, by = c("PreySpecies"))

#Calculate IRI
July.IRI <- July.3 %>% 
  group_by(FISHID) %>%
  mutate(NC = (PreyCount/AllPreyNo)* 100) %>%
  mutate(GC = (Prey_WW_g/AllPreyWt) * 100) %>%
  mutate(FO = (Freq/length(unique(July_Stomachs$FISHID)))) %>%
  mutate(IRI = (NC + GC) * FO) %>%
  mutate(percent_IRI = 100 * IRI/sum(IRI))


# Calculate PSIRI
PSIRI_Spp_July <- July.3 %>% 
  group_by(FISHID, PreySpecies) %>%  
  mutate(PN = ((PreyCount/AllPreyNo)*100)/Freq) %>% #Prey-specific counts
  mutate(PW = ((Prey_WW_g/AllPreyWt)*100)/Freq) %>%  #Prey-specific weight
  ungroup() %>%
  group_by(PreySpecies,FISHID,Freq, Year, Month, Heatwave) %>% 
  summarise(PN = sum(PN),PW = sum(PW)) %>%
  mutate(FO = (Freq/length(unique(Aug_Stomachs$FISHID)))) %>% 
  mutate(PSIRI = (FO * (PN+PW))/2) %>% 
  group_by(FISHID) %>%
  mutate(percent_PSIRI =  PSIRI/sum(PSIRI)) %>%
  arrange(-PSIRI)  


Wide_July_PSIRI <- PSIRI_Spp_July %>%
  select(FISHID, Year, Month, Heatwave, PreySpecies, percent_PSIRI) %>%
  pivot_wider(names_from = PreySpecies, values_from = percent_PSIRI) %>%
  replace(is.na(.), 0) 

ColumnMeans <- Wide_July_PSIRI %>%
  ungroup() %>%
  select( -Year, -Month, -Heatwave, -FISHID) %>%
  colMeans() 

ColumnMeans  <- as.data.frame(ColumnMeans) %>%
  mutate(Percent = ColumnMeans * 100) %>%
  format(scientific = F) 


#Add Taxa < 3.5% %PSIRI to Other
Wide_July_PSIRI_Greaterthan3.5 <- Wide_July_PSIRI %>%
  mutate(CALANOIDA = ACARTIALONGIREMIS + CALANUS +  EURYTEMORA + METRIDIA + UNKNOWNMEDCOPEPOD) %>%
  select(-ACARTIALONGIREMIS, -CALANUS, -EURYTEMORA, -METRIDIA, -UNKNOWNMEDCOPEPOD) %>%
  mutate(GAMMARIDAE = GAMMARIDAE + CORPHIIDAE ) %>%
  select(-CORPHIIDAE) %>%
  mutate(CLADOCERA = CLADOCERA + ClADOCERA) %>%
  select(-ClADOCERA) %>%
  mutate(OTHER = CIRRIPEDIA +  BarnNaup + BarnCyp  + NEVERITA + CUMACEA + OSTEICHTHYES + HYPERIIDAE  + TSPINIFERA  + POLYCHAETA + NEREISSP + NEPHTYSSP +SPIRONTOCARIS + UNKNOWNSHRIMP + INSECTA + MEGAYOLDIA + PTEROPODA + DECAPODA + SHRIMPLARVAE + ZOEA) %>%
  select(-CIRRIPEDIA ,- BarnNaup,-BarnCyp,- NEVERITA,-CUMACEA,-OSTEICHTHYES,- HYPERIIDAE,-TSPINIFERA,-POLYCHAETA,-NEREISSP,-NEPHTYSSP, -SPIRONTOCARIS,- UNKNOWNSHRIMP, -INSECTA, - MEGAYOLDIA, -PTEROPODA , -DECAPODA, -SHRIMPLARVAE , - ZOEA)
#write.csv(Wide_July_PSIRI_Greaterthan3.5, "JulyPSIRI_ROutput.csv")



#Look at prey taxa percentages for July overall
ColumnMeans_Trim <- Wide_July_PSIRI_Greaterthan3.5 %>%
  ungroup() %>%
  select( -Year, -Month, -Heatwave, -FISHID) %>%
  colMeans() 

ColumnMeans_Trim  <- as.data.frame(ColumnMeans_Trim) %>%
  mutate(Percent = ColumnMeans_Trim * 100) %>%
  format(scientific = F) 




#### Visualize diet composition using Mosaic Plots (Figures 2 and S4) ####

##### Mosaic Plots By Year (Figure S4) ####
July_PreyColumnAverages <- Wide_July_PSIRI_Greaterthan3.5 %>%
  group_by(Year) %>%
  reframe(CLADOCERA = mean(CLADOCERA), MYSIDAE = mean(MYSIDAE), GAMMARIDAE = mean(GAMMARIDAE), HARPACTICOIDA= mean(HARPACTICOIDA), CAPRELLIDAE = mean(CAPRELLIDAE), CALANOIDA = mean(CALANOIDA), OTHER = mean(OTHER)) 

July2006 <-data.frame(2006,0,0,0,0,0,0,0)
names(July2006)<-c("Year","CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "OTHER")

July2008 <-data.frame(2008,0,0,0,0,0,0,0)
names(July2008)<-c("Year","CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "OTHER")

July2011 <-data.frame(2011,0,0,0,0,0,0,0)
names(July2011)<-c("Year","CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "OTHER")

JulyYearlyPrey_wide <- rbind(July_PreyColumnAverages, July2006, July2008, July2011)

JulyMosaic <- JulyYearlyPrey_wide %>%
  pivot_longer(cols=c("CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "OTHER"), names_to='PreyType', values_to='PSIRI') %>%
  mutate(Heatwave = case_when(
    Year == 2006 | Year == 2007 | Year == 2008 | Year == 2009 | Year == 2010 | Year == 2011 | Year == 2012| Year == 2013 ~ "Before",
    Year == 2014 | Year == 2015 | Year == 2016 | Year == 2019~ "During",
    Year == 2017 | Year == 2018 ~ "Between")) %>%
  mutate(Year = as.factor(Year))

JulyMosaic$PreyType <-factor(JulyMosaic$PreyType, levels = c("CLADOCERA", "GAMMARIDAE", "HARPACTICOIDA", "CALANOIDA", "CAPRELLIDAE", "MYSIDAE",  "OTHER", "ANNELIDA"))

pal = pnw_palette("Bay", 8, type = "continuous")

labela <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), y = rep(1.04, 14),labels = c("n = 0", "23", "0", "27", "25", "0", "25", "25", "26", "9", "25", "26", "25", "25"))



July_PSIRI_Mosaic_byYear <- ggplot()  +
  geom_bar(aes(y = PSIRI, x = Year, fill = PreyType), data = JulyMosaic, stat="identity") +
  xlab("Year") + 
  ylab("Diet composition by % PSIRI") +
  theme_classic() +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())  +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.4)) +
  geom_text(data = labela, aes(x = x, y = y, label = labels), angle = 0, size = 5, colour = c("#332288","#332288","#332288","#332288","#332288","#332288","#332288","#332288", "#CC6677", "#CC6677", "#CC6677","#888888","#888888", "#CC6677")) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(colour = c("#332288","#332288","#332288","#332288","#332288","#332288","#332288","#332288", "#CC6677", "#CC6677", "#CC6677","#888888","#888888", "#CC6677"))) +
  ggtitle("July")+
  scale_fill_manual(values =pal) +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) 
print(July_PSIRI_Mosaic_byYear)



Aug_PreyColumnAverages <- Wide_Aug_PSIRI_Greaterthan3.5 %>%
  group_by(Year) %>%
  reframe(CLADOCERA = mean(CLADOCERA), MYSIDAE = mean(MYSIDAE), GAMMARIDAE = mean(GAMMARIDAE), HARPACTICOIDA= mean(HARPACTICOIDA), CAPRELLIDAE = mean(CAPRELLIDAE), CALANOIDA = mean(CALANOIDA), ANNELIDA = mean(ANNELIDA), OTHER = mean(OTHER)) 

Aug2011 <-data.frame(2011,0,0,0,0,0,0,0, 0)
names(Aug2011)<-c("Year","CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "ANNELIDA", "OTHER")

Aug2016 <-data.frame(2016,0,0,0,0,0,0,0, 0)
names(Aug2016)<-c("Year","CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "ANNELIDA", "OTHER")

AugYearlyPrey_wide <- rbind(Aug_PreyColumnAverages, Aug2011, Aug2016)

AugMosaic <- AugYearlyPrey_wide %>%
  pivot_longer(cols=c("CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "ANNELIDA", "OTHER"), names_to='PreyType', values_to='PSIRI') %>%
  mutate(Heatwave = case_when(
    Year == 2006 | Year == 2007 | Year == 2008 | Year == 2009 | Year == 2010 | Year == 2011 | Year == 2012| Year == 2013 ~ "Before",
    Year == 2014 | Year == 2015 | Year == 2016 | Year == 2019~ "During",
    Year == 2017 | Year == 2018 ~ "Between")) %>%
  mutate(Year = as.factor(Year))

AugMosaic$PreyType <-factor(AugMosaic$PreyType, levels = c("CLADOCERA", "GAMMARIDAE", "HARPACTICOIDA", "CALANOIDA", "CAPRELLIDAE", "MYSIDAE",  "OTHER", "ANNELIDA"))

pal = pnw_palette("Bay", 8, type = "continuous")

labelb <- data.frame(x = c(1:14), y = rep(1.04, 14),labels = c("n = 27", "21", "25", "26", "12", "0", "26", "25", "3", "24", "0", "25", "25", "26"))


Aug_PSIRI_Mosaic_byYear <- ggplot()  +
  geom_bar(aes(y = PSIRI, x = Year, fill = PreyType), data = AugMosaic, stat="identity") +
  xlab("Year") + 
  ylab("Diet composition by % PSIRI") +
  theme_classic() +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())  +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.4)) +
  geom_text(data = labelb, aes(x = x, y = y, label = labels), angle = 0, size = 5, colour  = c("#332288","#332288","#332288","#332288","#332288","#332288","#332288","#332288", "#CC6677", "#CC6677", "#CC6677","#888888","#888888", "#CC6677")) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(colour = c("#332288","#332288","#332288","#332288","#332288","#332288","#332288","#332288", "#CC6677", "#CC6677", "#CC6677","#888888","#888888", "#CC6677"))) +
  ggtitle("August")+
  scale_fill_manual(values =pal) +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) 
print(Aug_PSIRI_Mosaic_byYear)


### Combined Mosaic Plot by Year 
Mosaic_PSIRI_byYear <- July_PSIRI_Mosaic_byYear + Aug_PSIRI_Mosaic_byYear  +
  plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(Mosaic_PSIRI_byYear)
#dev.copy(jpeg,'Figure_S4.jpg', width=10, height=11, units='in', res=300)
#dev.off()

##### Mosaic Plots By Heatwave. (Figure 2) ####
July_PreyColumnAverages_HW <- Wide_July_PSIRI_Greaterthan3.5 %>%
  group_by(Heatwave) %>%
  reframe(CLADOCERA = mean(CLADOCERA), MYSIDAE = mean(MYSIDAE), GAMMARIDAE = mean(GAMMARIDAE), HARPACTICOIDA= mean(HARPACTICOIDA), CAPRELLIDAE = mean(CAPRELLIDAE), CALANOIDA = mean(CALANOIDA), OTHER = mean(OTHER)) 

JulyMosaic_HW <- July_PreyColumnAverages_HW  %>%
  pivot_longer(cols=c("CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "OTHER"), names_to='PreyType', values_to='PSIRI') 

JulyMosaic_HW$PreyType <-factor(JulyMosaic_HW$PreyType, levels = c("CLADOCERA", "GAMMARIDAE", "HARPACTICOIDA", "CALANOIDA", "CAPRELLIDAE", "MYSIDAE",  "OTHER", "ANNELIDA"))

pal = pnw_palette("Bay", 8, type = "continuous")

labelz <- data.frame(x = c(1, 2, 3), y = rep(1.04, 3),labels = c("n = 125", "51", "60"))

July_PSIRI_Mosaic_HW <- ggplot()  +
  geom_bar(aes(y = PSIRI, x = Heatwave, fill = PreyType), data = JulyMosaic_HW, stat="identity") +
  xlab("Heatwave") + 
  ylab("Diet composition by % PSIRI") +
  theme_classic() +
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())  +
  theme(legend.position="none") +
  geom_text(data = labelz, aes(x = x, y = y, label = labels), angle = 0, size = 5) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold")) +
  ggtitle("July")+
  scale_fill_manual(values =pal) +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) 
print(July_PSIRI_Mosaic_HW)

Aug_PreyColumnAverages_HW <- Wide_Aug_PSIRI_Greaterthan3.5 %>%
  group_by(Heatwave) %>%
  reframe(CLADOCERA = mean(CLADOCERA), MYSIDAE = mean(MYSIDAE), GAMMARIDAE = mean(GAMMARIDAE), HARPACTICOIDA= mean(HARPACTICOIDA), CAPRELLIDAE = mean(CAPRELLIDAE), CALANOIDA = mean(CALANOIDA), ANNELIDA = mean(ANNELIDA), OTHER = mean(OTHER)) 

AugMosaic_HW <- Aug_PreyColumnAverages_HW %>%
  pivot_longer(cols=c("CLADOCERA", "MYSIDAE", "GAMMARIDAE", "HARPACTICOIDA", "CAPRELLIDAE", "CALANOIDA", "ANNELIDA", "OTHER"), names_to='PreyType', values_to='PSIRI')

AugMosaic_HW$PreyType <-factor(AugMosaic_HW$PreyType, levels = c("CLADOCERA", "GAMMARIDAE", "HARPACTICOIDA", "CALANOIDA", "CAPRELLIDAE", "MYSIDAE",  "OTHER", "ANNELIDA"))

pal = pnw_palette("Bay", 8, type = "continuous")

labely <- data.frame(x = c(1:3), y = rep(1.04, 3),labels = c("n = 162", "50", "53"))

Aug_PSIRI_Mosaic_HW <- ggplot()  +
  geom_bar(aes(y = PSIRI, x = Heatwave, fill = PreyType), data = AugMosaic_HW, stat="identity") +
  xlab("Heatwave") + 
  ylab("Diet composition by % PSIRI") +
  theme_classic() +
  theme(legend.position="right", legend.direction="vertical", legend.title = element_blank())  +
  geom_text(data = labely, aes(x = x, y = y, label = labels), angle = 0, size = 5) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold")) +
  ggtitle("August")+
  scale_fill_manual(values =pal) +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) 
print(Aug_PSIRI_Mosaic_HW)


### Combined Mosaic Plot by HW
Mosaic_PSIRI_HW <- July_PSIRI_Mosaic_HW + Aug_PSIRI_Mosaic_HW  +
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16)) 
print(Mosaic_PSIRI_HW)
#dev.copy(jpeg,'Figure_2.jpg', width=12, height=7, units='in', res=300)
#dev.off()

#### Year-Specific Comparisons ####
##### 2019 #####

##Year Specific: Look at all the prey
### Just Look at 2019
ColumnMeans_July2019 <- Wide_July_PSIRI %>%
  ungroup() %>%
  filter(Year == "2019") %>%
  select( -Year, -Month, -Heatwave, -FISHID) %>%
  colMeans() 

ColumnMeans_July2019  <- as.data.frame(ColumnMeans_July2019) %>%
  mutate(Percent = ColumnMeans_July2019 * 100) %>%
  format(scientific = F) 

ColumnMeans_2019 <- Wide_Aug_PSIRI %>%
  ungroup() %>%
  filter(Year == "2019") %>%
  select( -Year, -Month, -Heatwave, -FISHID) %>%
  colMeans() 

ColumnMeans_2019  <- as.data.frame(ColumnMeans_2019) %>%
  mutate(Percent = ColumnMeans_2019 * 100) %>%
  format(scientific = F) 

Wide_Aug_PSIRI_2019 <- Wide_Aug_PSIRI %>%
  filter(Year == "2019") %>%
  mutate(CALANOIDA = ACARTIALONGIREMIS + CALANUS + COPEPODA + EURYTEMORA + METRIDIA + UNKNOWNSMALLCOPEPOD) %>%
  select(-ACARTIALONGIREMIS, -CALANUS, -COPEPODA, -EURYTEMORA, -METRIDIA, -UNKNOWNSMALLCOPEPOD) %>%
  mutate(GAMMARIDAE = GAMMARIDAE + CORPHIIDAE + UNKNOWNAMPHIPOD) %>%
  select(-CORPHIIDAE, - UNKNOWNAMPHIPOD) %>%
  mutate(SHRIMP = UNKNOWNPANDALID + CRANGONALASKENSIS + STRIPEDSHRIMP +  SPIRONTOCARIS ) %>%
  select(-UNKNOWNPANDALID, -CRANGONALASKENSIS, -STRIPEDSHRIMP, -SPIRONTOCARIS) %>%
  mutate(ANNELIDA = NEPHTYSSP  + SIPUNCULIDA + NEREISSP + POLYCHAETA ) %>%
  select( -NEPHTYSSP,  -SIPUNCULIDA, -NEREISSP, -POLYCHAETA) %>%
  mutate(CLADOCERA = CLADOCERA + ClADOCERA) %>%
  select(-ClADOCERA) %>%
  mutate(BARNACLELARVAE = BarnCyp + BarnNaup +  CIRRIPEDIA) %>%
  select(-BarnNaup, - BarnCyp, -CIRRIPEDIA) %>%
  select(-OSTEICHTHYES, -SHRIMPLARVAE, -INSECTA, -CLUPEAPALLASII, -MICROGADUSPROXIMUS, -DECAPODA, -MEGALOPE, -ISOPODA, -OITHONIA, -SALPS )
  

Wide_July_PSIRI_2019 <- Wide_July_PSIRI %>%
  filter(Year == "2019") %>%
  mutate(CALANOIDA = ACARTIALONGIREMIS + CALANUS  + EURYTEMORA + METRIDIA ) %>%
  select(-ACARTIALONGIREMIS, -CALANUS, -EURYTEMORA, -METRIDIA) %>%
  mutate(GAMMARIDAE = GAMMARIDAE + CORPHIIDAE) %>%
  select(-CORPHIIDAE) %>%
  mutate(SHRIMP =    SPIRONTOCARIS ) %>%
  select(   -SPIRONTOCARIS) %>%
  mutate(ANNELIDA = NEPHTYSSP + NEREISSP + POLYCHAETA ) %>%
  select( -NEPHTYSSP,   -NEREISSP, -POLYCHAETA) %>%
  mutate(CLADOCERA = CLADOCERA + ClADOCERA) %>%
  select(-ClADOCERA) %>%
  mutate(BARNACLELARVAE = BarnCyp + BarnNaup +  CIRRIPEDIA) %>%
  select(-BarnNaup, - BarnCyp, -CIRRIPEDIA) %>%
  select(-OSTEICHTHYES, -TSPINIFERA, -SHRIMPLARVAE, -UNKNOWNMEDCOPEPOD, -UNKNOWNSHRIMP, -DECAPODA, -INSECTA)

head(Wide_July_PSIRI_2019)

##### Mosaic Plots By Year (Figure S4) ####
July_PreyColumnAverages_2019 <- Wide_July_PSIRI_2019 %>%
  ungroup() %>%
  select(-FISHID, -Month, -Heatwave, -Year) %>%
  summarise_all(list(mean = ~mean(., na.rm = T))) %>%
  mutate(CLADOCERA = CLADOCERA_mean) %>%
  mutate(GAMMARIDAE = GAMMARIDAE_mean) %>%
  mutate(HARPACTICOIDA = HARPACTICOIDA_mean) %>%
  mutate(MYSIDAE = MYSIDAE_mean) %>% 
  mutate(PTEROPODA = PTEROPODA_mean) %>%
  mutate(CAPRELLIDAE = CAPRELLIDAE_mean) %>%
  mutate(CUMACEA = CUMACEA_mean) %>%
  mutate(NEVERITA = NEVERITA_mean) %>%
  mutate(HYPERIIDAE = HYPERIIDAE_mean) %>%
  mutate(MEGAYOLDIA = MEGAYOLDIA_mean) %>%
  mutate(ZOEA = ZOEA_mean) %>%
  mutate(CALANOIDA = CALANOIDA_mean) %>%
  mutate(SHRIMP = SHRIMP_mean) %>%
  mutate(ANNELIDA = ANNELIDA_mean) %>%
  mutate(BARNACLELARVAE = BARNACLELARVAE_mean) %>%
  select(-CLADOCERA_mean, -GAMMARIDAE_mean, -HARPACTICOIDA_mean, -MYSIDAE_mean, -PTEROPODA_mean, -CAPRELLIDAE_mean, -CUMACEA_mean, -NEVERITA_mean, -HYPERIIDAE_mean, -MEGAYOLDIA_mean, -ZOEA_mean, -CALANOIDA_mean, -SHRIMP_mean, -ANNELIDA_mean, -BARNACLELARVAE_mean)%>%
  mutate(Year = "2019") %>%
  mutate(Month = "July")

Aug_PreyColumnAverages_2019 <- Wide_Aug_PSIRI_2019 %>%
  ungroup() %>%
  select(-FISHID, -Month, -Heatwave, -Year) %>%
  summarise_all(list(mean = ~mean(., na.rm = T))) %>%
  mutate(CLADOCERA = CLADOCERA_mean) %>%
  mutate(GAMMARIDAE = GAMMARIDAE_mean) %>%
  mutate(HARPACTICOIDA = HARPACTICOIDA_mean) %>%
  mutate(MYSIDAE = MYSIDAE_mean) %>% 
  mutate(PTEROPODA = PTEROPODA_mean) %>%
  mutate(CAPRELLIDAE = CAPRELLIDAE_mean) %>%
  mutate(CUMACEA = CUMACEA_mean) %>%
  mutate(NEVERITA = NEVERITA_mean) %>%
  mutate(HYPERIIDAE = HYPERIIDAE_mean) %>%
  mutate(MEGAYOLDIA = MEGAYOLDIA_mean) %>%
  mutate(ZOEA = ZOEA_mean) %>%
  mutate(CALANOIDA = CALANOIDA_mean) %>%
  mutate(SHRIMP = SHRIMP_mean) %>%
  mutate(ANNELIDA = ANNELIDA_mean) %>%
  mutate(BARNACLELARVAE = BARNACLELARVAE_mean) %>%
  select(-CLADOCERA_mean, -GAMMARIDAE_mean, -HARPACTICOIDA_mean, -MYSIDAE_mean, -PTEROPODA_mean, -CAPRELLIDAE_mean, -CUMACEA_mean, -NEVERITA_mean, -HYPERIIDAE_mean, -MEGAYOLDIA_mean, -ZOEA_mean, -CALANOIDA_mean, -SHRIMP_mean, -ANNELIDA_mean, -BARNACLELARVAE_mean)%>%
  mutate(Year = "2019") %>%
  mutate(Month = "August")

Diet2019 <- rbind(July_PreyColumnAverages_2019, Aug_PreyColumnAverages_2019)
  
Mosaic2019 <- Diet2019 %>%
  pivot_longer(cols=c("CLADOCERA", "GAMMARIDAE", "HARPACTICOIDA", "MYSIDAE", "PTEROPODA", "CAPRELLIDAE", "CUMACEA", "NEVERITA", "HYPERIIDAE", "MEGAYOLDIA", "ZOEA", "CALANOIDA", "SHRIMP", "ANNELIDA", "BARNACLELARVAE"), names_to='PreyType', values_to='PSIRI') %>%
  mutate(Month = factor(Month, levels = c("July", "August")))
  

#library(pals)
DietPlot_2019 <- ggplot()  +
  geom_bar(aes(y = PSIRI, x = Month, fill = PreyType), data = Mosaic2019, stat="identity") +
  xlab("Month") + 
  ylab("Diet composition by % PSIRI") +
  theme_classic() +
  ggtitle("2019") + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_blank(), plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  scale_fill_manual( values=unname(trubetskoy())) 
print(DietPlot_2019)
#dev.copy(jpeg,'2019_CompleteDiet.jpg', width=8, height = 6, units='in', res=300)
#dev.off()
