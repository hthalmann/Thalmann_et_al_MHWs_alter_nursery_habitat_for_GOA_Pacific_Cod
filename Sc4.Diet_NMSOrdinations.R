## Kodiak Island Juvenile Pacific Cod Diet Non-metric Multidimensional Scaling Ordinations

## Date Created: 1/11/2024
## R version: 4.2.2

## Load Libraries
library(rstudioapi) #load R working directory
library(vegan) #Run ordinations
library(tidyverse)  #data wrangling and visualization
library(indicspecies) #indicator species analysis
library(goeveg) #scree plots
library(corrplot) #correlation matrices
library(patchwork) #figures


## Paired with:
#Prey-Specific Index of Relative Importance (PSIRI) data matrices for July and August
#Trident Temperature Data
#Size and Condition Data

##Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#Sets directory to folder that this R script is saved to

#### Data Wrangling ####
JulyPSIRI <- read.csv("D8.JulyPSIRI_ROutput.csv")
AugPSIRI <- read.csv("D9.AugPSIRI_ROutput.csv")
Size <- read.csv("D6.Size_Condition_ROutput.csv") %>%
  select(FISHID, SL_mm, WholeBodyWW_g, HSI, LWResiduals_All)
Sampling <- read.csv("D10.SamplingData_Master.csv") %>%
  select(FISHID, SeineID, SITEID)
Temp <- read.csv("D4.TridentBay_DailyTemperature.csv") %>%
  mutate(Date = ymd(Date)) %>%
  select(Date, TridentTemp)
Seines <- read.csv("D2.SeineData.csv") %>%
  mutate(Date = ymd(Date)) %>%
  select(SeineID, Date)
Sites <- read.csv("D1.SiteData.csv") %>%
  select(SITEID, Site, Bay, Habitat)

JulyOrd <- JulyPSIRI %>%
  left_join(Size) %>%
  left_join(Sampling) %>%
  left_join(Seines) %>%
  left_join(Temp) %>%
  left_join(Sites)

AugOrd <- AugPSIRI %>%
  left_join(Size) %>%
  left_join(Sampling) %>%
  left_join(Seines) %>%
  left_join(Temp) %>%
  left_join(Sites)

#### Run July Ordination ####

## Create Environmental and Species Matrices for July 

JulyPreyData <- JulyOrd %>%
  select(CALANOIDA, CAPRELLIDAE, CLADOCERA, GAMMARIDAE, HARPACTICOIDA, MYSIDAE,  OTHER)
head(JulyPreyData)

JulyEnvData <- JulyOrd %>%
  select(Month, Year, Heatwave, SL_mm, HSI, LWResiduals_All, WholeBodyWW_g, Site, Bay, Habitat, TridentTemp)


## Run a Scree Plot using the package "goeveg"
JulyScree <- dimcheckMDS(JulyPreyData,distance = "bray",k = 6,trymax = 200,autotransform = FALSE)

Stress <- as.vector(JulyScree)
Dimensions <- 1:6
JulyScreePlot <-as.data.frame(cbind(Stress, Dimensions))

JulyScree_Plot <- ggplot(JulyScreePlot, aes(x=Dimensions, y=Stress)) + 
  geom_line() +
  geom_point(size = 3, color = "red") +
  ggtitle("July: Stress value in tested dimensions")+
  theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) + 
  ylab("Stress") + 
  xlab("Dimensions")
print(JulyScree_Plot)

#3 dimensions seems reasonable

#### Run the July NMS Ordination ####
JulyMDS<-metaMDS(JulyPreyData, distance="bray", k=3, trymax=200, autotransform=FALSE)

#Stress
JulyMDS$stress

#Shepard Plot to evaluate Ordination Fit
JulyShepardPlot <- stressplot(object = JulyMDS, lwd = 5, main = "July Shepard Plot")

diss <- as.vector(JulyMDS$diss)
dist <- as.vector(JulyMDS$dist)
JulyShepard <- as.data.frame(cbind(diss, dist))

JulyShepard_Plot <- ggplot(JulyShepard, aes(x=diss, y=dist)) + 
  geom_point(size = 2, color = "navy") +
  ggtitle("July: Shepard Plot")+
  theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) + 
  ylab("Ordination Distance") + 
  xlab("Observed Dissimilarity")
print(JulyShepard_Plot)


#Rotate to align with fish size
RotateJuly <- MDSrotate(JulyMDS, JulyEnvData$SL_mm) 
plot(JulyMDS)

## Set up your three NMS Principal Axes 
NMDS1 <- RotateJuly$points[,1]
NMDS2 <- RotateJuly$points[,2]
NMDS3 <- RotateJuly$points[,3]


#Add them to the ordination data plot
Julypcod.plot<-cbind(JulyOrd, NMDS1, NMDS2, NMDS3)
head(Julypcod.plot)

#R Output of Ordination for Subsequent Analysis
#write.csv(Julypcod.plot, "JulyNMS_ROutput.csv") 


## Plot to assess proportion of variance explained by each point

gof <- goodness(object = RotateJuly)

July_GoodnessofFit <- ggplot(data=Julypcod.plot, aes(NMDS1, NMDS2))+
  geom_point(data=Julypcod.plot, aes(NMDS1, NMDS2, color=Heatwave), cex = (4 *gof/mean(gof)),  show.legend=F, position=position_jitter(.1))+
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_color_manual(values=c("#332288", "#888888", "#CC6677")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("July: Goodness of Fit")
print(July_GoodnessofFit)

##### Plot July NMS with Prey (Figure 3a) ####
#Take out Other 
JulyPrey_Trim <- JulyPreyData %>%
  select(-OTHER)
en_prey_July <- envfit(RotateJuly, JulyPrey_Trim, permutations = 999, na.rm = T)
prey_coord_July = as.data.frame(scores(en_prey_July, "vectors")) * ordiArrowMul(en_prey_July)

July_prey <- ggplot(data=Julypcod.plot, aes(NMDS1, NMDS2))+
  geom_point(data=Julypcod.plot, aes(NMDS1, NMDS2, color=Heatwave), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=Heatwave, color = Heatwave), alpha=.2,type='t',linewidth =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = prey_coord_July, linewidth = 1, alpha = 0.8, colour = "black", arrow = arrow()) +
  geom_label(data = prey_coord_July, aes(x=NMDS1, y = NMDS2), colour = "black", fontface = "bold", label = row.names(prey_coord_July), position=position_jitter(0.26),  alpha = 0.6, size = 3)+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("July Prey Vectors")
print(July_prey)


##### Plot July NMS with Environmental Variables #####

## Look at the effects of environmental variables on axes
JulyEnvData_ContinuousOnly <- JulyEnvData %>%
  select(SL_mm, HSI, LWResiduals_All, TridentTemp)
en_July_env <- envfit(RotateJuly, JulyEnvData_ContinuousOnly, permutations = 999, na.rm = T)
# Get the vectors the correct length
env_coord_July = as.data.frame(scores(en_July_env, "vectors")) * ordiArrowMul(en_July_env)

July_Env <- ggplot(data=Julypcod.plot, aes(NMDS1, NMDS2))+
  geom_point(data=Julypcod.plot, aes(NMDS1, NMDS2, color=Heatwave), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=Heatwave, color = Heatwave), alpha=.2,type='t',linewidth =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = env_coord_July, linewidth = 1, alpha = 0.8, colour = "black", arrow = arrow()) +
  geom_label(data =env_coord_July, aes(x=NMDS1, y = NMDS2), colour = "black", fontface = "bold", label = row.names(env_coord_July), position=position_jitter(0.20),  alpha = 0.6)+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("July Environmental Vectors")
print(July_Env)

##### Plot July Size as a Contour (Figure 3b) #####

#Quick and dirty visualization
ordisurf(RotateJuly,JulyEnvData$SL_mm,main="",col="forestgreen ", lwd = 2, font.lab = 2 )
orditorp(RotateJuly,display="species",col="grey30",air=0.1,cex=1)


species.scores1_July <- as.data.frame(scores(RotateJuly, "species"))
species.scores1_July$species <- rownames(species.scores1_July)
species.scores1_July$z <- NA
head(species.scores1_July)

lengthcontour1_July <- ordisurf(RotateJuly ~ JulyEnvData$SL_mm, plot = FALSE, scaling =3)
head(lengthcontour1_July)

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("NMDS1", "NMDS2", "z")
  return(xyz)
}

lengthcontour1_July$grid
length_contour.vals1_July <- extract.xyz(obj = lengthcontour1_July)
head(length_contour.vals1_July)


Julypcod.plot$z = rep("NA", length(Julypcod.plot$NMDS1))
head(Julypcod.plot)
head(species.scores1_July)

en_coord_cont2_July <- cbind(env_coord_July, rep("NA",4))
names(en_coord_cont2_July)[3]<-paste("z")

labelx <- data.frame(x = c( -0.62, 0, 0.7, 1.25, 1.60), y = c( 1.06, 1.38, 1.10, 1.10, 0.8), z = NA, labels = c("45 mm", "50 mm", "55 mm", "60 mm", "65 mm"))


July_Size <- ggplot(data = length_contour.vals1_July, aes(NMDS1, NMDS2, z = z)) +
  geom_point(data=Julypcod.plot, aes(NMDS1, NMDS2, color = Heatwave))+ 
  stat_ellipse(data = Julypcod.plot, aes(fill=Heatwave), alpha=.2,type='t',linewidth =1, geom="polygon")+ 
  stat_contour(data = length_contour.vals1_July, binwidth = 5, na.rm = T, show.legend = T, lwd = 1.25, color = "black" ) +
  geom_text(data = labelx, aes(x = x, y = y, label = labels), angle = -10, size = 4) +
  #geom_label(data = species.scores1_July, aes(NMDS1, NMDS2, label = species), colour = "black", size = 3, position = position_jitter(0.6), alpha = 0.7) +
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  theme(legend.position="none") +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("July Size Contours")
print(July_Size)

##### Plot July Temp as a Contour (Figure 3c) #######

ordisurf(RotateJuly,JulyEnvData$TridentTemp,main="",col="forestgreen ", lwd = 2, font.lab = 2 )
orditorp(RotateJuly,display="species",col="grey30",air=0.1,cex=1)

## Manually import contour lines into ggplot 
tempcontour_July <- ordisurf(RotateJuly ~ JulyEnvData_ContinuousOnly$TridentTemp, plot = FALSE, scaling =3)

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("NMDS1", "NMDS2", "z")
  return(xyz)
}

tempcontour_July$grid
contour.vals_July <- extract.xyz(obj = tempcontour_July)

Julypcod.plot$z = rep("NA", length(Julypcod.plot$NMDS1))

labelz <- data.frame(x = c( 0, -0.75, 0.1), y = c(-1.1, 0.25, 0.1), z = NA, labels = c( "10.5°C", "9.5°C", "10°C"))


July_Temp <- ggplot(data = contour.vals_July, aes(NMDS1, NMDS2, z = z)) +
  geom_point(data=Julypcod.plot, aes(NMDS1, NMDS2, color = Heatwave))+ 
  stat_ellipse(data = Julypcod.plot, aes(fill=Heatwave), alpha=.2,type='t',size =1, geom="polygon")+ 
  stat_contour(data = contour.vals_July, binwidth = .5, na.rm = T, show.legend = T, lwd = 1.25, color = "black", lty = 1) +
  geom_text(data = labelz, aes(x = x, y = y, label = labels), angle = 0, size = 4) +
  #geom_label(data = species.scores1_July, aes(NMDS1, NMDS2, label = species), colour = "black",  position=position_jitter(0.6), size =3, alpha = 0.9) +
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  theme(legend.position="none") +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("July Temperature Contours")
print(July_Temp)

##### July Correlation Matrices ##### 
#For Species

JulySpeciesCorr <- Julypcod.plot %>%
  select(NMDS1, NMDS2, NMDS3, CALANOIDA, CAPRELLIDAE, CLADOCERA, GAMMARIDAE, HARPACTICOIDA, MYSIDAE, OTHER)
head(JulySpeciesCorr)
JulySpeciesCorr.cor <- cor(JulySpeciesCorr)
corrplot(JulySpeciesCorr.cor)


#For Environmental Variables
JulyEnvCorr <- Julypcod.plot %>%
  select(NMDS1, NMDS2, NMDS3, TridentTemp,  SL_mm, HSI, LWResiduals_All)
JulyEnvCorr_NAOmit <- na.omit(JulyEnvCorr)
JulyEnvCorr_NAOmit.cor = cor(JulyEnvCorr_NAOmit)
corrplot(JulyEnvCorr_NAOmit.cor)

cor(JulyEnvCorr_NAOmit$SL_mm, JulyEnvCorr_NAOmit$NMDS1)
cor(JulyEnvCorr_NAOmit$LWResiduals_All, JulyEnvCorr_NAOmit$NMDS1)
cor(JulyEnvCorr_NAOmit$HSI, JulyEnvCorr_NAOmit$NMDS1)
cor(JulyEnvCorr_NAOmit$TridentTemp, JulyEnvCorr_NAOmit$NMDS1)

##### July MRPP By Heatwave ####
mrpp(JulyPreyData, JulyEnvData$Heatwave, distance = 'bray', weight = 3)

##### July ISA By Heatwave ####
July_indval = multipatt(JulyPreyData, JulyEnvData$Heatwave, control = how(nperm=999))
summary(July_indval)


#### Run August Ordination ####

## Create Environmental and Species Matrices for August

AugPreyData <- AugOrd %>%
  select(CALANOIDA, CAPRELLIDAE, CLADOCERA, GAMMARIDAE, HARPACTICOIDA, MYSIDAE, ANNELIDA, OTHER)
head(AugPreyData)

AugEnvData <- AugOrd %>%
  select(Month, Year, Heatwave, SL_mm, HSI, LWResiduals_All, WholeBodyWW_g, Site, Bay, Habitat, TridentTemp)


## Run a Scree Plot using the package "goeveg"
AugScree <- dimcheckMDS(AugPreyData,distance = "bray",k = 6,trymax = 200,autotransform = FALSE)
#3 dimensions seems reasonable

Stress <- as.vector(AugScree)
Dimensions <- 1:6
AugScreePlot <-as.data.frame(cbind(Stress, Dimensions))

AugScree_Plot <- ggplot(AugScreePlot, aes(x=Dimensions, y=Stress)) + 
  geom_line() +
  geom_point(size = 3, color = "red") +
  ggtitle("August: Stress value in tested dimensions")+
  theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) + 
  ylab("Stress") + 
  xlab("Dimensions")
print(AugScree_Plot)

#### Run the August NMS Ordination ####
AugMDS<-metaMDS(AugPreyData, distance="bray", k=3, trymax=200, autotransform=FALSE)



#Stress
AugMDS$stress
#0.123

#Shepard Plot to evaluate Ordination Fit
stressplot(object = AugMDS, lwd = 5)


diss <- as.vector(AugMDS$diss)
dist <- as.vector(AugMDS$dist)
AugShepard <- as.data.frame(cbind(diss, dist))

AugShepard_Plot <- ggplot(AugShepard, aes(x=diss, y=dist)) + 
  geom_point(size = 2, color = "navy") +
  ggtitle("August: Shepard Plot")+
  theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5)) + 
  ylab("Ordination Distance") + 
  xlab("Observed Dissimilarity")
print(AugShepard_Plot)


#Rotate to align with fish size
RotateAug <- MDSrotate(AugMDS, AugEnvData$SL_mm) 

## Set up your three NMS Principal Axes 
NMDS1 <- RotateAug$points[,1]
NMDS2 <- RotateAug$points[,2]
NMDS3 <- RotateAug$points[,3]

#Add them to the ordination data plot
Augpcod.plot<-cbind(AugOrd, NMDS1, NMDS2, NMDS3)
head(Augpcod.plot)
#write.csv(Augpcod.plot, "AugNMS_ROutput.csv") 


## Plot to assess proportion of variance explained by each point

gof <- goodness(object = RotateAug)

Aug_GoodnessofFit <- ggplot(data=Augpcod.plot, aes(NMDS1, NMDS2))+
  geom_point(data=Augpcod.plot, aes(NMDS1, NMDS2, color=Heatwave), cex = (4 *gof/mean(gof)),  show.legend=F, position=position_jitter(.1))+
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_color_manual(values=c("#332288", "#888888", "#CC6677")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("August: Goodness of Fit")
print(Aug_GoodnessofFit)


##### Plot Aug NMS with Prey (Figure 3d) ####
#Take out Other 
AugPrey_Trim <- AugPreyData %>%
  select(-OTHER, -CAPRELLIDAE)
en_prey_Aug <- envfit(RotateAug, AugPrey_Trim, permutations = 999, na.rm = T)
prey_coord_Aug = as.data.frame(scores(en_prey_Aug, "vectors")) * ordiArrowMul(en_prey_Aug)


Aug_prey <- ggplot(data=Augpcod.plot, aes(NMDS1, NMDS2))+
  geom_point(data=Augpcod.plot, aes(NMDS1, NMDS2, color=Heatwave), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=Heatwave, color = Heatwave), alpha=.2,type='t',linewidth =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = prey_coord_Aug, linewidth = 1, alpha = 0.8, colour = "black", arrow = arrow()) +
  geom_label(data = prey_coord_Aug, aes(x=NMDS1, y = NMDS2), colour = "black", fontface = "bold", label = row.names(prey_coord_Aug), position=position_jitter(0.26),  alpha = 0.6, size = 3)+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("August Prey Vectors")
print(Aug_prey)


##### Plot Aug NMS with Environmental Variables #####

## Look at the effects of environmental variables on axes
AugEnvData_ContinuousOnly <- AugEnvData %>%
  select(SL_mm, HSI, LWResiduals_All, TridentTemp)
en_Aug_env <- envfit(RotateAug, AugEnvData_ContinuousOnly, permutations = 999, na.rm = T)
# Get the vectors the correct length
env_coord_Aug = as.data.frame(scores(en_Aug_env, "vectors")) * ordiArrowMul(en_Aug_env)

Aug_Env <- ggplot(data=Augpcod.plot, aes(NMDS1, NMDS2))+
  geom_point(data=Augpcod.plot, aes(NMDS1, NMDS2, color=Heatwave), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=Heatwave, color = Heatwave), alpha=.2,type='t',linewidth =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), data = env_coord_Aug, linewidth = 1, alpha = 0.8, colour = "black", arrow = arrow()) +
  geom_label(data =env_coord_Aug, aes(x=NMDS1, y = NMDS2), colour = "black", fontface = "bold", label = row.names(env_coord_Aug), position=position_jitter(0.20),  alpha = 0.6)+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("August Environmental Vectors")
print(Aug_Env)



##### Plot August Size as a Contour (Figure 3e) #####

#Quick and dirty visualization
ordisurf(RotateAug,AugEnvData$SL_mm,main="",col="forestgreen ", lwd = 2, font.lab = 2 )
orditorp(RotateAug,display="species",col="grey30",air=0.1,cex=1)


species.scores1_Aug <- as.data.frame(scores(RotateAug, "species"))
species.scores1_Aug$species <- rownames(species.scores1_Aug)
species.scores1_Aug$z <- NA
head(species.scores1_Aug)

lengthcontour1_Aug <- ordisurf(RotateAug ~ AugEnvData$SL_mm, plot = FALSE, scaling =3)
head(lengthcontour1_Aug)

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("NMDS1", "NMDS2", "z")
  return(xyz)
}

lengthcontour1_Aug$grid
length_contour.vals1_Aug <- extract.xyz(obj = lengthcontour1_Aug)
head(length_contour.vals1_Aug)


Augpcod.plot$z = rep("NA", length(Augpcod.plot$NMDS1))
head(Augpcod.plot)
head(species.scores1_Aug)

en_coord_cont2_Aug <- cbind(env_coord_Aug, rep("NA",4))
names(en_coord_cont2_Aug)[3]<-paste("z")

labely <- data.frame(x = c( 1.6, 1.8, 0.9, 0.3, -0.15, -1.2, -1.45, -1.40, -1.50), y = c( -0.75, 0.35, 1.3, 1.3, 1.2, 0.9, 0.35, 0, -0.40), z = NA, labels = c("95 mm", "90 mm", "85 mm", "80 mm", "75 mm", "70 mm", "65 mm", "60 mm", "55 mm"))


Aug_Size <- ggplot(data = length_contour.vals1_Aug, aes(NMDS1, NMDS2, z = z)) +
  geom_point(data=Augpcod.plot, aes(NMDS1, NMDS2, color = Heatwave))+ 
  stat_ellipse(data = Augpcod.plot, aes(fill=Heatwave), alpha=.2,type='t',size =1, geom="polygon")+ 
  stat_contour(data = length_contour.vals1_Aug, binwidth = 5, na.rm = T, show.legend = T, lwd = 1.25, color = "black" ) +
  geom_text(data = labely, aes(x = x, y = y, label = labels),   size = 4) + 
 # geom_label(data = species.scores1_Aug, aes(NMDS1, NMDS2, label = species), colour = "black", size = 3, position = position_jitter(0.5), alpha = 0.7) +
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  theme(legend.position="none") +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("August Size Contours")
print(Aug_Size)


###### Plot August Temp as a Contour (Figure 3f) #######

ordisurf(RotateAug,AugEnvData$TridentTemp,main="",col="forestgreen ", lwd = 2, font.lab = 2 )
orditorp(RotateAug,display="species",col="grey30",air=0.1,cex=1)

## Manually import contour lines into ggplot 
tempcontour_Aug <- ordisurf(RotateAug ~ AugEnvData_ContinuousOnly$TridentTemp, plot = FALSE, scaling =3)

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("NMDS1", "NMDS2", "z")
  return(xyz)
}

tempcontour_Aug$grid
contour.vals_Aug <- extract.xyz(obj = tempcontour_Aug)

Augpcod.plot$z = rep("NA", length(Augpcod.plot$NMDS1))

labelw <- data.frame(x = c( 1.65, 1.85, 1.85, 1.55, 1.1, -0.3, -1.45), y = c(-0.5, 0, 0.25, 0.75, 1, 0.4, 0.9), z = NA, labels = c( "12°C", "11.5°C", "11°C", "10.5°C", "10°C", "10°C", "9.5°C"))


Aug_Temp <- ggplot(data = contour.vals_Aug, aes(NMDS1, NMDS2, z = z)) +
  geom_point(data=Augpcod.plot, aes(NMDS1, NMDS2, color = Heatwave))+ 
  stat_ellipse(data = Augpcod.plot, aes(fill=Heatwave), alpha=.2,type='t',size =1, geom="polygon")+ 
  stat_contour(data = contour.vals_Aug, binwidth = .5, na.rm = T, show.legend = T, lwd = 1.25, color = "black", lty = 1) +
  geom_text(data = labelw, aes(x = x, y = y, label = labels), angle = 0, size = 4) +
  #geom_label(data = species.scores1_Aug, aes(NMDS1, NMDS2, label = species), colour = "black",  position=position_jitter(0.04), size =3, alpha = 0.9) +
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  theme(legend.position="none") +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("August Temperature Contours")
print(Aug_Temp)


##### August Correlation Matrices ##### 
#For Species

AugustSpeciesCorr <- Augpcod.plot %>%
  select(NMDS1, NMDS2, NMDS3, CALANOIDA, CAPRELLIDAE, CLADOCERA, GAMMARIDAE, HARPACTICOIDA, MYSIDAE, OTHER)
head(AugustSpeciesCorr)
AugSpeciesCorr.cor <- cor(AugustSpeciesCorr)
corrplot(AugSpeciesCorr.cor)


#For Environmental Variables
AugEnvCorr <- Augpcod.plot %>%
  select(NMDS1, NMDS2, NMDS3, TridentTemp,  SL_mm, HSI, LWResiduals_All)
AugEnvCorr_NAOmit <- na.omit(AugEnvCorr)
AugEnvCorr_NAOmit.cor = cor(AugEnvCorr_NAOmit)
corrplot(AugEnvCorr_NAOmit.cor)

cor(AugEnvCorr_NAOmit$SL_mm, AugEnvCorr_NAOmit$NMDS1)
cor(AugEnvCorr_NAOmit$LWResiduals, AugEnvCorr_NAOmit$NMDS1)
cor(AugEnvCorr_NAOmit$HSI, AugEnvCorr_NAOmit$NMDS1)
cor(AugEnvCorr_NAOmit$TridentTemp, AugEnvCorr_NAOmit$NMDS1)

##### August MRPP By Heatwave ####
mrpp(AugPreyData, AugEnvData$Heatwave, distance = 'bray', weight = 3)

##### August ISA By Heatwave ####
Aug_indval = multipatt(AugPreyData, AugEnvData$Heatwave, control = how(nperm=999))
summary(Aug_indval)



#### Manuscript Figure 3 ####
#No prey on contour plots
DietOrdinations_PSIRI <- July_prey + July_Size + July_Temp + Aug_prey + Aug_Size + Aug_Temp  +
  plot_layout(nrow = 2) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(DietOrdinations_PSIRI )
#dev.copy(jpeg,'Figure_3.jpg', width=16, height=10, units='in', res=300)
#dev.off()


#### Ordination Validation Figure (Supplemental Figure S9) ####
OrdinationTests <- JulyScree_Plot+  AugScree_Plot + JulyShepard_Plot + AugShepard_Plot + July_GoodnessofFit  + Aug_GoodnessofFit  +
  plot_layout(nrow = 3) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(OrdinationTests )
#dev.copy(jpeg,'Figure_S8.jpg', width=12, height=12, units='in', res=300)
#dev.off()

#### Ordinations with fish of similar Sizes ####

##### July Reduced Ordination ####

##Restrict July Size to between 50 and 60 mm 

head(JulyOrd)
identify_rows_July <- JulyOrd$SL_mm >= 50 & JulyOrd$SL_mm <= 60
July_SizeRestricted <- JulyOrd[identify_rows_July, ]           
July_SizeRestricted

aggregate(July_SizeRestricted, by = list(July_SizeRestricted$Heatwave),FUN = length)

#Create Environmental and Species Matrices for July 
JulyPreyData_Restricted <- July_SizeRestricted %>%
  select(CLADOCERA, MYSIDAE, GAMMARIDAE, HARPACTICOIDA, CAPRELLIDAE, CALANOIDA, OTHER)
head(JulyPreyData_Restricted)

JulyEnvData_Restricted <- July_SizeRestricted %>%
  select(Month, Year, Heatwave, SL_mm, HSI, LWResiduals_All, WholeBodyWW_g, Site, Bay, Habitat, TridentTemp)
head(JulyEnvData_Restricted)


## Run a Scree Plot using the package "goeveg"
#JulyScree_Restricted <- dimcheckMDS(JulyPreyData_Restricted ,distance = "bray",k = 6,trymax = 200,autotransform = FALSE)
#3 dimensions seems reasonable

## Run the July NMS Ordination with Restricted Sizes
JulyMDS_Restricted <-metaMDS(JulyPreyData_Restricted, distance="bray", k=3, trymax=200, autotransform=FALSE)

#Stress
JulyMDS_Restricted$stress

#Shepard Plot to evaluate Ordination Fit
stressplot(object = JulyMDS_Restricted, lwd = 5, main = "July Shepard Plot")

#Rotate to align with fish size
RotateJuly_Restricted <- MDSrotate(JulyMDS_Restricted, JulyEnvData_Restricted$SL_mm) 

## Set up your three NMS Principal Axes 
NMDS1 <- RotateJuly_Restricted$points[,1]
NMDS2 <- RotateJuly_Restricted$points[,2]
NMDS3 <- RotateJuly_Restricted$points[,3]

#Add them to the ordination data plot
Julypcod.plot_Restricted <-cbind(July_SizeRestricted, NMDS1, NMDS2, NMDS3)
head(Julypcod.plot_Restricted)


## Plot Restricted July NMS with Prey
species.scores1_July_Restricted <- as.data.frame(scores(RotateJuly_Restricted, "species"))
species.scores1_July_Restricted$species <- rownames(species.scores1_July_Restricted)
species.scores1_July_Restricted$z <- NA
head(species.scores1_July_Restricted)


July_prey_Restricted <- ggplot(data=Julypcod.plot_Restricted, aes(NMDS1, NMDS2))+
  geom_point(data=Julypcod.plot_Restricted, aes(NMDS1, NMDS2, color=Heatwave), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=Heatwave, color = Heatwave), alpha=.2,type='t',linewidth =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  geom_label(data = species.scores1_July_Restricted, aes(NMDS1, NMDS2, label = species), colour = "black", size = 3.5, position = position_jitter(0.45), alpha = 0.7) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("July Restricted Size (50 - 60 mm)")
print(July_prey_Restricted)


##### August Reduced Ordination ####

#Restrict August Size to between 70 and 90 mm 
identify_rows_Aug <- AugOrd$SL_mm >= 80 & AugOrd$SL_mm <= 90
August_SizeRestricted <- AugOrd[identify_rows_Aug, ]           
head(August_SizeRestricted)

aggregate(August_SizeRestricted, by = list(August_SizeRestricted$Heatwave),FUN = length)

#Create Environmental and Species Matrices for August 

AugPreyData_Restricted <- August_SizeRestricted %>%
  select(CALANOIDA, CAPRELLIDAE, CLADOCERA, GAMMARIDAE, HARPACTICOIDA, MYSIDAE, ANNELIDA, OTHER)
head(AugPreyData)

AugEnvData_Restricted <- August_SizeRestricted %>%
  select(Month, Year, Heatwave, SL_mm, HSI, LWResiduals_All, WholeBodyWW_g, Site, Bay, Habitat, TridentTemp)


## Run a Scree Plot using the package "goeveg"
#AugScree_Restricted <- dimcheckMDS(AugPreyData_Restricted ,distance = "bray",k = 6,trymax = 200,autotransform = FALSE)
#3 dimensions seems reasonable

## Run the Aug NMS Ordination with Restricted Sizes
AugMDS_Restricted <-metaMDS(AugPreyData_Restricted, distance="bray", k=3, trymax=200, autotransform=FALSE)

#Stress
AugMDS_Restricted$stress

#Shepard Plot to evaluate Ordination Fit
stressplot(object = AugMDS_Restricted, lwd = 5, main = "July Shepard Plot")

#Rotate to align with fish size
RotateAug_Restricted <- MDSrotate(AugMDS_Restricted, AugEnvData_Restricted$SL_mm) 

## Set up your three NMS Principal Axes 
NMDS1 <- RotateAug_Restricted$points[,1]
NMDS2 <- RotateAug_Restricted$points[,2]
NMDS3 <- RotateAug_Restricted$points[,3]

#Add them to the ordination data plot
Augpcod.plot_Restricted <-cbind(August_SizeRestricted, NMDS1, NMDS2, NMDS3)
head(Augpcod.plot_Restricted)


## Plot Restricted July NMS with Prey 
species.scores1_Aug_Restricted <- as.data.frame(scores(RotateAug_Restricted, "species"))
species.scores1_Aug_Restricted$species <- rownames(species.scores1_Aug_Restricted)
species.scores1_Aug_Restricted$z <- NA
head(species.scores1_Aug_Restricted)


Aug_prey_Restricted <- ggplot(data=Augpcod.plot_Restricted, aes(NMDS1, NMDS2))+
  geom_point(data=Augpcod.plot_Restricted, aes(NMDS1, NMDS2, color=Heatwave), show.legend=F, position=position_jitter(.1))+
  stat_ellipse(aes(fill=Heatwave, color = Heatwave), alpha=.2,type='t',linewidth =1, geom="polygon")+ 
  theme_classic()+
  scale_linetype_manual(values = "solid") +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_color_manual(values=c("#332288", "#888888", "#661100", "#000000")) +
  geom_label(data = species.scores1_Aug_Restricted, aes(NMDS1, NMDS2, label = species), colour = "black", size = 3.5, position = position_jitter(0.25), alpha = 0.7) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ylab("NMS2") +
  xlab("NMS1") +
  ggtitle("August Restricted Size (70 - 90 mm)")
print(Aug_prey_Restricted)


#### Restricted Figure for Manuscript (Figure S5) ####

SizeRestricted_DietOrdinations_PSIRI <- July_prey_Restricted + Aug_prey_Restricted +
  plot_layout(nrow = 1) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(SizeRestricted_DietOrdinations_PSIRI )
#dev.copy(jpeg,'Figure_S5.jpg', width=12, height=7, units='in', res=300)
#dev.off()

#### Restricted MRPP and ISA ####
mrpp(JulyPreyData_Restricted, JulyEnvData_Restricted$Heatwave, distance = 'bray', weight = 3)

mrpp(AugPreyData_Restricted, AugEnvData_Restricted$Heatwave, distance = 'bray', weight = 3)

July_indval_Restricted = multipatt(JulyPreyData_Restricted, JulyEnvData_Restricted$Heatwave, control = how(nperm=999))
summary(July_indval_Restricted)


Aug_indval_Restricted = multipatt(AugPreyData_Restricted, AugEnvData_Restricted$Heatwave, control = how(nperm=999))
summary(Aug_indval_Restricted)

