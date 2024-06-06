#### Kodiak Island Juvenile Pacific Cod, Otolith Growth Back Calculations using a Biological Intercept 

## Date Created: 1/11/24
## R version: 4.2.2

## Load Libraries 
library(rstudioapi) #load R working directory
library(tidyverse)  #data wrangling and visualization
library(patchwork) #Data Visualization
library(ggpubr) #For scatter plots
library(viridis) #Color Palette


#Pairs with Final otolith increment measurements in wide format

##Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#Sets directory to folder that this R script is saved to

#### Read in files

wide_increments <- read.csv("D11.OtolithIncrements_Raw.csv")

### Pivot Long
AllData_Long <- wide_increments %>% pivot_longer(cols = c("increment21_edge", "increment20", "increment19", "increment18", "increment17", "increment16", "increment15", "increment14", "increment13", "increment12", "increment11", "increment10", "increment9", "increment8", "increment7", "increment6", "increment5", "increment4", "increment3", "increment2", "increment1_closertocore"), names_to = "Increment", values_to = "IncrementWidth")


AllData_Long <- AllData_Long %>%
  mutate(IncrementNumber = case_when(
   Increment == 'increment1_closertocore' ~ 1,
   Increment == 'increment2' ~ 2,
   Increment == 'increment3' ~ 3,
   Increment == 'increment4' ~ 4,
   Increment == 'increment5' ~ 5,
   Increment == 'increment6' ~ 6,
   Increment == 'increment7' ~ 7,
   Increment == 'increment8' ~ 8,
   Increment == 'increment9' ~ 9,
   Increment == 'increment10' ~ 10,
   Increment == 'increment11' ~ 11,
   Increment == 'increment12' ~ 12,
   Increment == 'increment13' ~ 13,
   Increment == 'increment14' ~ 14,
   Increment == 'increment15' ~ 15,
   Increment == 'increment16' ~ 16,
   Increment == 'increment17' ~ 17,
   Increment == 'increment18' ~ 18,
   Increment == 'increment19' ~ 19,
   Increment == 'increment20' ~ 20,
   Increment == 'increment21_edge' ~ 21
  ))
head(AllData_Long)

#### Calculate Radius at age for each of the 21 days ####
AllData_Long_Ri <- AllData_Long %>%
  group_by(FISHID) %>%
  mutate(cum_IW = cumsum(IncrementWidth)) %>%
  mutate(rad_at_age = OtoRadius_Avg - cum_IW)


#### Relationship between Standard Length (at capture) and Otolith Radius (at capture) ####

##Year as Factor and log transformed SL and Ri 

AllData_Long_Ri$Year <- factor(AllData_Long_Ri$Year, levels = c("2006", "2007", "2008", "2009", "2010", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))

AllData_Long_Ri <- AllData_Long_Ri %>%
  mutate(Year = as.factor(Year)) %>%
  mutate(log_sl_capture = log(SL_mm))%>%
  mutate(log_rad_at_age = log(rad_at_age)) %>%
  mutate(log_rad_capture  = log(OtoRadius_Avg))


July_Increments <- AllData_Long_Ri[AllData_Long_Ri$Month == 7,]
August_Increments <- AllData_Long_Ri[AllData_Long_Ri$Month == 8,]


## Plot log transform of standard length against otolith radius (Figure S10) ##

log_July_OR_SL <- ggscatter(July_Increments, x = "log_rad_capture", y = "log_sl_capture", 
                            color = "Year", add = "reg.line") +
  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"), color = Year)) +
  scale_color_viridis_d(option="turbo") +
  xlab(expression(paste(" log (Otolith Radius (", mu, "m))"))) +
  ylab("log (Standard Length (mm))") + 
  ggtitle("July") +
  theme(legend.position="right") +
  theme(axis.text.x=element_text(size=15, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title.y =element_text(size=14),  axis.title.x =element_text(size=14), legend.text = element_text(size = 12), legend.title = element_text(size = 14),   plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
print(log_July_OR_SL)

log_Aug_OR_SL <- ggscatter(August_Increments, x = "log_rad_capture", y = "log_sl_capture", 
                           color = "Year", add = "reg.line") +
  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"), color = Year)) +
  scale_color_viridis_d(option="turbo") +
  xlab(expression(paste(" log (Otolith Radius (", mu, "m))"))) +
  ylab("log (Standard Length (mm))") + 
  ggtitle("August") +
  theme(legend.position="right") +
  theme(axis.text.x=element_text(size=15, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title.y =element_text(size=14),  axis.title.x =element_text(size=14), legend.text = element_text(size = 12), legend.title = element_text(size = 14),   plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
print(log_Aug_OR_SL)

OR_SL_Relationship <-  log_July_OR_SL + log_Aug_OR_SL +
  plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(OR_SL_Relationship)
#dev.copy(jpeg,'Figure_S10.jpg', width=10, height=14, units='in', res=300)
#dev.off()


#### Biological Intercept Calculation ####
# Equation from Campana 1990:

# L(agex) = L(capture) + ((O(age) - O(capture))(L(capture) - L(age0))) / (O(capture) - O(age0))


#Biological intercepts from Narimatsu et al. 2007

#Fish TL at hatching = 4.1 mm - converted to SL below 
#Otolith nucleus check = 8.3 um


L0p.TL <- 4.1
R0 <- 8.3
L0 <- 0.952*L0p.TL

Lcpt <- AllData_Long_Ri$SL_mm
Rcpt <- AllData_Long_Ri$OtoRadius_Avg
Ri <- AllData_Long_Ri$rad_at_age

BIC_function <- function(L0, Lcpt, R0, Rcpt, Ri){
  result <- Lcpt + ((Ri - Rcpt)*(Lcpt - L0)) / (Rcpt - R0)
  return(result)
}

AllData_Long_Ri$Size_Initial_BIC <- BIC_function(L0, Lcpt, R0, Rcpt, Ri)
head(AllData_Long_Ri)

#### Calculate Absolute Growth Per Day 
AllData_Growth <- AllData_Long_Ri %>%
  group_by(FISHID) %>%
  mutate(totalgrowth_BIC = SL_mm - Size_Initial_BIC) %>%
  mutate(mm_day_short_BIC = c(NA , diff(-Size_Initial_BIC))) %>%
  mutate(mm_day_BIC = coalesce(mm_day_short_BIC, totalgrowth_BIC))

#### Calculate Relative Growth Per Day

BICGrowth <- AllData_Growth %>%
  mutate(mm_mm_day_BIC =mm_day_BIC /Size_Initial_BIC) %>%
  select(FISHID, Year, Month, Heatwave, SL_mm, OtoRadius_Avg, IncrementNumber, IncrementWidth, rad_at_age, log_sl_capture, log_rad_at_age, log_rad_capture, Size_Initial_BIC, mm_day_BIC, mm_mm_day_BIC)

#Write csv of back-calculated growth for subsequent analysis
#write.csv(BICGrowth, "BackCaculatedGrowth_ROutput.csv")

#### Subset July and August

AugDailyGrowth <- BICGrowth[BICGrowth$Month == 8,]
JulyDailyGrowth <-BICGrowth[BICGrowth$Month == 7,]

#### Plot Final 21 days of Relative Growth (Figure S11 ) ####

## Plot Relative Growth by year
July_RelativeGrowth_BIC <- ggplot(JulyDailyGrowth, aes(x=IncrementNumber, y = mm_mm_day_BIC)) +
  geom_smooth(aes(color = Year, fill = Year))+
  scale_color_viridis_d(option="turbo") +
  scale_fill_viridis_d(option="turbo") +
  coord_cartesian(ylim = c(0.005,0.025)) +
  theme_classic(base_size = 17, base_family = "")+
  theme(legend.position = "right")+
  labs(x = "Increment Number", y = "Relative Growth (mm/mm/day)")+
  ggtitle("July") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(July_RelativeGrowth_BIC)

Aug_RelativeGrowth_BIC <- ggplot(AugDailyGrowth, aes(x=IncrementNumber, y = mm_mm_day_BIC)) +
  geom_smooth(aes(color = Year, fill = Year))+
  scale_color_viridis_d(option="turbo") +
  scale_fill_viridis_d(option="turbo") +
  coord_cartesian(ylim = c(0.005,0.025)) +
  theme_classic(base_size = 17, base_family = "")+
  theme(legend.position = "right")+
  labs(x = "Increment Number", y = "Relative Growth (mm/mm/day)")+
  ggtitle("August") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(Aug_RelativeGrowth_BIC)

#Only plot years with both months included
AugDailyGrowth_Trim <- AugDailyGrowth %>%
  filter(!Year %in% 2006) %>%
  filter(!Year %in% 2008) %>%
  filter(!Year %in% 2016)

JulyDailyGrowth_Trim <- JulyDailyGrowth %>%
  filter(!Year %in% 2006) %>%
  filter(!Year %in% 2008) %>%
  filter(!Year %in% 2016)

#Plot years with trimmed intraannual data
July_RelativeGrowth_BIC_Trim <- ggplot(JulyDailyGrowth_Trim, aes(x=IncrementNumber, y = mm_mm_day_BIC)) +
  geom_smooth(aes(color = Year, fill = Year))+
  scale_color_viridis_d(option="turbo") +
  scale_fill_viridis_d(option="turbo") +
  coord_cartesian(ylim = c(0.005,0.025)) +
  theme_classic(base_size = 17, base_family = "")+
  theme(legend.position = "right")+
  labs(x = "Increment Number", y = "Relative Growth (mm/mm/day)")+
  ggtitle("July") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(July_RelativeGrowth_BIC_Trim)

Aug_RelativeGrowth_BIC_Trim <- ggplot(AugDailyGrowth_Trim, aes(x=IncrementNumber, y = mm_mm_day_BIC)) +
  geom_smooth(aes(color = Year, fill = Year))+
  coord_cartesian(ylim = c(0.005,0.025)) +
  scale_color_viridis_d(option="turbo") +
  scale_fill_viridis_d(option="turbo") +
  theme_classic(base_size = 17, base_family = "")+
  theme(legend.position = "right")+
  labs(x = "Increment Number", y = "Relative Growth (mm/mm/day)")+
  ggtitle("August") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(Aug_RelativeGrowth_BIC_Trim)


##Plot Relative Growth by HW, BIC
July_RelativeGrowth_HW_BIC <- ggplot(JulyDailyGrowth, aes(x=IncrementNumber, y = mm_mm_day_BIC)) +
  geom_smooth(aes(color = Heatwave, fill = Heatwave))+
  scale_color_manual(values=c("#402aaa", "grey1", "#CC6677")) +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  theme_classic(base_size = 17, base_family = "")+
  coord_cartesian(ylim = c(0.006,0.018)) +
  theme(legend.position = "right")+
  labs(x = "Increment Number", y = "Relative Growth (mm/mm/day)")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(July_RelativeGrowth_HW_BIC)

Aug_RelativeGrowth_HW_BIC <- ggplot(AugDailyGrowth, aes(x=IncrementNumber, y = mm_mm_day_BIC)) +
  geom_smooth( aes(color = Heatwave, fill = Heatwave))+
  scale_color_manual(values=c("#402aaa", "grey1", "#CC6677")) +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  theme_classic(base_size = 17, base_family = "")+
  coord_cartesian(ylim = c(0.006,0.018)) +
  theme(legend.position = "right")+
  labs(x = "Increment Number", y = "Relative Growth (mm/mm/day)")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
print(Aug_RelativeGrowth_HW_BIC)


Backcalculated_RelativeGrowth_Plots <- July_RelativeGrowth_BIC_Trim  + Aug_RelativeGrowth_BIC_Trim + July_RelativeGrowth_HW_BIC + Aug_RelativeGrowth_HW_BIC +
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16))
print(Backcalculated_RelativeGrowth_Plots)
#dev.copy(jpeg,'Figure_S11.jpg', width=14, height=10, units='in', res=300)
#dev.off()

