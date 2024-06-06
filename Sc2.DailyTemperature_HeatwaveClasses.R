#### Trident Bay, Kodiak Island Temperatures and Heatwave Classifications
### Date Created: 1/11/2024
# R v. 4.2.2

# Load libraries 
library(rstudioapi) #Load Working Directory
library(tidyverse) #For data wrangling and visualization
library(lubridate) #For working with dates
library(heatwaveR) #Heatwave Classifications
library(patchwork) #Data visualization
library(viridis) #Color Palette


#Paired with:
#Interpolated Trident Bay Data.  

##Set Working Directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#Sets directory to folder that this R script is saved to

#### Daily Temperature from Trident Bay ####
Trident <- read.csv("D4.TridentBay_DailyTemperature.csv") %>%
  mutate(Year = as.factor(Year)) %>%
  mutate(Date = as.Date(Date, "%Y-%m-%d")) %>%
  mutate(DATE_noyr = format(Date, "%b-%d")) %>%
  mutate(JulianDate = as.numeric(strftime(Date, format = "%j"))) 

#Subset Trident Data to years of cod collections 
TridentTrim <- subset(Trident,
                      Date >= '2006-01-01' &
                        Date <= '2019-12-31') %>%
  mutate(Heatwave = case_when(
    Year == '2006'| Year == '2007'|Year == '2008'|Year == '2009' | Year == '2010' | Year == '2011'| Year == '2012'| Year == '2013'~ "Before",
    Year == '2014' | Year == '2015' | Year =='2016'| Year == '2019' ~ "Heatwave",
    Year == "2017" | Year == "2018" ~ "Between"
  ))


##### Summary of Trident Temperatures (Table 1, Table S1) ####

Temp_summary <-  TridentTrim %>%
  filter(Month >= '7') %>%
  filter(Month <= '8') %>%
  group_by(Month,Year) %>%
  reframe(Temp_mean = mean(TridentTemp, na.rm = T), Temp_se = (sd(TridentTemp, na.rm = T)/sqrt(length(TridentTemp))))

Temp_summary_HW <-  TridentTrim %>%
  filter(Month >= '7') %>%
  filter(Month <= '8') %>%
  group_by(Month,Heatwave) %>%
  reframe(Temp_mean = mean(TridentTemp, na.rm = T), Temp_se = (sd(TridentTemp, na.rm = T)/sqrt(length(TridentTemp))))

##### Plot Temperature (Figure S1) ####
TridentTemps_byYear <- ggplot(data = TridentTrim, aes(as.Date(DATE_noyr, format = "%b-%d"), TridentTemp, color = Year)) + 
  geom_smooth(linewidth = 1.5) +
  scale_color_viridis_d(option="turbo") +
  xlab("Date") + 
  ylab("Trident Bay Temperature (\u00b0C)") +
  theme_classic(base_size = 17, base_family = "")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +   
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) 
print(TridentTemps_byYear)



TridentTemps_byHW <- ggplot(data = TridentTrim, aes(as.Date(DATE_noyr, format = "%b-%d"), TridentTemp, color = Heatwave)) +
  geom_smooth(linewidth = 1.5) +
  xlab("Date") + 
  ylab("Trident Bay Temperature (\u00b0C)") +
  theme_classic(base_size = 17, base_family = "")+
  scale_color_manual(values=c("#402aaa", "grey1", "#CC6677")) +
  scale_fill_manual(values=c("#332288", "#888888", "#CC6677")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=15, angle = 45, vjust = 0.5, hjust=.5), axis.text.y = element_text(size = 15),  axis.title=element_text(size=14,face="bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) 
print(TridentTemps_byHW)


TridentTemps <-  TridentTemps_byHW + TridentTemps_byYear +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold', size = 16), legend.position = "bottom")
print(TridentTemps)
#dev.copy(jpeg,'Figure_S1.jpg', width=12, height=7, units='in', res=300)
#dev.off()


#### Marine Heatwave Classifications for Trident Bay (Table S2) ####

TridentHW <- Trident %>%
  mutate(t = Date) %>%
  mutate(temp = TridentTemp) %>%
  select(t, temp)

#### Detect Heatwaves in the Trident Data
ts <- ts2clm(TridentHW, climatologyPeriod = c("1983-01-01", "2012-12-31"))
mhw <- detect_event(ts) 

mhw_cat <- category(mhw, S = FALSE) %>%
  mutate(date_peak = peak_date)

#MHW Classifications (Table S1)
TridentHW_Events <- mhw$event %>% 
  ungroup() %>%
  left_join(mhw_cat) %>%
  select(category, season, duration, date_start, date_peak, date_end,intensity_cumulative, intensity_mean, intensity_max) %>%
  filter(date_start >= '2006-01-01') %>%
  filter(date_start <='2019-12-31')
#write.csv(TridentHW_Events, "Table_S1.csv")


#Trident HW Anomaly Plot
event_line(mhw, spread = 4000, metric = "intensity_cumulative", start_date = "2006-01-01", end_date = "2020-12-31")

#Trident HW Category Plot
event_line(mhw, spread = 4000, metric = "intensity_cumulative", start_date = "2006-01-01", end_date = "2020-12-31", category = TRUE)

#Trident Lolli-plot
lolli_plot(mhw, metric = "intensity_cumulative")


