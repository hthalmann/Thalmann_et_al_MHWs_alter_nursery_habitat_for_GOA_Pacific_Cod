### Kodiak Island Juvenile Pacific Cod Sampling Map

#Load Libraries
library(tidyverse)  #data wrangling and visualization
library(PBSmapping) #Use nepacLLhigh dataset for Alaska coastline
library(cowplot) # Data Visualization

#Set Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in

# Zoomed in map of Kodiak with sites

data('nepacLLhigh') 
nepacLLhigh %>% 
  dplyr::select(group=PID, POS=POS,long=X, lat=Y) -> ak 

Sites_Zoom <- ggplot() + 
  geom_polygon(data = ak, aes(long, lat, group = group), fill=8, color='black') +
  theme_bw() +
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'W')))+
  coord_map(projection = "albers", lat0 = 35, lat1 = 72,
            xlim = c(-152.69, -152.22), ylim = c(57.75, 57.9)) +
  annotate(geom = "point", x=-152.633, y=57.87, size = 8, shape = 24, color="black", fill = "black") +
  annotate(geom = "point", x=-152.633, y=57.87, size = 6, shape = 24, color="black", fill = "dodgerblue") +
  annotate(geom = "point", x=-152.287, y=57.779, size = 8, shape = 24, color="black", fill = "black") +
  annotate(geom = "point", x=-152.287, y=57.779, size = 6, shape = 24, color="black", fill = "dodgerblue") +
  annotate(geom = "point", x=-152.392, y=57.778, size = 8, shape = 24, color="black", fill = "black") +
  annotate(geom = "point", x=-152.392, y=57.778, size = 6, shape = 24, color="black", fill = "dodgerblue") +
  annotate(geom = "label", x=-152.54, y=57.872, label="Anton Larsen Bay", alpha = 0.8,size =4,  label.size = 1, color="black") +
  annotate(geom = "label", x=-152.29, y=57.767, label="Cook Bay", alpha = 0.8, label.size = 1, size = 4, color="black") +
  annotate(geom = "label", x=-152.39, y=57.767, label="Trident Bay", alpha = 0.8, label.size = 1, size = 4, color="black") +
  
  annotate(geom = "point", x=-152.4, y=57.798, size = 5, shape = 21, color="black", fill = "black") +
  annotate(geom = "label", x=-152.45, y=57.798, label="Kodiak", alpha = 0.8, label.size = 1, size = 4, color="black") +
  annotate(geom = "label", x=-152.6, y=57.8, label="Kodiak\nIsland", alpha = 0, label.size = NA, size = 6, color="black", fontface = 2) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14),
        legend.text = element_text(size = 6), legend.title = element_blank()) +
  theme(panel.border = element_rect(colour = "red", fill=NA, size=4)) 

print(Sites_Zoom)


#Larger Gulf map with inset indicated

data('nepacLLhigh') 
nepacLLhigh %>% 
  dplyr::select(group=PID, POS=POS,long=X, lat=Y) -> CanadaAK 


Inset_AK <- ggplot() + 
  geom_polygon(data = CanadaAK, aes(long, lat, group = group), fill=8, color='black') +
  coord_map(projection = "albers", lat0 = 35, lat1 = 72,
            xlim = c(-158, -140), ylim = c(54, 61)) +
  theme_bw() +
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'W')))+
  annotate("rect", xmin = -152.22, xmax = -152.69, ymin = 57.75, ymax = 57.9, alpha = 0, size = 3, color = "black", fill = NA) +
  annotate("rect", xmin = -152.22, xmax = -152.69, ymin = 57.75, ymax = 57.9, alpha = 0, size = 0.5, color = "red", fill = NA) +
  annotate(geom = "label", x=-157.0, y=60, label="Alaska, USA", alpha = 0, label.size = NA, size = 7, color="black", fontface = 2) +
  annotate(geom = "label", x=-155.5, y=55, label="Gulf of Alaska", alpha = 0.8, label.size = NA, size = 7, color="black", fontface = 2) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 9), legend.title = element_blank())
print(Inset_AK)

# Sampling Map (Figure S7)
plot.with.inset <-
  ggdraw() +
  draw_plot(Inset_AK) +
  draw_plot(Sites_Zoom, x = 0.407, y = -0.05, width = .52, height = .8)

plot.with.inset
#Note, if trouble plotting, increase size of plotting window

#ggsave("Figure_S7.jpg", plot = plot.with.inset, device = "jpeg",dpi = 600)


