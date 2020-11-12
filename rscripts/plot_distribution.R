### Plotting diet data for Acanthophis antarcticus death adders

### Packages
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(cowplot)

### Set directories
setwd("C:/Users/jmcr/Documents/death-adder-lure-reflectance/")

# Load coordinates for A. antarcticus
spp_records <- read_csv(file = "data/distribution/Acanthophis_records_2020.csv", col_names = TRUE)

# Filter for a single species, remove outliers that are different species
Acan <- spp_records %>%
  filter(Species %in% "Acanthophis antarcticus", `State - parsed` != "Northern Territory",
         `IBRA 7 Regions` != "Mount Isa Inlier", `IBRA 7 Regions` != "Mitchell Grass Downs",
         Latitude <= -14) %>%
  select(Species, Latitude, Longitude, State = `State - parsed`) %>%
  drop_na(Latitude, Longitude) %>%
  mutate(Species_abbrev = gsub(x = Species, pattern = "canthophis", replacement = ".")) %>%
  drop_na(Latitude, Longitude, Species) %>%
  select(Species_abbrev, Latitude, Longitude)

## Try another  method
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")

## Plot study site population

# Smoky Bay co-ordinates
SB <- data.frame(Latitude = -32.3667, Longitude = 133.9333, Region = "Smoky Bay")
B <- data.frame(Latitude = -30.4520343, Longitude = 152.8966816, Region = "Belingen")

inset <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = Acan, aes(Longitude, Latitude), size = 2, colour = "black" ) +
  geom_point(data = SB, aes(Longitude, Latitude), size = 8, shape = 15, alpha = .90, colour = "#99CCFF") +
  coord_sf(xlim = c(129, 140), ylim = c(-38, -30.5), expand = FALSE)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

## Plot Australian distribution
Aus <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = Acan, aes(Longitude, Latitude), colour = "black") +
  geom_point(data = SB, aes(Longitude, Latitude), size = 6, shape = 15, alpha = .90, colour = "#99CCFF") +
  geom_point(data = B, aes(Longitude, Latitude), size = 6, shape = 15, alpha = .90, colour = "#FF9933") +
  coord_sf(xlim = c(110, 160), ylim = c(-45, -10), expand = FALSE) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

Aus

Acan_map <- ggdraw(Aus) +
  draw_plot(inset, x = 0.08, y = 0.6, 0.40, 0.4)

# Save plot
ggsave(filename = "plots/Acan_map.png",
       plot = Acan_map,
       device = "png",
       width = 297, 
       height = 210, 
       units = "mm")
##