### Plotting diet data for Acanthophis antarcticus death adders

### Packages
library(tidyverse)
library(magrittr)
library(RColorBrewer)

### Set directories
dir_data <- "~/death-adder-lure-reflectance/data/diet/"
dir_plots <- "~/death-adder-lure-reflectance/plots/"

### Load diet data
prey <- read_csv("prey.csv", col_names = TRUE) %>%
  rename(SVL_cm = `SVL (cm)`, Prey_activity = `Prey activity`, Catalog = `catalog#`) %>%
  select(Catalog, species, State, Country, SVL_cm, Endotherm, Ectotherm, Prey_taxa_major, Prey_taxa_family, Prey_activity) %>%
  mutate(Prey_type = case_when(Endotherm == 1 ~ 1,
                               Ectotherm == 1 ~ -1)) %>%
  mutate(SVL_mm = SVL_cm*10) %>%
  mutate(SVL_cm = round(SVL_cm)) %>%
  drop_na(SVL_cm) %>%
  filter(species == "ant") %>% # filter for only A. antarcticus
  filter(Country == "AUS") # filter for only specimens collected in Australia

### Plot prey items recorded ~ svl of specimens

# Filter diet by endothermic/ectothermic
endo <- prey %>% filter(Endotherm == 1)
ecto <- prey %>% filter(Ectotherm == 1)

# Bar plot
pDiet <- ggplot() +
  geom_bar(data = endo, aes(x = SVL_cm, y = Prey_type, fill = Prey_taxa_major), 
           colour = "black", stat = "identity", width = 1) +
  geom_bar(data = ecto, aes(x = SVL_cm, y = Prey_type, fill = Prey_taxa_major), 
           colour = "black", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(15, 75, by = 10)) +
  scale_fill_brewer(palette = "Set2") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
pDiet

## Save plot
ggsave(filename = "plot_diet.pdf",
       device = "pdf",
       path = normalizePath(dir_plots),
       plot = pDiet,
       width = 500, 
       height = 180, 
       units = "mm")

###
sessionInfo()
###

