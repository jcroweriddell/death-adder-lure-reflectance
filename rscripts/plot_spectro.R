### Plot reflectance colour data for Acanthophis antarcticus tail skin 

### Install necessary packages if not already
if(! 'tidyverse' %in% rownames(installed.packages())){install.packages('tidyverse')}
if(! 'magrittr' %in% rownames(installed.packages())){install.packages('magrittr')}
if(! 'devtools' %in% rownames(installed.packages())){install.packages('devtools')}
if(! 'ggpomological' %in% rownames(installed.packages())){devtools::install_github('gadenbuie/ggpomological')}

### Load packages
library(tidyverse)
library(magrittr)
library(assertthat)
library(cowplot)
library(viridisLite)

# Palette for plots
pal <- c("Black1" = "#333333", 
         "Black2" = "#666666",
         "BlackVentral" = "#999999",  
         "White1" = "#CC9999",
         "White2" = "#FFCCCC",
         "YellowVentral" = "#FFCC99",
         "WhiteVentral" = "#CCCCFF",
         "SupralabialWhite" = "#FFCCCC",
         "SupralabialBlack" = "#333333",
         "PreoccularBrown" = "#CC9966")

### Function to plot spectro data
plotSpectro <- function(vec_samples = NULL, path_data, ext_files = NULL, n_commentLines = 14,
                        min_reflectance = 0, max_reflectance = 100, pal1 = pal,
                        save_plot = FALSE, path_plot = NULL, name_plot = NULL, image_type = 'pdf'){
  
  ## Checking inputs - error if conditions aren't met
  assert_that(dir.exists(path_data), msg = "Directory path to spectro data can't be resolved. Check the path exists")
  assert_that(!is.null(ext_files), msg = "File extension hasn't been provided")
  if(isTRUE(save_plot)){
    assert_that(!is.null(path_plot), msg = "Path to plot directory has not been set")
    assert_that(!is.null(name_plot), msg = "Name of plot file has not been set")
  }
  
  ## Read data
  rawSpectro <- list.files(path = path_data, pattern = ext_files, full.names = TRUE) %>%
    set_names(sub(ext_files, '', basename(.)))

  ## Subset samples for selection
  if(!is.null(vec_samples)){
    rawSpectro <- rawSpectro[str_detect(string = rawSpectro, pattern = vec_samples)]
  }

  ## Read in data
  rawSpectro <- rawSpectro %>%
    map(read_tsv,
        col_names = c("wavelength", "reflectance"),
        col_types = cols(),
        skip = n_commentLines) %>%
    bind_rows(.id = 'sample')

  ## Get mean reflectance for each wavelength value grouped by colour + sample
  meanSpectro <- rawSpectro %>%
    separate(col = sample, into = c("id", "colour", "replicate")) %>%
    group_by(id, colour, wavelength) %>%
    summarise(mean_reflectance = mean(reflectance))

  ## Plotting data
  
  # Plot
  p <- meanSpectro %>%
    filter(mean_reflectance > min_reflectance & mean_reflectance <= max_reflectance) %>%
    ggplot(aes(x = wavelength, y = mean_reflectance, colour = colour)) +
    geom_point(size = 1) + 
    scale_colour_manual(values = pal1, name = "Scale colour") +
    theme_bw() +
    labs(x = 'Wavelength', y = 'Mean Reflectance (%)') +
    xlim(300, 700) +
    ylim(0, 100) +
    geom_vline(xintercept = 400, linetype = "dotted", colour = "grey20", size = 0.8) +
    #annotate(geom="text", x = 350, y = 100, label = "UV", color = "grey20") +
    #geom_vline(xintercept = 740, linetype = "dotted", colour = "gray", size = 0.8) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          strip.placement = "inside",
          strip.background = element_blank(),
          strip.text.y = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=3)),
           text = element_text(size = 20)) +
    facet_wrap(~id)

  ## Save plot or just return it
  if(isTRUE(save_plot)){

    ## Create output directory if it doesn't exist already
    dir.create(path = normalizePath(path_plot), showWarnings = FALSE, recursive = TRUE)

    ## Build file name
    n <- paste(name_plot, image_type, sep = ".")

    ggsave(filename = n,
           device = image_type,
           path = normalizePath(path_plot),
           plot = p)
  } else {
    return(p)
  }
}

### Set directories
path_info <- "C:/Users/jmcr/Documents/death-adder-lure-reflectance/data/"
path_data <- "C:/Users/jmcr/Documents/death-adder-lure-reflectance/data/reflectance/Tail_measurements/"
path_data_h <- "C:/Users/jmcr/Documents/death-adder-lure-reflectance/data/reflectance/Head_measurements/"
path_plot <- "C:/Users/jmcr/Documents/death-adder-lure-reflectance/plots/ "

## Using function to plot the any sample
plotSpectro(vec_samples = c("JCR0019"), # select specimen numbers
            path_data = path_data,
            ext_files = '.txt' ,
            save_plot = FALSE, # TRUE = save plot, FALSE = display plot
            path_plot = paste0(path_plot, "uniform_tails"),
            name_plot = "JCR0019", 
            image_type = "png")


### Mean spectro across phenotypes

## Load data
rawSpectro <- list.files(path = path_data, pattern = ".txt", full.names = TRUE) %>%
  set_names(sub(".txt", '', basename(.))) %>%
  map(read_tsv,
      col_names = c("wavelength", "reflectance"),
      col_types = cols(),
      skip = 14) %>%
  bind_rows(.id = 'sample')


## Get mean reflectance for each wavelength value grouped by colour + sample
meanSpectro <- rawSpectro %>%
  separate(col = sample, into = c("id", "colour", "reflection", "replicate")) %>%
  group_by(id, colour, wavelength) %>%
  summarise(mean_reflectance = mean(reflectance))

# tail phenotype info
sample_info <- read_csv(paste0(path_info, "sample_info.csv"), col_names = TRUE) %>%
  rename(id = `JCR ID`, phenotype = tail_phenotype) %>%
  select(id, phenotype, sex_m_f, svl_mm, mass_g)

# join spectro data with tail phenotype info
meanSpectro <- full_join(meanSpectro, sample_info)

##
pheno <- sample_info %>% 
  group_by(phenotype) %>% 
  arrange(phenotype, svl_mm, mass_g) 
View(pheno)
##

## Mean reflectance for each colour and phenotype
pheno_meanSpectro <- meanSpectro %>%
  mutate(colour = gsub(x = colour, pattern = "1|2", replacement = "")) %>% 
  group_by(phenotype, colour, wavelength) %>%
  summarise(meanSum_reflectance = mean(mean_reflectance), 
            min_reflectance = min(mean_reflectance),
            max_reflectance = max(mean_reflectance),
            sd_reflectance = sd(mean_reflectance))

## Palette
pal2 <- c("Black" = "#333333", 
          "BlackVentral" = "#999999",  
          "White" = "#CC9999",
          "YellowVentral" = "#FFCC99",
          "WhiteVentral" = "#CCCCFF")

## filter by phenotype
spectro_plot <- pheno_meanSpectro %>%
  filter(phenotype %in% "Banded") %>% # change text to modify phenotype: "Banded", "Transition", "Striped", "Uniform"
  ggplot(aes(x = wavelength, y = meanSum_reflectance, colour = colour)) +
  geom_line(size = 5) +
  geom_errorbar(aes(ymin = meanSum_reflectance - sd_reflectance, ymax = meanSum_reflectance + sd_reflectance), 
                width=.5,position=position_dodge(.9), alpha = 0.1)  +
  scale_colour_manual(values = pal2) +
  labs(x = '', y = '') +
  theme_classic() +
  xlim(300, 700) +
  ylim(0, 100) +
  geom_vline(xintercept = 380, linetype = "dotted", colour = "gray", size = 1.5) +
  #annotate(geom="text", x = 350, y = 100, label = "UV", color = "grey20") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(colour = 'black', size = 2),
        axis.ticks.x = element_line(colour = "black", size = 2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.ticks.y = element_line(colour = "black", size = 2),
        legend.position = "none",
        text = element_text(size = 20)) 
spectro_plot


## Save plots
ggsave(filename = "banded_mean.png",
       device = "png",
       path = normalizePath(path_plot),
       plot = spectro_plot,
       width = 297, 
       height = 210, 
       units = "mm")
##

## Plot reflectance for 'yellow ventral' measurement only
yellow <- meanSpectro %>% 
  filter(grepl("Yellow", colour)) %>%
  ggplot(aes(x = wavelength, y = mean_reflectance, colour = id)) +
  geom_point(size = 0.5) +
  theme_bw() +
  labs(x = 'Wavelength', y = 'Mean Reflectance (%)') +
  xlim(300, 700) +
  ylim(0, 100) +
  scale_color_viridis_d(alpha=0.6, direction = -1, name = "Individual") +
  geom_vline(xintercept = 400, linetype = "dotted", colour = "grey20", size = 0.8) +
  annotate(geom="text", x = 350, y = 100, label = "UV", color = "grey20") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.placement = "inside",
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=3)),
         text = element_text(size = 20))

## Save plot
ggsave(filename = "yellow.png",
       device = "png",
       path = normalizePath(path_plot),
       plot = yellow,
       width = 297, 
       height = 210, 
       units = "mm")
##

## Plot head measurements

plotSpectro(vec_samples = c("JCR0018"), # select specimen numbers
            path_data = path_data_h,
            ext_files = '.txt' ,
            save_plot = FALSE, # TRUE = save plot, FALSE = display plot
            path_plot = path_plot,
            name_plot = "head_measurements", 
            image_type = "png")
##