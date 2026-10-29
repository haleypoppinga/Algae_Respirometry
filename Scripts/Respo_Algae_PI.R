###### Respo Code for PI Curves ####### 
### Created by: Nyssa Silbiger, Maya Powell, Haley Poppinga
#### Last updated on: 2026-09-27

############## Introduction to code/script ####################
# PI (PPFD) Curves using Platt 1980 Photoinhibition Model

### Install Packages #####
## if these packages are not yet installed, install them 
## great for updates or new users 
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented')
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix')
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra')
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') devtools::install_github('colin-olito/LoLinR')
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron')
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse')
if ("here" %in% rownames(installed.packages()) == 'FALSE') install.packages('here')
if ("patchwork" %in% rownames(installed.packages()) == 'FALSE') install.packages('patchwork')
if ("PNWColors" %in% rownames(installed.packages()) == 'FALSE') install.packages('PNWColors')
if ("ggtext" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggtext')
if ("ggforce" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggforce')

#Read in required libraries
##### Include Versions of libraries
library(segmented)
library(plotrix)
library(gridExtra)
library(LoLinR)
library(lubridate)
library(chron)
library(patchwork)
library(tidyverse)
library(here)
library(PNWColors)
library(ggrepel)
library(reshape2)
library(viridis)
library(car)
library(future)
library(furrr)
library(dplyr)
library(beepr)
library(nls.multstart)
library(ggpubr)
library(ggforce)


############################################################
### FINAL PI ANALYSIS ###
############################################################


### Read in data (AFDW-normalized from Respo_Algae_PI_Models R script)
######### Load processed data and metadata ###############
RespoR <- read_csv(here("Data","Respo_Files","PI","Respo_Algae_R.csv"))
RespoR_Normalized <- read_csv(here("Data","Respo_Files","PI","Respo_Algae_RNormalized_AllPIRates.csv"))


## Species polette
# set up species colors
species_palette <- c(av  = "#033500",  # green algae division
                 hd  = "#66A61E",  # green
                 cs  = "#20B2AA",  # green
                 dc  = "#239665",  # green
                 as  = "#E41A1C",  # red
                 gs  = "#FFA500",  # red
                 sf  = "#7C0A02",  # red
                 da  = "#A57A54",  # brown
                 ds  = "#A89F51")  # brown

# full species names df
sp_names <- read_csv(here("Data/species_names.csv"))
# join species names with the normalized data
RespoR_Normalized <- RespoR_Normalized %>% left_join(sp_names, by = "species")

# full_species list:
# av, as, gs, da, ds, dc, sf, cs, hd


### Define final Platt photoinhibition model
# Platt et al. (1980) Photoinhibition Model ###
platt_model <- function(light_value, Ps, AQY, beta, Rd) {
  Ps *(1 - exp(-(AQY*light_value / Ps))) *exp(-(beta * light_value / Ps)) - Rd
}

##### MODEL PARAMETERS ######
# Ps = photosynthetic scaling parameter, maximum light-saturated photosynthetic rate 
# Pmax = fitted maximum photosynthetic rate (aka Pm), actual peak photosynthetic rate
# AQY = alpha = initial low-light slope
# beta = photoinhibition parameter
# beta = 0 little/no photoinhibition
# larger beta = stronger photoinhibition (high light decline)
# Rd = dark respiration

# Calculate eventually:
# Ik = saturation 
# Ik = Ps/alpha = Pmax / AQY
# Iopt = actual fitted irradiance where photosynthesis peaks, PPFD where the curve peaks before photoinhibition decline
# if beta = 0, no Iopt b/c no decline, Iopt would be beyond range/approaching infinity


# lists to save outputs from each individual Platt Model fit
coef_list <- list() # fitted parameters for each individual
curve_list <- list() # smooth fitted curve for each 
raw_list <- list() # observed raw data for each


# make forloop
for (i in unique(RespoR_Normalized$algae_id)) {
  # subset one individual's PI curve
  data <- RespoR_Normalized |>
    filter(algae_id == i) |>
    arrange(light_value)
  
  # identify species and zone 
  sp <- data$species[1] 
  sp_name <- data$species_name[1] # full species name
  z <- data$zone[1] # zone
  cat("\nModel 2:", sp, z, i, "\n") # catenate
  # make forloop
  for (i in unique(RespoR_Normalized$algae_id)) {
    # subset one individual's PI curve
    data <- RespoR_Normalized |>
      filter(algae_id == i) |>
      arrange(light_value)
    
    # identify species and zone 
    sp <- data$species[1] 
    sp_name <- data$species_name[1] # full species name
    z <- data$zone[1] # zone
    cat("\nModel 2:", sp, z, i, "\n") # catenate
    
    # use nls_multstart to fit the photoinhibition curve
    # nonlinear models where we dont know best starting parameter values
    curve.platt <- tryCatch(nls_multstart(
      umol.g.hr ~ platt_model(light_value, Ps, AQY, beta, Rd),
      data = data,
      iter = 500, # try 500 different starting values within the range below
      start_lower = c(Ps = max(data$umol.g.hr, na.rm = TRUE) * 0.5, # ranges of starting values
                      AQY = 0.001, beta = 0.0001, Rd = abs(min(data$umol.g.hr, na.rm = TRUE)) * 0.5),
      start_upper = c(Ps = max(data$umol.g.hr, na.rm = TRUE) * 3, AQY = 0.5, beta = 0.5,
                      Rd = abs(min(data$umol.g.hr, na.rm = TRUE)) * 2),
      lower = c(Ps = 0.000001, # cannot be 0 b/c light value/Ps cannot be divided by 0
                AQY = 0.000001, beta = 0,Rd = 0),
      supp_errors = "Y"),
      
      # if one individual fails, continue to the next one
      error = function(e) {
        message("MODEL 2 FAILED FOR ", i, ": ", e$message)
        return(NULL)
      }
    ) 
    # skip to next individual if Model 2 failed
    if(is.null(curve.platt)) {
      next
    }
    
    # extract fitted coefficients
    coefs <- coef(curve.platt)
    
    # calculate Pmax value from Platt model
    Pmax_val <- coefs["Ps"]* (coefs["AQY"] / (coefs["AQY"] + coefs["beta"])) *
      (coefs["beta"] / (coefs["AQY"] + coefs["beta"])) ^(coefs["beta"] / coefs["AQY"]) # see equation in Platt 1980
    
    # calculate Ik
    Ik_val<- Pmax_val/coefs["AQY"] # Ik = Pmax/AQY
    
    # calculate Iopt b/c if beta =0, there is no finite photoinhibition optimum
    # so Platt model would go toward a more non-photoinhibited curve
    if(coefs["beta"] > 0.000001) {
      Iopt_val <- (coefs["Ps"] / coefs["AQY"]) *
        log((coefs["AQY"] + coefs["beta"]) / coefs["beta"]) # Iopt accounts for beta
      
    } else {
      Iopt_val <- NA
    }
    
    coef_print <- coefs |> # make df like before
      t() |>
      as.data.frame() |> 
      mutate(Pmax = Pmax_val, Ik = Ik_val, Iopt = Iopt_val, species = sp, species_name = sp_name, zone = z, algae_id = i,
             model = "Platt_photoinhibition", AIC = AIC(curve.platt))
    
    coef_list[[i]] <- coef_print # save Ps, alpha, beta, Rd, Iopt, AIC
    
    # make smooth light sequence for fitted curve
    curve_df <- tibble(light_value = seq(0, 900, by = 5)) |>
      mutate(pred_rate = platt_model(light_value, Ps = coefs["Ps"], AQY = coefs["AQY"],  # predicted
                                     beta = coefs["beta"], Rd = coefs["Rd"]),
             species = sp, species_name = sp_name, zone = z, algae_id = i)
    
    curve_list[[i]] <- curve_df # save curve
    
    raw_list[[i]] <- data |> # save raw observations
      select(species, species_name, zone, algae_id, light_value, umol.g.hr)
  }
}
  

### Combine PI Model outputs

# one row per individual with fitted parameters
PI_individuals <- bind_rows(coef_list)

# all individual fitted curves
PI_curve_data <- bind_rows(curve_list)

# raw data points
PI_raw_data <- bind_rows(raw_list)

View(PI_individuals)



### Fit individual PI curves
## Individual curves for each sample (140)
##### INDIVIDUAL PLOTS #####

PI_individuals_plot <- ggplot() +
  geom_point(data = PI_raw_data, # observed points
             aes(x = light_value, y = umol.g.hr, color = species),
             size = 2) +
  geom_line(data = PI_curve_data, # fitted Platt curve
            aes(x = light_value, y = pred_rate, color = species, group = algae_id),
            linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +# line at net photosynthesis = 0
  scale_color_manual(values = species_palette) + # color palette
  facet_wrap(~ algae_id, ncol = 4, scales = "free_y") + # put all plots on one 
  scale_x_continuous(breaks = seq(0, 900, by = 100)) + # x-axis scale
  labs(x = expression("PAR ("*mu*"mol photons "*m^-2*s^-1*")"), # or Photon Flux Density?
       y = expression("Net Photosynthesic Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_classic() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

### Plot PI Model individual curves together
## Plot all inshore individual curves into one PDF and one for offshore

### INSHORE ### 
PI_inshore_raw <- PI_raw_data |>
  filter(zone == "inshore")

PI_inshore_curve <- PI_curve_data |>
  filter(zone == "inshore")






### OFFSHORE ### 










### Bootstrap analysis




### PI curve figures ###

## Within-species: inshore vs offshore

## Within-zone: compare species

## INSHORE ##

## OFFSHORE ##





### Individual parameter tables
# Pmax (Ps)
# AQY
# beta
# Rd
# Ik
# Iopt

### Parameter comparison figures

# Ik
# Pmax
# AQY
# Rd
# beta
# Iopt


### Statistical comparisons
# Species 

# Zone

# Individual vs. bootstrap estimates



### Export final tables and figures



