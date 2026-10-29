###### Respo Code for PI Curves ####### 
### Created by: Nyssa Silbiger, Maya Powell, Haley Poppinga
#### Last updated on: 2026-10-28
# old update with ewt weights: 2026-09-14

############## Introduction to code/script ####################
# Respiration rate calculation
# PI Models
# AIC comparison

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

options(scipen = 999) # makes it so R displays tables not in scientific notation

############## Algae PI Analysis ###############
# Calculate metabolic rates from LoLinR output
# Blank correct and normalize metabolic rates
# Fit PI models and compare AIC

######### Load processed data and metadata ###############

# read Respo_Algae_R.csv file 
RespoR <- read_csv(here("Data","Respo_Files","PI","Respo_Algae_R.csv"))

#### read and join in metadata
# Algae measurements: Load your respiration data file, with all the times, water volumes(mL), wet weight, dry weight, AFDW, zone
BioData <- read_csv(here("Data","Respo_Files","PI","Algae_Measurements_PI.csv")) |> 
  mutate(zone = case_when(str_ends(sample_id, "_IN") ~ "inshore", #create new column for zone
                          str_ends(sample_id, "_OFF") ~ "offshore",
                          TRUE ~ NA_character_))

# PI metadata 
RespoMeta <- read_csv(here("Data","Respo_Files","PI","Algae_PI_meta.csv"))
#View(BioData)
#View(RespoMeta)


# join the data together
Sample_Info <- left_join(RespoMeta, BioData)
#View(Sample_Info)


######### Calculate Respiration rate ###############
############################################

# make sure light_level is the same data type in both dataframes
Sample_Info <- Sample_Info |>
  mutate(light_level = as.integer(light_level))


RespoR2 <- RespoR |> 
  mutate(light_level = as.integer(light_level)) |> # make sure light level is integer
  select(-sample_id, -run_block, -light_value) |> # avoid duplicated columns after join (run_block/light_value already exist in RespoR)
  left_join(Sample_Info, by = c("file_id_csv", "light_level")) |>  #join calculated LoLinR slopes back to sample metadata
# what this does: for each LoLinR slope, find the metadata row with the same oxygen file and the same experimental light level
 # calculate rate using chamber water volumes
  mutate(Ch.Volume.mL = volume_mL, # mL to L conversion
         Ch.Volume.L  = Ch.Volume.mL * 0.001, #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of algae size) and removes per Lite
         umol.sec = umol.L.sec * Ch.Volume.L)
  #mutate_if(sapply(., is.character), as.factor) #convert character columns to factors
      


#Account for blank rate by sample run Block (if we do at least one blank per block)
# 12 samples, 2 of these are blanks

#View(RespoR2)


####### normalize the respo rates to the blanks ##### 
RespoR_Normalized <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, run_block, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE)) %>% # if you have multiple blanks per run take the average
  ungroup() %>% 
  #dplyr::select(Light_level, run_block, blank.rate, date) %>% # this is what I will use to join the blanks back with the raw data
  right_join(RespoR2) %>% # join blanks with the respo data
  mutate(umol.sec.corr = umol.sec - blank.rate, # subtract the blank rates from the raw rates   
         umol.g.hr = (umol.sec.corr*3600)/afdw_g, # change to AFDW when obtained
         umol.g.hr_uncorr = (umol.sec*3600)/afdw_g) %>%
  
  # join blank rates back onto ALL data, then correct, then normalize by afdw weight
  filter(blank!=1) %>% # remove the Blank data
  #ungroup() %>%
  dplyr::select(date,species, full_species, zone, sample_id, algae_id, light_dark,
                run_block, afdw_g, chamber_channel, temp_c, light_level, light_value, # change to AFDW when obtained
                umol.sec,blank.rate, umol.sec.corr, umol.g.hr, umol.g.hr_uncorr) # keep only what we need

#### filter out bad samples here if needed

#######################
### making a df for just blank data for future use in plots ### 
Blank_only <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, light_value, run_block, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

# CHECK THE RAW DATA OF YOUR BLANKS
# if one is high and the other normal, take out the high one
# CHECK FOR OUTLIERS (what does the blank look like, what does the raw data use)


#############################
write_csv(RespoR_Normalized , here("Data","Respo_Files","PI","Respo_Algae_RNormalized_AllPIRates.csv"))  

RespoR_Normalized <- read_csv(here("Data","Respo_Files","PI","Respo_Algae_RNormalized_AllPIRates.csv"))

## Plot the blanks across treatments to make sure nothing is funky
blank_only_plot <- Blank_only %>%
  ungroup() |> 
  ggplot(aes(x = light_value, blank.rate)) +
  geom_point()
blank_only_plot
#ggsave(here("Output","PI","blank_only_plot.pdf"), blank_only_plot, width = 10, height = 10)

#  basic plot of rates versus light before you make the real PI curve 
basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id))+
  geom_point()+
  geom_line()+
  facet_wrap(~species, scales = "free")
basic_PI_plot
#ggsave(here("Output","PI","basic_PI_plot.pdf"), basic_PI_plot, width = 10, height = 10)

#set up species colors
species_palette <- c(av  = "#033500",  # green
                 hd  = "#66A61E",  # green
                 cs  = "#20B2AA",  # green
                 dc  = "#239665",  # green
                 as  = "#E41A1C",  # red
                 gs  = "#FFA500",  # red
                 sf  = "#7C0A02",  # red
                 da  = "#A57A54",  # brown
                 ds  = "#A89F51")  # brown

# another basic plot thats maybe nicer
line_labels <- RespoR_Normalized %>% # one label per line, take the max light point for each algae_id within each species
  group_by(species, algae_id) %>%
  slice_max(light_value, n = 1, with_ties = FALSE) %>% # add label to each id
  ungroup()

algae_basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(data = line_labels, aes(label = algae_id), show.legend = FALSE,
                           size = 3, min.segment.length = 0, segment.color = "black") +
  scale_color_manual(values = species_palette, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 900, by = 100)) +
  coord_cartesian(xlim = c(0, 900)) +
  facet_wrap(~species, scales = "free_y")
#ggsave(here("Output","PI","Algae_Basic_PI_plot.pdf"), algae_basic_PI_plot, width = 10, height = 10)



### run an nls model for PI curve and extract Ik for each species ###

##### PLOTTING CURVES #####
############################################
# Plot PI curves using 3 different PI models and them compare AIC

### Read in rest of data
# full species names df
sp_names <- read_csv(here("Data/species_names.csv"))
# join species names with the normalized data
RespoR_Normalized <- RespoR_Normalized %>% left_join(sp_names, by = "species")

# full_species list:
# av, as, gs, da, ds, dc, sf, cs, hd

###### DEFINE MODEL #######

### MODEL #1: Marshall & Biscoe (1980) Non-Rectangular Hyperbola ###
# Nonlinear Least Squares regression of a non-rectangular hyperbola 
# maybe derived from (Marshall & Biscoe, 1980), I think from Maya's MCAP/PCOM code
### Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#using nls multstart
#define nrh model & function to make sure all variables are not saved from previous runs
nrh_model <- function(light_value, Pmax, AQY, Rd, theta) {
  (1/(2*theta)) *(AQY*light_value + Pmax -sqrt((AQY*light_value + Pmax)^2 -4*AQY*theta*Pmax*light_value)) - Rd
}


##### MODEL PARAMETERS ######
#Pmax = "Am" = (max photosynthetic rate) 
    # should be positive based on range of umol.cm2.hr - here ~0-2
    # fitted maximum gross photosynthetic capacity, can be higher than highest observed net rate 
         # b/c model estimates gross capacity and then subtracts respiration

#AQY (apparent quantum yield) = alpha
    # should be positive and between 0-1 (often much smaller)
    # high alpha -> alga responds strongly to small increases in light)
    # initial slope of the PI curve at low irradiance
    # represents photochemical efficiency at low light

#Rd = dark respiration
    # should be subtracted and usually between -2-0 but can adjust if larger
    # magnitude of the negative dark metabolic rate
    # fitted magnitude of dark respiration, the equation ends in - Rd, 
    # predicted net rate in darkness

#theta = fitted curvature 
    # is between 0 and 1 
    # theta controls bend btwn initial light-limited increase & the saturated part of the curve)
    # closer to 1 --> sharper transition toward saturation (sharp from the low-light rising portion toward saturation)
    # lower theta --> more rounded transition


### make function used to plot Model 1 curve
make_nrh_fun <- function(Pmax, AQY, Rd, theta) {
  force(Pmax); force(AQY); force(Rd); force(theta)
  function(x) nrh_model(x, Pmax, AQY, Rd, theta)
}

# create plot list and a coefficient table list for Model 1
# lists to save outputs from each individual model fit
coef_list_Model1 <- list()  # parameter estimates for each individual
curve_list_Model1 <- list() # smooth fitted curve for each 
raw_list_Model1 <- list()   # observed PI points for each

# MODEL 1 forloop: to calculate ik and plot
# fit nrh model to each species x zone population

### make a forloop
#Set i = "AS" or whatever you want to look at for each species_name as you go through
for (i in unique(RespoR_Normalized$algae_id)) {
  data <- RespoR_Normalized |> 
    filter(algae_id == i) |>  # subset individual's PI curve
    arrange(light_value) # not light level b/c light value is the actual measured irradiance used by the model
    
  
  # identify the species and zone for individual sample
  sp <-data$species[1]
  sp_name <- data$species_name[1]
  z <- data$zone[1]
  cat("\nModel 1:", sp, z, i, "\n")
  
  
  #use nls_multstart to calculate the curve
  curve.nlslrc <- tryCatch(nls_multstart( # tryCatch for the error function later
    umol.g.hr ~ nrh_model(light_value, Pmax, AQY, Rd, theta), #@haley you can just swap out this model with different ones here that you define above
    data = data,
    iter = 500,
    #IF YOU ARE NORMALIZING TO G INSTEAD OF CM2 OR ANOTHER METRIC YOU SHOULD LOOK AT YOUR DATA FIRST TO SEE BOUNDS FOR Pmax, AQY AND RD
    start_lower = c(Pmax = max(data$umol.g.hr, na.rm = TRUE)*0.5,     AQY = 0.001,  Rd = abs(min(data$umol.g.hr, na.rm = TRUE))*0.5,  theta = 0.1), # half of Pmax and min Rd as lower start
    start_upper = c(Pmax = max(data$umol.g.hr, na.rm = TRUE)*3, AQY = 0.5,  Rd = abs(min(data$umol.g.hr, na.rm = TRUE))*2,  theta = 1), # triple of Pmax and doubld of min Rd as upper start
    lower = c(Pmax = 0, AQY = 0, Rd = 0, theta = 0.001),
    supp_errors = 'Y'
  ),
  
  # need to add this error part for the species that go into photoinhibition
  # if one individual fails, continue to next one
  error = function(e) {
    message("MODEL 1 FAILED FOR ", i, ": ", e$message)
    return(NULL)
  }
  )

# skip to next individual if model failed
if(is.null(curve.nlslrc)) {
  next
}

# extract fitted coefficients
  coefs <- coef(curve.nlslrc)
  coef_print <- coefs %>% # coefficient estimates for the individual
    t() %>%
    as.data.frame() %>%
    rename(Pmax.gross = Pmax) %>%
    mutate(Ik = Pmax.gross/AQY,
           Ic = Rd/AQY,
           Pmax.net = Pmax.gross - Rd,
           species = sp,
           species_name = sp_name,
           zone = z,
           algae_id = i,
           model = "NRH",
           AIC = AIC(curve.nlslrc)) # for comparison between models
  
  #get Ik value (and force so that it updates each loop)
  ik_val <- coef_print$Ik
  force(ik_val)
  
  # save coefficient table for this individual
  coef_list_Model1[[i]] <- coef_print
  
  # fitted function for plotting
  fitted_fun <- make_nrh_fun(coefs["Pmax"], coefs["AQY"], coefs["Rd"], coefs["theta"])
  
  # make a smooth light sequence for the fitted curve (so that all curves can be plotted on one plot)
  curve_df <- tibble(light_value = seq(0, 900, by = 5)) |>
    mutate(pred_rate = fitted_fun(light_value), # photosynthetic rate predicted by a sample's fitted NRH model at that light value
           species = sp, species_name = sp_name, zone = z, algae_id = i)
  
  # save predicted fitted curve values
  curve_list_Model1[[i]] <- curve_df
  
  # save raw observed data for this individual
  raw_list_Model1[[i]] <- data |>
    select(species, species_name, zone, algae_id, light_value, umol.g.hr)
  
}

# Combine all individual Model 1 results into one df
Model1_individuals <- bind_rows(coef_list_Model1)

# fitted curve values in one df
Model1_curve_data <- bind_rows(curve_list_Model1)

# raw observed points in one df
Model1_raw_data <- bind_rows(raw_list_Model1)

# table of fitted parameters for each individual
Model1_Ik_table <- Model1_individuals |>
  select(species, species_name, zone, algae_id, 
         Ik, Pmax.gross, AQY, Rd, theta, AIC) |>
  arrange(species, zone, algae_id)

View(Model1_Ik_table)

### Plot Model 1 individual curves together
# get the end of each fitted curve so it can be labeled
Model1_curve_labels <- Model1_curve_data |>
  group_by(species, species_name, zone, algae_id) |>
  slice_max(light_value,n = 1, with_ties = FALSE) |>
  ungroup()

Model1_individualcurves_plot <- ggplot() +
  geom_point(data = Model1_raw_data, # actual observed PI points from each individual
             aes(x = light_value, y = umol.g.hr, color = species),
             size = 1.5, alpha = 0.6) +
  geom_line(data = Model1_curve_data, # fitted PI curve for each individual
            aes(x = light_value, y = pred_rate, color = species, group = algae_id),
            linewidth = 0.8) +
  ggrepel::geom_text_repel(data = Model1_curve_labels, # add individual labels to curves
                           aes(x = light_value, y = pred_rate, label = algae_id, color = "black"),
                           size = 2.5, nudge_x = 25, direction = "y", hjust = 0, show.legend = FALSE) +
  coord_cartesian(xlim = c(0, 1000), clip = "off") + # space after 900
  facet_wrap(vars(species_name, zone), # one panel for each species x zone
             scales = "free_y") +
  scale_color_manual(values = species_palette) + # species color pallette
  scale_x_continuous(breaks = seq(0, 900, by = 100)) + # light value breaks
  labs(title = "Model 1: Individual PI Curves",
       x = expression("PAR ("*mu*"mol photons "*m^-2*s^-1*")"),
       y = expression("Photosynthetic Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_classic() +
  theme(legend.position = "none", strip.text = element_text(face = "italic"))

Model1_individualcurves_plot
#ggsave(here("Output","PI","Model1_individualcurves_plot.pdf"), Model1_individualcurves_plot, width = 20, height = 10)



############################################################
############################################################
############################################################

###### DEFINE MODEL #######

### MODEL 2: Platt et al. (1980) Photoinhibition Model ###
platt_model <- function(light_value, Ps, AQY, beta, Rd) {
  Ps *(1 - exp(-(AQY*light_value / Ps))) *exp(-(beta * light_value / Ps)) - Rd
} #exp is the photoinhibition


# Another possible model but wont try
### MODEL 3: Jassby & Platt (1976) Hyperbolic Tangent ###
#jp_model <- function(light_value, Pmax, AQY, Rd) {
#  Pmax * tanh(AQY * light_value / Pmax) - Rd
#}

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

# lists to save outputs from each individual Model 2 fit
coef_list_Model2 <- list()  # fitted parameters for each individual
curve_list_Model2 <- list() # smooth fitted curve for each 
raw_list_Model2 <- list()   # observed raw data for each


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
  
  coef_list_Model2[[i]] <- coef_print # save Ps, alpha, beta, Rd, Iopt, AIC
 
   # make smooth light sequence for fitted curve
  curve_df <- tibble(light_value = seq(0, 900, by = 5)) |>
    mutate(pred_rate = platt_model(light_value, Ps = coefs["Ps"], AQY = coefs["AQY"],  # predicted
                                   beta = coefs["beta"], Rd = coefs["Rd"]),
           species = sp, species_name = sp_name, zone = z, algae_id = i)
  
  curve_list_Model2[[i]] <- curve_df # save curve

  raw_list_Model2[[i]] <- data |> # save raw observations
    select(species, species_name, zone, algae_id, light_value, umol.g.hr)
}


### Combine Model 2 outputs

# one row per individual with fitted parameters
Model2_individuals <- bind_rows(coef_list_Model2)

# all individual fitted curves
Model2_curve_data <- bind_rows(curve_list_Model2)

# raw data points
Model2_raw_data <- bind_rows(raw_list_Model2)

View(Model2_individuals)
  

### Plot Model 2 individual curves together
Model2_curve_labels <- Model2_curve_data |> # end of each fitted curve so it can be labeled
  group_by(species, species_name, zone, algae_id) |>
  slice_max(light_value,n = 1, with_ties = FALSE) |>
  ungroup()

Model2_individualcurves_plot <- ggplot() +
  geom_point(data = Model2_raw_data, # actual observed PI points from each individual
             aes(x = light_value, y = umol.g.hr, color = species),
             size = 1.5, alpha = 0.6) +
  geom_line(data = Model2_curve_data, # fitted PI curve for each individual
            aes(x = light_value, y = pred_rate, color = species, group = algae_id),
            linewidth = 0.8, alpha = 0.8) +
  ggrepel::geom_text_repel(data = Model2_curve_labels, # add individual labels to curves
                           aes(x = light_value, y = pred_rate, label = algae_id, color = "black"),
                           size = 2.5, nudge_x = 25, direction = "y", hjust = 0, show.legend = FALSE) +
  coord_cartesian(xlim = c(0, 1000), clip = "off") + # space after 900
  facet_wrap(vars(species_name, zone), # one panel for each species and zone
             scales = "free_y") +
  scale_color_manual(values = species_palette) + # species color palette
  scale_x_continuous(breaks = seq(0, 900, by = 100)) + # light value breaks
  labs(title = "Model 2: Individual Photoinhibition PI Curves",
       x = expression("PAR ("*mu*"mol photons "*m^-2*s^-1*")"),
       y = expression("Net Photosynthetic Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_classic() +
  theme(legend.position = "none", strip.text = element_text(face = "italic"))

Model2_individualcurves_plot
#ggsave(here("Output","PI","Model2_individualcurves_plot.pdf"), Model2_individualcurves_plot, width = 20, height = 10)

# make a table of model parameters 
Model2_parameter_table <- Model2_individuals |>
  select(species, species_name, zone, algae_id, Ps, Pmax, AQY, beta, Rd, Ik, Iopt, AIC) |>
  arrange(species, zone, algae_id)

View(Model2_parameter_table)

# NEXT STEPS: make boxplots of Ik, Pmax, AQY, and Rd
# run ANOVA or similar tests




########################################################################
### Now do this for population data (bootstrapped? error?)
# Model 2 only
# compare confidence intervals around parameters

## Prep PI population data
# calc means and SE for each species x zone at each experimental light level
PI_summary <- RespoR_Normalized |> 
  filter(!is.na(zone)) |> 
  group_by(species, species_name, zone, light_level) |>  
  summarise(light_value = mean(light_value, na.rm = TRUE),
            mean_rate = mean(umol.g.hr, na.rm = TRUE),
            n = sum(!is.na(umol.g.hr)),
            se_rate = sd(umol.g.hr, na.rm = TRUE) / sqrt(n),
            n_individuals = n_distinct(algae_id),
            .groups = "drop")


# df for all populations: list of every species x zone population
PI_populations <- RespoR_Normalized |>
  filter(!is.na(zone)) |>
  distinct(species, species_name, zone)
PI_populations

# save actual population Model 2 fits
coef_list_Model2_pop <- list()
curve_list_Model2_pop <- list()

# save bootstrap results
coef_list_Model2_boot <- list()
curve_list_Model2_boot <- list()

# make forloop... again
for (i in 1:nrow(PI_populations)) {
  
  # identify species x zone population
  sp <- PI_populations$species[i]
  sp_name <- PI_populations$species_name[i]
  z <- PI_populations$zone[i]
  cat("\nPopulation Model 2:", sp, z, "\n")

  # actual mean PI curve for this species and zone
  pop_mean <- PI_summary |>
    filter(species == sp, zone == z) |>
    arrange(light_level)
  
  # all individual PI curves from population
  pop_data <- RespoR_Normalized |>
    filter(species == sp, zone == z, !is.na(algae_id))
  
  # fit Model 2 to actual population mean curve
  fit_pop <- tryCatch(nls_multstart(
      mean_rate ~ platt_model(light_value, Ps, AQY, beta, Rd),
      data = pop_mean,
      iter = 500,
      start_lower = c(Ps = max(pop_mean$mean_rate, na.rm = TRUE)*0.5,
                      AQY = 0.001, beta = 0.0001, Rd = abs(min(pop_mean$mean_rate, na.rm = TRUE))*0.5),
      start_upper = c(Ps = max(pop_mean$mean_rate, na.rm = TRUE)*3, AQY = 0.5,
                      beta = 0.5, Rd = abs(min(pop_mean$mean_rate, na.rm = TRUE))*2),
      lower = c(Ps = 0.000001, AQY = 0.000001, beta = 0, Rd = 0),
      supp_errors = "Y"),
    
    error = function(e) {
      message("POPULATION MODEL 2 FAILED FOR ", sp, " ", z)
      return(NULL)
    }
  )
  if(is.null(fit_pop)) {
    next
  }
  
  pop_coefs <- coef(fit_pop)
  
  # calculate optimal irradiance (Iopt) only if beta > 0
  if(pop_coefs["beta"] > 0.000001) {
    Iopt_pop <- (pop_coefs["Ps"] / pop_coefs["AQY"]) *
      log((pop_coefs["AQY"] + pop_coefs["beta"]) /pop_coefs["beta"]
      )
  } else {
    Iopt_pop <- NA
  }
  
  pop_name <- paste(sp, z, sep = "_")
  
  coef_list_Model2_pop[[pop_name]] <- tibble(
    species = sp, species_name = sp_name, zone = z, Ps = pop_coefs["Ps"],
    AQY = pop_coefs["AQY"], beta = pop_coefs["beta"], Rd = pop_coefs["Rd"], 
    Ik = pop_coefs["Ps"]/pop_coefs["AQY"], Iopt = Iopt_pop, AIC = AIC(fit_pop))
  
  # make curve
  curve_list_Model2_pop[[pop_name]] <- tibble(
    light_value = seq(0, 900, by = 5)) |>
    mutate(pred_rate = platt_model(light_value, Ps = pop_coefs["Ps"], AQY = pop_coefs["AQY"],
                                   beta = pop_coefs["beta"], Rd = pop_coefs["Rd"]),
           species = sp,species_name = sp_name, zone = z)
  
  ### add bootstrap error loop
  set.seed(123) # reproduce same boostrap results
  n_boot <- 100 # for each population repeat bootstrap 1000 times
  
  # list of unique biological replicates in this population
  individual_ids <- unique(pop_data$algae_id)
  
  for (b in 1:n_boot) {
    sampled_ids <- sample(individual_ids, # bootstrap resampling
                          size = length(individual_ids),replace = TRUE)
    # df duplicated sampled individuals
    boot_ids <- tibble(boot_id = seq_along(sampled_ids), algae_id = sampled_ids)
    boot_data <- boot_ids |> 
      left_join(pop_data, by = "algae_id")
  
    # recalculate population mean for bootstrap
    boot_individual <- boot_data |>
      group_by(boot_id, light_level) |>
      summarise(light_value = mean(light_value, na.rm = TRUE),
                rate = mean(umol.g.hr, na.rm = TRUE), .groups = "drop")
    boot_mean <- boot_individual |> # average resampled individuals
      group_by(light_level) |>
      summarise(light_value = mean(light_value, na.rm = TRUE), mean_rate = mean(rate, na.rm = TRUE), .groups = "drop")
  
    # fit Platt (Model 2) to boostrapped mean curve
    fit_boot <- tryCatch(nls_multstart(
      mean_rate ~ platt_model(light_value, Ps, AQY, beta, Rd),
        data = boot_mean,
        # bootstrap means should be smoother, fewer multistarts usually enough
        iter = 100,
        start_lower = c(Ps = max(boot_mean$mean_rate, na.rm = TRUE)*0.5,
                        AQY = 0.001,beta = 0.0001, Rd = abs(min(boot_mean$mean_rate, na.rm = TRUE))*0.5),
        start_upper = c(Ps = max(boot_mean$mean_rate, na.rm = TRUE)*3,
                        AQY = 0.5, beta = 0.5, Rd = abs(min(boot_mean$mean_rate, na.rm = TRUE)) * 2),
        lower = c(Ps = 0.000001, AQY = 0.000001,beta = 0, Rd = 0),
        supp_errors = "Y"),
      
      error = function(e) { # so if it fails, loop still trys for next ones
        return(NULL)
      }
    )
    
    if(is.null(fit_boot)) {
      next
    }
    
    boot_coefs <- coef(fit_boot) # save boostrap fitted curve
    boot_name <- paste(sp, z, b, sep = "_")
    
    curve_list_Model2_boot[[boot_name]] <- tibble( # save curve
      light_value = seq(0, 900, by = 5)) |>
      mutate(pred_rate = platt_model(light_value, Ps = boot_coefs["Ps"],
                                     AQY = boot_coefs["AQY"], beta = boot_coefs["beta"],
                                     Rd = boot_coefs["Rd"]),
             species = sp, species_name = sp_name, zone = z, bootstrap = b)
  
    coef_list_Model2_boot[[boot_name]] <- tibble(
      species = sp, species_name = sp_name, zone = z,bootstrap = b,  Ps = boot_coefs["Ps"],
      AQY= boot_coefs["AQY"], beta = boot_coefs["beta"],  Rd = boot_coefs["Rd"])
    
  } # end bootstrap loop
  
} # end population loop

### Combine everything
# actual population fits
Model2_population_results <- bind_rows(coef_list_Model2_pop)
Model2_population_curve <- bind_rows(curve_list_Model2_pop)

# bootstrap fits
Model2_boot_parameters <- bind_rows(coef_list_Model2_boot)
Model2_boot_curves <- bind_rows(curve_list_Model2_boot)

# calculate 95% bootstrap confidence 
Model2_bootstrap <- Model2_boot_curves |>
  group_by(species, species_name, zone, light_value) |>
  summarise(lower = quantile(pred_rate, 0.025, na.rm = TRUE),
            median = quantile(pred_rate, 0.50, na.rm = TRUE),
            upper = quantile(pred_rate, 0.975, na.rm = TRUE),
            n_boot = n_distinct(bootstrap), .groups = "drop")
  

# make plot
Model2_Ik_labels <- Model2_population_results |> # add Ik labels
  mutate(Ik_label = paste0("Ik = ", round(Ik)))

Model2_population_plot <- ggplot() +
  geom_ribbon(data = Model2_bootstrap,  # 95% bootstrap confidence interval
              aes(x = light_value, ymin = lower, ymax = upper, fill = species), alpha = 0.2) +
  geom_line(data = Model2_population_curve, # actual fitted population Model 2 curve
            aes(x = light_value, y = pred_rate, color = species), linewidth = 1) +
  geom_point(data = PI_summary,  # actual observed population means
             aes(x = light_value, y = mean_rate, color = species), size = 2) +
  geom_errorbar(data = PI_summary, # standard error of observed population means
                aes(x = light_value, ymin = mean_rate - se_rate, ymax = mean_rate + se_rate, color = species), width = 10) +
  geom_vline(data = Model2_population_results, # adds Ik dashed line
             aes(xintercept = Ik, color = species), linetype = 2, linewidth = 0.7, show.legend = FALSE) +
  geom_text(data = Model2_Ik_labels, # adds Ik label
            aes(x = Ik, y = -Inf, label = Ik_label, color = species),
            hjust = -0.1, vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(values = species_palette) + # species palette
  scale_fill_manual(values = species_palette) +
  facet_wrap(vars(species_name, zone), scales = "free_y",
             labeller = labeller(species_name = function(x) paste0("<i>", x, "</i>"),
                                 zone = function(x) x, .multi_line = TRUE)) +
  scale_x_continuous(breaks = seq(0, 900, by = 100)) +
  labs(title = "Population PI Curves with Bootstrap Confidence Intervals",
       subtitle = "Platt et al. (1980) Photoinhibition Model",
       x = expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
       y = expression("Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_classic() +
  theme(legend.position = "none", strip.text = ggtext::element_markdown(face = "plain"))

Model2_population_plot
#ggsave(here("Output","PI","Model2_population_plot.pdf"), Model2_population_plot, width = 20, height = 12)





