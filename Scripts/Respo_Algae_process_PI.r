###### Respo Code for Algae PI Curve ####### 
### Created by: Haley Poppinga, Maya Powell, Nyssa Silbiger
#### Last updated on: 2026-09-02
#### PI Update: 2026-06-01, 2026-08-10, 2026-08-18

############## Introduction to code/script ####################
## this script will help us process the raw data gathered during respirometry runs. 
## need to change for specific project/experimental variables 

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
if ("nls.multstart" %in% rownames(installed.packages()) == 'FALSE') install.packages('nls.multstart')
if ("rTPC" %in% rownames(installed.packages()) == 'FALSE') install.packages('rTPC')

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

############# now it's time to code ############
################################################
# get the file path

#set the path to all of the raw oxygen datasheets
## these are saved onto the computer in whatever file path/naming scheme you saved things to 
path.p<-here("Data","Respo_Files","PI","RawO2") #the location of all your respirometry files
#you can change to individual run folders if needed

# bring in all of the individual files
filenames_final<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders

#basename above removes the subdirectory name from the file, re-name as file.names.full
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#empty chamber volume
ch.vol <- 482 #mL #of small chambers 12

######### Load and tidy files ###############
############################################
#Load your respiration data file, with all the times, water volumes(mL), #not doing dry weight just SA
#RespoMeta <- read_csv(here("Data","RespoFiles","Respo_Metadata_SGDDilutions_Cabral_Varari.csv"))
BioData <- read_csv(here("Data","Respo_Files","PI","Algae_Measurements_PI.csv")) |> 
  mutate(zone = case_when(str_ends(sample_id, "_IN") ~ "inshore", #create new column for zone
                          str_ends(sample_id, "_OFF") ~ "offshore",
                          TRUE ~ NA_character_))

RespoMeta <- read_csv(here("Data","Respo_Files","PI","Algae_PI_meta.csv"))
#View(BioData)
#View(RespoMeta)
## try first with prelim fake data to make sure script runs
## then switch to real calculated data after getting volumes and weight and surface area


# join the data together
Sample_Info <- left_join(RespoMeta, BioData) %>% #, by = c("sample_id", "algae_id"))
  mutate(light_level = as.integer(light_level)) %>% 
  #select(file_id_csv, light_level, date, start_time, stop_time, species, algae_id, sample_id, light_dark, 
  #blank, chamber_channel, run_block, light_value, volume_mL, wetweight_g) # only what we need
  #select(file_id_csv, light_level, date, start_time, stop_time, species, algae_id, light_value, light_dark, blank, run_block, chamber_channel, volume_mL, wetweight_g)
  select(file_id_csv, light_level, date, start_time, stop_time,
         species, full_species, algae_id, zone, light_value, light_dark, 
         blank, run_block, chamber_channel, volume_mL, wetweight_g)
# intentionally NOT including sample_id, temp_c already in RespoR


#View(Sample_Info)



##### Make sure times are consistent ####
# make start and stop times real times, so that we can join the respo output and sample_info data frames
Sample_Info <- Sample_Info %>% 
  #drop_na(sample_ID) %>% 
  unite(date,start_time,col="start_time",remove=F, sep=" ") %>% 
  unite(date,stop_time,col="stop_time",remove=F, sep=" ") %>%
  mutate(start_time = ymd_hms(start_time)) %>% 
  mutate(stop_time = ymd_hms(stop_time)) %>% 
  mutate(date = ymd(date)) %>% 
  mutate(light_level = as.integer(light_level))

#view(Sample_Info)

#generate a 4 column dataframe with specific column names
# data is in umol.L.sec

n_light_levels<-10 # number of unique light levels

RespoR <- tibble(.rows = length(filenames_final)*n_light_levels,
                 file_id_csv = NA_character_,
                 sample_id = NA_character_,
                 intercept = NA_real_,
                 umol.L.sec = NA_real_,
                 temp_c = NA_real_,
                 light_level = NA_integer_,
                 light_value = NA_real_,
                 run_block = NA_character_)



# create directory for output folder
out_dir <- here("Output", "PI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



######### Create a for loop! ###############
############################################

###forloop#####
for(i in 1:length(filenames_final)) { # loop each raw file?
  FRow <- as.numeric(which(Sample_Info$file_id_csv==filenames_final[i])) # stringsplit this renames our file
  
  Respo.Data1 <- read_csv(skip=1,file.path(path.p, paste0(file.names.full[i]))) %>% # reads in each file in list
    dplyr::select(Date, Time, Value, Temp) %>% # keep only what we need: Time stamp per 1sec, Raw O2 value per 1sec, in situ temp per 1sec
    unite(Date,Time,col="Time",remove=T, sep = " ") %>%
    drop_na() %>% 
    mutate(Time = mdy_hms(Time), Value = as.numeric(Value), Temp  = as.numeric(Temp))
  #filter(Time >= start_time & Time <= stop_time)
  #mutate(Time = ymd_hms(Time)) #%>% # convert time
  #mutate(help = i) ##if stuck in forloop with error from filter, can check RespoR and see at what row the forloop stopped working  
  
  ## cut the data by start and stop times from metadata
  #Use start time of each light step from the metadata to separate data by light stop
  
  oxy_subsets <- Sample_Info[FRow,] %>%
    pmap(function(light_level, start_time, stop_time, ...) {
      data <- Respo.Data1  %>%
        filter(Time >= start_time & Time <= stop_time) %>%
        arrange(Time) %>%
        mutate(t_sec = as.numeric(difftime(Time, first(Time), units = "secs"))) %>% #keep everything in seconds
        mutate(light_level = light_level) %>%
        filter(t_sec > 120) %>%                          # drop first 2 min (120 s)
        filter(row_number() %% 10 == 0)                  # keep every 10th row
    }) 
  
  
  # Combine into one long dataframe with ID labels
  combined_oxy <- bind_rows(oxy_subsets)
  
  # Get the filename without the .csv
  rename<- sub("_O2.csv","", filenames_final[i])
  
  ### plot and export the thinned data ####
  p1<- ggplot(combined_oxy, aes(x = t_sec, y = Value)) +
    geom_point(color = "dodgerblue") +
    labs(
      x = 'Time (seconds)',
      y = expression(paste(' O'[2],' (',mu,'mol/L)')),
      title = "original"
    )+
    facet_wrap(~light_level, scales = "free_y")
  
  
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "e "pc")
  
  # Define function for fitting LoLinR regressions to be applied to all intervals for all samples
  fit_reg <- function(data) {
    rankLocReg(xall = data$t_sec, yall = data$Value, 
               alpha = 0.2, method = "pc", verbose = FALSE)
  }
  
  # Setup for parallel processing
  future::plan(multisession)
  
  # Map LoLinR function onto all intervals of each sample's thinned dataset
  df <- combined_oxy %>%
    select(t_sec, Value, light_level, Temp)%>%
    mutate(t_sec = as.numeric(t_sec))%>%
    nest_by(light_level) %>%
    ungroup()%>%
    mutate(regs = furrr::future_map(data, fit_reg), # run the LOLinR fit in parallel
           temp_c = map_dbl(map(data, "Temp"), mean),# get the mean temperature
           RegStats =map(regs, function(x){ # extract the intercept and slope for the parameters
             x$allRegs %>%
               slice(1) %>%
               select(intercept = b0,
                      umol.L.sec = b1)
           }) )
  
  
  #  Plot regression diagnostics
  
  for(j in 1:length(df$light_level)){
    pdf(file.path(out_dir, paste0(rename, "_", j, ".pdf")))
    plot(df$regs[[j]])
    dev.off() 
  }
  
  
  # attach meta for this file (by light_level)
  meta_steps <- Sample_Info[FRow,] %>% select(light_level, light_value, run_block)
  
  df <- df %>%
    select(light_level, temp_c, RegStats) %>%
    unnest(RegStats) %>%
    left_join(meta_steps, by="light_level") %>%
    mutate(sample_id = rename)   # keep your original sample_id style
  
  
  
  
  ################################
  # index for the 10 light steps of this file
  idx <- ((i - 1) * n_light_levels + 1):(i * n_light_levels)
  
  # store file id for join later
  RespoR[idx, "file_id_csv"] <- filenames_final[i]
  
  
  # fill in all the O2 consumption and rate data
  RespoR[idx,"temp_c"]      <- df$temp_c
  RespoR[idx,"sample_id"]   <- df$sample_id
  RespoR[idx,"intercept"]   <- df$intercept
  RespoR[idx,"umol.L.sec"]  <- df$umol.L.sec
  RespoR[idx,"light_level"] <- df$light_level
  RespoR[idx,"light_value"] <- df$light_value
  RespoR[idx,"run_block"]   <- df$run_block
  
}  

######### end of for loop - celebrate victory of getting through that ###############
############################################




#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
write_csv(RespoR, here("Data","Respo_Files","PI","Respo_Algae_R.csv"))  

##### 

RespoR <- read_csv(here("Data","Respo_Files","PI","Respo_Algae_R.csv"), show_col_types = FALSE) #%>%
#select(file_id_csv, light_level, intercept, umol.L.sec, temp_c) %>%
#mutate(light_level = as.integer(light_level))




######### Calculate Respiration rate ###############

#avoid duplicated columns after join (run_block/light_value already exist in RespoR)
RespoR2 <- RespoR %>%
  mutate(light_level = as.integer(light_level)) %>% 
  left_join(Sample_Info) %>% #mutate(light_level = as.integer(light_level)),
  #by = c("file_id_csv", "light_level")) %>%
  mutate(Ch.Volume.mL = volume_mL,
         Ch.Volume.L  = Ch.Volume.mL * 0.001,
         umol.sec = umol.L.sec * Ch.Volume.L)

#Account for blank rate by sample run Block (if we do at least one blank per block)
#View(RespoR)



####### Normalize the respo rates to the blanks ########

blank_rates <-RespoR2 %>% # compute blank rates first then left join instead here
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, run_block, light_dark) %>%
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

RespoR_Normalized <- RespoR2 %>% 
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  #summarise(blank.rate = mean(umol.sec, na.rm = TRUE)) %>% # if you have multiple blanks per run take the average
  #ungroup() %>% 
  #dplyr::select(Light_level, run_block, blank.rate, date) %>% # this is what I will use to join the blanks back with the raw data
  #right_join(RespoR2) %>% # join blanks with the respo data
  
  # filter out all bad samples
  #filter(!(sample_id == "AV02_PI_RUN1"),
 #        !(sample_id == "AV03_PI_RUN1" & light_level == 3),
 #        !(sample_id == "AV05_PI_RUN3"),
 #        !(sample_id == "DA01_PI_RUN1" & light_level == 5),
  #       !(sample_id == "DA02_PI_RUN3"),
 #        !(sample_id == "GS02_PI_RUN1"),
  #       !(sample_id == "HD02_PI_RUN3" & light_level == 5),
  #       !(sample_id == "SF05_PI_RUN3" & light_level == 3 & light_level == 4)) %>% 
  
  # join blank rates back onto ALL data, then correct, then normalize by wet weight
  filter(blank != 1) %>%  # remove blank rows from the "sample" dataset, # remove the Blank data
  left_join(blank_rates, by = c("light_level", "run_block", "light_dark")) %>% 
  mutate(umol.sec.corr   = umol.sec - blank.rate, # subtract the blank rates from the raw rates
         umol.g.hr = (umol.sec.corr * 3600) / wetweight_g,
         umol.g.hr_uncorr= (umol.sec * 3600) / wetweight_g) %>%
  dplyr::select(date,species, full_species, zone, sample_id, algae_id, light_dark,
                run_block, wetweight_g, chamber_channel, temp_c, light_level, light_value,
                umol.sec,blank.rate, umol.sec.corr, umol.g.hr, umol.g.hr_uncorr)
  #dplyr::select(date, species, sample_id, algae_id, light_dark, run_block, wetweight_g, chamber_channel, #FIX .x AND .y COLUMNS
                #temp_c, light_level, light_value, umol.sec, blank.rate, umol.sec.corr, 
                #umol.g.hr, umol.g.hr_uncorr) #keep only what we need

###### Filter out bad samples #############
#AS0 (stressed)
#AV0
#AV0
#AV0
#DA0
#DA0
#DA0 (stressed)
#GS0
#GS0 (stressed)
#HD0
#SF0 (stressed)
#SF0 (stressed)
#SF05_3


#### making a df for just blank data for future use in plots #### 
Blank_only <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, light_value, run_block, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

write_csv(RespoR_Normalized , here("Data","Respo_Files","PI","Respo_Algae_RNormalized_AllPIRates.csv"))  


## Plot the blanks across treatments to make sure nothing is funky
Blank_only %>%
  ggplot(aes(x = light_value, blank.rate, group = interaction(run_block, light_dark))) +
  geom_point() +
  geom_line() +
  facet_wrap(~run_block)

#  basic plot of rates versus light before you make the real PI curve 
#basic_PI_plot <- RespoR_Normalized %>%
#  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id)) +
#  geom_point()+
#  geom_line()+
#  facet_wrap(~species, scales = "free_y")


# one label per line, take the max light point for each algae_id within each species
line_labels <- RespoR_Normalized %>%
  group_by(species, algae_id) %>%
  slice_max(light_value, n = 1, with_ties = FALSE) %>% # add label to each id
  ungroup()

# create color palette to match algae classification
# other hex colors: #2E8B57, #9FE2BF, #FB6A4A, #380000, 
#5C4827, #01FEC6, #610C04, #696006, #0B3F33, #8C510A, #B5A71F

species_pal <- c(av  = "#033500",  # green
                 hd  = "#66A61E",  # green
                 cs  = "#20B2AA",  # green
                 dc  = "#239665",  # green
                 as  = "#E41A1C",  # red
                 gs  = "#FFA500",  # red
                 sf  = "#7C0A02",  # red
                 da  = "#A57A54",  # brown
                 ds  = "#A89F51")  # brown

algae_basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(data = line_labels, aes(label = algae_id), show.legend = FALSE,
                           size = 3, min.segment.length = 0, segment.color = "black") +
  scale_color_manual(values = species_pal, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 900, by = 100)) +
  coord_cartesian(xlim = c(0, 900)) +
  facet_wrap(~species, scales = "free_y")

ggsave(here("Output","PI","Algae_Basic_PI_plot.pdf"), algae_basic_PI_plot, width = 10,
       height = 10)

### run an nls model for PI curve and extract Ik for each species ###


############################################
#### TESTING PI CURVE MODELS ####


#Plot PI curves using 3 different PI models and them compare AIC

#algae.data <- read.table("Data/Respo_Files/PI/Respo_Algae_RNormalized_AllPIRates.csv", header=TRUE, sep=",")

# means and se for each population inshore and offshore for each light level
PI_summary <- RespoR_Normalized |> 
  filter(!is.na(zone)) |> #filter inshore/offshore
  group_by(species, full_species, zone, light_level) |>  
  summarise(light_value = mean(light_value, na.rm = TRUE), #calc mean light
            mean_rate = mean(umol.g.hr, na.rm = TRUE), #calc mean rates
            n = sum(!is.na(umol.g.hr)),
            se_rate = sd(umol.g.hr, na.rm = TRUE) / sqrt(n), #error
            n_individuals = n_distinct(algae_id),
            .groups = "drop")

# Test one species x zone first
# Acantophora spicifera inshore population
as_in <- PI_summary |> 
  filter(species == "as", zone == "inshore") |> 
  arrange(light_level) # put experimental light treatments in order

# View(as_in) # look at data for this PI curve
as_in # see tibble in console

# Plot this test population data only
as_in |> 
  ggplot(aes(x = light_value, y = mean_rate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, # plot error
                    ymax = mean_rate + se_rate), width = 10) +
  geom_line() + 
  scale_x_continuous(breaks = seq(0, 900, by = 100)) +
  labs(title = "Acanthophora spicifera Inshore Population",
       x = expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
       y = expression("Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_classic()

### MODEL #1: NON-RECTANGULAR HYPERBOLA ##
# Nonlinear Least Squares regression of a non-rectangular hyperbola 
  # maybe derived from (Marshall & Biscoe, 1980), I think from Maya's MCAP/PCOM code

# set irradiance and photosynthetic rate for the test population
PAR <- as.numeric(as_in$light_value) # PAR = irradiance values = I
Pc <- as.numeric(as_in$mean_rate)  # Pc = mean metabolic rates (umol.g.hr)

# Set starting values for Model #1

# starting ma photosynthetic rate
Pmax_start<- max(Pc)

# starting quantum yield/alpha
AQY_start<-0.011

# starting value for dark respiration
Rd_start<-abs(min(Pc))

# starting value for curvature parameter
theta_start<-0.95

# Pmax_start = 19.15551 (start search for maximum gross photosynthetic rate near the highest observed rate)
# AQY_start = 0.011 (starting guess for alpha / apparent quantum yield; model will estimate final alpha
                    # high alpha -> alga responds strongly to small increases in light)
# Rd_start = 4.665729 (starting guess based on magnitude of the negative dark metabolic rate, Rdark = -Rd)
# theta_start = 0.95 (starting guess for curvature; model will estimate final theta)
  # theta controls bend btwn initial light-limited increase & the saturated part of the curve)
  # θ closer to 1 → sharper transition toward saturation, lower theta = more rounded transition

# Fit Model #1 to AS inshore
# fit the non-rectangular hyperbola model (maybe Marshall & Biscoe (1980))
as_in_M1 <- nls(Pc ~ (1/(2*theta)) * 
                  (AQY*PAR + Pmax -sqrt((AQY*PAR + Pmax)^2 
                                      - 4*AQY*theta*Pmax*PAR)) -Rd, #M1 = model 1 equation
# model is predicting what Pc (mean rate) should be

  # give nls the starting guesses defined above
  start = list(
    Pmax = Pmax_start, # starting guess for max gross photosynthetic rate
    AQY = AQY_start, # starting guess for alpha / apparent quantum yield
    Rd = Rd_start, # starting guess for dark respiration
    theta = theta_start)) # starting guess for curvature

summary(as_in_M1) # model results
#Parameters:
#Estimate Std. Error t value Pr(>|t|)    
# Pmax    24.190699   1.441037  16.787 2.85e-06 *** 
    # fitted maximum gross photosynthetic capacity, can be higher than 
    # highest observed net rate b/c model estimates gross capacity and then subtracts respiration

# AQY    0.046120   0.003962  11.640 2.42e-05 ***
    # alpha = initial slope of the PI curve at low irradiance
    # represents photochemical efficiency at low light

# Rd     4.320391   0.567994   7.606 0.000269 ***
    # fitted magnitude of dark respiration, the equation ends in - Rd, 
    # so the predicted net rate in darkness is about -4.32

# theta  0.962645   0.030635  31.423 6.90e-08 ***
    # fitted curvature parameter, this is close to 1 so fitted curve has 
    # sharp transition from the low-light rising portion toward saturation

coef(as_in_M1) # fitted parameter estimates
# Pmax         AQY          Rd       theta 
# 24.19069904  0.04612025  4.32039110  0.96264452 


## Extract fitted parameters from Model #1
coef_M1 <- coef(as_in_M1) # save model coefficients

# extract each fitted parameter
Pmax_M1 <- coef_M1["Pmax"] # maximum gross photosynthetic capacity

alpha_M1 <- coef_M1["AQY"] # alpha / apparent quantum yield

Rd_M1 <- coef_M1["Rd"] # dark respiration magnitude

theta_M1 <- coef_M1["theta"] # curvature parameter

# Pmax_M1 = 24.1907 
# alpha_M1 = 0.04612025 
# Rd_M1 = 4.320391 
# theta_M1 = 0.9626445 

## Plot fitted Model #1
# make a sequence of irradiance values to draw a smooth fitted curve
M1_curve <- tibble(PAR = seq(min(PAR), max(PAR), length.out = 200)) |> 
  # what does the model predict across the irradiance range

  # calculate predicted metabolic rate using fitted Model #1 parameters
  mutate(fitted_rate = (1/(2*theta_M1)) * 
           (alpha_M1*PAR + Pmax_M1 -sqrt((alpha_M1*PAR + Pmax_M1)^2 -
                                           4*alpha_M1*theta_M1*Pmax_M1*PAR)) -Rd_M1)

# plot observed population means and fitted Model #1 curve
ggplot() +
  geom_point(data = as_in, aes(x = light_value, y = mean_rate), size = 2) + # observed mean metabolic rates
  geom_errorbar(data = as_in,
                aes(x = light_value, ymin = mean_rate - se_rate, 
                    ymax = mean_rate + se_rate), width = 10) + # standard error around observed means
  geom_line(data = M1_curve,
            aes(x = PAR, y = fitted_rate), linewidth = 1) +  # fitted Model #1 curve
  scale_x_continuous(breaks = seq(0, 900, by = 100)) +
  labs(title = "Model 1: Acanthophora spicifera Inshore", 
       x = expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
       y = expression("Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_bw()

### Calculate AIC for Model 1
AIC_M1 <- AIC(as_in_M1)

AIC_M1
# [1] 25.01498

#AIC considers how well model fits the data and how complicated model is

# Save Model 1 results
M1_results <- tibble(species = "as", zone = "inshore",model = "non_rectangular_hyperbola",
                     Pmax = Pmax_M1, alpha = alpha_M1, Rd = Rd_M1,theta = theta_M1, AIC = AIC_M1)

M1_results
# species zone    model                      Pmax  alpha    Rd theta   AIC
#<chr>   <chr>   <chr>                     <dbl>  <dbl> <dbl> <dbl> <dbl>
# as     inshore non_rectangular_hyperbola  24.2 0.0461  4.32 0.963  25.0














##### PLOTTING CURVES #####

# diagnostic plot
PI_summary |> 
  ggplot(aes(x = light_value, y = mean_rate, color = zone, group = zone)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
                width = 10) +
  facet_wrap(~species, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 900, by = 100)) +
  labs(x = expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
       y = expression("Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"), color = "Zone") +
  theme_bw()

#shows mean response of all individuals within a species × zone at each experimental light step


#### organize data for individual PI curve fitting
PI_individual <- RespoR_Normalized |>
  filter(!is.na(zone)) |>
  select(species, full_species, zone, algae_id, light_level, light_value, umol.g.hr) |>
  arrange(species, zone, algae_id, light_value)

# diagnostic plot of the actual individual curves
PI_individual |>
  ggplot(aes(x = light_value, y = umol.g.hr, group = algae_id, color = zone)) +
  geom_point() +
  geom_line() +
  facet_wrap(~species, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 900, by = 100)) +
  labs(x = expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
       y = expression("Rate ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"), color = "Zone") +
  theme_bw()



   
  
# empty list to store output parameters
PI_outputs <- list()


#for loop to plot curves for all species
for(i in 1:nrow(PI_settings)) {
  
# get species and its starting values
  sp <- PI_settings$species[i]
  sp_name <- PI_settings$species_names[i]
  AQY_start <- PI_settings$AQY_start[i]
  theta_start <- PI_settings$theta_start[i]
  
  cat("\nFitting species:", sp, "\n")
 
  
# subset this species
  sp_resp <- RespoR_Normalized %>%
    filter(species == sp)
  
  
# calculate mean rate at each light value
  #sp_mean <- aggregate(umol.g.hr ~ light_value, data = sp_resp, FUN = mean)
  sp_mean<- sp_resp |> 
    group_by(light_value) |># group all measurments with the same irradiance
    summarise(umol.g.hr = mean(umol.g.hr, na.rm = TRUE), # calculate mean metabolic rate
              .groups = "drop")
  
  
# set PAR and Pc
  PAR <- as.numeric(sp_mean$light_value) # PAR = irradiance values
  Pc <- as.numeric(sp_mean$umol.g.hr)  # Pc = metabolic rates
  
  
# fit PI model (Marshall & Biscoe 1980)
# fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
  curve.nlslrc <- tryCatch(
    nls(Pc ~ (1/(2*theta)) *(AQY*PAR + Am - sqrt((AQY*PAR + Am)^2 - 4*AQY*theta*Am*PAR)) - Rd,
        start = list(
          Am = max(Pc), # Am = maximum gross photosynthetic rate
          AQY = AQY_start, # AQY = apparent quantum yield, or alpha,
          Rd = abs(min(Pc)), # Rd = dark respiration
          theta = theta_start)), # theta = curvature parameter
    
    error = function(e) {
      message("MODEL FAILED FOR ", sp, ": ", e$message)
      return(NULL)
    }
  )
  
  
  #### if model failed, skip to next species ####
  if(is.null(curve.nlslrc)) {
    next
  }
  
  
# extract model coefficients
  my.fit <- summary(curve.nlslrc)
  coef_fit <- coef(curve.nlslrc)
  
  Pmax.gross <- my.fit$parameters["Am", "Estimate"]
  AQY <- my.fit$parameters["AQY", "Estimate"]
  Rd <- my.fit$parameters["Rd", "Estimate"]
  theta <- my.fit$parameters["theta", "Estimate"]
  
  Ik <- Pmax.gross / AQY  
  Ic <- Rd / AQY
  Pmax.net <- Pmax.gross - Rd
  
  
  # save parameters 
  PI_outputs[[sp]] <- tibble(
    species = sp,
    Pg.max = Pmax.gross, # Max gross photosytnthetic rate
    Pn.max = Pmax.net, # Net photosynthetic rates
    Rdark = -Rd, # dark respiration
    alpha = AQY, # AQY (apparent quantum yield) alpha, initial slope/AQY
    theta = theta, # theta (curvature parameter)
    Ik = Ik, # Ik light saturation point
    Ic = Ic) # Ic light compensation point
  
  
  # make PDF plot
  pdf(here("Output", "PI", paste0(sp, "_PI_curve.pdf")),
      width = 7,height = 5)
  
  plot(PAR, Pc, xlab = "", ylab = "",
       xlim = c(0, max(PAR)), ylim = c(min(Pc) * 1.1, max(Pc) * 1.1),
       cex.lab = 0.8, cex.axis = 0.8, cex = 1,
       main = sp_name, font.main = 3, adj = 0.05)
  
  # fitted curve
  curve(
    (1/(2*coef_fit["theta"])) *
      (coef_fit["AQY"]*x + coef_fit["Am"] -
         sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 -
                4*coef_fit["AQY"]*
                coef_fit["theta"]*
                coef_fit["Am"]*x)) - coef_fit["Rd"],
    from = 0, to = max(PAR), lwd = 2, col = species_pal[sp], add = TRUE)
  
  # Ik line
  abline(v = Ik, col = species_pal[sp], lty = 2, lwd = 2)
  
  text(x = Ik + 75, y = 0, labels = paste0("Ik = ", round(Ik)))
  
  mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
        side = 1, line = 3.3, cex = 1)
  mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),
        side = 2, line = 2, cex = 1)
  
  dev.off()
}

# results combined into one table
PI_results <- bind_rows(PI_outputs)
print(PI_results, digits = 7)
   

 
  




