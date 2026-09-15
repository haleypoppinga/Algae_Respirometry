##################################################################################################
##################################################################################################
##################################################################################################
## CODE I TRIED DEVELOPING FOR COMPARING PI MODELS ####
## REFERENCE HERE TO USE SOME OF IT #####
##################################################################################################
##################################################################################################
##################################################################################################

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
ch.vol <- 482 #mL #of small chambers 12, volume of blank

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

RespoR <- tibble(.rows = length(filenames_final)*n_light_levels, #empty dataframe
                 file_id_csv = NA_character_,
                 sample_id = NA_character_,
                 intercept = NA_real_,
                 umol.L.sec = NA_real_,
                 temp_c = NA_real_,
                 light_level = NA_integer_,
                 light_value = NA_real_,
                 run_block = NA_character_,
                 zone = NA_character_) # fills in as forloop goes



# create directory for output folder
out_dir <- here("Output", "PI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



######### Create a for loop! ###############
############################################

###forloop#####

# to diagnose problems: set i = 1 and run through the loop line by line
# if run through loop and it fails: type "i" in the console to see what step (sample) forloop is failing on
for(i in 1:length(filenames_final)) { # loop each raw file?
  FRow <- as.numeric(which(Sample_Info$file_id_csv==filenames_final[i])) # stringsplit this renames our file
  # find row in meta data and match with filenames
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
  
  oxy_subsets <- Sample_Info[FRow,] %>% # set of dataframes, only used if had a file continuously used for multiple levels
    pmap(function(light_level, start_time, stop_time, ...) {
      data <- Respo.Data1  %>%
        filter(Time >= start_time & Time <= stop_time) %>%
        arrange(Time) %>%   # t_sec is difference between start and stop times
        mutate(t_sec = as.numeric(difftime(Time, first(Time), units = "secs"))) %>% #keep everything in seconds
        mutate(light_level = light_level) %>%
        filter(t_sec > 120) %>%     # drop first 2 min (120 s) (need at least 10 good minutes of data so dont cut out too much)
        filter(row_number() %% 10 == 0)    # keep every 10th row (why 10? this was because computer couldnt do every 1 second)
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
      y = expression(paste(' O'[2],' (',mu,'Viemol/L)')),
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
                      umol.L.sec = b1) # slope: y = mx +b done many times
           }) )
  # LoLinR (Local Linear Regressions) --> bootstrap method to go through data and objectively 
  
  
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
  left_join(Sample_Info) %>% #mutate(light_level = as.integer(light_level)),  # join with sample info which has all the meta data
  #by = c("file_id_csv", "light_level")) %>%
  mutate(Ch.Volume.mL = volume_mL, 
         Ch.Volume.L  = Ch.Volume.mL * 0.001, # convert from L to mL
         umol.sec = umol.L.sec * Ch.Volume.L) 

#Account for blank rate by sample run Block (if we do at least one blank per block)
# 12 samples, 2 of these are blanks
#View(RespoR)



####### Normalize the respo rates to the blanks ########

blank_rates <-RespoR2 %>% # compute blank rates first then left join instead here
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, run_block, light_dark) %>% # assuming 1 blank per run
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
         umol.g.hr = (umol.sec.corr * 3600) / wetweight_g, # normalize to wet weights (change when using AFDW)
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


# CHECK THE RAW DATA OF YOUR BLANKS
# if one is high and the other normal, take out the high one

#  basic plot of rates versus light before you make the real PI curve 
#basic_PI_plot <- RespoR_Normalized %>%
#  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id)) + # group by species, zone, etc.
#  geom_point()+
#  geom_line()+
#  facet_wrap(~species, scales = "free_y")

# CHECK FOR OUTLIERS (what does the blank look like, what does the raw data use)


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

# starting quantum yield/alpha (between 0-1)
AQY_start<-0.011

# starting value for dark respiration
Rd_start<-abs(min(Pc))

# starting value for curvature parameter (between 0-1)
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
# NOTES
# fitted function creates a line over raw data with error = "fitted_fun"
# add 95% confidence intervals around line of best fit
# IK is lowest possible light level of beginning of saturation
# not increasing photosynthesis, staying the same/stable
# find highest Ik and photoinhibition
# NEED ERROR FOR Ik
# fit it to each individual Ik for each individual or bootstrap 
# model table is important for my experiment (e.g. nls_model_coef_all)



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



##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
### CODE USED FOR FEBRUARY PI CURVES and WORKING THROUGH CODE IN AUGUST #####
##################################################################################################
##################################################################################################
##################################################################################################
###### Respo Code for Algae PI Curve ####### 
### Created by: Haley Poppinga, Maya Powell, Nyssa Silbiger
#### Last updated on: 2026-02-09

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
BioData <- read_csv(here("Data","Respo_Files","PI","Algae_Measurements_PI.csv"))

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
  select(file_id_csv, light_level, date, start_time, stop_time, species, algae_id, light_value, light_dark, blank, run_block, chamber_channel, volume_mL, wetweight_g)
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

n_light_levels<-8 # number of unique light levels

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
      #   data <- Respo.Data1  %>%
      #   filter(Time >= start_time & Time <= stop_time) %>%
      #   mutate(sec = row_number()) %>%# add an id for each row to help remove the first few mins
      #   mutate(Light_level = Light_level,
      #          sec = sec) %>%
      #   filter(sec > 60)  %>%# delete the first 2 mins of data assuming freq of 0.5 Hz
      #   mutate(row_number = row_number()) %>%
      #   filter(row_number %% 10 == 0) %>%  # keep every 10th row only to thin the data
      #   select(-row_number) %>%
      #   mutate(sec2 = row_number())  #update the row numbers
      # #return(subset)
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
  # index for the 8 light steps of this file
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

RespoR_Normalized_old <- RespoR2 %>% 
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  #summarise(blank.rate = mean(umol.sec, na.rm = TRUE)) %>% # if you have multiple blanks per run take the average
  #ungroup() %>% 
  #dplyr::select(Light_level, run_block, blank.rate, date) %>% # this is what I will use to join the blanks back with the raw data
  #right_join(RespoR2) %>% # join blanks with the respo data
  
  # join blank rates back onto ALL data, then correct, then normalize by wet weight
  filter(blank != 1) %>%  # remove blank rows from the "sample" dataset, # remove the Blank data
  left_join(blank_rates, by = c("light_level", "run_block", "light_dark")) %>% 
  mutate(umol.sec.corr   = umol.sec - blank.rate, # subtract the blank rates from the raw rates
         umol.g.hr = (umol.sec.corr * 3600) / wetweight_g,
         umol.g.hr_uncorr= (umol.sec * 3600) / wetweight_g) %>%
  dplyr::select(date, species, sample_id, algae_id, light_dark, run_block, wetweight_g, chamber_channel, #FIX .x AND .y COLUMNS
                temp_c, light_level, light_value, umol.sec, blank.rate, umol.sec.corr, 
                umol.g.hr, umol.g.hr_uncorr) #keep only what we need



#### making a df for just blank data for future use in plots #### 
Blank_only_old <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, light_value, run_block, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

#write_csv(RespoR_Normalized , here("Data","Respo_Files","PI","Respo_Algae_RNormalized_AllPIRates.csv"))  


## Plot the blanks across treatments to make sure nothing is funky
Blank_only_old %>%
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


# one label per line: take the max light point for each algae_id within each species
line_labels_old <- RespoR_Normalized_old %>%
  group_by(species, algae_id) %>%
  slice_max(light_value, n = 1, with_ties = FALSE) %>%
  ungroup()

basic_PI_plot <- RespoR_Normalized_old %>%
  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(data = line_labels_old, aes(label = algae_id), show.legend = FALSE,
                           size = 3, min.segment.length = 0) +
  facet_wrap(~species, scales = "free_y")
ggsave(here("Output","PI","basic_PI_plot.pdf"), basic_PI_plot, 
       width = 8, height = 8)

### run an nls model for PI curve and extract Ik for each species ###
beep(sound = 8, expr = NULL)








###############################################################################
###############################################################################
##### OLD CODE THAT WORKED WITH MAYA BUT HAVE EDITED SINCE ############ DON'T RUN
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
BioData <- read_csv(here("Data","Respo_Files","PI","Algae_Measurements_PI.csv"))

RespoMeta <- read_csv(here("Data","Respo_Files","PI","Algae_PI_meta.csv"))
#View(BioData)
#View(RespoMeta)
## try first with prelim fake data to make sure script runs
## then switch to real calculated data after getting volumes and weight and surface area


# join the data together
Sample_Info <- left_join(RespoMeta, BioData, by = c("sample_id", "algae_id"))
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

n_light_levels<-8 # number of unique light levels

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
      #   data <- Respo.Data1  %>%
      #   filter(Time >= start_time & Time <= stop_time) %>%
      #   mutate(sec = row_number()) %>%# add an id for each row to help remove the first few mins
      #   mutate(Light_level = Light_level,
      #          sec = sec) %>%
      #   filter(sec > 60)  %>%# delete the first 2 mins of data assuming freq of 0.5 Hz
      #   mutate(row_number = row_number()) %>%
      #   filter(row_number %% 10 == 0) %>%  # keep every 10th row only to thin the data
      #   select(-row_number) %>%
      #   mutate(sec2 = row_number())  #update the row numbers
      # #return(subset)
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
  # index for the 8 light steps of this file
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

RespoR <- read_csv(here("Data","Respo_Files","PI","Respo_Algae_R.csv"), show_col_types = FALSE) %>%
  #select(file_id_csv, light_level, intercept, umol.L.sec, temp_c) %>%
  mutate(light_level = as.integer(light_level))




######### Calculate Respiration rate ###############
############################################

#avoid duplicated columns after join (run_block/light_value already exist in RespoR)
RespoR2 <- RespoR %>%
  left_join(Sample_Info %>% mutate(light_level = as.integer(light_level)),
            by = c("file_id_csv", "light_level")) %>%
  mutate(Ch.Volume.mL = volume_mL,
         Ch.Volume.L  = Ch.Volume.mL * 0.001,
         umol.sec = umol.L.sec * Ch.Volume.L)

#Account for blank rate by sample run Block (if we do at least one blank per block)

#View(RespoR)

####### normalize the respo rates to the blanks ##### 
######### Blank correction (explicit join; keeps run_block intact) #########
# 1) calculate blank rate per (light_level x run_block x light_dark)
#blank_rates <- RespoR2 %>%
#  filter(blank == 1) %>%
#  group_by(light_level, run_block, light_dark) %>%
#  summarise(blank.rate = mean(umol.sec, na.rm = TRUE), .groups = "drop")

#RespoR2 <- RespoR %>%
  #drop_na(FileID_csv) %>% # drop NAs
#  left_join(Sample_Info) %>% # Join the raw respo calculations with the metadata
#  mutate(Ch.Volume.mL = volume_mL) %>% # 
#  mutate(Ch.Volume.L = Ch.Volume.mL * 0.001) %>% # mL to L conversion
#  mutate(umol.sec = umol.L.sec*Ch.Volume.L) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
#  mutate_if(sapply(., is.character), as.factor)  #convert character columns to factors

#Account for blank rate by sample run Block (if we do at least one blank per block)

#View(RespoR)

####### normalize the respo rates to the blanks ##### 

RespoR_Normalized <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, run_block.x, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE)) %>% # if you have multiple blanks per run take the average
  ungroup() %>% 
# join blank rates back onto ALL data, then correct, then normalize by wet weight
#RespoR_Normalized <- RespoR2 %>%
#  left_join(blank_rates, by = c("run_block", "light_dark")) %>%
  #dplyr::select(Light_level, run_block, blank.rate, date) %>% # this is what I will use to join the blanks back with the raw data
  right_join(RespoR2) %>% # join blanks with the respo data
  mutate(umol.sec.corr   = umol.sec - blank.rate, # subtract the blank rates from the raw rates
         umol.g.hr = (umol.sec.corr * 3600) / wetweight_g,
         umol.g.hr_uncorr= (umol.sec * 3600) / wetweight_g) %>%
  filter(blank != 1) %>%  # remove blank rows from the "sample" dataset, # remove the Blank data
  dplyr::select(date, species.x, sample_id.x, algae_id, light_dark, run_block.x, wetweight_g, chamber_channel, #FIX .x AND .y COLUMNS
                temp_c, light_level, light_value.x, umol.sec, blank.rate, umol.sec.corr, 
                umol.g.hr, umol.g.hr_uncorr) #keep only what we need


#######################
### making a df for just blank data for future use in plots ### 
#Blank_only <- blank_rates %>% 
#  left_join(Sample_Info %>% distinct(light_level, run_block, light_dark, light_value),
#            by = c("Light_level" = "light_level","run_block","light_dark"))
### making a df for just blank data for future use in plots ### 

Blank_only <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, light_value.x, run_block.x, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

write_csv(RespoR_Normalized , here("Data","Respo_Files","PI","Respo_Algae_RNormalized_AllPIRates.csv"))  


## Plot the blanks across treatments to make sure nothing is funky
Blank_only %>%
  ggplot(aes(x = light_value.x, blank.rate, group = interaction(run_block.x, light_dark))) +
  geom_point() +
  geom_line() +
  facet_wrap(~run_block.x)

#  basic plot of rates versus light before you make the real PI curve 
basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value.x, y = umol.g.hr, color = species.x, group = algae_id)) +
  geom_point()+
  geom_line()+
  facet_wrap(~species.x, scales = "free_y")
ggsave(here("Output","PI","basic_PI_plot.pdf"), basic_PI_plot)


### run an nls model for PI curve and extract Ik for each species ###
beep(sound = 8, expr = NULL)







###############################################################################
# Original PI Curve code last used June 2026

###### Respo Code for Algae PI Curve ####### 
### Created by: Haley Poppinga, Maya Powell, Nyssa Silbiger
#### Last updated on: 2026-
#### PI Update: 2026-06-01

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
BioData <- read_csv(here("Data","Respo_Files","PI","Algae_Measurements_PI.csv"))

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
  select(file_id_csv, light_level, date, start_time, stop_time, species, algae_id, light_value, light_dark, blank, run_block, chamber_channel, volume_mL, wetweight_g)
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

n_light_levels<-8 # number of unique light levels

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
      #   data <- Respo.Data1  %>%
      #   filter(Time >= start_time & Time <= stop_time) %>%
      #   mutate(sec = row_number()) %>%# add an id for each row to help remove the first few mins
      #   mutate(Light_level = Light_level,
      #          sec = sec) %>%
      #   filter(sec > 60)  %>%# delete the first 2 mins of data assuming freq of 0.5 Hz
      #   mutate(row_number = row_number()) %>%
      #   filter(row_number %% 10 == 0) %>%  # keep every 10th row only to thin the data
      #   select(-row_number) %>%
      #   mutate(sec2 = row_number())  #update the row numbers
      # #return(subset)
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
  # index for the 8 light steps of this file
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
  filter(!(sample_id == "AV02_PI_RUN1"),
         !(sample_id == "AV03_PI_RUN1" & light_level == 3),
         !(sample_id == "AV05_PI_RUN3"),
         !(sample_id == "DA01_PI_RUN1" & light_level == 5),
         !(sample_id == "DA02_PI_RUN3"),
         !(sample_id == "GS02_PI_RUN1"),
         !(sample_id == "HD02_PI_RUN3" & light_level == 5),
         !(sample_id == "SF05_PI_RUN3" & light_level == 3 & light_level == 4)) %>% 
  
  # join blank rates back onto ALL data, then correct, then normalize by wet weight
  filter(blank != 1) %>%  # remove blank rows from the "sample" dataset, # remove the Blank data
  left_join(blank_rates, by = c("light_level", "run_block", "light_dark")) %>% 
  mutate(umol.sec.corr   = umol.sec - blank.rate, # subtract the blank rates from the raw rates
         umol.g.hr = (umol.sec.corr * 3600) / wetweight_g,
         umol.g.hr_uncorr= (umol.sec * 3600) / wetweight_g) %>%
  dplyr::select(date, species, sample_id, algae_id, light_dark, run_block, wetweight_g, chamber_channel, #FIX .x AND .y COLUMNS
                temp_c, light_level, light_value, umol.sec, blank.rate, umol.sec.corr, 
                umol.g.hr, umol.g.hr_uncorr) #keep only what we need

###### Filter out bad samples #############
#AS02_PI_RUN2_8 (stressed)
#AV02_PI_RUN1
#AV03_PI_RUN1_3 AND 4 LOOK THE SAME
#AV05_PI_RUN3
#DA01_PI_RUN1_5
#DA02_PI_RUN3
#DA04_PI_RUN1_7 (stressed)
#GS02_PI_RUN1 BEFORE 7
#GS04_PI_RUN1_8 (stressed)
#HD02_PI_RUN3_5
#SF01_PI_RUN1_8 (stressed)
#SF03_PI_RUN2_7 AND 8 (stressed)
#SF05_PI_RUN3_3 AND 4


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
species_pal <- c(av  = "#1B9E77",  # green
                 hd  = "#66A61E",  # green
                 cs  = "#A6D854",  # green
                 as  = "#E41A1C",  # red
                 gs  = "#FB6A4A",  # red
                 sf  = "#380000",  # red
                 da  = "#8C510A")  # brown

algae_basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(data = line_labels, aes(label = algae_id), show.legend = FALSE,
                           size = 3, min.segment.length = 0, segment.color = "black") +
  scale_color_manual(values = species_pal, drop = FALSE) +
  facet_wrap(~species, scales = "free_y")

ggsave(here("Output","PI","Algae_Basic_PI_plot.pdf"), algae_basic_PI_plot, width = 8,
       height = 8)

### run an nls model for PI curve and extract Ik for each species ###

##### PLOTTING CURVES #####
##### Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)

#Plot curves
av_resp <- RespoR_Normalized %>% filter(species == "av")
as_resp <- RespoR_Normalized %>% filter(species == "as")
gs_resp <- RespoR_Normalized %>% filter(species == "gs")
cs_resp <- RespoR_Normalized %>% filter(species == "cs")
sf_resp <- RespoR_Normalized %>% filter(species == "sf")
da_resp <- RespoR_Normalized %>% filter(species == "da")

#algae.data <- read.table("Data/Respo_Files/PI/Respo_Algae_RNormalized_AllPIRates.csv", header=TRUE, sep=",")

# means and se for each species
av.mean <- aggregate(umol.g.hr ~ light_value, data = av_resp, FUN=mean)
av.se <- aggregate(umol.g.hr ~ light_value, data = av_resp, FUN=std.error)

as.mean <- aggregate(umol.g.hr ~ light_value, data = as_resp, FUN=mean)
as.se <- aggregate(umol.g.hr ~ light_value, data = as_resp, FUN=std.error)

gs.mean <- aggregate(umol.g.hr ~ light_value, data = gs_resp, FUN=mean)
gs.se <- aggregate(umol.g.hr ~ light_value, data = gs_resp, FUN=std.error)

cs.mean <- aggregate(umol.g.hr ~ light_value, data = cs_resp, FUN=mean)
cs.se <- aggregate(umol.g.hr ~ light_value, data = cs_resp, FUN=std.error)

sf.mean <- aggregate(umol.g.hr ~ light_value, data = sf_resp, FUN=mean)
sf.se <- aggregate(umol.g.hr ~ light_value, data = sf_resp, FUN=std.error)

da.mean <- aggregate(umol.g.hr ~ light_value, data = da_resp, FUN=mean)
da.se <- aggregate(umol.g.hr ~ light_value, data = da_resp, FUN=std.error)


### av Data ###
PAR <- as.numeric(av.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(av.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "av_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Avrainvillea lacerata", font.main = 3, adj = 0.05) #set plot info
abline(v = 369.02563713, col = "#06402B", lty = 2, lwd = 2) #green dashed line at Ik value
text(x =450, y = 0, labels = "Ik = 369") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.005, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for corals
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

#Amax (max gross photosytnthetic rate)
#Pmax.gross <- my.fit$parameters[1]
Pmax.gross <- my.fit$parameters["Am", "Estimate"]

#AQY (apparent quantum yield) alpha
#AQY <- my.fit$parameters[2]
AQY <- my.fit$parameters["AQY", "Estimate"]

#Rd (dark respiration)
#Rd <- my.fit$parameters[3]
Rd<- my.fit$parameters["Rd", "Estimate"]
theta <- my.fit$parameters["theta", "Estimate"]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
av.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(av.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(av.PI.Output)

#Data for av
# Pg.max   7.11017429
# Pn.max   5.76583832
# Rdark   -1.34433597
# alpha    0.01926743
# Ik     369.02563713
# Ic      69.77247223

#dev.off()
#dev.new()

#### as Data #####
PAR <- as.numeric(as.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(as.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "as_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Acanthophora spicifera", font.main = 3, adj = 0.05) #set plot info
abline(v = 519.35793869, col = "red", lty = 2, lwd = 2) #Red dashed line at Ik value
text(x =600, y = 0, labels = "Ik = 519") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.015, Rd=abs(min(Pc)),theta=0.95)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
as.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(as.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(as.PI.Output)

#Data for as
# Pg.max  33.61612753
# Pn.max  29.70310777
# Rdark   -3.91301976
# alpha    0.06472632
# Ik     519.35793869
# Ic      60.45484796

#### gs Data #####
PAR <- as.numeric(gs.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(gs.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "gs_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Gracilaria salicornia", font.main = 3, adj = 0.05) #set plot info
abline(v = 416.83242652, col = "#FFBF00", lty = 2, lwd = 2) #yellow dashed line at Ik value
text(x =550, y = 0, labels = "Ik = 417") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.02, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
gs.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(gs.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(gs.PI.Output)

#Data for gs
# Pg.max  12.49767214
# Pn.max  11.13203615
# Rdark   -1.36563599
# alpha    0.02998249
# Ik     416.83242652
# Ic      45.54779152



#### cs Data #####
PAR <- as.numeric(cs.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(cs.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "cs_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Caulerpa sertularioides", font.main = 3, adj = 0.05) #set plot info
abline(v = 245.2334301, col = "#32CD32", lty = 2, lwd = 2) #light green dashed line at Ik value
text(x =350, y = 0, labels = "Ik = 245") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.025, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
cs.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(cs.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(cs.PI.Output)

#Data for cs
# Pg.max  39.5763171
# Pn.max  35.4103641
# Rdark   -4.1659530
# alpha    0.1613822
# Ik     245.2334301
# Ic      25.8141996


#### sf Data #####
PAR <- as.numeric(sf.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(sf.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "sf_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Spyridia filamentosa", font.main = 3, adj = 0.05) #set plot info
abline(v = 356.2966178, col = "#950606", lty = 2, lwd = 2) #yellow dashed line at Ik value
text(x =450, y = 0, labels = "Ik = 356") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   #start=list(Am=(max(Pc)), AQY=0.01, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
                   start = list(Am =(max(Pc)), AQY = 0.0098, Rd = 3.156, theta = 0.88))
my.fit <- summary(curve.nlslrc) #summary of model fit


#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)

dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
sf.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(sf.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(sf.PI.Output)

#Data for sf
# Pg.max  28.0200573
# Pn.max  24.8640286
# Rdark   -3.1560287
# alpha    0.0786425
# Ik     356.2966178
# Ic      40.1313362



#### da Data #####
PAR <- as.numeric(da.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(da.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "da_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Dictyota acutiloba", font.main = 3, adj = 0.05) #set plot info
abline(v = 284.1703581, col = "#9B7A01", lty = 2, lwd = 2) #yellow dashed line at Ik value
text(x =350, y = 0, labels = "Ik = 284") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.08, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
da.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(da.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(da.PI.Output)

#Data for da
# Pg.max  86.1895094
# Pn.max  74.9572540
# Rdark  -11.2322554
# alpha    0.3033023
# Ik     284.1703581
# Ic      37.0332081


beep(sound = 8, expr = NULL)



###############################################################################
###############################################################################
###############################################################################
###############################################################################
################### old code to run forloop of models on all species ###########
# Code from 2026-08-10 run with two different models using a for loop instead of above
# insert after aggregate code

# starting values for each species
PI_settings <- tibble(
  species = c("av", "as", "gs", "cs", "sf", "da", "ds", "dc", "hd"),
  species_names = c(
    av = "Avrainvillea lacerata",
    as = "Acanthophora spicifera",
    gs = "Gracilaria salicornia",
    cs = "Caulerpa sertularioides",
    sf = "Spyridia filamentosa",
    da = "Dictyota acutiloba",
    ds = "Dictyota sandvicensis",
    dc = "Dictyosphaeria cavernosa",
    hd = "Halimeda discoidea"),
  
  # starting slope value 
  AQY_start = c( 
    0.004,  # av (good)
    0.014,  # as (good)
    0.020,  # gs (good)
    0.025,  # cs (good)
    0.020,  # sf (showing inhibition, cant fit)
    0.070,  # da (showing inhibition, cant fit)
    0.020,  # ds (good)
    0.010,  # dc (good)
    0.004), # hd (good)
  
  # theta values closer to 1 make transition toward saturation sharper, lower values more rounded
  theta_start = c( 
    0.90,   # av (good)
    0.95,   # as (good)
    0.90,   # gs (good)
    0.90,   # cs (good)
    0.60,   # sf (showing inhibition, cant fit
    0.60,   # da (showing inhibition, cant fit)
    0.90,   # ds (good)
    0.90,   # dc (good)
    0.90))  # hd (good)

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
  
  
  #### subset this species ####
  sp_resp <- RespoR_Normalized %>%
    filter(species == sp)
  
  
  #### calculate mean rate at each light value ####
  sp_mean <- aggregate(
    umol.g.hr ~ light_value,
    data = sp_resp,
    FUN = mean
  )
  
  
  #### set PAR and Pc ####
  PAR <- as.numeric(sp_mean$light_value)
  Pc <- as.numeric(sp_mean$umol.g.hr)
  
  
  #### fit PI model ####
  curve.nlslrc <- tryCatch(
    
    nls(Pc ~ (1/(2*theta)) *
          (AQY*PAR + Am -
             sqrt((AQY*PAR + Am)^2 -
                    4*AQY*theta*Am*PAR)) - Rd,
        start = list(
          Am = max(Pc),
          AQY = AQY_start,
          Rd = abs(min(Pc)),
          theta = theta_start)),
    
    error = function(e) {
      message("MODEL FAILED FOR ", sp, ": ", e$message)
      return(NULL)
    }
  )
  
  
  #### if model failed, skip to next species ####
  if(is.null(curve.nlslrc)) {
    next
  }
  
  
  #### extract model coefficients ####
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
    Pg.max = Pmax.gross,
    Pn.max = Pmax.net,
    Rdark = -Rd,
    alpha = AQY,
    theta = theta,
    Ik = Ik,
    Ic = Ic)
  
  
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

#species Pg.max Pn.max  Rdark  alpha theta    Ik    Ic
#<chr>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
#1 av        4.07   2.89 -1.18  0.0111 0.906  367. 107. 
#2 as       24.1   20.5  -3.53  0.0417 0.975  578.  84.8
#3 gs       11.5    9.72 -1.79  0.0233 0.886  495.  76.8
#4 cs       33.9   28.7  -5.25  0.112  0.738  302.  46.8
#5 ds       45.6   39.2  -6.34  0.160  0.942  285.  39.6
#6 dc        6.10   5.11 -0.984 0.0220 0.872  277.  44.8
#7 hd        6.82   5.97 -0.848 0.0130 0.989  523.  65.0

#beep(sound = 8, expr = NULL)



#### sf photoinhibition model ####
# Platt, Gallegos & Harrison (1980) photoinhibition model
# specifically developed to describe PI curves from the initial light-limited increase through the photoinhibited range
# Platt formulation is P=Ps(1−e^(−αI/Ps))e^−βI/P --> 2nd exponential allows curve to turn downward

#Jasby Platt? compare AIC (lowest AIC is best fit), calculate Ik

#Ps     = photosynthetic capacity parameter
#alpha  = initial low-light slope (your AQY-type parameter)
#beta   = photoinhibition parameter
#Rd     = dark respiration

PAR <- as.numeric(sf.mean$light_value)
Pc  <- as.numeric(sf.mean$umol.g.hr)

sf_dat <- data.frame(PAR, Pc)

# rough values from the data to define reasonable search ranges
Rd_guess <- abs(Pc[which.min(PAR)])
Ps_guess <- max(Pc) + Rd_guess

curve.sf <- nls.multstart::nls_multstart(
  Pc ~ Ps *  (1 - exp(-alpha * PAR / Ps)) *exp(-beta * PAR / Ps) - Rd,
  data = sf_dat, iter = 500,
  start_lower = c(Ps    = Ps_guess * 0.5,
                  alpha = 0.001,beta  = 0.0001, Rd= Rd_guess * 0.5),
  start_upper = c(Ps    = Ps_guess * 3,
                  alpha = 0.5, beta  = 0.5, Rd = Rd_guess * 2),
  lower = c(Ps = 0, alpha = 0,beta  = 0,Rd    = 0),
  supp_errors = "Y")

summary(curve.sf)
coef(curve.sf)
coef_sf <- coef(curve.sf)

pdf(here("Output", "PI", "sf_PI_curve_photoinhibition.pdf"),
    width = 7,height = 5)

plot(PAR, Pc, xlab = "",ylab = "", 
     xlim = c(0, 900), ylim = c(min(Pc) * 1.1, max(Pc) * 1.1),
     main = "Spyridia filamentosa", font.main = 3,adj = 0.05)

curve(coef_sf["Ps"] * (1 - exp(-coef_sf["alpha"] * x / coef_sf["Ps"])) *
        exp(-coef_sf["beta"] * x / coef_sf["Ps"]) -coef_sf["Rd"],
      from = 0,to = 900,lwd = 2,col = species_pal["sf"],add = TRUE)

mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
      side = 1, line = 3.3)
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side = 2, line = 2)
dev.off()


#### da photoinhibition model ####
# Platt, Gallegos & Harrison (1980) photoinhibition model
# specifically developed to describe PI curves from the initial light-limited increase through the photoinhibited range
# Platt formulation is P=Ps(1−e^(−αI/Ps))e^−βI/P --> 2nd exponential allows curve to turn downward

#Ps     = photosynthetic capacity parameter
#alpha  = initial low-light slope (your AQY-type parameter)
#beta   = photoinhibition parameter
#Rd     = dark respiration

PAR <- as.numeric(da.mean$light_value)
Pc  <- as.numeric(da.mean$umol.g.hr)

da_dat <- data.frame(PAR, Pc)

# rough values from the data to define reasonable search ranges
Rd_guess <- abs(Pc[which.min(PAR)])
Ps_guess <- max(Pc) + Rd_guess

curve.da <- nls.multstart::nls_multstart(
  Pc ~ Ps *  (1 - exp(-alpha * PAR / Ps)) *exp(-beta * PAR / Ps) - Rd,
  data = da_dat, iter = 500,
  start_lower = c(Ps    = Ps_guess * 0.5,
                  alpha = 0.001,beta  = 0.0001, Rd= Rd_guess * 0.5),
  start_upper = c(Ps    = Ps_guess * 3,
                  alpha = 0.5, beta  = 0.5, Rd = Rd_guess * 2),
  lower = c(Ps = 0, alpha = 0,beta  = 0,Rd    = 0),
  supp_errors = "Y")

summary(curve.da)
coef(curve.da)
coef_da <- coef(curve.da)

pdf(here("Output", "PI", "da_PI_curve_photoinhibition.pdf"),
    width = 7,height = 5)

plot(PAR, Pc, xlab = "",ylab = "", 
     xlim = c(0, 900), ylim = c(min(Pc) * 1.1, max(Pc) * 1.1),
     main = "Dictyota acutiloba", font.main = 3,adj = 0.05)

curve(coef_da["Ps"] * (1 - exp(-coef_da["alpha"] * x / coef_da["Ps"])) *
        exp(-coef_da["beta"] * x / coef_da["Ps"]) -coef_da["Rd"],
      from = 0,to = 900,lwd = 2,col = species_pal["da"],add = TRUE)

mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
      side = 1, line = 3.3)
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side = 2, line = 2)
dev.off()


###############################################################################
###############################################################################
###############################################################################
########## same code as above but slightly updated incase i need to go back###
#Plot curves
#algae.data <- read.table("Data/Respo_Files/PI/Respo_Algae_RNormalized_AllPIRates.csv", header=TRUE, sep=",")

# means and se for each species inshore and offshore for each light level
PI_summary <- RespoR_Normalized |> 
  filter(!is.na(zone)) |> #filter inshore/offshore
  group_by(species, full_species, zone, light_level) |>  
  summarise(light_value = mean(light_value, na.rm = TRUE),
            mean_rate = mean(umol.g.hr, na.rm = TRUE), #calc mean rates
            n = sum(!is.na(umol.g.hr)),
            se_rate = sd(umol.g.hr, na.rm = TRUE) / sqrt(n), #error
            n_individuals = n_distinct(algae_id),
            .groups = "drop")

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


# starting values for each species
PI_settings <- tibble(species = c("av", "as", "gs", "cs", "sf", "da", "ds", "dc", "hd"),
                      species_names = c(
                        av = "Avrainvillea lacerata",
                        as = "Acanthophora spicifera",
                        gs = "Gracilaria salicornia",
                        cs = "Caulerpa sertularioides",
                        sf = "Spyridia filamentosa",
                        da = "Dictyota acutiloba",
                        ds = "Dictyota sandvicensis",
                        dc = "Dictyosphaeria cavernosa",
                        hd = "Halimeda discoidea"),
                      
                      # starting slope value           
                      AQY_start = c(
                        0.004,  # av (good)
                        0.014,  # as (good)
                        0.020,  # gs (good)
                        0.025,  # cs (good)
                        0.020,  # sf (showing inhibition, cant fit)
                        0.070,  # da (showing inhibition, cant fit)
                        0.020,  # ds (good)
                        0.010,  # dc (good)
                        0.004), # hd (good)
                      
                      # theta values closer to 1 make transition toward saturation sharper, lower values more rounded
                      theta_start = c( 
                        0.90,   # av (good)
                        0.95,   # as (good)
                        0.90,   # gs (good)
                        0.90,   # cs (good)
                        0.60,   # sf (showing inhibition, cant fit
                        0.60,   # da (showing inhibition, cant fit)
                        0.90,   # ds (good)
                        0.90,   # dc (good)
                        0.90))  # hd (good)

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


###############################################################################
###############################################################################
###############################################################################
################### old code to calculate each individually ###################
### av Data ###
PAR <- as.numeric(av.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(av.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "av_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Avrainvillea lacerata", font.main = 3, adj = 0.05) #set plot info
abline(v = 369.02563713, col = "#06402B", lty = 2, lwd = 2) #green dashed line at Ik value
text(x =450, y = 0, labels = "Ik = 369") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.005, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for corals
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

#Amax (max gross photosytnthetic rate)
#Pmax.gross <- my.fit$parameters[1]
Pmax.gross <- my.fit$parameters["Am", "Estimate"]

#AQY (apparent quantum yield) alpha
#AQY <- my.fit$parameters[2]
AQY <- my.fit$parameters["AQY", "Estimate"]

#Rd (dark respiration)
#Rd <- my.fit$parameters[3]
Rd<- my.fit$parameters["Rd", "Estimate"]
theta <- my.fit$parameters["theta", "Estimate"]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
av.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(av.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(av.PI.Output)

#Data for av
# Pg.max   7.11017429
# Pn.max   5.76583832
# Rdark   -1.34433597
# alpha    0.01926743
# Ik     369.02563713
# Ic      69.77247223

#dev.off()
#dev.new()

#### as Data #####
PAR <- as.numeric(as.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(as.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "as_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Acanthophora spicifera", font.main = 3, adj = 0.05) #set plot info
abline(v = 519.35793869, col = "red", lty = 2, lwd = 2) #Red dashed line at Ik value
text(x =600, y = 0, labels = "Ik = 519") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.015, Rd=abs(min(Pc)),theta=0.95)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
as.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(as.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(as.PI.Output)

#Data for as
# Pg.max  33.61612753
# Pn.max  29.70310777
# Rdark   -3.91301976
# alpha    0.06472632
# Ik     519.35793869
# Ic      60.45484796

#### gs Data #####
PAR <- as.numeric(gs.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(gs.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "gs_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Gracilaria salicornia", font.main = 3, adj = 0.05) #set plot info
abline(v = 416.83242652, col = "#FFBF00", lty = 2, lwd = 2) #yellow dashed line at Ik value
text(x =550, y = 0, labels = "Ik = 417") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.02, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
gs.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(gs.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(gs.PI.Output)

#Data for gs
# Pg.max  12.49767214
# Pn.max  11.13203615
# Rdark   -1.36563599
# alpha    0.02998249
# Ik     416.83242652
# Ic      45.54779152



#### cs Data #####
PAR <- as.numeric(cs.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(cs.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "cs_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Caulerpa sertularioides", font.main = 3, adj = 0.05) #set plot info
abline(v = 245.2334301, col = "#32CD32", lty = 2, lwd = 2) #light green dashed line at Ik value
text(x =350, y = 0, labels = "Ik = 245") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.025, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
cs.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(cs.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(cs.PI.Output)

#Data for cs
# Pg.max  39.5763171
# Pn.max  35.4103641
# Rdark   -4.1659530
# alpha    0.1613822
# Ik     245.2334301
# Ic      25.8141996


#### sf Data #####
PAR <- as.numeric(sf.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(sf.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "sf_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Spyridia filamentosa", font.main = 3, adj = 0.05) #set plot info
abline(v = 356.2966178, col = "#950606", lty = 2, lwd = 2) #yellow dashed line at Ik value
text(x =450, y = 0, labels = "Ik = 356") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   #start=list(Am=(max(Pc)), AQY=0.01, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
                   start = list(Am =(max(Pc)), AQY = 0.0098, Rd = 3.156, theta = 0.88))
my.fit <- summary(curve.nlslrc) #summary of model fit


#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)

dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
sf.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(sf.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(sf.PI.Output)

#Data for sf
# Pg.max  28.0200573
# Pn.max  24.8640286
# Rdark   -3.1560287
# alpha    0.0786425
# Ik     356.2966178
# Ic      40.1313362



#### da Data #####
PAR <- as.numeric(da.mean$light_value) #PAR = irradiance values
Pc <- as.numeric(da.mean$umol.g.hr) #Pc = metabolic rates
pdf(here("Output", "PI", "da_PI_curve.pdf"), width = 7, height = 5)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(min(Pc) * 1.1, max(Pc) * 1.1), 
     cex.lab=0.8,cex.axis=0.8,cex=1, main="Dictyota acutiloba", font.main = 3, adj = 0.05) #set plot info
abline(v = 284.1703581, col = "#9B7A01", lty = 2, lwd = 2) #yellow dashed line at Ik value
text(x =350, y = 0, labels = "Ik = 284") #add Ik label to line
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
#Am = maximum gross photosynthetic rate, AQY = apparent quantum yield, or alpha, Rd = dark respiration, theta = curvature parameter
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                   start=list(Am=(max(Pc)), AQY=0.08, Rd=abs(min(Pc)),theta=0.9)) #changed to smaller AQY value and higher theta than for previous AV species
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
coef_fit <- coef(curve.nlslrc) #extract coefficients first
curve((1/(2*coef_fit["theta"])) *(coef_fit["AQY"]*x + coef_fit["Am"] -
                                    sqrt((coef_fit["AQY"]*x + coef_fit["Am"])^2 - 4*coef_fit["AQY"]*coef_fit["theta"]*coef_fit["Am"]*x)) - coef_fit["Rd"],
      from = 0, to = max(PAR), lwd = 2, col = "blue",add = TRUE)
dev.off()

Pmax.gross <- my.fit$parameters["Am", "Estimate"] #Amax (max gross photosynthetic rate) #instead of by row
AQY <- my.fit$parameters["AQY", "Estimate"] #AQY (apparent quantum yield) alpha
Rd<- my.fit$parameters["Rd", "Estimate"] #Rd (dark respiration)
theta <- my.fit$parameters["theta", "Estimate"]
Ik <- Pmax.gross/AQY # Ik light saturation point
Ic <- Rd/AQY # Ic light compensation point
Pmax.net <- Pmax.gross - Rd # Net photosynthetic rates

#output parameters into a table
da.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(da.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
print(da.PI.Output)

#Data for da
# Pg.max  86.1895094
# Pn.max  74.9572540
# Rdark  -11.2322554
# alpha    0.3033023
# Ik     284.1703581
# Ic      37.0332081






#################################################################
#################################################################
#### Respo code on 9/13/26 before editing forloop again ####
#################################################################
#################################################################
#################################################################
###### Respo Code for PI Curves ####### 
### Created by: Nyssa Silbiger, Maya Powell, Haley Poppinga
#### Last updated on: 2026-09-14

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
         umol.g.hr = (umol.sec.corr*3600)/wetweight_g, # change to AFDW when obtained
         umol.g.hr_uncorr = (umol.sec*3600)/wetweight_g) %>%
  
  # join blank rates back onto ALL data, then correct, then normalize by wet weight
  filter(blank!=1) %>% # remove the Blank data
  #ungroup() %>%
  dplyr::select(date,species, full_species, zone, sample_id, algae_id, light_dark,
                run_block, wetweight_g, chamber_channel, temp_c, light_level, light_value,
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
Blank_only %>%
  ggplot(aes(x = light_value, blank.rate)) +
  geom_point()
#ggsave(here("Output","PI","blank_only.pdf"), Blank_only)

#  basic plot of rates versus light before you make the real PI curve 
basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id))+
  geom_point()+
  geom_line()+
  facet_wrap(~species, scales = "free")
basic_PI_plot
ggsave(here("Output","PI","basic_PI_plot.pdf"), basic_PI_plot, width = 10, height = 10)

#set up species colors
species_pal <- c(av  = "#033500",  # green
                 hd  = "#66A61E",  # green
                 cs  = "#20B2AA",  # green
                 dc  = "#239665",  # green
                 as  = "#E41A1C",  # red
                 gs  = "#FFA500",  # red
                 sf  = "#7C0A02",  # red
                 da  = "#A57A54",  # brown
                 ds  = "#A89F51")  # brown

# another asic plot thats maybe nicer
# one label per line, take the max light point for each algae_id within each species
line_labels <- RespoR_Normalized %>%
  group_by(species, algae_id) %>%
  slice_max(light_value, n = 1, with_ties = FALSE) %>% # add label to each id
  ungroup()

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

### Prep individual PI data 
# df for every individual PI curve
PI_individuals <- RespoR_Normalized |>
  filter(!is.na(zone), !is.na(algae_id)) |> # filter by zone and algae individual
  distinct(species, species_name, zone, algae_id)

PI_individuals

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
#Pmax = "Am" = (max photosynthetic rate) should be positive based on range of umol.cm2.hr - here ~0-2
#AQY (apparent quantum yield) should be positive and between 0-1 (often much smaller)
#Rd should be subtracted and usually between -2-0 but can adjust if larger
#theta is between 0 and 1 


### make function uded to plot Model 1 curve
make_nrh_fun <- function(Pmax, AQY, Rd, theta) {
  force(Pmax); force(AQY); force(Rd); force(theta)
  function(x) nrh_model(x, Pmax, AQY, Rd, theta)
}

# create plot list and a coefficient table list for Model 1
plot_list_Model1 <- list() # for individual plots
coef_list_Model1 <- list() # coefficient table for each individual
curve_list_Model1 <- list() # predicted curve values for each individual
raw_list_Model1 <- list() # raw observed points for each individual

# MODEL 1 forloop: to calculate ik and plot
# fit nrh model to each species x zone population

### make a forloop
#Set i = "AS" or whatever you want to look at for each species_name as you go through
for (i in unique(RespoR_Normalized$algae_id)) {
  data <- RespoR_Normalized |> 
    filter(algae_id == i) |>  # subset individual's PI curve
    arrange(light_value) # or light level?
  
  
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
  ))
  
  # need to add this error part for the species that go into photoinhibition
  # if one individual fails, continue to next one
  error = function(e) {
    message("MODEL 1 FAILED FOR ", i, ": ", e$message)
    return(NULL)
  }
  
  
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
  
  # make a smooth light sequence for the fitted curve
  curve_df <- tibble(light_value = seq(0, 900, by = 5)) |>
    mutate(pred_rate = fitted_fun(light_value), # predicted rate
           species = sp, species_name = sp_name, zone = z, algae_id = i)
  
  # save predicted curve values
  curve_list_Model1[[i]] <- curve_df
  
  # save raw observed data for this individual
  raw_list_Model1[[i]] <- data |>
    select(species, species_name, zone, algae_id, light_value, umol.g.hr)
  
  # plot individual sample's PI curve
  basic_PI_plot <- data %>%
    ggplot(aes(x = light_value, y = umol.g.hr)) +
    geom_point(color = species_pal[sp], size = 2) +
    #geom_line()+
    stat_function(fun = fitted_fun, linewidth = 1, color = species_pal[sp]) + #fitted function
    geom_vline(xintercept = ik_val, color = species_pal[sp], lty = 2, linewidth = 1) + # adds Ik line
    annotate("text", x = ik_val + 40, y = 0, label = paste0("Ik = ", round(ik_val, 1))) + # labels Ik line
    #you can change ik plus however much you need to get the spacing correct
    scale_x_continuous(breaks = seq(0, 900, by = 100)) + # light value breaks
    labs(x = expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),
         y = expression(Rate*" ("*mu*"mol "*O[2]*" "*g^-1*h^-1*")"),
         title = paste0(sp_name, " - ", z, " - ", i))+ # add zone to species name
    theme_classic(base_size = 15) +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))
  
  plot_list_Model1[[i]] <- basic_PI_plot # save individual sample's plot and coefficients
  
}

# combine all individual Model 1 results into one dataframe
Model1_individuals <- bind_rows(coef_list_Model1)

# combine all fitted curve values into one dataframe
Model1_curve_data <- bind_rows(curve_list_Model1)

# combine all raw observed points into one dataframe
Model1_raw_data <- bind_rows(raw_list_Model1)

# look at coefficient output
View(Model1_individuals)

# save coefficient table
write.csv(Model1_individuals,
          here("Output/PI/Model1_individual_coefficients.csv"),
          row.names = FALSE)


all_species_Model1_plots <- ggarrange(
  plotlist = species_plot_list_Model1,
  ncol = 3, nrow = 3
)

all_species_Model1_plots

ggsave(here("Output/PI/Model1_all_species_grouped_plots.pdf"),
       all_species_Model1_plots,
       width = 16, height = 12)
# what is this below?
all_nls_PI_plots <- ggarrange(plotlist = plot_list_Model1, ncol = 5, nrow = 2)
all_nls_PI_plots
ggsave(here("Output/PI/ikonik_publishable_plot.pdf"), all_nls_PI_plots, width = 20, height = 8)
View(Model1_individuals)
write.csv(Model1_individuals, here("Output/PI/Model1_individual_coefficients.csv"), row.names = FALSE)





############################################################

###### DEFINE MODEL #######

### MODEL 2: Platt et al. (1980) Photoinhibition Model ###
platt_model <- function(light_value, Ps, AQY, beta, Rd) {
  Ps *(1 - exp(-AQY * light_value / Ps)) *exp(-beta * light_value / Ps) - Rd
}

# Another possible model but wont try
### MODEL 3: Jassby & Platt (1976) Hyperbolic Tangent ###
#jp_model <- function(light_value, Pmax, AQY, Rd) {
#  Pmax * tanh(AQY * light_value / Pmax) - Rd
#}

##### MODEL PARAMETERS ######
# Ps = photosynthetic scaling parameter
# AQY = alpha = initial low-light slope
# beta = photoinhibition parameter
# Rd = dark respiration magnitude
# beta = photoinhibition parameter (beta ≈ 0 little/no photoinhibition, larger beta = stronger high-light decline)


# empty list to save individual Model 2 results
M2_individual_list <- list()

# make forloop






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




