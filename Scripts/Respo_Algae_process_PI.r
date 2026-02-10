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
Sample_Info <- left_join(RespoMeta, BioData)
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
  
  
  # attach meta for this file (by light_level) - minimal and reliable
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
  select(file_id_csv, light_level, intercept, umol.L.sec, temp_c) %>%
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

RespoR2 <- RespoR %>%
  #drop_na(FileID_csv) %>% # drop NAs
  left_join(Sample_Info) %>% # Join the raw respo calculations with the metadata
  mutate(Ch.Volume.mL = volume_mL) %>% # 
  mutate(Ch.Volume.L = Ch.Volume.mL * 0.001) %>% # mL to L conversion
  mutate(umol.sec = umol.L.sec*Ch.Volume.L) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
  mutate_if(sapply(., is.character), as.factor)  #convert character columns to factors

#Account for blank rate by sample run Block (if we do at least one blank per block)

#View(RespoR)

####### normalize the respo rates to the blanks ##### 

RespoR_Normalized <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, run_block,light_dark) %>%
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
  dplyr::select(date, species, sample_id, algae_id, light_dark, run_block, wetweight_g, chamber_channel,
                temp_c, light_level, light_value,umol.sec, blank.rate, umol.sec.corr, 
                umol.g.hr, umol.g.hr_uncorr) #keep only what we need


#######################
### making a df for just blank data for future use in plots ### 
#Blank_only <- blank_rates %>% 
#  left_join(Sample_Info %>% distinct(light_level, run_block, light_dark, light_value),
#            by = c("Light_level" = "light_level","run_block","light_dark"))
### making a df for just blank data for future use in plots ### 

Blank_only <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, light_value, run_block,light_dark) %>%
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
basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value, y = umol.g.hr, color = species, group = algae_id)) +
  geom_point()+
  geom_line()+
  facet_wrap(~species, scales = "free_y")
ggsave(here("Output","PI","basic_PI_plot.pdf"), basic_PI_plot)


### run an nls model for PI curve and extract Ik for each species ###
beep(sound = 8, expr = NULL)
