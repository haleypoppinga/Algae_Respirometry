###### Respo Code for PI Curves ####### 
### Created by: Nyssa Silbiger, Maya Powell, Haley Poppinga
#### Last updated on: 2026-09-14

############## Introduction to code/script ####################
## this script will help us process the raw data gathered during respirometry runs. 

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
ch.vol <- 482 #mL #of small chambers x12

######### Load and tidy files ###############
############################################
#Load your respiration data file, with all the times, water volumes(mL), wet weight, dry weight, AFDW, zone
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
  mutate(date = ymd(date))

#view(Sample_Info)

#generate a 4 column dataframe with specific column names
# data is in umol.L.sec

n_light_levels<-10 # number of unique light levels (up to 900 PAR)

RespoR <- tibble(.rows =length(filenames_final)*n_light_levels, #empty dataframe
                 file_id_csv = NA_character_, # file ID used to join back to metadata
                 sample_id = NA,
                 intercept = NA,
                 umol.L.sec = NA,
                 temp_c = NA,
                 light_level = NA,
                 light_value = NA,
                 run_block = NA) # fills in as forloop goes


######### Create a for loop! ###############
############################################

###forloop#####
#to diagnose - you can set i = 1 and run through the loop line by line
#if you run through the loop and it fails:
#type "i" into the console and it will tell you what sample it failed on

for(i in 1:length(filenames_final)) {
  FRow <- as.numeric(which(Sample_Info$file_id_csv==filenames_final[i])) # stringsplit this renames our file
  
  Respo.Data1 <- read_csv(skip=1,file.path(path.p, paste0(file.names.full[i]))) %>% # reads in each file in list
    dplyr::select(Date, Time, Value, Temp) %>% # keep only what we need: Time stamp per 1sec, Raw O2 value per 1sec, in situ temp per 1sec
    unite(Date,Time,col="Time",remove=T, sep = " ") %>%
    drop_na() %>% 
    mutate(Time = mdy_hms(Time)) #%>% # convert time
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
      #   mutate(light_level = light_level,
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
    pdf(paste0(here("Output","PI"),"/",rename,"_",j,".pdf" ))
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
  # make sure the light value (or whatever other metadata you want) is in the final DF
  
  ################################
  # fill in all the O2 consumption and rate data
  
  RespoR[FRow,"file_id_csv"] <- filenames_final[i] # file ID for metadata join
  RespoR[FRow,"temp_c"]<-df$temp_c
  RespoR[FRow,"sample_id"]<-df$sample_id
  RespoR[FRow,"intercept"]<-df$intercept
  RespoR[FRow,"umol.L.sec"]<-df$umol.L.sec
  RespoR[FRow,"light_level"]<-df$light_level
  RespoR[FRow,"light_value"]<-df$light_value
  RespoR[FRow,"run_block"]<-df$run_block
  
}  


######### end of for loop - celebrate victory of getting through that ###############
############################################

###Get temperature data across all runs
median(Respo.Data1$Temp) #25.015
mean(Respo.Data1$Temp) #25.00991
max(Respo.Data1$Temp) #25.17
min(Respo.Data1$Temp) #24.749


#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
write_csv(RespoR, here("Data","Respo_Files","PI","Respo_Algae_R.csv"))  

##### 

RespoR <- read_csv(here("Data","Respo_Files","PI","Respo_Algae_R.csv"))


beep(sound = 8, expr = NULL)
