###### Respo Code for Algae Light and Dark Runs ####### 
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

############# now it's time to code ############
################################################
# get the file path

#set the path to all of the raw oxygen datasheets
## these are saved onto the computer in whatever file path/naming scheme you saved things to 
path.p<-here("Data","Respo_Files","PI") #the location of all your respirometry files
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
#View(Sample.Info)


# join the data together
Sample_Info <- left_join(RespoMeta, BioData)
#View(Sample.Info)

##### Make sure times are consistent ####
# make start and stop times real times, so that we can join the respo output and sample_info data frames
Sample_Info <- Sample_Info %>% 
  #drop_na(sample_ID) %>% 
  unite(date,start_time,col="start_time",remove=F, sep=" ") %>% 
  unite(date,stop_time,col="stop_time",remove=F, sep=" ") %>%
  mutate(start_time = mdy_hms(start_time)) %>% 
  mutate(stop_time = mdy_hms(stop_time)) %>% 
  mutate(date = mdy(date))

#view(Sample_Info)

#generate a 4 column dataframe with specific column names
# data is in umol.L.sec

n_light_levels<-8 # number of unique light levels

RespoR <- tibble(.rows =length(filenames_final)*n_light_levels,
                 sample_id = NA,
                 intercept = NA,
                 umol.L.sec = NA,
                 temp_c = NA,
                 light_level = NA,
                 light_value = NA,
                 run_block = NA)

######### Create a for loop! ###############
############################################

###forloop#####
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
           Temp.C = map_dbl(map(data, "Temp"), mean),# get the mean temperature
           RegStats =map(regs, function(x){ # extract the intercept and slope for the parameters
             x$allRegs %>%
               slice(1) %>%
               select(Intercept = b0,
                      umol.L.sec = b1)
           }) )
  
  
  #  Plot regression diagnostics
  
  for(j in 1:length(df$light_level)){
    pdf(paste0(here("Output","PI"),"/",rename,"_",j,".pdf" ))
    plot(df$regs[[j]])
    dev.off() 
  }
  
  
  df<-df %>%
    select(light_level, temp_c, RegStats ) %>%
    unnest(RegStats) %>%
    mutate(sample_id = rename) %>%
    left_join(Sample_Info[FRow,] %>%
                select(sample_id, light_level, light_value, run_block)) # make sure the light value (or whatever other metadata you want) is in the final DF
  
  ################################
  # fill in all the O2 consumption and rate data
  
  RespoR[FRow,"temp_c"]<-df$Temp.C
  RespoR[FRow,"sample_id"]<-df$sample_ID
  RespoR[FRow,"intercept"]<-df$Intercept
  RespoR[FRow,"umol.L.sec"]<-df$umol.L.sec
  RespoR[FRow,"light_level"]<-df$light_level
  RespoR[FRow,"light_value"]<-df$light_value
  RespoR[FRow,"run_block"]<-df$run_block
  
}  


######### end of for loop - celebrate victory of getting through that ###############
############################################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
write_csv(RespoR, here("Data","RespoFiles","PI","Respo_Algae_R.csv"))  

##### 

RespoR <- read_csv(here("Data","RespoFiles","PI","Respo_Algae_R.csv"))

######### Calculate Respiration rate ###############
############################################

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
  #dplyr::select(Light_level, run_block, blank.rate, date) %>% # this is what I will use to join the blanks back with the raw data
  right_join(RespoR2) %>% # join blanks with the respo data
  mutate(umol.sec.corr = umol.sec - blank.rate, # subtract the blank rates from the raw rates   
         umol.cm2.hr = (umol.sec.corr*3600)/wetweight_g,
         umol.cm2.hr_uncorr = (umol.sec*3600)/wetweight_g) %>%
  # mmol.cm2.hr = 0.001*(umol.sec.corr*3600)/SA_cm2, # convert to mmol cm-2 hr-1
  # mmol.cm2.hr_uncorr = 0.001*(umol.sec*3600)/SA_cm2) %>% 
  filter(blank!=1) %>% # remove the Blank data
  #ungroup() %>%
  dplyr::select(date, species, sample_id, algae_id, light_dark, run_block, wetweight_g, run_block, umol.cm2.hr, chamber_channel, 
                temp_c, umol.cm2.hr_uncorr, light_level, light_value, light_dark) #keep only what we need


#######################
### making a df for just blank data for future use in plots ### 
Blank_only <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(light_level, light_value, run_block,light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

#############################
write_csv(RespoR_Normalized , here("Data","RespoFiles","PI","Respo_Algae_RNormalized_AllPIRates.csv"))  


## Plot the blanks across treatments to make sure nothing is funky
Blank_only %>%
  ggplot(aes(x = light_value, blank.rate)) +
  geom_point()

#  basic plot of rates versus light before you make the real PI curve 
basic_PI_plot <- RespoR_Normalized %>%
  ggplot(aes(x = light_value, y = umol.cm2.hr, color = species, group = algae_id))+
  geom_point()+
  geom_line()+
  facet_wrap(~species, scales = "free")
ggsave(here("Output","PI","basic_PI_plot.pdf"), basic_PI_plot)

### run an nls model for PI curve and extract Ik for each species ###