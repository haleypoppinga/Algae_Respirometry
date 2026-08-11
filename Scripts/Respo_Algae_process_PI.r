###### Respo Code for Algae PI Curve ####### 
### Created by: Haley Poppinga, Maya Powell, Nyssa Silbiger
#### Last updated on: 2026-08-10
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

##### PLOTTING CURVES #####
##### Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)

#Plot curves
av_resp <- RespoR_Normalized %>% filter(species == "av")
as_resp <- RespoR_Normalized %>% filter(species == "as")
gs_resp <- RespoR_Normalized %>% filter(species == "gs")
cs_resp <- RespoR_Normalized %>% filter(species == "cs")
sf_resp <- RespoR_Normalized %>% filter(species == "sf")
da_resp <- RespoR_Normalized %>% filter(species == "da")
ds_resp <- RespoR_Normalized %>% filter(species == "ds")
dc_resp <- RespoR_Normalized %>% filter(species == "dc")
hd_resp <- RespoR_Normalized %>% filter(species == "hd")

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

ds.mean <- aggregate(umol.g.hr ~ light_value, data = ds_resp, FUN=mean)
ds.se <- aggregate(umol.g.hr ~ light_value, data = ds_resp, FUN=std.error)

dc.mean <- aggregate(umol.g.hr ~ light_value, data = dc_resp, FUN=mean)
dc.se <- aggregate(umol.g.hr ~ light_value, data = dc_resp, FUN=std.error)

hd.mean <- aggregate(umol.g.hr ~ light_value, data = hd_resp, FUN=mean)
hd.se <- aggregate(umol.g.hr ~ light_value, data = hd_resp, FUN=std.error)


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
