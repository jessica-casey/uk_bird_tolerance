# FILTERING AND COMPILING UK EBIRD DATA

#### PREPARATIONS ####

# Necessary packages
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(raster)
library(data.table)
library(ebirdst)
library(mgcv)
library(scam)
library(purrr)

# 'Continent' specifications

#need table with column titles: Country, ISO_code, Band1, Band2. Band columns are 0 or 1
uk <- read.csv("C:/Data1/UK/uk_bands.csv")
#Gets which countries (and their codes) are in each band
countries <- sort(uk$Country)
iso2 <- sort(uk$ISO_code)
band1 <- sort(uk[uk$Band1 == 1, "ISO_code"])
band2 <- sort(uk[uk$Band2 == 1, "ISO_code"])
band3 <- sort(uk[uk$Band3 == 1, "ISO_code"])
band4 <- sort(uk[uk$Band4 == 1, "ISO_code"])
band5 <- sort(uk[uk$Band5 == 1, "ISO_code"])

#Bounding boxes for each latitudinal band
#c(lng_min, lat_min, lng_max, lat_max)
band1_bbox <- c(-9, 55, 2, 60)
band2_bbox <- c(-9, 53.75, 2, 55)
band3_bbox <- c(-9, 52.5, 2, 53.75)
band4_bbox <- c(-9, 51.25, 2, 52.5)
band5_bbox <- c(-9, 50, 2, 51.25)

#Breeding season dates for each band
band1_dates <- c("*-04-01", "*-07-31")
band2_dates <- c("*-04-01", "*-07-31")
band3_dates <- c("*-04-01", "*-07-31")
band4_dates <- c("*-04-01", "*-07-31")
band5_dates <- c("*-04-01", "*-07-31")

# Set the working folder and import the data
auk::auk_set_ebd_path("C:/Data1/", overwrite = TRUE)

# Import Human Footprint Index data
hfi_full <- raster("C:/Data1/hfp2013_merisINT.tif")

#### Filter band 1 data ####

# Create a path for the filtered dataset
data_dir <- "C:/Data1/UK/Band1" 

for (country_code in band1){
  ebd_sampling <- auk_ebd(file = paste("C:/Data1/UK/ebd_", country_code ,"_smp_relApr-2024.txt", sep = ""),
                          file_sampling = "C:/Data1/UK/ebd_GB_smp_relApr-2024_sampling.txt")
  #Determine auk filters
  ebd_filters <- ebd_sampling %>%
    auk_country(country = country_code) %>%
    auk_bbox(bbox = band1_bbox) %>%
    auk_year(year = c(2018:2023)) %>%
    auk_date(date = band1_dates) %>%
    auk_protocol(protocol = c("Stationary", "Traveling")) %>%
    auk_duration(duration = c(30, 600)) %>%
    auk_distance(distance = c(0,8)) %>%
    auk_complete()
  
  f_sampling <- file.path(data_dir, paste("ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  f_ebd <- file.path(data_dir, paste("ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Filter the full dataset and save file in the path 
  # TIME CONSUMING STEP
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE) 
  
  sampling_events_filtered <- read_sampling(paste("C:/Data1/UK/Band1/ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  sp_filtered <- read_ebd(paste("C:/Data1/UK/Band1/ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Make observation_date into separate year, month, week and day columns
  sampling_events_filtered$Week <- as.numeric(data.table::week(sampling_events_filtered$observation_date))
  sampling_events_filtered$observation_date_to_break <- sampling_events_filtered$observation_date
  sampling_events_filtered <- sampling_events_filtered %>% separate(observation_date_to_break, c("Year","Month", "Day"), sep = "-")
  sampling_events_filtered$Year <- as.numeric(sampling_events_filtered$Year)
  sampling_events_filtered$Month <- as.numeric(sampling_events_filtered$Month)
  sampling_events_filtered$Day <- as.numeric(sampling_events_filtered$Day)
  sampling_events_filtered$Hour <- hour(as.POSIXct(sampling_events_filtered$time_observations_started, format = "%H:%M:%S"))
  
  # Change names of some columns for the next step
  colnames(sampling_events_filtered)[c(16:17)] <- c("lat", "lon")
  sampling_events_filtered$date <- yday(sampling_events_filtered$observation_date)
  
  # Filter the data spatially so that only one observation per week per 3x3 kilometer grid cell across years is included
  sampling_events_filtered <- grid_sample(sampling_events_filtered, coords = c("lon", "lat", "date"), res = c(3000, 3000, 7), jitter_grid = FALSE)
  sp_filtered <- sp_filtered[sp_filtered$checklist_id %in% sampling_events_filtered$checklist_id, ]
  
  # Include only approved observations
  sp_filtered <- sp_filtered[sp_filtered$approved == TRUE, ]
  
  # Make data into presence-absence data
  sp_filtered$observation_count <- 1
  
  # Select only those checklists from the sampling events data that are included in the final species data
  sampling_events_filtered <- sampling_events_filtered[sampling_events_filtered$checklist_id %in% sp_filtered$checklist_id, ]
  
  # Remove unnecessary columns from species data and save it
  sp_filtered <- sp_filtered[, c(1, 7:8, 10, 28:34, 37:38, 40)]
  save(sp_filtered, file = paste(data_dir, "/", country_code, "_sp_filtered.RData", sep = ""))
  rm(sp_filtered)
  
  # Make the sampling event data spatial
  sampling_events_filtered_spatial <- sf::st_as_sf(sampling_events_filtered, coords = c("lon", "lat"), crs = 4326)
  sampling_events_filtered_spatial <- st_transform(sampling_events_filtered_spatial, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  # Build a 3 kilometer buffer (1,5 km radius) around each checklist point and calculate the mean of the human footprint index inside the buffer
  # TIME CONSUMING STEP
  sampling_events_filtered$hfi_mean <- raster::extract(x = hfi_full, y = sampling_events_filtered_spatial, buffer = 1500, fun = "mean")
  
  # Remove rows with missing HFI values
  sampling_events_filtered <- sampling_events_filtered[complete.cases(sampling_events_filtered[ , "hfi_mean"]),]
  
  # Remove unnecessary columns from sampling event data and save it
  sampling_events_filtered <- sampling_events_filtered[, c(1, 16:22, 25:26, 28, 32:37, 38)]
  save(sampling_events_filtered, file = paste(data_dir, "/", country_code, "_sampling_events_filtered.RData", sep = ""))
}

# Combine datasets within the latitudinal band
load("C:/Data1/Uk/Band1/GB_sampling_events_filtered.RData")
sampling_events_filtered_Band1 <- sampling_events_filtered

load("C:/Data1/Uk/Band1/GB_sp_filtered.RData")
sp_filtered_Band1 <- sp_filtered

save(sampling_events_filtered_Band1, file = "C:/Data1/UK/Band1/sampling_events_filtered_Band1.RData")
save(sp_filtered_Band1, file = "C:/Data1/UK/Band1/sp_filtered_Band1.RData")

#### Filter band 2 data ####

# Create a path for the filtered dataset
data_dir <- "C:/Data1/UK/Band2" 

for (country_code in band2){
  ebd_sampling <- auk_ebd(file = paste("C:/Data1/UK/ebd_", country_code ,"_smp_relApr-2024.txt", sep = ""),
                          file_sampling = "C:/Data1/UK/ebd_GB_smp_relApr-2024_sampling.txt")
  #Determine auk filters
  ebd_filters <- ebd_sampling %>%
    auk_country(country = country_code) %>%
    auk_bbox(bbox = band2_bbox) %>%
    auk_year(year = c(2018:2023)) %>%
    auk_date(date = band2_dates) %>%
    auk_protocol(protocol = c("Stationary", "Traveling")) %>%
    auk_duration(duration = c(30, 600)) %>%
    auk_distance(distance = c(0,8)) %>%
    auk_complete()
  
  f_sampling <- file.path(data_dir, paste("ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  f_ebd <- file.path(data_dir, paste("ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Filter the full dataset and save file in the path 
  # TIME CONSUMING STEP
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE) 
  
  sampling_events_filtered <- read_sampling(paste("C:/Data1/UK/Band2/ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  sp_filtered <- read_ebd(paste("C:/Data1/UK/Band2/ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Make observation_date into separate year, month, week and day columns
  sampling_events_filtered$Week <- as.numeric(data.table::week(sampling_events_filtered$observation_date))
  sampling_events_filtered$observation_date_to_break <- sampling_events_filtered$observation_date
  sampling_events_filtered <- sampling_events_filtered %>% separate(observation_date_to_break, c("Year","Month", "Day"), sep = "-")
  sampling_events_filtered$Year <- as.numeric(sampling_events_filtered$Year)
  sampling_events_filtered$Month <- as.numeric(sampling_events_filtered$Month)
  sampling_events_filtered$Day <- as.numeric(sampling_events_filtered$Day)
  sampling_events_filtered$Hour <- hour(as.POSIXct(sampling_events_filtered$time_observations_started, format = "%H:%M:%S"))
  
  # Change names of some columns for the next step
  colnames(sampling_events_filtered)[c(16:17)] <- c("lat", "lon")
  sampling_events_filtered$date <- yday(sampling_events_filtered$observation_date)
  
  # Filter the data spatially so that only one observation per week per 3x3 kilometer grid cell across years is included
  sampling_events_filtered <- grid_sample(sampling_events_filtered, coords = c("lon", "lat", "date"), res = c(3000, 3000, 7), jitter_grid = FALSE)
  sp_filtered <- sp_filtered[sp_filtered$checklist_id %in% sampling_events_filtered$checklist_id, ]
  
  # Include only approved observations
  sp_filtered <- sp_filtered[sp_filtered$approved == TRUE, ]
  
  # Make data into presence-absence data
  sp_filtered$observation_count <- 1
  
  # Select only those checklists from the sampling events data that are included in the final species data
  sampling_events_filtered <- sampling_events_filtered[sampling_events_filtered$checklist_id %in% sp_filtered$checklist_id, ]
  
  # Remove unnecessary columns from species data and save it
  sp_filtered <- sp_filtered[, c(1, 7:8, 10, 28:34, 37:38, 40)]
  save(sp_filtered, file = paste(data_dir, "/", country_code, "_sp_filtered.RData", sep = ""))
  rm(sp_filtered)
  
  # Make the sampling event data spatial
  sampling_events_filtered_spatial <- sf::st_as_sf(sampling_events_filtered, coords = c("lon", "lat"), crs = 4326)
  sampling_events_filtered_spatial <- st_transform(sampling_events_filtered_spatial, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  # Build a 3 kilometer buffer (1,5 km radius) around each checklist point and calculate the mean of the human footprint index inside the buffer
  # TIME CONSUMING STEP
  sampling_events_filtered$hfi_mean <- raster::extract(x = hfi_full, y = sampling_events_filtered_spatial, buffer = 1500, fun = "mean")
  
  # Remove rows with missing HFI values
  sampling_events_filtered <- sampling_events_filtered[complete.cases(sampling_events_filtered[ , "hfi_mean"]),]
  
  # Remove unnecessary columns from sampling event data and save it
  sampling_events_filtered <- sampling_events_filtered[, c(1, 16:22, 25:26, 28, 32:37, 38)]
  save(sampling_events_filtered, file = paste(data_dir, "/", country_code, "_sampling_events_filtered.RData", sep = ""))
}

# Combine datasets within the latitudinal band
load("C:/Data1/Uk/Band2/GB_sampling_events_filtered.RData")
sampling_events_filtered_band2 <- sampling_events_filtered

load("C:/Data1/Uk/Band2/GB_sp_filtered.RData")
sp_filtered_band2 <- sp_filtered

save(sampling_events_filtered_band2, file = "C:/Data1/UK/Band2/sampling_events_filtered_Band2.RData")
save(sp_filtered_band2, file = "C:/Data1/UK/Band2/sp_filtered_Band2.RData")


#### Filter band 3 data ####

# Create a path for the filtered dataset
data_dir <- "C:/Data1/UK/Band3" 

for (country_code in band3){
  ebd_sampling <- auk_ebd(file = paste("C:/Data1/UK/ebd_", country_code ,"_smp_relApr-2024.txt", sep = ""),
                          file_sampling = "C:/Data1/UK/ebd_GB_smp_relApr-2024_sampling.txt")
  #Determine auk filters
  ebd_filters <- ebd_sampling %>%
    auk_country(country = country_code) %>%
    auk_bbox(bbox = band3_bbox) %>%
    auk_year(year = c(2018:2023)) %>%
    auk_date(date = band3_dates) %>%
    auk_protocol(protocol = c("Stationary", "Traveling")) %>%
    auk_duration(duration = c(30, 600)) %>%
    auk_distance(distance = c(0,8)) %>%
    auk_complete()
  
  f_sampling <- file.path(data_dir, paste("ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  f_ebd <- file.path(data_dir, paste("ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Filter the full dataset and save file in the path 
  # TIME CONSUMING STEP
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE) 
  
  sampling_events_filtered <- read_sampling(paste("C:/Data1/UK/Band3/ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  sp_filtered <- read_ebd(paste("C:/Data1/UK/Band3/ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Make observation_date into separate year, month, week and day columns
  sampling_events_filtered$Week <- as.numeric(data.table::week(sampling_events_filtered$observation_date))
  sampling_events_filtered$observation_date_to_break <- sampling_events_filtered$observation_date
  sampling_events_filtered <- sampling_events_filtered %>% separate(observation_date_to_break, c("Year","Month", "Day"), sep = "-")
  sampling_events_filtered$Year <- as.numeric(sampling_events_filtered$Year)
  sampling_events_filtered$Month <- as.numeric(sampling_events_filtered$Month)
  sampling_events_filtered$Day <- as.numeric(sampling_events_filtered$Day)
  sampling_events_filtered$Hour <- hour(as.POSIXct(sampling_events_filtered$time_observations_started, format = "%H:%M:%S"))
  
  # Change names of some columns for the next step
  colnames(sampling_events_filtered)[c(16:17)] <- c("lat", "lon")
  sampling_events_filtered$date <- yday(sampling_events_filtered$observation_date)
  
  # Filter the data spatially so that only one observation per week per 3x3 kilometer grid cell across years is included
  sampling_events_filtered <- grid_sample(sampling_events_filtered, coords = c("lon", "lat", "date"), res = c(3000, 3000, 7), jitter_grid = FALSE)
  sp_filtered <- sp_filtered[sp_filtered$checklist_id %in% sampling_events_filtered$checklist_id, ]
  
  # Include only approved observations
  sp_filtered <- sp_filtered[sp_filtered$approved == TRUE, ]
  
  # Make data into presence-absence data
  sp_filtered$observation_count <- 1
  
  # Select only those checklists from the sampling events data that are included in the final species data
  sampling_events_filtered <- sampling_events_filtered[sampling_events_filtered$checklist_id %in% sp_filtered$checklist_id, ]
  
  # Remove unnecessary columns from species data and save it
  sp_filtered <- sp_filtered[, c(1, 7:8, 10, 28:34, 37:38, 40)]
  save(sp_filtered, file = paste(data_dir, "/", country_code, "_sp_filtered.RData", sep = ""))
  rm(sp_filtered)
  
  # Make the sampling event data spatial
  sampling_events_filtered_spatial <- sf::st_as_sf(sampling_events_filtered, coords = c("lon", "lat"), crs = 4326)
  sampling_events_filtered_spatial <- st_transform(sampling_events_filtered_spatial, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  # Build a 3 kilometer buffer (1,5 km radius) around each checklist point and calculate the mean of the human footprint index inside the buffer
  # TIME CONSUMING STEP
  sampling_events_filtered$hfi_mean <- raster::extract(x = hfi_full, y = sampling_events_filtered_spatial, buffer = 1500, fun = "mean")
  
  # Remove rows with missing HFI values
  sampling_events_filtered <- sampling_events_filtered[complete.cases(sampling_events_filtered[ , "hfi_mean"]),]
  
  # Remove unnecessary columns from sampling event data and save it
  sampling_events_filtered <- sampling_events_filtered[, c(1, 16:22, 25:26, 28, 32:37, 38)]
  save(sampling_events_filtered, file = paste(data_dir, "/", country_code, "_sampling_events_filtered.RData", sep = ""))
}

# Combine datasets within the latitudinal band
load("C:/Data1/Uk/Band3/GB_sampling_events_filtered.RData")
sampling_events_filtered_band3 <- sampling_events_filtered

load("C:/Data1/Uk/Band3/GB_sp_filtered.RData")
sp_filtered_band3 <- sp_filtered

save(sampling_events_filtered_band3, file = "C:/Data1/UK/Band3/sampling_events_filtered_Band3.RData")
save(sp_filtered_band3, file = "C:/Data1/UK/Band3/sp_filtered_Band3.RData")


#### Filter band 4 data ####

# Create a path for the filtered dataset
data_dir <- "C:/Data1/UK/Band4" 

for (country_code in band4){
  ebd_sampling <- auk_ebd(file = paste("C:/Data1/UK/ebd_", country_code ,"_smp_relApr-2024.txt", sep = ""),
                          file_sampling = "C:/Data1/UK/ebd_GB_smp_relApr-2024_sampling.txt")
  #Determine auk filters
  ebd_filters <- ebd_sampling %>%
    auk_country(country = country_code) %>%
    auk_bbox(bbox = band4_bbox) %>%
    auk_year(year = c(2018:2023)) %>%
    auk_date(date = band4_dates) %>%
    auk_protocol(protocol = c("Stationary", "Traveling")) %>%
    auk_duration(duration = c(30, 600)) %>%
    auk_distance(distance = c(0,8)) %>%
    auk_complete()
  
  f_sampling <- file.path(data_dir, paste("ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  f_ebd <- file.path(data_dir, paste("ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Filter the full dataset and save file in the path 
  # TIME CONSUMING STEP
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE) 
  
  sampling_events_filtered <- read_sampling(paste("C:/Data1/UK/Band4/ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  sp_filtered <- read_ebd(paste("C:/Data1/UK/Band4/ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Make observation_date into separate year, month, week and day columns
  sampling_events_filtered$Week <- as.numeric(data.table::week(sampling_events_filtered$observation_date))
  sampling_events_filtered$observation_date_to_break <- sampling_events_filtered$observation_date
  sampling_events_filtered <- sampling_events_filtered %>% separate(observation_date_to_break, c("Year","Month", "Day"), sep = "-")
  sampling_events_filtered$Year <- as.numeric(sampling_events_filtered$Year)
  sampling_events_filtered$Month <- as.numeric(sampling_events_filtered$Month)
  sampling_events_filtered$Day <- as.numeric(sampling_events_filtered$Day)
  sampling_events_filtered$Hour <- hour(as.POSIXct(sampling_events_filtered$time_observations_started, format = "%H:%M:%S"))
  
  # Change names of some columns for the next step
  colnames(sampling_events_filtered)[c(16:17)] <- c("lat", "lon")
  sampling_events_filtered$date <- yday(sampling_events_filtered$observation_date)
  
  # Filter the data spatially so that only one observation per week per 3x3 kilometer grid cell across years is included
  sampling_events_filtered <- grid_sample(sampling_events_filtered, coords = c("lon", "lat", "date"), res = c(3000, 3000, 7), jitter_grid = FALSE)
  sp_filtered <- sp_filtered[sp_filtered$checklist_id %in% sampling_events_filtered$checklist_id, ]
  
  # Include only approved observations
  sp_filtered <- sp_filtered[sp_filtered$approved == TRUE, ]
  
  # Make data into presence-absence data
  sp_filtered$observation_count <- 1
  
  # Select only those checklists from the sampling events data that are included in the final species data
  sampling_events_filtered <- sampling_events_filtered[sampling_events_filtered$checklist_id %in% sp_filtered$checklist_id, ]
  
  # Remove unnecessary columns from species data and save it
  sp_filtered <- sp_filtered[, c(1, 7:8, 10, 28:34, 37:38, 40)]
  save(sp_filtered, file = paste(data_dir, "/", country_code, "_sp_filtered.RData", sep = ""))
  rm(sp_filtered)
  
  # Make the sampling event data spatial
  sampling_events_filtered_spatial <- sf::st_as_sf(sampling_events_filtered, coords = c("lon", "lat"), crs = 4326)
  sampling_events_filtered_spatial <- st_transform(sampling_events_filtered_spatial, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  # Build a 3 kilometer buffer (1,5 km radius) around each checklist point and calculate the mean of the human footprint index inside the buffer
  # TIME CONSUMING STEP
  sampling_events_filtered$hfi_mean <- raster::extract(x = hfi_full, y = sampling_events_filtered_spatial, buffer = 1500, fun = "mean")
  
  # Remove rows with missing HFI values
  sampling_events_filtered <- sampling_events_filtered[complete.cases(sampling_events_filtered[ , "hfi_mean"]),]
  
  # Remove unnecessary columns from sampling event data and save it
  sampling_events_filtered <- sampling_events_filtered[, c(1, 16:22, 25:26, 28, 32:37, 38)]
  save(sampling_events_filtered, file = paste(data_dir, "/", country_code, "_sampling_events_filtered.RData", sep = ""))
}

# Combine datasets within the latitudinal band
load("C:/Data1/Uk/Band4/GB_sampling_events_filtered.RData")
sampling_events_filtered_band4 <- sampling_events_filtered

load("C:/Data1/Uk/Band4/GB_sp_filtered.RData")
sp_filtered_band4 <- sp_filtered

save(sampling_events_filtered_band4, file = "C:/Data1/UK/Band4/sampling_events_filtered_Band4.RData")
save(sp_filtered_band4, file = "C:/Data1/UK/Band4/sp_filtered_Band4.RData")

#### Filter band 5 data ####

# Create a path for the filtered dataset
data_dir <- "C:/Data1/UK/Band5" 

for (country_code in band4){
  ebd_sampling <- auk_ebd(file = paste("C:/Data1/UK/ebd_", country_code ,"_smp_relApr-2024.txt", sep = ""),
                          file_sampling = "C:/Data1/UK/ebd_GB_smp_relApr-2024_sampling.txt")
  #Determine auk filters
  ebd_filters <- ebd_sampling %>%
    auk_country(country = country_code) %>%
    auk_bbox(bbox = band5_bbox) %>%
    auk_year(year = c(2018:2023)) %>%
    auk_date(date = band5_dates) %>%
    auk_protocol(protocol = c("Stationary", "Traveling")) %>%
    auk_duration(duration = c(30, 600)) %>%
    auk_distance(distance = c(0,8)) %>%
    auk_complete()
  
  f_sampling <- file.path(data_dir, paste("ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  f_ebd <- file.path(data_dir, paste("ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Filter the full dataset and save file in the path 
  # TIME CONSUMING STEP
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE) 
  
  sampling_events_filtered <- read_sampling(paste("C:/Data1/UK/Band5/ebd_checklists_", country_code, "_filtered.txt", sep = ""))
  sp_filtered <- read_ebd(paste("C:/Data1/UK/Band5/ebd_", country_code, "_filtered.txt", sep = ""))
  
  # Make observation_date into separate year, month, week and day columns
  sampling_events_filtered$Week <- as.numeric(data.table::week(sampling_events_filtered$observation_date))
  sampling_events_filtered$observation_date_to_break <- sampling_events_filtered$observation_date
  sampling_events_filtered <- sampling_events_filtered %>% separate(observation_date_to_break, c("Year","Month", "Day"), sep = "-")
  sampling_events_filtered$Year <- as.numeric(sampling_events_filtered$Year)
  sampling_events_filtered$Month <- as.numeric(sampling_events_filtered$Month)
  sampling_events_filtered$Day <- as.numeric(sampling_events_filtered$Day)
  sampling_events_filtered$Hour <- hour(as.POSIXct(sampling_events_filtered$time_observations_started, format = "%H:%M:%S"))
  
  # Change names of some columns for the next step
  colnames(sampling_events_filtered)[c(16:17)] <- c("lat", "lon")
  sampling_events_filtered$date <- yday(sampling_events_filtered$observation_date)
  
  # Filter the data spatially so that only one observation per week per 3x3 kilometer grid cell across years is included
  sampling_events_filtered <- grid_sample(sampling_events_filtered, coords = c("lon", "lat", "date"), res = c(3000, 3000, 7), jitter_grid = FALSE)
  sp_filtered <- sp_filtered[sp_filtered$checklist_id %in% sampling_events_filtered$checklist_id, ]
  
  # Include only approved observations
  sp_filtered <- sp_filtered[sp_filtered$approved == TRUE, ]
  
  # Make data into presence-absence data
  sp_filtered$observation_count <- 1
  
  # Select only those checklists from the sampling events data that are included in the final species data
  sampling_events_filtered <- sampling_events_filtered[sampling_events_filtered$checklist_id %in% sp_filtered$checklist_id, ]
  
  # Remove unnecessary columns from species data and save it
  sp_filtered <- sp_filtered[, c(1, 7:8, 10, 28:34, 37:38, 40)]
  save(sp_filtered, file = paste(data_dir, "/", country_code, "_sp_filtered.RData", sep = ""))
  rm(sp_filtered)
  
  # Make the sampling event data spatial
  sampling_events_filtered_spatial <- sf::st_as_sf(sampling_events_filtered, coords = c("lon", "lat"), crs = 4326)
  sampling_events_filtered_spatial <- st_transform(sampling_events_filtered_spatial, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  # Build a 3 kilometer buffer (1,5 km radius) around each checklist point and calculate the mean of the human footprint index inside the buffer
  # TIME CONSUMING STEP
  sampling_events_filtered$hfi_mean <- raster::extract(x = hfi_full, y = sampling_events_filtered_spatial, buffer = 1500, fun = "mean")
  
  # Remove rows with missing HFI values
  sampling_events_filtered <- sampling_events_filtered[complete.cases(sampling_events_filtered[ , "hfi_mean"]),]
  
  # Remove unnecessary columns from sampling event data and save it
  sampling_events_filtered <- sampling_events_filtered[, c(1, 16:22, 25:26, 28, 32:37, 38)]
  save(sampling_events_filtered, file = paste(data_dir, "/", country_code, "_sampling_events_filtered.RData", sep = ""))
}

# Combine datasets within the latitudinal band
load("C:/Data1/Uk/Band5/GB_sampling_events_filtered.RData")
sampling_events_filtered_band5 <- sampling_events_filtered

load("C:/Data1/Uk/Band5/GB_sp_filtered.RData")
sp_filtered_band5 <- sp_filtered

save(sampling_events_filtered_band5, file = "C:/Data1/UK/Band5/sampling_events_filtered_Band5.RData")
save(sp_filtered_band5, file = "C:/Data1/UK/Band5/sp_filtered_Band5.RData")
