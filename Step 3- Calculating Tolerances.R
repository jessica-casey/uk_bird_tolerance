# TOLERANCES OF UK SPECIES

#### PREPARATIONS ####

# Necessary packages
library(sf)
library(dplyr)
library(scam)
library(purrr)
library(terra)

# Import the data: sampling events, continent-scale species observations, resolved taxonomy
load("C:/Data1/UK/uk_sampling_events_filtered.RData")
load("C:/Data1/UK/uk_sp_filtered.RData")
load("C:/Data1/UK/uk_resolved_taxonomy.RData")
load("C:/Data1/UK/tolerances.RData")
uk_full_name_list <- uk_full_name_list[order(uk_full_name_list$scientific_ebird), ]

# Remove oceanic species
oceanic_species <- read.csv("C:/Data1/UK/oceanic_bird_list_birdlife.csv")
oceanic_species <- oceanic_species$Scientific.name
uk_full_name_list <- uk_full_name_list[-c(which(uk_full_name_list$scientific_ebird %in% oceanic_species)), ]
rm(oceanic_species)

# Add a small constant to HFI values that are zero to allow including them into the categories
uk_sampling_events_filtered$hfi_mean <- ifelse(uk_sampling_events_filtered$hfi_mean == 0, uk_sampling_events_filtered$hfi_mean + 0.00001, uk_sampling_events_filtered$hfi_mean)

#Remove unnecessary columns from sampling event dataframe
uk_sampling_events_filtered <- uk_sampling_events_filtered[, -c(4:7, 12, 14:15, 17)]

# Remove unnecessary columns from species dataframe
uk_sp_filtered <- uk_sp_filtered[, -c(2, 5:14)]

sp_full <- unique(uk_sp_filtered$scientific_name)
uk_sp_filtered <- uk_sp_filtered[uk_sp_filtered$scientific_name %in% uk_full_name_list$scientific_ebird, ]


#### QUANTIFY ANTHROPOGENIC DISTURBANCE TOLERANCE FOR SINGLE SPECIES ####

new_data <- data.frame(hfi_mean = seq(0, 50, by = 0.1), protocol_type = rep("Traveling", 501), duration_minutes = rep(60, 501), effort_distance_km = rep(1, 501), number_observers = rep(1, 501), Hour = rep(7, 501))

n_bootstrap <- 1
n_tolerances <- 5
oversampling_tolerance <- 1.3

gamma_value <- 1.4

# Repeat everything for each species
for (bird_species in uk_full_name_list$scientific_ebird){
  
  setwd("C:/Data1/UK")
  # Create a matrix to store the tolerance values in
  all_tolerances <- as.data.frame(matrix(nrow = 1, ncol = n_bootstrap * n_tolerances))
  rownames(all_tolerances) <- bird_species
  colnames(all_tolerances) <- c(c(paste0("peak", 1:n_bootstrap)), c(paste0("high_hfi_0.5 ", 1:n_bootstrap)), c(paste0("high_hfi_0.1 ", 1:n_bootstrap)), c(paste0("low_hfi_0.5 ", 1:n_bootstrap)), c(paste0("low_hfi_0.1 ", 1:n_bootstrap)))
  
  # Select species-specific observations from the full dataset
  # Select only those checklists that are done within the species' breeding or resident range
  
  sp_data <- uk_sp_filtered[uk_sp_filtered$scientific_name == bird_species, ]
  sp_data <- as.data.frame(merge(sp_data, uk_sampling_events_filtered, by = "checklist_id", all.y = TRUE))
  sp_data$observation_count <- replace(sp_data$observation_count, is.na(sp_data$observation_count), values = 0)
  sp_data$effort_distance_km <- replace(sp_data$effort_distance_km, is.na(sp_data$effort_distance_km), values = 0)
  sp_data <- cbind(sf::st_as_sf(sp_data, coords = c("lon", "lat"), crs = 4326), sp_data$lat, sp_data$lon)
  sp_data <- st_transform(sp_data, crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  #sp_data$select <- lengths(st_intersects(sp_data, sp_range_union))
  #sp_data <- sp_data[sp_data$select > 0, ]
  #sp_data <- st_drop_geometry(sp_data)
  sp_data$select <- NULL
  #rm(sp_range_union, sp_range)
  
  if(nrow(sp_data) < 60 | sum(sp_data$observation_count) < 5) { next }
  
  else {
    
    for (j in 1:n_bootstrap){
      message('Processing species ', bird_species, ' bootstrap ', j)
      
      # Take 75% of the species-specific data for the next steps within the bootstrap
      # Sample observations depending on the HFI values to obtain a balanced sample of observations per species
      # Create slices for which the sampling probabilities are set
      
      sp_data_sample <- sp_data[sample(nrow(sp_data), size = round(nrow(sp_data)*0.75)), ]
      sp_data_sample$slice <- cut(sp_data_sample$hfi_mean, breaks = seq(0, 50, by = 5))
      sp_data_sample$slice <- as.numeric(sp_data_sample$slice)
      
      # Calculate the proportion of presences out of total number of checklists done within each slice of HFI values
      # If the maximum proportion of presences out of the total number of checklists within a slice is less than half, scale the proportions up
      # This accounts for class imbalance in sampling and doesn't skew the subsampling towards the more abundant absences across the full dataset
      # Determine the target number of observations to sample from the full species data
      
      average_proportion <- sp_data_sample %>% dplyr::select(observation_count, slice) %>% group_by(slice, .drop = FALSE) %>% summarise(tot_check = n(), prop_obs = mean(observation_count))
      average_proportion$n_pos <- average_proportion$tot_check * average_proportion$prop_obs
      average_proportion$scaled <- ifelse(rep(max(average_proportion$prop_obs) < 0.5, times = nrow(average_proportion)), average_proportion$prop_obs * (0.5 / max(average_proportion$prop_obs)), average_proportion$prop_obs)
      average_proportion$target_all <- ifelse(average_proportion$tot_check < 100, average_proportion$tot_check, 100)
      average_proportion$target_pos <- round(average_proportion$target_all * average_proportion$scaled, digits = 0)
      average_proportion$target_neg <- average_proportion$target_all - average_proportion$target_pos
      average_proportion$positive_sampling_ratio <- average_proportion$target_pos / (average_proportion$tot_check * average_proportion$prop_obs)
      average_proportion$negative_sampling_ratio <- average_proportion$target_neg / (average_proportion$tot_check * (1 - average_proportion$prop_obs))
      average_proportion[!is.finite(average_proportion$positive_sampling_ratio), "positive_sampling_ratio"] <- 0
      average_proportion[!is.finite(average_proportion$negative_sampling_ratio), "negative_sampling_ratio"] <- 0
      average_proportion$target_pos_adjusted <- ifelse(average_proportion$positive_sampling_ratio > oversampling_tolerance | average_proportion$negative_sampling_ratio > oversampling_tolerance, round(average_proportion$n_pos * oversampling_tolerance), average_proportion$target_pos)
      average_proportion$target_all_adjusted <- ifelse(average_proportion$positive_sampling_ratio > oversampling_tolerance | average_proportion$negative_sampling_ratio > oversampling_tolerance, round(average_proportion$target_pos_adjusted / average_proportion$scaled), average_proportion$target_all)
      average_proportion$target_neg_adjusted <- ifelse(average_proportion$positive_sampling_ratio > oversampling_tolerance | average_proportion$negative_sampling_ratio > oversampling_tolerance, average_proportion$target_all_adjusted - average_proportion$target_pos_adjusted, average_proportion$target_neg)
      average_proportion$target_pos_adjusted <- replace(average_proportion$target_pos_adjusted, average_proportion$target_pos_adjusted == 0, 0.00001)
      average_proportion$target_neg_adjusted <- replace(average_proportion$target_neg_adjusted, average_proportion$target_neg_adjusted == 0, 0.00001)
      
      average_proportion <- as.data.frame(average_proportion)
      sp_data_sample <- as.data.frame(sp_data_sample)
      
      # Don't calculate tolerance for species that has fewer than 20 checklists making presence or absence observations or for species that has fewer than two checklists made in more than five HFI slices
      if (sum(average_proportion$target_pos_adjusted) < 30 | sum(average_proportion$target_neg_adjusted) < 30 | sum(average_proportion$tot_check < 10) > 5) { next }
      
      else {
        
        # Create a sampling probability column in the sp_data dataframe based on the slice-specific probabilities. Divide the species data back into presences and absences.
        sp_data_sample <- merge(sp_data_sample, average_proportion, by = "slice")
        sp_data_sample$slice <- as.factor(sp_data_sample$slice)
        sp_data_presence <- sp_data_sample[sp_data_sample$observation_count == 1, ]
        sp_data_absence <- sp_data_sample[sp_data_sample$observation_count == 0, ]
        
        # Sample presences and absences according to the target numbers
        sp_datapres <- sp_data_presence %>% group_split(slice, .drop = FALSE) %>% map2_dfr(as.matrix(round(average_proportion[average_proportion$target_pos_adjusted != 0, "target_pos_adjusted"])), ~ slice_sample(.x, n = .y, replace = TRUE))
        sp_dataabs <- sp_data_absence %>% group_split(slice, .drop = FALSE) %>% map2_dfr(as.matrix(average_proportion[average_proportion$target_neg_adjusted != 0, "target_neg_adjusted"]), ~ slice_sample(.x, n = .y, replace = TRUE))
        sp_data_for_model <- as.data.frame(rbind(sp_datapres, sp_dataabs))
        
        if (which.max(table(as.numeric(sp_datapres$slice))) == 10){
          
          all_tolerances[bird_species, j] <- 50 
          all_tolerances[bird_species, j+n_bootstrap] <- 50
          all_tolerances[bird_species, j+n_bootstrap*2] <- 50 
          
        } else {
          # Model
          message('Starting model for ', bird_species, ' bootstrap ', j)
          skip_to_next <- FALSE
          sp_model <- tryCatch(scam(observation_count ~ s(hfi_mean, k = 9, bs = "cv") + protocol_type + duration_minutes + effort_distance_km + Year + number_observers + s(Hour, bs = "cc", k = 9), family = binomial, gamma = gamma_value, data = sp_data_for_model), error = function(e) { skip_to_next <<- TRUE})
          message('Finished model for ', bird_species, ' bootstrap ', j)
          
          if (skip_to_next) { next }
          
          # Predict occurrences across the range of HFI values
          new_data1 <- cbind(new_data, sp_data.lat = rep(mean(sp_data_for_model$sp_data.lat), 501), sp_data.lon = rep(mean(sp_data_for_model$sp_data.lon), 501), Year = rep(max(sp_data_for_model$Year), 501))
          new_data1$predicted_sp <- as.numeric(predict(sp_model, newdata = new_data1, type = "response"))
          
          plot(new_data1$hfi_mean, new_data1$predicted_sp)
          
          # Calculate the three tolerance measures if model has converged properly
          # Measure 1: Peak tolerance = level of HFI where predicted occurrence probability is at the maximum
          
          if (length(new_data1[new_data1$predicted_sp == max(new_data1$predicted_sp), "hfi_mean"]) > 1){ next }
          
          else {
            all_tolerances[bird_species, j] <- new_data1[new_data1$predicted_sp == max(new_data1$predicted_sp), "hfi_mean"]
            
            if (all_tolerances[bird_species, j] != 50 & !is.na(all_tolerances[bird_species, j])){
              # Measure 2: High HFI 0.5 tolerance = level of HFI where predicted occurrence probability is 50% of the maximum predicted occurrence probability
              # Measure 3: High HFI 0.1 tolerance = level of HFI where predicted occurrence probability is 10% of the maximum predicted occurrence probability
              # To address the higher HFI end of the curve, include only those rows where occurrence probabilities have been predicted for HFI values higher than Measure 1 HFI level
              new_data2 <- new_data1[new_data1$hfi_mean > all_tolerances[bird_species, j], ]
              all_tolerances[bird_species, j+n_bootstrap]  <- new_data2[which.min(abs(new_data2$predicted_sp-max(new_data2$predicted_sp)*0.5)), "hfi_mean"]
              all_tolerances[bird_species, j+n_bootstrap*2] <- new_data2[which.min(abs(new_data2$predicted_sp-max(new_data2$predicted_sp)*0.1)), "hfi_mean"]
              
            } else {
              all_tolerances[bird_species, j+n_bootstrap]  <- all_tolerances[bird_species, j]
              all_tolerances[bird_species, j+n_bootstrap*2] <- all_tolerances[bird_species, j]
            }
            
            if (all_tolerances[bird_species, j] != 0 & !is.na(all_tolerances[bird_species, j])){
              # Measure 4: Low HFI 0.5 tolerance = level of HFI where predicted occurrence probability is 50% of the maximum predicted occurrence probability, on the left side
              # Measure 5: Low HFI 0.1 tolerance = level of HFI where predicted occurrence probability is 10% of the maximum predicted occurrence probability, on the left side
              # To address the lower HFI end of the curve, include only those rows where occurrence probabilities have been predicted for HFI values lower than Measure 1 HFI level
              new_data3 <- new_data1[new_data1$hfi_mean < all_tolerances[bird_species, j], ]
              all_tolerances[bird_species, j+n_bootstrap*3] <- new_data3[which.min(abs(new_data3$predicted_sp-max(new_data3$predicted_sp)*0.5)), "hfi_mean"]
              all_tolerances[bird_species, j+n_bootstrap*4] <- new_data3[which.min(abs(new_data3$predicted_sp-max(new_data3$predicted_sp)*0.1)), "hfi_mean"]
              
            } else {
              all_tolerances[bird_species, j+n_bootstrap*3] <- all_tolerances[bird_species, j]
              all_tolerances[bird_species, j+n_bootstrap*4] <- all_tolerances[bird_species, j]
            }
          }
        }
      } 
    } # bootstrap closed
    #Full_tolerances <- rbind(Full_tolerances, all_tolerances)
    #save(Full_tolerances, file = ("C:/Data1/UK/tolerances.RData"))
    gc()
  }
  gc()          
} # species for loop closed

rm("average_proportion", "new_data1", "new_data2", "new_data3", "sp_data_absence", "sp_data_for_model", "sp_data_presence", "sp_data_sample", "sp_dataabs", "sp_datapres", "sp_model")

save(Full_tolerances, file = "C:/Data1/UK/tolerances_uk.RData")
