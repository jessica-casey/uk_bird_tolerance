#### UK tolerance value compilation

#### PREPARATIONS ####

library(ggplot2)
library(paletteer)
library(ggpubr)
library(raster)
library(sf)
library(readxl)
library(ggridges)

setwd("C:/")

# Load the full file
load("C:/Data1/UK/tolerances.RData")
load("C:/Data1/UK/tolerances_with_metrics_uk.RData")

#### CALCULATE FINAL TOLERANCE METRICS FOR ALL SPECIES FROM BOOTSTRAPS ####

Full_tolerances$n_missing <- apply(is.na(Full_tolerances), MARGIN = 1, FUN = sum)
Full_tolerances$n_missing <- Full_tolerances$n_missing / 5 # divided by 5 because there are 5 tolerances

Full_tolerances <- Full_tolerances[Full_tolerances$n_missing != 50, ]
#hist(Full_tolerances$n_missing)


# Calculating breadth for each bootstrap
species_number <- nrow(Full_tolerances)
species_list <- rownames(Full_tolerances)
breadths <- as.data.frame(matrix(nrow = species_number, ncol =  100))
rownames(breadths) <- species_list
colnames(breadths) <- c(c(paste0("breadth_0.5 ", 1:50)), c(paste0("breadth_0.1 ", 1:50)))

for (i in 1:50){
  breadths[,i] <- Full_tolerances[, 50 + i] - Full_tolerances[, 150 + i]
}
for (i in 1:50){
  breadths[, 50 + i] <- Full_tolerances[, 100 + i] - Full_tolerances[, 200 + i]
}

Full_tolerances <- cbind(Full_tolerances, breadths)

# Calculating mean and sd for each tolerance category
Full_tolerances$mean_peak <- rowMeans(Full_tolerances[, 1:50], na.rm = TRUE)
Full_tolerances$sd_peak <- apply(Full_tolerances[, 1:50], 1, sd, na.rm = TRUE)  

Full_tolerances$mean_high_hfi_0.5 <- rowMeans(Full_tolerances[, 51:100], na.rm = TRUE)
Full_tolerances$sd_high_hfi_0.5 <- apply(Full_tolerances[, 51:100], 1, sd, na.rm = TRUE)  

Full_tolerances$mean_high_hfi_0.1 <- rowMeans(Full_tolerances[, 101:150], na.rm = TRUE)
Full_tolerances$sd_high_hfi_0.1 <- apply(Full_tolerances[, 101:150], 1, sd, na.rm = TRUE) 

Full_tolerances$mean_low_hfi_0.5 <- rowMeans(Full_tolerances[, 151:200], na.rm = TRUE)
Full_tolerances$sd_low_hfi_0.5 <- apply(Full_tolerances[, 151:200], 1, sd, na.rm = TRUE)  

Full_tolerances$mean_low_hfi_0.1 <- rowMeans(Full_tolerances[, 201:250], na.rm = TRUE)
Full_tolerances$sd_low_hfi_0.1 <- apply(Full_tolerances[, 201:250], 1, sd, na.rm = TRUE) 

Full_tolerances$mean_breadth_0.5 <- Full_tolerances$mean_high_hfi_0.5 - Full_tolerances$mean_low_hfi_0.5
Full_tolerances$sd_breadth_0.5 <- apply(Full_tolerances[, 252:301], 1, sd, na.rm = TRUE)

Full_tolerances$mean_breadth_0.1 <- Full_tolerances$mean_high_hfi_0.1 - Full_tolerances$mean_low_hfi_0.1
Full_tolerances$sd_breadth_0.1 <- apply(Full_tolerances[, 302:351], 1, sd, na.rm = TRUE)

Full_tolerances$mean_right_proportion_breadth_0.5 <- (Full_tolerances$mean_high_hfi_0.5 - Full_tolerances$mean_peak)/Full_tolerances$mean_breadth_0.5
#Full_tolerances$sd_right_proportion_breadth_0.5 <- (Full_tolerances$sd_high_hfi_0.5 - Full_tolerances$sd_peak)/Full_tolerances$sd_breadth_0.5

Full_tolerances$mean_right_proportion_breadth_0.1 <- (Full_tolerances$mean_high_hfi_0.1 - Full_tolerances$mean_peak)/Full_tolerances$mean_breadth_0.1
#Full_tolerances$sd_right_proportion_breadth_0.1 <- (Full_tolerances$sd_high_hfi_0.1 - Full_tolerances$sd_peak)/Full_tolerances$sd_breadth_0.1

# code to get names of species in a certain band of the histogram
#row.names(Full_tolerances)[which(Full_tolerances$mean_breadth_0.5 >= 45)]

# Number of species that have no missing bootstraps
sum(Full_tolerances$n_missing == 0) # no missing bootstraps for 142 species

Full_tolerances$margin_peak <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_peak / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$lowerinterval_peak <- Full_tolerances$mean_peak - Full_tolerances$margin_peak
Full_tolerances$upperinterval_peak <- Full_tolerances$mean_peak + Full_tolerances$margin_peak

Full_tolerances$margin_high_hfi_0.5 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_high_hfi_0.5 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$upperinterval_high_hfi_0.5 <- Full_tolerances$mean_high_hfi_0.5 + Full_tolerances$margin_high_hfi_0.5
Full_tolerances$lowerinterval_high_hfi_0.5 <- Full_tolerances$mean_high_hfi_0.5 - Full_tolerances$margin_high_hfi_0.5

Full_tolerances$margin_high_hfi_0.1 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_high_hfi_0.1 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$upperinterval_high_hfi_0.1 <- Full_tolerances$mean_high_hfi_0.1 + Full_tolerances$margin_high_hfi_0.1
Full_tolerances$lowerinterval_high_hfi_0.1 <- Full_tolerances$mean_high_hfi_0.1 - Full_tolerances$margin_high_hfi_0.1

Full_tolerances$margin_low_hfi_0.5 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_low_hfi_0.5 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$upperinterval_low_hfi_0.5 <- Full_tolerances$mean_low_hfi_0.5 + Full_tolerances$margin_low_hfi_0.5
Full_tolerances$lowerinterval_low_hfi_0.5 <- Full_tolerances$mean_low_hfi_0.5 - Full_tolerances$margin_low_hfi_0.5

Full_tolerances$margin_low_hfi_0.1 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_low_hfi_0.1 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$lowerinterval_low_hfi_0.1 <- Full_tolerances$mean_low_hfi_0.1 - Full_tolerances$margin_low_hfi_0.1
Full_tolerances$upperinterval_low_hfi_0.1 <- Full_tolerances$mean_low_hfi_0.1 + Full_tolerances$margin_low_hfi_0.1

Full_tolerances$margin_breadth_0.5 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_breadth_0.5 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$lower_interval_breadth_0.5 <- Full_tolerances$mean_breadth_0.5 - Full_tolerances$margin_breadth_0.5
Full_tolerances$upper_interval_breadth_0.5 <- Full_tolerances$mean_breadth_0.5 + Full_tolerances$margin_breadth_0.5

Full_tolerances$margin_breadth_0.1 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_breadth_0.1 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$lower_interval_breadth_0.1 <- Full_tolerances$mean_breadth_0.1 - Full_tolerances$margin_breadth_0.1
Full_tolerances$upper_interval_breadth_0.1 <- Full_tolerances$mean_breadth_0.1 + Full_tolerances$margin_breadth_0.1

# Skipped for now until definitely needed
Full_tolerances$margin_right_proportion_breadth_0.5 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_right_proportion_breadth_0.5 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$lower_interval_right_proportion_breadth_0.5 <- Full_tolerances$mean_right_proportion_breadth_0.5 - Full_tolerances$margin_right_proportion_breadth_0.5
Full_tolerances$upper_interval_right_proportion_breadth_0.5 <- Full_tolerances$mean_right_proportion_breadth_0.5 + Full_tolerances$margin_right_proportion_breadth_0.5

Full_tolerances$margin_right_proportion_breadth_0.1 <- qt(0.9, df = (50 - Full_tolerances$n_missing) - 1) * Full_tolerances$sd_right_proportion_breadth_0.1 / sqrt((50 - Full_tolerances$n_missing))
Full_tolerances$lower_interval_right_proportion_breadth_0.1 <- Full_tolerances$mean_right_proportion_breadth_0.1 - Full_tolerances$margin_right_proportion_breadth_0.1
Full_tolerances$upper_interval_right_proportion_breadth_0.1 <- Full_tolerances$mean_right_proportion_breadth_0.1 + Full_tolerances$margin_right_proportion_breadth_0.1

Full_tolerances_with_metrics <- Full_tolerances
save(Full_tolerances_with_metrics, file = "C:/Data1/UK/tolerances_with_metrics_uk.RData")

# Final result table format
tolerance_output <- Full_tolerances_with_metrics[, c(251,352:ncol(Full_tolerances_with_metrics))]
tolerance_output$Species <- rownames(tolerance_output)
tolerance_output$Continent <- "UK"

# Set the limits to confidence intervals
tolerance_output$lowerinterval_peak <- ifelse(tolerance_output$lowerinterval_peak < 0, 0, tolerance_output$lowerinterval_peak)
tolerance_output$upperinterval_peak <- ifelse(tolerance_output$upperinterval_peak > 50, 50, tolerance_output$upperinterval_peak)
tolerance_output$upperinterval_high_hfi_0.5 <- ifelse(tolerance_output$upperinterval_high_hfi_0.5 > 50, 50, tolerance_output$upperinterval_high_hfi_0.5)
tolerance_output$lowerinterval_high_hfi_0.5 <- ifelse(tolerance_output$lowerinterval_high_hfi_0.5 < 0, 0, tolerance_output$lowerinterval_high_hfi_0.5)
tolerance_output$upperinterval_high_hfi_0.1 <- ifelse(tolerance_output$upperinterval_high_hfi_0.1 > 50, 50, tolerance_output$upperinterval_high_hfi_0.1)
tolerance_output$lowerinterval_high_hfi_0.1 <- ifelse(tolerance_output$lowerinterval_high_hfi_0.1 < 0, 0, tolerance_output$lowerinterval_high_hfi_0.1)
tolerance_output$upperinterval_low_hfi_0.5 <- ifelse(tolerance_output$upperinterval_low_hfi_0.5 > 50, 50, tolerance_output$upperinterval_low_hfi_0.5)
tolerance_output$lowerinterval_low_hfi_0.5 <- ifelse(tolerance_output$lowerinterval_low_hfi_0.5 < 0, 0, tolerance_output$lowerinterval_low_hfi_0.5)
tolerance_output$lowerinterval_low_hfi_0.1 <- ifelse(tolerance_output$lowerinterval_low_hfi_0.1 < 0, 0, tolerance_output$lowerinterval_low_hfi_0.1)
tolerance_output$upperinterval_low_hfi_0.1 <- ifelse(tolerance_output$upperinterval_low_hfi_0.1 > 50, 50, tolerance_output$upperinterval_low_hfi_0.1)
tolerance_output$lower_interval_breadth_0.5 <- ifelse(tolerance_output$lower_interval_breadth_0.5 < 0, 0, tolerance_output$lower_interval_breadth_0.5)
tolerance_output$upper_interval_breadth_0.5 <- ifelse(tolerance_output$upper_interval_breadth_0.5 > 50, 50, tolerance_output$upper_interval_breadth_0.5)
tolerance_output$lower_interval_breadth_0.1 <- ifelse(tolerance_output$lower_interval_breadth_0.1 < 0, 0, tolerance_output$lower_interval_breadth_0.1)
tolerance_output$upper_interval_breadth_0.1 <- ifelse(tolerance_output$upper_interval_breadth_0.1 > 50, 50, tolerance_output$upper_interval_breadth_0.1)

tolerance_output$ci_range_peak <- tolerance_output$upperinterval_peak - tolerance_output$lowerinterval_peak
tolerance_output$ci_range_high_hfi_0.5 <- tolerance_output$upperinterval_high_hfi_0.5 - tolerance_output$lowerinterval_high_hfi_0.5
tolerance_output$ci_range_high_hfi_0.1 <- tolerance_output$upperinterval_high_hfi_0.1 - tolerance_output$upperinterval_high_hfi_0.1
tolerance_output$ci_range_low_hfi_0.5 <- tolerance_output$upperinterval_low_hfi_0.5 - tolerance_output$lowerinterval_low_hfi_0.5
tolerance_output$ci_range_low_hfi_0.1 <- tolerance_output$upperinterval_low_hfi_0.1 - tolerance_output$lowerinterval_low_hfi_0.1
tolerance_output$ci_range_breadth_0.5 <- tolerance_output$upper_interval_breadth_0.5 - tolerance_output$lower_interval_breadth_0.5
tolerance_output$ci_range_breadth_0.1 <- tolerance_output$upper_interval_breadth_0.1 - tolerance_output$lower_interval_breadth_0.1

# Exclude species with many missing bootstraps
tolerance_output <- tolerance_output[tolerance_output$n_missing == 0, ]

save(tolerance_output, file = "C:/Data1/UK/tolerance_output.RData")

# creating conservation status version
tolerance_output_with_conservation <- tolerance_output

load("C:/Data1/UK/red_species.RData")
load("C:/Data1/UK/amber_species.RData")
load("C:/Data1/UK/green_species.RData")
load("C:/Data1/UK/na_species.RData")

species_list <- rownames(tolerance_output_with_conservation)
tolerance_output_with_conservation$uk_conservation_status <- 0

for (i in 1:length(species_list)) {
  if (species_list[i] %in% red_species){
    tolerance_output_with_conservation$uk_conservation_status[i] <- 1
  }
  else if (species_list[i] %in% amber_species){
    tolerance_output_with_conservation$uk_conservation_status[i] <- 2
  }
  else if (species_list[i] %in% green_species){
    tolerance_output_with_conservation$uk_conservation_status[i] <- 3
  }
  else {
    tolerance_output_with_conservation$uk_conservation_status[i] <- 4
  }
}

load("C:/Data1/UK/bto_species_my_data.RData")

tolerance_output_with_conservation$population_change_27 <- NA
for (i in 1:length(species_list)) {
  if (species_list[i] %in% bto_species$species.scientific){
    tolerance_output_with_conservation$population_change_27[i] <- bto_species$change.27[i]
  }
}

tolerance_output_with_conservation$population_change_10 <- NA
for (i in 1:length(species_list)) {
  if (species_list[i] %in% bto_species$species.scientific){
    tolerance_output_with_conservation$population_change_10[i] <- bto_species$change.10[i]
  }
}

tolerance_output_with_conservation$population_change_5 <- NA
for (i in 1:length(species_list)) {
  if (species_list[i] %in% bto_species$species.scientific){
    tolerance_output_with_conservation$population_change_5[i] <- bto_species$change.5[i]
  }
}

tolerance_output_with_conservation$population_change_1 <- NA
for (i in 1:length(species_list)) {
  if (species_list[i] %in% bto_species$species.scientific){
    tolerance_output_with_conservation$population_change_1[i] <- bto_species$change.1[i]
  }
}

save(tolerance_output_with_conservation, file = "C:/Data1/UK/tolerance_output_with_conservation.RData")

#### CALCULATE SUMMARY STATISTICS ####

# Numbers of species for which calculations were made 
nrow(tolerance_output) # 142
nrow(Full_tolerances[Full_tolerances$n_missing != 50, ]) # 174

# Mean tolerances across all species 
mean(tolerance_output$mean_peak) # 26.68
mean(tolerance_output$mean_high_hfi_0.1) # 46.19
mean(tolerance_output$mean_high_hfi_0.5) # 39.63
mean(tolerance_output$mean_low_hfi_0.1) # 6.93
mean(tolerance_output$mean_low_hfi_0.5) # 12.53
mean(tolerance_output$mean_breadth_0.1) # 39.26
mean(tolerance_output$mean_breadth_0.5) # 27.10
mean(tolerance_output$mean_right_proportion_breadth_0.1) # 0.57
mean(tolerance_output$mean_right_proportion_breadth_0.5) # 0.57

# Correlations among tolerance measures 
cor(tolerance_output$mean_majority, tolerance_output$mean_conservative) 
cor(tolerance_output$mean_majority, tolerance_output$mean_maximum) 
cor(tolerance_output$mean_conservative, tolerance_output$mean_maximum)

# Average breadths of 80% confidence intervals across species 
mean(tolerance_output$ci_range_peak) # 0.67
mean(tolerance_output$ci_range_low_hfi_0.5) # 0.41
mean(tolerance_output$ci_range_low_hfi_0.1) # 0.36

# Number of species with tolerance only to wilderness areas (HFI < 1) 
nrow(tolerance_output[tolerance_output$mean_peak < 1, ]) # 7 
nrow(tolerance_output[tolerance_output$mean_high_hfi_0.1 < 1, ]) # 0
nrow(tolerance_output[tolerance_output$mean_high_hfi_0.5 < 1, ]) # 0

# Number of species with tolerance only to intact areas (HFI < 4)
nrow(tolerance_output[tolerance_output$mean_peak < 4, ]) # 9
nrow(tolerance_output[tolerance_output$mean_high_hfi_0.1 < 4, ]) # 0
nrow(tolerance_output[tolerance_output$mean_high_hfi_0.5 < 4, ]) # 0

# Number of species with tolerance to very high anthropogenic disturbances (HFI > 40) 
nrow(tolerance_output[tolerance_output$mean_peak > 40, ]) # 11
nrow(tolerance_output[tolerance_output$mean_high_hfi_0.1 > 40, ]) # 122
nrow(tolerance_output[tolerance_output$mean_high_hfi_0.5 > 40, ]) # 82
