library(ggplot2)
library(dplyr)

#### Big Graph ####
data_for_0.1_plot <- Full_tolerances[, c("mean_low_hfi_0.1", "mean_high_hfi_0.1")]
data_for_0.1_plot <- data_for_0.1_plot[order(data_for_0.1_plot$mean_low_hfi_0.1), ]
data_for_0.1_plot <- as.data.frame.table(data_for_0.1_plot)
data_for_0.1_plot <- data_for_0.1_plot[, -c(2)]
colnames(data_for_0.1_plot) <- c("Species", "Mean_Low_HFI", "Mean_High_HFI")

ggplot(data_for_0.1_plot) +
  geom_segment( aes(x=Species, xend=Species, y=Mean_Low_HFI, yend=Mean_High_HFI), color="black") +
  geom_point( aes(x=Species, y=Mean_Low_HFI, color=("Mean_Low_HFI")), size=2 ) +
  geom_point( aes(x=Species, y=Mean_High_HFI, color=("Mean_High_HFI")), size=2 ) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("orange", "deepskyblue4"), name   = "Key")+
  ylab("HFI") +
  ggtitle("Mean HFI breadth 0.1")


data_for_0.5_plot <- Full_tolerances[, c("mean_low_hfi_0.5", "mean_high_hfi_0.5")]
data_for_0.5_plot <- data_for_0.5_plot[order(data_for_0.5_plot$mean_low_hfi_0.5), ]
data_for_0.5_plot <- as.data.frame.table(data_for_0.5_plot)
data_for_0.5_plot <- data_for_0.5_plot[, -c(2)]
colnames(data_for_0.5_plot) <- c("Species", "Mean_Low_HFI", "Mean_High_HFI")

ggplot(data_for_0.5_plot) +
  geom_segment( aes(x=Species, xend=Species, y=Mean_Low_HFI, yend=Mean_High_HFI), color="black") +
  geom_point( aes(x=Species, y=Mean_Low_HFI, color=("Mean_Low_HFI")), size=2 ) +
  geom_point( aes(x=Species, y=Mean_High_HFI, color=("Mean_High_HFI")), size=2 ) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("orange", "deepskyblue4"), name   = "Key")+
  ylab("HFI") +
  ggtitle("Mean HFI breadth 0.5")


uk_full_name_list <- uk_full_name_list[-c(which(uk_full_name_list$scientific_ebird %in% oceanic_species)), ]
rm(oceanic_species)



#### Initial graphs ####
Full_tolerances <- Full_tolerances_with_metrics
species_conservation_categories <- read.csv("C:/Data1/UK/birds_of_conservation_concern_uk.csv")
red_species <- species_conservation_categories$Red.scientific
amber_species <- species_conservation_categories$Amber.scientific
red_and_amber_species <- c(red_species, amber_species)
green_species <- species_conservation_categories$Green.scientific
na_species <- species_conservation_categories$NA.scientific

# Plotting mean peak against 0.5 breadth comparing red/amber and green species
data_for_red_graph <- Full_tolerances[c(which(rownames(Full_tolerances) %in% red_species)), ]
plot(
  x = data_for_red_graph$mean_peak,
  y = data_for_red_graph$mean_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean breadth",
  main = "Red species, 0.5"
)

data_for_amber_graph <- Full_tolerances[c(which(rownames(Full_tolerances) %in% amber_species)), ]
plot(
  x = data_for_amber_graph$mean_peak,
  y = data_for_amber_graph$mean_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean breadth",
  main = "Amber species, 0.5"
)

data_for_green_graph <- Full_tolerances[c(which(rownames(Full_tolerances) %in% green_species)), ]
plot(
  x = data_for_green_graph$mean_peak,
  y = data_for_green_graph$mean_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean breadth",
  main = "Green species, 0.5"
)

data_for_ra_graph <- Full_tolerances[c(which(rownames(Full_tolerances) %in% red_and_amber_species)), ]
plot(
  x = data_for_ra_graph$mean_peak,
  y = data_for_ra_graph$mean_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean breadth",
  main = "Red and Amber species, 0.5",
)

data_for_na_graph <- Full_tolerances[c(which(rownames(Full_tolerances) %in% na_species)), ]
plot(
  x = data_for_na_graph$mean_peak,
  y = data_for_na_graph$mean_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean breadth",
  main = "NA species, 0.5"
)

# Looking at right proportion
plot(
  x = data_for_ra_graph$mean_high_hfi_0.5,
  y = data_for_ra_graph$mean_right_proportion_breadth_0.5,
  xlab = "Mean high HFI",
  ylab = "Mean right proportion breadth",
  main = "Red and amber species, 0.5"
)

plot(
  x = data_for_rest_graph$mean_peak,
  y = data_for_rest_graph$mean_right_proportion_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean right proportion breadth",
  main = "Green and NA species, 0.5"
)

plot(
  x = data_for_ra_graph$mean_peak,
  y = data_for_ra_graph$mean_right_proportion_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean right proportion breadth",
  main = "Red and amber species, 0.5",
  abline(a = 1, b = -0.02)
)

plot(
  x = data_for_green_graph$mean_peak,
  y = data_for_green_graph$mean_right_proportion_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean right proportion breadth",
  main = "Green species, 0.5",
  abline(a = 1, b = -0.02)
)

# Looking at bto species graphs
plot(
  x = bto_species$change.1,
  y = bto_species$mean_peak,
  xlab = "Population change over 1 year",
  ylab = "Mean HFI peak",
  main = "Mean peak by population change"
)

plot(
  x = bto_species$change.1,
  y = bto_species$mean_breadth_0.5,
  xlab = "Population change over 1 year",
  ylab = "Mean breadth 0.5",
  main = "Mean breadth (0.5) by population change"
)

plot(
  x = bto_species$change.1,
  y = bto_species$mean_right_proportion_breadth_0.5,
  xlab = "Population change over 1 year",
  ylab = "Mean right proportion breadth 0.5",
  main = "Mean right proportion of breadth (0.5) by population change"
)
# Comparing declining and increasing population trends
plot(
  x = bto_5_year_declining$change.5,
  y = bto_5_year_declining$mean_right_proportion_breadth_0.5,
  xlab = "Population change over 5 years",
  ylab = "Mean right proportion breadth 0.5",
  main = "Declining species"
)

plot(
  x = bto_5_year_increasing$change.5,
  y = bto_5_year_increasing$mean_right_proportion_breadth_0.5,
  xlab = "Population change over 5 years",
  ylab = "Mean right proportion breadth 0.5",
  main = "Increasing species"
)



#### Combined coloured graph ####

# Plotting red, amber, green on the same graph in different colours to compare

colors <- c("#CA0020",
            "#E66101",
            "#4DAC26",
            "#BABABA")
plot(
  x = test$mean_peak,
  y = test$mean_breadth_0.5,
  xlab = "Mean peak",
  ylab = "Mean breadth",
  main = "All species, 0.5",
  pch = 19,
  col = colors[factor(test$uk_conservation_status)]
)
