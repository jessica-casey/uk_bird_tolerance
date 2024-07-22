# RESOLVING TAXONOMIES
#Skip resolving taxonomies for now.

# Load eBird name files
load("C:/Data1/UK/Band1/sp_filtered_Band1.RData")
load("C:/Data1/UK/Band2/sp_filtered_Band2.RData")
load("C:/Data1/UK/Band3/sp_filtered_Band3.RData")
load("C:/Data1/UK/Band4/sp_filtered_Band4.RData")
load("C:/Data1/UK/Band5/sp_filtered_Band5.RData")

load("C:/Data1/UK/Band1/sampling_events_filtered_Band1.RData")
load("C:/Data1/UK/Band2/sampling_events_filtered_Band2.RData")
load("C:/Data1/UK/Band3/sampling_events_filtered_Band3.RData")
load("C:/Data1/UK/Band4/sampling_events_filtered_Band4.RData")
load("C:/Data1/UK/Band5/sampling_events_filtered_Band5.RData")

#Combining datasets
uk_sp_filtered <- rbind(sp_filtered_Band1, sp_filtered_band2, sp_filtered_band3, sp_filtered_band4, sp_filtered_band5)
save(uk_sp_filtered, file = "C:/Data1/UK/uk_sp_filtered.RData")
rm(uk_sp_filtered, sp_filtered_Band1, sp_filtered_band2, sp_filtered_band3, sp_filtered_band4, sp_filtered_band5)

uk_sampling_events_filtered <- rbind(sampling_events_filtered_Band1, sampling_events_filtered_band2, sampling_events_filtered_band3, sampling_events_filtered_band4, sampling_events_filtered_band5)
save(uk_sampling_events_filtered, file = "C:/Data1/UK/uk_sampling_events_filtered.RData")
rm(uk_sampling_events_filtered, sampling_events_filtered_Band1, sampling_events_filtered_band2, sampling_events_filtered_band3, sampling_events_filtered_band4, sampling_events_filtered_band5)

#Getting list of species (would be comparing them with birdlife but here they are all assumed to be fine)
load("C:/Data1/UK/uk_sp_filtered.RData")
species_list_uk <- unique(uk_sp_filtered[, c("scientific_name", "common_name")])
colnames(species_list_uk) <- c("Scientific", "English")

species_list_uk$scientific_birdlife <- species_list_uk$Scientific
colnames(species_list_uk) <- c("scientific_ebird", "english", "scientific_birdlife")

uk_full_name_list <- as.data.frame(species_list_uk)
uk_full_name_list <- uk_full_name_list[, c(2,1,3)]

save(uk_full_name_list, file = "C:/Data1/UK/uk_resolved_taxonomy.RData")