##########################################
## Time series analysis investigating the potntial expansion of Alexandrium pseudogonyaulax in northern European waters 
## by Kristof Möller, Jacob Carstensen, Hans Jakobsen, Annette Engesmo and Bengt Karlson 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 06/24
## Alfred-Wegener-Institute Bremerhaven / International Atomic Energy Agency Monaco
##########################################
####### Load custom functions ####### 
script_dir <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(paste0(script_dir, "/", "Time_series_analysis_custom_functions.R"))

# Install needed packages
install_packages()

####### Load NORWAY Data #######
norway_counts <- # including phytoplankton data
 read_delim(
  paste0(script_dir, "/", "norway", "/", "phyto_aug24.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
 )

# Load physical and chemical parameters of the water column
norway_waterquality <-
 read_delim(
  paste0(script_dir, "/", "norway", "/", "nutrients.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
 )

# Load CTD data
norway_ctd <-
  read_delim(
    paste0(script_dir, "/", "norway", "/", "norway_ctd_full.txt"),
    delim = "\t",
    col_names = TRUE,
    col_types = cols(.default = "c")
  )

####### General data transformations #######
norway_ctd <- norway_ctd %>%
  mutate(across(all_of(c(37:38, 40:46)), ~ as.numeric(gsub(",", ".", .))))

norway_ctd <- 
 pivot_wider(
  norway_ctd[, c(9, 36:46)],
  values_from = verdinum,
  names_from  = Parameter_navn,
  values_fn   = mean
 )

norway_ctd <- norway_ctd %>%
 dplyr::select(
  station     = Vannlokalitet,
  depth_start = depthH,
  depth_end   = depthL,
  date        = dato,
  year,
  month,
  day,
  doy,
  lat         = latitude,
  lon         = longitude,
  temp        = Temperatur,
  sal         = Salinitet,
  -Oksygenmetning,
  -Oksygen,
  -Tetthet
 )

# calculate the density column
norway_ctd <- norway_ctd %>% mutate(temp = as.numeric(temp), sal = as.numeric(sal)) %>%
  mutate(density = calc_seawater_density(temp, sal))

# combine both water depth in one (One sometimes contains NAs)
norway_ctd <-
 norway_ctd %>% 
  mutate(depth = rowMeans(dplyr::select(., one_of(
  "depth_start", "depth_end"
 )), na.rm = TRUE)) %>% 
  dplyr::select(-depth_start, -depth_end)

# Introduce stratification key if density difference of first two metres vs last two metres exceeds 1 g/cm^3
norway_ctd <- strat_index(norway_ctd)

# Average physical parameters of CTD data over first 10m of the water column
norway_ctd <- norway_ctd %>%
 group_by(date, year, doy, day, month, strat, station, lat, lon) %>%
 filter(depth < 10) %>%
 dplyr::summarise(
  sal = mean(sal, na.rm = TRUE),
  temp = mean(temp, na.rm = TRUE),
  density = mean(density, na.rm = TRUE)
 ) %>%
 ungroup()

norway_counts <- norway_counts %>%
 dplyr::select(
  station  = lok_kode_new_tlj,
  date     = dato,
  year,
  month,
  day,
  species  = VitenskapligNavn,
  cells_L  = celletall,
  doy,
  lat      = latitude,
  lon      = longitude
 )

norway_waterquality <- left_join(
 norway_waterquality,
 norway_counts %>% 
   dplyr::select(station, lat, lon) %>% 
   unique(),
 by = c("Vannlokalitet" = "station")
)

# Classify into 11 harmful algal genera, assign toxin syndromes, and generate probability keys
norway_counts <- norway_counts %>%
  mutate(cells_L = as.numeric(cells_L)) %>%
  process_harmful_genera_comprehensive()

# Remove any absent entries on dates where a harmful genus was actually present
# stems from the operation before that changed all other species to "No harmful algae"
# thus on all present days "No harmful algae" exists as well and needs to be removed
harmful_algae_intermediate <- norway_counts %>%
  filter(!species %in% c("No harmful algae") & harmful_algae == TRUE)

norway_counts <- filter_harmful_algae(norway_counts, species_col = "species")

norway_counts <- full_join(norway_counts, harmful_algae_intermediate)

norway_waterquality <- norway_waterquality %>%
  mutate(across(all_of(c(2:11, 13:14)), ~ as.numeric(.))) %>%
 dplyr::rename(
  station  = Vannlokalitet,
  PO4      = phospho,
  NO3      = nitro,
  silicate = sil,
  date     = dato
 ) 

####### Merge all Norwegian data files ####### 
# Average abiotic parameter dataframes by date and station as they occasionally have two measurements
norway_waterquality_sub <- norway_waterquality %>% 
  group_by(station, date) %>%  
  dplyr::select(-doy, -day, -month, -year, -lat, -lon) %>%
  reframe(across(c(1:6), \(x) mean(x, na.rm = TRUE))) 

norway_ctd_sub <- norway_ctd %>% 
  group_by(station, date, strat) %>%  
  dplyr::select(-doy, -day, -month, -year, -lat, -lon) %>%
  reframe(across(c(1:3), \(x) mean(x, na.rm = TRUE))) 

# Full join both abiotic parameter dataframes
norway_abiotic <- full_join(norway_waterquality_sub,
                            norway_ctd_sub,
                            by = c("date", "station"))

####### Loop over all unique stations and join abiotic and phytoplankton dataframes with a 1 day threshold ####### 
all_stations_norway <- unique(c(unique(norway_counts$station), unique(norway_abiotic$station)))
norway_list <- list()

norway_counts <- norway_counts %>%
  mutate(across(all_of(c(3:5, 7:10)), ~ as.numeric(.)), date = as.Date(date)) 

norway_abiotic <- norway_abiotic %>% mutate(date = as.Date(date))

for(each_station in all_stations_norway){
  
  norway_counts_sub <- norway_counts %>% 
    filter(station == each_station) 
  
  norway_abiotic_sub <- norway_abiotic%>% 
    filter(station == each_station) 

  norway_combined <- difference_full_join(
    norway_counts_sub,
    norway_abiotic_sub,
    by = c("date" = "date"),
    max_dist = 1
  )  %>%
    mutate(
      station = coalesce(station.x, station.y),
      date = coalesce(date.x, date.y)
    ) %>%
    dplyr::select(-matches("\\.(x|y)$"))

  norway_list[[each_station]] <- norway_combined
}

norway_combined <- do.call(rbind, norway_list)

####### Finalize dataframe ####### 
unique_stations <- norway_combined %>% group_by(station) %>%
  dplyr::summarise(lat = mean(lat, na.rm = T), 
                   lon = mean(lon, na.rm = T))

norway_combined <- full_join(norway_combined %>% ungroup() %>% dplyr::select(-lat, -lon),
                             unique_stations,
                             by = "station")

norway_combined <- update_dates(norway_combined %>% drop_na(date))

# Recompute overall probability (1 = present, 0 = absent) after joining with abiotic data
norway_combined <- norway_combined %>%
  mutate(
    probability = case_when(
      harmful_algae == TRUE & !is.na(cells_L) & cells_L > 0 ~ 1L,
      is.na(species) | is.na(cells_L)                        ~ NA_integer_,
      TRUE                                                    ~ 0L
    ),
    month = month(date)
  )

# Step 1: Aggregate genus-level probabilities, cell concentrations and categorical fields per date/station
norway_combined_probs <- norway_combined %>%
  group_by(date, station) %>%
  summarise(
    # Overall: 1 if any row is present, 0 if all absent, NA if all missing
    probability = case_when(
      any(probability == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability))              ~ NA_integer_,
      TRUE                                 ~ 0L
    ),
    # Individual genus probabilities (1/0/NA)
    probability_Alexandrium = case_when(
      any(probability_Alexandrium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Alexandrium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Dinophysis = case_when(
      any(probability_Dinophysis == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Dinophysis))              ~ NA_integer_, TRUE ~ 0L),
    probability_Pseudonitzschia = case_when(
      any(probability_Pseudonitzschia == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Pseudonitzschia))              ~ NA_integer_, TRUE ~ 0L),
    probability_Azadinium = case_when(
      any(probability_Azadinium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Azadinium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Chrysochromulina = case_when(
      any(probability_Chrysochromulina == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Chrysochromulina))              ~ NA_integer_, TRUE ~ 0L),
    probability_Prymnesium = case_when(
      any(probability_Prymnesium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Prymnesium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Amphidinium = case_when(
      any(probability_Amphidinium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Amphidinium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Pseudochattonella = case_when(
      any(probability_Pseudochattonella == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Pseudochattonella))              ~ NA_integer_, TRUE ~ 0L),
    probability_Phaeocystis = case_when(
      any(probability_Phaeocystis == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Phaeocystis))              ~ NA_integer_, TRUE ~ 0L),
    probability_Karlodinium = case_when(
      any(probability_Karlodinium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Karlodinium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Cyanobacteria = case_when(
      any(probability_Cyanobacteria == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Cyanobacteria))              ~ NA_integer_, TRUE ~ 0L),
    # Genus-level cell concentration totals per date/station
    cells_L_Alexandrium      = sum(ifelse(genus == "Alexandrium",       cells_L, 0), na.rm = TRUE),
    cells_L_Dinophysis       = sum(ifelse(genus == "Dinophysis",        cells_L, 0), na.rm = TRUE),
    cells_L_Pseudonitzschia  = sum(ifelse(genus == "Pseudo-nitzschia",  cells_L, 0), na.rm = TRUE),
    cells_L_Azadinium        = sum(ifelse(genus == "Azadinium",         cells_L, 0), na.rm = TRUE),
    cells_L_Chrysochromulina = sum(ifelse(genus == "Chrysochromulina",  cells_L, 0), na.rm = TRUE),
    cells_L_Prymnesium       = sum(ifelse(genus == "Prymnesium",        cells_L, 0), na.rm = TRUE),
    cells_L_Amphidinium      = sum(ifelse(genus == "Amphidinium",       cells_L, 0), na.rm = TRUE),
    cells_L_Pseudochattonella = sum(ifelse(genus == "Pseudochattonella", cells_L, 0), na.rm = TRUE),
    cells_L_Phaeocystis      = sum(ifelse(genus == "Phaeocystis",       cells_L, 0), na.rm = TRUE),
    cells_L_Karlodinium      = sum(ifelse(genus == "Karlodinium",       cells_L, 0), na.rm = TRUE),
    cells_L_Cyanobacteria    = sum(ifelse(genus == "Cyanobacteria",     cells_L, 0), na.rm = TRUE),
    # Keep first value of categorical/date-component fields
    across(c(day, month, year, doy, strat, species, genus, toxin_syndrome, harmful_algae), first),
    .groups = "drop"
  ) %>%
  mutate(
    logistic = probability,  # alias kept for downstream compatibility
    country  = "Norway"
  )

# Step 2: Average all numeric environmental variables per date/station
norway_combined_num <- norway_combined %>%
  group_by(date, station) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::select(-any_of(c("day", "month", "year", "doy", "logistic",
    "cells_L_Alexandrium", "cells_L_Dinophysis", "cells_L_Pseudonitzschia",
    "cells_L_Azadinium", "cells_L_Chrysochromulina", "cells_L_Prymnesium",
    "cells_L_Amphidinium", "cells_L_Pseudochattonella", "cells_L_Phaeocystis",
    "cells_L_Karlodinium", "cells_L_Cyanobacteria")))

# Step 3: Join aggregated probabilities with numeric environmental summary
norway_combined <- left_join(norway_combined_probs, norway_combined_num,
                              by = c("date", "station"))

# Step 4: Compute risk classification levels from aggregated genus-specific cell concentrations
norway_combined <- norway_combined %>%
  mutate(
    risk_level_Alexandrium     = case_when(cells_L_Alexandrium     >= 10000 ~ "High", cells_L_Alexandrium     >= 1000 ~ "Moderate", cells_L_Alexandrium     > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Dinophysis      = case_when(cells_L_Dinophysis      >= 10000 ~ "High", cells_L_Dinophysis      >= 1000 ~ "Moderate", cells_L_Dinophysis      > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Pseudonitzschia = case_when(cells_L_Pseudonitzschia >= 10000 ~ "High", cells_L_Pseudonitzschia >= 1000 ~ "Moderate", cells_L_Pseudonitzschia > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Azadinium       = case_when(cells_L_Azadinium       >= 10000 ~ "High", cells_L_Azadinium       >= 1000 ~ "Moderate", cells_L_Azadinium       > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Chrysochromulina = case_when(cells_L_Chrysochromulina >= 10000 ~ "High", cells_L_Chrysochromulina >= 1000 ~ "Moderate", cells_L_Chrysochromulina > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Prymnesium      = case_when(cells_L_Prymnesium      >= 10000 ~ "High", cells_L_Prymnesium      >= 1000 ~ "Moderate", cells_L_Prymnesium      > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Amphidinium     = case_when(cells_L_Amphidinium     >= 10000 ~ "High", cells_L_Amphidinium     >= 1000 ~ "Moderate", cells_L_Amphidinium     > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Pseudochattonella = case_when(cells_L_Pseudochattonella >= 10000 ~ "High", cells_L_Pseudochattonella >= 1000 ~ "Moderate", cells_L_Pseudochattonella > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Phaeocystis     = case_when(cells_L_Phaeocystis     >= 10000 ~ "High", cells_L_Phaeocystis     >= 1000 ~ "Moderate", cells_L_Phaeocystis     > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Karlodinium     = case_when(cells_L_Karlodinium     >= 10000 ~ "High", cells_L_Karlodinium     >= 1000 ~ "Moderate", cells_L_Karlodinium     > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Cyanobacteria   = case_when(cells_L_Cyanobacteria   >= 10000 ~ "High", cells_L_Cyanobacteria   >= 1000 ~ "Moderate", cells_L_Cyanobacteria   > 0 ~ "Low", TRUE ~ "None")
  )

# Save norway_combined as a new txt.file
write.table(norway_combined,
      file = paste0(script_dir, "/", "norway_combined.txt"),
      sep = "\t",
      row.names = FALSE)

####### DANISH Data ####### 
# Phytoplankton counts data
denmark_counts <-
  read_delim(
    file = paste0(script_dir, "/", "denmark", "/", "DK_counts.txt"),
    delim = "\t",
    col_names = TRUE,
    col_types = cols(.default = "c")
  )

denmark_counts <- denmark_counts %>% 
 mutate(date = as.Date(paste(day, month, year, sep = "."), format = "%d.%m.%Y"))

# CTD data
denmark_ctd <-
 read_delim(
  file = paste0(script_dir, "/", "denmark", "/", "DK_CTD.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
 )

denmark_ctd <- denmark_ctd %>% 
 mutate(date = as.Date(date, format = "%d. %b %y")) %>% 
 mutate(day = day(date))

# Chemical and physical water properties
denmark_waterquality <-
 read_delim(
  file = paste0(script_dir, "/", "denmark", "/", "DK_water_quality.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
 )

denmark_waterquality <-
 denmark_waterquality %>% 
 mutate(date = as.Date(date, format = "%d. %b %y")) %>% 
 mutate(day = day(date))

# Secci depth and Kd values
denmark_secci_kd <-
 read_delim(
  file = paste0(script_dir, "/", "denmark", "/", "DK_secci_kd.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
 )

denmark_secci_kd <-
 denmark_secci_kd %>% 
 mutate(date = as.Date(date, format = "%d. %b %y")) %>% 
 mutate(day = day(date))

# remove the -000X from the Limfjord stations as it does not match the synthax of denmark_counts
denmark_secci_kd   <- remove_station_suffix(denmark_secci_kd)
denmark_waterquality <- remove_station_suffix(denmark_waterquality)
denmark_ctd     <- remove_station_suffix(denmark_ctd)

####### General data transformations ####### 
denmark_ctd <- denmark_ctd %>%
 dplyr::rename(
  depth = `depth (m)`,
  sal = `Salinity (ppt)`,
  temp = `°C`,
  Chl = `Chlorophyl (flourometer,µg Chl L-1)`
 )

denmark_counts <-
 denmark_counts %>% dplyr::rename(
  lon = Longitude,
  lat = Latitude,
  species = species_name,
  cells_L = cells_per_L,
 )

# Replace lat/lon of the stations in denmark_ctd and denmark_waterquality with lat/lon of the stations in denmark_counts
# get lat lon of all stations in denmark_counts
all_stations <-
 denmark_counts %>% 
 dplyr::select(station, lon, lat) %>% 
 unique() %>%
 separate_rows(station, sep = ", ") %>%
 distinct() %>%
 group_by(station) %>%
 mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>%
 dplyr::reframe(lat = mean(lat), lon = mean(lon))

# Apply to each denmark df
denmark_ctd          <- join_station_metadata(denmark_ctd, all_stations)
denmark_waterquality <- join_station_metadata(denmark_waterquality, all_stations)
denmark_secci_kd     <- join_station_metadata(denmark_secci_kd, all_stations)

# calculate density
denmark_ctd <- denmark_ctd %>%
  mutate(temp = as.numeric(temp), sal = as.numeric(sal)) %>%
  mutate(density = calc_seawater_density(temp, sal))

# Introduce stratification key if density difference exceeds 1 g/cm^3
denmark_ctd <- strat_index(denmark_ctd)

# Average physical parameters of CTD data over the whole water column (0-10m)
denmark_ctd <- denmark_ctd %>%
 group_by(date, year, month, day, strat, station, lat, lon) %>%
  mutate(across(all_of(c(1:5)), ~ as.numeric(.))) %>%
         filter(depth < 10) %>%
 dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>%
 ungroup()

# General data transformations: Microalgae counts data ####### 
# Replace underscores with spaces in species names and classify into 11 harmful genera
denmark_counts$species <- str_replace_all(denmark_counts$species, "_", " ")

denmark_counts <- denmark_counts %>%
  mutate(cells_L = as.numeric(cells_L)) %>%
  process_harmful_genera_comprehensive()

# Aggregate carbon data per station and date first
denmark_carbon <- denmark_counts %>%
  dplyr::select(station, date, `C_(µgC_L-1)`) %>%
  group_by(station, date) %>%
  mutate(`C_(µgC_L-1)` = as.numeric(`C_(µgC_L-1)`)) %>%
  dplyr::summarise(C = mean(`C_(µgC_L-1)`))

# Remove duplicate "No harmful algae" entries on dates where a harmful genus was actually present
harmful_algae_intermediate <- denmark_counts %>%
  filter(!species %in% c("No harmful algae") & harmful_algae == TRUE)

denmark_counts <- filter_harmful_algae(denmark_counts, species_col = "species")

denmark_counts <- full_join(denmark_counts, harmful_algae_intermediate)

denmark_counts <- left_join(denmark_counts, denmark_carbon, by = c("station", "date"))

# summarize all numeric columns
denmark_ctd <- denmark_ctd %>%
 group_by(station, date, strat) %>%
  mutate(year = as.numeric(year), month = as.numeric(month)) %>%
 dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

denmark_waterquality <- denmark_waterquality %>%
 group_by(station, date, lat, lon) %>%
  mutate(across(everything(), ~ as.numeric(.))) %>%
 dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

# Merge all data files #####
denmark_secci_kd <- denmark_secci_kd %>% mutate(across(all_of(c("year", "month", "Kd", "Secci", "day")), ~ as.numeric(.)))
denmark_counts <- denmark_counts %>% mutate(across(all_of(c("year", "month", "lon", "lat", "day", "Cell_vol(µm3)", "C_(µgC_L-1)")), ~ as.numeric(.)))

denmark_combined <- denmark_counts %>%
 full_join(denmark_waterquality) %>%
 full_join(denmark_ctd) %>%
 full_join(denmark_secci_kd)

####### Merge all Danish data files ####### 
# Average abiotic parameter dataframes by date and station as they occassionally have two measurements
denmark_counts_sub <- denmark_counts %>%
  group_by(station, date, species, harmful_algae, genus, toxin_syndrome, across(starts_with("probability"))) %>%
  reframe(across(where(is.numeric), mean, na.rm = TRUE))

denmark_waterquality_sub <- denmark_waterquality %>% ungroup() %>%
  dplyr::select(-day, -month, -year, -lat, -lon)

denmark_ctd_sub <- denmark_ctd %>% ungroup() %>%
  dplyr::select(-day, -month, -year, -lat, -lon) 

denmark_secci_kd_sub <- denmark_secci_kd %>% 
  group_by(station, date) %>%  
  dplyr::select(-day, -month, -year, -lat, -lon) %>%
  reframe(across(where(is.numeric), mean, na.rm = TRUE))

# Full join both abiotic parameter dataframes
denmark_abiotic <- full_join(
  denmark_waterquality_sub,
  full_join(
    denmark_ctd_sub,
    denmark_secci_kd_sub,
    by = c("date", "station")),
    by = c("date", "station")
  )
  
####### Loop over all unique stations and join abiotic and phytoplankton dataframes with a 1 day threshold ####### 
all_stations_denmark <- unique(c(unique(denmark_counts_sub$station), unique(denmark_abiotic$station)))
denmark_list <- list()

for(each_station in all_stations_denmark){
  
  denmark_counts_sub2 <- denmark_counts_sub %>% 
    filter(station == each_station) 
  
  denmark_abiotic_sub <- denmark_abiotic %>% 
    filter(station == each_station) 
  
  denmark_combined <- difference_full_join(
    denmark_counts_sub2,
    denmark_abiotic_sub,
    by = c("date" = "date"),
    max_dist = 1
  )  %>%
    mutate(
      station = coalesce(station.x, station.y),
      date = coalesce(date.x, date.y)
    ) %>%
    dplyr::select(-matches("\\.(x|y)$"))
  
  denmark_list[[each_station]] <- denmark_combined
}

denmark_combined <- do.call(rbind, denmark_list)

####### Finalize dataframe ####### 
unique_stations <- full_join(
  full_join(
    full_join(
      denmark_counts %>% group_by(station) %>%
        dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T)),
      denmark_waterquality %>% group_by(station) %>%
        dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))
    ),
    denmark_ctd %>% group_by(station) %>%
      dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))
  ),
  denmark_secci_kd %>% group_by(station) %>%
    dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))
) %>%
  group_by(station) %>%
  dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))

denmark_combined <- full_join(denmark_combined %>% ungroup() %>% dplyr::select(-lat, -lon),
                             unique_stations,
                             by = "station")

denmark_combined <- update_dates(denmark_combined %>% drop_na(date))

####### Read wind data files ####### 
denmark_wind <-
 read_delim(
  paste0(script_dir, "/", "denmark", "/", "wind_data.txt"),
  delim = "\t",
  col_names = TRUE
 )

denmark_wind_stations <-
 fread(
  paste0(script_dir, "/", "denmark", "/", "stationen_wind.txt"),
  fill = TRUE
 )

# General data transformations
denmark_wind <- denmark_wind %>% 
 mutate(
 date =
  as.Date(paste(day, month, year, sep = "-"), format = "%d-%m-%Y"),
 doy = yday(date),
 `wind_m/s` = as.numeric(gsub(",", ".", `wind_m/s`))
) %>% 
 dplyr::select(-release_date)

# remove last two digits of the station name in the denmark_wind dataset to match the format in the station file
denmark_wind <- denmark_wind %>%
 mutate_at(vars(station), ~ as.numeric(str_sub(., end = -3)))

# select columns of interest
denmark_wind_stations <- denmark_wind_stations %>%
 dplyr::select(station, Lat, Lon)

# Merge both wind dataframes
denmark_wind <- full_join(denmark_wind, denmark_wind_stations) %>%
 drop_na(Lat)

# Find the two closest stations between danish monitoring stations and meterological "wind" stations
closest_stations <-
 find_closest_stations(
  denmark_wind   %>% distinct(station, .keep_all = TRUE),
  denmark_counts %>% distinct(station, .keep_all = TRUE)
 )

# Merge denmark_combined with denmark_wind and the closest stations
colnames(closest_stations) <- c("wind_match", "station", "distance")

denmark_combined <- merge(denmark_combined, closest_stations, by = "station", all.x = T)

colnames(denmark_wind)[1] <- "wind_match"

####### Merge denmark_combined and denmark wind ####### 
denmark_combined <-
 merge(
  denmark_combined,
  denmark_wind %>% dplyr::select(wind_match, date, year, month, day, `wind_m/s`),
  by = c("wind_match", "date", "year", "month", "day"),
  all.x = T,
  all.y = T
 )

# Introduce overall probability (1/0/NA) for all harmful genera
denmark_combined <- denmark_combined %>%
  mutate(
    probability = case_when(
      harmful_algae == TRUE & !is.na(cells_L) & cells_L > 0 ~ 1L,
      is.na(species) | is.na(cells_L)                        ~ NA_integer_,
      TRUE                                                    ~ 0L
    ),
    logistic = probability,  # alias kept for downstream compatibility
    country  = "Denmark"
  )

# Rename columns
denmark_combined <- denmark_combined %>% 
  dplyr::rename(
 # C = "C_(µgC_L-1)",
 DIN = "DIN_(µgL-1)",
 NH4 = "NH4_(µgL-1)",
 NO3 = "nitrat+nitrit_(µgL-1)",
 TN = "Total_nitrogen_(µgL-1)",
 DIP = `DIP_(µgL-1)`,
 PO4 = `Total_phosphate_(µgL-1)`,
 chl = `Chlorophyll_a_(µgr_L-1)`,
 silicate = `Silicate_(µgL-1)`,
 chl_2 = "Chl",
 kd = "Kd",
 secci = "Secci",
 wind_ms = `wind_m/s`
) %>%
 dplyr::select(-'Cell_vol(µm3)', -distance, -wind_match)

# Convert ug/L to umol/L to match other datasets
denmark_combined <-
 denmark_combined %>% mutate(
  C = C / 12.0107,
  DIN = DIN / 14.006720,
  NH4 = NH4 / 18.0383,
  NO3 = NO3 / 62.0049,
  TN = TN / 14.006720,
  DIP = DIP / 30.973762 ,
  PO4 = PO4 / 94.9712
 ) %>% drop_na(station, date) 

write.table(denmark_combined,
      paste0(script_dir, "/", "denmark_combined.txt"),
      sep = "\t",
      row.names = FALSE)

####### SWEDISH DATA #######
sweden_phys <-
 read_delim(
  paste0(script_dir, "/", "sweden", "/", "sharkweb_phys2.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
 )

sweden_phys_2024 <-
  read_delim(
    paste0(script_dir, "/", "sweden", "/", "sharkweb_phys_2024.txt"),
    delim = "\t",
    col_names = TRUE,
    col_types = cols(.default = "c")
  )

sweden_phys <- full_join(sweden_phys, sweden_phys_2024)

sweden_phyto <-
 read_delim(
  paste0(script_dir, "/", "sweden", "/", "sharkweb_phyto_new.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
 )

sweden_phyto_2024 <- read_delim(
  paste0(script_dir, "/", "sweden", "/", "sharkweb_phyto_2024.txt"),
  delim = "\t",
  col_names = TRUE,
  col_types = cols(.default = "c")
)

sweden_phyto <- full_join(sweden_phyto, sweden_phyto_2024)

####### General Data Preparation ####### 
sweden_phyto <- sweden_phyto %>%
 dplyr::select(
  station   = reported_station_name,
  date      = sample_date,
  lat       = sample_latitude_dm,
  lon       = sample_longitude_dm,
  wind_ms   = wind_speed_ms,
  secci     = secchi_depth_m,
  min_depth = sample_min_depth_m,
  max_depth = sample_max_depth_m,
  name      = scientific_name,
  parameter,
  counts    = value,
  unit,
  counts_coefficient   = coefficient,
  sedimentation_volume = sedimentation_volume_ml
 )

sweden_phyto <- sweden_phyto %>%
 mutate(across(
  c(
   wind_ms,
   secci,
   min_depth,
   max_depth,
   name,
   parameter,
   counts,
   unit,
   counts_coefficient,
   sedimentation_volume
  ),
  ~ gsub(",", ".", .)
 ))

sweden_phyto <- sweden_phyto %>%
 mutate(counts = as.numeric(counts),
     counts_coefficient = as.numeric(counts_coefficient))

# Perform the multiplication conditionally
sweden_phyto <- sweden_phyto %>%
 mutate(counts = case_when(
  str_detect(parameter, regex("counted", ignore_case = TRUE)) &
   !is.na(counts_coefficient) ~ counts * counts_coefficient,
  TRUE ~ counts
 ))

sweden_phyto <- sweden_phyto %>%
 mutate(
  lat  = sapply(lat, dmm_to_dd),
  lon  = sapply(lon, dmm_to_dd),
  date = as.Date(date)
 )

# Average over the first 10m of the water column
sweden_phyto <-
 sweden_phyto %>% 
 dplyr::rename(cells_L = counts, species = name)

sweden_phyto[, c(5:8)] <-
 sapply(sweden_phyto[, c(5:8)], function(col)
  as.numeric(col))

sweden_phyto <-
 sweden_phyto %>% filter(max_depth <= 10) %>%
 group_by(station, date, species, lat, lon, parameter) %>%
 reframe(
  wind_ms = mean(wind_ms, na.rm = T),
  cells_L = mean(cells_L, na.rm = T),
  secci  = mean(secci, na.rm = T)
 ) %>%
 filter(str_detect(parameter, c("Abundance|counted"))) %>%
 filter(parameter != "Abundance class")

# Classify into 11 harmful genera, assign toxin syndromes, and generate probability keys
sweden_phyto <- sweden_phyto %>%
  mutate(cells_L = as.numeric(cells_L)) %>%
  process_harmful_genera_comprehensive()

# Remove duplicate "No harmful algae" entries on dates where a harmful genus was actually present
harmful_algae_intermediate <- sweden_phyto %>%
  filter(!species %in% c("No harmful algae") & harmful_algae == TRUE)

sweden_phyto <- filter_harmful_algae(sweden_phyto, species_col = "species")

sweden_phyto <- full_join(sweden_phyto, harmful_algae_intermediate)

sweden_phys <- sweden_phys %>%
 mutate(
  lat = sapply(`Provets latitud (DM)`, dmm_to_dd),
  lon = sapply(`Provets longitud (DM)`, dmm_to_dd),
  date = as.Date(Provtagningsdatum, format = "%Y-%m-%d")
 ) %>%
 dplyr::select(
  station        = `Stationsnamn`,
  date,
  lat,
  lon,
  wind_direction = `Vindriktning (kod)`,
  wind_speed     = `Vindhastighet (m/s)`,
  air_temp       = `Lufttemperatur (C)`,
  secci          = `Siktdjup (m)`,
  depth          = `Provtagningsdjup (m)`,
  pressure       = `Tryck CTD (dbar)`,
  temp           = `Temperatur vattenhämtare (C)`,
  temp2          = `Temperatur CTD (C)`,
  sal            = `Salinitet vattenhämtare (o/oo psu)`,
  sal2           = `Salinitet CTD (o/oo psu)`,
  PO4            = `Fosfatfosfor PO4-P (umol/l)`,
  TP             = `Total fosfor Tot-P (umol/l)`,
  NO2            = `Nitritkväve NO2-N (umol/l)`,
  NO3            = `Nitratkväve NO3-N (umol/l)`,
  `NO2+NO3`      = `Nitrit+Nitratkväve NO2+NO3-N (umol/l)`,
  NH4            = `Ammoniumkväve NH4-N (umol/l)`,
  TN             = `Total kväve Tot-N (umol/l)`,
  silicate       = `Silikat SiO3-Si (umol/l)`,
  chl            = `Klorofyll-a vattenhämtare (ug/l)`,
  DOC            = `Löst Organiskt Kol DOC (umol/l)`,
  POC            = `Partikulärt Organiskt Kol POC (umol/l)`,
  TOC            = `Total Organiskt Kol TOC (mg/l)`,
  PON            = `Partikulärt Organiskt Kväve PON (umol/l)`,
  cur_dir        = `Strömriktning (dekagrader)`,
  cur            = `Strömhastighet (m/s)`,
  CDOM           = `CDOM (1/m)`
 ) 

# Convert columns with commas to numeric and change , to .
sweden_phys[, c(3:28)] <-
 sapply(sweden_phys[, c(3:28)], function(col)
  as.numeric(gsub(",", ".", col)))

# calculate the mean of both salinity and temperature columns excluding NAs and remove the redundant ones
sweden_phys <- sweden_phys %>%
 ungroup() %>%
 mutate(
  temp      = rowMeans(dplyr::select(., one_of("temp", "temp2")), na.rm = TRUE),
  sal       = sal,
  `NO2+NO3` = coalesce(NO3, 0) + coalesce(NO2, 0)
 ) %>%
 dplyr::select(-temp2, -sal2) %>%
 mutate(`NO2+NO3` = replace(`NO2+NO3`, is.na(NO3) &
                is.na(NO2), NA)) %>%
 mutate(`NO2+NO3` = replace(`NO2+NO3`, NO2 == 0 &
                is.na(NO3), NA)) %>%
 mutate(DIN = NH4 + `NO2+NO3`) %>% # Calculate the DIN (Dissolved Inorganic Nitrogen) concentration
 dplyr::select(-NO2, -NO3) %>%
 dplyr::rename(NO3 = 'NO2+NO3')

# calculate the density column
sweden_phys <- sweden_phys %>%
 mutate(density = calc_seawater_density(temp, sal))

# Introduce stratification key if density difference exceeds 1 g/cm^3 (?!)
sweden_phys <- strat_index(sweden_phys)

# Average the physical parameters of interest over a sampling depth of 10m to mimic danish data set
# in chunks of the dataframe because my laptop is ******* slow
# Calculate number of chunks
chunk_size <- 25000
num_chunks <- ceiling(nrow(sweden_phys) / chunk_size)

# Initialize an empty list to store the summarized dataframes
summarized_dfs <- list()

# Loop through each chunk, summarize it, and store the result
for (i in 1:num_chunks) {
 start_row <- (i - 1) * chunk_size + 1
 end_row <- min(i * chunk_size, nrow(sweden_phys))
 
 chunk <- sweden_phys[start_row:end_row,]
 
 summarized_chunk <- chunk %>%
  drop_na(station) %>%
  group_by(station, date, strat) %>%
  filter(depth < 10) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
 
 # Store the summarized chunk in the list
 summarized_dfs[[i]] <- summarized_chunk
}

# Combine all summarized dataframes
sweden_phys <- do.call(rbind, summarized_dfs)

# Merge microalgae and physical parameters datasets
sweden_abiotic <- sweden_phys %>% 
  dplyr::select(-lat, -lon)

# Loop over all unique stations and join abiotic and phytoplankton dataframes with a 1 day threshold ####
# Detect number of cores and register parallel backend
# Get list of all unique stations
all_stations_sweden <- unique(c(unique(sweden_phyto$station), unique(sweden_abiotic$station)))

sweden_phyo <- sweden_phyto %>% mutate(date = as.Date(date))
sweden_abiotic <- sweden_abiotic %>% mutate(date = as.Date(date))

sweden_list <- list()
for(each_station in all_stations_sweden){
  index <- which(all_stations_sweden == each_station)
  print(paste("Processing station", each_station, "at index", index))
  
  sweden_phyto_sub <- sweden_phyto %>% 
    filter(station == each_station) 
  
  sweden_abiotic_sub <- sweden_abiotic %>% 
    filter(station == each_station) 
  
  sweden_combined <- difference_full_join(
    sweden_phyto_sub,
    sweden_abiotic_sub,
    by = c("date" = "date"),
    max_dist = 1
  ) %>%
    mutate(
      station = coalesce(station.x, station.y),
      date = coalesce(date.x, date.y)
    ) %>%
    select(-matches("\\.(x|y)$"))
  
  sweden_list[[each_station]] <- sweden_combined
}

sweden_combined <- do.call(rbind, sweden_list)

# Update lat, lon and day, doy, month, year ####
unique_stations <- full_join(
  sweden_phys %>% group_by(station) %>%
    dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T)),
  sweden_phyto %>% group_by(station) %>%
    dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))
) %>% group_by(station) %>%
  dplyr::summarise(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))

sweden_combined <- full_join(sweden_combined %>% ungroup() %>% dplyr::select(-lat, -lon),
                             unique_stations,
                             by = "station")

sweden_combined <- update_dates(sweden_combined %>% drop_na(date))

sweden_combined <- sweden_combined %>%
  mutate(
    wind_ms = rowMeans(dplyr::select(., starts_with("wind")), na.rm = TRUE),
    country = "Sweden"
  ) %>%
  dplyr::select(-wind_speed) %>%
  update_dates() %>%
  filter(is.na(parameter) | parameter != "Abundance")

# Save sweden_combined as a new txt.file
write.table(sweden_combined,
      file = paste0(script_dir, "/", "sweden_combined.txt"),
      sep = "\t",
      row.names = FALSE)

####### IOW-Odin data ####### 
germany_combined <- read_delim(
  file = paste0(script_dir, "/", "germany", "/", "odin2_2025-09-12_070935_red.txt"),
  delim = "\t",
  skip = 2,
  col_names = FALSE,
  col_types = cols(.default = "c")
)

# Read the first two lines
first_two_lines <- read_lines(
 file = paste0(script_dir, "/", "germany", "/", "odin2_2024-01-31_095801.txt"),
 n_max = 2
)

# Split the lines into vectors based on tabs
line1_vector <- strsplit(first_two_lines[1], "\t")[[1]]
line2_vector <- strsplit(first_two_lines[2], "\t")[[1]]

# Concatenate corresponding values from line1 and line2
combined_header <- paste(line1_vector, line2_vector, sep = "_")

# Set the combined header to the dataframe
colnames(germany_combined) <- combined_header

# Load additional physical parameters (wind data and secci depth) ####### 
station_details = read_delim(
 file      = paste0(script_dir, "/", "germany", "/", "Station_details.txt"),
 delim     = "\t",
 col_names = T
)

####### General data transformations ####### 
# remove all columns containing biomass or carbon in the name as we are only interested in cell counts
germany_combined <- germany_combined %>%
 dplyr::select(-starts_with("NA_Carbon_"), -starts_with("NA_Biomass_"))

# Rename columns
germany_combined <- germany_combined %>%
 rename_with(remove_na_prefix, starts_with("NA_"))

# Replace , with .
germany_combined[,] <- sapply(germany_combined[, ], gsub, pattern = ",", replacement = ".")

# Apply transformation to numeric for columns consisting of numbers
germany_combined[, 3:605] <- sapply(germany_combined[, 3:605], as.numeric)

germany_combined <- germany_combined %>% 
 mutate(`Time_(start)` = str_sub(`Time_(start)`, 1, 10)) %>%
 mutate(`Time_(start)` = as.Date(`Time_(start)`, format = "%Y-%m-%d")) %>%
 mutate(
  month    = month(`Time_(start)`),
  year     = year(`Time_(start)`),
  doy      = yday(`Time_(start)`),
  temp     = rowMeans(dplyr::select(., starts_with("Temp")), na.rm = TRUE),
  sal      = rowMeans(dplyr::select(., starts_with("Sal")), na.rm = TRUE)
 ) %>%
 dplyr::rename(
  station  = Name,
  lon      = `Longitude_(start)_[°]`,
  lat      = `Latitude_(start)_[°]`,
  date     = `Time_(start)`
 ) %>%
 dplyr::select(-temp1, -temp2, -temp3, -sal1, -sal2, -sal3) %>%
 pivot_longer(
  cols      = contains("1/m**3"),
  names_to  = "species",
  values_to = "cells_L"
 ) %>%
 mutate(
  species = str_replace_all(
   species,
   c(
    "Alexandrium_pseudogonyaulax.*" = "Alexandrium pseudogonyaulax",
    "Alexandrium_ostenfeldii.*"  = "Alexandrium ostenfeldii"
   )
  ),
  `NO3+NO2` = rowMeans(dplyr::select(., starts_with("NO")), na.rm = TRUE),
  temp      = as.numeric(temp),
  sal       = as.numeric(sal)
 ) 

# Split station names into shorter strings
germany_combined$station <- sapply(strsplit(germany_combined$station, "_"), function(x) x[1])

germany_combined <- germany_combined %>%
  process_harmful_genera_comprehensive()

# Only keep distinct rows as in the IOW dataset each species contains a row even though it was not detected or counted
germany_combined <- germany_combined %>% distinct()

# Keep all present entries of any harmful genus; retain one "No harmful algae" row per date/station
harmful_algae_intermediate <- germany_combined %>%
  filter(if_any(starts_with("probability"), ~ .x == 1L))

germany_combined <- germany_combined %>%
  filter(species %in% c("No harmful algae", NA)) %>%
  group_by(station, date) %>%
  arrange(cells_L) %>%
  slice_head() %>%
  ungroup()

germany_combined <- full_join(germany_combined, harmful_algae_intermediate)

# Include secci depth and wind speed
station_details <- station_details %>%
 dplyr::select(
  station = 'Station_name',
  date    = 'Time_(start)',
  lon     = 'Longitude_(start)_[°]',
  lat     = 'Latitude_(start)_[°]',
  wind_ms = 'Wind_velocity_[m/s]',
  secci   = 'Secchi_depth_[m]'
 )

# extract date from time column and convert to date
station_details$date <-
 str_sub(station_details$date,
     start = 1,
     end   = 10) %>% 
  as.Date(format = "%d.%m.%Y")

# Average all numeric columns of the same station and date
germany_combined <- germany_combined %>%
  ungroup() %>%
  group_by(
    station, date,
    probability,
    probability_Alexandrium, probability_Dinophysis, probability_Pseudonitzschia,
    probability_Azadinium, probability_Chrysochromulina, probability_Prymnesium,
    probability_Amphidinium, probability_Pseudochattonella, probability_Phaeocystis,
    probability_Karlodinium, probability_Cyanobacteria,
    harmful_algae, genus, toxin_syndrome,
    species, `Depth_(start)_[m]`, `Depth_(end)_[m]`
  ) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

# # join algae (germany_combined) and station_details dataframe
germany_combined <-
 left_join(germany_combined, station_details, by = c("station", "date", "lat", "lon"))

germany_combined <- germany_combined %>%
 mutate(
  density = calc_seawater_density(temp, sal),
  DIN     = NH4 + NO3,
  depth   = `Depth_(start)_[m]`,
  cells_L = cells_L / 1000
 ) %>%
 strat_index()

germany_combined <- germany_combined[, ] %>%
  group_by(
    station, date,
    probability,
    probability_Alexandrium, probability_Dinophysis, probability_Pseudonitzschia,
    probability_Azadinium, probability_Chrysochromulina, probability_Prymnesium,
    probability_Amphidinium, probability_Pseudochattonella, probability_Phaeocystis,
    probability_Karlodinium, probability_Cyanobacteria,
    harmful_algae, genus, toxin_syndrome,
    strat, species
  ) %>%
  filter(as.numeric(`Depth_(end)_[m]`) <= 10) %>%
  group_modify(~ dplyr::summarise(.x, across(where(is.numeric), mean, na.rm = TRUE))) %>%
  dplyr::select(
    -c(
      `Longitude_(end)_[°]`,
      `Latitude_(end)_[°]`,
      `Depth_(end)_[m]`,
      dens,
      NO3
    )
  ) %>%
  dplyr::rename(NO3 = `NO3+NO2`, chl = Chl) %>%
  mutate(day = day(date), country = "Germany")

# If a harmful algae entry on a given date/station combination exists, remove "No harmful algae" rows
germany_combined <- germany_combined %>%
  group_by(station, date) %>%
  filter(
    any(species != "No harmful algae") & species != "No harmful algae" |
      all(species == "No harmful algae")
  ) %>%
  ungroup()

# export data
write.table(germany_combined,
      paste0(script_dir, "/", "germany_combined.txt"),
      sep = "\t",
      row.names = FALSE)

####### MERGING OF ALL MONITORING PROGRAMS ####### 
germany_combined = read_delim(
 paste0(script_dir, "/", "germany_combined.txt"),
 delim = "\t",
 col_names = T
)
sweden_combined = read_delim(
 paste0(script_dir, "/", "sweden_combined.txt"),
 delim = "\t",
 col_names = T
)
norway_combined = read_delim(
 paste0(script_dir, "/", "norway_combined.txt"),
 delim = "\t",
 col_names = T
)
denmark_combined = read_delim(
 paste0(script_dir, "/", "denmark_combined.txt"),
 delim = "\t",
 col_names = T
)

# Convert data frames to data tables
setDT(denmark_combined)
setDT(sweden_combined)
setDT(germany_combined)

# Update dates for all data 
data_frames <- list(norway_combined,
                    denmark_combined,
                    sweden_combined,
                    germany_combined)

updated_data_frames <- lapply(data_frames, update_dates)

all_data <- as.data.frame(reduce(updated_data_frames, full_join)) 

# combine stations in close proximity
unique_stations <- all_data %>%
 ungroup() %>%
 dplyr::select(station, lat, lon) %>%
 group_by(station) %>%
 dplyr::summarise(
  lat = mean(lat, na.rm = TRUE),
  lon = mean(lon, na.rm = TRUE),
  .groups = "drop"
 )

all_data <- all_data %>%
 group_by(station) %>%
 mutate(
  lat = mean(lat, na.rm = TRUE),
  lon = mean(lon, na.rm = TRUE)
 ) %>%
 ungroup()

close_stations <- combine_close_stations_dbscan(all_data)

# Update all_data with new stations, lat and lons
all_data <-
 full_join(
  all_data %>% dplyr::select(-lat, -lon),
  close_stations,
  by = c("station" = "old_station")
 )

# overwrite "limiting_conditions"
all_data <- all_data %>% mutate(limiting_conditions = as.factor(
 case_when(
  is.na(DIN) | is.na(PO4) ~ NA,
  DIN < 2 & PO4 < 0.2 ~ "yes",
  DIN > 2 | PO4 > 0.2 ~ "no"
 )
))

# Step 1: Aggregate genus-level probabilities, cell concentrations and categorical fields per combined_station/date
all_data_probs <- all_data %>%
  group_by(combined_station, date) %>%
  summarise(
    # Overall: 1 if any row is present, 0 if all absent, NA if all missing
    probability = case_when(
      any(probability == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability))              ~ NA_integer_,
      TRUE                                 ~ 0L
    ),
    # Individual genus probabilities (1/0/NA)
    probability_Alexandrium = case_when(
      any(probability_Alexandrium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Alexandrium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Dinophysis = case_when(
      any(probability_Dinophysis == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Dinophysis))              ~ NA_integer_, TRUE ~ 0L),
    probability_Pseudonitzschia = case_when(
      any(probability_Pseudonitzschia == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Pseudonitzschia))              ~ NA_integer_, TRUE ~ 0L),
    probability_Azadinium = case_when(
      any(probability_Azadinium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Azadinium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Chrysochromulina = case_when(
      any(probability_Chrysochromulina == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Chrysochromulina))              ~ NA_integer_, TRUE ~ 0L),
    probability_Prymnesium = case_when(
      any(probability_Prymnesium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Prymnesium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Amphidinium = case_when(
      any(probability_Amphidinium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Amphidinium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Pseudochattonella = case_when(
      any(probability_Pseudochattonella == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Pseudochattonella))              ~ NA_integer_, TRUE ~ 0L),
    probability_Phaeocystis = case_when(
      any(probability_Phaeocystis == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Phaeocystis))              ~ NA_integer_, TRUE ~ 0L),
    probability_Karlodinium = case_when(
      any(probability_Karlodinium == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Karlodinium))              ~ NA_integer_, TRUE ~ 0L),
    probability_Cyanobacteria = case_when(
      any(probability_Cyanobacteria == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability_Cyanobacteria))              ~ NA_integer_, TRUE ~ 0L),
    # Genus-level cell concentration totals
    cells_L_Alexandrium      = sum(ifelse(genus == "Alexandrium",        cells_L, 0), na.rm = TRUE),
    cells_L_Dinophysis       = sum(ifelse(genus == "Dinophysis",         cells_L, 0), na.rm = TRUE),
    cells_L_Pseudonitzschia  = sum(ifelse(genus == "Pseudo-nitzschia",   cells_L, 0), na.rm = TRUE),
    cells_L_Azadinium        = sum(ifelse(genus == "Azadinium",          cells_L, 0), na.rm = TRUE),
    cells_L_Chrysochromulina = sum(ifelse(genus == "Chrysochromulina",   cells_L, 0), na.rm = TRUE),
    cells_L_Prymnesium       = sum(ifelse(genus == "Prymnesium",         cells_L, 0), na.rm = TRUE),
    cells_L_Amphidinium      = sum(ifelse(genus == "Amphidinium",        cells_L, 0), na.rm = TRUE),
    cells_L_Pseudochattonella = sum(ifelse(genus == "Pseudochattonella", cells_L, 0), na.rm = TRUE),
    cells_L_Phaeocystis      = sum(ifelse(genus == "Phaeocystis",        cells_L, 0), na.rm = TRUE),
    cells_L_Karlodinium      = sum(ifelse(genus == "Karlodinium",        cells_L, 0), na.rm = TRUE),
    cells_L_Cyanobacteria    = sum(ifelse(genus == "Cyanobacteria",      cells_L, 0), na.rm = TRUE),
    # Categorical / metadata fields
    strat               = combine_strat(strat),
    limiting_conditions = combine_limiting_conditions(limiting_conditions),
    across(c(day, month, year, doy, species, genus, toxin_syndrome, country), first),
    .groups = "drop"
  ) %>%
  mutate(logistic = probability)  # alias kept for downstream compatibility

# Step 2: Average all numeric environmental variables per combined_station/date
all_data_num <- all_data %>%
  group_by(combined_station, date) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::select(-any_of(c("day", "month", "year", "doy", "logistic",
    "cells_L_Alexandrium", "cells_L_Dinophysis", "cells_L_Pseudonitzschia",
    "cells_L_Azadinium", "cells_L_Chrysochromulina", "cells_L_Prymnesium",
    "cells_L_Amphidinium", "cells_L_Pseudochattonella", "cells_L_Phaeocystis",
    "cells_L_Karlodinium", "cells_L_Cyanobacteria")))

# Step 3: Join and compute risk classification levels
all_data <- left_join(all_data_probs, all_data_num, by = c("combined_station", "date")) %>%
  mutate(
    risk_level_Alexandrium      = case_when(cells_L_Alexandrium      >= 10000 ~ "High", cells_L_Alexandrium      >= 1000 ~ "Moderate", cells_L_Alexandrium      > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Dinophysis       = case_when(cells_L_Dinophysis       >= 10000 ~ "High", cells_L_Dinophysis       >= 1000 ~ "Moderate", cells_L_Dinophysis       > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Pseudonitzschia  = case_when(cells_L_Pseudonitzschia  >= 10000 ~ "High", cells_L_Pseudonitzschia  >= 1000 ~ "Moderate", cells_L_Pseudonitzschia  > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Azadinium        = case_when(cells_L_Azadinium        >= 10000 ~ "High", cells_L_Azadinium        >= 1000 ~ "Moderate", cells_L_Azadinium        > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Chrysochromulina = case_when(cells_L_Chrysochromulina >= 10000 ~ "High", cells_L_Chrysochromulina >= 1000 ~ "Moderate", cells_L_Chrysochromulina > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Prymnesium       = case_when(cells_L_Prymnesium       >= 10000 ~ "High", cells_L_Prymnesium       >= 1000 ~ "Moderate", cells_L_Prymnesium       > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Amphidinium      = case_when(cells_L_Amphidinium      >= 10000 ~ "High", cells_L_Amphidinium      >= 1000 ~ "Moderate", cells_L_Amphidinium      > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Pseudochattonella = case_when(cells_L_Pseudochattonella >= 10000 ~ "High", cells_L_Pseudochattonella >= 1000 ~ "Moderate", cells_L_Pseudochattonella > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Phaeocystis      = case_when(cells_L_Phaeocystis      >= 10000 ~ "High", cells_L_Phaeocystis      >= 1000 ~ "Moderate", cells_L_Phaeocystis      > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Karlodinium      = case_when(cells_L_Karlodinium      >= 10000 ~ "High", cells_L_Karlodinium      >= 1000 ~ "Moderate", cells_L_Karlodinium      > 0 ~ "Low", TRUE ~ "None"),
    risk_level_Cyanobacteria    = case_when(cells_L_Cyanobacteria    >= 10000 ~ "High", cells_L_Cyanobacteria    >= 1000 ~ "Moderate", cells_L_Cyanobacteria    > 0 ~ "Low", TRUE ~ "None")
  ) %>%
  update_dates()

# Save all_data as a new txt.file
write.table(all_data,
      file = paste0(script_dir, "/", "all_data.txt"),
      sep = "\t",
      row.names = FALSE)

gc()
rm(list = ls())
.rs.restartR()

####### CELL DENSITY AND STATION PLOTS ####### 
script_dir <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(paste0(script_dir, "/", "Time_series_analysis_custom_functions.R"))
install_packages()

####### read in all_data to start here ####### 
all_data = read_delim(
 file = paste0(script_dir, "/", "all_data.txt"),
 delim = "\t",
 col_names = T
)

# Introduce 5 year and 10 year column
breaks <- c(1997, 2007, 2012, 2017, 2022) 
labels <- c("1997-2006", "2007-2011", "2012-2016", "2017-2021")

all_data <- all_data %>%
 mutate(year2 = cut(
  year,
  breaks      = breaks,
  labels      = labels,
  right       = F
 ))

# Find stations that meet the following criteria
# 1. sampling in the years of 2010:2020 with maximum two years missing;
# 2. sampling in the month 4-10 with maximum of three unique month
# 3. minimum of 20 observations per year) ######
stations_to_keep <- all_data %>%
 filter(year %in% 2010:2020, month %in% 5:9) %>%
 drop_na(probability) %>% 
 group_by(combined_station) %>%
 dplyr::summarise(unique_years = n_distinct(year),
          unique_month = n_distinct(month)) %>% 
 filter(unique_years >= 8 &
      unique_month >= 3) %>% 
 pull(combined_station)

stations_to_keep2 <- all_data %>%
 filter(year %in% 2010:2020) %>%
 drop_na(probability) %>% 
 group_by(combined_station, year) %>%
 dplyr::summarise(unique_doy = n_distinct(doy)) %>% 
 filter(unique_doy >= 10) %>% 
 ungroup() %>%
 pull(combined_station) %>%
 unique()

# Find stations that meet both filtering criteria and remove stations from df that do not fit the criteria
common_stations <- intersect(stations_to_keep, stations_to_keep2)

filtered_data <- all_data %>%
 filter(combined_station %in% common_stations)

# Create nutrient and wind lagged columns 
# Define parameters of interest
parameters_of_interest <-
  c("NH4",
    "NO3",
    "sal",
    "temp",
    "chl",
    "wind_ms",
    "TN",
    "PO4",
    "silicate",
    "N_P")

# convert chlorophyll a from mass/L to mol/L
filtered_data$chl <- filtered_data$chl/893.51

# Introduce nutrient ratios
filtered_data <-
 filtered_data %>% mutate(
  'N_P' = DIN/PO4,
  'Si_N'    = silicate/DIN,
  DON       = TN - DIN,
  'C_chl'   = C/chl
 ) %>%
 mutate('DON_DIN'      = DON/DIN) %>% 
 dplyr::rename(station = combined_station)

# Loop through parameters and lag weeks
n_lags <- 1:5
for (param in parameters_of_interest) {
  for (num_days in n_lags) {
    filtered_data <- create_lagged_columns(filtered_data, param, num_days)
  }
}

# Average latitude and longitude of stations with the same name
station_means <- filtered_data %>%
 group_by(station) %>%
 dplyr::summarise(
  lat = mean(lat, na.rm = TRUE),
  lon = mean(lon, na.rm = TRUE),
  .groups = "drop"
 )

filtered_data <- filtered_data %>%
 dplyr::select(-lat, -lon) %>%
 left_join(station_means, by = "station")

# Change absent/present in probability columns to 0/1 to match the format requirements of the GLMs and GAMs
filtered_data <- filtered_data %>%
  mutate_at(vars(starts_with("probability")),
            ~ case_when(
              is.na(.) ~ NA_real_,
              . == "absent" ~ 0,
              TRUE ~ 1
            ) %>% factor(levels = c(0, 1))
  )

# Group by station and year, then count the number of "present" observations in probability column
counts_overview <- all_data %>%
 filter(probability == "present" &
      combined_station %in% common_stations) %>%
 group_by(combined_station) %>%
 dplyr::summarise(
  present_count = n(),
  lat = as.numeric(sprintf("%.2f", first(lat))),
  lon = as.numeric(sprintf("%.2f", first(lon))),
  .groups = "drop" # Override grouped output
 )

# Group by station and year, then get station characteristics (cells_L, temp, sal...) of present observations
counts_overview2 <- all_data %>%
  filter(probability == "present", combined_station %in% common_stations) %>%
  group_by(combined_station) %>%
  dplyr::summarise(
    max_cells_L = sprintf("%.2e", max(cells_L, na.rm = TRUE)),
    temp_range  = format_range(temp),
    sal_range   = format_range(sal),
    year_range  = paste0(min(year), "-", max(year)),
    month_range = paste0(min(month), "-", max(month)),
    .groups     = "drop"
  ) 

# Group by station and get the total sampling range (years) of ALL observations, including absent observations
counts_overview3 <- all_data %>%
 filter(combined_station %in% common_stations) %>%
 group_by(combined_station) %>%
 dplyr::summarise(
  sampling_range = paste0(min(year), "-", max(year)),
  lat            = as.numeric(sprintf("%.2f", first(lat))),
  lon            = as.numeric(sprintf("%.2f", first(lon)))
 )

# Join all station characteristics dataframes together and arrange by descending latitude
station_characteristics <-
 as.data.frame(reduce(
  list(counts_overview, counts_overview2, counts_overview3),
  full_join
 )) %>% arrange(desc(lat))

# reorder column order
station_characteristics <-
 station_characteristics %>% 
 relocate(combined_station, lat, lon, sampling_range)

# Only keep the first station name of combined stations as otherwise the station names are very long
filtered_data$station <-
 sapply(strsplit(filtered_data$station, ","), function(x) x[1])

# Save filtered_data as a new txt.file
write.table(filtered_data,
      file = paste0(script_dir, "/", "filtered_data.txt"),
      sep = "\t",
      row.names = FALSE)

####### Get dataframe of unique combination of station / lat / lon and introduce station numbering ####### 
unique_stations <-
 filtered_data %>% 
  ungroup() %>%
 dplyr::select(station, lat, lon, probability) %>% 
 unique()

unique_stations <- unique_stations %>%
 arrange(station, desc(probability)) %>%
 group_by(station) %>%
 slice_head(n = 1) %>%
 ungroup() %>%
 arrange(desc(lat))

# Define named vector (mapping)
station_map <- c(
  "Kjempebakken"                   = "NW4",
  "Møkland"                        = "NW3",
  "RA2"                            = "B1",
  "A13"                            = "B2",
  "B7"                             = "B3",
  "Korsfjorden"                    = "NW2",
  "GA1"                            = "B4",
  "C3"                             = "B5",
  "Bjørnafjorden"                  = "NW1",       
  "Indre Oslofjord"                = "S12",     
  "Håøyfjorden"                    = "S13",        
  "Kosterfjorden (NR16)"           = "S11",   
  "SLV Bottnefjorden"              = "S10",     
  "Arendal"                        = "S14",
  "Stretudden"                     = "S9",
  "Havstensfjord"                  = "S8",
  "SLV Saltöfjorden"               = "S7",
  "Å17"                            = "S6",
  "SLÄGGÖ"                         = "S4",
  "SLV Havstensfjorden-Ljungskile" = "S5",
  "Koljöfjord"                     = "S3",
  "SLV Lyresund-Stigfjorden"       = "S2",  
  "Åstol"                          = "S1",
  "DANAFJORD"                      = "D18",
  "BY15 GOTLANDSDJ"                = "B6",
  "N7 OST Nidingen"                = "D17",
  "VIB3708"                        = "L2",
  "N14 Falkenberg"                 = "D15",
  "NOR409"                         = "D16",
  "ANHOLT E"                       = "D13",
  "NOR5503"                        = "D14",
  "VIB3727"                        = "L1",
  "L9  LAHOLMSBUKTEN"              = "D12",
  "9E REF M1V1"                    = "B7",
  "ARH170006"                      = "D10",
  "VSJ20925"                       = "D11",
  "RKB1"                           = "N2",
  "KBH431"                         = "D6",
  "ROS60"                          = "D7",
  "VEJ0006870"                     = "D9",
  "FYN6900017"                     = "D8",
  "RIB1510007"                     = "N1",
  "BRKBMPK2"                       = "B8",
  "FYN6300043"                     = "D5",
  "BY2 ARKONA"                     = "B9",
  "TF0360"                         = "D4",
  "TF0046"                         = "D1",
  "MECKLENBURGER BUCHT"            = "D3",
  "Heiligendamm"                   = "D2"
)

unique_stations <- unique_stations %>%
  mutate(
    station_number = station_map[station],
    label = paste0(station, " (", station_number, ")")
  ) %>%
  drop_na(station_number)

# export station characteristics
tab <-
  station_characteristics %>% mutate(station_number = station_map[combined_station]) %>%
  flextable() %>%
  autofit() %>%
  save_as_docx(path = paste0(script_dir, "/", "station_characteristics.docx"))


####### Plot time series stations for the manuscript ####### 
# Coordinates for two maps
coordinates1 <- data.frame(lon = c(3.5, 3.5, 23, 23),
                           lat = c(54, 73.5, 73.5, 54))

coordinates2 <- data.frame(lon = c(7.5, 7.5, 17, 17), 
                           lat = c(54, 60, 60, 54))

# Create maps using the function
stations_meeting_criteria <- make_station_map(
  coordinates1,
  point_size = 2,
  land_col = "grey75",
  legend_visible = TRUE
)

stations_meeting_criteria2 <- make_station_map(
  coordinates2,
  point_size = 2,
  land_col = "grey75",
  legend_visible = FALSE,
  north_arrow_pad = 0.1, 
  scale_bar_pad = 0.1
)

# Combine both maps side by side and export
P_combined <- ggarrange(stations_meeting_criteria, stations_meeting_criteria2, ncol = 2)

ggsave(
 "stations_meeting_criteria_both_without_label.png",
 P_combined,
 path = paste0(script_dir, "/", "maps"),
 dpi = 300,
 width = 8,
 height = 5,
 units = "in"
)

# For overview plot in supplemental:
# Define sub-basins as ribbons (you can adjust these boundaries)
stations_supplemental <- make_station_map(
  coordinates1,
  point_size = 1,
  land_col = "grey75",
  legend_visible = TRUE
)

ggsave(
  "stations_supplemental.png",
  stations_supplemental,
  path = paste0(script_dir, "/", "maps"),
  dpi = 300,
  width = 12,
  height = 7.5,
  units = "in"
)

# set up coordinates for the first map
coordinates <-
 data.frame(lon = c(8, 8, 15, 15), 
            lat = c(53.5, 60, 60, 53.5))

manuscript_stations <- c(
  "Arendal",
  "ARH170006",
  "VIB3708",
  "NOR409",
  "ANHOLT E",
  "Heiligendamm",
  "Indre Oslofjord",
  "SLV Havstensfjorden-Ljungskile",
  "TF0360",
  "L9  LAHOLMSBUKTEN"
)

# construct map
stations_yearly_probability <-
  basemap(
    data = coordinates,
    bathymetry = F,
    legends = F,
    land.col = "grey75",
    rotate = T
  ) +
  geom_spatial_point(
    data = unique_stations
    %>% filter(station %in% manuscript_stations),
    aes(x = lon, y = lat, col = probability),
    size = 1.5
  ) +
  labs(title = "Yearly analysed stations", x = "Longitude (°E)", y = "Latitude (°N)") +
  scale_color_manual(
    values = c("0" = "#E69F00", "1" = "#009E73"),
    labels = c("0" = "absent", "1" = "present")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    axis.ticks = element_line(),
    strip.background = element_blank(),
    axis.text = element_markdown(size = 12),
    axis.title = element_markdown(size = 12),
    legend.text = element_markdown(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    legend.position = "none",
    legend.justification = c(1, 1)
  ) +
  scale_x_continuous(breaks = c(8, 10, 12, 14),
                     labels = c("8", "10", "12", "14"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = c(54, 56, 58, 60),
                     labels = c("54", "56", "58", "60"),
                     expand = c(0, 0)) +
  ggspatial::annotation_north_arrow(
    location = "tr",
    pad_x = unit(0.05, "in"),
    pad_y = unit(0.1, "in"),
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering
  ) 

# Export stations in yearly probability figure map
ggsave( 
 "stations_yearly_probability.png",
 stations_yearly_probability,
 path = paste0(script_dir, "/", "maps"),
 dpi = 300,
 width = 11,
 height = 7,
 units = "cm"
)

# construct 10 year plot with mean cell densities of each species including the amount of present observations ####### 
mean_cell_dens <- filtered_data %>%
 filter(species != "No Alexandrium" &
      species != "Alexandrium minutum") %>%
 drop_na(cells_L) %>%
 group_by(station, year2, species, lat, lon) %>%
 drop_na(cells_L) %>%
 dplyr::summarise(cells_L = mean(cells_L, na.rm = TRUE), n = n()) %>%
 drop_na(year2)

mean_cell_dens <-
 full_join(mean_cell_dens, unique_stations) %>% 
 mutate(log_cells_L = log(cells_L))


####### Cell density plot preparation ####### 
# Specific breaks and their log-transformed values
breaks <- c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5)
log_breaks <- log(breaks)

custom_labels <- scales::scientific_format()(breaks)

custom_breaks <-
 c(1, 5, 10, 20, 40, 65) 
custom_sizes <-
 c(1, 5, 10, 20, 40, 65) 

coordinates2 <-
 data.frame(lon = c(7.5, 7.5, 20.5, 20.5), 
            lat = c(53.5, 60.5, 60.5, 53.5))

custom_colors <-
 c(viridis(256)[1:200], colorRampPalette(c("orange", "red"))(56))

species_vec <- c(
 "Alexandrium pseudogonyaulax",
 "Alexandrium ostenfeldii",
 "Alexandrium tamarense",
 "Alexandrium spp."
)

plot_list <- list()

####### Cell density plot ####### 
walk(species_vec, function(species_name) {
 
 species_data <- mean_cell_dens %>%
  filter(species == species_name) %>%
  drop_na(n)
 
 plot2 <- basemap(
  data = coordinates2,
  bathymetry = FALSE,
  legends = FALSE,
  land.col = "grey75",
  rotate = TRUE
 ) +
  ggspatial::geom_spatial_point(
   data = species_data,
   aes(x = lon, y = lat, col = log_cells_L, size = n)
  ) +
  facet_grid(. ~ year2) +
  scale_color_gradientn(
   colors = viridisLite::plasma(200),
   limits = log(c(1, 10^5)),
   breaks = log_breaks,
   labels = custom_labels,
   name = stringr::str_wrap("Cell density <br> (cells L<sup>-1</sup>)", width = 2)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  scale_size_area(
   limits = c(0, 65),
   breaks = custom_breaks,
   max_size = 6,
   labels = custom_sizes,
   name = stringr::str_wrap("Number of <br> present <br> observations", width = 2)
  ) +
  theme(
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.background = element_rect(fill = 'white'),
   legend.title = element_markdown(hjust = 0.5, size = 12, margin = margin(b = 10)),
   strip.background = element_blank(),
   axis.text = element_markdown(size = 12),
   axis.title = element_markdown(size = 12),
   plot.title = element_blank(),
   axis.text.x = element_blank(),
   axis.ticks.x = element_blank(),
   axis.title.x = element_blank(),
   strip.text = element_blank(),
   axis.ticks.y = element_blank(),
   axis.text.y = element_blank()
  ) +
  scale_x_continuous(breaks = c(8, 12, 16, 20), labels = c("8", "12", "16", "20"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = c(54, 56, 58, 60), labels = c("54", "56", "58", "60"),
                     expand = c(0, 0))
 
 if (species_name %in% c("Alexandrium pseudogonyaulax", "Alexandrium spp.", "Alexandrium ostenfeldii")) {
  plot2 <- plot2 + theme(legend.position = "none")
 }
 if (species_name == "Alexandrium pseudogonyaulax") {
  plot2 <- plot2 +
   theme(strip.text = element_markdown(size = 12, vjust = 1),
      plot.margin = unit(c(0, 0, -2, 0), "cm"))
 }
 if (species_name %in% c("Alexandrium ostenfeldii", "Alexandrium tamarense")) {
  plot2 <- plot2 + theme(plot.margin = unit(c(-2, 0, -2, 0), "cm"))
 }
 if (species_name == "Alexandrium spp.") {
  plot2 <- plot2 +
   labs(title = "", x = "Longitude (°E)", y = "Latitude (°N)") +
   theme(
    axis.text.x = element_markdown(size = 12),
    axis.ticks.x = element_line(),
    axis.title.x = element_markdown(size = 12),
    axis.ticks.y = element_line(),
    axis.text.y = element_markdown(size = 12),
    axis.title.y = element_markdown(size = 12),
    plot.margin = unit(c(-2, 0, 0, 0), "cm")
   )
 }
 plot_list[[species_name]] <<- plot2
})

# Combine the plots using patchwork
Alex_dens <-
 wrap_plots(plot_list) + plot_layout(axis_titles = "collect", ncol = 1) +
 guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))

ggsave(
 "AP_cell_densities.png",
 Alex_dens,
 path = paste0(script_dir, "/", "maps"),
 dpi = 300,
 width = 10,
 height = 6,
 units = "in"
)

####### Cell density plot 2 ####### 
mean_cell_dens <- filtered_data %>%
 filter(species %in% c("Alexandrium pseudogonyaulax", "No Alexandrium")) %>%
 drop_na(cells_L) %>%
 group_by(station, year, species, lat, lon) %>%
 dplyr::summarise(cells_L = mean(cells_L, na.rm = TRUE), n = n()) 

mean_cell_dens <-
 full_join(mean_cell_dens, unique_stations) %>% 
 mutate(log_cells_L = log(cells_L))

species_data <-
  mean_cell_dens %>% 
  filter(year >= 2006 &
           year <= 2011 ) %>% 
  drop_na(n, probability)

breaks <- c(10 ^ 1, 10 ^ 2, 10 ^ 3, 10 ^ 4, 10 ^ 5)
log_breaks <- log(breaks)
custom_labels <- scales::scientific_format()(breaks)

custom_breaks2 <- c(2, 6, 10) 
custom_sizes2 <-  c(2, 6, 10)  

species_data <- species_data %>%
  group_by(year, station) %>%
  filter(!("Alexandrium pseudogonyaulax" %in% species & species == "No Alexandrium")) %>%
  ungroup()

cell_dens_plot <-
  basemap(
    data = coordinates2,
    bathymetry = F,
    legends = F,
    land.col = "grey75",
    rotate = T
  ) +
  ggspatial::geom_spatial_point(
    data = subset(species_data, species == "Alexandrium pseudogonyaulax"),
    aes(
      x = lon,
      y = lat,
      col = log_cells_L,
      size = n
    )
  ) +
  ggspatial::geom_spatial_point(
    data = subset(species_data, species == "No Alexandrium"),
    aes(
      x = lon,
      y = lat
    ),
    shape = 1, 
    size = 1, 
    color = "black",
    stroke = 0.6        
  ) +
  facet_wrap( ~ year) +
  scale_color_gradientn(
    colors = viridisLite::plasma(200),
    limits = log(c(10, 10^5)),
    breaks = log_breaks,
    labels = custom_labels,
    name = stringr::str_wrap("Cell density <br> (cells L<sup>-1</sup>)", width = 2)
  ) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  theme_minimal() +
  scale_size_area(
    limits = c(0, 10),
    breaks = custom_breaks2,
    max_size = 6,
    labels = custom_sizes2,
    name = stringr::str_wrap("Number of <br> present <br> observations", width = 2)
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_markdown(
      hjust = 0.5,
      size = 10,
      margin = margin(b = 10)
    ),
    strip.background = element_blank(),
    axis.text = element_markdown(size = 12),
    axis.title = element_markdown(size = 12),
    plot.title = element_markdown(
      size = 12,
      hjust = 0,
      vjust = 1,
      colour = "white",
      face = "bold"
    ),
    axis.text.x = element_markdown(size = 12),
    axis.ticks.x = element_line(),
    axis.title.x = element_markdown(size = 12),
    axis.ticks.y = element_line(),
    axis.text.y = element_markdown(size = 12),
    axis.title.y = element_markdown(size = 12),
    strip.text = element_markdown(size = 12),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_x_continuous(breaks = c(8, 12, 16, 20),
                     labels = c("8", "12", "16", "20"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = c(54, 56, 58, 60),
                     labels = c("54", "56", "58", "60"),
                     expand = c(0, 0)) +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))


ggsave(
 "AP_cell_densities2.png",
 cell_dens_plot,
 path = paste0(script_dir, "/", "maps"),
 dpi = 300,
 width = 6,
 height = 4,
 units = "in"
)

####### TIME SERIES ANALYSIS ####### 
script_dir <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(paste0(script_dir, "/", "Time_series_analysis_custom_functions.R"))
install_packages()

# read in filtered_data to start here:
filtered_data = read_delim(
 paste0(script_dir, "/", "filtered_data.txt"),
 delim = "\t",
 col_names = T
)

# change the limiting conditions and stratification key to a factor as it gets read in as character
filtered_data <- filtered_data %>%
 mutate(limiting_conditions = as.factor(limiting_conditions),
     strat = as.factor(strat))

# find stations first that contain both absent and present entries
# which reduces calculation time of the following for-loops
stations_matching_condition <- filtered_data %>%
 group_by(station) %>%
 dplyr::summarise(n_prob = n_distinct(probability)) %>%
 filter(n_prob >= 2) %>%
 pull(station)

# apply check_and_fit function that checks for GLM requirements and
# then calculates the doyly, monthly and yearly probabilities including 95% CI intervals of A. pseudogonyaulax presence
# and fits a logistic regression to the yearly probability patterns
if (!exists("global_models_list", envir = .GlobalEnv)) {
  assign("global_models_list", list(), envir = .GlobalEnv)
}
global_models_list <- list()

for (each_station in stations_matching_condition) {
 check_and_fit(
  filtered_data,
  each_station,
  "probparm_station_monthly",
  "probparm_station_yearly",
  "predicted_probs_station_yearly",
  "results_list_",
  "probparm_station_doyly",
  "probparm_station_doyly_fit"
 )
}

# Filter the working directory to get a list of all non-empty dataframes
non_empty_dfs <-
 lapply(ls(envir = .GlobalEnv, pattern = "^results_list_"), function(list_name) {
  results_list <- get(list_name, envir = .GlobalEnv)
  lapply(results_list, function(df) {
   if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
    return(df)
   }
   return(NULL)
  })
 })

non_empty_dfs <- unlist(non_empty_dfs, recursive = FALSE)
list2env(non_empty_dfs, envir = .GlobalEnv)

parameters_of_interest <-
 c(
  "NH4",
  "NO3",
  "sal",
  "temp",
  "TN",
  "PO4",
  "silicate",
  "wind_ms",
  "DIN",
  "DON",
  "DIP",
  "N_P",
  "DON_DIN",
  "chl",
  "C_chl",
  "strat",
  "limiting_conditions",
  "DIN_DIP"
 )

# Loop through column names of all_data to include lagged parameter columns into parameters_of_interest
matching_columns <- c()
for (col_name in colnames(filtered_data)) {
 # Check if col_name contains any parameter in parameters_of_interest using grepl
 if (any(grepl(paste(parameters_of_interest, collapse = "|"), col_name))) {
  matching_columns <- c(matching_columns, col_name)
 }
}

# update parameters_of_interest and remove NO3_NO2
parameters_of_interest <- matching_columns
parameters_of_interest <- setdiff(parameters_of_interest, "NO3+NO2")

Coef_stations <- data.frame() 

# instead of using years_of_presence, rather use all years since the year of the first observation
# reasoning: after first observation, following years with only absent observations should also be included in the time series analysis!
# find the years for each Alexandrium species
probability_columns <- colnames(filtered_data)[str_detect(colnames(filtered_data), "probability")]
years_since_first_observation <-
 find_years_since_first_observation(filtered_data, probability_columns)

stations <- unique(filtered_data$station)
total_stations <- length(stations)

set.seed(1234)
filtered_data$DIN_DIP <- filtered_data$DIN / filtered_data$DIP

station_parameters_coefficients <-
  lapply(seq_along(stations), function(i) {
    each_station <- stations[i]
    cat("Station being analysed:", each_station, "- Station", i, "of", total_stations, "\n")
    tryCatch({
      station_data <- filtered_data %>% filter(station == each_station)
      years_station <-
        years_since_first_observation %>%
        filter(station == each_station) %>%
        ungroup() %>%
        dplyr::select(year, species)
      process_data(
        station_data,
        years_station,
        c("NO3", "PO4", "temp", "sal", "silicate", "N_P"),
        each_station,
        "probability"
      )
    }, error = function(e) {
      message(sprintf("Error processing station %s: %s", each_station, e$message))
      return(NULL)
    })
  })

# Combine the results into data frames and change the format
Coef_stations <- do.call(rbind, lapply(station_parameters_coefficients, function(station_results) {
  if (is.null(station_results)) return(NULL)
  do.call(rbind, lapply(station_results, function(probability_results) {
    if (is.null(probability_results)) return(NULL)
    do.call(rbind, lapply(probability_results, function(parameter_result) {
      if (is.null(parameter_result)) return(NULL)
      return(parameter_result$coefficients)
    }))
  }))
}))

Bootstrap_distributions <- do.call(rbind, lapply(station_parameters_coefficients, function(station_results) {
  if (is.null(station_results)) return(NULL)
  do.call(rbind, lapply(station_results, function(probability_results) {
    if (is.null(probability_results)) return(NULL)
    do.call(rbind, lapply(probability_results, function(parameter_result) {
      if (is.null(parameter_result)) return(NULL)
      return(parameter_result$bootstrap_distribution)
    }))
  }))
}))

# Introduce confidence intervals for monthly/yearly probability distributions with 1/0 as upper and lower limit, respectively
filtered_dataframes <- grep("^probparm_", ls(), value = TRUE) %>%
 keep( ~ !grepl("doy", .))

updated_dfs <- lapply(mget(filtered_dataframes), calculate_upr_lwr)

for (i in seq_along(filtered_dataframes)) {
 assign(filtered_dataframes[i], updated_dfs[[i]])
}

# Create an empty list to store the results
probparm_station_doyly_fit_characteristics <- list()

# Find doys at which the probability exceeds 0.05 and drops below
probparm_station_doyly_fit_characteristics <-
 probparm_station_doyly_fit %>%
 group_by(station) %>%
 dplyr::reframe(
  "t1"    = min(doy[data >= 0.1 & !is.na(data)], na.rm = TRUE),
  "t2"    = max(doy[data >= 0.1 & !is.na(data)], na.rm = TRUE),
  "p_max" = doy[which.max(data[!is.na(data)])]
 ) %>%
 dplyr::ungroup()

# Introduce water bodies to the filtered stations
water_bodies <- c(
  "estuary",
  "coastal",
  "estuary",
  "open",
  "coastal",
  "estuary",
  "estuary",
  "open",
  "estuary",
  "estuary",
  "coastal",
  "coastal",
  "estuary",
  "coastal",
  "estuary",
  "estuary",
  "estuary",
  "open",
  "estuary",
  "estuary",
  "estuary",
  "estuary",
  "coastal",
  "coastal",
  "open",
  "coastal",
  "estuary",
  "coastal",
  "open",
  "open",
  "estuary",
  "estuary",
  "coastal",
  "coastal",
  "coastal",
  "open",
  "estuary",
  "coastal",
  "estuary",
  "coastal",
  "estuary",
  "coastal",
  "open",
  "open",
  "open",
  "open",
  "open",
  "open",
  "coastal"
)

station_and_waterbody <-
  data.frame(station = unique(filtered_data$station)) %>%
  full_join(unique_stations) %>%
  drop_na(station_number) %>%
  arrange(desc(lat)) %>%
  mutate(water_body = water_bodies)

# Define the base directory path where plots will be saved
script_dir <- getwd() 
base_plot_dir <- file.path(script_dir, "figures", "abiotic_parameters", "deviation")

# Create the base directory if it doesn't exist
if (!dir.exists(base_plot_dir)) {
  dir.create(base_plot_dir, recursive = TRUE)
}

Coef_stations <- left_join(Coef_stations, station_and_waterbody %>% dplyr::select(station, station_number)) %>%
  mutate(
    station_number = factor(
      station_number,
      levels = unique(station_number)[order(
        substr(unique(station_number), 1, 1),
        as.numeric(sub("^[A-Z]+", "", unique(station_number)))
      )]
    )
  )

percent_parameters <- c("NO3", "PO4", "silicate", "N_P")

plot_list <- list()

# Loop through each parameter and create a plot
for (param in unique(Coef_stations$parameter)) {
  param_data <- Coef_stations %>% filter(parameter == param & species == "probability")  %>%
    mutate(sig_color = ifelse(lower_ci > 0 | upper_ci < 0, "darkred", "black"))
  
  y_label_function <- if (param %in% percent_parameters) {
    scales::label_number(suffix = "%", accuracy = 1)
  } else {
    scales::label_number(accuracy = 0.1)
  }
  
  plot <- ggplot(param_data, aes(x = station_number, y = relative_deviation)) +
    geom_point(aes(color = sig_color), size = 1.5) + 
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, color = sig_color), width = 0) +
    scale_color_identity() +  # <-- Use actual color values without legend
    labs(x = "Station", y = "Deviation when <i>A. pseudogonyaulax</i> is present") +
    theme_classic() +
    theme(axis.text.x = element_markdown(angle = 90, size = 9, hjust = 1, vjust = 0.5),
          axis.text.y = element_markdown(size = 9),
          plot.title = element_markdown(
            hjust = 0,
            vjust = 0,
            size = 12,
            face = "bold"
          ),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          legend.title = element_blank(),
          legend.text = element_markdown(size = 12),
          strip.background = element_blank(),
          axis.title.y = element_markdown(size = 12),
          axis.ticks.length.x = unit(0.35, "cm")) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2, minor.ticks = TRUE)) +
    scale_y_continuous(expand = c(0, 0), labels = y_label_function)
  
  plot_list[[param]] <- plot
  
  if (grepl("lag", param, ignore.case = TRUE)) {
    plot_dir <- file.path(base_plot_dir, "lag")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
  } else {
    plot_dir <- base_plot_dir
  }
  
  plot_path <- file.path(plot_dir, paste0(param, ".png"))
  
  ggsave(plot_path, plot, width = 10, height = 8, dpi = 300)
  
}

plots <- list(add_common_theme(plot_list$NO3, label = "a)",  x_label = "Station", y_label = "Deviation when <i>A. pseudogonyaulax</i> is present"), 
              add_common_theme(plot_list$PO4, label = "b)", x_label = "Station", y_label = "Deviation when <i>A. pseudogonyaulax</i> is present"),
              add_common_theme(plot_list$'N_P', label = "c)", x_label = "Station", y_label = "Deviation when <i>A. pseudogonyaulax</i> is present"),
              add_common_theme(plot_list$silicate, label = "d)", x_label = "Station", y_label = "Deviation when <i>A. pseudogonyaulax</i> is present"),
              add_common_theme(plot_list$temp, label =  "e)", x_label = "Station", y_label = "Deviation when <i>A. pseudogonyaulax</i> is present"),
              add_common_theme(plot_list$sal, label = "f)", x_label = "Station", y_label = "Deviation when <i>A. pseudogonyaulax</i> is present"))
              
all_plots <- wrap_plots(plots, ncol = 3, nrow = 2, guides = "collect") +
  plot_layout(axis_titles = "collect")    

ggsave(
  "fig6_manuscript4.png",
  all_plots,
  path = paste0(script_dir, "/", "figures"),
  dpi = 300,
  width = 6,
  height = 5,
  units = "in"
)

# Bootstrap GAM results
n_boot <- 10000
threshold <- 0.1
boot_results <- list()

doy_grid <- seq(0, 365, length.out = 366)

stations_sub <- filtered_data %>%
  group_by(station) %>%
  drop_na(probability) %>%
  dplyr::reframe(
    n_total = n(),
    n_presence = sum(probability == 1, na.rm = TRUE),
    n_unique_doy = n_distinct(doy)
  ) %>% 
  filter(n_presence >= 10)

set.seed(123)  

for (s in unique(stations_sub$station)) {
  dat_station <- filtered_data %>% 
    filter(station == s) %>%
    drop_na(probability)
  
  boot_summary <- replicate(n_boot, {
    dat_boot <- dat_station %>%
      group_by(probability) %>%
      group_modify(~ slice_sample(.x, n = nrow(.x), replace = TRUE)) %>%
      ungroup() %>%
      mutate(
        station = as.factor(as.character(station)),
        doy = as.numeric(doy)
      )
    
    fit <- tryCatch(
      mgcv::gam(probability ~ s(doy, bs = "cp"), 
                data = dat_boot, family = binomial,
                knots = list(doy = c(0, 365))),
      error = function(e) return(NULL)
    )
    
    if (is.null(fit)) return(rep(NA, 3))  
    
    pred <- predict(fit, newdata = data.frame(doy = doy_grid), type = "response")
    
    t1 <- doy_grid[which(pred >= threshold)[1]]
    t2 <- doy_grid[rev(which(pred >= threshold))[1]]
    p_max <- doy_grid[which.max(pred)]
    
    c(t1, t2, p_max)
  }, simplify = "matrix")
  
  boot_results[[s]] <- as.data.frame(t(boot_summary))
  colnames(boot_results[[s]]) <- c("t1", "t2", "p_max")
  boot_results[[s]]$station <- s
}

boot_all <- bind_rows(boot_results)

ci_summary <- boot_all %>%
  pivot_longer(cols = c("t1", "t2", "p_max")) %>%
  group_by(station, name) %>%
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = name, values_from = c(mean, lower, upper))

bootstrap_results <- left_join(boot_all, station_and_waterbody)

bootstrap_results$iteration <- rep(1:10000, times = length(unique(bootstrap_results$station)))

group_diffs <- bootstrap_results %>%
  group_by(iteration, water_body) %>%
 dplyr::summarise(mean_t1 = mean(t1, na.rm = T), .groups = "drop") %>%
  pivot_wider(names_from = water_body, values_from = mean_t1) %>%
  mutate(
    estuary_vs_coastal = estuary - coastal,
    estuary_vs_open = estuary - open,
    coastal_vs_open = coastal - open
  )

group_diffs %>%
  dplyr::summarise(
    mean_diff = mean(estuary_vs_coastal, na.rm = T),
    lower_CI = quantile(estuary_vs_coastal, 0.025, na.rm = T),
    upper_CI = quantile(estuary_vs_coastal, 0.975, na.rm = T)
  )

group_diffs %>%
  dplyr::summarise(
    mean_diff = mean(estuary_vs_open, na.rm = T),
    lower_CI = quantile(estuary_vs_open, 0.025, na.rm = T),
    upper_CI = quantile(estuary_vs_open, 0.975, na.rm = T)
  )

group_diffs %>%
  dplyr::summarise(
    mean_diff = mean(coastal_vs_open, na.rm = T),
    lower_CI = quantile(coastal_vs_open, 0.025, na.rm = T),
    upper_CI = quantile(coastal_vs_open, 0.975, na.rm = T)
  )


coastal_t1 <- bootstrap_results %>% filter(water_body == "coastal") %>% pull(t1)
estuary_t1 <- bootstrap_results %>% filter(water_body == "estuary") %>% pull(t1)
open_t1 <- bootstrap_results %>% filter(water_body == "open") %>% pull(t1)

mean_diff <- mean(estuary_t1, na.rm = TRUE) - mean(coastal_t1, na.rm = TRUE)
pooled_sd <- sqrt((sd(estuary_t1, na.rm = TRUE)^2 + sd(coastal_t1, na.rm = TRUE)^2) / 2)

cohens_d <- mean_diff / pooled_sd

mean_diff <- mean(estuary_t1, na.rm = TRUE) - mean(open_t1, na.rm = TRUE)
pooled_sd <- sqrt((sd(estuary_t1, na.rm = TRUE)^2 + sd(open_t1, na.rm = TRUE)^2) / 2)

cohens_d <- mean_diff / pooled_sd

group_diffs <- bootstrap_results %>%
  group_by(iteration, water_body) %>%
  dplyr::summarise(mean_t1 = mean(t2, na.rm = T), .groups = "drop") %>%
  pivot_wider(names_from = water_body, values_from = mean_t1) %>%
  mutate(
    estuary_vs_coastal = estuary - coastal,
    estuary_vs_open = estuary - open
  )

group_diffs %>%
  dplyr::summarise(
    mean_diff = mean(estuary_vs_coastal, na.rm = T),
    lower_CI = quantile(estuary_vs_coastal, 0.025, na.rm = T),
    upper_CI = quantile(estuary_vs_coastal, 0.975, na.rm = T)
  )

group_diffs %>%
  dplyr::summarise(
    mean_diff = mean(estuary_vs_open, na.rm = T),
    lower_CI = quantile(estuary_vs_open, 0.025, na.rm = T),
    upper_CI = quantile(estuary_vs_open, 0.975, na.rm = T)
  )

coastal_t2 <- bootstrap_results %>% filter(water_body == "coastal") %>% pull(t2)
estuary_t2 <- bootstrap_results %>% filter(water_body == "estuary") %>% pull(t2)
open_t2 <- bootstrap_results %>% filter(water_body == "open") %>% pull(t2)

mean_diff <- mean(estuary_t2, na.rm = TRUE) - mean(coastal_t2, na.rm = TRUE)
pooled_sd <- sqrt((sd(estuary_t2, na.rm = TRUE)^2 + sd(coastal_t2, na.rm = TRUE)^2) / 2)

cohens_d <- mean_diff / pooled_sd
cohens_d

mean_diff <- mean(estuary_t2, na.rm = TRUE) - mean(open_t2, na.rm = TRUE)
pooled_sd <- sqrt((sd(estuary_t2, na.rm = TRUE)^2 + sd(open_t2, na.rm = TRUE)^2) / 2)

cohens_d <- mean_diff / pooled_sd
cohens_d

group_diffs <- bootstrap_results %>%
  group_by(iteration, water_body) %>%
  dplyr::summarise(mean_t1 = mean(p_max, na.rm = T), .groups = "drop") %>%
  pivot_wider(names_from = water_body, values_from = mean_t1) %>%
  mutate(
    estuary_vs_coastal = estuary - coastal,
    estuary_vs_open = estuary - open
  )

group_diffs %>%
  dplyr::summarise(
    mean_diff = mean(estuary_vs_coastal, na.rm = T),
    lower_CI = quantile(estuary_vs_coastal, 0.025, na.rm = T),
    upper_CI = quantile(estuary_vs_coastal, 0.975, na.rm = T)
  )

group_diffs %>%
  dplyr::summarise(
    mean_diff = mean(estuary_vs_open, na.rm = T),
    lower_CI = quantile(estuary_vs_open, 0.025, na.rm = T),
    upper_CI = quantile(estuary_vs_open, 0.975, na.rm = T)
  )

# merge with probparm_station_doyly_fit_characteristics_df containing t1, t2, pmax
probparm_station_doyly_fit_characteristics <-
 probparm_station_doyly_fit_characteristics %>%
 full_join(station_and_waterbody) %>%
 filter(!is.infinite(t1) & !is.na(t1))

# Calculate seasonal means
filtered_data <- filtered_data %>% 
  mutate(across(starts_with("prob"), as.factor)) %>% 
  convert_as_factor(strat, limiting_conditions)

seasonal_means <- calculate_seasonal_mean(filtered_data, "probability")

prob_only <- seasonal_means %>%
 filter(parameter == "probability") %>%
 dplyr::select(-parameter) %>%
 dplyr::rename(probability = data)

seasonal_means <- seasonal_means %>%
 filter(!parameter == "probability") %>% 
 left_join(prob_only, by = c("time", "station"))

facet_parameters <- c("TN",
                      "DIN",
                      "chl",
                      "temp",
                      "limiting_conditions",
                      "strat",
                      "PO4",
                      "sal")

sa_plot <- seasonal_means %>%
 filter(
  parameter %in% facet_parameters) %>%
 mutate(parameter = factor(parameter, levels = facet_parameters)) %>%
 group_by(parameter, station) %>%
 dplyr::summarise(data = mean(data, na.rm = TRUE),
          probability  = mean(probability, na.rm = TRUE),
          .groups      = "drop")

plots <- Map(create_seasonal_means_plot, facet_parameters, letters[1:8])

facet_plots_grid <- wrap_plots(plots, ncol = 4, nrow = 2, guides = "collect") +
 plot_layout(axis_titles = "collect")

ggsave(
 "fig3_manuscript_defense3.png",
 facet_plots_grid,
 path = paste0(script_dir, "/", "figures"),
 dpi = 300,
 width = 9,
 height = 5,
 units = "in"
)

# Statistical analysis and pairwise comparisons of station characteristics, i.e., the days after which
# the probability of presence exceeds or falls below 10% and the DOY with the maximum probability
# List of variables to perform operations on
variables <- c("t1", "t2", "p_max")

perform_aov <- function(var) {
 formula <- as.formula(paste(var, "~ water_body"))

 aov_result <- aov(data = probparm_station_doyly_fit_characteristics, formula)

 print(var)
 print(summary(aov_result))
 print(TukeyHSD(aov_result))
}

aov_results_doyly_characteristics <- lapply(variables, perform_aov)

# add station name and latitude to the dataframes
probparm_station_doyly <-
 full_join(probparm_station_doyly,
      probparm_station_doyly_fit_characteristics %>% 
       dplyr::select(lat, lon, station))

probparm_station_doyly_fit <-
 full_join(probparm_station_doyly_fit,
      probparm_station_doyly_fit_characteristics %>% 
        dplyr::select(lat, lon, station))

tab <-
 probparm_station_doyly_fit_characteristics %>%
 arrange(water_body) %>%
 flextable() %>%
 autofit() %>%
 save_as_docx(path = paste0(script_dir, "/", "probparm_station_doyly_fit_characteristics.docx"))

# Extract the amount of presence observations of A. ostenfeldii and A. pseudogonyaulax and join them together
amount_of_observations_all <- bind_rows(
 filtered_data %>%
  filter(probability == 1) %>%
  mutate(species = "A. pseudogonyaulax"),
 
 filtered_data %>%
  filter(probability_AO == 1) %>%
  mutate(species = "A. ostenfeldii")
) %>%
 group_by(station, year, species) %>%
 dplyr::summarize(n_presence = n(), .groups = "drop") %>%
 left_join(unique_stations %>% dplyr::select(station, station_number), by = "station") %>%
 arrange(desc(station_number)) %>%
 mutate(
  station_number = factor(station_number, levels = unique(station_number)),
  species = factor(species, levels = c("A. pseudogonyaulax", "A. ostenfeldii"))
 )

# Export the dataframe as a table in a word-document
tab <-
  amount_of_observations_all %>%
 flextable() %>%
 autofit() %>%
 save_as_docx(path = paste0(script_dir, "/", "amount_of_observations.docx"))

####### Heatmap of the amount of present observations of A. ostenfeldii and A. pseudogonyaulax ####### 
reordered_levels <- amount_of_observations_all %>%
  dplyr::select(station_number) %>%
  unique() %>%  
  mutate(station_number = as.character(station_number)) %>%
  as_tibble() %>%
  mutate(
    group = str_extract(station_number, "^[A-Z]+"),
    number = as.numeric(str_extract(station_number, "\\d+"))
  ) %>%
  arrange(
    factor(group, levels = c("B", "D", "S", "L", "N", "NW")),
    number
  ) %>%
  pull(station_number)

amount_of_observations_all$station_number <- factor(amount_of_observations_all$station_number,
                                                    levels = rev(reordered_levels))

amount_of_observations_plot <-
  ggplot(amount_of_observations_all %>% filter(year >= 1997),
         aes(x = year, y = station_number, fill = n_presence)) +
  geom_tile() +
  scale_fill_viridis_c(
    limits = c(0, 13),
    breaks = seq(0, 13, by = 3),
    labels = scales::label_number(accuracy = 1),
    name = stringr::str_wrap("Number of <br> present <br> observations", width = 2) 
  ) +
  labs(title = "", x = "Time (year)", y = "Station") +
  theme_classic() +
  facet_wrap(
    ~ species,
    ncol = 1,
    labeller = labeller(species = c("A. ostenfeldii" = "<b>b)</b>",
                                    "A. pseudogonyaulax" = "<b>a)</b>")),
    scales = "free_x"
  ) +
  theme(
    plot.title = element_markdown(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_markdown(
      hjust = 0.5,
      size = 12,
      margin = margin(b = 10)
    ),
    axis.text.y = element_markdown(),
    strip.background = element_blank(),
    axis.text = element_markdown(),
    strip.placement = "outside",
    strip.text = element_markdown(hjust = 0, face = "bold")
  ) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2), expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(1996, 2025, by = 4),
    limits = c(1996, 2025),
    expand = c(0, 0),
    labels = function(x) {
      sprintf("%02d", x %% 100)
    }
  )

ggsave(
 filename = "amount_of_observations.png",
 plot     = amount_of_observations_plot,
 path     = paste0(script_dir, "/", "figures"),
 units    = "in",
 height   = 7,
 width    = 4.5
)

####### PLOT PROBABILITY PATTERNS OVER TIME ####### 
for (each_station in unique(probparm_station_monthly$station)) {
 p.subset <-
  probparm_station_monthly %>% 
  filter(station == each_station)
 filename <- paste0(each_station, "_CI", ".png")
 create_plot(
  p.subset,
  NULL,
  "station",
  "month",
  paste0(each_station, " monthly"),
  paste0(script_dir, "/", "figures", "/", "stations_monthly"),
  filename 
 )
}

# Monthly probability patterns of estuaries, coastal regions and open water stations
# add station_numbers and water_body to the filtered_data
filtered_data <-
 filtered_data %>% 
 left_join(station_and_waterbody %>% 
       dplyr::select(station, station_number, water_body))

result_all <-
 do.call(rbind, lapply(probability_columns, function(prob_col) {
  calculate_seasonal_probability(filtered_data %>% 
                   filter(station %in% unique_stations$station),
                  prob_col)
 }))

result_all <- result_all %>%
  mutate(
    water_body = if_else(water_body == "open", "open waters", water_body),
    water_body = factor(water_body, levels = c("estuary", "coastal", "open waters"))
  ) %>%
  filter(species %in% c("probability", "probability_AO"))

all_plots_monthly <- list()

for (i in unique(result_all$water_body)) {
 P <-
  ggplot(
   result_all %>% filter(water_body == i),
   aes(x = month, y = data, group = species)
  ) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.75) +
  ylab("Probability of presence of <i>Alexandrium</i>") +
  xlab("") +
  ggtitle(i) +
  theme_classic() +
  theme(
   plot.title = element_markdown(
    hjust = 0,
    vjust = 0,
    size = 12,
    face = "bold"
   ),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.background = element_rect(fill = 'white'),
   legend.title = element_blank(),
   axis.text.x = if (i != "open waters") element_blank() else element_markdown(angle = 90, size = 10),
   axis.text.y = element_markdown(size = 10),
   legend.text = element_markdown(size = 12),
   strip.background = element_blank(),
   axis.title.y = element_markdown(size = 12),
   axis.title.x = if (i != "open waters") element_blank() else element_markdown(size = 12)
  ) +
  scale_x_discrete(breaks = 1:12,
           labels = month.abb[1:12],
           limits = factor(1:12),
           expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0.0, 0.5, 0.1),
            labels = seq(0.0, 0.5, 0.1), expand = c(0, 0)) +
  geom_ribbon(aes(
   ymin = lwr,
   ymax = upr,
   fill = species
  ),
  alpha = 0.25) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_fill_manual(
   values = c("#377eb8", "#ff7f00"),
   labels = c("<i>A. pseudogonyaulax</i>", "<i>A. ostenfeldii</i>"),
   guide = guide_legend(
    nrow = 1,
    keyheight = unit(0.2, "cm"),
    keywidth = unit(0.2, "cm")
   )
  )
 all_plots_monthly[[i]] <- P
}

all_plots <- 
  all_plots_monthly[[2]] +
  all_plots_monthly[[1]] +
  all_plots_monthly[[3]] +
  plot_layout(axes = "collect_y", ncol = 1, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box.margin = margin(t = -5)))

ggsave(
 filename = "all_stations_monthly.png",
 plot = all_plots,
 path = paste0(script_dir, "/", "figures"),
 units = "in",
 width = 4,
 height = 6
)

# plot doyly data of each station and export
for (each_station in unique(probparm_station_doyly$station)) {
 p.subset <-
  probparm_station_doyly %>% filter(station == each_station) %>%
  mutate(doy = as.numeric(doy))
 p.subset2 <-
  probparm_station_doyly_fit %>% filter(station == each_station)
 filename <- paste0(each_station, "_CI", ".png")
 create_plot(
  p.subset,
  p.subset2,
  "station",
  "doy",
  paste0(each_station, " doyly"),
  paste0(script_dir, "/", "figures", "/", "station_doyly"),
  paste0(each_station, ".png")
 )
}

# plot yearly data of each station and export
for (each_station in unique(probparm_station_yearly$station)) {
 p.subset <-
  probparm_station_yearly %>% filter(station == each_station) %>% distinct()
 p.subset2 <-
  predicted_probs_station_yearly %>% filter(station == each_station)
 filename <- paste0(each_station, "_CI", ".png")
 create_plot(
  p.subset,
  p.subset2,
  "station",
  "year",
  paste0(each_station),
  paste0(script_dir, "/", "figures", "/", "station_yearly"),
  filename
 )
}

# Modelling the probability of presence as a function of temperature for each station
filtered_data_temp <- filtered_data %>% 
 drop_na(temp, probability)

results_list <- list()

gam_temp <- filtered_data_temp %>% 
filter(month <= 10 & month >= 5 & year >= 2008)

M1a_gam <- mgcv::gam(data = gam_temp %>% filter(temp < 21.5 & temp >= 10),
       probability ~ s(temp, bs = "tp", k = 3),
       family = binomial)
print(summary(M1a_gam))
gam.check(M1a_gam)
temp_gam_results <- data.frame(temp = seq(10, 21.5, length.out = 1000))
 
pred <- predict(M1a_gam, temp_gam_results, type = "response", se.fit = TRUE)
   
temp_gam_results$probability <- pred$fit
temp_gam_results$se <- pred$se.fit
temp_gam_results$lower <- temp_gam_results$probability - 1.96 * temp_gam_results$se
temp_gam_results$upper <- temp_gam_results$probability + 1.96 * temp_gam_results$se

temp_gam_plot <- ggplot(temp_gam_results, aes(x = temp, y = probability)) +
 geom_point(
  data = filtered_data_temp %>% filter(temp >= 10 &
                      temp < 21.5 &
                      year >= 2008 &
                      month <= 10 &
                      month >= 5) %>%
   mutate(
    temp        = round(temp, digits = 0),
    probability = as.numeric(as.character(probability))
   ) %>%
   aggregate(probability ~ temp, FUN = "mean"),
  aes(x = temp, y = probability)
 ) +
 geom_line(linewidth = 1) +
 geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
 labs(title = "a)", x = "Temperature (°C)", y = "Probability of presence<br> of <i>A. pseudogonyaulax</i>") +
 theme_classic() +
 theme(
  plot.title = element_markdown(face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = 'white'),
  legend.title = element_blank(),
  legend.text = element_markdown(),
  axis.title.y = element_markdown(size = 12),
  strip.background = element_blank(),
  axis.text = element_markdown(size = 10),
  strip.text = element_markdown(),
  strip.placement = "outside",
  legend.position = "top"
 ) +
  scale_x_continuous(breaks = seq(6, 32, by = 2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-0.01, 0.35), breaks = seq(0, 0.35, 0.05), expand = c(0.01, 0.01))

ggsave(
 filename = "temp_Ap.png",
 plot     = temp_gam_plot,
 path     = paste0(script_dir, "/", "figures"),
 units    = "in",
 width    = 4,
 height   = 3.5
)

# Modelling the probability of presence as a function of temperature for each station
filtered_data_sal <- filtered_data %>% drop_na(sal, probability)

results_df <- NULL
results_list <- NULL

gam_sal <-
 filtered_data_sal %>% filter(sal <= 32 & sal >= 5 &
                 year >= 2008 & month <= 10 & month >= 5)

M1a_gam <-
 mgcv::gam(data = gam_sal,
      probability ~ s(sal, k = 4, bs = "tp"),
      family = binomial)

print(summary(M1a_gam))
gam.check(M1a_gam)
sal_gam_results <- data.frame(sal = seq(5, 32, length.out = 10000))

pred <- predict(M1a_gam, sal_gam_results, type = "response", se.fit = TRUE)

sal_gam_results$probability <- pred$fit
sal_gam_results$se <- pred$se.fit
sal_gam_results$lower <- sal_gam_results$probability - 1.96 * sal_gam_results$se
sal_gam_results$upper <- sal_gam_results$probability + 1.96 * sal_gam_results$se

sal_gam_plot <- ggplot(sal_gam_results , aes(x = sal, y = probability)) +
 geom_line(linewidth = 1) +
 geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
 geom_point(
  data = filtered_data_sal %>% filter(sal <= 32 &
                     sal >= 5 &
                     year >= 2008 &
                     month <= 10 &
                     month >= 5) %>% mutate(
                      sal = round(sal, digits = 0),
                      probability = as.numeric(as.character(probability))
                     ) %>% aggregate(probability ~ sal, FUN = "mean"),
  aes(x = sal, y = probability)
 ) +
 labs(title = "b)", x = "Salinity", y = "Probability of presence<br> of <i>A. pseudogonyaulax</i>") +
 theme_classic() +
 theme(
  plot.title = element_markdown(face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = 'white'),
  legend.title = element_blank(),
  legend.text = element_markdown(),
  axis.title.y = element_markdown(size = 12),
  strip.background = element_blank(),
  axis.text = element_markdown(size = 10),
  strip.text = element_markdown(),
  strip.placement = "outside",
  legend.position = "top"
 ) +
 scale_x_continuous(breaks = seq(6, 32, by = 2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-0.01, 0.35), breaks = seq(0, 0.35, 0.05), expand = c(0.01, 0.01))

ggsave(
 filename = "sal_Ap.png",
 plot = sal_gam_plot,
 path = paste0(script_dir, "/", "figures"),
 units = "in",
 width = 4,
 height = 3.5
)

P_combined <- temp_gam_plot + sal_gam_plot + plot_layout(ncol = 2, axis_titles = "collect")

ggsave(
 filename = "sal_temp_Ap.png",
 plot = P_combined,
 path = paste0(script_dir, "/", "figures"),
 units = "in",
 width = 8,
 height = 3.5
)

# Create method figure as part of figure 2 in manuscript:
probparm_station_doyly <-
 left_join(probparm_station_doyly, station_and_waterbody)

probparm_station_doyly_fit <-
 left_join(probparm_station_doyly_fit, station_and_waterbody)

plot_list <- list()

for (each_water_body in unique(probparm_station_doyly$water_body)) {
 p.subset <-
  probparm_station_doyly %>% 
  filter(water_body == each_water_body & station != "STRETUDDEN") %>%
  mutate(doy = as.numeric(doy)) %>% 
  drop_na(water_body)
 
 p.subset2 <-
  probparm_station_doyly_fit %>% 
  filter(water_body == each_water_body) %>% 
  drop_na(water_body)
 
 plot <- create_plot(
  p.subset,
  p.subset2,
  "station",
  "doy",
  paste0(each_water_body, " doyly"),
  paste0(script_dir, "/", "figures"),
  paste0(each_water_body, ".png") # Specify the complete filename here
 )

 plot_list <- c(plot_list, plot)
}

gaussian_distribution <- function(x, mean, sd) {
  exp(-(x - mean) ^ 2 / (2 * sd ^ 2)) / (sd * sqrt(2 * pi))
}

x_values <- seq(-5, 5, length.out = 1000)
y_values <- gaussian_distribution(x_values, mean = 0, sd = 1.25)

y_values <- y_values / max(y_values) * 0.5

t1 <- min(x_values[y_values > 0.1])
t2 <- max(x_values[y_values > 0.1])

tmax <- x_values[which.max(y_values)]
pmax <- max(y_values)

df <- data.frame(x = x_values, y = y_values)

plot_segments <- data.frame(
 x = c(t1, t1, t2),
 xend = c(t2, t1, t2),
 y = c(0.1, 0, 0),
 yend = c(0.1, 0.1, 0.1)
)

doy_method <- ggplot(df, aes(x, y)) +
 geom_line() + 
 geom_segment(data = plot_segments,
        aes(
         x = x,
         xend = xend,
         y = y,
         yend = yend
        ),
        linetype = "dashed") +
 geom_richtext(
  x = tmax,
  y = 0.15,
  label = "D = t<sub>2</sub> - t<sub>1</sub>",
  fill = NA,
  label.color = NA,
  size = 4
 ) +
 geom_richtext(
  x = tmax,
  y = 0.55,
  label = "p<sub>max</sub>",
  fill = NA,
  label.color = NA,
  size = 4
 ) +
 scale_x_continuous(
  breaks = c(tmax, t1, t2),
  labels = c("t<sub>max</sub>", "t<sub>1</sub>", "t<sub>2</sub>"), expand = c(0, 0)
 ) +
 scale_y_continuous(breaks = seq(0, 0.6, by = 0.1),
           limits = c(0, 0.6), expand = c(0, 0)) + 
 labs(x = "", y = "Probability of presence of <i>A. pseudogonyaulax</i>") +
 ggtitle("a)") +
 theme_classic() +
 theme(
  plot.title = element_markdown(face = "bold", size = 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = 'white'),
  legend.title = element_blank(),
  axis.text.x = element_markdown(size = 12),
  axis.title.y = element_markdown(size = 12),
  axis.title.x = element_markdown(size = 12),
  strip.text = element_markdown(),
  strip.background = element_blank(),
  plot.margin = margin(
   t = 0.25,
   r = 0,
   b = 0.25,
   l = 0,
   unit = "cm"
  ),
  strip.text.x = element_markdown(hjust = 0, margin = margin(l = 0))
 )

plots <- list(add_common_theme(doy_method, label = "a)", x_label = "", y_label = "Probability of presence of <i>A. pseudogonyaulax</i>"), 
              add_common_theme(plot_list$`estuary doyly`, label = "b)", x_label = "", y_label = "Probability of presence of <i>A. pseudogonyaulax</i>"),
              add_common_theme(plot_list$`coastal doyly`, label = "c)", x_label = "Time (doy)", y_label = "Probability of presence of <i>A. pseudogonyaulax</i>"),
              add_common_theme(plot_list$`open doyly`, label = "d)", x_label = "Time (doy)", y_label = "Probability of presence of <i>A. pseudogonyaulax</i>"))

fig2_manuscript <- wrap_plots(plots, ncol = 2, nrow = 2) +
  plot_layout(axis_titles = "collect") 

ggsave(
 "fig2_manuscript.png",
 fig2_manuscript,
 path = paste0(script_dir, "/", "figures"),
 dpi = 300,
 width = 6,
 height = 6,
 units = "in"
)

manuscript_stations <-
 probparm_station_yearly %>% 
 filter(station %in% c(
   "Arendal",
   "ARH170006",
   "VIB3708",
   "NOR409",
   "ANHOLT E",
   "Heiligendamm",
   "Indre Oslofjord",
   "SLV Havstensfjorden-Ljungskile",
   "TF0360",
   "L9  LAHOLMSBUKTEN"
  ))

manuscript_stations <-
 calculate_upr_lwr(manuscript_stations) %>% 
 left_join(unique_stations, by = "station")

plots <- list()
plots2 <- list()
manuscript_subset <- data.frame()

for (each_station in unique(manuscript_stations$station)) {
 manuscript_subset <-
  manuscript_stations %>% 
  filter(station == each_station)
 
 p.subset2 <-
  predicted_probs_station_yearly %>% 
  filter(station == each_station)
 
 station_number <-
  manuscript_subset$station_number[1] 
 
 gg <- create_plot2(manuscript_subset,
          p.subset2,
          "station",
          "year",
          station_number,
          each_station)
 
 plots[[length(plots) + 1]]   <- ggplotGrob(gg)
 plots2[[each_station]] <- gg
 
}

for(i in 1:length(plots2)){
 name <- unique(plots2[[i]]$data$station)
 ggsave(
 paste0(name, ".png"),
 plots2[[i]],
 path = paste0(script_dir, "/", "figures", "/", "probability"),
 dpi = 300,
 width = 5,
 height = 2.5,
 units = "cm"
) 
}

Fig_manuscript <- plots2$Arendal + plots2$VIB3708 + plots2$NOR409 + plots2$ARH170006 + plot_layout(ncol = 1)
Fig_manuscript2 <- plots2$`SLV Havstensfjorden-Ljungskile` + plots2$`ANHOLT E` + plots2$`L9  LAHOLMSBUKTEN` + plots2$Heiligendamm + plot_layout(ncol = 1)
Fig_manuscript3 <- plots2$`Indre Oslofjord` + plots2$TF0360 + plots2$`Indre Oslofjord` + plots2$`Indre Oslofjord`  + plot_layout(ncol = 1)
Fig_manuscript4 <- plots2$TF0360 + plots2$VIB3708 + plots2$VIB3708 + plots2$VIB3708  + plot_layout(ncol = 1)

ggsave(
  "Fig_prediction1.png",
  Fig_manuscript,
  path = paste0(script_dir, "/", "figures", "/", "probability"),
  dpi = 300,
  width = 2,
  height = 6,
  units = "in"
)

ggsave(
  "Fig_prediction2.png",
  Fig_manuscript2,
  path = paste0(script_dir, "/", "figures", "/", "probability"),
  dpi = 300,
  width = 2,
  height = 6,
  units = "in"
)

ggsave(
  "Fig_prediction3.png",
  Fig_manuscript3,
  path = paste0(script_dir, "/", "figures", "/", "probability"),
  dpi = 300,
  width = 2,
  height = 6,
  units = "in"
)

ggsave(
  "Fig_prediction4.png",
  Fig_manuscript4,
  path = paste0(script_dir, "/", "figures", "/", "probability"),
  dpi = 300,
  width = 2,
  height = 6,
  units = "in"
)

gc()
rm(list = ls())
.rs.restartR()