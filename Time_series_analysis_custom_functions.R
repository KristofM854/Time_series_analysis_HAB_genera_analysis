##########################################
## Time series analysis investigating the potntial expansion of Alexandrium pseudogonyaulax in northern European waters 
## File to store all custom functions of the Time Series Analysis
## by Kristof Möller, Jacob Carstensen, Hans Jakobsen, Annette Engesmo and Bengt Karlson 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 06/24
## Alfred-Wegener-Institute Bremerhaven / International Atomic Energy Agency Monaco
##########################################

install_packages <- function() {
  if (!require("pacman"))
    install.packages("pacman")
  
  pacman::p_load(
    data.table,
    dbscan,
    dplyr,
    flextable,
    fishmethods,
    fuzzyjoin,
    furrr,
    geosphere,
    ggOceanMaps,
    ggplot2,
    ggpubr,
    ggtext,
    ggthemes,
    ggspatial,
    mgcv,
    officer,
    patchwork,
    purrr,
    readr,
    readxl,
    rstatix,
    rstudioapi,
    scales,
    stats,
    stringr,
    tidyverse,
    viridisLite
  )
}

####### Genera configuration: canonical genera vector, display names, and column-name helpers #######

genera_of_interest <- c(
  "Alexandrium", "Dinophysis", "Pseudonitzschia", "Azadinium",
  "Chrysochromulina", "Prymnesium", "Amphidinium", "Pseudochattonella",
  "Phaeocystis", "Karlodinium", "Cyanobacteria"
)

# Italic markdown names used in plot labels
genus_display_names <- c(
  Alexandrium       = "*Alexandrium* spp.",
  Dinophysis        = "*Dinophysis* spp.",
  Pseudonitzschia   = "*Pseudo-nitzschia* spp.",
  Azadinium         = "*Azadinium* spp.",
  Chrysochromulina  = "*Chrysochromulina* spp.",
  Prymnesium        = "*Prymnesium* spp.",
  Amphidinium       = "*Amphidinium* spp.",
  Pseudochattonella = "*Pseudochattonella* spp.",
  Phaeocystis       = "*Phaeocystis* spp.",
  Karlodinium       = "*Karlodinium* spp.",
  Cyanobacteria     = "Cyanobacteria"
)

# Column-name helper functions (single source of truth for naming conventions)
prob_col_name  <- function(genus) paste0("probability_", genus)
cells_col_name <- function(genus) paste0("cells_L_", genus)
risk_col_name  <- function(genus) paste0("risk_level_", genus)

# Mapping from genera_of_interest keys to the actual value stored in the 'genus' column.
# Required because "Pseudo-nitzschia" in the data != "Pseudonitzschia" column-name suffix.
genus_data_names <- c(
  Alexandrium       = "Alexandrium",
  Dinophysis        = "Dinophysis",
  Pseudonitzschia   = "Pseudo-nitzschia",
  Azadinium         = "Azadinium",
  Chrysochromulina  = "Chrysochromulina",
  Prymnesium        = "Prymnesium",
  Amphidinium       = "Amphidinium",
  Pseudochattonella = "Pseudochattonella",
  Phaeocystis       = "Phaeocystis",
  Karlodinium       = "Karlodinium",
  Cyanobacteria     = "Cyanobacteria"
)

# Apply the standard 4-level risk classification (None/Low/Moderate/High) to all
# genera in-place, deriving risk_level_<Genus> from cells_L_<Genus>.
apply_risk_levels <- function(data) {
  for (g in genera_of_interest) {
    cc <- cells_col_name(g)
    rc <- risk_col_name(g)
    data[[rc]] <- dplyr::case_when(
      data[[cc]] >= 10000 ~ "High",
      data[[cc]] >= 1000  ~ "Moderate",
      data[[cc]] > 0      ~ "Low",
      TRUE                ~ "None"
    )
  }
  data
}

####### Helper function for qualifying station selection (used in Phase 1 station filtering) #######
get_qualifying_stations <- function(data, prob_col,
                                    min_years        = 8,
                                    min_months       = 3,
                                    min_doy_per_year = 10,
                                    year_range       = 2010:2020,
                                    month_range      = 5:9) {
  stk1 <- data %>%
    filter(year %in% year_range, month %in% month_range) %>%
    drop_na(!!sym(prob_col)) %>%
    group_by(combined_station) %>%
    dplyr::summarise(uy = n_distinct(year), um = n_distinct(month), .groups = "drop") %>%
    filter(uy >= min_years & um >= min_months) %>%
    pull(combined_station)

  stk2 <- data %>%
    filter(year %in% year_range) %>%
    drop_na(!!sym(prob_col)) %>%
    group_by(combined_station, year) %>%
    dplyr::summarise(ud = n_distinct(doy), .groups = "drop") %>%
    filter(ud >= min_doy_per_year) %>%
    pull(combined_station) %>%
    unique()

  intersect(stk1, stk2)
}

####### Function to calculate seawater density #######
calc_seawater_density <- function(temp, sal) {
  999.842594 +
    6.793952e-2 * temp -
    9.095290e-3 * temp^2 +
    1.001685e-4 * temp^3 -
    1.120083e-6 * temp^4 +
    6.536332e-9 * temp^5 +
    (8.24493e-1 - 4.0899e-3 * temp + 7.6438e-5 * temp^2 -
       8.2467e-7 * temp^3 + 5.3875e-9 * temp^4) * sal +
    (-5.72466e-3 + 1.0227e-4 * temp - 1.6546e-6 * temp^2) * sal^1.5 +
    4.8314e-4 * sal^2
}

####### Function to introduce stratification index if difference between upper and lower two meters ####### 
# of the water colum exceed 1 gm cm^-3
strat_index <- function(df) {
  df %>%
  group_by(station, date) %>%
    mutate(density = as.numeric(density), depth = as.numeric(depth)) %>%
  dplyr::mutate(strat = as.factor(ifelse(
    abs(mean(c(density[depth >= max(depth) - 2], density[max(depth)]), na.rm = TRUE) -
          mean(c(density[depth <= min(depth) + 2], density[min(depth)]), na.rm = TRUE)) >= 1,
    "stratified",
    "not stratified"
  ))) %>%
  ungroup()
}

####### Function for comprehensive genus classification and probability key generation (11 harmful genera) #######
process_harmful_genera_comprehensive <- function(df) {
  df <- df %>%
    mutate(
      # Genus classification: comprehensive pattern matching
      # Pseudochattonella placed before Pseudo-nitzschia to prevent "pseud" prefix false-positives
      genus = case_when(
        str_detect(tolower(species), "alexandrium|alex min|alex ost|alex pse|alex tam|alexandz|alex exc") ~ "Alexandrium",
        str_detect(tolower(species), "dinophysis")                                               ~ "Dinophysis",
        str_detect(tolower(species), "pseudochattonella|pseudochat|pseudo chat")                 ~ "Pseudochattonella",
        str_detect(tolower(species), "pseudo-nitzschia|pseudonitzschia|pseudo nitzschia")        ~ "Pseudo-nitzschia",
        str_detect(tolower(species), "azadinium")                                                ~ "Azadinium",
        str_detect(tolower(species), "chrysochromulina|chryso chromulina")                       ~ "Chrysochromulina",
        str_detect(tolower(species), "prymnesium")                                               ~ "Prymnesium",
        str_detect(tolower(species), "amphidinium")                                              ~ "Amphidinium",
        str_detect(tolower(species), "phaeocystis")                                              ~ "Phaeocystis",
        str_detect(tolower(species), "karlodinium")                                              ~ "Karlodinium",
        str_detect(tolower(species), "nodularia|aphanizomenon|dolichospermum|anabaena|microcystis|planktothrix|woronichinia|cuspidothrix") ~ "Cyanobacteria",
        TRUE ~ "Other"
      ),
      
      # Toxin syndrome classification
      toxin_syndrome = case_when(
        genus == "Alexandrium"                                                                    ~ "PSP",
        genus == "Dinophysis"                                                                     ~ "DSP",
        genus == "Pseudo-nitzschia"                                                               ~ "ASP",
        genus == "Azadinium"                                                                      ~ "AZP",
        genus %in% c("Chrysochromulina", "Prymnesium", "Pseudochattonella", "Karlodinium")        ~ "Fish_killer",
        genus == "Amphidinium"                                                                    ~ "Amphidinol",
        genus == "Phaeocystis"                                                                    ~ "Ecosystem_disruptor",
        genus == "Cyanobacteria"                                                                  ~ "Multiple_toxins",
        TRUE ~ "Other"
      ),
      
      # Harmful algae flag
      harmful_algae = genus %in% c(
        "Alexandrium", "Dinophysis", "Pseudo-nitzschia", "Azadinium",
        "Chrysochromulina", "Prymnesium", "Amphidinium", "Pseudochattonella",
        "Phaeocystis", "Karlodinium", "Cyanobacteria"
      )
    ) %>%
    mutate(
      # Overall harmful algae probability: 1 = present, 0 = absent, NA = missing data
      probability = case_when(
        harmful_algae == TRUE & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L)                        ~ NA_integer_,
        TRUE                                                    ~ 0L
      ),

      # Binary 1/0 per genus (1 = present, 0 = absent, NA = missing data)
      probability_Alexandrium = case_when(
        genus == "Alexandrium"       & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Dinophysis = case_when(
        genus == "Dinophysis"        & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Pseudonitzschia = case_when(
        genus == "Pseudo-nitzschia"  & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Azadinium = case_when(
        genus == "Azadinium"         & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Chrysochromulina = case_when(
        genus == "Chrysochromulina"  & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Prymnesium = case_when(
        genus == "Prymnesium"        & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Amphidinium = case_when(
        genus == "Amphidinium"       & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Pseudochattonella = case_when(
        genus == "Pseudochattonella" & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Phaeocystis = case_when(
        genus == "Phaeocystis"       & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Karlodinium = case_when(
        genus == "Karlodinium"       & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      probability_Cyanobacteria = case_when(
        genus == "Cyanobacteria"     & !is.na(cells_L) & cells_L > 0 ~ 1L,
        is.na(species) | is.na(cells_L) ~ NA_integer_,
        TRUE ~ 0L
      ),
      
      # Rename non-HAB species for downstream filtering
      species = case_when(
        harmful_algae == TRUE ~ species,
        is.na(species)        ~ NA_character_,
        TRUE                  ~ "No harmful algae"
      )
    )
  return(df)
}

####### Function retaining only a single 'No harmful algae' entry per day per station #######
filter_harmful_algae <- function(df, species_col = "species") {
  harmful_species <- df %>%
    filter(!.data[[species_col]] %in% c("No harmful algae")) %>%
    filter(harmful_algae == TRUE)

  no_harmful_species <- df %>%
    filter(.data[[species_col]] == "No harmful algae") %>%
    distinct(date, station, .keep_all = TRUE)

  final_df <- bind_rows(harmful_species, no_harmful_species) %>%
    arrange(date, station)

  return(final_df)
}

####### Function to update day, month, year and doy from date ####### 
update_dates <- function(df) {
  df <- df %>% mutate(date = as.Date(date))
  df %>%
    mutate(
      day   = day(date),
      doy   = yday(date),
      month = month(date),
      year  = year(date))
}

####### Function to remove station suffix from danish datasets ####### 
remove_station_suffix <- function(df) {
  df$station <- sub("-.*", "", df$station)
  df
}

####### Function to parse scientific number format ####### 
scientific_10 <- function(x) {
  parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

####### Function to harmonize metadata of danish datasets ####### 
join_station_metadata <- function(df, station_meta) {
  left_join(df, station_meta, by = "station")
}

####### Function to find the closest stations of wind stations and phytoplankton stations (used for Danish data set) ####### 
find_closest_stations <- function(df1, df2) {
  closest_stations <-
    data.frame(station1 = character(),
               station2 = character(),
               distance = numeric())
  
  for (j in 1:nrow(df2)) {
    station2 <- df2[j,]
    closest_station1 <- NULL
    min_distance <- Inf
    
    for (i in 1:nrow(df1)) {
      station1 <- df1[i,]
      
      distance <-
        distHaversine(station1[c("Lon", "Lat")], station2[c("lon", "lat")])
      
      if (distance < min_distance) {
        closest_station1 <- station1$station
        min_distance <- distance
      }
    }
    
    closest_stations <-
      rbind(
        closest_stations,
        data.frame(
          station1 = closest_station1,
          station2 = station2$station,
          distance = min_distance
        )
      )
  }
  
  return(closest_stations)
}

####### Function to change DMM to DM format ####### 
dmm_to_dd <- function(dmm_coordinates) {
  dmm_parts <-
    strsplit(dmm_coordinates, " ") 
  degrees <-
    as.numeric(dmm_parts[[1]][1])    
  minutes <-
    as.numeric(dmm_parts[[1]][2]) 
  dd_coordinates <-
    degrees + (minutes / 60)
  return(dd_coordinates)
}

####### Function to remove "NA_" prefix ####### 
remove_na_prefix <- function(col_name) {
  gsub("^NA_", "", col_name)
}

####### Function to combine stations in close proximity with a threshold of 1km ####### 
combine_close_stations_dbscan <- function(df, threshold_km = 1) {

  station_coords <- df %>%
    dplyr::select(station, lat, lon) %>%
    dplyr::distinct() %>%
    tidyr::drop_na()
  
  dist_matrix <- geosphere::distm(station_coords[, c("lon", "lat")], fun = geosphere::distHaversine)
  
  db <- dbscan::dbscan(dist_matrix, eps = threshold_km * 1000, minPts = 1) 
  station_coords$cluster_id <- db$cluster
  
  grouped_stations <- station_coords %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(
      combined_station = paste(sort(unique(station)), collapse = ", "),
      combined_lat = mean(lat, na.rm = T),
      combined_lon = mean(lon, na.rm = T),
      .groups = "drop"
    )

    station_mapping <- station_coords %>%
    dplyr::left_join(grouped_stations, by = "cluster_id") %>%
    dplyr::transmute(
      combined_station,
      lat = combined_lat,
      lon = combined_lon,
      old_station = station
    )
  
  return(station_mapping)
}

####### Functions to combine 'strat' and 'limiting_condition' values ####### 
combine_strat <- function(strat_values) {
  if ("stratified" %in% strat_values &
      "not stratified" %in% strat_values) {
    return("stratified/not stratified")
  } else if ("stratified" %in% strat_values) {
    return("stratified")
  } else if ("not stratified" %in% strat_values) {
    return("not stratified")
  } else {
    return(NA)
  }
}

combine_limiting_conditions <-
  function(limiting_conditions_values) {
    if ("yes" %in% limiting_conditions_values &
        "no" %in% limiting_conditions_values) {
      return("yes/no")
    } else if ("yes" %in% limiting_conditions_values) {
      return("yes")
    } else if ("no" %in% limiting_conditions_values) {
      return("no")
    } else {
      return(NA)
    }
  }

####### Function to add lag columns for parameters of interest ####### 
create_lagged_columns <- function(data, column_name, num_days) {
  data <- data[order(data$date), ]
  lagged_values <- rep(NA_real_, nrow(data))  
  
  for (i in 1:nrow(data)) {
    current_date <- data$date[i]
    
    target_date <- current_date - num_days
    match_index <- which(data$date == target_date)
    
    if (length(match_index) == 1) {
      lagged_values[i] <- data[[column_name]][match_index]
    }
  }
  
  col_name <- paste0(column_name, "_lag_", num_days)
  data[[col_name]] <- lagged_values
  
  return(data)
}

####### Log transformation function ####### 
log_transform <- function(x) {
  result <- ifelse(x != 0 & !is.infinite(x), log(x), NA)
  return(result)
}

####### Function to perform probability of presence GLMs #######
# prob_col: name of the binary 0/1 probability column to model (default: "probability")
# stations_mc: vector of stations that have both 0 and 1 for this prob_col (per-genus)
check_and_fit <- function(data,
                          each_group,
                          result_month_name,
                          result_year_name,
                          predicted_probs_name,
                          results_list_name,
                          result_doy_name,
                          result_doy_gam_name,
                          prob_col  = "probability",
                          stations_mc = NULL) {
  if (!exists(results_list_name, envir = .GlobalEnv)) {
    assign(
      results_list_name,
      list(
        result_month    = data.frame(),
        result_year     = data.frame(),
        predicted_probs = data.frame()
      ),
      envir = .GlobalEnv
    )
  }

  # Use a genus-namespaced global_models_list to avoid collisions across genera
  gml_name <- paste0("global_models_list_", gsub("probability_?", "", prob_col))
  if (!exists(gml_name, envir = .GlobalEnv)) {
    assign(gml_name, list(), envir = .GlobalEnv)
  }
  global_models_list_local <- get(gml_name, envir = .GlobalEnv)
  if (!each_group %in% names(global_models_list_local)) {
    global_models_list_local[[each_group]] <- list()
  }

  result_month       <- data.frame()
  result_year        <- data.frame()
  result_doy         <- data.frame()
  result_doy_gam     <- data.frame()
  predicted_probs_df <- data.frame()

  results_list <- get(results_list_name, envir = .GlobalEnv)
  grouping_var <- "station"

  # Determine qualifying stations for this prob_col
  if (is.null(stations_mc)) {
    stations_mc <- unique(data$station)
  }

  data <- data %>%
    filter(!is.na(.data[[prob_col]]) &
             station %in% stations_mc)

  present_years_only <- data %>%
    filter(!!rlang::sym(grouping_var) == each_group) %>%
    mutate(limiting_conditions = as.factor(limiting_conditions),
           strat = as.factor(strat))

  presence_years <-
    present_years_only %>%
    filter(year >= 2007 & .data[[prob_col]] == 1L) %>%
    pull(year) %>%
    unique()

  if (length(presence_years) == 0) return(invisible(NULL))
  first_year_for_model <- min(presence_years)

  present_years_only2 <- present_years_only %>%
    filter(year >= first_year_for_model)

  if (sum(present_years_only[[prob_col]] == 0L, na.rm = TRUE) >= 2 &
      sum(present_years_only[[prob_col]] == 1L, na.rm = TRUE) >= 2 &
      length(unique(present_years_only$year))  > 3 &
      length(unique(present_years_only$month)) > 3) {

    M1 <- glm(
      as.formula(paste0("factor(", prob_col, ") ~ factor(month) + 0")),
      data   = present_years_only %>% drop_na(!!sym(prob_col)),
      family = "binomial"
    )
    global_models_list_local[[each_group]]$probability_month <- M1

    M1a <- glm(
      as.formula(paste0("factor(", prob_col, ") ~ factor(doy) + 0")),
      data   = present_years_only %>% drop_na(!!sym(prob_col)),
      family = "binomial"
    )
    global_models_list_local[[each_group]]$probability_doy <- M1a

    if (nlevels(as.factor(as.character(present_years_only$limiting_conditions))) >= 2) {
      M1_lim <- glm(
        factor(limiting_conditions) ~ factor(month) + 0,
        data   = present_years_only %>% drop_na(limiting_conditions),
        family = "binomial"
      )
      global_models_list_local[[each_group]]$limiting_conditions_month <- M1_lim
    }

    if (nlevels(as.factor(as.character(present_years_only$strat))) >= 2) {
      M1_strat <- glm(
        factor(strat) ~ factor(month) + 0,
        data   = present_years_only %>% drop_na(strat),
        family = "binomial"
      )
      global_models_list_local[[each_group]]$strat_month <- M1_strat
    }

    M1a_gam <- mgcv::gam(
      as.formula(paste0(prob_col, " ~ s(doy, bs = 'cp')")),
      data   = present_years_only %>% drop_na(!!sym(prob_col)),
      family = binomial,
      knots  = list(x = c(0, 365))
    )
    global_models_list_local[[each_group]]$probability_doy_gam <- M1a_gam

    M1a_gam_predict <- predict(M1a_gam, se.fit = TRUE)

    M2 <- glm(
      as.formula(paste0("factor(", prob_col, ") ~ factor(year) + 0")),
      data   = present_years_only %>% drop_na(!!sym(prob_col)),
      family = "binomial"
    )
    global_models_list_local[[each_group]]$probability_year <- M2

    present_years_only$month <- as.character(present_years_only$month)

    CI_M1 <- confint(M1, level = 0.95)[, 1]

    result_month <- data.frame(
      data  = exp(coef(M1)) / (1 + exp(coef(M1))),
      month = unlist(str_extract_all(names(coef(M1)), "\\d+\\.?\\d*|\\.\\d+")),
      CI    = exp(CI_M1) / (1 + exp(CI_M1))
    )

    if (exists("M1_lim")) {
      limiting_conditions_df <- data.frame(
        month              = unlist(str_extract_all(names(coef(M1_lim)), "\\d+\\.?\\d*|\\.\\d+")),
        limiting_conditions = exp(coef(M1_lim)) / (1 + exp(coef(M1_lim)))
      )
    }

    if (exists("M1_strat")) {
      strat_df <- data.frame(
        month = unlist(str_extract_all(names(coef(M1_strat)), "\\d+\\.?\\d*|\\.\\d+")),
        strat = exp(coef(M1_strat)) / (1 + exp(coef(M1_strat)))
      )
    }

    if (exists("limiting_conditions_df")) {
      result_month <- merge(result_month, limiting_conditions_df, by = "month", all = TRUE)
    }

    if (exists("strat_df")) {
      result_month <- merge(result_month, strat_df, by = "month", all = TRUE)
    }

    result_month[[grouping_var]] <- as.factor(each_group)

    CI_M2 <- confint(M2, level = 0.95)[, 1]

    result_year <- data.frame(
      data = exp(coef(M2)) / (1 + exp(coef(M2))),
      year = as.numeric(unlist(str_extract_all(names(coef(M2)), "\\d+\\.?\\d*|\\.\\d+"))),
      CI   = exp(CI_M2) / (1 + exp(CI_M2))
    )

    result_year[[grouping_var]] <- as.factor(each_group)

    result_doy <- bind_rows(
      result_doy,
      data.frame(
        data = exp(coef(M1a)) / (1 + exp(coef(M1a))),
        doy  = unlist(str_extract_all(names(coef(M1a)), "\\d+\\.?\\d*|\\.\\d+"))
      )
    )

    result_doy[[grouping_var]] <- as.factor(each_group)

    result_doy_gam <- data.frame(
      data  = NA,
      doy   = present_years_only$doy,
      CI    = NA,
      p_val = NA
    )

    result_doy_gam$p_val <- summary(M1a_gam)$p.pv[[1]]

    if (summary(M1a_gam)$p.pv < 0.1) {
      result_doy_gam$data <- exp(M1a_gam_predict$fit)  / (1 + exp(M1a_gam_predict$fit))
      result_doy_gam$CI   <- exp(M1a_gam_predict$se.fit) / (1 + exp(M1a_gam_predict$se.fit))
    }

    result_doy_gam[[grouping_var]] <- as.factor(each_group)

    p_vec <- as.integer(as.character(present_years_only2[[prob_col]]))
    model <- glm(
      cbind(p_vec, 1L - p_vec) ~ year,
      data   = present_years_only2 %>% drop_na(!!sym(prob_col)),
      family = binomial
    )

    global_models_list_local[[each_group]]$probability_year_log_reg <- model
    assign(gml_name, global_models_list_local, envir = .GlobalEnv)

    max_year <- max(present_years_only2$year)
    years    <- seq(first_year_for_model, max_year, by = 1)

    predicted_probs <- predict(model, newdata = data.frame(year = years), type = "response")

    predicted_probs_df <- bind_rows(
      predicted_probs_df,
      data.frame(
        year                  = years,
        predicted_probability = predicted_probs,
        p_val = rep(summary(model)$coefficients[1, 4], length.out = length(years))
      )
    )

    predicted_probs_df[[grouping_var]] <- as.factor(each_group)

    results_list[[result_month_name]]  <- bind_rows(results_list[[result_month_name]],  result_month)
    results_list[[result_year_name]]   <- bind_rows(results_list[[result_year_name]],   result_year)
    results_list[[predicted_probs_name]] <- bind_rows(results_list[[predicted_probs_name]], predicted_probs_df)
    results_list[[result_doy_name]]    <- bind_rows(results_list[[result_doy_name]],    result_doy)
    results_list[[result_doy_gam_name]] <- bind_rows(results_list[[result_doy_gam_name]], result_doy_gam)

    assign(results_list_name, results_list, envir = .GlobalEnv)
    return(results_list_name)
  }
}

####### Functions to calculate seasonal means of abiotic parameters #######  
# Helper function to set up the binomial glm of the stratification and limiting conditions index
fit_binomial_glm <- function(df, response, station) {
  if (nlevels(droplevels(df[[response]])) < 2) return(NULL)
  mod <- as.formula(sprintf("factor(%s) ~ factor(month) + 0", response))
  glm_fit <- glm(mod, data = df, family = "binomial")
  coef_exp <- exp(coef(glm_fit))
  data.frame(data = coef_exp / (1 + coef_exp)) %>%
    tibble::rownames_to_column("time") %>%
    filter(str_detect(time, "month")) %>%
    mutate(
      time = readr::parse_number(time),
      station = station,
      parameter = response
    )
}

# Helper function to perform the binomial glm of the stratification and limiting conditions index ####
process_prob_col <- function(station_data, prob_col, station) {
  min_year <- station_data %>%
    filter(.data[[prob_col]] == 1) %>%
    pull(year) %>%
    min(na.rm = TRUE)
  
  min_year <- ifelse(is.finite(min_year) && min_year < 2000, 2000, min_year)
  
  present_years_only <- station_data %>%
    filter(year >= 2008, month >= 5, month <= 10)
  
  absent_entries <- sum(present_years_only[[prob_col]] == 0, na.rm = TRUE)
  present_entries <- sum(present_years_only[[prob_col]] == 1, na.rm = TRUE)
  
  if (absent_entries >= 25 & present_entries >= 6) {
    responses <- c(prob_col, "limiting_conditions", "strat")
    map_dfr(responses, ~fit_binomial_glm(present_years_only, .x, station))
  } else {
    NULL
  }
}

# Helper function to perform the glms for the abiotic parameters ####
process_parameter_model <- function(station_data, param, station) {

    param_data <- station_data %>%
    filter(year >= 2008, month >= 5, month <= 10) %>%
    drop_na(!!sym(param)) %>%
    filter(!is.infinite(!!sym(param))) 

  if (nrow(param_data) <= 25) return(NULL)
  
  mod <- as.formula(sprintf("%s ~ factor(month) + 0", param))
  
  # Parameters that should NOT use log-link
  no_log <- c("temp", "wind_ms", "sal", "chl", "C_chl")
  use_log <- !(param %in% no_log)
  
  # If log-link is to be used, check for negatives and add small constant if needed
  if (use_log) {
    if (any(param_data[[param]] < 0, na.rm = TRUE)) return(NULL)
    min_positive <- min(param_data[[param]][param_data[[param]] > 0], na.rm = TRUE)
    if (is.finite(min_positive)) {
      param_data[[param]] <- param_data[[param]] + min_positive / 10
    }
    
    glm_fit <- tryCatch(
      glm(mod, data = param_data, family = gaussian(link = "log")),
      error = function(e) NULL
    )
  } else {
    glm_fit <- tryCatch(
      glm(mod, data = param_data, family = gaussian),
      error = function(e) NULL
    )
  }
  
  if (is.null(glm_fit)) return(NULL)
  
  coefs <- coef(glm_fit)
  coefs <- if (use_log) exp(coefs) else coefs
  
  tibble::enframe(coefs, name = "month_factor", value = "data") %>%
    filter(stringr::str_detect(month_factor, "month")) %>%
    mutate(
      time = readr::parse_number(month_factor),
      station = station,
      parameter = param
    ) %>%
    dplyr::select(time, station, parameter, data)
}

####### Function to calculate the seasonal means using the three upper helper functions ####### 
calculate_seasonal_mean <- function(data, probability_columns) {
  parameters <- c("NH4", "NO3", "sal", "temp", "TN", "PO4", "chl",
                  "DIP", "DIN", "N_P", "C_chl", "silicate", "Si_N")
  
  data %>%
    split(.$station) %>%
    map_dfr(function(station_data) {
      station <- unique(station_data$station)
      
      # Defensive: skip if empty or missing year column
      if (nrow(station_data) == 0 || !"year" %in% names(station_data)) {
        return(tibble())
      }
      
      prob_results <- map_dfr(probability_columns, ~process_prob_col(station_data, .x, station))
      
      param_results <- map_dfr(parameters, ~process_parameter_model(station_data, .x, station))
      
      bind_rows(prob_results, param_results)
    })
}

####### Custom labeller function for abiotic parameter plots ####### 
custom_labeller <- function(variable) {
  titles <- c(
    "NH4" = "NH<sub>4</sub><sup>+</sup> (<i>&mu;</i>mol L<sup>-1</sup>)",
    "NO3" = "NO<sub>3</sub><sup>-</sup> (<i>&mu;</i>mol L<sup>-1</sup>)",
    "sal" = "Salinity",
    "silicate" = "Silicate (<i>&mu;</i>mol L<sup>-1</sup>)",
    "temp" = "Temperature (°C)",
    "TN" = "TDN (<i>&mu;</i>mol L<sup>-1</sup>)",
    "PO4" = "PO<sub>4</sub><sup>3-</sup> (<i>&mu;</i>mol L<sup>-1</sup>)",
    "limiting_conditions" = "Probability of nutrient<br>limiting conditions",
    "strat" = "Probability of<br>stratification",
    "chl" = "Chl-a (<i>&mu;</i>mol L<sup>-1</sup>)",
    "wind_ms" = "Wind speed (m s<sup>-1</sup>)",
    "C" = "Carbon (&mu;mol L<sup>-1</sup>)",
    "DIP" = "DIP (<i>&mu;</i>mol L<sup>-1</sup>)",
    "DIN" = "DIN (<i>&mu;</i>mol L<sup>-1</sup>)",
    "NO3_PO4" = "DIN:PO<sub>4</sub><sup>3-</sup>",
    "Si_N" = "Si:N",
    "DON" = "DON (<i>&mu;</i>g L<sup>-1</sup>)",
    "DON_DIN" = "DON:DIN",
    "C_chl" = "C:Chl"
  )
  
  return(titles[variable])
}

####### Function to construct seasonal means plot #######
# genus: optional genus name key used for y-axis label (e.g. "Alexandrium").
#        If NULL, uses generic "any harmful algae" label.
create_seasonal_means_plot <- function(facet_name, label, genus = NULL) {

  y_label <- if (!is.null(genus) && genus %in% names(genus_display_names)) {
    paste0("Probability of presence<br> of ", genus_display_names[[genus]])
  } else {
    "Probability of presence<br> of any harmful algae"
  }

  p_subset <- sa_plot %>%
    filter(parameter == facet_name & probability > 0) %>%
    drop_na(probability, data)

  p <- ggplot(p_subset, aes(x = data, y = probability)) +
    geom_point(size = 1.5) +
    facet_wrap(
      ~ parameter,
      scales       = "free_x",
      labeller     = labeller(parameter = custom_labeller),
      strip.position = "bottom"
    ) +
    scale_y_continuous(
      limits = c(0, 0.5),
      breaks = seq(0, 0.5, 0.1),
      labels = seq(0, 0.5, 0.1),
      expand = c(0, 0)
    ) +
    scale_x_continuous(limits = c(min(p_subset$data) * 0.95, max(p_subset$data) * 1.05)) +
    theme_classic() +
    theme(
      plot.title       = element_markdown(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      legend.title     = element_blank(),
      legend.text      = element_markdown(),
      axis.title.y     = element_markdown(size = 12),
      strip.background = element_blank(),
      axis.text        = element_markdown(size = 10),
      strip.text       = element_markdown(),
      strip.placement  = "outside",
      legend.position  = "top"
    ) +
    xlab("") +
    ylab(y_label) +
    ggtitle(paste0(label, ")"))
  return(p)
}

####### Function to find all years with minimum one present observation of respective Alexandrium ####### 
find_years_of_presence <- function(data, probability_columns) {
  years_of_presence <- list()
  
  for (prob_col in probability_columns) {
    years_of_presence_df <- data.frame(station = " ", year = " ")
    
    for (each_station in unique(data$station)) {
      station_years <- data %>%
        filter(station == each_station) %>%
        dplyr::select(year, station, !!rlang::sym(prob_col)) %>%
        filter(!!rlang::sym(prob_col) == 1) %>%
        as.data.frame(station = station, year = year) %>%
        dplyr::select(station, year) %>%
        distinct()
      years_of_presence_df <-
        rbind(years_of_presence_df, station_years)
    }
    
    years_of_presence_df <- years_of_presence_df[-1,]
    
    years_of_presence_df$year <-
      as.numeric(years_of_presence_df$year)
    
    years_of_presence[[prob_col]] <- years_of_presence_df
  }
  
  return(years_of_presence)
}

####### Function to find all years with minimum one present observation of respective Alexandrium ####### 
find_years_since_first_observation <-
  function(data, probability_columns) {
    years_since_first_observation <- list()
    
    for (prob_col in probability_columns) {
      years_df <- data %>%
        filter(!!rlang::sym(prob_col) == 1) %>%
        group_by(station) %>%
        mutate(min_year_present = min(year)) %>%
        filter(year >= min_year_present) %>%
        dplyr::select(station, year, min_year_present) %>%
        distinct() %>%
        mutate(
          year_since_first_observation = year - min_year_present + 1,
          species = prob_col
        )
      
      years_since_first_observation[[prob_col]] <- years_df
    }
    
    years_since_first_observation_df <-
      do.call(rbind, years_since_first_observation)
    years_since_first_observation_df <-
      years_since_first_observation_df %>%
      dplyr::select(station,
                    year,
                    year_since_first_observation,
                    min_year_present,
                    species)
    
    return(years_since_first_observation_df)
  }

####### Function to calculate mean abiotic parameter when respective Alexandrium is present or absent ####### 
# Seasonal variation of the respective abiotic parameter is approximated via a sinus- and cosinusoidal function

sin_function <- function(doy, amplitude) {
  amplitude * sin(2 * pi * doy / 365)
}
cos_function <- function(doy, amplitude) {
  amplitude * cos(2 * pi * doy / 365)
}

# Helper function to fit the model
fit_model <- function(data, parameter) {
  amplitude_init <- 2

  nls_formula <- as.formula(
    sprintf(
      "%s ~ sin_function(doy, amplitude) + cos_function(doy, amplitude)",
      parameter
    )
  )
  
  nls_fit <- nls(
    nls_formula,
    data = data,
    start = c(amplitude = amplitude_init)
  )
  
  return(coef(nls_fit))
}

process_data <- function(data,
                         years_of_presence,
                         parameters_of_interest,
                         each_station,
                         probability_columns,
                         n_bootstraps = 10000) {
  
  parameter_results <- list()
  
  for (prob_col in probability_columns) {
    parameter_results[[prob_col]] <- lapply(parameters_of_interest, function(parameter) {
      
      cat("Parameter being analysed:", parameter, "\n")
      
      param_data <- data %>%
        drop_na(!!rlang::sym(parameter), !!rlang::sym(prob_col))
      
      years_of_presence2 <- years_of_presence %>%
        filter(species == prob_col) %>%
        pull(year)
      
      if (!is.null(years_of_presence2)) {
        param_data <- param_data %>%
          filter(year >= 2008,
                 month >= 5, month <= 10)
        
        # Restrict to months where both presence (1) and absence (0) exist
        valid_months <- param_data %>%
          group_by(month) %>%
          dplyr::summarise(
            has_absent = any(.data[[prob_col]] == 0),
            has_present = any(.data[[prob_col]] == 1)
          ) %>%
          filter(has_absent & has_present) %>%
          pull(month)
        
        param_data <- param_data %>%
          filter(month %in% valid_months)
      }
      
      cat("Months with both presence and absence for", prob_col, ":", valid_months, "\n")
      
      if(length(valid_months) >= 3){
      
      if (!(parameter %in% c(
        parameters_of_interest[grep("temp", parameters_of_interest)],
        "sal",
        "strat",
        "limiting_conditions"
      ))) {
        param_data[[parameter]] <- log_transform(as.numeric(param_data[[parameter]]))
      }
      
      if (sum(param_data[[prob_col]] == 1, na.rm = TRUE) >= 2 &
          nrow(param_data %>% filter(!!rlang::sym(prob_col) == 0)) >= 5 &
          nrow(param_data %>% filter(!!rlang::sym(prob_col) == 1)) >= 5 &
          length(unique(param_data[[parameter]])) >= 2 &
          !(parameter %in% c("strat", "limiting_conditions"))) {
        
        param_data[[prob_col]] <- as.factor(as.character(param_data[[prob_col]]))

        # Pre-compute group row indices once — avoids rebuilding them inside
        # every bootstrap iteration. param_data is already drop_na'd above.
        abs_idx  <- which(param_data[[prob_col]] == "0")
        pres_idx <- which(param_data[[prob_col]] == "1")

        bootstrap_log_devs <- numeric(n_bootstraps)   # differences on log-scale
        bootstrap_absent   <- numeric(n_bootstraps)
        bootstrap_present  <- numeric(n_bootstraps)

        for (i in 1:n_bootstraps) {
          # Base-R stratified sampling: ~10x faster than group_modify(slice_sample)
          boot_idx    <- c(sample(abs_idx,  length(abs_idx),  replace = TRUE),
                           sample(pres_idx, length(pres_idx), replace = TRUE))
          boot_sample <- param_data[boot_idx, , drop = FALSE]

          result <- tryCatch(fit_model(boot_sample, parameter),
                             error = function(e) NULL)
          if (is.null(result)) next
          
          mod <- as.formula(
            sprintf(
              "%s ~ sin_function(doy, %f) + cos_function(doy, %f) + factor(%s) + 0",
              parameter,
              result[["amplitude"]],
              result[["amplitude"]],
              prob_col
            )
          )
          
          M <- tryCatch({
            glm(data = boot_sample, formula = mod, family = gaussian)
          }, error = function(e) {
            message("Error fitting GLM: ", e$message)
            return(NULL)
          })
          
          if (is.null(M)) next
          
          if (parameter %in% c("strat", "limiting_conditions") &
              length(na.omit(unique(boot_sample[[parameter]]))) >= 2) {
            M <- glm(
              data    = boot_sample,
              formula = as.formula(paste0(parameter, " ~ year + month + factor(", prob_col, ") + 0")),
              family  = "binomial"
            )
          }
          
          coef_summary <- summary(M)$coefficients
          if (nrow(coef_summary) >= 4 && !is.na(coef_summary[3, 1]) && !is.na(coef_summary[4, 1])) {
            bootstrap_absent[i] <- coef_summary[3, 1]
            bootstrap_present[i] <- coef_summary[4, 1]
            bootstrap_log_devs[i] <- bootstrap_present[i] - bootstrap_absent[i]
          }
        }
        
        log_transformed <- !(parameter %in% c(
          parameters_of_interest[grep("temp", parameters_of_interest)],
          "sal",
          "strat",
          "limiting_conditions"
        ))
        
        if (log_transformed) {
          ci_log_dev <- quantile(bootstrap_log_devs, probs = c(0.025, 0.975), na.rm = TRUE)
          median_log_dev <- median(bootstrap_log_devs, na.rm = TRUE)
          
          relative_deviation <- (exp(median_log_dev) - 1) * 100
          lower_ci <- (exp(ci_log_dev[1]) - 1) * 100
          upper_ci <- (exp(ci_log_dev[2]) - 1) * 100
          
          median_absent <- exp(median(bootstrap_absent, na.rm = TRUE))
          median_present <- exp(median(bootstrap_present, na.rm = TRUE))
          
          mean_absent <- median_absent
          mean_present <- median_present
        } else {
          diff_samples <- bootstrap_present - bootstrap_absent
          
          median_dev <- median(diff_samples, na.rm = TRUE)
          ci_dev <- quantile(diff_samples, probs = c(0.025, 0.975), na.rm = TRUE)
          
          relative_deviation <- median_dev
          lower_ci <- ci_dev[1]
          upper_ci <- ci_dev[2]
          
          median_absent <- median(bootstrap_absent, na.rm = TRUE)
          median_present <- median(bootstrap_present, na.rm = TRUE)
          
          mean_absent <- median_absent
          mean_present <- median_present
        }
        
        coefficients_data <- data.frame(
          parameter = parameter,
          relative_deviation = relative_deviation,  
          present = mean_present,
          absent = mean_absent,
          lower_ci = lower_ci,
          upper_ci = upper_ci,
          species = prob_col,
          station = each_station
        )
        
        bootstrap_distribution <- data.frame(
          parameter = parameter,
          species = prob_col,
          station = each_station,
          bootstrap_iteration = seq_along(bootstrap_log_devs),
          log_difference = if (log_transformed) bootstrap_log_devs else bootstrap_present - bootstrap_absent,
          present = if (log_transformed) exp(bootstrap_present) else bootstrap_present,
          absent = if (log_transformed) exp(bootstrap_absent) else bootstrap_absent,
          transformed = log_transformed
        )
        
        return(list(
          coefficients = coefficients_data,
          bootstrap_distribution = bootstrap_distribution
        ))
        }} else {
        return(NULL)
      }
    })
  }
  return(parameter_results)
}

####### Function to backtransform log transformed coefficients ####### 
backtransformation_of_log_transformed_coefficients <-
  function(Coef_stations_back,
           parameters_of_interest) {
    Coef_stations_back3 <- data.frame()
      filter_col <- "station"

    for (grouping_value in unique(Coef_stations_back[[filter_col]])) {
      for (each_parameter in unique(Coef_stations_back$parameter)) {
        parameter_subset <- Coef_stations_back %>%
          filter(!!sym(filter_col) == grouping_value,
                 parameter == each_parameter)
        
        unique_species <- unique(parameter_subset$species)

        for (each_species in unique_species) {

            Coef_stations_back2_sub <- parameter_subset %>%
              filter(species == each_species)
            
            if (nrow(Coef_stations_back2_sub) > 0) {

              mean_val <- Coef_stations_back2_sub$mean_deviation
              lower_ci <- Coef_stations_back2_sub$lower_ci
              upper_ci <- Coef_stations_back2_sub$upper_ci
              station <- Coef_stations_back2_sub$station
              absent <- Coef_stations_back2_sub$absent
              present <- Coef_stations_back2_sub$present
              
              if (each_parameter %in% c(
                unique(Coef_stations_back$parameter)[grep("temp", unique(Coef_stations_back$parameter))],
                "strat",
                "limiting_conditions"
              )) {
                new_row <- data.frame(
                  parameter = each_parameter,
                  mean_val = mean_val,
                  lower_ci = lower_ci,
                  upper_ci = upper_ci,
                  mean_absent = absent,
                  mean_present = present,
                  species = each_species,
                  station = station
                )

                Coef_stations_back3 <-
                  bind_rows(Coef_stations_back3, new_row)
              } else if (!(each_parameter %in% c("temp", "strat", "limiting_conditions")) &&
                         !is.na(mean_val)) {

                new_row <- data.frame(
                  parameter = each_parameter,
                  mean_val =  exp(mean_val) / (1 + exp(mean_val)),
                  lower_ci = exp(lower_ci) / (1 + exp(lower_ci)),
                  upper_ci = exp(upper_ci) / (1 + exp(upper_ci)),
                  mean_absent = exp(absent) / (1 + exp(absent)),
                  mean_present = exp(present) / (1 + exp(present)),
                  species = each_species,
                  station = station
                )
                
                Coef_stations_back3 <-
                  bind_rows(Coef_stations_back3, new_row)
              }
            }
          }
        }
      }
    
    return(Coef_stations_back3)
  }

####### Function to calculate upr and lwr limit of the confidence interval of monthly and yearly probabilities ####### 
calculate_upr_lwr <- function(df) {
  if (nrow(df) > 0 && !is.null(df$CI)) {
    df$upr <- pmin(df$data + df$CI, 1)
    df$lwr <- pmax(df$data - df$CI, 0)
  } else {
    df$upr <- NA
    df$lwr <- NA
  }
  return(df)
}

####### Function to calculate seasonal probability of each water body ####### 
calculate_seasonal_probability <- function(df, prob_col) {
  result_all <- data.frame()
  
  for (water_body_level in unique(df$water_body)) {
    present_years_only <- df %>%
      filter(water_body == water_body_level) %>%
      filter(station %in% unique_stations$station) %>%
      filter(year >= (
        df %>%
          filter(!!rlang::sym(prob_col) == 1, water_body == water_body_level) %>%
          pull(year) %>%
          min()
      )) %>%
      drop_na(prob_col)
    
    if (sum(present_years_only[[prob_col]] == 0) >= 2 &
        sum(present_years_only[[prob_col]] == 1) >= 2 &
        length(unique(present_years_only$year)) > 3 &
        length(unique(present_years_only$month)) > 3) {
      
      M1 <-
        glm(formula = as.formula(sprintf("%s ~ factor(month) + 0", prob_col)),
            data = present_years_only,
            family = "binomial")
      
      CI_M1 = confint(M1, level = 0.95)[, 1]
      
      result_month <- data.frame(
        data = exp(coef(M1)) / (1 + exp(coef(M1))),
        month = unlist(str_extract_all(
          names(coef(M1)), "\\d+\\.?\\d*|\\.\\d+"
        )),
        CI = exp(CI_M1) / (1 + exp(CI_M1))
      )
      
      result_month <- calculate_upr_lwr(result_month) %>%
        mutate(station = "all stations",
               species = prob_col,
               water_body = water_body_level)
      
      result_all <- rbind(result_all, result_month)
    }
  }
  
  return(result_all)
}

format_range <- function(x, digits = 1) {
  if (all(is.na(x)))
    return(NA_character_)
  paste0(round(min(x, na.rm = TRUE), digits), "-", round(max(x, na.rm = TRUE), digits))
}

####### Helper function to add common theme and labels ####### 
add_common_theme <- function(plot, label = "", x_label = "time (doy)", y_label = "", ylim_range = NULL) {
  plot <- plot +
    ggtitle(label) +
    theme(plot.title = element_markdown(hjust = 0, size = 12, face = "bold"))
  
  if (!is.null(ylim_range)) {
    plot <- plot + ylim(ylim_range)
  }
  
  plot <- plot + labs(x = x_label, y = y_label) +
    theme(axis.title.y = element_markdown(size = 12))
  
  return(plot)
}


####### Functions to generate maps ####### 
make_station_map <- function(coordinates,
                             point_size = 3,
                             land_col = "grey90",
                             legend_pos = c(0.225, 0.775),
                             legend_visible = TRUE,
                             north_arrow_pad = 0.1,
                             scale_bar_pad = 0.1) {
  basemap(
    data = coordinates,
    bathymetry = FALSE,
    legends = FALSE,
    land.col = land_col,
    rotate = TRUE
  ) +
    geom_spatial_point(data = unique_stations,
                       aes(x = lon, y = lat, col = probability),
                       size = point_size) +
    labs(title = "Yearly analysed stations", x = "Longitude (°)", y = "Latitude (°)") +
    scale_color_manual(
      values = c("0" = "#E69F00", "1" = "#009E73"),
      labels = c("0" = "absent", "1" = "present")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = 'white'),
      legend.title = element_blank(),
      axis.ticks = element_line(),
      axis.text = element_markdown(size = 12),
      axis.title = element_markdown(size = 12),
      legend.text = element_markdown(size = 12, face = "bold", colour = ifelse(legend_visible, "black", "white")),
      legend.position = if (legend_visible) legend_pos else "none"
    ) +
    ggspatial::annotation_north_arrow(
      location = "tl",
      pad_x = unit(north_arrow_pad, "in"),
      pad_y = unit(north_arrow_pad, "in"),
      style = north_arrow_fancy_orienteering
    ) + 
    scale_x_continuous(
      name = "Longitude (°)",
      labels = function(x) paste0(x),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "Latitude (°)",
      labels = function(x) paste0(x),
      expand = c(0, 0)
    ) +
    ggspatial::annotation_scale(
      location = "tr",
      pad_x = unit(scale_bar_pad, "in"),
      pad_y = unit(scale_bar_pad, "in")
    )
}


create_plot_abiotic <- function(facet_name, label) {
  y_label <- custom_labeller(facet_name)
  
  p_subset <- Coef_stations_sub %>%
    filter(parameter == facet_name, species == "probability") %>%
    drop_na(mean)
  ggplot(p_subset,
         aes(
           x = station_number,
           y = mean,
           fill = treat,
           label = significance
         )) +
    scale_fill_colorblind() +
    ggtitle(paste0(label, ")")) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_linerange(
      aes(ymin = mean, ymax = mean + sd),
      position = position_dodge(width = 0.9),
      linewidth = 0.5
    ) +
    ylab(y_label) +
    xlab("Station") +
    theme_classic() +
    theme(
      plot.title = element_markdown(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = 'white'),
      legend.title = element_blank(),
      axis.text.x = element_markdown(angle = 90, size = 10),
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
      strip.text.x = element_markdown(hjust = 0, margin = margin(l = 0))# Adjust overall margins
    )
}

create_plot <-
  function(data,
           model,
           group_var,
           x_var,
           title_var,
           save_path,
           filename,
           genus = NULL) {

    plots <- list()

    y_lab <- if (!is.null(genus) && genus %in% names(genus_display_names)) {
      paste0("Probability of presence of ", genus_display_names[[genus]])
    } else {
      "Probability of presence"
    }

    gg <-
      ggplot(data, aes(
        x = !!sym(x_var),
        y = data,
        group = !!sym(group_var)
      )) +
      ylab(y_lab) +
      xlab(x_var) +
      ggtitle(title_var) +
      theme_classic() +
      theme(
        plot.title = element_markdown(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.title = element_blank(),
        axis.text.x = element_markdown(angle = 90),
        strip.background = element_blank(),
        axis.text = element_markdown()
      )
    
    if (x_var == "month") {
      gg <-
        gg + scale_x_discrete(breaks = 1:12,
                              labels = month.abb[1:12],
                              limits = factor(1:12)) + geom_line(aes(y = upr), col = "red", expand = c(0, 0)) +
        geom_line(aes(y = lwr), col = "red") +
        geom_ribbon(aes(ymin = lwr, ymax = upr),
                    fill = "grey",
                    alpha = 0.25) +
        geom_line() +
        geom_point()
      
    }
    
    if (x_var == "year" && any(grepl("month", filename))) {
      gg <-
        gg + scale_x_discrete(breaks = 1:12,
                              labels = month.abb[1:12],
                              limits = factor(1:12)) + geom_line(aes(y = upr), col = "red", expand = c(0, 0)) +
        geom_line(aes(y = lwr), col = "red") +
        geom_ribbon(aes(ymin = lwr, ymax = upr),
                    fill = "grey",
                    alpha = 0.25) +
        geom_line() +
        geom_point()
      
    }
    
    if (x_var == "year" &&
        any(grepl("year", save_path)) && !is.null(model)) {
      gg <-
        gg + geom_line(data = model,
                       aes(x = year, y = predicted_probability),
                       col = "blue") + geom_line(aes(y = upr), col = "red") +
        geom_line(aes(y = lwr), col = "red") +
        geom_ribbon(aes(ymin = lwr, ymax = upr),
                    fill = "grey",
                    alpha = 0.25) +
        geom_line() +
        geom_point()+
        annotate(
          "text",
          x = 2005,
          y = 0.2,
          label = round(unique(model$p_val), digits = 4),
          fontface = "bold",
          size = 6
        ) 
    }
    
    if (x_var == "doy" & !(grepl("estuary|open|coastal", title_var)))
    {
      gg <-
        gg + geom_line(data = model, aes(x = doy, y = data), col = "blue") +
        geom_point() +
        xlab(paste0("time (", x_var, ")")) +
        scale_x_continuous(limits = c(0, 365), 
                           breaks = seq(0, 365, 30), labels = seq(0, 365, 30), expand = c(0, 0))  +
        scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1), expand = c(0, 0))
    }
    
    if (grepl("estuary|open|coastal", title_var)) {
      gg <-
        gg + geom_line(data = model, aes(x = doy, y = data), col = "blue") +
        scale_x_continuous(limits = c(0, 365), 
                           breaks = seq(0, 365, 30), labels = seq(0, 365, 30), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1), expand = c(0, 0))
      plots[[title_var]] <- gg
    }
    
    if (exists("manuscript_subset")) {
      if (identical(data, manuscript_subset)) {
        gg <- gg +
          scale_x_continuous(breaks = seq(1996, 2020, by = 2),
                             limits = c(1996, 2020), expand = c(0, 0)) +
          scale_y_continuous(limits = c(0, 0.7),
                             breaks = seq(0, 0.7, by = 0.1), expand = c(0, 0)) +
          theme(plot.title = element_blank()) +
          geom_text(aes(x = 2005, y = 0.6, label = title_var), size = 5) +
          geom_point()
        
      }
    }
    
    ggsave(filename = filename,
           plot = gg,
           path = save_path)
    
    return(plots)
  }

create_plot2 <-
  function(data,
           model,
           group_var,
           x_var,
           station_number,
           each_station,
           genus = NULL) {
    gg <-
      ggplot(data, aes(
        x = !!sym(x_var),
        y = data,
        group = !!sym(group_var)
      )) +
      geom_point(size = 1)  +
      geom_line(
        data = model,
        aes(x = year, y = predicted_probability),
        col = "blue",
        linewidth = 0.75,
        linetype = if (unique(model$p_val) <= 0.05) "solid" else "dashed"
      ) +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 10),
        axis.text.x = element_markdown(size = 8),
        axis.text.y.left = element_markdown(size = 8, margin = margin(r = 10)),
        panel.background = element_rect(fill = 'white'),
        legend.position = "none",
        plot.title.position = "plot"
      ) +
      labs(x = NULL) +
      scale_y_continuous(limits = c(0, 0.65),
                         breaks = seq(0, 0.5, by = 0.25),
                         expand = c(0.01, 0.01))  +
      geom_line(aes(y = upr), col = "red", linewidth = 0.75) +
      geom_line(aes(y = lwr), col = "red", linewidth = 0.75) +
      geom_ribbon(aes(ymin = lwr, ymax = upr),
                  fill = "grey",
                  alpha = 0.25) +
      geom_line() +
      scale_x_continuous(
        name = "Time (Year)",
        breaks = seq(2007, 2023, by = 4),
        labels = sprintf("%02d", seq(7, 24, by = 4)),
        limits = c(2007, 2024),
        expand = c(0.01, 0.01)
      ) +
      xlab(if (each_station %in% c("Heiligendamm","TF0360", "NOR409", "L9  LAHOLMSBUKTEN", "ARH170006"))
        x_var
        else
          NULL)
    
    if (!(each_station %in% c("Heiligendamm", "TF0360", "ARH170006"))) {
      gg <- gg +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
    }
    
    if ((each_station %in% c("Heiligendamm", "VEJ0006870"))) {
      gg <- gg +
        theme(axis.text.x = element_markdown(size = 8)) +
        scale_y_continuous(
          limits = c(0, 0.65),
          breaks = seq(0, 0.5, by = 0.25),
          position = "right",
          expand = c(0.01, 0.01)
        )
    }
    
    if ((
      each_station %in% c(
        "ANHOLT E",
        "Pricken",
        "Heiligendamm",
        "L9  LAHOLMSBUKTEN",
        "SLV Lyresund-Stigfjorden",
        "DMU444",
        "SLV Havstensfjorden-Ljungskile"
      )
    )) {
      gg <- gg +
        scale_y_continuous(
          limits = c(0, 0.65),
          breaks = seq(0, 0.5, by = 0.25),
          position = "right",
          expand = c(0.01, 0.01)
        ) +
        theme(axis.text.y.right = element_markdown(size = 8, margin = margin(l = 10)))
    }
    if (each_station == "Indre Oslofjord") {
      gg <- gg +
        scale_y_continuous(limits = c(0, 0.8),
                           breaks = seq(0, 0.75, by = 0.25),
                           expand = c(0.01, 0.01))
    }
    
    return(gg)  
  }
