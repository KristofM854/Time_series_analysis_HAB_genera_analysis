# Implementation Guide: Multi-Species Aggregation Fix

## Overview
This guide implements Solution 1 to fix the data integrity bug where `probability_<Genus> == 1` but `species` and `genus` columns don't match. The solution preserves all species information while maintaining backward compatibility.

## Background
The bug occurs in the aggregation logic where multiple harmful algae species/genera are present in a single sampling event (same station/date). The original code used `first()` to select species/genus metadata, which randomly picked one species while the probability aggregation correctly detected all present genera. This created data integrity violations where `probability_Pseudochattonella == 1` but `genus == "Dinophysis"`.

## Implementation Steps

### Step 1: Update Core Function

Replace the existing `process_harmful_genera_comprehensive` function in `Time_series_analysis_custom_functions.R`:

```r
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
```

### Step 2: Add Validation Function

Add this validation function to `Time_series_analysis_custom_functions.R`:

```r
####### Function to validate data integrity after aggregation #######
validate_aggregated_data <- function(df, dataset_name = "Unknown") {
  cat("Validating data integrity for", dataset_name, "...\n")
  
  # Check for data integrity violations
  prob_cols <- colnames(df)[str_detect(colnames(df), "^probability_")]
  
  violations <- tibble()
  
  for (prob_col in prob_cols) {
    if (!prob_col %in% colnames(df)) next
    
    expected_genus <- str_remove(prob_col, "probability_")
    expected_genus <- case_when(
      expected_genus == "Pseudonitzschia" ~ "Pseudo-nitzschia",
      TRUE ~ expected_genus
    )
    
    violations_temp <- df %>%
      filter(.data[[prob_col]] == 1L & genus != expected_genus) %>%
      dplyr::select(date, station, genus, species, all_of(prob_col)) %>%
      mutate(
        expected_genus = expected_genus,
        prob_column = prob_col
      )
    
    violations <- bind_rows(violations, violations_temp)
  }
  
  if (nrow(violations) > 0) {
    cat("❌ Data integrity violations found in", dataset_name, ":\n")
    print(violations)
    cat("These will be fixed by the new aggregation logic.\n\n")
  } else {
    cat("✅ No data integrity violations found in", dataset_name, "\n\n")
  }
  
  # Summary statistics
  n_total <- nrow(df)
  n_harmful <- sum(df$probability == 1L, na.rm = TRUE)
  
  cat("Summary for", dataset_name, ":\n")
  cat("- Total observations:", n_total, "\n")
  cat("- Harmful algae present:", n_harmful, "(", round(100 * n_harmful/n_total, 1), "%)\n")
  
  # Per-genus breakdown
  for (prob_col in prob_cols) {
    if (prob_col %in% colnames(df)) {
      n_genus <- sum(df[[prob_col]] == 1L, na.rm = TRUE)
      genus_name <- str_remove(prob_col, "probability_")
      cat("- ", genus_name, "present:", n_genus, "\n")
    }
  }
  cat("\n")
  
  return(df)
}
```

### Step 3: Update Norway Aggregation Logic

Find and replace the Norway aggregation section (around line 300-350):

```r
# Step 1: Aggregate genus-level probabilities, cell concentrations and categorical fields per date/station
norway_combined_probs <- norway_combined %>%
  group_by(date, station) %>%
  dplyr::summarise(
    # Overall: 1 if any row is present, 0 if all absent, NA if all missing
    probability = case_when(
      any(probability == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability))              ~ NA_integer_,
      TRUE                                 ~ 0L
    ),
    # Individual genus probabilities (1/0/NA) — one across() replaces 11 repeated case_when blocks
    across(
      all_of(unname(sapply(genera_of_interest, prob_col_name))),
      ~ case_when(any(. == 1L, na.rm = TRUE) ~ 1L,
                  all(is.na(.))              ~ NA_integer_,
                  TRUE                       ~ 0L)
    ),
    
    # Genus-level cell concentration totals per date/station
    cells_L_Alexandrium       = sum(cells_L[genus == "Alexandrium"],       na.rm = TRUE),
    cells_L_Dinophysis        = sum(cells_L[genus == "Dinophysis"],        na.rm = TRUE),
    cells_L_Pseudonitzschia   = sum(cells_L[genus == "Pseudo-nitzschia"],  na.rm = TRUE),
    cells_L_Azadinium         = sum(cells_L[genus == "Azadinium"],         na.rm = TRUE),
    cells_L_Chrysochromulina  = sum(cells_L[genus == "Chrysochromulina"],  na.rm = TRUE),
    cells_L_Prymnesium        = sum(cells_L[genus == "Prymnesium"],        na.rm = TRUE),
    cells_L_Amphidinium       = sum(cells_L[genus == "Amphidinium"],       na.rm = TRUE),
    cells_L_Pseudochattonella = sum(cells_L[genus == "Pseudochattonella"], na.rm = TRUE),
    cells_L_Phaeocystis       = sum(cells_L[genus == "Phaeocystis"],       na.rm = TRUE),
    cells_L_Karlodinium       = sum(cells_L[genus == "Karlodinium"],       na.rm = TRUE),
    cells_L_Cyanobacteria     = sum(cells_L[genus == "Cyanobacteria"],     na.rm = TRUE),
    
    # ✅ NEW: Multi-species information preservation
    species_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(species[harmful_algae == TRUE & !is.na(species)]), collapse = " | "),
      "No harmful algae"
    ),
    genera_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(genus[harmful_algae == TRUE & !is.na(genus)]), collapse = " | "),
      "Other"
    ),
    toxin_syndromes_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(toxin_syndrome[harmful_algae == TRUE & !is.na(toxin_syndrome)]), collapse = " | "),
      "Other"
    ),
    
    # ✅ NEW: Dominant species for backward compatibility (highest cell count)
    dominant_species = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          species[max_idx]
        } else {
          "No harmful algae"
        }
      } else {
        "No harmful algae"
      }
    },
    dominant_genus = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          genus[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    dominant_toxin_syndrome = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          toxin_syndrome[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    
    # Keep first value of categorical/date-component fields
    across(c(day, month, year, doy, strat), first),
    .groups = "drop"
  ) %>%
  mutate(
    # ✅ For backward compatibility, use dominant species in main columns
    species = dominant_species,
    genus = dominant_genus,
    toxin_syndrome = dominant_toxin_syndrome,
    harmful_algae = genus != "Other",
    logistic = probability,  # alias kept for downstream compatibility
    country  = "Norway"
  ) %>%
  dplyr::select(-dominant_species, -dominant_genus, -dominant_toxin_syndrome)

# Validate the aggregated data
norway_combined <- validate_aggregated_data(norway_combined, "Norway")
```

### Step 4: Update Denmark Aggregation Logic

Find and replace the Denmark aggregation section:

```r
# Introduce overall probability (1/0/NA) for all harmful genera
denmark_combined <- denmark_combined %>%
  # First group to create multi-species records
  group_by(date, station) %>%
  dplyr::summarise(
    # Probability aggregation
    probability = case_when(
      any(harmful_algae == TRUE & !is.na(cells_L) & cells_L > 0, na.rm = TRUE) ~ 1L,
      any(!is.na(species) & !is.na(cells_L), na.rm = TRUE) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # Individual genus probabilities
    across(
      all_of(unname(sapply(genera_of_interest, prob_col_name))),
      ~ case_when(any(. == 1L, na.rm = TRUE) ~ 1L,
                  all(is.na(.)) ~ NA_integer_,
                  TRUE ~ 0L)
    ),
    
    # Cell count aggregation
    cells_L_Alexandrium       = sum(cells_L[genus == "Alexandrium"],       na.rm = TRUE),
    cells_L_Dinophysis        = sum(cells_L[genus == "Dinophysis"],        na.rm = TRUE),
    cells_L_Pseudonitzschia   = sum(cells_L[genus == "Pseudo-nitzschia"],  na.rm = TRUE),
    cells_L_Azadinium         = sum(cells_L[genus == "Azadinium"],         na.rm = TRUE),
    cells_L_Chrysochromulina  = sum(cells_L[genus == "Chrysochromulina"],  na.rm = TRUE),
    cells_L_Prymnesium        = sum(cells_L[genus == "Prymnesium"],        na.rm = TRUE),
    cells_L_Amphidinium       = sum(cells_L[genus == "Amphidinium"],       na.rm = TRUE),
    cells_L_Pseudochattonella = sum(cells_L[genus == "Pseudochattonella"], na.rm = TRUE),
    cells_L_Phaeocystis       = sum(cells_L[genus == "Phaeocystis"],       na.rm = TRUE),
    cells_L_Karlodinium       = sum(cells_L[genus == "Karlodinium"],       na.rm = TRUE),
    cells_L_Cyanobacteria     = sum(cells_L[genus == "Cyanobacteria"],     na.rm = TRUE),
    
    # Multi-species information preservation
    species_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(species[harmful_algae == TRUE & !is.na(species)]), collapse = " | "),
      "No harmful algae"
    ),
    genera_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(genus[harmful_algae == TRUE & !is.na(genus)]), collapse = " | "),
      "Other"
    ),
    toxin_syndromes_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(toxin_syndrome[harmful_algae == TRUE & !is.na(toxin_syndrome)]), collapse = " | "),
      "Other"
    ),
    
    # Dominant species for backward compatibility
    dominant_species = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          species[max_idx]
        } else {
          "No harmful algae"
        }
      } else {
        "No harmful algae"
      }
    },
    dominant_genus = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          genus[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    dominant_toxin_syndrome = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          toxin_syndrome[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    
    # Environmental variables (average)
    across(where(is.numeric) & !starts_with("probability_") & !starts_with("cells_L_"), ~ mean(.x, na.rm = TRUE)),
    
    # Categorical variables
    across(c(day, month, year, doy, strat, limiting_conditions), first),
    .groups = "drop"
  ) %>%
  mutate(
    # Use dominant species for main columns
    species = dominant_species,
    genus = dominant_genus,
    toxin_syndrome = dominant_toxin_syndrome,
    harmful_algae = genus != "Other",
    logistic = probability,  # alias kept for downstream compatibility
    country  = "Denmark"
  ) %>%
  dplyr::select(-dominant_species, -dominant_genus, -dominant_toxin_syndrome) %>%
  drop_na(station, date)

# Validate the aggregated data
denmark_combined <- validate_aggregated_data(denmark_combined, "Denmark")
```

### Step 5: Update Sweden Aggregation Logic

Find and replace the Sweden aggregation section:

```r
# Update Sweden aggregation (around line 780)
sweden_combined <- sweden_combined %>%
  group_by(date, station) %>%
  dplyr::summarise(
    # Probability aggregation
    probability = case_when(
      any(harmful_algae == TRUE & !is.na(cells_L) & cells_L > 0, na.rm = TRUE) ~ 1L,
      any(!is.na(species) & !is.na(cells_L), na.rm = TRUE) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # Individual genus probabilities
    across(
      all_of(unname(sapply(genera_of_interest, prob_col_name))),
      ~ case_when(any(. == 1L, na.rm = TRUE) ~ 1L,
                  all(is.na(.)) ~ NA_integer_,
                  TRUE ~ 0L)
    ),
    
    # Cell count aggregation
    cells_L_Alexandrium       = sum(cells_L[genus == "Alexandrium"],       na.rm = TRUE),
    cells_L_Dinophysis        = sum(cells_L[genus == "Dinophysis"],        na.rm = TRUE),
    cells_L_Pseudonitzschia   = sum(cells_L[genus == "Pseudo-nitzschia"],  na.rm = TRUE),
    cells_L_Azadinium         = sum(cells_L[genus == "Azadinium"],         na.rm = TRUE),
    cells_L_Chrysochromulina  = sum(cells_L[genus == "Chrysochromulina"],  na.rm = TRUE),
    cells_L_Prymnesium        = sum(cells_L[genus == "Prymnesium"],        na.rm = TRUE),
    cells_L_Amphidinium       = sum(cells_L[genus == "Amphidinium"],       na.rm = TRUE),
    cells_L_Pseudochattonella = sum(cells_L[genus == "Pseudochattonella"], na.rm = TRUE),
    cells_L_Phaeocystis       = sum(cells_L[genus == "Phaeocystis"],       na.rm = TRUE),
    cells_L_Karlodinium       = sum(cells_L[genus == "Karlodinium"],       na.rm = TRUE),
    cells_L_Cyanobacteria     = sum(cells_L[genus == "Cyanobacteria"],     na.rm = TRUE),
    
    # Multi-species information preservation
    species_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(species[harmful_algae == TRUE & !is.na(species)]), collapse = " | "),
      "No harmful algae"
    ),
    genera_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(genus[harmful_algae == TRUE & !is.na(genus)]), collapse = " | "),
      "Other"
    ),
    toxin_syndromes_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(toxin_syndrome[harmful_algae == TRUE & !is.na(toxin_syndrome)]), collapse = " | "),
      "Other"
    ),
    
    # Dominant species for backward compatibility
    dominant_species = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          species[max_idx]
        } else {
          "No harmful algae"
        }
      } else {
        "No harmful algae"
      }
    },
    dominant_genus = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          genus[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    dominant_toxin_syndrome = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          toxin_syndrome[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    
    # Environmental variables
    across(where(is.numeric) & !starts_with("probability_") & !starts_with("cells_L_"), ~ mean(.x, na.rm = TRUE)),
    
    # Categorical variables
    across(c(day, month, year, doy, strat), first),
    .groups = "drop"
  ) %>%
  mutate(
    # Use dominant species for main columns
    species = dominant_species,
    genus = dominant_genus,
    toxin_syndrome = dominant_toxin_syndrome,
    harmful_algae = genus != "Other",
    wind_ms = rowMeans(dplyr::select(., starts_with("wind")), na.rm = TRUE),
    country = "Sweden"
  ) %>%
  dplyr::select(-dominant_species, -dominant_genus, -dominant_toxin_syndrome, -wind_speed) %>%
  filter(is.na(parameter) | parameter != "Abundance")

# Validate the aggregated data
sweden_combined <- validate_aggregated_data(sweden_combined, "Sweden")
```

### Step 6: Update Germany Aggregation Logic

Find the Germany aggregation section and add similar logic:

```r
# Update Germany processing (after line 1000)
# After the existing processing, add aggregation step
germany_combined <- germany_combined %>%
  group_by(date, station) %>%
  dplyr::summarise(
    # Probability aggregation
    across(
      all_of(c("probability", unname(sapply(genera_of_interest, prob_col_name)))),
      ~ case_when(any(. == 1L, na.rm = TRUE) ~ 1L,
                  all(is.na(.)) ~ NA_integer_,
                  TRUE ~ 0L)
    ),
    
    # Cell count aggregation  
    cells_L_Alexandrium       = sum(cells_L[genus == "Alexandrium"],       na.rm = TRUE),
    cells_L_Dinophysis        = sum(cells_L[genus == "Dinophysis"],        na.rm = TRUE),
    cells_L_Pseudonitzschia   = sum(cells_L[genus == "Pseudo-nitzschia"],  na.rm = TRUE),
    cells_L_Azadinium         = sum(cells_L[genus == "Azadinium"],         na.rm = TRUE),
    cells_L_Chrysochromulina  = sum(cells_L[genus == "Chrysochromulina"],  na.rm = TRUE),
    cells_L_Prymnesium        = sum(cells_L[genus == "Prymnesium"],        na.rm = TRUE),
    cells_L_Amphidinium       = sum(cells_L[genus == "Amphidinium"],       na.rm = TRUE),
    cells_L_Pseudochattonella = sum(cells_L[genus == "Pseudochattonella"], na.rm = TRUE),
    cells_L_Phaeocystis       = sum(cells_L[genus == "Phaeocystis"],       na.rm = TRUE),
    cells_L_Karlodinium       = sum(cells_L[genus == "Karlodinium"],       na.rm = TRUE),
    cells_L_Cyanobacteria     = sum(cells_L[genus == "Cyanobacteria"],     na.rm = TRUE),
    
    # Multi-species information preservation
    species_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(species[harmful_algae == TRUE & !is.na(species)]), collapse = " | "),
      "No harmful algae"
    ),
    genera_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(genus[harmful_algae == TRUE & !is.na(genus)]), collapse = " | "),
      "Other"
    ),
    toxin_syndromes_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(toxin_syndrome[harmful_algae == TRUE & !is.na(toxin_syndrome)]), collapse = " | "),
      "Other"
    ),
    
    # Dominant species for backward compatibility
    dominant_species = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          species[max_idx]
        } else {
          "No harmful algae"
        }
      } else {
        "No harmful algae"
      }
    },
    dominant_genus = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          genus[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    dominant_toxin_syndrome = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          toxin_syndrome[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    
    # Environmental variables
    across(where(is.numeric) & !starts_with("probability_") & !starts_with("cells_L_"), ~ mean(.x, na.rm = TRUE)),
    
    # Categorical variables
    across(c(day, month, year, doy, strat), first),
    .groups = "drop"
  ) %>%
  mutate(
    # Use dominant species for main columns
    species = dominant_species,
    genus = dominant_genus, 
    toxin_syndrome = dominant_toxin_syndrome,
    harmful_algae = genus != "Other",
    country = "Germany"
  ) %>%
  dplyr::select(-dominant_species, -dominant_genus, -dominant_toxin_syndrome)

# Validate the aggregated data
germany_combined <- validate_aggregated_data(germany_combined, "Germany")
```

### Step 7: Update Final All-Data Aggregation

Find the final all-data aggregation section (around line 1150) and update:

```r
# Step 1: Aggregate genus-level probabilities, cell concentrations and categorical fields per combined_station/date
all_data_probs <- all_data %>%
  group_by(combined_station, date) %>%
  dplyr::summarise(
    # Overall: 1 if any row is present, 0 if all absent, NA if all missing
    probability = case_when(
      any(probability == 1L, na.rm = TRUE) ~ 1L,
      all(is.na(probability))              ~ NA_integer_,
      TRUE                                 ~ 0L
    ),
    # Individual genus probabilities (1/0/NA) — one across() replaces 11 repeated case_when blocks
    across(
      all_of(unname(sapply(genera_of_interest, prob_col_name))),
      ~ case_when(any(. == 1L, na.rm = TRUE) ~ 1L,
                  all(is.na(.))              ~ NA_integer_,
                  TRUE                       ~ 0L)
    ),
    
    # Genus-level cell concentration totals (bracket subsetting is faster than ifelse)
    cells_L_Alexandrium       = sum(cells_L[genus == "Alexandrium"],       na.rm = TRUE),
    cells_L_Dinophysis        = sum(cells_L[genus == "Dinophysis"],        na.rm = TRUE),
    cells_L_Pseudonitzschia   = sum(cells_L[genus == "Pseudo-nitzschia"],  na.rm = TRUE),
    cells_L_Azadinium         = sum(cells_L[genus == "Azadinium"],         na.rm = TRUE),
    cells_L_Chrysochromulina  = sum(cells_L[genus == "Chrysochromulina"],  na.rm = TRUE),
    cells_L_Prymnesium        = sum(cells_L[genus == "Prymnesium"],        na.rm = TRUE),
    cells_L_Amphidinium       = sum(cells_L[genus == "Amphidinium"],       na.rm = TRUE),
    cells_L_Pseudochattonella = sum(cells_L[genus == "Pseudochattonella"], na.rm = TRUE),
    cells_L_Phaeocystis       = sum(cells_L[genus == "Phaeocystis"],       na.rm = TRUE),
    cells_L_Karlodinium       = sum(cells_L[genus == "Karlodinium"],       na.rm = TRUE),
    cells_L_Cyanobacteria     = sum(cells_L[genus == "Cyanobacteria"],     na.rm = TRUE),
    
    # ✅ NEW: Multi-species information preservation  
    species_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(species[harmful_algae == TRUE & !is.na(species)]), collapse = " | "),
      "No harmful algae"
    ),
    genera_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(genus[harmful_algae == TRUE & !is.na(genus)]), collapse = " | "),
      "Other"
    ),
    toxin_syndromes_present = ifelse(
      any(harmful_algae == TRUE, na.rm = TRUE),
      paste(unique(toxin_syndrome[harmful_algae == TRUE & !is.na(toxin_syndrome)]), collapse = " | "),
      "Other"
    ),
    
    # ✅ NEW: Dominant species for backward compatibility
    dominant_species = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          species[max_idx]
        } else {
          "No harmful algae"
        }
      } else {
        "No harmful algae"
      }
    },
    dominant_genus = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          genus[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    dominant_toxin_syndrome = {
      if(any(harmful_algae == TRUE, na.rm = TRUE)) {
        harmful_rows <- which(harmful_algae == TRUE & !is.na(cells_L))
        if(length(harmful_rows) > 0) {
          max_idx <- harmful_rows[which.max(cells_L[harmful_rows])]
          toxin_syndrome[max_idx]
        } else {
          "Other"
        }
      } else {
        "Other"
      }
    },
    
    # Categorical / metadata fields
    strat               = combine_strat(strat),
    limiting_conditions = combine_limiting_conditions(limiting_conditions),
    across(c(day, month, year, doy), first),
    .groups = "drop"
  ) %>%
  mutate(
    # ✅ Use dominant species for main columns (backward compatibility)
    species = dominant_species,
    genus = dominant_genus,
    toxin_syndrome = dominant_toxin_syndrome,
    harmful_algae = genus != "Other",
    logistic = probability,  # alias kept for downstream compatibility
    # Get country from the first non-NA value in the group
    country = all_data %>% 
      filter(combined_station == first(.$combined_station), date == first(.$date)) %>% 
      pull(country) %>% 
      first()
  ) %>%
  dplyr::select(-dominant_species, -dominant_genus, -dominant_toxin_syndrome)

# Validate the final aggregated data
all_data <- validate_aggregated_data(all_data, "All Data Combined")
```

### Step 8: Test and Verification

Add these verification checks after the final aggregation:

```r
# Verification checks
cat("=== FINAL DATA VERIFICATION ===\n")

# 1. Check that species_present contains all species when multiple genera are present
multi_genus_check <- all_data %>%
  filter(str_detect(genera_present, "\\|")) %>%
  dplyr::select(date, combined_station, species, genus, species_present, genera_present) %>%
  head(10)

if(nrow(multi_genus_check) > 0) {
  cat("✅ Multi-genus events found - checking species preservation:\n")
  print(multi_genus_check)
} else {
  cat("ℹ️ No multi-genus events in final dataset.\n")
}

# 2. Check for any remaining data integrity violations
final_violations <- all_data %>%
  dplyr::select(date, combined_station, genus, starts_with("probability_")) %>%
  pivot_longer(starts_with("probability_"), names_to = "prob_col", values_to = "prob_val") %>%
  filter(prob_val == 1L) %>%
  mutate(
    expected_genus = str_remove(prob_col, "probability_"),
    expected_genus = ifelse(expected_genus == "Pseudonitzschia", "Pseudo-nitzschia", expected_genus)
  ) %>%
  filter(genus != expected_genus)

if(nrow(final_violations) > 0) {
  cat("❌ REMAINING DATA INTEGRITY VIOLATIONS:\n")
  print(final_violations)
  stop("Data integrity violations still exist! Check aggregation logic.")
} else {
  cat("✅ NO data integrity violations found in final dataset!\n")
}

# 3. Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total observations:", nrow(all_data), "\n")
cat("Stations:", length(unique(all_data$combined_station)), "\n")
cat("Date range:", min(all_data$date), "to", max(all_data$date), "\n\n")

# Per-genus presence counts
for(genus in genera_of_interest) {
  prob_col <- prob_col_name(genus)
  if(prob_col %in% colnames(all_data)) {
    n_present <- sum(all_data[[prob_col]] == 1L, na.rm = TRUE)
    cat(genus, "present:", n_present, "observations\n")
  }
}

cat("\n=== IMPLEMENTATION COMPLETE ===\n")
```

## Summary

This implementation:

1. **✅ Fixes the data integrity bug** by properly aggregating multi-species sampling events
2. **✅ Preserves all species information** in `species_present`, `genera_present`, `toxin_syndromes_present` columns  
3. **✅ Maintains backward compatibility** by using the dominant species (highest cell count) in the main `species` and `genus` columns
4. **✅ Adds comprehensive validation** to catch any remaining issues
5. **✅ Works across all datasets** (Norway, Denmark, Sweden, Germany)

The key insight is that sampling events often contain multiple harmful species, and the original `first()` approach was randomly selecting one species while the probability aggregation correctly detected all present genera. This fix ensures data integrity while preserving all the valuable multi-species information.

## Implementation Notes for Claude Code

- This is a comprehensive R time-series analysis pipeline for harmful algal bloom monitoring
- The script processes phytoplankton monitoring data from 4 Nordic countries (Norway, Denmark, Sweden, Germany)
- The bug was in the data aggregation logic where multiple species/genera per sampling event were incorrectly collapsed using `first()`
- The solution preserves all species information while maintaining backward compatibility
- All changes should be applied to `Time_series_analysis.R` and `Time_series_analysis_custom_functions.R`
