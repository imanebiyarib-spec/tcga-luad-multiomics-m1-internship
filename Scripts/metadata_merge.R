# ==============================================================================
# Script Name: Clinical Data Enrichment and Visualization
# Description: This script loads raw sample, clinical, and exposure data (e.g., 
#              from TCGA-LUAD). It deduplicates patient metadata to prevent 
#              unintended row expansion, merges the information into a single 
#              enriched sample sheet, and standardizes key clinical variables. 
#              Finally, it generates 100% stacked bar charts to compare 
#              demographics (Smoking Status, Gender) across Healthy vs. Tumor tissues.
# Author:      Imane BIYAR
# Date:        July 03, 2026
#
# Inputs:      - sample_sheet.tsv
#              - clinical.tsv
#              - exposure.tsv
# Outputs:     - metadata.csv (Enriched sample sheet)
#              - BarPlot_Smoking_Status.png
#              - BarPlot_Gender.png
# ==============================================================================

# Load the libraries
library(dplyr)
library(ggplot2)
library(stringr)

# Load raw data files
samples  <- read.delim("Data/Metadata_Clinical/sample_sheet.tsv", stringsAsFactors = FALSE)
clinical <- read.delim("Data/Metadata_Clinical/clinical.tsv", stringsAsFactors = FALSE)
exposure <- read.delim("Data/Metadata_Clinical/exposure.tsv", stringsAsFactors = FALSE)

# Deduplicate clinical and exposure data
# We use distinct() to ensure there is only one row per patient. 
# This prevents the creation of duplicate sample rows during the join.
clin_minimal <- clinical %>% 
  select(cases.submitter_id, demographic.gender) %>%
  distinct(cases.submitter_id, .keep_all = TRUE)

expo_minimal <- exposure %>% 
  select(cases.submitter_id, exposures.tobacco_smoking_status) %>%
  distinct(cases.submitter_id, .keep_all = TRUE)

# Merge datasets
# Joining the cleaned clinical and exposure data to the sample sheet.
# Because of the deduplication above, this will maintain the correct row count (601 rows).
sample_sheet_enrichi <- samples %>%
  left_join(clin_minimal, by = c("Case.ID" = "cases.submitter_id")) %>%
  left_join(expo_minimal, by = c("Case.ID" = "cases.submitter_id"))

# Standardize and categorize smoking status
sample_sheet_enrichi <- sample_sheet_enrichi %>%
  mutate(
    exposures.tobacco_smoking_status = case_when(
      grepl("Lifelong Non-Smoker", exposures.tobacco_smoking_status, ignore.case = TRUE) ~ "Non_Smoker",
      grepl("Current Smoker", exposures.tobacco_smoking_status, ignore.case = TRUE) ~ "Smoker_(Current)",
      grepl("Current Reformed Smoker", exposures.tobacco_smoking_status, ignore.case = TRUE) ~ "Smoker_(Reformed)",
      
      # Catch-all for empty values, NAs, or "--"
      TRUE ~ "Unknown"
    )
  )

# Prepare data for visualization
plot_data <- sample_sheet_enrichi %>%
  mutate(
    # Create condition labels (Healthy vs Tumor) for the X-axis
    Condition = case_when(
      grepl("Normal", Tissue.Type, ignore.case = TRUE) ~ "Healthy",
      grepl("Tumor", Tissue.Type, ignore.case = TRUE) ~ "Tumor",
      TRUE ~ "Other"
    ),
    
    # Clean and standardize gender labels
    Gender = case_when(
      is.na(demographic.gender) | demographic.gender == "'--" ~ "Unknown",
      TRUE ~ str_to_title(demographic.gender)
    )
  ) %>%
  # Filter out any 'Other' tissue types to only plot Healthy and Tumor samples
  filter(Condition %in% c("Healthy", "Tumor"))


# ==============================================================================
# Plot 1: Smoking Status
# ==============================================================================
smoking_counts <- plot_data %>%
  count(Condition, exposures.tobacco_smoking_status) %>%
  group_by(Condition) %>%
  mutate(Percentage = n / sum(n) * 100)

plot_smoking <- ggplot(smoking_counts, aes(x = Condition, y = Percentage, fill = exposures.tobacco_smoking_status)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%\n(n=", n, ")"), "")), 
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 4) +
  # Assign specific colors to the 4 newly created smoking categories
  scale_fill_manual(values = c("Non_Smoker" = "darkgreen", 
                               "Smoker_(Current)" = "orange",
                               "Smoker_(Reformed)" = "red", 
                               "Unknown" = "grey")) +
  theme_minimal(base_size = 14) +
  labs(title = "Statut Tabagique : Tissu Sain vs Tumeur", 
       x = NULL, y = "Pourcentage (%)", fill = "Statut :") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("Results/Metadata/BarPlot_Smoking_Status.png", plot = plot_smoking, width = 7, height = 6, dpi = 300)



#### Overall plot 
overall_smoking <- sample_sheet_enrichi %>%
  distinct(Case.ID, .keep_all = TRUE) %>%
  count(exposures.tobacco_smoking_status) %>%
  mutate(
    Percentage = n / sum(n) * 100,
    # Create the exact text label seen in the image: n=XXX \n (YY.Y%)
    label_text = paste0("n=", n, "\n(", round(Percentage, 1), "%)")
  ) %>%
  # Force the specific order from the image
  mutate(exposures.tobacco_smoking_status = factor(
    exposures.tobacco_smoking_status, 
    levels = c("Non_Smoker", "Smoker_(Current)", "Smoker_(Reformed)", "Unknown")
  ))

# Build the plot
plot_smoking <- ggplot(overall_smoking, aes(x = exposures.tobacco_smoking_status, y = n, fill = exposures.tobacco_smoking_status)) +
  # Use geom_bar with black outlines
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, width = 0.6) +
  # Add the text labels slightly above the bars
  geom_text(aes(label = label_text), vjust = -0.5, fontface = "bold", size = 3.5) +
  # Colors matching the image
  scale_fill_manual(values = c(
    "Non_Smoker" = "darkgreen", 
    "Smoker_(Current)" = "orange", 
    "Smoker_(Reformed)" = "red", 
    "Unknown" = "grey"
  )) +
  # Clean up the X-axis labels to match the image text
  scale_x_discrete(labels = c(
    "Non_Smoker" = "Non-Smoker", 
    "Smoker_(Current)" = "Smoker (Current)", 
    "Smoker_(Reformed)" = "Smoker (Reformed)", 
    "Unknown" = "Unknown/Not Reported"
  )) +
  # Expand the Y-axis slightly so the top labels don't get cut off
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  # Use theme_classic for the L-shaped axis lines and white background
  theme_classic(base_size = 14) +
  labs(
    title = "Overall Smoking Demographics (Entire TCGA-LUAD Cohort)", 
    x = NULL, 
    y = "Total Number of Patients"
  ) +
  theme(
    legend.position = "none", # Remove legend
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 20)),
    # Angle the x-axis text just like the image
    axis.text.x = element_text(angle = 15, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(color = "black")
  )

ggsave("Results/Metadata/BarPlot_Smoking_Overall.png", plot = plot_smoking, width = 9, height = 6, dpi = 300)

# ==============================================================================
# Plot 2: Gender Distribution
# ==============================================================================
gender_counts <- plot_data %>%
  count(Condition, Gender) %>%
  group_by(Condition) %>%
  mutate(Percentage = n / sum(n) * 100)

plot_gender <- ggplot(gender_counts, aes(x = Condition, y = Percentage, fill = Gender)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%\n(n=", n, ")"), "")), 
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Male" = 'blue', "Female" = "red", "Unknown" = "grey")) +
  theme_minimal(base_size = 14) +
  labs(title = "Répartition par Genre : Tissu Sain vs Tumeur", 
       x = NULL, y = "Pourcentage (%)", fill = "Genre :") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("Results/Metadata/BarPlot_Gender.png", plot = plot_gender, width = 7, height = 6, dpi = 300)

# Save the enriched metadata to CSV
write.csv(sample_sheet_enrichi, "Results/Metadata/metadata.csv", row.names = FALSE)