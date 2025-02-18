# my_qPCR_program_Exp_25
# ZMP
# 08/30/2024

# Description --------------------------------
# This is a program designed to quickly interpret qPCR data via the Livak Method 
# and produce graphical representations of each step of the analysis
# make sure to set your current directory as your working directory!

# to do
# facet graphs
# p values
# stars
# for loop to cycle through probes, put probe targets in a list - different colors per graph
# combine functions into one
# input their output locations for figure - automatically make a folder?
# combine functions into one
# change Sample.X into Sample - track development of name throughout


# Initialize ---------------------------------

install.packages("tidyr")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggpubr")
install.packages("ggbreak")
install.packages("ggforce")

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggbreak) 
library(ggforce)

# Load Functions -----------------------------

# this script has functions that are usable for any qPCR analysis
source("my_qPCR_functions_v2.R") 

# Global Variables ---------------------------

# (CHANGE THESE TO SUIT YOUR NEEDS)
# make sure your qPCR_data and your R scripts are in the same folder (directory) 
# and that this is set as your working directory
skip_meta = 20 # number of metadata rows to skip
qpcr_data <- read.csv("~/Desktop/CU_coding/qPCR/raw_data/Exp28_ZP_Exp6_CM_9.20.24_U4Cs_real_20240920_122151_Results_20240920 135718.csv", 
                      skip = skip_meta, 
                      stringsAsFactors = FALSE) # put in path to YOUR data

my_cols <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000") # should be sample length as probe_list

# Load necessary libraries
library(dplyr)

# Load necessary libraries
library(dplyr)

# Local variables
probe = "RSAD2"  # The probe of interest
reference = "GAPDH"  # The housekeeping gene
ctrl_sample = "pTRIP JAK1 WT 0 u/mL a2b"  # Control sample for relative calculations
plot_title = "mRNA"  # Plot title
fill_color = "#DC267F"  # Fill color for plots

# Step 1: Prepare the data
tidy_qpcr_data <- qpcr_data %>%
  filter(!is.na(Target) & Target != "") %>%
  select(Well.Position, Sample, Target, Cq) %>%
  mutate(Cq = as.numeric(Cq))  # Ensure Cq values are numeric

# Step 2: Create a data table with only Cq values from the probe being tested (e.g., IFI27)
test_data <- tidy_qpcr_data %>%
  filter(Target == probe) %>%
  select(Sample, Well.Position, Cq)  # Keep individual replicates

# Step 3: Create a data table with only Cq values from the reference gene (e.g., GAPDH)
ref_data <- tidy_qpcr_data %>%
  filter(Target == reference) %>%
  select(Sample, Well.Position, Cq) %>%
  rename(ref_Cq = Cq)  # Rename for clarity

# Step 4: Join the test_data and ref_data data frames without aggregating replicates
delta_cq_data <- inner_join(test_data, ref_data, by = c("Sample", "Well.Position")) %>%
  mutate(delta_Cq = Cq - ref_Cq)  # Calculate delta Cq for each technical replicate

# Filter data for "pTRIP JAK1 WT" samples
pTRIP_JAK1_WT <- delta_cq_data %>%
  filter(grepl("pTRIP JAK1 WT 0 u/mL a2b", Sample))

# Step 6: Calculate average delta_Cq for the control sample, handling NA values
average_delta_Cq_ref <- mean(pTRIP_JAK1_WT$delta_Cq, na.rm = TRUE)

# Step 7: Calculate delta_delta_Cq for all replicates
delta_cq_data <- delta_cq_data %>%
  mutate(delta_delta_cq = delta_Cq - average_delta_Cq_ref)

# Step 8: Find the relative concentrations of the target being tested for every sample
rel_conc_data <- delta_cq_data %>%
  mutate(rel_conc = 2^(-delta_delta_cq))  # Calculate relative concentrations

rel_conc_data <- rel_conc_data %>%
  filter(Sample != "pTRIP Luciferase 0 u/mL IFN a2b") %>%
  filter(Sample != "pTRIP Luciferase 10 u/mL a2b") %>%
  filter(Sample != "pTRIP Luciferase 100 u/mL a2b")
# filter(Sample != "pTRIP JAK1 S703I 1") %>%
# filter(Sample != "pTRIP JAK1 S703I 2") %>%
# filter(Sample != "pTRIP JAK1 S703I 3")

# Output final data with individual replicates
rel_conc_data


#############################

# Local variables
probe = "RSAD2"  # The probe of interest
reference = "JAK1"  # The reference gene
ctrl_sample = "pTRIP JAK1 WT 0 u/mL a2b"  # Control sample for relative calculations
plot_title = "mRNA"  # Plot title
fill_color = "#DC267F"  # Fill color for plots

# Step 1: Prepare the data
tidy_qpcr_data <- qpcr_data %>%
  filter(!is.na(Target) & Target != "") %>%
  select(Well.Position, Sample, Target, Cq) %>%
  mutate(Cq = as.numeric(Cq))  # Ensure Cq values are numeric

# Step 2: Create a data table with only Cq values from the probe being tested (IFI27)
test_data <- tidy_qpcr_data %>%
  filter(Target == probe) %>%
  select(Sample, Cq)  # Keep individual replicates, removing Well.Position

# Step 3: Create a data table with only Cq values from the reference gene (JAK1)
ref_data <- tidy_qpcr_data %>%
  filter(Target == reference) %>%
  select(Sample, Cq) %>%
  rename(ref_Cq = Cq)  # Rename for clarity

# Step 4: Join the test_data and ref_data data frames without aggregating replicates
# Here we need to summarize the reference data per sample since they are in different wells
ref_data_summary <- ref_data %>%
  group_by(Sample) %>%
  summarize(ref_Cq = mean(ref_Cq, na.rm = TRUE))  # Average the reference Cq values for each sample

delta_cq_data <- inner_join(test_data, ref_data_summary, by = "Sample") %>%
  mutate(delta_Cq = Cq - ref_Cq)  # Calculate delta Cq for each technical replicate

# Step 5: Filter data for "pTRIP JAK1 WT" samples
pTRIP_JAK1_WT <- delta_cq_data %>%
  filter(grepl(ctrl_sample, Sample))  # Use ctrl_sample for filtering

# Step 6: Calculate average delta_Cq for the control sample, handling NA values
average_delta_Cq_ref <- mean(pTRIP_JAK1_WT$delta_Cq, na.rm = TRUE)

# Step 7: Calculate delta_delta_Cq for all replicates
delta_cq_data <- delta_cq_data %>%
  mutate(delta_delta_cq = delta_Cq - average_delta_Cq_ref)

# Step 8: Find the relative concentrations of the target being tested for every sample
rel_conc_data <- delta_cq_data %>%
  mutate(rel_conc = 2^(-delta_delta_cq))  # Calculate relative concentrations

# Step 9: Filter out specific samples
rel_conc_data <- rel_conc_data %>%
  filter(!Sample %in% c("pTRIP Luciferase 0 u/mL IFN a2b",
                        "pTRIP Luciferase 10 u/mL a2b",
                        "pTRIP Luciferase 100 u/mL a2b"))

# Output final data with individual replicates
rel_conc_data

