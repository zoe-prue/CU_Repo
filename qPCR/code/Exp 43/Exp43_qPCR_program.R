# my_qPCR_program_Exp_43_JAK1_levels_or_ISG_GAPDH
# ZMP
# 10/2/2024

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
install.packages("readxl")

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggbreak) 
library(ggforce)
library(tidyverse)
library(readxl)

# Load Functions -----------------------------

# this script has functions that are usable for any qPCR analysis
source("my_qPCR_functions_v2.R") 

# Global Variables ---------------------------

# (CHANGE THESE TO SUIT YOUR NEEDS)
# make sure your qPCR_data and your R scripts are in the same folder (directory) 
# and that this is set as your working directory
skip_meta = 20 # number of metadata rows to skip
qpcr_data <- read_excel("~/Desktop/CU_coding/qPCR/raw_data/Exp43_qPCR_data_combined.xlsx", 
                      skip = skip_meta) # put in path to YOUR data


# View the result
print(qpcr_data)

my_cols <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000") # should be sample length as probe_list

# Single Graph -------------------------------
# you're welcome to just make one plot and change it up each time!
# here is some simple code to make one plot
# Local variables:
probe = "IFI27" # list of probes of interest
reference = "B-ACTIN" # probe that you're "normalizing" to
ctrl_sample = "JAK1 WT 0 IFN" # sample you're making everything "relative" to
plot_title = " mRNA"
fill_color = "#648FFF"

##

# avg_cq_values <- qpcr_data %>%
#   mutate(Sample = str_trim(Sample)) %>%  # Trim whitespace
#   mutate(Sample = ifelse(Sample == "293T Clone 15", "293T Clone 22", Sample)) %>%
#   filter(!Sample %in% c("293T WT", "293T Clone 19", "293T Clone 14", "293T Clone 22")) %>%
#   group_by(Sample, Target, Well.Position) %>%
#   mutate(Cq = as.numeric(Cq)) %>%
#   summarize(Cq = mean(Cq, na.rm = TRUE))
# 
# # select only columns from qpcr_data that are needed and ensure that values are class numeric
# tidy_qpcr_data <- avg_cq_values %>%
#   # filter(!str_detect(Sample, "^pTRIP JAK1 E139K 1$")) %>%
#   filter(Target != is.na(Target) & Target != "") %>%
#   select(Well.Position, Sample, Target, Cq) %>%
#   group_by(Well.Position) %>% 
#   mutate(Cq = as.numeric(Cq))
# # filter(Cq <= 34)


#### FOR JAK1 Analysis
# tidy_qpcr_data <- qpcr_data %>%
#   mutate(Sample = str_trim(Sample)) %>%
#   # mutate(Sample = ifelse(Sample == "293T Clone 15", "293T Clone 22", Sample)) %>%
#   # filter(!Sample %in% c("293T WT", "293T Clone 19", "293T Clone 14", "293T Clone 22")) %>%
#   group_by(Sample, Target, Well.Position) %>%
#   mutate(Cq = as.numeric(Cq)) %>%
#   select(Well.Position, Sample, Target, Cq)


##### for ISG Analysis by IFN stim type

tidy_qpcr_data <- qpcr_data %>%
  mutate(Sample = str_trim(Sample)) %>%
  # mutate(Sample = ifelse(Sample == "293T Clone 15", "293T Clone 22", Sample)) %>%
  # filter(!Sample %in% c("293T WT", "293T Clone 19", "293T Clone 14", "293T Clone 22")) %>%
  filter(str_ends(Sample, "IFNL2") | str_ends(Sample, "0 IFN")) %>%  # Filter targets that end with "IFN a2b" or "0 IFN"
  group_by(Sample, Target, Well.Position) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  select(Well.Position, Sample, Target, Cq)


# create a data table with only Cq values from target being tested (e.g. JAK1)
test_data <- tidy_qpcr_data %>%
  filter(Target == probe)

# create a data table with only Cq values from housekeeping gene (e.g. GAPDH)
ref_data <- tidy_qpcr_data %>%
  filter(Target == reference) %>%
  rename("ref_target" = Target) %>%
  rename("ref_Cq" = Cq)

# Group by Sample and calculate the average Cq
test_data <- test_data %>%
  group_by(Sample, Target) #%>%
# summarize(Cq = mean(Cq))
# # Remove the last number from Sample column to group the data correctly
# test_data <- test_data %>%
#   mutate(Sample = gsub("\\s[0-9]$", "", Sample))

# Step 1: Calculate the average Cq for the reference gene (JAK1) across technical replicates for each Sample
ref_data_avg <- ref_data %>%
  group_by(Sample) %>%
  summarize(avg_ref_Cq = mean(ref_Cq, na.rm = TRUE))

# Step 2: Join the test data (target gene, e.g., IFI27) with the averaged reference data
delta_cq_data <- inner_join(test_data, ref_data_avg, by = "Sample") %>%
  mutate(delta_Cq = Cq - avg_ref_Cq)  # Calculate delta_Cq using the average reference Cq

# Check the result
print(delta_cq_data)

# # Merge the two data frames based on the Well.Position column
# delta_cq_data <- merge(test_data, ref_data, by = "Well.Position")
# # Calculate the difference between Cq and ref_Cq
# delta_cq_data <- delta_cq_data %>%
#   mutate(delta_Cq = Cq - ref_Cq)
# print(delta_cq_data)

# Filter data for "pTRIP JAK1 WT" samples
pTRIP_JAK1_WT <- delta_cq_data %>%
  filter(grepl("JAK1 WT 0 IFN", Sample))

pTRIP_JAK1_WT

# Calculate average delta_Cq, handling NA values
average_delta_Cq_ref <- mean(pTRIP_JAK1_WT$delta_Cq, na.rm = TRUE)

# Remove numbers at the end of Sample names
# pTRIP_JAK1_WT$Sample <- gsub("\\s[0-9]+$", "", pTRIP_JAK1_WT$Sample)

# Output the vector of average delta_Cq
average_delta_Cq_ref

# # find the mean delta_cq for your ctrl_sample (to which every other sample is "relative")
# mean_ctrl_sample_delta <- mean(delta_cq_data$delta_Cq[delta_cq_data$Sample == ctrl_sample],
#                                na.rm = TRUE)

# subtract mean delta_cq for your ctrl_sample from all delta_cq values
delta_delta_cq_data <- delta_cq_data %>%
  mutate(delta_delta_cq = delta_Cq - average_delta_Cq_ref)

# find the relative concentrations of the target being tested for every sample
rel_conc_data <- delta_delta_cq_data %>%
  mutate(rel_conc = 2^(-(delta_delta_cq)))


# graph_rel_conc <- function(rel_conc_data, fill_color, probe) {
#   
#   # Filter out rows with NA values in 'rel_conc'
#   rel_conc_data <- rel_conc_data[complete.cases(rel_conc_data$rel_conc), ]
#   
#   # Extract group names without numbers at the end
#   rel_conc_data$group <- gsub("\\s[0-9]+$", "", rel_conc_data$Sample)
#   
#   # Aggregate data to calculate mean and standard deviation for each group of three points
#   graphable_data <- rel_conc_data %>%
#     group_by(group) %>%
#     summarise(mean_rel_conc = mean(rel_conc, na.rm = TRUE),
#               sde_rel_conc = sd(rel_conc, na.rm = TRUE) / sqrt(n()))
#   
#   # Make a bar graph of your rel_conc_data with error bars and data points
#   rel_conc_bar <- ggplot(graphable_data, aes(x = group, y = mean_rel_conc)) +
#     theme_classic() +
#     geom_bar(stat = "identity", fill = fill_color, width = 0.5, color = "black") +
#     geom_errorbar(aes(ymin = mean_rel_conc - sde_rel_conc, ymax = mean_rel_conc + sde_rel_conc),
#                   width = 0.2, position = position_dodge(width = 0.5)) +
#     geom_point(data = rel_conc_data, aes(x = group, y = rel_conc), position = position_jitter(width = 0.2), alpha = 0.8, color = "black", size = 1) +
#     theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,), axis.title.x = element_blank(), plot.title=element_text(hjust=0.5)) +
#     xlab("JAK1") +
#     scale_y_continuous(trans = "log10") +
#     ylab(paste0(probe, "/", "GAPDH", " Relative Units")) +
#     theme(plot.title = element_text(size = 16),
#           axis.title.x=element_text(size=14),
#           axis.title.y=element_text(size=14))
#   
#   return(rel_conc_bar)
# }
# 
# # one graph with Path #1:
# # rel_conc_data <- rel_conc_p1(avg_cq_values, probe, reference, ctrl_sample)
# rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color, probe) + ggtitle(paste0(probe, plot_title))
# rel_conc_bar

graph_rel_conc <- function(rel_conc_data, fill_color, probe) {
  
  # Filter out rows with NA values in 'rel_conc'
  rel_conc_data <- rel_conc_data[complete.cases(rel_conc_data$rel_conc), ]
  
  # Extract group names without numbers at the end
  rel_conc_data$group <- gsub("\\s[0-9]+$", "", rel_conc_data$Sample)
  
  # Make a bar graph with individual data points for each technical replicate
  rel_conc_bar <- ggplot(rel_conc_data, aes(x = group, y = rel_conc)) +
    theme_classic() +
    geom_bar(stat = "identity", fill = fill_color, width = 0.5, color = "black", position = position_dodge(width = 0.7)) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.8, color = "black", size = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    xlab(probe) +
    scale_y_continuous(trans = "log10") +
    ylab(paste0(probe, "/", "B-ACTIN", " Relative Units")) +
    theme(plot.title = element_text(size = 16), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
  
  return(rel_conc_bar)
}

# Generate the plot
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color, probe) + ggtitle(paste0(probe, " mRNA"))
rel_conc_bar

