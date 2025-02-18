# my_qPCR_program_Exp_30
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
qpcr_data <- read.csv("~/Desktop/CU_coding/qPCR/raw_data/Exp31_11.4.2024_U4Cs_20241104_164931_Results_20241104 182327.csv", 
                      skip = skip_meta, 
                      stringsAsFactors = FALSE) # put in path to YOUR data

my_cols <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000") # should be sample length as probe_list

# Single Graph -------------------------------
# you're welcome to just make one plot and change it up each time!
# here is some simple code to make one plot
# Local variables:
probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "pTRIP JAK1 WT" # sample you're making everything "relative" to
plot_title = " mRNA"
fill_color = "#DC267F"

avg_cq_values <- qpcr_data %>%
  group_by(Sample, Target, Well.Position) %>%
  summarize(Cq = mean(Cq, na.rm = TRUE))
# mutate(Sample = gsub("\\s\\d+$", "", Sample)) %>%
# filter(Cq <= 34)

# select only columns from qpcr_data that are needed and ensure that values are class numeric
tidy_qpcr_data <- qpcr_data %>%
  filter(Target != is.na(Target) & Target != "") %>%
  select(Well.Position, Sample, Target, Cq) %>%
  group_by(Well.Position) %>% 
  mutate(Cq = as.numeric(Cq))
# filter(Cq <= 34)

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
  group_by(Sample, Target) %>%
  summarize(Cq = mean(Cq))
# # Remove the last number from Sample column to group the data correctly
# test_data <- test_data %>%
#   mutate(Sample = gsub("\\s[0-9]$", "", Sample))

# Group by Sample and calculate the average Cq
ref_data <- ref_data %>%
  group_by(Sample, ref_target) %>%
  summarize(ref_Cq = mean(ref_Cq))

# join the test_data and ref_data data frames and form a new column with the calculated delta_cq for each well
delta_cq_data <- inner_join(test_data,
                            ref_data,
                            join_by("Sample"),
                            relationship = "one-to-one") %>%
  group_by(Sample) %>%
  mutate(delta_Cq = Cq - ref_Cq)
# select(-ends_with(".y")) %>%
# rename(Sample = Sample.x)

# # Merge the two data frames based on the Well.Position column
# delta_cq_data <- merge(test_data, ref_data, by = "Well.Position")
# # Calculate the difference between Cq and ref_Cq
# delta_cq_data <- delta_cq_data %>%
#   mutate(delta_Cq = Cq - ref_Cq)
# print(delta_cq_data)

# Filter data for "pTRIP JAK1 WT" samples
pTRIP_JAK1_WT <- delta_cq_data %>%
  filter(grepl("pTRIP JAK1 WT", Sample))

# Calculate average delta_Cq, handling NA values
average_delta_Cq_ref <- mean(pTRIP_JAK1_WT$delta_Cq, na.rm = TRUE)

# Remove numbers at the end of Sample names
pTRIP_JAK1_WT$Sample <- gsub("\\s[0-9]+$", "", pTRIP_JAK1_WT$Sample)

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

rel_conc_data <- rel_conc_data %>%
  filter(Sample != "pTRIP Luciferase 1") %>%
  filter(Sample != "pTRIP Luciferase 2") %>%
  filter(Sample != "pTRIP Luciferase 3")
# filter(Sample != "pTRIP JAK1 S703I 1") %>%
# filter(Sample != "pTRIP JAK1 S703I 2") %>%
# filter(Sample != "pTRIP JAK1 S703I 3")

graph_rel_conc <- function(rel_conc_data, fill_color, probe) {
  
  # Filter out rows with NA values in 'rel_conc'
  rel_conc_data <- rel_conc_data[complete.cases(rel_conc_data$rel_conc), ]
  
  # Extract group names without numbers at the end
  rel_conc_data$group <- gsub("\\s[0-9]+$", "", rel_conc_data$Sample)
  
  # Aggregate data to calculate mean and standard deviation for each group of three points
  graphable_data <- rel_conc_data %>%
    group_by(group) %>%
    summarise(mean_rel_conc = mean(rel_conc, na.rm = TRUE),
              sde_rel_conc = sd(rel_conc, na.rm = TRUE) / sqrt(n()))
  
  # Make a bar graph of your rel_conc_data with error bars and data points
  rel_conc_bar <- ggplot(graphable_data, aes(x = group, y = mean_rel_conc)) +
    theme_classic() +
    geom_bar(stat = "identity", fill = fill_color, width = 0.5, color = "black") +
    geom_errorbar(aes(ymin = mean_rel_conc - sde_rel_conc, ymax = mean_rel_conc + sde_rel_conc),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(data = rel_conc_data, aes(x = group, y = rel_conc), position = position_jitter(width = 0.2), alpha = 0.8, color = "black", size = 1) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,), axis.title.x = element_blank(), plot.title=element_text(hjust=0.5)) +
    xlab("JAK1") +
    scale_y_continuous(trans = "log10") +
    ylab(paste0(probe, "/", reference, " Relative Units")) +
    theme(plot.title = element_text(size = 16),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))
  
  return(rel_conc_bar)
}


# one graph with Path #1:
# rel_conc_data <- rel_conc_p1(avg_cq_values, probe, reference, ctrl_sample)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color, probe) + ggtitle(paste0(probe, plot_title))
rel_conc_bar

desired_order=(c("pTRIP JAK1 WT",
                 "pTRIP JAK1 S703I", 
                 "pTRIP JAK1 E139K",
                 "pTRIP JAK1 S153T",
                 "pTRIP JAK1 S512L",
                 "pTRIP JAK1 R826S",
                 "pTRIP JAK1 Q1016H"
))

rel_conc_bar_order <- rel_conc_bar + 
  scale_x_discrete(labels = c(
    "pTRIP JAK1 WT" = "WT", 
    "pTRIP JAK1 S703I" = "S703I",
    "pTRIP JAK1 E139K" = "E139K",
    "pTRIP JAK1 S153T" = "S153T",
    "pTRIP JAK1 S512L" = "S512L",
    "pTRIP JAK1 R826S" = "R826S",
    "pTRIP JAK1 Q1016H" = "Q1016H"
  ), 
  limits = desired_order
  ) +
  # scale_y_break(c(15, 75)) +
  theme(axis.text.x = element_text(size = 14)) 
# geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +  # Add dotted line
# geom_hline(yintercept = 1.1518898, linetype = 'dashed', color = 'black') +  # Add dotted line
# geom_hline(yintercept = 10.6, linetype = 'dashed', color = 'black')  # Add dotted line
rel_conc_bar_order

