# my_qPCR_program_Exp_11
# ZMP
# 04/08/2024

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
qpcr_data <- read.csv("~/Desktop/CU_coding/qPCR/raw_data/exp11_4.19.2024_u4cs_real_Copy_20240419_160053_Results_20240419 173509.csv", 
                      skip = skip_meta, 
                      stringsAsFactors = FALSE) # put in path to YOUR data

# There are two scenarios; either you use GAPDH vic in the same well, or 18S fam in separate wells
# let's call these Path #1 and #2, respectively

# Path #1 -------------------------------

# Local variables 
# (CHANGE THESE TO SUIT YOUR NEEDS)
probe_list = c("JAK1") # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "pTRIP Luciferase" # sample you're making everything "relative" to
my_cols <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000") # should be sample length as probe_list
file_name <- "exp10" # name of file you want to make for plots
plot_title = " Baseline Expression in U4Cs"
my_plots <- list() # make an empty list to put plots into

# # create list of probes aside from your reference probe
# probe_list <- qpcr_data %>%
#   filter(Target != reference) %>%
#   distinct(Target) %>%
#   pull(Target)

# make multiple plots at a time and save them into a folder
if (!file.exists(file_name)) { # this creates a folder to put your plots in
  dir.create(file_name) # as long as a folder with that name does not already exist
}

for (probe in probe_list) { # cycling through each probe in the list of probes of interest
  rel_conc_data <- rel_conc_p1(qpcr_data, probe, reference, ctrl_sample) # calculate the relative concentrations of your probe in each sample
  color_index <- which(probe_list == probe)  # cycle through a different color for each new probe
  # now graph the relative concentrations to see baseline expression of gene of interest, using a specific color
  rel_conc_bar <- graph_rel_conc(rel_conc_data, my_cols[color_index]) + ggtitle(paste0(probe, plot_title))
  # rel_conc_bar <- stat_compare_means(method = "anova", label = "p.format", hide.ns = TRUE)x
  my_plots[[probe]] <- rel_conc_bar # add each new plot to a list
  ggsave(filename = paste0(file_name, "/plot_", probe, ".png"), plot = rel_conc_bar, width = 10, height = 6) # save each new plot to a shared folder
}
rel_conc_bar

# Path #2 -------------------------------

# Local variables 
# (CHANGE THESE TO SUIT YOUR NEEDS)
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 Luc/0" # sample you're making everything "relative" to
my_cols <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000") # should be sample length as probe_list
file_name <- "Exp8_RSAD2_norm_JAK1" # name of file you want to make for plots
plot_title = " Baseline Expression in U4Cs"
my_plots <- list() # make an empty list to put plots into

# create list of probes aside from your reference probe
probe_list <- qpcr_data %>%
  filter(Target != reference) %>%
  distinct(Target) %>%
  pull(Target)

# make multiple plots at a time and save them into a folder
if (!file.exists(file_name)) { # this creates a folder to put your plots in
  dir.create(file_name) # as long as a folder with that name does not already exist
}

for (probe in probe_list) { # cycling through each probe in the list of probes of interest
  
  rel_conc_data <- rel_conc_p2(qpcr_data, probe, reference, ctrl_sample) # calculate the relative concentrations of your probe in each sample
  color_index <- which(probe_list == probe)  # cycle through a different color for each new probe
  # now graph the relative concentrations to see baseline expression of gene of interest
  rel_conc_bar <- graph_rel_conc(rel_conc_data, my_cols[color_index]) + ggtitle(paste0(probe, plot_title))
  # rel_conc_bar <- stat_compare_means(method = "anova", label = "p.format", hide.ns = TRUE)x
  my_plots[[probe]] <- rel_conc_bar
  ggsave(filename = paste0(file_name, "/plot_", probe, ".png"), plot = rel_conc_bar, width = 10, height = 6)
  
}

# Single Graph -------------------------------
# you're welcome to just make one plot and change it up each time!
# here is some simple code to make one plot
# Local variables:
probe = "IFI27" # list of probes of interest
reference = "JAK1" # probe that you're "normalizing" to
ctrl_sample = "pTRIP JAK1 WT" # sample you're making everything "relative" to
plot_title = " mRNA"
fill_color = "#785EF0"

# data pre-processing - averaging technical replicates
# Group the data by the "Target" column and calculate the mean Cq for each group
qpcr_data$Cq <- as.numeric(qpcr_data$Cq)
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
  filter(grepl("pTRIP JAK1 WT", Sample)) ####CHANGE THIS each time
# Calculate average delta_Cq
average_delta_Cq_ref <- mean(pTRIP_JAK1_WT$delta_Cq)
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
    summarise(mean_rel_conc = mean(rel_conc),
              sde_rel_conc = sd(rel_conc) / sqrt(n()))
  
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

# rel_conc_bar <- rel_conc_bar + scale_x_discrete(labels = c(
#   "pTRIP JAK1 WT" = "WT", 
#   "pTRIP JAK1 S703I" = "S703I", 
#   "pTRIP JAK1 D82G" = "D82G"))
# rel_conc_bar

desired_order=(c("pTRIP JAK1 WT", 
                 "pTRIP JAK1 S703I", 
                 "pTRIP JAK1 M1085I"))

rel_conc_bar_order <- rel_conc_bar + 
  scale_x_discrete(labels = c(
    "pTRIP JAK1 WT" = "WT", 
    "pTRIP JAK1 S703I" = "S703I", 
    "pTRIP JAK1 M1085I" = "M1085I"), 
    limits = desired_order
  ) +
  # scale_y_break(c(2.5, 10)) +
  theme(axis.text.x = element_text(size = 14))
rel_conc_bar_order
# Reorder x-axis labels

  


# rel_conc_data <- rel_conc_p2(qpcr_data, probe, reference, ctrl_sample)
# rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color) + ggtitle(paste0(probe, plot_title))
# rel_conc_bar


