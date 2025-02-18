# my_qPCR_program_v2
# ZMP
# 02/12/2024

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

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Load Functions -----------------------------

# this script has functions that are usable for any qPCR analysis
source("my_qPCR_functions_v2.R") 

# Global Variables ---------------------------

# (CHANGE THESE TO SUIT YOUR NEEDS)
# make sure your qPCR_data and your R scripts are in the same folder (directory) 
# and that this is set as your working directory
skip_meta = 23 # number of metadata rows to skip
qpcr_data <- read.csv("~/Desktop/CU_coding/qPCR/raw_data/Exp9_real_4.4.24_20240404_122458_Results_20240404_141702.csv", 
                      skip = skip_meta, 
                      stringsAsFactors = FALSE) # put in path to YOUR data

# There are two scenarios; either you use GAPDH vic in the same well, or 18S fam in separate wells
# let's call these Path #1 and #2, respectively

# Path #1 -------------------------------

# Local variables 
# (CHANGE THESE TO SUIT YOUR NEEDS)
# probe_list = c("JAK1", "RSAD2", "IFI27") # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "pTRIP Luciferase" # sample you're making everything "relative" to
my_cols <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000") # should be sample length as probe_list
file_name <- "Exp8_plots_norm_JAK1_500s" # name of file you want to make for plots
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
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 Luc/0" # sample you're making everything "relative" to
plot_title = " Baseline Expression in U4Cs after 48h of Transfection"
fill_color = "#648FFF"


# one graph with Path #1:
rel_conc_data <- rel_conc_p1(qpcr_data, probe, reference, ctrl_sample)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color) + ggtitle(paste0(probe, plot_title))
rel_conc_bar

# one graph with Path #2:
rel_conc_data <- rel_conc_p2(qpcr_data, probe, reference, ctrl_sample)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color) + ggtitle(paste0(probe, plot_title))
rel_conc_bar

