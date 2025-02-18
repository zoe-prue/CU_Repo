# summer_student_qPCR_program
# ZMP
# 02/12/2024

# Description --------------------------------
# This is a program designed to quickly interpret qPCR data via the Livak Method 
# and produce graphical representations of each step of the analysis
# make sure to set your current directory as your working directory!

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
source("summer_student_qPCR_functions.R") 

# Global Variables ---------------------------

# (CHANGE THESE TO SUIT YOUR NEEDS)
# make sure your qPCR_data and your R scripts are in the same folder (directory) 
# and that this is set as your working directory
skip_meta = 23 # number of metadata rows to skip
qpcr_data <- read.csv("~/Desktop/CU_coding/qPCR/raw_data/Exp9_real_4.4.24_20240404_122458_Results_20240404_141702.csv", 
                      skip = skip_meta, 
                      stringsAsFactors = FALSE) # put in path to YOUR data

# There are two scenarios; either you use a housekeeping VIC-TAMRA in the same well, or housekeeping FAM in separate wells
# let's call these Path #1 and #2, respectively

# Single Graph -------------------------------
# here is some simple code to make one plot
# Local variables:
probe = "JAK1" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "pTRIP Luciferase" # sample you're making everything "relative" to
plot_title = " Baseline Expression in U4Cs after 48h of Transfection"
fill_color = "#648FFF"
# and here are some colors to use!
c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")


# one graph with Path #1:
rel_conc_data <- rel_conc_p1(qpcr_data, probe, reference, ctrl_sample)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color, probe) + ggtitle(paste0(probe, plot_title))
rel_conc_bar

# one graph with Path #2:
rel_conc_data <- rel_conc_p2(qpcr_data, probe, reference, ctrl_sample)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color) + ggtitle(paste0(probe, plot_title))
rel_conc_bar

