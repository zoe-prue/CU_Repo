# exp_8_qpcr_program_v2
# ZMP
# 02/15/2024

# Description --------------------------------
# Quickly interpret qPCR data via the Livak Method (ddCq) method
# and produce graphs of each step of the analysis
# make sure to set your current directory as your working directory!

# Initialize ---------------------------------

library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Load Functions -----------------------------

source("my_qPCR_functions_v2.R") 

# Global Variables ---------------------------
skip_meta = 20 # number of metadata rows to skip
qpcr_data <- read.csv("Exp8_02.09.2024_real_20240209_121450_Results_20240209 134841.csv", 
                      skip = skip_meta, 
                      stringsAsFactors = FALSE) # put in path to YOUR data
head(qpcr_data)

# Program Body JAK1 -------------------------------

# tidy_data <- well_positions(qpcr_data) # yields tidy_data organized by well positions
# head(tidy_data)
# plate <- diagram_plate(tidy_data) # yields diagram of tidy_data
# plate

#SET THESE VARIABLES
probe = "JAK1" # probe of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "0/0" # sample you're making everything "relative" to
fill_color = "palevioletred2"

summarized_data <- summarize_data(qpcr_data, probe, reference) # yields summarized_data
test_data <- target_cq(summarized_data, probe) # yields test_data
ref_data <- reference_cq(summarized_data, reference)
delta_cq_data <- delta_cq(test_data, ref_data) # change output to delta_cq_data

# delta_cq_bar <- graph_delta_cq(delta_cq_data, probe, fill_color)
# delta_cq_bar


delta_delta_cq_data <- delta_delta_cq(delta_cq_data, ctrl_sample)
# delta_delta_cq_bar <- graph_delta_delta_cq(delta_delta_cq_data, probe, fill_color)
# delta_delta_cq_bar


rel_conc_data <- rel_conc(delta_delta_cq_data)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color)
rel_conc_bar
rel_conc_bar + ggtitle (paste0(probe, " mRNA in U4Cs after 48h of Transfection"))


### playing around
# Rename the column Sample.x to Sample_x
colnames(rel_conc_data)[which(names(rel_conc_data) == "Sample.x")] <- "Sample_x"

rel_conc_data$Sample_x <- trimws(rel_conc_data$Sample_x)
has_50 <- any(startsWith(rel_conc_data$Sample_x, "50/"))
print(has_50)
rel_conc_data$Sample_x <- as.character(rel_conc_data$Sample_x)

# Create subsets for each category
subset_0 <- subset(rel_conc_data, grepl("^0/", Sample_x))
subset_50 <- subset(rel_conc_data, grepl("^50/", Sample_x))
subset_500 <- subset(rel_conc_data, grepl("^500 ", Sample_x))

# Convert Sample_x to character type to avoid any potential issues
rel_conc_data$Sample_x <- as.character(rel_conc_data$Sample_x)

# Extract the first three characters of each string
first_three_chars <- substr(rel_conc_data$Sample_x, 1, 3)

# Check if any of the strings start with "50/"
has_50 <- any(first_three_chars == "50/")

# Print the result
print(has_50)

# Set up the layout for the plots
par(mfrow=c(3, 1)) # 3 rows, 1 column

# Plot for "0/" category
barplot(subset_0$rel_conc, names.arg = subset_0$Sample_x,
        xlab = "Categories", ylab = "rel_conc",
        main = "rel_conc for Categories starting with 0/",
        col = "skyblue", ylim = c(0, max(subset_0$rel_conc) * 1.2))

# Plot for "50/" category
barplot(subset_50$rel_conc, names.arg = subset_50$Sample_x,
        xlab = "Categories", ylab = "rel_conc",
        main = "rel_conc for Categories starting with 50/",
        col = "skyblue", ylim = c(0, max(subset_50$rel_conc) * 1.2))

# Plot for "500/" category
barplot(subset_500$rel_conc, names.arg = subset_500$Sample_x,
        xlab = "Categories", ylab = "rel_conc",
        main = "rel_conc for Categories starting with 500/",
        col = "skyblue", ylim = c(0, max(subset_500$rel_conc) * 1.2))


###

rel_conc_bar + ggtitle (paste0(probe, " mRNA in U4Cs after 48h of Transfection"))
# adjust to the title you want
# for JAK1 as probe, need to compare all transfected samples to the non-transfected sample, respectively
# for ISGs as probe, need to compare between 100ng and 1000 ng transfection treatments of each type, with WT as ctrl_sample

# Program Body RSAD2 50 ng -------------------------------

# tidy_data <- well_positions(qpcr_data) # yields tidy_data organized by well positions
# head(tidy_data)
# plate <- diagram_plate(tidy_data) # yields diagram of tidy_data
# plate

#SET THESE VARIABLES
probe = "RSAD2" # probe of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 Luc/0" # sample you're making everything "relative" to
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000", 
                        "50 Luc/0", "50 Luc/10", "50 Luc/100", "50 Luc/1000", 
                        "50 WT/0", "50 WT/10", "50 WT/100", "50 WT/1000", 
                        "50 S703I/0", "50 S703I/10", "50 S703I/100", "50 S703I/1000", 
                        "50 R326K/0", "50 R326K/10", "50 R326K/100", "50 R326K/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

summarized_data <- summarize_data(subset_data, probe, reference) # yields summarized_data
test_data <- target_cq(summarized_data, probe) # yields test_data
ref_data <- reference_cq(summarized_data, reference)
delta_cq_data <- delta_cq(test_data, ref_data) # change output to delta_cq_data

delta_cq_bar <- graph_delta_cq(delta_cq_data, probe, fill_color)
delta_cq_bar


delta_delta_cq_data <- delta_delta_cq(delta_cq_data, ctrl_sample)
delta_delta_cq_bar <- graph_delta_delta_cq(delta_delta_cq_data, probe, fill_color)
delta_delta_cq_bar


rel_conc_data <- rel_conc(delta_delta_cq_data)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color)
rel_conc_bar

# Exp 8
rel_conc_bar + ggtitle (paste0(probe, " Baseline Expression in U4Cs after 48h of Transfection")) # adjust to the title you want
# for JAK1 as probe, need to compare all transfected samples to the non-transfected sample, respectively
# for ISGs as probe, need to compare between 100ng and 1000 ng transfection treatments of each type, with WT as ctrl_sample


# Program Body RSAD2 500 ng -------------------------------

# tidy_data <- well_positions(qpcr_data) # yields tidy_data organized by well positions
# head(tidy_data)
# plate <- diagram_plate(tidy_data) # yields diagram of tidy_data
# plate

#SET THESE VARIABLES
probe = "RSAD2" # probe of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 Luc/0" # sample you're making everything "relative" to
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000", 
                        "500 Luc/0", "500 Luc/10", "500 Luc/100", "500 Luc/1000", 
                        "500 WT/0", "500 WT/10", "500 WT/100", "500 WT/1000", 
                        "500 S703I/0", "500 S703I/10", "500 S703I/100", "500 S703I/1000", 
                        "500 R326K/0", "500 R326K/10", "500 R326K/100", "500 R326K/1000")
fill_color = "skyblue3"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

summarized_data <- summarize_data(subset_data, probe, reference) # yields summarized_data
test_data <- target_cq(summarized_data, probe) # yields test_data
ref_data <- reference_cq(summarized_data, reference)
delta_cq_data <- delta_cq(test_data, ref_data) # change output to delta_cq_data

delta_cq_bar <- graph_delta_cq(delta_cq_data, probe, fill_color)
delta_cq_bar


delta_delta_cq_data <- delta_delta_cq(delta_cq_data, ctrl_sample)
delta_delta_cq_bar <- graph_delta_delta_cq(delta_delta_cq_data, probe, fill_color)
delta_delta_cq_bar


rel_conc_data <- rel_conc(delta_delta_cq_data)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color)
rel_conc_bar

# Exp 8
rel_conc_bar + ggtitle (paste0(probe, "Baseline Expression in U4Cs after 48h of Transfection")) # adjust to the title you want
# for JAK1 as probe, need to compare all transfected samples to the non-transfected sample, respectively
# for ISGs as probe, need to compare between 100ng and 1000 ng transfection treatments of each type, with WT as ctrl_sample

# Program Body IFI27 50 ng -------------------------------

# tidy_data <- well_positions(qpcr_data) # yields tidy_data organized by well positions
# head(tidy_data)
# plate <- diagram_plate(tidy_data) # yields diagram of tidy_data
# plate

#SET THESE VARIABLES
probe = "IFI27" # probe of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 Luc/0" # sample you're making everything "relative" to
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000", 
                        "50 Luc/0", "50 Luc/10", "50 Luc/100", "50 Luc/1000", 
                        "50 WT/0", "50 WT/10", "50 WT/100", "50 WT/1000", 
                        "50 S703I/0", "50 S703I/10", "50 S703I/100", "50 S703I/1000", 
                        "50 R326K/0", "50 R326K/10", "50 R326K/100", "50 R326K/1000")

fill_color = "darkseagreen3"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

summarized_data <- summarize_data(subset_data, probe, reference) # yields summarized_data
test_data <- target_cq(summarized_data, probe) # yields test_data
ref_data <- reference_cq(summarized_data, reference)
delta_cq_data <- delta_cq(test_data, ref_data) # change output to delta_cq_data

delta_cq_bar <- graph_delta_cq(delta_cq_data, probe, fill_color)
delta_cq_bar


delta_delta_cq_data <- delta_delta_cq(delta_cq_data, ctrl_sample)
delta_delta_cq_bar <- graph_delta_delta_cq(delta_delta_cq_data, probe, fill_color)
delta_delta_cq_bar


rel_conc_data <- rel_conc(delta_delta_cq_data)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color)
rel_conc_bar

# Exp 8
rel_conc_bar + ggtitle (paste0(probe, " Baseline Expression in U4Cs after 48h of Transfection")) # adjust to the title you want
# for JAK1 as probe, need to compare all transfected samples to the non-transfected sample, respectively
# for ISGs as probe, need to compare between 100ng and 1000 ng transfection treatments of each type, with WT as ctrl_sample

# Program Body IFI27 500 ng -------------------------------

# tidy_data <- well_positions(qpcr_data) # yields tidy_data organized by well positions
# head(tidy_data)
# plate <- diagram_plate(tidy_data) # yields diagram of tidy_data
# plate

#SET THESE VARIABLES
probe = "IFI27" # probe of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 Luc/0" # sample you're making everything "relative" to
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000", 
                        "500 Luc/0", "500 Luc/10", "500 Luc/100", "500 Luc/1000", 
                        "500 WT/0", "500 WT/10", "500 WT/100", "500 WT/1000", 
                        "500 S703I/0", "500 S703I/10", "500 S703I/100", "500 S703I/1000", 
                        "500 R326K/0", "500 R326K/10", "500 R326K/100", "500 R326K/1000")
fill_color = "darkseagreen4"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

summarized_data <- summarize_data(subset_data, probe, reference) # yields summarized_data
test_data <- target_cq(summarized_data, probe) # yields test_data
ref_data <- reference_cq(summarized_data, reference)
delta_cq_data <- delta_cq(test_data, ref_data) # change output to delta_cq_data

delta_cq_bar <- graph_delta_cq(delta_cq_data, probe, fill_color)
delta_cq_bar


delta_delta_cq_data <- delta_delta_cq(delta_cq_data, ctrl_sample)
delta_delta_cq_bar <- graph_delta_delta_cq(delta_delta_cq_data, probe, fill_color)
delta_delta_cq_bar


rel_conc_data <- rel_conc(delta_delta_cq_data)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color)
rel_conc_bar

# Exp 8
rel_conc_bar + ggtitle (paste0(probe, " Baseline Expression in U4Cs after 48h of Transfection")) # adjust to the title you want
# for JAK1 as probe, need to compare all transfected samples to the non-transfected sample, respectively
# for ISGs as probe, need to compare between 100ng and 1000 ng transfection treatments of each type, with WT as ctrl_sample


################## new code


# "0/0", "0/10", "0/100", "0/1000",
# "500 Luc/0", "500 Luc/10", "500 Luc/100", "500 Luc/1000")
# "500 WT/0", "500 WT/10", "500 WT/100", "500 WT/1000",
# "500 S703I/0", "500 S703I/10", "500 S703I/100", "500 S703I/1000",
# "500 R326K/0", "500 R326K/10", "500 R326K/100", "500 R326K/1000"



probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "0/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data1 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_0 <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_0

#####

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 Luc/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 Luc/0", "50 Luc/10", "50 Luc/100", "50 Luc/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data2 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_luc <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_luc

###

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 WT/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 WT/0", "50 WT/10", "50 WT/100", "50 WT/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data3 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_wt <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_wt

####

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 S703I/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 S703I/0", "50 S703I/10", "50 S703I/100", "50 S703I/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data4 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_s703i <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_s703i


####

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 R326K/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 R326K/0", "50 R326K/10", "50 R326K/100", "50 R326K/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data5 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)

###

combined_data <- bind_rows(
  rel_conc_data5 %>% mutate(dataset_id = "Dataset 5"),
  rel_conc_data4 %>% mutate(dataset_id = "Dataset 4"),
  rel_conc_data3 %>% mutate(dataset_id = "Dataset 3"),
  rel_conc_data2 %>% mutate(dataset_id = "Dataset 2"),
  rel_conc_data1 %>% mutate(dataset_id = "Dataset 1")
)

# Calculate means and standard deviations per sample category
summary_data <- combined_data %>%
  group_by(Sample, dataset_id) %>%
  summarise(mean = mean(rel_conc), sd = sd(rel_conc))

# Define custom colors
custom_colors <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000")

# Create a bar plot with error bars
p <- ggplot(summary_data, aes(x = Sample, y = mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           aes(fill = dataset_id), color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.4, position = position_dodge(width = 0.8)) +
  geom_point(data = combined_data, aes(x = Sample, y = rel_conc), 
             position = position_jitter(width = 0.2), 
             color = "black", size = 1) + # Add jittered points
  labs(x = "Sample", y = "RSAD2 Relative Concentration", 
       title = "RSAD2 Cytokine-Stimulated Expression in U4Cs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels if needed
  scale_y_continuous(trans = "log10") + # Set log scale on Y-axis
  guides(fill = FALSE) + # Remove legend for fill color
  scale_fill_manual(values = custom_colors) + # Set custom fill colors for bars
  scale_color_manual(values = "black") # Set black color for points

# Display the plot
print(p)

### IFI27####
probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "0/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data1 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_0 <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_0

#####

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 Luc/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 Luc/0", "50 Luc/10", "50 Luc/100", "50 Luc/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data2 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_luc <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_luc

###

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 WT/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 WT/0", "50 WT/10", "50 WT/100", "50 WT/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data3 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_wt <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_wt

####

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 S703I/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 S703I/0", "50 S703I/10", "50 S703I/100", "50 S703I/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data4 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_s703i <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_s703i


####

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 R326K/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("50 R326K/0", "50 R326K/10", "50 R326K/100", "50 R326K/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data5 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)

###

combined_data <- bind_rows(
  rel_conc_data5 %>% mutate(dataset_id = "Dataset 5"),
  rel_conc_data4 %>% mutate(dataset_id = "Dataset 4"),
  rel_conc_data3 %>% mutate(dataset_id = "Dataset 3"),
  rel_conc_data2 %>% mutate(dataset_id = "Dataset 2"),
  rel_conc_data1 %>% mutate(dataset_id = "Dataset 1")
)

# Calculate means and standard deviations per sample category
summary_data <- combined_data %>%
  group_by(Sample, dataset_id) %>%
  summarise(mean = mean(rel_conc), sd = sd(rel_conc))

# Define custom colors
custom_colors <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000")

# Create a bar plot with error bars
p <- ggplot(summary_data, aes(x = Sample, y = mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           aes(fill = dataset_id), color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.4, position = position_dodge(width = 0.8)) +
  geom_point(data = combined_data, aes(x = Sample, y = rel_conc), 
             position = position_jitter(width = 0.2), 
             color = "black", size = 1) + # Add jittered points
  labs(x = "Sample", y = "IFI27 Relative Concentration", 
       title = "IFI27 Cytokine-Stimulated Expression in U4Cs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels if needed
  scale_y_continuous(trans = "log10") + # Set log scale on Y-axis
  guides(fill = FALSE) + # Remove legend for fill color
  scale_fill_manual(values = custom_colors) + # Set custom fill colors for bars
  scale_color_manual(values = "black") # Set black color for points

# Display the plot
print(p)



### 500 RSAD2#####

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "0/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data1 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_0 <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_0

#####

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 Luc/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 Luc/0", "500 Luc/10", "500 Luc/100", "500 Luc/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data2 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_luc <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_luc

###*

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 WT/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 WT/0", "500 WT/10", "500 WT/100", "500 WT/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data3 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_wt <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_wt

####

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 S703I/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 S703I/0", "500 S703I/10", "500 S703I/100", "500 S703I/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data4 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_s703i <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_s703i


####

probe = "RSAD2" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 R326K/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 R326K/0", "500 R326K/10", "500 R326K/100", "500 R326K/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data5 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)

###

combined_data <- bind_rows(
  rel_conc_data5 %>% mutate(dataset_id = "Dataset 5"),
  rel_conc_data4 %>% mutate(dataset_id = "Dataset 4"),
  rel_conc_data3 %>% mutate(dataset_id = "Dataset 3"),
  rel_conc_data2 %>% mutate(dataset_id = "Dataset 2"),
  rel_conc_data1 %>% mutate(dataset_id = "Dataset 1")
) %>% 
  na.omit()  # Remove rows with NAs

# Calculate means and standard deviations per sample category
summary_data <- combined_data %>%
  group_by(Sample, dataset_id) %>%
  summarise(mean = mean(rel_conc), sd = sd(rel_conc))

# Define custom colors
custom_colors <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000")

# Create a bar plot with error bars
p <- ggplot(summary_data, aes(x = Sample, y = mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           aes(fill = dataset_id), color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.4, position = position_dodge(width = 0.8)) +
  geom_point(data = combined_data, aes(x = Sample, y = rel_conc), 
             position = position_jitter(width = 0.2), 
             color = "black", size = 1) + # Add jittered points
  labs(x = "Sample", y = "RSAD2 Relative Concentration", 
       title = "RSAD2 Cytokine-Stimulated Expression in U4Cs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels if needed
  scale_y_continuous(trans = "log10") + # Set log scale on Y-axis
  guides(fill = FALSE) + # Remove legend for fill color
  scale_fill_manual(values = custom_colors) + # Set custom fill colors for bars
  scale_color_manual(values = "black") # Set black color for points

# Display the plot
print(p)


### 500 IFI27#####

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "0/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("0/0", "0/10", "0/100", "0/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data1 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_0 <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_0

#####

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 Luc/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 Luc/0", "500 Luc/10", "500 Luc/100", "500 Luc/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data2 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_luc <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_luc

###*

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 WT/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 WT/0", "500 WT/10", "500 WT/100", "500 WT/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data3 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_wt <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_wt

####

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 S703I/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 S703I/0", "500 S703I/10", "500 S703I/100", "500 S703I/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data4 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)
# fifty_rsad2_s703i <- graph_rel_conc(rel_conc_data, fill_color)
# fifty_rsad2_s703i


####

probe = "IFI27" # list of probes of interest
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "500 R326K/0" # sample you're making everything "relative" to
plot_title = " Cytokine Stimulated Expression in U4Cs after 48h of Transfection"
samples_to_include <- c("500 R326K/0", "500 R326K/10", "500 R326K/100", "500 R326K/1000")
fill_color = "skyblue1"

# Subset the data frame using dplyr
subset_data <- qpcr_data %>%
  filter(Sample %in% samples_to_include)

# one graph with Path #2:
rel_conc_data5 <- rel_conc_p2(qpcr_data=subset_data, probe, reference, ctrl_sample)

###

combined_data <- bind_rows(
  rel_conc_data5 %>% mutate(dataset_id = "Dataset 5"),
  rel_conc_data4 %>% mutate(dataset_id = "Dataset 4"),
  rel_conc_data3 %>% mutate(dataset_id = "Dataset 3"),
  rel_conc_data2 %>% mutate(dataset_id = "Dataset 2"),
  rel_conc_data1 %>% mutate(dataset_id = "Dataset 1")
) %>% 
  na.omit()  # Remove rows with NAs

# Calculate means and standard deviations per sample category
summary_data <- combined_data %>%
  group_by(Sample, dataset_id) %>%
  summarise(mean = mean(rel_conc), sd = sd(rel_conc))

# Define custom colors
custom_colors <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000")

# Create a bar plot with error bars
p <- ggplot(summary_data, aes(x = Sample, y = mean)) + # create barplot
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           aes(fill = dataset_id), color = "black") + # format the barplot
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), # add error bars with standard deviation
                width = 0.4, position = position_dodge(width = 0.8)) +
  geom_point(data = combined_data, aes(x = Sample, y = rel_conc), # add the individual data points
             position = position_jitter(width = 0.2),# make sure the points don't overlap
             color = "black", size = 1) + 
  labs(x = "Sample", y = "IFI27 Relative Concentration", # add labels
       title = "IFI27 Cytokine-Stimulated Expression in U4Cs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) + # Rotate x-axis labels if needed
  scale_y_continuous(trans = "log10") + # Set log10 scale on Y-axis
  guides(fill = FALSE) + # Remove legend for fill color
  scale_fill_manual(values = custom_colors) + # Set custom fill colors for bars
  scale_color_manual(values = "black") # Set black color for points
# remove graph background

# Display the plot
print(p)





### baseline expression with significance tests


probe_list = c("JAK1", "RSAD2", "IFI27") # list of probes of interest
probe = "IFI27"
reference = "GAPDH" # probe that you're "normalizing" to
ctrl_sample = "50 Luc/0" # sample you're making everything "relative" to
my_cols <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#0072B2", "#FFB000") # should be sample length as probe_list
file_name <- "Exp8_plots_4" # name of file you want to make for plots
plot_title = " Baseline Expression in U4Cs after 48h of Transfection"
my_plots <- list() # make an empty list to put plots into

# new function to aciheive this
# group by stim level, compare stim levels with t-test - except WT-S703I relationship
graph_rel_conc <- function(rel_conc_data, fill_color, probe) {
  
  # Extract the numeric portion from the Sample column
  rel_conc_data <- rel_conc_data %>%
    mutate(Sample_Group = case_when(
      grepl("/0", Sample) ~ "/0",
      grepl("/10", Sample) ~ "/10",
      grepl("/100", Sample) ~ "/100",
      grepl("/1000", Sample) ~ "/1000",
      TRUE ~ "Other"
    ))
  
  # Calculate mean and standard deviation for each sample group
  summary_data <- rel_conc_data %>%
    group_by(Sample_Group) %>%
    summarise(mean_rel_conc = mean(rel_conc, na.rm = TRUE),
              sd_rel_conc = sd(rel_conc, na.rm = TRUE))
  
  # Plotting
  rel_conc_bar <- ggplot(summary_data, aes(x = Sample_Group, y = mean_rel_conc, fill = Sample_Group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean_rel_conc - sd_rel_conc, ymax = mean_rel_conc + sd_rel_conc),
                  position = position_dodge(width = 0.9), width = 0.2) +
    theme_classic() +
    scale_y_continuous(trans = "log10") +
    ylab(paste0("Relative Concentrations of ", probe)) +
    xlab("Sample Type") +
    ggtitle("Relative Concentration Distribution by Sample Type")
  
  return(rel_conc_bar)
}

rel_conc_data <- rel_conc_p2(qpcr_data, probe, reference, ctrl_sample)
rel_conc_bar <- graph_rel_conc(rel_conc_data, fill_color, probe) + ggtitle(paste0("IFI27", plot_title))
rel_conc_bar


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

