# summer_student_qPCR_functions_v2
# ZMP
# 02/12/2024

##################################################
# FUNCTION: rel_conc_p1
# packages: dplyr
# purpose: apply delta delta cq method to data, 
# where reference vic is in the same wells as samples
# input: qpcr_data
# output: rel_conc_data
# ------------------------------------------------
rel_conc_p1 <- function(qpcr_data, probe, reference, ctrl_sample) {
  
  # Check if the dataset has only the column "Well" and rename it to "Well.Position" if true
  if ("Well" %in% colnames(qpcr_data) & !("Well.Position" %in% colnames(qpcr_data))) {
    colnames(qpcr_data)[colnames(qpcr_data) == "Well"] <- "Well.Position"
  }
  
  # Data pre-processing - averaging technical replicates by Well.Position, Sample, and Target
  qpcr_data$Cq <- as.numeric(qpcr_data$Cq)
  avg_cq_values <- qpcr_data %>%
    group_by(Well.Position, Sample, Target) %>%
    summarize(Cq = mean(Cq, na.rm = TRUE), .groups = 'drop')
  
  # Filter and select relevant columns
  tidy_qpcr_data <- avg_cq_values %>%
    filter(!is.na(Target) & Target != "") %>%
    select(Well.Position, Sample, Target, Cq)
  
  # Create data tables for the target and reference genes
  test_data <- tidy_qpcr_data %>%
    filter(Target == probe)
  
  ref_data <- tidy_qpcr_data %>%
    filter(Target == reference) %>%
    rename(ref_target = Target, ref_Cq = Cq)
  
  # Calculate average reference Cq for each well
  avg_ref_Cq <- ref_data %>%
    group_by(Well.Position) %>%
    summarize(avg_ref_Cq = mean(ref_Cq, na.rm = TRUE), .groups = 'drop')
  
  # Calculate delta Cq (ΔCq) values for each well
  delta_cq_data <- test_data %>%
    left_join(avg_ref_Cq, by = "Well.Position") %>%
    mutate(delta_Cq = Cq - avg_ref_Cq)
  
  # Calculate mean ΔCq for the control sample
  mean_ctrl_sample_delta <- mean(delta_cq_data$delta_Cq[delta_cq_data$Sample == ctrl_sample], 
                                 na.rm = TRUE)
  
  # Calculate delta-delta Cq (ΔΔCq) values
  delta_delta_cq_data <- delta_cq_data %>%
    mutate(delta_delta_cq = delta_Cq - mean_ctrl_sample_delta)
  
  # Calculate relative concentrations
  rel_conc_data <- delta_delta_cq_data %>%
    mutate(rel_conc = 2^(-delta_delta_cq))
  
  return(rel_conc_data)
}
##################################################

##################################################
# FUNCTION: rel_conc_p2
# packages: dplyr
# purpose: apply delta delta cq method to data, 
# where reference fam is in the separate wells as samples
# input: qpcr_data
# output: rel_conc_data
# ------------------------------------------------
rel_conc_p2 <- function(qpcr_data, probe, reference, ctrl_sample) {
  
  # Check if the dataset has only the column "Well" and rename it to "Well.Position" if true
  if ("Well" %in% colnames(qpcr_data) & !("Well.Position" %in% colnames(qpcr_data))) {
    colnames(qpcr_data)[colnames(qpcr_data) == "Well"] <- "Well.Position"
  }
  
  # select only columns from qpcr_data that are needed and ensure that values are class numeric
  tidy_qpcr_data <- qpcr_data %>%
    filter(Target != is.na(Target) & Target != "") %>%
    select(Well.Position, Sample, Target, Cq) %>%
    group_by(Well.Position) %>% 
    mutate(Cq = as.numeric(Cq))
  
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
    group_by(Sample.x) %>% 
    mutate(delta_Cq = Cq - ref_Cq) %>%
    mutate(delta_Cq = as.numeric(delta_Cq)) %>%
    # select(-ends_with(".y")) %>%
    # rename(Sample = Sample.x)
  
  # find the mean delta_cq for your ctrl_sample (to which every other sample is "relative")
  mean_ctrl_sample_delta <- mean(delta_cq_data$delta_Cq[delta_cq_data$Sample == ctrl_sample], 
                                 na.rm = TRUE)
  
  # subtract mean delta_cq for your ctrl_sample from all delta_cq values
  delta_delta_cq_data <- delta_cq_data %>%
    mutate(delta_delta_cq = delta_Cq - mean_ctrl_sample_delta)
  
  # find the relative concentrations of the target being tested for every sample
  rel_conc_data <- delta_delta_cq_data %>%
    mutate(rel_conc = 2^(-(delta_delta_cq)))
  
  return(rel_conc_data)
}
##################################################


##################################################
# FUNCTION: graph_rel_conc
# packages: ggplot2
# purpose: plot relative mRNA concentrations
# input: rel_conc_data
# output: rel_conc_bar
# ------------------------------------------------
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
    geom_point(data = rel_conc_data, aes(x = group, y = rel_conc), 
               position = position_jitter(width = 0.2), 
               alpha = 0.8, color = "black", 
               size = 1) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,), 
          axis.title.x = element_blank(), 
          plot.title=element_text(hjust=0.5)) +
    # xlab("JAK1") +
    scale_y_continuous(trans = "log10") +
    ylab(paste0(probe, "/", reference, " Relative Units")) +
    theme(plot.title = element_text(size = 16),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))
  
  return(rel_conc_bar)
}
###########################################

