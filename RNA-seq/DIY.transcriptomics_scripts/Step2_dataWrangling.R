# Introduction to the Step 2 script ----
# Now that you've read your transcript-level or gene-level data into R, you're ready to begin working with your data.

# Step 2 Learning Objectives:
# 1 - Filter and normalize your data
# 2 - use ggplot2 to visualize the impact of filtering and normalization on our data.
# 3 - understand why gene expression data is messy, and how to make it 'tidy'

# Notes:
# recall that your abundance data are RPM, while the counts are read counts mapping to each gene or transcript
# do bioinformatics (filtering, normalization) with RNA-seq data
# also need to do plotting (data science)
# tidyverse: a composite universe of 9 different packages! (dplyr, tidyr, ggplot2, etc.)

# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # lets us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure

# Examine your data up to this point ----
# any list object can be examined using $, like a data frame
myTPM <- Txi_gene$abundance # transcripts per million
myCounts <- Txi_gene$counts
colSums(myTPM) # proportion of a million - should be a million
colSums(myCounts) # my counts (total reads mapping to reference per sample)

# capture sample labels from the study design file that you worked with and saved as 'targets' in step 1
targets # study design file is a tibble
sampleLabels <- targets$sample # grabs sample names (good for throwing on a graph later on)

# Generate summary stats for your data ----
# 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
# then we use the 'rowSds', 'rowMeans' and 'rowMedians' functions from the matrixStats package
# calculated new columns based on existing columns
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM)) 

# look at what you created
head(myTPM.stats)

# Create your first plot using ggplot2 ----
# produce a scatter plot of the transformed data
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=25, size=3)
# Experiment with point shape and size in the plot above
# Experiment with other plot types (e.g. 'geom_hex' instead of 'geom_point')
# Add a theme to your ggplot code above.  Try 'theme_bw()'
# How would these graphs change if you log2 converted the data?

# Let's expand on the plot above a bit more and take a look at each 'layer' of the ggplot code
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=16, size=2) + # shape = number, correlates to different symbols
  geom_smooth(method=lm) +
  geom_hex(show.legend = FALSE) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  theme_classic() + # themes are complex - lots of arguments in complete ggplot themes
  theme_dark() + 
  theme_bw()

# Make a DGElist from your counts, and plot ----
# DGEList of raw counts
myDGEList <- DGEList(myCounts)
# take a look at the DGEList object 
myDGEList
#DEGList objects are a good R data file to consider saving to you working directory
save(myDGEList, file = "myDGEList")
#Saved DGEList objects can be easily shared and loaded into an R environment
load(file = "myDGEList")

# use the 'cpm' function from EdgeR to get counts per million
cpm <- edgeR::cpm(myDGEList) # counts per million
colSums(cpm)
log2.cpm <- edgeR::cpm(myDGEList, log=TRUE) # log=TRUE converts to log2 automatically, can do log10

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") # drops rownames (thinks it should be another column in data) so we get it back
log2.cpm.df
# add your sample names to this dataframe (we lost these when we read our data in with tximport) top row!
colnames(log2.cpm.df) <- c("geneID", sampleLabels) # see how sampleLabels is a useful variable?
# Tidy data philosophy:
# each variable forms a column
# each observation forms a row
# each type of observational unit forms a table
# e.g. if every variable formed a column (treatment 1, treatment 2, etc.) then the data frame is very long
# pivot longer so there is a "treatment" column with the treatment numbers underneath as the rows
# use the tidy package to 'pivot' your dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)
# tidy data is ideal for plotting

# let's look at the impact of pivoting the data
log2.cpm.df.pivot

# not it is easy to plot this pivoted data
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() + 
  coord_flip()
# just the distribution of data across samples
# unusual if disrtibutions looks different between samples
  
# what do you think of the distribution of this data?
# Try using coord_flip() at the end of the ggplot code

# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(myDGEList$counts==0)==10) # how many samples had zero counts 
# FALSE  TRUE 
# 29278  6014 
# breaking down the line above is a little tricky.  Let's try:
# 1st - 'myDGEList$counts==0' returns a new 'logical matrix' where each observation (gene) is evaluated (TRUE/FALSE) for each variable (sample) as to whether it has zero counts
# 2nd - passing this logical matrix to 'rowsums' allows you to sum the total number of times an observation was 'TRUE' across all samples
# 3rd - adding the '==10' is a simple way of asking how many of the rowsums equaled 10. In other words, how many genes had 0 counts (TRUE) for all samples in our dataset
# 4th - passing all this to the 'table' function just provides a handy way to summarize the large logical produced in the previous step

# now set some cut-off to get rid of genes/transcripts with low counts
# again using rowSums to tally up the 'TRUE' results of a simple evaluation
# how many genes had more than 1 CPM (TRUE) in at least 3 samples

# The line below is important! This is where the filtering starts
# Be sure to adjust this cutoff for the number of samples in the smallest group of comparison.
keepers <- rowSums(cpm>1)>=5 # how many genes had counts per million that were greater than 1 and greater than or equal to 5 samples
# 5 samples bc theres 5 healthy and 5 sick patients
# now use base R's simple subsetting method to filter your DGEList based on the logical produced above
# [rows,columns]
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered) # dimensions

log2.cpm.filtered <- edgeR::cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)


p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Normalize your data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
# take a look at this new DGEList object...how has it changed?

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
# pivot this NORMALIZED data, just as you did earlier
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = -1, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)


p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# what if we wanted to put all three violin plots together?
# go back and assign each plot to a variable (rather than printing to the plots viewer)
# here we assigned the last 3 plots to p1, p2 and p3
# we'll use the 'plot_grid' function from the cowplot package to put these together in a figure
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
print("Step 2 complete!")

# the essentials ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = -1, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = -1, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

