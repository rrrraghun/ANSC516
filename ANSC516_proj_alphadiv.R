setwd("C:/Users/rajsr/OneDrive/Documents/ANSC_516/ANSC516-proj")
list.files()

# Modified from the original online version available at 
# http://rpubs.com/dillmcfarlan/R_microbiotaSOP

# and Tutorial: Integrating QIIME2 and R for data visualization 
# and analysis using qiime2R
# https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

# Data manipulation
## Load Packages

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
sessionInfo()
library(tidyverse)
library(qiime2R)
library(ggpubr)
library(ggplot2)

##Load Data

meta<-read_q2metadata("proj-metadata.tsv")
str(meta) #structure of metadata
colnames(meta)[2] <- "Bread" #changes column headings[#colno] 
colnames(meta)[3] <- "process"
colnames(meta)[5] <- "Ethnicity"
colnames(meta)[17] <- "PersonNumber1"
str(meta)

evenness = read_qza("evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged
colnames(shannon)[2] <- "shannon"

faith_pd = read_qza("faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
faith_pd = faith_pd[,-1]
colnames(faith_pd)[1] <- "SampleID"
colnames(faith_pd)[2] <- "faith_pd"

## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"

# You always need to check the data types in your tables to make 
# sure they are what you want. We will now change some data types 
# in the meta now

str(meta)
#observed_features$observed_features_num <- lapply(observed_features$observed_features, as.numeric)
#observed_features$observed_features <- as.numeric(observed_features$observed_features)
str(observed_features)

###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID
#meta = meta[,-1]
str(meta)

#Plots
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(meta$shannon, title = "Shannon")
ggqqplot(meta$faith_pd, title = "Faith PD")
ggqqplot(meta$pielou_e, title = "Evenness")
ggqqplot(meta$observed_features, title = "Observed Features")

# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta$shannon)
shapiro.test(meta$faith_pd)
shapiro.test(meta$pielou_e)
shapiro.test(meta$observed_features)

# The null hypothesis of these tests is that “sample distribution 
# is normal”. If the test is significant, the distribution is non-normal.

#Run the ANOVA and save it as an object - for process
aov.evenness.process = aov(pielou_evenness ~ process, data=meta)
aov.faith_pd.process = aov(faith_pd ~ process, data=meta)
aov.shannon.process = aov(shannon ~ process, data=meta)
aov.observed_features.process = aov(observed_features ~ process, data=meta)

#Call for the summary of that ANOVA, which will include P-values - for process
summary(aov.evenness.process)
summary(aov.faith_pd.process)
summary(aov.shannon.process)
summary(aov.observed_features.process)

#To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

TukeyHSD(aov.evenness.process)
TukeyHSD(aov.faith_pd.process)
TukeyHSD(aov.shannon.process)
TukeyHSD(aov.observed_features.process)

levels(meta$process)
#Re-order the groups because the default is alphabetical order
meta$process.ord = factor(meta$process, c("boiled", "extruded", "sour", "yeast", "NoYeast", "none"))
#process.ord is the name of the column in metadata file which is rearranged
levels(meta$process.ord)

#Plot
boxplot(pielou_evenness ~ process.ord, data=meta, ylab="Pielou evenness")

evenness_boxplot <- ggplot(meta, aes(process.ord, pielou_evenness)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/evenness-process.png", evenness_boxplot, height = 5, width = 3)

evenness_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(process.ord) %>%   # the grouping variable
  summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness = n(),  # calculates the sample size per group
            se_evenness = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

evenness_se <- ggplot(evenness_summary, aes(process.ord, mean_evenness, fill = process.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's evenness  ± s.e.", x = "") 
ggsave("output/evenness_se-process.png", evenness_se, height = 4, width = 3)

############FAITH_PD##################

boxplot(faith_pd ~ process.ord, data=meta, ylab="Faith's Phylogenetic Diversity")

faith_pd_boxplot <- ggplot(meta, aes(process.ord, faith_pd)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/faith_pd-process.png", faith_pd_boxplot, height = 4, width = 3)

faith_pd_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(process.ord) %>%   # the grouping variable
  summarise(mean_faith_pd = mean(faith_pd),  # calculates the mean of each group
            sd_faith_pd = sd(faith_pd), # calculates the standard deviation of each group
            n_faith_pd = n(),  # calculates the sample size per group
            se_faith_pd = sd(faith_pd)/sqrt(n())) # calculates the standard error of each group

faith_pd_se <- ggplot(faith_pd_summary, aes(process.ord, mean_faith_pd, fill = process.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_faith_pd - se_faith_pd, ymax = mean_faith_pd + se_faith_pd), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith's Phylogenetic Diversity  ± s.e.", x = "") 
ggsave("output/faith_pd_se-process.png", faith_pd_se, height = 4, width = 3)

#########SHANNON###############

boxplot(shannon ~ process.ord, data=meta, ylab="Shannon Diversity")

shannon_boxplot <- ggplot(meta, aes(process.ord, shannon)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/shannon-process.png", shannon_boxplot, height = 3, width = 3)

shannon_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(process.ord) %>%   # the grouping variable
  summarise(mean_shannon = mean(shannon),  # calculates the mean of each group
            sd_shannon = sd(shannon), # calculates the standard deviation of each group
            n_shannon = n(),  # calculates the sample size per group
            se_shannon = sd(shannon)/sqrt(n())) # calculates the standard error of each group

shannon_se <- ggplot(shannon_summary, aes(process.ord, mean_shannon, fill = process.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon Diversity  ± s.e.", x = "") 
ggsave("output/shannon_se-process.png", shannon_se, height = 4, width = 3)

#############OBSERVED_FEATURES#############

boxplot(observed_features ~ process.ord, data=meta, ylab="Observed Features")

observed_features_boxplot <- ggplot(meta, aes(process.ord, observed_features)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/observed_features-process.png", observed_features_boxplot, height = 3, width = 3)

# Now, the above graph is kind of not correct. Our test and our graphic do not exactly match. ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. Unfortunately plotting the average and standard deviation is a little complicated.

observed_features_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(process.ord) %>%   # the grouping variable
  summarise(mean_observed_features = mean(observed_features),  # calculates the mean of each group
            sd_observed_features = sd(observed_features), # calculates the standard deviation of each group
            n_observed_features = n(),  # calculates the sample size per group
            se_observed_features = sd(observed_features)/sqrt(n())) # calculates the standard error of each group

observed_features_se <- ggplot(observed_features_summary, aes(process.ord, mean_observed_features, fill = process.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_observed_features - se_observed_features, ymax = mean_observed_features + se_observed_features), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features  ± s.e.", x = "") 
ggsave("output/observed_features_se-process.png", observed_features_se, height = 4, width = 3)

## **Non-normally distributed metrics**

# Since process is categorical, we use Kruskal-Wallis (non-parametric equivalent of 
# ANOVA). If we have only two levels, we would run Wilcoxon rank sum 
# test (non-parametric equivalent of t-test)

#kruskal.test(pielou_evenness ~ process.ord, data=meta)
#kruskal.test(faith_pd ~ process.ord, data=meta)
#kruskal.test(shannon ~ process.ord, data=meta)
#kruskal.test(observed_features ~ process.ord, data=meta)

# We can test pairwise within the age groups with Wilcoxon Rank Sum 
# Tests. This test has a slightly different syntax than our other tests

#pairwise.wilcox.test(meta$faith_pd, meta$process.ord, p.adjust.method="BH")
#pairwise.wilcox.test(meta$pielou_evenness, meta$process.ord, p.adjust.method="BH")
#pairwise.wilcox.test(meta$shannon, meta$process.ord, p.adjust.method="BH")
#pairwise.wilcox.test(meta$observed_features, meta$process.ord, p.adjust.method="BH")