
#Load the packages
installed.packages()
sessionInfo()

install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")

#remove.packages("rlang")
#install.packages("rlang")

library(tidyverse)
library(vegan)
library(qiime2R)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#proj-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza

##############################################

###Set your working directory
setwd("C:/Users/rajsr/OneDrive/Documents/ANSC_516/ANSC516-proj")

list.files()

if(!dir.exists("output"))
  dir.create("output")

#How to load a file into R
metadata2 <- read.delim("proj-metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata2 <- metadata2[-1,]

#Now the qiime2R method
metadata<-read_q2metadata("proj-metadata.tsv")
str(metadata)
levels(metadata$`process`)
colnames(metadata)[2] <- "Bread"
colnames(metadata)[3] <- "process"
colnames(metadata)[5] <- "Ethnicity"
str(metadata)

row.names(metadata) <- metadata[,1]
metadata.without.col1 <- metadata[,-1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("bray_curtis_pcoa_results.qza")
wUF <- read_qza("weighted_unifrac_pcoa_results.qza")
jaccard_PCoA<-read_qza("jaccard_pcoa_results.qza")
uwUF <- read_qza("unweighted_unifrac_pcoa_results.qza")

body_colors <- c("Red", "Blue", "Green", "Gray", "Yellow", "Pink")

#levels(metadata$PersonNumber1)
#Re-order the groups because the default is alphabetical order
#metadata$process.ord = factor(metadata$process, c("no1", "no2", "no3", "no4", "no5", "no6", "no7","no8", "no9", "no10"))
#process.ord is the name of the column in metadata file which is rearranged
#levels(metadata$process.ord)

#####BRAYCURTIS#############

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))
#rlang::last_trace()
#rlang::last_trace(drop = FALSE)

# Now we are going to make an ordination plot
ggplot(data = bc_meta, aes(x=PC1, y=PC2, color=process)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=c("Red", "Blue", "Green", "Gray", "Yellow", "Pink"), name = "process")

# Now we are going to make our code a little more re-usable
my_column <- "process"
#my_column <- "Bread"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + 
  #geom_point(aes(shape= Bread), size = 3) +
  theme_q2r() +
  #facet_grid(~PersonNumber1) +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-basic-", my_column,".tiff"), height=5, width=25, device="tiff") # save a PDF 3 inches by 4 inches

#rlang::last_trace()
#rlang::last_trace(drop = FALSE)

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "process"

################WEIGHTED UNIFRAC##############

Wuni_PCoA<-read_qza("weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)
colnames(centroids)[1] <- "process"

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + 
  #geom_point(aes(shape= Bread), size = 3)
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #facet_grid(~PersonNumber1)
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Process")
ggsave(paste0("output/Wuni-ellipse-Ethnicity", my_column,".png"), height=3, width=4.5, device="png") # save a PDF 3 inches by 4 inches

###############UNWEIGHTED UNIFRAC###########

uwUF_PCoA<-read_qza("unweighted_unifrac_pcoa_results.qza")

uwUF_meta <- uwUF_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)
colnames(centroids)[1] <- "process"

ggplot(uwUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + 
  #geom_point(aes(shape= Bread), size = 3) + 
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*uwUF_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*uwUF_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Process")
ggsave(paste0("output/uwUF-ellipse-process", my_column,".png"), height=3, width=4.5, device="png") # save a PDF 3 inches by 4 inches

#######JACCARD##############

jaccard_PCoA<-read_qza("weighted_unifrac_pcoa_results.qza")

jaccard_meta <- jaccard_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),jaccard_meta,mean)
colnames(centroids)[1] <- "process"

ggplot(jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + 
  #geom_point(aes(shape= Bread), size = 3) + 
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Process")
ggsave(paste0("output/jaccard-ellipse-process", my_column,".png"), height=3, width=4.5, device="png") # save a PDF 3 inches by 4 inches

bc_dist_mat<-read_qza("bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(bc_dm ~ process*PersonNumber1, data = metadata_sub)

write.table(PERMANOVA_out,"output/process_bc_Adonis_overall.csv",sep=",", row.names = TRUE) 

Wuni_dist_mat<-read_qza("weighted_unifrac_distance_matrix.qza")
Wuni_dm <- as.matrix(Wuni_dist_mat$data) 
rownames(Wuni_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(Wuni_dm),metadata$SampleID),]
rownames(Wuni_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(Wuni_dm ~ process*PersonNumber1, data = metadata_sub)
write.table(PERMANOVA_out,"output/process_Wuni_Adonis_overall.csv",sep=",", row.names = TRUE)

uwUF_dist_mat<-read_qza("unweighted_unifrac_distance_matrix.qza")
uwUF_dm <- as.matrix(uwUF_dist_mat$data) 
rownames(uwUF_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(uwUF_dm),metadata$SampleID),]
rownames(uwUF_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(uwUF_dm ~ process*PersonNumber1, data = metadata_sub)
write.table(PERMANOVA_out,"output/process_uwUF_Adonis_overall.csv",sep=",", row.names = TRUE)

uwUF_dist_mat<-read_qza("unweighted_unifrac_distance_matrix.qza")
uwUF_dm <- as.matrix(uwUF_dist_mat$data) 
rownames(uwUF_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(uwUF_dm),metadata$SampleID),]
rownames(uwUF_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(uwUF_dm ~ process*PersonNumber1, data = metadata_sub)
write.table(PERMANOVA_out,"output/process_uwUF_Adonis_overall.csv",sep=",", row.names = TRUE)

jaccard_dist_mat<-read_qza("jaccard_distance_matrix.qza")
jaccard_dm <- as.matrix(jaccard_dist_mat$data) 
rownames(jaccard_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(jaccard_dm),metadata$SampleID),]
rownames(jaccard_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(jaccard_dm ~ process*PersonNumber1, data = metadata_sub)
write.table(PERMANOVA_out,"output/process_jaccard_Adonis_overall.csv",sep=",", row.names = TRUE)

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

process_Pair <- pairwise.adonis2(bc_dm ~ process, data = metadata_sub)
write.table(process_Pair,"output/process_bc_Adonis_pairwise.csv",sep=",", row.names = TRUE) 

#####

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

process_Pair <- pairwise.adonis2(Wuni_dm ~ process, data = metadata_sub)
write.table(process_Pair,"output/process_Wuni_Adonis_pairwise.csv",sep=",", row.names = TRUE)

####
pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

process_Pair <- pairwise.adonis2(uwUF_dm ~ process, data = metadata_sub)
write.table(process_Pair,"output/process_uwUF_Adonis_pairwise.csv",sep=",", row.names = TRUE)

####
pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

process_Pair <- pairwise.adonis2(jaccard_dm ~ process, data = metadata_sub)
write.table(process_Pair,"output/process_jaccard_Adonis_pairwise.csv",sep=",", row.names = TRUE)

