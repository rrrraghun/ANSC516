# This is a demo for running the co-occurrence analysis 

#make sure you have these libraries
library(Hmisc)
library(plyr)
library(reshape2)
library(qiime2R)
library(ggplot2)
#library(igraph)
#library(fdrtool)

setwd("C:/Users/rajsr/OneDrive/Documents/ANSC_516/ANSC516-proj")

ASVs <- read_qza("table.qza")
ASV_table <- as.data.frame(ASVs$data)

#####################################################################
##I added in this chuck. What does it do and why?

ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table)) #pastes ASVnos from ASVs onto ASV_table
ASV_table$ASVstring <- rownames(ASV_table) #assign each ASV codename to ASVnos - asv string refers to the long names
rownames(ASV_table) <- ASV_table$ASVnos #changes row names of ASV_table to ASVnos
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #combine ASV rownames with ASVnos into new file ASV_key
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)] #remove ASV rows and ASVnos from ASV_table
######################################################################

dataset <- as.data.frame(t(ASV_table))


# we are going to create a network per treatment
head(dataset[,1:10]) #are we transposing?
metadata<-read_q2metadata("proj-metadata.tsv")
str(metadata)
colnames(metadata)[3] = "process"

dataset <- merge(metadata, dataset, by.x = "SampleID", by.y = 0)
treatments<-as.vector(unique(dataset$process)) #need to choose column, i.e. how to subset
datasetn<-dataset
datasetn[datasetn==0]<-NA #any cells with 0 values are renamed NA

#i = 1

summary(metadata$process)

my_column <- "process"
n1 <- 6
n2 <- 6
n3 <- 6
n4 <- 6
n5 <- 6
n6 <- 6

# above lines indicate that ASVs must be in at least 6 samples to be co occurring

num_metadata_columns <- 21

q_cutoff <- 0.05

final_results<-data.frame()

for(i in 1:length(treatments)){
  #subset the data for a particular treatment YOU MUST ENTER THE HEADER OF THE COLUMN THAT HAS THE DIFFERENT TREATMENTS IN THIS CASE “Foaming_Status”
  print(paste("reading ",treatments[i],sep=""))
  temp<-subset(dataset, get(my_column)==treatments[i])
  tempn<-subset(datasetn, get(my_column)==treatments[i])
  print(paste("finished reading ",treatments[i],sep=""))
  results<-rcorr(as.matrix(temp[,-c(1:num_metadata_columns)]),type="spearman")
  resultsn<-rcorr(as.matrix(tempn[,-c(1:num_metadata_columns)]),type="spearman")
  
  #make two seperate objects for p-value and correlation coefficients
  rhos<-results$r
  ps<-results$P
  ns<-resultsn$n
  # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
  ps_melt<-na.omit(melt(ps))
  #creating a qvalue based on FDR
  ps_melt$qval<-p.adjust(ps_melt$value, method = "BH")
  #making column names more relevant
  
  names(ps_melt)[3]<-"pval"
  # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
  ps_sub<-subset(ps_melt, qval < q_cutoff)
  
  # now melting the rhos, note the similarity between ps_melt and rhos_melt
  rhos_melt<-na.omit(melt(rhos))
  names(rhos_melt)[3]<-"rho"
  
  # now melting the ns
  ns_melt<-(melt(ns))
  names(ns_melt)[3]<-"n"
  
  #merging together and remove negative rhos
  merged<-merge(ps_sub,rhos_melt,by=c("Var1","Var2"))
  if (treatments[i]==treatments[1]) {
    merged<-merge(merged,subset(ns_melt, n > n1),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[2]) {
    merged<-merge(merged,subset(ns_melt, n > n2),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[3]) {
    merged<-merge(merged,subset(ns_melt, n > n3),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[4]) {
    merged<-merge(merged,subset(ns_melt, n > n4),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[5]) {
    merged<-merge(merged,subset(ns_melt, n > n5),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[6]) {
    merged<-merge(merged,subset(ns_melt, n > n6),by=c("Var1","Var2"))
  }   else
    print("Somethings wrong with your treatment designations. Please Check!!")
  
  if (nrow(merged) > 0) {
    merged$trt<-treatments[i]
    final_results<-rbind(final_results, merged)
  }   else {
    print("no correlations for this variable")
  }
  
  print(paste("finished ",treatments[i],sep=""))
}

strong_results<-subset(final_results, rho >= 0.9)
strong_results<-subset(final_results, rho >= 0.7)
strong_results<-subset(final_results, rho >= 0.75)
colnames

###############################################################
# If you want to see the correlation scatterplot of 
# two significant ASVs
###############################################################

gut_ASVs<-subset(dataset, get(my_column)==treatments[1])

colnames(gut_ASVs[1:10])
head(final_results)
ggplot(gut_ASVs, aes(x = ASV3, y = ASV55)) +
  geom_point()

###############################################################
# If you want to see the the taxonomic assignment of  
# these significant ASVs
###############################################################

taxonomy<-read_qza("taxonomy.qza")
head(taxonomy$data)

tax.clean<-parse_taxonomy(taxonomy$data)
head(tax.clean)

#All this is OK except that in future use of the taxonomy table, 
#these ASVs will be ignored because they are not classified. Why 
#are ASVs not classified? Its because there is not a close enough 
#match in the database. Just because there is not a good match in 
#the database does not mean they don’t exist, so I wanted to make 
#sure this data was not lost. So in my new code, 
#I make it so that ASVs that are unclassified at any level are 
#classified as the lowest taxonomic level for which there is a 
#classification.
#Next, all these `NA` classifications with the last level that was 
#classified

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("unclassified_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("unclassified_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("unclassified_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("unclassified_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("unclassified_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("unclassified_",tax.clean$Genus[i], sep = "_")
  }
}

strong_results_taxa <- merge(strong_results, ASVkey, by.x = "Var1", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, ASVkey, by.x = "Var2", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.x", by.y = 0)
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.y", by.y = 0)


write.csv(strong_results_taxa, "output/proj-strong-results-taxa.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="yeast"), "output/proj-strong-results-taxa-yeast.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="boiled"), "output/proj-strong-results-taxa-boiled.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="none"), "output/proj-strong-results-taxa-none.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="sour"), "output/proj-strong-results-taxa-sour.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="NoYeast"), "output/proj-strong-results-taxa-NoYeast.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="extruded"), "output/proj-strong-results-taxa-extruded.csv", row.names = F)

############################################################
# Now you can import these into cytoscape to visualize the network
#Open Cytoscape
#Import Network from File system
#Anything with a .x indicate as a "source node attribute"
#Anything with a .y indicate as a "target node attribute"
#"Var2" indicate as "Target Node"
#"Var1" indicate as "Source Node"
#the rest leave as "edge attribute"
#Edit>Remove Duplicated Edges, click "Ignore edge direction"
#To manually assign aesthetics, use "Discrete Mapping"
#######################################################################

