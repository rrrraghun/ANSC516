!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20                    ##number of requested CPUs
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name qiime_classifier ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /depot/microbiome/data/ANSC516/Rajsri
mkdir proj_slurm_output
cd proj_slurm_output
mkdir raw
mkdir seq
ls
pwd

prefetch -O raw –option-file SraAccList.csv
ls
pwd
cp raw/SRR*/*.sra seq/
cd seq

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2022.8"
module load Qiime/2-2022.8

#step 1, cd into the proper directory. This directory must already exist

cd /depot/microbiome/data/ANSC516/Rajsri/data_download_proj
mkdir raw
mkdir seq
ls
pwd
prefetch -O raw –option-file SraAccList.csv
ls 
pwd
cp raw/SRR*/*.sra seq/
cd seq 

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2022.8"
module load Qiime/2-2022.8

#step 3, Run the code you want

echo "initiating qiime2 moving pictures tutorial pipeline"

echo "downloading necessary files: metadata, sequence files and classifier"

#Get the databases to classify the sequences
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2022.8/common/gg-13-8-99-515-806-nb-classifier.qza"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path moving-pictures-manifest.tsv \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data single-end-demux.qza \
  --o-visualization single-end-demux.qzv

echo "dada2"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 190 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats-dada2.qza\



#Convert .qza to .qzv format

qiime metadata tabulate \
--m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

#Create a .qzv from the output of denoise

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file proj-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

echo "Create a phylogenetic tree. (Phylogenetic method)"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#Calculate diversity metrics. You will need to designate a sampling depth. 

qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth 16254 \
--m-metadata-file proj-metadata.tsv \
--output-dir core-metrics-results

#alpha diversity significance
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
--m-metadata-file proj-metadata.tsv \
--o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/evenness_vector.qza \
--m-metadata-file proj-metadata.tsv \
--o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file proj-metadata.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

#beta diversity significance
qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata.tsv \
--m-metadata-column process \
--o-visualization core-metrics-results/unweighted-unifrac-process-significance.qzv \
--p-pairwise

qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file proj-metadata.tsv \
--m-metadata-column Bread \
--o-visualization core-metrics-results/unweighted-unifrac-Bread-significance.qzv \
--p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column process \
  --o-visualization core-metrics-results/bray-curtis-process-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column Bread \
  --o-visualization core-metrics-results/bray-curtis-Bread-significance.qzv \
  --p-pairwise

#emperor plots
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-custom-axes PersonNumber-1 \
  --o-visualization core-metrics-results/unweighted-unifrac-PersonNumber.qzv


qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-custom-axes PersonNumber-1 \
  --o-visualization core-metrics-results/bray-curtis-emperor-PersonNumber.qzv

#taxabarplots
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file proj-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

#ancom-boiled
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-where "[process]='boiled'" \
  --o-filtered-table boiled-table.qza

qiime composition add-pseudocount \
  --i-table boiled-table.qza \
  --o-composition-table comp-boiled-table.qza

qiime composition ancom \
  --i-table comp-boiled-table.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization ancom-boiled-PersonNumber-1.qzv

qiime taxa collapse \
  --i-table boiled-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table boiled-table-l6.qza

qiime composition add-pseudocount \
  --i-table boiled-table-l6.qza \
  --o-composition-table comp-boiled-table-l6.qza

qiime composition ancom \
  --i-table comp-boiled-table-l6.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization l6-ancom-boiled-PersonNumber-1.qzv

#ancom-extruded
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-where "[process]='extruded'" \
  --o-filtered-table extruded-table.qza

qiime composition add-pseudocount \
  --i-table extruded-table.qza \
  --o-composition-table comp-extruded-table.qza

qiime composition ancom \
  --i-table comp-extruded-table.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization ancom-extruded-PersonNumber-1.qzv

qiime taxa collapse \
  --i-table extruded-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table extruded-table-l6.qza

qiime composition add-pseudocount \
  --i-table extruded-table-l6.qza \
  --o-composition-table comp-extruded-table-l6.qza

qiime composition ancom \
  --i-table comp-extruded-table-l6.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization l6-ancom-extruded-PersonNumber-1.qzv
#ancom-none
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-where "[process]='none'" \
  --o-filtered-table none-table.qza

qiime composition add-pseudocount \
  --i-table none-table.qza \
  --o-composition-table comp-none-table.qza

qiime composition ancom \
  --i-table comp-none-table.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization ancom-none-PersonNumber-1.qzv

qiime taxa collapse \
  --i-table none-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table none-table-l6.qza

qiime composition add-pseudocount \
  --i-table none-table-l6.qza \
  --o-composition-table comp-none-table-l6.qza

qiime composition ancom \
  --i-table comp-none-table-l6.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization l6-ancom-none-PersonNumber-1.qzv

#ancom-NoYeast
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-where "[process]='NoYeast'" \
  --o-filtered-table NoYeast-table.qza

qiime composition add-pseudocount \
  --i-table NoYeast-table.qza \
  --o-composition-table comp-NoYeast-table.qza

qiime composition ancom \
  --i-table comp-NoYeast-table.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization ancom-NoYeast-PersonNumber-1.qzv

qiime taxa collapse \
  --i-table NoYeast-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table NoYeast-table-l6.qza

qiime composition add-pseudocount \
  --i-table NoYeast-table-l6.qza \
  --o-composition-table comp-NoYeast-table-l6.qza

qiime composition ancom \
  --i-table comp-NoYeast-table-l6.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization l6-ancom-NoYeast-PersonNumber-1.qzv

#ancom-yeast
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-where "[process]='yeast'" \
  --o-filtered-table yeast-table.qza

qiime composition add-pseudocount \
  --i-table yeast-table.qza \
  --o-composition-table comp-yeast-table.qza

qiime composition ancom \
  --i-table comp-yeast-table.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization ancom-yeast-PersonNumber-1.qzv

qiime taxa collapse \
  --i-table yeast-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table yeast-table-l6.qza

qiime composition add-pseudocount \
  --i-table yeast-table-l6.qza \
  --o-composition-table comp-yeast-table-l6.qza

qiime composition ancom \
  --i-table comp-yeast-table-l6.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization l6-ancom-yeast-PersonNumber-1.qzv
	
#ancom-sour
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file proj-metadata.tsv \
  --p-where "[process]='sour'" \
  --o-filtered-table sour-table.qza

qiime composition add-pseudocount \
  --i-table yeast-table.qza \
  --o-composition-table comp-sour-table.qza

qiime composition ancom \
  --i-table comp-sour-table.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization ancom-sour-PersonNumber-1.qzv

qiime taxa collapse \
  --i-table sour-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table sour-table-l6.qza

qiime composition add-pseudocount \
  --i-table sour-table-l6.qza \
  --o-composition-table comp-sour-table-l6.qza

qiime composition ancom \
  --i-table comp-sour-table-l6.qza \
  --m-metadata-file proj-metadata.tsv \
  --m-metadata-column PersonNumber-1 \
  --o-visualization l6-ancom-sour-PersonNumber-1.qzv

#Correlation Analysis
qiime tools import \
  --input-path feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table-biom-filtered.qza

qiime SCNIC sparcc-filter \
  --i-table feature-table-biom-filtered.qza \
  --o-table-filtered feature-table-filtered.qza

qiime SCNIC calculate-correlations \
  --i-table feature-table-filtered.qza \
  --p-method sparcc \
  --o-correlation-table feature-table_correls.qza

qiime SCNIC build-correlation-network-r \
  --i-correlation-table feature-table_correls.qza \
  --p-min-val .35 \
  --o-correlation-network feature_net.qza

qiime SCNIC make-modules-on-correlations \
  --i-correlation-table feature-table_correls.qza \
  --i-feature-table feature-table-biom-filtered.qza \
  --p-min-r .35 \
  --o-collapsed-table feature-table_collapsed.qza \
  --o-correlation-network feature_net.modules.qza \
  --o-module-membership feature_membership.qza

qiime metadata tabulate \
  --m-input-file feature_membership.qza \
  --o-visualization feature_membership.qzv

