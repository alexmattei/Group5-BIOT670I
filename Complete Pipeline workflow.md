#######  Complete Pipeline workflow. Run each command sequentially, except when indicated otherwise.  #######

## Setup: Run the installation commands in the installs file to set up the environments (sraenv and qiime) and install needed programs
# Make a directory for the analysis.
mkdir metagenome_dir

# enter the directory
cd metagenome_dir

# Using NCBI Run Selector, pick the desired samples from your study and download the metadata and accession list files
# They will probably be named SraRunTable.txt and SRR_Acc_List.txt

### STEP 1: Download fastq data
# Run the prefetch_fasterq-dump program
# may need to make executable first, if so, run chmod +x prefetch_fasterq-dump.sh
conda activate sraenv

./prefetch_fasterq-dump SRR_Acc_List.txt

# once done, enter the reads directory. Display all of the file names with ls command
ls
### WARNING: If your data are on the reverse strand, run this loop to reverse complement them to the forward strand.
# Do not run these commands otherwise.
for FILE in *.fastq; do cat $FILE | reformat.sh in=$FILE out=$FILE.rev rcomp=t; done
# move the reversed files to another directory with the mv command (recommend keeping separate)
mkdir reversed
mv *.rev /new/directory/for/trimmed/

#rename them back to fastq
for FILE in *.fastq.rev; do mv "$FILE" "${FILE%.fastq.rev}.fastq"; done

### STEP 2: Quality visualization:
# generate fastqc reports for all files as so, in the directory with fastq files:
fastqc *.fastq


# view files with firefox <file.html> or unzip it.
firefox <filename>.html


### STEP 3: Quality Control:
# BBDUK informed command: The V35 region is about 548 Thus,
# the force trim right should be about this long (e.g. 540). 
# make sure the adapters file is in the same folder as your fastq.
for FILE in *.fastq; do cat $FILE | bbduk.sh in=$FILE out=$FILE.trim ref=adapters.fa  k=23 ktrim=r mink=11 edist=1 qtrim=rl trimq=28 forcetrimright=540; done
# make directory for trimmed files
mkdir trimmed
cd trimmed

mv *.trim /current/dir/trimmed/

# rename ones with .trim to fastq again
for FILE in *.fastq.trim; do mv "$FILE" "${FILE%.fastq.trim}.fastq"; done

## do another fastqc report
fastqc <filenames>

# view with firefox as previously
# turn off sraenv environment
conda deactivate

### STEP 4: Import data into Qiime2
# Manifest may be made manually, or you can use the included interactive python script manifestmakeq2.py
# make sure the SRA metadata file is in the same folder as manifestmakeq2, then run:
python manifestmakeq2.py

# Also recommend creating Qiime2 compatible metadata file now. Run the provided script name metamakeq2.py
python metamakeq2.py


# activate your qiime environment
conda activate qiime2-amplicon-2024.5
  
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest-q2.tsv \
  --output-path seqs-trim.qza \
  --input-format SingleEndFastqManifestPhred33V2

mkdir qiimestart
# copy seqs-trim.qza to qiimeclustering with cp command
cd qiimestart

### STEP 5: v search command to dereplicate 
qiime vsearch dereplicate-sequences \
  --i-sequences seqs-trim.qza \
  --o-dereplicated-table table.qza \
  --o-dereplicated-sequences rep-seqs.qza
  
### STEP 6: Clustering section: Uses closed reference OTU clustering from vsearch.
# This method will generally identify fewer clusters than open reference or de novo clustering,
# but the sequences identified will have a high level of confidence.
  
# Recommend closed ref with 99% threshold. The database used is found at: http://ftp.microbio.me/greengenes_release/2022.10/
# recommend using the fna full length backbone one, in .qza. Must download a taxonomy backbone too: file name: 2022.10.backbone.tax.qza 

qiime vsearch cluster-features-closed-reference \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --i-reference-sequences 2022.10.backbone.full-length.fna.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table table-cr-99.qza \
  --o-clustered-sequences rep-seqs-cr-99.qza \
  --o-unmatched-sequences unmatched-cr-99.qza \
  --p-threads 6

# make new folder and move into
mkdir clustered99

  
### STEP 7 OPTIONAL: Pre-rarefaction diversity and taxonomy: Qualitative, more useful for seeing if 
# rarefaction distorts the group and sample-level picture of diversity. 

mkdir initialmetrics
cd initialmetrics
# Also needs backbone taxonomy, just copy or move with file manager

qiime diversity alpha \
  --i-table table-cr-99.qza \
  --p-metric observed_features \
  --o-alpha-diversity observed_features_vector.qza
  
qiime diversity alpha-group-significance \
  --i-alpha-diversity observed_features_vector.qza \
  --m-metadata-file metadataq2.tsv \
  --o-visualization alpha-group-sig-obs-feats.qzv
  
qiime taxa collapse \
    --i-table table-cr-99.qza \
    --i-taxonomy 2022.10.backbone.tax.qza \
    --p-level 2 \
    --o-collapsed-table prerare-table-l2.qza
  
qiime feature-table summarize \
   --i-table prerare-table-l2.qza \
   --o-visualization viz-prerare-table-l2.qzv

qiime taxa barplot \
  --i-table prerare-table-l2.qza \
  --m-metadata-file metadataq2.tsv \
  --o-visualization prerare-phylum-bar-plots.qzv
  

  
### STEP 8: Alpha rarefaction: This is done to determine at what sampling depth diversity "levels off". The level off
# point is a good point to set the depth since few features are likely to be found after it.

# Start with alignment to get a rooted tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-cr-99.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
# Use the highest sampling depth amongst the samples in your dataset for creating rarefaction curve.
# Recommend looking at a visualization of frequency data table (viz-table-cr99.qzv), this will show
# the sample with the highest depth. Use this to generate the plot. For example 1124 here
qiime feature-table summarize \
   --i-table table-cr-99.qza \
   --o-visualization viz-table-cr-99.qzv
   
qiime diversity alpha-rarefaction \
  --i-table table-cr-99.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 1124 \
  --m-metadata-file metadataq2.tsv \
  --o-visualization alpha-rarefaction.qzv
  
# After looking at rarefaction (alpha-rarefaction.qzv), find the point where feature diversity levels off
# The diversity is likely captured at this point. However, be careful not to lose too many samples that do not reach
# that depth since they will be excluded (See Moving Pictures Qiime2 tutorial for more).
 

### STEP 9: Creates the essential files for downstream diversity analysis. Some of these can be viewed as well (.qzv).
# run core metrics phylogenetic:
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-cr-99.qza \
  --p-sampling-depth 500 \
  --m-metadata-file metadataq2.tsv \
  --output-dir core-metrics-results
  
cd core-metrics-results.

### STEP 10: Alpha Diversity Analysis
## Phylogenetic alpha diversity

qiime diversity alpha-group-significance \
  --i-alpha-diversity faith_pd_vector.qza \
  --m-metadata-file metadataq2.tsv \
  --o-visualization faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity evenness_vector.qza \
  --m-metadata-file metadataq2.tsv \
  --o-visualization evenness-group-significance.qzv

mkdir alphadiv
# move evenness-group-significance.qzv and faith-pd-group-significance.qzv to alphadiv

### STEP 11: Beta Diversity Analysis

## phylogenetic beta diversity
# move back to the core-metrics-results folder in terminal


qiime diversity beta-group-significance \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadataq2.tsv  \
  --m-metadata-column group \
  --o-visualization unweighted-unifrac-group-significance.qzv \
  --p-pairwise

# Now do beta group significance for male vs female, if applicable. Change column name as needed.
qiime diversity beta-group-significance \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadataq2.tsv  \
  --m-metadata-column sex \
  --o-visualization unweighted-unifrac-sex-significance.qzv \
  --p-pairwise
 
 mkdir betadiv
 # place output files into betadiv foler
  
### STEP 12: Taxa collapse: 
# Essential step. In general, OTUs are mapped to a taxon, but lack actual taxonomic annotation.
# Taxa collapse provides this at the desired level, summing so the total feature count stays the same
# regardless of the level (e.g. phylum, class, species)

# Post rarefaction, using the rarefied_table from core-metrics-diversity:
# Collapse to given taxonomic levels using the Greengenes taxonomic backbone file.
# Here, l2, l6, and l7 are phylum, genus, and species, respectively.
# Bring the backbone.tax file into the core-metrics folder.

qiime taxa collapse \
    --i-table rarefied_table.qza \
    --i-taxonomy 2022.10.backbone.tax.qza \
    --p-level 7 \
    --o-collapsed-table rarefied-table-l7.qza
    
qiime taxa collapse \
    --i-table rarefied_table.qza \
    --i-taxonomy 2022.10.backbone.tax.qza \
    --p-level 2 \
    --o-collapsed-table rarefied-table-l2.qza
    
qiime taxa collapse \
    --i-table rarefied_table.qza \
    --i-taxonomy 2022.10.backbone.tax.qza \
    --p-level 6 \
    --o-collapsed-table rarefied-table-l6.qza
 
mkdir filtered
# move all of the rarefied tables with an l7, etc number to filtered folder.

### STEP 13: Filtering. Per Bokulich et al. 2013, recommend filtering features with abundance less than 0.005%. To do so, first
# summarize the data from table-l7 to get an idea of composition. For example, With 22,326 features, 0.005% is 1.1163
# so filter singletons at least. 
# Rarefied version has 10150 features, so singleton still applies.


# Recommend filtering singletons first (taxa occurring once in only one sample). 
qiime feature-table filter-features \
 --i-table rarefied-table-l7.qza \
 --p-min-frequency 2 \
 --o-filtered-table singleton-filtered-rarefied-table-l7.qza
 
# for genus
qiime feature-table filter-features \
  --i-table rarefied-table-l6.qza \
  --p-min-frequency 2 \
  --o-filtered-table singleton-filtered-rarefied-table-l6.qza
  
# for phylum
qiime feature-table filter-features \
  --i-table rarefied-table-l2.qza \
  --p-min-frequency 2 \
  --o-filtered-table singleton-filtered-rarefied-table-l2.qza


# move the  singleton-filtered-rarefied-table files to the "filtered" folder.

### STEP 14: Qualitative analysis with Barplots
# Requires the metadata file used earlier for taxa collapse.

qiime taxa barplot \
  --i-table singleton-filtered-rarefied-table-l7.qza \
  --m-metadata-file metadataq2.tsv \
  --o-visualization species-bar-plots.qzv
  
qiime taxa barplot \
  --i-table singleton-filtered-rarefied-table-l6.qza \
  --m-metadata-file metadataq2.tsv \
  --o-visualization genus-bar-plots.qzv
  
qiime taxa barplot \
  --i-table singleton-filtered-rarefied-table-l2.qza \
  --m-metadata-file metadataq2.tsv \
  --o-visualization phylum-bar-plots.qzv

### STEP 15:  differential abundance using ANCOM-BC Qiime2 plugin.
 
# Run the ancombc function on case vs control column, here called group
# Note that no significant species level differences were found at 0.05 or lower

# species
qiime composition ancombc \
  --i-table singleton-filtered-rarefied-table-l7.qza \
  --m-metadata-file metadataq2.tsv \
  --p-formula 'group' \
  --p-reference-levels group::None \
  --o-differentials ancombc-group-l7.qza
  
# genus
qiime composition ancombc \
  --i-table singleton-filtered-rarefied-table-l6.qza \
  --m-metadata-file metadataq2.tsv \
  --p-formula 'group' \
  --p-reference-levels group::None \
  --o-differentials ancombc-group-l6.qza
  
#phlyum
qiime composition ancombc \
  --i-table singleton-filtered-rarefied-table-l2.qza \
  --m-metadata-file metadataq2.tsv \
  --p-formula 'group' \
  --p-reference-levels group::None \
  --o-differentials ancombc-group-l2.qza

### STEP 16: Create ANCOM-BC barplots with certain sig threshold, can be altered if needed.
# On project dataset with p=0.001, differences were found at the phylum level.

mkdir ancombc
# move the ancombc-group-l(n).qza files to ancombc folder
cd ancombc

 qiime composition da-barplot \
  --i-data ancombc-group-l7.qza \
  --p-significance-threshold 0.001 \
  --p-level-delimiter ';' \
  --o-visualization da-barplot-group-l7.qzv
  
# for genus
 qiime composition da-barplot \
  --i-data ancombc-group-l6.qza \
  --p-significance-threshold 0.001 \
  --p-level-delimiter ';' \
  --o-visualization da-barplot-group-l6.qzv
   
# For phylum
 qiime composition da-barplot \
  --i-data ancombc-group-l2.qza \
  --p-significance-threshold 0.001 \
  --p-level-delimiter ';' \
  --o-visualization da-barplot-group-l2.qzv
  
  
### STEP 17: Picrust2 Qiime2 plugin workflow

# Recommend using the rarefied table from core metrics and the rep-seqs file from OTU clustering outputs,
# though it is optional and the non-rarefied table (table-cr-99.qza
# This command may take some time to run, at least 10-15 minutes. Different placement method (epa-ng) recommended
# If your computer has very good RAM, i.e. more than 16 GB.
qiime picrust2 full-pipeline \
   --i-table rarefied_table.qza \
   --i-seq rep-seqs-cr-99.qza \
   --output-dir picrustouts \
   --p-placement-tool sepp \
   --p-threads 6 \
   --p-hsp-method mp \
   --p-max-nsti 2 \
   --p-edge-exponent 0 \
   --verbose

# move to the newly created folder:
cd picrustouts

# Filter picrust results

qiime feature-table filter-features \
  --i-table pathway_abundance.qza \
  --p-min-frequency 2 \
  --o-filtered-table singleton-filtered-pathway_abundance.qza
  
qiime feature-table filter-features \
  --i-table ec_metagenome.qza \
  --p-min-frequency 2 \
  --o-filtered-table singleton-filtered-ec_metagenome.qza
  
qiime feature-table filter-features \
  --i-table ko_metagenome.qza \
  --p-min-frequency 2 \
  --o-filtered-table singleton-filtered-ko_metagenome.qza
  
# STEP 18: Optional: visualize the results, of filtered picrust2 outputs:

qiime feature-table summarize \
   --i-table singleton-filtered-pathway_abundance.qza \
   --o-visualization viz_pathway_abundance.qzv
   
qiime feature-table summarize \
   --i-table singleton-filtered-ec_metagenome.qza \
   --o-visualization viz_ec_metagenome.qzv
   
qiime feature-table summarize \
   --i-table singleton-filtered-ko_metagenome.qza \
   --o-visualization viz_ko_metagenome.qzv
   

# One picrust option is to due diversity analysis on the predictions.
qiime diversity core-metrics \
   --i-table singleton-filtered-pathway_abundance.qza \
   --p-sampling-depth n \
   --m-metadata-file metadataq2.tsv \
   --output-dir pathabun_core_metrics_out \
   --p-n-jobs 6

### STEP 19: Export For use in the ggpicrust2 R package, qza files must be converted to biom and then tsv

# metacyc
qiime tools export \
   --input-path singleton-filtered-pathway_abundance.qza \
   --output-path metacyc-pathabun_exported
  
cd metacyc-pathabun_exported

biom convert \
   -i feature-table.biom \
   -o metacyc-pathabun.biom.tsv \
   --to-tsv
   
# KO
# return to picrustouts folder

qiime tools export \
   --input-path singleton-filtered-ko_metagenome.qza \
   --output-path ko-metagenome_exported
   
cd ko-metagenome_exported

biom convert \
   -i feature-table.biom \
   -o ko-metagenome.biom.tsv \
   --to-tsv

# Now move the .tsv files to a directory in R that contains the ggpicrust2 script for the pipeline.
# The ggpicrust2 portion can be run in Windows if desired using R Studio, simply
# email the metacyc and other files to yourself.
# With a text editor, the metacyc pathway and KO file will need slight changes: remove the top row that starts with a #, and 
# then change the name of the first column form OTU ID to pathway (for the metacyc file).
# The next pipeline steps are to be completed in R, following the instructions in the scripts.

  




##########
### END
##########

  
  



