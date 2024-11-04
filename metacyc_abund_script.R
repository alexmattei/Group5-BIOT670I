
### This script will produce metacyc pathway outputs. It should be run in order, top to bottom. 
### The needed files, namely qiime_metadata.tsv and metacyc_pathway_abundance.tsv should be in the same folder.
### All of these commands are modified or original versions of the ggpicrust2 commands listed on the github page:
### https://github.com/cafferychen777/ggpicrust2
### Also makes use of ggfun By: D Wang, G Chen, L Li, S Wen, Z Xie, X Luo, L Zhan, S Xu, J Li, R Wang, Q Wang, G Yu. Reducing language barriers, promoting information absorption, and communication
#### using fanyi. Chinese Medical Journal. 2024, 137(16):1950-1956. doi: 10.1097/CM9.0000000000003242


# Library load commands
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(tidyverse)
library(ggh4x)
library(ggfun)

# Load pathway abundance file. Make sure it is in the same directory as this script, or provide a path to it.
abundance_file <- "metacyc-pathabun.biom.tsv"
metacyc_abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Load metadata. Ensure it is in the same directory as this script or provide path.
metadata <- read_delim(
  "metadataq2.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Run with metacyc data in this case:
results_file_input <- ggpicrust2(data = metacyc_abundance_data,
                                 metadata = metadata,
                                 group = "group",
                                 pathway = "MetaCyc",
                                 daa_method = "LinDA",
                                 ko_to_kegg = FALSE,
                                 order = "group",
                                 p_values_bar = TRUE,
                                 x_lab = "description")

results_file_input[[1]]$plot
results_file_input[[1]]$results


#### If no significant differences (at given p value) it will throw an error about the array size. At that point the 
#### data most likely do not differ, and no plots will be created. The developer of ggpicrust2 recommends increasing
#### the sample size if possible. Also check for errors in the metadata file. 
#### Also, if there are too many significant differences (over 30) an error will result as well. See below:
# Steps if an error results from too many significant differences: try a highly significant p value instead.

# the low p feature object will list the most significantly different pathways/KO's, top 20 as shown here.
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "group", daa_method = "LinDA")
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)

metacyc_daa_annotated_results_df <- metacyc_daa_annotated_results_df[!is.na(metacyc_daa_annotated_results_df$feature),]
low_p_feature <- metacyc_daa_annotated_results_df[order(metacyc_daa_annotated_results_df$p_adjust), ]$feature[1:20]
# like low_p_feature but the whole frame
lowest_metacyc_p <- metacyc_daa_annotated_results_df[order(metacyc_daa_annotated_results_df$p_adjust), ][ 1:20, c(1,2,3,4,5,6,7,8)]#[1:20]
metacyc_daa_annotated_results_df$p_adjust <- format(metacyc_daa_annotated_results_df$p_adjust, digits=5)


p <- pathway_errorbar(abundance = metacyc_abundance_data %>% column_to_rownames("pathway"),
                      daa_results_df = lowest_metacyc_p,
                      Group = metadata$group,
                      ko_to_kegg = FALSE,
                      p_values_threshold = 10,
                      order = "group",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "description")


# Click run next to p to make the plot (it is a patchwork)
p


pathway_errorbar(abundance = metacyc_abundance_data %>% column_to_rownames("pathway"), daa_results_df = metacyc_daa_annotated_results_df, Group = metadata$group, ko_to_kegg = FALSE, p_values_threshold = 0.00000000005, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")


#### Heatmaps for pathway


# Perform differential abundance analysis
metacyc_daa_results_df <- pathway_daa(
  abundance = metacyc_abundance_data %>% column_to_rownames("pathway"),
  metadata = metadata,
  group = "group",
  daa_method = "LinDA"
)

# Annotate the results
annotated_metacyc_daa_results_df <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = metacyc_daa_results_df,
  ko_to_kegg = FALSE
)

# Filter features with p < 0.05. Can be made more stringent as well.
feature_with_p_0.05 <- metacyc_daa_results_df %>% 
  filter(p_adjust < 0.05)

### Attention: only run one of the heatmap commands below, whichever is appropriate to your case.

# Now to Create the heatmap (use this command if no errors have occurred)
pathway_heatmap(
  abundance = metacyc_abundance_data %>% 
    right_join(
      annotated_metacyc_daa_results_df %>% select(all_of(c("feature","description"))),
      by = c("pathway" = "feature")
    ) %>% 
    filter(pathway %in% feature_with_p_0.05$feature) %>% 
    select(-"pathway") %>% 
    column_to_rownames("description"),
  metadata = metadata, 
  group = "group"
)

# Run this command if you had a high number of features (0ver 20). Refers back to that low_p_feature object.
pathway_heatmap(
  abundance = metacyc_abundance_data %>% 
    right_join(
      annotated_metacyc_daa_results_df %>% select(all_of(c("feature","description"))),
      by = c("pathway" = "feature")
    ) %>% 
    filter(pathway %in% low_p_feature) %>% 
    select(-"pathway") %>% 
    column_to_rownames("description"),
  metadata = metadata, 
  group = "group"
)



#### MetaCyc Pathway PCA plot (works regardless of significant feature count):

pathway_pca(abundance = metacyc_abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "group")


### optional commands to export to a text file if desired:

write.table(lowest_metacyc_p, "Metacyc_lowest_p_out.txt")

