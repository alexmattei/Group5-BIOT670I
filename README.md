# Group5-BIOT670I
Repository for capstone project.
Group 5 members are Holly Winbigler, Teresa Acevedo, and Alex Mattei.

This project aims to profile the functional capabilities of gut microbiomes in IBD patients and healthy controls using metagenomic data.

Microbial functions and pathways associated with IBD will be identified and compared.

A bioinformatics pipeline will be developed to perform data quality control and pre-processing, annotate metagenomic data, and visualize the functional differences between patient and control groups. This pipeline requires a relatively recent (~2022 or higher) Linux distribution to run, and also requires an up to date R installation (4.3 or higher). R Studio is highly recommended. 

The key files needed to use the pipeline are visible on the main page. These include the dependencies installation file (installs_EZ.txt), prefetch_fasterq-dump.sh by Dr. Asad Prodhan, the Complete Pipeline Commands File, two python files for making manifest and metadata files written by Group 5, and the R script files needed for library installation and plot generation of metabolic pathways (citation in file). In the Additional Files folder, there is a file from BBMAP and used by BBDuk2 called adapters.fa. The folder also contains example Metadata and SRA accession list files for use with the pipeline.
