#!/bin/bash

# Description: This script automatically downloads NCBI SRA reads
# Author: Asad Prodhan PhD
# Email: prodhan82@gmail.com
# Date: 2024-07-01
# Version: 1.0


# File containing the list of SRA accession numbers
SRA_LIST="SRR_Acc_List.txt"

# Loop through each accession number in the list
while IFS= read -r accession; do
    echo "Processing $accession"
    prefetch.3.1.1 $accession && fasterq-dump.3.1.1 $accession --outdir reads
done < "$SRA_LIST"

# The end
