#!/bin/bash

phenotype=("Epilepsy" "GGE")

subnetworks="/mnt/isilon/projects/isbsequencing/epi_rare/net_06032025/cluster/"

for v in "${phenotype[@]}"; do
    echo "phenotype: $v"
    ls ${subnetworks}*.tsv | parallel -j5  Rscript /mnt/isilon/projects/isbsequencing/epi_rare/network_29092025/GeneOverlap.R {} ${v}
done


