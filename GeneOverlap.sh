#!/bin/bash

phenotype=("Epilepsy" "GGE")

subnetworks="/cluster/"

for v in "${phenotype[@]}"; do
    echo "phenotype: $v"
    ls ${subnetworks}*.tsv | parallel -j5  Rscript GeneOverlap.R {} ${v}
done


