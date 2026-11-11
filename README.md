
# Functional Genomics of L. fermentum B1 28

This repository contains the R scripts used to analyze the genome of *Limosilactobacillus fermentum* B1 28 (ATCC 14931) for the manuscript: "Functional genomics of the type strain..."

## Requirements
- R (version 4.0+)
- Packages: `tidyverse`, `readxl`, `Biostrings`, `GenomicRanges`, `rtracklayer`, `GO.db`

## Data Inputs
To run this script, you must provide:
1. `GCF_005341425.1_ASM534142v1_genomic.fna` (Download from NCBI)
2. `genomic.gff` (Download from NCBI)
3. `emapper_annotations.xlsx` (Generated via EggNOG-mapper)

## Outputs
Running the script generates:
- `figures/`: All PDF/PNG plots included in the manuscript.
- `tables/`: CSV files for all supplementary data (Enzymes, Anticancer candidates, GO, KEGG, COG).
