# FBCseq

This package enables differential abundance analysis of rabbit HCDR3 sequences, identified via the FBC-seq whole cell panning method (NGS); publication pending.

Package download via devtools in R:
```
library(devtools)
install_github("hedawils/FBCseq")
```

Reads should span the HCDR3 in a reverse orientation, starting from CH1 or constant sequence immediately downstream of VH domain.  Package includes primer design capability for around-the-horn based recovery of single clones from mixed selection outputs via PCR.

Published raw FASTQ files from publication can be downloaded from GEO under accession [GSE222897](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222897).

BASH script for downloading raw FASTQs via CURL from the SRA can be run as following:
```
bash GSE222897_fastq_download.sh
```


---
Pipeline for differential HCDR3 analysis and primer design
---

```
## navigate to folder where files are located via setwd()

## initiate differential abundance analysis
H929_results <- differential_HCDR3(c("314_1_1.fastq", "314_2_1.fastq", "314_3_1.fastq"), c("314_7_1.fastq", "314_8_1.fastq", "314_9_1.fastq"), primer = "GCCCTTGGTGGAGGC")

H929_df <- H929_results[[1]]

## select sequences for recovery based upon 8-fold enriched, pvalue < 1e-4, and absolute abundance > 1e-4
H929_primer_set <- H929_df[H929_df$mean_targ_normalized > 1e-4 & H929_df$pvalue < 1e-4 & H929_df$log2FoldChange > 3,]

### design primers for inverse PCR recovery
H929_primer_set$F_primer <- sapply(H929_primer_set$F_primer_template, get_primer(annealing = 62))
H929_primer_set$R_primer <- sapply(H929_primer_set$R_primer_template, get_primer(annealing = 62))
```

Modification of the following parameters to the differential_HCDR3 function could enable extension of this differential abundance pipeline to other protein/peptide engineering applications, including antibody repertoires from other species:
upstream_junction, downstream_junction, tail_pattern, split_pattern 





