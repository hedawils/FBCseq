# FBCseq

This package enables differential abundance analysis of rabbit HCDR3 sequences, identified via the FBC-seq whole cell panning method (NGS); publication pending.

Package download via devtools in R:
```
library(devtools)
install_github("hedawils/FBCseq")
```

Reads should span the HCDR3 in a reverse orientation, starting from CH1 or constant sequence immediately downstream of VH domain.  Package includes primer design capability for around-the-horn based recovery of single clones from mixed selection outputs via PCR.

Published raw FASTQ files from publication can be downloaded here: XXXXXX

---
#Pipeline for differential HCDR3 analysis and primer design
---

```
## navigate to folder where files are located

## initiate differential abundance analysis
H929_results <- differential_HCDR3(c("314_1_1.fastq", "314_2_1.fastq", "314_3_1.fastq"), c("314_7_1.fastq", "314_8_1.fastq", "314_9_1.fastq"), primer = "GCCCTTGGTGGAGGC")

## generate volcano plot
H929_plot <- volcano_HCDR3s(H929_results, )

## select sequences for recovery based upon 8-fold enriched, pvalue < 1e-4, and absolute abundance > 1e-4
H929_primer_set <- H929_results[H929_results$mean_targ_normalized > 1e-4 & H929_results$pvalue < 1e-4 & H929_results$log2FoldChange > 3)
H929_primer_set$F_primer <- sapply(H929_primer_set$F_primer_template, get_primer(annealing = 62))
H929_primer_set$R_primer <- sapply(H929_primer_set$R_primer_template, get_primer(annealing = 62))
```

Modifications of the following parameters to the differential_HCDR3 function could enable generalizable extension of this differential abundance pipeline to other protein engineering applications:
upstream_junction, downstream_junction, tail_pattern, split_pattern 





