# FBCseq

R package for differential HCDR3 abundance analysis of rabbit antibody repertoires selected by **FBC-seq** (whole-cell bead/FACS panning coupled to next-generation sequencing). Given paired experimental (selected) and control (unselected/non-target) FASTQ libraries, FBCseq identifies antibody heavy-chain CDR3 (HCDR3) sequences enriched by selection and designs primers to recover the corresponding clones by inverse ("around-the-horn") PCR.

Used in Cyr et al., *J Mol Biol* 2023 — see [Citation](#citation).

## Contents

- [How it works](#how-it-works)
- [Installation](#installation)
- [Data availability](#data-availability)
- [Quick start](#quick-start)
- [Function reference](#function-reference)
- [Adapting to other repertoires](#adapting-to-other-repertoires)
- [Citation](#citation)
- [License](#license)

## How it works

Reads are single-end MiSeq reads that span the HCDR3 in reverse orientation, starting from the CH1/constant region immediately downstream of the VH domain. For each input FASTQ, `differential_HCDR3()`:

1. Reverse-complements the raw reads and discards any without the expected CH1 reverse primer (`filter_embedded_primer()`).
2. Frame-corrects and translates reads to protein.
3. Splices out the HCDR3 loop using flanking junction motifs and discards sequences that are out of frame, contain stop codons, or lack the expected FW4 tail (`splice_rabbit_HCDR3s()`, `extract_HCDR3s()`).
4. Builds a count matrix of unique HCDR3s shared across all experimental and all control replicates, and runs **DESeq2** to test for differential abundance between the two groups.
5. Returns per-HCDR3 fold-change, p-value, and normalized abundance in each group, along with template DNA sequences for primer design.

Enriched clones can then be recovered from the pooled selection output by inverse PCR: `get_primer()` trims a DNA template from its 3' end until the calculated melting temperature (`tm_calc()`) drops below a target annealing temperature, producing forward/reverse primers for that clone.

## Installation

```r
library(devtools)
install_github("hedawils/FBCseq")
```

**Dependencies:** tidyverse, magrittr, seqinr, R.utils, ShortRead, Biostrings, DESeq2

## Data availability

- Published, processed outputs: GEO accession [GSE222897](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222897)
- Raw FASTQ reads: SRA project [PRJNA923794](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA923794&o=acc_s%3Aa)

Raw FASTQs can be downloaded from the ENA with the included script (Linux or macOS):

```bash
bash PRJNA923794_fastq_download.sh
```

## Quick start

```r
## set working directory to the folder containing the fastq files
setwd("path/to/fastqs")

library(FBCseq)

## run the differential HCDR3 abundance pipeline
H929_results <- differential_HCDR3(
  experimental_fastq_files = c("SRR23080982.fastq.gz", "SRR23080981.fastq.gz", "SRR23080980.fastq.gz"),
  control_fastq_files      = c("SRR23080976.fastq.gz", "SRR23080975.fastq.gz", "SRR23080974.fastq.gz"),
  primer = "GCCCTTGGTGGAGGC"
)

H929_df <- H929_results[[1]]

## select clones for recovery: 8-fold enriched, p < 1e-4, abundance > 1e-4
H929_primer_set <- H929_df[tidyr::replace_na(
  H929_df$mean_targ_normalized > 1e-4 & H929_df$pvalue < 1e-4 & H929_df$log2FoldChange > 3,
  FALSE
), ]

## design inverse-PCR recovery primers for the selected clones
H929_primer_set$F_primer <- sapply(H929_primer_set$F_primer_template, get_primer, annealing = 62)
H929_primer_set$R_primer <- sapply(H929_primer_set$R_primer_template, get_primer, annealing = 62)
```

`differential_HCDR3()` returns a list of two elements:

1. `results_ordered` — one row per HCDR3, with DESeq2 statistics (`log2FoldChange`, `pvalue`, `padj`, `baseMean`), normalized abundance in target (`mean_targ_normalized`) and non-target (`mean_non_targ_normalized`) libraries, and DNA templates (`DNA`, `F_primer_template`, `R_primer_template`) for primer design.
2. `HCDR3_libs` — per-library HCDR3 count tables, one per input FASTQ file.

## Function reference

| Function | Purpose |
|---|---|
| `differential_HCDR3()` | End-to-end pipeline: filter reads, extract HCDR3s, run DESeq2 differential abundance analysis |
| `filter_embedded_primer()` | Keep reads containing the expected CH1 reverse primer and trim to the correct reading frame |
| `splice_rabbit_HCDR3s()` | Extract the HCDR3 substring from a translated sequence using flanking junction motifs |
| `extract_HCDR3s()` | Clean spliced HCDR3s (remove out-of-frame, stop-codon, or mistranslated sequences) |
| `get_primer()` | Trim a DNA template to a primer of a target annealing temperature |
| `tm_calc()` | Calculate the melting temperature of a DNA oligo |

Full argument documentation is available via `?differential_HCDR3`, etc., after installation.

## Adapting to other repertoires

The rabbit VH junction motifs are configurable, so the pipeline can be repointed at other species' antibody repertoires or other protein/peptide engineering panning workflows by adjusting these `differential_HCDR3()` arguments:

- `upstream_junction` / `downstream_junction` — protein-sequence regexes flanking the region of interest
- `tail_pattern` — fixed downstream sequence used to discard incompletely extracted reads
- `split_pattern` — sequence used to locate the primer-design split point

## Citation

Cyr MG, Wilson HD, Spierling AL, Chang J, Peng H, Steinberger P, Rader C. Concerted Antibody and Antigen Discovery by Differential Whole-cell Phage Display Selections and Multi-omic Target Deconvolution. *J Mol Biol.* 2023;435(10):168085. doi: [10.1016/j.jmb.2023.168085](https://doi.org/10.1016/j.jmb.2023.168085). PMID: [37019174](https://pubmed.ncbi.nlm.nih.gov/37019174/).

## License

MIT — see [LICENSE](LICENSE).
