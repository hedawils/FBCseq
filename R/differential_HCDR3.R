#' Calculate depmap score mean for depmap cell lines, and N per gene
#'
#' Function to conduct differential HCDR3 abundance analysis of single end MiSeq reads
#' from FBCseq experiment
#'
#' @param experimental_fastq_files vector of names of experimental fastq files for analysis. Ordered to be paired with control_fastq_files
#' @param control_fastq_files vector of names of control fastq files for analysis. Ordered to be paired with experimental_fastq_files
#' @param primer CH1 hybridizing reverse primer. 
#'
#'
#' @export
#'
#' @example results <- differential_HCDR3(experimental_fastq_files = c(), control_fastq_files = c(), primer = "GCCCTTGGTGGAGGC")
#' ### download example fastqs from ENA
#' 
#' 



differential_HCDR3 <- function(experimental_fastq_files, control_fastq_files, primer = "GCCCTTGGTGGAGGC"){
  
  library(tidyverse)
  library(magrittr)
  library(seqinr)
  library(R.utils)
  library(ShortRead)
  library(Biostrings)
  library(DESeq2)
  
  library_names <- c(experimental_fastq_files, control_fastq_files)
  
  if (!all(library_names %in% list.files(pattern = ".fastq"))){
    stop("not all fastq files detected in directory. need to change directory or ensure input parameters are correct")
  }
  
  HCDR3_libs <- list()
  DNA_libs <- list()
  HCDR3s <- list()
  
  for (i in 1:length(library_names)){
    DNA <- Biostrings::reverseComplement(ShortRead::sread(ShortRead::readFastq(library_names[i])))
    DNA <- filter_embedded_primer(DNA, primer)
    aa <- Biostrings::translate(Biostrings::DNAStringSet(DNA), if.fuzzy.codon = "X")
    HCDR3 <- extract_HCDR3s(aa)
    HCDR3_libs[[i]] <- as.data.frame(table(HCDR3)) %>% dplyr::mutate_at(vars(HCDR3), as.character)
    names(HCDR3_libs[[i]])[2] <- paste0("Count_", library_names[i]) %>% stringr::str_remove_all(".fastq|.fastq.gz")
    HCDR3s[[i]] <- HCDR3_libs[[i]]$HCDR3
    DNA_libs[[i]] <- as.data.frame(table(DNA)) %>%dplyr:: mutate_at(vars(DNA), as.character)
    
    names(DNA_libs[[i]])[2] <- paste0("Count_", library_names[i]) %>% stringr::str_remove_all(".fastq|.fastq.gz")
    
    if ((i == 1) | (i == length(experimental_fastq_files) + 1)) {
      DNA_libs[[i]]$AA <- DNAStringSet(DNA_libs[[i]]$DNA) %>% translate(if.fuzzy.codon = "X") %>% as.character()
      DNA_libs[[i]]$HCDR3 <- sapply(strsplit(sapply(strsplit(as.character(DNA_libs[[i]]$AA), "WG.G"), "[", 1), 
                                             ".TYFC|A.YFC|AT.FC|ATY.C"), "[", 2)
      names(DNA_libs[[i]])[2] <- "counts"
    }
    
    
  }
 
  ### obtain control matrix - keep sequences common to either all experimental or control outputs
  control_HCDR3s <- Reduce(intersect, HCDR3s[1:length(control_fastq_files)])
  experimental_HCDR3s <- Reduce(intersect, HCDR3s[(1+length(control_fastq_files)):length(HCDR3s)])
  
  unique_HCDR3s <- c(control_HCDR3s, experimental_HCDR3s) %>% unique
  unique_HCDR3s <- unique_HCDR3s[nchar(unique_HCDR3s) > 1]
  
  ## filter down HCDR3_libs
  for (i in 1:length(HCDR3_libs)){
    HCDR3_libs[[i]] <- HCDR3_libs[[i]][HCDR3_libs[[i]]$HCDR3 %in% unique_HCDR3s,]
  }
  
  matrix <- HCDR3_libs %>% reduce(dplyr::full_join, by = "HCDR3") 
  row.names(matrix) <- matrix$HCDR3
  matrix$HCDR3 <- NULL
  matrix <- matrix[order(matrix %>% rowSums(na.rm=TRUE), decreasing = TRUE),]
  matrix[is.na(matrix)] <- 0
  
  countdata <- as.matrix(matrix)
  #Condition for deseq2 paired analysis
  #Levels of substrate factor correspond to user-defined target and non-target
  #substrate names.
  #Technical_replicate stores pairing information of replicates consistent with
  #library_names variable.
  substrate <- factor(c(rep("target", length(experimental_fastq_files)), 
                        rep("non.target", length(control_fastq_files))),
                        levels = c("non.target", "target"))
  #technical_replicate <- factor(c(rep(as.character(c(1:3)), 2)))
  ##### The above line is included for paired analyses, which requires a
  #corresponding modification of the design paramater below.
  coldata <- data.frame(row.names=colnames(countdata), substrate)
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=
                                  ~substrate)
  dds <- DESeq(dds, fitType='local')
  results <- results(dds)
  #save(dds, file = paste0("dds_", exp, "_", analysis))
  #Append results df with HCDR3 and mean abundance in target and non-target outputs
  results$HCDR3 <- as.character(row.names(matrix))
  results$mean_non_targ_normalized <- rowMeans(as.matrix(counts(dds,
                                                                normalized=TRUE))[,which(as.vector(substrate) == "non.target")])
  results$mean_targ_normalized <- rowMeans(as.matrix(counts(dds,
                                                            normalized=TRUE))[,which(as.vector(substrate) == "target")])
  
  results$mean_non_targ_normalized <- results$mean_non_targ_normalized/sum(results$mean_non_targ_normalized)
  results$mean_targ_normalized <- results$mean_targ_normalized/sum(results$mean_targ_normalized)
  
  #Order results by pvalue
  results_ordered <- as.data.frame(results[order(-results$baseMean),])
  
  #Extract DNA sequences from DNA_libs
  DNA_data <-rbind(DNA_libs[[1]], DNA_libs[[length(experimental_fastq_files) + 1]])
  DNA_data <- DNA_data %>% arrange(desc(counts))
  
  results_ordered$DNA_sequence <- DNA_data$DNA[match(results_ordered$HCDR3, DNA_data$HCDR3)]

  read_length_no_primer <- results_ordered$DNA_sequence %>% nchar %>% unique
  
  ### forward primers start from second residue in FW4
  results_ordered$FW4_primer <- results_ordered$DNA_sequence %>% 
    substr(read_length_no_primer-29,read_length_no_primer-11)
  
  ## reverse primer template to be trimmed to appropriate length after 
  ## subsetting of sequences for mAb recovery
  results_ordered$primer_template <- results_ordered$DNA_sequence %>% 
    substr(1,read_length_no_primer-30) %>% 
    DNAStringSet %>% reverseComplement %>% as.character
  
  
  return(results_ordered)
  
}
c

