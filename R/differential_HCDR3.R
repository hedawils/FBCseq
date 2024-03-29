#' Differential HCDR3 analysis of rabbit Ab NGS outputs
#'
#' Function to conduct differential HCDR3 abundance analysis of single end MiSeq reads
#' from FBCseq whole cell phage selections
#'
#' @param experimental_fastq_files vector of names of experimental fastq files for analysis. Ordered to be paired with control_fastq_files
#' @param control_fastq_files vector of names of control fastq files for analysis. Ordered to be paired with experimental_fastq_files
#' @param primer CH1 hybridizing reverse primer used for NGS.
#' @param split_pattern DNA sequence used for sequence trimming and primer design. Forward primer hybridizes in FW4, reverse in HCDR3. Split point at protein sequence XXXXXXW^GPCTLVTVSS
#' @param upstream_junction protein sequence REGEX upstream of rabbit HCDR3
#' @param downstream_junction protein sequence REGEX downstream of rabbit HCDR3
#' @param tail_pattern fixed sequence in FW4 of rabbit VH domains. sequences without this are removed
#'
#' @export
#'
#' @examples H929_results <- differential_HCDR3(experimental_fastq_files = c("SRR23080982.fastq.gz", "SRR23080981.fastq.gz", "SRR23080980.fastq.gz"), control_fastq_files = c("SRR23080976.fastq.gz", "SRR23080975.fastq.gz", "SRR23080974.fastq.gz"), primer = "GCCCTTGGTGGAGGC")
#' U266_results <- differential_HCDR3(experimental_fastq_files = c("SRR23080979.fastq.gz", "SRR23080978.fastq.gz", "SRR23080977.fastq.gz"), control_fastq_files = c("SRR23080976.fastq.gz", "SRR23080975.fastq.gz", "SRR23080974.fastq.gz"), primer = "GCCCTTGGTGGAGGC")


differential_HCDR3 <- function(experimental_fastq_files, control_fastq_files, 
                               primer = "GCCCTTGGTGGAGGC", split_pattern = "GGCCCAGG|GGCCAGGG", 
                               upstream_junction = ".TYFC|A.YFC|AT.FC|ATY.C", downstream_junction = "WG.G", 
                               tail_pattern = "VT.SS"){
  
  library(tidyverse)
  library(magrittr)
  library(seqinr)
  library(R.utils)
  library(ShortRead)
  library(Biostrings)
  library(DESeq2)
  library(FBCseq)
  
  library_names <- c(experimental_fastq_files, control_fastq_files)
  
  if (!all(library_names %in% list.files(pattern = ".fastq"))){
    stop("not all fastq files detected in directory. need to change directory or ensure input parameters are correct")
  }
  
  
  HCDR3_libs <- list()
  DNA_libs <- list()
  HCDR3s <- list()
  

  
  for (i in 1:length(library_names)){
    cat(paste0("processing ",library_names[i],"\n"))
    DNA <- library_names[i] %>% readFastq %>% sread %>% reverseComplement
    
    if (i == 1){
      ### detect split point for recovery primers - forward primer hybridizes in
      ### FW4, reverse in HCDR3. Split point at XXXXXXW^GPCTLVTVSS
      split_point <- strsplit(as.character(DNA), split_pattern) %>% 
        sapply(`[`, 1) %>% nchar %>% table %>% which.max %>% names %>% as.numeric
      read_length <- nchar(DNA %>% as.character) %>% table %>% which.max %>% names %>% as.numeric
    }
    
    ## filter reads with correct primer sequence
    DNA <- filter_embedded_primer(DNA, primer, read_length)
    
    # trim DNA from 5' end such that translated sequence is in frame
    DNA <- substr(DNA, 1 + split_point %% 3, read_length)
    
    ### translate and extract HCDR3s
    aa <- Biostrings::translate(Biostrings::DNAStringSet(DNA), if.fuzzy.codon = "X")
    HCDR3 <- extract_HCDR3s(aa, tail_pattern, upstream_junction, downstream_junction)
    HCDR3_libs[[i]] <- as.data.frame(table(HCDR3)) %>% 
      dplyr::mutate_at(vars(HCDR3), as.character)
    names(HCDR3_libs[[i]])[2] <- paste0("Count_", library_names[i]) %>% 
                                    stringr::str_remove_all(".fastq|.fastq.gz")
    HCDR3s[[i]] <- HCDR3_libs[[i]]$HCDR3
    DNA_libs[[i]] <- as.data.frame(table(DNA)) %>% 
      dplyr::mutate_at(vars(DNA), as.character)
    
    names(DNA_libs[[i]])[2] <- paste0("Count_", library_names[i]) %>% 
                                    stringr::str_remove_all(".fastq|.fastq.gz")
    
    ### for one of each experimental and control, calculate primer templates
    ### from DNA sequence
    if ((i == 1) | (i == length(experimental_fastq_files) + 1)) {
      DNA_libs[[i]]$AA <- Biostrings::DNAStringSet(DNA_libs[[i]]$DNA) %>% 
        Biostrings::translate(if.fuzzy.codon = "X") %>% as.character()
      DNA_libs[[i]]$HCDR3 <- splice_rabbit_HCDR3s(DNA_libs[[i]]$AA, upstream_junction, downstream_junction)
      names(DNA_libs[[i]])[2] <- "counts"
      
      DNA_libs[[i]]$F_primer_template <- DNA_libs[[i]]$DNA %>% 
        substr(split_point - (split_point %% 3) + 1, read_length)
      
      DNA_libs[[i]]$R_primer_template <- DNA_libs[[i]]$DNA %>% 
        substr(1, split_point - (split_point %% 3)) %>% DNAStringSet %>% 
        reverseComplement %>% as.character
      
    }
    
    
  }
 
  
  cat("Finding unique HCDR3s in experimental and control fastq files\n")
  
  ### obtain control matrix - keep sequences common to either all experimental 
  ### or all control outputs
  control_HCDR3s <- Reduce(intersect, HCDR3s[1:length(experimental_fastq_files)])
  experimental_HCDR3s <- Reduce(intersect, HCDR3s[(1+length(experimental_fastq_files)):length(HCDR3s)])
  
  unique_HCDR3s <- c(control_HCDR3s, experimental_HCDR3s) %>% unique
  unique_HCDR3s <- unique_HCDR3s[nchar(unique_HCDR3s) > 1]
  
  ## filter down HCDR3_libs to only those unique HCDR3s
  for (i in 1:length(HCDR3_libs)){
    HCDR3_libs[[i]] <- HCDR3_libs[[i]][HCDR3_libs[[i]]$HCDR3 %in% unique_HCDR3s,]
  }
  cat("Running differential abundance analysis\n")
  ## create count matrix
  matrix <- HCDR3_libs %>% purrr::reduce(dplyr::full_join, by = "HCDR3") 
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
  
  results$mean_non_targ_normalized <- results$mean_non_targ_normalized/
                                      sum(results$mean_non_targ_normalized)
  results$mean_targ_normalized <- results$mean_targ_normalized/
                                  sum(results$mean_targ_normalized)
  
  #Order results by pvalue
  results_ordered <- as.data.frame(results[order(-results$baseMean),])
  
  #Extract DNA sequences from DNA_libs
  DNA_data <- rbind(DNA_libs[[1]], DNA_libs[[length(experimental_fastq_files) + 1]])
  DNA_data <- DNA_data %>% arrange(desc(counts))
  DNA_data <- DNA_data[!duplicated(DNA_data$HCDR3),]
  DNA_data <- DNA_data[!is.na(DNA_data$HCDR3),]
  
  DNA_data <- DNA_data[match(results_ordered$HCDR3, DNA_data$HCDR3),c("DNA", "F_primer_template", "R_primer_template") ]
  
  results_ordered <- cbind(results_ordered, DNA_data)
  
  cat("Analysis complete\n")
  
  return(list(results_ordered, HCDR3_libs))
  
}

