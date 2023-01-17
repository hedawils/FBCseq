library(tidyverse)
library(FBCseq)

U266_results <- differential_HCDR3(c("314_4_1.fastq", "314_5_1.fastq", "314_6_1.fastq"), c("314_7_1.fastq", "314_8_1.fastq", "314_9_1.fastq"), primer = "GCCCTTGGTGGAGGC")
H929_results <- differential_HCDR3(c("314_1_1.fastq", "314_2_1.fastq", "314_3_1.fastq"), c("314_7_1.fastq", "314_8_1.fastq", "314_9_1.fastq"), primer = "GCCCTTGGTGGAGGC")
HEK_FBC_results <- differential_HCDR3(c("322_4_1_Total.fastq", "322_5_1_Total.fastq", "322_6_1_Total.fastq"), c("322_1_1_Total.fastq", "322_2_1_Total.fastq", "322_3_1_Total.fastq"), primer = "GCCCTTGGTGGAGGC")
HEK_conv_results <- differential_HCDR3(c("322_10_1_Total.fastq", "322_11_1_Total.fastq", "322_12_1_Total.fastq"), c("322_7_1_Total.fastq", "322_8_1_Total.fastq", "322_9_1_Total.fastq"), primer = "GCCCTTGGTGGAGGC")


write_outputs <- function(results, title){
  HCDR3_libs <- results[[2]]
  results_ordered <- results[[1]]
  
  
  for (lib in 1:length(HCDR3_libs)){
    data <- HCDR3_libs[[lib]]
    name <- stringr::str_remove(names(data)[2], "Count_")
    names(data)[2] <- "count"
    data <- data %>% dplyr::arrange(desc(count))
    write_tsv(data, paste0(name, "_counts.txt"), quote = "none")
    cat(paste0("wrote ", name, "\n"))
  }

  write_tsv(results_ordered, paste0(title,"_processed_output.txt"))
  }

write_outputs(U266_results, title = "U266_v_PBMCs")
write_outputs(H929_results, title = "H929_v_PBMCs")
write_outputs(HEK_FBC_results, title = "HEK293_ROR1pos_v_ROR1neg_FBC")
write_outputs(HEK_conv_results , title = "HEK293_ROR1pos_v_ROR1neg_conventional")
