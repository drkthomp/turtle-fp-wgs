library(cn.mops)
library(tidyverse) # for string split

#### 0. Prepare Reference Scaffold ####
reference.scaffolds <- "/orange/cschnitzler/NS1865-MMartindale_HFHKVDSXY-Lane1-4_DNASeq2020Data/NS1865-MMartindale_HFHKVDSXY-Lane1-4/merged_reads/chelonia_mydas_new_referencescaffolds.txt" 
reference.scaffolds <- read.table(reference.scaffolds, header = T, check.names = F, sep = "\t")
reference.scaffolds <- reference.scaffolds[-c(29:98), ]
reference.scaffolds <- as.character(reference.scaffolds[, "RefSeq-Accn"])

#### 1. Indicate Setttings ####
samples <- list.files(path = "/orange/cschnitzler/NS1865-MMartindale_HFHKVDSXY-Lane1-4_DNASeq2020Data/NS1865-MMartindale_HFHKVDSXY-Lane1-4/merged_reads/Platypus_input_files/", pattern = "\\.bam$", full.names = TRUE)
windowsize <- 1000
output.path <- "/blue/cschnitzler/drewthompson/turtle-fp-wgs/CNV Analysis/readcounts/"
cat("\n Counting reads in specified tumour and host .bam files... \n")
for (sample in samples) {
  cat(paste("\n", sample, "\n"))
  tumor <- str_split_fixed(str_split_fixed(sample, "/", 8)[, 8], "_", 4)[, 2]
count_file <-  paste0(output.path, "readcounts_", windowsize, "BP_", tumor, ".txt")
if(file.exists(count_file)){
cat("\n Read count file already exists, skipping. \n")
}
else{
#### 2. Call Read Counts with cn.mops ####
  readcounts <- getReadCountsFromBAM(
    BAMFiles = sample,
    sampleNames = paste(sample),
    WL = windowsize,
    refSeqName = reference.scaffolds
  )
#### 3. Format Read Counts ####
  cat("\n Formatting read count outputs... \n")
  counts.df <- data.frame(
    seqnames = seqnames(readcounts),
    starts = start(readcounts),
    ends = end(readcounts),
    scores = elementMetadata(readcounts)
  )
  counts.mt <- as.matrix(counts.df)
  colnames(counts.mt)[1:4] <- c("SCAFFOLD", "START", "END", "SAMPLE COUNTS")
#### 4. Save Read Counts ####  
cat("\n Writing read count matrix to output... \n")
  write.table(counts.mt, # backup in case the next write errors 
    file = "temp.txt",
    quote = F, col.names = TRUE, row.names = F, sep = "\t"
  )
  write.table(counts.mt,
    file = count_file,
    quote = F, col.names = TRUE, row.names = F, sep = "\t"
  )
  # save(counts.mt, file = paste0(output.path, "readcounts_", windowsize, "BP_yuRKTW1Fdna.Rdata"))
}
}
