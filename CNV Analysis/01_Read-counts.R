library(cn.mops)
samples <- list.files(path = "/orange/cschnitzler/NS1865-MMartindale_HFHKVDSXY-Lane1-4_DNASeq2020Data/NS1865-MMartindale_HFHKVDSXY-Lane1-4/merged_reads/Platypus_input_files/", pattern = "\\.bam$", full.names = TRUE)
reference.scaffolds <- "/orange/cschnitzler/NS1865-MMartindale_HFHKVDSXY-Lane1-4/NS1865-MMartindale_HFHKVDSXY-Lane1-4/merged_reads/chelonia_mydas_new_referencescaffolds.txt"
output.path <- "/blue/cschnitzler/drewthompson/turtle-fp-wgs/readcounts/"
reference.scaffolds <- read.table(reference.scaffolds, header = T, check.names = F, sep = "\t")
reference.scaffolds <- reference.scaffolds[-c(29:98), ]
reference.scaffolds <- as.character(reference.scaffolds[, "RefSeq-Accn"])
cat("\n Counting reads in specified tumour and host .bam files... \n")
windowsize <- 5000
for (sample in samples) {
  cat(paste("\n", sample, "\n"))
  readcounts <- getReadCountsFromBAM(
    BAMFiles = sample,
    sampleNames = paste(sample),
    WL = windowsize,
    refSeqName = reference.scaffolds
  )
  cat("\n Formatting read count outputs... \n")
  counts.df <- data.frame(
    seqnames = seqnames(readcounts),
    starts = start(readcounts),
    ends = end(readcounts),
    scores = elementMetadata(readcounts)
  )
  counts.mt <- as.matrix(counts.df)
  colnames(counts.mt)[1:4] <- c("SCAFFOLD", "START", "END", "SAMPLE COUNTS")
  cat("\n Writing read count matrix to output... \n")
  tumor <- str_split_fixed(str_split_fixed(sample, "/", 8)[, 8], "_", 4)[, 2]
  write.table(counts.mt,
    file = paste0(output.path, "readcounts_", windowsize, "BP_", sample, ".txt"),
    quote = F, col.names = TRUE, row.names = F, sep = "\t"
  )
  # save(counts.mt, file = paste0(output.path, "readcounts_", windowsize, "BP_yuRKTW1Fdna.Rdata"))
}
