####################################################
#        Sea Turtle FP Tumour CNV Analysis         #
#                                                  #
# 26.01.2020                                       #
# TCG Cambridge, Department of Veterinary Medicine #
# mrs72@cam.ac.uk                                  #
####################################################

# Libraries and Dependencies
library(cn.mops)
library(DNAcopy)
library(GenomicRanges)
library(tidyverse)
library(janitor)
library(GenomicFeatures)
library(conflicted)
library(BSgenome.Cmydas.NCBI.rCheMyd1)
library(CopyNumberPlots)

conflict_prefer("desc", "IRanges")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

load("partials/sample_data")

##########################
# 1. Load Data Imports
##########################

dataloc <- "readcounts/" # note where the readcounts are
bp <- 5000 # then bp we're using
# and the final file name before the tumor:
file_prefix <- paste0(dataloc, "readcounts_", bp, "BP_")

# list of valid tumor files to be filled


# get relevant reference data
txdb <- loadDb("~/2021-REU/CNV Analysis/rCheMyd1.sqlite")
ref <- getChromInfoFromNCBI("GCF_015237465.1", assembled.molecules.only = TRUE, assembly.units = "Primary Assembly")
ref_seq <- getChromInfoFromNCBI("GCF_015237465.1", as.Seqinfo = TRUE, assembled.molecules.only = TRUE, assembly.units = "Primary Assembly")
seqnames(ref_seq) <- ref$RefSeqAccn
chr_ref <-
  mutate(ref, CHR = parse_number(SequenceName)) %>% select("RefSeqAccn", "CHR")

# only doing for yucca
normal_names <- unique(grep("yu", sample_data$normal, fixed = TRUE, value = TRUE))
tumor_names <- unique(grep("yu", sample_data$sample_name, fixed = TRUE, value = TRUE))

normals <- read_in_all(normal_names)
tumors <- read_in_all(tumor_names)



read_in_counts <- function(x) {
  result <- read.table(
    x,
    header = T,
    strip.white = T,
    check.names = F,
    sep = "\t"
  )
  type_convert(result)

  # dropping the mitochondria mapping here
  result <-
    result %>%
    filter(SCAFFOLD %in% ref$RefSeqAccn) %>%
    mutate(SCAFFOLD = as.factor(SCAFFOLD))
  return(result)
}
read_in_all <- function(names) {
  files <- paste0(paste0(dataloc, "readcounts_", bp, "BP_"), names, ".txt")
  counts <- lapply(files, read_in_counts)
  names(counts) <- normal_names
  result <- counts %>% purrr::reduce(left_join,
    by = c("SCAFFOLD", "START", "END")
  )
  colnames(result)[1:length(names) + 3] <- names
  final <-
    makeGRangesFromDataFrame(result,
      seqnames.field = "SCAFFOLD",
      keep.extra.columns = TRUE
    )
  seqinfo(final) <- ref_seq # update the seqinfo with the relevant genome info
  return(final)
}


segmented <- referencecn.mops(cases = tumors, controls = normals, normType = "quant", segAlgorithm = "fast")

counts <- calcIntegerCopyNumbers(segmented)
segplot(counts, lwd = 1, ylim = c(-5, 5), plot.type = "w", pt.cols = c("#999999", "#909090"), pt.cex = 1)
