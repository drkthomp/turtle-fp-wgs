# BiocManager::install("STRINGdb")
library(cn.mops)
library(DNAcopy)
library(GenomicRanges)
library(tidyverse)
library(janitor)
library(GenomicFeatures)
library(conflicted)
library(BSgenome.Cmydas.NCBI.rCheMyd1)
library(CopyNumberPlots)
library(STRINGdb)

conflict_prefer("filter", "dplyr")

setwd("~/2021-REU/CNV Analysis")
load("partials/sample_data")

txdb <- loadDb("~/2021-REU/CNV Analysis/rCheMyd1.sqlite")
ref <-
  getChromInfoFromNCBI(
    "GCF_015237465.1",
    assembled.molecules.only = TRUE,
    assembly.units = "Primary Assembly"
  )
ref_seq <-
  getChromInfoFromNCBI(
    "GCF_015237465.1",
    as.Seqinfo = TRUE,
    assembled.molecules.only = TRUE,
    assembly.units = "Primary Assembly"
  )
seqnames(ref_seq) <- ref$RefSeqAccn
chr_ref <-
  mutate(ref, CHR = parse_number(SequenceName)) %>% dplyr::select("RefSeqAccn", "CHR")

dna <- extractTranscriptSeqs(BSgenome.Cmydas.NCBI.rCheMyd1, txdb,
  use.names = TRUE
)
aa <- suppressWarnings(translate(dna))

saveSeq <- function(cnv, seqs, save_file) {
  # get the actual protein sequences
  sequences <- getSeq(seqs, as.character(unlist(cnv$tx_name)))

  # and save them for lookup
  writeLines(paste0("> ", names(sequences), "\n", sequences), save_file)
}

annotateTranscripts <- function(cnv, txdb, tum) {
  # by Martin Morgan
  # from https://support.bioconductor.org/p/96427/
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb")) # error checking
  anno <- transcripts(txdb) # get the genes from the ref
  olaps <- findOverlaps(cnv, anno, ignore.strand = TRUE, type = "within") # find the overlaps
  mcols(olaps)$tx_name <- anno$tx_name[subjectHits(olaps)]
  cnv_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  cnv$tx_name <- splitAsList(mcols(olaps)$tx_name, cnv_factor) # and add them to the cnvs
  return(cnv)
}

annotateGenes <- function(cnv, txdb) {
  # by Martin Morgan
  # from https://support.bioconductor.org/p/96427/
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb")) # error checking
  anno <- genes(txdb) # get the genes from the ref
  olaps <- findOverlaps(cnv, anno, ignore.strand = TRUE, type = "within") # find the overlaps
  mcols(olaps)$gene_id <- anno$gene_id[subjectHits(olaps)]
  cnv_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  cnv$gene_id <- splitAsList(mcols(olaps)$gene_id, cnv_factor) # and add them to the cnvs
  return(cnv)
}

cnvDownstream <- function(x) {
  tum <- x[names(x) == "sample_name"]
  normal <- x[names(x) == "normal"]
  tum_name <- str_replace(tum, "_", "-")
  # load in the called data
  load(paste0("partials/normalized/", tum, "-", normal, ".gz"))
  load(paste0("partials/segmented/", tum, "-", normal, ".gz"))


  # cn.mops style chromosome graph
  png(
    file = paste0(
      "plots/",
      tum,
      " FP - chromosome segplot.png"
    ),
    width = 3000,
    height = 1500
  )
  segplot(
    segmented,
    lwd = 1,
    ylim = c(-5, 5),
    plot.type = "s"
  )
  dev.off()

  # cn.mops style segmentation graph
  png(
    file = paste0(
      "plots/",
      tum,
      " FP - segplot.png"
    ),
    width = 3000,
    height = 1500
  )
  segplot(
    segmented,
    lwd = 1,
    ylim = c(-5, 5),
    plot.type = "w"
  )
  dev.off()


  # add annotations
  print(paste0("=== ANNOTATING FOR TUMOR ", tum))
  cnv <- cnvr(segmented)
  cnv <- annotateGenes(cnv, txdb)
  cnv <- annotateTranscripts(cnv, txdb, tum)
  write(cnv$gene_id %>% unlist(), paste0("data/", tum, "_cnvs.txt"))
  write.table(
    cnv,
    file = paste0("data/", tum, "_cnvs.tsv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  saveSeq(nv, aa, paste0("data/", tum, "_protein_sequences.txt"))
  saveSeq(nv, dna, paste0("data/", tum, "_dna_sequences.txt"))

  print(paste0("=== PLOTTING FOR TUMOR ", tum))
  title <- paste(
    # turtle tissue FP - segmented CNVs

    str_replace(x[names(x) == "turtle.x"], "_kidney", ""),
    filter(sample_data, sample_name == tum_name)$tumor_sample_location,
    "FP - segmented CNVs"
  )
  file_title <- paste0( # tum_normal FP - segmented CNVs(.pdf)
    tum,
    "_",
    normal,
    " FP - segmented CNVs"
  )

  # karyoPlot(segmented,title,file_title)
  # Plotting
  # nice seg plot takes a while, sometimes want to not rerun

  nice.seg.plot(
    x.segmented = segmented,
    x.normalised.counts = norm,
    indicateChrom = FALSE,
    # don't indicate the chromosomes (busy)
    colorChrom = FALSE,
    # def don't color the chromosomes (busy, no legend)
    annotate = TRUE,
    # do add small annotations
    title = title,
    file_title = file_title
  )
}


run <- sample_data
# and go and do it for everything!
apply(run, 1, cnvDownstream)
