library(BSgenome.Cmydas.NCBI.rCheMyd1)
library(GenomicRanges)

cds_seqs <- extractTranscriptSeqs(BSgenome.Cmydas.NCBI.rCheMyd1, cdsBy(txdb, by = "tx", use.names = TRUE, ))
translated <- translate(cds_seqs)

seqlevels(cnv) <- seqlevels(BSgenome.Cmydas.NCBI.rCheMyd1)
seqinfo(cnv) <- seqinfo(BSgenome.Cmydas.NCBI.rCheMyd1)

test <- transcriptsBy(txdb, by = "gene")

cnv_trim <- trim(cnv)

# olaps <- extractAt(translated, ranges(cnv_trim))
olaps <- findOverlaps(cnv, test, ignore.strand = TRUE, type = "within") # find the overlaps
mcols(olaps)$gene_id <- anno$gene_id[subjectHits(olaps)]
cnv_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
cnv$gene_id <- splitAsList(mcols(olaps)$gene_id, cnv_factor) # and add them to the cnvs
return(cnv)
