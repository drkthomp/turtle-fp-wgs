setwd("~/2021-REU/bsgenome")

# download fna.gz from ftp
library(tidyverse)
library(Biostrings)

library(DNAcopy)
library(GenomicRanges)
library(GenomicFeatures)
dna <- readDNAStringSet("~/2021-REU/somatypus/input/GCF_015237465.1_rCheMyd1.pri_genomic.fna")

### Check the seqnames.
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
#library(GenomeInfoDb)

#chrominfo <- getChromInfoFromNCBI("GCF_015237465.1")%>% filter(RefSeqAccn %in% current_GenBankAccn)
#ref_seq <- getChromInfoFromNCBI("GCF_015237465.1",as.Seqinfo=TRUE)

# re ordered just in case
#chrom_ordered <- chrominfo[match(current_RefSeqAccn, chrominfo$RefSeqAccn), ]        # Reorder data frame

### Rename the sequences numerically?
names(dna) <- current_RefSeqAccn


### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "rCheMyd1.pri.sorted.2bit")

library(BSgenome)
forgeBSgenomeDataPkg("BSgenome.Cmydas.NCBI.rCheMyd1-seed",destdir="output")


