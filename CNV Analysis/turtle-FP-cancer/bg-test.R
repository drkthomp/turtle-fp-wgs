### Download GCA_000317375.1_MicOch1.0_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/317/375/GCA_000317375.1_MicOch1.0/
library(tidyverse)
library(BSgenome.Cmydas.NCBI.rCheMyd1)
library(CopyNumberPlots)

kp <- plotKaryotype(genome = "BSgenome.Cmydas.NCBI.rCheMyd1", method = "NCBI", chromosomes = 1:28)
cn.calls <- loadCopyNumberCalls(cnv, cn.col = "TUMOUR")

# adjust to numeric chromosome
seqlevels(cn.calls) <- as.character(1:length(levels(seqnames(cn.calls))))
genome(cn.calls) <- rep("rCheMyd1.pri", length(genome(cn.calls)))
# and change cn calls to numeric
cn.calls$cn <- parse_number(cn.calls$cn)
plotCopyNumberCalls(kp, cn.calls = cn.calls)


chr_ref <- as.data.frame(SequenceName = ref$SequenceName %>% str_remove_all("SUPER_"), RefSeqAccn = ref$RefSeqAccn)

unmatch(ref[1, "RefSeqAccn"], ref)
