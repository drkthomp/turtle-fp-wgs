
library(cn.mops)
library(GenomicRanges)
library(tidyverse)
library(janitor)
library(GenomicFeatures)
library(BSgenome.Cmydas.NCBI.rCheMyd1) # this was custom generated and should be provided in the GitHub

#### 0. Define Options ####
saveloc <- "02-Annotation/"
reuse <- TRUE

#### 1. Load Genome Data ####
sample_data <- read.csv("~/2021-REU/CNV Analysis/partials/sample_data.tsv", sep="")

dna <- extractTranscriptSeqs(BSgenome.Cmydas.NCBI.rCheMyd1, txdb,
                      use.names=TRUE)
aa <- suppressWarnings(translate(dna))

txdb <- loadDb("rCheMyd1.sqlite") # generated from rCheMyd1
ref_file <- "partials/ref.Rdata"
if (file.exists(ref_file)) {
  load(ref_file)
} else {
  ref <- getChromInfoFromNCBI("GCF_015237465.1",
                              assembled.molecules.only = TRUE,
                              assembly.units = "Primary Assembly"
  )
  save(ref, file = ref_file)
}

ref_seq_file <- "partials/ref_seq.Rdata"
if (file.exists(ref_seq_file)) {
  load(ref_seq_file)
} else {
  ref_seq <- getChromInfoFromNCBI("GCF_015237465.1",
                                  as.Seqinfo = TRUE,
                                  assembled.molecules.only = TRUE,
                                  assembly.units = "Primary Assembly"
  )
  seqnames(ref_seq) <- ref$RefSeqAccn
  save(ref_seq, file = ref_seq_file)
}

protein_table <- makeGRangesFromDataFrame(read.delim("~/2021-REU/CNV Analysis/ProteinTable_13308_1483792_Cm_NEW_NCBI.txt", quote="") %>% clean_names() %>% rename(seqnames=accession) %>% dplyr::select(-x_name),keep.extra.columns=TRUE)

#### 2. Annotation Functions ####
# The functions to annotate both for transcripts and genes:
saveSeq <- function(cnv, seqs, save_file){
  # get the actual protein sequences
  sequences <- getSeq(seqs, as.character(unlist(cnv$tx_name)))

  # and save them for lookup
  writeLines(paste0("> ", names(sequences), "\n", sequences), save_file)

}

annotateAll <- function(cnv, protein_table){
  # Modifed from overlap by Martin Morgan
  # from https://support.bioconductor.org/p/96427/
  stopifnot(is(cnv, "GRanges"), is(protein_table, "GRanges"))
  olaps = findOverlaps(cnv, protein_table, ignore.strand = TRUE, type = type)
  mcols(olaps) <-  as.data.frame(protein_table[subjectHits(olaps)]) %>% mutate(protein_name_short = str_squish( str_replace_all(str_remove_all(protein_name, '-like|isoform|LOW QUALITY PROTEIN|"'), "F10", "F")))
  names(mcols(olaps))[1:5] <- paste0("gene_",names(mcols(olaps))[1:5])
  cnv_factor = factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  mcols(cnv) <- cbind(mcols(cnv), apply(mcols(olaps), MARGIN=2, function (x){splitAsList(x, cnv_factor)}))
  return(cnv)
}

annotateTranscripts <- function(cnv, txdb, tum)
{
  # by Martin Morgan
  # from https://support.bioconductor.org/p/96427/
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb")) # error checking
  anno = transcripts(txdb) # get the genes from the ref
  olaps = findOverlaps(cnv, anno, ignore.strand = TRUE, type = type) # find the overlaps
  mcols(olaps)$tx_name = anno$tx_name[subjectHits(olaps)]
  cnv_factor = factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  cnv$tx_name = splitAsList(mcols(olaps)$tx_name, cnv_factor) # and add them to the cnvs
  return(cnv)
}

annotateGenes <- function(cnv, txdb)
{
  # by Martin Morgan
  # from https://support.bioconductor.org/p/96427/
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb")) # error checking
  anno = genes(txdb) # get the genes from the ref
  olaps = findOverlaps(cnv, anno, ignore.strand = TRUE, type = type) # find the overlaps
  mcols(olaps)$gene_id = anno$gene_id[subjectHits(olaps)]
  cnv_factor = factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  cnv$gene_id = splitAsList(mcols(olaps)$gene_id, cnv_factor) # and add them to the cnvs
  return(cnv)
}


# filtering out Yucca's messy chromosomes
throwout <- 'yu'
throwout_percent <- .05

toThrowout <- function(start){
  loc <- as.numeric(start)/as.numeric(ref[ref$RefSeqAccn == "NC_051241.1",]$SequenceLength)
  return(loc<throwout_percent | loc > (1-throwout_percent))
}

# load all the possible segmentation files
seg_files <- list.files(path="partials/segmented",recursive = TRUE,full.names = FALSE)
options <-  data.frame(bp=str_split_fixed(seg_files, "/",3)[,1],normType=str_split_fixed(seg_files, "/",3)[,2],
                       sizeFactor=str_split_fixed(str_split_fixed(seg_files, "/",3)[,3], "-",4)[,1],
                       minReadCount=as.numeric(str_split_fixed(str_split_fixed(seg_files, "/",3)[,3], "-",4)[,2]),
                       tum=str_split_fixed(str_split_fixed(seg_files, "/",3)[,3],  "-",4)[,3],
                       normal=str_remove_all(str_split_fixed(str_split_fixed(seg_files, "/",3)[,3], "-",4)[,4], ".gz"),
                       files=seg_files) %>% na.omit()

# makes groups for each possible option run
options <- options %>% group_by(bp, normType, sizeFactor, minReadCount)
message("Iterating through ", count(group_keys(options)), " options.")

apply(group_keys(options),MARGIN=1,function(x){

  seg_files <- (options %>% filter(bp==x[1], normType==x[2],
                                   sizeFactor==x[3],minReadCount==x[4])) # get all the files matching the conditions
  seg_data <- lapply(seg_files$files, function(x){ # then load them
    load(paste0("partials/segmented/", x))
    return(segmented)
    })
  names(seg_data) <- seg_files$tum # name them by the tumor sample

  # ==== FILTRATION
  # 1. Filter genome areas with normalized read counts
  # significantly below or above the sample-specific normalized coverage
  # distribution across all tumours and normals



  message(paste(x,collapse=" "))
})





cnvDownstream <- function(x){
  tryCatch({
  #bp <- str_split_fixed(x, "/",3)[1]
  #normType <-  str_split_fixed(x, "/",3)[2]
  #sizeFactor <-  str_split_fixed(str_split_fixed(x, "/",3)[3], "-",4)[1]
  #minReadCount <- str_split_fixed(str_split_fixed(x, "/",3)[3], "-",4)[2]
  tum <- str_split_fixed(str_split_fixed(x, "/",3)[3], "-",4)[3]
  normal <- str_split_fixed(str_split_fixed(x, "/",3)[3], "-",4)[4]

  for(loc in c(saveloc, "partials/cnv")){
    # ensure that the directories are created
    dir.create(paste(loc,type,sep="/"))
    dir.create(paste(loc,type, str_split_fixed(x, "/",3)[1], sep="/"))
    dir.create(paste(loc,type, str_split_fixed(x, "/",3)[1],str_split_fixed(x, "/",3)[2], sep="/"))
  }

# define where everything will be saved
cnv_file <- paste0("partials/cnv/",type,"/", x)
  gene_file <- paste0(saveloc,type,"/", str_remove(x, ".gz"), "_cnvs.txt")
  table_file <- paste0(saveloc,type,"/", str_remove(x, ".gz"), "_cnvs.tsv")
  protein_file <- paste0(saveloc,type,"/", str_remove(x, ".gz"), "_protein_sequences.txt")
  dna_file <- paste0(saveloc,type,"/", str_remove(x, ".gz"), "_dna_sequences.txt")


  if(reuse & file.exists(cnv_file) & file.exists(gene_file) & file.exists(table_file) & file.exists(protein_file) & file.exists(dna_file)){
    message("All files exist for ", x) # don't redo anything if everything exists, just move on
  }
  else{
  if(reuse & file.exists(cnv_file)){
    message("Loading annotate cnv file ", cnv_file)
    load(cnv_file)
  }
  else{
  print(paste0("=== ANNOTATING FOR TUMOR ", tum))
  # load in the called data
  message("Loading segmentation file ", x)
  load(paste0("partials/segmented/", x))

  # add annotations

  cnv <- cnvr(segmented) # should be loaded in
  #cnv <- annotateGenes(cnv, txdb)
  cnv <- annotateTranscripts(cnv, txdb,tum)
  cnv <- annotateAll(cnv, protein_table)

  # and then save the annotations
  message("Saving annotation file ", cnv_file)
  save(cnv, file=cnv_file,compress=TRUE)
  }


  cnv <- as.data.frame(cnv) %>% mutate(throwout = toThrowout(start) | toThrowout(end))
  if(grepl(throwout, x)){
    print(paste("Threw out", count(cnv %>% filter(throwout)), "regions"))
    print(cnv%>% filter(throwout))
    cnv <- cnv %>% filter(!throwout)
  }

  if(!(reuse & file.exists(gene_file))){
    message("saving gene file at ", gene_file)
    write(cnv$gene_id %>% unlist(), gene_file) # only gene information
  }
  if(!(reuse & file.exists(table_file))){
        message("saving cnv table at ", table_file)
  write.table( # more full data, allows recreation
    cnv,
    file = table_file,
    quote = FALSE,
    sep = '\t',
    row.names = FALSE
  )
  }
  # and then save the protein and dna sequences based on the annotated genes
  if(!(reuse & file.exists(protein_file))){
    message("Saving cnv protein sequences at ", protein_file)
    saveSeq(cnv,aa,protein_file)
  }
  if(!(reuse & file.exists(dna_file))){
    message("Saving cnv dna sequences at ", dna_file)
    saveSeq(cnv,dna, dna_file)
  }
  }}, error=function(cond){
    message("CNV annotation failed for ", x)
    message(cond, "\n")
    })
}

for(type in c('within', 'any')){
lapply(seg_files, cnvDownstream)
}



# do some cleanup
rm(dna)
rm(aa)
gc()
```



