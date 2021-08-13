library(cn.mops)
library(GenomicRanges)
library(tidyverse)
library(janitor)
library(GenomicFeatures)
library(BSgenome.Cmydas.NCBI.rCheMyd1)

# Calls CNVs and performs basic filtering of the read counts

# setwd("~/2021-REU/CNV Analysis")
bp <- 5000
#### 0. Define Settings ####
reuse <- TRUE # whether to reuse files
dataloc <- "readcounts/" # note where the readcounts are
# and the final file name before the tumor:
file_prefix <- paste0(dataloc, "readcounts_", bp, "BP_")
sample_data <- read.csv("partials/sample_data.tsv", sep = "")
protein_table <- makeGRangesFromDataFrame(read.delim("~/2021-REU/CNV Analysis/ProteinTable_13308_1483792_Cm_NEW_NCBI.txt", quote = "") %>% clean_names() %>% dplyr::rename("seqnames" = "accession") %>% dplyr::select(-x_name), keep.extra.columns = TRUE)

annoloc <- "02-Annotation/"

#### 2. Load Reference Data ####
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

dna <- extractTranscriptSeqs(BSgenome.Cmydas.NCBI.rCheMyd1, txdb,
                             use.names = TRUE
)
aa <- suppressWarnings(translate(dna))

#### 3. Define Functions ####
##### 3.1 Count Functions #####

###### 3.1.1 Read in the counts based on a row from sample_data ######
readCounts <- function(x) {
  tumor_file <-
    paste0(file_prefix, x[names(x) == "sample_name"], ".txt")
  normal_file <-
    paste0(file_prefix, x[names(x) == "normal"], ".txt")

  message(paste("Reading in", tumor_file, "and", normal_file))

  # make sure that we've actually downloaded the files
  if (file.exists(tumor_file) && file.exists(normal_file)) {
    tumor <-
      read.table(
        tumor_file,
        header = T,
        strip.white = T,
        check.names = F,
        sep = "\t"
      )
    normal <-
      read.table(
        normal_file,
        header = T,
        strip.white = T,
        check.names = F,
        sep = "\t"
      )

    # and now append host
    result <- left_join(
      tumor,
      normal,
      suffix = c("_TUMOR", "_NORMAL"),
      by = c("SCAFFOLD", "START", "END")
    )
    colnames(result)[4:5] <- c("TUMOUR", "HOST")
    type_convert(result)

    # dropping the mitochondria mapping here
    result <-
      result %>%
      filter(SCAFFOLD %in% ref$RefSeqAccn) %>%
      mutate(SCAFFOLD = as.factor(SCAFFOLD))

    final <-
      makeGRangesFromDataFrame(result,
        seqnames.field = "SCAFFOLD",
        keep.extra.columns = TRUE
      )
    seqinfo(final) <- ref_seq # update the seqinfo with the relevant genome info
    return(final)
  }
}

###### 3.1.2 Normalize the counts based on a row from sample_data file ######
normCounts <- function(x, bp = 5000, normType = "poisson", sizeFactor = "mean") {
  tum <- x[names(x) == "sample_name"]
  gr <- readCounts(x)
  if (!is.null(gr)) { # ensure that we can proceed
    message(paste("=== NORMALIZING FOR TUMOR", tum))
    norm <- normalizeGenome(gr, sizeFactor = sizeFactor, normType = normType)
    return(norm)
  }
  return(NULL)
}


###### 3.1.3 Filter Counts Across Samples ######
## Filter or mask genome areas with normalised read counts
##  significantly below or above the sample-specific normalised
##  coverage distribution across all tumours and normals

filterCountsAcrossSamplesHelper <- function(norm_data_chrom, filter = 0.01) {
  ## Shouldn't be accessed outside of filterCountsAcrossSamples
  lower_norm_range <- norm_data_chrom %>% summarise(across(where(is.numeric), quantile, probs = c(filter)))
  upper_norm_range <- norm_data_chrom %>% summarise(across(where(is.numeric), quantile, probs = c(1 - filter)))
  norm_data_chrom <- norm_data_chrom %>% mutate(pass = case_when(
    if_all(where(is.numeric), ~ .x < lower_norm_range[[cur_column()]]) ~ FALSE,
    if_all(where(is.numeric), ~ .x > upper_norm_range[[cur_column()]]) ~ FALSE,
    TRUE ~ TRUE
  ))
  return(norm_data_chrom$pass)
}

filterCountsAcrossSamples <- function(norm_grange, filter = 0.01,
                                      byChrom = TRUE, savefile = NULL, reuse = TRUE) {
  if (reuse & file.exists(savefile)) {
    load(savefile)
  } else {
    norm_data <- lapply(names(norm_grange), function(x) {
      grange <- norm_grange[[x]]
      if (!is.null(grange)) { # guard against there being nothing, in which case can't call mcols
        names(mcols(grange)) <- c(x, sample_data[sample_data$sample_name == x, ]$normal)
        return(cbind(as.data.frame(seqnames(grange)), as.data.frame(mcols(grange))))
      }
      return(NULL)
    })
    if (is.null(norm_data)) { # exit early in case of complete null
      return(NULL)
    } # otherwise continue
    norm_data <- norm_data %>% purrr::reduce(cbind)
    norm_data <- norm_data[, !duplicated(colnames(norm_data))] %>%
      rename(chr = value)
    if (byChrom) {
      pass <- lapply(norm_data %>%
        group_by(chr) %>%
        group_split(), filterCountsAcrossSamplesHelper, filter = filter) %>% unlist()
    } else {
      pass <- filterCountsAcrossSamplesHelper(norm_data, filter = filter)
    }
    if (!is.null(savefile)) {
      save(pass, file = savefile)
    }
    message(paste0("Dropping ", count(norm_data %>% filter(!pass)), "/", count(norm_data), " (", round((count(norm_data %>% filter(!pass)) / count(norm_data)) * 100, 2), "%) read count regions.\n"))
  }

  norm_grange <- lapply(norm_grange, function(x) {
    x[pass]
  })

  return(norm_grange)
}

###### 3.1.4. Read In Counts ######
prepCounts <- function(sample_data, bp = 5000, normType = "poisson",
                       sizeFactor = "mean", filter = 0.01,
                       byChrom = TRUE, reuse = TRUE) {
  base_loc <- paste0(
    "partials/normalized/", bp, "/", normType, "/",
    sizeFactor, "-", filter, "-", byChrom
  )
  raw_file <- paste0(base_loc, "_RAW.gz")
  pass_file <- paste0(base_loc, "_PASS.gz")
  count_file <- paste0(base_loc, "_COUNTS.gz")
  if (reuse & file.exists(raw_file) & file.exists(pass_file) & file.exists(count_file)) {
    message("Loading counts at ", count_file)
    load(count_file)
  } else {
    # Read in and normalize the counts
    if (file.exists(raw_file)) {
      message("Raw file exists at ", raw_file)
      load(raw_file)
    } else {
      norm_grange <- apply(sample_data, 1, normCounts,
        bp = bp,
        normType = normType, sizeFactor = sizeFactor
      )
      names(norm_grange) <- sample_data$sample_name
      save(norm_grange, file = raw_file)
      message("Saving raw normalization counts at ", raw_file)
    }

    norm_grange_filtered <- filterCountsAcrossSamples(norm_grange,
      filter = filter, byChrom = byChrom,
      savefile = pass_file, reuse = reuse
    )
    message("Saving filtered normalized data at", count_file)
    save(norm_grange_filtered, file = count_file)
  }
  return(norm_grange_filtered)
}

##### 3.2 Segmentation Function #####
###### 3.2.1 Calling Segmentation ######
callCNVfromNorm <- function(name, norm_full, minReadCount = 5, parallel = 2) {
  message("Starting segmentation call for ", name)
  file <- paste0("partials/segmented/temp/", minReadCount, "-", name, ".gz")
  if (file.exists(file)) {
    message("Temporary seg call exists at ", file)
    load(file)
  } else {
    norm <- norm_full[[name]]
    segmented <- referencecn.mops(
      cases = norm[, "TUMOUR"],
      controls = norm[, "HOST"],
      minReadCount = minReadCount,
      segAlgorithm = "fast",
      norm = 0, # since we're inputting already normalized read counts
      parallel = parallel
    )
    segmented <- calcIntegerCopyNumbers(segmented)
    message("Saved temporary file at ", file)
    save(segmented, file = file)
  }
  return(segmented)
}
###### 3.2.2 Filtering ######
###### TODO make functional. currently obselete
filterCNVHelper <- function(cnv_single, quant_single) {
  if (length(cnv_single) > 0) { # make sure there are CNVs in the chromosome
    pass <- lapply(cnv_single, function(x) {
      return(x < quant_single[1] | x > quant_single[2])
    }) %>% unlist()
    return(pass)
  }
}
####### 3.2.2.1 By Signifigance
# For the remainder list of CNV segments, go through each segment/tumour pair
# and remove the call if a tumour’s normalised median read count within
# it is not significantly above or below its background distribution
# across the whole genome/chromosome
filterCNVSignifigance <- function(segmented, filter = 0.05, byChrom = TRUE) {
  cnvs <- cnvs(segmented)
  if (byChrom) {
    cnvs <- splitAsList(cnvs, seqnames(cnvs))
    quantile <- lapply(splitAsList(normalizedData(segmented), seqnames(normalizedData(segmented))), function(x) {
      return(quantile(x$TUMOUR, probs = c(.5 - filter, .5 + filter)))
    })

    pass <- lapply(names(cnvs), function(x, cnvs, quantile) {
      median_data <- mcols(cnvs[[x]])$median
      if (length(median_data) > 0) {
        return(filterCNVHelper(median_data, quantile[[x]]))
      }
    }, cnvs = cnvs, quantile = quantile) %>% unlist()
  } else {
    quantile <- quantile(normalizedData(segmented)$TUMOUR, probs = c(.5 - filter, .5 + filter))
    pass <- filterCNVHelper(mcols(cnvs)$median, quantile)
  }
  return(cnvs(segmented)[pass])
}

####### 3.2.2.2 By Shared Deviation
# remove the call if both tumour AND matched normal(s)
# deviate significantly from their background distributions
# (in this case it’s likely not a tumour-specific event)
filterCNVTumour <- function(norm, cnvs, filter = 0.05) {
  if (byChrom) {
    cnvs <- splitAsList(cnvs, seqnames(cnvs))
    norm <- splitAsList(norm, seqnames(norm))

    # TODO fix repeats here
    tumour_quantile <- lapply(norm, function(x) {
      return(quantile(x$TUMOUR, probs = c(filter, 1 - filter)))
    })
    host_quantile <- lapply(norm, function(x) {
      return(quantile(x$HOST, probs = c(filter, 1 - filter)))
    })

    tumour_pass <- lapply(names(norm), function(x, norm, quantile) {
      if (length(norm[[x]]) > 0) {
        return(filterCNVHelper(norm[[x]]$TUMOUR, quantile[[x]]))
      }
    }, norm = norm, quantile = tumour_quantile) %>% unlist()
    host_pass <- lapply(names(norm), function(x, norm, quantile) {
      if (length(norm[[x]]) > 0) {
        return(filterCNVHelper(norm[[x]]$HOST, quantile[[x]]))
      }
    }, norm = norm, quantile = host_quantile) %>% unlist()
  } else {
    quantile <- quantile(normalizedData(segmented)$TUMOUR, probs = c(.5 - filter, .5 + filter))
    pass <- filterCNVHelper(mcols(cnvs)$median, quantile)
  }
  return(cnvs(segmented)[pass])
}


###### 3.2.3 Visualization ######
plotClear <- function(){
  if((!is.null(dev.list())) & length(dev.list()) > 0){
    for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
      tryCatch(dev.off(),error=message)
    }}
}
plotSingleCNVs <- function(which, tum, segmented_single, base_loc, margin = 10) {
  tryCatch({
  file <- paste0(base_loc, paste(tum, cnvr(segmented_single)[which], mcols(cnvr(segmented_single))[which, ], sep = "-"), ".png")
  },error=function(cond){
    message(paste("Failed plotting CNV", which, tum))
    message(cond)
    return(NA) # exit before trying to start plotting
  })
  message(paste("Plotting single CNV", which,"at", file))
  png(filename = file, width = 1000, height = 965)
  tryCatch(plot(segmented_single, which = which, margin = c(margin, margin)),error=print)
  dev.off()
  plotClear()

}

plotCNV <- function(tum, segmented, base_loc, lim = 5) {
  chromplot_file <- paste0(base_loc, tum, " FP - chromosome segplot.png")
  segplot_file <- paste0(base_loc, tum, " FP - segplot.pdf")
  size_file <- paste0(base_loc, tum, " FP - Size Distribution.png")
  type_file <- paste0(base_loc, tum, " FP - Type Distribution.png")
  if (!(reuse && file.exists(chromplot_file))) {
    message("Generating chromosome plot at ", chromplot_file)

    png(
      filename = chromplot_file,
      width = 3000,
      height = 1500
    )
    tryCatch(
      segplot(
        segmented[[tum]],
        ylim = c(-lim, lim),
        plot.type = "s"
      ),error=message
    )
    dev.off()
    plotClear()
  }
  if (!(reuse && file.exists(segplot_file))) {
    message("Generating chromosome plot at ", segplot_file)

    pdf(file=segplot_file,
      width = 18,
      height
      = 9
    )
    tryCatch(
      segplot(
        segmented[[tum]],
        ylim = c(-lim, lim),
        plot.type = "w",
        pt.cols = c("#333333", "#000000")
      ),error=message
    )
    dev.off()
    plotClear()
  }
  data <- data.frame(width=width(cnvr(segmented[[tum]])),type=mcols(cnvr(segmented[[tum]]))$TUMOUR)
  if(!(reuse && file.exists(size_file))){
    ggplot(data, aes(x=width,fill=type)) + geom_histogram(bins=100) + labs(title=paste(tum, "Width Distribution"),x="Width of CNV Region") + theme_minimal()
    ggsave(filename = size_file)
  }
  if(!(reuse && file.exists(type_file))){
    ggplot(data, aes(x=parse_number(type))) + geom_histogram(stat="count",binwidth=1) + labs(title=paste(tum, "Type Distribution"),x="CNV Type") + theme_minimal()
    ggsave(filename = type_file)
  }

  # just in case
  plotClear()
 # if(length(cnvr(segmented[[tum]])) > 0){
 #   lapply(1:length(cnvr(segmented[[tum]])), plotSingleCNVs, tum = tum, segmented = segmented[[tum]], base_loc = base_loc)
 # }
 # plotClear()
}


##### 3.3 Annotation Functions #####

# The functions to annotate both for transcripts and genes:
saveSeq <- function(cnv, seqs, save_file) {
  # get the actual protein sequences
  sequences <- getSeq(seqs, as.character(unlist(cnv$tx_name)))

  # and save them for lookup
  writeLines(paste0("> ", names(sequences), "\n", sequences), save_file)
}

annotateAll <- function(cnv, protein_table, type = "any") {
  # Modifed from overlap by Martin Morgan
  # from https://support.bioconductor.org/p/96427/
  stopifnot(is(cnv, "GRanges"), is(protein_table, "GRanges")) # required for find overlaps
  olaps <- GenomicRanges::findOverlaps(cnv, protein_table, ignore.strand = TRUE, type = type)
  mcols(olaps) <- as.data.frame(protein_table[subjectHits(olaps)]) %>% mutate(protein_name_short = str_squish(str_replace_all(str_remove_all(protein_name, '-like|isoform|LOW QUALITY PROTEIN|"'), "F10", "F")))
  names(mcols(olaps))[1:5] <- paste0("gene_", names(mcols(olaps))[1:5])
  cnv_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  mcols(cnv) <- cbind(mcols(cnv), apply(mcols(olaps), MARGIN = 2, function(x) {
    splitAsList(x, cnv_factor)
  }))
  return(cnv)
}

annotateTranscripts <- function(cnv, txdb, tum, type = "any") {
  # by Martin Morgan
  # from https://support.bioconductor.org/p/96427/
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb")) # error checking
  anno <- transcripts(txdb) # get the genes from the ref
  olaps <- findOverlaps(cnv, anno, ignore.strand = TRUE, type = type) # find the overlaps
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
  olaps <- findOverlaps(cnv, anno, ignore.strand = TRUE, type = type) # find the overlaps
  mcols(olaps)$gene_id <- anno$gene_id[subjectHits(olaps)]
  cnv_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  cnv$gene_id <- splitAsList(mcols(olaps)$gene_id, cnv_factor) # and add them to the cnvs
  return(cnv)
}

throwout <- "yu"
throwout_percent <- .05

toThrowout <- function(start) {
  loc <- as.numeric(start) / as.numeric(ref[ref$RefSeqAccn == "NC_051241.1", ]$SequenceLength)
  return(loc < throwout_percent | loc > (1 - throwout_percent))
}


cnvDownstream <- function(tum, src, segmented, type = "any") {
  # bp <- str_split_fixed(x, "/",3)[1]
  # normType <-  str_split_fixed(x, "/",3)[2]
  # sizeFactor <-  str_split_fixed(str_split_fixed(x, "/",3)[3], "-",4)[1]
  # minReadCount <- str_split_fixed(str_split_fixed(x, "/",3)[3], "-",4)[2]


  for (loc in c(annoloc, "partials/cnv")) {
    # ensure that the directories are created
    dir.create(paste(loc, type, sep = "/"))
    dir.create(paste(loc, type, str_split_fixed(src, "/", 3)[1], sep = "/"))
    dir.create(paste(loc, type, str_split_fixed(src, "/", 3)[1], str_split_fixed(src, "/", 3)[2], sep = "/"))
  }

  # define where everything will be saved
  files <- list(
    cnv = paste0("partials/cnv/", type, "/", str_remove(src, "\\.gz"), "-", tum, "_cnvs.gz"),
    raw_protein = paste0(annoloc, type, "/", str_remove(src, ".gz"), "-", tum, "_proteins.txt"),
    gene  =  paste0(annoloc, type, "/", str_remove(src, ".gz"), "-", tum, "_cnvs.txt"),
    table = paste0(annoloc, type, "/", str_remove(src, ".gz"), "-", tum, "_cnvs.tsv"),
    protein = paste0(annoloc, type, "/", str_remove(src, ".gz"), "-", tum, "_protein_sequences.txt"),
    dna = paste0(annoloc, type, "/", str_remove(src, ".gz"), "-", tum, "_dna_sequences.txt"))


  if (reuse & all(lapply(files, file.exists) %>% unlist())) {
    message("All files exist for ", tum) # don't redo anything if everything exists, just move on
  } else {
    if (reuse & file.exists(files$cnv)) {
      message("Loading annotate cnv file ", files$cnv)
      load(files$cnv)
    } else {
      print(paste0("=== ANNOTATING FOR TUMOR ", tum))
      # add annotations

      cnv <- cnvr(segmented[[tum]]) # should be loaded in
      # cnv <- annotateGenes(cnv, txdb)
      # cnv <- annotateTranscripts(cnv, txdb, tum)
      cnv <- annotateAll(cnv, protein_table, type = type)

      # and then save the annotations
      message("Saving annotation file ", files$cnv)
      save(cnv, file = files$cnv, compress = TRUE)
    }


    cnv <- as.data.frame(cnv) %>% mutate(throwout = toThrowout(start) | toThrowout(end))
    if (grepl(throwout, tum)) {
      print(paste("Threw out", count(cnv %>% filter(throwout)), "regions"))
      print(cnv %>% filter(throwout))
      cnv <- cnv %>% filter(!throwout)
    }

    if (!(reuse & file.exists(files$gene))) {
      message("saving gene file at ", files$gene)
      write(cnv$locus %>% unlist(), files$gene) # only gene information
    }
    if (!(reuse & file.exists(files$raw_protein))) {
      message("saving raw proteins names for STRINGDB at ", files$raw_protein)
      write(cnv$protein_name_short %>% unlist(), files$raw_protein) # only gene information
    }
    if (!(reuse & file.exists(files$table))) {
      message("saving cnv table at ", files$table)
      write.table( # more full data, allows recreation
        cnv,
        file = files$table,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
    }
    # and then save the protein and dna sequences based on the annotated genes
    if (!(reuse & file.exists(files$protein))) {
      message("Saving cnv protein sequences at ", files$protein)
      saveSeq(cnv, aa, files$protein)
    }
    if (!(reuse & file.exists(files$dna))) {
      message("Saving cnv dna sequences at ", files$dna)
      saveSeq(cnv, dna, files$dna)
    }
  }

}
#### 4. Run Pipleline ####
partloc <- paste0("partials/segmented/", bp, "/")
plotloc <- paste0("01_CNV-Plots/", bp)

dir.create(partloc, showWarnings = FALSE)
dir.create(plotloc, showWarnings = FALSE)

# for (bp in c(5000, 1000)) {
# here the point is to run multiple norm options
for (normType in c("poisson", "mode", "mean", "median")) {
  dir.create(paste0("partials/segmented", bp, normType, sep = "/"), showWarnings = FALSE)
  dir.create(paste0("partials/normalized", bp, normType, sep = "/"), showWarnings = FALSE)
  dir.create(paste(plotloc, normType, sep = "/"), showWarnings = FALSE)

  # with multiple possible size factors
  for (sizeFactor in c("mean", "median", "mode")) {
    dir.create(paste0("partials/segmented", bp, normType, sizeFactor, sep = "/"), showWarnings = FALSE)
    dir.create(paste0("partials/normalized", bp, normType, sizeFactor, sep = "/"), showWarnings = FALSE)
    dir.create(paste(plotloc, normType, sizeFactor, sep = "/"), showWarnings = FALSE)
    for (filter in c(0, 0.01, 0.05)) { # and different filters for across samples
      dir.create(paste(plotloc,normType,  sizeFactor, filter, sep = "/"), showWarnings = FALSE)
      for (byChrom in c(TRUE, FALSE)) { # either by chrom or not by chrom
        # then we can normalize and save each option
        dir.create(paste(plotloc, normType, sizeFactor, filter, byChrom, sep = "/"), showWarnings = FALSE)
        norm_grange <- prepCounts(sample_data,
          bp = bp, normType = normType,
          sizeFactor = sizeFactor, filter = filter,
          byChrom = byChrom
        )
        # then the only seg option, minReadCount
        for (minReadCount in c(4, 8, 10, 16, 20, 32)) {
          dir.create(paste(plotloc, normType, sizeFactor, filter, byChrom, minReadCount, sep = "/"), showWarnings = FALSE)
          partial_footer <- paste0(bp, "/", normType, "/",
          sizeFactor, "-", filter, "-", byChrom, "-", minReadCount, ".gz")
          seg_file <- paste0(
            "partials/segmented/",partial_footer
          )
          if (reuse & file.exists(seg_file)) {
            message("Complete segmentation file exists for minReadCount ", minReadCount, ".")
          } else {
            message("Running segmentation for minReadCount ", minReadCount)
            dir.create("partials/segmented/temp")
            segmented <- lapply(names(norm_grange), callCNVfromNorm, norm_full = norm_grange, minReadCount = minReadCount)
            unlink("partials/segmented/temp/", recursive = TRUE)
            names(segmented) <- names(norm_grange)
            save(segmented, file = seg_file)
          }
          # ensure directories exist


          tum <- names(norm_grange)[length(norm_grange)]
        base_loc <- paste(plotloc, normType, sizeFactor, filter, byChrom, minReadCount, "", sep = "/")
          seg_output_files <-  list(chromplot = paste0(base_loc, tum, " FP - chromosome segplot.png"),
                                    segplot = paste0(base_loc, tum, " FP - segplot.pdf"),
                                    size = paste0(base_loc, tum, " FP - Size Distribution.png"),
                                    type = paste0(base_loc, tum, " FP - Type Distribution.png"),
                                    dna = paste0(annoloc, "any", "/", str_remove(partial_footer, ".gz"), "-", names(norm_grange)[length(names(norm_grange))], "_dna_sequences.txt"))

          if(reuse & all(lapply(seg_output_files, file.exists) %>% unlist())){
            message("All output files exists, skipping annotation and plotting.")
          } else {
            message("Loading segmentation")
            

            load(seg_file)
            message("Plotting segmentation")
          lapply(names(segmented), plotCNV,
            segmented = segmented,
            base_loc = base_loc
          )
          message(paste("Starting annotation for within on", paste(names(segmented), collapse = " ")))
          lapply(names(segmented), cnvDownstream, src = partial_footer, segmented = segmented, type = "within")
          message(paste("Starting annotation for any on", paste(names(segmented), collapse = " ")))
          lapply(names(segmented), cnvDownstream, src = partial_footer, segmented = segmented, type = "any")
          }
        }


        gc()
      }
    }
  }
}
# }
