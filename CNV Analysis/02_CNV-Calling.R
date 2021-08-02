library(cn.mops)
library(GenomicRanges)
library(tidyverse)
library(janitor)
library(GenomicFeatures)

# Calls CNVs and performs basic filtering of the read counts

setwd("~/2021-REU/CNV Analysis")

#### 0. Define Settings ####
reuse <- TRUE # whether to reuse files
dataloc <- "readcounts/" # note where the readcounts are
bp <- 5000 # then bp we're using
# and the final file name before the tumor:
file_prefix <- paste0(dataloc, "readcounts_", bp, "BP_")

#### 1. Load Sample Data ####
loadSampleData <- function() {
  if (file.exists("partials/sample_data")) {
    message("Sample data file exists, loading.")
    # if we've already done this, just load it
    load("partials/sample_data")
  } else {
    message("Generating sample data file.")
    sample_data <-
      read.delim("~/2021-REU/Genomic_Samples_Table.xlsx - Sheet1.tsv") %>%
      clean_names() %>%
      mutate(sample_name = str_replace_all(sample_name, " - HiSeq", "")) %>%
      filter(type == "New")

    # correct column names
    names(sample_data)[names(sample_data) == "internal_external"] <-
      "tumor_src"
    names(sample_data)[names(sample_data) == "tumor_non_tumor"] <-
      "tissue_status"

    # and STATICALLY correct Yucca's kidneys to be analyzed seperately
    sample_data <-
      sample_data %>% mutate(turtle = case_when(
        (turtle == "Yucca" &
          tumor_sample_location == "kidney ") ~ "Yucca_kidney",
        TRUE ~ turtle
      ))


    # create only the match data for turtle, tumor, and non-tumor
    match_data <- sample_data %>%
      pivot_wider(
        names_from = tissue_status,
        values_from = sample_name,
        id_cols = turtle
      ) %>%
      unnest(c(tumor, "non-tumor"))
    # and correct the non-tumor to a more tidy normal
    names(match_data)[names(match_data) == "non-tumor"] <- "normal"

    # then combine the data for easy iteration
    sample_data <-
      left_join(sample_data %>% filter(tissue_status == "tumor"),
        match_data,
        by = c("sample_name" = "tumor")
      ) %>%
      select(-"turtle.y") %>%
      distinct()

    # and save it for use in other scripts for simplicity
    save(sample_data, "partials/sample_data")
  }
  return(sample_data)
}
sample_data <- loadSampleData()


#### 2. Load Reference Data ####
txdb <- loadDb("rCheMyd1.sqlite") # generated from rCheMyd1
ref <- getChromInfoFromNCBI("GCF_015237465.1",
  assembled.molecules.only = TRUE,
  assembly.units = "Primary Assembly"
)
ref_seq <- getChromInfoFromNCBI("GCF_015237465.1",
  as.Seqinfo = TRUE,
  assembled.molecules.only = TRUE,
  assembly.units = "Primary Assembly"
)
seqnames(ref_seq) <- ref$RefSeqAccn


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
      names(mcols(grange)) <- c(x, sample_data[sample_data$sample_name == x, ]$normal)
      return(cbind(as.data.frame(seqnames(grange)), as.data.frame(mcols(grange))))
    }) %>% purrr::reduce(cbind)
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
    if(file.exists(raw_file)){
      message("Raw file exists at ", raw_file)
      load(raw_file)
    }else{
    norm_grange <- apply(sample_data, 1, normCounts,
      bp = bp,
      normType = normType, sizeFactor = sizeFactor
    )
    names(norm_grange) <- sample_data$sample_name
    save(norm_grange, file=raw_file)
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
callCNVfromNorm <- function(norm, minReadCount = 5, parallel = 2) {
  segmented <- referencecn.mops(
    cases = norm[, "TUMOUR"],
    controls = norm[, "HOST"],
    minReadCount = minReadCount,
    segAlgorithm = "fast",
    norm = 0, # since we're inputting already normalized read counts
    parallel = parallel
  )
  segmented <- calcIntegerCopyNumbers(segmented)
  return(segmented)
}
###### 3.2.2 Filtering ######
# For the remainder list of CNV segments, go through each segment/tumour pair
# and remove the call if a tumour’s normalised median read count within
# it is not significantly above or below its background distribution
# across the whole genome/chromosome
filterCNVSignifigance <- function(segmented){
}

# remove the call if both tumour AND matched normal(s)
# deviate significantly from their background distributions
# (in this case it’s likely not a tumour-specific event)
filterCNVTumour <- function(segmented){

}
###### 3.2.3 Visualization ######
plotCNV <- function(tum, segmented, base_loc, lim = 5) {
  chromplot_file <- paste0(base_loc, "-", tum, " FP - chromosome segplot.png")
  segplot_file <- paste0(base_loc, "-", tum, " FP - segplot.pdf")
  if (!(reuse && file.exists(chromplot_file))) {
    message("Generating chromosome plot at ", chromplot_file)
    try({
    png(
      file = chromplot_file,
      width = 3000,
      height = 1500
    )
    segplot(
      segmented[[tum]],
      ylim = c(-lim, lim),
      plot.type = "s"
    )
    dev.off()
    })
  }
  if (!(reuse && file.exists(segplot_file))) {
    message("Generating chromosome plot at ", segplot_file)
    try({
    pdf(segplot_file,
      width = 18,
      height
      = 9
    )
    segplot(
      segmented[[tum]],
      ylim = c(-lim, lim),
      plot.type = "w",
      pt.cols = c("#333333", "#000000")
    )
    dev.off()
    })
  }
}


#### 4. Run Pipleline ####
partloc <- paste0("partials/segmented/", bp, "/")
plotloc <- paste0("01_CNV-Plots/", bp, "/")

dir.create(partloc)
dir.create(plotloc)

# here the point is to run multiple norm options
for (normType in c("poisson", "mode", "mean", "min", "median", "quant")) {
  dir.create(paste0("partials/segmented/", bp, "/", normType))
  dir.create(paste0("partials/normalized/", bp, "/", normType))

  # with multiple possible size factors
  for (sizeFactor in c("mean", "median", "mode", "quant")) {
    for (filter in c(0, 0.01, 0.05)) { # and different filters for across samples
      for (byChrom in c(TRUE, FALSE)) { # either by chrom or not by chrom
        # then we can normalize and save each option
        norm_grange <- prepCounts(sample_data,
          bp = bp, normType = normType,
          sizeFactor = sizeFactor, filter = filter,
          byChrom = byChrom
        )
        # then the only seg option, minReadCount
        for (minReadCount in c(4, 8, 10, 16, 20, 32)) {
          seg_file <-paste0(
            "partials/segmented/", bp, "/", normType, "/",
            sizeFactor, "-", filter, "-", byChrom, "-", minReadCount, ".gz"
          )
          if(reuse & file.exists(seg_file)){
            message("Segmentation file exists for minReadCount ", minReadCount)
            load(seg_file)
          } else{
            message("Running segmentation for minReadCount ", minReadCount)
          segmented <- lapply(norm_grange, callCNVfromNorm, minReadCount = minReadCount)
          save(segmented, file = seg_file)
          }

          message("Plotting segmentation")
          lapply(names(segmented), plotCNV,
            segmented = segmented,
            base_loc = paste0(
              plotloc, "/", normType, "/",
              sizeFactor, "-", filter, "-",
              byChrom, "-", minReadCount
            )
          )
          rm(segmented)
          gc()
        }
      }
    }
  }
}
