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

# setwd("/Users/ms37/Desktop/Labwork/Side Projects/Turtles/src")

##########################
# 0. Load Host Information
##########################
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
  left_join(sample_data %>% filter(tissue_status == "tumor"), match_data, by = c("sample_name" = "tumor")) %>%
  select(-"turtle.y") %>%
  distinct()

# and save it for use in other scripts for simplicity
save(sample_data, file = "partials/sample_data")

##########################
# 1. Load Data Imports
##########################

dataloc <- "readcounts/" # note where the readcounts are
bp <- 5000 # then bp we're using
# and the final file name before the tumor:
file_prefix <- paste0(dataloc, "readcounts_", bp, "BP_")

# list of valid tumor files to be filled
tumors <<- c()

# get relevant reference data
txdb <- loadDb("~/2021-REU/CNV Analysis/rCheMyd1.sqlite")
ref <- getChromInfoFromNCBI("GCF_015237465.1", assembled.molecules.only = TRUE, assembly.units = "Primary Assembly")
ref_seq <- getChromInfoFromNCBI("GCF_015237465.1", as.Seqinfo = TRUE, assembled.molecules.only = TRUE, assembly.units = "Primary Assembly")
seqnames(ref_seq) <- ref$RefSeqAccn
chr_ref <-
  mutate(ref, CHR = parse_number(SequenceName)) %>% select("RefSeqAccn", "CHR")

##########################
# 2. Define functions
##########################

# function to read in paired tumor and normal read count files
read_in <- function(x) {
  tumor_file <-
    paste0(file_prefix, x[names(x) == "sample_name"], ".txt")
  normal_file <-
    paste0(file_prefix, x[names(x) == "normal"], ".txt")

  print(paste("Reading in", tumor_file, "and", normal_file))

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
    tumors <<-
      c(tumors, as.character(str_replace_all(x[names(x) == "sample_name"], "-", "_")))

    return(final)
  }
}

# wrapper function to call segmentation
callCounts <- function(norm) {
  segmented <- referencecn.mops(
    cases = norm[, "TUMOUR"],
    controls = norm[, "HOST"],
    norm = 0,
    normType = "mode",
    minReadCount = 3,
    segAlgorithm = "fast"
  )
  return(calcIntegerCopyNumbers(segmented))
}


# segmentation plot
nice.seg.plot <-
  function(x.segmented,
           x.normalised.counts,
           title = "",
           file_title = "",
           annotate = FALSE, # try to tag the annotated genes (visually busy)
           adjustForChrom = TRUE, # adjust the positions to not overlap by chrom
           indicateChrom = FALSE, # must adjust to indicate
           colorChrom = FALSE) { # color chrom not recommending for display
    ratio <-
      elementMetadata(x.normalised.counts)[, "TUMOUR"] / elementMetadata(x.normalised.counts)[, "HOST"]
    ratio <- ratio * 2

    pdf(paste0("plots/", file_title, ".pdf"),
      width = 14,
      height
      = 7
    )
    mar.default <- c(5, 4, 4, 2) + 0.1
    par(mar = mar.default + c(0, 3, 0, 0))

    x_data <- as.matrix(ranges(x.normalised.counts))[, 1]
    # if we're converting the x by chromosome
    if (adjustForChrom) {
      # get absolute lengths
      to_adjust <- as.data.frame(x.normalised.counts)
      adjusts <- list()
      adj <- 0
      for (j in length(levels(seqnames(norm))):1) {
        adj <- adj + max(0, max((to_adjust %>% filter(as.numeric(seqnames) == j + 1))$start))
        adjusts[j] <- adj
      }
      # and then apply them to the data
      to_adjust <- mutate(to_adjust, start = start + as.numeric(adjusts[as.numeric(seqnames)]))
      to_adjust <- mutate(to_adjust, start = abs(start - max(start)))

      to_adjust <- mutate(to_adjust, end = start + width)


      ranges(x.normalised.counts) <- IRanges(
        names = to_adjust$seqnames,
        start = to_adjust$start,
        end = to_adjust$end,
        # width=to_adjust$width,
        strand = to_adjust$strand,
        TUMOUR = to_adjust$TUMOUR,
        HOST = to_adjust$HOST
      )
      x_data <- to_adjust$start
    }


    colors <- "grey"

    if (colorChrom) {
      colors <- seqnames(norm)
      # 28 most distinct colors. not pretty but does the job
      levels(colors) <- c(
        "#696969", "#8b4513", "#006400", "#808000", "#483d8b", "#008b8b",
        "#00008b", "#7f007f", "#8fbc8f", "#b03060", "#ff4500", "#ffa500",
        "#ffff00", "#00ff00", "#8a2be2", "#dc143c", "#00ffff", "#0000ff",
        "#adff2f", "#da70d6", "#ff00ff", "#1e90ff", "#f0e68c", "#fa8072",
        "#90ee90", "#87ceeb", "#ff1493", "#ffc0cb"
      )
    }

    # plot the background points
    plot(
      x = x_data,
      y = ratio,
      type = "p",
      xlab = "Genome Position",
      ylab = "Copy Number",
      pch = 16,
      cex = 0.1,
      col = as.character(colors),
      ylim = c(0, 5),
      main = title,
      cex.lab = 2,
      cex.main = 3,
    )

    # if indicating, add lines for chromosomes
    if (adjustForChrom & indicateChrom) {
      pos <- 0
      for (i in seq_along(levels(seqnames(x.normalised.counts)))) {
        prev <- pos
        pos <- pos + seqlengths(x.normalised.counts)[i]
        lines(x = rep(pos, 2), y = c(-5, 5), col = "lightgrey", lwd = 1)
        text(y = 4.5, x = (prev + pos) / 2, labels = i)
      }
    }

    # add segmentation lines in red and label
    segmented_segs <- as.data.frame(segmentation(x.segmented))
    segs <- as.data.frame(cnvs(x.segmented))

    if (adjustForChrom) {
      segs <- mutate(segs, start = start + as.numeric(adjusts[as.numeric(seqnames)]))
      segs <- mutate(segs, start = abs(start - max(start)))

      segs <- mutate(segs, end = start + width)
    }

    segs[, "CN"] <-
      as.numeric(str_split_fixed(segs[, "CN"], "N", 2)[, 2])
    segs <- arrange(segs, start)

    for (i in 1:nrow(segs)) {
      lines(
        x = c(segs[i, "start"], segs[i, "end"]),
        y = rep(segs[i, "CN"], 2),
        col = "red",
        lwd = 2
      )

      if ((annotate) & as.numeric(segs[i, "CN"]) != 2) {
        annotated_cnvs <- annotateCnvs(makeGRangesFromDataFrame(segmented_segs[i, ]), txdb)$gene_id@unlistData
        annotated_cnvs <-
          grep("^LOC",
            annotated_cnvs,
            invert = TRUE,
            value = TRUE
          )
        if (length(annotated_cnvs) > 0) {
          print(paste("Plotting for real at", segs[i, "start"]))
          text(
            x = segs[i, "start"],
            y = segs[i, "CN"],
            pos = 3,
            adj = c(0, 0),
            cex = 0.2,
            srt = 90,
            labels = paste0(annotated_cnvs, collapse = ", ")
          )
        }
      }
    }

    dev.off()
  }



karyoPlot <- function(segmented, title, file_title) {
  cn.calls <- loadCopyNumberCalls(segmentation(segmented), cn.col = "TUMOUR")

  # adjust to numeric chromosome
  numbers <- as.numeric(names(sort(seqlengths(BSgenome.Cmydas.NCBI.rCheMyd1), decreasing = TRUE)))
  numbers <- numbers[numbers < 29]
  seqlevels(cn.calls) <- as.character(numbers)

  seqlevels(cn.calls) <- seqlevels(BSgenome.Cmydas.NCBI.rCheMyd1)
  seqinfo(cn.calls) <- seqinfo(BSgenome.Cmydas.NCBI.rCheMyd1)
  # and change cn calls to numeric
  cn.calls$cn <- parse_number(cn.calls$CN)

  pdf(paste0("plots/", file_title, "_karyo.pdf"),
    width = 14,
    height
    = 7
  )
  kp <- plotKaryotype(genome = "BSgenome.Cmydas.NCBI.rCheMyd1", chromosomes = 1:28)
  plotCopyNumberCalls(kp, cn.calls, cn.colors = "red_blue", r1 = 0.8, label.cex = 0.2)
  cn.cols <- getCopyNumberColors(colors = "red_blue")
  title(main = title, sub = "note chr24-28 are mislabeled", adj = 1, line = -4.5)
  legend("right", legend = names(cn.cols), fill = cn.cols, ncol = length(cn.cols))
  dev.off()
}

# requires an input with
# a sample_name (matching read count file)
# a normal (matching normal read count file)
# a turtle name and a tumor_sample_location
cnvAnalysis <- function(x) {
  tum <- x[names(x) == "sample_name"]
  normal <- x[names(x) == "normal"]
  tum_name <- str_replace(tum, "_", "-")

  # Read in
  gr <- read_in(x)
  if (!is.null(gr)) { # ensure that we can proceed
    print(paste0("=== NORMALIZING FOR TUMOR ", tum))
    norm <- normalizeGenome(gr)
    # and save the call for future use
    save(norm, file = paste0("partials/normalized/", tum, "-", normal, ".gz"), compress = TRUE)

    # Call counts
    print(paste0("=== CALLING FOR TUMOR ", tum))
    segmented <- callCounts(norm)
    # and save the call for future use
    save(segmented, file = paste0("partials/segmented/", tum, "-", normal, ".gz"), compress = TRUE)


    # cn.mops style chromosome graph
    png(
      file = paste0(
        "plots/",
        tum,
        " FP - segplot.png"
      ),
      width = 3000,
      height = 1500
    )
    segplot(segmented, lwd = 1, ylim = c(-5, 5), plot.type = "w")
    dev.off()



    print(paste0("=== ANNOTATING FOR TUMOR ", tum))
    cnv <- cnvr(segmented)
    cnv <- annotateGenes(cnv, txdb)
    cnv <- annotateTranscripts(cnv, txdb)
    write(cnv$gene_id %>% unlist(), paste0("data/", tum, "_cnvs.txt"))
    write.table(cnv, file = paste0("data/", tum, "_cnvs.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

    print(paste0("=== PLOTTING FOR TUMOR ", tum))
    title <- paste( # turtle tissue FP - segmented CNVs
      str_replace(filter(match_data, tumor == tum_name)$turtle, "_kidney", ""),
      filter(sample_data, sample_name == tum_name)$tumor_sample_location,
      "FP - segmented CNVs"
    )
    file_title <- paste0( # tum_normal FP - segmented CNVs(.pdf)
      tum,
      "_",
      filter(match_data, tumor == tum_name)$normal,
      " FP - segmented CNVs"
    )
    # karyoPlot(segmented,title,file_title)
    # Plotting
    if ((!exists("doNiceSegPlot")) | doNiceSegPlot) { # nice seg plot takes a while, sometimes want to not rerun

      nice.seg.plot(
        x.segmented = segmented,
        x.normalised.counts = norm,
        indicateChrom = FALSE, # don't indicate the chromosomes (busy)
        colorChrom = FALSE, # def don't color the chromosomes (busy, no legend)
        annotate = TRUE, # do add small annotations
        title = title,
        file_title = file_title
      )
    }
  }
}

doNiceSegPlot <- FALSE

run <- sample_data
# and go and do it for everything!
apply(run, 1, cnvAnalysis)
