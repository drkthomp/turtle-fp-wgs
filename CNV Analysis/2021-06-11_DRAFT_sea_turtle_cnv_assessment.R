####################################################
# cnMOPS-based preliminary copy number analysis    #
#                                                  #
# TCG Cambridge, Department of Veterinary Medicine #
# 03.02.2021                                       #
# Max Stammnitz, mrs72@cam.ac.uk                   #
####################################################

# Libraries and Dependencies
library(cn.mops)
library(tidyverse)
library(janitor)

# 0. Load Host Information
##########################
sample_data <-
  read.delim("~/2021-REU/Genomic_Samples_Table.xlsx - Sheet1.tsv") %>%
  clean_names() %>%
  mutate(sample_name = str_replace_all(sample_name, " - HiSeq", "")) %>%
  filter(type == "New")

names(sample_data)[names(sample_data) == "internal_external"] <-
  "tumor_src"
names(sample_data)[names(sample_data) == "tumor_non_tumor"] <-
  "tissue_status"

sample_data <- sample_data %>% mutate(turtle = case_when((turtle == "Yucca" & tumor_sample_location == "kidney ") ~ "Yucca_kidney", TRUE ~ turtle))


names(sample_data)[names(sample_data) == "non-tumor"] <- "normal"

tumor_data <- sample_data %>% filter(tissue_status == "tumor")

match_data <- sample_data %>%
  pivot_wider(names_from = tissue_status, values_from = sample_name, id_cols = turtle) %>%
  unnest(c(tumor, `non-tumor`))

names(match_data)[names(match_data) == "non-tumor"] <- "normal"

sample_data <- left_join(tumor_data, match_data, by = c("sample_name" = "tumor")) %>% distinct()

# Arguments
windowsize <- 5000
reference.scaffolds <- "chelonia_mynas_new_referencescaffolds.txt"
output.path <- "plots/"

###############
# 1. Read Counts
#################

#  Specify scaffolds
reference.scaffolds <- read.table(reference.scaffolds, header = T, check.names = F, sep = "\t")
reference.scaffolds <- reference.scaffolds[-c(29:98), ]
reference.scaffolds <- as.character(reference.scaffolds[, "RefSeq-Accn"])


#  Read in count matrix
get_matched <- function(x) {
  x <- as_tibble_row(x, .name_repair = "universal")
  print(x)
  if (!is_null(x$sample_name) && !is_null(x$normal)) {
    normal_file <-
      paste0(
        "~/2021-REU/CNV Analysis/readcounts/readcounts_",
        windowsize,
        "BP_",
        x$normal,
        ".Rdata"
      )
    if (file.exists(normal_file)) {
      print(normal_file)
      load(normal_file)

      assign(x$normal, counts.mt, envir = .GlobalEnv)
      for (tumor in x$sample_name) {
        tumor_file <-
          paste0(
            "~/2021-REU/CNV Analysis/readcounts/readcounts_",
            windowsize,
            "BP_",
            tumor,
            ".Rdata"
          )
        if (file.exists(tumor_file)) {
          print(tumor_file)
          load(tumor_file)
          assign(tumor, counts.mt, envir = .GlobalEnv)
        }
      }
    }
    # load(paste0("~/2021-REU/CNV Analysis/readcounts/readcounts_",windowsize,"BP_", x$tumor, ".Rdata"))
    # assign(counts.mt, x$tumor)
  }
  # load("~/2021-REU/CNV Analysis/readcounts/readcounts_5000BP_flSCBLdna.Rdata")
}

apply(sample_data, 1, get_matched)


##########
# 2. Plots
##########

# DOC plots

plot_DLOG <- function(paired, tumor, normal) {
  if (length(tumor) > 0 &&
    exists(tumor) && length(normal) > 0 && exists(normal)) {
    counts.tumor <- as.data.frame(get(tumor))
    counts.normal <- as.data.frame(get(normal))


    counts <-
      left_join(
        counts.tumor,
        counts.normal,
        suffix = c("_TUMOR", "_NORMAL"),
        by = c("SCAFFOLD", "START", "END")
      ) %>%
      clean_names() %>%
      mutate_at(c("sample_counts_tumor", "sample_counts_normal", "start", "end"), as.numeric)

    counts <- counts %>%
      mutate(sample_counts = log(sample_counts_tumor / sample_counts_normal, base = 2)) %>%
      drop_na(sample_counts) %>%
      filter(sample_counts != Inf)
    assign(paste0(tumor, ".", normal), counts, envir = .GlobalEnv)
    png(
      paste0(
        output.path,
        "doc_plot_",
        windowsize,
        "BP",
        paired,
        "_",
        tumor,
        "_",
        normal,
        ".png"
      ),
      width = 1400,
      height = 1000
    )
    mar.default <- c(5, 4, 4, 2) + 0.1
    par(mar = mar.default + c(0, 2, 0, 0))

    ## "ylim" parameter sets an upper threshold to visible reads per location;
    ## "cex" parameter can increase/decrease point size
    ## both need to be adjusted, depending on the chosen windowsize and sequencing coverage
    plot(
      counts$sample_counts,
      ylim = c(-4, 4),
      pch = 16,
      cex = 0.2,
      col = "darkgreen",
      ylab = "Log2(Tumor/Host)",
      cex.lab = 1.3,
      main  = paste0(str_replace_all(paired, "_", " "), " read counts in bins of ", windowsize, " BP"),
      cex.main = 2,
      xlab = "Bin index along Chelonia mynas scaffolds"
    )
    dev.off()
  }
}
plot_DOC <- function(x) {
  if (length(x) > 0 && exists(x)) {
    counts.mt <- get(x)
    png(
      paste0(output.path, "doc_plot_", windowsize, "BP", x, ".png"),
      width = 1400,
      height = 1000
    )
    mar.default <- c(5, 4, 4, 2) + 0.1
    par(mar = mar.default + c(0, 2, 0, 0))

    ## "ylim" parameter sets an upper threshold to visible reads per location;
    ## "cex" parameter can increase/decrease point size
    ## both need to be adjusted, depending on the chosen windowsize and sequencing coverage
    print(max(counts.mt[, "SAMPLE COUNTS"]))
    plot(
      counts.mt[, "SAMPLE COUNTS"],
      ylim = c(0, 1e6),
      pch = 16,
      cex = 0.2,
      col = "darkblue",
      ylab = "Read counts",
      cex.lab = 1.3,
      main  = paste0(x, " read counts in bins of ", windowsize, " BP"),
      cex.main = 2,
      xlab = "Bin index along Chelonia mynas scaffolds"
    )
    dev.off()
  }
}


plot_DOCs <- function(x) {
  x <- as_tibble_row(x, .name_repair = "universal")
  for (tumor in x$sample_name) {
    tum_name <- str_replace(tumor, "_", "-")
    plot_DOC(as.character(tumor))
    plot_DLOG(
      paste(str_replace(filter(match_data, tumor == tum_name)$turtle, "_kidney", ""), filter(sample_data, sample_name == tum_name)$tumor_sample_location),
      as.character(tumor),
      as.character(x$normal)
    )
  }
  plot_DOC(as.character(x$normal))
}

apply(sample_data, 1, plot_DOCs)
