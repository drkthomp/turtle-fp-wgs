# SOMATYPUS: A PLATYPUS-BASED VARIANT CALLING PIPELINE FOR CANCER DATA
# Adrian Baez-Ortega, Transmissible Cancer Group, University of Cambridge
# 03/06/2016

# BasicManipulation.R
# Loads variant data extracted with ExtractVcfData.py, computes VAF, generates tumour-only
# and host (non-tumour) variant sets, and produces some basic exploratory plots.
# For this script to work, tumour samples must contain a 'T' in their sample names.
# Only vars are used in this example.
library(tidyverse)
library(janitor)


# Read data extracted with ExtractVcfData.py
setwd("~/2021-REU/SNV")
what <- "Indels"
header <- scan(file=paste0("Somatypus_", what, "_final_NR.txt"),nlines=1,what=character())
vars.metadata <- read.table(paste0("Somatypus_",what,"_final_Metadata.txt"), header=T,nrows = 212450)
var.ids <- paste0(vars.metadata$CHROM, ":", vars.metadata$POS, vars.metadata$REF,">", vars.metadata$ALT)
vars.nr <- matrix(scan(file=paste0("Somatypus_", what, "_final_NR.txt"),skip=1,what=numeric()),ncol = length(header),dimnames = list(var.ids, header))
vars.nv <- matrix(scan(file=paste0("Somatypus_", what, "_final_NV.txt"),skip=1,what=numeric()),ncol = length(header),dimnames = list(var.ids, header))


# 1) Compute VAF (nv/nr)
vars.vaf <- vars.nv / vars.nr
rm(vars.nr)
rm(vars.nv)
vars.vaf[is.na(vars.vaf)] <- 0
gc()

# 2) Tidy the data
data <- cbind(vars.vaf, vars.metadata)
names <- colnames(vars.vaf)
sample_names <- str_remove_all(names, "bwa|_|CheMyd|sorted|HiSeq")
vcfcols <- colnames(vars.metadata)
#rm(vars.vaf)
rm(vars.metadata)
gc()




load("~/2021-REU/CNV Analysis/partials/sample_data")
sample_data_all <-
  read.csv("~/2021-REU/samplelist.csv") %>% clean_names()%>% mutate(sample_name = str_replace_all(sample_name, " - HiSeq", "")) %>% mutate(type = grepl("FP", tissue)) %>% mutate(sample_name = str_split(sample_name, " ",n = 1))

tumours <- grepl(paste((sample_data_all %>% filter(type))$sample_name, sample_data$sample_name, sep="|",  collapse="|"), names,ignore.case = TRUE) & !grepl("HiSeq", names, ignore.case=TRUE)
hosts <- grepl(paste((sample_data_all %>% filter(!type))$sample_name, sample_data$normal, sep="|", collapse="|"), names,ignore.case=TRUE) & !grepl("HiSeq", names, ignore.case=TRUE)
names[!(tumours | hosts)] # make sure we've labeled everything except HiSeq

# select out the tumour and host variants
tumor_variants <- data %>% select(c(all_of(vcfcols),names[tumours]))
host_variants <- data %>% select(c(all_of(vcfcols),names[hosts]))
rm(data)

vafAdjust <- function(x){
  normal <- sample_data[str_to_upper(sample_data$sample_name)  == str_to_upper(sample_names[x]),]$normal
  if(length(normal) > 0){
    normal_file <- names[grepl(normal, names,ignore.case=TRUE)][1]
    tumor_vaf <- tumor_variants[names[x]]
    normal_vaf <- host_variants[normal_file]
    adjusted <- tumor_vaf - normal_vaf
    return(adjusted)
  }
  return(NULL)
}

vaf_adjust <- lapply(seq_along(sample_names)[tumours], vafAdjust)
rm(host_variants)
names(vaf_adjust) <- sample_names[tumours]
gc()
vaf_adjusted <- bind_cols(vaf_adjust)
#save(vaf_adjusted, file="vaf_adjusted")
#load("vaf_adjusted")
final <- cbind(vaf_adjusted, tumor_variants %>% select(vcfcols))

final$mean <- final %>% select(starts_with("bwa")) %>% abs() %>% mutate(mean = rowMeans(.,na.rm = TRUE)) %>% select(mean)
final$total <- final %>% select(starts_with("bwa")) %>% abs() %>% adorn_totals("col") %>% select(Total)
final <- final %>% mutate(id = paste0(CHROM, ":", POS, REF, ">", ALT))

# only actually signifigant ones
#final_filter <- final %>% filter(if_any(starts_with("bwa"), ~ abs(.) > 0.25 ))

#final_filter_tidy <- final_filter %>% pivot_longer(starts_with("bwa")) %>% mutate(name = str_remove_all(name, "bwa|_|CheMyd|sorted|HiSeq"))
#final_filter_tidy <- final_filter_tidy %>% mutate(id = paste0(CHROM, ":", POS, REF, ">", ALT))
#save(final_filter_tidy, file="final_filter")
#load("final_filter")
#final_filter_tidy <- final_filter_tidy %>% filter(abs(value) > 0.25) %>% arrange(-abs(value))
#ggplot(final_filter_tidy %>% filter(id %in% head(final, 100)$id), aes(x=name,y=reorder(id,value),fill=value)) + geom_tile() +
#  theme_void() + scale_fill_viridis_c() +
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank(),
#  axis.text.x = element_text(angle = 45, hjust=0.9,vjust=1))

library(devtools)
library(ComplexHeatmap)
colnames(vaf_adjusted) <-  str_remove_all(colnames(vaf_adjusted), "bwa|_|CheMyd|sorted|HiSeq")
rownames(vaf_adjusted) <- rownames(final)
#save(vaf_adjusted, file="vaf_adjusted")
rows <- unlist(head(rownames(final %>% arrange(-mean)), 100000))
library(magick)

pdf(file=paste0(what,"_10k_Heatmap.pdf"),width=14,height=17)
Heatmap((vaf_adjusted[rownames(vaf_adjusted) %in% head(rows,50000),]),name = "Tumor - Host VAF", cluster_rows = FALSE,show_row_names = FALSE,column_km = 6,use_raster=FALSE,row_title=paste("Top 50k (by mean)", what,"in CHR1"),column_dend_height=unit(20,"mm"))
dev.off()
pdf(file=paste0(what,"_100k_Heatmap.pdf"),width=14,height=17)
Heatmap((vaf_adjusted[rownames(vaf_adjusted) %in% rows,]),name = "Tumor - Host VAF", cluster_rows = FALSE,show_row_names = FALSE,column_km = 6,use_raster=FALSE,row_title=paste("Top 100k (by mean)", what,"in CHR1"),column_dend_height=unit(20,"mm"))
dev.off()



