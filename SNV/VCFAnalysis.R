setwd("~/2021-REU/SNV")
#BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
library(GenomicRanges)
library(tidyverse)
library(janitor)
library(GenomicFeatures)
load("~/2021-REU/CNV Analysis/partials/sample_data")
sample_data_all <-
  read.csv("~/2021-REU/samplelist.csv") %>% clean_names()%>% mutate(sample_name = str_replace_all(sample_name, " - HiSeq", "")) %>% mutate(type = grepl("FP", tissue)) %>% mutate(sample_name = str_split(sample_name, " ",n = 1))

#ref<- getChromInfoFromNCBI("GCF_015237465.1")$RefSeqAccn
genes <- genes(loadDb("~/2021-REU/CNV Analysis/rCheMyd1.sqlite"))
#genes_by_chr <- splitAsList(genes, seqnames(genes))
#rm(genes)
fl <- "inputs/Somatypus_Indels_final.vcf.gz"

# need to load one chromosome at a time or its too big! + currently subsetting to exons
# only need the NR and NV data from GENO, nothing from INFO

convertToMatrix <- function(x){
  rows <- rownames(x)
  cols <- colnames(x)
  data <- as.numeric(x)
  dim(data) <- c(length(rows),length(cols))
  dimnames(data) <- list(rows, cols)
  return(data)
}

### PLOTTING FUNCTIONS
getX <- function(ind, tumor=TRUE){ # getX for VAF
  variants <- names(vars.vaf[tumour.only.idx == tumor, ind])
  x <- ranges[names(ranges) %in% variants]
  x <- x[match(variants, names(x)),]
  return(start(x))
}

plotVAF <- function(ind,ylab=NA){ # plot the vaf
  # plot the normals
  plot(getX(ind, tumor=FALSE), vars.vaf[!tumour.only.idx, ind],
       ylab=ylab, xlab=names[ind], pch=20, col="gray72", cex=0.6,xlim=xlim)
  # then the tumor in red
  points(getX(ind, tumor=TRUE), vars.vaf[tumour.only.idx, ind[1]], pch=18, col="red", cex=0.7)
}

plotCoverage <- function(ind,ylab=NA){
  # Coverage (log10): all variant
  x <- ranges
  y <- log10(vars.nr[,ind])
  plot(start(x[match(names(y), names(x)),]), y, ylab=ylab,
       xlab=names[ind], pch=20, col="gray72", cex=0.7,ylim=c(0,2),xlim=xlim)
}

plotFull <- function(toPlot,name){
  i <- seq_along(samples)[toPlot]
  # Create PNG file
  png(paste0(length(vars.vaf[tumour.only.idx,]), "_", names(gene), "_", name, ".png"), 2400, 1000)
  par(mfrow=c(2,length(samples[toPlot])))
  plotVAF(i[1],ylab="VAF")
  mtext(paste("VAF for gene", names(gene), "in sample"), outer = TRUE, side=3,line=-2)
  for (j in i[2:length(i)]) {
    plotVAF(j)
  }

  plotCoverage(i[1],ylab="log10(Coverage)")

  for (j in i[2:length(i)]) {
    plotCoverage(j)
    # Close PNG
  }
  mtext(paste("Coverage for gene", names(gene), "in sample"),outer = TRUE, side=3,line=-55)
  dev.off()
}

start_gene <- 304
display_percent <- 100
for(g in seq_along(genes)[start_gene:length(genes)]){
  gene <- genes[g]
  curr_percent <-round(g/length(genes) * 100,digits=2)
  if(display_percent != curr_percent){
  message(names(gene), " - ",g, "/",length(genes), "(", curr_percent, "%)")
    display_percent <- curr_percent
  }
  gc()
  setwd("~/2021-REU/SNV")
  tryCatch({
    svp <- ScanVcfParam(geno=c("NR","NV"), info=NA, which=gene)
    vcf <- readVcf(fl,as.character(seqnames(gene)), svp)

    vars.nv <- convertToMatrix(geno(vcf)$NV)
    if(length(vars.nv) > 1){
      vars.nr <- convertToMatrix(geno(vcf)$NR)

      vars.vaf <- vars.nv/vars.nr

      ### Get samples
      samples <- colnames(vars.nv)
      names <- str_remove_all(samples, "bwa_|_CheMyd|_sorted")
      tumours <- grepl(paste((sample_data_all %>% filter(type))$sample_name, sample_data$sample_name, sep="|",  collapse="|"), samples,ignore.case = TRUE) & !grepl("HiSeq", samples, ignore.case=TRUE)
      hosts <- grepl(paste((sample_data_all %>% filter(!type))$sample_name, sample_data$normal, sep="|", collapse="|"), samples,ignore.case=TRUE) & !grepl("HiSeq", samples, ignore.case=TRUE)
      samples[!(tumours | hosts)] # make sure we've labeled everything except HiSeq


      ### BASIC MANIPULATION
      MIN.VAF <- 0.25
      not.hosts <- apply(vars.vaf[,hosts], 1, function(vaf) {
        all(vaf < MIN.VAF)
      })

      # Second, take the variants that are found with at least 3 supporting reads in at least one tumour
      MIN.READS <- 3
      one.tumour <- apply(vars.nv[,tumours], 1, function(nv) {
        any(nv >= MIN.READS)
      })

      # Extract tumour-only variants and host variants
      tumour.only.idx <- not.hosts & one.tumour
      tumour.only.idx[is.na(tumour.only.idx)] <- FALSE # when in doubt, mark false (host)

      tumour.only.snvs <- rowRanges(vcf)[tumour.only.idx,]
      host.snvs <- rowRanges(vcf)[!tumour.only.idx,]
      #rm(vcf)

      # 3) Take histogram showing number of samples the T-only variants are shared between

      # Get number of samples in which each variant is present (≥3 reads)
      #num.tumours <- apply(vars.nv[tumour.only.idx, tumours], 1, function(nv) {
      #  sum(nv >= MIN.READS)
      #})
      #hist(num.tumours, breaks=100, col="cornflowerblue", border="white", main="Number of tumours containing each T-only variant (≥3 reads)", xlab=NULL)


      setwd("outputs/BasicManipulation/")
      xlim <- c(start(range(gene)), end(range(gene)))
      ranges <-ranges(rowRanges(vcf))

      plotFull(tumours, "tumor")
      plotFull(hosts,"host")
      rm(vars.nr)
      rm(vars.nv)
      rm(vars.vaf)
    }
  },
  error=function(cond) {
    print(paste("ERROR: Problem with gene", names(gene), "number", g))
    print(cond)
    list <- dev.list()
    tryCatch(dev.off(),error=function(w){print("Nothing to close.")})
  })
}
